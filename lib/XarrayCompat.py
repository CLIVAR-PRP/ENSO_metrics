# -*- coding:UTF-8 -*-
"""
XarrayCompat.py
===============
Compatibility-layer replacements for the retired CDAT/UV-CDAT
``cdms2.TransientVariable`` interface and related axis/grid objects used by
ENSO_metrics.

This module is intended as a compatibility shim for CDAT-style diagnostic code.
It is most reliable for decoded, rectilinear gridded fields such as common
CMIP, ERA5, GPCP, and E3SM post-processed data with dimensions like
``(time, lat, lon)`` or ``(time, lev, lat, lon)``.

Important scope note
--------------------
This is not a full replacement for xarray, xESMF, ESMPy, or native
grid-aware analysis tools. Native unstructured or curvilinear grids, such as
E3SM ``ncol`` grids, MPAS grids, and ocean-model grids, should usually be
regridded or handled with grid-aware tools before conversion to this CDAT-like
object model.

All internal numerical computation is done with ``numpy.ma``. Coordinate
metadata are stored alongside the data so callers that use CDAT-style
introspection methods such as ``.getAxisList()``, ``.getGrid()``,
``.getTime().asComponentTime()``, and related accessors can continue to work
with minimal changes.

Exported public API
-------------------
CDATVariable   - CDAT-like variable object used in place of
                 ``cdms2.TransientVariable`` within the refactored workflow
_Axis          - CDAT-like coordinate axis object
_TimeAxis      - CDAT-like time-axis object with component-time support
_Grid          - CDAT-like rectilinear grid object

Factory helpers
---------------
create_axis(values, id='', units='', attributes=None, axis_type=None)
create_uniform_lat_axis(start, n, delta)
create_uniform_lon_axis(start, n, delta)
create_rect_grid(lat_axis, lon_axis, order='yx', grid_type='generic', mask=None)
create_variable(data, axes=None, grid=None, mask=None, id='', attributes=None)

Conversion helpers
------------------
da_to_cdat(da, varname=None)   - xr.DataArray  -> CDATVariable
cdat_to_da(var)                - CDATVariable  -> xr.DataArray
"""

from __future__ import annotations

import datetime as _datetime
from typing import List, Optional

import numpy as np
import numpy.ma as ma
import xarray as xr


# ---------------------------------------------------------------------------
# Module-level behaviour flags
# ---------------------------------------------------------------------------

# When True, _validate_grid raises ValueError instead of warning when lon
# precedes lat in the axis list (non-standard xy ordering).
#
#   import lib.XarrayCompat as XC; XC.STRICT_GRID = True
#
STRICT_GRID: bool = False

# When True, _detect_axis_type refuses to return a type for the time axis
# based on name heuristics alone (requires CF axis/standard_name/units).
# This prevents silent misclassification of e.g. latitude → T.
#
#   import lib.XarrayCompat as XC; XC.STRICT_AXIS_DETECTION = False
#
STRICT_AXIS_DETECTION: bool = True


__all__ = [
    "CDATVariable",
    "_Axis",
    "_TimeAxis",
    "_Grid",
    "_clean_attrs",
    "create_axis",
    "create_uniform_lat_axis",
    "create_uniform_lon_axis",
    "create_rect_grid",
    "create_variable",
    "da_to_cdat",
    "cdat_to_da",
    "validate_cdat_variable",
]


# ---------------------------------------------------------------------------
# Axis identification
# ---------------------------------------------------------------------------

_LAT_IDS = {"lat", "latitude", "j", "y", "Y", "yt_ocean", "yu_ocean"}
_LON_IDS = {"lon", "longitude", "i", "x", "X", "xt_ocean", "xu_ocean"}
_TIME_IDS = {"time", "t", "T"}
_LEV_IDS = {"lev", "level", "depth", "plev", "z", "Z", "st_ocean", "sw_ocean"}


def _detect_axis_type(ax_id: str) -> str:
    """
    Infer a CDAT-style axis code from a coordinate/dimension name alone.

    **Only used as a last resort by _dim_to_axis_type when CF metadata are
    absent.**  When ``STRICT_AXIS_DETECTION`` is True (the default) this
    function never returns ``"T"`` — time detection from names alone is too
    unreliable and has historically caused lat→T misclassification.
    """
    ax_id = str(ax_id)
    low = ax_id.lower()
    # Time: only returned when STRICT_AXIS_DETECTION is off, because
    # dimension names like 'lat_bnds', 'time_of_day', or single-letter
    # aliases can falsely match the 'time' substring check.
    if not STRICT_AXIS_DETECTION:
        if ax_id in _TIME_IDS or "time" in low:
            return "T"
    if ax_id in _LAT_IDS or "lat" in low:
        return "Y"
    if ax_id in _LON_IDS or "lon" in low:
        return "X"
    if ax_id in _LEV_IDS or "lev" in low or "depth" in low:
        return "Z"
    # Strict: time names require CF metadata; fall through to "-"
    if STRICT_AXIS_DETECTION:
        if ax_id in _TIME_IDS or "time" in low:
            import warnings
            warnings.warn(
                f"_detect_axis_type: dimension {ax_id!r} looks like a time axis "
                "by name but has no CF axis/standard_name/units metadata. "
                "Returning '-' to avoid silent misclassification. "
                "Add 'axis: T' or 'standard_name: time' metadata to the coordinate.",
                stacklevel=3,
            )
            return "-"
    return "-"


def _dim_to_axis_type(dim_name: str, coord) -> str:
    """Infer axis type, prioritising CF metadata over name heuristics."""
    # CF axis/standard_name/units are checked first so that dimensions like
    # "dim_0" or "y" with explicit metadata are classified correctly even when
    # the name alone would be ambiguous or unrecognised.
    if coord is not None:
        cf_axis = str(coord.attrs.get("axis", "")).upper()
        if cf_axis in {"T", "Y", "X", "Z"}:
            return cf_axis

        standard_name = str(coord.attrs.get("standard_name", "")).lower()
        if standard_name in {"latitude", "grid_latitude",
                             "projection_y_coordinate", "rotated_latitude"}:
            return "Y"
        if standard_name in {"longitude", "grid_longitude",
                             "projection_x_coordinate", "rotated_longitude"}:
            return "X"
        if standard_name == "time":
            return "T"
        if standard_name in {"air_pressure", "altitude", "depth", "height",
                             "ocean_sigma_coordinate", "sigma",
                             "height_above_geopotential_datum",
                             "height_above_mean_sea_level"}:
            return "Z"

        units = str(coord.attrs.get("units", "")).lower()
        if units in {"degrees_north", "degree_north", "degrees_n", "degree_n"}:
            return "Y"
        if units in {"degrees_east", "degree_east", "degrees_e", "degree_e"}:
            return "X"
        if "since" in units:
            return "T"
        if units in {"pa", "hpa", "mb", "mbar"}:
            return "Z"

        # Detect decoded datetime-like coordinates by value type —
        # cftime.datetime and np.datetime64 are unambiguously time axes even
        # when CF units/axis/standard_name metadata are absent (e.g. after
        # xESMF regridding strips coordinate attributes).
        try:
            first = coord.values.flat[0] if coord.size > 0 else None
            if first is not None and _is_datetime_like(first):
                return "T"
            if first is not None and np.issubdtype(type(first), np.datetime64):
                return "T"
        except Exception:
            pass

    # Last resort: name-based heuristic
    return _detect_axis_type(dim_name)


# ---------------------------------------------------------------------------
# Private time-handling helpers
# ---------------------------------------------------------------------------


def _is_datetime_like(obj) -> bool:
    """Duck-type check for datetime-like objects: cftime, datetime, Timestamp."""
    return (
        hasattr(obj, "year")
        and hasattr(obj, "month")
        and hasattr(obj, "day")
        and not isinstance(obj, (int, float, np.integer, np.floating))
    )


def _get_time_coder():
    """Return xarray's CFDatetimeCoder, handling API changes across versions."""
    try:
        return xr.coders.CFDatetimeCoder(use_cftime=True)
    except AttributeError:
        return xr.coding.times.CFDatetimeCoder(use_cftime=True)


def _decode_times_safe(arr: np.ndarray, units: str, calendar: str) -> list:
    """
    Decode numeric times to datetime-like objects through xarray's CF coder.

    If a decoder hits a leap-second-like value, retry by subtracting one second
    from the raw numeric offset. This keeps the behavior robust across common
    CF calendars used by CMIP, E3SM, ERA5-derived products, and observations.
    """
    coder = _get_time_coder()

    try:
        var = xr.Variable("time", arr, {"units": units, "calendar": calendar})
        decoded = coder.decode(var, name="time").values
        result = []
        for dt in decoded:
            if getattr(dt, "second", 0) == 60:
                try:
                    dt = dt.replace(second=59)
                except Exception:
                    pass
            result.append(dt)
        return result
    except Exception:
        pass

    result = []
    for v in np.asarray(arr, dtype=float):
        v_adj = float(v)
        dt = None
        for _ in range(2):
            try:
                var = xr.Variable(
                    "time",
                    np.array([v_adj], dtype=float),
                    {"units": units, "calendar": calendar},
                )
                decoded_v = coder.decode(var, name="time").values[0]
                if getattr(decoded_v, "second", 0) == 60:
                    v_adj -= 1.0 / 86400.0
                    continue
                dt = decoded_v
                break
            except ValueError:
                v_adj -= 1.0 / 86400.0
            except Exception:
                break

        if dt is None:
            # Time decoding failed for this value.  The 2000-01-01 sentinel
            # was removed because silently substituting a fake date causes
            # downstream seasonal averages, time-slicing, and climatologies to
            # produce wrong results without any visible error.
            #
            # We raise to surface the problem immediately so it can be fixed at
            # the source (bad time:units, unsupported calendar, leap-second
            # artefact, etc.).
            raise RuntimeError(
                f"Time decoding failed for numeric value {float(v)!r} "
                f"(units={units!r}, calendar={calendar!r}). "
                "Check for bad time coordinate values in the source file. "
                "Possible causes: unsupported calendar, invalid leap-second "
                "offset, or missing/incorrect time:units attribute."
            )
        result.append(dt)
    return result


def _encode_times(comp, units: str, calendar: str) -> np.ndarray:
    """Encode datetime-like values to numeric time using cftime lazily."""
    import cftime as _cft

    return np.asarray(_cft.date2num(comp, units, calendar=calendar), dtype=float)


def _parse_datetime_string(value: str):
    """Parse common ISO-like datetime/date strings to Python datetime."""
    text = value.strip().replace("T", " ")
    if text.endswith("Z"):
        text = text[:-1]
    text = text.split(".")[0]

    for fmt in (
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%d %H:%M",
        "%Y-%m-%d",
        "%Y/%m/%d %H:%M:%S",
        "%Y/%m/%d",
        "%Y-%m",
        "%Y",
    ):
        try:
            dt = _datetime.datetime.strptime(text, fmt)
            if fmt == "%Y-%m":
                dt = dt.replace(day=1)
            return dt
        except ValueError:
            continue

    # pandas/xarray environments often parse more formats, but keep pandas
    # optional by not importing it directly.
    try:
        return np.datetime64(text).astype("datetime64[us]").astype(_datetime.datetime)
    except Exception as exc:
        raise ValueError(f"Could not parse time bound {value!r}") from exc


def _coerce_time_bound(bound, ref):
    """Convert a time-selection bound to a comparable object using ref type."""
    if _is_datetime_like(bound):
        return bound

    if isinstance(bound, np.datetime64):
        bound = bound.astype("datetime64[us]").astype(_datetime.datetime)
        return bound

    if isinstance(bound, str):
        py_dt = _parse_datetime_string(bound)
        if _is_datetime_like(ref) and type(ref).__module__.startswith("cftime"):
            # Construct the same cftime calendar class as the axis values.
            try:
                return type(ref)(
                    py_dt.year,
                    py_dt.month,
                    py_dt.day,
                    py_dt.hour,
                    py_dt.minute,
                    py_dt.second,
                )
            except Exception:
                return py_dt
        return py_dt

    return bound


def _time_in_range(t, lo, hi) -> bool:
    """Robust inclusive datetime comparison with string fallback."""
    try:
        return lo <= t <= hi
    except TypeError:
        # Cross-calendar cftime comparisons may fail. ISO-style strings are a
        # fallback only after true datetime comparison has failed.
        return str(lo) <= str(t) <= str(hi)


# ---------------------------------------------------------------------------
# Metadata validation and array helpers
# ---------------------------------------------------------------------------


def _coerce_axes(axes, ndim: int) -> list:
    """
    Ensure every element in an axes list is an _Axis or None.

    Raw numpy arrays, such as from ``axis[:]``, are wrapped in a generic _Axis.
    """
    if axes is None:
        return []

    result = []
    for i, a in enumerate(axes):
        if a is None or isinstance(a, _Axis):
            result.append(a)
        elif isinstance(a, np.ndarray):
            result.append(_Axis(f"dim_{i}", a))
        else:
            try:
                result.append(_Axis(f"dim_{i}", np.asarray(a)))
            except Exception:
                result.append(None)
    return result


def _validate_axes_shape(data, axes, context: str = "CDATVariable"):
    """
    Minimal metadata safety check.

    It catches the most dangerous silent failure mode: data shape and coordinate
    axis metadata no longer match.
    """
    if axes is None or len(axes) == 0:
        return

    if len(axes) != data.ndim:
        raise ValueError(
            f"{context}: number of axes ({len(axes)}) does not match "
            f"data ndim ({data.ndim})"
        )

    for i, ax in enumerate(axes):
        if ax is None:
            continue
        if len(ax) != data.shape[i]:
            raise ValueError(
                f"{context}: axis {i} ({ax.id}) length {len(ax)} does not "
                f"match data shape {data.shape[i]}"
            )


def _maybe_mask_invalid_numeric(arr: ma.MaskedArray) -> ma.MaskedArray:
    """Mask NaN/Inf for floating-point arrays while preserving existing masks.

    Integer arrays are returned unchanged — NaN/Inf cannot appear in them and
    ``arr.filled(np.nan)`` raises a TypeError on integer dtypes.
    """
    if np.issubdtype(arr.dtype, np.integer):
        return arr
    if np.issubdtype(arr.dtype, np.number):
        with np.errstate(invalid="ignore"):
            invalid = ~np.isfinite(arr.filled(np.nan).astype(float))
        if invalid.shape == arr.shape:
            arr = ma.array(arr.data, mask=ma.getmaskarray(arr) | invalid, fill_value=arr.fill_value)
    return arr


def _make_masked_array(data, mask=None, fill_value=1e20, attributes: Optional[dict] = None):
    """Create a masked array while respecting common missing metadata.

    Integer-typed data is kept as integer; float fill_value (1e20) is only
    applied to floating-point arrays to avoid the OverflowError that numpy
    raises when trying to store 1e20 into an int64 fill_value slot.
    """
    attributes = attributes or {}

    if isinstance(data, CDATVariable):
        raw = data._data.copy()
    elif isinstance(data, ma.MaskedArray):
        raw = data.copy()
    else:
        arr = np.asarray(data)
        if np.issubdtype(arr.dtype, np.integer):
            raw = ma.array(arr)           # keep dtype; no float fill_value
        else:
            raw = ma.array(arr.astype(float), fill_value=fill_value)

    # Resolve an appropriate fill_value for the actual dtype of raw.
    if np.issubdtype(raw.dtype, np.integer):
        _fv = int(np.iinfo(raw.dtype).max)
    else:
        _fv = fill_value

    combined_mask = ma.getmaskarray(raw)

    for key in ("_FillValue", "missing_value"):
        if key in attributes:
            mv = attributes[key]
            try:
                if np.ndim(mv) == 0:
                    combined_mask |= np.asarray(raw.data == mv)
                else:
                    for one_mv in np.ravel(mv):
                        combined_mask |= np.asarray(raw.data == one_mv)
            except Exception:
                pass

    if mask is not None:
        combined_mask |= np.asarray(mask, dtype=bool)

    raw = ma.array(raw.data, mask=combined_mask, fill_value=_fv)
    raw = _maybe_mask_invalid_numeric(raw)
    return raw


def _clean_attrs(attrs: dict) -> dict:
    """
    Remove attributes whose values are None or otherwise not serialisable to
    netCDF (xarray rejects None, and also non-scalar non-array objects).

    Valid netCDF attribute types: str, Number, ndarray, list, tuple, bytes.
    """
    import numbers
    keep = (str, numbers.Number, np.ndarray, list, tuple, bytes)
    return {k: v for k, v in attrs.items() if v is not None and isinstance(v, keep)}


def _validate_grid(grid, axes, context: str = "CDATVariable"):
    """Light check that a rectilinear grid is compatible with available axes."""
    if grid is None or not axes:
        return
    lat = grid.getLatitude()
    lon = grid.getLongitude()
    lat_axes = [ax for ax in axes if ax is not None and ax.isLatitude()]
    lon_axes = [ax for ax in axes if ax is not None and ax.isLongitude()]
    if lat_axes and lat is not None and len(lat) != len(lat_axes[0]):
        raise ValueError(f"{context}: grid latitude length does not match latitude axis")
    if lon_axes and lon is not None and len(lon) != len(lon_axes[0]):
        raise ValueError(f"{context}: grid longitude length does not match longitude axis")
    # Guard against the degenerate case where the same axis object was used for
    # both lat and lon (e.g. a typo in the calling code).
    if lat_axes and lon_axes and lat_axes[0].id == lon_axes[0].id:
        raise ValueError(
            f"{context}: latitude and longitude axes have the same id "
            f"{lat_axes[0].id!r} — likely a metadata error."
        )
    # Warn when axes appear in (lon, lat) order rather than the standard CDAT
    # (lat, lon) / yx convention.  This does not raise because xy ordering is
    # valid in some downstream tools, but it is a common source of silent bugs.
    if lat_axes and lon_axes and axes:
        try:
            lat_idx = next(i for i, ax in enumerate(axes) if ax is lat_axes[0])
            lon_idx = next(i for i, ax in enumerate(axes) if ax is lon_axes[0])
            if lon_idx < lat_idx:
                import warnings
                msg = (
                    f"{context}: longitude axis ({lon_axes[0].id!r}, dim {lon_idx}) "
                    f"precedes latitude axis ({lat_axes[0].id!r}, dim {lat_idx}). "
                    "Standard CDAT convention is (lat, lon) / yx order."
                )
                if STRICT_GRID:
                    raise ValueError(msg)
                warnings.warn(msg, stacklevel=2)
        except StopIteration:
            pass  # axis not found in list (e.g. grid built from external axes)


def _build_grid_from_axes(axes) -> Optional["_Grid"]:
    """Build a rectilinear grid from available latitude/longitude axes."""
    if not axes:
        return None
    # Use ax.axis directly — _Axis.__init__ already resolves it from CF
    # metadata (axis attr, standard_name, units) and name heuristics.  Calling
    # _detect_axis_type(ax.id) again is redundant and emits spurious warnings
    # under STRICT_AXIS_DETECTION=True when the name looks like 'time'.
    lat_ax = next((ax for ax in axes if ax is not None and ax.axis == "Y"), None)
    lon_ax = next((ax for ax in axes if ax is not None and ax.axis == "X"), None)
    return _Grid(lat_ax, lon_ax) if (lat_ax is not None and lon_ax is not None) else None


def validate_cdat_variable(var, context: str = "CDATVariable"):
    """
    Validate that a CDATVariable carries axis metadata consistent with data shape.

    This catches the common failure mode that caused downstream ENSO metric
    crashes: data have been reduced/sliced while stale or missing axis metadata
    were propagated.
    """
    if not isinstance(var, CDATVariable):
        raise TypeError(f"{context}: expected CDATVariable, got {type(var)!r}")
    axes = var.getAxisList()
    _validate_axes_shape(var._data, axes, context=context)
    _validate_grid(var.getGrid(), axes, context=context)
    return var


def _coord_attrs_with_axis(coord, axis_type: str) -> dict:
    """Return coordinate attrs with CF axis metadata strengthened."""
    attrs = dict(getattr(coord, "attrs", {}) or {})
    if axis_type in {"T", "Y", "X", "Z"}:
        attrs.setdefault("axis", axis_type)
    if axis_type == "Y":
        attrs.setdefault("standard_name", "latitude")
        attrs.setdefault("units", "degrees_north")
    elif axis_type == "X":
        attrs.setdefault("standard_name", "longitude")
        attrs.setdefault("units", "degrees_east")
    elif axis_type == "T":
        attrs.setdefault("standard_name", "time")
        cal = attrs.get("calendar") or getattr(coord, "encoding", {}).get("calendar", None)
        if cal is not None:
            attrs.setdefault("calendar", cal)
    return attrs


# ---------------------------------------------------------------------------
# Axis
# ---------------------------------------------------------------------------


class _Axis:
    """Lightweight replacement for a cdms2 Axis object."""

    def __init__(
        self,
        id: str,
        values,
        *,
        units: str = "",
        attributes: Optional[dict] = None,
        axis_type: Optional[str] = None,
    ):
        self.id = id
        if values is None:
            self._values = np.array([])
        elif isinstance(values, np.ndarray):
            self._values = values
        else:
            try:
                self._values = np.asarray(values)
            except Exception:
                self._values = np.array(list(values), dtype=object)

        self.units = units
        self._attributes = dict(attributes or {})
        # Priority: explicit arg > CF attrs in _attributes > name heuristic.
        # This ensures _Axis objects built from xarray coords with standard CF
        # metadata (e.g. axis='Y', standard_name='latitude', units='degrees_north')
        # are typed correctly even when the dimension name is opaque (e.g. 'dim_0').
        if axis_type and axis_type != "-":
            self.axis = axis_type
        else:
            _cf = str(self._attributes.get("axis", "")).upper()
            if _cf in {"T", "Y", "X", "Z"}:
                self.axis = _cf
            else:
                _sn = str(self._attributes.get("standard_name", "")).lower()
                _u  = str(self._attributes.get("units", self.units or "")).lower()
                if _sn in {"latitude", "grid_latitude",
                           "projection_y_coordinate", "rotated_latitude"} \
                        or _u in {"degrees_north", "degree_north",
                                  "degrees_n", "degree_n"}:
                    self.axis = "Y"
                elif _sn in {"longitude", "grid_longitude",
                             "projection_x_coordinate", "rotated_longitude"} \
                        or _u in {"degrees_east", "degree_east",
                                  "degrees_e", "degree_e"}:
                    self.axis = "X"
                elif _sn == "time" or "since" in _u:
                    self.axis = "T"
                elif _sn in {"air_pressure", "altitude", "depth", "height",
                             "ocean_sigma_coordinate", "sigma",
                             "height_above_geopotential_datum",
                             "height_above_mean_sea_level"} \
                        or _u in {"pa", "hpa", "mb", "mbar"}:
                    self.axis = "Z"
                else:
                    self.axis = _detect_axis_type(id)
        self.long_name = self._attributes.get("long_name", id)
        self.regions: Optional[str] = None
        self.reference: Optional[str] = None
        self.calendar: Optional[str] = self._attributes.get("calendar", None)

    def isTime(self) -> bool:
        return self.axis == "T"

    def isLatitude(self) -> bool:
        return self.axis == "Y"

    def isLongitude(self) -> bool:
        return self.axis == "X"

    def isLevel(self) -> bool:
        return self.axis == "Z"

    @property
    def shape(self):
        return self._values.shape

    def __len__(self):
        return len(self._values)

    def __getitem__(self, key):
        return self._values[key]

    def __setitem__(self, key, value):
        self._values[key] = value

    def __array__(self, dtype=None):
        return np.array(self._values, dtype=dtype)

    def __repr__(self):
        return f"_Axis(id={self.id!r}, axis={self.axis!r}, len={len(self._values)})"

    def asComponentTime(self) -> list:
        """Return list of datetime-like objects, mirroring cdtime behavior."""
        vals = self._values
        if len(vals) == 0:
            return []

        # Handle numpy.datetime64 (from xarray-decoded time axes without cftime).
        # These fail _is_datetime_like (no .year attr) and cannot be cast to float
        # meaningfully, so they must be converted before the generic path.
        if isinstance(vals[0], np.datetime64):
            result = []
            for v in vals:
                try:
                    # .item() returns datetime.datetime for sub-day precision,
                    # datetime.date for day-only — both satisfy downstream year/month.
                    dt = v.item()
                    if not isinstance(dt, _datetime.datetime):
                        dt = _datetime.datetime(dt.year, dt.month, dt.day)
                except Exception:
                    try:
                        secs = int(v.astype("datetime64[s]").astype(np.int64))
                        dt = _datetime.datetime(1970, 1, 1) + _datetime.timedelta(seconds=secs)
                    except Exception:
                        dt = _datetime.datetime(2000, 1, 1)
                if getattr(dt, "second", 0) == 60:
                    try:
                        dt = dt.replace(second=59)
                    except Exception:
                        pass
                result.append(dt)
            return result

        if _is_datetime_like(vals[0]):
            result = []
            for dt in vals:
                if getattr(dt, "second", 0) == 60:
                    try:
                        dt = dt.replace(second=59)
                    except Exception:
                        pass
                result.append(dt)
            return result

        if self.units:
            try:
                cal = self.calendar or "standard"
                arr = np.asarray(vals, dtype=float)
                return _decode_times_safe(arr, self.units, cal)
            except Exception:
                pass

        # Year 0 is invalid in Python's datetime; clamp to 1.
        result = []
        for v in vals:
            yr = max(1, int(v))
            try:
                result.append(_datetime.datetime(yr, 1, 1))
            except Exception:
                result.append(_datetime.datetime(1, 1, 1))
        return result

    def toRelativeTime(self, units: str):
        """Convert datetime-like values to numeric relative time in-place."""
        try:
            cal = self.calendar or "standard"
            comp = self.asComponentTime()
            self._values = _encode_times(comp, units, cal)
            self.units = units
        except Exception:
            pass

    def copy(self) -> "_Axis":
        new = _Axis(
            self.id,
            self._values.copy(),
            units=self.units,
            attributes=dict(self._attributes),
            axis_type=self.axis,
        )
        new.regions = self.regions
        new.reference = self.reference
        new.calendar = self.calendar
        return new


class _TimeAxis(_Axis):
    """Time-specific axis, always typed 'T'."""

    def __init__(self, id: str = "time", values=None, **kwargs):
        kwargs.setdefault("axis_type", "T")
        super().__init__(id, values, **kwargs)


# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------


class _Grid:
    """Lightweight replacement for cdms2.RectGrid."""

    def __init__(self, lat_axis: _Axis, lon_axis: _Axis, id: str = ""):
        self._lat = lat_axis
        self._lon = lon_axis
        self.id = id
        self.type = "generic"

    @property
    def shape(self):
        return (len(self._lat), len(self._lon))

    def getLatitude(self) -> _Axis:
        return self._lat

    def getLongitude(self) -> _Axis:
        return self._lon

    def __repr__(self):
        return f"_Grid(shape={self.shape}, id={self.id!r})"


# ---------------------------------------------------------------------------
# CDATVariable
# ---------------------------------------------------------------------------


class CDATVariable:
    """
    CDAT-like variable object used by the ENSO_metrics compatibility layer.

    This class preserves the subset of the legacy ``cdms2.TransientVariable``
    interface needed by ENSO_metrics while storing data internally as
    ``numpy.ma.MaskedArray``. Coordinate metadata are carried as a list of
    ``_Axis`` objects plus an optional rectilinear ``_Grid``.

    It is intended for decoded rectilinear gridded fields and should not be
    treated as a full replacement for the complete CDAT/UV-CDAT variable API.
    """
    def __init__(
        self,
        data,
        axes: Optional[List[_Axis]] = None,
        grid: Optional[_Grid] = None,
        mask=None,
        id: str = "",
        attributes: Optional[dict] = None,
        fill_value=1e20,
    ):
        self._attributes: dict = dict(attributes or {})
        self._data: ma.MaskedArray = _make_masked_array(
            data,
            mask=mask,
            fill_value=fill_value,
            attributes=self._attributes,
        )

        self._axes: List[_Axis] = _coerce_axes(axes, self._data.ndim)
        _validate_axes_shape(self._data, self._axes, context=f"CDATVariable({id})")
        _validate_grid(grid, self._axes, context=f"CDATVariable({id})")

        self._grid: Optional[_Grid] = grid
        self.id: str = id
        self.name: str = id

    @property
    def shape(self):
        return self._data.shape

    @property
    def ndim(self):
        return self._data.ndim

    @property
    def size(self):
        return self._data.size

    @property
    def dtype(self):
        return self._data.dtype

    @property
    def mask(self):
        return ma.getmaskarray(self._data)

    @property
    def units(self) -> str:
        return self._attributes.get("units", "")

    @units.setter
    def units(self, value: str):
        self._attributes["units"] = value

    @property
    def attributes(self) -> dict:
        return self._attributes

    @attributes.setter
    def attributes(self, value: dict):
        self._attributes = dict(value)

    def __array__(self, dtype=None):
        return np.array(self._data, dtype=dtype)

    def __float__(self):
        return float(self._data)

    def __int__(self):
        return int(self._data)

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        for i in range(len(self._data)):
            yield self._wrap(self._data[i], axes=self._axes[1:])

    def __repr__(self):
        return f"CDATVariable(id={self.id!r}, shape={self.shape}, dtype={self.dtype})"

    def __str__(self):
        return str(self._data)

    # ------------------------------------------------------------------
    # Indexing / slicing
    # ------------------------------------------------------------------
    def __getitem__(self, key):
        result = self._data[key]
        if not isinstance(result, (np.ndarray, ma.MaskedArray)):
            return result
        new_axes = self._sliced_axes(key, result.shape)
        # Don't propagate a stale grid when axes were dropped (advanced indexing
        # fallback returns [] when axis mapping cannot be determined).
        grid = self._grid if new_axes else None
        return CDATVariable(
            result,
            axes=new_axes,
            grid=grid,
            id=self.id,
            attributes=dict(self._attributes),
        )

    def __setitem__(self, key, value):
        if isinstance(value, CDATVariable):
            self._data[key] = value._data
        else:
            self._data[key] = value

    def _sliced_axes(self, key, new_shape):
        """Best-effort axis update for NumPy-style slicing."""
        if not self._axes:
            return []

        if not isinstance(key, tuple):
            key = (key,)

        # Expand ellipsis and append missing full slices.
        # Use `any(k is Ellipsis ...)` to avoid the ambiguous truth value of
        # numpy arrays that may appear as advanced-indexing keys.
        key_list = list(key)
        if any(k is Ellipsis for k in key_list):
            ell_idx = next(i for i, k in enumerate(key_list) if k is Ellipsis)
            n_missing = self._data.ndim - (len(key_list) - 1)
            key_list = key_list[:ell_idx] + [slice(None)] * n_missing + key_list[ell_idx + 1:]
        if len(key_list) < self._data.ndim:
            key_list += [slice(None)] * (self._data.ndim - len(key_list))

        new_axes = []
        ax_idx = 0
        for k in key_list:
            if k is None:
                # np.newaxis has no original coordinate axis.
                new_axes.append(None)
                continue
            if ax_idx >= len(self._axes):
                break

            ax = self._axes[ax_idx]
            ax_idx += 1

            if isinstance(k, (int, np.integer)):
                continue

            if ax is None:
                new_axes.append(None)
                continue

            try:
                new_vals = ax._values[k]
            except Exception:
                # Advanced indexing can be complicated; preserve the axis only
                # if the length still matches after the resulting data check.
                new_vals = ax._values

            new_axes.append(
                _Axis(
                    ax.id,
                    new_vals,
                    units=ax.units,
                    attributes=dict(ax._attributes),
                    axis_type=ax.axis,
                )
            )

        # If advanced indexing changed the number of dimensions we cannot
        # reliably map old axes to new ones — drop all metadata.
        if len(new_axes) != len(new_shape):
            import warnings
            # Identify which CF axis types were lost so the warning is actionable.
            lost = [ax.axis for ax in self._axes if ax is not None and ax.axis in ("T", "Y", "X", "Z")]
            warnings.warn(
                f"CDATVariable[{self.id!r}]: advanced indexing changed number of "
                f"dimensions; axis metadata dropped to avoid mis-labelling. "
                f"Axes that were present before slicing: {lost}.",
                stacklevel=3,
            )
            return []
        # Per-axis length mismatch: replace only the broken axis with None so
        # that the remaining (still-valid) coordinate axes are preserved.
        _important_axes = {"T", "Y", "X"}
        for i, (ax, n) in enumerate(zip(new_axes, new_shape)):
            if ax is not None and len(ax) != n:
                import warnings
                if ax.axis in _important_axes:
                    warnings.warn(
                        f"CDATVariable[{self.id!r}]: slicing broke axis "
                        f"{ax.id!r} (type {ax.axis!r}, expected {n} values, "
                        f"got {len(ax)}); axis metadata dropped for that dimension. "
                        "Downstream calls to getTime()/getLatitude()/getLongitude() "
                        "may fail.",
                        stacklevel=3,
                    )
                else:
                    warnings.warn(
                        f"CDATVariable[{self.id!r}]: advanced indexing broke axis "
                        f"length for {ax.id!r} (expected {n}, got {len(ax)}); "
                        "axis metadata for that dimension dropped.",
                        stacklevel=3,
                    )
                new_axes[i] = None
        return new_axes

    # ------------------------------------------------------------------
    # Arithmetic operators
    # ------------------------------------------------------------------
    def _unpack(self, other):
        if isinstance(other, CDATVariable):
            return other._data
        return other

    def _binary_axes(self, result, other):
        """
        Preserve axes when the result has the same number of dimensions as self.

        * Exact shape match (common case, including scalar ops) → copy axes as-is.
        * Same ndim but a dimension was broadcast-expanded → keep axis metadata
          where the length still matches; replace broadcast-expanded dims with None
          rather than propagating wrong coordinate values.
        * ndim changed (e.g. outer product) → drop all axes.
        """
        if not isinstance(result, (np.ndarray, ma.MaskedArray)):
            return []
        if result.ndim != self.ndim:
            return []
        if tuple(result.shape) == tuple(self.shape):
            return [ax.copy() if ax is not None else None for ax in self._axes]
        # Same ndim, different shape: some dimension was broadcast-expanded.
        import warnings
        new_axes = []
        for ax, n in zip(self._axes, result.shape):
            if ax is None:
                new_axes.append(None)
            elif len(ax) == n:
                new_axes.append(ax.copy())
            else:
                # Broadcast expanded this dim — coordinate values no longer
                # correspond to the result shape; drop axis and warn.
                warnings.warn(
                    f"CDATVariable[{self.id!r}]: dimension {ax.id!r} was "
                    f"broadcast-expanded ({len(ax)} → {n}); "
                    "axis metadata dropped for that dimension.",
                    stacklevel=3,
                )
                new_axes.append(None)
        return new_axes

    def _wrap(self, data, axes=None) -> "CDATVariable":
        if not isinstance(data, (np.ndarray, ma.MaskedArray)):
            return data
        axes = self._axes if axes is None else axes
        axes = [ax.copy() if ax is not None else None for ax in axes]
        return CDATVariable(
            data,
            axes=axes,
            grid=_build_grid_from_axes(axes),
            id=self.id,
            attributes=dict(self._attributes),
        )

    def _wrap_binary(self, result, other):
        return self._wrap(result, axes=self._binary_axes(result, other))

    def __add__(self, other):       return self._wrap_binary(self._data + self._unpack(other), other)
    def __radd__(self, other):      return self._wrap_binary(self._unpack(other) + self._data, other)
    def __sub__(self, other):       return self._wrap_binary(self._data - self._unpack(other), other)
    def __rsub__(self, other):      return self._wrap_binary(self._unpack(other) - self._data, other)
    def __mul__(self, other):       return self._wrap_binary(self._data * self._unpack(other), other)
    def __rmul__(self, other):      return self._wrap_binary(self._unpack(other) * self._data, other)
    def __truediv__(self, other):   return self._wrap_binary(self._data / self._unpack(other), other)
    def __rtruediv__(self, other):  return self._wrap_binary(self._unpack(other) / self._data, other)
    def __neg__(self):              return self._wrap(-self._data)
    def __abs__(self):              return self._wrap(abs(self._data))
    def __pow__(self, exp):         return self._wrap(self._data ** exp)

    def __gt__(self, other): return self._data > self._unpack(other)
    def __lt__(self, other): return self._data < self._unpack(other)
    def __ge__(self, other): return self._data >= self._unpack(other)
    def __le__(self, other): return self._data <= self._unpack(other)
    def __eq__(self, other): return self._data == self._unpack(other)
    def __ne__(self, other): return self._data != self._unpack(other)

    # ------------------------------------------------------------------
    # numpy.ma delegation
    # ------------------------------------------------------------------
    def filled(self, fill_value=1e20):
        return self._data.filled(fill_value)

    def argmin(self, axis=None, fill_value=None, out=None):
        return self._data.argmin(axis=axis)

    def argmax(self, axis=None, fill_value=None, out=None):
        return self._data.argmax(axis=axis)

    def mean(self, axis=None):
        return self._data.mean(axis=axis)

    def sum(self, axis=None):
        return self._data.sum(axis=axis)

    def std(self, axis=None, ddof=0):
        return self._data.std(axis=axis, ddof=ddof)

    def min(self, axis=None):
        return self._data.min(axis=axis)

    def max(self, axis=None):
        return self._data.max(axis=axis)

    def fill(self, value):
        self._data.fill(value)

    def astype(self, dtype) -> "CDATVariable":
        return self._wrap(self._data.astype(dtype))

    def squeeze(self, axis=None) -> "CDATVariable":
        result = self._data.squeeze(axis=axis)
        if axis is None:
            # Drop exactly the singleton dimensions removed by numpy.squeeze.
            new_axes = [
                ax.copy() if ax is not None else None
                for i, ax in enumerate(self._axes)
                if self._data.shape[i] != 1
            ]
        else:
            # Accept int, numpy.integer, or tuple/list of ints. numpy.squeeze
            # already raises if a requested axis is not singleton.
            axes_to_drop = {int(axis)} if isinstance(axis, (int, np.integer)) else {int(a) for a in axis}
            axes_to_drop = {a if a >= 0 else self.ndim + a for a in axes_to_drop}
            new_axes = [
                ax.copy() if ax is not None else None
                for i, ax in enumerate(self._axes)
                if i not in axes_to_drop
            ]
        return CDATVariable(
            result,
            axes=new_axes,
            grid=_build_grid_from_axes(new_axes),
            id=self.id,
            attributes=dict(self._attributes),
        )

    def compress(self, condition, axis: int = 0) -> "CDATVariable":
        result = self._data.compress(condition, axis=axis)
        new_axes = [ax.copy() if ax is not None else None for ax in self._axes]
        axis = axis if axis >= 0 else self.ndim + axis
        if self._axes and axis < len(self._axes) and self._axes[axis] is not None:
            old_ax = self._axes[axis]
            new_vals = old_ax._values[np.asarray(condition, dtype=bool)]
            new_axes[axis] = _Axis(old_ax.id, new_vals, units=old_ax.units, attributes=dict(old_ax._attributes), axis_type=old_ax.axis)
        return CDATVariable(result, axes=new_axes, grid=_build_grid_from_axes(new_axes), id=self.id, attributes=dict(self._attributes))

    def copy(self) -> "CDATVariable":
        new_axes = [ax.copy() if ax is not None else None for ax in self._axes]
        return CDATVariable(
            self._data.copy(),
            axes=new_axes,
            # Rebuild the grid from the *copied* axes rather than carrying a
            # reference to the original grid (which points to old axis objects).
            grid=_build_grid_from_axes(new_axes),
            id=self.id,
            attributes=dict(self._attributes),
        )

    # ------------------------------------------------------------------
    # CDAT-compatible introspection methods
    # ------------------------------------------------------------------
    def getAxisList(self) -> List[_Axis]:
        return list(self._axes)

    def getAxis(self, n: int) -> Optional[_Axis]:
        if self._axes and 0 <= n < len(self._axes):
            return self._axes[n]
        return None

    def setAxis(self, n: int, ax: _Axis):
        while len(self._axes) <= n:
            self._axes.append(None)
        self._axes[n] = ax
        _validate_axes_shape(self._data, self._axes, context=f"CDATVariable({self.id}).setAxis")
        # Rebuild the rectilinear grid so getGrid() is never stale after a
        # setAxis call (e.g. after toRelativeTime modifies the time axis).
        self._grid = _build_grid_from_axes(self._axes)

    def setAxisList(self, axes: list):
        axes = _coerce_axes(axes, self._data.ndim)
        _validate_axes_shape(self._data, axes, context=f"CDATVariable({self.id}).setAxisList")
        self._axes = axes
        # Rebuild the rectilinear grid so getGrid() stays in sync after a
        # bulk axis replacement (mirrors the grid-rebuild in setAxis).
        self._grid = _build_grid_from_axes(self._axes)

    def getGrid(self) -> Optional[_Grid]:
        return self._grid

    def setGrid(self, grid: Optional[_Grid]):
        _validate_grid(grid, self._axes, context=f"CDATVariable({self.id}).setGrid")
        self._grid = grid

    def getTime(self) -> Optional[_Axis]:
        # ax.axis is resolved by _Axis.__init__ — no name-heuristic fallback needed.
        return next((ax for ax in self._axes if ax is not None and ax.axis == "T"), None)

    def getLatitude(self) -> Optional[_Axis]:
        return next((ax for ax in self._axes if ax is not None and ax.axis == "Y"), None)

    def getLongitude(self) -> Optional[_Axis]:
        return next((ax for ax in self._axes if ax is not None and ax.axis == "X"), None)

    def getLevel(self) -> Optional[_Axis]:
        return next((ax for ax in self._axes if ax is not None and ax.axis == "Z"), None)

    def getOrder(self) -> str:
        order = ""
        for ax in self._axes:
            if ax is None:
                order += "-"
            elif ax.axis == "T":
                order += "t"
            elif ax.axis == "Y":
                order += "y"
            elif ax.axis == "X":
                order += "x"
            elif ax.axis == "Z":
                order += "z"
            else:
                order += "-"
        return order

    def reorder(self, order: str) -> "CDATVariable":
        """
        Reorder axes.

        Accepts examples like 'tyx', 'txy', '10', '210', 't...', and '...t'.
        Invalid permutations raise ValueError instead of silently returning an
        unchanged variable.
        """
        ndim = self._data.ndim
        if ndim == 0:
            return self.copy()

        def _ax_type(ax):
            """Return the single-letter axis type; ax.axis is already resolved by __init__."""
            return ax.axis

        if order in ("t...", "T..."):
            t_n = next((i for i, ax in enumerate(self._axes) if ax is not None and _ax_type(ax) == "T"), 0)
            perm = [t_n] + [i for i in range(ndim) if i != t_n]
        elif order in ("...t", "...T"):
            t_n = next((i for i, ax in enumerate(self._axes) if ax is not None and _ax_type(ax) == "T"), ndim - 1)
            perm = [i for i in range(ndim) if i != t_n] + [t_n]
        elif order and all(c.isdigit() for c in order):
            perm = [int(c) for c in order]
            # Pad with any axes not yet listed, preserving their original order.
            # This matches CDAT behaviour: reorder("10") on a 4-D array → [1,0,2,3].
            used = set(perm)
            perm = perm + [i for i in range(ndim) if i not in used]
        else:
            char_map = {}
            for i, ax in enumerate(self._axes):
                if ax is None:
                    continue
                t = _ax_type(ax)
                if t == "T":
                    char_map["t"] = i
                elif t == "Y":
                    char_map["y"] = i
                elif t == "X":
                    char_map["x"] = i
                elif t == "Z":
                    char_map["z"] = i

            perm = []
            used = set()
            for c in order.lower():
                if c == ".":
                    continue
                if c in char_map and char_map[c] not in used:
                    perm.append(char_map[c])
                    used.add(char_map[c])
            for i in range(ndim):
                if i not in used:
                    perm.append(i)
                    used.add(i)

        if len(perm) != ndim or sorted(perm) != list(range(ndim)):
            raise ValueError(
                f"Invalid reorder specification {order!r}: resolved permutation "
                f"{perm!r} is not valid for {ndim} dimensions"
            )

        new_data = np.transpose(self._data, perm)
        # Guard: _axes may be empty (e.g. from MV2zeros) — build stub axes in
        # that case so reorder doesn't raise IndexError.
        safe_axes = self._axes if len(self._axes) == self._data.ndim else [None] * self._data.ndim
        new_axes = [safe_axes[i].copy() if safe_axes[i] is not None else None for i in perm]
        return CDATVariable(new_data, axes=new_axes, grid=_build_grid_from_axes(new_axes), id=self.id, attributes=dict(self._attributes))

    # ------------------------------------------------------------------
    # CDAT-style callable selection: var(time=..., latitude=..., longitude=...)
    # ------------------------------------------------------------------
    def __call__(self, *args, **kwargs) -> "CDATVariable":
        if not kwargs:
            return self.copy()

        result_data = self._data.copy()
        result_axes = [ax.copy() if ax is not None else None for ax in self._axes]

        if kwargs.get("squeeze"):
            result_data = result_data.squeeze()
            result_axes = [ax for ax in result_axes if ax is not None and len(ax) > 1]
            return CDATVariable(
                result_data, axes=result_axes,
                # Rebuild grid from remaining axes; self._grid may reference a
                # level axis that was just squeezed out.
                grid=_build_grid_from_axes(result_axes),
                id=self.id, attributes=dict(self._attributes))

        def _sel_axis(ax_idx, bounds):
            nonlocal result_data
            ax = result_axes[ax_idx]
            if ax is None:
                return

            if not isinstance(bounds, (tuple, list)) or len(bounds) != 2:
                raise ValueError("Selection bounds must be a 2-element tuple/list")

            if ax.axis == "T":
                t_vals = ax.asComponentTime()
                if not t_vals:
                    indices = []
                else:
                    t0 = _coerce_time_bound(bounds[0], t_vals[0])
                    t1 = _coerce_time_bound(bounds[1], t_vals[0])
                    try:
                        lo, hi = (t0, t1) if t0 <= t1 else (t1, t0)
                    except TypeError:
                        lo, hi = (t0, t1) if str(t0) <= str(t1) else (t1, t0)
                    indices = [i for i, t in enumerate(t_vals) if _time_in_range(t, lo, hi)]
            else:
                raw = ax._values.astype(float)
                lo, hi = min(float(bounds[0]), float(bounds[1])), max(float(bounds[0]), float(bounds[1]))
                indices = list(np.where((raw >= lo) & (raw <= hi))[0])

            slices = [slice(None)] * result_data.ndim
            slices[ax_idx] = indices
            result_data = result_data[tuple(slices)]

            old_ax = ax
            new_vals = old_ax._values[indices]
            result_axes[ax_idx] = _Axis(
                old_ax.id,
                new_vals,
                units=old_ax.units,
                attributes=dict(old_ax._attributes),
                axis_type=old_ax.axis,
            )

        if "time" in kwargs:
            t_idx = next((i for i, ax in enumerate(result_axes)
                          if ax is not None and ax.axis == "T"), None)
            # Fallback: if no axis is typed "T", look for any axis whose values
            # look like datetime objects (handles axes with axis=="-" when CF
            # metadata was absent during construction).
            if t_idx is None:
                for i, ax in enumerate(result_axes):
                    if ax is None:
                        continue
                    vals = ax._values
                    if len(vals) > 0 and (_is_datetime_like(vals[0])
                                          or isinstance(vals[0], np.datetime64)
                                          or "since" in str(ax.units).lower()):
                        t_idx = i
                        ax.axis = "T"  # re-tag so getTime() works afterwards
                        break
            if t_idx is not None:
                _sel_axis(t_idx, kwargs["time"])

        if "latitude" in kwargs:
            lat_idx = next((i for i, ax in enumerate(result_axes)
                            if ax is not None and ax.axis == "Y"), None)
            if lat_idx is None:
                raise ValueError(
                    f"CDATVariable.__call__: latitude= selection requested but no axis "
                    f"tagged as 'Y' was found in variable '{self.id}'. "
                    f"Axes: {[ax.id if ax is not None else None for ax in self._axes]}"
                )
            _sel_axis(lat_idx, kwargs["latitude"])

        if "longitude" in kwargs:
            lon_idx = next((i for i, ax in enumerate(result_axes)
                            if ax is not None and ax.axis == "X"), None)
            if lon_idx is None:
                raise ValueError(
                    f"CDATVariable.__call__: longitude= selection requested but no axis "
                    f"tagged as 'X' was found in variable '{self.id}'. "
                    f"Axes: {[ax.id if ax is not None else None for ax in self._axes]}"
                )
            _sel_axis(lon_idx, kwargs["longitude"])

        lat_ax = next((ax for ax in result_axes if ax is not None and ax.axis == "Y"), None)
        lon_ax = next((ax for ax in result_axes if ax is not None and ax.axis == "X"), None)
        new_grid = _Grid(lat_ax, lon_ax) if (lat_ax and lon_ax) else None

        return CDATVariable(result_data, axes=result_axes, grid=new_grid, id=self.id, attributes=dict(self._attributes))


# ---------------------------------------------------------------------------
# Factory helpers
# ---------------------------------------------------------------------------


def create_axis(values, id: str = "", units: str = "", 
                attributes: Optional[dict] = None,
                axis_type: Optional[str] = None) -> _Axis:
    """Create a CDAT-like axis object for the compatibility layer.

    This helper preserves the subset of the legacy ``cdms2.createAxis`` calling
    convention used by ENSO_metrics: ``values`` is the first positional
    argument and ``id`` is an optional keyword argument.

    Parameters
    ----------
    values : array-like
        Axis coordinate values.

    id : str, optional
        Axis identifier, for example ``"time"``, ``"lat"``, ``"lon"``, or
        ``"lev"``.

    units : str, optional
        Axis units. For CF-style spatial axes, typical values are
        ``"degrees_north"`` for latitude and ``"degrees_east"`` for longitude.

    attributes : dict, optional
        Additional metadata attributes to attach to the axis.

    axis_type : str, optional
        Explicit CF axis code, one of ``"T"``, ``"Y"``, ``"X"``, or ``"Z"``.
        Pass this when the axis id or units alone would not be enough to detect
        the type under ``STRICT_AXIS_DETECTION=True``; for example, a synthetic
        integer time axis created with ``id="time"`` but no CF time units yet.

    Returns
    -------
    _Axis
        CDAT-like axis object used by ``CDATVariable``.
    """
    return _Axis(id, values, units=units, attributes=attributes, axis_type=axis_type)


def create_uniform_lat_axis(start: float, n: int, delta: float) -> _Axis:
    """Create a uniformly spaced CDAT-like latitude axis.

    This helper preserves the subset of the legacy
    ``cdms2.createUniformLatitudeAxis`` behavior used by ENSO_metrics while
    returning an ``_Axis`` object for the xarray compatibility layer.

    Parameters
    ----------
    start : float
        First latitude coordinate value.

    n : int
        Number of latitude points.

    delta : float
        Latitude spacing in degrees.

    Returns
    -------
    _Axis
        Latitude axis with CF-style metadata, including ``axis="Y"``,
        ``standard_name="latitude"``, and ``units="degrees_north"``.
    """
    vals = np.array([start + i * delta for i in range(int(n))])
    return _Axis("lat", vals, units="degrees_north", axis_type="Y")


def create_uniform_lon_axis(start: float, n: int, delta: float) -> _Axis:
    """Create a uniformly spaced CDAT-like longitude axis.

    This helper preserves the subset of the legacy
    ``cdms2.createUniformLongitudeAxis`` behavior used by ENSO_metrics while
    returning an ``_Axis`` object for the xarray compatibility layer.

    Parameters
    ----------
    start : float
        First longitude coordinate value.

    n : int
        Number of longitude points.

    delta : float
        Longitude spacing in degrees.

    Returns
    -------
    _Axis
        Longitude axis with CF-style metadata, including ``axis="X"``,
        ``standard_name="longitude"``, and ``units="degrees_east"``.
    """
    vals = np.array([start + i * delta for i in range(int(n))])
    return _Axis("lon", vals, units="degrees_east", axis_type="X")


def create_rect_grid(lat_axis: _Axis, lon_axis: _Axis, order: str = "yx", grid_type: str = "generic", mask=None) -> _Grid:
    """Create a CDAT-like rectilinear grid for the compatibility layer.

    This helper preserves the subset of the legacy ``cdms2.createRectGrid``
    behavior used by ENSO_metrics while returning an ``_Grid`` object backed by
    the xarray compatibility layer.

    Parameters
    ----------
    lat_axis : _Axis
        Latitude axis.

    lon_axis : _Axis
        Longitude axis.

    order : str, optional
        Axis order for the grid. The default ``"yx"`` means latitude followed
        by longitude.
        default value = ``"yx"``

    grid_type : str, optional
        Descriptive grid type, for example ``"generic"``, ``"gaussian"``,
        ``"uniform"``, or ``"equalarea"``.
        default value = ``"generic"``

    mask : array-like, optional
        Optional grid mask to attach to the returned grid.
        default value = None

    Returns
    -------
    _Grid
        CDAT-like rectilinear grid object used by ``CDATVariable``.
    """
    g = _Grid(lat_axis, lon_axis)
    g.type = grid_type
    return g


def create_variable(data, axes: Optional[list] = None, grid: Optional[_Grid] = None, mask=None, id: str = "", attributes: Optional[dict] = None) -> CDATVariable:
    """Create a CDAT-like variable for the xarray compatibility layer.

    This helper preserves the subset of the legacy ``cdms2.createVariable``
    behavior used by ENSO_metrics while returning a ``CDATVariable`` backed by
    the modern compatibility layer.

    Parameters
    ----------
    data : array-like
        Input data values.

    axes : list, optional
        List of CDAT-like axes associated with the dimensions of ``data``.
        default value = None

    grid : _Grid, optional
        CDAT-like rectilinear grid associated with the variable.
        default value = None

    mask : array-like, optional
        Optional mask to apply to ``data``.
        default value = None

    id : str, optional
        Variable identifier/name.
        default value = ``""``

    attributes : dict, optional
        Metadata attributes to attach to the variable.
        default value = None

    Returns
    -------
    CDATVariable
        CDAT-like variable object used by the refactored ENSO_metrics workflow.
    """
    return CDATVariable(data, axes=axes, grid=grid, mask=mask, id=id, attributes=attributes)


# ---------------------------------------------------------------------------
# Conversion helpers
# ---------------------------------------------------------------------------


def _extract_aux_lat_lon_axes(da: xr.DataArray):
    """
    Best-effort support for auxiliary 1D lat/lon coordinates.

    This helps with some post-processed products, but does not make native
    unstructured grids equivalent to rectilinear CDAT grids.

    **Raises** ``NotImplementedError`` if any auxiliary coordinate has >1
    dimension, because 2D lat/lon arrays indicate a curvilinear or unstructured
    grid that is not supported by this shim.  The caller should regrid the data
    to a regular lat-lon grid before converting via da_to_cdat.
    """
    lat_ax = None
    lon_ax = None
    for cname, coord in da.coords.items():
        if coord.ndim > 1:
            ax_type = _dim_to_axis_type(cname, coord)
            if ax_type in ("Y", "X"):
                raise NotImplementedError(
                    f"da_to_cdat: coordinate {cname!r} has {coord.ndim} dimensions "
                    f"(shape {coord.shape}), indicating a curvilinear or unstructured "
                    "grid that is not supported by this CDAT compatibility shim. "
                    "Regrid to a regular rectilinear lat-lon grid first "
                    "(e.g. using xESMF or Regrid())."
                )
            # Non-lat/lon 2-D coordinates (e.g. time_bnds, vertices) are fine;
            # skip them silently.
            continue
        ax_type = _dim_to_axis_type(cname, coord)
        attrs = dict(coord.attrs)
        if ax_type == "Y" and lat_ax is None:
            lat_ax = _Axis(cname, coord.values, units=attrs.get("units", ""), attributes=attrs, axis_type="Y")
        elif ax_type == "X" and lon_ax is None:
            lon_ax = _Axis(cname, coord.values, units=attrs.get("units", ""), attributes=attrs, axis_type="X")
    return lat_ax, lon_ax


def da_to_cdat(da: xr.DataArray, varname: Optional[str] = None) -> CDATVariable:
    """
    Convert an ``xarray.DataArray`` to a ``CDATVariable`` with strengthened
    axis metadata.

    This routine is the main boundary between xarray and legacy CDAT-style
    code. It makes the axis structure explicit and validates it before
    returning. Dask-backed arrays are materialized because the shim uses eager
    ``numpy.ma`` storage by design.
    """
    name = varname or da.name or ""

    # Reject curvilinear/unstructured grids eagerly: if any coordinate that
    # looks like latitude or longitude has more than one dimension, the data
    # are on a non-rectilinear grid that this shim cannot represent correctly.
    # Dim-coordinates (1D) are fine; only non-dim 2D+ coordinates trigger this.
    for cname, coord in da.coords.items():
        if coord.ndim > 1:
            ax_type = _dim_to_axis_type(cname, coord)
            if ax_type in ("Y", "X"):
                raise NotImplementedError(
                    f"da_to_cdat({name!r}): coordinate {cname!r} has "
                    f"{coord.ndim} dimensions (shape {coord.shape}), indicating "
                    "a curvilinear or unstructured grid that is not supported by "
                    "this CDAT compatibility shim. Regrid to a regular "
                    "rectilinear lat-lon grid first (e.g. using xESMF or Regrid())."
                )

    axes = []
    for dim in da.dims:
        coord = da.coords.get(dim)
        ax_type = _dim_to_axis_type(dim, coord)
        if coord is None:
            vals = np.arange(da.sizes[dim])
            attrs = {"axis": ax_type} if ax_type in {"T", "Y", "X", "Z"} else {}
            units = attrs.get("units", "")
        else:
            vals = coord.values
            attrs = _coord_attrs_with_axis(coord, ax_type)
            units = attrs.get("units", "")
        ax = _Axis(dim, vals, units=units, attributes=attrs, axis_type=ax_type)
        if coord is not None:
            cal = attrs.get("calendar") or getattr(coord, "encoding", {}).get("calendar", None)
            if cal is not None:
                ax.calendar = cal
        axes.append(ax)

    lat_ax = next((ax for ax in axes if ax is not None and ax.isLatitude()), None)
    lon_ax = next((ax for ax in axes if ax is not None and ax.isLongitude()), None)

    # Best-effort fallback for simple 1D auxiliary lat/lon coordinates.
    if lat_ax is None or lon_ax is None:
        aux_lat, aux_lon = _extract_aux_lat_lon_axes(da)
        lat_ax = lat_ax or aux_lat
        lon_ax = lon_ax or aux_lon

    grid = _Grid(lat_ax, lon_ax) if (lat_ax is not None and lon_ax is not None) else None

    attrs = dict(da.attrs)
    fill_value = attrs.get("_FillValue", attrs.get("missing_value", 1e20))
    try:
        fill_value = float(np.ravel(fill_value)[0])
    except Exception:
        fill_value = 1e20

    try:
        data = da.to_masked_array(copy=False)
    except Exception:
        raw = da.values
        if np.issubdtype(np.asarray(raw).dtype, np.number):
            data = ma.masked_invalid(raw)
        else:
            data = ma.array(raw)

    var = CDATVariable(data, axes=axes, grid=grid, id=name, attributes=attrs, fill_value=fill_value)
    if "units" in attrs:
        var.units = attrs.get("units", "")
    return validate_cdat_variable(var, context=f"da_to_cdat({name})")


def cdat_to_da(var: CDATVariable, name: Optional[str] = None) -> xr.DataArray:
    """
    Convert a ``CDATVariable`` to an ``xarray.DataArray`` while preserving CF
    axis metadata on coordinates.

    Masked floating data are represented with NaN. Masked non-floating numeric
    data are promoted to float only when a mask is present, avoiding invalid
    integer fill values while keeping unmasked integer arrays integer.
    """
    name = name or var.id or "var"
    validate_cdat_variable(var, context=f"cdat_to_da({name}) input")
    axes = _coerce_axes(var._axes, var._data.ndim)
    _validate_axes_shape(var._data, axes, context=f"cdat_to_da({name})")

    dims = [ax.id if ax is not None else f"dim_{i}" for i, ax in enumerate(axes)]
    coords = {}
    for dim, ax in zip(dims, axes):
        if ax is None:
            continue
        vals = ax._values
        attrs = dict(ax._attributes)
        if "units" not in attrs and ax.units:
            attrs["units"] = ax.units

        if ax.axis in ("T", "Y", "X", "Z"):
            attrs["axis"] = ax.axis
        else:
            detected = _detect_axis_type(ax.id)
            if detected in ("T", "Y", "X", "Z"):
                attrs["axis"] = detected

        if attrs.get("axis") == "Y":
            attrs.setdefault("standard_name", "latitude")
            attrs.setdefault("units", "degrees_north")
        elif attrs.get("axis") == "X":
            attrs.setdefault("standard_name", "longitude")
            attrs.setdefault("units", "degrees_east")
        elif attrs.get("axis") == "T":
            attrs.setdefault("standard_name", "time")
            if ax.calendar is not None:
                attrs.setdefault("calendar", ax.calendar)
            # Decode numeric time values to cftime objects so that xcdat's
            # add_missing_bounds (and other tools) see proper datetime-like
            # coordinates.  After toRelativeTime() the axis holds float64
            # "days since ..." values; we decode them here at the xarray
            # boundary so callers never receive a raw-numeric time coord.
            try:
                if np.issubdtype(np.asarray(vals).dtype, np.number):
                    units_t = attrs.get("units", ax.units or "")
                    cal_t = ax.calendar or attrs.get("calendar", "standard")
                    if units_t and "since" in units_t:
                        decoded = _decode_times_safe(
                            np.asarray(vals, dtype=float), units_t, cal_t
                        )
                        if decoded:
                            vals = np.array(decoded, dtype=object)
                            # Keep units/calendar in encoding, not attrs,
                            # so xarray does not try to re-encode them.
                            attrs.pop("units", None)
            except Exception:
                pass  # fallback: leave vals as-is (float), emit no warning

        clean = _clean_attrs(attrs)
        try:
            coords[dim] = xr.Variable(dim, vals, attrs=clean)
        except Exception:
            coords[dim] = xr.Variable(dim, np.arange(len(vals)), attrs=clean)

    data = var._data
    if ma.isMaskedArray(data) and np.any(ma.getmaskarray(data)):
        if np.issubdtype(data.dtype, np.floating):
            raw = data.filled(np.nan)
        elif np.issubdtype(data.dtype, np.number):
            raw = data.astype(float).filled(np.nan)
        else:
            raw = data.filled(None)
    else:
        raw = np.asarray(data)

    return xr.DataArray(raw, dims=dims, coords=coords, name=name, attrs=_clean_attrs(var._attributes))
