# -*- coding:UTF-8 -*-
"""
XarrayCompat.py
===============
Backward-compatible drop-in replacements for the retired CDAT/UV-CDAT
``cdms2.TransientVariable`` and associated axis/grid objects.

All computation is done with ``numpy.ma``; coordinate metadata is stored
alongside the data so callers that use CDAT-style introspection
(``.getAxisList()``, ``.getGrid()``, ``.getTime().asComponentTime()``, …)
continue to work without modification.

Exported public API
-------------------
CDATVariable   – replaces cdms2.TransientVariable
_Axis          – replaces cdms2.Axis
_TimeAxis      – replaces cdms2.Axis (time-specific)
_Grid          – replaces cdms2.RectGrid

Factory helpers (replacements for cdms2 creation functions)
-----------------------------------------------------------
create_axis(id, values, units='', attributes=None)
create_uniform_lat_axis(start, n, delta)
create_uniform_lon_axis(start, n, delta)
create_rect_grid(lat_axis, lon_axis, order='yx', grid_type='generic')
create_variable(data, axes=None, grid=None, mask=None, id='',
                attributes=None)

Conversion helpers
------------------
da_to_cdat(da, varname=None)   – xr.DataArray  → CDATVariable
cdat_to_da(var)                – CDATVariable  → xr.DataArray
"""

from __future__ import annotations

import copy
import re
from typing import List, Optional, Sequence, Union

import cftime
import numpy as np
import numpy.ma as ma
import xarray as xr

__all__ = [
    "CDATVariable",
    "_Axis",
    "_TimeAxis",
    "_Grid",
    "create_axis",
    "create_uniform_lat_axis",
    "create_uniform_lon_axis",
    "create_rect_grid",
    "create_variable",
    "da_to_cdat",
    "cdat_to_da",
]

# ---------------------------------------------------------------------------
# Axis
# ---------------------------------------------------------------------------

_LAT_IDS  = {"lat", "latitude", "j", "y", "Y", "yt_ocean", "yu_ocean"}
_LON_IDS  = {"lon", "longitude", "i", "x", "X", "xt_ocean", "xu_ocean"}
_TIME_IDS = {"time", "t", "T"}
_LEV_IDS  = {"lev", "level", "depth", "plev", "z", "Z", "st_ocean",
             "sw_ocean"}


def _detect_axis_type(ax_id: str) -> str:
    if ax_id in _TIME_IDS or "time" in ax_id.lower():
        return "T"
    if ax_id in _LAT_IDS or "lat" in ax_id.lower():
        return "Y"
    if ax_id in _LON_IDS or "lon" in ax_id.lower():
        return "X"
    if ax_id in _LEV_IDS or "lev" in ax_id.lower() or "depth" in ax_id.lower():
        return "Z"
    return "-"


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
        # Store raw values; may be cftime objects for time axes
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
        self.axis = axis_type or _detect_axis_type(id)
        self.long_name = self._attributes.get("long_name", id)
        # Optional CDAT-style extra attributes (set by callers)
        self.regions: Optional[str] = None
        self.reference: Optional[str] = None
        self.calendar: Optional[str] = self._attributes.get("calendar", None)

    # ------------------------------------------------------------------
    # Type predicates
    # ------------------------------------------------------------------
    def isTime(self) -> bool:      return self.axis == "T"
    def isLatitude(self) -> bool:  return self.axis == "Y"
    def isLongitude(self) -> bool: return self.axis == "X"
    def isLevel(self) -> bool:     return self.axis == "Z"

    # ------------------------------------------------------------------
    # Array-like protocol
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # cdms2-compatible methods
    # ------------------------------------------------------------------
    def asComponentTime(self) -> list:
        """Return list of cftime datetime objects (mirrors cdtime behaviour)."""
        vals = self._values
        if len(vals) == 0:
            return []
        if isinstance(vals[0], (cftime.datetime,)):
            # Already decoded — clamp any second=60 that somehow slipped through
            result = []
            for dt in vals:
                if hasattr(dt, 'second') and dt.second == 60:
                    dt = dt.replace(second=59)
                result.append(dt)
            return result
        # Decode numeric values element-by-element so that a single leap-second
        # (second=60) does not abort the entire array.  cftime raises ValueError
        # when constructing e.g. DatetimeNoLeap(... second=60); we catch it
        # per-element and subtract 1 second from the raw numeric value so the
        # re-decode produces second=59.  This handles every calendar type
        # (noleap, gregorian, proleptic_gregorian, …).
        if self.units:
            try:
                cal = self.calendar or "standard"
                arr = np.asarray(vals, dtype=float)
                result = []
                for v in arr:
                    try:
                        dt = cftime.num2date(v, self.units, calendar=cal)
                        if getattr(dt, 'second', 0) == 60:
                            dt = cftime.num2date(
                                v - 1.0 / 86400, self.units, calendar=cal)
                    except ValueError:
                        # second=60: subtract 1 s and re-decode
                        dt = cftime.num2date(
                            v - 1.0 / 86400, self.units, calendar=cal)
                    result.append(dt)
                return result
            except Exception:
                pass
        # Fallback: wrap floats as years
        return [cftime.datetime(int(v), 1, 1) for v in vals]

    def toRelativeTime(self, units: str):
        """Convert cftime values to numeric relative time in-place."""
        try:
            cal = self.calendar or "standard"
            comp = self.asComponentTime()
            self._values = cftime.date2num(comp, units, calendar=cal)
            self.units = units
        except Exception:
            pass

    def copy(self) -> "_Axis":
        return _Axis(
            self.id,
            self._values.copy(),
            units=self.units,
            attributes=dict(self._attributes),
            axis_type=self.axis,
        )


class _TimeAxis(_Axis):
    """Time-specific axis – identical to _Axis but always typed 'T'."""

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
    Drop-in replacement for ``cdms2.TransientVariable``.

    Stores data as a ``numpy.ma.MaskedArray`` and carries coordinate
    metadata (list of ``_Axis`` objects + optional ``_Grid``) so that
    all CDAT-style introspection methods continue to work.
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
        # ---- data --------------------------------------------------------
        if isinstance(data, CDATVariable):
            raw = data._data.copy()
        elif isinstance(data, ma.MaskedArray):
            raw = data.copy()
        else:
            raw = ma.array(np.asarray(data, dtype=float), fill_value=fill_value)
        if mask is not None:
            raw = ma.array(raw.data, mask=mask, fill_value=fill_value)
        self._data: ma.MaskedArray = raw

        # ---- metadata ----------------------------------------------------
        self._axes: List[_Axis] = _coerce_axes(axes, raw.ndim)
        self._grid: Optional[_Grid] = grid
        self.id: str = id
        self.name: str = id
        self._attributes: dict = dict(attributes or {})

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # Numpy / masked-array protocol
    # ------------------------------------------------------------------
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
            yield self._wrap(self._data[i])

    def __repr__(self):
        return (
            f"CDATVariable(id={self.id!r}, shape={self.shape}, "
            f"dtype={self.dtype})"
        )

    def __str__(self):
        return str(self._data)

    # ------------------------------------------------------------------
    # Indexing / slicing
    # ------------------------------------------------------------------
    def __getitem__(self, key):
        result = self._data[key]
        if not isinstance(result, np.ndarray) and not isinstance(result, ma.MaskedArray):
            return result  # scalar
        new_axes = self._sliced_axes(key, result.shape)
        return CDATVariable(
            result,
            axes=new_axes,
            grid=self._grid,
            id=self.id,
            attributes=dict(self._attributes),
        )

    def __setitem__(self, key, value):
        if isinstance(value, CDATVariable):
            self._data[key] = value._data
        else:
            self._data[key] = value

    def _sliced_axes(self, key, new_shape):
        """Best-effort: return axes for a sliced result."""
        if not self._axes:
            return []
        if isinstance(key, tuple):
            new_axes = []
            ax_idx = 0
            for k in key:
                if ax_idx >= len(self._axes):
                    break
                ax = self._axes[ax_idx]
                if isinstance(k, int):
                    ax_idx += 1
                    continue  # dimension collapsed
                elif isinstance(k, slice):
                    new_vals = ax._values[k]
                    new_axes.append(
                        _Axis(ax.id, new_vals, units=ax.units,
                              attributes=ax._attributes, axis_type=ax.axis)
                    )
                    ax_idx += 1
                elif isinstance(k, (list, np.ndarray)):
                    new_vals = ax._values[k]
                    new_axes.append(
                        _Axis(ax.id, new_vals, units=ax.units,
                              attributes=ax._attributes, axis_type=ax.axis)
                    )
                    ax_idx += 1
            return new_axes
        # Single index/slice
        if isinstance(key, int):
            return self._axes[1:] if len(self._axes) > 1 else []
        if isinstance(key, (slice, list, np.ndarray)):
            ax = self._axes[0]
            new_vals = ax._values[key]
            new_ax = _Axis(ax.id, new_vals, units=ax.units,
                           attributes=ax._attributes, axis_type=ax.axis)
            return [new_ax] + self._axes[1:]
        return self._axes

    # ------------------------------------------------------------------
    # Arithmetic operators (all return CDATVariable)
    # ------------------------------------------------------------------
    def _unpack(self, other):
        if isinstance(other, CDATVariable):
            return other._data
        return other

    def _wrap(self, data) -> "CDATVariable":
        if not isinstance(data, (np.ndarray, ma.MaskedArray)):
            return data
        return CDATVariable(
            data,
            axes=list(self._axes),
            grid=self._grid,
            id=self.id,
            attributes=dict(self._attributes),
        )

    def __add__(self, other):       return self._wrap(self._data + self._unpack(other))
    def __radd__(self, other):      return self._wrap(self._unpack(other) + self._data)
    def __sub__(self, other):       return self._wrap(self._data - self._unpack(other))
    def __rsub__(self, other):      return self._wrap(self._unpack(other) - self._data)
    def __mul__(self, other):       return self._wrap(self._data * self._unpack(other))
    def __rmul__(self, other):      return self._wrap(self._unpack(other) * self._data)
    def __truediv__(self, other):   return self._wrap(self._data / self._unpack(other))
    def __rtruediv__(self, other):  return self._wrap(self._unpack(other) / self._data)
    def __neg__(self):              return self._wrap(-self._data)
    def __abs__(self):              return self._wrap(abs(self._data))
    def __pow__(self, exp):         return self._wrap(self._data ** exp)

    def __gt__(self, other):  return self._data > self._unpack(other)
    def __lt__(self, other):  return self._data < self._unpack(other)
    def __ge__(self, other):  return self._data >= self._unpack(other)
    def __le__(self, other):  return self._data <= self._unpack(other)
    def __eq__(self, other):  return self._data == self._unpack(other)
    def __ne__(self, other):  return self._data != self._unpack(other)

    # ------------------------------------------------------------------
    # numpy.ma delegation
    # ------------------------------------------------------------------
    def filled(self, fill_value=1e20):
        return self._data.filled(fill_value)

    def fill(self, value):
        self._data.fill(value)

    def astype(self, dtype) -> "CDATVariable":
        return self._wrap(self._data.astype(dtype))

    def squeeze(self, axis=None) -> "CDATVariable":
        result = self._data.squeeze(axis=axis)
        new_axes = [ax for ax in self._axes if len(ax) > 1]
        return CDATVariable(result, axes=new_axes, grid=self._grid,
                            id=self.id, attributes=dict(self._attributes))

    def compress(self, condition, axis: int = 0) -> "CDATVariable":
        result = self._data.compress(condition, axis=axis)
        new_axes = list(self._axes)
        if self._axes and axis < len(self._axes):
            old_ax = self._axes[axis]
            new_vals = old_ax._values[np.asarray(condition, dtype=bool)]
            new_axes[axis] = _Axis(old_ax.id, new_vals, units=old_ax.units,
                                   attributes=old_ax._attributes,
                                   axis_type=old_ax.axis)
        return CDATVariable(result, axes=new_axes, grid=self._grid,
                            id=self.id, attributes=dict(self._attributes))

    def copy(self) -> "CDATVariable":
        return CDATVariable(
            self._data.copy(),
            axes=[ax.copy() for ax in self._axes],
            grid=self._grid,
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

    def setAxisList(self, axes: list):
        self._axes = list(axes)

    def getGrid(self) -> Optional[_Grid]:
        return self._grid

    def setGrid(self, grid: Optional[_Grid]):
        self._grid = grid

    def getTime(self) -> Optional[_Axis]:
        for ax in self._axes:
            if ax is not None and ax.isTime():
                return ax
        return None

    def getLatitude(self) -> Optional[_Axis]:
        for ax in self._axes:
            if ax is not None and ax.isLatitude():
                return ax
        return None

    def getLongitude(self) -> Optional[_Axis]:
        for ax in self._axes:
            if ax is not None and ax.isLongitude():
                return ax
        return None

    def getLevel(self) -> Optional[_Axis]:
        for ax in self._axes:
            if ax is not None and ax.isLevel():
                return ax
        return None

    def getOrder(self) -> str:
        order = ""
        for ax in self._axes:
            if ax is None:
                order += "-"
            elif ax.isTime():
                order += "t"
            elif ax.isLatitude():
                order += "y"
            elif ax.isLongitude():
                order += "x"
            elif ax.isLevel():
                order += "z"
            else:
                order += "-"
        return order

    def reorder(self, order: str) -> "CDATVariable":
        """
        Reorder axes.  Accepts:
        * full permutation string: 'tyx', 'txy', '10', '210', …
        * 't...'  – move time to first position
        * '...t'  – move time to last position
        * '10'    – numeric index permutation
        """
        ndim = self._data.ndim
        if ndim == 0:
            return self.copy()

        # ---- resolve permutation ----------------------------------------
        if order in ("t...", "T..."):
            t_n = next((i for i, ax in enumerate(self._axes) if ax is not None and ax.isTime()), 0)
            perm = [t_n] + [i for i in range(ndim) if i != t_n]
        elif order in ("...t", "...T"):
            t_n = next((i for i, ax in enumerate(self._axes) if ax is not None and ax.isTime()), ndim - 1)
            perm = [i for i in range(ndim) if i != t_n] + [t_n]
        elif all(c.isdigit() for c in order):
            perm = [int(c) for c in order]
        else:
            char_map: dict = {}
            for i, ax in enumerate(self._axes):
                if ax is None:
                    continue
                if ax.isTime():
                    char_map["t"] = i
                elif ax.isLatitude():
                    char_map["y"] = i
                elif ax.isLongitude():
                    char_map["x"] = i
                elif ax.isLevel():
                    char_map["z"] = i
            perm = []
            used = set()
            remaining = list(range(ndim))
            for c in order.lower():
                if c == ".":
                    continue
                if c in char_map and char_map[c] not in used:
                    perm.append(char_map[c])
                    used.add(char_map[c])
            # Fill in remaining dims for '...' parts
            for i in remaining:
                if i not in used:
                    perm.append(i)
                    used.add(i)
            if len(perm) != ndim:
                return self.copy()

        new_data = np.transpose(self._data, perm)
        new_axes = [self._axes[i] if i < len(self._axes) else None for i in perm]

        # Rebuild grid if lat/lon moved
        lat_ax = next((ax for ax in new_axes if ax is not None and ax.isLatitude()), None)
        lon_ax = next((ax for ax in new_axes if ax is not None and ax.isLongitude()), None)
        new_grid = _Grid(lat_ax, lon_ax) if (lat_ax and lon_ax) else self._grid

        return CDATVariable(
            new_data,
            axes=new_axes,
            grid=new_grid,
            id=self.id,
            attributes=dict(self._attributes),
        )

    # ------------------------------------------------------------------
    # CDAT-style callable selection  var(time=..., latitude=..., longitude=...)
    # ------------------------------------------------------------------
    def __call__(self, *args, **kwargs) -> "CDATVariable":
        if not kwargs:
            return self.copy()

        result_data = self._data.copy()
        result_axes = [ax.copy() if ax is not None else None for ax in self._axes]

        if "squeeze" in kwargs and kwargs["squeeze"]:
            result_data = result_data.squeeze()
            result_axes = [ax for ax in result_axes if ax is not None and len(ax) > 1]
            return CDATVariable(result_data, axes=result_axes, grid=self._grid,
                                id=self.id, attributes=dict(self._attributes))

        def _sel_axis(ax_idx, bounds):
            nonlocal result_data, result_axes
            ax = result_axes[ax_idx]
            if ax is None:
                return
            if ax.isTime():
                t_vals = ax.asComponentTime()
                t0_str = str(bounds[0]).split(".")[0]
                t1_str = str(bounds[1]).split(".")[0]
                indices = [i for i, t in enumerate(t_vals)
                           if t0_str <= str(t).split(".")[0] <= t1_str]
            else:
                raw = ax._values.astype(float)
                lo, hi = min(float(bounds[0]), float(bounds[1])), max(float(bounds[0]), float(bounds[1]))
                indices = list(np.where((raw >= lo) & (raw <= hi))[0])
            if not indices:
                return
            slices = [slice(None)] * result_data.ndim
            slices[ax_idx] = indices
            result_data = result_data[tuple(slices)]
            old_ax = ax
            new_vals = old_ax._values[indices]
            result_axes[ax_idx] = _Axis(old_ax.id, new_vals, units=old_ax.units,
                                        attributes=old_ax._attributes,
                                        axis_type=old_ax.axis)

        if "time" in kwargs:
            t_idx = next((i for i, ax in enumerate(result_axes)
                          if ax is not None and ax.isTime()), None)
            if t_idx is not None:
                _sel_axis(t_idx, kwargs["time"])

        if "latitude" in kwargs:
            lat_idx = next((i for i, ax in enumerate(result_axes)
                            if ax is not None and ax.isLatitude()), None)
            if lat_idx is not None:
                _sel_axis(lat_idx, kwargs["latitude"])

        if "longitude" in kwargs:
            lon_idx = next((i for i, ax in enumerate(result_axes)
                            if ax is not None and ax.isLongitude()), None)
            if lon_idx is not None:
                _sel_axis(lon_idx, kwargs["longitude"])

        lat_ax = next((ax for ax in result_axes if ax is not None and ax.isLatitude()), None)
        lon_ax = next((ax for ax in result_axes if ax is not None and ax.isLongitude()), None)
        new_grid = _Grid(lat_ax, lon_ax) if (lat_ax and lon_ax) else self._grid

        return CDATVariable(result_data, axes=result_axes, grid=new_grid,
                            id=self.id, attributes=dict(self._attributes))


# ---------------------------------------------------------------------------
# Factory helpers  (replacements for cdms2 factory functions)
# ---------------------------------------------------------------------------

def create_axis(
    id: str,
    values,
    units: str = "",
    attributes: Optional[dict] = None,
) -> _Axis:
    """Replacement for ``cdms2.createAxis``."""
    return _Axis(id, values, units=units, attributes=attributes)


def create_uniform_lat_axis(start: float, n: int, delta: float) -> _Axis:
    """Replacement for ``cdms2.createUniformLatitudeAxis``."""
    vals = np.array([start + i * delta for i in range(n)])
    return _Axis("lat", vals, units="degrees_north", axis_type="Y")


def create_uniform_lon_axis(start: float, n: int, delta: float) -> _Axis:
    """Replacement for ``cdms2.createUniformLongitudeAxis``."""
    vals = np.array([start + i * delta for i in range(n)])
    return _Axis("lon", vals, units="degrees_east", axis_type="X")


def create_rect_grid(
    lat_axis: _Axis,
    lon_axis: _Axis,
    order: str = "yx",
    grid_type: str = "generic",
    mask=None,
) -> _Grid:
    """Replacement for ``cdms2.createRectGrid``."""
    g = _Grid(lat_axis, lon_axis)
    g.type = grid_type
    return g


def create_variable(
    data,
    axes: Optional[list] = None,
    grid: Optional[_Grid] = None,
    mask=None,
    id: str = "",
    attributes: Optional[dict] = None,
) -> CDATVariable:
    """Replacement for ``cdms2.createVariable``."""
    return CDATVariable(data, axes=axes, grid=grid, mask=mask,
                        id=id, attributes=attributes)


# ---------------------------------------------------------------------------
# Conversion helpers
# ---------------------------------------------------------------------------

def _dim_to_axis_type(dim_name: str, coord) -> str:
    typ = _detect_axis_type(dim_name)
    if typ != "-":
        return typ
    if coord is not None:
        cf_axis = coord.attrs.get("axis", "")
        if cf_axis.upper() == "T":
            return "T"
        if cf_axis.upper() == "Y":
            return "Y"
        if cf_axis.upper() == "X":
            return "X"
        if cf_axis.upper() == "Z":
            return "Z"
    return "-"


def da_to_cdat(da: xr.DataArray, varname: Optional[str] = None) -> CDATVariable:
    """
    Convert an ``xarray.DataArray`` to a ``CDATVariable``.

    Time coordinates containing ``cftime`` objects are preserved as-is.
    """
    name = varname or da.name or ""

    # ---- build axes ---------------------------------------------------
    axes = []
    for dim in da.dims:
        coord = da.coords.get(dim)
        ax_type = _dim_to_axis_type(dim, coord)
        if coord is None:
            ax = _Axis(dim, np.arange(da.sizes[dim]), axis_type=ax_type)
        else:
            vals = coord.values
            units = coord.attrs.get("units", "")
            cal = coord.attrs.get("calendar", None)
            attrs = dict(coord.attrs)
            ax = _Axis(dim, vals, units=units,
                       attributes=attrs, axis_type=ax_type)
            if cal:
                ax.calendar = cal
        axes.append(ax)

    # ---- build grid ---------------------------------------------------
    lat_ax = next((ax for ax in axes if ax.isLatitude()), None)
    lon_ax = next((ax for ax in axes if ax.isLongitude()), None)
    grid = _Grid(lat_ax, lon_ax) if (lat_ax and lon_ax) else None

    # ---- build data ---------------------------------------------------
    raw = da.values
    if hasattr(raw, "mask"):
        data = ma.array(raw.data, mask=raw.mask, fill_value=1e20)
    else:
        data = ma.array(np.asarray(raw, dtype=float),
                        mask=np.isnan(raw.astype(float)),
                        fill_value=1e20)

    var = CDATVariable(data, axes=axes, grid=grid,
                       id=name, attributes=dict(da.attrs))
    var.units = da.attrs.get("units", "")
    return var


def _coerce_axes(axes, ndim: int) -> list:
    """
    Ensure every element in an axes list is an _Axis (or None).
    Raw numpy arrays (e.g. from ``axis[:]``) are wrapped in a generic _Axis.
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


def cdat_to_da(var: CDATVariable, name: Optional[str] = None) -> xr.DataArray:
    """
    Convert a ``CDATVariable`` to an ``xarray.DataArray``.

    Numeric time axes are left as-is; cftime-valued time axes are placed
    into the DataArray coords directly so xarray can handle them.
    """
    name = name or var.id or "var"
    # Ensure axes are proper _Axis objects (guard against any raw arrays that
    # slipped through before _coerce_axes was introduced)
    axes = _coerce_axes(var._axes, var._data.ndim)
    dims = [ax.id if ax is not None else f"dim_{i}"
            for i, ax in enumerate(axes)]
    coords = {}
    for i, (dim, ax) in enumerate(zip(dims, axes)):
        if ax is None:
            continue
        vals = ax._values
        attrs = dict(ax._attributes)
        if "units" not in attrs and ax.units:
            attrs["units"] = ax.units
        # Fill in canonical CF units when still missing, verified by value range
        if "units" not in attrs and ax.axis in ("Y", "X"):
            try:
                _v = np.asarray(vals, dtype=float)
                _lo, _hi = float(_v.min()), float(_v.max())
                if ax.axis == "Y" and -90.0 <= _lo and _hi <= 90.0:
                    attrs["units"] = "degrees_north"
                elif ax.axis == "X" and -180.0 <= _lo and _hi <= 360.0:
                    attrs["units"] = "degrees_east"
            except Exception:
                pass
        if ax.axis in ("T", "Y", "X", "Z"):
            attrs["axis"] = ax.axis
        try:
            coords[dim] = xr.Variable(dim, vals, attrs=attrs)
        except Exception:
            coords[dim] = xr.Variable(dim, np.arange(len(vals)), attrs=attrs)

    raw = np.ma.filled(var._data, fill_value=np.nan)
    da = xr.DataArray(raw, dims=dims, coords=coords,
                      name=name, attrs=dict(var._attributes))
    return da
