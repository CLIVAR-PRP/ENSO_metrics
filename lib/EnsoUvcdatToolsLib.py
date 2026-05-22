# -*- coding:UTF-8 -*-
from calendar import monthrange
import copy
import warnings
from datetime import date
from inspect import stack as INSPECTstack
from packaging.version import Version
import ntpath
import os as _os_sn
import tempfile as _tmp_sn
import shutil as _shutil_sn
import numpy
from numpy import array as NParray
from numpy import exp as NPexp
from numpy import histogram as NPhistogram
from numpy import isnan as NPisnan
from numpy import nan as NPnan
from numpy import nonzero as NPnonzero
from numpy import ones as NPones

if Version(numpy.__version__) < Version('1.25.0'):
    from numpy import product as NPproduct
else:
    from numpy import prod as NPproduct

from numpy import where as NPwhere
from numpy.ma.core import MaskedArray as NPma__core__MaskedArray
from os.path import isdir as OSpath_isdir
from os.path import isfile as OSpath__isfile
from os.path import join as OSpath__join
from os.path import split as OSpath__split
from scipy.signal import detrend as SCIPYsignal_detrend  # kept for inline use
from scipy.stats import skew as SCIPYstats__skew
from sys import prefix as SYS_prefix

# ENSO_metrics package functions:
from .EnsoCollectionsLib import CmipVariables
from .EnsoCollectionsLib import ReferenceObservations
from .EnsoCollectionsLib import ReferenceRegions
from . import EnsoErrorsWarnings
from .EnsoToolsLib import add_up_errors, find_xy_min_max, string_in_dict

# ---------------------------------------------------------------------------
# New-stack imports  (replaces retired CDAT/UV-CDAT packages)
# ---------------------------------------------------------------------------
import numpy as np                                 # replaces MV2 numeric ops
import numpy.ma as ma                              # replaces MV2 masked ops
import xarray as xr                                # replaces cdms2 variable/axis
try:
    import xcdat as _xcdat                         # registers xr.Dataset.spatial / .temporal
except ImportError:
    _xcdat = None
from scipy.stats import linregress as _linregress  # replaces genutil.linearregression
# scipy.signal.detrend / scipy.stats.skew imported below as SCIPYsignal_detrend / SCIPYstats__skew
try:
    import regionmask as _regionmask               # replaces cdutil.generateLandSeaMask
    _HAS_REGIONMASK = True
except ImportError:
    _HAS_REGIONMASK = False
try:
    import xesmf as _xesmf                         # replaces regrid2 / cdms2.regrid
    _HAS_XESMF = True
except (ImportError, OSError):
    # OSError can occur when esmpy's libesmf_fullylinked.so is missing/mislinked
    _HAS_XESMF = False

import cftime as _cft  # avoids module-level dependency

# When True (default), _guess_dim raises ValueError for dimensions whose axis
# type cannot be determined with any confidence (score <= 0). Set to False to
# demote the error to a warning and return the best-guess dim — useful for
# legacy observational datasets that lack CF axis/standard_name/units metadata.
#
#   import lib.EnsoUvcdatToolsLib as E; E.STRICT_DIM_GUESS = False
#
# **User Note:**
# If you work with legacy or poorly-annotated datasets, set STRICT_DIM_GUESS = False
# to avoid errors when axis metadata is missing. This will emit warnings instead.
# For best results, add CF-compliant axis/standard_name/units metadata to your data.
STRICT_DIM_GUESS: bool = True

# Compatibility shim: provides CDATVariable, a CDAT-like variable object
# used in place of the legacy cdms2.TransientVariable interface by ENSO_metrics.
from .XarrayCompat import (
    CDATVariable,
    _Axis,
    _build_grid_from_axes,
    _clean_attrs,
    _get_time_coder,
    _dim_to_axis_type,
    create_axis,
    create_uniform_lat_axis,
    create_uniform_lon_axis,
    create_rect_grid,
    create_variable,
    da_to_cdat,
    cdat_to_da,
    validate_cdat_variable,
)

def open_file(path, mode="r"):
    """Open a NetCDF file; returns an _XcDatasetHandle wrapping xcdat."""
    return _XcDatasetHandle(path, mode)

def CDTIMEcomptime(
            year, month=1, day=1, hour=0, minute=0, second=0.0,
            calendar="standard"
    ):
    """Replacement for cdtime.comptime()."""
    # Clamp leap-second (second=60) to 59 — cftime rejects second=60
    return _cft.datetime(
        year, month, day, hour, minute,
        min(int(second), 59), calendar=calendar
    )

def _coord_attrs_lower(coord):
    return {
        k.lower(): str(v).lower()
        for k, v in coord.attrs.items()
    }

def _detect_axis(da: xr.DataArray, axis_type: str) -> str:
    """
    Best-effort axis detection.

    Returns dim name or "".
    Never raises or warns.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            return _guess_dim(da, axis_type, strict=False)
        except Exception:
            return ""

def _validate_axis(da: xr.DataArray, axis_type: str, *, context: str = "") -> str:
    """
    Scientific validation layer.

    Requires at least heuristic confidence, i.e. score >= 1.
    Rejects value-range-only detection.
    """
    dim, score = _guess_dim(
        da,
        axis_type,
        strict=False,
        return_score=True,
    )

    if not dim or score <= 0:
        raise ValueError(
            f"{context}: Could not safely validate axis {axis_type!r}. "
            f"Detected dim={dim!r}, score={score}, "
            f"dims={list(da.dims)}, coords={list(da.coords)}. "
            "Add CF axis/standard_name/units metadata."
        )

    if axis_type.upper() == "T":
        coord = da.coords.get(dim)

        if coord is None:
            raise ValueError(
                f"{context}: Unsafe time-axis detection for dim {dim!r}. "
                "Time axis must have datetime-like values or CF time metadata."
            )

        attrs = _coord_attrs_lower(coord)

        if not (
            _is_datetime_like_time(coord)
            or "since" in attrs.get("units", "")
            or attrs.get("axis", "").upper() == "T"
            or attrs.get("standard_name", "") == "time"
        ):
            raise ValueError(
                f"{context}: Unsafe time-axis detection for dim {dim!r}. "
                "Time axis must have datetime-like values or CF time metadata."
            )

    return dim

def _require_axis(da: xr.DataArray, axis_type: str, *, context: str = "") -> str:
    """
    Hard enforcement layer.

    Requires metadata-supported detection, i.e. score >= 2.
    Use for time selection, regridding, EOFs, and strict reductions.
    """
    dim, score = _guess_dim(
        da,
        axis_type,
        strict=True,
        return_score=True,
    )

    if not dim or score < 2:
        raise ValueError(
            f"{context}: Required metadata-supported axis {axis_type!r} "
            f"not found. Detected dim={dim!r}, score={score}, "
            f"dims={list(da.dims)}, coords={list(da.coords)}. "
            "Add CF axis/standard_name/units metadata."
        )

    return dim

# ---------------------------------------------------------------------------
# MV2 aliases  → numpy.ma equivalents
#
# These functions retain the original CDAT MV2 names so that the 27 000-line
# EnsoMetricsLib.py function bodies work without modification.  They are
# purely internal; no user-facing API exposes these names.
# ---------------------------------------------------------------------------
def MV2add(a, b):            return _mv_wrap(ma.add(_mv(a), _mv(b)), a)
def MV2arange(*args):        return ma.array(np.arange(*args))
def MV2array(data, **kw):    return CDATVariable(ma.array(data, **kw), id="")
def MV2average(a, axis=None, weights=None):
    result = ma.average(_mv(a), axis=axis,
                        weights=_mv(weights) if weights is not None else None)
    if isinstance(a, CDATVariable) and isinstance(result, (np.ndarray, ma.MaskedArray)):
        if axis is None:
            new_axes = []
        else:
            ax = _axis_to_int(a, axis)
            drop = {int(i) for i in ax} if isinstance(ax, tuple) else {int(ax)}
            drop = {i if i >= 0 else a.ndim + i for i in drop}
            new_axes = [
                old_ax.copy() if old_ax is not None else None
                for i, old_ax in enumerate(a._axes)
                if i not in drop
            ]
        # Post-check: catch accidental time-axis loss only when time was
        # not in the set of axes that were explicitly reduced.
        _t_ax = a.getTime()
        _t_idx = next(
            (i for i, ax in enumerate(a._axes) if ax is not None and ax is _t_ax), None
        ) if _t_ax is not None else None
        _time_explicitly_dropped = (_t_idx is not None and _t_idx in drop)
        if (_t_ax is not None and result.ndim >= 1
                and not any(ax is not None and ax.axis == "T" for ax in new_axes)
                and not _time_explicitly_dropped):
            raise RuntimeError(
                f"Time axis unexpectedly dropped in MV2average: "
                f"id={a.id!r}, original_axes={[ax.id if ax else None for ax in a._axes]}, "
                f"axis={axis!r}, drop={drop}, result_shape={result.shape}"
            )
        return CDATVariable(
            result,
            axes=new_axes,
            grid=None if len(new_axes) < 2 else a._grid,
            id=a.id,
            attributes=dict(a._attributes),
        )
    return result
def MV2compress(condition, a, axis=0):
    if isinstance(a, CDATVariable):
        return a.compress(condition, axis=axis)
    cond = np.asarray(condition, dtype=bool)
    raw = _mv(a)
    return ma.array(np.compress(cond, raw, axis=axis),
                    mask=np.compress(cond, ma.getmaskarray(raw), axis=axis))
def MV2concatenate(seq, axis=0):
    seq = list(seq)
    arrs  = [_mv(x) for x in seq]
    masks = [ma.getmaskarray(x) for x in arrs]
    result = ma.array(
        np.concatenate(arrs, axis=axis),
        mask=np.concatenate(masks, axis=axis)
        )
    if seq and all(isinstance(x, CDATVariable) for x in seq):
        tmpl = seq[0]
        new_axes = list(tmpl._axes)
        ax_int = axis if isinstance(axis, (int, np.integer)) else 0
        if (0 <= ax_int < len(new_axes) and
                all(ax_int < len(x._axes) and x._axes[ax_int] is not None
                    for x in seq)):
            old_ax = tmpl._axes[ax_int]
            cat_vals = np.concatenate([x._axes[ax_int]._values for x in seq])
            new_ax = _Axis(
                old_ax.id, cat_vals,
                units=old_ax.units,
                attributes=dict(old_ax._attributes),
                axis_type=old_ax.axis
            )
            new_ax.calendar = old_ax.calendar
            new_axes[ax_int] = new_ax
        return CDATVariable(result, axes=new_axes, grid=tmpl._grid,
                            id=tmpl.id, attributes=dict(tmpl._attributes))
    return result
def MV2divide(a, b):         return _mv_wrap(ma.divide(_mv(a), _mv(b)), a)
def MV2masked_where(cond, a):return _mv_wrap(ma.masked_where(np.asarray(cond, dtype=bool), _mv(a)), a)
def MV2maximum(a):           return float(ma.max(_mv(a)))
def MV2minimum(a):           return float(ma.min(_mv(a)))
def MV2multiply(a, b):       return _mv_wrap(ma.multiply(_mv(a), _mv(b)), a)
def MV2ones(shape):          return CDATVariable(ma.ones(shape), id="")
def MV2subtract(a, b):       return _mv_wrap(ma.subtract(_mv(a), _mv(b)), a)
def MV2sum(a, axis=None, fill_value=0, dtype=None):
    result = ma.sum(_mv(a), axis=axis, dtype=dtype)
    if isinstance(a, CDATVariable) and isinstance(result, (np.ndarray, ma.MaskedArray)):
        if axis is None:
            new_axes = []
        else:
            ax = _axis_to_int(a, axis)
            drop = {int(i) for i in ax} if isinstance(ax, tuple) else {int(ax)}
            drop = {i if i >= 0 else a.ndim + i for i in drop}
            new_axes = [
                old_ax.copy() if old_ax is not None else None
                for i, old_ax in enumerate(a._axes)
                if i not in drop
            ]
        # Post-check: catch accidental time-axis loss only when time was
        # not in the set of axes that were explicitly reduced.
        _t_ax = a.getTime()
        _t_idx = next(
            (i for i, ax in enumerate(a._axes) if ax is not None and ax is _t_ax), None
        ) if _t_ax is not None else None
        _time_explicitly_dropped = (_t_idx is not None and _t_idx in drop)
        if (_t_ax is not None and result.ndim >= 1
                and not any(ax is not None and ax.axis == "T" for ax in new_axes)
                and not _time_explicitly_dropped):
            raise RuntimeError(
                f"Time axis unexpectedly dropped in MV2sum: "
                f"id={a.id!r}, original_axes={[ax.id if ax else None for ax in a._axes]}, "
                f"axis={axis!r}, drop={drop}, result_shape={result.shape}"
            )
        return CDATVariable(
            result,
            axes=new_axes,
            grid=None if len(new_axes) < 2 else a._grid,
            id=a.id,
            attributes=dict(a._attributes),
        )
    return result
def MV2take(a, indices, axis=0):
    raw = _mv(a)
    taken_data = np.take(np.asarray(raw), indices, axis=axis)
    taken_mask = np.take(ma.getmaskarray(raw), indices, axis=axis)
    result = ma.array(taken_data, mask=taken_mask)
    if isinstance(a, CDATVariable):
        new_axes = list(a._axes)
        ax_int = axis if isinstance(axis, (int, np.integer)) else 0
        if 0 <= ax_int < len(new_axes) and new_axes[ax_int] is not None:
            old_ax = new_axes[ax_int]
            idx = np.asarray(indices)
            new_ax = _Axis(
                old_ax.id, old_ax._values[idx],
                units=old_ax.units,
                attributes=dict(old_ax._attributes),
                axis_type=old_ax.axis
            )
            new_ax.calendar = old_ax.calendar
            new_axes[ax_int] = new_ax
        return CDATVariable(result, axes=new_axes, grid=a._grid,
                            id=a.id, attributes=dict(a._attributes))
    return result
def MV2where(condition, x, y):
    return _mv_wrap(ma.where(np.asarray(condition, dtype=bool), _mv(x), _mv(y)), x)
def MV2zeros(shape):         return CDATVariable(ma.zeros(shape), id="")

# Internal helpers for MV2 aliases
def _mv(x):
    """Extract numpy.ma array from CDATVariable or return as-is."""
    if isinstance(x, CDATVariable):
        return x._data
    if isinstance(x, ma.MaskedArray):
        return x
    return np.asarray(x)

def _mv_wrap(result, template):
    """Re-wrap a numpy.ma result in a CDATVariable if template was one."""
    if isinstance(template, CDATVariable):
        return CDATVariable(result, axes=list(template._axes),
                            grid=template._grid, id=template.id,
                            attributes=dict(template._attributes))
    return result

def _nan_majority_mask(data, axis):
    """
    Return True where a reduction should be masked because at least half of the
    contributing cells are missing.
    """
    arr = ma.masked_invalid(ma.asarray(data))
    mask = ma.getmaskarray(arr)
    axis = tuple(axis) if isinstance(axis, (tuple, list)) else axis
    invalid = np.sum(mask, axis=axis)
    total = np.sum(np.ones(mask.shape, dtype=np.int64), axis=axis)
    return np.asarray(invalid * 2 >= total)

def _apply_nan_majority_policy(result, source, axis):
    """Mask *result* wherever *source* has >=50% invalid values along *axis*."""
    src = _to_cdat(source)
    ax = _axis_to_int(src, axis)
    bad = _nan_majority_mask(_mv(src), ax)
    res_data = ma.array(_mv(result), copy=True)
    res_data = ma.array(
        res_data,
        mask=ma.getmaskarray(res_data) | np.broadcast_to(bad, res_data.shape),
    )
    if isinstance(result, CDATVariable):
        result = result.copy()
        result._data = res_data
        return result
    return res_data

def _valid_weight_sum(source, weights, axis):
    """Sum weights over *axis*, excluding cells where *source* is missing."""
    src = _to_cdat(source)
    weights = _to_cdat(weights)
    data = ma.masked_invalid(_mv(src))
    area = ma.masked_invalid(_mv(weights))
    ax = _axis_to_int(src, axis)

    if area.shape != data.shape:
        shape = [1] * data.ndim
        for ii, area_axis in enumerate(weights.getAxisList()):
            axis_type = getattr(area_axis, "axis", None)
            if axis_type == "Y":
                shape[_axis_to_int(src, "y")] = area.shape[ii]
            elif axis_type == "X":
                shape[_axis_to_int(src, "x")] = area.shape[ii]
            elif axis_type == "T":
                shape[_axis_to_int(src, "t")] = area.shape[ii]
            elif ii < data.ndim and area.shape[ii] == data.shape[ii]:
                shape[ii] = area.shape[ii]
        area = area.reshape(shape)

    area = ma.array(
        np.broadcast_to(area, data.shape),
        mask=np.broadcast_to(ma.getmaskarray(area), data.shape),
    )
    area = ma.array(area, mask=ma.getmaskarray(area) | ma.getmaskarray(data))
    return ma.sum(area, axis=ax)

def _to_cdat(x):
    """
    Ensure *x* is a CDATVariable.

    If it is already a CDATVariable it is returned unchanged.
    If it is a plain numpy.ma.MaskedArray (or any ndarray) it is wrapped in a
    CDATVariable so that CDAT-style methods like .getTime() work on it.
    Any existing axes/metadata on a CDATVariable are preserved.
    """
    if isinstance(x, CDATVariable):
        return x
    raw = _mv(x)
    return CDATVariable(raw, id="")

def _get_lat_weights(var, axis=None):
    """Return area weights broadcast to *var*'s shape, or None.

    Priority:
    1. ``cell_area`` attribute on *var* — exact area weights for stretched /
        regionally-refined (RRM) grids (E3SM, MPAS-A regular outputs).
    2. Cosine of latitude — standard approximation for regular grids.
    """
    var = _to_cdat(var)

    var_attrs = getattr(var, "attributes", {}) or {}
    var_label = (
        getattr(var, "id", None)
        or getattr(var, "name", None)
        or var_attrs.get("id", None)
        or var_attrs.get("short_name", None)
        or var_attrs.get("standard_name", None)
        or var_attrs.get("long_name", None)
        or "unknown"
    )

    # Priority 1: explicit cell_area (exact for stretched / RRM grids)
    cell_area = getattr(var, 'cell_area', None)
    if cell_area is not None:
        try:
            data = _mv(var)
            w = ma.masked_invalid(np.asarray(cell_area, dtype=float))
            if w.ndim > data.ndim:
                raise ValueError(
                    f"cell_area has {w.ndim} dims but data has only {data.ndim}; "
                    "cannot broadcast."
                )
            # Broadcast to data shape (handles time-invariant (y,x) areas)
            wbc = ma.array(np.broadcast_to(w, data.shape), copy=False)
            # Combine weight mask with data mask so excluded data points
            # are also excluded from the weight denominator.
            combined_mask = ma.getmaskarray(wbc) | ma.getmaskarray(data)
            wbc = ma.array(wbc, mask=combined_mask)
            return wbc

        except Exception as _cell_area_exc:
            warnings.warn(
                f"_get_lat_weights: cell_area weight computation failed for "
                f"variable {var_label!r} ({_cell_area_exc}); "
                "falling back to cosine-latitude weights. "
                "Results may be incorrect for stretched/RRM grids.",
                stacklevel=2,
            )

    # Priority 2: cos(lat) — standard approximation for regular rectilinear
    # grids.  This is accurate only when grid cells have uniform zonal width;
    # for stretched, RRM, or curvilinear grids the caller should attach
    # cell_area to avoid biased spatial averages.
    warnings.warn(
        f"_get_lat_weights: no cell_area found on variable {var_label!r}; "
        "falling back to cosine-latitude weights. "
        "Attach var.cell_area, such as areacella for atmospheric fields or areacello "
        "for ocean fields, for accurate spatial averages on non-uniform grids "
        "(E3SM-RRM, MPAS, stretched CMIP grids).",
        stacklevel=3,
    )

    lat = var.getLatitude()
    if lat is None:
        raise ValueError(
            f"_get_lat_weights: latitude axis not found for variable "
            f"{var_label!r}. Cannot compute spatial weights. "
            "Ensure the CDATVariable has a latitude axis with CF axis='Y' metadata."
        )

    lat_vals = np.asarray(lat[:], dtype=float)
    w = np.cos(np.deg2rad(lat_vals))
    w = ma.masked_invalid(w)

    lat_axis = _axis_to_int(var, "y")
    if lat_axis is None:
        raise ValueError(
            f"_get_lat_weights: latitude axis index cannot be determined for variable "
            f"{var_label!r}. Attach CF axis/standard_name metadata."
        )

    shape = [1] * var.ndim
    shape[lat_axis] = len(w)
    return w.reshape(shape)

def _weighted_spatial_average(tab, axes=("Y", "X")):
    """
    Compute a latitude-weighted spatial average using the ENSO_metrics
    compatibility layer.

    This helper preserves the subset of legacy ``cdutil.averager`` behavior
    needed by ENSO_metrics while using ``numpy.ma`` operations and CDAT-like
    compatibility metadata.

    Parameters
    ----------
    tab : CDATVariable or masked-array-like
        Input field with latitude/longitude axes when spatial weighting is
        requested.

    axes : tuple of str, optional
        CDAT-style axis-type strings to reduce over, for example ``("Y", "X")``,
        ``("Y",)``, or ``("X",)``.
        default value = ``("Y", "X")``

    Returns
    -------
    numpy.ma.MaskedArray
        Spatially averaged values with the requested dimensions collapsed.
    """
    tab = _to_cdat(tab)
    data = _mv(tab)
    do_y = "Y" in axes
    do_x = "X" in axes
    lat_w = _get_lat_weights(tab) if do_y else None
    # When latitude weighting is needed but weights are unavailable, raise
    # immediately — a silent unweighted fallback produces scientifically
    # incorrect ENSO metrics.
    if lat_w is None and do_y:
        raise RuntimeError(
            "_weighted_spatial_average: latitude weights missing — cannot perform "
            f"weighted meridional average for variable {getattr(tab, 'id', '?')!r}. "
            "Ensure the CDATVariable has a latitude axis with CF metadata."
        )
    if not do_y:
        # Zonal-only average: no latitude weighting needed, unweighted is correct.
        if do_x:
            ax_int = _axis_to_int(tab, "x")
            if ax_int is None:
                raise ValueError(
                    f"_weighted_spatial_average: cannot determine X averaging axis "
                    f"for axes={axes!r}; attach CF axis metadata."
                )
            result = ma.mean(data, axis=ax_int)
            result = ma.array(
                result,
                mask=ma.getmaskarray(result) | _nan_majority_mask(data, ax_int),
            )
            return result
        return data  # nothing to reduce
    mask = ma.getmaskarray(data)
    wm = ma.array(np.broadcast_to(lat_w, data.shape), mask=mask)
    # Re-apply the same mask to data so that data*wm uses identical masking
    weighted = ma.array(data, mask=mask) * wm
    if do_x and do_y:
        ax_int = _axis_to_int(tab, "xy")
    elif do_y:
        ax_int = _axis_to_int(tab, "y")
    else:  # do_x only
        ax_int = _axis_to_int(tab, "x")
    if ax_int is None:
        raise ValueError(
            f"_weighted_spatial_average: cannot determine averaging axis "
            f"for axes={axes!r}; attach CF axis metadata or pass an explicit integer."
        )
    num = ma.sum(weighted, axis=ax_int)
    den = ma.sum(wm, axis=ax_int)
    # Guard against fully-masked latitude bands: use clean masked division
    den_safe = ma.where(den == 0, ma.masked, den)
    result = num / den_safe
    result = ma.array(
        result,
        mask=ma.getmaskarray(result) | _nan_majority_mask(data, ax_int),
    )
    return result


def _is_unstructured_grid(da):
    """Return True when *da* looks like an unstructured (MPAS/ICON/SE) grid."""
    _unstructured_dims = {"ncol", "ncell", "ncells", "nvertices", "nedges"}
    if hasattr(da, 'dims'):             # xarray DataArray / Dataset
        dims = {str(d).lower() for d in da.dims}
    elif isinstance(da, CDATVariable):
        dims = {ax.id.lower() for ax in da._axes if ax is not None}
    else:
        return False
    return bool(dims & _unstructured_dims)


def _coord_values_degrees(coord):
    """Return coordinate values in degrees, accepting common MPAS radian coords."""
    vals = np.asarray(coord.values, dtype=float)
    units = str(coord.attrs.get("units", "")).lower()
    if units in {"radian", "radians", "rad"}:
        vals = np.rad2deg(vals)
    return vals


def _unstructured_lat_lon_coords(da: xr.DataArray, ds: xr.Dataset | None = None):
    """Return MPAS-like 1-D lat/lon coordinate arrays and their shared cell dim."""
    candidates = [
        ("latCell", "lonCell"),
        ("lat_cell", "lon_cell"),
        ("latitude", "longitude"),
        ("lat", "lon"),
    ]
    for lat_name, lon_name in candidates:
        lat = da.coords.get(lat_name)
        lon = da.coords.get(lon_name)
        if (lat is None or lon is None) and ds is not None:
            lat = ds[lat_name] if lat_name in ds else lat
            lon = ds[lon_name] if lon_name in ds else lon
        if lat is None or lon is None or lat.ndim != 1 or lon.ndim != 1:
            continue
        if lat.dims != lon.dims:
            continue
        cell_dim = lat.dims[0]
        if cell_dim in da.dims:
            return lat_name, lon_name, cell_dim, lat, lon
    return None, None, None, None, None


def _subset_unstructured_lat_lon(da: xr.DataArray, ds: xr.Dataset | None = None, latitude=None, longitude=None):
    """Subset MPAS-like unstructured cell data using 1-D cell lat/lon coords."""
    lat_name, lon_name, cell_dim, lat, lon = _unstructured_lat_lon_coords(da, ds)
    if cell_dim is None:
        return da, False
    lat_vals = _coord_values_degrees(lat)
    lon_vals = np.mod(_coord_values_degrees(lon), 360.0)
    mask = np.ones(lat_vals.shape, dtype=bool)
    if latitude is not None:
        lo, hi = sorted([float(v) for v in latitude])
        mask &= (lat_vals >= lo) & (lat_vals <= hi)
    if longitude is not None:
        mask &= _longitude_mask(lon_vals, longitude)
    keep = xr.DataArray(mask, dims=cell_dim)
    out = da.where(keep, drop=True)
    out = out.assign_coords(
        {
            lat_name: xr.DataArray(
                lat_vals[mask],
                dims=cell_dim,
                attrs={"axis": "Y", "standard_name": "latitude", "units": "degrees_north"},
            ),
            lon_name: xr.DataArray(
                lon_vals[mask],
                dims=cell_dim,
                attrs={"axis": "X", "standard_name": "longitude", "units": "degrees_east"},
            ),
        }
    )
    return out, True


def _axis_to_int(arr, axis):
    """
    Convert a CDAT-style axis spec to an integer or tuple of integers
    that numpy.ma can accept.

    Handles:
    - int / None / tuple  → returned as-is
    - "0", "1"            → int
    - "01", "10", "12"   → tuple of ints (one digit per axis)
    - "t" / "time"        → time axis index
    - "y" / "lat"         → latitude axis index
    - "x" / "lon"         → longitude axis index
    - "z" / "lev"         → level axis index
    - "xy", "yx"          → (lat_idx, lon_idx)

    Uses CDATVariable axis metadata when available, falling back to
    common layout heuristics (t=0, y=ndim-2, x=ndim-1).
    """
    if axis is None or isinstance(axis, (int, np.integer)):
        return axis
    # Coerce plain Python float or numpy floating to int to avoid
    # "'float' object cannot be interpreted as an integer" from numpy.
    if isinstance(axis, (float, np.floating)):
        return int(round(float(axis)))
    if isinstance(axis, (tuple, list)):
        return tuple(int(a) for a in axis)
    axis_s = str(axis).strip()
    # Pure-digit string: "0"→0, "01"→(0,1), "10"→(1,0)
    if axis_s.isdigit():
        digits = tuple(int(c) for c in axis_s)
        return digits[0] if len(digits) == 1 else digits
    ndim = getattr(arr, 'ndim', 1)
    # Build named-axis index map from CDATVariable metadata when available
    ax_map = {}
    src = arr if isinstance(arr, CDATVariable) else None
    if src is not None:
        for i, ax in enumerate(src._axes):
            if ax is not None and ax.axis in ("T", "Y", "X", "Z"):
                ax_map[ax.axis] = i
    # Limited positional fallback heuristics.
    #
    # Applied only when:
    #   - no CF axis metadata exists at all, and
    #   - ndim <= 2.
    #
    # This supports simple legacy arrays while avoiding dangerous silent
    # axis inference for higher-dimensional CMIP/E3SM datasets.
    #
    # High-dimensional arrays without identifiable axis metadata are unsafe
    # to interpret positionally. Respect STRICT_DIM_GUESS so legacy workflows
    # can demote this to a warning when needed.
    if not ax_map and ndim > 3:
        msg = (
            f"Cannot infer axes for {ndim}-D array without CF axis metadata; "
            "attach axis types or pass an explicit integer axis."
        )

        if STRICT_DIM_GUESS:
            raise ValueError(msg)

        warnings.warn(msg, stacklevel=2)

        # Explicitly stop positional guessing for high-dimensional arrays.
        return axis

    # Positional heuristics: only applied when *no* CF axis metadata was found
    # at all (ax_map is empty after the loop above).  Mixing positional guesses
    # with partial real metadata can silently reduce the wrong dimension on
    # non-standard layouts (CMIP ensemble dim, E3SM extra dims, etc.).
    if not ax_map and ndim <= 2:
        if ndim == 2:
            ax_map.setdefault("Y", 0)
            ax_map.setdefault("X", 1)
        elif ndim == 1:
            ax_map.setdefault("X", 0)

    axis_l = axis_s.lower()
    if axis_l in ("t", "time"):
        if "T" not in ax_map:
            raise ValueError(
                "Cannot determine time axis from metadata. "
                "Attach CF axis='T' or standard_name='time'."
            )
        return ax_map["T"]
    if axis_l in ("y", "lat", "latitude"):
        if "Y" not in ax_map:
            raise ValueError(
                "Cannot determine latitude axis from metadata."
            )
        return ax_map["Y"]
    if axis_l in ("x", "lon", "longitude"):
        if "X" not in ax_map:
            raise ValueError(
                "Cannot determine longitude axis from metadata."
            )
        return ax_map["X"]
    if axis_l in ("z", "lev", "level", "depth"):
        if "Z" not in ax_map:
            raise ValueError("Cannot determine vertical axis from metadata.")
        return ax_map["Z"]
    if axis_l in ("xy", "yx"):
        if "Y" not in ax_map or "X" not in ax_map:
            raise ValueError(
                "Cannot reliably determine lat/lon axes from metadata; "
                "ensure the variable has explicit 'Y'/'X' axis types."
            )
        return (ax_map["Y"], ax_map["X"])
    # Multi-character: scan char by char
    indices, seen = [], set()
    for c in axis_l:
        idx = None
        if c == "t":
            idx = ax_map.get("T")
        elif c == "y":
            idx = ax_map.get("Y")
        elif c == "x":
            idx = ax_map.get("X")
        elif c == "z":
            idx = ax_map.get("Z")
        if idx is not None and idx not in seen:
            indices.append(idx)
            seen.add(idx)
    if indices:
        return indices[0] if len(indices) == 1 else tuple(indices)
    return axis  # unknown string → pass through; numpy will raise a clear error


# ----------------------------------------------------------------------------------
# Legacy cdutil-style season helpers implemented with xarray/xcdat-compatible logic
# ----------------------------------------------------------------------------------
class _SeasonHelper:
    """Thin wrapper that mimics cdutil.JAN, cdutil.DJF, etc."""
    def __init__(self, months: list):
        self._months = months

    def __call__(self, tab):
        return _seasonal_mean(tab, self._months, compute_anom=False)

    def departures(self, tab):
        return _seasonal_mean(tab, self._months, compute_anom=True)


# ---------------------------------------------------------------------------
# Time-axis safety helpers (single authoritative block)
# ---------------------------------------------------------------------------
def _require_time_axis(var, context=''):
    """Return the time axis of *var*, raising ValueError with a clear message if absent."""
    ax = var.getTime()
    if ax is None:
        label = getattr(var, 'id', '') or context or 'variable'
        raise ValueError(
            f"No time axis found on '{label}'. "
            "Ensure the input data has a recognised time dimension."
        )
    return ax


def _is_datetime_like_time(coord):
    try:
        if np.issubdtype(coord.dtype, np.datetime64):
            return True
    except Exception:
        pass
    try:
        vals = coord.values
        return len(vals) > 0 and hasattr(vals[0], "year")
    except Exception:
        return False


def _ensure_time_encoding(ds: xr.Dataset, path: str = "") -> xr.Dataset:
    """
    STRICT: Ensure time is decoded. Never silently pass numeric time.
    """
    if "time" not in ds.coords:
        return ds

    tc = ds["time"]
    cal = tc.attrs.get("calendar") or tc.encoding.get("calendar") or "standard"

    if _is_datetime_like_time(tc):
        ds["time"].encoding.setdefault("calendar", cal)
        return ds

    # numeric → must decode
    try:
        is_numeric = np.issubdtype(tc.dtype, np.number)
    except Exception:
        is_numeric = False

    if is_numeric:
        units = tc.attrs.get("units") or tc.encoding.get("units") or ""

        if "since" not in units:
            # Non-standard or missing time units — cannot decode; fail early
            # so the caller gets a clear diagnostic rather than a silent
            # mis-identification later in the pipeline.
            raise RuntimeError(
                f"Time coordinate has numeric dtype but unrecognisable units "
                f"{units!r} (expected CF 'X since Y' format). "
                f"File: {path!r}. Fix the source file's time:units attribute."
            )

        try:
            coder = _get_time_coder()
            decoded = coder.decode(
                xr.Variable("time", tc.values, {"units": units, "calendar": cal}),
                name="time",
            )

            ds = ds.assign_coords(
                time=xr.DataArray(
                    decoded.values,
                    dims=tc.dims,
                    attrs={**tc.attrs, "calendar": cal},
                )
            )

        except Exception as e:
            raise RuntimeError(
                f"Time decoding failed for {path}: {e}"
            ) from e

    if not _is_datetime_like_time(ds["time"]):
        raise RuntimeError(
            f"Time still not decoded: dtype={ds['time'].dtype}, file={path}"
        )

    ds["time"].encoding.setdefault("calendar", cal)
    return ds


def _has_time_axis(var):
    return _to_cdat(var).getTime() is not None


def _get_time_axis_index(var, context=""):
    var = _to_cdat(var)
    time_ax = _require_time_axis(var, context)
    for i, ax in enumerate(var.getAxisList()):
        if ax is time_ax or getattr(ax, "axis", None) == "T":
            return i
    raise ValueError(f"Cannot determine time axis index in {context!r}")


def _get_component_time(var, context=""):
    var = _to_cdat(var)
    comp = _require_time_axis(var, context).asComponentTime()
    if len(comp) == 0:
        raise ValueError(f"Empty time axis in {context}")
    return comp


def _safe_time_bounds_for_debug(tab):
    try:
        comp = _get_component_time(tab, "TimeBounds")
        return str(comp[0]), str(comp[-1])
    except Exception:
        return None, None


def _seasonal_mean(tab, month_list, compute_anom=False):
    """
    Compute seasonal mean or departures using xcdat on the given CDATVariable.
    Returns a CDATVariable.
    """
    tab = _to_cdat(tab)
    _require_time_axis(tab, "_seasonal_mean")

    da = cdat_to_da(tab)
    varname = tab.id or "var"
    ds = da.to_dataset(name=varname)
    ds = _ensure_time_encoding(ds, path=f"in-memory:{varname}")

    try:
        if "time" in ds.coords and _is_datetime_like_time(ds["time"]):
            ds = ds.bounds.add_missing_bounds(axes=["T"])
    except Exception:
        pass

    _month_abbr = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ]
    month_names = [_month_abbr[m - 1] for m in month_list]
    season_cfg = {"custom_seasons": [month_names]}

    try:
        if compute_anom:
            result_ds = ds.temporal.departures(
                varname,
                freq="season",
                weighted=True,
                season_config=season_cfg,
            )
        else:
            result_ds = ds.temporal.group_average(
                varname,
                freq="season",
                weighted=True,
                season_config=season_cfg,
            )
        return _finalize_cdat(result_ds[varname], varname=tab.id,
                              context="_seasonal_mean", require_time=True)

    except Exception:
        comp = _get_component_time(tab, "_seasonal_mean fallback")

        months_arr = np.array([t.month for t in comp])
        years_arr = np.array([t.year for t in comp])
        mask_months = np.isin(months_arr, month_list)
        idx = np.where(mask_months)[0]

        raw = _mv(tab)[idx]
        years_sel = years_arr[idx]
        months_sel = months_arr[idx]

        unique_years = sorted(set(years_sel))
        season_data = []
        season_years = []

        for yr in unique_years:
            yr_mask = years_sel == yr
            present = sorted(months_sel[yr_mask].tolist())
            # Strict match: require all season months to be present.
            # Partial seasons (e.g. DJF with only JF) are dropped rather
            # than included as biased means — consistent with CDAT behavior.
            if present == sorted(month_list):
                season_data.append(ma.mean(raw[yr_mask], axis=0))
                season_years.append(yr)

        if not season_data:
            return tab.copy()

        result_raw = ma.array(season_data)

        if compute_anom:
            result_raw = result_raw - ma.mean(result_raw, axis=0, keepdims=True)

        time_new = _Axis(
            "time",
            np.array(season_years, dtype="int32"),
            units="years since 0001-01-01",
            axis_type="T",
        )

        new_axes = [time_new] + (
            tab.getAxisList()[1:] if len(tab.shape) > 1 else []
        )

        return _finalize_existing_cdat(
            CDATVariable(
                result_raw,
                axes=new_axes,
                grid=tab.getGrid(),
                id=tab.id,
                attributes=dict(tab.attributes),
            ),
            context="_seasonal_mean:fallback",
            require_time=True,
        )


# Build the sea_dict equivalent (populated after _SeasonHelper is defined)
_MONTH_MAP = {
    "JAN": [1], "FEB": [2], "MAR": [3], "APR": [4], "MAY": [5], "JUN": [6],
    "JUL": [7], "AUG": [8], "SEP": [9], "OCT": [10], "NOV": [11], "DEC": [12],
    "JF": [1, 2], "FM": [2, 3], "MA": [3, 4], "AM": [4, 5], "MJ": [5, 6],
    "JJ": [6, 7], "JA": [7, 8], "AS": [8, 9], "SO": [9, 10], "ON": [10, 11],
    "ND": [11, 12], "DJ": [12, 1],
    "JFM": [1, 2, 3], "FMA": [2, 3, 4], "MAM": [3, 4, 5], "AMJ": [4, 5, 6],
    "MJJ": [5, 6, 7], "JJA": [6, 7, 8], "JAS": [7, 8, 9], "ASO": [8, 9, 10],
    "SON": [9, 10, 11], "OND": [10, 11, 12], "NDJ": [11, 12, 1], "DJF": [12, 1, 2],
    "JFMA": [1, 2, 3, 4], "FMAM": [2, 3, 4, 5], "MAMJ": [3, 4, 5, 6],
    "AMJJ": [4, 5, 6, 7], "MJJA": [5, 6, 7, 8], "JJAS": [6, 7, 8, 9],
    "JASO": [7, 8, 9, 10], "ASON": [8, 9, 10, 11], "SOND": [9, 10, 11, 12],
    "ONDJ": [10, 11, 12, 1], "NDJF": [11, 12, 1, 2], "DJFM": [12, 1, 2, 3],
}
sea_dict = {k: _SeasonHelper(v) for k, v in _MONTH_MAP.items()}

# ---------------------------------------------------------------------------
# Grid consistency guard
# ---------------------------------------------------------------------------
def check_grid_consistency(a, b, context: str = "", regrid_to: str = "b") -> tuple:
    """Check whether two fields are on compatible horizontal grids.

    The function verifies that both fields have matching latitude/longitude
    dimensions. If coordinate sizes match but latitude or longitude coordinate
    values differ, it emits a warning so callers can detect same-shape but
    different-grid cases.

    Parameters
    ----------
    a, b : CDATVariable or array-like
        Fields to compare.
    context : str, optional
        Label used in error/warning messages, for example the calling function
        name.
    regrid_to : str, optional
        Retained for API compatibility. This function does not perform
        automatic regridding; callers should use ``Regrid()`` before metric
        calculation when grids differ.

    Returns
    -------
    tuple
        ``(a, b)`` unchanged.

    Raises
    ------
    ValueError
        If latitude/longitude sizes differ.
    """
    a_cdat = _to_cdat(a) if not isinstance(a, CDATVariable) else a
    b_cdat = _to_cdat(b) if not isinstance(b, CDATVariable) else b
    a_lat = a_cdat.getLatitude()
    b_lat = b_cdat.getLatitude()
    a_lon = a_cdat.getLongitude()
    b_lon = b_cdat.getLongitude()

    # Shape mismatch check (lat/lon sizes must match)
    a_lat_n = len(a_lat) if a_lat is not None else None
    a_lon_n = len(a_lon) if a_lon is not None else None
    b_lat_n = len(b_lat) if b_lat is not None else None
    b_lon_n = len(b_lon) if b_lon is not None else None

    if (a_lat_n is not None and b_lat_n is not None and a_lat_n != b_lat_n) or \
        (a_lon_n is not None and b_lon_n is not None and a_lon_n != b_lon_n):
        raise ValueError(
            f"{context}: spatial grid mismatch — "
            f"a has ({a_lat_n}, {a_lon_n}) lat/lon points, "
            f"b has ({b_lat_n}, {b_lon_n}). "
            "Regrid both fields to a common grid before computing metrics "
            "(use Regrid() with a target grid or set newgrid_name='generic_1x1deg')."
        )

    # Value-level mismatch: same size but different coordinate values
    if (a_lat is not None and b_lat is not None and a_lat_n == b_lat_n):
        if not np.allclose(
            np.asarray(a_lat[:], dtype=float),
            np.asarray(b_lat[:], dtype=float), atol=1e-4
            ):
            warnings.warn(
                f"{context}: latitude coordinate values differ between a and b "
                "despite matching sizes — verify both are on the same grid.",
                stacklevel=3,
            )

    if (a_lon is not None and b_lon is not None and a_lon_n == b_lon_n):
        if not np.allclose(
            np.asarray(a_lon[:], dtype=float),
            np.asarray(b_lon[:], dtype=float),
            atol=1e-4,
        ):
            warnings.warn(
                f"{context}: longitude coordinate values differ between a and b "
                "despite matching sizes — verify both are on the same grid.",
                stacklevel=3,
            )
    return a, b


# ---------------------------------------------------------------------------
# genutil.statistics replacements
# ---------------------------------------------------------------------------
def GENUTILcorrelation(a, b, weights=None, axis=0, centered=1, biased=1):
    a, b = check_grid_consistency(a, b, context="GENUTILcorrelation")
    x = ma.masked_invalid(_mv(a))
    y = ma.masked_invalid(_mv(b))
    axis = _axis_to_int(a, axis)
    if (
        x.shape == y.shape
        and np.array_equal(ma.getmaskarray(x), ma.getmaskarray(y))
        and np.allclose(
            x.astype(float).filled(np.nan),
            y.astype(float).filled(np.nan),
            equal_nan=True,
            atol=1e-12,
            rtol=1e-12,
        )
    ):
        valid_count = x.count(axis=axis)
        result = ma.array(np.ones(np.shape(valid_count), dtype=float))
        result.mask = (np.asarray(valid_count) < 2) | _nan_majority_mask(x, axis)
        return result
    if centered:
        x = x - ma.mean(x, axis=axis, keepdims=True)
        y = y - ma.mean(y, axis=axis, keepdims=True)
    if isinstance(weights, str) and weights.lower() == "weighted":
        w = _get_lat_weights(a, axis=axis)
    elif weights is not None:
        w = _mv(weights)
    else:
        w = None
    if w is not None:
        # Broadcast weights to data shape; use manual sum/sum for multi-axis support
        wbc = np.broadcast_to(np.asarray(w), x.shape)
        # Single unified mask: exclude points where *either* x or y is missing
        # so that cov, varx, and vary all share the same effective sample.
        wm   = ma.array(wbc, mask=ma.getmaskarray(x) | ma.getmaskarray(y))
        den  = ma.sum(wm, axis=axis)
        den_safe = ma.where(den == 0, ma.masked, den)
        cov  = ma.sum(x * y * wm, axis=axis) / den_safe
        varx = ma.sum(x ** 2 * wm, axis=axis) / den_safe
        vary = ma.sum(y ** 2 * wm, axis=axis) / den_safe
    else:
        n   = x.count(axis=axis) if biased else (x.count(axis=axis) - 1)
        if np.any(np.asarray(n) == 0):
            return ma.masked
        cov  = ma.sum(x * y,  axis=axis) / n
        varx = ma.sum(x ** 2, axis=axis) / n
        vary = ma.sum(y ** 2, axis=axis) / n
    return cov / ma.sqrt(varx * vary)

def GENUTILrms(a, b, weights=None, axis=0, centered=0, biased=1):
    a, b = check_grid_consistency(a, b, context="GENUTILrms")
    x = ma.masked_invalid(_mv(a))
    y = ma.masked_invalid(_mv(b))
    axis = _axis_to_int(a, axis)
    diff = x - y
    if centered:
        diff = diff - ma.mean(diff, axis=axis, keepdims=True)
    if isinstance(weights, str) and weights.lower() == "weighted":
        w = _get_lat_weights(a, axis=axis)
    elif weights is not None:
        w = _mv(weights)
    else:
        w = None
    if w is not None:
        # Broadcast weights to data shape; use manual sum/sum for multi-axis support
        wbc = np.broadcast_to(np.asarray(w), diff.shape)
        # Use union mask of x and y so the effective sample size matches
        # GENUTILcorrelation — avoids inconsistency when masking is asymmetric.
        wm  = ma.array(wbc, mask=ma.getmaskarray(x) | ma.getmaskarray(y))
        den = ma.sum(wm, axis=axis)
        den_safe = ma.where(den == 0, ma.masked, den)
        return ma.sqrt(ma.sum(diff ** 2 * wm, axis=axis) / den_safe)
    n = diff.count(axis=axis) if biased else (diff.count(axis=axis) - 1)
    if np.any(np.asarray(n) == 0):
        return ma.masked
    return ma.sqrt(ma.sum(diff ** 2, axis=axis) / n)

def GENUTILstd(a, weights=None, axis=0, centered=1, biased=1):
    x = ma.masked_invalid(_mv(a))
    axis = _axis_to_int(a, axis)
    if centered:
        x = x - ma.mean(x, axis=axis, keepdims=True)
    if isinstance(weights, str) and weights.lower() == "weighted":
        w = _get_lat_weights(a, axis=axis)
    elif weights is not None:
        w = _mv(weights)
    else:
        w = None
    if w is not None:
        # Broadcast weights to data shape; use manual sum/sum for multi-axis support
        wbc = np.broadcast_to(np.asarray(w), x.shape)
        wm  = ma.array(wbc, mask=ma.getmaskarray(x))
        den = ma.sum(wm, axis=axis)
        den_safe = ma.where(den == 0, ma.masked, den)
        result = ma.sqrt(ma.sum(x ** 2 * wm, axis=axis) / den_safe)
    else:
        ddof = 0 if biased else 1
        result = ma.std(x, axis=axis, ddof=ddof)
    if isinstance(a, CDATVariable):
        drop = {axis} if isinstance(axis, (int, np.integer)) else set(axis)
        new_axes = [ax for i, ax in enumerate(a._axes) if i not in drop]
        return CDATVariable(result, axes=new_axes, grid=a._grid,
                            id=a.id, attributes=dict(a._attributes))
    return result

def GENUTILlinearregression(y, x=None, error=1, nointercept=None):
    """
    Flattened replacement for genutil.statistics.linearregression.

    This helper collapses all dimensions with ravel() before regression and
    therefore returns one global scalar slope/intercept pair. It is appropriate
    for scalar/global regressions, but it should not be used when the caller
    expects axis-preserving regression fields such as:

        (time, lon)        -> (lon)
        (year, month, lon) -> (month, lon)

    For axis-preserving regression, use CustomLinearRegression(), which applies
    CustomLinearRegression1d pointwise along the first axis.
    """
    yy = np.ma.masked_invalid(np.ma.asarray(_mv(y)).ravel())

    if x is None:
        xx = np.ma.array(np.arange(yy.size, dtype=float), mask=np.ma.getmaskarray(yy))
    else:
        xx = np.ma.masked_invalid(np.ma.asarray(_mv(x)).ravel())

        if xx.size != yy.size:
            min_len = min(xx.size, yy.size)
            xx = xx[:min_len]
            yy = yy[:min_len]

    # Joint mask: keep only pairs where both x and y are valid
    joint_mask = np.ma.getmaskarray(xx) | np.ma.getmaskarray(yy)
    xf = np.asarray(np.ma.array(xx, mask=joint_mask).compressed(), dtype=float)
    yf = np.asarray(np.ma.array(yy, mask=joint_mask).compressed(), dtype=float)

    if len(xf) < 2:
        slope_int = np.array([[np.nan, np.nan]])
        stderr = np.array([[np.nan, np.nan]])
        return (slope_int, stderr) if error else slope_int

    same_series = np.allclose(xf, yf, equal_nan=True, atol=1e-12, rtol=1e-12)

    if nointercept == 1:
        denom = np.dot(xf, xf)
        if denom == 0:
            slope = np.nan
            se = np.nan
        elif same_series:
            slope = 1.0
            se = 0.0
        else:
            slope = float(np.dot(xf, yf) / denom)
            resid = yf - slope * xf
            se = float(
                np.sqrt(np.sum(resid ** 2) / max(len(xf) - 1, 1)) /
                np.sqrt(denom)
            )

        slope_int = np.array([[slope, 0.0]])
        stderr = np.array([[se, 0.0]])

    else:
        if len(np.unique(xf)) < 2:
            slope_int = np.array([[np.nan, np.nan]])
            stderr = np.array([[np.nan, np.nan]])
        elif same_series:
            slope_int = np.array([[1.0, 0.0]])
            stderr = np.array([[0.0, 0.0]])
        else:
            res = _linregress(xf, yf)
            slope_int = np.array([[res.slope, res.intercept]])
            stderr = np.array([[res.stderr, res.intercept_stderr]])

    if error:
        return slope_int, stderr
    return slope_int

def _add_cf_units_to_ds(ds: xr.Dataset) -> xr.Dataset:
    """
    Add CF-standard units to lat/lon coordinates that are missing a 'units'
    attribute, inferred from value ranges.  This avoids xcdat warnings when
    opening files that lack coordinate units.
    """
    for name, coord in ds.coords.items():
        if "units" in coord.attrs:
            continue
        # Only infer units for 1-D coordinates with numeric values;
        # skip 2-D bounds, projected coords, and rotated-grid auxiliaries.
        if coord.ndim != 1:
            continue
        try:
            vals = coord.values.astype(float)
        except Exception:
            continue
        lo, hi = float(vals.min()), float(vals.max())
        # Latitude: strict [-90, 90]
        if -90.0 <= lo and hi <= 90.0 and name.lower() in (
                "lat", "latitude", "y"):
            ds[name].attrs["units"] = "degrees_north"
        # Longitude: [-180, 360]
        elif -180.0 <= lo and hi <= 360.0 and name.lower() in (
                "lon", "longitude", "x"):
            ds[name].attrs["units"] = "degrees_east"
    return ds


def _fix_leap_seconds_in_raw(ds: xr.Dataset) -> xr.Dataset:
    """
    Given a dataset opened with ``decode_times=False``, find any time
    coordinate or time_bnds values that encode a leap second (would decode
    to second=60 in cftime) and subtract 1/86400 day so they decode cleanly
    to second=59.  Must be called *before* any cftime decoding step.
    """
    if "time" not in ds.coords:
        return ds
    time_coord = ds.coords["time"]
    units = time_coord.attrs.get("units", "")
    calendar = time_coord.attrs.get("calendar", "standard") or "standard"
    if not units or "since" not in units:
        return ds

    # Fix the time coordinate and any associated bounds variable
    candidates = ["time"]
    bounds_name = time_coord.attrs.get("bounds", "")
    if bounds_name and bounds_name in ds:
        candidates.append(bounds_name)
    for alt in ("time_bnds", "time_bounds"):
        if alt in ds and alt not in candidates:
            candidates.append(alt)

    updated = {}
    for vname in candidates:
        da = ds[vname] if vname in ds else ds.coords.get(vname)
        if da is None:
            continue
        raw = np.asarray(da.values, dtype=float)
        flat = raw.ravel()
        # Short-circuit: leap seconds encode as fractional-day offsets whose
        # sub-day remainder > 59 s.  Skip the per-value scan when none exist.
        sub_day_seconds = (flat % 1.0) * 86400.0
        if not np.any(sub_day_seconds > 59.0):
            continue
        fixed = flat.copy()
        changed = False
        for i, v in enumerate(flat):
            try:
                dt = _cft.num2date(v, units, calendar)
                if getattr(dt, "second", 0) == 60:
                    fixed[i] = v - 1.0 / 86400
                    changed = True
            except ValueError:
                # cftime raises ValueError for second=60
                fixed[i] = v - 1.0 / 86400
                changed = True
        if changed:
            updated[vname] = fixed.reshape(raw.shape)

    if not updated:
        return ds

    ds = ds.copy()
    for vname, vals in updated.items():
        if vname in ds.coords:
            ds = ds.assign_coords(
                {vname: xr.DataArray(
                    vals, dims=ds[vname].dims,
                    attrs=ds[vname].attrs
                )}
            )
        else:
            ds[vname] = xr.DataArray(vals, dims=ds[vname].dims,
                                     attrs=ds[vname].attrs)
    return ds


def _sanitize_time_bound(t) -> str:
    """
    Convert a CDAT-style time bound string that may contain second=60
    (e.g. '2015-12-31 23:59:60.0') to a cftime/xarray-valid string
    ('2015-12-31 23:59:59').

    This is needed because downstream PCMDI code constructs string time
    bounds with second=60, and xarray's .sel() passes them directly to
    cftime which rejects second=60 with ValueError.
    """
    s = str(t).strip().replace("  ", " ")
    s = s.split(".")[0]        # drop sub-second part
    if s.endswith(":60"):
        s = s[:-3] + ":59"
    return s


def _to_sel_bound(raw_bound, t_dim: str, da: "xr.DataArray"):
    """
    Convert a time bound string to a cftime object whose type matches the
    coordinate in *da*, then pass that to xarray's ``.sel()``.

    **Why this matters**: when the time coordinate contains
    ``cftime.DatetimeNoLeap`` (or any other ``cftime.datetime`` sub-class)
    and the slice bounds are plain Python strings, xarray internally calls
    ``datetime.datetime.strptime(...)`` which returns a standard
    ``datetime.datetime`` object.  Comparing ``datetime.datetime`` against
    ``cftime.datetime`` raises a ``TypeError`` in Python 3; xarray and
    pandas catch this silently and report *zero* matching indices — i.e.
    an empty DataArray — even when the requested range spans the entire
    file.  Providing a ``cftime.datetime`` of the *correct calendar*
    prevents the cross-type comparison entirely.
    """
    s = _sanitize_time_bound(raw_bound)
    coord = da.coords.get(t_dim) if t_dim in da.coords else None
    if coord is None or len(coord) == 0:
        return s                        # nothing to match against; keep string
    first = coord.values[0]
    if not hasattr(first, "calendar"):
        return s                        # numpy datetime64 — strings work fine
    try:
        clean = s.replace("T", " ").strip()
        parts = clean.split()
        d = parts[0].split("-")
        tt = parts[1].split(":") if len(parts) > 1 else ["0", "0", "0"]
        return _cft.datetime(
            int(d[0]), int(d[1]), int(d[2]),
            int(tt[0]), int(tt[1]),
            min(int(float(tt[2])), 59) if len(tt) > 2 else 0,
            calendar=first.calendar,
        )
    except Exception:
        return s                        # parsing failed; fall back to string


# ---------------------------------------------------------------------------
# Data-construction / metadata finalization helpers
# ---------------------------------------------------------------------------
def _coord_is_datetime_like(da: xr.DataArray, dim: str) -> bool:
    """True if *dim* has a datetime-like coordinate on *da*."""
    if dim not in da.coords:
        return False
    return _is_datetime_like_time(da.coords[dim])

def _standardize_da_axes(da: xr.DataArray, *, context: str = "") -> xr.DataArray:
    """
    Attach CF-style axis metadata to every dimension before conversion to CDATVariable.

    Iterates each dimension directly using _dim_to_axis_type (CF attrs →
    standard_name → units → name heuristic) rather than calling _guess_dim
    (which searches for a dim given an axis type and emits "low-confidence
    fallback" warnings even when detection is successful).  This eliminates
    spurious warnings at the read boundary while still stamping axis='T/Y/X/Z'
    before da_to_cdat sees the coordinate.

    Additional fallbacks handled here (not in _dim_to_axis_type):
    - Datetime-like coord values → axis='T' (covers xcdat-decoded time coords
      that lack CF attrs).
    - Numeric coord with units="… since …" → axis='T' (unencoded numeric time).
    - Dim without a coordinate: creates a minimal arange coord with CF attrs so
      da_to_cdat can attach the correct axis type to the resulting _Axis object.
    """
    da = da.copy()

    for dim in da.dims:
        coord = da.coords.get(dim)

        # --- Pre-classify before calling _dim_to_axis_type ---
        # Check CF attrs and value-based signals first so that we never reach
        # _detect_axis_type (which warns for time-named dims under
        # STRICT_AXIS_DETECTION=True) when the axis type is unambiguous.
        ax_type = "-"
        if coord is not None:
            attrs = _coord_attrs_lower(coord)
            cf_axis = attrs.get("axis", "").upper()
            if cf_axis in {"T", "Y", "X", "Z"}:
                ax_type = cf_axis
            else:
                sn = attrs.get("standard_name", "")
                units_str = attrs.get("units", "")
                if (sn == "time"
                        or _is_datetime_like_time(coord)
                        or "since" in units_str):
                    ax_type = "T"
        # Fall through to _dim_to_axis_type for standard_name/units of
        # Y/X/Z and name heuristics that don't risk a spurious warning.
        if ax_type == "-":
            ax_type = _dim_to_axis_type(dim, coord)

        # Extra fallbacks when _dim_to_axis_type returns "-" (no CF metadata,
        # name heuristic also failed — e.g. STRICT_AXIS_DETECTION=True blocks
        # time-by-name detection in _detect_axis_type).
        if ax_type == "-":
            if coord is not None:
                attrs = _coord_attrs_lower(coord)
                if _is_datetime_like_time(coord):
                    ax_type = "T"
                elif "since" in attrs.get("units", ""):
                    ax_type = "T"

        if ax_type not in ("T", "Y", "X", "Z"):
            continue

        # Build the updated attrs dict — setdefault preserves existing values.
        attrs = dict(coord.attrs) if coord is not None else {}
        attrs["axis"] = ax_type
        if ax_type == "T":
            attrs.setdefault("standard_name", "time")
            attrs.setdefault("long_name", "time")
            if coord is not None:
                cal = (
                    coord.attrs.get("calendar")
                    or coord.encoding.get("calendar", None)
                    or "standard"
                )
                attrs.setdefault("calendar", cal)
        elif ax_type == "Y":
            attrs.setdefault("standard_name", "latitude")
            attrs.setdefault("units", "degrees_north")
        elif ax_type == "X":
            attrs.setdefault("standard_name", "longitude")
            attrs.setdefault("units", "degrees_east")
        elif ax_type == "Z":
            if coord is not None:
                attrs.setdefault("positive", coord.attrs.get("positive", "up"))

        if coord is not None:
            da = da.assign_coords(
                {dim: xr.DataArray(coord.values, dims=coord.dims, attrs=attrs)}
            )
        else:
            # Pure dimension (no coordinate) — create a minimal index so
            # da_to_cdat can attach the axis type to the _Axis object.
            # Use NaN-filled float for T (we do not invent time values) and
            # arange for spatial/vertical dims.
            if ax_type == "T":
                vals = np.full(da.sizes[dim], np.nan)
            else:
                vals = np.arange(da.sizes[dim], dtype=float)
            da = da.assign_coords(
                {dim: xr.DataArray(vals, dims=dim, attrs=attrs)}
            )

    if da.name is None:
        da.name = "var"
    return da


def _validate_cdat_axes(var: CDATVariable, *, context: str = "", require_time: bool = False) -> CDATVariable:
    """Fail early when data shape and CDAT-like axes diverge."""
    axes = var.getAxisList()
    if var.ndim != len(axes):
        raise ValueError(
            f"Axis/data mismatch after {context or 'conversion'}: "
            f"id={getattr(var, 'id', '')!r}, shape={var.shape}, axes={axes}"
        )
    if require_time and var.getTime() is None:
        raise ValueError(
            f"No time axis after {context or 'conversion'}: "
            f"id={getattr(var, 'id', '')!r}, shape={var.shape}, "
            f"axes={[getattr(ax, 'id', None) for ax in axes]}"
        )
    return var


def _finalize_cdat(
        da: xr.DataArray, varname: str | None = None, *,
        context: str = "", require_time: bool = False
    ) -> CDATVariable:
    """
    Convert a DataArray to CDATVariable after normalizing axis metadata.

    Use this at file-read boundaries and after xarray/xcdat/xesmf operations.
    It prevents later failures such as "no time axis", empty axis lists, and
    ambiguous dimension handling.
    """
    name = varname or da.name or "var"
    da = _standardize_da_axes(da.rename(name), context=context)
    out = da_to_cdat(da, varname=name)
    return _validate_cdat_axes(out, context=context, require_time=require_time)


def _finalize_existing_cdat(var, *, context: str = "", require_time: bool = False):
    """Validate an existing CDATVariable and refresh its internal grid.

    Must be called after any operation that can mutate or replace the
    CDATVariable without rebuilding the grid (squeeze, setAxis, MV2masked_where,
    arithmetic sign flip, slicing).  Rebuilding from the current axes ensures
    getGrid().shape always matches getLatitude()/getLongitude() axis lengths.
    """
    out = _to_cdat(var)
    # Rebuild the grid from whatever axes are currently attached so that
    # getGrid() is never stale after squeeze / setAxis / masking mutations.
    if isinstance(out, CDATVariable):
        out._grid = _build_grid_from_axes(out._axes)
    return _validate_cdat_axes(out, context=context, require_time=require_time)


def _curvilinear_lat_lon_coords(da: xr.DataArray):
    """Return 2-D latitude/longitude coordinate names for logically rectangular grids."""
    lat_name = None
    lon_name = None
    for cname, coord in da.coords.items():
        if coord.ndim != 2:
            continue
        ax_type = _dim_to_axis_type(cname, coord)
        if ax_type == "Y" and lat_name is None:
            lat_name = cname
        elif ax_type == "X" and lon_name is None:
            lon_name = cname
    if lat_name is None or lon_name is None:
        return None, None
    lat = da.coords[lat_name]
    lon = da.coords[lon_name]
    if lat.dims != lon.dims:
        return None, None
    if not all(dim in da.dims for dim in lat.dims):
        return None, None
    return lat_name, lon_name


def _longitude_mask(lon_values, lon_bounds):
    """Build a longitude mask after normalizing both data and bounds to 0..360."""
    lon = np.mod(np.asarray(lon_values, dtype=float), 360.0)
    lo, hi = [float(v) for v in lon_bounds]
    lo = lo % 360.0
    hi = hi % 360.0
    if lo <= hi:
        return (lon >= lo) & (lon <= hi)
    return (lon >= lo) | (lon <= hi)


def _subset_curvilinear_lat_lon(da: xr.DataArray, latitude=None, longitude=None):
    """
    Apply lat/lon region selection for 2-D curvilinear coordinates.

    xarray cannot use ``.sel(lat=slice(...), lon=slice(...))`` when latitude
    and longitude are auxiliary 2-D coordinates.  This masks by coordinate
    values and drops fully outside rows/columns, preserving a compact logical
    rectangle for downstream legacy operations.
    """
    lat_name, lon_name = _curvilinear_lat_lon_coords(da)
    if lat_name is None or lon_name is None:
        return da, False

    lat = da.coords[lat_name]
    lon = da.coords[lon_name]
    mask = xr.DataArray(
        np.ones(lat.shape, dtype=bool),
        dims=lat.dims,
        coords={dim: da.coords[dim] for dim in lat.dims if dim in da.coords},
    )
    if latitude is not None:
        lo, hi = min(latitude), max(latitude)
        mask = mask & ((lat >= lo) & (lat <= hi))
    if longitude is not None:
        mask = mask & xr.DataArray(
            _longitude_mask(lon.values, longitude),
            dims=lon.dims,
            coords={dim: da.coords[dim] for dim in lon.dims if dim in da.coords},
        )
    return da.where(mask, drop=True), True


def _target_axis_values(bounds, resolution=1.0, *, longitude=False):
    """Create regular target cell centers inside the requested bounds."""
    if bounds is None:
        return None
    lo, hi = [float(v) for v in bounds]
    if longitude:
        lo = lo % 360.0
        hi = hi % 360.0
        if hi <= lo:
            hi += 360.0
    else:
        lo, hi = sorted([lo, hi])
    npts = int(round((hi - lo) / resolution))
    if npts <= 0:
        npts = 1
    vals = lo + resolution / 2.0 + np.arange(npts) * resolution
    if longitude:
        vals = np.mod(vals, 360.0)
    return vals


def _regrid_curvilinear_to_rectilinear(
    da: xr.DataArray,
    *,
    latitude=None,
    longitude=None,
    resolution=1.0,
) -> xr.DataArray:
    """
    Regrid 2-D curvilinear lat/lon data to a regular lat-lon grid.

    This is used at the file-read boundary before converting to the CDAT-like
    compatibility object.  It handles logically rectangular ocean variables such
    as ``zos(time, j, i)`` with auxiliary ``lat(j, i)`` / ``lon(j, i)`` coords.
    """
    lat_name, lon_name = _curvilinear_lat_lon_coords(da)
    if lat_name is None or lon_name is None:
        return da

    lat = da.coords[lat_name]
    lon = da.coords[lon_name]
    ydim, xdim = lat.dims
    lat_vals_src = np.asarray(lat.values, dtype=float)
    lon_vals_src = np.mod(np.asarray(lon.values, dtype=float), 360.0)

    target_lat = _target_axis_values(
        latitude
        if latitude is not None
        else (np.nanmin(lat_vals_src), np.nanmax(lat_vals_src)),
        resolution,
    )
    target_lon = _target_axis_values(
        longitude
        if longitude is not None
        else (np.nanmin(lon_vals_src), np.nanmax(lon_vals_src)),
        resolution,
        longitude=True,
    )

    lat_attrs = dict(lat.attrs)
    lat_attrs["axis"] = "Y"
    lat_attrs.setdefault("standard_name", "latitude")
    lat_attrs.setdefault("units", "degrees_north")
    lon_attrs = dict(lon.attrs)
    lon_attrs["axis"] = "X"
    lon_attrs.setdefault("standard_name", "longitude")
    lon_attrs.setdefault("units", "degrees_east")

    target_ds = xr.Dataset(
        coords={
            "lat": xr.DataArray(target_lat, dims="lat", attrs=lat_attrs),
            "lon": xr.DataArray(target_lon, dims="lon", attrs=lon_attrs),
        }
    )

    if _HAS_XESMF:
        try:
            regridder = _xesmf.Regridder(
                da.to_dataset(name=da.name or "var"),
                target_ds,
                method="bilinear",
                extrap_method="nearest_s2d",
                reuse_weights=False,
            )
            result = regridder(da)
            result = result.assign_coords(
                {
                    "lat": target_ds["lat"],
                    "lon": target_ds["lon"],
                }
            )
            for dim in da.dims:
                if dim in (ydim, xdim, "lat", "lon"):
                    continue
                if dim in result.coords and dim in da.coords:
                    result = result.assign_coords({dim: da.coords[dim]})
            return result
        except Exception as e:
            warnings.warn(
                "xESMF curvilinear regridding failed; falling back to "
                f"scipy.interpolate.griddata. Original error: {type(e).__name__}: {e}",
                stacklevel=2,
            )

    target_lon_for_interp = target_lon.astype(float)
    source_lon_for_interp = lon_vals_src.copy()
    if longitude is not None:
        lo, hi = [float(v) % 360.0 for v in longitude]
        if hi <= lo:
            source_lon_for_interp = np.where(
                source_lon_for_interp < lo,
                source_lon_for_interp + 360.0,
                source_lon_for_interp,
            )
            target_lon_for_interp = np.where(
                target_lon_for_interp < lo,
                target_lon_for_interp + 360.0,
                target_lon_for_interp,
            )

    try:
        from scipy.interpolate import griddata as _scipy_griddata
    except ImportError as e:
        raise ImportError(
            "scipy is required to regrid curvilinear ocean variables "
            "to a regular lat-lon grid before ENSO metric calculation."
        ) from e

    src_points = np.column_stack(
        [source_lon_for_interp.ravel(), lat_vals_src.ravel()]
    )
    dst_lon2d, dst_lat2d = np.meshgrid(target_lon_for_interp, target_lat)
    dst_points = (dst_lon2d, dst_lat2d)

    data = np.asarray(da.values, dtype=float)
    spatial_shape = lat_vals_src.shape
    if data.shape[-2:] != spatial_shape:
        axis_order = [dim for dim in da.dims if dim not in (ydim, xdim)] + [ydim, xdim]
        da = da.transpose(*axis_order)
        data = np.asarray(da.values, dtype=float)

    leading_shape = data.shape[:-2]
    flat_data = data.reshape((-1,) + spatial_shape)
    out = np.full((flat_data.shape[0], len(target_lat), len(target_lon)), np.nan, dtype=float)
    finite_points = np.isfinite(src_points).all(axis=1)

    for idx, field in enumerate(flat_data):
        values = field.ravel()
        valid = finite_points & np.isfinite(values)
        if np.count_nonzero(valid) < 3:
            continue
        out_field = _scipy_griddata(
            src_points[valid],
            values[valid],
            dst_points,
            method="linear",
        )
        missing = ~np.isfinite(out_field)
        if np.any(missing):
            nearest = _scipy_griddata(
                src_points[valid],
                values[valid],
                dst_points,
                method="nearest",
            )
            out_field = np.where(missing, nearest, out_field)
        out[idx] = out_field

    out = out.reshape(leading_shape + (len(target_lat), len(target_lon)))

    drop_names = [name for name in (lat_name, lon_name) if name not in da.dims]
    if drop_names:
        da = da.drop_vars(drop_names)

    dims = tuple([dim for dim in da.dims if dim not in (ydim, xdim)] + ["lat", "lon"])
    coords = {}
    for dim in dims:
        if dim == "lat":
            coords[dim] = xr.DataArray(target_lat, dims=dim, attrs=lat_attrs)
        elif dim == "lon":
            coords[dim] = xr.DataArray(target_lon, dims=dim, attrs=lon_attrs)
        elif dim in da.coords:
            coords[dim] = da.coords[dim]

    return xr.DataArray(
        out,
        dims=dims,
        coords=coords,
        attrs=dict(da.attrs),
        name=da.name,
    )


def _regrid_unstructured_to_rectilinear(
    da: xr.DataArray,
    ds: xr.Dataset | None = None,
    *,
    latitude=None,
    longitude=None,
    resolution=1.0,
) -> xr.DataArray:
    """
    Regrid MPAS-like unstructured cell data to a regular lat-lon grid.

    This handles variables such as ``tos(time, nCells)`` or ``zos(time, nCells)``
    with 1-D ``latCell(nCells)`` / ``lonCell(nCells)`` coordinates.
    """
    lat_name, lon_name, cell_dim, lat, lon = _unstructured_lat_lon_coords(da, ds)
    if cell_dim is None:
        raise NotImplementedError(
            "Unstructured/MPAS-like data were detected, but no 1-D cell "
            "latitude/longitude coordinates were found. Expected coordinates "
            "such as latCell/lonCell on the same nCells dimension. Remap to "
            "regular lat-lon before running ENSO_metrics, or provide CF-style "
            "cell latitude/longitude coordinates."
        )

    lat_vals_src = _coord_values_degrees(lat)
    lon_vals_src = np.mod(_coord_values_degrees(lon), 360.0)

    target_lat = _target_axis_values(
        latitude
        if latitude is not None
        else (np.nanmin(lat_vals_src), np.nanmax(lat_vals_src)),
        resolution,
    )
    target_lon = _target_axis_values(
        longitude
        if longitude is not None
        else (np.nanmin(lon_vals_src), np.nanmax(lon_vals_src)),
        resolution,
        longitude=True,
    )

    target_lon_for_interp = target_lon.astype(float)
    source_lon_for_interp = lon_vals_src.copy()
    if longitude is not None:
        lo, hi = [float(v) % 360.0 for v in longitude]
        if hi <= lo:
            source_lon_for_interp = np.where(
                source_lon_for_interp < lo,
                source_lon_for_interp + 360.0,
                source_lon_for_interp,
            )
            target_lon_for_interp = np.where(
                target_lon_for_interp < lo,
                target_lon_for_interp + 360.0,
                target_lon_for_interp,
            )

    try:
        from scipy.interpolate import griddata as _scipy_griddata
    except ImportError as e:
        raise ImportError(
            "scipy is required to regrid MPAS-like unstructured ocean variables "
            "to a regular lat-lon grid before ENSO metric calculation."
        ) from e

    src_points = np.column_stack([source_lon_for_interp, lat_vals_src])
    dst_lon2d, dst_lat2d = np.meshgrid(target_lon_for_interp, target_lat)
    dst_points = (dst_lon2d, dst_lat2d)

    if da.dims[-1] != cell_dim:
        axis_order = [dim for dim in da.dims if dim != cell_dim] + [cell_dim]
        da = da.transpose(*axis_order)
    data = np.asarray(da.values, dtype=float)
    leading_shape = data.shape[:-1]
    flat_data = data.reshape((-1, data.shape[-1]))
    out = np.full((flat_data.shape[0], len(target_lat), len(target_lon)), np.nan, dtype=float)
    finite_points = np.isfinite(src_points).all(axis=1)

    for idx, field in enumerate(flat_data):
        valid = finite_points & np.isfinite(field)
        if np.count_nonzero(valid) < 3:
            continue
        out_field = _scipy_griddata(
            src_points[valid],
            field[valid],
            dst_points,
            method="linear",
        )
        missing = ~np.isfinite(out_field)
        if np.any(missing):
            nearest = _scipy_griddata(
                src_points[valid],
                field[valid],
                dst_points,
                method="nearest",
            )
            out_field = np.where(missing, nearest, out_field)
        out[idx] = out_field

    out = out.reshape(leading_shape + (len(target_lat), len(target_lon)))
    lat_attrs = {"axis": "Y", "standard_name": "latitude", "units": "degrees_north"}
    lon_attrs = {"axis": "X", "standard_name": "longitude", "units": "degrees_east"}
    dims = tuple([dim for dim in da.dims if dim != cell_dim] + ["lat", "lon"])
    coords = {}
    for dim in dims:
        if dim == "lat":
            coords[dim] = xr.DataArray(target_lat, dims=dim, attrs=lat_attrs)
        elif dim == "lon":
            coords[dim] = xr.DataArray(target_lon, dims=dim, attrs=lon_attrs)
        elif dim in da.coords:
            coords[dim] = da.coords[dim]

    return xr.DataArray(
        out,
        dims=dims,
        coords=coords,
        attrs=dict(da.attrs),
        name=da.name,
    )


# ---------------------------------------------------------------------------
# Revised file handle
# ---------------------------------------------------------------------------
class _XcDatasetHandle:
    """
    Compatibility wrapper that mimics the subset of the legacy ``cdms2.open()``
    file-handle interface used by ENSO_metrics.

    Internally, this class uses xarray/xcdat-compatible reading and writing,
    while preserving the legacy call pattern ``handle(varname, ...)`` and basic
    write-mode behavior expected by downstream ENSO_metrics routines.
    """

    def __init__(self, path: str, mode: str = "r"):
        self._path = path
        self._mode = mode
        self._ds = None
        self._open_error = None
        self._write_vars = {}
        self._global_attrs = {}

        if mode in ("r", "", "a"):
            try:
                ds_raw = xr.open_dataset(path, decode_times=False)
                ds_raw = _add_cf_units_to_ds(ds_raw)
                ds_raw = _fix_leap_seconds_in_raw(ds_raw)

                try:
                    coder = xr.coders.CFDatetimeCoder(use_cftime=True)
                    try:
                        ds = xr.decode_cf(ds_raw, decode_times=coder)
                    except TypeError:
                        ds = xr.decode_cf(ds_raw, use_cftime=True)
                except AttributeError:
                    ds = xr.decode_cf(ds_raw, use_cftime=True)

                self._ds = _ensure_time_encoding(ds, path=path)

            except Exception as e:
                self._open_error = e
                self._ds = None

    def __call__(self, varname: str, **kwargs):
        """Read variable, optionally subset by time/latitude/longitude."""
        if self._ds is None:
            raise IOError(
                f"File not open for reading: {self._path}; "
                f"original error: {self._open_error}"
            )

        da = self._ds[varname]
        curvilinear_latlon = False
        unstructured_latlon = False

        if _is_unstructured_grid(da):
            da, unstructured_latlon = _subset_unstructured_lat_lon(
                da,
                self._ds,
                latitude=kwargs.get("latitude"),
                longitude=kwargs.get("longitude"),
            )

        if "latitude" in kwargs and not unstructured_latlon:
            lat_bnds = kwargs["latitude"]
            da, curvilinear_latlon = _subset_curvilinear_lat_lon(
                da, latitude=lat_bnds, longitude=kwargs.get("longitude")
            )
            if not curvilinear_latlon:
                lo, hi = min(lat_bnds), max(lat_bnds)
                lat_dim = _validate_axis(da, "Y", context="latitude selection")
                if lat_dim:
                    da = da.sel({lat_dim: slice(lo, hi)})

        if "longitude" in kwargs and "latitude" not in kwargs and not unstructured_latlon:
            lon_bnds = kwargs["longitude"]
            da, curvilinear_latlon = _subset_curvilinear_lat_lon(
                da, longitude=lon_bnds
            )

        if "longitude" in kwargs and not curvilinear_latlon and not unstructured_latlon:
            lon_bnds = kwargs["longitude"]
            lo, hi = min(lon_bnds), max(lon_bnds)
            lon_dim = _validate_axis(da, "X", context="longitude selection")
            if lon_dim:
                da = da.sel({lon_dim: slice(lo, hi)})

        if "time" in kwargs:
            t_bnds = kwargs["time"]
            t_dim = _require_axis(da, "T", context="time selection")
            if t_dim:
                if not _is_datetime_like_time(da[t_dim] if t_dim in da.coords else self._ds[t_dim]):
                    raise RuntimeError(
                        f"time axis not decoded before selection in {self._path}"
                    )
                da = da.sel(
                    {
                        t_dim: slice(
                            _to_sel_bound(t_bnds[0], t_dim, da),
                            _to_sel_bound(t_bnds[1], t_dim, da),
                        )
                    }
                )

        if kwargs.get("squeeze"):
            da = da.squeeze()

        if unstructured_latlon or _is_unstructured_grid(da):
            da = _regrid_unstructured_to_rectilinear(
                da,
                self._ds,
                latitude=kwargs.get("latitude"),
                longitude=kwargs.get("longitude"),
            )
        elif curvilinear_latlon or _curvilinear_lat_lon_coords(da) != (None, None):
            da = _regrid_curvilinear_to_rectilinear(
                da,
                latitude=kwargs.get("latitude"),
                longitude=kwargs.get("longitude"),
            )

        # Build a robust CDAT-like object at the read boundary.  This keeps
        # T/Y/X/Z metadata attached before downstream PCMDI code sees it.
        require_time = any(_coord_is_datetime_like(da, dim) for dim in da.dims)
        return _finalize_cdat(
            da, varname=varname,
            context=f"read:{self._path}:{varname}",
            require_time=require_time
        )

    @staticmethod
    def _clean_dataset_for_write(ds):
        """
        Clean dataset attributes and encodings before NetCDF writing.

        This reduces failures caused by stale xarray/backend encodings or
        non-NetCDF-safe attributes when emulating legacy CDAT append behavior.
        """
        ds = ds.copy()

        ds.attrs = _clean_attrs(getattr(ds, "attrs", {}))

        for name in list(ds.variables):
            ds[name].attrs = _clean_attrs(ds[name].attrs)

            # Drop stale backend-specific encodings that often cause write
            # conflicts when rewriting merged NetCDF files.
            keep_encoding = {}
            for key in ("_FillValue", "dtype", "zlib", "complevel", "chunksizes"):
                if key in ds[name].encoding:
                    keep_encoding[key] = ds[name].encoding[key]
            ds[name].encoding = keep_encoding

        return ds

    @staticmethod
    def _coords_compatible(coord_old, coord_new):
        """
        Return True if two coordinates are safely compatible for merge.

        Numeric coordinates are compared with allclose; non-numeric coordinates
        are compared with exact equality.
        """
        try:
            if coord_old.shape != coord_new.shape:
                return False

            old_vals = np.asarray(coord_old.values)
            new_vals = np.asarray(coord_new.values)

            try:
                return np.allclose(
                    old_vals.astype(float),
                    new_vals.astype(float),
                    equal_nan=True,
                    atol=1e-8,
                    rtol=1e-8,
                )
            except Exception:
                return np.array_equal(old_vals, new_vals)

        except Exception:
            return False

    @staticmethod
    def _merge_datasets_like_append(ds_old, ds_new):
        """
        Merge ``ds_new`` into ``ds_old`` while preserving old variables.

        This emulates the legacy CDAT append/write behavior used by
        ``CDMS2open(path, "a")``: new variables replace same-named old
        variables, while unrelated existing variables are preserved.
        """
        ds_old = _XcDatasetHandle._clean_dataset_for_write(ds_old)
        ds_new = _XcDatasetHandle._clean_dataset_for_write(ds_new)

        replace_names = [
            v for v in ds_new.data_vars
            if v in ds_old.data_vars
        ]

        ds_old_keep = ds_old.drop_vars(
            replace_names,
            errors="ignore",
        )

        ds_new_safe = ds_new.copy()
        for vname in list(ds_new_safe.data_vars):
            da = ds_new_safe[vname]
            rename_dims = {}
            safe_vname = (
                str(vname)
                .replace("/", "_")
                .replace(" ", "_")
                .replace(":", "_")
            )

            for cname in da.dims:
                if cname not in da.coords or cname not in ds_old_keep.coords:
                    continue
                if _XcDatasetHandle._coords_compatible(ds_old_keep[cname], da[cname]):
                    continue

                new_cname = f"{cname}_{safe_vname}"
                if new_cname in ds_old_keep.coords or new_cname in ds_old_keep.dims:
                    suffix = 1
                    base_name = new_cname
                    while (
                        new_cname in ds_old_keep.coords
                        or new_cname in ds_old_keep.dims
                        or new_cname in ds_new_safe.coords
                        or new_cname in ds_new_safe.dims
                    ):
                        suffix += 1
                        new_cname = f"{base_name}_{suffix}"

                rename_dims[cname] = new_cname

            if rename_dims:
                ds_new_safe = ds_new_safe.drop_vars(vname)
                ds_new_safe[vname] = da.rename(rename_dims)

        used_dims = {
            dim
            for da in ds_new_safe.data_vars.values()
            for dim in da.dims
        }
        unused_dim_coords = [
            cname for cname in ds_new_safe.coords
            if cname in ds_new_safe.dims and cname not in used_dims
        ]
        if unused_dim_coords:
            ds_new_safe = ds_new_safe.drop_vars(unused_dim_coords, errors="ignore")

        ds_new = ds_new_safe

        try:
            ds_merged = xr.merge(
                [ds_old_keep, ds_new],
                compat="override",
                join="outer",
            )

        except Exception:
            # Fallback: append variables one by one so one coordinate conflict
            # does not cause all previously written variables to be lost.
            ds_merged = ds_old_keep.copy()

            for vname in ds_new.data_vars:
                da = ds_new[vname]

                for cname in list(da.coords):
                    if cname not in ds_merged.coords:
                        ds_merged = ds_merged.assign_coords({cname: da[cname]})
                        continue

                    if _XcDatasetHandle._coords_compatible(
                        ds_merged[cname],
                        da[cname],
                    ):
                        continue

                    # If a dimension coordinate conflicts, rename it to a
                    # variable-specific coordinate. This preserves the new
                    # variable without corrupting old variables that already
                    # depend on the existing coordinate.
                    if cname in da.dims:
                        safe_vname = (
                            str(vname)
                            .replace("/", "_")
                            .replace(" ", "_")
                            .replace(":", "_")
                        )
                        new_cname = f"{cname}_{safe_vname}"

                        # Avoid accidental collision if the generated name
                        # already exists in the merged dataset.
                        if new_cname in ds_merged.coords or new_cname in ds_merged.dims:
                            suffix = 1
                            base_name = new_cname
                            while (
                                new_cname in ds_merged.coords
                                or new_cname in ds_merged.dims
                            ):
                                suffix += 1
                                new_cname = f"{base_name}_{suffix}"

                        da = da.rename({cname: new_cname})

                    else:
                        # Non-dimension coordinate conflict. Drop this coordinate
                        # from the new variable rather than corrupting an existing
                        # coordinate used by old variables.
                        da = da.drop_vars(cname, errors="ignore")

                ds_merged[vname] = da

        merged_attrs = dict(getattr(ds_old, "attrs", {}))
        merged_attrs.update(_clean_attrs(getattr(ds_new, "attrs", {})))
        ds_merged.attrs = _clean_attrs(merged_attrs)

        return _XcDatasetHandle._clean_dataset_for_write(ds_merged)

    def write(self, var, attributes=None, dtype="float32", id=None):
        """Buffer a variable for writing."""
        name = id or (var.id if isinstance(var, CDATVariable) else "var")
        if isinstance(var, CDATVariable):
            da = cdat_to_da(var, name=name)
        elif isinstance(var, xr.DataArray):
            da = var.rename(name)
        else:
            # Scalar float/int or plain numpy array — wrap in a DataArray.
            da = xr.DataArray(np.asarray(var), name=name)

        if attributes:
            da.attrs.update(_clean_attrs(attributes))

        da = da.rename(name).astype(dtype)
        rename_dims = {}
        safe_name = str(name).replace("/", "_").replace(" ", "_").replace(":", "_")
        for existing in self._write_vars.values():
            for cname in da.dims:
                if cname in rename_dims:
                    continue
                if cname not in da.coords or cname not in existing.coords:
                    continue
                if self._coords_compatible(existing[cname], da[cname]):
                    continue
                new_cname = f"{cname}_{safe_name}"
                suffix = 1
                base_name = new_cname
                used = set(da.coords) | set(da.dims)
                for old_da in self._write_vars.values():
                    used |= set(old_da.coords) | set(old_da.dims)
                while new_cname in used:
                    suffix += 1
                    new_cname = f"{base_name}_{suffix}"
                rename_dims[cname] = new_cname

        if rename_dims:
            da = da.rename(rename_dims)

        self._write_vars[name] = da

    def __setattr__(self, key, value):
        if key.startswith("_") or key in (
            "_path",
            "_mode",
            "_ds",
            "_open_error",
            "_write_vars",
            "_global_attrs",
        ):
            super().__setattr__(key, value)
        else:
            try:
                self._global_attrs[key] = value
            except AttributeError:
                super().__setattr__(key, value)

    def close(self):
        if self._mode in ("w", "w+", "a") and self._write_vars:
            _os_sn.makedirs(_os_sn.path.dirname(self._path) or ".", exist_ok=True)

            # If this handle was opened in append mode, __init__ may have opened
            # the existing file for reading. Close it before reopening/replacing.
            if self._ds is not None:
                self._ds.close()
                self._ds = None

            ds_new = xr.Dataset(
                self._write_vars,
                attrs=_clean_attrs(self._global_attrs),
            ).load()

            ds_new = self._clean_dataset_for_write(ds_new)

            ds_merged = None

            if _os_sn.path.exists(self._path):
                try:
                    with xr.open_dataset(
                        self._path,
                        engine="netcdf4",
                        chunks={},
                        decode_times=False,
                    ) as ds_old_open:
                        ds_old = ds_old_open.load()

                    ds_merged = self._merge_datasets_like_append(ds_old, ds_new)

                except Exception as e:
                    try:
                        ds_new.close()
                    except Exception:
                        pass

                    raise RuntimeError(
                        f"Failed to append variables to existing NetCDF file: {self._path!r}. "
                        "The existing file was not overwritten. "
                        f"New variables were: {list(self._write_vars.keys())}. "
                        f"Original error: {type(e).__name__}: {e}"
                    ) from e

            else:
                ds_merged = ds_new

            _fd, tmpfile = _tmp_sn.mkstemp(
                prefix=ntpath.basename(self._path) + ".",
                suffix=".tmp",
                dir=ntpath.dirname(self._path) or ".",
            )
            _os_sn.close(_fd)

            try:
                ds_merged.to_netcdf(tmpfile, mode="w", format="NETCDF4")
                _os_sn.replace(tmpfile, self._path)

            except Exception as e:
                if _os_sn.path.exists(tmpfile):
                    _os_sn.remove(tmpfile)
                raise RuntimeError(
                    f"Failed to write merged NetCDF file atomically: {self._path!r}. "
                    f"Original error: {type(e).__name__}: {e}"
                ) from e

            finally:
                if ds_merged is not None:
                    ds_merged.close()
                if ds_new is not ds_merged:
                    ds_new.close()

        if self._ds is not None:
            self._ds.close()
            self._ds = None

# ---------------------------------------------------------------------------
# Module-level lookup tables for _guess_dim — defined once, not per-call.
# ---------------------------------------------------------------------------
_AXIS_STANDARD_NAMES: dict[str, set[str]] = {
    "Y": {"latitude", "grid_latitude", "projection_y_coordinate","rotated_latitude"},
    "X": {"longitude", "grid_longitude", "projection_x_coordinate","rotated_longitude"},
    "T": {"time"},
    "Z": {
        "air_pressure", "altitude", "depth", "height",
        "ocean_sigma_coordinate", "sigma", "eta",
        "height_above_geopotential_datum",
        "height_above_mean_sea_level"
    },
}
_AXIS_UNITS_PATTERNS: dict[str, tuple[str, ...]] = {
    "Y": ("degrees_north", "degree_north", "degrees_n", "degree_n", "degreesnorth", "degreen"),
    "X": ("degrees_east", "degree_east", "degrees_e", "degree_e", "degreeseast", "degreee"),
    "T": ("since",),   # matches "days since …", "hours since …", etc.
    "Z": ("pa", "hpa", "mb", "mbar", "meter", "m", "sigma", "hybrid"),
}
_AXIS_NAME_HINTS: dict[str, tuple[str, ...]] = {
    "Y": ("lat", "latitude", "nav_lat", "rlat", "y", "j", "nlat", "y_1", "y_2"),
    "X": ("lon", "longitude", "nav_lon", "rlon", "x", "i", "nlon", "x_1", "x_2"),
    "T": ("time", "t"),
    "Z": ("lev", "level", "plev", "depth", "sigma", "eta", "z", "k", "nlev"),
}

def _safe_guess_dim(da: xr.DataArray, axis_type: str) -> str:
    """Best-effort axis lookup used only for metadata finalization."""
    return _detect_axis(da, axis_type)


def _guess_dim(
    da: xr.DataArray,
    axis_type: str,
    *,
    strict: bool | None = None,
    return_score: bool = False,
) -> str | tuple[str, int]:
    """Return the dimension name for a given axis type (T/Y/X/Z).

    Scoring lookup:
        score 4 — CF axis attribute
        score 3 — CF standard_name attribute
        score 2 — CF units pattern
        score 1 — token/endswith name heuristic
        score 0 — coordinate value-range
    """
    at = axis_type.upper()
    sn_set = _AXIS_STANDARD_NAMES.get(at, set())
    u_pats = _AXIS_UNITS_PATTERNS.get(at, ())
    hints = _AXIS_NAME_HINTS.get(at, ())

    scores: list[tuple[int, str]] = []

    for dim in da.dims:
        score = -1

        if dim in da.coords:
            coord = da.coords[dim]
            attrs = _coord_attrs_lower(coord)

            if attrs.get("axis", "").upper() == at:
                score = max(score, 4)

            if attrs.get("standard_name", "") in sn_set:
                score = max(score, 3)

            u = attrs.get("units", "")
            if any(p in u for p in u_pats):
                score = max(score, 2)

            if score < 0 and at in ("Y", "X"):
                try:
                    raw = np.asarray(coord.values)

                    if not (
                        np.issubdtype(raw.dtype, np.integer)
                        and raw.size >= 1
                        and int(raw.min()) == 0
                        and int(raw.max()) == raw.size - 1
                    ):
                        vals = raw.astype(float)
                        vals = vals[np.isfinite(vals)]

                        if vals.size == 0:
                            continue

                        vmin = float(vals.min())
                        vmax = float(vals.max())

                        if at == "Y" and -90.0 <= vmin <= vmax <= 90.0:
                            score = max(score, 0)

                        if at == "X" and (
                            (-180.0 <= vmin <= vmax <= 180.0)
                            or (0.0 <= vmin <= vmax <= 360.0)
                        ):
                            score = max(score, 0)

                except Exception:
                    pass

        if score < 0:
            dl = dim.lower()
            tokens = set(dl.replace("-", "_").split("_"))

            if (
                dl in hints
                or any(h in tokens for h in hints)
                or any(dl.endswith(h) for h in hints)
            ):
                score = max(score, 1)

        if score >= 0:
            scores.append((score, dim))

    if not scores:
        _strict = STRICT_DIM_GUESS if strict is None else strict
        msg = (
            f"_guess_dim cannot determine axis {axis_type!r} from dims "
            f"{list(da.dims)}: no candidate matched. "
            "Add CF axis/standard_name/units metadata (e.g., 'axis', 'standard_name', or 'units' attributes) to your coordinates. "
            "If working with legacy data, set STRICT_DIM_GUESS = False to allow fallback heuristics."
        )

        if _strict:
            raise ValueError(msg)

        warnings.warn(msg, stacklevel=2)
        return ("", -1) if return_score else ""

    dim_order = {d: i for i, d in enumerate(da.dims)}
    best_score, best_dim = max(
        scores,
        key=lambda t: (t[0], -dim_order.get(t[1], 0)),
    )

    if best_score <= 0:
        _strict = STRICT_DIM_GUESS if strict is None else strict
        msg = (
            f"_guess_dim low-confidence value-range fallback for axis "
            f"{axis_type!r}: chosen={best_dim!r}, score={best_score}. "
            "Only coordinate-value heuristics matched; results may be "
            "scientifically unreliable. Add CF metadata."
        )

        if _strict:
            raise ValueError(msg)

        warnings.warn(msg, stacklevel=2)
        return (best_dim, best_score) if return_score else best_dim

    if best_score == 1:
        msg = (
            f"_guess_dim name-heuristic fallback: axis={axis_type!r}, "
            f"chosen={best_dim!r}, score={best_score}, dims={list(da.dims)}. "
            "Add CF axis/standard_name/units metadata to avoid incorrect axis mapping."
        )

        if at == "T":
            warnings.warn(msg, stacklevel=2)
            return ("", best_score) if return_score else ""

        _strict = STRICT_DIM_GUESS if strict is None else strict
        if _strict:
            raise ValueError(msg)

        warnings.warn(msg, stacklevel=2)
        return (best_dim, best_score) if return_score else best_dim

    return (best_dim, best_score) if return_score else best_dim

# cdutil averager stub (used inside this module)
class _CdutilAverager:
    """Wraps xcdat spatial/temporal averaging to mimic cdutil.averager()."""
    @staticmethod
    def averager(tab, axis="xy", weights="weighted", action="average"):
        if weights is None:
            warnings.warn(
                "cdutil.averager called with weights=None (equal weights); "
                "CDAT default is cosine-latitude weighting (weights='weighted'). "
                "Results may differ from CDAT.",
                stacklevel=2,
            )
        da = cdat_to_da(tab)
        varname = getattr(tab, 'id', None) or "var"
        # Reject unstructured grids early — they need a different code path
        if _is_unstructured_grid(da):
            raise NotImplementedError(
                "Native unstructured/MPAS grids are not supported by this CDAT replacement path. "
                "Please remap to regular lat-lon first."
            )
        ds = da.to_dataset(name=varname)
        # add_missing_bounds: only include "T" when time is actually decoded so
        # xcdat does not emit "Bounds cannot be created for 'time'" warnings.
        _amb_axes = ["X", "Y"]
        if "time" in ds.coords and _is_datetime_like_time(ds["time"]):
            _amb_axes.append("T")
        try:
            ds = ds.bounds.add_missing_bounds(axes=_amb_axes)
        except Exception:
            try:
                ds = ds.bounds.add_missing_bounds(axes=["X", "Y"])
            except Exception:
                pass  # proceed without bounds; spatial.average will still work
        axis_s = axis.lower().replace(" ", "")
        xcdat_axes = []
        do_time = False
        if "x" in axis_s or "lon" in axis_s:
            xcdat_axes.append("X")
        if "y" in axis_s or "lat" in axis_s:
            xcdat_axes.append("Y")
        if "t" in axis_s or "time" in axis_s:
            do_time = True
        # Numeric axis string e.g. "01", "10", "12" — map via CDATVariable metadata
        if not xcdat_axes and not do_time and axis_s.isdigit():
            ndim = len(tab.shape)
            for c in axis_s:
                idx = int(c)
                ax_type = None
                if isinstance(tab, CDATVariable) and idx < len(tab._axes):
                    ax = tab._axes[idx]
                    ax_type = ax.axis if ax is not None else None
                if ax_type == "T":
                    do_time = True
                elif ax_type == "Y":
                    xcdat_axes.append("Y")
                elif ax_type == "X":
                    xcdat_axes.append("X")
                else:
                    # NOTE: fallback assumes (t, y, x) layout; may not hold
                    # for staggered grids or E3SM/MPAS unconventional ordering.
                    if idx == ndim - 2:
                        xcdat_axes.append("Y")
                    elif idx == ndim - 1:
                        xcdat_axes.append("X")
                    elif idx == 0 and ndim >= 3:
                        do_time = True
        if do_time and not xcdat_axes:
            t_dim = _require_axis(da, "T", context="time selection") or "time"
            result = _finalize_cdat(
                ds[varname].mean(dim=t_dim), varname=varname,
                context="cdutil.averager:time"
            )
            return _apply_nan_majority_policy(result, tab, "t")
        if xcdat_axes:
            # Cosine-latitude weighted average — deterministic: always use manual
            # implementation when weights="weighted" and Y is in the reduction
            # axes.  Only fall through to xcdat when lat weighting is explicitly
            # not requested or latitude axis is unavailable.
            if weights == "weighted" and "Y" in xcdat_axes:
                try:
                    result_raw = _weighted_spatial_average(tab, axes=tuple(xcdat_axes))
                    if do_time:
                        spatial_result = ma.array(result_raw, copy=True)
                        result_raw = ma.mean(result_raw, axis=0)
                        result_raw = ma.array(
                            result_raw,
                            mask=(
                                ma.getmaskarray(result_raw)
                                | _nan_majority_mask(spatial_result, 0)
                            ),
                        )
                    if isinstance(tab, CDATVariable) and isinstance(result_raw, (np.ndarray, ma.MaskedArray)):
                        reduce_types = set(xcdat_axes) | ({"T"} if do_time else set())
                        surviving = [
                            ax.copy() if ax is not None else None
                            for ax in tab._axes
                            if ax is None or ax.axis not in reduce_types
                        ]
                        return _finalize_existing_cdat(
                            CDATVariable(
                                result_raw,
                                axes=surviving,
                                grid=None,
                                id=varname,
                                attributes=dict(tab._attributes),
                            ),
                            context="cdutil.averager:weighted",
                        )
                    return _finalize_existing_cdat(
                        CDATVariable(result_raw, id=varname),
                        context="cdutil.averager:weighted",
                    )
                except Exception as _e:
                    raise RuntimeError(
                        "_weighted_spatial_average failed — cannot guarantee "
                        "correctness: " + str(_e)
                    ) from _e
            result_ds = ds.spatial.average(varname, axis=xcdat_axes)
            result = result_ds[varname]
            if do_time:
                t_dim = _require_axis(da, "T", context="time selection") or "time"
                if t_dim in result.dims:
                    result = result.mean(dim=t_dim)
            result_cdat = _finalize_cdat(
                result, varname=varname,
                context="cdutil.averager:spatial",
                require_time=("T" not in set(xcdat_axes) and _has_time_axis(tab))
            )
            reduce_axis = "".join([ax.lower() for ax in xcdat_axes])
            if do_time:
                reduce_axis += "t"
            return _apply_nan_majority_policy(result_cdat, tab, reduce_axis)
        return tab.copy()

    @staticmethod
    def setTimeBoundsMonthly(tab):
        # Intentional no-op: xcdat adds its own time bounds internally via
        # ds.bounds.add_missing_bounds(axes=["T"]) before every temporal
        # average.  Storing explicit bounds on the _Axis object is not
        # required by the xcdat-based pipeline.
        pass
    @staticmethod
    def setTimeBoundsDaily(tab):
        pass  # see setTimeBoundsMonthly
    @staticmethod
    def setTimeBoundsYearly(tab):
        pass  # see setTimeBoundsMonthly

    class ANNUALCYCLE:
        @staticmethod
        def departures(tab):
            da = cdat_to_da(tab)
            varname = tab.id or "var"
            ds = da.to_dataset(name=varname)
            ds = _ensure_time_encoding(ds, path="in-memory:ANNUALCYCLE")
            try:
                if "time" in ds.coords and _is_datetime_like_time(ds["time"]):
                    ds = ds.bounds.add_missing_bounds(axes=["T"])
            except (KeyError, Exception):
                pass  # proceed without time bounds
            result = ds.temporal.departures(varname, freq="month", weighted=True)
            return _finalize_cdat(
                result[varname], varname=varname,
                context="ANNUALCYCLE.departures",
                require_time=True
                )


    @staticmethod
    def generateLandSeaMask(d, debug=False):
        """
        Generate an estimated land-sea mask for the ENSO_metrics compatibility layer.

        This method preserves the subset of the legacy ``cdutil.generateLandSeaMask``
        behavior needed by ENSO_metrics while using ``regionmask`` and Natural Earth
        land polygons internally.

        Parameters
        ----------
        d : CDATVariable or xarray.DataArray
            Input field whose latitude/longitude coordinates define the target grid.
            Singleton time dimensions are dropped before mask generation.

        debug : boolean, optional
            If True, print diagnostic information during mask generation.
            default value = False

        Returns
        -------
        CDATVariable
            Estimated land-fraction mask on the input horizontal grid, with
            ``1.0`` over land and ``0.0`` over ocean. The returned mask preserves
            the same latitude/longitude coordinates and longitude convention as the
            input field. Callers may multiply by ``100`` if downstream code expects
            ``sftlf`` percent units.
        """
        if not _HAS_REGIONMASK:
            raise RuntimeError(
                "Land-sea mask generation failed because regionmask is not available. "
                "Install regionmask with: conda install -c conda-forge regionmask. "
                "Alternatively, provide an explicit sftlf file."
            )

        da = cdat_to_da(d) if isinstance(d, CDATVariable) else d

        # Drop non-spatial singleton dimensions, e.g. time=1.
        for dim in list(da.dims):
            if dim not in da.coords:
                continue
            ax = str(da[dim].attrs.get("axis", "")).upper()
            if ax == "T" and da.sizes.get(dim, 1) == 1:
                da = da.isel({dim: 0}, drop=True)

        try:
            lat_name = _validate_axis(da, "Y", context="generateLandSeaMask")
            lon_name = _validate_axis(da, "X", context="generateLandSeaMask")

            for nm, atype in [(lat_name, "Y"), (lon_name, "X")]:
                coord = da[nm]
                attrs = _coord_attrs_lower(coord)

                has_cf = (
                    attrs.get("axis", "").upper() == atype
                    or attrs.get("standard_name", "") in _AXIS_STANDARD_NAMES[atype]
                )

                has_units = any(
                    p in attrs.get("units", "")
                    for p in _AXIS_UNITS_PATTERNS[atype]
                )

                if not (has_cf or has_units):
                    warnings.warn(
                        f"Low-confidence {atype} axis detection for {nm!r}; "
                        "mask generation may be unreliable.",
                        stacklevel=2,
                    )

            lat = da[lat_name]
            lon = da[lon_name]

            lat_vals = lat.values.astype(float)
            lon_vals = lon.values.astype(float)

            if np.all(~np.isfinite(lat_vals)):
                raise ValueError("Latitude coordinate contains no finite values")

            if np.all(~np.isfinite(lon_vals)):
                raise ValueError("Longitude coordinate contains no finite values")

            lat_min = float(np.nanmin(lat_vals))
            lat_max = float(np.nanmax(lat_vals))

            if lat_min < -90.0 or lat_max > 90.0:
                raise ValueError(
                    "Latitude coordinate is outside the valid geographic range "
                    f"[-90, 90]: min={lat_min}, max={lat_max}. "
                    "This likely indicates incorrect latitude-axis detection or "
                    "non-geographic coordinates."
                )

            if debug:
                print("[DEBUG] generateLandSeaMask: lat shape:", lat.shape, "lon shape:", lon.shape)
                print("[DEBUG] generateLandSeaMask: lat min/max:", np.nanmin(lat_vals), np.nanmax(lat_vals))
                print("[DEBUG] generateLandSeaMask: lon min/max:", np.nanmin(lon_vals), np.nanmax(lon_vals))

            # Convert longitude only for regionmask polygon lookup.
            # The returned mask is restored to the original input longitude convention.
            lon_vals_for_mask = lon_vals.copy()

            if np.nanmax(lon_vals_for_mask) > 180.0:
                lon_vals_for_mask = np.where(
                    lon_vals_for_mask > 180.0,
                    lon_vals_for_mask - 360.0,
                    lon_vals_for_mask,
                )
                if debug:
                    print("[DEBUG] generateLandSeaMask: converted lon to -180/180 for regionmask lookup")

            lon_for_mask = xr.DataArray(
                lon_vals_for_mask,
                dims=lon.dims,
                coords=lon.coords,
                attrs=lon.attrs,
            )

            try:
                land = _regionmask.defined_regions.natural_earth.land_110
            except AttributeError:
                land = _regionmask.defined_regions.natural_earth_v5_0_0.land_110

            if lat.ndim == 1 and lon.ndim == 1:
                if debug:
                    print("[DEBUG] generateLandSeaMask: applying regionmask (1D)")

                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        "ignore",
                        message="No gridpoint belongs to any region.*",
                        category=UserWarning,
                    )
                    raw_mask = land.mask(lon_for_mask, lat)

                expected_shape = (lat.size, lon.size)

            elif lat.ndim == 2 and lon.ndim == 2:
                if debug:
                    print("[DEBUG] generateLandSeaMask: applying regionmask (2D)")

                with warnings.catch_warnings():
                    warnings.filterwarnings(
                        "ignore",
                        message="No gridpoint belongs to any region.*",
                        category=UserWarning,
                    )
                    raw_mask = land.mask(lon_for_mask, lat)

                expected_shape = lat.shape

            else:
                raise ValueError(
                    f"Unsupported lat/lon dimensionality: "
                    f"{lat_name}.ndim={lat.ndim}, "
                    f"{lon_name}.ndim={lon.ndim}"
                )

            if debug:
                raw_vals = np.asarray(raw_mask.values, dtype=float)
                print(
                    "[DEBUG] generateLandSeaMask: raw_mask shape:",
                    raw_mask.shape,
                    "finite count:",
                    np.isfinite(raw_vals).sum(),
                )
                if np.any(np.isfinite(raw_vals)):
                    print(
                        "[DEBUG] generateLandSeaMask: raw_mask min/max:",
                        np.nanmin(raw_vals),
                        np.nanmax(raw_vals),
                    )

            # regionmask convention:
            #   finite value = land polygon ID
            #   NaN          = ocean / outside land polygons
            #
            # ENSO_metrics final convention:
            #   1.0 = land
            #   0.0 = ocean
            land01 = xr.where(np.isfinite(raw_mask), 1.0, 0.0).rename("sftlf")

            # Restore original input coordinates.
            # lon_for_mask may use -180–180 only for regionmask polygon lookup;
            # the returned mask must stay on the same grid/convention as the input field.
            land01 = land01.assign_coords({lat_name: lat, lon_name: lon})

            vals = np.asarray(land01.values, dtype=float)

            if np.all(~np.isfinite(vals)):
                raise RuntimeError(
                    "Generated land-sea mask is all NaN after conversion. "
                    "This should not happen; check regionmask conversion logic."
                )

            if debug:
                print(
                    "[DEBUG] generateLandSeaMask: returned lon min/max:",
                    np.nanmin(land01[lon_name].values),
                    np.nanmax(land01[lon_name].values),
                )
                print(
                    "[DEBUG] generateLandSeaMask: land01 shape:",
                    land01.shape,
                    "min:",
                    np.nanmin(vals),
                    "max:",
                    np.nanmax(vals),
                )

            if np.nanmax(vals) <= 0:
                if debug:
                    print("[DEBUG] generateLandSeaMask: WARNING - all ocean mask")
                warnings.warn(
                    "Generated land mask is entirely ocean for this region. "
                    "This may be expected if the region is all ocean.",
                    stacklevel=2,
                )

            if np.nanmin(vals) >= 1:
                if debug:
                    print("[DEBUG] generateLandSeaMask: WARNING - all land mask")
                warnings.warn(
                    "Generated land mask is entirely land for this region. "
                    "This may be expected if the region is land-only; otherwise check "
                    "lat/lon detection and regionmask behavior.",
                    stacklevel=2,
                )

            if land01.shape != expected_shape:
                if debug:
                    print("[DEBUG] generateLandSeaMask: ERROR - mask shape mismatch")
                raise RuntimeError(
                    f"Generated mask shape mismatch: "
                    f"mask={land01.shape}, "
                    f"expected={expected_shape}, "
                    f"data={da.shape}"
                )

            land01.attrs.update(
                {
                    "long_name": "estimated land fraction",
                    "standard_name": "land_area_fraction",
                    "units": "1",
                    "comment": (
                        "Estimated from Natural Earth land polygons "
                        "using regionmask; "
                        "1=land, 0=ocean. "
                        "Longitude was converted only internally for polygon lookup; "
                        "returned coordinates preserve the input grid convention. "
                        "Prefer native sftlf when available."
                    ),
                }
            )

            return _finalize_cdat(
                land01,
                varname="sftlf",
                context="generateLandSeaMask",
            )

        except Exception as e:
            raise RuntimeError(
                "Land-sea mask generation failed inside regionmask. "
                f"Original error: {type(e).__name__}: {e}. "
                f"dims={getattr(da, 'dims', None)}, "
                f"coords={list(getattr(da, 'coords', []))}."
            ) from e

    class times:
        @staticmethod
        def Seasons(season_str: str):
            return _SeasonHelper(_MONTH_MAP.get(season_str, []))

# Attach season constants to the stub
for _s in [
    "JAN","FEB","MAR","APR","MAY","JUN",
    "JUL","AUG","SEP","OCT","NOV","DEC",
    "MAM","JJA","SON","DJF"
    ]:
    setattr(_CdutilAverager, _s, _SeasonHelper(_MONTH_MAP[_s]))

cdutil = _CdutilAverager()


# Normalize legacy ENSO_metrics/CDAT-style regridding method names
# to the xESMF method names used by the modern backend.
_REGRID_METHOD_MAP = {
    "linear":       "bilinear",
    "bilinear":     "bilinear",
    "conserve":     "conservative",
    "conservative": "conservative",
    "nearest":      "nearest_s2d",
    "nearest_s2d":  "nearest_s2d",
    "patch":        "patch",
}


class REGRID2horizontal__Horizontal:
    """
    Compatibility replacement for the legacy ``regrid2.horizontal.Horizontal``
    interface used by ENSO_metrics.

    The modern implementation uses xESMF/ESMF when available, with a scipy-based
    fallback for supported rectilinear-grid cases.
    """
    def __init__(self, src_grid, dst_grid, method="bilinear"):
        self._src = src_grid
        self._dst = dst_grid
        self._method = _REGRID_METHOD_MAP.get(str(method).lower(), "bilinear")

    def __call__(self, tab):
        if not _HAS_XESMF:
            raise ImportError(
                "xesmf is required for regrid2-style horizontal regridding. "
                "Install with:  conda install -c conda-forge xesmf"
            )
        da = cdat_to_da(tab, name=getattr(tab, 'id', 'var'))
        if _is_unstructured_grid(da):
            raise NotImplementedError(
                "REGRID2horizontal: native unstructured/MPAS grids are not supported. "
                "Remap to a regular lat-lon grid first."
            )
        dst_lat = np.asarray(self._dst.getLatitude()[:])
        dst_lon = np.asarray(self._dst.getLongitude()[:])
        # Identity check — skip regridding if source and target grids are identical
        src_lat_dim = _validate_axis(da, "Y", context="latitude selection")
        src_lon_dim = _validate_axis(da, "X", context="longitude selection")
        if (src_lat_dim and src_lon_dim
                and np.array_equal(np.asarray(da[src_lat_dim]), dst_lat)
                and np.array_equal(np.asarray(da[src_lon_dim]), dst_lon)):
            return tab
        target_ds = xr.Dataset(coords={"lat": dst_lat, "lon": dst_lon})
        src_da = da.rename({src_lat_dim: 'lat', src_lon_dim: 'lon'}) if src_lat_dim and src_lon_dim else da
        # reuse_weights requires a pre-computed weight file; without one
        # xesmf raises "To reuse weights, you need to provide either filename
        # or weights".  Always recompute weights — they are fast for the small
        # output grids used in ENSO metrics.
        regridder = _xesmf.Regridder(
            src_da.to_dataset(name='var'),
            target_ds,
            method=self._method,
            extrap_method="nearest_s2d",
            reuse_weights=False,
        )
        result = regridder(da)
        valid_src = xr.where(np.isfinite(src_da), 1.0, 0.0)
        valid_regridder = _xesmf.Regridder(
            valid_src.to_dataset(name='valid'),
            target_ds,
            method=self._method,
            reuse_weights=False,
        )
        valid_fraction = valid_regridder(valid_src)
        result = result.where(np.isfinite(valid_fraction) & (valid_fraction >= 1.0 - 1e-12))
        # xESMF strips coordinate attributes from non-spatial dimensions (time, lev,
        # etc.).  Restore them from the source DataArray so that _finalize_cdat /
        # da_to_cdat can still detect axis types.
        for _dim in list(src_da.dims):
            if _dim in ("lat", "lon"):
                continue  # spatial dims were replaced by target grid — skip
            if _dim in result.coords and _dim in src_da.coords:
                _orig_attrs = src_da.coords[_dim].attrs
                if _orig_attrs and not result.coords[_dim].attrs:
                    result = result.assign_coords(
                        {_dim: result.coords[_dim].assign_attrs(_orig_attrs)}
                    )
        # Constant-field preservation: if source is spatially uniform, fill result to
        # that constant to avoid interpolation artefacts / numerical drift
        data_flat = _mv(tab)
        if data_flat.ndim >= 2 and not np.any(ma.getmaskarray(data_flat)):
            spatial = data_flat.reshape(data_flat.shape[:-2] + (-1,))
            if np.allclose(spatial, spatial[..., :1], atol=1e-8):
                result_cdat = _finalize_cdat(
                    result, varname=getattr(tab, 'id', 'var'),
                    context="REGRID2horizontal:constant",
                    require_time=_has_time_axis(tab)
                )
                result_cdat._data[:] = data_flat.flat[0]
                return result_cdat
        return _finalize_cdat(
            result, varname=getattr(tab, 'id', 'var'),
            context="REGRID2horizontal",
            require_time=_has_time_axis(tab)
            )


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of simple uvcdat functions used in EnsoMetricsLib.py
#
def ArrayOnes(tab, id='new_variable_ones'):
    """
    #################################################################################
    Description:
    Create a masked_array filled with ones with the same properties as tab (shape, axes, grid, mask)
    #################################################################################
    """
    tab = _to_cdat(tab)
    return create_variable(MV2ones(tab.shape), axes=tab.getAxisList(), grid=tab.getGrid(), mask=tab.mask, id=id)


def ArrayZeros(tab, id='new_variable_zeros'):
    """
    #################################################################################
    Description:
    Create a masked_array filled with zeros with the same properties as tab (shape, axes, grid, mask)
    #################################################################################
    """
    tab = _to_cdat(tab)
    return create_variable(MV2zeros(tab.shape), axes=tab.getAxisList(), grid=tab.getGrid(), mask=tab.mask, id=id)


def _variable_label(tab):
    attrs = getattr(tab, "attributes", {}) or {}
    for candidate in [
        getattr(tab, "name", None),
        getattr(tab, "id", None),
        attrs.get("id"),
        attrs.get("name"),
        attrs.get("variable_id"),
        attrs.get("short_name"),
        attrs.get("standard_name"),
        attrs.get("long_name"),
        attrs.get("source_id"),
    ]:
        if candidate is None:
            continue
        candidate = str(candidate)
        if candidate and candidate not in ["?", "var", "variable"]:
            return candidate
    return "unknown"


def _make_coslat_areacell(tab):
    """Build a cosine-latitude area-weight CDATVariable matching *tab*'s grid.

    Used as a guaranteed fallback when the caller did not supply an explicit
    areacell (e.g. observation files without areacella).  This is more robust
    than passing areacell=None to cdutil.averager, which can silently return
    None for certain grid configurations.

    Returns a CDATVariable with the same lat/lon axes as *tab* and values
    proportional to cos(lat), broadcast to match the full spatial shape.
    Returns None if *tab* has no latitude axis.
    """
    tab = _to_cdat(tab)
    lat_ax = tab.getLatitude()
    if lat_ax is None:
        return None
    lon_ax = tab.getLongitude()
    lat_vals = np.asarray(lat_ax[:], dtype=float)
    w_lat = ma.masked_invalid(np.cos(np.deg2rad(lat_vals)))  # (nlat,)
    if lon_ax is not None:
        # broadcast (nlat,) -> (nlat, nlon)
        w_2d = np.tile(w_lat[:, np.newaxis], (1, len(lon_ax)))
        axes = [lat_ax, lon_ax]
    else:
        w_2d = w_lat
        axes = [lat_ax]
    _var_name = _variable_label(tab)
    warnings.warn(
        f"areacell is None for variable {_var_name!r}; "
        "synthesising cosine-latitude weights. "
        "Provide areacella/sftlf for accurate spatial averages.",
        stacklevel=3,
    )
    return create_variable(
        ma.masked_invalid(w_2d.astype(float)),
        axes=axes,
        id="areacell_coslat",
    )


def AverageHorizontal(tab, areacell=None, region=None, **kwargs):
    """
    #################################################################################
    Description:
    Averages along 'xy' axis
    #################################################################################
    """
    keyerror = None
    tab = _to_cdat(tab)
    lat_num = get_num_axis(tab, "latitude")
    lon_num = get_num_axis(tab, "longitude")
    # Guard: int() ensures str() never produces float repr (e.g. "1.02.0")
    # if a numpy scalar ever leaks through get_num_axis.
    snum = str(int(lat_num)) + str(int(lon_num))
    _tab_grid = tab.getGrid()
    _area_grid = areacell.getGrid() if areacell is not None else None
    if areacell is None or _tab_grid is None or _area_grid is None or _tab_grid.shape != _area_grid.shape:
        if areacell is not None and _tab_grid is not None and _area_grid is not None \
                and _tab_grid.shape != _area_grid.shape:
            print("\033[93m" + str().ljust(15) + "EnsoUvcdatToolsLib AverageHorizontal" + "\033[0m")
            print(
                "\033[93m" + str().ljust(25) + "tab.grid " + str(_tab_grid.shape) +
                " is not the same as areacell.grid " + str(_area_grid.shape) + " \033[0m"
            )
        areacell = _make_coslat_areacell(tab)
    if areacell is not None:
        averaged_tab = MV2multiply(tab, areacell)
        # Sum in reverse index order so removing a higher-indexed axis first
        # does not shift the remaining lower index before it is used.
        for ax in sorted([int(lat_num), int(lon_num)], reverse=True):
            averaged_tab = MV2sum(averaged_tab, axis=ax)
        averaged_tab = averaged_tab / _valid_weight_sum(
            tab, areacell, (int(lat_num), int(lon_num))
        )
        averaged_tab = _apply_nan_majority_policy(averaged_tab, tab, (int(lat_num), int(lon_num)))
    else:
        # No latitude axis at all — last resort fallback
        try:
            averaged_tab = cdutil.averager(tab, axis="xy", weights="weighted", action="average")
        except Exception:
            try:
                averaged_tab = cdutil.averager(tab, axis=snum, weights="weighted", action="average")
            except Exception:
                keyerror = "cannot perform horizontal average: no latitude axis and cdutil fallback failed"
                averaged_tab = None
                list_strings = [
                    "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": horizontal average",
                    str().ljust(5) + keyerror]
                EnsoErrorsWarnings.my_warning(list_strings)
    # Fail-fast: if result is still None, set a keyerror so the caller knows
    if averaged_tab is None and keyerror is None:
        keyerror = "AverageHorizontal returned None — check grid, weights, and axis metadata"
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": horizontal average",
            str().ljust(5) + keyerror]
        EnsoErrorsWarnings.my_warning(list_strings)
    if averaged_tab is not None:
        averaged_tab = _finalize_existing_cdat(
            averaged_tab, context="AverageHorizontal",
            require_time=_has_time_axis(tab),
        )
    return averaged_tab, keyerror


def AverageMeridional(tab, areacell=None, region=None, **kwargs):
    """
    #################################################################################
    Description:
    Average along the latitude / meridional axis.

    This modernized implementation preserves the legacy ENSO_metrics interface
    while using CDAT-like compatibility objects. If an areacell field is missing
    or incompatible with the input grid, a cosine-latitude area proxy is
    synthesized so the meridional average can still be computed on regular
    rectilinear grids.
    #################################################################################
    """
    keyerror = None

    tab = _to_cdat(tab)

    var_label = _variable_label(tab)

    try:
        lat_num = get_num_axis(tab, "latitude")
        lon_num = get_num_axis(tab, "longitude")
        snum = str(int(lat_num))
    except Exception as e:
        keyerror = (
            f"cannot determine latitude/longitude axis for meridional average "
            f"on variable {var_label!r}: {e}"
        )
        averaged_tab = None
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": meridional average",
            str().ljust(5) + keyerror,
        ]
        EnsoErrorsWarnings.my_warning(list_strings)
        return averaged_tab, keyerror

    _tab_grid = tab.getGrid()
    _area_grid = areacell.getGrid() if areacell is not None else None

    # Synthesise cosine-latitude weights when areacell is absent or on a
    # different grid. This avoids silent failures and preserves legacy behavior
    # for regular rectilinear grids.
    if (
        areacell is None
        or _tab_grid is None
        or _area_grid is None
        or _tab_grid.shape != _area_grid.shape
    ):
        if (
            areacell is not None
            and _tab_grid is not None
            and _area_grid is not None
            and _tab_grid.shape != _area_grid.shape
        ):
            print("\033[93m" + str().ljust(15) + "EnsoUvcdatToolsLib AverageMeridional" + "\033[0m")
            print(
                "\033[93m" + str().ljust(25)
                + "tab.grid " + str(_tab_grid.shape)
                + " is not the same as areacell.grid " + str(_area_grid.shape)
                + " \033[0m"
            )

        try:
            areacell = _make_coslat_areacell(tab)
        except Exception as e:
            areacell = None
            warnings.warn(
                f"AverageMeridional: failed to synthesize cosine-latitude "
                f"areacell for variable {var_label!r} ({e}); "
                "falling back to cdutil-style weighted average if available.",
                stacklevel=2,
            )

    if areacell is not None:
        try:
            averaged_tab = MV2multiply(tab, areacell)
            averaged_tab = (
                MV2sum(averaged_tab, axis=int(lat_num))
                / _valid_weight_sum(tab, areacell, int(lat_num))
            )
            averaged_tab = _apply_nan_majority_policy(averaged_tab, tab, int(lat_num))

        except Exception as e:
            keyerror = (
                f"cannot perform meridional average with areacell for variable "
                f"{var_label!r}: {e}"
            )
            averaged_tab = None
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": meridional average",
                str().ljust(5) + keyerror,
            ]
            EnsoErrorsWarnings.my_warning(list_strings)

    else:
        try:
            averaged_tab = cdutil.averager(
                tab,
                axis="y",
                weights="weighted",
                action="average",
            )
        except Exception:
            try:
                averaged_tab = cdutil.averager(
                    tab,
                    axis=snum,
                    weights="weighted",
                    action="average",
                )
            except Exception:
                keyerror = (
                    f"cannot perform meridional average for variable "
                    f"{var_label!r}: no compatible areacell and cdutil fallback failed"
                )
                averaged_tab = None
                list_strings = [
                    "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": meridional average",
                    str().ljust(5) + keyerror,
                ]
                EnsoErrorsWarnings.my_warning(list_strings)

    # Fail-fast: if result is still None, set a keyerror.
    if averaged_tab is None and keyerror is None:
        keyerror = (
            f"AverageMeridional returned None for variable {var_label!r}; "
            "check grid, weights, and axis metadata"
        )
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": meridional average",
            str().ljust(5) + keyerror,
        ]
        EnsoErrorsWarnings.my_warning(list_strings)

    # Preserve a 1-D longitude axis after averaging over latitude, matching the
    # legacy CDAT behavior for curvilinear / 2-D longitude inputs.
    if averaged_tab is not None:
        lon = tab.getLongitude()
        if lon is not None and len(lon.shape) > 1:
            lonn = create_axis(
                MV2array(lon[0, :]),
                id="longitude",
                units=getattr(lon, "units", "degrees_east"),
                attributes={
                    "axis": "X",
                    "standard_name": "longitude",
                },
            )

            lon_num = get_num_axis(tab, "longitude")
            try:
                averaged_tab.setAxis(lon_num, lonn)
            except Exception:
                averaged_tab.setAxis(lon_num - 1, lonn)

    if averaged_tab is not None:
        averaged_tab = _finalize_existing_cdat(
            averaged_tab,
            context="AverageMeridional",
            require_time=_has_time_axis(tab),
        )

    return averaged_tab, keyerror


def AverageTemporal(tab, areacell=None, **kwargs):
    """
    Averages along the time axis.
    """
    keyerror = None
    tab = _to_cdat(tab)

    if not _has_time_axis(tab):
        keyerror = "cannot perform temporal average: no time axis"
        averaged_tab = None
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": temporal average",
            str().ljust(5) + keyerror,
        ]
        EnsoErrorsWarnings.my_warning(list_strings)
        return averaged_tab, keyerror

    try:
        averaged_tab = cdutil.averager(tab, axis="t")
    except Exception:
        try:
            time_num = _get_time_axis_index(tab, "AverageTemporal")
            averaged_tab = cdutil.averager(tab, axis=str(time_num))
        except Exception:
            keyerror = "cannot perform temporal average"
            averaged_tab = None
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": temporal average",
                str().ljust(5) + keyerror,
            ]
            EnsoErrorsWarnings.my_warning(list_strings)

    if averaged_tab is not None:
        averaged_tab = _apply_nan_majority_policy(averaged_tab, tab, "t")
        # Temporal averaging intentionally removes time, so do not require T.
        averaged_tab = _finalize_existing_cdat(averaged_tab, context="AverageTemporal")
    return averaged_tab, keyerror


def AverageZonal(tab, areacell=None, region=None, **kwargs):
    """
    #################################################################################
    Description:
    Averages along 'x' axis
    #################################################################################
    """
    keyerror = None
    lat_num = get_num_axis(tab, "latitude")
    lon_num = get_num_axis(tab, "longitude")
    snum = str(int(lon_num))
    tab = _to_cdat(tab)
    _tab_grid = tab.getGrid()
    _area_grid = areacell.getGrid() if areacell is not None else None
    # Synthesise cosine-latitude weights when areacell is absent or on a
    # different grid — prevents the silent None return from cdutil.averager.
    if areacell is None or _tab_grid is None or _area_grid is None or _tab_grid.shape != _area_grid.shape:
        if areacell is not None and _tab_grid is not None and _area_grid is not None \
                and _tab_grid.shape != _area_grid.shape:
            print("\033[93m" + str().ljust(15) + "EnsoUvcdatToolsLib AverageZonal" + "\033[0m")
            print(
                "\033[93m" + str().ljust(25) + "tab.grid " + str(_tab_grid.shape) +
                " is not the same as areacell.grid " + str(_area_grid.shape) + " \033[0m"
            )
        areacell = _make_coslat_areacell(tab)
    if areacell is not None:
        averaged_tab = MV2multiply(tab, areacell)
        averaged_tab = MV2sum(averaged_tab, axis=int(lon_num)) / _valid_weight_sum(
            tab, areacell, int(lon_num)
        )
        averaged_tab = _apply_nan_majority_policy(averaged_tab, tab, int(lon_num))
    else:
        try:
            averaged_tab = cdutil.averager(tab, axis="x", weights="weighted", action="average")
        except Exception:
            try:
                averaged_tab = cdutil.averager(tab, axis=snum, weights="weighted", action="average")
            except Exception:
                keyerror = "cannot perform zonal average: no latitude axis and cdutil fallback failed"
                averaged_tab = None
                list_strings = [
                    "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": zonal average",
                    str().ljust(5) + keyerror
                ]
                EnsoErrorsWarnings.my_warning(list_strings)
    # Fail-fast: if result is still None, set a keyerror
    if averaged_tab is None and keyerror is None:
        keyerror = "AverageZonal returned None — check grid, weights, and axis metadata"
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": zonal average",
                        str().ljust(5) + keyerror]
        EnsoErrorsWarnings.my_warning(list_strings)
    if averaged_tab is not None:
        lat = tab.getLatitude()
        if lat is not None and len(lat.shape) > 1:
            latn = create_axis(MV2array(lat[:, 0]), id="latitude")
            latn.units = lat.units
            lat_num = get_num_axis(tab, "latitude")
            try:
                averaged_tab.setAxis(lat_num, latn)
            except Exception:
                averaged_tab.setAxis(lat_num - 1, latn)
    if averaged_tab is not None:
        averaged_tab = _finalize_existing_cdat(
            averaged_tab, context="AverageZonal",
            require_time=_has_time_axis(tab),
        )
    return averaged_tab, keyerror


# Dictionary of averaging methods
dict_average = {'horizontal': AverageHorizontal, 'meridional': AverageMeridional, 'time': AverageTemporal,
                'zonal': AverageZonal}


def Concatenate(tab1, tab2, events1=[], events2=[]):
    my_events = events1 + events2
    if len(my_events) > 0:
        my_events_sort = sorted(my_events)
        tab_out = None
        for yy in my_events_sort:
            if tab_out is None:
                if yy in events1:
                    tab_out = MV2array([tab1[events1.index(yy)]])
                else:
                    tab_out = MV2array([tab2[events2.index(yy)]])
            else:
                if yy in events1:
                    tab_out = MV2concatenate((tab_out, MV2array([tab1[events1.index(yy)]])))
                else:
                    tab_out = MV2concatenate((tab_out, MV2array([tab2[events2.index(yy)]])))
        axes = create_axis(MV2array(my_events_sort, dtype="int32"), id="years")
        if len(events1):
            tmp = copy.copy(tab1)
        else:
            tmp = copy.copy(tab2)
        att = tmp.attributes
        if len(tmp.shape) > 1:
            mask = tmp[0].mask
            mask2 = MV2zeros(tab_out.shape)
            mask2[:] = mask
            dictvar = {"axes": [axes] + _to_cdat(tab1[0]).getAxisList(), "mask": mask2, "grid": tmp.getGrid(), "attributes": att}
        else:
            dictvar = {"axes": [axes], "attributes": att}
        tab_out = create_variable(tab_out, **dictvar)
    else:
        tab_out = MyEmpty(tab1[:5, 0], time=True, time_id="years")
    return tab_out


def closest_grid(region, nlat, nlon):
    res = [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75]
    region_ref = ReferenceRegions(region)
    lats = region_ref["latitude"]
    dy = float(abs(max(lats) - min(lats))) / nlat
    lyy = [abs(dy - ii) for ii in res]
    lyy = res[lyy.index(min(lyy))]
    lons = region_ref["longitude"]
    dx = float(abs(max(lons) - min(lons))) / nlon
    lxx = [abs(dx - ii) for ii in res]
    lxx = res[lxx.index(min(lxx))]
    if lxx == lyy:
        grid = "generic_" + str(lxx) + "x" + str(lxx) + "deg"
    else:
        dx = abs(lxx + lyy) / 2.
        lxx = [abs(dx - ii) for ii in res]
        lxx = res[lxx.index(min(lxx))]
        grid = "generic_" + str(lxx) + "x" + str(lxx) + "deg"
    return grid


def ComputeInterannualAnomalies(tab):
    """
    #################################################################################
    Description:
    Computes interannual anomalies
    #################################################################################
    """
    return cdutil.ANNUALCYCLE.departures(tab)


def Correlation(tab, ref, weights=None, axis=0, centered=1, biased=1):
    """
    #################################################################################
    Description:
    Computes correlation
    #################################################################################
    """
    return GENUTILcorrelation(tab, ref, weights=weights, axis=axis, centered=centered, biased=biased)


def OperationAdd(tab, number_or_tab):
    """
    #################################################################################
    Description:
    Adds every elements of 'tab' by 'number_or_tab'
    If 'number_or_tab' is an array it must have the same shape as tab
    #################################################################################
    """
    if not isinstance(number_or_tab, int) and not isinstance(number_or_tab, float):
        if tab.shape != number_or_tab.shape:
            EnsoErrorsWarnings.mismatch_shapes_error(tab, number_or_tab, INSPECTstack())
    return MV2add(tab, number_or_tab)


def OperationDivide(tab, number_or_tab):
    """
    #################################################################################
    Description:
    Divides every elements of 'tab' by 'number_or_tab'
    #################################################################################
    """
    if not isinstance(number_or_tab, int) and not isinstance(number_or_tab, float):
        if tab.shape != number_or_tab.shape:
            EnsoErrorsWarnings.mismatch_shapes_error(tab, number_or_tab, INSPECTstack())
    return MV2divide(tab, number_or_tab)


def OperationMultiply(tab, number_or_tab):
    """
    #################################################################################
    Description:
    Multiplies every elements of 'tab' by 'number_or_tab'
    #################################################################################
    """
    tab = _to_cdat(tab)
    if not isinstance(number_or_tab, int) and not isinstance(number_or_tab, float):
        if tab.shape != number_or_tab.shape:
            EnsoErrorsWarnings.mismatch_shapes_error(tab, number_or_tab, INSPECTstack())
    tab_out = MV2multiply(tab, number_or_tab)
    axes = tab.getAxisList()
    att = tab.attributes
    if len(tab.shape) > 1:
        dictvar = {"axes": axes, "mask": tab.mask, "grid": tab.getGrid(), "attributes": att}
    else:
        dictvar = {"axes": axes, "attributes": att}
    tab_out = create_variable(tab_out, **dictvar)
    return tab_out


def OperationSubtract(tab, number_or_tab):
    """
    #################################################################################
    Description:
    Subtracts every elements of 'tab' by 'number_or_tab'
    #################################################################################
    """
    if not isinstance(number_or_tab, int) and not isinstance(number_or_tab, float):
        if tab.shape != number_or_tab.shape:
            EnsoErrorsWarnings.mismatch_shapes_error(tab, number_or_tab, INSPECTstack())
    return MV2subtract(tab, number_or_tab)


# Dictionary of operations
dict_operations = {
    "divide": OperationDivide,
    "minus": OperationSubtract,
    "multiply": OperationMultiply,
    "plus": OperationAdd
    }


def RmsAxis(tab, ref, weights=None, axis=0, centered=0, biased=1):
    """
    #################################################################################
    Description:
    Compute the root-mean-square difference between ``tab`` and ``ref`` along
    a selected axis.

    This function preserves the legacy ENSO_metrics interface while using the
    modern numpy/compatibility-layer statistics pathway.
    #################################################################################

    :param tab: CDATVariable or masked-array-like
        Input field, usually the model field.

    :param ref: CDATVariable or masked-array-like
        Reference field, usually the observational field.

    :param weights: array-like or string, optional
        Weights used in the RMS calculation. If ``"weighted"``, latitude-based
        weights are used where appropriate.
        default value = None

    :param axis: int or string, optional
        Axis over which to compute the RMS difference. May be an integer axis
        index or a CDAT-style axis specifier such as ``"x"``, ``"y"``,
        ``"t"``, or ``"xy"``.
        default value = 0

    :param centered: integer, optional
        Legacy centered flag passed to the RMS calculation.
        ``0`` means the mean difference is not removed before computing RMS;
        ``1`` means the mean difference is removed first.
        default value = 0

    :param biased: integer, optional
        Legacy normalization flag passed to the RMS calculation.
        ``1`` uses biased normalization; ``0`` uses unbiased normalization.
        default value = 1

    :return rmse, keyerror:
        RMS difference along the requested axis and any accumulated keyerror
        message.
    """
    tab = _to_cdat(tab)
    ref = _to_cdat(ref)
    keyerror = None
    # Computes the root mean square difference
    try:
        rmse = GENUTILrms(tab, ref, weights=weights, axis=axis, centered=centered, biased=biased)
    except Exception:
        keyerror = "cannot perform RMS along given axis: tab (" + str(tab.shape) + ") and ref (" + str(ref.shape) +\
            ") are not on the same grid"
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": RMS over axis " + str(axis),
                        str().ljust(5) + "cannot perform RMS along given axis",
                        str().ljust(10) + "axes may not be in the same order in 'ref' and 'tab'",
                        str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                        str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList())]
        EnsoErrorsWarnings.my_warning(list_strings)
    try:
        rmse = float(rmse)
    except Exception:
        rmse = None
    return rmse, keyerror


def RmsHorizontal(tab, ref, centered=0, biased=1):
    """
    #################################################################################
    Description:
    Compute the horizontal root-mean-square difference between ``tab`` and
    ``ref``.

    The RMS is computed over latitude and longitude dimensions using the modern
    numpy/compatibility-layer statistics pathway, while preserving the legacy
    ENSO_metrics integer-flag interface.
    #################################################################################

    :param tab: CDATVariable or masked-array-like
        Input field, usually the model field, with latitude/longitude axes.

    :param ref: CDATVariable or masked-array-like
        Reference field, usually the observational field, with
        latitude/longitude axes.

    :param centered: integer, optional
        Legacy centered flag passed to the RMS calculation.
        ``0`` means the mean difference is not removed before computing RMS;
        ``1`` means the mean difference is removed first.
        default value = 0

    :param biased: integer, optional
        Legacy normalization flag passed to the RMS calculation.
        ``1`` uses biased normalization; ``0`` uses unbiased normalization.
        default value = 1

    :return rmse, keyerror:
        Horizontal RMS difference and any accumulated keyerror message.
    """
    tab = _to_cdat(tab)
    ref = _to_cdat(ref)
    keyerror = None
    # Computes the root mean square difference
    try:
        rmse = GENUTILrms(tab, ref, weights="weighted", axis="xy", centered=centered, biased=biased)
    except Exception:
        lat_num = get_num_axis(tab, "latitude")
        lon_num = get_num_axis(tab, "longitude")
        try:
            rmse = GENUTILrms(
                tab, ref, weights="weighted", axis=str(lat_num)+str(lon_num), centered=centered,
                biased=biased
            )
        except Exception:
            keyerror = "cannot perform horizontal RMS (x=" + str(lon_num) + ", y=" + str(lat_num) + "): tab (" +\
                str(tab.shape) + ") and ref (" + str(ref.shape) + ") are not on the same grid"
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": horizontal RMS",
                str().ljust(5) + "cannot perform horizontal RMS",
                str().ljust(10) + "either lat and lon cannot be found in 'ref' / 'tab'",
                str().ljust(10) + "or lat and lon are not in the same order in 'ref' and 'tab'",
                str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList())]
            EnsoErrorsWarnings.my_warning(list_strings)
    try:
        rmse = float(rmse)
    except Exception:
        rmse = None
    return rmse, keyerror


def RmsMeridional(tab, ref, centered=0, biased=1):
    """
    #################################################################################
    Description:
    Compute the meridional root-mean-square difference between ``tab`` and
    ``ref``.

    The RMS is computed over the latitude dimension using the modern
    numpy/compatibility-layer statistics pathway, while preserving the legacy
    ENSO_metrics integer-flag interface.
    #################################################################################

    :param tab: CDATVariable or masked-array-like
        Input field, usually the model field, with a latitude axis.

    :param ref: CDATVariable or masked-array-like
        Reference field, usually the observational field, with a latitude axis.

    :param centered: integer, optional
        Legacy centered flag passed to the RMS calculation.
        ``0`` means the mean difference is not removed before computing RMS;
        ``1`` means the mean difference is removed first.
        default value = 0

    :param biased: integer, optional
        Legacy normalization flag passed to the RMS calculation.
        ``1`` uses biased normalization; ``0`` uses unbiased normalization.
        default value = 1

    :return rmse, keyerror:
        Meridional RMS difference and any accumulated keyerror message.
    """
    tab = _to_cdat(tab)
    ref = _to_cdat(ref)
    keyerror = None
    # Computes the root mean square difference
    try:
        rmse = GENUTILrms(tab, ref, axis="y", centered=centered, biased=biased)
    except Exception:
        lat_num = get_num_axis(tab, "latitude")
        try:
            rmse = GENUTILrms(tab, ref, axis=str(lat_num), centered=centered, biased=biased)
        except Exception:
            keyerror = "cannot perform meridional RMS (y=" + str(lat_num) + "): tab (" + str(tab.shape) +\
                ") and ref (" + str(ref.shape) + ") are not on the same grid"
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": meridional RMS",
                str().ljust(5) + "cannot perform meridional RMS",
                str().ljust(10) + "lat cannot be found in 'ref' / 'tab'",
                str().ljust(10) + "or lat is not in the same order in 'ref' and 'tab'",
                str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList())]
            EnsoErrorsWarnings.my_warning(list_strings)
    try:
        rmse = float(rmse)
    except Exception:
        rmse = None
    return rmse, keyerror


def RmsTemporal(tab, ref, centered=0, biased=1):
    """
    #################################################################################
    Description:
    Compute the temporal root-mean-square difference between ``tab`` and
    ``ref``.

    The RMS is computed over the time dimension using the modern
    numpy/compatibility-layer statistics pathway, while preserving the legacy
    ENSO_metrics integer-flag interface.
    #################################################################################

    :param tab: CDATVariable or masked-array-like
        Input time series or time-dependent field, usually the model field.

    :param ref: CDATVariable or masked-array-like
        Reference time series or time-dependent field, usually the
        observational field.

    :param centered: integer, optional
        Legacy centered flag passed to the RMS calculation.
        ``0`` means the mean difference is not removed before computing RMS;
        ``1`` means the mean difference is removed first.
        default value = 0

    :param biased: integer, optional
        Legacy normalization flag passed to the RMS calculation.
        ``1`` uses biased normalization; ``0`` uses unbiased normalization.
        default value = 1

    :return rmse, keyerror:
        Temporal RMS difference and any accumulated keyerror message.
    """
    tab = _to_cdat(tab)
    ref = _to_cdat(ref)
    keyerror = None
    rmse = None
    # Computes the root mean square difference
    try:
        rmse = GENUTILrms(tab, ref, axis="t", centered=centered, biased=biased)
    except Exception:
        try:
            time_num = _get_time_axis_index(tab, "RmsTemporal")
        except Exception:
            keyerror = "cannot determine time axis for temporal RMS: tab (" + str(tab.shape) + ")"
            return None, keyerror
        try:
            rmse = GENUTILrms(tab, ref, axis=str(time_num), centered=centered, biased=biased)
        except Exception:
            keyerror = "cannot perform temporal RMS (t=" + str(time_num) + "): tab (" + str(tab.shape) + \
                ") and ref (" + str(ref.shape) + ") are not on the same grid"
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": temporal RMS",
                str().ljust(5) + "cannot perform temporal RMS",
                str().ljust(10) + "time cannot be found in 'ref' / 'tab'",
                str().ljust(10) + "or time is not in the same order in 'ref' and 'tab'",
                str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList())]
            EnsoErrorsWarnings.my_warning(list_strings)
    try:
        rmse = float(rmse)
    except Exception:
        rmse = None
    return rmse, keyerror


def RmsZonal(tab, ref, centered=0, biased=1):
    """
    #################################################################################
    Description:
    Compute the zonal root-mean-square difference between ``tab`` and ``ref``.

    The RMS is computed over the longitude dimension using the modern
    numpy/compatibility-layer statistics pathway, while preserving the legacy
    ENSO_metrics integer-flag interface.
    #################################################################################

    :param tab: CDATVariable or masked-array-like
        Input field, usually the model field, with a longitude axis.

    :param ref: CDATVariable or masked-array-like
        Reference field, usually the observational field, with a longitude axis.

    :param centered: integer, optional
        Legacy centered flag passed to the RMS calculation.
        ``0`` means the mean difference is not removed before computing RMS;
        ``1`` means the mean difference is removed first.
        default value = 0

    :param biased: integer, optional
        Legacy normalization flag passed to the RMS calculation.
        ``1`` uses biased normalization; ``0`` uses unbiased normalization.
        default value = 1

    :return rmse, keyerror:
        Zonal RMS difference and any accumulated keyerror message.
    """
    tab = _to_cdat(tab)
    ref = _to_cdat(ref)
    keyerror = None
    # Computes the root mean square difference
    try:
        rmse = GENUTILrms(tab, ref, axis="x", centered=centered, biased=biased)
    except Exception:
        lon_num = get_num_axis(tab, "longitude")
        try:
            rmse = GENUTILrms(tab, ref, axis=str(lon_num), centered=centered, biased=biased)
        except Exception:
            keyerror = "cannot perform zonal RMS (t=" + str(lon_num) + "): tab (" + str(tab.shape) + \
                ") and ref (" + str(ref.shape) + ") are not on the same grid"
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": zonal RMS",
                str().ljust(5) + "cannot perform zonal RMS",
                str().ljust(10) + "lon cannot be found in 'ref' / 'tab'",
                str().ljust(10) + "or lon is not in the same order in 'ref' and 'tab'",
                str().ljust(15) + "order: ref = " + str(ref.getOrder()) + ", tab = " + str(tab.getOrder()),
                str().ljust(15) + "axes: ref = " + str(ref.getAxisList()) + ", tab = " + str(tab.getAxisList())]
            EnsoErrorsWarnings.my_warning(list_strings)
    try:
        rmse = float(rmse)
    except Exception:
        rmse = None
    return rmse, keyerror


# Dictionary of RMS methods
dict_rms = {"axis": RmsAxis, "horizontal": RmsHorizontal, "meridional": RmsMeridional, "time": RmsTemporal,
            "zonal": RmsZonal}


def Std(tab, weights=None, axis=0, centered=1, biased=1):
    """
    #################################################################################
    Description:
    Computes standard deviation
    #################################################################################
    """
    tmp = GENUTILstd(tab, weights=weights, axis=axis, centered=centered, biased=biased)
    try:
        tmp.setGrid(tab.getGrid())
    except Exception:
        pass
    return tmp


def SumAxis(tab, axis=None, fill_value=0, dtype=None):
    """
    #################################################################################
    Description:
    Sum values along the requested axis.

    This helper preserves the legacy ENSO_metrics interface while using the
    modern numpy/compatibility-layer masked-array pathway.
    #################################################################################

    :param tab: CDATVariable or masked-array-like
        Input field.

    :param axis: int or string, optional
        Axis over which to sum. May be an integer axis index or a CDAT-style
        axis specifier.
        default value = None

    :param fill_value: number, optional
        Fill value used for missing values during summation.
        default value = 0

    :param dtype: data-type, optional
        Optional output dtype.
        default value = None

    :return:
        Sum along the requested axis, following the existing ENSO_metrics
        return convention.
    """
    tab = _to_cdat(tab)
    keyerror = None
    try:
        sum_along_axis = MV2sum(tab, axis=axis, fill_value=fill_value, dtype=dtype)
    except Exception:
        keyerror = "cannot sum along given axis (" + str(axis) + "): tab (" + str(tab.shape) + ")"
        sum_along_axis = None
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": sum over axis " + str(axis),
            str().ljust(5) + "cannot perform sum along given axis",
            str().ljust(10) + "axes: " + str(tab.getAxisList()),

            str().ljust(15) + "axis = " + str(axis) + " ; fill_value = " + str(fill_value) + " ; dtype = " + str(dtype)]
        EnsoErrorsWarnings.my_warning(list_strings)
    return sum_along_axis, keyerror


def TimeBounds(tab):
    """
    Finds first and last dates of tab's time axis.

    Safe for debug use: returns (None, None) if no time axis exists.
    """
    return _safe_time_bounds_for_debug(tab)
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of more complex functions (based on uvcdat) used in EnsoMetricsLib.py
#
def annualcycle(tab):
    """
    Computes the annual cycle: climatological value of each calendar month.
    """
    tab = _to_cdat(tab)
    _require_time_axis(tab, "annualcycle")

    initorder = tab.getOrder()
    tab = tab.reorder("t...")
    axes = tab.getAxisList()
    time_ax = _get_component_time(tab, "annualcycle")

    months = MV2array([tt.month for tt in time_ax])
    cyc = []

    for ii in range(12):
        ids = MV2compress(months == (ii + 1), list(range(len(tab))))
        tmp = MV2take(tab, ids, axis=0)
        tmp = MV2average(tmp, axis=0)
        cyc.append(tmp)
        del tmp

    time = create_axis(list(range(12)), id="time", axis_type="T")
    moy = create_variable(
        MV2array(cyc),
        axes=[time] + axes[1:],
        grid=tab.getGrid(),
        attributes=tab.attributes,
    )
    moy = moy.reorder(initorder)

    time = create_axis(list(range(12)), id="months", axis_type="T")
    moy.setAxis(get_num_axis(moy, "time"), time)
    return moy


def ApplyLandmask(tab, landmask, maskland=True, maskocean=False):
    """
    #################################################################################
    Description:
    Applies the landmask on the given tab
        if maskland is True, mask where landmask==100
        if maskocean is True, mask where landmask==0
    #################################################################################

    :param tab: masked_array
    :param landmask: masked_array
    :param maskland: boolean, optional
        masks land points
        default value is True
    :param maskocean: boolean, optional
        masks ocean points
        default value is False

    :return: tab: masked_array
        masked_array where land points and/or ocean points are masked
    """
    keyerror = None
    if maskland is True or maskocean is True:
        tab = _to_cdat(tab)
        landmask = _to_cdat(landmask)
        _tg, _lg = tab.getGrid(), landmask.getGrid()
        if _tg is None or _lg is None or _tg.shape != _lg.shape:
            keyerror = "tab (" + str(_tg.shape if _tg is not None else None) + ") and landmask (" + \
                str(_lg.shape if _lg is not None else None) + ") are not on the same grid"
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": applying landmask",
                str().ljust(5) + keyerror, str().ljust(5) + "cannot apply landmask",
                str().ljust(5) + "this metric will be skipped"
            ]
            EnsoErrorsWarnings.my_warning(list_strings)
        else:
            landmask_nd = MV2zeros(tab.shape)
            if landmask_nd.shape == landmask.shape:
                landmask_nd = copy.copy(landmask)
            else:
                try:
                    landmask_nd[:] = landmask
                except Exception:
                    try:
                        landmask_nd[:, :] = landmask
                    except Exception:
                        keyerror = "ApplyLandmask: tab must be more than 4D and this is not taken into account yet (" +\
                            str(tab.shape) + ") and landmask (" + str(landmask.shape) + ")"
                        list_strings = [
                            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": landmask shape",
                            str().ljust(5) + keyerror, str().ljust(5) + "cannot reshape landmask"
                        ]
                        EnsoErrorsWarnings.my_warning(list_strings)
            if keyerror is None:
                tab = MV2masked_where(landmask_nd.mask, tab)
                # if land = 100 instead of 1, divides landmask by 100
                if MV2minimum(landmask_nd) == 0 and MV2maximum(landmask_nd) == 100:
                    landmask_nd = landmask_nd / 100.
                if maskland is True:
                    tab = MV2masked_where(landmask_nd != 0, tab)
                if maskocean is True:
                    tab = MV2masked_where(landmask_nd != 1, tab)
    return tab, keyerror


def ApplyLandmaskToArea(area, landmask, maskland=True, maskocean=False):
    """
    #################################################################################
    Description:
    Applies the landmask on the given tab
        if maskland is True, mask where landmask==1 and area=area*(1-landmask) (to weight island and coastal points)
        if maskocean is True, mask where landmask==0 and area=area*landmask (to weight island and coastal points)
    #################################################################################

    :param area: masked_array
        areacell
    :param landmask: masked_array
    :param maskland: boolean, optional
        masks land points and weights island and coastal points
        default value is True
    :param maskocean: boolean, optional
        masks ocean points and weights island and coastal points
        default value is False

    :return: tab: masked_array
        masked_array where land points and/or ocean points are masked
    """
    keyerror = None
    if maskland is True or maskocean is True:
        if area.getGrid().shape != landmask.getGrid().shape:
            keyerror = "ApplyLandmaskToArea: area (" + str(area.getGrid().shape) + ") and landmask (" +\
                str(landmask.getGrid().shape) + ") are not on the same grid"
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": applying landmask to areacell",
                str().ljust(5) + keyerror, str().ljust(5) + "cannot apply landmask to areacell"
            ]
            EnsoErrorsWarnings.my_warning(list_strings)
        if keyerror is None:
            # if land = 100 instead of 1, divides landmask by 100
            if MV2minimum(landmask) == 0 and MV2maximum(landmask) == 100:
                landmask = landmask / 100.
            area = MV2masked_where(landmask.mask, area)
            if maskland:
                area = MV2masked_where(landmask != 0, area)
                area = MV2multiply(area, 1-landmask)
            if maskocean:
                area = MV2masked_where(landmask != 1, area)
                area = MV2multiply(area, landmask)
    return area, keyerror


def ArrayListAx(tab, list1, ax_name_ax="", ax_long_name="", ax_ref=""):
    tab_out = MV2array(tab)
    ax = create_axis(list(range(len(list1))), id=ax_name_ax)
    ax.regions = str(list1)
    if len(ax_long_name) > 0:
        ax.long_name = ax_long_name
    if len(ax_ref) > 0:
        ax.reference = ax_ref
    tab_out.setAxis(0, ax)
    return tab_out


def ArrayToList(tab):
    mask_val = tab.mask
    tmp_mask = mask_val if hasattr(mask_val, '__len__') and not isinstance(mask_val, np.bool_) else [mask_val]
    if all(ii is False for ii in tmp_mask) is True or all(ii == False for ii in tmp_mask) == True:
        tmp = NParray(tab)
    else:
        tmp = NParray(MV2where(tab.mask, 1e20, tab))
    if len(tab.shape) == 1:
        tab_out = list(tmp)
    elif len(tab.shape) == 2:
        tab_out = [list(tmp[ii]) for ii in list(range(len(tab)))]
    else:
        tab_out = [None]
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": bad shape",
            str().ljust(5) + "cannot transform this array to a list",
            str().ljust(10) + "the length (" + str(len(tab.shape)) + ") of the shape (" + str(tab.shape) +
            ") is too large",
            str().ljust(10) + "it is not programed yet"
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    return tab_out


def BasinMask(
        tab_in, region_mask, box=None, lat1=None, lat2=None,
        latkey='', lon1=None, lon2=None, lonkey='',
        debug=False
    ):
    keyerror = None
    tab_in = _to_cdat(tab_in)
    keys = ["between", "outside"]
    # open file
    this_dir, this_filename = OSpath__split(__file__)
    # check basin file
    basin_generic_ncfile = OSpath__join(this_dir, "../share/EnsoMetrics/basin_generic_1x1deg.nc")
    if not OSpath__isfile(basin_generic_ncfile):
        basin_generic_ncfile = OSpath__join(SYS_prefix, "share", "EnsoMetrics", "basin_generic_1x1deg.nc")
    if debug is True:
        dict_debug = {
            "line1": "(path) " + str(this_dir), "line2": "(file) " + str(this_filename),
            "line3": "(basin) " + str(basin_generic_ncfile)
        }
        EnsoErrorsWarnings.debug_mode("\033[93m", "OSpath__split", 20, **dict_debug)
    ff = open_file(basin_generic_ncfile)
    # read basins
    if box is not None:
        region_ref = ReferenceRegions(box)
        basin = ff("basin", latitude=region_ref["latitude"], longitude=region_ref["longitude"])
    else:
        basin = ff("basin")
    if debug is True:
        dict_debug = {
            "axes1": str([ax.id for ax in basin.getAxisList()]), "shape1": str(basin.shape),
            "line1": "order = " + str(basin.getOrder())
        }
        EnsoErrorsWarnings.debug_mode("\033[93m", "in BasinMask", 20, **dict_debug)
    # choose basin
    keybasin = {"atlantic": 1, "pacific": 2, "indian": 3, "antarctic": 10, "arctic": 11}
    mask = MV2zeros(basin.shape)
    if region_mask.lower() not in list(keybasin.keys()):
        keyerror = "unknown region: " + region_mask + " (basin_generic_1x1deg.nc, regrided file from NOAA NODC" + \
            "WOA09 Masks basin Data Files)"
        list_strings = [
            "WARNING" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": region",
            str().ljust(5) + keyerror,
            str().ljust(5) + "https://iridl.ldeo.columbia.edu/SOURCES/.NOAA/.NODC/.WOA09/.Masks/.basin/"
            + "datafiles.html"
        ]
        EnsoErrorsWarnings.my_warning(list_strings)
    else:
        mask = MV2where(basin == keybasin[region_mask], 1, mask)
        mask = MV2where(basin.mask, 0, mask)
    # basin mask is selected only between or outside lat1 and lat2
    if latkey in keys and lat1 is not None and lat2 is not None:
        lat2d = MV2zeros(basin.shape)
        lat2d = lat2d.reorder("10")
        lat2d[:] = basin.getLatitude()
        lat2d = lat2d.reorder("10")
        tmp = MV2where(lat2d > lat1, 1, 0) + MV2where(lat2d < lat2, 1, 0)
        if latkey == "between":
            mask = MV2where(tmp != 2, 0, mask)
        else:
            mask = MV2where(tmp == 2, 0, mask)
    # basin mask is selected only between or outside lon1 and lon2
    if lonkey in keys and lat1 is not None and lat2 is not None:
        lon2d = MV2zeros(basin.shape)
        lon2d[:] = basin.getLongitude()
        tmp = MV2where(lon2d > lon1, 1, 0) + MV2where(lon2d < lon2, 1, 0)
        if latkey == "between":
            mask = MV2where(tmp != 2, 0, mask)
        else:
            mask = MV2where(tmp == 2, 0, mask)
    # apply mask
    tab_out = MV2masked_where(mask == 1, tab_in)
    tab_out = create_variable(
        tab_out,
        axes=tab_in.getAxisList(),
        grid=tab_in.getGrid(),
        mask=tab_in.mask,
        attributes=tab_in.attributes,
        id=tab_in.id
    )
    return tab_out, keyerror


def CheckTime(
        tab1, tab2,
        frequency="monthly",
        min_time_steps=None,
        metric_name="",
        debug=False,
        **kwargs
    ):
    """
    #################################################################################
    Description:
    Checks if tab1 and tab2 cover the same time period and adjust if not
    Checks if the minimum_length of the time period criterion if fulfilled
    #################################################################################

    :param tab1: masked_array
    :param tab2: masked_array
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :param min_time_steps: int, optional
        minimum number of time steps for the metric to make sens
        e.g., for 30 years of monthly data mintimesteps=360
    :param metric_name: string, optional
        name of the metric calling the function
    :return:
    """
    if debug is True:
        # dict_debug = {"shape1": "tab1.shape = " + str(tab1.shape), "shape2": "tab2.shape = " + str(tab2.shape),
        #               "time1": "tab1.time = " + str(TimeBounds(tab1)),
        #               "time2": "tab2.time = " + str(TimeBounds(tab2))}
        dict_debug = {"shape1": "tab1.shape = " + str(tab1.shape), "shape2": "tab2.shape = " + str(tab2.shape)}
        EnsoErrorsWarnings.debug_mode("\033[93m", "in CheckTime (input)", 20, **dict_debug)
    tab1 = _to_cdat(tab1)
    tab2 = _to_cdat(tab2)
    # gets dates of the first and last the time steps of tab1
    stime1 = _require_time_axis(tab1, 'CheckTime').asComponentTime()[0]
    etime1 = _require_time_axis(tab1, 'CheckTime').asComponentTime()[-1]

    # gets dates of the first and last the time steps of tab2
    stime2 = _require_time_axis(tab2, 'CheckTime').asComponentTime()[0]
    etime2 = _require_time_axis(tab2, 'CheckTime').asComponentTime()[-1]

    # retains only the latest start date and the earliest end date
    if stime1.year > stime2.year:
        stime = stime1
    elif stime1.year < stime2.year:
        stime = stime2
    else:
        if stime1.month > stime2.month:
            stime = stime1
        elif stime1.month < stime2.month:
            stime = stime2
        else:
            if stime1.day > stime2.day:
                stime = stime1
            elif stime1.day < stime2.day:
                stime = stime2
            else:
                stime = max(stime1, stime2)
    if etime1.year < etime2.year:
        etime = etime1
    elif etime1.year > etime2.year:
        etime = etime2
    else:
        if etime1.month < etime2.month:
            etime = etime1
        elif etime1.month > etime2.month:
            etime = etime2
        else:
            if etime1.day < etime2.day:
                etime = etime1
            elif etime1.day > etime2.day:
                etime = etime2
            else:
                etime = min(etime1, etime2)
    # stime = max(stime1, stime2)
    # etime = min(etime1, etime2)

    # defines the period between the two dates
    if frequency == "daily":
        stime_adjust = CDTIMEcomptime(stime.year, stime.month, stime.day, 0, 0, 0.0)
        etime_adjust = CDTIMEcomptime(etime.year, etime.month, etime.day, 23, 59, 0)
    elif frequency == "monthly":
        etime_day = monthrange(etime.year, etime.month)[-1]
        stime_adjust = CDTIMEcomptime(stime.year, stime.month, 1, 0, 0, 0.0)
        etime_adjust = CDTIMEcomptime(etime.year, etime.month, etime_day, 23, 59, 0)
    elif frequency == "yearly":
        stime_adjust = CDTIMEcomptime(stime.year, 1, 1, 0, 0, 0.0)
        etime_adjust = CDTIMEcomptime(etime.year, 12, 31, 23, 59, 0)
    else:
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())

    # retains only the time-period common to both tab1 and tab2
    tab1_sliced = tab1(time=(stime_adjust, etime_adjust))
    tab2_sliced = tab2(time=(stime_adjust, etime_adjust))
    if debug is True:
        # dict_debug = {"shape1": "tab1.shape = " + str(tab1_sliced.shape),
        #               "shape2": "tab2.shape = " + str(tab2_sliced.shape),
        #               "time1": "tab1.time = " + str(TimeBounds(tab1_sliced)),
        #               "time2": "tab1.time = " + str(tab1_sliced.getTime().asComponentTime()[:]),
        #               "time3": "tab2.time = " + str(TimeBounds(tab2_sliced)),
        #               "time4": "tab2.time = " + str(tab2_sliced.getTime().asComponentTime()[:])}
        dict_debug = {
            "shape1": "tab1.shape = " + str(tab1_sliced.shape),
            "shape2": "tab2.shape = " + str(tab2_sliced.shape)
        }
        EnsoErrorsWarnings.debug_mode("\033[93m", "in CheckTime (output)", 20, **dict_debug)
    if len(tab1_sliced.getTime()[:]) != len(tab2_sliced.getTime()[:]):
        keyerror1 = "missing time step within the given period"
    else:
        keyerror1 = None
        tab2_sliced.setAxis(0, tab1_sliced.getTime())

    # checks if the remaining time-period fulfills the minimum length criterion
    if min_time_steps is not None:
        if len(tab1_sliced) < min_time_steps or len(tab2_sliced) < min_time_steps:
            shortest = min(len(tab1_sliced), len(tab2_sliced))
            EnsoErrorsWarnings.too_short_time_period(metric_name, shortest, min_time_steps, INSPECTstack())
            keyerror2 = "too short time period (variable1:" + str(len(tab1_sliced)) + " ; variable2:" +\
                        str(len(tab2_sliced)) + ")"
        else:
            keyerror2 = None
    else:
        keyerror2 = None

    # errors
    if keyerror1 is not None or keyerror2 is not None:
        keyerror = add_up_errors([keyerror1, keyerror2])
    else:
        keyerror = None
    return tab1_sliced, tab2_sliced, keyerror


def CheckUnits(tab, var_name, name_in_file, units, return_tab_only=True, **kwargs):
    """
    #################################################################################
    Description:
    Checks the units of the variable and changes it if necessary
    Works for current/wind velocities, depth, heat flux, precipitation, pressure, temperature, wind stress

    Uses MV2 (uvcdat) to find the minimum value, to multiply and to subtract
    #################################################################################

    :param tab: array
        array containing 'var_name'
    :param var_name: string
        name of the variable included in 'tab'
    :param name_in_file: string
        name of the variable in the file (usually the short_name)
    :param units: string
        units of the variable included in 'tab'
    :param return_tab_only: boolean, optional
        default value = True, only the tab is returned
        True if you want only the tab, if you want the new units also pass anything but true
    :return tab: array
        array with new units (if applicable)
    """
    keyerror = None
    if var_name in ["temperature"]:
        if units in [
                "K", "Kelvin", "Kelvins", "degree K", "degree Kelvin", "degree Kelvins", "degree_K",
                "degree_Kelvin", "degree_Kelvins", "degreeK", "degreeKelvin", "degreeKelvins", "degrees K",
                "degrees Kelvin", "degrees Kelvins", "degrees_K", "degrees_Kelvin", "degrees_Kelvins", "degreesK",
                "degreesKelvin", "degreesKelvins", "deg K", "deg Kelvin", "deg Kelvins", "deg_K", "deg_Kelvin",
                "deg_Kelvins", "degK", "degKelvin", "degKelvins", "deg. K", "deg. Kelvin", "deg. Kelvins"
            ]:
            # check if the temperature units is really K
            if float(MV2minimum(tab)) > 150:
                # unit change of the temperature: from K to degC
                tab = dict_operations["minus"](tab, 273.15)
            else:
                minmax = [MV2minimum(tab), MV2maximum(tab)]
                EnsoErrorsWarnings.unlikely_units(var_name, name_in_file, units, minmax, INSPECTstack())
                keyerror = "unlikely units: " + str(units) + "(" + str(minmax) + ")"
        elif units in [
                "C", "celsius", "Celsius", "degree C", "degree celsius", "degree Celsius", "degree_C",
                "degree_celsius", "degree_Celsius", "degreeC", "degreecelsius", "degreeCelsius", "degrees C",
                "degrees celsius", "degrees Celsius", "degrees_C", "degrees_celsius", "degrees_Celsius",
                "degreesC", "degreescelsius", "degreesCelsius", "deg C", "deg celsius", "deg Celsius", "deg_C",
                "deg_celsius", "deg_Celsius", "degC", "degcelsius", "degCelsius", "deg. C", "deg. celsius",
                "deg. Celsius"
            ]:
            # check if the temperature units is really degC
            if float(MV2minimum(tab)) > 50:
                minmax = [MV2minimum(tab), MV2maximum(tab)]
                EnsoErrorsWarnings.unlikely_units(var_name, name_in_file, minmax, units, INSPECTstack())
                keyerror = "unlikely units: " + str(units) + "(" + str(minmax) + ")"
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "degC"
    elif var_name in ["precipitations"]:
        if units in [
                "kg/m2/s", "kg/m^2/s", "kg/m**2/s", "kg m-2 s-1", "kg m^-2 s^-1", "kg m**-2 s**-1", "Kg/m2/s",
                "Kg/m^2/s", "Kg/m**2/s", "Kg m-2 s-1", "Kg m^-2 s^-1", "Kg m**-2 s**-1"
            ]:
            # changes units of the precipitation flux: from kg/(m2.s) to mm/day
            # it must be divided by the density of water = 1000 kg/m3
            #     and multiplied by 1000 (m to mm) and by 60*60*24 (s to day)
            tab = dict_operations["multiply"](tab, 86400)
        elif units in ["mm/day", "mm day-1", "mm day^-1", "mm day**-1", "mm/d", "mm d-1", "mm d^-1", "mm d**-1"]:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "mm/day"
    elif var_name in ["wind stress"]:
        if units not in [
                "N/m2", "N/m^2", "N/m**2", "N m-2", "N m^-2", "N m**-2",
                "Pa", "pascal", "pascals", "Pascal", "Pascals"
            ]:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "N/m2"
    elif var_name in ["velocity"]:
        if units in ["cm/s", "cm s-1", "cm s^-1", "cm s**-1", "cm/sec", "cm sec-1", "cm sec^-1", "cm sec**-1"]:
            # unit change of the velocity: from cm/s to m/s
            tab = dict_operations["multiply"](tab, 1e-2)
        elif units in ["m/s", "m s-1", "m s^-1", "m s**-1", "m/sec", "m sec-1", "m sec^-1", "m sec**-1"]:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "m/s"
    elif var_name in ["heat flux"]:
        if units in [
                "W/m2", "W/m^2", "W/m**2", "W m-2", "W m^-2", "W m**-2", "Watt/m2", "Watt/m^2", "Watt/m**2",
                "Watt m-2", "Watt m^-2", "Watt m**-2", "Watts/m2", "Watts/m^2", "Watts/m**2", "Watts m-2",
                "Watts m^-2", "Watts m**-2"
            ]:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "W/m2"
    elif var_name in ["pressure"]:
        if units in [
                "N/m2", "N/m^2", "N/m**2", "N m-2", "N m^-2", "N m**-2",
                "Pa", "pascal", "pascals", "Pascal", "Pascals"
            ]:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "Pa"
    elif var_name in ["depth", "sea surface height"]:
        if units in ["cm", "centimeter", "centimeters"]:
            # unit change of the sea surface height: from cm to m
            tab = dict_operations["multiply"](tab, 1e-2)
        elif units in ["m", "meter", "meters"]:
            pass
        else:
            EnsoErrorsWarnings.unknown_units(var_name, name_in_file, units, INSPECTstack())
            keyerror = "unknown units: " + str(units) + "(as " + str(var_name) + ")"
        units = "m"
    else:
        list_strings = ["WARNING" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": variable name",
                        str().ljust(5) + "unknown variable name: " + var_name + " (" + name_in_file + ")"]
        EnsoErrorsWarnings.my_warning(list_strings)
    if return_tab_only is True:
        return tab
    else:
        return tab, units, keyerror


def Event_selection(tab, frequency, nbr_years_window=None, list_event_years=[]):
    tab = _to_cdat(tab)
    if frequency not in ["daily", "monthly", "yearly"]:
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())
    if len(list_event_years) == 0:
        tax = _require_time_axis(tab, 'Event_selection').asComponentTime()
        list_event_years = sorted(list(set([tax[ii].year for ii in list(range(len(tax)))])))
    else:
        list_event_years = sorted(list_event_years)
    # function to fill array with masked value where the data is not available
    def fill_array(tab, units, freq):
        # Robustly get the time axis: prefer getTime() but fall back to the
        # first axis when getTime() returns None (e.g. axis typed "-" instead
        # of "T" after a slice that didn't preserve axis metadata).
        def _get_tax(v):
            ax = v.getTime() if isinstance(v, CDATVariable) else None
            if ax is None and isinstance(v, CDATVariable):
                axes = v.getAxisList()
                ax = axes[0] if axes else None
            if ax is None:
                raise ValueError(
                    f"No time axis in fill_array input "
                    f"(id={getattr(v, 'id', '?')!r}, shape={getattr(v, 'shape', None)})"
                )
            return ax
        _tab_comp = _get_tax(tab).asComponentTime()
        if not _tab_comp:
            raise ValueError(
                f"fill_array: input time axis is empty "
                f"(id={getattr(tab, 'id', '?')!r}, shape={getattr(tab, 'shape', None)})"
            )
        y1, m1, d1 = _tab_comp[0].year, _tab_comp[0].month, _tab_comp[0].day
        if len(tab.shape) == 1:
            raw_out = MV2zeros(nbr_years_window * 12)
        elif len(tab.shape) == 2:
            raw_out = MV2zeros((nbr_years_window * 12, tab.shape[1]))
        else:
            raw_out = MV2zeros((nbr_years_window * 12, tab.shape[1], tab.shape[2]))
        raw_out = MV2masked_where(raw_out == 0, raw_out)
        time_axis = create_axis(list(range(len(raw_out))), id="time", units=units, axis_type="T")
        other_axes = tab.getAxisList()[1:] if isinstance(tab, CDATVariable) and len(tab.shape) > 1 else []
        tab_out = create_variable(
            raw_out, axes=[time_axis] + other_axes,
            grid=tab.getGrid() if isinstance(tab, CDATVariable) else None,
            id=getattr(tab, 'id', '')
        )
        _tab_out_comp = _get_tax(tab_out).asComponentTime()
        for ii in list(range(len(tab))):
            y2 = _tab_out_comp[ii].year
            m2 = _tab_out_comp[ii].month
            d2 = _tab_out_comp[ii].day
            if freq == "yearly":
                if y2 == y1:
                    tab_out[ii:ii + len(tab)] = copy.copy(tab)
                    break
            elif freq == "monthly":
                if y2 == y1 and m2 == m1:
                    tab_out[ii:ii + len(tab)] = copy.copy(tab)
                    break
            elif freq == "daily":
                if y2 == y1 and m2 == m1 and d2 == d1:
                    tab_out[ii:ii + len(tab)] = copy.copy(tab)
                    break
        return tab_out
    # compute composite
    if nbr_years_window is not None:
        composite = list()
        for yy in list_event_years:
            # first and last years of the window
            yy1, yy2 = yy + 1 - nbr_years_window // 2, yy + nbr_years_window // 2
            # create time bounds from "first and last years of the window"
            timebnds = (str(yy1) + "-01-01 00:00:00.0", str(yy2) + "-12-31 23:59:59.0")
            # select the right time period in the given tab
            tmp1 = tab(time=timebnds)
            # sometimes there is some errors with "time=timebnds"
            # if the time slice selected has the right length: do nothing
            # else: fill the beginning / end of the time series by masked values (done by the function "fill_array")
            if frequency == "yearly":
                length = nbr_years_window
                units = "years since " + timebnds[0]
                units_out = "years since 0001-07-02 12:00:00"
            elif frequency == "monthly":
                length = nbr_years_window * 12
                units = "months since " + timebnds[0]
                units_out = "months since 0001-01-15 12:00:00"
            elif frequency == "daily":
                date1 = date(yy1, 1, 1)
                date2 = date(yy2, 12, 31)
                length = (date2 - date1).days
                units = "days since " + timebnds[0]
                units_out = "days since 0001-01-01 12:00:00"
            if len(tmp1) == length:
                tmp2 = copy.copy(tmp1)
            else:
                tmp2 = fill_array(tmp1, units, frequency)
            # save the selected time slice
            composite.append(tmp2)
        composite = MV2array(composite)
        # axis list
        axis0 = create_axis(MV2array(list_event_years, dtype="int32"), id="years")
        axis1 = create_axis(list(range(len(composite[0]))), id="months")
        axis1.units = units_out
        axes = [axis0, axis1]
        if len(tab.shape) > 1:
            axes = axes + tab.getAxisList()[1:]
        composite.setAxisList(axes)
    else:
        time_ax = _require_time_axis(tab, 'Event_selection').asComponentTime()  # gets component time of tab
        list_years = [yy.year for yy in time_ax[:]]  # listing years in tab (from component time)
        indices = MV2arange(tab.size)
        # creates a tab of "condition" where True is set when the event is found, False otherwise
        try:
            condition = [True if yy in list_event_years else False for yy in list_years]
        except Exception:
            list_event_years = [str(yy) for yy in list_event_years]
            condition = [True if str(yy) in list_event_years else False for yy in list_years]
        ids = MV2compress(condition, indices)  # gets indices of events
        composite = MV2take(tab, ids, axis=0)  # gets events
        axis0 = create_axis(MV2array(list_event_years, dtype="int32"), id="years")
        composite.setAxis(0, axis0)
    return composite


def Composite(tab, list_event_years, frequency, nbr_years_window=None):
    return MV2average(
        Event_selection(tab, frequency, nbr_years_window=nbr_years_window, list_event_years=list_event_years), axis=0)


def DetectEvents(tab, season, threshold, normalization=False, nino=True, compute_season=True, duration=1):
    """
    #################################################################################
    Description:
    Detects Nina or Nino events
    These events are detected when 'tab' anomalies during 'season' are above (less) then 'threshold'
    The anomalies can be normalized
    Uses MV2 (uvcdat) to create an empty array, to create an array of indices, to define conditions, to select the
    indices depending on the conditions and to select the years of the events
    #################################################################################
    :param tab: masked_array
        masked_array containing a variable from which the events are detected. Most likely SST
    :param season: string
        one month (e.g, 'DEC'), two months (e.g., 'DJ'), three months (e.g., 'NDJ'), four months (e.g., 'NDJF'), period
        when the events are detected
    :param threshold: float
        threshold to define the events (e.g., 0.75 for El Nino, -0.75 for La Nina)
    :param normalization: boolean, optional
        True if events are detected based on the standard deviation, if not pass anything but True
    :param nino: boolean, optional
        True if events are detected if above threshold (El Nino like), if not pass anything but True (La Nina like)
    :param compute_season: boolean, optional
        True if the the seasonal mean anomalies during the given season needs to be computed
    :param duration: integer, optional
        number of consecutive months/seasons when given threshold must be met to define an ENSO event
    :return list_of_years: list
        list of years including a detected event
    """
    if duration == 1:
        # Seasonal mean and anomalies
        if compute_season is True:
            tab = SeasonalMean(tab, season, compute_anom=True)
        tab = _to_cdat(tab)
        # Normalization ?
        if normalization is True:
            threshold = threshold * float(GENUTILstd(tab, axis=0, centered=1, biased=1))
        # Initialization
        tab_threshold = MV2zeros(tab.shape)
        tab_threshold.fill(threshold)
        list_years = sorted(list(set([_require_time_axis(tab, 'DetectEvents').asComponentTime()[yy].year for yy in range(len(tab))])))
        indices = MV2arange(len(list_years))
        # Conditions
        if nino is True:
            condition = MV2where(tab > tab_threshold, True, False)
        else:
            condition = MV2where(tab < tab_threshold, True, False)
        # Indices of the events
        ids = MV2compress(condition, indices)
        # Events years
        events = list(MV2take(list_years, ids, axis=0))
    else:
        if season == "DEC":
            lseasons = ["NOV", season, "JAN"]
            if duration >= 3:
                lseasons = ["OCT"] + lseasons + ["FEB"]
            if duration >= 4:
                lseasons = ["SEP"] + lseasons + ["MAR"]
            if duration >= 5:
                lseasons = ["AUG"] + lseasons + ["APR"]
            if duration >= 6:
                lseasons = ["JUL"] + lseasons + ["MAY"]
        elif season == "NDJ":
            lseasons = ["OND", season, "DJF"]
            if duration >= 3:
                lseasons = ["SON"] + lseasons + ["JFM"]
            if duration >= 4:
                lseasons = ["ASO"] + lseasons + ["FMA"]
            if duration >= 5:
                lseasons = ["JAS"] + lseasons + ["MAM"]
            if duration >= 6:
                lseasons = ["JJA"] + lseasons + ["AMJ"]
        else:
            lseasons = list()
            list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": season",
                            str().ljust(5) + "unknown season for ENSO event detection: " + str(season)]
            EnsoErrorsWarnings.my_error(list_strings)
        # Main seasonal mean and anomalies
        enso = SeasonalMean(tab, season, compute_anom=True)
        list_years = [_require_time_axis(enso, 'DetectEvents').asComponentTime()[yy].year for yy in range(len(enso))]
        indices = MV2arange(len(list_years))
        y0 = list_years[0]
        enso_by_sea = list()
        for sea in lseasons:
            # Seasonal mean and anomalies
            tmp = SeasonalMean(tab, sea, compute_anom=True)
            y1 = _require_time_axis(tmp, 'DetectEvents').asComponentTime()[0].year
            if y1 == y0:
                tmp = tmp[:len(enso)]
            elif y1+1 == y0:
                tmp = tmp[1:len(enso)+1]
            if sea == "DEC" or (season == "NDJ" and sea == "DJF"):
                y0 += 1
            enso_by_sea.append(tmp)
            del tmp, y1
        enso_by_sea = MV2array(enso_by_sea)
        # Normalization ?
        if normalization is True:
            thr = threshold * float(GENUTILstd(enso, axis=0, centered=1, biased=1))
        else:
            thr = copy.deepcopy(threshold)
        # Initialization
        tab_threshold = MV2zeros(enso_by_sea.shape)
        tab_threshold.fill(thr)
        # Conditions
        if nino is True:
            condition = MV2where(enso_by_sea > tab_threshold, 1, 0)
        else:
            condition = MV2where(enso_by_sea < tab_threshold, 1, 0)
        # sum by duration window
        d_window = MV2zeros(enso_by_sea.shape)
        d_window = d_window[:len(lseasons) - duration + 1]
        for ii in range(len(d_window)):
            d_window[ii] = MV2sum(condition[ii: ii + duration], axis=0)
        # test if threshold met during at least duration
        tab_threshold = MV2zeros(d_window.shape)
        tab_threshold.fill(duration)
        condition = MV2where(d_window >= tab_threshold, 1, 0)
        condition = MV2sum(condition, axis=0)
        condition = MV2where(condition >= 1, True, False)
        # Indices of the events
        ids = MV2compress(condition, indices)
        # Events years
        events = list(MV2take(list_years, ids, axis=0))
    return events


def Detrend(tab, info, axis=0, method="linear", bp=0):
    """
    #################################################################################
    Description:
    Removes trend along 'axis' from 'tab'
    #################################################################################

    :param tab: array
        tab of data to detrend
    :param info: string
        information about what is done to 'tab'
    :param axis: int, optional
        axis along which to detrend the data
        default value is the first axis (0)
    :param method: string, optional
        detrending method:
        'constant': only the mean of 'tab' is subtracted
        'linear':   the result of a linear least-squares fit to 'tab' is subtracted from 'tab'
    :param bp: array of integer, optional
        a sequence of break points. If given, an individual linear fit is performed for each part of 'tab' between two
        break points
        break points are specified as indices into 'tab'
    :return new_tab: array
        detrended data
    """
    if method not in ["linear", "constant"]:
        new_tab = None
        keyerror = "cannot detrend: unknown method"
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": method",
                        str().ljust(5) + "unknown method: " + str(method)]
        EnsoErrorsWarnings.my_warning(list_strings)
    else:
        tab = _to_cdat(tab)
        axes = tab.getAxisList()
        grid = tab.getGrid()
        mask = tab.mask
        mean, keyerror = AverageTemporal(tab)
        # Fill masked values before passing to scipy.signal.detrend so the
        # computation does not receive NaN/inf (masked positions are restored
        # afterwards via MV2masked_where).
        _raw_for_detrend = _mv(tab).filled(0.0) if ma.isMaskedArray(_mv(tab)) else np.asarray(tab)
        new_tab = MV2array(SCIPYsignal_detrend(_raw_for_detrend, axis=axis, type=method, bp=bp))
        new_tab = new_tab + mean
        new_tab = MV2masked_where(mask, new_tab)
        new_tab.setAxisList(axes)
        new_tab.setGrid(grid)
        if method == "linear":
            info = info + ", time series are linearly detrended"
        else:
            info = info + ", the mean value of the time series is subtracted"
    return new_tab, info, keyerror


def DurationAllEvent(tab, threshold, nino=True, debug=False):
    """
    #################################################################################
    Description:
    Duration of Nina or Nino events
    The duration is the number of consecutive timestep when tab < threshold for La Nina and tab > threshold for El Nino

    Uses xarray/XarrayCompat for axes
    #################################################################################

    :param tab: masked_array
        masked_array containing a variable from which the events are detected. Most likely SST
    :param threshold: float
        threshold to define the events (e.g., 0.75 for El Nino, -0.75 for La Nina)
    :param nino: boolean, optional
        True if events are detected if above threshold (El Nino like), if not pass anything but True (La Nina like)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :return list_of_years: list
        list of years including a detected event
    """
    tmp = MV2array([DurationEvent(tab[tt], threshold, nino=nino, debug=debug) for tt in list(range(len(tab)))])
    tmp.setAxis(0, tab.getAxis(0))
    return tmp


def DurationEvent(tab, threshold, nino=True, debug=False):
    """
    #################################################################################
    Description:
    Duration of Nina or Nino events
    The duration is the number of consecutive timestep when tab < threshold for La Nina and tab > threshold for El Nino

    Uses MV2 (uvcdat)
    #################################################################################

    :param tab: masked_array
        masked_array containing a variable from which the events are detected. Most likely SST
    :param threshold: float
        threshold to define the events (e.g., 0.75 for El Nino, -0.75 for La Nina)
    :param nino: boolean, optional
        True if events are detected if above threshold (El Nino like), if not pass anything but True (La Nina like)
    :param debug: bolean, optional
        default value = False debug mode not activated
        If you want to activate the debug mode set it to True (prints regularly to see the progress of the calculation)
    :return list_of_years: list
        list of years including a detected event
    """
    mask_val = tab.mask
    mask = mask_val if hasattr(mask_val, '__len__') and not isinstance(mask_val, np.bool_) else [mask_val]
    # if debug is True:
    #     dict_debug = {'line1': 'threshold = ' + str(threshold) + '  ;  nino = ' + str(nino)
    #                            + '  ;  len(tab) = ' + str(len(tab)),
    #                   'line2': 'tab = ' + str(tab) + '\nmask = ' + str(mask)
    #                   }
    #     EnsoErrorsWarnings.DebugMode('\033[93m', 'in DurationEvent', 20, **dict_debug)
    if all(ii is False for ii in mask) is True:
        pass
    else:
        if nino is True:
            tab = MV2where(tab.mask, -9999, tab)
        else:
            tab = MV2where(tab.mask, 9999, tab)
    # if debug is True:
    #     dict_debug = {'line1': 'after unmasking',
    #                   'line2': 'tab = ' + str(tab)}
    #     EnsoErrorsWarnings.DebugMode('\033[93m', 'in DurationEvent', 20, **dict_debug)
    tmp1 = list(reversed(tab[: len(tab) // 2]))
    tmp2 = list(tab[len(tab) // 2:])
    if nino is True:
        try:
            nbr_before = next(x[0] for x in enumerate(tmp1) if x[1] <= threshold)
        except Exception:
            if all(ii == -9999 for ii in tmp1):
                nbr_before = 0
            elif all(ii > threshold for ii in tmp1):
                nbr_before = len(tmp1)
        try:
            nbr_after = next(x[0] for x in enumerate(tmp2) if x[1] <= threshold)
        except Exception:
            if all(ii == -9999 for ii in tmp2):
                nbr_after = 0
            elif all(ii > threshold for ii in tmp2):
                nbr_after = len(tmp2)
    else:
        try:
            nbr_before = next(x[0] for x in enumerate(tmp1) if x[1] >= threshold)
        except Exception:
            if all(ii == 9999 for ii in tmp1):
                nbr_before = 0
            elif all(ii < threshold for ii in tmp1):
                nbr_before = len(tmp1)
        try:
            nbr_after = next(x[0] for x in enumerate(tmp2) if x[1] >= threshold)
        except Exception:
            if all(ii == 9999 for ii in tmp2):
                nbr_after = 0
            elif all(ii < threshold for ii in tmp2):
                nbr_after = len(tmp2)
    duration = nbr_before + nbr_after
    # if debug is True:
    #     dict_debug = {'line1': 'duration of the event = ' + str(duration)}
    #     EnsoErrorsWarnings.DebugMode('\033[93m', 'in DurationEvent', 20, **dict_debug)
    return duration


def get_num_axis(tab, name_axis):
    """
    #################################################################################
    Description:
    Finds the number of the axis named "name_axis"
    #################################################################################

    :param tab: array
        tab of data to normalize by the standard deviation
    :param name_axis: string
        name of an axis
        e.g., name_axis='latitude'
    :return number: int
        position of the axis named "name_axis"
    """
    num = None
    if name_axis == "depth":
        axis_nick = "lev"
        axis_nicks = ["z", "Z", "st_ocean", "sw_ocean"]
    if name_axis == "latitude":
        axis_nick = "lat"
        axis_nicks = ["j", "y", "Y", "yt_ocean", "yu_ocean"]
    elif name_axis == "longitude":
        axis_nick = "lon"
        axis_nicks = ["i", "x", "X", "xt_ocean", "xu_ocean"]
    elif name_axis == "time":
        axis_nick = "time"
        axis_nicks = ["t", "T"]
    # Fast path: use CF axis type metadata when available — also avoids
    # IndexError when tab._axes is shorter than tab.shape (empty axes list).
    # Only activated when the target axis type is actually present in _axes;
    # pure-heuristic lookups can return wrong indices for 1-D time series.
    _ax_type_map = {"latitude": "y", "longitude": "x", "time": "t", "depth": "z"}
    _ax_cf_map   = {"latitude": "Y", "longitude": "X", "time": "T", "depth": "Z"}
    if name_axis in _ax_type_map and isinstance(tab, CDATVariable):
        if any(ax is not None and ax.axis == _ax_cf_map[name_axis] for ax in tab._axes):
            idx = _axis_to_int(tab, _ax_type_map[name_axis])
            if isinstance(idx, (int, np.integer)):
                return int(idx)
    # String-matching fallback — guard against _axes being shorter than shape.
    tab = _to_cdat(tab)
    axlist = tab.getAxisList()
    for nn in list(range(len(tab.shape))):
        if nn >= len(axlist):
            break
        if axis_nick in axlist[nn].id:
            num = nn
            break
    if num is None:
        for nn in list(range(len(tab.shape))):
            if nn >= len(axlist):
                break
            for ax in axis_nicks:
                if ax == axlist[nn].id:
                    num = nn
                    break
    if num is None:
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": axis",
                        str().ljust(5) + "cannot find axis named: " + str(name_axis),
                        str().ljust(5) + "axes: " + str(tab.getAxisList())]
        EnsoErrorsWarnings.my_error(list_strings)
    return num


def get_year_by_year(tab, frequency="monthly"):
    """
    #################################################################################
    Description:
    Reshape array to a year by year array
    #################################################################################

    :param tab: masked_array
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :return: tab: array
        array of the year by year values
    """
    tab = _to_cdat(tab)
    tab = tab.reorder("t...")
    time_ax = _require_time_axis(tab, 'Reshape').asComponentTime()
    myshape = [1] + [ss for ss in tab.shape[1:]]
    zeros = MV2zeros(myshape)
    zeros = MV2masked_where(zeros == 0, zeros)
    if frequency == "daily":
        days = MV2array(list(tt.day for tt in time_ax))
        months = MV2array(list(tt.month for tt in time_ax))
        months = MV2array([(mm*100)+dd for dd, mm in zip(days, months)])
        tmm = create_axis(list(range(365)), id="days")
        m1 = time_ax[0].day
        m2 = time_ax[-1].day
        t2 = 365
    elif frequency == "monthly":
        months = MV2array(list(tt.month for tt in time_ax))
        tmm = create_axis(list(range(12)), id="months")
        m1 = time_ax[0].month
        m2 = time_ax[-1].month
        t2 = 12
    else:
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())
    years = sorted(set(MV2array(list(tt.year for tt in time_ax))))
    tyy = create_axis(MV2array(years, dtype="int32"), id="years")
    axes = [tyy] + [tmm]
    val = sorted(set(months))
    tab_out = list()
    if frequency == "daily":
        val.remove(129)
    for ii in val:
        tmp = tab.compress(months == ii, axis=0)
        if m1 != 1 and len(tmp) != len(years):
            tmp = MV2concatenate((zeros, tmp))
        if m2 != t2 and len(tmp) != len(years):
            tmp = MV2concatenate((tmp, zeros))
        tab_out.append(tmp)
    tab_out = MV2array(tab_out)
    tab_out = MV2masked_where(tab_out == 0, tab_out)
    tab_out = tab_out.reorder("10")
    if len(tab.shape) == 1:
        tab_out = create_variable(tab_out, axes=axes, attributes=tab.attributes, id=tab.id)
    else:
        axes = axes + tab.getAxisList()[1:]
        grid = tab[0].getGrid()
        mask = tab[0].mask
        mask_out = MV2zeros(tab_out.shape)
        mask_out[:, :] = mask
        tab_out = create_variable(
            tab_out,
            axes=axes,
            grid=grid,
            mask=mask_out,
            attributes=tab.attributes,
            id=tab.id
        )
    return tab_out


def MinMax(tab):
    return [float(MV2minimum(tab)), float(MV2maximum(tab))]


def MyEmpty(tab, time=True, time_id=''):
    tab = _to_cdat(tab)
    tab_out = ArrayZeros(tab)
    if time is True:
        axis = create_axis(MV2array(len(tab_out), dtype="int32"), id=time_id)
        axes = [axis] + tab.getAxisList()[1:]
    else:
        axes = tab.getAxisList()
    tab_out.setAxisList(axes)
    return tab_out


def Normalize(tab, frequency):
    """
    #################################################################################
    Description:
    Removes trend along 'axis' from 'tab'
    #################################################################################

    :param tab: array
        tab of data to normalize by the standard deviation
    :param frequency: string, optional
        time frequency of the datasets
        e.g., frequency='monthly'
    :return tab: masked_array
        normalized data
    """
    tab = _to_cdat(tab)
    keyerror = None
    axes = tab.getAxisList()
    if frequency == "daily":
        time_steps_per_year = 365
    elif frequency == "monthly":
        time_steps_per_year = 12
    elif frequency == "yearly":
        time_steps_per_year = 1
    else:
        keyerror = "unknown frequency"
        time_steps_per_year = None
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())
    if time_steps_per_year is not None:
        if len(tab) % time_steps_per_year != 0:
            tab_out = None
            keyerror = "cannot perform normalization: the function can only handle full years (len(tab) = " +\
                str(len(tab)) + ")"
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": data length",
                str().ljust(5) + "the normalization function can only handle full years: " +
                str(len(tab) // time_steps_per_year) + " years " + str(len(tab) % time_steps_per_year),
                str().ljust(10) + "frequency: " + str(frequency) + " (time steps per year = " +
                str(time_steps_per_year) + "), len(dataset) = " + str(len(tab)) + ", so " +
                str(len(tab) / float(time_steps_per_year)) + " years",
            ]
            EnsoErrorsWarnings.my_warning(list_strings)
        else:
            # reshape tab like [yy,nb]
            new_tab = list()
            for yy in list(range(len(tab) // time_steps_per_year)):
                new_tab.append(tab[yy * time_steps_per_year:(yy + 1) * time_steps_per_year])
            new_tab = MV2array(new_tab)
            std = MV2zeros(new_tab[0].shape)
            for dd in list(range(time_steps_per_year)):
                std[dd] = float(GENUTILstd(new_tab[:,dd], weights=None, axis=0, centered=1, biased=1))
            tab_out = copy.copy(tab)
            for yy in list(range(len(tab) // time_steps_per_year)):
                tab_out[yy * time_steps_per_year:(yy + 1) * time_steps_per_year] = \
                    tab_out[yy * time_steps_per_year:(yy + 1) * time_steps_per_year] / std
            if len(tab.shape) == 1:
                tab_out = create_variable(tab_out, axes=axes, attributes=tab.attributes, id=tab.id)
            else:
                grid = tab.getGrid()
                mask = tab.mask
                tab_out = create_variable(tab_out, axes=axes, grid=grid, mask=mask, attributes=tab.attributes, id=tab.id)
    else:
        tab_out = None  # unknown frequency; keyerror already set above
    return tab_out, keyerror


def ReadAndSelectRegion(filename, varname, box=None, time_bounds=None, frequency=None, **kwargs):
    """
    #################################################################################
    Description:
    Read a variable from one or more input files and optionally select a time
    period and geographic region.

    This function preserves the legacy ENSO_metrics interface while using the
    modern xarray/xcdat-compatible file-reading pathway and CDAT-like
    compatibility objects provided by the refactored workflow.
    #################################################################################

    :param filename: string or list
        Path or paths to the input NetCDF file(s).

    :param varname: string
        Name of the variable to read from ``filename``.

    :param box: string, optional
        Region name to select. The region must be defined in
        ``EnsoCollectionsLib.ReferenceRegions``.
        default value = None

    :param time_bounds: tuple, optional
        First and last dates to extract from the file(s), for example
        ``('1979-01-01T00:00:00', '2017-01-01T00:00:00')``.
        default value = None

    :param frequency: string, optional
        Time frequency of the dataset, for example ``"monthly"``.
        default value = None

    usual kwargs:
        Additional options passed through the existing ENSO_metrics reading,
        selection, and preprocessing workflow.

    :return tab: CDATVariable or masked-array-like
        Variable read from ``filename`` and, if requested, selected over
        ``box`` and ``time_bounds``.
    """
    fi = open_file(filename)
    if box is None:  # no box given
        if time_bounds is None:  # no time period given
            # read file
            tab = fi(varname)
        else:  # time period given by the user
            # read file
            tab = fi(varname, time=time_bounds)
    else:  # box given by the user
        # define box
        region_ref = ReferenceRegions(box)
        if time_bounds is None:  # no time period given
            #  read file
            tab = fi(varname, latitude=region_ref["latitude"], longitude=region_ref["longitude"])
        else:
            # read file
            tab = fi(varname, time=time_bounds, latitude=region_ref["latitude"], longitude=region_ref["longitude"])
    # sign correction
    try:
        att1 = tab.attributes["standard_name"].lower().replace(" ", "_")
    except Exception:
        att1 = ''
    try:
        att2 = tab.attributes["long_name"].lower().replace(" ", "_")
    except Exception:
        att2 = ''
    reversed_sign = False
    if "latent_heat" in att1 or "latent_heat" in att2 or "sensible_heat" in att1 or "sensible_heat" in att2 or\
            (varname in ["tauu", "tauuo", "tauv", "tauvo", "taux", "tauy", "uflx", "vflx"]):
        if "upward" in att1 or "upward" in att2 or\
                (varname in ["tauu", "tauuo", "tauv", "tauvo", "taux", "tauy", "uflx", "vflx"] and
                ("in_air" in att1 or "in_air" in att2)):
            # I need to be in the ocean point of view so the heat fluxes must be downwards
            print("\033[93m" + str().ljust(15) + "EnsoUvcdatToolsLib ReadAndSelectRegion" + "\033[0m")
            print("\033[93m" + str().ljust(25) + varname + " sign reversed" + "\033[0m")
            print(
                "\033[93m" + str().ljust(5) + "range old = " + "{0:+.2f}".format(round(MV2minimum(tab), 2)) + " to " +
                "{0:+.2f}".format(round(MV2maximum(tab), 2)) + "\033[0m"
            )
            tab = -1 * tab
            print(
                "\033[93m" + str().ljust(5) + "range new = " + "{0:+.2f}".format(round(MV2minimum(tab), 2)) + " to " +
                "{0:+.2f}".format(round(MV2maximum(tab), 2)) + "\033[0m"
            )
            reversed_sign = True
    # CDATVariable arithmetic (-1 * tab) preserves axes via __rmul__/_wrap_binary.
    # _to_cdat is a no-op here but kept as a safety net for any branch that
    # didn't go through _finalize_cdat (e.g. a future fallback read path).
    tab = _to_cdat(tab)
    if time_bounds is not None:
        # sometimes the time boundaries are wrong, even with 'time=time_bounds'
        # this section checks if one time step has not been included by error at the beginning or the end of the time
        # series
        if isinstance(time_bounds[0], str):
            _comp = _require_time_axis(tab, 'ReadAndSelectRegion').asComponentTime()
            if _comp and str(_comp[0]) < time_bounds[0]:
                tab = tab[1:]
            _comp = _require_time_axis(tab, 'ReadAndSelectRegion').asComponentTime()
            if _comp and str(_comp[-1]) > time_bounds[1]:
                tab = tab[:-1]
    time_ax = _require_time_axis(tab, 'ReadAndSelectRegion')
    # Find the actual position of the time axis rather than assuming index 0.
    # CMIP6 files are always (time, lat, lon) but defensive coding ensures
    # non-standard dim orders do not corrupt the axis metadata.
    _t_ax_idx = next(
        (i for i, ax in enumerate(tab.getAxisList()) if ax is not None and ax.axis == "T"),
        0,  # safe fallback: CMIP6/ERA5 files always have T at 0
    )
    _comp0 = time_ax.asComponentTime()
    if not _comp0:
        raise ValueError(
            f"ReadAndSelectRegion: empty time axis after reading {varname!r} from {filename!r}. "
            "Check that time_bounds overlaps the data range or that the file has valid time data."
        )
    time_units = "days since " + str(_comp0[0].year) + "-01-01 12:00:00"
    time_ax.id = "time"
    time_ax.toRelativeTime(time_units)
    tab.setAxis(_t_ax_idx, time_ax)
    if frequency is None:  # no frequency given
        pass
    elif frequency == "daily":
        cdutil.setTimeBoundsDaily(tab)
    elif frequency == "monthly":
        cdutil.setTimeBoundsMonthly(tab)
    elif frequency == "yearly":
        cdutil.setTimeBoundsYearly(tab)
    else:
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())
    # remove axis 'level' if its length is 1
    if tab.getLevel():
        if len(tab.getLevel()) == 1:
            tab = tab(squeeze=1)
    # HadISST has -1000 values... mask them
    if "HadISST" in filename or "hadisst" in filename:
        tab = MV2masked_where(tab == -1000, tab)
    # Force the spatial mask to be constant through time (original CDAT invariant):
    # any grid point that is masked at ANY time step is masked for ALL time steps.
    # Using pure numpy operations here avoids fragile CDATVariable intermediate
    # steps on the mask array (ma.getmaskarray always returns a full-shape bool
    # array, never a scalar False, so edge cases are handled cleanly).
    _raw_mask = ma.getmaskarray(tab._data)  # full-shape bool array, dtype bool
    # Use the actual T-axis index so non-standard dim orders (unlikely but
    # possible for observational files) are handled correctly.
    _t_mask_idx = next(
        (i for i, ax in enumerate(tab.getAxisList()) if ax is not None and ax.axis == "T"),
        0,
    )
    if _raw_mask.ndim >= 2 and _raw_mask.shape[_t_mask_idx] > 1:
        _spatial_mask = np.any(_raw_mask, axis=_t_mask_idx)   # True = masked at >= 1 t
        # np.any collapses the T dimension; take a T=0 slice for comparison
        _t0_mask = np.take(_raw_mask, 0, axis=_t_mask_idx)
        if np.any(_spatial_mask != _t0_mask):   # mask not yet constant through time
            _full_mask = np.broadcast_to(
                np.expand_dims(_spatial_mask, axis=_t_mask_idx),
                _raw_mask.shape,
            ).copy()
            tab = MV2masked_where(_full_mask, tab)
    # check taux sign
    if varname in ["taux", "tauu", "tauuo", "uflx"] and reversed_sign is False:
        # define box
        region_ref = ReferenceRegions("nino4")
        if time_bounds is None:  # no time period given
            #  read file
            taux = fi(varname, latitude=region_ref["latitude"], longitude=region_ref["longitude"])
        else:
            # read file
            taux = fi(varname, time=time_bounds, latitude=region_ref["latitude"], longitude=region_ref["longitude"])
        # horizontal average
        taux, keyerror = AverageHorizontal(taux, region="nino4")
        if keyerror is None:
            taux, keyerror = AverageTemporal(taux)
            if keyerror is None and float(taux) > 0:
                print(
                    "\033[93m" + str().ljust(25) + "NOTE: taux sign reversed by the code (mean nino4 = " +
                    str(float(taux)) + ")" + "\033[0m"
                )
                tab = -1 * tab
    fi.close()
    # Re-run finalization after all in-place mutations (sign flip, slicing,
    # setAxis, toRelativeTime, squeeze, masking) to guarantee the returned
    # CDATVariable has valid, CF-typed axes.
    tab = _finalize_existing_cdat(tab, context=f"ReadAndSelectRegion:{varname}", require_time=True)
    return tab


def ReadAreaSelectRegion(filename, areaname='', box=None, **kwargs):
    """
    #################################################################################
    Description:
    Read an area-cell field from an input file and optionally select a
    geographic region.

    This function preserves the legacy ENSO_metrics interface while using the
    modern xarray/xcdat-compatible file-reading pathway and CDAT-like
    compatibility objects provided by the refactored workflow.
    #################################################################################

    :param filename: string
        Path to the input NetCDF file.

    :param areaname: string, optional
        Name of the area-cell variable in ``filename``, for example
        ``"areacella"`` or ``"areacello"``.
        default value = ``''``

    :param box: string, optional
        Region name to select. The region must be defined in
        ``EnsoCollectionsLib.ReferenceRegions``.
        default value = None

    usual kwargs:
        Additional options passed through the existing ENSO_metrics reading
        and regional-selection workflow.

    :return area: CDATVariable or masked-array-like
        Area-cell field read from ``filename`` and, if requested, selected over
        ``box``.
    """
    fi = open_file(filename)
    if box is None:  # no box given
        # read file
        try:
            areacell = fi(areaname)
        except Exception:
            try:
                areacell = fi('areacell')
            except Exception:
                try:
                    areacell = fi('areacella')
                except Exception:
                    try:
                        areacell = fi('areacello')
                    except Exception:
                        areacell = None
    else:  # box given by the user
        # define box
        region_ref = ReferenceRegions(box)
        # read file
        try:
            areacell = fi(areaname, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
        except Exception:
            try:
                areacell = fi('areacell', latitude=region_ref['latitude'], longitude=region_ref['longitude'])
            except Exception:
                try:
                    areacell = fi('areacella', latitude=region_ref['latitude'], longitude=region_ref['longitude'])
                except Exception:
                    try:
                        areacell = fi('areacello', latitude=region_ref['latitude'], longitude=region_ref['longitude'])
                    except Exception:
                        areacell = None
    fi.close()
    # Ensure areacell has valid spatial axes after read (no time axis expected).
    if areacell is not None:
        areacell = _finalize_existing_cdat(areacell, context=f"ReadAreaSelectRegion:{areaname}")
    return areacell

def ReadLandmaskSelectRegion(tab, filename, landmaskname='', box=None, **kwargs):
    """
    #################################################################################
    Description:
    Read a landmask field from an input file and optionally select a geographic
    region.

    This function preserves the legacy ENSO_metrics interface while using the
    modern xarray/xcdat-compatible file-reading pathway and CDAT-like
    compatibility objects provided by the refactored workflow. If no native
    landmask is provided or available, the workflow may fall back to an
    estimated landmask based on the input field.
    #################################################################################

    :param tab: CDATVariable or masked-array-like
        Input field used as a reference for grid, region, or fallback landmask
        generation.

    :param filename: string
        Path to the input NetCDF file.

    :param landmaskname: string, optional
        Name of the landmask variable in ``filename``, for example ``"sftlf"``,
        ``"lsmask"``, or ``"landmask"``.
        default value = ``''``

    :param box: string, optional
        Region name to select. The region must be defined in
        ``EnsoCollectionsLib.ReferenceRegions``.
        default value = None

    usual kwargs:
        Additional options passed through the existing ENSO_metrics reading,
        regional-selection, and fallback landmask workflow.

    :return landmask: CDATVariable or masked-array-like
        Landmask field read from ``filename`` and, if requested, selected over
        ``box``.
    """
    # Get landmask
    if OSpath__isfile(filename):
        # Open file
        fi = open_file(filename)

        if box is None:  # no box given
            # read file
            try:
                landmask = fi(landmaskname)
            except Exception:
                try:
                    landmask = fi('landmask')
                except Exception:
                    try:
                        landmask = fi('lsmask')
                    except Exception:
                        try:
                            landmask = fi('sftlf')
                        except Exception:
                            landmask = None
        else:  # box given by the user
            # define box
            region_ref = ReferenceRegions(box)

            # read file
            try:
                landmask = fi(landmaskname, latitude=region_ref['latitude'], longitude=region_ref['longitude'])
            except Exception:
                try:
                    landmask = fi('landmask', latitude=region_ref['latitude'], longitude=region_ref['longitude'])
                except Exception:
                    try:
                        landmask = fi('lsmask', latitude=region_ref['latitude'], longitude=region_ref['longitude'])
                    except Exception:
                        try:
                            landmask = fi('sftlf', latitude=region_ref['latitude'], longitude=region_ref['longitude'])
                        except Exception:
                            landmask = None

        fi.close()
    else:
        landmask = None

    _tg = tab.getGrid()
    _lg = landmask.getGrid() if landmask is not None else None

    if OSpath__isfile(filename) is False or landmask is None or _tg is None or _lg is None or _tg.shape != _lg.shape:
        # Estimate landmask
        landmask = EstimateLandmask(tab)

        if landmask is not None:
            try:
                check_grid_consistency(tab, landmask, context="ReadLandmaskSelectRegion")
            except Exception:
                landmask = None

        if landmask is not None and box is not None:
            # define box
            region_ref = ReferenceRegions(box)

            # subset
            landmask = landmask(latitude=region_ref['latitude'], longitude=region_ref['longitude'])

    # Return
    # Ensure landmask has valid spatial axes; no time axis is expected.
    if landmask is not None:
        landmask = _finalize_existing_cdat(landmask, context="ReadLandmaskSelectRegion")

    return landmask

def EstimateLandmask(d):
    """
    #################################################################################
    Description:
    Estimate a land-sea mask when no native landmask is provided.

    The refactored implementation uses the ENSO_metrics compatibility layer.
    The land mask is generated through the modern regionmask-based
    ``generateLandSeaMask`` pathway where available, while preserving the legacy
    ENSO_metrics land-fraction convention expected downstream.
    #################################################################################

    :param d: CDATVariable or masked-array-like
        Input model variable with latitude/longitude axes.

    :return landmask: CDATVariable
        Estimated land fraction mask on the input horizontal grid, with
        ``id='sftlf'``.
    """
    print('\033[93m' + str().ljust(25) + 'NOTE: Estimated landmask applied' + '\033[0m')
    n = 1
    sft = cdutil.generateLandSeaMask(d(*(slice(0, 1),) * n),debug=False) * 100.0
    sft[:] = sft.filled(100.0)
    lmsk = sft
    lmsk.setAxis(0, d.getAxis(1))
    lmsk.setAxis(1, d.getAxis(2))
    lmsk.id = 'sftlf'
    return lmsk

def Regrid(tab_to_regrid, newgrid, missing=None, order=None, mask=None,
           regridder='xesmf', regridTool='esmf', regridMethod='bilinear', **kwargs):
    """
    #################################################################################
    Description:
    Regrid ``tab_to_regrid`` to ``newgrid`` using the modern xarray/xESMF/ESMF
    regridding backend through the ENSO_metrics compatibility layer.

    This function preserves the legacy ENSO_metrics calling interface, but the
    active backend is now xESMF/ESMF rather than CDAT/cdms2.
    ``ReferenceRegions`` is expected to define longitude bounds in 0–360
    convention. When ``newgrid`` is constructed from ``newgrid_name`` and
    ``region``, the target longitude axis is kept monotonic and within
    ``[0, 360)``.
    #################################################################################

    :param tab_to_regrid: CDATVariable or masked-array-like
        Input field to regrid. In the refactored implementation this is expected
        to be a CDAT-like compatibility object, usually
        ``XarrayCompat.CDATVariable``, with valid latitude/longitude axes and a
        rectilinear grid.

    :param newgrid: grid-like object, string, or None
        Destination rectilinear grid. If ``newgrid`` is ``None`` or a string,
        the destination grid is constructed from ``newgrid_name`` and
        ``region``.

    :param missing: float, optional
        Missing-data value, if any. Retained for API compatibility.

    :param order: string, optional
        Axis order, for example ``"tzyx"`` or ``"tyx"``. Retained for API
        compatibility.

    :param mask: array of booleans, optional
        Optional mask for the destination grid. The mask may be 2-D or have the
        same shape as the target field, depending on the downstream use.

    :param regridder: string, optional
        Regridding backend name. The modern supported/default value is
        ``"xesmf"``.
        default value is ``"xesmf"``.

    :param regridTool: string, optional
        Regridding tool identifier. Retained for compatibility with existing
        ENSO_metrics configuration dictionaries. The modern expected value is
        ``"esmf"``, accessed through xESMF when available.
        default value is ``"esmf"``.

    :param regridMethod: string, optional
        Regridding method. The modern default is ``"bilinear"``.
        Supported method names are normalized internally; legacy ``"linear"``
        is treated as an alias for ``"bilinear"`` where supported.
        Common values include ``"bilinear"``, ``"conservative"``, and
        ``"nearest_s2d"``, depending on the backend.
        default value is ``"bilinear"``.

    usual kwargs:
    :param newgrid_name: string, optional
        Name used to construct the destination grid when ``newgrid`` is not
        provided explicitly. The name should specify a supported grid type and
        resolution, for example ``"generic_1x1deg"`` or ``"generic 1x1deg"``.
        default value is ``"generic_1x1deg"``.

    :param region: string, optional
        Region/domain used when constructing the destination grid. The region
        name must be defined in ``EnsoCollectionsLib.ReferenceRegions``. In this
        workflow, region longitude bounds are expected to be in 0–360
        convention.

    :return new_tab: CDATVariable
        ``tab_to_regrid`` regridded onto ``newgrid``, with CDAT-like metadata
        preserved for downstream ENSO metric calculations.
    """
    known_args = {"newgrid_name", "region"}
    extra_args = set(kwargs) - known_args
    if extra_args:
        EnsoErrorsWarnings.unknown_key_arg(extra_args, INSPECTstack())
    # test given arguments
    known_regridder = ["xesmf"]
    if regridder not in known_regridder:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": regridder",
            str().ljust(5) + "unknown regridder: " + str(regridder),
            str().ljust(10) + "known regridder: " + str(known_regridder),
        ]
        EnsoErrorsWarnings.my_error(list_strings)
    # test the given 'newgrid'
    if isinstance(newgrid, str) or newgrid is None:
        #
        # newgrid is not a grid, so a grid will be created
        # to do this, kwargs['newgrid_name'] and kwargs['region'] must be defined
        #
        if isinstance(newgrid, str) or newgrid is None:
            if "newgrid_name" not in kwargs or "region" not in kwargs:
                raise ValueError(
                    "Regrid: newgrid is None or a string, so both 'newgrid_name' and "
                    "'region' must be provided to construct the target grid."
                )

        # define the grid type
        for gtype in ["equalarea", "gaussian", "generic", "uniform"]:
            if gtype in kwargs['newgrid_name']:
                GridType = gtype
                break
            else:
                GridType = "generic"
        # define resolution (same resolution in lon and lat)
        for res in [
                "0.25x0.25deg", "0.5x0.5deg", "0.75x0.75deg", "1x1deg",
                "1.25x1.25deg", "1.5x1.5deg","1.75x1.75deg", "2x2deg",
                "2.25x2.25deg", "2.5x2.5deg", "2.75x2.75deg"
            ]:
            if res in kwargs['newgrid_name']:
                if res == "0.25x0.25deg":
                    GridRes = 0.25
                elif res == "0.5x0.5deg":
                    GridRes = 0.5
                elif res == "0.75x0.75deg":
                    GridRes = 0.75
                elif res == "1x1deg":
                    GridRes = 1.
                elif res == "1.25x1.25deg":
                    GridRes = 1.25
                elif res == "1.5x1.5deg":
                    GridRes = 1.5
                elif res == "1.75x1.75deg":
                    GridRes = 1.75
                elif res == "2x2deg":
                    GridRes = 2.
                elif res == "2.25x2.25deg":
                    GridRes = 2.25
                elif res == "2.5x2.5deg":
                    GridRes = 2.5
                else:
                    GridRes = 2.75
                break
        else:
            GridRes = 1.

        # Define bounds of 'region'
        region_ref = ReferenceRegions(kwargs["region"])
        lat1, lat2 = region_ref["latitude"][0], region_ref["latitude"][1]
        lon1, lon2 = region_ref["longitude"][0], region_ref["longitude"][1]

        # Ensure increasing latitude bounds
        lat_min, lat_max = sorted([float(lat1), float(lat2)])

        # Longitude bounds should already be in 0–360 convention from ReferenceRegions
        lon1 = float(lon1)
        lon2 = float(lon2)

        if lon1 < 0.0 or lon2 < 0.0 or lon1 >= 360.0 or lon2 > 360.0:
            raise ValueError(
                f"Regrid: ReferenceRegions must use 0–360 longitude bounds. "
                f"Got region={kwargs.get('region')!r}, lon=({region_ref['longitude'][0]}, {region_ref['longitude'][1]})."
            )

        if lon2 <= lon1:
            raise ValueError(
                f"Regrid: region={kwargs.get('region')!r} has non-increasing longitude bounds "
                f"lon=({region_ref['longitude'][0]}, {region_ref['longitude'][1]}). "
                "ReferenceRegions should define regions as non-crossing 0–360 intervals "
                "for target-grid construction."
            )

        nlat = int(round((lat_max - lat_min) / GridRes))
        nlon = int(round((lon2 - lon1) / GridRes))

        if nlat <= 0 or nlon <= 0:
            raise ValueError(
                f"Regrid: invalid target grid from newgrid_name={kwargs.get('newgrid_name')!r}, "
                f"region={kwargs.get('region')!r}, "
                f"lat=({lat1}, {lat2}), lon=({region_ref['longitude'][0]}, {region_ref['longitude'][1]}), "
                f"GridRes={GridRes}. Computed nlat={nlat}, nlon={nlon}."
            )

        # Create target-axis values directly
        lat_vals = lat_min + (GridRes / 2.0) + np.arange(nlat) * GridRes
        lon_vals = lon1 + (GridRes / 2.0) + np.arange(nlon) * GridRes

        if np.nanmin(lon_vals) < 0.0 or np.nanmax(lon_vals) >= 360.0:
            raise ValueError(
                f"Regrid: constructed longitude centers are outside [0, 360): "
                f"{np.nanmin(lon_vals)} to {np.nanmax(lon_vals)} for "
                f"region={kwargs.get('region')!r}."
            )

        # Create axes with explicit CF metadata
        lat = create_axis(
            lat_vals,
            id="lat",
            units="degrees_north",
            attributes={"axis": "Y", "standard_name": "latitude"},
        )

        lon = create_axis(
            lon_vals,
            id="lon",
            units="degrees_east",
            attributes={"axis": "X", "standard_name": "longitude"},
        )

        # Create grid
        newgrid = create_rect_grid(lat, lon, "yx", grid_type=GridType, mask=None)
        newgrid.id = kwargs["newgrid_name"]

    #
    # regrid
    #
    # Map regridMethod to the xESMF equivalent once, reused by both paths
    xesmf_method = _REGRID_METHOD_MAP.get(str(regridMethod).lower(), "bilinear") if regridMethod else "bilinear"
    regridFCT = REGRID2horizontal__Horizontal(
        tab_to_regrid.getGrid(),
        newgrid,
        method=xesmf_method,
    )
    new_tab = regridFCT(tab_to_regrid)
    if mask is not None:
        target_mask = np.asarray(mask, dtype=bool)
        new_data = ma.array(_mv(new_tab), copy=True)
        if target_mask.shape != new_data.shape:
            target_mask = np.broadcast_to(target_mask, new_data.shape)
        new_data = ma.array(new_data, mask=ma.getmaskarray(new_data) | target_mask)
        new_tab = new_tab.copy()
        new_tab._data = new_data

    return new_tab


def SaveNetcdf(
        netcdf_name, var1=None, var1_attributes={}, var1_name='', var1_time_name=None, var2=None,
        var2_attributes={}, var2_name='', var2_time_name=None, var3=None, var3_attributes={}, var3_name='',
        var3_time_name=None, var4=None, var4_attributes={}, var4_name='', var4_time_name=None, var5=None,
        var5_attributes={}, var5_name='', var5_time_name=None, var6=None, var6_attributes={}, var6_name='',
        var6_time_name=None, var7=None, var7_attributes={}, var7_name='', var7_time_name=None, var8=None,
        var8_attributes={}, var8_name='', var8_time_name=None, var9=None, var9_attributes={}, var9_name='',
        var9_time_name=None, var10=None, var10_attributes={}, var10_name='', var10_time_name=None, var11=None,
        var11_attributes={}, var11_name='', var11_time_name=None, var12=None, var12_attributes={}, var12_name='',
        var12_time_name=None, frequency="monthly", global_attributes={}, **kwargs
    ):

    _out_dir = ntpath.dirname(netcdf_name) or "."
    _os_sn.makedirs(_out_dir, exist_ok=True)

    if OSpath__isfile(netcdf_name) is True:
        o = open_file(netcdf_name, "a")
    else:
        o = open_file(netcdf_name, "w+")

    try:
        if var1 is not None:
            if var1_name == '':
                var1_name = var1.id
            if var1_time_name is not None:
                var1 = TimeButNotTime(var1, var1_time_name, frequency)
            o.write(var1, attributes=var1_attributes, dtype="float32", id=var1_name)

        if var2 is not None:
            if var2_name == '':
                var2_name = var2.id
            if var2_time_name is not None:
                var2 = TimeButNotTime(var2, var2_time_name, frequency)
            o.write(var2, attributes=var2_attributes, dtype="float32", id=var2_name)

        if var3 is not None:
            if var3_name == '':
                var3_name = var3.id
            if var3_time_name is not None:
                var3 = TimeButNotTime(var3, var3_time_name, frequency)
            o.write(var3, attributes=var3_attributes, dtype="float32", id=var3_name)

        if var4 is not None:
            if var4_name == '':
                var4_name = var4.id
            if var4_time_name is not None:
                var4 = TimeButNotTime(var4, var4_time_name, frequency)
            o.write(var4, attributes=var4_attributes, dtype="float32", id=var4_name)

        if var5 is not None:
            if var5_name == '':
                var5_name = var5.id
            if var5_time_name is not None:
                var5 = TimeButNotTime(var5, var5_time_name, frequency)
            o.write(var5, attributes=var5_attributes, dtype="float32", id=var5_name)

        if var6 is not None:
            if var6_name == '':
                var6_name = var6.id
            if var6_time_name is not None:
                var6 = TimeButNotTime(var6, var6_time_name, frequency)
            o.write(var6, attributes=var6_attributes, dtype="float32", id=var6_name)

        if var7 is not None:
            if var7_name == '':
                var7_name = var7.id
            if var7_time_name is not None:
                var7 = TimeButNotTime(var7, var7_time_name, frequency)
            o.write(var7, attributes=var7_attributes, dtype="float32", id=var7_name)

        if var8 is not None:
            if var8_name == '':
                var8_name = var8.id
            if var8_time_name is not None:
                var8 = TimeButNotTime(var8, var8_time_name, frequency)
            o.write(var8, attributes=var8_attributes, dtype="float32", id=var8_name)

        if var9 is not None:
            if var9_name == '':
                var9_name = var9.id
            if var9_time_name is not None:
                var9 = TimeButNotTime(var9, var9_time_name, frequency)
            o.write(var9, attributes=var9_attributes, dtype="float32", id=var9_name)

        if var10 is not None:
            if var10_name == '':
                var10_name = var10.id
            if var10_time_name is not None:
                var10 = TimeButNotTime(var10, var10_time_name, frequency)
            o.write(var10, attributes=var10_attributes, dtype="float32", id=var10_name)

        if var11 is not None:
            if var11_name == '':
                var11_name = var11.id
            if var11_time_name is not None:
                var11 = TimeButNotTime(var11, var11_time_name, frequency)
            o.write(var11, attributes=var11_attributes, dtype="float32", id=var11_name)

        if var12 is not None:
            if var12_name == '':
                var12_name = var12.id
            if var12_time_name is not None:
                var12 = TimeButNotTime(var12, var12_time_name, frequency)
            o.write(var12, attributes=var12_attributes, dtype="float32", id=var12_name)

        my_keys = sorted(
            [
                key for key in list(kwargs.keys())
                if "var" in key and str(key.replace("var", "")).isdigit() is True
            ],
            key=lambda v: v.upper(),
        )

        for key in my_keys:
            if kwargs[key] is not None:
                if key + "_name" not in list(kwargs.keys()) or \
                        (key + "_name" in list(kwargs.keys()) and kwargs[key + "_name"] == ''):
                    kwargs[key + "_name"] = kwargs[key].id

                if key + "_time_name" in list(kwargs.keys()) and kwargs[key + "_time_name"] is not None:
                    kwargs[key] = TimeButNotTime(kwargs[key], kwargs[key + "_time_name"], frequency)

                if key + "_attributes" not in list(kwargs.keys()):
                    kwargs[key + "_attributes"] = {}

                o.write(
                    kwargs[key],
                    attributes=kwargs[key + "_attributes"],
                    dtype="float32",
                    id=kwargs[key + "_name"],
                )

        for att in sorted(list(global_attributes.keys()), key=lambda v: v.upper()):
            o.__setattr__(att, global_attributes[att])

        o.close()

    except Exception:
        try:
            o.close()
        except Exception:
            pass
        raise

    return

def SkewnessTemporal(tab):
    """
    #################################################################################
    Description:
    Computes the skewness along the time axis
    #################################################################################

    :param tab: masked_array
    :return: tab: array
        array of the temporal skewness
    """
    tab = _to_cdat(tab)
    tab = tab.reorder('t...')
    if len(tab.shape) > 4:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": too many dimensions",
            str().ljust(5) + "tab.shape = " + str(tab.shape)]
        EnsoErrorsWarnings.my_error(list_strings)
    if len(tab.shape) == 1:
        valid = ma.compressed(ma.masked_invalid(_mv(tab)))
        skew = float(SCIPYstats__skew(valid)) if len(valid) > 0 else float('nan')
    else:
        if len(tab.shape) == 2:
            raw2d = ma.masked_invalid(_mv(tab))
            skew = SCIPYstats__skew(raw2d.filled(np.nan), axis=0, nan_policy='omit')
        else:
            # switch to numpy
            dataset = ma.masked_invalid(_mv(tab))
            # masked values -> nan
            dataset = dataset.filled(fill_value=NPnan)
            # Store information about the shape/size of the input data
            time_ax = dataset.shape[0]
            spac_ax = dataset.shape[1:]
            channels = NPproduct(spac_ax)
            # Reshape to two dimensions (time, space) creating the design matrix
            dataset = dataset.reshape([time_ax, channels])
            # Find the indices of values that are not missing in one row. All the rows will have missing values in the
            # same places provided the array was centered. If it wasn't then it is possible that some missing values
            # will be missed and the singular value decomposition will produce not a number for everything.
            nonMissingIndex = NPwhere(~NPisnan(dataset[0]))[0]
            # Remove missing values from the design matrix.
            dataNoMissing = dataset[:, nonMissingIndex]
            new_dataset = SCIPYstats__skew(dataNoMissing, axis=0)
            flatE = NPones([channels], dtype=dataset.dtype) * NPnan
            flatE = flatE.astype(dataset.dtype)
            flatE[nonMissingIndex] = new_dataset
            skew = flatE.reshape(spac_ax)
            skew = MV2masked_where(NPisnan(skew), skew)
        skew = create_variable(
            MV2array(skew),
            axes=tab.getAxisList()[1:],
            grid=tab.getGrid(),
            mask=tab[0].mask,
            attributes=tab.attributes,
            id='skewness'
        )
    return skew


def SmoothGaussian(tab, axis=0, window=5):
    """
    #################################################################################
    Description:
    Smooth 'tab' along 'axis' using gaussian moving window average
    #################################################################################

    :param tab: masked_array
        masked_array to smooth
    :param axis: integer, optional
        axis along which to smooth the data
        default value is the first axis (0)
    :param window: odd integer, optional
        number of points used for the triangle moving window average
        default value is 5
    :return smoothed_tab: masked_array
        smoothed data
    """
    if window % 2 == 0:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": smoothing window (running mean)",
            str().ljust(5) + "the window of smoothing must be an odd number: " + str(window)]
        EnsoErrorsWarnings.my_error(list_strings)
    if axis > len(tab.shape) - 1:
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": axis",
                        str().ljust(5) + "axis number too big: " + str(axis)]
        EnsoErrorsWarnings.my_error(list_strings)
    # Reorder tab in order to put 'axis' in first position
    tab = _to_cdat(tab)
    indices = list(range(len(tab.shape)))
    indices.remove(axis)
    newOrder = str(axis)
    for ii in indices:
        newOrder = newOrder + str(ii)
    new_tab = tab.reorder(newOrder)

    # degree
    degree = window // 2

    # Create the gaussian weight array
    weight = list()
    for ii in list(range(window)):
        ii = ii - degree + 1
        frac = ii / float(window)
        gauss = float(1. / (NPexp((4 * frac) ** 2)))
        ww = MV2zeros(new_tab.shape[1:])
        ww.fill(gauss)
        weight.append(ww)
        del frac, gauss, ww
    weight = MV2array(weight)

    # Smoothing
    smoothed_tab = MV2zeros(new_tab.shape)
    smoothed_tab = smoothed_tab[:len(new_tab) - window + 1]
    for ii in range(len(smoothed_tab)):
        tmp1 = MV2array(new_tab[ii: ii + window])
        tmp2 = MV2masked_where(tmp1.mask, weight)
        tmp1 = MV2sum(tmp1 * tmp2, axis=0) / MV2sum(tmp2, axis=0)
        tmp1 = MV2masked_where(MV2sum(tmp2.mask.astype("f"), axis=0) / window > 0.5, tmp1)
        smoothed_tab[ii] = tmp1
        del tmp1, tmp2

    # Axes list
    axes0 = new_tab[degree: len(new_tab) - degree].getAxisList()[0]
    if len(tab.shape) > 1:
        axes = [axes0] + new_tab.getAxisList()[1:]
    else:
        axes = [axes0]
    smoothed_tab.setAxisList(axes)
    if tab.getGrid():
        try:
            smoothed_tab.setGrid(tab.getGrid())
        except Exception:
            pass

    # Reorder to the input order
    for ii in range(axis):
        smoothed_tab = smoothed_tab.reorder(newOrder)
    return smoothed_tab


def SmoothSquare(tab, axis=0, window=5):
    """
    #################################################################################
    Description:
    Smooth 'tab' along 'axis' using square moving window average
    #################################################################################
    :param tab: masked_array
        masked_array to smooth
    :param axis: integer, optional
        axis along which to smooth the data
        default value is the first axis (0)
    :param window: odd integer, optional
        number of points used for the square moving window average
        default value is 5
    :return smoothed_tab: masked_array
        smoothed data
    """
    if window % 2 == 0:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": smoothing window (running mean)",
            str().ljust(5) + "the window of smoothing must be an odd number: " + str(window)]
        EnsoErrorsWarnings.my_error(list_strings)
    if axis > len(tab.shape)-1:
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": axis",
                        str().ljust(5) + "axis number too big: " + str(axis)]
        EnsoErrorsWarnings.my_error(list_strings)

    # Reorder tab in order to put 'axis' in first position
    tab = _to_cdat(tab)
    indices = list(range(len(tab.shape)))
    indices.remove(axis)
    newOrder = str(axis)
    for ii in indices:
        newOrder = newOrder+str(ii)
    new_tab = tab.reorder(newOrder)

    # degree
    degree = window // 2

    # Create the weight array (uniform)
    weight = MV2ones((window,) + new_tab.shape[1:])

    # Smoothing
    smoothed_tab = MV2zeros(new_tab.shape)
    smoothed_tab = smoothed_tab[:len(new_tab) - window + 1]
    for ii in range(len(smoothed_tab)):
        tmp1 = MV2array(new_tab[ii: ii + window])
        tmp2 = MV2masked_where(tmp1.mask, weight)
        tmp1 = MV2sum(tmp1, axis=0) / MV2sum(tmp2, axis=0)
        tmp1 = MV2masked_where(MV2sum(tmp2.mask.astype("f"), axis=0)/window>0.5, tmp1)
        smoothed_tab[ii] = tmp1
        del tmp1, tmp2

    # Axes list
    axes0 = new_tab[degree: len(new_tab) - degree].getAxisList()[0]
    if len(tab.shape) > 1:
        axes = [axes0] + new_tab.getAxisList()[1:]
    else:
        axes = [axes0]
    smoothed_tab.setAxisList(axes)
    if tab.getGrid():
        try:
            smoothed_tab.setGrid(tab.getGrid())
        except Exception:
            pass

    # Reorder to the input order
    for ii in range(axis):
        smoothed_tab = smoothed_tab.reorder(newOrder)
    return smoothed_tab


def SmoothTriangle(tab, axis=0, window=5):
    """
    #################################################################################
    Description:
    Smooth 'tab' along 'axis' using triangle moving window average
    #################################################################################
    :param tab: masked_array
        masked_array to smooth
    :param axis: integer, optional
        axis along which to smooth the data
        default value is the first axis (0)
    :param window: odd integer, optional
        number of points used for the triangle moving window average
        default value is 5
    :return smoothed_tab: masked_array
        smoothed data
    """
    if window % 2 == 0:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": smoothing window (running mean)",
            str().ljust(5) + "the window of smoothing must be an odd number: " + str(window)]
        EnsoErrorsWarnings.my_error(list_strings)
    if axis > len(tab.shape)-1:
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": axis",
                        str().ljust(5) + "axis number too big: " + str(axis)]
        EnsoErrorsWarnings.my_error(list_strings)

    # Reorder tab in order to put 'axis' in first position
    tab = _to_cdat(tab)
    indices = list(range(len(tab.shape)))
    indices.remove(axis)
    newOrder = str(axis)
    for ii in indices:
        newOrder = newOrder+str(ii)
    new_tab = tab.reorder(newOrder)

    # degree
    degree = window // 2

    # Create the weight array (triangle)
    weight = list()
    for ii in range(0, (2 * degree)+1):
        ww = MV2zeros(new_tab.shape[1:])
        ww.fill(float(1 + degree - abs(degree - ii)))
        weight.append(ww)
        del ww
    weight = MV2array(weight)

    # Smoothing
    smoothed_tab = MV2zeros(new_tab.shape)
    smoothed_tab = smoothed_tab[:len(new_tab) - window + 1]
    for ii in range(len(smoothed_tab)):
        tmp1 = MV2array(new_tab[ii: ii + window])
        tmp2 = MV2masked_where(tmp1.mask, weight)
        tmp1 = MV2sum(tmp1 * tmp2, axis=0) / MV2sum(tmp2, axis=0)
        tmp1 = MV2masked_where(MV2sum(tmp2.mask.astype("f"), axis=0) / window > 0.5, tmp1)
        smoothed_tab[ii] = tmp1
        del tmp1, tmp2

    # Axes list
    axes0 = new_tab[degree: len(new_tab) - degree].getAxisList()[0]
    if len(tab.shape) > 1:
        axes = [axes0] + new_tab.getAxisList()[1:]
    else:
        axes = [axes0]
    smoothed_tab.setAxisList(axes)
    if tab.getGrid():
        try:
            smoothed_tab.setGrid(tab.getGrid())
        except Exception:
            pass

    # Reorder to the input order
    for ii in range(axis):
        smoothed_tab = smoothed_tab.reorder(newOrder)
    return smoothed_tab


# Dictionary of seasons
# NOTE: sea_dict is already built from _MONTH_MAP above; this block
# is intentionally removed to avoid shadowing that definition.


def SeasonalMean(tab, season, compute_anom=False):
    """
    Computes seasonal mean or seasonal anomaly.
    """
    tab = _to_cdat(tab)
    _require_time_axis(tab, "SeasonalMean")

    if season in list(sea_dict.keys()):
        if compute_anom is True:
            tab = sea_dict[season].departures(tab)
        else:
            tab = sea_dict[season](tab)
    else:
        EnsoErrorsWarnings.unknown_key_arg("season", season, sorted(list(sea_dict.keys())), INSPECTstack())

    if season == "DJF":
        tab = _to_cdat(tab)
        time_ax = _require_time_axis(tab, "SeasonalMean DJF")
        time_num = _get_time_axis_index(tab, "SeasonalMean DJF")

        if len(time_ax) > 1:
            time_ax[:] = time_ax[:] - (time_ax[1] - time_ax[0])
            tab.setAxis(time_num, time_ax)

    return tab


# Dictionary of smoothing methods
dict_smooth = {'gaussian': SmoothGaussian, 'square': SmoothSquare, 'triangle': SmoothTriangle}


def Smoothing(tab, info, axis=0, window=5, method='triangle'):
    """
    #################################################################################
    Description:
    Smooth 'tab' along 'axis' using moving window average based on 'method'
    #################################################################################

    :param tab: masked_array
        masked_array to smooth
    :param info: string
        information about what was done on tab
    :param axis: integer, optional
        axis along which to smooth the data
        default value is the first axis (0)
    :param window: odd integer, optional
        number of points used for the moving window average
        default value is 5
    :param method: string, optional
        smoothing method:
            'gaussian': gaussian shaped window
            'square':   square shaped window
            'triangle': triangle shaped window
    :return: smoothed_tab: masked_array
        smoothed data
    """
    try: dict_smooth[method]
    except Exception:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": smoothing method (running mean)",
            str().ljust(5) + "unkwown smoothing method: " + str(method),
            str().ljust(10) + "known smoothing method: " + str(
                sorted(list(dict_smooth.keys()), key=lambda v: v.upper()))]
        EnsoErrorsWarnings.my_error(list_strings)
        return None, info
    info = info + ', smoothing using a ' + str(method) + ' shaped window of ' + str(window) + ' points'
    return dict_smooth[method](tab, axis=axis, window=window), info


def SkewMonthly(tab):
    """
    Computes monthly skewness for each calendar month.
    """
    tab = _to_cdat(tab)
    _require_time_axis(tab, "SkewMonthly")

    initorder = tab.getOrder()
    tab = tab.reorder("t...")
    axes = tab.getAxisList()
    time_ax = _get_component_time(tab, "SkewMonthly")

    months = MV2array([tt.month for tt in time_ax])
    cyc = []

    for ii in range(12):
        tmp = tab.compress(months == (ii + 1), axis=0)
        tmp = SCIPYstats__skew(tmp)
        cyc.append(tmp)
        del tmp

    time = create_axis(list(range(12)), id="time", axis_type="T")
    skew = create_variable(
        MV2array(cyc),
        axes=[time] + axes[1:],
        grid=tab.getGrid(),
        attributes=tab.attributes,
    )
    skew = skew.reorder(initorder)

    time = create_axis(list(range(12)), id="months", axis_type="T")
    skew.setAxis(get_num_axis(skew, "time"), time)
    return skew


def StdMonthly(tab):
    """
    Computes monthly standard deviation for each calendar month.
    """
    tab = _to_cdat(tab)
    _require_time_axis(tab, "StdMonthly")

    initorder = tab.getOrder()
    tab = tab.reorder("t...")
    axes = tab.getAxisList()
    time_ax = _get_component_time(tab, "StdMonthly")

    months = MV2array([tt.month for tt in time_ax])
    cyc = []

    for ii in range(12):
        tmp = tab.compress(months == (ii + 1), axis=0)
        tmp = Std(tmp, axis=0)
        cyc.append(tmp)
        del tmp

    time = create_axis(list(range(12)), id="time", axis_type="T")
    std = create_variable(
        MV2array(cyc),
        axes=[time] + axes[1:],
        grid=tab.getGrid(),
        attributes=tab.attributes,
    )
    std = std.reorder(initorder)

    time = create_axis(list(range(12)), id="months", axis_type="T")
    std.setAxis(get_num_axis(std, "time"), time)
    return std


def TimeButNotTime(tab, new_time_name, frequency):
    """
    Replace the time axis with a non-time axis while preserving length.
    """
    tab_out = copy.copy(_to_cdat(tab))

    time_num = _get_time_axis_index(tab_out, "TimeButNotTime")
    timeax = _get_component_time(tab_out, "TimeButNotTime")

    year1, month1, day1 = timeax[0].year, timeax[0].month, timeax[0].day

    if frequency == "daily":
        freq = "days"
    elif frequency == "monthly":
        freq = "months"
    elif frequency == "yearly":
        freq = "years"
    else:
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())
        freq = frequency

    axis = create_axis(list(range(len(tab_out))), id=new_time_name, axis_type="-")
    axis.units = f"{freq} since {year1}-{month1}-{day1}"
    axis.axis = freq

    tab_out.setAxis(time_num, axis)
    return tab_out
# ---------------------------------------------------------------------------------------------------------------------#


# ---------------------------------------------------------------------------------------------------------------------#
#
# Set of often used combinations of previous functions
#
def ComputePDF(tab, nbr_bins=10, interval=None, axis_name='axis'):
    """
    #################################################################################
    Description:
    Compute a probability density function or histogram from the input values.

    Missing values are excluded before binning. This function preserves the
    legacy ENSO_metrics interface while using the modern numpy-based statistics
    pathway.
    #################################################################################

    :param tab: CDATVariable or masked-array-like
        Input values used to compute the distribution.

    :param nbr_bin: integer, optional
        Number of bins to use.
        default value = None

    :param interval: list or tuple, optional
        Lower and upper bounds for the histogram interval.
        default value = None

    :param normalized: boolean, optional
        If True, normalize the histogram as a probability density.
        default value = True

    :return pdf, bins:
        Computed distribution and bin information, following the existing
        ENSO_metrics return convention.
    """
    tmp = NPhistogram(tab, bins=nbr_bins, range=interval)
    axis = [(tmp[1][ii] + tmp[1][ii + 1]) / 2. for ii in list(range(len(tmp[1]) - 1))]
    pdf = MV2array(tmp[0]) / float(len(tab))
    axis = create_axis(MV2array(axis, dtype='f'), id=axis_name)
    pdf.setAxis(0, axis)
    return pdf


def CustomLinearRegression(y, x, sign_x=0, return_stderr=True, return_intercept=True):
    """
    #################################################################################
    Description:
    Custom linear-regression helper used by ENSO_metrics.

    This function preserves the legacy ENSO_metrics interface while using the
    modern compatibility-layer replacement for genutil.linearregression. It can
    compute the regression using all x values, or only values with x >= 0 or
    x <= 0 through ``sign_x``.

    Regression is applied along the first axis. All remaining axes are
    preserved, which is required for longitude-profile and Hovmoeller outputs.

    The function accepts CDAT-like compatibility variables produced by the
    refactored xarray/XarrayCompat workflow, as well as masked-array-like inputs
    supported by the ENSO_metrics compatibility layer.
    #################################################################################

    :param y: CDATVariable or masked-array-like
        Dependent variable. In the refactored implementation this is typically
        an ``XarrayCompat.CDATVariable`` or compatible masked-array-like object
        with metadata preserved through the compatibility layer.

    :param x: CDATVariable or masked-array-like
        Independent variable used for the regression. Must be shape-compatible
        with ``y``.

    :param sign_x: int, optional
        Default value = 0. If 0, computes the regression of y over x using all
        valid points. If 1, computes the regression using x > 0. If -1,
        computes the regression using x < 0, following ``CustomLinearRegression1d``.

    :param return_stderr: boolean, optional
        Default value = True. If True, returns the unadjusted standard error of
        the regression slope.

    :param return_intercept: boolean, optional
        Default value = True. If True, returns the regression intercept.

    :return tab: float, CDATVariable, or list
        If both ``return_stderr`` and ``return_intercept`` are False, returns
        only the slope. Otherwise returns a list containing the slope, and
        optionally the standard error and intercept.
    """
    if sign_x not in [-1, 0, 1]:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": sign_x",
            str().ljust(5) + "unknown sign_x " + str(sign_x),
            str().ljust(5) + "known values are -1, 0, 1",
        ]
        EnsoErrorsWarnings.my_error(list_strings)

    try:
        len(y[0])
    except Exception:
        # 1-D input: direct regression.
        slope, intercept, stderr = CustomLinearRegression1d(
            y,
            x,
            sign_x=sign_x,
        )

    else:
        # Multi-dimensional input: perform pointwise regression along the first
        # axis and preserve all remaining axes.
        if x.shape != y.shape:
            list_strings = [
                "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": array shape",
                str().ljust(5) + "different array shape for x " + str(x.shape) + " and y " + str(y.shape),
            ]
            EnsoErrorsWarnings.my_error(list_strings)

        slope = MV2zeros(y[0].shape)
        intercept = MV2zeros(y[0].shape)
        stderr = MV2zeros(y[0].shape)

        for ii in list(range(len(y[0]))):
            try:
                len(y[0, ii])
            except Exception:
                slope[ii], intercept[ii], stderr[ii] = CustomLinearRegression1d(
                    y[:, ii],
                    x[:, ii],
                    sign_x=sign_x,
                )

            else:
                for jj in list(range(len(y[0, ii]))):
                    try:
                        len(y[0, ii, jj])
                    except Exception:
                        slope[ii, jj], intercept[ii, jj], stderr[ii, jj] = CustomLinearRegression1d(
                            y[:, ii, jj],
                            x[:, ii, jj],
                            sign_x=sign_x,
                        )

                    else:
                        for kk in list(range(len(y[0, ii, jj]))):
                            try:
                                len(y[0, ii, jj, kk])
                            except Exception:
                                slope[ii, jj, kk], intercept[ii, jj, kk], stderr[ii, jj, kk] = (
                                    CustomLinearRegression1d(
                                        y[:, ii, jj, kk],
                                        x[:, ii, jj, kk],
                                        sign_x=sign_x,
                                    )
                                )

                            else:
                                list_strings = [
                                    "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) +
                                    ": array shape",
                                    str().ljust(5) + str(x.shape) + " too many dimensions (not programmed)",
                                    str().ljust(5) + "Please check and modify the program if needed",
                                ]
                                EnsoErrorsWarnings.my_error(list_strings)

    try:
        len(slope)
    except Exception:
        # Scalar result: keep legacy scalar behavior.
        pass

    else:
        # Preserve the axes/grid/mask of y after removing the regression axis.
        y0 = _to_cdat(y[0])
        axes = y0.getAxisList()
        grid = y0.getGrid()
        mask = ma.getmaskarray(_mv(y0))

        slope = create_variable(
            _mv(slope),
            mask=mask,
            grid=grid,
            axes=axes,
            id='slope',
        )
        stderr = create_variable(
            _mv(stderr),
            mask=mask,
            grid=grid,
            axes=axes,
            id='standard_error',
        )
        intercept = create_variable(
            _mv(intercept),
            mask=mask,
            grid=grid,
            axes=axes,
            id='intercept',
        )

    if return_stderr is False and return_intercept is False:
        tab = copy.copy(slope)
    else:
        tab = [slope]
        if return_stderr is True:
            tab.append(stderr)
        if return_intercept is True:
            tab.append(intercept)

    return tab

def CustomLinearRegression1d(y, x, sign_x=1):
    x = np.ma.masked_invalid(NParray(x))
    y = np.ma.masked_invalid(NParray(y))
    if sign_x == 1:
        idx = NPnonzero((x > 0.) & ~np.ma.getmaskarray(x) & ~np.ma.getmaskarray(y))
    elif sign_x == -1:
        idx = NPnonzero((x < 0.) & ~np.ma.getmaskarray(x) & ~np.ma.getmaskarray(y))
    else:
        idx = NPnonzero(~np.ma.getmaskarray(x) & ~np.ma.getmaskarray(y))
    if len(idx[0]) == 0:
        slope, intercept, stderr = np.nan, np.nan, np.nan
    else:
        results = GENUTILlinearregression(y[idx], x=x[idx], error=1, nointercept=None)
        slope, intercept, stderr = float(results[0][0][0]), float(results[0][0][1]), float(results[1][0][0])
    return slope, intercept, stderr


def fill_dict_teleconnection(
        tab1, tab2, dataset1, dataset2, timebounds1, timebounds2, nyear1, nyear2, nbr, var_name,
        add_name, units, centered_rmse=0, biased_rmse=1, dict_metric={}, dict_nc={}, ev_name=None,
        events1=None, events2=None
    ):
    # Metric 1
    rmse_dive, keyerror = RmsAxis(tab1, tab2, axis="xy", centered=centered_rmse, biased=biased_rmse)
    rmse_error_dive = None
    # Metric 2
    corr_dive = float(Correlation(tab1, tab2, axis="xy", centered=1, biased=1))
    corr_error_dive = None
    # Metric 3
    std_mod_dive = float(Std(tab1, weights=None, axis="xy", centered=1, biased=1))
    std_obs_dive = float(Std(tab2, weights=None, axis="xy", centered=1, biased=1))
    std_dive = std_mod_dive / std_obs_dive
    std_error_dive = None
    list_met_name = [
        "RMSE_" + dataset2, "RMSE_error_" + dataset2,
        "CORR_" + dataset2, "CORR_error_" + dataset2,
        "STD_" + dataset2, "STD_error_" + dataset2
    ]
    list_metric_value = [float(rmse_dive), rmse_error_dive, corr_dive, corr_error_dive, std_dive, std_error_dive]
    for tmp1, tmp2 in zip(list_met_name, list_metric_value):
        dict_metric[tmp1 + "_" + add_name] = tmp2
    dict_nc["var" + str(nbr)] = tab1
    dict_dive = {
        "units": units,
        "number_of_years_used": nyear1,
        "time_period": str(timebounds1),
        "spatialSTD_" + dataset1: std_mod_dive
    }
    if isinstance(events1, list) is True:
        dict_dive[ev_name + "_years"] = str(events1)
    dict_nc["var" + str(nbr) + "_attributes"] = dict_dive
    dict_nc["var" + str(nbr) + "_name"] = var_name + dataset1
    dict_dive = {
        "units": units, "number_of_years_used": nyear2, "time_period": str(timebounds2),
        "spatialSTD_" + dataset2: std_obs_dive
    }
    if isinstance(events2, list) is True:
        dict_dive[ev_name + "_years"] = str(events2)
    dict_nc["var" + str(nbr + 1)] = tab2
    dict_nc["var" + str(nbr + 1) + "_attributes"] = dict_dive
    dict_nc["var" + str(nbr + 1) + "_name"] = var_name + dataset2
    return dict_metric, dict_nc


def FindXYMinMaxInTs(tab, return_val='both', smooth=False, axis=0, window=5, method='triangle'):
    """
    #################################################################################
    Description:
    Find the spatial locations of minimum and maximum values in a
    time-dependent field.

    This helper preserves the legacy ENSO_metrics interface while using
    CDAT-like compatibility objects provided by the refactored xarray workflow.
    #################################################################################

    :param tab: CDATVariable or masked-array-like
        Time-dependent field with latitude/longitude axes.

    :param metric_name: string
        Metric name used for diagnostic messages and context.

    :return:
        Coordinates or values associated with the spatial minima and maxima,
        following the existing ENSO_metrics return convention.
    """
    tab_ts = list()
    for tt in list(range(len(tab))):
        if smooth is True:
            tmp, unneeded = Smoothing(tab[tt], '', axis=axis, window=window, method=method)
        else:
            tmp = copy.copy(tab[tt])
        tab_ts.append(find_xy_min_max(tmp, return_val=return_val))
    tab_ts = MV2array(tab_ts)
    tab_ts.setAxis(0, tab.getAxis(0))
    return tab_ts


def MyDerive(project, internal_variable_name, dict_var):
    # get dictionary of observations
    dict_obs = ReferenceObservations()
    # test input parameters
    keyerror1, keyerror2, keyerror3, keyerror4 = None, None, None, None
    if not isinstance(project, str):
        keyerror1 = "project is not well defined (" + str(project) + ")"
        EnsoErrorsWarnings.object_type_error('project', 'string', type(project), INSPECTstack())
    if not isinstance(internal_variable_name, str):
        keyerror2 = "internal_variable_name is not well defined (" + str(internal_variable_name) + ")"
        EnsoErrorsWarnings.object_type_error(
            'internal_variable_name', 'string', type(internal_variable_name),
            INSPECTstack()
        )
    if not isinstance(dict_var, dict):
        keyerror3 = "dictionary of variable is not well defined (" + str(internal_variable_name) + ")"
        EnsoErrorsWarnings.object_type_error('project', 'dictionary', type(dict_var), INSPECTstack())
    # wrong project?
    if 'CMIP' not in project and project not in list(dict_obs.keys()) and keyerror1 is None:
        keyerror4 = "project is not well defined (" + str(project) + ")"
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": project",
            str().ljust(5) + "unknown 'project' (or observations dataset): " + str(project),
            str().ljust(10) + "it must be either a 'CMIP' project or an observations dataset defined in " +
            "EnsoCollectionsLib.ReferenceObservations",
            str().ljust(10) + "known observations dataset: " + str(
                sorted(list(dict_obs.keys()), key=lambda v: v.upper()))]
        EnsoErrorsWarnings.my_warning(list_strings)

    precomputed_aliases = {
        "thf": ["thf", "netflux", "hfds", "thflx", "sohefldo"],
    }
    if internal_variable_name in precomputed_aliases and isinstance(dict_var, dict):
        for alias in precomputed_aliases[internal_variable_name]:
            if alias in dict_var:
                return dict_var[alias], None

    if (keyerror1 is not None or keyerror2 is not None or keyerror3 is not None or keyerror4 is not None):
        outvar = None
        keyerror = add_up_errors([keyerror1, keyerror2, keyerror3, keyerror4])
    else:
        # compute 'internal_variable_name' in 'CMIP' case
        if 'CMIP' in project:
            # get dictionary of CMIP
            dict_CMIP = CmipVariables()['variable_name_in_file']
            # test if 'internal_variable_name' is defined in EnsoCollectionsLib.CmipVariables
            if internal_variable_name in list(dict_CMIP.keys()):
                list_var = dict_CMIP[internal_variable_name]['var_name']
                outvar, keyerror = MyDeriveCompute(
                    list_var, dict_var, dict_att=dict_CMIP, variable=internal_variable_name, project=project)
            else:
                outvar = None
                keyerror = "variable (" + str(internal_variable_name) + ") not defined in CmipVariables"
        # compute 'internal_variable_name' in 'obs' case
        else:
            # 'project' is defined in EnsoCollectionsLib.ReferenceObservations
            dict_obs_var = dict_obs[project]['variable_name_in_file']
            # test if 'internal_variable_name' is defined for this observations dataset
            if internal_variable_name in list(dict_obs_var.keys()):
                list_var = dict_obs_var[internal_variable_name]['var_name']
                outvar, keyerror = MyDeriveCompute(
                    list_var, dict_var, dict_att=dict_obs_var, variable=internal_variable_name, isObs=True
                )
            else:
                outvar = None
                keyerror = "variable (" + str(internal_variable_name) + ") not defined in ReferenceObservations[" +\
                    str(project) + "]"
    return outvar, keyerror


def MyDeriveCompute(list_var, dict_var, dict_att={}, variable='', isObs=False, project=''):
    # test if keys in list_var are in 'dict_var'
    string_in_dict(list_var, dict_var, INSPECTstack())
    if isinstance(list_var, str):
        # this 'internal_variable_name' is based on one variable
        keyerror = None
        outvar = dict_var[list_var]
    else:
        # this 'internal_variable_name' is based on several variables
        list_operator = dict_att[variable]['algebric_calculation']
        if len(list_operator) != len(list_var):
            outvar = None
            keyerror = str(len(list_var)) + " variables are needed to compute " + str(variable) + " but " + \
                    str(len(list_operator)) + " operator(s) are given"
            if isObs is True:
                list_strings = [
                    "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) +
                    ": variable definition in EnsoCollectionsLib.ReferenceObservations(" + str(project) + ")",
                    str().ljust(5) + str(len(list_var)) + " variables are needed to compute " +
                    str(variable) + " but " + str(len(list_operator)) + " operator(s) are given"
                ]
            else:
                list_strings = [
                    "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) +
                    ": variable definition in EnsoCollectionsLib.CmipVariables",
                    str().ljust(5) + str(len(list_var)) + " variables are needed to compute " +
                    str(variable) + " but " + str(len(list_operator)) + " operator(s) are given"
                ]
            EnsoErrorsWarnings.my_warning(list_strings)
        else:
            keyerror = None
            # compute the output variable
            if list_operator[0] == 'minus':
                outvar = -1 * dict_var[list_var[0]]
            else:
                outvar = dict_var[list_var[0]]
            for ii in list(range(1, len(list_var))):
                outvar = dict_operations[list_operator[ii]](outvar, dict_var[list_var[ii]])
            outvar.setAxisList(dict_var[list_var[0]].getAxisList())
            outvar = MV2masked_where(dict_var[list_var[0]].mask, outvar)
            outvar.setGrid(dict_var[list_var[0]].getGrid())
    return outvar, keyerror


def LinearRegressionAndNonlinearity(y, x, return_stderr=True, return_intercept=True):
    """
    #################################################################################
    Description:
    Compute linear regressions of ``y`` against ``x`` for all values, positive
    values of ``x``, and negative values of ``x``.

    This function preserves the legacy ENSO_metrics interface while using the
    modern numpy/scipy-based regression utilities and CDAT-like compatibility
    objects provided by the refactored xarray workflow.
    #################################################################################

    :param y: CDATVariable or masked-array-like
        Dependent variable.

    :param x: CDATVariable or masked-array-like
        Independent variable used to split the regression into all, positive,
        and negative branches.

    :param return_stderr: boolean, optional
        If True, return the unadjusted standard error of the regression slope.
        default value = True

    :param return_intercept: boolean, optional
        If True, return the regression intercept.
        default value = True

    :return:
        Regression results for all, positive, and negative values of ``x``.
        Each result contains slope, optional standard error, and optional
        intercept depending on ``return_stderr`` and ``return_intercept``.
    """
    # all points
    all_values = CustomLinearRegression(y, x, 0, return_stderr=return_stderr, return_intercept=return_intercept)
    # positive SSTA = El Nino
    positive_values = CustomLinearRegression(y, x, 1, return_stderr=return_stderr, return_intercept=return_intercept)
    # negative SSTA = La Nina
    negative_values = CustomLinearRegression(y, x, -1, return_stderr=return_stderr, return_intercept=return_intercept)
    return all_values, positive_values, negative_values


def _linear_regression_nointercept_axis0(y, x):
    y_data = ma.masked_invalid(ma.asarray(_mv(y)))
    x_data = ma.masked_invalid(ma.asarray(_mv(x)))
    if x_data.shape != y_data.shape:
        x_data = ma.array(
            np.broadcast_to(ma.getdata(x_data), y_data.shape),
            mask=np.broadcast_to(ma.getmaskarray(x_data), y_data.shape),
        )
    joint_mask = ma.getmaskarray(y_data) | ma.getmaskarray(x_data)
    y_data = ma.array(y_data, mask=joint_mask)
    x_data = ma.array(x_data, mask=joint_mask)

    n = y_data.count(axis=0)
    den = ma.sum(x_data ** 2, axis=0)
    num = ma.sum(x_data * y_data, axis=0)
    valid = (n >= 2) & ~ma.getmaskarray(den) & (ma.filled(den, 0) != 0)
    slope_data = ma.masked_where(~valid, ma.divide(num, den))

    same_series = ma.max(ma.abs(y_data - x_data), axis=0) <= 1e-12
    same_series = ma.filled(same_series, False) & valid
    slope_data = ma.where(same_series, 1.0, slope_data)

    y_hat = x_data * slope_data
    resid = y_data - y_hat
    mse = ma.masked_where(~valid, ma.divide(ma.sum(resid ** 2, axis=0), n - 1))
    stderr_data = ma.sqrt(ma.divide(mse, den))
    stderr_data = ma.where(same_series, 0.0, stderr_data)
    return slope_data, stderr_data


def LinearRegressionTsAgainstMap(y, x, return_stderr=True):
    """
    #################################################################################
    Description:
    Compute the linear regression of a time series against a time-dependent map.

    This function preserves the legacy ENSO_metrics interface while using the
    modern numpy/scipy-based regression utilities and CDAT-like compatibility
    objects provided by the refactored xarray workflow.
    #################################################################################

    :param y: CDATVariable or masked-array-like
        Time-dependent map or field used as the dependent variable.

    :param x: CDATVariable or masked-array-like
        Independent time series used for regression.

    :param return_stderr: boolean, optional
        If True, return the unadjusted standard error of the regression slope.
        default value = True

    :return slope, stderr:
        Regression slope map and, if requested, the unadjusted standard-error map.
    """
    y = _to_cdat(y)
    x = _to_cdat(x)
    x_data = ma.masked_invalid(_mv(x))
    x_1d = x_data.ravel()
    expand = (slice(None),) + (np.newaxis,) * (ma.asarray(_mv(y)).ndim - 1)
    x_bc = x_1d[expand]
    slope_data, stderr_data = _linear_regression_nointercept_axis0(y, x_bc)
    # Attach spatial axes from y (drop the leading time axis)
    spatial_axes = [ax.copy() if ax is not None else None
                    for ax in y.getAxisList()[1:]]
    slope = CDATVariable(
        slope_data,
        axes=spatial_axes,
        grid=y.getGrid(),
        id=getattr(x, 'id', '')
        )
    if return_stderr:
        stderr = CDATVariable(
            stderr_data,
            axes=spatial_axes,
            grid=y.getGrid(),
            id=getattr(x, 'id', '')
        )
        return slope, stderr
    return slope


def LinearRegressionTsAgainstTs(y, x, nbr_years_window, return_stderr=True, frequency=None, debug=False):
    """
    #################################################################################
    Description:
    Compute the lead-lag linear regression of a time series against another
    time series.

    This function preserves the legacy ENSO_metrics interface while using the
    modern numpy/scipy-based regression utilities and CDAT-like compatibility
    objects provided by the refactored xarray workflow.
    #################################################################################

    :param y: CDATVariable or masked-array-like
        Dependent time series.

    :param x: CDATVariable or masked-array-like
        Independent time series used for lead-lag regression.

    :param nbr_years_window: integer
        Number of years used to compute the lead-lag regression window.

    :param return_stderr: boolean, optional
        If True, return the unadjusted standard error of the regression slope.
        default value = True

    :param frequency: string, optional
        Time frequency of the datasets, for example ``"monthly"``.
        default value = None

    :param debug: boolean, optional
        If True, print diagnostic information during calculation.
        default value = False

    :return slope, stderr:
        Regression slope and, if requested, the unadjusted standard error.
    """
    if frequency == 'daily':
        nbr_timestep = nbr_years_window * 365
    elif frequency == 'monthly':
        nbr_timestep = nbr_years_window * 12
    elif frequency == 'yearly':
        nbr_timestep = nbr_years_window
    else:
        EnsoErrorsWarnings.unknown_frequency(frequency, INSPECTstack())
    y = _to_cdat(y)
    x = _to_cdat(x)
    tab_yy_mm = Event_selection(y, frequency, nbr_years_window=nbr_years_window)
    myshape = [nbr_timestep] + [ss for ss in y.shape[1:]]
    tmp_ax = create_axis(list(range(nbr_timestep)), id='months')
    slope_out = MV2zeros(myshape)
    slope_out.setAxisList([tmp_ax] + y.getAxisList()[1:])
    stderr_out = MV2zeros(myshape)
    stderr_out.setAxisList([tmp_ax] + y.getAxisList()[1:])
    for ii in list(range(nbr_timestep)):
        tmp1 = tab_yy_mm[:, ii]
        tmp2 = copy.copy(x)
        yy1 = tab_yy_mm.getAxis(0)[0]
        yy2 = _require_time_axis(tmp2, 'LinearRegressionAndNonlinearity').asComponentTime()[0].year
        if yy1 == yy2:
            tmp1 = tmp1[:len(tmp2)]
        elif yy1 < yy2:
            tmp1 = tmp1[yy2 - yy1:len(tmp2)]
        else:
            tmp2 = tmp2[yy2 - yy1:]
            tmp1 = tmp1[:len(x)]
        if len(tmp2) > len(tmp1):
            tmp2 = tmp2[:len(tmp1)]
        # if debug is True:
        #     yy1 = tmp1.getAxis(0)[0]
        #     yy2 = tmp2.getTime().asComponentTime()[0].year
        #     EnsoErrorsWarnings.DebugMode('\033[93m', "EnsoUvcdatToolsLib LinearRegressionTsAgainstTs", 20)
        #     dict_debug = {'axes1': str([ax.id for ax in tmp1.getAxisList()]), 'shape1': str(tmp1.shape),
        #                   'line1': "first year is " + str(yy1),
        #                   'axes2': str([ax.id for ax in tmp2.getAxisList()]), 'shape2': str(tmp2.shape),
        #                   'line2': "first year is " + str(yy2)}
        #     EnsoErrorsWarnings.DebugMode('\033[93m', str(x.id) + " regressed against " + str(y.id), 25, **dict_debug)
        if tmp2.shape == tmp1.shape:
            tmp3 = copy.copy(tmp2)
        else:
            tmp3 = MV2zeros(tmp1.shape)
            for jj in list(range(len(tmp3))):
                tmp3[jj].fill(tmp2[jj])
        tmp3_mask = ma.getmaskarray(_mv(tmp3)) | ma.getmaskarray(_mv(tmp1))
        tmp3 = create_variable(tmp3, mask=tmp3_mask, grid=tmp1.getGrid(), axes=tmp1.getAxisList(), id=x.id)
        slope, stderr = _linear_regression_nointercept_axis0(tmp1, tmp3)
        slope_out[ii] = slope
        stderr_out[ii] = stderr
        del slope, stderr, tmp1, tmp2, tmp3, tmp3_mask, yy1, yy2
    if return_stderr:
        return slope_out, stderr_out
    else:
        return slope_out


def PreProcessTS(tab, info, areacell=None, average=False, compute_anom=False, compute_sea_cycle=False, debug=False,
                 region=None, **kwargs):
    keyerror = None
    # removes annual cycle (anomalies with respect to the annual cycle)
    if compute_anom is True:
        tab = ComputeInterannualAnomalies(tab)
    # Normalization of the anomalies
    if kwargs['normalization']:
        if kwargs['frequency'] is not None:
            tab, keyerror = Normalize(tab, kwargs['frequency'])
            info = info + ', normalized'
    # Removing linear trend
    if keyerror is None:
        if isinstance(kwargs['detrending'], dict):
            known_args = {'axis', 'method', 'bp'}
            extra_args = set(kwargs['detrending']) - known_args
            if extra_args:
                EnsoErrorsWarnings.unknown_key_arg(extra_args, INSPECTstack())
            tab, info, keyerror = Detrend(tab, info, **kwargs['detrending'])
    if keyerror is None:
        # Smoothing time series
        if isinstance(kwargs['smoothing'], dict):
            known_args = {'axis', 'method', 'window'}
            extra_args = set(kwargs['smoothing']) - known_args
            if extra_args:
                EnsoErrorsWarnings.unknown_key_arg(extra_args, INSPECTstack())
            tab, info = Smoothing(tab, info, **kwargs['smoothing'])
        # computes mean annual cycle
        if compute_sea_cycle is True:
            tab = annualcycle(tab)
        # average
        if average is not False:
            if debug is True:
                EnsoErrorsWarnings.debug_mode('\033[93m', "EnsoUvcdatToolsLib PreProcessTS", 20)
                dict_debug = {'axes1':  str([ax.id for ax in tab.getAxisList()]), 'shape1': str(tab.shape)}
                EnsoErrorsWarnings.debug_mode('\033[93m', "averaging to perform: " + str(average), 25, **dict_debug)
            if isinstance(average, str):
                if average not in dict_average:
                    EnsoErrorsWarnings.unknown_averaging(average, list(dict_average.keys()), INSPECTstack())
                else:
                    tab, keyerror = dict_average[average](tab, areacell, region=region, **kwargs)
                    if keyerror is None:
                        if debug is True:
                            dict_debug = {'axes1': str([ax.id for ax in tab.getAxisList()]), 'shape1': str(tab.shape)}
                            EnsoErrorsWarnings.debug_mode('\033[93m', "performed " + str(average), 25, **dict_debug)
            elif isinstance(average, list):
                for av in average:
                    if av not in dict_average:
                        EnsoErrorsWarnings.unknown_averaging(average, list(dict_average.keys()), INSPECTstack())
                    else:
                        tab, keyerror = dict_average[av](tab, areacell, region=region, **kwargs)
                        if keyerror is None:
                            if debug is True:
                                dict_debug = {
                                    'axes1': str([ax.id for ax in tab.getAxisList()]),
                                    'shape1': str(tab.shape)
                                    }
                                EnsoErrorsWarnings.debug_mode('\033[93m', "performed " + str(av), 25, **dict_debug)
                        else:
                            break
            else:
                EnsoErrorsWarnings.unknown_averaging(average, list(dict_average.keys()), INSPECTstack())
    else:
        tab = None
    return tab, info, keyerror


def ReadSelectRegionCheckUnits(filename, varname, varfamily, box=None, time_bounds=None, frequency=None, **keyarg):
    """
    #################################################################################
    Description:
    Read a variable, select the requested region and time period, check or
    harmonize units, and return the processed field for metric calculation.

    This function preserves the legacy ENSO_metrics interface while using the
    modern xarray/xcdat-compatible file-reading pathway and CDAT-like
    compatibility objects provided by the refactored workflow.
    #################################################################################

    :param filename: string or list
        Path or paths to the input NetCDF file(s).

    :param varname: string
        Name of the variable to read from ``filename``.

    :param varfamily: string
        Variable family used for unit checks, for example temperature,
        precipitation, velocity, or heat flux.

    :param box: string, optional
        Region name to select. The region must be defined in
        ``EnsoCollectionsLib.ReferenceRegions``.
        default value = None

    :param time_bounds: tuple, optional
        First and last dates to extract from the file(s), for example
        ``('1979-01-01T00:00:00', '2017-01-01T00:00:00')``.
        default value = None

    :param frequency: string, optional
        Temporal frequency of the input data, for example ``"monthly"``.
        default value = None

    usual kwargs:
        Additional options passed through the existing ENSO_metrics reading,
        selection, and preprocessing workflow.

    :return tab, keyerror:
        Processed field and any accumulated keyerror message.
    """
    tab = ReadAndSelectRegion(filename, varname, box=box, time_bounds=time_bounds, frequency=frequency)
    tab, units, keyerror = CheckUnits(tab, varfamily, varname, tab.units, return_tab_only=False)
    tab.name = varname
    tab.units = units
    # Final gate: validate the CDATVariable structure before returning to metric
    # computation.  Catches axis/shape mismatches introduced by CheckUnits or
    # any upstream mutation and gives a clear error rather than a silent None.
    validate_cdat_variable(tab, context=f"ReadSelectRegionCheckUnits:{varname}")
    return tab, keyerror


def Read_data_mask_area(file_data, name_data, type_data, metric, region, file_area='', name_area='', file_mask='',
                        name_mask='', maskland=False, maskocean=False, time_bounds=None, debug=False, **kwargs):
    keyerror1, keyerror2, keyerror3 = None, None, None
    # Read variable
    if debug is True:
        dict_debug = {'file1': '(' + type_data + ') ' + str(file_data), 'var1': '(' + type_data + ') ' + str(name_data)}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'Files', 20, **dict_debug)
    variable, keyerror1 = ReadSelectRegionCheckUnits(file_data, name_data, type_data, box=region,
                                                     time_bounds=time_bounds, **kwargs)
    if debug is True:
        dict_debug = {
            'axes1': '(' + type_data + ') ' + str([ax.id for ax in variable.getAxisList()]),
            'shape1': '(' + type_data + ') ' + str(variable.shape),
            'time1': '(' + type_data + ') ' + str(TimeBounds(variable))
            }
        EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadSelectRegionCheckUnits', 20, **dict_debug)
    # checks if the time-period fulfills the minimum length criterion
    if isinstance(kwargs['min_time_steps'], int):
        if len(variable) < kwargs['min_time_steps']:
            EnsoErrorsWarnings.too_short_time_period(metric, len(variable), kwargs['min_time_steps'], INSPECTstack())
            keyerror2 = "too short time period (" + str(len(variable)) + ")"
    # Read areacell & mask
    variable, areacell, keyerror3 = Read_mask_area(
        variable, name_data, file_data, type_data, region, file_area=file_area, name_area=name_area,
        file_mask=file_mask, name_mask=name_mask, maskland=maskland, maskocean=maskocean, debug=debug, **kwargs)
    if keyerror1 is not None or keyerror2 is not None or keyerror3 is not None:
        keyerror = add_up_errors([keyerror1, keyerror2, keyerror3])
    else:
        keyerror = None
    return variable, areacell, keyerror


def Read_data_mask_area_multifile(
                file_data, name_data, type_data, variable, metric, region, file_area='', name_area='',
                file_mask='', name_mask='', maskland=False, maskocean=False, debug=False,
                interpreter='', **kwargs
    ):
    dict_area, dict_keye, dict_var = dict(), dict(), dict()
    def safe_get(seq, idx):
        try:
            return seq[idx]
        except Exception:
            return ''
    if isinstance(file_data, str):
        tab, areacell, keyerror = \
            Read_data_mask_area(file_data, name_data, type_data, metric, region, file_area=file_area,
                                name_area=name_area, file_mask=file_mask, name_mask=name_mask, maskland=maskland,
                                maskocean=maskocean, debug=debug, **kwargs)
        dict_area[name_data], dict_keye[name_data], dict_var[name_data] = areacell, keyerror, tab
    else:
        for ii in range(len(file_data)):
            ff1 = safe_get(file_data, ii)
            nn1 = safe_get(name_data, ii)
            fa1 = safe_get(file_area, ii)
            an1 = safe_get(name_area, ii)
            fl1 = safe_get(file_mask, ii)
            ln1 = safe_get(name_mask, ii)
            tab, areacell, keyerror = Read_data_mask_area(
                ff1, nn1, type_data, metric, region,
                file_area=fa1, name_area=an1, file_mask=fl1,
                name_mask=ln1, maskland=maskland, maskocean=maskocean, debug=debug, **kwargs
            )
            dict_area[nn1], dict_keye[nn1], dict_var[nn1] = areacell, keyerror, tab
    keyerror = add_up_errors([dict_keye[ii] for ii in list(dict_keye.keys())])
    if keyerror is None:
        list_var = sorted(list(dict_var.keys()))
        if len(list_var) > 1:
            for ii in list(range(2)):
                for var in list_var[1:]:
                    dict_var[list_var[0]], dict_var[var], keyerror =\
                        CheckTime(dict_var[list_var[0]], dict_var[var], metric_name=metric, **kwargs)
                    if keyerror is not None:
                        break
    if keyerror is not None:
        tab, areacell = None, None
    else:
        tab, keyerror = MyDerive(kwargs[interpreter], variable, dict_var)
        areacell = dict_area[list(dict_area.keys())[0]]
    return tab, areacell, keyerror


def Read_mask_area(tab, name_data, file_data, type_data, region, file_area='', name_area='', file_mask='', name_mask='',
                   maskland=False, maskocean=False, debug=False, **kwargs):
    tab_out = copy.copy(tab)
    keyerror1, keyerror2 = None, None
    # Read areacell
    if file_area:
        areacell = ReadAreaSelectRegion(file_area, areaname=name_area, box=region, **kwargs)
    else:
        areacell = ReadAreaSelectRegion(file_data, areaname=name_area, box=region, **kwargs)
    if areacell is not None and tab.getGrid().shape != areacell.getGrid().shape:
        areacell = None
    if debug is True:
        if areacell is not None:
            dict_debug = {
                'axes1': '(' + type_data + ') ' + str([ax.id for ax in areacell.getAxisList()]),
                'shape1': '(' + type_data + ') ' + str(areacell.shape)
                }
            EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadAreaSelectRegion', 20, **dict_debug)
        else:
            dict_debug = {'line1': 'areacell is None '}
            EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadAreaSelectRegion', 20, **dict_debug)
    # Read landmask
    lvari = [
        "latent_heatflux", "lhf", "lwr", "meridional_wind_stress", "msla", "net_heating",
        "net_longwave_heatflux_downwards", "net_shortwave_heatflux_downwards", "net_surface_heatflux_downwards",
        "netflux", "sea_surface_height", "sea_surface_temperature", "sensible_heatflux", "shf", "sla", "sohefldo",
        "sometauy", "sossheig", "sosstsst", "sozotaux", "ssh", "sshg", "sst", "swr", "tauuo", "tauvo", "taux",
        "tauy", "thf", "thflx", "tmpsf", "tos", "zonal_wind_stress", "zos"
        ]
    if (name_data.lower() in lvari and "_Amon_" not in file_data) or \
        (name_data.lower() in ["pr", "slp"] and "_Omon_" in file_data):
        landmask = None
    elif file_mask:
        landmask = ReadLandmaskSelectRegion(tab, file_mask, landmaskname=name_mask, box=region, **kwargs)
    else:
        landmask = ReadLandmaskSelectRegion(tab, file_data, landmaskname=name_mask, box=region, **kwargs)
    if debug is True:
        if landmask is not None:
            dict_debug = {
                'axes1': '(' + type_data + ') ' + str([ax.id for ax in landmask.getAxisList()]),
                'shape1': '(' + type_data + ') ' + str(landmask.shape)
                }
            EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadLandmaskSelectRegion', 20, **dict_debug)
        else:
            dict_debug = {'line1': 'landmask is None '}
            EnsoErrorsWarnings.debug_mode('\033[93m', 'after ReadLandmaskSelectRegion', 20, **dict_debug)
    # Apply landmask
    if landmask is not None:
        tab_out, keyerror1 = ApplyLandmask(tab_out, landmask, maskland=maskland, maskocean=maskocean)
        if keyerror1 is None:
            if areacell is None:
                areacell = ArrayOnes(landmask, id='areacell')
            areacell, keyerror2 = ApplyLandmaskToArea(areacell, landmask, maskland=maskland, maskocean=maskocean)
    # When neither file_area nor landmask provided a valid areacell, synthesise
    # cosine-latitude weights so downstream Average* functions never receive
    # areacell=None (which would trigger silent None returns).
    if areacell is None:
        areacell = _make_coslat_areacell(tab_out)
        if debug is True:
            dict_debug = {'line1': 'areacell synthesised from cosine-latitude (no areacella/landmask available)'}
            EnsoErrorsWarnings.debug_mode('\033[93m', 'after areacell fallback', 20, **dict_debug)
    if keyerror1 is not None or keyerror2 is not None:
        keyerror = add_up_errors([keyerror1, keyerror2])
    else:
        keyerror = None
    return tab_out, areacell, keyerror


def SlabOcean(tab1, tab2, month1, month2, events, frequency=None, tmin=0.1, debug=False):
    """
    #################################################################################
    Description:
    Compute a simple slab-ocean estimate by integrating total heat-flux
    anomalies over time.

    Based on:
    Bayr, T., C. Wengel, M. Latif, D. Dommenget, J. Lübbecke, W. Park (2018)
    Error compensation of ENSO atmospheric feedbacks in climate models and its
    influence on simulated ENSO dynamics. Climate Dynamics.
    doi:10.1007/s00382-018-4575-7

    This function preserves the legacy ENSO_metrics interface while using the
    modern numpy/compatibility-layer masked-array pathway.
    #################################################################################

    :param tab1: CDATVariable or masked-array-like
        Sea-surface-temperature anomaly field, usually SSTA, with metadata and
        axes preserved through the compatibility layer.

    :param tab2: CDATVariable or masked-array-like
        Total heat-flux anomaly field, usually THFA, with metadata and axes
        preserved through the compatibility layer.

    :param month1: string
        First month of integration, for example ``"JUN"``.

    :param month2: string
        Last month of integration, for example ``"DEC"``.

    :param events: list of integer
        Years considered as ENSO events to be selected.

    :param frequency: string, optional
        Time frequency of the datasets, for example ``"monthly"``.
        default value = None

    :param tmin: float, optional
        Minimum temperature threshold for events. This is useful when
        Hovmöller diagnostics are provided and the minimum ENSO temperature
        is not reached everywhere.
        default value = 0.1

    :param debug: boolean, optional
        If True, print diagnostic information during the calculation.
        default value = False

    :return dSST, dSSTthf, dSSToce: CDATVariable or masked-array-like
        ``dSST`` is the normalized cumulative SST change from 0 to 1, in C/C;
        ``dSSTthf`` is the normalized cumulative heat-flux-driven SST change,
        in C/C; and ``dSSToce`` is the normalized cumulative SST change driven
        by anomalous ocean circulation, in C/C.
    """
    if debug is True:
        EnsoErrorsWarnings.debug_mode('\033[93m', "EnsoUvcdatToolsLib SlabOcean", 20)
    # months and associated position
    list_months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
    mm1, mm2 = 0, 0
    if month1 in list_months and month2 in list_months:
        mm1 = list_months.index(month1)
        mm2 = list_months.index(month2)
    else:
        list_strings = ["ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": month"]
        if month1 not in list_months:
            list_strings.append(str().ljust(5) + "unknown month1 : " + str(month1))
        if month2 not in list_months:
            list_strings.append(str().ljust(5) + "unknown month2 : " + str(month2))
        EnsoErrorsWarnings.my_error(list_strings)
    # sea water constants
    cp = 4000   # J/(kg * K) (specific heat capacity at constant pressure of sea water)
    rho = 1024  # kg/m3      (average density of sea water)
    H = 50      # m          (depth of the slab ocean)
    fraction = 60 * 60 * 24 * 30.42 / (cp * rho * H)  # W/m2 to C
    # selecting events
    sstA = Event_selection(tab1, frequency, nbr_years_window=2, list_event_years=events)
    thfA = Event_selection(tab2, frequency, nbr_years_window=2, list_event_years=events)
    if debug is True:
        dict_debug = {
            'axes1': '(sst) ' + str([ax.id for ax in sstA.getAxisList()]),
            'axes2': '(thf) ' + str([ax.id for ax in thfA.getAxisList()]),
            'shape1': '(sst) ' + str(sstA.shape), 'shape2': '(thf) ' + str(thfA.shape)
            }
        EnsoErrorsWarnings.debug_mode('\033[93m', 'after Event_selection', 25, **dict_debug)
    # cumulative anomalies
    myshape = [len(events), mm2-mm1+1] + [ss for ss in tab1.shape[1:]]
    dSST = MV2zeros(myshape)
    dSSTthf = MV2zeros(myshape)
    for ii in range(mm1, mm2):
        dSST[:, ii - mm1 + 1] = dSST[:, ii - mm1] + sstA[:, ii + 1] - sstA[:, ii]
        dSSTthf[:, ii - mm1 + 1] = dSSTthf[:, ii - mm1] + thfA[:, ii + 1]
    if debug is True:
        dict_debug = {'shape1': '(dSST) ' + str(dSST.shape), 'shape2': '(dSSTthf) ' + str(dSSTthf.shape)}
        EnsoErrorsWarnings.debug_mode('\033[93m', 'after cumulative_anomalies', 25, **dict_debug)
    # normalized heat flux-driven SST change
    dt = MV2zeros(dSSTthf.shape)
    dt = dt.reorder('10')
    dt[:] = dSST[:, -1]
    dt = dt.reorder('10')
    dt = MV2masked_where(abs(dt) < tmin, dt)
    dSSTthf[:] = fraction * dSSTthf[:] / dt
    # normalized SST change
    dSST[:] = dSST[:] / dt
    # normalized SST change by an anomalous ocean circulation
    dSSToce = dSST - dSSTthf
    # averaging across events
    dSST = MV2average(dSST, axis=0)
    dSSTthf = MV2average(dSSTthf, axis=0)
    dSSToce = MV2average(dSSToce, axis=0)
    # axes
    axes = [create_axis(MV2array(list(range(12-len(dSST), 12))), id='months')]
    if debug is True:
        dict_debug = {
            'axes1': 'axes ' + str(axes[0]), 'axes2': 'axes[:] ' + str(axes[0][:]),
            'shape1': '(dSST) ' + str(dSST.shape), 'shape2': '(dSSTthf) ' + str(dSSTthf.shape),
            'shape3': '(dSSToce) ' + str(dSSToce.shape)
            }
        EnsoErrorsWarnings.debug_mode('\033[93m', 'after mean dSST', 25, **dict_debug)
    if len(tab1.shape) > 1:
        axes = axes + tab1.getAxisList()[1:]
    dSST.setAxisList(axes)
    dSSTthf.setAxisList(axes)
    dSSToce.setAxisList(axes)
    if debug is True:
        dict_debug = {
            'axes1': '(dSST) ' + str([ax.id for ax in dSST.getAxisList()]),
            'axes2': '(dSSTthf) ' + str([ax.id for ax in dSSTthf.getAxisList()]),
            'axes3': '(dSSToce) ' + str([ax.id for ax in dSSToce.getAxisList()]),
            'shape1': '(dSST) ' + str(dSST.shape), 'shape2': '(dSSTthf) ' + str(dSSTthf.shape),
            'shape3': '(dSSToce) ' + str(dSSToce.shape)
            }
        EnsoErrorsWarnings.debug_mode('\033[93m', 'output', 25, **dict_debug)
    return dSST, dSSTthf, dSSToce


def TimeAnomaliesLinearRegressionAndNonlinearity(tab2, tab1, return_stderr=True):
    """
    #################################################################################
    Description:
    Compute linear regressions between interannual anomalies of two input fields
    after horizontal averaging and annual-cycle removal.

    The regression of ``tab2`` on ``tab1`` is computed for all values of
    ``tab1``, positive values of ``tab1``, and negative values of ``tab1``.
    This function preserves the legacy ENSO_metrics interface while using the
    modern compatibility-layer regression pathway.
    #################################################################################

    :param tab2: CDATVariable or masked-array-like
        Dependent variable.

    :param tab1: CDATVariable or masked-array-like
        Independent variable used to define all, positive, and negative
        regression branches.

    :param return_stderr: boolean, optional
        If True, return the unadjusted standard error of the regression slope.
        default value = True

    :return lr, lrpos, lrneg, keyerror:
        Regression results for all, positive, and negative values of ``tab1``,
        plus any accumulated keyerror message.
    """
    # horizontal average
    tab1, keyerror1 = dict_average['horizontal'](tab1)
    tab2, keyerror2 = dict_average['horizontal'](tab2)
    if keyerror1 is not None or keyerror2 is not None:
        lr, lrpos, lrneg = None, None, None
        keyerror = add_up_errors([keyerror1, keyerror2])
    else:
        keyerror = None
        # removes annual cycle (anomalies with respect to the annual cycle)
        tab1 = ComputeInterannualAnomalies(tab1)
        tab2 = ComputeInterannualAnomalies(tab2)
        # computes linear regression of tab2 on tab1 for all values of tab1, for values of tab1>=0,
        # for values of tab1<=0
        lr, lrpos, lrneg = LinearRegressionAndNonlinearity(tab2, tab1, return_stderr=return_stderr)
    return lr, lrpos, lrneg, keyerror


def TimeAnomaliesStd(tab):
    """
    #################################################################################
    Description:
    Compute the standard deviation of a spatially averaged anomaly time series.

    This function first computes a horizontal average, removes the annual cycle,
    and then calculates the temporal standard deviation. It preserves the legacy
    ENSO_metrics interface while using the modern numpy/scipy-based utilities and
    CDAT-like compatibility objects provided by the refactored xarray workflow.
    #################################################################################

    :param tab: CDATVariable or masked-array-like
        Input variable with latitude/longitude axes and time metadata.

    :return std, keyerror:
        Standard deviation of the spatially averaged time anomalies after the
        annual cycle has been removed, plus any accumulated keyerror message.
    """
    # horizontal average
    tab, keyerror = dict_average['horizontal'](tab)
    if keyerror is not None:
        std = None
    else:
        # computes standard deviation
        std = float(GENUTILstd(tab, weights=None, axis=0, centered=1, biased=1))
    return std, keyerror


def TsToMap(tab, map_ref):
    """
    #################################################################################
    Description:
    Expand a one-dimensional time series to match the shape and metadata of a
    reference map or field.

    The time-series values are inserted or broadcast into the spatial structure
    of ``map_ref`` so the output has the same dimensional layout as the
    reference field. This preserves the legacy ENSO_metrics interface while
    using CDAT-like compatibility objects provided by the refactored xarray
    workflow.
    #################################################################################

    :param tab: CDATVariable or masked-array-like
        One-dimensional input time series.

    :param map_ref: CDATVariable or masked-array-like
        Reference field whose shape, axes, and metadata are used to construct
        the output map.

    :return map_out: CDATVariable or masked-array-like
        Field with ``tab`` values expanded to the shape of ``map_ref``.
    """
    if len(map_ref.shape) > 6:
        list_strings = [
            "ERROR" + EnsoErrorsWarnings.message_formating(INSPECTstack()) + ": too many dimensions",
            str().ljust(5) + "map_ref.shape = " + str(map_ref.shape)]
        EnsoErrorsWarnings.my_error(list_strings)
    map_out = MV2zeros(map_ref.shape)
    map_out = create_variable(
        map_out,
        axes=map_ref.getAxisList(),
        grid=map_ref.getGrid(),
        mask=map_ref.mask,
        attributes=map_ref.attributes,
        id=tab.id
    )
    initorder = map_out.getOrder()
    map_out = map_out.reorder('...t')
    if len(map_ref.shape) == 2:
        map_out[:] = tab
    elif len(map_ref.shape) == 3:
        map_out[:, :] = tab
    elif len(map_ref.shape) == 4:
        map_out[:, :, :] = tab
    elif len(map_ref.shape) == 5:
        map_out[:, :, :, :] = tab
    else:
        map_out[:, :, :, :, :] = tab
    map_out = map_out.reorder(initorder)
    return map_out


def TwoVarRegrid(model, obs, info, region=None, model_orand_obs=0, newgrid=None, **keyarg):
    """
    #################################################################################
    Description:
    Regrid ``model``, ``obs``, or both so that the two fields are on a common
    horizontal grid before metric calculation.

    This function preserves the legacy ENSO_metrics calling interface, but the
    active regridding backend is now the xarray/xESMF/ESMF-compatible pathway
    implemented through ``Regrid()`` and the ENSO_metrics compatibility layer.
    #################################################################################

    :param model: CDATVariable or masked-array-like
        Model field with valid latitude/longitude axes and a rectilinear grid.

    :param obs: CDATVariable or masked-array-like
        Observational field with valid latitude/longitude axes and a
        rectilinear grid.

    :param info: string
        Description string updated to record what regridding operation was
        applied.

    :param region: string, optional
        Region name used when constructing a target grid from ``newgrid_name``.
        The region must be defined in ``EnsoCollectionsLib.ReferenceRegions``.
        In this workflow, ``ReferenceRegions`` is expected to use 0–360
        longitude bounds.

    :param model_orand_obs: integer, optional
        Controls which field is regridded:
            0: regrid model data onto the observations grid
            1: regrid observations data onto the model grid
            2: regrid both model and observations data onto ``newgrid`` or onto
               a grid constructed by ``Regrid()`` from ``newgrid_name`` and
               ``region``
        default value = 0

    :param newgrid: grid-like object, optional
        Target grid used when ``model_orand_obs=2``. If ``newgrid`` is ``None``
        or a string, ``Regrid()`` constructs the target grid from
        ``newgrid_name`` and ``region``.

    usual kwargs:
    :param newgrid_name: string, optional
        Name used by ``Regrid()`` to construct the destination grid when
        ``newgrid`` is not provided explicitly. The name should specify a
        supported grid type and resolution, for example ``"generic_1x1deg"``.
        default value = ``"generic_1x1deg"``

    :param regridder: string, optional
        Regridding backend name. The modern supported/default value is
        ``"xesmf"``.

    :param regridTool: string, optional
        Regridding tool identifier retained for compatibility with existing
        ENSO_metrics configuration dictionaries. The modern expected value is
        ``"esmf"``.

    :param regridMethod: string, optional
        Regridding method. The modern default is ``"bilinear"``. Legacy
        ``"linear"`` is normalized internally to ``"bilinear"`` by ``Regrid()``
        where supported.

    :return: model, obs, info
        Model and observation fields on a common grid, plus an updated
        description of the regridding operation.
    """
    debug = keyarg.pop('debug', False)
    known_args = {'missing', 'order', 'mask', 'newgrid_name', 'regridder', 'regridTool', 'regridMethod'}
    extra_args = set(keyarg) - known_args
    if extra_args:
        EnsoErrorsWarnings.unknown_key_arg(extra_args, INSPECTstack())
    grid_obs = obs.getGrid()
    grid_model = model.getGrid()
    # select case:
    if model_orand_obs == 0:
        model = Regrid(model, grid_obs, **keyarg)
        info = info + ', model regridded to observations'
    elif model_orand_obs == 1:
        obs = Regrid(obs, grid_model, **keyarg)
        info = info + ', observations regridded to model'
    elif model_orand_obs == 2:
        if debug:
            if newgrid is not None and not isinstance(newgrid, str):
                print("DEBUG: [TwoVarRegrid] newgrid lat:", newgrid.getLatitude()[:])
                print("DEBUG: [TwoVarRegrid] newgrid lon:", newgrid.getLongitude()[:])
            else:
                print(
                    "DEBUG: [TwoVarRegrid] newgrid will be constructed by Regrid() "
                    f"from newgrid_name={keyarg.get('newgrid_name')!r}, region={region!r}"
                )
            print("DEBUG: [TwoVarRegrid] model original shape:", model.shape, "obs original shape:", obs.shape)
        model = Regrid(model, newgrid, region=region, **keyarg)
        obs = Regrid(obs, newgrid, region=region, **keyarg)
        if debug:
            print("DEBUG: [TwoVarRegrid] model regridded shape:", model.shape, "obs regridded shape:", obs.shape)
        try:
            grid_name = newgrid.id
        except Exception:
            try:
                grid_name = newgrid.name
            except Exception:
                try:
                    grid_name = keyarg['newgrid_name']
                except Exception:
                    grid_name = 'newgrid'
        info = info + ', observations and model regridded to ' + str(grid_name)
    else:
        raise ValueError(
            f"TwoVarRegrid: unknown model_orand_obs={model_orand_obs!r}; "
            "expected 0, 1, or 2."
        )
    # Validate that regridding produced consistent spatial shapes.
    # Only documented regridding modes 0, 1, and 2 are supported.
    if model_orand_obs in (0, 1, 2):
        m_spatial = model.shape[-2:] if model.ndim >= 2 else model.shape
        o_spatial = obs.shape[-2:] if obs.ndim >= 2 else obs.shape
        if m_spatial != o_spatial:
            raise ValueError(
                f"TwoVarRegrid: spatial shape mismatch after regridding "
                f"(model_orand_obs={model_orand_obs}): "
                f"model {model.shape} vs obs {obs.shape}. "
                "Check newgrid or regridding parameters."
            )
    if model.shape == obs.shape:
        if model.mask.shape != ():
            mask = model.mask
            if obs.mask.shape != ():
                mask = MV2where(obs.mask, obs.mask, mask)
        else:
            if obs.mask.shape != ():
                mask = obs.mask
            else:
                mask = MV2where(MV2zeros(model.shape)==0, False, True)
        model = MV2masked_where(mask, model)
        obs = MV2masked_where(mask, obs)
    else:
        if obs[0].mask.shape != ():
            tab = MV2zeros(model.shape)
            for tt in list(range(len(tab))):
                tab[tt] = MV2masked_where(obs[0].mask, tab[tt])
            model = MV2masked_where(tab.mask, model)
        if model[0].mask.shape != ():
            tab = MV2zeros(obs.shape)
            for tt in list(range(len(tab))):
                tab[tt] = MV2masked_where(model[0].mask, tab[tt])
            obs = MV2masked_where(tab.mask, obs)
    return model, obs, info
