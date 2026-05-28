#!/usr/bin/env python
"""End-to-end verification: package works without cdms2, with xcdat/numpy."""
import sys, os, tempfile
import numpy as np
import numpy.ma as ma
import cftime
import xarray as xr

PASS, FAIL = [], []
def ok(msg):   PASS.append(msg); print(f"OK   {msg}")
def fail(msg): FAIL.append(msg); print(f"FAIL {msg}")

# ── 1. Environment ────────────────────────────────────────────────────────
print("\n=== 1. Environment ===")
try:
    import cdms2
    fail("cdms2 is importable — it must NOT be present")
except ImportError:
    ok("cdms2 absent")

import xcdat;      ok(f"xcdat {xcdat.__version__}")
import numpy;      ok(f"numpy {numpy.__version__}")
import scipy;      ok(f"scipy {scipy.__version__}")
import cftime as _cf; ok(f"cftime {_cf.__version__}")
try:
    import regionmask; ok(f"regionmask {regionmask.__version__}")
except ImportError:
    print("WARN regionmask not installed (land mask fallback active)")
try:
    import xesmf; ok(f"xesmf {xesmf.__version__}")
except (ImportError, OSError) as _e:
    print(f"WARN xesmf unavailable ({_e.__class__.__name__}: {_e}) — REGRID2 will raise on use")

# ── 2. Package imports ────────────────────────────────────────────────────
print("\n=== 2. Package imports ===")
try:
    from EnsoMetrics.XarrayCompat import (
        CDATVariable, _Axis, _Grid,
        create_axis, create_uniform_lat_axis, create_uniform_lon_axis,
        create_rect_grid, create_variable, da_to_cdat, cdat_to_da)
    ok("XarrayCompat")
except Exception as e:
    fail(f"XarrayCompat: {e}"); sys.exit(1)

try:
    from EnsoMetrics.EnsoUvcdatToolsLib import (
        create_variable, create_axis, open_file, CDTIMEcomptime,
        MV2add, MV2subtract, MV2multiply, MV2divide, MV2array, MV2zeros,
        MV2ones, MV2average, MV2sum, MV2compress, MV2concatenate,
        MV2where, MV2masked_where, MV2take, MV2maximum, MV2minimum,
        GENUTILrms, GENUTILstd, GENUTILcorrelation, GENUTILlinearregression,
        SCIPYsignal_detrend, SCIPYstats__skew, cdutil, sea_dict,
        _axis_to_int, TimeBounds)
    ok("EnsoUvcdatToolsLib (all public names)")
except Exception as e:
    fail(f"EnsoUvcdatToolsLib: {e}"); sys.exit(1)

try:
    from EnsoMetrics import EnsoMetricsLib
    ok("EnsoMetricsLib (27 k lines)")
except Exception as e:
    fail(f"EnsoMetricsLib: {e}")

try:
    from EnsoMetrics import EnsoComputeMetricsLib
    ok("EnsoComputeMetricsLib")
except Exception as e:
    fail(f"EnsoComputeMetricsLib: {e}")

# ── 3. CDATVariable ───────────────────────────────────────────────────────
print("\n=== 3. CDATVariable ===")
lat  = create_uniform_lat_axis(-89.5, 180, 1.0)
lon  = create_uniform_lon_axis(0.5,  360, 1.0)
time_vals = np.array([cftime.datetime(1980+i//12, i%12+1, 15) for i in range(24)],
                     dtype=object)
time_ax = _Axis("time", time_vals, units="days since 1980-01-01", axis_type="T")
grid = create_rect_grid(lat, lon)
data = ma.array(np.random.randn(24, 180, 360))
var  = create_variable(data, axes=[time_ax, lat, lon], grid=grid, id="sst")

assert var.shape == (24, 180, 360);      ok(f"shape {var.shape}")
assert var.getGrid().shape == (180,360); ok(f"grid {var.getGrid().shape}")
assert var.getTime() is not None;        ok("getTime()")
assert var.getLatitude() is not None;    ok("getLatitude()")
assert var.getOrder() == "tyx";          ok(f"getOrder() = {var.getOrder()!r}")

# slicing
sliced = var[0:2]
assert sliced.shape == (2,180,360);      ok(f"slice [0:2] = {sliced.shape}")

# callable selection
sub = var(latitude=(-30,30), longitude=(150,250))
assert sub.shape[1] < 180;              ok(f"var(lat,lon) subset = {sub.shape}")

# reorder
reordered = var.reorder("tyx")
assert reordered.shape == (24,180,360); ok(f"reorder('tyx') = {reordered.shape}")

ct = time_ax.asComponentTime()
assert len(ct) == 24;                   ok(f"asComponentTime() len={len(ct)}")

# ── 4. MV2 aliases ────────────────────────────────────────────────────────
print("\n=== 4. MV2 aliases ===")
a = MV2array(np.arange(12, dtype=float).reshape(3,4))
b = MV2ones((3,4))
assert MV2add(a, b).shape == (3,4);         ok("MV2add")
assert MV2subtract(a, b).shape == (3,4);    ok("MV2subtract")
assert MV2multiply(a, b).shape == (3,4);    ok("MV2multiply")
assert MV2divide(a, b+1).shape == (3,4);    ok("MV2divide")
assert float(MV2sum(a)) == 66.0;            ok("MV2sum")
assert float(MV2average(a)) == 5.5;         ok("MV2average")
assert float(MV2maximum(a)) == 11.0;        ok("MV2maximum")
assert float(MV2minimum(a)) == 0.0;         ok("MV2minimum")
c1 = MV2array([[1.,2.],[3.,4.]])
c2 = MV2array([[5.,6.],[7.,8.]])
assert MV2concatenate([c1,c2]).shape == (4,2); ok("MV2concatenate")
cond = np.array([True,False]*6)
assert MV2compress(cond, a.ravel()).shape == (6,); ok("MV2compress")

# ── 5. CDTIMEcomptime ────────────────────────────────────────────────────
print("\n=== 5. CDTIMEcomptime ===")
t = CDTIMEcomptime(2000, 6, 15, 12, 0, 0.0)
assert t.year == 2000 and t.month == 6; ok(f"{t}")

# ── 6. GENUTIL replacements ───────────────────────────────────────────────
print("\n=== 6. GENUTIL replacements ===")
x = ma.array(np.random.randn(12,4,5))
y = ma.array(np.random.randn(12,4,5))

r0 = GENUTILrms(x, y, axis=0)
assert r0.shape == (4,5);               ok(f"GENUTILrms(axis=0) -> {r0.shape}")

r_xy = GENUTILrms(x, y, axis="xy")
assert r_xy.shape == (12,);             ok(f"GENUTILrms(axis='xy') -> {r_xy.shape}")

r_w = GENUTILrms(x, y, weights="weighted", axis="xy")
assert r_w.shape == (12,);              ok(f"GENUTILrms(weights='weighted') -> {r_w.shape}")

r_t = GENUTILrms(x[:,0,0], y[:,0,0], axis="t")
ok(f"GENUTILrms(axis='t') = {float(r_t):.4f}")

r_y  = GENUTILrms(x[:,0,:], y[:,0,:], axis="x")
ok(f"GENUTILrms(axis='x') -> {r_y.shape if hasattr(r_y,'shape') else r_y}")

s = GENUTILstd(x, axis="xy")
assert s.shape == (12,);                ok(f"GENUTILstd(axis='xy') -> {s.shape}")

s0 = GENUTILstd(x, axis=0)
assert s0.shape == (4,5);               ok(f"GENUTILstd(axis=0) -> {s0.shape}")

corr = GENUTILcorrelation(x[:,0,0], y[:,0,0], axis=0)
ok(f"GENUTILcorrelation = {float(corr):.4f}")

slope_int, stderr = GENUTILlinearregression(x[:,0,0])
ok(f"GENUTILlinearregression slope={slope_int[0,0]:.4f} stderr={stderr[0,0]:.4f}")

# ── 7. cdutil.averager (all axis forms) ──────────────────────────────────
print("\n=== 7. cdutil.averager ===")
lat2 = create_uniform_lat_axis(-45, 10, 9.)
lon2 = create_uniform_lon_axis(120, 20, 9.)
var2d = create_variable(
    ma.array(np.random.randn(24,10,20)),
    axes=[time_ax, lat2, lon2],
    grid=create_rect_grid(lat2, lon2), id="ts")

for ax_spec, expected_shape in [
    ("t",  (10,20)),
    ("xy", (24,)),
    ("y",  (24,20)),
    ("x",  (24,10)),
    ("1",  (24,20)),   # numeric: lat axis index
    ("2",  (24,10)),   # numeric: lon axis index
]:
    result = cdutil.averager(var2d, axis=ax_spec)
    assert result.shape == expected_shape, \
        f"axis={ax_spec!r}: expected {expected_shape}, got {result.shape}"
    ok(f"cdutil.averager(axis={ax_spec!r}) -> {result.shape}")

# cdutil.ANNUALCYCLE.departures
try:
    anom = cdutil.ANNUALCYCLE.departures(var2d)
    ok(f"cdutil.ANNUALCYCLE.departures() shape={anom.shape}")
except Exception as e:
    fail(f"ANNUALCYCLE.departures: {e}")

# ── 8. sea_dict ──────────────────────────────────────────────────────────
print("\n=== 8. sea_dict ===")
for k in ["JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC",
          "DJF","MAM","JJA","SON","NDJ","NDJF","DJFM"]:
    assert k in sea_dict, f"missing {k}"
ok(f"{len(sea_dict)} season keys present")

# ── 9. cdat_to_da / da_to_cdat round-trip ────────────────────────────────
print("\n=== 9. CDATVariable <-> xr.DataArray round-trip ===")
da = cdat_to_da(var2d)
assert isinstance(da, xr.DataArray);    ok(f"cdat_to_da -> DataArray {da.shape}")
back = da_to_cdat(da, varname="ts")
assert back.shape == var2d.shape;       ok(f"da_to_cdat -> CDATVariable {back.shape}")
assert back.getLatitude() is not None;  ok("lat axis preserved through round-trip")

# ── 10. open_file NetCDF I/O ──────────────────────────────────────────────
print("\n=== 10. open_file (NetCDF I/O) ===")
with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as f:
    tmp = f.name
ds_out = xr.Dataset({"sst": xr.DataArray(
    np.random.randn(12,4,5).astype("float32"),
    dims=["time","lat","lon"],
    coords={"time": [cftime.datetime(2000,m,15) for m in range(1,13)],
            "lat":  np.linspace(-15,15,4),
            "lon":  np.linspace(150,270,5)})})
ds_out.to_netcdf(tmp)

fh = open_file(tmp)
sst = fh("sst")
assert sst.shape == (12,4,5);           ok(f"read shape={sst.shape}")
sst_sub = fh("sst", latitude=(-10,10), longitude=(160,250))
ok(f"spatial subset shape={sst_sub.shape}")
t0, t1 = TimeBounds(sst)
ok(f"TimeBounds = {t0!r} .. {t1!r}")
fh.close()
os.unlink(tmp)

# ── 11. _axis_to_int edge cases ──────────────────────────────────────────
print("\n=== 11. _axis_to_int ===")
arr = ma.zeros((3,4,5))
cases = [
    (None,  None),
    (1,     1),
    ("xy",  (1,2)),
    ("yx",  (1,2)),
    ("x",   2),
    ("y",   1),
    ("t",   0),
    ("01",  (0,1)),
    ("10",  (1,0)),
    ("12",  (1,2)),
]
for inp, expected in cases:
    result = _axis_to_int(arr, inp)
    assert result == expected, f"_axis_to_int(arr, {inp!r}) = {result!r}, expected {expected!r}"
ok(f"all {len(cases)} axis-string cases correct")

# ── Summary ───────────────────────────────────────────────────────────────
print("\n" + "="*60)
if FAIL:
    print(f"RESULT: {len(FAIL)} FAILURES, {len(PASS)} passed")
    for f in FAIL:
        print(f"  FAIL  {f}")
    sys.exit(1)
else:
    print(f"RESULT: ALL {len(PASS)} CHECKS PASSED")
    print("Package works without cdms2 (xcdat {xcdat.__version__} / numpy {numpy.__version__})")
    print("="*60)
