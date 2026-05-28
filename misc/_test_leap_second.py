"""Test that leap-second (second=60) timestamps are handled correctly."""
import cftime, numpy as np, xarray as xr, tempfile, os, netCDF4 as nc4

# Build a NetCDF file with second=60 encoded as a numeric offset
ds = xr.Dataset({'sst': xr.DataArray(
    np.ones((2,)), dims=['time'],
    coords={'time': [cftime.DatetimeNoLeap(2015, 12, 31, 23, 59, 59),
                     cftime.DatetimeNoLeap(2015, 12, 31, 23, 59, 59)]})})
with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as f:
    tmp = f.name
ds.to_netcdf(tmp)

# Nudge the last time step by +1 second to create second=60
with nc4.Dataset(tmp, 'r+') as nf:
    t = nf.variables['time']
    t[1] = t[1] + 1.0 / 86400

# Open via our wrapper and call TimeBounds (exercises _clamp_leap_seconds +
# asComponentTime paths)
from EnsoMetrics.EnsoUvcdatToolsLib import open_file, TimeBounds
fh = open_file(tmp)
sst = fh('sst')
t0, t1 = TimeBounds(sst)
print('TimeBounds:', t0, '..', t1)
assert '60' not in t1, f'leap second NOT clamped: {t1}'
print('OK - leap second clamped correctly')
os.unlink(tmp)
