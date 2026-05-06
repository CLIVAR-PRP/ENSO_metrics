from setuptools import setup
import glob
import re

# Read the version from lib/version.py — that file is the single source of
# truth and is tracked by git.  setup.py no longer writes it.
with open("lib/version.py") as _vf:
    _src = _vf.read()
Version = re.search(r"__version__\s*=\s*'([^']+)'", _src).group(1)

# data_files = (
#               ('share/EnsoMetrics', ['share/EnsoMetrics/basin_generic_1x1deg.nc']),
#              )

data_files = (
    (
        'share/EnsoMetrics',
        [
            'share/EnsoMetrics/basin_generic_1x1deg.nc',
            'share/EnsoMetrics/cmip5_historical_ENSO_perf_v20200427_allModels_allRuns.json',
            'share/EnsoMetrics/cmip5_historical_ENSO_proc_v20200427_allModels_allRuns.json',
            'share/EnsoMetrics/cmip5_historical_ENSO_tel_v20200427_allModels_allRuns.json',
            'share/EnsoMetrics/cmip6_historical_ENSO_perf_v20200427_allModels_allRuns.json',
            'share/EnsoMetrics/cmip6_historical_ENSO_proc_v20200427_allModels_allRuns.json',
            'share/EnsoMetrics/cmip6_historical_ENSO_tel_v20200427_allModels_allRuns.json',
            'share/EnsoMetrics/obs2obs_historical_ENSO_perf_v20201231_allObservations.json',
            'share/EnsoMetrics/obs2obs_historical_ENSO_proc_v20201231_allObservations.json',
            'share/EnsoMetrics/obs2obs_historical_ENSO_tel_v20201231_allObservations.json'],
    ),
)

setup(name="EnsoMetrics",
      author="Eric Guilyardi",
      version=Version,
      description = "Library for ENSO Metrics",
      url="https://github.com/CLIVAR-PRP/ENSO_metrics",
      packages=['EnsoMetrics', 'EnsoPlots'],
      package_dir={'EnsoMetrics': 'lib', 'EnsoPlots': 'plots'},
      scripts=glob.glob("scripts/*.py"),
      data_files=data_files,
      install_requires=[
          'numpy',
          'scipy',
          'packaging',
          'xarray',
          'xcdat',
          'cftime',
          'regionmask',
          'xesmf',
      ])
