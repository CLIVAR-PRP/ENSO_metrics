from distutils.core import setup
import subprocess
import glob

Version = "1.1.3"

p = subprocess.Popen(
    ("git",
     "describe",
     "--tags"),
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
try:
    descr = p.stdout.readlines()[0].strip()
    Version = "-".join(descr.split("-")[:-2])
    if Version == "":
        Version = descr
except:
    descr = Version

p = subprocess.Popen(
    ("git",
     "log",
     "-n1",
     "--pretty=short"),
    stdin=subprocess.PIPE,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE)
try:
    commit = str(p.stdout.readlines()[0].split()[1])
except:
    commit = "unknown"
f = open("lib/version.py", "w")
print("__version__ = '%s'" % Version, file=f)
print("__git_tag_describe__ = '%s'" % descr, file=f)
print("__git_sha1__ = %s" % commit, file=f)
f.close()

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
      data_files=data_files, requires=['numpy'])
