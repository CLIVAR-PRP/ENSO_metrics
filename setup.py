
from distutils.core import setup
import subprocess
import glob

Version="0.1"

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
    commit = p.stdout.readlines()[0].split()[1]
except:
    commit = ""
f = open("lib/version.py", "w")
print("__version__ = '%s'" % Version, file=f)
print("__git_tag_describe__ = '%s'" % descr, file=f)
print("__git_sha1__ = '%s'" % commit, file=f)
f.close()

data_files = (
              ('share/EnsoMetrics', ['share/EnsoMetrics/basin_generic_1x1deg.nc']),
             )

setup (name="EnsoMetrics",
       author="Eric Guilyardi",
       version=Version,
       description = "Library for ENSO Metrics",
       url="https://github.com/eguil/ENSO_metrics",
       packages=['EnsoMetrics'],
       package_dir={'EnsoMetrics': 'lib'},
       scripts=glob.glob("scripts/*.py"),
       data_files=data_files,
      )

