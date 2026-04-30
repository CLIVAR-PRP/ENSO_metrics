"""Quick symbol consistency check across the refactored package."""
import re, sys

# 1. create_* factory functions reachable from EnsoUvcdatToolsLib
from EnsoMetrics.EnsoUvcdatToolsLib import (
    create_axis, create_variable,
    create_uniform_lat_axis, create_uniform_lon_axis,
    create_rect_grid, open_file)
print("OK  create_* factory functions present in EnsoUvcdatToolsLib")

# 2. XarrayCompat.__all__ is self-consistent
from EnsoMetrics import XarrayCompat as xc
bad = [n for n in xc.__all__ if not hasattr(xc, n)]
print("XarrayCompat __all__ issues:", bad or "none")

# 3. Every name EnsoMetricsLib imports from EnsoUvcdatToolsLib exists
with open("lib/EnsoMetricsLib.py") as f:
    src = f.read()
block = re.search(
    r"from .EnsoUvcdatToolsLib import(.*?)(?=\nfrom |\nimport |\ndef |\nclass )",
    src, re.DOTALL)
if block:
    raw = block.group(1).replace("\\", " ").replace(",", " ")
    names = [t for t in raw.split() if re.match(r"^[A-Za-z_]\w*$", t)]
    import EnsoMetrics.EnsoUvcdatToolsLib as uv
    missing = [n for n in names if not hasattr(uv, n)]
    print(f"EnsoMetricsLib needs {len(names)} names; missing: {missing or 'none'}")

# 4. Every name EnsoComputeMetricsLib imports from EnsoMetricsLib exists
with open("lib/EnsoComputeMetricsLib.py") as f:
    src2 = f.read()
block2 = re.search(
    r"from .EnsoMetricsLib import(.*?)(?=\nfrom |\nimport |\ndef |\nclass |\Z)",
    src2, re.DOTALL)
if block2:
    raw2 = block2.group(1).replace("\\", " ").replace(",", " ")
    names2 = [t for t in raw2.split() if re.match(r"^[A-Za-z_]\w*$", t)]
    import EnsoMetrics.EnsoMetricsLib as em
    missing2 = [n for n in names2 if not hasattr(em, n)]
    print(f"EnsoComputeMetricsLib needs {len(names2)} names from EnsoMetricsLib; missing: {missing2 or 'none'}")

print("DONE")
