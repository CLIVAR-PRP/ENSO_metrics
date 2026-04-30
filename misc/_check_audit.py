"""
Package-wide consistency audit for ENSO_metrics.
Run with:  python misc/_check_audit.py
"""
import os
import re
import ast

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

ACTIVE_PY = [
    "lib/XarrayCompat.py",
    "lib/EnsoUvcdatToolsLib.py",
    "lib/EnsoMetricsLib.py",
    "lib/EnsoComputeMetricsLib.py",
    "lib/EnsoCollectionsLib.py",
    "lib/EnsoToolsLib.py",
    "lib/EnsoErrorsWarnings.py",
    "lib/EnsoPlotLib.py",
    "lib/KeyArgLib.py",
    "lib/version.py",
    "lib/__init__.py",
    "plots/__init__.py",
    "plots/EnsoMetricPlot.py",
    "plots/EnsoPlotTemplate.py",
    "plots/EnsoPlotToolsLib.py",
    "scripts/driver_tools_lib.py",
    "scripts/driver_metric_collection.py",
    "scripts/enso_metrics_driver_example.py",
    "setup.py",
]

ACTIVE_OTHER = [
    "conda/meta.yaml",
]

OLD_API_NAMES = [
    "CDMS2open", "CDMS2createAxis", "CDMS2createVariable",
    "CDMS2createUniformLat", "CDMS2createUniformLon", "CDMS2createRect",
    "CDMS2setAutoBounds",
]

RETIRED_PKGS = ("cdms2", "cdutil", "genutil", "regrid2", "udunits2")

issues = []
ok = []


def check(path_rel):
    path = os.path.join(ROOT, path_rel)
    if not os.path.exists(path):
        issues.append(f"MISSING FILE: {path_rel}")
        return None, None
    with open(path) as f:
        src = f.read()
    return src, src.splitlines()


# -----------------------------------------------------------------------
# 1. Syntax check all Python files
# -----------------------------------------------------------------------
for rel in ACTIVE_PY:
    src, lines = check(rel)
    if src is None:
        continue
    try:
        ast.parse(src)
        ok.append(f"syntax OK: {rel}")
    except SyntaxError as e:
        issues.append(f"SYNTAX ERROR {rel}:{e.lineno}: {e.msg}")

# -----------------------------------------------------------------------
# 2. Live (non-docstring) retired CDAT imports
# -----------------------------------------------------------------------
for rel in ACTIVE_PY:
    src, lines = check(rel)
    if src is None or lines is None:
        continue
    for i, raw in enumerate(lines, 1):
        stripped = raw.strip()
        if stripped.startswith("#"):
            continue
        prev = lines[i - 2].strip() if i >= 2 else ""
        if "for more information" in prev or ">>>" in prev:
            continue
        for pkg in RETIRED_PKGS:
            if re.match(rf"^\s*(import {pkg}|from {pkg})\b", raw):
                issues.append(f"LIVE CDAT IMPORT {rel}:{i}: {stripped}")

# -----------------------------------------------------------------------
# 3. Module-level cftime import (XarrayCompat must have none)
# -----------------------------------------------------------------------
src, lines = check("lib/XarrayCompat.py")
if src and lines:
    for i, raw in enumerate(lines[:70], 1):
        stripped = raw.strip()
        if re.match(r"^import cftime|^from cftime", stripped):
            issues.append(f"MODULE-LEVEL cftime XarrayCompat.py:{i}: {stripped}")

# -----------------------------------------------------------------------
# 4. Old public API names (outside backup files)
# -----------------------------------------------------------------------
for rel in ACTIVE_PY:
    src, lines = check(rel)
    if src is None or lines is None:
        continue
    for i, raw in enumerate(lines, 1):
        if raw.strip().startswith("#"):
            continue
        for name in OLD_API_NAMES:
            if name in raw:
                issues.append(f"OLD API NAME {rel}:{i}: {raw.strip()}")

# -----------------------------------------------------------------------
# 5. lib/__init__.py health
# -----------------------------------------------------------------------
src, _ = check("lib/__init__.py")
if src:
    if "XarrayCompat" not in src:
        issues.append("lib/__init__.py: missing 'import XarrayCompat'")
    for bad in ("__git_tag_describe__", "__git_sha1__"):
        if bad in src:
            issues.append(f"lib/__init__.py: still imports removed symbol {bad}")
    ok.append("lib/__init__.py checked")

# -----------------------------------------------------------------------
# 6. setup.py version plumbing
# -----------------------------------------------------------------------
src, _ = check("setup.py")
if src:
    if "version=Version" not in src:
        issues.append("setup.py: does not use 'version=Version' — version may be hardcoded")
    if re.search(r"open\(.*version\.py.*['\"]w['\"]", src):
        issues.append("setup.py: still WRITES lib/version.py (should only read it)")
    ok.append("setup.py checked")

# -----------------------------------------------------------------------
# 7. Version consistency: lib/version.py vs conda/meta.yaml
# -----------------------------------------------------------------------
ver_lib = ver_conda = None
src, _ = check("lib/version.py")
if src:
    m = re.search(r"__version__\s*=\s*['\"]([^'\"]+)['\"]", src)
    ver_lib = m.group(1) if m else None
    if ver_lib is None:
        issues.append("lib/version.py: cannot parse __version__")

src_yaml, _ = check("conda/meta.yaml")
if src_yaml:
    # conda meta.yaml uses Jinja2: {% set version = "1.2.0" %}
    m = re.search(r'{%\s*set\s+version\s*=\s*["\']([^"\']+)["\']', src_yaml)
    ver_conda = m.group(1) if m else None
    if ver_conda is None:
        issues.append("conda/meta.yaml: cannot parse version")

if ver_lib and ver_conda and ver_lib != ver_conda:
    issues.append(f"VERSION MISMATCH: lib/version.py={ver_lib!r}  conda/meta.yaml={ver_conda!r}")
elif ver_lib and ver_conda:
    ok.append(f"Version consistent: {ver_lib}")

# -----------------------------------------------------------------------
# 8. scripts: open_file used (not old CDMS2open)
# -----------------------------------------------------------------------
for rel in ["scripts/driver_tools_lib.py"]:
    src, _ = check(rel)
    if src:
        if "CDMS2open" in src:
            issues.append(f"{rel}: still uses old CDMS2open")
        if "open_file" not in src:
            issues.append(f"{rel}: does not call open_file")

# -----------------------------------------------------------------------
# 9. conda/meta.yaml: no retired packages in run deps
# -----------------------------------------------------------------------
if src_yaml:
    for pkg in RETIRED_PKGS:
        # look only in 'run:' section
        m = re.search(r"run:(.*?)(?=\ntest:|\nbuild:|\nabout:|\Z)", src_yaml, re.DOTALL)
        if m and re.search(rf"\b{pkg}\b", m.group(1)):
            issues.append(f"conda/meta.yaml run deps: still lists retired package '{pkg}'")
    ok.append("conda/meta.yaml checked")

# -----------------------------------------------------------------------
# 10. build/ directory should not contain stale .py sources
# -----------------------------------------------------------------------
build_lib = os.path.join(ROOT, "build", "lib", "EnsoMetrics")
if os.path.isdir(build_lib):
    issues.append(
        f"STALE BUILD ARTIFACTS: {build_lib} exists — "
        "run 'python setup.py clean --all' or remove build/ to avoid stale imports"
    )

# -----------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------
print(f"\n{'='*60}")
if issues:
    print(f"ISSUES FOUND: {len(issues)}")
    print("="*60)
    for iss in issues:
        print(" ", iss)
else:
    print("ALL CONSISTENCY CHECKS PASSED")
print("="*60)
print(f"  {len(ok)} checks OK,  {len(issues)} issues")
