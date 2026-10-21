#!/usr/bin/env python3
"""Report what this machine can run: packages, raw data, derived data.

    python scripts/check_environment.py
"""
from __future__ import annotations
import importlib, os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PINNED = {"numpy": "2.3.5", "pandas": "2.3.3", "scipy": "1.16.3", "statsmodels": "0.14.5",
          "matplotlib": "3.10.6", "wfdb": "4.3.1"}
OPTIONAL = {"yaml": "configs/analysis.yaml reader (GI stages)", "seaborn": "some figures",
            "pytest": "tests/", "pyarrow": "parquet annotation cache"}


def version(name):
    try:
        return importlib.import_module(name).__version__
    except ImportError:
        return None
    except AttributeError:
        return "?"


print(f"python {sys.version.split()[0]}")
ok = True
for name, pinned in PINNED.items():
    v = version(name)
    state = "MISSING" if v is None else ("ok" if v == pinned else f"differs (pinned {pinned})")
    print(f"  {name:<12} {v or '-':<10} {state}")
    ok &= v is not None
for name, why in OPTIONAL.items():
    v = version(name)
    print(f"  {name:<12} {v or '-':<10} {'ok' if v else 'optional, missing'}  ({why})")

raw = Path(os.environ.get("EHG_RAW_DIR", ROOT / "data" / "raw")).expanduser()
has_raw = (raw / "RECORDS").exists()
print(f"\nraw data     : {raw}  ->  {'FOUND' if has_raw else 'not found (only needed for the raw stages)'}")
needed = ["work/windows/window_manifest_60s.csv", "work/features/window_features_60s_EHG9.csv",
          "work/primary/event_unit_analysis_60s_EHG9.csv", "work/gi_features_60s_EHG9.csv"]
has_derived = all((ROOT / p).exists() for p in needed)
print(f"derived data : {'FOUND' if has_derived else 'not restored  ->  python scripts/restore_derived.py'}")
print("\nrunnable now :")
print("  tests, conditional analyses, S7/S8 tables      :", "yes" if version("scipy") and version("pandas") else "no")
print("  mixed-model stages (need statsmodels + derived) :", "yes" if version("statsmodels") and has_derived else "no")
print("  raw-data stages (need wfdb + raw data)          :", "yes" if version("wfdb") and has_raw else "no")
sys.exit(0 if ok else 1)
