#!/usr/bin/env python3
"""Unpack data/derived/ehg_derived_inputs.zip into work/ and verify every file's SHA-256.

    python scripts/restore_derived.py [path/to/ehg_derived_inputs.zip]

The zip is not stored in git history (it is ~35 MB): download it from the repository's release
page / Zenodo record named in README.md and place it in data/derived/, or pass its path.
"""
import hashlib, sys, zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# What the zip contains (relative to the repository root).
DERIVED_DIRS = ["work/annotations", "work/audit", "work/windowing", "work/windows", "work/features",
                "work/feature_qc", "work/primary"]
DERIVED_FILES = ["work/gi_features_60s_EHG9.csv", "work/cross_channel_gi_features.csv",
                 "work/invariance_window_results.csv", "work/gain_resolution_window_results.csv",
                 "work/resolution_ablation_features.csv"]


def main():
    zpath = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "data" / "derived" / "ehg_derived_inputs.zip"
    if not zpath.exists():
        sys.exit(f"{zpath} not found. Download ehg_derived_inputs.zip (see README, 'Getting the derived data').")
    manifest = ROOT / "data" / "derived" / "MANIFEST.sha256"
    want = {}
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            h, name = line.split("  ", 1)
            want[name] = h
        h = hashlib.sha256(zpath.read_bytes()).hexdigest()
        if want.get(zpath.name) not in (None, h):
            sys.exit(f"checksum mismatch for {zpath.name}: expected {want[zpath.name]}, got {h}")
    with zipfile.ZipFile(zpath) as z:
        z.extractall(ROOT)
        bad = [n for n in z.namelist() if n in want and hashlib.sha256((ROOT / n).read_bytes()).hexdigest() != want[n]]
    if bad:
        sys.exit(f"{len(bad)} files failed checksum, e.g. {bad[:3]}")
    print(f"restored {len(z.namelist())} files into {ROOT/'work'} (checksums OK)")


if __name__ == "__main__":
    main()
