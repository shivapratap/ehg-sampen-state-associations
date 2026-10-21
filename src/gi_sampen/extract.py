#!/usr/bin/env python3
"""Stage A2 -- compute the Paper 2 feature block over the frozen analysis set.

Runs over the 4,083 EHG9/FWH analysis units. Two modes matter:

    --parity N   Step 0. Recompute the FROZEN summaries on N windows and check
                 them against the cached Stage 04 columns. Nothing downstream
                 is trustworthy until this passes. Run it first.

    --probe N    Time N windows and extrapolate, so you know what the full run
                 costs before committing to it.

Then the full run with no flags.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
import pandas as pd

from . import paths as _common
from . import gi_profile

PARITY_ATOL = 1e-10
PARITY_RTOL = 1e-8


def load_grid() -> np.ndarray:
    path = _common.RESULTS_DIR / "grid_definition.json"
    if not path.exists():
        raise RuntimeError(
            f"{path} not found. Run select_grid.py before extract_gi.py -- the "
            "grid must be fixed from baseline windows before UC/FM features "
            "are computed.")
    definition = json.loads(path.read_text())
    print(f"grid: tau in [{definition['tau_min']:.3f}, "
          f"{definition['tau_max']:.3f}] sd, "
          f"{definition['n_grid_points']} points "
          f"(fixed {definition['created_utc']})", flush=True)
    return np.asarray(definition["grid"], dtype=np.float64)


# ===========================================================================
# Step 0 -- parity against the frozen cached values
# ===========================================================================

def run_parity(n_windows: int) -> int:
    frame = _common.load_analysis_windows()
    rng = np.random.default_rng(1)
    picked = np.sort(rng.choice(len(frame), size=min(n_windows, len(frame)),
                                replace=False))
    sample = frame.iloc[picked].copy()

    cached = pd.read_csv(
        _common.FEATURES_PATH,
        usecols=["window_id"] + list(_common.FROZEN_COLUMNS.values()),
    ).set_index("window_id")

    rows = []
    for row, signal in _common.iter_windows(sample):
        computed = gi_profile.frozen_summaries(signal)
        reference = cached.loc[row.window_id]
        record = {"window_id": row.window_id}
        for name, column in _common.FROZEN_COLUMNS.items():
            got = computed[f"frozen_{name}"]
            want = float(reference[column])
            record[f"{name}_computed"] = got
            record[f"{name}_frozen"] = want
            record[f"{name}_abs_diff"] = abs(got - want)
            record[f"{name}_match"] = bool(
                np.isclose(got, want, atol=PARITY_ATOL, rtol=PARITY_RTOL)
                or (np.isnan(got) and np.isnan(want)))
        rows.append(record)

    report = pd.DataFrame(rows)
    _common.write_csv(report, _common.RESULTS_DIR / "parity_check.csv",
                      label="parity check")

    print()
    print("PARITY AGAINST FROZEN STAGE 04")
    failed = 0
    for name in _common.FROZEN_COLUMNS:
        matches = int(report[f"{name}_match"].sum())
        worst = float(np.nanmax(report[f"{name}_abs_diff"]))
        status = "OK  " if matches == len(report) else "FAIL"
        if matches != len(report):
            failed += 1
        print(f"  {status} {name:<10} {matches}/{len(report)} match, "
              f"max |diff| = {worst:.3e}")
    print()
    if failed:
        print("PARITY FAILED. Stop here -- window reconstruction or the "
              "reimplementation differs from the frozen pipeline. Do not "
              "interpret any downstream number.")
        return 1
    print("PARITY PASSED. Window reconstruction and the reimplementation "
          "match the frozen pipeline. Safe to proceed.")
    return 0


# ===========================================================================
# Stage A2 -- the feature block
# ===========================================================================

def run_extract(limit: int | None, probe: int | None, quiet: bool) -> int:
    taus = load_grid()
    n_requested = probe if probe is not None else limit
    frame = _common.load_analysis_windows(limit=n_requested)

    started = time.perf_counter()
    rows = []
    for row, signal in _common.iter_windows(frame, progress=not quiet):
        features = gi_profile.all_features(signal, taus)
        features.update({
            "window_id": row.window_id,
            "recording_id": row.recording_id,
            "patient_id": row.patient_id,
            "window_class": row.window_class,
            "window_start_s": row.window_start_s,
        })
        rows.append(features)
    elapsed = time.perf_counter() - started

    result = pd.DataFrame(rows)
    lead = ["window_id", "recording_id", "patient_id", "window_class",
            "window_start_s"]
    result = result[lead + [c for c in result.columns if c not in lead]]

    if probe is not None:
        total = len(_common.load_analysis_windows())
        per_window = elapsed / max(len(result), 1)
        print()
        print(f"PROBE: {len(result)} windows in {elapsed:.1f}s "
              f"({per_window*1000:.0f} ms/window)")
        print(f"  projected full run over {total} windows: "
              f"{per_window*total/60:.1f} min")
        print("  (the projection is optimistic -- the probe reuses the first "
              "few recordings, so per-recording read/filter cost is amortised "
              "over more windows than average)")
        return 0

    _common.write_csv(result, _common.WORK_DIR / "gi_features_60s_EHG9.csv",
                      label="GI features")
    print(f"elapsed: {elapsed/60:.1f} min")

    finite = result.filter(like="gi_").notna().mean()
    print()
    print("finite fraction by feature:")
    for name, value in finite.items():
        flag = "" if value > 0.98 else "   <-- check"
        print(f"  {name:<28} {value:.3f}{flag}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parity", type=int, metavar="N",
                        help="Step 0: check N windows against frozen values")
    parser.add_argument("--probe", type=int, metavar="N",
                        help="time N windows and project the full run")
    parser.add_argument("--limit", type=int, metavar="N",
                        help="process only the first N windows")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    if args.parity:
        return run_parity(args.parity)
    return run_extract(args.limit, args.probe, args.quiet)


if __name__ == "__main__":
    raise SystemExit(main())
