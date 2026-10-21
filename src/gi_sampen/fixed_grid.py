#!/usr/bin/env python3
"""Stage A1 -- choose the fixed dimensionless tolerance grid. BASELINE ONLY.

Why this is a separate stage with its own output file
-----------------------------------------------------
The obvious reviewer objection to a fixed-grid representation is that the grid
was tuned to produce the desired result. This stage removes that objection by
construction:

  * the grid is chosen from IN (basal inactivity) windows ONLY -- never from
    the UC or FM windows that the gate is evaluated on;
  * the rule is fixed in advance and stated here: keep the longest contiguous
    run of tolerances at which at least MIN_COVERAGE of sampled baseline
    windows yield a finite, retained profile value;
  * the result is written to results/grid_definition.json with a UTC timestamp
    before any UC or FM number is computed, so the ordering is auditable.

Run this before extract_gi.py. It is cheap (a few minutes).
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from . import paths as _common
from . import gi_profile

# --- Prespecified selection parameters. Do not tune after seeing results. ---
CANDIDATE_MIN = 0.01          # sd units
CANDIDATE_MAX = 4.00          # sd units
CANDIDATE_STEP = 0.01
MIN_COVERAGE = 0.95           # fraction of baseline windows that must be valid
N_BASELINE_SAMPLE = 300
RANDOM_SEED = 20260920
FINAL_GRID_POINTS = 200       # equally spaced points on the selected range

RULE_TEXT = (
    "Candidate tolerances tau in [0.01, 4.00] sd, step 0.01. A tolerance is "
    "covered if the sd-normalised entropy profile yields a finite retained "
    "value there. The selected range is the longest contiguous run of "
    "candidate tolerances covered in at least 95% of a fixed random sample of "
    "300 baseline (IN) windows (seed 20260920). The analysis grid is 200 "
    "equally spaced points on that range. Chosen from baseline windows only, "
    "before any UC or FM feature was computed."
)


def longest_contiguous_run(mask: np.ndarray) -> tuple[int, int]:
    """Return (start, stop_exclusive) of the longest run of True."""
    best_start = best_len = current_start = current_len = 0
    for index, flag in enumerate(mask):
        if flag:
            if current_len == 0:
                current_start = index
            current_len += 1
            if current_len > best_len:
                best_len, best_start = current_len, current_start
        else:
            current_len = 0
    if best_len == 0:
        raise RuntimeError(
            "No tolerance reached the coverage threshold. Widen the candidate "
            "range or lower MIN_COVERAGE -- but record that you did so.")
    return best_start, best_start + best_len


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-baseline", type=int, default=N_BASELINE_SAMPLE)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    frame = _common.load_analysis_windows()
    baseline = frame.loc[frame["window_class"] == "baseline"].copy()
    if len(baseline) < args.n_baseline:
        raise RuntimeError(
            f"only {len(baseline)} baseline windows available; "
            f"{args.n_baseline} requested")

    rng = np.random.default_rng(RANDOM_SEED)
    picked = rng.choice(len(baseline), size=args.n_baseline, replace=False)
    sample = baseline.iloc[np.sort(picked)].copy()

    candidates = np.arange(
        CANDIDATE_MIN, CANDIDATE_MAX + CANDIDATE_STEP / 2, CANDIDATE_STEP)
    coverage_counts = np.zeros(candidates.size, dtype=np.int64)
    n_windows = 0

    print(f"Scanning {len(sample)} baseline windows over {candidates.size} "
          f"candidate tolerances...", flush=True)
    for _row, signal in _common.iter_windows(sample, progress=not args.quiet):
        _profile, valid = gi_profile.gi_profile_on_grid(signal, candidates)
        coverage_counts += valid.astype(np.int64)
        n_windows += 1

    coverage = coverage_counts / float(n_windows)
    start, stop = longest_contiguous_run(coverage >= MIN_COVERAGE)
    tau_min = float(candidates[start])
    tau_max = float(candidates[stop - 1])
    grid = np.linspace(tau_min, tau_max, FINAL_GRID_POINTS)

    window_ids = ",".join(sorted(sample["window_id"].astype(str)))
    definition = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "rule": RULE_TEXT,
        "channel": _common.CHANNEL,
        "band": _common.BAND,
        "min_coverage": MIN_COVERAGE,
        "candidate_min": CANDIDATE_MIN,
        "candidate_max": CANDIDATE_MAX,
        "candidate_step": CANDIDATE_STEP,
        "random_seed": RANDOM_SEED,
        "n_baseline_windows_used": n_windows,
        "baseline_window_id_sha256": hashlib.sha256(
            window_ids.encode("utf-8")).hexdigest(),
        "tau_min": tau_min,
        "tau_max": tau_max,
        "n_grid_points": FINAL_GRID_POINTS,
        "grid": grid.tolist(),
    }
    out_path = _common.RESULTS_DIR / "grid_definition.json"
    out_path.write_text(json.dumps(definition, indent=2))

    coverage_frame = pd.DataFrame({
        "tau": candidates,
        "coverage_fraction": coverage,
        "selected": (np.arange(candidates.size) >= start)
                    & (np.arange(candidates.size) < stop),
    })
    _common.write_csv(coverage_frame,
                      _common.RESULTS_DIR / "grid_coverage_baseline.csv",
                      label="grid coverage")

    print()
    print(f"SELECTED GRID: tau in [{tau_min:.3f}, {tau_max:.3f}] sd, "
          f"{FINAL_GRID_POINTS} points")
    print(f"  from {n_windows} baseline windows, coverage >= {MIN_COVERAGE:.0%}")
    print(f"  written to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
