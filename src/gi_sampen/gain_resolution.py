#!/usr/bin/env python3
"""Section 7 -- gain sensitivity versus numerical (quantization) resolution.

Answers: does finer absolute rounding reduce the frozen construction's gain
sensitivity, and does GI profiling stay invariant regardless of resolution?

Reuses, with no new SampEn implementation:
  * the exact 60-window gain-stress subset from invariance_check.select_windows
    (imported, not re-implemented) -- same windows as Stage B and Section 5;
  * gi_profile.frozen_summaries / gi_profile.gi_summaries_rounded, called at
    the SAME decimals settings already defined and named by
    resolution_ablation.py's CONDITIONS (gi_d1=1, gi_d2=2, gi_d3=3 [the
    default GI_DISTANCE_DECIMALS], frozen_d3=3 [the default
    PROFILE_DISTANCE_DECIMALS], frozen_d4=4, frozen_d5=5) -- these are the
    exact settings that produced the ~60/~600/~5,000 and ~80/~800/~7,000 bin
    counts already reported in Paper2/README.md's Stage E table and in
    resolution_ablation.py's own docstring. No new resolution settings are
    invented here.

What is new here (not present anywhere else in Paper2): looping the SAME
six resolution conditions over the SAME six non-unit gain factors used by
Stage B (0.25, 0.5, 1.3, 2, 4, 7.1), for Total/Average/SD SampEn only
(Kurtosis and Skewness already failed the resolution ablation in Stage E and
are out of scope per the task prompt). resolution_ablation.py itself only
ever runs at gain = 1x (its job is the FM-vs-IN mixed model across the full
4,083-window analysis set, not a gain-stress test), so this loop does not
exist anywhere else and has to be added.

Relative-difference and invariance-pass definitions are copied verbatim from
invariance_check.py / avg_sampen_gain_report.py (percent relative change,
1e-10 floor on the denominator; np.isclose atol=1e-10, rtol=1e-8).
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from . import paths as _common
from . import gi_profile
from .gain_stress import select_windows

# (label, function, kwargs, output-key prefix used by that function's return dict)
RESOLUTIONS = (
    ("frozen_d3", gi_profile.frozen_summaries, dict(decimals=3), "frozen"),
    ("frozen_d4", gi_profile.frozen_summaries, dict(decimals=4), "frozen"),
    ("frozen_d5", gi_profile.frozen_summaries, dict(decimals=5), "frozen"),
    ("gi_d1", gi_profile.gi_summaries_rounded, dict(decimals=1), "gi"),
    ("gi_d2", gi_profile.gi_summaries_rounded, dict(decimals=2), "gi"),
    ("gi_d3", gi_profile.gi_summaries_rounded, dict(decimals=3), "gi"),
)
# Approximate median occupied-bin counts already documented in README.md's
# Stage E table -- reported alongside the measured median here for
# cross-checking, not assumed.
APPROX_BIN_COUNT = {
    "frozen_d3": 80, "frozen_d4": 800, "frozen_d5": 7000,
    "gi_d1": 60, "gi_d2": 600, "gi_d3": 5000,
}
SUMMARIES = ("total", "average", "sd")
GAINS = (0.25, 0.5, 1.3, 2.0, 4.0, 7.1)
INVARIANCE_ATOL = 1e-10
INVARIANCE_RTOL = 1e-8
NEAR_ZERO_DENOMINATOR_FLOOR = 1e-10


def relabel(raw: dict, new_prefix: str, old_prefix: str) -> dict:
    cut = len(old_prefix) + 1  # strip "prefix_"
    return {f"{new_prefix}_{key[cut:]}": value for key, value in raw.items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probe", type=int, metavar="N",
                        help="time N windows (all gains x resolutions) and stop")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    frame = _common.load_analysis_windows()
    selected = select_windows(frame)
    if args.probe is not None:
        selected = selected.head(args.probe)

    print(f"{len(selected)} windows x {len(GAINS) + 1} gains "
          f"x {len(RESOLUTIONS)} resolutions", flush=True)

    started = time.perf_counter()
    rows = []
    for row, signal in _common.iter_windows(selected, progress=not args.quiet):
        for gain in (1.0,) + GAINS:
            scaled = gain * signal
            record = {
                "window_id": row.window_id, "state": row.window_class,
                "patient_id": row.patient_id, "gain": gain,
            }
            for label, fn, kwargs, old_prefix in RESOLUTIONS:
                raw = fn(scaled, **kwargs)
                record.update(relabel(raw, label, old_prefix))
            rows.append(record)
    elapsed = time.perf_counter() - started

    long = pd.DataFrame(rows)

    if args.probe is not None:
        per_window = elapsed / max(len(selected), 1)
        print(f"PROBE: {len(selected)} windows in {elapsed:.1f}s "
              f"({per_window*1000:.0f} ms/window, all gains+resolutions)")
        print(f"  projected full 60-window run: {per_window*60:.1f}s")
        return 0

    _common.write_csv(long, _common.WORK_DIR / "gain_resolution_window_results.csv",
                      label="gain x resolution per-window results")

    base = long.loc[long["gain"] == 1.0].set_index("window_id")
    summary_rows = []
    for label, _fn, _kwargs, old_prefix in RESOLUTIONS:
        mode = "GI" if label.startswith("gi") else "frozen"
        measured_bins = float(long.loc[long["gain"] == 1.0, f"{label}_n_retained"].median())
        for summary in SUMMARIES:
            column = f"{label}_{summary}"
            baseline = base[column]
            for gain in GAINS:
                subset = long.loc[long["gain"] == gain].set_index("window_id")
                want = baseline.reindex(subset.index)
                got = subset[column]
                with np.errstate(divide="ignore", invalid="ignore"):
                    relative = 100.0 * (got - want) / want.abs().where(
                        want.abs() > NEAR_ZERO_DENOMINATOR_FLOOR)
                passed = np.isclose(got, want, atol=INVARIANCE_ATOL, rtol=INVARIANCE_RTOL)
                rel_valid = relative.dropna()
                summary_rows.append({
                    "resolution": label,
                    "approx_bin_count_documented": APPROX_BIN_COUNT[label],
                    "approx_bin_count_measured_median": measured_bins,
                    "mode": mode,
                    "feature": summary,
                    "gain": gain,
                    "median_relative_difference_percent": float(rel_valid.median()) if len(rel_valid) else float("nan"),
                    "max_abs_relative_difference_percent": float(rel_valid.abs().max()) if len(rel_valid) else float("nan"),
                    "pass_count": int(passed.sum()),
                    "total_count": int(len(passed)),
                })

    summary = pd.DataFrame(summary_rows)
    _common.write_csv(summary, _common.RESULTS_DIR / "gain_resolution_sensitivity.csv",
                      label="gain x resolution sensitivity summary")

    print(f"elapsed: {elapsed/60:.1f} min")
    print("Next: python gain_resolution_sensitivity_report.py")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
