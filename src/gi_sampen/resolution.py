#!/usr/bin/env python3
"""Stage E1 -- resolution ablation: separate invariance from bin count.

The problem
-----------
Normalising by sd does two things, not one: it removes gain-dependence AND it
changes the effective resolution of the 3-decimal rounding grid (EHG sd is of
order 1e-2, so a 0.001-sd grid is roughly 60x finer than the frozen
0.001-signal-unit grid). The gate result -- GI-SD survives, GI-Total flips
sign -- could be a genuine invariance effect, or it could be a resolution
artefact riding along with it. This stage tells them apart.

Method
------
Hold each property fixed in turn and vary the other:

    frozen_d3   (Paper1 baseline, already cached)   ~80 bins     NOT invariant
    frozen_d4   (new, this stage)                    ~800 bins    NOT invariant
    frozen_d5   (new, this stage)                    ~8,000 bins  NOT invariant
    gi_d1       (new, this stage)                     ~60 bins     invariant
    gi_d2       (new, this stage)                     ~600 bins    invariant
    gi_d3       (already computed, Stage A2)          ~5,000 bins  invariant

Read together with Stage E2's fitted models:

  * If a feature's significance holds across frozen_d3/d4/d5 despite gain-
    sensitivity being unchanged, resolution alone can rescue a signal --
    concerning, but distinguishable from what the GI features show.
  * If a feature's significance holds across gi_d1/d2/d3 regardless of bin
    count, the effect is resolution-independent and the invariance framing
    is supported.
  * If a feature only reaches significance at the finest resolution in BOTH
    families, that specific feature's "signal" is resolution-driven, not an
    invariance effect, whichever family it appears in.

Run after Stage A2 (extract_gi.py) and before deciding the paper's feature
set. Single pass over the 4,083 analysis windows; no new WFDB reads beyond
this stage's own pass.
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

from . import paths as _common
from . import gi_profile

CONDITIONS = (
    ("gi_d1", gi_profile.gi_summaries_rounded, dict(decimals=1), "gi"),
    ("gi_d2", gi_profile.gi_summaries_rounded, dict(decimals=2), "gi"),
    ("frozen_d4", gi_profile.frozen_summaries, dict(decimals=4), "frozen"),
    ("frozen_d5", gi_profile.frozen_summaries, dict(decimals=5), "frozen"),
)


def relabel(raw: dict, new_prefix: str, old_prefix: str) -> dict:
    cut = len(old_prefix) + 1  # strip "prefix_"
    return {f"{new_prefix}_{key[cut:]}": value for key, value in raw.items()}


def reference_bin_counts() -> None:
    """Print median n_retained for the two already-computed baselines, for
    context alongside the new conditions. Best-effort -- missing files are
    skipped rather than treated as an error."""
    gi_d3_path = _common.WORK_DIR / "gi_features_60s_EHG9.csv"
    if gi_d3_path.exists():
        gi_d3 = pd.read_csv(gi_d3_path, usecols=["gi_n_retained"])
        print(f"  gi_d3        median n_retained = "
              f"{gi_d3['gi_n_retained'].median():.0f}  (all 4,083 windows, Stage A2)")

    scale_path = (_common.PAPER1_ROOT / "work" / "entropy"
                 / "scale_sensitivity_window_results.csv")
    if scale_path.exists():
        scale = pd.read_csv(scale_path, usecols=["scale_factor", "n_retained_profile_values"])
        at_1x = scale.loc[scale["scale_factor"] == 1.0, "n_retained_profile_values"]
        if len(at_1x):
            print(f"  frozen_d3    median n_retained = {at_1x.median():.0f}  "
                  f"(60-window scale-sensitivity subset, Paper1)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, metavar="N")
    parser.add_argument("--probe", type=int, metavar="N")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    n_requested = args.probe if args.probe is not None else args.limit
    frame = _common.load_analysis_windows(limit=n_requested)

    started = time.perf_counter()
    rows = []
    for row, signal in _common.iter_windows(frame, progress=not args.quiet):
        record = {
            "window_id": row.window_id, "recording_id": row.recording_id,
            "patient_id": row.patient_id, "window_class": row.window_class,
        }
        for label, fn, kwargs, old_prefix in CONDITIONS:
            raw = fn(signal, **kwargs)
            record.update(relabel(raw, label, old_prefix))
        rows.append(record)
    elapsed = time.perf_counter() - started

    result = pd.DataFrame(rows)

    if args.probe is not None:
        total = len(_common.load_analysis_windows())
        per_window = elapsed / max(len(result), 1)
        print()
        print(f"PROBE: {len(result)} windows in {elapsed:.1f}s "
              f"({per_window*1000:.0f} ms/window)")
        print(f"  projected full run: {per_window*total/60:.1f} min")
        return 0

    _common.write_csv(result,
                      _common.WORK_DIR / "resolution_ablation_features.csv",
                      label="resolution ablation features")

    print()
    print("median n_retained by condition (resolution check):")
    for label, _fn, _kwargs, _old in CONDITIONS:
        column = f"{label}_n_retained"
        print(f"  {label:<12} median n_retained = {result[column].median():.0f}")
    reference_bin_counts()
    print(f"elapsed: {elapsed/60:.1f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
