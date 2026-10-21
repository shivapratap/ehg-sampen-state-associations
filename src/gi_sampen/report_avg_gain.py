#!/usr/bin/env python3
"""Section 5 (ICASSP reviewer-objection prompt) -- AvgSampEn gain-sensitivity
audit, derived entirely from Stage B's existing output.

Reuses, with no new SampEn computation whatsoever:
  * the same 60-window gain-stress subset (Paper 1's frozen 20 UC / 20 FM /
    20 IN selection, reproduced by invariance_check.select_windows);
  * the same six frozen/GI SampEn summaries (frozen_total/average/sd,
    gi_total/average/sd), already computed by gi_profile.py and already
    written by invariance_check.py to work/invariance_window_results.csv;
  * the identical relative-difference and invariance-pass definitions used
    in invariance_check.py (percent relative change with a 1e-10 floor on
    the denominator; np.isclose with atol=1e-10, rtol=1e-8). Nothing about
    that definition is changed here.

Precondition: invariance_check.py must have been run WITH --extra-scales,
so that work/invariance_window_results.csv contains all six non-unit gains
required here (0.25, 0.5, 1.3, 2, 4, 7.1), not just the five dyadic ones.

    python invariance_check.py --extra-scales

This script only reshapes/derives; it does not read raw EHG data and does
not touch the frozen Paper 1 pipeline or any cached Paper 1 value.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import paths as _common

# feature name -> (frozen column, GI column) in invariance_window_results.csv
FEATURES = {
    "total": ("frozen_total", "gi_total"),
    "average": ("frozen_average", "gi_average"),
    "sd": ("frozen_sd", "gi_sd"),
}
REQUIRED_GAINS = (0.25, 0.5, 1.3, 2.0, 4.0, 7.1)
# Same tolerances as invariance_check.py -- do not change independently of it.
INVARIANCE_ATOL = 1e-10
INVARIANCE_RTOL = 1e-8
NEAR_ZERO_DENOMINATOR_FLOOR = 1e-10


def main() -> int:
    source = _common.WORK_DIR / "invariance_window_results.csv"
    if not source.exists():
        raise RuntimeError(
            f"{source} not found. Run invariance_check.py --extra-scales first.")
    long = pd.read_csv(source)

    present = set(long["scale_factor"].unique())
    missing = sorted(set(REQUIRED_GAINS) - present)
    if missing:
        raise RuntimeError(
            f"{source} is missing gain factor(s) {missing}. "
            "Re-run: python invariance_check.py --extra-scales")

    meta_cols = ["window_id", "class_label", "patient_id"]
    base = long.loc[long["scale_factor"] == 1.0].set_index("window_id")

    detail_frames = []
    for summary, (frozen_col, gi_col) in FEATURES.items():
        for mode, column in (("frozen", frozen_col), ("GI", gi_col)):
            baseline = base[column]
            for gain in REQUIRED_GAINS:
                subset = long.loc[long["scale_factor"] == gain].set_index("window_id")
                want = baseline.reindex(subset.index)
                got = subset[column]
                abs_diff = (got - want).abs()
                with np.errstate(divide="ignore", invalid="ignore"):
                    relative = 100.0 * (got - want) / want.abs().where(
                        want.abs() > NEAR_ZERO_DENOMINATOR_FLOOR)
                passed = np.isclose(got, want, atol=INVARIANCE_ATOL,
                                    rtol=INVARIANCE_RTOL)
                block = pd.DataFrame({
                    "window_id": subset.index,
                    "state": subset["class_label"].values,
                    "gain": gain,
                    "feature": summary,
                    "mode": mode,
                    "baseline_value": want.values,
                    "scaled_value": got.values,
                    "abs_difference": abs_diff.values,
                    "relative_difference": relative.values,
                    "invariance_pass": passed,
                })
                detail_frames.append(block)

    detail = pd.concat(detail_frames, ignore_index=True)
    _common.write_csv(detail, _common.RESULTS_DIR / "avg_sampen_gain_test.csv",
                      label="AvgSampEn gain test (per window)")

    summary_rows = []
    for (feature, mode, gain), block in detail.groupby(["feature", "mode", "gain"]):
        rel = block["relative_difference"].dropna()
        summary_rows.append({
            "feature": feature,
            "mode": mode,
            "gain": gain,
            "n_windows": int(len(block)),
            "median_relative_difference_percent": float(rel.median()) if len(rel) else float("nan"),
            "max_abs_relative_difference_percent": float(rel.abs().max()) if len(rel) else float("nan"),
            "n_invariance_pass": int(block["invariance_pass"].sum()),
            "proportion_invariance_pass": float(block["invariance_pass"].mean()),
        })
    summary = pd.DataFrame(summary_rows).sort_values(["feature", "mode", "gain"])
    _common.write_csv(summary, _common.RESULTS_DIR / "avg_sampen_gain_summary.csv",
                      label="AvgSampEn gain summary")

    print()
    print("=" * 78)
    print("AvgSampEn GAIN-SENSITIVITY AUDIT (frozen vs GI, Total/Average/SD)")
    print("=" * 78)
    for feature in FEATURES:
        print()
        print(f"-- {feature} --")
        print(f"{'mode':<8}{'gain':>7}{'median rel.%':>14}{'max |rel.%|':>13}{'invariant':>12}")
        for mode in ("frozen", "GI"):
            block = summary.loc[(summary["feature"] == feature) & (summary["mode"] == mode)]
            for _, r in block.iterrows():
                print(f"{mode:<8}{r['gain']:>7.2f}{r['median_relative_difference_percent']:>14.4f}"
                      f"{r['max_abs_relative_difference_percent']:>13.4f}"
                      f"{r['proportion_invariance_pass']:>12.2f}")

    # Q1/Q2: is frozen AvgSampEn invariant, and how does it compare to frozen Total?
    def worst_case(feature, mode):
        block = summary.loc[(summary["feature"] == feature) & (summary["mode"] == mode)]
        return float(block["max_abs_relative_difference_percent"].max())

    print()
    print("-" * 78)
    print("Worst-case |relative difference| across all tested gains (percent):")
    for feature in FEATURES:
        print(f"  frozen_{feature:<9} {worst_case(feature, 'frozen'):>10.4f}%   "
              f"GI_{feature:<9} {worst_case(feature, 'GI'):>10.4f}%")
    print()
    print("(Interpretation against the five reviewer questions in Section 5 of")
    print(" the task prompt is left to a human/Claude reading pass over these")
    print(" numbers -- this script deliberately does not pre-judge the result.)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
