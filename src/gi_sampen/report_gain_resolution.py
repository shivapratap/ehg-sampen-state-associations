#!/usr/bin/env python3
"""Section 7 (report stage) -- turn gain_resolution_sensitivity.csv into a
compact aggregate-by-resolution table and a short markdown report.

Run after gain_resolution_sensitivity.py. Pure aggregation; no new SampEn
computation, no new data read.
"""
from __future__ import annotations

import pandas as pd

from . import paths as _common

RESOLUTION_ORDER = ["frozen_d3", "frozen_d4", "frozen_d5", "gi_d1", "gi_d2", "gi_d3"]


def main() -> int:
    path = _common.RESULTS_DIR / "gain_resolution_sensitivity.csv"
    if not path.exists():
        raise RuntimeError(f"{path} not found. Run gain_resolution_sensitivity.py first.")
    frame = pd.read_csv(path)

    # Aggregate by resolution: worst-case (max) relative deviation and overall
    # pass rate, pooled across the six gains and three summaries, per mode.
    agg = frame.groupby(["resolution", "mode"], sort=False).agg(
        approx_bin_count_documented=("approx_bin_count_documented", "first"),
        approx_bin_count_measured_median=("approx_bin_count_measured_median", "first"),
        worst_median_relative_difference_percent=("median_relative_difference_percent", lambda s: s.abs().max()),
        worst_max_relative_difference_percent=("max_abs_relative_difference_percent", "max"),
        pass_count=("pass_count", "sum"),
        total_count=("total_count", "sum"),
    ).reset_index()
    agg["proportion_pass"] = agg["pass_count"] / agg["total_count"]
    agg["order"] = agg["resolution"].map({r: i for i, r in enumerate(RESOLUTION_ORDER)})
    agg = agg.sort_values("order").drop(columns="order")

    lines = []
    lines.append("# Gain sensitivity versus numerical resolution\n")
    lines.append(
        "Same 60-window gain-stress subset as Stage B and the AvgSampEn gain "
        "audit; same six non-unit gains (0.25, 0.5, 1.3, 2, 4, 7.1); Total, "
        "Average and SD SampEn only. Frozen and GI profiling evaluated at the "
        "same decimals settings used in Stage E's resolution ablation.\n")
    lines.append("## Aggregate by resolution (worst case across gain x summary)\n")
    lines.append("| resolution | mode | approx bins (doc) | approx bins (measured) | "
                 "worst median rel.diff % | worst max rel.diff % | pass rate |")
    lines.append("|---|---|---:|---:|---:|---:|---:|")
    for _, r in agg.iterrows():
        lines.append(
            f"| {r['resolution']} | {r['mode']} | {r['approx_bin_count_documented']:.0f} | "
            f"{r['approx_bin_count_measured_median']:.0f} | "
            f"{r['worst_median_relative_difference_percent']:.4f} | "
            f"{r['worst_max_relative_difference_percent']:.4f} | "
            f"{r['proportion_pass']:.3f} |")

    lines.append("\n## Per feature x resolution x gain\n")
    lines.append("| resolution | mode | feature | gain | median rel.diff % | "
                 "max rel.diff % | pass/total |")
    lines.append("|---|---|---|---:|---:|---:|---:|")
    detail = frame.copy()
    detail["order"] = detail["resolution"].map({r: i for i, r in enumerate(RESOLUTION_ORDER)})
    detail = detail.sort_values(["order", "feature", "gain"])
    for _, r in detail.iterrows():
        lines.append(
            f"| {r['resolution']} | {r['mode']} | {r['feature']} | {r['gain']:.2f} | "
            f"{r['median_relative_difference_percent']:.4f} | "
            f"{r['max_abs_relative_difference_percent']:.4f} | "
            f"{int(r['pass_count'])}/{int(r['total_count'])} |")

    lines.append("\n## Questions (answer by reading the tables above, not asserted here)\n")
    lines.append("1. Does finer absolute resolution (frozen_d3 -> d4 -> d5) reduce "
                 "frozen gain sensitivity, and does it reach zero at any tested resolution?")
    lines.append("2. Is GI profiling's worst-case relative deviation at or below the "
                 "numerical-invariance threshold (rtol 1e-8) at every tested resolution "
                 "(gi_d1, gi_d2, gi_d3), independent of bin count?")
    lines.append("3. Is the frozen-vs-GI gap consistent with gain dependence arising from "
                 "the interaction of signal scale and FIXED ABSOLUTE quantization, as "
                 "opposed to quantization in general?")
    lines.append("4. Is the frozen trend monotonic in resolution? Report what the numbers "
                 "show; do not force a monotonic reading if they are not monotonic.")

    report_path = _common.AUDIT_DIR / "gain_resolution_sensitivity_report.md"
    report_path.write_text("\n".join(lines) + "\n")
    _common.write_csv(agg.drop(columns=["pass_count", "total_count"]),
                      _common.RESULTS_DIR / "gain_resolution_sensitivity_aggregate.csv",
                      label="gain x resolution aggregate")
    print(f"wrote {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
