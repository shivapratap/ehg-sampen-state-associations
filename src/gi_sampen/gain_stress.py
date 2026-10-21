#!/usr/bin/env python3
"""Stage B -- gain-invariance correctness test.

This is a TEST ON THE IMPLEMENTATION, not a scientific result. The GI features
are invariant by construction; if they are not invariant here, there is a bug
(most likely sd computed at the wrong point) and nothing downstream should be
read.

It reuses Paper 1's exact 60-window selection (20 UC / 20 FM / 20 IN, same
deterministic rule as _sampen_scale_impl.select_windows) and the same five
scale factors, so the output table is directly comparable to manuscript
Table S16 and can be dropped into the ICASSP paper as the side-by-side.

Note on the scale factors: 0.25, 0.5, 2 and 4 are exact in binary floating
point, so the GI features should agree to the last bit. --extra-scales adds
1.3 and 7.1, which are not, as a harder test.
"""
from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from . import paths as _common
from . import gi_profile

SCALES = (0.25, 0.5, 1.0, 2.0, 4.0)
EXTRA_SCALES = (1.3, 7.1)
N_PER_CLASS = 20
CLASS_MAP = {"contraction": "UC", "baseline": "IN", "fetal_movement": "FM"}
CLASS_ORDER = ("UC", "IN", "FM")
INVARIANCE_ATOL = 1e-10
INVARIANCE_RTOL = 1e-8


def select_windows(frame: pd.DataFrame) -> pd.DataFrame:
    """Reproduces Paper 1 `_sampen_scale_impl.select_windows` exactly."""
    eligible = frame.loc[frame["window_class"].isin(CLASS_MAP)].copy()
    eligible["class_label"] = eligible["window_class"].map(CLASS_MAP)
    selected = []
    for label in CLASS_ORDER:
        group = eligible.loc[eligible["class_label"] == label].sort_values(
            ["patient_id", "recording_id", "window_start_s", "window_id"],
            kind="mergesort")
        group["within_patient_rank"] = group.groupby(
            "patient_id", sort=False).cumcount()
        group = group.sort_values(
            ["within_patient_rank", "patient_id", "recording_id",
             "window_start_s", "window_id"], kind="mergesort").head(N_PER_CLASS)
        if len(group) != N_PER_CLASS:
            raise RuntimeError(
                f"only {len(group)} eligible {label} windows; "
                f"{N_PER_CLASS} required")
        selected.append(group)
    return pd.concat(selected, ignore_index=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extra-scales", action="store_true",
                        help="also test non-dyadic gains (1.3, 7.1)")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    scales = SCALES + (EXTRA_SCALES if args.extra_scales else ())
    taus = np.asarray(
        __import__("json").loads(
            (_common.RESULTS_DIR / "grid_definition.json").read_text())["grid"],
        dtype=np.float64)

    frame = _common.load_analysis_windows()
    selected = select_windows(frame)
    print(f"{len(selected)} windows x {len(scales)} scales", flush=True)

    rows = []
    for row, signal in _common.iter_windows(selected, progress=not args.quiet):
        label = CLASS_MAP[row.window_class]
        for scale in scales:
            features = gi_profile.all_features(
                scale * signal, taus, include_frozen=True)
            features.update({
                "window_id": row.window_id, "class_label": label,
                "patient_id": row.patient_id, "scale_factor": scale,
            })
            rows.append(features)

    long = pd.DataFrame(rows)
    feature_columns = [c for c in long.columns if c not in (
        "window_id", "class_label", "patient_id", "scale_factor")
        and not c.endswith("_n_retained")
        and not c.endswith("_max_rounded_distance")
        and not c.endswith("_valid_fraction")]

    reference = long.loc[long["scale_factor"] == 1.0].set_index("window_id")
    records = []
    for column in feature_columns:
        base = reference[column]
        for scale in scales:
            subset = long.loc[long["scale_factor"] == scale].set_index("window_id")
            got, want = subset[column], base.reindex(subset.index)
            with np.errstate(divide="ignore", invalid="ignore"):
                relative = 100.0 * (got - want) / want.abs().where(
                    want.abs() > 1e-10)
            invariant = np.isclose(got, want, atol=INVARIANCE_ATOL,
                                   rtol=INVARIANCE_RTOL)
            records.append({
                "feature": column,
                "scale_factor": scale,
                "n_windows": int(got.notna().sum()),
                "median_relative_change_percent": float(np.nanmedian(relative)),
                "q1_relative_change_percent": float(np.nanpercentile(relative.dropna(), 25))
                    if relative.notna().any() else float("nan"),
                "q3_relative_change_percent": float(np.nanpercentile(relative.dropna(), 75))
                    if relative.notna().any() else float("nan"),
                "max_absolute_difference": float(np.nanmax(np.abs(got - want))),
                "n_numerically_invariant": int(np.sum(invariant)),
                "proportion_invariant": float(np.mean(invariant)),
            })

    summary = pd.DataFrame(records)
    _common.write_csv(long, _common.WORK_DIR / "invariance_window_results.csv",
                      label="per-window invariance")
    _common.write_csv(summary, _common.RESULTS_DIR / "invariance_summary.csv",
                      label="invariance summary")

    print()
    print("GAIN-INVARIANCE CORRECTNESS TEST")
    print(f"{'feature':<26}{'worst |diff|':>14}{'invariant':>12}")
    failures = []
    for column in feature_columns:
        block = summary.loc[summary["feature"] == column]
        worst = float(block["max_absolute_difference"].max())
        proportion = float(block["proportion_invariant"].min())
        is_gi = column.startswith("gi")
        status = "PASS" if proportion == 1.0 else "----"
        if is_gi and proportion < 1.0:
            status = "FAIL"
            failures.append(column)
        print(f"{column:<26}{worst:>14.3e}{proportion:>10.2f}  {status}")

    print()
    if failures:
        print("FAILED -- these GI features are not invariant: "
              + ", ".join(failures))
        print("This is an implementation bug. Fix before reading Stage C.")
        return 1
    print("PASSED -- every GI feature is numerically invariant at every scale.")
    print("(Frozen features are expected to move; that is the Table S16 result.)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
