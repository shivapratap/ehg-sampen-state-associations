#!/usr/bin/env python3
"""Section 6 -- CONDITIONAL AvgSampEn extension of Stage E2.

Do NOT run this until the AvgSampEn gain-stress result (avg_sampen_gain_*)
has been inspected. Only proceeds if frozen AvgSampEn turns out to be
approximately invariant in practice -- if it is clearly gain-sensitive, this
script is unnecessary (the gain-stress control already settles the point)
and should not be run.

This is a SEPARATE script from resolution_ablation_report.py (Stage E2),
deliberately -- Stage E2's output (results/resolution_ablation_summary.csv)
is already referenced elsewhere (facts_pack.md, README.md) and is left
untouched. This script reuses the exact same model specification and reuses
work/resolution_ablation_features.csv, which ALREADY CONTAINS the average
columns at every resolution (gi_d1_average, gi_d2_average, frozen_d4_average,
frozen_d5_average were computed by the existing resolution_ablation.py --
frozen_summaries/gi_summaries_rounded return all six summaries, average
included, so no new feature extraction is needed here, only new model
fitting). gi_d3 average and frozen_d3 average are already fit and cached
(results/gi_model_results.csv's gi_average rows; Paper1's frozen
mixed_model_results_60s_fwh_EHG9.csv's fwh_entropy_averagesampen rows) and
are pulled in the same way Stage E2 pulls its d3 baselines.

Writes results/resolution_ablation_summary_average.csv -- a NEW file, not a
modification of Stage E2's resolution_ablation_summary.csv.
"""
from __future__ import annotations

import pandas as pd

from . import paths as _common
from .statistics import fit_one, apply_fdr

CONDITION_ORDER = {
    "frozen_d3": 0, "frozen_d4": 1, "frozen_d5": 2,
    "gi_d1": 3, "gi_d2": 4, "gi_d3": 5,
}
NEW_CONDITIONS = ("gi_d1", "gi_d2", "frozen_d4", "frozen_d5")


def fit_new_conditions(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for condition in NEW_CONDITIONS:
        column = f"{condition}_average"
        if column not in frame.columns:
            print(f"  skipping {column}: not in feature table")
            continue
        for record in fit_one(frame, column):
            record["condition"] = condition
            record["summary"] = "average"
            rows.append(record)
    fitted = apply_fdr(pd.DataFrame(rows))
    return fitted.loc[fitted["contrast"] == "fetal_movement_vs_baseline"].copy()


def load_gi_d3_average() -> pd.DataFrame:
    path = _common.RESULTS_DIR / "gi_model_results.csv"
    gi = pd.read_csv(path)
    gi = gi.loc[(gi["contrast"] == "fetal_movement_vs_baseline")
               & (gi["feature"] == "gi_average")].copy()
    gi["condition"] = "gi_d3"
    gi["summary"] = "average"
    return gi


def load_frozen_d3_average() -> pd.DataFrame:
    path = (_common.PAPER1_ROOT / "results" / "primary"
           / "mixed_model_results_60s_fwh_EHG9.csv")
    frozen = pd.read_csv(path)
    frozen = frozen.loc[frozen["channel"] == _common.CHANNEL]
    match = frozen.loc[
        (frozen["feature"] == "fwh_entropy_averagesampen")
        & (frozen["contrast"] == "fetal_movement_vs_baseline")]
    if match.empty:
        return pd.DataFrame()
    record = match.iloc[0]
    return pd.DataFrame([{
        "feature": "fwh_entropy_averagesampen", "contrast": "fetal_movement_vs_baseline",
        "estimate": float(record["estimate"]), "q_value": float(record["q_value"]),
        "n": int(record["n"]), "n_women": int(record["n_women"]),
        "condition": "frozen_d3", "summary": "average",
    }])


def main() -> int:
    features_path = _common.WORK_DIR / "resolution_ablation_features.csv"
    if not features_path.exists():
        raise RuntimeError(
            f"{features_path} not found. Run resolution_ablation.py first "
            "(it already exists from Stage E; no re-run should be needed).")
    frame = pd.read_csv(features_path)

    print("fitting AvgSampEn at the 4 new resolutions...", flush=True)
    new_results = fit_new_conditions(frame)
    gi_d3 = load_gi_d3_average()
    frozen_d3 = load_frozen_d3_average()

    keep = ["feature", "contrast", "estimate", "q_value", "condition", "summary", "n", "n_women"]
    parts = [p[keep] for p in (new_results, gi_d3, frozen_d3) if not p.empty]
    combined = pd.concat(parts, ignore_index=True)
    combined["order"] = combined["condition"].map(CONDITION_ORDER)
    combined = combined.sort_values("order").reset_index(drop=True)

    print()
    print("=" * 78)
    print("AvgSampEn ACROSS RESOLUTION -- FM versus IN, EHG9/FWH")
    print("=" * 78)
    print(f"{'condition':<12}{'invariant?':<12}{'estimate':>10}{'q_value':>12}{'sig (q<.05)':>13}")
    gi_sig, frozen_sig = [], []
    for _, r in combined.iterrows():
        invariant = r["condition"].startswith("gi")
        sig = bool(pd.notna(r["q_value"]) and r["q_value"] < 0.05)
        (gi_sig if invariant else frozen_sig).append(sig)
        print(f"{r['condition']:<12}{'yes' if invariant else 'no':<12}"
              f"{r['estimate']:>10.3f}{r['q_value']:>12.4f}{'yes' if sig else 'no':>13}")

    if all(gi_sig) and not any(frozen_sig):
        verdict = "INVARIANCE-DRIVEN -- like GI-Total, treat as a potentially story-changing result"
    elif not any(gi_sig):
        verdict = "NOT SUPPORTED at any GI resolution"
    elif gi_sig and gi_sig[-1] and not gi_sig[0]:
        verdict = "RESOLUTION-SENSITIVE -- only significant at the finest GI resolution"
    elif any(frozen_sig):
        verdict = "AMBIGUOUS -- resolution alone also reaches significance somewhere"
    else:
        verdict = "MIXED -- inspect the table directly"
    print()
    print(f"verdict: {verdict}")

    combined.drop(columns=["order"]).to_csv(
        _common.RESULTS_DIR / "resolution_ablation_summary_average.csv", index=False)
    print()
    print("wrote results/resolution_ablation_summary_average.csv "
          f"({len(combined)} rows) -- send this one back")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
