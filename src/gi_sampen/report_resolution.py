#!/usr/bin/env python3
"""Stage E2 -- fit the resolution ladder and report whether significance
tracks invariance or bin count.

Run after resolution_ablation.py. Reuses the exact model specification from
refit_models.py (fit_one, apply_fdr) so every condition is fit identically to
the gate. Only the FM-versus-IN contrast is reported; that is the contrast
the gate and the whole paper turn on.

This ablation's BH-FDR family is its own family: 4 conditions x 4 summaries
= 16 fits, q-values computed within that set. The frozen_d3 and gi_d3 rows
shown alongside it are pulled in from their own already-computed families
(Paper1's frozen model and Stage C's GI model) for the side-by-side and are
NOT re-included in this family's FDR correction.
"""
from __future__ import annotations

import pandas as pd

from . import paths as _common
from .statistics import fit_one, apply_fdr

SUMMARIES = ("sd", "kurtosis", "skewness", "total")
NEW_CONDITIONS = ("gi_d1", "gi_d2", "frozen_d4", "frozen_d5")
CONDITION_ORDER = {
    "frozen_d3": 0, "frozen_d4": 1, "frozen_d5": 2,
    "gi_d1": 3, "gi_d2": 4, "gi_d3": 5,
}
FROZEN_COLUMN = {
    "sd": "fwh_entropy_sdsampen", "kurtosis": "fwh_entropy_kurtosissampen",
    "skewness": "fwh_entropy_skewnesssampen", "total": "fwh_entropy_totalsampen",
}


def fit_new_conditions(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for condition in NEW_CONDITIONS:
        for summary in SUMMARIES:
            column = f"{condition}_{summary}"
            if column not in frame.columns:
                print(f"  skipping {column}: not in feature table")
                continue
            for record in fit_one(frame, column):
                record["condition"] = condition
                record["summary"] = summary
                rows.append(record)
    fitted = apply_fdr(pd.DataFrame(rows))
    return fitted.loc[fitted["contrast"] == "fetal_movement_vs_baseline"].copy()


def load_gi_d3() -> pd.DataFrame:
    path = _common.RESULTS_DIR / "gi_model_results.csv"
    gi_d3 = pd.read_csv(path)
    gi_d3 = gi_d3.loc[
        (gi_d3["contrast"] == "fetal_movement_vs_baseline")
        & gi_d3["feature"].isin([f"gi_{s}" for s in SUMMARIES])
    ].copy()
    gi_d3["condition"] = "gi_d3"
    gi_d3["summary"] = gi_d3["feature"].str.replace("gi_", "", regex=False)
    return gi_d3


def load_frozen_d3() -> pd.DataFrame:
    path = (_common.PAPER1_ROOT / "results" / "primary"
           / "mixed_model_results_60s_fwh_EHG9.csv")
    frozen = pd.read_csv(path)
    frozen = frozen.loc[frozen["channel"] == _common.CHANNEL]
    rows = []
    for summary, column in FROZEN_COLUMN.items():
        match = frozen.loc[
            (frozen["feature"] == column)
            & (frozen["contrast"] == "fetal_movement_vs_baseline")]
        if match.empty:
            continue
        record = match.iloc[0]
        rows.append({
            "feature": column, "contrast": "fetal_movement_vs_baseline",
            "estimate": float(record["estimate"]),
            "q_value": float(record["q_value"]),
            "n": int(record["n"]), "n_women": int(record["n_women"]),
            "condition": "frozen_d3", "summary": summary,
        })
    return pd.DataFrame(rows)


def main() -> int:
    features_path = _common.WORK_DIR / "resolution_ablation_features.csv"
    if not features_path.exists():
        raise RuntimeError(
            f"{features_path} not found. Run resolution_ablation.py first.")
    frame = pd.read_csv(features_path)

    print("fitting the 4 new conditions x 4 summaries...", flush=True)
    new_results = fit_new_conditions(frame)
    gi_d3 = load_gi_d3()
    frozen_d3 = load_frozen_d3()

    keep = ["feature", "contrast", "estimate", "q_value", "condition",
            "summary", "n", "n_women"]
    combined = pd.concat(
        [new_results[keep], gi_d3[keep], frozen_d3[keep]], ignore_index=True)
    combined["order"] = combined["condition"].map(CONDITION_ORDER)
    combined = combined.sort_values(["summary", "order"]).reset_index(drop=True)

    print()
    print("=" * 78)
    print("RESOLUTION ABLATION -- FM versus IN, EHG9/FWH")
    print("does significance track INVARIANCE or BIN COUNT?")
    print("=" * 78)
    verdicts = {}
    for summary in SUMMARIES:
        block = combined.loc[combined["summary"] == summary]
        print()
        print(f"-- {summary} --")
        print(f"{'condition':<12}{'invariant?':<12}{'estimate':>10}"
              f"{'q_value':>12}{'sig (q<.05)':>13}")
        gi_sig, frozen_sig = [], []
        for _, r in block.iterrows():
            invariant = r["condition"].startswith("gi")
            sig = bool(pd.notna(r["q_value"]) and r["q_value"] < 0.05)
            (gi_sig if invariant else frozen_sig).append(sig)
            print(f"{r['condition']:<12}{'yes' if invariant else 'no':<12}"
                  f"{r['estimate']:>10.3f}{r['q_value']:>12.4f}"
                  f"{'yes' if sig else 'no':>13}")
        if all(gi_sig) and not any(frozen_sig):
            verdict = "INVARIANCE-DRIVEN: significant at every GI resolution, never without invariance"
        elif not any(gi_sig):
            verdict = "NOT SUPPORTED: not significant at any GI resolution"
        elif gi_sig[-1] and not gi_sig[0]:
            verdict = "RESOLUTION-SENSITIVE: only significant at the finest GI resolution -- treat with caution"
        elif any(frozen_sig):
            verdict = "AMBIGUOUS: resolution alone (without invariance) also reaches significance somewhere"
        else:
            verdict = "MIXED: inspect the table directly"
        verdicts[summary] = verdict
        print(f"  verdict: {verdict}")

    print()
    print("=" * 78)
    print("SUMMARY")
    for summary, verdict in verdicts.items():
        print(f"  {summary:<10} {verdict}")
    print("=" * 78)

    combined.drop(columns=["order"]).to_csv(
        _common.RESULTS_DIR / "resolution_ablation_summary.csv", index=False)
    print()
    print(f"wrote results/resolution_ablation_summary.csv "
          f"({len(combined)} rows) -- send this one back")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
