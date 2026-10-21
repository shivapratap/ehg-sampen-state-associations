#!/usr/bin/env python3
"""Stage F2 -- fit EHG10-12, pool with EHG9, report cross-channel replication.

Decision rule, fixed before running
------------------------------------
For each feature (GI-SD, GI-Total), FM-versus-IN contrast:

  SAME SIGN in all four channels (EHG9-12)        -- required
  FDR q < 0.05 in at least 3 of 4 channels          -- "replicated"
  otherwise                                          -- "direction only" or "not replicated"

This mirrors Paper 1's own cross-channel bar (README: "seven highlighted
features retained the same UC-versus-IN direction across EHG9-EHG12 ...
five ... FDR-supported"): same-direction across channels is the baseline
claim, FDR support in most channels is the stronger one.

FDR family
----------
All 16 tests (4 channels x 2 features x 2 contrasts) are corrected together
as ONE family. EHG9's raw p-values are pulled from Stage C's fit (already on
file) and re-corrected here alongside EHG10-12 -- its Stage-C q-value, which
was computed within the larger 13-feature GI family, is NOT reused, so the
cross-channel claim rests on its own clean family rather than borrowing
significance from a differently-scoped correction.
"""
from __future__ import annotations

import pandas as pd

from . import paths as _common
from .statistics import fit_one, apply_fdr

FEATURES = ("gi_sd", "gi_total")
REPLICATION_CHANNELS = ("EHG10", "EHG11", "EHG12")
ALL_CHANNELS = ("EHG9",) + REPLICATION_CHANNELS
FROZEN_COLUMN = {"gi_sd": "fwh_entropy_sdsampen", "gi_total": "fwh_entropy_totalsampen"}


def fit_replication_channels(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for channel in REPLICATION_CHANNELS:
        subset = frame.loc[frame["channel"] == channel]
        for feature in FEATURES:
            for record in fit_one(subset, feature):
                record["channel"] = channel
                rows.append(record)
    return pd.DataFrame(rows)


def load_ehg9_raw_fits() -> pd.DataFrame:
    """EHG9's fit rows from Stage C, with q_value/validated stripped so they
    can be re-corrected inside this stage's own family."""
    path = _common.RESULTS_DIR / "gi_model_results.csv"
    frame = pd.read_csv(path)
    frame = frame.loc[frame["feature"].isin(FEATURES)].copy()
    frame["channel"] = "EHG9"
    return frame.drop(columns=["q_value", "validated", "frozen_estimate",
                               "frozen_q_value"], errors="ignore")


def load_frozen_reference() -> pd.DataFrame:
    rows = []
    for channel in ALL_CHANNELS:
        path = _common.FROZEN_MODEL_RESULTS_PATH[channel]
        frame = pd.read_csv(path)
        for feature, column in FROZEN_COLUMN.items():
            match = frame.loc[
                (frame["feature"] == column)
                & (frame["contrast"] == "fetal_movement_vs_baseline")]
            if match.empty:
                continue
            record = match.iloc[0]
            rows.append({
                "channel": channel, "feature": feature,
                "frozen_estimate": float(record["estimate"]),
                "frozen_q_value": float(record["q_value"]),
            })
    return pd.DataFrame(rows)


def main() -> int:
    features_path = _common.WORK_DIR / "cross_channel_gi_features.csv"
    if not features_path.exists():
        raise RuntimeError(
            f"{features_path} not found. Run cross_channel_replication.py first.")
    frame = pd.read_csv(features_path)

    print("fitting EHG10-12...", flush=True)
    new_fits = fit_replication_channels(frame)
    ehg9_fits = load_ehg9_raw_fits()

    combined = pd.concat([new_fits, ehg9_fits], ignore_index=True)
    combined = apply_fdr(combined)  # single family, 4 channels x 2 features x 2 contrasts
    fm = combined.loc[combined["contrast"] == "fetal_movement_vs_baseline"].copy()

    frozen = load_frozen_reference()
    fm = fm.merge(frozen, on=["channel", "feature"], how="left")

    channel_order = {c: i for i, c in enumerate(ALL_CHANNELS)}
    fm["order"] = fm["channel"].map(channel_order)
    fm = fm.sort_values(["feature", "order"])

    print()
    print("=" * 78)
    print("CROSS-CHANNEL REPLICATION -- FM versus IN, GI-SD and GI-Total")
    print("=" * 78)
    verdicts = {}
    for feature in FEATURES:
        block = fm.loc[fm["feature"] == feature]
        print()
        print(f"-- {feature} --")
        print(f"{'channel':<8}{'GI est':>9}{'GI q':>10}{'sig':>6}"
              f"{'frozen est':>12}{'frozen q':>10}")
        signs, sigs = [], []
        for _, r in block.iterrows():
            sig = bool(pd.notna(r["q_value"]) and r["q_value"] < 0.05)
            signs.append(1 if r["estimate"] > 0 else -1)
            sigs.append(sig)
            print(f"{r['channel']:<8}{r['estimate']:>9.3f}{r['q_value']:>10.4f}"
                  f"{'yes' if sig else 'no':>6}"
                  f"{r['frozen_estimate']:>12.3f}{r['frozen_q_value']:>10.4f}")
        same_sign = len(set(signs)) == 1
        n_sig = sum(sigs)
        if same_sign and n_sig >= 3:
            verdict = f"REPLICATED -- same sign in all 4 channels, FDR-supported in {n_sig}/4"
        elif same_sign:
            verdict = f"DIRECTION ONLY -- same sign in all 4 channels, FDR-supported in only {n_sig}/4"
        else:
            verdict = "NOT REPLICATED -- sign is not consistent across channels"
        verdicts[feature] = verdict
        print(f"  verdict: {verdict}")

    print()
    print("=" * 78)
    print("SUMMARY")
    for feature, verdict in verdicts.items():
        print(f"  {feature:<10} {verdict}")
    print("=" * 78)

    fm.drop(columns=["order"]).to_csv(
        _common.RESULTS_DIR / "cross_channel_replication_summary.csv", index=False)
    print()
    print(f"wrote results/cross_channel_replication_summary.csv "
          f"({len(fm)} rows) -- send this one back")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
