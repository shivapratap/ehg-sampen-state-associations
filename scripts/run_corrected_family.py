#!/usr/bin/env python3
"""Corrected multiplicity families for the GI cross-channel result (Table F5 / T2 of the manuscript).

Refits the GI Total and GI SD models (state contrasts, woman + recording random effects) in all
four channels and applies Benjamini-Hochberg FDR over
  * the 16-test family: 4 channels x 2 GI summaries x 2 contrasts (UC-IN, FM-IN), and
  * the 2-summary family inside each contrast for EHG9.
Needs statsmodels. Inputs : work/gi_features_60s_EHG9.csv, work/cross_channel_gi_features.csv,
                            results/gi/tables/cross_channel_replication_summary.csv (earlier 8-test q, for comparison)
Outputs: results/manuscript/corrected_family/{corrected_cross_channel_all_16.csv,
         corrected_cross_channel_missing_uc_in.csv, corrected_two_summary_family.csv}
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from gi_sampen.statistics import fit_one  # noqa: E402

OUT = ROOT / "results" / "manuscript" / "corrected_family"
OUT.mkdir(parents=True, exist_ok=True)
gi = pd.read_csv(ROOT / "work" / "gi_features_60s_EHG9.csv")
cc = pd.read_csv(ROOT / "work" / "cross_channel_gi_features.csv")
saved = pd.read_csv(ROOT / "results" / "gi" / "tables" / "cross_channel_replication_summary.csv")

rows = []
frames = [("EHG9", gi)] + [(ch, cc[cc.channel == ch].copy()) for ch in ["EHG10", "EHG11", "EHG12"]]
for ch, frame in frames:
    for feat in ["gi_total", "gi_sd"]:
        for r in fit_one(frame, feat):
            r["channel"] = ch
            rows.append(r)
allr = pd.DataFrame(rows)
allr["q_16test_family"] = np.nan
m = allr.p_value.notna()
allr.loc[m, "q_16test_family"] = false_discovery_control(allr.loc[m, "p_value"].to_numpy(), method="bh")
old = saved[["channel", "feature", "contrast", "q_value"]].rename(columns={"q_value": "saved_8test_q"})
out = allr.merge(old, on=["channel", "feature", "contrast"], how="left")
out.to_csv(OUT / "corrected_cross_channel_all_16.csv", index=False)
out[out.contrast.eq("contraction_vs_baseline")].to_csv(OUT / "corrected_cross_channel_missing_uc_in.csv", index=False)

g = allr[allr.channel.eq("EHG9")][["feature", "contrast", "p_value", "estimate", "std_error",
                                    "ci_low", "ci_high", "channel"]].copy()
g["q_2summary_family"] = np.nan
for c in g.contrast.dropna().unique():
    mm = (g.contrast == c) & g.p_value.notna()
    g.loc[mm, "q_2summary_family"] = false_discovery_control(g.loc[mm, "p_value"].to_numpy(), method="bh")
g.to_csv(OUT / "corrected_two_summary_family.csv", index=False)
print(out[["channel", "feature", "contrast", "estimate", "p_value", "q_16test_family", "saved_8test_q"]]
      .to_string(index=False, float_format=lambda v: f"{v:.4g}"))
