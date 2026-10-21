#!/usr/bin/env python3
"""Compare the key numbers in results/ and work/ with the values reported in the manuscript.

    python scripts/check_checkpoints.py

Needs only pandas. Run it after restoring the derived data (and again after re-running any stage):
a PASS means the file on disk holds the number the paper reports; it does not re-fit anything.
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
rows = []


def check(name, ok, detail=""):
    rows.append((bool(ok), name, detail))


def rd(p):
    return pd.read_csv(ROOT / p)


def near(a, b, tol):
    return abs(float(a) - float(b)) <= tol


try:
    u = rd("work/primary/event_unit_analysis_60s_EHG9.csv")
    vc = u.window_class.value_counts()
    check("analysed windows = 4,083", len(u) == 4083, len(u))
    check("IN / UC / FM = 3,358 / 333 / 392", (vc["baseline"], vc["contraction"], vc["fetal_movement"]) == (3358, 333, 392), {k: int(v) for k, v in vc.items()})
    check("42 women, 112 recordings", (u.patient_id.nunique(), u.recording_id.nunique()) == (42, 112))
    s = rd("work/primary/mixed_model_summary_60s_fwh_all_channels.csv").query("channel == 'EHG9'").iloc[0]
    check("EHG9 UC-vs-IN FDR-supported features = 11 of 34", (s.contraction_vs_baseline_validated, s.features_requested) == (11, 34))
    m = rd("work/windows/window_manifest_60s.csv")
    check("candidate windows = 4,283", len(m) == 4283, len(m))
except FileNotFoundError as e:
    check("derived data present (run scripts/restore_derived.py)", False, e.filename)

try:
    g = rd("results/manuscript/corrected_family/corrected_cross_channel_all_16.csv")
    r = g[(g.contrast == "contraction_vs_baseline") & (g.channel == "EHG9") & (g.feature == "gi_total")].iloc[0]
    check("GI Total EHG9 UC-IN = -0.113", near(r.estimate, -0.113, 0.0005), round(r.estimate, 4))
    check("GI Total EHG9 q (16-test family) = 0.021", near(r.q_16test_family, 0.0206, 0.0005), round(r.q_16test_family, 4))
    n = int(((g.contrast == "contraction_vs_baseline") & (g.q_16test_family < 0.05)).sum())
    check("UC-IN tests with q<0.05 in the 16-test family = 3 (EHG9 Total and SD, EHG10 SD)", n == 3, n)
except FileNotFoundError as e:
    check("corrected family present", False, e.filename)

c = rd("results/conditional/refit_gi_total_rms_results.csv")
uc = c[c.term == "UC_minus_IN"].set_index(["outcome", "adjusted_log_RMS"]).estimate
check("GI Total UC: -0.113 -> +0.275 after log RMS", near(uc[("gi_total", False)], -0.113, 5e-4) and near(uc[("gi_total", True)], 0.275, 5e-4),
      f"{uc[('gi_total', False)]:+.3f} -> {uc[('gi_total', True)]:+.3f}")
check("Frozen Total UC: 0.409 -> 0.020 after log RMS", near(uc[("frozen_total", False)], 0.409, 5e-4) and near(uc[("frozen_total", True)], 0.020, 5e-4),
      f"{uc[('frozen_total', False)]:+.3f} -> {uc[('frozen_total', True)]:+.3f}")
lw = rd("results/conditional/lowo_results.csv")
check("leave-one-woman-out: 42 refits", len(lw) == 42, len(lw))
check("LOWO GI Total unadjusted always negative, adjusted always positive",
      (lw.gi_total_unadj_UC < 0).all() and (lw.gi_total_adj_UC > 0).all(),
      f"{lw.gi_total_unadj_UC.min():+.3f}..{lw.gi_total_unadj_UC.max():+.3f}; {lw.gi_total_adj_UC.min():+.3f}..{lw.gi_total_adj_UC.max():+.3f}")
sv = rd("results/conditional/statsmodels_verification.csv")
check("statsmodels cross-check within tolerance (coef <= 5e-3, SE <= 5e-3)", sv.within_tolerance.all(),
      f"max |coef diff| {sv.abs_diff_estimate.max():.1e}, max |SE diff| {sv.abs_diff_se.max():.1e}")

lag = rd("results/robustness/lag1/lag1_residual_autocorrelation.csv")
o11 = set(rd("results/robustness/dependence/original_11_stability_summary.csv").feature_name)
l11 = lag[lag.feature.isin(o11)]
check("lag-1: 55 combinations, 25 with |r|>=0.30, max 0.56 (original-11 features)",
      len(l11) == 55 and int((l11.lag1_autocorr.abs() >= 0.3).sum()) == 25 and near(l11.lag1_autocorr.abs().max(), 0.560, 5e-4),
      f"{len(l11)}, {int((l11.lag1_autocorr.abs() >= 0.3).sum())}, {l11.lag1_autocorr.abs().max():.3f}")
chk = (ROOT / "results/manuscript/supplement/Supplement_S7_S8_checks.txt").read_text().splitlines()
check("Supplement S7/S8: every quoted number matches its source", sum(l.startswith("PASS") for l in chk) == 36 and not any(l.startswith("FAIL") for l in chk),
      f"{sum(l.startswith('PASS') for l in chk)} PASS, {sum(l.startswith('FAIL') for l in chk)} FAIL")

for ok, name, detail in rows:
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")
print(f"\n{sum(o for o, *_ in rows)}/{len(rows)} checkpoints pass")
sys.exit(0 if all(o for o, *_ in rows) else 1)
