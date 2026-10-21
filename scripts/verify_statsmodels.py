#!/usr/bin/env python3
"""Independent statsmodels check of the custom-REML EHG9 fits (frozen/GI Total, with and without log-RMS).

Run from anywhere (needs statsmodels):
    python scripts/verify_statsmodels.py

Reads  results/conditional/windows_gi_rms_frozen_full.csv   (4,083 windows, saved)
       results/conditional/refit_gi_total_rms_results.csv   (custom-REML results, saved)
Writes results/conditional/statsmodels_verification.{csv,md}
Expected: OVERALL: ALL WITHIN TOLERANCE (max coefficient difference ~3e-6, max SE difference ~3e-3).

Model (same as the custom fit): outcome standardised over all 4,083 windows (ddof=1);
fixed effects: intercept (IN), UC, FM [, standardised log RMS];
random effects: woman intercept + recording-within-woman variance component; REML.
"""
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import statsmodels
import statsmodels.formula.api as smf

TOL_COEF = 5e-3   # absolute tolerance on coefficients
TOL_SE = 5e-3     # absolute tolerance on standard errors

COND = Path(__file__).resolve().parents[1] / "results" / "conditional"
df = pd.read_csv(COND / "windows_gi_rms_frozen_full.csv")
df = df[df.window_class.isin(["baseline", "contraction", "fetal_movement"])].copy()
df["uc"] = (df.window_class == "contraction").astype(float)
df["fm"] = (df.window_class == "fetal_movement").astype(float)
lr = np.log(df.rms)
df["lrms"] = (lr - lr.mean()) / lr.std(ddof=1)
print("n", len(df), "women", df.patient_id.nunique(), "recordings", df.recording_id.nunique(),
      df.window_class.value_counts().to_dict())
print("statsmodels", statsmodels.__version__)

saved = pd.read_csv(COND / "refit_gi_total_rms_results.csv")
rows = []
for outcome in ["frozen_total", "gi_total"]:
    df["z"] = (df[outcome] - df[outcome].mean()) / df[outcome].std(ddof=1)
    for adj in [False, True]:
        formula = "z ~ uc + fm" + (" + lrms" if adj else "")
        best = None
        # try several optimisers / starts; keep the best REML log-likelihood
        for method in (["lbfgs"], ["bfgs"], ["powell"], ["lbfgs", "powell"]):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    m = smf.mixedlm(formula, df, groups=df["patient_id"],
                                    re_formula="1",
                                    vc_formula={"recording": "0 + C(recording_id)"})
                    r = m.fit(reml=True, method=method, maxiter=500)
                if best is None or r.llf > best.llf:
                    best = r
            except Exception as e:
                print("  optimiser", method, "failed:", e)
        r = best
        names = {"uc": "UC_minus_IN", "fm": "FM_minus_IN", "lrms": "log_RMS_std"}
        for term, tname in names.items():
            if term not in r.params.index:
                continue
            s = saved[(saved.outcome == outcome) & (saved.adjusted_log_RMS == adj) & (saved.term == tname)].iloc[0]
            dc = abs(r.params[term] - s.estimate)
            ds = abs(r.bse[term] - s.se)
            rows.append(dict(outcome=outcome, adjusted_log_RMS=adj, term=tname,
                             statsmodels_estimate=r.params[term], custom_estimate=s.estimate, abs_diff_estimate=dc,
                             statsmodels_se=r.bse[term], custom_se=s.se, abs_diff_se=ds,
                             statsmodels_ci_low=r.conf_int().loc[term, 0], statsmodels_ci_high=r.conf_int().loc[term, 1],
                             sm_var_woman=float(r.cov_re.iloc[0, 0]), sm_var_recording=float(r.vcomp[0]), sm_var_resid=float(r.scale),
                             custom_var_woman=s.var_woman, custom_var_recording=s.var_recording, custom_var_resid=s.var_resid,
                             sm_converged=bool(r.converged),
                             pass_=bool(dc <= TOL_COEF and ds <= TOL_SE)))

out = pd.DataFrame(rows).rename(columns={"pass_": "within_tolerance"})
out.to_csv(COND / "statsmodels_verification.csv", index=False)
ok = bool(out.within_tolerance.all())
with open(COND / "statsmodels_verification.md", "w") as fh:
    fh.write(f"# statsmodels cross-check of the EHG9 conditional fits\n\n"
             f"statsmodels {statsmodels.__version__}; python {sys.version.split()[0]}; MixedLM, REML, woman intercept + recording variance component.\n"
             f"Tolerance: |coefficient difference| <= {TOL_COEF}, |SE difference| <= {TOL_SE}.\n\n"
             f"**Overall: {'ALL WITHIN TOLERANCE' if ok else 'DIFFERENCES FOUND - inspect statsmodels_verification.csv'}**\n\n"
             + "```\n" + out[["outcome", "adjusted_log_RMS", "term", "statsmodels_estimate", "custom_estimate", "abs_diff_estimate",
                    "statsmodels_se", "custom_se", "abs_diff_se", "sm_converged", "within_tolerance"]]
             .to_string(index=False, float_format=lambda v: f"{v:.5f}") + "\n```\n")
print(out[["outcome", "adjusted_log_RMS", "term", "statsmodels_estimate", "custom_estimate",
           "abs_diff_estimate", "statsmodels_se", "custom_se", "abs_diff_se", "sm_converged", "within_tolerance"]]
      .to_string(float_format=lambda v: f"{v:.5f}"))
print("\nOVERALL:", "ALL WITHIN TOLERANCE" if ok else "DIFFERENCES FOUND")
