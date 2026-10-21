#!/usr/bin/env python3
"""Lag-1 residual autocorrelation diagnostic (EHG9) -> Supplementary Table S7A.

For the EHG9 primary-validated features, refit the primary mixed model
(state contrast, woman + recording-within-woman random effects), take the
residuals, and for the 5 recordings with the most windows compute the lag-1
autocorrelation of residuals ordered by window start time.

Needs statsmodels. Inputs : work/primary/event_unit_analysis_60s_EHG9.csv
                            results/primary/mixed_model_results_60s_fwh_EHG9.csv
Outputs: results/robustness/lag1/lag1_residual_autocorrelation.csv
         results/robustness/lag1/lag1_report.md, fig_lag1_residual_order_<rec>.png
Expected: 25 of 55 combinations with |r1| >= 0.30 for the original-11 subset, max |r1| 0.560.
"""
from __future__ import annotations
import sys, warnings
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import config  # noqa: E402
import models  # noqa: E402

CHANNEL, N_TOP, FLAG = "EHG9", 5, 0.30
OUT = config.ROBUSTNESS_DIR / "lag1"


def fit_with_residuals(units, feature):
    d = units[[feature, "window_class", "patient_id", "recording_id", "window_start_s"]].copy()
    d = d.replace([np.inf, -np.inf], np.nan).dropna()
    d["window_class"] = pd.Categorical(d["window_class"], categories=["baseline", "fetal_movement", "contraction"])
    d = d.dropna(subset=["window_class"]).reset_index(drop=True)
    sd = float(d[feature].std())
    if d.empty or not np.isfinite(sd) or sd == 0:
        return None, None
    d["_z"] = (d[feature] - d[feature].mean()) / sd
    formula = '_z ~ C(window_class, Treatment(reference="baseline"))'
    warns: list = []
    att = models._try_tier(formula, d, True, "woman_recording", warns) \
        or models._try_tier(formula, d, False, "woman_only", warns)
    return (att.fit, d) if att is not None else (None, None)


def lag1(x):
    if len(x) < 3 or np.std(x[:-1]) == 0 or np.std(x[1:]) == 0:
        return float("nan")
    return float(np.corrcoef(x[:-1], x[1:])[0, 1])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res = pd.read_csv(config.PRIMARY_DIR / f"mixed_model_results_60s_fwh_{CHANNEL}.csv")
    feats = res.loc[res.contrast.eq("contraction_vs_baseline") & res.primary_validated, "feature"].unique().tolist()
    units = pd.read_csv(config.UNITS_DIR / f"event_unit_analysis_60s_{CHANNEL}.csv")
    top = units.groupby("recording_id").size().sort_values(ascending=False).head(N_TOP).index.tolist()
    print(f"{len(feats)} primary-validated features; top recordings: {top}")

    rows, cache = [], {}
    for f in feats:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit, d = fit_with_residuals(units, f)
        if fit is None:
            continue
        try:
            resid, rtype = fit.resid.to_numpy(), "conditional"
        except Exception:  # singular random-effects covariance: use marginal residuals
            resid = np.asarray(fit.model.endog) - np.asarray(fit.model.exog) @ np.asarray(fit.fe_params)
            rtype = "marginal_fallback"
        d = d.assign(_res=resid)
        cache[f] = d
        for rec in top:
            g = d[d.recording_id == rec].sort_values("window_start_s")
            if len(g) >= 3:
                rows.append(dict(feature=f, recording_id=rec, n_windows=len(g),
                                 lag1_autocorr=lag1(g["_res"].to_numpy()), residual_type=rtype))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "lag1_residual_autocorrelation.csv", index=False)

    v = df.dropna(subset=["lag1_autocorr"])
    flagged = int((v.lag1_autocorr.abs() >= FLAG).sum())
    summary = (f"{len(v)} feature x recording combinations; {flagged} with |r1|>={FLAG}; "
               f"max |r1| = {v.lag1_autocorr.abs().max():.3f}; mean |r1| = {v.lag1_autocorr.abs().mean():.3f}; "
               f"{int((v.residual_type == 'marginal_fallback').sum())} rows use marginal residuals.")
    (OUT / "lag1_report.md").write_text(
        "# Lag-1 residual autocorrelation, EHG9\n\n" + summary + "\n\n```\n"
        + v.sort_values(["feature", "recording_id"]).to_string(index=False, float_format=lambda x: f"{x:.3f}") + "\n```\n")
    # Subset reported in the manuscript: the 11 features FDR-supported in the original 34-feature family.
    f11 = config.ROBUSTNESS_DIR / "dependence" / "original_11_stability_summary.csv"
    if f11.exists():
        o11 = set(pd.read_csv(f11)["feature_name"])
        s11 = v[v.feature.isin(o11)]
        sub = (f"Original-11 subset: {len(s11)} combinations; {int((s11.lag1_autocorr.abs() >= FLAG).sum())} "
               f"with |r1|>={FLAG}; max |r1| = {s11.lag1_autocorr.abs().max():.3f}; "
               f"{int((s11.residual_type == 'marginal_fallback').sum())} use marginal residuals.")
        summary += "\n\n" + sub
        with open(OUT / "lag1_report.md", "a") as fh:
            fh.write("\n" + sub + "\n")
    print(summary)

    plot_f = [f for f in feats if any(t in f for t in ("totalsampen", "root_mean_square", "spectral_edge"))][:3] or feats[:3]
    fig, axes = plt.subplots(len(plot_f), 1, figsize=(7, 2.3 * len(plot_f)), sharex=True, squeeze=False)
    for ax, f in zip(axes[:, 0], plot_f):
        g = cache[f][cache[f].recording_id == top[0]].sort_values("window_start_s")
        ax.plot(range(len(g)), g["_res"], marker="o", ms=3, lw=0.8)
        ax.axhline(0, color="grey", lw=0.7, ls="--")
        ax.set_title(f, fontsize=9)
    axes[-1, 0].set_xlabel("within-recording window order")
    fig.tight_layout()
    fig.savefig(OUT / f"fig_lag1_residual_order_{top[0]}.png", dpi=200)


if __name__ == "__main__":
    main()
