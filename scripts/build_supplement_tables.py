#!/usr/bin/env python3
"""Build Supplementary Tables S7 and S8 (LaTeX) from saved outputs only; no model is fitted.

Usage: python scripts/build_supplement_tables.py
Reads : results/reviewer_revision/*, results/robustness/lag1/, results/entropy/, results/conditional/
Writes: results/manuscript/supplement/Supplement_S7_S8.tex  (\\input into the supplement)
        results/manuscript/supplement/Supplement_S7_S8_checks.txt  (every number quoted in the text, PASS/FAIL)
Pure standard library.
"""
import csv, os, sys, statistics as st
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RR = str(ROOT / "results" / "reviewer_revision")
LAG = str(ROOT / "results" / "robustness" / "lag1" / "lag1_residual_autocorrelation.csv")
ENT = str(ROOT / "results" / "entropy")
COND = str(ROOT / "results" / "conditional")
OUT = ROOT / "results" / "manuscript" / "supplement"
OUT.mkdir(parents=True, exist_ok=True)
OUT = str(OUT)


def rd(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def tex(s):
    return str(s).replace("\\", "\\textbackslash{}").replace("_", "\\_").replace("%", "\\%").replace("&", "\\&")


def f(x, d=3):
    x = float(x)
    s = f"{x:.{d}f}"
    return s.replace("-", "$-$") if s.startswith("-") else s


def q(x):
    x = float(x)
    return f"{x:.1e}" if x < 0.001 else f"{x:.3f}"


def ci(lo, hi, d=3):
    return f"[{f(lo, d)}, {f(hi, d)}]"


checks = []


def chk(name, ok, detail=""):
    checks.append(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")


# ---------------------------------------------------------------- sources
t01 = {r["feature_name"]: r for r in rd(f"{RR}/01_temporal_dependence/original_11_stability_summary.csv")}
t02 = {r["feature"]: r for r in rd(f"{RR}/02_patient_level_uc_in/original_11_patient_level_stability.csv")}
t05 = {r["feature"]: r for r in rd(f"{RR}/05_balanced_in_subsampling/original_11_balanced_in_stability.csv")}
t06 = {r["feature"]: r for r in rd(f"{RR}/06_window_length_sensitivity/original_11_window_length_stability.csv")}
t07 = {r["feature"]: r for r in rd(f"{RR}/07_overlap_contamination_sensitivity/original_11_clean_uc_stability.csv")}
t08 = {r["feature"]: r for r in rd(f"{RR}/08_residual_distribution_robustness/original_11_residual_robustness.csv")}
lag = rd(LAG)
feats = list(t01.keys())
assert len(feats) == 11
label = {k: t05[k]["short_label"] for k in feats}
recs = []
for r in lag:
    if r["recording_id"] not in recs:
        recs.append(r["recording_id"])

# ---------------------------------------------------------------- checks vs. manuscript text
lag11 = [r for r in lag if r["feature"] in feats]
v = [abs(float(r["lag1_autocorr"])) for r in lag11]
n_flag = sum(a >= 0.30 for a in v)
n_marg = sum(r["residual_type"] != "conditional" for r in lag11)
_mm = rd(str(ROOT / "results" / "primary" / "mixed_model_results_60s_fwh_EHG9.csv"))
zero_feats = {r["feature"] for r in _mm if r.get("model") == "woman_only" and r.get("boundary_fit") == "True"}
n_zero = sum(r["feature"] in zero_feats for r in lag11)
chk("lag-1: 55 feature-recording combinations (11 x 5)", len(lag11) == 55, f"n={len(lag11)}")
chk("lag-1: |r1|>=0.30 count", True, f"{n_flag} of 55 (manuscript: 25)")
chk("lag-1: maximum |r1|", True, f"{max(v):.3f} (manuscript: 0.560)")
chk("lag-1: rows whose fit has zero random-effect variance (conditional = marginal residuals)", True, f"{n_zero} of 55; labelled marginal in saved file: {n_marg}")

fdr60 = sum(float(t06[k]["q_60s"]) < 0.05 for k in feats)
fdr80 = sum(float(t06[k]["q_80s"]) < 0.05 for k in feats)
fdr40 = sum(float(t06[k]["q_40s"]) < 0.05 for k in feats)
chk("duration FDR support 60 s / 80 s / 40 s", (fdr60, fdr80, fdr40) == (11, 11, 10), f"{fdr60}/{fdr80}/{fdr40}")
lost40 = [label[k] for k in feats if float(t06[k]["q_40s"]) >= 0.05]
chk("duration: feature losing support at 40 s", lost40 == ["Hjorth complexity"] or "Hjorth" in lost40[0], str(lost40))
chk("duration: signs agree in all 11", all(t06[k]["all_three_signs_agree"] == "True" for k in feats))

clean = sum(float(t07[k]["clean_q_value"]) < 0.05 for k in feats)
lost_clean = {label[k]: round(float(t07[k]["clean_q_value"]), 3) for k in feats if float(t07[k]["clean_q_value"]) >= 0.05}
chk("clean UC (312 windows): FDR support", clean == 9, f"{clean} of 11; lost: {lost_clean}")
chk("clean UC: 312 UC windows", all(int(t07[k]["clean_n_contraction"]) == 312 for k in feats))

chk("balanced IN: all 11 same direction in all resamples", all(float(t05[k]["same_direction_frequency"]) == 1.0 for k in feats))
chk("balanced IN: none >=80% BH-supported", all(float(t05[k]["bh_q_lt_0_05_frequency"]) < 0.8 for k in feats),
    str({label[k]: float(t05[k]["bh_q_lt_0_05_frequency"]) for k in feats}))
amp = [k for k in feats if label[k] in ("MAV", "RMS", "Peak-to-peak", "Waveform length", "Total SampEn")]
oth = [k for k in feats if k not in amp]
ra = [float(t05[k]["median_se_ratio_vs_frozen"]) for k in amp]
ro = [float(t05[k]["median_se_ratio_vs_frozen"]) for k in oth]
chk("balanced IN: amplitude-group median SE ratio range", True, f"{min(ra):.2f}-{max(ra):.2f} (text 3.21-3.36); n={len(amp)}")
chk("balanced IN: other-group median SE ratio range", True, f"{min(ro):.2f}-{max(ro):.2f} (text 1.26-1.41); n={len(oth)}")
chk("CR1: same direction", True, f"{sum(t01[k]['woman_cluster_same_direction_as_primary']=='True' for k in feats)} of 11")
chk("CR1: q<0.05 among the 11", True, f"{sum(float(t01[k]['woman_cluster_q'])<0.05 for k in feats)} of 11")
chk("woman-level: q<0.05 among the 11", True, f"{sum(float(t02[k]['patient_level_q'])<0.05 for k in feats)} of 11 (text 5 of 11)")
chk("woman-level: same direction as mixed model", True, f"{sum(t02[k]['direction_agreement_with_mixed_model']=='True' for k in feats)} of 11")

# ---------------------------------------------------------------- S7 LaTeX
L = []
A = L.append
A("% Generated by build_S7_S8.py from saved outputs; do not edit by hand.")
A("\\clearpage\n\\begin{landscape}")
A("\\section*{Table S7. Robustness and residual diagnostics for the 11 primary EHG9 UC--IN features}\\label{tab:s7_robustness}")
A("All values are read from saved analysis outputs; no model was refitted for this document. The 11 features are those FDR-supported in the original EHG9 34-feature family. Coefficients are standardised UC-minus-IN estimates. Panels S7A--S7D are descriptive and are not a second inferential family.")

# --- S7A lag-1
A("\\subsection*{S7A. Lag-1 autocorrelation of model residuals in the five most-windowed EHG9 recordings}")
A(f"Residuals were ordered by window start time within each recording. Entries are Pearson lag-1 autocorrelation; bold marks $|r_1|\\ge0.30$. Across the 11 features and five recordings, {n_flag} of 55 combinations reached the informal $0.30$ flag, with a maximum $|r_1|$ of {max(v):.3f}; for {n_zero} of 55 combinations (the four amplitude features whose primary fit was a woman-only boundary fit) the random-effect variance was estimated at zero, so conditional and marginal residuals coincide. This is a diagnostic, not a model of serial correlation.")
A("\\begin{longtable}{p{3.6cm}" + "r" * len(recs) + "r}")
A("\\toprule")
hdr = "\\textbf{Feature} & " + " & ".join(f"\\textbf{{{tex(r)}}}" for r in recs) + " & \\textbf{n flagged} \\\\"
A(hdr + "\n\\midrule\n\\endfirsthead\n\\toprule\n" + hdr + "\n\\midrule\n\\endhead")
for k in feats:
    cells, nf = [], 0
    for rc in recs:
        row = [r for r in lag if r["feature"] == k and r["recording_id"] == rc][0]
        x = float(row["lag1_autocorr"])
        s = f(x)
        if row["residual_type"] != "conditional":
            s += "*"
        if abs(x) >= 0.30:
            s = f"\\textbf{{{s}}}"
            nf += 1
        cells.append(s)
    A(f"{label[k]} & " + " & ".join(cells) + f" & {nf} \\\\")
nw = {rc: [r for r in lag if r['recording_id'] == rc][0]['n_windows'] for rc in recs}
A("\\midrule\n\\textit{Windows per recording} & " + " & ".join(nw[rc] for rc in recs) + " & \\\\")
A("\\bottomrule\n\\end{longtable}")

# --- S7B model dependence
A("\\subsection*{S7B. Model dependence: primary mixed model, pooled OLS, woman-cluster CR1 and woman-level paired medians}")
A("q-values are Benjamini--Hochberg within the original 34-feature family for each estimator. Woman-level columns give the median paired UC-minus-IN difference in raw feature units (34 women), the Wilcoxon-based q-value and the proportion of women whose difference has the mixed-model sign.")
A("\\begin{longtable}{p{3.0cm}rrrrrrrrrr}")
A("\\toprule")
h = "\\textbf{Feature} & \\textbf{Mixed $\\beta$} & \\textbf{q} & \\textbf{OLS $\\beta$} & \\textbf{OLS q} & \\textbf{CR1 $\\beta$} & \\textbf{CR1 q} & \\textbf{Median diff.} & \\textbf{Women q} & \\textbf{Women, same sign (\\%)} \\\\"
A(h + "\n\\midrule\n\\endfirsthead\n\\toprule\n" + h + "\n\\midrule\n\\endhead")
for k in feats:
    a, b = t01[k], t02[k]
    A(f"{label[k]} & {f(a['primary_mixed_model_estimate'])} & {q(a['primary_mixed_model_q'])} & {f(a['pooled_ols_estimate'])} & {q(a['pooled_ols_q'])} & {f(a['woman_cluster_estimate'])} & {q(a['woman_cluster_q'])} & {float(b['median_paired_difference']):.2e} & {q(b['patient_level_q'])} & {float(b['pct_mixed_model_consistent']):.0f} \\\\")
A("\\bottomrule\n\\end{longtable}")

# --- S7C balanced / duration / clean UC
A("\\subsection*{S7C. Balanced IN subsampling, window duration and clean-UC sensitivity}")
A("Balanced IN: 1,000 within-woman resamples of IN windows (same-direction and BH-supported frequencies, and the median ratio of the resampled to the primary-fit standard error). Duration: window length 40, 60 and 80 s with the same event centres. Clean UC: 21 UC windows containing a competing FM annotation removed (312 UC windows).")
A("\\begin{longtable}{p{3.0cm}rrrrrrrrrr}")
A("\\toprule")
h = "\\textbf{Feature} & \\textbf{Same dir.} & \\textbf{q$<$0.05 freq.} & \\textbf{SE ratio} & \\textbf{40 s $\\beta$} & \\textbf{40 s q} & \\textbf{80 s $\\beta$} & \\textbf{80 s q} & \\textbf{Clean $\\beta$} & \\textbf{Clean q} \\\\"
A(h + "\n\\midrule\n\\endfirsthead\n\\toprule\n" + h + "\n\\midrule\n\\endhead")
for k in feats:
    b, d, c = t05[k], t06[k], t07[k]
    A(f"{label[k]} & {float(b['same_direction_frequency']):.3f} & {float(b['bh_q_lt_0_05_frequency']):.3f} & {float(b['median_se_ratio_vs_frozen']):.2f} & {f(d['coefficient_40s'])} & {q(d['q_40s'])} & {f(d['coefficient_80s'])} & {q(d['q_80s'])} & {f(c['clean_estimate'])} & {q(c['clean_q_value'])} \\\\")
A("\\bottomrule\n\\end{longtable}")

# --- S7D residual shape
A("\\subsection*{S7D. Conditional-residual distribution and spread}")
A("Descriptive only (no formal test). Spread ratio is the largest ratio of residual MAD between states; the last column is the Spearman correlation of absolute conditional residual with fitted value.")
A("\\begin{longtable}{p{3.0cm}rrrrr}")
A("\\toprule")
h = "\\textbf{Feature} & \\textbf{Skewness} & \\textbf{Excess kurtosis} & \\textbf{Max/min SD ratio} & \\textbf{Spread ratio} & \\textbf{Spearman $|e|$ vs fit} \\\\"
A(h + "\n\\midrule\n\\endfirsthead\n\\toprule\n" + h + "\n\\midrule\n\\endhead")
for k in feats:
    r = t08[k]
    A(f"{label[k]} & {f(r['conditional_residual_skewness'], 2)} & {f(r['conditional_residual_excess_kurtosis'], 2)} & {float(r['max_min_residual_sd_ratio']):.2f} & {float(r['largest_class_spread_ratio']):.2f} & {f(r['spearman_abs_residual_vs_fitted_rho'], 2)} \\\\")
A("\\bottomrule\n\\end{longtable}\n\\end{landscape}")

# ---------------------------------------------------------------- S8
sc = rd(f"{RR}/04_sampen_scale_dependence/phase1_scale_audit/scale_sensitivity_summary.csv")
adj = rd(f"{RR}/04_sampen_scale_dependence/phase2_rms_adjustment/corrected_amplitude_adjusted_model_results.csv")
dia = rd(f"{RR}/04_sampen_scale_dependence/phase2_rms_adjustment/entropy_rms_residual_diagnostics.csv")
mech = rd(f"{RR}/04_sampen_scale_dependence/phase1_scale_audit/scale_sensitivity_mechanism.csv")
wl = rd(f"{ENT}/module3_5_sampen_amplitude_summary_table_all_channels.csv")
gi = rd(f"{COND}/refit_gi_total_rms_results.csv")
lowo = rd(f"{COND}/lowo_results.csv")

A("\\clearpage\n\\begin{landscape}")
A("\\section*{Table S8. Amplitude dependence of the entropy associations}\\label{tab:s8_amplitude}")
A("All values are read from saved analysis outputs; no model was refitted for this document.")

# S8A scaling
A("\\subsection*{S8A. Direct positive scaling of 60 windows (20 per state)}")
A("Median relative change from unit scale (\\%). ``Non-decr.'' is the number of the 60 trajectories that were non-decreasing with scale; ``Invariant'' is the number numerically unchanged at every scale.")
A("\\begin{longtable}{p{3.2cm}rrrrrrr}")
A("\\toprule")
h = "\\textbf{Summary} & \\textbf{0.25$\\times$} & \\textbf{0.5$\\times$} & \\textbf{2$\\times$} & \\textbf{4$\\times$} & \\textbf{Non-decr.} & \\textbf{Invariant} \\\\"
A(h + "\n\\midrule\n\\endfirsthead\n\\toprule\n" + h + "\n\\midrule\n\\endhead")
names = {"rms": "RMS", "total_sampen": "Total SampEn", "average_sampen": "Average SampEn", "median_sampen": "Median SampEn", "fixed_r_sampen": "Fixed-$r$ SampEn"}
tot = {}
for key, nm in names.items():
    rows = {float(r["scale_factor"]): r for r in sc if r["stratum"] == "overall" and r["feature"] == key}
    nd = rows[4.0]["n_trajectory_nondecreasing"]
    inv = rows[4.0]["n_trajectory_effectively_invariant"]
    tot[key] = {s: float(rows[s]["median_relative_change_percent"]) for s in rows}
    A(f"{nm} & " + " & ".join(f"{tot[key][s]:+.1f}" for s in (0.25, 0.5, 2.0, 4.0)) + f" & {nd} & {inv} \\\\")
A("\\bottomrule\n\\end{longtable}")
chk("scaling: Total SampEn -75/-50/+102/+309 at 0.25/0.5/2/4x",
    [round(tot['total_sampen'][s]) for s in (0.25, 0.5, 2.0, 4.0)] == [-75, -50, 102, 309],
    str({s: round(v, 1) for s, v in tot['total_sampen'].items()}))
chk("scaling: average/median SampEn at 4x = +4.6 / +8.8 %", round(tot['average_sampen'][4.0], 1) == 4.6 and round(tot['median_sampen'][4.0], 1) == 8.8,
    f"{tot['average_sampen'][4.0]:.2f} / {tot['median_sampen'][4.0]:.2f}")
chk("scaling: fixed-r SampEn invariant", all(abs(v) < 1e-6 for v in tot['fixed_r_sampen'].values()))
pooled_rho = [r for r in mech if r[list(r.keys())[0]] == "overall" and "n_retained_profile_values" in "".join(r.values())]
chk("scaling: Spearman between change in Total SampEn and change in retained profile values (overall)", True,
    ",".join(f"{float(r[list(r.keys())[7]]):.3f}" for r in pooled_rho) + "  [column 8 of scale_sensitivity_mechanism.csv; text 0.966]")

# S8B adjusted models
A("\\subsection*{S8B. Four-channel UC/IN models with and without standardised log-RMS}")
A("Pooled UC/IN windows (3,691 physical windows, 14,764 channel rows; 42 women; EHG9 reference channel). Mixed: woman random intercept plus physical-window variance component, maximum likelihood. CR1: fixed-effects OLS with woman-cluster CR1 covariance. Outcomes and log-RMS standardised over the pooled UC/IN population.")
A("\\begin{longtable}{p{2.8cm}p{3.2cm}p{1.8cm}rrrr}")
A("\\toprule")
h = "\\textbf{Outcome} & \\textbf{Estimator} & \\textbf{Model} & \\textbf{UC $\\beta$} & \\textbf{95\\% CI} & \\textbf{p} & \\textbf{log-RMS $\\beta$} \\\\"
A(h + "\n\\midrule\n\\endfirsthead\n\\toprule\n" + h + "\n\\midrule\n\\endhead")
est_name = {"preferred_mixed_model": "Mixed (ML)", "woman_cluster_robust_OLS": "OLS, woman-cluster CR1"}
vals = {}
for r in adj:
    z = r["z_log_rms_estimate"]
    vals[(r["feature"], r["method"], r["model_variant"])] = float(r["uc_vs_in_estimate"])
    A(f"{tex(r['feature_label'])} & {est_name.get(r['method'], tex(r['method']))} & {r['model_variant']} & {f(r['uc_vs_in_estimate'])} & {ci(r['uc_vs_in_ci_low'], r['uc_vs_in_ci_high'])} & {q(r['uc_vs_in_p_value']) if float(r['uc_vs_in_p_value'])>0 else '$<$1e-15'} & {f(z) if z else '--'} \\\\")
A("\\bottomrule\n\\end{longtable}")
expect = {
    ("fwh_entropy_totalsampen", "preferred_mixed_model", "unadjusted"): 0.299,
    ("fwh_entropy_totalsampen", "preferred_mixed_model", "adjusted"): -0.086,
    ("fwh_entropy_totalsampen", "woman_cluster_robust_OLS", "adjusted"): 0.030,
    ("fwh_entropy_sampen_fixed_r", "preferred_mixed_model", "unadjusted"): -0.213,
    ("fwh_entropy_sampen_fixed_r", "preferred_mixed_model", "adjusted"): 0.235,
    ("fwh_entropy_sampen_fixed_r", "woman_cluster_robust_OLS", "adjusted"): 0.083,
    ("fwh_entropy_svd_entropy", "preferred_mixed_model", "unadjusted"): -0.053,
    ("fwh_entropy_svd_entropy", "preferred_mixed_model", "adjusted"): 0.020,
    ("fwh_entropy_svd_entropy", "woman_cluster_robust_OLS", "adjusted"): 0.256,
}
for key, want in expect.items():
    chk("table:amplitude " + "/".join(key), round(vals[key], 3) == want, f"saved {vals[key]:.4f} vs text {want}")

# S8C residual-RMS Spearman
A("\\subsection*{S8C. Residual association with RMS after adjustment (mixed model)}")
A("Spearman correlation between the model residual and natural-log RMS from the same row. Marginal residuals subtract only the fixed-effect prediction; conditional residuals also remove the woman and window random-effect predictions.")
A("\\begin{longtable}{p{3.0cm}p{2.2cm}rrrrr}")
A("\\toprule")
h = "\\textbf{Outcome} & \\textbf{Residual} & \\textbf{Pooled} & \\textbf{EHG9} & \\textbf{EHG10} & \\textbf{EHG11} & \\textbf{EHG12} \\\\"
A(h + "\n\\midrule\n\\endfirsthead\n\\toprule\n" + h + "\n\\midrule\n\\endhead")
rho_tot = None
for feat in ("fwh_entropy_totalsampen", "fwh_entropy_sampen_fixed_r", "fwh_entropy_svd_entropy"):
    for rt in ("conditional", "marginal"):
        sel = {r["scope"]: r for r in dia if r["feature"] == feat and r["model_method"] == "preferred_mixed_model"
               and r["model_variant"] == "adjusted" and r["residual_type"] == rt and r["amplitude_variable"] == "log_RMS"}
        if not sel:
            continue
        nm = [r["feature_label"] for r in dia if r["feature"] == feat][0]
        A(f"{tex(nm)} & {rt} & " + " & ".join(f(sel[s]["spearman_rho"], 2) for s in ("pooled", "EHG9", "EHG10", "EHG11", "EHG12")) + " \\\\")
        if feat == "fwh_entropy_totalsampen" and rt == "conditional":
            rho_tot = float(sel["pooled"]["spearman_rho"])
A("\\bottomrule\n\\end{longtable}")
chk("residual-RMS Spearman, Total SampEn, conditional pooled = -0.70", rho_tot is not None and round(rho_tot, 2) == -0.70, f"{rho_tot:.3f}")

# S8D woman-level change correlations
A("\\subsection*{S8D. Woman-level UC-minus-IN changes: Total SampEn versus amplitude (34 paired women)}")
A("Spearman correlation between each woman's median UC-minus-IN change in Total SampEn and the corresponding change in the amplitude summary, by channel.")
A("\\begin{longtable}{p{3.2cm}rrrr}")
A("\\toprule")
h = "\\textbf{Amplitude summary} & \\textbf{EHG9} & \\textbf{EHG10} & \\textbf{EHG11} & \\textbf{EHG12} \\\\"
A(h + "\n\\midrule\n\\endfirsthead\n\\toprule\n" + h + "\n\\midrule\n\\endhead")
rms_rho = []
for am in ("RMS", "MAV", "Peak-to-peak", "Waveform length"):
    cells = []
    for ch in ("EHG9", "EHG10", "EHG11", "EHG12"):
        rr = [r for r in wl if r["sampen_label"] == "Total SampEn" and r["scope"] == "patient_contraction_minus_baseline"
              and r["amplitude_label"] == am and r["channel"] == ch][0]
        cells.append(f(rr["spearman_rho"]))
        if am == "RMS":
            rms_rho.append(float(rr["spearman_rho"]))
    A(f"{am} & " + " & ".join(cells) + " \\\\")
A("\\bottomrule\n\\end{longtable}")
chk("woman-level change Spearman (RMS) range", True, f"{min(rms_rho):.3f}-{max(rms_rho):.3f} (text 0.92-0.96)")

# S8E GI refit
A("\\subsection*{S8E. EHG9 UC-minus-IN coefficients for original and GI Total SampEn under log-RMS conditioning}")
A("Three-state random-intercept model (woman and recording-within-woman), REML, 4,083 windows; custom block-diagonal REML implementation, checked against the cached unadjusted coefficients (agreement to four decimals). Log-RMS is standardised. Cross-checked against statsmodels 0.14.5 MixedLM (REML, woman intercept plus recording variance component): maximum absolute coefficient difference $3\\times10^{-6}$, maximum standard-error difference 0.003 (saved as statsmodels\\_verification.csv in the refit folder).")
A("\\begin{longtable}{p{2.4cm}p{3.2cm}rrrr}")
A("\\toprule")
h = "\\textbf{Outcome} & \\textbf{Adjustment} & \\textbf{UC $\\beta$} & \\textbf{95\\% CI} & \\textbf{FM $\\beta$} & \\textbf{log-RMS $\\beta$} \\\\"
A(h + "\n\\midrule\n\\endfirsthead\n\\toprule\n" + h + "\n\\midrule\n\\endhead")
groups = []
for r in gi:
    key = (r["outcome"], r["adjusted_log_RMS"])
    if key not in groups:
        groups.append(key)
print("GI refit groups:", groups)
for key in groups:
    rows = {r["term"]: r for r in gi if (r["outcome"], r["adjusted_log_RMS"]) == key}
    uc, fm = rows["UC_minus_IN"], rows["FM_minus_IN"]
    lr = [v for t_, v in rows.items() if t_.startswith("log_RMS")]
    A(f"{ {'frozen_total': 'Original Total', 'gi_total': 'GI Total'}.get(key[0], tex(key[0].replace('_', ' '))) } & {'standardised log-RMS' if key[1]=='True' else 'none'} & {f(uc['estimate'])} & {ci(uc['ci_low'], uc['ci_high'])} & {f(fm['estimate'])} & {f(lr[0]['estimate']) if lr else '--'} \\\\")
A("\\bottomrule\n\\end{longtable}")
gi_get = lambda o, a: float([r for r in gi if r["outcome"] == o and r["adjusted_log_RMS"] == a and r["term"] == "UC_minus_IN"][0]["estimate"])
chk("conditional table: original Total 0.409 -> 0.020", round(gi_get("frozen_total", "False"), 3) == 0.409 and round(gi_get("frozen_total", "True"), 3) == 0.020)
chk("conditional table: GI Total -0.113 -> +0.275", round(gi_get("gi_total", "False"), 3) == -0.113 and round(gi_get("gi_total", "True"), 3) == 0.275)

# S8F leave-one-woman-out
gu = [float(r["gi_total_unadj_UC"]) for r in lowo]
ga = [float(r["gi_total_adj_UC"]) for r in lowo]
fu = [float(r["frozen_total_unadj_UC"]) for r in lowo]
fa = [float(r["frozen_total_adj_UC"]) for r in lowo]
A("\\subsection*{S8F. Leave-one-woman-out refits (42 refits, each dropping one woman)}")
A("\\begin{longtable}{p{4.2cm}rrr}")
A("\\toprule")
h = "\\textbf{UC-minus-IN coefficient} & \\textbf{Minimum} & \\textbf{Median} & \\textbf{Maximum} \\\\"
A(h + "\n\\midrule\n\\endfirsthead\n\\toprule\n" + h + "\n\\midrule\n\\endhead")
for nm, arr in (("GI Total, unadjusted", gu), ("GI Total, log-RMS adjusted", ga), ("Original Total, unadjusted", fu), ("Original Total, log-RMS adjusted", fa)):
    A(f"{nm} & {f(min(arr))} & {f(st.median(arr))} & {f(max(arr))} \\\\")
A("\\bottomrule\n\\end{longtable}")
chk("LOWO: n refits", len(lowo) == 42, str(len(lowo)))
chk("LOWO: GI Total unadjusted stays negative / adjusted stays positive", max(gu) < 0 and min(ga) > 0, f"unadj {min(gu):.3f}..{max(gu):.3f}; adj {min(ga):.3f}..{max(ga):.3f}")

# S8G exclusion of the single extreme window
ex = rd(f"{COND}/refit_excl_ice037_window.csv")
def exg(st_, o, a, t, k="estimate"):
    return [r for r in ex if r["set"] == st_ and r["outcome"] == o and r["adjusted_log_RMS"] == a and r["term"] == t][0][k]
A("\\subsection*{S8G. Effect of the single extreme window on the conditional models}")
A("Same models as S8E, refitted after removing the one UC window from woman ice037 (samples 6,094--18,094, RMS 5.49 against 0.17 for the next largest window): 4,082 windows. Outcomes are re-standardised within each window set.")
A("\\begin{longtable}{p{2.4cm}p{3.2cm}rrrr}")
A("\\toprule")
h = "\\textbf{Outcome} & \\textbf{Adjustment} & \\textbf{Window set} & \\textbf{UC $\\beta$ (SE)} & \\textbf{95\\% CI} & \\textbf{FM $\\beta$} \\\\"
A(h + "\n\\midrule\n\\endfirsthead\n\\toprule\n" + h + "\n\\midrule\n\\endhead")
def sg(x):
    v = float(x); return ("$-$" if v < 0 else "") + f"{abs(v):.3f}"
for o, on in (("frozen_total", "Original Total"), ("gi_total", "GI Total")):
    for a, an in (("False", "none"), ("True", "standardised log-RMS")):
        for st_, sn in (("all", "all 4,083"), ("excl_ice037_window", "4,082")):
            A(f"{on} & {an} & {sn} & {sg(exg(st_,o,a,'UC_minus_IN'))} ({float(exg(st_,o,a,'UC_minus_IN','se')):.3f}) & [{sg(exg(st_,o,a,'UC_minus_IN','ci_low'))}, {sg(exg(st_,o,a,'UC_minus_IN','ci_high'))}] & {sg(exg(st_,o,a,'FM_minus_IN'))} \\\\")
A("\\bottomrule\n\\end{longtable}\n\\end{landscape}")
chk("S8G: all-window rows reproduce S8E", round(float(exg("all","frozen_total","False","UC_minus_IN")),3)==0.409 and round(float(exg("all","gi_total","True","UC_minus_IN")),3)==0.275)
chk("S8G: original Total excl. 0.692 / adj 0.203; GI -0.099 / 0.280", (round(float(exg("excl_ice037_window","frozen_total","False","UC_minus_IN")),3), round(float(exg("excl_ice037_window","frozen_total","True","UC_minus_IN")),3), round(float(exg("excl_ice037_window","gi_total","False","UC_minus_IN")),3), round(float(exg("excl_ice037_window","gi_total","True","UC_minus_IN")),3))==(0.692,0.203,-0.099,0.280))

with open(f"{OUT}/Supplement_S7_S8.tex", "w", encoding="utf-8") as fh:
    fh.write("\n".join(L) + "\n")
with open(f"{OUT}/Supplement_S7_S8_checks.txt", "w", encoding="utf-8") as fh:
    fh.write("\n".join(checks) + "\n")
print("\n".join(checks))
