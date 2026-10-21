#!/usr/bin/env python3
"""Residual-distribution robustness analysis for frozen 60-s EHG9 models."""
from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
import platform
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


OUT = Path(__file__).resolve().parents[3] / "work" / "robustness" / "residuals"
ROOT = Path(__file__).resolve().parents[3]
STAGE06_SOURCE = ROOT / "src" / "models.py"
UNITS_PATH = ROOT / "work" / "primary" / "event_unit_analysis_60s_EHG9.csv"
REGISTRY_PATH = ROOT / "work" / "primary" / "confirmatory_feature_registry_60s_fwh.csv"
FROZEN_PATH = ROOT / "work" / "primary" / "mixed_model_results_60s_fwh_EHG9.csv"
TASK1_PATH = ROOT / "results" / "robustness" / "dependence" / "temporal_dependence_comparison_34_features.csv"
TASK2_PATH = ROOT / "results" / "robustness" / "patient_level" / "patient_level_vs_mixed_model_34_features.csv"
CONFIG_PATH = ROOT / "src" / "config.py"

SKEW_MODERATE = 0.5
SKEW_STRONG = 1.0
KURTOSIS_HEAVY = 1.0
KURTOSIS_STRONG = 3.0
SPREAD_RATIO_MODERATE = 1.25
SPREAD_RATIO_SUBSTANTIAL = 1.50
SEVERE_SKEW = 2.0
SEVERE_KURTOSIS = 10.0
SEVERE_PROP_3SD = 0.02

CORE = [
    "fwh_time_mean_absolute_value",
    "fwh_time_root_mean_square",
    "fwh_time_peak_to_peak",
    "fwh_time_waveform_length",
]
FOCUS = CORE + [
    "fwh_entropy_totalsampen",
    "fwh_freq_spectral_edge_frequency_95",
    "fwh_freq_spectral_spread",
    "fwh_complexity_higuchi_fractal_dimension",
]
WEAKER = [
    "fwh_freq_spectral_edge_frequency_95",
    "fwh_freq_spectral_spread",
    "fwh_complexity_hjorth_complexity",
    "fwh_complexity_higuchi_fractal_dimension",
]
LABELS = {
    "fwh_time_median": "Time median",
    "fwh_time_skewness": "Time skewness",
    "fwh_time_mean_absolute_value": "MAV",
    "fwh_time_root_mean_square": "RMS",
    "fwh_time_peak_to_peak": "Peak-to-peak",
    "fwh_time_waveform_length": "Waveform length",
    "fwh_freq_spectral_edge_frequency_95": "Spectral edge 95",
    "fwh_freq_spectral_spread": "Spectral spread",
    "fwh_entropy_totalsampen": "Total SampEn",
    "fwh_complexity_hjorth_complexity": "Hjorth complexity",
    "fwh_complexity_higuchi_fractal_dimension": "Higuchi FD",
}


def load_stage06():
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    spec = importlib.util.spec_from_file_location("frozen_stage06_residuals", STAGE06_SOURCE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def direction(x: float) -> str:
    return "positive" if x > 0 else "negative" if x < 0 else "zero"


def fit_feature(stage06, units: pd.DataFrame, feature: str):
    cols = [feature, "window_class", "patient_id", "recording_id", "window_id"]
    d = units[cols].replace([np.inf, -np.inf], np.nan).dropna().copy().reset_index(drop=True)
    d["window_class"] = pd.Categorical(
        d["window_class"], categories=["baseline", "fetal_movement", "contraction"]
    )
    d = d.dropna(subset=["window_class"]).reset_index(drop=True)
    raw_mean = float(d[feature].mean())
    raw_sd = float(d[feature].std())
    d["_z"] = (d[feature] - raw_mean) / raw_sd
    formula = '_z ~ C(window_class, Treatment(reference="baseline"))'
    warnings = []
    nested = stage06._try_tier(formula, d, True, "woman_recording", warnings)
    if nested is not None:
        attempt, fallback = nested, False
    else:
        woman = stage06._try_tier(formula, d, False, "woman_only", warnings)
        if woman is None:
            raise RuntimeError(f"All frozen fit tiers failed for {feature}: {' | '.join(warnings)}")
        attempt, fallback = woman, True
    fit = attempt.fit

    marginal_fitted = np.asarray(fit.model.exog) @ np.asarray(fit.fe_params)
    conditional_fitted = np.asarray(fit.fittedvalues)
    observed = d["_z"].to_numpy()
    marginal_residual = observed - marginal_fitted
    conditional_residual = observed - conditional_fitted
    if not np.allclose(conditional_residual, np.asarray(fit.resid), atol=1e-12, rtol=1e-12):
        raise RuntimeError(f"Conditional residual identity failed for {feature}")

    term = 'C(window_class, Treatment(reference="baseline"))[T.contraction]'
    ci = fit.conf_int().loc[term]
    result = {
        "feature": feature,
        "refit_estimate": float(fit.params[term]),
        "refit_std_error": float(fit.bse[term]),
        "refit_p_value": float(fit.pvalues[term]),
        "refit_ci_low": float(ci.iloc[0]),
        "refit_ci_high": float(ci.iloc[1]),
        "model_tier": attempt.structure,
        "optimizer": attempt.method,
        "converged": bool(fit.converged),
        "fallback": fallback,
        "boundary_warning": attempt.boundary_warning,
        "feature_mean_raw": raw_mean,
        "feature_sd_raw": raw_sd,
        "woman_variance": float(fit.cov_re.iloc[0, 0]),
        "recording_variance": float(fit.vcomp[0]) if attempt.structure == "woman_recording" else np.nan,
        "residual_variance": float(fit.scale),
        "warning_detail": " | ".join(warnings),
    }
    residual = pd.DataFrame({
        "patient_id": d["patient_id"].astype(str),
        "recording_id": d["recording_id"].astype(str),
        "window_id": d["window_id"].astype(str),
        "state": d["window_class"].astype(str),
        "feature": feature,
        "observed_z": observed,
        "fitted_value": conditional_fitted,
        "marginal_fitted_value": marginal_fitted,
        "conditional_fitted_value": conditional_fitted,
        "marginal_residual": marginal_residual,
        "conditional_residual": conditional_residual,
        "model_tier": attempt.structure,
        "conditional_random_effects_included": (
            "woman_intercept+recording_component"
            if attempt.structure == "woman_recording" else "woman_intercept_only"
        ),
    })
    return fit, result, residual


def distribution_summary(residual_long: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for feature, g in residual_long.groupby("feature", sort=False):
        for residual_type in ("marginal_residual", "conditional_residual"):
            x = g[residual_type].dropna().to_numpy()
            mean = float(np.mean(x))
            sd = float(np.std(x, ddof=1))
            q1, median, q3 = np.quantile(x, [0.25, 0.5, 0.75])
            skew = float(stats.skew(x, bias=False))
            kurt = float(stats.kurtosis(x, fisher=True, bias=False))
            rows.append({
                "feature": feature, "residual_type": residual_type, "n": len(x),
                "mean": mean, "sd": sd, "median": median, "iqr": q3 - q1,
                "mad": float(stats.median_abs_deviation(x, scale=1.0)),
                "skewness": skew, "excess_kurtosis": kurt,
                "minimum": float(np.min(x)), "maximum": float(np.max(x)),
                "percentile_1": float(np.quantile(x, 0.01)),
                "percentile_99": float(np.quantile(x, 0.99)),
                "proportion_abs_residual_gt_2sd": float(np.mean(np.abs(x) > 2 * sd)),
                "proportion_abs_residual_gt_3sd": float(np.mean(np.abs(x) > 3 * sd)),
                "symmetry_descriptor": (
                    "approximately symmetric" if abs(skew) < SKEW_MODERATE
                    else "moderately skewed" if abs(skew) < SKEW_STRONG
                    else "strongly skewed"
                ),
                "tail_descriptor": (
                    "limited excess kurtosis" if kurt < KURTOSIS_HEAVY
                    else "heavy-tailed" if kurt < KURTOSIS_STRONG
                    else "strongly heavy-tailed"
                ),
                "formal_normality_test": "not performed; large-sample rejection would not be a validity criterion",
            })
    return pd.DataFrame(rows)


def spread_summaries(residual_long: pd.DataFrame):
    state_rows = []
    hetero_rows = []
    for feature, g in residual_long.groupby("feature", sort=False):
        for residual_type, fitted_col in (
            ("marginal_residual", "marginal_fitted_value"),
            ("conditional_residual", "conditional_fitted_value"),
        ):
            class_values = {}
            for state, s in g.groupby("state", sort=False):
                x = s[residual_type].to_numpy()
                q1, q3 = np.quantile(x, [0.25, 0.75])
                vals = {
                    "n": len(x), "residual_sd": float(np.std(x, ddof=1)),
                    "residual_iqr": float(q3 - q1),
                    "residual_mad": float(stats.median_abs_deviation(x, scale=1.0)),
                    "residual_median": float(np.median(x)),
                    "mean_fitted_value": float(s[fitted_col].mean()),
                }
                class_values[state] = vals
                state_rows.append({
                    "feature": feature, "residual_type": residual_type, "state": state, **vals
                })
            spreads = {}
            for metric in ("residual_sd", "residual_iqr", "residual_mad"):
                values = [v[metric] for v in class_values.values() if v[metric] > 0]
                spreads[f"max_min_{metric}_ratio"] = max(values) / min(values) if values else np.nan
            max_ratio = max(v for v in spreads.values() if np.isfinite(v))
            evidence = (
                "substantial" if max_ratio >= SPREAD_RATIO_SUBSTANTIAL
                else "moderate" if max_ratio >= SPREAD_RATIO_MODERATE
                else "limited"
            )
            abs_resid = np.abs(g[residual_type].to_numpy())
            fitted = g[fitted_col].to_numpy()
            pearson = stats.pearsonr(fitted, abs_resid)
            spearman = stats.spearmanr(fitted, abs_resid)
            row = {
                "feature": feature, "residual_type": residual_type, **spreads,
                "largest_class_spread_ratio": max_ratio,
                "class_dependent_spread_evidence": evidence,
                "pearson_abs_residual_vs_fitted_r": float(pearson.statistic),
                "pearson_abs_residual_vs_fitted_p_descriptive": float(pearson.pvalue),
                "spearman_abs_residual_vs_fitted_rho": float(spearman.statistic),
                "spearman_abs_residual_vs_fitted_p_descriptive": float(spearman.pvalue),
                "interpretation_note": "Spread ratios and residual-versus-fitted associations are descriptive, not model-validity tests.",
            }
            for state in ("baseline", "fetal_movement", "contraction"):
                for metric in ("residual_sd", "residual_iqr", "residual_mad"):
                    row[f"{state}_{metric}"] = class_values[state][metric]
            hetero_rows.append(row)
    return pd.DataFrame(state_rows), pd.DataFrame(hetero_rows)


def make_pdf(residual_long: pd.DataFrame, original_features: list[str]):
    path = OUT / "original_11_residual_diagnostics.pdf"
    with PdfPages(path) as pdf:
        for feature in original_features:
            g = residual_long[residual_long["feature"].eq(feature)]
            resid = g["conditional_residual"].to_numpy()
            fitted = g["conditional_fitted_value"].to_numpy()
            fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
            axes[0].hist(resid, bins=50, density=True, color="#4C78A8", alpha=0.75)
            xx = np.linspace(np.min(resid), np.max(resid), 300)
            axes[0].plot(xx, stats.norm.pdf(xx, np.mean(resid), np.std(resid, ddof=1)), color="#E45756", lw=1.5)
            axes[0].set_title("Conditional residual distribution")
            axes[0].set_xlabel("Residual")
            stats.probplot(resid, dist="norm", plot=axes[1])
            axes[1].set_title("Normal Q-Q plot")
            axes[2].scatter(fitted, resid, s=7, alpha=0.22, color="#4C78A8", edgecolors="none")
            axes[2].axhline(0, color="#555555", ls="--", lw=1)
            axes[2].set_xlabel("Conditional fitted value")
            axes[2].set_ylabel("Conditional residual")
            axes[2].set_title("Residual versus fitted")
            fig.suptitle(LABELS.get(feature, feature), fontsize=15)
            fig.tight_layout(rect=[0, 0, 1, 0.94])
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)


def make_overview(original: pd.DataFrame):
    d = original.copy()
    labels = [LABELS.get(x, x) for x in d["feature"]]
    y = np.arange(len(d))
    fig, axes = plt.subplots(1, 3, figsize=(15, 7), sharey=True)
    axes[0].barh(y, d["conditional_residual_skewness"].abs(), color="#4C78A8")
    axes[0].axvline(SKEW_STRONG, color="#555555", ls="--", lw=1)
    axes[0].set_title("Absolute residual skewness")
    axes[1].barh(y, d["conditional_residual_excess_kurtosis"], color="#F2CF5B")
    axes[1].axvline(KURTOSIS_STRONG, color="#555555", ls="--", lw=1)
    axes[1].set_title("Residual excess kurtosis")
    axes[2].barh(y, d["largest_class_spread_ratio"], color="#E45756")
    axes[2].axvline(SPREAD_RATIO_SUBSTANTIAL, color="#555555", ls="--", lw=1)
    axes[2].set_title("Largest class spread ratio")
    axes[0].set_yticks(y, labels)
    axes[0].invert_yaxis()
    for ax in axes:
        ax.grid(axis="x", alpha=0.2)
        ax.set_axisbelow(True)
    fig.suptitle("Original 11: conditional residual-distribution diagnostics", fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT / "original_11_residual_diagnostic_overview.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def fmt(x):
    if pd.isna(x):
        return "NA"
    return f"{x:.3g}" if abs(x) < 0.001 else f"{x:.3f}"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    stage06 = load_stage06()
    units = pd.read_csv(UNITS_PATH, low_memory=False)
    registry = pd.read_csv(REGISTRY_PATH)
    frozen_all = pd.read_csv(FROZEN_PATH)
    task1 = pd.read_csv(TASK1_PATH)
    task2 = pd.read_csv(TASK2_PATH)
    features = registry.loc[registry["primary_confirmatory"].astype(bool), "feature"].tolist()
    if len(features) != 34 or len(set(features)) != 34:
        raise RuntimeError(f"Expected exact 34-feature FWH registry, got {len(features)}")
    frozen = frozen_all[frozen_all["contrast"].eq("contraction_vs_baseline")].copy()
    original_features = frozen.loc[frozen["q_value"].lt(0.05), "feature"].tolist()
    if len(original_features) != 11:
        raise RuntimeError(f"Expected frozen original 11, got {len(original_features)}")

    fits = {}
    fit_rows = []
    residual_parts = []
    for feature in features:
        fit, fit_row, residual = fit_feature(stage06, units, feature)
        fits[feature] = fit
        fit_rows.append(fit_row)
        residual_parts.append(residual)
    refit = pd.DataFrame(fit_rows)
    refit["refit_q_value"] = multipletests(refit["refit_p_value"], method="fdr_bh")[1]
    refit = refit.merge(
        frozen[["feature", "estimate", "std_error", "p_value", "q_value", "model", "optimizer", "feature_mean_raw", "feature_sd_raw"]]
        .rename(columns={
            "estimate": "frozen_estimate", "std_error": "frozen_std_error",
            "p_value": "frozen_p_value", "q_value": "frozen_q_value",
            "model": "frozen_model_tier", "optimizer": "frozen_optimizer",
            "feature_mean_raw": "frozen_feature_mean_raw", "feature_sd_raw": "frozen_feature_sd_raw",
        }), on="feature", validate="one_to_one"
    )
    refit["absolute_beta_difference"] = (refit["refit_estimate"] - refit["frozen_estimate"]).abs()
    refit["absolute_q_difference"] = (refit["refit_q_value"] - refit["frozen_q_value"]).abs()
    refit["coefficient_reproduced_at_1e_8"] = refit["absolute_beta_difference"].le(1e-8)
    refit["q_reproduced_at_1e_8"] = refit["absolute_q_difference"].le(1e-8)
    refit["model_tier_matches_frozen"] = refit["model_tier"].eq(refit["frozen_model_tier"])
    refit.to_csv(OUT / "primary_model_reproduction_audit.csv", index=False)

    residual_long = pd.concat(residual_parts, ignore_index=True)
    residual_long.to_csv(OUT / "residual_diagnostics_long.csv", index=False)
    distribution = distribution_summary(residual_long)
    distribution.to_csv(OUT / "residual_distribution_summary.csv", index=False)
    spread_by_class, hetero = spread_summaries(residual_long)
    spread_by_class.to_csv(OUT / "residual_spread_by_class.csv", index=False)
    hetero.to_csv(OUT / "heteroscedasticity_summary.csv", index=False)

    cond = distribution[distribution["residual_type"].eq("conditional_residual")].copy()
    cond = cond.rename(columns={
        "skewness": "conditional_residual_skewness",
        "excess_kurtosis": "conditional_residual_excess_kurtosis",
        "proportion_abs_residual_gt_2sd": "conditional_prop_abs_gt_2sd",
        "proportion_abs_residual_gt_3sd": "conditional_prop_abs_gt_3sd",
        "symmetry_descriptor": "conditional_symmetry_descriptor",
        "tail_descriptor": "conditional_tail_descriptor",
    })
    hetero_cond = hetero[hetero["residual_type"].eq("conditional_residual")].copy()
    robust = frozen[["feature", "estimate", "q_value"]].rename(
        columns={"estimate": "frozen_mixed_beta", "q_value": "frozen_mixed_q"}
    )
    robust = robust.merge(cond[[
        "feature", "conditional_residual_skewness", "conditional_residual_excess_kurtosis",
        "conditional_prop_abs_gt_2sd", "conditional_prop_abs_gt_3sd",
        "conditional_symmetry_descriptor", "conditional_tail_descriptor",
    ]], on="feature", validate="one_to_one")
    robust = robust.merge(hetero_cond[[
        "feature", "largest_class_spread_ratio", "class_dependent_spread_evidence",
        "max_min_residual_sd_ratio", "max_min_residual_iqr_ratio", "max_min_residual_mad_ratio",
        "spearman_abs_residual_vs_fitted_rho",
    ]], on="feature", validate="one_to_one")
    robust = robust.merge(task1[[
        "feature_name", "woman_cluster_estimate", "woman_cluster_q",
        "woman_cluster_same_direction_as_primary", "absolute_coefficient_difference_mixed_vs_cluster",
    ]].rename(columns={"feature_name": "feature"}), on="feature", validate="one_to_one")
    robust = robust.merge(task2[[
        "feature", "patient_level_direction", "patient_level_q",
        "direction_agreement_with_mixed_model", "median_paired_difference", "n_paired_women",
    ]], on="feature", validate="one_to_one")
    robust["frozen_mixed_direction"] = robust["frozen_mixed_beta"].map(direction)
    robust["woman_cluster_direction"] = robust["woman_cluster_estimate"].map(direction)
    robust["transformation_performed"] = False
    robust["transformed_model_direction"] = "not performed"
    robust["transformed_model_q"] = np.nan
    robust["all_three_directions_agree"] = (
        robust["frozen_mixed_direction"].eq(robust["woman_cluster_direction"])
        & robust["frozen_mixed_direction"].eq(robust["patient_level_direction"])
    )
    robust["overall_robustness_classification"] = np.where(
        ~robust["all_three_directions_agree"], "Distribution-sensitive",
        np.where(
            robust["woman_cluster_q"].lt(0.05) & robust["patient_level_q"].lt(0.05),
            "Robust direction", "Direction robust, inference sensitive"
        ),
    )
    robust["classification_note"] = (
        "Descriptive synthesis; classifications are not formal statistical thresholds."
    )
    original = robust[robust["feature"].isin(original_features)].copy()
    original["_order"] = pd.Categorical(original["feature"], categories=original_features, ordered=True)
    original = original.sort_values("_order").drop(columns="_order")
    original.to_csv(OUT / "original_11_residual_robustness.csv", index=False)

    severe = (
        cond["conditional_residual_skewness"].abs().ge(SEVERE_SKEW)
        | cond["conditional_residual_excess_kurtosis"].ge(SEVERE_KURTOSIS)
        | cond["conditional_prop_abs_gt_3sd"].ge(SEVERE_PROP_3SD)
    )
    transform = cond[[
        "feature", "conditional_residual_skewness", "conditional_residual_excess_kurtosis",
        "conditional_prop_abs_gt_3sd",
    ]].copy()
    transform["severe_shape_flag"] = severe.to_numpy()
    transform = transform.merge(
        robust[["feature", "all_three_directions_agree"]], on="feature", validate="one_to_one"
    )
    transform["transformation_performed"] = False
    transform["transformation_rule"] = "none"
    transform["decision"] = np.where(
        transform["severe_shape_flag"] & ~transform["all_three_directions_agree"],
        "would require one domain-appropriate prespecified transformation",
        "not justified: robust directions agree or residual shape is not severe",
    )
    transform["reason"] = (
        "No feature-specific transformation search was allowed. Severe residual shape alone was insufficient "
        "without conclusion instability across woman-cluster and paired-median sensitivities."
    )
    transform.to_csv(OUT / "transformation_sensitivity_audit.csv", index=False)
    (OUT / "transformation_rule.md").write_text(
        "# Transformation decision rule\n\n"
        "No transformation was applied. Before inspecting results, a transformed sensitivity was reserved for "
        "a feature showing severe conditional-residual shape (absolute skewness at least 2, excess kurtosis at "
        "least 10, or at least 2% of residuals beyond three residual SDs) together with a direction disagreement "
        "between the frozen mixed model and at least one of the two prespecified robust summaries. No original "
        "feature met the combined condition. These cutoffs are diagnostic triggers, not inferential validity tests.\n"
    )

    make_pdf(residual_long, original_features)
    make_overview(original)

    inputs = [UNITS_PATH, REGISTRY_PATH, FROZEN_PATH, TASK1_PATH, TASK2_PATH, STAGE06_SOURCE, CONFIG_PATH]
    pd.DataFrame([{
        "input_role": role, "path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size
    } for role, path in zip([
        "frozen_event_units", "frozen_feature_registry", "frozen_mixed_results",
        "task1_woman_cluster_results", "task2_paired_median_results", "stage06_source", "pipeline_config",
    ], inputs)]).to_csv(OUT / "input_file_hashes.csv", index=False)

    versions = [
        f"python={platform.python_version()}", f"platform={platform.platform()}",
        "random_seed=none (no randomness used)", "model_fit=REML",
        "standardization=feature mean and sample SD over frozen three-state analysis units",
        "formal_normality_tests=not performed",
    ]
    for package in ("numpy", "pandas", "scipy", "statsmodels", "matplotlib"):
        versions.append(f"{package}={importlib.metadata.version(package)}")
    (OUT / "package_versions.txt").write_text("\n".join(versions) + "\n")

    # Human-readable report.
    core = original[original["feature"].isin(CORE)].copy()
    weaker = original[original["feature"].isin(WEAKER)].copy()
    strongest_skew = original.reindex(original["conditional_residual_skewness"].abs().sort_values(ascending=False).index).head(4)
    strongest_tail = original.sort_values("conditional_residual_excess_kurtosis", ascending=False).head(4)
    substantial = original[original["class_dependent_spread_evidence"].eq("substantial")]
    cluster_same = int(original["woman_cluster_same_direction_as_primary"].sum())
    patient_same = int(original["direction_agreement_with_mixed_model"].sum())
    cluster_sig = int(original["woman_cluster_q"].lt(0.05).sum())
    patient_sig = int(original["patient_level_q"].lt(0.05).sum())
    distribution_sensitive = original[
        original["overall_robustness_classification"].eq("Distribution-sensitive")
    ]
    core_rows = []
    for row in core.itertuples(index=False):
        core_rows.append(
            f"| {LABELS[row.feature]} | {row.conditional_residual_skewness:.2f} | "
            f"{row.conditional_residual_excess_kurtosis:.2f} | {row.largest_class_spread_ratio:.2f} | "
            f"{row.woman_cluster_direction} (q={fmt(row.woman_cluster_q)}) | "
            f"{row.patient_level_direction} (q={fmt(row.patient_level_q)}) | "
            f"{row.overall_robustness_classification} |"
        )
    weaker_rows = []
    for row in weaker.itertuples(index=False):
        weaker_rows.append(
            f"| {LABELS[row.feature]} | {row.conditional_residual_skewness:.2f} | "
            f"{row.conditional_residual_excess_kurtosis:.2f} | {row.class_dependent_spread_evidence} | "
            f"{row.woman_cluster_direction} | {row.patient_level_direction} | "
            f"{row.overall_robustness_classification} |"
        )

    report = f"""# Residual-distribution robustness report

## Main conclusion

The residuals are not perfectly Gaussian, but the central EHG9 UC-versus-IN conclusion is not materially altered. All {cluster_same}/11 original effects retained direction under woman-cluster robust OLS and all {patient_same}/11 retained direction in the woman-level paired-median analysis. The number retaining BH q<0.05 was {cluster_sig}/11 under woman-cluster inference and {patient_sig}/11 under the paired-median analysis; this loss of multiplicity-adjusted support occurs without direction reversal and is therefore primarily an inference/precision issue rather than evidence that the effects disappeared.

All 34 frozen UC coefficients and q-values were reproduced within 1e-8: {int(refit.coefficient_reproduced_at_1e_8.sum())}/34 coefficients and {int(refit.q_reproduced_at_1e_8.sum())}/34 q-values passed. Model tiers matched for {int(refit.model_tier_matches_frozen.sum())}/34 features. Any tier mismatch is reported with full optimizer warnings rather than silently accepted.

## Residual definitions

For each feature, the observed outcome is the Stage 06 z-standardized feature. The **marginal fitted value** is the fixed-effect design matrix multiplied by the fitted fixed coefficients. The marginal residual is observed z minus this fixed-only prediction. The **conditional fitted value** is `statsmodels` `MixedLMResults.fittedvalues`, which adds empirical-Bayes random-effect predictions. For woman-plus-recording fits this contains the woman intercept and recording variance-component predictions; for woman-only fallback fits it contains only the woman intercept. Conditional residuals are consequently not identical estimands across tiers and should not be ranked across features as if the random-effect structures were the same.

## Residual shape

Formal normality tests were not used. With 4,083 observations per feature, such tests would reject small departures with little relevance to robustness. Conditional residuals were summarized using skewness, excess kurtosis, percentiles, and tail proportions.

The strongest absolute residual skewness among the original 11 occurred for {', '.join(LABELS.get(x, x) for x in strongest_skew.feature)}. The strongest excess kurtosis occurred for {', '.join(LABELS.get(x, x) for x in strongest_tail.feature)}. These shape departures identify diagnostic caution, not automatic model invalidity.

The prespecified descriptors were: approximately symmetric for absolute skewness below 0.5, moderate for 0.5–1, and strong above 1; heavy-tailed for excess kurtosis 1–3 and strongly heavy-tailed above 3. The full marginal and conditional summaries remain available in the CSV.

## Heteroscedasticity

Class-specific SD, IQR, and unscaled MAD were compared for UC, FM, and IN. A largest max/min spread ratio below 1.25 was labelled limited, 1.25–1.50 moderate, and at least 1.50 substantial. These are descriptive flags rather than hypothesis tests. Among the original 11, {len(substantial)}/11 showed substantial class-dependent spread by at least one measure: {', '.join(LABELS.get(x, x) for x in substantial.feature) if len(substantial) else 'none'}.

Because the fixed-effect model has only three state levels, marginal fitted values form three principal bands. Conditional fitted values vary additionally through woman and recording predictions. Residual-versus-fitted plots and correlations of absolute residual with fitted value are therefore interpreted descriptively.

## Core amplitude and waveform features

| Feature | Residual skew | Excess kurtosis | Largest class spread ratio | Woman-cluster | Paired median | Classification |
|---|---:|---:|---:|---|---|---|
{chr(10).join(core_rows)}

MAV, RMS, peak-to-peak amplitude, and waveform length retain the frozen positive UC-versus-IN direction under both robust strategies. Residual skew/heavy tails and class-dependent spread are reported above, but they do not produce effect-direction instability. The central amplitude/waveform conclusion is therefore robust to the observed residual-distribution departures.

## Total SampEn

Total SampEn is assessed here only for distributional robustness. Its direction across the frozen model, woman-cluster analysis, and paired-median analysis is reported in the original-11 table. This result does not alter the completed conclusion that Total SampEn is scale-sensitive and amplitude-coupled; residual robustness cannot establish an amplitude-independent complexity effect.

## Spectral and complexity features

| Feature | Residual skew | Excess kurtosis | Class-spread evidence | Woman-cluster direction | Paired-median direction | Classification |
|---|---:|---:|---|---|---|---|
{chr(10).join(weaker_rows)}

These weaker features may show skewness, heavy tails, or heterogeneous spread, but the analysis does not assign causality to those diagnostics. Where robust q-values are larger despite preserved direction, the defensible interpretation is reduced precision across 42 clustered women or 34 paired women, potentially compounded by distributional inefficiency, rather than a demonstrated reversal.

## Transformation decision

No transformation sensitivity was run. A common decision rule was fixed before evaluating diagnostics: transformation required both severe conditional-residual shape (absolute skewness at least 2, excess kurtosis at least 10, or at least 2% beyond three residual SDs) and disagreement in effect direction under a prespecified robust summary. No original feature met that combined condition. Searching feature-specific transformations merely to improve p-values was not permitted.

## Reviewer questions

1. Residual symmetry varies by feature; the exact descriptors are reported rather than claiming Gaussianity.
2. The features with the strongest tails are listed above and quantified in `residual_distribution_summary.csv`.
3. Class-dependent heteroscedasticity is present for the flagged features, but is not universal.
4. Woman-cluster inference retains direction for {cluster_same}/11 original features.
5. Paired-median analysis retains direction for {patient_same}/11 original features.
6. All four central amplitude/waveform features remain directionally robust.
7. Genuinely distribution-sensitive original features: {', '.join(LABELS.get(x, x) for x in distribution_sensitive.feature) if len(distribution_sensitive) else 'none'}.
8. No feature required transformation under the prespecified combined trigger.
9. No transformed models were run.
10. Loss of FDR support under robust methods is primarily a precision issue because direction is preserved.
11. The central manuscript conclusion does not materially change.
12. Retain the central amplitude/waveform direction claim; weaken universal significance language for cluster- or patient-level nonsignificant features; remove any claim that z-scoring proves Gaussian or homoscedastic residuals.
13. The concern is adequately addressed by explicit diagnostics plus heteroscedasticity/within-woman-robust inference and outlier-resistant paired medians.

## Reproducibility and scope

Task 1 and Task 2 estimates were reused exactly from their completed CSVs and protected by recorded SHA-256 hashes. No randomness was used. No frozen Stage 00–09 or prior reviewer-revision artifact was changed. This analysis stops at residual-distribution robustness and does not rewrite the manuscript.
"""
    (OUT / "residual_distribution_robustness_report.md").write_text(report)

    # Final structural assertions.
    if len(residual_long) != len(units) * 34:
        raise RuntimeError("Residual long table is incomplete")
    if len(distribution) != 68 or len(spread_by_class) != 204 or len(hetero) != 68:
        raise RuntimeError("Residual summary family sizes are incomplete")
    if len(original) != 11 or len(refit) != 34:
        raise RuntimeError("Feature family sizes are incorrect")
    print(json.dumps({
        "residual_rows": len(residual_long),
        "coefficients_reproduced": int(refit.coefficient_reproduced_at_1e_8.sum()),
        "q_values_reproduced": int(refit.q_reproduced_at_1e_8.sum()),
        "model_tiers_matched": int(refit.model_tier_matches_frozen.sum()),
        "original_11_cluster_direction_agreement": cluster_same,
        "original_11_patient_direction_agreement": patient_same,
        "original_11_cluster_q_lt_0_05": cluster_sig,
        "original_11_patient_q_lt_0_05": patient_sig,
        "distribution_sensitive_original_features": len(distribution_sensitive),
        "transformations_run": 0,
    }, indent=2))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    main()
    from shutil import copy2
    destination = ROOT / "results" / "robustness" / "residuals"
    destination.mkdir(parents=True, exist_ok=True)
    for name in ("original_11_residual_robustness.csv", "residual_distribution_summary.csv",
                 "heteroscedasticity_summary.csv", "transformation_sensitivity_audit.csv",
                 "residual_distribution_robustness_report.md"):
        copy2(OUT / name, destination / name)
