#!/usr/bin/env python3
"""RMS-adjusted entropy analysis with corrected physical-window dependence.

Primary corrected model (locked to cross-channel Phase 1/2):
    z_entropy ~ state + channel [+ z_log_RMS]
    ML MixedLM with woman random intercept and physical-window variance component.

Required sensitivity:
    identical fixed effects with woman-cluster-robust OLS covariance (CR1,
    finite-cluster t inference).

The script also refits the frozen older three-state woman+recording model in
memory to reproduce Table 8 and audit its residual/RMS calculation. It writes
only inside this Phase-2 directory.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import json
import platform
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm, pearsonr, spearmanr, t
import statsmodels
import statsmodels.formula.api as smf


SCRIPT_PATH = Path(__file__).resolve()
OUTPUT_DIR = Path(__file__).resolve().parents[2] / "results" / "entropy"
PIPELINE_ROOT = SCRIPT_PATH.parents[2]
STAGE06_DIR = PIPELINE_ROOT / "work" / "primary"
STAGE05_LONG = PIPELINE_ROOT / "work" / "feature_qc" / "window_features_60s_long.csv"
CHANNELS = ("EHG9", "EHG10", "EHG11", "EHG12")
FEATURES = (
    "fwh_entropy_totalsampen",
    "fwh_entropy_sampen_fixed_r",
    "fwh_entropy_svd_entropy",
)
LABELS = {
    "fwh_entropy_totalsampen": "Total SampEn",
    "fwh_entropy_sampen_fixed_r": "fixed-r SampEn",
    "fwh_entropy_svd_entropy": "SVD entropy",
}
RMS_FEATURE = "fwh_time_root_mean_square"
STATE_TERM = 'C(state, Treatment(reference="IN"))[T.UC]'
STATE_FACTOR = 'C(state, Treatment(reference="IN"))'
CHANNEL_FACTOR = 'C(channel, Treatment(reference="EHG9"))'
UNADJUSTED_FORMULA = f"z_entropy ~ {STATE_FACTOR} + {CHANNEL_FACTOR}"
ADJUSTED_FORMULA = UNADJUSTED_FORMULA + " + z_log_RMS"
OPTIMIZERS = ("lbfgs", "bfgs")
ALPHA = 0.05
REPRO_ATOL = 1e-10
REPRO_RTOL = 1e-8
RANDOM_SEED = None


@dataclass
class MixedAttempt:
    optimizer: str
    model_variant: str
    fit: object | None
    converged: bool
    finite_llf: bool
    fixed_covariance_finite: bool
    warnings: list[str]
    exception: str

    @property
    def usable(self) -> bool:
        return bool(self.fit is not None and self.converged and self.finite_llf and self.fixed_covariance_finite)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def load_corrected_population() -> tuple[pd.DataFrame, dict]:
    frames = []
    for channel in CHANNELS:
        frame = pd.read_csv(STAGE06_DIR / f"event_unit_analysis_60s_{channel}.csv", low_memory=False)
        frame["channel"] = channel
        frames.append(frame)
    pooled = pd.concat(frames, ignore_index=True)
    data = pooled.loc[pooled["window_class"].isin(["baseline", "contraction"])].copy()
    data["state"] = data["window_class"].map({"baseline": "IN", "contraction": "UC"})
    data["state"] = pd.Categorical(data["state"], categories=["IN", "UC"])
    data["channel"] = pd.Categorical(data["channel"], categories=list(CHANNELS))
    data["physical_window_key"] = data["window_id"].astype(str)
    per_window = data.groupby("physical_window_key", observed=True).agg(
        rows=("channel", "size"), channels=("channel", "nunique"), women=("patient_id", "nunique"),
        recordings=("recording_id", "nunique"), states=("state", "nunique"),
    )
    if not (
        per_window["rows"].eq(4).all() and per_window["channels"].eq(4).all()
        and per_window["women"].eq(1).all() and per_window["recordings"].eq(1).all()
        and per_window["states"].eq(1).all()
    ):
        raise RuntimeError("Corrected UC/IN physical-window structure is incomplete or inconsistent")
    unique_windows = data.drop_duplicates("physical_window_key")
    population = {
        "women": int(data["patient_id"].nunique()),
        "recordings": int(data["recording_id"].nunique()),
        "physical_windows": int(data["physical_window_key"].nunique()),
        "rows": int(len(data)),
        "uc_windows": int((unique_windows["state"] == "UC").sum()),
        "in_windows": int((unique_windows["state"] == "IN").sum()),
        "channels": int(data["channel"].nunique()),
        "rows_per_window": 4,
    }
    return data, population


def prepare_feature(population: pd.DataFrame, feature: str) -> tuple[pd.DataFrame, dict]:
    columns = [
        feature, RMS_FEATURE, "patient_id", "recording_id", "window_id", "physical_window_key",
        "state", "channel", "window_class",
    ]
    data = population[columns].replace([np.inf, -np.inf], np.nan).dropna().reset_index(drop=True).copy()
    if (data[RMS_FEATURE] <= 0).any():
        raise RuntimeError(f"{feature}: non-positive RMS prevents log transformation")
    entropy_mean = float(data[feature].mean())
    entropy_sd = float(data[feature].std(ddof=1))
    log_rms = np.log(data[RMS_FEATURE].to_numpy(dtype=float))
    log_rms_mean = float(np.mean(log_rms))
    log_rms_sd = float(np.std(log_rms, ddof=1))
    if entropy_sd <= 0 or log_rms_sd <= 0:
        raise RuntimeError(f"{feature}: invalid standardization SD")
    data["z_entropy"] = (data[feature] - entropy_mean) / entropy_sd
    data["log_RMS"] = log_rms
    data["z_log_RMS"] = (log_rms - log_rms_mean) / log_rms_sd
    stats = {
        "feature": feature, "feature_label": LABELS[feature], "n_rows": len(data),
        "n_women": data["patient_id"].nunique(), "n_recordings": data["recording_id"].nunique(),
        "n_physical_windows": data["physical_window_key"].nunique(),
        "n_uc_windows": data.loc[data["state"] == "UC", "physical_window_key"].nunique(),
        "n_in_windows": data.loc[data["state"] == "IN", "physical_window_key"].nunique(),
        "entropy_mean": entropy_mean, "entropy_sd_ddof1": entropy_sd,
        "raw_rms_mean": float(data[RMS_FEATURE].mean()), "raw_rms_sd_ddof1": float(data[RMS_FEATURE].std(ddof=1)),
        "log_rms_mean": log_rms_mean, "log_rms_sd_ddof1": log_rms_sd,
        "rms_transform": "natural logarithm of same-row, same-channel FWH RMS",
        "standardization": "pooled UC/IN rows across EHG9-EHG12; sample SD (ddof=1)",
    }
    return data, stats


def try_mixed(formula: str, data: pd.DataFrame, optimizer: str, variant: str) -> MixedAttempt:
    fit = None
    messages: list[str] = []
    exception = ""
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            model = smf.mixedlm(
                formula, data, groups=data["patient_id"], re_formula="1",
                vc_formula={"physical_window": "0 + C(physical_window_key)"},
            )
            fit = model.fit(reml=False, method=optimizer, maxiter=1000, disp=False)
        messages = [str(item.message) for item in caught]
    except Exception as exc:
        exception = f"{type(exc).__name__}: {exc}"
    converged = bool(getattr(fit, "converged", False)) if fit is not None else False
    finite_llf = bool(np.isfinite(float(getattr(fit, "llf", np.nan)))) if fit is not None else False
    covariance_finite = False
    if fit is not None:
        try:
            names = list(fit.fe_params.index)
            covariance_finite = bool(np.isfinite(fit.cov_params().loc[names, names].to_numpy(dtype=float)).all())
        except Exception as exc:
            exception = " | ".join(filter(None, [exception, f"covariance: {type(exc).__name__}: {exc}"]))
    return MixedAttempt(optimizer, variant, fit, converged, finite_llf, covariance_finite, messages, exception)


def fit_preferred(formula: str, data: pd.DataFrame, variant: str) -> tuple[MixedAttempt | None, list[MixedAttempt]]:
    attempts = []
    for optimizer in OPTIMIZERS:
        attempt = try_mixed(formula, data, optimizer, variant)
        attempts.append(attempt)
        if attempt.usable:
            return attempt, attempts
    return None, attempts


def fit_robust(formula: str, data: pd.DataFrame):
    ordinary = smf.ols(formula, data=data).fit()
    return ordinary.get_robustcov_results(
        cov_type="cluster", groups=data["patient_id"], use_correction=True,
        df_correction=True, use_t=True,
    )


def mixed_term(fit, term: str) -> dict:
    estimate = float(fit.fe_params[term])
    se = float(fit.bse_fe[term])
    z_value = estimate / se
    critical = float(norm.ppf(0.975))
    return {"estimate": estimate, "se": se, "ci_low": estimate - critical * se, "ci_high": estimate + critical * se, "statistic": z_value, "df": np.nan, "p_value": float(2 * norm.sf(abs(z_value)))}


def robust_term(fit, term: str) -> dict:
    names = list(fit.model.exog_names)
    index = names.index(term)
    estimate = float(fit.params[index])
    se = float(fit.bse[index])
    statistic = estimate / se
    df = float(fit.df_resid_inference)
    critical = float(t.ppf(0.975, df))
    return {"estimate": estimate, "se": se, "ci_low": estimate - critical * se, "ci_high": estimate + critical * se, "statistic": statistic, "df": df, "p_value": float(2 * t.sf(abs(statistic), df))}


def mixed_variances(fit) -> dict:
    vc_names = list(getattr(fit.model.exog_vc, "names", []))
    vc_values = np.asarray(getattr(fit, "vcomp", []), dtype=float)
    vc_map = dict(zip(vc_names, vc_values))
    return {
        "woman_variance": float(fit.cov_re.iloc[0, 0]),
        "physical_window_variance": float(vc_map.get("physical_window", np.nan)),
        "residual_variance": float(fit.scale),
    }


def corrected_models(population: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict, pd.DataFrame]:
    result_rows = []
    convergence_rows = []
    selected_adjusted: dict[str, tuple[pd.DataFrame, object, str]] = {}
    scaling_rows = []
    for feature in FEATURES:
        print(f"Corrected models: {LABELS[feature]}", flush=True)
        data, scaling = prepare_feature(population, feature)
        scaling_rows.append(scaling)
        fitted = {}
        for variant, formula in (("unadjusted", UNADJUSTED_FORMULA), ("adjusted", ADJUSTED_FORMULA)):
            preferred, attempts = fit_preferred(formula, data, variant)
            robust = fit_robust(formula, data)
            fitted[(variant, "preferred_mixed_model")] = preferred.fit if preferred else None
            fitted[(variant, "woman_cluster_robust_OLS")] = robust
            for attempt_number, attempt in enumerate(attempts, 1):
                fit = attempt.fit
                variances = mixed_variances(fit) if fit is not None else {"woman_variance": np.nan, "physical_window_variance": np.nan, "residual_variance": np.nan}
                boundary = any(any(token in message.lower() for token in ("boundary", "singular", "hessian", "positive definite")) for message in attempt.warnings)
                convergence_rows.append({
                    "analysis_phase": "corrected_UC_IN", "feature": feature, "feature_label": LABELS[feature],
                    "model_variant": variant, "method": "preferred_mixed_model",
                    "model_specification": "ML MixedLM; woman random intercept + physical-window variance component",
                    "attempt_number": attempt_number, "optimizer": attempt.optimizer,
                    "selected_attempt": bool(preferred is attempt), "usable": attempt.usable,
                    "converged": attempt.converged, "finite_llf": attempt.finite_llf,
                    "llf": float(getattr(fit, "llf", np.nan)) if fit is not None else np.nan,
                    "fixed_covariance_finite": attempt.fixed_covariance_finite,
                    "boundary_or_singularity_warning": boundary,
                    "warning_messages_json": json.dumps(attempt.warnings, ensure_ascii=False),
                    "exception": attempt.exception, **variances,
                    "n_rows": len(data), "n_women": data["patient_id"].nunique(),
                    "n_recordings": data["recording_id"].nunique(),
                    "n_physical_windows": data["physical_window_key"].nunique(),
                })
            convergence_rows.append({
                "analysis_phase": "corrected_UC_IN", "feature": feature, "feature_label": LABELS[feature],
                "model_variant": variant, "method": "woman_cluster_robust_OLS",
                "model_specification": "OLS with CR1 woman-cluster covariance, df correction, t reference",
                "attempt_number": 1, "optimizer": "closed_form_OLS", "selected_attempt": preferred is None,
                "usable": True, "converged": True, "finite_llf": bool(np.isfinite(robust.llf)),
                "llf": float(robust.llf), "fixed_covariance_finite": bool(np.isfinite(robust.cov_params()).all()),
                "boundary_or_singularity_warning": False, "warning_messages_json": "[]", "exception": "",
                "woman_variance": np.nan, "physical_window_variance": np.nan,
                "residual_variance": float(np.mean(np.square(robust.resid))),
                "n_rows": len(data), "n_women": data["patient_id"].nunique(),
                "n_recordings": data["recording_id"].nunique(),
                "n_physical_windows": data["physical_window_key"].nunique(),
            })
        primary_method = "preferred_mixed_model" if all(fitted[(v, "preferred_mixed_model")] is not None for v in ("unadjusted", "adjusted")) else "woman_cluster_robust_OLS"
        selected_adjusted[feature] = (data, fitted[("adjusted", primary_method)], primary_method)
        for method in ("preferred_mixed_model", "woman_cluster_robust_OLS"):
            if fitted[("unadjusted", method)] is None or fitted[("adjusted", method)] is None:
                continue
            term_func = mixed_term if method == "preferred_mixed_model" else robust_term
            unadjusted_state = term_func(fitted[("unadjusted", method)], STATE_TERM)
            adjusted_state = term_func(fitted[("adjusted", method)], STATE_TERM)
            amplitude = term_func(fitted[("adjusted", method)], "z_log_RMS")
            absolute_change = abs(adjusted_state["estimate"] - unadjusted_state["estimate"])
            for variant, state in (("unadjusted", unadjusted_state), ("adjusted", adjusted_state)):
                result_rows.append({
                    "feature": feature, "feature_label": LABELS[feature], "method": method,
                    "method_role": "primary" if method == primary_method else "required_sensitivity",
                    "model_variant": variant, "formula": UNADJUSTED_FORMULA if variant == "unadjusted" else ADJUSTED_FORMULA,
                    "uc_vs_in_estimate": state["estimate"], "uc_vs_in_se": state["se"],
                    "uc_vs_in_ci_low": state["ci_low"], "uc_vs_in_ci_high": state["ci_high"],
                    "uc_vs_in_statistic": state["statistic"], "reference_df": state["df"],
                    "uc_vs_in_p_value": state["p_value"],
                    "absolute_uc_coefficient_change_after_adjustment": absolute_change,
                    "z_log_rms_estimate": amplitude["estimate"] if variant == "adjusted" else np.nan,
                    "z_log_rms_se": amplitude["se"] if variant == "adjusted" else np.nan,
                    "z_log_rms_ci_low": amplitude["ci_low"] if variant == "adjusted" else np.nan,
                    "z_log_rms_ci_high": amplitude["ci_high"] if variant == "adjusted" else np.nan,
                    "z_log_rms_p_value": amplitude["p_value"] if variant == "adjusted" else np.nan,
                    "n_rows": len(data), "n_women": data["patient_id"].nunique(),
                    "n_recordings": data["recording_id"].nunique(),
                    "n_physical_windows": data["physical_window_key"].nunique(),
                    "outcome_standardization": "pooled UC/IN EHG9-EHG12; sample SD ddof=1",
                    "rms_transformation": "natural log then pooled UC/IN EHG9-EHG12 z-score; sample SD ddof=1",
                })
    return pd.DataFrame(result_rows), pd.DataFrame(convergence_rows), selected_adjusted, pd.DataFrame(scaling_rows)


def residual_vectors(data: pd.DataFrame, fit, method: str) -> dict[str, np.ndarray]:
    if method == "preferred_mixed_model":
        fixed_prediction = np.asarray(fit.model.exog @ fit.fe_params.to_numpy(dtype=float), dtype=float)
        marginal = data["z_entropy"].to_numpy(dtype=float) - fixed_prediction
        conditional = np.asarray(fit.resid, dtype=float)
        return {"marginal": marginal, "conditional": conditional}
    prediction = np.asarray(fit.model.exog @ np.asarray(fit.params), dtype=float)
    marginal = data["z_entropy"].to_numpy(dtype=float) - prediction
    return {"marginal": marginal}


def correlation_pair(y: np.ndarray, x: np.ndarray) -> dict:
    finite = np.isfinite(y) & np.isfinite(x)
    y, x = y[finite], x[finite]
    pearson = pearsonr(y, x)
    spearman = spearmanr(y, x)
    return {
        "n": len(y), "pearson_r": float(pearson.statistic), "pearson_p": float(pearson.pvalue),
        "spearman_rho": float(spearman.statistic), "spearman_p": float(spearman.pvalue),
    }


def corrected_residual_diagnostics(selected_adjusted: dict) -> pd.DataFrame:
    rows = []
    for feature, (data, fit, method) in selected_adjusted.items():
        vectors = residual_vectors(data, fit, method)
        amplitudes = {
            "raw_RMS": data[RMS_FEATURE].to_numpy(dtype=float),
            "log_RMS": data["log_RMS"].to_numpy(dtype=float),
            "z_log_RMS": data["z_log_RMS"].to_numpy(dtype=float),
        }
        scopes = [("pooled", np.ones(len(data), dtype=bool))] + [
            (channel, data["channel"].astype(str).to_numpy() == channel) for channel in CHANNELS
        ]
        for residual_type, residual in vectors.items():
            for amplitude_name, amplitude in amplitudes.items():
                for scope, mask in scopes:
                    rows.append({
                        "feature": feature, "feature_label": LABELS[feature], "model_method": method,
                        "model_variant": "adjusted", "residual_type": residual_type,
                        "residual_definition": (
                            "observed z-entropy minus fixed-effect prediction" if residual_type == "marginal"
                            else "observed z-entropy minus fixed-plus-random-effect prediction"
                        ),
                        "amplitude_variable": amplitude_name, "scope": scope,
                        **correlation_pair(residual[mask], amplitude[mask]),
                    })
    return pd.DataFrame(rows)


def evidence_classification(results: pd.DataFrame, residuals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for feature in FEATURES:
        preferred = results.loc[(results["feature"] == feature) & (results["method_role"] == "primary")]
        robust = results.loc[(results["feature"] == feature) & (results["method"] == "woman_cluster_robust_OLS")]
        p_un = preferred.loc[preferred["model_variant"] == "unadjusted"].iloc[0]
        p_ad = preferred.loc[preferred["model_variant"] == "adjusted"].iloc[0]
        r_ad = robust.loc[robust["model_variant"] == "adjusted"].iloc[0]
        cond = residuals.loc[
            (residuals["feature"] == feature) & (residuals["residual_type"] == "conditional")
            & (residuals["amplitude_variable"] == "log_RMS") & (residuals["scope"] == "pooled")
        ].iloc[0]
        if feature == "fwh_entropy_totalsampen":
            scale = "demonstrated scale-sensitive by construction in Phase 1"
        elif feature == "fwh_entropy_sampen_fixed_r":
            scale = "demonstrated numerically invariant to positive scaling in Phase 1"
        else:
            scale = "algebraically scale-invariant after singular-value normalization; not empirically rescaled in Phase 1"
        if p_ad["uc_vs_in_p_value"] < ALPHA and r_ad["uc_vs_in_p_value"] < ALPHA:
            state = "supported after adjustment by both corrected mixed and woman-cluster analyses"
        elif p_ad["uc_vs_in_p_value"] < ALPHA or r_ad["uc_vs_in_p_value"] < ALPHA:
            state = "model-dependent after adjustment"
        else:
            state = "not supported after adjustment by either corrected analysis"
        if p_ad["z_log_rms_p_value"] < ALPHA and r_ad["z_log_rms_p_value"] < ALPHA:
            amplitude = "direct log-RMS association supported by both corrected analyses"
        elif p_ad["z_log_rms_p_value"] < ALPHA or r_ad["z_log_rms_p_value"] < ALPHA:
            amplitude = "direct log-RMS association is model-dependent"
        else:
            amplitude = "no direct log-RMS association supported"
        if np.sign(p_un["uc_vs_in_estimate"]) != np.sign(p_ad["uc_vs_in_estimate"]):
            attenuation = "sign reversal after RMS adjustment"
        elif abs(p_ad["uc_vs_in_estimate"]) < abs(p_un["uc_vs_in_estimate"]):
            attenuation = "attenuated toward zero after RMS adjustment"
        else:
            attenuation = "magnitude increased away from zero after RMS adjustment"
        residual = (
            f"conditional pooled log-RMS Pearson r={cond['pearson_r']:.3f}; "
            f"Spearman rho={cond['spearman_rho']:.3f}"
        )
        rows.append({
            "feature": feature, "feature_label": LABELS[feature],
            "state_association": state, "direct_amplitude_association": amplitude,
            "mathematical_scale_sensitivity": scale, "attenuation_after_rms_adjustment": attenuation,
            "absolute_uc_coefficient_change": p_ad["absolute_uc_coefficient_change_after_adjustment"],
            "adjusted_primary_uc_estimate": p_ad["uc_vs_in_estimate"],
            "adjusted_primary_uc_p": p_ad["uc_vs_in_p_value"],
            "adjusted_woman_cluster_uc_estimate": r_ad["uc_vs_in_estimate"],
            "adjusted_woman_cluster_uc_p": r_ad["uc_vs_in_p_value"],
            "residual_association_with_rms": residual,
            "conditional_pearson_log_rms": cond["pearson_r"],
            "conditional_spearman_log_rms": cond["spearman_rho"],
        })
    return pd.DataFrame(rows)


def plot_coefficients(results: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharey=True, constrained_layout=True)
    methods = (("preferred_mixed_model", "Window-aware mixed model"), ("woman_cluster_robust_OLS", "Woman-cluster OLS"))
    colors = {"unadjusted": "#4C78A8", "adjusted": "#E45756"}
    for ax, (method, title) in zip(axes, methods):
        subset = results.loc[results["method"] == method]
        for i, feature in enumerate(FEATURES):
            for variant, offset in (("unadjusted", -0.10), ("adjusted", 0.10)):
                row = subset.loc[(subset["feature"] == feature) & (subset["model_variant"] == variant)].iloc[0]
                ax.errorbar(row["uc_vs_in_estimate"], i + offset,
                            xerr=[[row["uc_vs_in_estimate"] - row["uc_vs_in_ci_low"]], [row["uc_vs_in_ci_high"] - row["uc_vs_in_estimate"]]],
                            fmt="o", color=colors[variant], capsize=3, label=variant.title() if i == 0 else None)
        ax.axvline(0, color="0.35", linewidth=1)
        ax.set_title(title)
        ax.set_xlabel("UC minus IN coefficient (outcome SD units)")
        ax.grid(axis="x", alpha=0.25)
        ax.legend(frameon=False)
    axes[0].set_yticks(range(len(FEATURES)), [LABELS[f] for f in FEATURES])
    axes[0].invert_yaxis()
    fig.savefig(OUTPUT_DIR / "uc_in_coefficients_before_after_rms_adjustment.png", dpi=180)
    plt.close(fig)


def plot_corrected_residuals(residuals: pd.DataFrame) -> None:
    pooled = residuals.loc[(residuals["scope"] == "pooled") & (residuals["amplitude_variable"].isin(["raw_RMS", "log_RMS"]))].copy()
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    for ax, metric, label in ((axes[0], "pearson_r", "Pearson r"), (axes[1], "spearman_rho", "Spearman rho")):
        x = np.arange(len(FEATURES))
        width = 0.18
        combos = [("marginal", "raw_RMS"), ("marginal", "log_RMS"), ("conditional", "raw_RMS"), ("conditional", "log_RMS")]
        colors = ["#4C78A8", "#72B7B2", "#E45756", "#F2CF5B"]
        for j, ((resid_type, amp), color) in enumerate(zip(combos, colors)):
            vals = []
            for feature in FEATURES:
                vals.append(float(pooled.loc[(pooled["feature"] == feature) & (pooled["residual_type"] == resid_type) & (pooled["amplitude_variable"] == amp), metric].iloc[0]))
            ax.bar(x + (j - 1.5) * width, vals, width, color=color, label=f"{resid_type}, {amp}")
        ax.axhline(0, color="0.35", linewidth=1)
        ax.set_xticks(x, [LABELS[f] for f in FEATURES], rotation=18, ha="right")
        ax.set_ylabel(label)
        ax.set_title(f"Adjusted residual association: {label}")
        ax.grid(axis="y", alpha=0.25)
    axes[1].legend(frameon=False, fontsize=8, loc="best")
    fig.savefig(OUTPUT_DIR / "corrected_residual_correlations.png", dpi=180)
    plt.close(fig)


def write_versions() -> None:
    packages = ("numpy", "pandas", "scipy", "matplotlib", "statsmodels")
    lines = [f"python={platform.python_version()}", f"platform={platform.platform()}"]
    lines += [f"{package}={importlib.metadata.version(package)}" for package in packages]
    lines += [f"random_seed={RANDOM_SEED}", "corrected_mixed_fit=ML", "old_table8_fit=REML", "cluster_correction=CR1_with_df_correction_and_t_reference"]
    (OUTPUT_DIR / "package_versions.txt").write_text("\n".join(lines) + "\n")


def fmt(value: float, digits: int = 3) -> str:
    return "NA" if not np.isfinite(value) else f"{value:.{digits}g}"


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    population_data, population = load_corrected_population()
    results, convergence, selected_adjusted, scaling = corrected_models(population_data)
    residuals = corrected_residual_diagnostics(selected_adjusted)
    evidence = evidence_classification(results, residuals)
    results.to_csv(OUTPUT_DIR / "corrected_amplitude_adjusted_model_results.csv", index=False)
    residuals.to_csv(OUTPUT_DIR / "entropy_rms_residual_diagnostics.csv", index=False)
    convergence.to_csv(OUTPUT_DIR / "model_convergence_diagnostics.csv", index=False)
    evidence.to_csv(OUTPUT_DIR / "entropy_feature_evidence_classification.csv", index=False)
    scaling.to_csv(OUTPUT_DIR / "analysis_population_and_scaling.csv", index=False)
    plot_coefficients(results)
    plot_corrected_residuals(residuals)
    write_versions()
    report = ["# Corrected RMS-adjusted entropy analysis", "",
              "UC-versus-IN models account for woman and shared physical-window dependence.",
              "The historical raw-RMS residual correlation (r=0.938) is documented in RESULTS_MAP.md; it is not an input to this analysis.", "",
              f"Population: {population}", "", "## Corrected model output", "",
              "```\n" + results.to_string(index=False) + "\n```", "", "## Evidence classification", "",
              "```\n" + evidence.to_string(index=False) + "\n```"]
    (OUTPUT_DIR / "rms_adjusted_entropy_report.md").write_text("\n".join(report) + "\n")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
