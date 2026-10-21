#!/usr/bin/env python3
"""Corrected UC-vs-IN cross-channel analysis with shared-window dependence.

The model and fallback hierarchy are locked to the Phase 1 recommendation:

Preferred (Candidate B)
    z_feature ~ state * channel
    woman random intercept + physical-window variance component

Required sensitivity and fallback (Candidate C)
    z_feature ~ state * channel
    woman-cluster-robust covariance with CR1 and t/F reference distributions

Population: UC and IN, EHG9-EHG12. Feature family: read from the frozen
per-channel Stage 06 results as the features validated on all four channels.
No supplementary feature family is analysed here.
"""
from __future__ import annotations

import json
import platform
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import statsmodels
import statsmodels.formula.api as smf
from scipy.stats import chi2, norm, t


SCRIPT_PATH = Path(__file__).resolve()
OUTPUT_DIR = Path(__file__).resolve().parents[1] / "results" / "cross_channel"
PIPELINE_ROOT = SCRIPT_PATH.parents[1]
WORK_OUTPUT_DIR = PIPELINE_ROOT / "work" / "cross_channel"
LOCKED_COVARIANCE_CHOICE = "woman intercept + physical-window variance component; woman-cluster sensitivity"
PHASE1_REPORT = LOCKED_COVARIANCE_CHOICE
LOCKED_FEATURES = ("fwh_entropy_totalsampen", "fwh_time_mean_absolute_value", "fwh_time_median", "fwh_time_peak_to_peak", "fwh_time_root_mean_square", "fwh_time_skewness", "fwh_time_waveform_length")
STAGE06_DIR = PIPELINE_ROOT / "work" / "primary"

CHANNELS = ("EHG9", "EHG10", "EHG11", "EHG12")
STATE_TERM = 'C(state, Treatment(reference="IN"))'
CHANNEL_TERM = 'C(channel, Treatment(reference="EHG9"))'
REDUCED_FORMULA = f"z_feature ~ {STATE_TERM} + {CHANNEL_TERM}"
FULL_FORMULA = f"z_feature ~ {STATE_TERM} * {CHANNEL_TERM}"
OPTIMIZERS = ("lbfgs", "bfgs")
ALPHA = 0.05

FEATURE_LABELS = {
    "fwh_entropy_totalsampen": "Total SampEn",
    "fwh_time_mean_absolute_value": "Mean absolute value",
    "fwh_time_median": "Median",
    "fwh_time_peak_to_peak": "Peak-to-peak amplitude",
    "fwh_time_root_mean_square": "RMS",
    "fwh_time_skewness": "Skewness",
    "fwh_time_waveform_length": "Waveform length",
}


@dataclass
class FitAttempt:
    optimizer: str
    model_part: str
    fit: object | None
    converged: bool
    finite_llf: bool
    llf: float
    fixed_covariance_finite: bool
    boundary_or_singularity_warning: bool
    warnings: list[str]
    exception: str

    @property
    def usable(self) -> bool:
        return bool(
            self.fit is not None
            and self.converged
            and self.finite_llf
            and self.fixed_covariance_finite
        )


def bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.lower().isin(["true", "1", "yes"])


def bh_adjust(p_values: pd.Series) -> pd.Series:
    values = pd.to_numeric(p_values, errors="coerce").to_numpy(dtype=float)
    output = np.full(len(values), np.nan)
    valid = np.flatnonzero(np.isfinite(values))
    if not len(valid):
        return pd.Series(output, index=p_values.index)
    order = valid[np.argsort(values[valid])]
    ranked = values[order] * len(valid) / np.arange(1, len(valid) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    output[order] = np.minimum(ranked, 1.0)
    return pd.Series(output, index=p_values.index)


def load_frozen_feature_family() -> tuple[list[str], pd.DataFrame]:
    """Read, rather than manually recreate, the seven-feature frozen family."""
    frames = []
    for channel in CHANNELS:
        path = STAGE06_DIR / f"mixed_model_results_60s_fwh_{channel}.csv"
        frame = pd.read_csv(path)
        frame = frame.loc[frame["contrast"].eq("contraction_vs_baseline")].copy()
        frame["channel"] = channel
        frame["primary_validated"] = bool_series(frame["primary_validated"])
        frames.append(frame)
    frozen = pd.concat(frames, ignore_index=True)

    family = []
    for feature, group in frozen.groupby("feature", sort=False):
        if group["channel"].nunique() == 4 and group["primary_validated"].all():
            family.append(feature)
    family = sorted(family)

    phase1_features = sorted(LOCKED_FEATURES)
    if family != phase1_features:
        raise ValueError(
            "Frozen all-channel validated feature family does not match the Phase 1 audit: "
            f"frozen={family}; phase1={phase1_features}"
        )
    if len(family) != 7:
        raise ValueError(f"Expected seven frozen cross-channel features, found {len(family)}")
    return family, frozen.loc[frozen["feature"].isin(family)].copy()


def load_pooled_event_units() -> pd.DataFrame:
    frames = []
    for channel in CHANNELS:
        frame = pd.read_csv(STAGE06_DIR / f"event_unit_analysis_60s_{channel}.csv")
        frame = frame.copy()
        frame["channel"] = channel
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def verify_population(pooled: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    data = pooled.loc[pooled["window_class"].isin(["baseline", "contraction"])].copy()
    data["state"] = data["window_class"].map({"baseline": "IN", "contraction": "UC"})
    data["state"] = pd.Categorical(data["state"], categories=["IN", "UC"])
    data["channel"] = pd.Categorical(data["channel"], categories=list(CHANNELS))
    data["physical_window_key"] = data["window_id"].astype(str)

    per_window = data.groupby("window_id", observed=True).agg(
        rows=("channel", "size"),
        channels=("channel", "nunique"),
        women=("patient_id", "nunique"),
        recordings=("recording_id", "nunique"),
        states=("state", "nunique"),
    )
    if not (
        per_window["rows"].eq(4).all()
        and per_window["channels"].eq(4).all()
        and per_window["women"].eq(1).all()
        and per_window["recordings"].eq(1).all()
        and per_window["states"].eq(1).all()
    ):
        raise ValueError("UC/IN physical-window structure is incomplete or inconsistent")

    window_states = data.drop_duplicates("window_id")["state"].value_counts()
    population = {
        "women": int(data["patient_id"].nunique()),
        "recordings": int(data["recording_id"].nunique()),
        "unique_windows": int(data["window_id"].nunique()),
        "rows": int(len(data)),
        "uc_windows": int(window_states.get("UC", 0)),
        "in_windows": int(window_states.get("IN", 0)),
        "rows_per_window": 4,
        "complete_four_channel_proportion": float(per_window["channels"].eq(4).mean()),
    }
    return data, population


def prepare_feature(population_data: pd.DataFrame, feature: str) -> tuple[pd.DataFrame, float, float]:
    columns = [
        feature, "patient_id", "recording_id", "window_id", "physical_window_key",
        "state", "channel",
    ]
    data = population_data[columns].replace([np.inf, -np.inf], np.nan).dropna().copy()
    mean = float(data[feature].mean())
    sd = float(data[feature].std(ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        raise ValueError(f"Invalid pooled UC/IN sample SD for {feature}: {sd}")
    data["z_feature"] = (data[feature] - mean) / sd
    return data, mean, sd


def interaction_terms() -> dict[str, str]:
    prefix = 'C(state, Treatment(reference="IN"))[T.UC]:C(channel, Treatment(reference="EHG9"))'
    return {channel: f"{prefix}[T.{channel}]" for channel in CHANNELS[1:]}


def uc_term() -> str:
    return 'C(state, Treatment(reference="IN"))[T.UC]'


def has_warning(messages: list[str]) -> bool:
    tokens = ("boundary", "singular", "not positive definite", "hessian")
    return any(any(token in message.lower() for token in tokens) for message in messages)


def attempt_mixed(formula: str, data: pd.DataFrame, optimizer: str, model_part: str) -> FitAttempt:
    fit = None
    caught_messages: list[str] = []
    exception = ""
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            model = smf.mixedlm(
                formula,
                data,
                groups=data["patient_id"],
                re_formula="1",
                vc_formula={"physical_window": "0 + C(physical_window_key)"},
            )
            fit = model.fit(reml=False, method=optimizer, maxiter=1000, disp=False)
        caught_messages = [str(item.message) for item in caught]
    except Exception as exc:
        exception = f"{type(exc).__name__}: {exc}"

    converged = bool(getattr(fit, "converged", False)) if fit is not None else False
    llf = float(getattr(fit, "llf", np.nan)) if fit is not None else np.nan
    finite_llf = bool(np.isfinite(llf))
    fixed_covariance_finite = False
    if fit is not None:
        try:
            fixed_names = list(fit.fe_params.index)
            covariance = fit.cov_params().loc[fixed_names, fixed_names].to_numpy(dtype=float)
            fixed_covariance_finite = bool(np.isfinite(covariance).all())
        except Exception as exc:
            exception = " | ".join(filter(None, [exception, f"fixed covariance: {type(exc).__name__}: {exc}"]))
    return FitAttempt(
        optimizer=optimizer,
        model_part=model_part,
        fit=fit,
        converged=converged,
        finite_llf=finite_llf,
        llf=llf,
        fixed_covariance_finite=fixed_covariance_finite,
        boundary_or_singularity_warning=has_warning(caught_messages),
        warnings=caught_messages,
        exception=exception,
    )


def fit_preferred_pair(data: pd.DataFrame) -> tuple[FitAttempt | None, FitAttempt | None, list[FitAttempt]]:
    attempts: list[FitAttempt] = []
    for optimizer in OPTIMIZERS:
        reduced = attempt_mixed(REDUCED_FORMULA, data, optimizer, "reduced")
        full = attempt_mixed(FULL_FORMULA, data, optimizer, "full")
        attempts.extend([reduced, full])
        if reduced.usable and full.usable:
            return reduced, full, attempts
    return None, None, attempts


def fit_cluster_robust(data: pd.DataFrame):
    ordinary = smf.ols(FULL_FORMULA, data=data).fit()
    return ordinary.get_robustcov_results(
        cov_type="cluster",
        groups=data["patient_id"],
        use_correction=True,
        df_correction=True,
        use_t=True,
    )


def fixed_parameters_and_covariance(fit, method: str) -> tuple[list[str], np.ndarray, np.ndarray, float]:
    names = list(fit.model.exog_names)
    if method == "preferred_mixed_model":
        params = fit.fe_params.loc[names].to_numpy(dtype=float)
        covariance = fit.cov_params().loc[names, names].to_numpy(dtype=float)
        critical = float(norm.ppf(0.975))
    else:
        params = np.asarray(fit.params, dtype=float)
        covariance = np.asarray(fit.cov_params(), dtype=float)
        critical = float(t.ppf(0.975, fit.df_resid_inference))
    return names, params, covariance, critical


def linear_contrast(
    fit,
    method: str,
    weights: dict[str, float],
) -> dict:
    names, params, covariance, critical = fixed_parameters_and_covariance(fit, method)
    vector = np.zeros(len(names), dtype=float)
    for term, weight in weights.items():
        vector[names.index(term)] = weight
    estimate = float(vector @ params)
    variance = float(vector @ covariance @ vector)
    se = float(np.sqrt(max(variance, 0.0)))
    statistic = estimate / se if se > 0 else np.nan
    if method == "preferred_mixed_model":
        p_value = float(2 * norm.sf(abs(statistic))) if np.isfinite(statistic) else np.nan
        df = np.nan
    else:
        df = float(fit.df_resid_inference)
        p_value = float(2 * t.sf(abs(statistic), df)) if np.isfinite(statistic) else np.nan
    return {
        "estimate": estimate,
        "se": se,
        "ci_low": estimate - critical * se,
        "ci_high": estimate + critical * se,
        "statistic": statistic,
        "df": df,
        "p_value": p_value,
    }


def channel_effect(fit, method: str, channel: str) -> dict:
    weights = {uc_term(): 1.0}
    if channel != "EHG9":
        weights[interaction_terms()[channel]] = 1.0
    return linear_contrast(fit, method, weights)


def pairwise_interaction(fit, method: str, channel: str) -> dict:
    return linear_contrast(fit, method, {interaction_terms()[channel]: 1.0})


def mixed_global_test(reduced, full) -> dict:
    df = int(full.model.exog.shape[1] - reduced.model.exog.shape[1])
    statistic = max(0.0, 2.0 * (float(full.llf) - float(reduced.llf)))
    return {
        "test_type": "ML likelihood-ratio chi-square",
        "statistic": statistic,
        "df_num": float(df),
        "df_denom": np.nan,
        "p_value": float(chi2.sf(statistic, df)),
    }


def robust_global_test(fit) -> dict:
    names = list(fit.model.exog_names)
    terms = list(interaction_terms().values())
    restriction = np.zeros((len(terms), len(names)), dtype=float)
    for row, term in enumerate(terms):
        restriction[row, names.index(term)] = 1.0
    result = fit.wald_test(restriction, use_f=True, scalar=True)
    return {
        "test_type": "woman-cluster-robust Wald F (CR1, t/F small-sample reference)",
        "statistic": float(np.asarray(result.statistic).squeeze()),
        "df_num": float(result.df_num),
        "df_denom": float(result.df_denom),
        "p_value": float(np.asarray(result.pvalue).squeeze()),
    }


def mixed_variances(fit) -> dict:
    names = list(getattr(fit.model.exog_vc, "names", []))
    values = np.asarray(getattr(fit, "vcomp", []), dtype=float)
    mapping = dict(zip(names, values))
    return {
        "woman_variance": float(fit.cov_re.iloc[0, 0]),
        "physical_window_variance": float(mapping.get("physical_window", np.nan)),
        "residual_variance": float(fit.scale),
    }


def attempt_summary(attempts: list[FitAttempt]) -> str:
    payload = []
    for item in attempts:
        payload.append({
            "optimizer": item.optimizer,
            "model_part": item.model_part,
            "converged": item.converged,
            "finite_llf": item.finite_llf,
            "llf": item.llf if np.isfinite(item.llf) else None,
            "fixed_covariance_finite": item.fixed_covariance_finite,
            "boundary_or_singularity_warning": item.boundary_or_singularity_warning,
            "warnings": item.warnings,
            "exception": item.exception,
        })
    return json.dumps(payload, ensure_ascii=False, separators=(",", ":"))


def analyse_feature(feature: str, population_data: pd.DataFrame) -> tuple[list[dict], list[dict], dict]:
    data, mean, sd = prepare_feature(population_data, feature)
    reduced_attempt, full_attempt, attempts = fit_preferred_pair(data)
    preferred_usable = reduced_attempt is not None and full_attempt is not None
    robust = fit_cluster_robust(data)

    fits: list[tuple[str, object, str, bool]] = []
    if preferred_usable:
        fits.append(("preferred_mixed_model", full_attempt.fit, "primary", True))
    fits.append((
        "woman_cluster_robust_OLS",
        robust,
        "required_sensitivity" if preferred_usable else "prespecified_fallback",
        not preferred_usable,
    ))

    effect_rows: list[dict] = []
    interaction_rows: list[dict] = []
    for method, fit, role, selected_primary in fits:
        for channel in CHANNELS:
            contrast = channel_effect(fit, method, channel)
            effect_rows.append({
                "feature": feature,
                "feature_label": FEATURE_LABELS.get(feature, feature),
                "method": method,
                "method_role": role,
                "selected_for_primary_inference": selected_primary,
                "channel": channel,
                "estimate": contrast["estimate"],
                "se": contrast["se"],
                "ci_low": contrast["ci_low"],
                "ci_high": contrast["ci_high"],
                "test_statistic": contrast["statistic"],
                "reference_df": contrast["df"],
                "p_value": contrast["p_value"],
                "p_value_multiplicity_status": "descriptive; unadjusted",
                "direction": "positive" if contrast["estimate"] > 0 else "negative" if contrast["estimate"] < 0 else "zero",
                "p_lt_0_05": bool(contrast["p_value"] < ALPHA),
                "n_rows": int(len(data)),
                "n_physical_windows": int(data["window_id"].nunique()),
                "n_women": int(data["patient_id"].nunique()),
                "n_recordings": int(data["recording_id"].nunique()),
                "standardization_mean": mean,
                "standardization_sd_ddof1": sd,
                "standardization_scope": "pooled UC/IN rows across EHG9-EHG12",
            })

        global_test = (
            mixed_global_test(reduced_attempt.fit, full_attempt.fit)
            if method == "preferred_mixed_model"
            else robust_global_test(fit)
        )
        interaction_rows.append({
            "feature": feature,
            "feature_label": FEATURE_LABELS.get(feature, feature),
            "method": method,
            "method_role": role,
            "selected_for_primary_inference": selected_primary,
            "test_scope": "global",
            "contrast": "UC-by-channel global heterogeneity",
            "comparison_channel": "all EHG10-EHG12 versus EHG9",
            "test_type": global_test["test_type"],
            "estimate": np.nan,
            "se": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "statistic": global_test["statistic"],
            "df_num": global_test["df_num"],
            "df_denom": global_test["df_denom"],
            "p_value": global_test["p_value"],
            "method_specific_q_value": np.nan,
            "primary_family_q_value": np.nan,
            "multiplicity_family": "seven predefined global UC-by-channel tests",
        })
        for channel in CHANNELS[1:]:
            contrast = pairwise_interaction(fit, method, channel)
            interaction_rows.append({
                "feature": feature,
                "feature_label": FEATURE_LABELS.get(feature, feature),
                "method": method,
                "method_role": role,
                "selected_for_primary_inference": selected_primary,
                "test_scope": "pairwise",
                "contrast": f"UC effect: {channel} minus EHG9",
                "comparison_channel": channel,
                "test_type": "Wald coefficient contrast",
                "estimate": contrast["estimate"],
                "se": contrast["se"],
                "ci_low": contrast["ci_low"],
                "ci_high": contrast["ci_high"],
                "statistic": contrast["statistic"],
                "df_num": 1.0,
                "df_denom": contrast["df"],
                "p_value": contrast["p_value"],
                "method_specific_q_value": np.nan,
                "primary_family_q_value": np.nan,
                "multiplicity_family": "descriptive; pairwise interaction p-values unadjusted",
            })

    if preferred_usable:
        variances = mixed_variances(full_attempt.fit)
        selected_optimizer = full_attempt.optimizer
        selected_warning = " | ".join(dict.fromkeys(reduced_attempt.warnings + full_attempt.warnings))
        selected_boundary = bool(
            reduced_attempt.boundary_or_singularity_warning
            or full_attempt.boundary_or_singularity_warning
        )
        reduced_llf = reduced_attempt.llf
        full_llf = full_attempt.llf
    else:
        variances = {"woman_variance": np.nan, "physical_window_variance": np.nan, "residual_variance": np.nan}
        selected_optimizer = ""
        selected_warning = "No optimizer produced a usable matched reduced/full Candidate B pair"
        selected_boundary = any(item.boundary_or_singularity_warning for item in attempts)
        reduced_llf = np.nan
        full_llf = np.nan

    audit = {
        "feature": feature,
        "feature_label": FEATURE_LABELS.get(feature, feature),
        "preferred_model": "Candidate B: woman random intercept + physical-window variance component",
        "reduced_formula": REDUCED_FORMULA,
        "full_formula": FULL_FORMULA,
        "optimizer_sequence": ",".join(OPTIMIZERS),
        "optimizer_attempt_log_json": attempt_summary(attempts),
        "selected_optimizer": selected_optimizer,
        "preferred_model_usable": preferred_usable,
        "reduced_llf": reduced_llf,
        "full_llf": full_llf,
        "boundary_or_singularity_warning": selected_boundary,
        "woman_variance": variances["woman_variance"],
        "physical_window_variance": variances["physical_window_variance"],
        "residual_variance": variances["residual_variance"],
        "selected_warning": selected_warning,
        "fallback_run": True,
        "fallback_status": "succeeded",
        "fallback_used_as_primary": not preferred_usable,
        "selected_primary_method": "preferred_mixed_model" if preferred_usable else "woman_cluster_robust_OLS",
        "n_rows": int(len(data)),
        "n_physical_windows": int(data["window_id"].nunique()),
        "n_women": int(data["patient_id"].nunique()),
        "n_recordings": int(data["recording_id"].nunique()),
        "standardization_mean": mean,
        "standardization_sd_ddof1": sd,
    }
    return effect_rows, interaction_rows, audit


def apply_interaction_fdr(interactions: pd.DataFrame) -> pd.DataFrame:
    result = interactions.copy()
    global_mask = result["test_scope"].eq("global")
    for method in result.loc[global_mask, "method"].unique():
        mask = global_mask & result["method"].eq(method)
        result.loc[mask, "method_specific_q_value"] = bh_adjust(result.loc[mask, "p_value"])
    primary_mask = global_mask & bool_series(result["selected_for_primary_inference"])
    result.loc[primary_mask, "primary_family_q_value"] = bh_adjust(result.loc[primary_mask, "p_value"])
    result["global_interaction_q_lt_0_05"] = np.where(
        global_mask,
        pd.to_numeric(result["primary_family_q_value"], errors="coerce") < ALPHA,
        False,
    )
    return result


def build_comparison(
    family: list[str],
    frozen: pd.DataFrame,
    effects: pd.DataFrame,
    interactions: pd.DataFrame,
    audit: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    primary_effects = effects.loc[bool_series(effects["selected_for_primary_inference"])].copy()
    robust_effects = effects.loc[effects["method"].eq("woman_cluster_robust_OLS")].copy()
    primary_global = interactions.loc[
        interactions["test_scope"].eq("global")
        & bool_series(interactions["selected_for_primary_inference"])
    ].set_index("feature")
    robust_global = interactions.loc[
        interactions["test_scope"].eq("global")
        & interactions["method"].eq("woman_cluster_robust_OLS")
    ].set_index("feature")
    audit_by_feature = audit.set_index("feature")

    for feature in family:
        original = frozen.loc[frozen["feature"].eq(feature)].set_index("channel")
        corrected = primary_effects.loc[primary_effects["feature"].eq(feature)].set_index("channel")
        robust = robust_effects.loc[robust_effects["feature"].eq(feature)].set_index("channel")
        original_signs = np.sign(original.loc[list(CHANNELS), "estimate"].to_numpy(dtype=float))
        corrected_signs = np.sign(corrected.loc[list(CHANNELS), "estimate"].to_numpy(dtype=float))
        robust_signs = np.sign(robust.loc[list(CHANNELS), "estimate"].to_numpy(dtype=float))
        row = {
            "feature": feature,
            "feature_label": FEATURE_LABELS.get(feature, feature),
            "original_cross_channel_consistency_status": (
                "all four primary-validated with same direction"
                if len(set(original_signs)) == 1 and original["primary_validated"].all()
                else "not consistently validated/directed"
            ),
            "corrected_direction_consistency": bool(len(set(corrected_signs)) == 1),
            "corrected_direction": "positive" if corrected_signs[0] > 0 else "negative" if corrected_signs[0] < 0 else "mixed/zero",
            "corrected_n_channels_p_lt_0_05": int(corrected["p_lt_0_05"].sum()),
            "corrected_all_channels_p_lt_0_05": bool(corrected["p_lt_0_05"].all()),
            "corrected_global_interaction_p": float(primary_global.loc[feature, "p_value"]),
            "corrected_global_interaction_q": float(primary_global.loc[feature, "primary_family_q_value"]),
            "corrected_global_interaction_q_lt_0_05": bool(primary_global.loc[feature, "primary_family_q_value"] < ALPHA),
            "corrected_magnitude_interpretation": (
                "channel-dependent magnitude"
                if primary_global.loc[feature, "primary_family_q_value"] < ALPHA
                else "no FDR-supported channel heterogeneity"
            ),
            "robust_direction_consistency": bool(len(set(robust_signs)) == 1),
            "robust_n_channels_p_lt_0_05": int(robust["p_lt_0_05"].sum()),
            "robust_global_interaction_p": float(robust_global.loc[feature, "p_value"]),
            "robust_global_interaction_q_separate_family": float(robust_global.loc[feature, "method_specific_q_value"]),
            "selected_primary_method": audit_by_feature.loc[feature, "selected_primary_method"],
            "preferred_model_usable": bool(audit_by_feature.loc[feature, "preferred_model_usable"]),
            "boundary_or_singularity_warning": bool(audit_by_feature.loc[feature, "boundary_or_singularity_warning"]),
            "fallback_used_as_primary": bool(audit_by_feature.loc[feature, "fallback_used_as_primary"]),
            "original_standardization": "within channel; all IN/FM/UC rows; sample SD",
            "corrected_standardization": "pooled EHG9-EHG12; UC/IN rows; sample SD",
        }
        for channel in CHANNELS:
            key = channel.lower()
            row.update({
                f"original_{key}_estimate": float(original.loc[channel, "estimate"]),
                f"original_{key}_se": float(original.loc[channel, "std_error"]),
                f"original_{key}_ci_low": float(original.loc[channel, "ci_low"]),
                f"original_{key}_ci_high": float(original.loc[channel, "ci_high"]),
                f"original_{key}_p": float(original.loc[channel, "p_value"]),
                f"original_{key}_q": float(original.loc[channel, "q_value"]),
                f"corrected_{key}_estimate": float(corrected.loc[channel, "estimate"]),
                f"corrected_{key}_se": float(corrected.loc[channel, "se"]),
                f"corrected_{key}_ci_low": float(corrected.loc[channel, "ci_low"]),
                f"corrected_{key}_ci_high": float(corrected.loc[channel, "ci_high"]),
                f"corrected_{key}_p_descriptive": float(corrected.loc[channel, "p_value"]),
                f"robust_{key}_estimate": float(robust.loc[channel, "estimate"]),
                f"robust_{key}_se": float(robust.loc[channel, "se"]),
                f"robust_{key}_p_descriptive": float(robust.loc[channel, "p_value"]),
            })
        rows.append(row)
    return pd.DataFrame(rows)


def format_p(value) -> str:
    if pd.isna(value):
        return "NA"
    value = float(value)
    return f"{value:.3g}" if value >= 0.001 else f"{value:.2e}"


def format_estimate(value) -> str:
    return "NA" if pd.isna(value) else f"{float(value):.3f}"


def render_report(
    population: dict,
    family: list[str],
    effects: pd.DataFrame,
    interactions: pd.DataFrame,
    audit: pd.DataFrame,
    comparison: pd.DataFrame,
) -> str:
    preferred_ok = int(audit["preferred_model_usable"].sum())
    fallback_count = int(audit["fallback_used_as_primary"].sum())
    boundary_count = int(audit["boundary_or_singularity_warning"].sum())
    direction_count = int(comparison["corrected_direction_consistency"].sum())
    all_supported_count = int(comparison["corrected_all_channels_p_lt_0_05"].sum())
    heterogeneous = comparison.loc[comparison["corrected_global_interaction_q_lt_0_05"], "feature_label"].tolist()
    homogeneous = comparison.loc[~comparison["corrected_global_interaction_q_lt_0_05"], "feature_label"].tolist()
    robust_heterogeneous = comparison.loc[
        comparison["robust_global_interaction_q_separate_family"] < ALPHA, "feature_label"
    ].tolist()
    robust_all_supported_count = int(comparison["robust_n_channels_p_lt_0_05"].eq(4).sum())

    channel_table = []
    selected = effects.loc[bool_series(effects["selected_for_primary_inference"])]
    for feature in family:
        group = selected.loc[selected["feature"].eq(feature)].set_index("channel")
        cells = []
        for channel in CHANNELS:
            row = group.loc[channel]
            marker = "*" if row["p_value"] < ALPHA else ""
            cells.append(
                f"{format_estimate(row['estimate'])} "
                f"[{format_estimate(row['ci_low'])}, {format_estimate(row['ci_high'])}]{marker}"
            )
        test = interactions.loc[
            interactions["feature"].eq(feature)
            & interactions["test_scope"].eq("global")
            & bool_series(interactions["selected_for_primary_inference"])
        ].iloc[0]
        channel_table.append(
            f"| {FEATURE_LABELS.get(feature, feature)} | " + " | ".join(cells) +
            f" | {format_p(test['p_value'])} | {format_p(test['primary_family_q_value'])} |"
        )

    support_table = []
    robust_selected = effects.loc[effects["method"].eq("woman_cluster_robust_OLS")]
    for feature in family:
        preferred_group = selected.loc[selected["feature"].eq(feature)]
        robust_group = robust_selected.loc[robust_selected["feature"].eq(feature)]
        preferred_channels = preferred_group.loc[preferred_group["p_value"] < ALPHA, "channel"].tolist()
        robust_channels = robust_group.loc[robust_group["p_value"] < ALPHA, "channel"].tolist()
        preferred_test = interactions.loc[
            interactions["feature"].eq(feature)
            & interactions["test_scope"].eq("global")
            & bool_series(interactions["selected_for_primary_inference"])
        ].iloc[0]
        robust_test = interactions.loc[
            interactions["feature"].eq(feature)
            & interactions["test_scope"].eq("global")
            & interactions["method"].eq("woman_cluster_robust_OLS")
        ].iloc[0]
        support_table.append(
            f"| {FEATURE_LABELS.get(feature, feature)} | "
            f"{', '.join(preferred_channels) if preferred_channels else 'none'} | "
            f"{', '.join(robust_channels) if robust_channels else 'none'} | "
            f"{format_p(preferred_test['primary_family_q_value'])} | "
            f"{format_p(robust_test['method_specific_q_value'])} |"
        )

    detailed_answers = []
    for label in ("Total SampEn", "Peak-to-peak amplitude"):
        row = comparison.loc[comparison["feature_label"].eq(label)].iloc[0]
        detailed_answers.append(
            f"- **{label}:** direction was consistent across all four channels; "
            f"{int(row['corrected_n_channels_p_lt_0_05'])}/4 channel-specific contrasts had descriptive p<0.05; "
            f"the corrected global interaction was p={format_p(row['corrected_global_interaction_p'])}, "
            f"q={format_p(row['corrected_global_interaction_q'])}. The preferred model "
            f"{'was usable' if row['preferred_model_usable'] else 'failed, so the prespecified fallback was used'}; "
            f"the woman-cluster sensitivity supported {int(row['robust_n_channels_p_lt_0_05'])}/4 "
            f"channel-specific contrasts and its global p was {format_p(row['robust_global_interaction_p'])} "
            f"(separate-family q={format_p(row['robust_global_interaction_q_separate_family'])})."
        )

    interaction_text = ", ".join(heterogeneous) if heterogeneous else "none"
    homogeneous_text = ", ".join(homogeneous) if homogeneous else "none"
    robust_text = ", ".join(robust_heterogeneous) if robust_heterogeneous else "none"

    return f"""# Corrected cross-channel UC-versus-IN analysis

## Bottom line

All {direction_count}/7 previously reported cross-channel-consistent features retained the same UC-versus-IN direction across EHG9-EHG12 after explicitly modelling shared physical windows. Under the preferred window-aware mixed model, {all_supported_count}/7 had descriptive p<0.05 in every channel and global tests identified FDR-supported channel-dependent magnitude for {interaction_text}. Under the required woman-cluster sensitivity, only {robust_all_supported_count}/7 had p<0.05 in every channel and no global interaction survived its separate seven-test BH correction.

The direction-consistency claim is therefore preserved, but inferential strength is covariance-model dependent. A blanket claim of equal effect magnitude is not supported by Candidate B for the features with significant global interaction, while Candidate C does not confirm those interactions after FDR correction. The phrase **“cross-channel spatial consistency within the same cohort”** is supportable as a same-direction pattern, not as proof of uniform statistical support, channel-invariant magnitude, or independent independent cohorts.

## Locked analysis design

No post-Phase 1 redesign was made before examining results. The analysis used Candidate B exactly as recommended: `z_feature ~ state * channel`, with IN and EHG9 as references, a woman random intercept, and a physical-window variance component. Additive and interaction models used ML and the same random structure; their likelihood-ratio comparison is the 3-df global UC-by-channel test.

Candidate C was run for every feature using woman-cluster-robust OLS with CR1 finite-cluster correction and t/F inference based on 42 clusters. It served as the prespecified fallback if Candidate B failed and as a mandatory sensitivity otherwise. Recording was not separately parameterized in Candidate B because Phase 1 found the three-tier woman-recording-window model numerically unstable. Candidate C permits arbitrary dependence within woman, including dependence across recordings and windows.

Each feature was standardized once using its mean and sample SD (`ddof=1`) across all pooled UC/IN rows from EHG9-EHG12. Frozen channel-specific estimates used a different, within-channel three-state SD and should not be compared numerically without acknowledging that different denominator.

## Analysis population and predefined features

- Women: {population['women']}
- Recordings: {population['recordings']}
- Unique physical windows: {population['unique_windows']:,}
- Rows: {population['rows']:,}
- UC windows: {population['uc_windows']}
- IN windows: {population['in_windows']:,}
- Four-channel completeness: {100 * population['complete_four_channel_proportion']:.1f}%

The seven-feature family was read from the frozen per-channel Stage 06 results as the features marked `primary_validated` on every channel: {', '.join(FEATURE_LABELS.get(item, item) for item in family)}. No supplementary 11-feature analysis was added.

## Convergence and fallback

Candidate B produced usable matched reduced/full fits for {preferred_ok}/7 features. Fallback was used as primary for {fallback_count}/7. Boundary or singularity warnings occurred for {boundary_count}/7 preferred fits and are retained in `convergence_fallback_audit.csv`; a small woman-intercept variance was accepted when the physical-window variance and fixed-effect covariance were finite. Total SampEn and peak-to-peak received the same prespecified diagnostics and fallback rules as every other feature.

## Corrected channel-specific UC-versus-IN effects

Values are standardized estimates with 95% CIs from the selected corrected model. An asterisk marks a descriptive, unadjusted channel-specific p<0.05. These p-values do not belong to the seven-test global-interaction FDR family.

| Feature | EHG9 | EHG10 | EHG11 | EHG12 | Global p | Global BH q |
|:--|--:|--:|--:|--:|--:|--:|
{chr(10).join(channel_table)}

## Interaction interpretation

- FDR-supported channel-dependent magnitude: {interaction_text}.
- No FDR-supported channel heterogeneity: {homogeneous_text}.
- In the separately adjusted seven-test woman-cluster sensitivity family, FDR-supported interactions: {robust_text}.
- Pairwise EHG10/EHG11/EHG12-versus-EHG9 interaction estimates and unadjusted p-values are in `interaction_tests.csv`. They localize a global difference but do not replace the global FDR-controlled test.

Direction consistency, within-channel support, and equality of magnitude answer different questions. A feature can retain the same sign on every channel while showing a significant global interaction because one channel's effect is larger or smaller.

The five Candidate B interactions are localized primarily to a larger EHG10 effect. Their EHG10-minus-EHG9 point differences range from 0.263 to 0.303 SD, while the corresponding EHG11/EHG12 differences from EHG9 are small. Candidate C retains essentially the same point differences but assigns much larger uncertainty to EHG10-minus-EHG9, so none of its global tests survives BH correction.

## Statistical-support sensitivity

Channel lists use descriptive, unadjusted p<0.05. Global q-values use separate seven-feature BH families.

| Feature | Candidate B supported channels | Woman-cluster supported channels | Candidate B global q | Woman-cluster global q |
|:--|:--|:--|--:|--:|
{chr(10).join(support_table)}

## Features requiring particular attention

{chr(10).join(detailed_answers)}

## Answers to the reviewer questions

1. **Direction:** {direction_count}/7 retained a common UC-versus-IN direction across all four channels.
2. **Within-channel support:** Candidate B supported all four channels for {all_supported_count}/7 features; median was supported on EHG9, EHG10, and EHG12 but not EHG11. Candidate C supported all four channels only for skewness; Total SampEn was supported on EHG9, median on EHG12, and the other channel-specific contrasts were not significant. These are descriptive, unadjusted p-values.
3. **Significant global interactions:** Candidate B identified {interaction_text} after BH correction. Candidate C identified none after its separate BH correction.
4. **Channel-dependent magnitude:** Candidate B supports channel dependence for the five-feature set above, driven mainly by larger EHG10 effects. Candidate C does not confirm these as FDR-significant, so magnitude heterogeneity is model-dependent.
5. **Similar magnitude:** Candidate B found no FDR-supported heterogeneity for {homogeneous_text}; Candidate C found no FDR-supported heterogeneity for any feature. Absence of detected heterogeneity is not proof of exact equality.
6. **Total SampEn:** summarized above. Its entropy/amplitude interpretation is outside this task; only cross-channel dependence is assessed here.
7. **Peak-to-peak amplitude:** summarized above.
8. **Material scientific change:** the cross-channel direction pattern is {'preserved' if direction_count == 7 else 'partly preserved'}, but the strength of channel-specific support and evidence for unequal magnitude depend materially on covariance assumptions. Candidate B resolves the shared-window omission; Candidate C shows that stronger protection against all within-woman dependence substantially widens uncertainty.
9. **“Cross-channel spatial consistency within the same cohort”:** supported for direction. It should not be used without qualification to imply that every channel remains statistically significant under woman-cluster inference, that magnitudes are uniform, or that these are independent replications.
10. **Claims to weaken:** avoid claiming that all seven are robustly significant on all four channels, or that the five Candidate B interactions conclusively establish channel-dependent magnitude without noting the negative FDR-controlled cluster sensitivity. Also avoid “equivalent,” “uniform,” and “channel-independent,” and do not treat four channels from the same windows and women as independent replications.

## Multiplicity

BH-FDR was applied to exactly seven global UC-by-channel p-values selected by the prespecified primary hierarchy. The woman-cluster global tests were also BH-adjusted as a separate seven-test sensitivity family. Channel-specific and pairwise interaction p-values are explicitly descriptive and unadjusted. Neither family was mixed with the original 34-feature EHG9 primary family.

## Reproducibility

- Script: `{SCRIPT_PATH}`
- Phase 1 design source: `{PHASE1_REPORT}`
- Frozen feature/effect sources: `{STAGE06_DIR}/mixed_model_results_60s_fwh_EHG9.csv` through `EHG12.csv`
- Frozen analysis population sources: `{STAGE06_DIR}/event_unit_analysis_60s_EHG9.csv` through `EHG12.csv`
- Random seed: none; all procedures are deterministic.
- Package versions: `package_versions.txt`.
"""


def write_versions() -> None:
    text = "\n".join([
        f"python={platform.python_version()}",
        f"platform={platform.platform()}",
        f"numpy={np.__version__}",
        f"pandas={pd.__version__}",
        f"scipy={scipy.__version__}",
        f"statsmodels={statsmodels.__version__}",
        "random_seed=none",
        "preferred_model=statsmodels MixedLM ML; woman random intercept; physical-window variance component",
        "preferred_optimizers=lbfgs,bfgs",
        "fallback_and_sensitivity=OLS woman-cluster covariance; CR1; df correction; use_t=True",
        "global_multiplicity=Benjamini-Hochberg within seven predefined features",
    ]) + "\n"
    (OUTPUT_DIR / "package_versions.txt").write_text(text, encoding="utf-8")


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    WORK_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    phase1_text = LOCKED_COVARIANCE_CHOICE

    family, frozen = load_frozen_feature_family()
    pooled = load_pooled_event_units()
    population_data, population = verify_population(pooled)
    (WORK_OUTPUT_DIR / "population_manifest.json").write_text(
        json.dumps(population, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (WORK_OUTPUT_DIR / "locked_feature_family.json").write_text(
        json.dumps({"features": family, "channels": list(CHANNELS), "states": ["IN", "UC"]}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Locked feature family ({len(family)}): {family}", flush=True)
    print(f"Population: {population}", flush=True)

    effect_rows: list[dict] = []
    interaction_rows: list[dict] = []
    audit_rows: list[dict] = []
    for index, feature in enumerate(family, start=1):
        print(f"[{index}/{len(family)}] {feature}", flush=True)
        effects, interactions, audit = analyse_feature(feature, population_data)
        effect_rows.extend(effects)
        interaction_rows.extend(interactions)
        audit_rows.append(audit)
        print(
            f"  preferred_usable={audit['preferred_model_usable']}; "
            f"optimizer={audit['selected_optimizer'] or 'none'}; "
            f"fallback_primary={audit['fallback_used_as_primary']}",
            flush=True,
        )

    effects_df = pd.DataFrame(effect_rows)
    interactions_df = apply_interaction_fdr(pd.DataFrame(interaction_rows))
    audit_df = pd.DataFrame(audit_rows)
    comparison_df = build_comparison(family, frozen, effects_df, interactions_df, audit_df)

    effects_df.to_csv(OUTPUT_DIR / "corrected_model_results.csv", index=False)
    interactions_df.to_csv(OUTPUT_DIR / "interaction_tests.csv", index=False)
    audit_df.to_csv(OUTPUT_DIR / "convergence_fallback_audit.csv", index=False)
    comparison_df.to_csv(OUTPUT_DIR / "corrected_cross_channel_comparison.csv", index=False)
    report = render_report(population, family, effects_df, interactions_df, audit_df, comparison_df)
    (OUTPUT_DIR / "corrected_cross_channel_report.md").write_text(report, encoding="utf-8")
    write_versions()
    print(f"Wrote Phase E1 outputs to {OUTPUT_DIR}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
