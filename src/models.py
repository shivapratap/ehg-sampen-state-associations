"""Stage 06 -- Mixed-effects models, UC/FM vs IN, per channel.

Standalone port of pipeline/module2_5/{prepare_analysis.py, run_models.py}
into this project's numbered-stage layout, extended across all four Paper 1
channels (EHG9, EHG10, EHG11, EHG12). Logic is ported, not imported: the
feature-exclusion registry, event-unit selection, model formula, optimizer
fallback hierarchy, standardization, and FDR-correction-per-contrast-family
are unchanged from the frozen module2_5 design. See that module's own
docstring (module2_5/run_models.py) for the full statistical rationale
(reproduced in condensed form below); this file focuses on what is
Paper-1-specific: running the same per-feature model independently for
each of the four channels rather than only EHG9.

Primary model (per feature, per channel)
-----------------------------------------
    Feature(z) ~ Class + (1 | Woman)         [+ recording-level variance
                                                 component, see below]

    groups     = patient_id
    vc_formula = {"recording": "0 + C(recording_id)"}   (nested tier only)

Class is a 3-level categorical (baseline/fetal_movement/contraction,
"baseline" = reference), giving three contrasts per fit:
    contraction_vs_baseline        (PRIMARY confirmatory contrast)
    fetal_movement_vs_baseline     (secondary, descriptive)
    contraction_vs_fetal_movement  (secondary, specificity; derived post-hoc
                                     from the same fit's coefficient
                                     covariance, not a separate model)

Each feature is z-scored (its own per-feature mean/SD) before fitting;
estimates are reported both standardized and back-transformed to raw
units. Optimizer/structure fallback: nested+lbfgs -> nested+bfgs ->
woman-only+lbfgs -> woman-only+bfgs (first clean fit wins; boundary fits
are retained, not discarded, and flagged). Benjamini-Hochberg FDR is
applied separately per contrast family PER CHANNEL (this project's
extension of the original single-channel design: each channel's four
extra recording sites are a separate signal-acquisition condition, so
correcting across channels as one pooled family would conflate four
different physical measurements into one hypothesis-testing family).
`validated` (q < 0.05) is computed symmetrically for ALL THREE contrasts
using this same per-family FDR procedure -- see Paper1_Additional_Analysis_Plan.md
Step 1. `primary_validated` remains a narrower, backward-compatible column
(True only for contraction_vs_baseline AND validated), since it is read by
name elsewhere (the 06-A reproducibility gate check, Stage 08's
load_validated_features()); new code needing fetal_movement_vs_baseline or
contraction_vs_fetal_movement validation status should read `validated`.

Multi-channel extension is intentionally NOT a modeling change: channel
is never a fixed effect or grouping factor in the formula above (per the
fixed protocol -- "No classification or detection models" and the
migration plan's explicit single-channel-per-run design). It is simply
run four times, once per channel, on that channel's own event-unit table.

Input:  results/paper1/04_feature_extraction/window_features_60s_{channel}.csv
        (read directly per channel; does not require Stage 05's long
        table, though it reads the same underlying Stage 04 output)
Output: results/paper1/06_mixed_models/
    confirmatory_feature_registry_60s_fwh.csv       (shared; feature-value-
                                                       independent, identical
                                                       across channels)
    event_unit_analysis_60s_{channel}.csv           (per-channel analysis units)
    event_unit_selection_60s_{channel}.csv          (per-channel selection audit)
    module2_5_preparation_summary_60s_{channel}.csv
    mixed_model_results_60s_fwh_{channel}.csv       (per-channel model results)
    mixed_model_fit_audit_60s_fwh_{channel}.csv
    mixed_model_summary_60s_fwh_{channel}.csv
    mixed_model_results_60s_fwh_all_channels.csv    (all 4 channels stacked)
"""
from __future__ import annotations

import argparse
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
_SRC_DIR = Path(__file__).resolve().parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(1, str(_PROJECT_ROOT))

import config as project_config  # noqa: E402


# =============================================================================
# Ported from pipeline/module2_5/prepare_analysis.py (unchanged)
# =============================================================================

EXCLUDED_SUFFIXES = {
    "time_sum_value": ("time_mean", "fixed window length: sum is a multiple of mean"),
    "time_integrated_absolute_value": ("time_mean_absolute_value", "fixed window length: IAV is a multiple of MAV"),
    "time_standard_deviation": ("time_root_mean_square", "near-zero-mean bandpass signal: SD is equivalent to RMS"),
    "time_variance": ("time_root_mean_square", "variance is SD squared and redundant with RMS"),
    "complexity_petrosian_fractal_dimension": ("time_slope_sign_change_count", "monotone transform of derivative sign changes at fixed N"),
    "time_mean": ("", "QC/descriptive only; near-zero after bandpass"),
    "time_minimum": ("", "QC/descriptive only; polarity/extreme-sensitive"),
    "time_maximum": ("", "QC/descriptive only; polarity/extreme-sensitive"),
}


def build_registry(features: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for band in ("fwh", "broad"):
        cols = [c for c in features.columns if c.startswith(f"{band}_")]
        for col in cols:
            suffix = col[len(band) + 1:]
            if suffix in EXCLUDED_SUFFIXES:
                rep, reason = EXCLUDED_SUFFIXES[suffix]
                status = "excluded"
            else:
                rep, reason, status = "", "eligible feature", "eligible"
            rows.append({"feature": col, "band": band, "suffix": suffix,
                         "status": status, "retained_representative": f"{band}_{rep}" if rep else "",
                         "reason": reason, "primary_confirmatory": bool(band == "fwh" and status == "eligible")})
    for col in [c for c in features.columns if c.startswith("freq_band_power_")]:
        rows.append({"feature": col, "band": "global", "suffix": col[15:], "status": "excluded",
                     "retained_representative": "", "reason": "absolute band power redundant with band-specific RMS",
                     "primary_confirmatory": False})
    return pd.DataFrame(rows)


def _select(group: pd.DataFrame, midpoint, registry: list) -> tuple:
    g = group.copy()
    g["_complete"] = g[registry].notna().all(axis=1)
    eligible = g[g["_complete"]]
    if eligible.empty:
        return None, "no_complete_registry_window"
    if midpoint is not None:
        eligible = eligible.assign(_distance=(eligible["window_start_s"] + eligible["window_end_s"]) / 2 - midpoint)
        eligible["_distance"] = eligible["_distance"].abs()
        best = eligible.sort_values(["_distance", "window_start_s", "window_id"], kind="mergesort").iloc[0]
    else:
        best = eligible.sort_values(["window_start_s", "window_id"], kind="mergesort").iloc[0]
    return best, "selected"


def build_event_units(features: pd.DataFrame, registry: list) -> tuple:
    selected, audit = [], []
    events = features[features["kind"].eq("event") & features["window_class"].isin(["contraction", "fetal_movement"])]
    for (woman, recording, ann, cls), g in events.groupby(
        ["patient_id", "recording_id", "annotation_index", "window_class"], dropna=False, sort=True
    ):
        midpoint = float(g["event_time_s"].dropna().iloc[0]) if g["event_time_s"].notna().any() else None
        row, reason = _select(g, midpoint, registry)
        unit = f"event:{recording}:{int(ann)}"
        for _, r in g.iterrows():
            audit.append({
                "analysis_unit": unit, "unit_type": "event", "patient_id": woman, "recording_id": recording,
                "annotation_index": ann, "class": cls, "window_id": r.window_id,
                "selection_status": "selected" if row is not None and r.window_id == row.window_id else "excluded",
                "exclusion_reason": "" if row is not None and r.window_id == row.window_id else reason,
            })
        if row is not None:
            selected.append(row)
    base = features[features["kind"].eq("baseline")].copy()
    for (woman, recording), g in base.groupby(["patient_id", "recording_id"], sort=True):
        g = g.sort_values(["window_start_s", "window_id"], kind="mergesort")
        for i, (_, r) in enumerate(g.iterrows(), 1):
            unit = f"baseline:{recording}:{i:04d}"
            audit.append({
                "analysis_unit": unit, "unit_type": "baseline", "patient_id": woman, "recording_id": recording,
                "annotation_index": -1, "class": "baseline", "window_id": r.window_id,
                "selection_status": "selected", "exclusion_reason": "",
            })
            selected.append(r)
    out = pd.DataFrame(selected).reset_index(drop=True)
    aud = pd.DataFrame(audit)
    if not out.empty:
        out.insert(0, "analysis_unit", aud[aud.selection_status.eq("selected")].analysis_unit.to_numpy())
        out.insert(1, "unit_type", aud[aud.selection_status.eq("selected")].unit_type.to_numpy())
    return out, aud


# =============================================================================
# Ported from pipeline/module2_5/run_models.py (unchanged)
# =============================================================================

STRUCTURE_TIERS = (("woman_recording", True), ("woman_only", False))
OPTIMIZER_METHODS = ("lbfgs", "bfgs")
VARIANCE_BOUNDARY_THRESHOLD = 1e-8
MIN_BASELINE_N = 5
MIN_CONTRACTION_N = 5


@dataclass
class FitAttempt:
    fit: Any
    structure: str
    method: str
    messages: list = field(default_factory=list)
    boundary_warning: bool = False
    iterations: float = float("nan")
    optimizer_message: str = None
    optimizer_status: float = float("nan")


def _safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _extract_optimizer_metadata(fit: Any) -> tuple:
    try:
        for candidate in (
            getattr(fit, "mle_retvals", None),
            getattr(getattr(fit, "_results", None), "mle_retvals", None),
        ):
            if isinstance(candidate, dict):
                iterations = _safe_float(candidate.get("iterations", candidate.get("nit")))
                status = _safe_float(candidate.get("warnflag", candidate.get("status")))
                message = candidate.get("message")
                return iterations, (str(message) if message is not None else None), status
    except Exception:
        pass
    return float("nan"), None, float("nan")


def _attempt_fit(formula: str, d: pd.DataFrame, method: str, use_vc: bool, structure: str) -> FitAttempt:
    kwargs: dict = {"groups": d["patient_id"], "re_formula": "1"}
    if use_vc:
        kwargs["vc_formula"] = {"recording": "0 + C(recording_id)"}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = smf.mixedlm(formula, d, **kwargs).fit(reml=True, method=method, maxiter=1000, disp=False)
        messages = [str(w.message) for w in caught]
    if not bool(getattr(fit, "converged", False)):
        raise RuntimeError(f"{method} optimizer did not report convergence")
    boundary_warning = any(
        any(token in m.lower() for token in ("boundary", "random effects covariance is singular", "singular covariance"))
        for m in messages
    )
    iterations, optimizer_message, optimizer_status = _extract_optimizer_metadata(fit)
    return FitAttempt(fit=fit, structure=structure, method=method, messages=messages,
                       boundary_warning=boundary_warning, iterations=iterations,
                       optimizer_message=optimizer_message, optimizer_status=optimizer_status)


def _try_tier(formula: str, d: pd.DataFrame, use_vc: bool, structure: str, all_warnings: list) -> FitAttempt:
    first_boundary = None
    for method in OPTIMIZER_METHODS:
        try:
            attempt = _attempt_fit(formula, d, method, use_vc, structure)
        except Exception as exc:
            all_warnings.append(f"{structure}/{method}: {exc}")
            continue
        all_warnings.extend(f"{structure}/{method}: {m}" for m in attempt.messages)
        if not attempt.boundary_warning:
            return attempt
        if first_boundary is None:
            first_boundary = attempt
    return first_boundary


def _validate_class_counts(n_baseline: int, n_contraction: int):
    baseline_short = n_baseline < MIN_BASELINE_N
    contraction_short = n_contraction < MIN_CONTRACTION_N
    if baseline_short and contraction_short:
        return (f"insufficient baseline (n={n_baseline} < {MIN_BASELINE_N}) and "
                f"contraction (n={n_contraction} < {MIN_CONTRACTION_N}) observations")
    if baseline_short:
        return f"insufficient baseline observations (n={n_baseline} < {MIN_BASELINE_N})"
    if contraction_short:
        return f"insufficient contraction observations (n={n_contraction} < {MIN_CONTRACTION_N})"
    return None


def _fit_one(data: pd.DataFrame, feature: str) -> tuple:
    d = data[[feature, "window_class", "patient_id", "recording_id"]].copy()
    d = d.replace([np.inf, -np.inf], np.nan).dropna()
    d["window_class"] = pd.Categorical(d["window_class"], categories=["baseline", "fetal_movement", "contraction"])
    d = d.dropna(subset=["window_class"])

    base: dict = {
        "feature": feature, "n": len(d), "n_women": d.patient_id.nunique(), "n_recordings": d.recording_id.nunique(),
        "n_baseline": int((d.window_class == "baseline").sum()),
        "n_contraction": int((d.window_class == "contraction").sum()),
        "n_fetal_movement": int((d.window_class == "fetal_movement").sum()),
        "model": "", "optimizer": "", "iterations": float("nan"), "optimizer_message": None,
        "optimizer_status": float("nan"), "converged": False, "fallback": False,
        "boundary_fit": False, "boundary_reason": "", "standardized": True,
        "feature_mean_raw": float("nan"), "feature_sd_raw": float("nan"),
        "woman_variance": float("nan"), "recording_variance": float("nan"), "residual_variance": float("nan"),
        "warning": "", "warning_summary": "", "warning_detail": "",
    }

    skip_reason = _validate_class_counts(base["n_baseline"], base["n_contraction"])
    if skip_reason is not None:
        base["warning"] = base["warning_detail"] = skip_reason
        base["warning_summary"] = "insufficient_data"
        return [], base

    raw_mean = float(d[feature].mean())
    raw_sd = float(d[feature].std())
    if not np.isfinite(raw_sd) or raw_sd == 0:
        msg = "zero or non-finite variance in feature; cannot standardize/fit"
        base["warning"] = base["warning_detail"] = msg
        base["warning_summary"] = "zero_variance"
        return [], base
    base["feature_mean_raw"], base["feature_sd_raw"] = raw_mean, raw_sd
    d["_z"] = (d[feature] - raw_mean) / raw_sd

    formula = '_z ~ C(window_class, Treatment(reference="baseline"))'
    all_warnings: list = []

    nested = _try_tier(formula, d, True, "woman_recording", all_warnings)
    if nested is not None:
        attempt, fallback = nested, False
    else:
        woman_only = _try_tier(formula, d, False, "woman_only", all_warnings)
        if woman_only is None:
            msg = " | ".join(all_warnings)
            base["warning"] = base["warning_detail"] = msg
            base["warning_summary"] = "optimizer_failure"
            return [], base
        attempt, fallback = woman_only, True

    fit = attempt.fit
    base.update({
        "model": attempt.structure, "optimizer": attempt.method, "iterations": attempt.iterations,
        "optimizer_message": attempt.optimizer_message, "optimizer_status": attempt.optimizer_status,
        "converged": True, "fallback": fallback,
    })

    try:
        base["woman_variance"] = float(fit.cov_re.iloc[0, 0])
    except Exception:
        pass
    if attempt.structure == "woman_recording":
        try:
            base["recording_variance"] = float(fit.vcomp[0])
        except Exception:
            pass
    try:
        base["residual_variance"] = float(fit.scale)
    except Exception:
        pass

    boundary_reasons = []
    if attempt.boundary_warning:
        warning_lower = " ".join(attempt.messages).lower()
        if "random effects covariance is singular" in warning_lower or "singular covariance" in warning_lower:
            boundary_reasons.append("singular_random_effects_warning")
        if "boundary" in warning_lower:
            boundary_reasons.append("optimizer_boundary_warning")
    if np.isfinite(base["woman_variance"]) and base["woman_variance"] < VARIANCE_BOUNDARY_THRESHOLD:
        boundary_reasons.append("near_zero_woman_variance")
    if (attempt.structure == "woman_recording" and np.isfinite(base["recording_variance"])
            and base["recording_variance"] < VARIANCE_BOUNDARY_THRESHOLD):
        boundary_reasons.append("near_zero_recording_variance")
    base["boundary_fit"] = bool(boundary_reasons)
    base["boundary_reason"] = ",".join(boundary_reasons)

    warning_text = " | ".join(all_warnings)
    base["warning"] = base["warning_detail"] = warning_text
    if attempt.structure == "woman_recording":
        base["warning_summary"] = "nested_boundary" if base["boundary_fit"] else "nested_clean"
    else:
        base["warning_summary"] = "woman_only_boundary" if base["boundary_fit"] else "woman_only_clean"

    rows = []
    for contrast, key in [
        ("fetal_movement_vs_baseline", 'C(window_class, Treatment(reference="baseline"))[T.fetal_movement]'),
        ("contraction_vs_baseline", 'C(window_class, Treatment(reference="baseline"))[T.contraction]'),
    ]:
        if key not in fit.params:
            continue
        ci = fit.conf_int().loc[key]
        estimate_z = float(fit.params[key])
        ci_low_z, ci_high_z = float(ci.iloc[0]), float(ci.iloc[1])
        rows.append({
            **base, "contrast": contrast, "term": key,
            "estimate": estimate_z, "estimate_raw_units": estimate_z * raw_sd,
            "std_error": float(fit.bse[key]), "statistic": float(fit.tvalues[key]), "p_value": float(fit.pvalues[key]),
            "ci_low": ci_low_z, "ci_high": ci_high_z,
            "ci_low_raw": ci_low_z * raw_sd, "ci_high_raw": ci_high_z * raw_sd,
        })

    c_key = 'C(window_class, Treatment(reference="baseline"))[T.contraction]'
    fm_key = 'C(window_class, Treatment(reference="baseline"))[T.fetal_movement]'
    if c_key in fit.params and fm_key in fit.params:
        estimate_z = float(fit.params[c_key] - fit.params[fm_key])
        covariance = float(fit.cov_params().loc[c_key, fm_key])
        variance = float(fit.cov_params().loc[c_key, c_key] + fit.cov_params().loc[fm_key, fm_key] - 2 * covariance)
        if np.isfinite(variance) and variance > 0:
            se = float(np.sqrt(variance))
            statistic = estimate_z / se
            from scipy.stats import norm
            p_value = float(2 * norm.sf(abs(statistic)))
            ci_low_z, ci_high_z = estimate_z - 1.96 * se, estimate_z + 1.96 * se
            rows.append({
                **base, "contrast": "contraction_vs_fetal_movement", "term": f"{c_key} - {fm_key}",
                "estimate": estimate_z, "estimate_raw_units": estimate_z * raw_sd,
                "std_error": se, "statistic": statistic, "p_value": p_value,
                "ci_low": ci_low_z, "ci_high": ci_high_z,
                "ci_low_raw": ci_low_z * raw_sd, "ci_high_raw": ci_high_z * raw_sd,
            })
    return rows, base


def _apply_fdr(result_df: pd.DataFrame) -> pd.DataFrame:
    """Benjamini-Hochberg FDR, applied separately per contrast family (three
    independent families: contraction_vs_baseline, fetal_movement_vs_baseline,
    contraction_vs_fetal_movement -- never pooled into one flat correction).

    `validated` (q < 0.05) is now computed symmetrically for ALL THREE
    contrasts, using the identical per-family BH-FDR procedure -- this is
    Paper1_Additional_Analysis_Plan.md Step 1, resolved as the "extend flat
    FDR symmetrically" option rather than adding a domain-level omnibus
    tier: the existing per-contrast-family flat FDR (already validated
    against the frozen EHG9 reference) is retained unchanged, just applied
    to all three contrasts instead of one.

    `primary_validated` is kept as a separate, narrower column (True only
    for contraction_vs_baseline AND validated) for backward compatibility:
    it is read by name in validate_reproducibility.py's 06-A gate check and
    by 08_subject_level_analysis.py's load_validated_features(), and
    changing its meaning would silently break both. New code that needs
    fetal_movement_vs_baseline or contraction_vs_fetal_movement validation
    status should read `validated`, not `primary_validated`.
    """
    result_df = result_df.copy()
    result_df["q_value"] = np.nan
    for contrast in ("contraction_vs_baseline", "fetal_movement_vs_baseline", "contraction_vs_fetal_movement"):
        mask = result_df["contrast"].eq(contrast) & result_df["p_value"].notna()
        if mask.any():
            result_df.loc[mask, "q_value"] = multipletests(result_df.loc[mask, "p_value"], method="fdr_bh")[1]
    result_df["validated"] = result_df["q_value"] < 0.05
    result_df["primary_validated"] = result_df["contrast"].eq("contraction_vs_baseline") & result_df["validated"]
    result_df["fit_quality"] = np.select(
        [
            (result_df["model"] == "woman_recording") & (~result_df["boundary_fit"]),
            (result_df["model"] == "woman_recording") & (result_df["boundary_fit"]),
            result_df["model"] == "woman_only",
        ],
        ["nested_clean", "nested_boundary", "woman_only_fallback"],
        default="unknown",
    )
    return result_df


def run_channel(channel: str, features_dir: Path, results_dir: Path, registry_features: list) -> dict:
    """Run Stage 06 for a single channel; write per-channel outputs; return
    a summary dict for the cross-channel roll-up."""
    path = features_dir / f"window_features_60s_{channel}.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Stage 04 output missing for {channel}: {path}")
    df = pd.read_csv(path)

    units, audit = build_event_units(df, registry_features)
    units.to_csv(results_dir / f"event_unit_analysis_60s_{channel}.csv", index=False)
    audit.to_csv(results_dir / f"event_unit_selection_60s_{channel}.csv", index=False)
    prep_summary = pd.DataFrame([{
        "channel": channel, "input_rows": len(df),
        "primary_fwh_features": len(registry_features), "analysis_units": len(units),
        "event_units": int((units.unit_type == "event").sum()) if not units.empty else 0,
        "baseline_units": int((units.unit_type == "baseline").sum()) if not units.empty else 0,
        "excluded_windows": int((audit.selection_status == "excluded").sum()),
    }])
    prep_summary.to_csv(results_dir / f"module2_5_preparation_summary_60s_{channel}.csv", index=False)

    results, audits = [], []
    for feature in registry_features:
        rows, feat_audit = _fit_one(units, feature)
        results.extend(rows)
        audits.append(feat_audit)

    result_df = pd.DataFrame(results)
    if not result_df.empty:
        result_df = _apply_fdr(result_df)
        result_df.insert(0, "channel", channel)
    result_df.to_csv(results_dir / f"mixed_model_results_60s_fwh_{channel}.csv", index=False)

    audit_df = pd.DataFrame(audits)
    audit_df.insert(0, "channel", channel)
    audit_df.to_csv(results_dir / f"mixed_model_fit_audit_60s_fwh_{channel}.csv", index=False)

    nested_clean = int(((audit_df.model == "woman_recording") & (~audit_df.boundary_fit)).sum())
    nested_boundary = int(((audit_df.model == "woman_recording") & (audit_df.boundary_fit)).sum())
    woman_only_fallback = int((audit_df.model == "woman_only").sum())
    failed = int((~audit_df.converged).sum())
    def _validated_count(contrast: str) -> int:
        if result_df.empty:
            return 0
        return int(result_df.loc[result_df.contrast.eq(contrast), "validated"].sum())

    primary_validated = _validated_count("contraction_vs_baseline")
    fm_validated = _validated_count("fetal_movement_vs_baseline")
    ucfm_validated = _validated_count("contraction_vs_fetal_movement")
    summary_row = {
        "channel": channel, "features_requested": len(registry_features),
        "features_with_results": int(result_df.feature.nunique()) if not result_df.empty else 0,
        "nested_clean": nested_clean, "nested_boundary": nested_boundary,
        "woman_only_fallback": woman_only_fallback, "failed": failed,
        "primary_validated_features": primary_validated,
        "contraction_vs_baseline_validated": primary_validated,
        "fetal_movement_vs_baseline_validated": fm_validated,
        "contraction_vs_fetal_movement_validated": ucfm_validated,
        "n_event_units": int((units.unit_type == "event").sum()) if not units.empty else 0,
        "n_baseline_units": int((units.unit_type == "baseline").sum()) if not units.empty else 0,
    }
    pd.DataFrame([summary_row]).to_csv(results_dir / f"mixed_model_summary_60s_fwh_{channel}.csv", index=False)
    return summary_row, result_df


def run(
    channels: tuple = project_config.CHANNELS,
    features_dir: Path = None,
    results_dir: Path = None,
) -> int:
    features_dir = features_dir or (project_config.FEATURES_DIR)
    results_dir = results_dir or (project_config.UNITS_DIR)
    results_dir.mkdir(parents=True, exist_ok=True)

    # Registry is feature-value-independent (depends only on which columns
    # exist), so build it once from any one channel's file and reuse for all.
    first_df = pd.read_csv(features_dir / f"window_features_60s_{channels[0]}.csv")
    registry = build_registry(first_df)
    registry.to_csv(results_dir / "confirmatory_feature_registry_60s_fwh.csv", index=False)
    registry_features = registry.loc[registry.primary_confirmatory, "feature"].tolist()

    summaries, all_results = [], []
    for channel in channels:
        print(f"=== channel={channel} ===")
        summary_row, result_df = run_channel(channel, features_dir, results_dir, registry_features)
        summaries.append(summary_row)
        if not result_df.empty:
            all_results.append(result_df)
        print(pd.DataFrame([summary_row]).to_string(index=False))

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(results_dir / "mixed_model_results_60s_fwh_all_channels.csv", index=False)

    pd.DataFrame(summaries).to_csv(results_dir / "mixed_model_summary_60s_fwh_all_channels.csv", index=False)
    # Publish compact manuscript tables; keep event-unit and fit intermediates in work/.
    from shutil import copy2
    final_dir = project_config.PRIMARY_DIR
    final_dir.mkdir(parents=True, exist_ok=True)
    selected = ["confirmatory_feature_registry_60s_fwh.csv",
                "mixed_model_results_60s_fwh_all_channels.csv",
                "mixed_model_fit_audit_60s_fwh_EHG9.csv"]
    selected += [f"mixed_model_results_60s_fwh_{channel}.csv" for channel in channels]
    for name in selected:
        path = results_dir / name
        if path.exists():
            copy2(path, final_dir / name)
    print(f"done. outputs -> {results_dir}")
    return 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--channels", type=str, default=None, help="Comma-separated, e.g. EHG9,EHG10")
    parser.add_argument("--features-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    channels = tuple(c.strip() for c in args.channels.split(",")) if args.channels else project_config.CHANNELS
    sys.exit(run(channels=channels, features_dir=args.features_dir, results_dir=args.results_dir))
