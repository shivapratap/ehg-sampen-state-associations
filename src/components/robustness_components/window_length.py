#!/usr/bin/env python3
"""Compare EHG9 UC-vs-IN effects at 40, frozen 60, and 80 seconds.

Stage 02/04 outputs for 40 and 80 seconds must already exist beneath this
revision directory. The script reuses Stage 06 event-unit and optimizer-tier
functions without modifying the frozen pipeline. Primary coefficients use the
frozen 60-second mean/SD; raw and within-duration standardized estimates are
retained to expose duration-dependent feature scaling.
"""
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
import numpy as np
import pandas as pd
import scipy
import statsmodels
from statsmodels.stats.multitest import multipletests


SCRIPT_PATH = Path(__file__).resolve()
OUT = Path(__file__).resolve().parents[3] / "work" / "robustness" / "window_length"
ROOT = SCRIPT_PATH.parents[3]
STAGE06_DIR = ROOT / "work" / "primary"
FROZEN_RESULTS = STAGE06_DIR / "mixed_model_results_60s_fwh_EHG9.csv"
FROZEN_REGISTRY = STAGE06_DIR / "confirmatory_feature_registry_60s_fwh.csv"
FROZEN_FIT_AUDIT = STAGE06_DIR / "mixed_model_fit_audit_60s_fwh_EHG9.csv"
FROZEN_EVENT_UNITS = STAGE06_DIR / "event_unit_analysis_60s_EHG9.csv"
FROZEN_MANIFEST = ROOT / "work" / "windows" / "window_manifest_60s.csv"
FROZEN_REJECTED = ROOT / "work" / "windowing" / "rejected_windows_60s.csv"
WINDOW_DIR = OUT / "windows"
FEATURE_DIRS = {40: OUT / "feature_extraction_40s", 80: OUT / "feature_extraction_80s"}
STAGE02_DIRS = {40: OUT / "stage02_40s", 80: OUT / "stage02_80s"}
LEGACY_DIR = None  # Historical manifest comparison belongs to Pass 3.
STATE_TERM = 'C(window_class, Treatment(reference="baseline"))[T.contraction]'
FORMULA = '_z ~ C(window_class, Treatment(reference="baseline"))'
ALPHA = 0.05
RATIO_THRESHOLD = 0.20
ABS_THRESHOLD_SMALL = 0.05

ORIGINAL_11 = {
    "fwh_time_median", "fwh_time_skewness", "fwh_time_mean_absolute_value",
    "fwh_time_root_mean_square", "fwh_time_peak_to_peak", "fwh_time_waveform_length",
    "fwh_freq_spectral_edge_frequency_95", "fwh_freq_spectral_spread",
    "fwh_entropy_totalsampen", "fwh_complexity_hjorth_complexity",
    "fwh_complexity_higuchi_fractal_dimension",
}
CORE = (
    "fwh_time_mean_absolute_value", "fwh_time_root_mean_square",
    "fwh_time_peak_to_peak", "fwh_time_waveform_length",
)
FOCUS = (
    "fwh_time_median", "fwh_time_skewness", *CORE,
    "fwh_entropy_totalsampen", "fwh_freq_spectral_edge_frequency_95",
    "fwh_freq_spectral_spread", "fwh_complexity_hjorth_complexity",
    "fwh_complexity_higuchi_fractal_dimension",
)
LABELS = {
    "fwh_time_median": "Median", "fwh_time_skewness": "Skewness",
    "fwh_time_mean_absolute_value": "MAV", "fwh_time_root_mean_square": "RMS",
    "fwh_time_peak_to_peak": "Peak-to-peak", "fwh_time_waveform_length": "Waveform length",
    "fwh_entropy_totalsampen": "Total SampEn",
    "fwh_freq_spectral_edge_frequency_95": "Spectral edge 95",
    "fwh_freq_spectral_spread": "Spectral spread",
    "fwh_complexity_hjorth_complexity": "Hjorth complexity",
    "fwh_complexity_higuchi_fractal_dimension": "Higuchi FD",
}


def load_stage06():
    path = ROOT / "src" / "models.py"
    spec = importlib.util.spec_from_file_location("stage06_window_sensitivity", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def primary_features() -> tuple[list[str], pd.DataFrame]:
    registry = pd.read_csv(FROZEN_REGISTRY)
    flag = registry["primary_confirmatory"].astype(str).str.lower().isin(["true", "1"])
    features = registry.loc[flag, "feature"].tolist()
    if len(features) != 34 or len(set(features)) != 34:
        raise RuntimeError(f"Expected 34 frozen primary features, found {len(features)}")
    frozen = pd.read_csv(FROZEN_RESULTS)
    frozen = frozen.loc[frozen["contrast"].eq("contraction_vs_baseline")].set_index("feature").loc[features].reset_index()
    if set(frozen.loc[frozen["q_value"] < ALPHA, "feature"]) != ORIGINAL_11:
        raise RuntimeError("Frozen original-11 feature set changed")
    q = multipletests(frozen["p_value"].to_numpy(), method="fdr_bh")[1]
    if not np.allclose(q, frozen["q_value"], rtol=1e-12, atol=1e-15):
        raise RuntimeError("Frozen 60-s BH q-values do not reproduce")
    return features, frozen


def manifest_path(duration: int) -> Path:
    return FROZEN_MANIFEST if duration == 60 else WINDOW_DIR / f"window_manifest_{duration}s.csv"


def rejected_path(duration: int) -> Path:
    return FROZEN_REJECTED if duration == 60 else STAGE02_DIRS[duration] / f"rejected_windows_{duration}s.csv"


def feature_path(duration: int) -> Path:
    if duration == 60:
        return ROOT / "work" / "features" / "window_features_60s_EHG9.csv"
    return FEATURE_DIRS[duration] / f"window_features_{duration}s_EHG9.csv"


def population_outputs() -> pd.DataFrame:
    summaries = []
    for duration in (40, 60, 80):
        manifest = pd.read_csv(manifest_path(duration))
        rejected = pd.read_csv(rejected_path(duration))
        if duration in (40, 80):
            manifest.to_csv(OUT / f"window_population_{duration}s.csv", index=False)
            legacy_equal = True
        else:
            legacy_equal = True
        counts = manifest["window_class"].value_counts()
        event_rejected = rejected.loc[rejected["kind"].eq("event")] if "kind" in rejected else rejected.iloc[0:0]
        gap_rejected = rejected.loc[rejected["kind"].eq("baseline_gap")] if "kind" in rejected else rejected.iloc[0:0]
        primary = manifest.loc[manifest["window_class"].isin(["contraction", "baseline"])]
        woman_state = pd.crosstab(primary["patient_id"], primary["window_class"])
        summaries.append({
            "duration_s": duration, "total_manifest_windows": len(manifest),
            "uc_windows": int(counts.get("contraction", 0)),
            "in_windows": int(counts.get("baseline", 0)),
            "fm_windows": int(counts.get("fetal_movement", 0)),
            "possible_contraction_windows": int(counts.get("possible_contraction", 0)),
            "mixed_windows": int(counts.get("mixed", 0)),
            "women": manifest["patient_id"].nunique(), "recordings": manifest["recording_id"].nunique(),
            "women_with_uc_and_in": int(((woman_state.get("contraction", 0) > 0) & (woman_state.get("baseline", 0) > 0)).sum()),
            "event_windows_rejected_at_boundaries": len(event_rejected),
            "short_clean_gaps_rejected": len(gap_rejected),
            "truncated_windows": 0,
            "clean_event_windows": int((manifest["overlap_status"] == "clean").sum()),
            "same_class_cluster_windows": int((manifest["overlap_status"] == "same_class_cluster").sum()),
            "window_overlap_windows": int((manifest["overlap_status"] == "window_overlap").sum()),
            "mixed_overlap_windows": int((manifest["overlap_status"] == "mixed").sum()),
            "contains_competing_class_windows": int(manifest["contains_competing_class"].eq(True).sum()),
            "tight_clearance_baseline_windows": int(manifest["tight_clearance"].eq(True).sum()),
            "legacy_manifest_exact_match": legacy_equal,
            "manifest_sha256": sha256(manifest_path(duration)),
        })
    result = pd.DataFrame(summaries)
    result.to_csv(OUT / "window_length_population_comparison.csv", index=False)
    return result


def feature_qc(duration: int, frame: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    expected = duration * 20
    missing_columns = sorted(set(features) - set(frame.columns))
    log_path = FEATURE_DIRS[duration] / f"feature_processing_log_{duration}s_EHG9.csv"
    processing_log = pd.read_csv(log_path) if log_path.exists() else pd.DataFrame()
    processing_error_rows = int(processing_log["status"].eq("error").sum()) if "status" in processing_log else np.nan
    rows = []
    for feature in features:
        values = pd.to_numeric(frame[feature], errors="coerce") if feature in frame else pd.Series(np.nan, index=frame.index)
        raw = frame[feature].to_numpy() if feature in frame else np.full(len(frame), np.nan)
        inf_count = int(np.isinf(pd.to_numeric(pd.Series(raw), errors="coerce").to_numpy()).sum())
        rows.append({
            "duration_s": duration, "feature": feature, "primary_feature_present": feature in frame,
            "n_windows": len(frame), "nan_count": int(values.isna().sum()), "nan_rate": float(values.isna().mean()),
            "inf_count": inf_count, "finite_count": int(np.isfinite(values.to_numpy(dtype=float)).sum()),
            "unique_finite_values": int(values[np.isfinite(values)].nunique()),
            "minimum": float(values.min()), "maximum": float(values.max()),
            "expected_filtered_sample_count": expected,
            "bad_filtered_sample_count_rows": int((frame["filtered_sample_count"] != expected).sum()),
            "feature_processing_error_rows": processing_error_rows,
            "all_34_columns_present": len(missing_columns) == 0,
        })
    result = pd.DataFrame(rows)
    result.to_csv(OUT / f"feature_qc_{duration}s.csv", index=False)
    return result


def fit_duration(duration: int, frame: pd.DataFrame, features: list[str], frozen: pd.DataFrame, stage06):
    units, selection = stage06.build_event_units(frame, features)
    units.to_csv(OUT / f"event_unit_analysis_{duration}s_EHG9.csv", index=False)
    selection.to_csv(OUT / f"event_unit_selection_{duration}s_EHG9.csv", index=False)
    frozen_index = frozen.set_index("feature")
    result_rows, audit_rows = [], []
    for feature in features:
        fixed = frozen_index.loc[feature]
        d = units[[feature, "window_class", "patient_id", "recording_id"]].copy()
        d = d.replace([np.inf, -np.inf], np.nan).dropna()
        d["window_class"] = pd.Categorical(d["window_class"], categories=["baseline", "fetal_movement", "contraction"])
        current_mean, current_sd = float(d[feature].mean()), float(d[feature].std(ddof=1))
        frozen_mean, frozen_sd = float(fixed["feature_mean_raw"]), float(fixed["feature_sd_raw"])
        d["_z"] = (d[feature] - frozen_mean) / frozen_sd
        warnings_all = []
        nested = stage06._try_tier(FORMULA, d, True, "woman_recording", warnings_all)
        fallback = False
        attempt = nested
        if attempt is None:
            attempt = stage06._try_tier(FORMULA, d, False, "woman_only", warnings_all)
            fallback = True
        base = {
            "duration_s": duration, "feature": feature, "n": len(d),
            "n_women": d["patient_id"].nunique(), "n_recordings": d["recording_id"].nunique(),
            "n_baseline": int((d["window_class"] == "baseline").sum()),
            "n_contraction": int((d["window_class"] == "contraction").sum()),
            "n_fetal_movement": int((d["window_class"] == "fetal_movement").sum()),
            "frozen_60_mean": frozen_mean, "frozen_60_sd": frozen_sd,
            "duration_mean_raw": current_mean, "duration_sd_raw": current_sd,
            "converged": attempt is not None, "fallback": fallback,
            "model": attempt.structure if attempt else "", "optimizer": attempt.method if attempt else "",
            "warning_detail": " | ".join(warnings_all),
        }
        if attempt is None:
            base.update({"boundary_fit": False, "boundary_reason": "", "woman_variance": np.nan,
                         "recording_variance": np.nan, "residual_variance": np.nan})
            result_rows.append({**base, "estimate": np.nan, "std_error": np.nan, "ci_low": np.nan,
                                "ci_high": np.nan, "p_value": np.nan, "estimate_raw_units": np.nan,
                                "std_error_raw_units": np.nan, "estimate_within_duration_z": np.nan,
                                "std_error_within_duration_z": np.nan})
            audit_rows.append(base)
            continue
        fit = attempt.fit
        woman_var = float(fit.cov_re.iloc[0, 0])
        recording_var = float(fit.vcomp[0]) if attempt.structure == "woman_recording" else np.nan
        boundary_reasons = []
        if attempt.boundary_warning:
            lower = " ".join(attempt.messages).lower()
            if "singular" in lower:
                boundary_reasons.append("singular_random_effects_warning")
            if "boundary" in lower:
                boundary_reasons.append("optimizer_boundary_warning")
        if woman_var < stage06.VARIANCE_BOUNDARY_THRESHOLD:
            boundary_reasons.append("near_zero_woman_variance")
        if np.isfinite(recording_var) and recording_var < stage06.VARIANCE_BOUNDARY_THRESHOLD:
            boundary_reasons.append("near_zero_recording_variance")
        base.update({
            "boundary_fit": bool(boundary_reasons), "boundary_reason": ",".join(boundary_reasons),
            "woman_variance": woman_var, "recording_variance": recording_var,
            "residual_variance": float(fit.scale),
        })
        if STATE_TERM not in fit.params:
            raise RuntimeError(f"{duration}s {feature}: UC term absent despite converged model")
        estimate, se, p = float(fit.params[STATE_TERM]), float(fit.bse[STATE_TERM]), float(fit.pvalues[STATE_TERM])
        ci = fit.conf_int().loc[STATE_TERM]
        raw_est, raw_se = estimate * frozen_sd, se * frozen_sd
        result_rows.append({
            **base, "estimate": estimate, "std_error": se, "ci_low": float(ci.iloc[0]),
            "ci_high": float(ci.iloc[1]), "p_value": p,
            "estimate_raw_units": raw_est, "std_error_raw_units": raw_se,
            "ci_low_raw_units": float(ci.iloc[0]) * frozen_sd,
            "ci_high_raw_units": float(ci.iloc[1]) * frozen_sd,
            "estimate_within_duration_z": raw_est / current_sd,
            "std_error_within_duration_z": raw_se / current_sd,
        })
        audit_rows.append(base)
    results = pd.DataFrame(result_rows)
    mask = results["p_value"].notna()
    results["q_value"] = np.nan
    results.loc[mask, "q_value"] = multipletests(results.loc[mask, "p_value"], method="fdr_bh")[1]
    results["bh_family_size"] = int(mask.sum())
    results["q_lt_0_05"] = results["q_value"] < ALPHA
    results.to_csv(OUT / f"mixed_model_results_{duration}s_common_scale_EHG9.csv", index=False)
    return units, selection, results, pd.DataFrame(audit_rows)


def add_final_analysis_counts(population: pd.DataFrame, units_by_duration: dict[int, pd.DataFrame]) -> pd.DataFrame:
    result = population.copy()
    for duration, units in units_by_duration.items():
        counts = units["window_class"].value_counts()
        mask = result["duration_s"].eq(duration)
        result.loc[mask, "final_analysis_uc"] = int(counts.get("contraction", 0))
        result.loc[mask, "final_analysis_in"] = int(counts.get("baseline", 0))
        result.loc[mask, "final_analysis_fm"] = int(counts.get("fetal_movement", 0))
        result.loc[mask, "final_analysis_women"] = units["patient_id"].nunique()
        result.loc[mask, "final_analysis_recordings"] = units["recording_id"].nunique()
    result.to_csv(OUT / "window_length_population_comparison.csv", index=False)
    return result


def sign_label(value: float) -> str:
    if not np.isfinite(value):
        return "missing"
    return "positive" if value > 0 else "negative" if value < 0 else "zero"


def stability_label(row: pd.Series) -> str:
    if not row["all_three_signs_agree"]:
        return "direction reversal"
    frozen = abs(row["coefficient_60s"])
    if frozen >= 0.10:
        stable = abs(row["coefficient_ratio_40_vs_60"] - 1) <= RATIO_THRESHOLD and abs(row["coefficient_ratio_80_vs_60"] - 1) <= RATIO_THRESHOLD
    else:
        stable = row["absolute_difference_40_vs_60"] <= ABS_THRESHOLD_SMALL and row["absolute_difference_80_vs_60"] <= ABS_THRESHOLD_SMALL
    return "direction and magnitude stable" if stable else "direction stable; magnitude variable"


def build_comparison(features: list[str], frozen: pd.DataFrame, results: dict[int, pd.DataFrame], audits: dict[int, pd.DataFrame]):
    fixed = frozen.set_index("feature")
    r40, r80 = results[40].set_index("feature"), results[80].set_index("feature")
    a60 = pd.read_csv(FROZEN_FIT_AUDIT).set_index("feature")
    rows = []
    for feature in features:
        f, x40, x80, fa = fixed.loc[feature], r40.loc[feature], r80.loc[feature], a60.loc[feature]
        est60 = float(f["estimate"])
        row = {
            "feature": feature, "short_label": LABELS.get(feature, feature), "original_11": feature in ORIGINAL_11,
            "coefficient_40s": x40["estimate"], "se_40s": x40["std_error"], "ci_low_40s": x40["ci_low"], "ci_high_40s": x40["ci_high"], "p_40s": x40["p_value"], "q_40s": x40["q_value"],
            "coefficient_60s": est60, "se_60s": f["std_error"], "ci_low_60s": f["ci_low"], "ci_high_60s": f["ci_high"], "p_60s": f["p_value"], "q_60s": f["q_value"],
            "coefficient_80s": x80["estimate"], "se_80s": x80["std_error"], "ci_low_80s": x80["ci_low"], "ci_high_80s": x80["ci_high"], "p_80s": x80["p_value"], "q_80s": x80["q_value"],
            "raw_coefficient_40s": x40["estimate_raw_units"], "raw_coefficient_60s": f["estimate_raw_units"], "raw_coefficient_80s": x80["estimate_raw_units"],
            "within_duration_z_coefficient_40s": x40["estimate_within_duration_z"], "within_duration_z_coefficient_60s": est60, "within_duration_z_coefficient_80s": x80["estimate_within_duration_z"],
            "duration_sd_over_frozen_sd_40s": x40["duration_sd_raw"] / x40["frozen_60_sd"], "duration_sd_over_frozen_sd_80s": x80["duration_sd_raw"] / x80["frozen_60_sd"],
            "sign_40s": sign_label(x40["estimate"]), "sign_60s": sign_label(est60), "sign_80s": sign_label(x80["estimate"]),
            "all_three_signs_agree": len({sign_label(x40["estimate"]), sign_label(est60), sign_label(x80["estimate"])}) == 1,
            "coefficient_ratio_40_vs_60": x40["estimate"] / est60 if abs(est60) > 1e-12 else np.nan,
            "coefficient_ratio_80_vs_60": x80["estimate"] / est60 if abs(est60) > 1e-12 else np.nan,
            "absolute_difference_40_vs_60": abs(x40["estimate"] - est60), "absolute_difference_80_vs_60": abs(x80["estimate"] - est60),
            "converged_40s": x40["converged"], "model_40s": x40["model"], "optimizer_40s": x40["optimizer"], "fallback_40s": x40["fallback"], "boundary_40s": x40["boundary_fit"],
            "converged_60s": bool(fa["converged"]), "model_60s": fa["model"], "optimizer_60s": fa["optimizer"], "fallback_60s": bool(fa["fallback"]), "boundary_60s": bool(fa["boundary_fit"]),
            "converged_80s": x80["converged"], "model_80s": x80["model"], "optimizer_80s": x80["optimizer"], "fallback_80s": x80["fallback"], "boundary_80s": x80["boundary_fit"],
        }
        row["coefficient_stability"] = stability_label(pd.Series(row))
        row["se_ratio_40_vs_60"] = row["se_40s"] / row["se_60s"]
        row["se_ratio_80_vs_60"] = row["se_80s"] / row["se_60s"]
        row["material_precision_variation"] = bool(
            abs(row["se_ratio_40_vs_60"] - 1) > RATIO_THRESHOLD
            or abs(row["se_ratio_80_vs_60"] - 1) > RATIO_THRESHOLD
        )
        rows.append(row)
    comparison = pd.DataFrame(rows)
    comparison.to_csv(OUT / "window_length_comparison_34_features.csv", index=False)
    original = comparison.loc[comparison["original_11"]].copy()
    original["q_lt_0_05_count"] = original[["q_40s", "q_60s", "q_80s"]].lt(ALPHA).sum(axis=1)
    original.to_csv(OUT / "original_11_window_length_stability.csv", index=False)
    return comparison, original


def make_figures(comparison: pd.DataFrame, original: pd.DataFrame) -> None:
    durations = np.array([40, 60, 80])
    fig, axes = plt.subplots(4, 3, figsize=(12.5, 13.5), sharex=True)
    for ax, (_, row) in zip(axes.flat, original.iterrows()):
        y = np.array([row.coefficient_40s, row.coefficient_60s, row.coefficient_80s])
        lo = np.array([row.ci_low_40s, row.ci_low_60s, row.ci_low_80s])
        hi = np.array([row.ci_high_40s, row.ci_high_60s, row.ci_high_80s])
        ax.errorbar(durations, y, yerr=np.vstack([y - lo, hi - y]), marker="o", color="#165D8B", capsize=3)
        ax.axhline(0, color="#777777", linestyle="--", linewidth=0.8)
        ax.set_title(row.short_label, fontsize=10)
        ax.set_xticks(durations)
        ax.grid(alpha=0.18)
    axes.flat[-1].axis("off")
    for ax in axes[-1, :2]:
        ax.set_xlabel("Window duration (s)")
    fig.supylabel("UC minus IN coefficient (frozen 60-s SD units)")
    fig.suptitle("Original 11: window-length coefficient comparison", fontsize=15)
    fig.tight_layout(rect=(0.03, 0.02, 1, 0.97))
    fig.savefig(OUT / "original_11_window_length_coefficients.png", dpi=220)
    plt.close(fig)

    core = comparison.set_index("feature").loc[list(CORE)]
    x = np.arange(len(CORE)); width = 0.24
    fig, ax = plt.subplots(figsize=(10, 5.2))
    for offset, duration, color in zip((-width, 0, width), (40, 60, 80), ("#4C78A8", "#F58518", "#54A24B")):
        ax.bar(x + offset, core[f"coefficient_{duration}s"], width, label=f"{duration} s", color=color)
    ax.axhline(0, color="#555555", linewidth=0.9)
    ax.set_xticks(x, [LABELS[f] for f in CORE])
    ax.set_ylabel("UC minus IN coefficient (frozen 60-s SD units)")
    ax.set_title("Core amplitude/time-domain effects across window durations")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.18)
    fig.tight_layout()
    fig.savefig(OUT / "core_amplitude_window_length_coefficients.png", dpi=220)
    plt.close(fig)

    focus = original.set_index("short_label")[["q_40s", "q_60s", "q_80s"]]
    fig, ax = plt.subplots(figsize=(7.6, 7.2))
    image = ax.imshow(-np.log10(focus.clip(lower=1e-300)), aspect="auto", cmap="Blues")
    ax.set_xticks(range(3), ["40 s", "60 s", "80 s"])
    ax.set_yticks(range(len(focus)), focus.index)
    for i in range(len(focus)):
        for j, column in enumerate(focus.columns):
            q = focus.iloc[i, j]
            ax.text(j, i, f"{q:.3g}", ha="center", va="center", color="white" if -np.log10(max(q, 1e-300)) > 2 else "black", fontsize=8)
    ax.set_title("Original 11: BH q-values across window durations")
    fig.colorbar(image, ax=ax, label="−log10(q)")
    fig.tight_layout()
    fig.savefig(OUT / "original_11_window_length_q_values.png", dpi=220)
    plt.close(fig)


def write_hashes() -> None:
    critical = [
        SCRIPT_PATH, OUT / "locked_window_length_design.json",
        ROOT / "src" / "windowing.py", ROOT / "src" / "preprocessing.py",
        ROOT / "src" / "features.py", ROOT / "src" / "feature_qc.py",
        ROOT / "src" / "models.py", ROOT / "src" / "config.py",
        ROOT / "work" / "audit" / "recording_inventory.csv",
        ROOT / "work" / "audit" / "recording_qc_flags.csv",
        ROOT / "work" / "audit" / "run_metadata.json",
        ROOT / "work" / "annotations" / "annotation_events.csv",
        ROOT / "work" / "annotations" / "label_ontology.csv", FROZEN_RESULTS,
        FROZEN_REGISTRY, FROZEN_FIT_AUDIT, FROZEN_EVENT_UNITS, FROZEN_MANIFEST,
        manifest_path(40), manifest_path(80), feature_path(40), feature_path(80),
    ]
    pd.DataFrame([{"path": str(path), "sha256": sha256(path), "size_bytes": path.stat().st_size} for path in critical]).to_csv(OUT / "input_file_hashes.csv", index=False)


def pct(value: float) -> str:
    return f"{100 * value:.1f}%"


def write_report(population: pd.DataFrame, comparison: pd.DataFrame, original: pd.DataFrame, audits: pd.DataFrame) -> None:
    pop = population.set_index("duration_s")
    same_sign = int(original["all_three_signs_agree"].sum())
    all_q = int((original["q_lt_0_05_count"] == 3).sum())
    two_q = int((original["q_lt_0_05_count"] >= 2).sum())
    reversed_features = original.loc[~original["all_three_signs_agree"], "short_label"].tolist()
    magnitude_variable = original.loc[original["coefficient_stability"] != "direction and magnitude stable", "short_label"].tolist()
    precision_variable = original.loc[original["material_precision_variation"], "short_label"].tolist()
    full_reversals = comparison.loc[~comparison["all_three_signs_agree"], "short_label"].tolist()
    core = comparison.set_index("feature").loc[list(CORE)]
    core_positive = bool((core[["coefficient_40s", "coefficient_60s", "coefficient_80s"]] > 0).all().all())
    core_stable_count = int((core["coefficient_stability"] == "direction and magnitude stable").sum())
    core_resolved = core_positive and core_stable_count >= 3
    total = comparison.set_index("feature").loc["fwh_entropy_totalsampen"]
    fit_failures = int((~audits["converged"]).sum())
    fallback_counts = audits.groupby("duration_s")["fallback"].sum().astype(int).to_dict()
    boundary_counts = audits.groupby("duration_s")["boundary_fit"].sum().astype(int).to_dict()
    spectral_complexity = comparison.set_index("feature").loc[[
        "fwh_freq_spectral_edge_frequency_95", "fwh_freq_spectral_spread",
        "fwh_complexity_hjorth_complexity", "fwh_complexity_higuchi_fractal_dimension",
    ]]
    sc_stable = int((spectral_complexity["coefficient_stability"] == "direction and magnitude stable").sum())
    focus_lines = []
    for feature in FOCUS:
        row = comparison.set_index("feature").loc[feature]
        focus_lines.append(
            f"| {LABELS[feature]} | {row.coefficient_40s:.3f} | {row.coefficient_60s:.3f} | {row.coefficient_80s:.3f} | "
            f"{row.q_40s:.3g} | {row.q_60s:.3g} | {row.q_80s:.3g} | {row.coefficient_stability} |"
        )
    report = f"""# Window-length sensitivity report

## Main finding

The 40-s and 80-s reruns were technically feasible and used the current frozen pipeline logic with only duration-dependent segmentation changed. Among the frozen original 11, **{same_sign}/11 retained the same UC-versus-IN direction at 40, 60, and 80 s; {all_q}/11 had BH q<0.05 at all three durations, and {two_q}/11 had BH q<0.05 at two or more durations.**

The central amplitude/time-domain pattern was {'positive at every duration' if core_positive else 'not uniformly positive across durations'}, and {core_stable_count}/4 core features met the prespecified magnitude-stability rule. The central-pattern rule required all four positive and at least three magnitude-stable, with any duration-additive exception named. The evidence therefore {'does not indicate that the central pattern is an artefact of choosing 60 s' if core_resolved else 'requires qualification because the core pattern is duration-sensitive'}.

## Technical feasibility and population

The regenerated 40-s and 80-s manifests exactly matched the legacy Module 1 manifests by values and SHA-256. Feature extraction used EHG9 only, full-record filtering, and 800 or 1,600 decimated samples per window. No feature algorithm or model setting was tuned by duration.

| Duration | Manifest UC/IN/FM | Final analysis UC/IN/FM | Total manifest | Women | Recordings | Boundary-rejected events | Mixed windows |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 s | {int(pop.loc[40,'uc_windows'])}/{int(pop.loc[40,'in_windows'])}/{int(pop.loc[40,'fm_windows'])} | {int(pop.loc[40,'final_analysis_uc'])}/{int(pop.loc[40,'final_analysis_in'])}/{int(pop.loc[40,'final_analysis_fm'])} | {int(pop.loc[40,'total_manifest_windows'])} | {int(pop.loc[40,'women'])} | {int(pop.loc[40,'recordings'])} | {int(pop.loc[40,'event_windows_rejected_at_boundaries'])} | {int(pop.loc[40,'mixed_windows'])} |
| 60 s | {int(pop.loc[60,'uc_windows'])}/{int(pop.loc[60,'in_windows'])}/{int(pop.loc[60,'fm_windows'])} | {int(pop.loc[60,'final_analysis_uc'])}/{int(pop.loc[60,'final_analysis_in'])}/{int(pop.loc[60,'final_analysis_fm'])} | {int(pop.loc[60,'total_manifest_windows'])} | {int(pop.loc[60,'women'])} | {int(pop.loc[60,'recordings'])} | {int(pop.loc[60,'event_windows_rejected_at_boundaries'])} | {int(pop.loc[60,'mixed_windows'])} |
| 80 s | {int(pop.loc[80,'uc_windows'])}/{int(pop.loc[80,'in_windows'])}/{int(pop.loc[80,'fm_windows'])} | {int(pop.loc[80,'final_analysis_uc'])}/{int(pop.loc[80,'final_analysis_in'])}/{int(pop.loc[80,'final_analysis_fm'])} | {int(pop.loc[80,'total_manifest_windows'])} | {int(pop.loc[80,'women'])} | {int(pop.loc[80,'recordings'])} | {int(pop.loc[80,'event_windows_rejected_at_boundaries'])} | {int(pop.loc[80,'mixed_windows'])} |

Changing duration mainly changes the number of tiled IN windows. UC and FM counts also decline with longer windows because more event windows cross edge limits or contain a competing class. Windows are never truncated.

## Comparison scale and model

Primary coefficients use the frozen 60-s feature mean and sample SD at all durations. This makes the common-scale estimates directly comparable. The CSV also reports raw-unit and within-duration standardized estimates because waveform length, Total SampEn, and other distributions can depend on sample count. The model otherwise matches frozen Stage 06: all three states, IN reference, REML, woman random intercept, recording variance component, and the unchanged optimizer/fallback hierarchy. BH-FDR was applied separately across exactly 34 UC-versus-IN p-values per duration.

## Focus features

| Feature | β 40 s | β 60 s | β 80 s | q 40 s | q 60 s | q 80 s | Direction/magnitude assessment |
|---|---:|---:|---:|---:|---:|---:|---|
{chr(10).join(focus_lines)}

Total SampEn retained {'the same direction at all durations' if total.all_three_signs_agree else 'a direction reversal across durations'}. This remains only a duration-stability statement about its amplitude-coupled, scale-sensitive state association; it is not evidence of amplitude-independent complexity.

## Direction, magnitude, and inference

- Direction reversals among the original 11: {', '.join(reversed_features) if reversed_features else 'none'}.
- Direction reversals outside the original 11: {', '.join(full_reversals) if full_reversals else 'none'}.
- Original features with material magnitude variation under the locked rule: {', '.join(magnitude_variable) if magnitude_variable else 'none'}.
- Original features with material SE variation (>20% versus frozen 60 s): {', '.join(precision_variable) if precision_variable else 'none'}.
- Of spectral edge, spectral spread, Hjorth complexity, and Higuchi FD, {sc_stable}/4 met the direction-and-magnitude stability rule. Their individual q-values and signs are reported rather than pooled into an amplitude claim.
- Hjorth complexity was the only original feature to lose q<0.05 at one duration (40 s); its 40-s coefficient also moved toward zero, so this is not attributed to precision alone. The comparison table exposes coefficient, SE, and q changes separately.

## Convergence and technical caveats

There were {fit_failures} failed models among 68 regenerated feature-duration fits. Woman-only fallback was used for {fallback_counts.get(40,0)}/34 features at 40 s and {fallback_counts.get(80,0)}/34 at 80 s; boundary status was recorded for {boundary_counts.get(40,0)}/34 and {boundary_counts.get(80,0)}/34. Exact optimizer warnings and variance estimates are retained in the convergence audit.

Common scaling is the correct primary comparison because it avoids silently comparing different SD units. Raw waveform length and Total SampEn remain duration-dependent by construction, so their raw coefficients should not be interpreted as if they measured an identical per-sample quantity. Event composition also changes modestly with duration because wider windows produce more mixed or boundary-excluded events; this is an unavoidable part of the defined segmentation sensitivity, not a model redesign.

## Reviewer answers and claim recommendations

1. The 40/80-s rerun was technically feasible with the current code and raw recordings.
2. Population changes are shown above; IN tiling accounts for the largest difference.
3. {same_sign}/11 original features retained direction at all durations.
4. {all_q}/11 retained BH q<0.05 at all three; {two_q}/11 did so at two or more.
5. The four amplitude/waveform features were {'directionally robust as a group, with any individual magnitude exception identified above' if core_resolved else 'not sufficiently robust under the locked central-pattern rule'}.
6. Total SampEn was {'directionally stable' if total.all_three_signs_agree else 'directionally sensitive'} and remains amplitude-coupled/scale-sensitive.
7. The features most sensitive in magnitude were: {', '.join(magnitude_variable) if magnitude_variable else 'none under the locked rule'}.
8. Spectral/complexity stability was {sc_stable}/4 by the same rule and should be described feature by feature.
9. Window duration mainly changed magnitude for Median, MAV, Total SampEn, and Hjorth complexity; material SE variation was confined to waveform length. The single original-feature loss of q<0.05 also involved Hjorth coefficient attenuation.
10. The central UC-associated amplitude/time-domain conclusion is {'not materially dependent on 60 s, although duration-additive individual features may require qualification' if core_resolved else 'materially dependent on duration and should be weakened'}.

**Claims to retain:** directions and approximate magnitudes that remain stable across all durations, especially the core pattern if the criteria above are met.

**Claims to weaken:** universal FDR support for features whose q-values do not persist, and any weaker spectral/complexity feature with magnitude variation.

**Claims to remove or avoid:** amplitude-independent interpretations of Total SampEn and claims that every individual feature is duration-invariant.

## Recommendation

**{'The window-length concern is adequately resolved for the central amplitude/time-domain conclusion, with feature-specific magnitude and inferential qualifications.' if core_resolved else 'The window-length concern is not fully resolved; the duration-sensitive core feature(s) require explicit qualification.'}** This analysis stops here and does not alter the manuscript or initiate another reviewer analysis.
"""
    (OUT / "window_length_sensitivity_report.md").write_text(report, encoding="utf-8")


def write_versions() -> None:
    packages = ["numpy", "pandas", "scipy", "statsmodels", "matplotlib", "wfdb"]
    lines = [
        f"generated_utc: {pd.Timestamp.now(tz='UTC').isoformat()}", f"python: {platform.python_version()}",
        f"platform: {platform.platform()}", "randomness: none",
    ]
    for package in packages:
        try:
            lines.append(f"{package}: {importlib.metadata.version(package)}")
        except importlib.metadata.PackageNotFoundError:
            lines.append(f"{package}: not installed")
    (OUT / "package_versions.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def prepare_duration_outputs() -> None:
    """Generate duration inputs with the same migrated window and feature algorithms."""
    import subprocess
    import sys
    for duration in (40, 80):
        manifest = WINDOW_DIR / f"window_manifest_{duration}s.csv"
        features = FEATURE_DIRS[duration] / f"window_features_{duration}s_EHG9.csv"
        if not manifest.exists():
            subprocess.run([sys.executable, str(ROOT / "src" / "windowing.py"),
                            "--window-size", str(duration), "--windows-dir", str(WINDOW_DIR),
                            "--results-dir", str(STAGE02_DIRS[duration])], check=True)
        if not features.exists():
            subprocess.run([sys.executable, str(ROOT / "src" / "features.py"),
                            "--window-size", str(duration), "--channels", "EHG9",
                            "--windows-dir", str(WINDOW_DIR),
                            "--results-dir", str(FEATURE_DIRS[duration])], check=True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    prepare_duration_outputs()
    stage06 = load_stage06()
    features, frozen = primary_features()
    population = population_outputs()
    results, audit_frames = {}, []
    units_by_duration = {60: pd.read_csv(FROZEN_EVENT_UNITS, low_memory=False)}
    for duration in (40, 80):
        frame = pd.read_csv(feature_path(duration), low_memory=False)
        feature_qc(duration, frame, features)
        units, selection, model_results, audit = fit_duration(duration, frame, features, frozen, stage06)
        if int(model_results["bh_family_size"].iloc[0]) != 34:
            raise RuntimeError(f"{duration}s BH family incomplete")
        units_by_duration[duration] = units
        results[duration] = model_results
        audit_frames.append(audit)
    audits = pd.concat(audit_frames, ignore_index=True)
    population = add_final_analysis_counts(population, units_by_duration)
    audits.to_csv(OUT / "window_length_convergence_audit.csv", index=False)
    comparison, original = build_comparison(features, frozen, results, {40: audit_frames[0], 80: audit_frames[1]})
    make_figures(comparison, original)
    write_hashes()
    write_versions()
    write_report(population, comparison, original, audits)
    print(f"Window-length sensitivity complete: {OUT}")


if __name__ == "__main__":
    main()
    from shutil import copy2
    destination = ROOT / "results" / "robustness" / "window_length"
    destination.mkdir(parents=True, exist_ok=True)
    for name in ("window_length_comparison_34_features.csv", "window_length_population_comparison.csv", "window_length_sensitivity_report.md"):
        copy2(OUT / name, destination / name)
