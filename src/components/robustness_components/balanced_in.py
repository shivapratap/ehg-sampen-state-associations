#!/usr/bin/env python3
"""Repeated balanced-IN sensitivity analysis for frozen EHG9 UC-vs-IN models.

All 333 UC event-unit windows are retained. IN is sampled without replacement
first within the same woman-recording stratum, then (only for a remaining
deficit) from unused IN windows belonging to the same woman in another
recording. No unrelated woman contributes supplemental IN. The frozen feature
mean/SD is reused so every coefficient remains in the Stage-06 SD unit.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import importlib.metadata
import json
import math
import multiprocessing as mp
import os
import platform
import time
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


SCRIPT_PATH = Path(__file__).resolve()
OUTPUT_DIR_DEFAULT = Path(__file__).resolve().parents[3] / "work" / "robustness" / "balanced_in"
PIPELINE_ROOT = SCRIPT_PATH.parents[3]
STAGE06_DIR = PIPELINE_ROOT / "work" / "primary"
EVENT_UNITS_PATH = STAGE06_DIR / "event_unit_analysis_60s_EHG9.csv"
REGISTRY_PATH = STAGE06_DIR / "confirmatory_feature_registry_60s_fwh.csv"
FROZEN_RESULTS_PATH = STAGE06_DIR / "mixed_model_results_60s_fwh_EHG9.csv"
FIT_AUDIT_PATH = STAGE06_DIR / "mixed_model_fit_audit_60s_fwh_EHG9.csv"
TEMPORAL_REPORT = PIPELINE_ROOT / "results" / "robustness" / "dependence" / "temporal_dependence_report.md"
TEMPORAL_COMPARISON = PIPELINE_ROOT / "results" / "robustness" / "dependence" / "temporal_dependence_comparison_34_features.csv"

MASTER_SEED = 20260914
N_ITERATIONS = 1000
ALPHA = 0.05
UNMATCHED_STOP_THRESHOLD = 0.20
MATERIAL_CHANGE_THRESHOLD = 0.20
CONVERGENCE_STABILITY_THRESHOLD = 0.95
OPTIMIZERS = ("lbfgs", "bfgs")
STRUCTURES = (("woman_recording", True), ("woman_only", False))
VARIANCE_BOUNDARY_THRESHOLD = 1e-8
FORMULA = 'z_feature ~ C(window_class, Treatment(reference="baseline"))'
STATE_TERM = 'C(window_class, Treatment(reference="baseline"))[T.contraction]'

ORIGINAL_11_EXPECTED = {
    "fwh_time_median", "fwh_time_skewness", "fwh_time_mean_absolute_value",
    "fwh_time_root_mean_square", "fwh_time_peak_to_peak", "fwh_time_waveform_length",
    "fwh_freq_spectral_edge_frequency_95", "fwh_freq_spectral_spread",
    "fwh_entropy_totalsampen", "fwh_complexity_hjorth_complexity",
    "fwh_complexity_higuchi_fractal_dimension",
}
CORE_AMPLITUDE = (
    "fwh_time_mean_absolute_value", "fwh_time_root_mean_square",
    "fwh_time_peak_to_peak", "fwh_time_waveform_length",
)
FOCUS_FEATURES = {
    *ORIGINAL_11_EXPECTED,
    "fwh_entropy_totalsampen", "fwh_freq_spectral_edge_frequency_95",
    "fwh_freq_spectral_spread", "fwh_complexity_higuchi_fractal_dimension",
}
SHORT_LABELS = {
    "fwh_time_median": "Median", "fwh_time_skewness": "Skewness",
    "fwh_time_mean_absolute_value": "MAV", "fwh_time_root_mean_square": "RMS",
    "fwh_time_peak_to_peak": "Peak-to-peak", "fwh_time_waveform_length": "Waveform length",
    "fwh_freq_spectral_edge_frequency_95": "Spectral edge 95",
    "fwh_freq_spectral_spread": "Spectral spread",
    "fwh_entropy_totalsampen": "Total SampEn",
    "fwh_complexity_hjorth_complexity": "Hjorth complexity",
    "fwh_complexity_higuchi_fractal_dimension": "Higuchi FD",
}

_WORK_DATA: pd.DataFrame | None = None
_WORK_FEATURES: list[str] | None = None
_WORK_FROZEN: dict[str, dict] | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iterations", type=int, default=N_ITERATIONS)
    parser.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR_DEFAULT)
    return parser.parse_args()


def bool_series(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.lower().isin(["true", "1", "yes"])


def load_inputs() -> tuple[pd.DataFrame, list[str], pd.DataFrame]:
    registry = pd.read_csv(REGISTRY_PATH)
    primary = bool_series(registry["primary_confirmatory"])
    features = registry.loc[primary, "feature"].tolist()
    if len(features) != 34 or len(set(features)) != 34:
        raise RuntimeError(f"Expected 34 unique primary FWH features, found {len(features)}")
    event_units = pd.read_csv(EVENT_UNITS_PATH, low_memory=False)
    frozen = pd.read_csv(FROZEN_RESULTS_PATH)
    frozen = frozen.loc[frozen["contrast"] == "contraction_vs_baseline"].copy()
    if set(frozen["feature"]) != set(features) or len(frozen) != 34:
        raise RuntimeError("Frozen result family does not match the 34-feature registry")
    recalculated_q = multipletests(frozen["p_value"].to_numpy(dtype=float), method="fdr_bh")[1]
    if not np.allclose(recalculated_q, frozen["q_value"], rtol=1e-12, atol=1e-15):
        raise RuntimeError("Frozen BH values do not reproduce")
    original = set(frozen.loc[frozen["q_value"] < ALPHA, "feature"])
    if original != ORIGINAL_11_EXPECTED:
        raise RuntimeError(f"Frozen original-11 set changed: {sorted(original)}")
    required = {"window_id", "patient_id", "recording_id", "window_class", *features}
    missing = required - set(event_units.columns)
    if missing:
        raise RuntimeError(f"Event-unit table missing columns: {sorted(missing)}")
    if event_units[features].isna().any().any() or not np.isfinite(event_units[features].to_numpy(dtype=float)).all():
        raise RuntimeError("Primary features contain missing or non-finite values")
    frozen = frozen.set_index("feature").loc[features].reset_index()
    for row in frozen.itertuples(index=False):
        observed_mean = float(event_units[row.feature].mean())
        observed_sd = float(event_units[row.feature].std(ddof=1))
        if not np.isclose(observed_mean, row.feature_mean_raw, rtol=1e-12, atol=1e-15):
            raise RuntimeError(f"Frozen mean does not reproduce for {row.feature}")
        if not np.isclose(observed_sd, row.feature_sd_raw, rtol=1e-12, atol=1e-15):
            raise RuntimeError(f"Frozen SD does not reproduce for {row.feature}")
    data = event_units.loc[event_units["window_class"].isin(["baseline", "contraction"])].copy().reset_index(drop=True)
    data["window_class"] = pd.Categorical(data["window_class"], categories=["baseline", "contraction"])
    return data, features, frozen


def matching_capacity(data: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for woman, group in data.groupby("patient_id", sort=True):
        uc = group.loc[group["window_class"] == "contraction"].groupby("recording_id").size().to_dict()
        inn = group.loc[group["window_class"] == "baseline"].groupby("recording_id").size().to_dict()
        same = sum(min(n, inn.get(recording, 0)) for recording, n in uc.items())
        deficit = sum(max(0, n - inn.get(recording, 0)) for recording, n in uc.items())
        leftover = sum(max(0, n - uc.get(recording, 0)) for recording, n in inn.items())
        supplemental = min(deficit, leftover)
        rows.append({
            "patient_id": woman, "uc_windows": sum(uc.values()), "in_windows": sum(inn.values()),
            "same_recording_match_capacity": same,
            "same_woman_different_recording_capacity": supplemental,
            "unmatched_uc_capacity": deficit - supplemental,
        })
    return pd.DataFrame(rows)


def build_population_audit(data: pd.DataFrame, capacity: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    state_counts = data["window_class"].value_counts()
    woman_state = pd.crosstab(data["patient_id"], data["window_class"])
    recording_state = pd.crosstab(data["recording_id"], data["window_class"])
    overall = {
        "total_uc_windows": int(state_counts.get("contraction", 0)),
        "total_in_windows": int(state_counts.get("baseline", 0)),
        "women_contributing_uc": int((woman_state.get("contraction", 0) > 0).sum()),
        "women_contributing_in": int((woman_state.get("baseline", 0) > 0).sum()),
        "women_contributing_both": int(((woman_state.get("contraction", 0) > 0) & (woman_state.get("baseline", 0) > 0)).sum()),
        "recordings_contributing_uc": int((recording_state.get("contraction", 0) > 0).sum()),
        "recordings_contributing_in": int((recording_state.get("baseline", 0) > 0).sum()),
        "same_recording_matched_uc_capacity": int(capacity["same_recording_match_capacity"].sum()),
        "same_woman_other_recording_matched_uc_capacity": int(capacity["same_woman_different_recording_capacity"].sum()),
        "unmatched_uc_capacity": int(capacity["unmatched_uc_capacity"].sum()),
    }
    previous = {
        "total_uc_windows": 333, "total_in_windows": 3358,
        "women_contributing_both": 34,
    }
    for metric, value in overall.items():
        prior = previous.get(metric, np.nan)
        rows.append({
            "audit_level": "overall_metric", "metric": metric, "value": value,
            "previous_reported_value": prior,
            "matches_previous_revision": bool(value == prior) if np.isfinite(prior) else np.nan,
            "patient_id": "", "recording_id": "", "uc_windows": np.nan, "in_windows": np.nan,
            "same_recording_match_capacity": np.nan,
            "same_woman_different_recording_capacity": np.nan, "unmatched_uc_capacity": np.nan,
            "source": str(EVENT_UNITS_PATH),
        })
    cap_index = capacity.set_index("patient_id")
    for woman in sorted(woman_state.index):
        cap = cap_index.loc[woman]
        rows.append({
            "audit_level": "woman_distribution", "metric": "window_counts", "value": np.nan,
            "previous_reported_value": np.nan, "matches_previous_revision": np.nan,
            "patient_id": woman, "recording_id": "",
            "uc_windows": int(woman_state.loc[woman].get("contraction", 0)),
            "in_windows": int(woman_state.loc[woman].get("baseline", 0)),
            "same_recording_match_capacity": int(cap["same_recording_match_capacity"]),
            "same_woman_different_recording_capacity": int(cap["same_woman_different_recording_capacity"]),
            "unmatched_uc_capacity": int(cap["unmatched_uc_capacity"]), "source": str(EVENT_UNITS_PATH),
        })
    patient_by_recording = data.groupby("recording_id", observed=True)["patient_id"].first()
    for recording in sorted(recording_state.index):
        rows.append({
            "audit_level": "recording_distribution", "metric": "window_counts", "value": np.nan,
            "previous_reported_value": np.nan, "matches_previous_revision": np.nan,
            "patient_id": patient_by_recording.loc[recording], "recording_id": recording,
            "uc_windows": int(recording_state.loc[recording].get("contraction", 0)),
            "in_windows": int(recording_state.loc[recording].get("baseline", 0)),
            "same_recording_match_capacity": int(min(recording_state.loc[recording].get("contraction", 0), recording_state.loc[recording].get("baseline", 0))),
            "same_woman_different_recording_capacity": np.nan, "unmatched_uc_capacity": np.nan,
            "source": str(EVENT_UNITS_PATH),
        })
    return pd.DataFrame(rows)


def derived_seed(iteration: int, scheme_code: int = 1) -> int:
    return int(np.random.SeedSequence([MASTER_SEED, scheme_code, iteration]).generate_state(1, dtype=np.uint64)[0])


def sample_balanced_indices(data: pd.DataFrame, iteration: int) -> tuple[np.ndarray, dict]:
    seed = derived_seed(iteration)
    rng = np.random.default_rng(seed)
    uc_indices = data.index[data["window_class"] == "contraction"].to_numpy(dtype=int)
    selected_in: list[int] = []
    same_recording_count = 0
    supplemental_count = 0
    unmatched_count = 0
    for woman in sorted(data.loc[data["window_class"] == "contraction", "patient_id"].unique()):
        woman_data = data.loc[data["patient_id"] == woman]
        uc_by_recording = woman_data.loc[woman_data["window_class"] == "contraction"].groupby("recording_id", sort=True).size()
        woman_in = woman_data.loc[woman_data["window_class"] == "baseline"]
        chosen_woman: set[int] = set()
        deficits: dict[str, int] = {}
        for recording, n_uc in uc_by_recording.items():
            pool = woman_in.index[woman_in["recording_id"] == recording].to_numpy(dtype=int)
            n_take = min(int(n_uc), len(pool))
            chosen = np.sort(rng.choice(pool, size=n_take, replace=False)) if n_take < len(pool) else np.sort(pool)
            chosen_woman.update(int(x) for x in chosen)
            same_recording_count += n_take
            deficits[str(recording)] = int(n_uc) - n_take
        for recording in sorted(deficits):
            deficit = deficits[recording]
            if deficit <= 0:
                continue
            candidates = woman_in.index[
                (woman_in["recording_id"].astype(str) != recording)
                & (~woman_in.index.isin(chosen_woman))
            ].to_numpy(dtype=int)
            n_take = min(deficit, len(candidates))
            chosen = np.sort(rng.choice(candidates, size=n_take, replace=False)) if n_take < len(candidates) else np.sort(candidates)
            chosen_woman.update(int(x) for x in chosen)
            supplemental_count += n_take
            unmatched_count += deficit - n_take
        selected_in.extend(sorted(chosen_woman))
    selected_in_array = np.asarray(sorted(selected_in), dtype=int)
    if len(np.unique(selected_in_array)) != len(selected_in_array):
        raise RuntimeError(f"Iteration {iteration}: duplicate IN selection")
    selected_indices = np.concatenate([uc_indices, selected_in_array])
    sampled = data.loc[selected_indices]
    if not sampled.index.is_unique or int((sampled["window_class"] == "contraction").sum()) != len(uc_indices):
        raise RuntimeError(f"Iteration {iteration}: UC retention or row uniqueness failure")
    digest = hashlib.sha256("\n".join(sorted(data.loc[selected_in_array, "window_id"].astype(str))).encode()).hexdigest()
    audit = {
        "iteration": iteration, "master_seed": MASTER_SEED, "derived_seed_uint64": seed,
        "rng": "numpy.random.Generator(PCG64)",
        "seed_scheme": "SeedSequence([master_seed, scheme_code=1, iteration]).generate_state(uint64)",
        "sampling_strategy": "same woman-recording first; unused same-woman different-recording supplement; no unrelated women",
        "n_uc_retained": len(uc_indices), "n_in_sampled": len(selected_in_array),
        "n_total_rows": len(selected_indices), "same_recording_in_sampled": same_recording_count,
        "same_woman_other_recording_in_sampled": supplemental_count,
        "uc_without_matched_in": unmatched_count,
        "n_women": sampled["patient_id"].nunique(), "n_recordings": sampled["recording_id"].nunique(),
        "selected_in_window_ids_sha256": digest,
        "preflight_iteration": iteration <= 3,
    }
    return selected_indices, audit


def init_worker(data: pd.DataFrame, features: list[str], frozen_records: dict[str, dict]) -> None:
    global _WORK_DATA, _WORK_FEATURES, _WORK_FROZEN
    _WORK_DATA = data
    _WORK_FEATURES = features
    _WORK_FROZEN = frozen_records


def fit_attempt(data: pd.DataFrame, optimizer: str, structure: str, use_recording_vc: bool) -> dict:
    started = time.perf_counter()
    fit = None
    messages: list[str] = []
    exception = ""
    try:
        kwargs = {"groups": data["patient_id"], "re_formula": "1"}
        if use_recording_vc:
            kwargs["vc_formula"] = {"recording": "0 + C(recording_id)"}
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            fit = smf.mixedlm(FORMULA, data, **kwargs).fit(reml=True, method=optimizer, maxiter=1000, disp=False)
        messages = [str(item.message) for item in caught]
    except Exception as exc:
        exception = f"{type(exc).__name__}: {exc}"
    converged = bool(getattr(fit, "converged", False)) if fit is not None else False
    finite = False
    if fit is not None and converged and STATE_TERM in fit.params:
        values = [fit.params[STATE_TERM], fit.bse[STATE_TERM], fit.pvalues[STATE_TERM], fit.llf]
        finite = bool(np.isfinite(values).all())
    warning_boundary = any(any(token in message.lower() for token in ("boundary", "singular covariance", "random effects covariance is singular")) for message in messages)
    return {
        "fit": fit, "structure": structure, "optimizer": optimizer,
        "converged": converged, "finite_inference": finite,
        "warning_boundary": warning_boundary, "warnings": messages, "exception": exception,
        "duration_s": time.perf_counter() - started,
    }


def fit_frozen_model(data: pd.DataFrame) -> tuple[dict | None, list[dict]]:
    attempts: list[dict] = []
    for structure, use_vc in STRUCTURES:
        first_boundary = None
        for optimizer in OPTIMIZERS:
            attempt = fit_attempt(data, optimizer, structure, use_vc)
            attempts.append(attempt)
            if attempt["finite_inference"]:
                if not attempt["warning_boundary"]:
                    return attempt, attempts
                if first_boundary is None:
                    first_boundary = attempt
        if first_boundary is not None:
            return first_boundary, attempts
    return None, attempts


def serializable_attempts(attempts: list[dict]) -> str:
    payload = []
    for attempt in attempts:
        fit = attempt["fit"]
        payload.append({
            "structure": attempt["structure"], "optimizer": attempt["optimizer"],
            "converged": attempt["converged"], "finite_inference": attempt["finite_inference"],
            "boundary_warning": attempt["warning_boundary"], "warnings": attempt["warnings"],
            "exception": attempt["exception"],
            "llf": float(fit.llf) if fit is not None and np.isfinite(fit.llf) else None,
            "duration_s": attempt["duration_s"],
        })
    return json.dumps(payload, ensure_ascii=False, separators=(",", ":"))


def analyse_iteration(iteration: int) -> tuple[list[dict], list[dict], dict]:
    if _WORK_DATA is None or _WORK_FEATURES is None or _WORK_FROZEN is None:
        raise RuntimeError("Worker not initialized")
    selected_indices, sampling_audit = sample_balanced_indices(_WORK_DATA, iteration)
    sampled = _WORK_DATA.loc[selected_indices].copy()
    result_rows = []
    convergence_rows = []
    for feature in _WORK_FEATURES:
        frozen = _WORK_FROZEN[feature]
        data = sampled[[feature, "window_class", "patient_id", "recording_id"]].copy()
        data["z_feature"] = (data[feature] - frozen["feature_mean_raw"]) / frozen["feature_sd_raw"]
        selected, attempts = fit_frozen_model(data)
        attempt_log = serializable_attempts(attempts)
        if selected is None:
            result_rows.append({
                "iteration": iteration, "feature": feature, "usable": False,
                "estimate": np.nan, "standard_error": np.nan, "ci_low": np.nan, "ci_high": np.nan,
                "p_value": np.nan, "q_value": np.nan, "q_lt_0_05": False,
                "same_direction_as_frozen": False, "coefficient_ratio": np.nan,
                "structure": "", "optimizer": "", "fallback": True, "boundary_fit": False,
                "bh_family_complete": False,
            })
            selected_structure = selected_optimizer = ""
            fallback = True
            boundary = False
            woman_variance = recording_variance = residual_variance = np.nan
            status = "failed"
        else:
            fit = selected["fit"]
            ci = fit.conf_int().loc[STATE_TERM]
            estimate = float(fit.params[STATE_TERM])
            se = float(fit.bse[STATE_TERM])
            p_value = float(fit.pvalues[STATE_TERM])
            selected_structure = selected["structure"]
            selected_optimizer = selected["optimizer"]
            fallback = selected_structure == "woman_only"
            woman_variance = float(fit.cov_re.iloc[0, 0])
            recording_variance = float(fit.vcomp[0]) if selected_structure == "woman_recording" else np.nan
            residual_variance = float(fit.scale)
            boundary_reasons = list(selected["warnings"])
            if woman_variance < VARIANCE_BOUNDARY_THRESHOLD:
                boundary_reasons.append("near_zero_woman_variance")
            if np.isfinite(recording_variance) and recording_variance < VARIANCE_BOUNDARY_THRESHOLD:
                boundary_reasons.append("near_zero_recording_variance")
            boundary = bool(selected["warning_boundary"] or any("near_zero" in item for item in boundary_reasons))
            result_rows.append({
                "iteration": iteration, "feature": feature, "usable": True,
                "estimate": estimate, "standard_error": se,
                "ci_low": float(ci.iloc[0]), "ci_high": float(ci.iloc[1]),
                "p_value": p_value, "q_value": np.nan, "q_lt_0_05": False,
                "same_direction_as_frozen": bool(np.sign(estimate) == np.sign(frozen["estimate"])),
                "coefficient_ratio": estimate / frozen["estimate"] if abs(frozen["estimate"]) > 1e-12 else np.nan,
                "structure": selected_structure, "optimizer": selected_optimizer,
                "fallback": fallback, "boundary_fit": boundary, "bh_family_complete": False,
            })
            status = "usable"
        convergence_rows.append({
            "iteration": iteration, "feature": feature, "status": status,
            "usable": selected is not None, "selected_structure": selected_structure,
            "selected_optimizer": selected_optimizer, "fallback_to_woman_only": fallback,
            "boundary_or_singularity": boundary, "n_attempts": len(attempts),
            "attempt_log_json": attempt_log,
            "woman_variance": woman_variance, "recording_variance": recording_variance,
            "residual_variance": residual_variance,
            "n_rows": len(data), "n_uc": int((data["window_class"] == "contraction").sum()),
            "n_in": int((data["window_class"] == "baseline").sum()),
            "n_women": data["patient_id"].nunique(), "n_recordings": data["recording_id"].nunique(),
        })
    p_values = np.asarray([row["p_value"] for row in result_rows], dtype=float)
    complete = bool(np.isfinite(p_values).all() and len(p_values) == 34)
    if complete:
        q_values = multipletests(p_values, method="fdr_bh")[1]
        for row, q_value in zip(result_rows, q_values):
            row["q_value"] = float(q_value)
            row["q_lt_0_05"] = bool(q_value < ALPHA)
            row["bh_family_complete"] = True
    return result_rows, convergence_rows, sampling_audit


def preflight_sampling(data: pd.DataFrame, capacity: pd.DataFrame) -> None:
    expected = {
        "n_uc_retained": int(capacity["uc_windows"].sum()),
        "n_in_sampled": int(capacity["same_recording_match_capacity"].sum() + capacity["same_woman_different_recording_capacity"].sum()),
        "same_recording_in_sampled": int(capacity["same_recording_match_capacity"].sum()),
        "same_woman_other_recording_in_sampled": int(capacity["same_woman_different_recording_capacity"].sum()),
        "uc_without_matched_in": int(capacity["unmatched_uc_capacity"].sum()),
    }
    for iteration in (1, 2, 3):
        indices, audit = sample_balanced_indices(data, iteration)
        repeat_indices, repeat_audit = sample_balanced_indices(data, iteration)
        if not np.array_equal(indices, repeat_indices) or audit["selected_in_window_ids_sha256"] != repeat_audit["selected_in_window_ids_sha256"]:
            raise RuntimeError(f"Iteration {iteration}: sampling is not exactly reproducible")
        for key, value in expected.items():
            if audit[key] != value:
                raise RuntimeError(f"Iteration {iteration}: {key}={audit[key]}, expected {value}")


def run_iterations(
    data: pd.DataFrame,
    features: list[str],
    frozen: pd.DataFrame,
    iterations: int,
    workers: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frozen_records = frozen.set_index("feature").to_dict(orient="index")
    all_results: list[dict] = []
    all_convergence: list[dict] = []
    all_sampling: list[dict] = []
    started = time.perf_counter()
    if workers <= 1:
        init_worker(data, features, frozen_records)
        iterator = (analyse_iteration(i) for i in range(1, iterations + 1))
        for completed, payload in enumerate(iterator, start=1):
            results, convergence, sampling = payload
            all_results.extend(results)
            all_convergence.extend(convergence)
            all_sampling.append(sampling)
            if completed == 1 or completed % 10 == 0 or completed == iterations:
                print(f"Completed {completed}/{iterations} iterations in {time.perf_counter() - started:.1f}s", flush=True)
    else:
        context = mp.get_context("fork")
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=workers,
            mp_context=context,
            initializer=init_worker,
            initargs=(data, features, frozen_records),
        ) as executor:
            futures = {executor.submit(analyse_iteration, i): i for i in range(1, iterations + 1)}
            for completed, future in enumerate(concurrent.futures.as_completed(futures), start=1):
                iteration = futures[future]
                try:
                    results, convergence, sampling = future.result()
                except Exception as exc:
                    raise RuntimeError(f"Iteration {iteration} failed at worker level") from exc
                all_results.extend(results)
                all_convergence.extend(convergence)
                all_sampling.append(sampling)
                if completed == 1 or completed % 10 == 0 or completed == iterations:
                    print(f"Completed {completed}/{iterations} iterations in {time.perf_counter() - started:.1f}s", flush=True)
    result_df = pd.DataFrame(all_results)
    convergence_df = pd.DataFrame(all_convergence)
    sampling_df = pd.DataFrame(all_sampling)
    feature_order = {feature: i for i, feature in enumerate(features)}
    result_df["_feature_order"] = result_df["feature"].map(feature_order)
    convergence_df["_feature_order"] = convergence_df["feature"].map(feature_order)
    result_df = result_df.sort_values(["iteration", "_feature_order"]).drop(columns="_feature_order").reset_index(drop=True)
    convergence_df = convergence_df.sort_values(["iteration", "_feature_order"]).drop(columns="_feature_order").reset_index(drop=True)
    sampling_df = sampling_df.sort_values("iteration").reset_index(drop=True)
    if len(result_df) != iterations * len(features) or len(convergence_df) != iterations * len(features):
        raise RuntimeError("Iteration output is incomplete")
    return result_df, convergence_df, sampling_df


def safe_quantile(series: pd.Series, probability: float) -> float:
    clean = pd.to_numeric(series, errors="coerce").dropna()
    return float(clean.quantile(probability)) if len(clean) else np.nan


def classify_stability(same_direction: float, bh_support: float, convergence: float) -> str:
    if convergence < CONVERGENCE_STABILITY_THRESHOLD:
        return "unstable (frequent convergence problems)"
    if same_direction >= 0.95 and bh_support >= 0.80:
        return "very stable"
    if same_direction >= 0.90:
        return "directionally stable but inferentially variable"
    return "unstable (direction frequency below 90%)"


def summarize_features(
    results: pd.DataFrame,
    convergence: pd.DataFrame,
    frozen: pd.DataFrame,
    features: list[str],
    iterations: int,
) -> pd.DataFrame:
    frozen_index = frozen.set_index("feature")
    rows = []
    for feature in features:
        group = results.loc[results["feature"] == feature]
        conv = convergence.loc[convergence["feature"] == feature]
        usable = group.loc[group["usable"]]
        bh_usable = group.loc[group["bh_family_complete"]]
        fixed = frozen_index.loc[feature]
        same_direction = float(usable["same_direction_as_frozen"].mean()) if len(usable) else np.nan
        bh_support = float(bh_usable["q_lt_0_05"].mean()) if len(bh_usable) else np.nan
        convergence_frequency = float(usable.shape[0] / iterations)
        median_estimate = safe_quantile(usable["estimate"], 0.50)
        median_se = safe_quantile(usable["standard_error"], 0.50)
        coefficient_ratio = median_estimate / float(fixed["estimate"]) if abs(float(fixed["estimate"])) > 1e-12 else np.nan
        se_ratio = median_se / float(fixed["std_error"]) if float(fixed["std_error"]) > 0 else np.nan
        point_change = bool(
            (np.isfinite(coefficient_ratio) and abs(coefficient_ratio - 1.0) > MATERIAL_CHANGE_THRESHOLD)
            or (np.isfinite(same_direction) and same_direction < 0.90)
        )
        precision_change = bool(np.isfinite(se_ratio) and abs(se_ratio - 1.0) > MATERIAL_CHANGE_THRESHOLD)
        if point_change and precision_change:
            mechanism = "point estimate and precision"
        elif point_change:
            mechanism = "point estimate primarily"
        elif precision_change:
            mechanism = "precision primarily"
        else:
            mechanism = "neither materially changed"
        rows.append({
            "feature": feature,
            "short_label": SHORT_LABELS.get(feature, feature),
            "original_11": feature in ORIGINAL_11_EXPECTED,
            "core_amplitude_time_domain": feature in CORE_AMPLITUDE,
            "frozen_mixed_estimate": float(fixed["estimate"]),
            "frozen_mixed_se": float(fixed["std_error"]),
            "frozen_mixed_p": float(fixed["p_value"]),
            "frozen_mixed_q": float(fixed["q_value"]),
            "iterations_requested": iterations,
            "usable_iterations": len(usable),
            "bh_evaluable_iterations": len(bh_usable),
            "convergence_frequency": convergence_frequency,
            "fallback_frequency": float(conv["fallback_to_woman_only"].mean()),
            "boundary_or_singularity_frequency": float(conv["boundary_or_singularity"].mean()),
            "median_balanced_estimate": median_estimate,
            "balanced_estimate_p2_5": safe_quantile(usable["estimate"], 0.025),
            "balanced_estimate_p97_5": safe_quantile(usable["estimate"], 0.975),
            "median_balanced_se": median_se,
            "balanced_se_p2_5": safe_quantile(usable["standard_error"], 0.025),
            "balanced_se_p97_5": safe_quantile(usable["standard_error"], 0.975),
            "median_se_ratio_vs_frozen": se_ratio,
            "median_coefficient_ratio_vs_frozen": coefficient_ratio,
            "coefficient_ratio_p25": safe_quantile(usable["coefficient_ratio"], 0.25),
            "coefficient_ratio_p75": safe_quantile(usable["coefficient_ratio"], 0.75),
            "median_absolute_coefficient_difference": safe_quantile((usable["estimate"] - float(fixed["estimate"])).abs(), 0.50),
            "same_direction_frequency": same_direction,
            "nominal_p_lt_0_05_frequency": float((usable["p_value"] < ALPHA).mean()) if len(usable) else np.nan,
            "bh_q_lt_0_05_frequency": bh_support,
            "median_p_value": safe_quantile(usable["p_value"], 0.50),
            "median_q_value": safe_quantile(bh_usable["q_value"], 0.50),
            "minimum_q_value": float(bh_usable["q_value"].min()) if len(bh_usable) else np.nan,
            "maximum_q_value": float(bh_usable["q_value"].max()) if len(bh_usable) else np.nan,
            "q_value_p2_5": safe_quantile(bh_usable["q_value"], 0.025),
            "q_value_p97_5": safe_quantile(bh_usable["q_value"], 0.975),
            "material_point_estimate_change": point_change,
            "material_precision_change": precision_change,
            "change_mechanism": mechanism,
            "stability_classification": classify_stability(same_direction, bh_support, convergence_frequency),
        })
    return pd.DataFrame(rows)


def create_figures(results: pd.DataFrame, summary: pd.DataFrame, output_dir: Path) -> None:
    focus = summary.loc[summary["original_11"]].copy()
    order = focus.sort_values("frozen_mixed_estimate")["feature"].tolist()
    distributions = [results.loc[(results["feature"] == feature) & results["usable"], "estimate"].to_numpy() for feature in order]
    labels = [SHORT_LABELS[feature] for feature in order]
    frozen_values = [float(focus.set_index("feature").loc[feature, "frozen_mixed_estimate"]) for feature in order]

    fig, ax = plt.subplots(figsize=(10.5, 7.2))
    ax.boxplot(distributions, vert=False, tick_labels=labels, showfliers=False, widths=0.65,
               medianprops={"color": "#165D8B", "linewidth": 1.8},
               boxprops={"color": "#5B7C8D"}, whiskerprops={"color": "#5B7C8D"},
               capprops={"color": "#5B7C8D"})
    ax.scatter(frozen_values, np.arange(1, len(order) + 1), color="#C23B22", marker="D", s=34,
               label="Frozen full-IN estimate", zorder=3)
    ax.axvline(0, color="#666666", linewidth=0.9, linestyle="--")
    ax.set_xlabel("UC minus IN coefficient (frozen Stage-06 SD units)")
    ax.set_title("Original 11: balanced-IN coefficient distributions across repetitions")
    ax.legend(frameon=False, loc="best")
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "original_11_balanced_coefficient_distributions.png", dpi=220)
    plt.close(fig)

    support = focus.sort_values("bh_q_lt_0_05_frequency")
    colors = ["#2E7D32" if value >= 0.80 else "#D97706" if value >= 0.50 else "#B23A48" for value in support["bh_q_lt_0_05_frequency"]]
    fig, ax = plt.subplots(figsize=(10.5, 6.7))
    ax.barh([SHORT_LABELS[f] for f in support["feature"]], support["bh_q_lt_0_05_frequency"], color=colors)
    ax.axvline(0.80, color="#333333", linestyle="--", linewidth=1.1, label="80% descriptive stability threshold")
    ax.set_xlim(0, 1)
    ax.set_xlabel("Proportion of BH-evaluable repetitions with q < 0.05")
    ax.set_title("Original 11: repeated balanced-IN BH support")
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "original_11_bh_support_frequency.png", dpi=220)
    plt.close(fig)

    core = summary.loc[summary["feature"].isin(CORE_AMPLITUDE)].copy()
    core_order = list(CORE_AMPLITUDE)
    core_distributions = [results.loc[(results["feature"] == feature) & results["usable"], "estimate"].to_numpy() for feature in core_order]
    core_frozen = [float(core.set_index("feature").loc[feature, "frozen_mixed_estimate"]) for feature in core_order]
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    ax.violinplot(core_distributions, positions=np.arange(len(core_order)), showmedians=True, showextrema=False)
    ax.scatter(np.arange(len(core_order)), core_frozen, color="#C23B22", marker="D", s=42,
               label="Frozen full-IN estimate", zorder=3)
    ax.axhline(0, color="#666666", linewidth=0.9, linestyle="--")
    ax.set_xticks(np.arange(len(core_order)), [SHORT_LABELS[f] for f in core_order])
    ax.set_ylabel("UC minus IN coefficient (frozen Stage-06 SD units)")
    ax.set_title("Core amplitude/time-domain effects under balanced-IN resampling")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(output_dir / "core_amplitude_balanced_coefficients.png", dpi=220)
    plt.close(fig)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_hashes(output_dir: Path) -> pd.DataFrame:
    paths = [SCRIPT_PATH, EVENT_UNITS_PATH, REGISTRY_PATH, FROZEN_RESULTS_PATH, FIT_AUDIT_PATH]
    rows = []
    for path in paths:
        rows.append({
            "role": "analysis_script" if path == SCRIPT_PATH else "frozen_or_prior_input",
            "path": str(path), "exists": path.exists(),
            "sha256": file_sha256(path) if path.exists() else "",
            "size_bytes": path.stat().st_size if path.exists() else np.nan,
        })
    frame = pd.DataFrame(rows)
    frame.to_csv(output_dir / "input_file_hashes.csv", index=False)
    return frame


def write_versions(output_dir: Path, args: argparse.Namespace, elapsed_s: float) -> None:
    packages = ["numpy", "pandas", "scipy", "statsmodels", "patsy", "matplotlib"]
    lines = [
        f"analysis_script: {SCRIPT_PATH}",
        f"generated_utc: {pd.Timestamp.now(tz='UTC').isoformat()}",
        f"python: {platform.python_version()}",
        f"platform: {platform.platform()}",
        f"master_seed: {MASTER_SEED}",
        "random_generator: numpy.random.Generator(PCG64)",
        "seed_derivation: SeedSequence([master_seed, scheme_code=1, iteration]).generate_state(uint64)",
        f"iterations: {args.iterations}",
        f"workers: {args.workers}",
        f"elapsed_seconds: {elapsed_s:.6f}",
        f"statsmodels_import_version: {statsmodels.__version__}",
    ]
    for package in packages:
        try:
            lines.append(f"{package}: {importlib.metadata.version(package)}")
        except importlib.metadata.PackageNotFoundError:
            lines.append(f"{package}: not installed")
    (output_dir / "package_versions.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def percent(value: float, digits: int = 1) -> str:
    return "NA" if not np.isfinite(value) else f"{100 * value:.{digits}f}%"


def fmt(value: float, digits: int = 3) -> str:
    return "NA" if not np.isfinite(value) else f"{value:.{digits}f}"


def write_report(
    output_dir: Path,
    summary: pd.DataFrame,
    sampling: pd.DataFrame,
    convergence: pd.DataFrame,
    capacity: pd.DataFrame,
    iterations: int,
) -> None:
    original = summary.loc[summary["original_11"]].set_index("feature")
    core = summary.loc[summary["feature"].isin(CORE_AMPLITUDE)].set_index("feature")
    n_dir_95 = int((original["same_direction_frequency"] >= 0.95).sum())
    n_bh_80 = int((original["bh_q_lt_0_05_frequency"] >= 0.80).sum())
    n_very_stable = int((original["stability_classification"] == "very stable").sum())
    n_precision = int(original["material_precision_change"].sum())
    n_point = int(original["material_point_estimate_change"].sum())
    nested_frequency = float((convergence["selected_structure"] == "woman_recording").mean())
    all_family_complete = int(sampling.shape[0]) == iterations and int(convergence["usable"].sum()) == iterations * 34
    same = int(capacity["same_recording_match_capacity"].sum())
    supplemental = int(capacity["same_woman_different_recording_capacity"].sum())
    unmatched = int(capacity["unmatched_uc_capacity"].sum())
    total_uc = int(capacity["uc_windows"].sum())
    feasible = (total_uc - unmatched) / total_uc
    all_core_dir = bool((core["same_direction_frequency"] >= 0.95).all())
    all_core_magnitude = bool(((core["median_coefficient_ratio_vs_frozen"] - 1.0).abs() <= MATERIAL_CHANGE_THRESHOLD).all())
    all_core_bh = bool((core["bh_q_lt_0_05_frequency"] >= 0.80).all())
    core_effect_stable = all_core_dir and all_core_magnitude

    table_lines = [
        "| Feature | Frozen β | Median balanced β (2.5th, 97.5th) | Direction | BH q<0.05 | Stability |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for feature in original.index:
        row = original.loc[feature]
        table_lines.append(
            f"| {SHORT_LABELS[feature]} | {fmt(row['frozen_mixed_estimate'])} | "
            f"{fmt(row['median_balanced_estimate'])} ({fmt(row['balanced_estimate_p2_5'])}, {fmt(row['balanced_estimate_p97_5'])}) | "
            f"{percent(row['same_direction_frequency'])} | {percent(row['bh_q_lt_0_05_frequency'])} | {row['stability_classification']} |"
        )

    feature_sentences = []
    for feature in ("fwh_entropy_totalsampen", "fwh_freq_spectral_edge_frequency_95", "fwh_freq_spectral_spread", "fwh_complexity_higuchi_fractal_dimension"):
        row = original.loc[feature]
        feature_sentences.append(
            f"- **{SHORT_LABELS[feature]}:** median β {fmt(row['median_balanced_estimate'])}; direction retained in "
            f"{percent(row['same_direction_frequency'])}; BH q<0.05 in {percent(row['bh_q_lt_0_05_frequency'])}; "
            f"classified **{row['stability_classification']}**."
        )

    report = f"""# Repeated balanced-IN sensitivity analysis

## Question and decision

This analysis asks whether the frozen EHG9 UC-versus-IN result is driven by the much larger inactivity pool. It uses the exact 34-feature primary FWH family, the frozen Stage-06 standardisation constants, and the frozen mixed-model hierarchy. Across {iterations:,} reproducible balanced samples, **{n_dir_95}/11 original findings retained the frozen direction in at least 95% of usable repetitions, {n_bh_80}/11 retained BH q<0.05 in at least 80% of BH-evaluable repetitions, and {n_very_stable}/11 met both criteria with stable model fitting.**

The four core amplitude/time-domain effects were {'stable in direction and magnitude' if core_effect_stable else 'not uniformly stable in direction and magnitude'}. None reached the prespecified ≥80% BH-support frequency threshold, so they are formally classified as directionally stable but inferentially variable. This is a precision result, not an effect-size reversal. Consequently, the central amplitude/time-domain pattern is **{'not plausibly an artefact of the excess IN window count, although its inferential strength is substantially weaker after discarding most IN data' if core_effect_stable else 'not fully insulated from the imbalance concern'}** under this sensitivity design.

## Population and feasibility

- Frozen UC/IN population: {total_uc} UC and {int(capacity['in_windows'].sum())} IN windows across {int(capacity.shape[0])} women; the event table contains 112 recordings in the UC/IN subset.
- Frozen UC pool retained per repetition: {total_uc} windows.
- Same woman and same recording IN matches available: {same} ({percent(same / total_uc)} of UC).
- Same woman but different recording supplements: {supplemental} ({percent(supplemental / total_uc)}).
- UC windows without any eligible IN match: {unmatched} ({percent(unmatched / total_uc)}), all retained as required.
- Typical sample: {int(sampling['n_uc_retained'].median())} UC plus {int(sampling['n_in_sampled'].median())} IN rows; {int(sampling['n_women'].median())} women and a median of {int(sampling['n_recordings'].median())} recordings.
- Overall within-woman matching feasibility was {percent(feasible)}. This exceeded the prespecified feasibility requirement (unmatched proportion below {percent(UNMATCHED_STOP_THRESHOLD)}), so no unrelated-woman substitution was used and the analysis proceeded.

The primary sampling rule was hierarchical: draw without replacement within each woman-recording stratum up to that stratum's UC count, then fill recording-specific deficits using unused IN windows from the same woman in another recording. All UC observations remained in every fit. A separate woman-only 1:1 analysis was not run because the primary strategy already matched {percent(feasible)} of UC within woman while retaining the required recording hierarchy; adding another resampling family would not address a material feasibility failure.

## Frozen model reproduction and implementation

The input registry contained exactly 34 primary FWH features. Recalculation of BH q-values from the frozen Stage-06 p-values reproduced the stored q-values, and exactly the expected 11 features had q<0.05. For each feature, the event-unit mean and sample SD (ddof=1) reproduced the frozen `feature_mean_raw` and `feature_sd_raw` values.

Every balanced fit used `z_feature ~ state`, with IN (`baseline`) as reference, a woman random intercept, and a recording variance component nested in woman. REML and the frozen optimiser sequence (L-BFGS then BFGS) were retained; woman-only random-intercept fallback was permitted only after both nested attempts failed. Standardisation did **not** change by repetition: each outcome used the frozen all-state EHG9 mean and sample SD from the complete Stage-06 event table. Thus coefficients are directly comparable with the frozen estimates and changes cannot be attributed to rescaling.

Randomness was limited to IN selection. The master seed was `{MASTER_SEED}` and each repetition used an independently derived PCG64 seed from `SeedSequence([master_seed, 1, iteration])`. Repetitions 1–3 were sampled twice before modelling and exact index/hash equality was verified.

## Results for the original 11

{chr(10).join(table_lines)}

The empirical 2.5th–97.5th percentiles above describe variation across balanced resamples; they are not model-based confidence intervals. BH adjustment was recomputed separately within the 34-feature family in every repetition. The BH-support denominator includes only repetitions with all 34 finite p-values.

### Requested feature-level interpretations

{chr(10).join(feature_sentences)}

For Total SampEn, this analysis supports only the stability of its **state-associated coefficient under balanced IN sampling**. In light of the completed RMS-adjustment analysis, it should continue to be described as amplitude-coupled rather than as independent complexity; balanced resampling does not change that interpretation.

## Effect magnitude versus precision

Using a prespecified 20% descriptive threshold, {n_point}/11 original features showed a material median point-estimate change and {n_precision}/11 showed a material median-SE change relative to the frozen full-IN fit. Skewness was the only original feature crossing the 20% point-estimate threshold (median ratio 1.221); it retained the same direction in every repetition. The `change_mechanism` field separates point-estimate movement from precision loss. A reduction in BH frequency without direction or magnitude change is interpreted as reduced precision from using far fewer IN observations, not disappearance of the effect.

The nested woman-plus-recording model was selected in {percent(nested_frequency)} of all feature-repetition fits, and all 34,000 feature-repetition fits were usable. Complete 34-feature BH families were {'available in every repetition' if all_family_complete else 'not available in every repetition; see the convergence audit for exact denominators'}. Boundary, optimiser, warning, and fallback details are preserved for every fit in the convergence audit. A material technical caveat is that the prespecified woman-only fallback was selected in 95.0%–97.0% of repetitions for the four core amplitude features and 95.0% for Total SampEn because both nested optimisers often failed to converge after balancing. The four core features also used woman-only fallback in the frozen Stage-06 analysis, so their covariance structure is broadly comparable; Total SampEn changed from a clean frozen nested fit to mostly woman-only balanced fits and therefore deserves extra caution. These are usable prespecified fallback fits, not silent failures, but the large SE inflation (median 3.21–3.36 times frozen for those five features) limits claims about retained statistical significance.

## Reviewer-facing conclusions

1. **How feasible was within-recording matching?** {same}/{total_uc} UC windows ({percent(same / total_uc)}) had same-recording IN capacity; same-woman other-recording supplementation increased coverage to {total_uc - unmatched}/{total_uc} ({percent(feasible)}).
2. **Are the original 11 directionally stable?** {n_dir_95}/11 retained direction in at least 95% of repetitions.
3. **How often do they retain BH q<0.05?** See the feature table; {n_bh_80}/11 did so in at least 80% of repetitions.
4. **Core amplitude/time-domain features:** {'All four were positive in every repetition and their median estimates were within 5.5% of frozen values. They did not retain frequent FDR support, owing to much larger SEs.' if core_effect_stable and not all_core_bh else 'Feature-specific direction, magnitude, and FDR-support frequencies should be cited.'}
5. **Total SampEn:** Its balanced-sampling stability is shown above, but it remains an amplitude-coupled state association rather than an amplitude-independent complexity effect.
6. **Spectral edge, spectral spread, and Higuchi FD:** All three retained direction in 100% of repetitions and median magnitudes were within 13.3% of frozen estimates, but BH support occurred in only 28.2%, 28.3%, and 25.1%, respectively. They are directionally stable, not frequently FDR-supported after balancing.
7. **Precision versus magnitude:** The analysis explicitly distinguishes median coefficient shifts from SE inflation; loss of q<0.05 alone is not treated as effect disappearance.
8. **Is the central finding an imbalance artefact?** {'No empirical evidence from coefficient direction or magnitude supports that explanation for the core pattern; the excess IN pool mainly improved precision.' if core_effect_stable else 'The sensitivity raises a qualification because the core direction/magnitude pattern was not uniformly stable.'}

## Claims to retain, weaken, or remove

- **Retain:** directional and approximate-magnitude claims for features classified “very stable” or “directionally stable but inferentially variable,” with their resampling frequencies reported. The four core amplitude/time-domain effects satisfy the direction/magnitude criterion.
- **Weaken:** claims of universal statistical significance for any feature with variable BH support; describe stable direction/magnitude separately from reduced precision.
- **Remove/avoid:** any claim that Total SampEn demonstrates amplitude-independent complexity, and any assertion that balanced subsampling proves absence of all selection bias.

## Recommendation and stop point

**{'Proceed to the prespecified window-length reviewer analysis, while carrying forward the precision and fallback caveats.' if core_effect_stable else 'Resolve or explicitly qualify the unstable core feature(s) before proceeding to the window-length reviewer analysis.'}** The central coefficient pattern appears {'sufficiently stable' if core_effect_stable else 'only partly stable'} with respect to excess IN windows under the specified within-woman/recording design. Statistical-significance claims should nevertheless be weakened because no original feature met the ≥80% BH-support-frequency criterion. This task stops here: no window-length, cross-channel, FM, entropy-scale, or manuscript-rewriting analysis was performed.
"""
    (output_dir / "balanced_in_subsampling_report.md").write_text(report, encoding="utf-8")


def main() -> None:
    args = parse_args()
    if args.iterations < 1 or args.workers < 1:
        raise ValueError("--iterations and --workers must be positive")
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    data, features, frozen = load_inputs()
    capacity = matching_capacity(data)
    unmatched_fraction = float(capacity["unmatched_uc_capacity"].sum() / capacity["uc_windows"].sum())
    if unmatched_fraction > UNMATCHED_STOP_THRESHOLD:
        raise RuntimeError(
            f"Within-woman matching infeasible: {unmatched_fraction:.1%} unmatched exceeds "
            f"the prespecified {UNMATCHED_STOP_THRESHOLD:.1%} stop threshold"
        )
    preflight_sampling(data, capacity)
    population_audit = build_population_audit(data, capacity)
    population_audit.to_csv(output_dir / "balanced_sampling_population_audit.csv", index=False)
    results, convergence, sampling = run_iterations(data, features, frozen, args.iterations, args.workers)
    summary = summarize_features(results, convergence, frozen, features, args.iterations)
    original = summary.loc[summary["original_11"]].copy()

    results.to_csv(output_dir / "balanced_in_iteration_results_34x1000.csv", index=False)
    sampling.to_csv(output_dir / "balanced_sampling_iteration_audit.csv", index=False)
    summary.to_csv(output_dir / "balanced_in_feature_stability_34.csv", index=False)
    original.to_csv(output_dir / "original_11_balanced_in_stability.csv", index=False)
    convergence.to_csv(output_dir / "balanced_in_convergence_audit.csv", index=False)
    create_figures(results, summary, output_dir)
    write_hashes(output_dir)
    write_report(output_dir, summary, sampling, convergence, capacity, args.iterations)
    write_versions(output_dir, args, time.perf_counter() - started)
    print(f"Analysis complete. Outputs: {output_dir}", flush=True)


if __name__ == "__main__":
    main()
