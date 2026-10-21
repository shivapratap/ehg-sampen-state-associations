#!/usr/bin/env python3
"""Paper 1 reproducible pipeline -- Stage 04: Feature Extraction.

Standalone port of pipeline/module2/{run_module2,features}.py's logic
(frozen reference; not imported at runtime -- see Paper1/MIGRATION_PLAN.md).
Extracts feature registry (time/frequency/entropy/complexity, FWH + Broad bands,
plus shared band-power ratios) for every window in the Stage 02 manifest, across all
four Paper 1 channels (EHG9, EHG10, EHG11, EHG12).

Updated to include extended Sample Entropy Profile statistics:
  - MedianSampEn
  - SDSampEn
  - KurtosisSampEn
  - SkewnessSampEn

Usage:
    python src/features.py
    python src/features.py --channels EHG9
    python src/features.py --force
    python src/features.py --max-windows 20    # fast smoke test
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import logging
import os
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from scipy.stats import kurtosis, skew

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
_SRC_DIR = Path(__file__).resolve().parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(1, str(_PROJECT_ROOT))

import config as project_config  # noqa: E402



def _load_sibling_module(name: str, filename: str):
    """Loads a module by file path and registers it in sys.modules first."""
    path = Path(__file__).resolve().parent / filename
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


preprocessing = _load_sibling_module("stage03_preprocessing", "preprocessing.py")

from amrita_biosignal_feature_engine import (  # noqa: E402
    BandPowerRequest,
    ExtractorConfig,
    FeatureExtractor,
    WelchPSDConfig,
    compute_psd,
)

import wfdb  # noqa: E402

LOGGER = logging.getLogger("paper1.stage04")


# =============================================================================
# Section 1 -- Memory-safe entropy adaptations
# =============================================================================

MEMORY_SAFE_ENTROPY_FEATURES = ("approximate_entropy", "fuzzy_entropy", "distribution_entropy")

# Registered new Sample Entropy profile scalar statistics
LOCAL_ENTROPY_FEATURES = (
    *MEMORY_SAFE_ENTROPY_FEATURES,
    "sampen_fixed_r",
    "totalsampen",
    "averagesampen",
    "mediansampen",      # <-- NEW
    "sdsampen",          # <-- NEW
    "kurtosissampen",    # <-- NEW
    "skewnesssampen",    # <-- NEW
)

ENGINE_FEATURES = (
    *project_config.TIME_FEATURES,
    *tuple(name for name in project_config.FREQUENCY_FEATURES if name != "spectral_spread"),
    *tuple(name for name in project_config.ENTROPY_FEATURES if name not in LOCAL_ENTROPY_FEATURES),
    *project_config.COMPLEXITY_FEATURES,
)


@dataclass(frozen=True)
class BandExtraction:
    values: dict
    band_power: float
    diagnostics: tuple


def _aligned_templates(data: np.ndarray, dimension: int, count: int) -> np.ndarray:
    return np.column_stack([data[offset:count + offset] for offset in range(dimension)])


def approximate_entropy_memory_safe(signal: Any, *, order: int = 2, tolerance_fraction: float = 0.2) -> float:
    data = preprocessing.validate_raw_signal(signal)
    if order < 1 or data.size < order + 2:
        return float("nan")
    standard_deviation = float(np.std(data, ddof=0))
    if standard_deviation == 0:
        return float("nan")
    tolerance = tolerance_fraction * standard_deviation

    def phi(dimension: int) -> float:
        count = data.size - dimension + 1
        templates = _aligned_templates(data, dimension, count)
        tree = cKDTree(templates)
        matches = tree.query_ball_point(templates, tolerance, p=np.inf, return_length=True, workers=1)
        probabilities = np.asarray(matches, dtype=np.float64) / count
        return float(np.mean(np.log(probabilities)))

    return float(phi(order) - phi(order + 1))


def _centered_templates(data: np.ndarray, dimension: int, template_count: int) -> np.ndarray:
    templates = _aligned_templates(data, dimension, template_count)
    return np.asarray(templates - np.mean(templates, axis=1, keepdims=True), dtype=np.float64)


def _pairwise_similarity_mean(templates: np.ndarray, tolerance: float, exponent: float, *, block_size: int = 512) -> float:
    count = templates.shape[0]
    similarity_sum = 0.0
    pair_count = 0
    scale = -np.log(2.0)
    for left_start in range(0, count, block_size):
        left = templates[left_start:left_start + block_size]
        for right_start in range(left_start, count, block_size):
            right = templates[right_start:right_start + block_size]
            distances = cdist(left, right, metric="chebyshev")
            with np.errstate(over="ignore"):
                similarities = np.exp(scale * np.power(distances / tolerance, exponent))
            if left_start == right_start:
                indices = np.triu_indices(similarities.shape[0], k=1)
                similarity_sum += float(np.sum(similarities[indices]))
                pair_count += int(indices[0].size)
            else:
                similarity_sum += float(np.sum(similarities))
                pair_count += int(similarities.size)
    return similarity_sum / pair_count if pair_count else float("nan")


def fuzzy_entropy_memory_safe(signal: Any, *, order: int = 2, delay: int = 1,
                                tolerance_fraction: float = 0.2, exponent: float = 2.0) -> float:
    data = preprocessing.validate_raw_signal(signal)
    if order < 1 or delay < 1 or data.size < order * delay + 2:
        return float("nan")
    if delay != 1:
        raise ValueError("Stage 04 memory-safe fuzzy entropy requires delay=1")
    standard_deviation = float(np.std(data, ddof=0))
    if standard_deviation == 0:
        return float("nan")
    tolerance = tolerance_fraction * standard_deviation
    template_count = data.size - order * delay
    phi_order = _pairwise_similarity_mean(_centered_templates(data, order, template_count), tolerance, exponent)
    phi_next = _pairwise_similarity_mean(_centered_templates(data, order + 1, template_count), tolerance, exponent)
    if phi_order == 0 or phi_next == 0 or not np.isfinite(phi_order + phi_next):
        return float("nan")
    return float(np.log(phi_order) - np.log(phi_next))


def distribution_entropy_memory_safe(signal: Any, *, order: int = 2, delay: int = 1,
                                       n_bins: int = 500, normalize: bool = True) -> float:
    data = preprocessing.validate_raw_signal(signal)
    if order < 1 or delay < 1 or n_bins < 2:
        raise ValueError("invalid distribution entropy parameters")
    template_count = data.size - (order - 1) * delay
    if template_count < 2:
        return float("nan")
    templates = np.column_stack([data[offset * delay: offset * delay + template_count] for offset in range(order)])
    tree = cKDTree(templates)
    nearest = tree.query(templates, k=2, p=np.inf, workers=1)[0][:, 1]
    minimum_distance = float(np.min(nearest))
    maximum_distance = float(np.max(np.ptp(templates, axis=0)))
    if minimum_distance == maximum_distance:
        return float("nan")

    edges = np.linspace(minimum_distance, maximum_distance, n_bins + 1)
    thresholds = np.nextafter(edges[1:], -np.inf)
    thresholds[-1] = edges[-1]
    ordered_shell_counts = np.asarray(
        tree.count_neighbors(tree, thresholds, p=np.inf, cumulative=False), dtype=np.int64,
    )
    ordered_shell_counts[0] -= template_count
    counts = ordered_shell_counts.astype(np.float64) / 2.0
    probabilities = counts[counts > 0] / np.sum(counts)
    result = -float(np.sum(probabilities * np.log(probabilities)))
    if normalize:
        result /= float(np.log(n_bins))
    return result


def sample_entropy_fixed_r(signal: Any, *, m: int = project_config.SAMPEN_EMBEDDING_DIMENSION,
                              r_fraction: float = project_config.SAMPEN_R_FRACTION) -> float:
    """Richman-Moorman SampEn with memory-safe KD-tree approach."""
    data = np.asarray(signal, dtype=np.float64).reshape(-1)
    n_samples = data.size
    if m < 1:
        raise ValueError("m must be positive")
    if not np.isfinite(r_fraction) or r_fraction <= 0:
        raise ValueError("r_fraction must be finite and positive")
    if n_samples < m + 3 or not np.all(np.isfinite(data)):
        return float("nan")
    standard_deviation = float(np.std(data, ddof=0))
    if standard_deviation == 0:
        return float("nan")
    tolerance = r_fraction * standard_deviation

    templates_m = np.lib.stride_tricks.sliding_window_view(data, m)
    templates_m1 = np.lib.stride_tricks.sliding_window_view(data, m + 1)
    aligned_count = templates_m1.shape[0]
    templates_m = np.asarray(templates_m[:aligned_count], dtype=np.float64)
    templates_m1 = np.asarray(templates_m1, dtype=np.float64)

    def ordered_nonself_matches(templates: np.ndarray) -> float:
        tree = cKDTree(templates)
        return float(tree.count_neighbors(tree, tolerance, p=np.inf) - aligned_count)

    denominator = float(aligned_count * (aligned_count - 1))
    b_probability = ordered_nonself_matches(templates_m) / denominator
    a_probability = ordered_nonself_matches(templates_m1) / denominator
    if b_probability <= 0 or a_probability <= 0:
        return float("nan")
    return float(-np.log(a_probability / b_probability))


def _rounded_distance_shell_counts(templates: np.ndarray, maximum_bin: int, *, decimals: int) -> np.ndarray:
    scale = float(10 ** decimals)
    bins = np.arange(maximum_bin + 1, dtype=np.float64)
    upper_boundaries = (bins + 0.5) / scale
    odd_bins = (np.arange(maximum_bin + 1) % 2) == 1
    upper_boundaries[odd_bins] = np.nextafter(upper_boundaries[odd_bins], -np.inf)
    tree = cKDTree(templates)
    counts = np.asarray(
        tree.count_neighbors(tree, upper_boundaries, p=np.inf, cumulative=False), dtype=np.int64,
    )
    counts[0] -= templates.shape[0]
    if np.any(counts < 0):
        raise RuntimeError("invalid sample-entropy profile pair counts")
    return counts


def sample_entropy_profile_summaries(
    signal: Any, *, order: int = project_config.SAMPEN_EMBEDDING_DIMENSION,
    decimals: int = project_config.SAMPEN_PROFILE_DISTANCE_DECIMALS,
    epsilon: float = project_config.SAMPEN_PROFILE_EPSILON,
) -> tuple[float, float, float, float, float, float]:
    """Computes TotalSampEn, AvgSampEn, MedianSampEn, SDSampEn, KurtosisSampEn, SkewnessSampEn."""
    data = preprocessing.validate_raw_signal(signal)
    
    # 6 NaN fallback tuple for guard clauses
    nan_tuple = (float("nan"), float("nan"), float("nan"), float("nan"), float("nan"), float("nan"))
    
    if order < 1 or data.size < order + 2:
        return nan_tuple

    template_count = data.size - order
    templates_m = _aligned_templates(data, order, template_count)
    templates_m1 = _aligned_templates(data, order + 1, template_count)

    maximum_distance = max(
        float(np.max(np.ptp(templates_m, axis=0))), float(np.max(np.ptp(templates_m1, axis=0))),
    )
    maximum_rounded = float(np.round(maximum_distance, decimals))
    maximum_bin = int(np.rint(maximum_rounded * (10 ** decimals)))

    shell_m = _rounded_distance_shell_counts(templates_m, maximum_bin, decimals=decimals)
    shell_m1 = _rounded_distance_shell_counts(templates_m1, maximum_bin, decimals=decimals)

    observed = (shell_m > 0) | (shell_m1 > 0)
    denominator = float(template_count * (template_count - 1))
    b = np.cumsum(shell_m, dtype=np.float64) / denominator
    a = np.cumsum(shell_m1, dtype=np.float64) / denominator

    profile = np.log((b + epsilon) / (a + epsilon))
    valid = observed & np.isfinite(profile) & (b > epsilon) & (a > epsilon)
    values = profile[valid]

    if values.size == 0:
        return nan_tuple

    # Calculate profile summary statistics
    total_sampen = float(np.sum(values))
    avg_sampen = float(np.mean(values))
    median_sampen = float(np.median(values))
    sd_sampen = float(np.std(values))
    kurtosis_sampen = float(kurtosis(values))
    skewness_sampen = float(skew(values))

    return total_sampen, avg_sampen, median_sampen, sd_sampen, kurtosis_sampen, skewness_sampen


# =============================================================================
# Section 2 -- Band feature extraction (vendored engine + local additions)
# =============================================================================

def _extractor() -> FeatureExtractor:
    psd_config = WelchPSDConfig(
        window_length=project_config.WELCH_WINDOW_LENGTH_S, overlap=project_config.WELCH_OVERLAP_S,
        detrend=project_config.WELCH_DETREND, scaling=project_config.WELCH_SCALING, window=project_config.WELCH_WINDOW,
    )
    return FeatureExtractor(ExtractorConfig(sampling_frequency=project_config.OUTPUT_SAMPLING_FREQUENCY_HZ, psd=psd_config))


def _spectral_spread(signal: np.ndarray) -> float:
    psd = compute_psd(
        signal, project_config.OUTPUT_SAMPLING_FREQUENCY_HZ,
        WelchPSDConfig(
            window_length=project_config.WELCH_WINDOW_LENGTH_S, overlap=project_config.WELCH_OVERLAP_S,
            detrend=project_config.WELCH_DETREND, scaling=project_config.WELCH_SCALING, window=project_config.WELCH_WINDOW,
        ),
    )
    frequencies, power = psd.frequencies, psd.values
    total = float(trapezoid(power, frequencies))
    if total <= 0:
        return float("nan")
    centroid = float(trapezoid(frequencies * power, frequencies) / total)
    variance = float(trapezoid(((frequencies - centroid) ** 2) * power, frequencies) / total)
    return float(np.sqrt(max(variance, 0.0)))


def extract_band_features(signal: Any, band_name: str) -> BandExtraction:
    """Extract scalar features including local memory-safe entropy metrics."""
    if band_name not in project_config.BANDS_HZ:
        raise ValueError(f"unknown band: {band_name!r}")
    data = preprocessing.validate_raw_signal(signal)
    low_hz, high_hz = project_config.BANDS_HZ[band_name]
    power_output_name = "_absolute_band_power"
    requests = (*ENGINE_FEATURES, BandPowerRequest(output_name=power_output_name, band=(low_hz, high_hz), relative=False))
    result = _extractor().extract(data, features=requests)

    diagnostics = tuple(
        f"{item.feature_name}:{item.severity.value}:{item.code.value}:{item.message}" for item in result.diagnostics
    )
    values = {name: float(result.values[name]) for name in ENGINE_FEATURES}

    # Spectral Spread
    try:
        values["spectral_spread"] = _spectral_spread(data)
    except Exception as exc:
        values["spectral_spread"] = float("nan")
        diagnostics = (*diagnostics, f"spectral_spread:error:{type(exc).__name__}:{exc}")

    # Fixed-r SampEn
    try:
        values["sampen_fixed_r"] = sample_entropy_fixed_r(data)
    except Exception as exc:
        values["sampen_fixed_r"] = float("nan")
        diagnostics = (*diagnostics, f"sampen_fixed_r:error:{type(exc).__name__}:{exc}")

    # Memory-safe entropy features
    for feature_name, function in (
        ("approximate_entropy", approximate_entropy_memory_safe),
        ("fuzzy_entropy", fuzzy_entropy_memory_safe),
        ("distribution_entropy", distribution_entropy_memory_safe),
    ):
        try:
            values[feature_name] = function(data)
        except Exception as exc:
            values[feature_name] = float("nan")
            diagnostics = (*diagnostics, f"{feature_name}:error:{type(exc).__name__}:{exc}")

    # SampEn Profile statistics
    try:
        (
            total_sampen,
            average_sampen,
            median_sampen,
            sd_sampen,
            kurtosis_sampen,
            skewness_sampen,
        ) = sample_entropy_profile_summaries(data)

        values["totalsampen"] = total_sampen
        values["averagesampen"] = average_sampen
        values["mediansampen"] = median_sampen
        values["sdsampen"] = sd_sampen
        values["kurtosissampen"] = kurtosis_sampen
        values["skewnesssampen"] = skewness_sampen
    except Exception as exc:
        values["totalsampen"] = float("nan")
        values["averagesampen"] = float("nan")
        values["mediansampen"] = float("nan")
        values["sdsampen"] = float("nan")
        values["kurtosissampen"] = float("nan")
        values["skewnesssampen"] = float("nan")
        diagnostics = (*diagnostics, f"sample_entropy_profile:error:{type(exc).__name__}:{exc}")

    return BandExtraction(values=values, band_power=float(result.values[power_output_name]), diagnostics=diagnostics)


def namespace_band_features(band_name: str, values: dict) -> dict:
    output = {}
    for domain, names in (
        ("time", project_config.TIME_FEATURES), ("freq", project_config.FREQUENCY_FEATURES),
        ("entropy", project_config.ENTROPY_FEATURES), ("complexity", project_config.COMPLEXITY_FEATURES),
    ):
        for name in names:
            output[f"{band_name}_{domain}_{name}"] = float(values.get(name, np.nan))
    return output


def _window_features(fwh_window: np.ndarray, broad_window: np.ndarray) -> tuple:
    if fwh_window.size != broad_window.size or fwh_window.size < 2:
        raise ValueError("filtered window arrays must be aligned and nonempty")
    fwh = extract_band_features(fwh_window, "fwh")
    broad = extract_band_features(broad_window, "broad")
    values = {
        **namespace_band_features("fwh", fwh.values), **namespace_band_features("broad", broad.values),
        "freq_band_power_fwh": fwh.band_power, "freq_band_power_broad": broad.band_power,
        "freq_power_ratio_fwh_over_broad": (
            fwh.band_power / broad.band_power if np.isfinite(broad.band_power) and broad.band_power > 0 else float("nan")
        ),
    }
    if tuple(values) != project_config.FEATURE_COLUMNS:
        raise RuntimeError("extracted feature order does not match the registry")
    diagnostics = tuple(f"fwh:{item}" for item in fwh.diagnostics) + tuple(f"broad:{item}" for item in broad.diagnostics)
    return values, diagnostics


def _nan_feature_row() -> dict:
    return {name: float("nan") for name in project_config.FEATURE_COLUMNS}


# =============================================================================
# Section 3 -- Orchestration (resumable, per (channel, recording_id))
# =============================================================================

def _configure_logging(results_dir: Path, suffix: str) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.setLevel(logging.INFO)
    LOGGER.handlers.clear()
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    file_handler = logging.FileHandler(results_dir / f"stage04_{suffix}.log", mode="a")
    file_handler.setFormatter(formatter)
    LOGGER.addHandler(stream_handler)
    LOGGER.addHandler(file_handler)


def _validate_manifest(manifest: pd.DataFrame, window_size_s: int) -> None:
    required = {"window_id", "recording_id", "window_size_s", "start_sample", "end_sample", "window_class", "patient_id"}
    missing = sorted(required - set(manifest.columns))
    if missing:
        raise ValueError(f"manifest is missing required columns: {missing}")
    if manifest.empty:
        raise ValueError("manifest contains no windows")
    if manifest["window_id"].duplicated().any():
        raise ValueError("manifest contains duplicate window_id values")
    if not (manifest["window_size_s"] == window_size_s).all():
        raise ValueError("manifest window_size_s does not match requested size")
    if manifest[["start_sample", "end_sample"]].isna().any().any():
        raise ValueError("manifest contains missing sample bounds")
    expected_raw_samples = int(round(window_size_s * project_config.INPUT_SAMPLING_FREQUENCY_HZ))
    observed = manifest["end_sample"].astype(int) - manifest["start_sample"].astype(int)
    if not (observed == expected_raw_samples).all():
        raise ValueError(f"manifest windows must each contain {expected_raw_samples} raw samples")


def _read_channel(recording_id: str, channel: str) -> tuple:
    record_path = project_config.RAW_DIR / recording_id
    record = wfdb.rdrecord(str(record_path), channel_names=[channel])
    signal = np.asarray(record.p_signal[:, 0], dtype=np.float64)
    sampling_frequency_hz = float(record.fs)
    if not np.isclose(sampling_frequency_hz, project_config.INPUT_SAMPLING_FREQUENCY_HZ, rtol=0, atol=1e-9):
        raise ValueError(
            f"{recording_id}: expected {project_config.INPUT_SAMPLING_FREQUENCY_HZ} Hz, observed {sampling_frequency_hz} Hz"
        )
    return signal, sampling_frequency_hz


def _process_recording_channel(recording_id: str, channel: str, recording_windows: pd.DataFrame,
                                 expected_filtered_samples: int) -> tuple:
    compact_rows, log_rows = [], []
    record_error = ""
    try:
        raw_signal, fs_hz = _read_channel(recording_id, channel)
        filtered = preprocessing.filter_recording(raw_signal, fs_hz)
    except Exception as exc:
        record_error = f"{type(exc).__name__}: {exc}"
        LOGGER.exception("%s/%s: recording load/filter failed", recording_id, channel)
        filtered = {}

    for _, window in recording_windows.iterrows():
        metadata = window.to_dict()
        base = {
            **metadata, "signal_channel": channel,
            "input_sampling_frequency_hz": project_config.INPUT_SAMPLING_FREQUENCY_HZ,
            "filtered_sampling_frequency_hz": project_config.OUTPUT_SAMPLING_FREQUENCY_HZ,
            "filtered_sample_count": expected_filtered_samples,
        }
        error = record_error
        values = _nan_feature_row()
        diagnostics: tuple = ()
        if not error:
            try:
                fwh_window = preprocessing.filtered_window(filtered["fwh"], int(window["start_sample"]), int(window["end_sample"]))
                broad_window = preprocessing.filtered_window(filtered["broad"], int(window["start_sample"]), int(window["end_sample"]))
                if fwh_window.size != expected_filtered_samples or broad_window.size != expected_filtered_samples:
                    raise ValueError(
                        f"expected {expected_filtered_samples} filtered samples, got FWH={fwh_window.size}, Broad={broad_window.size}"
                    )
                values, diagnostics = _window_features(fwh_window, broad_window)
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"
                LOGGER.exception("%s/%s: window %s extraction failed", recording_id, channel, window["window_id"])

        compact_rows.append({**base, **values})
        log_rows.append({
            "window_id": window["window_id"], "recording_id": recording_id, "channel": channel,
            "status": "error" if error else "ok", "error": error,
            "diagnostic_count": len(diagnostics), "diagnostics": json.dumps(diagnostics, separators=(",", ":")),
            "feature_nan_count": int(sum(pd.isna(values[name]) for name in project_config.FEATURE_COLUMNS)),
        })

    return compact_rows, log_rows, record_error


def run(
    window_size_s: int = project_config.WINDOW_SIZE_S,
    channels: tuple = project_config.CHANNELS,
    windows_dir: Path = None,
    results_dir: Path = None,
    max_windows: int | None = None,
    force: bool = False,
) -> int:
    windows_dir = windows_dir or project_config.WINDOWS_DIR
    results_dir = results_dir or (project_config.FEATURES_DIR)
    results_dir.mkdir(parents=True, exist_ok=True)

    suffix = f"{window_size_s}s"
    _configure_logging(results_dir, suffix)

    manifest_path = windows_dir / f"window_manifest_{window_size_s}s.csv"
    LOGGER.info("Loading manifest: %s", manifest_path)
    manifest = pd.read_csv(manifest_path)
    _validate_manifest(manifest, window_size_s)
    if max_windows is not None:
        manifest = manifest.groupby("recording_id", group_keys=False).apply(lambda g: g.head(max_windows))
    manifest = manifest.reset_index(drop=True)

    expected_filtered_samples = int(round(window_size_s * project_config.OUTPUT_SAMPLING_FREQUENCY_HZ))
    run_log_path = results_dir / f"run_log_{suffix}.jsonl"

    overall_status = 0
    for channel in channels:
        LOGGER.info("=== channel=%s ===", channel)
        partial_dir = results_dir / "_partial" / channel
        partial_dir.mkdir(parents=True, exist_ok=True)

        recordings = list(manifest.groupby("recording_id", sort=False))
        for i, (recording_id, recording_windows) in enumerate(recordings, 1):
            partial_path = partial_dir / f"{recording_id}.csv"
            if partial_path.is_file() and not force:
                LOGGER.info("%d/%d %s/%s: resumed (partial exists)", i, len(recordings), channel, recording_id)
                _append_run_log(run_log_path, channel, recording_id, "resumed_skip", 0, "")
                continue

            started = datetime.now(timezone.utc)
            compact_rows, log_rows, record_error = _process_recording_channel(
                recording_id, channel, recording_windows, expected_filtered_samples,
            )
            duration_ms = int((datetime.now(timezone.utc) - started).total_seconds() * 1000)

            if record_error:
                log_dir = results_dir / "_partial_logs" / channel
                log_dir.mkdir(parents=True, exist_ok=True)
                pd.DataFrame(log_rows).to_csv(log_dir / f"{recording_id}.FAILED.csv", index=False)
                _append_run_log(run_log_path, channel, recording_id, "recording_load_failed", duration_ms, record_error)
                LOGGER.error("%d/%d %s/%s: RECORDING-LEVEL FAILURE, will retry next run -- %s",
                             i, len(recordings), channel, recording_id, record_error)
                overall_status = 1
                continue

            pd.DataFrame(compact_rows).to_csv(partial_path, index=False)
            log_dir = results_dir / "_partial_logs" / channel
            log_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(log_rows).to_csv(log_dir / f"{recording_id}.csv", index=False)

            n_window_errors = sum(1 for r in log_rows if r["status"] == "error")
            status = "ok_with_window_errors" if n_window_errors else "ok"
            _append_run_log(run_log_path, channel, recording_id, status, duration_ms, "")
            LOGGER.info("%d/%d %s/%s: %s (%d windows, %d window-errors, %dms)",
                        i, len(recordings), channel, recording_id, status, len(compact_rows), n_window_errors, duration_ms)
            if n_window_errors:
                overall_status = 1

        missing_recordings = [rid for rid, _ in recordings if not (partial_dir / f"{rid}.csv").is_file()]
        if missing_recordings:
            LOGGER.warning(
                "channel=%s: %d/%d recording(s) still have no successful output: %s",
                channel, len(missing_recordings), len(recordings), missing_recordings,
            )
            overall_status = 1

        compact_frames = [pd.read_csv(partial_dir / f"{rid}.csv") for rid, _ in recordings if (partial_dir / f"{rid}.csv").is_file()]
        compact_df = pd.concat(compact_frames, ignore_index=True) if compact_frames else pd.DataFrame()
        log_dir = results_dir / "_partial_logs" / channel
        log_frames = [pd.read_csv(log_dir / f"{rid}.csv") for rid, _ in recordings if (log_dir / f"{rid}.csv").is_file()]
        log_df = pd.concat(log_frames, ignore_index=True) if log_frames else pd.DataFrame()

        feature_path = results_dir / f"window_features_{suffix}_{channel}.csv"
        log_path = results_dir / f"feature_processing_log_{suffix}_{channel}.csv"
        compact_df.to_csv(feature_path, index=False)
        log_df.to_csv(log_path, index=False)
        LOGGER.info(
            "channel=%s: wrote %s (%d rows from %d/%d recordings), %s",
            channel, feature_path.name, len(compact_df), len(recordings) - len(missing_recordings),
            len(recordings), log_path.name,
        )

    _write_run_metadata(results_dir, suffix, window_size_s, channels, manifest_path, max_windows)
    LOGGER.info("done. outputs -> %s", results_dir)
    return overall_status


def _append_run_log(path: Path, channel: str, recording_id: str, status: str, duration_ms: int, error: str) -> None:
    entry = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(), "channel": channel,
        "recording_id": recording_id, "status": status, "duration_ms": duration_ms, "error": error,
    }
    with open(path, "a") as f:
        f.write(json.dumps(entry) + "\n")


def _write_run_metadata(results_dir: Path, suffix: str, window_size_s: int, channels: tuple,
                          manifest_path: Path, max_windows) -> None:
    metadata = {
        "stage": "04_feature_extraction",
        "scope": {"window_size_s": window_size_s, "channels": list(channels), "max_windows": max_windows},
        "signal_processing": {
            "input_sampling_frequency_hz": project_config.INPUT_SAMPLING_FREQUENCY_HZ,
            "analysis_sampling_frequency_hz": project_config.ANALYSIS_SAMPLING_FREQUENCY_HZ,
            "decimation_factor": project_config.DECIMATION_FACTOR,
            "filter": "Butterworth SOS, zero-phase sosfiltfilt", "filter_order": project_config.FILTER_ORDER,
            "bands_hz": project_config.BANDS_HZ, "filter_scope": "full recording before manifest-window slicing",
        },
        "feature_count": len(project_config.FEATURE_COLUMNS),
        "feature_columns": list(project_config.FEATURE_COLUMNS),
        "input_manifest": str(manifest_path),
        "feature_engine_commit": project_config.ENGINE_COMMIT,
        "execution_timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    (results_dir / f"run_metadata_{suffix}.json").write_text(json.dumps(metadata, indent=2, default=str))


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--window-size", type=int, default=project_config.WINDOW_SIZE_S)
    parser.add_argument("--channels", type=str, default=None, help="Comma-separated, e.g. EHG9,EHG10")
    parser.add_argument("--windows-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--max-windows", type=int, default=None, help="Per-recording cap, for fast smoke tests")
    parser.add_argument("--force", action="store_true", help="Recompute even if a partial result exists")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    channels = tuple(c.strip() for c in args.channels.split(",")) if args.channels else project_config.CHANNELS
    sys.exit(run(
        window_size_s=args.window_size, channels=channels, windows_dir=args.windows_dir,
        results_dir=args.results_dir, max_windows=args.max_windows, force=args.force,
    ))
