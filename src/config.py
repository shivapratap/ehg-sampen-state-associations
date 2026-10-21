"""Shared configuration for the Paper 1 reproducible pipeline.

Standalone project: does not import anything from ../../pipeline/ at
runtime. The FROZEN_* paths below point at the old pipeline's already
-computed outputs and are used ONLY by the reproducibility-gate validation
checks (results/paper1/validation/), never by the numbered pipeline stages
themselves while they compute their own results.

Fixed protocol (see Paper1/MIGRATION_PLAN.md):
    channels    EHG9 (primary), EHG10, EHG11, EHG12
    window size 60 seconds
    band        FWH (0.3-1.0 Hz) primary; Broad (0.1-3.0 Hz) retained as
                secondary/sensitivity, same as the old pipeline
    classes     UC (uterine contraction), FM (fetal movement), IN (inactivity)
"""
from __future__ import annotations

import os
from pathlib import Path

# Repository paths. Raw WFDB data stay outside Git; all large derived data use work/.
PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT_ROOT / "data"
WORK_DIR = PROJECT_ROOT / "work"
RESULTS_DIR = PROJECT_ROOT / "results"
RAW_DIR_DEFAULT = DATA_DIR / "raw"
ANNOTATIONS_DIR = WORK_DIR / "annotations"
WINDOWS_DIR = WORK_DIR / "windows"
AUDIT_DIR = WORK_DIR / "audit"
WINDOWING_DIR = WORK_DIR / "windowing"
FEATURES_DIR = WORK_DIR / "features"
QC_DIR = WORK_DIR / "feature_qc"
UNITS_DIR = WORK_DIR / "primary"
PRIMARY_DIR = RESULTS_DIR / "primary"
ROBUSTNESS_DIR = RESULTS_DIR / "robustness"
ENGINE_COMMIT = "4dc5025b265f300d7e3375118f7794bfb5fca711"

def resolve_raw_dir() -> Path:
    """Raw WFDB directory: $EHG_RAW_DIR if set, else <repo>/data/raw.

    Not required to exist at import time, so the statistics stages can run from
    the derived tables without the raw database. Stages that read recordings
    check for it themselves (see data/README.md).
    """
    value = os.environ.get("EHG_RAW_DIR")
    return (Path(value).expanduser() if value else RAW_DIR_DEFAULT).resolve()

RAW_DIR = resolve_raw_dir()
RECORDS_FILE = RAW_DIR / "RECORDS"
SHA256SUMS_FILE = RAW_DIR / "SHA256SUMS.txt"

# --- Dataset expectations (ported from pipeline/module0/config.py) --------
# Comparison-only values; never assumed correct, always recomputed and
# reported against these for drift-detection (same design note as the old
# module0/config.py).
PUBLISHED_DOC_N_RECORDINGS = 122
PUBLISHED_DOC_N_ANNOTATED = 111
EXPECTED_N_PATIENTS = 45
EXPECTED_N_CHANNELS = 16
EXPECTED_FS_HZ = 200.0

# ice002_p_1of3 / ice039_p_6aof6 style stems; trailing-letter suffix is a
# genuine dataset naming pattern (mid-session split recordings), not a typo.
RECORD_STEM_RE = r"^(ice\d{3})_(p|l)_(\d+[a-z]?)of(\d+)$"

ANNOTATION_LABELS = {
    "C": "contraction",
    "(c)": "possible_contraction",
    "fm": "fetal_movement",
    "pm": "participant_movement",
    "em": "equipment_manipulation",
    "pos": "position_change",
}

QC_SEVERITY_ORDER = ("info", "warning", "error", "fatal")
RECOMMENDED_ACTIONS = ("retain", "inspect", "quarantine")
RECORDING_EXTENSIONS = (".hea", ".dat", ".jpg")

# --- Module 1 windowing constants (ported from pipeline/module1/config.py) -
GUARD_INTERVAL_S = 120.0
EDGE_TRIM_S = 15.0
TIGHT_CLEARANCE_THRESHOLD_S = 30.0
BOUNDARY_OVERLAP_THRESHOLD_S = 1.0
GUARD_ANNOTATION_TYPES = (
    "contraction", "possible_contraction", "fetal_movement",
    "participant_movement", "equipment_manipulation", "position_change",
)
PRIMARY_EVENT_TYPES = ("contraction", "possible_contraction", "fetal_movement")
SKIP_RECOMMENDED_ACTIONS = ("quarantine",)

# --- Paper 1 fixed protocol -------------------------------------------------
CHANNELS = ("EHG9", "EHG10", "EHG11", "EHG12")
PRIMARY_CHANNEL = "EHG9"
WINDOW_SIZE_S = 60
PRIMARY_BAND = "fwh"
SECONDARY_BAND = "broad"

# window_class (Module 1 vocabulary) -> Paper 1 protocol class_label.
# Anything not in this map is out of the primary UC/FM/IN protocol -- see
# EXCLUDED_WINDOW_CLASSES and src/01_label_ontology.py for the reasoning.
CLASS_MAP = {
    "contraction": "UC",
    "fetal_movement": "FM",
    "baseline": "IN",
}
EXCLUDED_WINDOW_CLASSES = ("possible_contraction", "mixed")

# ABFE is installed from the exact Git revision in requirements.txt.
# --- Stage 03/04 signal-processing constants (ported from module2/config.py)
INPUT_SAMPLING_FREQUENCY_HZ = 200.0
ANALYSIS_SAMPLING_FREQUENCY_HZ = 20.0
OUTPUT_SAMPLING_FREQUENCY_HZ = ANALYSIS_SAMPLING_FREQUENCY_HZ
DECIMATION_FACTOR = int(INPUT_SAMPLING_FREQUENCY_HZ / ANALYSIS_SAMPLING_FREQUENCY_HZ)
if INPUT_SAMPLING_FREQUENCY_HZ / ANALYSIS_SAMPLING_FREQUENCY_HZ != DECIMATION_FACTOR:
    raise RuntimeError("Stage 03/04 requires an integer decimation factor")
FILTER_ORDER = 4

BANDS_HZ = {
    "fwh": (0.3, 1.0),
    "broad": (0.1, 3.0),
}

# A 10-second Welch segment gives 0.1 Hz bin spacing, resolving the lower
# Broad-band edge and retaining multiple segments even for 40-second windows.
WELCH_WINDOW_LENGTH_S = 10.0
WELCH_OVERLAP_S = 5.0
WELCH_WINDOW = "hann"
WELCH_DETREND = "constant"
WELCH_SCALING = "density"

SAMPEN_EMBEDDING_DIMENSION = 2
SAMPEN_R_FRACTION = 0.2
SAMPEN_PROFILE_DISTANCE_DECIMALS = 3
SAMPEN_PROFILE_EPSILON = 1e-12

# Names exposed by the vendored engine, excluding the request-only largest
# Lyapunov exponent. Spectral spread and fixed-r SampEn are added locally
# (memory-safe adaptations of core/sampen_fixed_r.py and
# core/sampleEntropy_Gayathri.py -- see src/04_feature_extraction.py).
TIME_FEATURES = (
    "minimum", "maximum", "sum_value", "mean", "median", "standard_deviation",
    "variance", "kurtosis", "skewness", "mean_absolute_value", "root_mean_square",
    "peak_to_peak", "integrated_absolute_value", "waveform_length",
    "zero_crossing_count", "slope_sign_change_count",
)
FREQUENCY_FEATURES = (
    "peak_frequency", "mean_frequency", "median_frequency",
    "spectral_edge_frequency_95", "spectral_entropy", "spectral_spread",
)
ENTROPY_FEATURES = (
    "approximate_entropy", "permutation_entropy", "fuzzy_entropy",
    "distribution_entropy", "svd_entropy", "sampen_fixed_r",
    "totalsampen", "averagesampen", "mediansampen", "sdsampen",
    "kurtosissampen", "skewnesssampen",
)
COMPLEXITY_FEATURES = (
    "lempel_ziv_complexity", "hjorth_mobility", "hjorth_complexity",
    "fisher_information", "petrosian_fractal_dimension", "katz_fractal_dimension",
    "higuchi_fractal_dimension", "detrended_fluctuation_analysis",
)
SHARED_POWER_FEATURES = (
    "freq_band_power_fwh", "freq_band_power_broad", "freq_power_ratio_fwh_over_broad",
)


def _band_feature_names(band: str) -> tuple:
    return (
        *(f"{band}_time_{name}" for name in TIME_FEATURES),
        *(f"{band}_freq_{name}" for name in FREQUENCY_FEATURES),
        *(f"{band}_entropy_{name}" for name in ENTROPY_FEATURES),
        *(f"{band}_complexity_{name}" for name in COMPLEXITY_FEATURES),
    )


FEATURE_COLUMNS = (
    *_band_feature_names("fwh"),
    *_band_feature_names("broad"),
    *SHARED_POWER_FEATURES,
)

if len(FEATURE_COLUMNS) != 87 or len(set(FEATURE_COLUMNS)) != 87:
    raise RuntimeError("Stage 04 feature registry must contain 87 unique columns")


# Features with documented cross-platform floating-point sensitivity.
NUMERIC_SENSITIVITY_FEATURES = (
    "fwh_time_median", "fwh_time_skewness", "fwh_time_kurtosis", "fwh_time_sum_value",
    "broad_time_median", "broad_time_skewness", "broad_time_kurtosis", "broad_time_sum_value",
    "fwh_entropy_distribution_entropy", "fwh_entropy_totalsampen", "fwh_entropy_averagesampen",
    "broad_entropy_distribution_entropy", "broad_entropy_totalsampen", "broad_entropy_averagesampen",
    "broad_freq_median_frequency",
)
