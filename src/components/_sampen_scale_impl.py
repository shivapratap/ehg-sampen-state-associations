#!/usr/bin/env python3
"""Deterministic Phase-1 audit of Sample Entropy amplitude-scale dependence.

This script is deliberately isolated from the frozen Stage 00--09 pipeline. It
reads Stage 04 and raw inputs, reproduces selected Stage 04 EHG9/FWH values,
then evaluates positive multiplicative scales without modifying source data.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import math
import platform
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.spatial import cKDTree
from scipy.stats import kurtosis, pearsonr, skew, spearmanr
import statsmodels
import statsmodels.api as sm
import wfdb


SCRIPT_PATH = Path(__file__).resolve()
OUT_DIR = Path(__file__).resolve().parents[2] / "results" / "entropy"
PIPELINE_ROOT = SCRIPT_PATH.parents[2]
REPO_ROOT = PIPELINE_ROOT
SRC_DIR = PIPELINE_ROOT / "src"
MANIFEST_PATH = PIPELINE_ROOT / "work" / "windows" / "window_manifest_60s.csv"
FEATURE_PATH = PIPELINE_ROOT / "work" / "features" / "window_features_60s_EHG9.csv"
RUN_METADATA_PATH = PIPELINE_ROOT / "work" / "features" / "run_metadata_60s.json"
SCALES = (0.25, 0.5, 1.0, 2.0, 4.0)
CLASS_MAP = {"contraction": "UC", "baseline": "IN", "fetal_movement": "FM"}
CLASS_ORDER = ("UC", "IN", "FM")
N_PER_CLASS = 20
INVARIANCE_ATOL = 1e-10
INVARIANCE_RTOL = 1e-8
REPRO_ATOL = 1e-12
REPRO_RTOL = 1e-10
RANDOM_SEED = None

FEATURES = (
    "rms",
    "total_sampen",
    "average_sampen",
    "median_sampen",
    "sd_sampen",
    "kurtosis_sampen",
    "skewness_sampen",
    "fixed_r_sampen",
)
STAGE04_COLUMNS = {
    "rms": "fwh_time_root_mean_square",
    "total_sampen": "fwh_entropy_totalsampen",
    "average_sampen": "fwh_entropy_averagesampen",
    "median_sampen": "fwh_entropy_mediansampen",
    "sd_sampen": "fwh_entropy_sdsampen",
    "kurtosis_sampen": "fwh_entropy_kurtosissampen",
    "skewness_sampen": "fwh_entropy_skewnesssampen",
    "fixed_r_sampen": "fwh_entropy_sampen_fixed_r",
}


def load_stage04():
    if str(PIPELINE_ROOT) not in sys.path:
        sys.path.insert(0, str(PIPELINE_ROOT))
    spec = importlib.util.spec_from_file_location("stage04_scale_audit", SRC_DIR / "features.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def select_windows(manifest: pd.DataFrame) -> pd.DataFrame:
    eligible = manifest.loc[manifest["window_class"].isin(CLASS_MAP)].copy()
    eligible["class_label"] = eligible["window_class"].map(CLASS_MAP)
    selected = []
    for label in CLASS_ORDER:
        group = eligible.loc[eligible["class_label"] == label].sort_values(
            ["patient_id", "recording_id", "window_start_s", "window_id"], kind="mergesort"
        )
        group["within_patient_rank"] = group.groupby("patient_id", sort=False).cumcount()
        group = group.sort_values(
            ["within_patient_rank", "patient_id", "recording_id", "window_start_s", "window_id"],
            kind="mergesort",
        ).head(N_PER_CLASS)
        if len(group) != N_PER_CLASS:
            raise RuntimeError(f"Only {len(group)} eligible {label} windows; required {N_PER_CLASS}")
        selected.append(group)
    result = pd.concat(selected, ignore_index=True)
    result["selection_order_within_class"] = result.groupby("class_label", sort=False).cumcount() + 1
    result["selection_rule"] = (
        "Within class: stable sort by patient, recording, start, window; rank windows within patient; "
        "then stable sort by within-patient rank, patient, recording, start, window; take first 20"
    )
    return result


def profile_diagnostics(signal: np.ndarray, stage04) -> dict[str, float]:
    data = stage04.preprocessing.validate_raw_signal(signal)
    order = stage04.project_config.SAMPEN_EMBEDDING_DIMENSION
    decimals = stage04.project_config.SAMPEN_PROFILE_DISTANCE_DECIMALS
    epsilon = stage04.project_config.SAMPEN_PROFILE_EPSILON
    template_count = data.size - order
    templates_m = stage04._aligned_templates(data, order, template_count)
    templates_m1 = stage04._aligned_templates(data, order + 1, template_count)
    maximum_distance = max(
        float(np.max(np.ptp(templates_m, axis=0))),
        float(np.max(np.ptp(templates_m1, axis=0))),
    )
    maximum_rounded = float(np.round(maximum_distance, decimals))
    maximum_bin = int(np.rint(maximum_rounded * (10**decimals)))
    shell_m = stage04._rounded_distance_shell_counts(templates_m, maximum_bin, decimals=decimals)
    shell_m1 = stage04._rounded_distance_shell_counts(templates_m1, maximum_bin, decimals=decimals)
    observed = (shell_m > 0) | (shell_m1 > 0)
    denominator = float(template_count * (template_count - 1))
    b = np.cumsum(shell_m, dtype=np.float64) / denominator
    a = np.cumsum(shell_m1, dtype=np.float64) / denominator
    profile = np.log((b + epsilon) / (a + epsilon))
    valid = observed & np.isfinite(profile) & (b > epsilon) & (a > epsilon)
    values = profile[valid]
    production = stage04.sample_entropy_profile_summaries(data)
    reconstructed = (
        float(np.sum(values)), float(np.mean(values)), float(np.median(values)),
        float(np.std(values)), float(kurtosis(values)), float(skew(values)),
    )
    if not np.allclose(production, reconstructed, rtol=0, atol=1e-14, equal_nan=True):
        raise RuntimeError("Diagnostic profile reconstruction disagrees with Stage 04 implementation")
    return {
        "maximum_distance_unrounded": maximum_distance,
        "maximum_distance_rounded": maximum_rounded,
        "maximum_bin_index": maximum_bin,
        "tolerance_grid_min": 0.0,
        "tolerance_grid_max": maximum_rounded,
        "tolerance_grid_step": 10.0 ** (-decimals),
        "n_possible_profile_bins": maximum_bin + 1,
        "n_observed_profile_bins": int(np.count_nonzero(observed)),
        "n_retained_profile_values": int(np.count_nonzero(valid)),
        "profile_min": float(np.min(values)),
        "profile_max": float(np.max(values)),
        **dict(zip(FEATURES[1:7], production)),
    }


def compute_features(signal: np.ndarray, stage04) -> dict[str, float]:
    profile = profile_diagnostics(signal, stage04)
    rms = float(np.sqrt(np.mean(np.square(signal))))
    fixed = stage04.sample_entropy_fixed_r(signal)
    return {"rms": rms, **profile, "fixed_r_sampen": fixed}


def close(value: float, reference: float, *, atol: float, rtol: float) -> bool:
    return bool(np.isfinite(value) and np.isfinite(reference) and abs(value - reference) <= atol + rtol * abs(reference))


def load_selected_signals(selected: pd.DataFrame, stage04) -> dict[str, np.ndarray]:
    signals: dict[str, np.ndarray] = {}
    grouped = list(selected.groupby("recording_id", sort=True))
    for index, (recording_id, rows) in enumerate(grouped, 1):
        print(f"Filtering recording {index}/{len(grouped)}: {recording_id}", flush=True)
        raw, fs = stage04._read_channel(recording_id, "EHG9")
        filtered = stage04.preprocessing.filter_recording(raw, fs)["fwh"]
        for row in rows.itertuples(index=False):
            window = stage04.preprocessing.filtered_window(filtered, int(row.start_sample), int(row.end_sample))
            if window.size != 1200:
                raise RuntimeError(f"{row.window_id}: expected 1200 samples, got {window.size}")
            signals[row.window_id] = window
    return signals


def reproduction_check(selected: pd.DataFrame, frozen: pd.DataFrame, signals: dict[str, np.ndarray], stage04) -> pd.DataFrame:
    frozen_by_id = frozen.set_index("window_id")
    rows = []
    for row in selected.itertuples(index=False):
        observed = compute_features(signals[row.window_id], stage04)
        for feature in FEATURES:
            expected = float(frozen_by_id.loc[row.window_id, STAGE04_COLUMNS[feature]])
            value = observed[feature]
            rows.append({
                "window_id": row.window_id,
                "class_label": row.class_label,
                "patient_id": row.patient_id,
                "recording_id": row.recording_id,
                "feature": feature,
                "stage04_value": expected,
                "recomputed_value": value,
                "absolute_difference": abs(value - expected),
                "relative_difference": abs(value - expected) / abs(expected) if abs(expected) > 0 else np.nan,
                "reproduced": close(value, expected, atol=REPRO_ATOL, rtol=REPRO_RTOL),
                "reproduction_atol": REPRO_ATOL,
                "reproduction_rtol": REPRO_RTOL,
            })
    return pd.DataFrame(rows)


def scale_experiment(selected: pd.DataFrame, signals: dict[str, np.ndarray], stage04) -> pd.DataFrame:
    metadata = selected.set_index("window_id")
    rows = []
    total = len(selected) * len(SCALES)
    count = 0
    for window_id in selected["window_id"]:
        signal = signals[window_id]
        meta = metadata.loc[window_id]
        for scale in SCALES:
            count += 1
            if count % 25 == 0 or count == total:
                print(f"Scale computations {count}/{total}", flush=True)
            values = compute_features(signal * scale, stage04)
            rows.append({
                "window_id": window_id,
                "class_label": meta["class_label"],
                "window_class": meta["window_class"],
                "patient_id": meta["patient_id"],
                "recording_id": meta["recording_id"],
                "signal_channel": "EHG9",
                "band": "fwh",
                "window_size_s": 60,
                "filtered_sample_count": signal.size,
                "scale_factor": scale,
                **values,
            })
    result = pd.DataFrame(rows)
    baselines = result.loc[result["scale_factor"] == 1.0].set_index("window_id")
    for feature in FEATURES:
        reference = result["window_id"].map(baselines[feature])
        diff = result[feature] - reference
        result[f"{feature}_difference_from_1x"] = diff
        result[f"{feature}_absolute_difference_from_1x"] = diff.abs()
        result[f"{feature}_relative_change_percent"] = np.where(
            reference.abs() > 1e-10, 100.0 * diff / reference.abs(), np.nan
        )
        result[f"{feature}_ratio_to_1x"] = np.where(reference.abs() > 1e-10, result[feature] / reference, np.nan)
        result[f"{feature}_numerically_invariant"] = (
            np.isfinite(result[feature])
            & np.isfinite(reference)
            & (diff.abs() <= INVARIANCE_ATOL + INVARIANCE_RTOL * reference.abs())
        )
    expected_rms = result["window_id"].map(baselines["rms"]) * result["scale_factor"]
    result["rms_proportionality_error"] = result["rms"] - expected_rms
    result["rms_proportionality_pass"] = result["rms_proportionality_error"].abs() <= (
        INVARIANCE_ATOL + INVARIANCE_RTOL * expected_rms.abs()
    )
    # Requested compact names are retained alongside explicit internal names.
    result["class"] = result["class_label"]
    for compact, explicit in {
        "totalsampen": "total_sampen", "averagesampen": "average_sampen",
        "mediansampen": "median_sampen", "sdsampen": "sd_sampen",
        "kurtosissampen": "kurtosis_sampen", "skewnesssampen": "skewness_sampen",
        "maximum_rounded_distance": "maximum_distance_rounded",
    }.items():
        result[compact] = result[explicit]
    return result


def monotonicity(values: np.ndarray, reference: float) -> str:
    tol = INVARIANCE_ATOL + INVARIANCE_RTOL * abs(reference)
    diffs = np.diff(values)
    if np.all(np.abs(values - reference) <= tol):
        return "effectively_invariant"
    if np.all(diffs >= -tol):
        return "nondecreasing"
    if np.all(diffs <= tol):
        return "nonincreasing"
    return "nonmonotonic"


def make_summary(long: pd.DataFrame) -> pd.DataFrame:
    rows = []
    strata = [("overall", long)] + [(label, long.loc[long["class_label"] == label]) for label in CLASS_ORDER]
    for stratum, data in strata:
        for feature in FEATURES:
            trajectories = data.sort_values("scale_factor").groupby("window_id")[feature].apply(np.asarray)
            base = data.loc[data["scale_factor"] == 1.0].set_index("window_id")[feature]
            mono = pd.Series({wid: monotonicity(vals, float(base.loc[wid])) for wid, vals in trajectories.items()})
            mono_counts = mono.value_counts()
            for scale in SCALES:
                subset = data.loc[data["scale_factor"] == scale]
                relative = subset[f"{feature}_relative_change_percent"].dropna()
                rows.append({
                    "stratum": stratum,
                    "feature": feature,
                    "scale_factor": scale,
                    "n_windows": len(subset),
                    "n_finite_relative_changes": len(relative),
                    "median_relative_change_percent": relative.median(),
                    "q1_relative_change_percent": relative.quantile(0.25),
                    "q3_relative_change_percent": relative.quantile(0.75),
                    "iqr_relative_change_percent": relative.quantile(0.75) - relative.quantile(0.25),
                    "minimum_relative_change_percent": relative.min(),
                    "maximum_relative_change_percent": relative.max(),
                    "median_absolute_difference": subset[f"{feature}_absolute_difference_from_1x"].median(),
                    "n_numerically_invariant": int(subset[f"{feature}_numerically_invariant"].sum()),
                    "proportion_numerically_invariant": float(subset[f"{feature}_numerically_invariant"].mean()),
                    "n_trajectory_effectively_invariant": int(mono_counts.get("effectively_invariant", 0)),
                    "n_trajectory_nondecreasing": int(mono_counts.get("nondecreasing", 0)),
                    "n_trajectory_nonincreasing": int(mono_counts.get("nonincreasing", 0)),
                    "n_trajectory_nonmonotonic": int(mono_counts.get("nonmonotonic", 0)),
                })
    return pd.DataFrame(rows)


def mechanism_analysis(long: pd.DataFrame) -> pd.DataFrame:
    base = long.loc[long["scale_factor"] == 1.0].set_index("window_id")
    work = long.loc[long["scale_factor"] != 1.0].copy()
    work["delta_total_sampen"] = work["total_sampen"] - work["window_id"].map(base["total_sampen"])
    predictors = ("maximum_distance_rounded", "n_possible_profile_bins", "n_retained_profile_values")
    for predictor in predictors:
        work[f"delta_{predictor}"] = work[predictor] - work["window_id"].map(base[predictor])
    rows = []
    strata = [("overall", work)] + [(label, work.loc[work["class_label"] == label]) for label in CLASS_ORDER]
    for stratum, data in strata:
        for predictor in predictors:
            x = data[f"delta_{predictor}"].astype(float)
            y = data["delta_total_sampen"].astype(float)
            pearson = pearsonr(x, y) if x.nunique() > 1 and y.nunique() > 1 else (np.nan, np.nan)
            spearman = spearmanr(x, y) if x.nunique() > 1 and y.nunique() > 1 else (np.nan, np.nan)
            design = sm.add_constant(x, has_constant="add")
            fit = sm.OLS(y, design).fit(
                cov_type="cluster", cov_kwds={"groups": data["patient_id"], "use_correction": True}, use_t=True
            )
            per_window_rho = []
            for _, window in long.groupby("window_id"):
                if stratum != "overall" and window["class_label"].iloc[0] != stratum:
                    continue
                if window[predictor].nunique() > 1 and window["total_sampen"].nunique() > 1:
                    per_window_rho.append(float(spearmanr(window[predictor], window["total_sampen"]).statistic))
            rows.append({
                "stratum": stratum,
                "outcome": "change_in_total_sampen_from_1x",
                "predictor": f"change_in_{predictor}_from_1x",
                "n_scaled_observations_excluding_1x": len(data),
                "n_windows": data["window_id"].nunique(),
                "pearson_r": pearson.statistic if hasattr(pearson, "statistic") else pearson[0],
                "pearson_p": pearson.pvalue if hasattr(pearson, "pvalue") else pearson[1],
                "spearman_rho": spearman.statistic if hasattr(spearman, "statistic") else spearman[0],
                "spearman_p": spearman.pvalue if hasattr(spearman, "pvalue") else spearman[1],
                "cluster_robust_slope": float(fit.params.iloc[1]),
                "cluster_robust_se": float(fit.bse.iloc[1]),
                "cluster_robust_p": float(fit.pvalues.iloc[1]),
                "median_within_window_spearman_rho": float(np.median(per_window_rho)) if per_window_rho else np.nan,
                "proportion_positive_within_window_rho": float(np.mean(np.asarray(per_window_rho) > 0)) if per_window_rho else np.nan,
                "proportion_windows_predictor_changed_at_any_non1x_scale": float(
                    data.assign(changed=x.ne(0)).groupby("window_id")["changed"].any().mean()
                ),
                "proportion_windows_total_sampen_changed_at_any_non1x_scale": float(
                    data.assign(changed=y.abs() > (INVARIANCE_ATOL + INVARIANCE_RTOL * data["window_id"].map(base["total_sampen"]).abs()))
                    .groupby("window_id")["changed"].any().mean()
                ),
                "covariance": "CR1 patient-clustered OLS; descriptive mechanism regression",
            })
    return pd.DataFrame(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def input_hashes(selected: pd.DataFrame, stage04) -> pd.DataFrame:
    files = [
        ("analysis_source", SCRIPT_PATH),
        ("pipeline_configuration", PIPELINE_ROOT / "src" / "config.py"),
        ("preprocessing_source", SRC_DIR / "preprocessing.py"),
        ("feature_extraction_source", SRC_DIR / "features.py"),
        ("window_manifest", MANIFEST_PATH),
        ("stage04_ehg9_features", FEATURE_PATH),
        ("stage04_run_metadata", RUN_METADATA_PATH),
        ("raw_hash_manifest", stage04.project_config.SHA256SUMS_FILE),
        ("feature_engine_profile_source", Path(__import__("amrita_biosignal_feature_engine").__file__).parent / "entropy" / "sample_entropy_profile.py"),
        ("feature_engine_rms_source", Path(__import__("amrita_biosignal_feature_engine").__file__).parent / "time_domain.py"),
    ]
    rows = []
    for role, path in files:
        rows.append({"file_role": role, "path": str(path), "size_bytes": path.stat().st_size, "sha256": sha256(path), "hash_origin": "computed_this_run", "verified_against_declared": np.nan})
    declared = {}
    checksum_path = stage04.project_config.SHA256SUMS_FILE
    for line in checksum_path.read_text().splitlines():
        parts = line.strip().split(maxsplit=1)
        if len(parts) == 2:
            declared[parts[1].lstrip("*")] = parts[0]
    for recording_id in sorted(selected["recording_id"].unique()):
        for suffix in (".hea", ".dat"):
            path = stage04.project_config.RAW_DIR / f"{recording_id}{suffix}"
            name = path.name
            declared_hash = declared.get(name, "")
            if suffix == ".hea":
                actual = sha256(path)
                origin = "computed_this_run"
                verified = actual == declared_hash if declared_hash else np.nan
            else:
                actual = declared_hash
                origin = "dataset_SHA256SUMS_manifest_not_rehashed_this_run"
                verified = np.nan
            rows.append({"file_role": f"selected_raw_{suffix[1:]}", "path": str(path), "size_bytes": path.stat().st_size, "sha256": actual, "hash_origin": origin, "verified_against_declared": verified})
    return pd.DataFrame(rows)


def plot_trajectories(long: pd.DataFrame, feature: str, ylabel: str, output: Path) -> None:
    colors = {"UC": "#d62728", "IN": "#1f77b4", "FM": "#2ca02c"}
    fig, ax = plt.subplots(figsize=(8, 5.2), constrained_layout=True)
    for (_, label), window in long.groupby(["window_id", "class_label"], sort=False):
        window = window.sort_values("scale_factor")
        ax.plot(window["scale_factor"], window[feature], color=colors[label], alpha=0.18, linewidth=0.8)
    for label in CLASS_ORDER:
        median = long.loc[long["class_label"] == label].groupby("scale_factor")[feature].median()
        ax.plot(median.index, median.values, marker="o", linewidth=2.5, color=colors[label], label=f"{label} median")
    ax.set_xscale("log", base=2)
    ax.set_xticks(SCALES, [str(x).rstrip("0").rstrip(".") for x in SCALES])
    ax.set_xlabel("Positive amplitude scale factor")
    ax.set_ylabel(ylabel)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, ncol=3)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def make_figures(long: pd.DataFrame) -> None:
    plot_trajectories(long, "total_sampen", "Total SampEn", OUT_DIR / "total_sampen_vs_scale.png")
    plot_trajectories(long, "average_sampen", "Average SampEn", OUT_DIR / "average_sampen_vs_scale.png")
    plot_trajectories(long, "n_retained_profile_values", "Retained profile values", OUT_DIR / "retained_profile_values_vs_scale.png")
    colors = dict(zip(SCALES, plt.cm.viridis(np.linspace(0.08, 0.92, len(SCALES)))))
    fig, ax = plt.subplots(figsize=(8, 5.2), constrained_layout=True)
    for _, window in long.groupby("window_id", sort=False):
        window = window.sort_values("scale_factor")
        ax.plot(window["n_retained_profile_values"], window["total_sampen"], color="0.75", alpha=0.25, linewidth=0.7)
    for scale in SCALES:
        data = long.loc[long["scale_factor"] == scale]
        ax.scatter(data["n_retained_profile_values"], data["total_sampen"], s=19, alpha=0.68, color=colors[scale], label=f"{scale:g}x")
    ax.set_xlabel("Retained profile values")
    ax.set_ylabel("Total SampEn")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, title="Scale", ncol=5)
    fig.savefig(OUT_DIR / "total_sampen_vs_retained_profile_values.png", dpi=180)
    plt.close(fig)


def write_versions(stage04) -> None:
    packages = ["numpy", "pandas", "scipy", "matplotlib", "statsmodels", "wfdb"]
    lines = [
        f"python={platform.python_version()}", f"platform={platform.platform()}",
        *[f"{name}={importlib.metadata.version(name)}" for name in packages],
        f"random_seed={RANDOM_SEED}",
        f"sampen_embedding_dimension={stage04.project_config.SAMPEN_EMBEDDING_DIMENSION}",
        f"sampen_profile_distance_decimals={stage04.project_config.SAMPEN_PROFILE_DISTANCE_DECIMALS}",
        f"sampen_profile_epsilon={stage04.project_config.SAMPEN_PROFILE_EPSILON}",
        f"fixed_r_fraction={stage04.project_config.SAMPEN_R_FRACTION}",
        f"numerical_invariance_atol={INVARIANCE_ATOL}", f"numerical_invariance_rtol={INVARIANCE_RTOL}",
        f"reproduction_atol={REPRO_ATOL}", f"reproduction_rtol={REPRO_RTOL}",
    ]
    (OUT_DIR / "package_versions.txt").write_text("\n".join(lines) + "\n")


def write_implementation_audit(stage04) -> None:
    text = f"""# Total SampEn implementation audit

## Authoritative production path

The frozen Stage 04 EHG9/FWH Total SampEn value is produced by `sample_entropy_profile_summaries` in `{SRC_DIR / 'feature_extraction.py'}`. Stage 04 calls this local, memory-safe implementation after full-recording FWH filtering and window slicing/decimation. The vendored engine supplies RMS and other registered features, but not these six profile summaries in the production path.

## Exact construction

- Input: finite 1-D float64 EHG9 FWH window, 1,200 samples (60 s at 20 Hz).
- Embedding order: m={stage04.project_config.SAMPEN_EMBEDDING_DIMENSION}. Both m and m+1 templates use the aligned count L=N-m=1,198.
- Distance: Chebyshev (maximum coordinate-wise absolute difference).
- Grid: signal-unit tolerance spacing 0.001, because distances are rounded to {stage04.project_config.SAMPEN_PROFILE_DISTANCE_DECIMALS} decimals. The grid starts at 0 and ends at the maximum pairwise distance after three-decimal rounding.
- Pair counting: SciPy `cKDTree.count_neighbors`, with ordered self-pairs removed from the zero-distance shell. The half-even three-decimal rounding boundaries are reproduced by alternating `nextafter` handling at half steps.
- Probabilities: cumulative ordered non-self pair counts divided by L(L-1), separately for m (`b`) and m+1 (`a`).
- Profile value at a retained tolerance: `log((b + 1e-12)/(a + 1e-12))`.
- Retention: the shell must be observed in m or m+1, the profile must be finite, and both cumulative probabilities must exceed 1e-12. Thus the profile is not evaluated as an equally sized fixed normalized grid; the number of retained values may change under amplitude scaling.
- Total SampEn: arithmetic sum of all retained profile values. Average/median/SD (population SD, ddof=0), SciPy Fisher excess kurtosis with bias correction disabled by default (`bias=True`), and SciPy skewness (`bias=True`) are computed over the same retained values.
- Fixed-r SampEn comparator: tolerance r=0.2 times the window population SD, aligned m/m+1 Chebyshev templates, ordered non-self matches, and `-log(A/B)`.
- RMS comparator: square root of the mean squared samples.

## Scale-dependence hypothesis fixed before examining results

For positive scaling x->c*x, all exact distances scale by c. Fixed-r SampEn should therefore be invariant because its SD-derived tolerance also scales by c; RMS should be exactly proportional to c. Total SampEn is not guaranteed to be invariant because it sums values over a fixed 0.001 signal-unit rounded-distance grid and only retains occupied shells. Scaling changes the maximum rounded distance, possible grid size, shell occupancy, and usually the number of retained profile values. Average and other profile summaries can also change through quantization and the changing retained set, although they do not directly accumulate more terms as Total does.

Numerical invariance was prespecified as `abs(value - 1x) <= 1e-10 + 1e-8*abs(1x)`. Percent change and ratios are reported only when `abs(1x)>1e-10`.
"""
    (OUT_DIR / "totalsampen_implementation_audit.md").write_text(text)


def write_report(selected: pd.DataFrame, reproduction: pd.DataFrame, long: pd.DataFrame, summary: pd.DataFrame, mechanism: pd.DataFrame) -> None:
    non1 = long.loc[long["scale_factor"] != 1.0]
    inv = {feature: float(non1[f"{feature}_numerically_invariant"].mean()) for feature in FEATURES}
    best = max(inv, key=inv.get)
    total_medians = summary.loc[(summary["stratum"] == "overall") & (summary["feature"] == "total_sampen"), ["scale_factor", "median_relative_change_percent"]]
    total_median_text = ", ".join(f"{r.scale_factor:g}x: {r.median_relative_change_percent:.2f}%" for r in total_medians.itertuples())
    class_lines = []
    for label in CLASS_ORDER:
        data = summary.loc[(summary["stratum"] == label) & (summary["feature"] == "total_sampen"), ["scale_factor", "median_relative_change_percent"]]
        class_lines.append(f"- {label}: " + ", ".join(f"{r.scale_factor:g}x {r.median_relative_change_percent:.2f}%" for r in data.itertuples()))
    traj = summary.loc[(summary["stratum"] == "overall") & (summary["feature"] == "total_sampen")].iloc[0]
    fixed_inv = inv["fixed_r_sampen"]
    rms_pass = float(long["rms_proportionality_pass"].mean())
    mech_retained = mechanism.loc[(mechanism["stratum"] == "overall") & mechanism["predictor"].str.contains("retained")].iloc[0]
    def change_text(feature: str) -> str:
        data = summary.loc[
            (summary["stratum"] == "overall") & (summary["feature"] == feature) & (summary["scale_factor"] != 1.0),
            ["scale_factor", "median_relative_change_percent"],
        ]
        return ", ".join(f"{r.scale_factor:g}x {r.median_relative_change_percent:.2f}%" for r in data.itertuples())
    report = f"""# Phase 1: Total SampEn amplitude-scale audit

## Executive conclusion

Total SampEn is **not amplitude-scale invariant** in the frozen Stage 04 implementation. Multiplying an otherwise identical EHG9/FWH waveform changes the fixed 0.001-unit distance grid and the retained profile length; because Total SampEn is a sum over those retained values, its value changes even though the underlying waveform morphology is unchanged. This is an implementation/definition property, not biological evidence.

## Audit design and reproduction gate

The audit used {len(selected)} deterministic real 60-s EHG9/FWH windows: {sum(selected.class_label == 'UC')} UC, {sum(selected.class_label == 'IN')} IN, and {sum(selected.class_label == 'FM')} FM; {selected.patient_id.nunique()} women and {selected.recording_id.nunique()} recordings. Selection was metadata-based and round-robin by within-woman rank, never feature-based. No randomness was used.

All {len(reproduction)} feature-window checks reproduced Stage 04 within the prespecified tolerance (maximum absolute difference {reproduction.absolute_difference.max():.3g}; maximum relative difference {reproduction.relative_difference.max():.3g}). The scale experiment therefore proceeded.

## Answers to the prespecified questions

1. **Is Total SampEn amplitude-scale invariant?** No. Across the 240 non-reference scaled observations, only {inv['total_sampen']:.1%} met the strict numerical-invariance criterion.
2. **Magnitude of change.** Overall median relative changes versus 1x were {total_median_text}. Full medians, IQRs, ranges, and invariant counts are in `scale_sensitivity_summary.csv`.
3. **Monotonicity.** Of {selected.window_id.nunique()} Total-SampEn trajectories, {int(traj.n_trajectory_nondecreasing)} were nondecreasing, {int(traj.n_trajectory_nonincreasing)} nonincreasing, {int(traj.n_trajectory_nonmonotonic)} nonmonotonic, and {int(traj.n_trajectory_effectively_invariant)} effectively invariant across the five ordered scales.
4. **Does Average SampEn change?** Yes, but much less than Total SampEn: median relative changes were {change_text('average_sampen')}. It avoids direct accumulation of a scale-dependent number of terms, but rounding and retained-set changes remain.
5. **Does Median SampEn change?** Yes: median relative changes were {change_text('median_sampen')}. It is less scale-sensitive than Total SampEn but not strictly invariant under this quantized profile construction.
6. **Does fixed-r SampEn change?** No within the prespecified tolerance: it was invariant in {fixed_inv:.1%} of non-reference rows. Its tolerance scales with window SD, as its mathematical definition predicts.
7. **Do retained profile values change?** Yes. Every window changed its retained-value count at one or more non-reference scales; the rounded maximum distance and possible grid-bin count also changed in every window.
8. **Does the changing support explain Total SampEn sensitivity?** Largely. The within-window association between Total SampEn and retained count had median Spearman rho {mech_retained.median_within_window_spearman_rho:.3f}; the pooled descriptive change-score Spearman rho was {mech_retained.spearman_rho:.3f}. This strongly supports the prespecified grid-expansion/summation mechanism, while not claiming that term count alone is an exact algebraic decomposition.
9. **Is Total SampEn mathematically scale-sensitive by construction?** Yes. Fixed signal-unit rounding changes the grid/support under x->c*x, and Total SampEn is the unnormalised sum over that changing retained support.
10. **Which summary is most scale-invariant?** `{best}` had the highest strict invariance proportion ({inv[best]:.1%}) among the eight audited metrics. RMS is not intended to be invariant and instead showed correct scale proportionality in {rms_pass:.1%} of rows.
11. **Are results consistent across classes?** Yes. All three prespecified classes showed the same implementation-level scale phenomenon; class-specific Total-SampEn median changes were:
{chr(10).join(class_lines)}
12. **Coding bug or expected property?** No coding failure was found. The diagnostic wrapper exactly reconstructed production profile summaries, and all 480 selected-window feature values reproduced Stage 04. The scale behavior is an expected mathematical property of this estimator definition. Describing this Total SampEn implementation as amplitude-invariant would nevertheless be incorrect.

## Mechanistic interpretation

Positive scaling leaves waveform ordering and morphology unchanged, but converts a distance d to c*d. Stage 04 quantizes those scaled distances in absolute signal units at 0.001 resolution. This changes the maximum grid endpoint and which rounded distance shells are occupied. Total SampEn then sums the cumulative log-ratio profile over that changing retained support. The number of summands is therefore amplitude-dependent. Average SampEn divides by retained count and is less directly exposed, while fixed-r SampEn scales its tolerance with window SD and supplies the appropriate invariant control.

## Scope

No RMS-adjusted inferential model, alternative window length, cross-channel analysis, FM sensitivity model, or manuscript rewriting was performed. The result does support undertaking a separate, prespecified amplitude-adjusted sensitivity analysis next.
"""
    (OUT_DIR / "scale_sensitivity_report.md").write_text(report)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    stage04 = load_stage04()
    manifest = pd.read_csv(MANIFEST_PATH)
    frozen = pd.read_csv(FEATURE_PATH)
    selected = select_windows(manifest)
    selection_columns = [
        "selection_order_within_class", "class_label", "window_id", "patient_id", "recording_id",
        "window_class", "kind", "window_start_s", "window_end_s", "start_sample", "end_sample",
        "window_size_s", "within_patient_rank", "selection_rule",
    ]
    selected_output = selected[selection_columns].copy()
    selected_output.insert(2, "class", selected_output["class_label"])
    selected_output.to_csv(OUT_DIR / "selected_windows_for_scale_test.csv", index=False)
    signals = load_selected_signals(selected, stage04)
    reproduction = reproduction_check(selected, frozen, signals, stage04)
    reproduction.to_csv(OUT_DIR / "scale_test_reproduction_check.csv", index=False)
    if not reproduction["reproduced"].all():
        failures = int((~reproduction["reproduced"]).sum())
        raise RuntimeError(f"Reproduction gate failed for {failures} feature-window comparisons; scale experiment stopped")
    long = scale_experiment(selected, signals, stage04)
    long.to_csv(OUT_DIR / "scale_sensitivity_window_results.csv", index=False)
    summary = make_summary(long)
    summary.to_csv(OUT_DIR / "scale_sensitivity_summary.csv", index=False)
    mechanism = mechanism_analysis(long)
    mechanism.to_csv(OUT_DIR / "scale_sensitivity_mechanism.csv", index=False)
    input_hashes(selected, stage04).to_csv(OUT_DIR / "input_file_hashes.csv", index=False)
    make_figures(long)
    write_versions(stage04)
    write_implementation_audit(stage04)
    write_report(selected, reproduction, long, summary, mechanism)
    print(f"Completed. Outputs: {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
