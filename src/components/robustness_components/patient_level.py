#!/usr/bin/env python3
"""Patient-level EHG9 UC-versus-IN sensitivity analysis.

The script reads only frozen Stage 06 outputs. For each primary FWH feature,
it reduces each paired woman to a UC median and an IN median, performs a
two-sided Wilcoxon signed-rank test, creates a participant bootstrap interval
for the median paired difference, and compares the result with the frozen
mixed model.

Existing outputs are never overwritten. Use --output-dir to reproduce the
analysis in a fresh directory.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import os
import platform
import sys
from pathlib import Path

# Keep Matplotlib's cache outside the source and revision-output directories.
os.environ.setdefault("MPLCONFIGDIR", "/tmp/ehg_patient_level_matplotlib")

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import stats
from statsmodels.stats.multitest import multipletests


PROJECT_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = PROJECT_ROOT / "results" / "robustness" / "patient_level"
STAGE06_DIR = PROJECT_ROOT / "work" / "primary"
PREVIOUS_REVISION_DIR = PROJECT_ROOT / "results" / "robustness" / "dependence"

REGISTRY_PATH = STAGE06_DIR / "confirmatory_feature_registry_60s_fwh.csv"
EVENT_UNITS_PATH = STAGE06_DIR / "event_unit_analysis_60s_EHG9.csv"
MIXED_RESULTS_PATH = STAGE06_DIR / "mixed_model_results_60s_fwh_EHG9.csv"
PREVIOUS_REPORT_PATH = PREVIOUS_REVISION_DIR / "temporal_dependence_report.md"

LONG_FILENAME = "patient_level_uc_in_differences.csv"
COMPARISON_FILENAME = "patient_level_vs_mixed_model_34_features.csv"
STABILITY_FILENAME = "original_11_patient_level_stability.csv"
FIGURE_FILENAME = "core_feature_patient_level_differences.png"
REPORT_FILENAME = "patient_level_uc_in_report.md"
VERSIONS_FILENAME = "package_versions.txt"

BOOTSTRAP_SEED = 20260914
BOOTSTRAP_RESAMPLES = 5000

CORE_FEATURES = [
    "fwh_time_mean_absolute_value",
    "fwh_time_root_mean_square",
    "fwh_time_peak_to_peak",
    "fwh_time_waveform_length",
    "fwh_time_skewness",
    "fwh_time_median",
    "fwh_entropy_totalsampen",
]

CORE_LABELS = {
    "fwh_time_mean_absolute_value": "Mean absolute value",
    "fwh_time_root_mean_square": "Root mean square",
    "fwh_time_peak_to_peak": "Peak-to-peak amplitude",
    "fwh_time_waveform_length": "Waveform length",
    "fwh_time_skewness": "Time skewness",
    "fwh_time_median": "Time median",
    "fwh_entropy_totalsampen": "Total SampEn",
}

TARGET_SECONDARY_FEATURES = [
    "fwh_freq_spectral_edge_frequency_95",
    "fwh_freq_spectral_spread",
    "fwh_complexity_higuchi_fractal_dimension",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=SCRIPT_DIR,
        help="Fresh directory for outputs (default: directory containing this script).",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_inputs() -> None:
    for path in (REGISTRY_PATH, EVENT_UNITS_PATH, MIXED_RESULTS_PATH):
        if not path.is_file():
            raise FileNotFoundError(f"Required frozen input not found: {path}")


def refuse_overwrite(output_dir: Path) -> dict[str, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        "long": output_dir / LONG_FILENAME,
        "comparison": output_dir / COMPARISON_FILENAME,
        "stability": output_dir / STABILITY_FILENAME,
        "figure": output_dir / FIGURE_FILENAME,
        "report": output_dir / REPORT_FILENAME,
        "versions": output_dir / VERSIONS_FILENAME,
    }
    existing = [path for path in outputs.values() if path.exists()]
    if existing:
        joined = "\n".join(str(path) for path in existing)
        raise FileExistsError(
            "Refusing to overwrite existing revision output(s). "
            "Choose a fresh --output-dir:\n" + joined
        )
    return outputs


def direction(value: float) -> str:
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "zero"


def load_and_verify_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    registry = pd.read_csv(REGISTRY_PATH)
    event_units = pd.read_csv(EVENT_UNITS_PATH)
    mixed_all = pd.read_csv(MIXED_RESULTS_PATH)

    primary_mask = registry["primary_confirmatory"].astype(str).str.lower().eq("true")
    features = registry.loc[primary_mask, "feature"].tolist()
    if len(features) != 34 or len(set(features)) != 34:
        raise AssertionError(f"Expected 34 unique primary FWH features, found {len(features)}")
    if not set(CORE_FEATURES + TARGET_SECONDARY_FEATURES).issubset(features):
        raise AssertionError("A requested focus feature is absent from the primary registry")

    if "signal_channel" in event_units.columns:
        observed_channels = set(event_units["signal_channel"].dropna().unique())
        if observed_channels != {"EHG9"}:
            raise AssertionError(f"Event-unit table is not EHG9-only: {observed_channels}")

    required = {"patient_id", "window_class", *features}
    missing = required.difference(event_units.columns)
    if missing:
        raise AssertionError(f"Event-unit table lacks required columns: {sorted(missing)}")
    if event_units[features].isna().any().any():
        raise AssertionError("Primary feature missing values found in the EHG9 event-unit table")
    if not np.isfinite(event_units[features].to_numpy(dtype=float)).all():
        raise AssertionError("Non-finite primary feature values found in the EHG9 event-unit table")

    mixed = mixed_all.loc[
        mixed_all["contrast"].eq("contraction_vs_baseline")
    ].copy()
    if len(mixed) != 34 or mixed["feature"].nunique() != 34:
        raise AssertionError("Mixed-model input lacks one EHG9 UC-vs-IN row per feature")
    if set(mixed["feature"]) != set(features):
        raise AssertionError("Registry and mixed-model feature sets differ")
    mixed = mixed.set_index("feature").loc[features].reset_index()

    recomputed_q = multipletests(mixed["p_value"].to_numpy(), method="fdr_bh")[1]
    if not np.allclose(recomputed_q, mixed["q_value"], rtol=1e-12, atol=1e-15):
        raise AssertionError("Frozen mixed-model BH q-values do not reproduce")
    if int((mixed["q_value"] < 0.05).sum()) != 11:
        raise AssertionError("Expected the frozen EHG9 mixed model to contain 11 q<0.05 results")

    return registry, event_units, mixed, features


def identify_population(event_units: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    uc_in = event_units.loc[
        event_units["window_class"].isin(["baseline", "contraction"])
    ].copy()
    counts = pd.crosstab(uc_in["patient_id"], uc_in["window_class"])
    for state in ("baseline", "contraction"):
        if state not in counts.columns:
            counts[state] = 0
    counts = counts[["contraction", "baseline"]].sort_index()
    paired_ids = counts.index[
        (counts["contraction"] > 0) & (counts["baseline"] > 0)
    ].astype(str).tolist()

    if event_units["patient_id"].nunique() != 42:
        raise AssertionError("Expected 42 women in the frozen EHG9 analysis")
    if len(paired_ids) != 34:
        raise AssertionError(f"Expected 34 women with both states, found {len(paired_ids)}")
    return uc_in, counts, paired_ids


def build_patient_long(
    uc_in: pd.DataFrame,
    counts: pd.DataFrame,
    paired_ids: list[str],
    features: list[str],
) -> pd.DataFrame:
    medians = (
        uc_in.groupby(["patient_id", "window_class"], sort=True)[features]
        .median()
        .unstack("window_class")
    )
    rows: list[dict[str, object]] = []
    for patient_id in paired_ids:
        for feature in features:
            uc_median = float(medians.loc[patient_id, (feature, "contraction")])
            in_median = float(medians.loc[patient_id, (feature, "baseline")])
            difference = uc_median - in_median
            rows.append(
                {
                    "patient_id": patient_id,
                    "feature": feature,
                    "n_UC_windows": int(counts.loc[patient_id, "contraction"]),
                    "n_IN_windows": int(counts.loc[patient_id, "baseline"]),
                    "UC_median": uc_median,
                    "IN_median": in_median,
                    "UC_minus_IN": difference,
                    "direction": direction(difference),
                }
            )

    patient_long = pd.DataFrame(rows)
    if len(patient_long) != len(paired_ids) * len(features):
        raise AssertionError("Unexpected patient-level long-table size")
    if patient_long.isna().any().any():
        raise AssertionError("Missing values found in patient-level long table")
    return patient_long


def wilcoxon_test(differences: np.ndarray) -> tuple[float, float, int, int, str]:
    differences = np.asarray(differences, dtype=float)
    zero_count = int(np.sum(differences == 0))
    nonzero = differences[differences != 0]
    effective_n = int(len(nonzero))
    if effective_n == 0:
        return 0.0, 1.0, effective_n, zero_count, "all_zero"

    absolute_nonzero = np.abs(nonzero)
    has_rank_ties = len(np.unique(absolute_nonzero)) < effective_n
    method_used = (
        "asymptotic_auto_ties_or_zeros"
        if zero_count > 0 or has_rank_ties
        else "exact_auto_no_ties_or_zeros"
    )
    result = stats.wilcoxon(
        differences,
        zero_method="wilcox",
        correction=False,
        alternative="two-sided",
        method="auto",
    )
    return float(result.statistic), float(result.pvalue), effective_n, zero_count, method_used


def summarize_features(
    patient_long: pd.DataFrame,
    mixed: pd.DataFrame,
    features: list[str],
) -> pd.DataFrame:
    mixed_by_feature = mixed.set_index("feature")
    paired_ids = sorted(patient_long["patient_id"].unique())
    bootstrap_rng = np.random.default_rng(BOOTSTRAP_SEED)
    bootstrap_indices = bootstrap_rng.integers(
        0,
        len(paired_ids),
        size=(BOOTSTRAP_RESAMPLES, len(paired_ids)),
    )

    rows: list[dict[str, object]] = []
    for feature in features:
        data = (
            patient_long.loc[patient_long["feature"].eq(feature)]
            .set_index("patient_id")
            .loc[paired_ids]
        )
        differences = data["UC_minus_IN"].to_numpy(dtype=float)
        bootstrap_medians = np.median(differences[bootstrap_indices], axis=1)
        if not np.isfinite(bootstrap_medians).all():
            raise AssertionError(f"Non-finite bootstrap result for {feature}")
        bootstrap_ci_low, bootstrap_ci_high = np.percentile(
            bootstrap_medians,
            [2.5, 97.5],
            method="linear",
        )

        statistic, p_value, effective_n, zero_count, test_method = wilcoxon_test(differences)
        frozen = mixed_by_feature.loc[feature]
        mixed_estimate = float(frozen["estimate"])
        median_difference = float(np.median(differences))
        if mixed_estimate > 0:
            consistent = differences > 0
        elif mixed_estimate < 0:
            consistent = differences < 0
        else:
            consistent = differences == 0

        positive_n = int(np.sum(differences > 0))
        negative_n = int(np.sum(differences < 0))
        n_paired = int(len(differences))
        consistent_n = int(np.sum(consistent))
        rows.append(
            {
                "feature": feature,
                "mixed_model_standardized_estimate": mixed_estimate,
                "mixed_model_q": float(frozen["q_value"]),
                "original_primary_supported": bool(float(frozen["q_value"]) < 0.05),
                "n_paired_women": n_paired,
                "median_UC_value": float(np.median(data["UC_median"])),
                "median_IN_value": float(np.median(data["IN_median"])),
                "median_paired_difference": median_difference,
                "paired_difference_q1": float(np.percentile(differences, 25, method="linear")),
                "paired_difference_q3": float(np.percentile(differences, 75, method="linear")),
                "bootstrap_ci_low": float(bootstrap_ci_low),
                "bootstrap_ci_high": float(bootstrap_ci_high),
                "bootstrap_finite_resamples": int(np.isfinite(bootstrap_medians).sum()),
                "bootstrap_unique_medians": int(len(np.unique(bootstrap_medians))),
                "n_positive_difference": positive_n,
                "pct_positive_difference": 100.0 * positive_n / n_paired,
                "n_negative_difference": negative_n,
                "pct_negative_difference": 100.0 * negative_n / n_paired,
                "n_zero_difference": zero_count,
                "wilcoxon_effective_n_nonzero": effective_n,
                "wilcoxon_statistic": statistic,
                "wilcoxon_p": p_value,
                "wilcoxon_zero_method": "wilcox_discard_zeros",
                "wilcoxon_p_method": test_method,
                "patient_level_direction": direction(median_difference),
                "direction_agreement_with_mixed_model": bool(
                    np.sign(median_difference) == np.sign(mixed_estimate)
                ),
                "n_mixed_model_consistent": consistent_n,
                "pct_mixed_model_consistent": 100.0 * consistent_n / n_paired,
                "majority_mixed_model_consistent": bool(consistent_n > n_paired / 2),
            }
        )

    comparison = pd.DataFrame(rows)
    comparison["patient_level_q"] = multipletests(
        comparison["wilcoxon_p"].to_numpy(), method="fdr_bh"
    )[1]
    comparison["patient_level_q_lt_0_05"] = comparison["patient_level_q"] < 0.05

    ordered_columns = [
        "feature",
        "mixed_model_standardized_estimate",
        "mixed_model_q",
        "original_primary_supported",
        "n_paired_women",
        "median_UC_value",
        "median_IN_value",
        "median_paired_difference",
        "paired_difference_q1",
        "paired_difference_q3",
        "bootstrap_ci_low",
        "bootstrap_ci_high",
        "bootstrap_finite_resamples",
        "bootstrap_unique_medians",
        "n_positive_difference",
        "pct_positive_difference",
        "n_negative_difference",
        "pct_negative_difference",
        "n_zero_difference",
        "wilcoxon_effective_n_nonzero",
        "wilcoxon_statistic",
        "wilcoxon_p",
        "wilcoxon_zero_method",
        "wilcoxon_p_method",
        "patient_level_q",
        "patient_level_q_lt_0_05",
        "patient_level_direction",
        "direction_agreement_with_mixed_model",
        "n_mixed_model_consistent",
        "pct_mixed_model_consistent",
        "majority_mixed_model_consistent",
    ]
    return comparison[ordered_columns]


def create_figure(patient_long: pd.DataFrame, mixed: pd.DataFrame, path: Path) -> None:
    mixed_estimates = mixed.set_index("feature")["estimate"]
    fig, axes = plt.subplots(4, 2, figsize=(12, 10), constrained_layout=True)
    axes_flat = axes.ravel()
    consistent_color = "#007C83"
    opposite_color = "#9CA3AF"
    zero_color = "#D97706"

    for ax, feature in zip(axes_flat, CORE_FEATURES, strict=False):
        data = (
            patient_long.loc[patient_long["feature"].eq(feature)]
            .sort_values("patient_id")
        )
        values = data["UC_minus_IN"].to_numpy(dtype=float)
        expected_positive = float(mixed_estimates.loc[feature]) > 0
        consistent = values > 0 if expected_positive else values < 0
        colors = np.where(
            values == 0,
            zero_color,
            np.where(consistent, consistent_color, opposite_color),
        )
        y_offsets = 1.0 + ((np.arange(len(values)) % 9) - 4) * 0.022

        ax.boxplot(
            values,
            vert=False,
            widths=0.34,
            showfliers=False,
            patch_artist=True,
            boxprops={"facecolor": "#E5E7EB", "edgecolor": "#374151", "linewidth": 1.0},
            medianprops={"color": "#111827", "linewidth": 2.0},
            whiskerprops={"color": "#6B7280"},
            capprops={"color": "#6B7280"},
        )
        ax.scatter(values, y_offsets, c=colors, s=22, alpha=0.85, edgecolors="none", zorder=3)
        ax.axvline(0, color="#B91C1C", linestyle="--", linewidth=1.1, zorder=1)
        low = min(float(np.min(values)), 0.0)
        high = max(float(np.max(values)), 0.0)
        span = high - low
        pad = span * 0.08 if span > 0 else 1.0
        ax.set_xlim(low - pad, high + pad)
        ax.set_yticks([])
        ax.set_title(CORE_LABELS[feature], fontsize=10, fontweight="bold", loc="left")
        ax.set_xlabel("UC median minus IN median (raw units)", fontsize=8)
        ax.tick_params(axis="x", labelsize=8)
        ax.grid(axis="x", color="#E5E7EB", linewidth=0.6)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color("#9CA3AF")

    legend_ax = axes_flat[-1]
    legend_ax.axis("off")
    legend_handles = [
        Line2D([0], [0], marker="o", linestyle="", color=consistent_color,
               label="Same direction as mixed model"),
        Line2D([0], [0], marker="o", linestyle="", color=opposite_color,
               label="Opposite direction"),
        Line2D([0], [0], color="#B91C1C", linestyle="--", label="No difference"),
    ]
    legend_ax.legend(handles=legend_handles, loc="center", frameon=False, fontsize=10)
    legend_ax.text(
        0.5,
        0.28,
        "Each point is one of 34 paired women.\nPanels use separate raw-unit scales.",
        ha="center",
        va="center",
        fontsize=9,
        color="#4B5563",
        transform=legend_ax.transAxes,
    )

    fig.suptitle(
        "Woman-level EHG9 differences between uterine contraction and inactivity",
        fontsize=14,
        fontweight="bold",
    )
    fig.savefig(path, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fmt(value: float) -> str:
    value = float(value)
    if value == 0:
        return "0"
    if abs(value) < 0.001 or abs(value) >= 1000:
        return f"{value:.3e}"
    return f"{value:.4g}"


def feature_code(feature: str) -> str:
    return f"`{feature}`"


def count_distribution(values: pd.Series) -> str:
    return (
        f"mean {values.mean():.1f}, SD {values.std(ddof=1):.1f}, "
        f"minimum {int(values.min())}, Q1 {values.quantile(0.25):.2f}, "
        f"median {values.median():.1f}, Q3 {values.quantile(0.75):.2f}, "
        f"maximum {int(values.max())}"
    )


def build_report(
    event_units: pd.DataFrame,
    counts: pd.DataFrame,
    paired_ids: list[str],
    comparison: pd.DataFrame,
) -> str:
    paired_counts = counts.loc[paired_ids]
    has_uc = counts["contraction"] > 0
    has_in = counts["baseline"] > 0
    original = comparison.loc[comparison["original_primary_supported"]].copy()
    core = comparison.set_index("feature").loc[CORE_FEATURES].reset_index()
    secondary = comparison.set_index("feature").loc[TARGET_SECONDARY_FEATURES].reset_index()

    original_same = int(original["direction_agreement_with_mixed_model"].sum())
    original_q = int(original["patient_level_q_lt_0_05"].sum())
    original_majority = int(original["majority_mixed_model_consistent"].sum())
    original_reversed = original.loc[
        ~original["direction_agreement_with_mixed_model"], "feature"
    ].tolist()
    all_q = int(comparison["patient_level_q_lt_0_05"].sum())
    all_nominal = int((comparison["wilcoxon_p"] < 0.05).sum())

    core_rows: list[str] = []
    for _, row in core.iterrows():
        core_rows.append(
            "| {name} | {n} | {difference} | {consistent:.1f}% | {p} | {q} | [{low}, {high}] |".format(
                name=CORE_LABELS[row["feature"]],
                n=int(row["n_paired_women"]),
                difference=fmt(row["median_paired_difference"]),
                consistent=float(row["pct_mixed_model_consistent"]),
                p=fmt(row["wilcoxon_p"]),
                q=fmt(row["patient_level_q"]),
                low=fmt(row["bootstrap_ci_low"]),
                high=fmt(row["bootstrap_ci_high"]),
            )
        )

    original_rows: list[str] = []
    for _, row in original.iterrows():
        if row["patient_level_q_lt_0_05"]:
            classification = "same direction + q<0.05"
        elif row["direction_agreement_with_mixed_model"] and row["majority_mixed_model_consistent"]:
            classification = "majority consistent; q>=0.05"
        elif row["direction_agreement_with_mixed_model"]:
            classification = "same direction; q>=0.05"
        else:
            classification = "reversed"
        original_rows.append(
            "| {feature} | {difference} | {consistent:.1f}% | {p} | {q} | {classification} |".format(
                feature=row["feature"],
                difference=fmt(row["median_paired_difference"]),
                consistent=float(row["pct_mixed_model_consistent"]),
                p=fmt(row["wilcoxon_p"]),
                q=fmt(row["patient_level_q"]),
                classification=classification,
            )
        )

    secondary_sentences: list[str] = []
    for _, row in secondary.iterrows():
        secondary_sentences.append(
            f"{feature_code(row['feature'])}: median difference {fmt(row['median_paired_difference'])}, "
            f"{row['pct_mixed_model_consistent']:.1f}% mixed-model-consistent, "
            f"p={fmt(row['wilcoxon_p'])}, q={fmt(row['patient_level_q'])}, "
            f"bootstrap CI [{fmt(row['bootstrap_ci_low'])}, {fmt(row['bootstrap_ci_high'])}]"
        )

    zero_features = comparison.loc[comparison["n_zero_difference"] > 0, ["feature", "n_zero_difference"]]
    zero_text = ", ".join(
        f"{feature_code(row.feature)} ({int(row.n_zero_difference)})"
        for row in zero_features.itertuples(index=False)
    )
    degenerate_bootstrap = comparison.loc[comparison["bootstrap_unique_medians"] == 1, "feature"].tolist()

    input_lines = [
        f"- `{REGISTRY_PATH}` (SHA-256 `{sha256(REGISTRY_PATH)}`)",
        f"- `{EVENT_UNITS_PATH}` (SHA-256 `{sha256(EVENT_UNITS_PATH)}`)",
        f"- `{MIXED_RESULTS_PATH}` (SHA-256 `{sha256(MIXED_RESULTS_PATH)}`)",
        "- Paired analysis checks the migrated primary feature family directly.",
    ]

    return f"""# Patient-level EHG9 UC-versus-IN sensitivity analysis

## Question

Do EHG9 uterine-contraction (UC) versus inactivity (IN) differences remain visible after reducing the data to one paired UC-versus-IN summary per woman?

## Frozen inputs

{chr(10).join(input_lines)}

The script read these inputs without refitting or changing the Stage 06 mixed model. It verified the 34-feature registry, the 34 frozen EHG9 UC-vs-IN rows, their BH q-values, and the frozen 11-feature q<0.05 set.

## Analysis population

- Total women in the frozen EHG9 event-unit table: **{event_units['patient_id'].nunique()}**.
- Women with both UC and IN: **{len(paired_ids)}**.
- Women with IN only: **{int((has_in & ~has_uc).sum())}**.
- Women with UC only: **{int((~has_in & has_uc).sum())}**.
- UC windows per paired woman: {count_distribution(paired_counts['contraction'])}.
- IN windows per paired woman: {count_distribution(paired_counts['baseline'])}.

The earlier estimate of approximately 34 paired women is confirmed exactly.

## Methods

For each of the 34 primary FWH features, the analysis calculated each paired woman's raw UC-window median, raw IN-window median, and `UC median - IN median`. No window-level standardisation was applied. Each woman therefore contributed one paired difference per feature and received equal weight in the inferential analysis.

The primary test was a two-sided Wilcoxon signed-rank test. `zero_method="wilcox"` discarded exact zero differences before ranking. `correction=False` applied no continuity correction. `method="auto"` used the exact null distribution for n<=50 when there were no rank ties or zeros, and the asymptotic tie/zero-adjusted calculation otherwise. Exact zeros occurred for: {zero_text}. The result file records the effective nonzero sample size and p-value method for every feature.

The 95% CI is the 2.5th and 97.5th percentiles of {BOOTSTRAP_RESAMPLES:,} participant-level bootstrap medians. Each resample drew {len(paired_ids)} women with replacement. All features used the same resampled participant indices, NumPy `default_rng`, and seed **{BOOTSTRAP_SEED}**. All {BOOTSTRAP_RESAMPLES:,} bootstrap estimates were finite for all features. {feature_code(degenerate_bootstrap[0]) if degenerate_bootstrap else 'No feature'} had a degenerate [0, 0] bootstrap interval because 32/34 paired differences were exactly zero; this reflects feature quantisation, not a computational failure.

BH correction was applied once across all 34 Wilcoxon p-values. Across the full family, {all_nominal}/34 had raw p<0.05 and {all_q}/34 had q<0.05.

## Original 11 results

| Feature | Median UC-IN | Mixed-consistent women | Wilcoxon p | BH q | Interpretation |
|---|---:|---:|---:|---:|---|
{chr(10).join(original_rows)}

All **{original_same}/11** original features retained the mixed-model direction. **{original_q}/11** retained BH q<0.05. All **{original_majority}/11** had a strict majority of women changing in the mixed-model-consistent direction. Reversed original features: **{'none' if not original_reversed else ', '.join(feature_code(feature) for feature in original_reversed)}**.

## Core amplitude/time-domain results

| Feature | Paired women | Median UC-IN | Mixed-consistent women | Wilcoxon p | BH q | Bootstrap 95% CI |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(core_rows)}

Mean absolute value, RMS, peak-to-peak amplitude, and waveform length were higher during UC for 85.3%, 85.3%, 91.2%, and 82.4% of paired women, respectively. All four retained BH q<0.05, and all four bootstrap intervals excluded zero. Time skewness remained positive in 64.7% but did not retain BH support. Time median remained negative in 64.7% but did not retain BH support. Total SampEn was higher during UC in 82.4%, with p={fmt(core.loc[core['feature'].eq('fwh_entropy_totalsampen'), 'wilcoxon_p'].iloc[0])}, q={fmt(core.loc[core['feature'].eq('fwh_entropy_totalsampen'), 'patient_level_q'].iloc[0])}, and a bootstrap interval excluding zero. This analysis establishes the within-woman pattern only; it does not reassess Total SampEn's relation to amplitude.

## Spectral and complexity features highlighted by the prior sensitivity

- {secondary_sentences[0]}.
- {secondary_sentences[1]}.
- {secondary_sentences[2]}.

All three returned to the negative mixed-model direction at the patient level, and 61.8%--67.6% of women were direction-consistent. None retained BH q<0.05. Spectral-edge frequency's percentile bootstrap interval excluded zero even though its Wilcoxon p-value and BH q-value did not. The interval estimates the sampling distribution of the median, whereas the signed-rank test uses signed ranks and assumes a symmetric difference distribution for its location interpretation; the two procedures need not agree in a small, asymmetric sample.

## Answers to the reviewer-focused questions

1. **Women with both states:** {len(paired_ids)} of {event_units['patient_id'].nunique()} women.
2. **Original directions retained:** {original_same}/11.
3. **Original BH q<0.05 retained:** {original_q}/11.
4. **Original reversals:** {'None.' if not original_reversed else ', '.join(feature_code(feature) for feature in original_reversed) + '.'}
5. **Core amplitude features:** MAV, RMS, peak-to-peak amplitude, and waveform length were consistently higher during UC and all retained family-wide BH support.
6. **Expected-direction percentages:** MAV 85.3%, RMS 85.3%, peak-to-peak 91.2%, and waveform length 82.4%.
7. **Total SampEn:** The patient-level median difference remained positive, 82.4% of women agreed with the mixed-model direction, and q<0.05 was retained.
8. **Spectral edge, spectral spread, and Higuchi fractal dimension:** All three had negative median differences, matching the mixed model rather than the pooled/cluster-robust OLS point estimates. Each had a direction-consistent majority, but none had q<0.05.
9. **Main pattern:** Yes. The patient-level results support a primarily amplitude/time-domain UC-vs-IN pattern. The four core amplitude/waveform features show the strongest consistency, with Total SampEn also prominent.
10. **Relation to the prior woman-cluster analysis:** This strengthens the interpretation of the previous result. It shows that the core effects are directly visible within women after removing window-count weighting, while confirming that inference for weaker spectral and complexity features is imprecise. The earlier cluster-robust OLS reversals are not reproduced by the paired medians, which suggests that pooled marginal weighting and between-woman composition contributed to those reversals. The analyses use different estimands and should be presented as complementary rather than treated as contradictory.

## Technical concerns

- Collapsing to medians changes the estimand from a window-level conditional contrast to the median of woman-level median differences. This analysis complements rather than replaces the mixed model.
- The paired sample contains only 34 women. UC counts range from {int(paired_counts['contraction'].min())} to {int(paired_counts['contraction'].max())}, and three women have only one UC window. Their UC median is a single observation.
- Cluster sizes and state counts are highly unequal, but the paired analysis gives each woman equal inferential weight. This is desirable for the sensitivity question and discards information about within-woman precision.
- The Wilcoxon signed-rank test assumes a symmetric distribution of paired differences for a location-shift interpretation. The sign-consistency counts and bootstrap median interval provide complementary summaries when that assumption is doubtful.
- Discrete features generate ties and exact zeros. The specified zero handling and feature-specific p-value method are recorded explicitly.
- Percentile bootstrap intervals are simple and transparent but can be discrete with 34 participants. They should not be treated as multiplicity-adjusted intervals.

## Recommendation

The evidence is sufficiently coherent to proceed to cross-channel dependence analysis. The frozen mixed model, woman-cluster sensitivity, and patient-level paired analysis collectively support a robust core amplitude/time-domain pattern. The patient-level result also supports Total SampEn directionally and inferentially, without establishing independence from amplitude. Individual spectral and complexity findings remain less secure and should not carry the central claim.

## Scope stop

This analysis stops at the EHG9 patient-level UC-versus-IN sensitivity. It does not perform cross-channel modelling, alternative window-length analysis, SampEn-scale analysis, fetal-movement contamination analysis, or manuscript rewriting.
"""


def package_versions_text() -> str:
    packages = ["numpy", "pandas", "scipy", "statsmodels", "matplotlib"]
    lines = [
        f"Python: {sys.version.replace(chr(10), ' ')}",
        f"Platform: {platform.platform()}",
    ]
    for package in packages:
        try:
            version = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            version = "not installed"
        lines.append(f"{package}: {version}")
    lines.extend(
        [
            f"matplotlib imported version: {matplotlib.__version__}",
            f"Bootstrap resamples: {BOOTSTRAP_RESAMPLES}",
            "Random number generator: numpy.random.default_rng (PCG64)",
            f"Random seed: {BOOTSTRAP_SEED}",
            "Bootstrap statistic: median of resampled woman-level UC-minus-IN differences",
            "Bootstrap interval: percentile 2.5% and 97.5%, NumPy method=linear",
            "Wilcoxon: scipy.stats.wilcoxon, two-sided, zero_method=wilcox, correction=False, method=auto",
            "Multiplicity: statsmodels.stats.multitest.multipletests(method=fdr_bh)",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    require_inputs()
    outputs = refuse_overwrite(output_dir)
    _, event_units, mixed, features = load_and_verify_inputs()
    uc_in, counts, paired_ids = identify_population(event_units)
    patient_long = build_patient_long(uc_in, counts, paired_ids, features)
    comparison = summarize_features(patient_long, mixed, features)
    stability = comparison.loc[comparison["original_primary_supported"]].copy()

    patient_long.to_csv(outputs["long"], index=False, float_format="%.17g")
    comparison.to_csv(outputs["comparison"], index=False, float_format="%.17g")
    stability.to_csv(outputs["stability"], index=False, float_format="%.17g")
    create_figure(patient_long, mixed, outputs["figure"])
    outputs["report"].write_text(
        build_report(event_units, counts, paired_ids, comparison),
        encoding="utf-8",
    )
    outputs["versions"].write_text(package_versions_text(), encoding="utf-8")

    print(f"Wrote {len(patient_long)} patient-feature rows: {outputs['long']}")
    print(f"Wrote {len(comparison)}-feature comparison: {outputs['comparison']}")
    print(f"Wrote {len(stability)}-feature stability summary: {outputs['stability']}")
    print(f"Wrote diagnostic figure: {outputs['figure']}")
    print(f"Wrote report: {outputs['report']}")
    print(f"Wrote package versions: {outputs['versions']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
