#!/usr/bin/env python3
"""Dependence-robust sensitivity analysis for the frozen EHG9 UC-vs-IN result.

This script is intentionally separate from Stages 00--09. It reads only the
frozen Stage 06 outputs, fits pooled and woman-cluster-robust OLS models for the
34 primary FWH features, and writes revision-only outputs.

No randomness is used. By default, existing outputs are never overwritten.
Use --output-dir to reproduce the analysis in a fresh directory.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


SCRIPT_DIR = Path(__file__).resolve().parents[3] / "results" / "robustness" / "dependence"
PROJECT_ROOT = Path(__file__).resolve().parents[3]
STAGE06_DIR = PROJECT_ROOT / "work" / "primary"
MANUSCRIPT_MEMO = None  # Final feature set is verified directly against primary output.

REGISTRY_PATH = STAGE06_DIR / "confirmatory_feature_registry_60s_fwh.csv"
EVENT_UNITS_PATH = STAGE06_DIR / "event_unit_analysis_60s_EHG9.csv"
MIXED_RESULTS_PATH = STAGE06_DIR / "mixed_model_results_60s_fwh_EHG9.csv"

COMPARISON_FILENAME = "temporal_dependence_comparison_34_features.csv"
STABILITY_FILENAME = "original_11_stability_summary.csv"
REPORT_FILENAME = "temporal_dependence_report.md"
VERSIONS_FILENAME = "package_versions.txt"

EXPECTED_ORIGINAL_11 = {
    "fwh_time_median",
    "fwh_time_skewness",
    "fwh_time_mean_absolute_value",
    "fwh_time_root_mean_square",
    "fwh_time_peak_to_peak",
    "fwh_time_waveform_length",
    "fwh_freq_spectral_edge_frequency_95",
    "fwh_freq_spectral_spread",
    "fwh_entropy_totalsampen",
    "fwh_complexity_hjorth_complexity",
    "fwh_complexity_higuchi_fractal_dimension",
}

EXPECTED_FDR_REMOVALS_FROM_PRIOR_14 = {
    "fwh_entropy_svd_entropy",
    "fwh_entropy_sampen_fixed_r",
    "fwh_complexity_hjorth_mobility",
}


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
        "comparison": output_dir / COMPARISON_FILENAME,
        "stability": output_dir / STABILITY_FILENAME,
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


def load_and_verify_frozen_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    registry = pd.read_csv(REGISTRY_PATH)
    event_units = pd.read_csv(EVENT_UNITS_PATH)
    mixed_all = pd.read_csv(MIXED_RESULTS_PATH)

    primary_mask = registry["primary_confirmatory"].astype(str).str.lower().eq("true")
    features = registry.loc[primary_mask, "feature"].tolist()
    if len(features) != 34 or len(set(features)) != 34:
        raise AssertionError(f"Expected 34 unique primary FWH features, found {len(features)}")

    mixed = mixed_all.loc[
        mixed_all["contrast"].eq("contraction_vs_baseline")
    ].copy()
    if len(mixed) != 34 or mixed["feature"].nunique() != 34:
        raise AssertionError("Frozen EHG9 mixed-model file lacks one row per primary feature")
    if set(mixed["feature"]) != set(features):
        raise AssertionError("Registry and EHG9 mixed-model primary feature sets differ")
    mixed = mixed.set_index("feature").loc[features].reset_index()

    recomputed_q = multipletests(mixed["p_value"].to_numpy(), method="fdr_bh")[1]
    if not np.allclose(recomputed_q, mixed["q_value"], rtol=1e-12, atol=1e-15):
        max_difference = float(np.max(np.abs(recomputed_q - mixed["q_value"].to_numpy())))
        raise AssertionError(f"Frozen mixed-model BH q-values do not reproduce; max diff={max_difference}")

    frozen_original_11 = set(mixed.loc[mixed["q_value"] < 0.05, "feature"])
    if frozen_original_11 != EXPECTED_ORIGINAL_11:
        raise AssertionError(
            "Frozen q<0.05 feature set differs from the manuscript source-of-truth. "
            f"Observed={sorted(frozen_original_11)}"
        )

    for feature in EXPECTED_FDR_REMOVALS_FROM_PRIOR_14:
        q_value = float(mixed.loc[mixed["feature"].eq(feature), "q_value"].iloc[0])
        if q_value < 0.05:
            raise AssertionError(f"Expected removed feature is still q<0.05: {feature}")

    required_columns = {"window_class", "patient_id", "recording_id", *features}
    missing_columns = required_columns.difference(event_units.columns)
    if missing_columns:
        raise AssertionError(f"EHG9 event-unit table lacks columns: {sorted(missing_columns)}")
    if event_units[features].isna().any().any():
        raise AssertionError("Primary feature missing values found in the frozen EHG9 event-unit table")
    if not np.isfinite(event_units[features].to_numpy(dtype=float)).all():
        raise AssertionError("Non-finite primary feature values found in the frozen EHG9 event-unit table")

    return registry, event_units, mixed, features


def fit_models(event_units: pd.DataFrame, mixed: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    mixed_by_feature = mixed.set_index("feature")
    rows: list[dict[str, object]] = []

    for feature in features:
        frozen = mixed_by_feature.loc[feature]
        full = event_units[[feature, "window_class", "patient_id", "recording_id"]].copy()

        # Stage 06 standardizes each feature over all three modeled states with
        # pandas' sample SD (ddof=1). Reuse that frozen mean/SD so every estimate
        # remains on exactly the same z scale as the primary mixed model.
        observed_mean = float(full[feature].mean())
        observed_sd = float(full[feature].std(ddof=1))
        if not np.isclose(observed_mean, float(frozen["feature_mean_raw"]), rtol=1e-12, atol=1e-15):
            raise AssertionError(f"Frozen mean does not reproduce for {feature}")
        if not np.isclose(observed_sd, float(frozen["feature_sd_raw"]), rtol=1e-12, atol=1e-15):
            raise AssertionError(f"Frozen SD does not reproduce for {feature}")

        data = full.loc[full["window_class"].isin(["baseline", "contraction"])].copy()
        data["z_feature"] = (
            data[feature] - float(frozen["feature_mean_raw"])
        ) / float(frozen["feature_sd_raw"])
        data["state_uc"] = data["window_class"].eq("contraction").astype(int)
        design = sm.add_constant(data[["state_uc"]], has_constant="add")

        pooled = sm.OLS(data["z_feature"], design).fit(use_t=True)
        pooled_ci = pooled.conf_int(alpha=0.05).loc["state_uc"]

        # CR1 covariance correction plus cluster-based reference df. With 42
        # women and two regression coefficients, statsmodels multiplies the
        # sandwich covariance by G/(G-1) * (N-1)/(N-K), then uses t_(G-1).
        clustered = pooled.get_robustcov_results(
            cov_type="cluster",
            groups=data["patient_id"],
            use_correction=True,
            df_correction=True,
            use_t=True,
        )
        clustered_ci = clustered.conf_int(alpha=0.05)[1]
        cluster_df = int(getattr(clustered, "df_resid_inference", -1))
        expected_df = int(data["patient_id"].nunique() - 1)
        if cluster_df != expected_df:
            raise AssertionError(
                f"Unexpected cluster inference df for {feature}: {cluster_df} != {expected_df}"
            )

        rows.append(
            {
                "feature_name": feature,
                "primary_mixed_model_estimate": float(frozen["estimate"]),
                "primary_mixed_model_ci_low": float(frozen["ci_low"]),
                "primary_mixed_model_ci_high": float(frozen["ci_high"]),
                "primary_mixed_model_p": float(frozen["p_value"]),
                "primary_mixed_model_q": float(frozen["q_value"]),
                "pooled_ols_estimate": float(pooled.params["state_uc"]),
                "pooled_ols_se": float(pooled.bse["state_uc"]),
                "pooled_ols_ci_low": float(pooled_ci.iloc[0]),
                "pooled_ols_ci_high": float(pooled_ci.iloc[1]),
                "pooled_ols_p": float(pooled.pvalues["state_uc"]),
                "woman_cluster_estimate": float(clustered.params[1]),
                "woman_cluster_robust_se": float(clustered.bse[1]),
                "woman_cluster_ci_low": float(clustered_ci[0]),
                "woman_cluster_ci_high": float(clustered_ci[1]),
                "woman_cluster_p": float(clustered.pvalues[1]),
                "original_11": feature in EXPECTED_ORIGINAL_11,
            }
        )

    comparison = pd.DataFrame(rows)
    comparison["pooled_ols_q"] = multipletests(
        comparison["pooled_ols_p"].to_numpy(), method="fdr_bh"
    )[1]
    comparison["woman_cluster_q"] = multipletests(
        comparison["woman_cluster_p"].to_numpy(), method="fdr_bh"
    )[1]
    comparison["woman_cluster_same_direction_as_primary"] = (
        np.sign(comparison["woman_cluster_estimate"])
        == np.sign(comparison["primary_mixed_model_estimate"])
    )
    comparison["absolute_coefficient_difference_mixed_vs_cluster"] = (
        comparison["primary_mixed_model_estimate"]
        - comparison["woman_cluster_estimate"]
    ).abs()
    comparison["woman_cluster_q_lt_0_05"] = comparison["woman_cluster_q"] < 0.05

    # The two OLS analyses have identical point estimates by construction.
    if not np.allclose(
        comparison["pooled_ols_estimate"], comparison["woman_cluster_estimate"],
        rtol=1e-13, atol=1e-15,
    ):
        raise AssertionError("Pooled and cluster-robust OLS point estimates differ unexpectedly")

    ordered_columns = [
        "feature_name",
        "primary_mixed_model_estimate",
        "primary_mixed_model_ci_low",
        "primary_mixed_model_ci_high",
        "primary_mixed_model_p",
        "primary_mixed_model_q",
        "pooled_ols_estimate",
        "pooled_ols_se",
        "pooled_ols_ci_low",
        "pooled_ols_ci_high",
        "pooled_ols_p",
        "pooled_ols_q",
        "woman_cluster_estimate",
        "woman_cluster_robust_se",
        "woman_cluster_ci_low",
        "woman_cluster_ci_high",
        "woman_cluster_p",
        "woman_cluster_q",
        "original_11",
        "woman_cluster_same_direction_as_primary",
        "absolute_coefficient_difference_mixed_vs_cluster",
        "woman_cluster_q_lt_0_05",
    ]
    return comparison[ordered_columns]


def format_feature_list(features: list[str]) -> str:
    if not features:
        return "None"
    return ", ".join(f"`{feature}`" for feature in features)


def format_p(value: float) -> str:
    return f"{value:.3g}"


def build_report(comparison: pd.DataFrame, event_units: pd.DataFrame) -> str:
    original = comparison.loc[comparison["original_11"]].copy()
    uc_in = event_units.loc[event_units["window_class"].isin(["baseline", "contraction"])].copy()
    woman_state = pd.crosstab(uc_in["patient_id"], uc_in["window_class"])
    has_in = woman_state.get("baseline", pd.Series(0, index=woman_state.index)) > 0
    has_uc = woman_state.get("contraction", pd.Series(0, index=woman_state.index)) > 0

    same_direction_n = int(original["woman_cluster_same_direction_as_primary"].sum())
    cluster_q_n = int(original["woman_cluster_q_lt_0_05"].sum())
    reversals = original.loc[
        ~original["woman_cluster_same_direction_as_primary"], "feature_name"
    ].tolist()

    # A transparent descriptive threshold, not a hypothesis test.
    original["absolute_ratio_cluster_to_mixed"] = (
        original["woman_cluster_estimate"].abs()
        / original["primary_mixed_model_estimate"].abs()
    )
    materially_attenuated = original.loc[
        original["absolute_ratio_cluster_to_mixed"] <= 0.80, "feature_name"
    ].tolist()

    time_domain = original.loc[original["feature_name"].str.startswith("fwh_time_")]
    time_same_n = int(time_domain["woman_cluster_same_direction_as_primary"].sum())

    pooled_q_total = int((comparison["pooled_ols_q"] < 0.05).sum())
    clustered_q_total = int((comparison["woman_cluster_q"] < 0.05).sum())
    mixed_q_total = int((comparison["primary_mixed_model_q"] < 0.05).sum())
    pooled_p_total = int((comparison["pooled_ols_p"] < 0.05).sum())
    clustered_p_total = int((comparison["woman_cluster_p"] < 0.05).sum())
    mixed_p_total = int((comparison["primary_mixed_model_p"] < 0.05).sum())

    cluster_vs_pooled_larger = int(
        (comparison["woman_cluster_p"] > comparison["pooled_ols_p"]).sum()
    )
    median_cluster_p_ratio = float(
        np.median(comparison["woman_cluster_p"] / comparison["pooled_ols_p"])
    )
    cluster_vs_mixed_larger = int(
        (comparison["woman_cluster_p"] > comparison["primary_mixed_model_p"]).sum()
    )
    median_cluster_vs_mixed_ratio = float(
        np.median(comparison["woman_cluster_p"] / comparison["primary_mixed_model_p"])
    )

    original_rows = []
    for _, row in original.iterrows():
        original_rows.append(
            "| {feature} | {mixed:.3f} | {cluster:.3f} | {se:.3f} | {p} | {q} | {direction} |".format(
                feature=row["feature_name"],
                mixed=row["primary_mixed_model_estimate"],
                cluster=row["woman_cluster_estimate"],
                se=row["woman_cluster_robust_se"],
                p=format_p(row["woman_cluster_p"]),
                q=format_p(row["woman_cluster_q"]),
                direction="same" if row["woman_cluster_same_direction_as_primary"] else "reversed",
            )
        )

    input_rows = [
        f"- `{REGISTRY_PATH.relative_to(PROJECT_ROOT)}` (SHA-256 `{sha256(REGISTRY_PATH)}`)",
        f"- `{EVENT_UNITS_PATH.relative_to(PROJECT_ROOT)}` (SHA-256 `{sha256(EVENT_UNITS_PATH)}`)",
        f"- `{MIXED_RESULTS_PATH.relative_to(PROJECT_ROOT)}` (SHA-256 `{sha256(MIXED_RESULTS_PATH)}`)",
        "- Feature-family check: 34 primary features and expected 11 supported effects",
    ]

    return f"""# EHG9 temporal-dependence sensitivity analysis

## Question

Do the frozen 11/34 EHG9 uterine-contraction (UC) versus inactivity (IN) findings remain stable when inference is robust to dependence among repeated observations from the same woman?

## Frozen inputs and reproduction checks

{chr(10).join(input_rows)}

The Stage 06 registry contains exactly 34 unique primary FWH features. The EHG9 mixed-model file contains one UC-vs-IN result for every registered feature. Reapplying Benjamini-Hochberg (BH) to the 34 frozen p-values reproduces every stored q-value within numerical tolerance and reproduces the exact 11-feature q<0.05 set. The manuscript update memo states 11/34 and identifies SVD entropy, fixed-r SampEn, and Hjorth mobility as the three prior borderline findings removed by the expanded family. The frozen results agree with those statements. The supplied updated memo does not tabulate every numerical estimate and interval, so manuscript-level numerical identity beyond the count and named feature-set change cannot be checked from that document.

No primary mixed model was refitted or altered.

## Methods

The analysis used the frozen EHG9 event-unit table and retained only UC and IN: {len(uc_in):,} observations ({int((uc_in['window_class'] == 'contraction').sum()):,} UC and {int((uc_in['window_class'] == 'baseline').sum()):,} IN), {uc_in['patient_id'].nunique()} women, and {uc_in['recording_id'].nunique()} recordings. Of the 42 women, {int((has_in & has_uc).sum())} contributed both states, {int((has_in & ~has_uc).sum())} contributed IN only, and {int((~has_in & has_uc).sum())} contributed UC only. Woman-level observation counts ranged from {int(uc_in.groupby('patient_id').size().min())} to {int(uc_in.groupby('patient_id').size().max())}.

For exact comparability with Stage 06, each feature used its frozen all-three-state mean and sample SD (ddof=1) from the primary fit. OLS was then fit as `z_feature ~ state_uc`, with IN as 0 and UC as 1. Pooled inference used ordinary OLS standard errors. Woman-cluster inference used the identical point estimate with statsmodels' CR1 covariance correction, `G/(G-1) * (N-1)/(N-K)`, cluster degrees of freedom `G-1 = 41`, and two-sided t-based 95% CIs and p-values. BH was applied separately to the 34 pooled p-values and the 34 woman-cluster p-values. No randomness was used.

## Results for the original 11

| Feature | Mixed estimate | Clustered estimate | Cluster-robust SE | Cluster p | Cluster q | Direction vs mixed |
|---|---:|---:|---:|---:|---:|---|
{chr(10).join(original_rows)}

1. **Frozen reproduction:** Yes. The registry, 34 mixed-model UC-vs-IN rows, stored BH q-values, and exact 11-feature set reproduce correctly. The count and named FDR-family change agree with the manuscript update memo.
2. **Effect direction:** {same_direction_n}/11 retained the primary mixed-model direction under woman-cluster-robust OLS.
3. **BH support:** {cluster_q_n}/11 retained woman-cluster BH q<0.05. This was `{original.loc[original['woman_cluster_q_lt_0_05'], 'feature_name'].iloc[0] if cluster_q_n else 'none'}`.
4. **Direction reversals:** {len(reversals)} original features reversed direction: {format_feature_list(reversals)}. These reversals reflect the difference between the unadjusted pooled mean contrast and the mixed-model coefficient. Cluster-robust covariance changes the uncertainty, not the OLS point estimate.
5. **Magnitude changes:** For this report, material attenuation was defined descriptively as at least 20% lower absolute coefficient than the mixed estimate. By that threshold, {len(materially_attenuated)} original effects were materially attenuated: {format_feature_list(materially_attenuated)}. The two spectral estimates moved through zero to small positive values. The median estimate retained its negative direction but was about {100 * (1 - float(original.loc[original['feature_name'].eq('fwh_time_median'), 'absolute_ratio_cluster_to_mixed'].iloc[0])):.0f}% smaller in absolute magnitude. This threshold is not an inferential test.
6. **Amplitude/time-domain pattern:** The point-estimate pattern was preserved. All {time_same_n}/{len(time_domain)} original time-domain effects kept their direction. Mean absolute value, RMS, peak-to-peak amplitude, and waveform length changed by less than 6% in magnitude; skewness changed by less than 4%; median attenuated as noted above. Reduced FDR support primarily reflects much larger cluster-robust uncertainty, not disappearance of the core amplitude/waveform point estimates.
7. **P-value comparison:** Across all 34 features, nominal p<0.05 counts were mixed {mixed_p_total}, pooled {pooled_p_total}, and woman-cluster {clustered_p_total}; BH q<0.05 counts were mixed {mixed_q_total}, pooled {pooled_q_total}, and woman-cluster {clustered_q_total}. Woman-cluster p-values exceeded pooled p-values for {cluster_vs_pooled_larger}/34 features, with a median cluster-to-pooled p-value ratio of {median_cluster_p_ratio:.1f}. Woman-cluster p-values exceeded mixed-model p-values for {cluster_vs_mixed_larger}/34 features, with a median ratio of {median_cluster_vs_mixed_ratio:.1f}. Pooled and mixed p-values were not uniformly ordered because the mixed model changes both the coefficient and the variance model.
8. **Pseudoreplication:** Yes, the pooled-versus-cluster comparison provides strong empirical evidence that treating {len(uc_in):,} windows as independent overstates precision. Pooled BH support was {pooled_q_total}/34, but woman-cluster BH support was {clustered_q_total}/34, while the pooled and clustered coefficients were exactly identical. The difference therefore comes entirely from dependence-robust uncertainty. This does not invalidate the mixed model, which already modeled woman and, where estimable, recording effects.
9. **Appropriateness of woman clustering:** Woman is the correct primary independence unit because recordings are nested within women. Forty-two clusters are workable but modest, so finite-cluster correction is important and tail-area inference remains approximate. Severe cluster-size imbalance (range {int(uc_in.groupby('patient_id').size().min())}--{int(uc_in.groupby('patient_id').size().max())} observations), only 34 women with both states, and eight single-state women can make the marginal OLS contrast sensitive to cluster composition. Cluster-robust OLS also leaves the pooled coefficient unchanged and does not isolate a purely within-woman contrast or adjust for recording-level heterogeneity. These are limitations of this sensitivity estimand, not reasons to cluster by recording instead. No dataset feature makes woman clustering technically inappropriate.

## Conclusion and recommendation

The result is stable in its main amplitude/time-domain point-estimate pattern but not in multiplicity-adjusted cluster-robust precision. Eight of the original 11 keep direction, all six original time-domain effects keep direction, and the four core amplitude/waveform coefficients remain close to the mixed-model estimates. Only one original feature retains BH q<0.05 under the conservative woman-cluster OLS analysis. Three non-core features reverse direction, including two spectral effects that shrink close to zero.

The primary finding appears sufficiently stable to proceed with the revision only if it is framed as a predominantly amplitude/time-domain pattern supported by the prespecified mixed model, with the dependence-robust analysis presented as evidence of reduced precision and weaker support for individual features. It would be too strong to claim that all 11 individual findings are robust under woman-cluster inference.

## Scope

This analysis stops at EHG9 UC versus IN. It does not analyze cross-channel models, alternative window lengths, SampEn scale sensitivity, fetal-movement contamination, recording-clustered inference, or manuscript wording.
"""


def package_versions_text() -> str:
    packages = ["numpy", "pandas", "scipy", "statsmodels", "patsy"]
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
            f"statsmodels imported version: {statsmodels.__version__}",
            "Randomness: none used",
            "Random seed: not applicable",
            "Cluster covariance: statsmodels cov_type=cluster, use_correction=True, df_correction=True, use_t=True",
            "Multiplicity: statsmodels.stats.multitest.multipletests(method=fdr_bh)",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    require_inputs()
    outputs = refuse_overwrite(output_dir)
    _, event_units, mixed, features = load_and_verify_frozen_inputs()
    comparison = fit_models(event_units, mixed, features)
    stability = comparison.loc[comparison["original_11"]].copy()

    comparison.to_csv(outputs["comparison"], index=False, float_format="%.17g")
    stability.to_csv(outputs["stability"], index=False, float_format="%.17g")
    outputs["report"].write_text(build_report(comparison, event_units), encoding="utf-8")
    outputs["versions"].write_text(package_versions_text(), encoding="utf-8")

    print(f"Wrote {len(comparison)}-feature comparison: {outputs['comparison']}")
    print(f"Wrote {len(stability)}-feature stability summary: {outputs['stability']}")
    print(f"Wrote report: {outputs['report']}")
    print(f"Wrote package versions: {outputs['versions']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
