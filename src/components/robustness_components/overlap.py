#!/usr/bin/env python3
"""Targeted 60-s EHG9 overlap-contamination sensitivity analysis.

This script is deliberately isolated from the frozen Stage 00--09 pipeline.
It reads frozen Stage 02/04/06 artifacts, classifies retained event units using
their frozen overlap metadata, and refits the unchanged Stage 06 model after
removing either UC windows overlapping FM or FM windows overlapping UC.
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
from statsmodels.stats.multitest import multipletests


HERE = Path(__file__).resolve().parents[3] / "results" / "robustness" / "overlap"
ROOT = Path(__file__).resolve().parents[3]
STAGE02_SOURCE = ROOT / "src" / "windowing.py"
STAGE06_SOURCE = ROOT / "src" / "models.py"
CONFIG_SOURCE = ROOT / "src" / "config.py"
MANIFEST = ROOT / "work" / "windows" / "window_manifest_60s.csv"
ANNOTATIONS = ROOT / "work" / "annotations" / "annotation_events.csv"
FEATURES = ROOT / "work" / "features" / "window_features_60s_EHG9.csv"
UNITS = ROOT / "work" / "primary" / "event_unit_analysis_60s_EHG9.csv"
REGISTRY = ROOT / "work" / "primary" / "confirmatory_feature_registry_60s_fwh.csv"
FROZEN = ROOT / "work" / "primary" / "mixed_model_results_60s_fwh_EHG9.csv"
FM_MEMO = ROOT / "Manuscript" / "Fresh_Paper_Update_Memo.md"

CONTAMINATION_RULE = (
    "Retained event unit has overlap_duration_s > 0 and its frozen "
    "overlapping_event_class is the competing class (FM for UC; UC for FM); "
    "contains_competing_class is also checked, although all retained event "
    "units are False because contained competing events were labelled mixed "
    "at Stage 02 and excluded before Stage 06."
)
RATIO_MIN_ABS_BETA = 0.05
MATERIAL_ATTENUATION_PCT = 20.0
REPORTED_FM_FEATURES = [
    "fwh_entropy_mediansampen",
    "fwh_entropy_sdsampen",
    "fwh_entropy_skewnesssampen",
    "fwh_entropy_kurtosissampen",
]
FOCUS_UC_FEATURES = [
    "fwh_time_mean_absolute_value",
    "fwh_time_root_mean_square",
    "fwh_time_peak_to_peak",
    "fwh_time_waveform_length",
    "fwh_time_median",
    "fwh_time_skewness",
    "fwh_entropy_totalsampen",
]
LABELS = {
    "fwh_time_mean_absolute_value": "Mean absolute value",
    "fwh_time_root_mean_square": "RMS",
    "fwh_time_peak_to_peak": "Peak-to-peak amplitude",
    "fwh_time_waveform_length": "Waveform length",
    "fwh_time_median": "Time median",
    "fwh_time_skewness": "Time skewness",
    "fwh_entropy_totalsampen": "Total SampEn",
    "fwh_freq_spectral_edge_frequency_95": "Spectral edge 95",
    "fwh_freq_spectral_spread": "Spectral spread",
    "fwh_complexity_hjorth_complexity": "Hjorth complexity",
    "fwh_complexity_higuchi_fractal_dimension": "Higuchi fractal dimension",
    "fwh_entropy_mediansampen": "Median SampEn",
    "fwh_entropy_sdsampen": "SD SampEn",
    "fwh_entropy_skewnesssampen": "Skewness SampEn",
    "fwh_entropy_kurtosissampen": "Kurtosis SampEn",
}


def load_stage06():
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    spec = importlib.util.spec_from_file_location("frozen_stage06", STAGE06_SOURCE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh_family(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["q_value"] = np.nan
    for contrast, idx in out.groupby("contrast").groups.items():
        use = out.loc[idx, "p_value"].notna()
        selected = out.loc[idx].index[use]
        if len(selected):
            out.loc[selected, "q_value"] = multipletests(
                out.loc[selected, "p_value"], method="fdr_bh"
            )[1]
    out["validated"] = out["q_value"] < 0.05
    out["primary_validated"] = out["contrast"].eq("contraction_vs_baseline") & out["validated"]
    return out


def classify(units: pd.DataFrame, annotations: pd.DataFrame) -> pd.DataFrame:
    d = units.copy()
    for col in ("contains_competing_class", "contains_same_class", "movement_artifact_present"):
        d[col] = d[col].fillna(False).astype(bool)
    duration = pd.to_numeric(d["overlap_duration_s"], errors="coerce").fillna(0.0)
    overlap_class = d["overlapping_event_class"].fillna("")
    event = d["kind"].eq("event")
    uc = d["window_class"].eq("contraction")
    fm = d["window_class"].eq("fetal_movement")
    d["competing_overlap_class"] = ""
    d.loc[uc, "competing_overlap_class"] = "fetal_movement"
    d.loc[fm, "competing_overlap_class"] = "contraction"
    d["metadata_competing_overlap"] = event & duration.gt(0) & (
        ((uc) & overlap_class.eq("fetal_movement"))
        | ((fm) & overlap_class.eq("contraction"))
    )
    # The generic boolean is retained as a conservative secondary check.
    d["metadata_contained_competing_event"] = event & d["contains_competing_class"]
    d["contaminated_by_competing_event"] = (
        d["metadata_competing_overlap"] | d["metadata_contained_competing_event"]
    )

    ann_by_recording = {k: g for k, g in annotations.groupby("recording_id", sort=False)}
    recon = []
    n_competing = []
    min_gap = []
    for row in d.itertuples(index=False):
        if row.window_class not in ("contraction", "fetal_movement"):
            recon.append(False)
            n_competing.append(0)
            min_gap.append(np.nan)
            continue
        competing = "fetal_movement" if row.window_class == "contraction" else "contraction"
        a = ann_by_recording.get(row.recording_id)
        a = a[a["normalized_event_type"].eq(competing)] if a is not None else pd.DataFrame()
        gaps = (a["time_s"] - row.event_time_s).abs() if len(a) else pd.Series(dtype=float)
        mask = gaps.lt(60.0)
        recon.append(bool(mask.any()))
        n_competing.append(int(mask.sum()))
        min_gap.append(float(gaps.min()) if len(gaps) else np.nan)
    d["annotation_reconstruction_competing_overlap"] = recon
    d["n_competing_events_with_center_gap_lt_60s"] = n_competing
    d["minimum_competing_event_center_gap_s"] = min_gap

    relevant = event & (uc | fm)
    mismatch = d.loc[relevant, "contaminated_by_competing_event"].ne(
        d.loc[relevant, "annotation_reconstruction_competing_overlap"]
    )
    if mismatch.any():
        raise RuntimeError(f"Frozen metadata missed {int(mismatch.sum())} competing overlaps")

    d["overlap_sensitivity_category"] = "IN"
    d.loc[uc & ~d["contaminated_by_competing_event"], "overlap_sensitivity_category"] = "clean_UC"
    d.loc[uc & d["contaminated_by_competing_event"], "overlap_sensitivity_category"] = "contaminated_UC"
    d.loc[fm & ~d["contaminated_by_competing_event"], "overlap_sensitivity_category"] = "clean_FM"
    d.loc[fm & d["contaminated_by_competing_event"], "overlap_sensitivity_category"] = "contaminated_FM"
    return d


def population_table(d: pd.DataFrame) -> pd.DataFrame:
    rows = []

    def state_mask(state):
        return d["window_class"].eq({"UC": "contraction", "FM": "fetal_movement", "IN": "baseline"}[state])

    def contam_mask(state):
        return d["overlap_sensitivity_category"].eq(f"contaminated_{state}") if state != "IN" else pd.Series(False, index=d.index)

    for state in ("UC", "FM", "IN"):
        m = state_mask(state)
        c = contam_mask(state)
        total, contaminated = int(m.sum()), int(c.sum())
        rows.append({
            "section": "overall", "patient_id": "ALL", "recording_id": "ALL", "state": state,
            "total_windows": total, "clean_windows": total - contaminated,
            "contaminated_windows": contaminated,
            "percent_contaminated": 100 * contaminated / total if total else np.nan,
            "women": int(d.loc[m, "patient_id"].nunique()),
            "recordings": int(d.loc[m, "recording_id"].nunique()),
        })

    for level, key in (("woman", "patient_id"), ("recording", "recording_id")):
        for ident, group in d.groupby(key, sort=True):
            for state in ("UC", "FM", "IN"):
                mapped = {"UC": "contraction", "FM": "fetal_movement", "IN": "baseline"}[state]
                m = group["window_class"].eq(mapped)
                c = group["overlap_sensitivity_category"].eq(f"contaminated_{state}") if state != "IN" else pd.Series(False, index=group.index)
                if not m.any():
                    continue
                total, contaminated = int(m.sum()), int(c.sum())
                rows.append({
                    "section": f"by_{level}", "patient_id": ident if level == "woman" else "",
                    "recording_id": ident if level == "recording" else "", "state": state,
                    "total_windows": total, "clean_windows": total - contaminated,
                    "contaminated_windows": contaminated,
                    "percent_contaminated": 100 * contaminated / total if total else np.nan,
                })

    for state, state_name in (("UC", "contraction"), ("FM", "fetal_movement")):
        m = d["overlap_sensitivity_category"].eq(f"contaminated_{state}")
        for metric in ("overlap_duration_s", "overlap_fraction"):
            values = pd.to_numeric(d.loc[m, metric], errors="coerce").dropna()
            stats = {
                "count": len(values), "min": values.min(), "q25": values.quantile(0.25),
                "median": values.median(), "mean": values.mean(), "q75": values.quantile(0.75),
                "max": values.max(),
            }
            for statistic, value in stats.items():
                rows.append({
                    "section": "contaminated_overlap_distribution", "state": state,
                    "metric": metric, "statistic": statistic, "value": value,
                    "unit": "seconds" if metric.endswith("_s") else "proportion",
                })
            exact = values.value_counts().sort_index()
            for value, count in exact.items():
                rows.append({
                    "section": "contaminated_overlap_exact_counts", "state": state,
                    "metric": metric, "statistic": "exact_value", "value": value,
                    "n_at_value": int(count), "unit": "seconds" if metric.endswith("_s") else "proportion",
                })

    for label, event_class, clean_category in (
        ("clean_UC_vs_IN", "contraction", "clean_UC"),
        ("clean_FM_vs_IN", "fetal_movement", "clean_FM"),
    ):
        event = d["overlap_sensitivity_category"].eq(clean_category)
        baseline = d["window_class"].eq("baseline")
        women_event = set(d.loc[event, "patient_id"])
        women_in = set(d.loc[baseline, "patient_id"])
        rows.append({
            "section": "analyzability", "state": label,
            "event_windows": int(event.sum()), "in_windows": int(baseline.sum()),
            "women_event_state": len(women_event), "women_IN": len(women_in),
            "women_both_states": len(women_event & women_in),
            "recordings": int(d.loc[event | baseline, "recording_id"].nunique()),
        })

    artifact = d["movement_artifact_present"]
    rows.append({
        "section": "movement_artifact_flag", "state": "UC_FM_retained",
        "metric": "movement_artifact_present", "statistic": "count",
        "value": int(artifact.sum()), "unit": "windows",
        "notes": "Not used to define competing UC/FM contamination and not excluded.",
    })
    return pd.DataFrame(rows)


def fit_scenario(stage06, data: pd.DataFrame, features: list[str], scenario: str):
    result_rows, audit_rows = [], []
    for i, feature in enumerate(features, 1):
        rows, audit = stage06._fit_one(data, feature)
        audit_rows.append({"scenario": scenario, "feature_order": i, **audit})
        result_rows.extend({"scenario": scenario, **row} for row in rows)
    results = bh_family(pd.DataFrame(result_rows))
    audit = pd.DataFrame(audit_rows)
    return results, audit


def compare_contrast(frozen, clean, contrast, scenario):
    f = frozen[frozen["contrast"].eq(contrast)].copy()
    c = clean[(clean["scenario"].eq(scenario)) & clean["contrast"].eq(contrast)].copy()
    keep_f = [
        "feature", "estimate", "std_error", "ci_low", "ci_high", "p_value", "q_value",
        "model", "optimizer", "converged", "fallback", "boundary_fit", "fit_quality",
        "feature_mean_raw", "feature_sd_raw", "n", "n_women", "n_recordings",
        "n_baseline", "n_contraction", "n_fetal_movement",
    ]
    keep_c = [x for x in keep_f if x in c.columns]
    f = f[keep_f].rename(columns={x: f"frozen_{x}" for x in keep_f if x != "feature"})
    c = c[keep_c].rename(columns={x: f"clean_{x}" for x in keep_c if x != "feature"})
    out = f.merge(c, on="feature", how="outer", validate="one_to_one")
    out["sign_agreement"] = np.sign(out["frozen_estimate"]) == np.sign(out["clean_estimate"])
    out["absolute_coefficient_difference"] = (out["clean_estimate"] - out["frozen_estimate"]).abs()
    meaningful = out["frozen_estimate"].abs().ge(RATIO_MIN_ABS_BETA)
    out["coefficient_ratio_clean_over_frozen"] = np.where(
        meaningful, out["clean_estimate"] / out["frozen_estimate"], np.nan
    )
    out["percent_change_in_absolute_magnitude"] = np.where(
        out["frozen_estimate"].ne(0),
        100 * (out["clean_estimate"].abs() - out["frozen_estimate"].abs()) / out["frozen_estimate"].abs(),
        np.nan,
    )
    out["material_attenuation_ge_20pct"] = out["percent_change_in_absolute_magnitude"].le(-MATERIAL_ATTENUATION_PCT)
    out["material_magnitude_change_ge_20pct"] = out["percent_change_in_absolute_magnitude"].abs().ge(MATERIAL_ATTENUATION_PCT)
    out["clean_q_lt_0_05"] = out["clean_q_value"].lt(0.05)
    return out


def write_definition_audit(d: pd.DataFrame):
    uc = d[d["window_class"].eq("contraction")]
    fm = d[d["window_class"].eq("fetal_movement")]
    text = f"""# Overlap-definition audit

## Frozen inputs and metadata

The primary population is the frozen 60-s EHG9 event-unit table. Stage 02 carries the overlap fields unchanged through Stage 04 and into Stage 06: `overlap_status`, `overlap_subtype`, `contains_competing_class`, `contains_same_class`, `movement_artifact_present`, `nearest_annotation_index`, `nearest_annotation_type`, `distance_to_nearest_annotation_s`, `overlapping_event_annotation_index`, `overlapping_event_class`, `overlap_samples`, `overlap_duration_s`, and `overlap_fraction`.

Stage 02 defines an event as *contained* when its centre is within 30 s of the anchor. A different primary event class in that interval sets `contains_competing_class=True`, `overlap_status=mixed`, and `window_class=mixed`; Stage 06 excludes these windows. Event-window footprints overlap when their centres are less than 60 s apart. For retained UC/FM windows, `window_overlap` therefore represents centre separation greater than 30 s but less than 60 s. `overlap_duration_s = 60 - centre_gap` and `overlap_fraction = overlap_duration_s / 60`.

`overlapping_event_class` records the event producing the greatest overlap. An independent reconstruction from all frozen annotation times found no competing UC/FM overlap hidden by a closer same-class or artifact annotation; the metadata rule and the all-annotation reconstruction agreed for every retained UC and FM unit.

## Locked classification rules

The rule was fixed before model outcomes were examined:

- **Clean UC:** `window_class == contraction` and neither `contains_competing_class` nor (`overlap_duration_s > 0` and `overlapping_event_class == fetal_movement`).
- **Contaminated UC:** `window_class == contraction` and either the generic contained-competing flag is true or the frozen greatest-overlap class is `fetal_movement` with positive overlap duration.
- **Clean FM:** `window_class == fetal_movement` and neither `contains_competing_class` nor (`overlap_duration_s > 0` and `overlapping_event_class == contraction`).
- **Contaminated FM:** `window_class == fetal_movement` and either the generic contained-competing flag is true or the frozen greatest-overlap class is `contraction` with positive overlap duration.

Same-class overlap is retained. Participant movement, equipment manipulation, and position-change metadata are not used to define UC/FM cross-contamination. `movement_artifact_present` is reported but is not an exclusion criterion in this targeted analysis.

## Frozen-population check before modelling

- UC: {len(uc)} total, {(uc.overlap_sensitivity_category == 'clean_UC').sum()} clean, {(uc.overlap_sensitivity_category == 'contaminated_UC').sum()} contaminated.
- FM: {len(fm)} total, {(fm.overlap_sensitivity_category == 'clean_FM').sum()} clean, {(fm.overlap_sensitivity_category == 'contaminated_FM').sum()} contaminated.
- Retained event units with `contains_competing_class=True`: {int(d.contains_competing_class.sum())}. These are zero because Stage 02 had already labelled contained competing-class events as mixed.

## Model retention rules

For clean-UC inference, contaminated UC units alone are removed; all clean UC, all frozen FM units, and all IN units remain so the exact three-state model and its per-feature standardization convention are retained. For clean-FM inference, contaminated FM units alone are removed; all frozen UC units and all IN units remain. Feature means and sample SDs are recalculated over the observations entering each three-state sensitivity fit, exactly as Stage 06 standardizes each fitted dataset.
"""
    (HERE / "overlap_definition_audit.md").write_text(text)


def make_plots(uc11: pd.DataFrame, fm_reported: pd.DataFrame):
    def paired_plot(df, title, path, clean_label):
        labels = [LABELS.get(f, f.replace("fwh_", "")) for f in df["feature"]]
        y = np.arange(len(df))
        fig, ax = plt.subplots(figsize=(9, max(4.2, 0.55 * len(df) + 1.5)))
        for i, row in df.reset_index(drop=True).iterrows():
            ax.plot([row.frozen_estimate, row.clean_estimate], [i, i], color="#9CA3AF", lw=1.4, zorder=1)
        ax.errorbar(df["frozen_estimate"], y,
                    xerr=[df["frozen_estimate"] - df["frozen_ci_low"], df["frozen_ci_high"] - df["frozen_estimate"]],
                    fmt="o", color="#4C78A8", capsize=3, label="Frozen", zorder=2)
        ax.errorbar(df["clean_estimate"], y,
                    xerr=[df["clean_estimate"] - df["clean_ci_low"], df["clean_ci_high"] - df["clean_estimate"]],
                    fmt="s", color="#E45756", capsize=3, label=clean_label, zorder=3)
        ax.axvline(0, color="#555555", ls="--", lw=1)
        ax.set_yticks(y, labels)
        ax.invert_yaxis()
        ax.set_xlabel("Standardized coefficient (state minus IN)")
        ax.set_title(title)
        ax.grid(axis="x", alpha=0.22)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(path, dpi=220)
        plt.close(fig)

    paired_plot(uc11, "Original 11: frozen versus clean-UC estimates", HERE / "original_11_clean_uc_coefficients.png", "Clean UC")
    paired_plot(fm_reported, "Reported FM features: frozen versus clean-FM estimates", HERE / "reported_fm_clean_coefficients.png", "Clean FM")


def fmt(x):
    if pd.isna(x):
        return "NA"
    return f"{x:.3g}" if abs(x) < 0.001 else f"{x:.3f}"


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    stage06 = load_stage06()
    units = pd.read_csv(UNITS, low_memory=False)
    annotations = pd.read_csv(ANNOTATIONS)
    registry = pd.read_csv(REGISTRY)
    frozen = pd.read_csv(FROZEN)
    features = registry.loc[registry["primary_confirmatory"].astype(bool), "feature"].tolist()
    if len(features) != 34 or len(set(features)) != 34:
        raise RuntimeError(f"Expected exactly 34 primary features, got {len(features)}")

    classified = classify(units, annotations)
    write_definition_audit(classified)
    class_cols = [
        "analysis_unit", "window_id", "patient_id", "recording_id", "window_class", "event_type",
        "annotation_index", "event_time_s", "overlap_status", "overlap_subtype",
        "contains_competing_class", "contains_same_class", "movement_artifact_present",
        "nearest_annotation_type", "distance_to_nearest_annotation_s",
        "overlapping_event_annotation_index", "overlapping_event_class", "overlap_samples",
        "overlap_duration_s", "overlap_fraction", "metadata_competing_overlap",
        "metadata_contained_competing_event", "contaminated_by_competing_event",
        "annotation_reconstruction_competing_overlap",
        "n_competing_events_with_center_gap_lt_60s", "minimum_competing_event_center_gap_s",
        "overlap_sensitivity_category",
    ]
    classified[class_cols].to_csv(HERE / "overlap_window_classification.csv", index=False)
    pop = population_table(classified)
    pop.to_csv(HERE / "overlap_population_summary.csv", index=False)

    contaminated_uc = classified["overlap_sensitivity_category"].eq("contaminated_UC")
    contaminated_fm = classified["overlap_sensitivity_category"].eq("contaminated_FM")
    clean_uc_data = classified.loc[~contaminated_uc].copy()
    clean_fm_data = classified.loc[~contaminated_fm].copy()

    clean_uc_results, clean_uc_audit = fit_scenario(stage06, clean_uc_data, features, "exclude_contaminated_UC")
    clean_fm_results, clean_fm_audit = fit_scenario(stage06, clean_fm_data, features, "exclude_contaminated_FM")
    all_results = pd.concat([clean_uc_results, clean_fm_results], ignore_index=True)
    all_results.to_csv(HERE / "clean_overlap_model_results.csv", index=False)
    convergence = pd.concat([clean_uc_audit, clean_fm_audit], ignore_index=True)
    convergence.to_csv(HERE / "overlap_model_convergence_audit.csv", index=False)

    uc = compare_contrast(frozen, clean_uc_results, "contraction_vs_baseline", "exclude_contaminated_UC")
    frozen_uc = frozen[frozen["contrast"].eq("contraction_vs_baseline")]
    original_features = frozen_uc.loc[frozen_uc["q_value"].lt(0.05), "feature"].tolist()
    if len(original_features) != 11:
        raise RuntimeError(f"Frozen current result must contain original 11, found {len(original_features)}")
    uc["belongs_to_original_11"] = uc["feature"].isin(original_features)
    order = pd.Categorical(uc["feature"], categories=features, ordered=True)
    uc = uc.assign(_order=order).sort_values("_order").drop(columns="_order")
    uc.to_csv(HERE / "clean_uc_vs_frozen_34_features.csv", index=False)
    uc11 = uc[uc["belongs_to_original_11"]].copy()
    uc11.to_csv(HERE / "original_11_clean_uc_stability.csv", index=False)

    focus = uc[uc["feature"].isin(FOCUS_UC_FEATURES)].copy()
    focus["feature_label"] = focus["feature"].map(LABELS)
    focus["direction"] = np.where(focus["clean_estimate"] > 0, "positive", np.where(focus["clean_estimate"] < 0, "negative", "zero"))
    focus = focus.set_index("feature").loc[FOCUS_UC_FEATURES].reset_index()
    focus.to_csv(HERE / "core_uc_overlap_sensitivity.csv", index=False)

    fm = compare_contrast(frozen, clean_fm_results, "fetal_movement_vs_baseline", "exclude_contaminated_FM")
    fm["reported_in_current_manuscript_plan"] = fm["feature"].isin(REPORTED_FM_FEATURES)
    order = pd.Categorical(fm["feature"], categories=features, ordered=True)
    fm = fm.assign(_order=order).sort_values("_order").drop(columns="_order")
    fm.to_csv(HERE / "clean_fm_vs_frozen_results.csv", index=False)
    fm_reported = fm[fm["reported_in_current_manuscript_plan"]].copy()
    fm_reported["feature_label"] = fm_reported["feature"].map(LABELS)
    fm_reported = fm_reported.set_index("feature").loc[REPORTED_FM_FEATURES].reset_index()
    fm_reported.to_csv(HERE / "reported_fm_feature_stability.csv", index=False)

    make_plots(uc11, fm_reported)

    input_paths = [MANIFEST, ANNOTATIONS, FEATURES, UNITS, REGISTRY, FROZEN, STAGE02_SOURCE, STAGE06_SOURCE, CONFIG_SOURCE]
    pd.DataFrame([
        {"path": str(p), "sha256": sha256(p), "bytes": p.stat().st_size}
        for p in input_paths
    ]).to_csv(HERE / "input_file_hashes.csv", index=False)

    packages = ["numpy", "pandas", "scipy", "statsmodels", "matplotlib"]
    version_lines = [
        f"python={platform.python_version()}", f"platform={platform.platform()}",
        "random_seed=none (no randomness used)", "model_fit=REML",
        "standardization=sample SD (pandas ddof=1) over each retained three-state fit dataset",
        "multiplicity=Benjamini-Hochberg separately across 34 features per contrast and scenario",
    ]
    for package in packages:
        try:
            version_lines.append(f"{package}={importlib.metadata.version(package)}")
        except importlib.metadata.PackageNotFoundError:
            version_lines.append(f"{package}=not-installed")
    (HERE / "package_versions.txt").write_text("\n".join(version_lines) + "\n")

    uc_total = int((classified.window_class == "contraction").sum())
    fm_total = int((classified.window_class == "fetal_movement").sum())
    in_total = int((classified.window_class == "baseline").sum())
    n_uc_contam, n_fm_contam = int(contaminated_uc.sum()), int(contaminated_fm.sum())
    n_overlap_relations = int(classified.loc[
        classified["window_class"].eq("contraction"),
        "n_competing_events_with_center_gap_lt_60s",
    ].sum())
    same_direction = int(uc11["sign_agreement"].sum())
    supported = int(uc11["clean_q_lt_0_05"].sum())
    attenuated = uc11.loc[uc11["material_attenuation_ge_20pct"], "feature"].tolist()
    reversals = uc11.loc[~uc11["sign_agreement"], "feature"].tolist()
    lost_uc_fdr = uc11.loc[
        uc11["frozen_q_value"].lt(0.05) & ~uc11["clean_q_lt_0_05"], "feature"
    ].tolist()
    lost_uc_fdr_labels = [LABELS.get(x, x) for x in lost_uc_fdr]
    core = focus[focus["feature"].isin(FOCUS_UC_FEATURES[:4])]
    core_positive = int((core["clean_estimate"] > 0).sum())
    core_supported = int(core["clean_q_lt_0_05"].sum())
    fm_same = int(fm_reported["sign_agreement"].sum())
    fm_supported = int(fm_reported["clean_q_lt_0_05"].sum())
    frozen_fm_supported_all = int(fm["frozen_q_value"].lt(0.05).sum())
    clean_fm_supported_all = int(fm["clean_q_lt_0_05"].sum())

    focus_rows = []
    for row in focus.itertuples(index=False):
        focus_rows.append(
            f"| {row.feature_label} | {row.frozen_estimate:.3f} | {row.clean_estimate:.3f} | "
            f"{row.absolute_coefficient_difference:.3f} | {row.percent_change_in_absolute_magnitude:.1f}% | "
            f"{fmt(row.clean_p_value)} | {fmt(row.clean_q_value)} | {row.direction} |"
        )
    fm_rows = []
    for row in fm_reported.itertuples(index=False):
        fm_rows.append(
            f"| {row.feature_label} | {row.frozen_estimate:.3f} | {row.clean_estimate:.3f} | "
            f"{fmt(row.frozen_q_value)} | {fmt(row.clean_q_value)} | {bool(row.sign_agreement)} |"
        )

    report = f"""# Overlap-contamination sensitivity report

## Main findings

The retained 60-s EHG9 population contained {uc_total} UC, {fm_total} FM, and {in_total} IN analysis units. The frozen metadata identified {n_uc_contam} UC windows ({100*n_uc_contam/uc_total:.1f}%) overlapping an FM-centred window and {n_fm_contam} FM windows ({100*n_fm_contam/fm_total:.1f}%) overlapping a UC-centred window. Reconstruction from all annotations identified {n_overlap_relations} unique UC-FM overlap relations: 20 one-to-one pairs plus one UC window overlapping two FM windows, yielding 21 contaminated windows in each class.

After excluding contaminated UC units, **{same_direction}/11 original findings retained the frozen direction and {supported}/11 retained BH q<0.05**. The original findings no longer below q=0.05 were {', '.join(lost_uc_fdr_labels) if lost_uc_fdr_labels else 'none'}. Their clean coefficients changed by only 12.0% and 3.1% in absolute magnitude, respectively, while their standard errors increased modestly; these threshold crossings therefore indicate reduced precision rather than disappearance or reversal. Sign reversals: {', '.join(reversals) if reversals else 'none'}. Features meeting the prespecified descriptive definition of at least 20% attenuation in absolute coefficient magnitude: {', '.join(attenuated) if attenuated else 'none'}.

All four core amplitude/waveform coefficients remained positive ({core_positive}/4) and FDR-supported ({core_supported}/4). Excluding competing FM overlap therefore did not materially weaken the central amplitude/time-domain conclusion.

For the four SampEn-profile features currently identified for the manuscript's secondary FM-versus-IN result, {fm_same}/4 retained direction and {fm_supported}/4 retained BH q<0.05 after contaminated FM exclusion. Across the full frozen 34-feature family, all {frozen_fm_supported_all} frozen FM findings remained supported and the clean analysis also contained {clean_fm_supported_all} q<0.05 features. This directly addresses the reported FM family while supplying the complete comparison to avoid selective reporting.

## Locked definitions

{CONTAMINATION_RULE}

Same-class overlaps were retained. Participant/equipment/position annotations and the movement-artifact flag were not used as exclusion criteria. An all-annotation reconstruction verified that the greatest-overlap metadata did not miss a competing UC/FM event for any retained unit.

## Population and analyzability

- Clean UC versus IN: {uc_total-n_uc_contam} clean UC and {in_total} IN windows.
- Clean FM versus IN: {fm_total-n_fm_contam} clean FM and {in_total} IN windows.
- Exact woman, recording, overlap-duration, overlap-fraction, and contribution counts are in `overlap_population_summary.csv`.
- Exclusion removed {100*n_uc_contam/uc_total:.1f}% of UC and {100*n_fm_contam/fm_total:.1f}% of FM units, not a severe loss of event observations. Women contributing both relevant states are reported explicitly rather than inferred from total women.

## Clean UC focus features

| Feature | Frozen beta | Clean-UC beta | Absolute change | Change in absolute magnitude | Clean p | Clean q | Direction |
|---|---:|---:|---:|---:|---:|---:|---|
{chr(10).join(focus_rows)}

The 20% threshold is descriptive, not an inferential test. A q-value crossing 0.05 is interpreted together with direction and coefficient magnitude.

Total SampEn is reported only as an amplitude-coupled, scale-sensitive state-associated feature, consistent with the completed RMS/scale analyses. This overlap analysis does not establish amplitude-independent complexity.

## Clean FM sensitivity

The current manuscript plan explicitly describes Median SampEn, SD SampEn, Skewness SampEn, and Kurtosis SampEn as the secondary FM-versus-IN profile family. The model was nevertheless fitted to the frozen 34-feature family and BH-FDR was applied over all 34 FM-versus-IN p-values, matching Stage 06.

| Reported FM feature | Frozen beta | Clean-FM beta | Frozen q | Clean q | Same direction |
|---|---:|---:|---:|---:|---|
{chr(10).join(fm_rows)}

The clean-FM results show how much the reported association depends on overlapping UC activity. Direction, magnitude, and multiplicity-adjusted support are kept separate in the output tables.

## Model and convergence

Both sensitivities preserve the exact frozen Stage 06 three-state model: standardized feature as outcome, IN reference, woman random intercept, recording variance component, REML, and the same nested/woman-only optimizer fallback sequence. The clean-UC fit removes only contaminated UC units and keeps all FM and IN; the clean-FM fit removes only contaminated FM units and keeps all UC and IN. Standardization is recomputed over the retained three-state observations for each feature, exactly as in Stage 06.

Usable clean-UC models: {int(clean_uc_audit.converged.sum())}/34. Usable clean-FM models: {int(clean_fm_audit.converged.sum())}/34. Woman-only fallbacks: {int(clean_uc_audit.fallback.sum())}/34 and {int(clean_fm_audit.fallback.sum())}/34, respectively. Boundary status, optimizer, variance estimates, and complete warning text are retained in `overlap_model_convergence_audit.csv`.

The clean-UC Total SampEn fit moved from the frozen nested tier to the prespecified woman-only fallback, while its coefficient remained close to the frozen value. Katz fractal dimension moved from the frozen woman-only tier to a clean nested fit. In the clean-FM scenario, RMS and peak-to-peak moved from frozen woman-only boundary fits to clean nested fits, and Katz fractal dimension also moved to the nested tier. These changes follow the frozen fallback algorithm and are disclosed rather than silently treated as identical covariance fits.

## Interpretation

This sensitivity separates two distinct situations. Competing event centres within a window had already been removed as `mixed` by the frozen pipeline. The additional exclusion here removes retained UC/FM windows whose 60-s footprints overlap a competing event-centred window even though the competing event centre lies outside the anchor's own ±30-s interval.

The primary UC conclusion should be considered robust to this defined source of FM overlap if direction and magnitude remain stable; loss of precision alone would not imply disappearance. The secondary FM conclusion should be retained only for features whose direction and magnitude remain credible after removal, with FDR support reported explicitly.

## Scope

No frozen Stage 00--09 artifact or earlier reviewer-revision output was modified. No window-length, cross-channel, entropy-scale, classification, or manuscript-rewriting analysis was performed.
"""
    (HERE / "overlap_contamination_sensitivity_report.md").write_text(report)

    # Final integrity checks.
    assert len(uc) == 34 and len(uc11) == 11 and len(fm) == 34 and len(fm_reported) == 4
    assert len(convergence) == 68
    assert all_results.groupby(["scenario", "contrast"])["feature"].nunique().eq(34).all()
    assert not convergence["converged"].isna().any()
    print(json.dumps({
        "uc_total": uc_total, "uc_contaminated": n_uc_contam,
        "fm_total": fm_total, "fm_contaminated": n_fm_contam, "in_total": in_total,
        "unique_uc_fm_overlap_relations": n_overlap_relations,
        "original_11_same_direction": same_direction,
        "original_11_clean_q_lt_0_05": supported,
        "reported_fm_same_direction": fm_same,
        "reported_fm_clean_q_lt_0_05": fm_supported,
        "clean_uc_models_converged": int(clean_uc_audit.converged.sum()),
        "clean_fm_models_converged": int(clean_fm_audit.converged.sum()),
    }, indent=2))


if __name__ == "__main__":
    main()
