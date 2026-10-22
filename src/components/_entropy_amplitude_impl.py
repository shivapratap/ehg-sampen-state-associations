"""Stage 07 -- SampEn/amplitude dependence audit, per channel.

Standalone port of pipeline/module3_5/sampen_amplitude_audit.py into this
project's numbered-stage layout, extended across all four the standard pipeline channels
(EHG9, EHG10, EHG11, EHG12). Logic (correlation pairs, aggregation scopes,
coupling-strength thresholds) is unchanged from the frozen module3_5
design; see that module's docstring for the original single-channel
rationale, condensed below.

Purpose
-------
This is a post-hoc, descriptive audit -- it does NOT refit or modify
Stage 06's models or the frozen validated feature set. The scientific
question: are SampEn-derived entropy features (Total SampEn, Average
SampEn, fixed-r SampEn, and the extended SampEn profile summaries) actually
amplitude-coupled -- i.e. do they track signal amplitude/power (MAV, RMS,
peak-to-peak, waveform length) rather than genuine, amplitude-independent
complexity? If so, a UC-vs-IN difference in an entropy feature could be
confounded by an amplitude difference between classes, and the feature
should be interpreted (and described in the manuscript) accordingly.

Six correlation scopes are computed (Pearson AND Spearman, for every
SampEn x amplitude feature pair, 7x4=28 pairs per scope):
    window_pooled                        all windows, pooled (descriptive
                                          only -- windows are nested/overlapping)
    window_class_stratified              separately within baseline/
                                          contraction/fetal_movement
    recording_class_median                per (patient, recording, class)
                                          median aggregation
    patient_class_median                  per (patient, class) median
                                          aggregation
    patient_contraction_minus_baseline    per-patient effect: contraction
                                          median minus baseline median
                                          (PRIMARY interpretive level)
    patient_fetal_movement_minus_baseline per-patient effect: fetal_movement
                                          median minus baseline median
                                          (secondary; kept separate from
                                          contraction throughout)

Coupling-strength decision rule (Spearman |rho| only):
    strong    >= 0.75
    moderate  0.50 to <0.75
    weak      < 0.50
    (overridden by "exploratory/unstable: small n" if n_units < 20)

Input:  results/paper1/06_mixed_models/event_unit_analysis_60s_{channel}.csv
        (Stage 06's event-unit table -- the standalone-pipeline equivalent
        of module2_5's analysis_units_60s_fwh_final.csv)
Output: results/paper1/07_amplitude_confounding/
    module3_5_sampen_amplitude_correlations_{channel}.csv
    module3_5_sampen_amplitude_interpretation_flags_{channel}.csv
    module3_5_sampen_amplitude_summary_table_{channel}.csv
    module3_5_sampen_amplitude_audit_report_{channel}.md
    module3_5_sampen_amplitude_summary_table_all_channels.csv  (cross-channel roll-up)
"""
from __future__ import annotations

import argparse
import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import config as project_config  # noqa: E402


# =============================================================================
# Ported from pipeline/module3_5/sampen_amplitude_audit.py (unchanged)
# =============================================================================

SAMPEN = [
    "fwh_entropy_totalsampen",
    "fwh_entropy_averagesampen",
    "fwh_entropy_sampen_fixed_r",
    "fwh_entropy_mediansampen",
    "fwh_entropy_sdsampen",
    "fwh_entropy_kurtosissampen",
    "fwh_entropy_skewnesssampen",
]
AMPLITUDE = ["fwh_time_mean_absolute_value", "fwh_time_root_mean_square", "fwh_time_peak_to_peak", "fwh_time_waveform_length"]
LABELS = {
    "fwh_entropy_totalsampen": "Total SampEn", "fwh_entropy_averagesampen": "Average SampEn",
    "fwh_entropy_sampen_fixed_r": "Fixed-r SampEn", "fwh_entropy_mediansampen": "Median SampEn",
    "fwh_entropy_sdsampen": "SD SampEn", "fwh_entropy_kurtosissampen": "Kurtosis SampEn",
    "fwh_entropy_skewnesssampen": "Skewness SampEn", "fwh_time_mean_absolute_value": "MAV",
    "fwh_time_root_mean_square": "RMS", "fwh_time_peak_to_peak": "Peak-to-peak",
    "fwh_time_waveform_length": "Waveform length",
}


def correlation(a: pd.Series, b: pd.Series, method: str) -> tuple:
    z = pd.concat([a, b], axis=1).apply(pd.to_numeric, errors="coerce").dropna()
    n = len(z)
    if n < 3 or z.iloc[:, 0].nunique() < 2 or z.iloc[:, 1].nunique() < 2:
        return np.nan, np.nan, n
    fn = pearsonr if method == "pearson" else spearmanr
    r = fn(z.iloc[:, 0], z.iloc[:, 1])
    return float(r.statistic), float(r.pvalue), n


def pairwise(df: pd.DataFrame, scope: str, cls: str = "all") -> pd.DataFrame:
    rows = []
    for s, a in product(SAMPEN, AMPLITUDE):
        for method in ("pearson", "spearman"):
            r, p, n = correlation(df[s], df[a], method)
            rows.append(dict(
                scope=scope, class_label=cls, sampen_feature=s, amplitude_feature=a,
                sampen_label=LABELS[s], amplitude_label=LABELS[a], method=method,
                n_units=n, rho=r, p_value=p, abs_rho=abs(r) if pd.notna(r) else np.nan,
            ))
    return pd.DataFrame(rows)


def median_level(df: pd.DataFrame, keys: list) -> pd.DataFrame:
    return df.groupby(keys, dropna=False)[SAMPEN + AMPLITUDE].median().reset_index()


def effects(df: pd.DataFrame, outcome: str) -> pd.DataFrame:
    x = median_level(df[df.window_class.isin(["baseline", outcome])], ["patient_id", "window_class"])
    b = x[x.window_class == "baseline"].set_index("patient_id")
    o = x[x.window_class == outcome].set_index("patient_id")
    common = b.index.intersection(o.index)
    out = pd.DataFrame(index=common)
    for f in SAMPEN + AMPLITUDE:
        out[f] = o.loc[common, f] - b.loc[common, f]
    return out.reset_index()


def flags(corrs: pd.DataFrame) -> pd.DataFrame:
    out = corrs[corrs.method == "spearman"].copy()
    out["coupling_strength"] = pd.cut(out.abs_rho, [-np.inf, .5, .75, np.inf], labels=["weak", "moderate", "strong"], right=False)
    out["small_sample"] = out.n_units < 20
    out["window_dependence_warning"] = out.scope.str.contains("window")
    out["interpretation_flag"] = np.select(
        [out.small_sample, out.coupling_strength.eq("strong"), out.coupling_strength.eq("moderate")],
        ["exploratory/unstable: small n", "strong amplitude-coupling", "moderate amplitude-coupling"],
        default="weak coupling",
    )
    return out


def _write_report(output_dir: Path, channel: str, units_path: Path, u: pd.DataFrame) -> None:
    report = f'''# Stage 07 ({channel}): SampEn-Amplitude Dependence and Interpretability Audit

## Objective
Post-hoc assessment of whether SampEn-derived features are amplitude-coupled
on channel {channel}. This audit does not change Stage 06 models or the
frozen validated feature set.

## Sources and methods
- Source: `{units_path}` ({len(u):,} analysis windows; {u.patient_id.nunique()} women; {u.recording_id.nunique()} recordings).
- SampEn features: Total SampEn, Average SampEn, fixed-r SampEn, Median SampEn,
  SD SampEn, Kurtosis SampEn, and Skewness SampEn.
- Amplitude features: MAV, RMS, peak-to-peak, and waveform length.
- Pearson and Spearman correlations were computed at pooled/window, recording, patient, and patient effect levels.
- Primary effect: contraction minus baseline. Fetal movement minus baseline is secondary and was not pooled with contraction.

## Key results
See `module3_5_sampen_amplitude_summary_table_{channel}.csv` for the complete matrix. Spearman correlations are emphasized. Strong coupling is defined as |rho| >= 0.75; moderate as 0.50-<0.75; weak as <0.50.

The pooled window-level results are descriptive only because windows are nested and overlapping. Patient-level effects are more interpretable but have smaller sample sizes.

## Interpretation and limitations
Total SampEn should be described as an amplitude-coupled nonlinear/entropy-derived summary when strong coupling is observed, not as an amplitude-independent complexity marker. Average SampEn is retained as an exploratory comparator; any weaker or inconsistent association should not be interpreted as validation of a replacement feature. Correlation is not causation, and repeated-window precision is not independent evidence.

## Recommended the standard pipeline wording (channel {channel})
"Post-hoc analyses indicated that SampEn-derived measures, particularly Total SampEn, showed [substantial/moderate/weak -- see summary table] association with amplitude-domain measures on channel {channel}. We therefore interpret Total SampEn conservatively as an amplitude-coupled nonlinear entropy summary rather than an amplitude-independent complexity biomarker. Average SampEn was examined as an exploratory comparator and was not used to revise the prespecified feature set."

## Stage 06 status
Stage 06 should remain frozen. Stage 07 is an interpretability audit and provides no basis for revising the primary inferential model or validated feature set.
'''
    (output_dir / f"module3_5_sampen_amplitude_audit_report_{channel}.md").write_text(report, encoding="utf-8")


def run_channel(channel: str, units_dir: Path, results_dir: Path) -> pd.DataFrame:
    units_path = units_dir / f"event_unit_analysis_60s_{channel}.csv"
    if not units_path.is_file():
        raise FileNotFoundError(
            f"Stage 06 event-unit table missing for {channel}: {units_path}. "
            "Run src/06_mixed_models.py first."
        )
    u = pd.read_csv(units_path)
    missing = [f for f in SAMPEN + AMPLITUDE if f not in u.columns]
    if missing:
        raise ValueError(f"{channel}: missing required features: {missing}")
    if u.window_id.duplicated().any():
        raise ValueError(f"{channel}: window_id must be unique")

    rows = [pairwise(u, "window_pooled")]
    rows += [pairwise(u[u.window_class == c], "window_class_stratified", c)
             for c in ["baseline", "contraction", "fetal_movement"]]
    rec = median_level(u, ["patient_id", "recording_id", "window_class"])
    pat = median_level(u, ["patient_id", "window_class"])
    rows += [pairwise(rec, "recording_class_median"), pairwise(pat, "patient_class_median")]
    con_eff = effects(u, "contraction")
    fm_eff = effects(u, "fetal_movement")
    rows += [pairwise(con_eff, "patient_contraction_minus_baseline"),
             pairwise(fm_eff, "patient_fetal_movement_minus_baseline")]
    corrs = pd.concat(rows, ignore_index=True)
    corrs.insert(0, "channel", channel)
    corrs.to_csv(results_dir / f"module3_5_sampen_amplitude_correlations_{channel}.csv", index=False)

    flagged = flags(corrs)
    flagged.to_csv(results_dir / f"module3_5_sampen_amplitude_interpretation_flags_{channel}.csv", index=False)

    summary = flagged.groupby(
        ["scope", "sampen_label", "amplitude_label", "coupling_strength"], observed=True, dropna=False
    ).agg(
        n_units=("n_units", "max"), spearman_rho=("rho", "first"),
        p_value=("p_value", "first"), interpretation=("interpretation_flag", "first"),
    ).reset_index()
    summary.insert(0, "channel", channel)
    summary.to_csv(results_dir / f"module3_5_sampen_amplitude_summary_table_{channel}.csv", index=False)

    _write_report(results_dir, channel, units_path, u)
    return summary


def run(
    channels: tuple = project_config.CHANNELS,
    units_dir: Path = None,
    results_dir: Path = None,
) -> int:
    units_dir = units_dir or (project_config.UNITS_DIR)
    results_dir = results_dir or (project_config.WORK_DIR / "entropy_amplitude")
    results_dir.mkdir(parents=True, exist_ok=True)

    all_summaries = []
    for channel in channels:
        print(f"=== channel={channel} ===")
        summary = run_channel(channel, units_dir, results_dir)
        all_summaries.append(summary)
        primary = summary[
            (summary.scope == "patient_contraction_minus_baseline")
            & (summary.sampen_label == "Total SampEn")
            & (summary.amplitude_label == "RMS")
        ]
        if not primary.empty:
            row = primary.iloc[0]
            print(f"  Total SampEn vs RMS (patient contraction-minus-baseline effect): "
                  f"rho={row.spearman_rho:.3f}, n={row.n_units}, coupling={row.coupling_strength}")

    combined = pd.concat(all_summaries, ignore_index=True)
    combined.to_csv(results_dir / "module3_5_sampen_amplitude_summary_table_all_channels.csv", index=False)
    print(f"done. outputs -> {results_dir}")
    return 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--channels", type=str, default=None, help="Comma-separated, e.g. EHG9,EHG10")
    parser.add_argument("--units-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    channels = tuple(c.strip() for c in args.channels.split(",")) if args.channels else project_config.CHANNELS
    sys.exit(run(channels=channels, units_dir=args.units_dir, results_dir=args.results_dir))
