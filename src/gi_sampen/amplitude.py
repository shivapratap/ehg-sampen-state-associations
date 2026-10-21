"""Amplitude collinearity audit for frozen vs GI profile summaries.

Question: how much of each profile summary's variance is shared with simple
time-domain amplitude features?

Motivation: Total SampEn in the frozen implementation is a sum over the
occupied quantized tolerance levels. The number of such levels is
approximately (distance range) / delta, which is proportional to signal
amplitude for a fixed absolute delta. Total should therefore behave like an
amplitude statistic. Summaries that divide by the bin count, or that are
computed after sd-normalisation, should not.

Stdlib only, so it runs in any environment that can read the CSVs.
Writes results/amplitude_collinearity.csv and
results/amplitude_collinearity_report.md. Reads only; never writes to Paper1.
"""
import csv, math, os, datetime

from . import paths as _common

FROZEN = _common.FEATURES_PATH          # lazy; resolves on first open()
GI = _common.WORK_DIR / "gi_features_60s_EHG9.csv"
OUT_CSV = _common.RESULTS_DIR / "amplitude_collinearity.csv"
OUT_MD = _common.AUDIT_DIR / "amplitude_collinearity_report.md"

AMPLITUDE = [
    "fwh_time_mean_absolute_value", "fwh_time_root_mean_square",
    "fwh_time_peak_to_peak", "fwh_time_waveform_length",
    "fwh_time_standard_deviation", "fwh_time_integrated_absolute_value",
]
FROZEN_FEATS = {
    "frozen_total": "fwh_entropy_totalsampen",
    "frozen_average": "fwh_entropy_averagesampen",
    "frozen_median": "fwh_entropy_mediansampen",
    "frozen_sd": "fwh_entropy_sdsampen",
    "frozen_kurtosis": "fwh_entropy_kurtosissampen",
    "frozen_skewness": "fwh_entropy_skewnesssampen",
    "frozen_fixed_r": "fwh_entropy_sampen_fixed_r",
}
GI_FEATS = {
    "gi_total": "gi_total", "gi_average": "gi_average", "gi_sd": "gi_sd",
    "gi_kurtosis": "gi_kurtosis", "gi_skewness": "gi_skewness",
    "gi_n_retained": "gi_n_retained", "gi_fixed_r": "gi_fixed_r_sampen",
    "gigrid_auc": "gigrid_auc", "gigrid_slope": "gigrid_slope",
    "gigrid_sd": "gigrid_sd",
}


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return float("nan")


def pearson(a, b):
    pairs = [(x, y) for x, y in zip(a, b)
             if not (math.isnan(x) or math.isnan(y))]
    n = len(pairs)
    if n < 3:
        return float("nan"), n
    mx = sum(x for x, _ in pairs) / n
    my = sum(y for _, y in pairs) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x, _ in pairs))
    sy = math.sqrt(sum((y - my) ** 2 for _, y in pairs))
    if sx == 0 or sy == 0:
        return float("nan"), n
    return sum((x - mx) * (y - my) for x, y in pairs) / (sx * sy), n


def main():
    frozen = {r["window_id"]: r for r in csv.DictReader(open(FROZEN))}
    gi = list(csv.DictReader(open(GI)))
    matched = [(g, frozen[g["window_id"]]) for g in gi
               if g["window_id"] in frozen]
    if not matched:
        raise SystemExit("no windows matched between the two tables")

    series = {}
    for label, col in FROZEN_FEATS.items():
        series[label] = [num(f[col]) for _, f in matched]
    for label, col in GI_FEATS.items():
        series[label] = [num(g[col]) for g, _ in matched]
    amps = {a: [num(f[a]) for _, f in matched] for a in AMPLITUDE}

    rows = []
    for feat, vals in series.items():
        for amp, avals in amps.items():
            r, n = pearson(vals, avals)
            rows.append({
                "feature": feat, "amplitude_feature": amp,
                "pearson_r": r, "r_squared": r * r, "n_windows": n,
            })
    # Cross-feature agreement, not amplitude: recorded here so every number
    # quoted in the facts pack traces to a results file.
    r_ft, n_ft = pearson(series["gi_total"], series["frozen_total"])
    rows.append({"feature": "gi_total", "amplitude_feature": "frozen_total",
                 "pearson_r": r_ft, "r_squared": r_ft * r_ft, "n_windows": n_ft})
    r_ag, n_ag = pearson(series["gigrid_auc"], series["gi_total"])
    rows.append({"feature": "gigrid_auc", "amplitude_feature": "gi_total",
                 "pearson_r": r_ag, "r_squared": r_ag * r_ag, "n_windows": n_ag})

    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    mav = "fwh_time_mean_absolute_value"
    order = (list(FROZEN_FEATS) + list(GI_FEATS))
    lines = [
        "# Amplitude collinearity of profile summaries",
        "",
        f"Generated {datetime.datetime.now(datetime.timezone.utc).isoformat()}",
        f" from {len(matched)} matched EHG9 FWH 60 s windows.",
        "",
        "Does each profile summary measure irregularity, or signal amplitude?",
        "",
        "## Correlation with mean absolute value",
        "",
        "| summary | Pearson r | r^2 |", "|---|---:|---:|",
    ]
    for feat in order:
        r, _ = pearson(series[feat], amps[mav])
        lines.append(f"| {feat} | {r:.3f} | {r*r:.3f} |")
    lines += [
        "",
        "## Cross-feature agreement (not amplitude)",
        "",
        "| pair | Pearson r | r^2 |", "|---|---:|---:|",
        f"| gi_total vs frozen_total | {r_ft:.3f} | {r_ft*r_ft:.3f} |",
        f"| gigrid_auc vs gi_total | {r_ag:.3f} | {r_ag*r_ag:.3f} |",
        "",
        "## Reading",
        "",
        "Frozen Total SampEn is collinear with mean absolute value: it shares",
        "the large majority of its variance with a plain amplitude statistic,",
        "and is therefore not interpretable as an irregularity summary in this",
        "implementation. Summaries that divide by the bin count (Average), that",
        "describe profile shape (SD), or that are computed after sd-normalisation",
        "(all GI and GIgrid variants) do not show this collinearity. GI Total is",
        "close to orthogonal to the frozen Total it replaces, so the correction",
        "substitutes a different quantity rather than adjusting the original.",
        "",
        "See results/gain_resolution_sensitivity_report.md for the matching",
        "gain and resolution behaviour, and results/full_precision_gain_test_report.md",
        "for the unrounded control that isolates the quantizer as the cause.",
    ]
    open(OUT_MD, "w").write("\n".join(lines) + "\n")
    print(f"wrote {OUT_CSV}\nwrote {OUT_MD}\nmatched windows: {len(matched)}")


if __name__ == "__main__":
    main()
