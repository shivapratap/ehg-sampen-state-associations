#!/usr/bin/env python3
"""Figure 1 -- mechanism and consequence.

(a) Amplitude collinearity: standard Total SampEn and GI-Total SampEn against
    time-domain mean absolute value (MAV), across the matched EHG9 windows.
(b) Cross-channel consequence: FM-versus-IN standardised effect estimates for
    SD and Total SampEn, standard versus GI, EHG9-EHG12.

This is a VISUALISATION script. It runs no models, recomputes no FDR, and
writes nothing outside results/gi/figures/. Panel A is a deterministic join of
two existing feature tables from the standard and GI implementations on window_id; panel B is read verbatim from
the cross-channel result table.

Two layouts, both reproducible from existing outputs:

    collinearity  (a) standard/GI Total against amplitude, (b) cross-channel
                  consequence.
    mechanism     (a) tolerance-support scaling against gain -- the cause;
                  (b) gain sensitivity of profile summaries, including the
                  AvgSampEn and unrounded controls -- the cheap fixes; (c)
                  cross-channel FM-versus-IN consequence.

Run from anywhere:

    python src/gi_sampen/figure1.py
    python src/gi_sampen/figure1.py --layout mechanism
    python src/gi_sampen/figure1.py --display zlinear   # the literal
                                                          # z-scored spec
    python src/gi_sampen/figure1.py --ci                # add 95% CIs to (b)

Display modes for panel (a)
---------------------------
loglog   (default) raw units on log10 axes; shared x range, independent y
         range (the two Totals are not the same quantity on a common scale).
zlog     z-scores of log10(value); identical x and y limits in both panels.
zlinear  z-scores of the raw values; identical limits. This is the literal
         "standardised axes" specification. It is NOT recommended here: MAV is
         a positive amplitude spanning three decades, so on linear z-scores
         about 99.9% of windows fall inside z < 1 while the largest sits near
         z = 63, and the cloud collapses into one corner. Provided so the
         effect can be seen rather than taken on trust.

In every mode the annotated Pearson r is computed on the ORIGINAL,
untransformed feature values, exactly as reported in the facts pack. Spearman
rho is annotated alongside it because it is invariant to the display
transform, and because it is the honest answer to "is that r driven by the
amplitude outliers?" (it is not).
"""
from __future__ import annotations

import argparse
import datetime as _dt
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager, ticker
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# --------------------------------------------------------------------------
# Paths -- mirrors src/gi_sampen/paths.py so this script stays in step with the
# rest of the pipeline, without importing it (that module pulls in the standard pipeline's
# preprocessing and requires EHG_RAW_DIR, which a figure does not need).
# --------------------------------------------------------------------------
from . import paths as _common

REPO_ROOT = _common.REPO_ROOT

GI_FEATURES_PATH = _common.WORK_DIR / "gi_features_60s_EHG9.csv"
CROSS_CHANNEL_PATH = _common.RESULTS_DIR / "cross_channel_replication_summary.csv"
FIGURE_DIR = _common.FIGURE_DIR


def _frozen_feature_candidates():
    """Cached standard finite-resolution feature table, canonical location
    first and the published-results copy second."""
    root = _common.upstream_dir()
    return (
        root / "work" / "features" / "window_features_60s_EHG9.csv",
        root / "results" / "paper1" / "04_feature_extraction"
        / "window_features_60s_EHG9.csv",
    )


# Column names, read from the repository rather than guessed.
COL_WINDOW_ID = "window_id"
COL_MAV = "fwh_time_mean_absolute_value"
COL_FROZEN_TOTAL = "fwh_entropy_totalsampen"
COL_GI_TOTAL = "gi_total"

CONTRAST = "fetal_movement_vs_baseline"
PANEL_B_FEATURES = ("gi_sd", "gi_total")
FEATURE_TITLES = {"gi_sd": "SD SampEn", "gi_total": "Total SampEn"}
CHANNEL_ORDER = ("EHG9", "EHG10", "EHG11", "EHG12")

INVARIANCE_WINDOWS_PATH = _common.WORK_DIR / "invariance_window_results.csv"
FULL_PRECISION_PATH = _common.RESULTS_DIR / "full_precision_gain_test.csv"
AVG_GAIN_SUMMARY_PATH = _common.RESULTS_DIR / "avg_sampen_gain_summary.csv"
COL_SCALE = "scale_factor"
BASELINE_GAIN = 1.0

# Validated values (facts pack SS5.8). Reproduced, never imposed.
EXPECTED = {
    "n_windows": 4083,
    "frozen_total_vs_mav_r": 0.984,
    "gi_total_vs_mav_r": -0.102,
}
R_TOLERANCE = 0.002

# --------------------------------------------------------------------------
# Visual system
# --------------------------------------------------------------------------
# The standard implementation is deliberately achromatic: it is the uncorrected baseline, GI is the
# result under discussion. #949494 rather than #7A7A7A -- the latter sits at
# OKLab dE 14.0 from #0072B2 for normal colour vision, below the dE 15 floor at
# which two adjacent categorical marks stop being reliably separable. #949494
# gives dE 18.8 (normal) / 14.8 (protanopia) against the same blue. Standard bars
# additionally carry a hatch, so identity survives greyscale printing.
COLOR_GI = "#0072B2"
# Third series: the unrounded (full-precision) control. Vermillion rather than
# green -- green/blue collapses to dE 8.6 under tritanopia, vermillion/blue
# holds at 14.8. Line style and marker carry identity as well as hue.
COLOR_UNROUNDED = "#D55E00"
COLOR_FROZEN = "#949494"
COLOR_FROZEN_LINE = "#6E6E6E"  # darker, for the trend line over grey points
INK = "#1A1A1A"
INK_MUTED = "#555555"
ZERO_LINE = "#333333"

FS_PANEL = 9.5
FS_TITLE = 8.5
FS_AXIS = 8.5
FS_TICK = 7.5
FS_LEGEND = 7.5
FS_ANNOT = 7.5

FIG_W, FIG_H = 7.1, 4.7


def _pick_font() -> str:
    """First available of the preferred families. Installs nothing."""
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in ("Arial", "Helvetica", "Helvetica Neue", "Liberation Sans",
                 "Nimbus Sans", "DejaVu Sans"):
        if name in available:
            return name
    return "DejaVu Sans"


def _apply_style(font: str) -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [font, "DejaVu Sans"],
        "font.size": FS_TICK,
        "axes.labelsize": FS_AXIS,
        "axes.titlesize": FS_TITLE,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND,
        "axes.edgecolor": "#444444",
        "axes.linewidth": 0.7,
        "axes.labelcolor": INK,
        "text.color": INK,
        "xtick.color": "#444444",
        "ytick.color": "#444444",
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 2.6,
        "ytick.major.size": 2.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "grid.color": "#D8D8D8",
        "grid.linewidth": 0.5,
        "figure.dpi": 150,
        "savefig.bbox": None,
        "pdf.fonttype": 42,   # TrueType, editable text in the vector PDF
        "ps.fonttype": 42,
    })


# --------------------------------------------------------------------------
# Data
# --------------------------------------------------------------------------
def _resolve_frozen_features() -> Path:
    for candidate in _frozen_feature_candidates():
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "No standard-implementation EHG9 feature table found. Looked in:\n  "
        + "\n  ".join(str(p) for p in _frozen_feature_candidates()))


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.corrcoef(a, b)[0, 1])


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra = pd.Series(a).rank().to_numpy()
    rb = pd.Series(b).rank().to_numpy()
    return _pearson(ra, rb)


def load_panel_a() -> tuple[pd.DataFrame, dict]:
    """Deterministic join of the standard and GI feature tables on window_id."""
    frozen_path = _resolve_frozen_features()
    frozen = pd.read_csv(frozen_path, low_memory=False)
    gi = pd.read_csv(GI_FEATURES_PATH, low_memory=False)

    for col, frame, name in ((COL_MAV, frozen, frozen_path.name),
                             (COL_FROZEN_TOTAL, frozen, frozen_path.name),
                             (COL_GI_TOTAL, gi, GI_FEATURES_PATH.name)):
        if col not in frame.columns:
            raise KeyError(f"column {col!r} not present in {name}")

    checks = {}
    checks["frozen_rows"] = len(frozen)
    checks["gi_rows"] = len(gi)
    checks["frozen_duplicate_ids"] = int(frozen[COL_WINDOW_ID].duplicated().sum())
    checks["gi_duplicate_ids"] = int(gi[COL_WINDOW_ID].duplicated().sum())
    if checks["frozen_duplicate_ids"] or checks["gi_duplicate_ids"]:
        raise ValueError(
            "duplicate window_id values present; the join would not be "
            f"one-to-one (frozen={checks['frozen_duplicate_ids']}, "
            f"gi={checks['gi_duplicate_ids']})")

    merged = gi[[COL_WINDOW_ID, COL_GI_TOTAL]].merge(
        frozen[[COL_WINDOW_ID, COL_MAV, COL_FROZEN_TOTAL]],
        on=COL_WINDOW_ID, how="inner", validate="one_to_one")
    merged = merged.dropna(subset=[COL_MAV, COL_FROZEN_TOTAL, COL_GI_TOTAL])

    checks["matched_windows"] = len(merged)
    checks["join_one_to_one"] = True
    if len(merged) != EXPECTED["n_windows"]:
        raise ValueError(
            f"expected {EXPECTED['n_windows']} matched windows, got {len(merged)}. "
            "Stopping rather than drawing a figure from a different sample.")

    mav = merged[COL_MAV].to_numpy(float)
    checks["frozen_r"] = _pearson(merged[COL_FROZEN_TOTAL].to_numpy(float), mav)
    checks["gi_r"] = _pearson(merged[COL_GI_TOTAL].to_numpy(float), mav)
    checks["frozen_rho"] = _spearman(merged[COL_FROZEN_TOTAL].to_numpy(float), mav)
    checks["gi_rho"] = _spearman(merged[COL_GI_TOTAL].to_numpy(float), mav)
    checks["frozen_features_path"] = frozen_path

    for key, expected in (("frozen_r", EXPECTED["frozen_total_vs_mav_r"]),
                          ("gi_r", EXPECTED["gi_total_vs_mav_r"])):
        if abs(checks[key] - expected) > R_TOLERANCE:
            raise ValueError(
                f"{key} = {checks[key]:+.4f} does not reproduce the validated "
                f"{expected:+.3f} (tolerance {R_TOLERANCE}). Stopping; do not "
                "draw the figure until this is explained.")
    return merged, checks


def load_support_vs_gain(common_subset: bool = False) -> tuple[pd.DataFrame, dict]:
    """Median retained tolerance-support size at each gain, per construction.

    Standard and GI come from the 60-window gain-stress table; the unrounded
    control from the 12-window full-precision table (a strict subset of the
    60). With common_subset=True, standard and GI are restricted to those same
    12 windows so all three series rest on identical data. Median over
    windows only -- no model, no estimation.
    """
    inv = pd.read_csv(INVARIANCE_WINDOWS_PATH, low_memory=False)
    for col in (COL_SCALE, "frozen_n_retained", "gi_n_retained"):
        if col not in inv.columns:
            raise KeyError(f"{col!r} missing from {INVARIANCE_WINDOWS_PATH.name}")

    fp = pd.read_csv(FULL_PRECISION_PATH, low_memory=False)
    if "fp_n_retained" not in fp.columns:
        raise KeyError(f"'fp_n_retained' missing from {FULL_PRECISION_PATH.name}")

    fp_ids = set(fp[COL_WINDOW_ID])
    inv_ids = set(inv[COL_WINDOW_ID])
    counts = {"n_frozen_gi_available": len(inv_ids),
              "n_unrounded": len(fp_ids),
              "unrounded_is_subset": fp_ids <= inv_ids,
              "common_subset": bool(common_subset)}
    if not counts["unrounded_is_subset"]:
        raise ValueError(
            "the full-precision windows are not a subset of the gain-stress "
            "windows; the series would not be comparable at all")

    if common_subset:
        inv = inv[inv[COL_WINDOW_ID].isin(fp_ids)]
    counts["n_frozen_gi_used"] = inv[COL_WINDOW_ID].nunique()

    med = (inv.groupby(COL_SCALE)[["frozen_n_retained", "gi_n_retained"]]
              .median().reset_index().rename(columns={COL_SCALE: "gain"}))
    fp_med = (fp.groupby("gain")["fp_n_retained"].median()
                .reset_index().rename(columns={"fp_n_retained": "unrounded_n_retained"}))

    out = med.merge(fp_med, on="gain", how="outer").sort_values("gain")
    if BASELINE_GAIN not in set(out["gain"]):
        raise ValueError("no unit-gain baseline row in the gain-stress tables")
    return out.reset_index(drop=True), counts


def _median_abs_deviation(df, col, gain_col):
    """Median |relative deviation| from each window's own unit-gain value."""
    base = df[df[gain_col] == BASELINE_GAIN].set_index(COL_WINDOW_ID)[col]
    rows = []
    for gain, block in df[df[gain_col] != BASELINE_GAIN].groupby(gain_col):
        ref = base.reindex(block[COL_WINDOW_ID]).to_numpy(float)
        got = block[col].to_numpy(float)
        rows.append({"gain": float(gain),
                     "abs_dev": float(np.nanmedian(np.abs(got - ref) / np.abs(ref)) * 100.0)})
    return pd.DataFrame(rows).sort_values("gain").reset_index(drop=True)


def load_deviation_vs_gain(common_subset: bool = False) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Median |relative deviation| from the unit-gain baseline, by gain.

    By default the standard and GI series are read verbatim from the standard
    AvgSampEn gain summary (60 windows). With common_subset=True they are
    reduced from the per-window gain-stress table restricted to the 12
    full-precision windows; that reduction is first validated against the
    published summary on the full 60, so the restricted curve is trustworthy.
    The unrounded series always comes from the full-precision table and is
    checked against its published verdict of exactly zero.
    """
    summ = pd.read_csv(AVG_GAIN_SUMMARY_PATH)
    need = {"feature", "mode", "gain", "median_relative_difference_percent"}
    if not need <= set(summ.columns):
        raise KeyError(f"{AVG_GAIN_SUMMARY_PATH.name} is missing {need - set(summ.columns)}")
    summ = summ.copy()
    summ["abs_dev"] = summ["median_relative_difference_percent"].abs()

    fp = pd.read_csv(FULL_PRECISION_PATH, low_memory=False)
    fp_dev = _median_abs_deviation(fp, "fp_total", "gain")
    worst = float(fp_dev["abs_dev"].max())
    if worst > 1e-6:
        raise ValueError(
            f"unrounded control deviates by {worst:.3g}% -- the published result "
            "is exactly zero. Stopping rather than drawing a figure that "
            "disagrees with the standard-implementation report.")

    notes = {"frozen_gi_source": _rel(AVG_GAIN_SUMMARY_PATH),
             "n_unrounded": fp[COL_WINDOW_ID].nunique(),
             "n_frozen_gi": int(summ["n_windows"].iloc[0]) if "n_windows" in summ else 60}

    if common_subset:
        inv = pd.read_csv(INVARIANCE_WINDOWS_PATH, low_memory=False)
        # Validate the reduction against the published summary on the full set
        # before trusting it on the restricted one.
        check = _median_abs_deviation(inv, "frozen_total", COL_SCALE).set_index("gain")["abs_dev"]
        pub = (summ[(summ["feature"] == "total") & (summ["mode"] == "frozen")]
               .set_index("gain")["abs_dev"])
        shared = check.index.intersection(pub.index)
        drift = float((check[shared] - pub[shared]).abs().max())
        if drift > 0.05:
            raise ValueError(
                f"recomputed standard-implementation deviation differs from the published summary "
                f"by up to {drift:.3g} percentage points; not substituting it")
        notes["reduction_validated_drift_pp"] = drift

        fp_ids = set(fp[COL_WINDOW_ID])
        inv = inv[inv[COL_WINDOW_ID].isin(fp_ids)]
        notes["n_frozen_gi"] = inv[COL_WINDOW_ID].nunique()
        notes["frozen_gi_source"] = _rel(INVARIANCE_WINDOWS_PATH) + " (restricted)"
        rebuilt = []
        for feature, col, mode in (("total", "frozen_total", "frozen"),
                                   ("average", "frozen_average", "frozen"),
                                   ("total", "gi_total", "GI")):
            d = _median_abs_deviation(inv, col, COL_SCALE)
            d["feature"], d["mode"] = feature, mode
            rebuilt.append(d)
        summ = pd.concat(rebuilt, ignore_index=True)
        summ["median_relative_difference_percent"] = np.nan

    return summ, fp_dev, notes


def load_panel_b() -> pd.DataFrame:
    frame = pd.read_csv(CROSS_CHANNEL_PATH)
    sel = frame[(frame["contrast"] == CONTRAST)
                & (frame["feature"].isin(PANEL_B_FEATURES))].copy()
    if sel.empty:
        raise ValueError(f"no rows for contrast {CONTRAST!r} in {CROSS_CHANNEL_PATH}")
    missing = set(CHANNEL_ORDER) - set(sel["channel"])
    if missing:
        raise ValueError(f"channels missing from the standard result table: {sorted(missing)}")
    sel["channel"] = pd.Categorical(sel["channel"], CHANNEL_ORDER, ordered=True)
    keep = ["feature", "channel", "estimate", "q_value",
            "frozen_estimate", "frozen_q_value", "ci_low", "ci_high", "n", "n_women"]
    keep = [c for c in keep if c in sel.columns]
    return sel[keep].sort_values(["feature", "channel"]).reset_index(drop=True)


# --------------------------------------------------------------------------
# Panel A
# --------------------------------------------------------------------------
def _transform(values: np.ndarray, mode: str) -> np.ndarray:
    if mode == "loglog":
        return values
    if mode == "zlog":
        v = np.log10(values)
    elif mode == "zlinear":
        v = values.astype(float)
    else:
        raise ValueError(f"unknown display mode {mode!r}")
    return (v - v.mean()) / v.std(ddof=0)


def _thin_log_axis(axis, values) -> None:
    """Label intermediate decades when the data spans a short log range."""
    decades = np.log10(np.nanmax(values)) - np.log10(np.nanmin(values))
    if decades < 1.5:
        axis.set_major_locator(ticker.LogLocator(base=10, subs=(1, 2, 3, 5), numticks=12))
        axis.set_major_formatter(ticker.FuncFormatter(
            lambda v, _pos: f"{v:g}" if v >= 1 else f"{v:g}"))
        axis.set_minor_formatter(ticker.NullFormatter())


def _draw_scatter(ax, x, y, *, color, line_color, title, r, rho, mode,
                  xlabel, ylabel, show_ylabel):
    ax.scatter(x, y, s=3.2, c=color, alpha=0.16, linewidths=0,
               rasterized=True, zorder=2)

    # Least-squares guide, fitted in the plotted coordinates.
    if mode == "loglog":
        lx, ly = np.log10(x), np.log10(y)
        slope, intercept = np.polyfit(lx, ly, 1)
        lo, hi = np.percentile(lx, [1, 99])
        gx = np.linspace(lo, hi, 100)
        ax.plot(10 ** gx, 10 ** (slope * gx + intercept),
                color=line_color, lw=1.1, zorder=3, solid_capstyle="round")
        ax.set_xscale("log")
        ax.set_yscale("log")
        _thin_log_axis(ax.yaxis, y)
    else:
        slope, intercept = np.polyfit(x, y, 1)
        lo, hi = np.percentile(x, [1, 99])
        gx = np.linspace(lo, hi, 100)
        ax.plot(gx, slope * gx + intercept,
                color=line_color, lw=1.1, zorder=3, solid_capstyle="round")

    ax.set_title(title, fontsize=FS_TITLE, color=INK, pad=4)
    ax.set_xlabel(xlabel, labelpad=2)
    if show_ylabel:
        ax.set_ylabel(ylabel, labelpad=2)
    ax.annotate(f"$r$ = {r:+.3f}   $r^2$ = {r * r:.3f}\n"
                fr"$\rho$ = {rho:+.3f}",
                xy=(0.965, 0.035), xycoords="axes fraction",
                va="bottom", ha="right", fontsize=FS_ANNOT, color=INK_MUTED,
                linespacing=1.35, zorder=6,
                bbox=dict(boxstyle="square,pad=0.22", facecolor="white",
                          edgecolor="none", alpha=0.85))
    ax.tick_params(length=2.6, pad=1.8)


def draw_panel_a(ax_left, ax_right, data: pd.DataFrame, checks: dict, mode: str):
    mav_raw = data[COL_MAV].to_numpy(float)
    fro_raw = data[COL_FROZEN_TOTAL].to_numpy(float)
    gi_raw = data[COL_GI_TOTAL].to_numpy(float)

    x = _transform(mav_raw, mode)
    y_fro = _transform(fro_raw, mode)
    y_gi = _transform(gi_raw, mode)

    if mode == "loglog":
        xlabel, ylabel = "MAV (mV, log scale)", "Total SampEn (log scale)"
    elif mode == "zlog":
        xlabel, ylabel = "MAV (standardised $\\log_{10}$)", "Total SampEn (standardised $\\log_{10}$)"
    else:
        xlabel, ylabel = "MAV (standardised)", "Total SampEn (standardised)"

    _draw_scatter(ax_left, x, y_fro, color=COLOR_FROZEN,
                  line_color=COLOR_FROZEN_LINE, title="Standard Total SampEn",
                  r=checks["frozen_r"], rho=checks["frozen_rho"], mode=mode,
                  xlabel=xlabel, ylabel=ylabel, show_ylabel=True)
    _draw_scatter(ax_right, x, y_gi, color=COLOR_GI, line_color=COLOR_GI,
                  title="GI-Total SampEn",
                  r=checks["gi_r"], rho=checks["gi_rho"], mode=mode,
                  xlabel=xlabel, ylabel=ylabel, show_ylabel=False)

    # Identical x limits always; identical y limits except in loglog, where the
    # two quantities have genuinely different ranges and a shared axis would
    # read as a variance claim rather than an association claim.
    xlo = min(ax_left.get_xlim()[0], ax_right.get_xlim()[0])
    xhi = max(ax_left.get_xlim()[1], ax_right.get_xlim()[1])
    for ax in (ax_left, ax_right):
        ax.set_xlim(xlo, xhi)
    if mode != "loglog":
        ylo = min(ax_left.get_ylim()[0], ax_right.get_ylim()[0])
        yhi = max(ax_left.get_ylim()[1], ax_right.get_ylim()[1])
        for ax in (ax_left, ax_right):
            ax.set_ylim(ylo, yhi)
    ax_right.tick_params(labelleft=True)


# --------------------------------------------------------------------------
# Mechanism panels
# --------------------------------------------------------------------------
SERIES = (
    # key, label, colour, linestyle, marker
    ("frozen", "Standard", COLOR_FROZEN_LINE, "-", "o"),
    ("gi", "GI", COLOR_GI, "-", "s"),
    ("unrounded", "Unrounded profile", COLOR_UNROUNDED, "--", "^"),
)


def draw_support_panel(ax, support: pd.DataFrame, counts: dict):
    """(a) Retained tolerance-support size against signal gain."""
    gains = support["gain"].to_numpy(float)
    cols = {"frozen": "frozen_n_retained", "gi": "gi_n_retained",
            "unrounded": "unrounded_n_retained"}

    # Measured log-log slope of the standard support against gain. A slope of 1
    # is exact proportionality; annotating the fitted value is more use than a
    # reference line, which would sit invisibly underneath the standard series.
    base = float(support.loc[support["gain"] == BASELINE_GAIN,
                             "frozen_n_retained"].iloc[0])
    fz = support["frozen_n_retained"].to_numpy(float)
    ok = np.isfinite(fz) & (fz > 0)
    slope = float(np.polyfit(np.log10(gains[ok]), np.log10(fz[ok]), 1)[0])
    ax.annotate(f"Standard slope = {slope:.2f}", xy=(0.97, 0.035),
                xycoords="axes fraction", ha="right", va="bottom",
                fontsize=FS_ANNOT, color=COLOR_FROZEN_LINE)

    n_fg, n_un = counts["n_frozen_gi_used"], counts["n_unrounded"]
    n_by_key = {"frozen": n_fg, "gi": n_fg, "unrounded": n_un}
    for key, label, colour, ls, marker in SERIES:
        label = f"{label} ($n$ = {n_by_key[key]})"
        y = support[cols[key]].to_numpy(float)
        ok = np.isfinite(y)
        ax.plot(gains[ok], y[ok], ls=ls, lw=1.2, color=colour, marker=marker,
                ms=3.4, mew=0.6, mfc="white", mec=colour, label=label, zorder=3)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Signal gain ($\\times$)", labelpad=2)
    ax.set_ylabel("Retained tolerance levels", labelpad=2)
    ax.set_title("Tolerance-support scaling", fontsize=FS_TITLE, color=INK, pad=4)
    ax.set_xticks(list(gains))
    ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, lw=0.5, color="#E3E3E3")
    ax.tick_params(length=2.6, pad=1.8)
    # Placed in the clear band between the GI and unrounded series.
    ax.legend(loc="upper left", bbox_to_anchor=(0.02, 0.80), frameon=False,
              handlelength=2.0, handletextpad=0.5, borderpad=0.2,
              labelspacing=0.32, fontsize=FS_LEGEND)


def draw_deviation_panel(ax, summ: pd.DataFrame, fp_dev: pd.DataFrame, notes: dict):
    """(b) Median |deviation| of each summary from its unit-gain value."""
    frozen = summ[summ["mode"].str.lower() == "frozen"]
    gi = summ[summ["mode"].str.upper() == "GI"]

    n_fg, n_un = notes["n_frozen_gi"], notes["n_unrounded"]
    spec = [("total", f"Standard Total ($n$ = {n_fg})", COLOR_FROZEN_LINE, "-", "o"),
            ("average", f"Standard Average (AvgSampEn) ($n$ = {n_fg})",
             COLOR_FROZEN_LINE, "--", "s")]
    for feature, label, colour, ls, marker in spec:
        block = frozen[frozen["feature"] == feature].sort_values("gain")
        ax.plot(block["gain"], block["abs_dev"], ls=ls, lw=1.2, color=colour,
                marker=marker, ms=3.4, mew=0.6, mfc="white", mec=colour,
                label=label, zorder=3)

    gi_total = gi[gi["feature"] == "total"].sort_values("gain")
    ax.plot(gi_total["gain"], gi_total["abs_dev"], ls="-", lw=1.2,
            color=COLOR_GI, marker="s", ms=3.4, mew=0.6, mfc="white",
            mec=COLOR_GI, label=f"GI Total ($n$ = {n_fg})", zorder=4)
    ax.plot(fp_dev["gain"], fp_dev["abs_dev"], ls="--", lw=1.2,
            color=COLOR_UNROUNDED, marker="^", ms=3.4, mew=0.6, mfc="white",
            mec=COLOR_UNROUNDED,
            label=f"Unrounded Total ($n$ = {n_un})", zorder=4)

    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1.0, linscale=0.45)
    ax.set_ylim(-0.35, 1000)
    ax.set_yticks([0, 1, 10, 100, 1000])
    ax.get_yaxis().set_major_formatter(ticker.FuncFormatter(lambda v, _p: f"{v:g}"))
    gains = sorted(summ["gain"].unique())
    ax.set_xticks(gains)
    ax.get_xaxis().set_major_formatter(ticker.FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.set_xlabel("Signal gain ($\\times$)", labelpad=2)
    ax.set_ylabel("Absolute deviation from 1$\\times$ (%)", labelpad=2)
    ax.set_title("Gain sensitivity of profile summaries", fontsize=FS_TITLE,
                 color=INK, pad=4)
    ax.axhline(0, color=ZERO_LINE, lw=0.7, zorder=2)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, lw=0.5, color="#E3E3E3")
    ax.tick_params(length=2.6, pad=1.8)
    ax.annotate("symlog; linear below 1%", xy=(0.03, 0.50),
                xycoords="axes fraction", ha="left", va="top",
                fontsize=FS_ANNOT - 0.5, color=INK_MUTED, style="italic")
    ax.legend(loc="upper left", frameon=False, handlelength=2.0,
              handletextpad=0.5, borderpad=0.2, labelspacing=0.32,
              fontsize=FS_LEGEND - 0.3)


# --------------------------------------------------------------------------
# Panel B
# --------------------------------------------------------------------------
def draw_panel_b(ax_left, ax_right, panel_b: pd.DataFrame, show_ci: bool):
    width = 0.34
    idx = np.arange(len(CHANNEL_ORDER))
    axes = {"gi_sd": ax_left, "gi_total": ax_right}

    est_cols = panel_b[["estimate", "frozen_estimate"]].to_numpy(float)
    lo, hi = np.nanmin(est_cols), np.nanmax(est_cols)
    if show_ci and {"ci_low", "ci_high"} <= set(panel_b.columns):
        lo = min(lo, float(panel_b["ci_low"].min()))
        hi = max(hi, float(panel_b["ci_high"].max()))
    span = hi - lo
    ylim = (lo - 0.16 * span, hi + 0.42 * span)

    for feature, ax in axes.items():
        sub = panel_b[panel_b["feature"] == feature].set_index("channel")
        sub = sub.reindex(CHANNEL_ORDER)

        frozen = sub["frozen_estimate"].to_numpy(float)
        gi = sub["estimate"].to_numpy(float)

        ax.bar(idx - width / 2, frozen, width, color=COLOR_FROZEN,
               edgecolor="#3F3F3F", linewidth=0.55, hatch="////",
               label="Standard", zorder=2)
        ax.bar(idx + width / 2, gi, width, color=COLOR_GI,
               edgecolor="#12384F", linewidth=0.55, label="GI", zorder=2)

        if show_ci and {"ci_low", "ci_high"} <= set(sub.columns):
            ax.errorbar(idx + width / 2, gi,
                        yerr=[gi - sub["ci_low"].to_numpy(float),
                              sub["ci_high"].to_numpy(float) - gi],
                        fmt="none", ecolor="#12384F", elinewidth=0.7,
                        capsize=1.8, capthick=0.7, zorder=4)

        # One asterisk = BH-FDR q < 0.05 in the 16-test cross-channel family.
        for pos, value, q in ((idx - width / 2, frozen, sub["frozen_q_value"].to_numpy(float)),
                              (idx + width / 2, gi, sub["q_value"].to_numpy(float))):
            for xi, vi, qi in zip(pos, value, q):
                if np.isfinite(qi) and qi < 0.05:
                    off = -0.022 * span if vi < 0 else 0.022 * span
                    ax.annotate("*", xy=(xi, vi + off),
                                ha="center",
                                va="top" if vi < 0 else "bottom",
                                fontsize=FS_ANNOT + 1.5, color=INK, zorder=5)

        ax.axhline(0, color=ZERO_LINE, lw=0.8, zorder=3)
        ax.set_xticks(idx)
        ax.set_xticklabels(CHANNEL_ORDER)
        ax.set_ylim(*ylim)
        ax.set_title(FEATURE_TITLES[feature], fontsize=FS_TITLE, color=INK, pad=4)
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, lw=0.5, color="#DCDCDC")
        ax.xaxis.grid(False)
        ax.tick_params(length=2.6, pad=1.8)

    ax_left.set_ylabel("FM-vs-IN effect\n(standardised)", labelpad=2)
    ax_right.tick_params(labelleft=True)

    handles = [
        Patch(facecolor=COLOR_FROZEN, edgecolor="#3F3F3F", hatch="////",
              linewidth=0.55, label="Standard"),
        Patch(facecolor=COLOR_GI, edgecolor="#12384F", linewidth=0.55, label="GI"),
        Line2D([], [], linestyle="none", marker="$*$", color=INK,
               markersize=6, label="BH-FDR $q$ < 0.05"),
    ]
    ax_left.legend(handles=handles, loc="upper left", frameon=False,
                   handlelength=1.5, handletextpad=0.5, borderpad=0.2,
                   labelspacing=0.3, fontsize=FS_LEGEND, ncol=3,
                   columnspacing=1.0)


# --------------------------------------------------------------------------
# Assembly
# --------------------------------------------------------------------------
def _new_figure():
    fig = plt.figure(figsize=(FIG_W, FIG_H), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.035, h_pad=0.035,
                                    wspace=0.045, hspace=0.10)
    return fig


def _label(ax, text):
    ax.text(-0.215, 1.10, text, transform=ax.transAxes, fontsize=FS_PANEL,
            fontweight="bold", color=INK, ha="left", va="bottom", clip_on=False)


def build_figure(data, checks, panel_b, mode, show_ci):
    """Collinearity layout: (a) amplitude association, (b) consequence."""
    fig = _new_figure()
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.82])

    ax_a_l = fig.add_subplot(gs[0, 0])
    ax_a_r = fig.add_subplot(gs[0, 1])
    ax_b_l = fig.add_subplot(gs[1, 0])
    ax_b_r = fig.add_subplot(gs[1, 1])

    draw_panel_a(ax_a_l, ax_a_r, data, checks, mode)
    draw_panel_b(ax_b_l, ax_b_r, panel_b, show_ci)

    _label(ax_a_l, "(a)")
    _label(ax_b_l, "(b)")
    return fig


def build_mechanism_figure(support, counts, summ, fp_dev, notes, panel_b, show_ci):
    """Mechanism layout: (a) cause, (b) cheap fixes fail, (c) consequence."""
    fig = _new_figure()
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.82])

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c_l = fig.add_subplot(gs[1, 0])
    ax_c_r = fig.add_subplot(gs[1, 1])

    draw_support_panel(ax_a, support, counts)
    draw_deviation_panel(ax_b, summ, fp_dev, notes)
    draw_panel_b(ax_c_l, ax_c_r, panel_b, show_ci)

    _label(ax_a, "(a)")
    _label(ax_b, "(b)")
    _label(ax_c_l, "(c)")
    return fig


def _rel(path: Path) -> str:
    """Repo-relative where possible, so the audit is not machine-specific."""
    root = REPO_ROOT
    try:
        return str(Path(path).resolve().relative_to(root))
    except ValueError:
        return str(path)


def write_audit(path, checks, panel_b, mode, outputs, show_ci, *,
                layout="collinearity", support=None, fp_dev=None,
                counts=None, notes=None):
    now = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    mode_note = {
        "loglog": "raw units on log10 axes (no standardisation); independent "
                  "y range per panel, shared x range",
        "zlog": "z-scores of log10(value); identical x and y limits",
        "zlinear": "z-scores of raw values; identical x and y limits",
    }[mode]
    if layout == "mechanism":
        base = float(support.loc[support["gain"] == BASELINE_GAIN,
                                 "frozen_n_retained"].iloc[0])
        lines = [
            "# Figure 1 audit -- mechanism layout", "",
            f"Generated {now} by `src/make_figure1.py --layout mechanism`.",
            "",
            "## Source files", "",
            f"- Panel (a) standard/GI support: `{_rel(INVARIANCE_WINDOWS_PATH)}`",
            f"- Panel (a) unrounded support, panel (b) unrounded deviation: "
            f"`{_rel(FULL_PRECISION_PATH)}`",
            f"- Panel (b) standard/GI deviation: `{_rel(AVG_GAIN_SUMMARY_PATH)}`",
            f"- Panel (c): `{_rel(CROSS_CHANNEL_PATH)}`",
            "",
            "## Columns used", "",
            f"- gain: `{COL_SCALE}` (window table) / `gain` (full-precision table)",
            "- support size: `frozen_n_retained`, `gi_n_retained`, `fp_n_retained`",
            "- deviation: `median_relative_difference_percent` (standard/GI), "
            "`fp_total` (unrounded)",
            "- Panel (c): `estimate`, `frozen_estimate`, `q_value`, `frozen_q_value` "
            "(legacy source-column names retained internally)",
            "",
            "## Window sets behind each series", "",
            f"- standard and GI: n = {counts['n_frozen_gi_used']} windows "
            f"(of {counts['n_frozen_gi_available']} in the gain-stress set)",
            f"- unrounded control: n = {counts['n_unrounded']} windows",
            f"- the unrounded windows are a strict subset of the gain-stress "
            f"windows: {counts['unrounded_is_subset']}",
            f"- common-subset mode: {counts['common_subset']}",
            "",
            ("All three series rest on the same 12 windows."
             if counts["common_subset"] else
             "The series therefore do NOT rest on identical window sets. "
             "Restricting standard and GI to the same 12 windows "
             "(`--common-subset`) shifts the support medians by under 3%, the "
             "fitted standard slope from 0.976 to 0.981, and the standard Total "
             "deviation at 7.1x from 630.6% to 627.4%; the GI and unrounded "
             "series are exactly zero either way. The comparison is therefore "
             "insensitive to the choice, but the n for each series is stated "
             "in the legend and should be stated in the caption."),
            "",
            f"- panel (b) standard/GI source: {notes['frozen_gi_source']}",
            ]
        if "reduction_validated_drift_pp" in notes:
            lines.append(
                f"- restricted reduction validated against the published "
                f"summary on the full set; worst drift "
                f"{notes['reduction_validated_drift_pp']:.3g} percentage points")
        lines += [
            "",
            "## Panel (a) -- retained tolerance-support size, median over windows",
            "",
            "| gain | standard | standard / 1x | GI | unrounded |",
            "|---:|---:|---:|---:|---:|",
        ]
        for _, row in support.iterrows():
            lines.append(
                f"| {row['gain']:g} | {row['frozen_n_retained']:,.0f} | "
                f"{row['frozen_n_retained'] / base:.2f} | "
                f"{row['gi_n_retained']:,.0f} | "
                f"{row['unrounded_n_retained']:,.0f} |")
        lines += [
            "",
            "The standard support tracks gain almost exactly one-for-one; the GI "
            "and unrounded supports are constant. Only medians over the "
            "existing gain-stress windows are computed -- no model, no "
            "estimation, no new statistic.",
            "",
            "## Panel (b) -- median |deviation| from the unit-gain value", "",
            f"- standard and GI series read verbatim from "
            f"`{_rel(AVG_GAIN_SUMMARY_PATH)}` (absolute value taken for display).",
            "- unrounded series reduced from the full-precision per-window table "
            "and checked against its published verdict of exactly zero; the "
            f"script stops if it exceeds 1e-6%. Observed worst case: "
            f"{fp_dev['abs_dev'].max():.3g}%.",
            "- y axis is symlog (linear below 1%, log above) so an exact zero "
            "and a 600% deviation can share one axis.",
            "",
            "## Panel (c) values (read verbatim; no model or FDR was recomputed)",
            "",
            "| feature | channel | GI estimate | GI q | standard estimate | standard q |",
            "|---|---|---:|---:|---:|---:|",
        ]
        for _, row in panel_b.iterrows():
            lines.append(
                f"| {row['feature']} | {row['channel']} | {row['estimate']:+.4f} | "
                f"{row['q_value']:.3g} | {row['frozen_estimate']:+.4f} | "
                f"{row['frozen_q_value']:.3g} |")
        lines += [
            "",
            "q-values are BH-FDR within the 16-test cross-channel family "
            "(4 channels x 2 features x 2 contrasts), as computed upstream. The "
            "four channels are spatially corresponding measurements from the "
            "same women and the same physical windows: cross-channel "
            "consistency, not independent replication.",
            "",
            "## Display-only transformations", "",
            "- Panels (a) and (b): log x; log y in (a), symlog y in (b).",
            "- Absolute value taken in (b) so standard Total's sign change across "
            "gain does not split the series; signed values are in the source CSV.",
            "- No transformation of any kind is applied to panel (c) values.",
            "",
            "## Outputs", "",
        ]
        lines += [f"- `{_rel(p)}`" for p in outputs]
        lines.append("")
        path.write_text("\n".join(lines))
        return

    lines = [
        "# Figure 1 audit -- collinearity layout", "",
        f"Generated {now} by `src/make_figure1.py` (display mode: `{mode}`).",
        "",
        "## Source files", "",
        f"- Panel A, standard-implementation features: `{_rel(checks['frozen_features_path'])}`",
        f"- Panel A, GI features: `{_rel(GI_FEATURES_PATH)}`",
        f"- Panel B: `{_rel(CROSS_CHANNEL_PATH)}`",
        "",
        "## Columns used", "",
        f"- join key: `{COL_WINDOW_ID}`",
        f"- amplitude: `{COL_MAV}`",
        f"- standard Total SampEn (source column): `{COL_FROZEN_TOTAL}`",
        f"- GI Total SampEn: `{COL_GI_TOTAL}`",
        "- Panel B: `estimate`, `frozen_estimate`, `q_value`, `frozen_q_value` "
        "(legacy source-column names retained internally)"
        + (", `ci_low`, `ci_high`" if show_ci else ""),
        "",
        "## Panel A join", "",
        f"- standard-implementation table rows: {checks['frozen_rows']}",
        f"- GI table rows: {checks['gi_rows']}",
        f"- duplicate window_id: standard {checks['frozen_duplicate_ids']}, "
        f"GI {checks['gi_duplicate_ids']}",
        f"- one-to-one join verified: {checks['join_one_to_one']} "
        "(pandas `validate=\"one_to_one\"`)",
        f"- matched windows plotted: {checks['matched_windows']} "
        f"(expected {EXPECTED['n_windows']})",
        "",
        "## Panel A correlations (computed on ORIGINAL untransformed values)", "",
        "| pair | Pearson r | r^2 | Spearman rho |",
        "|---|---:|---:|---:|",
        f"| standard Total vs MAV | {checks['frozen_r']:+.4f} | "
        f"{checks['frozen_r']**2:.4f} | {checks['frozen_rho']:+.4f} |",
        f"| GI Total vs MAV | {checks['gi_r']:+.4f} | "
        f"{checks['gi_r']**2:.4f} | {checks['gi_rho']:+.4f} |",
        "",
        f"Validated targets reproduced within {R_TOLERANCE}: "
        f"standard {EXPECTED['frozen_total_vs_mav_r']:+.3f}, "
        f"GI {EXPECTED['gi_total_vs_mav_r']:+.3f}. The script raises and stops "
        "if either check fails.",
        "",
        "Spearman rho is reported because it is invariant to the display "
        "transform and shows the association is not an artefact of the "
        "amplitude tail.",
        "",
        "## Panel B values (read verbatim; no model or FDR was recomputed)", "",
        "| feature | channel | GI estimate | GI q | standard estimate | standard q |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for _, row in panel_b.iterrows():
        lines.append(
            f"| {row['feature']} | {row['channel']} | {row['estimate']:+.4f} | "
            f"{row['q_value']:.3g} | {row['frozen_estimate']:+.4f} | "
            f"{row['frozen_q_value']:.3g} |")
    lines += [
        "",
        "q-values are BH-FDR within the 16-test cross-channel family "
        "(4 channels x 2 features x 2 contrasts), as computed upstream. The "
        "four channels are spatially corresponding measurements from the same "
        "women and the same physical windows: cross-channel consistency, not "
        "independent replication.",
        "",
        "## Display-only transformations", "",
        f"- Panel A axes: {mode_note}.",
        "- No transformation of any kind is applied to Panel B values.",
        "- Scatter layers are rasterised in the vector outputs; all text, "
        "lines and bars remain vector.",
        "",
        "## Outputs", "",
    ]
    lines += [f"- `{_rel(p)}`" for p in outputs]
    lines.append("")
    path.write_text("\n".join(lines))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--display", choices=("loglog", "zlog", "zlinear"),
                        default="loglog",
                        help="panel (a) axis treatment (default: loglog)")
    parser.add_argument("--ci", action="store_true",
                        help="add 95%% CIs to the GI bars in panel (b)")
    parser.add_argument("--layout", choices=("collinearity", "mechanism"),
                        default="collinearity",
                        help="collinearity: amplitude association + consequence "
                             "(default). mechanism: tolerance support + summary "
                             "sensitivity + consequence.")
    parser.add_argument("--common-subset", action="store_true",
                        help="mechanism layout: restrict the standard and GI "
                             "series to the 12 windows the unrounded control "
                             "uses, so all three curves rest on identical data")
    parser.add_argument("--stem", default=None,
                        help="output basename (defaults per layout)")
    args = parser.parse_args()
    if args.stem is None:
        args.stem = ("Figure1_mechanism_consequence" if args.layout == "collinearity"
                     else "Figure1_support_mechanism"
                          + ("_commonsubset" if args.common_subset else ""))

    font = _pick_font()
    _apply_style(font)

    panel_b = load_panel_b()
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    data = checks = support = summ = fp_dev = counts = notes = None
    if args.layout == "collinearity":
        data, checks = load_panel_a()
        fig = build_figure(data, checks, panel_b, args.display, args.ci)
    else:
        support, counts = load_support_vs_gain(args.common_subset)
        summ, fp_dev, notes = load_deviation_vs_gain(args.common_subset)
        fig = build_mechanism_figure(support, counts, summ, fp_dev, notes,
                                     panel_b, args.ci)

    outputs = []
    for ext, kwargs in (("pdf", {}), ("png", {"dpi": 600}), ("svg", {})):
        out = FIGURE_DIR / f"{args.stem}.{ext}"
        fig.savefig(out, **kwargs)
        outputs.append(out)
    plt.close(fig)

    if args.layout == "collinearity":
        a_src = FIGURE_DIR / "Figure1_panelA_source.csv"
        data.to_csv(a_src, index=False)
        outputs.append(a_src)
    else:
        a_src = FIGURE_DIR / "Figure1_support_source.csv"
        support.to_csv(a_src, index=False)
        outputs.append(a_src)
        b_dev = FIGURE_DIR / "Figure1_deviation_source.csv"
        keep = summ[["feature", "mode", "gain", "median_relative_difference_percent",
                     "abs_dev"]].copy()
        keep["series"] = keep["mode"] + " " + keep["feature"]
        extra = fp_dev.assign(feature="total", mode="unrounded",
                              median_relative_difference_percent=np.nan,
                              series="unrounded total")
        pd.concat([keep, extra], ignore_index=True).to_csv(b_dev, index=False)
        outputs.append(b_dev)

    b_src = FIGURE_DIR / "Figure1_panelB_source.csv"
    panel_b.to_csv(b_src, index=False)
    outputs.append(b_src)

    audit = FIGURE_DIR / f"Figure1_audit_{args.layout}.md"
    write_audit(audit, checks, panel_b, args.display, outputs, args.ci,
                layout=args.layout, support=support, fp_dev=fp_dev,
                counts=counts, notes=notes)
    outputs.append(audit)

    print(f"font: {font}   layout: {args.layout}")
    if args.layout == "collinearity":
        print(f"display mode: {args.display}")
        print(f"matched windows: {checks['matched_windows']}")
        print(f"standard Total vs MAV: r = {checks['frozen_r']:+.4f}  "
              f"r^2 = {checks['frozen_r']**2:.4f}  rho = {checks['frozen_rho']:+.4f}")
        print(f"GI Total vs MAV:     r = {checks['gi_r']:+.4f}  "
              f"r^2 = {checks['gi_r']**2:.4f}  rho = {checks['gi_rho']:+.4f}")
    else:
        base = float(support.loc[support["gain"] == BASELINE_GAIN,
                                 "frozen_n_retained"].iloc[0])
        print("support size vs gain (median over windows):")
        for _, row in support.iterrows():
            print(f"  gain {row['gain']:5.2f}: standard {row['frozen_n_retained']:8.0f} "
                  f"({row['frozen_n_retained']/base:5.2f}x)  "
                  f"GI {row['gi_n_retained']:8.0f}  "
                  f"unrounded {row['unrounded_n_retained']:10.0f}")
        print(f"windows: standard/GI n = {counts['n_frozen_gi_used']} "
              f"(of {counts['n_frozen_gi_available']} available), "
              f"unrounded n = {counts['n_unrounded']}"
              + ("  [common subset]" if counts["common_subset"] else ""))
        print(f"unrounded Total worst |deviation|: {fp_dev['abs_dev'].max():.3g}%")
    for out in outputs:
        print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
