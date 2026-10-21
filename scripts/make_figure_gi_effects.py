#!/usr/bin/env python3
"""Figure 4: GI Total and GI profile SD, UC-minus-IN and FM-minus-IN, EHG9-EHG12.

Reads results/manuscript/corrected_family/corrected_cross_channel_all_16.csv (no model is fitted).
Filled markers: q < 0.05 in the 16-test Benjamini-Hochberg family (column q_16test_family);
open markers: q >= 0.05. Bars: model-based Wald 95% CI.

Usage: python scripts/make_figure_gi_effects.py [--out DIR]
Writes: F5_GI_state_associations.pdf (vector), .png (600 dpi) and .csv (plotted values).
"""
import argparse
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CH = ["EHG9", "EHG10", "EHG11", "EHG12"]
COL = {"EHG9": "#0b5d8f", "EHG10": "#c4510b", "EHG11": "#12875f", "EHG12": "#b45a95"}
MRK = {"EHG9": "o", "EHG10": "s", "EHG11": "^", "EHG12": "D"}
SUMMARY = [("gi_total", "GI Total"), ("gi_sd", "GI profile SD")]
CONTRAST = [("contraction_vs_baseline", "(a) UC versus IN"), ("fetal_movement_vs_baseline", "(b) FM versus IN")]

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.facecolor": "white"})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    ap.add_argument("--src", type=Path, default=root / "results/manuscript/corrected_family/corrected_cross_channel_all_16.csv")
    ap.add_argument("--out", type=Path, default=root / "results/manuscript/figures")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(a.src, newline="", encoding="utf-8")))
    assert len(rows) == 16, len(rows)
    get = {(r["feature"], r["contrast"], r["channel"]): r for r in rows}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), sharey=True)
    fig.subplots_adjust(left=0.115, right=0.985, top=0.845, bottom=0.30, wspace=0.07)
    out_rows = []
    off = {c: (1.5 - i) * 0.17 for i, c in enumerate(CH)}  # EHG9 at the top of each group
    for ax, (cn, ct) in zip(axes, CONTRAST):
        for gi, (fk, fl) in enumerate(SUMMARY):
            yc = 1 - gi
            for c in CH:
                r = get[(fk, cn, c)]
                est, lo, hi, q = (float(r[k]) for k in ("estimate", "ci_low", "ci_high", "q_16test_family"))
                y = yc + off[c]
                sig = q < 0.05
                ax.plot([lo, hi], [y, y], color=COL[c], lw=1.3, solid_capstyle="butt", zorder=2)
                ax.plot(est, y, marker=MRK[c], ms=5.6, mfc=COL[c] if sig else "white", mec=COL[c], mew=1.2, zorder=3)
                out_rows.append([ct[4:], fl, c, f"{est:.4f}", f"{lo:.4f}", f"{hi:.4f}", f"{q:.3g}", "filled" if sig else "open"])
        ax.axvline(0, color="#555555", lw=0.8, ls=(0, (3, 2)), zorder=1)
        ax.set_title(ct, loc="left", fontsize=8, weight="bold", pad=4)
        ax.set_ylim(-0.62, 1.62)
        ax.grid(axis="x", color="#e6e6e6", lw=0.6, zorder=0)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=6.5, length=3)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.set_xlabel("state-minus-IN coefficient (model-standardised units)", fontsize=6.8)
    axes[0].set_yticks([1, 0]); axes[0].set_yticklabels([s[1] for s in SUMMARY], fontsize=7.2)
    axes[0].tick_params(axis="y", length=0); axes[1].tick_params(axis="y", length=0)
    fig.suptitle("Gain-invariant (GI) entropy summaries: state-minus-IN associations on four simultaneous channels",
                 fontsize=8, weight="bold", x=0.115, ha="left", y=0.985)
    h1 = [Line2D([], [], color=COL[c], marker=MRK[c], ms=5.6, mfc=COL[c], mec=COL[c], lw=1.3, label=c + (" (primary)" if c == "EHG9" else "")) for c in CH]
    h2 = [Line2D([], [], color="#333333", marker="o", ms=5.6, mfc="#333333", mec="#333333", lw=0, label="filled marker: q < 0.05 (16-test BH family)"),
          Line2D([], [], color="#333333", marker="o", ms=5.6, mfc="white", mec="#333333", mew=1.2, lw=0, label="open marker: q \u2265 0.05"),
          Line2D([], [], color="#333333", lw=1.3, label="bar: 95% CI (model-based Wald)")]
    kw = dict(loc="lower center", fontsize=6.8, frameon=False, columnspacing=2.0, handletextpad=0.5)
    fig.legend(handles=h1, ncol=4, bbox_to_anchor=(0.55, 0.075), title="Channel", title_fontsize=6.8, alignment="center", **kw)
    fig.legend(handles=h2, ncol=3, bbox_to_anchor=(0.55, 0.0), **kw)
    fig.savefig(a.out / "F5_GI_state_associations.pdf")
    fig.savefig(a.out / "F5_GI_state_associations.png", dpi=600)
    with open(a.out / "F5_GI_state_associations.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["contrast", "summary", "channel", "estimate", "ci_low", "ci_high", "q_16test", "marker"]); w.writerows(out_rows)
    print("wrote", a.out / "F5_GI_state_associations.pdf")


if __name__ == "__main__":
    main()
