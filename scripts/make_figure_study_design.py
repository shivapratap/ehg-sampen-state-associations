#!/usr/bin/env python3
"""Figure 1 (study design and feature pipeline) as vector artwork.

Schematic only: it reads no data and fits no model. The two occupied-bin figures in
panel (b) (about 80 for the fixed 10^-3 grid in signal units, about 5,000 for SD units)
are the documented bin counts in Supplementary Table S3D.
Electrode facts (PhysioNet EHGDB v1.0.0 description): 16 monopolar electrodes in a 4x4 grid,
17.5 mm between centres; electrodes 9-12 form the third vertical line on the median axis.

Usage: python scripts/make_figure_study_design.py [--out DIR]
Writes: F1_study_design.pdf (vector) and F1_study_design.png (600 dpi preview).
"""
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, Circle, Polygon

UC, FM, IN_ = "#d9483b", "#4f55c2", "#4f9a5b"
TEAL, AMBER, SLATE, GREY = "#0b5d5b", "#d98a00", "#33414d", "#8a8f94"
FS, FSS, FST = 6.2, 5.6, 7.6

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42,
                     "mathtext.default": "regular", "savefig.facecolor": "white"})


def box(ax, x, y, w, h, text, fc="white", ec=SLATE, lw=0.8, fs=FS, color="black", ls="-", weight="normal", r=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, ls=ls))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=color, weight=weight, linespacing=1.15)


def arrow(ax, x0, y0, x1, y1, color=SLATE, lw=0.8, ls="-"):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, ls=ls, shrinkA=0, shrinkB=0, mutation_scale=6))


def panel(fig, rect, title, xlim=(0, 100), ylim=(0, 100), equal=False):
    ax = fig.add_axes(rect)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim); ax.axis("off")
    if equal:
        ax.set_aspect("equal")
    ax.text(0, 1.015, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=FST, weight="bold")
    return ax


def panel_a(fig):
    ax = panel(fig, [0.012, 0.445, 0.265, 0.50], "(a) 16-electrode array", xlim=(0, 72), ylim=(0, 100), equal=True)
    xs, ys, r = [11, 26, 41, 56], [76, 64, 52, 40], 4.4
    n = 1
    for cx in xs:
        for cy in ys:
            med = cx == 41
            ax.add_patch(Circle((cx, cy), r, fc=TEAL if med else "white", ec=TEAL if med else GREY, lw=1.0 if med else 0.9))
            ax.text(cx, cy, str(n), ha="center", va="center", fontsize=FS, weight="bold" if med else "normal",
                    color="white" if med else "#333333")
            n += 1
    ax.add_patch(Circle((41, 76), r + 1.6, fc="none", ec=TEAL, lw=0.9))
    ax.text(47.6, 81.5, "\u2605", color=TEAL, fontsize=8, ha="left", va="center")
    ax.plot([41, 41], [87, 30], color=TEAL, lw=0.8, ls=(0, (4, 3)), zorder=0)
    ax.text(41, 94, "median axis", ha="center", fontsize=FSS, color=TEAL)
    ax.text(36, 25, "fundus \u2191   pubic symphysis \u2193", ha="center", fontsize=FSS - 0.4, color="#555555")
    ax.add_patch(Circle((4, 17), 2.1, fc=TEAL, ec=TEAL)); ax.add_patch(Circle((4, 17), 3.2, fc="none", ec=TEAL, lw=0.7))
    ax.text(9.5, 17, "EHG9, primary channel", va="center", fontsize=FSS)
    ax.add_patch(Circle((4, 10.5), 2.1, fc=TEAL, ec=TEAL)); ax.text(9.5, 10.5, "EHG10\u2013EHG12, supportive", va="center", fontsize=FSS)
    ax.add_patch(Circle((4, 4), 2.1, fc="white", ec=GREY)); ax.text(9.5, 4, "not analysed", va="center", fontsize=FSS)
    ax.text(71, 4, "centre-to-centre\n17.5 mm", va="center", ha="right", fontsize=FSS - 0.6, color="#555555", linespacing=1.1)


def panel_b(fig):
    ax = panel(fig, [0.30, 0.445, 0.69, 0.50], "(b) Feature pipeline and the gain-sensitive step")
    f2 = FS - 0.6
    cw = 14.5
    box(ax, 0.5, 70, cw, 22, "monopolar EHG\n200 Hz,\n16 channels", fc="#f4f6f7", fs=f2)
    box(ax, 0.5, 40, cw, 22, "band\n0.3\u20131.0 Hz,\nthen 20 Hz", fc="#f4f6f7", fs=f2)
    box(ax, 0.5, 8, cw, 24, "60-s window x\nUC, FM or IN\n(EHG9\u2013EHG12)", fc="#f4f6f7", fs=f2)
    arrow(ax, 7.7, 70, 7.7, 62); arrow(ax, 7.7, 40, 7.7, 32)
    colx = [19.5, 36.2, 52.9, 69.6]
    w = 14.2
    yo, yg, h = 66, 26, 21
    ax.text(19.5, 95, "Original profile: fixed absolute grid", fontsize=FS, weight="bold", color=SLATE, va="center")
    ax.text(19.5, 55.5, "GI profile: SD-scaled", fontsize=FS, weight="bold", color=TEAL, va="center")
    box(ax, colx[0], yo, w, h, "x unchanged\n(signal units,\nany gain g)", fc="#f4f6f7", ls=(0, (3, 2)), fs=f2)
    box(ax, colx[1], yo, w, h, "template\nChebyshev\ndistances\n(Eq. 2)", fs=f2)
    box(ax, colx[2], yo, w, h, "round to\n$10^{-3}$\nsignal units", fc="#fdf1d6", ec=AMBER, lw=1.5, weight="bold", fs=f2)
    box(ax, colx[3], yo, w, h, "profile:\nTotal,\nprofile SD", fs=f2)
    box(ax, colx[0], yg, w, h, "$\\tilde{x}=x/\\mathrm{SD}(x)$\npopulation SD", fc="#e4f1f0", ec=TEAL, fs=f2)
    box(ax, colx[1], yg, w, h, "template\nChebyshev\ndistances of\n$\\tilde{x}$ (Eq. 2)", fs=f2)
    box(ax, colx[2], yg, w, h, "round to\n$10^{-3}$\nSD units", fc="#e4f1f0", ec=TEAL, lw=1.5, weight="bold", fs=f2)
    box(ax, colx[3], yg, w, h, "GI profile:\nGI Total,\nGI profile SD", fs=f2)
    for row_y in (yo, yg):
        for i in range(3):
            arrow(ax, colx[i] + w, row_y + h / 2, colx[i + 1], row_y + h / 2)
    ax.text(colx[2] + w / 2, yo - 2.5, "gain enters here:\n\u224880 occupied bins", ha="center", va="top", fontsize=FSS, color=AMBER, weight="bold", linespacing=1.1)
    ax.text(colx[2] + w / 2, yg - 2.5, "gain removed:\n\u22485,000 occupied bins", ha="center", va="top", fontsize=FSS, color=TEAL, weight="bold", linespacing=1.1)
    ya, yb = yo + h / 2, yg + h / 2
    ax.plot([15, 17.2, 17.2], [20, 20, ya], color=SLATE, lw=0.8)
    arrow(ax, 17.2, ya, colx[0], ya); ax.plot([17.2, 17.2], [20, yb], color=SLATE, lw=0.8); arrow(ax, 17.2, yb, colx[0], yb)
    ex = 88.5
    box(ax, ex, yg, 11, yo + h - yg, "state\ncontrasts:\n\nUC vs IN\nFM vs IN\n\n(mixed\nmodels)", fc="#f4f6f7", fs=f2)
    arrow(ax, colx[3] + w, ya, ex, ya); arrow(ax, colx[3] + w, yb, ex, yb)
    ax.text(colx[0] + 2, 2.5, "same women, windows and model hierarchy for both profiles", ha="left", va="center", fontsize=FSS, style="italic", color="#444444")


def panel_c(fig):
    ax = panel(fig, [0.012, 0.025, 0.45, 0.335], "(c) Repeated-measures hierarchy and inference")
    labels = ["Woman", "Recording", "Event or interval\n(UC, FM or IN)", "Analysis window\n(60 s)", "Channel features\n(EHG9\u2013EHG12)"]
    hh, gap = 13.0, 4.8
    y0 = 100 - hh - 2
    for i, a in enumerate(labels):
        y = y0 - i * (hh + gap)
        box(ax, 1, y, 40, hh, a, fs=FS - 0.4, weight="bold")
        if i < 4:
            arrow(ax, 21, y, 21, y - gap)
    box(ax, 50, 62, 49, 36, "Inference tiers\n\u2022 mixed model: woman and\n   recording random effects\n\u2022 pooled OLS, woman-cluster CR1\n\u2022 paired woman-level comparison", fc="#f4f6f7", fs=FSS)
    box(ax, 50, 27, 49, 31, "Multiplicity (BH, within family)\n\u2022 34 features per channel, contrast\n\u2022 16-test GI family, chosen\n   post hoc", fc="#f4f6f7", fs=FSS)
    box(ax, 50, 1, 49, 22, "EHG9\u2013EHG12 share one physical\nwindow: spatial evidence,\nnot independent replication", fc="white", ec=GREY, ls=(0, (3, 2)), fs=FSS, color="#444444")
    ax.plot([41, 45.5, 45.5], [91, 91, 80], color=TEAL, lw=0.8, ls=(0, (3, 2)))
    arrow(ax, 45.5, 80, 50, 80, color=TEAL)


def panel_d(fig):
    ax = panel(fig, [0.488, 0.025, 0.50, 0.335], "(d) Window construction (60 s; 40/80 s inset)")
    ax.add_patch(Rectangle((1, 92.5), 3, 5, fc=UC, ec="none")); ax.text(5.5, 95, "UC  (C)", fontsize=FSS, va="center")
    ax.add_patch(Rectangle((21, 92.5), 3, 5, fc=FM, ec="none")); ax.text(25.5, 95, "FM  (fm)", fontsize=FSS, va="center")
    ax.add_patch(Polygon([[44, 95], [46.2, 97.5], [48.4, 95], [46.2, 92.5]], fc=GREY, ec="none")); ax.text(50, 95, "other annotation", fontsize=FSS, va="center")
    ax.add_patch(Rectangle((78, 92.5), 3, 5, fc=IN_, ec="none")); ax.text(82.5, 95, "IN window", fontsize=FSS, va="center")
    arrow(ax, 1, 70, 99, 70, lw=1.0); ax.text(99, 66, "time", ha="right", fontsize=FSS)
    ex = {"uc": 13, "oth": 31, "fm": 49}
    ax.add_patch(Rectangle((ex["uc"] - 2.2, 76), 4.4, 8, fc=UC, ec="none"))
    ax.add_patch(Polygon([[ex["oth"] - 2.6, 80], [ex["oth"], 84.5], [ex["oth"] + 2.6, 80], [ex["oth"], 75.5]], fc=GREY, ec="none"))
    ax.add_patch(Rectangle((ex["fm"] - 2.2, 76), 4.4, 8, fc=FM, ec="none"))
    for k in ex.values():
        ax.plot([k, k], [75, 52], color=GREY, lw=0.6, ls=(0, (3, 3)))
    ax.add_patch(Rectangle((3, 59), 58, 6, fc="white", ec=GREY, hatch="////", lw=0.6)); ax.text(32, 62, " \u00b1120 s guard ", ha="center", va="center", fontsize=FSS - 0.2, bbox=dict(fc="white", ec="none", pad=0.6))
    wy, wh = 40, 12
    box(ax, 6.2, wy, 13.6, wh, "60 s", fc="#f2b1ab", ec=UC, fs=FSS, r=0.4)
    box(ax, 24.2, wy, 13.6, wh, "60 s", fc="#e6e6e6", ec=GREY, fs=FSS, ls=(0, (3, 2)), color=GREY, r=0.4)
    box(ax, 42.2, wy, 13.6, wh, "60 s", fc="#adaee3", ec=FM, fs=FSS, r=0.4)
    ax.text(13, 37, "UC window,\ncentred", ha="center", va="top", fontsize=FSS - 0.4, color=UC, linespacing=1.1)
    ax.text(31, 37, "excluded", ha="center", va="top", fontsize=FSS - 0.4, color=GREY)
    ax.text(49, 37, "FM window,\ncentred", ha="center", va="top", fontsize=FSS - 0.4, color=FM, linespacing=1.1)
    for i in range(3):
        box(ax, 64 + i * 9.4, wy, 9.4, wh, "60 s", fc="#a9d2b0", ec=IN_, fs=FSS - 0.3, r=0.3)
    box(ax, 64 + 3 * 9.4, wy, 7, wh, "<60", fc="white", ec=GREY, ls=(0, (1, 1.5)), fs=FSS - 0.8, color=GREY, r=0.3)
    ax.text(82, 37, "IN windows tiled from\nguard-cleared interval;\nshort remainder unused", ha="center", va="top", fontsize=FSS - 0.4, color=IN_, linespacing=1.1)
    for k, (wd, lab) in enumerate([(26, "80 s"), (19.5, "60 s"), (13, "40 s")]):
        ax.add_patch(Rectangle((2 + (26 - wd) / 2, 14.5 - k * 4.6), wd, 3.6, fc=["#f7c9c5", "#ec7f76", "#f7c9c5"][k], ec="none"))
        ax.text(31, 16.3 - k * 4.6, lab, va="center", fontsize=FSS - 0.8, color="#444444")
    ax.text(40, 10, "sensitivity windows (UC, FM): same centre;\nschematic, not to scale", va="center", fontsize=FSS - 0.6, color="#444444", linespacing=1.1)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    root = Path(__file__).resolve().parents[1]
    ap.add_argument("--out", type=Path, default=root / "results" / "manuscript" / "figures")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(7.2, 5.1))
    panel_a(fig); panel_b(fig); panel_c(fig); panel_d(fig)
    for rect in ([0.004, 0.42, 0.283, 0.576], [0.29, 0.42, 0.706, 0.576], [0.004, 0.004, 0.465, 0.405], [0.478, 0.004, 0.518, 0.405]):
        fig.add_artist(Rectangle(rect[:2], rect[2], rect[3], transform=fig.transFigure, fc="none", ec="#c9ced2", lw=0.6))
    fig.savefig(a.out / "F1_study_design.pdf")
    fig.savefig(a.out / "F1_study_design.png", dpi=600)
    print("wrote", a.out / "F1_study_design.pdf")


if __name__ == "__main__":
    main()
