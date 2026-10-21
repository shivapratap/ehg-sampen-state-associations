#!/usr/bin/env python3
"""Conditional log-RMS analysis of GI Total and frozen Total SampEn (EHG9) -> Supplementary S8.

Stages (run one, or `all`):
  extract  needs RAW data + wfdb. One pass over the 4,083 EHG9 windows; writes
             results/conditional/windows_gi_rms_frozen_full.csv   (GI/frozen Total, RMS, ...)
             results/conditional/window_shape_descriptors.csv     (crest factor, excess kurtosis)
  refit    UC/FM-vs-IN models for GI and frozen Total, with and without log RMS
             -> refit_gi_total_rms_results.csv
  within   same, with within-recording-centred log RMS (+ Spearman checks) -> within_recording_results.txt
  shape    GI Total with crest factor / kurtosis / log RMS covariates -> shape_models_output.txt
  lowo     leave-one-woman-out refits (42 x 4 fits) -> lowo_results.csv

Everything except `extract` uses only the saved window table (numpy/pandas/scipy; no raw data).
Expected: GI Total UC coefficient -0.113 -> +0.275 after log-RMS adjustment;
          frozen Total UC 0.409 -> 0.020.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import kurtosis, spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from reml import fit_table, load_windows, reml_fit, zscore  # noqa: E402

OUT = ROOT / "results" / "conditional"
WINDOWS = OUT / "windows_gi_rms_frozen_full.csv"
SHAPE = OUT / "window_shape_descriptors.csv"


def extract():
    from gi_sampen import gi_profile as gp, paths
    frame = paths.load_analysis_windows()
    rows, shape = [], []
    for row, w in paths.iter_windows(frame):
        g, f = gp.gi_summaries_rounded(w), gp.frozen_summaries(w)
        rows.append({"window_id": row.window_id, "recording_id": row.recording_id,
                     "window_class": row.window_class, "patient_id": row.patient_id,
                     "window_start_s": row.window_start_s, "rms": gp.rms(w),
                     "gi_total": g["gi_total"], "gi_sd": g["gi_sd"], "gi_average": g["gi_average"],
                     "gi_median": g["gi_median"], "gi_kurtosis": g["gi_kurtosis"],
                     "gi_skewness": g["gi_skewness"], "gi_n_retained": g["gi_n_retained"],
                     "frozen_total": f["frozen_total"], "frozen_sd": f["frozen_sd"]})
        rms, sd, mean = np.sqrt(np.mean(w ** 2)), np.std(w), np.mean(w)
        shape.append({"window_id": row.window_id, "crest": np.max(np.abs(w)) / rms,
                      "excess_kurtosis": kurtosis(w, fisher=True), "abs_mean_ratio": abs(mean) / sd,
                      "rms_sq_check": rms ** 2 - (sd ** 2 + mean ** 2)})
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(WINDOWS, index=False)
    pd.DataFrame(shape).to_csv(SHAPE, index=False)
    print("wrote", WINDOWS, SHAPE)


def refit(out_name="refit_gi_total_rms_results.csv"):
    df = load_windows(WINDOWS)
    print("n", len(df), "women", df.patient_id.nunique(), "recordings", df.recording_id.nunique(),
          df.window_class.value_counts().to_dict())
    tabs = [fit_table(df, o, a) for o in ["frozen_total", "gi_total"] for a in (False, True)]
    out = pd.concat(tabs)
    out.to_csv(OUT / out_name, index=False)
    print(out.to_string(float_format=lambda v: f"{v:.5f}"))
    return out


def within():
    df = load_windows(WINDOWS)
    df["lrms"] = np.log(df.rms)
    df["lrms_within"] = df.lrms - df.groupby("recording_id").lrms.transform("mean")
    lines = [f"between-recording SD of log RMS: {df.groupby('recording_id').lrms.mean().std():.3f}; "
             f"within-recording SD: {df.lrms_within.std():.3f}"]
    for o in ["gi_total", "frozen_total"]:
        for cov in ["lrms", "lrms_within"]:
            X = np.column_stack([np.ones(len(df)), df.uc, df.fm, zscore(df[cov])])
            b, se, _, _ = reml_fit(zscore(df[o]), X, pd.factorize(df.patient_id)[0], pd.factorize(df.recording_id)[0])
            lines.append(f"{o:12s} cov={cov:12s} UC={b[1]:+.4f} (SE {se[1]:.4f}, CI {b[1]-1.96*se[1]:+.3f} to "
                         f"{b[1]+1.96*se[1]:+.3f}) FM={b[2]:+.4f} (SE {se[2]:.4f}) cov_coef={b[3]:+.4f}")
    inn = df[df.window_class == "baseline"]
    lines.append(f"IN windows Spearman GI Total vs within-recording log RMS: {spearmanr(inn.gi_total, inn.lrms_within)[0]:.3f}")
    lines.append(f"all-window Spearman GI Total vs within-recording log RMS: {spearmanr(df.gi_total, df.lrms_within)[0]:.3f}")
    (OUT / "within_recording_results.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def shape():
    df = load_windows(WINDOWS).merge(pd.read_csv(SHAPE), on="window_id")
    df["lr"], df["crest_z"], df["kurt_z"] = (zscore(np.log(df.rms)), zscore(df.crest), zscore(df.excess_kurtosis))
    inn = df[df.window_class == "baseline"]
    lines = ["IN Spearman: GI~RMS %.3f, GI~crest %.3f, GI~kurt %.3f, RMS~crest %.3f, RMS~kurt %.3f" % (
        spearmanr(inn.gi_total, inn.rms)[0], spearmanr(inn.gi_total, inn.crest)[0],
        spearmanr(inn.gi_total, inn.excess_kurtosis)[0], spearmanr(inn.rms, inn.crest)[0],
        spearmanr(inn.rms, inn.excess_kurtosis)[0])]
    y = zscore(df.gi_total)
    for name, cov in {"none": [], "lr": ["lr"], "crest": ["crest_z"], "kurt": ["kurt_z"],
                      "lr+crest+kurt": ["lr", "crest_z", "kurt_z"]}.items():
        X = np.column_stack([np.ones(len(df)), df.uc, df.fm] + [df[c].values for c in cov])
        b, se, _, _ = reml_fit(y, X, pd.factorize(df.patient_id)[0], pd.factorize(df.recording_id)[0])
        extra = " ".join(f"{c}={b[3+i]:+.3f}" for i, c in enumerate(cov))
        lines.append(f"GI_total {name:15s} UC={b[1]:+.3f} (SE {se[1]:.3f})  {extra}")
    (OUT / "shape_models_output.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def lowo(limit=None):
    df = load_windows(WINDOWS)
    lrz = zscore(np.log(df.rms))
    Z = {o: zscore(df[o]) for o in ["gi_total", "frozen_total"]}  # standardisation fixed to the full population
    out = []
    women = sorted(df.patient_id.unique())[:limit]   # limit=None -> all 42
    for w in women:
        keep = (df.patient_id != w).values
        d = df[keep]
        row = {"dropped_woman": w, "n_windows": int(len(d)), "n_uc_dropped": int(df[~keep].uc.sum())}
        for o in Z:
            for adj in (False, True):
                X = np.column_stack([np.ones(len(d)), d.uc.values, d.fm.values] + ([lrz[keep]] if adj else []))
                b, *_ = reml_fit(Z[o][keep], X, pd.factorize(d.patient_id)[0], pd.factorize(d.recording_id)[0])
                row[f"{o}_{'adj' if adj else 'unadj'}_UC"] = b[1]
        out.append(row)
        print(w, {k: round(float(v), 3) for k, v in row.items() if k.endswith("_UC")}, flush=True)
    res = pd.DataFrame(out)
    name = "lowo_results.csv" if limit is None else f"lowo_results_first{limit}_smoketest.csv"
    res.to_csv(OUT / name, index=False)   # a limited run never overwrites the full result
    for c in [c for c in res.columns if c.endswith("_UC")]:
        print(f"{c}: {res[c].min():+.3f} to {res[c].max():+.3f}")


STAGES = {"extract": extract, "refit": refit, "within": within, "shape": shape, "lowo": lowo}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=[*STAGES, "all"], help="'all' = refit, within, shape, lowo (not extract)")
    ap.add_argument("--limit", type=int, help="lowo only: refit just the first N women (quick smoke test; writes a separate file)")
    a = ap.parse_args()
    for s in (["refit", "within", "shape", "lowo"] if a.stage == "all" else [a.stage]):
        print(f"\n=== {s} ===")
        STAGES[s](a.limit) if s == "lowo" else STAGES[s]()
