#!/usr/bin/env python3
"""Run the analysis pipeline, stage by stage.

    python scripts/run_all.py --list                  # show every stage, what it needs, what it writes
    python scripts/run_all.py                         # default: all stages that need NO raw data
    python scripts/run_all.py --only primary,lag1     # chosen stages
    python scripts/run_all.py --from robustness       # resume
    python scripts/run_all.py --with-raw              # also raw-data and feature-engine stages (needs EHG_RAW_DIR)
    python scripts/run_all.py --dry-run

Two groups of stages:
  derived  start from the derived tables in work/ (restore them first with scripts/restore_derived.py,
           or regenerate them with the raw stages). Need numpy/pandas/scipy; most also need statsmodels.
  engine   audits that recompute the 34 features from packaged windows (need the pinned feature engine; run with --with-raw or --only)
  raw      start from the PhysioNet recordings (set EHG_RAW_DIR, needs wfdb). These produce the
           derived tables. `features` additionally needs the pinned amrita-biosignal-feature-engine.

Each stage is an ordinary script; running one by hand is equivalent. A failing stage stops the run.
"""
from __future__ import annotations
import argparse, os, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable
S, SC = "src", "scripts"

# name, group, command (relative to ROOT), needs, ~time, writes
STAGES = [
    # ---------------- raw -> derived ----------------
    ("data", "raw", [S + "/data.py"], "raw", "~2 min", "work/audit/, work/annotations/ (inventory)"),
    ("annotations", "raw", [S + "/annotations.py"], "raw", "<1 min", "work/annotations/"),
    ("windowing", "raw", [S + "/windowing.py", "--window-size", "60"], "raw", "<1 min", "work/windows/, work/windowing/"),
    ("features", "raw", [S + "/features.py", "--window-size", "60", "--channels", "EHG9,EHG10,EHG11,EHG12"],
     "raw + feature-engine", "hours", "work/features/"),
    ("feature_qc", "raw", [S + "/feature_qc.py", "--window-size", "60"], "raw", "~1 min", "work/feature_qc/"),
    ("gi_pipeline", "raw", [SC + "/reproduce_paper.py"], "raw", "~3 h (cross-channel extraction dominates)",
     "work/gi_*.csv, work/cross_channel_gi_features.csv, results/gi/"),
    ("conditional_extract", "raw", [SC + "/run_conditional_amplitude.py", "extract"], "raw", "~40 min",
     "results/conditional/windows_gi_rms_frozen_full.csv, window_shape_descriptors.csv"),
    ("waveforms", "raw", [SC + "/manuscript_audits/export_waveform_package.py"], "raw", "~2 min",
     "results/manuscript/audits/waveforms/ (60 audited windows)"),
    # ---------------- audits that recompute the 34 features: need the pinned feature engine ----------------
    ("gain_audit", "engine", [SC + "/manuscript_audits/run_D_gain.py"], "feature engine", "~10 min",
     "results/manuscript/audits/D/ (gain audit on the 60 packaged windows)"),
    ("gain1_gate", "engine", [SC + "/manuscript_audits/run_D_gain1_gate.py"], "feature engine", "~2 min",
     "results/manuscript/audits/D/gain1_gate_*"),
    # ---------------- derived -> results ----------------
    ("primary", "derived", [S + "/models.py"], "statsmodels", "~10 min", "work/primary/, results/primary/"),
    ("robustness", "derived", [S + "/robustness.py"], "statsmodels", "~15 min", "results/robustness/{dependence,patient_level,overlap,residuals}"),
    ("cross_channel", "derived", [S + "/cross_channel.py"], "statsmodels", "~10 min", "results/cross_channel/"),
    ("entropy", "derived", [S + "/entropy_analysis.py"], "statsmodels", "~5 min", "results/entropy/"),
    ("fm", "derived", [S + "/fm_analysis.py"], "-", "<1 min", "results/fm/"),
    ("gi_statistics", "derived", [SC + "/run_statistics.py"], "statsmodels", "~5 min", "results/gi/tables/gi_model_results.csv"),
    ("gi_reports", "derived", [SC + "/run_gate_report.py"], "-", "<1 min", "results/gi/tables/gate_summary.csv"),
    ("corrected_family", "derived", [SC + "/run_corrected_family.py"], "statsmodels", "~5 min", "results/manuscript/corrected_family/"),
    ("lag1", "derived", [SC + "/run_lag1_diagnostic.py"], "statsmodels", "~5 min", "results/robustness/lag1/"),
    ("conditional", "derived", [SC + "/run_conditional_amplitude.py", "all"], "scipy only", "~25 min (lowo dominates)", "results/conditional/"),
    ("grid_check", "derived", [SC + "/manuscript_audits/run_C3_grid.py"], "scipy only", "<1 min", "results/manuscript/audits/C3/"),
    ("exact_gate", "derived", [SC + "/manuscript_audits/run_exact_gate.py"], "statsmodels", "~10 min", "results/manuscript/audits/exact_gate_*"),
    ("annotation_sensitivity", "derived", [SC + "/manuscript_audits/run_B_sensitivity.py"], "statsmodels", "~5 min", "results/manuscript/audits/B/"),
    ("verify_statsmodels", "derived", [SC + "/verify_statsmodels.py"], "statsmodels", "~1 min", "results/conditional/statsmodels_verification.*"),
    ("figure_inputs", "derived", [SC + "/build_figure_inputs.py"], "pandas", "<1 min", "results/manuscript/adapters/"),
    ("manuscript_tables", "derived", [SC + "/build_manuscript_tables.py"], "stdlib", "<1 s", "results/manuscript/tables/ (T2, S1-S5, F5 sources)"),
    ("manuscript_figures", "derived", [SC + "/make_manuscript_figures.py", "--overwrite"], "matplotlib", "<1 min", "results/manuscript/figures/ (Fig. 2-4, S1, Table 1)"),
    ("tables", "derived", [S + "/tables.py"], "-", "<1 min", "results/tables/"),
    ("figures", "derived", [S + "/figures.py"], "-", "<1 min", "results/figures/"),
    ("supplement", "derived", [SC + "/build_supplement_tables.py"], "stdlib", "<1 s", "results/manuscript/supplement/ (36 checks)"),
]
# Optional heavy analyses (not in the default run)
HEAVY = [("robustness_heavy", "heavy", [S + "/robustness.py", "--include-heavy"], "statsmodels", "hours",
          "balanced-IN (1,000 resamples) and 40/80-s window sensitivity")]
NAMES = [s[0] for s in STAGES + HEAVY]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--only", help="comma-separated stage names")
    ap.add_argument("--from", dest="frm", choices=NAMES, help="start at this stage")
    ap.add_argument("--with-raw", action="store_true", help="include the raw-data stages")
    ap.add_argument("--heavy", action="store_true", help="include the heavy robustness stage")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    allst = STAGES + (HEAVY if (a.heavy or (a.only and "robustness_heavy" in a.only)) else [])
    if a.list:
        for n, g, cmd, needs, t, w in STAGES + HEAVY:
            print(f"{n:<24}{g:<8}{needs:<22}{t:<26}{w}")
        return
    if a.only:
        want = a.only.split(",")
        bad = [w for w in want if w not in NAMES]
        if bad:
            sys.exit(f"unknown stage(s): {bad}. See --list.")
        plan = [s for s in allst if s[0] in want]
    else:
        plan = [s for s in allst if a.with_raw or s[1] not in ("raw", "engine")]
        if a.frm:
            plan = plan[[s[0] for s in plan].index(a.frm):]
    env = dict(os.environ, PYTHONPATH=str(ROOT / S) + os.pathsep + os.environ.get("PYTHONPATH", ""))
    t0 = time.time()
    for n, g, cmd, needs, t, w in plan:
        print(f"\n=== {n}  [{g}; {needs}; {t}] -> {w}", flush=True)
        c = [PY] + [str(ROOT / cmd[0])] + cmd[1:]
        print("$", " ".join(c), flush=True)
        if a.dry_run:
            continue
        if subprocess.run(c, cwd=ROOT, env=env).returncode:
            sys.exit(f"\nstage '{n}' failed; fix it and resume with --from {n}")
    print(f"\nfinished {len(plan)} stage(s) in {(time.time()-t0)/60:.1f} min" + (" (dry run)" if a.dry_run else ""))


if __name__ == "__main__":
    main()
