#!/usr/bin/env python3
"""Run the deterministic paper-level pipeline end to end.

Requires ``EHG_RAW_DIR`` (or ``configs/analysis.yaml``) pointing at the
PhysioNet recordings, and ``GI_SAMPEN_UPSTREAM_DIR`` pointing at the prior
study's tree, which supplies the cached parity target and window manifest.
See docs/reproducibility.md before running this.

Usage::

    python scripts/reproduce_paper.py            # everything
    python scripts/reproduce_paper.py --stage parity
    python scripts/reproduce_paper.py --from gain_stress
    python scripts/reproduce_paper.py --dry-run   # print the plan, run nothing

Each stage is the corresponding ``scripts/run_*.py`` / ``scripts/report_*.py``
wrapper, run in a subprocess with this interpreter, in the order the
manuscript results depend on them. A failing stage stops the run; already
tracked results are never overwritten by a stage that did not reach its
``write_csv`` call, so a partial run is safe to inspect and safe to resume
with ``--from``.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent

# (stage name, script, extra args, ~cost note, what it produces)
STAGES = [
    ("parity", "run_parity_check.py", [],
     "~1 min", "results/gi/tables/parity_check.csv -- Step 0, must pass first"),
    ("select_grid", "run_select_grid.py", [],
     "~3-5 min", "results/gi/tables/grid_definition.json"),
    ("extract_gi", "run_extract_gi.py", [],
     "~30-45 min", "work/gi_features_60s_EHG9.csv"),
    ("gain_stress", "run_gain_stress.py", ["--extra-scales"],
     "~2 min", "results/gi/tables/invariance_summary.csv, "
               "work/invariance_window_results.csv"),
    ("avg_gain_report", "report_avg_sampen_gain.py", [],
     "<1 min", "results/gi/tables/avg_sampen_gain_summary.csv (needs gain_stress)"),
    ("full_precision", "run_full_precision_control.py", [],
     "~1 min (12-window subset)", "results/gi/tables/full_precision_gain_test.csv"),
    ("resolution_ablation", "run_resolution_ablation.py", [],
     "~5-10 min", "work/resolution_ablation_features.csv"),
    ("resolution_report", "report_resolution_ablation.py", [],
     "<1 min", "results/gi/tables/resolution_ablation_summary.csv"),
    ("gain_resolution", "run_gain_resolution.py", [],
     "~5-10 min", "results/gi/tables/gain_resolution_sensitivity.csv"),
    ("gain_resolution_report", "report_gain_resolution.py", [],
     "<1 min", "results/gi/audits/gain_resolution_sensitivity_report.md"),
    ("statistics", "run_statistics.py", [],
     "~5 min", "results/gi/tables/gi_model_results.csv -- needs statsmodels"),
    ("gate", "run_gate_report.py", [],
     "<1 min", "results/gi/tables/gate_summary.csv -- the decision"),
    ("cross_channel", "run_cross_channel_analysis.py", [],
     "~2 hours (EHG9-EHG12)", "work/cross_channel_gi_features.csv"),
    ("cross_channel_report", "report_cross_channel.py", [],
     "<1 min", "results/gi/tables/cross_channel_replication_summary.csv -- Table 1"),
    ("amplitude", "run_amplitude_collinearity.py", [],
     "<1 min", "results/gi/tables/amplitude_collinearity.csv"),
    ("figure1", "make_figure1.py", [],
     "<1 min", "results/gi/figures/Figure1_support_mechanism.{pdf,png,svg}"),
]

STAGE_NAMES = [name for name, *_ in STAGES]


def _run(script: str, extra_args: list[str], *, dry_run: bool) -> None:
    cmd = [sys.executable, str(SCRIPTS_DIR / script), *extra_args]
    print(f"$ {' '.join(cmd)}", flush=True)
    if dry_run:
        return
    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise SystemExit(
            f"stage failed: {script} (exit {result.returncode}). "
            "Fix the underlying problem and resume with "
            f"--from {script.replace('run_', '').replace('.py', '')}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stage", choices=STAGE_NAMES,
                        help="Run exactly one stage.")
    parser.add_argument("--from", dest="from_stage", choices=STAGE_NAMES,
                        help="Resume from this stage onward (inclusive).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the commands that would run, without running them.")
    args = parser.parse_args()

    if args.stage:
        plan = [s for s in STAGES if s[0] == args.stage]
    elif args.from_stage:
        start = STAGE_NAMES.index(args.from_stage)
        plan = STAGES[start:]
    else:
        plan = STAGES

    print("Plan:")
    for name, script, extra, cost, produces in plan:
        flag = f" {' '.join(extra)}" if extra else ""
        print(f"  {name:<24} scripts/{script}{flag}")
        print(f"  {'':<24} {cost} -> {produces}")
    print()

    started = time.monotonic()
    for name, script, extra, _cost, _produces in plan:
        print(f"=== {name} ===")
        _run(script, extra, dry_run=args.dry_run)
    elapsed = time.monotonic() - started

    if args.dry_run:
        print("Dry run only; nothing was executed.")
    else:
        print(f"Done in {elapsed / 60:.1f} min.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
