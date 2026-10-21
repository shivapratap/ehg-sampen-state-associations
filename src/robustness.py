"""Canonical Phase D reviewer robustness entry point.

The component scripts are local to this paper pipeline; no external repository
paths are required at runtime.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
COMPONENTS = Path(__file__).resolve().parent / "components" / "robustness_components"


def _run(name: str, args: list[str] | None = None) -> None:
    subprocess.run([sys.executable, str(COMPONENTS / name)] + (args or []), cwd=ROOT, check=True)


def _run_fresh(name: str, final: Path) -> None:
    """dependence.py and patient_level.py refuse to overwrite existing outputs.
    Run them into a fresh scratch directory, then copy the files over `final`
    (the committed copies stay recoverable through git)."""
    scratch = ROOT / "work" / "robustness" / "_fresh" / final.name
    shutil.rmtree(scratch, ignore_errors=True)
    _run(name, ["--output-dir", str(scratch)])
    final.mkdir(parents=True, exist_ok=True)
    for path in scratch.iterdir():
        shutil.copy2(path, final / path.name)


def run_dependence() -> None:
    _run_fresh("dependence.py", ROOT / "results" / "robustness" / "dependence")


def run_patient_level() -> None:
    _run_fresh("patient_level.py", ROOT / "results" / "robustness" / "patient_level")


def run_overlap() -> None:
    _run("overlap.py")


def run_residuals() -> None:
    _run("residuals.py")


def run_balanced_in(iterations: int = 1000, workers: int = 8) -> None:
    _run("balanced_in.py", ["--iterations", str(iterations), "--workers", str(workers)])


def run_window_duration() -> None:
    _run("window_length.py")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--include-heavy", action="store_true")
    args = parser.parse_args()
    for function in (run_dependence, run_patient_level, run_overlap, run_residuals):
        function()
    if args.include_heavy:
        run_balanced_in()
        run_window_duration()


if __name__ == "__main__":
    main()
