#!/usr/bin/env python3
"""Canonical Phase E2 entropy, scale, and RMS-adjustment analyses."""
from __future__ import annotations

import importlib.util
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
WORK = ROOT / "work" / "entropy"
RESULTS = ROOT / "results" / "entropy"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def run_descriptive_entropy_amplitude(channels=None):
    """Run the descriptive entropy-amplitude coupling audit."""
    mod = _load("canonical_entropy_amplitude", SRC / "components" / "_entropy_amplitude_impl.py")
    channels = tuple(channels or ("EHG9", "EHG10", "EHG11", "EHG12"))
    out = RESULTS
    out.mkdir(parents=True, exist_ok=True)
    return mod.run(channels=channels, units_dir=ROOT / "work" / "primary", results_dir=out)


def run_sampen_scale_experiment():
    """Run the locked five-factor Total SampEn scale experiment."""
    mod = _load("canonical_sampen_scale", SRC / "components" / "_sampen_scale_impl.py")
    result = mod.main()
    WORK.mkdir(parents=True, exist_ok=True)
    expanded = RESULTS / "scale_sensitivity_window_results.csv"
    if expanded.exists():
        shutil.copy2(expanded, WORK / expanded.name)
    return result


def run_rms_adjusted_entropy():
    """Run corrected physical-window-aware RMS-adjusted entropy models."""
    mod = _load("canonical_rms_adjusted_entropy", SRC / "components" / "_rms_adjusted_entropy_impl.py")
    result = mod.main()
    source = RESULTS / "corrected_amplitude_adjusted_model_results.csv"
    if source.exists():
        shutil.copy2(source, RESULTS / "rms_adjusted_entropy_results.csv")
    source = RESULTS / "entropy_feature_evidence_classification.csv"
    if source.exists():
        shutil.copy2(source, RESULTS / "entropy_evidence_summary.csv")
    return result


def main() -> int:
    WORK.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(parents=True, exist_ok=True)
    run_descriptive_entropy_amplitude()
    raw = Path(__import__("os").environ.get("EHG_RAW_DIR", ROOT / "data" / "raw")).expanduser()
    if (raw / "RECORDS").exists():
        run_sampen_scale_experiment()
    else:
        print("SKIPPED the Total SampEn scale experiment: it re-reads 60 raw windows (needs EHG_RAW_DIR and the "
              "feature engine). Saved results are in results/entropy/scale_sensitivity_*.csv.")
    run_rms_adjusted_entropy()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
