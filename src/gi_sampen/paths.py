#!/usr/bin/env python3
"""Path resolution and the validated window-loading path.

This module is the single choke point for every path the analysis touches. It
was migrated from ``src/gi_sampen/paths.py``; the loading sequence, the filter
call order and the frozen column names are unchanged, because the parity check
in ``tests/test_standard_parity.py`` depends on them being identical.

Window reconstruction, unchanged from the upstream study:

    read one channel from WFDB
    -> zero-phase bandpass the WHOLE recording (never the window)
    -> slice the manifest interval
    -> decimate 200 Hz -> 20 Hz

The upstream study's ``preprocessing.py`` is imported by path so that the
filter design, band edges and decimation factor cannot drift from the
definitions used to produce the cached feature tables. Its
``feature_extraction.py`` is deliberately NOT imported: it pulls in the
vendored ``amrita-biosignal-feature-engine``, which this analysis does not
need.

External locations
------------------
One input lives outside the repository and is never committed (the second is
the repository itself, see below):

``raw_dir``
    The PhysioNet Icelandic 16-electrode EHG recordings (WFDB). See
    ``data/README.md``. Override with ``EHG_RAW_DIR``.

``upstream_dir``
    Where ``src/preprocessing.py`` and the standard pipeline's ``work/`` and
    ``results/`` tables live (window manifest, cached finite-resolution
    feature table = the parity target, event-unit tables, per-channel model
    results). In this repository that is the repository root, which is the
    default; ``GI_SAMPEN_UPSTREAM_DIR`` can point elsewhere.

Both default from ``configs/analysis.yaml``. No absolute path is hard-coded in
this repository; the defaults in that file are relative to the repository root
and are expected to be edited or overridden on any other machine.
"""
from __future__ import annotations

import importlib.util
import os
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

# --- Repository layout -----------------------------------------------------
# src/gi_sampen/paths.py -> parents[0]=gi_sampen, [1]=src, [2]=repo root
REPO_ROOT = Path(__file__).resolve().parents[2]

# Historical names kept so older call sites keep working. In this repository the
# standard pipeline (formerly "the standard pipeline") and the GI stages (formerly "the GI stages")
# share one tree, so all of these are the repository root.
PAPER1_ROOT = REPO_ROOT
PAPER2_ROOT = REPO_ROOT

CONFIG_PATH = REPO_ROOT / "configs" / "analysis.yaml"

WORK_DIR = REPO_ROOT / "work"
# GI outputs live under results/gi/ so they cannot collide with the standard
# pipeline's results/{primary,robustness,cross_channel,entropy,fm,tables}.
RESULTS_DIR = REPO_ROOT / "results" / "gi" / "tables"
AUDIT_DIR = REPO_ROOT / "results" / "gi" / "audits"
FIGURE_DIR = REPO_ROOT / "results" / "gi" / "figures"
MANIFEST_DIR = REPO_ROOT / "data" / "manifests"

for _directory in (WORK_DIR, RESULTS_DIR, AUDIT_DIR, FIGURE_DIR):
    _directory.mkdir(parents=True, exist_ok=True)


# --- Configuration ---------------------------------------------------------

_DEFAULTS = {
    "raw_dir": "data/raw",
    "upstream_dir": ".",
    "channel": "EHG9",
    "cross_channels": ["EHG9", "EHG10", "EHG11", "EHG12"],
    "band": "fwh",
    "expected_samples": 1200,
}


@lru_cache(maxsize=1)
def config() -> dict:
    """Read ``configs/analysis.yaml``, falling back to the built-in defaults.

    PyYAML is optional: if it is absent, or the file is missing, the defaults
    above apply and can still be overridden by environment variables.
    """
    values = dict(_DEFAULTS)
    if CONFIG_PATH.is_file():
        try:
            import yaml  # noqa: PLC0415 - optional dependency

            loaded = yaml.safe_load(CONFIG_PATH.read_text()) or {}
        except ImportError:
            loaded = {}
        paths_section = loaded.get("paths", {}) or {}
        analysis_section = loaded.get("analysis", {}) or {}
        values.update({k: v for k, v in paths_section.items() if v is not None})
        values.update({k: v for k, v in analysis_section.items() if v is not None})
    return values


def _resolve(value: str) -> Path:
    path = Path(str(value)).expanduser()
    return path if path.is_absolute() else (REPO_ROOT / path).resolve()


def upstream_dir() -> Path:
    """Directory of the prior study that supplies the cached inputs."""
    override = os.environ.get("GI_SAMPEN_UPSTREAM_DIR")
    path = _resolve(override) if override else _resolve(config()["upstream_dir"])
    if not path.is_dir():
        raise RuntimeError(
            f"Upstream study directory not found at {path}. Set "
            "GI_SAMPEN_UPSTREAM_DIR or edit configs/analysis.yaml. See "
            "docs/reproducibility.md."
        )
    return path


def raw_dir() -> Path:
    """Directory holding the PhysioNet WFDB recordings."""
    override = os.environ.get("EHG_RAW_DIR")
    path = _resolve(override) if override else _resolve(config()["raw_dir"])
    if not path.is_dir():
        raise RuntimeError(
            f"Raw recording directory not found at {path}. Set EHG_RAW_DIR or "
            "edit configs/analysis.yaml. See data/README.md for how to obtain "
            "the database."
        )
    return path


# Analysis constants, unchanged.
CHANNEL = _DEFAULTS["channel"]
CROSS_CHANNELS = tuple(_DEFAULTS["cross_channels"])
BAND = _DEFAULTS["band"]
EXPECTED_SAMPLES = _DEFAULTS["expected_samples"]


# --- Derived upstream paths ------------------------------------------------
# Callables rather than constants because they must not be resolved at import
# time: importing this module has to succeed on a machine that has neither the
# raw database nor the upstream tree, so that the synthetic-data tests run.

def manifest_path() -> Path:
    return upstream_dir() / "work" / "windows" / "window_manifest_60s.csv"


def features_path(channel: str = CHANNEL) -> Path:
    """Cached finite-resolution feature table; the parity target."""
    return upstream_dir() / "work" / "features" / f"window_features_60s_{channel}.csv"


def units_path(channel: str = CHANNEL) -> Path:
    return upstream_dir() / "work" / "primary" / f"event_unit_analysis_60s_{channel}.csv"


def frozen_model_results_path(channel: str = CHANNEL) -> Path:
    """Per-channel model results from the prior study.

    The ``frozen_`` name is retained only because it appears in committed CSV
    column names (``frozen_estimate``, ``frozen_q_value``); in documentation
    and prose these are the *standard finite-resolution* results. See
    docs/provenance.md.
    """
    return (upstream_dir() / "work" / "primary"
            / f"mixed_model_results_60s_fwh_{channel}.csv")


class _LazyPath:
    """Deferred path, so ``paths.MANIFEST_PATH`` keeps working as a constant
    at call sites while still failing only when actually used."""

    def __init__(self, factory):
        self._factory = factory

    def __fspath__(self) -> str:
        return str(self._factory())

    def __truediv__(self, other):
        return self._factory() / other

    def __str__(self) -> str:
        return str(self._factory())

    def __repr__(self) -> str:
        return f"<lazy {self._factory()}>"


MANIFEST_PATH = _LazyPath(manifest_path)
FEATURES_PATH = _LazyPath(features_path)
UNITS_PATH = _LazyPath(units_path)
FROZEN_MODEL_RESULTS_PATH = {
    channel: _LazyPath(lambda c=channel: frozen_model_results_path(c))
    for channel in CROSS_CHANNELS
}

# Frozen columns used as the parity target. These are column names in a cached
# CSV produced by the prior study; they are not renamed. See docs/provenance.md.
FROZEN_COLUMNS = {
    "total": "fwh_entropy_totalsampen",
    "average": "fwh_entropy_averagesampen",
    "median": "fwh_entropy_mediansampen",
    "sd": "fwh_entropy_sdsampen",
    "kurtosis": "fwh_entropy_kurtosissampen",
    "skewness": "fwh_entropy_skewnesssampen",
}


# --- Signal loading --------------------------------------------------------

def load_paper1_preprocessing():
    """Import the upstream study's preprocessing module by path."""
    root = upstream_dir()
    src = root / "src"
    os.environ.setdefault("EHG_RAW_DIR", str(raw_dir()))
    for entry in (str(root), str(src)):
        if entry not in sys.path:
            sys.path.insert(0, entry)
    spec = importlib.util.spec_from_file_location(
        "paper1_preprocessing", src / "preprocessing.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def read_channel(recording_id: str, channel: str = CHANNEL):
    """One channel of one recording, at the native 200 Hz."""
    import wfdb

    record = wfdb.rdrecord(str(raw_dir() / recording_id), channel_names=[channel])
    signal = np.asarray(record.p_signal[:, 0], dtype=np.float64)
    fs = float(record.fs)
    if not np.isclose(fs, 200.0, rtol=0, atol=1e-9):
        raise ValueError(f"{recording_id}: expected 200 Hz, observed {fs} Hz")
    return signal, fs


def load_analysis_windows(limit: int | None = None, *,
                          channel: str = CHANNEL) -> pd.DataFrame:
    """The 4,083 analysis units, with the manifest columns needed to
    reconstruct each window.

    Restricting to the analysis units (rather than all 4,283 manifest windows)
    keeps this exactly aligned with the prior study's models.

    ``channel`` selects which channel's event-unit file defines the window
    list. In practice this makes no difference -- the window_id sets are
    identical across EHG9-EHG12 -- but each channel's own file is used so that
    fact is never silently assumed.
    """
    manifest = pd.read_csv(manifest_path())
    units = pd.read_csv(units_path(channel),
                        usecols=["window_id", "window_class",
                                 "patient_id", "recording_id"])
    keep = manifest["window_id"].isin(set(units["window_id"]))
    frame = manifest.loc[keep].copy()
    frame = frame.sort_values(["recording_id", "window_start_s", "window_id"],
                              kind="mergesort").reset_index(drop=True)
    if limit is not None:
        frame = frame.head(limit).copy()
    return frame


def iter_windows(frame: pd.DataFrame, *, progress: bool = True,
                 channel: str = CHANNEL):
    """Yield (row, signal) grouped by recording so each record is read and
    filtered once. This is the dominant cost of the whole stage."""
    preprocessing = load_paper1_preprocessing()
    grouped = list(frame.groupby("recording_id", sort=True))
    for index, (recording_id, rows) in enumerate(grouped, 1):
        if progress:
            print(f"[{index:>3}/{len(grouped)}] {recording_id} "
                  f"({len(rows)} windows, channel {channel})", flush=True)
        raw, fs = read_channel(recording_id, channel)
        filtered = preprocessing.filter_recording(raw, fs)[BAND]
        for row in rows.itertuples(index=False):
            window = preprocessing.filtered_window(
                filtered, int(row.start_sample), int(row.end_sample))
            if window.size != EXPECTED_SAMPLES:
                raise RuntimeError(
                    f"{row.window_id}: expected {EXPECTED_SAMPLES} samples, "
                    f"got {window.size}")
            yield row, window


def write_csv(frame: pd.DataFrame, path, *, label: str = "") -> None:
    path = Path(os.fspath(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    print(f"wrote {label or path.name}: {len(frame)} rows -> {path}", flush=True)
