#!/usr/bin/env python3
"""EHG sample-entropy pipeline -- Stage 03: Preprocessing.

Standalone port of the filtering half of pipeline/module2/features.py's
logic (frozen reference; not imported at runtime -- see
the design notes). The old pipeline folds filtering into its
feature-extraction run; this project separates it into its own stage for
clarity, but keeps the identical algorithm: a full-recording, zero-phase
Butterworth SOS bandpass filter (order 4) for the FWH (0.3-1.0 Hz) and
Broad (0.1-3.0 Hz) bands, followed by integer decimation (200 Hz -> 20 Hz,
factor 10) applied AFTER filtering, and only when a window is sliced out.

Filtering the full recording before slicing (rather than filtering each
window independently) is a deliberate, load-bearing design choice carried
over unchanged: it means every manifest window inherits the filter's
transient-free context and does not acquire its own per-window edge
effects. This is also why Stage 02's 15s edge-trim margin exists.

This is a **function library**, not a script with its own persisted
output: Stage 04
imports `filter_recording` and `filtered_window` directly and calls them
once per (recording, channel), exactly mirroring the old pipeline's
approach and avoiding an unnecessary large intermediate artifact (filtered
waveforms for all four channels would be a substantial disk footprint with
no independent use).

Design source (logic ported, not imported): pipeline/module2/features.py
(`filter_recording`, `filtered_window`).
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.signal import butter, sosfiltfilt

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import config as project_config  # noqa: E402


def validate_raw_signal(signal: Any) -> np.ndarray:
    """Return a finite, one-dimensional float64 signal or raise."""
    data = np.asarray(signal, dtype=np.float64)
    if data.ndim != 1:
        raise ValueError("raw signal must be one-dimensional")
    if data.size < 2:
        raise ValueError("raw signal must contain at least two samples")
    if not np.all(np.isfinite(data)):
        raise ValueError("raw signal contains NaN or infinite samples")
    return data


def filter_recording(signal: Any, sampling_frequency_hz: float) -> dict[str, np.ndarray]:
    """Zero-phase filter a full recording into FWH and Broad bands.

    Filtering occurs before window slicing, so every manifest window
    inherits the context of the recording and does not acquire per-window
    edge effects.
    """
    data = validate_raw_signal(signal)
    fs = float(sampling_frequency_hz)
    if not np.isfinite(fs) or fs <= 0:
        raise ValueError("sampling_frequency_hz must be finite and positive")

    filtered: dict[str, np.ndarray] = {}
    for band_name, (low_hz, high_hz) in project_config.BANDS_HZ.items():
        if high_hz >= fs / 2:
            raise ValueError(f"{band_name} upper edge must be below Nyquist")
        sos = butter(
            project_config.FILTER_ORDER, (low_hz, high_hz),
            btype="bandpass", fs=fs, output="sos",
        )
        filtered[band_name] = np.asarray(sosfiltfilt(sos, data), dtype=np.float64)
    return filtered


def filtered_window(
    filtered_recording: np.ndarray,
    start_sample: int,
    end_sample: int,
    *,
    decimation_factor: int = project_config.DECIMATION_FACTOR,
) -> np.ndarray:
    """Slice a Stage 02 manifest interval and decimate its already
    band-limited data."""
    start = int(start_sample)
    end = int(end_sample)
    if start < 0 or end <= start or end > filtered_recording.size:
        raise ValueError(
            f"invalid window bounds [{start}, {end}) for {filtered_recording.size} samples"
        )
    if decimation_factor < 1:
        raise ValueError("decimation_factor must be positive")
    return np.asarray(
        filtered_recording[start:end:decimation_factor], dtype=np.float64,
    )


if __name__ == "__main__":
    print(__doc__)
    print(
        "This stage is a function library consumed by src/04_feature_extraction.py "
        "(filter_recording, filtered_window) -- it has no standalone CLI output. "
        "Run src/04_feature_extraction.py instead."
    )
    sys.exit(0)
