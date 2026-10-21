"""Standard finite-resolution parity.

Two independent checks:

1. ``frozen_summaries`` agrees with the released ``sampen-profile`` package
   (a verification dependency only -- see docs/provenance.md) on
   deterministic synthetic signals, to floating-point tolerance. This test
   needs no external data and always runs.

2. ``frozen_summaries`` agrees with the prior study's CACHED feature table on
   real analysis windows. This is the manuscript's own parity check (200
   windows, max abs deviation ~1e-15) and needs the upstream tree and the raw
   recordings, so it is marked ``requires_data`` and skipped unless those are
   present.
"""
from __future__ import annotations

import numpy as np
import pytest

from gi_sampen import gi_profile, paths

pytest.importorskip("sampen_profile", reason=(
    "sampen-profile is a test-only verification dependency; install the "
    "extra: pip install git+https://github.com/shivapratap/sampen-profile.git@263b72f7b906208241291c0579e8ec103fb52d79"
))
from sampen_profile import sample_entropy_profile as upstream_profile  # noqa: E402

_KEYMAP = {
    "total": "TotalSampEn",
    "average": "AvgSampEn",
    "median": "MedianSampEn",
    "sd": "StdSampEn",
    "kurtosis": "KurtosisSampEn",
    "skewness": "SkewnessSampEn",
}


def _synthetic_signals():
    rng = np.random.default_rng(0)
    return {
        "ehg_like_cumsum": np.cumsum(rng.standard_normal(1200)) * 1e-2,
        "white_noise": rng.standard_normal(600),
        "sine_plus_noise": (
            np.sin(np.linspace(0, 40, 800)) + 0.1 * rng.standard_normal(800)
        ),
    }


@pytest.mark.parametrize("name,signal", list(_synthetic_signals().items()))
def test_parity_against_released_upstream_package(name, signal):
    """gi_profile.frozen_summaries matches sampen-profile on synthetic data."""
    local = gi_profile.frozen_summaries(signal)
    upstream = upstream_profile(signal, m=gi_profile.EMBEDDING_DIMENSION)

    assert local["frozen_n_retained"] == upstream["n_r_points"], name

    for local_key, upstream_key in _KEYMAP.items():
        local_value = local[f"frozen_{local_key}"]
        upstream_value = upstream[upstream_key]
        assert local_value == pytest.approx(upstream_value, rel=1e-9, abs=1e-12), (
            f"{name}: {local_key} local={local_value!r} upstream={upstream_value!r}"
        )


def _upstream_available() -> bool:
    try:
        paths.upstream_dir()
        paths.raw_dir()
        return paths.manifest_path().is_file() and paths.features_path().is_file()
    except RuntimeError:
        return False


@pytest.mark.requires_data
@pytest.mark.skipif(not _upstream_available(), reason=(
    "Needs the prior study's tree and the PhysioNet raw recordings. See "
    "docs/reproducibility.md."
))
def test_parity_against_cached_manuscript_values():
    """Reproduces the manuscript's own 200-window parity check.

    Matches the tolerance reported in the facts pack: max absolute deviation
    on the order of floating-point noise (~1e-15), not an approximate bound.
    """
    frame = paths.load_analysis_windows(limit=200)
    cached = pd_read_cached_features(frame["window_id"].tolist())

    max_abs_diff = 0.0
    for row, window in paths.iter_windows(frame, progress=False):
        computed = gi_profile.frozen_summaries(window)
        reference = cached.loc[row.window_id]
        for key, column in paths.FROZEN_COLUMNS.items():
            diff = abs(computed[f"frozen_{key}"] - reference[column])
            max_abs_diff = max(max_abs_diff, diff)

    assert max_abs_diff < 1e-9, f"parity drifted: max abs diff {max_abs_diff:.3e}"


def pd_read_cached_features(window_ids):
    import pandas as pd

    frame = pd.read_csv(paths.features_path())
    frame = frame.set_index("window_id")
    return frame.loc[window_ids]
