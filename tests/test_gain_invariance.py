"""GI invariance: for positive gains, GI profile summaries are unchanged.

Covers the exact gains named in the paper's gain-stress design: 0.25, 0.5,
1.3, 2, 4, 7.1 (the last two are non-power-of-two, deliberately harder).

This is a correctness test on ``gi_profile.py``, not a scientific result --
see the module docstring there. It needs no external data.
"""
from __future__ import annotations

import numpy as np
import pytest

from gi_sampen import gi_profile

GAINS = (0.25, 0.5, 1.3, 2.0, 4.0, 7.1)

# The manuscript's own invariance check uses < 1e-8 relative tolerance
# (facts_pack.md section 5.2: "every GI feature was invariant to < 1e-8
# relative tolerance at every scale"). Matched here, not loosened.
RELATIVE_TOLERANCE = 1e-8


def _synthetic_windows():
    rng = np.random.default_rng(1)
    return {
        "ehg_like_cumsum_n1200": np.cumsum(rng.standard_normal(1200)) * 1e-2,
        "white_noise_n600": rng.standard_normal(600),
        "sine_plus_noise_n800": (
            np.sin(np.linspace(0, 40, 800)) + 0.1 * rng.standard_normal(800)
        ),
    }


@pytest.mark.parametrize("name,signal", list(_synthetic_windows().items()))
@pytest.mark.parametrize("gain", GAINS)
def test_gi_summaries_rounded_invariant_to_gain(name, signal, gain):
    """gi_summaries_rounded(c * x) == gi_summaries_rounded(x) for c > 0."""
    base = gi_profile.gi_summaries_rounded(signal)
    scaled = gi_profile.gi_summaries_rounded(gain * signal)

    for key, base_value in base.items():
        if not np.isfinite(base_value):
            continue
        scaled_value = scaled[key]
        relative = abs(scaled_value - base_value) / max(abs(base_value), 1e-12)
        assert relative < RELATIVE_TOLERANCE, (
            f"{name} gain={gain}: {key} base={base_value!r} "
            f"scaled={scaled_value!r} relative={relative:.3e}"
        )


@pytest.mark.parametrize("name,signal", list(_synthetic_windows().items()))
@pytest.mark.parametrize("gain", GAINS)
def test_gi_summaries_grid_invariant_to_gain(name, signal, gain):
    """gi_summaries_grid(c * x) == gi_summaries_grid(x) for c > 0."""
    taus = np.arange(0.05, 2.0 + 1e-12, 0.01)
    base = gi_profile.gi_summaries_grid(signal, taus)
    scaled = gi_profile.gi_summaries_grid(gain * signal, taus)

    assert base["gigrid_valid_fraction"] > 0.90
    assert scaled["gigrid_valid_fraction"] > 0.90

    for key in ("gigrid_auc", "gigrid_slope", "gigrid_centroid",
               "gigrid_sd", "gigrid_kurtosis", "gigrid_skewness"):
        base_value = base[key]
        if not np.isfinite(base_value):
            continue
        scaled_value = scaled[key]
        relative = abs(scaled_value - base_value) / max(abs(base_value), 1e-12)
        assert relative < RELATIVE_TOLERANCE, (
            f"{name} gain={gain}: {key} base={base_value!r} "
            f"scaled={scaled_value!r} relative={relative:.3e}"
        )


def test_frozen_summaries_are_not_gain_invariant():
    """Sanity check on the premise: the STANDARD implementation is gain
    sensitive. If this ever starts passing, something upstream changed and
    the paper's motivating problem no longer holds as stated."""
    rng = np.random.default_rng(2)
    signal = np.cumsum(rng.standard_normal(1200)) * 1e-2
    base = gi_profile.frozen_summaries(signal)
    scaled = gi_profile.frozen_summaries(4.0 * signal)

    relative = abs(scaled["frozen_total"] - base["frozen_total"]) / max(
        abs(base["frozen_total"]), 1e-12)
    assert relative > 1e-3, (
        "frozen_summaries appears gain-invariant at 4x gain; this "
        "contradicts the paper's stated motivation and should be "
        "investigated before trusting any other result here"
    )
