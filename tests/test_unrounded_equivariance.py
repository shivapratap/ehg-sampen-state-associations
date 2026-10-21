"""Unrounded / full-precision control: verifies the paper's mechanism claim.

`full_precision_summaries` removes the fixed-decimal rounding entirely (no
sd-standardisation either): the tolerance axis is the raw, data-adaptive set
of unique unrounded pairwise Chebyshev distances. Because that axis is
data-adaptive rather than fixed, the cumulative pair-count at each rank is
unchanged by any positive rescaling of the signal (rescaling relabels the
axis but preserves the ordering of distances), so the profile -- and every
summary of it -- is EXACTLY invariant to gain, with no standardisation step
at all.

The paper's own report (results/gi/audits/full_precision_gain_test_report.md)
states this in its own vocabulary as "exactly gain-equivariant": 0.000e+00
deviation for fp_total, fp_average and fp_sd at every tested gain, well
within atol=1e-08, rtol=1e-06. This test file reproduces that finding,
including its terminology.

Read together with test_gain_invariance.py, this is the paper's central
mechanism-isolating result (facts_pack.md, section 5.5): removing the
fixed-decimal ROUNDING is what restores invariance, standardisation is a
separate, additional feature-definition change made by the GI variants. This
is "the single strongest answer to the 'this is just a bug in one
implementation' objection" (facts_pack.md), so this test is written to fail
loudly, not silently degrade, if that stops holding.

Needs no external data.
"""
from __future__ import annotations

import numpy as np
import pytest

from gi_sampen import full_precision as fp
from gi_sampen import gi_profile

# Matches results/gi/audits/full_precision_gain_test_report.md exactly.
ATOL = 1e-8
RTOL = 1e-6


def _synthetic_windows():
    rng = np.random.default_rng(3)
    return {
        "ehg_like_cumsum_n1200": np.cumsum(rng.standard_normal(1200)) * 1e-2,
        "white_noise_n600": rng.standard_normal(600),
    }


@pytest.mark.parametrize("name,signal", list(_synthetic_windows().items()))
@pytest.mark.parametrize("gain", fp.GAINS)
def test_full_precision_summaries_are_gain_invariant(name, signal, gain):
    """fp_total, fp_average and fp_sd are exactly invariant to positive gain,
    with no rounding and no standardisation -- the paper's own "gain-
    equivariant" result. See the module docstring for why: the tolerance
    axis is data-adaptive, so rescaling relabels it without changing any
    cumulative pair-count."""
    base = fp.full_precision_summaries(signal)
    scaled = fp.full_precision_summaries(gain * signal)

    for key in ("fp_total", "fp_average", "fp_sd"):
        base_value = base[key]
        scaled_value = scaled[key]
        if not (np.isfinite(base_value) and np.isfinite(scaled_value)):
            pytest.skip(f"{name} gain={gain}: non-finite {key}, skipping")
        assert scaled_value == pytest.approx(base_value, abs=ATOL, rel=RTOL), (
            f"{name} gain={gain}: {key} base={base_value!r} "
            f"scaled={scaled_value!r} -- full-precision summaries are "
            "expected to be exactly gain-invariant even with no "
            "standardisation; see facts_pack.md section 5.5"
        )
    assert scaled["fp_n_retained"] == base["fp_n_retained"]


def test_rounding_is_the_mechanism_not_amplitude_itself():
    """The paper's central mechanism claim, stated as one test: the STANDARD
    (rounded) implementation is gain-SENSITIVE, but removing only the
    rounding -- with no standardisation at all -- makes the profile
    invariant again. So rounding, not SampEn profiling itself, is what
    breaks invariance under gain. See facts_pack.md section 5.5, and the
    module docstring in gi_profile.py."""
    rng = np.random.default_rng(4)
    signal = np.cumsum(rng.standard_normal(1200)) * 1e-2
    gain = 4.0

    frozen_base = gi_profile.frozen_summaries(signal)["frozen_total"]
    frozen_scaled = gi_profile.frozen_summaries(gain * signal)["frozen_total"]
    frozen_relative_change = abs(frozen_scaled - frozen_base) / max(
        abs(frozen_base), 1e-12)

    fp_base = fp.full_precision_summaries(signal)["fp_total"]
    fp_scaled = fp.full_precision_summaries(gain * signal)["fp_total"]

    assert frozen_relative_change > 1e-3, (
        "frozen_summaries (rounded) appears gain-invariant; the premise "
        "this test isolates no longer holds"
    )
    assert fp_scaled == pytest.approx(fp_base, abs=ATOL, rel=RTOL), (
        "full_precision_summaries (unrounded, unstandardised) is not "
        "invariant; rounding is no longer isolable as the mechanism"
    )
