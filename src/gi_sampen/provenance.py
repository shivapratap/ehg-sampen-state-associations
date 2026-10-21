#!/usr/bin/env python3
"""Machine-readable provenance for the implementations in this repository.

The prose version is ``docs/provenance.md``. This module exists so that the
same identifiers can be asserted in tests and stamped into generated outputs,
rather than living only in documentation that can drift.

Nothing here changes any numerical result.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from typing import Optional

__all__ = [
    "Implementation",
    "REFERENCE_MATLAB",
    "FEATURE_ENGINE",
    "UPSTREAM_PACKAGE",
    "THIS_REPOSITORY",
    "PROFILE_DISTANCE_DECIMALS",
    "chain",
    "as_dict",
]


@dataclass(frozen=True)
class Implementation:
    """One link in the standard finite-resolution implementation chain."""

    name: str
    url: Optional[str]
    commit: Optional[str]
    version: Optional[str]
    language: str
    role: str
    rounding: str
    notes: str


#: Number of decimals to which pairwise Chebyshev distances are rounded before
#: the tolerance support is formed. This is a property of the *released
#: implementation*, not of the published sample-entropy-profile formulation,
#: which contains no rounding step. See docs/provenance.md.
PROFILE_DISTANCE_DECIMALS = 3


REFERENCE_MATLAB = Implementation(
    name="Sample entropy profile, original MATLAB release",
    url=None,
    commit="2a1d496c8bcf1e859ae9887d9911d168dd9a17dc",
    version=None,
    language="MATLAB",
    role=(
        "Canonical source of the standard finite-resolution behaviour. The "
        "rounding enters at `round(d,3)` / `round(d1,3)` in "
        "`Sample entropy/CHM.m`, introduced by this commit (2020-09-24)."
    ),
    rounding="three-decimal, absolute, in signal units",
    notes=(
        "The same release's approximate-entropy profile "
        "(`Approximate entropy/CHMforApEn.m`) takes unique() over the raw "
        "distances with no rounding, so the quantizer is a deliberate, "
        "method-specific choice rather than an incidental one."
    ),
)

FEATURE_ENGINE = Implementation(
    name="amrita-biosignal-feature-engine",
    url="https://github.com/shivapratap/amrita-biosignal-feature-engine",
    commit="4dc5025b265f300d7e3375118f7794bfb5fca711",
    version="0.2.0.dev0",
    language="Python",
    role=(
        "Produced the cached finite-resolution feature table "
        "(`fwh_entropy_*sampen` columns) that the prior study analysed and "
        "that this repository uses as its parity target."
    ),
    rounding="three-decimal, absolute, in signal units",
    notes=(
        "Not a runtime dependency here. The relevant module is "
        "`entropy/sample_entropy_profile.py`."
    ),
)

UPSTREAM_PACKAGE = Implementation(
    name="sampen-profile",
    url="https://github.com/shivapratap/sampen-profile",
    commit="263b72f7b906208241291c0579e8ec103fb52d79",
    version="0.1.0",
    language="Python",
    role=(
        "Released, installable implementation of the same method. Published "
        "2026-07-20, i.e. AFTER the analysis in this repository was run, so "
        "it did not produce any manuscript number. It is pinned here as a "
        "test-only dependency: `tests/test_standard_parity.py` asserts that "
        "`gi_sampen.gi_profile.frozen_summaries` agrees with "
        "`sampen_profile.sample_entropy_profile` on deterministic signals."
    ),
    rounding="three-decimal, absolute, in signal units",
    notes=(
        "Relevant functions: `sampen_profile.reference.sample_entropy_profile` "
        "(readable backend) and `sampen_profile.optimized.sample_entropy_profile`. "
        "Both round pairwise Chebyshev distances to three decimals and then "
        "take the distinct occupied values as the tolerance support."
    ),
)

THIS_REPOSITORY = Implementation(
    name="gi_sampen.gi_profile",
    url=None,
    commit=None,
    version=None,
    language="Python",
    role=(
        "Implementation of record for every manuscript number. "
        "`frozen_summaries` reimplements the standard finite-resolution "
        "construction; `gi_summaries_rounded` is gain-invariant (GI) SampEn "
        "profiling, which applies the otherwise unchanged construction to "
        "z = x / sd(x); `gi_summaries_grid` evaluates the profile on a fixed "
        "dimensionless grid with no rounding."
    ),
    rounding=(
        "frozen_summaries: three decimals in signal units. "
        "gi_summaries_rounded: three decimals in units of the window's own "
        "standard deviation. gi_summaries_grid: none. "
        "full_precision: none, and no standardisation either."
    ),
    notes=(
        "Deliberately standalone rather than a call into the feature engine, "
        "so that parity against the cached values can be asserted numerically "
        "before anything new is interpreted."
    ),
)


def chain() -> tuple[Implementation, ...]:
    """The implementations in provenance order."""
    return (REFERENCE_MATLAB, FEATURE_ENGINE, UPSTREAM_PACKAGE, THIS_REPOSITORY)


def as_dict() -> dict:
    """Serialisable form, for stamping into generated outputs."""
    return {
        "profile_distance_decimals": PROFILE_DISTANCE_DECIMALS,
        "implementations": [asdict(item) for item in chain()],
    }


if __name__ == "__main__":
    import json

    print(json.dumps(as_dict(), indent=2))
