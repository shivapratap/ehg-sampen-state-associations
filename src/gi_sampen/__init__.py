"""Gain-invariant sample entropy profiling for electrohysterography.

The numerical core is :mod:`gi_sampen.gi_profile`:

``frozen_summaries``
    The standard finite-resolution construction -- pairwise Chebyshev
    distances rounded to three decimals in signal units, the distinct occupied
    values taken as the tolerance support.

``gi_summaries_rounded``
    Gain-invariant (GI) SampEn profiling. Identical to the above except that
    the window is divided by its own standard deviation first, so the
    tolerance support is expressed in dimensionless units.

``gi_summaries_grid``
    The profile of the standardised window on a fixed dimensionless tolerance
    grid, with no rounding, which makes profiles comparable across windows.

``gi_sampen.full_precision``
    The unrounded / full-precision control: rounding removed, no
    standardisation. This isolates the mechanism and is not a proposed
    feature.

:mod:`gi_sampen.provenance` records where the standard implementation comes
from. See ``docs/provenance.md``.
"""

__all__ = ["__version__"]

__version__ = "1.0.0"
