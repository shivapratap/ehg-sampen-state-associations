# Methods

## The problem

The standard finite-resolution SampEn-profile implementation quantizes each
pairwise Chebyshev distance to a fixed absolute resolution -- three decimals,
i.e. 0.001 signal units -- before forming the tolerance support: the profile
is evaluated at the *distinct occupied quantized values*, not a sweep of a
regular grid. That quantizer does not commute with a multiplicative gain.
Under gain `c > 0`, every pairwise distance scales by exactly `c`, but
distances that shared a quantized level at unit gain can separate under
scaling, and distinct levels can merge -- so the profile is evaluated at a
support that is not the scaled image of the original one. Total SampEn is the
most exposed summary, because the support's cardinality grows roughly in
proportion to signal amplitude divided by the fixed resolution.

This is a property of the *released implementation*
(`Sample entropy/CHM.m`, commit `2a1d496c8bcf1e859ae9887d9911d168dd9a17dc`,
2020-09-24), not of the published SampEn-profile formulation, which contains
no rounding step. See `docs/provenance.md` for the full attribution chain.

## The method: Gain-Invariant (GI) SampEn profiling

Divide each analysis window by its own standard deviation before an
otherwise unchanged profiling stage:

```
z = x / sd(x)
```

then apply the standard construction -- same Chebyshev distances, same
three-decimal rounding, same retention rule, same six summaries -- to `z`
instead of `x`. This is `gi_summaries_rounded` in
`src/gi_sampen/gi_profile.py`. For `c > 0`, `sd(c*x) = c*sd(x)`, so `z` is
unchanged and the construction is invariant to gain up to floating point.
`tests/test_gain_invariance.py` asserts this at six gains (0.25, 0.5, 1.3, 2,
4, 7.1) to < 1e-8 relative tolerance -- a correctness test on the
implementation, not a scientific result.

**Standardisation also refines the tolerance resolution as a side effect.**
EHG windows in this dataset have `sd` of order 1e-2 signal units, so a 0.001
signal-unit grid gives roughly 1e2 bins, whereas a 0.001-sd grid gives
roughly 1e3-1e4. The GI features are therefore not merely "the standard
features made invariant" -- they are computed on a finer profile. This is a
feature-definition change, not a bug fix, and the manuscript states it as
such (see the design notes

**Gain-invariance is not cross-window comparability.** The GI *rounded*
variant still runs out to each window's own maximum distance, so its support
size varies between windows even though it is invariant to gain within one
window. `gi_summaries_grid` addresses this: it evaluates the standardised
profile on a *fixed dimensionless tolerance grid*, shared by every window,
with no rounding at all. GI-AUC (the mean profile height over the grid),
GI-slope and GI-centroid are defined on this grid and are the proposed
cross-window-comparable summaries.

## The mechanism-isolating control

A third variant, `full_precision_summaries`
(`src/gi_sampen/full_precision.py`), removes the fixed-decimal rounding
entirely but applies **no standardisation**. Because its tolerance axis is
the raw, data-adaptive set of unique unrounded pairwise distances, the
cumulative pair-count at each rank is unchanged by any positive rescaling of
the signal -- rescaling relabels the axis but preserves the ordering of
distances. The result is *exactly* gain-invariant (0.000e+00 deviation; see
`results/gi/audits/full_precision_gain_test_report.md`), with no
standardisation step at all. This isolates fixed-decimal quantization,
rather than SampEn profiling itself, as the mechanism responsible for the
standard implementation's gain sensitivity -- the strongest available answer
to "this is just amplitude sensitivity, not a quantization artifact."

## Comparators

- **Fixed-r SampEn** (Richman-Moorman, `r = 0.2 * sd`): already exactly
  gain-invariant by construction, carried through as the reference point
  against which the profile representations must justify themselves.
- **The standard (rounded) implementation itself**
  (`frozen_summaries`): the parity target, and the thing whose gain
  sensitivity motivates this work. Not proposed as a feature.

## Scope of the invariance claim

Constant **positive multiplicative** gain only. No claim is made about
additive offsets, non-constant (time-varying) gain, or negative scaling.

## Terminology used throughout this repository

| Term | Meaning |
|---|---|
| Standard finite-resolution implementation | The released three-decimal SampEn-profile implementation (`frozen_summaries` here) |
| Gain-Invariant (GI) SampEn profiling | The proposed SD-standardised variant |
| Unrounded / full-precision control | The mechanism-isolation variant: rounding removed, no standardisation |

"Frozen" is avoided in every public-facing README, filename, figure label,
table heading and prose document in this repository. It survives only inside
committed CSV column names inherited from the prior study
(`frozen_estimate`, `frozen_total`, etc.), where renaming would break the
join back to that study's own cached files; `docs/provenance.md` documents
what they mean.

## What is not claimed

- No claim of priority for gain-invariant / amplitude-robust entropy
  estimation in general.
- No claim that SampEn profiling is inherently gain-sensitive -- only that
  the specific released finite-resolution implementation, as characterised
  above, is.
- No claim that GI-SampEn removes all amplitude effects.
