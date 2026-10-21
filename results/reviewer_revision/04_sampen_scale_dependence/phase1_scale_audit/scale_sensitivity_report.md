# Phase 1: Total SampEn amplitude-scale audit

## Executive conclusion

Total SampEn is **not amplitude-scale invariant** in the frozen Stage 04 implementation. Multiplying an otherwise identical EHG9/FWH waveform changes the fixed 0.001-unit distance grid and the retained profile length; because Total SampEn is a sum over those retained values, its value changes even though the underlying waveform morphology is unchanged. This is an implementation/definition property, not biological evidence.

## Audit design and reproduction gate

The audit used 60 deterministic real 60-s EHG9/FWH windows: 20 UC, 20 IN, and 20 FM; 25 women and 32 recordings. Selection was metadata-based and round-robin by within-woman rank, never feature-based. No randomness was used.

All 480 feature-window checks reproduced Stage 04 within the prespecified tolerance (maximum absolute difference 1.78e-15; maximum relative difference 4.53e-14). The scale experiment therefore proceeded.

## Answers to the prespecified questions

1. **Is Total SampEn amplitude-scale invariant?** No. Across the 240 non-reference scaled observations, only 0.0% met the strict numerical-invariance criterion.
2. **Magnitude of change.** Overall median relative changes versus 1x were 0.25x: -75.44%, 0.5x: -50.39%, 1x: 0.00%, 2x: 102.27%, 4x: 309.18%. Full medians, IQRs, ranges, and invariant counts are in `scale_sensitivity_summary.csv`.
3. **Monotonicity.** Of 60 Total-SampEn trajectories, 60 were nondecreasing, 0 nonincreasing, 0 nonmonotonic, and 0 effectively invariant across the five ordered scales.
4. **Does Average SampEn change?** Yes, but much less than Total SampEn: median relative changes were 0.25x -5.56%, 0.5x -1.71%, 2x 2.12%, 4x 4.62%. It avoids direct accumulation of a scale-dependent number of terms, but rounding and retained-set changes remain.
5. **Does Median SampEn change?** Yes: median relative changes were 0.25x -12.08%, 0.5x -4.19%, 2x 4.18%, 4x 8.77%. It is less scale-sensitive than Total SampEn but not strictly invariant under this quantized profile construction.
6. **Does fixed-r SampEn change?** No within the prespecified tolerance: it was invariant in 100.0% of non-reference rows. Its tolerance scales with window SD, as its mathematical definition predicts.
7. **Do retained profile values change?** Yes. Every window changed its retained-value count at one or more non-reference scales; the rounded maximum distance and possible grid-bin count also changed in every window.
8. **Does the changing support explain Total SampEn sensitivity?** Largely. The within-window association between Total SampEn and retained count had median Spearman rho 1.000; the pooled descriptive change-score Spearman rho was 0.966. This strongly supports the prespecified grid-expansion/summation mechanism, while not claiming that term count alone is an exact algebraic decomposition.
9. **Is Total SampEn mathematically scale-sensitive by construction?** Yes. Fixed signal-unit rounding changes the grid/support under x->c*x, and Total SampEn is the unnormalised sum over that changing retained support.
10. **Which summary is most scale-invariant?** `fixed_r_sampen` had the highest strict invariance proportion (100.0%) among the eight audited metrics. RMS is not intended to be invariant and instead showed correct scale proportionality in 100.0% of rows.
11. **Are results consistent across classes?** Yes. All three prespecified classes showed the same implementation-level scale phenomenon; class-specific Total-SampEn median changes were:
- UC: 0.25x -75.42%, 0.5x -50.52%, 1x 0.00%, 2x 102.26%, 4x 308.36%
- IN: 0.25x -75.51%, 0.5x -50.31%, 1x 0.00%, 2x 102.21%, 4x 309.03%
- FM: 0.25x -75.43%, 0.5x -50.45%, 1x 0.00%, 2x 102.71%, 4x 310.09%
12. **Coding bug or expected property?** No coding failure was found. The diagnostic wrapper exactly reconstructed production profile summaries, and all 480 selected-window feature values reproduced Stage 04. The scale behavior is an expected mathematical property of this estimator definition. Describing this Total SampEn implementation as amplitude-invariant would nevertheless be incorrect.

## Mechanistic interpretation

Positive scaling leaves waveform ordering and morphology unchanged, but converts a distance d to c*d. Stage 04 quantizes those scaled distances in absolute signal units at 0.001 resolution. This changes the maximum grid endpoint and which rounded distance shells are occupied. Total SampEn then sums the cumulative log-ratio profile over that changing retained support. The number of summands is therefore amplitude-dependent. Average SampEn divides by retained count and is less directly exposed, while fixed-r SampEn scales its tolerance with window SD and supplies the appropriate invariant control.

## Scope

No RMS-adjusted inferential model, alternative window length, cross-channel analysis, FM sensitivity model, or manuscript rewriting was performed. The result does support undertaking a separate, prespecified amplitude-adjusted sensitivity analysis next.
