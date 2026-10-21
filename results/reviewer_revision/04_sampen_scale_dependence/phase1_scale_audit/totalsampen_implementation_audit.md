# Total SampEn implementation audit

## Authoritative production path

The frozen Stage 04 EHG9/FWH Total SampEn value is produced by `sample_entropy_profile_summaries` in `<repo>/src/04_feature_extraction.py`. Stage 04 calls this local, memory-safe implementation after full-recording FWH filtering and window slicing/decimation. The vendored engine supplies RMS and other registered features, but not these six profile summaries in the production path.

## Exact construction

- Input: finite 1-D float64 EHG9 FWH window, 1,200 samples (60 s at 20 Hz).
- Embedding order: m=2. Both m and m+1 templates use the aligned count L=N-m=1,198.
- Distance: Chebyshev (maximum coordinate-wise absolute difference).
- Grid: signal-unit tolerance spacing 0.001, because distances are rounded to 3 decimals. The grid starts at 0 and ends at the maximum pairwise distance after three-decimal rounding.
- Pair counting: SciPy `cKDTree.count_neighbors`, with ordered self-pairs removed from the zero-distance shell. The half-even three-decimal rounding boundaries are reproduced by alternating `nextafter` handling at half steps.
- Probabilities: cumulative ordered non-self pair counts divided by L(L-1), separately for m (`b`) and m+1 (`a`).
- Profile value at a retained tolerance: `log((b + 1e-12)/(a + 1e-12))`.
- Retention: the shell must be observed in m or m+1, the profile must be finite, and both cumulative probabilities must exceed 1e-12. Thus the profile is not evaluated as an equally sized fixed normalized grid; the number of retained values may change under amplitude scaling.
- Total SampEn: arithmetic sum of all retained profile values. Average/median/SD (population SD, ddof=0), SciPy Fisher excess kurtosis with bias correction disabled by default (`bias=True`), and SciPy skewness (`bias=True`) are computed over the same retained values.
- Fixed-r SampEn comparator: tolerance r=0.2 times the window population SD, aligned m/m+1 Chebyshev templates, ordered non-self matches, and `-log(A/B)`.
- RMS comparator: square root of the mean squared samples.

## Scale-dependence hypothesis fixed before examining results

For positive scaling x->c*x, all exact distances scale by c. Fixed-r SampEn should therefore be invariant because its SD-derived tolerance also scales by c; RMS should be exactly proportional to c. Total SampEn is not guaranteed to be invariant because it sums values over a fixed 0.001 signal-unit rounded-distance grid and only retains occupied shells. Scaling changes the maximum rounded distance, possible grid size, shell occupancy, and usually the number of retained profile values. Average and other profile summaries can also change through quantization and the changing retained set, although they do not directly accumulate more terms as Total does.

Numerical invariance was prespecified as `abs(value - 1x) <= 1e-10 + 1e-8*abs(1x)`. Percent change and ratios are reported only when `abs(1x)>1e-10`.
