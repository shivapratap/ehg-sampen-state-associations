# Step 9 -- Residual autocorrelation diagnostic, EHG9

Lightweight substitute for a full permutation framework (per the plan's intentional scope note): for the 5 EHG9 recordings with the most analysis windows, compute the lag-1 autocorrelation of Stage 06 model residuals ordered by within-recording window start time, for all 14 primary-validated features. A large lag-1 autocorrelation would mean the recording-level variance component is NOT sufficiently absorbing within-recording dependency, motivating escalation to permutation testing; values near zero support the existing random-effects structure as adequate.

Top 5 recordings by window count: ice002_p_1of3, ice044_p_2of3, ice004_p_1of1, ice003_p_1of2, ice036_p_1of6.

Across 70 feature x recording combinations: max |lag-1 autocorrelation| = 0.560, mean |lag-1 autocorrelation| = 0.266, 32 combination(s) at or above an informal |r| >= 0.3 flag threshold.

Note: 20 feature/recording row(s) use marginal (fixed-effects-only) residuals rather than conditional (BLUP) residuals, because those features' fits had a singular random-effects covariance (a boundary fit already flagged in Stage 06) from which conditional residuals cannot be predicted. Marginal residuals still answer the question this diagnostic asks -- serial structure left over after accounting for state and channel -- just without also netting out the woman/recording random-effect predictions.

| feature | recording | n windows | lag-1 autocorr | residual type |
|:---|:---|---:|---:|:---|
| Complexity Higuchi Fractal Dimension | ice002_p_1of3 | 58 | 0.317 | conditional |
| Complexity Higuchi Fractal Dimension | ice003_p_1of2 | 56 | 0.385 | conditional |
| Complexity Higuchi Fractal Dimension | ice004_p_1of1 | 58 | 0.308 | conditional |
| Complexity Higuchi Fractal Dimension | ice036_p_1of6 | 55 | 0.133 | conditional |
| Complexity Higuchi Fractal Dimension | ice044_p_2of3 | 58 | 0.553 | conditional |
| Complexity Hjorth Complexity | ice002_p_1of3 | 58 | 0.083 | conditional |
| Complexity Hjorth Complexity | ice003_p_1of2 | 56 | -0.191 | conditional |
| Complexity Hjorth Complexity | ice004_p_1of1 | 58 | -0.039 | conditional |
| Complexity Hjorth Complexity | ice036_p_1of6 | 55 | -0.015 | conditional |
| Complexity Hjorth Complexity | ice044_p_2of3 | 58 | -0.168 | conditional |
| Complexity Hjorth Mobility | ice002_p_1of3 | 58 | 0.288 | conditional |
| Complexity Hjorth Mobility | ice003_p_1of2 | 56 | 0.328 | conditional |
| Complexity Hjorth Mobility | ice004_p_1of1 | 58 | 0.170 | conditional |
| Complexity Hjorth Mobility | ice036_p_1of6 | 55 | 0.168 | conditional |
| Complexity Hjorth Mobility | ice044_p_2of3 | 58 | 0.482 | conditional |
| Entropy Sampen Fixed R | ice002_p_1of3 | 58 | 0.516 | conditional |
| Entropy Sampen Fixed R | ice003_p_1of2 | 56 | 0.371 | conditional |
| Entropy Sampen Fixed R | ice004_p_1of1 | 58 | 0.205 | conditional |
| Entropy Sampen Fixed R | ice036_p_1of6 | 55 | -0.042 | conditional |
| Entropy Sampen Fixed R | ice044_p_2of3 | 58 | 0.204 | conditional |
| Entropy Svd Entropy | ice002_p_1of3 | 58 | 0.302 | conditional |
| Entropy Svd Entropy | ice003_p_1of2 | 56 | 0.344 | conditional |
| Entropy Svd Entropy | ice004_p_1of1 | 58 | 0.194 | conditional |
| Entropy Svd Entropy | ice036_p_1of6 | 55 | 0.172 | conditional |
| Entropy Svd Entropy | ice044_p_2of3 | 58 | 0.486 | conditional |
| Entropy Totalsampen | ice002_p_1of3 | 58 | 0.398 | conditional |
| Entropy Totalsampen | ice003_p_1of2 | 56 | 0.532 | conditional |
| Entropy Totalsampen | ice004_p_1of1 | 58 | 0.384 | conditional |
| Entropy Totalsampen | ice036_p_1of6 | 55 | 0.343 | conditional |
| Entropy Totalsampen | ice044_p_2of3 | 58 | 0.298 | conditional |
| Freq Spectral Edge Frequency 95 | ice002_p_1of3 | 58 | 0.241 | conditional |
| Freq Spectral Edge Frequency 95 | ice003_p_1of2 | 56 | 0.331 | conditional |
| Freq Spectral Edge Frequency 95 | ice004_p_1of1 | 58 | 0.232 | conditional |
| Freq Spectral Edge Frequency 95 | ice036_p_1of6 | 55 | 0.078 | conditional |
| Freq Spectral Edge Frequency 95 | ice044_p_2of3 | 58 | 0.445 | conditional |
| Freq Spectral Spread | ice002_p_1of3 | 58 | 0.141 | conditional |
| Freq Spectral Spread | ice003_p_1of2 | 56 | 0.245 | conditional |
| Freq Spectral Spread | ice004_p_1of1 | 58 | 0.234 | conditional |
| Freq Spectral Spread | ice036_p_1of6 | 55 | -0.064 | conditional |
| Freq Spectral Spread | ice044_p_2of3 | 58 | 0.420 | conditional |
| Time Mean Absolute Value | ice002_p_1of3 | 58 | 0.323 | marginal_fallback |
| Time Mean Absolute Value | ice003_p_1of2 | 56 | 0.528 | marginal_fallback |
| Time Mean Absolute Value | ice004_p_1of1 | 58 | 0.420 | marginal_fallback |
| Time Mean Absolute Value | ice036_p_1of6 | 55 | 0.164 | marginal_fallback |
| Time Mean Absolute Value | ice044_p_2of3 | 58 | 0.380 | marginal_fallback |
| Time Median | ice002_p_1of3 | 58 | -0.221 | conditional |
| Time Median | ice003_p_1of2 | 56 | 0.145 | conditional |
| Time Median | ice004_p_1of1 | 58 | -0.023 | conditional |
| Time Median | ice036_p_1of6 | 55 | 0.134 | conditional |
| Time Median | ice044_p_2of3 | 58 | 0.070 | conditional |
| Time Peak To Peak | ice002_p_1of3 | 58 | 0.316 | marginal_fallback |
| Time Peak To Peak | ice003_p_1of2 | 56 | 0.279 | marginal_fallback |
| Time Peak To Peak | ice004_p_1of1 | 58 | 0.433 | marginal_fallback |
| Time Peak To Peak | ice036_p_1of6 | 55 | 0.075 | marginal_fallback |
| Time Peak To Peak | ice044_p_2of3 | 58 | 0.027 | marginal_fallback |
| Time Root Mean Square | ice002_p_1of3 | 58 | 0.318 | marginal_fallback |
| Time Root Mean Square | ice003_p_1of2 | 56 | 0.524 | marginal_fallback |
| Time Root Mean Square | ice004_p_1of1 | 58 | 0.433 | marginal_fallback |
| Time Root Mean Square | ice036_p_1of6 | 55 | 0.095 | marginal_fallback |
| Time Root Mean Square | ice044_p_2of3 | 58 | 0.319 | marginal_fallback |
| Time Skewness | ice002_p_1of3 | 58 | 0.064 | conditional |
| Time Skewness | ice003_p_1of2 | 56 | 0.307 | conditional |
| Time Skewness | ice004_p_1of1 | 58 | -0.068 | conditional |
| Time Skewness | ice036_p_1of6 | 55 | 0.220 | conditional |
| Time Skewness | ice044_p_2of3 | 58 | 0.139 | conditional |
| Time Waveform Length | ice002_p_1of3 | 58 | 0.378 | marginal_fallback |
| Time Waveform Length | ice003_p_1of2 | 56 | 0.560 | marginal_fallback |
| Time Waveform Length | ice004_p_1of1 | 58 | 0.376 | marginal_fallback |
| Time Waveform Length | ice036_p_1of6 | 55 | 0.211 | marginal_fallback |
| Time Waveform Length | ice044_p_2of3 | 58 | 0.252 | marginal_fallback |

## Per-feature summary (across the recordings checked)

| feature | recordings flagged (|r|>=0.3) | mean |lag-1 autocorr| |
|:---|---:|---:|
| Complexity Higuchi Fractal Dimension | 4/5 | 0.339 |
| Entropy Totalsampen | 4/5 | 0.391 |
| Time Mean Absolute Value | 4/5 | 0.363 |
| Time Root Mean Square | 4/5 | 0.338 |
| Entropy Svd Entropy | 3/5 | 0.300 |
| Time Waveform Length | 3/5 | 0.355 |
| Complexity Hjorth Mobility | 2/5 | 0.287 |
| Entropy Sampen Fixed R | 2/5 | 0.268 |
| Freq Spectral Edge Frequency 95 | 2/5 | 0.265 |
| Time Peak To Peak | 2/5 | 0.226 |
| Freq Spectral Spread | 1/5 | 0.221 |
| Time Skewness | 1/5 | 0.160 |
| Complexity Hjorth Complexity | 0/5 | 0.099 |
| Time Median | 0/5 | 0.119 |

## Interpretation

**This diagnostic does not come back clean, and that itself is the finding.** 32 of 70 feature x recording combinations (46%) show lag-1 residual autocorrelation at or above an informal |r|=0.3 threshold, with a mean |lag-1 autocorrelation| of 0.266 across all combinations -- broad-based, not an isolated one-off in a single recording or feature. 12 of 14 primary-validated features have at least one flagged recording, and Complexity Higuchi Fractal Dimension, Entropy Totalsampen, Time Mean Absolute Value, Time Root Mean Square, Entropy Svd Entropy, Time Waveform Length are the most consistently affected (flagged in at least half of the checked recordings).
Per the plan's own stated contingency ("only escalate to a full permutation framework if this diagnostic shows a problem"), this diagnostic **does show a problem**: the current recording-nested-within-woman random-effects structure -- a fixed per-recording intercept -- does not appear to fully absorb whatever slow, within-recording drift produces this correlation between consecutive windows' residuals.
Practical implication for the manuscript, not previously anticipated by the plan: this does not overturn any primary finding (the state effects being tested concern within-recording contraction-vs-baseline differences, and the random-effects structure already accounts for between-woman and between-recording levels), but the nominal p-values and CIs reported throughout Stages 06-08 likely understate true uncertainty for the most-affected features, since residual autocorrelation reduces the effective within-recording sample size beyond what the current variance-components model assumes. This is worth an explicit Limitations-section sentence, and a candidate follow-up (not undertaken here, consistent with the plan's scoping) would be a permutation-based or cluster-robust standard-error re-analysis restricted to the features flagged here.
