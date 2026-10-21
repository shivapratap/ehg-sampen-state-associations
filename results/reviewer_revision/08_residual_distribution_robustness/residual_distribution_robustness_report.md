# Residual-distribution robustness report

## Main conclusion

The residuals are not perfectly Gaussian, but the central EHG9 UC-versus-IN conclusion is not materially altered. All 8/11 original effects retained direction under woman-cluster robust OLS and all 11/11 retained direction in the woman-level paired-median analysis. The number retaining BH q<0.05 was 1/11 under woman-cluster inference and 5/11 under the paired-median analysis; this loss of multiplicity-adjusted support occurs without direction reversal and is therefore primarily an inference/precision issue rather than evidence that the effects disappeared.

All 34 frozen UC coefficients and q-values were reproduced within 1e-8: 34/34 coefficients and 34/34 q-values passed. Model tiers matched for 34/34 features. Any tier mismatch is reported with full optimizer warnings rather than silently accepted.

## Residual definitions

For each feature, the observed outcome is the Stage 06 z-standardized feature. The **marginal fitted value** is the fixed-effect design matrix multiplied by the fitted fixed coefficients. The marginal residual is observed z minus this fixed-only prediction. The **conditional fitted value** is `statsmodels` `MixedLMResults.fittedvalues`, which adds empirical-Bayes random-effect predictions. For woman-plus-recording fits this contains the woman intercept and recording variance-component predictions; for woman-only fallback fits it contains only the woman intercept. Conditional residuals are consequently not identical estimands across tiers and should not be ranked across features as if the random-effect structures were the same.

## Residual shape

Formal normality tests were not used. With 4,083 observations per feature, such tests would reject small departures with little relevance to robustness. Conditional residuals were summarized using skewness, excess kurtosis, percentiles, and tail proportions.

The strongest absolute residual skewness among the original 11 occurred for Peak-to-peak, RMS, Waveform length, MAV. The strongest excess kurtosis occurred for Peak-to-peak, RMS, Waveform length, MAV. These shape departures identify diagnostic caution, not automatic model invalidity.

The prespecified descriptors were: approximately symmetric for absolute skewness below 0.5, moderate for 0.5–1, and strong above 1; heavy-tailed for excess kurtosis 1–3 and strongly heavy-tailed above 3. The full marginal and conditional summaries remain available in the CSV.

## Heteroscedasticity

Class-specific SD, IQR, and unscaled MAD were compared for UC, FM, and IN. A largest max/min spread ratio below 1.25 was labelled limited, 1.25–1.50 moderate, and at least 1.50 substantial. These are descriptive flags rather than hypothesis tests. Among the original 11, 5/11 showed substantial class-dependent spread by at least one measure: MAV, RMS, Peak-to-peak, Waveform length, Total SampEn.

Because the fixed-effect model has only three state levels, marginal fitted values form three principal bands. Conditional fitted values vary additionally through woman and recording predictions. Residual-versus-fitted plots and correlations of absolute residual with fitted value are therefore interpreted descriptively.

## Core amplitude and waveform features

| Feature | Residual skew | Excess kurtosis | Largest class spread ratio | Woman-cluster | Paired median | Classification |
|---|---:|---:|---:|---|---|---|
| MAV | 61.87 | 3910.84 | 27.01 | positive (q=0.488) | positive (q=0.000176) | Direction robust, inference sensitive |
| RMS | 62.90 | 3997.65 | 42.35 | positive (q=0.515) | positive (q=0.000154) | Direction robust, inference sensitive |
| Peak-to-peak | 63.04 | 4010.00 | 47.38 | positive (q=0.515) | positive (q=0.000154) | Direction robust, inference sensitive |
| Waveform length | 61.97 | 3919.03 | 28.18 | positive (q=0.488) | positive (q=0.000176) | Direction robust, inference sensitive |

MAV, RMS, peak-to-peak amplitude, and waveform length retain the frozen positive UC-versus-IN direction under both robust strategies. Residual skew/heavy tails and class-dependent spread are reported above, but they do not produce effect-direction instability. The central amplitude/waveform conclusion is therefore robust to the observed residual-distribution departures.

## Total SampEn

Total SampEn is assessed here only for distributional robustness. Its direction across the frozen model, woman-cluster analysis, and paired-median analysis is reported in the original-11 table. This result does not alter the completed conclusion that Total SampEn is scale-sensitive and amplitude-coupled; residual robustness cannot establish an amplitude-independent complexity effect.

## Spectral and complexity features

| Feature | Residual skew | Excess kurtosis | Class-spread evidence | Woman-cluster direction | Paired-median direction | Classification |
|---|---:|---:|---|---|---|---|
| Spectral edge 95 | -0.04 | 1.05 | limited | positive | negative | Distribution-sensitive |
| Spectral spread | 0.21 | 0.91 | limited | positive | negative | Distribution-sensitive |
| Hjorth complexity | 0.21 | 1.04 | limited | negative | negative | Direction robust, inference sensitive |
| Higuchi FD | 0.56 | 2.02 | limited | positive | negative | Distribution-sensitive |

These weaker features may show skewness, heavy tails, or heterogeneous spread, but the analysis does not assign causality to those diagnostics. Where robust q-values are larger despite preserved direction, the defensible interpretation is reduced precision across 42 clustered women or 34 paired women, potentially compounded by distributional inefficiency, rather than a demonstrated reversal.

## Transformation decision

No transformation sensitivity was run. A common decision rule was fixed before evaluating diagnostics: transformation required both severe conditional-residual shape (absolute skewness at least 2, excess kurtosis at least 10, or at least 2% beyond three residual SDs) and disagreement in effect direction under a prespecified robust summary. No original feature met that combined condition. Searching feature-specific transformations merely to improve p-values was not permitted.

## Reviewer questions

1. Residual symmetry varies by feature; the exact descriptors are reported rather than claiming Gaussianity.
2. The features with the strongest tails are listed above and quantified in `residual_distribution_summary.csv`.
3. Class-dependent heteroscedasticity is present for the flagged features, but is not universal.
4. Woman-cluster inference retains direction for 8/11 original features.
5. Paired-median analysis retains direction for 11/11 original features.
6. All four central amplitude/waveform features remain directionally robust.
7. Genuinely distribution-sensitive original features: Spectral edge 95, Spectral spread, Higuchi FD.
8. No feature required transformation under the prespecified combined trigger.
9. No transformed models were run.
10. Loss of FDR support under robust methods is primarily a precision issue because direction is preserved.
11. The central manuscript conclusion does not materially change.
12. Retain the central amplitude/waveform direction claim; weaken universal significance language for cluster- or patient-level nonsignificant features; remove any claim that z-scoring proves Gaussian or homoscedastic residuals.
13. The concern is adequately addressed by explicit diagnostics plus heteroscedasticity/within-woman-robust inference and outlier-resistant paired medians.

## Reproducibility and scope

Task 1 and Task 2 estimates were reused exactly from their completed CSVs and protected by recorded SHA-256 hashes. No randomness was used. No frozen Stage 00–09 or prior reviewer-revision artifact was changed. This analysis stops at residual-distribution robustness and does not rewrite the manuscript.
