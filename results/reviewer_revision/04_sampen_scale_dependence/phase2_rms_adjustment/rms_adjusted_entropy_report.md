# RMS-adjusted entropy analysis with corrected dependence

## Bottom line

After adjusting for same-row log RMS and explicitly modelling shared physical windows and women, Total SampEn does not retain a stable positive independent UC-versus-IN association. The corrected mixed model reverses to a small negative coefficient, whereas the woman-cluster sensitivity is near zero and positive. Their direction and inferential support disagree, so there is no robust adjusted state effect. In this implementation, Total SampEn is best described as an amplitude-coupled, scale-sensitive entropy-derived summary rather than evidence of complexity independent of amplitude.

## Population and locked model

The corrected analysis used UC and IN only: 3691 physical windows (333 UC, 3358 IN), 14764 channel rows, 42 women, 112 recordings, and four complete channels per window. Each outcome was standardized over the pooled UC/IN EHG9-EHG12 rows with sample SD (`ddof=1`). RMS was transformed by the natural logarithm using the same row/channel's FWH RMS and then standardized over the identical population.

Primary corrected models used ML MixedLM with a woman random intercept and a physical-window variance component, exactly following the cross-channel model recommendation. Required sensitivity models used the identical fixed effects with CR1 woman-cluster-robust covariance, degrees-of-freedom correction, and t inference (42 clusters). The fixed-effect formulas were:

- Unadjusted: `z_entropy ~ state + channel`
- Adjusted: `z_entropy ~ state + channel + z_log_RMS`

## Corrected UC-versus-IN results

| Feature | Window-aware unadjusted | Window-aware RMS-adjusted | Absolute coefficient change | Woman-cluster RMS-adjusted |
|:--|:--|:--|--:|:--|
| Total SampEn | 0.299 (0.216, 0.383); p=2.26e-12 | -0.086 (-0.166, -0.005); p=0.0363 | 0.385 | 0.030 (-0.075, 0.134); p=0.57 |
| fixed-r SampEn | -0.213 (-0.307, -0.119); p=8.35e-06 | 0.235 (0.150, 0.321); p=7.38e-08 | 0.448 | 0.083 (-0.129, 0.295); p=0.435 |
| SVD entropy | -0.053 (-0.135, 0.030); p=0.209 | 0.020 (-0.061, 0.101); p=0.63 | 0.073 | 0.256 (0.027, 0.485); p=0.0294 |

P-values are feature-specific and unadjusted; no multiplicity family was added because the request predefined exactly these three diagnostic outcomes and did not specify FDR inference.

## Reproduction and explanation of the old 0.938 result

The older Table 8 code used all IN, FM, and UC observations across all four channels (16332 rows), standardized entropy outcomes over that pool, transformed RMS with the natural logarithm and standardized it over the same pool, and fitted REML models with woman random intercepts plus recording variance components. It did **not** model physical-window dependence. Its reported residual was the conditional residual, `observed - (fixed + random prediction)`.

The old values reproduced: Total SampEn unadjusted beta 0.311097, adjusted beta -0.015530, and conditional-residual Pearson correlation with **raw RMS** 0.937639 (frozen 0.937639). The same conditional residual had Pearson correlation -1.77e-16 with the actually modelled log-RMS covariate, while its Spearman correlation was -0.682 with either raw or log RMS.

Therefore, the 0.938 number is **correctly reproducible numerically but based on a mismatched RMS transformation**: the model adjusted linearly for standardized log RMS, then the diagnostic switched back to raw RMS for Pearson correlation. Raw RMS was extremely right-skewed (median 0.0105, 99.9th percentile 0.0787, maximum 15.583); only 5 rows exceeded RMS=1. Removing those rows changes the raw-RMS Pearson correlation to -0.435. Thus a few extreme raw-scale observations create most of the positive 0.938 Pearson coefficient. It should not have been interpreted as proof that the fitted log-RMS adjustment failed. The non-zero Spearman result still shows a strong nonlinear/rank residual pattern, so the corrected conclusion should not claim complete amplitude independence from the near-zero Pearson log-RMS value alone.

## Feature-specific evidence classification

### Total SampEn

- State association: model-dependent after adjustment.
- Direct amplitude association: direct log-RMS association supported by both corrected analyses.
- Mathematical scale sensitivity: demonstrated in Phase 1; the fixed signal-unit profile grid and unnormalised sum make it amplitude-dependent by construction.
- Adjustment: sign reversal after RMS adjustment.
- Residual association: conditional pooled log-RMS Pearson r=-0.000; Spearman rho=-0.698.

Total SampEn is strongly associated with RMS statistically and is scale-sensitive mathematically; this observational model cannot determine how much of that amplitude is physiological rather than measurement-related. Because the adjusted estimates are neither directionally concordant nor jointly supported across covariance strategies, “amplitude-coupled” is the defensible interpretation; “independent complexity” is not.

### fixed-r SampEn

- State association: model-dependent after adjustment.
- Direct amplitude association: direct log-RMS association supported by both corrected analyses.
- Mathematical scale sensitivity: Phase 1 found exact numerical invariance under positive scaling because its tolerance scales with signal SD.
- Adjustment: sign reversal after RMS adjustment.
- Residual association: conditional pooled log-RMS Pearson r=-0.000; Spearman rho=0.023.

### SVD entropy

- State association: model-dependent after adjustment.
- Direct amplitude association: direct log-RMS association is model-dependent.
- Mathematical scale sensitivity: the implementation normalizes singular values by their sum, so positive amplitude scaling cancels algebraically; SVD entropy was not part of the Phase-1 rescaling experiment.
- Adjustment: sign reversal after RMS adjustment.
- Residual association: conditional pooled log-RMS Pearson r=0.000; Spearman rho=-0.070.

## Residual diagnostics

`entropy_rms_residual_diagnostics.csv` contains both required residual definitions, each correlated with raw RMS, log RMS, and standardized log RMS by Pearson and Spearman, pooled and separately for EHG9-EHG12. Marginal residuals exclude all random-effect predictions. Conditional residuals subtract both fixed and estimated woman/window random effects. Pearson correlation with log RMS tests remaining linear association on the modelled scale; Spearman correlation tests monotone residual association and is identical for raw and log RMS because the logarithm is monotone.

## Convergence and limitations

All six corrected preferred fits (three features by two variants) produced usable selected optimizer results: True. Boundary/singularity warnings and all attempted optimizers are retained in `model_convergence_diagnostics.csv`. A random physical-window intercept models exchangeable within-window correlation rather than an arbitrary 4x4 channel covariance; the woman-cluster sensitivity protects inference against broader within-woman dependence. Neither model makes the amplitude coefficient causal.

## Scope

No other entropy feature, FM inferential model, alternative scale correction, or manuscript revision was performed.

## Inputs

- Corrected event-unit tables: `<repo>/results/paper1/06_mixed_models/event_unit_analysis_60s_EHG9.csv` through `EHG12.csv`
- Old analysis source: `<repo>/src/additional_analysis/step06_amplitude_adjusted_model.py`
- Frozen old result: `<repo>/results/paper1/additional_analysis/phase3/step06_amplitude_adjusted_results.csv`
- Phase-1 scale audit: `<repo>/reviewer_revision/04_sampen_scale_dependence/phase1_scale_audit`
- Corrected dependence recommendation: `<repo>/reviewer_revision/03_cross_channel_dependence/phase1_model_audit/cross_channel_model_recommendation.md`
