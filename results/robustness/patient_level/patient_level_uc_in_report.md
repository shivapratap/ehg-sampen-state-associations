# Patient-level EHG9 UC-versus-IN sensitivity analysis

## Question

Do EHG9 uterine-contraction (UC) versus inactivity (IN) differences remain visible after reducing the data to one paired UC-versus-IN summary per woman?

## Frozen inputs

- `/Users/shivangayathri/Work/EHGProcessing/ehg-sampen-state-associations/work/primary/confirmatory_feature_registry_60s_fwh.csv` (SHA-256 `bdbd34522dfd5d64fa852fb90ab734cca03eec5cb3aa11771b95b5997450ae12`)
- `/Users/shivangayathri/Work/EHGProcessing/ehg-sampen-state-associations/work/primary/event_unit_analysis_60s_EHG9.csv` (SHA-256 `b9d383031470aa718c63210ae4f0f9ab7d3fd13a5a009e1a405a3d5e16e67d68`)
- `/Users/shivangayathri/Work/EHGProcessing/ehg-sampen-state-associations/work/primary/mixed_model_results_60s_fwh_EHG9.csv` (SHA-256 `dcfc060c1b0311b72eed8f81adc6c0a57312787566368a2449d9562df7ad8638`)
- Paired analysis checks the migrated primary feature family directly.

The script read these inputs without refitting or changing the Stage 06 mixed model. It verified the 34-feature registry, the 34 frozen EHG9 UC-vs-IN rows, their BH q-values, and the frozen 11-feature q<0.05 set.

## Analysis population

- Total women in the frozen EHG9 event-unit table: **42**.
- Women with both UC and IN: **34**.
- Women with IN only: **7**.
- Women with UC only: **1**.
- UC windows per paired woman: mean 9.6, SD 7.9, minimum 1, Q1 4.00, median 6.5, Q3 13.00, maximum 32.
- IN windows per paired woman: mean 88.6, SD 50.2, minimum 15, Q1 52.75, median 81.5, Q3 112.75, maximum 228.

The earlier estimate of approximately 34 paired women is confirmed exactly.

## Methods

For each of the 34 primary FWH features, the analysis calculated each paired woman's raw UC-window median, raw IN-window median, and `UC median - IN median`. No window-level standardisation was applied. Each woman therefore contributed one paired difference per feature and received equal weight in the inferential analysis.

The primary test was a two-sided Wilcoxon signed-rank test. `zero_method="wilcox"` discarded exact zero differences before ranking. `correction=False` applied no continuity correction. `method="auto"` used the exact null distribution for n<=50 when there were no rank ties or zeros, and the asymptotic tie/zero-adjusted calculation otherwise. Exact zeros occurred for: `fwh_time_zero_crossing_count` (4), `fwh_time_slope_sign_change_count` (4), `fwh_freq_peak_frequency` (32), `fwh_complexity_lempel_ziv_complexity` (2). The result file records the effective nonzero sample size and p-value method for every feature.

The 95% CI is the 2.5th and 97.5th percentiles of 5,000 participant-level bootstrap medians. Each resample drew 34 women with replacement. All features used the same resampled participant indices, NumPy `default_rng`, and seed **20260914**. All 5,000 bootstrap estimates were finite for all features. `fwh_freq_peak_frequency` had a degenerate [0, 0] bootstrap interval because 32/34 paired differences were exactly zero; this reflects feature quantisation, not a computational failure.

BH correction was applied once across all 34 Wilcoxon p-values. Across the full family, 10/34 had raw p<0.05 and 5/34 had q<0.05.

## Original 11 results

| Feature | Median UC-IN | Mixed-consistent women | Wilcoxon p | BH q | Interpretation |
|---|---:|---:|---:|---:|---|
| fwh_time_median | -7.598e-05 | 64.7% | 0.01868 | 0.1058 | majority consistent; q>=0.05 |
| fwh_time_skewness | 0.03203 | 64.7% | 0.06433 | 0.1988 | majority consistent; q>=0.05 |
| fwh_time_mean_absolute_value | 0.00268 | 85.3% | 1.911e-05 | 1.763e-04 | same direction + q<0.05 |
| fwh_time_root_mean_square | 0.00388 | 85.3% | 9.052e-06 | 1.539e-04 | same direction + q<0.05 |
| fwh_time_peak_to_peak | 0.02209 | 91.2% | 5.138e-06 | 1.539e-04 | same direction + q<0.05 |
| fwh_time_waveform_length | 0.5114 | 82.4% | 2.118e-05 | 1.763e-04 | same direction + q<0.05 |
| fwh_freq_spectral_edge_frequency_95 | -0.01227 | 67.6% | 0.1629 | 0.3307 | majority consistent; q>=0.05 |
| fwh_freq_spectral_spread | -0.003775 | 61.8% | 0.08745 | 0.2286 | majority consistent; q>=0.05 |
| fwh_entropy_totalsampen | 0.9821 | 82.4% | 2.593e-05 | 1.763e-04 | same direction + q<0.05 |
| fwh_complexity_hjorth_complexity | -0.007742 | 67.6% | 0.09075 | 0.2286 | majority consistent; q>=0.05 |
| fwh_complexity_higuchi_fractal_dimension | -0.001184 | 61.8% | 0.1848 | 0.3307 | majority consistent; q>=0.05 |

All **11/11** original features retained the mixed-model direction. **5/11** retained BH q<0.05. All **11/11** had a strict majority of women changing in the mixed-model-consistent direction. Reversed original features: **none**.

## Core amplitude/time-domain results

| Feature | Paired women | Median UC-IN | Mixed-consistent women | Wilcoxon p | BH q | Bootstrap 95% CI |
|---|---:|---:|---:|---:|---:|---:|
| Mean absolute value | 34 | 0.00268 | 85.3% | 1.911e-05 | 1.763e-04 | [0.001887, 0.003518] |
| Root mean square | 34 | 0.00388 | 85.3% | 9.052e-06 | 1.539e-04 | [0.002126, 0.004908] |
| Peak-to-peak amplitude | 34 | 0.02209 | 91.2% | 5.138e-06 | 1.539e-04 | [0.01073, 0.02966] |
| Waveform length | 34 | 0.5114 | 82.4% | 2.118e-05 | 1.763e-04 | [0.3152, 0.6865] |
| Time skewness | 34 | 0.03203 | 64.7% | 0.06433 | 0.1988 | [-0.001526, 0.05798] |
| Time median | 34 | -7.598e-05 | 64.7% | 0.01868 | 0.1058 | [-2.188e-04, 2.199e-05] |
| Total SampEn | 34 | 0.9821 | 82.4% | 2.593e-05 | 1.763e-04 | [0.7865, 1.548] |

Mean absolute value, RMS, peak-to-peak amplitude, and waveform length were higher during UC for 85.3%, 85.3%, 91.2%, and 82.4% of paired women, respectively. All four retained BH q<0.05, and all four bootstrap intervals excluded zero. Time skewness remained positive in 64.7% but did not retain BH support. Time median remained negative in 64.7% but did not retain BH support. Total SampEn was higher during UC in 82.4%, with p=2.593e-05, q=1.763e-04, and a bootstrap interval excluding zero. This analysis establishes the within-woman pattern only; it does not reassess Total SampEn's relation to amplitude.

## Spectral and complexity features highlighted by the prior sensitivity

- `fwh_freq_spectral_edge_frequency_95`: median difference -0.01227, 67.6% mixed-model-consistent, p=0.1629, q=0.3307, bootstrap CI [-0.02649, -0.001491].
- `fwh_freq_spectral_spread`: median difference -0.003775, 61.8% mixed-model-consistent, p=0.08745, q=0.2286, bootstrap CI [-0.01035, 0.001351].
- `fwh_complexity_higuchi_fractal_dimension`: median difference -0.001184, 61.8% mixed-model-consistent, p=0.1848, q=0.3307, bootstrap CI [-0.004389, 3.323e-04].

All three returned to the negative mixed-model direction at the patient level, and 61.8%--67.6% of women were direction-consistent. None retained BH q<0.05. Spectral-edge frequency's percentile bootstrap interval excluded zero even though its Wilcoxon p-value and BH q-value did not. The interval estimates the sampling distribution of the median, whereas the signed-rank test uses signed ranks and assumes a symmetric difference distribution for its location interpretation; the two procedures need not agree in a small, asymmetric sample.

## Answers to the reviewer-focused questions

1. **Women with both states:** 34 of 42 women.
2. **Original directions retained:** 11/11.
3. **Original BH q<0.05 retained:** 5/11.
4. **Original reversals:** None.
5. **Core amplitude features:** MAV, RMS, peak-to-peak amplitude, and waveform length were consistently higher during UC and all retained family-wide BH support.
6. **Expected-direction percentages:** MAV 85.3%, RMS 85.3%, peak-to-peak 91.2%, and waveform length 82.4%.
7. **Total SampEn:** The patient-level median difference remained positive, 82.4% of women agreed with the mixed-model direction, and q<0.05 was retained.
8. **Spectral edge, spectral spread, and Higuchi fractal dimension:** All three had negative median differences, matching the mixed model rather than the pooled/cluster-robust OLS point estimates. Each had a direction-consistent majority, but none had q<0.05.
9. **Main pattern:** Yes. The patient-level results support a primarily amplitude/time-domain UC-vs-IN pattern. The four core amplitude/waveform features show the strongest consistency, with Total SampEn also prominent.
10. **Relation to the prior woman-cluster analysis:** This strengthens the interpretation of the previous result. It shows that the core effects are directly visible within women after removing window-count weighting, while confirming that inference for weaker spectral and complexity features is imprecise. The earlier cluster-robust OLS reversals are not reproduced by the paired medians, which suggests that pooled marginal weighting and between-woman composition contributed to those reversals. The analyses use different estimands and should be presented as complementary rather than treated as contradictory.

## Technical concerns

- Collapsing to medians changes the estimand from a window-level conditional contrast to the median of woman-level median differences. This analysis complements rather than replaces the mixed model.
- The paired sample contains only 34 women. UC counts range from 1 to 32, and three women have only one UC window. Their UC median is a single observation.
- Cluster sizes and state counts are highly unequal, but the paired analysis gives each woman equal inferential weight. This is desirable for the sensitivity question and discards information about within-woman precision.
- The Wilcoxon signed-rank test assumes a symmetric distribution of paired differences for a location-shift interpretation. The sign-consistency counts and bootstrap median interval provide complementary summaries when that assumption is doubtful.
- Discrete features generate ties and exact zeros. The specified zero handling and feature-specific p-value method are recorded explicitly.
- Percentile bootstrap intervals are simple and transparent but can be discrete with 34 participants. They should not be treated as multiplicity-adjusted intervals.

## Recommendation

The evidence is sufficiently coherent to proceed to cross-channel dependence analysis. The frozen mixed model, woman-cluster sensitivity, and patient-level paired analysis collectively support a robust core amplitude/time-domain pattern. The patient-level result also supports Total SampEn directionally and inferentially, without establishing independence from amplitude. Individual spectral and complexity findings remain less secure and should not carry the central claim.

## Scope stop

This analysis stops at the EHG9 patient-level UC-versus-IN sensitivity. It does not perform cross-channel modelling, alternative window-length analysis, SampEn-scale analysis, fetal-movement contamination analysis, or manuscript rewriting.
