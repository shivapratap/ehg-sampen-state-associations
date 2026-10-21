# EHG9 temporal-dependence sensitivity analysis

## Question

Do the frozen 11/34 EHG9 uterine-contraction (UC) versus inactivity (IN) findings remain stable when inference is robust to dependence among repeated observations from the same woman?

## Frozen inputs and reproduction checks

- `work/primary/confirmatory_feature_registry_60s_fwh.csv` (SHA-256 `bdbd34522dfd5d64fa852fb90ab734cca03eec5cb3aa11771b95b5997450ae12`)
- `work/primary/event_unit_analysis_60s_EHG9.csv` (SHA-256 `b9d383031470aa718c63210ae4f0f9ab7d3fd13a5a009e1a405a3d5e16e67d68`)
- `work/primary/mixed_model_results_60s_fwh_EHG9.csv` (SHA-256 `dcfc060c1b0311b72eed8f81adc6c0a57312787566368a2449d9562df7ad8638`)
- Feature-family check: 34 primary features and expected 11 supported effects

The Stage 06 registry contains exactly 34 unique primary FWH features. The EHG9 mixed-model file contains one UC-vs-IN result for every registered feature. Reapplying Benjamini-Hochberg (BH) to the 34 frozen p-values reproduces every stored q-value within numerical tolerance and reproduces the exact 11-feature q<0.05 set. The manuscript update memo states 11/34 and identifies SVD entropy, fixed-r SampEn, and Hjorth mobility as the three prior borderline findings removed by the expanded family. The frozen results agree with those statements. The supplied updated memo does not tabulate every numerical estimate and interval, so manuscript-level numerical identity beyond the count and named feature-set change cannot be checked from that document.

No primary mixed model was refitted or altered.

## Methods

The analysis used the frozen EHG9 event-unit table and retained only UC and IN: 3,691 observations (333 UC and 3,358 IN), 42 women, and 112 recordings. Of the 42 women, 34 contributed both states, 7 contributed IN only, and 1 contributed UC only. Woman-level observation counts ranged from 7 to 233.

For exact comparability with Stage 06, each feature used its frozen all-three-state mean and sample SD (ddof=1) from the primary fit. OLS was then fit as `z_feature ~ state_uc`, with IN as 0 and UC as 1. Pooled inference used ordinary OLS standard errors. Woman-cluster inference used the identical point estimate with statsmodels' CR1 covariance correction, `G/(G-1) * (N-1)/(N-K)`, cluster degrees of freedom `G-1 = 41`, and two-sided t-based 95% CIs and p-values. BH was applied separately to the 34 pooled p-values and the 34 woman-cluster p-values. No randomness was used.

## Results for the original 11

| Feature | Mixed estimate | Clustered estimate | Cluster-robust SE | Cluster p | Cluster q | Direction vs mixed |
|---|---:|---:|---:|---:|---:|---|
| fwh_time_median | -0.219 | -0.146 | 0.074 | 0.0563 | 0.319 | same |
| fwh_time_skewness | 0.188 | 0.181 | 0.052 | 0.00118 | 0.0402 | same |
| fwh_time_mean_absolute_value | 0.310 | 0.296 | 0.191 | 0.129 | 0.488 | same |
| fwh_time_root_mean_square | 0.269 | 0.256 | 0.192 | 0.19 | 0.515 | same |
| fwh_time_peak_to_peak | 0.248 | 0.235 | 0.192 | 0.227 | 0.515 | same |
| fwh_time_waveform_length | 0.302 | 0.298 | 0.191 | 0.126 | 0.488 | same |
| fwh_freq_spectral_edge_frequency_95 | -0.117 | 0.029 | 0.150 | 0.849 | 0.88 | reversed |
| fwh_freq_spectral_spread | -0.126 | 0.014 | 0.168 | 0.934 | 0.934 | reversed |
| fwh_entropy_totalsampen | 0.409 | 0.412 | 0.186 | 0.0323 | 0.275 | same |
| fwh_complexity_hjorth_complexity | -0.144 | -0.203 | 0.082 | 0.0179 | 0.217 | same |
| fwh_complexity_higuchi_fractal_dimension | -0.090 | 0.110 | 0.170 | 0.521 | 0.674 | reversed |

1. **Frozen reproduction:** Yes. The registry, 34 mixed-model UC-vs-IN rows, stored BH q-values, and exact 11-feature set reproduce correctly. The count and named FDR-family change agree with the manuscript update memo.
2. **Effect direction:** 8/11 retained the primary mixed-model direction under woman-cluster-robust OLS.
3. **BH support:** 1/11 retained woman-cluster BH q<0.05. This was `fwh_time_skewness`.
4. **Direction reversals:** 3 original features reversed direction: `fwh_freq_spectral_edge_frequency_95`, `fwh_freq_spectral_spread`, `fwh_complexity_higuchi_fractal_dimension`. These reversals reflect the difference between the unadjusted pooled mean contrast and the mixed-model coefficient. Cluster-robust covariance changes the uncertainty, not the OLS point estimate.
5. **Magnitude changes:** For this report, material attenuation was defined descriptively as at least 20% lower absolute coefficient than the mixed estimate. By that threshold, 3 original effects were materially attenuated: `fwh_time_median`, `fwh_freq_spectral_edge_frequency_95`, `fwh_freq_spectral_spread`. The two spectral estimates moved through zero to small positive values. The median estimate retained its negative direction but was about 34% smaller in absolute magnitude. This threshold is not an inferential test.
6. **Amplitude/time-domain pattern:** The point-estimate pattern was preserved. All 6/6 original time-domain effects kept their direction. Mean absolute value, RMS, peak-to-peak amplitude, and waveform length changed by less than 6% in magnitude; skewness changed by less than 4%; median attenuated as noted above. Reduced FDR support primarily reflects much larger cluster-robust uncertainty, not disappearance of the core amplitude/waveform point estimates.
7. **P-value comparison:** Across all 34 features, nominal p<0.05 counts were mixed 17, pooled 22, and woman-cluster 5; BH q<0.05 counts were mixed 11, pooled 21, and woman-cluster 1. Woman-cluster p-values exceeded pooled p-values for 33/34 features, with a median cluster-to-pooled p-value ratio of 15.1. Woman-cluster p-values exceeded mixed-model p-values for 25/34 features, with a median ratio of 4.4. Pooled and mixed p-values were not uniformly ordered because the mixed model changes both the coefficient and the variance model.
8. **Pseudoreplication:** Yes, the pooled-versus-cluster comparison provides strong empirical evidence that treating 3,691 windows as independent overstates precision. Pooled BH support was 21/34, but woman-cluster BH support was 1/34, while the pooled and clustered coefficients were exactly identical. The difference therefore comes entirely from dependence-robust uncertainty. This does not invalidate the mixed model, which already modeled woman and, where estimable, recording effects.
9. **Appropriateness of woman clustering:** Woman is the correct primary independence unit because recordings are nested within women. Forty-two clusters are workable but modest, so finite-cluster correction is important and tail-area inference remains approximate. Severe cluster-size imbalance (range 7--233 observations), only 34 women with both states, and eight single-state women can make the marginal OLS contrast sensitive to cluster composition. Cluster-robust OLS also leaves the pooled coefficient unchanged and does not isolate a purely within-woman contrast or adjust for recording-level heterogeneity. These are limitations of this sensitivity estimand, not reasons to cluster by recording instead. No dataset feature makes woman clustering technically inappropriate.

## Conclusion and recommendation

The result is stable in its main amplitude/time-domain point-estimate pattern but not in multiplicity-adjusted cluster-robust precision. Eight of the original 11 keep direction, all six original time-domain effects keep direction, and the four core amplitude/waveform coefficients remain close to the mixed-model estimates. Only one original feature retains BH q<0.05 under the conservative woman-cluster OLS analysis. Three non-core features reverse direction, including two spectral effects that shrink close to zero.

The primary finding appears sufficiently stable to proceed with the revision only if it is framed as a predominantly amplitude/time-domain pattern supported by the prespecified mixed model, with the dependence-robust analysis presented as evidence of reduced precision and weaker support for individual features. It would be too strong to claim that all 11 individual findings are robust under woman-cluster inference.

## Scope

This analysis stops at EHG9 UC versus IN. It does not analyze cross-channel models, alternative window lengths, SampEn scale sensitivity, fetal-movement contamination, recording-clustered inference, or manuscript wording.
