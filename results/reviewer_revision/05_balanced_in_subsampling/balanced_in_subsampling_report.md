# Repeated balanced-IN sensitivity analysis

## Question and decision

This analysis asks whether the frozen EHG9 UC-versus-IN result is driven by the much larger inactivity pool. It uses the exact 34-feature primary FWH family, the frozen Stage-06 standardisation constants, and the frozen mixed-model hierarchy. Across 1,000 reproducible balanced samples, **11/11 original findings retained the frozen direction in at least 95% of usable repetitions, 0/11 retained BH q<0.05 in at least 80% of BH-evaluable repetitions, and 0/11 met both criteria with stable model fitting.**

The four core amplitude/time-domain effects were stable in direction and magnitude. None reached the prespecified ≥80% BH-support frequency threshold, so they are formally classified as directionally stable but inferentially variable. This is a precision result, not an effect-size reversal. Consequently, the central amplitude/time-domain pattern is **not plausibly an artefact of the excess IN window count, although its inferential strength is substantially weaker after discarding most IN data** under this sensitivity design.

## Population and feasibility

- Frozen UC/IN population: 333 UC and 3358 IN windows across 42 women; the event table contains 112 recordings in the UC/IN subset.
- Frozen UC pool retained per repetition: 333 windows.
- Same woman and same recording IN matches available: 301 (90.4% of UC).
- Same woman but different recording supplements: 25 (7.5%).
- UC windows without any eligible IN match: 7 (2.1%), all retained as required.
- Typical sample: 333 UC plus 326 IN rows; 35 women and a median of 87 recordings.
- Overall within-woman matching feasibility was 97.9%. This exceeded the prespecified feasibility requirement (unmatched proportion below 20.0%), so no unrelated-woman substitution was used and the analysis proceeded.

The primary sampling rule was hierarchical: draw without replacement within each woman-recording stratum up to that stratum's UC count, then fill recording-specific deficits using unused IN windows from the same woman in another recording. All UC observations remained in every fit. A separate woman-only 1:1 analysis was not run because the primary strategy already matched 97.9% of UC within woman while retaining the required recording hierarchy; adding another resampling family would not address a material feasibility failure.

## Frozen model reproduction and implementation

The input registry contained exactly 34 primary FWH features. Recalculation of BH q-values from the frozen Stage-06 p-values reproduced the stored q-values, and exactly the expected 11 features had q<0.05. For each feature, the event-unit mean and sample SD (ddof=1) reproduced the frozen `feature_mean_raw` and `feature_sd_raw` values.

Every balanced fit used `z_feature ~ state`, with IN (`baseline`) as reference, a woman random intercept, and a recording variance component nested in woman. REML and the frozen optimiser sequence (L-BFGS then BFGS) were retained; woman-only random-intercept fallback was permitted only after both nested attempts failed. Standardisation did **not** change by repetition: each outcome used the frozen all-state EHG9 mean and sample SD from the complete Stage-06 event table. Thus coefficients are directly comparable with the frozen estimates and changes cannot be attributed to rescaling.

Randomness was limited to IN selection. The master seed was `20260914` and each repetition used an independently derived PCG64 seed from `SeedSequence([master_seed, 1, iteration])`. Repetitions 1–3 were sampled twice before modelling and exact index/hash equality was verified.

## Results for the original 11

| Feature | Frozen β | Median balanced β (2.5th, 97.5th) | Direction | BH q<0.05 | Stability |
|---|---:|---:|---:|---:|---|
| Median | -0.219 | -0.198 (-0.272, -0.122) | 100.0% | 24.5% | directionally stable but inferentially variable |
| Skewness | 0.188 | 0.230 (0.131, 0.313) | 100.0% | 48.8% | directionally stable but inferentially variable |
| MAV | 0.310 | 0.298 (0.290, 0.308) | 100.0% | 0.0% | directionally stable but inferentially variable |
| RMS | 0.269 | 0.256 (0.251, 0.261) | 100.0% | 0.0% | directionally stable but inferentially variable |
| Peak-to-peak | 0.248 | 0.234 (0.228, 0.240) | 100.0% | 0.0% | directionally stable but inferentially variable |
| Waveform length | 0.302 | 0.291 (0.284, 0.298) | 100.0% | 0.0% | directionally stable but inferentially variable |
| Spectral edge 95 | -0.117 | -0.132 (-0.196, -0.070) | 100.0% | 28.2% | directionally stable but inferentially variable |
| Spectral spread | -0.126 | -0.131 (-0.191, -0.070) | 100.0% | 28.3% | directionally stable but inferentially variable |
| Total SampEn | 0.409 | 0.404 (0.390, 0.419) | 100.0% | 0.2% | directionally stable but inferentially variable |
| Hjorth complexity | -0.144 | -0.132 (-0.204, -0.060) | 100.0% | 7.5% | directionally stable but inferentially variable |
| Higuchi FD | -0.090 | -0.100 (-0.150, -0.052) | 100.0% | 25.1% | directionally stable but inferentially variable |

The empirical 2.5th–97.5th percentiles above describe variation across balanced resamples; they are not model-based confidence intervals. BH adjustment was recomputed separately within the 34-feature family in every repetition. The BH-support denominator includes only repetitions with all 34 finite p-values.

### Requested feature-level interpretations

- **Total SampEn:** median β 0.404; direction retained in 100.0%; BH q<0.05 in 0.2%; classified **directionally stable but inferentially variable**.
- **Spectral edge 95:** median β -0.132; direction retained in 100.0%; BH q<0.05 in 28.2%; classified **directionally stable but inferentially variable**.
- **Spectral spread:** median β -0.131; direction retained in 100.0%; BH q<0.05 in 28.3%; classified **directionally stable but inferentially variable**.
- **Higuchi FD:** median β -0.100; direction retained in 100.0%; BH q<0.05 in 25.1%; classified **directionally stable but inferentially variable**.

For Total SampEn, this analysis supports only the stability of its **state-associated coefficient under balanced IN sampling**. In light of the completed RMS-adjustment analysis, it should continue to be described as amplitude-coupled rather than as independent complexity; balanced resampling does not change that interpretation.

## Effect magnitude versus precision

Using a prespecified 20% descriptive threshold, 1/11 original features showed a material median point-estimate change and 11/11 showed a material median-SE change relative to the frozen full-IN fit. Skewness was the only original feature crossing the 20% point-estimate threshold (median ratio 1.221); it retained the same direction in every repetition. The `change_mechanism` field separates point-estimate movement from precision loss. A reduction in BH frequency without direction or magnitude change is interpreted as reduced precision from using far fewer IN observations, not disappearance of the effect.

The nested woman-plus-recording model was selected in 82.5% of all feature-repetition fits, and all 34,000 feature-repetition fits were usable. Complete 34-feature BH families were available in every repetition. Boundary, optimiser, warning, and fallback details are preserved for every fit in the convergence audit. A material technical caveat is that the prespecified woman-only fallback was selected in 95.0%–97.0% of repetitions for the four core amplitude features and 95.0% for Total SampEn because both nested optimisers often failed to converge after balancing. The four core features also used woman-only fallback in the frozen Stage-06 analysis, so their covariance structure is broadly comparable; Total SampEn changed from a clean frozen nested fit to mostly woman-only balanced fits and therefore deserves extra caution. These are usable prespecified fallback fits, not silent failures, but the large SE inflation (median 3.21–3.36 times frozen for those five features) limits claims about retained statistical significance.

## Reviewer-facing conclusions

1. **How feasible was within-recording matching?** 301/333 UC windows (90.4%) had same-recording IN capacity; same-woman other-recording supplementation increased coverage to 326/333 (97.9%).
2. **Are the original 11 directionally stable?** 11/11 retained direction in at least 95% of repetitions.
3. **How often do they retain BH q<0.05?** See the feature table; 0/11 did so in at least 80% of repetitions.
4. **Core amplitude/time-domain features:** All four were positive in every repetition and their median estimates were within 5.5% of frozen values. They did not retain frequent FDR support, owing to much larger SEs.
5. **Total SampEn:** Its balanced-sampling stability is shown above, but it remains an amplitude-coupled state association rather than an amplitude-independent complexity effect.
6. **Spectral edge, spectral spread, and Higuchi FD:** All three retained direction in 100% of repetitions and median magnitudes were within 13.3% of frozen estimates, but BH support occurred in only 28.2%, 28.3%, and 25.1%, respectively. They are directionally stable, not frequently FDR-supported after balancing.
7. **Precision versus magnitude:** The analysis explicitly distinguishes median coefficient shifts from SE inflation; loss of q<0.05 alone is not treated as effect disappearance.
8. **Is the central finding an imbalance artefact?** No empirical evidence from coefficient direction or magnitude supports that explanation for the core pattern; the excess IN pool mainly improved precision.

## Claims to retain, weaken, or remove

- **Retain:** directional and approximate-magnitude claims for features classified “very stable” or “directionally stable but inferentially variable,” with their resampling frequencies reported. The four core amplitude/time-domain effects satisfy the direction/magnitude criterion.
- **Weaken:** claims of universal statistical significance for any feature with variable BH support; describe stable direction/magnitude separately from reduced precision.
- **Remove/avoid:** any claim that Total SampEn demonstrates amplitude-independent complexity, and any assertion that balanced subsampling proves absence of all selection bias.

## Recommendation and stop point

**Proceed to the prespecified window-length reviewer analysis, while carrying forward the precision and fallback caveats.** The central coefficient pattern appears sufficiently stable with respect to excess IN windows under the specified within-woman/recording design. Statistical-significance claims should nevertheless be weakened because no original feature met the ≥80% BH-support-frequency criterion. This task stops here: no window-length, cross-channel, FM, entropy-scale, or manuscript-rewriting analysis was performed.
