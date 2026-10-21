# Corrected cross-channel UC-versus-IN analysis

## Bottom line

All 7/7 previously reported cross-channel-consistent features retained the same UC-versus-IN direction across EHG9-EHG12 after explicitly modelling shared physical windows. Under the preferred window-aware mixed model, 6/7 had descriptive p<0.05 in every channel and global tests identified FDR-supported channel-dependent magnitude for Total SampEn, Mean absolute value, Peak-to-peak amplitude, RMS, Waveform length. Under the required woman-cluster sensitivity, only 1/7 had p<0.05 in every channel and no global interaction survived its separate seven-test BH correction.

The direction-consistency claim is therefore preserved, but inferential strength is covariance-model dependent. A blanket claim of equal effect magnitude is not supported by Candidate B for the features with significant global interaction, while Candidate C does not confirm those interactions after FDR correction. The phrase **“cross-channel spatial consistency within the same cohort”** is supportable as a same-direction pattern, not as proof of uniform statistical support, channel-invariant magnitude, or independent independent cohorts.

## Locked analysis design

No post-Phase 1 redesign was made before examining results. The analysis used Candidate B exactly as recommended: `z_feature ~ state * channel`, with IN and EHG9 as references, a woman random intercept, and a physical-window variance component. Additive and interaction models used ML and the same random structure; their likelihood-ratio comparison is the 3-df global UC-by-channel test.

Candidate C was run for every feature using woman-cluster-robust OLS with CR1 finite-cluster correction and t/F inference based on 42 clusters. It served as the prespecified fallback if Candidate B failed and as a mandatory sensitivity otherwise. Recording was not separately parameterized in Candidate B because Phase 1 found the three-tier woman-recording-window model numerically unstable. Candidate C permits arbitrary dependence within woman, including dependence across recordings and windows.

Each feature was standardized once using its mean and sample SD (`ddof=1`) across all pooled UC/IN rows from EHG9-EHG12. Frozen channel-specific estimates used a different, within-channel three-state SD and should not be compared numerically without acknowledging that different denominator.

## Analysis population and predefined features

- Women: 42
- Recordings: 112
- Unique physical windows: 3,691
- Rows: 14,764
- UC windows: 333
- IN windows: 3,358
- Four-channel completeness: 100.0%

The seven-feature family was read from the frozen per-channel Stage 06 results as the features marked `primary_validated` on every channel: Total SampEn, Mean absolute value, Median, Peak-to-peak amplitude, RMS, Skewness, Waveform length. No supplementary 11-feature analysis was added.

## Convergence and fallback

Candidate B produced usable matched reduced/full fits for 7/7 features. Fallback was used as primary for 0/7. Boundary or singularity warnings occurred for 5/7 preferred fits and are retained in `convergence_fallback_audit.csv`; a small woman-intercept variance was accepted when the physical-window variance and fixed-effect covariance were finite. Total SampEn and peak-to-peak received the same prespecified diagnostics and fallback rules as every other feature.

## Corrected channel-specific UC-versus-IN effects

Values are standardized estimates with 95% CIs from the selected corrected model. An asterisk marks a descriptive, unadjusted channel-specific p<0.05. These p-values do not belong to the seven-test global-interaction FDR family.

| Feature | EHG9 | EHG10 | EHG11 | EHG12 | Global p | Global BH q |
|:--|--:|--:|--:|--:|--:|--:|
| Total SampEn | 0.256 [0.143, 0.368]* | 0.519 [0.406, 0.631]* | 0.211 [0.099, 0.324]* | 0.211 [0.099, 0.324]* | 3.17e-07 | 1.33e-06 |
| Mean absolute value | 0.171 [0.059, 0.284]* | 0.452 [0.340, 0.565]* | 0.147 [0.035, 0.260]* | 0.148 [0.035, 0.260]* | 7.51e-07 | 1.33e-06 |
| Median | -0.143 [-0.254, -0.032]* | -0.122 [-0.233, -0.011]* | -0.109 [-0.220, 0.001] | -0.142 [-0.253, -0.031]* | 0.815 | 0.815 |
| Peak-to-peak amplitude | 0.130 [0.018, 0.243]* | 0.422 [0.309, 0.534]* | 0.118 [0.005, 0.230]* | 0.116 [0.004, 0.229]* | 9.47e-07 | 1.33e-06 |
| RMS | 0.141 [0.028, 0.253]* | 0.432 [0.319, 0.544]* | 0.125 [0.013, 0.238]* | 0.124 [0.012, 0.237]* | 8.83e-07 | 1.33e-06 |
| Skewness | 0.188 [0.076, 0.300]* | 0.246 [0.133, 0.358]* | 0.223 [0.111, 0.335]* | 0.208 [0.096, 0.320]* | 0.555 | 0.648 |
| Waveform length | 0.150 [0.037, 0.262]* | 0.453 [0.341, 0.566]* | 0.131 [0.019, 0.244]* | 0.130 [0.018, 0.243]* | 7.43e-07 | 1.33e-06 |

## Interaction interpretation

- FDR-supported channel-dependent magnitude: Total SampEn, Mean absolute value, Peak-to-peak amplitude, RMS, Waveform length.
- No FDR-supported channel heterogeneity: Median, Skewness.
- In the separately adjusted seven-test woman-cluster sensitivity family, FDR-supported interactions: none.
- Pairwise EHG10/EHG11/EHG12-versus-EHG9 interaction estimates and unadjusted p-values are in `interaction_tests.csv`. They localize a global difference but do not replace the global FDR-controlled test.

Direction consistency, within-channel support, and equality of magnitude answer different questions. A feature can retain the same sign on every channel while showing a significant global interaction because one channel's effect is larger or smaller.

The five Candidate B interactions are localized primarily to a larger EHG10 effect. Their EHG10-minus-EHG9 point differences range from 0.263 to 0.303 SD, while the corresponding EHG11/EHG12 differences from EHG9 are small. Candidate C retains essentially the same point differences but assigns much larger uncertainty to EHG10-minus-EHG9, so none of its global tests survives BH correction.

## Statistical-support sensitivity

Channel lists use descriptive, unadjusted p<0.05. Global q-values use separate seven-feature BH families.

| Feature | Candidate B supported channels | Woman-cluster supported channels | Candidate B global q | Woman-cluster global q |
|:--|:--|:--|--:|--:|
| Total SampEn | EHG9, EHG10, EHG11, EHG12 | EHG9 | 1.33e-06 | 0.133 |
| Mean absolute value | EHG9, EHG10, EHG11, EHG12 | none | 1.33e-06 | 0.133 |
| Median | EHG9, EHG10, EHG12 | EHG12 | 0.815 | 0.771 |
| Peak-to-peak amplitude | EHG9, EHG10, EHG11, EHG12 | none | 1.33e-06 | 0.133 |
| RMS | EHG9, EHG10, EHG11, EHG12 | none | 1.33e-06 | 0.133 |
| Skewness | EHG9, EHG10, EHG11, EHG12 | EHG9, EHG10, EHG11, EHG12 | 0.648 | 0.466 |
| Waveform length | EHG9, EHG10, EHG11, EHG12 | none | 1.33e-06 | 0.133 |

## Features requiring particular attention

- **Total SampEn:** direction was consistent across all four channels; 4/4 channel-specific contrasts had descriptive p<0.05; the corrected global interaction was p=3.17e-07, q=1.33e-06. The preferred model was usable; the woman-cluster sensitivity supported 1/4 channel-specific contrasts and its global p was 0.0896 (separate-family q=0.133).
- **Peak-to-peak amplitude:** direction was consistent across all four channels; 4/4 channel-specific contrasts had descriptive p<0.05; the corrected global interaction was p=9.47e-07, q=1.33e-06. The preferred model was usable; the woman-cluster sensitivity supported 0/4 channel-specific contrasts and its global p was 0.0423 (separate-family q=0.133).

## Answers to the reviewer questions

1. **Direction:** 7/7 retained a common UC-versus-IN direction across all four channels.
2. **Within-channel support:** Candidate B supported all four channels for 6/7 features; median was supported on EHG9, EHG10, and EHG12 but not EHG11. Candidate C supported all four channels only for skewness; Total SampEn was supported on EHG9, median on EHG12, and the other channel-specific contrasts were not significant. These are descriptive, unadjusted p-values.
3. **Significant global interactions:** Candidate B identified Total SampEn, Mean absolute value, Peak-to-peak amplitude, RMS, Waveform length after BH correction. Candidate C identified none after its separate BH correction.
4. **Channel-dependent magnitude:** Candidate B supports channel dependence for the five-feature set above, driven mainly by larger EHG10 effects. Candidate C does not confirm these as FDR-significant, so magnitude heterogeneity is model-dependent.
5. **Similar magnitude:** Candidate B found no FDR-supported heterogeneity for Median, Skewness; Candidate C found no FDR-supported heterogeneity for any feature. Absence of detected heterogeneity is not proof of exact equality.
6. **Total SampEn:** summarized above. Its entropy/amplitude interpretation is outside this task; only cross-channel dependence is assessed here.
7. **Peak-to-peak amplitude:** summarized above.
8. **Material scientific change:** the cross-channel direction pattern is preserved, but the strength of channel-specific support and evidence for unequal magnitude depend materially on covariance assumptions. Candidate B resolves the shared-window omission; Candidate C shows that stronger protection against all within-woman dependence substantially widens uncertainty.
9. **“Cross-channel spatial consistency within the same cohort”:** supported for direction. It should not be used without qualification to imply that every channel remains statistically significant under woman-cluster inference, that magnitudes are uniform, or that these are independent replications.
10. **Claims to weaken:** avoid claiming that all seven are robustly significant on all four channels, or that the five Candidate B interactions conclusively establish channel-dependent magnitude without noting the negative FDR-controlled cluster sensitivity. Also avoid “equivalent,” “uniform,” and “channel-independent,” and do not treat four channels from the same windows and women as independent replications.

## Multiplicity

BH-FDR was applied to exactly seven global UC-by-channel p-values selected by the prespecified primary hierarchy. The woman-cluster global tests were also BH-adjusted as a separate seven-test sensitivity family. Channel-specific and pairwise interaction p-values are explicitly descriptive and unadjusted. Neither family was mixed with the original 34-feature EHG9 primary family.

## Reproducibility

- Script: `/Users/shivangayathri/Work/EHGProcessing/ehg-sampen-state-associations/src/cross_channel.py`
- Phase 1 design source: `woman intercept + physical-window variance component; woman-cluster sensitivity`
- Frozen feature/effect sources: `/Users/shivangayathri/Work/EHGProcessing/ehg-sampen-state-associations/work/primary/mixed_model_results_60s_fwh_EHG9.csv` through `EHG12.csv`
- Frozen analysis population sources: `/Users/shivangayathri/Work/EHGProcessing/ehg-sampen-state-associations/work/primary/event_unit_analysis_60s_EHG9.csv` through `EHG12.csv`
- Random seed: none; all procedures are deterministic.
- Package versions: `package_versions.txt`.
