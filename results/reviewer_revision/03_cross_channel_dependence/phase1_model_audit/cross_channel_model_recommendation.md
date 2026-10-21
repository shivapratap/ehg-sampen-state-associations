# Cross-channel dependence: Phase 1 model audit and feasibility assessment

## Bottom line

Candidate B is the recommended mixed-model specification because Candidate A was not stable across all three feasibility features. The revised primary estimand should be the UC-versus-IN channel heterogeneity test, not the frozen UC/FM omnibus. Phase 2 is a **GO**, conditional on using the prespecified failure hierarchy below and retaining all convergence diagnostics.

## 1. What the frozen analysis did

The frozen Step 07 pooled four Stage 06 event-unit tables and fitted standardized feature values using state, channel, and state-by-channel fixed effects, with EHG9 and baseline inactivity as references. Standardization used the feature-specific pooled mean and sample SD over all three states and all four channels. Woman was the top-level random-intercept grouping variable; recording was included as a variance component when that tier fitted. `window_id` was present in the source data but absent from the model. Thus, the model did not explicitly represent the correlation among EHG9-EHG12 measurements from the same physiological 60-s interval.

The frozen likelihood-ratio test had 6 df because it jointly tested both FM-by-channel and UC-by-channel terms. That is broader than the reviewer's direct UC-versus-IN question. Of seven features, five produced finite tests; Total SampEn and peak-to-peak failed with non-finite likelihoods after fallback (`fwh_entropy_totalsampen, fwh_time_peak_to_peak`).

The frozen model was incomplete as a covariance model, not wrong in its construction of fixed effects. Its point estimand can still be meaningful, but its standard errors and likelihood-ratio test rely on a residual-independence assumption that is false within physical windows.

## 2. Verified multichannel structure

- Source: the four frozen Stage 06 `event_unit_analysis_60s_EHG9`-`EHG12` tables.
- Rows: 16,332; unique physical windows: 4,083.
- Rows per physical window: 4 to 4 (expected 4).
- Windows with all four channels: 100.0%.
- Women: 42; recordings: 112.
- Window counts: UC 333, IN 3,358, FM 392.
- Channel availability is complete and balanced: True.
- Each `window_id` maps to exactly one woman, one recording, and one state. It is globally unique in these canonical tables, so no composite window key is required. The script nevertheless uses an explicit `physical_window_key` copy and a woman-recording composite key to make nesting unambiguous.

## 3. Reviewer-response estimand and coding

Restrict the primary revised model to UC and IN. For each feature, standardize over the pooled UC/IN four-channel analysis rows using the sample standard deviation (`ddof=1`), then fit:

`z_feature ~ C(state, Treatment(reference="IN")) * C(channel, Treatment(reference="EHG9"))`

The UC coefficient is the UC-minus-IN contrast on EHG9. Each UC-by-channel coefficient is how the corresponding channel's UC-minus-IN contrast differs from EHG9. The primary heterogeneity test is a joint 3-df test that the EHG10, EHG11, and EHG12 interaction coefficients are all zero. This restriction changes the question from the original UC/FM omnibus to the reviewer-requested UC-specific estimand; it does not otherwise change the scientific contrast.

This UC/IN-only standardization differs from the frozen three-state scaling population. It changes the numerical SD unit but not the unstandardized contrast or the null hypothesis of no channel heterogeneity; the final report should not compare standardized coefficients across the two scalings as though they used the same denominator.

For the corrected mixed models, compare the additive reduced model with the interaction model by ML (`reml=False`) while holding the random structure fixed. Report the 3-df likelihood-ratio test plus the three interaction estimates and confidence intervals. REML must not be used for this fixed-effects likelihood-ratio comparison.

## 4. Candidate assessment

### Candidate A: woman + recording-within-woman + physical window

Specification: woman random intercept through `groups=patient_id`; recording-within-woman variance component `0 + C(recording_key)`; physical-window variance component `0 + C(physical_window_key)`. This is statistically coherent because windows are nested in recordings, recordings in women, and the four channels repeat within each window. The design contains the replication needed in principle to separate window from recording, and the formulation is accepted by the installed statsmodels `MixedLM`; however, numerical estimability must be judged from the actual fits rather than from design alone.

Feasibility: 0/3 matched reduced/full pairs converged to finite likelihoods under L-BFGS or BFGS. Because none yielded a usable pair, the feasibility run could not confirm stable numerical separation of recording and physical-window variance. The counts of retained positive physical-window (0/3) and recording (0/3) estimates are therefore not evidence that either component is absent; no failed fit was interpreted. Candidate A should not be the routine Phase 2 specification in the current environment.

### Candidate B: woman + physical window

Specification: woman random intercept plus physical-window variance component, omitting recording. This preserves the fixed-effect estimand and directly handles shared-window correlation, but it no longer separately parameterizes between-recording heterogeneity within a woman. Feasibility: 3/3 matched pairs converged, and each returned a positive physical-window variance. All 3/3 emitted a boundary warning, driven by a very small woman-intercept variance; this requires transparent reporting and makes the woman-cluster sensitivity essential. Given Candidate A's systematic failure, Candidate B is the recommended Phase 2 mixed model.

### Candidate C: woman-cluster-robust OLS

Specification: the same fixed-effects formula with CR1 finite-cluster correction and reference-distribution adjustment based on 42 woman clusters. It allows arbitrary correlation within each woman, covering channel, recording, window, and residual temporal dependence without separately estimating their variance components. It is statistically coherent at 42 clusters, but less efficient and less explanatory about the covariance hierarchy. It should be a prespecified fallback and sensitivity analysis, not automatically the primary model. Feasibility: 3/3 fits succeeded and all UC-by-channel terms were extractable.

## 5. Minimal feasibility results (not final scientific results)

No multiplicity correction was applied. The test column is an LR chi-square for Candidates A/B and a cluster-robust F statistic for Candidate C.

| Feature | Candidate | Status | Window variance | Recording variance | Test statistic | p | Interactions extracted |
|:--|:--|:--|--:|--:|--:|--:|:--|
| RMS | A_woman_recording_window_MixedLM | failed | NA | NA | NA | NA | no |
| RMS | B_woman_window_MixedLM | converged | 0.3245 | NA | 30.92 | 8.833e-07 | yes |
| RMS | C_woman_cluster_robust_OLS | converged | NA | NA | 2.636 | 0.06247 | yes |
| Total SampEn | A_woman_recording_window_MixedLM | failed | NA | NA | NA | NA | no |
| Total SampEn | B_woman_window_MixedLM | converged | 0.384 | NA | 33.03 | 3.169e-07 | yes |
| Total SampEn | C_woman_cluster_robust_OLS | converged | NA | NA | 2.319 | 0.08958 | yes |
| Peak-to-peak amplitude | A_woman_recording_window_MixedLM | failed | NA | NA | NA | NA | no |
| Peak-to-peak amplitude | B_woman_window_MixedLM | converged | 0.3295 | NA | 30.78 | 9.468e-07 | yes |
| Peak-to-peak amplitude | C_woman_cluster_robust_OLS | converged | NA | NA | 2.981 | 0.04231 | yes |

These three-feature results establish computational feasibility only. They must not be presented as the final cross-channel analysis or used to select features/models by significance.

The Candidate B likelihood-ratio p-values were much smaller than the woman-cluster-robust p-values in all three examples (Candidate C ranged from 0.042 to 0.090). This does not decide the scientific result at Phase 1, but it demonstrates that inference is sensitive to the covariance assumptions. Phase 2 must therefore report the Candidate C sensitivity alongside Candidate B and distinguish coefficient stability from changes in standard errors or reference distributions.

## 6. Final model recommendation

1. **Primary model:** Candidate B in the current Python/statsmodels environment: woman random intercept plus a physical-window variance component. Candidate A is scientifically attractive but failed for 3/3 representative features with both prespecified optimizers, so it is not recommended as the routine final specification unless a different validated implementation is introduced before Phase 2.
2. **Population:** UC and IN only for the primary reviewer-response analysis. FM may remain a separate secondary question, but it must not be included in this primary 3-df test.
3. **References:** IN for state and EHG9 for channel.
4. **Test:** a single joint 3-df ML likelihood-ratio test per feature for the three UC-by-channel coefficients, accompanied by all three coefficient estimates, 95% CIs, and channel-specific UC-minus-IN contrasts.
5. **Failure hierarchy:** fit Candidate B with a documented optimizer sequence. If either member of the matched reduced/full pair is non-finite, non-convergent, or singular in a way that invalidates inference, use Candidate C with woman-cluster robust covariance and a small-sample correction for that feature. Never mix random structures inside one LR comparison. Candidate A may be reported as an attempted, unsuccessful richer structure rather than silently dropped.
6. **Sensitivity:** run Candidate C for every final feature even when Candidate B succeeds, so conclusions can be checked against a covariance estimator robust to all within-woman dependence. If Candidate B and Candidate C differ, emphasize estimates and uncertainty under both rather than selecting the smaller p-value.
7. **Feature family for Phase 2:** analyse the same seven prespecified cross-channel features from the frozen Step 07 analysis: Total SampEn, mean absolute value, median, peak-to-peak, RMS, skewness, and waveform length. Do not expand or select this family based on the three feasibility fits. Apply the planned family-wise BH procedure only in Phase 2, not here.

## 7. Technical concerns and go/no-go

- A random intercept for physical window imposes exchangeable residual correlation across the four channels after the higher tiers. It directly addresses the reviewer concern, but it does not model an arbitrary 4x4 within-window residual covariance. The woman-cluster sensitivity protects against misspecification of that exchangeable form.
- Variance components approached the boundary in Candidate B because the woman-intercept variance was small after adding the window tier. This is not a reason to delete the window component; it is a reason to retain diagnostics and the woman-cluster sensitivity.
- Total SampEn and peak-to-peak previously produced infinite likelihoods; their corrected fits and warnings should receive feature-specific scrutiny in Phase 2. A successful feasibility fit does not guarantee every optimizer or all seven features will be stable.
- With 42 women, woman-cluster robust inference is reasonable but not asymptotically luxurious; retain the CR1/df correction and avoid clustering by recording as the primary fallback because recordings are nested within women.
- All four channels are present for every physical window, so missing-channel imbalance is not a concern in these frozen inputs.

**Phase 2 recommendation: GO, with conditions.** The physical-window unit is verified, the intended UC-specific estimand is identifiable, and Candidate B plus Candidate C are operational for every feasibility feature. The final run should use Candidate B as the mixed model, require Candidate C for all seven features, and stop for manual review if a mixed fit is non-finite or diagnostically singular rather than silently accepting it. The large difference between model-based and cluster-robust feasibility p-values precludes an unconditional GO based on Candidate B alone.

## Reproducibility

- Script: `<repo>/reviewer_revision/03_cross_channel_dependence/phase1_model_audit/cross_channel_model_audit.py`
- Frozen current-analysis sources: `<repo>/src/additional_analysis/step07_channel_state_interaction.py`, `<repo>/results/paper1/additional_analysis/phase3/step07_channel_interaction_results.csv`, `<repo>/results/paper1/additional_analysis/phase3/step07_channel_interaction_report.md`
- Manuscript/planning sources inspected: `<repo>/manuscript/Paper1_Manuscript_Blueprint.md`, `<repo>/manuscript/Paper1_Additional_Analysis_Plan.md`
- Random seed: none; all procedures are deterministic.
- Package versions are recorded in `package_versions.txt`.
