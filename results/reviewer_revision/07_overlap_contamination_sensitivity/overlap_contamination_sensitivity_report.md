# Overlap-contamination sensitivity report

## Main findings

The retained 60-s EHG9 population contained 333 UC, 392 FM, and 3358 IN analysis units. The frozen metadata identified 21 UC windows (6.3%) overlapping an FM-centred window and 21 FM windows (5.4%) overlapping a UC-centred window. Reconstruction from all annotations identified 22 unique UC-FM overlap relations: 20 one-to-one pairs plus one UC window overlapping two FM windows, yielding 21 contaminated windows in each class.

After excluding contaminated UC units, **11/11 original findings retained the frozen direction and 9/11 retained BH q<0.05**. The original findings no longer below q=0.05 were Spectral edge 95, Higuchi fractal dimension. Their clean coefficients changed by only 12.0% and 3.1% in absolute magnitude, respectively, while their standard errors increased modestly; these threshold crossings therefore indicate reduced precision rather than disappearance or reversal. Sign reversals: none. Features meeting the prespecified descriptive definition of at least 20% attenuation in absolute coefficient magnitude: none.

All four core amplitude/waveform coefficients remained positive (4/4) and FDR-supported (4/4). Excluding competing FM overlap therefore did not materially weaken the central amplitude/time-domain conclusion.

For the four SampEn-profile features currently identified for the manuscript's secondary FM-versus-IN result, 4/4 retained direction and 4/4 retained BH q<0.05 after contaminated FM exclusion. Across the full frozen 34-feature family, all 11 frozen FM findings remained supported and the clean analysis also contained 11 q<0.05 features. This directly addresses the reported FM family while supplying the complete comparison to avoid selective reporting.

## Locked definitions

Retained event unit has overlap_duration_s > 0 and its frozen overlapping_event_class is the competing class (FM for UC; UC for FM); contains_competing_class is also checked, although all retained event units are False because contained competing events were labelled mixed at Stage 02 and excluded before Stage 06.

Same-class overlaps were retained. Participant/equipment/position annotations and the movement-artifact flag were not used as exclusion criteria. An all-annotation reconstruction verified that the greatest-overlap metadata did not miss a competing UC/FM event for any retained unit.

## Population and analyzability

- Clean UC versus IN: 312 clean UC and 3358 IN windows.
- Clean FM versus IN: 371 clean FM and 3358 IN windows.
- Exact woman, recording, overlap-duration, overlap-fraction, and contribution counts are in `overlap_population_summary.csv`.
- Exclusion removed 6.3% of UC and 5.4% of FM units, not a severe loss of event observations. Women contributing both relevant states are reported explicitly rather than inferred from total women.

## Clean UC focus features

| Feature | Frozen beta | Clean-UC beta | Absolute change | Change in absolute magnitude | Clean p | Clean q | Direction |
|---|---:|---:|---:|---:|---:|---:|---|
| Mean absolute value | 0.310 | 0.323 | 0.013 | 4.1% | 4.79e-08 | 8.15e-07 | positive |
| RMS | 0.269 | 0.283 | 0.014 | 5.1% | 1.76e-06 | 1.49e-05 | positive |
| Peak-to-peak amplitude | 0.248 | 0.262 | 0.014 | 5.6% | 9.8e-06 | 6.66e-05 | positive |
| Waveform length | 0.302 | 0.316 | 0.014 | 4.6% | 8.99e-08 | 1.02e-06 | positive |
| Time median | -0.219 | -0.219 | 0.000 | -0.2% | 0.000147 | 0.00083 | negative |
| Time skewness | 0.188 | 0.180 | 0.008 | -4.4% | 0.002 | 0.012 | positive |
| Total SampEn | 0.409 | 0.422 | 0.013 | 3.2% | 7.95e-13 | 2.7e-11 | positive |

The 20% threshold is descriptive, not an inferential test. A q-value crossing 0.05 is interpreted together with direction and coefficient magnitude.

Total SampEn is reported only as an amplitude-coupled, scale-sensitive state-associated feature, consistent with the completed RMS/scale analyses. This overlap analysis does not establish amplitude-independent complexity.

## Clean FM sensitivity

The current manuscript plan explicitly describes Median SampEn, SD SampEn, Skewness SampEn, and Kurtosis SampEn as the secondary FM-versus-IN profile family. The model was nevertheless fitted to the frozen 34-feature family and BH-FDR was applied over all 34 FM-versus-IN p-values, matching Stage 06.

| Reported FM feature | Frozen beta | Clean-FM beta | Frozen q | Clean q | Same direction |
|---|---:|---:|---:|---:|---|
| Median SampEn | -0.154 | -0.157 | 0.015 | 0.019 | True |
| SD SampEn | -0.209 | -0.198 | 0.000194 | 0.000818 | True |
| Skewness SampEn | 0.219 | 0.207 | 0.000449 | 0.002 | True |
| Kurtosis SampEn | 0.184 | 0.171 | 0.005 | 0.015 | True |

The clean-FM results show how much the reported association depends on overlapping UC activity. Direction, magnitude, and multiplicity-adjusted support are kept separate in the output tables.

## Model and convergence

Both sensitivities preserve the exact frozen Stage 06 three-state model: standardized feature as outcome, IN reference, woman random intercept, recording variance component, REML, and the same nested/woman-only optimizer fallback sequence. The clean-UC fit removes only contaminated UC units and keeps all FM and IN; the clean-FM fit removes only contaminated FM units and keeps all UC and IN. Standardization is recomputed over the retained three-state observations for each feature, exactly as in Stage 06.

Usable clean-UC models: 34/34. Usable clean-FM models: 34/34. Woman-only fallbacks: 5/34 and 2/34, respectively. Boundary status, optimizer, variance estimates, and complete warning text are retained in `overlap_model_convergence_audit.csv`.

The clean-UC Total SampEn fit moved from the frozen nested tier to the prespecified woman-only fallback, while its coefficient remained close to the frozen value. Katz fractal dimension moved from the frozen woman-only tier to a clean nested fit. In the clean-FM scenario, RMS and peak-to-peak moved from frozen woman-only boundary fits to clean nested fits, and Katz fractal dimension also moved to the nested tier. These changes follow the frozen fallback algorithm and are disclosed rather than silently treated as identical covariance fits.

## Interpretation

This sensitivity separates two distinct situations. Competing event centres within a window had already been removed as `mixed` by the frozen pipeline. The additional exclusion here removes retained UC/FM windows whose 60-s footprints overlap a competing event-centred window even though the competing event centre lies outside the anchor's own ±30-s interval.

The primary UC conclusion should be considered robust to this defined source of FM overlap if direction and magnitude remain stable; loss of precision alone would not imply disappearance. The secondary FM conclusion should be retained only for features whose direction and magnitude remain credible after removal, with FDR support reported explicitly.

## Scope

No frozen Stage 00--09 artifact or earlier reviewer-revision output was modified. No window-length, cross-channel, entropy-scale, classification, or manuscript-rewriting analysis was performed.
