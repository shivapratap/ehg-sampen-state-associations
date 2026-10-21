# Window-length sensitivity report

## Main finding

The 40-s and 80-s reruns were technically feasible and used the current frozen pipeline logic with only duration-dependent segmentation changed. Among the frozen original 11, **11/11 retained the same UC-versus-IN direction at 40, 60, and 80 s; 10/11 had BH q<0.05 at all three durations, and 11/11 had BH q<0.05 at two or more durations.**

The central amplitude/time-domain pattern was positive at every duration, and 3/4 core features met the prespecified magnitude-stability rule. The central-pattern rule required all four positive and at least three magnitude-stable, with any duration-additive exception named. The evidence therefore does not indicate that the central pattern is an artefact of choosing 60 s.

## Technical feasibility and population

The regenerated 40-s and 80-s manifests exactly matched the legacy Module 1 manifests by values and SHA-256. Feature extraction used EHG9 only, full-record filtering, and 800 or 1,600 decimated samples per window. No feature algorithm or model setting was tuned by duration.

| Duration | Manifest UC/IN/FM | Final analysis UC/IN/FM | Total manifest | Women | Recordings | Boundary-rejected events | Mixed windows |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 s | 344/5200/406 | 344/5200/406 | 6127 | 42 | 112 | 5 | 108 |
| 60 s | 333/3358/392 | 333/3358/392 | 4283 | 42 | 112 | 7 | 132 |
| 80 s | 328/2470/382 | 328/2470/382 | 3391 | 42 | 112 | 11 | 145 |

Changing duration mainly changes the number of tiled IN windows. UC and FM counts also decline with longer windows because more event windows cross edge limits or contain a competing class. Windows are never truncated.

## Comparison scale and model

Primary coefficients use the frozen 60-s feature mean and sample SD at all durations. This makes the common-scale estimates directly comparable. The CSV also reports raw-unit and within-duration standardized estimates because waveform length, Total SampEn, and other distributions can depend on sample count. The model otherwise matches frozen Stage 06: all three states, IN reference, REML, woman random intercept, recording variance component, and the unchanged optimizer/fallback hierarchy. BH-FDR was applied separately across exactly 34 UC-versus-IN p-values per duration.

## Focus features

| Feature | β 40 s | β 60 s | β 80 s | q 40 s | q 60 s | q 80 s | Direction/magnitude assessment |
|---|---:|---:|---:|---:|---:|---:|---|
| Median | -0.288 | -0.219 | -0.203 | 6.67e-05 | 0.0005 | 0.000247 | direction stable; magnitude variable |
| Skewness | 0.188 | 0.188 | 0.182 | 0.0029 | 0.0055 | 0.00997 | direction and magnitude stable |
| MAV | 0.392 | 0.310 | 0.268 | 2.03e-08 | 1.12e-06 | 1.59e-05 | direction stable; magnitude variable |
| RMS | 0.308 | 0.269 | 0.245 | 8.28e-07 | 2.35e-05 | 0.000212 | direction and magnitude stable |
| Peak-to-peak | 0.231 | 0.248 | 0.258 | 3.83e-06 | 0.000109 | 0.00061 | direction and magnitude stable |
| Waveform length | 0.253 | 0.302 | 0.346 | 8.39e-08 | 1.55e-06 | 1.92e-05 | direction and magnitude stable |
| Total SampEn | 0.523 | 0.409 | 0.352 | 3.79e-13 | 5.02e-11 | 8.88e-09 | direction stable; magnitude variable |
| Spectral edge 95 | -0.111 | -0.117 | -0.120 | 0.0466 | 0.0269 | 0.0127 | direction and magnitude stable |
| Spectral spread | -0.106 | -0.126 | -0.128 | 0.0466 | 0.00991 | 0.00548 | direction and magnitude stable |
| Hjorth complexity | -0.083 | -0.144 | -0.145 | 0.206 | 0.0151 | 0.00704 | direction stable; magnitude variable |
| Higuchi FD | -0.114 | -0.090 | -0.093 | 0.0133 | 0.0308 | 0.0142 | direction and magnitude stable |

Total SampEn retained the same direction at all durations. This remains only a duration-stability statement about its amplitude-coupled, scale-sensitive state association; it is not evidence of amplitude-independent complexity.

## Direction, magnitude, and inference

- Direction reversals among the original 11: none.
- Direction reversals outside the original 11: fwh_complexity_lempel_ziv_complexity, fwh_complexity_katz_fractal_dimension.
- Original features with material magnitude variation under the locked rule: Median, MAV, Total SampEn, Hjorth complexity.
- Original features with material SE variation (>20% versus frozen 60 s): Waveform length.
- Of spectral edge, spectral spread, Hjorth complexity, and Higuchi FD, 3/4 met the direction-and-magnitude stability rule. Their individual q-values and signs are reported rather than pooled into an amplitude claim.
- Hjorth complexity was the only original feature to lose q<0.05 at one duration (40 s); its 40-s coefficient also moved toward zero, so this is not attributed to precision alone. The comparison table exposes coefficient, SE, and q changes separately.

## Convergence and technical caveats

There were 0 failed models among 68 regenerated feature-duration fits. Woman-only fallback was used for 2/34 features at 40 s and 4/34 at 80 s; boundary status was recorded for 5/34 and 6/34. Exact optimizer warnings and variance estimates are retained in the convergence audit.

Common scaling is the correct primary comparison because it avoids silently comparing different SD units. Raw waveform length and Total SampEn remain duration-dependent by construction, so their raw coefficients should not be interpreted as if they measured an identical per-sample quantity. Event composition also changes modestly with duration because wider windows produce more mixed or boundary-excluded events; this is an unavoidable part of the defined segmentation sensitivity, not a model redesign.

## Reviewer answers and claim recommendations

1. The 40/80-s rerun was technically feasible with the current code and raw recordings.
2. Population changes are shown above; IN tiling accounts for the largest difference.
3. 11/11 original features retained direction at all durations.
4. 10/11 retained BH q<0.05 at all three; 11/11 did so at two or more.
5. The four amplitude/waveform features were directionally robust as a group, with any individual magnitude exception identified above.
6. Total SampEn was directionally stable and remains amplitude-coupled/scale-sensitive.
7. The features most sensitive in magnitude were: Median, MAV, Total SampEn, Hjorth complexity.
8. Spectral/complexity stability was 3/4 by the same rule and should be described feature by feature.
9. Window duration mainly changed magnitude for Median, MAV, Total SampEn, and Hjorth complexity; material SE variation was confined to waveform length. The single original-feature loss of q<0.05 also involved Hjorth coefficient attenuation.
10. The central UC-associated amplitude/time-domain conclusion is not materially dependent on 60 s, although duration-additive individual features may require qualification.

**Claims to retain:** directions and approximate magnitudes that remain stable across all durations, especially the core pattern if the criteria above are met.

**Claims to weaken:** universal FDR support for features whose q-values do not persist, and any weaker spectral/complexity feature with magnitude variation.

**Claims to remove or avoid:** amplitude-independent interpretations of Total SampEn and claims that every individual feature is duration-invariant.

## Recommendation

**The window-length concern is adequately resolved for the central amplitude/time-domain conclusion, with feature-specific magnitude and inferential qualifications.** This analysis stops here and does not alter the manuscript or initiate another reviewer analysis.
