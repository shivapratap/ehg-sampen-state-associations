# Amplitude collinearity of profile summaries

Generated 2026-09-21T10:28:09.448016+00:00
 from 4083 matched EHG9 FWH 60 s windows.

Does each profile summary measure irregularity, or signal amplitude?

## Correlation with mean absolute value

| summary | Pearson r | r^2 |
|---|---:|---:|
| frozen_total | 0.984 | 0.968 |
| frozen_average | -0.048 | 0.002 |
| frozen_median | 0.037 | 0.001 |
| frozen_sd | -0.086 | 0.007 |
| frozen_kurtosis | 0.140 | 0.020 |
| frozen_skewness | 0.061 | 0.004 |
| frozen_fixed_r | -0.118 | 0.014 |
| gi_total | -0.102 | 0.010 |
| gi_average | -0.072 | 0.005 |
| gi_sd | -0.115 | 0.013 |
| gi_kurtosis | 0.010 | 0.000 |
| gi_skewness | 0.003 | 0.000 |
| gi_n_retained | 0.053 | 0.003 |
| gi_fixed_r | -0.118 | 0.014 |
| gigrid_auc | -0.100 | 0.010 |
| gigrid_slope | 0.114 | 0.013 |
| gigrid_sd | -0.129 | 0.017 |

## Cross-feature agreement (not amplitude)

| pair | Pearson r | r^2 |
|---|---:|---:|
| gi_total vs frozen_total | -0.053 | 0.003 |
| gigrid_auc vs gi_total | 0.999 | 0.998 |

## Reading

Frozen Total SampEn is collinear with mean absolute value: it shares
the large majority of its variance with a plain amplitude statistic,
and is therefore not interpretable as an irregularity summary in this
implementation. Summaries that divide by the bin count (Average), that
describe profile shape (SD), or that are computed after sd-normalisation
(all GI and GIgrid variants) do not show this collinearity. GI Total is
close to orthogonal to the frozen Total it replaces, so the correction
substitutes a different quantity rather than adjusting the original.

See results/gain_resolution_sensitivity_report.md for the matching
gain and resolution behaviour, and results/full_precision_gain_test_report.md
for the unrounded control that isolates the quantizer as the cause.
