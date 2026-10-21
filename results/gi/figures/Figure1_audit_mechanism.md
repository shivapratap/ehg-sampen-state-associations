# Figure 1 audit -- mechanism layout

Generated 2026-09-22T05:03:10+00:00 by `src/make_figure1.py --layout mechanism`.

## Source files

- Panel (a) standard/GI support: `Paper2/work/invariance_window_results.csv`
- Panel (a) unrounded support, panel (b) unrounded deviation: `Paper2/results/full_precision_gain_test.csv`
- Panel (b) standard/GI deviation: `Paper2/results/avg_sampen_gain_summary.csv`
- Panel (c): `Paper2/results/cross_channel_replication_summary.csv`

## Columns used

- gain: `scale_factor` (window table) / `gain` (full-precision table)
- support size: `frozen_n_retained`, `gi_n_retained`, `fp_n_retained`
- deviation: `median_relative_difference_percent` (standard/GI), `fp_total` (unrounded)
- Panel (c): `estimate`, `frozen_estimate`, `q_value`, `frozen_q_value` (legacy source-column names retained internally)

## Window sets behind each series

- standard and GI: n = 60 windows (of 60 in the gain-stress set)
- unrounded control: n = 12 windows
- the unrounded windows are a strict subset of the gain-stress windows: True
- common-subset mode: False

The series therefore do NOT rest on identical window sets. Restricting standard and GI to the same 12 windows (`--common-subset`) shifts the support medians by under 3%, the fitted standard slope from 0.976 to 0.981, and the standard Total deviation at 7.1x from 630.6% to 627.4%; the GI and unrounded series are exactly zero either way. The comparison is therefore insensitive to the choice, but the n for each series is stated in the legend and should be stated in the caption.

- panel (b) standard/GI source: Paper2/results/avg_sampen_gain_summary.csv

## Panel (a) -- retained tolerance-support size, median over windows

| gain | standard | standard / 1x | GI | unrounded |
|---:|---:|---:|---:|---:|
| 0.25 | 22 | 0.26 | 5,258 | 680,924 |
| 0.5 | 42 | 0.51 | 5,258 | 680,924 |
| 1 | 84 | 1.00 | 5,258 | 680,924 |
| 1.3 | 109 | 1.30 | 5,258 | 680,924 |
| 2 | 166 | 1.98 | 5,258 | 680,924 |
| 4 | 329 | 3.92 | 5,258 | 680,924 |
| 7.1 | 571 | 6.80 | 5,258 | 680,924 |

The standard support tracks gain almost exactly one-for-one; the GI and unrounded supports are constant. Only medians over the existing gain-stress windows are computed -- no model, no estimation, no new statistic.

## Panel (b) -- median |deviation| from the unit-gain value

- standard and GI series read verbatim from `Paper2/results/avg_sampen_gain_summary.csv` (absolute value taken for display).
- unrounded series reduced from the full-precision per-window table and checked against its published verdict of exactly zero; the script stops if it exceeds 1e-6%. Observed worst case: 0%.
- y axis is symlog (linear below 1%, log above) so an exact zero and a 600% deviation can share one axis.

## Panel (c) values (read verbatim; no model or FDR was recomputed)

| feature | channel | GI estimate | GI q | standard estimate | standard q |
|---|---|---:|---:|---:|---:|
| gi_sd | EHG9 | -0.2460 | 1.88e-07 | -0.2087 | 0.000194 |
| gi_sd | EHG10 | -0.3452 | 1.63e-12 | -0.2927 | 1.14e-08 |
| gi_sd | EHG11 | -0.3323 | 3.59e-11 | -0.3086 | 5.28e-09 |
| gi_sd | EHG12 | -0.2787 | 2.52e-08 | -0.2626 | 1.96e-06 |
| gi_total | EHG9 | -0.2003 | 7.16e-06 | +0.0777 | 0.459 |
| gi_total | EHG10 | -0.2736 | 1.17e-09 | +0.0193 | 0.962 |
| gi_total | EHG11 | -0.2937 | 4.03e-10 | +0.0739 | 0.487 |
| gi_total | EHG12 | -0.2545 | 4.76e-08 | +0.0668 | 0.554 |

q-values are BH-FDR within the 16-test cross-channel family (4 channels x 2 features x 2 contrasts), as computed upstream. The four channels are spatially corresponding measurements from the same women and the same physical windows: cross-channel consistency, not independent replication.

## Display-only transformations

- Panels (a) and (b): log x; log y in (a), symlog y in (b).
- Absolute value taken in (b) so standard Total's sign change across gain does not split the series; signed values are in the source CSV.
- No transformation of any kind is applied to panel (c) values.

## Outputs

- `Paper2/results/figures/Figure1_support_mechanism.pdf`
- `Paper2/results/figures/Figure1_support_mechanism.png`
- `Paper2/results/figures/Figure1_support_mechanism.svg`
- `Paper2/results/figures/Figure1_support_source.csv`
- `Paper2/results/figures/Figure1_deviation_source.csv`
- `Paper2/results/figures/Figure1_panelB_source.csv`
