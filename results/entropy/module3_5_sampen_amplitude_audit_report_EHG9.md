# Stage 07 (EHG9): SampEn-Amplitude Dependence and Interpretability Audit

## Objective
Post-hoc assessment of whether SampEn-derived features are amplitude-coupled
on channel EHG9. This audit does not change Stage 06 models or the
frozen validated feature set.

## Sources and methods
- Source: `/Users/shivangayathri/Work/EHGProcessing/ehg-sampen-state-associations/work/primary/event_unit_analysis_60s_EHG9.csv` (4,083 analysis windows; 42 women; 112 recordings).
- SampEn features: Total SampEn, Average SampEn, fixed-r SampEn, Median SampEn,
  SD SampEn, Kurtosis SampEn, and Skewness SampEn.
- Amplitude features: MAV, RMS, peak-to-peak, and waveform length.
- Pearson and Spearman correlations were computed at pooled/window, recording, patient, and patient effect levels.
- Primary effect: contraction minus baseline. Fetal movement minus baseline is secondary and was not pooled with contraction.

## Key results
See `module3_5_sampen_amplitude_summary_table_EHG9.csv` for the complete matrix. Spearman correlations are emphasized. Strong coupling is defined as |rho| >= 0.75; moderate as 0.50-<0.75; weak as <0.50.

The pooled window-level results are descriptive only because windows are nested and overlapping. Patient-level effects are more interpretable but have smaller sample sizes.

## Interpretation and limitations
Total SampEn should be described as an amplitude-coupled nonlinear/entropy-derived summary when strong coupling is observed, not as an amplitude-independent complexity marker. Average SampEn is retained as an exploratory comparator; any weaker or inconsistent association should not be interpreted as validation of a replacement feature. Correlation is not causation, and repeated-window precision is not independent evidence.

## Recommended Paper 1 wording (channel EHG9)
"Post-hoc analyses indicated that SampEn-derived measures, particularly Total SampEn, showed [substantial/moderate/weak -- see summary table] association with amplitude-domain measures on channel EHG9. We therefore interpret Total SampEn conservatively as an amplitude-coupled nonlinear entropy summary rather than an amplitude-independent complexity biomarker. Average SampEn was examined as an exploratory comparator and was not used to revise the prespecified feature set."

## Stage 06 status
Stage 06 should remain frozen. Stage 07 is an interpretability audit and provides no basis for revising the primary inferential model or validated feature set.
