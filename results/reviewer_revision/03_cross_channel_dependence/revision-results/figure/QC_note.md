# Figure quality-control note

- The plotting table contains the seven expected features in the requested order: Total SampEn, mean absolute value, root mean square, peak-to-peak amplitude, waveform length, time median, and time skewness.
- All four channels (EHG9, EHG10, EHG11, EHG12) are present for every feature (28 estimates total).
- Plotted estimates and confidence limits are copied from the `preferred_mixed_model` primary rows of `phase2_final_analysis/corrected_model_results.csv`; the plotting CSV preserves these values.
- The vertical dashed reference line is drawn at coefficient 0 (no UC-versus-IN difference).
- No values were taken from the obsolete cross-channel interaction analysis or from the original per-channel standardisation.
- The PNG was rendered at 600 dpi and the PDF is vector output. Visual inspection confirmed legible labels, distinguishable channels, white background, and correctly positioned confidence intervals.

Caption:

> Cross-channel spatial consistency of UC-versus-IN effects for seven predefined features across EHG9–EHG12. Points show standardised UC-versus-IN estimates and horizontal lines show 95% confidence intervals from the corrected physical-window-aware analysis. The vertical dashed line indicates no state difference. Agreement in effect direction across channels is interpreted as spatial consistency within the same cohort rather than independent replication. Channel-dependent differences in effect magnitude are evaluated separately using the UC-by-channel interaction models.
