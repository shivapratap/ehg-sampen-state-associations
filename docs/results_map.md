# Where each manuscript element comes from

Display items of the Results in Engineering manuscript, the script that produces their numbers, and the file to look at.
"derived" = runs from `data/derived` (no raw data). "raw" = needs the PhysioNet recordings.

| Manuscript element | Produced by | Look at |
|---|---|---|
| Fig. 1 study design and feature pipeline | `scripts/make_figure_study_design.py` (schematic; no data read; the two occupied-bin counts are those of Supplementary Table S3, Panel D) | `results/manuscript/figures/F1_study_design.pdf` |
| Table 1 cohort / population (4,283 candidate → 4,083 analysed) | `build_figure_inputs.py` → `make_manuscript_figures.py` (T1) ; counts from `src/windowing.py` (raw) | `results/manuscript/figures/T1_cohort.tex`, `results/manuscript/cohort_reconciliation.csv`, `work/windowing/cohort_summary_60s.md` |
| Fig. 2 all 34 EHG9 UC–IN coefficients | `src/models.py` (derived) → `make_manuscript_figures.py` (F2) | `results/primary/mixed_model_results_60s_fwh_EHG9.csv`, `results/manuscript/figures/F2_*.pdf` |
| Table 2 robustness (temporal dependence, woman-level, balanced IN, window length, overlap, residuals) | `src/robustness.py` (derived; `--include-heavy` for balanced IN and 40/80 s) | `results/robustness/{dependence,patient_level,overlap,residuals,window_length}/`; frozen snapshot `results/reviewer_revision/` |
| Fig. 3 gain response (60 windows) | `manuscript_audits/export_waveform_package.py` (raw) → `run_D_gain.py` (from the packaged windows; needs the feature engine) → `make_manuscript_figures.py` (F3) | `results/manuscript/audits/D/gain_audit_window_results.csv`, `results/manuscript/figures/F3_*.pdf` |
| Table 3 RMS adjustment of entropy features | `src/entropy_analysis.py` (derived) | `results/entropy/corrected_amplitude_adjusted_model_results.csv`, `rms_adjusted_entropy_results.csv` |
| Table 4 matched original vs GI associations | `scripts/run_statistics.py`, `run_corrected_family.py`, `build_manuscript_tables.py` | `results/manuscript/tables/T2_matched_original_corrected_UC_FM.csv`, `results/gi/tables/gi_model_results.csv` |
| Table 5 numerical validation controls | `reproduce_paper.py` stages (parity, gain stress, full-precision control, resolution ablation, fixed grid; raw) + `run_exact_gate.py`, `run_D_gain1_gate.py`, `run_C3_grid.py` | `results/gi/tables/*`, `results/manuscript/audits/{D,C3}/`, `results/manuscript/tables/S2*`, `S3*` |
| Fig. 4 GI Total and GI profile SD, UC–IN and FM–IN, EHG9–EHG12 (16-test BH family) | `scripts/run_corrected_family.py` (derived) → `scripts/make_figure_gi_effects.py` | `results/manuscript/corrected_family/corrected_cross_channel_all_16.csv`, `results/manuscript/figures/F5_GI_state_associations.pdf` |
| Table 6 conditional (log-RMS) analysis | `scripts/run_conditional_amplitude.py` (`extract` raw; the rest derived) | `results/conditional/refit_gi_total_rms_results.csv`, `within_recording_results.txt`, `shape_models_output.txt`, `lowo_results.csv` |
| Sensitivity to annotation exclusions / whole-woman influence (Table S5) | `manuscript_audits/run_B_sensitivity.py`; `run_conditional_amplitude.py lowo` | `results/manuscript/audits/B/`, `results/manuscript/tables/S4_annotation_sensitivity_actual.csv`, `results/conditional/lowo_results.csv` |
| Supp. Table S1 features | `src/features.py`, `src/models.py` registry | `results/primary/confirmatory_feature_registry_60s_fwh.csv` |
| Supp. Table S2 model dependence | `build_figure_inputs.py` | `results/manuscript/tables/S1_model_dependence_actual.csv` |
| Supp. Table S3 validation controls | as Table 5 | `results/manuscript/tables/S2*.csv`, `S3*.csv` |
| Supp. Table S4 complete corrected associations | `run_corrected_family.py` | `results/manuscript/corrected_family/` |
| Supp. Table S6 analysis and multiplicity map | `build_manuscript_tables.py` | `results/manuscript/tables/S5_GI_family_roles.csv` |
| Compact Supplement (Tables S1–S7, as submitted); complete listings with old numbering are in `results/manuscript/supplement/full/` | `scripts/build_compact_supplement.py` (18 checks) | `results/manuscript/supplement/Supplement_compact.tex`, `Supplement_compact.pdf`, `Supplement_compact_checks.txt` |
| Complete listings: Supp. Tables S7 (robustness panels A–D, lag-1) and S8 (RMS adjustment, conditional) | `scripts/build_supplement_tables.py`; lag-1 from `run_lag1_diagnostic.py` | `results/manuscript/supplement/Supplement_S7_S8.tex` and `…_checks.txt` (36 checks) |
| Supp. Fig. S1 independent grid agreement | `run_C3_grid.py` → `make_manuscript_figures.py` (F4) | `results/manuscript/audits/C3/independent_grid_agreement.csv`, `F4_*.pdf` |
| Statistical-software cross-check (Methods) | `scripts/verify_statsmodels.py` | `results/conditional/statsmodels_verification.{csv,md}` |

The `results/reviewer_revision/` folder is a frozen copy of the outputs behind Supplement S7/S8 (they were produced by the standalone
scripts that `src/components/robustness_components/` now wraps); the S7/S8 builder reads it so the supplement can be rebuilt without refitting anything.
