# Sample-entropy gain sensitivity and state associations in electrohysterography (EHG9)

Code, derived tables and saved results for a study of uterine-contraction (UC), fetal-movement (FM) and inactivity (IN) windows in the public Icelandic 16-electrode EHG database. The study:

1. cuts 60-s windows around annotated events (and from guard-cleared quiet intervals), filters the whole recording at 0.3–1.0 Hz (the "FWH" band), decimates to 20 Hz, and computes a 34-feature family per window (amplitude, frequency, entropy, Hjorth, fractal, complexity) on four median-axis channels, EHG9 (primary) to EHG12;
2. tests UC–IN and FM–IN differences with woman- and recording-nested mixed models (REML, Benjamini–Hochberg FDR), plus robustness analyses (temporal dependence, woman-level analysis, balanced IN, window length, overlap, residual diagnostics);
3. shows that finite-resolution sample-entropy-profile summaries (Total SampEn) are sensitive to signal gain, and evaluates a gain-invariant alternative (**GI-SampEn**: divide each window by its SD before profiling), including how its state association changes when conditioned on amplitude (log RMS).

Raw recordings are **not** in this repository: get them from PhysioNet (<https://physionet.org/content/ehgdb/>, see [data/README.md](data/README.md)). A 40 MB zip of derived tables is included, so everything downstream of feature extraction can be rerun without them.

This repository combines the standard 34-feature analysis, the gain-invariant (GI) sample-entropy stages and the manuscript-stage scripts in one tree. A glossary note: "frozen" in file and column names (for example `frozen_total`) denotes the original (fixed-grid rounding) definition used as the comparator in the paper. [docs/known_issues.md](docs/known_issues.md) says exactly what has and has not been re-run since the clean-up. Please read it.

---

## 1. Quick start (no raw data needed)

```bash
git clone <this repository> && cd ehg-sampen-state-associations
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt            # numpy, pandas, scipy, statsmodels, matplotlib, ...
python scripts/restore_derived.py          # unpack data/derived/ehg_derived_inputs.zip -> work/ (checksums verified)
python scripts/check_environment.py        # what can this machine run?
python scripts/check_checkpoints.py        # 15 key numbers on disk vs the manuscript: expect 15/15
pytest                                     # unit tests (no data needed)
python scripts/run_all.py --list           # every stage, what it needs, what it writes
python scripts/run_all.py                  # re-run everything that does not need raw data
git diff --stat results/                   # anything that changed is a number the code did not reproduce
```

The repository ships the results of the authors' runs in `results/`. Re-running overwrites them in place, so the cleanest reproducibility test is to commit first, run, then look at `git diff`. Small last-digit differences across platforms are expected.

Python 3.13 was used (3.10+ should work). `requirements.txt` pins the versions that produced the numbers.

## 2. Reproducing from the raw recordings (optional)

```bash
pip install -r requirements-raw.txt        # adds wfdb and the pinned feature-engine commit
export EHG_RAW_DIR=/path/to/ehgdb          # folder with the .dat/.hea/.atr files, RECORDS, SHA256SUMS.txt
python scripts/run_all.py --with-raw       # raw stages first, then everything else
```

Raw stages, in order: `data` (inventory and checksums), `annotations`, `windowing`, `features` (34 features × 4 channels; the slow one), `feature_qc`, `gi_pipeline` (GI feature extraction, gain stress tests, resolution ablations, cross-channel extraction; wraps `scripts/reproduce_paper.py`), `conditional_extract` (per-window GI/frozen Total, RMS, crest factor, kurtosis for EHG9), `waveforms` (packs the 60 audited windows). The raw stages regenerate the contents of `data/derived/` and the window-level conditional table.

## 3. Repository layout

```
README.md                    this file
requirements.txt             core environment (pinned)
requirements-raw.txt         + wfdb and feature engine, for the raw stages
CITATION.cff, LICENSE        provisional citation (TODOs inside), MIT
conftest.py                  makes src/ importable for pytest
configs/
  analysis.yaml              paths (relative) and analysis constants for the GI stages
  figures.json               specification of the manuscript figure/table renderer
data/
  README.md                  where to get the raw data, what the derived zip contains
  raw/                       (empty, git-ignored) put or link the PhysioNet files here
  derived/ehg_derived_inputs.zip   derived tables (40 MB) and MANIFEST.sha256
src/                         library code, one module per pipeline stage
  config.py                  all fixed protocol constants and directory locations
  data.py annotations.py windowing.py preprocessing.py features.py feature_qc.py
                             raw -> windows -> features (the standard 34-feature pipeline)
  models.py                  primary mixed-effects models (34 features x 4 channels), FDR
  robustness.py              runs src/components/robustness_components/*
  cross_channel.py           UC-IN across the four channels, shared-window dependence
  entropy_analysis.py        entropy-amplitude coupling, scale experiment, RMS adjustment
  fm_analysis.py tables.py figures.py   FM results, summary tables, original figures
  reml.py                    small REML fitter for the conditional (log-RMS) analyses
  components/                implementation modules wrapped by the files above
  gi_sampen/                 GI-SampEn library: profile summaries (gi_profile.py), paths,
                             gain stress, full-precision control, fixed grid, statistics, reports
scripts/
  run_all.py                 orchestrator (use --list)
  check_environment.py       packages / raw data / derived data present?
  check_checkpoints.py       key numbers on disk vs the manuscript
  restore_derived.py  pack_derived.py   unpack / (maintainers) rebuild the derived zip
  reproduce_paper.py         the GI pipeline end to end (wrapped by run_all 'gi_pipeline')
  run_*.py report_*.py make_figure1.py   individual GI stages
  run_corrected_family.py    16-test BH family for the GI cross-channel result
  run_lag1_diagnostic.py     residual lag-1 autocorrelation (Supplement S6A)
  run_conditional_amplitude.py   log-RMS, within-recording, shape, leave-one-woman-out (Supplement S7)
  verify_statsmodels.py      statsmodels cross-check of the custom REML fitter
  manuscript_audits/         waveform-based audits: gain audit, gain-1 gate, grid check,
                             exact-model gate, annotation sensitivity, waveform packaging
  build_figure_inputs.py  build_manuscript_tables.py  make_manuscript_figures.py
                             manuscript display items from saved CSVs (fit nothing)
  build_supplement_tables.py full-listing Supplement S7/S8 LaTeX and its 36-number check list (complete listings in results/manuscript/supplement/full/)
tests/                       gain invariance, unrounded equivariance, parity, result fingerprints
docs/                        methods.md, provenance.md, results_map.md, known_issues.md
work/                        (git-ignored) intermediate tables; restored from data/derived
results/                     saved outputs, committed
  primary/ robustness/ cross_channel/ entropy/ fm/ tables/ figures/   standard analysis
  gi/{tables,audits,figures}/   GI-SampEn outputs
  conditional/               log-RMS / shape / leave-one-woman-out outputs
  manuscript/                corrected family, audits, adapters, tables, figures, supplement
  reviewer_revision/         frozen robustness outputs behind Supplement S6/S7
```

## 4. Stages

`python scripts/run_all.py --list` prints this table. Times are rough estimates, not measurements.

| Stage | Group | Needs | Writes |
|---|---|---|---|
| data, annotations, windowing | raw | raw data, wfdb | `work/audit`, `work/annotations`, `work/windows`, `work/windowing` |
| features | raw | raw data, feature engine; hours | `work/features/window_features_60s_EHG{9..12}.csv` |
| feature_qc | raw | | `work/feature_qc/` |
| gi_pipeline | raw | raw data; ~3 h | `work/gi_*.csv`, `work/cross_channel_gi_features.csv`, `results/gi/` |
| conditional_extract | raw | raw data | `results/conditional/windows_gi_rms_frozen_full.csv`, `window_shape_descriptors.csv` |
| waveforms | raw | raw data | `results/manuscript/audits/waveforms/` |
| gain_audit, gain1_gate | engine | the 60 packaged windows + feature engine | `results/manuscript/audits/D/` |
| primary | derived | statsmodels | `work/primary/`, `results/primary/` |
| robustness | derived | statsmodels | `results/robustness/*` |
| cross_channel, entropy, fm | derived | statsmodels | `results/cross_channel`, `results/entropy`, `results/fm` |
| gi_statistics, gi_reports | derived | statsmodels | `results/gi/tables/gi_model_results.csv`, `gate_summary.csv` |
| corrected_family | derived | statsmodels | `results/manuscript/corrected_family/` |
| lag1 | derived | statsmodels | `results/robustness/lag1/` |
| conditional | derived | scipy only; roughly 2 h in the build environment (about 3 min per woman), almost all of it the leave-one-woman-out refits (`refit`, `within`, `shape` take a few minutes) | `results/conditional/` |
| grid_check | derived | scipy only | `results/manuscript/audits/C3/` |
| exact_gate, annotation_sensitivity | derived | statsmodels | `results/manuscript/audits/` |
| verify_statsmodels | derived | statsmodels | `results/conditional/statsmodels_verification.*` |
| figure_inputs, manuscript_tables, manuscript_figures | derived | pandas, matplotlib | `results/manuscript/{adapters,tables,figures}` |
| tables, figures | derived | | `results/tables`, `results/figures` |
| supplement | derived | stdlib | `results/manuscript/supplement/` |
| robustness_heavy (`--heavy`) | derived | statsmodels; hours | balanced-IN (1,000 resamples), 40/80-s windows |

Run one stage: `python scripts/run_all.py --only lag1`. Resume: `--from robustness`. Each stage is also an ordinary script you can run directly (`python scripts/run_lag1_diagnostic.py`).

**Stages that need only numpy/pandas/scipy** (so they run anywhere): `conditional`, `grid_check`, `supplement`, the figure/table builders, and `check_checkpoints.py`.

## 5. What a correct run looks like

Fixed checkpoints (all verified by `scripts/check_checkpoints.py`):

| Quantity | Value |
|---|---|
| Candidate / analysed 60-s windows | 4,283 / 4,083 (IN 3,358, UC 333, FM 392) |
| Women / recordings analysed | 42 / 112 |
| EHG9 UC-vs-IN features FDR-supported (34-feature family) | 11 of 34 |
| GI Total, EHG9, UC − IN (standardised) | −0.113, q = 0.021 in the 16-test family |
| Total SampEn UC − IN, adjusting for log RMS | frozen 0.409 → 0.020; GI −0.113 → +0.275 |
| Leave-one-woman-out (42 refits), GI Total | unadjusted −0.134…−0.088, adjusted +0.244…+0.305 |
| Lag-1 residual autocorrelation, original-11 features | 25 of 55 with \|r\| ≥ 0.30, max 0.560 |
| Custom REML vs statsmodels | coefficients ≤ 2.7e-6, SEs ≤ 2.8e-3 |
| Supplement S7/S8 number checks | 36 of 36 PASS |
| Compact Supplement (S1–S7) checks | 18 of 18 PASS (`scripts/build_compact_supplement.py`) |

`tests/test_reproducibility.py` also fingerprints the GI result files in `results/gi/`.

## 6. Method notes worth knowing

- **Filtering before slicing.** The full recording is band-passed (4th-order zero-phase Butterworth, SOS) and then windows are cut and decimated 200 → 20 Hz, so windows have no independent filter-edge transients. `src/preprocessing.py` is the single definition used everywhere.
- **GI-SampEn.** `src/gi_sampen/gi_profile.py`: `frozen_summaries` is the standard finite-resolution profile (distances rounded to 3 decimals in signal units); `gi_summaries_rounded` divides the window by its SD first; `gi_summaries_grid` evaluates the standardised profile on a fixed dimensionless grid. See [docs/methods.md](docs/methods.md).
- **Mixed models.** `src/models.py` uses statsmodels `MixedLM` (woman random intercept + recording variance component, REML, optimiser fallbacks, woman-only fallback). `src/reml.py` is a small independent block-diagonal REML implementation used for the conditional analyses and the leave-one-woman-out refits (42 × 4 fits are slow in statsmodels); its agreement with statsmodels is recorded in `results/conditional/statsmodels_verification.md`.
- **Multiplicity.** BH within defined families: 34 features per contrast for the standard analysis; for the GI result, the declared 16-test family (4 channels × 2 summaries × 2 contrasts). GI summaries were chosen after development; the manuscript says so.
- **Randomness.** The analyses are deterministic. Balanced-IN subsampling (heavy stage) is seeded inside its component.

## 7. Data, ethics, citation

The recordings are public (PhysioNet; see `data/README.md` for the licence and the citation to use). The derived tables contain window-level features and the database's own record identifiers, no signal samples. If you use this code, cite the manuscript (provisional entry in `CITATION.cff`, still to be completed) and the database.

Code is MIT-licensed (see `LICENSE`).
