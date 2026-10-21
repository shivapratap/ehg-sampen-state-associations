# Known issues and decisions

Read this before trusting a green run.

## What was and was not executed when this repository was assembled

Executed in the build environment (numpy/pandas/scipy/matplotlib, no statsmodels, no wfdb, no raw data) and compared with saved values:

| Item | Result |
|---|---|
| `run_conditional_amplitude.py refit / within / shape` | reproduces saved tables (coefficients ≤ 4e-8, SEs ≤ 2e-7); shape output identical |
| `run_conditional_amplitude.py lowo` | slow (about 3 min per woman in the build environment). The first 19 of 42 women were refitted and agree with the saved `lowo_results.csv` at the 3 decimals printed; the rest was not waited for. The saved file holds all 42 |
| `manuscript_audits/run_C3_grid.py` | reproduces `independent_grid_agreement.csv` to 1.8e-15 |
| `build_supplement_tables.py` | 36/36 checks PASS; the LaTeX is byte-identical to the one used in the manuscript |
| `build_figure_inputs.py`, `build_manuscript_tables.py`, `make_manuscript_figures.py` | run; outputs equal the saved ones (differences ≤ 4e-15 or only in a source-path text column) |
| `restore_derived.py` | extracts 62 files, all SHA-256 verified |
| `check_checkpoints.py` | 15/15 |
| `tests/` | 69 tests pass and 1 skips (needs the optional `sampen-profile` package). They were run with a minimal stand-in for pytest because pytest could not be installed in the build environment, so run `pytest` yourself once |

Not executed here (static syntax check only), because they need statsmodels, wfdb, the raw recordings or the feature-engine package:
`src/models.py`, `robustness`, `cross_channel`, `entropy_analysis`, `fm_analysis`, `tables`, `figures`, `run_statistics`, `run_corrected_family`, `run_lag1_diagnostic`,
`run_exact_gate`, `run_B_sensitivity`, `run_D_gain*`, and every **raw** stage. They are the original scripts with paths changed to be repository-relative
(`run_corrected_family.py` and `run_lag1_diagnostic.py` are ports of one-off scripts). Their *saved outputs* are in `results/` and are what the manuscript reports.
The statsmodels cross-check (`verify_statsmodels.py`) was run on the authors' machine (statsmodels 0.14.5, Python 3.13.9): all differences within tolerance.

**Recommended first check on your machine:** commit the repository as shipped, run `python scripts/run_all.py`, then `git diff --stat results/`.
Anything that moved is a result the code did not reproduce (tiny last-digit differences across platforms are expected for floating-point output).

## Decisions

- **One REML fitter.** The conditional analyses use `src/reml.py` (variance components bounded below by exp(−12)). Earlier hand-run copies used −30 in one place and −12 in another; −30 makes one covariance block numerically non-positive-definite in the shape models. Coefficients are the same either way (checked: ≤ 4e-8).
- **Robustness outputs are refreshed through a scratch directory.** Two of the robustness components refuse to overwrite existing files; `src/robustness.py` runs them into `work/robustness/_fresh/` and copies the result over `results/robustness/`.
- **Stale outputs left out.** The one-off driver that produced the 16-test family also wrote "blocked" placeholder folders and a stale edge-case report; those were not carried over. The bootstrap pilot output (`A`) and the early `C4` edge-case report are omitted. The gain audit (`D`) keeps only the files that match the final run.
- **Absolute paths** inside saved provenance files (`input_hashes.csv`, run metadata, some reports) were replaced by `<repo>`; the recorded hashes are unchanged.
- **Timing figures** in the README are rough estimates, not measurements, apart from the GI extraction times quoted in `scripts/reproduce_paper.py`.

## Things to be aware of

- `results/primary` etc. are overwritten in place by their stages; use git to compare.
- `features` (hours) needs the pinned `amrita-biosignal-feature-engine` commit; restoring `data/derived` avoids it.
- Balanced-IN (1,000 resamples) and the 40/80-s window analyses are not part of the default run (`--heavy`). Their small summary files are in `results/reviewer_revision/`; the 34×1,000 iteration table (7.7 MB) and the convergence audit (22 MB) were not copied; the heavy stage regenerates them under `work/robustness/balanced_in/`.
- Figure 1 (study design) is artwork, not code output.
- `docs/methods.md` and `docs/provenance.md` are carried over from the earlier GI repository; they still cite that repository's internal "facts pack", which is not distributed.
