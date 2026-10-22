# Provenance

> **Note on references.** Some comments and file names use the labels "standard pipeline" (the 34-feature analysis) and "GI stages" (the gain-invariant analysis). Both share one tree in this repository (see README.md); GI outputs live under `results/gi/`. References to design notes point to working notes that are not distributed.


This document maps every manuscript claim to the script, input files and
output files that produce it, and records where the "standard" implementation
it is compared against actually comes from. `src/gi_sampen/provenance.py` is
the machine-readable form of the second half of this document.

## 1. The standard finite-resolution implementation: full chain

The brief for this repository asked that the standard implementation be
pinned to the exact upstream commit used for the paper. Auditing the actual
history showed that chain is four layers deep, not one, and the most recent
layer -- the released `sampen-profile` package -- **postdates the analysis
and produced no manuscript number.** Recording that honestly, rather than
collapsing it into a single citation, is the point of this section.

| # | Implementation | Language | Commit / version | Role |
|---|---|---|---|---|
| 1 | `Sample entropy/CHM.m`, original release | MATLAB | `2a1d496c8bcf1e859ae9887d9911d168dd9a17dc` (2020-09-24) | Canonical source of the rounding behaviour. `round(d,3)` / `round(d1,3)`, introduced by this commit. |
| 2 | `amrita-biosignal-feature-engine` | Python | `4dc5025b265f300d7e3375118f7794bfb5fca711` (0.2.0.dev0) | Produced the *cached* `fwh_entropy_*sampen` feature table that the prior study analysed, and that this repository uses as its parity target. Not a runtime dependency here. |
| 3 | `sampen-profile` (github.com/shivapratap/sampen-profile) | Python | `263b72f7b906208241291c0579e8ec103fb52d79` (0.1.0), released 2026-07-20 | Independently released, installable implementation of the same method. **Published after this analysis was run.** Produced no manuscript number. Pinned here as a **test-only** dependency (the `test` extra in `pyproject.toml`) so that `tests/test_standard_parity.py` can assert numerical agreement. |
| 4 | `src/gi_sampen/gi_profile.py::frozen_summaries` | Python | this repository | **Implementation of record for every manuscript number.** A standalone reimplementation, chosen deliberately so parity against layer 2's cached values could be asserted numerically (200 windows, max abs deviation ~1e-15) before anything new was interpreted. |

Every manuscript sentence asserting the rounding behaviour cites layer 1
(the implementation) alongside the parent method paper, never the parent
paper alone -- the published formulation contains no rounding step.

**Why `sampen-profile` is not "the implementation used for the paper".**
Its two commits are dated 2026-07-20; the cached feature table this
repository verifies against (layer 2) was produced on 2026-07-23 by a
different, separately vendored package with its own commit history. There is
no evidence the two share code, only that they implement the same rounding
rule. Declaring `sampen-profile` as *the* pinned dependency, as an early
draft of this reorganisation was asked to do, would have invented provenance.
What is true, and useful, is that it independently reproduces layer 1's
behaviour, which is exactly what a parity test can check without asserting
anything about *how* it was produced.

**Parity result** (see `tests/test_standard_parity.py`, which runs this on
every commit): on three deterministic synthetic signals -- including an
EHG-like `sd ~ 1e-2` cumulative-sum signal -- `frozen_summaries` and
`sampen_profile.sample_entropy_profile` agree to floating-point noise
(worst observed relative deviation 2.0e-13, from median tie-breaking; five of
six summaries agree to <5e-16).

**A specific, documented divergence.** The rounding is specific to the
SampEn-profile implementation: the same release's ApEn-profile code
(`Approximate entropy/CHMforApEn.m`) takes `unique()` over raw distances with
no rounding step. This project's own earlier Python port
(`core/sampleEntropy_Gayathri.py`, not part of this repository) reproduces the
same documentation-versus-implementation divergence: its docstring promises
unrounded r-values while its code rounds to three decimals. Both facts
support treating the rounding as a deliberate, method-specific implementation
choice, not an incidental slip -- see the design notes for
the manuscript's use of this.

## 2. What each implementation in this repository computes

See the module docstring in `src/gi_sampen/gi_profile.py` for the full
mathematical description. Summary:

| Function | Rounding | Standardisation | Role |
|---|---|---|---|
| `frozen_summaries` | 3 decimals, signal units | none | Standard finite-resolution construction; parity target only, not a proposed feature |
| `gi_summaries_rounded` | 3 decimals, units of window SD | divide by window SD | **Gain-Invariant (GI) SampEn profiling.** The gate features. |
| `gi_summaries_grid` | none (fixed dimensionless grid) | divide by window SD | Cross-window-comparable representation (GI-AUC etc.) |
| `full_precision_summaries` (`src/gi_sampen/full_precision.py`) | none | none | Mechanism-isolating control -- exactly gain-invariant with no rounding AND no standardisation; isolates rounding, not amplitude sensitivity per se, as the mechanism. See `results/gi/audits/full_precision_gain_test_report.md`. |
| `fixed_r_sampen` | n/a (fixed tolerance = 0.2 x SD) | implicit | Richman-Moorman comparator, already exactly gain-invariant by construction |

**On terminology.** "Frozen" does not appear in this document or in any
other public-facing file; the manuscript's committed CSV column names
(`frozen_estimate`, `frozen_q_value`, `frozen_total`, etc.) are retained
verbatim for provenance reasons -- renaming them would break the join back to
the prior study's own cached files -- and are documented here instead. Read
`frozen_*` as "standard finite-resolution" throughout.

## 3. Claim -> script -> input(s) -> output(s)

| Manuscript claim | Script | Input(s) | Output(s) |
|---|---|---|---|
| 200-window parity check | `scripts/run_parity_check.py` (`gi_sampen.extract`, `--parity 200`) | upstream window manifest, event-unit table, cached feature table, raw recordings | `results/gi/tables/parity_check.csv` |
| 60-window gain-stress test (Table S16 design + GI) | `scripts/run_gain_stress.py` (`gi_sampen.gain_stress`, `--extra-scales`) | upstream event-unit table, raw recordings | `results/gi/tables/invariance_summary.csv`, `work/invariance_window_results.csv` |
| 12-window full-precision control | `scripts/run_full_precision_control.py` (`gi_sampen.full_precision`) | upstream event-unit table, raw recordings | `results/gi/tables/full_precision_gain_test.csv`, `results/gi/audits/full_precision_gain_test_report.md` |
| AvgSampEn gain sensitivity | `scripts/report_avg_sampen_gain.py` (`gi_sampen.report_avg_gain`) | `work/invariance_window_results.csv` (from gain-stress) | `results/gi/tables/avg_sampen_gain_summary.csv`, `results/gi/tables/avg_sampen_gain_test.csv` |
| Resolution sensitivity | `scripts/run_gain_resolution.py` + `scripts/report_gain_resolution.py` | raw recordings, upstream event-unit table | `results/gi/tables/gain_resolution_sensitivity*.csv`, `results/gi/audits/gain_resolution_sensitivity_report.md` |
| GI resolution ablation | `scripts/run_resolution_ablation.py` + `scripts/report_resolution_ablation.py` | `work/gi_features_60s_EHG9.csv`, `results/gi/tables/gi_model_results.csv` | `work/resolution_ablation_features.csv`, `results/gi/tables/resolution_ablation_summary.csv` |
| Dimensionless fixed-grid check | `scripts/run_select_grid.py` (`gi_sampen.fixed_grid`) | upstream event-unit table (baseline windows), raw recordings | `results/gi/tables/grid_definition.json`, `results/gi/tables/grid_coverage_baseline.csv` |
| EHG9 primary statistical analysis | `scripts/run_statistics.py` (`gi_sampen.statistics`) | `work/gi_features_60s_EHG9.csv`, upstream cached feature table | `results/gi/tables/gi_model_results.csv` |
| Go/no-go decision | `scripts/run_gate_report.py` (`gi_sampen.gate`) | `results/gi/tables/invariance_summary.csv`, `results/gi/tables/gi_model_results.csv`, `results/gi/tables/grid_definition.json` | `results/gi/tables/gate_summary.csv` |
| EHG9-EHG12 cross-channel analysis (**Table 1**) | `scripts/run_cross_channel_analysis.py` + `scripts/report_cross_channel.py` | upstream event-unit tables (4 channels), upstream per-channel model results, raw recordings | `work/cross_channel_gi_features.csv`, `results/gi/tables/cross_channel_replication_summary.csv` |
| Amplitude collinearity | `scripts/run_amplitude_collinearity.py` (`gi_sampen.amplitude`) | upstream cached feature table, `work/gi_features_60s_EHG9.csv` | `results/gi/tables/amplitude_collinearity.csv`, `results/gi/audits/amplitude_collinearity_report.md` |
| **Figure 1** | `scripts/make_figure1.py` (`gi_sampen.figure1`) | `work/gi_features_60s_EHG9.csv`, `results/gi/tables/cross_channel_replication_summary.csv`, `work/invariance_window_results.csv`, `results/gi/tables/full_precision_gain_test.csv`, `results/gi/tables/avg_sampen_gain_summary.csv`, upstream cached feature table | `results/gi/figures/Figure1_support_mechanism.{pdf,png,svg}` + 3 source CSVs |
| **Table 1 or its source CSV** | (same as cross-channel analysis, above) | | `results/gi/tables/cross_channel_replication_summary.csv` |

`scripts/reproduce_paper.py` runs every row above in dependency order; see
`docs/reproducibility.md`.

## 4. What is NOT claimed

Carried over verbatim from the facts pack (the design notes,
section 12) because these boundaries are binding on this document too:

- Not a claim of priority for gain-invariant / amplitude-robust entropy
  estimation, which is well established (Richman-Moorman fixed-r SampEn is
  already exactly gain-invariant, and is carried through here as a
  comparator, not a novel contribution).
- Not a claim that SampEn profiling is *inherently* gain-sensitive -- only
  that the released finite-resolution implementation, as specified above, is.
- Not a claim that GI-SampEn removes all amplitude effects.
- The invariance claim is scoped to **constant positive multiplicative gain
  only**.
