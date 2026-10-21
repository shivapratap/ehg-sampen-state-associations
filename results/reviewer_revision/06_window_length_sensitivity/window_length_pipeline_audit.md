# Window-length pipeline audit

## Scope and conclusion

The current Paper 1 code supports arbitrary integer window durations without changing the frozen Stage 00–09 products. Stage 02 and Stage 04 expose `--window-size` and isolated output-directory arguments. Stage 03 is a function library consumed by Stage 04. Stage 05 exposes a duration argument but hard-codes its input filename pattern to `60s`; Stage 06 similarly hard-codes `60s` output names. The revision analysis therefore reuses their audited functions from a separate wrapper rather than editing either frozen stage.

The current Stage 02 code regenerated the 40-s and 80-s manifests in the revision directory. Both regenerated files match the corresponding legacy Module 1 manifests exactly, including row order, values, columns, and SHA-256 hashes. Legacy manifests are treated only as validation references. No legacy 40-s or 80-s feature table exists, so feature values must be regenerated with the current Stage 03/04 implementation.

## Event-centred windows

- Contraction, possible-contraction, and fetal-movement annotations are primary window anchors.
- A window is centred on the annotation time with start `event_time - duration/2` and end `event_time + duration/2`.
- UC and FM therefore use the same centring rule. Possible contractions are constructed for audit completeness but excluded from the primary UC/FM/IN analysis.
- The usable recording interval is `[15 s, recording_duration - 15 s]`. An event window extending outside it is rejected as `outside_usable_range_after_edge_trim`; event windows are never truncated or shifted.
- Sample indices are rounded after multiplying the exact window boundaries by the recording sampling frequency.

## Overlap handling

- All annotations are examined when assigning nearest-event and overlap diagnostics.
- If another primary event lies within half a window of the anchor, the anchor window contains that event. A different primary class produces `mixed`; the same primary class produces `same_class_cluster`.
- Any pair of event-centred windows with annotation separation shorter than the window duration overlaps. When there is no contained primary event, overlap of at most 1 s is labelled `boundary_overlap`; larger overlap is `substantive_overlap`.
- Participant movement, equipment manipulation, and position change do not anchor event windows but set the movement-artifact diagnostic when contained.
- Mixed windows and possible-contraction windows remain in the manifest but are excluded by the frozen Stage 06 event-unit rules.

## Inactivity construction and baseline guards

- Every contraction, possible contraction, fetal movement, participant movement, equipment manipulation, and position change annotation creates a ±120-s forbidden interval.
- Overlapping forbidden intervals are merged. Their complement within the 15-s edge-trim limits defines clean baseline gaps.
- Each clean gap is tiled from its left boundary with non-overlapping windows of the requested duration. The number is `floor(gap_duration/window_duration)`; an unused remainder shorter than one window is recorded as a rejected clean gap.
- Baseline windows are constructed only for recordings with an annotation file. No baseline is inferred from unannotated recordings.
- Baseline windows are not truncated. Guard clearance and tight-clearance status (<30 s beyond the guard boundary) are retained as diagnostics.

## Preprocessing order

- Raw EHG9 is validated as finite one-dimensional data sampled at 200 Hz.
- The complete recording is filtered before any window is sliced.
- Filtering uses a fourth-order zero-phase Butterworth SOS bandpass for FWH (0.3–1.0 Hz) and Broad (0.1–3.0 Hz).
- The already filtered recording is sliced at manifest sample bounds and decimated by taking every tenth sample, producing 20-Hz analysis data.
- Therefore 40-s, 60-s, and 80-s windows contain 800, 1,200, and 1,600 filtered samples, respectively. No isolated-window filtering or duration-specific filter retuning occurs.

## Feature extraction assumptions

- The same 87-column Stage 04 registry is computed for every duration; the primary family is the frozen set of 34 non-redundant FWH features.
- Welch PSD settings remain a 10-s Hann segment with 5-s overlap, constant detrending, and density scaling. All durations contain multiple Welch segments.
- Fixed-r SampEn uses embedding dimension 2 and `r = 0.2 × within-window population SD`.
- Sample-entropy-profile features use embedding dimension 2, three-decimal distance shells, and epsilon `1e-12`. Several entropy functions require at least `order + 2` samples; all three durations exceed this by a wide margin.
- Approximate entropy, fuzzy entropy, and distribution entropy contain explicit short-vector/zero-variance guards. Unexpected NaN or infinite output is audited rather than imputed.
- Duration-additive features such as waveform length and Total SampEn can change raw scale when sample count changes. Their algorithms are held fixed; this duration dependence is part of the sensitivity question and is not normalized away during extraction.

## Event-unit and model rules

- Stage 06 retains contraction, fetal-movement, and baseline classes and applies the exact frozen 34-feature completeness requirement.
- Event units are grouped by woman, recording, annotation, and class. If more than one eligible row exists, the complete row nearest the annotation midpoint is selected deterministically. All complete baseline windows are retained.
- The frozen model is REML `z_feature ~ C(window_class, Treatment(reference="baseline"))`, with a woman random intercept and recording variance component. The frozen fallback hierarchy is nested L-BFGS, nested BFGS, woman-only L-BFGS, then woman-only BFGS; the first clean fit is used, with a boundary fit retained within a tier when no clean fit exists.
- BH-FDR is applied separately to the 34 UC-versus-IN p-values at each duration.

## Locked comparison scale

The primary 40/60/80 coefficient comparison uses each feature's frozen 60-s Stage 06 mean and sample SD for all three durations. This common positive scaling preserves signs and p-values while expressing every coefficient in the same frozen 60-s SD unit. Raw-unit coefficients and within-duration standardized coefficients are also reported because some feature distributions are inherently duration-dependent. This choice was locked before inspecting 40-s or 80-s outcomes.

