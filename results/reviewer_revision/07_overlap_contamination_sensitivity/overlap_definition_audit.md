# Overlap-definition audit

## Frozen inputs and metadata

The primary population is the frozen 60-s EHG9 event-unit table. Stage 02 carries the overlap fields unchanged through Stage 04 and into Stage 06: `overlap_status`, `overlap_subtype`, `contains_competing_class`, `contains_same_class`, `movement_artifact_present`, `nearest_annotation_index`, `nearest_annotation_type`, `distance_to_nearest_annotation_s`, `overlapping_event_annotation_index`, `overlapping_event_class`, `overlap_samples`, `overlap_duration_s`, and `overlap_fraction`.

Stage 02 defines an event as *contained* when its centre is within 30 s of the anchor. A different primary event class in that interval sets `contains_competing_class=True`, `overlap_status=mixed`, and `window_class=mixed`; Stage 06 excludes these windows. Event-window footprints overlap when their centres are less than 60 s apart. For retained UC/FM windows, `window_overlap` therefore represents centre separation greater than 30 s but less than 60 s. `overlap_duration_s = 60 - centre_gap` and `overlap_fraction = overlap_duration_s / 60`.

`overlapping_event_class` records the event producing the greatest overlap. An independent reconstruction from all frozen annotation times found no competing UC/FM overlap hidden by a closer same-class or artifact annotation; the metadata rule and the all-annotation reconstruction agreed for every retained UC and FM unit.

## Locked classification rules

The rule was fixed before model outcomes were examined:

- **Clean UC:** `window_class == contraction` and neither `contains_competing_class` nor (`overlap_duration_s > 0` and `overlapping_event_class == fetal_movement`).
- **Contaminated UC:** `window_class == contraction` and either the generic contained-competing flag is true or the frozen greatest-overlap class is `fetal_movement` with positive overlap duration.
- **Clean FM:** `window_class == fetal_movement` and neither `contains_competing_class` nor (`overlap_duration_s > 0` and `overlapping_event_class == contraction`).
- **Contaminated FM:** `window_class == fetal_movement` and either the generic contained-competing flag is true or the frozen greatest-overlap class is `contraction` with positive overlap duration.

Same-class overlap is retained. Participant movement, equipment manipulation, and position-change metadata are not used to define UC/FM cross-contamination. `movement_artifact_present` is reported but is not an exclusion criterion in this targeted analysis.

## Frozen-population check before modelling

- UC: 333 total, 312 clean, 21 contaminated.
- FM: 392 total, 371 clean, 21 contaminated.
- Retained event units with `contains_competing_class=True`: 0. These are zero because Stage 02 had already labelled contained competing-class events as mixed.

## Model retention rules

For clean-UC inference, contaminated UC units alone are removed; all clean UC, all frozen FM units, and all IN units remain so the exact three-state model and its per-feature standardization convention are retained. For clean-FM inference, contaminated FM units alone are removed; all frozen UC units and all IN units remain. Feature means and sample SDs are recalculated over the observations entering each three-state sensitivity fit, exactly as Stage 06 standardizes each fitted dataset.
