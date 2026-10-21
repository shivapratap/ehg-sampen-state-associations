#!/usr/bin/env python3
"""Paper 1 reproducible pipeline -- Stage 02: Windowing / Segmentation.

Standalone port of pipeline/module1/*'s logic (frozen reference; not
imported at runtime -- see Paper1/MIGRATION_PLAN.md). Builds the 60-second
window manifest: event windows (contraction / possible_contraction /
fetal_movement -centered, purity/overlap classified) and baseline windows
(guard-interval-tiled quiet intervals).

**Scope note**: the old pipeline builds manifests for four window sizes
(40/60/80/100s). Per the fixed Paper 1 protocol (60s windows only), this
stage builds the 60s manifest exclusively -- see Paper1/MIGRATION_PLAN.md
Section 9, open decision 2. The 40/80/100s sensitivity sizes remain
available in the old (frozen) pipeline/module1/outputs/ if needed later.

Design source (logic ported, not imported): pipeline/module1/{run_module1,
windows,config}.py. This stage's window_manifest_60s.csv is validated
against the frozen pipeline/module1/outputs/window_manifest_60s.csv by the
reproducibility gate before Stage 03 onward is trusted.

Input:
    results/paper1/00_data_audit/recording_inventory.csv
    results/paper1/00_data_audit/recording_qc_flags.csv
    data/annotations/annotation_events.csv
    data/annotations/label_ontology.csv   (consistency check only)

Output:
    data/windows/window_manifest_60s.csv                          (primary, reusable)
    results/paper1/02_windowing/rejected_windows_60s.csv
    results/paper1/02_windowing/recording_processing_log_60s.csv
    results/paper1/02_windowing/cohort_summary_60s.md
    results/paper1/02_windowing/run_metadata_60s.json

Usage:
    python src/02_windowing.py
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
_SRC_DIR = Path(__file__).resolve().parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(1, str(_PROJECT_ROOT))

import config as project_config  # noqa: E402


# =============================================================================
# Section 1 -- Window construction primitives
# (ported from pipeline/module1/windows.py)
# =============================================================================

def usable_range(duration_s: float, edge_trim_s: float = project_config.EDGE_TRIM_S):
    return edge_trim_s, duration_s - edge_trim_s


@dataclass
class EventWindowResult:
    recording_id: str
    window_size_s: float
    annotation_index: int
    event_time_s: float
    event_type: str
    window_start_s: float
    window_end_s: float
    start_sample: int
    end_sample: int
    accepted: bool
    reject_reason: str = ""
    overlap_status: str = "clean"
    overlap_subtype: str = ""
    contains_competing_class: bool = False
    contains_same_class: bool = False
    movement_artifact_present: bool = False
    nearest_annotation_index: int = -1
    nearest_annotation_type: str = ""
    distance_to_nearest_annotation_s: float = float("inf")
    overlapping_event_annotation_index: int = -1
    overlapping_event_class: str = ""
    overlap_samples: int = 0
    overlap_duration_s: float = 0.0
    overlap_fraction: float = 0.0
    window_class: str = ""


def build_event_windows(events: list, recording_id: str, window_size_s: float,
                         duration_s: float, fs_hz: float) -> list:
    """events: list of dicts with keys annotation_index, time_s,
    normalized_event_type, for ONE recording, any order. Returns one
    EventWindowResult per PRIMARY event-type annotation (contraction /
    possible_contraction / fetal_movement) -- pm/em/pos never anchor their
    own window, they only contribute to the purity/overlap checks below.
    """
    half = window_size_s / 2.0
    u_start, u_end = usable_range(duration_s)
    results = []

    primary_events = [e for e in events if e["normalized_event_type"] in project_config.PRIMARY_EVENT_TYPES]

    for e in primary_events:
        t = e["time_s"]
        w_start, w_end = t - half, t + half
        start_sample = int(round(w_start * fs_hz))
        end_sample = int(round(w_end * fs_hz))

        r = EventWindowResult(
            recording_id=recording_id, window_size_s=window_size_s,
            annotation_index=e["annotation_index"], event_time_s=t,
            event_type=e["normalized_event_type"],
            window_start_s=w_start, window_end_s=w_end,
            start_sample=start_sample, end_sample=end_sample, accepted=True,
        )

        if w_start < u_start or w_end > u_end:
            r.accepted = False
            r.reject_reason = "outside_usable_range_after_edge_trim"
            results.append(r)
            continue

        best_overlap_dur = 0.0
        best_overlap_other = None
        nearest_gap = float("inf")
        nearest_other = None

        for o in events:
            if o is e:
                continue
            gap = abs(o["time_s"] - t)
            if gap < nearest_gap:
                nearest_gap = gap
                nearest_other = o

            contains = gap <= half
            windows_overlap = gap < window_size_s
            if windows_overlap:
                overlap_dur = window_size_s - gap
                if overlap_dur > best_overlap_dur:
                    best_overlap_dur = overlap_dur
                    best_overlap_other = o

            if contains:
                o_type = o["normalized_event_type"]
                if o_type in project_config.PRIMARY_EVENT_TYPES:
                    if o_type == r.event_type:
                        r.contains_same_class = True
                    else:
                        r.contains_competing_class = True
                else:
                    r.movement_artifact_present = True

        if nearest_other is not None:
            r.nearest_annotation_index = nearest_other["annotation_index"]
            r.nearest_annotation_type = nearest_other["normalized_event_type"]
            r.distance_to_nearest_annotation_s = nearest_gap

        if best_overlap_other is not None:
            r.overlapping_event_annotation_index = best_overlap_other["annotation_index"]
            r.overlapping_event_class = best_overlap_other["normalized_event_type"]
            r.overlap_duration_s = best_overlap_dur
            r.overlap_samples = int(round(best_overlap_dur * fs_hz))
            r.overlap_fraction = best_overlap_dur / window_size_s

        if r.contains_competing_class:
            r.overlap_status = "mixed"
        elif r.contains_same_class:
            r.overlap_status = "same_class_cluster"
        elif best_overlap_other is not None:
            r.overlap_status = "window_overlap"
            r.overlap_subtype = (
                "boundary_overlap" if best_overlap_dur <= project_config.BOUNDARY_OVERLAP_THRESHOLD_S
                else "substantive_overlap"
            )
        else:
            r.overlap_status = "clean"

        r.window_class = "mixed" if r.overlap_status == "mixed" else r.event_type
        results.append(r)

    return results


@dataclass
class BaselineWindowResult:
    recording_id: str
    window_size_s: float
    window_index: int
    window_start_s: float
    window_end_s: float
    start_sample: int
    end_sample: int
    nearest_annotation_index: int
    nearest_annotation_type: str
    distance_to_nearest_annotation_s: float
    guard_clearance_s: float
    tight_clearance: bool
    window_class: str = "baseline"


def _merge_intervals(intervals: list) -> list:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for s, e in intervals[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def compute_guard_intervals(events: list, guard_s: float = project_config.GUARD_INTERVAL_S) -> list:
    raw = [
        (e["time_s"] - guard_s, e["time_s"] + guard_s)
        for e in events
        if e["normalized_event_type"] in project_config.GUARD_ANNOTATION_TYPES
    ]
    return _merge_intervals(raw)


def _clean_gaps(u_start: float, u_end: float, forbidden: list) -> list:
    gaps = []
    cursor = u_start
    for f_start, f_end in forbidden:
        if f_start > cursor:
            gaps.append((cursor, min(f_start, u_end)))
        cursor = max(cursor, f_end)
        if cursor >= u_end:
            break
    if cursor < u_end:
        gaps.append((cursor, u_end))
    return [(s, e) for s, e in gaps if e > s]


def _nearest_annotation(t: float, events: list):
    best = None
    best_gap = float("inf")
    for e in events:
        gap = abs(e["time_s"] - t)
        if gap < best_gap:
            best_gap = gap
            best = e
    return best, best_gap


def build_baseline_windows(events: list, recording_id: str, window_size_s: float,
                            duration_s: float, fs_hz: float,
                            guard_s: float = project_config.GUARD_INTERVAL_S):
    """Returns (accepted_windows, rejected_gaps)."""
    u_start, u_end = usable_range(duration_s)
    forbidden = compute_guard_intervals(events, guard_s)
    gaps = _clean_gaps(u_start, u_end, forbidden)
    guard_events = [e for e in events if e["normalized_event_type"] in project_config.GUARD_ANNOTATION_TYPES]

    accepted, rejected_gaps = [], []
    idx = 0
    for gap_start, gap_end in gaps:
        gap_len = gap_end - gap_start
        if gap_len < window_size_s:
            rejected_gaps.append({
                "recording_id": recording_id, "window_size_s": window_size_s,
                "gap_start_s": gap_start, "gap_end_s": gap_end,
                "gap_duration_s": gap_len,
                "reject_reason": "clean_gap_shorter_than_window_size",
            })
            continue

        n_windows = int(gap_len // window_size_s)
        for i in range(n_windows):
            w_start = gap_start + i * window_size_s
            w_end = w_start + window_size_s
            start_sample = int(round(w_start * fs_hz))
            end_sample = int(round(w_end * fs_hz))

            nearest, _ = _nearest_annotation((w_start + w_end) / 2.0, guard_events)
            if nearest is not None:
                t = nearest["time_s"]
                if t < w_start:
                    edge_gap = w_start - t
                elif t > w_end:
                    edge_gap = t - w_end
                else:
                    edge_gap = 0.0
            else:
                edge_gap = float("inf")

            clearance = edge_gap - guard_s
            accepted.append(BaselineWindowResult(
                recording_id=recording_id, window_size_s=window_size_s, window_index=idx,
                window_start_s=w_start, window_end_s=w_end,
                start_sample=start_sample, end_sample=end_sample,
                nearest_annotation_index=nearest["annotation_index"] if nearest else -1,
                nearest_annotation_type=nearest["normalized_event_type"] if nearest else "",
                distance_to_nearest_annotation_s=edge_gap,
                guard_clearance_s=clearance,
                tight_clearance=clearance < project_config.TIGHT_CLEARANCE_THRESHOLD_S,
            ))
            idx += 1

    return accepted, rejected_gaps


# =============================================================================
# Section 2 -- Per-recording orchestration
# (ported from pipeline/module1/run_module1.py)
# =============================================================================

def _window_id(recording_id: str, window_size_s: float, start_sample: int, end_sample: int) -> str:
    return f"{recording_id}|W{window_size_s:g}|{start_sample}-{end_sample}"


def _process_recording_for_size(recording_id: str, duration_s: float, fs_hz: float,
                                  has_atr: bool, events_df: pd.DataFrame, window_size_s: float):
    events = events_df[events_df["recording_id"] == recording_id][
        ["annotation_index", "time_s", "normalized_event_type"]
    ].to_dict("records")

    event_results = build_event_windows(events, recording_id, window_size_s, duration_s, fs_hz)

    baseline_included = has_atr  # deliberate, flagged default -- same as old Module 1
    baseline_results, rejected_gaps = ([], [])
    if baseline_included:
        baseline_results, rejected_gaps = build_baseline_windows(
            events, recording_id, window_size_s, duration_s, fs_hz,
        )

    manifest_rows, rejected_rows = [], []

    for r in event_results:
        row = {
            "window_id": _window_id(recording_id, window_size_s, r.start_sample, r.end_sample),
            "recording_id": recording_id, "window_size_s": window_size_s,
            "kind": "event", "window_class": r.window_class, "event_type": r.event_type,
            "annotation_index": r.annotation_index, "event_time_s": r.event_time_s,
            "window_start_s": r.window_start_s, "window_end_s": r.window_end_s,
            "start_sample": r.start_sample, "end_sample": r.end_sample,
            "overlap_status": r.overlap_status, "overlap_subtype": r.overlap_subtype,
            "contains_competing_class": r.contains_competing_class,
            "contains_same_class": r.contains_same_class,
            "movement_artifact_present": r.movement_artifact_present,
            "nearest_annotation_index": r.nearest_annotation_index,
            "nearest_annotation_type": r.nearest_annotation_type,
            "distance_to_nearest_annotation_s": r.distance_to_nearest_annotation_s,
            "overlapping_event_annotation_index": r.overlapping_event_annotation_index,
            "overlapping_event_class": r.overlapping_event_class,
            "overlap_samples": r.overlap_samples, "overlap_duration_s": r.overlap_duration_s,
            "overlap_fraction": r.overlap_fraction,
            "guard_clearance_s": None, "tight_clearance": None,
        }
        if r.accepted:
            manifest_rows.append(row)
        else:
            rejected_rows.append({**row, "reject_reason": r.reject_reason})

    for w in baseline_results:
        manifest_rows.append({
            "window_id": _window_id(recording_id, window_size_s, w.start_sample, w.end_sample),
            "recording_id": recording_id, "window_size_s": window_size_s,
            "kind": "baseline", "window_class": "baseline", "event_type": "",
            "annotation_index": -1, "event_time_s": None,
            "window_start_s": w.window_start_s, "window_end_s": w.window_end_s,
            "start_sample": w.start_sample, "end_sample": w.end_sample,
            "overlap_status": "", "overlap_subtype": "",
            "contains_competing_class": False, "contains_same_class": False,
            "movement_artifact_present": False,
            "nearest_annotation_index": w.nearest_annotation_index,
            "nearest_annotation_type": w.nearest_annotation_type,
            "distance_to_nearest_annotation_s": w.distance_to_nearest_annotation_s,
            "overlapping_event_annotation_index": -1, "overlapping_event_class": "",
            "overlap_samples": 0, "overlap_duration_s": 0.0, "overlap_fraction": 0.0,
            "guard_clearance_s": w.guard_clearance_s, "tight_clearance": w.tight_clearance,
        })

    for g in rejected_gaps:
        rejected_rows.append({
            "window_id": "", "recording_id": recording_id, "window_size_s": window_size_s,
            "kind": "baseline_gap", "window_class": "", "event_type": "",
            "annotation_index": -1, "event_time_s": None,
            "window_start_s": g["gap_start_s"], "window_end_s": g["gap_end_s"],
            "start_sample": None, "end_sample": None,
            "overlap_status": "", "overlap_subtype": "",
            "contains_competing_class": False, "contains_same_class": False,
            "movement_artifact_present": False,
            "nearest_annotation_index": -1, "nearest_annotation_type": "",
            "distance_to_nearest_annotation_s": None,
            "overlapping_event_annotation_index": -1, "overlapping_event_class": "",
            "overlap_samples": 0, "overlap_duration_s": 0.0, "overlap_fraction": 0.0,
            "guard_clearance_s": None, "tight_clearance": None,
            "reject_reason": g["reject_reason"],
        })

    log_row = {
        "recording_id": recording_id, "window_size_s": window_size_s,
        "n_events_total": len(events),
        "n_event_windows_accepted": sum(1 for r in event_results if r.accepted),
        "n_event_windows_rejected": sum(1 for r in event_results if not r.accepted),
        "n_mixed": sum(1 for r in event_results if r.overlap_status == "mixed"),
        "n_same_class_cluster": sum(1 for r in event_results if r.overlap_status == "same_class_cluster"),
        "n_window_overlap": sum(1 for r in event_results if r.overlap_status == "window_overlap"),
        "baseline_included": baseline_included,
        "n_baseline_windows": len(baseline_results),
        "n_baseline_tight_clearance": sum(1 for w in baseline_results if w.tight_clearance),
        "n_rejected_baseline_gaps": len(rejected_gaps),
    }

    return manifest_rows, rejected_rows, log_row


# =============================================================================
# Section 3 -- Ontology consistency check
# (new -- cross-checks Stage 01's label_ontology.csv against the manifest's
#  actual window_class values, so the two artifacts cannot silently drift)
# =============================================================================

def _check_window_class_against_ontology(manifest_df: pd.DataFrame, ontology_df: pd.DataFrame) -> list:
    if manifest_df.empty:
        return []
    known_event_types = set(ontology_df["normalized_event_type"])
    observed = set(manifest_df["window_class"].dropna().unique()) - {""}
    unknown = sorted(observed - known_event_types)
    return unknown


# =============================================================================
# Section 4 -- Orchestration
# =============================================================================

def _load_inputs(audit_dir: Path, annotations_dir: Path):
    inventory = pd.read_csv(audit_dir / "recording_inventory.csv")
    qc = pd.read_csv(audit_dir / "recording_qc_flags.csv")
    events_df = pd.read_csv(annotations_dir / "annotation_events.csv")
    ontology_df = pd.read_csv(annotations_dir / "label_ontology.csv")
    run_metadata_path = audit_dir / "run_metadata.json"
    m0_run_metadata = json.loads(run_metadata_path.read_text()) if run_metadata_path.is_file() else {}
    return inventory, qc, events_df, ontology_df, m0_run_metadata


def run(window_size_s: int = project_config.WINDOW_SIZE_S,
        audit_dir: Path = None, annotations_dir: Path = None,
        windows_dir: Path = None, results_dir: Path = None) -> int:
    audit_dir = audit_dir or (project_config.AUDIT_DIR)
    annotations_dir = annotations_dir or project_config.ANNOTATIONS_DIR
    windows_dir = windows_dir or project_config.WINDOWS_DIR
    results_dir = results_dir or (project_config.WINDOWING_DIR)
    windows_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    inventory, qc, events_df, ontology_df, m0_run_metadata = _load_inputs(audit_dir, annotations_dir)

    merged = inventory.merge(qc[["recording_id", "recommended_action"]], on="recording_id", how="left")
    skip = merged["recommended_action"].isin(project_config.SKIP_RECOMMENDED_ACTIONS)
    n_skipped = int(skip.sum())
    active = merged[~skip].copy()

    print(f"[02_windowing] {len(active)}/{len(merged)} recordings active "
          f"({n_skipped} skipped by Stage 00 QC recommended_action)")
    print(f"[02_windowing] === window_size_s={window_size_s} ===")

    manifest_rows, rejected_rows, log_rows = [], [], []
    for i, rec in enumerate(active.itertuples(), 1):
        m_rows, r_rows, log_row = _process_recording_for_size(
            rec.recording_id, rec.duration_s, rec.fs_hz, bool(rec.has_atr),
            events_df, float(window_size_s),
        )
        manifest_rows.extend(m_rows)
        rejected_rows.extend(r_rows)
        log_rows.append(log_row)
        if i % 30 == 0 or i == len(active):
            print(f"[02_windowing]   {i}/{len(active)} recordings processed")

    manifest_df = pd.DataFrame(manifest_rows)
    rejected_df = pd.DataFrame(rejected_rows)
    log_df = pd.DataFrame(log_rows)

    manifest_df = manifest_df.merge(
        active[["recording_id", "patient_id", "record_type"]], on="recording_id", how="left",
    )

    unknown_classes = _check_window_class_against_ontology(manifest_df, ontology_df)
    if unknown_classes:
        print(f"[02_windowing] WARNING: window_class values with no ontology row: {unknown_classes}", file=sys.stderr)

    suffix = f"{int(window_size_s)}s"
    manifest_path = windows_dir / f"window_manifest_{suffix}.csv"
    manifest_df.to_csv(manifest_path, index=False)
    rejected_df.to_csv(results_dir / f"rejected_windows_{suffix}.csv", index=False)
    log_df.to_csv(results_dir / f"recording_processing_log_{suffix}.csv", index=False)

    print(f"[02_windowing]   -> {len(manifest_df)} manifest rows, "
          f"{len(rejected_df)} rejected rows written for {suffix}")

    _write_cohort_summary(results_dir, suffix, manifest_df, active, n_skipped, unknown_classes)
    _write_run_metadata(results_dir, suffix, window_size_s, m0_run_metadata, active, n_skipped)

    print(f"[02_windowing] done. window manifest -> {manifest_path}")
    print(f"[02_windowing] audit outputs -> {results_dir}")
    return 1 if unknown_classes else 0


def _write_cohort_summary(results_dir: Path, suffix: str, manifest_df: pd.DataFrame,
                            active: pd.DataFrame, n_skipped: int, unknown_classes: list) -> None:
    lines = [
        "# Paper 1 Reproducible Pipeline -- Stage 02 Cohort Summary",
        "",
        f"Recordings active: {len(active)} (skipped by Stage 00 QC: {n_skipped}).",
        f"Recordings with no annotation file, excluded from baseline generation "
        f"by default (flagged, not silently decided): {int((~active['has_atr']).sum())}.",
        "",
        "All counts below are computed from this run's manifest -- none are hand-typed.",
        "",
        f"## Window size: {suffix} (PRIMARY -- fixed protocol; only size built by this stage)",
        "",
    ]
    class_counts = manifest_df["window_class"].value_counts()
    for cls, count in class_counts.items():
        lines.append(f"- `{cls}`: {count}")
    n_tight = int((manifest_df["tight_clearance"] == True).sum())  # noqa: E712
    lines.append(f"- baseline windows flagged `tight_clearance`: {n_tight}")
    n_patients_with_windows = manifest_df["patient_id"].nunique() if not manifest_df.empty else 0
    lines.append(f"- distinct patients contributing at least one window: {n_patients_with_windows}")
    lines += ["", "## Ontology consistency check", ""]
    if unknown_classes:
        lines.append(f"- **FAIL**: window_class values with no Stage 01 ontology row: {unknown_classes}")
    else:
        lines.append("- **PASS**: every window_class value in the manifest has a Stage 01 ontology row.")
    lines += [
        "",
        "## Reproducibility note",
        "",
        "This stage's window_manifest_60s.csv is validated against the frozen "
        "pipeline/module1/outputs/window_manifest_60s.csv by the "
        "reproducibility gate before Stage 03 onward is trusted. See "
        "results/paper1/validation/.",
    ]
    (results_dir / f"cohort_summary_{suffix}.md").write_text("\n".join(lines) + "\n")


def _get_git_commit(repo_root: Path):
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(repo_root),
            capture_output=True, text=True, timeout=10,
        )
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:  # noqa: BLE001
        pass
    return None


def _write_run_metadata(results_dir: Path, suffix: str, window_size_s: int,
                          m0_run_metadata: dict, active: pd.DataFrame, n_skipped: int) -> None:
    metadata = {
        "stage": "02_windowing",
        "stage00_run_metadata_hashes": m0_run_metadata.get("content_hashes", {}),
        "config": {
            "window_size_s": window_size_s,
            "primary_band": project_config.PRIMARY_BAND,
            "guard_interval_s": project_config.GUARD_INTERVAL_S,
            "edge_trim_s": project_config.EDGE_TRIM_S,
            "tight_clearance_threshold_s": project_config.TIGHT_CLEARANCE_THRESHOLD_S,
            "boundary_overlap_threshold_s": project_config.BOUNDARY_OVERLAP_THRESHOLD_S,
            "scope_note": "60s only, per fixed Paper 1 protocol -- see MIGRATION_PLAN.md",
        },
        "code_version": {"git_commit": _get_git_commit(project_config.PROJECT_ROOT)},
        "execution_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "n_recordings_active": len(active),
        "n_recordings_skipped_by_stage00_qc": n_skipped,
    }
    (results_dir / f"run_metadata_{suffix}.json").write_text(json.dumps(metadata, indent=2, default=str))


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--window-size", type=int, default=project_config.WINDOW_SIZE_S)
    parser.add_argument("--audit-dir", type=Path, default=None)
    parser.add_argument("--annotations-dir", type=Path, default=None)
    parser.add_argument("--windows-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    sys.exit(run(
        window_size_s=args.window_size, audit_dir=args.audit_dir,
        annotations_dir=args.annotations_dir, windows_dir=args.windows_dir,
        results_dir=args.results_dir,
    ))
