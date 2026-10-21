"""Stage 05 -- Feature QC across all four channels (EHG9/10/11/12).

Standalone stage: reads the four per-channel Stage 04 outputs
(results/paper1/04_feature_extraction/window_features_60s_{channel}.csv),
never touches pipeline/ or re-derives features.

Per MIGRATION_PLAN.md Section 3 (05_feature_qc.py row), this stage performs:
  - per-channel missingness (NaN count/rate per feature column)
  - duplicate (window_id, signal_channel) detection
  - cross-channel row-count parity (all four channels must share the same
    4283 window_ids -- window structure is channel-independent, only
    feature values differ)
  - sample-count integrity vs. the Stage 02 window manifest
    (filtered_sample_count vs. expected window duration at the analysis
    sampling rate)
  - feature-column-set drift check against config.FEATURE_COLUMNS (the
    87-column registry)
  - distribution/NaN summary per channel (lightweight text/CSV form; no
    plotting dependency required)

It also builds the canonical long-format table used by later stages
(one row per woman x recording x event x window x channel): this is the
union of the four per-channel wide CSVs with a `channel` column added,
stacked long-wise. Not explicitly named in the migration plan as a
separate artifact, but required by the original Paper 1 design brief
("one canonical long-format feature table ... one row per woman x
recording x event x window x channel"). Written as
results/paper1/05_feature_qc/window_features_60s_long.csv.

Outputs (results/paper1/05_feature_qc/):
  qc_missingness_60s.csv         per (channel, feature) NaN count/rate
  qc_duplicates_60s.csv          any duplicate (window_id, channel) rows found (empty if none)
  qc_row_count_parity_60s.csv    per-channel row count + window_id-set match vs EHG9
  qc_integrity_report_60s.md     sample-count integrity + feature-column-set drift summary
  qc_report_60s.md               human-readable roll-up of all checks, pass/fail per check
  window_features_60s_long.csv   canonical long-format table (all 4 channels stacked)
  run_metadata_60s.json          run provenance
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
_SRC_DIR = Path(__file__).resolve().parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(1, str(_PROJECT_ROOT))

import config as project_config  # noqa: E402

LOGGER = logging.getLogger("stage05_feature_qc")


# --- Identifier / merge columns kept in the long-format table --------------
# (everything up to but not including the feature columns; feature columns
# are appended per config.FEATURE_COLUMNS)
ID_COLUMNS = (
    "window_id", "recording_id", "patient_id", "record_type", "signal_channel",
    "window_size_s", "kind", "window_class", "event_type", "annotation_index",
    "event_time_s", "window_start_s", "window_end_s", "start_sample", "end_sample",
    "overlap_status", "overlap_subtype", "contains_competing_class", "contains_same_class",
    "movement_artifact_present", "nearest_annotation_index", "nearest_annotation_type",
    "distance_to_nearest_annotation_s", "overlapping_event_annotation_index",
    "overlapping_event_class", "overlap_samples", "overlap_duration_s", "overlap_fraction",
    "guard_clearance_s", "tight_clearance",
    "input_sampling_frequency_hz", "filtered_sampling_frequency_hz", "filtered_sample_count",
)


def _configure_logging(results_dir: Path) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    fmt = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(fmt)
    file_handler = logging.FileHandler(results_dir / "stage05_60s.log", mode="a")
    file_handler.setFormatter(fmt)
    root.handlers = [stream_handler, file_handler]


def _load_channel(features_dir: Path, channel: str) -> pd.DataFrame:
    path = features_dir / f"window_features_60s_{channel}.csv"
    if not path.is_file():
        raise FileNotFoundError(f"Stage 04 output missing for {channel}: {path}")
    df = pd.read_csv(path)
    df["channel"] = channel
    return df


def check_missingness(frames: dict) -> pd.DataFrame:
    rows = []
    for channel, df in frames.items():
        n = len(df)
        for col in project_config.FEATURE_COLUMNS:
            nan_count = int(df[col].isna().sum()) if col in df.columns else n
            rows.append({
                "channel": channel, "feature": col, "n_rows": n,
                "nan_count": nan_count,
                "nan_rate": (nan_count / n) if n else float("nan"),
                "is_fragile_feature": col in project_config.NUMERIC_SENSITIVITY_FEATURES,
            })
    return pd.DataFrame(rows)


def check_duplicates(frames: dict) -> pd.DataFrame:
    dup_frames = []
    for channel, df in frames.items():
        dup_mask = df.duplicated(subset=["window_id", "signal_channel"], keep=False)
        if dup_mask.any():
            dup_frames.append(df.loc[dup_mask, ["window_id", "recording_id", "signal_channel"]])
    if dup_frames:
        return pd.concat(dup_frames, ignore_index=True)
    return pd.DataFrame(columns=["window_id", "recording_id", "signal_channel"])


def check_row_count_parity(frames: dict, reference_channel: str) -> pd.DataFrame:
    ref_ids = set(frames[reference_channel]["window_id"]) if reference_channel in frames else set()
    rows = []
    for channel, df in frames.items():
        ids = set(df["window_id"])
        rows.append({
            "channel": channel,
            "n_rows": len(df),
            "n_unique_window_ids": len(ids),
            "matches_reference_window_id_set": ids == ref_ids,
            "reference_channel": reference_channel,
            "missing_vs_reference": len(ref_ids - ids),
            "extra_vs_reference": len(ids - ref_ids),
        })
    return pd.DataFrame(rows)


def check_sample_count_integrity(frames: dict, window_size_s: int) -> list:
    """filtered_sample_count should equal window_size_s * analysis sampling rate
    for every row of every channel (all Stage 02 windows are fixed-length by
    construction; only the leading/trailing recording windows near
    guard/edge trims should ever legitimately differ, and Stage 02 already
    excludes those by construction, so none are expected here)."""
    expected = int(round(window_size_s * project_config.ANALYSIS_SAMPLING_FREQUENCY_HZ))
    findings = []
    for channel, df in frames.items():
        if "filtered_sample_count" not in df.columns:
            findings.append(f"- **{channel}**: `filtered_sample_count` column missing -- cannot check.")
            continue
        bad = df[df["filtered_sample_count"] != expected]
        if bad.empty:
            findings.append(f"- **{channel}**: PASS -- all {len(df)} rows have `filtered_sample_count == {expected}`.")
        else:
            findings.append(
                f"- **{channel}**: FAIL -- {len(bad)}/{len(df)} rows have `filtered_sample_count != {expected}` "
                f"(e.g. window_id={bad['window_id'].iloc[0]!r} count={bad['filtered_sample_count'].iloc[0]})."
            )
    return findings


def check_feature_column_drift(frames: dict) -> list:
    expected = set(project_config.FEATURE_COLUMNS)
    findings = []
    for channel, df in frames.items():
        present = set(df.columns) & expected
        missing = expected - present
        extra_feature_like = {c for c in df.columns if c not in ID_COLUMNS and c not in expected and c != "channel"}
        if not missing and not extra_feature_like:
            findings.append(
                f"- **{channel}**: PASS -- all {len(project_config.FEATURE_COLUMNS)} registry columns present, "
                "no unexpected feature-like columns."
            )
        else:
            findings.append(
                f"- **{channel}**: FAIL -- missing={sorted(missing)}, unexpected={sorted(extra_feature_like)}."
            )
    return findings


def build_long_table(frames: dict) -> pd.DataFrame:
    parts = []
    for channel, df in frames.items():
        cols = [c for c in ID_COLUMNS if c in df.columns] + ["channel"] + list(project_config.FEATURE_COLUMNS)
        cols = [c for c in cols if c in df.columns or c == "channel"]
        parts.append(df[cols])
    long_df = pd.concat(parts, ignore_index=True)
    long_df = long_df.sort_values(["patient_id", "recording_id", "window_id", "channel"]).reset_index(drop=True)
    return long_df


def _write_qc_report(
    results_dir: Path, missingness: pd.DataFrame, duplicates: pd.DataFrame,
    parity: pd.DataFrame, integrity_findings: list, drift_findings: list,
    long_df: pd.DataFrame, channels: tuple,
) -> str:
    lines = ["# Paper 1 Reproducible Pipeline -- Stage 05 Feature QC Report", ""]
    lines.append(f"Channels checked: {', '.join(channels)}")
    lines.append(f"Canonical long-format table: {len(long_df)} rows "
                 f"({long_df['recording_id'].nunique()} recordings x {len(channels)} channels).")
    lines.append("")

    overall_pass = True

    lines.append("## 1. Row-count parity (window_id sets identical across channels)")
    lines.append("")
    lines.append(parity.to_markdown(index=False))
    if not parity["matches_reference_window_id_set"].all():
        overall_pass = False
        lines.append("\n**FAIL** -- at least one channel's window_id set does not match the reference.")
    else:
        lines.append("\n**PASS** -- all channels share an identical window_id set.")
    lines.append("")

    lines.append("## 2. Duplicate (window_id, channel) detection")
    lines.append("")
    if duplicates.empty:
        lines.append("**PASS** -- no duplicate (window_id, channel) pairs found in any channel.")
    else:
        overall_pass = False
        lines.append(f"**FAIL** -- {len(duplicates)} duplicate row(s) found. See qc_duplicates_60s.csv.")
    lines.append("")

    lines.append("## 3. Sample-count integrity vs. manifest")
    lines.append("")
    lines.extend(integrity_findings)
    if any("FAIL" in f for f in integrity_findings):
        overall_pass = False
    lines.append("")

    lines.append(f"## 4. Feature-column-set drift vs. config.FEATURE_COLUMNS ({len(project_config.FEATURE_COLUMNS)} columns)")
    lines.append("")
    lines.extend(drift_findings)
    if any("FAIL" in f for f in drift_findings):
        overall_pass = False
    lines.append("")

    lines.append("## 5. Missingness summary")
    lines.append("")
    per_channel_nan = missingness.groupby("channel")["nan_count"].sum()
    fragile_nan = missingness[missingness["is_fragile_feature"]].groupby("channel")["nan_count"].sum()
    for channel in channels:
        total_nan = int(per_channel_nan.get(channel, 0))
        frag_nan = int(fragile_nan.get(channel, 0))
        lines.append(f"- **{channel}**: {total_nan} total NaN cells across {len(project_config.FEATURE_COLUMNS)} feature columns "
                      f"({frag_nan} in fragile-feature columns).")
    if per_channel_nan.sum() == 0:
        lines.append("\n**PASS** -- zero missingness across all channels and feature columns.")
    else:
        lines.append("\n**INFO** -- non-zero missingness present; see qc_missingness_60s.csv for the full "
                      "per-(channel, feature) breakdown before deciding whether this blocks downstream analysis.")
    lines.append("")

    lines.insert(2, f"**Overall QC outcome: {'PASS' if overall_pass else 'FAIL'}**")
    lines.insert(3, "")

    report = "\n".join(lines)
    (results_dir / "qc_report_60s.md").write_text(report)
    return "PASS" if overall_pass else "FAIL"


def run(
    channels: tuple = project_config.CHANNELS,
    features_dir: Path = None,
    results_dir: Path = None,
    window_size_s: int = project_config.WINDOW_SIZE_S,
) -> int:
    features_dir = features_dir or (project_config.FEATURES_DIR)
    results_dir = results_dir or (project_config.QC_DIR)
    results_dir.mkdir(parents=True, exist_ok=True)
    _configure_logging(results_dir)

    LOGGER.info("Stage 05: loading Stage 04 outputs for channels=%s from %s", channels, features_dir)
    frames = {channel: _load_channel(features_dir, channel) for channel in channels}
    for channel, df in frames.items():
        LOGGER.info("  %s: %d rows", channel, len(df))

    missingness = check_missingness(frames)
    duplicates = check_duplicates(frames)
    parity = check_row_count_parity(frames, reference_channel=project_config.PRIMARY_CHANNEL)
    integrity_findings = check_sample_count_integrity(frames, window_size_s)
    drift_findings = check_feature_column_drift(frames)
    long_df = build_long_table(frames)

    missingness.to_csv(results_dir / "qc_missingness_60s.csv", index=False)
    duplicates.to_csv(results_dir / "qc_duplicates_60s.csv", index=False)
    parity.to_csv(results_dir / "qc_row_count_parity_60s.csv", index=False)
    long_df.to_csv(results_dir / "window_features_60s_long.csv", index=False)

    integrity_report = "\n".join(
        ["# Stage 05 -- Sample-count integrity & feature-column-set drift", "",
         "## Sample-count integrity", ""] + integrity_findings +
        ["", "## Feature-column-set drift", ""] + drift_findings
    )
    (results_dir / "qc_integrity_report_60s.md").write_text(integrity_report)

    outcome = _write_qc_report(
        results_dir, missingness, duplicates, parity, integrity_findings, drift_findings,
        long_df, channels,
    )

    metadata = {
        "stage": "05_feature_qc", "window_size_s": window_size_s, "channels": list(channels),
        "n_rows_per_channel": {c: len(df) for c, df in frames.items()},
        "n_long_table_rows": len(long_df),
        "overall_outcome": outcome,
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }
    (results_dir / "run_metadata_60s.json").write_text(json.dumps(metadata, indent=2, default=str))

    LOGGER.info("Stage 05 done. outcome=%s outputs -> %s", outcome, results_dir)
    return 0 if outcome == "PASS" else 1


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--channels", type=str, default=None, help="Comma-separated, e.g. EHG9,EHG10")
    parser.add_argument("--features-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--window-size", type=int, default=project_config.WINDOW_SIZE_S)
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    channels = tuple(c.strip() for c in args.channels.split(",")) if args.channels else project_config.CHANNELS
    sys.exit(run(
        channels=channels, features_dir=args.features_dir, results_dir=args.results_dir,
        window_size_s=args.window_size,
    ))
