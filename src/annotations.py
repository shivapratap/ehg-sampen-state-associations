#!/usr/bin/env python3
"""EHG sample-entropy pipeline -- Stage 01: Label Ontology Consolidation.

NEW script -- no direct 1:1 predecessor in the old pipeline. Makes the
wfdb-annotation-symbol -> class_label mapping used throughout this project
one explicit, auditable table, instead of leaving it implicit across
pipeline/module0/config.py's ANNOTATION_LABELS dict and the class-assignment
branches inside pipeline/module1/windows.py (window_class = event_type, or
"mixed" if overlap_status == "mixed").

Fixed the standard pipeline protocol classes:
    UC  uterine contraction   (wfdb symbol "C")
    FM  fetal movement        (wfdb symbol "fm")
    IN  inactivity / baseline (unannotated, guard-cleared windows)

`possible_contraction` ("(c)"), the movement/equipment/position markers
(pm/em/pos), and "mixed" event windows are explicitly excluded from the
primary UC/FM/IN protocol here, with the reason recorded -- not silently
dropped downstream. See the design notes.

Input:  data/annotations/annotation_events.csv (Stage 00 output)
Output: data/annotations/label_ontology.csv
        results/paper1/01_label_ontology/label_ontology_report.md
        results/paper1/01_label_ontology/annotation_events_with_class.csv

Usage:
    python src/01_label_ontology.py
"""
from __future__ import annotations

import argparse
import sys
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
# The ontology table itself -- the one place this mapping is defined.
# =============================================================================

def build_ontology_table() -> pd.DataFrame:
    """One row per wfdb annotation symbol (plus a synthetic 'baseline' row
    for unannotated, guard-cleared windows), mapping to the fixed the standard pipeline
    protocol class label, or explicitly documenting exclusion.
    """
    rows = [
        {
            "wfdb_symbol": "C",
            "normalized_event_type": "contraction",
            "class_label": "UC",
            "in_primary_protocol": True,
            "exclusion_reason": "",
            "notes": "Primary confirmatory class -- uterine contraction.",
        },
        {
            "wfdb_symbol": "(c)",
            "normalized_event_type": "possible_contraction",
            "class_label": None,
            "in_primary_protocol": False,
            "exclusion_reason": "ambiguous/weak annotation; not part of the fixed UC/FM/IN protocol",
            "notes": "Retained in annotation_events.csv and used as a guard-interval "
                     "annotation type in windowing, but does not anchor its own "
                     "class in this pipeline's primary analysis.",
        },
        {
            "wfdb_symbol": "fm",
            "normalized_event_type": "fetal_movement",
            "class_label": "FM",
            "in_primary_protocol": True,
            "exclusion_reason": "",
            "notes": "Secondary confirmatory class -- fetal movement, modeled "
                     "separately from UC and never pooled with it.",
        },
        {
            "wfdb_symbol": "pm",
            "normalized_event_type": "participant_movement",
            "class_label": None,
            "in_primary_protocol": False,
            "exclusion_reason": "movement-artifact marker, not an analysis class",
            "notes": "Contributes to guard intervals and event-window purity "
                     "checks only.",
        },
        {
            "wfdb_symbol": "em",
            "normalized_event_type": "equipment_manipulation",
            "class_label": None,
            "in_primary_protocol": False,
            "exclusion_reason": "equipment-artifact marker, not an analysis class",
            "notes": "Contributes to guard intervals and event-window purity "
                     "checks only.",
        },
        {
            "wfdb_symbol": "pos",
            "normalized_event_type": "position_change",
            "class_label": None,
            "in_primary_protocol": False,
            "exclusion_reason": "position-change marker, not an analysis class",
            "notes": "Contributes to guard intervals and event-window purity "
                     "checks only.",
        },
        {
            "wfdb_symbol": "",
            "normalized_event_type": "baseline",
            "class_label": "IN",
            "in_primary_protocol": True,
            "exclusion_reason": "",
            "notes": "Synthetic row -- unannotated, guard-cleared interval "
                     "(Stage 02 windowing rule), not a wfdb annotation symbol.",
        },
        {
            "wfdb_symbol": "<mixed window>",
            "normalized_event_type": "mixed",
            "class_label": None,
            "in_primary_protocol": False,
            "exclusion_reason": "window contains a competing-class annotation; "
                                 "not a clean single-class unit",
            "notes": "Synthetic row -- window_class assigned by Stage 02 "
                     "windowing when a candidate window's overlap_status is "
                     "'mixed', not a wfdb annotation symbol.",
        },
    ]
    df = pd.DataFrame(rows)
    # Cross-check against config.CLASS_MAP / EXCLUDED_WINDOW_CLASSES so this
    # table and config.py cannot silently drift apart.
    for event_type, class_label in project_config.CLASS_MAP.items():
        match = df.loc[df["normalized_event_type"] == event_type, "class_label"]
        assert not match.empty, f"config.CLASS_MAP has {event_type!r} not present in ontology table"
        assert match.iloc[0] == class_label, (
            f"config.CLASS_MAP maps {event_type!r} -> {class_label!r} but ontology "
            f"table has {match.iloc[0]!r}"
        )
    for event_type in project_config.EXCLUDED_WINDOW_CLASSES:
        match = df.loc[df["normalized_event_type"] == event_type, "in_primary_protocol"]
        assert not match.empty, f"config.EXCLUDED_WINDOW_CLASSES has {event_type!r} not present in ontology table"
        assert match.iloc[0] == False, (  # noqa: E712
            f"config.EXCLUDED_WINDOW_CLASSES lists {event_type!r} but ontology "
            "table marks it in_primary_protocol=True"
        )
    return df


# =============================================================================
# Coverage / consistency checks
# =============================================================================

def load_annotation_events(annotations_dir: Path) -> pd.DataFrame:
    events_path = annotations_dir / "annotation_events.csv"
    if not events_path.is_file():
        raise FileNotFoundError(
            f"{events_path} not found. Run src/00_data_audit.py first -- "
            "Stage 01 consumes its annotation_events.csv output."
        )
    return pd.read_csv(events_path)


def check_symbol_coverage(events_df: pd.DataFrame, ontology_df: pd.DataFrame) -> list:
    """Every wfdb_symbol actually present in the dataset's annotations must
    have exactly one ontology row. Returns a list of any symbols found in
    the data but absent from the ontology table (expected: empty list)."""
    if events_df.empty:
        return []
    observed_symbols = set(events_df["wfdb_symbol"].dropna().unique())
    known_symbols = set(ontology_df.loc[ontology_df["wfdb_symbol"] != "", "wfdb_symbol"])
    # "<mixed window>" and "" (baseline) are synthetic rows, not real wfdb
    # symbols, so they're excluded from the "known" comparison set above.
    known_symbols.discard("<mixed window>")
    unmapped = sorted(observed_symbols - known_symbols)
    return unmapped


def apply_ontology(events_df: pd.DataFrame, ontology_df: pd.DataFrame) -> pd.DataFrame:
    """Left-join annotation_events.csv against the ontology table on
    wfdb_symbol, attaching class_label and in_primary_protocol. Rows for
    unmapped symbols (should not occur if check_symbol_coverage is clean)
    get class_label=NaN, in_primary_protocol=NaN rather than being dropped.
    """
    lookup = ontology_df.loc[ontology_df["wfdb_symbol"] != ""].set_index("wfdb_symbol")
    merged = events_df.copy()
    merged["class_label"] = merged["wfdb_symbol"].map(lookup["class_label"])
    merged["in_primary_protocol"] = merged["wfdb_symbol"].map(lookup["in_primary_protocol"])
    return merged


# =============================================================================
# Report
# =============================================================================

def _write_report(path: Path, ontology_df: pd.DataFrame, events_df: pd.DataFrame,
                   unmapped_symbols: list, class_counts: pd.Series) -> None:
    lines = [
        "# EHG Sample-Entropy Pipeline -- Stage 01 Label Ontology Report",
        "",
        "This table is the single source of truth for wfdb-annotation-symbol "
        "to the standard pipeline protocol class (UC/FM/IN) mapping. It is new relative to "
        "the old pipeline, which left this mapping implicit -- see "
        "the design notes.",
        "",
        "## Ontology table",
        "",
        "| wfdb_symbol | normalized_event_type | class_label | in_primary_protocol | exclusion_reason |",
        "|---|---|---|---|---|",
    ]
    for _, r in ontology_df.iterrows():
        sym = r["wfdb_symbol"] if r["wfdb_symbol"] else "*(none -- baseline)*"
        cls = r["class_label"] if pd.notna(r["class_label"]) else "--"
        lines.append(f"| `{sym}` | {r['normalized_event_type']} | {cls} | {r['in_primary_protocol']} | {r['exclusion_reason'] or '--'} |")

    lines += ["", "## Symbol coverage check", ""]
    if events_df.empty:
        lines.append(
            "- `annotation_events.csv` is empty (no raw data was available "
            "when Stage 00 ran) -- coverage check skipped."
        )
    elif unmapped_symbols:
        lines.append(f"- **FAIL**: {len(unmapped_symbols)} symbol(s) found in the data with no ontology row: {unmapped_symbols}")
    else:
        lines.append(f"- **PASS**: every distinct wfdb_symbol in `annotation_events.csv` ({events_df['wfdb_symbol'].nunique()} distinct symbols) has an ontology row.")

    lines += ["", "## Class distribution (annotation-anchored events only; excludes synthetic baseline/mixed rows)", ""]
    if class_counts.empty:
        lines.append("- No annotation events available to summarize.")
    else:
        for cls, count in class_counts.items():
            lines.append(f"- {cls}: {count}")

    lines += [
        "",
        "## Downstream usage",
        "",
        "- Stage 02 (`02_windowing.py`) uses this table to decide which "
        "annotation types anchor event windows (`in_primary_protocol == True` "
        "types, plus `possible_contraction` and movement/equipment/position "
        "markers still contribute to guard intervals and purity checks per "
        "the original Module 1 windowing rules).",
        "- Stage 06 (`06_mixed_models.py`) uses this table to confirm no "
        "non-protocol class (`possible_contraction`, `mixed`) leaks into the "
        "confirmatory UC/FM/IN analysis.",
    ]
    path.write_text("\n".join(lines) + "\n")


# =============================================================================
# Orchestration
# =============================================================================

def run(annotations_dir: Path = None, results_dir: Path = None) -> int:
    annotations_dir = annotations_dir or project_config.ANNOTATIONS_DIR
    results_dir = results_dir or (project_config.ANNOTATIONS_DIR)
    annotations_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    print("[01_label_ontology] building ontology table")
    ontology_df = build_ontology_table()
    ontology_df.to_csv(annotations_dir / "label_ontology.csv", index=False)

    print("[01_label_ontology] loading annotation events from Stage 00")
    events_df = load_annotation_events(annotations_dir)

    unmapped = check_symbol_coverage(events_df, ontology_df)
    if unmapped:
        print(f"[01_label_ontology] WARNING: unmapped symbols found: {unmapped}", file=sys.stderr)

    merged = apply_ontology(events_df, ontology_df)
    merged.to_csv(results_dir / "annotation_events_with_class.csv", index=False)

    class_counts = (
        merged.loc[merged["class_label"].notna(), "class_label"].value_counts()
        if not merged.empty else pd.Series(dtype=int)
    )

    _write_report(results_dir / "label_ontology_report.md", ontology_df, events_df, unmapped, class_counts)

    print(f"[01_label_ontology] ontology table written to {annotations_dir / 'label_ontology.csv'}")
    print(f"[01_label_ontology] report written to {results_dir / 'label_ontology_report.md'}")
    if unmapped:
        print("[01_label_ontology] done WITH warnings (unmapped symbols present).")
        return 1
    print("[01_label_ontology] done. All annotation symbols covered.")
    return 0


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--annotations-dir", type=Path, default=None)
    parser.add_argument("--results-dir", type=Path, default=None)
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    sys.exit(run(annotations_dir=args.annotations_dir, results_dir=args.results_dir))
