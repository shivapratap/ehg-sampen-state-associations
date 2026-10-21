#!/usr/bin/env python3
"""Paper 1 reproducible pipeline -- Stage 00: Data & Annotation Audit.

Standalone port of pipeline/module0/*'s logic (frozen reference; not
imported at runtime -- see Paper1/MIGRATION_PLAN.md). Validates raw
PhysioNet recordings against their own RECORDS/SHA256SUMS manifests, parses
headers/annotations/clinical text, builds the raw-position -> canonical
-channel mapping, and aggregates per-recording QC severity into a single
audit. This underwrites every later stage: channel_mapping.csv decides
which raw column is "EHG9"/"EHG10"/"EHG11"/"EHG12", and annotation_events.csv
seeds Stage 01's label ontology and Stage 02's windowing.

Design source (logic ported, not imported): pipeline/module0/{run_module0,
headers,annotations,channel_mapping,checksums,qc,toco,file_inventory,
provenance}.py. This script's EHG9-relevant outputs (recording_inventory.csv,
channel_mapping.csv, annotation_events.csv) are validated against that
module's frozen outputs by the reproducibility gate before any of the later
Paper 1 stages are trusted -- see results/paper1/validation/.

Two-level failure policy (unchanged from the old pipeline):
  * Dataset-level fatal errors (RECORDS unparseable, SHA256SUMS unparseable,
    duplicate recording identifiers) abort the entire run immediately.
  * Recording-level structural errors (bad header, missing .dat, annotation
    decode failure, etc.) are caught, the recording is flagged and
    quarantined, and the run continues through the rest of the cohort.

Usage:
    python src/00_data_audit.py
    EHG_RAW_DIR=/path/to/raw python src/00_data_audit.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import traceback
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy
import pandas as pd
import wfdb
from PIL import Image, UnidentifiedImageError

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
_SRC_DIR = Path(__file__).resolve().parent
if str(_SRC_DIR) not in sys.path:
    sys.path.insert(0, str(_SRC_DIR))
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(1, str(_PROJECT_ROOT))

import config as project_config  # noqa: E402
from src.common.hashing import (  # noqa: E402
    hash_directory_listing,
    sha256_of_bytes,
    sha256_of_file,
)


# =============================================================================
# Section 1 -- File inventory / RECORDS reconciliation
# (ported from pipeline/module0/file_inventory.py)
# =============================================================================

_RECORDING_STEM_PREFIX_RE = re.compile(r"^ice\d{3}_")


@dataclass
class ReconciliationResult:
    records_listed: set
    hea_found: set
    dat_found: set
    atr_found: set
    jpg_found: set

    missing_hea_for_records: set = field(init=False)
    hea_not_in_records: set = field(init=False)
    duplicate_in_records_file: list = field(init=False, default_factory=list)
    incomplete_groups: dict = field(init=False)

    def __post_init__(self):
        self.missing_hea_for_records = self.records_listed - self.hea_found
        self.hea_not_in_records = self.hea_found - self.records_listed
        self.incomplete_groups = {}
        for stem in sorted(self.records_listed | self.hea_found):
            missing = []
            if stem not in self.hea_found:
                missing.append(".hea")
            if stem not in self.dat_found:
                missing.append(".dat")
            if stem not in self.jpg_found:
                missing.append(".jpg")
            # .atr intentionally not required -- absence is normal (~91% coverage).
            if missing:
                self.incomplete_groups[stem] = missing

    def summary(self) -> dict:
        return {
            "n_records_listed": len(self.records_listed),
            "n_hea_found": len(self.hea_found),
            "n_dat_found": len(self.dat_found),
            "n_atr_found": len(self.atr_found),
            "n_jpg_found": len(self.jpg_found),
            "n_missing_hea_for_listed_records": len(self.missing_hea_for_records),
            "n_hea_not_in_records_file": len(self.hea_not_in_records),
            "n_incomplete_file_groups": len(self.incomplete_groups),
            "published_doc_n_recordings": project_config.PUBLISHED_DOC_N_RECORDINGS,
            "published_doc_n_annotated": project_config.PUBLISHED_DOC_N_ANNOTATED,
            "recordings_vs_published_doc_delta": len(self.hea_found) - project_config.PUBLISHED_DOC_N_RECORDINGS,
            "annotated_vs_published_doc_delta": len(self.atr_found) - project_config.PUBLISHED_DOC_N_ANNOTATED,
        }


def read_records_file(raw_dir: Path) -> list:
    records_path = raw_dir / "RECORDS"
    if not records_path.is_file():
        raise FileNotFoundError(
            f"RECORDS manifest not found at {records_path}. This is a "
            "dataset-fatal error -- Stage 00 will not substitute a glob of "
            "*.hea files for the dataset's own authoritative record list."
        )
    lines = [ln.strip() for ln in records_path.read_text().splitlines()]
    return [ln for ln in lines if ln]


def _find_duplicates(items: list) -> list:
    seen, dupes = set(), []
    for item in items:
        if item in seen:
            dupes.append(item)
        seen.add(item)
    return dupes


def discover_local_stems(raw_dir: Path, ext: str, recordings_only: bool = False) -> set:
    stems = {p.stem for p in raw_dir.glob(f"*{ext}")}
    if recordings_only:
        stems = {s for s in stems if _RECORDING_STEM_PREFIX_RE.match(s)}
    return stems


def reconcile(raw_dir: Path) -> ReconciliationResult:
    records_list = read_records_file(raw_dir)
    duplicates = _find_duplicates(records_list)
    if duplicates:
        raise ValueError(
            f"RECORDS manifest contains duplicate entries: {duplicates}. "
            "This is dataset-fatal -- duplicate identifiers prevent a "
            "unique recording-to-file mapping."
        )
    result = ReconciliationResult(
        records_listed=set(records_list),
        hea_found=discover_local_stems(raw_dir, ".hea"),
        dat_found=discover_local_stems(raw_dir, ".dat"),
        atr_found=discover_local_stems(raw_dir, ".atr"),
        jpg_found=discover_local_stems(raw_dir, ".jpg", recordings_only=True),
    )
    result.duplicate_in_records_file = duplicates
    return result


def parse_record_stem(stem: str) -> dict:
    m = re.match(project_config.RECORD_STEM_RE, stem)
    if not m:
        return {
            "recording_id": stem, "patient_id": None, "record_type": None,
            "record_number": None, "record_number_int": None,
            "record_number_suffix": None, "record_count": None, "parse_ok": False,
        }
    patient_id, rtype, num, count = m.groups()
    num_match = re.match(r"^(\d+)([a-z]?)$", num)
    num_int, num_suffix = num_match.groups()
    return {
        "recording_id": stem, "patient_id": patient_id,
        "record_type": {"p": "pregnancy", "l": "labour"}.get(rtype, rtype),
        "record_number": num, "record_number_int": int(num_int),
        "record_number_suffix": num_suffix or None, "record_count": int(count),
        "parse_ok": True,
    }


# =============================================================================
# Section 2 -- Checksums / .dat structural check
# (ported from pipeline/module0/checksums.py)
# =============================================================================

CHECKSUM_VERIFIED = "verified"
CHECKSUM_MISMATCH = "mismatch"
CHECKSUM_NOT_LISTED = "not_listed"


def parse_sha256sums(path: Path) -> dict:
    mapping = {}
    if not path.is_file():
        return mapping
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            raise ValueError(f"Unparseable line in {path.name}: {line!r}")
        digest, filename = parts
        mapping[filename.strip()] = digest.strip().lower()
    return mapping


def verify_checksum(path: Path, expected_hashes: dict) -> str:
    expected = expected_hashes.get(path.name)
    if expected is None:
        return CHECKSUM_NOT_LISTED
    actual = sha256_of_file(path)
    return CHECKSUM_VERIFIED if actual == expected else CHECKSUM_MISMATCH


def check_dat_structure(dat_path: Path, n_samples: int, n_signals: int, fmt: list) -> list:
    flags = []
    if not dat_path.is_file():
        flags.append("dat_file_missing")
        return flags
    bytes_per_sample = 2 if all(f == "16" for f in fmt) else None
    if bytes_per_sample is not None:
        expected_bytes = n_samples * n_signals * bytes_per_sample
        actual_bytes = dat_path.stat().st_size
        if actual_bytes != expected_bytes:
            flags.append(f"dat_size_mismatch:expected={expected_bytes}:actual={actual_bytes}")
    else:
        flags.append(f"unhandled_dat_format_for_size_check:{sorted(set(fmt))}")

    rec_stem = str(dat_path.with_suffix(""))
    try:
        wfdb.rdrecord(rec_stem, sampfrom=0, sampto=min(10, n_samples))
    except Exception as exc:  # noqa: BLE001 -- deliberately broad QC probe
        flags.append(f"dat_decode_failed_at_start:{exc}")
    try:
        start = max(0, n_samples - 10)
        wfdb.rdrecord(rec_stem, sampfrom=start, sampto=n_samples)
    except Exception as exc:  # noqa: BLE001
        flags.append(f"dat_decode_failed_at_end:{exc}")
    return flags


# =============================================================================
# Section 3 -- Header + clinical-info parsing
# (ported from pipeline/module0/headers.py)
# =============================================================================

_KEY_NORMALIZE = {
    "id": "patient_id",
    "record type": "record_type",
    "record number": "record_number_raw",
    "age(years)": "age_years",
    "bmi before pregnancy": "bmi_before_pregnancy",
    "bmi at recording": "bmi_at_recording",
    "gravidity": "gravidity",
    "parity": "parity",
    "previous caesarean": "previous_caesarean",
    "placental position": "placental_position",
    "gestational age at recording(w/d)": "gestational_age_at_recording_wd",
    "gestational age at delivery": "gestational_age_at_delivery_wd",
    "gestational age at delivery(w/d)": "gestational_age_at_delivery_wd",
    "mode of delivery": "mode_of_delivery",
    "synthetic oxytocin use in labour": "oxytocin_use",
    "epidural during labour": "epidural_use",
    "comments for recording": "comments_recording",
    "comments for delivery": "comments_delivery",
}
_INT_FIELDS = {"gravidity", "parity"}
_FLOAT_FIELDS = {"age_years", "bmi_before_pregnancy", "bmi_at_recording"}
_WD_FIELDS = {"gestational_age_at_recording_wd", "gestational_age_at_delivery_wd"}


def read_header(rec_path: Path) -> dict:
    record = wfdb.rdheader(str(rec_path))
    return {
        "fs_hz": float(record.fs),
        "n_samples": int(record.sig_len),
        "n_signals": int(record.n_sig),
        "sig_name_raw_order": list(record.sig_name),
        "fmt": list(record.fmt),
        "units": list(record.units),
        "adc_gain": list(record.adc_gain),
        "comments": list(record.comments),
        "duration_s": record.sig_len / record.fs if record.fs else None,
    }


def _parse_weeks_days(raw: str):
    if not raw:
        return None, None, None
    m = re.match(r"^\s*(\d+)\s*/\s*(\d+)\s*$", raw)
    if not m:
        return None, None, None
    weeks, days = int(m.group(1)), int(m.group(2))
    return weeks, days, weeks * 7 + days


def _safe_int(value: str):
    try:
        return int(value.strip())
    except (ValueError, AttributeError):
        return None


def _safe_float(value: str):
    try:
        return float(value.strip())
    except (ValueError, AttributeError):
        return None


def parse_clinical_info(comments: list) -> dict:
    raw_kv: dict = {}
    current_key = None
    flags = []
    for line in comments:
        if line.strip().lower() == "info:":
            continue
        key_candidate = None
        if ":" in line:
            candidate, _ = line.split(":", 1)
            key_candidate = candidate.strip().lower()
        if key_candidate is not None and key_candidate in _KEY_NORMALIZE:
            _, value = line.split(":", 1)
            current_key = key_candidate
            raw_kv[key_candidate] = value.strip()
            continue
        if current_key is None:
            flags.append(f"unattached_continuation_line:{line!r}")
            continue
        if ":" in line:
            flags.append(f"colon_in_continuation_line:{current_key}")
        existing = raw_kv.get(current_key, "")
        raw_kv[current_key] = (existing + " " + line.strip()).strip()

    out = {}
    for key_norm, value in raw_kv.items():
        field_name = _KEY_NORMALIZE.get(key_norm)
        if field_name is None:
            out.setdefault("_unrecognized_fields", {})[key_norm] = value
            flags.append(f"unrecognized_clinical_field:{key_norm}")
            continue
        out[f"{field_name}_raw"] = value
        if field_name in _INT_FIELDS:
            out[field_name] = _safe_int(value)
        elif field_name in _FLOAT_FIELDS:
            out[field_name] = _safe_float(value)
        elif field_name in _WD_FIELDS:
            weeks, days, total_days = _parse_weeks_days(value)
            out[f"{field_name}_weeks"] = weeks
            out[f"{field_name}_days_component"] = days
            out[f"{field_name}_total_days"] = total_days
        else:
            out[field_name] = value
    out["clinical_parse_flags"] = flags
    return out


# =============================================================================
# Section 4 -- Channel mapping (raw header order -> canonical EHG1..16)
# (ported from pipeline/module0/channel_mapping.py)
# =============================================================================

_LABEL_RE = re.compile(r"^EHG(\d+)$")


def build_channel_mapping(recording_id: str, sig_name_raw_order: list) -> pd.DataFrame:
    rows = []
    for raw_position, label in enumerate(sig_name_raw_order):
        m = _LABEL_RE.match(label)
        electrode_number = int(m.group(1)) if m else None
        rows.append({
            "recording_id": recording_id, "raw_position": raw_position,
            "header_label": label, "electrode_number": electrode_number,
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("electrode_number", na_position="last").reset_index(drop=True)
        df["canonical_position"] = range(len(df))
    return df


def validate_channel_mapping(mapping: pd.DataFrame) -> list:
    flags = []
    n = len(mapping)
    if n != project_config.EXPECTED_N_CHANNELS:
        flags.append(f"unexpected_channel_count:{n}")
    electrode_numbers = mapping["electrode_number"].tolist()
    if any(e is None for e in electrode_numbers):
        flags.append("unparseable_channel_label_present")
    valid_numbers = [e for e in electrode_numbers if e is not None]
    expected = set(range(1, project_config.EXPECTED_N_CHANNELS + 1))
    found = set(valid_numbers)
    missing = expected - found
    if missing:
        flags.append(f"missing_electrodes:{sorted(missing)}")
    extra = found - expected
    if extra:
        flags.append(f"unexpected_electrode_numbers:{sorted(extra)}")
    duplicates = [e for e in valid_numbers if valid_numbers.count(e) > 1]
    if duplicates:
        flags.append(f"duplicate_electrode_numbers:{sorted(set(duplicates))}")
    return flags


# =============================================================================
# Section 5 -- Annotation parsing
# (ported from pipeline/module0/annotations.py)
# =============================================================================

def parse_annotations(rec_path, recording_id: str, n_samples: int, fs_hz: float) -> pd.DataFrame:
    ann = wfdb.rdann(str(rec_path), "atr")
    samples = numpy.asarray(ann.sample)
    symbols = list(ann.symbol)
    aux_notes = list(ann.aux_note) if ann.aux_note is not None else [""] * len(symbols)

    n = len(samples)
    rows = []
    diffs = numpy.diff(samples) if n > 1 else numpy.array([])
    decreasing_at = set(numpy.where(diffs < 0)[0] + 1) if diffs.size else set()
    tie_at = set(numpy.where(diffs == 0)[0] + 1) if diffs.size else set()

    seen_pairs: dict = {}
    for i in range(n):
        sample = int(samples[i])
        symbol = symbols[i]
        aux_note = aux_notes[i]
        time_s = sample / fs_hz if fs_hz else None

        flags = []
        if i in decreasing_at:
            flags.append("timestamp_decreased_from_previous")
        if i in tie_at:
            flags.append("timestamp_tied_with_previous")

        pair_key = (sample, symbol)
        if pair_key in seen_pairs:
            flags.append(f"duplicate_of_index_{seen_pairs[pair_key]}")
        else:
            seen_pairs[pair_key] = i

        if n_samples is not None and (sample < 0 or sample >= n_samples):
            flags.append("out_of_recording_bounds")

        normalized = project_config.ANNOTATION_LABELS.get(symbol)
        is_recognized = normalized is not None
        if not is_recognized:
            flags.append("unrecognized_symbol")

        rows.append({
            "recording_id": recording_id, "annotation_index": i, "sample": sample,
            "time_s": time_s, "wfdb_symbol": symbol, "aux_note": aux_note,
            "normalized_event_type": normalized or "other",
            "is_recognized": is_recognized, "validation_flags": ";".join(flags),
        })

    return pd.DataFrame(rows, columns=[
        "recording_id", "annotation_index", "sample", "time_s", "wfdb_symbol",
        "aux_note", "normalized_event_type", "is_recognized", "validation_flags",
    ])


def summarize_annotation_flags(events: pd.DataFrame) -> list:
    if events.empty:
        return []
    flags = set()
    for f in events["validation_flags"]:
        if f:
            flags.update(f.split(";"))
    rolled = set()
    for f in flags:
        if f.startswith("duplicate_of_index_"):
            rolled.add("duplicate_annotation_present")
        elif f == "timestamp_decreased_from_previous":
            rolled.add("annotation_timestamps_decreasing")
        elif f == "timestamp_tied_with_previous":
            rolled.add("annotation_timestamps_tied")
        elif f == "out_of_recording_bounds":
            rolled.add("annotation_out_of_bounds")
        elif f == "unrecognized_symbol":
            rolled.add("unrecognized_annotation_symbol_present")
    return sorted(rolled)


# =============================================================================
# Section 6 -- Tocograph scan inspection
# (ported from pipeline/module0/toco.py)
# =============================================================================

_BLANK_STD_THRESHOLD = 5.0


def inspect_toco_image(jpg_path: Path) -> dict:
    result = {
        "toco_image_exists": jpg_path.is_file(), "toco_image_readable": False,
        "toco_image_width": None, "toco_image_height": None,
        "toco_image_filesize_bytes": None, "toco_image_pixel_std": None,
        "toco_image_blank_suspect": None, "toco_image_flags": [],
    }
    if not result["toco_image_exists"]:
        return result
    result["toco_image_filesize_bytes"] = jpg_path.stat().st_size
    try:
        with Image.open(jpg_path) as img:
            img.verify()
        with Image.open(jpg_path) as img:
            width, height = img.size
            arr = numpy.asarray(img.convert("L"), dtype=numpy.float32)
        result["toco_image_readable"] = True
        result["toco_image_width"] = width
        result["toco_image_height"] = height
        pixel_std = float(arr.std())
        result["toco_image_pixel_std"] = pixel_std
        result["toco_image_blank_suspect"] = pixel_std < _BLANK_STD_THRESHOLD
        if result["toco_image_blank_suspect"]:
            result["toco_image_flags"].append("low_pixel_variance_possible_blank_scan")
    except (UnidentifiedImageError, OSError) as exc:
        result["toco_image_flags"].append(f"unreadable_image:{exc}")
    return result


# =============================================================================
# Section 7 -- QC severity aggregation
# (ported from pipeline/module0/qc.py)
# =============================================================================

SEVERITY_RANK = {"info": 0, "warning": 1, "error": 2, "fatal": 3}

FLAG_SEVERITY = {
    "missing_hea": ("fatal", "quarantine"),
    "missing_dat": ("fatal", "quarantine"),
    "missing_jpg": ("warning", "inspect"),
    "dat_size_mismatch": ("error", "quarantine"),
    "dat_decode_failed_at_start": ("fatal", "quarantine"),
    "dat_decode_failed_at_end": ("error", "quarantine"),
    "unhandled_dat_format_for_size_check": ("warning", "inspect"),
    "checksum_mismatch_hea": ("error", "quarantine"),
    "checksum_mismatch_dat": ("error", "quarantine"),
    "checksum_mismatch_atr": ("error", "quarantine"),
    "checksum_mismatch_jpg": ("warning", "inspect"),
    "checksum_not_listed": ("info", "retain"),
    "record_stem_unparseable": ("error", "inspect"),
    "unexpected_channel_count": ("fatal", "quarantine"),
    "unparseable_channel_label_present": ("error", "quarantine"),
    "missing_electrodes": ("fatal", "quarantine"),
    "unexpected_electrode_numbers": ("error", "quarantine"),
    "duplicate_electrode_numbers": ("fatal", "quarantine"),
    "unexpected_fs": ("warning", "inspect"),
    "header_parse_failed": ("fatal", "quarantine"),
    "record_number_mismatch_header_vs_filename": ("warning", "inspect"),
    "clinical_parse_incomplete": ("info", "retain"),
    "unrecognized_clinical_field": ("warning", "inspect"),
    "colon_in_continuation_line": ("info", "retain"),
    "unattached_continuation_line": ("warning", "inspect"),
    "annotation_parse_failed": ("error", "quarantine"),
    "duplicate_annotation_present": ("warning", "inspect"),
    "annotation_timestamps_decreasing": ("error", "inspect"),
    "annotation_timestamps_tied": ("info", "retain"),
    "annotation_out_of_bounds": ("error", "inspect"),
    "unrecognized_annotation_symbol_present": ("warning", "inspect"),
    "unreadable_image": ("warning", "inspect"),
    "low_pixel_variance_possible_blank_scan": ("info", "inspect"),
}


@dataclass
class RecordingQC:
    recording_id: str
    flags: list = field(default_factory=list)

    def _flag_key(self, flag: str) -> str:
        return flag.split(":", 1)[0]

    def severity(self) -> str:
        if not self.flags:
            return "info"
        ranks = [SEVERITY_RANK[self._lookup(f)[0]] for f in self.flags]
        worst = max(ranks)
        return [k for k, v in SEVERITY_RANK.items() if v == worst][0]

    def recommended_action(self) -> str:
        if not self.flags:
            return "retain"
        actions = [self._lookup(f)[1] for f in self.flags]
        if "quarantine" in actions:
            return "quarantine"
        if "inspect" in actions:
            return "inspect"
        return "retain"

    def _lookup(self, flag: str):
        return FLAG_SEVERITY.get(self._flag_key(flag), ("warning", "inspect"))

    def to_row(self) -> dict:
        return {
            "recording_id": self.recording_id, "qc_severity": self.severity(),
            "recommended_action": self.recommended_action(),
            "n_flags": len(self.flags), "flags": ";".join(self.flags),
        }


# =============================================================================
# Section 8 -- Provenance / run metadata
# (ported from pipeline/module0/provenance.py)
# =============================================================================

def get_git_commit(repo_root: Path):
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


def build_run_metadata(raw_dir: Path, project_root: Path, counts: dict) -> dict:
    return {
        "dataset_name": "Icelandic 16-electrode Electrohysterogram Database (ice16ehgdb)",
        "dataset_source": "https://physionet.org/content/ehgdb/1.0.0/",
        "pipeline": "paper1_reproducible_pipeline (standalone; see Paper1/MIGRATION_PLAN.md)",
        "content_hashes": {
            "records_file_sha256": sha256_of_file(project_config.RECORDS_FILE) if project_config.RECORDS_FILE.is_file() else None,
            "sha256sums_file_sha256": sha256_of_file(project_config.SHA256SUMS_FILE) if project_config.SHA256SUMS_FILE.is_file() else None,
            "raw_dir_listing_hash": hash_directory_listing(raw_dir),
        },
        "code_version": {
            "git_commit": get_git_commit(project_root),
            "stage00_schema_version": "1.0.0",
        },
        "library_versions": {
            "wfdb": wfdb.__version__,
            "pandas": pd.__version__,
            "numpy": numpy.__version__,
        },
        "execution_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "counts": counts,
    }


# =============================================================================
# Section 9 -- Per-recording orchestration
# (ported from pipeline/module0/run_module0.py)
# =============================================================================

def _process_one_recording(stem: str, reconciliation, sha_map: dict, raw_dir: Path):
    qc = RecordingQC(recording_id=stem)
    stem_info = parse_record_stem(stem)
    if not stem_info["parse_ok"]:
        qc.flags.append("record_stem_unparseable")

    has_hea = stem in reconciliation.hea_found
    has_dat = stem in reconciliation.dat_found
    has_atr = stem in reconciliation.atr_found
    has_jpg = stem in reconciliation.jpg_found

    if not has_hea:
        qc.flags.append("missing_hea")
    if not has_dat:
        qc.flags.append("missing_dat")
    if not has_jpg:
        qc.flags.append("missing_jpg")

    rec_path = raw_dir / stem
    header = None
    channel_mapping_df = pd.DataFrame()
    clinical = {}
    events_df = pd.DataFrame()

    if has_hea:
        try:
            header = read_header(rec_path)
        except Exception:  # noqa: BLE001
            qc.flags.append(f"header_parse_failed:{traceback.format_exc(limit=1).strip()[-200:]}")

    if header is not None:
        if abs(header["fs_hz"] - project_config.EXPECTED_FS_HZ) > 1e-6:
            qc.flags.append(f"unexpected_fs:{header['fs_hz']}")

        channel_mapping_df = build_channel_mapping(stem, header["sig_name_raw_order"])
        qc.flags.extend(validate_channel_mapping(channel_mapping_df))

        clinical = parse_clinical_info(header["comments"])
        if clinical.get("clinical_parse_flags"):
            qc.flags.extend(clinical["clinical_parse_flags"])

        rn_raw = clinical.get("record_number_raw_raw")
        if rn_raw and "/" in rn_raw and stem_info["parse_ok"]:
            num_part, _, count_part = rn_raw.partition("/")
            if num_part.strip() != stem_info["record_number"] or count_part.strip() != str(stem_info["record_count"]):
                qc.flags.append("record_number_mismatch_header_vs_filename")

        if has_dat:
            qc.flags.extend(check_dat_structure(
                rec_path.with_suffix(".dat"), header["n_samples"], header["n_signals"], header["fmt"],
            ))

        if has_atr:
            try:
                events_df = parse_annotations(rec_path, stem, header["n_samples"], header["fs_hz"])
                qc.flags.extend(summarize_annotation_flags(events_df))
            except Exception:  # noqa: BLE001
                qc.flags.append(f"annotation_parse_failed:{traceback.format_exc(limit=1).strip()[-200:]}")

    for ext, present in (("hea", has_hea), ("dat", has_dat), ("atr", has_atr), ("jpg", has_jpg)):
        if not present:
            continue
        path = raw_dir / f"{stem}.{ext}"
        status = verify_checksum(path, sha_map)
        if status == CHECKSUM_MISMATCH:
            qc.flags.append(f"checksum_mismatch_{ext}")
        elif status == CHECKSUM_NOT_LISTED:
            qc.flags.append("checksum_not_listed")

    toco_info = inspect_toco_image(raw_dir / f"{stem}.jpg")
    for f in toco_info["toco_image_flags"]:
        qc.flags.append(f)

    inventory_row = {
        **stem_info, "has_hea": has_hea, "has_dat": has_dat, "has_atr": has_atr, "has_jpg": has_jpg,
        "fs_hz": header["fs_hz"] if header else None,
        "n_samples": header["n_samples"] if header else None,
        "n_signals": header["n_signals"] if header else None,
        "duration_s": header["duration_s"] if header else None,
        **{k: v for k, v in toco_info.items() if k != "toco_image_flags"},
    }
    clinical_row = {"recording_id": stem, **{k: v for k, v in clinical.items() if k != "clinical_parse_flags"}}

    return {
        "inventory_row": inventory_row, "clinical_row": clinical_row,
        "channel_mapping_df": channel_mapping_df, "events_df": events_df, "qc": qc,
    }


def run(raw_dir: Path = None, annotations_dir: Path = None, results_dir: Path = None) -> int:
    raw_dir = raw_dir or project_config.RAW_DIR
    annotations_dir = annotations_dir or project_config.ANNOTATIONS_DIR
    results_dir = results_dir or (project_config.AUDIT_DIR)
    annotations_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    print(f"[00_data_audit] raw_dir={raw_dir}")
    if not raw_dir.is_dir():
        print(f"[00_data_audit] FATAL: raw_dir does not exist: {raw_dir}", file=sys.stderr)
        print(
            "[00_data_audit] See data/raw/README.md -- set EHG_RAW_DIR or "
            "place the PhysioNet recordings under data/raw/.",
            file=sys.stderr,
        )
        return 2

    reconciliation = fi_reconcile_or_raise(raw_dir)
    sha_map = parse_sha256sums(project_config.SHA256SUMS_FILE)

    all_stems = sorted(reconciliation.records_listed | reconciliation.hea_found)
    print(f"[00_data_audit] auditing {len(all_stems)} recording stems...")

    inventory_rows, clinical_rows, qc_rows = [], [], []
    channel_mapping_frames, events_frames = [], []

    for i, stem in enumerate(all_stems, 1):
        result = _process_one_recording(stem, reconciliation, sha_map, raw_dir)
        inventory_rows.append(result["inventory_row"])
        clinical_rows.append(result["clinical_row"])
        qc_rows.append(result["qc"].to_row())
        if not result["channel_mapping_df"].empty:
            channel_mapping_frames.append(result["channel_mapping_df"])
        if not result["events_df"].empty:
            events_frames.append(result["events_df"])
        if i % 20 == 0 or i == len(all_stems):
            print(f"[00_data_audit]   {i}/{len(all_stems)} recordings processed")

    recording_inventory = pd.DataFrame(inventory_rows)
    clinical_metadata = pd.DataFrame(clinical_rows)
    recording_qc_flags = pd.DataFrame(qc_rows)
    channel_mapping = pd.concat(channel_mapping_frames, ignore_index=True) if channel_mapping_frames else pd.DataFrame()
    annotation_events = pd.concat(events_frames, ignore_index=True) if events_frames else pd.DataFrame()

    annotation_counts = pd.DataFrame()
    if not annotation_events.empty:
        annotation_counts = (
            annotation_events.groupby(["recording_id", "normalized_event_type"])
            .size().reset_index(name="count")
        )

    # --- Write outputs -------------------------------------------------
    # data/annotations/ -- reusable inputs for later stages
    channel_mapping.to_csv(annotations_dir / "channel_mapping.csv", index=False)
    annotation_counts.to_csv(annotations_dir / "annotation_counts.csv", index=False)
    if not annotation_events.empty:
        annotation_events.to_csv(annotations_dir / "annotation_events.csv", index=False)
        try:
            annotation_events.to_parquet(annotations_dir / "annotation_events.parquet", index=False)
        except ImportError:
            print("[00_data_audit] WARNING: pyarrow not installed, skipping .parquet output", file=sys.stderr)
    else:
        # Always write the file (even empty) so downstream stages have a
        # stable path to read from and fail loudly/clearly if genuinely empty.
        pd.DataFrame(columns=[
            "recording_id", "annotation_index", "sample", "time_s", "wfdb_symbol",
            "aux_note", "normalized_event_type", "is_recognized", "validation_flags",
        ]).to_csv(annotations_dir / "annotation_events.csv", index=False)

    # results/paper1/00_data_audit/ -- audit-specific outputs
    recording_inventory.to_csv(results_dir / "recording_inventory.csv", index=False)
    clinical_metadata.to_csv(results_dir / "clinical_metadata.csv", index=False)
    recording_qc_flags.to_csv(results_dir / "recording_qc_flags.csv", index=False)

    valid_patients = recording_inventory.loc[recording_inventory["parse_ok"], "patient_id"]
    severity_counts = recording_qc_flags["qc_severity"].value_counts().to_dict() if not recording_qc_flags.empty else {}
    action_counts = recording_qc_flags["recommended_action"].value_counts().to_dict() if not recording_qc_flags.empty else {}
    symbol_counts = annotation_events["wfdb_symbol"].value_counts().to_dict() if not annotation_events.empty else {}

    dataset_summary = {
        "reconciliation": reconciliation.summary(),
        "n_recordings_audited": len(all_stems),
        "n_patients": int(valid_patients.nunique()),
        "n_recordings_with_annotations": int(recording_inventory["has_atr"].sum()),
        "annotation_symbol_counts": symbol_counts,
        "qc_severity_counts": severity_counts,
        "qc_recommended_action_counts": action_counts,
        "expected_n_patients": project_config.EXPECTED_N_PATIENTS,
        "patients_vs_expected_delta": int(valid_patients.nunique()) - project_config.EXPECTED_N_PATIENTS,
    }
    (results_dir / "dataset_summary.json").write_text(json.dumps(dataset_summary, indent=2, default=str))

    run_metadata = build_run_metadata(raw_dir, project_config.PROJECT_ROOT, dataset_summary)
    (results_dir / "run_metadata.json").write_text(json.dumps(run_metadata, indent=2, default=str))

    _write_cohort_summary(results_dir / "cohort_summary.md", dataset_summary, reconciliation)

    n_fatal_or_error = sum(1 for r in qc_rows if r["qc_severity"] in ("fatal", "error"))
    print(f"[00_data_audit] done. {n_fatal_or_error} recording(s) flagged error/fatal severity.")
    print(f"[00_data_audit] data outputs written to {annotations_dir}")
    print(f"[00_data_audit] audit outputs written to {results_dir}")
    return 1 if n_fatal_or_error > 0 else 0


def fi_reconcile_or_raise(raw_dir: Path) -> ReconciliationResult:
    """Thin wrapper name kept distinct from `reconcile` so the dataset-fatal
    vs. recording-level failure boundary reads clearly at the call site."""
    return reconcile(raw_dir)


def _write_cohort_summary(path: Path, summary: dict, reconciliation) -> None:
    r = summary["reconciliation"]
    lines = [
        "# Paper 1 Reproducible Pipeline -- Stage 00 Cohort Summary",
        "",
        "All counts below are computed from the raw dataset by this run -- "
        "none are hardcoded.",
        "",
        f"- Recordings listed in RECORDS manifest: {r['n_records_listed']}",
        f"- Recordings with a local .hea file: {r['n_hea_found']}",
        f"- Recordings with a local .dat file: {r['n_dat_found']}",
        f"- Recordings with a local .atr (annotation) file: {r['n_atr_found']}",
        f"- Recordings with a local .jpg (tocograph scan): {r['n_jpg_found']}",
        f"- Patients (parsed from recording IDs): {summary['n_patients']} "
        f"(expected {summary['expected_n_patients']}, delta {summary['patients_vs_expected_delta']})",
        "",
        "## Local files vs. PhysioNet's published project-page text",
        "",
        f"- Published doc states {r['published_doc_n_recordings']} recordings; "
        f"local RECORDS/.hea count is {r['n_hea_found']} "
        f"(delta {r['recordings_vs_published_doc_delta']:+d}).",
        f"- Published doc states {r['published_doc_n_annotated']} annotated recordings; "
        f"local .atr count is {r['n_atr_found']} "
        f"(delta {r['annotated_vs_published_doc_delta']:+d}).",
        "",
        "## QC severity distribution",
        "",
    ]
    for sev in ("info", "warning", "error", "fatal"):
        lines.append(f"- {sev}: {summary['qc_severity_counts'].get(sev, 0)}")
    lines += ["", "## Recommended-action distribution", ""]
    for act in ("retain", "inspect", "quarantine"):
        lines.append(f"- {act}: {summary['qc_recommended_action_counts'].get(act, 0)}")
    lines += ["", "## Annotation symbol counts (raw, across all recordings)", ""]
    for sym, count in sorted(summary["annotation_symbol_counts"].items(), key=lambda kv: -kv[1]):
        lines.append(f"- `{sym}`: {count}")
    if r["n_missing_hea_for_listed_records"] or r["n_hea_not_in_records_file"]:
        lines += ["", "## File-inventory anomalies", ""]
        if r["n_missing_hea_for_listed_records"]:
            lines.append(
                f"- {r['n_missing_hea_for_listed_records']} recording(s) listed in RECORDS "
                "have no local .hea file."
            )
        if r["n_hea_not_in_records_file"]:
            lines.append(
                f"- {r['n_hea_not_in_records_file']} local .hea file(s) are not listed in RECORDS."
            )
    lines += [
        "",
        "## Reproducibility note",
        "",
        "This stage's EHG9-relevant outputs (recording_inventory.csv, "
        "channel_mapping.csv, annotation_events.csv) are validated against "
        "the frozen pipeline/module0/outputs/ files by the reproducibility "
        "gate before Stage 01 onward is trusted. See "
        "results/paper1/validation/.",
    ]
    path.write_text("\n".join(lines) + "\n")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--raw-dir", type=Path, default=None, help="Override raw data directory (else EHG_RAW_DIR env var, else data/raw/)")
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    raw_dir = args.raw_dir.expanduser().resolve() if args.raw_dir else None
    sys.exit(run(raw_dir=raw_dir))
