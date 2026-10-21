"""Regression fingerprints for the committed manuscript results.

These do NOT recompute anything -- they need no raw data and run in CI. They
exist so that an unrelated code change cannot silently alter a tracked
results file without a test failing. If a change is supposed to alter a
results file, regenerate ``tests/fixtures/results_fingerprint.json`` (see the
comment at the bottom of this file) as a deliberate, reviewed step -- never as
a side effect of making this test pass.

Figure 1 and Table 1 get an additional, human-readable numeric spot-check
(a handful of specific values, not a hash) so a failure here says what
actually moved.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from gi_sampen import paths

FIXTURES_DIR = Path(__file__).resolve().parent / "fixtures"
FINGERPRINT_PATH = FIXTURES_DIR / "results_fingerprint.json"


def _fingerprint() -> dict:
    return json.loads(FINGERPRINT_PATH.read_text())


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.mark.parametrize("relative_path,expected", list(_fingerprint().items()))
def test_committed_results_unchanged(relative_path, expected):
    path = paths.REPO_ROOT / relative_path
    assert path.is_file(), f"missing tracked result: {relative_path}"
    actual_sha256 = _sha256(path)
    assert actual_sha256 == expected["sha256"], (
        f"{relative_path} content changed.\n"
        f"  expected sha256: {expected['sha256']}\n"
        f"  actual sha256:   {actual_sha256}\n"
        "If this change is intentional, regenerate "
        "tests/fixtures/results_fingerprint.json as a deliberate step -- see "
        "the comment at the bottom of this file -- and say so in the commit "
        "message."
    )
    if "n_rows" in expected:
        with open(path, newline="") as fh:
            rows = list(csv.reader(fh))
        n_rows = len(rows) - 1 if rows else 0
        assert n_rows == expected["n_rows"], (
            f"{relative_path}: row count changed "
            f"({expected['n_rows']} -> {n_rows})"
        )
        header = rows[0] if rows else []
        assert header == expected["columns"], (
            f"{relative_path}: columns changed"
        )


def test_table1_source_spot_check():
    """Table 1 (the four-channel FM-vs-baseline cross-channel result):
    the primary feature's EHG9 estimate and q-value, read from the tracked
    CSV, not recomputed."""
    path = paths.REPO_ROOT / "results" / "gi" / "tables" / "cross_channel_replication_summary.csv"
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))

    row = next(
        r for r in rows
        if r["feature"] == "gi_sd"
        and r["channel"] == "EHG9"
        and r["contrast"] == "fetal_movement_vs_baseline"
    )
    assert float(row["estimate"]) == pytest.approx(-0.2460081397964411, rel=1e-9)
    assert float(row["q_value"]) == pytest.approx(1.8834442379185916e-07, rel=1e-6)
    assert row["validated"] == "True"


def test_gate_summary_decision_spot_check():
    """The go/no-go gate's decision and its three gate features, read from
    the tracked CSV. See docs/provenance.md for the decision rule."""
    path = paths.REPO_ROOT / "results" / "gi" / "tables" / "gate_summary.csv"
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))

    fm_rows = {
        r["feature"]: r for r in rows
        if r["contrast"] == "fetal_movement_vs_baseline"
    }
    gate_features = ("gi_sd", "gi_kurtosis", "gi_skewness")
    expected_signs = {"gi_sd": -1, "gi_kurtosis": 1, "gi_skewness": 1}

    for feature in gate_features:
        assert feature in fm_rows, f"gate feature {feature} missing from gate_summary.csv"
        row = fm_rows[feature]
        assert row["decision"] == "GO"
        estimate = float(row["estimate"])
        assert (estimate > 0) == (expected_signs[feature] > 0), (
            f"{feature}: sign flipped (estimate={estimate})"
        )


# ---------------------------------------------------------------------------
# To regenerate the fingerprint file after a deliberate, reviewed change to a
# tracked results file:
#
#   python - <<'PY'
#   import csv, hashlib, json, pathlib
#   root = pathlib.Path("results/gi")
#   manifest = {}
#   for path in sorted(root.rglob("*.csv")):
#       with open(path, newline="") as fh:
#           rows = list(csv.reader(fh))
#       header = rows[0] if rows else []
#       manifest[str(path)] = {
#           "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
#           "n_rows": len(rows) - 1 if rows else 0,
#           "n_cols": len(header),
#           "columns": header,
#       }
#   for path in sorted(root.rglob("*.json")):
#       manifest[str(path)] = {"sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
#   out = pathlib.Path("tests/fixtures/results_fingerprint.json")
#   out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
#   PY
# ---------------------------------------------------------------------------
