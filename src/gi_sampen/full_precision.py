#!/usr/bin/env python3
"""Section 8 (optional) -- full-precision (unrounded) profile sanity check.

Tests the mathematical inference from the provenance audit: that an
UNROUNDED, data-derived tolerance axis is exactly gain-equivariant by
construction, and that gain-dependence is introduced specifically by
fixed-absolute-decimal rounding, not by the profiling method itself.

This is a SEPARATE experimental code path. It does not replace or modify
frozen_summaries or gi_summaries_rounded (gi_profile.py is untouched). The
only machinery reused from gi_profile.py is _validate/_aligned_templates
(pure helpers, no rounding inside them) and the module-level constants
(embedding order, epsilon). Pair counting uses the same cKDTree
count_neighbors approach as the rest of Paper 2, but at the window's own
UNIQUE RAW (unrounded, float64) pairwise Chebyshev distances as the
tolerance axis -- i.e. what core/sampleEntropy_Gayathri.py's own docstring
describes ("r_range: unique sorted r-values derived from the data") before
that file's actual code rounds to 3 decimals.

Cost note: O(N^2) pairwise distances per window (pdist), N ~ 1198 templates
for a 60s/20Hz window -> ~1.4M pairs, cheap. To keep total runtime low this
defaults to a SMALL subset of the existing 60-window gain-stress set (first
4 windows per class = 12 windows); override with --n-windows. Use --probe
first to time it before committing to a larger subset.
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist

from . import paths as _common
from . import gi_profile
from .gain_stress import select_windows

GAINS = (0.25, 0.5, 1.3, 2.0, 4.0, 7.1)
DEFAULT_N_PER_CLASS = 4  # -> 12 windows total across UC/IN/FM
EPSILON = gi_profile.PROFILE_EPSILON
ORDER = gi_profile.EMBEDDING_DIMENSION
# Looser than the rounded-variant invariance test (1e-10 / 1e-8): cdist +
# independent per-scale unique-value extraction accumulates more
# floating-point noise than the rounded/fixed-grid variants, which never
# form an explicit distance matrix. This threshold is for this script only
# and does not alter invariance_check.py's definition.
INVARIANCE_ATOL = 1e-8
INVARIANCE_RTOL = 1e-6


def _unique_raw_distances(templates: np.ndarray) -> np.ndarray:
    """Unrounded, unique, off-diagonal pairwise Chebyshev distances."""
    condensed = pdist(templates, metric="chebyshev")
    return np.unique(condensed)


def full_precision_summaries(signal, *, order: int = ORDER, epsilon: float = EPSILON) -> dict:
    data = gi_profile._validate(signal)
    nan_out = {"fp_total": float("nan"), "fp_average": float("nan"),
              "fp_sd": float("nan"), "fp_n_retained": float("nan")}
    if data.size == 0 or data.size < order + 2:
        return nan_out

    template_count = data.size - order
    templates_m = gi_profile._aligned_templates(data, order, template_count)
    templates_m1 = gi_profile._aligned_templates(data, order + 1, template_count)

    d_m = _unique_raw_distances(templates_m)
    d_m1 = _unique_raw_distances(templates_m1)
    taus = np.unique(np.concatenate([d_m, d_m1]))  # union axis, as chm() does

    tree_m = cKDTree(templates_m)
    tree_m1 = cKDTree(templates_m1)

    def cumulative_counts(tree: cKDTree) -> np.ndarray:
        shells = np.asarray(
            tree.count_neighbors(tree, taus, p=np.inf, cumulative=False),
            dtype=np.float64)
        shells[0] -= template_count  # remove ordered self-pairs at distance 0
        return np.cumsum(shells)

    denominator = float(template_count * (template_count - 1))
    b = cumulative_counts(tree_m) / denominator
    a = cumulative_counts(tree_m1) / denominator

    profile = np.log((b + epsilon) / (a + epsilon))
    valid = np.isfinite(profile) & (b > epsilon) & (a > epsilon)
    values = profile[valid]
    if values.size == 0:
        return {**nan_out, "fp_n_retained": 0.0}
    return {
        "fp_total": float(np.sum(values)),
        "fp_average": float(np.mean(values)),
        "fp_sd": float(np.std(values)),
        "fp_n_retained": float(values.size),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-per-class", type=int, default=DEFAULT_N_PER_CLASS)
    parser.add_argument("--probe", action="store_true",
                        help="time 1 window across all 7 gains and stop")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    frame = _common.load_analysis_windows()
    selected = select_windows(frame)
    per_class = selected.groupby("class_label", group_keys=False).apply(
        lambda g: g.head(args.n_per_class))
    if args.probe:
        per_class = per_class.head(1)

    print(f"{len(per_class)} windows x {len(GAINS) + 1} gains "
          f"(full-precision, no rounding)", flush=True)

    started = time.perf_counter()
    rows = []
    for row, signal in _common.iter_windows(per_class, progress=not args.quiet):
        for gain in (1.0,) + GAINS:
            summaries = full_precision_summaries(gain * signal)
            rows.append({
                "window_id": row.window_id, "state": row.window_class,
                "patient_id": row.patient_id, "gain": gain, **summaries,
            })
    elapsed = time.perf_counter() - started

    if args.probe:
        print(f"PROBE: 1 window x {len(GAINS)+1} gains in {elapsed:.1f}s")
        print(f"  projected for {len(selected)} windows (full 60-window set): "
              f"{elapsed*len(selected)/60:.1f} min")
        print(f"  projected for the default {DEFAULT_N_PER_CLASS*3}-window subset: "
              f"{elapsed*DEFAULT_N_PER_CLASS*3/60:.1f} min")
        return 0

    long = pd.DataFrame(rows)
    _common.write_csv(long, _common.RESULTS_DIR / "full_precision_gain_test.csv",
                      label="full-precision gain test")

    base = long.loc[long["gain"] == 1.0].set_index("window_id")
    lines = ["# Full-precision (unrounded) profile -- gain sanity check\n",
             f"{len(per_class)} windows (subset of the 60-window gain-stress set), "
             f"gains {GAINS}, no rounding at any step.\n",
             "| feature | gain | median |diff| | max |diff| | max relative |diff| |",
             "|---|---:|---:|---:|---:|"]
    verdict_ok = True
    for feature in ("fp_total", "fp_average", "fp_sd"):
        baseline = base[feature]
        for gain in GAINS:
            subset = long.loc[long["gain"] == gain].set_index("window_id")
            want = baseline.reindex(subset.index)
            got = subset[feature]
            diff = (got - want).abs()
            with np.errstate(divide="ignore", invalid="ignore"):
                rel = (diff / want.abs().where(want.abs() > 1e-10))
            passed = np.isclose(got, want, atol=INVARIANCE_ATOL, rtol=INVARIANCE_RTOL)
            if not bool(passed.all()):
                verdict_ok = False
            lines.append(f"| {feature} | {gain:.2f} | {diff.median():.3e} | "
                        f"{diff.max():.3e} | {rel.max():.3e} |")

    lines.append("")
    lines.append(f"Verdict: {'PASS' if verdict_ok else 'FAIL'} -- "
                f"{'all' if verdict_ok else 'not all'} full-precision summaries "
                f"agreed with the 1x baseline within atol={INVARIANCE_ATOL}, "
                f"rtol={INVARIANCE_RTOL} at every tested gain.")
    lines.append("")
    lines.append("This is a mathematical-inference check on a small subset, not a "
                "scientific claim on its own -- read alongside Section 4/5 of the "
                "provenance audit and the gain-resolution sensitivity report.")

    report_path = _common.AUDIT_DIR / "full_precision_gain_test_report.md"
    report_path.write_text("\n".join(lines) + "\n")
    print(f"elapsed: {elapsed/60:.1f} min")
    print(f"wrote {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
