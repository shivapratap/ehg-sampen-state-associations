#!/usr/bin/env python3
"""Stage D -- apply the decision rule and print the one table to paste back.

Decision rule, fixed before any result was seen
-----------------------------------------------
Gate features : gi_sd, gi_kurtosis, gi_skewness   (FM vs IN, EHG9/FWH)
Expected signs: negative, positive, positive      (frozen directions)

  GO     at least 2 of the 3 have q < 0.05 AND the same sign as frozen
  NO-GO  0 or 1 do

Precondition: invariance_check.py must have passed. An invariance failure is a
bug, not a result, and voids the gate.

Recorded but NOT gating: UC vs IN on the GI features; the GI-AUC/shape
descriptors on the fixed grid; frozen-vs-GI effect-size deltas.

Output: results/gate_summary.csv -- small by design, ~30 rows. That file is
the only thing that needs to leave this machine.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import paths as _common
from .statistics import GATE_FEATURES

EXPECTED_SIGN = {"gi_sd": -1.0, "gi_kurtosis": 1.0, "gi_skewness": 1.0}


def check_invariance() -> tuple[bool, str]:
    path = _common.RESULTS_DIR / "invariance_summary.csv"
    if not path.exists():
        return False, "invariance_summary.csv not found -- Stage B not run"
    summary = pd.read_csv(path)
    gi = summary.loc[summary["feature"].str.startswith("gi")]
    if gi.empty:
        return False, "no GI features in the invariance summary"
    worst = float(gi["proportion_invariant"].min())
    if worst < 1.0:
        failing = sorted(set(
            gi.loc[gi["proportion_invariant"] < 1.0, "feature"]))
        return False, f"not invariant: {', '.join(failing)}"
    return True, "all GI features numerically invariant at every scale"


def main() -> int:
    results = pd.read_csv(_common.RESULTS_DIR / "gi_model_results.csv")
    invariance_ok, invariance_note = check_invariance()

    grid_note = ""
    grid_path = _common.RESULTS_DIR / "grid_definition.json"
    if grid_path.exists():
        definition = json.loads(grid_path.read_text())
        grid_note = (f"tau in [{definition['tau_min']:.3f}, "
                     f"{definition['tau_max']:.3f}] sd, "
                     f"{definition['n_grid_points']} points")

    fm = results.loc[results["contrast"] == "fetal_movement_vs_baseline"]
    uc = results.loc[results["contrast"] == "contraction_vs_baseline"]

    print("=" * 76)
    print("ICASSP GO/NO-GO GATE -- EHG9 / FWH / 60 s")
    print("=" * 76)
    print(f"invariance precondition : "
          f"{'PASS' if invariance_ok else 'FAIL'}  ({invariance_note})")
    if grid_note:
        print(f"fixed grid              : {grid_note}")
    print()

    print("GATE -- FM versus IN")
    print(f"{'feature':<16}{'GI est':>9}{'GI q':>10}{'frozen est':>12}"
          f"{'frozen q':>10}{'sign':>7}{'q<.05':>7}")
    passes = 0
    for feature in GATE_FEATURES:
        row = fm.loc[fm["feature"] == feature]
        if row.empty:
            print(f"{feature:<16}{'missing':>9}")
            continue
        row = row.iloc[0]
        estimate = float(row["estimate"])
        q = float(row["q_value"])
        sign_ok = np.isfinite(estimate) and np.sign(estimate) == EXPECTED_SIGN[feature]
        q_ok = np.isfinite(q) and q < 0.05
        if sign_ok and q_ok:
            passes += 1
        print(f"{feature:<16}{estimate:>9.3f}{q:>10.4f}"
              f"{float(row['frozen_estimate']):>12.3f}"
              f"{float(row['frozen_q_value']):>10.4f}"
              f"{'ok' if sign_ok else 'FLIP':>7}"
              f"{'yes' if q_ok else 'no':>7}")

    print()
    print("REST OF THE GI FAMILY -- FM versus IN (recorded, not gating)")
    print(f"{'feature':<18}{'est':>9}{'q':>10}{'model':>18}")
    for _, row in fm.iterrows():
        if row["feature"] in GATE_FEATURES:
            continue
        print(f"{row['feature']:<18}{float(row['estimate']):>9.3f}"
              f"{float(row['q_value']):>10.4f}{str(row['model']):>18}")

    print()
    print("UC versus IN (the contrast -- weakening here is the expected result)")
    print(f"{'feature':<18}{'GI est':>9}{'GI q':>10}{'frozen est':>12}{'frozen q':>10}")
    for _, row in uc.iterrows():
        print(f"{row['feature']:<18}{float(row['estimate']):>9.3f}"
              f"{float(row['q_value']):>10.4f}"
              f"{float(row['frozen_estimate']):>12.3f}"
              f"{float(row['frozen_q_value']):>10.4f}")

    decision = "GO" if (invariance_ok and passes >= 2) else "NO-GO"
    print()
    print("=" * 76)
    print(f"{passes}/3 gate features retained direction and FDR support")
    print(f"DECISION: {decision}")
    if decision == "NO-GO" and invariance_ok:
        print("  The FM entropy signature does not survive gain normalisation.")
        print("  Stop the ICASSP submission. This belongs in the Access")
        print("  supplement as a limitation -- it is a real finding.")
    elif decision == "NO-GO":
        print("  Blocked by the invariance precondition: fix the bug, re-run.")
    else:
        print("  The FM signature survives. Write the paper.")
    print("=" * 76)

    export = pd.concat([fm, uc])[[
        "contrast", "feature", "n", "n_women", "estimate", "std_error",
        "ci_low", "ci_high", "p_value", "q_value", "frozen_estimate",
        "frozen_q_value", "model", "fallback", "boundary_fit",
        "warning_summary"]].copy()
    export.insert(0, "decision", decision)
    export.insert(1, "invariance_pass", invariance_ok)
    _common.write_csv(export, _common.RESULTS_DIR / "gate_summary.csv",
                      label="gate summary (paste this one back)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
