#!/usr/bin/env python3
"""Stage C -- patient-aware mixed models on the GI features.

The model specification is copied from Paper 1 `src/models.py::_fit_one` and
must not drift from it, because the whole point of the gate is a like-for-like
comparison against the frozen estimates:

    z-standardised feature
    ~ C(window_class, Treatment(reference="baseline"))
    groups = patient_id, re_formula = "1"
    vc_formula = {"recording": "0 + C(recording_id)"}      (nested tier)
    REML, optimizers tried in order: lbfgs, bfgs, maxiter 1000
    fall back to woman-only random effects if the nested tier fails
    boundary/singularity warnings recorded, not silently swallowed

BH-FDR is applied WITHIN THE GI FAMILY and within each contrast. The GI family
is its own family and is deliberately NOT merged into Paper 1's frozen 34, for
two reasons: these are new features, and merging would retrospectively change
the frozen q-values.
"""
from __future__ import annotations

import argparse
import warnings
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

from . import paths as _common

# --- Frozen Paper 1 model constants ----------------------------------------
OPTIMIZER_METHODS = ("lbfgs", "bfgs")
VARIANCE_BOUNDARY_THRESHOLD = 1e-8
MIN_BASELINE_N = 5
MIN_CONTRACTION_N = 5
FORMULA = '_z ~ C(window_class, Treatment(reference="baseline"))'
FM_KEY = 'C(window_class, Treatment(reference="baseline"))[T.fetal_movement]'
UC_KEY = 'C(window_class, Treatment(reference="baseline"))[T.contraction]'

# --- The GI family ---------------------------------------------------------
# Order matters only for readability. The three gate features come first.
GI_FAMILY = (
    "gi_sd", "gi_kurtosis", "gi_skewness",          # <-- the gate
    "gi_average", "gi_median", "gi_total",
    "gi_fixed_r_sampen",
    "gigrid_auc", "gigrid_slope", "gigrid_centroid",
    "gigrid_sd", "gigrid_kurtosis", "gigrid_skewness",
)
GATE_FEATURES = ("gi_sd", "gi_kurtosis", "gi_skewness")

# Frozen counterparts, for the side-by-side column.
FROZEN_COUNTERPART = {
    "gi_sd": "fwh_entropy_sdsampen",
    "gi_kurtosis": "fwh_entropy_kurtosissampen",
    "gi_skewness": "fwh_entropy_skewnesssampen",
    "gi_average": "fwh_entropy_averagesampen",
    "gi_median": "fwh_entropy_mediansampen",
    "gi_total": "fwh_entropy_totalsampen",
    "gi_fixed_r_sampen": "fwh_entropy_sampen_fixed_r",
}


@dataclass
class FitAttempt:
    fit: Any
    structure: str
    method: str
    messages: list
    boundary_warning: bool


def _attempt_fit(data: pd.DataFrame, method: str, use_vc: bool,
                 structure: str) -> FitAttempt:
    kwargs: dict = {"groups": data["patient_id"], "re_formula": "1"}
    if use_vc:
        kwargs["vc_formula"] = {"recording": "0 + C(recording_id)"}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = smf.mixedlm(FORMULA, data, **kwargs).fit(
            reml=True, method=method, maxiter=1000, disp=False)
        messages = [str(w.message) for w in caught]
    if not bool(getattr(fit, "converged", False)):
        raise RuntimeError(f"{method} optimizer did not report convergence")
    boundary = any(
        any(token in m.lower() for token in (
            "boundary", "random effects covariance is singular",
            "singular covariance"))
        for m in messages)
    return FitAttempt(fit, structure, method, messages, boundary)


def _try_tier(data: pd.DataFrame, use_vc: bool, structure: str,
              all_warnings: list):
    first_boundary = None
    for method in OPTIMIZER_METHODS:
        try:
            attempt = _attempt_fit(data, method, use_vc, structure)
        except Exception as exc:
            all_warnings.append(f"{structure}/{method}: {exc}")
            continue
        all_warnings.extend(f"{structure}/{method}: {m}" for m in attempt.messages)
        if not attempt.boundary_warning:
            return attempt
        if first_boundary is None:
            first_boundary = attempt
    return first_boundary


def fit_one(frame: pd.DataFrame, feature: str) -> list:
    data = frame[[feature, "window_class", "patient_id", "recording_id"]].copy()
    data = data.replace([np.inf, -np.inf], np.nan).dropna()
    data["window_class"] = pd.Categorical(
        data["window_class"],
        categories=["baseline", "fetal_movement", "contraction"])
    data = data.dropna(subset=["window_class"])

    base = {
        "feature": feature,
        "n": len(data),
        "n_women": int(data.patient_id.nunique()),
        "n_recordings": int(data.recording_id.nunique()),
        "n_baseline": int((data.window_class == "baseline").sum()),
        "n_contraction": int((data.window_class == "contraction").sum()),
        "n_fetal_movement": int((data.window_class == "fetal_movement").sum()),
        "model": "", "optimizer": "", "converged": False, "fallback": False,
        "boundary_fit": False, "boundary_reason": "", "warning_summary": "",
    }

    if base["n_baseline"] < MIN_BASELINE_N or base["n_contraction"] < MIN_CONTRACTION_N:
        base["warning_summary"] = "insufficient_data"
        return [{**base, "contrast": c, "estimate": float("nan"),
                 "std_error": float("nan"), "p_value": float("nan"),
                 "ci_low": float("nan"), "ci_high": float("nan")}
                for c in ("fetal_movement_vs_baseline", "contraction_vs_baseline")]

    raw_sd = float(data[feature].std())
    if not np.isfinite(raw_sd) or raw_sd == 0:
        base["warning_summary"] = "zero_variance"
        return [{**base, "contrast": c, "estimate": float("nan"),
                 "std_error": float("nan"), "p_value": float("nan"),
                 "ci_low": float("nan"), "ci_high": float("nan")}
                for c in ("fetal_movement_vs_baseline", "contraction_vs_baseline")]
    data["_z"] = (data[feature] - float(data[feature].mean())) / raw_sd

    all_warnings: list = []
    attempt = _try_tier(data, True, "woman_recording", all_warnings)
    fallback = False
    if attempt is None:
        attempt = _try_tier(data, False, "woman_only", all_warnings)
        fallback = True
    if attempt is None:
        base["warning_summary"] = "optimizer_failure"
        return [{**base, "contrast": c, "estimate": float("nan"),
                 "std_error": float("nan"), "p_value": float("nan"),
                 "ci_low": float("nan"), "ci_high": float("nan")}
                for c in ("fetal_movement_vs_baseline", "contraction_vs_baseline")]

    fit = attempt.fit
    base.update({"model": attempt.structure, "optimizer": attempt.method,
                 "converged": True, "fallback": fallback})

    woman_variance = float("nan")
    recording_variance = float("nan")
    try:
        woman_variance = float(fit.cov_re.iloc[0, 0])
    except Exception:
        pass
    if attempt.structure == "woman_recording":
        try:
            recording_variance = float(fit.vcomp[0])
        except Exception:
            pass

    reasons = []
    if attempt.boundary_warning:
        reasons.append("optimizer_boundary_warning")
    if np.isfinite(woman_variance) and woman_variance < VARIANCE_BOUNDARY_THRESHOLD:
        reasons.append("near_zero_woman_variance")
    if (attempt.structure == "woman_recording" and np.isfinite(recording_variance)
            and recording_variance < VARIANCE_BOUNDARY_THRESHOLD):
        reasons.append("near_zero_recording_variance")
    base["boundary_fit"] = bool(reasons)
    base["boundary_reason"] = ",".join(reasons)
    base["warning_summary"] = (
        f"{attempt.structure}_" + ("boundary" if reasons else "clean"))

    rows = []
    for contrast, key in (("fetal_movement_vs_baseline", FM_KEY),
                          ("contraction_vs_baseline", UC_KEY)):
        if key not in fit.params:
            continue
        ci = fit.conf_int().loc[key]
        rows.append({
            **base, "contrast": contrast, "term": key,
            "estimate": float(fit.params[key]),
            "std_error": float(fit.bse[key]),
            "statistic": float(fit.tvalues[key]),
            "p_value": float(fit.pvalues[key]),
            "ci_low": float(ci.iloc[0]), "ci_high": float(ci.iloc[1]),
        })
    return rows


def apply_fdr(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["q_value"] = np.nan
    for contrast in frame["contrast"].unique():
        mask = (frame["contrast"] == contrast) & frame["p_value"].notna()
        if mask.sum() == 0:
            continue
        frame.loc[mask, "q_value"] = multipletests(
            frame.loc[mask, "p_value"], method="fdr_bh")[1]
    frame["validated"] = frame["q_value"] < 0.05
    return frame


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    gi_path = _common.WORK_DIR / "gi_features_60s_EHG9.csv"
    if not gi_path.exists():
        raise RuntimeError(f"{gi_path} not found. Run extract_gi.py first.")
    frame = pd.read_csv(gi_path)

    available = [f for f in GI_FAMILY if f in frame.columns]
    missing = [f for f in GI_FAMILY if f not in frame.columns]
    if missing:
        print(f"warning: not in feature table, skipped: {missing}")

    rows = []
    for feature in available:
        if not args.quiet:
            print(f"fitting {feature}...", flush=True)
        rows.extend(fit_one(frame, feature))

    results = apply_fdr(pd.DataFrame(rows))

    # Attach the frozen counterpart estimates for the side-by-side.
    # results/primary/mixed_model_results_60s_fwh_EHG9.csv carries all three
    # contrasts; results/fm/fm_vs_in_34_features.csv carries only FM-vs-IN.
    frozen = pd.read_csv(
        _common.PAPER1_ROOT / "results" / "primary"
        / "mixed_model_results_60s_fwh_EHG9.csv")
    frozen = frozen.loc[frozen["channel"] == _common.CHANNEL]
    frozen_map = {}
    for _, row in frozen.iterrows():
        frozen_map[(row["feature"], row["contrast"])] = (
            row["estimate"], row["q_value"])

    def lookup(feature, contrast, index):
        column = FROZEN_COUNTERPART.get(feature)
        if column is None:
            return float("nan")
        value = frozen_map.get((column, contrast))
        return float(value[index]) if value is not None else float("nan")

    results["frozen_estimate"] = [
        lookup(f, c, 0) for f, c in zip(results["feature"], results["contrast"])]
    results["frozen_q_value"] = [
        lookup(f, c, 1) for f, c in zip(results["feature"], results["contrast"])]

    _common.write_csv(results, _common.RESULTS_DIR / "gi_model_results.csv",
                      label="GI model results")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
