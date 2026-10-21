#!/usr/bin/env python3
"""Gain-invariant sample-entropy profiling -- Paper 2 core numerics.

Design notes
------------
1.  The frozen Paper 1 algorithm is RE-IMPLEMENTED here rather than imported.
    That is deliberate. It lets `parity_check` assert, numerically, that this
    file reproduces the frozen Stage 04 `fwh_entropy_*sampen` columns before
    anything new is interpreted. It also keeps Paper 1 untouched.

2.  Two gain-invariant variants are provided, and they answer different
    questions.

    (A) `gi_summaries_rounded` -- the FROZEN algorithm applied to
        z = x / sd(x).  Structurally identical to Paper 1: same Chebyshev
        distances, same 3-decimal rounding, same retention rule, same six
        summaries.  The ONLY change is that distances are expressed in units
        of the window's own standard deviation.  These are the features the
        go/no-go gate is evaluated on, because they are the direct,
        like-for-like analogues of the frozen ones.

        Consequence to state in the paper: normalising also REFINES the
        tolerance resolution.  EHG FWH windows have sd of order 1e-2 signal
        units, so a 0.001 signal-unit grid gives only ~1e2 bins, whereas a
        0.001-sd grid gives ~1e3-1e4.  The GI analogues are therefore not
        merely "the frozen features made invariant"; they are computed on a
        finer profile.  This is a feature-definition change, not a bug fix.

    (B) `gi_summaries_grid` -- the PROPOSED representation.  The profile is
        evaluated on a FIXED dimensionless tolerance grid shared by every
        window, with no rounding at all.  This is what makes profiles
        comparable ACROSS windows: variant (A), like the frozen feature, still
        runs out to each window's own maximum distance, so its support varies
        between windows even though it is invariant to gain within a window.
        GI-AUC, GI-slope and GI-centroid are defined on this fixed grid.

3.  Exact gain invariance.  For c > 0, sd(cx) = c sd(x) and every Chebyshev
    distance scales by c, so z is unchanged and BOTH variants are invariant up
    to floating point.  `invariance_check.py` asserts this empirically; it is
    a correctness test on this file, not a scientific result.

Author: drafted for the ICASSP go/no-go gate. Not part of the frozen Paper 1
pipeline.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import kurtosis, skew

# --- Frozen Paper 1 constants (src/config.py). Do not change. --------------
EMBEDDING_DIMENSION = 2
R_FRACTION = 0.2
PROFILE_DISTANCE_DECIMALS = 3
PROFILE_EPSILON = 1e-12

# --- Paper 2 additions -----------------------------------------------------
# Rounding resolution for the GI rounded variant, in units of sd. Kept at 3 to
# mirror the frozen algorithm exactly; see design note 2(A).
GI_DISTANCE_DECIMALS = 3

SUMMARY_NAMES = ("total", "average", "median", "sd", "kurtosis", "skewness")

# numpy renamed trapz -> trapezoid in 2.0; support both so the script runs
# under whatever numpy the Paper 1 conda environment pins.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz


# ===========================================================================
# Shared primitives
# ===========================================================================

def _validate(signal) -> np.ndarray:
    data = np.asarray(signal, dtype=np.float64).reshape(-1)
    if data.size == 0 or not np.all(np.isfinite(data)):
        return np.empty(0, dtype=np.float64)
    return data


def _aligned_templates(data: np.ndarray, dimension: int, count: int) -> np.ndarray:
    """Frozen Paper 1 template construction (feature_extraction._aligned_templates)."""
    return np.column_stack([data[offset:count + offset] for offset in range(dimension)])


def _rounded_shell_counts(templates: np.ndarray, maximum_bin: int, *, decimals: int) -> np.ndarray:
    """Frozen Paper 1 shell counting, including the half-even rounding boundaries.

    Reproduces feature_extraction._rounded_distance_shell_counts exactly.
    """
    scale = float(10 ** decimals)
    bins = np.arange(maximum_bin + 1, dtype=np.float64)
    upper_boundaries = (bins + 0.5) / scale
    odd_bins = (np.arange(maximum_bin + 1) % 2) == 1
    upper_boundaries[odd_bins] = np.nextafter(upper_boundaries[odd_bins], -np.inf)
    tree = cKDTree(templates)
    counts = np.asarray(
        tree.count_neighbors(tree, upper_boundaries, p=np.inf, cumulative=False),
        dtype=np.int64,
    )
    counts[0] -= templates.shape[0]
    if np.any(counts < 0):
        raise RuntimeError("invalid sample-entropy profile pair counts")
    return counts


def _nan_summaries(prefix: str) -> dict:
    out = {f"{prefix}_{name}": float("nan") for name in SUMMARY_NAMES}
    out[f"{prefix}_n_retained"] = float("nan")
    out[f"{prefix}_max_rounded_distance"] = float("nan")
    return out


def _summarise(values: np.ndarray, prefix: str, max_rounded: float) -> dict:
    """Frozen summary set: sum, mean, median, population sd, Fisher excess
    kurtosis (bias=True), skewness (bias=True) -- scipy defaults, as Paper 1."""
    return {
        f"{prefix}_total": float(np.sum(values)),
        f"{prefix}_average": float(np.mean(values)),
        f"{prefix}_median": float(np.median(values)),
        f"{prefix}_sd": float(np.std(values)),
        f"{prefix}_kurtosis": float(kurtosis(values)),
        f"{prefix}_skewness": float(skew(values)),
        f"{prefix}_n_retained": float(values.size),
        f"{prefix}_max_rounded_distance": float(max_rounded),
    }


def _rounded_profile_values(data: np.ndarray, *, order: int, decimals: int,
                            epsilon: float) -> tuple[np.ndarray, float]:
    """Retained profile values under the frozen rounded-grid construction."""
    template_count = data.size - order
    templates_m = _aligned_templates(data, order, template_count)
    templates_m1 = _aligned_templates(data, order + 1, template_count)

    maximum_distance = max(
        float(np.max(np.ptp(templates_m, axis=0))),
        float(np.max(np.ptp(templates_m1, axis=0))),
    )
    maximum_rounded = float(np.round(maximum_distance, decimals))
    maximum_bin = int(np.rint(maximum_rounded * (10 ** decimals)))

    shell_m = _rounded_shell_counts(templates_m, maximum_bin, decimals=decimals)
    shell_m1 = _rounded_shell_counts(templates_m1, maximum_bin, decimals=decimals)

    observed = (shell_m > 0) | (shell_m1 > 0)
    denominator = float(template_count * (template_count - 1))
    b = np.cumsum(shell_m, dtype=np.float64) / denominator
    a = np.cumsum(shell_m1, dtype=np.float64) / denominator

    profile = np.log((b + epsilon) / (a + epsilon))
    valid = observed & np.isfinite(profile) & (b > epsilon) & (a > epsilon)
    return profile[valid], maximum_rounded


# ===========================================================================
# (0) Frozen reimplementation -- parity target
# ===========================================================================

def frozen_summaries(signal, *, order: int = EMBEDDING_DIMENSION,
                     decimals: int = PROFILE_DISTANCE_DECIMALS,
                     epsilon: float = PROFILE_EPSILON) -> dict:
    """Reimplementation of Paper 1 `sample_entropy_profile_summaries`.

    Must reproduce the cached `fwh_entropy_{total,average,median,sd,kurtosis,
    skewness}sampen` columns. Used only by the parity check.
    """
    data = _validate(signal)
    if data.size == 0 or order < 1 or data.size < order + 2:
        return _nan_summaries("frozen")
    values, max_rounded = _rounded_profile_values(
        data, order=order, decimals=decimals, epsilon=epsilon)
    if values.size == 0:
        return _nan_summaries("frozen")
    return _summarise(values, "frozen", max_rounded)


# ===========================================================================
# (A) Gain-invariant analogues -- frozen algorithm on sd-normalised data
# ===========================================================================

def gi_summaries_rounded(signal, *, order: int = EMBEDDING_DIMENSION,
                         decimals: int = GI_DISTANCE_DECIMALS,
                         epsilon: float = PROFILE_EPSILON) -> dict:
    """Frozen construction applied to z = x / sd(x).

    Exactly one line differs from `frozen_summaries`: the normalisation.
    These are the gate features.
    """
    data = _validate(signal)
    if data.size == 0 or order < 1 or data.size < order + 2:
        return _nan_summaries("gi")
    sd = float(np.std(data, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        return _nan_summaries("gi")

    z = data / sd                                   # <-- the only change

    values, max_rounded = _rounded_profile_values(
        z, order=order, decimals=decimals, epsilon=epsilon)
    if values.size == 0:
        return _nan_summaries("gi")
    return _summarise(values, "gi", max_rounded)


# ===========================================================================
# (B) Proposed representation -- fixed dimensionless grid, no rounding
# ===========================================================================

def gi_profile_on_grid(signal, taus: np.ndarray, *,
                       order: int = EMBEDDING_DIMENSION,
                       epsilon: float = PROFILE_EPSILON):
    """Entropy profile of z = x / sd(x) evaluated at the tolerances `taus`.

    No distance rounding is involved: pair counts are taken directly at each
    tolerance, so the profile is exact at every grid point and identical for
    every window's grid. Returns (profile_values, valid_mask).
    """
    data = _validate(signal)
    taus = np.asarray(taus, dtype=np.float64).reshape(-1)
    if data.size == 0 or order < 1 or data.size < order + 2 or taus.size == 0:
        return np.full(taus.size, np.nan), np.zeros(taus.size, dtype=bool)
    sd = float(np.std(data, ddof=0))
    if not np.isfinite(sd) or sd <= 0:
        return np.full(taus.size, np.nan), np.zeros(taus.size, dtype=bool)

    z = data / sd
    template_count = z.size - order
    tree_m = cKDTree(_aligned_templates(z, order, template_count))
    tree_m1 = cKDTree(_aligned_templates(z, order + 1, template_count))

    # cumulative=False does one traversal over sorted radii; cumulative=True
    # is markedly slower here. Shell counts are accumulated manually, and the
    # ordered self-pairs (distance 0) are removed from the first shell, which
    # is what the frozen implementation does.
    def cumulative_counts(tree: cKDTree) -> np.ndarray:
        shells = np.asarray(
            tree.count_neighbors(tree, taus, p=np.inf, cumulative=False),
            dtype=np.float64)
        shells[0] -= template_count
        return np.cumsum(shells)

    counts_m = cumulative_counts(tree_m)
    counts_m1 = cumulative_counts(tree_m1)

    denominator = float(template_count * (template_count - 1))
    b = counts_m / denominator
    a = counts_m1 / denominator

    profile = np.log((b + epsilon) / (a + epsilon))
    valid = np.isfinite(profile) & (b > epsilon) & (a > epsilon)
    profile = np.where(valid, profile, np.nan)
    return profile, valid


def gi_summaries_grid(signal, taus: np.ndarray, *,
                      order: int = EMBEDDING_DIMENSION,
                      epsilon: float = PROFILE_EPSILON,
                      min_valid_fraction: float = 0.90) -> dict:
    """GI-AUC and shape descriptors on the fixed dimensionless grid.

    GI-AUC is the mean profile height over the grid (the normalised area,
    i.e. the integral divided by the grid width). Slope is the OLS slope of
    the profile against tau. Centroid is the tau-weighted centre of mass of
    the profile; profile values are non-negative by construction (b >= a), so
    this is well defined.
    """
    taus = np.asarray(taus, dtype=np.float64).reshape(-1)
    nan_out = {
        "gigrid_auc": float("nan"), "gigrid_slope": float("nan"),
        "gigrid_centroid": float("nan"), "gigrid_sd": float("nan"),
        "gigrid_kurtosis": float("nan"), "gigrid_skewness": float("nan"),
        "gigrid_valid_fraction": float("nan"),
    }
    profile, valid = gi_profile_on_grid(signal, taus, order=order, epsilon=epsilon)
    valid_fraction = float(valid.mean()) if valid.size else 0.0
    if valid_fraction < min_valid_fraction:
        nan_out["gigrid_valid_fraction"] = valid_fraction
        return nan_out

    t = taus[valid]
    e = profile[valid]
    if t.size < 3:
        nan_out["gigrid_valid_fraction"] = valid_fraction
        return nan_out

    width = float(t[-1] - t[0])
    auc = float(_trapezoid(e, t) / width) if width > 0 else float("nan")
    slope = float(np.polyfit(t, e, 1)[0])
    mass = float(np.sum(e))
    centroid = float(np.sum(t * e) / mass) if mass > 0 else float("nan")

    return {
        "gigrid_auc": auc,
        "gigrid_slope": slope,
        "gigrid_centroid": centroid,
        "gigrid_sd": float(np.std(e)),
        "gigrid_kurtosis": float(kurtosis(e)),
        "gigrid_skewness": float(skew(e)),
        "gigrid_valid_fraction": valid_fraction,
    }


# ===========================================================================
# Comparators
# ===========================================================================

def fixed_r_sampen(signal, *, m: int = EMBEDDING_DIMENSION,
                   r_fraction: float = R_FRACTION) -> float:
    """Richman-Moorman SampEn, r = 0.2 sd. Frozen Paper 1 comparator.

    Already exactly gain-invariant; carried through as the reference point
    against which the profile representations must justify themselves.
    """
    data = _validate(signal)
    if data.size < m + 3:
        return float("nan")
    sd = float(np.std(data, ddof=0))
    if sd == 0:
        return float("nan")
    tolerance = r_fraction * sd

    templates_m = np.lib.stride_tricks.sliding_window_view(data, m)
    templates_m1 = np.lib.stride_tricks.sliding_window_view(data, m + 1)
    aligned_count = templates_m1.shape[0]
    templates_m = np.asarray(templates_m[:aligned_count], dtype=np.float64)
    templates_m1 = np.asarray(templates_m1, dtype=np.float64)

    def ordered_nonself(templates: np.ndarray) -> float:
        tree = cKDTree(templates)
        return float(tree.count_neighbors(tree, tolerance, p=np.inf) - aligned_count)

    denominator = float(aligned_count * (aligned_count - 1))
    b = ordered_nonself(templates_m) / denominator
    a = ordered_nonself(templates_m1) / denominator
    if b <= 0 or a <= 0:
        return float("nan")
    return float(-np.log(a / b))


def rms(signal) -> float:
    data = _validate(signal)
    if data.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(np.square(data))))


# ===========================================================================
# Convenience: everything for one window
# ===========================================================================

def all_features(signal, taus: np.ndarray, *, include_frozen: bool = False) -> dict:
    """Compute the Paper 2 feature block for a single window."""
    out: dict = {}
    if include_frozen:
        out.update(frozen_summaries(signal))
    out.update(gi_summaries_rounded(signal))
    out.update(gi_summaries_grid(signal, taus))
    out["gi_fixed_r_sampen"] = fixed_r_sampen(signal)
    out["rms"] = rms(signal)
    return out


if __name__ == "__main__":
    # Smoke test on synthetic data: assert exact gain invariance.
    rng = np.random.default_rng(0)
    x = np.cumsum(rng.standard_normal(1200)) * 1e-2
    taus = np.arange(0.05, 2.0 + 1e-12, 0.01)
    for c in (0.25, 0.5, 2.0, 4.0):
        base = all_features(x, taus)
        scaled = all_features(c * x, taus)
        worst = max(
            abs(scaled[k] - base[k]) / max(abs(base[k]), 1e-12)
            for k in base
            if k not in ("rms",) and np.isfinite(base[k]) and np.isfinite(scaled[k])
        )
        print(f"scale {c:>5}: worst relative deviation = {worst:.3e}")
