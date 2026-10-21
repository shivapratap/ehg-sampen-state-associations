"""Small REML mixed-model fitter used for the conditional (log-RMS) analyses.

Model: y = X b + woman intercept + recording-within-woman intercept + noise,
fitted by REML with a block-diagonal covariance (one block per woman).
Cross-checked against statsmodels MixedLM (scripts/verify_statsmodels.py).
Variance components are log-parameterised and bounded below by exp(-30).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.linalg import LinAlgError, cho_factor, cho_solve
from scipy.optimize import minimize
from scipy.stats import norm

LOG_VAR_BOUNDS = (-30.0, 6.0)
BAD = 1e12   # objective returned when a trial point gives a numerically non-PD covariance
STARTS = [(-2, -2, -1), (-4, -3, -1), (0, 0, 0), (-6, -6, -1)]


def reml_fit(y, X, woman, rec):
    """Return (beta, se, (var_woman, var_recording, var_resid), scipy result)."""
    blocks = []
    for w in np.unique(woman):
        idx = np.where(woman == w)[0]
        r = rec[idx]
        blocks.append((idx, np.ones((len(idx), len(idx))), (r[:, None] == r[None, :]).astype(float)))

    def parts(theta):
        su2, sv2, se2 = np.exp(theta)
        A = np.zeros((X.shape[1],) * 2)
        b = np.zeros(X.shape[1])
        logdet, cache = 0.0, []
        for idx, J, S in blocks:
            V = su2 * J + sv2 * S + se2 * np.eye(len(idx))
            try:
                c = cho_factor(V, lower=True)
            except LinAlgError:
                # rounding error on some BLAS builds; retry with a negligible ridge, else reject the trial point
                try:
                    c = cho_factor(V + 1e-9 * np.mean(np.diag(V)) * np.eye(len(idx)), lower=True)
                except LinAlgError:
                    return None
            logdet += 2.0 * np.sum(np.log(np.diag(c[0])))
            A += X[idx].T @ cho_solve(c, X[idx])
            b += X[idx].T @ cho_solve(c, y[idx])
            cache.append((idx, c))
        beta = np.linalg.solve(A, b)
        quad = sum((y[i] - X[i] @ beta) @ cho_solve(c, y[i] - X[i] @ beta) for i, c in cache)
        return A, beta, logdet, quad

    def nll(theta):
        out = parts(theta)
        if out is None:
            return BAD
        A, _, logdet, quad = out
        return 0.5 * (logdet + np.linalg.slogdet(A)[1] + quad)

    best = None
    for s in STARTS:
        res = minimize(nll, np.array(s, float), method="L-BFGS-B", bounds=[LOG_VAR_BOUNDS] * 3)
        if best is None or res.fun < best.fun:
            best = res
    out = parts(best.x)
    if out is None:
        raise RuntimeError("REML optimum is numerically non-positive-definite")
    A, beta, _, _ = out
    return beta, np.sqrt(np.diag(np.linalg.inv(A))), np.exp(best.x), best


def zscore(s):
    s = np.asarray(s, float)
    return (s - s.mean()) / s.std(ddof=1)


def load_windows(path):
    """Window-level table restricted to the three states, with uc/fm indicators."""
    df = pd.read_csv(path)
    df = df[df.window_class.isin(["baseline", "contraction", "fetal_movement"])].copy()
    df["uc"] = (df.window_class == "contraction").astype(float)
    df["fm"] = (df.window_class == "fetal_movement").astype(float)
    return df


def fit_table(df, outcome, adjust=False):
    """UC/FM coefficients (vs IN) for one outcome, optionally adjusted for standardised log RMS."""
    cols = [np.ones(len(df)), df.uc.values, df.fm.values]
    names = ["intercept_IN", "UC_minus_IN", "FM_minus_IN"]
    if adjust:
        cols.append(zscore(np.log(df.rms.values)))
        names.append("log_RMS_std")
    beta, se, var, res = reml_fit(zscore(df[outcome]), np.column_stack(cols),
                                  pd.factorize(df.patient_id)[0], pd.factorize(df.recording_id)[0])
    out = pd.DataFrame({"term": names, "estimate": beta, "se": se})
    out["ci_low"] = out.estimate - norm.ppf(0.975) * out.se
    out["ci_high"] = out.estimate + norm.ppf(0.975) * out.se
    out["p"] = 2 * norm.sf(np.abs(out.estimate / out.se))
    out.insert(0, "outcome", outcome)
    out.insert(1, "adjusted_log_RMS", adjust)
    out["var_woman"], out["var_recording"], out["var_resid"] = var
    out["converged"] = bool(res.success)
    return out
