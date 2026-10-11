"""Per-tail fitting of drying laws (library): the candidate models, least squares, BIC / AIC with an effective
sample size, and the population statistics used on real and simulated tails alike.

Candidates (t = minutes since t0_trace; every candidate has a free start W0 and a free offset c; A = W + c):
    F1  exponential         W = W0 exp(-t / tau)                                      k = 3 (W0, tau, c)
    F4  two layers          W = W0 ((1 - phi) exp(-t / tau_f) + phi exp(-t / tau_s)), tau_s > tau_f      k = 5
    F3  nonlinear           dW/dt = -k W^n, n in [1, 3]: W = (W0^(1-n) + (n - 1) k t)^(1 / (1 - n))      k = 4
    F5  linear + quadratic  dW/dt = -a W - b W^2: W = a W0 e^(-a t) / (a + b W0 (1 - e^(-a t)))           k = 4
Scores: the slow noise makes neighbouring minutes dependent, so the information in n minutes is less than n. Every
tail gets ONE effective sample size, n_eff = n (1 - rho) / (1 + rho), rho = lag-1 autocorrelation of the residuals of
the reference candidate (default F4, the most flexible). Then, for every candidate,
    BIC = n_eff ln(RSS / n) + k ln(n_eff),    AIC = n_eff ln(RSS / n) + 2 k.
The candidate with the lowest score wins the tail; summed over the tails of a dataset, the lowest sum wins the dataset.
"""
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

CANDIDATES = {"F1": ["W0", "tau", "c"], "F4": ["W0", "tau_f", "tau_s", "phi", "c"], "F3": ["W0", "k", "n", "c"],
              "F5": ["W0", "a", "b", "c"]}
N_BOUNDS = (1.0, 3.0)                                               # F3 exponent bounds (set by the caller)


def curve(cand, x, t):
    """Model A(t) for the internal parameter vector x of a candidate."""
    if cand == "F1":
        w0, tau, c = x
        return w0 * np.exp(-t / tau) + c
    if cand == "F4":
        w0, tau_f, dtau, phi, c = x                                 # tau_s = tau_f + dtau keeps tau_s > tau_f
        return w0 * ((1 - phi) * np.exp(-t / tau_f) + phi * np.exp(-t / (tau_f + dtau))) + c
    if cand == "F3":
        w0, k, n, c = x
        if abs(n - 1) < 1e-4:
            return w0 * np.exp(-k * t) + c
        return (w0 ** (1 - n) + (n - 1) * k * t) ** (1 / (1 - n)) + c
    if cand == "F5":
        w0, a, b, c = x
        e = np.exp(-a * t)
        return a * w0 * e / (a + b * w0 * (1 - e)) + c
    raise ValueError(cand)


def named(cand, x):
    """Internal vector -> named parameters (F4: tau_s = tau_f + dtau)."""
    if cand == "F4":
        return {"W0": x[0], "tau_f": x[1], "tau_s": x[1] + x[2], "phi": x[3], "c": x[4]}
    return dict(zip(CANDIDATES[cand], x))


def _bounds(cand):
    lo_n, hi_n = N_BOUNDS
    return {"F1": ([0, 0.5, -5], [30, 1e4, 5]),
            "F4": ([0, 0.5, 0.5, 0, -5], [30, 1e4, 1e5, 1, 5]),
            "F3": ([1e-3, 1e-7, lo_n, -5], [30, 10, hi_n, 5]),
            "F5": ([1e-3, 1e-6, 0, -5], [30, 2, 10, 5])}[cand]


def _starts(cand, w0, tau, c0):
    if cand == "F1":
        return [[w0, tau, c0]]
    if cand == "F4":
        return [[w0, tau / 2, tau * 4, 0.2, c0], [w0, tau / 3, tau * 20, 0.1, c0]]
    if cand == "F3":
        out = []
        for n in (1.3, 2.2):                                        # k from the half time ~ tau ln 2
            n = min(max(n, N_BOUNDS[0] + 1e-3), N_BOUNDS[1])
            out.append([w0, (2 ** (n - 1) - 1) / ((n - 1) * w0 ** (n - 1) * tau * np.log(2)), n, c0])
        return out
    return [[w0, 1 / tau, 0.01, c0], [w0, 0.2 / tau, 0.8 / (tau * w0), c0]]


def _clip(x, lo, hi):
    return np.minimum(np.maximum(x, np.asarray(lo) + 1e-9), np.asarray(hi) - 1e-9)


def lag1(r):
    """Lag-1 autocorrelation of residuals over consecutive valid minutes (r has NaN where A was NaN)."""
    a, b = r[:-1], r[1:]
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5:
        return 0.0
    a, b = a[m] - np.mean(r[np.isfinite(r)]), b[m] - np.mean(r[np.isfinite(r)])
    den = np.sqrt(np.sum(a * a) * np.sum(b * b))
    return float(np.sum(a * b) / den) if den > 0 else 0.0


def fit_tail(t, a, candidates=("F1", "F3", "F4", "F5"), n_eff_from="F4", n_eff_min=6):
    """Fit every candidate to one tail. t: minutes, a: 1 min A (NaN allowed). Returns a list of dicts."""
    ok = np.isfinite(a)
    tt, aa = t[ok].astype(float), a[ok].astype(float)
    n = len(aa)
    if n < 20:
        return []
    c0 = float(np.median(aa[-60:])) if n >= 120 else 0.0
    w0 = max(float(np.mean(aa[:5])) - c0, 0.1)
    sm = pd.Series(aa).rolling(5, min_periods=1).mean().to_numpy()
    hit = np.nonzero(sm - c0 < w0 / np.e)[0]
    tau = float(np.clip(tt[hit[0]] if len(hit) else tt[-1] / 3, 1, 1000))
    fits = {}
    for cand in candidates:
        lo, hi = _bounds(cand)
        best = None
        for x0 in _starts(cand, w0, tau, c0):
            try:
                res = least_squares(lambda x: curve(cand, x, tt) - aa, _clip(np.asarray(x0, float), lo, hi),
                                    bounds=(lo, hi), method="trf", x_scale="jac", max_nfev=400)
            except (ValueError, FloatingPointError):
                continue
            rss = float(np.sum(res.fun ** 2))
            if np.isfinite(rss) and (best is None or rss < best[1]):
                best = (res.x, rss, res.status > 0)
        if best is not None:
            fits[cand] = best
    if not fits:
        return []
    ref = n_eff_from if n_eff_from in fits else min(fits, key=lambda k: fits[k][1])
    r_full = np.full(len(a), np.nan)
    r_full[ok] = curve(ref, fits[ref][0], tt) - aa
    rho = float(np.clip(lag1(r_full), 0.0, 0.99))
    n_eff = float(np.clip(n * (1 - rho) / (1 + rho), n_eff_min, n))
    out = []
    for cand, (x, rss, conv) in fits.items():
        k = len(CANDIDATES[cand])
        ll = n_eff * np.log(max(rss, 1e-12) / n)
        out.append({"candidate": cand, "n": n, "rho": rho, "n_eff": n_eff, "k": k, "rss": rss,
                    "bic": ll + k * np.log(n_eff), "aic": ll + 2 * k, "converged": conv, **named(cand, x)})
    return out


def fit_many(job):
    """Worker entry (picklable): job = (n_bounds, candidates, n_eff_from, n_eff_min, [(tail_id, t, a), ...])."""
    global N_BOUNDS
    n_bounds, cands, ref, nmin, tails = job
    N_BOUNDS = tuple(n_bounds)
    rows = []
    for tid, t, a in tails:
        for r in fit_tail(np.asarray(t), np.asarray(a), cands, ref, nmin):
            rows.append({"tail_id": tid, **r})
    return rows


def winners(fits, score="bic"):
    """Per tail: the candidate with the lowest score. fits: DataFrame with tail_id, candidate, score columns."""
    i = fits.groupby("tail_id")[score].idxmin()
    return fits.loc[i, ["tail_id", "candidate"]].set_index("tail_id").candidate


def t_half(df, a_start):
    """Per tail: first minute the 5 min trailing mean of A < A_start / 2 (NaN = not reached)."""
    out = {}
    for tid, g in df.groupby("tail_id"):
        a = g.sort_values("minutes_since_t0_trace").A_db.to_numpy(float)
        sm = pd.Series(a).rolling(5, min_periods=1).mean().to_numpy()
        hit = np.nonzero(sm < a_start[tid] / 2)[0]
        out[tid] = int(hit[0]) if len(hit) else np.nan
    return pd.Series(out)


def dadt_curve(df, smooth, edges, min_count):
    """Median dA/dt (dB/min) per bin of smoothed A, pooled over tails. A is smoothed by a centred rolling mean of
    `smooth` minutes per tail; the derivative is the central difference of the smoothed A."""
    parts = []
    for _, g in df.groupby("tail_id"):
        a = g.sort_values("minutes_since_t0_trace").A_db
        s = a.rolling(smooth, center=True, min_periods=int(smooth * 2 / 3)).mean().to_numpy()
        d = np.full(len(s), np.nan)
        d[1:-1] = (s[2:] - s[:-2]) / 2
        parts.append(pd.DataFrame({"A": s, "dadt": d}))
    x = pd.concat(parts).dropna()
    b = pd.cut(x.A, edges)
    g = x.groupby(b, observed=False).dadt
    out = pd.DataFrame({"median_dadt": g.median(), "count": g.size()})
    out.loc[out["count"] < min_count, "median_dadt"] = np.nan
    out.index = [f"{iv.left:g}..{iv.right:g}" for iv in out.index]
    out["A_mid"] = [(lo + hi) / 2 for lo, hi in zip(edges[:-1], edges[1:])]
    return out
