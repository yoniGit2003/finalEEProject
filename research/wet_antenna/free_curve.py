"""Free drying curve for the pooled fit, and the equation of a curve (library, step 3b-3).

Jump 1, the free curve. The pooled fit of pooled_fit.py with a law that is NOT taken from a menu:
    dW/dt = -f_i W k(W),      k(W) > 0   ("drying rate per dB of water", 1 / min)
    log k = natural cubic spline in x = log W with K knots; the parameters are the K values of log k at the knots.
    Outside the knot range log k continues as a straight line in log W (that is what "natural" means), and k is
    kept inside [1e-4, 1] / min so that a wild extrapolation cannot break the time table.
    FLAT k = exponential;  RISING k = nonlinear in W (m3, m5);  a STEP in k = state switch (m6).
Everything else is the pooled likelihood of step 3b-1b (offset prior, ARMA(1,1) noise, per-tail start profiled,
speed factor f_i = exp(s g), g ~ N(0, 1), so the median f is 1, marginalised). The fit maximises
    log L(theta, s)  -  lambda * integral (d^2 log k / dx^2)^2 dx  -  sum_k (theta_k - prior_mean)^2 / (2 prior_sd^2)
(smoothness penalty; and a weak prior that only keeps knots without data from drifting away).

Jump 2, the equation of a curve. h(W) = W k(W) sampled on a grid is regressed on a small library
(W, W^2, W^3, max(W - A_c, 0)) with the thresholded backward elimination of sindy.py; no constant, h(0) = 0.
"""
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize

from research.wet_antenna import pooled_fit as wp
from research.wet_antenna import sindy as sd


class Curve:
    """The spline of log k on fixed knots (dB). theta = log k at the knots."""

    def __init__(self, knots_db, table_range_db=(0.01, 30.0), table_n=300, prior_mean=np.log(1 / 30), prior_sd=3.0):
        self.knots = np.asarray(knots_db, float)
        self.xk = np.log(self.knots)
        self.K = len(self.knots)
        self.xt = np.linspace(np.log(table_range_db[0]), np.log(table_range_db[1]), table_n)   # table grid (log W)
        self.Bt = self.basis(self.xt)
        fine = np.linspace(self.xk[0], self.xk[-1], 400)
        d2 = CubicSpline(self.xk, np.eye(self.K), bc_type="natural")(fine, 2)
        self.omega = d2.T @ d2 * (fine[1] - fine[0])               # theta' omega theta = integral (log k)''^2 dx
        self.prior_mean, self.prior_sd = prior_mean, prior_sd
        self.fit_noise = False                                     # step 4b: set to True to fit the noise-variance
                                                                   #   multiplier too; x is then [theta, log s, log m]

    def basis(self, x):
        """Matrix B with log k(x) = B theta (x = log W; straight-line continuation outside the knots)."""
        x = np.atleast_1d(np.asarray(x, float))
        cs = CubicSpline(self.xk, np.eye(self.K), bc_type="natural")
        B = cs(np.clip(x, self.xk[0], self.xk[-1]))
        for end, out in ((self.xk[0], x < self.xk[0]), (self.xk[-1], x > self.xk[-1])):
            B[out] += (x[out] - end)[:, None] * cs(end, 1)[None, :]
        return B

    def k(self, theta, w_db):
        return np.exp(self.basis(np.log(np.asarray(w_db, float))) @ np.asarray(theta, float))

    def table(self, theta):
        """The table that the law FC of pooled_fit.py reads: [n, x_0, dx, T, k]."""
        return self._pack_table(np.exp(np.clip(self.Bt @ np.asarray(theta, float), np.log(1e-4), 0.0)))

    def _pack_table(self, kt):
        """k on the table grid -> [n, x_0, dx, T, k] with T(x) = integral of dx / k (trapezoid rule)."""
        dx = self.xt[1] - self.xt[0]
        T = np.concatenate([[0.0], np.cumsum(0.5 * dx * (1 / kt[1:] + 1 / kt[:-1]))])
        return np.concatenate([[len(self.xt), self.xt[0], dx], T, kt])

    def penalty(self, theta, lam):
        th = np.asarray(theta, float)
        return lam * float(th @ self.omega @ th) + float(np.sum((th - self.prior_mean) ** 2)) / (2 * self.prior_sd ** 2)

    def loglik(self, x, D):
        """Pooled marginal log-likelihood (no penalty). x = [theta_1 .. theta_K, log s] (+ [log m] with fit_noise)."""
        if self.fit_noise:
            ll = wp.pooled_loglik("FC", self.table(x[:-2]), float(np.exp(x[-2])), D, noise_mult=float(np.exp(x[-1])))
        else:
            ll = wp.pooled_loglik("FC", self.table(x[:-1]), float(np.exp(x[-1])), D)
        return ll if np.isfinite(ll) else -1e18

    def cost(self, x, D, lam):
        return -self.loglik(x, D) + self.penalty(x[:-2] if self.fit_noise else x[:-1], lam)

    def fit(self, D, x0, lam, maxfev, step=0.4, fatol=0.1):
        """Nelder-Mead from x0 (as pooled_fit.fit: explicit simplex, `step` in every coordinate)."""
        x0 = np.asarray(x0, float)
        simplex = np.vstack([x0] + [x0 + step * np.eye(len(x0))[i] for i in range(len(x0))])
        res = minimize(self.cost, x0, args=(D, lam), method="Nelder-Mead",
                       options={"maxfev": maxfev, "initial_simplex": simplex, "xatol": 1e-3, "fatol": fatol})
        return {"x": [float(v) for v in res.x], "cost": float(res.fun), "ll": float(self.loglik(res.x, D)),
                "nfev": int(res.nfev), "success": bool(res.success)}

    def covariance(self, x, D, lam, h=0.15):
        """Covariance of x from the curvature of the penalised cost (central finite differences, step h).
        Returns None if the curvature matrix is not positive definite (the optimum is then not well determined)."""
        x = np.asarray(x, float)
        p = len(x)
        f0 = self.cost(x, D, lam)
        e = np.eye(p) * h
        fp = np.array([self.cost(x + e[i], D, lam) for i in range(p)])
        fm = np.array([self.cost(x - e[i], D, lam) for i in range(p)])
        H = np.zeros((p, p))
        for i in range(p):
            H[i, i] = (fp[i] - 2 * f0 + fm[i]) / h ** 2
            for j in range(i):
                fpp = self.cost(x + e[i] + e[j], D, lam)
                fmm = self.cost(x - e[i] - e[j], D, lam)
                H[i, j] = H[j, i] = (fpp - fp[i] - fp[j] + 2 * f0 - fm[i] - fm[j] + fmm) / (2 * h ** 2)
        if not np.all(np.isfinite(H)) or np.linalg.eigvalsh(H).min() <= 1e-8:
            return None
        return np.linalg.inv(H)


class StepCurve(Curve):
    """The free curve with an optional smooth step (step 3b-3c). p = [theta_1 .. theta_K, u_c, u_d]:
        log k(x) = g(x) spline(x) + (1 - g(x)) (spline(x_c) - d),   g = logistic((x - x_c) / width),  x = log W
    so k = k_low = spline(W_c) exp(-d) (a constant) below W_c and the spline above it. x_c = log W_c is kept inside
    `wc_range_db` by a logistic transform of u_c; the step height d = exp(u_d) >= 0 (log units). d -> 0 with W_c
    below the data gives back the plain curve only approximately (k is constant below W_c instead of a straight
    line in log-log), so the comparison with the plain curve is by BIC, not by a nested test."""

    def __init__(self, knots_db, wc_range_db, width, **kw):
        super().__init__(knots_db, **kw)
        self.xc_lo, self.xc_hi = np.log(wc_range_db[0]), np.log(wc_range_db[1])
        self.width = width

    def u_of(self, wc_db, height):
        q = (np.log(wc_db) - self.xc_lo) / (self.xc_hi - self.xc_lo)
        return [float(np.log(q / (1 - q))), float(np.log(height))]

    def step(self, p):
        """(W_c in dB, step height d in log units) of the parameter vector p."""
        xc = self.xc_lo + (self.xc_hi - self.xc_lo) / (1 + np.exp(-np.clip(p[self.K], -30, 30)))
        return float(np.exp(xc)), float(np.exp(np.clip(p[self.K + 1], -8, 3)))

    def _logk(self, p, x, B):
        th = np.asarray(p[:self.K], float)
        wc, d = self.step(p)
        g = 1 / (1 + np.exp(-np.clip((x - np.log(wc)) / self.width, -50, 50)))
        return g * (B @ th) + (1 - g) * (float(self.basis(np.log(wc))[0] @ th) - d)

    def k(self, p, w_db):
        x = np.log(np.asarray(w_db, float))
        return np.exp(self._logk(p, x, self.basis(x)))

    def table(self, p):
        return self._pack_table(np.exp(np.clip(self._logk(p, self.xt, self.Bt), np.log(1e-4), 0.0)))

    def penalty(self, p, lam):
        return super().penalty(p[:self.K], lam)


class FormulaLaw(Curve):
    """A named formula h(W) = sum_j c_j phi_j(W) as a law of the pooled model (jump 3): k(W) = h(W) / W through the
    same table as the free curve. Parameters u_j = c_j / c0_j (c0 = the coefficients jump 2 found, so the start is
    u = 1). Terms: 'W', 'W^2', 'W^3', 'sw<x>' = max(W - x dB, 0). k is kept inside [1e-4, 1] / min."""

    def __init__(self, terms, c0, table_range_db=(0.01, 30.0), table_n=300):
        self.terms, self.c0, self.K = list(terms), np.asarray(c0, float), len(terms)
        self.xt = np.linspace(np.log(table_range_db[0]), np.log(table_range_db[1]), table_n)
        wt = np.exp(self.xt)
        self.Gt = term_matrix(wt, self.terms) / wt[:, None]

    def coef(self, u):
        return np.asarray(u, float) * self.c0

    def k(self, u, w_db):
        w = np.asarray(w_db, float)
        return np.clip(term_matrix(w, self.terms) @ self.coef(u) / w, 1e-4, 1.0)

    def table(self, u):
        return self._pack_table(np.clip(self.Gt @ self.coef(u), 1e-4, 1.0))

    def penalty(self, u, lam):
        return 0.0


def slope_weights(x):
    """Weights a with (least-squares slope of y on x) = a @ y."""
    xc = np.asarray(x, float) - np.mean(x)
    return xc / np.sum(xc ** 2)


# ------------------------------------------------------------------------------------------------- jump 2
def curve_terms(switch_grid):
    return ["W", "W^2", "W^3"] + [f"sw{g:g}" for g in switch_grid]


def curve_library(w, switch_grid):
    w = np.asarray(w, float)
    return np.column_stack([w, w ** 2, w ** 3] + [np.maximum(w - g, 0.0) for g in switch_grid])


def name_curves(w, H, switch_grid, lambdas):
    """Name the equation of every curve. w: grid (dB); H: curves x grid, values of h(W) = W k(W).
    Weights = 1 / variance of h over the curves (floor: 1 % of the mean curve). The threshold is chosen by
    leave-one-curve-out cross-validation (terms and coefficients from the mean of the other curves, error on the
    curve left out; one-standard-error rule: the largest threshold within one SE of the best).
    Returns (lambda, selected: curves x terms (bool), coefficients: curves x terms)."""
    H = np.asarray(H, float)
    G = curve_library(w, switch_grid)
    mean = H.mean(axis=0)
    wt = 1.0 / np.sqrt(H.var(axis=0, ddof=1) + (0.01 * mean) ** 2)
    Gw = G * wt[:, None]
    err = np.zeros((len(lambdas), len(H)))
    for i, lam in enumerate(lambdas):
        for r in range(len(H)):
            others = np.delete(H, r, axis=0).mean(axis=0) * wt
            xi = sd.select(Gw, others, lam)
            res = H[r] * wt - Gw @ xi
            err[i, r] = float(res @ res) / float((H[r] * wt) @ (H[r] * wt))
    m, se = err.mean(axis=1), err.std(axis=1, ddof=1) / np.sqrt(len(H))
    best = int(np.argmin(m))
    lam = float(np.asarray(lambdas)[np.flatnonzero(m <= m[best] + se[best]).max()])
    coef = np.array([sd.select(Gw, H[r] * wt, lam) for r in range(len(H))])
    return lam, coef != 0, coef



# ------------------------------------------------------------------------------------------------- jump 2, one curve
def term_matrix(w, terms):
    """Library columns for the named terms: 'W', 'W^2', 'W^3', 'sw<x>' = max(W - x dB, 0)."""
    w = np.asarray(w, float)
    cols = {"W": lambda: w, "W^2": lambda: w ** 2, "W^3": lambda: w ** 3}
    return np.column_stack([cols[t]() if t in cols else np.maximum(w - float(t[2:]), 0.0) for t in terms])


def best_subset(w, h, terms, max_terms=3):
    """Every combination of up to `max_terms` library terms, least squares on RELATIVE residuals (each point divided
    by h, so a 5 % miss counts the same at every W). Returns one row per combination, best BIC first:
    {terms, coef, rss, bic, size}; BIC = n log(RSS / n) + (number of terms) log n."""
    from itertools import combinations
    w, h = np.asarray(w, float), np.asarray(h, float)
    G = term_matrix(w, terms) / h[:, None]
    y = np.ones(len(w))
    out = []
    for size in range(1, max_terms + 1):
        for ix in combinations(range(len(terms)), size):
            c, *_ = np.linalg.lstsq(G[:, ix], y, rcond=None)
            r = y - G[:, ix] @ c
            rss = max(float(r @ r), 1e-300)
            out.append({"terms": [terms[i] for i in ix], "coef": [float(v) for v in c], "rss": rss, "size": size,
                        "bic": len(w) * np.log(rss / len(w)) + size * np.log(len(w))})
    return sorted(out, key=lambda r: r["bic"])


def bagged_select(w, h, terms, lam, n_boot, rng, method="waa_sindy"):
    """Ensemble SINDy on one curve: the points are resampled with replacement `n_boot` times, a sparse regression
    is run on every resample (relative residuals, unit-norm columns), and a term is reported by its inclusion
    probability and the median of its coefficients over the resamples where it was included.
    method 'waa_sindy': thresholded backward elimination of waa_sindy.select (relative threshold lam);
    method 'pysindy': STLSQ of the pysindy package (threshold lam * ||y|| on the unit-norm columns).
    Returns (inclusion probability, median coefficient) per term."""
    w, h = np.asarray(w, float), np.asarray(h, float)
    G = term_matrix(w, terms) / h[:, None]
    norm = np.linalg.norm(G, axis=0)
    Gn, y = G / norm, np.ones(len(w))
    coefs = np.zeros((n_boot, len(terms)))
    if method == "pysindy":
        from pysindy.optimizers import STLSQ
    for b in range(n_boot):
        ix = rng.integers(0, len(w), len(w))
        if method == "pysindy":
            opt = STLSQ(threshold=lam * np.linalg.norm(y[ix]), alpha=1e-8)
            opt.fit(Gn[ix], y[ix].reshape(-1, 1))
            coefs[b] = np.asarray(opt.coef_).ravel() / norm
        else:
            coefs[b] = sd.select(Gn[ix], y[ix], lam) / norm
    inc = (coefs != 0).mean(axis=0)
    med = np.array([np.median(coefs[coefs[:, j] != 0, j]) if inc[j] > 0 else 0.0 for j in range(len(terms))])
    return inc, med


def true_k(truth, par, w):
    """k(W) = -(dW/dt) / W of the base law of a simulated truth (speed factor 1); None for m4 (two states: not a
    function of W). par: {truth: its parameters}."""
    w = np.asarray(w, float)
    p = par[truth]
    if truth == "m1":
        return np.full_like(w, 1 / p["tau"])
    if truth == "m3":
        return p["c"] * w ** (p["n"] - 1)
    if truth == "m5":
        return p["a"] + p["b"] * w
    if truth == "m6":
        return np.where(w > p["A_c"], p["k_evap"] + p["k_run"] * (1 - p["A_c"] / w), p["k_evap"])
    return None


def start_theta(spec, knots):
    """Start values of log k at the knots from a start description of the config: a straight line in log-log
    ({k1, slope}) or a step ({at_db, low, high})."""
    if "slope" in spec:
        return np.log(spec["k1"]) + spec["slope"] * np.log(knots)
    return np.log(np.where(knots <= spec["at_db"], spec["low"], spec["high"]))
