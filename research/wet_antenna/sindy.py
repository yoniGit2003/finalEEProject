"""Minimal weak-form SINDy for ONE state variable, pooled over many tails (library, step 3b-2).

pysindy is not installed in this environment, so the method is written out here (numpy only). It looks for a sparse
drying law
    dA/dt = sum_j xi_j theta_j(A),      theta = 1, A, A^2, A^3, max(A - A_c, 0) for a grid of A_c ("switch" terms)

Weak form (no numerical derivative of the noisy data): multiply the equation by a smooth test function phi(t) that is
zero at both ends of a short window and integrate by parts:
    - integral phi'(t) A(t) dt  =  sum_j xi_j  integral phi(t) theta_j(A(t)) dt
     \\________ b ___________/                 \\__________ G_j ____________/
One window gives one row (b, G); every tail gives several windows; all tails of a dataset are stacked and the
coefficients xi are found by sparse regression.
Test function: phi(t) = (1 - tau^2)^p, tau = (t - centre) / half_width in [-1, 1] (p = 2: phi and phi' are 0 at the ends).

Sparse regression: thresholded backward elimination (function `select`). The usual SINDy regression (STLSQ: drop
the terms with small coefficients) fails here even on noise-free curves, because library terms are nearly copies of
each other (while A > 0.1 dB, max(A - 0.1, 0) is A minus a constant): the fit splits the coefficient between the
copies and none is small. So a term is judged by how much the fit worsens without it, and among the terms that can
go, the least simple goes first. lambda (the threshold) is chosen by cross-validation over GROUPS of rows (events,
never single tails or windows): the largest lambda whose held-out error is within one standard error of the best.
Stability selection: the same regression on bootstrap resamples of the groups; the result is how often each term is
selected.
"""
import numpy as np


def term_names(switch_grid):
    return ["1", "A", "A^2", "A^3"] + [f"sw{g:g}" for g in switch_grid]


def library(a, switch_grid):
    """theta(A): one column per term (rows = samples)."""
    a = np.asarray(a, float)
    return np.column_stack([np.ones_like(a), a, a ** 2, a ** 3] + [np.maximum(a - g, 0.0) for g in switch_grid])


def test_function(half, power=2):
    """phi and dphi/dt on the 2 * half + 1 samples of a window (sample step = 1 time unit)."""
    tau = np.arange(-half, half + 1) / half
    return (1 - tau ** 2) ** power, -2 * power * tau * (1 - tau ** 2) ** (power - 1) / half


def fill_gaps(a, max_gap):
    """Linear interpolation over NaN runs of at most max_gap samples that have data on both sides."""
    a = np.array(a, float)
    bad = np.isnan(a)
    if not bad.any():
        return a
    idx = np.arange(len(a))
    good = idx[~bad]
    if len(good) < 2:
        return a
    filled = np.interp(idx, good, a[good])
    edges = np.flatnonzero(np.diff(np.r_[0, bad.astype(int), 0]))      # starts and ends of the NaN runs
    for s, e in zip(edges[::2], edges[1::2]):
        if e - s <= max_gap and s > 0 and e < len(a):
            a[s:e] = filled[s:e]
    return a


def weak_rows(a, half, stride, switch_grid, power=2, max_gap=3):
    """Weak-form rows of one tail. a: samples on a regular grid (step 1 minute, NaN allowed).
    Returns (b, G): one row per window that has no NaN left after gap filling."""
    a = fill_gaps(a, max_gap)
    phi, dphi = test_function(half, power)
    theta = library(a, switch_grid)
    b, G = [], []
    for c in range(half, len(a) - half, stride):
        seg = slice(c - half, c + half + 1)
        if np.isnan(a[seg]).any():
            continue
        b.append(-np.sum(dphi * a[seg]))                   # phi and phi' are 0 at the ends: the trapezoid rule is a sum
        G.append(phi @ theta[seg])
    return np.asarray(b), np.asarray(G).reshape(len(b), theta.shape[1])


def select(G, b, lam, ridge=1e-8):
    """Thresholded backward elimination, simplest terms kept first. Returns xi (0 for dropped terms), OLS on the rest.
    cost of a term = sqrt(increase of the squared residual when it alone is removed) / ||b||  (for independent
    columns this is |xi_j| ||G_j|| / ||b||, the usual size of a term; for a term that another one can replace it is
    about 0). While some term costs less than lambda, the LAST such term in the library order (1, A, A^2, A^3,
    switch terms from small to large A_c) is removed and the costs are computed again."""
    norm = np.linalg.norm(G, axis=0)
    active = norm > 1e-12 * max(norm.max(), 1e-300)
    bn = max(np.linalg.norm(b), 1e-300)
    while active.any():
        idx = np.flatnonzero(active)
        Gs = G[:, idx] / norm[idx]
        inv = np.linalg.inv(Gs.T @ Gs + ridge * np.eye(len(idx)))
        z = inv @ (Gs.T @ b)
        cost = np.abs(z) / np.sqrt(np.maximum(np.diag(inv), 1e-300)) / bn
        cheap = np.flatnonzero(cost < lam)
        if len(cheap) == 0:
            break
        active[idx[cheap.max()]] = False
    return refit(G, b, active)


def refit(G, b, mask):
    """Ordinary least squares on the selected terms only (no ridge, no shrinkage)."""
    xi = np.zeros(G.shape[1])
    if mask.any():
        xi[mask] = np.linalg.lstsq(G[:, mask], b, rcond=None)[0]
    return xi


def cv_lambda(G, b, groups, lambdas, folds, rng, ridge=1e-8):
    """Choose lambda by K-fold cross-validation over groups. Error of a fold = ||b - G xi||^2 / ||b||^2 on the
    held-out groups (xi refitted by OLS on the selected terms of the training groups). Returns (lambda, table):
    the largest lambda with a mean error <= best mean + its standard error (one-SE rule)."""
    ug = rng.permutation(np.unique(groups))
    fold_of = {g: i % folds for i, g in enumerate(ug)}
    f = np.array([fold_of[g] for g in groups])
    err = np.full((len(lambdas), folds), np.nan)
    for k in range(folds):
        tr, te = f != k, f == k
        if te.sum() == 0 or tr.sum() <= G.shape[1]:
            continue
        den = max(float(b[te] @ b[te]), 1e-300)
        for i, lam in enumerate(lambdas):
            xi = refit(G[tr], b[tr], select(G[tr], b[tr], lam, ridge) != 0)
            r = b[te] - G[te] @ xi
            err[i, k] = float(r @ r) / den
    mean = np.nanmean(err, axis=1)
    se = np.nanstd(err, axis=1, ddof=1) / np.sqrt(np.sum(np.isfinite(err), axis=1))
    best = int(np.nanargmin(mean))
    ok = np.flatnonzero(mean <= mean[best] + se[best])
    return float(lambdas[ok.max()]), {"lambda": list(map(float, lambdas)), "cv_mean": mean.tolist(), "cv_se": se.tolist()}


def stability(G, b, groups, lam, resamples, rng, ridge=1e-8):
    """Selection frequency of every term over bootstrap resamples of the groups (drawn with replacement)."""
    ug = np.unique(groups)
    rows_of = {g: np.flatnonzero(groups == g) for g in ug}
    hits = np.zeros(G.shape[1])
    for _ in range(resamples):
        ix = np.concatenate([rows_of[g] for g in rng.choice(ug, len(ug), replace=True)])
        hits += select(G[ix], b[ix], lam, ridge) != 0
    return hits / resamples


def _selftest():
    """Noise-free curves of known laws: the method must return the true terms and coefficients."""
    grid = [0.1, 0.2, 0.3, 0.5, 1.0]
    names = term_names(grid)
    t = np.arange(60.0)
    a0s = np.linspace(1.0, 4.0, 40)
    laws = {"exponential, dA/dt = -A / 30": lambda a0: a0 * np.exp(-t / 30),
            "quadratic, dA/dt = -0.05 A^2": lambda a0: 1 / (1 / a0 + 0.05 * t),
            "switch, dA/dt = -0.005 A - 0.04 max(A - 0.2, 0)": None}
    for name, fn in laws.items():
        bs, Gs = [], []
        for a0 in a0s:
            if fn is None:                                  # integrate the switch law with small Euler steps
                a, x = [], a0
                for _ in t:
                    a.append(x)
                    for _ in range(100):
                        x += 0.01 * (-0.005 * x - 0.04 * max(x - 0.2, 0))
                a = np.array(a)
            else:
                a = fn(a0)
            b_, G_ = weak_rows(a, 10, 5, grid)
            bs.append(b_), Gs.append(G_)
        b_, G_ = np.concatenate(bs), np.vstack(Gs)
        for lam in (0.02, 0.1):
            xi = select(G_, b_, lam)
            print(f"{name}, lambda {lam}: " + ", ".join(f"{n} {v:+.4f}" for n, v in zip(names, xi) if v != 0))


if __name__ == "__main__":
    _selftest()
