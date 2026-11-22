"""Library of step 8: closed forms for the drying rate k_dry(W) (decisions.md, "Step 8: naming the drying law").

    Wp, w = points()                                  the A4 drying-domain W bins (centres) and their density weights (mean 1)
    y = k1_lnk(x, Wp)                                 ln k_dry of an A7 K1 parameter vector at those points
    f = physics()["sat"]                              a form: k = f.k(W, theta); theta = its fitted numbers
    theta, wrss = fit_curve(f, Wp, y, w)              weighted least squares of ln k(form) - y, best of the form's starts
    aicc(wrss, n, p)                                  the AICc of a weighted Gaussian fit (sigma counted)
    fm = FormModel(prob, f, wet)                      the A7 model with k_dry = f(W; theta) in place of the K1 spline
    z = fm.fit(fm.start(theta, x_a7), use)            refit theta and the F, P extensions (wetting held)
    S, N = held_out(fm, z, use, fold, score, prob)    5-fold held-out Huber sums per event on the scored values

How a form enters the model: ln k(W; theta) is tabulated on a fine geometric grid of W and handed to the UNCHANGED kernel of
`ksurf` as the knots of a K1 drying rate (piecewise linear in W, held above the last knot). The kernel's exact gradient in
the knot values is chained to theta through the derivative of the table.
"""
import numpy as np
import pandas as pd
import sympy as sp
from scipy.optimize import least_squares, minimize

from research.waa_ksurface import data, ksurf, ladder
from research.waa_ksurface import paths as P

C8 = data.branch_config()["a8"]
TINY = 1e-12                                                    # a form with k <= 0 gets ln k = ln(TINY) (it loses the fit)
GRID = np.r_[0.0, np.geomspace(C8["grid_lo_db"], C8["grid_hi_db"], C8["grid_points"])]     # the W knots of a tabulated form
GRID_EVAL = np.r_[C8["grid_w0_db"], GRID[1:]]                   # the value at W = 0 is the form at a tiny W
W = sp.Symbol("W")


# =================================================================================================== the points
def points():
    """(W at the bin centres, weights scaled to mean 1, minutes) of the frozen A4 drying domain, pooled over the bands."""
    d = pd.read_csv(P.A4_DOMAIN)
    d = d[(d.band == "pooled") & (d.regime == "drying") & d.in_domain]
    g = d.groupby(["wbin", "W_lo", "W_hi"]).minutes_L2.sum().reset_index().sort_values("wbin")
    m = g.minutes_L2.to_numpy(float)
    return 0.5 * (g.W_lo.to_numpy() + g.W_hi.to_numpy()), m / m.mean(), m


def k1_model(nb=2):
    """A K1 drying model shell for reading k_dry from a stored parameter vector (Model.k does not read the data)."""
    return ksurf.Model({"nb": nb}, "K0", "K1")


def k1_lnk(x, Wp, nb=2):
    m = k1_model(nb)
    assert m.n == len(x), "the parameter vector does not fit a K0 / K1 model with this number of bands"
    return np.log(m.k(np.asarray(x, float), "drying", 0.0, np.asarray(Wp, float)))


# =================================================================================================== the forms
class Form:
    """A closed form k(W; theta). expr = a sympy expression of k in W and the symbols `params` (= theta, in the order
    given); natural(theta) -> {name: value} the numbers as written in the formula; starts = starting thetas."""

    def __init__(self, name, expr, params, starts, text, natural=None, bounds=None, source="physics"):
        self.name, self.expr, self.params, self.text, self.source = name, expr, list(params), text, source
        self.p = len(self.params)
        self.starts = [np.asarray(s, float) for s in starts]
        self.natural = natural or (lambda t: {str(s): float(v) for s, v in zip(self.params, t)})
        self.bounds = bounds or [(None, None)] * self.p
        self._f = sp.lambdify([W] + self.params, expr, "numpy")

    def k(self, Wv, theta):
        with np.errstate(all="ignore"):
            v = np.asarray(self._f(np.asarray(Wv, float), *np.asarray(theta, float)), float) * np.ones(np.shape(Wv))
        return v

    def lnk(self, Wv, theta):
        v = self.k(Wv, theta)
        return np.log(np.where(np.isfinite(v) & (v > TINY), v, TINY))

    def formula(self, theta):
        """The formula with its numbers (natural parameters), as text."""
        nat = self.natural(theta)
        return self.text + "  with  " + ", ".join(f"{k} = {v:.4g}" for k, v in nat.items())


def physics():
    """The four physics forms of step 8 (k in 1 / min, W in dB). Positive parameters through their logarithm."""
    t = sp.symbols("t0:3")
    e = sp.exp
    L = np.log
    return {
        "linear": Form("linear", e(t[0]) + t[1] * W, t[:2], [(L(a), b) for a in (0.01, 0.05, 0.2) for b in (-0.02, 0.0, 0.05)],
                       "k = a + b*W", lambda th: {"a (1/min)": np.exp(th[0]), "b (1/(min dB))": th[1]}, [(-30, 5), (-10, 10)]),
        "sat": Form("sat", e(t[0]) * W / (W + e(t[1])), t[:2], [(L(a), L(h)) for a in (0.05, 0.2, 1.0) for h in (0.1, 0.5, 2.0, 5.0)],
                    "k = k_max*W/(W + W_h)", lambda th: {"k_max (1/min)": np.exp(th[0]), "W_h (dB)": np.exp(th[1])}, [(-30, 5), (-30, 10)]),
        "power": Form("power", e(t[0]) * W ** t[1], t[:2], [(L(a), q) for a in (0.02, 0.1, 0.3) for q in (0.2, 0.5, 1.0, 2.0)],
                      "k = a*W^p", lambda th: {"a (1/(min dB^p))": np.exp(th[0]), "p": th[1]}, [(-30, 5), (-10, 10)]),
        "sat_floor": Form("sat_floor", e(t[0]) + e(t[1]) * W / (W + e(t[2])), t,
                          [(L(k0), L(a), L(h)) for k0 in (0.002, 0.01, 0.03) for a in (0.05, 0.2, 1.0) for h in (0.1, 0.5, 2.0, 5.0)],
                          "k = k_0 + k_max*W/(W + W_h)", lambda th: {"k_0 (1/min)": np.exp(th[0]), "k_max (1/min)": np.exp(th[1]), "W_h (dB)": np.exp(th[2])},
                          [(-30, 5), (-30, 5), (-30, 10)]),
    }


def from_sympy(name, expr_text):
    """A PySR expression (text in W) as a Form: every floating constant becomes a parameter, started at its value."""
    ex = sp.sympify(expr_text, locals={"W": W})
    consts = sorted({a for a in ex.atoms(sp.Float)}, key=lambda a: float(a))
    syms = sp.symbols(f"c0:{len(consts)}") if consts else []
    ex2 = ex.xreplace(dict(zip(consts, syms)))
    return Form(name, ex2, syms, [[float(c) for c in consts]], f"k = {sp.sstr(ex2)}", source="PySR")


# =================================================================================================== the curve fit and AICc
def fit_curve(form, Wp, y, w, starts=None):
    """(theta, weighted residual sum of squares) of ln k(form) - y, the best over the starts."""
    best = None
    sw = np.sqrt(np.asarray(w, float))
    for s in (starts if starts is not None else form.starts):
        if form.p == 0:
            r = np.sum(np.asarray(w) * (form.lnk(Wp, []) - y) ** 2)
            return np.zeros(0), float(r)
        lo = np.array([b[0] if b[0] is not None else -np.inf for b in form.bounds])
        hi = np.array([b[1] if b[1] is not None else np.inf for b in form.bounds])
        s = np.clip(np.asarray(s, float), lo + 1e-9, hi - 1e-9)
        try:
            r = least_squares(lambda t: sw * (form.lnk(Wp, t) - y), s, bounds=(lo, hi), method="trf", max_nfev=C8["lsq_max_nfev"])
        except ValueError:                                      # a start where the form is not finite
            continue
        if best is None or r.cost < best.cost:
            best = r
    return best.x, float(2.0 * best.cost)


def aicc(wrss, n, p):
    K = p + 1                                                   # sigma is a parameter too
    return float(n * np.log(wrss / n) + 2 * K + 2 * K * (K + 1) / (n - K - 1))


def rank(tab):
    """tab: one row per candidate with columns name, p, aicc_gauge, aicc_radar. Adds the deltas, Delta_max, the chosen form,
    the Pareto front of (p, Delta_max) and the trial order (decisions.md, step 8 (2))."""
    t = tab.copy()
    for e in ("gauge", "radar"):
        t[f"daicc_{e}"] = t[f"aicc_{e}"] - t[f"aicc_{e}"].min()
    t["dmax"] = t[["daicc_gauge", "daicc_radar"]].max(axis=1)
    t["aicc_sum"] = t.aicc_gauge + t.aicc_radar
    best = t.sort_values(["dmax", "aicc_sum"]).index[0]
    near = t[t.dmax <= t.dmax[best] + C8["aicc_within"]]
    chosen = near.sort_values(["p", "dmax"]).index[0]
    t["best"], t["chosen"] = t.index == best, t.index == chosen
    t["front"] = [not ((t.p <= r.p) & (t.dmax < r.dmax) & (t.index != i)).any() for i, r in t.iterrows()]
    rest = t[t.front & (t.index != chosen) & (t.dmax <= C8["dmax_try"])].sort_values(["dmax", "p"])
    order = [chosen] + list(rest.index)
    order = order[:C8["max_trials"]]
    t["trial"] = [order.index(i) if i in order else -1 for i in t.index]
    return t, order


# =================================================================================================== the plug-back model
class FormModel:
    """The A7 model (wetting K0 held at `wet` = (gw, hw)) with the drying rate k_dry = form(W; theta). The free vector is
    z = [theta | F increments | P increments]; F, P as in ksurf.Model (measured parts fixed, extensions free)."""

    def __init__(self, prob, form, wet, cfg=None):
        c = dict(cfg or ksurf.config())
        c["w_knots_dry_db"] = [float(v) for v in GRID]
        self.m = ksurf.Model(prob, "K0", "K1", 0.0, cfg=c)      # no smoothness penalty: a form has no knots
        self.form, self.cfg = form, c
        x = self.m.start()
        x[self.m.sl["gw"]], x[self.m.sl["hw"]] = wet
        self.base = x
        self.nF, self.nP = self.m.sl["F"].stop - self.m.sl["F"].start, self.m.sl["P"].stop - self.m.sl["P"].start
        self.bounds = list(form.bounds) + list(zip(self.m.lo[self.m.sl["F"]], self.m.hi[self.m.sl["F"]])) + list(zip(self.m.lo[self.m.sl["P"]], self.m.hi[self.m.sl["P"]]))

    # ---- parameters ------------------------------------------------------------------------------------------------
    def start(self, theta, x_ext):
        """z from theta and the F, P extension of a ksurf parameter vector (e.g. the A7 fit at this estimate)."""
        src = k1_model(self.m.nb)
        x_ext = np.asarray(x_ext, float)
        return np.r_[np.asarray(theta, float), x_ext[src.sl["F"]], x_ext[src.sl["P"]]]

    def table(self, theta):
        """ln k at the grid, clipped to the kernel's bounds of ln k."""
        return np.clip(self.form.lnk(GRID_EVAL, theta), *self.m.lnk)

    def xfull(self, z):
        p = self.form.p
        L = self.table(z[:p])
        x = self.base.copy()
        x[self.m.sl["gd"]] = L[0]
        x[self.m.sl["hd"]] = L[1:] - L[0]
        x[self.m.sl["F"]] = z[p:p + self.nF]
        x[self.m.sl["P"]] = z[p + self.nF:]
        return x

    def jac_table(self, theta):
        """d table / d theta (grid x p), central differences of the (cheap) form."""
        theta = np.asarray(theta, float)
        J = np.zeros((len(GRID), len(theta)))
        for j in range(len(theta)):
            h = C8["fd_step"] * max(1.0, abs(theta[j]))
            a, b = theta.copy(), theta.copy()
            a[j] += h
            b[j] -= h
            J[:, j] = (self.table(a) - self.table(b)) / (2 * h)
        return J

    # ---- loss and fit ------------------------------------------------------------------------------------------------
    def loss(self, z, use=None):
        z = np.asarray(z, float)
        p = self.form.p
        val, g = self.m.loss(self.xfull(z), use, penalty=False)
        gL = np.r_[g[self.m.sl["gd"]][0] - g[self.m.sl["hd"]].sum(), g[self.m.sl["hd"]]]
        gt = self.jac_table(z[:p]).T @ gL if p else np.zeros(0)
        return val, np.r_[gt, g[self.m.sl["F"]], g[self.m.sl["P"]]]

    def fit(self, z0, use=None, maxiter=None):
        r = minimize(lambda z: self.loss(z, use), np.asarray(z0, float), jac=True, method="L-BFGS-B", bounds=self.bounds,
                     options={"maxiter": maxiter or self.cfg["maxiter"], "ftol": 1e-10, "gtol": 1e-7})
        self.last = {"loss": float(r.fun), "nit": int(r.nit), "nfev": int(r.nfev), "converged": bool(r.success)}
        return r.x

    def simulate(self, z):
        return self.m.simulate(self.xfull(np.asarray(z, float)))

    def per_link_event(self, z):
        return self.m.per_link_event(self.xfull(np.asarray(z, float)))


# =================================================================================================== held-out sums
def held_out(model, x, use, fold, score, prob, maxiter=None):
    """Per-event (Huber sum, count) of the held-out scored values, each fold refitted from x (the code of a7_review.stage_w).
    model = a ksurf.Model (with its held parts) or a FormModel."""
    bin_le = np.repeat(np.arange(len(prob["le_n"])), prob["le_n"] // 5)
    ev = prob["event"]
    n_ev = int(ev.max()) + 1
    S, N = np.zeros(n_ev), np.zeros(n_ev)
    use = np.asarray(use, float)
    for f in range(fold.max() + 1):
        out = (fold == f) & (use > 0)
        xf = model.fit(x, use * (fold != f), maxiter=maxiter or ladder.C["maxiter_fold"])
        b = out[bin_le] & score
        h = ladder._huber(model.simulate(xf)[1][b] - prob["target"][b], ksurf.config()["huber_delta_db"])
        S += np.bincount(ev[bin_le[b]], weights=h, minlength=n_ev)
        N += np.bincount(ev[bin_le[b]], minlength=n_ev).astype(float)
    return S, N
