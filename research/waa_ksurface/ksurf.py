"""Library of steps A5-A7: the one-state WAA model with a learned rate k(R, W), and its fit.

The model, per link-event, one step per minute, W = 0 at the start of the window:

    v      = F_full(R(t))                          the steady WAA at the rain of this minute (per band)
    regime = wetting if v > W, else drying
    k      = exp( g_regime(R(t)) + h_regime(W) )   always positive
    W      <- v + (W - v) * exp(-k)                one exact exponential step of 1 min
    A - c  = W + L * P_full(R(t - d))

F_full, P_full = the MEASURED curves of step E4 inside their valid range, and above it a free, monotone, continuous
extension: linear in R between fixed knots, parameters = the increments from knot to knot (each >= 0).
g and h are piecewise linear between knots (g in ln(1 + R), h in W) and held outside them. The ladder:

    K0  k constant        K1  k(W)        K2  k(R)        K3  k = (function of R) * (function of W)

The loss is the mean Huber loss (0.3 dB = one quantization step) of the 5 min means of the model output against the 5 min
means of the measured A - c, over whole events, plus a smoothness penalty on the knot values. The gradient is exact (a
reverse pass through the minutes). Everything is defined in docs/waa_ksurface/decisions.md, step A5; settings: config `a5`.

    prob = prepare(D, band_idx, curves, rain, d)          the arrays the kernel reads (rain = the 1 min rain input, padded)
    m = Model(prob, "K0", "K1", lam)                      the wetting and the drying rung
    x = m.fit(m.start())                                  -> the parameter vector
    m.k(x, "drying", R, W)                                the learned rate (1 / min)
"""
import numpy as np
from numba import njit, prange
from scipy.optimize import minimize

from research.waa_ksurface import data
from research.wet_antenna.system_id import PADM                # minutes of rain stored before and after every window (data layer)

RUNGS = {"K0": (False, False), "K1": (False, True), "K2": (True, False), "K3": (True, True)}     # (depends on R, depends on W)
ORDER = ["K0", "K1", "K2", "K3"]                               # from the simplest


def config():
    return data.branch_config()["a5"]


# =================================================================================================== the curves F_full, P_full
def ext_coords(r, r_top, knots, hold):
    """Where a rain value sits in the extension above r_top: (segment j, fraction f). j = -1 inside the valid range.
    The extension value is  sum(increments before j) + f * increment[j].  hold = True: held above the last knot (F);
    False: f grows past 1 there, so the curve continues with its last slope (P)."""
    x = np.r_[r_top, np.asarray(knots, float)]
    r = np.asarray(r, float)
    j = np.clip(np.searchsorted(x, r, side="right") - 1, 0, len(x) - 2)
    f = (r - x[j]) / (x[j + 1] - x[j])
    if hold:
        f = np.minimum(f, 1.0)
    return np.where(r > r_top, j, -1).astype(np.int64), np.where(r > r_top, f, 0.0)


def full_curve(curve, kind, inc, r, cfg=None):
    """F_full or P_full of one band at the rain values r: the measured curve inside its valid range, the extension above.
    kind = "F" or "P"; inc = the increments of the extension."""
    cfg = cfg or config()
    r = np.asarray(r, float)
    top = curve.F_max if kind == "F" else curve.P_max
    base = (curve.F if kind == "F" else curve.P)(np.minimum(r, top))
    j, f = ext_coords(r, top, cfg["f_knots_mmh"] if kind == "F" else cfg["p_knots_mmh"], kind == "F")
    cum = np.r_[0.0, np.cumsum(inc)]
    jj = np.maximum(j, 0)
    return base + np.where(j >= 0, cum[jj] + f * np.asarray(inc, float)[jj], 0.0)


def increments_of(curve, kind, values_at_knots):
    """The increments that put the extension through the given values at its knots."""
    top = (curve.F if kind == "F" else curve.P)(np.array([curve.F_max if kind == "F" else curve.P_max]))[0]
    return np.diff(np.r_[top, np.asarray(values_at_knots, float)])


# =================================================================================================== the arrays of a fit
def prepare(D, band_idx, curves, rain, d=0, cfg=None):
    """The arrays the kernel reads. D = a dataset of link-events (data layer layout); band_idx = the band number of every
    link-event (the position of its band in `curves`); curves = [Curves per band]; rain = the 1 min rain input of every
    link-event with PADM minutes before and after its window (the layout of D["rain"]); d = the lag of the path term (min)."""
    cfg = cfg or config()
    n = D["le_n"]
    idx = np.concatenate([o + PADM + np.arange(k) for o, k in zip(D["le_roff"], n)])
    rep = np.repeat(np.arange(len(n)), n)
    r_w = np.nan_to_num(np.asarray(rain, float)[idx])                       # the rain that drives W
    r_p = np.nan_to_num(np.asarray(rain, float)[idx - int(d)])              # the rain of the path term, d minutes earlier
    b_min = np.asarray(band_idx)[rep]
    Ffix, Pfix = np.empty(len(idx)), np.empty(len(idx))
    jF, fF = np.empty(len(idx), np.int64), np.empty(len(idx))
    jP, fP = np.empty(len(idx), np.int64), np.empty(len(idx))
    for b, c in enumerate(curves):
        m = b_min == b
        Ffix[m] = c.F(np.minimum(r_w[m], c.F_max))
        Pfix[m] = c.P(np.minimum(r_p[m], c.P_max))
        jF[m], fF[m] = ext_coords(r_w[m], c.F_max, cfg["f_knots_mmh"], True)
        jP[m], fP[m] = ext_coords(r_p[m], c.P_max, cfg["p_knots_mmh"], False)
    uk = np.log1p(np.asarray(cfg["r_knots_mmh"], float))
    u = np.log1p(r_w)
    iR = np.clip(np.searchsorted(uk, u, side="right") - 1, 0, len(uk) - 2)
    fR = np.clip((u - uk[iR]) / (uk[iR + 1] - uk[iR]), 0.0, 1.0)
    rep5 = np.repeat(np.arange(len(n)), n // 5)
    return {"le_n": n.astype(np.int64), "boff1": D["le_boff1"].astype(np.int64), "boff5": D["le_boff5"].astype(np.int64),
            "band": np.asarray(band_idx, np.int64), "L": D["le_L"].astype(float), "event": D["le_event"].astype(np.int64),
            "Ffix": Ffix, "jF": jF, "fF": fF, "Pfix": Pfix, "jP": jP, "fP": fP, "iR": iR.astype(np.int64), "fR": fR,
            "target": (D["y5"] - D["le_c"][rep5]).astype(float), "w5": (D["w5"] > 0).astype(float), "nb": len(curves), "curves": curves, "d": int(d)}


# =================================================================================================== the kernel
@njit(cache=True, inline="always")
def _h(W, wk, h):
    """h(W): piecewise linear through (wk[0], 0), (wk[1], h[0]), ...; held outside the knots.
    Returns (value, segment i, fraction f, slope dh/dW)."""
    if len(h) == 0:
        return 0.0, 0, 0.0, 0.0
    last = len(wk) - 1
    if W <= wk[0]:
        return 0.0, 0, 0.0, 0.0
    if W >= wk[last]:
        return h[last - 1], last - 1, 1.0, 0.0
    i = 0
    while W >= wk[i + 1]:
        i += 1
    lo = 0.0 if i == 0 else h[i - 1]
    f = (W - wk[i]) / (wk[i + 1] - wk[i])
    return lo + f * (h[i] - lo), i, f, (h[i] - lo) / (wk[i + 1] - wk[i])


@njit(cache=True, parallel=True)
def _kernel(gw, hw, gd, hd, incF, incP, wkw, wkd, delta, lnk_lo, lnk_hi, want_grad, want_sim,
            le_n, boff1, boff5, band, L, Ffix, jF, fF, Pfix, jP, fP, iR, fR, target, w5, use):
    """The model on every link-event with use[e] > 0. Returns (loss sum per link-event, number of 5 min values per
    link-event, gradient of the loss sum per link-event, W per minute, model output per 5 min bin)."""
    nle = len(le_n)
    nGw, nHw, nGd, nHd = len(gw), len(hw), len(gd), len(hd)
    nb, nFk = incF.shape
    nPk = incP.shape[1]
    oHw, oGd, oHd = nGw, nGw + nHw, nGw + nHw + nGd
    oF = oHd + nHd
    oP = oF + nb * nFk
    npar = oP + nb * nPk
    cumF, cumP = np.zeros((nb, nFk + 1)), np.zeros((nb, nPk + 1))
    for b in range(nb):
        for j in range(nFk):
            cumF[b, j + 1] = cumF[b, j] + incF[b, j]
        for j in range(nPk):
            cumP[b, j + 1] = cumP[b, j] + incP[b, j]
    loss, cnt = np.zeros(nle), np.zeros(nle)
    G = np.zeros((nle if want_grad else 1, npar))
    W1 = np.zeros(len(Ffix) if want_sim else 1)
    out5 = np.zeros(len(target) if want_sim else 1)
    for e in prange(nle):
        if use[e] <= 0:
            continue
        n, o1, o5, b, Le = le_n[e], boff1[e], boff5[e], band[e], L[e]
        Ws, vs, ks, Es = np.empty(n + 1), np.empty(n), np.empty(n), np.empty(n)
        dres = np.zeros(n // 5)
        W, acc = 0.0, 0.0
        Ws[0] = 0.0
        for t in range(n):
            q = o1 + t
            v = Ffix[q]
            if jF[q] >= 0:
                v += cumF[b, jF[q]] + fF[q] * incF[b, jF[q]]
            if v > W:                                           # wetting
                lnk = (gw[iR[q]] + fR[q] * (gw[iR[q] + 1] - gw[iR[q]]) if nGw > 1 else gw[0]) + _h(W, wkw, hw)[0]
            else:                                               # drying
                lnk = (gd[iR[q]] + fR[q] * (gd[iR[q] + 1] - gd[iR[q]]) if nGd > 1 else gd[0]) + _h(W, wkd, hd)[0]
            k = np.exp(min(max(lnk, lnk_lo), lnk_hi))
            E = np.exp(-k)
            vs[t], ks[t], Es[t] = v, (k if (lnk > lnk_lo and lnk < lnk_hi) else 0.0), E      # k = 0 marks "clipped": no gradient through ln k
            W = v + (W - v) * E
            Ws[t + 1] = W
            p = Pfix[q]
            if jP[q] >= 0:
                p += cumP[b, jP[q]] + fP[q] * incP[b, jP[q]]
            acc += W + Le * p
            if want_sim:
                W1[q] = W
            if t % 5 == 4:
                jb = o5 + t // 5
                m = acc / 5.0
                acc = 0.0
                if want_sim:
                    out5[jb] = m
                if w5[jb] > 0:
                    r = m - target[jb]
                    ar = abs(r)
                    loss[e] += 0.5 * r * r if ar <= delta else delta * (ar - 0.5 * delta)
                    cnt[e] += 1.0
                    dres[t // 5] = (r if ar <= delta else delta * np.sign(r)) / 5.0
        if not want_grad:
            continue
        lam = 0.0                                               # d loss / d W after minute t, from the later minutes
        for t in range(n - 1, -1, -1):
            q = o1 + t
            gm = dres[t // 5]                                   # d loss / d (model output of minute t)
            if gm != 0.0 and jP[q] >= 0:                        # the P extension
                for j in range(jP[q]):
                    G[e, oP + b * nPk + j] += gm * Le
                G[e, oP + b * nPk + jP[q]] += gm * Le * fP[q]
            tot = lam + gm
            if tot == 0.0:
                lam = 0.0
                continue
            W, v, k, E = Ws[t], vs[t], ks[t], Es[t]
            dv = tot * (1.0 - E)
            dlnk = -tot * (W - v) * E * k
            wet = v > W
            if wet:
                _, i, f, slope = _h(W, wkw, hw)
            else:
                _, i, f, slope = _h(W, wkd, hd)
            lam = tot * E + dlnk * slope
            if dv != 0.0 and jF[q] >= 0:                        # the F extension
                for j in range(jF[q]):
                    G[e, oF + b * nFk + j] += dv
                G[e, oF + b * nFk + jF[q]] += dv * fF[q]
            if dlnk != 0.0:
                og, ng, oh, nh = (0, nGw, oHw, nHw) if wet else (oGd, nGd, oHd, nHd)
                if ng > 1:
                    G[e, og + iR[q]] += dlnk * (1.0 - fR[q])
                    G[e, og + iR[q] + 1] += dlnk * fR[q]
                else:
                    G[e, og] += dlnk
                if nh > 0:
                    if i >= 1:
                        G[e, oh + i - 1] += dlnk * (1.0 - f)
                    G[e, oh + i] += dlnk * f
    return loss, cnt, G, W1, out5


def _second_diff(n):
    """The matrix that takes n knot values to their second differences (empty if n < 3)."""
    M = np.zeros((max(n - 2, 0), n))
    for i in range(n - 2):
        M[i, i:i + 3] = (1.0, -2.0, 1.0)
    return M


# =================================================================================================== the model
class Model:
    """One rung for wetting and one for drying on a prepared dataset. The parameter vector is
    [g_wet | h_wet | g_dry | h_dry | F increments (band by band) | P increments (band by band)]."""

    def __init__(self, prob, rung_wet="K0", rung_dry="K0", lam=0.0, cfg=None):
        self.p, self.cfg, self.lam = prob, cfg or config(), float(lam)
        self.rung = {"wetting": rung_wet, "drying": rung_dry}
        c = self.cfg
        self.rk = np.asarray(c["r_knots_mmh"], float)
        self.wk = {"wetting": np.asarray(c["w_knots_wet_db"], float), "drying": np.asarray(c["w_knots_dry_db"], float)}
        self.nF, self.nP, self.nb = len(c["f_knots_mmh"]), len(c["p_knots_mmh"]), prob["nb"]
        sizes = []
        for reg in ("wetting", "drying"):
            on_r, on_w = RUNGS[self.rung[reg]]
            sizes += [len(self.rk) if on_r else 1, len(self.wk[reg]) - 1 if on_w else 0]
        sizes += [self.nb * self.nF, self.nb * self.nP]
        edges = np.r_[0, np.cumsum(sizes)]
        self.sl = {name: slice(edges[i], edges[i + 1]) for i, name in enumerate(("gw", "hw", "gd", "hd", "F", "P"))}
        self.n = int(edges[-1])
        lo_k, hi_k = np.log(1.0 / c["tau_max_min"]), np.log(1.0 / c["tau_min_min"])
        self.lnk = (lo_k, hi_k)
        self.lo, self.hi = np.empty(self.n), np.empty(self.n)
        for name, (a, b) in (("gw", (lo_k, hi_k)), ("gd", (lo_k, hi_k)), ("hw", (-c["h_bound"], c["h_bound"])), ("hd", (-c["h_bound"], c["h_bound"])),
                             ("F", (0.0, c["f_inc_max_db"])), ("P", (0.0, c["p_inc_max_db_km"]))):
            self.lo[self.sl[name]], self.hi[self.sl[name]] = a, b
        # the smoothness penalty: second differences of g and of (0, h) per regime
        self.pen = []
        for g, h in (("gw", "hw"), ("gd", "hd")):
            ng, nh = self.sl[g].stop - self.sl[g].start, self.sl[h].stop - self.sl[h].start
            if ng >= 3:
                self.pen.append((self.sl[g], _second_diff(ng)))
            if nh >= 2:
                self.pen.append((self.sl[h], _second_diff(nh + 1)[:, 1:]))      # the first value of (0, h) is the fixed 0
        self.fixed = np.zeros(self.n, bool)

    # ---- parameters --------------------------------------------------------------------------------------------------
    def start(self, inc_F=None, inc_P=None):
        """Starting values: 1 / k = k_start_min for both regimes, no dependence on R or W; the extension as given
        (default: F almost held, P with the last slope of the measured curve)."""
        x = np.zeros(self.n)
        x[self.sl["gw"]] = x[self.sl["gd"]] = np.log(1.0 / self.cfg["k_start_min"])
        if inc_F is None:
            inc_F = [np.full(self.nF, self.cfg["f_inc_start_db"]) for _ in range(self.nb)]
        if inc_P is None:
            inc_P = [increments_of(c, "P", c.P(np.asarray(self.cfg["p_knots_mmh"], float))) for c in self.p["curves"]]
        x[self.sl["F"]], x[self.sl["P"]] = np.concatenate(inc_F), np.concatenate(inc_P)
        return np.clip(x, self.lo, self.hi)

    def hold(self, *names):
        """Hold blocks of parameters at the values they have in the starting vector ("gw", "hw", "gd", "hd", "F", "P")."""
        for name in names:
            self.fixed[self.sl[name]] = True
        return self

    def parts(self, x):
        s = self.sl
        return (x[s["gw"]], x[s["hw"]], x[s["gd"]], x[s["hd"]], x[s["F"]].reshape(self.nb, self.nF), x[s["P"]].reshape(self.nb, self.nP))

    def carry(self, other, x_other):
        """A starting vector for this model from the fit of another rung (constants are spread over the knots; W parts that
        the other rung does not have start at 0; the extension is copied)."""
        x = self.start()
        for g, h in (("gw", "hw"), ("gd", "hd")):
            go, ho = x_other[other.sl[g]], x_other[other.sl[h]]
            ng, nh = self.sl[g].stop - self.sl[g].start, self.sl[h].stop - self.sl[h].start
            x[self.sl[g]] = go if len(go) == ng else (np.full(ng, go.mean() if len(go) == 1 else np.median(go)))
            x[self.sl[h]] = ho if len(ho) == nh else 0.0
        x[self.sl["F"]], x[self.sl["P"]] = x_other[other.sl["F"]], x_other[other.sl["P"]]
        return np.clip(x, self.lo, self.hi)

    # ---- loss --------------------------------------------------------------------------------------------------------
    def _run(self, x, use, grad, sim=False):
        p = self.p
        gw, hw, gd, hd, iF, iP = self.parts(np.asarray(x, float))
        return _kernel(gw, hw, gd, hd, np.ascontiguousarray(iF), np.ascontiguousarray(iP), self.wk["wetting"], self.wk["drying"], self.cfg["huber_delta_db"],
                       self.lnk[0], self.lnk[1], grad, sim, p["le_n"], p["boff1"], p["boff5"], p["band"], p["L"], p["Ffix"], p["jF"], p["fF"],
                       p["Pfix"], p["jP"], p["fP"], p["iR"], p["fR"], p["target"], p["w5"], use)

    def loss(self, x, use=None, penalty=True):
        """(loss, gradient): the mean Huber loss of the link-events with use > 0 (use = how often each counts), plus the penalty."""
        use = np.ones(len(self.p["le_n"])) if use is None else np.asarray(use, float)
        ls, cnt, G, _, _ = self._run(x, use, True)
        N = float((use * cnt).sum())
        val, g = float((use * ls).sum() / N), (use[:, None] * G).sum(axis=0) / N
        if penalty and self.lam > 0:
            for sl, M in self.pen:
                dd = M @ x[sl]
                val += self.lam * float(dd @ dd)
                g[sl] += 2.0 * self.lam * (M.T @ dd)
        return val, g

    def per_link_event(self, x):
        """(Huber loss sum, number of 5 min values) per link-event: for held-out scores and event bootstraps."""
        ls, cnt, _, _, _ = self._run(x, np.ones(len(self.p["le_n"])), False)
        return ls, cnt

    def simulate(self, x):
        """(W per minute, model output per 5 min bin) for every link-event."""
        _, _, _, W1, out5 = self._run(x, np.ones(len(self.p["le_n"])), False, sim=True)
        return W1, out5

    def fit(self, x0, use=None, maxiter=None):
        """Minimise the loss from x0 (L-BFGS-B inside the bounds; held parameters stay at their x0 value)."""
        lo, hi = np.where(self.fixed, x0, self.lo), np.where(self.fixed, x0, self.hi)
        r = minimize(lambda z: self.loss(z, use), np.clip(x0, lo, hi), jac=True, method="L-BFGS-B", bounds=list(zip(lo, hi)),
                     options={"maxiter": maxiter or self.cfg["maxiter"], "ftol": 1e-10, "gtol": 1e-7})
        self.last = {"loss": float(r.fun), "nit": int(r.nit), "nfev": int(r.nfev), "converged": bool(r.success)}
        return r.x

    # ---- what was learned ----------------------------------------------------------------------------------------------
    def k(self, x, regime, R, W):
        """The learned rate k (1 / min) of a regime at rain R (mm/h) and level W (dB)."""
        g = x[self.sl["gw" if regime == "wetting" else "gd"]]
        h = x[self.sl["hw" if regime == "wetting" else "hd"]]
        R, W = np.broadcast_arrays(np.asarray(R, float), np.asarray(W, float))
        lnk = np.interp(np.log1p(R), np.log1p(self.rk), g) if len(g) > 1 else np.full(R.shape, g[0])
        if len(h):
            lnk = lnk + np.interp(W, self.wk[regime], np.r_[0.0, h])
        return np.exp(np.clip(lnk, *self.lnk))

    def extension(self, x):
        """{("F" | "P", band index): the values of the extension at its knots}."""
        _, _, _, _, iF, iP = self.parts(x)
        out = {}
        for b, c in enumerate(self.p["curves"]):
            out[("F", b)] = full_curve(c, "F", iF[b], self.cfg["f_knots_mmh"], self.cfg)
            out[("P", b)] = full_curve(c, "P", iP[b], self.cfg["p_knots_mmh"], self.cfg)
        return out
