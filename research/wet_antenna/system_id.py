"""Library of steps S1 and S2: the wet-antenna models with rain as input (M0-M4; for S2 also a level drift and the
discovery rung M5), their fit, the link-event datasets and the leave-one-event-out judge.

The model of one link-event (one event on one sub-link), minute by minute:

    A_pred(t) = W(t) + L * P_band(R(t - d)) + c

    R        rain on the path (mm/h), held constant inside its bin; d = one lag in whole minutes for all links
    P_band   path attenuation per km: piecewise linear in R on fixed knots, 0 at R = 0, never falling
    F_band   the WAA level the rain drives toward (dB): the same kind of curve
    a_link   one amplitude per link that multiplies F (antennas differ)
    c        the offset of the link-event (measured before the rain, not fitted)
    W        M0: 0.   M1: a F(R).   M2: dW/dt = (a F - W) / tau.   M3: tau_wet while a F > W, tau_dry otherwise.
             M4: two stores, W = W_f + W_s: the fast one gets phi a F (with tau_f_wet / tau_f_dry), the slow one
             (1 - phi) a F with tau_s.
Inside one minute the input is constant, so every store is updated exactly: W <- u + (W - u) exp(-1 / tau).
The loss is the Huber loss of the residuals of the 5 min means (or of the 1 min values), summed over the points.
`kernel` computes the loss and, by going backwards through the minutes, its exact gradient (numba).
"""
import warnings

import numpy as np
import pandas as pd
from numba import njit
from scipy.optimize import minimize

from research.wet_antenna import paths as PATHS
from research.wet_antenna import sindy
from research.wet_antenna.simulator import ar1, quantize
from research.wet_antenna.tails import itu_p838

PADM = 10                                                  # minutes of rain stored before and after every window (= the largest lag)
MODELS = ["M0", "M1", "M2", "M3", "M4"]
MODE = {"M0": 0, "M1": 1, "M2": 2, "M3": 2, "M4": 2}       # 0 no WAA, 1 static, 2 with memory
DYN = {"M0": [], "M1": [], "M2": ["tau"], "M3": ["tau_wet", "tau_dry"], "M4": ["phi", "tau_f_wet", "tau_f_dry", "tau_s"]}


@njit(cache=True)
def kernel(mode, d, Finc, Pinc, a, tau, phi, wle, le_roff, le_n, le_boff, le_link, le_band, le_L, le_c, rain, y, w, agg,
           delta, knots, pred, le_loss, gF, gP, ga, gtau, gphi, want_grad):
    """Loss of all link-events (weights wle) and, if want_grad, its gradient. Fills `pred` (one value per output bin
    of `agg` minutes) and `le_loss` (the unweighted loss of every link-event). tau[band] = (fast wet, fast dry, slow)."""
    nb, nk = Finc.shape
    cF, cP = np.zeros((nb, nk)), np.zeros((nb, nk))        # value of the curve at the lower knot of each segment
    for b in range(nb):
        for j in range(1, nk):
            cF[b, j] = cF[b, j - 1] + Finc[b, j - 1]
            cP[b, j] = cP[b, j - 1] + Pinc[b, j - 1]
    maxn = le_n.max()
    Wf, Ws, fr = np.empty(maxn), np.empty(maxn), np.empty(maxn)
    jj, br = np.empty(maxn, np.int64), np.empty(maxn, np.int64)
    psi = np.empty(maxn)
    G0F, G1F, G0P, G1P = np.zeros((nb, nk)), np.zeros((nb, nk)), np.zeros((nb, nk)), np.zeros((nb, nk))
    total = 0.0
    for e in range(len(le_n)):
        n, b, link, L, c = le_n[e], le_band[e], le_link[e], le_L[e], le_c[e]
        al, ph = a[link], phi[b]
        ro, bo = le_roff[e] + PADM - d, le_boff[e]
        kfw, kfd, ks = np.exp(-1.0 / tau[b, 0]), np.exp(-1.0 / tau[b, 1]), np.exp(-1.0 / tau[b, 2])
        wf, ws, acc = 0.0, 0.0, 0.0
        for t in range(n):                                 # forward: simulate the minutes
            r = rain[ro + t]
            j, f_ = nk - 1, 1.0
            if r <= 0.0:
                j, f_ = 0, 0.0
            else:
                for q in range(nk):
                    if r < knots[q + 1]:
                        j, f_ = q, (r - knots[q]) / (knots[q + 1] - knots[q])
                        break
            jj[t], fr[t] = j, f_
            wtot = 0.0
            if mode == 1:
                wtot = al * (cF[b, j] + f_ * Finc[b, j])
            elif mode == 2:
                u = al * (cF[b, j] + f_ * Finc[b, j])
                uf, us = ph * u, (1.0 - ph) * u
                Wf[t], Ws[t] = wf, ws
                if uf > wf:
                    br[t] = 0
                    wf = uf + (wf - uf) * kfw
                else:
                    br[t] = 1
                    wf = uf + (wf - uf) * kfd
                ws = us + (ws - us) * ks
                wtot = wf + ws
            acc += wtot + L * (cP[b, j] + f_ * Pinc[b, j]) + c
            if (t + 1) % agg == 0:
                pred[bo + t // agg] = acc / agg
                acc = 0.0
        loss = 0.0
        for q in range(n // agg):                          # Huber loss of the output bins
            psi[q] = 0.0
            if w[bo + q] > 0:
                res = pred[bo + q] - y[bo + q]
                if abs(res) <= delta:
                    loss += 0.5 * res * res
                    psi[q] = res
                else:
                    loss += delta * (abs(res) - 0.5 * delta)
                    psi[q] = delta if res > 0 else -delta
        le_loss[e] = loss
        total += wle[e] * loss
        if want_grad and wle[e] != 0.0:                    # backward: the gradient
            lf, ls, gkfw, gkfd, gks = 0.0, 0.0, 0.0, 0.0, 0.0
            for t in range(n - 1, -1, -1):
                gp = wle[e] * psi[t // agg] / agg
                j, f_ = jj[t], fr[t]
                G0P[b, j] += gp * L
                G1P[b, j] += gp * L * f_
                if mode == 0:
                    continue
                fval = cF[b, j] + f_ * Finc[b, j]
                gu = gp
                if mode == 2:
                    u = al * fval
                    uf, us = ph * u, (1.0 - ph) * u
                    lf += gp
                    ls += gp
                    k = kfw if br[t] == 0 else kfd
                    guf = lf * (1.0 - k)
                    if br[t] == 0:
                        gkfw += lf * (Wf[t] - uf)
                    else:
                        gkfd += lf * (Wf[t] - uf)
                    lf *= k
                    gus = ls * (1.0 - ks)
                    gks += ls * (Ws[t] - us)
                    ls *= ks
                    gu = ph * guf + (1.0 - ph) * gus
                    gphi[b] += u * (guf - gus)
                ga[link] += gu * fval
                G0F[b, j] += gu * al
                G1F[b, j] += gu * al * f_
            gtau[b, 0] += gkfw * kfw / tau[b, 0] ** 2
            gtau[b, 1] += gkfd * kfd / tau[b, 1] ** 2
            gtau[b, 2] += gks * ks / tau[b, 2] ** 2
    if want_grad:                                          # a segment's increment also lifts every segment above it
        for b in range(nb):
            sF, sP = 0.0, 0.0
            for j in range(nk - 1, -1, -1):
                gF[b, j] = sF + G1F[b, j]
                gP[b, j] = sP + G1P[b, j]
                sF += G0F[b, j]
                sP += G0P[b, j]
    return total


@njit(cache=True)
def simulate_stores(u, phi, tau_f_wet, tau_f_dry, tau_s):
    """W per minute for a given drive u = a F(R) per minute (the truth of the fake data; same update as `kernel`)."""
    kfw, kfd, ks = np.exp(-1.0 / tau_f_wet), np.exp(-1.0 / tau_f_dry), np.exp(-1.0 / tau_s)
    out, wf, ws = np.empty(len(u)), 0.0, 0.0
    for t in range(len(u)):
        uf, us = phi * u[t], (1.0 - phi) * u[t]
        wf = uf + (wf - uf) * (kfw if uf > wf else kfd)
        ws = us + (ws - us) * ks
        out[t] = wf + ws
    return out


class Problem:
    """One model on one dataset: the parameter vector, its bounds, the loss and its gradient.
    D = the dataset (dict of arrays, see system_id_ladder.build_dataset); M = config `wet_antenna.system_id.model`.
    Options: agg (5 or 1 minutes per output bin), no_path (P = 0), tau_per_band (own time constants per band)."""

    def __init__(self, D, model, M, agg=5, no_path=False, tau_per_band=False):
        self.D, self.model, self.M, self.agg, self.no_path = D, model, M, agg, no_path or False
        self.mode, self.nb, self.nl, self.nk = MODE[model], int(D["nb"]), int(D["nl"]), len(M["knots_mmh"]) - 1
        self.knots = np.asarray(M["knots_mmh"], float)
        self.nt = self.nb if tau_per_band else 1           # copies of every dynamic parameter
        self.y, self.w, self.boff = (D["y5"], D["w5"], D["le_boff5"]) if agg == 5 else (D["y1"], D["w1"], D["le_boff1"])
        self.pen = 0.5 * (M["huber_delta_db"] / M["a_prior_sd"]) ** 2
        B = M["bounds"]
        names, lo, hi = [], [], []

        def add(name, n, bound):
            names.extend([name] * n), lo.extend([bound[0]] * n), hi.extend([bound[1]] * n)

        if not self.no_path:
            add("logP", self.nb * self.nk, B["log_increment"])
        if self.mode > 0:
            add("logF", self.nb * self.nk, B["log_increment"])
            add("loga", self.nl, B["log_a"])
        for p in DYN[model]:
            add(p, self.nt, B["logit_phi"] if p == "phi" else list(np.log(B["tau_min"])))
        self.names, self.bounds = np.array(names), list(zip(lo, hi))
        self.pred, self.le_loss = np.zeros(len(self.y)), np.zeros(len(D["le_n"]))

    def start(self, s):
        """A generic start: F and P straight lines up to f_top / p_top at the last knot, a = 1, the start's time constants."""
        x = np.zeros(len(self.names))
        width = np.diff(self.knots) / self.knots[-1]
        x[self.names == "logP"] = np.tile(np.log(s["p_top"] * width), self.nb)[:(self.names == "logP").sum()]
        x[self.names == "logF"] = np.tile(np.log(s["f_top"] * width), self.nb)[:(self.names == "logF").sum()]
        for p in DYN[self.model]:
            x[self.names == p] = np.log(s[p] / (1 - s[p])) if p == "phi" else np.log(s[p])
        return x

    def unpack(self, x):
        """Parameter vector -> (F increments, P increments, a per link, tau[band, 3], phi[band])."""
        get = lambda n: x[self.names == n]                 # noqa: E731
        Pinc = np.zeros((self.nb, self.nk)) if self.no_path else np.exp(get("logP")).reshape(self.nb, self.nk)
        Finc = np.exp(get("logF")).reshape(self.nb, self.nk) if self.mode > 0 else np.zeros((self.nb, self.nk))
        a = np.exp(get("loga")) if self.mode > 0 else np.ones(self.nl)
        tau, phi = np.ones((self.nb, 3)), np.ones(self.nb)
        band = lambda n: np.broadcast_to(np.exp(get(n)), (self.nb,))   # noqa: E731
        if self.model == "M2":
            tau[:, 0] = tau[:, 1] = band("tau")
        elif self.model == "M3":
            tau[:, 0], tau[:, 1] = band("tau_wet"), band("tau_dry")
        elif self.model == "M4":
            tau[:, 0], tau[:, 1], tau[:, 2] = band("tau_f_wet"), band("tau_f_dry"), band("tau_s")
            phi = np.broadcast_to(1 / (1 + np.exp(-get("phi"))), (self.nb,)).copy()
        return Finc, Pinc, a, tau, phi

    def run(self, x, d, wle, want_grad=True):
        """(loss + ridge, gradient) for lag d and link-event weights wle. Leaves self.pred and self.le_loss filled."""
        D = self.D
        Finc, Pinc, a, tau, phi = self.unpack(x)
        gF, gP, ga = np.zeros_like(Finc), np.zeros_like(Pinc), np.zeros(self.nl)
        gtau, gphi = np.zeros((self.nb, 3)), np.zeros(self.nb)
        loss = kernel(self.mode, int(d), Finc, Pinc, a, tau, phi, wle, D["le_roff"], D["le_n"], self.boff, D["le_link"], D["le_band"],
                      D["le_L"], D["le_c"], D["rain"], self.y, self.w, self.agg, self.M["huber_delta_db"], self.knots, self.pred, self.le_loss,
                      gF, gP, ga, gtau, gphi, want_grad)
        g = np.zeros(len(x))
        if self.mode > 0:
            la = x[self.names == "loga"]
            loss += self.pen * float(la @ la)
        if not want_grad:
            return loss, g
        if not self.no_path:
            g[self.names == "logP"] = (gP * Pinc).ravel()
        if self.mode > 0:
            g[self.names == "logF"] = (gF * Finc).ravel()
            g[self.names == "loga"] = ga * a + 2 * self.pen * la
        fold = lambda v: v if self.nt > 1 else np.array([v.sum()])     # noqa: E731  shared parameter: sum over the bands
        if self.model == "M2":
            g[self.names == "tau"] = fold((gtau[:, 0] + gtau[:, 1]) * tau[:, 0])
        elif self.model == "M3":
            g[self.names == "tau_wet"], g[self.names == "tau_dry"] = fold(gtau[:, 0] * tau[:, 0]), fold(gtau[:, 1] * tau[:, 1])
        elif self.model == "M4":
            for i, n in enumerate(("tau_f_wet", "tau_f_dry", "tau_s")):
                g[self.names == n] = fold(gtau[:, i] * tau[:, i])
            g[self.names == "phi"] = fold(gphi * phi * (1 - phi))
        return loss, g

    def fit(self, x0, d, wle, maxiter):
        x0 = np.clip(x0, [b[0] for b in self.bounds], [b[1] for b in self.bounds])
        r = minimize(self.run, x0, args=(d, wle), jac=True, method="L-BFGS-B", bounds=self.bounds,
                     options={"maxiter": maxiter, "maxfun": int(maxiter * 1.5), "ftol": 1e-10, "gtol": 1e-6})
        return {"x": r.x, "loss": float(r.fun), "nit": int(r.nit), "nfev": int(r.nfev), "d": int(d)}

    def predict(self, x, d):
        """Predicted A per output bin for every link-event."""
        self.run(x, d, np.ones(len(self.D["le_n"])), want_grad=False)
        return self.pred.copy()

    def params(self, x, d):
        """The readable parameters: d, the time constants, phi, the curves at the knots, the spread of a_link."""
        Finc, Pinc, a, tau, phi = self.unpack(x)
        out = {"d": int(d)}
        for b in range(self.nb if self.nt > 1 else 1):
            tag = f"[{self.D['bands'][b]}]" if self.nt > 1 else ""
            vals = {"tau": tau[b, 0], "tau_wet": tau[b, 0], "tau_dry": tau[b, 1], "tau_f_wet": tau[b, 0], "tau_f_dry": tau[b, 1],
                    "tau_s": tau[b, 2], "phi": phi[b]}
            out.update({p + tag: float(vals[p]) for p in DYN[self.model]})
        for b in range(self.nb):
            for j in range(self.nk):
                out[f"F[{self.D['bands'][b]}]({self.knots[j + 1]:g})"] = float(Finc[b, :j + 1].sum())
                out[f"P[{self.D['bands'][b]}]({self.knots[j + 1]:g})"] = float(Pinc[b, :j + 1].sum())
        if self.mode > 0:
            la = np.log(a)
            out.update({"a_log_sd": float(la.std()), "a_p10": float(np.exp(np.percentile(la, 10))), "a_p90": float(np.exp(np.percentile(la, 90)))})
        return out


def fit_all(P, starts, lags, maxiter, log=None):
    """The fit on a whole dataset (config `lag_candidates`, first part): every start at d = 0, the lag profile from the
    best one, every start again at the best d. Returns the solution per lag (x, loss, per-link-event loss), the best
    lag and the losses of the starts."""
    ones = np.ones(len(P.D["le_n"]))
    first = {s: P.fit(P.start(v), 0, ones, maxiter["start"]) for s, v in starts.items()}
    prof = {0: min(first.values(), key=lambda r: r["loss"])}
    for step in (1, -1):
        for d in range(step, (lags[1] if step > 0 else lags[0]) + step, step):
            prof[d] = P.fit(prof[d - step]["x"], d, ones, maxiter["lag"])
    d_best = min(prof, key=lambda d: prof[d]["loss"])
    again = {s: P.fit(P.start(v), d_best, ones, maxiter["start"]) for s, v in starts.items()} if d_best != 0 else first
    best = min(list(again.values()) + [prof[d_best]], key=lambda r: r["loss"])
    prof[d_best] = best
    for d, r in prof.items():                              # the loss of every link-event at every lag (for the folds)
        P.run(r["x"], d, ones, want_grad=False)
        r["le_loss"] = P.le_loss.copy()
    return {"profile": prof, "d": int(d_best), "x": best["x"], "loss": best["loss"],
            "starts_d0": {s: r["loss"] for s, r in first.items()}, "starts_dbest": {s: r["loss"] for s, r in again.items()},
            "nit": {s: r["nit"] for s, r in again.items()}}


def fit_weighted(P, allfit, wle, n_cand, maxiter):
    """A fit on a re-weighted training set (a leave-one-event-out fold, or a bootstrap draw): the training loss of the
    profile solutions decides which lags are refitted; the lowest training loss after the refit wins. Links without
    any training link-event get a_link = 1."""
    prof = allfit["profile"]
    pen = lambda x: P.pen * float(x[P.names == "loga"] @ x[P.names == "loga"]) if P.mode > 0 else 0.0   # noqa: E731
    train = {d: float(wle @ r["le_loss"]) + pen(r["x"]) for d, r in prof.items()}
    cands = sorted(train, key=train.get)[:n_cand]
    best = min((P.fit(prof[d]["x"], d, wle, maxiter) for d in cands), key=lambda r: r["loss"])
    if P.mode > 0:
        seen = np.bincount(P.D["le_link"], weights=wle, minlength=P.nl) > 0
        la = best["x"][P.names == "loga"]
        la[~seen] = 0.0
        best["x"][P.names == "loga"] = la
    return best


# =================================================================================================== the data (steps S1 and S2)
ARRAYS = ["le_roff", "le_n", "le_boff5", "le_boff1", "le_link", "le_band", "le_L", "le_c", "le_event", "le_last", "le_first",
          "rain", "y5", "w5", "ok5", "y1", "w1", "pre1", "wet5", "after5", "r5", "bin_le"]


def build_dataset(d, DATA, event_split="train", link_split="train", baselines=None):
    """The link-events of the short links (config section `data` of the step). Returns (D, le, drops): D = the arrays
    the models use, le = one row per link-event, drops = how many pairs each rule removed.
    event_split = the events used ("train" for S1 / S2 / S2b; step S3 also reads "val"; step S8 "test").
    link_split = the link split (the TRAIN links everywhere; step S8 also builds the TEST links).
    baselines = the baseline table (default: the frozen baselines.csv of the baseline module; step S8 passes its TEST baselines)."""
    BANDS = list(DATA["bands_ghz"])
    links = d.links()
    links = links[(links.split_link == link_split) & (links.length_m < DATA["max_length_km"] * 1000)]
    subs = d.sublinks().set_index("sublink_id")
    ev = d.events()
    train_ev = ev[ev.split_event == event_split].set_index("event_id")
    if baselines is None:
        baselines = pd.read_csv(PATHS.BASELINE_MODULE / "outputs" / "products" / "baselines.csv")
    base = baselines.query("split_event == @event_split")   # other splits dropped here
    base = base[base.cml_id.isin(links.cml_id)].merge(
        pd.read_csv(PATHS.BASELINE_UNCERTAINTY)[["sublink_id", "offset_std"]].rename(columns={"offset_std": "offset_std_unc"}), on="sublink_id", how="left")
    drops = {f"pairs of the short TRAIN links in {event_split.upper()} events": len(base)}
    for name, rule in (("no baseline", lambda b: ~b.missing), ("no wet bin of the link in the event", lambda b: b.link_has_wet_slot),
                       ("sub-link not good", lambda b: b.status == "good"), ("use_in_waa_main false", lambda b: b.use_in_waa_main),
                       (f"offset_std > {DATA['max_offset_std_db']} dB", lambda b: ~(b.offset_std_unc > DATA["max_offset_std_db"]))):
        keep = rule(base)                                  # the rules are applied one after the other
        drops[name] = int((~keep).sum())
        base = base[keep]
    base["frequency_ghz"] = base.sublink_id.map(subs.frequency_ghz)
    base["band"] = -1
    for i, (lo, hi) in enumerate(DATA["bands_ghz"].values()):
        base.loc[(base.frequency_ghz >= lo) & (base.frequency_ghz < hi), "band"] = i
    drops["frequency outside the two bands"] = int((base.band < 0).sum())
    base = base[base.band >= 0].sort_values(["event_id", "sublink_id"]).reset_index(drop=True)

    st = d.store("link_reference_state_5min").state.load()
    T5, state = pd.DatetimeIndex(st.time.values), st.to_numpy()
    s_col = {int(c): j for j, c in enumerate(st.cml_id.values)}
    ra = d.store("radar_along_links_5min").R_mean.load()
    assert list(ra.cml_id.values) == list(st.cml_id.values) and (pd.DatetimeIndex(ra.time.values) == T5).all()
    R = ra.to_numpy()
    sub_ids = sorted(base.sublink_id.unique())
    s1 = d.store("cml_1min")
    T1 = pd.DatetimeIndex(s1.time.values)
    assert T1[0] == T5[0]
    TL = s1.tl_mean.sel(sublink_id=sub_ids).load().to_numpy()
    STEP = s1.n_tsl_step.sel(sublink_id=sub_ids).load().to_numpy() > 0
    c1 = {s_: j for j, s_ in enumerate(sub_ids)}
    blank = np.zeros(len(T1), bool)                        # minutes inside an event of another split
    for e in ev[ev.split_event != event_split].itertuples():
        blank[max(T1.searchsorted(e.start), 0):T1.searchsorted(e.end)] = True
    lgm = d.table("link_gauge_matches")
    lgm = lgm[(lgm.network == "city") & (lgm.distance_to_path_m <= DATA["gauge_max_dist_m"])].sort_values("distance_to_path_m").drop_duplicates("cml_id").set_index("cml_id")
    gauge = d.store("gauges_city_1min").rate_mm_h.load().to_pandas()
    assert gauge.index[0] == T1[0]
    wet_bins = {c: np.nonzero(state[:, j] == 1)[0] for c, j in s_col.items() if c in set(base.cml_id)}
    nT, nb_before, nb_after = len(T5), DATA["before_min"] // 5, DATA["after_min"] // 5
    ow = DATA["offset_window_min"]
    length_km = links.set_index("cml_id").length_m / 1000
    rows, parts = [], {k: [] for k in ("rain", "rain_gauge", "y5", "ok5", "y1", "pre1", "wet5", "after5", "r5")}
    for name in ("cut at a TSL step", "TSL step before the rain", "too few valid 5 min bins", "too many 5 min bins without radar"):
        drops[name] = 0
    for r in base.itertuples():
        e = train_ev.loc[r.event_id]
        wb, j = wet_bins[r.cml_id], s_col[r.cml_id]
        mine = wb[(wb >= T5.searchsorted(e.start)) & (wb < T5.searchsorted(e.end))]
        first, last = int(mine[0]), int(mine[-1])
        k0, k1 = np.searchsorted(wb, first), np.searchsorted(wb, last)
        prev, nxt = (int(wb[k0 - 1]) if k0 > 0 else -1), (int(wb[k1 + 1]) if k1 + 1 < len(wb) else nT)
        b0, b1 = max(first - nb_before, prev + 1, 2), min(last + 1 + nb_after, nxt, nT - 2)
        step = np.nonzero(STEP[b0 * 5:b1 * 5, c1[r.sublink_id]])[0]
        if len(step):
            b1 = b0 + int(step[0]) // 5
            if b1 <= first:
                drops["TSL step before the rain"] += 1
                continue
            drops["cut at a TSL step"] += 1
            last = int(mine[mine < b1][-1])
        a1 = TL[b0 * 5:b1 * 5, c1[r.sublink_id]].astype(float) - r.baseline_db
        a1[blank[b0 * 5:b1 * 5]] = np.nan
        m5 = a1.reshape(-1, 5)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            y5 = np.where(np.isfinite(m5).sum(axis=1) >= 3, np.nanmean(m5, axis=1), np.nan)
        r5 = R[b0:b1, j].astype(float)
        if np.isfinite(y5).mean() < DATA["min_valid_share"]:
            drops["too few valid 5 min bins"] += 1
            continue
        if np.isnan(r5).mean() > DATA["max_radar_nan_share"]:
            drops["too many 5 min bins without radar"] += 1
            continue
        minute = np.arange(b0 * 5, b1 * 5)
        pre = (minute >= first * 5 - ow[0]) & (minute < first * 5 - ow[1]) & np.repeat(state[b0:b1, j] == 0, 5)
        n_pre = int((pre & np.isfinite(a1)).sum())
        c = float(np.nanmedian(a1[pre])) if n_pre >= DATA["offset_min_minutes"] else 0.0
        g_id = lgm.gauge_id.get(r.cml_id, None)
        g1 = gauge[g_id].to_numpy()[b0 * 5 - PADM:b1 * 5 + PADM].astype(float) if g_id else np.full(len(a1) + 2 * PADM, np.nan)
        parts["rain"].append(np.repeat(np.nan_to_num(R[b0 - 2:b1 + 2, j].astype(float)), 5))
        parts["rain_gauge"].append(g1)
        parts["y5"].append(y5), parts["y1"].append(a1), parts["pre1"].append(pre), parts["r5"].append(r5)
        parts["ok5"].append(np.isfinite(r5) & ~blank[b0 * 5:b1 * 5].reshape(-1, 5).any(axis=1))
        parts["wet5"].append(state[b0:b1, j] == 1)
        parts["after5"].append((np.arange(b0, b1) - last - 1) * 5)
        rows.append({"event_id": r.event_id, "sublink_id": r.sublink_id, "cml_id": r.cml_id, "L_km": length_km[r.cml_id],
                     "frequency_ghz": r.frequency_ghz, "polarization": subs.polarization[r.sublink_id], "band": int(r.band), "band_name": BANDS[int(r.band)],
                     "baseline_db": r.baseline_db, "offset_std": r.offset_std_unc, "significant": bool(r.significant),
                     "window_start": T5[b0], "window_end": T5[b1 - 1] + pd.Timedelta("5min"), "first_wet": T5[first], "last_wet": T5[last],
                     "n_bins": b1 - b0, "first_bin": first - b0, "last_bin": last - b0, "tsl_cut": bool(len(step)),
                     "offset_c_db": c, "offset_minutes": n_pre, "offset_flag": n_pre < DATA["offset_min_minutes"],
                     "gauge_id": g_id or "", "gauge_distance_m": float(lgm.distance_to_path_m.get(r.cml_id, np.nan))})
    le = pd.DataFrame(rows)
    drops["link-events kept"] = len(le)
    return pack(le, parts, BANDS), le, drops


def pack(le, parts, bands):
    """The arrays of a list of link-events as one dataset (everything concatenated, with the offsets per link-event)."""
    n5 = le.n_bins.to_numpy()
    D = {k: np.concatenate(v) for k, v in parts.items()}
    D.update({"le_n": n5 * 5, "le_boff5": np.r_[0, np.cumsum(n5)[:-1]], "le_boff1": np.r_[0, np.cumsum(n5 * 5)[:-1]],
              "le_roff": np.r_[0, np.cumsum(n5 * 5 + 2 * PADM)[:-1]], "le_link": pd.factorize(le.cml_id)[0], "le_band": le.band.to_numpy(),
              "le_L": le.L_km.to_numpy(float), "le_c": le.offset_c_db.to_numpy(float), "le_event": pd.factorize(le.event_id)[0],
              "le_first": le.first_bin.to_numpy(), "le_last": le.last_bin.to_numpy(), "bin_le": np.repeat(np.arange(len(le)), n5)})
    return finish(D, bands)


def finish(D, bands):
    """The masks that follow from y and ok5, and the types numba wants."""
    D["w5"] = (D["ok5"] & np.isfinite(D["y5"])).astype(float)
    D["w1"] = (np.repeat(D["ok5"], 5) & np.isfinite(D["y1"])).astype(float)
    D["y5"], D["y1"] = np.nan_to_num(D["y5"]), np.nan_to_num(D["y1"])      # the weights say which values count
    for k in ("le_roff", "le_n", "le_boff5", "le_boff1", "le_link", "le_band", "le_event", "le_first", "le_last", "bin_le", "after5"):
        D[k] = D[k].astype(np.int64)
    D["nb"], D["nl"], D["bands"] = len(bands), int(D["le_link"].max()) + 1, list(bands)
    return D


def subset(D, keep, renumber_events=True):
    """The dataset with only the link-events `keep` (boolean per link-event). Links (and events) are renumbered."""
    idx = np.nonzero(keep)[0]
    out = {}
    for k in ("y5", "ok5", "wet5", "after5", "r5"):
        out[k] = np.concatenate([D[k][D["le_boff5"][i]:D["le_boff5"][i] + D["le_n"][i] // 5] for i in idx])
    for k in ("y1", "pre1"):
        out[k] = np.concatenate([D[k][D["le_boff1"][i]:D["le_boff1"][i] + D["le_n"][i]] for i in idx])
    for k in ("rain", "rain_gauge"):
        out[k] = np.concatenate([D[k][D["le_roff"][i]:D["le_roff"][i] + D["le_n"][i] + 2 * PADM] for i in idx])
    out["y5"] = np.where(np.concatenate([D["w5"][D["le_boff5"][i]:D["le_boff5"][i] + D["le_n"][i] // 5] for i in idx]) > 0, out["y5"], np.nan)
    out["y1"] = np.where(np.concatenate([D["w1"][D["le_boff1"][i]:D["le_boff1"][i] + D["le_n"][i]] for i in idx]) > 0, out["y1"], np.nan)
    n5 = D["le_n"][idx] // 5
    out.update({"le_n": n5 * 5, "le_boff5": np.r_[0, np.cumsum(n5)[:-1]], "le_boff1": np.r_[0, np.cumsum(n5 * 5)[:-1]],
                "le_roff": np.r_[0, np.cumsum(n5 * 5 + 2 * PADM)[:-1]], "le_link": np.unique(D["le_link"][idx], return_inverse=True)[1],
                "le_band": D["le_band"][idx], "le_L": D["le_L"][idx], "le_c": D["le_c"][idx],
                "le_event": np.unique(D["le_event"][idx], return_inverse=True)[1] if renumber_events else D["le_event"][idx],
                "le_first": D["le_first"][idx], "le_last": D["le_last"][idx], "bin_le": np.repeat(np.arange(len(idx)), n5)})
    return finish(out, D["bands"])


def save_data(D, file):
    np.savez_compressed(file, **{k: D[k] for k in ARRAYS + ["rain_gauge"] if k in D})
    return str(file)


def load_data(file, bands):
    D = dict(np.load(file))
    D["nb"], D["nl"], D["bands"] = len(bands), int(D["le_link"].max()) + 1, list(bands)
    return D


# =================================================================================================== the judge (steps S1 and S2)
def event_errors(D, pred, J, w=None):
    """Per held-out event: the RMSE of A over its 5 min points (all, during rain, 0-60 and 60-180 min after the last
    wet bin), and the means over its link-events of the peak error and of the integral of A - c after the rain.
    J = config section `judge` of the step."""
    w = (D["w5"] > 0) if w is None else w
    ev_bin, n_ev = D["le_event"][D["bin_le"]], int(D["le_event"].max()) + 1
    res2 = np.where(w, (pred - D["y5"]) ** 2, 0.0)
    out = {}
    (a0, a1), (b0, b1) = J["after_windows_min"]
    masks = {"rmse": w, "rmse_rain": w & D["wet5"], f"rmse_{a0}_{a1}": w & (D["after5"] >= a0) & (D["after5"] < a1),
             f"rmse_{b0}_{b1}": w & (D["after5"] >= b0) & (D["after5"] < b1)}
    with np.errstate(all="ignore"):
        for name, m in masks.items():
            out[name] = np.sqrt(np.bincount(ev_bin[m], weights=res2[m], minlength=n_ev) / np.bincount(ev_bin[m], minlength=n_ev))
        le_tab = link_event_table(D, pred, J, w)
        g = le_tab.groupby(D["le_event"])
        out["peak_error"] = g.peak_error.mean().reindex(range(n_ev)).to_numpy()
        out["integral_obs"], out["integral_pred"] = g.integral_obs.mean().reindex(range(n_ev)).to_numpy(), g.integral_pred.mean().reindex(range(n_ev)).to_numpy()
    return pd.DataFrame(out)


def link_event_table(D, pred, J, w=None):
    """Per link-event: observed and predicted peak of A, the integral of A - c over 0-180 min after the last wet bin
    (dB min), and the residual statistics."""
    w = (D["w5"] > 0) if w is None else w
    after = w & (D["after5"] >= 0) & (D["after5"] < J["after_windows_min"][-1][1])
    c_bin = D["le_c"][D["bin_le"]]
    df = pd.DataFrame({"lev": D["bin_le"], "obs": np.where(w, D["y5"], np.nan), "pred": np.where(w, pred, np.nan), "res": np.where(w, pred - D["y5"], np.nan),
                       "io": np.where(after, (D["y5"] - c_bin) * 5, 0.0), "ip": np.where(after, (pred - c_bin) * 5, 0.0), "wet": w & D["wet5"]})
    g = df.groupby("lev")
    out = pd.DataFrame({"points": g.res.count(), "wet_points": g.wet.sum(), "peak_obs": g.obs.max(), "peak_pred": g.pred.max(),
                        "integral_obs": g.io.sum(), "integral_pred": g.ip.sum(), "bias": g.res.mean(), "rmse": np.sqrt((df.res ** 2).groupby(df.lev).mean())})
    out["peak_error"] = out.peak_pred - out.peak_obs
    return out.reindex(range(len(D["le_n"])))


def judge(errors, order, J):
    """The ladder table with event-bootstrap intervals and the stop rule. errors = {model: per-event table}."""
    n_ev = len(next(iter(errors.values())))
    rng = np.random.default_rng(J["seed"])
    draws = rng.integers(0, n_ev, (J["boot_reps"], n_ev))
    lo_hi = lambda v: np.nanpercentile(v, J["interval"])                 # noqa: E731
    mean = {m: e.mean() for m, e in errors.items()}
    boot = {m: {c: np.nanmean(e[c].to_numpy()[draws], axis=1) for c in e} for m, e in errors.items()}
    best = min(order, key=lambda m: mean[m]["rmse"])
    rows, chosen = [], None
    for m in order:
        dif = boot[m]["rmse"] - boot[best]["rmse"]
        d_lo, d_hi = lo_hi(dif)
        inside = bool(d_lo <= 0 <= d_hi) or m == best
        if chosen is None and inside:
            chosen = m
        skill = 1 - boot[m]["rmse"] / boot["M0"]["rmse"] if "M0" in boot else np.full(J["boot_reps"], np.nan)
        row = {"model": m, "rmse": mean[m]["rmse"], "rmse_lo": lo_hi(boot[m]["rmse"])[0], "rmse_hi": lo_hi(boot[m]["rmse"])[1],
               "minus_best": mean[m]["rmse"] - mean[best]["rmse"], "diff_lo": d_lo, "diff_hi": d_hi, "interval_holds_0": inside,
               "skill_vs_M0": 1 - mean[m]["rmse"] / mean["M0"]["rmse"] if "M0" in mean else np.nan, "skill_lo": lo_hi(skill)[0], "skill_hi": lo_hi(skill)[1]}
        for c in errors[m].columns[1:]:
            row[c] = mean[m][c]
            if c.startswith("rmse") or c.startswith("peak_error"):
                row[c + "_lo"], row[c + "_hi"] = lo_hi(boot[m][c])
        row["integral_pred_over_obs"] = mean[m]["integral_pred"] / mean[m]["integral_obs"]
        rows.append(row)
    return pd.DataFrame(rows), best, chosen


def showers(D, e, gap_min):
    """The showers of a link-event: runs of wet bins, gaps up to gap_min bridged. Returns [(first bin, last bin)]."""
    o, n = D["le_boff5"][e], D["le_n"][e] // 5
    wet = np.nonzero(D["wet5"][o:o + n])[0]
    if not len(wet):
        return []
    cut = np.nonzero(np.diff(wet) > gap_min // 5 + 1)[0]
    return list(zip(wet[np.r_[0, cut + 1]], wet[np.r_[cut, len(wet) - 1]]))


def shower_table(D, gap_min):
    """Every shower of every link-event: (link-event, event, first bin, last bin as positions in the 5 min arrays)."""
    rows = []
    for e in range(len(D["le_n"])):
        o = D["le_boff5"][e]
        rows += [(e, D["le_event"][e], o + s0, o + s1 + 1) for s0, s1 in showers(D, e, gap_min)]
    return np.array(rows, dtype=np.int64).reshape(-1, 4)


def shower_peaks(D, sh, pred, min_peak):
    """Predicted - observed peak of A for every shower whose observed peak of A - c is at least min_peak.
    Returns (event of the shower, peak error)."""
    w = D["w5"] > 0
    obs, prd = np.where(w, D["y5"], -np.inf), np.where(w, pred, -np.inf)
    ev, err = [], []
    for e, event, i0, i1 in sh:
        po = obs[i0:i1].max()
        if np.isfinite(po) and po - D["le_c"][e] >= min_peak:
            ev.append(event), err.append(prd[i0:i1].max() - po)
    return np.array(ev, dtype=np.int64), np.array(err)


# =================================================================================================== the models of step S2
# M1, M3, M4 as in step S1; plus a level drift per link-event (M3D) and the discovery rung M5:
#     dW/dt = s_band * a_link * g(R, W) - h(W)
# g = values on a grid of rain knots x W knots (bilinear, 0 at R = 0, never falling with R); h = piecewise linear in W,
# 0 at W = 0, never falling. M5 is integrated with explicit Euler steps inside each minute, W kept >= 0.
SPEC2 = {"M1": (1, [], False), "M3": (2, ["tau_wet", "tau_dry"], False), "M3D": (2, ["tau_wet", "tau_dry"], True),
         "M4": (2, ["phi", "tau_f_wet", "tau_f_dry", "tau_s"], False), "M5": (3, [], False), "M5D": (3, [], True),   # mode, dynamic parameters, drift
         "M4b": (2, ["phi", "tau_f_wet", "tau_f_dry", "tau_s"], False),   # S2b: M4 with tau_s within bounds.tau_s_min
         "M4c": (2, ["phi", "tau_f_wet", "tau_f_dry", "tau_s_wet", "tau_s_dry"], False)}   # S5a: the slow store fills with
                                                    #   tau_s_wet and empties with tau_s_dry (bounds.tau_s_split)


@njit(cache=True)
def kernel2(mode, d, Finc, Pinc, a, tau, phi, G, H, rk, wk, sb, nsub, bd, wle, le_roff, le_n, le_boff, le_link, le_band, le_L, le_c, le_t0,
            rain, y, w, agg, delta, knots, pred, le_loss, gF, gP, ga, gtau, gphi, gG, gH, gsb, gb, want_grad, split_s):
    """As `kernel`, for the models of step S2. mode 1 static, 2 stores (M3 = one store with phi = 1), 3 the free
    surfaces G[R knot, W knot] and H[W knot]. bd = drift slope per link-event (dB/h, 0 if the model has none), counted
    from the minute le_t0 of the window. Gradients w.r.t. the VALUES of G and H are returned in gG and gH.
    split_s (M4c, step S5a): the slow store fills with tau[b, 2] and empties with tau[b, 3]; otherwise it has the one
    time constant tau[b, 2] (and tau[b, 3] is never read)."""
    nb, nk = Pinc.shape
    nR, nW = G.shape
    cF, cP = np.zeros((nb, nk)), np.zeros((nb, nk))
    for b in range(nb):
        for j in range(1, nk):
            cF[b, j] = cF[b, j - 1] + Finc[b, j - 1]
            cP[b, j] = cP[b, j - 1] + Pinc[b, j - 1]
    maxn = le_n.max()
    Wf, Ws, fr = np.empty(maxn), np.empty(maxn), np.empty(maxn)
    jj, br, bs = np.empty(maxn, np.int64), np.empty(maxn, np.int64), np.empty(maxn, np.int64)
    psi = np.empty(maxn)
    jr, f2s = np.empty(maxn, np.int64), np.empty(maxn)
    Wu, qu = np.empty(maxn * nsub), np.empty(maxn * nsub)
    ku, cu = np.empty(maxn * nsub, np.int64), np.empty(maxn * nsub, np.int64)
    G0F, G1F, G0P, G1P = np.zeros((nb, nk)), np.zeros((nb, nk)), np.zeros((nb, nk)), np.zeros((nb, nk))
    dt = 1.0 / nsub
    total = 0.0
    for e in range(len(le_n)):
        n, b, link, L, c = le_n[e], le_band[e], le_link[e], le_L[e], le_c[e]
        al, ph, slope, t0 = a[link], phi[b], bd[e], le_t0[e]
        scale = sb[b] * al
        ro, bo = le_roff[e] + PADM - d, le_boff[e]
        kfw, kfd, ks = np.exp(-1.0 / tau[b, 0]), np.exp(-1.0 / tau[b, 1]), np.exp(-1.0 / tau[b, 2])
        ksd = np.exp(-1.0 / tau[b, 3]) if split_s else ks
        wf, ws, wu, acc = 0.0, 0.0, 0.0, 0.0
        for t in range(n):                                 # forward: simulate the minutes
            r = rain[ro + t]
            j, f_ = nk - 1, 1.0
            if r <= 0.0:
                j, f_ = 0, 0.0
            else:
                for q in range(nk):
                    if r < knots[q + 1]:
                        j, f_ = q, (r - knots[q]) / (knots[q + 1] - knots[q])
                        break
            jj[t], fr[t] = j, f_
            wtot = 0.0
            if mode == 1:
                wtot = al * (cF[b, j] + f_ * Finc[b, j])
            elif mode == 2:
                u = al * (cF[b, j] + f_ * Finc[b, j])
                uf, us = ph * u, (1.0 - ph) * u
                Wf[t], Ws[t] = wf, ws
                if uf > wf:
                    br[t] = 0
                    wf = uf + (wf - uf) * kfw
                else:
                    br[t] = 1
                    wf = uf + (wf - uf) * kfd
                if split_s and us <= ws:                   # M4c: the slow store empties with its own time constant
                    bs[t] = 1
                    ws = us + (ws - us) * ksd
                else:
                    bs[t] = 0
                    ws = us + (ws - us) * ks
                wtot = wf + ws
            elif mode == 3:
                j2, f2 = nR - 2, 1.0                       # the rain cell of the surface (constant above the last knot)
                if r <= 0.0:
                    j2, f2 = 0, 0.0
                else:
                    for q in range(nR - 1):
                        if r < rk[q + 1]:
                            j2, f2 = q, (r - rk[q]) / (rk[q + 1] - rk[q])
                            break
                jr[t], f2s[t] = j2, f2
                for s in range(nsub):
                    i = t * nsub + s
                    k2 = nW - 2
                    for q in range(nW - 1):
                        if wu < wk[q + 1]:
                            k2 = q
                            break
                    qh = (wu - wk[k2]) / (wk[k2 + 1] - wk[k2])     # above the last knot: h goes on with its last slope,
                    qg = qh if qh < 1.0 else 1.0                   #   g stays constant
                    g0 = G[j2, k2] + qg * (G[j2, k2 + 1] - G[j2, k2])
                    g1 = G[j2 + 1, k2] + qg * (G[j2 + 1, k2 + 1] - G[j2 + 1, k2])
                    hv = H[k2] + qh * (H[k2 + 1] - H[k2])
                    Wu[i], ku[i], qu[i] = wu, k2, qh
                    wu = wu + dt * (scale * (g0 + f2 * (g1 - g0)) - hv)
                    cu[i] = 0
                    if wu < 0.0:
                        wu, cu[i] = 0.0, 1
                wtot = wu
            acc += wtot + L * (cP[b, j] + f_ * Pinc[b, j]) + c + slope * (t - t0) / 60.0
            if (t + 1) % agg == 0:
                pred[bo + t // agg] = acc / agg
                acc = 0.0
        loss = 0.0
        for q in range(n // agg):                          # Huber loss of the output bins
            psi[q] = 0.0
            if w[bo + q] > 0:
                res = pred[bo + q] - y[bo + q]
                if abs(res) <= delta:
                    loss += 0.5 * res * res
                    psi[q] = res
                else:
                    loss += delta * (abs(res) - 0.5 * delta)
                    psi[q] = delta if res > 0 else -delta
        le_loss[e] = loss
        total += wle[e] * loss
        if want_grad and wle[e] != 0.0:                    # backward: the gradient
            lf, ls, gkfw, gkfd, gks, gksd, lam = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            for t in range(n - 1, -1, -1):
                gp = wle[e] * psi[t // agg] / agg
                gb[e] += gp * (t - t0) / 60.0
                j, f_ = jj[t], fr[t]
                G0P[b, j] += gp * L
                G1P[b, j] += gp * L * f_
                if mode == 3:
                    lam += gp
                    j2, f2 = jr[t], f2s[t]
                    for s in range(nsub - 1, -1, -1):
                        i = t * nsub + s
                        if cu[i] == 1:                     # W was set to 0 here: nothing before it matters
                            lam = 0.0
                            continue
                        k2, qh = ku[i], qu[i]
                        qg = qh if qh < 1.0 else 1.0
                        dW = wk[k2 + 1] - wk[k2]
                        g0 = G[j2, k2] + qg * (G[j2, k2 + 1] - G[j2, k2])
                        g1 = G[j2 + 1, k2] + qg * (G[j2 + 1, k2 + 1] - G[j2 + 1, k2])
                        gv = g0 + f2 * (g1 - g0)
                        co = dt * lam
                        ga[link] += co * sb[b] * gv
                        gsb[b] += co * al * gv
                        cs = co * scale
                        gG[j2, k2] += cs * (1.0 - f2) * (1.0 - qg)
                        gG[j2, k2 + 1] += cs * (1.0 - f2) * qg
                        gG[j2 + 1, k2] += cs * f2 * (1.0 - qg)
                        gG[j2 + 1, k2 + 1] += cs * f2 * qg
                        gH[k2] -= co * (1.0 - qh)
                        gH[k2 + 1] -= co * qh
                        dg = 0.0
                        if qh < 1.0:
                            dg = ((1.0 - f2) * (G[j2, k2 + 1] - G[j2, k2]) + f2 * (G[j2 + 1, k2 + 1] - G[j2 + 1, k2])) / dW
                        lam = lam * (1.0 + dt * (scale * dg - (H[k2 + 1] - H[k2]) / dW))
                    continue
                fval = cF[b, j] + f_ * Finc[b, j]
                gu = gp
                if mode == 2:
                    u = al * fval
                    uf, us = ph * u, (1.0 - ph) * u
                    lf += gp
                    ls += gp
                    k = kfw if br[t] == 0 else kfd
                    guf = lf * (1.0 - k)
                    if br[t] == 0:
                        gkfw += lf * (Wf[t] - uf)
                    else:
                        gkfd += lf * (Wf[t] - uf)
                    lf *= k
                    if split_s and bs[t] == 1:
                        gus = ls * (1.0 - ksd)
                        gksd += ls * (Ws[t] - us)
                        ls *= ksd
                    else:
                        gus = ls * (1.0 - ks)
                        gks += ls * (Ws[t] - us)
                        ls *= ks
                    gu = ph * guf + (1.0 - ph) * gus
                    gphi[b] += u * (guf - gus)
                ga[link] += gu * fval
                G0F[b, j] += gu * al
                G1F[b, j] += gu * al * f_
            gtau[b, 0] += gkfw * kfw / tau[b, 0] ** 2
            gtau[b, 1] += gkfd * kfd / tau[b, 1] ** 2
            gtau[b, 2] += gks * ks / tau[b, 2] ** 2
            if split_s:
                gtau[b, 3] += gksd * ksd / tau[b, 3] ** 2
    if want_grad:                                          # a segment's increment also lifts every segment above it
        for b in range(nb):
            sF, sP = 0.0, 0.0
            for j in range(nk - 1, -1, -1):
                gF[b, j] = sF + G1F[b, j]
                gP[b, j] = sP + G1P[b, j]
                sF += G0F[b, j]
                sP += G0P[b, j]
    return total


@njit(cache=True)
def surface_path(G, H, rk, wk, scale, r, w0, minutes, nsub):
    """W per minute under a constant rain r from W = w0, for the surfaces of M5 (the same Euler steps as `kernel2`)."""
    nR, nW = G.shape
    j2, f2 = nR - 2, 1.0
    if r <= 0.0:
        j2, f2 = 0, 0.0
    else:
        for q in range(nR - 1):
            if r < rk[q + 1]:
                j2, f2 = q, (r - rk[q]) / (rk[q + 1] - rk[q])
                break
    out, wu, dt = np.empty(minutes), w0, 1.0 / nsub
    for t in range(minutes):
        for s in range(nsub):
            k2 = nW - 2
            for q in range(nW - 1):
                if wu < wk[q + 1]:
                    k2 = q
                    break
            qh = (wu - wk[k2]) / (wk[k2 + 1] - wk[k2])
            qg = qh if qh < 1.0 else 1.0
            g0 = G[j2, k2] + qg * (G[j2, k2 + 1] - G[j2, k2])
            g1 = G[j2 + 1, k2] + qg * (G[j2 + 1, k2 + 1] - G[j2 + 1, k2])
            wu = max(wu + dt * (scale * (g0 + f2 * (g1 - g0)) - (H[k2] + qh * (H[k2 + 1] - H[k2]))), 0.0)
        out[t] = wu
    return out


def crossing_time(path, level):
    """Minutes until a path (one value per minute, from minute 0) first crosses `level`, interpolated inside the minute."""
    side = path >= level
    i = int(np.argmax(side != side[0]))
    if i == 0:
        return np.nan
    return float(i - 1 + (level - path[i - 1]) / (path[i] - path[i - 1]))


class Problem2:
    """One model of step S2 on one dataset (the same interface as `Problem`). M = config `system_id_s2.model`.
    lam = the smoothness weight of M5; s_b = the prior sd of the drift slopes (dB/h) of a model with drift.
    M["increments"] = "linear" (step S2b): the increments of P and F are plain parameters >= 0 (bounds.increment)
    instead of exp(parameter), so that an increment cannot get stuck near 0 (the trap found in S2)."""

    def __init__(self, D, model, M, agg=5, lam=0.0, s_b=None):
        self.D, self.model, self.M, self.agg = D, model, M, agg
        self.mode, self.dyn, self.drift = SPEC2[model]
        self.lin = M.get("increments") == "linear"
        self.pname, self.fname = ("Pinc", "Finc") if self.lin else ("logP", "logF")
        self.nb, self.nl, self.ne, self.nk = int(D["nb"]), int(D["nl"]), len(D["le_n"]), len(M["knots_mmh"]) - 1
        self.knots = np.asarray(M["knots_mmh"], float)
        self.rk, self.wk = np.asarray(M["m5"]["r_knots_mmh"], float), np.asarray(M["m5"]["w_knots_db"], float)
        self.nR, self.nW, self.nsub = len(self.rk), len(self.wk), int(M["m5"]["substeps_per_min"])
        self.y, self.w, self.boff = (D["y5"], D["w5"], D["le_boff5"]) if agg == 5 else (D["y1"], D["w1"], D["le_boff1"])
        self.pen = 0.5 * (M["huber_delta_db"] / M["a_prior_sd"]) ** 2
        self.s_b = M["drift"]["s_b_start_db_h"] if s_b is None else s_b
        self.lam = lam * float((D["w5"] > 0).sum()) / 1e5
        self.t0 = (D["le_first"] * 5 - M["drift"]["t_ref_min_before_rain"]).astype(float)
        self.tau_h = None                                  # centred time of every 5 min point in hours (made when needed)
        B = M["bounds"]
        names, lo, hi = [], [], []

        def add(name, n, bound):
            names.extend([name] * n), lo.extend(np.broadcast_to(bound[0], n)), hi.extend(np.broadcast_to(bound[1], n))

        inc_bound = B["increment"] if self.lin else B["log_increment"]
        add(self.pname, self.nb * self.nk, inc_bound)
        if self.mode in (1, 2):
            add(self.fname, self.nb * self.nk, inc_bound)
        if self.mode == 3:
            add("logG", (self.nR - 1) * self.nW, B["log_g_increment"])
            add("logH", self.nW - 1, (B["log_g_increment"][0], np.log(np.diff(self.wk) * B["h_max_slope_per_min"])))
            add("logs", self.nb - 1, B["log_band_scale"])
        add("loga", self.nl, B["log_a"])
        for p in self.dyn:
            tb = B["tau_s_min"] if (p == "tau_s" and model == "M4b") else (B["tau_s_split"] if p in ("tau_s_wet", "tau_s_dry") else B["tau_min"])
            add(p, 1, B["logit_phi"] if p == "phi" else list(np.log(tb)))
        if self.drift:
            add("b", self.ne, B["drift_db_h"])
        self.names, self.bounds = np.array(names), list(zip(lo, hi))
        self.pred, self.le_loss = np.zeros(len(self.y)), np.zeros(self.ne)

    def start(self, s, Pinc=None):
        """A generic start (config `starts`): F and P straight lines, a = 1, the start's time constants; for M5 the
        surface g = (f_top min(R, R_last) / R_last / tau_wet) (1 - 0.5 W / W_last) and h = W / tau_dry.
        Pinc (nb x nk, linear increments only) replaces the straight line of P (step S2b: P starts)."""
        x = np.zeros(len(self.names))
        width = np.diff(self.knots) / self.knots[-1]
        tr = (lambda v: v) if self.lin else np.log                                  # noqa: E731
        x[self.names == self.pname] = tr(np.tile(s["p_top"] * width, self.nb) if Pinc is None else np.asarray(Pinc, float).ravel())
        if self.mode in (1, 2):
            x[self.names == self.fname] = tr(np.tile(s["f_top"] * width, self.nb))
        if self.mode == 3:
            inc = np.outer(s["f_top"] * np.diff(self.rk) / self.rk[-1] / s["tau_wet"], 1 - 0.5 * self.wk / self.wk[-1])
            x[self.names == "logG"] = np.log(inc).ravel()
            x[self.names == "logH"] = np.log(np.diff(self.wk) / s["tau_dry"])
        for p in self.dyn:
            x[self.names == p] = np.log(s[p] / (1 - s[p])) if p == "phi" else np.log(s[p])
        return x

    def unpack(self, x):
        get = lambda n: x[self.names == n]                 # noqa: E731
        tr = (lambda v: v) if self.lin else np.exp                                  # noqa: E731
        Pinc = tr(get(self.pname)).reshape(self.nb, self.nk)
        Finc = tr(get(self.fname)).reshape(self.nb, self.nk) if self.mode in (1, 2) else np.zeros((self.nb, self.nk))
        a = np.exp(get("loga"))
        tau, phi = np.ones((self.nb, 4)), np.ones(self.nb)
        if "tau_wet" in self.dyn:
            tau[:, 0], tau[:, 1] = np.exp(get("tau_wet"))[0], np.exp(get("tau_dry"))[0]
        if "tau_s" in self.dyn:
            tau[:, 0], tau[:, 1], tau[:, 2] = np.exp(get("tau_f_wet"))[0], np.exp(get("tau_f_dry"))[0], np.exp(get("tau_s"))[0]
            phi[:] = 1 / (1 + np.exp(-get("phi")[0]))
        if "tau_s_wet" in self.dyn:
            tau[:, 0], tau[:, 1] = np.exp(get("tau_f_wet"))[0], np.exp(get("tau_f_dry"))[0]
            tau[:, 2], tau[:, 3] = np.exp(get("tau_s_wet"))[0], np.exp(get("tau_s_dry"))[0]
            phi[:] = 1 / (1 + np.exp(-get("phi")[0]))
        G, H, sb = np.zeros((self.nR, self.nW)), np.zeros(self.nW), np.ones(self.nb)
        if self.mode == 3:
            G[1:] = np.cumsum(np.exp(get("logG")).reshape(self.nR - 1, self.nW), axis=0)
            H[1:] = np.cumsum(np.exp(get("logH")))
            sb[1:] = np.exp(get("logs"))
        b = get("b") if self.drift else np.zeros(self.ne)
        return Finc, Pinc, a, tau, phi, G, H, sb, b

    def rough(self, G, H):
        """Smoothness penalty of M5 and its gradient w.r.t. the values of G and H: squared second differences on the knot index."""
        gG, gH, val = np.zeros_like(G), np.zeros_like(H), 0.0
        for axis in (0, 1):
            d2 = np.diff(G, 2, axis=axis)
            val += float((d2 ** 2).sum())
            sl = lambda i, j: tuple(slice(i, j) if ax == axis else slice(None) for ax in (0, 1))   # noqa: E731
            n = G.shape[axis]
            gG[sl(0, n - 2)] += 2 * d2
            gG[sl(1, n - 1)] -= 4 * d2
            gG[sl(2, n)] += 2 * d2
        d2 = np.diff(H, 2)
        val += float((d2 ** 2).sum())
        gH[:-2] += 2 * d2
        gH[1:-1] -= 4 * d2
        gH[2:] += 2 * d2
        return val, gG, gH

    def run(self, x, d, wle, want_grad=True):
        """(loss + penalties, gradient) for lag d and link-event weights wle. Leaves self.pred and self.le_loss filled."""
        D = self.D
        Finc, Pinc, a, tau, phi, G, H, sb, b = self.unpack(x)
        gF, gP, ga = np.zeros_like(Finc), np.zeros_like(Pinc), np.zeros(self.nl)
        gtau, gphi, gG, gH, gsb, gb = np.zeros((self.nb, 4)), np.zeros(self.nb), np.zeros_like(G), np.zeros_like(H), np.zeros(self.nb), np.zeros(self.ne)
        loss = kernel2(self.mode, int(d), Finc, Pinc, a, tau, phi, G, H, self.rk, self.wk, sb, self.nsub, b, wle, D["le_roff"], D["le_n"], self.boff,
                       D["le_link"], D["le_band"], D["le_L"], D["le_c"], self.t0, D["rain"], self.y, self.w, self.agg, self.M["huber_delta_db"],
                       self.knots, self.pred, self.le_loss, gF, gP, ga, gtau, gphi, gG, gH, gsb, gb, want_grad, "tau_s_wet" in self.dyn)
        g = np.zeros(len(x))
        la = x[self.names == "loga"]
        loss += self.pen * float(la @ la)
        pen_b = 0.5 * (self.M["huber_delta_db"] / self.s_b) ** 2
        if self.drift:
            loss += pen_b * float(wle @ b ** 2)
        if self.mode == 3 and self.lam > 0:
            val, rG, rH = self.rough(G, H)
            loss += self.lam * val
            gG, gH = gG + self.lam * rG, gH + self.lam * rH
        if not want_grad:
            return loss, g
        g[self.names == self.pname] = (gP if self.lin else gP * Pinc).ravel()
        g[self.names == "loga"] = ga * a + 2 * self.pen * la
        if self.mode in (1, 2):
            g[self.names == self.fname] = (gF if self.lin else gF * Finc).ravel()
        if self.mode == 3:                                 # an increment in R lifts every knot above it
            up = np.cumsum(gG[:0:-1], axis=0)[::-1]
            g[self.names == "logG"] = (up * np.diff(G, axis=0)).ravel()
            g[self.names == "logH"] = np.cumsum(gH[:0:-1])[::-1] * np.diff(H)
            g[self.names == "logs"] = gsb[1:] * sb[1:]
        if "tau_wet" in self.dyn:
            g[self.names == "tau_wet"], g[self.names == "tau_dry"] = (gtau[:, 0] * tau[:, 0]).sum(), (gtau[:, 1] * tau[:, 1]).sum()
        if "tau_s" in self.dyn:
            for i, n in enumerate(("tau_f_wet", "tau_f_dry", "tau_s")):
                g[self.names == n] = (gtau[:, i] * tau[:, i]).sum()
            g[self.names == "phi"] = (gphi * phi * (1 - phi)).sum()
        if "tau_s_wet" in self.dyn:
            for i, n in enumerate(("tau_f_wet", "tau_f_dry", "tau_s_wet", "tau_s_dry")):
                g[self.names == n] = (gtau[:, i] * tau[:, i]).sum()
            g[self.names == "phi"] = (gphi * phi * (1 - phi)).sum()
        if self.drift:
            g[self.names == "b"] = gb + 2 * pen_b * wle * b
        return loss, g

    def fit(self, x0, d, wle, maxiter):
        x0 = np.clip(x0, [b[0] for b in self.bounds], [b[1] for b in self.bounds])
        r = minimize(self.run, x0, args=(d, wle), jac=True, method="L-BFGS-B", bounds=self.bounds,
                     options={"maxiter": maxiter, "maxfun": int(maxiter * 1.5), "ftol": 1e-10, "gtol": 1e-6})
        return {"x": r.x, "loss": float(r.fun), "nit": int(r.nit), "nfev": int(r.nfev), "d": int(d), "converged": bool(r.success)}

    def predict(self, x, d):
        self.run(x, d, np.ones(self.ne), want_grad=False)
        return self.pred.copy()

    def drift_times(self):
        """(t - t_ref) in hours of every 5 min point (its middle)."""
        if self.tau_h is None:
            D = self.D
            pos = np.arange(len(D["y5"])) - D["le_boff5"][D["bin_le"]]
            self.tau_h = (pos * 5 + 2.0 - self.t0[D["bin_le"]]) / 60.0
        return self.tau_h

    def update_s_b(self, x):
        """One EM round for the prior sd of the drift slopes (config `drift`)."""
        D, delta = self.D, self.M["huber_delta_db"]
        b = x[self.names == "b"]
        s_tt = np.bincount(D["bin_le"], weights=D["w5"] * self.drift_times() ** 2, minlength=self.ne)
        self.s_b = float(np.sqrt(np.mean(b ** 2 + 1.0 / (s_tt / delta ** 2 + 1.0 / self.s_b ** 2))))
        return self.s_b

    def params(self, x, d, rain=(0.5, 2.0, 5.0)):
        """The readable parameters: d, the time constants, phi, the curves at the knots, the spread of a_link; for a
        model with drift s_b and the mean / sd of the slopes; for M5 the surfaces and the numbers derived from them
        (t_wet, t_dry, the level reached under constant rain, per band)."""
        Finc, Pinc, a, tau, phi, G, H, sb, b = self.unpack(x)
        present = np.bincount(self.D["le_band"], minlength=self.nb) > 0
        out = {"d": int(d)}
        vals = {"tau_wet": tau[0, 0], "tau_dry": tau[0, 1], "tau_f_wet": tau[0, 0], "tau_f_dry": tau[0, 1], "tau_s": tau[0, 2], "phi": phi[0],
                "tau_s_wet": tau[0, 2], "tau_s_dry": tau[0, 3]}
        out.update({p: float(vals[p]) for p in self.dyn})
        for bi in range(self.nb):
            name = self.D["bands"][bi]
            for j in range(self.nk):
                if self.mode in (1, 2):
                    out[f"F[{name}]({self.knots[j + 1]:g})"] = float(Finc[bi, :j + 1].sum()) if present[bi] else np.nan
                out[f"P[{name}]({self.knots[j + 1]:g})"] = float(Pinc[bi, :j + 1].sum()) if present[bi] else np.nan
        la = np.log(a)
        out.update({"a_log_sd": float(la.std()), "a_p10": float(np.exp(np.percentile(la, 10))), "a_p90": float(np.exp(np.percentile(la, 90)))})
        if self.drift:
            out.update({"s_b": float(self.s_b), "b_mean": float(b.mean()), "b_sd": float(b.std())})
        if self.mode == 3:
            for bi in range(self.nb):
                out[f"s[{self.D['bands'][bi]}]"] = float(sb[bi]) if present[bi] else np.nan
            for j in range(1, self.nR):
                for k in range(self.nW):
                    out[f"g({self.rk[j]:g},{self.wk[k]:g})"] = float(G[j, k])
            for k in range(1, self.nW):
                out[f"h({self.wk[k]:g})"] = float(H[k])
            for bi in range(self.nb):
                for r in rain:
                    out[f"Weq[{self.D['bands'][bi]}]({r:g})"] = float(surface_path(G, H, self.rk, self.wk, sb[bi], r, 0.0, 2000, self.nsub)[-1]) if present[bi] else np.nan
            b0 = int(np.argmax(present))
            up = surface_path(G, H, self.rk, self.wk, sb[b0], 2.0, 0.0, 2000, self.nsub)
            out["t_wet"] = crossing_time(np.r_[0.0, up], 0.632 * up[-1]) if up[-1] > 0 else np.nan
            out["t_dry"] = crossing_time(np.r_[1.0, surface_path(G, H, self.rk, self.wk, sb[b0], 0.0, 1.0, 2000, self.nsub)], 0.368)
        return out


def fit_all2(P, starts, lags, maxiter, em_rounds=0, state=None, checkpoint=None):
    """The fit of a step S2 model on a whole dataset (config `maxiter`): every start at d = 0, the lag profile from the
    best one, a polish at the best d, and for a model with drift the EM rounds for s_b.
    state (a dict, optional) keeps the finished starts and lag fits; checkpoint() is called after each one, so that an
    interrupted fit goes on from there (step S2b)."""
    ones = np.ones(P.ne)
    state = {} if state is None else state
    first = state.setdefault("first", {})
    for s, v in starts.items():                            # a start: config dict or x
        if s not in first:
            first[s] = P.fit(v if isinstance(v, np.ndarray) else P.start(v), 0, ones, maxiter["start"])
            if checkpoint:
                checkpoint()
    prof = state.setdefault("prof", {})
    prof[0] = min(first.values(), key=lambda r: r["loss"])
    for step in (1, -1):
        for d in range(step, (lags[1] if step > 0 else lags[0]) + step, step):
            if d not in prof:
                prof[d] = P.fit(prof[d - step]["x"], d, ones, maxiter["lag"])
                if checkpoint:
                    checkpoint()
    d_best = min(prof, key=lambda d: prof[d]["loss"])
    best = P.fit(prof[d_best]["x"], d_best, ones, maxiter["polish"])
    s_b = [P.s_b]
    if P.drift:
        for _ in range(em_rounds):
            s_b.append(P.update_s_b(best["x"]))
            best = P.fit(best["x"], d_best, ones, maxiter["em"])
    return {"d": int(d_best), "x": best["x"], "loss": best["loss"], "nit": best["nit"], "s_b": s_b,
            "starts_d0": {s: r["loss"] for s, r in first.items()}, "lag_profile": {int(d): r["loss"] for d, r in prof.items()},
            "converged": best.get("converged"), "starts_converged": {s: r.get("converged") for s, r in first.items()},
            "starts_nit": {s: r["nit"] for s, r in first.items()}}


def fit_weighted2(P, allfit, wle, maxiter):
    """A fit on a re-weighted training set (a leave-one-event-out fold or a bootstrap draw): warm start from the fit on
    all events, with its lag. Links without a training link-event get a_link = 1, such link-events a drift of 0."""
    best = P.fit(allfit["x"], allfit["d"], wle, maxiter)
    seen = np.bincount(P.D["le_link"], weights=wle, minlength=P.nl) > 0
    la = best["x"][P.names == "loga"]
    la[~seen] = 0.0
    best["x"][P.names == "loga"] = la
    if P.drift:
        b = best["x"][P.names == "b"]
        b[wle == 0] = 0.0
        best["x"][P.names == "b"] = b
    return best


def drift_allowed(D, pred, tau_h, s_b, delta):
    """The held-out prediction after one linear drift per link-event is allowed (config `judge.drift_allowed`): the
    ridge estimate of the slope from the residuals, prior sd s_b, residual sd delta. Returns (prediction, slopes)."""
    w = D["w5"]
    ne = len(D["le_n"])
    num = np.bincount(D["bin_le"], weights=w * tau_h * (D["y5"] - pred), minlength=ne)
    den = np.bincount(D["bin_le"], weights=w * tau_h ** 2, minlength=ne) + (delta / s_b) ** 2
    slope = num / den
    return pred + slope[D["bin_le"]] * tau_h, slope


# =================================================================================================== step S2: naming the M5 surfaces
def m5_states(P5, x, d):
    """The 1 min states of a fit of M5: the rain R(t - d) and the store W(t) of every minute of every link-event, from a
    1 min run of the same model and parameters (W = A_pred - c - L P(R))."""
    D = P5.D
    P1 = Problem2(D, P5.model, P5.M, agg=1)
    pred = P1.predict(x, d)
    Pinc = P1.unpack(x)[1]
    cP = np.concatenate([np.zeros((P1.nb, 1)), np.cumsum(Pinc, axis=1)], axis=1)
    R, W = [], []
    for e in range(P1.ne):
        n, o, ro = D["le_n"][e], D["le_boff1"][e], D["le_roff"][e] + PADM - d
        r = D["rain"][ro:ro + n]
        R.append(r)
        W.append(pred[o:o + n] - D["le_c"][e] - D["le_L"][e] * np.interp(r, P1.knots, cP[D["le_band"][e]]))
    return np.concatenate(R), np.concatenate(W)


def knot_visits(R, W, rk, wk):
    """How often the data visit each knot of the surfaces (bilinear shares of the 1 min states): for g the minutes with
    rain on the (R, W) grid, for h every minute on the W knots."""
    def cell(v, k):
        v = np.clip(v, k[0], k[-1])
        i = np.clip(np.searchsorted(k, v, side="right") - 1, 0, len(k) - 2)
        return i, (v - k[i]) / (k[i + 1] - k[i])

    iw, fw = cell(W, wk)
    hv = np.zeros(len(wk))
    np.add.at(hv, iw, 1 - fw)
    np.add.at(hv, iw + 1, fw)
    wet = R > 0
    ir, fr = cell(R[wet], rk)
    gv = np.zeros((len(rk), len(wk)))
    for i, wi in ((ir, 1 - fr), (ir + 1, fr)):
        for k, wk_ in ((iw[wet], 1 - fw[wet]), (iw[wet] + 1, fw[wet])):
            np.add.at(gv, (i, k), wi * wk_)
    return gv, hv


def naming_terms(NM):
    """The term libraries of config `naming` as (name, function) pairs: g(R, W) terms and h(W) terms.
    exp(-R / R0) of the config is used as 1 - exp(-R / R0), so that every R-term is 0 at R = 0 (g(0, W) = 0 in M5)."""
    r_terms = [("R", lambda R, W: R), ("R^2", lambda R, W: R ** 2)]
    r_terms += [(f"R/(1+R/{r0:g})", lambda R, W, r0=r0: R / (1 + R / r0)) for r0 in NM["r0_mmh"]]
    e_terms = [(f"(1-exp(-R/{r0:g}))", lambda R, W, r0=r0: 1 - np.exp(-R / r0)) for r0 in NM["r0_mmh"]]
    g = r_terms + [("W", lambda R, W: W + 0 * R), ("R W", lambda R, W: R * W)] + e_terms
    g += [(f"{n} (1-W/{wm:g})", lambda R, W, f=f, wm=wm: f(R, W) * (1 - W / wm)) for n, f in r_terms + e_terms for wm in NM["wm_db"]]
    h = [("W", lambda W: W), ("W^2", lambda W: W ** 2)] + [(f"max(W-{c:g},0)", lambda W, c=c: np.maximum(W - c, 0.0)) for c in NM["h_switch_db"]]
    return g, h


def esindy(theta, y, wts, cfg, rng):
    """Bagged sparse regression (E-SINDy) of samples y on the library columns theta, with the backward elimination of
    sindy.py. Rows are weighted by wts (sqrt of the share); a bag = a bootstrap sample of the rows. For every lambda of
    the grid: the selection frequency of each term over the bags, the terms kept (frequency >= keep_share), their
    coefficients refitted on all rows and the relative weighted error. The pick: the formula with the fewest terms whose
    error is within 10 % of the best error over the grid."""
    sw = np.sqrt(wts / wts.sum())
    G, b = theta * sw[:, None], y * sw
    rows = []
    for lam in cfg["lambda_grid"]:
        hits = np.zeros(theta.shape[1])
        for _ in range(cfg["bags"]):
            ix = rng.integers(0, len(y), len(y))
            hits += sindy.select(G[ix], b[ix], lam) != 0
        keep = hits / cfg["bags"] >= cfg["keep_share"]
        xi = sindy.refit(G, b, keep)
        rows.append({"lambda": float(lam), "freq": hits / cfg["bags"], "keep": keep, "xi": xi,
                     "terms": int(keep.sum()), "rel_error": float(np.linalg.norm(b - G @ xi) / max(np.linalg.norm(b), 1e-300))})
    best = min(r["rel_error"] for r in rows)
    ok = [r for r in rows if r["rel_error"] <= 1.1 * best and r["terms"] > 0] or rows
    return rows, min(ok, key=lambda r: (r["terms"], r["rel_error"]))


class Problem2n(Problem2):
    """M5n: M5 with the surfaces given by a named formula, g(R, W) = sum_i theta_i g_i(R, W) and h(W) = sum_j eta_j h_j(W)
    (the coefficients free; g is the formula at the knots with R > 0 and 0 at R = 0, bilinear in between, as M5).
    P, s_band, a_link and the lag as M5; no smoothness penalty. g_terms / h_terms = (name, function) pairs."""

    def __init__(self, D, M, g_terms, h_terms, agg=5):
        super().__init__(D, "M5", M, agg=agg)
        Rg, Wg = np.meshgrid(self.rk, self.wk, indexing="ij")
        self.Bg = np.zeros((len(g_terms), self.nR, self.nW))
        for i, (_, f) in enumerate(g_terms):
            self.Bg[i, 1:] = f(Rg[1:], Wg[1:])
        self.Bh = np.array([f(self.wk) for _, f in h_terms]).reshape(len(h_terms), self.nW)
        old = [(n, bd) for n, bd in zip(self.names, self.bounds) if n not in ("logG", "logH")]
        cut = sum(n == "logP" for n, _ in old)
        new = old[:cut] + [("cg", (-10.0, 10.0))] * len(g_terms) + [("ch", (-10.0, 10.0))] * len(h_terms) + old[cut:]
        self.names, self.bounds = np.array([n for n, _ in new]), [bd for _, bd in new]

    def unpack(self, x):
        get = lambda n: x[self.names == n]                 # noqa: E731
        Pinc = np.exp(get("logP")).reshape(self.nb, self.nk)
        G, H = np.tensordot(get("cg"), self.Bg, 1), get("ch") @ self.Bh
        sb = np.ones(self.nb)
        sb[1:] = np.exp(get("logs"))
        return np.zeros((self.nb, self.nk)), Pinc, np.exp(get("loga")), np.ones((self.nb, 3)), np.ones(self.nb), G, H, sb, np.zeros(self.ne)

    def run(self, x, d, wle, want_grad=True):
        D = self.D
        Finc, Pinc, a, tau, phi, G, H, sb, b = self.unpack(x)
        gF, gP, ga = np.zeros_like(Finc), np.zeros_like(Pinc), np.zeros(self.nl)
        gtau, gphi, gG, gH, gsb, gb = np.zeros((self.nb, 3)), np.zeros(self.nb), np.zeros_like(G), np.zeros_like(H), np.zeros(self.nb), np.zeros(self.ne)
        loss = kernel2(3, int(d), Finc, Pinc, a, tau, phi, G, H, self.rk, self.wk, sb, self.nsub, b, wle, D["le_roff"], D["le_n"], self.boff,
                       D["le_link"], D["le_band"], D["le_L"], D["le_c"], self.t0, D["rain"], self.y, self.w, self.agg, self.M["huber_delta_db"],
                       self.knots, self.pred, self.le_loss, gF, gP, ga, gtau, gphi, gG, gH, gsb, gb, want_grad, False)
        la = x[self.names == "loga"]
        loss += self.pen * float(la @ la)
        g = np.zeros(len(x))
        if not want_grad:
            return loss, g
        g[self.names == "logP"] = (gP * Pinc).ravel()
        g[self.names == "loga"] = ga * a + 2 * self.pen * la
        g[self.names == "cg"] = np.tensordot(self.Bg, gG, ([1, 2], [0, 1]))
        g[self.names == "ch"] = self.Bh @ gH
        g[self.names == "logs"] = gsb[1:] * sb[1:]
        return loss, g

    def start_from(self, P5, x5, cg, ch):
        """A start from a fit of M5 (its P, s_band, a_link) and the coefficients of the named formula."""
        x = np.zeros(len(self.names))
        for n in ("logP", "logs", "loga"):
            x[self.names == n] = x5[P5.names == n]
        x[self.names == "cg"], x[self.names == "ch"] = cg, ch
        return x


# =================================================================================================== data and fake data of S2 / S2b
def merged_bands(d, DATA):
    """Config `min_band_links`: a band with too few TRAIN links below max_length_km is merged with the nearest band."""
    links = d.links()
    links = links[(links.split_link == "train") & (links.length_m < DATA["max_length_km"] * 1000)]
    freq = d.sublinks().groupby("cml_id").frequency_ghz.mean().reindex(links.cml_id)
    bands = {k: list(v) for k, v in DATA["bands_ghz"].items()}
    notes = []
    while len(bands) > 1:
        n = {k: int(((freq >= lo) & (freq < hi)).sum()) for k, (lo, hi) in bands.items()}
        small = [k for k in bands if n[k] < DATA["min_band_links"]]
        if not small:
            break
        k = small[0]
        centre = {q: float(freq[(freq >= lo) & (freq < hi)].median()) if n[q] else np.mean([lo, hi]) for q, (lo, hi) in bands.items()}
        other = min((q for q in bands if q != k), key=lambda q: abs(centre[q] - centre[k]))
        bands[other] = [min(bands[other][0], bands[k][0]), max(bands[other][1], bands[k][1])]
        notes.append(f"{k} ({n[k]} links) merged with {other}")
        del bands[k]
    return bands, notes, {k: int(((freq >= lo) & (freq < hi)).sum()) for k, (lo, hi) in bands.items()}


def gauge_variants(D, has_g):
    """The link-events with a City gauge: (the dataset with the GAUGE rain as input, the dataset with the RADAR rain
    on the same link-events and the same 5 min points)."""
    out = []
    for use_gauge in (True, False):
        Dc = subset(D, has_g)
        g_nan = np.isnan(Dc["rain_gauge"])
        bad5 = np.concatenate([g_nan[o + PADM:o + PADM + n].reshape(-1, 5).any(axis=1) for o, n in zip(Dc["le_roff"], Dc["le_n"])])
        Dc["y5"], Dc["y1"] = np.where(Dc["w5"] > 0, Dc["y5"], np.nan), np.where(Dc["w1"] > 0, Dc["y1"], np.nan)
        Dc["ok5"] = Dc["ok5"] & ~bad5
        if use_gauge:
            Dc["rain"] = np.nan_to_num(Dc["rain_gauge"])
        out.append(finish(Dc, D["bands"]))
    return out


def noise_table():
    """The noise of every sub-link from the simulator table: {sublink: (white sd, slow sd, slow time), "median": ...}."""
    nz = pd.read_csv(PATHS.SIM_NOISE).set_index("sublink_id")
    noise = {s_: (r.sigma_white_db, r.sigma_slow_db, r.tau_noise_min) for s_, r in nz.iterrows()}
    noise["median"] = (float(nz.sigma_white_db.median()), float(nz.sigma_slow_db.median()), float(nz.tau_noise_min.median()))
    return noise


def truth_curves(le, bands, FK):
    """The truth of the fake data: F(R) (all bands) and the ITU-R P.838 power law per band (fake data only)."""
    freq = le.groupby("band").frequency_ghz.median()
    itu = {int(b): itu_p838(float(freq[b]), "v") for b in freq.index}
    F = lambda r: FK["f_top_db"] * (1 - np.exp(-np.asarray(r, float) / FK["f_scale_mmh"]))   # noqa: E731
    return F, itu


def make_fake(D, le, spec, seed, noise, FK, t_ref_min, offset_min_minutes):
    """A fake dataset on the link-events of D (rows of `le`): the truth `spec` (config `fake.truths`) driven by the
    rain of D, observed with noise, quantization and an unknown offset; c by the same rule as the real data.
    t_ref_min = model.drift.t_ref_min_before_rain, offset_min_minutes = data.offset_min_minutes."""
    rng = np.random.default_rng(seed)
    a = np.exp(rng.normal(0, FK["a_sd"], D["nl"]))
    F, itu = truth_curves(le, D["bands"], FK)
    y1, c = np.full(len(D["y1"]), np.nan), np.zeros(len(le))
    f_var = np.sqrt(FK["noise_variance_factor"])
    lag = spec["lag_min"]
    for e, r in enumerate(le.itertuples()):
        n, o1 = D["le_n"][e], D["le_boff1"][e]
        rain = D["rain"][D["le_roff"][e] + PADM - lag:D["le_roff"][e] + PADM - lag + n]
        u = a[D["le_link"][e]] * F(rain)
        if spec["model"] in ("M4", "M4b"):
            wv = simulate_stores(u, spec["phi"], spec["tau_f_wet"], spec["tau_f_dry"], spec["tau_s"])
        else:
            wv = simulate_stores(u, 1.0, spec["tau_wet"], spec["tau_dry"], 1.0)
        truth_a = wv + r.L_km * itu[r.band][0] * np.power(rain, itu[r.band][1])
        if spec.get("drift_sd_db_h"):
            t0 = D["le_first"][e] * 5 - t_ref_min
            truth_a = truth_a + rng.normal(0, spec["drift_sd_db_h"]) * (np.arange(n) - t0) / 60.0
        sw, ss, tn = noise.get(r.sublink_id, noise["median"])
        off = rng.standard_normal() * (r.offset_std if np.isfinite(r.offset_std) else 0.0)
        slow = np.repeat(ar1(rng.standard_normal(n), ss * f_var, tn, rng.standard_normal() * ss * f_var), 6) if ss > 0 else 0.0
        x = np.repeat(truth_a, 6) + off + slow + sw * f_var * rng.standard_normal(n * 6)
        obs = quantize(x, rng.uniform(0, FK["rsl_step_db"]), FK["rsl_step_db"]).reshape(n, 6).mean(axis=1)
        obs[D["w1"][o1:o1 + n] == 0] = np.nan                      # the real missing minutes
        y1[o1:o1 + n] = obs
        pre = D["pre1"][o1:o1 + n] & np.isfinite(obs)
        c[e] = float(np.median(obs[pre])) if pre.sum() >= offset_min_minutes else 0.0
    Fd = dict(D)
    Fd["y1"], Fd["le_c"] = y1, c
    m5 = y1.reshape(-1, 5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        Fd["y5"] = np.where(np.isfinite(m5).sum(axis=1) >= 3, np.nanmean(m5, axis=1), np.nan)
    return finish(Fd, D["bands"])


def fake_from_truth(D, le, truth1, seed, noise, FK, offset_min_minutes):
    """A fake dataset with a GIVEN 1 min truth of A (one value per minute of every window of D, c = 0), observed with the
    noise of `make_fake` (white + AR(1) per sub-link, variance x FK.noise_variance_factor, quantization of the 10 s values,
    a true offset, the real missing minutes); c by the data rule. Step S5a (truth = a model with data-learned curves)."""
    rng = np.random.default_rng(seed)
    y1, c = np.full(len(D["y1"]), np.nan), np.zeros(len(le))
    f_var = np.sqrt(FK["noise_variance_factor"])
    for e, r in enumerate(le.itertuples()):
        n, o1 = D["le_n"][e], D["le_boff1"][e]
        sw, ss, tn = noise.get(r.sublink_id, noise["median"])
        off = rng.standard_normal() * (r.offset_std if np.isfinite(r.offset_std) else 0.0)
        slow = np.repeat(ar1(rng.standard_normal(n), ss * f_var, tn, rng.standard_normal() * ss * f_var), 6) if ss > 0 else 0.0
        x = np.repeat(truth1[o1:o1 + n], 6) + off + slow + sw * f_var * rng.standard_normal(n * 6)
        obs = quantize(x, rng.uniform(0, FK["rsl_step_db"]), FK["rsl_step_db"]).reshape(n, 6).mean(axis=1)
        obs[D["w1"][o1:o1 + n] == 0] = np.nan
        y1[o1:o1 + n] = obs
        pre = D["pre1"][o1:o1 + n] & np.isfinite(obs)
        c[e] = float(np.median(obs[pre])) if pre.sum() >= offset_min_minutes else 0.0
    Fd = dict(D)
    Fd["y1"], Fd["le_c"] = y1, c
    m5 = y1.reshape(-1, 5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        Fd["y5"] = np.where(np.isfinite(m5).sum(axis=1) >= 3, np.nanmean(m5, axis=1), np.nan)
    return finish(Fd, D["bands"])


def fix_params(Pm, fixed):
    """Hold dynamic parameters of a Problem2 fixed (step S3b variant B, S5a): equal lower and upper bounds at the value
    (logit for phi, log for a time constant)."""
    for p, v in fixed.items():
        val = float(np.log(v / (1 - v))) if p == "phi" else float(np.log(v))
        for i in np.nonzero(Pm.names == p)[0]:
            Pm.bounds[i] = (val, val)


def x_into(Pm, Psrc, xsrc, dyn=None, curves=True):
    """A parameter vector of Psrc carried into Pm (another model with the same links): P, F (if curves; the same knots)
    and a_link copied by name, the dynamic parameters of Pm from `dyn` (readable values) or, where Psrc has them, from xsrc."""
    x = np.zeros(len(Pm.names))
    for n in ((Pm.pname, Pm.fname) if curves else ()) + ("loga",):
        if (Pm.names == n).any():
            x[Pm.names == n] = xsrc[Psrc.names == n]
    for p in Pm.dyn:
        if dyn is not None and p in dyn:
            v = dyn[p]
            x[Pm.names == p] = np.log(v / (1 - v)) if p == "phi" else np.log(v)
        else:
            x[Pm.names == p] = xsrc[Psrc.names == p]
    return x


# =================================================================================================== judge helpers of S2 / S2b
def paired(a, b, n_ev, J):
    """Mean and 90 % event-bootstrap interval of the paired difference a - b of two per-event error columns."""
    draws = np.random.default_rng(J["seed"]).integers(0, n_ev, (J["boot_reps"], n_ev))
    dif = np.nanmean(a[draws], axis=1) - np.nanmean(b[draws], axis=1)
    return float(np.nanmean(a) - np.nanmean(b)), *map(float, np.nanpercentile(dif, J["interval"]))


def shower_stat(ev, err, n_ev, J):
    """Median of the per-shower peak errors with a 90 % interval from resampling the events."""
    if len(err) == 0:
        return np.nan, np.nan, np.nan, 0
    by_ev = [err[ev == e] for e in range(n_ev)]
    draws = np.random.default_rng(J["seed"]).integers(0, n_ev, (J["boot_reps"], n_ev))
    med = [np.median(np.concatenate([by_ev[e] for e in dr])) for dr in draws]
    return float(np.median(err)), *map(float, np.percentile(med, J["interval"])), len(err)


def after_curves(D, preds, J, min_peak_db):
    """Median (A - c) / observed peak of the last shower, per 5 min after the last wet bin, observed and predicted,
    over the link-events whose last shower has an observed peak of at least min_peak_db."""
    n_after = J["after_windows_min"][-1][1] // 5
    curves = {k: [] for k in ["observed"] + list(preds)}
    for e in range(len(D["le_n"])):
        o, n, c = D["le_boff5"][e], D["le_n"][e] // 5, D["le_c"][e]
        w = D["w5"][o:o + n] > 0
        obs = np.where(w, D["y5"][o:o + n] - c, np.nan)
        sh = showers(D, e, J["shower_gap_min"])
        last = D["le_last"][e]
        if not sh or not np.isfinite(obs[sh[-1][0]:last + 1]).any():
            continue
        peak = float(np.nanmax(obs[sh[-1][0]:last + 1]))
        if peak < min_peak_db:
            continue
        pad = lambda v: np.r_[v[last + 1:last + 1 + n_after], np.full(n_after, np.nan)][:n_after] / peak   # noqa: E731
        curves["observed"].append(pad(obs))
        for k, p in preds.items():
            curves[k].append(pad(np.where(w, p[o:o + n] - c, np.nan)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return pd.DataFrame({"minutes_after_last_wet_bin": np.arange(n_after) * 5, "link_events": len(curves["observed"]),
                             **{k: np.nanmedian(np.array(v), axis=0) for k, v in curves.items()}})
