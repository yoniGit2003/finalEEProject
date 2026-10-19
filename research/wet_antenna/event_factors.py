"""Library of steps S4.1 and S4.3: per-event factors on a frozen S3b fit, the fake data for them, and the radar / gauge ratio
per event.

    A_pred = alpha_e * W[a_link F(R / k_e) dynamics] + L * P(R(t - d) / k_e) + c

k_e = the rain-input scale of event e, alpha_e = its antenna factor. W and L * P are simulated separately for every k on a
log grid (alpha only multiplies W), so the Huber loss of an event is known on the whole (k, alpha) grid for every link;
the estimate is the grid minimum with a weak log-ridge toward 1, and a link bootstrap inside the event gives the interval.
Settings: config `wet_antenna.conditions_s4` (factors, fake noise of system_id_s2.fake).
Moved unchanged from conditions_s4.py (step S4.3 uses the same code).
"""
import copy
import pickle

import numpy as np
import pandas as pd
from numba import njit, prange
from scipy import stats

from research.common.config import load_config
from research.wet_antenna import length_test as LTL
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.simulator import ar1, quantize

CFG = load_config()
WA = CFG["wet_antenna"]
S2, S2B, S3B, C = WA["system_id_s2"], WA["system_id_s2b"], WA["system_id_s3b"], WA["conditions_s4"]
M = copy.deepcopy(S2["model"])                                # the S2b model settings (as S3b)
M["increments"], M["maxiter"] = S2B["model"]["increments"], dict(S2B["model"]["maxiter"])
M["bounds"] = {**M["bounds"], "increment": S2B["model"]["increment_bounds"], "tau_s_min": S2B["model"]["tau_s_m4b_min"]}
FK = S2["fake"]
DELTA = M["huber_delta_db"]
FA = C["factors"]
FITS = {k: tuple(v) for k, v in S3B["fits"].items()}           # "M4b-B" -> ("M4b", "B")

# =================================================================================================== the loss grid
@njit(cache=True, parallel=True)
def grid_loss(W, Pp, base, link, agrid, delta, nl):
    """Huber loss summed per link on the (k, alpha) grid: residual = alpha W[k] + Pp[k] + base (base = c - y)."""
    nk, nb = W.shape
    na = len(agrid)
    out = np.zeros((nk, na, nl))
    for i in prange(nk):
        for j in range(nb):
            wv, pb, lk = W[i, j], Pp[i, j] + base[j], link[j]
            for q in range(na):
                r = agrid[q] * wv + pb
                ar = abs(r)
                out[i, q, lk] += 0.5 * r * r if ar <= delta else delta * (ar - 0.5 * delta)
    return out


def factor_grid(bounds):
    """The log-spaced grid over the bounds, plus the value 1 itself."""
    g = np.geomspace(bounds[0], bounds[1], FA["grid_points"])
    return np.unique(np.r_[g, 1.0])


KG, AG = factor_grid(FA["k_bounds"]), factor_grid(FA["alpha_bounds"])
LKG, LAG = np.log(KG), np.log(AG)
I1K, I1A = int(np.argmin(np.abs(KG - 1))), int(np.argmin(np.abs(AG - 1)))


def frozen(name, Dt, let):
    """A frozen S3b fit: x, lag d, a_link by cml_id (from the TRAIN all-length dataset it was fitted on)."""
    r = pickle.loads((P.SID3B_CACHE / f"fit_{name}.pkl").read_bytes())
    return {"x": r["all"]["x"], "d": int(r["all"]["d"]), "a": dict(zip(LTL.link_ids(Dt, let), r["a_link"])), "model": FITS[name][0], "params": r["params"]}


def bin_positions(D, idx):
    """The 5 min bin positions in D of the link-events idx (in order)."""
    return np.concatenate([np.arange(D["le_boff5"][i], D["le_boff5"][i] + D["le_n"][i] // 5) for i in idx])


def bin_times(D, le):
    """The start time of every 5 min bin of D."""
    n5 = D["le_n"] // 5
    start = np.repeat(pd.to_datetime(le.window_start).to_numpy(), n5)
    off = np.arange(len(D["y5"])) - np.repeat(D["le_boff5"], n5)
    return pd.DatetimeIndex(start + off * np.timedelta64(5, "m"))


class EventParts:
    """One event of a dataset under a frozen fit: its sub-dataset and W, L * P (c = 0) for every k of the grid."""

    def __init__(self, D, le, f, idx):
        self.idx = idx
        self.De = sid.subset(D, np.isin(np.arange(len(le)), idx))
        self.le = le.iloc[idx].reset_index(drop=True)
        self.D0 = dict(self.De)
        self.D0["le_c"] = np.zeros(len(idx))
        self.rain0 = self.De["rain"].copy()
        self.f = f
        self.nl = int(self.De["nl"])
        self.link_bin = self.De["le_link"][self.De["bin_le"]]
        Pm = sid.Problem2(self.D0, f["model"], M)
        x, _ = LTL.x_for(Pm, f, LTL.link_ids(self.De, self.le))
        self.xW, self.xP = x.copy(), x.copy()
        self.xW[Pm.names == Pm.pname] = 0.0                # W alone: no path term
        self.xP[Pm.names == Pm.fname] = 0.0                # L * P alone: no antenna term
        self.W, self.Pp = np.empty((len(KG), len(self.De["y5"]))), np.empty((len(KG), len(self.De["y5"])))
        for i, k in enumerate(KG):
            self.D0["rain"] = self.rain0 / k
            self.W[i], self.Pp[i] = Pm.predict(self.xW, f["d"]), Pm.predict(self.xP, f["d"])
        self.D0["rain"] = self.rain0

    def truth_1min(self, k, alpha):
        """alpha W + L P at 1 min resolution for the input R / k (c = 0): the truth of the fake data."""
        Pm = sid.Problem2(self.D0, self.f["model"], M, agg=1)
        self.D0["rain"] = self.rain0 / k
        out = alpha * Pm.predict(self.xW, self.f["d"]) + Pm.predict(self.xP, self.f["d"])
        self.D0["rain"] = self.rain0
        return out

    def grid(self, y5, w5, c_le):
        """The loss per link on the (k, alpha) grid for observations y5 (weights w5) and offsets c_le."""
        v = w5 > 0
        base = (c_le[self.De["bin_le"]] - y5)[v]
        return grid_loss(np.ascontiguousarray(self.W[:, v]), np.ascontiguousarray(self.Pp[:, v]), base, self.link_bin[v], AG, DELTA, self.nl)

    def pred(self, ik, ia, c_le):
        return AG[ia] * self.W[ik] + self.Pp[ik] + c_le[self.De["bin_le"]]


def estimate(G, lam, rng=None, boot=0):
    """Grid minimum of the summed link losses + the log-ridge; optional link bootstrap. Returns a dict of grid indices."""
    nk, na, nl = G.shape
    pen = lam * (LKG[:, None] ** 2 + LAG[None, :] ** 2)
    tot = G.sum(axis=2) + nl * pen
    ik, ia = np.unravel_index(int(np.argmin(tot)), tot.shape)
    out = {"ik": int(ik), "ia": int(ia),
           "ik_refit": int(np.argmin(G[:, I1A, :].sum(axis=1) + nl * lam * LKG ** 2)),       # k alone (alpha = 1)
           "ia_refit": int(np.argmin(G[I1K, :, :].sum(axis=1) + nl * lam * LAG ** 2))}       # alpha alone (k = 1)
    if boot:
        mult = np.stack([np.bincount(rng.integers(0, nl, nl), minlength=nl) for _ in range(boot)]).astype(float)
        tb = (G.reshape(-1, nl) @ mult.T).T.reshape(boot, nk, na) + nl * pen[None]
        flat = tb.reshape(boot, -1).argmin(axis=1)
        out["boot_ik"], out["boot_ia"] = flat // na, flat % na
    return out


def fake_obs(ev, truth1, rng, noise):
    """Fake 5 min observations of one event: the noise of system_id.make_fake (white + AR(1) per sub-link, variance x 3,
    0.3 dB quantization of the 10 s values, a true offset, the real missing minutes) and c by the data rule."""
    De, le = ev.De, ev.le
    f_var = np.sqrt(FK["noise_variance_factor"])
    y1, c = np.full(len(De["y1"]), np.nan), np.zeros(len(le))
    for e, r in enumerate(le.itertuples()):
        n, o1 = De["le_n"][e], De["le_boff1"][e]
        sw, ss, tn = noise.get(r.sublink_id, noise["median"])
        off = rng.standard_normal() * (r.offset_std if np.isfinite(r.offset_std) else 0.0)
        slow = np.repeat(ar1(rng.standard_normal(n), ss * f_var, tn, rng.standard_normal() * ss * f_var), 6) if ss > 0 else 0.0
        x = np.repeat(truth1[o1:o1 + n], 6) + off + slow + sw * f_var * rng.standard_normal(n * 6)
        obs = quantize(x, rng.uniform(0, FK["rsl_step_db"]), FK["rsl_step_db"]).reshape(n, 6).mean(axis=1)
        obs[De["w1"][o1:o1 + n] == 0] = np.nan
        y1[o1:o1 + n] = obs
        pre = De["pre1"][o1:o1 + n] & np.isfinite(obs)
        c[e] = float(np.median(obs[pre])) if pre.sum() >= S2["data"]["offset_min_minutes"] else 0.0
    m5 = y1.reshape(-1, 5)
    with np.errstate(all="ignore"):
        y5 = np.where(np.isfinite(m5).sum(axis=1) >= 3, np.nanmean(m5, axis=1), np.nan)
    w5 = (De["ok5"] & np.isfinite(y5)).astype(float)
    return np.nan_to_num(y5), w5, c


# =================================================================================================== radar / gauge, statistics
def rescaled(D, le, scale_ev):
    """The dataset with the radar rain R / g of its event, in the model input (1 min) and in the 5 min radar of the bins
    (moved unchanged from conditions_s4c.py; used by S4.3 and S5a)."""
    sc = le.event_id.map(scale_ev).to_numpy(float)
    D2 = dict(D)
    D2["rain"] = D["rain"] / np.repeat(sc, D["le_n"] + 2 * sid.PADM)
    D2["r5"] = D["r5"] / np.repeat(sc, D["le_n"] // 5)
    return D2


def radar_gauge_primary(D, le):
    """Per event: radar path-mean / City gauge over the gauge link-events (one sub-link per link), both valid minutes."""
    rows = []
    PADM = sid.PADM
    has = le[le.gauge_id.astype(str) != ""].drop_duplicates(["event_id", "cml_id"])
    for e_id, g in has.groupby("event_id", sort=False):
        sr, sg, n_l = 0.0, 0.0, 0
        for i in g.index:
            o, n, o5 = D["le_roff"][i], D["le_n"][i], D["le_boff5"][i]
            rad = D["rain"][o + PADM:o + PADM + n]
            gau = D["rain_gauge"][o + PADM:o + PADM + n]
            ok = np.repeat(np.isfinite(D["r5"][o5:o5 + n // 5]), 5) & np.isfinite(gau)
            sr += rad[ok].sum() / 60
            sg += gau[ok].sum() / 60
            n_l += 1
        rows.append({"event_id": e_id, "gauge_links": n_l, "radar_mm_sum": sr, "gauge_mm_sum": sg, "gauge_mm_per_link": sg / n_l,
                     "radar_over_gauge": sr / sg if sg > 0 else np.nan})
    return pd.DataFrame(rows)


def radar_gauge_pixel(d, ev_tab):
    """Per event: radar pixel at the 10 City gauges / the gauges (5 min), over [start, end), both valid (information)."""
    rad = d.store("radar_at_gauges_5min").R_pixel.load().to_pandas()
    gau = d.store("gauges_city_5min").rate_mm_h.load().to_pandas()
    cols = [c for c in gau.columns if c in rad.columns]
    rad, gau = rad[cols], gau[cols].reindex(rad.index)
    rows = []
    for e_id, r in ev_tab.iterrows():
        m = (rad.index >= r.start) & (rad.index < r.end)
        a, b = rad[m].to_numpy(float), gau[m].to_numpy(float)
        ok = np.isfinite(a) & np.isfinite(b)
        sa, sb = a[ok].sum() * 5 / 60, b[ok].sum() * 5 / 60
        rows.append({"event_id": e_id, "pixel_radar_mm": sa, "pixel_gauge_mm": sb, "pixel_radar_over_gauge": sa / sb if sb > 0 else np.nan})
    return pd.DataFrame(rows)


def spearman_boot(x, y, reps, rng):
    """Spearman rho, its p-value and a 90 % event-bootstrap interval."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 5:
        return np.nan, np.nan, np.nan, np.nan, len(x)
    rho, p = stats.spearmanr(x, y)
    b = []
    for _ in range(reps):
        i = rng.integers(0, len(x), len(x))
        if np.std(x[i]) > 0 and np.std(y[i]) > 0:
            b.append(stats.spearmanr(x[i], y[i])[0])
    lo, hi = np.nanpercentile(b, [5, 95])
    return float(rho), float(p), float(lo), float(hi), len(x)


def bh(p, q):
    """Benjamini-Hochberg: True where the test is a discovery at level q."""
    p = np.asarray(p, float)
    o = np.argsort(p)
    m = len(p)
    passed = p[o] <= q * np.arange(1, m + 1) / m
    k = np.nonzero(passed)[0].max() + 1 if passed.any() else 0
    out = np.zeros(m, bool)
    out[o[:k]] = True
    return out
