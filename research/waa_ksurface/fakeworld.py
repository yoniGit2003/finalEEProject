"""Library of the A3 rescue: the fake world in which a length test is tried before it is used on real data.

The fake world keeps what is real (the links, the storms, the rain of every link-event) and replaces the attenuation by a
KNOWN answer, so that a method can be judged by whether it returns that answer:

    true rain R_true      = the real rain input of the link-event (5 min)
    W(t)                  = one store relaxing to F_true(R_true), slower when it dries than when it wets (round numbers)
    A - c                 = W + L * P_true(R_true) + noise (measured in this branch on dry minutes)
    observed rain R_obs   = R_true * exp(error), the error model measured in step E1:
                            a scale part per link-event + a part shared by all links at that time + a part of the link's own

Nothing here comes from Part 1: the truth is generic, the noise and the rain error are this branch's own measurements.
Settings: config section `fake`.
"""
import numpy as np
import pandas as pd
from numba import njit

from research.waa_ksurface import data


def config():
    return data.branch_config()["fake"]


# =================================================================================================== the truth
def f_true(r, cfg):
    """The true steady WAA of the fake world (dB): saturating, round numbers."""
    r = np.maximum(np.nan_to_num(np.asarray(r, float)), 0.0)
    return cfg["truth"]["f_top_db"] * r / (r + cfg["truth"]["f_half_mmh"])


def p_true(band, r, cfg):
    """The true path loss per km of the fake world (dB/km): linear in R."""
    return cfg["truth"]["p_per_mmh"][band] * np.maximum(np.nan_to_num(np.asarray(r, float)), 0.0)


@njit(cache=True)
def _store(target5, boff, n5, k_wet, k_dry):
    """The 5 min mean of W for a target that is constant over each 5 min bin: one exact step per minute, W = 0 at the start."""
    out = np.empty(len(target5))
    for e in range(len(boff)):
        w = 0.0
        for t in range(n5[e]):
            v = target5[boff[e] + t]
            acc = 0.0
            for _ in range(5):
                w = v + (w - v) * (k_wet if v > w else k_dry)
                acc += w
            out[boff[e] + t] = acc / 5.0
    return out


@njit(cache=True)
def _ar1(z, rho):
    """A stationary AR(1) series with unit variance from standard normal numbers (in place)."""
    f = np.sqrt(1.0 - rho * rho)
    for i in range(1, len(z)):
        z[i] = rho * z[i - 1] + f * z[i]
    return z


# =================================================================================================== sizes measured on real data
def noise_sizes(D, le):
    """The noise of the real TRAIN link-events on dry minutes (no rain model): the white noise of a 5 min value and the
    offset error of a link-event.
    white = the median over link-events of the sd of the 5 min means of (A - c) inside the offset window;
    offset = the sd over link-events of the mean of (A - c) over the dry bins of the 30 min before the first wet bin."""
    white, off = [], []
    for e in range(len(le)):
        o1, n, o5 = D["le_boff1"][e], D["le_n"][e], D["le_boff5"][e]
        y1 = np.where(D["w1"][o1:o1 + n] > 0, D["y1"][o1:o1 + n], np.nan) - D["le_c"][e]
        m = np.where(D["pre1"][o1:o1 + n], y1, np.nan).reshape(-1, 5)
        ok = np.isfinite(m).sum(axis=1) >= 3
        if ok.sum() >= 4:
            with np.errstate(all="ignore"):
                white.append(np.nanstd(np.nanmean(m[ok], axis=1), ddof=1))
        first = D["le_first"][e]
        sl = slice(o5 + max(first - 6, 0), o5 + first)
        dry = (~D["wet5"][sl]) & (D["w5"][sl] > 0)
        if dry.sum() >= 3:
            off.append(np.mean(D["y5"][sl][dry]) - D["le_c"][e])
    return {"white_sd_db": float(np.median(white)), "white_link_events": len(white), "offset_sd_db": float(np.std(off, ddof=1)), "offset_link_events": len(off)}


# =================================================================================================== the time axis shared by the links
def time_index(le):
    """(tb, ev_code): the absolute 5 min bin number of every bin of every link-event, and the event number per link-event."""
    t0 = le.window_start.min()
    start = ((le.window_start - t0) / pd.Timedelta("5min")).round().astype(int).to_numpy()
    n5 = le.n_bins.to_numpy()
    tb = np.concatenate([s + np.arange(n) for s, n in zip(start, n5)])
    return tb, pd.factorize(le.event_id)[0]


def area_mean(r5, le):
    """Z: per event and time the mean rain over all its link-events (the same value for every link at that time)."""
    tb, ev = time_index(le)
    ev_bin = np.repeat(ev, le.n_bins.to_numpy())
    key = ev_bin.astype(np.int64) * (tb.max() + 1) + tb
    uniq, inv = np.unique(key, return_inverse=True)
    fin = np.isfinite(r5)
    s = np.bincount(inv, weights=np.where(fin, r5, 0.0), minlength=len(uniq))
    c = np.bincount(inv, weights=fin.astype(float), minlength=len(uniq))
    return (s / np.maximum(c, 1.0))[inv]


# =================================================================================================== the rain error
def rain_error(le, n_cells, sigma, cfg, rng, scale_sd=None):
    """ln(R_obs / R_true) for every bin (E1 model): s per link-event + sqrt(q) sigma c_t (shared by all links of the event
    at that time) + sqrt((1 - q) / n_cells) sigma e (the link's own); c and e are AR(1); the mean of exp(error) is one."""
    R = cfg["rain_error"]
    q, rho = R["shared_share"], R["rho"]
    s_sd = R["scale_sd"] if scale_sd is None else scale_sd
    n5 = le.n_bins.to_numpy()
    rep = np.repeat(np.arange(len(le)), n5)
    tb, ev = time_index(le)
    ev_bin = ev[rep]
    c = np.empty(len(tb))
    for e_ in np.unique(ev):                                    # one shared series per event over its time span
        m = ev_bin == e_
        lo = tb[m].min()
        series = _ar1(rng.standard_normal(tb[m].max() - lo + 1), rho)
        c[m] = series[tb[m] - lo]
    own = np.empty(len(tb))
    z = rng.standard_normal(len(tb))
    boff = np.r_[0, np.cumsum(n5)[:-1]]
    for i in range(len(le)):
        own[boff[i]:boff[i] + n5[i]] = _ar1(z[boff[i]:boff[i] + n5[i]], rho)
    var_own = (1 - q) * sigma ** 2 / n_cells
    var_tot = s_sd ** 2 + q * sigma ** 2 + var_own
    s = rng.normal(0, s_sd, len(le)) if s_sd > 0 else np.zeros(len(le))
    return (s - 0.5 * var_tot)[rep] + np.sqrt(q) * sigma * c + np.sqrt(var_own)[rep] * own


# =================================================================================================== a fake dataset
def make(D, le, n_cells, sizes, sigma, cfg, seed, memory=True, with_noise=True, scale_sd=None):
    """A fake dataset on the link-events of D: known F and P, one store with memory, measured noise, the E1 rain error.
    Returns the dataset as the length test reads it (r5 = OBSERVED rain, y5 = fake A, c = 0) and keeps the true rain.
    sigma = the per-bin sd of the rain error; scale_sd = the scale part per link-event (default: the config value; with
    sigma = 0 there is no rain error at all)."""
    rng = np.random.default_rng(seed)
    T = cfg["truth"]
    n5 = D["le_n"] // 5
    rep = np.repeat(np.arange(len(le)), n5)
    r_true = D["r5"].astype(float)
    band = le.band_name.to_numpy()
    target = f_true(r_true, cfg)
    if memory:
        w = _store(target, D["le_boff5"], n5, np.exp(-1.0 / T["tau_wet_min"]), np.exp(-1.0 / T["tau_dry_min"]))
    else:
        w = target
    p = np.array([T["p_per_mmh"][b] for b in band])[rep] * np.nan_to_num(r_true)
    y = w + le.L_km.to_numpy(float)[rep] * p
    if with_noise:
        link_off = rng.normal(0, sizes["link_sd_db"], D["nl"])[D["le_link"]]
        y = y + sizes["white_sd_db"] * rng.standard_normal(len(y)) + (rng.normal(0, sizes["offset_sd_db"], len(le)) + link_off)[rep]
    err = rain_error(le, n_cells, sigma, cfg, rng, scale_sd=scale_sd) if sigma > 0 else np.zeros(len(y))
    Dk = dict(D)
    Dk["r_true"] = r_true
    Dk["r5"] = r_true * np.exp(err)
    Dk["y5"] = y
    Dk["le_c"] = np.zeros(len(le))
    Dk["w5"] = ((D["w5"] > 0) & np.isfinite(r_true)).astype(float)
    return Dk


# =================================================================================================== the judge
def judge(cell_tables, cfg, min_reps=None, l_med=None):
    """Per cell over the replicates: the median error of I against F_true and of s against P_true, both at the cell's own
    R_cell. A cell is judged if it is eligible in at least half of the replicates. Returns (table, n_fail, passed).
    l_med (km) given: the slope is judged by the absolute rule of Deviation 3, |s - P_true| <= intercept tolerance / l_med
    (the intercept tolerance applied to the path term of a median link). l_med = None: the first rule (relative), kept so
    that the E3 diagnostic can be reproduced."""
    rows = []
    for rep, t in enumerate(cell_tables):
        e = t[t.eligible]
        for r in e.itertuples():
            ft, pt = float(f_true(r.R_cell, cfg)), float(p_true(r.band, r.R_cell, cfg))
            rows.append({"rep": rep, "band": r.band, "class": r._2, "cls": r.cls, "R_cell": r.R_cell, "I": r.I, "s": r.s, "F_true": ft, "P_true": pt,
                         "I_err": r.I - ft, "s_err_abs": r.s - pt, "s_err_rel": (r.s - pt) / pt, "links_short": r.links_short})
    if not rows:
        return pd.DataFrame(), 0, False
    a = pd.DataFrame(rows)
    g = a.groupby(["band", "cls", "class"], sort=True)
    out = g.agg(reps=("rep", "size"), R_cell=("R_cell", "median"), F_true=("F_true", "median"), I=("I", "median"), I_err=("I_err", "median"),
                P_true=("P_true", "median"), s=("s", "median"), s_err_abs=("s_err_abs", "median"), s_err_rel=("s_err_rel", "median"), links_short=("links_short", "median")).reset_index()
    out = out[out.reps >= (min_reps or int(np.ceil(len(cell_tables) / 2)))].reset_index(drop=True)
    out["I_ok"] = out.I_err.abs() <= cfg["pass"]["intercept_db"]
    out["s_ok"] = (out.s_err_rel.abs() <= cfg["pass"]["slope_rel"]) if l_med is None else (out.s_err_abs.abs() <= cfg["pass"]["intercept_db"] / l_med)
    n_fail = int((~(out.I_ok & out.s_ok)).sum())
    return out, n_fail, bool(len(out) > 0 and n_fail <= cfg["pass"]["cells_may_fail"])
