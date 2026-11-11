"""Library of the A3 rescue (steps E2-E4): the length test of this branch, with its own definitions.

Nothing here comes from Part 1. Every definition is derived in `docs/waa_ksurface/decisions.md` (step E3) and its numbers
are in `research/waa_ksurface/config.yaml` (section `ltest`).

The idea: in steady rain the attenuation of a link is (A - c) = W + L * P. W (the wet antenna) does not grow with the link
length L, the path loss does. So over the links of one frequency band, in one rain class, the line

    (A - c) = I + s * L

has the intercept I = the steady WAA F(R) and the slope s = the path loss per km P(R).

    t = steady_values(D, le, cfg)            one row per link-event and rain class (medians over its steady bins)
    cells, draws = cell_fits(t, cfg, ...)    the line per band and class, with event-bootstrap draws
    curve = monotone_curve(...)              a rising curve through (0, 0) and the cell values
"""
import numpy as np
import pandas as pd
from numba import njit
from scipy.interpolate import PchipInterpolator

from research.waa_ksurface import data


def config():
    """The length-test settings of this branch (config section `ltest`)."""
    return data.branch_config()["ltest"]


def class_names(cfg):
    e = cfg["class_edges_mmh"]
    return [f"{lo:g}-{hi:g}" for lo, hi in zip(e[:-1], e[1:])]


# =================================================================================================== steady rain
@njit(cache=True)
def _steady(sel, own_ok, boff, n5, edges, k_pre, q):
    """The rain class of every steady bin, -1 for the other bins. A bin is steady if its selecting rain `sel` is inside a
    class and each of the k_pre bins before it lies within [sel / q, sel * q] (a missing value fails)."""
    out = np.full(len(sel), -1)
    nc = len(edges) - 1
    for e in range(len(boff)):
        o, n = boff[e], n5[e]
        for t in range(k_pre, n):
            r = sel[o + t]
            if not (r >= edges[0] and r < edges[nc]) or not own_ok[o + t]:
                continue
            ok = True
            for j in range(1, k_pre + 1):
                rp = sel[o + t - j]
                if not (rp >= r / q and rp <= r * q):
                    ok = False
                    break
            if ok:
                c = 0
                while r >= edges[c + 1]:
                    c += 1
                out[o + t] = c
    return out


def steady_values(D, le, cfg, rule, sel=None, r_stat="median"):
    """One row per link-event and rain class: the median of (A - c) and the median (or mean) of the link's own R over its
    steady bins. rule = {"t_pre_min", "q"}. sel = the rain that selects the bins (default: the link's own rain, D["r5"])."""
    edges = np.asarray(cfg["class_edges_mmh"], float)
    own = D["r5"].astype(float)
    own_ok = (D["w5"] > 0) & np.isfinite(own)
    cls = _steady(own if sel is None else np.asarray(sel, float), own_ok, D["le_boff5"], D["le_n"] // 5, edges, int(rule["t_pre_min"]) // 5, float(rule["q"]))
    keep = cls >= 0
    le_idx = np.repeat(np.arange(len(le)), D["le_n"] // 5)[keep]
    b = pd.DataFrame({"le_idx": le_idx, "cls": cls[keep], "y": (D["y5"] - D["le_c"][np.repeat(np.arange(len(le)), D["le_n"] // 5)])[keep], "r": own[keep]})
    t = b.groupby(["le_idx", "cls"]).agg(A_med=("y", "median"), R_med=("r", r_stat), n_bins=("y", "size")).reset_index()
    info = le[["event_id", "cml_id", "L_km", "band_name"]].reset_index(drop=True)
    t = t.join(info, on="le_idx")
    t["ev"] = D["le_event"][t.le_idx.to_numpy()]
    return t


# =================================================================================================== the line
def huber_line(L, y, w, k=1.345, iterations=50):
    """Robust line y = I + s L (Huber, iteratively re-weighted least squares; scale = weighted MAD / 0.6745).
    w = how many times each value counts (the event bootstrap). Returns (I, s), NaN if there are fewer than 3 lengths."""
    keep = w > 0
    L, y, w = L[keep], y[keep], w[keep]
    if len(np.unique(L)) < 3:
        return np.array([np.nan, np.nan])
    X = np.c_[np.ones_like(L), L]
    ww = w.copy()
    beta = np.zeros(2)
    for i in range(iterations + 1):
        sw = np.sqrt(ww)
        new = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
        if i and np.allclose(new, beta, rtol=1e-7, atol=1e-9):
            return new
        beta = new
        res = np.abs(y - X @ beta)
        o = np.argsort(res)
        cw = np.cumsum(w[o])
        scale = res[o][np.searchsorted(cw, 0.5 * cw[-1])] / 0.6745          # weighted median of |residual|
        if scale <= 1e-9:
            return beta
        u = res / scale
        ww = w * np.where(u <= k, 1.0, k / np.maximum(u, 1e-12))
    return beta


def multiplicities(n_ev, reps, seed):
    """The event bootstrap: per draw, how many times every event is drawn (events resampled with replacement)."""
    rng = np.random.default_rng(seed)
    return [np.bincount(rng.integers(0, n_ev, n_ev), minlength=n_ev).astype(float) for _ in range(reps)]


def cell_fits(t, cfg, n_short, mult=None, interval=(5, 95), r_stat="median"):
    """Per band and rain class the line (A - c) = I + s L. Returns (table, draws).
    fitted: values from at least cfg.min_events events. eligible: fitted, and at least n_short distinct links below
    cfg.short_km and at least n_short distinct links of cfg.short_km or more. draws[(band, class)] = array (reps, 2)."""
    names = class_names(cfg)
    rows, draws = [], {}
    for (band, c), g in t.groupby(["band_name", "cls"]):
        L, y, ev = g.L_km.to_numpy(float), g.A_med.to_numpy(float), g.ev.to_numpy()
        short = g[g.L_km < cfg["short_km"]].cml_id.nunique()
        longer = g[g.L_km >= cfg["short_km"]].cml_id.nunique()
        row = {"band": band, "class": names[c], "cls": int(c), "link_events": len(g), "links": g.cml_id.nunique(), "events": g.event_id.nunique(),
               "bins": int(g.n_bins.sum()), "links_short": short, "links_long": longer,
               "R_cell": float(np.median(g.R_med) if r_stat == "median" else np.mean(g.R_med)), "I": np.nan, "s": np.nan}
        row["fitted"] = bool(row["events"] >= cfg["min_events"] and len(np.unique(L)) >= 3)
        if row["fitted"]:
            b = huber_line(L, y, np.ones_like(L), cfg["huber_k"])
            row["I"], row["s"] = float(b[0]), float(b[1])
            if mult:
                B = np.array([huber_line(L, y, m[ev], cfg["huber_k"]) for m in mult])
                draws[(band, names[c])] = B
                for k, name in enumerate(("I", "s")):
                    row[f"{name}_lo"], row[f"{name}_hi"] = np.nanpercentile(B[:, k], interval) if np.isfinite(B[:, k]).any() else (np.nan, np.nan)
                    row[f"{name}_se"] = float(np.nanstd(B[:, k]))
        row["eligible"] = bool(row["fitted"] and np.isfinite(row["I"]) and short >= n_short and longer >= n_short)
        rows.append(row)
    return pd.DataFrame(rows), draws


def residuals(t, cells):
    """The residual of every link-event value from its cell line (cells with a fitted line only)."""
    c = cells.set_index(["band", "cls"])
    key = list(zip(t.band_name, t.cls))
    I = np.array([c.I.get(k, np.nan) for k in key])
    s = np.array([c.s.get(k, np.nan) for k in key])
    return t.assign(resid=t.A_med - I - s * t.L_km)


def variance_components(v, group):
    """One-way random effects by the method of moments: (between-group sd, within-group sd, groups, values)."""
    df = pd.DataFrame({"v": v, "g": group}).dropna()
    n_i = df.groupby("g").v.size()
    k, N = len(n_i), len(df)
    if k < 2 or N <= k:
        return np.nan, np.nan, k, N
    means = df.groupby("g").v.mean()
    ss_w = float(((df.v - df.g.map(means)) ** 2).sum())
    ss_b = float((n_i * (means - df.v.mean()) ** 2).sum())
    ms_w, ms_b = ss_w / (N - k), ss_b / (k - 1)
    n0 = (N - (n_i ** 2).sum() / N) / (k - 1)
    return float(np.sqrt(max((ms_b - ms_w) / n0, 0.0))), float(np.sqrt(ms_w)), k, N


# =================================================================================================== the curves
def pava(v, w):
    """The closest non-decreasing sequence to v (weighted least squares): pool adjacent violators."""
    blocks = [[float(v[i]), float(w[i]), 1] for i in range(len(v))]          # value, weight, number of points
    i = 0
    while i < len(blocks) - 1:
        if blocks[i][0] > blocks[i + 1][0] + 1e-12:
            a, b = blocks[i], blocks[i + 1]
            blocks[i:i + 2] = [[(a[0] * a[1] + b[0] * b[1]) / (a[1] + b[1]), a[1] + b[1], a[2] + b[2]]]
            i = max(i - 1, 0)
        else:
            i += 1
    return np.concatenate([[b[0]] * b[2] for b in blocks])


def monotone_curve(r_cells, values, weights, kind):
    """A function of R through (0, 0) and the cell values: made non-decreasing (pava), never negative, monotone cubic
    (PCHIP) between the cells. Above the highest cell F is held constant and P continues with the slope of its last segment."""
    x = np.r_[0.0, np.asarray(r_cells, float)]
    v = np.r_[0.0, np.maximum(pava(values, weights), 0.0)]
    f = PchipInterpolator(x, v)
    last_slope = (v[-1] - v[-2]) / (x[-1] - x[-2])

    def curve(r):
        r = np.asarray(r, float)
        return f(np.clip(r, 0.0, x[-1])) + (np.maximum(r - x[-1], 0.0) * last_slope if kind == "P" else 0.0)
    return curve, v[1:]
