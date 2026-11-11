"""SUPERSEDED (Detachment correction, 2026-11-11): this library used the Part 1 steady-rain rule, rain classes, eligibility rule
and Huber code. It is kept only so that the superseded A3 outputs can be reproduced. The current library is `ltest.py`.

Library of step A3: the model-free length test with options, and the curves F(R), P(R) built from it.

The steady-rain rule and the Huber line are those of Part 1 (`research.wet_antenna.length_test`, read only). This library
adds what Part 1 does not have: the same-class duration of every steady bin, a bootstrap that keeps its draws (so that two
variants can be compared draw by draw, and curves can be built per draw), and the monotone curves.

    bins = steady_bins(D, le, classes)                    one row per steady 5 min bin
    t = reduce(bins, D, le, min_run=1)                    one row per link-event and class (medians)
    cells, draws = cell_fits(t, mult, classes)            the line (A - c) = I + s L per band and class, with the draws
"""
import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

from research.wet_antenna import length_test as LTL
from research.waa_ksurface import data
from research.waa_ksurface import paths as P

LT = data.LT                                                   # the Part 1 length test settings (read only)


def class_names(classes):
    return [f"{lo:g}-{hi:g}" for lo, hi in classes]


def steady_bins(D, le, classes):
    """Every steady-rain 5 min bin (the Part 1 rule of `length_test.steady_table`, unchanged), one row per bin, with
    `run` = the number of consecutive bins ending at this bin whose R is in the same rain class (5 * run = its duration)."""
    thr, (q_lo, q_hi) = LT["wet_threshold_mmh"], LT["prev_ratio"]
    k_wet, k_prev = LT["wet_before_min"] // 5, LT["prev_window_min"] // 5
    edges = np.array([c[0] for c in classes] + [classes[-1][1]], float)
    y_all = np.where(D["w5"] > 0, D["y5"], np.nan)
    parts = []
    for e in range(len(le)):
        o, n = D["le_boff5"][e], D["le_n"][e] // 5
        r = D["r5"][o:o + n].astype(float)
        fin = np.isfinite(r)
        wet = np.where(fin, r, -1.0) >= thr
        t = np.arange(n)
        cw, cf, cr = np.r_[0, np.cumsum(wet)], np.r_[0, np.cumsum(fin)], np.r_[0, np.cumsum(np.where(fin, r, 0.0))]
        t_w, t_p = np.maximum(t - k_wet, 0), np.maximum(t - k_prev, 0)
        ok_wet = (t >= k_wet) & (cw[t] - cw[t_w] == k_wet)
        prev = (cr[t] - cr[t_p]) / k_prev
        ok_prev = (t >= k_prev) & (cf[t] - cf[t_p] == k_prev) & (prev >= q_lo * r) & (prev <= q_hi * r)
        y = y_all[o:o + n] - D["le_c"][e]
        ok = ok_wet & ok_prev & np.isfinite(y) & fin & (r >= edges[0]) & (r < edges[-1])
        if not ok.any():
            continue
        cls_all = np.where(fin & (r >= edges[0]) & (r < edges[-1]), np.searchsorted(edges, np.where(fin, r, 0.0), side="right") - 1, -1)
        run = np.ones(n, int)                                   # the length of the run of equal class ending at each bin
        for i in range(1, n):
            if cls_all[i] >= 0 and cls_all[i] == cls_all[i - 1]:
                run[i] = run[i - 1] + 1
        parts.append(pd.DataFrame({"le_idx": e, "bin": t[ok], "cls": cls_all[ok], "y": y[ok], "r": r[ok], "run": run[ok]}))
    return pd.concat(parts, ignore_index=True)


def reduce(bins, D, le, classes, min_run=1):
    """One row per link-event and class: the median of (A - c) and of R over its steady bins with run >= min_run."""
    b = bins[bins.run >= min_run]
    t = b.groupby(["le_idx", "cls"]).agg(A_med=("y", "median"), R_med=("r", "median"), n_bins=("y", "size")).reset_index()
    info = le[["event_id", "sublink_id", "cml_id", "L_km", "band_name"]].reset_index(drop=True)
    t = t.join(info, on="le_idx")
    t["ev"] = D["le_event"][t.le_idx.to_numpy()]
    t["class"] = np.array(class_names(classes))[t.cls.to_numpy()]
    return t


def multiplicities(n_ev, reps, seed):
    """The event bootstrap: per draw, how many times every event is drawn (events resampled with replacement)."""
    rng = np.random.default_rng(seed)
    return [np.bincount(rng.integers(0, n_ev, n_ev), minlength=n_ev).astype(float) for _ in range(reps)]


def cell_fits(t, mult, classes, interval=(5, 95), ycol="A_med"):
    """Per band and class the Huber line (A - c) = I + s L (Part 1 `huber_fit`). Returns (table, draws): the table has I, s,
    their bootstrap intervals, the counts and `eligible` (the Part 1 rule); draws[(band, class)] = array (reps, 2) of (I, s)."""
    names, vd = class_names(classes), LT["verdict"]
    rows, draws = [], {}
    for (band, c), g in t.groupby(["band_name", "cls"]):
        L, y, ev = g.L_km.to_numpy(float), g[ycol].to_numpy(float), g.ev.to_numpy()
        row = {"band": band, "class": names[c], "cls": int(c), "link_events": len(g), "links": g.cml_id.nunique(), "events": g.event_id.nunique(),
               "bins": int(g.n_bins.sum()), "links_under_1km": g[g.L_km < vd["short_km"]].cml_id.nunique(), "R_cell": float(np.median(g.R_med)),
               "I": np.nan, "s": np.nan}
        fitted = len(g) >= LT["min_link_events"] and len(np.unique(L)) >= 3
        if fitted:
            X = np.c_[np.ones_like(L), L]
            b = LTL.huber_fit(X, y, np.ones_like(L), LT)
            row["I"], row["s"] = float(b[0]), float(b[1])
            if mult:
                B = np.array([LTL.huber_fit(X, y, m[ev], LT) for m in mult])
                draws[(band, names[c])] = B
                for k, name in enumerate(("I", "s")):
                    row[f"{name}_lo"], row[f"{name}_hi"] = np.nanpercentile(B[:, k], interval) if np.isfinite(B[:, k]).any() else (np.nan, np.nan)
                    row[f"{name}_se"] = float(np.nanstd(B[:, k]))
        row["eligible"] = bool(fitted and np.isfinite(row["I"]) and row["links_under_1km"] >= vd["min_short_links"])
        rows.append(row)
    return pd.DataFrame(rows), draws


# =================================================================================================== the curves
def pava(v, w):
    """The closest non-decreasing sequence to v (weighted least squares): pool adjacent violators."""
    v, w = list(map(float, v)), list(map(float, w))
    blocks = [[v[i], w[i], 1] for i in range(len(v))]          # value, weight, number of points
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
    (PCHIP) between the cells. Above the highest cell: F is held constant, P continues with the slope of its last segment."""
    r_cells = np.asarray(r_cells, float)
    v = np.maximum(pava(values, weights), 0.0)
    x, yv = np.r_[0.0, r_cells], np.r_[0.0, v]
    f = PchipInterpolator(x, yv)
    last_slope = (yv[-1] - yv[-2]) / (x[-1] - x[-2])

    def curve(r):
        r = np.asarray(r, float)
        inside = f(np.clip(r, 0.0, x[-1]))
        return inside + (np.maximum(r - x[-1], 0.0) * last_slope if kind == "P" else 0.0)
    return curve, v


class Curves:
    """F(R) and P(R) per band, read from the product table of step A3 (a dense R grid, linear between the grid points).

        C = Curves()                      # the central curves
        C.F("28 GHz", r), C.P("28 GHz", r)
        Curves(which="lo") / Curves(which="hi") / Curves(which="gauge")    # the sensitivity curves of rule A3f (F only)
    """

    def __init__(self, which="mid", table=None):
        t = table if table is not None else pd.read_csv(P.A3_CURVES)
        self.col = {"mid": "F", "lo": "F_lo", "hi": "F_hi", "gauge": "F_gauge_events"}[which]
        self.t = {b: g.sort_values("R") for b, g in t.groupby("band")}
        self.bands = list(self.t)

    def F(self, band, r):
        g = self.t[band]
        return np.interp(np.asarray(r, float), g.R.to_numpy(), g[self.col].to_numpy())

    def P(self, band, r):
        g = self.t[band]
        r = np.asarray(r, float)
        R, Pv = g.R.to_numpy(), g.P.to_numpy()
        top = (Pv[-1] - Pv[-2]) / (R[-1] - R[-2])              # beyond the grid the last slope continues
        return np.interp(r, R, Pv) + np.maximum(r - R[-1], 0.0) * top

    def valid_range(self, band):
        g = self.t[band]
        m = g[g.measured]
        return float(m.R.min()), float(m.R.max())
