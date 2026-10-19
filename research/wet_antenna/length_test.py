"""Library of steps S3 and S3b: the model-free length test and the helpers to run a fitted model on another dataset.

The length test (S3.1): in steady rain, per link-event and rain class the median of (A - c); then per band and class a
robust line over all link lengths

    (A - c) = I + s * L        I = the part that does not grow with L (the antenna), s = the part per km (the path)

with an event bootstrap. Nothing of a model enters the bins or the line: the same line can then be fitted to a model's
own prediction (A_pred - c) on the same bins ("the model through the test") and compared with the observed I and s.
`lt` = config section `length_test` of step S3 (classes, steady-rain rule, Huber settings, bootstrap, verdict).
"""
import pickle

import numpy as np
import pandas as pd

from research.wet_antenna import paths as PATHS
from research.wet_antenna import system_id as sid


def class_names(lt):
    return [f"{lo:g}-{hi:g}" for lo, hi in lt["rain_classes_mmh"]]


# =================================================================================================== a fitted model on another dataset
def link_ids(D, le):
    """cml_id of every link index of a dataset."""
    ids = np.empty(D["nl"], object)
    ids[D["le_link"]] = le.cml_id.to_numpy()
    return ids


def s2b_order():
    """The cml_id of every link index of the S2b datasets S<3 ("S3") and S<1 ("S1") (from the S2b product)."""
    prod = pd.read_csv(PATHS.SID2B_LINK_EVENTS, usecols=["cml_id", "in_S1"])
    return {"S3": list(pd.unique(prod.cml_id)), "S1": list(pd.unique(prod.cml_id[prod.in_S1]))}


def s2b_fit(set_, model, order):
    """One S2b fit (on all its TRAIN events, read only): x, lag d, a_link by cml_id, the readable parameters, its
    leave-one-event-out prediction. `order` = s2b_order()."""
    res = pickle.loads((PATHS.SID2B_CACHE / f"real_{set_}_{model}.pkl").read_bytes())
    assert len(order[set_]) == len(res["a_link"])
    return {"x": res["all"]["x"], "d": int(res["all"]["d"]), "a": dict(zip(order[set_], res["a_link"])), "a_med": float(np.median(res["a_link"])),
            "params": res["params"], "loo_pred": res.get("loo_pred"), "model": model, "set": set_}


def x_for(Pm, fz, ids):
    """The parameter vector of a fit (fz: x, a = {cml_id: a_link}) on a new dataset (Problem2 Pm, its links `ids`): every
    curve and time constant as fitted, a_link of a link from the fit if the link was in the fit's set, else 1.
    Returns (x, number of links with a_link = 1)."""
    is_a = Pm.names == "loga"
    first = int(np.argmax(is_a))
    nl_old = len(fz["a"])
    xo = fz["x"]
    assert np.allclose(np.exp(xo[first:first + nl_old]), list(fz["a"].values()))
    x = np.empty(len(Pm.names))
    x[~is_a] = np.r_[xo[:first], xo[first + nl_old:]]
    x[is_a] = [np.log(fz["a"][c]) if c in fz["a"] else 0.0 for c in ids]
    return x, int(sum(c not in fz["a"] for c in ids))


def curve(params, kind, band, r, knots):
    """F or P of a fit at rain r (linear between the knots, 0 at R = 0, constant above the last knot)."""
    v = np.r_[0.0, [params.get(f"{kind}[{band}]({k:g})", np.nan) for k in knots[1:]]]
    return float(np.interp(r, knots, v)) if np.isfinite(v).all() else np.nan


# =================================================================================================== the length test
def steady_table(D, le, lt, preds=None, r15=False, late=False):
    """The steady-rain 5 min bins of every link-event, reduced to one row per link-event and rain class: the median of
    (A - c), of R, and (if preds = {name: 5 min prediction}) of (A_pred - c) of every prediction on the same bins.
    r15: the radar R replaced by its centred 15 min mean; late: only bins whose previous robustness.late_bins_min are wet."""
    classes = [tuple(c) for c in lt["rain_classes_mmh"]]
    names = class_names(lt)
    thr, (q_lo, q_hi) = lt["wet_threshold_mmh"], lt["prev_ratio"]
    k_wet = (lt["robustness"]["late_bins_min"] if late else lt["wet_before_min"]) // 5
    k_prev = lt["prev_window_min"] // 5
    edges = np.array([c[0] for c in classes] + [classes[-1][1]], float)
    y_all = np.where(D["w5"] > 0, D["y5"], np.nan)
    preds = preds or {}
    rows = []
    for e in range(len(le)):
        o, n = D["le_boff5"][e], D["le_n"][e] // 5
        r = D["r5"][o:o + n].astype(float)
        if r15:
            rr = np.full(n, np.nan)
            rr[1:-1] = (r[:-2] + r[1:-1] + r[2:]) / 3
            r = rr
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
        cls = np.searchsorted(edges, r[ok], side="right") - 1
        yk, rk = y[ok], r[ok]
        pk = {k: (p[o:o + n] - D["le_c"][e])[ok] for k, p in preds.items()}
        for c in np.unique(cls):
            m = cls == c
            rows.append({"le_idx": e, "cls": int(c), "A_med": float(np.median(yk[m])), "R_med": float(np.median(rk[m])), "n_bins": int(m.sum()),
                         **{f"pred_{k}": float(np.median(v[m])) for k, v in pk.items()}})
    t = pd.DataFrame(rows)
    info = le[["event_id", "sublink_id", "cml_id", "L_km", "band_name", "frequency_ghz"]].reset_index(drop=True)
    t = t.join(info, on="le_idx")
    t["ev"] = D["le_event"][t.le_idx.to_numpy()]
    t["class"] = [names[c] for c in t.cls]
    return t


def wmedian(v, w):
    o = np.argsort(v)
    cw = np.cumsum(w[o])
    return float(v[o][np.searchsorted(cw, 0.5 * cw[-1])])


def huber_fit(X, y, w, lt):
    """Huber regression by iteratively re-weighted least squares: k = lt.huber_k, the scale = the weighted MAD of the
    residuals / 0.6745, re-estimated every iteration. w = multiplicities (event bootstrap)."""
    keep = w > 0
    X, y, w = X[keep], y[keep], w[keep]
    if len(np.unique(X[:, 1])) < 3:
        return np.full(X.shape[1], np.nan)
    sw = np.sqrt(w)
    beta = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
    for _ in range(lt["irls_iterations"]):
        r = y - X @ beta
        s = wmedian(np.abs(r), w) / 0.6745
        if s <= 1e-9:
            break
        u = np.abs(r) / s
        ww = w * np.where(u <= lt["huber_k"], 1.0, lt["huber_k"] / np.maximum(u, 1e-12))
        sw = np.sqrt(ww)
        new = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
        if np.allclose(new, beta, rtol=1e-7, atol=1e-9):
            beta = new
            break
        beta = new
    return beta


def length_fits(t, n_ev, variant, lt, ycol="A_med", boot=True):
    """Per band and class: I + s L (and I + s L + q L^2) with event-bootstrap intervals; L0 = I / s."""
    names, vd = class_names(lt), lt["verdict"]
    rng = np.random.default_rng(lt["seed"])
    mult = [np.bincount(rng.integers(0, n_ev, n_ev), minlength=n_ev).astype(float) for _ in range(lt["boot_reps"] if boot else 0)]
    rows = []
    for (band, c), g in t.groupby(["band_name", "cls"]):
        L, y, ev = g.L_km.to_numpy(float), g[ycol].to_numpy(float), g.ev.to_numpy()
        row = {"variant": variant, "band": band, "class": names[c], "cls": c, "link_events": len(g), "links": g.cml_id.nunique(),
               "events": g.event_id.nunique(), "bins": int(g.n_bins.sum()), "links_under_1km": g[g.L_km < vd["short_km"]].cml_id.nunique(),
               "L_max_km": float(L.max()), "R_cell": float(np.median(g.R_med))}
        if len(g) < lt["min_link_events"] or len(np.unique(L)) < 3:
            rows.append(row)
            continue
        X1, X2 = np.c_[np.ones_like(L), L], np.c_[np.ones_like(L), L, L ** 2]
        b1, b2 = huber_fit(X1, y, np.ones_like(L), lt), huber_fit(X2, y, np.ones_like(L), lt)
        row.update({"I": b1[0], "s": b1[1], "L0": b1[0] / b1[1] if b1[1] > 0 else np.nan, "I_q": b2[0], "s_q": b2[1], "q": b2[2]})
        if boot:
            B1 = np.array([huber_fit(X1, y, m[ev], lt) for m in mult])
            B2 = np.array([huber_fit(X2, y, m[ev], lt) for m in mult])
            l0 = np.where(B1[:, 1] > 0, B1[:, 0] / np.where(B1[:, 1] > 0, B1[:, 1], 1.0), np.nan)
            for name, v in (("I", B1[:, 0]), ("s", B1[:, 1]), ("L0", l0), ("I_q", B2[:, 0]), ("s_q", B2[:, 1]), ("q", B2[:, 2])):
                row[f"{name}_lo"], row[f"{name}_hi"] = np.nanpercentile(v, lt["interval"]) if np.isfinite(v).any() else (np.nan, np.nan)
            row["L0_draws_s_le_0"] = int((B1[:, 1] <= 0).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def verdict_table(cmp_, share_needed, cols=("F_inside", "P_inside")):
    """Per variant and fit: how many eligible cells have both quantities inside, and the verdict."""
    ver = []
    a, b = cols
    for (variant, key), g in cmp_.groupby(["variant", "fit"], sort=False):
        e = g[g.eligible]
        both = e[a] & e[b]
        share = both.mean() if len(e) else np.nan
        ver.append({"variant": variant, "fit": key, "eligible_cells": len(e), a: int(e[a].sum()), b: int(e[b].sum()), "both_inside": int(both.sum()),
                    "share": share, "verdict": ("consistent" if share >= share_needed else "NOT consistent") if len(e) else "no eligible cell"})
    return pd.DataFrame(ver)


def compare_models(fits, fz_all, lt, knots, through=None):
    """Per fitted cell and fit: a F and P at R_cell (equilibrium) against the intervals of I and s; eligibility; the verdict.
    through (optional) = the length fits of the models' own predictions: their I and s are added to each row."""
    vd = lt["verdict"]
    rows = []
    for r in fits.to_dict("records"):
        if not np.isfinite(r.get("I", np.nan)):
            continue
        for key, fz in fz_all.items():
            F = fz["a_med"] * curve(fz["params"], "F", r["band"], r["R_cell"], knots)
            Pv = curve(fz["params"], "P", r["band"], r["R_cell"], knots)
            has = np.isfinite(F) and np.isfinite(Pv)
            row = {"variant": r["variant"], "band": r["band"], "class": r["class"], "fit": key, "R_cell": r["R_cell"], "links_under_1km": r["links_under_1km"],
                   "I": r["I"], "I_lo": r["I_lo"], "I_hi": r["I_hi"], "aF_model": F, "s": r["s"], "s_lo": r["s_lo"], "s_hi": r["s_hi"], "P_model": Pv,
                   "eligible": bool(has and r["links_under_1km"] >= vd["min_short_links"]),
                   "F_inside": bool(has and r["I_lo"] <= F <= r["I_hi"]), "P_inside": bool(has and r["s_lo"] <= Pv <= r["s_hi"])}
            row["both_inside"] = row["F_inside"] and row["P_inside"]
            if through is not None:
                m = through[(through.band == r["band"]) & (through["class"] == r["class"]) & (through.fit == key)]
                row["I_through_model"] = float(m.I.iloc[0]) if (len(m) and has) else np.nan
                row["s_through_model"] = float(m.s.iloc[0]) if (len(m) and has) else np.nan
            rows.append(row)
    cmp_ = pd.DataFrame(rows)
    return cmp_, verdict_table(cmp_, vd["share"])


# =================================================================================================== error set of a prediction
def val_metrics(D_, pred, J):
    """The S3.2 error set of one prediction on one dataset (event-bootstrap 90 % intervals). J = config `judge` of
    the step (moved unchanged from system_id_s3b.py; used by S3b and S5a)."""
    n_ev = int(D_["le_event"].max()) + 1
    err = sid.event_errors(D_, pred, J)
    tab = sid.judge({"x": err}, ["x"], J)[0].iloc[0]
    sh = sid.shower_table(D_, J["shower_gap_min"])
    pk = sid.shower_stat(*sid.shower_peaks(D_, sh, pred, J["shower_min_peak_db"]), n_ev, J)
    draws = np.random.default_rng(J["seed"]).integers(0, n_ev, (J["boot_reps"], n_ev))
    io, ip = err.integral_obs.to_numpy(), err.integral_pred.to_numpy()
    rb = np.nanmean(ip[draws], axis=1) / np.nanmean(io[draws], axis=1)
    return {"events": n_ev, "rmse": tab.rmse, "rmse_lo": tab.rmse_lo, "rmse_hi": tab.rmse_hi, "rmse_rain": tab.rmse_rain, "rmse_rain_lo": tab.rmse_rain_lo,
            "rmse_rain_hi": tab.rmse_rain_hi, "rmse_0_60": tab.rmse_0_60, "rmse_60_180": tab.rmse_60_180, "shower_peak_error": pk[0], "peak_lo": pk[1], "peak_hi": pk[2],
            "showers": pk[3], "integral_pred_over_obs": float(np.nanmean(ip) / np.nanmean(io)), "integral_lo": float(np.nanpercentile(rb, 5)),
            "integral_hi": float(np.nanpercentile(rb, 95))}
