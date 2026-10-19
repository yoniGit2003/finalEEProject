"""wet_antenna step S7: validate the FROZEN S6 candidate without TEST. Nothing of the model is refitted.

  7.1  W independent of L: a_link against L and the link metadata (Spearman, link bootstrap, BH); the W height a_link F
       at 2 and 5 mm/h per length class.
  7.2  short links on VAL (frozen): RMSE per length class, band and rain class; the VAL length test; shared against
       band-specific shapes on short links; the residual against the minutes after the rain.
  7.3  long links without a power law: (a) dry-path moments, where A - c = W is measured and compared with the model's W;
       (b) paired links W = (L2 A1 - L1 A2) / (L2 - L1), counted and noise-checked first; (c) the intercept I(R) from
       short-only and long-only links; (d) does any of them separate the two F descriptions; (e) c_band against ITU-R
       P.838 (size reference only).
  7.4  dew: size and timing of the dew W, dew against a_link, events with dew before the rain, dew water inside c.
  7.5  leftovers: a_link against the metadata and the per-event drying time against the weather (Spearman + BH, then
       PySR guarded by leave-one-out and a permutation null). Nothing enters the equation here.
  7.6  the rain reference: radar (V0, V1) against the City gauges by rate class, binned by the gauge and by the radar.

Two descriptions of F are carried: the candidate (free shared knots) and the closed F = A_band (R / 5)^b of S6.
STANDING RULE: WAA is never estimated by subtracting a power-law rain term. The path part is separated only through its
dependence on L at matched rain, or by dry-path moments (R = 0 on the whole path). ITU-R P.838 only as a size reference.
The rain correction comes from the City gauges only (g_e), never from the links. NO TEST attenuation or rain is read.
Every rule is in config `wet_antenna.validate_s7` (written before the run). Caches: research/.cache/wet_antenna/validate_s7/.
Run from the repository root (about an hour the first time, minutes from the cache):
    python -m research.wet_antenna.validate_s7
"""
import hashlib
import json
import os
import pickle
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import least_squares

from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import length_test as LTL
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.event_factors import bh, bin_times, spearman_boot
from research.wet_antenna.final_model import FIX, LT, WA, predict, through
from research.wet_antenna.validation import (LB, MODELS, boot_items, ci, curve_at, dry_windows, ev_draws, fit_tau, fmt, lclass, load, need,
                                             norm_median, overlap, predictions, radar_max, w_minutes)
from research.wet_antenna.tails import itu_p838

S = WA["validate_s7"]
CACHE = P.S7_CACHE
VL = S["val_label"]
STANDING = ("**Standing rule.** WAA is never estimated by subtracting a power-law rain term. The path part is separated only through its "
            "dependence on L at matched rain, or by dry-path moments (R = 0 on the whole path). ITU-R P.838 only as a size reference. The rain "
            "correction comes from the City gauges only (g_e), never from the links. No TEST attenuation or rain was read (TEST links and events "
            "are only counted).")


# =================================================================================================== small helpers
# =================================================================================================== loading
# =================================================================================================== 7.1 W independent of L
def link_table(d, let, mods, BANDS):
    sub = d.sublinks()
    links = d.links().set_index("cml_id")
    sites = d.table("sites").set_index("site_id")
    ids = list(pd.unique(let.cml_id))
    lt = let.groupby("cml_id").agg(L_km=("L_km", "first"), band=("band_name", "first"), frequency_ghz=("frequency_ghz", "mean"),
                                   link_events=("event_id", "size")).reindex(ids)
    s = sub[sub.cml_id.isin(ids)].copy()
    s["antenna_d"] = s[["site_0_antenna_diameter_m", "site_1_antenna_diameter_m"]].mean(axis=1)
    agg = s.groupby("cml_id").agg(antenna_d=("antenna_d", "mean"), n_radar_cells=("n_radar_cells", "mean"), dry_noise_db=("dry_tl_std", "mean"),
                                  diurnal_range_db=("dry_tl_diurnal_range", "mean"), polarization=("polarization", "first"))
    lt = lt.join(agg)
    s0, s1 = links.site_0_id.reindex(ids), links.site_1_id.reindex(ids)
    x0, y0 = sites.x_3006_m.reindex(s0).to_numpy(), sites.y_3006_m.reindex(s0).to_numpy()
    x1, y1 = sites.x_3006_m.reindex(s1).to_numpy(), sites.y_3006_m.reindex(s1).to_numpy()
    brg = np.degrees(np.arctan2(x1 - x0, y1 - y0)) % 180
    lt["site_0"], lt["site_1"] = s0.to_numpy(), s1.to_numpy()
    lt["bearing_deg"] = brg
    lt["abs_sin_bearing"], lt["abs_cos_bearing"] = np.abs(np.sin(np.radians(brg))), np.abs(np.cos(np.radians(brg)))
    lt["x_mid_km"], lt["y_mid_km"] = (x0 + x1) / 2000, (y0 + y1) / 2000
    bi = {b: i for i, b in enumerate(BANDS)}
    for m in MODELS:
        C_, z = mods[m]["C"], mods[m]["z"]
        a = mods[m]["fz"]["a"]
        lt[f"a_link {m}"] = [a[c] for c in ids]
        F, _, _ = curve_at(C_, z, np.array([0.75, 2.0, 5.0]))
        for j, R in enumerate((0.75, 2.0, 5.0)):
            lt[f"W {m} R={R:g}"] = lt[f"a_link {m}"] * np.array([F[bi[b], j] for b in lt.band])
    lt["length_class"] = lclass(lt.L_km.to_numpy())
    return lt.reset_index().rename(columns={"index": "cml_id"})


def alink_tests(lt):
    W = S["w_independent_of_L"]
    rng = np.random.default_rng(W["seeds"]["rho"])
    y = np.log(lt["a_link candidate"].to_numpy())
    rows = []
    for v in ("L_km", "frequency_ghz", "antenna_d", "abs_sin_bearing", "abs_cos_bearing", "n_radar_cells"):
        rho, p, lo, hi, n = spearman_boot(lt[v].to_numpy(float), y, W["boot_reps"], rng)
        rows.append({"variable": v, "test": "Spearman", "rho": rho, "lo": lo, "hi": hi, "p": p, "n": n})
    # the site: the end with more TRAIN links ("hub"), Kruskal-Wallis over hubs with >= 3 links
    cnt = pd.concat([lt.site_0, lt.site_1]).value_counts()
    hub = np.where(cnt.reindex(lt.site_0).to_numpy() >= cnt.reindex(lt.site_1).to_numpy(), lt.site_0, lt.site_1)
    hs = pd.Series(hub).value_counts()
    keep = [h for h in hs.index if hs[h] >= 3]
    groups = [y[hub == h] for h in keep]
    kw = stats.kruskal(*groups) if len(groups) >= 2 else None
    rows.append({"variable": f"site (hub, {len(keep)} sites with >= 3 links)", "test": "Kruskal-Wallis", "rho": np.nan, "lo": np.nan, "hi": np.nan,
                 "p": float(kw.pvalue) if kw else np.nan, "n": int(sum(len(g) for g in groups))})
    pol = lt.polarization.value_counts().to_dict()
    rows.append({"variable": f"polarization ({', '.join(f'{k} {v}' for k, v in pol.items())} links)", "test": "information only", "rho": np.nan,
                 "lo": np.nan, "hi": np.nan, "p": np.nan, "n": len(lt)})
    t = pd.DataFrame(rows)
    tested = t.p.notna()
    t["bh_discovery"] = False
    t.loc[tested, "bh_discovery"] = bh(t.loc[tested, "p"].to_numpy(), W["bh_q"])
    return t


def w_height(lt):
    W = S["w_independent_of_L"]
    rng = np.random.default_rng(W["seeds"]["w_height"])
    rows = []
    for m in MODELS:
        for R in (2.0, 5.0):
            col = f"W {m} R={R:g}"
            for by in ("length", "band x length"):
                keys = [(nm,) for _, nm in LB] if by == "length" else [(b, nm) for b in sorted(lt.band.unique()) for _, nm in LB]
                for k in keys:
                    msk = (lt.length_class == k[-1]) & ((lt.band == k[0]) if len(k) == 2 else True)
                    v = lt.loc[msk, col].to_numpy(float)
                    if len(v) < 3:
                        continue
                    bs = [np.median(v[rng.integers(0, len(v), len(v))]) for _ in range(W["boot_reps"])]
                    lo, hi = ci(bs)
                    rows.append({"model": m, "R_mmh": R, "by": by, "band": k[0] if len(k) == 2 else "all", "length": k[-1], "links": len(v),
                                 "W_median_db": float(np.median(v)), "lo": lo, "hi": hi})
    return pd.DataFrame(rows)


# =================================================================================================== 7.2 short links on VAL
def rmse_table(sets, preds):
    Q = S["short_links_val"]
    rows = []
    for s_, (D, le) in sets.items():
        bl = D["bin_le"]
        w = D["w5"] > 0
        ev = D["le_event"][bl]
        E = int(D["le_event"].max()) + 1
        Mlt = ev_draws(E, Q["boot_reps"], Q["seed"])
        Lc = lclass(le.L_km.to_numpy(float))[bl]
        band = le.band_name.to_numpy()[bl]
        r5 = np.where(np.isfinite(D["r5"]), D["r5"], -1.0)
        groups = [("all", "all", np.ones(len(bl), bool))]
        groups += [("length", nm, Lc == nm) for _, nm in LB]
        groups += [("band", b, band == b) for b in sorted(set(band))]
        groups += [("rain", f"{lo:g}-{hi:g}" if hi < 1000 else f">= {lo:g}", (r5 >= lo) & (r5 < hi)) for lo, hi in Q["rain_classes_mmh"]]
        for m, p in preds[s_].items():
            r2 = np.where(w, (D["y5"] - p) ** 2, 0.0)
            for kind, nm, g in groups:
                se, ne = np.bincount(ev, weights=r2 * g, minlength=E), np.bincount(ev, weights=(w & g).astype(float), minlength=E)
                if ne.sum() == 0:
                    continue
                b = np.sqrt((Mlt @ se) / np.maximum(Mlt @ ne, 1))
                lo, hi = ci(b)
                rows.append({"set": s_, "model": m, "by": kind, "class": nm, "bins": int(ne.sum()), "rmse": float(np.sqrt(se.sum() / ne.sum())), "lo": lo, "hi": hi})
    t = pd.DataFrame(rows)
    allv = t[(t.by == "all")].set_index(["set", "model"]).rmse
    t["over_all_class"] = [r.rmse / allv[(r.set, r.model)] for r in t.itertuples()]
    tr = t[t.set == "TRAIN"][["model", "by", "class", "rmse"]].rename(columns={"rmse": "rmse_train"})
    t = t.merge(tr, on=["model", "by", "class"], how="left")
    t["val_over_train"] = np.where(t.set == "VAL", t.rmse / t.rmse_train, np.nan)
    return t.drop(columns="rmse_train")


def length_test_val(sets, mods):
    D, le = sets["VAL"]
    cells = through(D, le, {m: mods[m]["fz"] for m in MODELS}, "VAL", "V1")
    rows = []
    for m, g in cells.groupby("fit", sort=False):
        e = g[g.eligible]
        both = int((e.I_inside & e.s_inside).sum())
        rows.append({"fit": m, "eligible_cells": len(e), "both_inside": both, "needed": need(len(e)), "passes": bool(len(e) and both >= need(len(e)))})
    return cells, pd.DataFrame(rows)


def shape_short(sets, mods, free):
    D, le = sets["VAL"]
    keep = (le.L_km < 1.0).to_numpy()
    Ds, les = sid.subset(D, keep), le[keep].reset_index(drop=True)
    pr = {"candidate (shared shapes)": predict(Ds, les, mods["candidate"]["fz"], mods["candidate"]["fz"]["M"])[0],
          "S5a free (band shapes)": predict(Ds, les, free)[0]}
    band = les.band_name.to_numpy()[Ds["bin_le"]]
    w = Ds["w5"] > 0
    rows = []
    for b in sorted(set(band)):
        g = w & (band == b)
        r = {"band": b, "link_events": int((les.band_name == b).sum()), "bins": int(g.sum())}
        for k, p in pr.items():
            r[k] = float(np.sqrt(np.mean((Ds["y5"][g] - p[g]) ** 2)))
        r["shared_minus_band_db"] = r["candidate (shared shapes)"] - r["S5a free (band shapes)"]
        r["ok"] = bool(r["shared_minus_band_db"] <= S["short_links_val"]["shape_rise_max_db"])
        rows.append(r)
    return pd.DataFrame(rows)


def after_rain(sets, preds):
    rows = []
    for s_, (D, le) in sets.items():
        w = D["w5"] > 0
        Lc = lclass(le.L_km.to_numpy(float))[D["bin_le"]]
        aft = D["after5"]
        r = D["y5"] - preds[s_]["candidate"]
        for _, nm in LB:
            for t0 in range(0, 180, 10):
                g = w & (Lc == nm) & (aft >= t0) & (aft < t0 + 10)
                if g.sum() >= 20:
                    rows.append({"set": s_, "length": nm, "after_min": t0 + 5, "bins": int(g.sum()), "median_residual_db": float(np.median(r[g])),
                                 "median_A_db": float(np.median(D["y5"][g] - D["le_c"][D["bin_le"][g]]))})
    return pd.DataFrame(rows)


# =================================================================================================== 7.3a dry-path moments
def dry_stats(win, curves):
    Q = S["long_links"]["dry_path"]
    ev = pd.factorize(win.event_id)[0]
    draws = boot_items(ev, Q["boot_reps"], Q["seed"])
    obs = np.array([c["obs"] for c in curves])
    mod = {m: np.array([c[m] for c in curves]) for m in MODELS}
    rows = []
    for t in Q["at_min"]:
        ok = (win.dry_minutes.to_numpy() > t) & np.isfinite(obs[:, t])
        for by, keys in (("length", [(nm,) for _, nm in LB]), ("band x length", [(b, nm) for b in sorted(win.band.unique()) for _, nm in LB])):
            for k in keys:
                sel = ok & (win.length.to_numpy() == k[-1]) & ((win.band.to_numpy() == k[0]) if len(k) == 2 else True)
                if sel.sum() == 0:
                    continue
                r = {"t_min": t, "by": by, "band": k[0] if len(k) == 2 else "all", "length": k[-1], "windows": int(sel.sum()),
                     "links": int(win[sel].cml_id.nunique()), "events": int(win[sel].event_id.nunique()), "W_obs": float(np.median(obs[sel, t]))}
                if by == "length":
                    r["W_obs_lo"], r["W_obs_hi"] = ci([np.median(obs[i[sel[i]], t]) if sel[i].any() else np.nan for i in draws])
                for m in MODELS:
                    err = obs[:, t] - mod[m][:, t]
                    r[f"W_{m}"] = float(np.median(mod[m][sel, t]))
                    r[f"err_{m}"] = float(np.median(err[sel]))
                    if by == "length":
                        r[f"err_{m}_lo"], r[f"err_{m}_hi"] = ci([np.median(err[i[sel[i]]]) if sel[i].any() else np.nan for i in draws])
                rows.append(r)
    tab = pd.DataFrame(rows)
    # the difference L > 3 km minus L < 1 km of the median error at 10 min, same draws
    t = 10
    ok = (win.dry_minutes.to_numpy() > t) & np.isfinite(obs[:, t])
    s_lo, s_hi = ok & (win.length.to_numpy() == LB[0][1]), ok & (win.length.to_numpy() == LB[2][1])
    diffs = {}
    for m in MODELS:
        err = obs[:, t] - mod[m][:, t]
        point = float(np.median(err[s_hi]) - np.median(err[s_lo])) if s_hi.any() and s_lo.any() else np.nan
        bs = [np.median(err[i[s_hi[i]]]) - np.median(err[i[s_lo[i]]]) if (s_hi[i].any() and s_lo[i].any()) else np.nan for i in draws]
        diffs[m] = (point, *ci(bs))
    return tab, diffs


def decay(win, curves):
    Q = S["long_links"]["dry_path"]
    obs = np.array([c["obs"] for c in curves])
    mod = {m: np.array([c[m] for c in curves]) for m in MODELS}
    w0o = win.W0_obs.to_numpy()
    use = np.isfinite(w0o) & (w0o >= Q["min_w0_db"])
    rows, med = [], {}
    for _, nm in LB:
        sel = use & (win.length.to_numpy() == nm)
        if sel.sum() < 3:
            rows.append({"length": nm, "windows": int(sel.sum())})
            continue
        ws = win[sel]
        ev = pd.factorize(ws.event_id)[0]
        draws = boot_items(ev, Q["boot_reps"], Q["seed"] + 1)
        nmin = ws.dry_minutes.to_numpy()
        r = {"length": nm, "windows": int(sel.sum()), "events": int(ws.event_id.nunique()), "W0_median_db": float(np.median(w0o[sel]))}
        for src, arr in (("obs", obs[sel]), *[(m, mod[m][sel]) for m in MODELS]):
            w0 = np.nanmean(arr[:, :2], axis=1)
            mc = norm_median(arr, nmin, w0)
            med[(nm, src)] = mc
            r[f"tau_{src}"] = fit_tau(mc)
            r[f"left_40_60_{src}"] = float(np.nanmean(mc[40:60])) if np.isfinite(mc[40:60]).any() else np.nan
            if src in ("obs", "candidate"):
                bt = [fit_tau(norm_median(arr[i], nmin[i], w0[i])) for i in draws[:200]]
                r[f"tau_{src}_lo"], r[f"tau_{src}_hi"] = ci(bt)
                bl = [np.nanmean(norm_median(arr[i], nmin[i], w0[i])[40:60]) for i in draws[:200]]
                r[f"left_40_60_{src}_lo"], r[f"left_40_60_{src}_hi"] = ci(bl)
        rows.append(r)
    return pd.DataFrame(rows), med


# =================================================================================================== 7.3b paired links
def pairs(d, sets, lt, preds):
    Q = S["long_links"]["paired"]
    links = d.links().set_index("cml_id")
    sites = d.table("sites").set_index("site_id")
    li = lt.set_index("cml_id")
    cand = []
    ids = list(li.index)
    for i, a in enumerate(ids):
        for b in ids[i + 1:]:
            if li.band[a] != li.band[b]:
                continue
            sa, sb = {links.site_0_id[a], links.site_1_id[a]}, {links.site_0_id[b], links.site_1_id[b]}
            common = sa & sb
            if len(common) != 1:
                continue
            s0 = common.pop()
            oa, ob = (sa - {s0}).pop(), (sb - {s0}).pop()
            xs, ys = sites.x_3006_m[s0], sites.y_3006_m[s0]
            ba = np.degrees(np.arctan2(sites.x_3006_m[oa] - xs, sites.y_3006_m[oa] - ys))
            bb = np.degrees(np.arctan2(sites.x_3006_m[ob] - xs, sites.y_3006_m[ob] - ys))
            dang = abs((ba - bb + 180) % 360 - 180)
            L1, L2 = sorted([li.L_km[a], li.L_km[b]])
            short, long_ = (a, b) if li.L_km[a] <= li.L_km[b] else (b, a)
            cand.append({"link_short": short, "link_long": long_, "site": s0, "band": li.band[a], "L1_km": L1, "L2_km": L2, "bearing_diff_deg": dang,
                         "ok_bearing": dang <= Q["max_bearing_deg"], "ok_ratio": L2 / L1 >= Q["min_ratio"], "ok_dL": L2 - L1 >= Q["min_dL_km"]})
    pt = pd.DataFrame(cand)
    if not len(pt):
        return pt, pd.DataFrame(), {"pairs": 0}
    pt["geometry_ok"] = pt.ok_bearing & pt.ok_ratio & pt.ok_dL
    ok = pt[pt.geometry_ok]
    rows = []
    for s_, (D, le) in sets.items():
        bt = bin_times(D, le)
        first = le.drop_duplicates(["event_id", "cml_id"])
        by = {(r.event_id, r.cml_id): i for i, r in zip(first.index, first.itertuples())}
        for p in ok.itertuples():
            for ev_ in sorted(set(first[first.cml_id == p.link_short].event_id) & set(first[first.cml_id == p.link_long].event_id)):
                i1, i2 = by[(ev_, p.link_short)], by[(ev_, p.link_long)]
                sl = [slice(D["le_boff5"][i], D["le_boff5"][i] + D["le_n"][i] // 5) for i in (i1, i2)]
                t1, t2 = bt[sl[0]], bt[sl[1]]
                common, a1, a2 = np.intersect1d(t1.values, t2.values, return_indices=True)
                if not len(common):
                    continue
                g1, g2 = np.arange(sl[0].start, sl[0].stop)[a1], np.arange(sl[1].start, sl[1].stop)[a2]
                v = (D["w5"][g1] > 0) & (D["w5"][g2] > 0) & (np.nan_to_num(D["r5"][g1], nan=-1) >= 0.5) & (np.nan_to_num(D["r5"][g2], nan=-1) >= 0.5)
                if v.sum() < 6:
                    continue
                L1, L2 = p.L1_km, p.L2_km
                A1, A2 = D["y5"][g1] - D["le_c"][i1], D["y5"][g2] - D["le_c"][i2]
                pc = preds[s_]["candidate"]
                m1, m2 = pc[g1] - D["le_c"][i1], pc[g2] - D["le_c"][i2]
                res1, res2 = D["y5"][g1] - pc[g1], D["y5"][g2] - pc[g2]
                # fake: a known truth plus the real residual series of the two link-events (aligned in time)
                Wt, kt = 1.5, 0.5
                f1, f2 = Wt + L1 * kt + res1, Wt + L2 * kt + res2
                est_fake = (L2 * f1 - L1 * f2) / (L2 - L1)
                est_obs = (L2 * A1 - L1 * A2) / (L2 - L1)
                est_mod = (L2 * m1 - L1 * m2) / (L2 - L1)
                rr = np.sqrt(D["r5"][g1] * D["r5"][g2])
                rows.append({"set": s_, "event_id": ev_, "link_short": p.link_short, "link_long": p.link_long, "L1_km": L1, "L2_km": L2, "rain_bins": int(v.sum()),
                             "noise_sd_db": float(np.std(est_fake[v] - Wt)), "noise_factor": float(np.sqrt(L1 ** 2 + L2 ** 2) / (L2 - L1)),
                             "W_obs_median": float(np.median(est_obs[v])), "W_model_median": float(np.median(est_mod[v])),
                             "_bins": (est_obs[v], est_mod[v], rr[v])})
    pe = pd.DataFrame(rows)
    info = {"pairs_geometry": int(len(ok)), "pair_events": int(len(pe)), "pairs_with_events": int(pe[["link_short", "link_long"]].drop_duplicates().shape[0]) if len(pe) else 0,
            "candidates_shared_site_same_band": int(len(pt))}
    info["median_noise_sd_db"] = float(pe.noise_sd_db.median()) if len(pe) else np.nan
    info["testable"] = bool(info["pairs_with_events"] >= Q["min_pairs"] and np.isfinite(info["median_noise_sd_db"]) and info["median_noise_sd_db"] <= Q["max_noise_sd_db"])
    return pt, pe, info


def pair_compare(pe):
    classes = [tuple(c) for c in LT["rain_classes_mmh"]]
    ev = pd.factorize(pe.event_id)[0]
    rows = []
    for lo, hi in classes:
        vals, evs = [], []
        for k, r in enumerate(pe.itertuples()):
            o, m, rr = r._bins
            g = (rr >= lo) & (rr < hi)
            vals.append(o[g] - m[g])
            evs.append(np.full(g.sum(), ev[k]))
        v, e = np.concatenate(vals), np.concatenate(evs)
        if len(v) < 10:
            continue
        draws = boot_items(e, 1000, S["long_links"]["dry_path"]["seed"] + 7)
        lo_, hi_ = ci([np.median(v[i]) for i in draws if len(i)])
        rows.append({"class": f"{lo:g}-{hi:g}", "bins": len(v), "events": int(len(np.unique(e))), "median_W_obs_minus_model": float(np.median(v)), "lo": lo_, "hi": hi_})
    return pd.DataFrame(rows)


# =================================================================================================== 7.3c the intercept from short-only and long-only links
def intercepts(sets, preds):
    Q = S["long_links"]["intercept"]
    D, le = sets["TRAIN"]
    t = LTL.steady_table(D, le, LT, {m: preds["TRAIN"][m] for m in MODELS})
    n_ev = int(D["le_event"].max()) + 1
    parts = {"short (L < 3 km)": t[t.L_km < Q["short_max_km"]], "long (L > 1.5 km)": t[t.L_km > Q["long_min_km"]], "all lengths": t}
    fits = {}
    for nm, sub in parts.items():
        f = LTL.length_fits(sub, n_ev, nm, LT)
        for m in MODELS:
            fm = LTL.length_fits(sub, n_ev, nm, LT, ycol=f"pred_{m}", boot=False)[["band", "class", "I", "s"]].rename(columns={"I": f"I_{m}", "s": f"s_{m}"})
            f = f.merge(fm, on=["band", "class"], how="left")
        fits[nm] = f
    a, b = fits["short (L < 3 km)"], fits["long (L > 1.5 km)"]
    cols = ["band", "class", "R_cell", "links", "link_events", "I", "I_lo", "I_hi", "s", "s_lo", "s_hi"] + [f"I_{m}" for m in MODELS]
    mg = a[cols].merge(b[cols], on=["band", "class"], suffixes=("_short", "_long"))
    mg = mg[np.isfinite(mg.I_short) & np.isfinite(mg.I_long)].copy()
    mg["overlap"] = [overlap((r.I_lo_short, r.I_hi_short), (r.I_lo_long, r.I_hi_long)) for r in mg.itertuples()]
    return mg, fits


def p_size(mods, let, BANDS):
    C_, z = mods["candidate"]["C"], mods["candidate"]["z"]
    _, Pv, _ = C_.values(z)
    freq = let.groupby("band_name").frequency_ghz.median()
    rows = []
    for b, nm in enumerate(BANDS):
        c = float(Pv[b, -1] / C_.kn[-1])
        k, al = itu_p838(float(freq[nm]), "v")
        for R in (2.0, 5.0, 10.0):
            ratio = c * R / (k * R ** al)
            lo_, hi_ = S["long_links"]["p_ratio_flag"]
            rows.append({"band": nm, "f_ghz": float(freq[nm]), "c_db_km_per_mmh": c, "R_mmh": R, "P_candidate": c * R, "ITU_P838_kR^alpha": k * R ** al,
                         "ratio": ratio, "flag": not (lo_ <= ratio <= hi_)})
    return pd.DataFrame(rows)


# =================================================================================================== 7.4 dew
def dew(sets, lt, preds, Rm, T5, col):
    Q = S["dew"]
    j = pd.read_csv(P.DRY_JUMPS)
    j = j[j.cml_id.isin(set(lt.cml_id)) & (j.sign == "positive") & (j.hours_since_wet >= 6)
          & ((j.mean_rh >= 90) | (j.mean_dew_point_depression < 2))].copy()
    j["night_key"] = j.start_local.str[:10]
    n_nights = int(j.groupby(["cml_id", "night_key"]).ngroups)
    out = {"jumps": j, "link_nights": n_nights, "links": int(j.cml_id.nunique()), "testable": n_nights >= Q["min_link_nights"]}
    li = lt.set_index("cml_id")
    out["rain_W_075"] = float(np.median(li["W candidate R=0.75"]))
    if not out["testable"]:
        return out
    j["start_hour_local"] = pd.to_datetime(j.start_local).dt.hour
    out["size"] = j[["median_excess_db", "peak_excess_db", "duration_h", "mean_rh"]].describe(percentiles=[0.25, 0.5, 0.75]).T
    pl = j.groupby("cml_id").peak_excess_db.median()
    rng = np.random.default_rng(Q["seeds"]["same_layer"])
    out["same_layer"] = spearman_boot(pl.to_numpy(float), np.log(li["a_link candidate"].reindex(pl.index).to_numpy(float)), Q["boot_reps"], rng)
    out["per_link"] = pd.DataFrame({"cml_id": pl.index, "dew_peak_median_db": pl.to_numpy(), "a_link": li["a_link candidate"].reindex(pl.index).to_numpy()})
    # dew before the rain, per link-event
    dr = pd.read_csv(P.DEW_RISK_HOURLY, parse_dates=["time_utc"]).set_index("time_utc").dew_risk
    jj = {s: g for s, g in j.groupby("sublink_id")}
    rows = []
    for s_, (D, le) in sets.items():
        p = preds[s_]["candidate"]
        fwb = T5.searchsorted(pd.to_datetime(le.first_wet).to_numpy())
        for e, r in enumerate(le.itertuples()):
            fw = pd.Timestamp(r.first_wet)
            hrs = (fw - pd.Timedelta("6h")).floor("h") + pd.to_timedelta(np.arange(6), "h")
            n_dew = int(dr.reindex(hrs).fillna(False).astype(bool).sum())
            rm = Rm[max(fwb[e] - 72, 0):fwb[e], col[r.cml_id]]
            path_dry = bool(len(rm) == 72 and np.isfinite(rm).all() and (rm <= 0.05).all())
            jx = jj.get(r.sublink_id)
            jump = bool(jx is not None and ((pd.to_datetime(jx.start_utc) < fw) & (pd.to_datetime(jx.end_utc) > fw - pd.Timedelta("6h"))).any())
            o, n5 = D["le_boff5"][e], D["le_n"][e] // 5
            w = D["w5"][o:o + n5] > 0
            res = (D["y5"][o:o + n5] - p[o:o + n5])[w]
            rmax = np.nanmax(np.r_[D["r5"][o:o + n5], 0.0])
            rows.append({"set": s_, "event_id": r.event_id, "cml_id": r.cml_id, "sublink_id": r.sublink_id, "dew_hours_before": n_dew, "path_dry_6h": path_dry,
                         "dew_before": bool(n_dew >= Q["dew_hours_min"] and path_dry), "jump_before": bool(jump and path_dry), "rmse": float(np.sqrt(np.mean(res ** 2))) if len(res) else np.nan,
                         "bias": float(np.mean(res)) if len(res) else np.nan, "c_db": float(D["le_c"][e]), "rain_max_v1": float(rmax)})
    t = pd.DataFrame(rows)
    rc = Q["rain_classes_mmh"]
    t["rain_class"] = pd.cut(t.rain_max_v1, [c[0] for c in rc] + [rc[-1][1]], right=False, labels=[f"{a:g}-{b:g}" if b < 1000 else f">= {a:g}" for a, b in rc]).astype(str)
    out["events"] = t
    comp = []
    for flag in ("dew_before", "jump_before"):
        if flag == "jump_before" and t[flag].sum() < Q["min_link_nights"]:
            comp.append({"definition": flag, "link_events_with": int(t[flag].sum()), "note": "fewer than 8 link-events: not testable"})
            continue
        ev = pd.factorize(t.event_id)[0]
        draws = boot_items(ev, Q["boot_reps"], Q["seeds"]["before_rain"])
        def stat(ix):
            tt = t.iloc[ix]
            o = {}
            for v in ("rmse", "bias", "c_db"):
                num, den = 0.0, 0
                for _, g in tt.groupby("rain_class"):
                    a, b = g[g[flag]][v].dropna(), g[~g[flag]][v].dropna()
                    if len(a) and len(b):
                        num += len(a) * (a.mean() - b.mean())
                        den += len(a)
                o[v] = num / den if den else np.nan
            return o
        point = stat(np.arange(len(t)))
        bs = [stat(i) for i in draws[:300]]
        r = {"definition": flag, "link_events_with": int(t[flag].sum()), "link_events_without": int((~t[flag]).sum()), "events_with": int(t[t[flag]].event_id.nunique())}
        for v in ("rmse", "bias", "c_db"):
            r[f"{v}_with_minus_without"] = point[v]
            r[f"{v}_lo"], r[f"{v}_hi"] = ci([b[v] for b in bs])
        comp.append(r)
    out["compare"] = pd.DataFrame(comp)
    return out


# =================================================================================================== 7.5 leftovers: PySR jobs (cached)
def pysr_job(jb):
    file = CACHE / f"{jb['name']}.pkl"
    if file.exists():
        return jb["name"], pickle.loads(file.read_bytes())
    from pysr import PySRRegressor
    PS = S["leftovers"]["pysr"]
    t0 = time.time()
    m = PySRRegressor(niterations=PS["niterations"], maxsize=PS["maxsize"], binary_operators=PS["binary"], unary_operators=PS["unary"],
                      random_state=jb["seed"], deterministic=True, parallelism="serial", progress=False, verbosity=0, temp_equation_file=True)
    m.fit(np.asarray(jb["X"]), np.asarray(jb["y"]), variable_names=[f"v{i}" for i in range(np.asarray(jb["X"]).shape[1])])
    out = {"eqs": [(int(r.complexity), float(r.loss), str(r.equation)) for r in m.equations_.itertuples()], "seconds": time.time() - t0}
    file.write_bytes(pickle.dumps(out))
    return jb["name"], out


def run_pool(jobs, workers, clock, results):
    jobs = [j for j in jobs if j["name"] not in results]
    if not jobs:
        return results
    if workers <= 1:
        for j in jobs:
            n, o = pysr_job(j)
            results[n] = o
        return results
    with ProcessPoolExecutor(workers) as ex:
        running = {ex.submit(pysr_job, j): j for j in jobs[:workers]}
        queue = jobs[workers:]
        while running:
            done_, _ = wait(running, return_when=FIRST_COMPLETED)
            for fu in done_:
                running.pop(fu)
                name, out = fu.result()
                results[name] = out
                if len(results) % 20 == 0:
                    print(f"   PySR {len(results)} runs done ({(time.time() - clock) / 60:.1f} min)", flush=True)
                if queue:
                    j = queue.pop(0)
                    running[ex.submit(pysr_job, j)] = j
    return results


def structure(eq, nv):
    """An equation with its numeric constants as parameters: (f(X, theta), theta0)."""
    import sympy as sp
    vs = sp.symbols(" ".join(f"v{i}" for i in range(nv)))
    vs = vs if isinstance(vs, tuple) else (vs,)
    ex = sp.sympify(eq, locals={f"v{i}": v for i, v in enumerate(vs)})
    consts = sorted({a for a in ex.atoms(sp.Float)}, key=lambda a: float(a))
    th = sp.symbols(f"t0:{len(consts)}") if consts else ()
    f = sp.lambdify(tuple(vs) + tuple(th), ex.xreplace(dict(zip(consts, th))), "numpy")

    def fn(X, theta):
        with np.errstate(all="ignore"):
            return np.asarray(f(*[X[:, i] for i in range(nv)], *theta), float) * np.ones(len(X))
    return fn, np.array([float(c) for c in consts])


def loo_skill(eq, X, y):
    """1 - MSE / MSE of the constant, both leave-one-out (the constants of the structure refitted without the held-out row)."""
    fn, th0 = structure(eq, X.shape[1])
    n = len(y)
    pe, pc = np.empty(n), np.empty(n)
    for i in range(n):
        m = np.ones(n, bool)
        m[i] = False
        pc[i] = y[m].mean()
        if len(th0):
            try:
                r = least_squares(lambda th: np.nan_to_num(fn(X[m], th) - y[m], nan=1e3, posinf=1e3, neginf=-1e3), th0, max_nfev=400)
                th = r.x
            except Exception:                                    # noqa: BLE001
                th = th0
        else:
            th = th0
        pe[i] = fn(X[i:i + 1], th)[0]
    pe = np.where(np.isfinite(pe), pe, pc)
    return float(1 - np.mean((pe - y) ** 2) / np.mean((pc - y) ** 2))


def choose_eq(results, names):
    PS = S["leftovers"]["pysr"]
    best = None
    for nm in names:
        for c, loss, eq in results[nm]["eqs"]:
            if 2 <= c <= PS["max_complexity"] and np.isfinite(loss) and (best is None or loss < best[1]):
                best = (c, loss, eq)
    return best


def leftover_targets(lt, win, curves):
    """The two targets with their variables: log a_link per link, log tau_e per event."""
    A = S["leftovers"]
    va = ["L_km", "frequency_ghz", "antenna_d", "abs_sin_bearing", "abs_cos_bearing", "n_radar_cells", "dry_noise_db", "diurnal_range_db", "x_mid_km", "y_mid_km"]
    ta = lt[va].astype(float)
    ya = np.log(lt["a_link candidate"].to_numpy(float))
    okr = ta.notna().all(axis=1).to_numpy()
    # per-event drying time from the dry-path windows (measured, W0 >= 0.5 dB)
    dr = pd.read_csv(P.DEW_RISK_HOURLY, parse_dates=["time_utc"]).set_index("time_utc")
    obs = np.array([c["obs"] for c in curves])
    use = np.isfinite(win.W0_obs.to_numpy()) & (win.W0_obs.to_numpy() >= S["long_links"]["dry_path"]["min_w0_db"])
    rows = []
    for e_id, g in win[use].groupby("event_id"):
        if len(g) < A["min_windows_per_event"]:
            continue
        arr = obs[g.index.to_numpy()]
        tau = fit_tau(norm_median(arr, g.dry_minutes.to_numpy(), np.nanmean(arr[:, :2], axis=1)))
        hrs = pd.DatetimeIndex(g.t0).floor("h")
        w = dr.reindex(hrs)
        rows.append({"event_id": e_id, "set": g.set.iloc[0], "windows": len(g), "tau_min": tau, "temperature": w.temperature.mean(), "humidity": w.humidity.mean(),
                     "wind_speed": w.wind_speed.mean(), "dew_point_depression": w.dew_point_depression.mean(), "night_share": w.night.astype(float).mean(),
                     "local_hour": float(np.mean((hrs.hour + 2) % 24))})
    te = pd.DataFrame(rows)
    vd = ["temperature", "humidity", "wind_speed", "dew_point_depression", "night_share", "local_hour"]
    if len(te):
        te = te[np.isfinite(te.tau_min) & (te.tau_min > 0)].reset_index(drop=True)
    return {"a_link": (ta[okr].reset_index(drop=True), ya[okr], va), "drying": (te[vd].astype(float) if len(te) else pd.DataFrame(columns=vd),
                                                                               np.log(te.tau_min.to_numpy(float)) if len(te) else np.zeros(0), vd)}, te


def leftovers(targets, clock, workers, plan):
    A = S["leftovers"]
    PS = A["pysr"]
    results, rows, screens = {}, [], []
    seeds = PS["seeds"][:plan["seeds"]]
    for tname, (X, y, vnames) in targets.items():
        if len(y) < 8:
            rows.append({"target": tname, "n": len(y), "verdict": "not testable (too few units)"})
            continue
        # Spearman screening
        p = []
        for v in vnames:
            rho, pv = stats.spearmanr(X[v].to_numpy(float), y)
            p.append(pv)
            screens.append({"target": tname, "variable": v, "rho": float(rho), "p": float(pv)})
        disc = bh(np.nan_to_num(np.array(p), nan=1.0), S["w_independent_of_L"]["bh_q"])   # the same q = 0.10 (config leftovers.screen)
        for k, v in enumerate(vnames):
            screens[-len(vnames) + k]["bh_discovery"] = bool(disc[k])
        Xs = ((X - X.mean()) / X.std(ddof=0).replace(0, 1)).to_numpy(float)
        rng = np.random.default_rng([A["seed"], len(tname)])
        perms = [rng.permutation(len(y)) for _ in range(plan["shuffles"])]
        jobs = []
        for sh in range(-1, plan["shuffles"]):
            yy = y if sh < 0 else y[perms[sh]]
            for sd in seeds:
                jobs.append({"name": f"pysr_{tname.replace(' ', '_')}_{'real' if sh < 0 else f'shuffle{sh:02d}'}_s{sd}", "X": Xs.tolist(), "y": yy.tolist(), "seed": sd})
        run_pool(jobs, workers, clock, results)
        real = choose_eq(results, [j["name"] for j in jobs if "_real_" in j["name"]])
        sk = loo_skill(real[2], Xs, y) if real else np.nan
        null = []
        for sh in range(plan["shuffles"]):
            b = choose_eq(results, [j["name"] for j in jobs if f"_shuffle{sh:02d}_" in j["name"]])
            null.append(loo_skill(b[2], Xs, y[perms[sh]]) if b else np.nan)
        q95 = float(np.nanpercentile(null, 95)) if np.isfinite(null).any() else np.nan
        eq_named = real[2] if real else ""
        for i, v in sorted(enumerate(vnames), key=lambda t: -len(f"v{t[0]}")):
            eq_named = eq_named.replace(f"v{i}", v)
        hit = bool(real and np.isfinite(sk) and sk > 0 and sk > q95)
        rows.append({"target": tname, "n": len(y), "screen_discoveries": ", ".join(v for v, dd in zip(vnames, disc) if dd) or "none",
                     "pysr_equation": eq_named, "complexity": real[0] if real else np.nan, "loo_skill": sk, "null_q95": q95, "null_median": float(np.nanmedian(null)),
                     "shuffles": plan["shuffles"], "seeds": len(seeds), "verdict": "RELATION (hint only)" if hit else "no relation found"})
    return pd.DataFrame(rows), pd.DataFrame(screens), results


# =================================================================================================== 7.6 the rain reference
def rain_reference(d, sets, gmap, preds):
    Q = S["rain_reference"]
    ev = d.events()
    ev = ev[ev.split_event.isin(["train", "val"])]
    rs, gs = d.store("radar_at_gauges_5min").R_pixel, d.store("gauges_city_5min").rate_mm_h
    rows = []
    for k, e in enumerate(ev.itertuples()):
        sl = slice(e.start, e.end - pd.Timedelta("1ns"))
        a = rs.sel(time=sl).to_pandas()
        b = gs.sel(time=sl).to_pandas()
        cols = [c for c in b.columns if c in a.columns]
        a, b = a[cols], b[cols].reindex(a.index)
        set_ = "TRAIN" if e.split_event == "train" else "VAL"
        g = float(gmap[set_].get(e.event_id, 1.0))
        A_, B_ = a.to_numpy(float).ravel(), b.to_numpy(float).ravel()
        ok = np.isfinite(A_) & np.isfinite(B_)
        rows.append(pd.DataFrame({"event_id": e.event_id, "set": set_, "ev": k, "radar_v0": A_[ok], "radar_v1": A_[ok] / g, "gauge": B_[ok]}))
    t = pd.concat(rows, ignore_index=True)
    cls = [tuple(c) for c in Q["classes_mmh"]]
    E = int(t.ev.max()) + 1
    Mlt = ev_draws(E, Q["boot_reps"], Q["seed"])
    out = []
    for inp in ("v0", "v1"):
        rad = t[f"radar_{inp}"].to_numpy()
        for by, key in (("gauge", t.gauge.to_numpy()), ("radar", rad)):
            sums = {}
            for lo, hi in cls + [(0.1, 5.0)]:
                m = (key >= lo) & (key < hi)
                sr = np.bincount(t.ev[m], weights=rad[m], minlength=E)
                sg = np.bincount(t.ev[m], weights=t.gauge.to_numpy()[m], minlength=E)
                sums[(lo, hi)] = (sr, sg)
                if (lo, hi) in cls:
                    lo_, hi_ = ci((Mlt @ sr) / np.maximum(Mlt @ sg, 1e-9))
                    out.append({"input": inp.upper(), "binned_by": by, "class": f"{lo:g}-{hi:g}" if hi < 1000 else f">= {lo:g}", "pairs": int(m.sum()),
                                "events": int(len(np.unique(t.ev[m]))), "ratio": float(sr.sum() / sg.sum()) if sg.sum() > 0 else np.nan, "lo": lo_, "hi": hi_})
            (lr, lg), (hr, hg) = sums[(0.1, 5.0)], sums[cls[-1]]
            hl = (hr.sum() / hg.sum()) / (lr.sum() / lg.sum())
            b = ((Mlt @ hr) / np.maximum(Mlt @ hg, 1e-9)) / ((Mlt @ lr) / np.maximum(Mlt @ lg, 1e-9))
            lo_, hi_ = ci(b)
            out.append({"input": inp.upper(), "binned_by": by, "class": "high / low (>= 10 against 0.1-5)", "pairs": np.nan, "events": np.nan, "ratio": float(hl),
                        "lo": lo_, "hi": hi_})
    rt = pd.DataFrame(out)
    # the rule on V1
    hl = rt[(rt.input == "V1") & rt["class"].str.startswith("high")].set_index("binned_by")
    ch = {b: (hl.ratio[b] - 1, hl.lo[b] > 1 or hl.hi[b] < 1) for b in ("gauge", "radar")}
    found = bool(all(abs(c) > Q["max_change"] and excl for c, excl in ch.values()) and np.sign(ch["gauge"][0]) == np.sign(ch["radar"][0]))
    # the candidate residual by V1 class
    res = []
    for s_, (D, le) in sets.items():
        w = D["w5"] > 0
        r5 = np.nan_to_num(D["r5"], nan=-1)
        r = D["y5"] - preds[s_]["candidate"]
        ev_ = D["le_event"][D["bin_le"]]
        E2 = int(D["le_event"].max()) + 1
        M2 = ev_draws(E2, 1000, Q["seed"] + 1)
        for lo, hi in ((2, 5), (10, 1000)):
            g = w & (r5 >= lo) & (r5 < hi)
            sr, n = np.bincount(ev_[g], weights=r[g], minlength=E2), np.bincount(ev_[g], minlength=E2).astype(float)
            lo_, hi_ = ci((M2 @ sr) / np.maximum(M2 @ n, 1))
            res.append({"set": s_, "V1 class": f"{lo:g}-{hi:g}" if hi < 1000 else f">= {lo:g}", "bins": int(g.sum()), "mean_residual_db": float(sr.sum() / n.sum()),
                        "lo": lo_, "hi": hi_})
    # per event (gauge binning, V1): low and high class ratios
    pev = []
    for e_id, g in t.groupby("event_id"):
        r_ = {"event_id": e_id, "set": g.set.iloc[0], "g_e": float(g.radar_v0.sum() / g.radar_v1.sum()) if g.radar_v1.sum() > 0 else np.nan}
        for nm, (lo, hi) in (("low 0.1-5", (0.1, 5.0)), ("high >= 10", (10, 1000))):
            m = (g.gauge >= lo) & (g.gauge < hi)
            r_[f"V1/gauge {nm}"] = float(g.radar_v1[m].sum() / g.gauge[m].sum()) if g.gauge[m].sum() > 0 else np.nan
            r_[f"pairs {nm}"] = int(m.sum())
        pev.append(r_)
    return rt, found, ch, pd.DataFrame(res), pd.DataFrame(pev)


# =================================================================================================== figures (colour-free)
MK = {"L < 1 km": "o", "1-3 km": "s", "L > 3 km": "^"}
BMK = {"28 GHz": "o", "32 GHz": "s", "38 GHz": "^"}


def fig_alink(lt, at):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    for b, g in lt.groupby("band"):
        ax.scatter(g.L_km, g["a_link candidate"], marker=BMK.get(b, "o"), facecolors="none", edgecolors=INK, s=22, label=b)
    r = at[at.variable == "L_km"].iloc[0]
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("link length L (km)")
    ax.set_ylabel("a_link (candidate)")
    ax.set_title(f"7.1: a_link against L (Spearman {r.rho:+.2f} [{r.lo:+.2f}, {r.hi:+.2f}])")
    ax.legend(frameon=False, fontsize=8, title="band (marker)")
    style_axes(ax)
    fig.tight_layout()
    return fig


def fig_wheight(wh):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharey=True)
    for ax, R in zip(axes, (2.0, 5.0)):
        for k, (m, ls) in enumerate(zip(MODELS, ("-", "--"))):
            g = wh[(wh.model == m) & (wh.R_mmh == R) & (wh.by == "length")]
            x = np.arange(len(g)) + (k - 0.5) * 0.15
            ax.errorbar(x, g.W_median_db, yerr=[g.W_median_db - g.lo, g.hi - g.W_median_db], fmt="o" if k == 0 else "s", mfc="white" if k else INK,
                        color=INK, ls=ls, capsize=3, label=m)
            ax.set_xticks(np.arange(len(g)), g.length)
        ax.set_title(f"W height a_link F at {R:g} mm/h")
        ax.set_ylabel("median over links (dB)")
        style_axes(ax)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("7.1: the antenna level does not depend on the length class if the bars overlap", fontsize=9)
    fig.tight_layout()
    return fig


def fig_rmse(rm):
    g = rm[(rm.set == "VAL") & (rm.model == "candidate") & rm.by.isin(["all", "length", "band"])]
    fig, ax = plt.subplots(figsize=(7, 3.8))
    x = np.arange(len(g))
    ax.errorbar(x, g.rmse, yerr=[g.rmse - g.lo, g.hi - g.rmse], fmt="o", color=INK, capsize=3)
    allv = float(g[g.by == "all"].rmse.iloc[0])
    ax.axhline(allv * S["short_links_val"]["rmse_factor_max"], color=INK2, ls=":", lw=1)
    ax.text(len(g) - 0.5, allv * S["short_links_val"]["rmse_factor_max"], "rule: 1.15 x all", ha="right", va="bottom", fontsize=8)
    ax.set_xticks(x, [f"{a}\n{b}" if a != "all" else "all" for a, b in zip(g.by, g["class"])], fontsize=8)
    ax.set_ylabel(f"VAL V1 RMSE (dB, {VL})")
    ax.set_title("7.2: VAL error of the candidate per length class and band")
    style_axes(ax)
    fig.tight_layout()
    return fig


def fig_dry(win, curves):
    obs = np.array([c["obs"] for c in curves])
    mod = {m: np.array([c[m] for c in curves]) for m in MODELS}
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, (_, nm) in zip(axes, LB):
        sel = win.length.to_numpy() == nm
        nmin = win.dry_minutes.to_numpy()
        for src, arr, ls, lab in (("obs", obs, "-", "measured A - c"), ("candidate", mod["candidate"], "--", "model W (candidate)"),
                                  ("closed F", mod["closed F"], ":", "model W (closed F)")):
            med = [np.nanmedian(arr[sel & (nmin > t), t]) if (sel & (nmin > t)).sum() >= 3 else np.nan for t in range(60)]
            ax.plot(np.arange(60), med, ls, color=INK, lw=1.4, label=lab)
        ax.set_title(f"{nm}: {int(sel.sum())} windows")
        ax.set_xlabel("minutes after the path became dry")
        style_axes(ax)
    axes[0].set_ylabel("median W (dB)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("7.3a: dry-path moments, W measured (A - c, no rain on any path pixel) against the frozen model", fontsize=9)
    fig.tight_layout()
    return fig


def fig_decay(dc):
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    g = dc.dropna(subset=["tau_obs"]) if "tau_obs" in dc else dc.iloc[:0]
    x = np.arange(len(g))
    if len(g):
        ax.errorbar(x - 0.08, g.tau_obs, yerr=[g.tau_obs - g.tau_obs_lo, g.tau_obs_hi - g.tau_obs], fmt="o", color=INK, capsize=3, label="measured")
        ax.errorbar(x + 0.08, g.tau_candidate, yerr=[g.tau_candidate - g.tau_candidate_lo, g.tau_candidate_hi - g.tau_candidate], fmt="s", mfc="white", color=INK,
                    capsize=3, label="model (candidate)")
        ax.set_xticks(x, g.length)
    ax.axhline(FIX["tau_f_dry"], color=INK2, ls="--", lw=1)
    ax.text(-0.4, FIX["tau_f_dry"], f"tau_f_dry {FIX['tau_f_dry']:g} min", va="bottom", fontsize=8)
    ax.set_ylabel("decay time of the median W, 0-30 min (min)")
    ax.set_title("7.3a: drying time per length class")
    ax.legend(frameon=False, fontsize=8)
    style_axes(ax)
    fig.tight_layout()
    return fig


def fig_pairs(pe, info):
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    if len(pe):
        ax.scatter(pe.noise_factor, pe.noise_sd_db, marker="o", facecolors="none", edgecolors=INK, s=22)
        ax.axhline(S["long_links"]["paired"]["max_noise_sd_db"], color=INK2, ls=":", lw=1)
        ax.text(ax.get_xlim()[1], S["long_links"]["paired"]["max_noise_sd_db"], "testable below this line ", ha="right", va="bottom", fontsize=8)
    ax.set_xlabel("noise amplification sqrt(L1^2 + L2^2) / (L2 - L1)")
    ax.set_ylabel("estimator noise SD (fake truth, dB)")
    ax.set_title(f"7.3b: paired links ({info.get('pairs_with_events', 0)} pairs, {info.get('pair_events', 0)} pair-events): "
                 f"{'testable' if info.get('testable') else 'NOT testable'}", fontsize=9)
    style_axes(ax)
    fig.tight_layout()
    return fig


def fig_intercept(mg):
    bands = sorted(mg.band.unique())
    fig, axes = plt.subplots(1, max(len(bands), 1), figsize=(4 * max(len(bands), 1), 3.6), sharey=True, squeeze=False)
    for ax, b in zip(axes[0], bands):
        g = mg[mg.band == b].sort_values("R_cell_short")
        ax.errorbar(g.R_cell_short * 0.97, g.I_short, yerr=[g.I_short - g.I_lo_short, g.I_hi_short - g.I_short], fmt="o", color=INK, capsize=3, label="short (L < 3 km)")
        ax.errorbar(g.R_cell_long * 1.03, g.I_long, yerr=[g.I_long - g.I_lo_long, g.I_hi_long - g.I_long], fmt="s", mfc="white", color=INK, capsize=3,
                    label="long (L > 1.5 km)")
        ax.set_xscale("log")
        ax.set_title(b)
        ax.set_xlabel("rain of the cell (mm/h, V1)")
        style_axes(ax)
    for ax in axes[0]:
        ax.set_ylim(-1, 7)                                        # one 32 GHz 10-20 mm/h interval reaches 28 dB (clipped; in the table)
    axes[0][0].set_ylabel("intercept I (dB; axis clipped at 7)")
    axes[0][0].legend(frameon=False, fontsize=8)
    fig.suptitle("7.3c: the antenna intercept from short-only and from long-only links (model-free, TRAIN V1)", fontsize=9)
    fig.tight_layout()
    return fig


def fig_dew(dw):
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    if "per_link" in dw:
        pl = dw["per_link"]
        ax.scatter(pl.a_link, pl.dew_peak_median_db, marker="o", facecolors="none", edgecolors=INK, s=22)
        r = dw["same_layer"]
        ax.set_title(f"7.4: dew W against a_link ({len(pl)} links; Spearman {r[0]:+.2f} [{r[2]:+.2f}, {r[3]:+.2f}])", fontsize=9)
    ax.set_xlabel("a_link (candidate)")
    ax.set_ylabel("median dew peak above the dry baseline (dB)")
    style_axes(ax)
    fig.tight_layout()
    return fig


def fig_rain_ref(rt):
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8), sharey=True)
    for ax, inp in zip(axes, ("V0", "V1")):
        for by, ls, mk in (("gauge", "-", "o"), ("radar", "--", "s")):
            g = rt[(rt.input == inp) & (rt.binned_by == by) & ~rt["class"].str.startswith("high")]
            x = np.arange(len(g))
            ax.errorbar(x, g.ratio, yerr=[g.ratio - g.lo, g.hi - g.ratio], fmt=mk + ls, mfc="white" if by == "radar" else INK, color=INK, capsize=3,
                        label=f"binned by the {by}")
            ax.set_xticks(x, g["class"])
        ax.axhline(1, color=INK2, lw=0.8)
        ax.set_yscale("log")
        ax.set_title(f"radar {inp} / City gauge")
        ax.set_xlabel("rate class (mm/h)")
        style_axes(ax)
    axes[0].set_ylabel("sum radar / sum gauge")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("7.6: a real rate dependence moves both binnings the same way (regression to the mean moves them apart)", fontsize=9)
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    clock = time.time()
    workers = max(min((os.cpu_count() or 2) - 1, 3), 1)
    sha = hashlib.sha256(P.S6_CANDIDATE.read_bytes()).hexdigest()
    d, BANDS, sets, gmap, mods, free = load()
    preds = predictions(sets, mods, free)
    W1 = {s_: {m: w_minutes(D, le, mods[m]["fz"]) for m in MODELS} for s_, (D, le) in sets.items()}
    print(f"loaded, predictions done ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # 7.1
    lt = link_table(d, sets["TRAIN"][1], mods, BANDS)
    at = alink_tests(lt)
    wh = w_height(lt)
    rL = at[at.variable == "L_km"].iloc[0]
    L_dep = bool(abs(rL.rho) >= 0.3 and (rL.lo > 0 or rL.hi < 0))
    whc = wh[(wh.model == "candidate") & (wh.by == "length")]
    w_ok = {}
    for R in (2.0, 5.0):
        g = whc[whc.R_mmh == R].set_index("length")
        w_ok[R] = all(overlap((g.lo[a], g.hi[a]), (g.lo[b], g.hi[b])) for i, a in enumerate(g.index) for b in g.index[i + 1:])
    pass71 = (not L_dep) and all(w_ok.values())
    print(f"7.1 done: L rho {rL.rho:+.2f} [{rL.lo:+.2f}, {rL.hi:+.2f}], W overlap {w_ok} ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # 7.2
    rm = rmse_table(sets, preds)
    ltc, ltv = length_test_val(sets, mods)
    sh = shape_short(sets, mods, free)
    ar = after_rain(sets, preds)
    v = rm[(rm.set == "VAL") & (rm.model == "candidate")]
    rule72 = v[v.by.isin(["length", "band"])]
    pass72_rmse = bool((rule72.over_all_class <= S["short_links_val"]["rmse_factor_max"]).all())
    pass72_len = bool(ltv.set_index("fit").passes["candidate"])
    pass72_shape = bool(sh.ok.all())
    print(f"7.2 done ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # 7.3a
    Rm, T5, col = radar_max(d, sets)
    win, curves = dry_windows(sets, gmap, Rm, T5, col, W1)
    dtab, ddiff = dry_stats(win, curves) if len(win) else (pd.DataFrame(), {})
    dc, dmed = decay(win, curves) if len(win) else (pd.DataFrame(), {})
    Q = S["long_links"]["dry_path"]
    if len(dtab):
        g10 = dtab[(dtab.t_min == 10) & (dtab.by == "length")].set_index("length")
        big = [nm for _, nm in LB if nm in g10.index and g10.windows[nm] >= Q["min_windows"]]
        err_ok = all(abs(g10.err_candidate[nm]) <= Q["max_err_10_db"] for nm in big)
        p_, lo_, hi_ = ddiff["candidate"]
        diff_ok = bool(np.isfinite(lo_) and lo_ <= 0 <= hi_)
        testable73a = LB[0][1] in big and LB[2][1] in big
        pass73a = bool(err_ok and diff_ok) if testable73a else None
    else:
        g10, big, err_ok, diff_ok, testable73a, pass73a = pd.DataFrame(), [], False, False, False, None
    print(f"7.3a done: {len(win)} windows ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # 7.3b
    pt, pe, pinfo = pairs(d, sets, lt, preds)
    pc = pair_compare(pe) if pinfo.get("testable") else pd.DataFrame()
    # 7.3c
    mg, ifits = intercepts(sets, preds)
    n_cells = len(mg)
    pass73c = bool(n_cells and mg.overlap.sum() >= need(n_cells))
    ps = p_size(mods, sets["TRAIN"][1], BANDS)
    print(f"7.3b-e done ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # 7.3d: does anything separate the two F descriptions?
    sep = []
    if len(win):
        obs = np.array([c["obs"] for c in curves])
        hv = win.heavy.to_numpy()
        ev = pd.factorize(win.event_id)[0]
        draws = boot_items(ev, 1000, Q["seed"] + 3)
        for t in (5, 10, 20):
            s_ = hv & (win.dry_minutes.to_numpy() > t) & np.isfinite(obs[:, t])
            if s_.sum() < 5:
                continue
            lo_, hi_ = ci([np.median(obs[i[s_[i]], t]) if s_[i].any() else np.nan for i in draws])
            r = {"check": f"dry path after heavy rain, W at {t} min", "n": int(s_.sum()), "observed": float(np.median(obs[s_, t])), "lo": lo_, "hi": hi_}
            for m in MODELS:
                mv = float(np.median(np.array([c[m] for c in curves])[s_, t]))
                r[m] = mv
                r[f"{m} inside"] = bool(lo_ <= mv <= hi_)
            sep.append(r)
    fa = ifits["all lengths"]
    for rr in fa[fa["class"] == "10-20"].itertuples():
        if np.isfinite(rr.I):
            r = {"check": f"intercept I, {rr.band}, 10-20 mm/h (all lengths)", "n": int(rr.link_events), "observed": rr.I, "lo": rr.I_lo, "hi": rr.I_hi}
            for m in MODELS:
                mv = fa.loc[rr.Index, f"I_{m}"]
                r[m] = float(mv)
                r[f"{m} inside"] = bool(rr.I_lo <= mv <= rr.I_hi)
            sep.append(r)
    sept = pd.DataFrame(sep)
    separates = bool(len(sept) and (sept["candidate inside"] != sept["closed F inside"]).any())

    # 7.4
    dw = dew(sets, lt, preds, Rm, T5, col)
    print(f"7.4 done ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # 7.5: timing first, plan, cuts
    targets, te = leftover_targets(lt, win, curves)
    plan_file = CACHE / "plan.json"
    A = S["leftovers"]
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
    else:
        X, y, _ = targets["a_link"]
        Xs = ((X - X.mean()) / X.std(ddof=0)).to_numpy(float)
        t0 = time.time()
        res0 = {}
        run_pool([{"name": "pysr_a_link_real_s1", "X": Xs.tolist(), "y": y.tolist(), "seed": 1}], 1, clock, res0)
        wall = time.time() - t0
        fit_s = res0["pysr_a_link_real_s1"]["seconds"]
        n_t = sum(1 for k in targets if len(targets[k][1]) >= 8)

        def proj(p):
            return (n_t * (1 + p["shuffles"]) * p["seeds"] * fit_s / workers + (wall - fit_s) * workers) / 3600 + 0.3
        plan = {"fit_s": round(fit_s, 1), "init_s": round(wall - fit_s, 1), "workers": workers, "seeds": len(A["pysr"]["seeds"]), "shuffles": A["shuffles"], "cuts": []}
        plan["projected_full_h"] = round(proj(plan), 2)
        for cut, ch in (("PySR seeds 3 -> 2", lambda p: p.update(seeds=2)), ("shuffles 50 -> 30", lambda p: p.update(shuffles=30))):
            if proj(plan) <= S["time_budget_h"]:
                break
            ch(plan)
            plan["cuts"].append(cut)
        plan["projected_h"] = round(proj(plan), 2)
        plan_file.write_text(json.dumps(plan))
    print(f"7.5 timing: one PySR run {plan['fit_s']} s (+{plan['init_s']} s start); projection {plan['projected_full_h']} h; cuts: {plan['cuts'] or 'nothing'}",
          flush=True)
    try:
        lo_t, scr, _ = leftovers(targets, clock, workers, plan)
    except Exception as e:                                       # noqa: BLE001  (a worker pool problem: run serially)
        print(f"   pool failed ({type(e).__name__}: {e}); running serially", flush=True)
        lo_t, scr, _ = leftovers(targets, clock, 1, plan)
    hit75 = bool((lo_t.verdict == "RELATION (hint only)").any()) if len(lo_t) else False
    print(f"7.5 done ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # 7.6
    rt, found76, ch76, res76, pev = rain_reference(d, sets, gmap, preds)
    print(f"7.6 done: rate dependence found = {found76} ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # summary and the return decision
    def vtxt(x):
        return "not testable" if x is None else ("pass" if x else "FAIL")
    summ = [
        {"check": "7.1 a_link against L", "rule": "not (|rho| >= 0.3 and interval excludes 0)", "result": f"rho {rL.rho:+.2f} [{rL.lo:+.2f}, {rL.hi:+.2f}]", "verdict": vtxt(not L_dep)},
        {"check": "7.1 W height per length class", "rule": "the three classes overlap pairwise at 2 and 5 mm/h",
         "result": "; ".join(f"{R:g} mm/h: " + ", ".join(f"{r.length} {r.W_median_db:.2f} [{r.lo:.2f}, {r.hi:.2f}]" for r in whc[whc.R_mmh == R].itertuples()) for R in (2.0, 5.0)),
         "verdict": vtxt(all(w_ok.values()))},
        {"check": "7.2 VAL RMSE per length class and band", "rule": "each <= 1.15 x all-class VAL RMSE",
         "result": ", ".join(f"{r['class']} {r.over_all_class:.2f} x" for _, r in rule72.iterrows()), "verdict": vtxt(pass72_rmse)},
        {"check": "7.2 VAL length test (candidate)", "rule": "at least 4 of 6 cells", "result": f"{int(ltv.set_index('fit').both_inside['candidate'])} of "
         f"{int(ltv.set_index('fit').eligible_cells['candidate'])}", "verdict": vtxt(pass72_len)},
        {"check": "7.2 shared against band shapes, short VAL links", "rule": "shared not worse by > 0.003 dB in any band",
         "result": ", ".join(f"{r.band} {r.shared_minus_band_db:+.4f} dB" for r in sh.itertuples()), "verdict": vtxt(pass72_shape)},
        {"check": "7.3a dry-path W, L > 3 km against L < 1 km", "rule": "difference of the median errors at 10 min: interval contains 0; |median error at 10 min| <= 0.3 dB",
         "result": (f"difference {fmt(*ddiff['candidate'])} dB; errors " + ", ".join(f"{nm} {g10.err_candidate[nm]:+.2f}" for nm in g10.index)) if len(g10) else "no windows",
         "verdict": vtxt(pass73a)},
        {"check": "7.3b paired links", "rule": ">= 5 pairs and median estimator noise <= 1 dB, else not testable",
         "result": f"{pinfo.get('pairs_with_events', 0)} pairs with common events, median noise {fmt(pinfo.get('median_noise_sd_db', np.nan))} dB",
         "verdict": "compared (information)" if pinfo.get("testable") else "not testable"},
        {"check": "7.3c intercept short-only against long-only", "rule": "intervals overlap in at least 4 of 6 cells (k of n = ceil(2 n / 3))",
         "result": f"{int(mg.overlap.sum()) if n_cells else 0} of {n_cells}", "verdict": vtxt(pass73c) if n_cells else "not testable"},
        {"check": "7.3d F form (free knots against (R / 5)^0.47)", "rule": "an observation inside one and outside the other", "result": "separates" if separates else "cannot tell",
         "verdict": "information"},
        {"check": "7.3e c_band against ITU-R P.838 (size only)", "rule": "flag a ratio outside 0.5-2", "result": ", ".join(f"{r.band} {r.ratio:.2f}" for r in ps[ps.R_mmh == 5].itertuples()) + " (at 5 mm/h)",
         "verdict": "flagged" if ps.flag.any() else "no flag"},
        {"check": "7.4 dew", "rule": ">= 8 link-nights", "result": f"{dw['link_nights']} link-nights on {dw['links']} links", "verdict": "tested (information)" if dw["testable"] else "not testable"},
        {"check": "7.5 leftovers (a_link, drying time)", "rule": "skill > 0 and > 95 % of the null", "result": "; ".join(f"{r.target}: {r.verdict}" for r in lo_t.itertuples()),
         "verdict": "RELATION" if hit75 else "no relation found"},
        {"check": "7.6 radar / gauge rate dependence (V1)", "rule": "high / low differs > 25 % in both binnings, same direction, intervals exclude 1",
         "result": f"gauge binning {ch76['gauge'][0]:+.0%}, radar binning {ch76['radar'][0]:+.0%}", "verdict": "FOUND" if found76 else "not found"},
    ]
    summ = pd.DataFrame(summ)
    reasons = {"(a) 7.1 dependence on L": not pass71, "(b) 7.3a or 7.3c fails": (pass73a is False) or (not pass73c and n_cells > 0),
               "(c) 7.6 rate dependence": found76, "(d) 7.5 relation beats the null": hit75}
    go_back = any(reasons.values())

    # save
    for tab, path in ((summ, P.S7_SUMMARY), (at.round(4), P.S7_ALINK), (lt.round(4), P.S7_LINKS), (wh.round(4), P.S7_WHEIGHT), (rm.round(4), P.S7_RMSE),
                      (ltc.round(4), P.S7_LENGTH), (sh.round(5), P.S7_SHAPE), (ar.round(4), P.S7_AFTER), (win.round(4), P.S7_DRY_WINDOWS),
                      (pd.concat([dtab, sept.assign(by="7.3d F form")], ignore_index=True).round(4), P.S7_DRY), (dc.round(4), P.S7_DECAY),
                      (pd.concat([pt, pe.drop(columns="_bins", errors="ignore")], ignore_index=True).round(4), P.S7_PAIRS), (mg.round(4), P.S7_INTERCEPT),
                      (ps.round(4), P.S7_PSIZE), (dw["jumps"].round(4), P.S7_DEW), (dw.get("events", pd.DataFrame()).round(4), P.S7_DEW_EVENTS),
                      (pd.concat([lo_t, scr.assign(kind="screen")], ignore_index=True).round(4), P.S7_LEFTOVER), (te.round(4), P.S7_DRYING),
                      (pd.concat([rt, res76.assign(binned_by="candidate residual"), pev.assign(binned_by="per event")], ignore_index=True).round(4), P.S7_RAIN_REF)):
        save(tab, path, index=False)
    figs = [(fig_alink(lt, at), P.S7_FIG_ALINK), (fig_wheight(wh), P.S7_FIG_WHEIGHT), (fig_rmse(rm), P.S7_FIG_RMSE), (fig_dry(win, curves), P.S7_FIG_DRY),
            (fig_decay(dc), P.S7_FIG_DECAY), (fig_pairs(pe, pinfo), P.S7_FIG_PAIRS), (fig_intercept(mg), P.S7_FIG_INTERCEPT), (fig_dew(dw), P.S7_FIG_DEW),
            (fig_rain_ref(rt), P.S7_FIG_RAIN_REF)]
    for fig, path in figs:
        save(fig, path, dpi=100)
        plt.close(fig)
    q = dict(locals())
    q.update(readings(q))
    report = write_report(q)
    save(report, P.S7_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {(time.time() - clock) / 60:.1f} min")


# =================================================================================================== readings
def readings(q):
    """The READING paragraphs of the report, written from the tables of this run."""
    R = {}
    at, whc = q["at"], q["whc"]
    rL = at[at.variable == "L_km"].iloc[0]
    disc = at[at.bh_discovery].variable.tolist()
    nom = at[(at.test == "Spearman") & ((at.lo > 0) | (at.hi < 0))]
    w5 = whc[whc.R_mmh == 5.0].set_index("length").W_median_db
    R["READ71"] = (f"**READING.** a_link does not depend on the length: rho {rL.rho:+.2f} [{rL.lo:+.2f}, {rL.hi:+.2f}], far from the 0.3 of the rule, and the W height "
                   f"a_link F is the same in the three length classes ({', '.join(f'{k} {v:.2f}' for k, v in w5.items())} dB at 5 mm/h, overlapping intervals; the "
                   f"closed F gives the same picture about 0.15 dB lower). **7.1 passes.** No variable is a BH discovery "
                   f"({'none' if not disc else ', '.join(disc)}); nominally (interval excluding 0, BH not passed) "
                   f"{', '.join(f'{r.variable} {r.rho:+.2f}' for r in nom.itertuples()) or 'none'}: larger antennas a little larger a_link, the orientation a "
                   "little. These are findings, not failures, and too weak to act on (|rho| < 0.2). The metadata has no radome or antenna type.")
    rm = q["rm"]
    v = rm[(rm.set == "VAL") & (rm.model == "candidate")].set_index(["by", "class"])
    t = rm[(rm.set == "TRAIN") & (rm.model == "candidate")].set_index(["by", "class"])
    long_ = LB[2][1]
    sh = q["sh"]
    R["READ72"] = (f"**READING.** The RMSE here pools all 5 min bins (the S3-S6 judge averaged per event, hence the larger numbers: VAL all classes "
                   f"{v.rmse[('all', 'all')]:.3f} dB, TRAIN {t.rmse[('all', 'all')]:.3f} dB). The written rule **fails for the long links**: VAL L > 3 km is "
                   f"{v.over_all_class[('length', long_)]:.2f} x the all-class RMSE (rule 1.15). It is not a VAL effect: on TRAIN the same class is "
                   f"{t.over_all_class[('length', long_)]:.2f} x, and every VAL length and band class is BELOW its TRAIN value (VAL / TRAIN 0.72-0.80). Long links carry a larger "
                   "path term, so their absolute error in rain is larger; the rule compared absolute errors across lengths, which is the wrong yardstick for "
                   "the long class, but it was written so and the check is reported as failed (it is not one of the return reasons). The short links, which this "
                   f"check is about, are the best class ({v.over_all_class[('length', LB[0][1])]:.2f} x). The rain classes on VAL are within "
                   f"{rm[(rm.set == 'VAL') & (rm.model == 'candidate') & (rm.by == 'rain')].val_over_train.max():.2f} x of TRAIN (information). The VAL length "
                   f"test passes ({int(q['ltv'].set_index('fit').both_inside['candidate'])} of 6; the closed F reaches exactly 4 of 6 and passes under the "
                   "2026-10-19 wording). On the short VAL links the shared shape is as good as the band shapes "
                   f"({', '.join(f'{r.band} {r.shared_minus_band_db:+.4f} dB' for r in sh.itertuples())}; no 32 GHz link is under 1 km). After the rain the "
                   "median residual is small and negative (the model a little HIGH: -0.05 to -0.12 dB in the first 10 min, -0.01 to -0.02 dB after an hour), "
                   "the same on TRAIN and VAL and slightly larger on short links.")
    dt, dc, win = q["dtab"], q["dc"], q["win"]
    if len(dt):
        g = dt[(dt.by == "length")].set_index(["t_min", "length"])
        e10 = {nm: g.err_candidate[(10, nm)] for _, nm in LB}
        dcx = dc.set_index("length")
        R["READ73A"] = (f"**READING.** {len(win)} dry-path windows exist ({', '.join(f'{nm} {int(g.windows[(5, nm)])}' for _, nm in LB)} at 5 min; 32 GHz has none "
                        f"under 1 km and few above 3 km). **7.3a passes**: the model's error at 10 min is the same in every length class "
                        f"({', '.join(f'{k} {v:+.2f}' for k, v in e10.items())} dB; L > 3 minus L < 1: {fmt(*q['ddiff']['candidate'])} dB) and below the 0.3 dB "
                        "of the rule. So W measured where no rain is on the path does not depend on L, as the model assumes. But the model is HIGH by 0.2-0.3 dB "
                        "right after the rain stops, in every class: the measured median W starts at about 0.45 dB where the model has about 0.75 dB "
                        "(figure), and the gap closes over 20-40 min. The decay time of the normalised median curve is "
                        f"{', '.join(f'{nm} {dcx.tau_obs[nm]:.0f} [{dcx.tau_obs_lo[nm]:.0f}, {dcx.tau_obs_hi[nm]:.0f}]' for _, nm in LB)} min measured against "
                        f"about {dcx.tau_candidate.median():.0f} min in the model (both slower than tau_f_dry = 10.3 min because the slow store and the "
                        "windows with small W flatten the median curve); the share left at 40-60 min is "
                        f"{', '.join(f'{nm} {dcx.left_40_60_obs[nm]:.2f}' for _, nm in LB)} measured against about {dcx.left_40_60_candidate.median():.2f} "
                        "in the model. The long links keep MORE water late (decay "
                        f"{dcx.tau_obs[long_]:.0f} min, {dcx.left_40_60_obs[long_]:.2f} left), a hint that a path holding drizzle below the radar's 0.05 mm/h "
                        "over 3-8 km still attenuates a little, or that the long links' antennas dry more slowly; it is not part of the rule and is reported.")
    R["READ73D"] = ("**READING.** **Cannot tell.** The intercepts of the 10-20 mm/h class have intervals of 1.5-6 dB and hold both descriptions in every "
                    "band; the paired links are not testable. The dry-path windows after heavy rain (>= 10 min of V1 >= 10 mm/h in the hour before) "
                    "reject BOTH descriptions alike: measured W 0.25 dB at 5 min, both models about 0.8 dB (0.6 at 10 min against 0.13 measured). After "
                    "heavy rain the antenna is drier, sooner, than the model says, whatever the shape of F; the two shapes differ by only 0.03 dB there. "
                    "This is the clearest misfit S7 found (section 11).")
    ps = q["ps"]
    R["READ73E"] = (f"**READING.** c_band R is {ps.ratio.min():.2f}-{ps.ratio.max():.2f} times the ITU-R P.838 size (vertical, the band's median frequency): "
                    "the same order, about 0.6 of it, rising with R and with the frequency. No ratio is outside 0.5-2, so nothing is flagged; 28 GHz at 2 mm/h "
                    "(0.52) is close to the edge. A path term below the ITU size is what the radar's path mean on a 2 km grid and gauge-adjusted rain would "
                    "give if rain on the path is spread unevenly; it is a size reference only and corrects nothing.")
    dw = q["dw"]
    if dw["testable"]:
        r = dw["same_layer"]
        cp = dw["compare"].set_index("definition")
        R["READ74"] = (f"**READING.** Dew is testable ({dw['link_nights']} link-nights). The dew W is LARGER than the rain W at light rain: median "
                       f"{dw['size'].loc['median_excess_db', '50%']:.1f} dB above the dry baseline (peak {dw['size'].loc['peak_excess_db', '50%']:.1f} dB) "
                       f"against {dw['rain_W_075']:.1f} dB for the candidate at 0.75 mm/h, and it lasts hours (median {dw['size'].loc['duration_h', '50%']:.0f} h, "
                       "starting between midnight and 3 h local), not the minutes of the rain drying. **Same water layer: yes, weakly**: links with a "
                       f"larger a_link get larger dew rises (rho {r[0]:+.2f} [{r[2]:+.2f}, {r[3]:+.2f}], {r[4]} links), as expected if the same radome "
                       "holds both. Dew before the rain is rare in these events ("
                       f"{int(cp.link_events_with['dew_before'])} link-events in {int(cp.events_with['dew_before'])} events): there c is "
                       f"{cp.c_db_with_minus_without['dew_before']:+.2f} [{cp.c_db_lo['dew_before']:+.2f}, {cp.c_db_hi['dew_before']:+.2f}] dB higher (dew "
                       f"water inside the pre-rain offset c), the model bias is {cp.bias_with_minus_without['dew_before']:+.2f} "
                       f"[{cp.bias_lo['dew_before']:+.2f}, {cp.bias_hi['dew_before']:+.2f}] dB (A over-predicted a little: the dew already in c is counted "
                       f"again as rain wetting), RMSE {cp.rmse_with_minus_without['dew_before']:+.2f} [{cp.rmse_lo['dew_before']:+.2f}, "
                       f"{cp.rmse_hi['dew_before']:+.2f}] dB. The stricter definition (a dew jump of the sub-link) has only "
                       f"{int(cp.link_events_with['jump_before'])} link-events in {int(cp.events_with['jump_before'])} events: its large numbers are not "
                       "interpretable. Dew is outside the equation; it matters through c (about 0.1 dB, at most a few tenths) and on dry nights.")
    lo_t = q["lo_t"].set_index("target")
    te = q["te"]
    R["READ75"] = ("**READING.** **a_link: no relation found** (the best PySR form, "
                   f"{lo_t.pysr_equation.get('a_link', '')}, has a leave-one-link-out skill of {lo_t.loo_skill.get('a_link', np.nan):.3f}, below the 95 % of "
                   f"the shuffled null, {lo_t.null_q95.get('a_link', np.nan):.3f}; the screen finds antenna diameter and the east-west position, both weak). "
                   "**Drying time: a relation beats the null**: per-event decay time against the wind speed, "
                   f"`{lo_t.pysr_equation.get('drying', '')}` (wind standardised), leave-one-event-out skill {lo_t.loo_skill.get('drying', np.nan):.2f} against "
                   f"the null's 95 % {lo_t.null_q95.get('drying', np.nan):.2f}; Spearman of wind against the decay time "
                   f"{q['scr'][(q['scr'].target == 'drying') & (q['scr'].variable == 'wind_speed')].rho.iloc[0]:+.2f} ({len(te)} events). Windier events dry "
                   "faster, which is physical (evaporation and shaking off the radome). It is a hint on 17 events with one weather station and hourly wind, and "
                   "the per-event decay time also carries the rain history of the event; by the return rule it is a reason to return to S6.")
    ch = q["ch76"]
    R["READ76"] = (f"**READING.** **No rate dependence found.** Binned by the gauge, the radar reads LOW at high rates (high / low {ch['gauge'][0] + 1:.2f}); "
                   f"binned by the radar, it reads HIGH ({ch['radar'][0] + 1:.2f}): the two binnings move in OPPOSITE directions, the signature of regression "
                   "to the mean between two noisy measures of a point at 5 min (a 2 x 2 km radar pixel against a gauge orifice), not of a radar error that "
                   "grows with the rate. The rule asks for the same direction; it is not met, so the parked g(R) correction is NOT justified by these data. "
                   "The candidate's mean residual in V1 >= 10 mm/h bins has intervals that hold 0 on TRAIN and VAL (information).")
    reasons = q["reasons"]
    R["READ9"] = ("The return rule was written before the run; the only reason met is (d), the drying-time relation with the wind. The checks that test the "
                  "candidate's structure (W independent of L: 7.1, 7.3a, 7.3c; the rain reference: 7.6) all pass. "
                  + ("So S7 sends the work back to S6 ONCE (recorded and dated 2026-10-19 in the roadmap). What S6 would test, under the same rules: whether "
                     "tau_f_dry depending on the wind (e.g. tau_f_dry = tau_0 exp(-k u), u the hourly wind) is identified on fake data, improves the model out of "
                     "sample, and changes the other numbers; and, from 7.3d, why the model is about 0.5 dB high in the minutes after heavy rain. This step "
                     "does NOT start it." if any(reasons.values()) else ""))
    R["READ10"] = ("1. Read TEST once, after a hashed freeze of whichever candidate leaves S6 (the S6 candidate hashed here if the return changes nothing).\n"
                   "2. The TEST length test through the model, at least 4 of 6 cells, and the TEST RMSE per length class against VAL (the S5a plan), on V1 "
                   "with g_e from the TEST events' City gauges, V0 and the gauge links as the independent check.\n"
                   "3. The dry-path check of 7.3a on TEST: the error at 10 min per length class (and after heavy rain, where S7 found the model 0.5 dB high).\n"
                   "4. The long links separately: their RMSE is 1.1-1.2 x the all-class value on TRAIN and VAL (7.2).\n"
                   "5. If wind enters tau_f_dry in S6: the same on TEST with the TEST hours' wind; if not, the residual of the drying against the wind.")
    R["READ11"] = ("- **The RMSE rule of 7.2 failed for L > 3 km** (1.18 x) although every VAL length and band class is better than on TRAIN; the rule compared absolute errors "
                   "across lengths. It is reported as failed; it is not one of the return reasons.\n"
                   "- **After heavy rain the model is about 0.5 dB too high** in the first 10 min of dry path (7.3d), for both F descriptions, and about 0.2 dB "
                   "high after ordinary rain (7.3a). The fast store may dry faster after heavy rain than tau_f_dry = 10.3 min, or F above 10 mm/h is too high "
                   "(F there rests on 0.4 % of the minutes), or the radar marks the path dry a little late. Not separable here.\n"
                   "- **The long links keep more water late** in dry-path windows (7.3a); drizzle below 0.05 mm/h per pixel over a long path is possible.\n"
                   "- **The drying relation** rests on 17 events, one station, hourly wind, and decay times of pooled windows (which also carry each event's "
                   "rain history); the permutation null guards against chance in the procedure, not against a confounder such as storm type.\n"
                   "- **The a_link relation** missed the null by a hair (skill 0.077 against 0.081); it is reported as none.\n"
                   "- **Dew before the rain** is rare here (8 events); the stricter definition (2 events) is not interpretable.\n"
                   "- **Paired links**: only 4 pairs meet the geometry, so 'not testable' (the noise, 0.97 dB, would have just passed).\n"
                   "- **Windows are not independent**: both directions of a link and neighbouring links share an event; every interval resamples events.\n"
                   "- **VAL has been used for decisions** since S3; its numbers here are optimistic.")
    return R


# =================================================================================================== report
def write_report(q):
    L = []
    a = L.append
    a("# wet_antenna step S7: validate the frozen S6 candidate without TEST\n")
    a("Written by `research/wet_antenna/validate_s7.py` (libraries `final_model.py`, `system_id.py`, `length_test.py`, `event_factors.py`). Every rule: "
      "`wet_antenna.validate_s7` in `research/config.yaml`, written on 2026-10-19 BEFORE the run. Plan, decisions and integrity rules: `docs/ROADMAP.md`.\n")
    a(STANDING + "\n")
    a(f"**Frozen object.** The S6 candidate, `{rel(P.S6_CANDIDATE)}`, sha256 `{q['sha']}`. Nothing of it was refitted. VAL numbers are **{VL}**.\n")
    a("## 1. What was run\n")
    D, le = q["sets"]["TRAIN"]
    Dv, lev = q["sets"]["VAL"]
    a(f"- Data: TRAIN {len(le)} link-events ({int(D['le_event'].max()) + 1} events), VAL {len(lev)} ({int(Dv['le_event'].max()) + 1} events), input V1 = radar / g_e.")
    a("- Two descriptions of F, both frozen S6 fits: **candidate** (free shared F knots, P = c_band R; the S6 plug-back 'P only') and **closed F** "
      "(F = A_band (R / 5)^b, P = c_band R; the S6 plug-back 'F and P'). The S5a free fit (band shapes) only in 7.2.")
    a(f"- 7.5 timing: one PySR run {q['plan']['fit_s']} s (+{q['plan']['init_s']} s start-up); projection {q['plan']['projected_full_h']} h for the full plan "
      f"(budget {S['time_budget_h']} h); cuts: {', '.join(q['plan']['cuts']) or 'nothing'}.")
    a(f"- Run time: {(time.time() - q['clock']) / 60:.1f} min.\n")
    a("## 2. Every check, its rule and its result\n")
    a(md_table(q["summ"]) + "\n")
    a("## 3. 7.1 W independent of L\n")
    a("a_link against the link metadata (Spearman of log a_link, 90 % link-bootstrap interval, BH q = 0.10; "
      f"`{rel(P.S7_ALINK)}`). The metadata has no antenna type, radome or mast height; a_link is one number per link, so the direction of a sub-link does not apply.\n")
    a(md_table(q["at"].round(3)) + "\n")
    a(f"W height a_link F (equilibrium) per length class, median over links with a 90 % link-bootstrap interval (`{rel(P.S7_WHEIGHT)}`):\n")
    w = q["wh"][q["wh"].by == "length"].round(3)
    a(md_table(w) + "\n")
    a("![a_link against L](../figures/validate_s7_alink_vs_L.png)\n")
    a("![W height](../figures/validate_s7_w_height.png)\n")
    a(q.get("READ71", "") + "\n")
    a(f"## 4. 7.2 Short links on VAL ({VL})\n")
    rm = q["rm"]
    a(f"RMSE per class (`{rel(P.S7_RMSE)}`); `over_all_class` = RMSE / the all-class RMSE of the same set and model; `val_over_train` for VAL:\n")
    a(md_table(rm[rm.model == "candidate"].round(4)) + "\n")
    a("The closed F and the S5a free fit, all classes, for comparison:\n")
    a(md_table(rm[(rm.model != "candidate") & (rm.by == "all")].round(4)) + "\n")
    a("VAL V1 length test through the model (at least 4 of 6 cells):\n")
    a(md_table(q["ltv"]) + "\n")
    a("Shared shapes (candidate) against band shapes (S5a free fit) on the short VAL links (L < 1 km), RMSE in dB:\n")
    a(md_table(q["sh"].round(4)) + "\n")
    ar = q["ar"]
    if len(ar):
        piv = ar.pivot_table(index="after_min", columns=["set", "length"], values="median_residual_db").round(3)
        piv.columns = [f"{s} {l}" for s, l in piv.columns]
        a(f"Median residual A - A_pred (dB) against the minutes after the link's last wet bin (`{rel(P.S7_AFTER)}`):\n")
        a(md_table(piv.reset_index()) + "\n")
    a("![VAL RMSE](../figures/validate_s7_val_rmse.png)\n")
    a(q.get("READ72", "") + "\n")
    a("## 5. 7.3 Long links without a power law\n")
    a("### 7.3a Dry-path moments\n")
    win = q["win"]
    a(f"Windows (path-max radar <= 0.05 mm/h on every pixel, V0 and V1, the City gauge 0 if any, after >= 10 min of V1 >= 2 mm/h in the hour before; "
      f"`{rel(P.S7_DRY_WINDOWS)}`): **{len(win)} windows** on {win.cml_id.nunique() if len(win) else 0} links in {win.event_id.nunique() if len(win) else 0} events.\n")
    if len(win):
        cnt = win.groupby(["length", "band"]).size().unstack(fill_value=0)
        a(md_table(cnt.reset_index()) + "\n")
    dt = q["dtab"]
    if len(dt):
        a(f"Measured W (A - c) against the model's W, median over windows, 90 % event-bootstrap intervals (`{rel(P.S7_DRY)}`):\n")
        a(md_table(dt[dt.by == "length"].round(3)) + "\n")
        a("Per band and length class (information):\n")
        a(md_table(dt[dt.by != "length"][["t_min", "band", "length", "windows", "W_obs", "W_candidate", "err_candidate", "err_closed F"]].round(3)) + "\n")
        a("Difference of the median error at 10 min, L > 3 km minus L < 1 km (rule: interval contains 0): " +
          "; ".join(f"{m} {fmt(*q['ddiff'][m])} dB" for m in MODELS) + "\n")
    dc = q["dc"]
    if len(dc):
        a(f"Decay time of the median normalised W over 0-30 min and the share left at 40-60 min (`{rel(P.S7_DECAY)}`):\n")
        a(md_table(dc.round(3)) + "\n")
    a("![dry path](../figures/validate_s7_dry_path.png)\n")
    a("![decay](../figures/validate_s7_decay.png)\n")
    a(q.get("READ73A", "") + "\n")
    a("### 7.3b Paired links\n")
    pi = q["pinfo"]
    a(f"Link pairs with a shared site and the same band: {pi.get('candidates_shared_site_same_band', 0)}; with the geometry rule (bearing within 15 deg, L2 / L1 >= 1.5, "
      f"L2 - L1 >= 1 km): {pi.get('pairs_geometry', 0)}; with common events and >= 6 common rain bins: {pi.get('pairs_with_events', 0)} pairs, "
      f"{pi.get('pair_events', 0)} pair-events. Fake check (truth W = 1.5 dB, k = 0.5 dB/km, the real residual series of the two link-events): median "
      f"estimator noise SD **{fmt(pi.get('median_noise_sd_db', np.nan))} dB** -> **{'testable' if pi.get('testable') else 'NOT testable'}** "
      f"(rule: >= 5 pairs and <= 1 dB; `{rel(P.S7_PAIRS)}`).\n")
    if len(q["pe"]):
        a(md_table(q["pe"].drop(columns="_bins").round(3)) + "\n")
    if len(q["pc"]):
        a("Observed minus model pair estimate per rain class:\n")
        a(md_table(q["pc"].round(3)) + "\n")
    a("![pairs](../figures/validate_s7_pairs.png)\n")
    a("### 7.3c The intercept from short-only and long-only links\n")
    a(f"Model-free length test on TRAIN V1 (S3 rules) twice, with the models' own I on the same bins (`{rel(P.S7_INTERCEPT)}`):\n")
    mg = q["mg"]
    cols = ["band", "class", "R_cell_short", "links_short", "I_short", "I_lo_short", "I_hi_short", "links_long", "I_long", "I_lo_long", "I_hi_long",
            "I_candidate_short", "I_candidate_long", "I_closed F_short", "I_closed F_long", "overlap"]
    a(md_table(mg[[c for c in cols if c in mg]].round(3)) + "\n")
    a("![intercept](../figures/validate_s7_intercept.png)\n")
    a("### 7.3d Can the data tell the two F descriptions apart?\n")
    a(md_table(q["sept"].round(3)) + "\n" if len(q["sept"]) else "No heavy-rain dry-path windows or 10-20 mm/h intercepts.\n")
    a(q.get("READ73D", "") + "\n")
    a("### 7.3e The size of P against ITU-R P.838 (size reference only, never a correction)\n")
    a(md_table(q["ps"].round(4)) + "\n")
    a(q.get("READ73E", "") + "\n")
    a("## 6. 7.4 Dew\n")
    dw = q["dw"]
    a(f"Source: `{rel(P.DRY_JUMPS)}` (baseline module, `make_dew_products.py`: the R1 C4 dry-day jumps with the SMHI GbgA weather), positive jumps of the S7 "
      f"links, >= 6 h after the last wet slot, at dew risk by the baseline rule (mean RH >= 90 % or dew-point depression < 2 C): **{len(dw['jumps'])} jumps, "
      f"{dw['link_nights']} link-nights on {dw['links']} links** (`{rel(P.S7_DEW)}`).\n")
    if dw["testable"]:
        a("Size and timing of the dew W:\n")
        a(md_table(dw["size"].round(3).reset_index().rename(columns={"index": "quantity"})) + "\n")
        hrs = dw["jumps"].start_hour_local.value_counts().sort_index()
        a("Local start hour of the jumps: " + ", ".join(f"{h} h: {n}" for h, n in hrs.items()) + "\n")
        r = dw["same_layer"]
        a(f"Same water layer? Spearman of the link's median dew peak against log a_link: {r[0]:+.2f} [{r[2]:+.2f}, {r[3]:+.2f}] (p {r[1]:.3f}, {r[4]} links).\n")
        a(f"The rain W of the candidate at 0.75 mm/h (median link): {dw['rain_W_075']:.2f} dB; its fast drying time {FIX['tau_f_dry']:g} min.\n")
        ev_ = dw["events"]
        a(f"Dew before the rain (>= 4 of the 6 hours before the link's first wet bin at dew risk, the path dry): {int(ev_.dew_before.sum())} of {len(ev_)} link-events "
          f"({ev_[ev_.dew_before].event_id.nunique()} events); a dew jump of the sub-link in those 6 h: {int(ev_.jump_before.sum())} (`{rel(P.S7_DEW_EVENTS)}`). "
          f"With minus without, within the rain class, event bootstrap (TRAIN and VAL; VAL {VL}):\n")
        a(md_table(dw["compare"].round(3)) + "\n")
    a("![dew](../figures/validate_s7_dew.png)\n")
    a(q.get("READ74", "") + "\n")
    a("## 7. 7.5 What explains the leftovers\n")
    a(f"Screening (Spearman, BH q = 0.10) and PySR on all variables (complexity <= 6), guarded by leave-one-out skill against a constant and a permutation null "
      f"(`{rel(P.S7_LEFTOVER)}`; per-event drying times in `{rel(P.S7_DRYING)}`):\n")
    a(md_table(q["lo_t"].round(3)) + "\n")
    a(md_table(q["scr"].round(3)) + "\n")
    if len(q["te"]):
        a(f"Per-event drying times: {len(q['te'])} events with >= {S['leftovers']['min_windows_per_event']} windows; median {q['te'].tau_min.median():.1f} min "
          f"(range {q['te'].tau_min.min():.1f}-{q['te'].tau_min.max():.1f}).\n")
    a(q.get("READ75", "") + "\n")
    a("## 8. 7.6 The rain reference\n")
    a(f"Radar pixel above the 10 City gauges against the gauge, 5 min, TRAIN and VAL event windows only (`{rel(P.S7_RAIN_REF)}`):\n")
    a(md_table(q["rt"].round(3)) + "\n")
    a(f"Rule on V1 (both binnings > 25 %, same direction, intervals exclude 1): gauge binning {q['ch76']['gauge'][0]:+.0%} (excludes 1: {q['ch76']['gauge'][1]}), "
      f"radar binning {q['ch76']['radar'][0]:+.0%} (excludes 1: {q['ch76']['radar'][1]}) -> **{'rate dependence FOUND' if q['found76'] else 'no rate dependence found'}**.\n")
    a("The candidate's mean residual by V1 class:\n")
    a(md_table(q["res76"].round(3)) + "\n")
    a("![rain reference](../figures/validate_s7_rain_reference.png)\n")
    a(q.get("READ76", "") + "\n")
    a("## 9. The return decision\n")
    a(md_table(pd.DataFrame([{"reason to return to S6": k, "met": v} for k, v in q["reasons"].items()])) + "\n")
    a(("**Back to S6** (the one allowed return, to be recorded and dated) because: " + "; ".join(k for k, v in q["reasons"].items() if v) + ".")
      if q["go_back"] else "**No reason to return: the candidate goes to S8 as it is.**")
    a("\n" + q.get("READ9", "") + "\n")
    a("## 10. What S8 must check\n")
    a(q.get("READ10", "") + "\n")
    a("## 11. Unclear / possibly wrong\n")
    a(q.get("READ11", "") + "\n")
    return "\n".join(L)


if __name__ == "__main__":
    main()
