"""wet_antenna step S6b: the ONE allowed return to S6 (S7 reason (d): the drying time falls with the wind), plus the S7 misfit
after the rain. Every effect seen so far is put in a box: A in the equation, B identified but outside it, C unexplained.
The equation changes ONLY for what passes the box-A rule; not adding a term is a valid outcome.

  Part 1  diagnose the after-rain misfit cause by cause, with M0 (the S6 candidate) frozen, on the S7 dry-path windows:
          (a) the radar calls the path dry late (gauge end against radar end); (b) the antenna dries faster than the model in
          general; (c) drying depends on the weather or on the preceding rain (partial Spearman, event bootstrap); (d) F too
          high at heavy rain (residual during >= 10 mm/h on links < 1 km).
  Part 2  refinements fitted on TRAIN (everything else as S6): M0r control, M1 tau_f_dry free, M2 tau_f_dry = tau_0
          exp(-k u) with the hourly wind (M3 with the preceding rain only if Part 1 asks for it); a fake check of every
          extra number.
  Part 3  the box-A rule (identified, physical sign, fixes the after-rain error on VAL, does not hurt, not carried by a few
          events, the simpler one preferred); for a refinement that passes: cold bootstrap and identifiability again.
  Part 4  the three boxes (every effect of S2-S7) and the VAL error budget.
  Part 5  the candidate for S8 with its sha256.

STANDING RULE: WAA is never estimated by subtracting a power-law rain term; the path part is separated only through L at
matched rain or by dry-path moments; ITU-R P.838 only as a size reference; the rain correction comes from the City gauges
only (g_e); NO TEST attenuation or rain is read. Every rule is in config `wet_antenna.return_s6b` (written before the run).
Caches: research/.cache/wet_antenna/return_s6b/. Run from the repository root:
    python -m research.wet_antenna.return_s6b
"""
import hashlib
import json
import os
import pickle
import sys
import time
import zlib
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import length_test as LTL
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.final_model import DATA, FIX, FK, J, LT, M, SHAPES, WA, Curves, plan_start, s6_models
from research.wet_antenna.tails import itu_p838
from research.wet_antenna.validation import LB, PADM, boot_items, ci, dry_windows, error_budget, fmt, load, need, radar_max
from research.wet_antenna.wind_drying import ProblemW, wind_hourly, wind_minutes

S = WA["return_s6b"]
CACHE = P.S6B_CACHE
VL = S["val_label"]
KB = tuple(S["refinements"]["k_bounds"])
K5 = int(np.nonzero(np.asarray(M["knots_mmh"], float)[1:] == 5.0)[0][0]) + 1
STANDING = ("**Standing rule.** WAA is never estimated by subtracting a power-law rain term; the path part is separated only through L at matched "
            "rain or by dry-path moments; ITU-R P.838 only as a size reference; the rain correction comes from the City gauges only (g_e); no TEST "
            "attenuation or rain was read. The S6 candidate stays the reference (M0); every new fit is compared with it.")
REFS = {"M0r": [], "M1": ["tau_f_dry"], "M2": ["tau_f_dry", "k_wind"], "M3": ["tau_f_dry", "k_wind"]}
COV = {"M0r": "wind", "M1": "wind", "M2": "wind", "M3": "peak_rain"}      # the covariate u of k_wind (k = 0 for M0r and M1)


# =================================================================================================== the parametrisation
def spec_shared_lin():
    return {"F": "shared", "P": ("form", {**SHAPES["linear"]}, True)}


def covariate_raw(name, D, le, wind):
    """The covariate per minute in the layout of D["rain"]: the hourly wind (m/s), or the preceding peak rain (max 1 min V1 in
    the previous 60 min)."""
    if name == "wind":
        return wind_minutes(D, le, wind, 0.0, 1.0)
    out = np.zeros(len(D["rain"]))
    for e in range(len(le)):
        o, n = D["le_roff"][e], D["le_n"][e] + 2 * PADM
        out[o:o + n] = pd.Series(D["rain"][o:o + n]).rolling(60, min_periods=1).max().to_numpy()
    return out


_W = {}                                                           # per process: the TRAIN data and covariates


def worker_data(jb):
    if "Dt1" not in _W:
        _W["Dt1"] = sid.load_data(P.S5A_CACHE / "data_train_v1.npz", jb["bands"])
        _W["let"] = pickle.loads((P.SID3_CACHE / "le_train_all.pkl").read_bytes())[0].reset_index(drop=True)
        _W["wind"] = wind_hourly(P.DEW_RISK_HOURLY)
    key = ("cov", jb["cov"])
    if key not in _W:
        raw = covariate_raw(jb["cov"], _W["Dt1"], _W["let"], _W["wind"])
        _W[key] = (raw - jb["cov_stats"][0]) / jb["cov_stats"][1]
    return _W["Dt1"], _W["let"], _W[key]


def truth_minutes(D, le, x, d, uw):
    """W + L P at 1 min (c = 0) of a ProblemW parameter vector: the truth of fake data."""
    D0 = dict(D)
    D0["le_c"] = np.zeros(len(le))
    return ProblemW(D0, M, uw, KB, agg=1).predict(np.asarray(x), d)


def job(jb):
    file = CACHE / f"{jb['name']}.pkl"
    if file.exists():
        return jb["name"], pickle.loads(file.read_bytes())
    t0 = time.time()
    D, le, uw = worker_data(jb)
    if jb["kind"] in ("fake", "ident"):
        tr = truth_minutes(D, le, jb["x_truth"], jb["d"], uw)
        D = sid.fake_from_truth(D, le, tr, jb["seed"], sid.noise_table(), FK, DATA["offset_min_minutes"])
    wle = np.ones(len(le))
    if jb.get("drop_events"):
        wle[le.event_id.isin(jb["drop_events"]).to_numpy()] = 0.0
    if jb["kind"] == "boot":
        n_ev = int(D["le_event"].max()) + 1
        rng = np.random.default_rng(jb["seed"])
        mults = [np.bincount(rng.integers(0, n_ev, n_ev), minlength=n_ev).astype(float) for _ in range(jb["i"] + 1)]
        wle = mults[jb["i"]][D["le_event"]]
    Cv = Curves(D, spec_shared_lin(), np.asarray(jb["x_ref"]), jb["d"], wle=wle, Pm=ProblemW(D, M, uw, KB), free=jb["free"])
    z = np.asarray(jb["z0"])
    nit = 0
    for cap in jb["caps"]:
        f = Cv.fit(z, cap)
        z, nit = f["z"], nit + f["nit"]
    out = {"z": z, "loss": f["loss"], "converged": f["converged"], "nit": nit, "seconds": time.time() - t0}
    file.write_bytes(pickle.dumps(out))
    return jb["name"], out


def run_pool(jobs, workers, clock, results):
    jobs = [j for j in jobs if j["name"] not in results]
    if not jobs:
        return results
    with ProcessPoolExecutor(workers) as ex:
        running = {ex.submit(job, j): j for j in jobs[:workers]}
        queue = jobs[workers:]
        while running:
            done_, _ = wait(running, return_when=FIRST_COMPLETED)
            for fu in done_:
                running.pop(fu)
                name, out = fu.result()
                results[name] = out
                print(f"   {name} done ({(time.time() - clock) / 60:.1f} min)", flush=True)
                if queue:
                    j = queue.pop(0)
                    running[ex.submit(job, j)] = j
    return results


def readout(Cv, z, bands):
    F, Pv, a = Cv.values(z)
    out = {}
    for b, nm in enumerate(bands):
        out[f"A[{nm}]"] = float(F[b, K5])
        out[f"c[{nm}]"] = float(Pv[b, K5] / 5)
    for j, k in enumerate(Cv.kn[1:], 1):
        if k != 5:
            out[f"G({k:g})"] = float(F[0, j] / F[0, K5])
    la = np.log(a)
    out.update({"a_p10": float(np.exp(np.percentile(la, 10))), "a_p90": float(np.exp(np.percentile(la, 90)))})
    for i, p in enumerate(Cv.free):
        v = z[Cv.nc + Cv.nl + i]
        out["tau_0" if p == "tau_f_dry" else p] = float(np.exp(v) if p == "tau_f_dry" else v)
    return out


# =================================================================================================== evaluation of a model
class Model:
    """A fitted model on TRAIN and VAL: 5 min predictions, 1 min W (path term 0, c = 0), the VAL judge and the length test."""

    def __init__(self, name, fz, sets, uw_sets):
        self.name, self.fz, self.p5, self.W1 = name, fz, {}, {}
        for s_, (D, le) in sets.items():
            ids = LTL.link_ids(D, le)
            Pm = ProblemW(D, fz["M"], uw_sets[s_], KB)
            x, _ = LTL.x_for(Pm, fz, ids)
            self.p5[s_] = Pm.predict(x, fz["d"])
            D0 = dict(D)
            D0["le_c"] = np.zeros(len(le))
            P1 = ProblemW(D0, fz["M"], uw_sets[s_], KB, agg=1)
            x1, _ = LTL.x_for(P1, fz, ids)
            x1[P1.names == P1.pname] = 0.0
            self.W1[s_] = P1.predict(x1, fz["d"])


def fz_w(Cv, z, ids):
    """fz of a Curves fit with ProblemW (its x includes k_wind; LTL.x_for keeps it)."""
    return {"x": Cv.build(z), "d": Cv.d, "a": dict(zip(ids, np.exp(z[Cv.nc:Cv.nc + Cv.nl]))), "model": "M4b", "M": Cv.M}


def fz_m0(mods):
    """The frozen S6 candidate in the ProblemW layout (k_wind = 0): the same predictions as S7."""
    f = dict(mods["candidate"]["fz"])
    f["x"] = np.r_[f["x"], 0.0]
    return f


def length_test(D, le, preds):
    """The VAL length test through the model (as final_model.through) for given 5 min predictions."""
    t = LTL.steady_table(D, le, LT, preds)
    n_ev = int(D["le_event"].max()) + 1
    obs = LTL.length_fits(t, n_ev, "VAL", LT)
    rows = []
    for k in preds:
        mod = LTL.length_fits(t, n_ev, "VAL", LT, ycol=f"pred_{k}", boot=False).set_index(["band", "class"])
        n_ok, n_el = 0, 0
        for r in obs.dropna(subset=["I"]).to_dict("records"):
            key = (r["band"], r["class"])
            if key not in mod.index or not np.isfinite(mod.loc[key, "I"]) or r["links_under_1km"] < LT["verdict"]["min_short_links"]:
                continue
            q = mod.loc[key]
            n_el += 1
            n_ok += int(r["I_lo"] <= q.I <= r["I_hi"] and r["s_lo"] <= q.s <= r["s_hi"])
        rows.append({"model": k, "eligible_cells": n_el, "both_inside": n_ok, "passes": bool(n_el and n_ok >= need(n_el))})
    return pd.DataFrame(rows).set_index("model")


# =================================================================================================== Part 1: the diagnosis
def window_table(d, sets, win, curves, wind_h):
    """Per window: errors of M0, decay ratio, weather of the hour of t0, preceding peak rain, gauge end and lag."""
    w = win.copy()
    obs = np.array([c["obs"] for c in curves])
    m0 = np.array([c["M0"] for c in curves])
    for t in (0, 5, 10):
        ok = (w.dry_minutes.to_numpy() > t)
        w[f"err{t}"] = np.where(ok, obs[:, t] - m0[:, t], np.nan)
    # (b) apparent decay over 0-10 min: ln(W(9-11) / W(0-1)), measured and model, clipped to (0.02, 0.98)
    def frac(a):
        with np.errstate(all="ignore"):
            f = np.nanmean(a[:, 9:12], axis=1) / np.nanmean(a[:, 0:2], axis=1)
        return np.clip(f, 0.02, 0.98)
    ok10 = (w.dry_minutes.to_numpy() > 11) & (w.W0_obs.to_numpy() >= 0.5)
    fo, fm = frac(obs), frac(m0)
    w["tau10_obs"] = np.where(ok10, -10 / np.log(fo), np.nan)
    w["tau10_M0"] = np.where(ok10, -10 / np.log(fm), np.nan)
    w["decay_ratio"] = w.tau10_obs / w.tau10_M0
    # (c) weather and the preceding peak rain
    met = pd.read_csv(P.DEW_RISK_HOURLY, parse_dates=["time_utc"]).set_index("time_utc")
    hr = pd.DatetimeIndex(w.t0).floor("h")
    for v in ("wind_speed", "humidity", "temperature"):
        w[v] = met[v].reindex(hr).to_numpy(float)
    peak = []
    for r in w.itertuples():
        D, le = sets[r.set]
        o = D["le_roff"][r.le_idx] + PADM
        m0_ = int((pd.Timestamp(r.t0) - pd.Timestamp(le.window_start.iloc[r.le_idx])) / pd.Timedelta("1min"))
        peak.append(float(np.max(D["rain"][o + max(m0_ - 60, -PADM):o + m0_])) if m0_ > 0 else np.nan)
    w["peak_rain_60"] = peak
    # (a) the gauge end
    Q = S["diagnosis"]
    lgm = d.table("link_gauge_matches")
    lgm = lgm[(lgm.network == "city") & (lgm.distance_to_path_m <= Q["gauge_max_m"])].sort_values("distance_to_path_m").drop_duplicates("cml_id").set_index("cml_id")
    w["gauge_id"] = w.cml_id.map(lgm.gauge_id)
    w["gauge_dist_m"] = w.cml_id.map(lgm.distance_to_path_m)
    ev = d.events()
    test = ev[ev.split_event == "test"]
    gs = d.store("gauges_city_1min").rate_mm_h
    w["gauge_end"] = pd.NaT
    for e_id, g in w[w.gauge_id.notna()].groupby("event_id"):
        t_lo, t_hi = pd.Timestamp(g.t0.min()) - pd.Timedelta("90min"), pd.Timestamp(g.t0.max()) + pd.Timedelta("61min")
        fr = gs.sel(time=slice(t_lo, t_hi)).to_pandas()
        for te in test.itertuples():                              # never a TEST minute
            fr.loc[(fr.index >= te.start) & (fr.index < te.end)] = np.nan
        for i, r in g.iterrows():
            if r.gauge_id not in fr.columns:
                continue
            t0 = pd.Timestamp(r.t0)
            s_ = fr[r.gauge_id][(fr.index >= t0 - pd.Timedelta("90min")) & (fr.index < t0 + pd.Timedelta(minutes=int(r.dry_minutes)))]
            wet = s_[s_ > 0]
            if len(wet):
                w.at[i, "gauge_end"] = wet.index[-1] + pd.Timedelta("1min")
    w["lag_min"] = (pd.to_datetime(w.t0) - pd.to_datetime(w.gauge_end)) / pd.Timedelta("1min")
    return w


def pspear(x, y, z=None):
    """Spearman of x and y, partial on z (ranks, linear residuals) if z is given."""
    ok = np.isfinite(x) & np.isfinite(y) & (np.isfinite(z) if z is not None else True)
    if ok.sum() < 10:
        return np.nan
    rx, ry = stats.rankdata(x[ok]), stats.rankdata(y[ok])
    if z is not None:
        rz = stats.rankdata(z[ok])
        X = np.c_[np.ones(len(rz)), rz]
        rx = rx - X @ np.linalg.lstsq(X, rx, rcond=None)[0]
        ry = ry - X @ np.linalg.lstsq(X, ry, rcond=None)[0]
    return float(np.corrcoef(rx, ry)[0, 1])


def diagnose(wt, sets, m0_p5, curves):
    Q = S["diagnosis"]
    ev = pd.factorize(wt.event_id)[0]
    draws = boot_items(ev, Q["boot_reps"], Q["seed"])
    rows, out = [], {}
    # (a) gauge end
    g = wt[wt.lag_min.notna()]
    ig = g.index.to_numpy()
    agree = g[g.lag_min.abs() <= Q["agree_lag_min"]]
    ia = agree.index.to_numpy()
    rho_lag = pspear(g.lag_min.to_numpy(float), g.err10.to_numpy(float))
    bs = [pspear(wt.lag_min.to_numpy(float)[i], wt.err10.to_numpy(float)[i]) for i in draws[:400]]
    lo_l, hi_l = ci(bs)
    agree_e = {t: float(np.nanmedian(agree[f"err{t}"])) for t in (5, 10)}
    agree_ci = {t: ci([np.nanmedian(wt[f"err{t}"].to_numpy()[i[np.isin(i, ia)]]) if np.isin(i, ia).any() else np.nan for i in draws[:400]]) for t in (5, 10)}
    sup_a = bool(all(abs(v) < 0.1 for v in agree_e.values()) and np.isfinite(hi_l) and hi_l < 0)
    obs = np.array([c["obs"] for c in curves])
    m0 = np.array([c["M0"] for c in curves])
    since_g = []
    for tau in range(0, 61):
        vals = []
        for i in ig:
            t = int(round(tau - wt.lag_min[i]))
            if 0 <= t < wt.dry_minutes[i] and t < 60:
                vals.append(obs[i, t] - m0[i, t])
        since_g.append({"minutes": tau, "since_gauge_end": float(np.nanmedian(vals)) if len(vals) >= 10 else np.nan, "n_gauge": len(vals),
                        "since_radar_end": float(np.nanmedian(obs[ig, tau] - m0[ig, tau])) if tau < 60 else np.nan})
    out["since"] = pd.DataFrame(since_g)
    rows.append({"cause": "(a) the radar calls the path dry late", "windows": len(g), "events": int(g.event_id.nunique()),
                 "numbers": f"lag (radar end - gauge end) median {g.lag_min.median():+.0f} min (IQR {g.lag_min.quantile(0.25):+.0f} to {g.lag_min.quantile(0.75):+.0f}); "
                            f"ends agree (|lag| <= {Q['agree_lag_min']} min) in {len(agree)} windows: median error at 5 / 10 min "
                            f"{agree_e[5]:+.2f} {fmt(np.nan, *agree_ci[5])} / {agree_e[10]:+.2f} {fmt(np.nan, *agree_ci[10])} dB (all gauge windows "
                            f"{np.nanmedian(g.err5):+.2f} / {np.nanmedian(g.err10):+.2f}); Spearman error at 10 min vs lag {rho_lag:+.2f} [{lo_l:+.2f}, {hi_l:+.2f}]",
                 "verdict": "supported" if sup_a else "not supported"})
    # (b) faster drying in general
    rat = wt.decay_ratio.to_numpy(float)
    hv = wt.heavy.to_numpy(bool)
    res_b = {}
    for nm, sel in (("all", np.ones(len(wt), bool)), ("ordinary", ~hv), ("heavy", hv)):
        m = sel & np.isfinite(rat)
        med = float(np.median(rat[m])) if m.any() else np.nan
        lo_, hi_ = ci([np.nanmedian(rat[i[m[i]]]) if m[i].any() else np.nan for i in draws])
        res_b[nm] = (med, lo_, hi_, int(m.sum()))
    sup_b = bool(all(res_b[k][0] <= Q["fast_ratio_max"] and res_b[k][2] < 1 for k in ("all", "ordinary")))
    e0 = float(np.nanmedian(wt.err0))
    e0_ci = ci([np.nanmedian(wt.err0.to_numpy()[i]) for i in draws])
    rows.append({"cause": "(b) the antenna dries faster than the model in general", "windows": res_b["all"][3], "events": int(wt[np.isfinite(rat)].event_id.nunique()),
                 "numbers": "apparent decay time 0-10 min, measured / model (median): " + "; ".join(f"{k} {v[0]:.2f} [{v[1]:.2f}, {v[2]:.2f}] ({v[3]} windows)" for k, v in res_b.items())
                            + f"; the error AT t = 0 (before any drying): {e0:+.2f} [{e0_ci[0]:+.2f}, {e0_ci[1]:+.2f}] dB",
                 "verdict": "supported" if sup_b else "not supported"})
    out["decay"] = res_b
    out["e0"] = (e0, *e0_ci)
    # (c) weather and preceding rain
    tests = []
    sign = {"wind_speed": -1, "peak_rain_60": -1, "humidity": +1, "temperature": -1}
    for target in ("err10", "decay_ratio"):
        y = wt[target].to_numpy(float)
        for v, zc in (("wind_speed", None), ("peak_rain_60", None), ("humidity", None), ("temperature", None),
                      ("wind_speed", "peak_rain_60"), ("peak_rain_60", "wind_speed"), ("humidity", "peak_rain_60"), ("temperature", "peak_rain_60")):
            x = wt[v].to_numpy(float)
            z = wt[zc].to_numpy(float) if zc else None
            r = pspear(x, y, z)
            bs = [pspear(x[i], y[i], z[i] if z is not None else None) for i in draws[:400]]
            lo_, hi_ = ci(bs)
            ok_ = bool(np.isfinite(lo_) and (lo_ > 0 or hi_ < 0) and np.sign(r) == sign[v])
            tests.append({"target": target, "variable": v, "given": zc or "", "rho": r, "lo": lo_, "hi": hi_, "physical_and_excludes_0": ok_})
    tt = pd.DataFrame(tests)
    out["weather"] = tt

    def sig(v, given, target):
        q = tt[(tt.target == target) & (tt.variable == v) & (tt.given == given)]
        return bool(q.physical_and_excludes_0.iloc[0]) if len(q) else False
    verd = {}
    for v, other in (("wind_speed", "peak_rain_60"), ("peak_rain_60", "wind_speed"), ("humidity", "peak_rain_60"), ("temperature", "peak_rain_60")):
        alone = sig(v, "", "err10") or sig(v, "", "decay_ratio")
        part = sig(v, other, "err10") or sig(v, other, "decay_ratio")
        if part:
            verd[v] = "supported"
        elif alone and v in ("wind_speed", "peak_rain_60") and not (sig(other, v, "err10") or sig(other, v, "decay_ratio")):
            verd[v] = "not separable"
        else:
            verd[v] = "not supported"
    out["weather_verdict"] = verd
    rows.append({"cause": "(c) drying depends on the weather / the preceding rain", "windows": int(np.isfinite(wt.err10).sum()), "events": int(wt.event_id.nunique()),
                 "numbers": "; ".join(f"{r.variable}{' | ' + r.given if r.given else ''} vs {r.target}: {r.rho:+.2f} [{r.lo:+.2f}, {r.hi:+.2f}]" for r in tt.itertuples()),
                 "verdict": ", ".join(f"{k}: {v}" for k, v in verd.items())})
    # (d) F too high at heavy rain, links < 1 km
    res_d = {}
    for s_, (D, le) in sets.items():
        w5 = D["w5"] > 0
        Lb = le.L_km.to_numpy(float)[D["bin_le"]]
        r5 = np.nan_to_num(D["r5"], nan=-1)
        g_ = w5 & (Lb < 1) & (r5 >= 10)
        r = (D["y5"] - m0_p5[s_])[g_]
        evb = D["le_event"][D["bin_le"]][g_]
        dd = boot_items(evb, 1000, Q["seed"] + 1)
        lo_, hi_ = ci([np.mean(r[i]) for i in dd if len(i)])
        res_d[s_] = (float(np.mean(r)) if len(r) else np.nan, lo_, hi_, float(np.median(r)) if len(r) else np.nan, int(g_.sum()), int(len(np.unique(evb))))
    sup_d = bool(all(np.isfinite(v[0]) and v[0] <= Q["f_heavy_max_db"] and v[2] < 0 for v in res_d.values()))
    rows.append({"cause": "(d) F too high at heavy rain", "windows": sum(v[4] for v in res_d.values()), "events": sum(v[5] for v in res_d.values()),
                 "numbers": "; ".join(f"{k}: mean residual during V1 >= 10 mm/h on L < 1 km {v[0]:+.2f} [{v[1]:+.2f}, {v[2]:+.2f}] dB (median {v[3]:+.2f}; {v[4]} bins, {v[5]} events)"
                                      for k, v in res_d.items()),
                 "verdict": "supported" if sup_d else "not supported"})
    out["f_heavy"] = res_d
    out["table"] = pd.DataFrame(rows)
    out["supported"] = {"a": sup_a, "b": sup_b, "c": verd, "d": sup_d}
    return out


# =================================================================================================== Part 3 helpers
def dry_errors(win, curves, model, sets_=("VAL",), heavy=None):
    sel = win.set.isin(sets_).to_numpy()
    if heavy is not None:
        sel = sel & (win.heavy.to_numpy(bool) == heavy)
    obs = np.array([c["obs"] for c in curves])
    mod = np.array([c[model] for c in curves])
    return sel, obs, mod


def improvement(win, curves, mx, sets_=("VAL",), heavy=None, boot=True, seed=8401):
    """|median error M0| - |median error Mx| at 5 and 10 min (windows still dry at t), with an event-bootstrap interval."""
    sel, obs, m0 = dry_errors(win, curves, "M0", sets_, heavy)
    mxv = np.array([c[mx] for c in curves])
    nmin = win.dry_minutes.to_numpy()
    ev = pd.factorize(win.event_id)[0]
    draws = boot_items(ev[sel], 1000, seed) if boot else []
    idx = np.nonzero(sel)[0]
    out = {}
    for t in (5, 10):
        ok = idx[(nmin[idx] > t) & np.isfinite(obs[idx, t])]
        e0, ex = obs[ok, t] - m0[ok, t], obs[ok, t] - mxv[ok, t]
        point = abs(np.median(e0)) - abs(np.median(ex)) if len(ok) else np.nan
        bs = []
        for i in draws:
            ii = idx[i]
            ii = ii[(nmin[ii] > t) & np.isfinite(obs[ii, t])]
            if len(ii):
                bs.append(abs(np.median(obs[ii, t] - m0[ii, t])) - abs(np.median(obs[ii, t] - mxv[ii, t])))
        out[t] = (float(point), *ci(bs), float(np.median(e0)) if len(ok) else np.nan, float(np.median(ex)) if len(ok) else np.nan, len(ok))
    return out


def length_diff(win, curves, mx, seed=8402):
    """The 7.3a difference of the median error at 10 min, L > 3 km minus L < 1 km, TRAIN + VAL windows, event bootstrap."""
    obs = np.array([c["obs"] for c in curves])
    m = np.array([c[mx] for c in curves])
    t = 10
    ok = (win.dry_minutes.to_numpy() > t) & np.isfinite(obs[:, t])
    lo_s, hi_s = ok & (win.length.to_numpy() == LB[0][1]), ok & (win.length.to_numpy() == LB[2][1])
    err = obs[:, t] - m[:, t]
    draws = boot_items(pd.factorize(win.event_id)[0], 1000, seed)
    point = float(np.median(err[hi_s]) - np.median(err[lo_s]))
    bs = [np.median(err[i[hi_s[i]]]) - np.median(err[i[lo_s[i]]]) for i in draws if hi_s[i].any() and lo_s[i].any()]
    return (point, *ci(bs))


# =================================================================================================== Part 4: boxes and budget
def boxes(q):
    """Every effect seen in S2-S7 with its box (A / B / C), evidence, size, conditions and the reason."""
    dg = q["dg"]
    wind_box = q["wind_box"]
    sup = dg["supported"]
    e0 = dg["e0"]
    rows = [
        ("rain rate -> antenna level F(R)", "A", "S6: A = 2.74 [2.20, 3.37] / 2.56 [1.88, 3.39] / 2.96 [2.49, 3.53] dB at 5 mm/h; fake 62 / 62; VAL length test 5 / 6",
         "about 1 dB at 0.5-1 mm/h, 2.7 dB at 5 mm/h (equilibrium; 0.75 x after 30 min)", "all rain", "strong, identified, holds on VAL"),
        ("one shape of F and P for all bands", "A", "S6: VAL 0.3875 vs 0.3878 dB; S7.2: short links, shared vs band shapes +0.0005 / -0.0027 dB",
         "-", "28 / 32 / 38 GHz", "the band shapes do not generalise (TRAIN AICc prefers them, VAL does not)"),
        ("wetting / drying memory (fast store)", "A", "S2b: tau_wet 2.6, tau_dry 10.3 min (gauge input); S3.4 VAL 2.0 / 8.8 min"
         + (f"; S6b {q['chosen']}: tau_f_dry {q['rd'][q['chosen']]['ro']['tau_0']:.1f} min with the radar input (an effective time)" if q["chosen"] else ""),
         "the whole dynamics", "all events", "physical, replicated on VAL" + ("; the drying time of the candidate is the S6b refit, not the gauge timing"
                                                                         if q["chosen"] else "; fixed at the gauge timing")),
        ("slow store", "A (structure); numbers C", "S2b: tau_s at its 240 min bound; S5a: M4c tau_s_wet at the 600 min bound, phi 0.19-0.84",
         "1 - phi = 28 % of the level, slow", "long events", "the structure helps (after-rain integral 0.79-1.21 x), its numbers are fixed, not identified"),
        ("path term P = c_band R", "A (S6)" + ("; NOT identified with " + q["chosen"] if q["chosen"] and not q["boot_tab"].set_index("quantity").identified.filter(like="c[").all() else ""),
         "S6: c = 0.107 / 0.142 / 0.197 dB/km per mm/h, plug-back holds; S7.3e 0.5-0.7 x ITU-R P.838 (size only)"
         + (f"; S6b {q['chosen']}: c halves and its bootstrap interval reaches 0 in every band" if q["chosen"] else ""),
         "L c R: 0.5 dB/km at 5 mm/h (28 GHz)", "R >= 1 mm/h", "identified above 1 mm/h in S6, straight line chosen among physics and PySR forms"
         + ("; the faster fast store of the S6b refinement takes over part of it" if q["chosen"] else "")),
        ("band amplitudes A_band, c_band", "A (numbers); no band difference in A detected", "S6: A intervals overlap (2.20-3.37, 1.88-3.39, 2.49-3.53); c rises with f",
         "-", "3 bands", "numbers per band; 3 frequencies are too few for a frequency law"),
        ("link multiplier a_link", "A (multiplier); its cause C", "S6: 10-90 % 0.76-1.23; S7.1 / 7.5: no metadata relation (L rho +0.05; PySR below the null)",
         "+-25 % of W", "per link", "needed, but nothing in the metadata explains it; the F level assumes median a_link = 1"),
        ("W independent of L", "A (confirmed)", "S7.1 rho +0.05; W height 2.8 dB in every class; dry-path error -0.18 dB in every class; intercepts 13 / 15",
         "-", "0.3-8.7 km", "three independent checks"),
        ("drying faster with wind", wind_box, q["wind_evidence"], q["wind_size"], "hourly wind, one station", q["wind_why"]),
        ("after-rain misfit (model high)", q["misfit_box"], f"S7.3a / S6b Part 1: error at t = 0 {e0[0]:+.2f} [{e0[1]:+.2f}, {e0[2]:+.2f}] dB, -0.18 dB at 10 min; "
         f"causes: (a) {'supported' if sup['a'] else 'not supported'}, (b) {'supported' if sup['b'] else 'not supported'}, (c) "
         f"{', '.join(f'{k} {v}' for k, v in sup['c'].items())}, (d) {'supported' if sup['d'] else 'not supported'}",
         "0.2-0.3 dB at the end of rain, 0.5 dB after heavy rain", "the first 20-40 min after the rain", q["misfit_why"]),
        ("dew", "B", "S7.4: 125 link-nights, median 1.7 dB (peak 2.3) for about 4 h after midnight; rho with a_link +0.26 [+0.06, +0.45]; raises c by 0.13 [0.05, 0.35] dB",
         "1.7 dB, hours", "dry nights at high humidity", "real WAA of the same radome, but not driven by rain; outside the rain equation"),
        ("late water on long links", "C", "S7.3a: decay 30 [23, 51] min and 31 % left at 40-60 min on L > 3 km against 21 [16, 26] min, 20 % on L < 1 km",
         "about 0.05 dB at 40 min", "L > 3 km, dry path", "not separable from drizzle below the radar's 0.05 mm/h along a long path"),
        ("F shape above 10 mm/h", "C", "S6: 0.4 % of the minutes; S7.3d: free knots and (R / 5)^0.47 cannot be told apart", "+-1-2 dB at 20-40 mm/h", "R > 10 mm/h",
         "not enough heavy rain"),
        ("P below 1 mm/h", "C", "S6: P(0.5), P(1) intervals at 0 (8 numbers not identified)", "< 0.1 dB/km", "light rain", "the path term is too small to see"),
        ("rain-reference error between events", "C (the main source of the error bars)", f"S4.3 g_e 0.9-1.7 between sets; S6 bootstrap 1-13 x the profile; S6b budget: "
         f"between events {q['between_share']:.0%} of the squared error on VAL", "the event's whole level", "every event", "the radar / gauge ratio of an event is known only roughly"),
        ("rate-dependent radar error", "closed: not found", "S7.6: gauge binning -57 %, radar binning +88 % (opposite: regression to the mean)", "-", "-",
         "the parked g(R) is not justified"),
    ]
    return pd.DataFrame(rows, columns=["effect", "box", "evidence", "size", "conditions", "why"])


# =================================================================================================== figures
def fig_diag(dg, wt):
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    s = dg["since"]
    ax = axes[0]
    ax.plot(s.minutes, s.since_radar_end, "-", color=INK, label="time from the radar end")
    ax.plot(s.minutes, s.since_gauge_end, "--", color=INK, label="time from the gauge end")
    ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xlabel("minutes after the end of rain")
    ax.set_ylabel("median error W_obs - W_M0 (dB)")
    ax.set_title("(a) gauge end against radar end")
    ax.legend(frameon=False, fontsize=8)
    for ax, v, lab in ((axes[1], "wind_speed", "hourly wind (m/s)"), (axes[2], "peak_rain_60", "peak rain in the hour before (mm/h, V1)")):
        x = wt[v].to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(wt.err10.to_numpy(float))
        qs = np.nanpercentile(x[ok], [0, 20, 40, 60, 80, 100])
        for k, (sub, mk, ls) in enumerate(((~wt.heavy.to_numpy(bool), "o", "-"), (wt.heavy.to_numpy(bool), "s", "--"))):
            xm, ym = [], []
            for a, b in zip(qs[:-1], qs[1:]):
                m = ok & sub & (x >= a) & (x <= b)
                if m.sum() >= 10:
                    xm.append(np.median(x[m]))
                    ym.append(np.median(wt.err10.to_numpy(float)[m]))
            ax.plot(xm, ym, mk + ls, color=INK, mfc="white" if k else INK, label="after heavy rain" if k else "ordinary windows")
        ax.axhline(0, color=INK2, lw=0.8)
        ax.set_xlabel(lab)
        ax.set_title(f"(c) error at 10 min against {'the wind' if v == 'wind_speed' else 'the preceding rain'}")
        ax.legend(frameon=False, fontsize=8)
    for ax in axes:
        style_axes(ax)
    fig.suptitle("S6b Part 1: the after-rain error of M0 by cause (dry-path windows, TRAIN + VAL)", fontsize=9)
    fig.tight_layout()
    return fig


def fig_fake(fk):
    fig, ax = plt.subplots(figsize=(6, 4))
    for k, (nm, g) in enumerate(fk.groupby("quantity")):
        ax.scatter(np.full(len(g), k), g.estimate / g.truth, marker="os^dv"[k % 5], facecolors="none", edgecolors=INK, s=30)
    ax.axhspan(1 - S["refinements"]["fake_tol"], 1 + S["refinements"]["fake_tol"], facecolor="none", edgecolor=INK2, hatch="//", lw=0)
    ax.axhline(1, color=INK2, lw=0.8)
    ax.set_xticks(range(fk.quantity.nunique()), [f"{a}" for a in sorted(fk.quantity.unique())], fontsize=8, rotation=20)
    ax.set_ylabel("estimate / truth (3 replicates)")
    ax.set_title("Part 2: fake recovery of the extra numbers (hatched: +-20 %)", fontsize=9)
    style_axes(ax)
    fig.tight_layout()
    return fig


def fig_dry(win, curves, models):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
    styles = {"measured": "-", "M0": "--", "M1": ":", "M2": "-."}
    obs = np.array([c["obs"] for c in curves])
    nmin = win.dry_minutes.to_numpy()
    for ax, (_, nm) in zip(axes, LB):
        sel = win.length.to_numpy() == nm
        for lab in ["measured"] + [m for m in ("M0", "M1", "M2") if m in models]:
            arr = obs if lab == "measured" else np.array([c[lab] for c in curves])
            med = [np.nanmedian(arr[sel & (nmin > t), t]) if (sel & (nmin > t)).sum() >= 3 else np.nan for t in range(60)]
            ax.plot(range(60), med, styles[lab], color=INK, lw=1.4, label=lab if lab == "measured" else f"model {lab}")
        ax.set_title(f"{nm}: {int(sel.sum())} windows")
        ax.set_xlabel("minutes after the path became dry")
        style_axes(ax)
    axes[0].set_ylabel("median W (dB)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("S6b: dry-path W, measured against M0 (S6 candidate), M1 (tau_f_dry free), M2 (tau_f_dry with the wind); TRAIN + VAL", fontsize=9)
    fig.tight_layout()
    return fig


def fig_boxes(bx):
    fig, ax = plt.subplots(figsize=(13, 0.42 * len(bx) + 1))
    ax.axis("off")
    tab = ax.table(cellText=[[r.box.split(";")[0].split(" (")[0] + (" *" if "NOT identified" in r.box else ""), r.effect, r.size[:45], r.why[:80]]
                             for r in bx.itertuples()],
                   colLabels=["box", "effect", "size", "why in this box"],
                   cellLoc="left", colLoc="left", loc="center", colWidths=[0.14, 0.22, 0.24, 0.40])
    tab.auto_set_font_size(False)
    tab.set_fontsize(7.5)
    tab.scale(1, 1.35)
    for (r, c), cell in tab.get_celld().items():
        cell.set_edgecolor(INK2)
        if r == 0:
            cell.set_text_props(weight="bold")
    ax.set_title("What the model claims: A in the equation, B identified outside it, C unexplained", fontsize=9)
    if bx.box.str.contains("NOT identified").any():
        ax.text(0.0, 0.02, "* in box A in S6, but NOT identified with the S6b refinement (its bootstrap interval reaches 0)", transform=ax.transAxes, fontsize=8)
    fig.tight_layout()
    return fig


def fig_budget(bud):
    b = bud[~bud.part.isin(["all"])]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    hatches = ["", "//", "..", "xx", "\\\\"]
    for i, r in enumerate(b.itertuples()):
        ax.barh(i, r.share_of_squared_error, color="white", edgecolor=INK, hatch=hatches[i % len(hatches)])
        ax.text(r.share_of_squared_error + 0.01, i, f"{r.share_of_squared_error:.0%} (RMSE {r.rmse_db:.2f} dB)", va="center", fontsize=8)
    ax.set_yticks(range(len(b)), b.part, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.15)
    ax.set_xlabel(f"share of the squared error on VAL ({VL})")
    ax.set_title("The error budget of the candidate", fontsize=9)
    style_axes(ax)
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    clock = time.time()
    workers = max(min((os.cpu_count() or 2) - 1, 3), 1)
    sha6 = hashlib.sha256(P.S6_CANDIDATE.read_bytes()).hexdigest()
    d, BANDS, sets, gmap, mods, free = load()
    Dt1, let = sets["TRAIN"]
    _, (Ck, zk) = s6_models(Dt1, let, BANDS)
    ids = LTL.link_ids(Dt1, let)
    wind = wind_hourly(P.DEW_RISK_HOURLY)
    raw_t = covariate_raw("wind", Dt1, let, wind)
    inside = np.concatenate([np.arange(o + PADM, o + PADM + n) for o, n in zip(Dt1["le_roff"], Dt1["le_n"])])
    wstats = (float(raw_t[inside].mean()), float(raw_t[inside].std()))
    uw = {s_: (covariate_raw("wind", D, le, wind) - wstats[0]) / wstats[1] for s_, (D, le) in sets.items()}

    # ------------------------------------------------------------------------------------------- M0 frozen, the windows
    f0 = fz_m0(mods)
    m0 = Model("M0", f0, sets, uw)
    Rm, T5, col = radar_max(d, sets)
    W1 = {s_: {"M0": m0.W1[s_]} for s_ in sets}
    win, curves = dry_windows(sets, gmap, Rm, T5, col, W1)
    print(f"M0 and {len(win)} windows ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # ------------------------------------------------------------------------------------------- Part 1
    wt = window_table(d, sets, win, curves, wind)
    dg = diagnose(wt, sets, m0.p5, curves)
    print(dg["table"][["cause", "verdict"]].to_string(), flush=True)
    run_m3 = dg["supported"]["c"].get("peak_rain_60") == "supported" and dg["supported"]["c"].get("wind_speed") != "supported"

    # ------------------------------------------------------------------------------------------- Part 2: the fits
    spec = spec_shared_lin()
    x_ref = np.r_[Ck.build(zk), 0.0]
    Cm0 = mods["candidate"]["C"]
    zc = mods["candidate"]["z"]
    _, Pc, _ = Cm0.values(zc)
    c5 = Pc[:, -1] / Cm0.kn[-1] * 5
    nF = Ck.blocks[0][2]
    la0 = zc[Cm0.nc:Cm0.nc + Cm0.nl]
    z_start = {m: np.r_[zk[:nF], c5, la0, [np.log(FIX["tau_f_dry"]), 0.0][:len(REFS[m])]] for m in REFS}
    base = {"bands": BANDS, "x_ref": x_ref.tolist(), "d": Ck.d}

    def covj(m):
        if COV[m] == "wind":
            return {"cov": "wind", "cov_stats": list(wstats)}
        raw = covariate_raw("peak_rain", Dt1, let, wind)
        return {"cov": "peak_rain", "cov_stats": [float(raw[inside].mean()), float(raw[inside].std())]}
    plan_file = CACHE / "plan.json"
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
    else:
        Cx = Curves(Dt1, spec, x_ref, Ck.d, Pm=ProblemW(Dt1, M, uw["TRAIN"], KB), free=REFS["M2"])
        Cx.fun(z_start["M2"])
        t0 = time.time()
        for _ in range(3):
            Cx.fun(z_start["M2"])
        te = (time.time() - t0) / 3
        n_fold = len(win[win.set == "TRAIN"].event_id.unique())
        ev_ = {"fits": 3 * 1500, "fakes": 6 * 1500, "folds": 2 * n_fold * 360, "drop3": 2 * 1200, "boot": 40 * 1500, "ident": 3 * 1500, "m3": 1500 + 3 * 1500}

        def proj(p):
            return (ev_["fits"] + ev_["fakes"] + ev_["drop3"] + ev_["folds"] * p["folds"] / n_fold + p["boot"] / 40 * ev_["boot"] + ev_["ident"]
                    + (ev_["m3"] if p["m3"] else 0)) * te * 1.2 / workers / 3600 + 0.3
        plan = {"t_eval_s": round(te, 3), "workers": workers, "boot": S["box_a_rule"]["bootstrap"]["reps"], "m3": True, "folds": n_fold, "cuts": []}
        plan["projected_full_h"] = round(proj(plan), 2)
        for cut, ch in (("bootstrap 40 -> 30", lambda p: p.update(boot=S["box_a_rule"]["bootstrap"]["reps_if_cut"])), ("M3", lambda p: p.update(m3=False)),
                        ("leave-one-event-out folds -> 10", lambda p: p.update(folds=min(10, p["folds"])))):
            if proj(plan) <= S["time_budget_h"]:
                break
            ch(plan)
            plan["cuts"].append(cut)
        plan["projected_h"] = round(proj(plan), 2)
        plan_file.write_text(json.dumps(plan))
    print(f"timing: one evaluation {plan['t_eval_s'] * 1e3:.0f} ms; projection {plan['projected_full_h']} h (worst case); cuts: {plan['cuts'] or 'nothing'}",
          flush=True)
    refs = ["M0r", "M1", "M2"] + (["M3"] if run_m3 and plan["m3"] else [])
    caps = [S["refinements"]["caps"]["start"], S["refinements"]["caps"]["polish"]]
    results = {}
    run_pool([{**base, **covj(m), "name": f"fit_{m}", "kind": "fit", "free": REFS[m], "z0": z_start[m].tolist(), "caps": caps} for m in refs], workers, clock, results)
    fits, models = {}, {"M0": m0}
    for m in refs:
        Cv = Curves(Dt1, spec, x_ref, Ck.d, Pm=ProblemW(Dt1, M, uw["TRAIN"], KB), free=REFS[m])
        fits[m] = (Cv, np.asarray(results[f"fit_{m}"]["z"]))
    # the VAL covariate for M3 (preceding rain) differs from the wind
    if "M3" in refs:
        st3 = covj("M3")["cov_stats"]
        uw3 = {s_: (covariate_raw("peak_rain", D, le, wind) - st3[0]) / st3[1] for s_, (D, le) in sets.items()}
    for m in refs:
        Cv, z = fits[m]
        models[m] = Model(m, fz_w(Cv, z, ids), sets, uw3 if m == "M3" else uw)
    W1 = {s_: {m: models[m].W1[s_] for m in models} for s_ in sets}
    win, curves = dry_windows(sets, gmap, Rm, T5, col, W1)
    print(f"Part 2 fits done ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # fake check of the extra numbers
    FR = S["refinements"]
    fjobs = []
    for m in [r for r in refs if REFS[r]]:
        Cv, z = fits[m]
        for r in range(FR["fake_replicates"]):
            fjobs.append({**base, **covj(m), "name": f"fake_{m}_{r}", "kind": "fake", "free": REFS[m], "z0": z_start[m].tolist(), "caps": caps,
                          "x_truth": Cv.build(z).tolist(), "seed": [FR["fake_seed"], r, zlib.crc32(m.encode())]})
    run_pool(fjobs, workers, clock, results)
    fk = []
    for m in [r for r in refs if REFS[r]]:
        Cv, z = fits[m]
        tru = readout(Cv, z, BANDS)
        for r in range(FR["fake_replicates"]):
            est = readout(Cv, np.asarray(results[f"fake_{m}_{r}"]["z"]), BANDS)
            for qn in [("tau_0" if p == "tau_f_dry" else p) for p in REFS[m]]:
                t_, e_ = tru[qn], est[qn]
                ok_ = abs(e_ - t_) <= 0.01 if abs(t_) < 0.05 else abs(e_ / t_ - 1) <= FR["fake_tol"]
                i_b = [p for p in REFS[m]].index("tau_f_dry" if qn == "tau_0" else qn)
                zb = np.asarray(results[f"fake_{m}_{r}"]["z"])[Cv.nc + Cv.nl + i_b]
                lo_b, hi_b = Cv.lo[Cv.nc + Cv.nl + i_b], Cv.hi[Cv.nc + Cv.nl + i_b]
                at_b = bool(min(abs(zb - lo_b), abs(zb - hi_b)) <= 0.01 * (hi_b - lo_b))
                fk.append({"model": m, "quantity": f"{m} {qn}", "replicate": r, "truth": t_, "estimate": e_, "within": bool(ok_), "at_bound": at_b,
                           "converged": bool(results[f"fake_{m}_{r}"]["converged"])})
    fkt = pd.DataFrame(fk)
    ident = {m: bool(len(fkt[fkt.model == m]) and fkt[fkt.model == m].within.all() and not fkt[fkt.model == m].at_bound.any()) for m in refs if REFS[m]}
    print(f"fake: {ident} ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # ------------------------------------------------------------------------------------------- Part 3: the box-A rule
    BA = S["box_a_rule"]
    Dv, lev = sets["VAL"]
    rm0 = LTL.val_metrics(Dv, m0.p5["VAL"], J)["rmse"]
    lt_ = length_test(Dv, lev, {m: models[m].p5["VAL"] for m in models})
    fit_rows, rule_rows = [], []
    rd = {}
    for m in refs:
        Cv, z = fits[m]
        ro = readout(Cv, z, BANDS)
        rmse = LTL.val_metrics(Dv, models[m].p5["VAL"], J)["rmse"]
        imp = improvement(win, curves, m)
        imp_h = improvement(win, curves, m, heavy=True, boot=False)
        imp_ht = improvement(win, curves, m, sets_=("TRAIN",), heavy=True, boot=False)
        ld = length_diff(win, curves, m)
        fit_rows.append({"model": m, "free": ", ".join(REFS[m]) or "none", "train_loss": results[f"fit_{m}"]["loss"], "converged": results[f"fit_{m}"]["converged"],
                         "tau_0_min": ro.get("tau_0", FIX["tau_f_dry"]), "k": ro.get("k_wind", 0.0), **{k: v for k, v in ro.items() if k.startswith(("A[", "c["))},
                         "a_p10": ro["a_p10"], "a_p90": ro["a_p90"], f"VAL rmse ({VL})": rmse, "rise_vs_M0": rmse - rm0,
                         "length_test": f"{int(lt_.both_inside[m])} of {int(lt_.eligible_cells[m])}",
                         **{f"VAL err {t} min M0 / Mx": f"{imp[t][3]:+.3f} / {imp[t][4]:+.3f}" for t in (5, 10)}})
        if not REFS[m]:
            continue
        r1 = ident.get(m, False)
        r2 = True if m == "M1" else bool(ro.get("k_wind", 0) > 0)
        r3 = all(imp[t][0] >= BA["rule_3_min_db"] and imp[t][1] > 0 for t in (5, 10)) and all(imp_h[t][0] >= BA["rule_3_heavy_min_db"] for t in (5, 10))
        r4 = bool(rmse - rm0 <= BA["rule_4_rmse_rise_max_db"] and lt_.passes[m] and ld[1] <= 0 <= ld[2])
        rd[m] = {"r1": r1, "r2": r2, "r3": r3, "r4": r4, "imp": imp, "imp_h": imp_h, "imp_ht": imp_ht, "ld": ld, "ro": ro, "rmse": rmse}
    # rule 5 only for a refinement that passes rules 1-4 (else it cannot enter box A)
    train_ev = list(pd.unique(win[win.set == "TRAIN"].event_id))
    wwin = wt.assign(ev=wt.event_id)[wt.set == "TRAIN"].groupby("event_id").wind_speed.mean().sort_values(ascending=False)
    windiest = list(wwin.index[:3])
    fold_rows = []
    for m in rd:
        if not all(rd[m][k] for k in ("r1", "r2", "r3", "r4")):
            rd[m]["r5"] = None
            continue
        Cv, z = fits[m]
        fj = [{**base, **covj(m), "name": f"fold_{m}_{e_}", "kind": "fold", "free": REFS[m], "z0": z.tolist(), "caps": [BA["fold_cap"]], "drop_events": [e_]}
              for e_ in train_ev[:plan["folds"]]]
        fj.append({**base, **covj(m), "name": f"drop3_{m}", "kind": "fold", "free": REFS[m], "z0": z.tolist(), "caps": [1000], "drop_events": windiest})
        run_pool(fj, workers, clock, results)
        signs = []
        for e_ in train_ev[:plan["folds"]]:
            ro_ = readout(Cv, np.asarray(results[f"fold_{m}_{e_}"]["z"]), BANDS)
            v = ro_.get("k_wind") if m != "M1" else ro_["tau_0"] - FIX["tau_f_dry"]
            signs.append(v)
            fold_rows.append({"model": m, "left_out": e_, "k_or_tau0_minus_10.3": v})
        ref_v = rd[m]["ro"].get("k_wind") if m != "M1" else rd[m]["ro"]["tau_0"] - FIX["tau_f_dry"]
        same = bool(all(np.sign(v) == np.sign(ref_v) for v in signs))
        md = Model(f"{m} without the 3 windiest", fz_w(Cv, np.asarray(results[f"drop3_{m}"]["z"]), ids), sets, uw3 if m == "M3" else uw)
        W1x = {s_: {"M0": m0.W1[s_], "X": md.W1[s_]} for s_ in sets}
        _, cx = dry_windows(sets, gmap, Rm, T5, col, W1x)
        imp3 = improvement(win, [dict(c, X=cc["X"]) for c, cc in zip(curves, cx)], "X", boot=False)
        keep = all(imp3[t][0] >= 0.5 * rd[m]["imp"][t][0] for t in (5, 10))
        rd[m]["r5"] = bool(same and keep)
        rd[m]["r5_detail"] = f"sign kept in {sum(np.sign(v) == np.sign(ref_v) for v in signs)} of {len(signs)} folds; without {', '.join(windiest)}: improvement " \
                             f"{imp3[5][0]:+.3f} / {imp3[10][0]:+.3f} dB (full {rd[m]['imp'][5][0]:+.3f} / {rd[m]['imp'][10][0]:+.3f})"
        fold_rows.append({"model": m, "left_out": "3 windiest: " + ", ".join(windiest), "k_or_tau0_minus_10.3": readout(Cv, np.asarray(results[f"drop3_{m}"]["z"]), BANDS).get(
            "k_wind", np.nan) if m != "M1" else readout(Cv, np.asarray(results[f"drop3_{m}"]["z"]), BANDS)["tau_0"] - FIX["tau_f_dry"]})
    passed = [m for m in rd if all(rd[m].get(k) for k in ("r1", "r2", "r3", "r4", "r5"))]
    chosen = None
    if "M1" in passed:                                            # rule 6: M2 must add >= 0.02 dB beyond M1 at 5 and 10 min, else M1
        adds = [abs(rd["M1"]["imp"][t][4]) - abs(rd["M2"]["imp"][t][4]) for t in (5, 10)] if "M2" in passed else [0.0]
        chosen = "M2" if ("M2" in passed and min(adds) >= BA["rule_6_db"]) else "M1"
    elif passed:
        chosen = max(passed, key=lambda m: min(rd[m]["imp"][t][0] for t in (5, 10)))
    for m in rd:
        r = rd[m]
        box = "A" if m == chosen else ("B" if (r["r2"] and r["imp"][5][0] > 0 and r["imp"][10][0] > 0) else "C")
        rule_rows.append({"model": m, "1 identified": r["r1"], "2 physical sign": r["r2"],
                          "3 fixes the after-rain error": r["r3"], "3 VAL improvement 5 / 10 min": f"{fmt(*r['imp'][5][:3], nd=3)} / {fmt(*r['imp'][10][:3], nd=3)}",
                          "3 heavy VAL / TRAIN 5 / 10 min": f"{r['imp_h'][5][0]:+.3f} / {r['imp_h'][10][0]:+.3f}; {r['imp_ht'][5][0]:+.3f} / {r['imp_ht'][10][0]:+.3f}",
                          "4 does not hurt": r["r4"], "4 detail": f"rise {r['rmse'] - rm0:+.4f} dB; length {int(lt_.both_inside[m])} of {int(lt_.eligible_cells[m])}; "
                                                                  f"L>3 - L<1 at 10 min {fmt(*r['ld'], nd=3)}",
                          "5 not a few events": "not run (rules 1-4 not all met)" if r["r5"] is None else r["r5"], "5 detail": r.get("r5_detail", ""), "box": box})
    rule_tab = pd.DataFrame(rule_rows)
    print(f"Part 3: passed {passed}, chosen {chosen} ({(time.time() - clock) / 60:.1f} min)\n{rule_tab.to_string()}", flush=True)

    # ------------------------------------------------------------------------------------------- the new candidate: bootstrap and identifiability
    boot_tab, ident_tab = pd.DataFrame(), pd.DataFrame()
    final = models[chosen] if chosen else m0
    if chosen:
        Cv, z = fits[chosen]
        Pg = sid.Problem2(Dt1, "M4b", M)
        sid.fix_params(Pg, FIX)
        prod_freq = let.groupby("band_name").frequency_ghz.median()
        kn = np.asarray(M["knots_mmh"], float)
        itu_inc = np.array([np.diff(itu_p838(float(prod_freq[b]), "v")[0] * kn ** itu_p838(float(prod_freq[b]), "v")[1]) for b in BANDS])   # START only
        xg = plan_start(Pg, "g4+Phalf", itu_inc)
        zg_shared = Ck.z_from_x(np.where(Ck.Pm.names == "loga", 0.0, xg))
        Pg_inc = xg[Pg.names == Pg.pname].reshape(len(BANDS), -1)
        zg = np.r_[zg_shared[:nF], Pg_inc[:, :K5].sum(axis=1), np.zeros(Cv.nl), [np.log(FIX["tau_f_dry"]), 0.0][:len(REFS[chosen])]]
        bj = [{**base, **covj(chosen), "name": f"boot_{chosen}_{i:02d}", "kind": "boot", "free": REFS[chosen], "z0": zg.tolist(), "caps": caps, "i": i,
               "seed": [BA["bootstrap"]["seed"], zlib.crc32(chosen.encode())]} for i in range(plan["boot"])]
        ij = [{**base, **covj(chosen), "name": f"ident_{chosen}_{r}", "kind": "ident", "free": REFS[chosen], "z0": zg.tolist(), "caps": caps,
               "x_truth": Cv.build(z).tolist(), "seed": [BA["bootstrap"]["seed"] + 1, r]} for r in range(3)]
        run_pool(bj + ij, workers, clock, results)
        est = readout(Cv, z, BANDS)
        draws = pd.DataFrame([readout(Cv, np.asarray(results[j["name"]]["z"]), BANDS) for j in bj if results[j["name"]]["converged"]])
        boot_tab = pd.DataFrame([{"quantity": k, "value": v, "lo": ci(draws[k])[0], "hi": ci(draws[k])[1], "draws": len(draws)} for k, v in est.items()])
        ir = []
        for r in range(3):
            e_ = readout(Cv, np.asarray(results[f"ident_{chosen}_{r}"]["z"]), BANDS)
            for k, t_ in est.items():
                ir.append({"quantity": k, "replicate": r, "truth": t_, "estimate": e_[k],
                           "within": bool(abs(e_[k] - t_) <= 0.01 if abs(t_) < 0.05 else abs(e_[k] / t_ - 1) <= 0.2)})
        ident_tab = pd.DataFrame(ir)
        # the S6 identifiability rule: fake within +-20 % in 3 of 3 AND the bootstrap interval not at a bound
        fk3 = ident_tab.groupby("quantity").within.all()
        at_b = {r.quantity: bool((r.lo <= np.exp(-3) + 1e-9 or r.hi >= np.exp(3) - 1e-9) if r.quantity.startswith("a_") else r.lo <= 0.001)
                for r in boot_tab.itertuples()}
        boot_tab["fake_within_3_of_3"] = boot_tab.quantity.map(fk3)
        boot_tab["interval_at_bound"] = boot_tab.quantity.map(at_b)
        boot_tab["identified"] = boot_tab.fake_within_3_of_3 & ~boot_tab.interval_at_bound

    # ------------------------------------------------------------------------------------------- Part 4: boxes and budget
    bud = error_budget(Dv, lev, final.p5["VAL"])
    bud0 = error_budget(Dv, lev, m0.p5["VAL"])
    bud = bud.merge(bud0[["part", "rmse_db", "share_of_squared_error"]].rename(columns={"rmse_db": "rmse_db M0", "share_of_squared_error": "share M0"}), on="part")
    q = dict(locals())
    q["between_share"] = float(bud[bud.part.str.contains("between")].share_of_squared_error.iloc[0])
    m2 = rd.get("M2")
    if m2:
        kv = m2["ro"].get("k_wind", 0.0)
        per_ms = (1 - np.exp(-kv / wstats[1])) * 100
        q["wind_evidence"] = (f"S7.5 rho -0.76 (17 events, beats the null); S6b M2: k = {kv:+.3f} per SD of the wind (SD {wstats[1]:.2f} m/s), tau_0 "
                              f"{m2['ro'].get('tau_0', np.nan):.1f} min; VAL improvement at 5 / 10 min {m2['imp'][5][0]:+.3f} / {m2['imp'][10][0]:+.3f} dB; "
                              f"identified {m2['r1']}")
        q["wind_size"] = f"tau_f_dry {'shorter' if kv > 0 else 'longer'} by about {abs(per_ms):.0f} % per m/s"
        q["wind_box"] = "A" if chosen == "M2" else ("B" if (m2["r2"] and m2["imp"][5][0] > 0 and m2["imp"][10][0] > 0) else "C")
        q["wind_why"] = {"A": "passes every rule of box A", "B": "right sign and improves, but fails a box-A rule (see Part 3)",
                         "C": "no improvement with the right sign out of sample"}[q["wind_box"]]
    else:
        q.update(wind_evidence="", wind_size="", wind_box="C", wind_why="not fitted")
    sup = dg["supported"]
    q["misfit_box"] = "B" if (sup["a"] or sup["d"] or any(v == "supported" for v in sup["c"].values())) else "C"
    q["misfit_why"] = ("a cause is identified (Part 1) but no refinement fixes it cleanly" if q["misfit_box"] == "B" and not chosen else
                       (f"the level is removed by {chosen} (VAL error at 5 / 10 min {rd[chosen]['imp'][5][4]:+.2f} / {rd[chosen]['imp'][10][4]:+.2f} dB); "
                        "its link to the preceding rain and dry air stays outside" if chosen else "no cause identified"))
    bx = boxes(q)
    q["bx"] = bx

    # ------------------------------------------------------------------------------------------- Part 5: the candidate
    s6 = json.loads(P.S6_CANDIDATE.read_text())
    if chosen:
        cand = {"step": "S6b", "based_on": f"the S6 candidate (sha256 {sha6}) with the refinement {chosen}", "refinement": chosen,
                "refinement_text": S["refinements"][chosen], "numbers": {r.quantity: {"value": r.value, "boot_lo": r.lo, "boot_hi": r.hi} for r in boot_tab.itertuples()},
                "wind_standardisation_m_s": list(wstats), "dynamics_fixed": FIX, "lag_min": int(Ck.d)}
    else:
        cand = {**s6, "step": "S6b", "unchanged_from_S6": f"the S6 candidate unchanged (sha256 {sha6}): no refinement passed box A"}
    cand["boxes"] = {b: bx[bx.box.str.startswith(b)].effect.tolist() for b in ("A", "B", "C")}
    if chosen:
        cand["not_identified"] = boot_tab[~boot_tab.identified].quantity.tolist() + ["P below 1 mm/h (S6)", "phi", "tau_s"]
    cand["closed_not_found"] = bx[bx.box.str.startswith("closed")].effect.tolist()
    txt = json.dumps(cand, indent=1, default=float)
    save(txt, P.S6B_CANDIDATE)
    sha = hashlib.sha256(P.S6B_CANDIDATE.read_bytes()).hexdigest()
    q.update({"cand": cand, "sha": sha, "sha6": sha6, "chosen": chosen})

    # save, figures, report
    for tab, path in ((dg["table"], P.S6B_DIAG), (wt.round(4), P.S6B_DIAG_WINDOWS), (pd.DataFrame(fit_rows).round(4), P.S6B_FITS), (fkt.round(4), P.S6B_FAKE),
                      (rule_tab, P.S6B_RULE), (pd.DataFrame(fold_rows).round(4), P.S6B_FOLDS), (pd.concat([boot_tab.assign(kind="bootstrap"), ident_tab.assign(kind="identifiability")],
                                                                                                       ignore_index=True).round(4), P.S6B_BOOT),
                      (bx, P.S6B_BOXES), (bud.round(4), P.S6B_BUDGET)):
        save(tab, path, index=False)
    figs = [(fig_diag(dg, wt), P.S6B_FIG_DIAG), (fig_fake(fkt), P.S6B_FIG_FAKE), (fig_dry(win, curves, models), P.S6B_FIG_DRY), (fig_boxes(bx), P.S6B_FIG_BOXES),
            (fig_budget(bud), P.S6B_FIG_BUDGET)]
    for fig, path in figs:
        save(fig, path, dpi=100)
        plt.close(fig)
    q.update(dict(locals()))
    q.update(readings(q))
    report = write_report(q)
    save(report, P.S6B_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {(time.time() - clock) / 60:.1f} min")


# =================================================================================================== readings
def readings(q):
    """The READING paragraphs, the candidate box and sections 7-8, written from this run's numbers."""
    R = {}
    dg, wt, rd = q["dg"], q["wt"], q["rd"]
    g = wt[wt.lag_min.notna()]
    tt = dg["weather"].set_index(["target", "variable", "given"])

    def rr(t, v, gv=""):
        r = tt.loc[(t, v, gv)]
        return f"{r.rho:+.2f} [{r.lo:+.2f}, {r.hi:+.2f}]"
    dec, e0, fh = dg["decay"], dg["e0"], dg["f_heavy"]
    R["READ2"] = (
        "**READING.** (a) **Not the radar's timing.** The radar calls the path dry a median "
        f"{g.lag_min.median():.0f} min after the gauge (IQR {g.lag_min.quantile(0.25):.0f}-{g.lag_min.quantile(0.75):.0f} min), but the error does "
        f"not grow with that lag (Spearman at 10 min {pspear(g.lag_min.to_numpy(float), g.err10.to_numpy(float)):+.2f}) "
        "and it is just as large where the two ends agree. (b) **Not a faster drying in general.** Over the first 10 min the measured W decays at "
        f"about the model's rate (ratio {dec['all'][0]:.2f} [{dec['all'][1]:.2f}, {dec['all'][2]:.2f}]); the misfit is ALREADY there when the path "
        f"becomes dry: {e0[0]:+.2f} [{e0[1]:+.2f}, {e0[2]:+.2f}] dB at t = 0. It is a LEVEL at the end of the rain, not a slow drying. "
        f"(c) The end-of-rain error is larger after heavier preceding rain ({rr('err10', 'peak_rain_60', 'wind_speed')} with the wind held fixed) and in "
        f"drier air (humidity {rr('err10', 'humidity', 'peak_rain_60')} with the rain held fixed). The wind alone goes with it too "
        f"({rr('err10', 'wind_speed')}) but NOT once the preceding rain is held fixed ({rr('err10', 'wind_speed', 'peak_rain_60')}): windy hours are "
        "the hours after heavy showers, and the wind adds nothing of its own at the window level. Warm hours dry faster (decay ratio against "
        f"temperature {rr('decay_ratio', 'temperature', 'peak_rain_60')}). All these are small (|rho| <= 0.26). (d) **F is not too high at heavy "
        f"rain**: on links < 1 km during V1 >= 10 mm/h the model is not high (TRAIN {fh['TRAIN'][0]:+.2f} [{fh['TRAIN'][1]:+.2f}, {fh['TRAIN'][2]:+.2f}] dB); "
        f"on VAL it is LOW ({fh['VAL'][0]:+.2f} [{fh['VAL'][1]:+.2f}, {fh['VAL'][2]:+.2f}] dB, {fh['VAL'][5]} events). So the misfit is an end-of-rain "
        "level that grows with the preceding rain and in dry air; the radar timing, a general faster drying and F at heavy rain are ruled out. "
        "By the Part 2 rule M3 (drying with the preceding rain) was run.")
    m1, m2, m3 = rd.get("M1"), rd.get("M2"), rd.get("M3")
    bt = q["boot_tab"].set_index("quantity") if len(q["boot_tab"]) else None
    fr = pd.DataFrame(q["fit_rows"]).set_index("model")
    if m1 and bt is not None:
        notid = bt[~bt.identified].index.tolist()
        R["READ3"] = (
            f"**READING.** **M1 passes every box-A rule.** tau_f_dry = {m1['ro']['tau_0']:.1f} min [{bt.lo['tau_0']:.1f}, {bt.hi['tau_0']:.1f}] instead of "
            f"10.3: the VAL dry-path error goes from {fr.loc['M1', 'VAL err 5 min M0 / Mx']} dB at 5 min and {fr.loc['M1', 'VAL err 10 min M0 / Mx']} dB at "
            f"10 min (improvement {m1['imp'][5][0]:.2f} / {m1['imp'][10][0]:.2f} dB, after heavy rain {m1['imp_h'][5][0]:.2f} / {m1['imp_h'][10][0]:.2f}); "
            f"the VAL RMSE does not rise ({m1['rmse'] - q['rm0']:+.4f} dB); the length test drops from 5 to 4 of 6 (still passes); 23 of 23 folds and the fit "
            "without the 3 windiest events keep it. "
            f"**M2 (wind)**: k = {m2['ro']['k_wind']:+.2f} per SD of the wind (tau_f_dry about 42 % shorter per m/s), identified, the right sign, and it "
            f"improves as much as M1 at 5-10 min, but it fails rule 3 after heavy rain ({m2['imp_h'][5][0]:+.3f} dB at 5 min, rule 0.1) and rule 4 (the "
            f"length-class difference {fmt(*m2['ld'], nd=3)} excludes 0; VAL RMSE {m2['rmse'] - q['rm0']:+.4f} dB): **box B**. Beyond M1 it adds nothing "
            "(-0.062 / -0.058 against -0.065 / -0.051 dB). **M3 (preceding rain)**: k = "
            f"{m3['ro']['k_wind']:+.3f}, essentially 0 with the wrong sign: once tau_f_dry is free, the preceding rain does not change it: **box C**.\n\n"
            "**What M1 does to the rest of the equation (the repeated S6 checks).** With tau_f_dry = 3.7 min, close to tau_f_wet = 2.6 min, the fast store "
            "follows the rain within minutes, so on every link it looks like a path term; only L tells them apart. The refit moves attenuation from the "
            f"path to the antenna: c_band halves (28 GHz {fr.loc['M0r', 'c[28 GHz]']:.3f} -> {fr.loc['M1', 'c[28 GHz]']:.3f} dB/km per mm/h; the same in "
            f"every band) and A_band rises by about a third ({fr.loc['M0r', 'A[28 GHz]']:.2f} -> {fr.loc['M1', 'A[28 GHz]']:.2f} dB at 28 GHz). In the cold "
            "bootstrap 8-11 of the 40 draws put c of a band at exactly 0 (all three bands in 6 draws): by the S6 identifiability rule "
            f"**NOT identified in M1: {', '.join(notid) if notid else 'none'}**. The path term was identified in S6. And 3.7 min is not the antenna's "
            "drying time measured with a gauge (S2b 10.3 min, VAL 8.8 min): it is an effective time for the 5 min radar input, which lets the model drop "
            "W as the radar rain fades.")
    bud = q["bud"].set_index("part")
    R["READ5"] = (
        f"**READING.** (Candidate {q['chosen'] or 'M0'}; M0 next to it.) The dry minutes (no rain, outside the first hour after it) carry "
        f"{bud.share_of_squared_error.iloc[0]:.0%} of the squared error at {bud.rmse_db.iloc[0]:.2f} dB RMSE: the noise floor (baseline noise, "
        "quantization, dew - box C and B). During rain: "
        f"{bud.share_of_squared_error.iloc[1]:.0%} at {bud.rmse_db.iloc[1]:.2f} dB, and almost all of it WITHIN events "
        f"({bud.share_of_squared_error.iloc[4]:.0%}) rather than as one offset per event ({bud.share_of_squared_error.iloc[3]:.0%}): the rain-reference "
        "error that dominates the parameter intervals is a small part of the prediction error; the bulk is the 5 min, 2 km radar against a 1 min "
        "link inside an event (box C). After the rain (the dry path, first hour): "
        f"{bud.share_of_squared_error.iloc[2]:.0%} at {bud.rmse_db.iloc[2]:.2f} dB (M0 {bud['rmse_db M0'].iloc[2]:.2f} dB): the misfit S6b targets is "
        "real but a small part of the total.")
    if q["chosen"] and bt is not None:
        b = bt
        lines = ["```", f"WAA candidate for S8 (S6b): the S6 candidate with tau_f_dry free ({q['chosen']}); everything else as S6",
                 "A(t) = W_f(t) + W_s(t) + L * c_band * R(t - d) + c,   R = radar / g_e,   d = +1 min",
                 "dW_f/dt = (phi a_link A_band G(R) - W_f) / tau_f,   tau_f = 2.6 min while wetting (fixed), tau_f_dry while drying",
                 "dW_s/dt = ((1 - phi) a_link A_band G(R) - W_s) / tau_s,   phi 0.72, tau_s 240 min (fixed, not identified)",
                 f"CHANGED  tau_f_dry   {b.value['tau_0']:.2f} min [{b.lo['tau_0']:.2f}, {b.hi['tau_0']:.2f}]   (S6: 10.3 min fixed, the gauge timing)"]
        for nm in ("28 GHz", "32 GHz", "38 GHz"):
            lines.append(f"         A[{nm}] {b.value[f'A[{nm}]']:.2f} [{b.lo[f'A[{nm}]']:.2f}, {b.hi[f'A[{nm}]']:.2f}] dB    c[{nm}] {b.value[f'c[{nm}]']:.3f} "
                         f"[{b.lo[f'c[{nm}]']:.3f}, {b.hi[f'c[{nm}]']:.3f}] dB/km per mm/h{'   NOT identified' if not b.identified[f'c[{nm}]'] else ''}")
        lines.append("         G(0.5, 1, 2, 10, 20, 40) " + ", ".join(f"{b.value[f'G({k})']:.2f}" for k in ("0.5", "1", "2", "10", "20", "40")) + "   (G(5) = 1)")
        lines.append(f"         a_link 10-90 %  {b.value['a_p10']:.2f} - {b.value['a_p90']:.2f}")
        lines += ["[ ] = 90 % cold event bootstrap (40 draws)", "```"]
        bx = q["bx"]
        R["BOX"] = ("\n".join(lines) + "\n\n**Box A (in the equation):** " + "; ".join(bx[bx.box.str.startswith("A")].effect) + ".\n\n**Box B (identified, outside):** "
                    + "; ".join(bx[bx.box.str.startswith("B")].effect) + ".\n\n**Box C (unexplained, in the error budget):** " + "; ".join(bx[bx.box.str.startswith("C")].effect)
                    + ".\n\n**A DECISION IS NEEDED BEFORE S8.** By the rules written before the run, M1 enters box A and is the candidate written here. "
                    "But M1 makes the path term unidentified (c_band reaches 0 in the bootstrap), which the success definition of the roadmap ('box A numbers "
                    "identified and bounded') and integrity rule 3 do not allow for a result. The two ways forward, for the learners to choose and date: "
                    "(1) take M1 as written (c_band reported as not identified; W carries more of the attenuation); (2) keep the S6 candidate (M0, sha256 "
                    f"{q['sha6'][:8]}...) for S8 and put 'a faster effective drying with the radar input removes the end-of-rain misfit but unbinds the path "
                    "term' in box B. No further return to S6 is allowed either way.")
    R["READ7"] = ("1. Freeze (hash) the candidate the learners choose (section 6), then read TEST once.\n"
                  "2. TEST against VAL per length class (RMSE, the S5a plan), the TEST length test (at least 4 of 6 cells), V1 with g_e from the TEST events' "
                  "gauges, V0 and the gauge links as the independent check.\n"
                  "3. The dry-path error at 5 and 10 min on TEST, M0 and M1 side by side (the after-rain misfit).\n"
                  "4. If M1: whether c_band stays above 0 on TEST (the length dependence at matched rain is the only thing that pins it).")
    R["READ8"] = ("- **M1 passes the written rules but unbinds the path term** (section 3): the rules for a refinement did not ask that the OTHER numbers stay "
                  "identified; the repeated S6 check shows they do not. Reported, not decided here.\n"
                  "- **tau_f_dry = 3.7 min with the radar input against 10.3 / 8.8 min with a gauge**: an effective time, not the antenna's drying.\n"
                  "- **Wind**: S7 found it per event (17 events); at the window level it does not survive holding the preceding rain fixed, and M2 fails two "
                  "rules: box B, not A.\n"
                  "- **The deviation**: the projected cuts (bootstrap 40 -> 30, M3) were undone before Part 2 ran (config `return_s6b.deviation`).\n"
                  "- **Heavy rain on VAL** rests on 3 events (63 windows); rule 3's heavy part is a point estimate.\n"
                  "- **(d) on VAL**: the model is LOW during >= 10 mm/h on short VAL links (+1.2 dB, 4 events): F at heavy rain may be too LOW there, the "
                  "opposite of the S7 suspicion; 4 events.\n"
                  "- **a_link percentiles** of the fit lie outside their bootstrap interval (duplicated events in a draw), as in S5a and S6.\n"
                  "- **VAL is used (optimistic)**: every rule of Part 3 was judged on it.")
    return R


# =================================================================================================== report
def write_report(q):
    L = []
    a = L.append
    a("# wet_antenna step S6b: the one return to S6 (drying and wind; the after-rain misfit) and the three boxes\n")
    a("Written by `research/wet_antenna/return_s6b.py` (libraries `wind_drying.py`, `validation.py`, `final_model.py`, `system_id.py`, `length_test.py`). Every rule: "
      "`wet_antenna.return_s6b` in `research/config.yaml`, written on 2026-10-19 BEFORE the run. This is the ONE allowed return to S6 (S7 reason (d)); after it the "
      "candidate goes to S8 whatever the outcome.\n")
    a(STANDING + "\n")
    a(f"Reference M0 = the S6 candidate, sha256 `{q['sha6']}`. VAL numbers are **{VL}**.\n")
    a("## 1. What was run\n")
    pl = q["plan"]
    a(f"- Dry-path windows (the S7 rule): {len(q['win'])} (TRAIN {int((q['win'].set == 'TRAIN').sum())}, VAL {int((q['win'].set == 'VAL').sum())}).")
    a(f"- Fits on TRAIN V1 (kernel `wind_drying.kernel_w`, exact M4b at k = 0): {', '.join(q['refs'])}; fake checks {len(q['fjobs'])}; the box-A rule; "
      f"{'a cold bootstrap of ' + str(pl['boot']) + ' draws and 3 identifiability replicates for ' + q['chosen'] if q['chosen'] else 'no bootstrap (nothing passed box A)'}.")
    a(f"- Wind: hourly SMHI GbgA, standardised over the TRAIN link-event minutes (mean {q['wstats'][0]:.2f} m/s, SD {q['wstats'][1]:.2f} m/s).")
    a(f"- Timing: one evaluation {pl['t_eval_s'] * 1e3:.0f} ms on {pl['workers']} processes; projection {pl['projected_full_h']} h; cuts: {', '.join(pl['cuts']) or 'nothing'}; "
      f"M3 {'run' if 'M3' in q['refs'] else 'not run (Part 1 does not ask for it)' if pl['m3'] else 'cut'}. Run time {(time.time() - q['clock']) / 60:.0f} min.\n")
    a("## 2. The after-rain misfit: cause by cause (M0 frozen)\n")
    a(md_table(q["dg"]["table"]) + "\n")
    a("Weather and preceding rain, Spearman and partial Spearman with 90 % event-bootstrap intervals (`physical_and_excludes_0` = the rule):\n")
    a(md_table(q["dg"]["weather"].round(3)) + "\n")
    a("![causes](../figures/return_s6b_after_rain_causes.png)\n")
    a(q.get("READ2", "") + "\n")
    a("## 3. The refinements and the box-A rule\n")
    a(f"Fits (`{rel(P.S6B_FITS)}`):\n")
    a(md_table(pd.DataFrame(q["fit_rows"]).round(4)) + "\n")
    a(f"Fake recovery of the extra numbers (`{rel(P.S6B_FAKE)}`):\n")
    a(md_table(q["fkt"].round(4)) + "\n")
    a("![fake](../figures/return_s6b_fake.png)\n")
    a(f"The box-A rule (`{rel(P.S6B_RULE)}`):\n")
    a(md_table(q["rule_tab"]) + "\n")
    if len(q["fold_rows"]):
        a(md_table(pd.DataFrame(q["fold_rows"]).round(4)) + "\n")
    a("![dry path](../figures/return_s6b_dry_path.png)\n")
    if len(q["boot_tab"]):
        a(f"The new candidate's cold bootstrap ({q['plan']['boot']} draws) and identifiability (`{rel(P.S6B_BOOT)}`):\n")
        a(md_table(q["boot_tab"].round(4)) + "\n")
        it = q["ident_tab"].groupby("quantity").within.all()
        idn = q["boot_tab"].set_index("quantity").identified
        a(f"Identifiability: {int(it.sum())} of {len(it)} quantities within +-20 % in 3 of 3 on fake data; by the full S6 rule (also the bootstrap interval "
          f"not at a bound) {int(idn.sum())} of {len(idn)} are identified; NOT identified: {', '.join(idn[~idn].index) or 'none'}.\n")
    a(q.get("READ3", "") + "\n")
    a("## 4. What the model claims: the three boxes\n")
    a(md_table(q["bx"]) + "\n")
    a("![boxes](../figures/return_s6b_boxes.png)\n")
    a("## 5. The error budget (VAL, " + VL + ")\n")
    a(md_table(q["bud"].round(4)) + "\n")
    a("![budget](../figures/return_s6b_error_budget.png)\n")
    a(q.get("READ5", "") + "\n")
    a("## 6. The candidate for S8\n")
    a(q.get("BOX", "") + "\n")
    a(f"`{rel(P.S6B_CANDIDATE)}`, sha256 **`{q['sha']}`**.\n")
    a(q.get("READ6", "") + "\n")
    a("## 7. What S8 must check\n")
    a(q.get("READ7", "") + "\n")
    a("## 8. Unclear / possibly wrong\n")
    a(q.get("READ8", "") + "\n")
    return "\n".join(L)


if __name__ == "__main__":
    main()
