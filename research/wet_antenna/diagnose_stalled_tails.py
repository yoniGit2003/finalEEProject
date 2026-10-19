"""wet_antenna step 4b: why do some real tails stall or rise? Contamination or antenna physics?

Step 4 found no drying law: about one real tail in five is still near its start level after 10 min, and the drying
rate per dB, k(W) = -(1/W) dW/dt, FALLS with W. Two explanations are kept open and tested here by their DIFFERENT
predictions, mostly without fitting any model:
  H1 contamination     something still adds attenuation after t0_trace (rain or drizzle below the references, dew, a
                       wrong end of the rain, a baseline shift). It belongs to the weather after the rain.
  H2 antenna physics   a drying law of the antenna itself whose rate per dB falls with A (saturation, drops that stay,
                       a delayed start of evaporation). It belongs to the antenna and the water on it.
What the script does (TRAIN tails of the frozen set v1 only; the 123 VAL tails are never read):
  classes   every tail gets one of four classes from the shape of its first 60 min (rules in config, fixed first);
  part A    outside evidence of continued wetting per tail (radar on the path and nearby, gauges, humidity, dew,
            time of day, wind), the same classes in the other tails of the event, the two directions of a link, and
            whether the rise grows with the path length;
  part B    B1 transition level, B2 the same link in different events, B3 plateau length vs wetness, B4 do all tails
            lie on one curve shifted in time (master curve), B5 the model-free k(A) per subset;
  part C    the two fixes of pooled_fit.py found in step 4 (C1: the start of F1 / F4 bounded like the other laws,
            C2: one fitted multiplier on the noise variance), a check on 6 cached fake datasets, and the menu fits
            and free curve of step 4 again with both fixes on four sets.
Tails are never chosen by how well a law fits them. One fit is timed first and the projection decides the cuts.
Fits are cached in research/.cache/wet_antenna/diagnose_stalled_tails/ (delete the folder to fit again).
Run from the repository root (about 40 min the first time on 3 processes; about 9 min when the fits are cached,
most of it the event bootstraps of parts A and B):
    python -m research.wet_antenna.diagnose_stalled_tails
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "1")            # one thread per process: the processes run in parallel

from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cmlpipe.api import CMLData
from research.common.config import load_config
from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import datasets as dsets
from research.wet_antenna import fit_jobs as fj
from research.wet_antenna import free_curve as wc
from research.wet_antenna import paths as P
from research.wet_antenna import pooled_fit as wp
from research.wet_antenna.tails import classify as tails_classify, raw_rate, spearman

CFG = load_config()
W = CFG["wet_antenna"]
S = W["diagnose_stalled_tails"]
CL, PA, PB, PC, MC = S["classes"], S["part_a"], S["part_b"], S["part_c"], S["part_b"]["master"]
S4, S1C, FIT, S3B = W["fit_real_drying"], W["stress_pooled_fit"], W["check_free_curve"], W["confirm_free_curve"]
WIN, W_LO = S["window_min"], S3B["w_lo_db"]
CACHE = P.CACHE / "diagnose_stalled_tails"
CACHE_BEFORE = P.CACHE / "check_discovery_method"          # the menu fits of the fake datasets before the fixes
FAKE_DIR = P.CACHE / "confirm_free_curve"                  # the fake datasets of step 3b-3b
CLASSES = ["decays", "plateau-then-decays", "rises", "flat / slow"]
NOT_DECAYING = CLASSES[1:]
FAMILY = {"F1": "exponential", "F3": "nonlinear in W (F3 / F5)", "F5": "nonlinear in W (F3 / F5)",
          "F4": "two layers", "F6": "state switch"}


# =================================================================================================== classes
def classify(a):
    """Class of one tail from its first 60 minutes of A (1 min values, NaN allowed). Rules: config `classes`.
    The code is in tails.py, because step 4c uses the same rules."""
    return tails_classify(a, CL)


# =================================================================================================== small statistics
def event_boot(df, stat, rng, reps):
    """Event bootstrap: draw whole events with replacement, compute `stat` on the tails of the draw. Returns the
    interval of config part_a.interval (90 %)."""
    groups = [g for _, g in df.groupby("event_id")]
    vals = [stat(pd.concat([groups[i].assign(draw=k) for k, i in enumerate(rng.integers(0, len(groups), len(groups)))]))
            for _ in range(reps)]                          # `draw` tells two copies of the same event apart
    vals = np.asarray(vals, float)
    return (np.nan, np.nan) if np.isnan(vals).mean() > 0.5 else tuple(np.nanpercentile(vals, PA["interval"]))


def odds_ratio(df, col):
    """Odds of `col` being true in the non-decaying tails over the odds in the decaying tails (0.5 added to every
    cell, so an empty cell does not give 0 or infinity). Tails where `col` is unknown are left out."""
    g = df[df[col].notna()]
    x, nd = g[col].astype(bool), g.not_decaying
    a, b, c, e = (x & nd).sum(), (~x & nd).sum(), (x & ~nd).sum(), (~x & ~nd).sum()
    return (a + 0.5) * (e + 0.5) / ((b + 0.5) * (c + 0.5))


def rho_with_interval(df, x, y, rng):
    g = df[[x, y, "event_id"]].dropna()
    if len(g) < 8:
        return np.nan, np.nan, np.nan, len(g)
    lo, hi = event_boot(g, lambda q: spearman(q[x], q[y]) if len(q) > 3 else np.nan, rng, PA["boot_reps"])
    return float(spearman(g[x], g[y])), lo, hi, len(g)


# =================================================================================================== part A
def outside_evidence(d, ix):
    """One row per tail: what the references say from t0_trace to t0_trace + 60 min (and 30 min before)."""
    links = d.links().set_index("cml_id")
    sites = d.table("sites").set_index("site_id")
    mid = pd.DataFrame({"x": (sites.x_3006_m.reindex(links.site_0_id).to_numpy() + sites.x_3006_m.reindex(links.site_1_id).to_numpy()) / 2,
                        "y": (sites.y_3006_m.reindex(links.site_0_id).to_numpy() + sites.y_3006_m.reindex(links.site_1_id).to_numpy()) / 2},
                       index=links.index)
    ra = d.store("radar_along_links_5min").R_max.load().to_pandas()       # (time, cml_id), mm/h, 5 min, interval start
    T5, R, col = pd.DatetimeIndex(ra.index), ra.to_numpy(), {c: i for i, c in enumerate(ra.columns)}
    xy = mid.reindex(ra.columns).to_numpy()
    near = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1]) < PA["near_km"] * 1000
    np.fill_diagonal(near, False)
    city = d.store("gauges_city_1min").amount_mm.load().to_pandas()       # (time, gauge_id), mm per minute
    smhi = d.store("gauges_smhi_15min").amount_mm.load().to_pandas().iloc[:, 0]
    met = pd.read_csv(P.DEW_RISK_HOURLY, parse_dates=["time_utc"]).set_index("time_utc")
    hour, before = pd.Timedelta("1h"), pd.Timedelta(minutes=S["context_before_min"])
    out = []
    for r in ix.itertuples():
        t0, t1 = r.t0_trace, r.t0_trace + pd.Timedelta(minutes=WIN)
        i0, i1, ib = T5.searchsorted(t0), T5.searchsorted(t1), T5.searchsorted(t0 - before)
        path, pre = R[i0:i1, col[r.cml_id]], R[ib:i0, col[r.cml_id]]
        nb = R[i0:i1][:, near[col[r.cml_id]]]
        g_id = links.nearest_city_gauge_id[r.cml_id]
        g = city.loc[t0:t1 - pd.Timedelta("1s"), g_id]
        s_ = smhi.loc[t0 - pd.Timedelta("14min"):t1 - pd.Timedelta("1s")]
        wx = met.reindex([t0.floor("1h"), (t1 - pd.Timedelta("1min")).floor("1h")])
        local = t0 + pd.Timedelta(hours=2)                               # Gothenburg summer time = UTC + 2 h
        row = {"tail_id": r.tail_id,
               "radar_path_max_mmh": float(np.nanmax(path)) if np.isfinite(path).any() else np.nan,
               "radar_path_nan_bins": int(np.isnan(path).sum()), "radar_path_bins": len(path),
               "radar_before_max_mmh": float(np.nanmax(pre)) if np.isfinite(pre).any() else np.nan,
               "radar_near_links": int(near[col[r.cml_id]].sum()),
               "radar_near_max_mmh": float(np.nanmax(nb)) if np.isfinite(nb).any() else np.nan,
               "radar_near_wet_share": float(np.nanmean(nb > 0)) if nb.size else np.nan,
               "city_gauge_id": g_id, "city_gauge_km": round(links.nearest_city_gauge_distance_m[r.cml_id] / 1000, 2),
               "city_gauge_mm": float(g.sum()) if g.notna().any() else np.nan,
               "smhi_gauge_km": round(links.nearest_smhi_gauge_distance_m[r.cml_id] / 1000, 2),
               "smhi_gauge_mm": float(s_.sum()) if s_.notna().any() else np.nan,
               "humidity_pct": float(wx.humidity.mean()), "temperature_c": float(wx.temperature.mean()),
               "dew_point_depression_c": float(wx.dew_point_depression.min()),
               "wind_speed_ms": float(wx.wind_speed.mean()), "hour_local": local.hour + local.minute / 60,
               "night": bool(wx.night.iloc[0]) if pd.notna(wx.night.iloc[0]) else np.nan,
               "sun_elevation_deg": float(wx.sun_elevation_deg.iloc[0])}
        flag = lambda v, known=True: bool(v) if known else np.nan        # noqa: E731
        row["ev_radar_path"] = flag(row["radar_path_max_mmh"] > 0, np.isfinite(row["radar_path_max_mmh"]))
        row["ev_radar_near"] = flag(row["radar_near_max_mmh"] > 0, np.isfinite(row["radar_near_max_mmh"]))
        row["ev_gauge_city"] = flag(row["city_gauge_mm"] > 0, np.isfinite(row["city_gauge_mm"]))
        row["ev_gauge_smhi"] = flag(row["smhi_gauge_mm"] > 0, np.isfinite(row["smhi_gauge_mm"]))
        # the hourly product says "no dew risk" where the humidity is missing: that is "unknown" here, not "no"
        row["ev_dew_risk"] = flag(wx.dew_risk.fillna(False).astype(bool).any(), wx.humidity.notna().all())
        rain = [row[k] for k in ("ev_radar_path", "ev_radar_near", "ev_gauge_city", "ev_gauge_smhi")]
        row["ev_any_rain"] = any(v is True for v in rain) if any(v is True for v in rain) or all(v is False for v in rain) else np.nan
        both = [row["ev_any_rain"], row["ev_dew_risk"]]
        row["ev_any_moisture"] = any(v is True for v in both) if any(v is True for v in both) or all(v is False for v in both) else np.nan
        out.append(row)
    return pd.DataFrame(out)


def coherence(T):
    """Per tail: share of the event's OTHER tails (other links, start within overlap_min) that are 'rises' or
    'flat / slow'. NaN if there is no such tail."""
    share, n = np.full(len(T), np.nan), np.zeros(len(T), int)
    bad = T["class"].isin(["rises", "flat / slow"]).to_numpy()
    t0, link = T.t0_trace.to_numpy(), T.cml_id.to_numpy()
    for _, idx in T.groupby("event_id").indices.items():
        for i in idx:
            o = idx[(link[idx] != link[i]) & (np.abs(t0[idx] - t0[i]) <= np.timedelta64(PA["overlap_min"], "m"))]
            n[i] = len(o)
            if len(o):
                share[i] = bad[o].mean()
    return share, n


# =================================================================================================== part B
def link_repeat_test(T, rng):
    """B2. Unit = (link, event): y = share of the link's tails in that event that do NOT decay, minus the mean of the
    event (so a 'bad' event does not count as a property of its links). Statistic = mean product of y over the pairs
    of events of the same link (links with >= 2 events): positive when a link behaves the same way in different
    events. Null: the y values shuffled among the links INSIDE each event."""
    u = T.groupby(["cml_id", "event_id"]).not_decaying.mean().rename("y").reset_index()
    u["y"] = u.y - u.groupby("event_id").y.transform("mean")
    keep = u.groupby("cml_id").event_id.transform("size") >= PB["min_events_per_link"]

    def stat(y):
        g = pd.DataFrame({"l": u.cml_id[keep], "y": y[keep]}).groupby("l").y
        s1, s2, n = g.sum(), g.apply(lambda v: float((v ** 2).sum())), g.size()
        return float(((s1 ** 2 - s2) / 2).sum() / (n * (n - 1) / 2).sum())

    obs = stat(u.y)
    null = np.array([stat(u.groupby("event_id").y.transform(lambda v: rng.permutation(v.to_numpy()))) for _ in range(PB["shuffles"])])
    return obs, float((np.sum(null >= obs) + 1) / (len(null) + 1)), int(u[keep].cml_id.nunique()), int(keep.sum()), float(np.std(null))


class Master:
    """B4. A law of the form dA/dt = -f_i h(A) (nothing but A decides the rate) makes every tail the SAME curve,
    shifted in time and stretched by its speed: A_i(t) = M(f_i t + d_i). M is built as the median of the smoothed
    'decays' tails after shifting (and stretching) them onto each other; no law is fitted."""

    def __init__(self):
        self.tau = np.arange(MC["tau_range"][0], MC["tau_range"][1] + 1.0)
        self.fs = np.geomspace(MC["stretch_grid"][0], MC["stretch_grid"][1], int(MC["stretch_grid"][2]))
        self.ds = np.arange(MC["tau_range"][0], MC["tau_range"][1] + 1.0, MC["shift_step_min"])
        self.M = np.full(len(self.tau), np.nan)

    def errors(self, t, y, fs):
        """Mean squared misfit of a tail (minutes t, values y) for every stretch in fs and every shift: (fs, shifts)."""
        E = np.full((len(fs), len(self.ds)), np.inf)
        for i, f in enumerate(fs):
            q = f * t[None, :] + self.ds[:, None]
            m = np.interp(q.ravel(), self.tau, self.M, left=np.nan, right=np.nan).reshape(q.shape)
            ok = np.isfinite(m)
            n = ok.sum(axis=1)
            e = np.where(ok, (y[None, :] - np.where(ok, m, 0.0)) ** 2, 0.0).sum(axis=1) / np.maximum(n, 1)
            e[(n < MC["min_overlap_min"]) | (n < MC["min_cover"] * len(t))] = np.inf
            E[i] = e
        return E

    def place(self, t, y, fs=None):
        """Best (f, d, rms) of one tail on the master; (nan, nan, nan) if it cannot be placed."""
        fs = self.fs if fs is None else fs
        if len(t) < MC["min_overlap_min"]:
            return np.nan, np.nan, np.nan
        E = self.errors(t, y, fs)
        i, j = np.unravel_index(np.argmin(E), E.shape)
        return (float(fs[i]), float(self.ds[j]), float(np.sqrt(E[i, j]))) if np.isfinite(E[i, j]) else (np.nan, np.nan, np.nan)

    def rebuild(self, series, placed):
        """Median of the placed tails on the tau grid. series: {tail: (t, y)}; placed: {tail: (f, d)}."""
        rows = []
        for k, (f, d_) in placed.items():
            t, y = series[k]
            if np.isfinite(f) and len(t) > 1:
                rows.append(np.interp((self.tau - d_) / f, t, y, left=np.nan, right=np.nan))
        rows = np.array(rows)
        n = np.isfinite(rows).sum(axis=0)
        with np.errstate(all="ignore"):
            self.M = np.where(n >= MC["min_tails"], np.nanmedian(rows, axis=0), np.nan)
        self.n = n

    def value(self, tau):
        return np.interp(tau, self.tau, self.M, left=np.nan, right=np.nan)


def build_master(series, T, per_link):
    """The master curve from the 'decays' tails of the events NOT flagged rise-heavy. First alignment: the minute at
    which a tail falls below ref_db is tau = 0, no stretch. Then `iterations` rounds: place every tail on the master
    (one stretch per link, or one shared stretch = 1 if per_link is False; its own shift), rebuild the median."""
    M = Master()
    base = T[(T["class"] == "decays") & ~T.event_rise_heavy]
    placed = {}
    for r in base.itertuples():
        t, y = series[r.tail_id]
        below = np.nonzero(y < MC["ref_db"])[0]
        if len(below) and below[0] > 0:
            placed[r.tail_id] = (1.0, -float(t[below[0]]))
    M.rebuild(series, placed)
    first = len(placed)
    for _ in range(MC["iterations"]):
        placed = {}
        for _, g in base.groupby("cml_id"):
            ids = [k for k in g.tail_id if len(series[k][0]) >= MC["min_overlap_min"]]
            if not ids:
                continue
            fs = M.fs if per_link else np.array([1.0])
            E = {k: M.errors(*series[k], fs) for k in ids}
            tot = sum(np.where(np.isfinite(E[k].min(axis=1)), E[k].min(axis=1), 1e3) for k in ids)
            i = int(np.argmin(tot))
            for k in ids:
                j = int(np.argmin(E[k][i]))
                if np.isfinite(E[k][i, j]):
                    placed[k] = (float(fs[i]), float(M.ds[j]))
        g_ = float(np.exp(np.mean([np.log(f) for f, _ in placed.values()])))   # keep the median speed at 1
        placed = {k: (f / g_, d_ / g_) for k, (f, d_) in placed.items()}
        M.rebuild(series, placed)
    return M, placed, len(base), first


# =================================================================================================== part C: workers
def noise_tables():
    """Noise and offset tables of the fit (a sub-link without control tails would get the median; none in v1 main)."""
    return fj.noise_tables()


def pack_file(file, window, min_points):
    noise, offset_sd = noise_tables()
    tails = dsets.load_tails(file)
    D = wp.pack(tails, noise, offset_sd, window, min_points)
    used = set(D["tail_id"])
    return tails, D, len({(tl["cml_id"], tl["event_copy"]) for tl in tails if tl["tail_id"] in used})


def menu_starts(cand):
    return list(S1C["t2_starts"][cand]) + ([] if cand == "F1" else ["f1sol"])


def menu_file(name, cand, label):
    return CACHE / f"menu_{name}_{cand}_{label}.json"


def run_menu(job):
    """One candidate on one dataset from every start of step 4 (the generic starts of step 3b-1c + the start built
    from the F1 solution), with BOTH fixes on: bound_w0 (C1) and the fitted noise multiplier (C2)."""
    name, cand, data, out = job["name"], job["cand"], {}, []
    for label in menu_starts(cand):
        f_ = menu_file(name, cand, label)
        if f_.exists():
            out.append(json.loads(f_.read_text()))
            continue
        if not data:
            _, data["D"], data["n"] = pack_file(job["file"], job["window"], job["min_points"])
        if label == "f1sol":
            f1 = max((json.loads(menu_file(name, "F1", lb).read_text()) for lb in menu_starts("F1")), key=lambda r: r["ll"])
            theta0, s0 = wp.f1_start(cand, f1["theta"][0]), f1["s"]
        else:
            theta0, s0 = S1C["t2_starts"][cand][label]
        t0 = time.time()
        r = wp.fit(cand, data["D"], theta0, s0, S1C["maxfev_generic"], bound_w0=True, fit_noise=True)
        r.update({"dataset": name, "candidate": cand, "start": label, "n_tails": len(data["D"]["tail_id"]),
                  "n_independent": data["n"], "n_points": int(len(data["D"]["T"])), "k": len(wp.NAMES[cand]) + 2,
                  "seconds": round(time.time() - t0, 1)})
        fj.save_json(f_, r)
        out.append(r)
    return out


def run_curve(job):
    """The free-curve pooled fit of step 4 (b) on one dataset, with the fitted noise multiplier (the free curve always
    had the start bound). Same stages: 4 starts, the best continued, the other penalties, penalty by 2-fold
    cross-validation over events. x = [log k at the 5 knots, log s, log m]."""
    f_ = CACHE / f"curve_{job['name']}.json"
    res = json.loads(f_.read_text()) if f_.exists() else {}
    if res.get("done"):
        return res
    t0 = time.time()
    tails, D, n_ind = pack_file(job["file"], WIN, FIT["min_points"])
    noise, offset_sd = noise_tables()
    cv = wc.Curve(job["knots"], prior_sd=FIT["prior_sd_logk"])
    cv.fit_noise = True
    lam0, mf = FIT["lambda_start"], FIT["maxfev"]
    res.update({"dataset": job["name"], "knots": job["knots"], "w_hi": job["w_hi"], "n_tails": len(D["tail_id"]),
                "n_independent": n_ind, "n_points": int(len(D["T"]))})
    res.setdefault("starts", {})
    for name, spec in FIT["starts"].items():
        if name not in res["starts"]:
            res["starts"][name] = cv.fit(D, np.r_[wc.start_theta(spec, cv.knots), np.log(FIT["s_start"]), 0.0], lam0, mf["start"])
            fj.save_json(f_, res)
    res.setdefault("lam", {})
    if str(lam0) not in res["lam"]:
        best = min(res["starts"], key=lambda n: res["starts"][n]["cost"])
        res["lam"][str(lam0)] = {**cv.fit(D, res["starts"][best]["x"], lam0, mf["finish"]), "from": best}
        fj.save_json(f_, res)
    for lam in FIT["lambda_grid"]:
        if str(lam) not in res["lam"]:
            res["lam"][str(lam)] = cv.fit(D, res["lam"][str(lam0)]["x"], lam, mf["other_lambda"])
            fj.save_json(f_, res)
    res.setdefault("cv", {})
    evs = np.random.default_rng(PC["cv_seed"]).permutation(sorted({tl["event_copy"] for tl in tails}))
    half, Dh = [set(evs[::2].tolist()), set(evs[1::2].tolist())], None
    for lam in FIT["lambda_grid"]:
        if str(lam) not in res["cv"]:
            if Dh is None:
                Dh = [wp.pack([tl for tl in tails if tl["event_copy"] in h], noise, offset_sd, WIN, FIT["min_points"]) for h in half]
            res["cv"][str(lam)] = float(sum(cv.loglik(np.asarray(cv.fit(Dh[a], res["lam"][str(lam)]["x"], lam, mf["cv"])["x"]), Dh[b])
                                            for a, b in ((0, 1), (1, 0))))
            fj.save_json(f_, res)
    lam_star = max(FIT["lambda_grid"], key=lambda v: res["cv"][str(v)])
    fit = res["lam"][str(lam_star)]
    res.update({"lam_star": lam_star, "x": fit["x"], "ll": fit["ll"], "s": float(np.exp(fit["x"][-2])),
                "noise_mult": float(np.exp(fit["x"][-1])), "bic": -2 * fit["ll"] + (cv.K + 2) * np.log(n_ind),
                "seconds": round(res.get("seconds", 0.0) + time.time() - t0, 1), "done": True})
    fj.save_json(f_, res)
    return res


def run_any(job):
    return {"menu": run_menu, "curve": run_curve}[job["part"]](job)


def read_menu(best):
    """The reading rules of step 3b-1c on the best fit of every candidate (as step 4): F3 and F5 are one family; F4
    and F6 closer than f4_f6_bic are not separated."""
    b = best.set_index("candidate").bic.sort_values()
    win = b.index[0]
    fam, members = FAMILY[win], {c for c in b.index if FAMILY[c] == FAMILY[win]}
    if win in ("F4", "F6") and b["F6" if win == "F4" else "F4"] - b[win] < S4["f4_f6_bic"]:
        fam, members = "two layers or state switch (F4 / F6 not separated)", {"F4", "F6"}
    rest = b[[c for c in b.index if c not in members]]
    return {"winner": win, "family": fam, "runner_up": b.index[1], "margin": float(b.iloc[1] - b.iloc[0]),
            "next_family": rest.index[0], "margin_family": float(rest.iloc[0] - b.iloc[0]), "f1_minus_winner": float(b["F1"] - b.iloc[0])}


def read_slope(slope):
    R = S4["reading"]
    if abs(slope) < R["flat_abs_slope"]:
        return "exponential in this range"
    if slope > R["quadratic_slope"]:
        return "quadratic family in this range"
    return "mildly rising, law not identified by the curve" if slope > 0 else "falling (no fake truth looked like this)"


# =================================================================================================== figures
def gallery_figure(d, T, series_raw, picks):
    """4 tails per class: A (dots = 1 min, line = 5 min median) from 30 min before t0_trace, and under it the radar on
    the path and the nearest City gauge. Everything is labelled in words; line styles, not colours, tell the lines apart."""
    ra, city = d.store("radar_along_links_5min").R_max, d.store("gauges_city_1min").rate_mm_h
    n = S["gallery"]["per_class"]
    fig = plt.figure(figsize=(4.1 * n, 12.5))
    gs = fig.add_gridspec(len(CLASSES) * 2, n, height_ratios=[3, 1] * len(CLASSES), hspace=0.62, wspace=0.22,
                          top=0.93, bottom=0.05, left=0.05, right=0.97)
    before = pd.Timedelta(minutes=S["context_before_min"])
    for i, cls in enumerate(CLASSES):
        for j, tid in enumerate(picks[cls]):
            r = T.loc[T.tail_id == tid].iloc[0]
            t0, t1 = r.t0_trace - before, r.t0_trace + pd.Timedelta(minutes=WIN)
            ax, bx = fig.add_subplot(gs[2 * i, j]), fig.add_subplot(gs[2 * i + 1, j])
            pre = d.cml([r.sublink_id], t0, r.t0_trace, res="1min", variables=["tl_mean"]).tl_mean.to_pandas().iloc[:, 0] - r.baseline_db
            ax.plot((pre.index - r.t0_trace) / pd.Timedelta("1min"), pre, "x", ms=3, color="#9a9a9a", label="A before t0_trace (crosses)")
            a = series_raw[tid]
            ax.plot(np.arange(WIN), a, "o", ms=2.6, color="#8a8a8a", label="A, 1 min (dots)")
            ax.plot(np.arange(WIN), pd.Series(a).rolling(CL["smooth_min"], center=True, min_periods=3).median(), "-", color=INK, lw=1.8,
                    label="A, 5 min median (solid)")
            for lev, ls, txt in ((1.0, "--", "A_start"), (CL["plateau_level"], ":", "0.8 x"), (CL["half"], ":", "0.5 x")):
                ax.axhline(lev * r.A_start_db, ls=ls, color=INK2, lw=0.9)
                ax.text(WIN + 0.5, lev * r.A_start_db, txt, fontsize=7, va="center")
            ax.axvline(0, color=INK, lw=0.8)
            ax.set_xlim(-S["context_before_min"], WIN + 8)
            ax.set_title(f"{cls.upper()}\ntail {tid}, event {r.event_id[1:9]}, L = {r.L_km:.1f} km", fontsize=8)
            style_axes(ax)
            rr = ra.sel(cml_id=r.cml_id, time=slice(t0, t1)).to_pandas()
            gg = city.sel(gauge_id=r.city_gauge_id, time=slice(t0, t1)).to_pandas().rolling(5, min_periods=1).mean()
            bx.step((rr.index - r.t0_trace) / pd.Timedelta("1min"), rr.fillna(0), where="post", color=INK, lw=1.4, label="radar R_max on the path (solid steps)")
            bx.plot((gg.index - r.t0_trace) / pd.Timedelta("1min"), gg, "--", color=INK, lw=1.1,
                    label="nearest City gauge, 5 min mean (dashed)")
            bx.axvline(0, color=INK, lw=0.8)
            bx.set_xlim(-S["context_before_min"], WIN + 8)
            bx.set_ylim(0, max(1.0, float(np.nanmax(np.r_[rr.to_numpy(), gg.to_numpy(), 0.0])) * 1.15))
            bx.text(0.99, 0.9, f"rain, mm/h (gauge {r.city_gauge_km:.1f} km away)", transform=bx.transAxes, fontsize=7, ha="right", va="top")
            style_axes(bx)
            if j == 0:
                ax.set_ylabel("A (dB)")
            if i == 0 and j == 0:
                ax.legend(frameon=False, fontsize=6.5, loc="upper right")
                bx.legend(frameon=False, fontsize=6.5, loc="center right")
            if i == len(CLASSES) - 1:
                bx.set_xlabel("minutes since t0_trace")
    fig.suptitle(f"Step 4b: {n} random TRAIN tails per class (seed {S['gallery']['seed']}). Top: A; bottom: radar on the path and the "
                 "nearest City gauge. Vertical line: t0_trace.", fontsize=10, y=0.985)
    return fig


def master_figure(M, series, T, placed, plateau_fit):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4), sharex=True, sharey=True)
    ax = axes[0]
    for k, (f, d_) in placed.items():
        t, y = series[k]
        ax.plot(f * t + d_, y, "-", color="#b5b5b5", lw=0.6)
    ax.plot(M.tau, M.M, "-", color=INK, lw=2.8)
    ax.text(0.98, 0.95, f"thin grey: {len(placed)} 'decays' tails of the events\nnot flagged rise-heavy, shifted and stretched\n"
            "thick solid: their median = the master curve", transform=ax.transAxes, fontsize=8.5, ha="right", va="top")
    ax.set(xlabel="master time tau = f t + d (min)", ylabel="A, 5 min median (dB)", title="How the master curve is built")
    ax = axes[1]
    ax.plot(M.tau, M.M, "-", color=INK, lw=2.8)
    for k, v in plateau_fit.items():
        if not np.isfinite(v["f"]):
            continue
        t, y = series[k]
        p = t < v["plateau_min"]
        ax.plot(v["f"] * t[~p] + v["d"], y[~p], "-", color="#7a7a7a", lw=0.8)
        ax.plot(v["f"] * t[p] + v["d"], y[p], "o", ms=2.2, mfc="none", color=INK2)
    ax.text(0.98, 0.95, "thick solid: the master curve\nthin grey lines: the DECAYING part of each\n'plateau-then-decays' tail, placed on it\n"
            "open circles: the PLATEAU part of the same tail,\nwith the same shift and stretch", transform=ax.transAxes,
            fontsize=8.5, ha="right", va="top")
    ax.set(xlabel="master time tau = f t + d (min)", title="The plateau tails on the master curve")
    ok = np.isfinite(M.M)
    for a_ in axes:
        a_.set_xlim(M.tau[ok][0] - 5, M.tau[ok][-1] + 5)
        style_axes(a_)
    fig.suptitle("Step 4b, B4: if only A decided the drying rate, every tail would be the same curve shifted in time", fontsize=10)
    fig.tight_layout()
    return fig


def rate_figure(rate):
    fig, ax = plt.subplots(figsize=(8.8, 5.4))
    styles = [("-", "o"), ("--", "s"), ("-.", "^"), (":", "D"), ((0, (5, 1, 1, 1, 1, 1)), "v"), ((0, (1, 1)), "x"), ("-", "+")]
    for (name, g), (ls, mk) in zip(rate[~rate.subset.str.contains("minutes")].groupby("subset", sort=False), styles):
        g = g[g.used_for_slope]
        if len(g):
            ax.plot(g.median_A_db, g.k_from_mean_per_min, ls=ls, marker=mk, ms=5, color=INK, lw=1.5, mfc="none")
            ax.text(g.median_A_db.iloc[-1] * 1.04, g.k_from_mean_per_min.iloc[-1], f"{name} ({g.slope.iloc[0]:+.2f})", fontsize=7.5, va="center")
    ax.set(xscale="log", yscale="log", xlabel="A (dB), median of the bin", ylabel="k = -mean dA/dt / median A (1 / min)",
           xlim=(W_LO * 0.9, rate[rate.used_for_slope].median_A_db.max() * 3.2))
    style_axes(ax)
    ax.set_title("Step 4b, B5: the drying rate per dB at each level, without any model, per subset of tails\n"
                 "(label: subset and slope of log k vs log A; only bins with enough minutes)", fontsize=9.5)
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    clock = time.time()
    workers = max((os.cpu_count() or 2) - 1, 1)
    rng = np.random.default_rng(PA["seed"])

    # ------------------------------------------------------------------------------------------- data (TRAIN only)
    ix = pd.read_csv(P.TAILS_MODEL_INDEX, parse_dates=["t0", "t0_trace", "t_end"]).query("split_event == 'train'")
    ix = ix.sort_values("tail_id").reset_index(drop=True)
    rows = pq.read_table(P.TAILS_MODEL, filters=[("split_event", "=", "train")],
                         columns=["tail_id", "sublink_id", "cml_id", "event_id", "minutes_since_t0_trace", "A_db"]).to_pandas()
    assert set(rows.tail_id) == set(ix.tail_id) and (ix.split_event == "train").all()        # no VAL tail
    extra = pd.read_csv(P.TAILS_INDEX, usecols=["tail_id", "sublink_id", "event_id", "baseline_db", "wet_radar_R_mean_sum_mm", "wet_span_min"])
    ix = ix.merge(extra, on=["tail_id", "sublink_id", "event_id"], how="left", validate="1:1")
    assert ix.baseline_db.notna().all()
    first = rows[rows.minutes_since_t0_trace < WIN]
    A = first.pivot(index="tail_id", columns="minutes_since_t0_trace", values="A_db").reindex(index=ix.tail_id, columns=range(WIN)).to_numpy()
    print(f"TRAIN tails only: {len(ix)} tails, {ix.event_id.nunique()} events, {ix.cml_id.nunique()} links; {workers} processes")

    # ------------------------------------------------------------------------------------------- classes
    cls_rows, smooth = [], {}
    for tid, a in zip(ix.tail_id, A):
        c, sm = classify(a)
        cls_rows.append({"tail_id": tid, **c})
        smooth[tid] = sm
    keep_cols = ["tail_id", "sublink_id", "cml_id", "event_id", "event_rise_heavy", "level_shift", "undershoot", "t0_trace", "L_km",
                 "frequency_ghz", "A_start_trace_db", "baseline_db", "wet_radar_R_mean_sum_mm", "wet_span_min"]
    T = ix[keep_cols].merge(pd.DataFrame(cls_rows), on="tail_id")
    T["not_decaying"] = T["class"] != "decays"
    series = {tid: (np.nonzero(np.isfinite(sm))[0].astype(float), sm[np.isfinite(sm)]) for tid, sm in smooth.items()}
    count = lambda g: g["class"].value_counts().reindex(CLASSES, fill_value=0)   # noqa: E731
    counts = pd.concat([count(T).rename("all tails").to_frame().T,
                        T.groupby("event_rise_heavy").apply(count).rename(index={True: "events flagged rise-heavy", False: "other events"}),
                        T.groupby("event_id").apply(count)])
    counts.insert(0, "tails", counts.sum(axis=1))
    counts["share_not_decaying"] = (1 - counts["decays"] / counts.tails).round(3)
    counts = counts.rename_axis("group").reset_index()
    counts.insert(1, "rise_heavy", counts.group.map(T.groupby("event_id").event_rise_heavy.first()))
    save(counts, P.STALL_COUNTS, index=False)
    overlap = {"decays and plateau": int((T.meets_decays & T.meets_plateau).sum()),
               "plateau and rises": int((T.meets_plateau & T.meets_rises).sum()),
               "decays and rises": int((T.meets_decays & T.meets_rises).sum())}
    print("classes: " + ", ".join(f"{c} {n}" for c, n in count(T).items()), flush=True)

    # ------------------------------------------------------------------------------------------- time one fit, plan
    sets = {"main": ("v1 main", T.tail_id, {}),
            "no_rise_heavy": ("without event_rise_heavy", T.tail_id[~T.event_rise_heavy], {}),
            "decays": ("class 'decays' only", T.tail_id[T["class"] == "decays"], {}),
            "decays_plus": ("'decays' + 'plateau-then-decays' from the plateau's end",
                            T.tail_id[T["class"].isin(CLASSES[:2])],
                            T[T["class"] == "plateau-then-decays"].set_index("tail_id").min_to_80pct.to_dict())}
    ev_code = {e: i for i, e in enumerate(sorted(ix.event_id.unique()))}
    a_start_ix = ix.set_index("tail_id").A_start_trace_db
    meta = {}
    for name, (label, ids, restart) in sets.items():
        d_ = rows[rows.tail_id.isin(ids)].copy()
        d_["minutes_since_t0_trace"] = d_.minutes_since_t0_trace - d_.tail_id.map(restart).fillna(0).astype(int)   # restart = the new minute 0
        d_ = d_[(d_.minutes_since_t0_trace >= 0) & (d_.minutes_since_t0_trace < WIN)]
        d_["A_start_trace_db"] = d_.tail_id.map(a_start_ix)
        if restart:                                        # a restarted tail: its start level is measured at the restart
            new = d_[d_.tail_id.isin(restart) & (d_.minutes_since_t0_trace < CL["start_minutes"])].groupby("tail_id").A_db.mean()
            d_.loc[d_.tail_id.isin(restart), "A_start_trace_db"] = d_.tail_id.map(new)
        d_["event_copy"] = d_.event_id.map(ev_code)
        file = CACHE / f"data_{name}.parquet"
        if not file.exists():
            d_.drop(columns="event_id").to_parquet(file)
        a = d_.A_db[d_.A_db.notna()].to_numpy()
        w_hi = float(np.sort(a)[-FIT["support_min_minutes"]])           # the 200 tail-minutes rule of step 3b-3
        inside = np.log(a[(a >= W_LO) & (a <= w_hi)])
        q = np.quantile(inside, np.linspace(0, 1, S3B["knots"])[1:-1])
        sel = ix[ix.tail_id.isin(ids)]
        meta[name] = {"label": label, "file": str(file), "tails": len(sel), "events": int(sel.event_id.nunique()),
                      "links": int(sel.cml_id.nunique()), "w_hi": w_hi, "points": int(len(a)), "minutes_in_range": int(len(inside)),
                      "knots": [W_LO] + [float(v) for v in np.exp(q)] + [w_hi]}
    fakes = {f"fake_{t}_{r}": (t, r, FAKE_DIR / f"sim_main_{t}_{r}.parquet") for t in PC["fake_truths"] for r in PC["fake_reps"]}
    fakes = {k: v for k, v in fakes.items() if v[2].exists()}
    plan_file = CACHE / "plan.json"
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
        print(f"plan from {rel(plan_file)}: cut: {plan['cuts'] or 'nothing'}")
    else:
        _, D, _ = pack_file(meta["main"]["file"], WIN, FIT["min_points"])
        cv = wc.Curve(meta["main"]["knots"], prior_sd=FIT["prior_sd_logk"])
        cv.fit_noise = True
        x_flat = np.r_[wc.start_theta(FIT["starts"]["flat"], cv.knots), 0.0, 0.0]
        wp.pooled_loglik("F5", [0.03, 0.01], 1.0, D, True, 1.5), cv.loglik(x_flat, D)       # compile
        t0 = time.time()
        for _ in range(5):
            wp.pooled_loglik("F5", [0.03, 0.01], 1.0, D, True, 1.5)
        t_menu = (time.time() - t0) / 5
        t0 = time.time()
        for _ in range(5):
            cv.loglik(x_flat, D)
        t_curve = (time.time() - t0) / 5
        t0 = time.time()
        one = run_menu({"name": "main", "file": meta["main"]["file"], "cand": "F1", "window": WIN, "min_points": FIT["min_points"]})
        t_fit = (time.time() - t0) / len(one)
        print(f"timing on the main set: one likelihood evaluation {t_menu * 1e3:.0f} ms (menu law), {t_curve * 1e3:.0f} ms (free curve); "
              f"one F1 fit with both fixes {t_fit:.1f} s ({np.mean([r['nfev'] for r in one]):.0f} evaluations)", flush=True)
        n_menu = sum(len(menu_starts(c)) for c in S4["candidates"])
        ev_curve = PC["typical_nfev"]["curve_factor"] * (len(FIT["starts"]) * FIT["maxfev"]["start"] + FIT["maxfev"]["finish"]
                                                         + (len(FIT["lambda_grid"]) - 1) * FIT["maxfev"]["other_lambda"] + len(FIT["lambda_grid"]) * FIT["maxfev"]["cv"])
        size = {n: m["points"] / meta["main"]["points"] for n, m in meta.items()}

        def project(names):
            sec = sum(size[n] * (n_menu * PC["typical_nfev"]["menu"] * t_menu + ev_curve * t_curve) for n in names)
            return (sec + len(fakes) * 1.15 * n_menu * PC["typical_nfev"]["menu"] * t_menu) / workers

        plan = {"sets": list(sets), "b4_per_link": True, "wind": True, "cuts": [], "wall_s": 0.0,
                "t_menu_ms": round(t_menu * 1e3, 1), "t_curve_ms": round(t_curve * 1e3, 1), "t_f1_fit_s": round(t_fit, 1)}
        proj = project(plan["sets"])
        plan["projected_full_s"] = round(proj)
        print(f"projection: {proj / 3600:.2f} h on {workers} processes (budget {S['time_budget_h']} h)", flush=True)
        for cut, change in (("the two per-class sets of C3", {"sets": ["main", "no_rise_heavy"]}),
                            ("B4 with one shared stretch", {"b4_per_link": False}), ("the wind columns of part A", {"wind": False})):
            if proj <= S["time_budget_h"] * 3600:
                break
            plan.update(change)
            plan["cuts"].append(cut)
            proj = project(plan["sets"])
            print(f"   cut: {cut} -> projected {proj / 3600:.2f} h", flush=True)
        plan_file.write_text(json.dumps(plan))

    # ------------------------------------------------------------------------------------------- part A
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    T = T.merge(outside_evidence(d, ix), on="tail_id")
    if not plan["wind"]:
        T["wind_speed_ms"] = np.nan
    T["coherence_share"], T["coherence_n"] = coherence(T)
    EV = [("ev_radar_path", "radar R_max > 0 on the path of the link"), ("ev_radar_near", f"radar R_max > 0 on a link path within {PA['near_km']:g} km"),
          ("ev_gauge_city", "nearest City gauge > 0 mm"), ("ev_gauge_smhi", "SMHI gauge > 0 mm"),
          ("ev_any_rain", "ANY of the four rain indicators"), ("ev_dew_risk", "dew risk (hourly product of the baseline module)"),
          ("ev_any_moisture", "ANY moisture evidence (rain indicators or dew risk) = THE RULE"),
          ("ev_night", "night at t0_trace (NOT moisture evidence; information, added after the first look)")]
    T["ev_night"] = T.night
    a_rows = []
    for col, text in EV:
        g = T[T[col].notna()]
        lo, hi = event_boot(g, lambda q, c=col: odds_ratio(q, c), rng, PA["boot_reps"])
        orr = odds_ratio(g, col)
        if not g[col].astype(bool).any():                  # never true: there is nothing to compare
            orr = lo = hi = np.nan
        a_rows.append({"indicator": text, "column": col, "tails_known": len(g),
                       "share_not_decaying": round(float(g[g.not_decaying][col].astype(bool).mean()), 3),
                       "share_decays": round(float(g[~g.not_decaying][col].astype(bool).mean()), 3),
                       **{f"share_{c}": round(float(g[g["class"] == c][col].astype(bool).mean()), 3) for c in CLASSES[1:]},
                       "odds_ratio": round(orr, 2), "or_lo": round(lo, 2), "or_hi": round(hi, 2),
                       "rule_met": bool(orr > PA["or_min"] and lo > 1) if col == "ev_any_moisture" else ""})
    part_a = pd.DataFrame(a_rows)
    save(part_a, P.STALL_PART_A, index=False)
    coh = T[T.coherence_share.notna()]
    coh_ratio = lambda q: q[q.not_decaying].coherence_share.mean() / max(q[~q.not_decaying].coherence_share.mean(), 1e-9)   # noqa: E731
    coh_val, coh_ci = float(coh_ratio(coh)), event_boot(coh, coh_ratio, rng, PA["boot_reps"])
    coh_nd, coh_d = float(coh[coh.not_decaying].coherence_share.mean()), float(coh[~coh.not_decaying].coherence_share.mean())
    pairs = T.groupby(["cml_id", "event_id"]).filter(lambda g: len(g) == 2)
    pg = pairs.groupby(["cml_id", "event_id"])
    pair_same = pg.not_decaying.nunique() == 1
    pair_same4 = pg["class"].nunique() == 1
    p_ev = T.groupby("event_id").not_decaying.mean()
    pair_exp = float(np.mean([p_ev[e] ** 2 + (1 - p_ev[e]) ** 2 for (_, e) in pair_same.index]))
    T["short"] = T.L_km < PA["short_km"]
    rise_all = rho_with_interval(T, "L_km", "rise_db", rng)
    rise_nd = rho_with_interval(T[T["class"] == "rises"], "L_km", "rise_db", rng)
    by_len = T.groupby("short").agg(tails=("tail_id", "size"), median_rise_db=("rise_db", "median"),
                                    share_rises=("class", lambda c: (c == "rises").mean()), share_not_decaying=("not_decaying", "mean"))
    met_tab = T.groupby("class")[["humidity_pct", "temperature_c", "dew_point_depression_c", "wind_speed_ms", "hour_local",
                                  "sun_elevation_deg", "radar_near_wet_share", "city_gauge_km"]].median().reindex(CLASSES).round(2)
    met_tab.insert(0, "night_share", T.groupby("class").night.apply(lambda v: v.dropna().astype(bool).mean()).reindex(CLASSES).round(2))
    met_tab.insert(0, "tails", T["class"].value_counts().reindex(CLASSES))
    h1_or, h1_coh = bool(part_a.set_index("column").rule_met["ev_any_moisture"] is True), bool(coh_val > PA["coherence_ratio"])
    any_row = part_a.set_index("column").loc["ev_any_moisture"]
    print(f"part A: any-moisture odds ratio {any_row.odds_ratio} ({any_row.or_lo}..{any_row.or_hi}); "
          f"coherence ratio {coh_val:.2f}; H1 {'supported' if h1_or or h1_coh else 'not supported'} by the rule", flush=True)

    # ------------------------------------------------------------------------------------------- part B
    rb = np.random.default_rng(PB["seed"])
    pl = T[T["class"] == "plateau-then-decays"]
    iqr = lambda v: float(v.quantile(0.75) - v.quantile(0.25)) if len(v) else np.nan   # noqa: E731
    b1 = {"n": len(pl), "median_db": float(pl.transition_level_db.median()), "iqr_db": iqr(pl.transition_level_db),
          "median_frac": float((pl.transition_level_db / pl.A_start_db).median()), "iqr_frac": iqr(pl.transition_level_db / pl.A_start_db),
          "iqr_a_start_db": iqr(pl.A_start_db), "iqr_a_start_all_db": iqr(T.A_start_db)}
    per_link = pl.groupby("cml_id").transition_level_db.agg(["size", lambda v: v.max() - v.min()])
    per_link = per_link[per_link["size"] >= 2]
    b1["links_2plus"], b1["median_range_in_link_db"] = len(per_link), float(per_link.iloc[:, 1].median()) if len(per_link) else np.nan
    b2_stat, b2_p, b2_links, b2_units, b2_null_sd = link_repeat_test(T, rb)
    b3 = {y: rho_with_interval(pl, y, "min_to_80pct", rb) for y in ("A_start_db", "wet_radar_R_mean_sum_mm", "wet_span_min")}
    T["min_to_80pct_censored"] = T.min_to_80pct.fillna(WIN)
    b3_all = {y: rho_with_interval(T, y, "min_to_80pct_censored", rb) for y in ("A_start_db", "wet_radar_R_mean_sum_mm", "wet_span_min")}
    b3_met = any(b3[y][0] > PB["b3_rho"] and b3[y][1] > 0 for y in ("A_start_db", "wet_radar_R_mean_sum_mm"))

    # B4: master curve
    M, placed, n_base, n_first = build_master(series, T, plan["b4_per_link"])
    own = {tid: M.place(*series[tid]) for tid in T.tail_id[T["class"] == "decays"]}
    b4_rows, plateau_fit = [], {}
    for r in pl.itertuples():
        t, y = series[r.tail_id]
        dec = t >= r.min_to_80pct
        f, d_, rms = M.place(t[dec], y[dec])
        v = {"f": f, "d": d_, "rms": rms, "plateau_min": r.min_to_80pct, "resid": np.nan, "covered": 0.0}
        if np.isfinite(f):
            m = M.value(f * t[~dec] + d_)
            v["covered"] = float(np.isfinite(m).mean()) if len(m) else 0.0
            v["resid"] = float(np.nanmean(y[~dec] - m)) if np.isfinite(m).any() else np.nan
        plateau_fit[r.tail_id] = v
    grp = {"'decays', events not rise-heavy (they built the curve; own shift and stretch)":
           [own[k][2] for k in T.tail_id[(T["class"] == "decays") & ~T.event_rise_heavy]],
           "'decays', rise-heavy events (NOT used to build the curve)": [own[k][2] for k in T.tail_id[(T["class"] == "decays") & T.event_rise_heavy]],
           "'plateau-then-decays', decaying part only": [v["rms"] for v in plateau_fit.values()]}
    for name, v in grp.items():
        v = np.asarray(v, float)
        b4_rows.append({"group": name, "tails": len(v), "placed": int(np.isfinite(v).sum()),
                        "median_rms_db": round(float(np.nanmedian(v)), 3) if np.isfinite(v).any() else np.nan,
                        "p25_rms_db": round(float(np.nanpercentile(v, 25)), 3) if np.isfinite(v).any() else np.nan,
                        "p75_rms_db": round(float(np.nanpercentile(v, 75)), 3) if np.isfinite(v).any() else np.nan})
    ref_rms = b4_rows[0]["median_rms_db"]
    res_ = np.array([v["resid"] for v in plateau_fit.values()], float)
    b4_pl = {"placed": int(np.isfinite([v["f"] for v in plateau_fit.values()]).sum()), "with_master": int(np.isfinite(res_).sum()),
             "on": int((np.abs(res_) <= ref_rms).sum()), "below": int((res_ < -ref_rms).sum()), "above": int((res_ > ref_rms).sum()),
             "median_resid": float(np.nanmedian(res_)) if np.isfinite(res_).any() else np.nan}
    mt = pd.DataFrame({"tau_min": M.tau, "master_A_db": np.round(M.M, 4), "tails": M.n})
    mt = mt[mt.master_A_db.notna()]
    save(pd.concat([pd.DataFrame(b4_rows).assign(table="misfit"), mt.assign(table="master curve")], ignore_index=True), P.STALL_MASTER, index=False)

    # B5: model-free k(A) per subset
    RR = S4["raw_rate"]
    dec_ids = set(T.tail_id[T["class"] == "decays"])
    rh_ids, big_ids = set(T.tail_id[T.event_rise_heavy]), set(T.tail_id[T.A_start_db >= PB["a_start_split_db"]])
    all_ids = set(T.tail_id)
    subsets = {"decaying tails, events not rise-heavy (the rule)": dec_ids - rh_ids, "decaying tails, all events": dec_ids,
               "all tails": all_ids, "events not rise-heavy": all_ids - rh_ids, "rise-heavy events": rh_ids,
               f"not rise-heavy, A_start < {PB['a_start_split_db']} dB": all_ids - rh_ids - big_ids,
               f"not rise-heavy, A_start >= {PB['a_start_split_db']} dB": (all_ids - rh_ids) & big_ids}

    def rate_slope(dd):
        if "draw" in dd:                                   # event bootstrap: a copy of an event is a new set of tails
            dd = dd.assign(tail_id=dd.tail_id * 100 + dd.draw)
        r_ = raw_rate(dd, RR["edges_db"], RR["smooth_min"])
        r_["used_for_slope"] = (r_.tail_minutes >= PB["b5_min_minutes"]) & (r_.k_from_mean_per_min > 0)
        u = r_[r_.used_for_slope]
        return r_, (float(np.log(u.k_from_mean_per_min) @ wc.slope_weights(np.log(u.median_A_db))) if len(u) >= 3 else np.nan)

    rate_rows, b5 = [], {}
    for name, ids in subsets.items():
        dd = first[first.tail_id.isin(ids)]
        r_, sl = rate_slope(dd)
        lo, hi = event_boot(dd, lambda q: rate_slope(q)[1], rb, PB["b5_boot_reps"]) if name.startswith("decaying tails, events not") or name == "all tails" else (np.nan, np.nan)
        b5[name] = (sl, lo, hi, len(ids), int(r_.used_for_slope.sum()))
        rate_rows.append(r_.assign(subset=name, tails_in_subset=len(ids), slope=round(sl, 3)))
    # added after the first look (information, no rule): the same table for the first minutes of a tail and for the
    # rest. A law in which only A decides the rate gives the same k at the same level whenever the level is reached.
    split_rows = []
    for name in ("all tails", "events not rise-heavy", "rise-heavy events"):
        for lab, (m0, m1) in PB["time_split_min"].items():
            dd = first[first.tail_id.isin(subsets[name]) & first.minutes_since_t0_trace.between(m0, m1 - 1)]
            r_, sl = rate_slope(dd)
            rate_rows.append(r_.assign(subset=f"{name}, {lab}", tails_in_subset=len(subsets[name]), slope=round(sl, 3)))
            split_rows.append(r_[r_.used_for_slope].assign(subset=name, minutes=lab))
    split = pd.concat(split_rows).pivot(index="A_bin_db", columns=["subset", "minutes"], values="k_from_mean_per_min")
    split = split.reindex(rate_rows[0].A_bin_db).dropna(how="all")
    split.columns = [f"{a}, {b}" for a, b in split.columns]
    rate = pd.concat(rate_rows, ignore_index=True)
    save(rate, P.STALL_RATE, index=False)
    b5_rule = b5["decaying tails, events not rise-heavy (the rule)"]
    tests = [("B1", "IQR of the transition level of the plateau-then-decays tails (dB)", b1["iqr_db"], f"< {PB['b1_iqr_db']}", b1["iqr_db"] < PB["b1_iqr_db"]),
             ("B2", "the same link in different events: permutation p", b2_p, f"< {PB['b2_p']}", b2_p < PB["b2_p"]),
             ("B3", "Spearman of the plateau duration with A_start / with the rain amount",
              max(b3["A_start_db"][0], b3["wet_radar_R_mean_sum_mm"][0]), f"> {PB['b3_rho']}, interval excluding 0", b3_met),
             ("B5", "slope of log k vs log A, decaying tails of the events not rise-heavy", b5_rule[0], f"< {PB['b5_slope']}", b5_rule[0] < PB["b5_slope"])]
    part_b = pd.DataFrame([{"test": a, "what": b, "value": round(float(c), 3), "rule": e, "met": bool(f)} for a, b, c, e, f in tests])
    n_met = int(part_b.met.sum())
    h2 = n_met >= PB["need"]
    save(part_b, P.STALL_PART_B, index=False)
    print(f"part B: {n_met} of 4 rules met ({', '.join(part_b.test[part_b.met]) or 'none'}); H2 {'supported' if h2 else 'not supported'} by the rule", flush=True)

    # per-tail product (classes + outside evidence + B4 placement)
    T["b4_rms_db"] = T.tail_id.map({**{k: v[2] for k, v in own.items()}, **{k: v["rms"] for k, v in plateau_fit.items()}})
    T["b4_plateau_minus_master_db"] = T.tail_id.map({k: v["resid"] for k, v in plateau_fit.items()})
    save(T.drop(columns=["baseline_db", "short", "min_to_80pct_censored"]).round(4), P.TAIL_CLASSES, index=False)

    # ------------------------------------------------------------------------------------------- part C: fits
    def run(jobs, label):
        t_run, out = time.time(), []
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(run_any, j): j for j in jobs}
            for k, fu in enumerate(as_completed(futs)):
                out.append((futs[fu], fu.result()))
                j = futs[fu]
                print(f"   {label} {k + 1}/{len(jobs)}: {j['part']} {j['name']} {j.get('cand', '')} ({time.time() - t_run:.0f} s)", flush=True)
        return out, time.time() - t_run

    cells = [{"name": n, "file": meta[n]["file"], "window": WIN, "min_points": FIT["min_points"]} for n in plan["sets"]] + [
        {"name": n, "file": str(v[2]), "window": S1C["window_min"], "min_points": S1C["min_points"]} for n, v in fakes.items()]
    jobs1 = [{"part": "menu", "cand": "F1", **c} for c in cells] + [
        {"part": "curve", "name": n, "file": meta[n]["file"], "knots": meta[n]["knots"], "w_hi": meta[n]["w_hi"]} for n in plan["sets"]]
    out1, wall1 = run(jobs1, "phase 1")
    out2, wall2 = run([{"part": "menu", "cand": c, **cell} for cell in cells for c in S4["candidates"] if c != "F1"], "phase 2")
    if wall1 + wall2 > 60:
        plan["wall_s"] = round(plan["wall_s"] + wall1 + wall2, 1)
        plan_file.write_text(json.dumps(plan))
    curve = {r["dataset"]: r for j, r in out1 if j["part"] == "curve"}
    fits = pd.DataFrame([r for j, rs in out1 + out2 if j["part"] == "menu" for r in rs])
    fits["bic"] = -2 * fits.ll + fits.k * np.log(fits.n_independent)
    best = fits.loc[fits.groupby(["dataset", "candidate"]).ll.idxmax()].copy()
    spread = (2 * (fits.groupby(["dataset", "candidate"]).ll.max() - fits.groupby(["dataset", "candidate"]).ll.min())).round(1)
    best = best.merge(spread.rename("bic_spread_starts").reset_index(), on=["dataset", "candidate"])
    best["dbic"] = (best.bic - best.groupby("dataset").bic.transform("min")).round(2)
    best["parameters"] = [", ".join(f"{n_} {v:.4g}" for n_, v in zip(wp.NAMES[c], th)) for c, th in zip(best.candidate, best.theta)]

    # C2 check on the fake datasets: before (cached fits of step 3b-3c part C, old code) vs after (both fixes)
    f_rows = []
    for name, (truth, rep, _) in fakes.items():
        old = []
        for c in S4["candidates"]:
            for lb in menu_starts(c):
                f_ = CACHE_BEFORE / f"fit_c_{truth}_{rep}_{c}_{lb}.json"
                if f_.exists():
                    old.append(json.loads(f_.read_text()))
        if not old:
            continue
        old = pd.DataFrame(old)
        old["bic"] = -2 * old.ll + old.k * np.log(old.n_independent)
        ob = old.loc[old.groupby("candidate").ll.idxmax()].set_index("candidate")
        nb = best[best.dataset == name].set_index("candidate")
        own_c = PC["own"][truth]
        ratio = np.array(nb.theta[own_c]) / np.array(ob.theta[own_c])
        mult = float(nb.noise_mult[nb.bic.idxmin()])
        ok = bool(ob.bic.idxmin() == nb.bic.idxmin() and np.all(np.abs(ratio - 1) <= PC["tol"])
                  and PC["mult_range"][0] <= mult <= PC["mult_range"][1])
        f_rows.append({"truth": truth, "replicate": rep, "winner_before": ob.bic.idxmin(), "winner_after": nb.bic.idxmin(),
                       "margin_before_bic": round(float(ob.bic.sort_values().iloc[1] - ob.bic.min()), 1),
                       "margin_after_bic": round(float(nb.bic.sort_values().iloc[1] - nb.bic.min()), 1),
                       "own_candidate": own_c, "parameters_before": ", ".join(f"{n_} {v:.4g}" for n_, v in zip(wp.NAMES[own_c], ob.theta[own_c])),
                       "parameters_after": ", ".join(f"{n_} {v:.4g}" for n_, v in zip(wp.NAMES[own_c], nb.theta[own_c])),
                       "largest_change_pct": round(float(np.max(np.abs(ratio - 1)) * 100), 1),
                       "s_before": round(float(ob.s[own_c]), 3), "s_after": round(float(nb.s[own_c]), 3),
                       "noise_multiplier": round(mult, 3),
                       "multiplier_all_candidates": f"{nb.noise_mult.min():.2f}..{nb.noise_mult.max():.2f}", "passed": ok})
    fake = pd.DataFrame(f_rows)
    save(fake, P.STALL_FAKE, index=False)
    fake_ok = bool(len(fake) == len(PC["fake_truths"]) * len(PC["fake_reps"]) and fake.passed.all())
    print(f"part C fake check: {int(fake.passed.sum()) if len(fake) else 0} of {len(fake)} datasets pass", flush=True)

    # C3: the real sets with both fixes
    real = best[best.dataset.isin(plan["sets"])].copy()
    real["family"] = real.candidate.map(FAMILY)
    menu = real[["dataset", "candidate", "family", "n_tails", "n_independent", "k", "ll", "bic", "dbic", "parameters", "s", "noise_mult",
                 "start", "bic_spread_starts"]].sort_values(["dataset", "bic"]).round({"ll": 2, "bic": 2, "s": 3, "noise_mult": 3})
    menu["dataset"] = pd.Categorical(menu.dataset, plan["sets"])
    menu = menu.sort_values(["dataset", "bic"])
    save(menu, P.STALL_MENU, index=False)
    old4 = pd.read_csv(P.REAL_SUMMARY).set_index("dataset") if P.REAL_SUMMARY.exists() else pd.DataFrame()
    c_rows, slopes, kcurves = [], {}, {}
    for name in plan["sets"]:
        g = real[real.dataset == name].set_index("candidate")
        rd, c_ = read_menu(g.reset_index()[["candidate", "bic"]]), curve[name]
        w = np.geomspace(W_LO, c_["w_hi"], FIT["grid_points"])
        k = wc.Curve(c_["knots"]).k(np.asarray(c_["x"])[:-2], w)
        slopes[name] = float(np.log(k) @ wc.slope_weights(np.log(w)))
        kcurves[name] = (w, k)
        row = {"dataset": name, "what": meta[name]["label"], "tails": c_["n_tails"], "events": meta[name]["events"],
               "w_hi_db": round(meta[name]["w_hi"], 2), "menu_winner": rd["winner"], "menu_family": rd["family"],
               "runner_up": rd["runner_up"], "margin_bic": round(rd["margin"], 1), "next_family": rd["next_family"],
               "margin_family_bic": round(rd["margin_family"], 1), "f1_minus_winner_bic": round(rd["f1_minus_winner"], 1),
               "f1_tau_min": round(g.theta["F1"][0], 1), "f3_n": round(g.theta["F3"][1], 3), "f5_b": float(f"{g.theta['F5'][1]:.3g}"),
               "noise_multiplier_winner": round(float(g.noise_mult[rd["winner"]]), 2),
               "noise_multiplier_range": f"{g.noise_mult.min():.2f}..{g.noise_mult.max():.2f}",
               "s_winner": round(float(g.s[rd["winner"]]), 3), "curve_slope": round(slopes[name], 3), "curve_reading": read_slope(slopes[name]),
               "curve_noise_multiplier": round(c_["noise_mult"], 2), "curve_s": round(c_["s"], 3), "curve_penalty": c_["lam_star"],
               "curve_minus_winner_bic": round(c_["bic"] - g.bic[rd["winner"]], 1),
               "curve_start_cost_spread": round(max(v["cost"] for v in c_["starts"].values()) - min(v["cost"] for v in c_["starts"].values()), 1),
               "k_at_knots_per_min": " ".join(f"{np.exp(v):.4f}" for v in c_["x"][:-2]),
               "knots_db": " ".join(f"{v:.3g}" for v in c_["knots"])}
        if name in old4.index:
            o = old4.loc[name]
            row.update({"step4_winner_as_run": o.menu_winner, "step4_like_for_like_winner": o.like_for_like_winner,
                        "step4_like_for_like_margin_bic": o.like_for_like_margin_family_bic, "step4_curve_slope": o.curve_slope,
                        "step4_f1_tau_min": o.f1_tau_min})
        c_rows.append(row)
    summ = pd.DataFrame(c_rows)
    save(summ, P.STALL_SUMMARY, index=False)

    # ------------------------------------------------------------------------------------------- figures
    rg = np.random.default_rng(S["gallery"]["seed"])
    picks = {c: sorted(rg.choice(T.tail_id[T["class"] == c], min(S["gallery"]["per_class"], int((T["class"] == c).sum())), replace=False).tolist())
             for c in CLASSES}
    fig = gallery_figure(d, T, dict(zip(ix.tail_id, A)), picks)
    save(fig, P.STALL_FIG_GALLERY, dpi=105)
    plt.close(fig)
    fig = master_figure(M, series, T, placed, plateau_fit)
    save(fig, P.STALL_FIG_MASTER, dpi=115)
    plt.close(fig)
    fig = rate_figure(rate)
    save(fig, P.STALL_FIG_RATE, dpi=115)
    plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    C0 = count(T)
    pct = lambda v: f"{100 * v:.0f} %"                       # noqa: E731
    ci = lambda v: f"{v[0]:+.2f} ({v[1]:+.2f}..{v[2]:+.2f}, n = {v[3]})" if np.isfinite(v[0]) else f"not computed (n = {v[3]})"   # noqa: E731
    ev_tab = counts[~counts.group.isin(["all tails", "events flagged rise-heavy", "other events"])].sort_values("tails", ascending=False)
    b3_tab = pd.DataFrame([{"plateau duration against": lab, "plateau-then-decays tails only (the rule)": ci(b3[y]),
                            "all tails, minutes until A < 0.8 A_start (60 if never; information)": ci(b3_all[y])}
                           for y, lab in (("A_start_db", "A_start (dB)"), ("wet_radar_R_mean_sum_mm", "rain amount of the event on the path (radar, mm)"),
                                          ("wet_span_min", "rain duration on the link (min)"))])
    b5_tab = pd.DataFrame([{"subset": n, "tails": v[3], "bins used": v[4], "slope of log k vs log A": round(v[0], 2) if np.isfinite(v[0]) else np.nan,
                            "event-bootstrap 90 %": f"{v[1]:+.2f}..{v[2]:+.2f}" if np.isfinite(v[1]) else ""} for n, v in b5.items()])
    rate_wide = rate[rate.used_for_slope].pivot(index="A_bin_db", columns="subset", values="k_from_mean_per_min").reindex(
        index=rate.A_bin_db.drop_duplicates(), columns=list(subsets)).dropna(how="all").round(4).reset_index()
    cuts = plan["cuts"] or ["nothing"]
    nan_radar = int((T.radar_path_nan_bins > 0).sum())
    wall = plan["wall_s"]
    report = f"""# Step 4b: why do some real tails stall or rise? Contamination or antenna physics?

Script: `{rel(P.MODULE / "diagnose_stalled_tails.py")}`. Settings and every rule: `research/config.yaml`, section
`wet_antenna.diagnose_stalled_tails` (written before the run; nothing was re-run with changed rules).
Data: the {len(T)} TRAIN tails of the frozen set v1 ({T.event_id.nunique()} events, {T.cml_id.nunique()} links), first {WIN} min after
`t0_trace`. The 123 VAL tails were never read; no file of an earlier step in `outputs/products/` was changed.

**The question.** Step 4 found that about one tail in five does not dry in its first minutes and that the drying rate
per dB falls with the level. Either something still wets the antenna after the rain reference says "dry" (H1,
contamination: it belongs to the weather), or the antenna itself dries more slowly per dB when there is more water on
it (H2, antenna physics: it belongs to the drying law). The two explanations predict different things, and this step
looks for those differences.

## 1. What was run, what was cut, run time

- Classes of all {len(T)} tails; part A (outside evidence); part B (B1 to B5); part C (two code fixes, the check on
  {len(fake)} fake datasets, the menu fits and the free curve again on {len(plan['sets'])} real sets).
- Timing first: one likelihood evaluation {plan['t_menu_ms']:.0f} ms (menu law) and {plan['t_curve_ms']:.0f} ms (free curve);
  one F1 fit with both fixes {plan['t_f1_fit_s']} s. Projection {plan['projected_full_s'] / 3600:.2f} h on {workers} processes, budget
  {S['time_budget_h']} h. **Cut: {'; '.join(cuts)}.**
- Fits: {len(fits)} menu fits and {len(curve)} free-curve fits, {wall / 60:.0f} min of wall time for the fits; everything else
  takes about a minute. Cached in `{rel(CACHE)}/`.
- The spec named files and tools of the old layout (`roadmap.md`, `DECISIONS.md`, `save_output` with importance
  levels, `check_organization.py`, `.cache/fits4b/`). The repository was reorganised since: this step uses
  `docs/ROADMAP.md`, `METHOD.md`, `save`, the cache folder above, and the checks that exist now. There are no
  importance levels any more, so no output is "key", "product" or "diag"; the per-tail table is in `outputs/products/`
  because later steps will read it.

## 2. The classes

Rules (fixed first). A is smoothed with a centred {CL['smooth_min']} min median; A_start = mean of minutes 0 to {CL['start_minutes'] - 1}.
The classes are tried in this order, and a tail gets the first one it meets:

1. **decays**: A falls below {CL['half']} x A_start within {CL['half_within_min']} min and never climbs more than {CL['rise_db']} dB above its lowest level so far.
2. **plateau-then-decays**: A stays at or above {CL['plateau_level']} x A_start for at least {CL['plateau_min']} min, then falls below {CL['half']} x A_start before {WIN} min.
3. **rises**: A exceeds A_start + {CL['rise_db']} dB at some time.
4. **flat / slow**: everything else.

The order is my reading of the spec (it did not say what happens when two rules hold). Tails meeting two rules:
{', '.join(f'{k}: {v}' for k, v in overlap.items())}. The columns `meets_decays`, `meets_plateau`, `meets_rises` of the per-tail table
keep all three answers.

{md_table(counts[counts.group.isin(['all tails', 'events flagged rise-heavy', 'other events'])].drop(columns='rise_heavy'))}

Per event (largest first):

{md_table(ev_tab)}

## 3. Part A: outside evidence of continued wetting (H1)

Per tail, from `t0_trace` to `t0_trace` + {WIN} min. An indicator is "true" when the reference shows any water at all
(no threshold). Odds ratio = odds in the non-decaying tails (classes 2 to 4, {int(T.not_decaying.sum())} tails) over the odds in the
decaying tails ({int((~T.not_decaying).sum())}); the interval is the 90 % event bootstrap ({PA['boot_reps']} draws of whole events).

{md_table(part_a.drop(columns='column'))}

**Read this table with one fact in mind.** The frozen tail set v1 starts every tail at the LAST rain trace
(`t0_trace`) and cuts it at the next one, where "trace" = radar R_max > 0 on the path, a radar gap, or a gauge within
1 km > 0. All {len(T)} tails are longer than {WIN} min, so inside the window the radar on the path and the near gauges are
dry BY CONSTRUCTION (radar on the path > 0 in {int((T.ev_radar_path == True).sum())} tails; {nan_radar} tails with a missing radar bin).
The only rain references that CAN show something are the ones farther away: other link paths within {PA['near_km']:g} km,
the nearest City gauge when it is more than 1 km away (median distance {T.city_gauge_km.median():.1f} km), and the SMHI gauge
(median {T.smhi_gauge_km.median():.1f} km). Humidity exists for {pct(T.humidity_pct.notna().mean())} of the tails: it is missing in the
hourly SMHI record for the events {', '.join(sorted(T[T.humidity_pct.isna()].event_id.unique()))} ({int(T.humidity_pct.isna().sum())} tails, of which
{int((T.humidity_pct.isna() & T.not_decaying).sum())} do not decay and {int((T.humidity_pct.isna() & (T['class'] == 'rises')).sum())} of the {int((T['class'] == 'rises').sum())} "rises" tails). There the dew indicator is
"unknown" and those tails are left out of its row. The weather is ONE station for the whole city, hourly.

Weather and time of day per class (medians):

{md_table(met_tab, index=True)}

Other parts of part A:

- **Same class in the event's other tails at the same time** (other links, start within {PA['overlap_min']} min): the share of
  "rises" or "flat / slow" neighbours is {coh_nd:.2f} for the non-decaying tails and {coh_d:.2f} for the decaying ones,
  ratio **{coh_val:.2f}** (event bootstrap {coh_ci[0]:.2f}..{coh_ci[1]:.2f}; rule: > {PA['coherence_ratio']}). {int(T.coherence_share.isna().sum())} tails
  have no such neighbour.
- **The two directions of one link** ({len(pair_same)} link-events with both sub-links): same answer to "decays or not" in
  {pct(pair_same.mean())} (expected {pct(pair_exp)} if the tails of an event were independent); same class of four in {pct(pair_same4.mean())}.
  The two directions share both antennas and the path, so this shows the classes are not noise; it cannot separate H1 from H2.
- **Does the rise grow with the path length?** Rise = highest smoothed A minus A_start. Spearman with L: all tails
  {ci(rise_all)}; tails of class "rises" {ci(rise_nd)}.

{md_table(by_len.rename(index={True: f'short, L < {PA["short_km"]:g} km', False: f'long, L >= {PA["short_km"]:g} km'}).rename_axis("links").round(3), index=True)}

**H1 verdict by the rule: {'SUPPORTED' if h1_or or h1_coh else 'NOT SUPPORTED'}.** Odds-ratio part ("any moisture evidence" > {PA['or_min']} with the
interval excluding 1): {'met' if h1_or else 'not met'}. Coherence part (ratio > {PA['coherence_ratio']}): {'met' if h1_coh else 'not met'}.

## 4. Part B: tests of a law of the antenna (H2)

**B1, transition level.** {b1['n']} plateau-then-decays tails. Transition level: median {b1['median_db']:.2f} dB, IQR
**{b1['iqr_db']:.2f} dB**; as a fraction of A_start: median {b1['median_frac']:.2f}, IQR {b1['iqr_frac']:.2f}. For comparison the IQR of
A_start itself is {b1['iqr_a_start_db']:.2f} dB in these tails and {b1['iqr_a_start_all_db']:.2f} dB in all tails. Links with two or more such tails:
{b1['links_2plus']}{f", median spread inside a link {b1['median_range_in_link_db']:.2f} dB" if b1['links_2plus'] else ''}.
Warning: by the definition in the rules the transition level is the first A below {CL['plateau_level']} x A_start, so it is
about {CL['plateau_level']} x A_start for every tail. Its spread is therefore the spread of A_start, and the fraction says nothing.

**B2, the same link in different events.** {b2_links} links with tails in at least {PB['min_events_per_link']} events ({b2_units} link-events).
Statistic (mean product of "does not decay", event mean removed, over pairs of events of a link): {b2_stat:+.4f}; the
{PB['shuffles']} shuffles inside events give 0 +- {b2_null_sd:.4f}; **p = {b2_p:.3f}**.

**B3, plateau duration against wetness** (Spearman, 90 % event bootstrap):

{md_table(b3_tab)}

**B4, master curve.** Built from the {n_base} "decays" tails of the events not flagged rise-heavy ({n_first} of them cross
{MC['ref_db']} dB and gave the first alignment; {len(placed)} are in the final curve), {'one stretch per link' if plan['b4_per_link'] else 'ONE shared stretch (cut)'},
{MC['iterations']} rounds. The curve is defined from {mt.master_A_db.max():.2f} dB down to {mt.master_A_db.min():.2f} dB over {len(mt)} min.

{md_table(pd.DataFrame(b4_rows))}

The plateau part of the plateau-then-decays tails, with the shift and stretch found on their decaying part:
{b4_pl['placed']} of {len(plateau_fit)} tails could be placed, {b4_pl['with_master']} have the master curve defined under their plateau.
Mean (A - master) over the plateau: median {b4_pl['median_resid']:+.2f} dB; **on** the curve (within the median misfit
{ref_rms} dB of the first row) {b4_pl['on']}, **below** it {b4_pl['below']}, **above** it {b4_pl['above']}. "Below" means: going back in time
the master curve is higher than the tail was, so the tail sat at its level instead of coming down from above; that
is a delayed start, not a slow part of one common curve.

**B5, the model-free drying rate per level.** k = -mean dA/dt / median A per bin of A (A smoothed over
{RR['smooth_min']} min, central difference), bins with at least {PB['b5_min_minutes']} tail-minutes and k > 0:

{md_table(b5_tab)}

k per bin (1 / min):

{md_table(rate_wide)}

Added after the first look (information, no rule): the same k, split by WHEN in the tail the level is reached. If
only A decided the rate, the two columns of a subset would agree.

{md_table(split.round(4), index=True)}

**The four rules:**

{md_table(part_b)}

**H2 verdict by the rule: {'SUPPORTED' if h2 else 'NOT SUPPORTED'}** ({n_met} of 4 met, {PB['need']} needed). B4 has no rule; it is reported above.

## 5. Part C: the two code fixes, the fake check, the fits again

**C1.** `pooled_fit.py`, option `bound_w0`: the start W0 of F1 and F4 is kept inside [0.001, 30] dB, as it always was
for F3, F5, F6 and the free curve. **C2.** option `noise_mult` / `fit(..., fit_noise=True)`: one fitted number m
multiplies the noise variance of every tail of the dataset; the likelihood gets the term -1/2 n log m. Both options
are OFF by default, so steps 3b-1b to 4 give exactly the numbers they gave before (checked: step 4 re-run from its
cache, `python -m research compare`).

**Check on fake data** ({len(fake)} cached datasets of step 3b-3b; "before" = the cached menu fits of step 3b-3c part C).
Rule: same winner, the parameters of the truth's own candidate within +-{int(PC['tol'] * 100)} %, the multiplier inside
{PC['mult_range'][0]}..{PC['mult_range'][1]}. **{'PASSED' if fake_ok else 'NOT PASSED'}: {int(fake.passed.sum()) if len(fake) else 0} of {len(fake)}.**

{md_table(fake)}

**The real sets with both fixes** (first {WIN} min; every candidate now has the same start bound, so the BIC column IS
like for like; k of the BIC includes the multiplier):

{md_table(summ[['dataset', 'what', 'tails', 'events', 'w_hi_db', 'menu_winner', 'menu_family', 'runner_up', 'margin_bic', 'next_family', 'margin_family_bic', 'f1_minus_winner_bic', 'noise_multiplier_winner']])}

{md_table(summ[['dataset', 'f1_tau_min', 'f3_n', 'f5_b', 's_winner', 'curve_slope', 'curve_reading', 'curve_noise_multiplier', 'curve_minus_winner_bic', 'curve_start_cost_spread', 'k_at_knots_per_min', 'knots_db']])}

Against step 4 (old code):

{md_table(summ[[c for c in ['dataset', 'step4_winner_as_run', 'step4_like_for_like_winner', 'step4_like_for_like_margin_bic', 'menu_winner', 'margin_family_bic', 'step4_curve_slope', 'curve_slope', 'step4_f1_tau_min', 'f1_tau_min'] if c in summ.columns]].dropna())}

All menu fits:

{md_table(menu.drop(columns=['n_independent', 'start']))}

Reading rules as in step 4: |slope| < {S4['reading']['flat_abs_slope']} exponential; > {S4['reading']['quadratic_slope']} quadratic family; in between mildly
rising; negative = falling (no fake truth). F3 and F5 are one family; F4 against F6 by less than {S4['f4_f6_bic']:g} BIC is not decisive.
The two per-class sets are chosen by the SHAPE of the tails, not by a fit, but that shape is the thing the law
describes: "decays" keeps only tails that halve within {CL['half_within_min']} min. Their fits describe those tails and are biased
towards fast drying; they are not the drying law of all tails.

{PLAIN}

## Files

- Per-tail table (classes, outside evidence, B4 placement): `{rel(P.TAIL_CLASSES)}`
- Tables: `{rel(P.STALL_COUNTS)}`, `{rel(P.STALL_PART_A)}`, `{rel(P.STALL_PART_B)}`, `{rel(P.STALL_MASTER)}`,
  `{rel(P.STALL_RATE)}`, `{rel(P.STALL_FAKE)}`, `{rel(P.STALL_MENU)}`, `{rel(P.STALL_SUMMARY)}`
- Figures: `{rel(P.STALL_FIG_GALLERY)}` ({S['gallery']['per_class']} tails per class with radar and gauge under A),
  `{rel(P.STALL_FIG_MASTER)}` (B4), `{rel(P.STALL_FIG_RATE)}` (B5)
"""
    save(report, P.STALL_REPORT)
    print(f"report: {rel(P.STALL_REPORT)} ({time.time() - clock:.0f} s)")


# Sections 6 to 9 of the report: written by hand AFTER the run, from the numbers above (they quote this run).
PLAIN = """## 6. In plain words: which explanation the data support

**Neither H1 nor H2 passes its rule. What the data do show is a third picture: a tail stalls or rises because of
WHEN it happens (the event), not because of WHICH link it is, how long the link is, or how much water is on it.**

- **H2 as "the rate per dB falls with the level": not supported, and fairly sure.** In the tails that decay, the rate
  does not fall with the level: slope +0.28 (event bootstrap +0.05..+0.77) in the events not flagged rise-heavy, +0.10
  in all events. The falling curve of step 4 (-0.80 here, model-free) appears only when tails that wait are mixed in
  (rise-heavy events: -1.69). The same link does not repeat its class in another event (B2, p = 0.35), and the plateau
  is not longer when the tail starts higher or the event had more rain (B3, at most +0.16). After their plateau, the
  plateau tails follow the same curve as the decaying tails just as well (rms 0.117 dB against 0.124 dB), and their
  plateau lies BELOW that curve in 45 of 61 tails: the tail waited, then dried like the others. A wait is a delay, not
  a slow part of one common curve. Limit of this test: "decays" is defined by a fast fall, so a tail that is slow
  because it is high would not be in that class at all (see 9).
- **H1 as "the references show water after t0_trace": not supported, but the test is weak.** The tails that do not
  decay have LESS rain nearby than the decaying ones (any rain indicator: 78 % against 93 %, odds ratio 0.29, interval
  0.10..0.58), not more. But the radar on the path is dry in all 429 tails by the way the tail set was built, and the
  humidity is missing for the two events that hold 60 of the 189 non-decaying tails and 22 of the 36 rising tails.
  So drizzle seen by the references is unlikely; dew or condensation could not be tested where it matters most.
- **What is supported: the cause is at the scale of the event.** The share of non-decaying tails runs from 0 % to
  100 % between events (19 % and 24 % in the two biggest clean events, 83 % in E20150617_0705, 100 % in
  E20150701_0030). The other links of the same event at the same time show the same behaviour more often (ratio 1.85,
  interval 1.27..2.43; the rule asked for 2). The two directions of a link agree in 81 % (64 % expected by chance).
  The rise does not grow with the path length (Spearman 0.00), and short and long links stall equally often (44 %),
  so the extra attenuation is not rain along the path. Rising tails are night tails half of the time (decaying: 8 %).
- **For which tails.** 240 of 429 decay at once (56 %). 64 wait 10 min or more, then decay normally (15 %). 36 rise by
  more than 0.3 dB (8 %; 22 of them in two night events without humidity data). 89 are flat or slow (21 %).
- **How sure.** Sure: the stall is not a property of the level, the link or the path length. Not sure: what it is. A
  weather cause fits (something still wets the antennas without rain on the references, or nothing dries them: night,
  cool air), but with one weather station, hourly values and 24 events it is a description, not a proven cause.

## 7. What this means for the drying law

- **h(W) does not need a falling k.** Do not build "slower when wetter" into the drying law to explain the stalled tails.
- **Which tails should define h(W): the decaying phase.** Two candidate sets are ready in
  `outputs/products/tails_model_classes.csv`: class "decays" (240 tails), and "decays" plus the plateau tails from the end
  of their plateau (304 tails). On both the free curve above 0.35 dB is flat (slopes +0.01 and -0.15, "exponential in
  this range") and F1 gives tau = 10.2 and 10.8 min; on all tails tau is 19.3 min. Warning: these sets keep the tails
  that fall fast, so tau is biased short, and the fits were made with a noise option whose fake-data check did not
  pass (section 5).
- **The stalled and rising minutes need something that is not in h(W).** For the 36 rising tails that is an INPUT
  (the wetting side g), but the rain references show nothing there, so rain R cannot be that input; humidity, dew or
  time of day would have to be. For the tails that only wait, a start time of drying per tail (or a rate that depends
  on the weather) describes them as well as an input does; the data here cannot tell these apart. That question is
  Step 6 (what drying depends on), which should now include night / day, humidity and wind from the start.
- **A new doubt about the "decays" class.** These tails follow rain more closely (rain nearby in 93 %, the gauge wet
  in 41 % against 20 %), and several fall steeply in the first minutes (gallery). Part of that fall may be the last
  path rain fading, not the antenna drying (step 2 already could not prove the fast drop is WAA). If so, the short
  tau of 10 min is too short for the antenna.
- **The menu with both fixes says the same thing in all four sets**: F6 first, F4 second (11.6 to 23 BIC), F1 third (17
  to 37 BIC behind F6); F3 at n = 1 and F5 at b = 0. F6 places its switch at 0.09 to 0.12 dB and F4 gives 3 to 4 % of
  the water a very slow layer: both add a small floor far below 0.35 dB, the range where the fake tests found a false
  feature for every truth. On one of the two fake exponential datasets F6 beat F1 by about 10 BIC with the old code. So read it as
  "exponential, plus something near 0.1 dB that the method cannot name", not as a state switch.
- **The noise during drying is 2.6 to 4.3 times the dry-weather noise variance**, also in the cleanest tails (2.6 in
  "decays"). The dry-weather noise model is too optimistic for drying tails whatever set is chosen.

## 8. Decisions for the learners (not taken)

1. Which tails define h(W): all of v1; "decays" (240); "decays" + plateau tails from the plateau's end (304); or all
   tails with the start of drying as an unknown per tail. A new frozen set is needed for anything but the first.
2. Whether the model gets a per-tail start delay, a weather-dependent rate (Step 6), or an input term without rain.
3. C1 (start bound for F1 / F4): make it the default? It removes the fault of step 4; earlier fake-data steps would
   then have to be re-run to keep their numbers honest (their tails never rise, so little should change).
4. C2 (fitted noise multiplier): the fake check did NOT pass its rule (4 of 6). Accept the two failures (see 9), fix the
   rule, or keep the multiplier out?
5. Get humidity for E20150617_0705 and E20150701_0030 from another station (they hold most of the rising tails).
   Without it, dew cannot be tested.
6. Test whether the fast fall of the "decays" tails is the end of the path rain (compare by path length with enough
   power, or use the two directions and neighbours), before the 10 min tau is called a drying time.
7. Importance levels no longer exist in this repository, so there is no key figure to promote or demote. The figure
   that carries the result is `diagnose_stalled_tails_b5_rate_per_level.png`.

## 9. Unclear, or possibly wrong

- **The class order is mine.** The spec did not say which class wins when two rules hold; I tried them in the listed
  order. The overlap counts are in section 2.
- **B1 cannot work as defined**: the transition level is about 0.8 x A_start by construction, so its IQR (0.35 dB)
  measures the spread of A_start (0.52 dB in these tails). It fails its rule, but it would be uninformative either way.
- **B5 and B4 are biased against H2 by the classes.** "decays" requires half of A_start within 30 min. If the truth
  were "slow when high", high tails would be classed as non-decaying and removed from the very table that tests it.
  The split by A_start inside the clean events (-0.12 for A_start >= 1.5 dB, -0.54 below) does not settle it.
- **The early / late table (added after the first look) is not a clean test either.** k at the same level is 3 to 7
  times higher in minutes 0-14 than later, but a tail that is already low after a few minutes is a fast tail, so the
  spread of speeds between tails gives this pattern too.
- **The dew indicator was corrected during the work**: the hourly product says "no dew risk" where humidity is
  missing; the first version counted that as "no". Now it is "unknown" (69 tails left out). The rule and its limits
  were not changed. "Radar nearby within 5 km" is my addition to the indicators the spec listed (fixed before the run).
- **"Any moisture evidence" is true for 98 % of both groups**, so that rule could hardly have been met: with dew risk
  and far gauges included, almost every hour after rain counts as moist.
- **Coherence 1.85 with interval 1.27..2.43**: below the rule's 2, clearly above 1. I report "not met" as the rule
  says; a reader may fairly call it partial support for a shared, weather-scale cause.
- **Part C, fake check failed by the rule (4 of 6).** m1 replicate 1: the winner changed from F6 to F5, both wrong
  and both only about 10 BIC or less ahead of F1; the tau of F1 moved 1 %. m6 replicate 0: same winner, but k_evap moved
  from 0.0031 to 0.016 / min, while two starts with k_evap 0.0023 and 0.016 end 0.2 log-likelihood apart: the first
  hour does not determine k_evap (known since 3b-1c). The multiplier is 1.03 to 1.08 in all six, inside the range but
  always above 1. My reading: the fixes did not break the tool, but the rule says "not passed", so the C3 numbers are
  not validated.
- **The free curve is still unstable on real tails**: its starts end 55 to 355 apart in cost (fake data: about 8), the
  cross-validation picks a different penalty per set, and the slope without the rise-heavy events moved from -0.05
  (step 4) to -0.51 here with the noise multiplier. Only the sign pattern is worth reading: falling when stalled
  tails are in, flat when they are out. The curve beats every menu law by 1500 to 2300 BIC, and it does so below
  0.35 dB, where it drops to k = 0.002 to 0.010 / min at 0.05 to 0.1 dB; that range is not read.
- **The master curve ignores the per-tail baseline offset** and uses the 5 min median of A, so 0.1 dB of misfit is
  about the noise; "on the curve" was set to the median misfit (0.124 dB), my choice.
- **One weather station, hourly**, for a city 30 km wide: humidity, wind and "night" are the same for every link of
  an event, so they cannot explain differences inside an event.
- **Earlier steps are unchanged**: both options are off by default; step 4 re-run from its cache gives identical tables.
"""


if __name__ == "__main__":
    main()
