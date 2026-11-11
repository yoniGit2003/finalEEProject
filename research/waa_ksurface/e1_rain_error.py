"""Step E1 of the k-surface branch: how large is the radar's own rain error? (rain data only, TRAIN events)

    python -m research.waa_ksurface.e1_rain_error

The length test sorts 5 min values into rain classes by the radar rain. If the radar rain is noisy, the classes are mixed
and the measured F(R) comes out too flat (step A3d). How much depends on the size of the radar's OWN error, which is not the
radar-vs-gauge scatter: a gauge is a point, and it resolves a 5 min rate only in steps of 1.2 mm/h (weighing, 0.1 mm) or
2.4 mm/h (tipping bucket, 0.2 mm). This script splits the scatter without calling any sensor the truth:
(a) triple collocation (radar + two gauges at one place): the error variance of each of the three;
(b) gauge pairs against their separation: the gauge's own error and the point-vs-pixel part;
(c) the binning reversal (gauge binned by radar, radar binned by gauge) and the split it implies;
(d) the result: [sigma_lo, sigma_hi] of the per-bin radar error at 5 min, the per-event scale part, the shared-in-space part.
No attenuation data is read.
"""
import itertools

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import data
from research.waa_ksurface import paths as P

C = data.branch_config()["e1"]
THR = C["log_min_rain_mmh"]                                    # a log needs rain in every sensor: all >= this
EDGES = np.array(C["class_edges_mmh"], float)


# =================================================================================================== data
def load():
    """The rain series inside the TRAIN events (inherited data layer): radar pixel at the gauges and the gauges, at 5 and 15 min.
    Returns dict with r5, g5 (time x City gauge), r15, g15, s15 (SMHI), ev5 / ev15 (event number per bin, -1 outside TRAIN events)."""
    d = data.open_release()
    ev = d.events()
    tr = ev[ev.split_event == "train"].reset_index(drop=True)
    rad = d.store("radar_at_gauges_5min").R_pixel.load().to_pandas()
    c5, c15, s15 = d.store("gauges_city_5min"), d.store("gauges_city_15min"), d.store("gauges_smhi_15min")
    g5 = c5.rate_mm_h.load().to_pandas().where(c5.n_valid.load().to_pandas() == 5)
    g15 = c15.rate_mm_h.load().to_pandas().where(c15.n_valid.load().to_pandas() == 15)
    sm = s15.rate_mm_h.load().to_pandas().where(s15.n_valid.load().to_pandas() >= 1)["GbgA"]
    r15 = rad.resample("15min").mean().where(rad.notna().resample("15min").sum() == 3)
    assert (rad.index == g5.index).all() and (r15.index == g15.index).all() and (sm.index == g15.index).all()

    def event_number(index, width):
        out = np.full(len(index), -1)
        for i, e in enumerate(tr.itertuples()):
            out[(index >= e.start) & (index + pd.Timedelta(width) <= e.end)] = i
        return out
    gauges = d.table("gauges").set_index("gauge_id")
    return {"r5": rad, "g5": g5, "r15": r15, "g15": g15, "s15": sm, "ev5": event_number(rad.index, "5min"), "ev15": event_number(r15.index, "15min"),
            "events": tr, "gauges": gauges, "g_table": data.g_table()}


def centre(x, ev):
    """Remove the mean of every event (per column). Returns the centred values and the degrees-of-freedom factor n / (n - events)."""
    x = np.asarray(x, float).copy()
    ids = np.unique(ev)
    for e in ids:
        m = ev == e
        x[m] -= x[m].mean(axis=0)
    n = len(x)
    return x, n / max(n - len(ids), 1)


def boot_weights(ev, reps, seed):
    """Event bootstrap as a weight per row (how many times its event was drawn)."""
    ids, inv = np.unique(ev, return_inverse=True)
    rng = np.random.default_rng(seed)
    return [np.bincount(rng.integers(0, len(ids), len(ids)), minlength=len(ids))[inv].astype(float) for _ in range(reps)]


def wcov(X, w):
    """Weighted covariance matrix of the columns of X."""
    w = w / w.sum()
    Xc = X - (w[:, None] * X).sum(axis=0)
    return (w[:, None] * Xc).T @ Xc


# =================================================================================================== (a) triple collocation
def tc_errors(Cm):
    """The error variance of each of three sensors from their covariance matrix (errors independent of each other)."""
    return np.array([Cm[0, 0] - Cm[0, 1] * Cm[0, 2] / Cm[1, 2], Cm[1, 1] - Cm[0, 1] * Cm[1, 2] / Cm[0, 2], Cm[2, 2] - Cm[0, 2] * Cm[1, 2] / Cm[0, 1]])


def triple_collocation(S):
    rows = []
    for name, (rad_cols, a, b) in C["triplets"].items():
        rad = S["r15"][rad_cols].mean(axis=1, skipna=False)
        ga = S["s15"] if a == "GbgA" else S["g15"][a]
        gb = S["s15"] if b == "GbgA" else S["g15"][b]
        X = np.c_[rad.to_numpy(), ga.to_numpy(), gb.to_numpy()].astype(float)
        ev = S["ev15"]
        base = (ev >= 0) & np.isfinite(X).all(axis=1)
        for scale in ("log", "linear"):
            keep = base & ((X >= THR).all(axis=1) if scale == "log" else (X >= C["linear_min_any_mmh"]).any(axis=1))
            Y = np.log(X[keep]) if scale == "log" else X[keep]
            e = ev[keep]
            for centred in ((False, True) if scale == "log" else (False,)):
                Z, dof = centre(Y, e) if centred else (Y, 1.0)
                var = tc_errors(wcov(Z, np.ones(len(Z))) * dof)
                B = np.array([tc_errors(wcov(Z, w) * dof) for w in boot_weights(e, C["boot_reps"], C["seed"])])
                lo, hi = np.nanpercentile(B, C["interval"], axis=0)
                for k, sensor in enumerate(("radar", a, b)):
                    rows.append({"triplet": name, "sensors": f"radar@{'+'.join(rad_cols)}, {a}, {b}", "scale": scale, "event_centred": centred,
                                 "sensor": sensor, "bins": int(keep.sum()), "events": int(len(np.unique(e))), "error_var": var[k], "var_lo": lo[k], "var_hi": hi[k],
                                 "error_sd": np.sqrt(max(var[k], 0.0)), "sd_lo": np.sqrt(max(lo[k], 0.0)), "sd_hi": np.sqrt(max(hi[k], 0.0)),
                                 "sensor_var": float(np.var(Z[:, k]) * dof)})
    return pd.DataFrame(rows)


# =================================================================================================== (b) gauge pairs
def gauge_pairs(S):
    g = S["gauges"]
    city = list(S["g15"].columns)
    ev = S["ev15"]
    rows = []
    for a, b in itertools.combinations(city, 2):
        x, y = S["g15"][a].to_numpy(float), S["g15"][b].to_numpy(float)
        ok = (ev >= 0) & np.isfinite(x) & np.isfinite(y)
        both = ok & (x >= THR) & (y >= THR)
        dl = np.log(x[both]) - np.log(y[both])
        rows.append({"gauge_a": a, "gauge_b": b, "type_a": g.gauge_type[a], "type_b": g.gauge_type[b],
                     "separation_km": float(np.hypot(g.x[a] - g.x[b], g.y[a] - g.y[b]) / 1000), "bins": int(ok.sum()),
                     "corr_rate": float(np.corrcoef(x[ok], y[ok])[0, 1]), "bins_both_rain": int(both.sum()),
                     "semivar_log": float(0.5 * np.mean(dl ** 2)) if both.sum() >= C["pairs_min_bins"] else np.nan,
                     "corr_log": float(np.corrcoef(np.log(x[both]), np.log(y[both]))[0, 1]) if both.sum() >= C["pairs_min_bins"] else np.nan})
    t = pd.DataFrame(rows).sort_values("separation_km").reset_index(drop=True)
    m = t[(t.separation_km <= C["pairs_fit_max_km"]) & t.semivar_log.notna()]
    slope, nug = np.polyfit(m.separation_km, m.semivar_log, 1, w=np.sqrt(m.bins_both_rain))
    fit = {"pairs_in_fit": len(m), "nugget": float(nug), "slope_per_km": float(slope), "gauge_own_error_sd": float(np.sqrt(max(nug, 0.0))),
           "semivar_at_1km": float(nug + slope), "semivar_at_2km": float(nug + 2 * slope),
           "repr_point_vs_pixel_var": float(max(slope, 0.0) * C["pixel_mean_distance_km"]),
           "repr_point_vs_pixel_sd": float(np.sqrt(max(slope, 0.0) * C["pixel_mean_distance_km"]))}
    return t, fit


# =================================================================================================== (c) binning reversal
def reversal(S, res):
    """Conditional means both ways (the signature), and the unit-gain split of the event-centred log variances."""
    r, g, ev = S[f"r{res}"], S[f"g{res}"], S[f"ev{res}"]
    gauges = [c for c in g.columns if S["gauges"].gauge_type[c] == "weighing"]
    R = r[gauges].to_numpy(float)
    G = g[gauges].to_numpy(float)
    E = np.repeat(ev[:, None], len(gauges), axis=1)
    ok = (E >= 0) & np.isfinite(R) & np.isfinite(G)
    rr, gg, ee = R[ok], G[ok], E[ok]
    rows = []
    for lo, hi in zip(EDGES[:-1], EDGES[1:]):
        a, b = (rr >= lo) & (rr < hi), (gg >= lo) & (gg < hi)
        rows.append({"resolution_min": res, "class_mmh": f"{lo:g}-{hi:g}", "n_radar_in_class": int(a.sum()), "mean_radar": float(rr[a].mean()) if a.any() else np.nan,
                     "mean_gauge_given_radar": float(gg[a].mean()) if a.any() else np.nan, "n_gauge_in_class": int(b.sum()),
                     "mean_gauge": float(gg[b].mean()) if b.any() else np.nan, "mean_radar_given_gauge": float(rr[b].mean()) if b.any() else np.nan})
    t = pd.DataFrame(rows)
    t["gauge_over_radar_given_radar"] = t.mean_gauge_given_radar / t.mean_radar
    t["radar_over_gauge_given_gauge"] = t.mean_radar_given_gauge / t.mean_gauge
    both = (rr >= THR) & (gg >= THR)
    L = np.c_[np.log(rr[both]), np.log(gg[both])]
    e = ee[both]

    def split(w):
        Cm = wcov(Z, w) * dof
        return np.array([Cm[0, 0] - Cm[0, 1], Cm[1, 1] - Cm[0, 1], Cm[0, 0] + Cm[1, 1] - 2 * Cm[0, 1], Cm[0, 1]])
    Z, dof = centre(L, e)
    est = split(np.ones(len(Z)))
    B = np.array([split(w) for w in boot_weights(e, C["boot_reps"], C["seed"])])
    lo, hi = np.percentile(B, C["interval"], axis=0)
    names = ("radar_error_var", "gauge_error_var", "ratio_var", "shared_var")
    out = {"resolution_min": res, "bins": int(both.sum()), "events": int(len(np.unique(e)))}
    for k, nm in enumerate(names):
        out[nm], out[f"{nm}_lo"], out[f"{nm}_hi"] = est[k], lo[k], hi[k]
    # the lag-1 correlation of the event-centred ln(radar / gauge) (per gauge, consecutive bins both rainy)
    d = np.where((R >= THR) & (G >= THR) & (E >= 0), np.log(np.where(R > 0, R, 1.0)) - np.log(np.where(G > 0, G, 1.0)), np.nan)
    pairs = []
    for j in range(d.shape[1]):
        for ev_id in np.unique(ev[ev >= 0]):
            col = d[ev == ev_id, j]
            col = col - np.nanmean(col) if np.isfinite(col).any() else col
            m = np.isfinite(col[:-1]) & np.isfinite(col[1:])
            pairs.append(np.c_[col[:-1][m], col[1:][m]])
    pairs = np.concatenate(pairs)
    out["rho_lag1"], out["rho_pairs"] = float(np.corrcoef(pairs[:, 0], pairs[:, 1])[0, 1]), len(pairs)
    return t, out


# =================================================================================================== (d) scale part and shared part
def scale_part(S):
    """Per TRAIN event the pooled log ratio radar / gauge (the scale part), and per (event, gauge) its own log ratio."""
    r, g, ev = S["r5"], S["g5"], S["ev5"]
    city = list(g.columns)
    R, G = r[city].to_numpy(float), g[city].to_numpy(float)
    ok = np.isfinite(R) & np.isfinite(G)
    used = set(S["g_table"][S["g_table"].path_used].event_id)
    rows, per_gauge = [], []
    for i, e in enumerate(S["events"].itertuples()):
        m = (ev == i)[:, None] & ok
        rs, gs = np.where(m, R, 0.0).sum(axis=0) / 12, np.where(m, G, 0.0).sum(axis=0) / 12          # mm per gauge
        row = {"event_id": e.event_id, "gauge_ratio_used": e.event_id in used, "radar_mm_mean": float(rs.mean()), "gauge_mm_mean": float(gs.mean()),
               "ln_ratio": float(np.log(rs.sum() / gs.sum())) if rs.sum() > 0 and gs.sum() > 0 else np.nan}
        rows.append(row)
        for j, c in enumerate(city):
            if gs[j] >= C["scale_min_gauge_mm"] and rs[j] > 0 and np.isfinite(row["ln_ratio"]):
                per_gauge.append({"event_id": e.event_id, "gauge": c, "gauge_mm": gs[j], "ln_ratio_gauge_minus_event": float(np.log(rs[j] / gs[j]) - row["ln_ratio"])})
    return pd.DataFrame(rows), pd.DataFrame(per_gauge)


def shared_in_space(S):
    """cov(d_i, d_j) between City sites, d = event-centred ln(radar / gauge) at 15 min: the part of the radar error shared in space."""
    r, g, ev = S["r15"], S["g15"], S["ev15"]
    gz = S["gauges"]
    rows = []
    for a, b in itertools.combinations(list(g.columns), 2):
        X = np.c_[r[a], g[a], r[b], g[b]].astype(float)
        keep = (ev >= 0) & np.isfinite(X).all(axis=1) & (X >= THR).all(axis=1)
        if keep.sum() < C["pairs_min_bins"]:
            continue
        D = np.c_[np.log(X[keep, 0] / X[keep, 1]), np.log(X[keep, 2] / X[keep, 3])]
        Z, dof = centre(D, ev[keep])
        Cm = np.cov(Z.T) * dof
        rows.append({"gauge_a": a, "gauge_b": b, "separation_km": float(np.hypot(gz.x[a] - gz.x[b], gz.y[a] - gz.y[b]) / 1000), "bins": int(keep.sum()),
                     "cov_d": float(Cm[0, 1]), "var_d_a": float(Cm[0, 0]), "var_d_b": float(Cm[1, 1])})
    t = pd.DataFrame(rows)
    return t, float(np.average(t.cov_d, weights=t.bins)), float(np.average(0.5 * (t.var_d_a + t.var_d_b), weights=t.bins))


# =================================================================================================== figures
def fig_pairs(pairs, fit, shared):
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9))
    ax = axes[0]
    for typ, mk in (("weighing", "o"), ("tipping_bucket", "s")):
        m = pairs[(pairs.type_a == typ) & (pairs.type_b == typ)]
        ax.plot(m.separation_km, m.corr_rate, mk, color=INK, mfc="white", label=f"both {typ.replace('_', ' ')}")
    m = pairs[pairs.type_a != pairs.type_b]
    ax.plot(m.separation_km, m.corr_rate, "x", color=INK2, label="mixed types")
    ax.set_ylabel("correlation of the 15 min rates")
    ax.set_title("City gauge pairs: correlation", fontsize=9)
    ax = axes[1]
    ax.plot(pairs.separation_km, pairs.semivar_log, "o", color=INK, mfc="white", label="pairs")
    h = np.linspace(0, C["pairs_fit_max_km"], 20)
    ax.plot(h, fit["nugget"] + fit["slope_per_km"] * h, "--", color=INK, label=f"line up to {C['pairs_fit_max_km']:g} km: {fit['nugget']:.3f} + {fit['slope_per_km']:.3f} h")
    ax.axvspan(0, 2, color=INK2, alpha=0.12, hatch="//", lw=0)
    ax.text(1.0, ax.get_ylim()[1] * 0.92, "0-2 km:\nextrapolated", ha="center", va="top", fontsize=7)
    ax.set_ylabel("semivariance of ln(rate)")
    ax.set_title("City gauge pairs: ½ var(ln g_i − ln g_j), both >= 0.5 mm/h", fontsize=9)
    ax = axes[2]
    ax.plot(shared.separation_km, shared.cov_d, "o", color=INK, mfc="white", label="cov between sites")
    ax.axhline(np.average(shared.cov_d, weights=shared.bins), color=INK, ls="--", label="mean (weighted)")
    ax.axhline(0, color=INK2, lw=0.6)
    ax.set_ylabel("cov of ln(radar / gauge) between two sites")
    ax.set_title("radar error shared in space (15 min, event-centred)", fontsize=9)
    for ax in axes:
        ax.set_xlabel("separation (km)")
        ax.legend(fontsize=7)
        style_axes(ax)
    fig.tight_layout()
    return fig


def fig_reversal(rev5, rev15):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, t, res in ((axes[0], rev5, 5), (axes[1], rev15, 15)):
        ax.plot(t.mean_radar, t.mean_gauge_given_radar, "o-", color=INK, mfc="white", label="binned by RADAR: mean gauge")
        ax.plot(t.mean_radar_given_gauge, t.mean_gauge, "s--", color=INK, mfc="white", label="binned by GAUGE: mean radar (x) against gauge (y)")
        lim = [0.05, 30]
        ax.plot(lim, lim, ":", color=INK2, label="1 : 1")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("radar rain (mm/h)")
        ax.set_title(f"binning reversal, {res} min, weighing City gauges", fontsize=9)
        ax.legend(fontsize=7)
        style_axes(ax)
    axes[0].set_ylabel("gauge rain (mm/h)")
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    S = load()
    print(f"TRAIN events {len(S['events'])}; 5 min bins inside them {int((S['ev5'] >= 0).sum())}; 15 min bins {int((S['ev15'] >= 0).sum())}", flush=True)
    tc = triple_collocation(S)
    pairs, fit = gauge_pairs(S)
    rev5_t, rev5 = reversal(S, 5)
    rev15_t, rev15 = reversal(S, 15)
    ev_scale, gauge_scale = scale_part(S)
    shared, cov_shared, var_d15 = shared_in_space(S)

    # ---- (d) the synthesis, exactly as written in decisions.md before the run
    rho = rev5["rho_lag1"]
    f15 = (3 + 4 * rho + 2 * rho ** 2) / 9                      # variance of a 15 min mean of an AR(1) 5 min error, relative to one bin
    t = tc[(tc.scale == "log") & tc.event_centred & (tc.sensor == "radar")]
    tc_var15 = float(np.median(t.error_var))
    tc_lo15, tc_hi15 = float(np.median(t.var_lo)), float(np.median(t.var_hi))
    est_tc = np.sqrt(max(tc_var15, 0.0) / f15)
    tc_lo, tc_hi = np.sqrt(max(tc_lo15, 0.0) / f15), np.sqrt(max(tc_hi15, 0.0) / f15)
    est_rev = np.sqrt(max(rev5["radar_error_var"], 0.0))
    rev_lo, rev_hi = np.sqrt(max(rev5["radar_error_var_lo"], 0.0)), np.sqrt(max(rev5["radar_error_var_hi"], 0.0))
    small, large = ((est_tc, tc_lo, tc_hi), (est_rev, rev_lo, rev_hi)) if est_tc <= est_rev else ((est_rev, rev_lo, rev_hi), (est_tc, tc_lo, tc_hi))
    sigma_lo, sigma_hi = small[1], large[2]
    all_to_radar = float(np.sqrt(max(rev5["ratio_var"], 0.0)))
    es = ev_scale[(ev_scale.gauge_mm_mean >= C["scale_min_gauge_mm"]) & ev_scale.ln_ratio.notna()]
    sigma_s, sigma_s_res = float(es.ln_ratio.std(ddof=1)), float(gauge_scale.ln_ratio_gauge_minus_event.std(ddof=1))
    light = ev_scale[~ev_scale.gauge_ratio_used]
    syn = pd.DataFrame([
        {"quantity": "rho_lag1 of the 5 min ln(radar / gauge), event-centred", "value": rho, "lo": np.nan, "hi": np.nan},
        {"quantity": "factor var(15 min mean) / var(5 min bin)", "value": f15, "lo": np.nan, "hi": np.nan},
        {"quantity": "est_TC: per-bin radar error sd at 5 min (from 15 min triple collocation)", "value": est_tc, "lo": tc_lo, "hi": tc_hi},
        {"quantity": "est_rev: per-bin radar error sd at 5 min (unit-gain split)", "value": est_rev, "lo": rev_lo, "hi": rev_hi},
        {"quantity": "sigma_lo (per-bin part, 5 min, one pixel)", "value": sigma_lo, "lo": np.nan, "hi": np.nan},
        {"quantity": "sigma_hi (per-bin part, 5 min, one pixel)", "value": sigma_hi, "lo": np.nan, "hi": np.nan},
        {"quantity": "all-to-radar limit: sd of ln(radar / gauge) at 5 min, event-centred", "value": all_to_radar, "lo": np.nan, "hi": np.nan},
        {"quantity": "gauge error sd at 5 min (unit-gain split)", "value": np.sqrt(max(rev5["gauge_error_var"], 0.0)), "lo": np.nan, "hi": np.nan},
        {"quantity": "sigma_s: scale part, sd over events of ln(sum radar / sum gauge)", "value": sigma_s, "lo": np.nan, "hi": np.nan},
        {"quantity": "median over events of ln(sum radar / sum gauge)", "value": float(es.ln_ratio.median()), "lo": np.nan, "hi": np.nan},
        {"quantity": "events in sigma_s", "value": len(es), "lo": np.nan, "hi": np.nan},
        {"quantity": "sigma_s_res: scale part left after g_e (gauge's own event ratio minus the pooled one)", "value": sigma_s_res, "lo": np.nan, "hi": np.nan},
        {"quantity": "shared-in-space variance of the per-bin part (15 min, mean pair covariance)", "value": cov_shared, "lo": np.nan, "hi": np.nan},
        {"quantity": "variance of ln(radar / gauge) at 15 min on the same bins", "value": var_d15, "lo": np.nan, "hi": np.nan},
        {"quantity": "gauge own error sd at 15 min (pair variogram at h = 0)", "value": fit["gauge_own_error_sd"], "lo": np.nan, "hi": np.nan},
        {"quantity": "point-vs-pixel representativeness sd at 15 min (pair variogram)", "value": fit["repr_point_vs_pixel_sd"], "lo": np.nan, "hi": np.nan},
        {"quantity": "g_e = 1 events: number", "value": len(light), "lo": np.nan, "hi": np.nan},
        {"quantity": "g_e = 1 events: pooled radar mm per gauge", "value": float(light.radar_mm_mean.sum()), "lo": np.nan, "hi": np.nan},
        {"quantity": "g_e = 1 events: pooled gauge mm per gauge", "value": float(light.gauge_mm_mean.sum()), "lo": np.nan, "hi": np.nan},
    ])
    print(syn.round(4).to_string(), flush=True)
    print(tc[(tc.scale == "log")][["triplet", "event_centred", "sensor", "bins", "events", "error_sd", "sd_lo", "sd_hi"]].round(3).to_string(), flush=True)
    for tab, f in ((tc, P.E1_TC), (pairs, P.E1_PAIRS), (pd.concat([rev5_t, rev15_t]), P.E1_REVERSAL), (pd.DataFrame([rev5, rev15]), P.E1_SPLIT),
                   (ev_scale, P.E1_SCALE), (shared, P.E1_SHARED), (syn, P.E1_SYNTHESIS), (pd.DataFrame([fit]), P.E1_PAIR_FIT)):
        save(tab.round(5), f, index=False)
    for fig, f in ((fig_pairs(pairs, fit, shared), P.E1_FIG_PAIRS), (fig_reversal(rev5_t, rev15_t), P.E1_FIG_REVERSAL)):
        save(fig, f, dpi=120)
        plt.close(fig)

    a = ["# E1: the radar's own rain error (rain data only, TRAIN events)\n",
         f"Written by `research/waa_ksurface/e1_rain_error.py`. {len(S['events'])} TRAIN events, {int((S['ev5'] >= 0).sum())} five-minute bins. No attenuation data. "
         f"Errors are multiplicative: sd of ln(rate), on bins where every sensor of the comparison reads >= {THR} mm/h. 90 % intervals from {C['boot_reps']} event draws.\n",
         "## (a) Triple collocation (15 min)\n",
         "Error sd of ln(rate) of each sensor, with no sensor taken as the truth. `event_centred = True` is the per-bin part (event means removed).\n",
         md_table(tc[tc.scale == "log"][["triplet", "sensors", "event_centred", "sensor", "bins", "events", "error_sd", "sd_lo", "sd_hi"]].round(3)),
         "\nLinear version (error sd in mm/h, bins where any sensor reads >= 0.1 mm/h):\n",
         md_table(tc[tc.scale == "linear"][["triplet", "sensor", "bins", "error_sd", "sd_lo", "sd_hi"]].round(3)),
         f"\n(`{rel(P.E1_TC)}`)\n",
         "## (b) Gauge pairs against separation (15 min)\n",
         f"Line through the semivariance of ln(rate) of the {fit['pairs_in_fit']} City pairs up to {C['pairs_fit_max_km']:g} km: γ(h) = {fit['nugget']:.3f} + {fit['slope_per_km']:.3f}·h. "
         f"**Extrapolated to 0–2 km** (the closest pair is 1.5 km apart): γ(0) = {fit['nugget']:.3f} (gauge's own error sd {fit['gauge_own_error_sd']:.2f}), "
         f"γ(1 km) = {fit['semivar_at_1km']:.3f}, γ(2 km) = {fit['semivar_at_2km']:.3f}. Point-vs-pixel representativeness: sd {fit['repr_point_vs_pixel_sd']:.2f} "
         f"(variance {fit['repr_point_vs_pixel_var']:.3f}). (`{rel(P.E1_PAIRS)}`, `{rel(P.E1_FIG_PAIRS)}`)\n",
         md_table(pairs[["gauge_a", "gauge_b", "separation_km", "corr_rate", "bins_both_rain", "semivar_log"]].head(12).round(3)),
         "\n## (c) The binning reversal (weighing City gauges)\n",
         md_table(pd.concat([rev5_t, rev15_t])[["resolution_min", "class_mmh", "n_radar_in_class", "mean_radar", "mean_gauge_given_radar", "gauge_over_radar_given_radar",
                                                "n_gauge_in_class", "mean_gauge", "mean_radar_given_gauge", "radar_over_gauge_given_gauge"]].round(3)),
         "\nThe unit-gain split of the event-centred log variances (bins with both >= 0.5 mm/h):\n",
         md_table(pd.DataFrame([rev5, rev15])[["resolution_min", "bins", "events", "radar_error_var", "radar_error_var_lo", "radar_error_var_hi", "gauge_error_var",
                                               "ratio_var", "shared_var", "rho_lag1"]].round(4)),
         f"\n(`{rel(P.E1_REVERSAL)}`, `{rel(P.E1_SPLIT)}`, `{rel(P.E1_FIG_REVERSAL)}`)\n",
         "## (d) The result\n", md_table(syn.round(4)), f"\n(`{rel(P.E1_SYNTHESIS)}`, `{rel(P.E1_SCALE)}`, `{rel(P.E1_SHARED)}`)\n"]
    save("\n".join(a) + "\n", P.E1_REPORT)
    print(f"report: {rel(P.E1_REPORT)}", flush=True)


if __name__ == "__main__":
    main()
