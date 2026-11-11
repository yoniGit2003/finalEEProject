"""Step T1 of the k-surface branch: the rain error that the length test actually meets.

    python -m research.waa_ksurface.t1_conditional_error

The length test does not use every 5 min bin. It uses bins in steady rain (30 min within a factor 1.5) and labels them by a
median. So the error that matters is the error of that label in such windows, not the error of a single bin (step E1).
(a) radar label against gauge rain in steady 35 min windows at the City gauges (rain data only);
(b) the count of steady bins: how fast it falls when a per-bin error is added, and the largest fast error that a perfectly
    constant rain could carry and still give the real share of steady bins.
No attenuation data is read in (a); (b) reads only the rain of the TRAIN link-events.
"""
import numpy as np
import pandas as pd

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import data, fakeworld as fw, ltest
from research.waa_ksurface import paths as P

LT, FK = ltest.config(), fw.config()
C = data.branch_config()["t1"]
RULE = LT["rules"]["R2"]                                        # the adopted steadiness rule (30 min, factor 1.5)
EDGES = np.asarray(LT["class_edges_mmh"], float)
CAVEAT = ("the σ that matches the real count is an EFFECTIVE bound, biased low because the fake truth (the observed radar) "
          "already carries error")


# =================================================================================================== (a) steady windows at the gauges
def windows_at_gauges(steady_only=True):
    """One row per 35 min window at a City gauge inside the gauge-ratio TRAIN events: label (median radar / g_e), gauge mean."""
    d = data.open_release()
    g = data.g_table()
    g = g[g.path_used].set_index("event_id").g_path
    ev = d.events().set_index("event_id").loc[g.index]
    rad = d.store("radar_at_gauges_5min").R_pixel.load().to_pandas()
    c1 = d.store("gauges_city_1min")
    g1 = c1.rate_mm_h.load().to_pandas().where(c1.n_valid.load().to_pandas() >= 1)
    k = RULE["t_pre_min"] // 5
    rows = []
    for e in ev.itertuples():
        r5 = rad[(rad.index >= e.start) & (rad.index < e.end)] / g[e.Index]
        gg = g1[(g1.index >= e.start) & (g1.index < e.end)]
        g5 = gg.resample("5min").mean().where(gg.notna().resample("5min").sum() == 5).reindex(r5.index)
        for gauge in g1.columns:
            r = r5[gauge].to_numpy(float)
            gv = g5[gauge].to_numpy(float)
            n = len(r)
            cls = ltest._steady(r, np.ones(n, bool), np.array([0]), np.array([n]), EDGES, k, float(RULE["q"]))
            for t in range(k, n):
                win_r, win_g = r[t - k:t + 1], gv[t - k:t + 1]
                if not (np.isfinite(win_r).all() and np.isfinite(win_g).all()):
                    continue
                if steady_only and cls[t] < 0:
                    continue
                rows.append({"event_id": e.Index, "gauge": gauge, "time": r5.index[t], "steady": bool(cls[t] >= 0), "label": float(np.median(win_r)),
                             "gauge_mean": float(win_g.mean())})
    return pd.DataFrame(rows)


def components(w):
    """sd of d = ln(label / gauge) and its split: between events (scale part) and within events."""
    d = np.log(w.label / w.gauge_mean)
    between, within, k, n = ltest.variance_components(d.to_numpy(), w.event_id.to_numpy())
    return {"windows": n, "events": k, "sd_total": float(d.std(ddof=1)), "median_d": float(d.median()), "sigma_ev": between, "sigma_w": within,
            "sigma_cond": float(np.sqrt(1 - C["shared_share"]) * within), "sigma_shared": float(np.sqrt(C["shared_share"]) * within)}


def conditional_error():
    w_all = windows_at_gauges(steady_only=False)
    out, parts = {}, {}
    for name, w in (("steady windows", w_all[w_all.steady]), ("all windows (reference)", w_all)):
        big = w[w.label >= C["min_label_mmh"]]
        used = big[big.gauge_mean > 0]
        row = {"set": name, "windows_label_ge_min": len(big), "dropped_gauge_zero": int((big.gauge_mean <= 0).sum())}
        row.update(components(used) if len(used) >= 3 else {"windows": len(used), "events": used.event_id.nunique()})
        out[name], parts[name] = row, used
    used = parts["steady windows"]
    ok = len(used) >= C["min_windows"] and used.event_id.nunique() >= C["min_events"]
    lo = hi = np.nan
    if ok:
        rng = np.random.default_rng(C["seed"])
        ids = used.event_id.unique()
        grp = {e: used[used.event_id == e] for e in ids}
        draws = []
        for _ in range(C["boot_reps"]):
            pick = rng.choice(ids, len(ids))
            b = pd.concat([grp[e].assign(event_id=f"{e}_{i}") for i, e in enumerate(pick)])
            draws.append(components(b)["sigma_cond"])
        lo, hi = np.nanpercentile(draws, C["interval"])
    by_class = []
    for lo_c, hi_c in zip(EDGES[:-1], EDGES[1:]):
        m = used[(used.label >= lo_c) & (used.label < hi_c)]
        if len(m) >= 10:
            dd = np.log(m.label / m.gauge_mean)
            by_class.append({"label_class_mmh": f"{lo_c:g}-{hi_c:g}", "windows": len(m), "events": m.event_id.nunique(), "gauges": m.gauge.nunique(),
                             "median_d": float(dd.median()), "sd_d": float(dd.std(ddof=1))})
    return pd.DataFrame(list(out.values())), pd.DataFrame(by_class), used, ok, float(lo), float(hi)


# =================================================================================================== (b) the steady-window count
def count_check():
    D, le, n_cells = data.gauge_ratio_train()
    n5 = D["le_n"] // 5
    r = D["r5"].astype(float)
    ok = np.isfinite(r)
    k, q = RULE["t_pre_min"] // 5, float(RULE["q"])
    ev = fw.time_index(le)[1]
    ev_bin = np.repeat(ev, n5)

    def count(rr):
        cls = ltest._steady(rr, ok, D["le_boff5"], n5, EDGES, k, q)
        return np.bincount(ev_bin[cls >= 0], minlength=ev.max() + 1)
    real = count(r)
    rows = []
    for sig in C["count_grid"]:
        reps = []
        for rep in range(C["count_reps"]):
            rng = np.random.default_rng(C["seed"] + 10 + rep)
            err = fw.rain_error(le, n_cells, sig, FK, rng, scale_sd=0.0) if sig > 0 else np.zeros(len(r))
            reps.append(count(r * np.exp(err)))
        c = np.mean(reps, axis=0)
        rows.append({"added_per_bin_sd": sig, "steady_bins": float(c.sum()), "share_of_real": float(c.sum() / real.sum()),
                     "events_with_fewer_than_half": int((c < 0.5 * real).sum()), "median_event_share": float(np.median(c[real > 0] / real[real > 0]))})
    # second reading: a constant truth plus a per-bin AR(1) error; which sd gives the real share of steady bins among rainy bins?
    in_class = ok & (r >= EDGES[0]) & (r < EDGES[-1])
    cls = ltest._steady(r, ok, D["le_boff5"], n5, EDGES, k, q)
    pos = np.concatenate([np.arange(n) for n in n5])
    cand = in_class & (pos >= k)                                # rainy bins that have 30 min before them inside the window
    share_real = float((cls >= 0).sum() / cand.sum())
    rng = np.random.default_rng(C["seed"] + 99)
    n = 400_000
    e = fw._ar1(rng.standard_normal(n), FK["rain_error"]["rho"])
    grid = np.linspace(0.02, 0.6, 59)
    shares = []
    for s in grid:
        x = s * e
        okk = np.ones(n - k, bool)
        for j in range(1, k + 1):
            okk &= np.abs(x[k:] - x[k - j:n - j]) <= np.log(q)
        shares.append(okk.mean())
    shares = np.array(shares)
    bound = float(np.interp(share_real, shares[::-1], grid[::-1]))
    return pd.DataFrame(rows), {"real_steady_bins": int(real.sum()), "events": int((real > 0).sum()), "rainy_bins_with_history": int(cand.sum()),
                                "real_share_steady": share_real, "constant_truth_upper_bound_sd": bound}, real


def main():
    summ, by_class, used, ok, lo, hi = conditional_error()
    print(summ.round(4).to_string(), flush=True)
    print(by_class.round(3).to_string(), flush=True)
    counts, info, real = count_check()
    print(counts.round(3).to_string(), flush=True)
    print(info, flush=True)
    s = summ.iloc[0]
    if ok:
        sc_lo, sc_hi, source = lo, hi, "(a) steady windows at the gauges, event bootstrap"
    else:
        sc_lo, sc_hi, source = 0.0, info["constant_truth_upper_bound_sd"], "(c) fallback: the constant-truth bound of (b)"
    res = pd.DataFrame([{"quantity": "sigma_c_lo", "value": sc_lo}, {"quantity": "sigma_c_hi", "value": sc_hi},
                        {"quantity": "sigma_cond (point)", "value": s.get("sigma_cond", np.nan)}, {"quantity": "sigma_shared (point)", "value": s.get("sigma_shared", np.nan)},
                        {"quantity": "sigma_ev (point)", "value": s.get("sigma_ev", np.nan)}, {"quantity": "sigma_w (point)", "value": s.get("sigma_w", np.nan)},
                        {"quantity": "sd_total", "value": s.get("sd_total", np.nan)}, {"quantity": "median_d", "value": s.get("median_d", np.nan)},
                        {"quantity": "windows", "value": s["windows"]}, {"quantity": "events", "value": s["events"]},
                        {"quantity": "fast_error_upper_bound (constant truth)", "value": info["constant_truth_upper_bound_sd"]},
                        {"quantity": "source", "value": source}])
    for t, f in ((summ, P.T1_SUMMARY), (by_class, P.T1_BY_CLASS), (used, P.T1_WINDOWS), (counts, P.T1_COUNTS), (res, P.T1_RESULT)):
        save(t.round(5) if f is not P.T1_RESULT else t, f, index=False)
    print(res.to_string(), flush=True)

    a = ["# T1: the rain error in steady windows (the error the length test meets)\n",
         f"Written by `research/waa_ksurface/t1_conditional_error.py`. Steadiness rule: {RULE['t_pre_min']} min within a factor {RULE['q']:g}. Part (a) reads rain data only "
         "(radar pixel at the 10 City gauges divided by g_e, the 1 min gauges), inside the 14 TRAIN events with a gauge ratio. VAL and TEST were not opened.\n",
         "## (a) Radar label against gauge rain in steady 35 min windows\n",
         f"d = ln(label / gauge); label = the median of the 7 radar values of the window, gauge = its mean over the same 35 min. Used: label >= {C['min_label_mmh']:g} mm/h and "
         "gauge > 0. `sigma_ev` = between events (the scale part), `sigma_w` = within events, `sigma_cond` = √0.54·sigma_w (the link-independent part), "
         "`sigma_shared` = √0.46·sigma_w.\n",
         md_table(summ.round(3)),
         f"\n**[σc_lo, σc_hi] = [{sc_lo:.3f}, {sc_hi:.3f}]** (source: {source}). The gauge's own error and the point-vs-pixel difference are not subtracted. "
         f"**Not measured below {C['min_label_mmh']:g} mm/h**: there a 35 min window holds only a few gauge steps.\n",
         "By label class (steady windows):\n", md_table(by_class.round(3)) if len(by_class) else "too few windows per class",
         f"\n(`{rel(P.T1_SUMMARY)}`, `{rel(P.T1_WINDOWS)}`, `{rel(P.T1_RESULT)}`)\n",
         "## (b) The count of steady bins (TRAIN link-events, gauge-ratio events)\n",
         f"Real data: {info['real_steady_bins']} steady bins in {info['events']} events. Fake replicates with the observed radar as the true rain and a per-bin error added "
         "(46 % shared, lag-1 correlation 0.35):\n",
         md_table(counts.round(3)),
         f"\n**Caveat: {CAVEAT}.** With the observed radar as the truth, the count at added error 0 is the real count by construction, so the match is at 0. The table shows "
         "how fast the count falls: the real radar cannot carry much MORE fast error than it already has.\n",
         f"Second reading: {100 * info['real_share_steady']:.1f} % of the rainy bins with 30 min of history are steady in the real data ({info['rainy_bins_with_history']} bins). "
         f"A rain that is perfectly constant in time gives this share with a per-bin error of sd **{info['constant_truth_upper_bound_sd']:.3f}** (lag-1 correlation 0.35). "
         "This is an upper bound on the fast per-bin error of a link's own radar value: real changes of the rain use up part of the allowed unsteadiness. "
         f"(`{rel(P.T1_COUNTS)}`)\n"]
    save("\n".join(a) + "\n", P.T1_REPORT)
    print(f"report: {rel(P.T1_REPORT)}", flush=True)


if __name__ == "__main__":
    main()
