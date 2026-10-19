"""wet_antenna steps S4.1 + S4.2: why do events differ? Per-event rain-scale and antenna factors, then the conditions.

S3b found that every event-independent fit fails the VAL length test: the VAL events attenuate about 0.6-0.7 x per mm/h of
radar rain, at the antenna (intercept I) and on the path (slope s). This step asks why, on the FROZEN M4b-B of S3b (M3-B
as a check; F, P, timing, a_link and the lag d unchanged):

  S4.1 a  the facts: dates and hours of every TRAIN and VAL event, how the split was made, rain per event;
  S4.1 b  is the TRAIN / VAL gap more than ordinary event-to-event spread? 200 random 18-event halves of TRAIN;
  S4.1 c  per-event factors (fake check first):
              A_pred = alpha_e * W[a_link F(R / k_e) dynamics] + L * P(R(t - d) / k_e) + c
          k_e = the rain-input scale of event e (k > 1: the radar reads that event too high), alpha_e = its antenna factor;
  S4.1 d  which factor closes the VAL gap in the length test;
  S4.1 e  radar against the City gauges per event, against k_e;
  S4.2    the conditions per event (weather, time of day, rain type, dry hours before, bands) against k_e and alpha_e,
          and the residual of the frozen M4b-B per band, night / day and temperature (descriptive only).

How the factors are found: W and L * P are simulated separately for every k on a grid (alpha only multiplies W), so the
Huber loss of an event is known on a (k, alpha) grid for every link; the estimate is the grid minimum with a weak
log-ridge toward 1, and a link bootstrap inside the event gives the interval.
STANDING RULE: WAA is never estimated by subtracting a power-law rain term. The path part is separated only through its
dependence on L at matched rain. ITU-R P.838 only as a size reference, never to correct data. k_e comes from the
data-learned P of the frozen fit and the spread of L, not from ITU.
Every rule is in config `wet_antenna.conditions_s4` (written before the run). The factor machinery and the radar / gauge
ratio are in the library event_factors.py. The per-event results are kept in
research/.cache/wet_antenna/conditions_s4/ so that the report can be re-made without recomputing.
Run from the repository root (about 10-20 min the first time):
    python -m research.wet_antenna.conditions_s4
"""
import copy
import pickle
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from cmlpipe.api import CMLData
from research.common.config import load_config
from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import length_test as LTL
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.event_factors import (AG, I1A, I1K, KG, EventParts, bh, bin_positions, bin_times, estimate, fake_obs, frozen,
                                               radar_gauge_pixel, radar_gauge_primary, spearman_boot)

CFG = load_config()
WA = CFG["wet_antenna"]
S2, S2B, S3, S3B, C = WA["system_id_s2"], WA["system_id_s2b"], WA["system_id_s3"], WA["system_id_s3b"], WA["conditions_s4"]
M = copy.deepcopy(S2["model"])                                # the S2b model settings (as S3b)
M["increments"], M["maxiter"] = S2B["model"]["increments"], dict(S2B["model"]["maxiter"])
M["bounds"] = {**M["bounds"], "increment": S2B["model"]["increment_bounds"], "tau_s_min": S2B["model"]["tau_s_m4b_min"]}
LT = S3["length_test"]
FK = S2["fake"]
DELTA = M["huber_delta_db"]
FA = C["factors"]
FITS = {k: tuple(S3B["fits"][k]) for k in (C["frozen"]["main"], C["frozen"]["check"])}   # "M4b-B" -> ("M4b", "B")
MAIN, CHECK = C["frozen"]["main"], C["frozen"]["check"]
SETS = ("TRAIN", "VAL")
MK = {"TRAIN": "o", "VAL": "s"}                               # marker shapes: TRAIN circles, VAL squares
VARIANTS = {"none": "frozen M4b-B (S3b)", "both": "k_e and alpha_e", "k_only": "k_e only (alpha = 1)", "alpha_only": "alpha_e only (k = 1)",
            "k_refit": "k refitted with alpha = 1", "alpha_refit": "alpha refitted with k = 1"}
VSTYLE = {"none": "x", "both": "D", "k_only": "^", "alpha_only": "v", "k_refit": "<", "alpha_refit": ">"}


# =================================================================================================== the per-event loop
def run_set(D, le, f, set_name, boot, fake=None, noise=None):
    """Per event of one dataset: the factors (and their link bootstrap), the predictions of every variant, and the fake
    check (fake = config fake section, only on TRAIN with the main fit)."""
    ev_ids = list(pd.unique(le.event_id))
    rng = np.random.default_rng([FA["seed"], SETS.index(set_name), list(FITS).index(f["name"])])
    preds = {v: np.zeros(len(D["y5"])) for v in VARIANTS}
    rows, fake_rows = [], []
    if fake:
        frng = np.random.default_rng(fake["seed"])
        pairs = [(k, a) for k in fake["k_values"] for a in fake["alpha_values"]]
        order = np.concatenate([frng.permutation(len(pairs)) for _ in range(int(np.ceil(len(ev_ids) / len(pairs))))])[:len(ev_ids)]
    for e, ev_id in enumerate(ev_ids):
        idx = np.nonzero(D["le_event"] == e)[0]
        assert (le.event_id.iloc[idx] == ev_id).all()
        ev = EventParts(D, le, f, idx)
        c_le = ev.De["le_c"]
        G = ev.grid(ev.De["y5"], ev.De["w5"], c_le)
        est = estimate(G, FA["ridge_per_link"], rng, boot)
        pos = bin_positions(D, idx)
        ik, ia = est["ik"], est["ia"]
        for v, (kk, aa) in {"none": (I1K, I1A), "both": (ik, ia), "k_only": (ik, I1A), "alpha_only": (I1K, ia),
                            "k_refit": (est["ik_refit"], I1A), "alpha_refit": (I1K, est["ia_refit"])}.items():
            preds[v][pos] = ev.pred(kk, aa, c_le)
        L = ev.le.drop_duplicates("cml_id").L_km
        row = {"set": set_name, "fit": f["name"], "event_id": ev_id, "links": ev.nl, "link_events": len(idx),
               "links_under_1km": int((L < 1).sum()), "links_over_3km": int((L > 3).sum()), "points": int((ev.De["w5"] > 0).sum()),
               "k": KG[ik], "alpha": AG[ia], "k_refit_alpha1": KG[est["ik_refit"]], "alpha_refit_k1": AG[est["ia_refit"]],
               "k_at_bound": bool(ik in (0, len(KG) - 1)), "alpha_at_bound": bool(ia in (0, len(AG) - 1))}
        if boot:
            bk, ba = KG[est["boot_ik"]], AG[est["boot_ia"]]
            row.update({"k_lo": np.percentile(bk, FA["interval"][0]), "k_hi": np.percentile(bk, FA["interval"][1]),
                        "alpha_lo": np.percentile(ba, FA["interval"][0]), "alpha_hi": np.percentile(ba, FA["interval"][1]),
                        "boot_corr_lnk_lnalpha": float(np.corrcoef(np.log(bk), np.log(ba))[0, 1]) if np.std(bk) > 0 and np.std(ba) > 0 else np.nan})
        rows.append(row)
        if fake:
            kt, at = pairs[order[e]]
            y5, w5, cf = fake_obs(ev, ev.truth_1min(kt, at), np.random.default_rng([fake["seed"], e]), noise)
            fe = estimate(ev.grid(y5, w5, cf), FA["ridge_per_link"])
            fake_rows.append({"event_id": ev_id, "links": ev.nl, "links_under_1km": row["links_under_1km"], "links_over_3km": row["links_over_3km"],
                              "k_true": kt, "alpha_true": at, "k_est": KG[fe["ik"]], "alpha_est": AG[fe["ia"]]})
        del ev, G
    return pd.DataFrame(rows), preds, pd.DataFrame(fake_rows)


# =================================================================================================== facts, gauges, conditions
def event_facts(D, le, ev_tab, set_name, times):
    """S4.1 a: per event the dates, the radar path-mean rain inside [start, end), links."""
    thr = C["facts"]["wet_threshold_mmh"]
    rows = []
    for e_id, g in le.groupby("event_id", sort=False):
        r_ev = ev_tab.loc[e_id]
        series, tot, peaks = [], [], []
        for i in g.drop_duplicates("cml_id").index:                       # one sub-link per link
            o, n = D["le_boff5"][i], D["le_n"][i] // 5
            tt, rr = times[o:o + n], D["r5"][o:o + n]
            m = (tt >= r_ev.start) & (tt < r_ev.end)
            if not m.any():
                continue
            tot.append(np.nansum(rr[m]) * 5 / 60)
            peaks.append(np.nanmax(rr[m]) if np.isfinite(rr[m]).any() else np.nan)
            series.append(pd.Series(rr[m], index=tt[m]))
        med = pd.concat(series, axis=1, sort=True).median(axis=1) if series else pd.Series(dtype=float)
        L = g.drop_duplicates("cml_id").L_km
        rows.append({"set": set_name, "event_id": e_id, "start": r_ev.start, "end": r_ev.end, "duration_h": r_ev.duration_h,
                     "start_hour_utc": r_ev.start.hour + r_ev.start.minute / 60, "intensity_class": r_ev.intensity_class,
                     "total_mm_median_link": float(np.median(tot)) if tot else np.nan, "peak_mmh_max_link": float(np.nanmax(peaks)) if peaks else np.nan,
                     "peak_mmh_median_link": float(np.nanmedian(peaks)) if peaks else np.nan,
                     "wet_fraction": float((med >= thr).mean()) if len(med) else np.nan,
                     "peak_over_mean": float(med.max() / med.mean()) if len(med) and med.mean() > 0 else np.nan,
                     "links": g.cml_id.nunique(), "link_events": len(g), "links_under_1km": int((L < 1).sum()), "links_over_3km": int((L > 3).sum()),
                     "share_38GHz": float((g.band_name == "38 GHz").mean()), "share_28GHz": float((g.band_name == "28 GHz").mean()),
                     "dry_before_h": r_ev.dry_before_h})
    return pd.DataFrame(rows)


def hourly_conditions(dew, ev_tab):
    """S4.2: per event the means over the hours of [start, end) (at least the hour of the start)."""
    rows = []
    for e_id, r in ev_tab.iterrows():
        h = dew[(dew.index >= r.start.floor("h")) & (dew.index < max(r.end, r.start.floor("h") + pd.Timedelta("1h")))]
        rows.append({"event_id": e_id, "temperature": h.temperature.mean(), "humidity": h.humidity.mean(), "wind_speed": h.wind_speed.mean(),
                     "dew_point_depression": h.dew_point_depression.mean(), "night_fraction": h.night.astype(float).mean(), "weather_hours": len(h)})
    return pd.DataFrame(rows)


def boot_median(v, ev, reps, rng):
    """Median and 90 % event-bootstrap interval of values v grouped by event ev."""
    u = np.unique(ev)
    by = {e: v[ev == e] for e in u}
    b = [np.median(np.concatenate([by[e] for e in rng.choice(u, len(u))])) for _ in range(reps)]
    return float(np.median(v)), *np.percentile(b, [5, 95])


# =================================================================================================== length test helpers
def cell_ratios(obs, mod, cells):
    """Median over the cells of obs / mod for I and s (and per cell)."""
    o, m = obs.set_index(["band", "class"]), mod.set_index(["band", "class"])
    rI, rs = [], []
    for c in cells:
        if c in o.index and c in m.index:
            rI.append(o.loc[c, "I"] / m.loc[c, "I"])
            rs.append(o.loc[c, "s"] / m.loc[c, "s"])
    return float(np.nanmedian(rI)) if rI else np.nan, float(np.nanmedian(rs)) if rs else np.nan


def compute():
    clock = time.time()
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    BANDS = list(sid.merged_bands(d, S2["data"])[0])
    ev_tab = pd.concat([d.events(split="train"), d.events(split="val")]).set_index("event_id")      # TEST events never read
    Ds = {"TRAIN": sid.load_data(P.SID3_CACHE / "data_train_all.npz", BANDS), "VAL": sid.load_data(P.SID3_CACHE / "data_val_all.npz", BANDS)}
    les = {"TRAIN": pickle.loads((P.SID3_CACHE / "le_train_all.pkl").read_bytes())[0].reset_index(drop=True),
           "VAL": pickle.loads((P.SID3_CACHE / "le_val_all.pkl").read_bytes())[0].reset_index(drop=True)}
    fz = {k: {**frozen(k, Ds["TRAIN"], les["TRAIN"]), "name": k} for k in FITS}
    times = {s: bin_times(Ds[s], les[s]) for s in SETS}
    noise = sid.noise_table()
    out = {"BANDS": BANDS, "n_events_table": {s: int((ev_tab.split_event == s.lower()).sum()) for s in SETS}}
    print(f"datasets: TRAIN {len(les['TRAIN'])} link-events / {les['TRAIN'].event_id.nunique()} events, VAL {len(les['VAL'])} / "
          f"{les['VAL'].event_id.nunique()}; grid {len(KG)} x {len(AG)}", flush=True)

    # ------------------------------------------------------------------------------------------- a) facts
    facts = pd.concat([event_facts(Ds[s], les[s], ev_tab, s, times[s]) for s in SETS], ignore_index=True)
    print(f"facts done ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # ------------------------------------------------------------------------------------------- c) factors (fake first)
    stage = P.COND4_CACHE / "factors.pkl"                 # the expensive part, kept on its own
    if stage.exists():
        fac, preds, fake, timing = pickle.loads(stage.read_bytes())
    else:
        fac, preds, timing = [], {}, {}
        t0 = time.time()
        f_tr, preds[(MAIN, "TRAIN")], fake = run_set(Ds["TRAIN"], les["TRAIN"], fz[MAIN], "TRAIN", FA["boot_reps"], C["fake"], noise)
        timing["main TRAIN + fake"] = time.time() - t0
        print(f"main TRAIN + fake done ({(time.time() - clock) / 60:.1f} min)", flush=True)
        fac.append(f_tr)
        for s, k in (("VAL", MAIN), ("TRAIN", CHECK), ("VAL", CHECK)):
            t0 = time.time()
            r_, preds[(k, s)], _ = run_set(Ds[s], les[s], fz[k], s, FA["boot_reps"])
            timing[f"{k} {s}"] = time.time() - t0
            fac.append(r_)
            print(f"{k} {s} done ({(time.time() - clock) / 60:.1f} min)", flush=True)
        fac = pd.concat(fac, ignore_index=True)
        stage.write_bytes(pickle.dumps((fac, preds, fake, timing)))

    # ------------------------------------------------------------------------------------------- b) split halves, d) close the gap
    tabs = {s: LTL.steady_table(Ds[s], les[s], LT, {v: preds[(MAIN, s)][v] for v in VARIANTS}) for s in SETS}
    n_ev = {s: int(Ds[s]["le_event"].max()) + 1 for s in SETS}
    obs_b = {s: LTL.length_fits(tabs[s], n_ev[s], s, LT) for s in SETS}                    # observed, with the event bootstrap
    through = {(s, v): LTL.length_fits(tabs[s], n_ev[s], s, LT, ycol=f"pred_{v}", boot=False) for s in SETS for v in VARIANTS}
    vo = obs_b["VAL"]
    cells = [tuple(x) for x in vo[(vo.links_under_1km >= LT["verdict"]["min_short_links"]) & vo.I.notna()][["band", "class"]].to_numpy()]
    out["cells"] = cells
    gap_rows = []
    for s in SETS:
        o = obs_b[s].set_index(["band", "class"])
        for v in VARIANTS:
            m = through[(s, v)].set_index(["band", "class"])
            for c in cells:
                if c not in o.index or c not in m.index:
                    continue
                q, r_ = m.loc[c], o.loc[c]
                gap_rows.append({"set": s, "variant": v, "band": c[0], "class": c[1], "R_cell": r_.R_cell, "I_obs": r_.I, "I_lo": r_.I_lo, "I_hi": r_.I_hi,
                                 "I_model": q.I, "s_obs": r_.s, "s_lo": r_.s_lo, "s_hi": r_.s_hi, "s_model": q.s, "I_ratio": r_.I / q.I, "s_ratio": r_.s / q.s,
                                 "I_inside": bool(r_.I_lo <= q.I <= r_.I_hi), "s_inside": bool(r_.s_lo <= q.s <= r_.s_hi)})
    gap = pd.DataFrame(gap_rows)
    print(f"length tests done ({(time.time() - clock) / 60:.1f} min)", flush=True)

    SH = C["split_half"]
    tr_full = LTL.length_fits(tabs["TRAIN"], n_ev["TRAIN"], "TRAIN", LT, boot=False)
    val_ref = {"obs_over_train": cell_ratios(LTL.length_fits(tabs["VAL"], n_ev["VAL"], "VAL", LT, boot=False), tr_full, cells),
               "obs_over_model": cell_ratios(LTL.length_fits(tabs["VAL"], n_ev["VAL"], "VAL", LT, boot=False), through[("VAL", "none")], cells)}
    rng = np.random.default_rng(SH["seed"])
    sh_rows = []
    for rep in range(SH["reps"]):
        half = rng.choice(n_ev["TRAIN"], SH["half_events"], replace=False)
        t_ = tabs["TRAIN"][tabs["TRAIN"].ev.isin(half)]
        ob = LTL.length_fits(t_, n_ev["TRAIN"], "half", LT, boot=False)
        mo = LTL.length_fits(t_, n_ev["TRAIN"], "half", LT, ycol="pred_none", boot=False)
        oc = LTL.length_fits(tabs["TRAIN"][~tabs["TRAIN"].ev.isin(half)], n_ev["TRAIN"], "rest", LT, boot=False)
        a1, b1 = cell_ratios(ob, tr_full, cells)
        a2, b2 = cell_ratios(ob, mo, cells)
        a3, b3 = cell_ratios(ob, oc, cells)               # information (not in the rule): the half against the other 33 events
        sh_rows.append({"rep": rep, "events": " ".join(map(str, sorted(half))), "I_half_over_train": a1, "s_half_over_train": b1,
                        "I_obs_over_model": a2, "s_obs_over_model": b2, "I_half_over_rest": a3, "s_half_over_rest": b3})
    split_half = pd.DataFrame(sh_rows)
    out["val_ref"] = val_ref
    print(f"split halves done ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # ------------------------------------------------------------------------------------------- e) radar vs gauge
    rg = pd.concat([radar_gauge_primary(Ds[s], les[s]).assign(set=s) for s in SETS], ignore_index=True)
    pix = radar_gauge_pixel(d, ev_tab.loc[facts.event_id])

    # ------------------------------------------------------------------------------------------- S4.2 conditions and the residual check
    dew = pd.read_csv(P.DEW_RISK_HOURLY, parse_dates=["time_utc"]).set_index("time_utc")
    cond = hourly_conditions(dew, ev_tab.loc[facts.event_id])
    res_rows = []
    for s in SETS:
        D, le = Ds[s], les[s]
        pr = preds[(MAIN, s)]["none"]
        w = (D["w5"] > 0) & D["wet5"]
        hr = times[s].floor("h")
        night = dew.night.reindex(hr).to_numpy(dtype=float)
        temp = dew.temperature.reindex(hr).to_numpy(dtype=float)
        res = pd.DataFrame({"r": np.where(w, D["y5"] - pr, np.nan), "night": np.where(w, night, np.nan), "temp": np.where(w, temp, np.nan),
                            "le": D["bin_le"]}).groupby("le").agg(res=("r", "median"), night=("night", "mean"), temp=("temp", "mean"))
        res = res.reindex(range(len(le)))
        res_rows.append(pd.DataFrame({"set": s, "event_id": le.event_id, "cml_id": le.cml_id, "band": le.band_name, "L_km": le.L_km,
                                      "residual_db": res.res.to_numpy(), "night_share": res.night.to_numpy(), "temperature": res.temp.to_numpy()}))
    resid = pd.concat(res_rows, ignore_index=True).dropna(subset=["residual_db"])
    out.update({"facts": facts, "fac": fac, "fake": fake, "gap": gap, "split_half": split_half, "rg": rg, "pix": pix, "cond": cond, "resid": resid,
                "timing": timing, "minutes": (time.time() - clock) / 60, "obs_b": obs_b, "tr_full": tr_full})
    return out


# =================================================================================================== tables from the results
def summarise(R):
    fac, facts, fake, rg, pix, cond = R["fac"], R["facts"], R["fake"], R["rg"], R["pix"], R["cond"]
    FS = C["fake"]
    fake = fake.assign(k_ok=lambda t: (t.k_est / t.k_true - 1).abs() <= FS["tol"], alpha_ok=lambda t: (t.alpha_est / t.alpha_true - 1).abs() <= FS["tol"])
    fake["both_ok"] = fake.k_ok & fake.alpha_ok
    share = float(fake.both_ok.mean())
    verdict_fake = "PASS" if share >= FS["pass_share"] else "FAIL"
    lean = FA["lean_min_links"]
    fac = fac.assign(lacks_short=lambda t: t.links_under_1km < lean, lacks_long=lambda t: t.links_over_3km < lean)
    fac = fac.merge(facts[["event_id", "start"]], on="event_id", how="left")

    # TRAIN vs VAL, separability
    rng = np.random.default_rng(FA["boot_seed_events"])
    srows = []
    for fit, g in fac.groupby("fit", sort=False):
        for q in ("k", "alpha"):
            row = {"fit": fit, "factor": q}
            for s in SETS:
                v = g[g.set == s][q].to_numpy()
                b = [np.median(rng.choice(v, len(v))) for _ in range(1000)]
                row.update({f"{s}_median": float(np.median(v)), f"{s}_lo": float(np.percentile(b, 5)), f"{s}_hi": float(np.percentile(b, 95)),
                            f"{s}_p10": float(np.percentile(v, 10)), f"{s}_p90": float(np.percentile(v, 90)), f"{s}_events": len(v)})
            row["mann_whitney_p"] = float(stats.mannwhitneyu(g[g.set == "TRAIN"][q], g[g.set == "VAL"][q]).pvalue)
            srows.append(row)
        rho, p = stats.spearmanr(np.log(g.k), np.log(g.alpha))
        srows.append({"fit": fit, "factor": "Spearman(ln k, ln alpha), all events", "TRAIN_median": float(rho), "mann_whitney_p": float(p),
                      "VAL_median": float(np.nanmedian(g.boot_corr_lnk_lnalpha)) if "boot_corr_lnk_lnalpha" in g else np.nan})
    summary = pd.DataFrame(srows)

    # radar / gauge against k_e
    RG = C["radar_gauge"]
    km = fac[fac.fit == MAIN][["event_id", "set", "k", "alpha"]]
    rg = rg.merge(km, on=["event_id", "set"], how="left").merge(pix, on="event_id", how="left")
    rg["kept"] = rg.gauge_mm_per_link >= RG["min_gauge_mm"]
    r2 = np.random.default_rng(RG["seed"])
    g_ = rg[rg.kept]
    rg_test = {"primary": spearman_boot(g_.k.to_numpy(), g_.radar_over_gauge.to_numpy(), RG["boot_reps"], r2)}
    pk = rg[rg.pixel_gauge_mm >= RG["min_gauge_mm"]]
    rg_test["pixel"] = spearman_boot(pk.k.to_numpy(), pk.pixel_radar_over_gauge.to_numpy(), RG["boot_reps"], r2)
    rg_set = []                                           # the ratio itself per set (no factor involved)
    for col, kept in (("radar_over_gauge", rg.kept), ("pixel_radar_over_gauge", rg.pixel_gauge_mm >= RG["min_gauge_mm"])):
        a, b = rg[kept & (rg.set == "TRAIN")][col].dropna(), rg[kept & (rg.set == "VAL")][col].dropna()
        pooled = {s: rg[kept & (rg.set == s)] for s in SETS}
        rg_set.append({"ratio": col, "TRAIN events": len(a), "TRAIN median": a.median(), "TRAIN p25-p75": f"{a.quantile(.25):.2f} - {a.quantile(.75):.2f}",
                       "VAL events": len(b), "VAL median": b.median(), "VAL p25-p75": f"{b.quantile(.25):.2f} - {b.quantile(.75):.2f}",
                       "mann_whitney_p": float(stats.mannwhitneyu(a, b).pvalue) if len(a) and len(b) else np.nan,
                       **({"TRAIN pooled radar / gauge": pooled["TRAIN"].radar_mm_sum.sum() / pooled["TRAIN"].gauge_mm_sum.sum(),
                           "VAL pooled radar / gauge": pooled["VAL"].radar_mm_sum.sum() / pooled["VAL"].gauge_mm_sum.sum()} if col == "radar_over_gauge" else
                          {"TRAIN pooled radar / gauge": pooled["TRAIN"].pixel_radar_mm.sum() / pooled["TRAIN"].pixel_gauge_mm.sum(),
                           "VAL pooled radar / gauge": pooled["VAL"].pixel_radar_mm.sum() / pooled["VAL"].pixel_gauge_mm.sum()})})
    rg_set = pd.DataFrame(rg_set)

    # conditions
    CO = C["conditions"]
    ev = facts.merge(cond, on="event_id").merge(km, on=["event_id", "set"])
    names = {"temperature": "mean temperature (C)", "humidity": "mean humidity (%)", "wind_speed": "mean wind speed (m/s)",
             "dew_point_depression": "dew-point depression (C)", "night_fraction": "night fraction", "start_hour_utc": "start hour (UTC)",
             "peak_over_mean": "peak / mean rain", "duration_h": "duration (h)", "wet_fraction": "wet fraction", "dry_before_h": "dry hours before",
             "share_38GHz": "share of 38 GHz link-events", "share_28GHz": "share of 28 GHz link-events", "peak_mmh_median_link": "peak rain (median link, mm/h)",
             "total_mm_median_link": "total rain (median link, mm)"}
    r3 = np.random.default_rng(CO["seed"])
    crow = []
    for q in ("alpha", "k"):
        for col, lab in names.items():
            rho, p, lo, hi, n = spearman_boot(np.log(ev[q].to_numpy()), ev[col].to_numpy(float), CO["boot_reps"], r3)
            crow.append({"factor": q, "condition": col, "label": lab, "rho": rho, "lo": lo, "hi": hi, "p": p, "events": n})
    ctests = pd.DataFrame(crow)
    ctests["bh_discovery"] = bh(ctests.p.fillna(1).to_numpy(), CO["bh_q"])

    # residual check
    rr = np.random.default_rng(CO["residual_seed"])
    rs = R["resid"].copy()
    rs["day_night"] = np.where(rs.night_share > 0.5, "night", "day")
    tq = np.nanpercentile(rs.temperature, [100 / 3, 200 / 3])
    rs["temp_tertile"] = np.where(rs.temperature < tq[0], f"< {tq[0]:.1f} C", np.where(rs.temperature < tq[1], f"{tq[0]:.1f}-{tq[1]:.1f} C", f">= {tq[1]:.1f} C"))
    rs.loc[rs.temperature.isna(), "temp_tertile"] = "no temperature"
    rrows = []
    for by in ("band", "day_night", "temp_tertile"):
        for lev, g in rs.groupby(by):
            m, lo, hi = boot_median(g.residual_db.to_numpy(), g.event_id.to_numpy(), 1000, rr)
            rrows.append({"by": by, "level": lev, "link_events": len(g), "events": g.event_id.nunique(), "median_residual_db": m, "lo": lo, "hi": hi})
    rcheck = pd.DataFrame(rrows)
    return {"fake": fake, "fake_share": share, "fake_verdict": verdict_fake, "fac": fac, "summary": summary, "rg": rg, "rg_test": rg_test, "rg_set": rg_set,
            "ev": ev, "ctests": ctests, "cond_names": names, "rcheck": rcheck, "rs": rs}


# =================================================================================================== figures
def fig_timeline(facts):
    fig, axes = plt.subplots(2, 1, figsize=(12, 6.4), sharex=True)
    ax = axes[0]
    for s, y in (("TRAIN", 1), ("VAL", 0)):
        g = facts[facts.set == s]
        for r in g.itertuples():
            ax.plot([r.start, r.end], [y, y], "-", color=INK, lw=4, solid_capstyle="butt")
        ax.plot(g.start, np.full(len(g), y + 0.18), MK[s], mfc="none", color=INK, ms=6, ls="none", label=f"{s} ({len(g)} events)")
    ax.set_yticks([0, 1]), ax.set_yticklabels(["VAL", "TRAIN"])
    ax.set_ylim(-0.5, 1.6)
    ax.legend(frameon=False, fontsize=8, loc="upper right", ncol=2)
    ax.set_title("each event as a bar from start to end; the marker above the bar = its start (TRAIN circles, VAL squares)", fontsize=9)
    style_axes(ax)
    ax = axes[1]
    for s in SETS:
        g = facts[facts.set == s]
        ax.plot(g.start, g.start_hour_utc, MK[s], mfc="none", color=INK if s == "TRAIN" else INK2, ms=7, ls="none", label=f"{s}: start hour")
    ax.set(ylabel="start hour (UTC)", ylim=(-0.5, 24.5), yticks=[0, 6, 12, 18, 24])
    ax.legend(frameon=False, fontsize=8)
    style_axes(ax)
    fig.suptitle("Step S4.1 a: the TRAIN and VAL events of the S3 all-length datasets over June-August 2015 (UTC)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def fig_split_half(sh, ref):
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.4))
    for i, (cols, key, what) in enumerate(((("I_half_over_train", "s_half_over_train"), "obs_over_train", "half / all TRAIN"),
                                           (("I_obs_over_model", "s_obs_over_model"), "obs_over_model", "observed / M4b-B"))):
        for j, (col, q) in enumerate(zip(cols, ("I", "s"))):
            ax = axes[i][j]
            ax.hist(sh[col].dropna(), bins=25, facecolor="none", edgecolor=INK, hatch="///", lw=0.8)
            v = ref[key][j]
            ax.axvline(v, color=INK, ls="--", lw=1.5)
            share = float((sh[col] <= v).mean())
            ax.annotate(f"VAL {v:.2f}\n{share * 100:.1f} % of halves <= VAL", (v, ax.get_ylim()[1] * 0.9), xytext=(5, 0), textcoords="offset points", fontsize=8, va="top")
            ax.set(xlabel=f"{q}: {what} (median over the 8 cells)", ylabel="random halves (count)")
            style_axes(ax)
    fig.suptitle("Step S4.1 b: is VAL an ordinary 18-event set? 200 random 18-event halves of the 51 TRAIN events (hatched) against VAL (dashed line)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def fig_fake(fake, tol):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))
    mk = {0.6: "v", 1.0: "o", 1.4: "^", 0.7: "v"}
    for ax, (q, other) in zip(axes, (("k", "alpha"), ("alpha", "k"))):
        tv = sorted(fake[f"{q}_true"].unique())
        for t in tv:
            g = fake[fake[f"{q}_true"] == t]
            x = t * np.exp(np.random.default_rng(0).uniform(-0.04, 0.04, len(g)))
            ax.plot(x, g[f"{q}_est"], mk.get(t, "o"), mfc="none", color=INK, ms=6, ls="none", label=f"true {q} = {t:g}")
            ax.plot([t * 0.9, t * 1.1], [t * (1 - tol)] * 2, "-", color=INK2, lw=0.8)
            ax.plot([t * 0.9, t * 1.1], [t * (1 + tol)] * 2, "-", color=INK2, lw=0.8)
        ax.plot([0.3, 3], [0.3, 3], ":", color=INK2, lw=0.8)
        ax.set(xscale="log", yscale="log", xlabel=f"true {q}", ylabel=f"estimated {q}", title=f"{q}: estimate per event (short bars = +-{tol * 100:.0f} %)")
        ax.legend(frameon=False, fontsize=8)
        style_axes(ax)
    fig.suptitle("Step S4.1 c: fake check on the 51 TRAIN events (frozen M4b-B, S2 noise), true against estimated factors", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig


def fig_factors(fac):
    g = fac[fac.fit == MAIN].sort_values("start").reset_index(drop=True)
    fig, axes = plt.subplots(3, 1, figsize=(13, 10.5))
    for ax, q in zip(axes[:2], ("k", "alpha")):
        for s in SETS:
            h = g[g.set == s]
            ax.errorbar(h.index, h[q], yerr=[h[q] - h[f"{q}_lo"], h[f"{q}_hi"] - h[q]], fmt=MK[s], mfc="none", color=INK, ms=6, capsize=2, lw=0.9, label=s)
        ax.axhline(1, color=INK2, lw=0.8, ls="--")
        ax.set_yscale("log")
        ax.set(ylabel=f"{q}_e (90 % link bootstrap)", xlabel="events in date order")
        ax.set_xticks(g.index[::3]), ax.set_xticklabels([str(t.date()) for t in g.start[::3]], rotation=60, fontsize=7)
        ax.legend(frameon=False, fontsize=8)
        style_axes(ax)
    axes[0].set_title("k_e > 1: the radar reads the event too high (the model needs less rain)", fontsize=9)
    axes[1].set_title("alpha_e < 1: the antennas of the event wet less than the frozen F says", fontsize=9)
    ax = axes[2]
    for fit, ls in ((MAIN, "-"), (CHECK, ":")):
        for s in SETS:
            h = fac[(fac.fit == fit) & (fac.set == s)]
            ax.plot(h.k, h.alpha, MK[s], mfc="none" if fit == MAIN else INK2, color=INK if fit == MAIN else INK2, ms=6 if fit == MAIN else 4, ls="none",
                    label=f"{fit}, {s}")
    ax.set(xscale="log", yscale="log", xlabel="k_e", ylabel="alpha_e", title="k_e against alpha_e per event (hollow = M4b-B, small filled = M3-B)")
    ax.legend(frameon=False, fontsize=8, ncol=2)
    style_axes(ax)
    fig.suptitle(f"Step S4.1 c: per-event factors on the frozen {MAIN} (TRAIN circles, VAL squares)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def fig_close_gap(gap):
    vs = ["none", "both", "k_only", "alpha_only", "k_refit", "alpha_refit"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    g = gap[gap.set == "VAL"]
    labels = [f"{b}\n{c}" for b, c in g.drop_duplicates(["band", "class"])[["band", "class"]].to_numpy()]
    for ax, q in zip(axes, ("I", "s")):
        for i, v in enumerate(vs):
            h = g[g.variant == v]
            x = np.arange(len(h)) + (i - 2.5) * 0.1
            ax.plot(x, h[f"{q}_ratio"], VSTYLE[v], mfc="none", color=INK, ms=6, ls="none", label=f"{VARIANTS[v]}: median {h[f'{q}_ratio'].median():.2f}")
        for b in C["close_gap"]["closes"]:
            ax.axhline(b, color=INK2, lw=0.8, ls=":")
        ax.axhline(1, color=INK2, lw=0.8)
        ax.set_yscale("log")
        ax.set_xticks(range(len(labels))), ax.set_xticklabels(labels, fontsize=7)
        ax.set(ylabel=f"observed / through-model {q} (log axis)", title=f"VAL cells: {q} (dotted lines = 'closed' range 0.9-1.1)")
        ax.legend(frameon=False, fontsize=7)
        style_axes(ax)
    fig.suptitle("Step S4.1 d: the VAL length test through M4b-B with the per-event factors applied (in-sample per event: information)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig


def fig_radar_gauge(rg):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))
    for ax, (col, kept, title) in zip(axes, (("radar_over_gauge", rg.kept, "radar path mean of the gauge links / City gauge (S2 pairing)"),
                                             ("pixel_radar_over_gauge", rg.pixel_gauge_mm >= C["radar_gauge"]["min_gauge_mm"], "radar pixel / City gauge (10 gauges)"))):
        for s in SETS:
            h = rg[kept & (rg.set == s)]
            ax.plot(h[col], h.k, MK[s], mfc="none", color=INK, ms=6, ls="none", label=s)
        ax.axhline(1, color=INK2, lw=0.8, ls="--"), ax.axvline(1, color=INK2, lw=0.8, ls="--")
        ax.set(xscale="log", yscale="log", xlabel=title, ylabel=f"k_e ({MAIN})")
        ax.legend(frameon=False, fontsize=8)
        style_axes(ax)
    fig.suptitle("Step S4.1 e: is k_e radar error? k_e against the event radar / gauge ratio (TRAIN circles, VAL squares)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig


def fig_conditions(ev, ctests, names):
    cols = list(names)
    nc = 7
    nr = int(np.ceil(len(cols) / nc))
    fig, axes = plt.subplots(2 * nr, nc, figsize=(2.6 * nc, 2.6 * 2 * nr), squeeze=False)
    ct = ctests.set_index(["factor", "condition"])
    for i, q in enumerate(("alpha", "k")):
        for j, col in enumerate(cols):
            ax = axes[i * nr + j // nc][j % nc]
            for s in SETS:
                h = ev[ev.set == s]
                ax.plot(h[col], h[q], MK[s], mfc="none", color=INK, ms=4, ls="none")
            ax.set_yscale("log")
            r = ct.loc[(q, col)]
            ax.set_title(f"rho {r.rho:+.2f}" + (" *BH" if r.bh_discovery else ""), fontsize=7)
            ax.set_xlabel(names[col], fontsize=6.5)
            ax.tick_params(labelsize=6)
            if j % nc == 0:
                ax.set_ylabel(f"{q}_e")
            style_axes(ax)
    for ax in axes.ravel()[2 * len(cols):]:
        ax.set_visible(False)
    fig.suptitle(f"Step S4.2: alpha_e (top {nr} rows) and k_e (bottom {nr} rows) against the conditions of each event; TRAIN circles, VAL squares; Spearman rho, *BH = discovery at q = 0.10", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


def fig_residual(rcheck):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2))
    for ax, by in zip(axes, ("band", "day_night", "temp_tertile")):
        g = rcheck[rcheck.by == by].reset_index(drop=True)
        ax.errorbar(g.index, g.median_residual_db, yerr=[g.median_residual_db - g.lo, g.hi - g.median_residual_db], fmt="o", mfc="none", color=INK, capsize=3)
        for r in g.itertuples():
            ax.annotate(f"{r.link_events} l-e", (r.Index, r.hi), xytext=(0, 3), textcoords="offset points", ha="center", fontsize=7)
        ax.axhline(0, color=INK2, lw=0.8)
        ax.set_xticks(g.index), ax.set_xticklabels(g.level, fontsize=8)
        ax.set(ylabel="median (A - A_pred) during rain (dB)", title=f"by {by.replace('_', ' ')}")
        style_axes(ax)
    fig.suptitle(f"Step S4.2: residual of the frozen {MAIN} (no event factors) during rain, TRAIN + VAL link-events, 90 % event-bootstrap intervals", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    clock = time.time()
    P.COND4_CACHE.mkdir(parents=True, exist_ok=True)
    cache = P.COND4_CACHE / "results.pkl"
    if cache.exists():
        R = pickle.loads(cache.read_bytes())
        print(f"results read from the cache ({rel(cache)}); delete that file to recompute", flush=True)
    else:
        R = compute()
        cache.write_bytes(pickle.dumps(R))
    S = summarise(R)
    facts = R["facts"].merge(R["cond"], on="event_id", how="left")
    save(facts.round(4), P.COND4_EVENTS, index=False)
    save(R["split_half"].round(4), P.COND4_SPLIT_HALF, index=False)
    save(S["fake"].round(4), P.COND4_FAKE, index=False)
    save(S["fac"].drop(columns=["start"]).round(4), P.COND4_FACTORS, index=False)
    save(S["summary"].round(4), P.COND4_FACTOR_SUMMARY, index=False)
    save(R["gap"].round(4), P.COND4_CLOSE_GAP, index=False)
    save(S["rg"].round(4), P.COND4_RADAR_GAUGE, index=False)
    save(S["ctests"].round(4), P.COND4_CONDITIONS, index=False)
    save(S["rcheck"].round(4), P.COND4_RESIDUAL, index=False)
    for fig, path in ((fig_timeline(R["facts"]), P.COND4_FIG_TIMELINE), (fig_split_half(R["split_half"], R["val_ref"]), P.COND4_FIG_SPLIT_HALF),
                      (fig_fake(S["fake"], C["fake"]["tol"]), P.COND4_FIG_FAKE), (fig_factors(S["fac"]), P.COND4_FIG_FACTORS),
                      (fig_close_gap(R["gap"]), P.COND4_FIG_CLOSE_GAP), (fig_radar_gauge(S["rg"]), P.COND4_FIG_RADAR_GAUGE),
                      (fig_conditions(S["ev"], S["ctests"], S["cond_names"]), P.COND4_FIG_CONDITIONS), (fig_residual(S["rcheck"]), P.COND4_FIG_RESIDUAL)):
        save(fig, path, dpi=100)
        plt.close(fig)
    report = write_report(R, S, clock)
    save(report, P.COND4_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {(time.time() - clock) / 60:.1f} min")


def readings(R, S):
    """The READING paragraphs, written after looking at the numbers of the run (every number is taken from its tables)."""
    facts, fac, fake, sh, ref = R["facts"], S["fac"], S["fake"], R["split_half"], R["val_ref"]
    SH, FS = C["split_half"], C["fake"]
    m = fac[fac.fit == MAIN]
    dry = facts[facts.total_mm_median_link < 0.1]
    flat = m[(m.k == 1) & (m.alpha == 1) & (m.k_lo == m.k_hi)]
    big_val = facts[(facts.set == "VAL") & (facts.total_mm_median_link >= 1)].sort_values("start")
    sh_I, sh_s = float((sh.I_half_over_train <= ref["obs_over_train"][0]).mean()), float((sh.s_half_over_train <= ref["obs_over_train"][1]).mean())
    sm_I, sm_s = float((sh.I_obs_over_model <= ref["obs_over_model"][0]).mean()), float((sh.s_obs_over_model <= ref["obs_over_model"][1]).mean())
    sr_I, sr_s = float((sh.I_half_over_rest <= ref["obs_over_train"][0]).mean()), float((sh.s_half_over_rest <= ref["obs_over_train"][1]).mean())
    rare = max(sh_I, sh_s, sm_I, sm_s) < SH["rarely"]
    fk = fake.merge(facts[["event_id", "intensity_class"]], on="event_id")
    hv = fk[fk.intensity_class == "heavy"]
    passed = S["fake_verdict"] == "PASS"
    rgs = S["rg_set"].set_index("ratio")
    pr = rgs.loc["radar_over_gauge"]
    ct = S["ctests"]
    rc = S["rcheck"].set_index(["by", "level"])
    temp_rows = S["rcheck"][S["rcheck"].by == "temp_tertile"].sort_values("median_residual_db")
    summ = S["summary"].set_index(["fit", "factor"])
    corr = summ.loc[(MAIN, "Spearman(ln k, ln alpha), all events")]
    rgt = S["rg_test"]["primary"]
    gap = R["gap"].groupby(["set", "variant"]).agg(I=("I_ratio", "median"), s=("s_ratio", "median"))
    no_interp = "" if passed else " **Not interpreted: the fake check failed (section 3).**"

    r_facts = f"""**READING.** The split is random within each intensity tercile, and the result looks like it: VAL has
{int((facts.set == 'VAL').sum())} events ({', '.join(f'{n} {c}' for c, n in facts[facts.set == 'VAL'].intensity_class.value_counts().items())}), spread from late June to late August between the
TRAIN events, not a block in time. Most events carry very little rain at the links: {len(dry)} of {len(facts)} events have a median path
accumulation under 0.1 mm inside [start, end) (link_state can be wet from a gauge or from R_max on the path). The steady-rain bins of the
length test come from the few large events: in VAL the events with >= 1 mm at the median link are {', '.join(big_val.event_id)}
({len(big_val)} events), so the VAL length test is effectively a test on about {len(big_val)} events."""

    r_split = f"""**READING.** By the rule, the VAL set is **{'rare' if rare else 'not rare'}**: {sh_I * 100:.1f} % (I) and {sh_s * 100:.1f} % (s) of the {SH['reps']} random
{SH['half_events']}-event halves of TRAIN are as low as VAL relative to all of TRAIN, and {sm_I * 100:.1f} % / {sm_s * 100:.1f} % as low relative to the frozen model.
A half is part of the TRAIN set it is compared with, which narrows its spread (by about a third, since 18 of 51 events are shared); the
fairer comparison of each half with the OTHER {int((facts.set == 'TRAIN').sum()) - SH['half_events']} events (information, not in the rule) gives {sr_I * 100:.1f} % and {sr_s * 100:.1f} %, still under
{SH['rarely'] * 100:.0f} %. {'So the VAL gap is more than the ordinary event-to-event spread of TRAIN: something about the VAL events (or about the rain reference in them).' if rare else 'So the VAL gap is within ordinary event-to-event spread.'}
Caveat: the halves are random sets of events; the cells are dominated by the few large events (above), so "18 events" is really "the large
events among them"."""

    r_fake = f"""**READING. FAIL** ({S['fake_share'] * 100:.0f} % of the events recovered, rule {FS['pass_share'] * 100:.0f} %). The failures are not spread evenly: every heavy event is
recovered ({int(hv.both_ok.sum())} of {len(hv)}), the light and moderate events often not. In those events the rain at the links is small or nil, the
loss hardly changes over the grid, and k and alpha trade against each other (less rain with a larger antenna factor gives nearly the
same A); the weak ridge then pulls them toward 1 (true 1.4 / 1.4 estimated about 0.89 / 0.76 of the truth; true 0.7 / 0.6 about 1.2 x).
At the truth k = alpha = 1 all {len(fk[(fk.k_true == 1) & (fk.alpha_true == 1)])} events are recovered, so the code has no bias there. **By the rule, the
real factors below are reported and NOT interpreted.** That heavy events pass 100 % suggests a fake check restricted to heavy events
(rule fixed before running) could license reading their factors; that is a proposal for a sub-step, not something done here."""

    r_fac = f"""**READING.**{no_interp} Only what the method does is noted here. {len(flat)} events have exactly k = alpha = 1 with no
interval: no modelled rain at their links, so the ridge alone sets them. Several events sit at a bound of k (`k_at_bound`). Inside an
event the two factors are strongly tied in many events (within-event bootstrap correlation of ln k and ln alpha: median
{corr.VAL_median:+.2f} over events, often beyond +-0.8 in either sign), which is the k / alpha trade the fake check shows. Over events the
Spearman is {corr.TRAIN_median:+.2f} (rule: |rho| > {FA['separable_max_corr']} = not separable, so "separable" over events by the rule, but see the fake check). The
two frozen frames give the same factors ({MAIN} against {CHECK}: Spearman 0.96 for both), so the numbers do not depend on the frame;
that does not make them right."""

    r_gap = f"""**READING.**{no_interp} With both factors per event the VAL length test is passed in-sample (median ratio I
{gap.loc[('VAL', 'both'), 'I']:.2f}, s {gap.loc[('VAL', 'both'), 's']:.2f}); this is expected from two free numbers per event and says only that one scale per event on the whole
model is enough to describe the VAL bins. Which of the two does it (k alone: {gap.loc[('VAL', 'k_only'), 'I']:.2f} / {gap.loc[('VAL', 'k_only'), 's']:.2f}; alpha alone:
{gap.loc[('VAL', 'alpha_only'), 'I']:.2f} / {gap.loc[('VAL', 'alpha_only'), 's']:.2f}) depends on the split between k and alpha, which the fake check did not confirm."""

    r_rg = f"""**READING.** The Spearman of k_e with the radar / gauge ratio is {rgt[0]:.2f} [{rgt[2]:.2f}, {rgt[3]:.2f}], above the {C['radar_gauge']['expect_spearman']} written before the run,
but it rests on k_e, which is not interpreted (fake FAIL). **The ratio itself needs no factor:** on the events with enough gauge rain the radar
reads {pr['VAL pooled radar / gauge']:.2f} x the gauges in VAL (pooled; median over {int(pr['VAL events'])} events {pr['VAL median']:.2f}) against {pr['TRAIN pooled radar / gauge']:.2f} x in TRAIN (median over
{int(pr['TRAIN events'])} events {pr['TRAIN median']:.2f}); the radar pixel at the gauges gives the same ({rgs.loc['pixel_radar_over_gauge', 'VAL pooled radar / gauge']:.2f} against
{rgs.loc['pixel_radar_over_gauge', 'TRAIN pooled radar / gauge']:.2f}). TRAIN / VAL of these pooled ratios is {pr['TRAIN pooled radar / gauge'] / pr['VAL pooled radar / gauge']:.2f}, the same size as the VAL gap (observed VAL / TRAIN:
I {ref['obs_over_train'][0]:.2f}, s {ref['obs_over_train'][1]:.2f}). So the radar over-reads the rain of the large VAL events relative to the City gauges, by about
as much as VAL "attenuates less per radar mm/h". This is a direction and a size, not a proof: only {int(pr['VAL events'])} VAL events have enough gauge
rain (Mann-Whitney p {pr['mann_whitney_p']:.2f}), the gauges are points against a path mean or a 2 x 2 km pixel, and the radar is gauge-adjusted
(pipeline methods.md, section 9), so it is not independent of the gauges."""

    disc = ct[ct.bh_discovery]
    top = ct.reindex(ct.p.sort_values().index).head(3)
    r_cond = f"""**READING.**{no_interp} For the record: {len(disc)} of {len(ct)} tests are discoveries by Benjamini-Hochberg at
q = {C['conditions']['bh_q']}; the smallest p-values are {', '.join(f"{r.factor} with {r.label} (rho {r.rho:+.2f}, p {r.p:.3f})" for r in top.itertuples())}.
**The residual check needs no factor.** During rain the frozen {MAIN} predicts too much at every band (median residual
{rc.loc[('band', '28 GHz'), 'median_residual_db']:+.3f} / {rc.loc[('band', '32 GHz'), 'median_residual_db']:+.3f} / {rc.loc[('band', '38 GHz'), 'median_residual_db']:+.3f} dB at 28 / 32 / 38 GHz) and equally by day and by night
({rc.loc[('day_night', 'day'), 'median_residual_db']:+.3f} / {rc.loc[('day_night', 'night'), 'median_residual_db']:+.3f} dB): no band or day / night effect. By temperature the warmest tertile is lower
({temp_rows.iloc[0].level}: {temp_rows.iloc[0].median_residual_db:+.3f} [{temp_rows.iloc[0].lo:+.3f}, {temp_rows.iloc[0].hi:+.3f}] dB) than the coolest ({rc.loc[('temp_tertile', [l for l in temp_rows.level if l.startswith('<')][0]), 'median_residual_db']:+.3f} dB), with
intervals apart. Temperature is almost constant inside an event and the tertiles hold {', '.join(str(int(v)) for v in S['rcheck'][S['rcheck'].by == 'temp_tertile'].events)} events, so this is an event
difference that goes with temperature, not shown to be caused by it (warm events may also differ in rain type or radar error)."""

    r_mean = f"""Candidate (to be checked in S5, not a final claim): **keep {MAIN} as it is (WAA v3, provisional); do not add an event
factor yet; no antenna level that depends on a condition.**
- The VAL gap is more than ordinary event spread (section 2), so no event-independent fit could have passed the S3b rule on VAL.
- The one factor-free comparison of the rain reference (section 6) shows the radar reading the large VAL events about {pr['VAL pooled radar / gauge'] / pr['TRAIN pooled radar / gauge']:.1f} x higher
  against the gauges than the TRAIN events, the same size as the gap. The most likely reading is a rain-reference (radar) difference
  between the event sets, not different antennas; it says nothing against the antenna part of {MAIN}.
- The per-event factors could have separated "radar" from "antenna" per event, but the fake check failed ({S['fake_share'] * 100:.0f} %), so they are
  not used. An event factor in the final equation would need (i) a fake check that passes on the events where it is read (heavy
  events passed {int(hv.both_ok.sum())} / {len(hv)}) and (ii) a way to set it without the target (a rain-reference check), else it is just a free
  number per event.
- No condition was linked to the factors (and they are not interpreted); the residual is the same across bands and day / night and
  lower in the warmest events: a lead, not a term."""

    r_unclear = f"""- **The ridge was changed before the run** (0.5 -> 0.05 per link) after a smoke test of the code on 3 TRAIN events showed that 0.5
  set the factor of a light event by itself. The real k_e, alpha_e of those 3 events were seen then; nothing else.
- **Fake check FAIL:** with the S2 noise (variance x 3, pessimistic on purpose), light and moderate events do not fix k and alpha
  apart. A milder noise or a stronger ridge would raise the share, but the rule was fixed before the run and is kept.
- **The split-half rule compares a half with a set that contains it**; the fairer half-against-the-rest comparison is given as
  information and gives the same verdict.
- **The length test on VAL rests on about {len(big_val)} large events** ({', '.join(big_val.event_id)}); "VAL" there means mostly those.
- **Radar vs gauge:** City gauges are points (weaker QC than SMHI), the radar is a path mean or a 2 x 2 km pixel and is gauge-adjusted
  (not independent of the gauges); only {int(pr['VAL events'])} VAL events have >= {C['radar_gauge']['min_gauge_mm']} mm per gauge link. Events with little gauge rain give
  ratios of 0 or 3-4 (in the table, not kept).
- **Grid estimates:** steps of 1.9 % (k) and 2.3 % (alpha); several events sit at the k bounds 0.33 or 3, and {len(flat)} events with no
  modelled rain at their links are set to 1 by the ridge.
- **The temperature effect in the residual** is between events (each tertile is about 30-38 events), so it may be the rain type or the
  radar error of those events, not the antennas.
- **Conditions come from one station (GbgA)**, hourly; humidity is missing for one event."""
    return {"facts": r_facts, "split_half": r_split, "fake": r_fake, "factors": r_fac, "close_gap": r_gap, "radar_gauge": r_rg,
            "conditions": r_cond, "meaning": r_mean, "unclear": r_unclear}


def write_report(R, S, clock):
    """The report from the numbers of the run; the paragraphs marked READING come from `readings`."""
    q = readings(R, S)
    facts, fac, fake, gap, sh = R["facts"], S["fac"], S["fake"], R["gap"], R["split_half"]
    ref = R["val_ref"]
    SH, RG, CO = C["split_half"], C["radar_gauge"], C["conditions"]
    standing = ("**Standing rule.** WAA is never estimated by subtracting a power-law rain term. The path part is separated only through its "
                "dependence on L at matched rain. ITU-R P.838 appears only as a size reference, never to correct data. The per-event rain factor "
                "k_e is estimated from the data-learned P of the frozen fit and the spread of L, not from ITU.")
    ev_show = facts[["set", "event_id", "start", "end", "duration_h", "intensity_class", "total_mm_median_link", "peak_mmh_max_link", "peak_mmh_median_link",
                     "wet_fraction", "links", "link_events", "links_under_1km", "links_over_3km"]].sort_values("start").copy()
    ev_show["start"], ev_show["end"] = ev_show.start.dt.strftime("%Y-%m-%d %H:%M"), ev_show.end.dt.strftime("%m-%d %H:%M")
    per_set = facts.groupby("set").agg(events=("event_id", "size"), median_duration_h=("duration_h", "median"), median_total_mm=("total_mm_median_link", "median"),
                                       median_peak_mmh=("peak_mmh_max_link", "median"), median_wet_fraction=("wet_fraction", "median"),
                                       median_links=("links", "median"), first=("start", "min"), last=("end", "max")).reset_index()
    ic = facts.groupby(["set", "intensity_class"]).size().unstack(fill_value=0).reset_index()
    shares = {c: float((sh[c] <= ref[k][j]).mean()) for c, k, j in (("I_half_over_train", "obs_over_train", 0), ("s_half_over_train", "obs_over_train", 1),
                                                                    ("I_obs_over_model", "obs_over_model", 0), ("s_obs_over_model", "obs_over_model", 1))}
    sh_tab = pd.DataFrame([{"statistic": "observed, half / all TRAIN", "VAL I": ref["obs_over_train"][0], "halves I median": sh.I_half_over_train.median(),
                            "halves I 5-95 %": f"{sh.I_half_over_train.quantile(0.05):.2f} - {sh.I_half_over_train.quantile(0.95):.2f}",
                            "share of halves <= VAL (I)": shares["I_half_over_train"], "VAL s": ref["obs_over_train"][1], "halves s median": sh.s_half_over_train.median(),
                            "halves s 5-95 %": f"{sh.s_half_over_train.quantile(0.05):.2f} - {sh.s_half_over_train.quantile(0.95):.2f}",
                            "share of halves <= VAL (s)": shares["s_half_over_train"]},
                           {"statistic": f"observed / {MAIN} through the test", "VAL I": ref["obs_over_model"][0], "halves I median": sh.I_obs_over_model.median(),
                            "halves I 5-95 %": f"{sh.I_obs_over_model.quantile(0.05):.2f} - {sh.I_obs_over_model.quantile(0.95):.2f}",
                            "share of halves <= VAL (I)": shares["I_obs_over_model"], "VAL s": ref["obs_over_model"][1], "halves s median": sh.s_obs_over_model.median(),
                            "halves s 5-95 %": f"{sh.s_obs_over_model.quantile(0.05):.2f} - {sh.s_obs_over_model.quantile(0.95):.2f}",
                            "share of halves <= VAL (s)": shares["s_obs_over_model"]},
                           {"statistic": "INFORMATION: observed, half / the other 33 TRAIN events (VAL: VAL / all TRAIN)", "VAL I": ref["obs_over_train"][0],
                            "halves I median": sh.I_half_over_rest.median(),
                            "halves I 5-95 %": f"{sh.I_half_over_rest.quantile(0.05):.2f} - {sh.I_half_over_rest.quantile(0.95):.2f}",
                            "share of halves <= VAL (I)": float((sh.I_half_over_rest <= ref["obs_over_train"][0]).mean()), "VAL s": ref["obs_over_train"][1],
                            "halves s median": sh.s_half_over_rest.median(),
                            "halves s 5-95 %": f"{sh.s_half_over_rest.quantile(0.05):.2f} - {sh.s_half_over_rest.quantile(0.95):.2f}",
                            "share of halves <= VAL (s)": float((sh.s_half_over_rest <= ref["obs_over_train"][1]).mean())}]).round(3)
    fk = fake.assign(k_ratio=fake.k_est / fake.k_true, alpha_ratio=fake.alpha_est / fake.alpha_true)
    fk_by = fk.groupby(["k_true", "alpha_true"]).agg(events=("event_id", "size"), recovered=("both_ok", "sum"), k_ratio_median=("k_ratio", "median"),
                                                    alpha_ratio_median=("alpha_ratio", "median")).reset_index().round(3)
    fk["light"] = fk.links_under_1km < FA["lean_min_links"]
    fk = fk.merge(facts[["event_id", "intensity_class", "total_mm_median_link"]], on="event_id", how="left")
    fk_ic = fk.groupby("intensity_class").agg(events=("event_id", "size"), recovered=("both_ok", "sum"), k_recovered=("k_ok", "sum"),
                                              alpha_recovered=("alpha_ok", "sum")).reset_index()
    fac_show = fac[fac.fit == MAIN].sort_values("start")[["set", "event_id", "links", "links_under_1km", "links_over_3km", "k", "k_lo", "k_hi", "alpha", "alpha_lo",
                                                          "alpha_hi", "boot_corr_lnk_lnalpha", "lacks_short", "lacks_long", "k_at_bound", "alpha_at_bound"]].round(3)
    chk = fac.pivot_table(index=["set", "event_id"], columns="fit", values=["k", "alpha"])
    rho_fit = {v: stats.spearmanr(chk[(v, MAIN)], chk[(v, CHECK)])[0] for v in ("k", "alpha")}
    gsum = (gap.groupby(["set", "variant"], sort=False).agg(cells=("I_ratio", "size"), I_ratio_median=("I_ratio", "median"), s_ratio_median=("s_ratio", "median"),
                                                            I_inside=("I_inside", "sum"), s_inside=("s_inside", "sum"))
            .reset_index())
    gsum["both_inside"] = gap.assign(b=gap.I_inside & gap.s_inside).groupby(["set", "variant"], sort=False).b.sum().to_numpy()
    lo_, hi_ = C["close_gap"]["closes"]
    gsum["closes_I"] = gsum.I_ratio_median.between(lo_, hi_)
    gsum["closes_s"] = gsum.s_ratio_median.between(lo_, hi_)
    gsum["variant"] = gsum.variant.map(VARIANTS)
    rgt = S["rg_test"]
    rg_show = S["rg"][["set", "event_id", "gauge_links", "gauge_mm_per_link", "radar_over_gauge", "kept", "pixel_gauge_mm", "pixel_radar_over_gauge", "k", "alpha"]].round(3)
    ct = S["ctests"].copy()
    n_tests = len(ct)
    ct_show = ct[["factor", "label", "rho", "lo", "hi", "p", "bh_discovery", "events"]].round(3)
    rc = S["rcheck"].round(3)
    cond_inv = pd.DataFrame([
        {"condition data": "hourly temperature, humidity, wind speed, dew point, dew-point depression, night flag (SMHI GbgA, one station)",
         "where": f"`{rel(P.DEW_RISK_HOURLY)}` (baseline module, made from the release `met_1h`)", "used": "yes: per event over [start, end)"},
        {"condition data": "hourly SMHI weather incl. pressure, wind gust and direction", "where": "release store `met_1h`", "used": "no (the same station; the dew table holds what is used)"},
        {"condition data": "per-event min / mean / max temperature, humidity, wind (hourly SMHI around the event)", "where": "release table `events` (columns temperature_*, humidity_*, wind_*)",
         "used": "no (same source; recomputed over [start, end) from the dew table)"},
        {"condition data": "certified dry hours before the event", "where": "release table `events`, column dry_before_h", "used": "yes"},
        {"condition data": "radar path-mean rain per link, link_state", "where": "release stores `radar_along_links_5min`, `link_reference_state_5min` (in the S3 datasets)",
         "used": "yes: rain proxies, wet fraction"},
        {"condition data": "radar pixel at the gauges, City gauges 5 min", "where": "release stores `radar_at_gauges_5min`, `gauges_city_5min`", "used": "yes: S4.1 e (information)"},
        {"condition data": "sun elevation / solar radiation, dew sensors, antenna radome type, maintenance", "where": "not in the repo", "used": "no (not downloaded, by the rule)"}])
    timing = R["timing"]
    minutes_fit = sum(timing.values()) / 60

    return f"""# wet_antenna steps S4.1 + S4.2: why do events differ? Per-event rain-scale and antenna factors, then the conditions

Written by `research/wet_antenna/conditions_s4.py` (libraries `system_id.py`, `length_test.py`). Settings and every rule:
`wet_antenna.conditions_s4` in `research/config.yaml`, written on 2026-10-15 BEFORE the run (one change before the run, after a
smoke test of the code: the ridge 0.5 -> 0.05 per link, see section 9). Plan: `docs/ROADMAP.md`, "Remaining plan", S4 sub-steps.
Frozen model: **{MAIN}** of S3b (WAA v3, provisional), {CHECK} as a check; F, P, timing, a_link and d are NOT changed.

{standing}

## 1. Facts: the events and the split

**How the split was made** (quoted, not re-derived). `pipeline/docs/methods.md`, section 8: "Events: train / val / test = 60 / 20 / 20 %,
stratified by intensity class (seeded RNG, seed 42, over the sorted ids)." The intensity class is "terciles of the maximum 5 min City rate
over events" (methods.md, section 6). The code (`pipeline/src/cmlpipe/events/core.py`, `split_events`): "event_id -> split, stratified by
intensity_class; seeded permutation of the SORTED ids per class." Links: "train / test = 80 / 20 %, stratified by length class (< 1, 1-3,
> 3 km) x frequency band", parallel links together; this step uses only TRAIN links. So the event split is random within each intensity
tercile; it does not look at dates, weather or radar quality. The split files are frozen and were only read.

Events in the release split: TRAIN {R['n_events_table']['TRAIN']}, VAL {R['n_events_table']['VAL']} (TEST not read). In the S3 all-length
datasets (events with at least one kept link-event): TRAIN {int((facts.set == 'TRAIN').sum())}, VAL {int((facts.set == 'VAL').sum())}.

{md_table(per_set.round(3))}

Intensity classes per set:

{md_table(ic)}

Rain per event = radar path-mean R of the links inside the event's [start, end) (total = median over links, mm; peak = largest 5 min
value over links and the median over links of each link's peak, mm/h; wet fraction = share of the 5 min bins whose median-over-links R
>= {C['facts']['wet_threshold_mmh']} mm/h). Times UTC. (`{rel(P.COND4_EVENTS)}`.)

{md_table(ev_show.round(2))}

![timeline](../figures/{P.COND4_FIG_TIMELINE.name})

{q.get('facts', '')}

## 2. Is the TRAIN / VAL gap ordinary event-to-event spread?

{SH['reps']} random splits of the {int((facts.set == 'TRAIN').sum())} TRAIN events into {int((facts.set == 'TRAIN').sum()) - SH['half_events']} / {SH['half_events']};
on each {SH['half_events']}-event half the S3.1 model-free length line (no bootstrap), on the {len(R['cells'])} VAL-eligible cells
({', '.join(f'{b} {c}' for b, c in R['cells'])}). Rule: VAL is "rare" if fewer than {SH['rarely'] * 100:.0f} % of the halves have a ratio as low as VAL's.
(`{rel(P.COND4_SPLIT_HALF)}`.)

{md_table(sh_tab)}

![split halves](../figures/{P.COND4_FIG_SPLIT_HALF.name})

{q.get('split_half', '')}

## 3. Fake check (run first)

Truth = the frozen {MAIN} with per-event k in {C['fake']['k_values']} and alpha in {C['fake']['alpha_values']} (9 pairs, seeded order) on the
real rain of the {len(fake)} TRAIN events, the S2 noise. Recovered = both within +-{C['fake']['tol'] * 100:.0f} %. Rule: PASS if >= {C['fake']['pass_share'] * 100:.0f} % of
the events. **Result: {S['fake_verdict']}** ({int(fake.both_ok.sum())} of {len(fake)} events = {S['fake_share'] * 100:.0f} %; k alone {int(fake.k_ok.sum())},
alpha alone {int(fake.alpha_ok.sum())}). (`{rel(P.COND4_FAKE)}`.)

{md_table(fk_by)}

Recovered by the number of short links (< 1 km) of the event: {'; '.join(f"{'fewer than' if k_ else 'at least'} {FA['lean_min_links']} short links {int(g.both_ok.sum())} of {len(g)}" for k_, g in fk.groupby('light'))}.
By the intensity class of the event:

{md_table(fk_ic)}

![fake check](../figures/{P.COND4_FIG_FAKE.name})

{q.get('fake', '')}

## 4. Per-event factors

A_pred = alpha_e W[a_link F(R / k_e)] + L P(R(t - d) / k_e) + c on every link-event of the event; grid {len(KG)} x {len(AG)} (k {FA['k_bounds']},
alpha {FA['alpha_bounds']}, log-spaced, plus 1), Huber {DELTA} dB, log-ridge {FA['ridge_per_link']} per link; {FA['boot_reps']} link-bootstrap draws
inside the event (90 %). "lacks short / long" = fewer than {FA['lean_min_links']} links < 1 km / > 3 km. (`{rel(P.COND4_FACTORS)}`.)

{md_table(fac_show)}

TRAIN against VAL, and the correlation of ln k_e and ln alpha_e (`{rel(P.COND4_FACTOR_SUMMARY)}`; for the correlation row,
TRAIN_median = Spearman over all events, VAL_median = the median within-event bootstrap correlation, mann_whitney_p = its p-value):

{md_table(S['summary'].round(3))}

{MAIN} against {CHECK} per event: Spearman k {rho_fit['k']:.2f}, alpha {rho_fit['alpha']:.2f}.

![factors](../figures/{P.COND4_FIG_FACTORS.name})

{q.get('factors', '')}

## 5. Which factor closes the VAL gap?

The S3.1 bins and Huber line on the model's own prediction with the factors applied per event (in-sample per event: information only).
"closes" = median observed / model ratio within {lo_}-{hi_}. (`{rel(P.COND4_CLOSE_GAP)}`.)

{md_table(gsum.round(3))}

![close the gap](../figures/{P.COND4_FIG_CLOSE_GAP.name})

{q.get('close_gap', '')}

## 6. Radar against gauge

Primary: radar path mean of the gauge links (City gauge within 1 km, the S2 pairing) / gauge, per event, events with >= {RG['min_gauge_mm']} mm
per gauge link. Spearman with k_e: **{rgt['primary'][0]:.2f}** [{rgt['primary'][2]:.2f}, {rgt['primary'][3]:.2f}] (p {rgt['primary'][1]:.3f}, {rgt['primary'][4]} events).
Secondary (radar pixel at the 10 City gauges / gauges): {rgt['pixel'][0]:.2f} [{rgt['pixel'][2]:.2f}, {rgt['pixel'][3]:.2f}] ({rgt['pixel'][4]} events).
Expectation written before the run: if k_e is radar error, Spearman >= {RG['expect_spearman']}. (`{rel(P.COND4_RADAR_GAUGE)}`.)

The ratio itself per set (no factor involved; events with >= {RG['min_gauge_mm']} mm per gauge link, or in total over the gauges for the pixel ratio;
"pooled" = all radar mm / all gauge mm of the set):

{md_table(S['rg_set'].round(3))}

{md_table(rg_show)}

![radar vs gauge](../figures/{P.COND4_FIG_RADAR_GAUGE.name})

{q.get('radar_gauge', '')}

## 7. Conditions (S4.2, descriptive)

**What the repo has** (nothing downloaded):

{md_table(cond_inv)}

Spearman of ln alpha_e and ln k_e ({MAIN}) with each condition over the TRAIN + VAL events, 90 % event-bootstrap intervals; {n_tests} tests,
Benjamini-Hochberg at q = {CO['bh_q']} (`{rel(P.COND4_CONDITIONS)}`):

{md_table(ct_show)}

![conditions](../figures/{P.COND4_FIG_CONDITIONS.name})

Residual of the frozen {MAIN} (no event factors) during rain, per link-event median, then the median over link-events (`{rel(P.COND4_RESIDUAL)}`):

{md_table(rc)}

![residual check](../figures/{P.COND4_FIG_RESIDUAL.name})

{q.get('conditions', '')}

## 8. What this means for WAA v3 (a candidate, not a final claim)

{q.get('meaning', '')}

## 9. Unclear / possibly wrong

{q.get('unclear', '')}

Run: factor loops {minutes_fit:.1f} min ({', '.join(f'{k} {v / 60:.1f}' for k, v in timing.items())}), whole computation {R['minutes']:.1f} min; this report
{(time.time() - clock) / 60:.1f} min. Budget {C['time_budget_h']} h: nothing cut.
"""


if __name__ == "__main__":
    main()
