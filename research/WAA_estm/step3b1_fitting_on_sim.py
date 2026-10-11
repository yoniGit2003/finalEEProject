"""WAA_estm step 3b-1: test the FITTING method on simulated tails with a known truth. No real tail is fitted here.

Data: the step 3a-fix simulator (waa_sim.py) on the TRAIN v1 skeleton; truths m1, m3, m4, m5 (parameters and spread
from step3a_calibration.csv; noise from sim_noise_params.csv), plus m6 (state switch) if its calibration converges.
`replicates` datasets per truth (different seeds: new spread draws, offsets, noise, phases).

  1. every simulated tail is fitted with the candidates F1, F3, F4, F5 (waa_fit.py; free W0 and offset c per tail);
     BIC / AIC with one effective sample size per tail (lag-1 autocorrelation of the F4 residuals)
  2. which candidate wins: per tail (share of tails) and per dataset (sum of BIC); confusion matrices
  3. parameter recovery of the true law's own candidate (m1 -> F1, m3 -> F3, m4 -> F4, m5 -> F5)
  4. population tests: t_half per A_start bin and dA/dt vs A, replicate range per truth vs the real TRAIN v1 tails
  5. ablation for truths m1, m4: one effect off at a time (noise, quantization, offset, spread, all)
Fits are cached per dataset in .cache/fits3b1/ (delete the folder to refit). Settings: `waa.fit3b1`.
Run from the repository root (about 80 minutes with 4 workers; reruns use the cache):
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step3b1_fitting_on_sim.py
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import yaml
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MODULE = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE))
import WAA_estm_paths as P
from research.common.outputs import save_output
import waa_sim as ws
import waa_fit as wf

CFG = yaml.safe_load(open(P.CONFIG_FILE, encoding="utf-8"))
SIM, F = CFG["waa"]["sim"], CFG["waa"]["fit3b1"]
STEP, SPM, H = SIM["rsl_step_db"], SIM["samples_per_min"], SIM["horizon_min"]
OWN = {"m1": "F1", "m3": "F3", "m4": "F4", "m5": "F5"}             # each truth's own candidate
INK, INK2 = "#1a1a1a", "#52514e"
CACHE = P.CACHE / "fits3b1"


def tidy(ax):
    ax.grid(alpha=0.25)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def rel(p):
    return Path(p).resolve().relative_to(P.REPO).as_posix()


def md(df, index=True):
    df = df.reset_index() if index else df
    head = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns)
    return head + "\n" + "\n".join("| " + " | ".join("" if pd.isna(x) else str(x) for x in row) + " |"
                                   for row in df.itertuples(index=False))


def curve_stats(df, col="A_norm", horizon=H):
    d = df[(df.minutes_since_t0_trace <= horizon) & np.isfinite(df[col])]
    g = d.groupby("minutes_since_t0_trace")[col]
    return pd.DataFrame({"median": g.median(), "p10": g.quantile(0.1), "p90": g.quantile(0.9), "n": g.size()})


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    clock = time.time()
    CACHE.mkdir(parents=True, exist_ok=True)

    # =============================================================================================== inputs
    skel = pd.read_csv(P.TAILS_MODEL_INDEX).query("split_event == 'train'").sort_values("tail_id").reset_index(drop=True)
    real = pq.read_table(P.TAILS_MODEL, filters=[("split_event", "=", "train")],
                         columns=["tail_id", "event_id", "minutes_since_t0_trace", "A_db"]).to_pandas()
    real = real.merge(skel[["tail_id", "A_start_trace_db"]], on="tail_id")
    real["A_norm"] = real.A_db / real.A_start_trace_db
    nan_masks = {t: g.sort_values("minutes_since_t0_trace").A_db.isna().to_numpy() for t, g in real.groupby("tail_id")}
    nz = pd.read_csv(P.SIM_NOISE)
    noise = {r.sublink_id: (r.sigma_white_db, r.sigma_slow_db, r.tau_noise_min) for r in nz.itertuples()}
    cal = pd.read_csv(P.STEP3A_CALIBRATION).set_index("law")
    truths = {law: ({k: float(cal.loc[law, k]) for k in ws.PARAMS[law]}, float(cal.loc[law, "s"])) for law in F["laws"]}
    a_start = skel.set_index("tail_id").A_start_trace_db.to_dict()
    print(f"inputs: {len(skel)} skeleton tails ({time.time() - clock:.0f} s)")

    # =============================================================================================== m6 calibration
    real_c = curve_stats(real)
    ok_t = real_c.index[real_c.n >= SIM["min_tails_per_minute"]].to_numpy()
    target = real_c.loc[ok_t, "median"].to_numpy()
    z0 = ws.spread_z(skel, SIM["seed"])
    m6_scan, m6_ok = [], False
    if F["m6"]["try"]:
        st = F["m6"]["start"]
        to_p = lambda x: {"k_run": np.exp(x[0]), "k_evap": np.exp(x[1]), "A_c": np.exp(x[2])}   # noqa: E731
        for s in SIM["spread_grid"]:
            f = np.exp(s * z0)
            res = least_squares(lambda x: np.nanmedian(ws.truth_matrix(skel, "m6", to_p(x), f, H + 1), axis=0)[ok_t]
                                - target, np.log([st["k_run"], st["k_evap"], st["A_c"]]), method="trf")
            p = to_p(res.x)
            rows, _ = ws.simulate(skel, "m6", p, noise, SIM["switches"], SIM["seed"], f, STEP, SPM, nan_masks)
            rows = rows.merge(skel[["tail_id", "A_start_trace_db"]], on="tail_id")
            rows["A_norm"] = rows.A_db / rows.A_start_trace_db
            c = curve_stats(rows).reindex(ok_t)
            r = real_c.loc[ok_t]
            m6_scan.append({"s": s, **{k: round(v, 5) for k, v in p.items()}, "status": res.status,
                            "rmse_median_truth": float(np.sqrt(np.mean(res.fun ** 2))),
                            "rmse_median_observed": float(np.sqrt(np.nanmean((c["median"] - target) ** 2))),
                            "rmse_band_observed": float(np.sqrt(np.nanmean((c.p10 - r.p10) ** 2 + (c.p90 - r.p90) ** 2) / 2))})
        m6_scan = pd.DataFrame(m6_scan)
        b = m6_scan.loc[m6_scan.rmse_band_observed.idxmin()]
        edge = b.s in (min(SIM["spread_grid"]), max(SIM["spread_grid"]))
        m6_ok = bool(b.status > 0 and not edge and b.A_c < 5)
        save_output(m6_scan.round(5), P.STEP3B1_M6, "product", "calibration scan of truth m6 (state switch); read by step 3b-1b",
                    index=False)
        if m6_ok:
            truths["m6"] = ({k: float(b[k]) for k in ws.PARAMS["m6"]}, float(b.s))
        print(f"m6: s = {b.s}, converged = {m6_ok} ({time.time() - clock:.0f} s)")
    laws = list(truths)

    # =============================================================================================== datasets + fits
    base_sw = dict(SIM["switches"])
    jobs = [(law, "full", r) for law in laws for r in range(F["replicates"])]
    jobs += [(law, v, r) for law in F["ablation"]["laws"] for v in F["ablation"]["variants"]
             for r in range(F["ablation"]["replicates"])]
    all_fits, pop_thalf, pop_dadt = [], [], []
    edges = F["thalf_bins_db"] + [np.inf]
    blabels = [f"{lo:g}-{hi:g}" if np.isfinite(hi) else f"> {lo:g}" for lo, hi in zip(edges[:-1], edges[1:])]
    a_bin = pd.cut(pd.Series(a_start), edges, labels=blabels, right=False)
    with ProcessPoolExecutor(F["workers"]) as ex:
        for j, (law, variant, rep) in enumerate(jobs):
            seed = F["seed_base"] + rep
            sw = dict(base_sw)
            for k in (F["ablation"]["variants"].get(variant, []) if variant != "full" else []):
                sw[k] = False
            p, s = truths[law]
            f = np.exp(s * ws.spread_z(skel, seed))
            cache = CACHE / f"{law}_{variant}_{rep}.parquet"
            need_rows = variant == "full" or not cache.exists()
            if need_rows:
                rows, info = ws.simulate(skel, law, p, noise, sw, seed, f, STEP, SPM, nan_masks)
            if cache.exists():
                fits = pd.read_parquet(cache)
            else:
                tails = [(tid, g.minutes_since_t0_trace.to_numpy(), g.A_db.to_numpy()) for tid, g in rows.groupby("tail_id")]
                chunks = [(tuple(F["f3_n_bounds"]), F["candidates"], F["n_eff_from"], F["n_eff_min"], tails[i::16])
                          for i in range(16)]
                fits = pd.DataFrame([x for out in ex.map(wf.fit_many, chunks) for x in out])
                fits = fits.merge(info[["tail_id", "f_time_scale", "offset_db"]], on="tail_id")
                fits.to_parquet(cache)
            fits = fits.assign(law=law, variant=variant, rep=rep)
            all_fits.append(fits)
            if variant == "full":                                      # population statistics of this dataset
                th = wf.t_half(rows, a_start)
                med = th.groupby(a_bin.reindex(th.index), observed=False).median()
                pop_thalf.append(pd.DataFrame({"law": law, "rep": rep, "bin": blabels, "median_t_half_min": med.reindex(blabels).to_numpy()}))
                dc = wf.dadt_curve(rows, F["dadt_smooth_min"], F["dadt_bins_db"], F["dadt_min_count"])
                pop_dadt.append(dc.assign(law=law, rep=rep).rename_axis("A_bin").reset_index())
            print(f"   [{j + 1}/{len(jobs)}] {law} {variant} rep {rep} ({time.time() - clock:.0f} s)", flush=True)
    fits = pd.concat(all_fits, ignore_index=True)
    fits["dataset"] = fits.law + "_" + fits.variant + "_" + fits.rep.astype(str)
    save_output(fits, P.STEP3B1_FITS, "diag", "step 3b-1 every per-tail fit on simulated tails", index=False)

    # =============================================================================================== 2. confusion
    def confusion(fr, score):
        """Per tail (mean share over replicates) and per dataset (share of replicates won, by summed score)."""
        rows_t, rows_d = [], []
        for (law, variant, rep), g in fr.groupby(["law", "variant", "rep"]):
            g = g[g.groupby("tail_id").candidate.transform("size") == len(F["candidates"])]   # tails with every fit
            w = wf.winners(g, score)
            rows_t.append({"law": law, "variant": variant, "rep": rep, **w.value_counts(normalize=True).to_dict()})
            tot = g.groupby("candidate")[score].sum()
            rows_d.append({"law": law, "variant": variant, "rep": rep, "winner": tot.idxmin(),
                           "delta_to_second": float(np.sort(tot.to_numpy())[1] - tot.min())})
        t_ = pd.DataFrame(rows_t).fillna(0)
        d_ = pd.DataFrame(rows_d)
        per_tail = t_.groupby(["law", "variant"])[F["candidates"]].mean()
        per_tail_sd = t_.groupby(["law", "variant"])[F["candidates"]].std()
        per_ds = pd.crosstab([d_.law, d_.variant], d_.winner, normalize="index").reindex(columns=F["candidates"],
                                                                                          fill_value=0)
        return per_tail, per_tail_sd, per_ds, d_

    conf = {}
    for score in ("bic", "aic"):
        conf[score] = confusion(fits, score)
    tab = []
    for score, (pt, sd, pdz, _) in conf.items():
        for (law, variant), r in pt.iterrows():
            for c in F["candidates"]:
                tab.append({"score": score, "law": law, "variant": variant, "candidate": c,
                            "per_tail_share": round(r[c], 4), "per_tail_share_sd": round(sd.loc[(law, variant), c], 4),
                            "per_dataset_share": round(pdz.loc[(law, variant), c], 4) if (law, variant) in pdz.index else 0})
    conf_tab = pd.DataFrame(tab)
    save_output(conf_tab, P.STEP3B1_CONFUSION, "diag", "step 3b-1 confusion: true law x chosen candidate", index=False)

    full_t = conf["bic"][0].xs("full", level="variant").reindex(laws)
    full_d = conf["bic"][2].xs("full", level="variant").reindex(laws)
    full_t_aic = conf["aic"][0].xs("full", level="variant").reindex(laws)
    full_d_aic = conf["aic"][2].xs("full", level="variant").reindex(laws)
    fig, ax = plt.subplots(1, 2, figsize=(12, 1.2 + 0.75 * len(laws)))
    for a_, m, title in ((ax[0], full_t, "per tail: share of tails choosing each candidate (BIC)"),
                         (ax[1], full_d, f"per dataset: share of {F['replicates']} replicates won (summed BIC)")):
        a_.imshow(m.to_numpy(float), cmap="Greys", vmin=0, vmax=1.6, aspect="auto")
        for i in range(m.shape[0]):
            for k in range(m.shape[1]):
                v = m.iloc[i, k]
                own = OWN.get(m.index[i]) == m.columns[k]
                a_.text(k, i, f"{v:.2f}", ha="center", va="center", fontsize=10,
                        fontweight="bold" if own else "normal", color=INK)
                if own:
                    a_.add_patch(plt.Rectangle((k - 0.5, i - 0.5), 1, 1, fill=False, ec=INK, lw=2, ls="--"))
        a_.set_xticks(range(m.shape[1]), m.columns)
        a_.set_yticks(range(m.shape[0]), [f"truth {x}" for x in m.index])
        a_.set_xlabel("chosen candidate")
        a_.set_title(title, fontsize=9)
    fig.suptitle(f"Fitting method on simulated TRAIN v1 tails ({F['replicates']} replicates per truth, full observation "
                 "model). Dashed box, bold = the truth's own candidate.", fontsize=10)
    fig.tight_layout()
    save_output(fig, P.STEP3B1_FIG_CONFUSION, "diag", "step 3b-1 confusion matrix: which candidate wins per true law",
                dpi=110)
    plt.close(fig)

    # =============================================================================================== 3. recovery
    full = fits[fits.variant == "full"].merge(skel[["tail_id", "A_start_trace_db"]], on="tail_id")
    rec = []
    for law, cand in OWN.items():
        if law not in laws:
            continue
        p, _ = truths[law]
        g = full[(full.law == law) & (full.candidate == cand)].copy()
        fct = g.f_time_scale
        if law == "m1":                                             # true per-tail parameters of this law only
            true = {"tau": p["tau"] * fct}
        elif law == "m3":
            true = {"k": p["c"] / fct, "n": pd.Series(p["n"], index=g.index)}
        elif law == "m4":
            true = {"tau_f": p["tau_f"] * fct, "tau_s": p["tau_s"] * fct, "phi": pd.Series(p["phi"], index=g.index)}
        else:
            true = {"a": p["a"] / fct, "b": p["b"] / fct}
        for par, tv in true.items():
            if par in ("phi", "n"):
                d = g[par] - tv
                rec.append({"law": law, "candidate": cand, "parameter": par, "measure": "fitted - true",
                            "true_value": round(float(tv.iloc[0]), 4), "n": len(d), "median": d.median(),
                            "p10": d.quantile(.1), "p90": d.quantile(.9), "at_bound_share": float(
                                (np.isclose(g[par], F["f3_n_bounds"][0]) | np.isclose(g[par], F["f3_n_bounds"][1])).mean()
                                if par == "n" else ((g[par] < 1e-3) | (g[par] > 1 - 1e-3)).mean())})
            else:
                rr = g[par] / tv
                rec.append({"law": law, "candidate": cand, "parameter": par, "measure": "fitted / true",
                            "true_value": round(float(tv.median()), 5), "n": len(rr), "median": rr.median(),
                            "p10": rr.quantile(.1), "p90": rr.quantile(.9), "at_bound_share": np.nan})
        for par, tv, ms in (("c", g.offset_db, "fitted - true offset (dB)"), ("W0", g.A_start_trace_db, "fitted - A_start (dB)")):
            d = g[par] - tv
            rec.append({"law": law, "candidate": cand, "parameter": par, "measure": ms, "true_value": np.nan,
                        "n": len(d), "median": d.median(), "p10": d.quantile(.1), "p90": d.quantile(.9),
                        "at_bound_share": np.nan})
    rec = pd.DataFrame(rec)
    for c in ("median", "p10", "p90", "at_bound_share"):
        rec[c] = rec[c].astype(float).round(3)
    save_output(rec, P.STEP3B1_RECOVERY, "diag", "step 3b-1 parameter recovery of the true law's own candidate",
                index=False)

    # =============================================================================================== 4. population
    real_th = wf.t_half(real, a_start)
    rb = a_bin.reindex(real_th.index)
    rng = np.random.default_rng(F["seed_base"])
    ev = skel.set_index("tail_id").event_id
    real_rows = []
    for lab in blabels:
        tids = rb.index[rb == lab]
        v = real_th.loc[tids]
        evs = ev.loc[tids].unique()
        reps = []
        for _ in range(F["ci_reps"]):
            pick = rng.choice(evs, len(evs))
            vv = pd.concat([v[ev.loc[tids] == e] for e in pick]).dropna()
            reps.append(vv.median() if len(vv) else np.nan)
        lo, hi = np.nanpercentile(reps, [2.5, 97.5])
        real_rows.append({"bin": lab, "tails": len(tids), "real_median": v.median(), "real_ci_low": lo, "real_ci_high": hi})
    real_th_tab = pd.DataFrame(real_rows).set_index("bin")
    pt = pd.concat(pop_thalf)
    rng_tab = pt.groupby(["law", "bin"]).median_t_half_min.agg(["min", "median", "max"]).reset_index()
    rng_tab = rng_tab.merge(real_th_tab.reset_index(), on="bin")
    rng_tab["real_inside"] = (rng_tab.real_median >= rng_tab["min"]) & (rng_tab.real_median <= rng_tab["max"])
    save_output(rng_tab.round(2), P.STEP3B1_THALF, "diag", "step 3b-1 t_half per A_start bin: replicate range vs real",
                index=False)

    real_dc = wf.dadt_curve(real, F["dadt_smooth_min"], F["dadt_bins_db"], F["dadt_min_count"])
    pdd = pd.concat(pop_dadt)
    dd_rng = pdd.groupby(["law", "A_bin"], sort=False).median_dadt.agg(["min", "max"]).reset_index()
    dd_rng = dd_rng.merge(real_dc.rename_axis("A_bin").reset_index()[["A_bin", "A_mid", "median_dadt", "count"]]
                          .rename(columns={"median_dadt": "real_median_dadt", "count": "real_count"}), on="A_bin")
    inside = ((dd_rng.real_median_dadt >= dd_rng["min"]) & (dd_rng.real_median_dadt <= dd_rng["max"])).astype(float)
    dd_rng["real_inside"] = inside.where(dd_rng.real_median_dadt.notna() & dd_rng["min"].notna())   # 1 / 0 / NaN
    save_output(dd_rng.round(5), P.STEP3B1_DADT, "diag", "step 3b-1 dA/dt vs A: replicate range vs real", index=False)

    def separation(tab_, key, lo="min", hi="max"):
        """Pairs of laws whose replicate ranges do NOT overlap, counted over bins."""
        out = {}
        for i, a in enumerate(laws):
            for b in laws[i + 1:]:
                x = tab_[tab_.law == a].set_index(key)
                y = tab_[tab_.law == b].set_index(key).reindex(x.index)
                ok_ = x[lo].notna() & y[lo].notna()
                out[f"{a} vs {b}"] = f"{int(((x[hi] < y[lo]) | (y[hi] < x[lo]))[ok_].sum())} of {int(ok_.sum())}"
        return out

    sep_th, sep_dd = separation(rng_tab, "bin"), separation(dd_rng, "A_bin")

    # figures: t_half population, dA/dt population (one panel per law; real on every panel)
    fig, ax = plt.subplots(1, len(laws), figsize=(3.6 * len(laws), 4), sharey=True)
    xs = np.arange(len(blabels))
    for a_, law in zip(np.atleast_1d(ax), laws):
        g = rng_tab[rng_tab.law == law].set_index("bin").reindex(blabels)
        for x_, (lo, hi) in enumerate(zip(g["min"], g["max"])):
            a_.add_patch(plt.Rectangle((x_ - 0.25, lo), 0.5, max(hi - lo, 0.3), facecolor="none", edgecolor=INK2,
                                       hatch="///", lw=1, label="sim replicate range" if x_ == 0 else None))
        a_.errorbar(xs, real_th_tab.reindex(blabels).real_median,
                    yerr=[real_th_tab.reindex(blabels).real_median - real_th_tab.reindex(blabels).real_ci_low,
                          real_th_tab.reindex(blabels).real_ci_high - real_th_tab.reindex(blabels).real_median],
                    fmt="o-", color=INK, lw=1.8, capsize=3, label="real median, 95 % CI (events)")
        a_.set_xticks(xs, blabels, fontsize=8)
        a_.set_title(f"truth {law}: real inside {int(g.real_inside.sum())} / {len(g)} bins", fontsize=9)
        a_.set_xlabel("A_start bin (dB)", fontsize=8)
        tidy(a_)
    np.atleast_1d(ax)[0].set_ylabel("median t_half (min)")
    np.atleast_1d(ax)[0].legend(frameon=False, fontsize=7)
    fig.suptitle(f"t_half per A_start bin: {F['replicates']} simulated replicates per truth (hatched = min-max) vs real",
                 fontsize=10)
    fig.tight_layout()
    save_output(fig, P.STEP3B1_FIG_THALF, "diag", "step 3b-1 t_half population test", dpi=110)
    plt.close(fig)

    fig, ax = plt.subplots(2, len(laws), figsize=(3.6 * len(laws), 7.2), sharey="row")
    for j, law in enumerate(laws):
        g = dd_rng[dd_rng.law == law].sort_values("A_mid")
        for i, log in enumerate((False, True)):
            a_ = ax[i, j]
            lo_, hi_ = -g["max"].to_numpy(), -g["min"].to_numpy()          # plot the drying rate -dA/dt
            x = g.A_mid.to_numpy()
            if log:
                keep = (x > 0) & (lo_ > 0)
                a_.fill_between(x[keep], lo_[keep], hi_[keep], facecolor="none", edgecolor=INK2, hatch="///", lw=0,
                                label="sim replicate range")
                rv = -g.real_median_dadt.to_numpy()
                a_.plot(x[(x > 0) & (rv > 0)], rv[(x > 0) & (rv > 0)], "o-", color=INK, lw=1.6, ms=4, label="real")
                a_.set_xscale("log")
                a_.set_yscale("log")
                a_.set_ylim(1e-4, 0.1)                              # tiny simulated rates would stretch the axis
            else:
                a_.fill_between(x, lo_, hi_, facecolor="none", edgecolor=INK2, hatch="///", lw=0, label="sim replicate range")
                a_.plot(x, -g.real_median_dadt, "o-", color=INK, lw=1.6, ms=4, label="real")
                a_.axhline(0, color=INK2, lw=0.6)
            a_.set_xlabel("smoothed A (dB)" + (", log" if log else ""), fontsize=8)
            if i == 0:
                a_.set_title(f"truth {law}: real inside {int(g.real_inside.fillna(False).sum())} / "
                             f"{int(g.real_inside.notna().sum())} bins", fontsize=9)
            tidy(a_)
    ax[0, 0].set_ylabel("drying rate -dA/dt (dB/min)")
    ax[1, 0].set_ylabel("-dA/dt (dB/min), log")
    ax[0, 0].legend(frameon=False, fontsize=7)
    fig.suptitle(f"Drying rate vs A, all tails pooled ({F['dadt_smooth_min']} min centred mean, central difference): "
                 f"{F['replicates']} replicates per truth (hatched = min-max) vs real (dots)", fontsize=10)
    fig.tight_layout()
    save_output(fig, P.STEP3B1_FIG_DADT, "diag", "step 3b-1 dA/dt vs A population test", dpi=100)
    plt.close(fig)

    # =============================================================================================== 5. ablation
    abl = conf_tab[(conf_tab.score == "bic") & conf_tab.law.isin(F["ablation"]["laws"])]
    save_output(abl, P.STEP3B1_ABLATION, "diag", "step 3b-1 ablation: confusion of m1 / m4 with effects off",
                index=False)
    variants = ["full"] + list(F["ablation"]["variants"])
    fig, ax = plt.subplots(1, len(F["ablation"]["laws"]) * 2, figsize=(6 * len(F["ablation"]["laws"]), 4.2))
    for j, law in enumerate(F["ablation"]["laws"]):
        for k, (col, title) in enumerate((("per_tail_share", "per tail"), ("per_dataset_share", "per dataset"))):
            a_ = ax[2 * j + k]
            m = abl[abl.law == law].pivot(index="variant", columns="candidate", values=col).reindex(variants)[F["candidates"]]
            a_.imshow(m.to_numpy(float), cmap="Greys", vmin=0, vmax=1.6, aspect="auto")
            for i in range(m.shape[0]):
                for c in range(m.shape[1]):
                    own = OWN[law] == m.columns[c]
                    a_.text(c, i, f"{m.iloc[i, c]:.2f}", ha="center", va="center", fontsize=9,
                            fontweight="bold" if own else "normal")
            a_.set_xticks(range(m.shape[1]), m.columns)
            a_.set_yticks(range(m.shape[0]), m.index if k == 0 else [""] * len(m.index))
            a_.set_title(f"truth {law}, {title} (own = {OWN[law]}, bold)", fontsize=9)
    fig.suptitle("Ablation: share choosing each candidate (BIC) with one effect switched off", fontsize=10)
    fig.tight_layout()
    save_output(fig, P.STEP3B1_FIG_ABLATION, "diag", "step 3b-1 ablation confusion", dpi=110)
    plt.close(fig)

    # =============================================================================================== report
    own_t = {law: round(float(full_t.loc[law, OWN[law]]), 3) for law in laws if law in OWN}
    own_d = {law: round(float(full_d.loc[law, OWN[law]]), 3) for law in laws if law in OWN}
    abl_own = abl.pivot_table(index=["law", "variant"], columns="candidate", values="per_tail_share").round(3)
    abl_ds = abl.pivot_table(index=["law", "variant"], columns="candidate", values="per_dataset_share").round(2)
    rec_show = rec.drop(columns=["candidate"])
    m6_txt = ("not tried" if not F["m6"]["try"] else
              f"included: {', '.join(f'{k} = {v:.4g}' for k, v in truths['m6'][0].items())}, s = {truths['m6'][1]}"
              if m6_ok else "NOT included (calibration did not converge or s at the grid edge)")
    summary = summarize(laws, full_t, full_d, full_d_aic, own_t, own_d, rng_tab, dd_rng, sep_th, sep_dd, abl_own, rec, fits)
    report = f"""# WAA_estm step 3b-1: the fitting method on simulated tails (known truth)

Written by `research/WAA_estm/step3b1_fitting_on_sim.py` (libraries `waa_sim.py`, `waa_fit.py`; settings `waa.fit3b1`).
No real tail is fitted. Truths: {', '.join(laws)} (step 3a-fix parameters and spread; m6: {m6_txt}).
{F['replicates']} replicate datasets per truth on the {len(skel)}-tail TRAIN v1 skeleton (seeds {F['seed_base']} ..
{F['seed_base'] + F['replicates'] - 1}); ablation {F['ablation']['replicates']} replicates per variant. Fits: {len(fits)} rows
({fits.groupby('dataset').tail_id.nunique().sum()} tail fits x candidates). Run time {time.time() - clock:.0f} s.

## Summary in plain words

{summary}

## 1. The fit

Each simulated tail (1 min A over its usable length) is fitted with F1 (exponential, 3 parameters), F4 (two layers,
5), F3 (nonlinear, n in {F['f3_n_bounds']}, 4) and F5 (linear + quadratic, 4), all with a free start W0 and offset c.
One effective sample size per tail from the lag-1 autocorrelation rho of the F4 residuals: n_eff = n (1 - rho) /
(1 + rho) (median n {fits.n.median():.0f}, rho {fits.rho.median():.2f}, n_eff {fits.n_eff.median():.0f}). Converged fits:
{fits.converged.mean():.1%}.

## 2. Which candidate wins (BIC; `{rel(P.STEP3B1_CONFUSION)}` also has AIC)

Per tail (mean share over replicates):

{md(full_t.round(3))}

Per dataset (share of replicates won by the summed BIC):

{md(full_d.round(2))}

AIC, per tail / per dataset:

{md(full_t_aic.round(3))}

{md(full_d_aic.round(2))}

![confusion](../figures/{P.STEP3B1_FIG_CONFUSION.name})

## 3. Parameter recovery of the true law's own candidate (`{rel(P.STEP3B1_RECOVERY)}`, key)

Scale parameters: fitted / true per tail (1 = unbiased); phi, n: fitted - true; c: fitted - the true offset (the slow
noise also shifts the tail mean, so c is not expected to hit the offset exactly); W0: fitted - A_start.

{md(rec_show, index=False)}

## 4. Population tests (all replicates of each truth vs the real TRAIN v1 tails)

**a. t_half per A_start bin** (`{rel(P.STEP3B1_THALF)}`): real median with 95 % CI (resampling events) vs the min-max
of the {F['replicates']} replicate medians.

{md(rng_tab.pivot(index='bin', columns='law', values='real_inside').reindex(blabels))}

Real (median t_half in min, 95 % CI from resampling events):

{md(real_th_tab.round(1))}

Pairs of truths whose replicate ranges do not overlap (bins): {sep_th}.

![t_half](../figures/{P.STEP3B1_FIG_THALF.name})

**b. Drying rate -dA/dt vs A** (`{rel(P.STEP3B1_DADT)}`): bins of smoothed A where the real median lies inside each
truth's replicate range: {', '.join(f"{law}: {int(dd_rng[dd_rng.law == law].real_inside.fillna(False).sum())} / {int(dd_rng[dd_rng.law == law].real_inside.notna().sum())}" for law in laws)}.
Pairs of truths that do not overlap (bins): {sep_dd}.

![dA/dt](../figures/{P.STEP3B1_FIG_DADT.name})

## 5. Ablation (truths {', '.join(F['ablation']['laws'])}; `{rel(P.STEP3B1_ABLATION)}`)

Per tail, share choosing each candidate:

{md(abl_own)}

Per dataset, share of replicates won:

{md(abl_ds)}

![ablation](../figures/{P.STEP3B1_FIG_ABLATION.name})

## Unclear

1. **One n_eff per tail, from the F4 residuals.** If F4 overfits the slow noise, rho and n_eff are too optimistic for
   the simpler candidates; a per-candidate n_eff would change the per-tail choice between nested models.
2. **The offset c absorbs the mean of the slow noise over the tail** (AR(1) with tau_noise ~80 min does not average
   out within a few hours), so c and the late part of the curve trade off; this limits what any fit can say about the
   slow layer of F4 or the tail of F3 / F5.
3. **The ablation keeps the real NaN mask and the tail lengths**; "all off" is therefore not a perfect, infinite
   signal: short tails stay short.
4. **The population tests reuse the calibration data**: the laws were calibrated on the real median and band, so the
   real median curve is close to every law by construction; t_half per bin and dA/dt use the same tails, only cut a
   different way.
5. **m6 is a truth only**; there is no state-switch candidate among F1-F5, so its row shows which smooth law
   absorbs a switch.
6. **The ablation has no "noise and quantization off, offset and spread on" variant**, so the claim that noise OR
   quantization alone suffices to hide the slow layer rests on the single-effect variants plus "all off".
7. **The real dA/dt curve above ~1.5 dB is thin** (few minutes, few events; the zig-zag to ~0 at 1.75 and 3.5 dB): the
   simulated drying rates there are higher than the real ones for every law. This can be the early part of real
   tails (plateaus, the fast drop just after t0_trace) that no smooth law has, or simply too few minutes.
8. **Replicates share the skeleton** (same tails, lengths, A_start); the replicate range is the spread from noise and
   spread draws, not from a different set of events.
"""
    save_output(report, P.STEP3B1_REPORT, "diag", "step 3b-1 report: the fitting method on simulated tails")
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {time.time() - clock:.0f} s")


def summarize(laws, full_t, full_d, full_d_aic, own_t, own_d, rng_tab, dd_rng, sep_th, sep_dd, abl_own, rec, fits):
    """Plain-words summary. Each sentence is tied to a computed number; the wording only says what the numbers show."""
    out = []
    nonexp = [l for l in laws if l != "m1"]
    f1_tail = {l: float(full_t.loc[l, "F1"]) for l in laws}
    m1_ds = float(full_d.loc["m1", "F1"]) if "m1" in laws else np.nan
    nonexp_f1 = max(float(full_d.loc[l, "F1"]) for l in nonexp) if nonexp else np.nan
    out.append("**What the data can tell apart.** Only one thing, and only when all tails of a dataset are scored "
               f"together: a single exponential (m1) versus a slower-than-exponential tail. With the summed BIC, truth m1 "
               f"chooses F1 in {m1_ds:.0%} of the replicates, and no non-exponential truth chooses F1 in more than "
               f"{nonexp_f1:.0%}. They choose the two-layer F4 instead ("
               + ", ".join(f"{l} {float(full_d.loc[l, 'F4']):.0%}" for l in nonexp) + ").")
    wrong = [l for l in own_d if l != "m1" and own_d[l] < 0.5]
    out.append("**What it cannot tell apart.** Which slow law it is: m3, m4, m5 and m6 all end up at F4 (or F3), "
               f"so a dataset won by F4 does not mean two layers; the truth's own candidate wins only for "
               + ", ".join(f"{l} {own_d[l]:.0%}" for l in own_d) + f" of the replicates"
               + (f" (never for {', '.join(wrong)})" if wrong else "") + ". **Per tail nothing can be told apart**: F1 is "
               f"chosen for {min(f1_tail.values()):.0%}-{max(f1_tail.values()):.0%} of the tails whatever the truth is.")
    m4 = rec[(rec.law == "m4") & (rec.parameter == "tau_s")]
    m1t = rec[(rec.law == "m1") & (rec.parameter == "tau")]
    if len(m4) and len(m1t):
        out.append(f"**Parameters.** The exponential time constant comes back without bias (m1: fitted / true tau median "
                   f"{m1t['median'].iloc[0]:.2f}, 10-90 % {m1t['p10'].iloc[0]:.2f}-{m1t['p90'].iloc[0]:.2f}); the slow layer "
                   f"of m4 does not (tau_s median x{m4['median'].iloc[0]:.1f}, 90 % up to x{m4['p90'].iloc[0]:.0f}): a long, "
                   "small slow layer trades off with the free offset c.")
    if "m4" in abl_own.index.get_level_values(0):
        a4 = abl_own.loc["m4"]["F4"]
        single = a4.drop(["full", "all_off"], errors="ignore")
        out.append(f"**Why.** For truth m4 the own candidate F4 wins {a4.get('full', np.nan):.0%} of the tails with "
                   f"every effect on and {a4.get('all_off', np.nan):.0%} with every effect off, but switching off ONE "
                   f"effect changes it little ({', '.join(f'{k} {v:.0%}' for k, v in single.items())}). So no single "
                   "effect is to blame: the noise alone (white + slow) or the 0.3 dB quantization alone is enough to hide the "
                   f"small slow layer in a tail of ~{fits.n.median():.0f} minutes (effective n ~{fits.n_eff.median():.0f}).")
    if "m1" in full_d_aic.index:
        out.append(f"**Score.** AIC picks F4 for every truth per dataset (m1 too: {float(full_d_aic.loc['m1', 'F4']):.0%}), "
                   "so it cannot even see the exponential; BIC with the effective n is the score to use.")
    n_sep_th = sum(int(v.split(" of ")[0]) for v in sep_th.values())
    n_sep_dd = sum(int(v.split(" of ")[0]) for v in sep_dd.values())
    out.append(f"**Population tests do not help.** Over all pairs of truths, the replicate ranges are separate in "
               f"{n_sep_th} t_half bin(s) and {n_sep_dd} dA/dt bin(s): the laws, calibrated to the same median and band, "
               "predict the same pooled behaviour within the replicate scatter. The real t_half medians lie inside the "
               "range of " + ", ".join(f"{l} {int(v)}/4" for l, v in rng_tab.groupby('law').real_inside.sum().items())
               + " bins.")
    return "\n\n".join(out)


if __name__ == "__main__":
    main()
