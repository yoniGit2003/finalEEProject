"""WAA_estm step 3b-1b: a more powerful (pooled) fitting method, re-tested on simulated tails. No real tail is fitted.

Why: in step 3b-1 the per-tail fit with a free offset could only separate "plain exponential" from "slower than
exponential". Here every candidate is fitted to a WHOLE dataset at once (waa_pooled.py):
  - shared law parameters, a per-tail start W0, a per-tail speed factor exp(s z) (half per link, half per tail; s
    fitted; the speed factors are marginalised: Laplace over the per-tail part, quadrature over the per-link part);
  - the offset with a prior c ~ N(0, offset_std^2) of the sub-link (baseline_uncertainty.csv), MAP;
  - the noise in the likelihood: white + AR(1) per sub-link (sim_noise_params.csv) = ARMA(1,1) whitening;
  - candidates F1, F3, F4, F5 and the new F6 (state switch); BIC = -2 log L + k ln N, N = number of distinct
    (link, event draw) pairs (the two directions of a link in one event see the same rain: not independent).
Data: truths m1, m3, m4, m5 (step3a_calibration.csv) and m6 (step3b1_m6_calibration.csv), simulated with waa_sim.py
on 10 EVENT-BOOTSTRAP replicates of the TRAIN v1 skeleton (events drawn with replacement, all their tails taken; the
same draws for every truth). Three fitting windows: first 60 min, first 6 h, whole tail.
Settings: `waa.pooled3b1b` (the decision rule was fixed there before the run). Simulations and fits are cached in
.cache/pooled3b1b/ (delete to redo). Run from the repository root (several hours the first time):
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step3b1b_pooled_fit.py
"""
from pathlib import Path
import json
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MODULE = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE))
import WAA_estm_paths as P
from research.common.save import save
import waa_sim as ws
import waa_fit as wf
import waa_pooled as wp

CFG = yaml.safe_load(open(P.CONFIG_FILE, encoding="utf-8"))
SIM, B = CFG["waa"]["sim"], CFG["waa"]["pooled3b1b"]
STEP, SPM = SIM["rsl_step_db"], SIM["samples_per_min"]
OWN = {"m1": "F1", "m3": "F3", "m4": "F4", "m5": "F5", "m6": "F6"}
LAW_OF = {v: k for k, v in OWN.items()}                            # candidate -> the truth with the same form
KEY_PARAMS = {"m1": {"tau": "tau"}, "m3": {"k": "c", "n": "n"}, "m4": {"tau_f": "tau_f", "tau_s": "tau_s", "phi": "phi"},
              "m5": {"a": "a", "b": "b"}, "m6": {"k_run": "k_run", "k_evap": "k_evap", "A_c": "A_c"}}
INK, INK2 = "#1a1a1a", "#52514e"
CACHE = P.CACHE / "pooled3b1b"


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


def bootstrap_skeleton(skel, rng):
    """Draw events with replacement and take all their tails; new tail ids, event_copy = draw number."""
    evs = np.sort(skel.event_id.unique())
    parts = []
    for k, e in enumerate(rng.choice(evs, len(evs), replace=True)):
        g = skel[skel.event_id == e].copy()
        g["orig_tail_id"] = g.tail_id
        g["tail_id"] = (k + 1) * 100000 + g.tail_id
        g["event_copy"] = k
        parts.append(g)
    return pd.concat(parts, ignore_index=True)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    clock = time.time()
    CACHE.mkdir(parents=True, exist_ok=True)

    # =============================================================================================== inputs
    skel = pd.read_csv(P.TAILS_MODEL_INDEX).query("split_event == 'train'").sort_values("tail_id").reset_index(drop=True)
    real = pq.read_table(P.TAILS_MODEL, filters=[("split_event", "=", "train")],
                         columns=["tail_id", "minutes_since_t0_trace", "A_db"]).to_pandas()
    nan_real = {t: g.sort_values("minutes_since_t0_trace").A_db.isna().to_numpy() for t, g in real.groupby("tail_id")}
    nz = pd.read_csv(P.SIM_NOISE)
    noise = {r.sublink_id: (r.sigma_white_db, r.sigma_slow_db, r.tau_noise_min) for r in nz.itertuples()}
    offset_sd = pd.read_csv(P.BASELINE_UNCERTAINTY).set_index("sublink_id").offset_std.to_dict()
    cal = pd.read_csv(P.STEP3A_CALIBRATION).set_index("law")
    truths = {law: ({k: float(cal.loc[law, k]) for k in ws.PARAMS[law]}, float(cal.loc[law, "s"]))
              for law in B["truths"] if law in cal.index}
    m6 = pd.read_csv(P.STEP3B1_M6)
    b6 = m6.loc[m6.rmse_band_observed.idxmin()]
    truths["m6"] = ({k: float(b6[k]) for k in ws.PARAMS["m6"]}, float(b6.s))
    starts = {"F1": [truths["m1"][0]["tau"]], "F3": [truths["m3"][0]["c"], 1.5],
              "F4": [truths["m4"][0]["tau_f"], truths["m4"][0]["tau_s"], truths["m4"][0]["phi"]],
              "F5": [truths["m5"][0]["a"], truths["m5"][0]["b"]],
              "F6": [truths["m6"][0]["k_run"], truths["m6"][0]["k_evap"], truths["m6"][0]["A_c"]]}
    print(f"truths: {', '.join(truths)}; {len(skel)} skeleton tails, {skel.event_id.nunique()} events")

    # =============================================================================================== simulate
    sims = {}
    for rep in range(B["replicates"]):
        seed = B["seed_base"] + rep
        sk = bootstrap_skeleton(skel, np.random.default_rng(seed))
        masks = {r.tail_id: nan_real[r.orig_tail_id] for r in sk.itertuples()}
        z = ws.spread_z(sk, seed)
        for law, (p, s) in truths.items():
            f_ = CACHE / f"sim_{law}_{rep}.parquet"
            if f_.exists():
                rows = pd.read_parquet(f_)
            else:
                rows, info = ws.simulate(sk, law, p, noise, SIM["switches"], seed, np.exp(s * z), STEP, SPM, masks)
                rows = rows.merge(sk[["tail_id", "event_copy", "A_start_trace_db"]], on="tail_id")
                rows.to_parquet(f_)
            sims[(law, rep)] = rows
    print(f"simulated {len(sims)} datasets ({time.time() - clock:.0f} s)")

    # =============================================================================================== pooled fits
    tasks = [(law, rep, w, c) for law in truths for rep in range(B["replicates"]) for w in B["windows"]
             for c in B["candidates"]]
    fit_rows = []
    packs = {}
    for j, (law, rep, wname, cand) in enumerate(tasks):
        f_ = CACHE / f"fit_{law}_{rep}_{wname}_{cand}.json"
        if f_.exists():
            r = json.loads(f_.read_text())
        else:
            key = (law, rep, wname)
            if key not in packs:
                packs.clear()                                       # one dataset-window in memory at a time
                rows = sims[(law, rep)]
                tails = [{"tail_id": t, "sublink_id": g.sublink_id.iloc[0], "cml_id": g.cml_id.iloc[0],
                          "event_copy": g.event_copy.iloc[0], "A_start": g.A_start_trace_db.iloc[0],
                          "t": g.minutes_since_t0_trace.to_numpy(), "a": g.A_db.to_numpy()} for t, g in rows.groupby("tail_id")]
                D = wp.pack(tails, noise, offset_sd, B["windows"][wname], B["min_points"])
                used = set(D["tail_id"])
                n_ind = len({(tl["cml_id"], tl["event_copy"]) for tl in tails if tl["tail_id"] in used})
                packs[key] = (D, n_ind)
            D, n_ind = packs[key]
            t0 = time.time()
            r = wp.fit(cand, D, starts[cand], B["s_start"], B["maxfev"])
            r.update({"law": law, "rep": rep, "window": wname, "candidate": cand, "n_tails": len(D["tail_id"]),
                      "n_independent": n_ind, "k": len(wp.NAMES[cand]) + 1, "seconds": round(time.time() - t0, 1)})
            f_.write_text(json.dumps(r))
        fit_rows.append(r)
        if (j + 1) % 25 == 0:
            print(f"   fits {j + 1}/{len(tasks)} ({time.time() - clock:.0f} s)", flush=True)
    fits = pd.DataFrame(fit_rows)
    fits["bic"] = -2 * fits.ll + fits.k * np.log(fits.n_independent)
    for i, nm in enumerate(["p1", "p2", "p3"]):
        fits[nm] = fits.theta.apply(lambda th: th[i] if len(th) > i else np.nan)
    fits["param_names"] = fits.candidate.map(lambda c: ",".join(wp.NAMES[c]))
    save(fits.drop(columns=["theta"]).round(6), P.STEP3B1B_FITS, index=False)

    # =============================================================================================== a. confusion
    win = fits.loc[fits.groupby(["law", "rep", "window"]).bic.idxmin(), ["law", "rep", "window", "candidate"]]
    second = fits.sort_values("bic").groupby(["law", "rep", "window"]).bic.apply(lambda b: b.iloc[1] - b.iloc[0])
    conf = {w: pd.crosstab(g.law, g.candidate, normalize="index").reindex(index=list(truths), columns=B["candidates"],
                                                                           fill_value=0)
            for w, g in win.groupby("window")}
    wins = list(B["windows"])
    fig, ax = plt.subplots(1, len(wins), figsize=(5.2 * len(wins), 1.4 + 0.7 * len(truths)))
    for a_, w in zip(np.atleast_1d(ax), wins):
        m = conf[w]
        a_.imshow(m.to_numpy(float), cmap="Greys", vmin=0, vmax=1.6, aspect="auto")
        for i in range(m.shape[0]):
            for k in range(m.shape[1]):
                own = OWN[m.index[i]] == m.columns[k]
                a_.text(k, i, f"{m.iloc[i, k]:.1f}", ha="center", va="center", fontsize=10,
                        fontweight="bold" if own else "normal")
                if own:
                    a_.add_patch(plt.Rectangle((k - 0.5, i - 0.5), 1, 1, fill=False, ec=INK, lw=2, ls="--"))
        a_.set_xticks(range(m.shape[1]), m.columns)
        a_.set_yticks(range(m.shape[0]), [f"truth {x}" for x in m.index])
        a_.set_xlabel("chosen candidate (lowest pooled BIC)")
        a_.set_title(f"window: {w}", fontsize=10)
    fig.suptitle(f"Pooled fit on simulated data: share of {B['replicates']} event-bootstrap replicates won by each "
                 "candidate. Dashed box, bold = the truth's own candidate.", fontsize=10)
    fig.tight_layout()
    save(fig, P.STEP3B1B_FIG_CONFUSION, dpi=110)
    plt.close(fig)

    # =============================================================================================== b. recovery
    rec = []
    for law, (p, s_true) in truths.items():
        cand = OWN[law]
        for w in wins:
            g = fits[(fits.law == law) & (fits.window == w) & (fits.candidate == cand)]
            for i, (pname, tname) in enumerate(KEY_PARAMS[law].items()):
                true = p[tname]
                ratio = g[f"p{i + 1}"] / true
                rec.append({"law": law, "candidate": cand, "window": w, "parameter": pname, "true": round(true, 5),
                            "median_ratio": ratio.median(), "p10_ratio": ratio.quantile(.1),
                            "p90_ratio": ratio.quantile(.9), "replicates": len(ratio)})
            ratio = g.s / s_true
            rec.append({"law": law, "candidate": cand, "window": w, "parameter": "s", "true": s_true,
                        "median_ratio": ratio.median(), "p10_ratio": ratio.quantile(.1), "p90_ratio": ratio.quantile(.9),
                        "replicates": len(ratio)})
    rec = pd.DataFrame(rec)
    rec["within_tol"] = (rec.median_ratio - 1).abs() <= B["decision"]["param_tol"]
    for c in ("median_ratio", "p10_ratio", "p90_ratio"):
        rec[c] = rec[c].round(3)
    save(rec, P.STEP3B1B_RECOVERY,
                index=False)

    # =============================================================================================== d. decision rule
    dec = []
    for law in truths:
        for w in wins:
            share = float(conf[w].loc[law, OWN[law]])
            r_ = rec[(rec.law == law) & (rec.window == w) & (rec.parameter != "s")]
            params_ok = bool(r_.within_tol.all())
            others = conf[w].loc[law].drop(OWN[law])
            dec.append({"law": law, "window": w, "in_rule": law in B["decision"]["laws"], "own_win_share": share,
                        "win_ok": share >= B["decision"]["own_win_share"], "params_ok": params_ok,
                        "params_off": ", ".join(f"{x.parameter} x{x.median_ratio}" for x in r_.itertuples() if not x.within_tol),
                        "main_confusion": f"{others.idxmax()} {others.max():.0%}" if others.max() > 0 else "",
                        "met": (share >= B["decision"]["own_win_share"]) and params_ok})
    dec = pd.DataFrame(dec)
    save(dec, P.STEP3B1B_DECISION, index=False)
    rule = dec[dec.in_rule].groupby("window").met.all()

    # =============================================================================================== c. population
    edges = B["thalf_bins_db"] + [np.inf]
    blabels = [f"{lo:g}-{hi:g}" if np.isfinite(hi) else f"> {lo:g}" for lo, hi in zip(edges[:-1], edges[1:])]

    def per_tail(df, a0):
        """Per tail: t_half (min), A_start and A / A_start at k minutes (mean over [k, k + 5))."""
        out = pd.DataFrame({"t_half": wf.t_half(df, a0)})
        out["a0"] = pd.Series(a0).reindex(out.index)
        for k in B["shape_minutes"]:
            v = df[(df.minutes_since_t0_trace >= k) & (df.minutes_since_t0_trace < k + 5)].groupby("tail_id").A_db.mean()
            out[f"r{k}"] = v.reindex(out.index) / out.a0
        return out

    def stats(pt):
        """Population statistics of one dataset from its per-tail table: median t_half per A_start bin, median
        A / A_start at k minutes."""
        bins = pd.cut(pt.a0, edges, labels=blabels, right=False)
        out = {f"t_half {b} dB": pt.t_half[bins == b].median() for b in blabels}
        for k in B["shape_minutes"]:
            out[f"A/A_start at {k} min"] = pt[f"r{k}"].median()
        return out

    sim_stats = []
    for (law, rep), rows in sims.items():
        a0 = rows.drop_duplicates("tail_id").set_index("tail_id").A_start_trace_db.to_dict()
        sim_stats.append({"law": law, "rep": rep, **stats(per_tail(rows, a0))})
    sim_stats = pd.DataFrame(sim_stats)
    pt_real = per_tail(real, skel.set_index("tail_id").A_start_trace_db.to_dict())
    pt_real["event_id"] = skel.set_index("tail_id").event_id.reindex(pt_real.index)
    real_val = stats(pt_real)
    rng = np.random.default_rng(B["seed_base"])
    evs = np.sort(skel.event_id.unique())
    by_ev = {e: g for e, g in pt_real.groupby("event_id")}
    boot = pd.DataFrame([stats(pd.concat([by_ev[e] for e in rng.choice(evs, len(evs))]))
                         for _ in range(B["ci_reps"])])
    stat_names = list(real_val)
    pop = []
    for st in stat_names:
        lo, hi = np.nanpercentile(boot[st], [2.5, 97.5])
        row = {"statistic": st, "real": real_val[st], "real_ci_low": lo, "real_ci_high": hi}
        for law in truths:
            v = sim_stats.loc[sim_stats.law == law, st]
            row[f"{law}_min"], row[f"{law}_max"] = v.min(), v.max()
            row[f"{law}_real_inside"] = bool(v.min() <= real_val[st] <= v.max())
        sep = [f"{a} vs {b}" for i, a in enumerate(truths) for b in list(truths)[i + 1:]
               if row[f"{a}_max"] < row[f"{b}_min"] or row[f"{b}_max"] < row[f"{a}_min"]]
        row["separates"] = "; ".join(sep)
        pop.append(row)
    pop = pd.DataFrame(pop)
    save(pop.round(4), P.STEP3B1B_POPULATION,
                index=False)
    fig, ax = plt.subplots(2, len(stat_names) // 2, figsize=(4 * len(stat_names) // 2, 7.5))
    for a_, st in zip(ax.ravel(), stat_names):
        r_ = pop.set_index("statistic").loc[st]
        a_.axhspan(r_.real_ci_low, r_.real_ci_high, facecolor="none", edgecolor="#9a9a9a", hatch="///", lw=0,
                   label="real 95 % CI (events)")
        a_.axhline(r_.real, color=INK, lw=1.8, label="real")
        for i, law in enumerate(truths):
            a_.plot([i, i], [r_[f"{law}_min"], r_[f"{law}_max"]], color=INK2, lw=5, solid_capstyle="butt")
        a_.set_xticks(range(len(truths)), list(truths))
        a_.set_title(st + (f"\nseparates: {r_.separates}" if r_.separates else "\nseparates: none"), fontsize=8)
        tidy(a_)
    ax.ravel()[0].legend(frameon=False, fontsize=7)
    fig.suptitle(f"Population statistics: range over {B['replicates']} event-bootstrap replicates per truth (bars) vs "
                 "the real TRAIN v1 value (line) and its 95 % CI from resampling events (hatched)", fontsize=10)
    fig.tight_layout()
    save(fig, P.STEP3B1B_FIG_POPULATION, dpi=100)
    plt.close(fig)

    # =============================================================================================== report
    conf_txt = "\n\n".join(f"**{w}**\n\n{md(conf[w].round(2))}" for w in wins)
    marg = (f"median dBIC to the second candidate: " + ", ".join(
        f"{w} {second.xs(w, level='window').median():.1f}" for w in wins))
    dec_show = dec.assign(own_win_share=dec.own_win_share.round(2))
    rec_show = rec.drop(columns=["candidate"])
    sep_any = pop[pop.separates != ""][["statistic", "separates"]]
    inside = pop[[c for c in pop.columns if c.endswith("_real_inside")]].sum().to_dict()
    summary = plain_summary(truths, wins, conf, dec, rec, pop, rule)
    report = f"""# WAA_estm step 3b-1b: pooled fit, re-tested on simulated tails (known truth)

Written by `research/WAA_estm/step3b1b_pooled_fit.py` (libraries `waa_pooled.py`, `waa_sim.py`, `waa_fit.py`; settings
`waa.pooled3b1b`). No real tail is fitted. Run time {time.time() - clock:.0f} s (fits cached in `.cache/pooled3b1b/`).

## Summary in plain words

{summary}

## Method

- **Data**: truths {', '.join(f"{k} ({ws.LAWS[k]})" for k in truths)}; {B['replicates']} replicates, each an EVENT
  bootstrap of the {len(skel)}-tail TRAIN v1 skeleton ({skel.event_id.nunique()} events drawn with replacement, all
  their tails; the same draws for every truth; replicate sizes {fits.groupby(['law', 'rep']).n_tails.max().min()}-{fits.groupby(['law', 'rep']).n_tails.max().max()} tails).
- **Pooled model**: A = W(t / f_i; theta, W0_i) + c_i + noise. theta shared; W0_i per tail (profiled; exact for F1 / F4,
  two Gauss-Newton steps for F3 / F5 / F6); f_i = exp(s (u_link + v_tail)), u, v ~ N(0, 1/2), s fitted; c_i ~ N(0,
  offset_std^2) MAP (closed form); noise white + AR(1) per sub-link = ARMA(1,1), residuals whitened (Box-Jenkins,
  first minute with its stationary variance; NaN minutes dropped).
- **Marginalisation of the speed factors**: per tail, a LAPLACE (Gaussian) approximation of the profile
  log-likelihood in g = u + v (lattice of {len(wp.G_LATTICE)} points, peak refined, curvature from +-0.05); the
  per-tail part v is then integrated analytically and the per-link part u by quadrature on {len(wp.U_NODES)} nodes of
  N(0, 1/2). Not a full mixed model: W0 and c are profiled, not integrated.
- **Score**: BIC = -2 log L + k ln N with k = law parameters + 1 (s) and N = the number of distinct (link, event draw)
  pairs that have tails in the window (both directions of a link in one event count once; a link in two events or
  two draws counts twice). Median N: {fits.n_independent.median():.0f} (vs {fits.n_tails.median():.0f} tails).
- **Optimiser**: Nelder-Mead in an unconstrained space from the real-data calibration of each candidate's own truth,
  s = {B['s_start']}; converged: {fits.success.mean():.0%} ({int((~fits.success).sum())} fits stopped at maxfev
  {B['maxfev']}); median {fits.nfev.median():.0f} evaluations.

## a. Confusion (share of replicates won, lowest pooled BIC; key figure)

{conf_txt}

{marg}.

![confusion](../figures/{P.STEP3B1B_FIG_CONFUSION.name})

## b. Recovery of the own candidate's parameters (`{rel(P.STEP3B1B_RECOVERY)}`, key)

fitted / true per replicate: median and 10-90 % over the {B['replicates']} replicates; within_tol = |median - 1| <=
{B['decision']['param_tol']}.

{md(rec_show, index=False)}

## c. Population statistics (`{rel(P.STEP3B1B_POPULATION)}`)

Real value with its 95 % CI (resampling events) vs the min-max over the replicates of each truth. Statistics that
separate at least one pair of truths (non-overlapping ranges):

{md(sep_any, index=False) if len(sep_any) else "none"}

Real value inside each truth's range (number of the {len(stat_names)} statistics): {inside}.

{md(pop[['statistic', 'real', 'real_ci_low', 'real_ci_high'] + [c for c in pop.columns if c.endswith(('_min', '_max'))]].round(3), index=False)}

![population](../figures/{P.STEP3B1B_FIG_POPULATION.name})

## d. Decision rule (fixed in `waa.pooled3b1b.decision` before the run)

Useful if, per dataset, the truth's own candidate wins >= {B['decision']['own_win_share']:.0%} of the replicates for
{', '.join(B['decision']['laws'])}, and every key law parameter is recovered within +-{B['decision']['param_tol']:.0%}
(median). Met per window (all four truths): {rule.to_dict()}.

{md(dec_show, index=False)}

## Unclear

1. **The noise model in the fit is the simulator's own** (same sub-link parameters, same white + AR(1) form). On real
   tails the true noise has more than one time scale (3a-fix), so the real-data likelihood will be less sharp than here.
2. **Quantization is not in the likelihood**: it is folded into the white part ((sigma_white^2 + 0.3^2 / 12) / 6),
   which is only right when the 10 s noise dithers the 0.3 dB steps; for the quietest sub-links it is not.
3. **Laplace in g** assumes each tail's profile likelihood is near-Gaussian in log speed; for tails that barely decay
   in the window (first 60 min of slow tails) it is flat and one-sided.
4. **W0 and c are profiled, not integrated** (a few hundred nuisance parameters per dataset); they are the same for
   every candidate, so they do not bias the comparison, but log L is optimistic in absolute terms.
5. **BIC's N** counts (link, event draw) pairs; neighbouring links in one event share weather and are not
   independent either, so N is still an upper bound on the independent information.
6. **One optimiser start per candidate** (the real-data calibration); a few fits stop at maxfev. A worse local optimum
   for the true candidate would show up as a confusion, not as a crash.
7. **m3 is not in the decision rule** (as set); its row is reported for completeness.
8. **The real early decay is faster than every simulated truth**: the real median A / A_start at 30 min
   ({pop.set_index('statistic').loc['A/A_start at 30 min', 'real']:.3f}) lies below the replicate range of every
   truth (lowest minimum {min(pop.set_index('statistic').loc['A/A_start at 30 min', f'{l}_min'] for l in truths):.3f}),
   though inside its own wide CI ({pop.set_index('statistic').loc['A/A_start at 30 min', 'real_ci_low']:.3f}-{pop.set_index('statistic').loc['A/A_start at 30 min', 'real_ci_high']:.3f}).
   The same holds at 60 min (real {pop.set_index('statistic').loc['A/A_start at 60 min', 'real']:.3f}, lowest truth
   minimum {min(pop.set_index('statistic').loc['A/A_start at 60 min', f'{l}_min'] for l in truths):.3f}), and at 180 min only
   m1's range contains the real value. So the real tails fall faster early and keep less of a slow tail than every
   calibrated truth: a hint (the real CIs are wide), not a result. With the calibrated spread s = 1.25 the simulated
   first hour stays a little slow. The real fit (step 4) fits s and the law parameters on the real tails, so it is not
   bound to the calibration.
9. **Replicate sizes vary** (event bootstrap: {fits.groupby(['law', 'rep']).n_tails.max().min()}-{fits.groupby(['law', 'rep']).n_tails.max().max()} tails), so the power varies between
   replicates; the confusion shares are over these varying sizes.
"""
    save(report, P.STEP3B1B_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {time.time() - clock:.0f} s")


def plain_summary(truths, wins, conf, dec, rec, pop, rule):
    """Plain words, every statement tied to a computed number."""
    out = []
    best_w = max(wins, key=lambda w: np.mean([conf[w].loc[l, OWN[l]] for l in truths]))
    diag = {w: ", ".join(f"{l} {conf[w].loc[l, OWN[l]]:.0%}" for l in truths) for w in wins}
    out.append("**Which law wins when it is the truth** (share of the event-bootstrap replicates won by the truth's "
               "own candidate): " + "; ".join(f"{w}: {diag[w]}" for w in wins) + f". The best window is {best_w}.")
    met = [w for w, v in rule.items() if v]
    out.append(f"**Pre-set rule** (own candidate >= 80 % for m1, m4, m5, m6 and key parameters within +-30 %): "
               + (f"met in {', '.join(met)}." if met else "met in no window."))
    fails = dec[dec.in_rule & ~dec.met & (dec.window == best_w)]
    for r in fails.itertuples():
        why = []
        if not r.win_ok:
            why.append(f"own candidate wins {r.own_win_share:.0%}, main confusion {r.main_confusion}")
        if not r.params_ok:
            why.append(f"parameters off: {r.params_off}")
        out.append(f"- {r.law} ({best_w}): " + "; ".join(why) + ".")
    good = rec[(rec.window == best_w) & rec.within_tol & (rec.parameter != "s")]
    out.append(f"**Parameters recovered within +-30 % ({best_w})**: "
               + (", ".join(f"{r.law} {r.parameter} (x{r.median_ratio})" for r in good.itertuples()) or "none") + ".")
    sep = pop[pop.separates != ""]
    out.append("**Population statistics** separate "
               + ("; ".join(f"{r.statistic}: {r.separates}" for r in sep.itertuples()) if len(sep) else "no pair of truths")
               + " (non-overlapping replicate ranges, event bootstrap).")
    conf_m3 = {w: float(conf[w].loc["m3", "F5"]) for w in wins if "m3" in conf[w].index}
    if bool(rule.all()):
        out.append("**What the real tails will be able to tell us (e).** On simulated data with a known noise model, "
                   "the pooled fit tells the five drying forms apart and recovers their parameters (exponential, two "
                   "layers, linear + quadratic and the state switch reliably; the pure power law m3 is partly taken "
                   "for linear + quadratic: F5 wins " + ", ".join(f"{w} {v:.0%}" for w, v in conf_m3.items())
                   + " of the m3 replicates). On the real tails it can therefore say WHICH of the five forms "
                   "describes the drying best, with its parameters and an event-bootstrap uncertainty. It cannot "
                   "prove a mechanism (the true law may be none of the five; the winner is then the best "
                   "approximation). And here the fit used the simulator's own noise model, while real dry noise has "
                   "more than one time scale: before a real-data winner is trusted, a misspecified-noise test is "
                   "needed (truth with two noise time scales, fit with one), because unmodelled slow noise can favour "
                   "the more flexible F4 / F6.")
    else:
        out.append("**What the real tails will be able to tell us (e).** The pre-set rule is not met in every window, "
                   "so a real-data winner must be read with the confusions above in mind.")
    return "\n\n".join(out)


if __name__ == "__main__":
    main()
