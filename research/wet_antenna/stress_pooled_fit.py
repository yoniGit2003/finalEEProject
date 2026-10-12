"""wet_antenna step 3b-1c: stress test of the pooled fit (pooled_fit.py). Fake data only; no real tail is fitted.

Why: step 3b-1b met its decision rule under friendly conditions: the fit used the same noise model that made the fake
data, the optimiser started at the calibrated truth, and the real early decay is faster than every simulated truth.
Here we try to BREAK the pooled fit with three tests (window: first 60 min; truths m1, m4, m5, m6; candidates F1, F3,
F4, F5, F6; same likelihood, BIC and event-bootstrap replicates as 3b-1b):
  T1  misspecified noise: the data get 10 s white + TWO 1 min AR(1) wanders (tau 20 and 300 min) + quantization; the
      fit keeps the single AR(1) of sim_noise_params.csv.
  T2  honest starting points: data as in 3b-1b, but the optimiser starts from generic values (the same for every
      truth) and from the F1 solution; the best likelihood of the starts is kept.
  T3  faster truths: every truth is recalibrated to the real median A / A_start at 10, 30 and 60 min (spread s = 1.25),
      plus one variant with s = 1.5 for m1 and m4.
The decision rule is fixed in `wet_antenna.stress_pooled_fit.decision` before the run: per test and truth the own candidate wins
>= 4 of 5 replicates and every key parameter has a median within +-30 % of the truth.

The only uses of real data: the NaN minutes of the real tails (copied into the fake tails, as in 3a), the real median
A / A_start at 10 / 30 / 60 min (the T3 calibration target) and the dry-noise statistics already stored in
sim_noise_params.csv (the T1 check). The frozen tail set v1 is read, never changed.

The fits run in parallel processes (cores - 1), one dataset per job; every finished fit is cached in
.cache/stress_pooled_fit/ (key = test, truth, replicate, candidate, start), so a crash or a restart resumes. Delete that folder
to redo everything. Run from the repository root (about 1-2 hours the first time, 1 minute when cached):
    python -m research.wet_antenna.stress_pooled_fit
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "1")            # one thread per process: the processes run in parallel

from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import sys
import time
import warnings

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from research.wet_antenna import paths as P
from research.common.save import save
from research.common.config import load_config
from research.common.report import md_table, rel
from research.wet_antenna import datasets as dsets
from research.wet_antenna import fit_jobs as fj
from research.wet_antenna import simulator as ws
from research.wet_antenna import pooled_fit as wp

CFG = load_config()
SIM, S = CFG["wet_antenna"]["simulator"], CFG["wet_antenna"]["stress_pooled_fit"]
STEP, SPM, LAGS = SIM["rsl_step_db"], SIM["samples_per_min"], SIM["noise_acf_lags_min"]
OWN = {"m1": "F1", "m4": "F4", "m5": "F5", "m6": "F6"}     # truth -> the candidate with the same form
TESTS = {"T1": "misspecified noise", "T2": "honest starts", "T3": "faster truths (s = 1.25)",
         "T3s": f"faster truths (s = {S['t3']['s_variant']})"}
INK, INK2 = "#1a1a1a", "#52514e"
CACHE = P.CACHE / "stress_pooled_fit"


# =================================================================================================== helpers (main)
def free_names(law):
    """The parameters the T3 recalibration may change (the slow ones in `t3.keep` stay at their 3a-fix value)."""
    return [k for k in ws.PARAMS[law] if k not in S["t3"]["keep"].get(law, [])]


def to_p(law, x, base):
    """Free parameters (log scale; phi through a logistic) -> full parameter dict; the others are taken from `base`."""
    p = dict(base)
    for k, v in zip(free_names(law), x):
        p[k] = 1 / (1 + np.exp(-v)) if k == "phi" else np.exp(v)
    return p


def from_p(law, p):
    return [np.log(p[k] / (1 - p[k])) if k == "phi" else np.log(p[k]) for k in free_names(law)]


def win_figure(dec, replicates):
    """One panel per test: replicates won by the own candidate, per truth. Failing cells are hatched and labelled."""
    tests = [t for t in TESTS if t in set(dec.test)]
    fig, ax = plt.subplots(1, len(tests), figsize=(4.1 * len(tests), 4.6), sharey=True)
    need = S["decision"]["own_win_share"] * replicates
    for a_, t in zip(np.atleast_1d(ax), tests):
        g = dec[dec.test == t].reset_index(drop=True)
        for i, r in g.iterrows():
            a_.bar(i, r.own_wins, width=0.62, facecolor="#d9d9d9" if r.met else "white", edgecolor=INK, lw=1.3,
                   hatch=None if r.met else "///")
            why = "PASS" if r.met else "FAIL" + (f"\n{r.wins_instead}" if r.wins_instead else "") + (
                "\nparams off" if not r.params_ok else "")
            a_.text(i, r.own_wins + 0.12, f"{r.own_wins}/{r.replicates}\n{why}", ha="center", va="bottom", fontsize=9,
                    fontweight="bold" if not r.met else "normal")
        a_.axhline(need, color=INK2, ls="--", lw=1.2)
        a_.text(len(g) - 0.55, need, f"rule: {need:g} of {replicates}", fontsize=8, va="bottom", ha="right", color=INK2)
        a_.set_xticks(range(len(g)), [f"{r.truth}\n(own {r.own})" for r in g.itertuples()])
        a_.set_ylim(0, replicates + 2.2)
        a_.set_yticks(range(replicates + 1))
        a_.set_title(f"{t}: {TESTS[t]}", fontsize=10)
        a_.grid(axis="y", alpha=0.25)
        for sp in ("top", "right"):
            a_.spines[sp].set_visible(False)
    np.atleast_1d(ax)[0].set_ylabel("replicates won by the truth's own candidate")
    fig.suptitle("Stress test of the pooled fit (fake data, first 60 min): does the true law still win?\n"
                 "Grey bar = cell passes the pre-set rule; hatched white bar = fails (label: who wins instead)",
                 fontsize=10)
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    clock = time.time()
    CACHE.mkdir(parents=True, exist_ok=True)
    workers = max((os.cpu_count() or 2) - 1, 1)
    R0, MIN = S["replicates"], S["t3"]["shape_minutes"]

    # ------------------------------------------------------------------------------------------- inputs
    skel = pd.read_csv(P.TAILS_MODEL_INDEX).query("split_event == 'train'").sort_values("tail_id").reset_index(drop=True)
    real = pq.read_table(P.TAILS_MODEL, filters=[("split_event", "=", "train")],
                         columns=["tail_id", "minutes_since_t0_trace", "A_db"]).to_pandas()
    nan_real = {t: g.sort_values("minutes_since_t0_trace").A_db.isna().to_numpy() for t, g in real.groupby("tail_id")}
    real_shape = dsets.shape(real, skel.set_index("tail_id").A_start_trace_db.to_dict(), MIN)
    del real
    nz, noise, _ = dsets.load_noise()
    cal = pd.read_csv(P.SIM_CALIBRATION).set_index("law")
    old = {law: {k: float(cal.loc[law, k]) for k in ws.PARAMS[law]} for law in ("m1", "m3", "m4", "m5")}
    m6 = pd.read_csv(P.PER_TAIL_M6)
    b6 = m6.loc[m6.rmse_band_observed.idxmin()]
    old["m6"] = {k: float(b6[k]) for k in ws.PARAMS["m6"]}
    truths = S["truths"]
    s_true, s_var = S["s_true"], S["t3"]["s_variant"]
    # the calibration start of every candidate (as 3b-1b): the 3a-fix calibration of the truth with the same form
    calib = {"F1": [old["m1"]["tau"]], "F3": [old["m3"]["c"], 1.5],
             "F4": [old["m4"][k] for k in ws.PARAMS["m4"]], "F5": [old["m5"][k] for k in ws.PARAMS["m5"]],
             "F6": [old["m6"][k] for k in ws.PARAMS["m6"]]}
    print(f"{len(skel)} skeleton tails, {skel.event_id.nunique()} events; truths {', '.join(truths)}; "
          f"{workers} worker processes; real A/A_start at {MIN} min: "
          + ", ".join(f"{real_shape[k]:.3f}" for k in MIN))

    # ------------------------------------------------------------------------------------------- T1: the unknown noise
    block = int(round(nz[nz.source == "fit"].median_control_len_min.median() / 10) * 10)
    w_fast, noise_chk, noise_scores = dsets.t1_noise(nz, block)
    T1N = S["t1_noise"]
    noise_t1 = {r.sublink_id: (r.sigma_white_db, np.sqrt(w_fast) * r.sigma_slow_db, T1N["tau_fast_min"],
                               np.sqrt(1 - w_fast) * r.sigma_slow_db, T1N["tau_slow_min"]) for r in nz.itertuples()}
    slow_total = float(nz.sigma_slow_db.median())
    print(f"\nT1 noise: white + AR(1) tau {T1N['tau_fast_min']} min + AR(1) tau {T1N['tau_slow_min']} min; total slow "
          f"sigma per sub-link = its 3a-fix sigma_slow (median {slow_total:.3f} dB); share of the fast part "
          f"{w_fast:g} (sigma {np.sqrt(w_fast) * slow_total:.3f} / {np.sqrt(1 - w_fast) * slow_total:.3f} dB at the "
          f"median). Dry check on {block} min blocks:\n{noise_chk.to_string(index=False)}")
    save(noise_chk, P.STRESS_NOISE,
                index=False)

    # ------------------------------------------------------------------------------------------- T3: faster truths
    zs = [ws.spread_z(skel, 100 + j) for j in range(S["t3"]["z_draws"])]
    target = np.array([real_shape[k] for k in MIN])

    def truth_shape(law, p, s):
        """Noise-free median of W / A_start at the shape minutes (mean over the spread draws)."""
        out = []
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)       # all-NaN rows: tails that end before the minute
            for z in zs:
                m = ws.truth_matrix(skel, law, p, np.exp(s * z), max(MIN) + 5)
                out.append([np.nanmedian(np.nanmean(m[:, k:k + 5], axis=1)) for k in MIN])
        return np.mean(out, axis=0)

    def recalibrate(law, s):
        res = least_squares(lambda x: truth_shape(law, to_p(law, x, old[law]), s) / target - 1, from_p(law, old[law]),
                            method="trf", diff_step=1e-3)
        return {k: float(v) for k, v in to_p(law, res.x, old[law]).items()}

    p_true = {"T1": {t: (old[t], s_true) for t in truths}, "T2": {t: (old[t], s_true) for t in truths},
              "T3": {t: (recalibrate(t, s_true), s_true) for t in truths},
              "T3s": {t: (recalibrate(t, s_var), s_var) for t in S["t3"]["variant_truths"]}}
    print("\nT3 recalibration (target = real median A/A_start at 10 / 30 / 60 min):")
    for test in ("T3", "T3s"):
        for t, (p, s) in p_true[test].items():
            print(f"   {test} {t} s = {s}: " + ", ".join(f"{k} {old[t][k]:.4g} -> {v:.4g}" for k, v in p.items())
                  + "; noise-free A/A_start " + " / ".join(f"{v:.3f}" for v in truth_shape(t, p, s)))

    # ------------------------------------------------------------------------------------------- simulate
    sim_file = lambda test, t, rep: CACHE / f"sim_{test}_{t}_{rep}.parquet"   # noqa: E731
    sim_shape = []
    for rep in range(R0):
        seed = S["seed_base"] + rep
        sk = dsets.bootstrap_skeleton(skel, np.random.default_rng(seed))
        masks = {r.tail_id: nan_real[r.orig_tail_id] for r in sk.itertuples()}
        z = ws.spread_z(sk, seed)
        a0 = sk.set_index("tail_id").A_start_trace_db.to_dict()
        for test in TESTS:
            for t, (p, s) in p_true[test].items():
                f_ = sim_file(test, t, rep)
                if f_.exists():
                    rows = pd.read_parquet(f_, columns=["tail_id", "minutes_since_t0_trace", "A_db"])
                else:
                    rows, _ = ws.simulate(sk, t, p, noise_t1 if test == "T1" else noise, SIM["switches"], seed,
                                          np.exp(s * z), STEP, SPM, masks)
                    rows = rows.merge(sk[["tail_id", "event_copy", "A_start_trace_db"]], on="tail_id")
                    rows.to_parquet(f_)
                sim_shape.append({"test": test, "truth": t, "rep": rep, **dsets.shape(rows, a0, MIN)})
    sim_shape = pd.DataFrame(sim_shape)
    print(f"\nsimulated {len(sim_shape)} datasets ({time.time() - clock:.0f} s)")

    recal = []
    for test in ("T3", "T3s"):
        for t, (p, s) in p_true[test].items():
            g, g_old = sim_shape[(sim_shape.test == test) & (sim_shape.truth == t)], sim_shape[
                (sim_shape.test == "T2") & (sim_shape.truth == t)]
            new_free, old_free = truth_shape(t, p, s), truth_shape(t, old[t], s_true)
            recal.append({"test": test, "truth": t, "s": s,
                          "old_parameters": ", ".join(f"{k} {old[t][k]:.4g}" for k in p),
                          "new_parameters": ", ".join(f"{k} {v:.4g}" for k, v in p.items()),
                          **{f"real_{k}min": round(real_shape[k], 3) for k in MIN},
                          **{f"old_sim_{k}min": round(float(g_old[k].median()), 3) for k in MIN},
                          **{f"new_sim_{k}min": round(float(g[k].median()), 3) for k in MIN},
                          **{f"new_sim_min_{k}min": round(float(g[k].min()), 3) for k in MIN},
                          **{f"new_sim_max_{k}min": round(float(g[k].max()), 3) for k in MIN},
                          **{f"old_noise_free_{k}min": round(float(v), 3) for k, v in zip(MIN, old_free)},
                          **{f"new_noise_free_{k}min": round(float(v), 3) for k, v in zip(MIN, new_free)}})
    recal = pd.DataFrame(recal)
    save(recal, P.STRESS_RECAL, index=False)
    print(recal[["test", "truth"] + [c for c in recal.columns if c.startswith(("real_", "old_sim", "new_sim_1",
          "new_sim_3", "new_sim_6"))]].to_string(index=False))

    # ------------------------------------------------------------------------------------------- jobs
    def jobs_for(replicates, variants, drop_start):
        jobs = []
        for test in TESTS:
            if test == "T3s" and not variants:
                continue
            for t in p_true[test]:
                for rep in range(replicates):
                    starts = {}
                    for c in S["candidates"]:
                        st = {}
                        if test == "T2":                              # generic starts first, then the F1 solution
                            st = {k: [v[0], v[1], S["maxfev_generic"]] for k, v in S["t2_starts"][c].items()
                                  if k != drop_start}
                            if c != "F1":
                                st["f1sol"] = [None, None, S["maxfev_generic"]]
                        st["calib"] = [calib[c], S["s_start"], S["maxfev"]]   # T2: the reference (not in the best-of)
                        starts[c] = st
                    jobs.append({"test": test, "truth": t, "rep": rep, "sim": str(sim_file(test, t, rep)), "cache": str(CACHE),
                                 "starts": starts})
        return jobs

    n_fits = lambda jobs: sum(len(st) for j in jobs for st in j["starts"].values())   # noqa: E731

    # ------------------------------------------------------------------------------------------- time one replicate
    plan_file = CACHE / "plan.json"
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
        print(f"\nplan from {rel(plan_file)}: {plan['replicates']} replicates, variants {plan['variants']}, "
              f"dropped start {plan['drop_start']}; cut: {plan['cuts'] or 'nothing'}")
    else:
        print("\ntiming one complete replicate (T2, m4, replicate 0: all candidates x all starts) ...", flush=True)
        tiny = [{"tail_id": i, "sublink_id": nz.sublink_id.iloc[0], "cml_id": i, "A_start": 1.0,
                 "t": np.arange(20.0), "a": np.exp(-np.arange(20.0) / 30)} for i in range(2)]
        for c in S["candidates"]:                                   # compile the numba code before the clock starts
            wp.pooled_loglik(c, calib[c], 1.0, wp.pack(tiny, noise, {}, None, 10))
        t0 = time.time()
        timed = fj.menu_fits_job(next(j for j in jobs_for(1, False, None) if j["test"] == "T2" and j["truth"] == "m4"))
        t_all = time.time() - t0
        sec = {(r["candidate"], r["start"]): r["seconds"] for r in timed}
        t_calib = sum(v for (c, st), v in sec.items() if st == "calib")

        def projected(replicates, variants, drop_start):
            t_t2 = sum(v for (c, st), v in sec.items() if st != drop_start)
            n_other = sum(j["test"] != "T2" for j in jobs_for(replicates, variants, drop_start))
            return (len(truths) * replicates * t_t2 + n_other * t_calib) / workers

        plan = {"replicates": R0, "variants": True, "drop_start": None, "cuts": [], "wall_s": 0.0,
                "timed_job_s": round(t_all, 1)}
        proj = projected(R0, True, None)
        print(f"   {len(timed)} fits in {t_all:.0f} s ({t_calib:.0f} s for the 5 calibration starts); projected total "
              f"{proj / 3600:.2f} h with {workers} processes (budget {S['time_budget_h']} h)", flush=True)
        for cut, change in (("replicates 5 -> 3", {"replicates": 3}), ("s = 1.5 variants dropped", {"variants": False}),
                            ("T2 starts 4 -> 3 (g3 dropped)", {"drop_start": "g3"})):
            if proj <= S["time_budget_h"] * 3600:
                break
            plan.update(change)
            plan["cuts"].append(cut)
            proj = projected(plan["replicates"], plan["variants"], plan["drop_start"])
            print(f"   cut: {cut} -> projected {proj / 3600:.2f} h", flush=True)
        plan["projected_s"] = round(proj)
        plan_file.write_text(json.dumps(plan))
    R = plan["replicates"]
    jobs = jobs_for(R, plan["variants"], plan["drop_start"])
    total = n_fits(jobs)

    # ------------------------------------------------------------------------------------------- fits (parallel)
    cached0 = len(list(CACHE.glob("fit_*.json")))
    print(f"\n{len(jobs)} datasets, {total} fits ({cached0} already cached) on {workers} processes", flush=True)
    t_fit = time.time()
    fit_rows, done = [], 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = [pool.submit(fj.menu_fits_job, j) for j in jobs]
        for k, fu in enumerate(as_completed(futs)):
            rws = fu.result()
            fit_rows += rws
            done += len(rws)
            el = time.time() - t_fit
            print(f"   datasets {k + 1}/{len(jobs)}, fits {done}/{total}: last {rws[0]['test']} {rws[0]['truth']} "
                  f"rep {rws[0]['rep']} ({el:.0f} s elapsed, about {el / done * (total - done):.0f} s left)", flush=True)
    if len(list(CACHE.glob("fit_*.json"))) > cached0:                # this run did work: add its wall time
        plan["wall_s"] = round(plan["wall_s"] + time.time() - clock, 1)
        plan_file.write_text(json.dumps(plan))

    # ------------------------------------------------------------------------------------------- tables
    fits = pd.DataFrame(fit_rows)
    fits["bic"] = -2 * fits.ll + fits.k * np.log(fits.n_independent)
    for i, nm in enumerate(["p1", "p2", "p3"]):
        fits[nm] = fits.theta.apply(lambda th: th[i] if len(th) > i else np.nan)
    fits["param_names"] = fits.candidate.map(lambda c: ",".join(wp.NAMES[c]))
    fits = fits.sort_values(["test", "truth", "rep", "candidate", "start"]).reset_index(drop=True)
    save(fits.drop(columns=["theta", "theta0"]).round(6), P.STRESS_FITS, index=False)

    key = ["test", "truth", "rep"]
    # the fit that counts: T1 / T3 the calibration start; T2 the best of the generic starts (calib is the reference)
    counts = fits[(fits.test != "T2") | (fits.start != "calib")]
    chosen = counts.loc[counts.groupby(key + ["candidate"]).bic.idxmin()]
    winner = chosen.loc[chosen.groupby(key).bic.idxmin(), key + ["candidate"]].rename(columns={"candidate": "winner"})

    win = (winner.groupby(["test", "truth", "winner"]).size().rename("wins").reset_index()
           .rename(columns={"winner": "candidate"}))
    grid = pd.MultiIndex.from_tuples([(te, t, c) for te in TESTS for t in p_true[te] for c in S["candidates"]
                                      if te in set(fits.test)], names=["test", "truth", "candidate"])
    win = win.set_index(["test", "truth", "candidate"]).reindex(grid, fill_value=0).reset_index()
    win["replicates"] = R
    win["share"] = win.wins / R
    win["own"] = win.candidate == win.truth.map(OWN)

    rec, dec = [], []
    for test in [t for t in TESTS if t in set(fits.test)]:
        for t, (p, s) in p_true[test].items():
            own = OWN[t]
            ch = chosen[(chosen.test == test) & (chosen.truth == t)]
            g = ch[ch.candidate == own]
            for i, name in enumerate(wp.NAMES[own]):                 # candidate and truth share parameter names
                ratio = g[f"p{i + 1}"] / p[name]
                rec.append({"test": test, "truth": t, "candidate": own, "parameter": name, "true": round(p[name], 5),
                            "median_ratio": ratio.median(), "min_ratio": ratio.min(), "max_ratio": ratio.max(),
                            "in_rule": True})
            ratio = g.s / s
            rec.append({"test": test, "truth": t, "candidate": own, "parameter": "s", "true": s,
                        "median_ratio": ratio.median(), "min_ratio": ratio.min(), "max_ratio": ratio.max(),
                        "in_rule": False})
            # own candidate against the best other candidate, per replicate (negative = the own candidate wins)
            piv = ch.pivot(index="rep", columns="candidate", values="bic")
            others = piv.drop(columns=own)
            d = piv[own] - others.min(axis=1)
            lost = d > 0
            instead = others[lost].idxmin(axis=1).value_counts()
            dec.append({"test": test, "truth": t, "own": own, "own_wins": int((~lost).sum()), "replicates": len(d),
                        "wins_instead": ", ".join(f"{c} x{n}" for c, n in instead.items()),
                        "median_dbic_own_minus_best_other": round(float(d.median()), 1),
                        "median_dbic_where_lost": round(float(d[lost].median()), 1) if lost.any() else np.nan})
    rec = pd.DataFrame(rec)
    rec["within_tol"] = (rec.median_ratio - 1).abs() <= S["decision"]["param_tol"]
    for c in ("median_ratio", "min_ratio", "max_ratio"):
        rec[c] = rec[c].round(3)
    dec = pd.DataFrame(dec)
    dec["win_ok"] = dec.own_wins / dec.replicates >= S["decision"]["own_win_share"]
    off = rec[rec.in_rule & ~rec.within_tol].groupby(["test", "truth"]).apply(
        lambda g: ", ".join(f"{r.parameter} x{r.median_ratio}" for r in g.itertuples()), include_groups=False)
    off = off.to_dict() if len(off) else {}
    dec["params_off"] = [off.get((r.test, r.truth), "") for r in dec.itertuples()]
    dec["params_ok"] = dec.params_off == ""
    dec["met"] = dec.win_ok & dec.params_ok
    failed = sorted({"T3" if t == "T3s" else t for t in dec[~dec.met].test})
    verdict = "PASS" if not failed else ("FAIL" if len(failed) == 3 else "PARTIAL")
    save(win, P.STRESS_WIN,
                index=False)
    save(rec, P.STRESS_RECOVERY, index=False)
    save(dec, P.STRESS_DECISION, index=False)

    # ------------------------------------------------------------------------------------------- T2 start sensitivity
    t2 = fits[fits.test == "T2"]
    gen = t2[t2.start != "calib"]
    per = gen.groupby(["truth", "rep", "candidate"]).agg(best=("bic", "min"), worst=("bic", "max"),
                                                         stopped=("success", lambda v: int((~v).sum())),
                                                         n_starts=("bic", "size"))
    per["best_start"] = gen.loc[gen.groupby(["truth", "rep", "candidate"]).bic.idxmin()].set_index(
        ["truth", "rep", "candidate"]).start
    per["calib"] = t2[t2.start == "calib"].set_index(["truth", "rep", "candidate"]).bic
    per["spread"] = per.worst - per.best
    per["gap"] = per.best - per.calib                              # > 0: every generic start is worse than the truth start
    G = S["t2_bic_gap"]
    sens = per.reset_index().groupby(["truth", "candidate"]).apply(lambda g: pd.Series({
        "datasets": len(g), "starts_differ": int((g.spread > G).sum()), "median_bic_range": round(g.spread.median(), 2),
        "max_bic_range": round(g.spread.max(), 1),
        "generic_worse_than_calib": int((g.gap > G).sum()), "generic_better_than_calib": int((g.gap < -G).sum()),
        "fits_stopped_at_maxfev": int(g.stopped.sum()), "fits": int(g.n_starts.sum()),
        "best_start_counts": ", ".join(f"{k} {v}" for k, v in g.best_start.value_counts().items())}),
        include_groups=False).reset_index()
    save(sens, P.STRESS_STARTS, index=False)
    ref = t2[t2.start == "calib"]
    win_ref = ref.loc[ref.groupby(["truth", "rep"]).bic.idxmin(), ["truth", "rep", "candidate"]].rename(
        columns={"candidate": "winner_calib"})
    same = winner[winner.test == "T2"].merge(win_ref, on=["truth", "rep"])
    same["same"] = same.winner == same.winner_calib
    by_cand = per.reset_index().groupby("candidate").apply(lambda g: pd.Series({
        "datasets": len(g), "starts_differ": int((g.spread > G).sum()), "median_bic_range": round(g.spread.median(), 2),
        "max_bic_range": round(g.spread.max(), 1), "generic_worse_than_calib": int((g.gap > G).sum()),
        "generic_better_than_calib": int((g.gap < -G).sum())}), include_groups=False).reset_index()

    # ------------------------------------------------------------------------------------------- figure
    fig = win_figure(dec, R)
    save(fig, P.STRESS_FIG_WIN,
                dpi=110)
    plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    conf = win.pivot(index=["test", "truth"], columns="candidate", values="wins").reindex(
        columns=S["candidates"]).reset_index()
    conf = conf.merge(dec[["test", "truth", "own", "win_ok", "params_ok", "met"]], on=["test", "truth"])
    conf["verdict"] = np.where(conf.met, "pass", "FAIL")
    conf = conf.drop(columns="met")
    fails = dec[~dec.met]
    fail_txt = "\n".join(
        f"- **{r.test} {r.truth}**: own candidate {r.own} wins {r.own_wins}/{r.replicates}"
        + (f"; instead: {r.wins_instead} (median dBIC own - winner in those replicates {r.median_dbic_where_lost:+.1f})"
           if r.wins_instead else "")
        + (f"; parameters off: {r.params_off}" if r.params_off else "") + "."
        for r in fails.itertuples()) or "No cell fails."
    stopped = fits.groupby(["test", "start"]).success.apply(lambda v: f"{int((~v).sum())}/{len(v)}")
    recal_par = recal[["test", "truth", "s", "old_parameters", "new_parameters"]]
    recal_shape = recal[["test", "truth"] + [f"{a}_{k}min" for k in MIN for a in ("real", "old_sim", "new_sim")]]
    recal_free = recal[["test", "truth"] + [f"{a}_{k}min" for k in MIN for a in ("old_noise_free", "new_noise_free")]]
    starts_txt = "\n".join(f"- {c} ({', '.join(wp.NAMES[c])}): " + "; ".join(
        f"{k} = {v[0]}, s {v[1]}" for k, v in S["t2_starts"][c].items() if k != plan["drop_start"])
        for c in S["candidates"])
    starts_txt += ("\n- `f1sol` (every candidate but F1), from the best F1 fit of the same dataset (tau, s): F3 k = "
                   "1 / tau, n = 1.05; F4 tau_f = tau, tau_s = 5 tau, phi = 0.05; F5 a = 1 / tau, b = 0.1 / tau; "
                   "F6 k_run = 0.1 / tau, k_evap = 1 / tau, A_c = 1 dB; s of the F1 fit.")
    report = f"""# wet_antenna step 3b-1c: stress test of the pooled fit (fake data only)

Written by `research/wet_antenna/stress_pooled_fit.py` (libraries `pooled_fit.py`, `simulator.py`; settings
`wet_antenna.stress_pooled_fit`). No real tail is fitted and the frozen tail set v1 is not changed.

**Overall verdict: {verdict}**{"" if not failed else " (tests with a failing cell: " + ", ".join(failed) + ")"}.
PASS = every cell meets the rule; PARTIAL = one or two of the three tests have a failing cell; FAIL = all three do.

## 1. What was run

- Window: the first {S['window_min']} min of every tail. Truths {', '.join(truths)}; candidates
  {', '.join(S['candidates'])}. {R} replicates per truth and test, each an event bootstrap of the {len(skel)}-tail
  TRAIN v1 skeleton (seeds {S['seed_base']} + replicate: new draws, not those of 3b-1b; the same draws for every test,
  so the tests differ only in what they change). Replicate sizes {fits.n_tails.min()}-{fits.n_tails.max()} tails,
  N of the BIC {fits.n_independent.min()}-{fits.n_independent.max()}.
- Same pooled likelihood as 3b-1b (offset prior, ARMA(1,1) whitening with the sub-link's single AR(1), Laplace over
  the per-tail speed), same BIC (N = distinct (link, event draw) pairs).
- **T1** misspecified noise: truths and spread as in 3b-1b, noise with two slow time scales (section 2). Optimiser
  start = the 3a-fix calibration, which here is the truth (so T1 changes only the noise).
- **T2** honest starts: data made exactly as in 3b-1b. {len(S['t2_starts']['F4']) - (plan['drop_start'] is not None) + 1}
  starts per candidate, the same for every truth; the best likelihood counts. The calibration (= truth) start is
  fitted too, only as a reference.
- **T3** faster truths: recalibrated truths (section 3), 3a-fix noise. Optimiser start = the OLD 3a-fix calibration,
  which is no longer the truth. **T3s**: the same with s = {s_var} for {', '.join(S['t3']['variant_truths'])}.
- Nelder-Mead budget: {S['maxfev']} evaluations from the calibration start (as 3b-1b), {S['maxfev_generic']} from
  the T2 starts. Fits stopped at the budget (test, start: stopped / fits): {', '.join(f"{a} {b}: {v}" for (a, b), v in stopped.items())}.
- **Cut for time**: {', '.join(plan['cuts']) if plan['cuts'] else 'nothing'}. One complete replicate (T2, m4: 5
  candidates x all starts) took {plan['timed_job_s']:.0f} s; projected total {plan['projected_s'] / 3600:.2f} h on
  {workers} processes (budget {S['time_budget_h']} h).
- **Total run time**: {plan['wall_s'] / 60:.0f} min wall time for {len(fits)} fits ({fits.seconds.sum() / 3600:.1f}
  process-hours of fitting; cached in `.cache/stress_pooled_fit/`). This run: {time.time() - clock:.0f} s.

## 2. T1 noise check (`{rel(P.STRESS_NOISE)}`)

T1 noise per sub-link = 10 s white (the sub-link's 3a-fix sigma_white) + AR(1) with tau {T1N['tau_fast_min']} min +
AR(1) with tau {T1N['tau_slow_min']} min + 0.3 dB quantization. The two slow parts together have the sub-link's 3a-fix
slow variance (median sigma_slow {slow_total:.3f} dB, so the total slow variance is about 0.15^2 dB^2); the share of
the tau {T1N['tau_fast_min']} min part is {w_fast:g} for every sub-link (median sigmas
{np.sqrt(w_fast) * slow_total:.3f} and {np.sqrt(1 - w_fast) * slow_total:.3f} dB), the grid value that fits the real
statistics best. Dry series of {block} min ({T1N['blocks']} per sub-link, {int((nz.source == 'fit').sum())}
sub-links) against the real control tails; medians over sub-links and the median absolute error per sub-link:

{md_table(noise_chk)}

Score per share of the fast part (sum of squared scaled errors; lower is better):

{md_table(noise_scores)}

## 3. T3 recalibration (`{rel(P.STRESS_RECAL)}`)

Each truth is fitted (least squares on the three relative errors) so that its noise-free median W / A_start, with
its spread, matches the real median A / A_start at {', '.join(map(str, MIN))} min (mean over [k, k + 5) min; the
calibration averages over {S['t3']['z_draws']} spread draws on the full skeleton). The slow parameters
({', '.join(f"{k} {', '.join(v)}" for k, v in S['t3']['keep'].items())}) keep their 3a-fix value: three points inside
the first hour cannot fix a slow time constant (a first attempt with every parameter free gave m4 tau_s = 7e11 min
and m6 k_evap = 4e-20 / min, i.e. a constant floor; no fit was run on those truths).

{md_table(recal_par)}

Real vs simulated A / A_start (simulated = median over the {R0} replicates of the OBSERVED fake tails, with noise and
event bootstrap; old = the 3b-1b truths):

{md_table(recal_shape)}

Noise-free calibration values (full skeleton):

{md_table(recal_free)}

## 4. Win shares and the decision rule (`{rel(P.STRESS_WIN)}`, `{rel(P.STRESS_DECISION)}`)

Replicates (of {R}) won by each candidate (lowest pooled BIC); `win_ok` = own candidate wins >=
{S['decision']['own_win_share']:.0%}; `params_ok` = every key parameter within +-{S['decision']['param_tol']:.0%}.

{md_table(conf)}

Margins (dBIC = BIC of the own candidate minus the best other candidate; negative = the own candidate wins):

{md_table(dec[['test', 'truth', 'own', 'own_wins', 'wins_instead', 'median_dbic_own_minus_best_other', 'median_dbic_where_lost']])}

Failing cells:

{fail_txt}

![win share](../figures/{P.STRESS_FIG_WIN.name})

## 5. Parameter recovery (`{rel(P.STRESS_RECOVERY)}`)

fitted / true of the own candidate: median, minimum and maximum over the {R} replicates. `s` is shown but is not in
the rule.

{md_table(rec.drop(columns=['candidate']))}

## 6. Start sensitivity (T2; `{rel(P.STRESS_STARTS)}`)

Starts:

{starts_txt}

`starts_differ` = datasets where the best and the worst start differ by more than {G:g} BIC;
`generic_worse_than_calib` = datasets where even the best generic start is more than {G:g} BIC worse than the fit
started at the truth (the optimum was missed); `generic_better_than_calib` = the opposite (the truth start stopped
in a worse place).

Per candidate (all truths):

{md_table(by_cand)}

Per truth and candidate:

{md_table(sens)}

Winner with the generic starts vs winner with the truth start: the same in {int(same.same.sum())} of {len(same)}
datasets{"." if same.same.all() else ": different in " + "; ".join(f"{r.truth} rep {r.rep} ({r.winner_calib} -> {r.winner})" for r in same[~same.same].itertuples()) + "."}

## 7. Overall verdict and what it means for Step 4

**{verdict}.** {MEANING}

## 8. Unclear / possibly wrong

{UNCLEAR}
"""
    save(report, P.STRESS_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {time.time() - clock:.0f} s")


# Written after the run of 2026-10-11 (the numbers in the tables are computed on every run; these two texts are the
# reading of that run and quote its numbers: update them if a rerun with other settings changes the tables).
MEANING = """Two of the 14 cells fail, both for truth m5 (linear + quadratic); T2 passes in every cell.

- **What held.** Exponential (m1), two layers (m4) and the state switch (m6) keep winning as their own candidate in
  every test (4-5 of 5 replicates) and their key parameters stay within +-30 % (median), also m4's tau_s and m6's
  A_c: with noise the fit does not know (T1), with starts that are not the truth (T2), with the faster truths (T3)
  and with the wider spread (T3s). The margins of m4 and m6 are large (median dBIC -17 to -69); m1 wins by only
  about 5 BIC in every test.
- **What broke: F5 against F3.** With the unknown noise (T1) the m5 data are won by F3 (dW/dt = -k W^n) in 5 of 5
  replicates, by a median of 1.8 BIC. Both candidates have two parameters, so this is a near tie in likelihood, and
  the m5 truth is almost a pure quadratic (its linear part a = 0.003 / min is a 330 min time constant, invisible in
  60 min), which F3 with n near 2 copies. With the right noise the same data give F5 4 of 5, by 1.9 BIC (T2). In T3
  F5 wins 4 of 5 but its linear rate a is biased (x0.65, range 0.47-1.10). So the fit can say "the decay is
  nonlinear in W" but NOT whether it is a W^n law or linear + quadratic, and it cannot measure the linear rate of F5
  from the first hour.
- **Starts matter, mostly for F6.** The best and worst of the four starts differ by more than 2 BIC in 50 of 100
  (dataset, candidate) pairs: F1 5 / 20, F3 4 / 20, F4 10 / 20, F5 11 / 20, F6 20 / 20 (F6: median spread 59 BIC,
  up to 267). For the true m6, even the best of four starts was more than 2 BIC worse than the truth start in 2 of 5
  datasets (it still won). In one m4 dataset F6 wins from a generic start (by 3.8 BIC) although F4 won when every
  candidate started at its calibration: the single start of 3b-1b had left F6 in a poor optimum there, which flattered F4.

**For Step 4 (real tails).** (1) Go ahead with the pooled fit, but run every candidate from several starts (at least
the four used here; more for F6, e.g. a small grid of A_c) and report the spread over the starts. (2) Read F3 and F5
as ONE family ("nonlinear in W"): a win of one over the other by a few BIC means nothing. (3) A win of F4 over F6
or of F6 over F4 by less than about 10 BIC is not decisive (F6 beat a true m4 by 3.8 and 8.4 BIC in single
replicates); a win of F1 is always by a small margin. (4) The single-AR(1) noise model is good enough for the
exponential / two-layer / switch question; it is not for F3 against F5. (5) Parameters to trust from the first hour:
tau, tau_f, phi, tau_s (the range over replicates reaches 0.5-1.6), k_run, A_c, b; not F5's a.
All of this is on fake tails whose truth is one of the candidates; real tails may follow none of them."""
UNCLEAR = """1. **T3 recalibration is a choice of mine.** With every parameter free, matching the three real points drove m4's
   tau_s and m6's k_evap to a constant floor (tau_s = 7e11 min), which no 60 min window could recover. I stopped
   that run before any T3 fit and kept the slow parameters at their 3a-fix value (`wet_antenna.stress_pooled_fit.t3.keep`). The
   learners may prefer another rule (e.g. an upper bound on tau_s). This was a change after seeing the calibration,
   not after seeing fit results.
2. **The recalibrated truths do not match the real points exactly.** Noise-free they give 0.49-0.63 / 0.24-0.30 /
   0.10-0.14 against the real 0.568 / 0.229 / 0.142: m1 has one parameter and cannot follow the shape (too slow at
   10 and 30 min, too fast at 60), the others overshoot the 10 min drop. The observed fake medians (5 bootstrap
   replicates) are 0.25-0.29 at 30 min: closer to the real 0.229 than the old truths (0.29-0.32), still above it.
   The real early decay is therefore still a little faster than every T3 truth.
3. **Noise-free and observed medians differ** (m1 at 60 min: 0.095 noise-free, 0.143 observed). I believe this is
   the event bootstrap (3 events hold 246 of 429 tails, and links share half the spread) plus the offset, but I
   did not separate the causes.
4. **T1 noise matches the real dry statistics a little worse than the 3a-fix noise**: median absolute error per
   sub-link 0.018 vs 0.011 dB (within-tail std), 0.044 vs 0.032 (ACF 1 min), 0.057 vs 0.052 (10 min), 0.066 vs 0.045
   (60 min); the network medians match as well or better. Reason: one fast share (0.3) for the whole network
   instead of a fit per sub-link, and the total slow variance was kept, not refitted, so the within-tail std is
   about 10 % low (0.127 vs 0.141 dB). I read "total slow variance about 0.15^2" as the sub-link's own 3a-fix value
   (median 0.15 dB), not 0.15 dB for every sub-link.
5. **T1 is one kind of misspecification** (two AR(1) instead of one, same total variance, Gaussian). Real noise may
   also have steps, a daily cycle or TSL jumps; none is tested.
6. **Only 5 replicates**: 4 of 5 and 5 of 5 are not distinguishable statistically, and a median of 5 ratios is
   rough (several ranges reach 0.5-1.6). A cell at exactly 4 of 5 (T2 m1, m4, m5; T3s m4) passes the rule but is
   one replicate from failing.
7. **Replicate sizes vary a lot** (135-753 tails, event bootstrap), more than in 3b-1b (289-832), so power differs
   between replicates.
8. **The reference start in T2 is the truth only for the own candidate.** For the other candidates "calib" is just
   one more start, so `generic_better_than_calib` is expected there; the alarming column is
   `generic_worse_than_calib` (3 of 100, all F6).
9. **Budget**: the generic starts got 800 evaluations, the calibration start 400 (as in 3b-1b); 5 of 100 T2
   calibration-start fits stopped at 400, no generic one stopped at 800. Part of "generic better than calib" may
   be this budget and not the start.
10. **Seeds**: new event draws (3000 + replicate), not the 3b-1b draws, so T2 is not a rerun of the same datasets;
    "the winner with the truth start" is fitted here on the new draws.
"""


if __name__ == "__main__":
    main()
