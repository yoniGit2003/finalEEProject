"""WAA_estm step 3b-3b: confirmation run of "draw the curve, then name it" on FRESH fake tails. No real tail is fitted.

Step 3b-3 found, AFTER seeing its results, that the free drying curve k(W) and its equation are right above about
0.35 dB and noise-dominated below. A result found post hoc must be confirmed on new data with the rule written first.
This step does that: same fit as 3b-3 (waa_curve.py, waa_pooled.py: offset prior, ARMA(1,1) noise, speed spread
marginalised, 4 starts, penalty by cross-validation over event draws), but
  - new data: event bootstrap with seeds 4000 + replicate (3b-3: 3000), first 60 min only;
  - judged range [W_lo, W_hi] = 0.35 dB to the level with 200 tail-minutes above it; 5 knots inside that range; below
    W_lo the curve is still fitted (straight line in log-log) but nothing is judged or passed to jump 2;
  - band = min-max over the replicates; shape = slope of log k vs log W per replicate (no curvature interval);
  - jump 2 library W, W^2, W^3, max(W - x, 0) for x = 0.5, 0.75 dB.
Settings and the decision rule: `waa.curve3b3b` in research/config.yaml (fixed before the run); the fit settings are
read from `waa.curve3b3`, unchanged. Real data is used only as in 3b-3: the NaN minutes of the real tails and the
dry-noise statistics. The frozen tail set v1 is read, never changed.

Replicate 0 is run first and timed; the projection for all replicates decides the cut (T1 robustness first).
Every stage of every fit is cached in .cache/fits3b3b/, so a restart resumes. Delete that folder to redo.
Run from the repository root (about 40 min the first time, under a minute when cached):
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step3b3b_confirm.py
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "1")            # one thread per process: the processes run in parallel

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import json
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MODULE = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE))
import WAA_estm_paths as P
from research.common.outputs import save_output
import waa_sim as ws
import waa_pooled as wp
import waa_curve as wc
import step3b1c_stress_test as st                          # event bootstrap, noise tables, T1 noise
import step3b3_curve_then_equation as s3                   # load_tails, true_k, start_theta, menu_fit, md, rel

CFG = s3.CFG
SIM, S, FIT, S1C = CFG["waa"]["sim"], CFG["waa"]["curve3b3b"], CFG["waa"]["curve3b3"], CFG["waa"]["stress3b1c"]
STEP, SPM = SIM["rsl_step_db"], SIM["samples_per_min"]
DEC = S["decision"]
W_LO = S["w_lo_db"]
JGRID = S["j2_switch_grid_db"]
JTERMS = wc.curve_terms(JGRID)
KINDS = {"main": "spread 1.25, 3a-fix noise", "T1": "spread 1.25, T1 noise"}
INK, INK2 = "#1a1a1a", "#52514e"
CACHE = P.CACHE / "fits3b3b"
md, rel = s3.md, s3.rel


# =================================================================================================== one job (worker)
_INPUTS = {}


def run_job(job):
    """Free-curve pooled fit of ONE dataset (first 60 min). The stages are those of step 3b-3 (without the curvature
    interval, which is not used here); each is cached in one json file per job."""
    f_ = CACHE / f"job_{job['tag']}.json"
    res = json.loads(f_.read_text()) if f_.exists() else {}
    if res.get("done"):
        return res
    if not _INPUTS:
        _, _INPUTS["noise"], _INPUTS["offset_sd"] = st.load_noise()
    noise, offset_sd = _INPUTS["noise"], _INPUTS["offset_sd"]

    def save():
        tmp = f_.with_suffix(".tmp")
        tmp.write_text(json.dumps(res))
        tmp.replace(f_)

    t0 = time.time()
    tails = s3.load_tails(job["sim"])
    D = wp.pack(tails, noise, offset_sd, S["window_min"], FIT["min_points"])
    used = set(D["tail_id"])
    n_ind = len({(tl["cml_id"], tl["event_copy"]) for tl in tails if tl["tail_id"] in used})
    cv = wc.Curve(job["knots"], prior_sd=FIT["prior_sd_logk"])
    lam0, mf = FIT["lambda_start"], FIT["maxfev"]
    res.update({k: job[k] for k in ("tag", "truth", "kind", "rep", "knots", "w_hi", "minutes_in_range")})
    res.update({"n_tails": len(D["tail_id"]), "n_independent": n_ind, "n_points": int(len(D["T"]))})

    # 1. the starts with the middle penalty; the best one is continued
    res.setdefault("starts", {})
    for name, spec in FIT["starts"].items():
        if name not in res["starts"]:
            x0 = np.r_[s3.start_theta(spec, cv.knots), np.log(FIT["s_start"])]
            res["starts"][name] = cv.fit(D, x0, lam0, mf["start"])
            save()
    res.setdefault("lam", {})
    if str(lam0) not in res["lam"]:
        best = min(res["starts"], key=lambda n: res["starts"][n]["cost"])
        res["lam"][str(lam0)] = {**cv.fit(D, res["starts"][best]["x"], lam0, mf["finish"]), "from": best}
        save()
    x_mid = res["lam"][str(lam0)]["x"]
    # 2. the other penalties, from that solution
    for lam in FIT["lambda_grid"]:
        if str(lam) not in res["lam"]:
            res["lam"][str(lam)] = cv.fit(D, x_mid, lam, mf["other_lambda"])
            save()
    # 3. 2-fold cross-validation over event draws: fit on one half, log-likelihood of the other half
    res.setdefault("cv", {})
    evs = np.random.default_rng([S["seed_base"], job["rep"]]).permutation(sorted({tl["event_copy"] for tl in tails}))
    half = [set(evs[::2].tolist()), set(evs[1::2].tolist())]
    Dh = None
    for lam in FIT["lambda_grid"]:
        if str(lam) not in res["cv"]:
            if Dh is None:
                Dh = [wp.pack([tl for tl in tails if tl["event_copy"] in h], noise, offset_sd, S["window_min"],
                              FIT["min_points"]) for h in half]
            out = 0.0
            for a, b in ((0, 1), (1, 0)):
                r = cv.fit(Dh[a], res["lam"][str(lam)]["x"], lam, mf["cv"])
                out += cv.loglik(np.asarray(r["x"]), Dh[b])
            res["cv"][str(lam)] = out
            save()
    lam_star = max(FIT["lambda_grid"], key=lambda v: res["cv"][str(v)])
    fit = res["lam"][str(lam_star)]
    res.update({"lam_star": lam_star, "x": fit["x"], "ll": fit["ll"], "s": float(np.exp(fit["x"][-1])),
                "bic": -2 * fit["ll"] + (cv.K + 1) * np.log(n_ind),
                "nfev": sum(r["nfev"] for r in list(res["starts"].values()) + list(res["lam"].values()))})
    # 4. F4 menu fit (m4 only): for the "free curve worse than F4" warning sign
    if job["menu_f4"]:
        s3.menu_fit(job, D, n_ind, "F4", S["menu_f4_starts"], res, save)
    res["seconds"] = round(res.get("seconds", 0.0) + time.time() - t0, 1)
    res["done"] = True
    save()
    return res


# =================================================================================================== figure
def k_figure(truth, reps, w, K, par, cover):
    """True k(W) (dashed) vs learned k(W) (solid), band of the replicates (hatched) inside the judged range (shaded)."""
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    wf = np.geomspace(0.05, w[-1], 120)                    # drawn also below W_lo, to show what is NOT judged
    Kf = np.array([wc.Curve(r["knots"]).k(np.asarray(r["x"])[:-1], wf) for r in reps])
    ax.axvspan(w[0], w[-1], color="#e9e9e9", lw=0, label=f"judged range {w[0]:.2f}-{w[-1]:.2f} dB (shaded)")
    ax.fill_between(w, K.min(axis=0), K.max(axis=0), facecolor="none", edgecolor="#6f6f6f", hatch="///", lw=0,
                    label=f"band: min-max of the {len(reps)} replicates (hatched)")
    for k_ in Kf:
        ax.plot(wf, k_, "-", color=INK2, lw=0.7)
    med = np.exp(np.median(np.log(Kf), axis=0))
    ax.plot(wf, med, "-", color=INK, lw=2.4, label="learned k(W): median (thick solid), each replicate (thin solid)")
    ax.text(wf[-1], med[-1], " learned", fontsize=9, va="center")
    kt = s3.true_k(truth, par, wf)
    if kt is not None:
        ax.plot(wf, kt, "--", color=INK, lw=2.0, label="true k(W) (dashed)")
        ax.text(wf[-1], kt[-1] * 0.9, " true", fontsize=9, va="top")
    else:                                                  # m4: the rate depends on where the tail started
        for a0, ls, name in ((1.0, ":", "dotted"), (3.0, "-.", "dash-dotted")):
            t = np.linspace(0, 720, 4000)
            wt = ws.law_curve("m4", par["m4"], t, a0)
            ax.plot(wt[1:], -np.diff(wt) / np.diff(t) / wt[1:], ls, color=INK, lw=1.6,
                    label=f"true rate along a tail that starts at {a0:g} dB ({name})")
    ax.axvline(w[0], color=INK, lw=1.0)
    ax.set(xscale="log", yscale="log", xlim=(wf[0], wf[-1] * 1.4), xlabel="W (dB)",
           ylabel="k(W) = drying rate per dB of water (1 / min)")
    ax.text(w[0] * 0.95, ax.get_ylim()[1], "below: noise-dominated, not judged  ", fontsize=8.5, ha="right", va="top",
            rotation=0)
    ax.grid(alpha=0.25, which="both")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    ax.set_title(f"Truth {truth} ({ws.LAWS[truth]}), spread 1.25, 3a-fix noise, first 60 min, fresh data\n"
                 + (f"true k inside the band on {cover:.0%} of the judged range" if np.isfinite(cover) else
                    "no single true curve: the rate depends on the tail"), fontsize=10)
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    clock = time.time()
    CACHE.mkdir(parents=True, exist_ok=True)
    workers = max((os.cpu_count() or 2) - 1, 1)
    R = S["replicates"]

    # ------------------------------------------------------------------------------------------- inputs and truths
    skel = pd.read_csv(P.TAILS_MODEL_INDEX).query("split_event == 'train'").sort_values("tail_id").reset_index(drop=True)
    real = pq.read_table(P.TAILS_MODEL, filters=[("split_event", "=", "train")],
                         columns=["tail_id", "minutes_since_t0_trace", "A_db"]).to_pandas()
    nan_real = {t: g.sort_values("minutes_since_t0_trace").A_db.isna().to_numpy() for t, g in real.groupby("tail_id")}
    del real
    nz, noise, _ = st.load_noise()
    cal = pd.read_csv(P.STEP3A_CALIBRATION).set_index("law")
    par = {law: {k: float(cal.loc[law, k]) for k in ws.PARAMS[law]} for law in ("m1", "m3", "m4", "m5")}
    m6 = pd.read_csv(P.STEP3B1_M6)
    b6 = m6.loc[m6.rmse_band_observed.idxmin()]
    par["m6"] = {k: float(b6[k]) for k in ws.PARAMS["m6"]}
    true_coef = {"m1": {"W": 1 / par["m1"]["tau"]}, "m3": {"W^2": par["m3"]["c"]},
                 "m5": {"W": par["m5"]["a"], "W^2": par["m5"]["b"]}}
    datasets = [("main", t) for t in S["truths"]] + [("T1", t) for t in S["t1_truths"]]
    print(f"{len(skel)} skeleton tails; datasets: " + ", ".join(f"{k} {t}" for k, t in datasets) + f"; {workers} processes")

    # ------------------------------------------------------------------------------------------- simulate (cached)
    sim_file = lambda kind, t, rep: CACHE / f"sim_{kind}_{t}_{rep}.parquet"   # noqa: E731
    noise_t1 = None
    meta = {}
    for rep in range(R):
        seed = S["seed_base"] + rep
        sk = masks = z = None
        for kind, t in datasets:
            f_ = sim_file(kind, t, rep)
            if not f_.exists():
                if sk is None:
                    sk = st.bootstrap_skeleton(skel, np.random.default_rng(seed))
                    masks = {r.tail_id: nan_real[r.orig_tail_id] for r in sk.itertuples()}
                    z = ws.spread_z(sk, seed)
                if kind == "T1" and noise_t1 is None:      # the T1 noise of step 3b-1c (same function, same settings)
                    block = int(round(nz[nz.source == "fit"].median_control_len_min.median() / 10) * 10)
                    w_fast, _, _ = st.t1_noise(nz, block)
                    T1N = S1C["t1_noise"]
                    noise_t1 = {r.sublink_id: (r.sigma_white_db, np.sqrt(w_fast) * r.sigma_slow_db, T1N["tau_fast_min"],
                                               np.sqrt(1 - w_fast) * r.sigma_slow_db, T1N["tau_slow_min"])
                                for r in nz.itertuples()}
                rows, _ = ws.simulate(sk, t, par[t], noise_t1 if kind == "T1" else noise, SIM["switches"], seed,
                                      np.exp(S["s_true"] * z), STEP, SPM, masks)
                rows.merge(sk[["tail_id", "event_copy", "A_start_trace_db"]], on="tail_id").to_parquet(f_)
            # judged range and knots of this dataset (first 60 min)
            a = pd.read_parquet(f_, columns=["minutes_since_t0_trace", "A_db"])
            a = a.A_db[(a.minutes_since_t0_trace < S["window_min"]) & a.A_db.notna()].to_numpy()
            w_hi = float(np.sort(a)[-FIT["support_min_minutes"]])
            inside = np.log(a[(a >= W_LO) & (a <= w_hi)])
            q = np.quantile(inside, np.linspace(0, 1, S["knots"])[1:-1])
            meta[(kind, t, rep)] = {"w_hi": w_hi, "points": int(len(a)), "minutes_in_range": int(len(inside)),
                                    "knots": [W_LO] + [float(v) for v in np.exp(q)] + [w_hi]}
    print(f"datasets ready ({time.time() - clock:.0f} s); judged range {W_LO} dB to "
          f"{min(m['w_hi'] for m in meta.values()):.2f}-{max(m['w_hi'] for m in meta.values()):.2f} dB")

    def jobs_for(reps, kinds):
        jobs = [{"tag": f"{kind}_{t}_{rep}", "truth": t, "kind": kind, "rep": rep, "sim": str(sim_file(kind, t, rep)),
                 "knots": m["knots"], "w_hi": m["w_hi"], "points": m["points"], "minutes_in_range": m["minutes_in_range"],
                 "menu_f4": kind == "main" and t == "m4"}
                for (kind, t, rep), m in meta.items() if rep in reps and kind in kinds]
        return sorted(jobs, key=lambda j: (-j["menu_f4"], -j["points"]))   # the long jobs first

    def run(jobs, label):
        t_run, out = time.time(), []
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(run_job, j): j for j in jobs}
            for k, fu in enumerate(as_completed(futs)):
                out.append(fu.result())
                print(f"   {label} {k + 1}/{len(jobs)}: {futs[fu]['tag']} ({time.time() - t_run:.0f} s)", flush=True)
        return out, time.time() - t_run

    # ------------------------------------------------------------------------------------------- time replicate 0, plan
    plan_file = CACHE / "plan.json"
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
        print(f"plan from {rel(plan_file)}: cut: {plan['cuts'] or 'nothing'}")
    else:
        print(f"\ntiming replicate 0 ({len(datasets)} fits) ...", flush=True)
        r0, wall0 = run(jobs_for([0], list(KINDS)), "replicate 0")
        sec = {k: sum(r["seconds"] for r in r0 if r["kind"] == k) for k in KINDS}
        proj = R * wall0                                   # every replicate has the same datasets
        plan = {"kinds": list(KINDS), "cuts": [], "wall_rep0_s": round(wall0, 1), "process_s_rep0": sec,
                "wall_s": round(wall0, 1)}
        print(f"   replicate 0: {wall0 / 60:.1f} min wall ({sec['main'] / 60:.1f} process-min main, {sec['T1'] / 60:.1f} T1); "
              f"projected total {proj / 3600:.2f} h on {workers} processes (budget {S['time_budget_h']} h)", flush=True)
        if proj > S["time_budget_h"] * 3600:
            plan["kinds"] = ["main"]
            plan["cuts"].append("T1 robustness")
            proj = wall0 + (R - 1) * wall0 * sec["main"] / sum(sec.values())
            print(f"   cut: T1 robustness -> projected {proj / 3600:.2f} h", flush=True)
        plan["projected_s"] = round(proj)
        plan_file.write_text(json.dumps(plan))
    jobs = jobs_for(range(R), plan["kinds"])
    done0 = sum(json.loads((CACHE / f"job_{j['tag']}.json").read_text()).get("done", False)
                for j in jobs if (CACHE / f"job_{j['tag']}.json").exists())
    print(f"\n{len(jobs)} fits ({done0} already finished) on {workers} processes", flush=True)
    results, wall = run(jobs, "fits")
    if done0 < len(jobs):
        plan["wall_s"] = round(plan["wall_s"] + wall, 1)
        plan_file.write_text(json.dumps(plan))

    # ------------------------------------------------------------------------------------------- curves on the judged grid
    by = {(r["kind"], r["truth"], r["rep"]): r for r in results}
    cfgs = [(k, t) for k, t in datasets if k in plan["kinds"]]
    need = DEC["min_reps"]
    fit_rows, cur_rows, cs_rows, sel_rows, coef_rows, j2_rows, store = [], [], [], [], [], [], {}
    for kind, truth in cfgs:
        reps = [by[(kind, truth, r)] for r in range(R)]
        w = np.geomspace(W_LO, min(r["w_hi"] for r in reps), FIT["grid_points"])   # judged range common to the replicates
        aw = wc.slope_weights(np.log(w))
        K = np.array([wc.Curve(r["knots"]).k(np.asarray(r["x"])[:-1], w) for r in reps])
        slopes = np.log(K) @ aw
        kt = s3.true_k(truth, par, w)
        store[(kind, truth)] = (reps, w, K)
        for r, k_, sl in zip(reps, K, slopes):
            costs = {n: v["cost"] for n, v in r["starts"].items()}
            fit_rows.append({"kind": kind, "truth": truth, "rep": r["rep"], "n_tails": r["n_tails"],
                             "n_independent": r["n_independent"], "w_lo_db": W_LO, "w_hi_db": round(r["w_hi"], 3),
                             "minutes_in_range": r["minutes_in_range"],
                             "knots_db": " ".join(f"{v:.3g}" for v in r["knots"]), "lambda": r["lam_star"],
                             "s": round(r["s"], 4), "ll": round(r["ll"], 2), "bic": round(r["bic"], 2),
                             "best_start": r["lam"][str(FIT["lambda_start"])]["from"],
                             "start_cost_spread": round(max(costs.values()) - min(costs.values()), 2),
                             "nfev": r["nfev"], "seconds": r["seconds"], "slope": round(float(sl), 4),
                             **{f"k_knot{i + 1}": round(float(np.exp(v)), 5) for i, v in enumerate(r["x"][:-1])}})
            for j in range(len(w)):
                cur_rows.append({"kind": kind, "truth": truth, "rep": r["rep"], "w_db": w[j], "k": float(k_[j]),
                                 "k_true": np.nan if kt is None else float(kt[j])})
        # curve and shape
        row = {"kind": kind, "truth": truth, "w_lo_db": W_LO, "w_hi_db": round(w[-1], 3),
               "band_coverage": np.nan if kt is None else round(float(np.mean((kt >= K.min(axis=0)) & (kt <= K.max(axis=0)))), 3),
               "median_k_ratio": np.nan if kt is None else round(float(np.exp(np.median(np.log(K / kt)))), 3),
               "true_slope": np.nan if kt is None else round(float(np.log(kt) @ aw), 3),
               "slope_min": round(float(slopes.min()), 3), "slope_median": round(float(np.median(slopes)), 3),
               "slope_max": round(float(slopes.max()), 3),
               "slopes": " ".join(f"{v:.2f}" for v in slopes),
               "flat_reps": int((np.abs(slopes) < DEC["flat_abs_slope"]).sum()),
               "rising_reps": int((slopes > DEC["rising_slope"]).sum()), "shape_wanted": DEC["shape"].get(truth, "-")}
        if truth in DEC["shape"]:
            row["curve_ok"] = bool(row["band_coverage"] >= DEC["band_coverage"]
                                   and abs(row["median_k_ratio"] - 1) <= DEC["level_tol"])
            row["shape_ok"] = bool(row[f"{row['shape_wanted']}_reps"] >= need)
        cs_rows.append(row)
        # jump 2: the equation of h(W) = W k(W) on the judged grid
        lam, sel, coef = wc.name_curves(w, K * w, JGRID, FIT["j2_lambda_grid"])
        rule = DEC["j2_terms"].get(truth)
        ix = {t: i for i, t in enumerate(JTERMS)}
        for t in JTERMS:
            ro = "-" if rule is None else ("true" if t in rule["required"] else ("allowed" if t in rule["allowed"] else "spurious"))
            sel_rows.append({"kind": kind, "truth": truth, "term": t, "role": ro, "n_selected": int(sel[:, ix[t]].sum()),
                             "replicates": R, "threshold": lam})
        jr = {"kind": kind, "truth": truth, "threshold": lam,
              "returned": ", ".join(f"{k} x{v}" for k, v in pd.Series(
                  [" + ".join(t for t, s_ in zip(JTERMS, s) if s_) or "(nothing)" for s in sel]).value_counts().items()),
              "coefficients_median": ", ".join(f"{t} {np.median(coef[:, ix[t]]):.4f}" for t in JTERMS if sel[:, ix[t]].sum() * 2 > R)}
        if rule is not None:
            need_n = {t: int(sel[:, ix[t]].sum()) for t in rule["required"]}
            spur = {t: int(sel[:, ix[t]].sum()) for t in JTERMS if t not in rule["allowed"] and sel[:, ix[t]].any()}
            ratios = {}
            for c, tc in true_coef[truth].items():
                ratio = coef[:, ix[c]] / tc
                counted = c in rule["coefs"]
                coef_rows.append({"kind": kind, "truth": truth, "coefficient": c, "true": round(tc, 5),
                                  "median_ratio": round(float(np.median(ratio)), 3), "min_ratio": round(float(ratio.min()), 3),
                                  "max_ratio": round(float(ratio.max()), 3),
                                  "replicates_nonzero": int((ratio != 0).sum()), "counted": counted})
                if counted:
                    ratios[c] = float(np.median(ratio))
            jr.update({"true_selected": "; ".join(f"{t} {k}/{R}" for t, k in need_n.items()),
                       "spurious": ", ".join(f"{t} {v}/{R}" for t, v in spur.items()),
                       "coef_ratios": ", ".join(f"{c} x{v:.2f}" for c, v in ratios.items()),
                       "true_ok": all(k >= need for k in need_n.values()),
                       "spurious_ok": all(v <= DEC["max_spurious"] for v in spur.values()),
                       "coef_ok": all(abs(v - 1) <= DEC["coef_tol"] for v in ratios.values())})
            jr["equation_ok"] = bool(jr["true_ok"] and jr["spurious_ok"] and jr["coef_ok"])
        j2_rows.append(jr)
    fits, cur, cs = pd.DataFrame(fit_rows), pd.DataFrame(cur_rows), pd.DataFrame(cs_rows)
    sel, coefs, j2 = pd.DataFrame(sel_rows), pd.DataFrame(coef_rows), pd.DataFrame(j2_rows)
    save_output(fits, P.STEP3B3B_FITS, "diag", "step 3b-3b every free-curve pooled fit: judged range, knots, penalty, k at "
                "the knots, s, BIC, starts, slope over the judged range", index=False)
    save_output(cur.round(6), P.STEP3B3B_CURVES, "product", "step 3b-3b learned k(W) per replicate on the judged grid "
                "(0.35 dB to W_hi) and the true k(W)", index=False)
    save_output(cs, P.STEP3B3B_CURVE_SHAPE, "product", "step 3b-3b curve and shape per truth and data kind: band coverage, "
                "level, slope of log k vs log W per replicate", index=False)
    save_output(sel, P.STEP3B3B_J2_SEL, "product", "step 3b-3b jump 2: replicates in which each library term is selected",
                index=False)
    save_output(coefs, P.STEP3B3B_J2_COEF, "product", "step 3b-3b jump 2: fitted / true coefficients of the true terms",
                index=False)

    # ------------------------------------------------------------------------------------------- m4: free curve vs F4
    m4 = pd.DataFrame([{"rep": r, "bic_free_curve": round(by[("main", "m4", r)]["bic"], 1),
                        "bic_f4": round(by[("main", "m4", r)]["menu_F4"]["bic"], 1),
                        "f4_best_start": by[("main", "m4", r)]["menu_F4"]["start"],
                        "curve_minus_f4_bic": round(by[("main", "m4", r)]["bic"] - by[("main", "m4", r)]["menu_F4"]["bic"], 1)}
                       for r in range(R)])
    m4["warning"] = m4.curve_minus_f4_bic > DEC["m4_dbic_vs_f4"]
    save_output(m4, P.STEP3B3B_M4, "product", "step 3b-3b m4 per replicate: BIC(free curve) - BIC(F4), the warning sign",
                index=False)

    # ------------------------------------------------------------------------------------------- verdicts
    dec = []
    for truth in DEC["shape"]:
        a_ = cs[(cs.kind == "main") & (cs.truth == truth)].iloc[0]
        b_ = j2[(j2.kind == "main") & (j2.truth == truth)].iloc[0]
        oks = [bool(a_.curve_ok), bool(a_.shape_ok), bool(b_.equation_ok)]
        dec.append({"truth": truth, "curve_ok": oks[0], "shape_ok": oks[1], "equation_ok": oks[2],
                    "verdict": "PASS" if all(oks) else ("FAIL" if not any(oks) else "PARTIAL"),
                    "basis": f"coverage {a_.band_coverage:.0%}, level x{a_.median_k_ratio:.2f}; {a_.shape_wanted} in "
                             f"{a_[a_.shape_wanted + '_reps']} of {R}; equation: {b_.returned}; {b_.coef_ratios}"})
    n4 = int(m4.warning.sum())
    dec.append({"truth": "m4", "verdict": "PASS" if n4 >= need else "FAIL",
                "basis": f"BIC(free curve) - BIC(F4) > {DEC['m4_dbic_vs_f4']:g} in {n4} of {R}"})
    dec.append({"truth": "m6", "verdict": "not judged", "basis": "switch at 0.16 dB, below W_lo: cannot be named here by design"})
    dec = pd.DataFrame(dec)
    judged = dec[dec.verdict != "not judged"].verdict
    overall = "PASS" if (judged == "PASS").all() else ("FAIL" if (judged == "FAIL").all() else "PARTIAL")
    save_output(pd.concat([dec, pd.DataFrame([{"truth": "overall", "verdict": overall}])], ignore_index=True),
                P.STEP3B3B_DECISION, "product", "step 3b-3b verdict per truth and overall (rule fixed before the run)",
                index=False)

    # ------------------------------------------------------------------------------------------- figures
    for truth in S["truths"]:
        reps, w, K = store[("main", truth)]
        cover = cs[(cs.kind == "main") & (cs.truth == truth)].band_coverage.iloc[0]
        fig = k_figure(truth, reps, w, K, par, cover)
        save_output(fig, Path(str(P.STEP3B3B_FIG_K).format(truth=truth)), "diag",
                    f"step 3b-3b true vs learned k(W), replicate band and judged range, truth {truth}", dpi=110)
        plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    rng_tab = fits.groupby(["kind", "truth"]).agg(
        w_lo_db=("w_lo_db", "first"), w_hi_min=("w_hi_db", "min"), w_hi_max=("w_hi_db", "max"),
        tails=("n_tails", "median"), independent=("n_independent", "median"),
        minutes_in_range=("minutes_in_range", "median"), knots_rep0=("knots_db", "first")).reset_index()
    rng_tab.insert(4, "w_hi_judged", [cs[(cs.kind == r.kind) & (cs.truth == r.truth)].w_hi_db.iloc[0] for r in rng_tab.itertuples()])
    fit_sum = fits.groupby("kind").agg(
        fits=("rep", "size"), s_median=("s", "median"), start_spread_median=("start_cost_spread", "median"),
        minutes_median=("seconds", lambda v: round(v.median() / 60, 1))).reset_index()
    lam_tab = fits.groupby("kind")["lambda"].value_counts().unstack(fill_value=0).reset_index()
    start_tab = fits.groupby("kind").best_start.value_counts().unstack(fill_value=0).reset_index()
    curve_cols = ["truth", "w_hi_db", "band_coverage", "median_k_ratio", "curve_ok"]
    shape_cols = ["truth", "true_slope", "slopes", "slope_median", "flat_reps", "rising_reps", "shape_wanted", "shape_ok"]
    j2_cols = ["truth", "threshold", "returned", "true_selected", "spurious", "coef_ratios", "equation_ok"]
    main_cs, main_j2 = cs[cs.kind == "main"], j2[j2.kind == "main"]
    sel_piv = (sel.assign(cell=[f"{'**' if r.role == 'true' else ''}{r.n_selected}/{r.replicates}{'**' if r.role == 'true' else ''}"
                                for r in sel.itertuples()])
               .pivot(index=["kind", "truth"], columns="term", values="cell").reindex(columns=JTERMS).reset_index())
    six = main_cs[main_cs.truth == "m6"].iloc[0]
    figs = "\n\n".join(f"![k {t}](../figures/{Path(str(P.STEP3B3B_FIG_K).format(truth=t)).name})" for t in S["truths"])
    t1 = "T1" in plan["kinds"]
    report = f"""# WAA_estm step 3b-3b: confirmation run of "draw the curve, then name it" (rule fixed before the run)

Written by `research/WAA_estm/step3b3b_confirm.py` (libraries `waa_curve.py`, `waa_pooled.py`, `waa_sindy.py`,
`waa_sim.py`; settings `waa.curve3b3b`, fit settings `waa.curve3b3`). Fake data only: no real tail is fitted, the frozen
tail set v1 is not changed.

**Overall verdict: {overall}.** Per truth: {', '.join(f"{r.truth} {r.verdict}" for r in dec.itertuples())}.

## 1. What was run

- **Why**: step 3b-3 found, after seeing its results, that the free curve and its equation are right above about
  0.35 dB. This run tests exactly that on fresh data, with the range, the tests and the rule written before the run.
- **Data**: NEW event bootstraps of the {len(skel)}-tail TRAIN v1 skeleton (seeds {S['seed_base']} + replicate; 3b-3 used
  3000), {R} replicates, 3a-fix noise, speed spread s = {S['s_true']}. Truths {', '.join(S['truths'])}. Robustness: T1 noise
  (two slow time scales, the fit keeps one) for {', '.join(S['t1_truths'])}{'' if t1 else ' (CUT)'}. Window: first {S['window_min']} min.
- **Fit**: the free-curve pooled fit of 3b-3, unchanged ({len(FIT['starts'])} starts of {FIT['maxfev']['start']} evaluations,
  the best continued for {FIT['maxfev']['finish']}, penalties {FIT['lambda_grid']} chosen by 2-fold cross-validation over
  event draws). Only the knots moved: {S['knots']} knots at W_lo, W_hi and 3 equal quantiles of log A in between. Below
  W_lo the curve continues as a straight line in log-log; it is fitted there (the tail ends need it) but not judged.
  The curvature interval of 3b-3 is not computed.
- **Fits**: {len(results)} free-curve fits, {int(fits.nfev.sum())} likelihood evaluations in the main stages; plus the F4
  menu fit of the {R} m4 datasets (starts {', '.join(S['menu_f4_starts'])}, {S1C['maxfev_generic']} evaluations each).
- **Timing and cut**: replicate 0 (all {len(datasets)} datasets) was run first: {plan['wall_rep0_s'] / 60:.1f} min wall on
  {workers} processes; projected total {plan['projected_s'] / 3600:.2f} h (budget {S['time_budget_h']} h).
  Cut: {'; '.join(plan['cuts']) if plan['cuts'] else 'nothing'}.
- **Run time**: {plan['wall_s'] / 60:.0f} min wall ({fits.seconds.sum() / 3600:.1f} process-hours; cached in
  `.cache/fits3b3b/`). This run: {time.time() - clock:.0f} s.

{md(fit_sum)}

Chosen penalty (fits):

{md(lam_tab)}

Best start (fits; with all knots above 0.2 dB the "step" start is a flat start at k = {FIT['starts']['step']['high']} / min):

{md(start_tab)}

## 2. Judged range per dataset

W_lo = {W_LO} dB (fixed). W_hi = the level with {FIT['support_min_minutes']} tail-minutes of A above it in the first
{S['window_min']} min, per dataset (`w_hi_min`-`w_hi_max` over the replicates). Band, slopes and jump 2 use
[W_lo, `w_hi_judged`] = the range common to the {R} replicates, on {FIT['grid_points']} log-spaced points.
`minutes_in_range` = tail-minutes with A inside the range (median over replicates).

{md(rng_tab)}

## 3. Curve (`{rel(P.STEP3B3B_CURVE_SHAPE)}`)

Rule: the true k(W) inside the min-max band of the {R} replicates on >= {DEC['band_coverage']:.0%} of the judged range,
AND the median learned / true k within +-{DEC['level_tol']:.0%}.

{md(main_cs[main_cs.truth.isin(DEC['shape'])][curve_cols])}

{figs}

## 4. Shape

Slope of log k vs log W over the judged range, one number per replicate (least squares on the grid). Rule: m1 flat =
|slope| < {DEC['flat_abs_slope']} in >= {need} of {R}; m3, m5 rising = slope > {DEC['rising_slope']} in >= {need} of {R}.
m4 and m6 are listed for comparison (no rule).

{md(main_cs[shape_cols])}

## 5. Jump 2: terms and coefficients (`{rel(P.STEP3B3B_J2_SEL)}`, `{rel(P.STEP3B3B_J2_COEF)}`)

h(W) = W k(W) of every replicate on the judged grid; library {', '.join(JTERMS)} (sw x = max(W - x dB, 0)); the
thresholded backward elimination and threshold choice of 3b-3 (leave-one-replicate-out, one-SE rule). Replicates in
which a term is selected (**bold** = a term of the true law):

{md(sel_piv)}

Rule: true terms in >= {need} of {R}, no other term in more than {DEC['max_spurious']} of {R}, counted coefficient
(median of fitted / true, 0 where not selected) within +-{DEC['coef_tol']:.0%}.

{md(main_j2[main_j2.truth.isin(DEC['j2_terms'])][j2_cols])}

{md(coefs[coefs.kind == 'main'].drop(columns='kind'))}

## 6. m6 and m4

**m6 (state switch at {par['m6']['A_c']:.2f} dB, below W_lo: not judged).** Above A_c the true law is
k = k_evap + k_run (1 - A_c / W): it rises toward k_run. Over the judged range: true slope {six.true_slope}, learned
slopes {six.slopes} (median {six.slope_median}); true k inside the band on {six.band_coverage:.0%} of the range, level
x{six.median_k_ratio}. Jump 2 returns: {main_j2[main_j2.truth == 'm6'].returned.iloc[0]} (median coefficients:
{main_j2[main_j2.truth == 'm6'].coefficients_median.iloc[0]}).

**m4 (two layers: no single curve).** Rule: PASS if BIC(free curve) - BIC(F4) > {DEC['m4_dbic_vs_f4']:g} in >= {need} of {R}.
The free curve counts {S['knots'] + 1} parameters ({S['knots']} knots + s), F4 counts 4.

{md(m4)}

What the free curve draws for m4 over the judged range: slopes {main_cs[main_cs.truth == 'm4'].slopes.iloc[0]}; jump 2
returns {main_j2[main_j2.truth == 'm4'].returned.iloc[0]}.

## 7. T1 robustness (noise the fit does not know)

{('Same tests as the main data; reported, not part of the verdict.' + chr(10) + chr(10) + md(cs[cs.kind == 'T1'][curve_cols + shape_cols[1:]]) + chr(10) + chr(10) + md(j2[j2.kind == 'T1'][j2_cols])) if t1 else 'Not run (cut for time).'}

Fitted spread s (true {S['s_true']}), median: {', '.join(f"{k} {t}: {v:.2f}" for (k, t), v in fits.groupby(['kind', 'truth']).s.median().items())}.

## 8. Verdict and what it means for Step 4

{md(dec)}

**Overall: {overall}** (over m1, m3, m5, m4; PASS = every criterion of every truth).

{MEANING}

## 9. Unclear / possibly wrong

{UNCLEAR}
"""
    save_output(report, P.STEP3B3B_REPORT, "product", "step 3b-3b report: confirmation run of the free drying curve and "
                "its equation above 0.35 dB on fresh fake tails, rule fixed before the run")
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {time.time() - clock:.0f} s")


# Written after the run (the tables are computed on every run; these two texts are the reading of that run and quote
# its numbers: update them if a rerun changes the tables).
MEANING = """**In plain words: the post hoc finding of 3b-3 is confirmed on fresh data. Above 0.35 dB the free curve and its
equation are right for the three single-curve laws, on every criterion. The m4 warning sign is NOT confirmed (2 of 5),
which is why the overall verdict is PARTIAL.**

- **m1, m3, m5: PASS on curve, shape and equation.** The true k(W) is inside the replicate band on 100 % of the
  judged range, the level is within 3 % (x0.99, x0.97, x1.03). The slopes separate the laws cleanly: m1 between
  -0.01 and 0.05 (true 0), m3 between 0.99 and 1.12 (true 1), m5 between 0.93 and 1.03 (true 0.92). Jump 2 returns
  exactly one term in 5 of 5: "W" for m1 (coefficient x0.99), "W^2" for m3 (x0.97) and for m5 (x1.13; m5's small
  linear part, 0.003 / min, is never selected and goes into the W^2 coefficient). The mixing trap (speed spread
  making an exponential look nonlinear) does not appear: m1 gets no false term in any replicate.
- **T1 noise (noise the fit does not know)**: same slopes and the same one-term equations (W x1.05, W^2 x0.96). The
  only miss is the band coverage of m3 (72 %, level x0.97): the band of 5 replicates is narrow there.
- **m6 (switch at 0.16 dB, below W_lo).** Over the judged range the curve is right (inside the band on 92 %, level
  x0.97, slopes 0.20-0.32 against a true 0.23): a mild rise, clearly neither flat (m1) nor W^2 (m3, m5). Jump 2 gives
  a three-term fit (W^2 + W^3 + a switch term at 0.5 dB with negative coefficients) in 5 of 5, which is a
  description of a bent curve, not a law. A real-data user would conclude: "above 0.35 dB the rate rises slowly
  with W (slope about 0.3); not exponential, not quadratic; the equation is not a clean one-term law". That is
  true, but it does not say "switch": the switch itself is below the judged range.
- **m4 (two layers): FAIL of the warning sign.** The free curve is worse than F4 by more than 10 BIC in only 2 of 5
  replicates (2.7, 4.2, 38.7, 53.2, 3.3). In 3b-3 it was 4 of 5. So "free curve worse than F4" cannot be relied on
  to expose a hidden second state in the first 60 min. What m4 does show over the judged range: slopes 0.36-0.47
  and a three-term equation, i.e. the same picture as m6. The curve cannot tell m4 from m6.
- **The slope is the most useful single number.** About 0 = exponential; about 1 = quadratic (m3 / m5, which stay
  one family as in 3b-1c); 0.2-0.5 with a multi-term equation = "something else" (switch or two layers), to be
  settled by the menu fits F6 and F4, not by the curve.

**What this means for Step 4 (real tails).**
1. Run exactly this setup on the real tails: first 60 min, 5 knots from 0.35 dB to W_hi, band from an event
   bootstrap, slope of log k vs log W, jump 2 with this library.
2. **The statement we can make on real data holds for W between 0.35 dB and W_hi (about 2 dB here; it must be
   computed on the real tails) and for the first 60 min only**: "in this range the drying rate per dB of water is
   flat / rises with slope x, and h(W) is (or is not) a single-term law". If the real slope is near 0 with equation
   "W", the exponential is supported in that range; near 1 with "W^2", a quadratic law. This run shows both would
   be recovered with the level within about 5 % and the leading coefficient within about 15 %.
3. A slope of 0.2-0.5 with a multi-term equation must be reported as "mildly rising, law not identified by the
   curve". It is what both a switch and two layers look like. Nothing may be said about W < 0.35 dB (so nothing
   about a switch level near 0.16 dB) or about the late, slow part of the tails from this tool.
4. Keep the menu fits (F1, F3 / F5, F4, F6) next to the curve, with the 3b-1c rules. Do not use "free curve vs F4"
   as a pass / fail check of a hidden state; report the BIC difference only as information."""
UNCLEAR = """1. **The F4 fit behind the m4 test may be too weak, and I did not check.** I gave F4 three generic starts (g1, g2,
   g3; 800 evaluations each); 3b-3 took the best of four (it also had the start built from the F1 solution). The
   best start was g3 in all 5. If F4 did not reach its optimum, its BIC is too high and the difference too small,
   which would turn a real warning into a FAIL. Against that: all 5 fits stopped on tolerance after 157-322
   evaluations and landed on similar parameters (tau_f 16-18 min, tau_s 156-295 min, phi 0.16-0.20). The more
   likely reading is that the new knots (all above 0.35 dB, a free straight line below) let one curve imitate two
   layers better than the 3b-3 knots did. I did not refit: the rule says no re-running after a failed cell. A check
   with the F1-solution start on these 5 cached datasets would settle it and takes minutes.
2. **F4 start choice was mine**, not in the task text (which only says BIC(F4)); it was fixed in the config before
   the run.
3. **"Same thresholded method and threshold choice" for jump 2**: with the smaller library the cross-validation
   picked thresholds 0.3-0.45, the top of the grid. A high threshold favours one-term answers, which helps m1, m3,
   m5 and would hide a small second term (m5's linear part is never found; a weak real second term would not be
   either).
4. **m5 and m3 are not told apart** (both "W^2", slopes 0.96 vs 1.05). This is the F3 / F5 family of 3b-1c again; the
   rule counted only W^2 for m5, so it passes.
5. **The band is the min-max of 5 replicates of the SAME skeleton** (event bootstrap). 100 % coverage with 5
   curves is not a calibrated interval. On real data there is one dataset, so the band must come from an event
   bootstrap of the fit (5 fits of about 3 min each); its coverage was not tested here, since here the replicates
   also had new noise and new speed factors.
6. **The straight line below W_lo is fitted, not judged, but it is not harmless**: the tails spend most of their
   minutes below 0.35 dB, and those minutes pull on the slope at the first knot. For m1, m3, m5 the truth happens to
   be a straight line in log-log there, so the continuation is exactly right. For m6 (a step at 0.16 dB) and m4 it
   is wrong, and the curve just above W_lo may bend to compensate (m4: steep between 0.35 and 0.5 dB, then a
   plateau). Real tails may be of that kind.
7. **The fitted spread s is low** (median 1.14-1.21 against a true 1.25) in every dataset. The same sign as in the
   earlier steps; the curve level is nevertheless within 3 %.
8. **The search budget is unchanged from 3b-3** (150 evaluations per start). The starts still differ by a median
   of 8 in cost after the first stage (3b-3: 34-41), so the problem is much better conditioned, but I did not
   verify that the final fits are at the optimum.
9. **W_hi differs between replicates** (1.9-3.5 dB) because it follows the bootstrap draw; the judged range is the
   common part, so the top knot of most replicates lies above it. On real tails there is one W_hi.
10. **Timing**: replicate 0 was run and timed in full this time (7.7 min, projected 0.64 h); the run took 32 min,
    nothing was cut.
11. **Key figure**: the limit of 6 is reached, so the five `step3b3b_k_<truth>.png` are diag. If one should become
    key I would take `step3b3b_k_m1.png` (flat curve inside the band despite the spread) and demote
    `outputs/tables/step3b1b_recovery.csv`. The real-tail curve of Step 4 is the better candidate, so I would wait."""


if __name__ == "__main__":
    main()
