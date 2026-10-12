"""WAA_estm step 3b-3: "draw the curve, then name it" on fake tails (known truth). No real tail is fitted.

Two jumps instead of SINDy on raw tails (which failed in step 3b-2):
  Jump 1  the pooled fit of steps 3b-1b / 3b-1c (offset prior, ARMA(1,1) noise, per-tail speed marginalised), but the
          drying law is a FREE smooth curve learned from all tails together: dW/dt = -f_i W k(W), log k = natural
          cubic spline in log W (library waa_curve.py). k(W) is the main output: flat = exponential, rising =
          nonlinear in W, a step = state switch.
  Jump 2  sparse regression on that clean curve h(W) = W k(W) (not on raw tails) to name its equation.
Question: with a known truth, does jump 1 draw the right curve and does jump 2 name the right equation?

Data: the fake datasets of steps 3b-1c / 3b-2 (event bootstrap, seeds 3000 + replicate, 5 replicates): truths m1, m3,
m4, m5, m6 with spread s = 1.25 and 3a-fix noise (main test); extras: m1 without spread, and m1, m3 with the T1 noise
of step 3b-1c. Windows: first 60 min and first 6 h. The decision rule is fixed in `waa.curve3b3.decision` before the
run. The only uses of real data: the NaN minutes of the real tails (copied into the fake tails) and the dry-noise
statistics in sim_noise_params.csv (T1 noise, as in 3b-1c). The frozen tail set v1 is read, never changed.

Per dataset and window (one job; jobs run in parallel processes, cores - 1):
  1. four starts (flat, rising, falling, step) with the middle penalty; the best is continued;
  2. the other penalties, from that solution; 2-fold cross-validation over EVENT DRAWS picks the penalty;
  3. the curvature of the fit at the chosen solution gives a per-replicate uncertainty of the curve.
Every stage of every job is cached in .cache/fits3b3/, so a restart resumes. Delete that folder to redo.
Run from the repository root (about 2.5 hours the first time, 1 minute when cached):
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step3b3_curve_then_equation.py
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
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MODULE = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE))
import WAA_estm_paths as P
from research.common.save import save
import waa_sim as ws
import waa_pooled as wp
import waa_curve as wc
import step3b1c_stress_test as st                          # the same event bootstrap, noise tables, T1 noise, starts

CFG = yaml.safe_load(open(P.CONFIG_FILE, encoding="utf-8"))
SIM, S, S1C = CFG["waa"]["sim"], CFG["waa"]["curve3b3"], CFG["waa"]["stress3b1c"]
STEP, SPM = SIM["rsl_step_db"], SIM["samples_per_min"]
DEC = S["decision"]
JGRID = S["j2_switch_grid_db"]
JTERMS = wc.curve_terms(JGRID)
KINDS = {"main": "spread 1.25, 3a-fix noise", "nospread": "no spread, 3a-fix noise", "T1": "spread 1.25, T1 noise"}
INK, INK2 = "#1a1a1a", "#52514e"
CACHE = P.CACHE / "fits3b3"
MENU_CACHE = P.CACHE / "sindy3b2"                          # 60 min menu fits (best of 4 starts) of steps 3b-1c / 3b-2
N_EVAL_COV = 1 + 2 * (S["knots"] + 1) + 2 * (S["knots"] + 1) * S["knots"] // 2   # evaluations of Curve.covariance
N_EVAL = (len(S["starts"]) * S["maxfev"]["start"] + S["maxfev"]["finish"]
          + (len(S["lambda_grid"]) - 1) * S["maxfev"]["other_lambda"]
          + len(S["lambda_grid"]) * S["maxfev"]["cv"] + N_EVAL_COV)   # planned (the 2 CV halves = 1 full dataset)


def rel(p):
    return Path(p).resolve().relative_to(P.REPO).as_posix()


def md(df, index=False):
    df = df.reset_index() if index else df
    head = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns)
    return head + "\n" + "\n".join("| " + " | ".join("" if pd.isna(x) else str(x) for x in row) + " |"
                                   for row in df.itertuples(index=False))


def true_k(truth, par, w):
    """k(W) = -(dW/dt) / W of the base law (speed factor 1); None for m4 (two states: not a function of W)."""
    w = np.asarray(w, float)
    p = par[truth]
    if truth == "m1":
        return np.full_like(w, 1 / p["tau"])
    if truth == "m3":
        return p["c"] * w ** (p["n"] - 1)
    if truth == "m5":
        return p["a"] + p["b"] * w
    if truth == "m6":
        return np.where(w > p["A_c"], p["k_evap"] + p["k_run"] * (1 - p["A_c"] / w), p["k_evap"])
    return None


def start_theta(spec, knots):
    if "slope" in spec:
        return np.log(spec["k1"]) + spec["slope"] * np.log(knots)
    return np.log(np.where(knots <= spec["at_db"], spec["low"], spec["high"]))


# =================================================================================================== one job (worker)
_INPUTS = {}


def load_tails(sim):
    rows = pd.read_parquet(sim).sort_values(["tail_id", "minutes_since_t0_trace"])
    return [{"tail_id": t, "sublink_id": g.sublink_id.iloc[0], "cml_id": g.cml_id.iloc[0],
             "event_copy": int(g.event_copy.iloc[0]), "A_start": g.A_start_trace_db.iloc[0],
             "t": g.minutes_since_t0_trace.to_numpy(), "a": g.A_db.to_numpy(float)} for t, g in rows.groupby("tail_id")]


def menu_fit(job, D, n_ind, cand, labels, res, save):
    """Menu candidate `cand` from the given starts of step 3b-1c (best likelihood). Cached in the job's result."""
    key = f"menu_{cand}"
    if key not in res:
        best = None
        for lab in labels:
            th0, s0 = S1C["t2_starts"][cand][lab]
            r = wp.fit(cand, D, th0, s0, S1C["maxfev_generic"])
            if best is None or r["ll"] > best["ll"]:
                best = {**r, "start": lab}
        best["bic"] = -2 * best["ll"] + (len(wp.NAMES[cand]) + 1) * np.log(n_ind)
        res[key] = best
        save()
    return res[key]


def run_job(job):
    """Free-curve pooled fit of ONE dataset in ONE window. Stages are cached in one json file per job."""
    f_ = CACHE / f"job_{job['tag']}_{job['window']}.json"
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
    tails = load_tails(job["sim"])
    D = wp.pack(tails, noise, offset_sd, job["minutes"], S["min_points"])
    used = set(D["tail_id"])
    n_ind = len({(tl["cml_id"], tl["event_copy"]) for tl in tails if tl["tail_id"] in used})
    cv = wc.Curve(job["knots"], prior_sd=S["prior_sd_logk"])
    lam0, mf = S["lambda_start"], S["maxfev"]
    res.update({k: job[k] for k in ("tag", "truth", "kind", "rep", "window", "knots", "w_hi")})
    res.update({"n_tails": len(D["tail_id"]), "n_independent": n_ind, "n_points": int(len(D["T"]))})

    # 1. the starts with the middle penalty; the best one is continued
    res.setdefault("starts", {})
    for name, spec in S["starts"].items():
        if name not in res["starts"]:
            res["starts"][name] = cv.fit(D, np.r_[start_theta(spec, cv.knots), np.log(S["s_start"])], lam0, mf["start"])
            save()
    res.setdefault("lam", {})
    if str(lam0) not in res["lam"]:
        best = min(res["starts"], key=lambda n: res["starts"][n]["cost"])
        res["lam"][str(lam0)] = {**cv.fit(D, res["starts"][best]["x"], lam0, mf["finish"]), "from": best}
        save()
    x_mid = res["lam"][str(lam0)]["x"]
    # 2. the other penalties, from that solution
    for lam in S["lambda_grid"]:
        if str(lam) not in res["lam"]:
            res["lam"][str(lam)] = cv.fit(D, x_mid, lam, mf["other_lambda"])
            save()
    # 3. 2-fold cross-validation over event draws: fit on one half, log-likelihood of the other half
    res.setdefault("cv", {})
    evs = np.random.default_rng([S["seed_base"], job["rep"]]).permutation(sorted({tl["event_copy"] for tl in tails}))
    half = [set(evs[::2].tolist()), set(evs[1::2].tolist())]
    Dh = None
    for lam in S["lambda_grid"]:
        if str(lam) not in res["cv"]:
            if Dh is None:
                Dh = [wp.pack([tl for tl in tails if tl["event_copy"] in h], noise, offset_sd, job["minutes"],
                              S["min_points"]) for h in half]
            out = 0.0
            for a, b in ((0, 1), (1, 0)):
                r = cv.fit(Dh[a], res["lam"][str(lam)]["x"], lam, mf["cv"])
                out += cv.loglik(np.asarray(r["x"]), Dh[b])
            res["cv"][str(lam)] = out
            save()
    lam_star = max(S["lambda_grid"], key=lambda v: res["cv"][str(v)])
    fit = res["lam"][str(lam_star)]
    x = np.asarray(fit["x"])
    # 4. uncertainty of the chosen solution from the curvature of the penalised fit
    if "cov" not in res:
        C = cv.covariance(x, D, lam_star)
        res["cov"] = None if C is None else C.tolist()
        save()
    res.update({"lam_star": lam_star, "x": fit["x"], "ll": fit["ll"], "s": float(np.exp(x[-1])),
                "bic": -2 * fit["ll"] + (cv.K + 1) * np.log(n_ind),
                "nfev": sum(r["nfev"] for r in list(res["starts"].values()) + list(res["lam"].values()))})
    # 5. F4 menu fit in the 6 h window (m4 only): for the "free curve worse than F4" warning
    if job["menu_f4"]:
        menu_fit(job, D, n_ind, "F4", S["menu_6h_starts"], res, save)
    res["seconds"] = round(res.get("seconds", 0.0) + time.time() - t0, 1)
    res["done"] = True
    save()
    return res


# =================================================================================================== figure
def k_figure(truth, cur, bands, par, R):
    """True k(W) (dashed) vs learned k(W) (solid) with the band of the replicates (hatched), one panel per window."""
    wins = list(S["windows"])
    fig, ax = plt.subplots(1, len(wins), figsize=(6.6 * len(wins), 4.9), sharey=True)
    for a_, wname in zip(np.atleast_1d(ax), wins):
        g = cur[(cur.truth == truth) & (cur.kind == "main") & (cur.window == wname)]
        if g.empty:
            a_.set_title(f"{wname}: not run")
            continue
        piv = g.pivot(index="w_db", columns="rep", values="k")
        w = piv.index.to_numpy()
        a_.fill_between(w, piv.min(axis=1), piv.max(axis=1), facecolor="none", edgecolor="#8a8a8a", hatch="///", lw=0,
                        label=f"band: min-max of the {R} replicates")
        for r in piv.columns:
            a_.plot(w, piv[r], "-", color=INK2, lw=0.7)
        a_.plot(w, np.exp(np.log(piv).median(axis=1)), "-", color=INK, lw=2.4, label="learned k(W), median (thin: each replicate)")
        kt = true_k(truth, par, w)
        if kt is not None:
            a_.plot(w, kt, "--", color=INK, lw=2.0, label="true k(W)")
            a_.text(w[-1], kt[-1], " true", fontsize=9, va="center")
        else:                                              # m4: the rate depends on where the tail started
            for a0, ls in ((1.0, ":"), (3.0, "-.")):
                t = np.linspace(0, 720, 4000)
                wt = ws.law_curve("m4", par["m4"], t, a0)
                a_.plot(wt[1:], -np.diff(wt) / np.diff(t) / wt[1:], ls, color=INK, lw=1.6,
                        label=f"true rate along a tail that starts at {a0:g} dB")
        a_.text(w[-1], np.exp(np.log(piv).median(axis=1)).iloc[-1], " learned", fontsize=9, va="center")
        b = bands[(bands.truth == truth) & (bands.kind == "main") & (bands.window == wname)]
        a_.set(xscale="log", yscale="log", xlabel="W (dB), supported range", xlim=(w[0], w[-1] * 1.35),
               title=f"{wname}" + (f": true k inside the band on {b.band_coverage.iloc[0]:.0%} of the range"
                                   if len(b) and np.isfinite(b.band_coverage.iloc[0]) else ""))
        a_.grid(alpha=0.25, which="both")
        for sp in ("top", "right"):
            a_.spines[sp].set_visible(False)
    np.atleast_1d(ax)[0].set_ylabel("k(W) = drying rate per dB of water (1 / min)")
    np.atleast_1d(ax)[0].legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.suptitle(f"Truth {truth} ({ws.LAWS[truth]}), spread 1.25, 3a-fix noise: free-curve pooled fit.\n"
                 "Flat k = exponential, rising k = nonlinear in W, a step in k = state switch", fontsize=10)
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    clock = time.time()
    CACHE.mkdir(parents=True, exist_ok=True)
    workers = max((os.cpu_count() or 2) - 1, 1)
    R0 = S["replicates"]

    # ------------------------------------------------------------------------------------------- inputs and truths
    skel = pd.read_csv(P.TAILS_MODEL_INDEX).query("split_event == 'train'").sort_values("tail_id").reset_index(drop=True)
    real = pq.read_table(P.TAILS_MODEL, filters=[("split_event", "=", "train")],
                         columns=["tail_id", "minutes_since_t0_trace", "A_db"]).to_pandas()
    nan_real = {t: g.sort_values("minutes_since_t0_trace").A_db.isna().to_numpy() for t, g in real.groupby("tail_id")}
    del real
    nz, noise, offset_sd = st.load_noise()
    cal = pd.read_csv(P.STEP3A_CALIBRATION).set_index("law")
    par = {law: {k: float(cal.loc[law, k]) for k in ws.PARAMS[law]} for law in ("m1", "m3", "m4", "m5")}
    m6 = pd.read_csv(P.STEP3B1_M6)
    b6 = m6.loc[m6.rmse_band_observed.idxmin()]
    par["m6"] = {k: float(b6[k]) for k in ws.PARAMS["m6"]}
    true_coef = {"m1": {"W": 1 / par["m1"]["tau"]}, "m3": {"W^2": par["m3"]["c"]},
                 "m5": {"W": par["m5"]["a"], "W^2": par["m5"]["b"]},
                 "m6": {"W": par["m6"]["k_evap"], "switch": par["m6"]["k_run"]}}
    datasets = ([("main", t) for t in S["truths"]] + [("nospread", t) for t in S["nospread_truths"]]
                + [("T1", t) for t in S["t1_truths"]])
    print(f"{len(skel)} skeleton tails; datasets: " + ", ".join(f"{k} {t}" for k, t in datasets) + f"; {workers} processes")

    # ------------------------------------------------------------------------------------------- simulate (cached)
    sim_file = lambda kind, t, rep: CACHE / f"sim_{kind}_{t}_{rep}.parquet"   # noqa: E731
    noise_t1 = None
    meta = {}
    for rep in range(R0):
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
                f = np.ones(len(sk)) if kind == "nospread" else np.exp(S["s_true"] * z)
                rows, _ = ws.simulate(sk, t, par[t], noise_t1 if kind == "T1" else noise, SIM["switches"], seed, f,
                                      STEP, SPM, masks)
                rows.merge(sk[["tail_id", "event_copy", "A_start_trace_db"]], on="tail_id").to_parquet(f_)
            # supported range and knots from the first 60 min (the same knots for both windows of a dataset)
            a = pd.read_parquet(f_, columns=["minutes_since_t0_trace", "A_db"])
            n_by_win = {w: int(a.A_db[a.minutes_since_t0_trace < m].notna().sum()) for w, m in S["windows"].items()}
            a = a.A_db[(a.minutes_since_t0_trace < 60) & a.A_db.notna()].to_numpy()
            w_hi = float(np.sort(a)[-S["support_min_minutes"]])
            inside = np.log(a[(a >= S["w_floor_db"]) & (a <= w_hi)])
            q = np.quantile(inside, np.linspace(0, 1, S["knots"])[1:-1])
            meta[(kind, t, rep)] = {"w_hi": w_hi, "points": n_by_win,
                                    "knots": [S["w_floor_db"]] + [float(v) for v in np.exp(q)] + [w_hi]}
    print(f"datasets ready ({time.time() - clock:.0f} s); supported range {S['w_floor_db']} dB to "
          f"{min(m['w_hi'] for m in meta.values()):.2f}-{max(m['w_hi'] for m in meta.values()):.2f} dB")

    def jobs_for(replicates, kinds, extras_6h):
        jobs = []
        for (kind, t, rep), m in meta.items():
            if rep >= replicates or kind not in kinds:
                continue
            for wname, minutes in S["windows"].items():
                if kind != "main" and wname != "first_60min" and not extras_6h:
                    continue
                jobs.append({"tag": f"{kind}_{t}_{rep}", "truth": t, "kind": kind, "rep": rep, "window": wname,
                             "minutes": minutes, "sim": str(sim_file(kind, t, rep)), "knots": m["knots"],
                             "w_hi": m["w_hi"], "points": m["points"][wname],
                             "menu_f4": kind == "main" and t == "m4" and wname != "first_60min"})
        return sorted(jobs, key=lambda j: -j["points"])    # the long jobs first

    # ------------------------------------------------------------------------------------------- time, then plan
    plan_file = CACHE / "plan.json"
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
        print(f"plan from {rel(plan_file)}: cut: {plan['cuts'] or 'nothing'}")
    else:
        print("\ntiming one replicate (m6, replicate 0): likelihood evaluations in both windows ...", flush=True)
        tails = load_tails(sim_file("main", "m6", 0))
        cvt = wc.Curve(meta[("main", "m6", 0)]["knots"])
        x = np.r_[start_theta(S["starts"]["flat"], cvt.knots), 0.0]
        spp = []
        for wname, minutes in S["windows"].items():
            D = wp.pack(tails, noise, offset_sd, minutes, S["min_points"])
            cvt.loglik(x, D)                               # compile the numba code before the clock starts
            t0 = time.time()
            for i in range(8):
                cvt.loglik(x + 0.01 * i, D)
            spp.append((time.time() - t0) / 8 / len(D["T"]))
            print(f"   {wname}: {len(D['tail_id'])} tails, {len(D['T'])} points, {spp[-1] * len(D['T']):.3f} s per evaluation")
        spp = float(np.mean(spp))

        def projected(replicates, kinds, extras_6h):
            jobs = jobs_for(replicates, kinds, extras_6h)
            return sum(N_EVAL * j["points"] * spp * (1.5 if j["menu_f4"] else 1.0) for j in jobs) / workers

        plan = {"replicates": R0, "kinds": list(KINDS), "extras_6h": True, "cuts": [], "wall_s": 0.0,
                "sec_per_point": spp, "evals_per_job": N_EVAL}
        proj = projected(R0, plan["kinds"], True)
        m60, m6h = (meta[("main", "m6", 0)]["points"][w] * spp * N_EVAL for w in S["windows"])
        print(f"   planned {N_EVAL} evaluations per dataset and window -> this replicate {m60 / 60:.0f} + {m6h / 60:.0f} "
              f"min; projected total {proj / 3600:.2f} h with {workers} processes (budget {S['time_budget_h']} h)", flush=True)
        for cut, change in (("6 h window of the extras (m1 no spread, T1 noise): they keep the first 60 min", {"extras_6h": False}),
                            ("T1 robustness", {"kinds": ["main", "nospread"]}),
                            ("m1 no-spread", {"kinds": ["main"]}),
                            ("replicates 5 -> 3", {"replicates": 3})):
            if proj <= S["time_budget_h"] * 3600:
                break
            plan.update(change)
            plan["cuts"].append(cut)
            proj = projected(plan["replicates"], plan["kinds"], plan["extras_6h"])
            print(f"   cut: {cut} -> projected {proj / 3600:.2f} h", flush=True)
        plan["projected_s"] = round(proj)
        plan_file.write_text(json.dumps(plan))
    R = plan["replicates"]
    jobs = jobs_for(R, plan["kinds"], plan["extras_6h"])

    # ------------------------------------------------------------------------------------------- run (parallel)
    done0 = sum(json.loads(f.read_text()).get("done", False) for f in CACHE.glob("job_*.json"))
    print(f"\n{len(jobs)} jobs (dataset x window; {done0} already finished) on {workers} processes", flush=True)
    t_run, results = time.time(), []
    total_pts = sum(j["points"] for j in jobs)
    done_pts = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(run_job, j): j for j in jobs}
        for k, fu in enumerate(as_completed(futs)):
            results.append(fu.result())
            done_pts += futs[fu]["points"]
            el = time.time() - t_run
            print(f"   jobs {k + 1}/{len(jobs)}: last {futs[fu]['tag']} {futs[fu]['window']} ({el:.0f} s elapsed, about "
                  f"{el / done_pts * (total_pts - done_pts):.0f} s left)", flush=True)
    if len(results) > done0:
        plan["wall_s"] = round(plan["wall_s"] + time.time() - clock, 1)
        plan_file.write_text(json.dumps(plan))

    # ------------------------------------------------------------------------------------------- per-fit table and curves
    wins = [w for w in S["windows"] if any(r["window"] == w for r in results)]
    by = {(r["kind"], r["truth"], r["rep"], r["window"]): r for r in results}
    cfgs = sorted({(r["kind"], r["truth"], r["window"]) for r in results},
                  key=lambda c: (list(KINDS).index(c[0]), c[1], wins.index(c[2])))
    z = DEC["interval_z"]
    fit_rows, cur_rows, band_rows, grids = [], [], [], {}
    for kind, truth, wname in cfgs:
        reps = [by[(kind, truth, r, wname)] for r in range(R)]
        w = np.geomspace(S["w_floor_db"], min(r["w_hi"] for r in reps), S["grid_points"])   # common supported range
        grids[(kind, truth, wname)] = w
        kt = true_k(truth, par, w)
        aw = wc.slope_weights(np.log(w))
        K = []
        for r in reps:
            cv = wc.Curve(r["knots"], prior_sd=S["prior_sd_logk"])
            x = np.asarray(r["x"])
            B = cv.basis(np.log(w))
            logk = B @ x[:-1]
            C = None if r["cov"] is None else np.asarray(r["cov"])[:-1, :-1]
            sd_logk = np.full(len(w), np.nan) if C is None else np.sqrt(np.einsum("ij,jk,ik->i", B, C, B))
            slope = float(aw @ logk)
            sd_slope = np.nan if C is None else float(np.sqrt(aw @ B @ C @ B.T @ aw))
            costs = {n: v["cost"] for n, v in r["starts"].items()}
            row = {"kind": kind, "truth": truth, "rep": r["rep"], "window": wname, "n_tails": r["n_tails"],
                   "n_independent": r["n_independent"], "w_hi_db": round(r["w_hi"], 3),
                   "knots_db": " ".join(f"{v:.3g}" for v in r["knots"]), "lambda": r["lam_star"],
                   "s": round(r["s"], 4), "ll": round(r["ll"], 2), "bic": round(r["bic"], 2),
                   "best_start": r["lam"][str(S["lambda_start"])]["from"],
                   "start_cost_spread": round(max(costs.values()) - min(costs.values()), 2),
                   "starts_within_2": int(sum(v - min(costs.values()) <= 2 for v in costs.values())),
                   "converged": bool(r["lam"][str(r["lam_star"])]["success"]), "nfev": r["nfev"],
                   "seconds": r["seconds"], "has_interval": C is not None,
                   "slope": round(slope, 4), "slope_lo": round(slope - z * sd_slope, 4),
                   "slope_hi": round(slope + z * sd_slope, 4),
                   **{f"k_knot{i + 1}": round(float(np.exp(v)), 5) for i, v in enumerate(x[:-1])}}
            if truth == "m6":
                ac = par["m6"]["A_c"]
                row["step_ratio"] = round(float(np.exp(logk[w >= 2 * ac].mean() - logk[w <= ac].mean())), 3)
            fit_rows.append(row)
            K.append(np.exp(logk))
            for j in range(len(w)):
                cur_rows.append({"kind": kind, "truth": truth, "window": wname, "rep": r["rep"], "w_db": w[j],
                                 "k": float(np.exp(logk[j])), "sd_log_k": float(sd_logk[j]),
                                 "k_true": np.nan if kt is None else float(kt[j])})
        K = np.array(K)
        cover = np.nan if kt is None else float(np.mean((kt >= K.min(axis=0)) & (kt <= K.max(axis=0))))
        band_rows.append({"kind": kind, "truth": truth, "window": wname, "band_coverage": cover,
                          "w_lo_db": w[0], "w_hi_db": round(w[-1], 3),
                          "median_k_ratio": np.nan if kt is None else round(float(np.exp(np.median(np.log(K / kt)))), 3)})
    fits = pd.DataFrame(fit_rows)
    cur = pd.DataFrame(cur_rows)
    bands = pd.DataFrame(band_rows)
    save(fits, P.STEP3B3_FITS, index=False)
    save(cur.round(6), P.STEP3B3_CURVES, index=False)

    # ------------------------------------------------------------------------------------------- jump 1 rule
    need = DEC["min_reps"] / 5 * R
    j1 = []
    for b in bands.itertuples():
        g = fits[(fits.kind == b.kind) & (fits.truth == b.truth) & (fits.window == b.window)]
        row = {"kind": b.kind, "truth": b.truth, "window": b.window, "band_coverage": b.band_coverage,
               "curve_ok": bool(b.band_coverage >= DEC["band_coverage"]) if np.isfinite(b.band_coverage) else np.nan,
               "median_k_ratio": b.median_k_ratio, "w_hi_db": b.w_hi_db, "median_slope": round(float(g.slope.median()), 3),
               "with_interval": int(g.has_interval.sum()), "shape_wanted": DEC["shape"].get(b.truth, "-")}
        flat = int((g.has_interval & (g.slope_lo <= 0) & (g.slope_hi >= 0)).sum())
        rising = int((g.has_interval & (g.slope_lo > 0)).sum())
        row.update({"flat_reps": flat, "rising_reps": rising})
        if b.truth == "m6":
            row["step_reps"] = int((g.step_ratio > DEC["step_ratio"]).sum())
            row["median_step_ratio"] = round(float(g.step_ratio.median()), 2)
        got = {"flat": flat, "rising": rising, "step": row.get("step_reps", 0)}.get(row["shape_wanted"])
        row["shape_ok"] = np.nan if got is None else bool(got >= need)
        j1.append(row)
    j1 = pd.DataFrame(j1)
    save(j1, P.STEP3B3_J1, index=False)

    # ------------------------------------------------------------------------------------------- jump 2
    sel_rows, coef_rows, j2 = [], [], []
    for kind, truth, wname in cfgs:
        w = grids[(kind, truth, wname)]
        H = np.array([cur[(cur.kind == kind) & (cur.truth == truth) & (cur.window == wname) & (cur.rep == r)].k.to_numpy() * w
                      for r in range(R)])
        lam, sel, coef = wc.name_curves(w, H, JGRID, S["j2_lambda_grid"])
        rule = DEC["j2_terms"].get(truth)
        for i, term in enumerate(JTERMS):
            ro = "-" if rule is None else ("true" if any(term in g for g in rule["required"]) else
                                           ("allowed" if term in rule["allowed"] else "spurious"))
            sel_rows.append({"kind": kind, "truth": truth, "window": wname, "term": term, "role": ro,
                             "n_selected": int(sel[:, i].sum()), "replicates": R, "threshold": lam})
        row = {"kind": kind, "truth": truth, "window": wname, "threshold": lam,
               "returned": ", ".join(f"{k} x{v}" for k, v in pd.Series(
                   [" + ".join(t for t, s_ in zip(JTERMS, s) if s_) or "(nothing)" for s in sel]).value_counts().items())}
        if rule is not None:
            ix = {t: i for i, t in enumerate(JTERMS)}
            need_n = [int(sel[:, [ix[t] for t in g]].any(axis=1).sum()) for g in rule["required"]]
            spur = {t: int(sel[:, ix[t]].sum()) for t in JTERMS if t not in rule["allowed"] and sel[:, ix[t]].any()}
            ratios = {}
            for c, tc in true_coef[truth].items():
                est = coef[:, [ix[t] for t in rule["required"][-1]]].sum(axis=1) if c == "switch" else coef[:, ix[c]]
                ratio = est / tc
                counted = c in rule["coefs"]
                coef_rows.append({"kind": kind, "truth": truth, "window": wname, "coefficient": c, "true": round(tc, 5),
                                  "median_ratio": round(float(np.median(ratio)), 3), "min_ratio": round(float(ratio.min()), 3),
                                  "max_ratio": round(float(ratio.max()), 3), "replicates_nonzero": int((est != 0).sum()),
                                  "counted": counted})
                if counted:
                    ratios[c] = float(np.median(ratio))
            row.update({"true_selected": "; ".join(f"{'/'.join(g)} {k}/{R}" for g, k in zip(rule["required"], need_n)),
                        "spurious": ", ".join(f"{t} {v}/{R}" for t, v in spur.items()),
                        "coef_ratios": ", ".join(f"{c} x{v:.2f}" for c, v in ratios.items()),
                        "true_ok": all(k >= need for k in need_n),
                        "spurious_ok": all(v <= DEC["max_spurious"] for v in spur.values()),
                        "coef_ok": all(abs(v - 1) <= DEC["coef_tol"] for v in ratios.values())})
            row["equation_ok"] = bool(row["true_ok"] and row["spurious_ok"] and row["coef_ok"])
        j2.append(row)
    sel, coefs, j2 = pd.DataFrame(sel_rows), pd.DataFrame(coef_rows), pd.DataFrame(j2)
    save(sel, P.STEP3B3_J2_SEL, index=False)
    save(coefs, P.STEP3B3_J2_COEF,
                index=False)

    # post hoc diagnostic, added after the first complete run (it changes no verdict): the same fits judged only above
    # the second knot. Between the floor (0.05 dB) and the second knot (about 0.3-0.4 dB) there is no knot, and the
    # first knot turned out to be badly determined; the curvature intervals were missing for nearly every fit, so
    # the shape is described here by the spread of the slope over the replicates instead.
    post = []
    for kind, truth, wname in cfgs:
        reps = [by[(kind, truth, r, wname)] for r in range(R)]
        w = np.geomspace(max(r["knots"][1] for r in reps), min(r["w_hi"] for r in reps), S["grid_points"])
        aw = wc.slope_weights(np.log(w))
        K = np.array([wc.Curve(r["knots"]).k(np.asarray(r["x"])[:-1], w) for r in reps])
        slopes = np.log(K) @ aw
        kt = true_k(truth, par, w)
        lam, sel2, coef2 = wc.name_curves(w, K * w, JGRID, S["j2_lambda_grid"])
        row = {"kind": kind, "truth": truth, "window": wname, "w_from_db": round(w[0], 3), "w_to_db": round(w[-1], 3),
               "band_coverage": np.nan if kt is None else round(float(np.mean((kt >= K.min(axis=0)) & (kt <= K.max(axis=0)))), 3),
               "median_k_ratio": np.nan if kt is None else round(float(np.exp(np.median(np.log(K / kt)))), 3),
               "true_slope": np.nan if kt is None else round(float(np.log(kt) @ aw), 2),
               "slope_min": round(float(slopes.min()), 2), "slope_median": round(float(np.median(slopes)), 2),
               "slope_max": round(float(slopes.max()), 2),
               "equation": ", ".join(f"{k} x{v}" for k, v in pd.Series(
                   [" + ".join(t_ for t_, s_ in zip(JTERMS, s) if s_) or "(nothing)" for s in sel2]).value_counts().items())}
        rule = DEC["j2_terms"].get(truth)
        if rule is not None:
            ix = {t_: i for i, t_ in enumerate(JTERMS)}
            row["coef_ratios"] = ", ".join(
                f"{c} x{np.median((coef2[:, [ix[t_] for t_ in rule['required'][-1]]].sum(axis=1) if c == 'switch' else coef2[:, ix[c]]) / true_coef[truth][c]):.2f}"
                for c in rule["coefs"])
        post.append(row)
    post = pd.DataFrame(post)
    save(post, P.STEP3B3_POSTHOC, index=False)

    # ------------------------------------------------------------------------------------------- menu fit, m4 warning
    W_ = DEC["m4_warning"]
    menu = []
    for truth in S["truths"]:
        for rep in range(R):
            row = {"truth": truth, "rep": rep}
            best = {}
            for c in S1C["candidates"]:                    # 60 min menu fits of steps 3b-1c / 3b-2 (best of 4 starts)
                fs = [json.loads(f.read_text()) for f in MENU_CACHE.glob(f"pfit_{truth}_{rep}_{c}_*.json")]
                if fs:
                    b = max(fs, key=lambda r: r["ll"])
                    best[c] = -2 * b["ll"] + b["k"] * np.log(b["n_independent"])
            r60 = by.get(("main", truth, rep, "first_60min"))
            if best and r60:
                order = sorted(best, key=best.get)
                row.update({"menu_winner_60min": order[0], "menu_dbic_to_second": round(best[order[1]] - best[order[0]], 1),
                            "curve_minus_best_menu_bic_60min": round(r60["bic"] - best[order[0]], 1),
                            "curve_minus_f4_bic_60min": round(r60["bic"] - best["F4"], 1)})
            r6h = by.get(("main", truth, rep, "first_6h"))
            if r6h and "menu_F4" in r6h:
                row["curve_minus_f4_bic_6h"] = round(r6h["bic"] - r6h["menu_F4"]["bic"], 1)
            if r60 and r6h:                                # do the two windows draw the same curve?
                a_, b_ = (fits[(fits.kind == "main") & (fits.truth == truth) & (fits.rep == rep) & (fits.window == w_)].iloc[0]
                          for w_ in ("first_60min", "first_6h"))
                wg = np.geomspace(S["w_floor_db"], min(r60["w_hi"], r6h["w_hi"]), S["grid_points"])
                lo_hi = []
                for r_ in (r60, r6h):
                    cv = wc.Curve(r_["knots"], prior_sd=S["prior_sd_logk"])
                    B = cv.basis(np.log(wg))
                    lk = B @ np.asarray(r_["x"])[:-1]
                    sdk = (np.full(len(wg), np.nan) if r_["cov"] is None else
                           np.sqrt(np.einsum("ij,jk,ik->i", B, np.asarray(r_["cov"])[:-1, :-1], B)))
                    lo_hi.append((lk - z * sdk, lk + z * sdk, lk))
                apart = (lo_hi[0][0] > lo_hi[1][1]) | (lo_hi[1][0] > lo_hi[0][1])
                row.update({"both_intervals": bool(a_.has_interval and b_.has_interval),
                            "windows_disagree_share": round(float(np.mean(apart)), 3),
                            "median_abs_log_ratio_6h_60min": round(float(np.median(np.abs(lo_hi[1][2] - lo_hi[0][2]))), 3)})
                row["warn_windows"] = bool(row["windows_disagree_share"] > W_["disagree_share"])
            d = [row.get("curve_minus_f4_bic_60min", np.nan), row.get("curve_minus_f4_bic_6h", np.nan)]
            row["warn_f4"] = bool(np.nanmax(d) > W_["dbic_vs_f4"]) if np.isfinite(d).any() else False
            row["warning"] = bool(row.get("warn_windows", False) or row["warn_f4"])
            menu.append(row)
    menu = pd.DataFrame(menu)
    save(menu, P.STEP3B3_MENU, index=False)

    # ------------------------------------------------------------------------------------------- verdicts
    dec = []
    for truth in S["truths"]:
        if truth == "m4":
            n = int(menu[menu.truth == "m4"].warning.sum())
            dec.append({"truth": "m4", "window_used": "both", "curve_ok": np.nan, "shape_ok": np.nan, "equation_ok": np.nan,
                        "verdict": "PASS" if n >= need else "FAIL", "basis": f"warning sign in {n} of {R} replicates"})
            continue
        cand = []
        for wname in wins:
            a_ = j1[(j1.kind == "main") & (j1.truth == truth) & (j1.window == wname)]
            b_ = j2[(j2.kind == "main") & (j2.truth == truth) & (j2.window == wname)]
            if len(a_):
                oks = [bool(a_.curve_ok.iloc[0]), bool(a_.shape_ok.iloc[0]), bool(b_.equation_ok.iloc[0])]
                cand.append((sum(oks), wname, oks))
        n_ok, wname, oks = max(cand, key=lambda c: c[0])   # the better window (the first one on a tie)
        dec.append({"truth": truth, "window_used": wname, "curve_ok": oks[0], "shape_ok": oks[1], "equation_ok": oks[2],
                    "verdict": "PASS" if n_ok == 3 else ("FAIL" if n_ok == 0 else "PARTIAL"),
                    "basis": "; ".join(f"{w_}: curve {o[0]}, shape {o[1]}, equation {o[2]}" for _, w_, o in cand)})
    dec = pd.DataFrame(dec)
    overall = "PASS" if (dec.verdict == "PASS").all() else ("FAIL" if (dec.verdict == "FAIL").all() else "PARTIAL")
    save(pd.concat([dec, pd.DataFrame([{"truth": "overall", "verdict": overall}])], ignore_index=True),
                P.STEP3B3_DECISION, index=False)

    # ------------------------------------------------------------------------------------------- figures
    for truth in S["truths"]:
        fig = k_figure(truth, cur, bands, par, R)
        save(fig, Path(str(P.STEP3B3_FIG_K).format(truth=truth)), dpi=110)
        plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    def pivot_sel(kind):
        g = sel[sel.kind == kind].assign(cell=lambda d: [f"{'**' if r.role == 'true' else ''}{r.n_selected}/{r.replicates}"
                                                         f"{'**' if r.role == 'true' else ''}" for r in d.itertuples()])
        return g.pivot(index=["truth", "window"], columns="term", values="cell").reindex(columns=JTERMS).reset_index()

    fit_sum = fits.groupby(["kind", "window"]).agg(
        jobs=("rep", "size"), tails=("n_tails", "median"), s_median=("s", "median"),
        converged=("converged", "sum"), with_interval=("has_interval", "sum"),
        start_spread_median=("start_cost_spread", "median"), starts_within_2_median=("starts_within_2", "median"),
        minutes_median=("seconds", lambda v: round(v.median() / 60, 1))).reset_index()
    lam_tab = fits.groupby(["kind", "window"])["lambda"].value_counts().unstack(fill_value=0).reset_index()
    start_tab = fits.groupby(["kind", "window"]).best_start.value_counts().unstack(fill_value=0).reset_index()
    j1_show = j1[["kind", "truth", "window", "band_coverage", "curve_ok", "median_k_ratio", "median_slope", "with_interval",
                  "flat_reps", "rising_reps"] + [c for c in ("step_reps", "median_step_ratio") if c in j1] +
                 ["shape_wanted", "shape_ok"]].round(3)
    j2_cols = ["truth", "window", "threshold", "returned", "true_selected", "spurious", "coef_ratios", "equation_ok"]
    m4 = menu[menu.truth == "m4"]
    warn_tab = menu.groupby("truth").agg(
        replicates=("rep", "size"), menu_winners_60min=("menu_winner_60min", lambda v: ", ".join(
            f"{k} x{n}" for k, n in v.value_counts().items())),
        curve_minus_best_menu_bic=("curve_minus_best_menu_bic_60min", "median"),
        curve_minus_f4_bic_60min=("curve_minus_f4_bic_60min", "median"),
        windows_disagree=("warn_windows", "sum") if "warn_windows" in menu else ("rep", "size"),
        warn_f4=("warn_f4", "sum"), warning=("warning", "sum")).reset_index()
    figs = "\n\n".join(f"![k {t}](../figures/{Path(str(P.STEP3B3_FIG_K).format(truth=t)).name})" for t in S["truths"])
    report = f"""# WAA_estm step 3b-3: "draw the curve, then name it" on fake tails (known truth)

Written by `research/WAA_estm/step3b3_curve_then_equation.py` (libraries `waa_curve.py`, `waa_pooled.py`, `waa_sindy.py`,
`waa_sim.py`; settings `waa.curve3b3`). Fake data only: no real tail is fitted, the frozen tail set v1 is not changed.

**Overall verdict: {overall}.** Per truth: {', '.join(f"{r.truth} {r.verdict}" for r in dec.itertuples())}.

## 1. What was run

- **Data**: the fake datasets of steps 3b-1c / 3b-2 (event bootstrap of the {len(skel)}-tail TRAIN v1 skeleton, seeds
  {S['seed_base']} + replicate), {R} replicates. Main test: {', '.join(S['truths'])} with spread s = {S['s_true']} and 3a-fix
  noise. Extras: {', '.join(f"{KINDS[k]} ({', '.join(sorted({t for kk, t in datasets if kk == k}))})" for k in plan['kinds'] if k != 'main') or 'none'}.
  Windows: {', '.join(wins)}{'' if plan['extras_6h'] else ' (the extras: first 60 min only)'}.
- **Jobs**: {len(results)} free-curve pooled fits (dataset x window), each = {len(S['starts'])} starts + {len(S['lambda_grid'])}
  penalties + 2-fold cross-validation + curvature; {int(fits.nfev.sum())} likelihood evaluations in the main fits.
- **Cut for time**: {'; '.join(plan['cuts']) if plan['cuts'] else 'nothing'}. Timing: {plan['sec_per_point'] * 1e6:.2f}
  microseconds per tail-minute and evaluation on one replicate (m6), {plan['evals_per_job']} evaluations planned per job;
  projected total {plan['projected_s'] / 3600:.2f} h on {workers} processes (budget {S['time_budget_h']} h).
- **Run time**: {plan['wall_s'] / 60:.0f} min wall time ({fits.seconds.sum() / 3600:.1f} process-hours; cached in
  `.cache/fits3b3/`). This run: {time.time() - clock:.0f} s.

{md(fit_sum)}

## 2. Spline, penalty, supported range, starts

- **Law**: dW/dt = -f_i W k(W); log k = natural cubic spline in log W with {S['knots']} knots (parameters: log k at the
  knots), straight-line continuation outside the knots; f_i = exp(s g), g ~ N(0, 1) (median f = 1), marginalised as
  in 3b-1b; s fitted. Likelihood, offset prior and noise model: exactly those of the menu fits (`waa_pooled.py`;
  checked: a flat curve gives the same log-likelihood as F1 to 0.01).
- **Supported range of W**: from {S['w_floor_db']} dB (fixed: the smallest switch level of jump 2; below the 1 min noise
  sd of about 0.14 dB the curve is held by the offset prior and the early part of the tails, not seen directly) up to
  the level that has {S['support_min_minutes']} tail-minutes of A above it in the first 60 min: {fits.w_hi_db.min():.2f}-{fits.w_hi_db.max():.2f} dB
  over the datasets. Bands, tests and jump 2 use the range common to the replicates of a truth, on {S['grid_points']}
  log-spaced points.
- **Knots**: the two ends of that range and the {S['knots'] - 2} equal quantiles of log A in between (tail-minutes of the
  first 60 min); the same knots for both windows of a dataset. Example (m6, replicate 0): {fits[(fits.truth == 'm6') & (fits.kind == 'main')].knots_db.iloc[0]} dB.
- **Penalty**: lambda * integral of (d2 log k / d(log W)2)2, lambda from {S['lambda_grid']} by 2-fold cross-validation
  over EVENT DRAWS (fit on one half of the event draws, log-likelihood of the other half, both ways; highest sum).
  Plus a weak prior on log k at every knot (sd {S['prior_sd_logk']} around log(1/30 min)). Chosen penalties (jobs):

{md(lam_tab)}

- **Starts**: {', '.join(S['starts'])} (flat = exponential), {S['maxfev']['start']} evaluations each with lambda =
  {S['lambda_start']}; the best is continued ({S['maxfev']['finish']} more). `start_spread` = worst minus best penalised
  cost after the first stage; `starts_within_2` = starts within 2 of the best. Best start (jobs):

{md(start_tab)}

- **Uncertainty**: the band of a truth = min-max of its {R} replicates (each a new event bootstrap). Per replicate: a
  {z} sd interval from the curvature of the penalised fit (finite differences); `with_interval` = fits where that
  curvature was positive definite.

## 3. Jump 1: is the curve right? (`{rel(P.STEP3B3_J1)}`)

`band_coverage` = share of the supported range where the true k(W) lies inside the replicate band (rule: >=
{DEC['band_coverage']:.0%}); `median_k_ratio` = learned / true k (median over range and replicates); shape: slope of log k
vs log W over the range with its interval (flat = interval contains 0; rising = interval above 0), m6: mean k above
2 A_c / mean k below A_c > {DEC['step_ratio']} (rule: in >= {DEC['min_reps']} of 5).

{md(j1_show[j1_show.kind == 'main'].drop(columns='kind'))}

{figs}

**Post hoc diagnostic, added after the first complete run (it changes no verdict; `{rel(P.STEP3B3_POSTHOC)}`).**
Two things went wrong with the judging itself: the per-replicate interval exists for only
{int(fits.has_interval.sum())} of {len(fits)} fits (the curvature matrix was not positive definite), so every slope test
fails for lack of an interval; and the first knot (0.05 dB) is badly determined while the next knot is at
{fits.knots_db.str.split().str[1].astype(float).min():.2f}-{fits.knots_db.str.split().str[1].astype(float).max():.2f} dB. Slope of log k vs log W per replicate over the WHOLE range (min / median / max):

{md(fits.groupby(['kind', 'truth', 'window']).slope.agg(['min', 'median', 'max']).round(2).reset_index())}

The same fits judged only ABOVE THE SECOND KNOT (band coverage, level, slope over the replicates against the true
slope there, and the equation jump 2 gives on that range):

{md(post)}

## 4. Jump 2: is the equation right? (`{rel(P.STEP3B3_J2_SEL)}`, `{rel(P.STEP3B3_J2_COEF)}`)

h(W) = W k(W) of every replicate on the grid, weights 1 / variance over the replicates, library {', '.join(JTERMS)}
(sw x = max(W - x dB, 0)), thresholded backward elimination of `waa_sindy.py`; threshold by leave-one-replicate-out
cross-validation (one-SE rule). Replicates in which a term is selected (**bold** = a term of the true law):

{md(pivot_sel('main'))}

{md(j2[j2.kind == 'main'][[c for c in j2_cols if c in j2]])}

Coefficients, fitted / true (0 where not selected):

{md(coefs[coefs.kind == 'main'].drop(columns='kind')) if len(coefs) else ''}

## 5. m4 (two layers: not one curve): warning sign or not (`{rel(P.STEP3B3_MENU)}`)

Warning = the 60 min and 6 h curves of a replicate do not overlap (their {z} sd intervals) on more than
{W_['disagree_share']:.0%} of the range, OR BIC(free curve) - BIC(F4) > {W_['dbic_vs_f4']:g} in either window.

{md(m4.drop(columns=['truth']))}

The same signs for every truth (false alarms on the single-curve truths; F4 in the 6 h window was fitted for m4 only):

{md(warn_tab)}

## 6. Robustness: no spread, T1 noise

{md(j1_show[j1_show.kind != 'main']) if (j1_show.kind != 'main').any() else 'Not run (cut for time).'}

{md(j2[j2.kind != 'main'][[c for c in ['kind'] + j2_cols if c in j2]]) if (j2.kind != 'main').any() else ''}

Fitted spread s (median): {', '.join(f"{k} {t} {w_}: {v:.2f}" for (k, t, w_), v in fits.groupby(['kind', 'truth', 'window']).s.median().items())}.

## 7. Comparison with the menu fit (first 60 min; best of 4 starts, steps 3b-1c / 3b-2)

{md(menu[[c for c in ['truth', 'rep', 'menu_winner_60min', 'menu_dbic_to_second', 'curve_minus_best_menu_bic_60min'] if c in menu]].groupby('truth').agg(menu_winners=('menu_winner_60min', lambda v: ', '.join(f"{k} x{n}" for k, n in v.value_counts().items())), menu_dbic_to_second=('menu_dbic_to_second', 'median'), curve_minus_best_menu_bic=('curve_minus_best_menu_bic_60min', 'median')).reset_index())}

BIC of the free curve counts {S['knots'] + 1} parameters ({S['knots']} knots + s); positive = the menu law is preferred.

## 8. Verdict and what it means for Step 4

{md(dec)}

**Overall: {overall}.** PASS = all three criteria (curve, shape, equation) in the better window; FAIL = none; m4: PASS =
a warning sign in >= {DEC['min_reps']} of 5 replicates.

{MEANING}

## 9. Unclear / possibly wrong

{UNCLEAR}
"""
    save(report, P.STEP3B3_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {time.time() - clock:.0f} s")


# Written after the run of 2026-10-11 (the tables are computed on every run; these two texts are the reading of that
# run and quote its numbers: update them if a rerun with other settings changes the tables).
MEANING = """**In plain words: by the rule fixed before the run this is PARTIAL, and m1 fails. But the failure sits in one place:
below about 0.35 dB. Above it the free curve is right for every single-curve truth, also with the speed spread.**

- **By the rule (whole range from 0.05 dB).** m1 FAIL: the curve bends up at the lowest W (median slope -0.24 where
  0 is true), the true k is inside the band on 68 % / 48 % of the range, and the equation gets false terms (W^2, W^3,
  a switch at 0.05 dB) with the W coefficient x2.1. m3 and m5 PARTIAL: the equation is right in the first 60 min
  (W^2 in 5 of 5, coefficient x0.94 and x1.07; m5 also in 6 h, x0.97), the band covers the truth on 80 % / 93 % in
  the 6 h window, but the slope test cannot be passed. m6 PARTIAL: curve and step test pass in the first 60 min
  (k above 2 A_c / k below A_c = 7.7, in 4 of 5), but the equation is W^2 + W^3, never a switch term. m4 PASS: a
  warning sign in 4 of 5 (the free curve is worse than F4 by 20-36 BIC in the first 60 min).
- **Two things in the judging were broken, not the curve.** (1) The per-replicate interval from the curvature of
  the fit exists for 3 of 65 fits, so every slope test fails for lack of an interval, whatever the slope. (2) The
  first knot is at 0.05 dB and the next at 0.30-0.40 dB: a whole stretch with one free value that the data barely
  hold. For m1 with spread that value comes out 2-4 times too high in 9 of 10 fits (0.08-0.14 instead of 0.034
  per min); for m3 it is right (0.003) when the search started from "rising" and 7-15 times too high (0.02-0.04) when it
  started from "flat".
- **Post hoc, above the second knot only (0.35-1.7 dB; no verdict uses this).** First 60 min, with spread:
  m1 slope -0.02 (true 0), level x1.03, equation "W" in 5 of 5 with coefficient x1.03: **the mixing trap is gone
  there** (3b-2 gave x1.79 and a false A^2). m3 slope 1.08 (true 1.0), equation W^2 x0.94 in 5 of 5. m5 slope 0.96
  (true 0.92), W^2 x1.06 in 5 of 5. m6 slope 0.27 (true 0.26), level x1.02. The four laws give four different,
  correct slopes (0, 0.26, about 1, 1), and m4 gives 0.37. The first 60 min is the cleaner window; in 6 h the
  equation of m1 picks up false terms.
- **What the curve cannot show.** m6's step is at 0.16 dB, inside the stretch that is not determined: above it the
  switch law is a smooth, slowly rising curve, and jump 2 names it "W + W^2". So the free curve cannot tell a
  switch from a mild nonlinearity; only the menu fit (F6) can, through the level it uses below 0.3 dB. m4 looks
  like a clean rising curve too (equation W^2 + W^3 in 5 of 5); only the comparison with F4 gives it away.
- **Free curve against the menu fit (first 60 min).** The menu fit names the right family in 23 of 25 datasets. The
  free curve has a LOWER BIC than the true menu law for m1, m3, m5 (median -21, -32, -62): a flexible curve should
  not beat the true law by that much. The gain comes from the low end (my reading: the likelihood has no 0.3 dB
  quantization, so below about 0.15 dB a signal that sticks to one level looks like faster drying). Not proven.

**For Step 4 (real tails).**
1. Use the free-curve fit to DRAW the drying curve, in the first 60 min, and read it only above about 0.35 dB (more
   exactly: above the second knot). There it gave the right level (within 2-6 %) and the right slope on fake data
   with the real spread and noise. Do not interpret the curve below that.
2. Jump 2 on that upper range names exponential and quadratic laws correctly. Read "W + W^2" as "mildly rising",
   which is also what a switch law looks like from above.
3. Keep the menu fits next to it: the switch (F6) and the two layers (F4) are visible only there. If the free
   curve is worse than F4 by more than 10 BIC, do not trust the single curve (m4 check: 4 of 5; no false alarm
   on the other truths in 20 datasets).
4. Before Step 4, three repairs are needed and should be tested on these same cached fake datasets: an interval
   that works (bootstrap over events, or a profile, instead of the curvature); knots that do not leave
   0.05-0.35 dB to one free value (or a range that starts at the second knot); and a longer or better search
   (the starts still differ by 30-40 in cost when the first stage ends)."""
UNCLEAR = """1. **The low-W problem may be a property of the likelihood, not only of the spline.** The free curve beats the
   true menu law by 21-62 BIC, and in my first check a curve with k held constant below 0.05 dB already gained 7
   log-likelihood units over the exact F3 law on m3 data. I suspect the missing quantization model (noted as open
   in 3b-1b), but I did not test it. If true, it also touches the menu fits where they use the low end: F6's
   A_c = 0.16 dB and F4's slow layer.
2. **The supported range was my definition** (floor 0.05 dB because jump 2's smallest switch is there; top = 200
   tail-minutes above). A count of tail-minutes cannot find the floor, because late minutes with A near 0 are
   plentiful but carry noise only. The floor should come from the noise level or from the knots.
3. **Knots at quantiles of log A put none between 0.05 and 0.30-0.40 dB.** With 5 knots the step of m6 at 0.16 dB
   could not have been drawn even with perfect data.
4. **The search is not finished when it stops.** 150 evaluations per start, then 450 more for the best one
   (Nelder-Mead, 6 parameters). Every final fit reports that it met the stopping tolerance, but the starts
   differ by a median of 34-41 in cost after the first stage. `starts_within_2` is 1 everywhere, so the "spread
   over starts" says the starts had not converged, not that there are several optima. The best start was
   "falling" for m1 in 9 of 10 main fits.
5. **The curvature interval failed** (positive definite in 3 of 65). So "flat" and "rising" could not be tested as
   written, and the "60 min and 6 h curves disagree beyond their bands" warning could never fire (0 of 25). The
   m4 PASS rests on the BIC comparison with F4 alone.
6. **The step test for m6 is not specific**: k above 2 A_c over k below A_c > 2 is also true for any rising law
   (m3 would pass it). Its pass says "k rises", not "there is a step".
7. **The band is the min-max of only 5 replicates** and is narrow where the fits agree (above 0.35 dB), so a
   level error of 5 % already puts the truth outside it (m3, first 60 min: coverage 68 % with the level at x0.95).
   Coverage is a harsh number there and a generous one at the low end, where the band is a decade wide.
8. **The post hoc diagnostic was added after seeing the results.** It reuses the cached fits (no refit) and the
   decision table is byte-identical before and after. Its range starts at the largest second knot of the replicates.
9. **One cut that was not on the list**: before dropping T1 altogether I dropped only the 6 h window of the extras
   (they keep the first 60 min). Projection 3.61 h -> 2.55 h; the run took 2.1 h. The projection came from timing
   likelihood evaluations on one replicate (2,053 evaluations planned per job), not from running a whole replicate.
10. **Penalty by 2-fold cross-validation** picked the smallest and the largest value about equally often (main,
    60 min: 11 / 4 / 10), which suggests the held-out likelihood hardly depends on it with 5 knots.
11. **Jump 2 weights** (1 / variance over 5 replicates) are very uneven; in the 6 h window they made the m1
    equation unstable (coefficient x-4 on the upper range) although the curve itself is flat there.
12. **T1 noise** (first 60 min only): above the second knot m1 stays flat (slope -0.02, level x1.03) and m3 keeps
    W^2 x0.96, but the m1 equation picks up false terms in 4 of 5.
13. **Key figure**: the limit of 6 is reached, so the five `step3b3_k_<truth>.png` are diag. If one should become
    key I would take `step3b3_k_m1.png` only after the repairs above (it now shows the low-W artefact), and demote
    `outputs/tables/step3b1b_recovery.csv`.
"""


if __name__ == "__main__":
    main()
