"""WAA_estm step 3b-2: weak-form SINDy on fake tails (known truth). Fake data only; no real tail is fitted.

Question: what does a SINDy-style equation finder return when we know the true drying law?
  Q1  the true terms for m1 (exponential), m3 (-c W^2), m5 (-a W - b W^2), m6 (state switch)?
  Q2  what for m4 (two layers: a hidden state that a library in A alone cannot hold)?
  Q3  mixture trap: does pooling exponential tails of different speed (spread s = 1.25) create a false nonlinear term?
  Q4  does it agree with the pooled fit of step 3b-1b / 3b-1c on the same datasets?

Method (library waa_sindy.py; pysindy is not installed, so the weak form is written out there):
  - data: the fake datasets of step 3b-1c, test T2 (3a-fix noise, event-bootstrap replicates, seeds 3000 + replicate),
    truths m1, m3, m4, m5, m6 with spread s = 1.25, plus m1, m3, m5, m6 without spread (s = 0); two windows: the first
    60 min and the whole tail;
  - one law for all tails of a dataset: dA/dt = sum_j xi_j theta_j(A), theta = 1, A, A^2, A^3, max(A - A_c, 0);
    weak-form rows per tail, weighted by 1 / (1 min noise sd of the sub-link);
  - sparse regression with a threshold chosen by cross-validation over event draws; stability selection over
    event-bootstrap resamples (how often is each term selected?);
  - the unknown offset c of each tail (A = W + c), three ways: `naive` (ignore it), `handled` (a free constant per
    tail in dA/dt: exact for a linear law, because dA/dt = -k (A - c) = -k A + k c) and `oracle` (the true c removed;
    only possible on fake data, it shows the ceiling);
  - rescaling check: the time axis of each tail divided by its speed factor, estimated by the pooled fit (`fit`) or the
    true one (`true`).
The decision rule is fixed in `waa.sindy3b2.decision` before the run. The only uses of real data: the NaN minutes of
the real tails (copied into the fake tails, as in 3a). The frozen tail set v1 is read, never changed.

One dataset per job, jobs in parallel processes (cores - 1); every result is cached in .cache/sindy3b2/ (pooled
fits per candidate and start, SINDy results per dataset), so a restart resumes. Delete that folder to redo.
Run from the repository root (about 20 minutes the first time, 1 minute when cached):
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step3b2_sindy_fake.py
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
import waa_sindy as sd
import step3b1c_stress_test as st                          # the same event bootstrap, noise tables and T2 starts

CFG = yaml.safe_load(open(P.CONFIG_FILE, encoding="utf-8"))
SIM, S, S1C = CFG["waa"]["sim"], CFG["waa"]["sindy3b2"], CFG["waa"]["stress3b1c"]
STEP, SPM = SIM["rsl_step_db"], SIM["samples_per_min"]
GRID = S["switch_grid_db"]
TERMS = sd.term_names(GRID)
SW = [t for t in TERMS if t.startswith("sw")]
RULE = S["decision"]["true_terms"]
THR = S["stability"]["select_freq"]
FAMILY = {"F1": "exponential", "F3": "nonlinear", "F5": "nonlinear", "F6": "switch", "F4": "two layers"}
INK, INK2 = "#1a1a1a", "#52514e"
CACHE = P.CACHE / "sindy3b2"
STRESS_CACHE = P.CACHE / "fits3b1c"                        # T2 fits of step 3b-1c (reused when the data are identical)


def rel(p):
    return Path(p).resolve().relative_to(P.REPO).as_posix()


def md(df, index=False):
    df = df.reset_index() if index else df
    head = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns)
    return head + "\n" + "\n".join("| " + " | ".join("" if pd.isna(x) else str(x) for x in row) + " |"
                                   for row in df.itertuples(index=False))


def tag(job):
    return f"{job['truth']}_{'spread' if job['spread'] else 'nospread'}_{job['rep']}"


# =================================================================================================== one job (worker)
_INPUTS = {}


def pooled_part(job, tails, noise, offset_sd):
    """Pooled fit of every candidate from the T2 starts of step 3b-1c (best likelihood per candidate), the BIC
    winner, and the winner's speed factor per tail. Fits are cached; identical 3b-1c fits are reused."""
    D = wp.pack(tails, noise, offset_sd, S1C["window_min"], S["pooled_min_points"])
    used = set(D["tail_id"])
    n_ind = len({(tl["cml_id"], tl["event_copy"]) for tl in tails if tl["tail_id"] in used})
    fits, f1_best = [], None
    for cand in S1C["candidates"]:
        starts = {k: (v[0], v[1]) for k, v in S1C["t2_starts"][cand].items()}
        if cand != "F1":
            starts["f1sol"] = (None, None)
        for label, (theta0, s0) in starts.items():
            f_ = CACHE / f"pfit_{job['truth']}_{job['rep']}_{cand}_{label}.json"
            old = STRESS_CACHE / f"fit_T2_{job['truth']}_{job['rep']}_{cand}_{label}.json"
            if f_.exists():
                r = json.loads(f_.read_text())
            elif job["reuse"] and old.exists():
                r = {**json.loads(old.read_text()), "reused": True}
                f_.write_text(json.dumps(r))
            else:
                if label == "f1sol":
                    theta0, s0 = st.f1_start(cand, f1_best["theta"][0]), f1_best["s"]
                r = wp.fit(cand, D, theta0, s0, S1C["maxfev_generic"])
                r.update({"candidate": cand, "start": label, "n_independent": n_ind, "k": len(wp.NAMES[cand]) + 1,
                          "reused": False})
                f_.write_text(json.dumps(r))
            fits.append(r)
        if cand == "F1":
            f1_best = max((r for r in fits if r["candidate"] == "F1"), key=lambda r: r["ll"])
    best = {}
    for r in fits:                                         # best start per candidate, then the lowest BIC
        r["bic"] = -2 * r["ll"] + r["k"] * np.log(n_ind)
        if r["candidate"] not in best or r["ll"] > best[r["candidate"]]["ll"]:
            best[r["candidate"]] = r
    win = min(best.values(), key=lambda r: r["bic"])
    speeds = wp.tail_speed_factors(win["candidate"], win["theta"], win["s"], D)
    pooled = {"winner": win["candidate"], "winner_theta": win["theta"], "winner_s": win["s"],
              "bic": {c: round(r["bic"], 2) for c, r in best.items()}, "f1_tau": best["F1"]["theta"][0],
              "n_new_fits": sum(not r.get("reused") for r in fits)}
    return pooled, speeds


def design(tails, window, offset, sigma):
    """Stack the weak-form rows of all tails. Returns (b, G, event draw per row, tail position per row)."""
    B, G, grp, pos = [], [], [], []
    for i, tl in enumerate(tails):
        a = tl["a"][:window] if window else tl["a"]
        if offset == "oracle":
            a = a - tl["offset"]
        b, g = sd.weak_rows(a, S["half_width_min"], S["stride_min"], GRID, S["test_power"], S["max_gap_min"])
        if len(b) < (2 if offset == "handled" else 1):
            continue
        if offset == "handled":                            # a free constant per tail = remove the tail's mean of
            b, g = b - b.mean(), g - g.mean(axis=0)        # every column (the constant column is the same in every row)
        w = 1.0 / sigma[tl["sublink_id"]]
        B.append(b * w), G.append(g * w)
        grp += [tl["event_copy"]] * len(b)
        pos += [i] * len(b)
    return np.concatenate(B), np.vstack(G), np.asarray(grp), np.asarray(pos)


def sindy_run(b, G, groups, seed):
    rng = np.random.default_rng(seed)
    lam, cv = sd.cv_lambda(G, b, groups, S["lambda_grid"], S["cv_folds"], rng)
    freq = sd.stability(G, b, groups, lam, S["stability"]["resamples"], rng)
    xi = sd.refit(G, b, freq >= THR)                       # coefficients: OLS on the stable terms
    # diagnostic added after the first run (it does not enter the rule): the threshold with the LOWEST held-out error
    i_min = int(np.nanargmin(cv["cv_mean"]))
    lam_min = cv["lambda"][i_min]
    fmin = freq if lam_min == lam else sd.stability(G, b, groups, lam_min, S["stability"]["resamples"], rng)
    xmin = sd.refit(G, b, fmin >= THR)
    return {"lam": lam, "n_rows": int(len(b)), **{f"freq_{t}": float(v) for t, v in zip(TERMS, freq)},
            **{f"coef_{t}": float(v) for t, v in zip(TERMS, xi)},
            "lam_min": lam_min, "cv_err_min": cv["cv_mean"][i_min], "cv_err_chosen": cv["cv_mean"][cv["lambda"].index(lam)],
            "cv_err_largest": cv["cv_mean"][-1], "cv_se_min": cv["cv_se"][i_min],
            **{f"fmin_{t}": float(v) for t, v in zip(TERMS, fmin)}, **{f"cmin_{t}": float(v) for t, v in zip(TERMS, xmin)}}


def run_dataset(job):
    """Everything for ONE fake dataset: pooled fit (spread datasets), then SINDy for every window, offset mode and
    rescaling. The result is cached as one json file."""
    f_ = CACHE / f"res_{tag(job)}.json"
    if f_.exists():
        return json.loads(f_.read_text())
    if not _INPUTS:
        _, _INPUTS["noise"], _INPUTS["offset_sd"] = st.load_noise()
    noise, offset_sd = _INPUTS["noise"], _INPUTS["offset_sd"]
    sigma = {k: float(np.sqrt(v[1] ** 2 + (v[0] ** 2 + STEP ** 2 / 12) / SPM)) for k, v in noise.items()}   # 1 min noise sd
    rows = pd.read_parquet(job["sim"]).sort_values(["tail_id", "minutes_since_t0_trace"])
    tails = [{"tail_id": t, "sublink_id": g.sublink_id.iloc[0], "cml_id": g.cml_id.iloc[0],
              "event_copy": int(g.event_copy.iloc[0]), "A_start": g.A_start_trace_db.iloc[0],
              "offset": float(g.offset_db.iloc[0]), "f_true": float(g.f_time_scale.iloc[0]),
              "t": g.minutes_since_t0_trace.to_numpy(), "a": g.A_db.to_numpy(float)} for t, g in rows.groupby("tail_id")]
    out = {"truth": job["truth"], "spread": job["spread"], "rep": job["rep"], "n_tails": len(tails), "pooled": None}
    t0 = time.time()
    scales = {"none": None}
    if job["spread"]:
        out["pooled"], speeds = pooled_part(job, tails, noise, offset_sd)
        scales["fit"] = np.array([speeds.get(tl["tail_id"], 1.0) for tl in tails])
        scales["true"] = np.array([tl["f_true"] for tl in tails])
    out["pooled_seconds"] = round(time.time() - t0, 1)
    t0 = time.time()
    runs = []
    for k, (wname, window) in enumerate(job["windows"].items()):
        for offset in S["offsets"]:
            b, G, grp, pos = design(tails, window, offset, sigma)
            for rname, f in scales.items():
                if rname != "none" and offset != "handled":   # the rescaling check is run with the handled offset only
                    continue
                Gs = G if f is None else G / f[pos][:, None]   # f dA/dt = h(A): divide the right-hand side by f
                seed = [S["stability"]["seed"], job["rep"], k, S["offsets"].index(offset), list(scales).index(rname)]
                runs.append({"window": wname, "offset": offset, "rescale": rname, "n_tails_used": int(len(set(pos))),
                             **sindy_run(b, Gs, grp, seed)})
    out["sindy"], out["sindy_seconds"] = runs, round(time.time() - t0, 1)
    tmp = f_.with_suffix(".tmp")
    tmp.write_text(json.dumps(out))
    tmp.replace(f_)
    return out


# =================================================================================================== analysis helpers
def family(selected):
    """Name the kind of law a set of selected terms stands for."""
    s = set(selected)
    poly, sw = bool(s & {"A^2", "A^3"}), bool(s & set(SW))
    if not s:
        return "nothing"
    if poly and sw:
        return "mixed"
    if poly:
        return "nonlinear"
    if sw:
        return "switch"
    return "exponential" if "A" in s else "constant only"


def role(truth, term):
    if truth not in RULE:
        return "-"
    r = RULE[truth]
    return "true" if any(term in g for g in r["required"]) else ("allowed" if term in r["allowed"] else "spurious")


def eval_rule(g, truth, true_coef):
    """The pre-set rule for one truth in one configuration. g: one row per replicate (sel_<term>, coef_<term>)."""
    r, D, n = RULE[truth], S["decision"], len(g)
    need = [int(g[[f"sel_{t}" for t in grp]].any(axis=1).sum()) for grp in r["required"]]
    spur = {t: int(g[f"sel_{t}"].sum()) for t in TERMS if t not in r["allowed"] and g[f"sel_{t}"].sum() > 0}
    ratios = {}
    for c in r["coefs"]:
        est = g[["coef_sw0.1", "coef_sw0.2"]].sum(axis=1) if c == "switch" else g[f"coef_{c}"]
        ratios[c] = float((est / true_coef[c]).median())
    true_ok = all(k / n >= D["min_true"] for k in need)
    spur_ok = all(v <= D["max_spurious"] for v in spur.values())
    coef_ok = all(abs(v - 1) <= D["coef_tol"] for v in ratios.values())
    why = []
    if not true_ok:
        why.append("true " + ", ".join(f"{'/'.join(grp)} {k}/{n}" for grp, k in zip(r["required"], need)
                                       if k / n < D["min_true"]))
    if not spur_ok:
        why.append("spurious " + ", ".join(f"{t} {v}/{n}" for t, v in spur.items() if v > D["max_spurious"]))
    if not coef_ok:
        why.append("coef " + ", ".join(f"{c} x{v:.2f}" for c, v in ratios.items() if abs(v - 1) > D["coef_tol"]))
    return {"true_selected": "; ".join(f"{'/'.join(grp)} {k}/{n}" for grp, k in zip(r["required"], need)),
            "spurious": ", ".join(f"{t} {v}/{n}" for t, v in spur.items()),
            "coef_ratios": ", ".join(f"{c} x{v:.2f}" for c, v in ratios.items()),
            "true_ok": true_ok, "spurious_ok": spur_ok, "coef_ok": coef_ok, "met": true_ok and spur_ok and coef_ok,
            "why_not": "; ".join(why)}


def selection_figure(sel, panels, R):
    """One panel per configuration: truth (rows) x term (columns), text = replicates in which the term is selected.
    Dashed box = a true term; hatched = a term outside the truth selected in more than one replicate."""
    fig, ax = plt.subplots(1, len(panels), figsize=(6.3 * len(panels), 4.6))
    for a_, (title, q) in zip(np.atleast_1d(ax), panels.items()):
        g = sel.query(q)
        truths = [t for t in S["truths"] if t in set(g.truth)]
        m = g.pivot(index="truth", columns="term", values="n_selected").reindex(index=truths, columns=TERMS)
        a_.imshow(m.to_numpy(float), cmap="Greys", vmin=0, vmax=R * 1.8, aspect="auto")
        for i, tr in enumerate(truths):
            for k, term in enumerate(TERMS):
                n, ro = int(m.iloc[i, k]), role(tr, term)
                a_.text(k, i, f"{n}/{R}", ha="center", va="center", fontsize=10,
                        fontweight="bold" if ro == "true" else "normal")
                if ro == "true":
                    a_.add_patch(plt.Rectangle((k - 0.46, i - 0.46), 0.92, 0.92, fill=False, ec=INK, lw=2, ls="--"))
                elif ro == "spurious" and n > S["decision"]["max_spurious"]:
                    a_.add_patch(plt.Rectangle((k - 0.5, i - 0.5), 1, 1, fill=False, ec=INK2, lw=0, hatch="///"))
        a_.set_xticks(range(len(TERMS)), [t.replace("sw", "sw ") for t in TERMS], fontsize=8)
        a_.set_yticks(range(len(truths)), [f"truth {t}" for t in truths])
        a_.set_xlabel("library term (sw x = max(A - x dB, 0))")
        a_.set_title(title, fontsize=10)
    fig.suptitle(f"Weak-form SINDy on fake tails, first 60 min: replicates (of {R}) in which each term is selected.\n"
                 "Dashed box, bold = a term of the true law; hatched = a term outside the truth selected in more than "
                 "one replicate (m4: no single-variable truth)", fontsize=10)
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
    _, noise, _ = st.load_noise()
    cal = pd.read_csv(P.STEP3A_CALIBRATION).set_index("law")
    par = {law: {k: float(cal.loc[law, k]) for k in ws.PARAMS[law]} for law in ("m1", "m3", "m4", "m5")}
    m6 = pd.read_csv(P.STEP3B1_M6)
    b6 = m6.loc[m6.rmse_band_observed.idxmin()]
    par["m6"] = {k: float(b6[k]) for k in ws.PARAMS["m6"]}
    true_coef = {"m1": {"A": -1 / par["m1"]["tau"]}, "m3": {"A^2": -par["m3"]["c"]},
                 "m5": {"A": -par["m5"]["a"], "A^2": -par["m5"]["b"]},
                 "m6": {"A": -par["m6"]["k_evap"], "switch": -par["m6"]["k_run"]}}
    print(f"{len(skel)} skeleton tails; truths {', '.join(S['truths'])} (+ no spread: {', '.join(S['nospread_truths'])}); "
          f"{workers} worker processes; library: {', '.join(TERMS)}")

    # ------------------------------------------------------------------------------------------- simulate (cached)
    sim_file = lambda t, spread, rep: CACHE / f"sim_{t}_{'spread' if spread else 'nospread'}_{rep}.parquet"   # noqa: E731
    reuse = {}
    for rep in range(R0):
        seed = S["seed_base"] + rep
        sk = st.bootstrap_skeleton(skel, np.random.default_rng(seed))
        masks = {r.tail_id: nan_real[r.orig_tail_id] for r in sk.itertuples()}
        z = ws.spread_z(sk, seed)
        for spread in (True, False):
            for t in (S["truths"] if spread else S["nospread_truths"]):
                f_ = sim_file(t, spread, rep)
                if not f_.exists():
                    f = np.exp(S["s_true"] * z) if spread else np.ones(len(sk))
                    rows, info = ws.simulate(sk, t, par[t], noise, SIM["switches"], seed, f, STEP, SPM, masks)
                    rows = rows.merge(sk[["tail_id", "event_copy", "A_start_trace_db"]], on="tail_id").merge(
                        info[["tail_id", "offset_db", "f_time_scale"]], on="tail_id")
                    rows.to_parquet(f_)
                old = STRESS_CACHE / f"sim_T2_{t}_{rep}.parquet"     # the 3b-1c dataset: reuse its fits only if identical
                if spread and old.exists():
                    a = pd.read_parquet(f_, columns=["tail_id", "minutes_since_t0_trace", "A_db"])
                    o = pd.read_parquet(old, columns=["tail_id", "minutes_since_t0_trace", "A_db"])
                    reuse[(t, rep)] = len(a) == len(o) and bool(np.array_equal(a.to_numpy(), o.to_numpy(), equal_nan=True))
    print(f"datasets simulated ({time.time() - clock:.0f} s); identical to the step 3b-1c T2 datasets: "
          f"{sum(reuse.values())} of {len(reuse)} compared (their pooled fits are reused)")

    def jobs_for(replicates, windows, nospread_truths):
        return [{"truth": t, "spread": spread, "rep": rep, "sim": str(sim_file(t, spread, rep)), "windows": windows,
                 "reuse": reuse.get((t, rep), False)}
                for spread in (True, False) for t in (S["truths"] if spread else nospread_truths)
                for rep in range(replicates)]

    # ------------------------------------------------------------------------------------------- time one replicate
    plan_file = CACHE / "plan.json"
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
        print(f"plan from {rel(plan_file)}: cut: {plan['cuts'] or 'nothing'}")
    else:
        print("\ntiming one replicate (m3 with spread, replicate 0: pooled fits from all starts + every SINDy run) ...",
              flush=True)
        tiny = [{"tail_id": i, "sublink_id": next(iter(noise)), "cml_id": i, "A_start": 1.0, "t": np.arange(20.0),
                 "a": np.exp(-np.arange(20.0) / 30)} for i in range(2)]
        for c in S1C["candidates"]:                                 # compile the numba code before the clock starts
            wp.pooled_loglik(c, [1.0] * len(wp.NAMES[c]) if c != "F4" else [10.0, 100.0, 0.3], 1.0,
                             wp.pack(tiny, noise, {}, None, 10))
        timed = run_dataset(next(j for j in jobs_for(1, S["windows"], S["nospread_truths"])
                                 if j["truth"] == "m3" and j["spread"]))
        t_pool, t_sindy = timed["pooled_seconds"], timed["sindy_seconds"]
        n_win = len(S["windows"])

        def projected(replicates, windows, nospread_truths):
            jobs = jobs_for(replicates, windows, nospread_truths)
            n_fit = sum(j["spread"] and not j["reuse"] for j in jobs)
            n_sindy = sum(5 if j["spread"] else 3 for j in jobs) / 5      # in units of one spread dataset
            return (n_fit * t_pool + n_sindy * t_sindy * len(windows) / n_win) / workers

        plan = {"replicates": R0, "windows": S["windows"], "nospread_truths": S["nospread_truths"], "cuts": [],
                "wall_s": 0.0, "timed_pooled_s": t_pool, "timed_sindy_s": t_sindy}
        proj = projected(R0, plan["windows"], plan["nospread_truths"])
        print(f"   pooled fits {t_pool:.0f} s, SINDy runs {t_sindy:.0f} s; projected total {proj / 60:.0f} min with "
              f"{workers} processes (budget {S['time_budget_h']} h)", flush=True)
        for cut, change in (("replicates 5 -> 3", {"replicates": 3}),
                            ("whole-tail window dropped", {"windows": {"first_60min": S["windows"]["first_60min"]}}),
                            ("no-spread variant of m5 dropped",
                             {"nospread_truths": [t for t in S["nospread_truths"] if t != "m5"]})):
            if proj <= S["time_budget_h"] * 3600:
                break
            plan.update(change)
            plan["cuts"].append(cut)
            proj = projected(plan["replicates"], plan["windows"], plan["nospread_truths"])
            print(f"   cut: {cut} -> projected {proj / 60:.0f} min", flush=True)
        plan["projected_s"] = round(proj)
        plan_file.write_text(json.dumps(plan))
    R = plan["replicates"]
    jobs = jobs_for(R, plan["windows"], plan["nospread_truths"])

    # ------------------------------------------------------------------------------------------- run (parallel)
    cached0 = len(list(CACHE.glob("res_*.json")))
    print(f"\n{len(jobs)} datasets ({cached0} already cached) on {workers} processes", flush=True)
    t_run, results = time.time(), []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(run_dataset, j): j for j in jobs}
        for k, fu in enumerate(as_completed(futs)):
            results.append(fu.result())
            el = time.time() - t_run
            print(f"   datasets {k + 1}/{len(jobs)}: last {tag(futs[fu])} ({el:.0f} s elapsed, about "
                  f"{el / (k + 1) * (len(jobs) - k - 1):.0f} s left)", flush=True)
    if len(list(CACHE.glob("res_*.json"))) > cached0:
        plan["wall_s"] = round(plan["wall_s"] + time.time() - clock, 1)
        plan_file.write_text(json.dumps(plan))

    # ------------------------------------------------------------------------------------------- tables
    fits = pd.DataFrame([{"truth": r["truth"], "spread": r["spread"], "rep": r["rep"], **run}
                         for r in results for run in r["sindy"]])
    for t in TERMS:
        fits[f"sel_{t}"] = fits[f"freq_{t}"] >= THR
    fits["selected"] = fits.apply(lambda r: " + ".join(t for t in TERMS if r[f"sel_{t}"]) or "(nothing)", axis=1)
    fits["family"] = fits.apply(lambda r: family([t for t in TERMS if r[f"sel_{t}"]]), axis=1)
    cfg = ["spread", "window", "offset", "rescale"]
    fits = fits.sort_values(["truth"] + cfg + ["rep"]).reset_index(drop=True)
    save(fits.round(6), P.STEP3B2_FITS, index=False)

    sel = []
    for keys, g in fits.groupby(["truth"] + cfg):
        for t in TERMS:
            sel.append({**dict(zip(["truth"] + cfg, keys)), "term": t, "role": role(keys[0], t),
                        "n_selected": int(g[f"sel_{t}"].sum()), "replicates": len(g),
                        "mean_frequency": round(float(g[f"freq_{t}"].mean()), 3)})
    sel = pd.DataFrame(sel)
    save(sel, P.STEP3B2_SELECTION, index=False)

    coefs, rule = [], []
    for keys, g in fits[fits.truth.isin(RULE)].groupby(["truth"] + cfg):
        truth = keys[0]
        for c, tc in true_coef[truth].items():
            est = g[["coef_sw0.1", "coef_sw0.2"]].sum(axis=1) if c == "switch" else g[f"coef_{c}"]
            ratio = est / tc
            coefs.append({**dict(zip(["truth"] + cfg, keys)), "coefficient": c, "true": round(tc, 5),
                          "median_ratio": round(float(ratio.median()), 3), "min_ratio": round(float(ratio.min()), 3),
                          "max_ratio": round(float(ratio.max()), 3), "replicates_nonzero": int((est != 0).sum()),
                          "counted": c in RULE[truth]["coefs"]})
        rule.append({**dict(zip(["truth"] + cfg, keys)), **eval_rule(g, truth, true_coef[truth])})
    coefs, rule = pd.DataFrame(coefs), pd.DataFrame(rule)
    save(coefs, P.STEP3B2_COEFS, index=False)

    M = S["main"]
    is_main = (fits.window == M["window"]) & (fits.offset == M["offset"]) & (fits.rescale == M["rescale"])
    main_q = f"window == '{M['window']}' and offset == '{M['offset']}' and rescale == '{M['rescale']}'"
    rule["main"] = (rule.window == M["window"]) & (rule.offset == M["offset"]) & (rule.rescale == M["rescale"]) & rule.spread
    q1 = rule[rule.main]
    v1 = "PASS" if q1.met.all() else ("FAIL" if not q1.met.any() else "PARTIAL")

    # Q2 (m4) and Q3 (m1 with spread): what comes out in the main configuration
    D = S["decision"]
    m4 = fits[is_main & fits.spread & (fits.truth == "m4")]
    clean = m4.family.isin(["exponential", "nonlinear", "switch"])
    v2 = "FAIL" if clean.mean() >= D["q2_q3"]["clean_share"] else ("PASS" if clean.sum() <= 1 else "PARTIAL")
    m1s = fits[is_main & fits.spread & (fits.truth == "m1")]
    false_terms = {t: int(m1s[f"sel_{t}"].sum()) for t in TERMS if t not in ("1", "A") and m1s[f"sel_{t}"].sum() > 0}
    worst = max(false_terms.values(), default=0)
    v3 = "FAIL" if worst / len(m1s) >= D["q2_q3"]["clean_share"] else ("PASS" if worst <= 1 else "PARTIAL")

    # Q4: the family named by SINDy vs the family of the pooled-fit winner, per dataset with spread
    agree = []
    for r in results:
        if not r["spread"]:
            continue
        g = fits[(fits.truth == r["truth"]) & fits.spread & (fits.rep == r["rep"]) & (fits.window == M["window"])
                 & (fits.offset == M["offset"])].set_index("rescale")
        pw = r["pooled"]["winner"]
        bic = sorted(r["pooled"]["bic"].values())
        agree.append({"truth": r["truth"], "rep": r["rep"], "n_tails": r["n_tails"], "pooled_winner": pw,
                      "pooled_family": FAMILY[pw], "pooled_dbic_to_second": round(bic[1] - bic[0], 1),
                      "sindy_terms": g.loc["none", "selected"], "sindy_family": g.loc["none", "family"],
                      "sindy_rescaled_terms": g.loc["fit", "selected"], "sindy_rescaled_family": g.loc["fit", "family"],
                      "pooled_F1_rate": round(1 / r["pooled"]["f1_tau"], 4),
                      "sindy_A_coef": round(-g.loc["none", "coef_A"], 4),
                      "sindy_rescaled_A_coef": round(-g.loc["fit", "coef_A"], 4)})
    agree = pd.DataFrame(agree).sort_values(["truth", "rep"]).reset_index(drop=True)
    agree["agree"] = agree.pooled_family == agree.sindy_family
    agree["agree_rescaled"] = agree.pooled_family == agree.sindy_rescaled_family
    cmp_ = agree[agree.truth.isin(RULE)]                             # m4 has no single-variable family: not counted
    share4 = float(cmp_.agree.mean())
    v4 = "PASS" if share4 >= D["q4"]["pass"] else ("PARTIAL" if share4 >= D["q4"]["partial"] else "FAIL")
    save(agree, P.STEP3B2_AGREEMENT, index=False)

    verdicts = pd.DataFrame([
        {"question": "Q1", "verdict": v1, "basis": f"{int(q1.met.sum())} of {len(q1)} truths meet the rule: "
         + "; ".join(f"{r.truth} {'met' if r.met else r.why_not}" for r in q1.itertuples())},
        {"question": "Q2", "verdict": v2, "basis": "m4 returns: " + ", ".join(
            f"{k} x{v}" for k, v in m4.selected.value_counts().items())},
        {"question": "Q3", "verdict": v3, "basis": "m1 with spread returns: " + ", ".join(
            f"{k} x{v}" for k, v in m1s.selected.value_counts().items())},
        {"question": "Q4", "verdict": v4, "basis": f"same family in {int(cmp_.agree.sum())} of {len(cmp_)} datasets "
         f"(m1, m3, m5, m6); with time rescaling {int(cmp_.agree_rescaled.sum())} of {len(cmp_)}"}])
    dec = pd.concat([rule.assign(kind="rule"), verdicts.assign(kind="verdict")], ignore_index=True)
    save(dec, P.STEP3B2_DECISION, index=False)

    # post hoc diagnostic (not part of the rule): the same runs with the threshold at the lowest held-out error
    alt = fits.copy()
    for t in TERMS:
        alt[f"sel_{t}"], alt[f"coef_{t}"] = alt[f"fmin_{t}"] >= THR, alt[f"cmin_{t}"]
    alt["selected"] = alt.apply(lambda r: " + ".join(t for t in TERMS if r[f"sel_{t}"]) or "(nothing)", axis=1)
    alt = alt[(alt.window == M["window"]) & (alt.offset == M["offset"])]
    alt_rule = pd.DataFrame([{"spread": k[1], "rescale": k[2], "truth": k[0], **eval_rule(g, k[0], true_coef[k[0]])}
                             for k, g in alt[alt.truth.isin(RULE)].groupby(["truth", "spread", "rescale"])])
    cv_tab = fits[is_main].groupby(["spread", "truth"]).agg(
        threshold_rule=("lam", "median"), threshold_lowest_error=("lam_min", "median"),
        heldout_error_rule=("cv_err_chosen", "median"), heldout_error_lowest=("cv_err_min", "median"),
        heldout_error_one_term=("cv_err_largest", "median"), se_of_lowest=("cv_se_min", "median")).round(4).reset_index()

    # ------------------------------------------------------------------------------------------- figure
    w60 = f"window == '{M['window']}' and offset == '{M['offset']}'"
    fig = selection_figure(sel, {"with spread (the rule's configuration)": f"{w60} and rescale == 'none' and spread",
                                 "no spread (s = 0)": f"{w60} and rescale == 'none' and not spread",
                                 "with spread, time rescaled by the pooled-fit speed": f"{w60} and rescale == 'fit'"}, R)
    save(fig, P.STEP3B2_FIG_SELECTION, dpi=110)
    plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    def pivot(q):
        g = sel.query(q)
        g = g.assign(cell=[f"{'**' if r.role == 'true' else ''}{r.n_selected}/{r.replicates}{'**' if r.role == 'true' else ''}"
                           for r in g.itertuples()])
        return g.pivot(index="truth", columns="term", values="cell").reindex(columns=TERMS).reset_index()

    def returned(q):
        """Per truth: the term sets that came out and how often."""
        g = fits.query(q)
        return pd.DataFrame([{"truth": t, "returned (x replicates)": ", ".join(
            f"{k} x{v}" for k, v in x.selected.value_counts().items()),
            "median threshold": x.lam.median()} for t, x in g.groupby("truth")])

    def what(g, by):
        """Per group: the term sets that came out, with the number of replicates."""
        return g.groupby(by).selected.agg(lambda v: ", ".join(f"{k} x{n}" for k, n in v.value_counts().items())).rename(
            "returned (x replicates)").reset_index()

    def rule_grid(spread):  # noqa: E306
        g = rule[rule.spread == spread].assign(cell=lambda d: np.where(d.met, "met", d.why_not))
        return g.pivot(index=["window", "offset", "rescale"], columns="truth", values="cell").reset_index()

    coef_main = coefs[(coefs.window == M["window"])].drop(columns=["window"])
    coef_show = coef_main[(coef_main.offset == M["offset"])].drop(columns=["offset"])
    lam_tab = fits.groupby(["spread", "window", "offset", "rescale"]).lam.agg(["median", "min", "max"]).reset_index()
    n_new = sum(r["pooled"]["n_new_fits"] for r in results if r["pooled"])
    n_fits = sum(len(S1C["t2_starts"][c]) + (c != "F1") for c in S1C["candidates"]) * sum(r["spread"] for r in results)
    rows_tab = fits.groupby(["window", "offset"]).agg(rows=("n_rows", "median"), tails=("n_tails_used", "median")).reset_index()
    report = f"""# WAA_estm step 3b-2: weak-form SINDy on fake tails (known truth)

Written by `research/WAA_estm/step3b2_sindy_fake.py` (libraries `waa_sindy.py`, `waa_sim.py`, `waa_pooled.py`; settings
`waa.sindy3b2`). Fake data only: no real tail is fitted and the frozen tail set v1 is not changed.

**Verdicts: Q1 {v1}, Q2 {v2}, Q3 {v3}, Q4 {v4}** (section 8).

## 1. What was run

- **Software**: pysindy is not installed in this environment (neither is scikit-learn, which it needs), so the
  weak-form SINDy is written out in `waa_sindy.py` (numpy {np.__version__} only, about 150 lines): weak-form rows,
  thresholded sparse regression, cross-validation over events, stability selection.
- **Data**: the fake datasets of step 3b-1c, test T2 (3a-fix noise, event bootstrap of the {len(skel)}-tail TRAIN v1
  skeleton, seeds {S['seed_base']} + replicate), {R} replicates per truth: {', '.join(S['truths'])} with spread
  s = {S['s_true']}, and {', '.join(plan['nospread_truths'])} with no spread (s = 0; same seeds, so the same noise and
  offsets). {sum(reuse.values())} of the {len(reuse)} spread datasets that step 3b-1c also made are bit-identical to
  its datasets (m3 was not a truth there). Windows: {', '.join(plan['windows'])}.
- **Runs**: {len(fits)} SINDy runs on {len(results)} datasets (window x offset mode x rescaling). Pooled fits for the
  {sum(r['spread'] for r in results)} spread datasets (5 candidates, 4 starts each, first 60 min): {n_fits - n_new} reused from
  step 3b-1c, {n_new} new (m3).
- **Cut for time**: {', '.join(plan['cuts']) if plan['cuts'] else 'nothing'}. The timed replicate (m3 with spread) took
  {plan['timed_pooled_s']:.0f} s for the pooled fits and {plan['timed_sindy_s']:.0f} s for its SINDy runs; projected
  total {plan['projected_s'] / 60:.0f} min on {workers} processes (budget {S['time_budget_h']} h).
- **Run time**: {plan['wall_s'] / 60:.0f} min wall time (cached in `.cache/sindy3b2/`). This run: {time.time() - clock:.0f} s.

## 2. Library, settings and the offset

- **Law searched**: dA/dt = sum of xi_j theta_j(A), one law for all tails of a dataset. Library ({len(TERMS)} terms):
  1, A, A^2, A^3 and the switch terms max(A - A_c, 0) for A_c = {', '.join(map(str, GRID))} dB (written `sw0.1` ...).
- **Weak form**: test function phi = (1 - tau^2)^{S['test_power']} on windows of {2 * S['half_width_min'] + 1} min, centres every
  {S['stride_min']} min; row: -integral(phi' A) = sum xi_j integral(phi theta_j(A)). NaN runs up to {S['max_gap_min']} min
  are interpolated, a window with a NaN left is skipped. Rows of a tail are weighted by 1 / (1 min noise sd of its
  sub-link, from `sim_noise_params.csv`). Median rows and tails per dataset:

{md(rows_tab)}

- **Sparse regression**: thresholded backward elimination (`waa_sindy.select`). The usual STLSQ (drop small
  coefficients) failed a noise-free self-test: while A > 0.1 dB the term max(A - 0.1, 0) is A minus a constant, the
  fit splits the coefficient between such copies and none is small. So a term is judged by how much the fit worsens
  without it (sqrt of the increase of the squared residual / ||b||), and among the terms that cost less than the
  threshold the least simple one (last in the library order) is removed first.
- **Threshold**: grid {S['lambda_grid']}, chosen per run by {S['cv_folds']}-fold cross-validation over EVENT DRAWS
  (all tails of an event draw are held out together), one-standard-error rule (largest threshold within one SE of
  the best held-out error). Chosen thresholds (median, min, max):

{md(lam_tab)}

- **Stability selection**: {S['stability']['resamples']} bootstrap resamples of the event draws at the chosen threshold; a
  term counts as **selected in a replicate if its frequency is >= {THR}**. Coefficients: ordinary least squares on the
  selected terms.
- **Offset** (A = W + c, c unknown per tail), three ways:
  - `naive`: ignore c (use A as W).
  - `handled` (my way, the configuration of the rule): a free constant per tail in dA/dt, removed by subtracting each
    tail's mean from all its weak-form rows. For a linear law this is exact (dA/dt = -k (A - c) = -k A + k c). For
    a nonlinear law it is not (W^2 = A^2 - 2 c A + c^2: a per-tail slope error is left). Price: the library's
    constant term `1` cannot be seen any more, and the absolute level of A (where W = 0 is) is given up.
  - `oracle`: the true c subtracted. Impossible on real data; it shows the ceiling.
- **Rescaling check** (spread datasets, handled offset): right-hand side divided by the tail's speed factor f_i
  (f_i dA/dt = h(A)): `fit` = f_i estimated by the BIC-winning pooled fit of the same dataset (posterior mode per
  tail), `true` = the simulator's f_i.
- **True terms of the rule** (`waa.sindy3b2.decision.true_terms`): m1: A; m3: A^2; m5: A^2 (A allowed, its coefficient
  reported but not counted); m6: A and a switch term at 0.1 or 0.2 dB (true A_c = {par['m6']['A_c']:.2f} dB). Seen in
  the noise-free self-test before the run: above A_c the switch law IS a straight line, dW/dt = k_run A_c -
  (k_run + k_evap) W, so on tails that stay above A_c the terms `1` + `A` describe it exactly.

## 3. Selection frequency: truth x term (`{rel(P.STEP3B2_SELECTION)}`)

Replicates (of {R}) in which the term is selected; **bold** = a term of the true law. First 60 min, handled offset.

**With spread (the configuration of the rule):**

{md(pivot(main_q + " and spread"))}

{md(returned(main_q + " and spread"))}

**No spread (s = 0):**

{md(pivot(main_q + " and not spread"))}

{md(returned(main_q + " and not spread"))}

**With spread, time rescaled by the pooled-fit speed:**

{md(pivot(f"{w60} and rescale == 'fit'"))}

**With spread, time rescaled by the TRUE speed:**

{md(pivot(f"{w60} and rescale == 'true'"))}

![selection](../figures/{P.STEP3B2_FIG_SELECTION.name})

The pre-set rule in EVERY configuration (cell = `met`, or why not). Only the first-60-min, handled, not rescaled row
with spread counts for Q1; the others show what helps and what hurts.

With spread:

{md(rule_grid(True))}

No spread:

{md(rule_grid(False))}

**Post hoc diagnostic, added after the first complete run (it changes no verdict).** In the rule's configuration
the one-standard-error rule picked the largest threshold of the grid in most runs. Held-out error (share of the
squared left-hand side that is NOT explained; medians over the replicates) with the rule's threshold, with the
threshold of lowest error, and with the largest threshold of the grid:

{md(cv_tab)}

What comes out, first 60 min, handled offset, if the threshold of LOWEST held-out error is used instead:

{md(what(alt, ['spread', 'rescale', 'truth']))}

The rule evaluated on that variant (for information only):

{md(alt_rule[['spread', 'rescale', 'truth', 'true_selected', 'spurious', 'coef_ratios', 'met']])}

## 4. Coefficient recovery (`{rel(P.STEP3B2_COEFS)}`)

fitted / true (median, min, max over the replicates; 0 where the term is not selected), first 60 min, handled
offset. True = the base law (speed factor 1). `counted` = part of the rule.

{md(coef_show)}

## 5. m4 (hidden state) and the spread mixture

**m4 (two layers), first 60 min:**

{md(what(fits[(fits.truth == 'm4') & (fits.window == M['window'])], ['offset', 'rescale']))}

**m1 (exponential), first 60 min: no spread vs spread vs rescaled:**

{md(what(fits[(fits.truth == 'm1') & (fits.window == M['window'])], ['spread', 'offset', 'rescale']))}

## 6. Offset: naive vs handled vs oracle (first 60 min, not rescaled)

{md(what(fits[(fits.window == M['window']) & (fits.rescale == 'none')], ['spread', 'truth', 'offset']))}

Coefficients (fitted / true, median):

{md(coef_main[coef_main.rescale == 'none'].pivot(index=['spread', 'truth', 'coefficient'], columns='offset', values='median_ratio').reset_index())}

## 7. Agreement with the pooled fit (`{rel(P.STEP3B2_AGREEMENT)}`)

Per spread dataset: the BIC winner of the pooled fit (best of 4 starts, first 60 min) and the family named by the
SINDy terms (first 60 min, handled offset). Families: exponential (F1; only A), nonlinear (F3 / F5; A^2 or A^3),
switch (F6; a switch term), two layers (F4; no SINDy family). Counted for Q4: m1, m3, m5, m6.

{md(agree.groupby('truth').agg(datasets=('agree', 'size'), pooled_winners=('pooled_winner', lambda v: ', '.join(f"{k} x{n}" for k, n in v.value_counts().items())), sindy_families=('sindy_family', lambda v: ', '.join(f"{k} x{n}" for k, n in v.value_counts().items())), agree=('agree', 'sum'), sindy_rescaled_families=('sindy_rescaled_family', lambda v: ', '.join(f"{k} x{n}" for k, n in v.value_counts().items())), agree_rescaled=('agree_rescaled', 'sum')).reset_index())}

Linear rate (1 / min): 1 / tau of the pooled F1 fit vs the SINDy coefficient of A (0 = A not selected), medians per truth:

{md(agree.groupby('truth')[['pooled_F1_rate', 'sindy_A_coef', 'sindy_rescaled_A_coef']].median().round(4).reset_index())}

## 8. Verdicts and what they mean for Step 5

{md(verdicts)}

Rule per truth (first 60 min, handled offset, with spread, not rescaled):

{md(q1[['truth', 'true_selected', 'spurious', 'coef_ratios', 'true_ok', 'spurious_ok', 'coef_ok', 'met']])}

{MEANING}

## 9. Unclear / possibly wrong

{UNCLEAR}
"""
    save(report, P.STEP3B2_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {time.time() - clock:.0f} s")


# Written after the run of 2026-10-11 (the tables are computed on every run; these two texts are the reading of that
# run and quote its numbers: update them if a rerun with other settings changes the tables).
MEANING = """**In plain words: on tails like ours, this weak-form SINDy does not find the drying law.** It returns "A" (an
exponential) for almost every truth, and where it adds a term the term is wrong.

- **Q1 (true terms), FAIL: 0 of 4 truths.** In the rule's configuration m3 and m5 (quadratic laws) come out as a
  plain "A" in 5 of 5 replicates and A^2 is never selected. m6 (switch) comes out as "A" with a rate of about
  0.05 / min, which is k_run + k_evap: above A_c the switch law is a straight line, and the free constant per tail
  swallows its constant part; no switch term is ever selected. m1 keeps A (5/5) but gets a false A^2 in 3 of 5 and
  its rate is x1.79. Without spread only m1 is met (A in 5/5, rate x0.94); m3 and m5 still come out as "A" (4/5).
- **Why.** About half of the left-hand side is noise, and the laws differ very little in what they explain: the
  held-out error is 0.36-0.49 without spread and 0.53-0.58 with spread whatever the threshold, and the best
  threshold beats the one-term model by 0.00-0.06, at most about 1.5 standard errors (0.03-0.04). So cross-validation
  cannot tell the laws apart and keeps one term. With the threshold of lowest error instead (post hoc check), the
  result is not better but unstable: 1 to 8 terms, different in every replicate.
- **The offset matters.** With the true offset removed and no spread, the first hour gives A^2 for m3 and m5 in 5 of
  5 (coefficients x0.83 and x0.77) and "A" for m1 (x0.96): the method works when W is known. Ignoring the offset
  (naive) puts a false constant term into nearly every result (m1 without spread: "1" in 5/5, rate x0.69). The free
  constant per tail repairs the linear law only; for the quadratic laws it removes the very information that
  separates A^2 from A in the first hour.
- **Q2 (m4, hidden state), FAIL = misleading.** m4 comes out as a clean single-variable law in 5 of 5 ("A" x3,
  "A + A^2" x2), the same answers as m1 and m6. Nothing in the output says that a second state is missing.
- **Q3 (mixture trap), PARTIAL = the trap is real.** Exponential tails with different speeds return a false A^2 in
  3 of 5 replicates (0 of 5 without spread), and the pooled rate is x1.79 the median tail's. Rescaling time by the
  pooled-fit speed does not cure it (A^2 in 4 of 5, rate x0.48); even the TRUE speed leaves A^2 in 2 of 5 and a
  rate x0.28, so the rescaling as I built it has a problem of its own (section 9, point 3).
- **Q4 (agreement with the pooled fit), FAIL.** Same family in 1 of 20 datasets. The pooled fit names the true
  family in 19 of 20 of them (m1: F1 x4; m3: F3 / F5 x5; m5: F5 / F3 x5; m6: F6 x5), so the disagreement is SINDy
  being wrong. With time rescaling they agree in 11 of 20, but only because rescaled SINDy says "nonlinear" for
  every truth, m1, m4 and m6 included.

**For Step 5 (SINDy on real tails).**
1. Do not run it as planned (one pooled law dA/dt = h(A) on the real tails) and read the result as the drying law.
   On fake data the answer is "A" whatever the truth, so "A" on real tails would not mean exponential, and an A^2
   term could be the speed mixture (Q3), not physics.
2. The pooled fit of 3b-1b / 3b-1c is the tool that works here. It wins because it uses what this SINDy throws
   away: the absolute level (offset prior), a speed per tail inside the model, and the noise model.
3. If SINDy stays in the plan, it needs a redesign that is tested on these same fake datasets first: the speed per
   tail inside the regression (not a pooled law), an offset prior instead of a free constant, a library without
   switch terms below the noise level, and a treatment of the noise inside the library terms. Until then it can
   serve only as a cross-check of the pooled-fit winner (does the fitted law's h(W) lie inside the SINDy scatter?).
4. Step 7 (wetting, with rain as input, whole events) is a different problem: the signal there is the rise during
   rain, several dB in minutes, far above the noise. This result does not say SINDy will fail there; step 3c should
   test it the same way."""
UNCLEAR = """1. **The sparse regression is mine, not the textbook one.** STLSQ failed the noise-free self-test (duplicate-like
   switch terms), so I used thresholded backward elimination that removes the least simple removable term first.
   That prefers polynomials over switch terms on ties. It may have helped "A" to win, but the post hoc check with a
   low threshold shows no hidden correct answer either.
2. **Thresholds sit at the edge of the grid** (0.45 in most runs of the rule's configuration). The one-SE rule with
   5 folds is conservative; a larger grid would not change the outcome (one term is already the result), and the
   lowest-error threshold was checked post hoc. That check was added after seeing the first results; it reran the
   same seeds and the decision table is identical to the first run (compared file to file).
3. **The time rescaling is probably flawed as built.** Rows are divided by the speed factor f on the right-hand
   side, so very fast tails (f down to 0.05), whose signal is gone after a few minutes, enter with noise-only
   library terms multiplied by up to 20. I think this is why even the true speed gives rates x0.28-0.55 (noise in
   the regressors pulls coefficients toward zero). A fair rescaling test needs another weighting or resampling
   each tail on its own time axis. So "rescaling does not cure the mixture" is not established; "this rescaling
   does not" is.
4. **Noise inside the library terms** (A^2 of a noisy A is biased up by the noise variance; all terms carry the
   noise of A) is not treated at all. Rates are low by 4-25 % even in the best cases (m1 x0.94, m3 x0.83, m5 x0.77).
5. **"handled" is my choice of offset treatment.** An offset prior as in the pooled fit (c near 0 within
   offset_std) would keep the level information that the free constant gives up; I did not build it because the
   constant it implies in dA/dt depends on the unknown law. It is the first thing to try in a redesign.
6. **The rule for m6 could not be met by construction in the first hour**: seen before the run in the noise-free
   self-test (above A_c the law is "1 + A"); I kept the rule as written. In the naive mode m6 does come out as
   "1 + A" in 4 of 5 with spread, but so do m1, m3, m4 and m5.
7. **m5's A is "allowed, not required"** and its coefficient not counted (my reading of the rule). It is selected
   5 of 5 with a rate x17: it is standing in for the quadratic term, not measuring a.
8. **Q2 and Q3 verdict rules are mine** (fixed in the config before the run): Q2 FAIL = a clean single-variable law
   in >= 4 of 5; Q3 PASS = a false term in at most 1 of 5, FAIL = in >= 4 of 5. Q4 counts families, with F3 and F5
   as one, and leaves m4 out.
9. **One window width (21 min) and one test function.** Wider windows average more noise but blur the fast first
   minutes; not explored (no re-running after the result).
10. **Whole-tail window**: about 34,000 rows per dataset, most of them late minutes with A near 0 (noise only). It
    helps the no-spread quadratic laws (m5 met, m3 A^2 x0.65) and hurts m1 (false A^2 in 3 of 5 without spread).
    I did not study why.
11. **5 replicates**, of very different size (135-753 tails): counts such as 2/5 vs 3/5 are not distinguishable.
12. **Key figure**: the limit of 6 is reached, so `step3b2_selection_frequency.png` is diag. I do not think a
    negative method test must be key; if the learners want it, the candidate to demote is
    `outputs/tables/step3b1b_recovery.csv` (3b-1c's recovery table covers the same question under harder conditions).
"""


if __name__ == "__main__":
    main()
