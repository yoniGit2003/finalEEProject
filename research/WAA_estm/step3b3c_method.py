"""WAA_estm step 3b-3c: finish the discovery method on fake tails before real tails. No real tail is fitted.

Five parts, each with its own rule fixed before the run (`waa.curve3b3c` in research/config.yaml):
  A  Jump 2 alone: the TRUE h(W) = W k(W) of m1, m3, m5, m6 plus 2 % / 5 % noise; does best-subset + BIC (and
     ensemble SINDy) name the true terms?
  B  Jump 1 with an optional step: free curve x smooth step at a fitted W_c (waa_curve.StepCurve); kept only if BIC
     improves by > 6. Must find the switch of m6 and no step for m1, m3, m5.
  C  Jump 3: every candidate (the jump-2 formulas, the step law, the menu laws F1, F3, F4, F5, F6) is refitted in the
     pooled model on half of the event draws and scored by the log-likelihood of the other half (both ways).
  D  The band for real data: one dataset only, so the band comes from an event bootstrap of the fit. Does it cover
     the truth?
  E  F4 with the F1-solution start on the m4 datasets (information for the failed m4 warning of 3b-3b).
Data of B-E: the cached fake datasets and free-curve fits of step 3b-3b (.cache/fits3b3b/, seeds 4000 + replicate,
first 60 min). Real data is used only through those (NaN minutes, noise statistics). The frozen tail set v1 is not read.

Replicate 0 is run first and timed; the projection decides the cuts. Everything is cached in .cache/fits3b3c/, so a
restart resumes. Delete that folder to redo. Run from the repository root (needs step 3b-3b run first):
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step3b3c_method.py --pysr   # optional: part A with PySR (15 min)
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step3b3c_method.py          # everything else and the report
Part A uses the packages pysindy and (with --pysr) pysr, if installed; without them those rows are left out.
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
import step3b1c_stress_test as st                          # noise tables, menu fits from several starts
import step3b3_curve_then_equation as s3                   # load_tails, true_k, md, rel

CFG = s3.CFG
S, S3B, FIT, S1C = (CFG["waa"][k] for k in ("curve3b3c", "curve3b3b", "curve3b3", "stress3b1c"))
A, B, C, D_ = S["a"], S["b"], S["c"], S["d"]
W_LO = S3B["w_lo_db"]
TERMS_A = wc.curve_terms(A["switch_grid_db"])              # library of part A (and of jump 2 when a step is kept)
TERMS_HI = wc.curve_terms(S3B["j2_switch_grid_db"])        # library of 3b-3b (judged range above 0.35 dB)
INK, INK2 = "#1a1a1a", "#52514e"
CACHE = P.CACHE / "fits3b3c"
CACHE_3B = P.CACHE / "fits3b3b"
FAMILY = {"F1": "exponential", "F3": "quadratic", "F5": "quadratic", "F4": "two layers", "F6": "switch", "step": "switch"}
md, rel = s3.md, s3.rel


def formula_family(terms):
    """Family of a jump-2 formula: 'W' = exponential; 'W^2' or 'W + W^2' = quadratic; W + switch term(s) at <= 0.3 dB
    = switch; anything else = other."""
    t = set(terms)
    if t == {"W"}:
        return "exponential"
    if t in ({"W^2"}, {"W", "W^2"}):
        return "quadratic"
    sw = [x for x in t if x.startswith("sw")]
    if "W" in t and sw and len(sw) == len(t) - 1 and all(float(x[2:]) <= 0.3 for x in sw):
        return "switch"
    return "other"


# =================================================================================================== workers
_INPUTS = {}


def _load(job):
    if not _INPUTS:
        _, _INPUTS["noise"], _INPUTS["offset_sd"] = st.load_noise()
    tails = s3.load_tails(CACHE_3B / f"sim_main_{job['truth']}_{job['rep']}.parquet")
    r3 = json.loads((CACHE_3B / f"job_main_{job['truth']}_{job['rep']}.json").read_text())
    return tails, r3, _INPUTS["noise"], _INPUTS["offset_sd"]


def _saver(f_, res):
    def save():
        tmp = f_.with_suffix(".tmp")
        tmp.write_text(json.dumps(res))
        tmp.replace(f_)
    return save


def run_bc(job):
    """Parts B, C (and E) for ONE 3b-3b dataset: the step fit, the jump-2 formulas, and every candidate refitted on
    the whole dataset (BIC) and on each half of the event draws (held-out log-likelihood of the other half)."""
    f_ = CACHE / f"bc_{job['truth']}_{job['rep']}.json"
    res = json.loads(f_.read_text()) if f_.exists() else {}
    if res.get("done"):
        return res
    save = _saver(f_, res)
    t0 = time.time()
    tails, r3, noise, offset_sd = _load(job)
    Dfull = wp.pack(tails, noise, offset_sd, S3B["window_min"], FIT["min_points"])
    used = set(Dfull["tail_id"])
    n_ind = len({(tl["cml_id"], tl["event_copy"]) for tl in tails if tl["tail_id"] in used})
    # the same two halves of the event draws as the cross-validation of 3b-3b
    evs = np.random.default_rng([S3B["seed_base"], job["rep"]]).permutation(sorted({tl["event_copy"] for tl in tails}))
    Dh = [wp.pack([tl for tl in tails if tl["event_copy"] in h], noise, offset_sd, S3B["window_min"], FIT["min_points"])
          for h in (set(evs[::2].tolist()), set(evs[1::2].tolist()))]
    K, x3, lam = len(r3["knots"]), np.asarray(r3["x"]), r3["lam_star"]
    plain = wc.Curve(r3["knots"], prior_sd=FIT["prior_sd_logk"])
    sc = wc.StepCurve(r3["knots"], B["wc_range_db"], B["width_logw"], prior_sd=FIT["prior_sd_logk"])
    res.update({"truth": job["truth"], "rep": job["rep"], "n_independent": n_ind, "w_hi": r3["w_hi"]})

    def held_out(fit_half, score):
        return float(sum(score(fit_half(Dh[a]), Dh[b]) for a, b in ((0, 1), (1, 0))))

    # ---- B: the step, from the 3b-3b solution (same knots, same penalty), one start per W_c
    kept = False
    if job["step"]:
        res.setdefault("step_starts", {})
        for wc0 in B["wc_starts_db"]:
            if str(wc0) not in res["step_starts"]:
                p0 = np.r_[x3[:-1], sc.u_of(wc0, B["height_start"]), x3[-1]]
                res["step_starts"][str(wc0)] = sc.fit(Dfull, p0, lam, B["maxfev"]["start"])
                save()
        if "step" not in res:
            best = min(res["step_starts"], key=lambda n: res["step_starts"][n]["cost"])
            r = sc.fit(Dfull, res["step_starts"][best]["x"], lam, B["maxfev"]["finish"])
            wc_db, height = sc.step(r["x"][:-1])
            r.update({"from": best, "wc_db": wc_db, "height": height, "bic": -2 * r["ll"] + (K + 3) * np.log(n_ind)})
            res["step"] = r
            save()
        kept = bool(r3["bic"] - res["step"]["bic"] > B["bic_gain"])
    res["step_kept"] = kept
    res.setdefault("seconds_b", round(time.time() - t0, 1))

    # ---- jump 2 on the chosen curve of this dataset: the best subset of each size
    if "j2" not in res:
        if kept:
            w, terms = np.geomspace(B["step_range_lo_db"], r3["w_hi"], FIT["grid_points"]), TERMS_A
            k = sc.k(res["step"]["x"][:-1], w)
        else:
            w, terms = np.geomspace(W_LO, r3["w_hi"], FIT["grid_points"]), TERMS_HI
            k = plain.k(x3[:-1], w)
        subs = wc.best_subset(w, w * k, terms, max(C["j2_sizes"]))
        res["j2"] = {"range_db": [float(w[0]), float(w[-1])], "bic_pick": subs[0]["terms"],
                     "by_size": {str(n): next({"terms": s_["terms"], "coef": s_["coef"]} for s_ in subs if s_["size"] == n)
                                 for n in C["j2_sizes"]}}
        save()

    # ---- C: every candidate on the whole dataset and on the halves
    res.setdefault("cand", {})
    for n in C["j2_sizes"]:
        name = f"J2-{n}"
        if name not in res["cand"]:
            sub = res["j2"]["by_size"][str(n)]
            fl = wc.FormulaLaw(sub["terms"], sub["coef"])
            mf = C["maxfev_formula"]
            full = fl.fit(Dfull, np.r_[np.ones(n), x3[-1]], 0.0, mf["full"], step=0.3)
            held = held_out(lambda Da: np.asarray(fl.fit(Da, full["x"], 0.0, mf["half"], step=0.3)["x"]), fl.loglik)
            res["cand"][name] = {"what": " + ".join(sub["terms"]), "terms": sub["terms"],
                                 "coef": [float(v) for v in fl.coef(full["x"][:-1])], "s": float(np.exp(full["x"][-1])),
                                 "ll": full["ll"], "k": n + 1, "held": held}
            save()
    if kept and "step" not in res["cand"]:
        xs = res["step"]["x"]
        held = held_out(lambda Da: np.asarray(sc.fit(Da, xs, lam, C["maxfev_half"])["x"]), sc.loglik)
        res["cand"]["step"] = {"what": f"free curve with a step at {res['step']['wc_db']:.2f} dB", "ll": res["step"]["ll"],
                               "k": K + 3, "held": held}
        save()
    if "menu_done" not in res:
        st.CACHE = CACHE                                   # the menu fits of this step are cached in fits3b3c/
        starts = {}
        for c in C["menu"]:
            starts[c] = {k: [v[0], v[1], S1C["maxfev_generic"]] for k, v in S1C["t2_starts"][c].items()}
            if c != "F1":
                starts[c]["f1sol"] = [None, None, S1C["maxfev_generic"]]
        rows = st.fit_dataset({"test": "c", "truth": job["truth"], "rep": job["rep"], "starts": starts,
                               "sim": str(CACHE_3B / f"sim_main_{job['truth']}_{job['rep']}.parquet")})
        for c in C["menu"]:
            fits = [r for r in rows if r["candidate"] == c]
            best = max(fits, key=lambda r: r["ll"])
            gen = max((r for r in fits if r["start"] != "f1sol"), key=lambda r: r["ll"])
            held = held_out(lambda Da: wp.fit(c, Da, best["theta"], best["s"], C["maxfev_half"]),
                            lambda r, Db: wp.pooled_loglik(c, r["theta"], r["s"], Db))
            res["cand"][c] = {"what": ", ".join(f"{n_} {v:.4g}" for n_, v in zip(wp.NAMES[c], best["theta"])),
                              "ll": best["ll"], "k": len(wp.NAMES[c]) + 1, "held": held, "start": best["start"],
                              "ll_generic": gen["ll"], "start_generic": gen["start"],
                              "ll_by_start": {r["start"]: r["ll"] for r in fits}}
            save()
        res["menu_done"] = True
    res["seconds"] = round(res.get("seconds", 0.0) + time.time() - t0, 1)
    res["done"] = True
    save()
    return res


def run_d(job):
    """Part D for ONE 3b-3b dataset: the free curve refitted on event bootstraps of that dataset (same knots and
    penalty, start = the fit of the dataset itself)."""
    f_ = CACHE / f"d_{job['truth']}_{job['rep']}.json"
    res = json.loads(f_.read_text()) if f_.exists() else {"truth": job["truth"], "rep": job["rep"], "x": [], "seconds": 0.0}
    if len(res["x"]) >= job["bootstraps"]:
        return res
    save = _saver(f_, res)
    tails, r3, noise, offset_sd = _load(job)
    by_ev = {}
    for tl in tails:
        by_ev.setdefault(tl["event_copy"], []).append(tl)
    evs = sorted(by_ev)
    cv = wc.Curve(r3["knots"], prior_sd=FIT["prior_sd_logk"])
    for b in range(len(res["x"]), job["bootstraps"]):
        t0 = time.time()
        rng = np.random.default_rng([D_["seed"], S["d"]["truths"].index(job["truth"]), job["rep"], b])
        draw = [tl for e in rng.choice(evs, len(evs), replace=True) for tl in by_ev[e]]
        Db = wp.pack(draw, noise, offset_sd, S3B["window_min"], FIT["min_points"])
        res["x"].append(cv.fit(Db, r3["x"], r3["lam_star"], D_["maxfev"])["x"])
        res["seconds"] = round(res["seconds"] + time.time() - t0, 1)
        save()
    return res


def run_any(job):
    return run_d(job) if job["part"] == "D" else run_bc(job)


# =================================================================================================== part A
def part_a(par):
    """Jump 2 alone on the true h(W) plus noise. Returns (per-draw table, per-cell table, software note)."""
    w = np.geomspace(*A["w_range_db"], A["grid_points"])
    methods = ["best-subset + BIC", "E-SINDy (waa_sindy, bagged)"]
    try:
        import pysindy
        methods.append("E-SINDy (pysindy STLSQ, bagged)")
        note = f"pysindy {pysindy.__version__} installed (STLSQ used inside the bagging loop)"
    except Exception as e:                                 # noqa: BLE001
        note = f"pysindy not usable ({type(e).__name__}); E-SINDy with waa_sindy.py only"
    true = a_truths(par)
    sw6 = set(A["m6_switch_terms"])
    E = A["esindy"]
    rows = []
    for i, truth in enumerate(A["truths"]):
        h = w * s3.true_k(truth, par, w)
        for j, sig in enumerate(A["noise"]):
            for d in range(A["draws"]):
                rng = np.random.default_rng([A["seed"], i, j, d])
                hn = h * (1 + sig * rng.standard_normal(len(w)))
                for m in methods:
                    if m.startswith("best"):
                        b = wc.best_subset(w, hn, TERMS_A, A["max_terms"])[0]
                        got = dict(zip(b["terms"], b["coef"]))
                    else:
                        inc, med = wc.bagged_select(w, hn, TERMS_A, E["threshold"], E["bootstraps"], rng,
                                                    "pysindy" if "pysindy" in m else "waa_sindy")
                        got = {t: float(c) for t, c, p_ in zip(TERMS_A, med, inc) if p_ >= E["inclusion"]}
                    sel = set(got)
                    if truth == "m6":
                        exact = "W" in sel and bool(sel & sw6) and sel <= {"W"} | sw6
                        incl = "W" in sel and bool(sel & sw6)
                        coef = {"W": got.get("W", 0.0), "switch": sum(got.get(t, 0.0) for t in sw6)}
                    else:
                        exact, incl = sel == set(true[truth]), set(true[truth]) <= sel
                        coef = {t: got.get(t, 0.0) for t in true[truth]}
                    rows.append({"truth": truth, "noise": sig, "draw": d, "method": m,
                                 "returned": " + ".join(t for t in TERMS_A if t in sel) or "(nothing)",
                                 "exact": bool(exact), "true_included": bool(incl),
                                 **{f"ratio_{c}": coef[c] / true[truth][c] for c in coef}})
    extra, pysr_s = pysr_rows(par)
    if extra:
        methods.append("PySR")
        note += f"; PySR: {len(extra)} fits ({pysr_s / 60:.0f} min in a separate process)"
    else:
        note += "; PySR not run"
    draws = pd.DataFrame(rows + extra)
    cells = []
    for (truth, sig, m), g in draws.groupby(["truth", "noise", "method"], sort=False):
        ok = g[g.exact]
        ratios = {c: float(ok[f"ratio_{c}"].median()) if len(ok) else np.nan for c in true[truth]}
        cells.append({"truth": truth, "noise": sig, "method": m, "exact": int(g.exact.sum()),
                      "true_included": int(g.true_included.sum()), "draws": len(g),
                      "most_common": f"{g.returned.value_counts().index[0]} x{g.returned.value_counts().iloc[0]}",
                      "coef_ratios": ", ".join(f"{c} x{v:.2f}" for c, v in ratios.items()),
                      "ok": bool(g.exact.sum() >= A["min_draws"]
                                 and all(abs(v - 1) <= A["coef_tol"] for v in ratios.values()))})
    return draws, pd.DataFrame(cells), note, methods


def a_truths(par):
    return {"m1": {"W": 1 / par["m1"]["tau"]}, "m3": {"W^2": par["m3"]["c"]},
            "m5": {"W": par["m5"]["a"], "W^2": par["m5"]["b"]},
            "m6": {"W": par["m6"]["k_evap"], "switch": par["m6"]["k_run"]}}


def pysr_terms(expr, w, h):
    """Terms of a PySR formula: expanded with sympy; a term below `drop_share` of h everywhere is ignored; c W^p with p
    within `power_round` of 1, 2, 3 is named W, W^2, W^3; anything else keeps its own name. Returns {name: c}."""
    import sympy as sp
    Y = A["pysr"]
    W = sp.Symbol("W", positive=True)
    out = {}
    for term in sp.Add.make_args(sp.expand(sp.sympify(str(expr), locals={"W": W}))):
        val = np.broadcast_to(np.asarray(sp.lambdify(W, term, "numpy")(w), float), w.shape)
        if np.max(np.abs(val) / h) < Y["drop_share"]:
            continue
        c, rest = term.as_coeff_Mul()
        if rest == 1:
            name = "const"
        elif rest == W:
            name = "W"
        elif rest.is_Pow and rest.base == W and rest.exp.is_number:
            p_ = float(rest.exp)
            near = abs(p_ - round(p_)) <= Y["power_round"] and round(p_) in (1, 2, 3)
            name = {1: "W", 2: "W^2", 3: "W^3"}[round(p_)] if near else f"W^{p_:.2f}"
        else:
            name = str(rest)
        out[name] = out.get(name, 0.0) + float(c)
    return out


def part_a_pysr(par):
    """Part A with PySR (symbolic regression) on the same noise draws. Run as a separate process:
    `step3b3c_method.py --pysr`; every fit is cached in .cache/fits3b3c/pysr.json."""
    from pysr import PySRRegressor
    Y = A["pysr"]
    f_ = CACHE / "pysr.json"
    done = json.loads(f_.read_text()) if f_.exists() else {}
    w = np.geomspace(*A["w_range_db"], A["grid_points"])
    for i, truth in enumerate(A["truths"]):
        h = w * s3.true_k(truth, par, w)
        for j, sig in enumerate(A["noise"]):
            for d in range(A["draws"]):
                key = f"{truth}_{sig}_{d}"
                if key in done:
                    continue
                hn = h * (1 + sig * np.random.default_rng([A["seed"], i, j, d]).standard_normal(len(w)))
                t0 = time.time()
                model = PySRRegressor(niterations=Y["niterations"], maxsize=Y["maxsize"], timeout_in_seconds=Y["timeout_s"],
                                      binary_operators=["+", "-", "*", "/", "^"], unary_operators=["exp"],
                                      constraints={"^": (-1, 1)}, parallelism="serial", progress=False, verbosity=0,
                                      output_directory=str(P.CACHE / "pysr"))
                for attempt in range(3):                   # Windows sometimes locks PySR's checkpoint file: try again
                    try:
                        model.fit(w[:, None], hn, weights=1 / hn ** 2, variable_names=["W"])
                        break
                    except PermissionError:
                        if attempt == 2:
                            raise
                        time.sleep(2)
                done[key] = {"expr": str(model.sympy()), "seconds": round(time.time() - t0, 1)}
                tmp = f_.with_suffix(".tmp")
                tmp.write_text(json.dumps(done))
                tmp.replace(f_)
                print(f"   PySR {key}: {done[key]['expr']} ({done[key]['seconds']} s)", flush=True)


def pysr_rows(par):
    """The cached PySR answers as rows of the part A draw table (empty list if PySR was not run)."""
    f_ = CACHE / "pysr.json"
    if not f_.exists():
        return [], 0.0
    done, true = json.loads(f_.read_text()), a_truths(par)
    w = np.geomspace(*A["w_range_db"], A["grid_points"])
    rows = []
    for truth in A["truths"]:
        h = w * s3.true_k(truth, par, w)
        for sig in A["noise"]:
            for d in range(A["draws"]):
                r = done.get(f"{truth}_{sig}_{d}")
                if r is None:
                    continue
                try:
                    got = pysr_terms(r["expr"], w, h)
                except Exception:                          # noqa: BLE001  (a formula sympy cannot split)
                    got = {"unreadable": np.nan}
                exact = truth != "m6" and set(got) == set(true[truth])
                rows.append({"truth": truth, "noise": sig, "draw": d, "method": "PySR",
                             "returned": " + ".join(sorted(got)) or "(nothing)", "exact": bool(exact),
                             "true_included": bool(truth != "m6" and set(true[truth]) <= set(got)),
                             **{f"ratio_{c}": got.get(c, 0.0) / true[truth][c] for c in true[truth]}, "formula": r["expr"]})
    return rows, sum(r["seconds"] for r in done.values())


# =================================================================================================== figures
def fig_b(truth, rows, par):
    """Part B, one panel per replicate: free curve without the step (thin solid), with the step (thick solid), true
    k (dashed), the fitted W_c (vertical line with a text label)."""
    fig, ax = plt.subplots(1, len(rows), figsize=(4.1 * len(rows), 4.3), sharey=True)
    for a_, (res, r3) in zip(ax, rows):
        w = np.geomspace(0.05, r3["w_hi"], 200)
        plain = wc.Curve(r3["knots"])
        sc = wc.StepCurve(r3["knots"], B["wc_range_db"], B["width_logw"])
        a_.plot(w, plain.k(np.asarray(r3["x"])[:-1], w), "-", color=INK2, lw=1.0, label="without the step (thin solid)")
        a_.plot(w, sc.k(res["step"]["x"][:-1], w), "-", color=INK, lw=2.3, label="with the step (thick solid)")
        a_.plot(w, s3.true_k(truth, par, w), "--", color=INK, lw=1.6, label="true k(W) (dashed)")
        a_.axvline(res["step"]["wc_db"], color=INK, lw=0.8, ls=":")
        a_.text(res["step"]["wc_db"], 0.97, f" W_c = {res['step']['wc_db']:.2f} dB", transform=a_.get_xaxis_transform(),
                fontsize=8, va="top")
        gain = r3["bic"] - res["step"]["bic"]
        a_.set(xscale="log", yscale="log", xlabel="W (dB)",
               title=f"replicate {res['rep']}: BIC gain {gain:.1f}\n-> step {'KEPT' if res['step_kept'] else 'not kept'}")
        a_.title.set_fontsize(9)
        a_.grid(alpha=0.25, which="both")
        for sp in ("top", "right"):
            a_.spines[sp].set_visible(False)
    ax[0].set_ylabel("k(W) (1 / min)")
    ax[0].legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.suptitle(f"Step 3b-3c part B, truth {truth} ({ws.LAWS[truth]}): free curve with an optional step "
                 f"(kept if BIC improves by > {B['bic_gain']:g})", fontsize=10)
    fig.tight_layout()
    return fig


def fig_d(truth, g):
    """Part D, one panel per dataset: bootstrap band 5-95 % (hatched), min-max (thin dotted), the fit (solid), true k
    (dashed)."""
    reps = sorted(g.dataset.unique())
    fig, ax = plt.subplots(1, len(reps), figsize=(4.6 * len(reps), 4.3), sharey=True)
    for a_, rep in zip(np.atleast_1d(ax), reps):
        q = g[g.dataset == rep]
        a_.fill_between(q.w_db, q.k_lo, q.k_hi, facecolor="none", edgecolor="#6f6f6f", hatch="///", lw=0,
                        label=f"bootstrap band {D_['band'][0]}-{D_['band'][1]} % (hatched)")
        a_.plot(q.w_db, q.k_min, ":", color=INK2, lw=1.0, label="bootstrap min and max (dotted)")
        a_.plot(q.w_db, q.k_max, ":", color=INK2, lw=1.0)
        a_.plot(q.w_db, q.k_fit, "-", color=INK, lw=2.0, label="fit of the dataset (solid)")
        a_.plot(q.w_db, q.k_true, "--", color=INK, lw=1.6, label="true k(W) (dashed)")
        cov = float(np.mean((q.k_true >= q.k_lo) & (q.k_true <= q.k_hi)))
        a_.set(xscale="log", yscale="log", xlabel="W (dB), judged range", title=f"dataset {rep}: true k inside on {cov:.0%}")
        a_.title.set_fontsize(9)
        a_.grid(alpha=0.25, which="both")
        for sp in ("top", "right"):
            a_.spines[sp].set_visible(False)
    np.atleast_1d(ax)[0].set_ylabel("k(W) (1 / min)")
    np.atleast_1d(ax)[0].legend(frameon=False, fontsize=7.5, loc="best")
    fig.suptitle(f"Step 3b-3c part D, truth {truth} ({ws.LAWS[truth]}): band from an event bootstrap of ONE dataset",
                 fontsize=10)
    fig.tight_layout()
    return fig


def verdict(oks):
    return "PASS" if all(oks) else ("FAIL" if not any(oks) else "PARTIAL")


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    clock = time.time()
    CACHE.mkdir(parents=True, exist_ok=True)
    workers = max((os.cpu_count() or 2) - 1, 1)
    R = S3B["replicates"]
    cal = pd.read_csv(P.STEP3A_CALIBRATION).set_index("law")
    par = {law: {k: float(cal.loc[law, k]) for k in ws.PARAMS[law]} for law in ("m1", "m3", "m4", "m5")}
    m6 = pd.read_csv(P.STEP3B1_M6)
    b6 = m6.loc[m6.rmse_band_observed.idxmin()]
    par["m6"] = {k: float(b6[k]) for k in ws.PARAMS["m6"]}
    need = B["min_reps"]
    if "--pysr" in sys.argv:
        CACHE.mkdir(parents=True, exist_ok=True)
        return part_a_pysr(par)

    # ------------------------------------------------------------------------------------------- A (seconds)
    t0 = time.time()
    a_draws, a_cells, software, a_methods = part_a(par)
    a_sec = time.time() - t0
    print(f"part A done in {a_sec:.0f} s; {software}", flush=True)
    save_output(a_draws.round(4), P.STEP3B3C_A_DRAWS, "product", "step 3b-3c part A: what each jump-2 method returned on "
                "every noise draw of the true h(W)", index=False)
    save_output(a_cells, P.STEP3B3C_A, "product", "step 3b-3c part A: jump 2 alone, draws with exactly the true terms and "
                "coefficient ratios per truth, noise level and method", index=False)

    # ------------------------------------------------------------------------------------------- jobs, timing, plan
    def jobs_for(reps, plan):
        jobs = [{"part": "BC", "truth": t, "rep": r, "step": t in plan["b_truths"]} for t in C["truths"] for r in reps]
        jobs += [{"part": "D", "truth": t, "rep": r, "bootstraps": plan["bootstraps"]} for t in D_["truths"]
                 for r in D_["datasets"] if r in reps]
        return sorted(jobs, key=lambda j: (j["rep"], j["part"] != "D"))

    def run(jobs, label):
        t_run, out = time.time(), []
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(run_any, j): j for j in jobs}
            for k, fu in enumerate(as_completed(futs)):
                out.append((futs[fu], fu.result()))
                j = futs[fu]
                print(f"   {label} {k + 1}/{len(jobs)}: {j['part']} {j['truth']} {j['rep']} ({time.time() - t_run:.0f} s)", flush=True)
        return out, time.time() - t_run

    plan_file = CACHE / "plan.json"
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
        print(f"plan from {rel(plan_file)}: cut: {plan['cuts'] or 'nothing'}")
    else:
        plan = {"b_truths": B["truths"], "bootstraps": D_["bootstraps"], "cuts": []}
        print(f"\ntiming replicate 0 ({len(jobs_for([0], plan))} jobs) ...", flush=True)
        out0, wall0 = run(jobs_for([0], plan), "replicate 0")
        sec_bc = sum(r["seconds"] for j, r in out0 if j["part"] == "BC")
        sec_b = {j["truth"]: r["seconds_b"] for j, r in out0 if j["part"] == "BC"}
        sec_d = sum(r["seconds"] for j, r in out0 if j["part"] == "D")
        per_s = wall0 / (sec_bc + sec_d)                   # wall seconds per process-second in the pool
        proj = per_s * (R * sec_bc + len(D_["datasets"]) * sec_d)
        plan.update({"wall_rep0_s": round(wall0, 1), "process_s_rep0": {"BC": sec_bc, "D": sec_d}, "wall_s": round(wall0, 1)})
        print(f"   replicate 0: {wall0 / 60:.1f} min wall ({sec_bc / 60:.1f} process-min B + C, {sec_d / 60:.1f} D); "
              f"projected total {proj / 3600:.2f} h on {workers} processes (budget {S['time_budget_h']} h)", flush=True)
        if proj > S["time_budget_h"] * 3600:
            plan["cuts"].append("E: nothing to cut (its F4 fits are part of C)")
            plan["bootstraps"] = 12
            plan["cuts"].append("D: bootstraps 20 -> 12")
            proj = per_s * (R * sec_bc + len(D_["datasets"]) * sec_d * 0.6)
            print(f"   cut: D bootstraps 20 -> 12 -> projected {proj / 3600:.2f} h", flush=True)
        if proj > S["time_budget_h"] * 3600:
            plan["b_truths"] = [t for t in B["truths"] if t not in ("m3", "m5")]
            plan["cuts"].append("B: truths m3, m5")
            proj -= per_s * (R - 1) * (sec_b["m3"] + sec_b["m5"])
            print(f"   cut: B truths m3, m5 -> projected {proj / 3600:.2f} h", flush=True)
        plan["projected_s"] = round(proj)
        plan_file.write_text(json.dumps(plan))
    jobs = jobs_for(range(R), plan)
    print(f"\n{len(jobs)} jobs on {workers} processes", flush=True)
    out, wall = run(jobs, "jobs")
    if wall > 30:
        plan["wall_s"] = round(plan["wall_s"] + wall, 1)
        plan_file.write_text(json.dumps(plan))
    bc = {(r["truth"], r["rep"]): r for j, r in out if j["part"] == "BC"}
    dd = {(r["truth"], r["rep"]): r for j, r in out if j["part"] == "D"}
    r3 = {k: json.loads((CACHE_3B / f"job_main_{k[0]}_{k[1]}.json").read_text()) for k in bc}
    dec = [{"part": "A", "truth": t, "verdict": verdict(list(g.ok)),
            "basis": "; ".join(f"{r.noise:.0%}: {r.exact}/{r.draws} ({r.coef_ratios})" for r in g.itertuples())}
           for t, g in a_cells[a_cells.method == a_methods[0]].groupby("truth", sort=False)]

    # ------------------------------------------------------------------------------------------- B
    b_rows, b_sum = [], []
    for t in plan["b_truths"]:
        chosen = []
        for r in range(R):
            res, q = bc[(t, r)], r3[(t, r)]
            sc = wc.StepCurve(q["knots"], B["wc_range_db"], B["width_logw"])
            stp = res["step"]
            b_rows.append({"truth": t, "rep": r, "bic_no_step": round(q["bic"], 1), "bic_step": round(stp["bic"], 1),
                           "bic_gain": round(q["bic"] - stp["bic"], 1), "step_kept": res["step_kept"],
                           "wc_db": round(stp["wc_db"], 3), "height_ln": round(stp["height"], 2),
                           "k_low": round(float(sc.k(stp["x"][:-1], [B["wc_range_db"][0] / 2])[0]), 5),
                           "best_wc_start": stp["from"], "seconds": res["seconds_b"]})
            chosen.append((sc.k if res["step_kept"] else wc.Curve(q["knots"]).k,
                           (stp["x"] if res["step_kept"] else q["x"])[:-1]))
        g = pd.DataFrame([x for x in b_rows if x["truth"] == t])
        n_kept = int(g.step_kept.sum())
        row = {"truth": t, "step_kept": f"{n_kept}/{R}", "median_bic_gain": round(float(g.bic_gain.median()), 1)}
        for lo, name in ((W_LO, "above_0.35"), (B["step_range_lo_db"], "above_0.1")):
            w = np.geomspace(lo, min(r3[(t, r)]["w_hi"] for r in range(R)), FIT["grid_points"])
            Kc = np.array([f(x, w) for f, x in chosen])
            kt = s3.true_k(t, par, w)
            row[f"coverage_{name}"] = round(float(np.mean((kt >= Kc.min(axis=0)) & (kt <= Kc.max(axis=0)))), 3)
            row[f"k_ratio_{name}"] = round(float(np.exp(np.median(np.log(Kc / kt)))), 3)
        if t == "m6":
            wcs = g[g.step_kept].wc_db
            row["median_wc_kept"] = round(float(wcs.median()), 3) if len(wcs) else np.nan
            row["wc_within_tol"] = int((abs(wcs / par["m6"]["A_c"] - 1) <= B["wc_tol"]).sum())
            ok = bool(n_kept >= need and len(wcs) and abs(wcs.median() / par["m6"]["A_c"] - 1) <= B["wc_tol"])
            basis = (f"step kept in {n_kept} of {R}; median W_c {row['median_wc_kept']} dB (true {par['m6']['A_c']:.2f}; "
                     f"within +-{B['wc_tol']:.0%} in {row['wc_within_tol']})")
        else:
            ok = bool(R - n_kept >= need and row["coverage_above_0.35"] >= B["band_coverage"])
            basis = f"no step in {R - n_kept} of {R}; coverage above 0.35 dB {row['coverage_above_0.35']:.0%}"
        row["ok"] = ok
        b_sum.append(row)
        dec.append({"part": "B", "truth": t, "verdict": "PASS" if ok else "FAIL", "basis": basis})
    b_tab, b_sum = pd.DataFrame(b_rows), pd.DataFrame(b_sum)
    save_output(b_tab, P.STEP3B3C_B, "product", "step 3b-3c part B: per dataset, BIC of the free curve with and without "
                "the step, fitted W_c and height, step kept or not", index=False)

    # ------------------------------------------------------------------------------------------- C
    c_rows, win_rows = [], []
    for t in C["truths"]:
        for r in range(R):
            res, q = bc[(t, r)], r3[(t, r)]
            n = np.log(res["n_independent"])
            cands = {k: {**v, "family": FAMILY.get(k) or formula_family(v["terms"])} for k, v in res["cand"].items()}
            best = max(v["held"] for v in cands.values())
            fc_held = q["cv"][str(q["lam_star"])]
            for k, v in cands.items():
                c_rows.append({"truth": t, "rep": r, "candidate": k, "what": v["what"], "family": v["family"],
                               "parameters": v["k"], "held_out_ll": round(v["held"], 2),
                               "held_out_minus_best": round(v["held"] - best, 2), "ll": round(v["ll"], 2),
                               "bic": round(-2 * v["ll"] + v["k"] * n, 1)})
            c_rows.append({"truth": t, "rep": r, "candidate": "free curve (reference)", "what": "5 knots, not a candidate",
                           "family": "-", "parameters": len(q["knots"]) + 1, "held_out_ll": round(fc_held, 2),
                           "held_out_minus_best": round(fc_held - best, 2), "ll": round(q["ll"], 2), "bic": round(q["bic"], 1)})
            order = sorted(cands, key=lambda k: -cands[k]["held"])
            other = next((k for k in order if cands[k]["family"] != cands[order[0]]["family"]), None)
            near = [k for k in order if best - cands[k]["held"] <= 2]
            bicw = min(cands, key=lambda k: -2 * cands[k]["ll"] + cands[k]["k"] * n)
            win_rows.append({"truth": t, "rep": r, "winner": order[0], "winner_what": cands[order[0]]["what"],
                             "family": cands[order[0]]["family"],
                             "in_true_family": cands[order[0]]["family"] in C["families"][t],
                             "margin_to_other_family": round(best - cands[other]["held"], 2) if other else np.nan,
                             "best_other_family": f"{other} ({cands[other]['family']})" if other else "",
                             "simplest_within_2": min(near, key=lambda k: cands[k]["k"]),
                             "bic_winner": bicw, "bic_winner_family": cands[bicw]["family"],
                             "free_curve_minus_winner": round(fc_held - best, 2), "step_kept": res["step_kept"],
                             "j2_bic_pick": " + ".join(res["j2"]["bic_pick"])})
    c_tab, win = pd.DataFrame(c_rows), pd.DataFrame(win_rows)
    save_output(c_tab, P.STEP3B3C_C, "product", "step 3b-3c part C (jump 3): every candidate of every dataset, held-out "
                "log-likelihood (2-fold over event draws) and BIC", index=False)
    save_output(win, P.STEP3B3C_C_WIN, "product", "step 3b-3c part C (jump 3): winner per dataset, its family, margin",
                index=False)
    # post hoc, added after the first complete run (no refit, no verdict uses it): the same held-out scores with the
    # MENU LAWS ONLY, because the step law won for every truth (part B) and the jump-2 formulas were then read from
    # the curve with the step.
    mo = c_tab[c_tab.candidate.isin(C["menu"])]
    menu_only = []
    for (t, r), g in mo.groupby(["truth", "rep"], sort=False):
        g = g.sort_values("held_out_ll", ascending=False)
        other = g[g.family != g.family.iloc[0]]
        menu_only.append({"truth": t, "rep": r, "held_out_winner": g.candidate.iloc[0], "family": g.family.iloc[0],
                          "in_true_family": g.family.iloc[0] in C["families"][t],
                          "margin_to_other_family": round(float(g.held_out_ll.iloc[0] - other.held_out_ll.iloc[0]), 1),
                          "bic_winner": g.loc[g.bic.idxmin(), "candidate"],
                          "bic_in_true_family": g.loc[g.bic.idxmin(), "family"] in C["families"][t]})
    menu_only = pd.DataFrame(menu_only)
    mo_sum = menu_only.groupby("truth", sort=False).agg(
        held_out_winners=("held_out_winner", lambda v: ", ".join(f"{k} x{n_}" for k, n_ in v.value_counts().items())),
        in_true_family=("in_true_family", "sum"), median_margin=("margin_to_other_family", "median"),
        bic_winners=("bic_winner", lambda v: ", ".join(f"{k} x{n_}" for k, n_ in v.value_counts().items())),
        bic_in_true_family=("bic_in_true_family", "sum")).reset_index()
    c_sum = win.groupby("truth", sort=False).agg(
        winners=("winner", lambda v: ", ".join(f"{k} x{n_}" for k, n_ in v.value_counts().items())),
        families=("family", lambda v: ", ".join(f"{k} x{n_}" for k, n_ in v.value_counts().items())),
        in_true_family=("in_true_family", "sum"),
        bic_in_true_family=("bic_winner_family", lambda v: 0), median_margin=("margin_to_other_family", "median"),
        free_curve_minus_winner=("free_curve_minus_winner", "median")).reset_index()
    c_sum["bic_in_true_family"] = [int(sum(f in C["families"][t] for f in win[win.truth == t].bic_winner_family)) for t in c_sum.truth]
    for r in c_sum.itertuples():
        dec.append({"part": "C", "truth": r.truth, "verdict": "PASS" if r.in_true_family >= C["min_reps"] else "FAIL",
                    "basis": f"held-out winner in the true family ({', '.join(C['families'][r.truth])}) in {r.in_true_family} of {R}: {r.winners}"})

    # ------------------------------------------------------------------------------------------- D
    d_rows, d_cur = [], []
    for t in D_["truths"]:
        for r in D_["datasets"]:
            q = r3[(t, r)]
            cv = wc.Curve(q["knots"])
            w = np.geomspace(W_LO, q["w_hi"], FIT["grid_points"])
            Kb = np.array([cv.k(np.asarray(x)[:-1], w) for x in dd[(t, r)]["x"][:plan["bootstraps"]]])
            lo, hi = np.percentile(Kb, D_["band"], axis=0)
            kt, kf = s3.true_k(t, par, w), cv.k(np.asarray(q["x"])[:-1], w)
            aw = wc.slope_weights(np.log(w))
            sl = np.log(Kb) @ aw
            d_rows.append({"truth": t, "dataset": r, "w_hi_db": round(q["w_hi"], 3), "bootstraps": len(Kb),
                           "coverage_band": round(float(np.mean((kt >= lo) & (kt <= hi))), 3),
                           "coverage_minmax": round(float(np.mean((kt >= Kb.min(axis=0)) & (kt <= Kb.max(axis=0)))), 3),
                           "band_half_width_pct": round(float(100 * (np.exp(np.median(np.log(hi / lo)) / 2) - 1)), 1),
                           "minmax_half_width_pct": round(float(100 * (np.exp(np.median(np.log(Kb.max(axis=0) / Kb.min(axis=0))) / 2) - 1)), 1),
                           "fit_k_ratio": round(float(np.exp(np.median(np.log(kf / kt)))), 3),
                           "slope_fit": round(float(np.log(kf) @ aw), 2), "slope_lo": round(float(np.percentile(sl, D_["band"][0])), 2),
                           "slope_hi": round(float(np.percentile(sl, D_["band"][1])), 2), "true_slope": round(float(np.log(kt) @ aw), 2),
                           "seconds": dd[(t, r)]["seconds"]})
            d_rows[-1]["ok"] = bool(d_rows[-1]["coverage_band"] >= D_["coverage"])
            for i in range(len(w)):
                d_cur.append({"truth": t, "dataset": r, "w_db": w[i], "k_fit": kf[i], "k_lo": lo[i], "k_hi": hi[i],
                              "k_min": Kb[:, i].min(), "k_max": Kb[:, i].max(), "k_true": kt[i]})
        n_ok = sum(x["ok"] for x in d_rows if x["truth"] == t)
        dec.append({"part": "D", "truth": t, "verdict": "PASS" if n_ok >= D_["min_datasets"] else "FAIL",
                    "basis": "true k inside the bootstrap band on >= 80 % in " + f"{n_ok} of {len(D_['datasets'])} datasets (" +
                             ", ".join(f"{x['coverage_band']:.0%}" for x in d_rows if x["truth"] == t) + ")"})
    d_tab, d_cur = pd.DataFrame(d_rows), pd.DataFrame(d_cur)
    save_output(d_tab, P.STEP3B3C_D, "product", "step 3b-3c part D: coverage of the true k by the event-bootstrap band of "
                "one dataset, band width, slope interval", index=False)
    save_output(d_cur.round(6), P.STEP3B3C_D_CURVES, "product", "step 3b-3c part D: the bootstrap bands on the judged grid",
                index=False)

    # ------------------------------------------------------------------------------------------- E
    e_rows = []
    for r in range(R):
        res, q = bc[("m4", r)], r3[("m4", r)]
        f4, n = res["cand"]["F4"], np.log(res["n_independent"])
        bic = lambda ll: -2 * ll + 4 * n                   # noqa: E731
        e_rows.append({"rep": r, "bic_free_curve": round(q["bic"], 1), "bic_f4_3b3b": round(q["menu_F4"]["bic"], 1),
                       "bic_f4_generic_now": round(bic(f4["ll_generic"]), 1), "bic_f4_f1sol": round(bic(f4["ll_by_start"]["f1sol"]), 1),
                       "best_start": f4["start"], "curve_minus_f4_3b3b": round(q["bic"] - q["menu_F4"]["bic"], 1),
                       "curve_minus_f4_now": round(q["bic"] - bic(f4["ll"]), 1)})
    e_tab = pd.DataFrame(e_rows)
    e_tab["warning_now"] = e_tab.curve_minus_f4_now > S3B["decision"]["m4_dbic_vs_f4"]
    save_output(e_tab, P.STEP3B3C_E, "product", "step 3b-3c part E: m4, BIC(free curve) - BIC(F4) with the F1-solution "
                "start added (information; the 3b-3b verdict stays)", index=False)

    # ------------------------------------------------------------------------------------------- verdicts
    dec = pd.DataFrame(dec)
    parts = dec.groupby("part", sort=False).verdict.agg(lambda v: verdict([x == "PASS" for x in v])).reset_index()
    overall = verdict([x == "PASS" for x in dec.verdict])
    save_output(pd.concat([dec, parts.assign(truth="all parts of " + parts.part, basis=""),
                           pd.DataFrame([{"part": "overall", "truth": "", "verdict": overall, "basis": ""}])], ignore_index=True),
                P.STEP3B3C_DECISION, "product", "step 3b-3c verdict per part and truth (rules fixed before the run)", index=False)

    # ------------------------------------------------------------------------------------------- figures
    for t in plan["b_truths"]:
        fig = fig_b(t, [(bc[(t, r)], r3[(t, r)]) for r in range(R)], par)
        save_output(fig, Path(str(P.STEP3B3C_FIG_B).format(truth=t)), "diag",
                    f"step 3b-3c part B: free curve with and without the step, truth {t}", dpi=100)
        plt.close(fig)
    for t in D_["truths"]:
        fig = fig_d(t, d_cur[d_cur.truth == t])
        save_output(fig, Path(str(P.STEP3B3C_FIG_D).format(truth=t)), "diag",
                    f"step 3b-3c part D: event-bootstrap band of one dataset vs the true k, truth {t}", dpi=100)
        plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    pv = dict(zip(parts.part, parts.verdict))
    figs_b = "\n\n".join(f"![B {t}](../figures/{Path(str(P.STEP3B3C_FIG_B).format(truth=t)).name})" for t in plan["b_truths"])
    figs_d = "\n\n".join(f"![D {t}](../figures/{Path(str(P.STEP3B3C_FIG_D).format(truth=t)).name})" for t in D_["truths"])
    a_show = a_cells.assign(noise=[f"{v:.0%}" for v in a_cells.noise])
    sec_bc = sum(r["seconds"] for r in bc.values())
    sec_d = sum(r["seconds"] for r in dd.values())
    c_show = win[["truth", "rep", "winner", "winner_what", "family", "in_true_family", "margin_to_other_family",
                  "best_other_family", "simplest_within_2", "bic_winner", "free_curve_minus_winner", "step_kept", "j2_bic_pick"]]
    held_piv = (c_tab.assign(col=c_tab.candidate).pivot_table(index=["truth", "rep"], columns="col", values="held_out_minus_best")
                .round(1).reset_index())
    report = f"""# WAA_estm step 3b-3c: the discovery method finished on fake tails (rules fixed before the run)

Written by `research/WAA_estm/step3b3c_method.py` (libraries `waa_curve.py`, `waa_pooled.py`, `waa_sindy.py`; settings
`waa.curve3b3c`). Fake data only: no real tail is fitted, the frozen tail set v1 is not changed.

**Overall: {overall}.** Per part: {', '.join(f"{k} {v}" for k, v in pv.items())}; E is information only.

## 1. What was run

- **Data of B-E**: the cached fake datasets and free-curve fits of step 3b-3b (event bootstrap seeds
  {S3B['seed_base']} + replicate, {R} replicates, 3a-fix noise, spread s = {S3B['s_true']}, first {S3B['window_min']} min,
  knots from {W_LO} dB to W_hi). Nothing of 3b-3b was refitted.
- **A** {len(a_draws)} small regressions (seconds: {a_sec:.0f} s). **B + C** {len(bc)} datasets: the step fit
  ({len(B['wc_starts_db'])} starts), three jump-2 formulas, the step law if kept, and {len(C['menu'])} menu laws from 4
  starts each, every candidate also on the two halves of the event draws. **D** {len(dd)} datasets x
  {plan['bootstraps']} event bootstraps. **E** read from the F4 fits of C.
- **Timing and cut**: replicate 0 of every part was run first: {plan['wall_rep0_s'] / 60:.1f} min wall on {workers}
  processes; projected total {plan['projected_s'] / 3600:.2f} h (budget {S['time_budget_h']} h). Cut:
  {'; '.join(plan['cuts']) if plan['cuts'] else 'nothing'}.
- **Run time**: {plan['wall_s'] / 60:.0f} min wall ({sec_bc / 3600:.1f} process-hours B + C, {sec_d / 3600:.1f} D; cached in
  `.cache/fits3b3c/`). This run: {time.time() - clock:.0f} s.
- **Software**: {software}. {SOFTWARE}

## 2. Part A: jump 2 alone (`{rel(P.STEP3B3C_A)}`) - {pv['A']}

Input: the TRUE h(W) = W k(W) on {A['grid_points']} log-spaced points in [{A['w_range_db'][0]}, {A['w_range_db'][1]}] dB,
times (1 + noise), {A['draws']} draws per noise level. Library: {', '.join(TERMS_A)} (sw x = max(W - x dB, 0)). All
regressions on relative residuals (every point divided by h).
- **best-subset + BIC**: every combination of up to {A['max_terms']} terms; BIC = n log(RSS / n) + terms x log n.
- **E-SINDy**: {A['esindy']['bootstraps']} bootstraps of the points; a term is kept if included in >=
  {A['esindy']['inclusion']:.0%} of them; coefficient = median over those. Sparse step: the thresholded backward
  elimination of `waa_sindy.py`, or pysindy's STLSQ (both with relative threshold {A['esindy']['threshold']}).
Rule (judged on best-subset, the main tool; the others are cross-checks): `exact` = exactly the true terms, in >=
{A['min_draws']} of {A['draws']} draws, and the median coefficient within +-{A['coef_tol']:.0%}. m6: W plus one or both of
{', '.join(A['m6_switch_terms'])} (the true A_c = {par['m6']['A_c']:.2f} dB lies between them; their sum is compared with k_run).

{md(a_show)}

## 3. Part B: jump 1 with an optional step (`{rel(P.STEP3B3C_B)}`) - {pv['B']}

k(W) = k_low (a constant) below W_c and the 3b-3b spline above it, joined by a logistic of width {B['width_logw']} in
ln W; W_c free in [{B['wc_range_db'][0]}, {B['wc_range_db'][1]}] dB, step height >= 0 (ln units). Fitted from the 3b-3b
solution of the dataset (same knots, same penalty) with W_c starting at {B['wc_starts_db']} dB ({B['maxfev']['start']}
evaluations each, the best continued for {B['maxfev']['finish']}). The step costs 2 parameters; it is kept only if
BIC(no step) - BIC(step) > {B['bic_gain']:g}. Rule: m6 step kept in >= {need} of {R} and the median W_c of the kept steps
within +-{B['wc_tol']:.0%} of {par['m6']['A_c']:.2f} dB; m1, m3, m5 no step in >= {need} of {R} and the true k inside the band
of the chosen curves on >= {B['band_coverage']:.0%} of [{W_LO} dB, W_hi].

{md(b_sum)}

{md(b_tab)}

{figs_b}

## 4. Part C: jump 3, refit and choose (`{rel(P.STEP3B3C_C_WIN)}`, `{rel(P.STEP3B3C_C)}`) - {pv['C']}

Candidates per dataset: **J2-1, J2-2, J2-3** = the best subset of 1, 2 and 3 terms on the chosen free curve of that
dataset (no step: library {', '.join(TERMS_HI)} on [{W_LO} dB, W_hi]; step kept: the library of part A on
[{B['step_range_lo_db']} dB, W_hi]), refitted as a law h(W) of the pooled model; **step** = the step law of B if kept;
**F1, F3, F4, F5, F6** = the menu laws, best of 4 starts (3 generic + the F1 solution; F1: 4 generic). Every candidate
is fitted on one half of the event draws and scored by the log-likelihood of the other half, both ways (the halves
of 3b-3b); winner = highest sum. Families: F1 / "W" = exponential; F3 / F5 / "W^2" / "W + W^2" = quadratic; F6 / step /
"W + switch at <= 0.3 dB" = switch; F4 = two layers; other formulas = other. Rule: winner in the true family in >=
{C['min_reps']} of {R}.

{md(c_sum)}

`margin_to_other_family` = held-out log-likelihood of the winner minus the best candidate of another family;
`simplest_within_2` = the candidate with the fewest parameters within 2 of the winner; `bic_winner` = lowest BIC on
the whole dataset; `free_curve_minus_winner` = the free curve (reference, not a candidate) minus the winner;
`j2_bic_pick` = what BIC alone would pick on the curve.

{md(c_show)}

Held-out log-likelihood of every candidate minus the best one of its dataset:

{md(held_piv)}

**Post hoc, added after the first complete run (no refit; no verdict uses it): the menu laws only.** The step law
wins for every truth (see part B), and the jump-2 formulas were read from the curve with the step. The same held-out
scores, restricted to F1, F3, F4, F5, F6:

{md(mo_sum)}

## 5. Part D: the band for real data (`{rel(P.STEP3B3C_D)}`) - {pv['D']}

One dataset = replicate {', '.join(map(str, D_['datasets']))} of 3b-3b. Its event draws are resampled with replacement
{plan['bootstraps']} times (same noise, same tails per event), the free curve is refitted on each (same knots and penalty,
start = the fit of the dataset, {D_['maxfev']} evaluations). Band = {D_['band'][0]}-{D_['band'][1]} % of the refits at every W.
Rule: true k inside that band on >= {D_['coverage']:.0%} of [{W_LO} dB, W_hi] in >= {D_['min_datasets']} of
{len(D_['datasets'])} datasets per truth. `band_half_width_pct` = median half-width of the band (+- %);
`slope_lo`-`slope_hi` = the same percentiles of the slope of log k vs log W.

{md(d_tab)}

{figs_d}

## 6. Part E: F4 check for the m4 warning of 3b-3b (`{rel(P.STEP3B3C_E)}`) - information only

3b-3b fitted F4 from 3 generic starts and found BIC(free curve) - BIC(F4) > 10 in 2 of 5. Here F4 also starts from
the F1 solution of the dataset. The 3b-3b verdict stays as written.

{md(e_tab)}

With the extra start: > {S3B['decision']['m4_dbic_vs_f4']:g} in {int(e_tab.warning_now.sum())} of {R}.

## 7. The Step 4 recipe that follows

{md(dec)}

{RECIPE}

## 8. In plain words: what we will be able to say about the real drying law

{PLAIN}

## 9. Unclear / possibly wrong

{UNCLEAR}
"""
    save_output(report, P.STEP3B3C_REPORT, "product", "step 3b-3c report: jump 2 alone, jump 1 with a step, jump 3, the "
                "bootstrap band and the F4 check on fake tails; the Step 4 recipe")
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {time.time() - clock:.0f} s")


# Written after the run (the tables are computed on every run; these texts are the reading of that run and quote its
# numbers: update them if a rerun changes the tables).
SOFTWARE = """`pip install pysindy pysr sympy` worked (pysindy 2.1.0, PySR 2.7.0, scikit-learn 1.9.1, sympy 1.14.0 added
  to the miniconda environment). PySR's Julia backend (Julia 1.11.9, about 1 GB under `miniconda3/julia_env`) failed
  twice because drive C: was full; after C: was cleaned it installed in 11 min (inside the 15 min allowed). PySR
  stopped once on a Windows file lock after 36 fits; the loop now retries and resumed from its cache."""
RECIPE = """**Read the table with care: the PASS cells of m6 in B and C are not evidence.** The step is kept, and the step
law wins jump 3, for EVERY truth, also the exponential one. A tool that says "switch" for everything has not found
the switch of m6.

What each part showed:
- **A (jump 2 alone).** Given the true curve with 2-5 % noise, all three working tools find the true terms; they
  differ in how often they add one. Best-subset + BIC: true terms present in 20 of 20 everywhere, but one extra
  small term in 3 to 8 of 20 draws for m1, m3, m5 (so 5 of 8 cells miss the 16-of-20 rule); m6 exact in 20 of 20 with
  W x0.96 and switch x1.02. E-SINDy with the backward elimination of `waa_sindy.py`: exact in 20 of 20 in all 8
  cells, m6 included. PySR: exact in 20 of 20 for m1 and m3, 18 of 20 for m5, coefficients within 2 %; it cannot
  write m6 (no max in its operators) and returns a different smooth formula in every draw. pysindy's STLSQ inside
  the same bagging fails at 5 % noise (keeps all 8 terms): plain thresholding cannot handle the switch columns.
- **B (step).** The step is kept in 5 of 5 for m1, m3, m5 and m6, with BIC gains of 15 to 170. For m1 the fitted
  curve is flat below W_c, jumps up by 30-60 % in a narrow spike just above W_c (0.2-0.4 dB) and comes back to the
  true level. That is not in the truth. W_c sits near 0.2 or 0.4 dB for every truth (m6: median 0.21, true 0.16).
- **C (jump 3).** The step law has the best HELD-OUT log-likelihood in 18 of the 20 single-curve datasets, by 16 to
  65 units over the true menu law for m1, m3, m5. So this is not overfitting to one half of the events: the pooled
  likelihood itself prefers a feature near 0.2-0.4 dB that the truth does not have. m4: F4 wins in 3 of 5 (rule: 4).
  Post hoc, menu laws only: BIC on the whole dataset names the true family in 4, 5, 4, 5, 5 of 5 (m1, m3, m4,
  m5, m6), as in 3b-1c; the held-out score does the same except for m1 (F6 beats F1 in 3 of 5, by about 1 unit).
- **D (band from one dataset).** The 5-95 % event-bootstrap band is +-1 to 6 % wide and holds the true k on 0 to
  78 % of the range; no dataset reaches 80 %. The fits themselves are off by -7 to +23 % in level. The band is too
  narrow by a factor of about 3. The slope is much steadier: fitted slopes are within 0.01-0.11 of the truth.
- **E (F4 check).** The F1-solution start changes BIC(F4) by at most 1.2. F4 was at its optimum in 3b-3b; the m4
  warning stays at 2 of 5.

**The Step 4 recipe that follows (only what passed, here or in 3b-3b):**
1. **Jump 1**: the free-curve pooled fit of 3b-3b, unchanged: first 60 min, 5 knots from 0.35 dB to W_hi, read only
   for W >= 0.35 dB. **No step**: the step tool is not used on real tails.
2. **Shape**: the slope of log k vs log W over [0.35 dB, W_hi]. About 0 = exponential there, about 1 = quadratic,
   0.2-0.5 = "mildly rising, law not identified by the curve" (3b-3b).
3. **Band**: none passed. Show the event-bootstrap band, but say next to it that on fake data it was about 3
   times too narrow. For the level use the spread seen on fake replicates (about +-20 %); for the slope about
   +-0.1. Do not write "the law is inside the band".
4. **Jump 2** on the curve above 0.35 dB (library W, W^2, W^3, switch at 0.5 / 0.75 dB): E-SINDy with the
   `waa_sindy.py` elimination as the main tool, PySR as cross-check, best-subset + BIC as the list of "true terms
   plus at most one extra". Not pysindy's STLSQ.
5. **Jump 3 is not used to choose the law.** Held-out likelihood rewards whatever bends the curve near
   0.2-0.4 dB. For the comparison of named laws use the menu fits F1, F3, F4, F5, F6 with BIC on the whole
   dataset and the 3b-1c rules (several starts, F3 / F5 one family, F4 vs F6 under about 10 BIC not decisive).
6. **Before any statement about a switch or about W < 0.35 dB**: repair the likelihood at low W and repeat B, C
   and D on the cached fake datasets (minutes to an hour each). See unclear point 1."""
PLAIN = """- **For W between 0.35 dB and W_hi (about 2 dB on fake data; to be computed on the real tails), first 60 min:** we can
  say how the drying rate per dB of water changes with W, as a slope with an uncertainty of about +-0.1, and
  whether the curve there is a one-term law. "Exponential in this range" (slope near 0, equation W) and
  "quadratic in this range" (slope near 1, equation W^2) would both be recovered; this was confirmed on fresh
  fake data (3b-3b) and by jump 2 alone here.
- **The level of the rate (how fast, not which shape)** is known to about +-20 %, not to the +-5 % the bootstrap
  band suggests.
- **A mild rise (slope 0.2-0.5)** can be stated as such, but not explained: a switch and two layers look the same
  above 0.35 dB.
- **For W below 0.35 dB we can say nothing yet**, and in particular nothing about a switch level. The tool built
  for that draws a step for every truth, including a pure exponential. A step found on real tails would
  therefore mean nothing.
- **Which named law is best** can be answered by the menu fits with BIC (true family in 4-5 of 5 on fake data,
  F3 / F5 as one family), not by held-out likelihood and not by the free curve against F4.
- So the honest real-data sentence will be of the form: "between 0.35 and about 2 dB, in the first hour after
  rain, the drying rate per dB of water [is constant / rises with slope x +- 0.1]; among our candidate laws the
  data prefer family Y by BIC; below 0.35 dB the method cannot yet tell"."""
UNCLEAR = """1. **Why every truth gets a step: my reading is the 0.3 dB quantization, and I did not test it.** The simulator
   quantizes RSL in 0.3 dB steps; the likelihood treats the noise as Gaussian. The false feature sits at 0.2-0.4 dB,
   it has the shape "flat, then a fast passage" (what a signal crossing one quantization step looks like), it wins
   on held-out events, and 3b-3 already saw the free curve beat the true law at the low end. All of that fits, none
   of it is proof. Other suspects: the offset prior, the start-level floor (0.05 dB), the per-tail start profile. A
   direct test: simulate m1 without quantization and refit the step (about 15 min). Real RSL is quantized the
   same way, so the problem is not confined to fake data.
2. **If point 1 is right, the menu fits are touched too.** F6 has its switch at 0.16 dB and beats F1 on m1 data
   in 3 of 5 held-out comparisons (by about 1 unit; by 12 in one). With BIC the true family still wins in 4-5
   of 5, so the damage is small there, but F6's A_c on real tails should not be read as a physical level.
3. **The step model is harsher than it had to be**: a constant below W_c and a logistic only 0.02 wide in ln W
   (my reading of "width 0.02 dB in log W"), narrower than one cell of the internal table (0.027). The narrow
   spike above W_c is partly the spline reacting to that jump. A gentler step might behave differently; not tried.
4. **m6's W_c**: the median is 0.21 dB, inside +-30 % of 0.16 dB, but only 3 of the 5 replicates are inside it, and
   m1 gives 0.21 dB in 2 of 5. Judging the median of the kept steps was my choice; judged per replicate the cell
   would fail (3 of 5, rule 4). Either way the cell carries no information (see section 7).
5. **Jump 3 candidates were my construction**: instead of one jump-2 answer I sent the best subset of 1, 2 and 3
   terms, because BIC on a smooth fitted curve always takes the largest formula (column `j2_bic_pick`). And when
   the step was kept (always), the formulas were read from the curve WITH the step on [0.1 dB, W_hi], so they
   inherit the artefact. Jump 3 on formulas read above 0.35 dB only was therefore never tested.
6. **Half-fits start from the full-data solution** (menu laws and step: 250 evaluations; formulas: 200). The
   fitted half has thus "seen" the other half through its start. It affects all candidates alike, but the
   held-out scores are slightly optimistic.
7. **The bootstrap band may be too narrow for two different reasons**, and I cannot separate them: (a) resampling
   event draws of one dataset keeps the same noise and the same link speed factors, which is what a real
   bootstrap also does; (b) each refit starts at the dataset's solution and gets 400 evaluations, so it may not
   move far enough. (b) is a weakness of my setup; a refit from the four starts would cost 4 times more.
8. **The band was judged up to each dataset's own W_hi** (2.0-3.5 dB), where the top knot has only 200
   tail-minutes; 3b-3b judged up to the smallest W_hi of the replicates. That makes D somewhat harsher. The level
   errors (x0.93 to x1.23) are medians over the whole range, though.
9. **"+-20 % for the level, +-0.1 for the slope" in the recipe are read off 9 fits**, not derived. They are a
   placeholder until a band passes.
10. **Part A is an easy test**: independent noise on the true curve. A fitted curve has smooth, correlated
    errors, which is why best-subset + BIC cannot be used on it directly. E-SINDy's 20 of 20 used one fixed
    threshold (0.05) that I chose; no other value was tried.
11. **PySR**: judged by splitting its formula into terms with sympy (terms under 2 % of h dropped, exponents
    within 0.1 of an integer rounded). The m5 misses are formulas like W^1.78 exp(0.12 W), which fit as well
    but are not "W + W^2". 30 iterations, about 5 s per fit.
12. **The post hoc "menu laws only" table** was added after seeing the results; it reuses the stored scores and
    the decision table does not use it.
13. **Cut order**: E cost nothing extra (C already fits F4 with the F1-solution start), so the first real cut
    would have been the bootstraps of D (20 -> 12, my reading of "D-replicates 5 -> 3"). Nothing was cut:
    projected 1.3 h, run 1.6 h, plus 15 min of PySR in a fourth process.
14. **Installed software stays on the machine** (about 1.2 GB on C:, mostly Julia). A temporary folder of my smoke
    test is left in `.cache/` (`tmp*`, safe to delete).
15. **Key figures**: none proposed. The step and band figures show tools that failed."""


if __name__ == "__main__":
    main()
