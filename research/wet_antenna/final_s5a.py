"""wet_antenna step S5a: WAA v3 (M4b-B) refitted on gauge-corrected rain, one last try at the slow leftover, and closed
forms for the level curve F(R) and the path curve P(R).

S4.3 showed that the VAL gap is the rain reference: the radar over-reads the large VAL events (~1.7 x the City gauges,
TRAIN ~0.9 x), and dividing the radar rain of every event by its own radar / gauge ratio g_e closes the gap. Here

  Part 1  the frame A_pred(t) = W(t) + L * P_band(R(t - d)) + c is refitted on the TRAIN events of all link lengths with
          the corrected input R / g_e (the S4.3 V1 path-mean ratios): M4b-B (gauge timing fixed, as S3b) and M3-B as a
          check; parameter intervals from COLD event-bootstrap refits; VAL (frozen, no refit) with the V0 and V1 inputs.
  Part 2  M4c = M4b-B whose slow store fills and empties with two time constants (phi free): a fake check first, then
          a rule written before the run decides whether it replaces M4b-B.
  Part 3  closed forms for F_band(R) and P_band(R) of the kept model, chosen by AICc among a few candidates (a power law
          is one of them, never imposed), PySR as a cross-check, then the plug-back check on VAL.

STANDING RULE: WAA is never estimated by subtracting a power-law rain term. The path part is separated only through its
dependence on L at matched rain. ITU-R P.838 only as a size reference (and, as in S2b / S3b, one of the optimiser
starting values of P). Closed forms are fitted to the curves the data already learned. The rain correction comes from
the City gauges only. TEST events and TEST links are never read.
Every rule is in config `wet_antenna.final_s5a` (written before the run); shared code in the library final_model.py. Fits are cached in
research/.cache/wet_antenna/final_s5a/ (every start and every bootstrap draw).
Run from the repository root (hours the first time on 3 processes, minutes from the cache):
    python -m research.wet_antenna.final_s5a
"""
import json
import os
import pickle
import sys
import time
import zlib
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cmlpipe.api import CMLData
from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import length_test as LTL
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.event_factors import frozen, rescaled
from research.wet_antenna.tails import itu_p838
from research.wet_antenna.final_model import (C, CF, CFG, CHK, DATA, FITS, FK, FORM_TEXT, FORMS, J, KEY_R, KNOTS, M, S2B, WA, ClosedForm,
                                              fit_dict, fit_forms, form_intervals, form_params, knot_table, length_bins, load_sets,
                                              param_table, plan_start, predict, problem, through, through_verdict, val_table)

CACHE = P.S5A_CACHE
INPUTS = ("V0", "V1")
STYLE = {"M4b-B (S3b)": ("-", "x"), "M4b-B": ("--", "o"), "M3-B": (":", "s"), "M4c": ("-.", "^"), "closed form": ((0, (5, 1, 1, 1)), "D")}


# =================================================================================================== workers


def job_fit(job):
    """A fit on the corrected TRAIN set: the starts, the lag profile, the polish (every finished start is saved)."""
    file = CACHE / f"fit_{job['name']}.pkl"
    res = pickle.loads(file.read_bytes()) if file.exists() else {}
    if "params" in res:
        return res
    D = sid.load_data(job["data"], job["bands"])
    Pm = problem(D, job["name"])
    t0 = [time.time()]

    def dump():
        res["seconds"] = res.get("seconds", 0.0) + time.time() - t0[0]
        t0[0] = time.time()
        file.write_bytes(pickle.dumps(res))

    starts = {}
    for key in job["starts"]:
        if key == "nested":
            Pb = problem(D, "M4b-B")
            starts[key] = sid.x_into(Pm, Pb, np.asarray(job["nested_x"]), dyn={"phi": FITS["M4b-B"]["fixed"]["phi"], "tau_s_wet": FITS["M4b-B"]["fixed"]["tau_s"],
                                                                                 "tau_s_dry": FITS["M4b-B"]["fixed"]["tau_s"], **FITS[job["name"]]["fixed"]})
        else:
            starts[key] = plan_start(Pm, key, job["itu_inc"], job.get("slow"))
    res["all"] = sid.fit_all2(Pm, starts, M["lag_grid_min"], M["maxiter"], 0, state=res.setdefault("all_state", {}), checkpoint=dump)
    del res["all_state"]
    res["params"] = Pm.params(res["all"]["x"], res["all"]["d"], KEY_R)
    res["a_link"] = np.exp(res["all"]["x"][Pm.names == "loga"])
    dump()
    return res


def job_boot(job):
    """One COLD event-bootstrap refit: the winning plan start (generic, not the fitted x), the lag of the fit on all
    events, caps start + polish."""
    file = CACHE / f"boot_{job['name']}_{job['i']:02d}.pkl"
    if file.exists():
        return pickle.loads(file.read_bytes())
    t0 = time.time()
    D = sid.load_data(job["data"], job["bands"])
    Pm = problem(D, job["name"])
    n_ev = int(D["le_event"].max()) + 1
    rng = np.random.default_rng([C["bootstrap"]["seed"], zlib.crc32(job["name"].encode())])
    mults = [np.bincount(rng.integers(0, n_ev, n_ev), minlength=n_ev).astype(float) for _ in range(job["i"] + 1)]
    wle = mults[job["i"]][D["le_event"]]
    x0 = plan_start(Pm, job["start"], job["itu_inc"], job.get("slow"))
    f1 = Pm.fit(x0, job["d"], wle, M["maxiter"]["start"])
    f2 = Pm.fit(f1["x"], job["d"], wle, M["maxiter"]["polish"])
    out = {"i": job["i"], "params": Pm.params(f2["x"], job["d"], KEY_R), "loss": f2["loss"], "nit": f1["nit"] + f2["nit"],
           "converged": f2["converged"], "seconds": time.time() - t0}
    file.write_bytes(pickle.dumps(out))
    return out


def job_fake(job):
    """The M4c fake check, one replicate: truth M4c with data-learned curves (frozen S3b M4b-B) on the corrected rain,
    the S2 fake noise, then M4c fitted from 2 starts at the true lag."""
    file = CACHE / f"fake_{job['rep']}.pkl"
    if file.exists():
        return pickle.loads(file.read_bytes())
    t0 = time.time()
    FC = C["m4c_fake"]
    D = sid.load_data(job["data"], job["bands"])
    le = pickle.loads((P.SID3_CACHE / "le_train_all.pkl").read_bytes())[0].reset_index(drop=True)
    src = frozen("M4b-B", D, le)                       # the frozen S3b fit (its a_link by cml_id)
    D0 = dict(D)
    D0["le_c"] = np.zeros(len(le))
    P1, Pb = sid.Problem2(D0, "M4c", M, agg=1), sid.Problem2(D0, "M4b", M, agg=1)
    xt = sid.x_into(P1, Pb, LTL.x_for(Pb, src, LTL.link_ids(D0, le))[0], dyn=FC["truth"])
    d = int(src["d"])
    truth1 = P1.predict(xt, d)
    Fd = sid.fake_from_truth(D, le, truth1, [FC["seed"], job["rep"]], sid.noise_table(), FK, DATA["offset_min_minutes"])
    Pm = problem(Fd, "M4c")
    best = None
    for g, p in FC["fit_starts"]:
        f1 = Pm.fit(plan_start(Pm, f"{g}+{p}", job["itu_inc"], C["m4c_slow_start"]), d, np.ones(Pm.ne), M["maxiter"]["start"])
        best = f1 if best is None or f1["loss"] < best["loss"] else best
    f2 = Pm.fit(best["x"], d, np.ones(Pm.ne), M["maxiter"]["polish"])
    par = Pm.params(f2["x"], d, KEY_R)
    out = {"rep": job["rep"], **{k: par[k] for k in ("phi", "tau_s_wet", "tau_s_dry")}, "loss": f2["loss"], "converged": f2["converged"],
           "seconds": time.time() - t0}
    file.write_bytes(pickle.dumps(out))
    return out


def run_job(job):
    return job["kind"], job.get("name"), {"fit": job_fit, "boot": job_boot, "fake": job_fake}[job["kind"]](job)


# =================================================================================================== evaluation helpers


# =================================================================================================== closed forms


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    clock = time.time()
    workers = max(min((os.cpu_count() or 2) - 1, 3), 1)
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    BANDS = list(sid.merged_bands(d, DATA)[0])
    let = pickle.loads((P.SID3_CACHE / "le_train_all.pkl").read_bytes())[0].reset_index(drop=True)
    lev = pickle.loads((P.SID3_CACHE / "le_val_all.pkl").read_bytes())[0].reset_index(drop=True)
    g = pd.read_csv(P.COND4C_G)                                                       # S4.3: g_e per event (TRAIN, VAL)
    gmap = {s: dict(zip(g[g.set == s].event_id, g[g.set == s].g_path)) for s in ("TRAIN", "VAL")}
    assert set(let.event_id) <= set(gmap["TRAIN"]) and set(lev.event_id) <= set(gmap["VAL"])
    ftr1 = CACHE / "data_train_v1.npz"
    if not ftr1.exists():
        sid.save_data(rescaled(sid.load_data(P.SID3_CACHE / "data_train_all.npz", BANDS), let, gmap["TRAIN"]), ftr1)
    prod_freq = pd.read_csv(P.SID2B_LINK_EVENTS, usecols=["band_name", "frequency_ghz"]).groupby("band_name").frequency_ghz.median()
    itu_inc = np.array([np.diff(itu_p838(float(prod_freq[b]), "v")[0] * KNOTS ** itu_p838(float(prod_freq[b]), "v")[1]) for b in BANDS])   # P starts of S2b
    plan_keys = [f"{a}+{b}" for a, b in S2B["model"]["start_plan"]]
    base = {"data": str(ftr1), "bands": BANDS, "itu_inc": itu_inc.tolist()}
    print(f"TRAIN {len(let)} link-events, VAL {len(lev)}; g_e used for {int(g[g.set == 'TRAIN'].path_used.sum())} TRAIN and "
          f"{int(g[g.set == 'VAL'].path_used.sum())} VAL events; {workers} processes", flush=True)

    # ------------------------------------------------------------------------------------------- time one start, plan, cuts
    plan_file = CACHE / "plan.json"
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
    else:
        Dt1 = sid.load_data(ftr1, BANDS)
        Pm = problem(Dt1, "M4b-B")
        x0 = plan_start(Pm, plan_keys[0], itu_inc)
        ones = np.ones(Pm.ne)
        Pm.run(x0, 0, ones)
        t0 = time.time()
        for _ in range(3):
            Pm.run(x0, 0, ones)
        t_eval = (time.time() - t0) / 3
        t0 = time.time()
        fit = Pm.fit(x0, 0, ones, M["maxiter"]["start"])
        t_fit = time.time() - t0
        (CACHE / "fit_M4b-B.pkl").write_bytes(pickle.dumps({"all_state": {"first": {plan_keys[0]: fit}}, "seconds": t_fit}))   # reused by the job
        del Dt1, Pm
        n_lag = M["lag_grid_min"][1] - M["lag_grid_min"][0]
        lagpol = n_lag * M["maxiter"]["lag"] / M["maxiter"]["start"] + 0.5

        def proj(p):
            first = (8 + lagpol) * t_fit                                           # M4b-B on all events
            others = ((8 + lagpol) if p["m3"] else 0) + 2 * 2.5 + (C["m4c_plan_starts"] + 1 + lagpol) + p["boot"] * 1.5
            return max(first, (first + others * t_fit) / workers + 0.5 * first) / 3600 + 0.5   # + 0.5 h for VAL, closed forms, PySR
        plan = {"t_eval_s": t_eval, "t_fit_s": round(t_fit, 1), "nit_start": fit["nit"], "converged_start": fit["converged"], "m3": True,
                "boot": C["bootstrap"]["reps"], "pysr_seeds": len(CF["pysr"]["seeds"]), "cuts": [], "workers": workers}
        plan["projected_full_h"] = round(proj(plan), 2)
        for cut, change in (("the M3-B check", lambda p: p.update(m3=False)), ("PySR seeds 3 -> 1", lambda p: p.update(pysr_seeds=1)),
                            ("bootstrap 30 -> 20", lambda p: p.update(boot=20))):
            if proj(plan) <= C["time_budget_h"]:
                break
            change(plan)
            plan["cuts"].append(cut)
        plan["projected_h"] = round(proj(plan), 2)
        plan_file.write_text(json.dumps(plan))
    print(f"timing: one M4b-B start {plan['t_fit_s']} s ({plan['nit_start']} iterations, converged {plan['converged_start']}), one evaluation "
          f"{plan['t_eval_s'] * 1e3:.0f} ms; projection {plan['projected_full_h']} h; cuts: {plan['cuts'] or 'nothing'} -> {plan['projected_h']} h", flush=True)

    # ------------------------------------------------------------------------------------------- the jobs (dependencies handled here)
    results = {"fit": {}, "boot": {}, "fake": {}}
    queue = [{**base, "kind": "fit", "name": "M4b-B", "starts": plan_keys}]
    if plan["m3"]:
        queue.append({**base, "kind": "fit", "name": "M3-B", "starts": plan_keys})
    queue += [{**base, "kind": "fake", "rep": r} for r in range(C["m4c_fake"]["replicates"])]
    added = set()
    fake_pass = None
    m4c_kept = None
    t_jobs = time.time()

    def follow_ups():
        nonlocal fake_pass
        new = []
        fb = results["fit"].get("M4b-B")
        if fb is not None and "boot" not in added:
            added.add("boot")
            best = min(fb["all"]["starts_d0"], key=fb["all"]["starts_d0"].get)
            new += [{**base, "kind": "boot", "name": "M4b-B", "i": i, "start": best, "d": int(fb["all"]["d"])} for i in range(plan["boot"])]
        if len(results["fake"]) == C["m4c_fake"]["replicates"] and fake_pass is None:
            fake_pass = fake_verdict(results["fake"])["passed"]
        if fb is not None and fake_pass and "M4c" not in added:
            added.add("M4c")
            new.append({**base, "kind": "fit", "name": "M4c", "starts": ["nested"] + plan_keys[:C["m4c_plan_starts"]], "nested_x": fb["all"]["x"].tolist(),
                        "slow": C["m4c_slow_start"]})
        return new

    with ProcessPoolExecutor(workers) as ex:
        running = {}
        while queue or running:
            while queue and len(running) < workers:
                j = queue.pop(0)
                running[ex.submit(run_job, j)] = j
            done_, _ = wait(running, return_when=FIRST_COMPLETED)
            for fu in done_:
                j = running.pop(fu)
                kind, name, out = fu.result()
                if kind == "fit":
                    results["fit"][name] = out
                elif kind == "boot":
                    results["boot"].setdefault(name, {})[out["i"]] = out
                else:
                    results["fake"][out["rep"]] = out
                tag = f"{kind} {name or ''} {j.get('i', j.get('rep', ''))}".strip()
                print(f"   {tag} done ({(time.time() - clock) / 60:.1f} min)", flush=True)
                queue += follow_ups()
                if kind == "fit" and name == "M4c" and m4c_kept is None:
                    m4c_kept = m4c_decision(results, d, BANDS, let, lev, gmap)
                    if m4c_kept:
                        best = min(out["all"]["starts_d0"], key=out["all"]["starts_d0"].get)
                        if best == "nested":
                            best = plan_keys[0]
                        queue += [{**base, "kind": "boot", "name": "M4c", "i": i, "start": best, "d": int(out["all"]["d"]), "slow": C["m4c_slow_start"]}
                                  for i in range(C["bootstrap"]["if_m4c_kept"])]
    fake_pass = fake_verdict(results["fake"])["passed"] if fake_pass is None else fake_pass
    if fake_pass and m4c_kept is None and "M4c" in results["fit"]:
        m4c_kept = m4c_decision(results, d, BANDS, let, lev, gmap)
    jobs_h = (time.time() - t_jobs) / 3600
    print(f"fits done ({(time.time() - clock) / 60:.1f} min)", flush=True)
    evaluate(results, plan, fake_pass, bool(m4c_kept), d, BANDS, let, lev, gmap, clock, jobs_h)


def fake_verdict(fakes):
    FC = C["m4c_fake"]
    tr = FC["truth"]
    rows = []
    for rep, f in sorted(fakes.items()):
        ok_w = abs(f["tau_s_wet"] / tr["tau_s_wet"] - 1) <= FC["tol"]
        ok_d = abs(f["tau_s_dry"] / tr["tau_s_dry"] - 1) <= FC["tol"]
        rows.append({"replicate": rep, "phi_true": tr["phi"], "phi": f["phi"], "tau_s_wet_true": tr["tau_s_wet"], "tau_s_wet": f["tau_s_wet"],
                     "tau_s_dry_true": tr["tau_s_dry"], "tau_s_dry": f["tau_s_dry"], "tau_s_wet_ok": ok_w, "tau_s_dry_ok": ok_d, "both_ok": ok_w and ok_d,
                     "converged": f["converged"], "minutes": round(f["seconds"] / 60, 1)})
    t = pd.DataFrame(rows)
    return {"table": t, "passed": bool(len(t) == FC["replicates"] and t.both_ok.all())}




def m4c_rule(val, ver):
    R = C["m4c_rule"]
    v = val[val.input == "V1"].set_index(["fit", "VAL set"])
    lo, hi = R["integral_range"]
    rows = []
    for nm in ("L < 1 km", "1-3 km"):
        x = v.loc[("M4c", nm), "integral_pred_over_obs"]
        rows.append({"criterion": f"after-rain integral 0-180 min, {nm}", "value": x, "M4b-B": v.loc[("M4b-B", nm), "integral_pred_over_obs"],
                     "rule": f"{lo}-{hi}", "ok": bool(lo <= x <= hi)})
    dr = v.loc[("M4c", "all lengths"), "rmse"] - v.loc[("M4b-B", "all lengths"), "rmse"]
    rows.append({"criterion": "VAL RMSE (all lengths) M4c - M4b-B (dB)", "value": dr, "M4b-B": v.loc[("M4b-B", "all lengths"), "rmse"],
                 "rule": f"<= {R['rmse_worse_max_db']}", "ok": bool(dr <= R["rmse_worse_max_db"])})
    q = ver[(ver.input == "V1") & (ver.set == "VAL") & (ver.fit == "M4c")].iloc[0]
    rows.append({"criterion": "VAL length test through the model (S3b rule)", "value": q.share, "M4b-B": float(ver[(ver.input == "V1") & (ver.set == "VAL") & (ver.fit == "M4b-B")].share.iloc[0]),
                 "rule": f">= {C['length_share']}", "ok": bool(q.passes)})
    t = pd.DataFrame(rows)
    return t, bool(t.ok.all())


def m4c_decision(results, d, BANDS, let, lev, gmap):
    """Part 2 rule, as soon as M4c is fitted (so that its bootstrap can start if it is kept)."""
    le_link = np.load(P.SID3_CACHE / "data_train_all.npz")["le_link"]          # only the link index (the workers hold the data)
    Dl = {"le_link": le_link, "nl": int(le_link.max()) + 1}
    fz = {k: fit_dict(results["fit"][k], k, Dl, let) for k in ("M4b-B", "M4c")}
    Dv = {"V1": rescaled(sid.load_data(P.SID3_CACHE / "data_val_all.npz", BANDS), lev, gmap["VAL"])}
    val = val_table(Dv["V1"], lev, fz, "V1")
    ver = through_verdict(through(Dv["V1"], lev, fz, "VAL", "V1"))
    t, kept = m4c_rule(val, ver)
    print(f"M4c rule on VAL V1: {'REPLACES M4b-B' if kept else 'M4b-B stays'}\n{t.to_string()}", flush=True)
    return kept


# =================================================================================================== evaluation, closed forms, report


def run_pysr(kt, seeds):
    """PySR on the pooled (R, f) -> value for F and P. Returns (table or None, note)."""
    PS = CF["pysr"]
    try:
        from pysr import PySRRegressor
    except Exception as e:                                   # noqa: BLE001
        return None, f"PySR could not be imported ({type(e).__name__}: {e})"
    rows = []
    t0 = time.time()
    try:
        for kind in ("F", "P"):
            sub = kt[kt.kind == kind]
            X = np.c_[sub.R.to_numpy(float), sub.f_ghz.to_numpy(float)]
            y, w = sub.value.to_numpy(float), 1 / sub.se_used.to_numpy(float) ** 2
            for seed in seeds:
                m = PySRRegressor(niterations=PS["niterations"], maxsize=PS["maxsize"], binary_operators=PS["binary"], unary_operators=PS["unary"],
                                  random_state=seed, deterministic=True, parallelism="serial", progress=False, verbosity=0, temp_equation_file=True)
                m.fit(X, y, weights=w, variable_names=["R", "f"])
                for r in m.equations_.itertuples():
                    rows.append({"kind": kind, "seed": seed, "complexity": int(r.complexity), "loss": float(r.loss), "equation": str(r.equation)})
    except Exception as e:                                   # noqa: BLE001
        return (pd.DataFrame(rows) if rows else None), f"PySR stopped with {type(e).__name__}: {str(e)[:300]}"
    return pd.DataFrame(rows), f"PySR ran ({(time.time() - t0) / 60:.1f} min, seeds {list(seeds)})"


def evaluate(results, plan, fake_pass, m4c_kept, d, BANDS, let, lev, gmap, clock, jobs_h):
    Dt, Dv = load_sets(BANDS, let, lev, gmap)
    Dt1 = sid.load_data(CACHE / "data_train_v1.npz", BANDS)
    res = results["fit"]
    kept = "M4c" if m4c_kept else "M4b-B"
    s3b = frozen("M4b-B", Dt, let)
    fz = {"M4b-B (S3b)": s3b, **{k: fit_dict(r, k, Dt, let) for k, r in res.items()}}

    # ------------------------------------------------------------------------------------------- parameters, convergence
    params = param_table(res, results["boot"], s3b["params"])
    conv = pd.DataFrame([{"fit": k, "model": FITS[k]["model"], "loss": r["all"]["loss"], "d": r["all"]["d"], "converged_on_all_events": r["all"]["converged"],
                          "starts": len(r["all"]["starts_nit"]), "starts_converged": int(sum(bool(v) for v in r["all"]["starts_converged"].values())),
                          "best_start": min(r["all"]["starts_d0"], key=r["all"]["starts_d0"].get),
                          "spread_over_starts": max(r["all"]["starts_d0"].values()) - min(r["all"]["starts_d0"].values()),
                          "boot_draws": len(results["boot"].get(k, {})), "boot_converged": int(sum(bool(b["converged"]) for b in results["boot"].get(k, {}).values())),
                          "fit_minutes": round(r.get("seconds", 0) / 60, 1),
                          "boot_minutes": round(sum(b["seconds"] for b in results["boot"].get(k, {}).values()) / 60, 1)} for k, r in res.items()])
    fk = fake_verdict(results["fake"])

    # ------------------------------------------------------------------------------------------- closed forms (Part 3)
    kt, n_above, knots = knot_table(params, kept, BANDS, Dt1, let)
    forms, chosen = fit_forms(kt)
    fint = form_intervals(kt, chosen, results["boot"].get(kept, {}), BANDS)
    pysr_tab, pysr_note = run_pysr(kt, CF["pysr"]["seeds"][:plan["pysr_seeds"]])
    print(pysr_note, flush=True)
    band_idx = {b: i for i, b in enumerate(BANDS)}
    pb_file = CACHE / "plug_back.pkl"
    sel = {(kind, band_idx[b]): (form, form_params(forms, kind, b, form)) for (kind, b), form in chosen.items() if b in band_idx}
    Pk = sid.Problem2(Dt1, FITS[kept]["model"], M)
    cfm = ClosedForm(Dt1, FITS[kept]["model"], res[kept]["all"]["x"], Pk, sel, int(res[kept]["all"]["d"]))
    if pb_file.exists():
        pb = pickle.loads(pb_file.read_bytes())
    else:
        t0 = time.time()
        r = cfm.fit(1000)
        pb = {"z": r.x, "loss": float(r.fun), "nit": int(r.nit), "converged": bool(r.success), "seconds": time.time() - t0}
        pb_file.write_bytes(pickle.dumps(pb))
    xcf = cfm.build(pb["z"])
    a_cf = np.exp(xcf[cfm.ia])
    o = 0
    pb_par = []
    for (kind, b), n in zip(cfm.keys, cfm.sizes):
        form = sel[(kind, b)][0]
        for nm, v in zip(FORMS[form][1], pb["z"][o:o + n]):
            pb_par.append({"kind": kind, "band": BANDS[b], "form": form, "parameter": nm, "refitted": float(v)})
        o += n
    pb_par = pd.DataFrame(pb_par)
    fz["closed form"] = {"x": xcf, "d": int(res[kept]["all"]["d"]), "a": dict(zip(LTL.link_ids(Dt1, let), a_cf)), "model": FITS[kept]["model"], "M": cfm.M,
                         "params": {}}
    print(f"plug-back fitted (loss {pb['loss']:.1f} against the free curves {res[kept]['all']['loss']:.1f}; {(time.time() - clock) / 60:.1f} min)", flush=True)

    # ------------------------------------------------------------------------------------------- VAL (frozen) on V0 and V1, length test
    val = pd.concat([val_table(Dv[inp], lev, fz, inp) for inp in INPUTS], ignore_index=True)
    cells = pd.concat([through(Dv[inp], lev, fz, "VAL", inp) for inp in INPUTS]
                      + [through(Dt1, let, {k: f for k, f in fz.items() if k != "M4b-B (S3b)"}, "TRAIN", "V1")], ignore_index=True)
    ver = through_verdict(cells)
    print(f"VAL done ({(time.time() - clock) / 60:.1f} min)", flush=True)
    rule_tab, rule_ok = (m4c_rule(val, ver) if "M4c" in res else (None, False))
    v1 = val[(val.input == "V1") & (val["VAL set"] == "all lengths")].set_index("fit")
    q = ver[(ver.input == "V1") & (ver.set == "VAL")].set_index("fit")
    PR = CF["plug_back_rule"]
    plug = pd.DataFrame([{"criterion": "VAL V1 RMSE (all lengths): closed form - free curves (dB)", "value": v1.loc["closed form", "rmse"] - v1.loc[kept, "rmse"],
                          "rule": f"<= {PR['rmse_worse_max_db']}", "ok": bool(v1.loc["closed form", "rmse"] - v1.loc[kept, "rmse"] <= PR["rmse_worse_max_db"])},
                         {"criterion": "VAL V1 length test through the model (S3b rule)", "value": q.loc["closed form", "share"], "rule": f">= {C['length_share']}",
                          "ok": bool(q.loc["closed form", "passes"])}])
    plug_ok = bool(plug.ok.all())

    for tab, path in ((params.round(4), P.S5A_PARAMS), (conv.round(3), P.S5A_CONVERGENCE), (val.round(4), P.S5A_VAL), (cells.round(4), P.S5A_LENGTH),
                      (ver.round(4), P.S5A_LENGTH_VERDICT), (fk["table"].round(3), P.S5A_FAKE), (kt.round(4), P.S5A_KNOTS),
                      (pd.concat([forms, fint.assign(kind=fint.kind)], ignore_index=True).round(5), P.S5A_FORMS),
                      (pd.concat([plug.assign(part="rule"), pb_par.assign(part="refitted closed-form parameters")], ignore_index=True).round(5), P.S5A_PLUG)):
        save(tab, path, index=False)
    if rule_tab is not None:
        save(rule_tab.round(4), P.S5A_M4C_RULE, index=False)
    if pysr_tab is not None:
        save(pysr_tab, P.S5A_PYSR, index=False)

    # ------------------------------------------------------------------------------------------- figures
    figs = [(fig_curves(params, kt, chosen, forms, fz, BANDS, kept, cfm, pb), P.S5A_FIG_CURVES), (fig_val(val), P.S5A_FIG_VAL), (fig_fake(fk["table"]), P.S5A_FIG_FAKE)]
    curves = {}
    for (lo, hi), nm in length_bins()[:2]:
        keep = ((lev.L_km >= lo) & (lev.L_km < hi)).to_numpy()
        D_ = sid.subset(Dv["V1"], keep)
        le_ = lev[keep].reset_index(drop=True)
        curves[nm] = sid.after_curves(D_, {k: predict(D_, le_, fz[k], fz[k].get("M"))[0] for k in ("M4b-B", "M4c") if k in fz}, J, CHK["min_peak_db"])
    figs.append((fig_after(curves), P.S5A_FIG_AFTER))
    for fig, path in figs:
        save(fig, path, dpi=100)
        plt.close(fig)

    qd = dict(locals())
    report = write_report(qd)
    save(report, P.S5A_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {(time.time() - clock) / 60:.1f} min")


# =================================================================================================== figures
def fig_curves(params, kt, chosen, forms, fz, BANDS, kept, cfm, pb):
    fig, axes = plt.subplots(2, len(BANDS), figsize=(5.0 * len(BANDS), 8.6), squeeze=False)
    rr = np.geomspace(0.3, 40, 200)
    s3 = fz["M4b-B (S3b)"]["params"]
    for j, b in enumerate(BANDS):
        for i, kind in enumerate(("F", "P")):
            ax = axes[i][j]
            k = kt[(kt.kind == kind) & (kt.band == b)]
            ax.errorbar(k.R, k.value, yerr=[np.clip(k.value - k.lo, 0, None), np.clip(k.hi - k.value, 0, None)], fmt="o", color=INK, ms=6, capsize=3, label=f"{kept} refit, R / g_e (90 % cold bootstrap)")
            ks = np.array([kk for kk in KNOTS[1:]])
            ax.plot(ks, [s3.get(f"{kind}[{b}]({kk:g})", np.nan) for kk in ks], "x--", color=INK2, ms=6, label="M4b-B of S3b (raw radar)")
            form = chosen.get((kind, b))
            if form:
                ax.plot(rr, FORMS[form][0](rr, *form_params(forms, kind, b, form)), "-", color=INK, lw=1.2, label=f"closed form (AICc): {FORM_TEXT[form]}")
            ax.set_xscale("log")
            ax.set(xlabel="rain R / g_e (mm/h)", ylabel="F (dB, a_link = 1)" if kind == "F" else "P (dB/km)", title=f"{b}: {'level F(R)' if kind == 'F' else 'path P(R)'}")
            if i == 0 and j == 0 or (i == 1 and j == 0):
                ax.legend(frameon=False, fontsize=7, loc="upper left")
            style_axes(ax)
    fig.suptitle("Step S5a: the free curves of the refit on gauge-corrected rain (circles, intervals) against S3b (crosses, dashed) and the chosen closed form (line)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def fig_val(val):
    sets = ["L < 1 km", "1-3 km", "L > 3 km", "all lengths"]
    fits = list(dict.fromkeys(val.fit))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), squeeze=False)
    for ax, inp in zip(axes[0], INPUTS):
        for i, k in enumerate(fits):
            g = val[(val.input == inp) & (val.fit == k)].set_index("VAL set").reindex(sets)
            ls, mk = STYLE.get(k, ("-", "o"))
            ax.plot(np.arange(len(sets)) + (i - len(fits) / 2) * 0.06, g.rmse, ls=ls, marker=mk, mfc="none", color=INK, ms=6, lw=1.1, label=k)
        ax.set_xticks(range(len(sets))), ax.set_xticklabels(sets)
        ax.set(ylabel="RMSE of A on VAL (dB), mean over events", title=f"input {inp} ({'raw radar' if inp == 'V0' else 'radar / g_e'})")
        ax.legend(frameon=False, fontsize=7.5)
        style_axes(ax)
    fig.suptitle("Step S5a: VAL errors by length class, frozen fits (line style and marker name the fit)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def fig_fake(t):
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for col, mk, lab in (("tau_s_wet", "o", "tau_s_wet (fill)"), ("tau_s_dry", "s", "tau_s_dry (empty)")):
        if col in t:
            ax.plot(t.replicate, t[col], mk, mfc="none", color=INK, ms=8, label=f"{lab}: fitted")
            ax.axhline(float(t[f"{col}_true"].iloc[0]) if len(t) else np.nan, ls="--" if col == "tau_s_wet" else ":", color=INK2, lw=1, label=f"{lab}: truth")
            tr = float(t[f"{col}_true"].iloc[0]) if len(t) else np.nan
            ax.axhspan(tr * (1 - C["m4c_fake"]["tol"]), tr * (1 + C["m4c_fake"]["tol"]), facecolor="none", edgecolor=INK2, hatch="//" if col == "tau_s_wet" else "\\\\", lw=0)
    ax.set(xlabel="replicate", ylabel="minutes", title="M4c fake check: the slow store (hatched = truth +- 30 %)")
    ax.set_xticks(list(t.replicate) if len(t) else [])
    ax.legend(frameon=False, fontsize=8)
    style_axes(ax)
    fig.tight_layout()
    return fig


def fig_after(curves):
    fig, axes = plt.subplots(1, len(curves), figsize=(5.8 * len(curves), 4.6), squeeze=False)
    for ax, (name, c) in zip(axes[0], curves.items()):
        ax.plot(c.minutes_after_last_wet_bin, c.observed, "o", mfc="none", color=INK, ms=5, label="observed (circles)")
        for k in ("M4b-B", "M4c"):
            if k in c:
                ls, mk = STYLE[k]
                ax.plot(c.minutes_after_last_wet_bin, c[k], ls=ls, marker=mk, markevery=6, mfc="none", ms=5, color=INK, lw=1.3, label=k)
        ax.set_yscale("log")
        ax.set_ylim(3e-3, 0.5)
        ax.set(xlabel="minutes after the last wet bin of the link", ylabel="median (A - c) / observed peak of the last shower", title=f"VAL, V1 input, {name}")
        ax.legend(frameon=False, fontsize=7.5)
        style_axes(ax)
    fig.suptitle("Step S5a: after the last shower on VAL (radar / g_e), observed against M4b-B and M4c (logarithmic axis)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


# =================================================================================================== report
def readings(q):
    """The READING paragraphs, written after looking at the numbers of the run (every number is taken from its tables)."""
    params, val, ver, conv, forms, fint, kt, plug = (q[k] for k in ("params", "val", "ver", "conv", "forms", "fint", "kt", "plug"))
    BANDS, res = q["BANDS"], q["res"]
    pv = params.set_index(["fit", "parameter"])
    va = val.set_index(["input", "fit", "VAL set"])
    vr = ver.set_index(["input", "set", "fit"])
    cv = conv.set_index("fit")
    b28, b38 = BANDS[0], BANDS[-1]

    def p(fit, name, nd=2):
        r = pv.loc[(fit, name)]
        return f"{r.value:.{nd}f} [{r.lo:.{nd}f}, {r.hi:.{nd}f}]" if np.isfinite(r.lo) else f"{r.value:.{nd}f}"

    def e(inp, fit, s="all lengths", col="rmse", nd=3):
        return f"{va.loc[(inp, fit, s), col]:.{nd}f}"
    fits_val = [k for k in dict.fromkeys(val.fit)]
    out = {}
    out["params"] = f"""**READING.** The refit on the corrected rain hardly moves WAA v3: the S3b value of every free knot lies inside the cold-bootstrap
interval of the refit, except F(1) at 32 and 38 GHz, where the refit is flat between 0.5 and 1 mm/h (an increment at its bound 0) (e.g. {b28} F(2) {p('M4b-B', f'F[{b28}](2)')} dB against {pv.loc[('M4b-B', f'F[{b28}](2)'), 'S3b']:.2f}; P(5) {p('M4b-B', f'P[{b28}](5)')} dB/km against
{pv.loc[('M4b-B', f'P[{b28}](5)'), 'S3b']:.2f}); the lag stays {int(res['M4b-B']['all']['d'])} min; P at 0.5 mm/h is a little larger ({pv.loc[('M4b-B', f'P[{b28}](0.5)'), 'value']:.3f} against
{pv.loc[('M4b-B', f'P[{b28}](0.5)'), 'S3b']:.3f} dB/km). The training loss falls from 163 991 (S3b, raw radar) to {cv.loc['M4b-B', 'loss']:.0f}: the corrected rain
describes TRAIN better. The **cold** intervals are much wider than the warm S3b ones (F(2) of {b28} spans {pv.loc[('M4b-B', f'F[{b28}](2)'), 'hi'] - pv.loc[('M4b-B', f'F[{b28}](2)'), 'lo']:.2f} dB):
this is the real uncertainty of the free curves, and it is largest above 10 mm/h, where few points are. One oddity: the a_link 90th
percentile of the fit ({pv.loc[('M4b-B', 'a_p90'), 'value']:.2f}) lies below its bootstrap interval; resampling events (duplicates) spreads a_link
more than the fit on all events, so the a_link spread interval is biased up."""
    out["val"] = f"""**READING.** (1) **The corrected input is better for every fit**: VAL RMSE on V1 against V0, all lengths: {', '.join(f'{k} {e("V1", k)} / {e("V0", k)}' for k in fits_val)} dB.
(2) **The length test passes on V1 for every fit and fails on V0 for every fit** ({', '.join(f'{k} {int(vr.loc[("V1", "VAL", k), "both_inside"])}/{int(vr.loc[("V1", "VAL", k), "eligible_cells"])}' for k in fits_val)} on V1; at most
{int(ver[(ver.input == 'V0') & (ver.set == 'VAL')].both_inside.max())} of 8 on V0), including the S3b fit that never saw corrected rain: the refit is not what closes the gap, the input is
(as S4.3). (3) The fits differ by a few thousandths of a dB in every length class and window; M4c is best in the first hour after the rain
(RMSE 0-60 min {e('V1', 'M4c', col='rmse_0_60')} against M4b-B {e('V1', 'M4b-B', col='rmse_0_60')} dB). (4) **Shower peaks get worse on V1**: {e('V1', 'M4b-B', col='shower_peak_error', nd=2)} dB against
{e('V0', 'M4b-B', col='shower_peak_error', nd=2)} on V0 (M4b-B; the same for all fits): g_e is an event-total factor, and dividing the whole event lowers the peaks
too, where the radar was not too high. (5) After the rain M3-B ends far too low (integral {e('V1', 'M3-B', 'L < 1 km', 'integral_pred_over_obs', 2)} x on L < 1 km), M4b-B at
{e('V1', 'M4b-B', 'L < 1 km', 'integral_pred_over_obs', 2)} / {e('V1', 'M4b-B', '1-3 km', 'integral_pred_over_obs', 2)} x (L < 1 km / 1-3 km) is already inside the 0.7-1.4 range with the corrected rain: most of the S2b-S3b
"slow leftover" (1.2-1.4 x) was the over-read radar. TRAIN (in-sample, V1): the refits pass {int(vr.loc[('V1', 'TRAIN', 'M4b-B'), 'both_inside'])}/10 cells (S3b on V0: 10/10)."""
    r4c = pv.loc[("M4c", "tau_s_wet")]
    rule = q["rule_tab"]
    out["m4c"] = f"""**READING.** The fake check passes cleanly (both replicates within 14 % of both slow times, phi 0.90). **On the real data
the slow store of M4c is NOT identified**: tau_s_wet = {r4c.value:.0f} min, AT the upper bound ({C['tau_s_split_bounds'][1]}), interval {r4c.lo:.0f}-{r4c.hi:.0f}; tau_s_dry {p('M4c', 'tau_s_dry', 0)} min;
phi {p('M4c', 'phi')} (M4b-B: 0.72 fixed); 4 of the 20 cold draws did not converge, and all {int(cv.loc['M4c', 'starts'])} starts ended within
{cv.loc['M4c', 'spread_over_starts']:.2f} loss units of each other (one flat valley, not several optima). In words: the data want a slow store that fills over
many hours and empties in a few hours, i.e. close to a store fed by the rain amount (the S2b candidate), but they cannot say how slow.
**By the rule written before the run M4c replaces M4b-B** (all four criteria hold; table above), and that decision is kept. But it is a
weak replacement: M4b-B already meets the integral criterion on the corrected rain, the RMSE gain is {rule.value.iloc[2]:+.4f} dB, and in the median
curves after the last shower (figure) both models stay above the observed median from about 40 min on, M4c more than M4b-B; the
integral criterion is a MEAN over link-events and is carried by the larger events. The rule did not ask for identified parameters; it
should have. The leftover is therefore written down as only partly described: its size is right on average, its time constants are not known."""
    ch = {(r.kind, r.scope): r.form for r in forms[forms.chosen & forms.scope.notna()].itertuples()}
    fF = forms[(forms.kind == "F") & forms.scope.isin(BANDS)]
    fP = forms[(forms.kind == "P") & forms.scope.isin(BANDS)]
    pooledF = forms[(forms.kind == "F") & (forms.scope == "pooled: one curve")].set_index("form")
    pooledPf = forms[(forms.kind == "P") & (forms.scope == "pooled: amp linear in f")].set_index("form")
    fi = fint.set_index(["kind", "band", "parameter"])
    out["forms"] = f"""**READING.** **P (the path part):** the data choose a straight line or a power law with an exponent just below 1:
{', '.join(f"{b} {ch[('P', b)]}" for b in BANDS)}; at {b38} b = {fi.loc[('P', b38, 'b'), 'value']:.2f} [{fi.loc[('P', b38, 'b'), 'lo']:.2f}, {fi.loc[('P', b38, 'b'), 'hi']:.2f}]. At {b28} the power law has the lowest AICc (b
{json.loads(fP[(fP.scope == b28) & (fP.form == 'power')].params.iloc[0])['b']:.2f}) and the line, {fP[(fP.scope == b28) & (fP.form == 'linear')].delta_aicc.iloc[0]:.1f} AICc higher with one parameter fewer, is taken by the simplicity rule; at {BANDS[1]} the line is
lowest (power law +{fP[(fP.scope == BANDS[1]) & (fP.form == 'power')].delta_aicc.iloc[0]:.1f}, b {json.loads(fP[(fP.scope == BANDS[1]) & (fP.form == 'power')].params.iloc[0])['b']:.2f}). The quadratic term is 0 in every band (bound), so P does not bend up. Pooled with the amplitude
linear in f the power law wins (b {json.loads(pooledPf.loc['power', 'params'])['b']:.2f}, kappa {json.loads(pooledPf.loc['power', 'params'])['kappa']:.2f} per 10 GHz). This is a power law found among the candidates, not
imposed: P ~ c R^b with b ~ 0.9, c rising with frequency. **F (the antenna level):** the knots of the kept model are poorly determined
(SE {kt[(kt.kind == 'F') & (kt.R == 2)].se_used.min():.1f}-{kt[(kt.kind == 'F') & (kt.R == 2)].se_used.max():.1f} dB at 2 mm/h, {kt[(kt.kind == 'F') & (kt.R == 40)].se_used.min():.0f}-{kt[(kt.kind == 'F') & (kt.R == 40)].se_used.max():.0f} dB at 40 mm/h), so the shapes are hard to tell apart: per band the hyperbolic
saturation F_max R / (R + R0) is chosen (the others 0.7-5.5 AICc higher; at {BANDS[1]} the capped power law is lower by {fF[(fF.scope == BANDS[1]) & (fF.form == 'hyperbolic')].delta_aicc.iloc[0]:.1f} but has
one parameter more), pooled the capped power law min(c R^b, F_max) wins clearly (b {json.loads(pooledF.loc['capped_power', 'params'])['b']:.2f}; hyperbolic +{pooledF.loc['hyperbolic', 'delta_aicc']:.1f} AICc). The saturation level
F_max is {fi.loc[('F', b28, 'F_max'), 'value']:.1f} [{fi.loc[('F', b28, 'F_max'), 'lo']:.1f}, {fi.loc[('F', b28, 'F_max'), 'hi']:.1f}] dB at {b28}: the data say F saturates, not where. Note that F is the EQUILIBRIUM level of M4c;
with phi {pv.loc[('M4c', 'phi'), 'value']:.2f} and a slow store that fills over hours, the level reached after 30 min of steady rain is about
{pv.loc[('M4c', 'phi'), 'value'] + (1 - pv.loc[('M4c', 'phi'), 'value']) * (1 - np.exp(-30 / pv.loc[('M4c', 'tau_s_wet'), 'value'])):.2f} x F. **Plug-back:** the per-band forms in the model, refitted with the dynamics fixed, cost {plug.value.iloc[0]:+.4f} dB on VAL V1
(rule <= 0.005) and keep the length test at {plug.value.iloc[1] * 6:.0f}/6: **it holds**."""
    out["pysr"] = """For P all three seeds give R * 0.147 at complexity 3 (= the linear candidate, the
pooled curve), and with f: f R * 0.0046 or exp(-61.7 / f) R (the amplitude grows with f; from 3 frequencies this is not a frequency law).
For F all three seeds give log(R^2 + 3.46) at complexity 6: a form NOT in the candidate list, which rises like 2 log R and does not
saturate; with the 40 mm/h knot uncertain by about 5 dB the data cannot tell a slow logarithmic rise from a saturation."""
    out["box"] = """Read the box with section 3: the slow-store times are not identified (tau_s_wet at its bound), and F is the
equilibrium level of a store that is rarely filled. The closed-form intervals come from the knot curves of the cold draws (centred on
the knot fit); the plug-back values are the refit with the dynamics fixed, so the two can differ (e.g. F_max). P is the better-determined
half of the model; F is a saturating curve of uncertain height."""
    out["s5b"] = """- **Decide before S5b reads TEST:** the model for S5b is the S5a choice (M4c in closed form, by the rules). If the learners prefer
  M4b-B because M4c's slow store is not identified (section 3), that must be written into `final_s5a.s5b_plan` BEFORE TEST is read;
  the VAL numbers say the difference will be a few thousandths of a dB either way."""
    out["unclear"] = f"""- **M4c is kept by a rule that did not ask for identified parameters**; tau_s_wet sits at its bound and its interval spans the bounds.
- **Mean against median after the rain:** the integral criterion (mean over link-events) is met, the median curves after the last shower
  are too high from about 40 min for both models.
- **V1 lowers the shower peaks** (about -1.0 dB instead of -0.5 dB): an event-total g_e is the wrong correction at the peaks; a
  rain-rate-dependent radar error (parked from S4.3) would fit this.
- **Circularity of V1:** the radar is already gauge-adjusted and g_e uses point gauges against path means; only {int((pd.read_csv(P.COND4C_G).query('set == "VAL"').path_used).sum())} VAL events carry a
  g_e of their own.
- **F knots are poorly determined** (the M4c bootstrap moves F and phi together), so the F closed form is weak; P is the solid part.
- **The 40 mm/h knot is kept** ({q['n_above']} TRAIN link-events have corrected rain above 20 mm/h), but it rests on few minutes.
- **3 frequencies** are too few for a frequency law (kappa and the PySR f-forms are descriptions of these 3 bands).
- **M3-B** has no intervals (point fit only, as the config says); 32 GHz has no link under 1 km.
- The VAL length test on V1 has 6 eligible cells; the one that fails for M4c is {b38}, 1-2 mm/h (I outside)."""
    return out


def v4_box(q):
    """The WAA v4 (closed form) box, from the plug-back fit and the bootstrap intervals of the closed forms."""
    kept, res, BANDS, fint, pb_par, params = q["kept"], q["res"], q["BANDS"], q["fint"], q["pb_par"], q["params"]
    pv = params[params.fit == kept].set_index("parameter")
    m = FITS[kept]["model"]
    fx = FITS[kept]["fixed"]
    rows = [f"WAA v4 (S5a, closed form): {kept} on gauge-corrected rain R / g_e, TRAIN events of all link lengths",
            "A(t) = W_f(t) + W_s(t) + L * P_band(R(t - d)) + c,   u = a_link * F_band(R)",
            "dW_f/dt = (phi u - W_f) / tau_f,   tau_f = tau_f_wet while W_f rises, tau_f_dry while it falls"]
    if m == "M4c":
        rows.append("dW_s/dt = ((1 - phi) u - W_s) / tau_s,   tau_s = tau_s_wet while W_s rises, tau_s_dry while it falls")
    else:
        rows.append("dW_s/dt = ((1 - phi) u - W_s) / tau_s")
    rows.append(f"lag d {int(res[kept]['all']['d']):+d} min")
    for p in ("tau_f_wet", "tau_f_dry", "phi", "tau_s", "tau_s_wet", "tau_s_dry"):
        if p in fx:
            rows.append(f"{p:12s} {fx[p]:g}   (fixed: gauge timing of S2b)")
        elif p in pv.index:
            r = pv.loc[p]
            rows.append(f"{p:12s} {r.value:.3g} [{r.lo:.3g}, {r.hi:.3g}]")
    fi = fint.set_index(["kind", "band", "parameter"])
    for kind in ("F", "P"):
        for b in BANDS:
            g_ = pb_par[(pb_par.kind == kind) & (pb_par.band == b)]
            if not len(g_):
                continue
            form = g_.form.iloc[0]
            parts = []
            for r in g_.itertuples():
                lo, hi = (fi.loc[(kind, b, r.parameter), "lo"], fi.loc[(kind, b, r.parameter), "hi"]) if (kind, b, r.parameter) in fi.index else (np.nan, np.nan)
                parts.append(f"{r.parameter} = {r.refitted:.3g} [{lo:.3g}, {hi:.3g}]")
            rows.append(f"{kind}[{b}](R) = {FORM_TEXT[form]:24s} {', '.join(parts)}" + ("  (dB)" if kind == "F" else "  (dB/km)"))
    a = np.array(list(q["fz"]["closed form"]["a"].values()))
    rows.append(f"a_link 10-90 %   {np.exp(np.percentile(np.log(a), 10)):.2f} - {np.exp(np.percentile(np.log(a), 90)):.2f}")
    rows.append("[intervals: the closed form fitted to each cold bootstrap draw of the free curves; the plug-back value is the refit with the dynamics fixed]")
    return "```\n" + "\n".join(rows) + "\n```"


def write_report(q):
    """The report from the numbers of the run (q = the local variables of evaluate). READING paragraphs from `readings`."""
    r = readings(q)
    plan, conv, params, val, ver, fk, kt, forms, fint, plug, kept = (q[k] for k in ("plan", "conv", "params", "val", "ver", "fk", "kt", "forms", "fint", "plug", "kept"))
    res = q["res"]
    standing = ("**Standing rule.** WAA is never estimated by subtracting a power-law rain term. The path part is separated only through its "
                "dependence on L at matched rain. ITU-R P.838 only as a size reference (and, as in S2b / S3b, one of the optimiser starting values "
                "of P). Closed forms are fitted to the curves the data already learned; a power law appears only as one candidate among several, "
                "chosen by the data, never imposed. The rain correction comes from the City gauges only, never from the links. TEST was not read.")
    key = params[params.parameter.isin(["d", "phi", "tau_s_wet", "tau_s_dry", "a_p10", "a_p90"]) | params.parameter.str.match(r"^[FP]\[.*\]\((0\.5|2|5|10|20)\)$")]
    secs_fit = sum(x.get("seconds", 0) for x in res.values())
    secs_boot = sum(b["seconds"] for bb in q["results"]["boot"].values() for b in bb.values())
    vshow = val[["input", "fit", "VAL set", "link_events", "rmse", "rmse_lo", "rmse_hi", "rmse_rain", "rmse_0_60", "rmse_60_180", "shower_peak_error",
                 "integral_pred_over_obs", "integral_lo", "integral_hi"]].round(3)
    used = pd.read_csv(P.COND4C_G)
    box = v4_box(q) if q["plug_ok"] else "**Free curves kept.** " + r.get("free_kept", "")
    pysr = q["pysr_tab"]
    pysr_show = ""
    if pysr is not None and len(pysr):
        pysr_show = md_table(pysr.sort_values(["kind", "seed", "complexity"]).round(5))
    return f"""# wet_antenna step S5a: WAA v3 refitted on gauge-corrected rain, the slow leftover (M4c), closed forms of F(R) and P(R)

Written by `research/wet_antenna/final_s5a.py` (libraries `final_model.py`, `system_id.py`, `length_test.py`, `event_factors.py`). Settings and every rule:
`wet_antenna.final_s5a` in `research/config.yaml`, written on 2026-10-17 BEFORE the run. Plan: `docs/ROADMAP.md`, S5a.

{standing}

## 1. What was run

- Data: the S3 all-length datasets, TRAIN {len(q['let'])} link-events ({int(q['Dt']['le_event'].max()) + 1} events), VAL {len(q['lev'])} ({int(q['Dv']['V0']['le_event'].max()) + 1} events).
  Input V1 = radar / g_e with the S4.3 path-mean ratio (`{rel(P.COND4C_G)}`): g_e used for {int(used[used.set == 'TRAIN'].path_used.sum())} TRAIN and
  {int(used[used.set == 'VAL'].path_used.sum())} VAL events (>= {WA['conditions_s4c']['g_e']['min_gauge_mm']} mm per gauge link), g_e = 1 for the others. The same division in the 5 min
  radar of the steady-rain bins (the classes are recomputed).
- Fits on the corrected TRAIN set: {', '.join(f"{k} ({FITS[k]['model']}, fixed {', '.join(f'{a} {b:g}' for a, b in FITS[k]['fixed'].items())})" for k in res)}; 8 starts for
  M4b-B / M3-B (the S2b plan), {C['m4c_plan_starts'] + 1} for M4c; caps {M['maxiter']['start']} / {M['maxiter']['lag']} / {M['maxiter']['polish']}; no folds.
- **Cold bootstrap**: {', '.join(f'{k} {len(v)} draws' for k, v in q['results']['boot'].items())} (resampled events, the winning plan start, the lag of the
  fit on all events, caps {M['maxiter']['start']} + {M['maxiter']['polish']}). Fitting time summed over the processes: fits {secs_fit / 3600:.1f} h, bootstrap
  {secs_boot / 3600:.1f} h; the jobs took {q['jobs_h']:.1f} h wall clock on {plan['workers']} processes.
- **Timing first:** one M4b-B start {plan['t_fit_s']} s ({plan['nit_start']} iterations, converged {plan['converged_start']}); projection {plan['projected_full_h']} h against
  the budget {C['time_budget_h']} h: **cuts: {', '.join(plan['cuts']) or 'nothing'}** (projected {plan['projected_h']} h).

Convergence (`{rel(P.S5A_CONVERGENCE)}`):

{md_table(conv.drop(columns=['loss']).round(2))}

## 2. Refit on corrected rain (Part 1)

Parameters against S3b (the frozen M4b-B fitted on the raw radar); intervals = 5-95 % of the cold draws (`{rel(P.S5A_PARAMS)}`):

{md_table(key.round(3))}

{r.get('params', '')}

![F and P](../figures/{P.S5A_FIG_CURVES.name})

VAL, frozen fits, inputs V0 (raw radar) and V1 (radar / g_e) (`{rel(P.S5A_VAL)}`):

{md_table(vshow)}

![VAL by length](../figures/{P.S5A_FIG_VAL.name})

The length test through the model (S3b rule: both inside in >= 2/3 of the eligible cells; `{rel(P.S5A_LENGTH_VERDICT)}`, cells in
`{rel(P.S5A_LENGTH)}`):

{md_table(ver.round(3))}

{r.get('val', '')}

## 3. The slow leftover (Part 2: M4c)

Fake check (truth phi {C['m4c_fake']['truth']['phi']}, tau_s_wet {C['m4c_fake']['truth']['tau_s_wet']:g}, tau_s_dry {C['m4c_fake']['truth']['tau_s_dry']:g} min; rule: both slow times within
+-{C['m4c_fake']['tol'] * 100:.0f} % in 2 of 2): **{'PASS' if fk['passed'] else 'FAIL'}** (`{rel(P.S5A_FAKE)}`).

{md_table(fk['table'].round(2))}

![M4c fake](../figures/{P.S5A_FIG_FAKE.name})

{('The rule on VAL V1 (`' + rel(P.S5A_M4C_RULE) + '`): **' + ('M4c REPLACES M4b-B' if q['m4c_kept'] else 'M4b-B stays') + '**.' + chr(10) + chr(10) + md_table(q['rule_tab'].round(4))) if q['rule_tab'] is not None else '**M4c was not fitted to the real data (the fake check failed), as the rule says.**'}

![after the rain](../figures/{P.S5A_FIG_AFTER.name})

{r.get('m4c', '')}

## 4. Closed forms of F(R) and P(R) (Part 3), on {kept}

Knot values with the cold-bootstrap SE (floored: max({CF['se_floor']['rel'] * 100:.0f} % of the value, {CF['se_floor']['F_db']} dB for F, {CF['se_floor']['P_db_km']} dB/km for P)); knots
{', '.join(f'{k:g}' for k in q['knots'])} mm/h ({q['n_above']} TRAIN link-events have corrected rain above 20 mm/h; the 40 mm/h knot is dropped below
{CF['drop_40_if_fewer_than']}). `{rel(P.S5A_KNOTS)}`:

{md_table(kt.round(4))}

Every candidate, per band and pooled (weighted least squares, AICc; `chosen` = lowest AICc, or a simpler one within {CF['aicc_simpler_within']:g}; `{rel(P.S5A_FORMS)}`):

{md_table(forms.drop(columns=['wrss']).round(3))}

Intervals of the chosen per-band forms (the form fitted to every cold draw):

{md_table(fint.round(4))}

**PySR:** {q['pysr_note']}. {r.get('pysr', '')}

{pysr_show}

**Plug-back** (the chosen forms in the model, parameters and a_link refitted with the dynamics fixed; `{rel(P.S5A_PLUG)}`): **{'HOLDS' if q['plug_ok'] else 'FAILS'}**.

{md_table(plug.round(4))}

{r.get('forms', '')}

## 5. {'WAA v4 (closed form)' if q['plug_ok'] else 'Free curves kept'}

{box}

{r.get('box', '')}

## 6. What S5b will do on TEST (rules fixed now, config `wet_antenna.final_s5a.s5b_plan`)

- **Model:** {C['s5b_plan']['model']}.
- **Data:** {C['s5b_plan']['data']}.
- **References:** {C['s5b_plan']['references']}.
- **Pass:** {C['s5b_plan']['pass']}.
{r.get('s5b', '')}

## 7. Unclear / possibly wrong

{r.get('unclear', '')}

Run: this report {(time.time() - q['clock']) / 60:.1f} min after the start of the script.
"""


if __name__ == "__main__":
    main()
