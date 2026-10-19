"""wet_antenna step S2: the WAA law on links up to 3 km, the slow part, and the discovery rung M5.

Step S1 (links below 1 km, radar rain) left five things open: the path term P and the WAA level F are not well
separated by lengths of 0.14-0.97 km; the wetting speed depends on the rain input (radar 17 min, gauge 2.4 min);
the slow store of M4 is too high after the rain (water, or a slow drift of the level?); the peaks are too low; and
no model was learned freely from the data. This step answers them with the same model frame

    A_pred(t) = W(t) + L * P_band(R(t - d)) + c      (TRAIN events only; VAL is never read)

RULE: WAA is never estimated by subtracting a power-law rain term. P is a free curve, told apart from W only by the
path length L. ITU-R P.838 is used only to BUILD the fake data and as a size reference.
What the script does:
  data     the link-events of the TRAIN links below 3 km (S<3), with the nested sets S<1 and S<2; two rain inputs:
           the radar on the path (all links) and the City gauge at 1 min (links with a gauge within 1 km);
  fake     four fake checks first (the input explanation, drift against a slow store, the length sets, the discovery
           rung), with the same code;
  ladder   M1 static, M3 one store, M3D = M3 + a level drift per link-event, M4 two stores, M5 the discovery rung
           dW/dt = g(R, W) - h(W) with free surfaces; judged by leave-one-event-out prediction, per set and input;
  rules    the length rule (are the WAA numbers the same with longer links?), the gauge rule (wetting speed and
           peaks with the gauge), the slow part (water or drift?), all fixed in the config before the run;
  after    parameters with event-bootstrap intervals, figures, the product per link-event, the report.
Every rule is in config `wet_antenna.system_id_s2`. One fit is timed first and the projection decides the cuts. Fits
are cached in research/.cache/wet_antenna/system_id_s2/ (every finished stage and every few folds are saved, so an
interrupted run goes on where it stopped; delete the folder to fit again).
Run from the repository root (hours the first time on 3 processes, about 2 min from the cache):
    python -m research.wet_antenna.system_id_s2
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "1")            # one thread per process: the processes run in parallel

from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import pickle
import subprocess
import sys
import time
import zlib

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cmlpipe.api import CMLData
from research.common.config import load_config
from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.tails import spearman

CFG = load_config()
S = CFG["wet_antenna"]["system_id_s2"]
DATA, M, J, FK, LR, GR, SP, NM, CHK = S["data"], S["model"], S["judge"], S["fake"], S["length_rule"], S["gauge_rule"], S["slow_part"], S["naming"], S["checks"]
CACHE = P.CACHE / "system_id_s2"
SMOKE = bool(os.environ.get("S2_SMOKE"))                   # a code test: every 12th link-event, tiny iteration caps, its own
if SMOKE:                                                  #   cache folder, nothing written to outputs/ (S2_SMOKE=1)
    CACHE = P.CACHE / "system_id_s2_smoke"
    M["maxiter"] = {"start": 8, "lag": 3, "polish": 8, "em": 5, "fold": 3, "boot": 2}
    M["lag_grid_min"] = [-1, 2]
    J["boot_reps"], J["param_boot_reps"] = 50, 4
    FK["truths"]["Tc"]["boot_reps"] = 4

    def save(obj, path, **kw):                             # noqa: F811
        print(f"   [smoke] not saved: {path.name}")
PADM = sid.PADM
KNOTS = np.asarray(M["knots_mmh"], float)
LADDER = list(M["ladder"])
KEY_R = list(LR["key_rain_mmh"])
LINE = {"M1": ":", "M3": "--", "M3D": "-.", "M4": (0, (5, 1, 1, 1, 1, 1)), "M5": "-", "M5D": (0, (1, 1))}
MARK = {"M1": "v", "M3": "s", "M3D": "D", "M4": "^", "M5": "o", "M5D": "x"}
BAND_STYLE = [("-", "o"), ("-.", "^"), ("--", "s")]        # one line style and marker per band
SET_LABEL = {"S1": "S<1, radar", "S2": "S<2, radar", "S3": "S<3, radar", "GG": "gauge links, GAUGE input", "GR": "gauge links, radar input (same points)"}


# =================================================================================================== the data
def merged_bands(d):
    return sid.merged_bands(d, DATA)


def gauge_variants(D, has_g):
    return sid.gauge_variants(D, has_g)


def count_rows(name, D, le):
    """Counts of one dataset, all bands and per band."""
    w = D["w5"] > 0
    rows = []
    for band, g in [("all", le)] + [(b, le[le.band_name == b]) for b in D["bands"]]:
        if len(g) == 0:
            continue
        m = np.isin(D["bin_le"], g.index)
        rows.append({"set": SET_LABEL[name], "band": band, "links": g.cml_id.nunique(), "sub_links": g.sublink_id.nunique(), "link_events": len(g),
                     "events": g.event_id.nunique(), "points_5min": int(w[m].sum()), "wet_points": int((w & D["wet5"])[m].sum()),
                     "L_min_km": round(g.L_km.min(), 2), "L_median_km": round(g.L_km.median(), 2), "L_max_km": round(g.L_km.max(), 2),
                     "offset_flagged": int(g.offset_flag.sum())})
    return rows


# =================================================================================================== the fake data
def truth_curves(le, bands):
    return sid.truth_curves(le, bands, FK)


def make_fake(D, le, spec, seed, noise):
    return sid.make_fake(D, le, spec, seed, noise, FK, M["drift"]["t_ref_min_before_rain"], DATA["offset_min_minutes"])


def smoothed_input(D, sd, seed):
    """Fake check Ta: the same data, but the rain the FIT gets is the 5 min block mean of the true rain times a
    log-normal factor per 5 min bin (what a radar pixel gives: a time average with an error)."""
    rng = np.random.default_rng(seed)
    Fd = dict(D)
    blocks = D["rain"].reshape(-1, 5).mean(axis=1)          # every stored rain series has a multiple of 5 minutes
    Fd["rain"] = np.repeat(blocks * np.exp(rng.normal(0, sd, len(blocks))), 5)
    Fd["y5"], Fd["y1"] = np.where(D["w5"] > 0, D["y5"], np.nan), np.where(D["w1"] > 0, D["y1"], np.nan)
    return sid.finish(Fd, D["bands"])


# =================================================================================================== the jobs (worker processes)
def run_job(job):
    """One model on one dataset: the fit on all events, the leave-one-event-out predictions and (if asked) the
    event-bootstrap refits. Every stage, and every few folds or draws, is saved in the cache."""
    file = CACHE / f"{job['key']}.pkl"
    res = pickle.loads(file.read_bytes()) if file.exists() else {}
    D = sid.load_data(job["data"], job["bands"])
    if job["model"] == "M5n":
        Pm = sid.Problem2n(D, M, *named_terms(job["formula"]))
    else:
        Pm = sid.Problem2(D, job["model"], M, lam=job.get("lam", 0.0))
    t0 = [time.time()]

    def dump():
        res["seconds"] = res.get("seconds", 0.0) + time.time() - t0[0]
        t0[0] = time.time()
        file.write_bytes(pickle.dumps(res))

    if "all" not in res and job["model"] == "M5n":         # the named formula: one start (M5's P, s_band, a_link and the
        ones = np.ones(Pm.ne)                              #   E-SINDy coefficients), M5's lag, the polish iterations
        r = Pm.fit(np.asarray(job["x0"]), job["d"], ones, M["maxiter"]["polish"])
        res["all"] = {"d": job["d"], "x": r["x"], "loss": r["loss"], "nit": r["nit"], "s_b": [None],
                      "starts_d0": {"M5 + E-SINDy": r["loss"]}, "lag_profile": {}}
        dump()
    if "all" not in res:
        res["all"] = sid.fit_all2(Pm, M["starts"], M["lag_grid_min"], M["maxiter"], M["drift"]["em_rounds"])
        dump()
    fit = res["all"]
    if Pm.drift:
        Pm.s_b = fit["s_b"][-1]
    if "params" not in res:
        res["pred_all"], res["params"] = Pm.predict(fit["x"], fit["d"]), Pm.params(fit["x"], fit["d"], KEY_R)
        res["a_link"] = np.exp(fit["x"][Pm.names == "loga"])
        res["b"] = fit["x"][Pm.names == "b"] if Pm.drift else None
        dump()
    if job.get("loo") and "loo_pred" not in res:
        st = res.setdefault("loo_state", {"pred": np.full(len(D["y5"]), np.nan), "folds": []})
        done_ev = {f["event"] for f in st["folds"]}
        for e in np.unique(D["le_event"]):
            if int(e) in done_ev:
                continue
            fw = sid.fit_weighted2(Pm, fit, (D["le_event"] != e).astype(float), M["maxiter"]["fold"])
            held = D["le_event"][D["bin_le"]] == e
            st["pred"][held] = Pm.predict(fw["x"], fw["d"])[held]
            st["folds"].append({"event": int(e), "nit": fw["nit"], "loss": fw["loss"]})
            if len(st["folds"]) % 6 == 0:
                dump()
        res["loo_pred"], res["folds"] = st["pred"], st["folds"]
        del res["loo_state"]
        dump()
    if job.get("boot") and len(res.get("boot", [])) < job["boot"]:
        draws = res.setdefault("boot", [])
        n_ev = int(D["le_event"].max()) + 1
        rng = np.random.default_rng([J["seed"], zlib.crc32(job["key"].encode())])
        mults = [np.bincount(rng.integers(0, n_ev, n_ev), minlength=n_ev).astype(float) for _ in range(job["boot"])]
        for i in range(len(draws), job["boot"]):
            fw = sid.fit_weighted2(Pm, fit, mults[i][D["le_event"]], M["maxiter"]["boot"])
            draws.append(Pm.params(fw["x"], fw["d"], KEY_R))
            if len(draws) % 10 == 0:
                dump()
        dump()
    return job["key"]


def done(job):
    file = CACHE / f"{job['key']}.pkl"
    if not file.exists():
        return False
    res = pickle.loads(file.read_bytes())
    return "params" in res and (not job.get("loo") or "loo_pred" in res) and (not job.get("boot") or len(res.get("boot", [])) >= job["boot"])


def run_pool(jobs, workers, label):
    jobs = [j for j in jobs if not done(j)]
    if not jobs:
        return
    t0 = time.time()
    print(f"{label}: {len(jobs)} jobs", flush=True)
    with ProcessPoolExecutor(workers) as ex:
        for i, f in enumerate(as_completed([ex.submit(run_job, j) for j in jobs]), 1):
            print(f"   {label}: {f.result()} done ({i} of {len(jobs)}, {(time.time() - t0) / 60:.1f} min)", flush=True)


_RESULTS = {}


def result(key):
    if key not in _RESULTS:
        _RESULTS[key] = pickle.loads((CACHE / f"{key}.pkl").read_bytes())
    return _RESULTS[key]


def evaluations(job, n_ev):
    """Loss evaluations a job needs (the iteration caps are usually reached; 1.1 evaluations per iteration)."""
    mi, n_lag = M["maxiter"], M["lag_grid_min"][1] - M["lag_grid_min"][0]
    n = len(M["starts"]) * mi["start"] + n_lag * mi["lag"] + mi["polish"] + (M["drift"]["em_rounds"] * mi["em"] if job["model"].endswith("D") else 0)
    return 1.1 * (n + (n_ev * mi["fold"] if job.get("loo") else 0) + job.get("boot", 0) * mi["boot"])


def named_terms(formula):
    """The (name, function) pairs of the terms of a named formula {"g": [names], "h": [names]} (config `naming`)."""
    g, h = sid.naming_terms(NM)
    gd, hd = dict(g), dict(h)
    return [(n, gd[n]) for n in formula["g"]], [(n, hd[n]) for n in formula["h"]]


def formula_text(names, xi, var):
    """'0.12 R/(1+R/2) - 0.03 W' from term names and coefficients."""
    parts = [f"{c:+.4g} {n}" for n, c in zip(names, xi)]
    return (f"{var} = " + " ".join(parts).lstrip("+")) if parts else f"{var} = 0"


def run_pysr():
    """The PySR cross-check of the naming stage (config `naming.pysr`), in its own process (Julia): g(R, W) and h(W)
    on the samples of naming_samples.json, weighted by the visits. Writes naming_pysr.json."""
    from pysr import PySRRegressor
    Y = NM["pysr"]
    smp = json.loads((CACHE / "naming_samples.json").read_text())
    out = {}
    for key, variables in (("g", ["R", "W"]), ("h", ["W"])):
        X, y, w = np.asarray(smp[key]["X"], float), np.asarray(smp[key]["y"], float), np.asarray(smp[key]["w"], float)
        t0 = time.time()
        model = PySRRegressor(niterations=Y["niterations"], maxsize=Y["maxsize"], binary_operators=[o for o in Y["operators"] if o != "exp"],
                              unary_operators=["exp"], random_state=Y["seed"], deterministic=True, parallelism="serial",
                              progress=False, verbosity=0, output_directory=str(CACHE / "pysr"))
        for attempt in range(3):                           # Windows sometimes locks PySR's checkpoint file: try again
            try:
                model.fit(X, y, weights=w / w.sum() + 1e-12, variable_names=variables)
                break
            except PermissionError:
                if attempt == 2:
                    raise
                time.sleep(2)
        eq = model.equations_
        out[key] = {"best": str(model.sympy()), "seconds": round(time.time() - t0, 1),
                    "front": [{"complexity": int(c), "loss": float(l), "equation": str(e)} for c, l, e in zip(eq.complexity, eq.loss, eq.equation)]}
        print(f"   PySR {key}: {out[key]['best']} ({out[key]['seconds']} s)", flush=True)
    (CACHE / "naming_pysr.json").write_text(json.dumps(out, indent=1))


def name_m5(D3, lam_star, n_ev):
    """The naming stage (config `naming`): E-SINDy on the surfaces of M5 (S<3, radar), PySR as a cross-check, and the
    named formula refitted as the model M5n and judged against M5. Returns the tables and the report text."""
    r5 = result("real_S3_M5")
    P5 = sid.Problem2(D3, "M5", M, lam=lam_star)
    rk, wk = P5.rk, P5.wk
    vis_file = CACHE / "naming_visits.npz"
    if not vis_file.exists():                             # how often the data visit each knot (a 1 min run of the fit)
        R, W = sid.m5_states(P5, r5["all"]["x"], r5["all"]["d"])
        gv, hv = sid.knot_visits(R, W, rk, wk)
        np.savez(vis_file, gv=gv, hv=hv, w_wet_q=np.percentile(W[R > 0], [50, 90, 99]), w_all_q=np.percentile(W, [50, 90, 99]))
    v = np.load(vis_file)
    gv, hv = v["gv"], v["hv"]
    G, H = P5.unpack(r5["all"]["x"])[5:7]
    g_terms, h_terms = sid.naming_terms(NM)
    Rg, Wg = (a.ravel() for a in np.meshgrid(rk[1:], wk, indexing="ij"))
    samples = {"g": (np.column_stack([Rg, Wg]), G[1:].ravel(), gv[1:].ravel()), "h": (wk[1:, None], H[1:], hv[1:])}
    theta = {"g": np.column_stack([f(Rg, Wg) for _, f in g_terms]), "h": np.column_stack([f(wk[1:]) for _, f in h_terms])}
    rng = np.random.default_rng(NM["esindy"]["seed"])
    es, pick, term_rows = {}, {}, []
    for key, terms in (("g", g_terms), ("h", h_terms)):
        rows, p = sid.esindy(theta[key], samples[key][1], samples[key][2], NM["esindy"], rng)
        es[key], pick[key] = rows, p
        for r in rows:
            for i, (n, _) in enumerate(terms):
                if r["freq"][i] > 0:
                    term_rows.append({"surface": key, "lambda": r["lambda"], "term": n, "bag_frequency": r["freq"][i], "kept": bool(r["keep"][i]),
                                      "coefficient": r["xi"][i], "terms_in_formula": r["terms"], "rel_error": r["rel_error"], "picked": r is p})
    terms_tab = pd.DataFrame(term_rows)
    formula = {k: [n for (n, _), kp in zip(t, pick[k]["keep"]) if kp] for k, t in (("g", g_terms), ("h", h_terms))}
    coef = {k: pick[k]["xi"][pick[k]["keep"]] for k in ("g", "h")}
    text_g, text_h = formula_text(formula["g"], coef["g"], "g(R, W)"), formula_text(formula["h"], coef["h"], "h(W)")
    print(f"naming: E-SINDy {text_g};  {text_h}", flush=True)

    # PySR, its own process (the cross-check only)
    pysr_file = CACHE / "naming_pysr.json"
    if not pysr_file.exists():
        (CACHE / "naming_samples.json").write_text(json.dumps({k: {"X": X.tolist(), "y": y.tolist(), "w": w.tolist()} for k, (X, y, w) in samples.items()}))
        print("naming: PySR (separate process) ...", flush=True)
        try:
            subprocess.run([sys.executable, "-m", "research.wet_antenna.system_id_s2", "--pysr"], check=True, cwd=P.REPO, timeout=60 * 60)
        except Exception as e:                             # noqa: BLE001  (PySR / Julia missing: the cross-check is left out)
            print(f"   PySR failed: {type(e).__name__}: {e}")
    pysr = json.loads(pysr_file.read_text()) if pysr_file.exists() else {}

    # M5n: the named formula as a model, judged by the same judge against M5
    Pn = sid.Problem2n(D3, M, *named_terms(formula))
    job = {"key": "real_S3_M5n", "data": CACHE / "data_S3.npz", "model": "M5n", "bands": list(D3["bands"]), "set": "S3", "loo": True,
           "formula": formula, "d": int(r5["all"]["d"]), "x0": Pn.start_from(P5, r5["all"]["x"], coef["g"], coef["h"]).tolist()}
    run_pool([job], 1, "naming (M5n, leave-one-event-out)")
    rn = result("real_S3_M5n")
    e5, en = sid.event_errors(D3, r5["loo_pred"], J), sid.event_errors(D3, rn["loo_pred"], J)
    dif = paired(en.rmse.to_numpy(), e5.rmse.to_numpy(), n_ev)
    sh = sid.shower_table(D3, J["shower_gap_min"])
    judge_rows = []
    for m, r, e in (("M5", r5, e5), ("M5n", rn, en)):
        pk = shower_stat(*sid.shower_peaks(D3, sh, r["loo_pred"], J["shower_min_peak_db"]), n_ev)
        judge_rows.append({"model": m, "rmse": float(e.rmse.mean()), "rmse_60_180": float(e["rmse_60_180"].mean()), "shower_peak_error_median": pk[0],
                           "t_wet": r["params"]["t_wet"], "t_dry": r["params"]["t_dry"], **{f"level({x:g})": r["params"][f"Weq[{D3['bands'][0]}]({x:g})"] for x in KEY_R},
                           "free_numbers_of_g_and_h": (len(rk) - 1) * len(wk) + len(wk) - 1 if m == "M5" else len(formula["g"]) + len(formula["h"])})
    judge_tab = pd.DataFrame(judge_rows)
    counts = bool(dif[1] <= 0 <= dif[2] or dif[2] < 0)
    Gn, Hn = Pn.unpack(rn["all"]["x"])[5:7]
    G0, H0 = Pn.unpack(np.asarray(job["x0"]))[5:7]         # is the named formula a usable model? W under 2 mm/h for 2 h from 0
    a0, s0 = Pn.unpack(np.asarray(job["x0"]))[2], Pn.unpack(np.asarray(job["x0"]))[7]
    w_one = sid.surface_path(G0, H0, rk, wk, 1.0, 2.0, 0.0, 120, Pn.nsub)[-1]
    w_run = sid.surface_path(G0, H0, rk, wk, a0.max() * s0.max(), 2.0, 0.0, 120, Pn.nsub)[-1]
    stalled = rn["all"]["nit"] <= 1 or rn["all"]["loss"] > 10 * r5["all"]["loss"]
    stall_text = (f"\n\n**The refit could not start.** With the E-SINDy coefficients g rises with W (the term W) while h stops rising above "
                  f"3 dB (h = {H0[-3]:.3f}, {H0[-2]:.3f}, {H0[-1]:.3f} dB/min at W = {wk[-3]:g}, {wk[-2]:g}, {wk[-1]:g} dB). Under 2 mm/h for 2 h W reaches {w_one:.1f} dB with "
                  f"a_link s_band = 1, but {w_run:.0f} dB on the link with the largest a_link s_band of M5 ({a0.max() * s0.max():.1f}): there W runs away. The loss at "
                  f"the start is {rn['all']['loss']:.3g} (M5: {r5['all']['loss']:.0f}) and the optimiser stopped after {rn['all']['nit']} iteration(s). "
                  "The formula is not a bounded model of W.") if stalled else ""
    fitted_g = formula_text(formula["g"], rn["all"]["x"][Pn.names == "cg"], "g(R, W)")
    fitted_h = formula_text(formula["h"], rn["all"]["x"][Pn.names == "ch"], "h(W)")
    surf = pd.DataFrame([{"surface": "g", "R_mmh": r_, "W_db": w_, "M5_db_per_min": G[j, k], "visits_min": gv[j, k], "M5n_db_per_min": Gn[j, k]}
                         for j, r_ in enumerate(rk) if j > 0 for k, w_ in enumerate(wk)]
                        + [{"surface": "h", "R_mmh": np.nan, "W_db": w_, "M5_db_per_min": H[k], "visits_min": hv[k], "M5n_db_per_min": Hn[k]} for k, w_ in enumerate(wk) if k > 0])
    pysr_lines = "; ".join(f"{k}: `{v['best']}`" for k, v in pysr.items()) if pysr else "PySR did not run (see the log)."
    text = f"""E-SINDy (the backward elimination of `sindy.py`, {NM['esindy']['bags']} bags of the {len(Rg)} g-points and {len(wk) - 1} h-points,
weights = minutes the 1 min states of the fit visit each knot; a term is kept if selected in >= {NM['esindy']['keep_share']:.0%} of the bags;
lambda by the fewest terms within 10 % of the best error). The g library has {len(g_terms)} terms, the h library {len(h_terms)}
(`exp(-R/R0)` is used as `1 - exp(-R/R0)` so that g(0, W) = 0). Picked: lambda {pick['g']['lambda']:g} for g (relative weighted error
{pick['g']['rel_error']:.2f}), lambda {pick['h']['lambda']:g} for h (error {pick['h']['rel_error']:.2f}):

    {text_g}
    {text_h}

PySR cross-check (operators {', '.join(NM['pysr']['operators'])}, maxsize {NM['pysr']['maxsize']}, its 'best' formula): {pysr_lines}

M5n = this formula as a model (its {len(formula['g']) + len(formula['h'])} coefficients free; P, s_band, a_link as M5; M5's lag; one start from M5 and the
E-SINDy coefficients), refitted on all events and judged leave-one-event-out. Coefficients after the refit{' (unchanged: the refit stopped at once, see below)' if stalled else ''}:

    {fitted_g}
    {fitted_h}

{md_table(judge_tab.round(4))}

{stall_text}

RMSE(M5n) - RMSE(M5) = {dif[0]:+.4g} dB (90 % interval {dif[1]:+.4g} .. {dif[2]:+.4g}): **the named formula {'COUNTS' if counts else 'does NOT count'}**
(rule: it counts only if its paired difference to M5 holds 0{'; here it is even better' if dif[2] < 0 else ''}). Terms per lambda: `{rel(P.SID2_NAMING_TERMS)}`."""
    return {"qualifies": True, "done": True, "text": text, "surf": surf, "terms": terms_tab, "judge": judge_tab, "counts": counts, "diff": dif,
            "formula": formula, "g_text": fitted_g, "h_text": fitted_h, "esindy_g": text_g, "esindy_h": text_h, "pysr": pysr,
            "Gn": Gn, "Hn": Hn, "g_best_error": min(r["rel_error"] for r in es["g"]), "visits": (gv, hv, v["w_wet_q"], v["w_all_q"]), "rmse_M5n": float(en.rmse.mean()), "stalled": stalled, "w_run": w_run}


# =================================================================================================== the judge
def paired(a, b, n_ev):
    return sid.paired(a, b, n_ev, J)


def shower_stat(ev, err, n_ev):
    return sid.shower_stat(ev, err, n_ev, J)


def ladder(prefix, D, models, preds=None):
    """The judge on the leave-one-event-out predictions of one dataset. preds = {model: prediction} replaces the
    cached predictions (used for the drift-allowed judge). Returns (table, best, chosen, errors per model)."""
    preds = preds or {m: result(f"{prefix}_{m}")["loo_pred"] for m in models}
    n_ev = int(D["le_event"].max()) + 1
    errors = {m: sid.event_errors(D, preds[m], J) for m in models}
    tab, best, chosen = sid.judge(errors, [m for m in J["order"] + ["M5D"] if m in models], J)
    sh = sid.shower_table(D, J["shower_gap_min"])
    peaks = pd.DataFrame([dict(zip(("shower_peak_error_median", "shower_peak_error_lo", "shower_peak_error_hi", "showers"),
                                   shower_stat(*sid.shower_peaks(D, sh, preds[m], J["shower_min_peak_db"]), n_ev)), model=m) for m in models])
    tab = tab.drop(columns=["skill_vs_M0", "skill_lo", "skill_hi"]).merge(peaks, on="model")
    return tab, best, chosen, errors


def drift_ladder(prefix, D, models):
    """The drift-allowed judge (config `judge.drift_allowed`): every model keeps one linear drift per held-out
    link-event, with the prior sd s_b of M3D on this dataset. Returns (table, errors, slopes of M3D's fit on all events)."""
    r = result(f"{prefix}_M3D")
    tau_h = sid.Problem2(D, "M1", M).drift_times()
    preds = {m: sid.drift_allowed(D, result(f"{prefix}_{m}")["loo_pred"], tau_h, r["all"]["s_b"][-1], M["huber_delta_db"])[0] for m in models}
    tab, _, _, errors = ladder(prefix, D, models, preds)
    return tab, errors, r["b"]


def slow_verdict(prefix, D, models=("M3D", "M4")):
    """The rule of config `slow_part` on one dataset. Returns a row of numbers with the verdict."""
    n_ev = int(D["le_event"].max()) + 1
    _, errors, b = drift_ladder(prefix, D, list(models))
    dif = paired(errors["M3D"].rmse.to_numpy(), errors["M4"].rmse.to_numpy(), n_ev)
    plain = paired(sid.event_errors(D, result(f"{prefix}_M3D")["loo_pred"], J).rmse.to_numpy(),
                   sid.event_errors(D, result(f"{prefix}_M4")["loo_pred"], J).rmse.to_numpy(), n_ev)
    s_ev, n_evs = np.bincount(D["le_event"], weights=b, minlength=n_ev), np.bincount(D["le_event"], minlength=n_ev)
    draws = np.random.default_rng(J["seed"]).integers(0, n_ev, (J["boot_reps"], n_ev))
    b_lo, b_hi = np.percentile(s_ev[draws].sum(axis=1) / n_evs[draws].sum(axis=1), J["interval"])
    depth = np.bincount(D["bin_le"], weights=np.nan_to_num(D["r5"]) * D["wet5"] / 12.0, minlength=len(D["le_n"]))
    water = dif[1] > 0
    verdict = "WATER-LIKE" if water else ("DRIFT-LIKE" if b_lo <= 0 <= b_hi else "SLOW RISE, FORM UNDECIDED")
    return {"M3D_minus_M4_from_R_alone": plain[0], "from_R_lo": plain[1], "from_R_hi": plain[2],
            "M3D_minus_M4_drift_allowed": dif[0], "drift_allowed_lo": dif[1], "drift_allowed_hi": dif[2],
            "s_b_db_h": result(f"{prefix}_M3D")["all"]["s_b"][-1], "mean_slope_db_h": float(b.mean()), "slope_lo": float(b_lo), "slope_hi": float(b_hi),
            "spearman_slope_rain_depth": float(spearman(pd.Series(b), pd.Series(depth))), "tau_s_M4_min": result(f"{prefix}_M4")["params"]["tau_s"], "verdict": verdict}


def key_names(model, bands):
    """The key WAA numbers of the length rule for a model: (time constants, level names per band)."""
    taus = {"M3": ["tau_wet", "tau_dry"], "M3D": ["tau_wet", "tau_dry"], "M4": ["tau_f_wet", "tau_f_dry"], "M5": ["t_wet", "t_dry"], "M5D": ["t_wet", "t_dry"], "M1": []}[model]
    level = "Weq" if model.startswith("M5") else "F"
    return taus, [f"{level}[{b}]({r:g})" for b in bands for r in KEY_R]


def boot_interval(key, name):
    v = np.array([dr.get(name, np.nan) for dr in result(key).get("boot", [])], float)
    return tuple(np.nanpercentile(v, J["interval"])) if np.isfinite(v).any() else (np.nan, np.nan)   # a band absent from the set: no interval


def after_curves(D, preds):
    return sid.after_curves(D, preds, J, CHK["min_peak_db"])


def residuals_by_length(D, le, pred):
    """Mean held-out residual (predicted - observed) during rain per length bin, with event-bootstrap intervals."""
    n_ev = int(D["le_event"].max()) + 1
    draws = np.random.default_rng(J["seed"]).integers(0, n_ev, (J["boot_reps"], n_ev))
    edges = LR["residual_bins_km"]
    L_bin, ev_bin = le.L_km.to_numpy()[D["bin_le"]], D["le_event"][D["bin_le"]]
    res = pred - D["y5"]
    rows = []
    for label, extra in (("all wet points", True), (f"wet points with R >= {LR['residual_min_rain_mmh']:g} mm/h", np.nan_to_num(D["r5"]) >= LR["residual_min_rain_mmh"])):
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (D["w5"] > 0) & D["wet5"] & (L_bin >= lo) & (L_bin < hi) & extra
            if m.sum() == 0:
                continue
            s, n = np.bincount(ev_bin[m], weights=res[m], minlength=n_ev), np.bincount(ev_bin[m], minlength=n_ev)
            q = np.nanpercentile(s[draws].sum(axis=1) / np.maximum(n[draws].sum(axis=1), 1), J["interval"])
            rows.append({"points": label, "length_km": f"{lo:g}-{hi:g}", "n_points": int(m.sum()), "links": int(le.cml_id[(le.L_km >= lo) & (le.L_km < hi)].nunique()),
                         "mean_residual_db": float(res[m].mean()), "lo": float(q[0]), "hi": float(q[1]), "rmse_db": float(np.sqrt((res[m] ** 2).mean())),
                         "mean_observed_A_minus_c_db": float((D["y5"] - D["le_c"][D["bin_le"]])[m].mean())})
    return pd.DataFrame(rows)


# =================================================================================================== figures
def figure_ladder(ladders, chosen_of):
    names = list(ladders)
    fig, axes = plt.subplots(1, len(names), figsize=(3.5 * len(names), 4.6), squeeze=False)
    for ax, name in zip(axes[0], names):
        tab = ladders[name]
        x = np.arange(len(tab))
        ax.errorbar(x, tab.rmse, yerr=[tab.rmse - tab.rmse_lo, tab.rmse_hi - tab.rmse], fmt="o", color=INK, ms=6, capsize=4, lw=1.3)
        for i, r in enumerate(tab.itertuples()):
            ax.text(i, r.rmse_hi, f"{r.rmse:.3f}" + ("\nCHOSEN" if r.model == chosen_of[name] else ""), fontsize=7.5, ha="center", va="bottom",
                    fontweight="bold" if r.model == chosen_of[name] else "normal")
            if r.model == chosen_of[name]:
                ax.plot(i, r.rmse, "s", ms=14, mfc="none", mec=INK, mew=1.5)
        ax.set(xticks=x, xlim=(-0.6, len(tab) - 0.4))
        ax.set_xticklabels(tab.model, fontsize=8)
        ax.set_title(SET_LABEL[name], fontsize=9)
        ax.margins(y=0.25)
        style_axes(ax)
    axes[0][0].set_ylabel("held-out RMSE of A (dB), mean over events")
    fig.suptitle("Step S2: leave-one-event-out error per length set and rain input (dot = mean over the held-out events, whisker = 90 % event-bootstrap\n"
                 "interval; square = the model chosen by the stop rule on that set). The sets have different link-events: compare inside a panel.", fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    return fig


def figure_length(length_tab, models):
    keys = list(dict.fromkeys(length_tab.key))
    fig, axes = plt.subplots(1, len(models), figsize=(7.2 * len(models), 0.36 * len(keys) + 2.6), squeeze=False)
    sets = ["S1", "S2", "S3"]
    for ax, model in zip(axes[0], models):
        t = length_tab[length_tab.model == model].set_index("key")
        ks = [k for k in keys if k in t.index]
        for i, k in enumerate(ks):
            for s_, (mk, off) in zip(sets, (("o", 0.22), ("s", 0.0), ("^", -0.22))):
                v, lo, hi = t.loc[k, f"{s_}_value"], t.loc[k, f"{s_}_lo"], t.loc[k, f"{s_}_hi"]
                if np.isfinite(v):
                    ax.errorbar(v, i + off, xerr=[[max(v - lo, 0)], [max(hi - v, 0)]] if np.isfinite(lo) else None, fmt=mk, color=INK, ms=5,
                                mfc="none" if s_ != "S3" else INK, capsize=3, lw=1.0, label=SET_LABEL[s_] if i == 0 else None)
        ax.set(yticks=range(len(ks)), xscale="log", xlabel="value (minutes for the time constants, dB for the levels)")
        ax.set_yticklabels(ks, fontsize=8)
        ax.invert_yaxis()
        ax.set_title(f"{model}: circle = S<1, square = S<2, filled triangle = S<3", fontsize=9)
        style_axes(ax)
    fig.suptitle("Step S2, the length rule: the key WAA numbers on the three nested length sets (radar input), with 90 % event-bootstrap intervals", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


def figure_curves(model, bands):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0))
    level = "Weq" if model.startswith("M5") else "F"
    for ax, name, unit in zip(axes, (level, "P"), ("dB", "dB / km")):
        xs = np.array(KEY_R) if name == "Weq" else KNOTS[1:]
        for s_, lw, fill in (("S1", 1.0, "none"), ("S3", 2.0, INK)):
            par = result(f"real_{s_}_{model}")["params"]
            for (ls, mk), band in zip(BAND_STYLE, bands):
                v = np.array([par.get(f"{name}[{band}]({k:g})", np.nan) for k in xs])
                if not np.isfinite(v).any():
                    continue
                ax.plot(xs, v, ls=ls, marker=mk, ms=5, mfc=fill, color=INK, lw=lw, label=f"{band}, {SET_LABEL[s_][:3]} ({'thick, filled' if s_ == 'S3' else 'thin, open'})")
                q = np.array([boot_interval(f"real_{s_}_{model}", f"{name}[{band}]({k:g})") for k in xs])
                ax.vlines(xs * (1.04 if s_ == "S3" else 0.96), q[:, 0], q[:, 1], color=INK, lw=lw * 0.6)
        ax.set_xscale("symlog", linthresh=0.5)
        ax.set_xticks(KNOTS), ax.set_xticklabels([f"{k:g}" for k in KNOTS])
        ax.set(xlabel="rain R (mm/h)", ylabel=f"{'level reached' if name == 'Weq' else name + '(R)'} ({unit})", xlim=(0, 45))
        ax.legend(frameon=False, fontsize=7.5, loc="upper left")
        style_axes(ax)
    axes[0].set_title("the WAA level the rain drives toward (a_link = 1)", fontsize=9.5)
    axes[1].set_title("P: the path term per km, learned only through the path length", fontsize=9.5)
    fig.suptitle(f"Step S2: the two curves of {model} on S<1 (thin, open markers) and on S<3 (thick, filled markers), per band (line style and marker);\n"
                 "vertical bars = 90 % event-bootstrap interval at each knot", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def figure_gallery(D, le, picks, preds, chosen):
    fig = plt.figure(figsize=(17, 9.6))
    gs = fig.add_gridspec(4, 4, height_ratios=[3, 1.1, 3, 1.1], hspace=0.55, wspace=0.22, top=0.9, bottom=0.06, left=0.045, right=0.985)
    for j, e in enumerate(picks):
        ax, bx = fig.add_subplot(gs[2 * (j // 4), j % 4]), fig.add_subplot(gs[2 * (j // 4) + 1, j % 4])
        o, n, r = D["le_boff5"][e], D["le_n"][e] // 5, le.iloc[e]
        x = (np.arange(n) - r.last_bin - 1) * 5
        w = D["w5"][o:o + n] > 0
        ax.plot(x[w], D["y5"][o:o + n][w], "o", ms=2.6, color="#8a8a8a")
        ax.plot(x[w], preds[chosen][o:o + n][w], "-", color=INK, lw=1.5)
        if chosen != "M3":
            ax.plot(x[w], preds["M3"][o:o + n][w], "--", color=INK, lw=1.0)
        ax.axhline(r.offset_c_db, color=INK2, lw=0.6)
        rain1 = D["rain"][D["le_roff"][e] + PADM:D["le_roff"][e] + PADM + n * 5]
        bx.plot((np.arange(n * 5) - (r.last_bin + 1) * 5), rain1, "-", color=INK, lw=0.8)
        for a_ in (ax, bx):
            a_.axvline(0, color=INK, lw=0.8)
            a_.set_xlim(x[0], x[-1])
            style_axes(a_)
        ax.set_title(f"{r.sublink_id}, event {r.event_id[1:]}\nL = {r.L_km:.2f} km, {r.band_name}, gauge {r.gauge_id} at {r.gauge_distance_m:.0f} m", fontsize=7.5)
        if j % 4 == 0:
            ax.set_ylabel("A (dB)"), bx.set_ylabel("mm/h")
        if j // 4 == 1:
            bx.set_xlabel("minutes since the last wet bin of the link")
    fig.suptitle(f"Step S2: 8 held-out link-events with the GAUGE rain as input (seed {S['gallery']['seed']}, observed peak >= {CHK['min_peak_db']} dB). Top: observed A (dots, 5 min means), "
                 f"the chosen model {chosen} (solid){'' if chosen == 'M3' else ' and M3 (dashed)'},\npredicted from the gauge rain alone with the event left out of the fit; thin line = the offset c. "
                 "Bottom: the 1 min rain of the City gauge. Vertical line: the end of the rain on the link.", fontsize=9.5)
    return fig


def figure_after(curve, models):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    for ax, logy in zip(axes, (False, True)):
        ax.plot(curve.minutes_after_last_wet_bin, curve.observed, "o", ms=4.5, color=INK, mfc="none", label="observed (circles)")
        for m in models:
            ax.plot(curve.minutes_after_last_wet_bin, curve[m], ls=LINE[m], color=INK, lw=1.5, marker=MARK[m], ms=3.5, markevery=6, label=f"{m}, held out")
        ax.set(xlabel="minutes after the last wet bin of the link", title="the same, logarithmic" if logy else "median (A - c) / observed peak of the last shower")
        if logy:
            ax.set_yscale("log")
            ax.set_ylim(3e-3, 0.5)
        ax.legend(frameon=False, fontsize=8.5)
        style_axes(ax)
    fig.suptitle(f"Step S2 (S<3, radar input): after the last shower ({int(curve.link_events.iloc[0])} link-events with an observed peak >= {CHK['min_peak_db']} dB), observed against "
                 "predicted from the rain alone", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def figure_m5(par, lam, naming=None):
    """The surfaces of M5; if a formula was named, M5n at the same knots as small grey crosses (no line)."""
    rk, wk = np.asarray(M["m5"]["r_knots_mmh"], float), np.asarray(M["m5"]["w_knots_db"], float)
    G = np.array([[par[f"g({r:g},{w_:g})"] for w_ in wk] for r in rk[1:]])
    H = np.r_[0.0, [par[f"h({w_:g})"] for w_ in wk[1:]]]
    named = naming is not None and naming.get("done")
    styles = ["-", "--", "-.", ":", (0, (5, 1, 1, 1, 1, 1)), (0, (1, 1))]
    marks = ["o", "s", "^", "v", "D", "x"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.9))
    for k, w_ in enumerate(wk):
        axes[0].plot(rk, np.r_[0.0, G[:, k]], ls=styles[k], marker=marks[k], ms=4.5, mfc="none", color=INK, lw=1.4, label=f"W = {w_:g} dB")
    axes[0].set_xscale("symlog", linthresh=0.5)
    axes[0].set_xticks(rk), axes[0].set_xticklabels([f"{r:g}" for r in rk])
    axes[0].set(xlabel="rain R (mm/h)", ylabel="g(R, W) (dB / min)", title="g against the rain, one line per W knot")
    for j, r in enumerate(rk[1:]):
        axes[1].plot(wk, G[j], ls=styles[j % 6], marker=marks[j % 6], ms=4.5, mfc="none", color=INK, lw=1.4, label=f"R = {r:g} mm/h")
    axes[1].plot(wk, H, "-", color=INK, lw=3.0, alpha=0.45, label="h(W) (thick grey)")
    axes[1].set(xlabel="W (dB)", ylabel="dB / min", title="g against W, one line per rain knot; where a line crosses h(W) the level settles")
    axes[2].plot(wk, H, "-o", color=INK, lw=1.7, ms=5, label="M5")
    if named:
        for k in range(len(wk)):
            axes[0].plot(rk, naming["Gn"][:, k], "+", color=INK2, ms=7, label="M5n (named formula), at the knots" if k == 0 else None)
        axes[2].plot(wk, naming["Hn"], "--+", color=INK2, ms=8, lw=1.2, label="M5n: " + naming["h_text"])
        axes[2].legend(frameon=False, fontsize=7.5, loc="upper left")
    with np.errstate(all="ignore"):
        for w_, h_ in zip(wk[1:], H[1:]):
            axes[2].text(w_, h_, f"  W / h = {w_ / h_:.0f} min", fontsize=7.5, va="top")
    axes[2].set(xlabel="W (dB)", ylabel="h(W) (dB / min)", title="h: the drying rate without rain")
    for ax in axes:
        style_axes(ax)
    axes[0].legend(frameon=False, fontsize=8), axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle(f"Step S2: the surfaces of the discovery rung M5 (S<3, radar input, smoothness lambda = {lam:g}):  dW/dt = s_band a_link g(R, W) - h(W)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    if "--pysr" in sys.argv:                               # the PySR part of the naming stage (run by main in its own process)
        return run_pysr()
    CACHE.mkdir(parents=True, exist_ok=True)
    clock = time.time()
    workers = max((os.cpu_count() or 2) - 1, 1)

    # ------------------------------------------------------------------------------------------- data
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    bands_cfg, band_notes, band_links = merged_bands(d)
    BANDS = list(bands_cfg)
    D3, le, drops = sid.build_dataset(d, {**DATA, "bands_ghz": bands_cfg})
    if SMOKE:
        thin = np.arange(len(le)) % 12 == 0
        D3, le = sid.subset(D3, thin, renumber_events=False), le[thin].reset_index(drop=True)
    L = le.L_km.to_numpy()
    has_g = (le.gauge_id != "").to_numpy()
    keep = {"S1": L < DATA["length_sets_km"]["S<1"], "S2": L < DATA["length_sets_km"]["S<2"], "S3": np.ones(len(le), bool), "GG": has_g, "GR": has_g}
    data = {"S3": D3, "S1": sid.subset(D3, keep["S1"]), "S2": sid.subset(D3, keep["S2"])}
    data["GG"], data["GR"] = gauge_variants(D3, has_g)
    les = {k: le[keep[k]].reset_index(drop=True) for k in data}
    files = {k: sid.save_data(v, CACHE / f"data_{k}.npz") for k, v in data.items()}
    n_ev = {k: int(v["le_event"].max()) + 1 for k, v in data.items()}
    minutes = {k: int(v["le_n"].sum()) for k, v in data.items()}
    counts = pd.DataFrame([row for k in ("S1", "S2", "S3", "GG", "GR") for row in count_rows(k, data[k], les[k])])
    save(counts, P.SID2_COUNTS, index=False)
    drops_tab = pd.DataFrame({"rule": [k.replace("the two bands", "the bands") for k in drops], "pairs": list(drops.values())})
    print(f"TRAIN only: S<3 {len(le)} link-events, {le.cml_id.nunique()} links, {n_ev['S3']} events; S<1 {keep['S1'].sum()}, S<2 {keep['S2'].sum()}, gauge links "
          f"{has_g.sum()} link-events; bands {band_links}; {workers} processes; loaded in {time.time() - clock:.0f} s", flush=True)

    def job(prefix, name, model, data_file=None, lam=None, tag=True, **kw):
        """A job: `model` on the dataset `name` (or on a fake dataset of that size). key = the name of its cache file."""
        j = {"key": f"{prefix}_{model}" + (f"_lam{lam:g}" if lam is not None and tag else ""), "data": data_file or files[name], "model": model, "bands": BANDS, "set": name, **kw}
        if lam is not None:
            j["lam"] = lam
        return j

    # ------------------------------------------------------------------------------------------- the fake datasets
    noise = sid.noise_table()
    src = {"gauge": "GG", "S<1": "S1", "S<3": "S3"}
    fake_files = {}

    def fake_file(truth, rep, variant=""):
        """The data file of a fake replicate (made on first use). variant: '' the data, 'smooth' (Ta), 'S1' (Tc: the S<1 part)."""
        key = f"data_fake_{truth}_{rep}{'_' + variant if variant else ''}"
        file = CACHE / f"{key}.npz"
        if not file.exists():
            spec = FK["truths"][truth]
            name = src[spec["set"]]
            Fd = make_fake(data[name], les[name], spec, [FK["seed_base"], list(FK["truths"]).index(truth), rep], noise)
            if variant == "smooth":
                Fd = smoothed_input(Fd, spec["smooth_noise"], [FK["seed_base"], 77, rep])
            elif variant == "S1":
                Fd = sid.subset(Fd, keep["S1"])
            sid.save_data(Fd, file)
        fake_files[key] = str(file)
        return str(file)

    # ------------------------------------------------------------------------------------------- time one fit, plan
    plan_file = CACHE / "plan.json"
    lam_fake = M["m5"]["smooth_fake"]

    def all_jobs(plan, lam_star=None, chosen=None, m5d=False):
        """Every job of the run for a plan. Stage 1 = everything that does not need the M5 smoothness or the chosen model."""
        reps = range(plan["fake_reps"])
        sets_radar = ["S3", "S1"] + (["S2"] if plan["s2"] else [])
        stage1, stage2 = [], []
        for name in sets_radar + ["GG", "GR"]:
            for m in LADDER:
                if m == "M5":
                    continue
                stage1.append(job(f"real_{name}", name, m, loo=True, boot=J["param_boot_reps"] if m == "M3" else 0))
        stage1 += [job("real_S3", "S3", "M5", lam=lam, loo=True) for lam in plan["smooth_grid"]]
        for rep in reps:
            stage1 += [job(f"fake_Ta_{rep}_exact", "GG", "M3", fake_file("Ta", rep)), job(f"fake_Ta_{rep}_smooth", "GG", "M3", fake_file("Ta", rep, "smooth"))]
            for truth in ("Tb", "Tb2"):
                stage1 += [job(f"fake_{truth}_{rep}", "S1", m, fake_file(truth, rep), loo=True) for m in ("M3D", "M4")]
            stage1 += [job(f"fake_Tc_{rep}_S3", "S3", "M3", fake_file("Tc", rep), boot=FK["truths"]["Tc"]["boot_reps"]),
                       job(f"fake_Tc_{rep}_S1", "S1", "M3", fake_file("Tc", rep, "S1"), boot=FK["truths"]["Tc"]["boot_reps"])]
            if plan["td"]:
                stage1 += [job(f"fake_Td_{rep}", "S1", "M5", fake_file("Td", rep), lam=lam_fake, tag=False, loo=True),
                           job(f"fake_Td_{rep}", "S1", "M3", fake_file("Td", rep), loo=True)]
        if lam_star is not None:
            for name in ["S1"] + (["S2"] if plan["s2"] else []) + ["GG", "GR"]:
                stage2.append(job(f"real_{name}", name, "M5", lam=lam_star, tag=False, loo=True))
            if m5d:
                stage2.append(job("real_S3", "S3", "M5D", lam=lam_star, tag=False, loo=True))
        if chosen not in (None, "M3", "M1"):
            for name in sets_radar + ["GG", "GR"]:
                kw = {"lam": lam_star, "tag": name == "S3"} if chosen == "M5" else {}
                stage2.append(job(f"real_{name}", name, chosen, loo=True, boot=J["param_boot_reps"], **kw))
        return stage1, stage2

    def cost_h(jobs, t_eval):
        """Projected hours of a list of jobs on the workers (t_eval = seconds per evaluation of M3 and of M5 on S<3)."""
        tot = sum(evaluations(j, n_ev[j["set"]]) * t_eval["M5" if j["model"].startswith("M5") else "M3"] * minutes[j["set"]] / minutes["S3"] for j in jobs if not done(j))
        return tot / workers / 3600

    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
        print(f"plan from {rel(plan_file)}: cut: {plan['cuts'] or 'nothing'}")
    else:
        ones, t_eval = np.ones(len(le)), {}
        for m in ("M3", "M5"):
            Pm = sid.Problem2(D3, m, M, lam=lam_fake)
            x0 = Pm.start(M["starts"]["g1"])
            Pm.run(x0, 0, ones)                                                        # compile
            t0 = time.time()
            for _ in range(5):
                Pm.run(x0, 0, ones)
            t_eval[m] = (time.time() - t0) / 5
        Pm = sid.Problem2(D3, "M3", M)
        t0 = time.time()
        fit = Pm.fit(Pm.start(M["starts"]["g1"]), 0, ones, M["maxiter"]["start"])
        t_fit = time.time() - t0
        plan = {"t_eval_ms": {m: round(v * 1e3, 1) for m, v in t_eval.items()}, "t_fit_s": round(t_fit, 1), "nfev": fit["nfev"], "cuts": [],
                "fake_reps": FK["replicates"], "smooth_grid": list(M["m5"]["smooth_grid"]), "s2": True, "td": True}

        def project():
            s1_, s2_ = all_jobs(plan, lam_star=plan["smooth_grid"][0], chosen="M5")      # the dearest case: M5 is chosen
            return cost_h(s1_ + s2_, t_eval)

        proj = project()
        plan["projected_full_h"] = round(proj, 2)
        print(f"timing on S<3: one M3 fit from a generic start {t_fit:.0f} s ({fit['nfev']} evaluations); one evaluation M3 {t_eval['M3'] * 1e3:.0f} ms, M5 {t_eval['M5'] * 1e3:.0f} ms; "
              f"projection {proj:.2f} h on {workers} processes if M5 is chosen (budget {S['time_budget_h']} h)", flush=True)
        for cut, change in (("the M5 smoothness grid to 2 values", {"smooth_grid": list(M["m5"]["smooth_grid_cut"])}), ("S<2", {"s2": False}),
                            ("fake Td", {"td": False}), ("fake replicates 3 -> 2", {"fake_reps": 2})):
            if proj <= S["time_budget_h"]:
                break
            plan.update(change), plan["cuts"].append(cut)
            proj = project()
            print(f"   cut: {cut} -> projected {proj:.2f} h", flush=True)
        plan["projected_h"], plan["t_eval_s"] = round(proj, 2), t_eval
        plan_file.write_text(json.dumps(plan))
    sets_radar = ["S1"] + (["S2"] if plan["s2"] else []) + ["S3"]
    sets_all = sets_radar + ["GG", "GR"]

    # ------------------------------------------------------------------------------------------- stage 1: fake checks, ladders without the chosen M5
    stage1, _ = all_jobs(plan)
    heavy = lambda j: -evaluations(j, n_ev[j["set"]]) * minutes[j["set"]] * (2.5 if j["model"].startswith("M5") else 1)   # noqa: E731
    run_pool(sorted(stage1, key=heavy), workers, "stage 1")

    # the smoothness of M5: the lowest held-out RMSE on S<3 (radar)
    lam_rows = []
    for lam in plan["smooth_grid"]:
        r = result(f"real_S3_M5_lam{lam:g}")
        e = sid.event_errors(D3, r["loo_pred"], J)
        lam_rows.append({"lambda": lam, "held_out_rmse": float(e.rmse.mean()), "training_loss": r["all"]["loss"], "t_wet": r["params"]["t_wet"], "t_dry": r["params"]["t_dry"]})
    lam_tab = pd.DataFrame(lam_rows)
    lam_star = float(lam_tab.loc[lam_tab.held_out_rmse.idxmin(), "lambda"])

    def sync_m5():
        """The M5 fit with the chosen lambda also under the plain name, so that every set has a file "real_<set>_M5"."""
        (CACHE / "real_S3_M5.pkl").write_bytes((CACHE / f"real_S3_M5_lam{lam_star:g}.pkl").read_bytes())
        _RESULTS.pop("real_S3_M5", None)

    sync_m5()

    # slow water or drift, on S<3 (radar): decides whether M5 is also fitted with a drift
    slow_real = slow_verdict("real_S3", D3)
    _, stage2 = all_jobs(plan, lam_star=lam_star, m5d=SP["m5_with_drift"] and slow_real["verdict"] != "WATER-LIKE")
    spent_h = sum(pickle.loads(f.read_bytes()).get("seconds", 0.0) for f in CACHE.glob("*.pkl") if f.name != "real_S3_M5.pkl") / workers / 3600
    m5d_job = [j for j in stage2 if j["model"] == "M5D"]
    if m5d_job and not done(m5d_job[0]) and spent_h + cost_h(stage2, plan["t_eval_s"]) > S["time_budget_h"]:
        stage2 = [j for j in stage2 if j["model"] != "M5D"]            # the optional job goes first if the run is late
        plan["cuts"] = plan["cuts"] + ["M5D (optional)"]
    run_pool(sorted(stage2, key=heavy), workers, "stage 2")
    has_m5d = (CACHE / "real_S3_M5D.pkl").exists() and done({"key": "real_S3_M5D", "loo": True})

    # ------------------------------------------------------------------------------------------- the ladders and the chosen model
    ladders, chosen_of, best_of, errors_of = {}, {}, {}, {}
    for name in sets_all:
        ladders[name], best_of[name], chosen_of[name], errors_of[name] = ladder(f"real_{name}", data[name], LADDER)
    chosen = chosen_of["S3"]
    _, stage3 = all_jobs(plan, lam_star=lam_star, chosen=chosen)
    stage3 = [j for j in stage3 if j.get("boot")]
    run_pool(sorted(stage3, key=heavy), workers, "stage 3 (bootstrap of the chosen model)")
    sync_m5()
    _RESULTS.clear()
    ref = list(dict.fromkeys(["M3", chosen])) if chosen != "M1" else ["M3"]

    ladder_tab = pd.concat([t.assign(set=SET_LABEL[k], best=best_of[k], chosen=chosen_of[k]) for k, t in ladders.items()], ignore_index=True)
    ladder_tab.insert(0, "set", ladder_tab.pop("set"))
    save(ladder_tab.round(4), P.SID2_LADDER, index=False)
    a0, a1 = J["after_windows_min"]
    late = f"rmse_{a1[0]}_{a1[1]}"
    best_late = {k: t.model[t[late].idxmin()] for k, t in ladders.items()}
    best_peak = {k: t.model[t.shower_peak_error_median.abs().idxmin()] for k, t in ladders.items()}

    da_tabs = {}
    for name in sets_all:
        models = LADDER + (["M5D"] if has_m5d and name == "S3" else [])
        da_tabs[name] = drift_ladder(f"real_{name}", data[name], models)[0]
    da_tab = pd.concat([t.assign(set=SET_LABEL[k]) for k, t in da_tabs.items()], ignore_index=True)
    da_tab.insert(0, "set", da_tab.pop("set"))
    save(da_tab.round(4), P.SID2_LADDER_DRIFT, index=False)

    fit_rows = []
    for name in sets_all:
        for m in LADDER + (["M5D"] if has_m5d and name == "S3" else []):
            r = result(f"real_{name}_{m}")
            st = r["all"]["starts_d0"]
            fit_rows.append({"set": SET_LABEL[name], "model": m, "d": r["all"]["d"], "loss": r["all"]["loss"], "spread_over_starts_at_d0": max(st.values()) - min(st.values()),
                             "best_start": min(st, key=st.get), "s_b_db_h": r["all"]["s_b"][-1] if m.endswith("D") else np.nan,
                             "folds_at_iteration_cap": int(sum(f["nit"] >= M["maxiter"]["fold"] for f in r["folds"])), "folds": len(r["folds"]), "seconds": round(r["seconds"])})
    fits = pd.DataFrame(fit_rows)
    save(fits.round(3), P.SID2_FITS, index=False)
    print(f"ladder on S<3 (radar): best {best_of['S3']}, chosen {chosen}; M5 lambda {lam_star:g}; slow part: {slow_real['verdict']}", flush=True)

    # ------------------------------------------------------------------------------------------- the fake checks
    tol = FK["tol"]
    Ftrue, itu = truth_curves(le, BANDS)
    reps = range(plan["fake_reps"])
    fake_rows, verdicts = [], {}
    within = lambda v, truth: abs(v / truth - 1) <= tol          # noqa: E731
    # Ta: the input explanation
    spec = FK["truths"]["Ta"]
    ta = pd.DataFrame([{"replicate": rep, "input": inp, **{k: result(f"fake_Ta_{rep}_{inp}_M3")["params"][k] for k in ("tau_wet", "tau_dry", "d")}} for rep in reps for inp in ("exact", "smooth")])
    med = ta.groupby("input")[["tau_wet", "tau_dry"]].median()
    ok_exact = within(med.tau_wet["exact"], spec["tau_wet"]) and within(med.tau_dry["exact"], spec["tau_dry"])
    slow = med.tau_wet["smooth"] >= spec["slow_factor"] * spec["tau_wet"]
    verdicts["Ta"] = {"PASS": bool(ok_exact and slow), "text": f"exact input: tau_wet {med.tau_wet['exact']:.2f} (true {spec['tau_wet']}), tau_dry {med.tau_dry['exact']:.2f} (true {spec['tau_dry']}): "
                      f"{'recovered' if ok_exact else 'NOT recovered'}; smoothed input: tau_wet {med.tau_wet['smooth']:.2f} = {med.tau_wet['smooth'] / spec['tau_wet']:.1f} x the truth "
                      f"({'at least' if slow else 'less than'} the factor {spec['slow_factor']:g} of the rule), tau_dry {med.tau_dry['smooth']:.2f}"}
    fake_rows += [{"truth": "Ta", "replicate": r.replicate, "what": f"M3 fitted with the {r.input} input", "tau_wet": r.tau_wet, "tau_dry": r.tau_dry, "d": r.d} for r in ta.itertuples()]
    # Tb / Tb2: drift against a slow store
    slow_rows = [{"data": "REAL S<3, radar", **slow_real}]
    tb = {}
    for truth in ("Tb", "Tb2"):
        for rep in reps:
            v = slow_verdict(f"fake_{truth}_{rep}", sid.load_data(fake_file(truth, rep), BANDS))
            slow_rows.append({"data": f"fake {truth} ({'M3 + drift' if truth == 'Tb' else 'M4, slow store'}), replicate {rep}", **v})
            tb.setdefault(truth, []).append(v["verdict"])
            fake_rows.append({"truth": truth, "replicate": rep, "what": "slow-part rule", "verdict": v["verdict"], "tau_s": v["tau_s_M4_min"], "s_b": v["s_b_db_h"]})
    ok_b = sum(v == "DRIFT-LIKE" for v in tb["Tb"]) >= 2 and sum(v == "WATER-LIKE" for v in tb["Tb2"]) >= 2
    verdicts["Tb"] = {"PASS": bool(ok_b), "text": f"truth M3 + drift read as: {', '.join(tb['Tb'])}; truth M4 (slow store, tau_s {FK['truths']['Tb2']['tau_s']} min) read as: {', '.join(tb['Tb2'])}"}
    slow_tab = pd.DataFrame(slow_rows)
    save(slow_tab.round(4), P.SID2_SLOW, index=False)
    # Tc: the length sets
    spec = FK["truths"]["Tc"]
    tc_rows = []
    for rep in reps:
        for s_ in ("S1", "S3"):
            key = f"fake_Tc_{rep}_{s_}_M3"
            par = result(key)["params"]
            for bi, band in enumerate(BANDS):
                for kind, rr in (("F", KEY_R), ("P", [2, 5, 10])):
                    for r_ in rr:
                        name = f"{kind}[{band}]({r_:g})"
                        truth = float(Ftrue(r_)) if kind == "F" else itu[bi][0] * r_ ** itu[bi][1]
                        lo, hi = boot_interval(key, name)
                        tc_rows.append({"replicate": rep, "set": s_, "kind": kind, "name": name, "value": par.get(name, np.nan), "truth": truth, "rel_width": (hi - lo) / truth})
            fake_rows.append({"truth": "Tc", "replicate": rep, "what": f"M3 fitted on {SET_LABEL[s_][:3]}", "tau_wet": par["tau_wet"], "tau_dry": par["tau_dry"], "d": par["d"]})
    tc = pd.DataFrame(tc_rows)
    tc_med = tc.groupby(["set", "kind", "name"], as_index=False).agg(value=("value", "median"), truth=("truth", "first"), rel_width=("rel_width", "median"))
    tc_med["ratio"] = tc_med.value / tc_med.truth
    in_s1 = tc_med[(tc_med.set == "S1") & tc_med.value.notna()].name
    s3_ok = bool((tc_med[tc_med.set == "S3"].ratio.sub(1).abs() <= tol).all())
    wid = tc_med[tc_med.name.isin(in_s1)].groupby(["set", "kind"]).rel_width.mean().unstack(0)
    narrower = bool((wid["S3"] < wid["S1"]).all())
    worst = tc_med[tc_med.set == "S3"].assign(off=lambda t: t.ratio.sub(1).abs()).sort_values("off").iloc[-1]
    s1_off = tc_med[(tc_med.set == "S1") & tc_med.value.notna()].ratio.sub(1).abs()
    verdicts["Tc"] = {"PASS": bool(s3_ok and narrower), "text": f"on S<3 every F and P number is {'within' if s3_ok else 'NOT within'} +-{tol:.0%} of the truth (worst: {worst['name']} x{worst.ratio:.2f}); "
                      f"mean relative interval width F {wid.loc['F', 'S1']:.2f} (S<1) -> {wid.loc['F', 'S3']:.2f} (S<3), P {wid.loc['P', 'S1']:.2f} -> {wid.loc['P', 'S3']:.2f}; "
                      f"on S<1 {int((s1_off > tol).sum())} of {len(s1_off)} numbers are off by more than {tol:.0%} (farthest: x{tc_med.ratio[s1_off.idxmax()]:.3f}, {tc_med.name[s1_off.idxmax()]})"}
    # Td: the discovery rung on M3 data
    if plan["td"]:
        spec = FK["truths"]["Td"]
        td_rows, td_ok = [], 0
        for rep in reps:
            Dd = sid.load_data(fake_file("Td", rep), BANDS)
            e5, e3 = sid.event_errors(Dd, result(f"fake_Td_{rep}_M5")["loo_pred"], J).rmse.mean(), sid.event_errors(Dd, result(f"fake_Td_{rep}_M3")["loo_pred"], J).rmse.mean()
            par = result(f"fake_Td_{rep}_M5")["params"]
            b0 = BANDS[int(np.argmax(np.bincount(Dd["le_band"], minlength=len(BANDS)) > 0))]
            td_rows.append({"replicate": rep, "rmse_M5": e5, "rmse_M3": e3, "t_wet": par["t_wet"], "t_dry": par["t_dry"], **{f"level({r_:g})": par[f"Weq[{b0}]({r_:g})"] for r_ in KEY_R}})
            td_ok += e5 <= e3 * (1 + spec["rmse_margin"])
            fake_rows.append({"truth": "Td", "replicate": rep, "what": "M5 fitted on M3 data", "tau_wet": par["t_wet"], "tau_dry": par["t_dry"], "d": par["d"], "rmse_M5": e5, "rmse_M3": e3})
        td = pd.DataFrame(td_rows)
        tm = td.median()
        truths = {"t_wet": spec["tau_wet"], "t_dry": spec["tau_dry"], **{f"level({r_:g})": float(Ftrue(r_)) for r_ in KEY_R}}
        ratios = {k: tm[k] / v for k, v in truths.items()}
        par_ok = all(abs(r_ - 1) <= tol for r_ in ratios.values())
        verdicts["Td"] = {"PASS": bool(par_ok and td_ok >= 2), "text": "M5 / truth: " + ", ".join(f"{k} x{r_:.2f}" for k, r_ in ratios.items())
                          + f"; held-out RMSE of M5 within {spec['rmse_margin']:.0%} of M3 in {td_ok} of {len(td)} replicates (M5 {tm.rmse_M5:.4f}, M3 {tm.rmse_M3:.4f} dB)"}
    else:
        verdicts["Td"] = {"PASS": False, "text": "not run (time cut)"}
    fake = pd.DataFrame(fake_rows)
    for t_, v in verdicts.items():
        fake.loc[fake.truth.str.startswith(t_), "PASS"] = v["PASS"]
    save(fake.round(4), P.SID2_FAKE, index=False)
    print("fake checks: " + "; ".join(f"{t_} {'PASS' if v['PASS'] else 'FAIL'}" for t_, v in verdicts.items()), flush=True)
    tag = {"input": "" if verdicts["Ta"]["PASS"] else " [not validated: fake check Ta failed]", "slow": "" if verdicts["Tb"]["PASS"] else " [not validated: fake check Tb failed]",
           "length": "" if verdicts["Tc"]["PASS"] else " [not validated: fake check Tc failed]", "m5": "" if verdicts["Td"]["PASS"] else " [not validated: fake check Td failed]"}

    # ------------------------------------------------------------------------------------------- the length rule
    len_rows = []
    for m in ref:
        taus, levels = key_names(m, BANDS)
        for k in taus + levels:
            row = {"model": m, "key": k}
            for s_ in ("S1", "S2", "S3"):
                has = s_ in sets_radar
                row[f"{s_}_value"] = result(f"real_{s_}_{m}")["params"].get(k, np.nan) if has else np.nan
                row[f"{s_}_lo"], row[f"{s_}_hi"] = boot_interval(f"real_{s_}_{m}", k) if has else (np.nan, np.nan)
            row["in_rule"] = bool(np.isfinite(row["S1_value"]))
            row["S3_inside_S1_interval"] = bool(row["S1_lo"] <= row["S3_value"] <= row["S1_hi"]) if row["in_rule"] else None
            len_rows.append(row)
    length_tab = pd.DataFrame(len_rows)
    save(length_tab.round(4), P.SID2_LENGTH, index=False)
    rule_model = chosen if chosen != "M1" else "M3"       # M1 has no time constants: M3 decides
    in_rule = length_tab[(length_tab.model == rule_model) & length_tab.in_rule]
    agree = bool(in_rule.S3_inside_S1_interval.all())
    outside = in_rule[~in_rule.S3_inside_S1_interval.astype(bool)]
    main_set = "S3" if agree else "S1"
    resid = residuals_by_length(D3, le, result(f"real_S3_{chosen}")["loo_pred"])
    save(resid.round(4), P.SID2_RESID, index=False)
    r_all = resid[resid.points == "all wet points"]
    trend = bool(r_all.iloc[0].hi < r_all.iloc[-1].lo or r_all.iloc[-1].hi < r_all.iloc[0].lo)

    # ------------------------------------------------------------------------------------------- the gauge rule
    gauge_rows = []
    for m in ref:
        taus, levels = key_names(m, BANDS)
        for name in ("GG", "GR"):
            r, t = result(f"real_{name}_{m}"), ladders[name].set_index("model")
            row = {"model": m, "input": "GAUGE" if name == "GG" else "radar, same points", "held_out_rmse": t.rmse[m], "rmse_lo": t.rmse_lo[m], "rmse_hi": t.rmse_hi[m],
                   "rmse_rain": t.rmse_rain[m], "shower_peak_error_median": t.shower_peak_error_median[m], "peak_lo": t.shower_peak_error_lo[m], "peak_hi": t.shower_peak_error_hi[m], "d": r["params"]["d"]}
            for k, label in zip(taus, ("wetting_min", "drying_min")):
                lo, hi = boot_interval(f"real_{name}_{m}", k)
                row[label], row[label + "_lo"], row[label + "_hi"] = r["params"][k], lo, hi
            for k in levels:
                if np.isfinite(r["params"].get(k, np.nan)):
                    row[k] = r["params"][k]
            row["P_10_first_band"] = r["params"].get(f"P[{BANDS[0]}](10)", np.nan)
            gauge_rows.append(row)
    gauge_tab = pd.DataFrame(gauge_rows)
    save(gauge_tab.round(4), P.SID2_GAUGE, index=False)
    gt = gauge_tab.set_index(["model", "input"])
    g_model = rule_model
    dry_g, dry_r = gt.drying_min[(g_model, "GAUGE")], gt.drying_min[(g_model, "radar, same points")]
    dry_agree = bool(abs(dry_g / dry_r - 1) <= GR["dry_agree_tol"])

    # ------------------------------------------------------------------------------------------- parameters, the after-rain curve, the M5 surfaces
    par_rows = []
    for name in sets_all:
        for m in ref:
            r = result(f"real_{name}_{m}")
            for k, v in r["params"].items():
                if k.startswith(("g(", "h(")) and name != "S3":
                    continue
                lo, hi = boot_interval(f"real_{name}_{m}", k)
                par_rows.append({"set": SET_LABEL[name], "model": m, "parameter": k, "value": v, "lo": lo, "hi": hi})
    params = pd.DataFrame(par_rows)
    save(params.round(4), P.SID2_PARAMS, index=False)

    loo3 = {m: result(f"real_S3_{m}")["loo_pred"] for m in LADDER}
    curve = after_curves(D3, {m: loo3[m] for m in LADDER if m != "M1"})
    save(curve.round(4), P.SID2_AFTER, index=False)

    par5 = result("real_S3_M5")["params"]
    rk, wk = M["m5"]["r_knots_mmh"], M["m5"]["w_knots_db"]
    surf = pd.DataFrame([{"surface": "g", "R_mmh": r_, "W_db": w_, "M5_db_per_min": par5[f"g({r_:g},{w_:g})"]} for r_ in rk[1:] for w_ in wk]
                        + [{"surface": "h", "R_mmh": np.nan, "W_db": w_, "M5_db_per_min": par5[f"h({w_:g})"]} for w_ in wk[1:]])
    t3 = ladders["S3"].set_index("model")
    m5_qualifies = bool(chosen == "M5" or t3.interval_holds_0["M5"])
    naming = {"qualifies": m5_qualifies, "done": False}
    if m5_qualifies and not SMOKE:
        naming = name_m5(D3, lam_star, n_ev["S3"])
        surf = naming["surf"]
        save(naming["terms"].round(5), P.SID2_NAMING_TERMS, index=False)
    save(surf.round(5), P.SID2_NAMING, index=False)

    # ------------------------------------------------------------------------------------------- product and figures
    lt = sid.link_event_table(D3, loo3[chosen], J)
    product = pd.concat([le.drop(columns=["band"]), pd.DataFrame({"in_S1": keep["S1"], "in_S2": keep["S2"], "has_city_gauge_1km": has_g}),
                         lt.add_suffix(f"_{chosen}").rename(columns={f"points_{chosen}": "points", f"wet_points_{chosen}": "wet_points", f"peak_obs_{chosen}": "peak_obs_db",
                                                                     f"integral_obs_{chosen}": "integral_obs_db_min"})], axis=1)
    if chosen != "M3":
        product = pd.concat([product, sid.link_event_table(D3, loo3["M3"], J)[["peak_pred", "integral_pred", "bias", "rmse"]].add_suffix("_M3")], axis=1)
    product.insert(0, "chosen_model", chosen)
    num = product.select_dtypes("number").columns
    product[num] = product[num].round(4)
    save(product, P.SID2_LINK_EVENTS, index=False)

    fig = figure_ladder(ladders, chosen_of)
    save(fig, P.SID2_FIG_LADDER, dpi=110), plt.close(fig)
    fig = figure_length(length_tab, ref)
    save(fig, P.SID2_FIG_LENGTH, dpi=110), plt.close(fig)
    curve_model = chosen if chosen != "M1" else "M3"
    fig = figure_curves(curve_model, BANDS)
    save(fig, P.SID2_FIG_CURVES, dpi=110), plt.close(fig)
    Dg, leg = data["GG"], les["GG"]
    g_preds = {m: result(f"real_GG_{m}")["loo_pred"] for m in dict.fromkeys(ref + [chosen])}
    ltg = sid.link_event_table(Dg, g_preds[chosen], J)
    cand = np.nonzero((ltg.peak_obs - leg.offset_c_db >= CHK["min_peak_db"]).to_numpy())[0]
    picks = np.sort(np.random.default_rng(S["gallery"]["seed"]).choice(cand, S["gallery"]["link_events"], replace=False))
    fig = figure_gallery(Dg, leg, picks, g_preds, chosen)
    save(fig, P.SID2_FIG_GALLERY, dpi=100), plt.close(fig)
    fig = figure_after(curve, [m for m in LADDER if m != "M1"])
    save(fig, P.SID2_FIG_AFTER, dpi=110), plt.close(fig)
    fig = figure_m5(par5, lam_star, naming)
    save(fig, P.SID2_FIG_M5, dpi=110), plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    import numba
    import scipy
    show = ["set", "model", "rmse", "rmse_lo", "rmse_hi", "minus_best", "diff_lo", "diff_hi", "interval_holds_0"]
    more = ["set", "model", "rmse_rain", f"rmse_{a0[0]}_{a0[1]}", late, "shower_peak_error_median", "shower_peak_error_lo", "shower_peak_error_hi", "showers", "integral_pred_over_obs"]
    da_show = ["set", "model", "rmse", "minus_best", "diff_lo", "diff_hi", late, "integral_pred_over_obs"]
    seconds = sum(pickle.loads(f.read_bytes()).get("seconds", 0.0) for f in CACHE.glob("*.pkl") if f.name != "real_S3_M5.pkl")
    key_par = params[~params.parameter.str.startswith(("F[", "P[", "g(", "h(", "Weq[")) & params.set.isin([SET_LABEL[main_set], SET_LABEL["GG"]])].round(3)
    curve_par = params[params.parameter.str.startswith(("F[", "P[", "Weq[", "s[")) & (params.model == curve_model) & (params.set == SET_LABEL[main_set])].dropna(subset=["value"]).round(3)
    ctx = dict(locals())
    box, slow_text, plain, unclear = texts(ctx)
    cut_text = ", ".join(plan["cuts"]) or "nothing"
    not_val = [v for v in tag.values() if v]
    report = f"""# wet_antenna step S2: the WAA law on links up to 3 km, the slow part, and the discovery rung M5 (TRAIN only)

Written by `research/wet_antenna/system_id_s2.py` (library `system_id.py`). Settings and every rule:
`wet_antenna.system_id_s2` in `research/config.yaml`, written before the run. Plan: `docs/ROADMAP.md`, "Reboot (2026-10-13)".
{chr(10) + '**Results that carry the label "not validated" below: ' + '; '.join(v.strip(' []') for v in not_val) + '.**' + chr(10) if not_val else ''}
## 1. What was run

- Every model predicts A_pred(t) = W(t) + L * P_band(R(t - d)) + c at 1 min and is compared with the 5 min means of A.
  M1 static W = a_link F_band(R); M3 one store with a wetting and a drying time constant; M3D = M3 plus a linear level
  drift per link-event (not driven by rain, slopes ~ N(0, s_b^2), s_b fitted); M4 a fast and a slow store; M5 the
  discovery rung dW/dt = s_band a_link g(R, W) - h(W) with g a free surface on {len(rk)} x {len(wk)} knots and h a free rising
  curve. P and F are free piecewise-linear curves (0 at R = 0, never falling). **No power-law rain term is subtracted
  anywhere**; ITU-R P.838 only builds the fake data and is a size reference.
- Loss, starts, ridge on a_link as S1. The fit is cheaper than in S1 (config `maxiter`): lag grid {M['lag_grid_min'][0]}..{M['lag_grid_min'][1]} min, the
  folds keep the lag of the fit on all events.
- Judge: leave-one-event-out, {J['boot_reps']} event-bootstrap draws, the S1 stop rule on the order {' < '.join(J['order'])}.
- TRAIN events only; VAL is never read. Nothing earlier is changed (the data and judge functions of S1 moved into the
  library unchanged; the S1 outputs are identical).
- Timing: one evaluation on S<3 takes {plan['t_eval_ms']['M3']:.0f} ms (M3) and {plan['t_eval_ms']['M5']:.0f} ms (M5); projection {plan['projected_full_h']} h on {workers} processes
  against a budget of {S['time_budget_h']} h. Cut: {cut_text}. Fitting time summed over the processes: {seconds / 3600:.1f} h. Run time
  of this run: {(time.time() - clock) / 60:.1f} min (fits are cached).
- Software: numba {numba.__version__}, scipy {scipy.__version__} (L-BFGS-B, exact gradients), numpy {np.__version__}.

## 2. Counts

Bands (TRAIN links below {DATA['max_length_km']:g} km): {', '.join(f'{k} {v} links' for k, v in band_links.items())}. {'Merged: ' + '; '.join(band_notes) + '.' if band_notes else f"No band was merged (rule: at least {DATA['min_band_links']} links)."}

{md_table(counts)}

(`{rel(P.SID2_COUNTS)}`.) A link-event = one TRAIN event on one sub-link. How the pairs of S<3 were reduced (the S1 rules):

{md_table(drops_tab)}

## 3. The fake checks (first)

Truth of all fake data: F = {FK['f_top_db']} (1 - exp(-R / {FK['f_scale_mmh']:g})) dB, P = ITU-R P.838 size, a_link log-normal (sd {FK['a_sd']}), tau_wet 3 min,
tau_dry 8 min, noise variance x {FK['noise_variance_factor']:g}, {FK['rsl_step_db']} dB quantization, unknown offsets; {plan['fake_reps']} replicates per truth.

{chr(10).join(f"- **{t_} ({what}): {'PASS' if v['PASS'] else 'FAIL'}.** {v['text']}." for (t_, v), what in zip(verdicts.items(), ('the input explanation: gauge rain as the truth, a smoothed input for the fit', 'drift against a slow store', 'F and P on S<1 against S<3', 'the discovery rung on M3 data')))}

{md_table(fake.round(3))}

(`{rel(P.SID2_FAKE)}`.) Ta per input (median over the replicates):

{md_table(med.reset_index().round(2))}

Tc, median over the replicates (value / truth and the relative width of the 90 % interval):

{md_table(tc_med.round(3))}

## 4. The ladder per set and input

{md_table(ladder_tab[show].round(4))}

(`{rel(P.SID2_LADDER)}`; rmse = mean over the held-out events of the RMSE of A (dB), predicted from the rain alone;
minus_best with the 90 % interval of the paired difference. A held-out link-event of M3D has no drift.) More metrics
(shower peak error = predicted - observed peak, median over the showers with an observed peak >= {J['shower_min_peak_db']} dB):

{md_table(ladder_tab[more].round(3))}

M5: the smoothness lambda was chosen on S<3 (radar) by the held-out RMSE:

{md_table(lam_tab.round(4))}

The fits (`{rel(P.SID2_FITS)}`):

{md_table(fits.round(3))}

**Chosen by the stop rule on S<3 (radar): {chosen}** (best: {best_of['S3']}). Per set: {', '.join(f'{SET_LABEL[k]}: best {best_of[k]}, chosen {chosen_of[k]}' for k in sets_all)}.
Best on the 60-180 min window: {', '.join(f'{SET_LABEL[k]}: {v}' for k, v in best_late.items())}. Smallest shower peak error: {', '.join(f'{SET_LABEL[k]}: {v}' for k, v in best_peak.items())}.

![ladder](../figures/{P.SID2_FIG_LADDER.name})

## 5. The length rule{tag['length']}

Key numbers of {' and '.join(ref)} on the three nested sets (radar input; 90 % event-bootstrap intervals):

{md_table(length_tab.round(3))}

(`{rel(P.SID2_LENGTH)}`.) **The rule, read on {rule_model}: {'AGREE' if agree else 'DISAGREE'}.** {f'Every key number of S<3 lies inside the S<1 interval: S<3 is the main result.' if agree else f'{len(outside)} of {len(in_rule)} key numbers of S<3 lie outside the S<1 interval ({", ".join(outside.key)}): S<1 is the main result, and the difference is a finding for Stage 2.'}{tag['length']}

![length rule](../figures/{P.SID2_FIG_LENGTH.name})

![curves](../figures/{P.SID2_FIG_CURVES.name})

Held-out residuals of {chosen} on S<3 during rain, per length bin (predicted - observed):

{md_table(resid.round(3))}

(`{rel(P.SID2_RESID)}`.) **{'TREND with L: the intervals of the shortest and the longest bin do not overlap, so L * P(R) is not enough.' if trend else 'No trend by the rule: the intervals of the shortest and the longest bin overlap.'}**

## 6. The gauge rule{tag['input']}

{md_table(gauge_tab.round(3))}

(`{rel(P.SID2_GAUGE)}`. The same link-events and 5 min points for both inputs.) Drying ({g_model}): {dry_g:.1f} min with the gauge, {dry_r:.1f} min with
the radar: **{'agree' if dry_agree else 'DO NOT agree'}** (tolerance {GR['dry_agree_tol']:.0%}). The wetting speed and the peaks are read from the GAUGE run.

![gallery](../figures/{P.SID2_FIG_GALLERY.name})

## 7. The slow part: a water store or a drift?{tag['slow']}

{md_table(slow_tab.round(4))}

(`{rel(P.SID2_SLOW)}`.) The ladder again when every model may keep one linear drift per held-out link-event (prior sd s_b of M3D):

{md_table(da_tab[da_show].round(4))}

(`{rel(P.SID2_LADDER_DRIFT)}`.)

{slow_text}

![after the last shower](../figures/{P.SID2_FIG_AFTER.name})

Median (A - c) / peak after the last shower (`{rel(P.SID2_AFTER)}`):

{md_table(curve.iloc[[0, 1, 2, 3, 4, 6, 9, 12, 18, 24, 30, 35]].round(3))}

## 8. M5 and its name{tag['m5']}

The surfaces of M5 on S<3 (radar, lambda {lam_star:g}), in dB per minute (`{rel(P.SID2_NAMING)}`):

{md_table(pd.DataFrame([{'g: R (mm/h)': r_, **{f'W = {w_:g} dB': round(par5[f'g({r_:g},{w_:g})'], 4) for w_ in wk}} for r_ in rk[1:]]))}

{md_table(surf[surf.surface == 'h'].drop(columns=['surface', 'R_mmh']).round(4).rename(columns={'M5_db_per_min': 'h(W) of M5, dB/min', 'M5n_db_per_min': 'h(W) of M5n, dB/min'}))}

Numbers derived from the surfaces: t_wet {par5['t_wet']:.1f} min, t_dry {par5['t_dry']:.1f} min, level reached {', '.join(f"{par5[f'Weq[{BANDS[0]}]({r_:g})']:.2f} dB at {r_:g} mm/h" for r_ in KEY_R)} ({BANDS[0]}).

**M5 {'qualifies for naming' if m5_qualifies else 'does not qualify for naming'}** (rule: M5 is the chosen model, or its paired difference to the best holds 0).
{naming.get('text', 'No formula was named.' if not m5_qualifies else '(the naming is to be run)')}

![M5](../figures/{P.SID2_FIG_M5.name})

## 9. The current WAA equation

{box}

{md_table(key_par)}

The curves of {curve_model} on the main set ({SET_LABEL[main_set]}):

{md_table(curve_par)}

(`{rel(P.SID2_PARAMS)}`.)

## 10. In plain words

{plain}

## 11. Unclear / possibly wrong

{unclear}
"""
    save(report, P.SID2_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {(time.time() - clock) / 60:.1f} min")


def texts(q):
    """Sections 7 (the reading of the slow part), 9 (the equation box), 10 (plain words) and 11 (unclear) of the report,
    written from the numbers of the run (q = the local variables of main)."""
    params, chosen, main_set, sv, curve = q["params"], q["chosen"], q["main_set"], q["slow_real"], q["curve"]
    t3, naming, gt = q["ladders"]["S3"].set_index("model"), q["naming"], q["gt"]
    bands = q["BANDS"]

    def par(model, name, set_):
        r = params[(params.set == SET_LABEL[set_]) & (params.model == model) & (params.parameter == name)]
        return (np.nan, np.nan, np.nan) if r.empty else tuple(float(v) for v in r[["value", "lo", "hi"]].iloc[0])

    def show(model, name, set_, nd=1):
        v, lo, hi = par(model, name, set_)
        return f"{v:.{nd}f} [{lo:.{nd}f}, {hi:.{nd}f}]"

    def efold(t0, t1):
        c = curve.set_index("minutes_after_last_wet_bin").observed
        return (t1 - t0) / np.log(c[t0] / c[t1])

    obs = curve.set_index("minutes_after_last_wet_bin")
    integ = t3.integral_pred_over_obs

    # ---------------------------------------------------------------- section 7: the reading
    slow_text = f"""**Reading: {sv['verdict']}.** When every model may keep a free linear drift per held-out link-event, M4 is not
better than M3D (RMSE(M3D) - RMSE(M4) = {sv['M3D_minus_M4_drift_allowed']:+.4f} dB, interval {sv['drift_allowed_lo']:+.4f} .. {sv['drift_allowed_hi']:+.4f}):
a slow store driven by the rain explains no more than a straight ramp. But the ramps of M3D are not centred on 0: their
mean is {sv['mean_slope_db_h']:+.4f} dB/h (interval {sv['slope_lo']:+.4f} .. {sv['slope_hi']:+.4f}), so the level rises slowly through an event, and the slopes
hardly follow the rain depth of the link-event (Spearman {sv['spearman_slope_rain_depth']:+.2f}). From the rain alone M4 and M3D are equal within
the interval ({sv['M3D_minus_M4_from_R_alone']:+.4f} dB, {sv['from_R_lo']:+.4f} .. {sv['from_R_hi']:+.4f}). On fake data the same rule told a drift from a slow
store in 3 of 3 replicates each way (check Tb), so the rule works; the real data simply are neither of the two clean cases.

What the slow part looks like (the curve below, median over {int(obs.link_events.iloc[0])} link-events): the observed (A - c) / peak falls from
{obs.observed[30]:.3f} at 30 min to {obs.observed[90]:.3f} at 90 min and {obs.observed[175]:.3f} at 175 min after the last wet bin, an e-folding time of about
{efold(30, 90):.0f} min, then {efold(90, 175):.0f} min. M3 and M3D reach 0 within about an hour (integral over 0-180 min after the rain {integ['M3']:.2f} x and
{integ['M3D']:.2f} x the observed); the slow store of M4 (tau_s {sv['tau_s_M4_min']:.0f} min on S<3, {par('M4', 'tau_s', 'S1')[0]:.0f} min with the upper bound 600 inside the
interval on S<1) stays too high ({integ['M4']:.2f} x), and so does M5 ({integ['M5']:.2f} x). **So there is a slow part, small (2-4 % of the peak after
30-60 min) with a time constant of about 1-2 h: longer than M3 allows, much shorter than the slow store of M4.** Whether it
is water (a store) or a level change (a ramp) is not decided by these data."""

    # ---------------------------------------------------------------- section 9: the box
    m = chosen if chosen != "M1" else "M3"
    ms = SET_LABEL[main_set]
    if m == "M4":
        lines = [f"A(t) = W_f(t) + W_s(t) + L * P_band(R(t - d)) + c          model M4, {ms} (main set by the length rule), 90 % intervals",
                 "u(t) = a_link * F_band(R(t - d))                                the level the rain drives toward",
                 "W_f -> phi u       with tau_f_wet while W_f rises, tau_f_dry while it falls    (exact steps per minute)",
                 "W_s -> (1 - phi) u with tau_s",
                 f"d = {par('M4', 'd', main_set)[0]:+.0f} min;  phi = {show('M4', 'phi', main_set, 2)}",
                 f"tau_f_wet = {show('M4', 'tau_f_wet', main_set)} min with the radar  |  {show('M4', 'tau_f_wet', 'GG')} min with the GAUGE input (use this one)",
                 f"tau_f_dry = {show('M4', 'tau_f_dry', main_set)} min with the radar  |  {show('M4', 'tau_f_dry', 'GG')} min with the GAUGE input",
                 f"tau_s     = {show('M4', 'tau_s', main_set, 0)} min: NOT identified (the upper bound is 600); too high after the rain (section 7)"]
    else:
        lines = [f"A(t) = W(t) + L * P_band(R(t - d)) + c, dW/dt = (a_link F(R) - W) / tau, model {m}, {ms}",
                 f"tau_wet = {show(m, 'tau_wet', main_set)} min (radar) | {show(m, 'tau_wet', 'GG')} min (GAUGE)",
                 f"tau_dry = {show(m, 'tau_dry', main_set)} min (radar) | {show(m, 'tau_dry', 'GG')} min (GAUGE)"]
    for b in bands:
        if np.isfinite(par(m, f"F[{b}](2)", main_set)[0]):
            lines.append(f"F[{b}](0.5, 2, 5 mm/h) = " + ", ".join(show(m, f"F[{b}]({r:g})", main_set, 2) for r in (0.5, 2, 5)) + " dB")
    lines.append(f"a_link: 10-90 % range {par(m, 'a_p10', main_set)[0]:.2f} - {par(m, 'a_p90', main_set)[0]:.2f}")
    lines.append("P_band(R): free, never falling; " + ", ".join(f"P[{b}](10) = {show(m, f'P[{b}](10)', main_set, 2)} dB/km" for b in bands
                                                       if np.isfinite(par(m, f'P[{b}](10)', main_set)[0])))
    lines.append("P at R <= 2 mm/h and the split between F and P are NOT trusted (an optimiser trap, section 11)")
    box = "```\n" + "\n".join(lines) + "\n```"

    # ---------------------------------------------------------------- section 10: plain words
    g4 = gt.loc[("M4", "GAUGE")] if ("M4", "GAUGE") in gt.index else None
    rmse = t3.rmse
    plain = f"""- **The models are close.** On the links below 3 km the held-out RMSE is {rmse['M1']:.3f} dB without memory (M1), {rmse['M3']:.3f} with one
  store (M3), {rmse['M4']:.3f} with two stores (M4) and {rmse['M5']:.3f} with the free surfaces (M5). Memory clearly helps (M1 -> M3); every
  step after that is a few thousandths of a dB. The stop rule takes {chosen} on S<3 ({', '.join(f"{c} on {SET_LABEL[k]}" for k, c in q['chosen_of'].items() if k != 'S3')}). The WAA law is
  therefore known in its rough shape (a store that fills with rain and empties after it), not in its details.
- **Wetting is fast.** With the 1 min City gauge as input the antenna wets with a time constant of about
  {gt.wetting_min[('M4', 'GAUGE')] if g4 is not None else np.nan:.1f} min (M4; M3 {gt.wetting_min[('M3', 'GAUGE')]:.1f} min). With the radar it looks like 15-50 min, and the fake check Ta shows why:
  a 5 min average of the rain with a 30 % error makes the fitted wetting {q['med'].tau_wet['smooth'] / FK['truths']['Ta']['tau_wet']:.1f} x slower than the truth. The slow wetting of S1
  was the radar, not the antenna.
- **Drying of the fast part takes about 5-13 min** (M4: {par('M4', 'tau_f_dry', 'S1')[0]:.1f} min on S<1, {par('M4', 'tau_f_dry', 'S3')[0]:.1f} min on S<3, {gt.drying_min[('M4', 'GAUGE')]:.1f} min
  with the gauge). Gauge and radar disagree by more than the 30 % allowed, so the drying time is not yet one number.
- **Saturation.** With the radar the level F reaches about 3 dB by 2 mm/h and stays there up to 10 mm/h (times a_link,
  0.8-1.2 between links). With the gauge the level keeps rising up to 5 mm/h (F(0.5) {gt['F[28 GHz](0.5)'][('M4', 'GAUGE')]:.2f}, F(5) {gt['F[28 GHz](5)'][('M4', 'GAUGE')]:.2f} dB
  at 28 GHz): how fast the level saturates depends on the rain input.
- **Frequency.** 28, 32 and 38 GHz give the same F within the intervals (38 GHz a little higher); one law for all bands
  is good enough for now.
- **Length.** The numbers of M3 are the same on S<1, S<2 and S<3; for M4 the fast drying is shorter on S<3
  ({par('M4', 'tau_f_dry', 'S3')[0]:.1f} against {par('M4', 'tau_f_dry', 'S1')[0]:.1f} min), so by the rule S<1 stays the main result. The held-out residuals show no
  trend with the length. This is "not validated": the fake check Tc failed (section 11).
- **After the rain** a small slow part remains (section 7): about 1-2 h, not caught by any rung.
- **Peaks** are too low in every model (median {t3.shower_peak_error_median[chosen]:+.2f} dB for {chosen} with the radar, {gt.shower_peak_error_median[('M4', 'GAUGE')]:+.2f} dB with the gauge).
- **M5** (the free surfaces) is the best model on S<3 by a hair, but it did not pass its fake check (Td), its fits do not
  agree between starts, and {'its named formula ' + ('counts' if naming.get('counts') else 'does not count') if naming.get('done') else 'no formula was named'}; it is not a new law yet.

**What S3 must check:** (1) the fitting of P and F: a parametrisation without the trap of section 11, the fake check Tc
again, then the length rule again; (2) a slow part with a time constant of 1-2 h (a bounded slow store, or a store plus a
ramp), judged on the after-rain integral as well as the RMSE; (3) the gauge-input numbers as the wetting law, and why the
drying differs between gauge and radar; (4) why every model misses the peaks; (5) the M5 fit (more iterations, more
starts) before any formula is read from it. VAL is still untouched."""

    # ---------------------------------------------------------------- section 11: unclear
    fits = q["fits"]
    capped = fits[fits.folds_at_iteration_cap == fits.folds]
    nm = ""
    if naming.get("done"):
        nm = f"""
- **The naming.** No library formula describes the g surface well (the best relative weighted error over the lambda grid is
  {naming['g_best_error']:.2f}). The picked h has terms of both signs: as a continuous formula it is below 0 between 0 and about
  0.5 dB (at the knots, which M5n uses, it is not). The formula {'makes W run away in rain on the links with a large a_link (' + f"{naming['w_run']:.0f}" + ' dB after 2 h of 2 mm/h), so M5n could not even be refitted; ' if naming['stalled'] else ''}M5n
  {'counts' if naming['counts'] else 'does not count'} by the rule. E-SINDy fits the knot values without any condition that keeps W bounded (h rising,
  g not growing with W faster than h); a naming that is a usable model needs such conditions. And M5 itself failed Td,
  so a name read from its surfaces would say little about the antenna anyway."""
    unclear = f"""- **An optimiser trap in P (found after the run, one extra check, not in the config).** In the fake check Tc the fitted
  P(2 mm/h) is about 0 although the fit gets the true rain. P and F are sums of exp(parameter) increments, and an
  increment that goes low (exp(-8)) has almost no gradient and stays there. On fake replicate 0 (S<1) a refit that starts
  from the true P reaches a LOWER loss (13741 against 13792) and recovers P(2) = 0.42 (truth 0.38). So the real P at
  R <= 2 mm/h (0.002-0.3 dB/km) and therefore part of F are probably shaped by the optimiser, not only by the data. This is
  why Tc failed; the length rule and the F/P split are not validated.
- **Fake check Td failed:** M5 fitted on M3 data predicts as well as M3, but its derived numbers are wrong (levels
  0.36-0.52 x the truth, t_wet x1.56): the surfaces of M5 are not identified by these data. The real M5 surfaces
  (g(R, 0) about 0.005 dB/min, t_wet {q['par5']['t_wet']:.0f} min, the smallest smoothness lambda of the grid chosen) must not be read as physics.
- **Fits stop at the iteration caps.** The leave-one-event-out refits are warm starts capped at {M['maxiter']['fold']} iterations;
  {len(capped)} of {len(fits)} model-set fits reach the cap in every fold. The spread of the loss over the generic starts is large
  for M4, M3D and M5 (up to {fits.spread_over_starts_at_d0.max():.0f} loss units for M5): the fits on all events may not be the global
  minimum. Differences of a few thousandths of a dB between M3, M3D, M4 and M5 can change with longer fits.
- **tau_s of M4 is not identified** (S<1: {show('M4', 'tau_s', 'S1', 0)} min, the bound is 600).
- **The radar wetting time of M3** has a very wide interval ({show('M3', 'tau_wet', 'S1')} min on S<1): with the radar the
  wetting speed is not measured, only bounded from below.
- **Gauge input:** the gauge is up to 1 km from the path and measures at one point; the fitted lag is d = -3 min (A moves
  3 min before the gauge), and the drying times of gauge and radar runs do not agree.
- **The drift-allowed judge** uses the held-out A to fit one slope per link-event; it serves only the slow-part question.
- **Run time:** the computer slept for about 11.6 h during stage 2; the 'seconds' of real_S1_M5, real_S2_M5 and real_S3_M5D
  in the fits table include that sleep (each about 12 h of wall clock, about 30 min of computing).{nm}"""
    return box, slow_text, plain, unclear


if __name__ == "__main__":
    main()
