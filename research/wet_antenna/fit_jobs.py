"""Fitting jobs that more than one step script runs. Each job fits ONE dataset and saves every finished stage in a
json file, so that an interrupted run continues where it stopped and a second run only reads the files.

A job is a dict (it travels to a worker process). Besides its own fields it carries "cache": the folder of the step
that runs it (research/.cache/wet_antenna/<script>/). Settings are read from research/config.yaml:
`wet_antenna.stress_pooled_fit` (the menu fits), `wet_antenna.check_free_curve` (how a free curve is fitted) and
`wet_antenna.confirm_free_curve` (window and seeds of the confirmation run).
"""
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from research.common.config import load_config
from research.wet_antenna import datasets
from research.wet_antenna import free_curve as wc
from research.wet_antenna import pooled_fit as wp

CFG = load_config()
S1C = CFG["wet_antenna"]["stress_pooled_fit"]              # menu fits: candidates, window, starts
FIT = CFG["wet_antenna"]["check_free_curve"]               # free curve: starts, penalties, budgets, prior
S3B = CFG["wet_antenna"]["confirm_free_curve"]             # confirmation run: window, seeds, F4 starts

_INPUTS = {}                                               # per process: noise and offset tables (read once)


def noise_tables():
    """(noise per sub-link, offset std per sub-link), read once per process."""
    if not _INPUTS:
        _, _INPUTS["noise"], _INPUTS["offset_sd"] = datasets.load_noise()
    return _INPUTS["noise"], _INPUTS["offset_sd"]


def save_json(file, result):
    """Write `result` to `file` through a temporary file, so that a crash never leaves a half-written json."""
    tmp = Path(file).with_suffix(".tmp")
    tmp.write_text(json.dumps(result))
    tmp.replace(file)


def menu_fits_job(job):
    """Fit every menu candidate from every start to ONE simulated dataset; each fit is cached as a json file.
    job: test, truth, rep, sim (parquet file), cache (folder), starts {candidate: {label: [theta0, s0, maxfev]}}
    (theta0 = None for the label `f1sol`: built from the best generic F1 fit of this dataset).
    Returns the list of fit rows. (Step 3b-1c; step 3b-3c reuses it for its menu fits.)"""
    noise, offset_sd = noise_tables()
    cache = Path(job["cache"])
    data = {}

    def dataset():
        if not data:                                       # only loaded if a fit is missing from the cache
            rows = pd.read_parquet(job["sim"])
            tails = [{"tail_id": t, "sublink_id": g.sublink_id.iloc[0], "cml_id": g.cml_id.iloc[0],
                      "event_copy": g.event_copy.iloc[0], "A_start": g.A_start_trace_db.iloc[0],
                      "t": g.minutes_since_t0_trace.to_numpy(), "a": g.A_db.to_numpy()} for t, g in rows.groupby("tail_id")]
            D = wp.pack(tails, noise, offset_sd, S1C["window_min"], S1C["min_points"])
            used = set(D["tail_id"])
            data["D"] = D                                  # N of the BIC: distinct (link, event draw) pairs
            data["n_ind"] = len({(tl["cml_id"], tl["event_copy"]) for tl in tails if tl["tail_id"] in used})
        return data["D"], data["n_ind"]

    out, f1_best = [], None
    for cand in S1C["candidates"]:
        for label, (theta0, s0, maxfev) in job["starts"][cand].items():
            f_ = cache / f"fit_{job['test']}_{job['truth']}_{job['rep']}_{cand}_{label}.json"
            if f_.exists():
                r = json.loads(f_.read_text())
            else:
                D, n_ind = dataset()
                if label == "f1sol":
                    theta0, s0 = wp.f1_start(cand, f1_best["theta"][0]), f1_best["s"]
                t0 = time.time()
                r = wp.fit(cand, D, theta0, s0, maxfev)
                r.update({"test": job["test"], "truth": job["truth"], "rep": job["rep"], "candidate": cand,
                          "start": label, "theta0": [float(v) for v in theta0], "s0": float(s0),
                          "n_tails": len(D["tail_id"]), "n_independent": n_ind, "k": len(wp.NAMES[cand]) + 1,
                          "seconds": round(time.time() - t0, 1)})
                save_json(f_, r)
            out.append(r)
        if cand == "F1":                                   # the best generic F1 fit seeds the `f1sol` starts
            gen = [r for r in out if r["candidate"] == "F1" and r["start"] != "calib"]
            f1_best = max(gen, key=lambda r: r["ll"]) if gen else None
    return out


def menu_fit(D, n_ind, cand, labels, res, save):
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


def free_curve_job(job):
    """Free-curve pooled fit of ONE dataset (first 60 min), as in the confirmation run (step 3b-3b; step 3b-4
    reuses it). The stages are those of step 3b-3 without the curvature interval; each is cached in one json file
    per job: <cache>/job_<tag>.json."""
    f_ = Path(job["cache"]) / f"job_{job['tag']}.json"
    res = json.loads(f_.read_text()) if f_.exists() else {}
    if res.get("done"):
        return res
    noise, offset_sd = noise_tables()

    def save():
        save_json(f_, res)

    t0 = time.time()
    tails = datasets.load_tails(job["sim"])
    D = wp.pack(tails, noise, offset_sd, S3B["window_min"], FIT["min_points"])
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
            x0 = np.r_[wc.start_theta(spec, cv.knots), np.log(FIT["s_start"])]
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
    evs = np.random.default_rng([S3B["seed_base"], job["rep"]]).permutation(sorted({tl["event_copy"] for tl in tails}))
    half = [set(evs[::2].tolist()), set(evs[1::2].tolist())]
    Dh = None
    for lam in FIT["lambda_grid"]:
        if str(lam) not in res["cv"]:
            if Dh is None:
                Dh = [wp.pack([tl for tl in tails if tl["event_copy"] in h], noise, offset_sd, S3B["window_min"],
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
        menu_fit(D, n_ind, "F4", S3B["menu_f4_starts"], res, save)
    res["seconds"] = round(res.get("seconds", 0.0) + time.time() - t0, 1)
    res["done"] = True
    save()
    return res
