"""Step 3d: the A7 procedure with 1 / k_wet FIXED at 2.5, 5 and 10 min, on fake data (decisions.md, step 3d, written
before this code).

    python -m research.waa_ksurface.s3d_kwet_fixed              every replicate in order (a finished one is read from the cache)
    python -m research.waa_ksurface.s3d_kwet_fixed --report     write the report from the cache
    python -m research.waa_ksurface.s3d_kwet_fixed --one a_r0   one replicate in this process (used by the main run)

The procedure is that of step 3c (`s3c_a7_variant.py`): the tail flag before any fit, the 14 learning events only, drying
K0 against K1 on L-dryhigh with the margin of the current variant (0.00291 dB), extension free. What changes: stage 1
(gauge input, wetting K0) runs once and gives the start of the extension; then stage 2 runs three times, each time with
the wetting rate HELD at 1 / k_wet = 2.5, 5, 10 min. `fit_fixed_kwet` is the function A7 uses on the real data.

Truths: (a) constant k (k_wet = 1/5, k_dry = 1/15); (c) k_dry = (1 + W) / 30. Two replicates each, seeds 3840-3843.
PASS = (a) drying K0 with 1 / k_dry within 25 % and (c) drying K1, at every value; a rung flip is a STOP.
Fake data only. VAL and TEST are not opened.
"""
import pickle
import subprocess
import sys
import threading
import time

import numpy as np
import pandas as pd
import psutil

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a6_recovery as a6, data, fake6, ladder, ksurf, s1_separability as s1, s3_recovery as s3, s3b_recovery as s3b
from research.waa_ksurface import s3c_a7_variant as s3c
from research.waa_ksurface import paths as P
from research.wet_antenna.system_id import PADM

A6, S1, D3 = a6.A6, s1.S1, data.branch_config()["s3d"]
GB = s3b.GB
log = a6.log


# =================================================================================================== the procedure with k_wet fixed (A7 uses it)
def fit_fixed_kwet(Dk, le, band_idx, curves, gauge, use, stage2, inc_F=None, inc_P=None, tag=""):
    """Stage 1 once (gauge input, wetting K0: the start of the extension), then stage 2 (drying K0 vs K1, the selection of
    step 3c) with 1 / k_wet held at every value of `tau_wet_min`. Returns {tau: {"wet", "dry", "model", "x"}}."""
    rain_g, ok5 = ladder.gauge_input(Dk)
    prob_g = ksurf.prepare(Dk, band_idx, curves, rain_g, 0)
    prob_g["w5"] = prob_g["w5"] * ok5
    st1 = ladder.stage(prob_g, le, use * gauge, "wetting", "K1", None, inc_F, inc_P, False, log=log, tag=tag, rungs=s3.WET_RUNGS)
    m1, x1 = st1["model"], st1["x"]
    gw, hw = x1[m1.sl["gw"]].copy(), x1[m1.sl["hw"]].copy()
    ext = m1.parts(x1)[4:]
    log(f"   {tag}stage 1: 1 / k_wet estimated {1 / m1.k(x1, 'wetting', 0, 0):.2f} min (not used: it is held below)")
    del prob_g, rain_g, ok5
    prob_r = ksurf.prepare(Dk, band_idx, curves, Dk["rain"], 0)
    prob_r["w5"] = np.asarray(stage2["w5"], float)
    out = {}
    for tau in D3["tau_wet_min"]:
        wet = (np.full_like(gw, np.log(1.0 / tau)), hw)          # the wetting rate held at 1 / tau
        st2 = ladder.stage(prob_r, le, use, "drying", st1["rung"], wet, [e.copy() for e in ext[0]], [e.copy() for e in ext[1]], False, log=log,
                           tag=f"{tag}[1/k_wet {tau:g}] ", score=np.asarray(stage2["score"], bool), margin_db=stage2["margin_db"],
                           selector=stage2["selector"], rungs=stage2["rungs"])
        out[tau] = {"wet": st1, "dry": st2, "model": st2["model"], "x": st2["x"]}
    return out


def k_dry_at(model, x, W):
    """1 / k_dry (min) of a fit at the levels W (dB), after rain (R = 0)."""
    return 1.0 / model.k(x, "drying", 0.0, np.asarray(W, float))


# =================================================================================================== one replicate
def spec_seed(name):
    t, rep = name.split("_r")
    spec = A6["truths"]["a"] if t == "a" else A6["truths"]["c"]
    return spec, D3["seed"] + 2 * ("ac".index(t)) + int(rep)


def names():
    return [f"{t}_r{k}" for t in ("a", "c") for k in range(D3["replicates"])]


def run_one(name, S):
    """One fake dataset, the flag, then the fit at the three values of 1 / k_wet; judged at each."""
    spec, seed = spec_seed(name)
    D, le, learn = S["D"], S["le"], S["learn"]
    t0, c0 = time.time(), time.process_time()
    sc = np.repeat(np.where(learn, 1.0, S1["radar_scale_added"]), D["le_n"] + 2 * PADM)     # the fake world of step 3c
    Dt = dict(D)
    Dt["rain"] = np.asarray(D["rain"], float) / sc
    Dk, ex = fake6.make(Dt, le, S["n_cells"], spec, S["curves"], S["iF"], S["iP"], seed, S["step"], S["dist"],
                        ghost_pool=np.sort(le.event_id[learn].unique()))
    Dk["rain"] = Dk["rain"] * sc
    fl = a6.tail_flags(S, Dk)
    use = (learn & ~le.event_id.isin(set(fl.event_id[fl.flagged])).to_numpy()).astype(float)
    used = (use > 0)[S["bin_le"]]
    counted, dry, high, _ = s1.score_sets(S, Dk)
    stage2 = {"w5": counted & used, "score": high & used, "margin_db": s3c.margin(), "selector": s3.select_k0_k1, "rungs": s3.DRY_RUNGS}
    E = A6["extension"]
    start_F = [np.full(len(i), E["f_start_inc_db"]) for i in S["iF"]]
    start_P = [E["p_start_share"] * i for i in S["iP"]]
    fits = fit_fixed_kwet(Dk, le, ex["band_idx"], S["curves"], S["gauge"].astype(float), use, stage2, start_F, start_P, tag=f"({name}) ")
    dens = a6.state_density(S, Dk, ex)
    rows, cvs = [], []
    for tau, res in fits.items():
        j, _ = a6.judge(S, spec, res, dens, True)
        d = j.set_index("regime").loc["drying"]
        g, lo, hi = s3.k1_gain(res["dry"]["table"])
        rows.append({"replicate": name, "truth": name[0], "seed": seed, "1/k_wet (min)": tau, "true rung": d.true_rung, "drying rung": d.selected_rung,
                     "K1 gain (dB)": g, "gain lo": lo, "gain hi": hi, "1/k_dry fit (min)": d.tau_fit_min, "1/k_dry truth (min)": d.tau_true_min,
                     "share within 25 %": d.share_within_25pct, "k_ok": bool(d.k_ok) if d.k_ok is not None else None,
                     "1/k_dry at W = 0.5, 1, 2 dB": ", ".join(f"{v:.1f}" for v in k_dry_at(res["model"], res["x"], [0.5, 1.0, 2.0]))})
        cvs.append(res["dry"]["table"].assign(replicate=name, tau_wet=tau))
    return {"name": name, "rows": pd.DataFrame(rows), "cv": pd.concat(cvs, ignore_index=True), "flags": fl,
            "events_flagged": int(fl.flagged.sum()), "events_checked": int(len(fl)),
            "minutes": (time.time() - t0) / 60, "cpu_minutes": (time.process_time() - c0) / 60}


# =================================================================================================== verdicts
def row_verdict(r):
    if r["truth"] == "a":
        if r["drying rung"] != "K0":
            return "FLIP"
        return "ok" if r["k_ok"] else "FAIL (1/k_dry outside 25 %)"
    return "ok" if r["drying rung"] == "K1" else "FLIP"


def verdict(rows):
    v = rows.apply(row_verdict, axis=1)
    if (v == "FLIP").any():
        return "STOP (rung flip)"
    return "PASS" if (v == "ok").all() else "FAIL"


# =================================================================================================== report
def report(R, state, stopped):
    rows = pd.concat([r["rows"] for r in R.values()], ignore_index=True) if R else pd.DataFrame()
    if len(rows):
        rows["verdict"] = rows.apply(row_verdict, axis=1)
        save(rows, P.S3D_RUNS, index=False)
        save(pd.concat([r["cv"] for r in R.values()], ignore_index=True), P.S3D_CV, index=False)
    L = ["# Step 3d: the A7 procedure with 1 / k_wet fixed at 2.5, 5, 10 min (fake data)\n",
         "Written by `research/waa_ksurface/s3d_kwet_fixed.py` (procedure of `s3c_a7_variant.py`). Fake attenuation as in step 3c, fresh seeds "
         f"(base {D3['seed']}). The tail flag before any fit; the 14 learning events only; stage 1 (gauge input, wetting K0) once, then stage 2 "
         f"(drying K0 against K1 on L-dryhigh, margin {s3c.margin():.5f} dB, paired 90 % interval excluding zero) with the wetting rate HELD at "
         "1 / k_wet = 2.5, 5, 10 min. Truth (a): constant k, 1 / k_dry = 15 min; truth (c): k_dry = (1 + W) / 30 (ratio 3). PASS = (a) K0 with "
         "1 / k_dry within 25 % and (c) K1 at every value; a rung flip is a STOP. No real attenuation is read. VAL and TEST were not opened.\n"]
    if len(rows):
        show = rows.assign(**{"K1 gain (dB)": rows["K1 gain (dB)"].round(5), "1/k_dry fit (min)": rows["1/k_dry fit (min)"].round(1),
                              "share within 25 %": rows["share within 25 %"].round(2),
                              "interval": [f"{a:+.5f} to {b:+.5f}" for a, b in zip(rows["gain lo"], rows["gain hi"])]})
        L += ["## Every fit\n", md_table(show[["replicate", "seed", "1/k_wet (min)", "true rung", "drying rung", "K1 gain (dB)", "interval",
                                               "1/k_dry fit (min)", "1/k_dry truth (min)", "share within 25 %", "1/k_dry at W = 0.5, 1, 2 dB", "verdict"]]), ""]
        done = len(R) == len(names())
        L.append(f"**Verdict: {verdict(rows)}**" + ("" if done else f" (only {len(R)} of {len(names())} replicates)") + ".\n")
    if stopped:
        L.append(f"\n**STOP by the rule: {stopped}.** The 'rung not determined' reading would be expected on real data; the review decides before A7.\n")
    if R:
        ex = pd.DataFrame([{"replicate": n, "events flagged": f"{r['events_flagged']} of {r['events_checked']}", "minutes": round(r["minutes"], 1),
                            "cpu minutes": round(r["cpu_minutes"], 1), "peak memory (GB)": round(r.get("peak_gb", np.nan), 2)} for n, r in R.items()])
        L += ["\n## Execution\n", md_table(ex), ""]
    if state:
        L.append(f"Keep-awake request accepted by Windows: {state.get('awake')}; gaps of the memory watch above 3 min: {len(state['gaps'])}; "
                 f"lowest free memory seen {state['min_free']:.2f} GB; low-memory events: {len(state['low'])}.\n")
        if state.get("incomplete"):
            L.append(f"\n**NOT COMPLETE: {state['incomplete']}. Run again to go on from the cache.**\n")
    save("\n".join(L) + "\n", P.S3D_REPORT)
    log(f"report: {rel(P.S3D_REPORT)}")


# =================================================================================================== main
def run_children(module, cache, items, state, label, run_rule=None):
    """Run every item in its own child process (cached); stop cleanly below the memory floor or if a child dies.
    run_rule(name, result) -> a STOP text or "". Returns ({name: result}, stop text)."""
    R = {}
    for name in items:
        f = cache / f"{name}.pkl"
        if not f.exists():
            free = psutil.virtual_memory().available / GB
            if free < D3["min_free_gb"]:
                state["incomplete"] = f"free memory {free:.1f} GB is below {D3['min_free_gb']} GB: stopped before ({name})"
                log("=== " + state["incomplete"] + " ===")
                break
            log(f"=== ({name}) {label(name)} (free memory {free:.1f} GB) {time.strftime('%H:%M:%S')} ===")
            before = s3b.programs()
            rc = subprocess.run([sys.executable, "-u", "-m", module, "--one", name]).returncode
            if rc != 0:
                cur = s3b.programs()
                grew = sorted(((cur[n] - before.get(n, 0.0), n) for n in cur), reverse=True)[:5]
                state["incomplete"] = (f"({name}) died with exit code {rc}; free memory {psutil.virtual_memory().available / GB:.1f} GB; "
                                       "grew most during it: " + "; ".join(f"{n} {g:+.2f} GB" for g, n in grew))
                log("=== " + state["incomplete"] + " ===")
                break
        R[name] = pickle.loads(f.read_bytes())
        stop = run_rule(name, R[name]) if run_rule else ""
        if stop:
            return R, stop
    return R, ""


def start_state(what):
    state = {"min_free": np.inf, "gaps": [], "low": [], "peak": 0.0, "incomplete": ""}
    free = psutil.virtual_memory().available / GB
    log(f"=== {what} started {time.strftime('%Y-%m-%d %H:%M:%S')}; free memory {free:.1f} GB ===")
    if free < D3["min_free_gb"]:
        log(f"=== not started: free memory {free:.1f} GB is below {D3['min_free_gb']} GB ===")
        sys.exit(3)
    state["awake"] = s3b.keep_awake()
    log(f"=== keep-awake (ES_CONTINUOUS | ES_SYSTEM_REQUIRED): {'on' if state['awake'] else 'NOT accepted'} ===")
    state["min_free"] = free
    threading.Thread(target=s3b.watch_memory, args=(state,), daemon=True).start()
    return state


def main():
    cache = P.CACHE / "s3d_kwet_fixed"
    cache.mkdir(parents=True, exist_ok=True)
    if "--one" in sys.argv:
        name = sys.argv[sys.argv.index("--one") + 1]
        r = run_one(name, s3.setting())
        r["peak_gb"] = psutil.Process().memory_info().peak_wset / GB
        (cache / f"{name}.pkl").write_bytes(pickle.dumps(r))
        for _, x in r["rows"].iterrows():
            log(f"=== ({name}) 1/k_wet {x['1/k_wet (min)']:g} min: drying {x['drying rung']}; 1/k_dry {x['1/k_dry fit (min)']:.1f} min "
                f"(truth {x['1/k_dry truth (min)']:.0f}); K1 gain {x['K1 gain (dB)']:+.5f} dB; {row_verdict(x)} ===")
        log(f"=== ({name}) done: {r['minutes']:.1f} min wall, {r['cpu_minutes']:.1f} min processor; peak memory {r['peak_gb']:.2f} GB ===")
        return
    if "--report" in sys.argv:
        R = {n: pickle.loads((cache / f"{n}.pkl").read_bytes()) for n in names() if (cache / f"{n}.pkl").exists()}
        report(R, None, "")
        return
    state = start_state("step 3d")

    def rule(name, r):
        return f"({name}) rung flip" if (r["rows"].apply(row_verdict, axis=1) == "FLIP").any() else ""

    R, stopped = run_children("research.waa_ksurface.s3d_kwet_fixed", cache, names(), state, lambda n: f"seed {spec_seed(n)[1]}", rule)
    report(R, state, stopped)
    log(f"=== step 3d ended {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
    if stopped or state["incomplete"]:
        sys.exit(1 if stopped else 2)


if __name__ == "__main__":
    main()
