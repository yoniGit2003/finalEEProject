"""Step 3c: fake runs of the A7 procedure (decisions.md, Deviation 14 and step 3c, written before this code).

    python -m research.waa_ksurface.s3c_a7_variant              every run in order (a finished run is read from the cache)
    python -m research.waa_ksurface.s3c_a7_variant --report     write the report from the cache
    python -m research.waa_ksurface.s3c_a7_variant --one g_088  one run in this process (used by the main run)

The A7 procedure = step 3b with one change (Deviation 14): the drying fit uses the 14 learning events only (S1 "current"
variant). The tail flag is computed on the fake data before any fit; the link-events of the flagged events are left out.
The 37 added events take no part: not in stage 1, not in stage 2, not in the scored set L-dryhigh. Drying K0 against K1
on L-dryhigh with the margin of the current variant in S1 (0.00291 dB).

    (a) link spread    truth (g) (constant k, K0) with a level shift per link of sd 0.088 dB (measured, E4) and 0.15 dB.
                       PASS = drying K0 and 1 / k_dry within 25 % at both ends of k_wet. A false K1 at 0.088 is a STOP.
    (b) detectability  k_dry(W) at the ratios 1.5, 2, 3, five replicates each (no early stop); detected = K1 in 4 of 5.
    (c) flag coverage  how many of the 14 learning events the flag can check (event timing only, no attenuation).

One child process per run, one at a time; keep-awake and the memory rules of step 3b. Fake data only. No real
attenuation is read. VAL and TEST are not opened.
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
from research.waa_ksurface import a6_recovery as a6, data, fake6, ladder, s1_separability as s1, s3_recovery as s3, s3b_recovery as s3b
from research.waa_ksurface import paths as P
from research.wet_antenna.system_id import PADM

A6, S1, C = a6.A6, s1.S1, data.branch_config()["s3c"]
SPREAD = {"g_088": 0.088, "g_150": 0.15}                        # (a): run name -> level shift per link, sd (dB)
STOP_RUN = "g_088"                                              # a false K1 here is a STOP
GB = s3b.GB
log = a6.log


# =================================================================================================== names, seeds, margin
def margin():
    """The margin of the current variant (Deviation 14): the largest (a) gain on the adopted loss in S1, current variant."""
    t = pd.read_csv(P.S1_TABLE)
    r = t[(t.loss == C["loss"]) & (t.variant == C["variant"])].iloc[0]
    assert bool(r.separable)
    return max(float(r.a_gain_max_db), 0.0)


def det_names(a):
    return [a6.det_name("cw", a, k) for k in range(C["replicates"])]


def spec_seed(name):
    """(truth, seed) of a run (decisions.md, step 3c)."""
    if name in SPREAD:
        return {**A6["truths"]["g"], "link_sd_db": SPREAD[name]}, C["seed"] + list(SPREAD).index(name)
    _, a, rep = name.split("_")
    a = float(a)
    spec = {"kind": "dry_w", "k_wet": 0.2, "k_dry0": (1 / 15) / (1 + a), "w_scale_db": 1 / a, "size_a": a}     # as a6.spec_of
    return spec, C["seed"] + 10 + 10 * C["detect_sizes"].index(a) + int(rep[1:])


def truth_of(name):
    return "g" if name in SPREAD else "cw"


# =================================================================================================== one run
def run_one(name, S):
    """One run: fake data, the flag, the two stages on the unflagged learning events only, judged as in step 3b."""
    spec, seed = spec_seed(name)
    D, le, learn = S["D"], S["le"], S["learn"]
    t0, c0 = time.time(), time.process_time()
    sc = np.repeat(np.where(learn, 1.0, S1["radar_scale_added"]), D["le_n"] + 2 * PADM)     # the fake world as in step 3b
    Dt = dict(D)
    Dt["rain"] = np.asarray(D["rain"], float) / sc
    Dk, ex = fake6.make(Dt, le, S["n_cells"], spec, S["curves"], S["iF"], S["iP"], seed, S["step"], S["dist"],
                        ghost_pool=np.sort(le.event_id[learn].unique()))
    Dk["rain"] = Dk["rain"] * sc
    bi = ex["band_idx"]

    fl = a6.tail_flags(S, Dk)                                   # the flag, before any fit (learning events)
    flagged = set(fl.event_id[fl.flagged])
    use = (learn & ~le.event_id.isin(flagged).to_numpy()).astype(float)     # Deviation 14: learning events only
    used = (use > 0)[S["bin_le"]]

    counted, dry, high, _ = s1.score_sets(S, Dk)
    stage2 = {"w5": counted & used, "score": high & used, "margin_db": margin(), "selector": s3.select_k0_k1, "rungs": s3.DRY_RUNGS}
    E = A6["extension"]
    start_F = [np.full(len(i), E["f_start_inc_db"]) for i in S["iF"]]
    start_P = [E["p_start_share"] * i for i in S["iP"]]
    res = ladder.two_stage(Dk, le, bi, S["curves"], S["gauge"], inc_F=start_F, inc_P=start_P, hold_ext=False, use=use, log=log, tag=f"({name}) ",
                           stage2=stage2, wet_rungs=s3.WET_RUNGS)
    dens = a6.state_density(S, Dk, ex)
    out = {"name": name, "truth": truth_of(name), "spec": spec, "seed": seed, "flags": fl,
           "events_flagged": int(fl.flagged.sum()), "events_checked": int(len(fl)), "false_removals": int(fl.flagged.sum()),
           "scored_values": int((high & used).sum()), "scored_link_events": int(len(np.unique(S["bin_le"][high & used])))}
    out["cv"] = pd.concat([res["wet"]["table"], res["dry"]["table"]], ignore_index=True)
    out["tau_wet"] = res.get("tau_wet")
    out["judge"], _ = a6.judge(S, spec, res, dens, True)
    out["minutes"] = (time.time() - t0) / 60
    out["cpu_minutes"] = (time.process_time() - c0) / 60
    return out


# =================================================================================================== verdicts
def verdict(name, r):
    d = s3b.drying(r)
    k_ok = bool(d.k_ok) if d.k_ok is not None else False
    if name in SPREAD:
        return "PASS" if (d.selected_rung == "K0" and k_ok) else "FAIL"
    return "K1 returned" if d.selected_rung == "K1" else "K0 returned"


def flag_coverage(S):
    """(c): per learning event, the short-link link-events with the 20 bins after the end of rain that the flag needs
    (the rule of `a6_recovery.tail_flags`, from the event timing only: no attenuation is read)."""
    D, le, F = S["D"], S["le"], A6["flag"]
    rows = []
    for e in np.nonzero(S["short"])[0]:
        n, end = D["le_n"][e], (D["le_last"][e] + 1) * 5
        rows.append({"event_id": le.event_id.iloc[e], "ok": not (end + 20 > n or end < 5)})
    t = pd.DataFrame(rows).groupby("event_id").agg(short_link_events=("ok", "size"), with_tail=("ok", "sum")).reset_index()
    ev = pd.DataFrame({"event_id": np.sort(le.event_id[S["learn"]].unique())}).merge(t, how="left").fillna(0)
    ev["checkable"] = ev.with_tail >= F["min_link_events"]
    return ev


# =================================================================================================== report
def report(R, cov, state, stopped):
    rows = []
    for n, r in R.items():
        d = s3b.drying(r)
        g, lo, hi = s3.k1_gain(r["cv"])
        rows.append({"run": n, "seed": r["seed"], "true rung": d.true_rung, "drying rung": d.selected_rung, "K1 gain on L-dryhigh (dB)": round(g, 5),
                     "interval": f"{lo:+.5f} to {hi:+.5f}", "recovery (k_wet gauge / radar)": s3b.recovery_text(r),
                     "events flagged": f"{r['events_flagged']} of {r['events_checked']}", "scored values": r["scored_values"],
                     "scored link-events": r["scored_link_events"], "verdict": verdict(n, r), "minutes": round(r["minutes"], 1),
                     "cpu minutes": round(r["cpu_minutes"], 1), "peak memory (GB)": round(r.get("peak_gb", np.nan), 2)})
    runs = pd.DataFrame(rows)
    save(runs, P.S3C_RUNS, index=False)
    if R:
        save(pd.concat([r["cv"].assign(run=n) for n, r in R.items()], ignore_index=True), P.S3C_CV, index=False)
        save(pd.concat([r["flags"].assign(run=n) for n, r in R.items()], ignore_index=True), P.S3C_FLAGS, index=False)
    save(cov, P.S3C_COVERAGE, index=False)
    m = margin()
    L = ["# Step 3c: fake runs of the A7 procedure (Deviation 14: learning events only)\n",
         "Written by `research/waa_ksurface/s3c_a7_variant.py` (libraries `s3b_recovery.py`, `s3_recovery.py`, `a6_recovery.py`, `s1_separability.py`, "
         "`fake6.py`, `ladder.py`). Fake attenuation as in step 3b, fresh seeds (base "
         f"{C['seed']}). **The fit uses the 14 learning events only** (S1 current variant): the tail flag before any fit, the link-events of flagged events "
         f"left out; the 37 added events take no part. Drying K0 against K1 on L-dryhigh, margin {m:.5f} dB (current variant, S1) with a paired 90 % interval "
         "that excludes zero; wetting K0, not judged (R3). Where two numbers are given as `x / y`, the first is with k_wet at its gauge-input estimate, the "
         "second at its radar-input estimate. A returned drying K1 is named: *" + s3b.K1_NAME + "* (R5). No real attenuation is read. VAL and TEST were not opened.\n"]
    sp = runs[runs.run.isin(list(SPREAD))]
    if len(sp):
        L += ["## (a) The link spread (truth K0, 1 / k_dry = 15 min)\n",
              "Every fake truth of A6, S1, step 3 and step 3b already carried a level shift per link of sd 0.088 dB, the measured spread (E4), as part "
              "of the default noise; so `g_088` is truth (a) on a fresh seed under the A7 procedure, and every truth that passed so far passed at the "
              "measured link spread. It is run first because its STOP rule is new.\n",
              md_table(sp.assign(**{"link sd (dB)": sp.run.map(SPREAD)})[["run", "link sd (dB)", "seed", "drying rung", "K1 gain on L-dryhigh (dB)", "interval",
                                                                        "recovery (k_wet gauge / radar)", "events flagged", "scored values", "verdict"]]), ""]
    if stopped:
        L.append(f"\n**STOP by the rule: {stopped}.** The per-link question is reopened by measurement; nothing after it was run.\n")
    det_rows = []
    for a in C["detect_sizes"]:
        names = [n for n in det_names(a) if n in R]
        if names:
            hits = sum(s3b.drying(R[n]).selected_rung == "K1" for n in names)
            det_rows.append({"size_a": a, "ratio": 1 + 2 * a, "replicates_run": len(names), "hits (K1 returned)": hits,
                             "detected": bool(len(names) == C["replicates"] and hits >= C["replicates"] - 1), "runs": ", ".join(names)})
    if det_rows:
        det = pd.DataFrame(det_rows)
        save(det, P.S3C_DETECT, index=False)
        L += ["\n## (b) Detectability of k_dry(W) for the A7 variant\n",
              "k_dry = (1/15)·(1 + a·W) / (1 + a), k_wet = 1/5; ratio k(W = 2 dB) / k(W = 0) = 1 + 2a. Five replicates per ratio, no early stop. "
              "Detected = drying K1 returned in at least 4 of 5.\n", md_table(det), ""]
        d_ok = det[det.detected]
        L.append(f"**Smallest ratio detected in at least 4 of 5: {d_ok.ratio.min():g}.**\n" if len(d_ok) else "**No ratio detected in at least 4 of 5.**\n")
        L.append("Why step 3b showed 4 of 4 and not 5: there the replicates of a ratio stopped as soon as the verdict was decided (4 hits, or 2 misses; "
                 "Deviation 13), so after 4 hits the fifth replicate was never run. Here every ratio runs all 5.\n")
    L += ["\n## (c) The flag on the 14 learning events\n",
          f"Checkable = at least {A6['flag']['min_link_events']} short-link link-events with the 20 bins after the end of rain (event timing only; no "
          f"attenuation read). **{int(cov.checkable.sum())} of {len(cov)} learning events are checkable.** Removed on real data: not counted here (it needs "
          "the real attenuation at the end of the rain, which is not read before A7; it is the first number of A7).\n",
          md_table(cov), ""]
    if len(runs):
        L += ["\n## Execution\n", f"Runs: {len(R)}; wall time {runs.minutes.sum() / 60:.1f} h, processor time {runs['cpu minutes'].sum() / 60:.1f} h; "
              f"largest peak memory {runs['peak memory (GB)'].max():.2f} GB.\n",
              md_table(runs[["run", "minutes", "cpu minutes", "peak memory (GB)"]]), ""]
    if state:
        L.append(f"Keep-awake request accepted by Windows: {state.get('awake')}; gaps of the memory watch above 3 min: {len(state['gaps'])} "
                 f"({', '.join(f'{g:.0f} min' for g in state['gaps']) or 'none'}); lowest free memory seen {state['min_free']:.2f} GB; low-memory events: {len(state['low'])}.\n")
        if state.get("incomplete"):
            L.append(f"\n**NOT COMPLETE: {state['incomplete']}. Run again to go on from the cache.**\n")
    save("\n".join(L) + "\n", P.S3C_REPORT)
    log(f"report: {rel(P.S3C_REPORT)}")


# =================================================================================================== main
def main():
    cache = P.CACHE / "s3c_a7_variant"
    cache.mkdir(parents=True, exist_ok=True)
    if "--one" in sys.argv:
        name = sys.argv[sys.argv.index("--one") + 1]
        r = run_one(name, s3.setting())
        r["peak_gb"] = psutil.Process().memory_info().peak_wset / GB
        (cache / f"{name}.pkl").write_bytes(pickle.dumps(r))
        log(f"=== ({name}): {verdict(name, r)}; drying {s3b.drying(r).selected_rung}; {s3b.recovery_text(r)}; flagged {r['events_flagged']} of "
            f"{r['events_checked']}; {r['minutes']:.1f} min wall, {r['cpu_minutes']:.1f} min processor; peak memory {r['peak_gb']:.2f} GB ===")
        return
    only_report = "--report" in sys.argv
    state = {"min_free": np.inf, "gaps": [], "low": [], "peak": 0.0, "incomplete": ""}
    if not only_report:
        free = psutil.virtual_memory().available / GB
        log(f"=== step 3c started {time.strftime('%Y-%m-%d %H:%M:%S')}; free memory {free:.1f} GB ===")
        if free < C["min_free_start_gb"]:
            log(f"=== not started: free memory {free:.1f} GB is below {C['min_free_start_gb']} GB ===")
            sys.exit(3)
        state["awake"] = s3b.keep_awake()
        log(f"=== keep-awake (ES_CONTINUOUS | ES_SYSTEM_REQUIRED): {'on' if state['awake'] else 'NOT accepted'} ===")
        state["min_free"] = free
        threading.Thread(target=s3b.watch_memory, args=(state,), daemon=True).start()
    cov = flag_coverage(s3.setting())
    R = {}

    def get(name):
        f = cache / f"{name}.pkl"
        if not f.exists():
            if only_report or state["incomplete"]:
                return None
            free = psutil.virtual_memory().available / GB
            need = C["min_free_gb"]                             # execution setting (decisions.md): no peak term; a run that dies resumes from the cache
            if free < need:
                state["incomplete"] = f"free memory {free:.1f} GB is below the {need:.1f} GB a run needs: stopped before ({name})"
                log("=== " + state["incomplete"] + " ===")
                return None
            log(f"=== ({name}) seed {spec_seed(name)[1]} (free memory {free:.1f} GB) {time.strftime('%H:%M:%S')} ===")
            before = s3b.programs()
            rc = subprocess.run([sys.executable, "-u", "-m", "research.waa_ksurface.s3c_a7_variant", "--one", name]).returncode
            if rc != 0:                                         # the child died (e.g. out of memory): say which programs grew, stop cleanly
                cur = s3b.programs()
                grew = sorted(((cur[n] - before.get(n, 0.0), n) for n in cur), reverse=True)[:5]
                state["incomplete"] = (f"({name}) died with exit code {rc}; free memory {psutil.virtual_memory().available / GB:.1f} GB; "
                                       "grew most during it: " + "; ".join(f"{n} {g:+.2f} GB" for g, n in grew))
                log("=== " + state["incomplete"] + " ===")
                return None
        r = pickle.loads(f.read_bytes())
        state["peak"] = max(state["peak"], r.get("peak_gb", 0.0))
        R[name] = r
        return r

    stopped = ""
    r = get(STOP_RUN)                                           # (a) at the measured spread first: a false K1 is a STOP
    if r is not None and s3b.drying(r).selected_rung == "K1":
        stopped = f"({STOP_RUN}) false K1 at the measured link spread 0.088 dB"
    if r is not None and not stopped and get("g_150") is not None:
        for a in C["detect_sizes"]:                             # (b) five replicates per ratio, no early stop
            if any(get(n) is None for n in det_names(a)):
                break
    report(R, cov, None if only_report else state, stopped)
    if not only_report:
        log(f"=== step 3c ended {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
    if stopped or state["incomplete"]:
        sys.exit(1 if stopped else 2)


if __name__ == "__main__":
    main()
