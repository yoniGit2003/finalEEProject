"""Step 3b: the fake-data truths again with the tail flag in the procedure (Deviation 12), the history test B3'
calibrated (Deviation 13), then (g) and the detectability study of k_dry(W).

    python -m research.waa_ksurface.s3b_recovery                every run in order (a finished run is read from the cache)
    python -m research.waa_ksurface.s3b_recovery --report       write the report from the cache
    python -m research.waa_ksurface.s3b_recovery --one c_r2     one run in this process (used by the main run)

What differs from `s3_recovery.py` (decisions.md, Deviations 12 and 13, written before this code):

    the flag     the tail flag of A6 (`a6_recovery.tail_flags`, model-free, event-level, learning events only) is computed
                 on the fake data before any fit; the link-events of flagged events are left out of both stages and of
                 the scored set L-dryhigh. There is no fit without the exclusion any more.
    B3'          the after-rain residual at matched W, highest minus lowest third of preceding rain, from the drying K1 fit
                 (whatever rung is selected), over the short links of the unflagged learning events. Null = 5 replicates
                 of (c) (one state), signal = 5 replicates of (e) (two stores); separable if every (e) is above max (c).
    naming       a returned drying K1 is "non-exponential drying: the drying rate falls as the antenna dries".

Order: (a), (c), (d), (d-high), (e), (f) [a FAIL on (a), (c), (d), (d-high) or a spurious k_dry in (f) is a STOP];
then replicates 1-4 of (c) and (e) for B3'; then (g) and the detectability study at the ratios 3, 2, 1.5 (ratio 3 is
truth (c): its 5 replicates are those of B3'). One child process per run, one at a time; peak memory and processor time
are logged; a clean stop if the machine is short of memory. The parent keeps Windows awake while it runs.
Fake data only. No real attenuation is used as a target. VAL and TEST are not opened.
"""
import atexit
import ctypes
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
from research.waa_ksurface import a6_recovery as a6, data, fake6, ladder, s1_separability as s1, s3_recovery as s3
from research.waa_ksurface import paths as P
from research.wet_antenna.system_id import PADM

A6, S1, B = a6.A6, s1.S1, data.branch_config()["s3b"]
TRUTHS = ["a", "c", "d", "d_high", "e", "f"]                    # Deviation 12: every truth again, with the flag
RUNS = TRUTHS + ["g"]                                           # seed = base + position in this list
STOP_ON_FAIL = ["a", "c", "d", "d_high"]
REP_OFFSET = {"c": 20, "e": 30}                                 # B3' replicates 1-4: seed = base + offset + r - 1
K1_NAME = "non-exponential drying: the drying rate falls as the antenna dries"     # Deviation 13
GB = 2 ** 30
log = a6.log


# =================================================================================================== names and seeds
def rep_names(t):
    """The 5 B3' replicates of truth t ("c" or "e"): replicate 0 is the run of Deviation 12."""
    return [t] + [f"{t}_r{k}" for k in range(1, B["b3_replicates"])]


def det_name(a, rep):
    return a6.det_name("cw", a, rep)


def truth_of(name):
    """The truth a run draws: "a" ... "g", or "cw" for a replicate of the detectability study."""
    if name in RUNS:
        return name
    if name.startswith("cw_"):
        return "cw"
    return name.split("_r")[0]


def spec_seed(name):
    """(truth, seed) of a run (decisions.md, Deviation 12 reading 6 and Deviation 13)."""
    if name in RUNS:
        return A6["truths"][name], B["seed"] + RUNS.index(name)
    if name.startswith("cw_"):
        _, a, rep = name.split("_")
        return a6.spec_of(name)[0], B["seed"] + 100 + 10 * B["detect_sizes"].index(float(a)) + int(rep[1:])
    t, rep = name.split("_r")
    return A6["truths"][t], B["seed"] + REP_OFFSET[t] + int(rep) - 1


# =================================================================================================== one run
def run_one(name, S):
    """One run: fake data on all events, the flag, the two stages without the flagged events, judged as in step 3."""
    spec, seed = spec_seed(name)
    D, le, learn = S["D"], S["le"], S["learn"]
    t0, c0 = time.time(), time.process_time()
    sc = np.repeat(np.where(learn, 1.0, S1["radar_scale_added"]), D["le_n"] + 2 * PADM)
    Dt = dict(D)
    Dt["rain"] = np.asarray(D["rain"], float) / sc              # the true rain of an added event = its radar / 2.5
    Dk, ex = fake6.make(Dt, le, S["n_cells"], spec, S["curves"], S["iF"], S["iP"], seed, S["step"], S["dist"], ghost_pool=np.sort(le.event_id[learn].unique()))
    Dk["rain"] = Dk["rain"] * sc                                # the fit sees the radar value
    bi = ex["band_idx"]

    fl = a6.tail_flags(S, Dk)                                   # Deviation 12: the flag, before any fit
    fl["unseen_water_in_truth"] = fl.event_id.isin(ex["ghost_events"])
    flagged = set(fl.event_id[fl.flagged])
    use = (~le.event_id.isin(flagged)).to_numpy().astype(float)
    used = (use > 0)[S["bin_le"]]

    counted, dry, high, _ = s1.score_sets(S, Dk)
    w5 = np.where(learn[S["bin_le"]], counted, dry)             # stage 2: of an added link-event only its L-dry values count
    stage2 = {"w5": w5, "score": high & used, "margin_db": s3.margin(), "selector": s3.select_k0_k1, "rungs": s3.DRY_RUNGS}
    E = A6["extension"]
    start_F = [np.full(len(i), E["f_start_inc_db"]) for i in S["iF"]]
    start_P = [E["p_start_share"] * i for i in S["iP"]]
    res = ladder.two_stage(Dk, le, bi, S["curves"], S["gauge"], inc_F=start_F, inc_P=start_P, hold_ext=False, use=use, log=log, tag=f"({name}) ",
                           stage2=stage2, wet_rungs=s3.WET_RUNGS)
    dens = a6.state_density(S, Dk, ex)
    out = {"name": name, "truth": truth_of(name), "spec": spec, "seed": seed, "flags": fl,
           "events_flagged": int(fl.flagged.sum()), "events_checked": int(len(fl)),
           "flagged_with_unseen_water": int((fl.flagged & fl.unseen_water_in_truth).sum()),
           "false_removals": int((fl.flagged & ~fl.unseen_water_in_truth).sum()), "unseen_water_events": len(ex["ghost_events"]),
           "scored_values": int((high & used).sum()), "scored_link_events": int(len(np.unique(S["bin_le"][high & used])))}
    out["cv"] = pd.concat([res["wet"]["table"], res["dry"]["table"]], ignore_index=True)
    out["tau_wet"] = res.get("tau_wet")
    out["judge"], _ = a6.judge(S, spec, res, dens, True)
    Su = dict(S, short=S["short"] & (use > 0))                 # B3 and B3': the unflagged learning events only
    out["after_rain_selected"] = a6.after_rain(Su, Dk, res)     # old B3: from the selected fit
    m1, x1 = res["dry"]["fits"]["K1"]
    out["b3_prime"] = a6.after_rain(Su, Dk, {"model": m1, "x": x1})     # B3': from the K1 fit
    out["minutes"] = (time.time() - t0) / 60
    out["cpu_minutes"] = (time.process_time() - c0) / 60
    return out


# =================================================================================================== verdicts
def drying(r):
    return r["judge"].set_index("regime").loc["drying"]


def verdict(name, r):
    """The verdict of a run (Deviation 12 reading 4): only drying is judged, on the fit of the procedure."""
    t, d = truth_of(name), drying(r)
    sel, k_ok = d.selected_rung, bool(d.k_ok) if d.k_ok is not None else False
    if t == "e":
        return "report"
    if t in ("c", "cw"):
        return "PASS" if (sel == "K1" and k_ok) else "FAIL"
    if t == "d":
        ms = A6["judge"]["min_share"]
        ok = sel == "K0" and d.share_within_25pct >= ms and (not np.isfinite(d.share_alt) or d.share_alt >= ms)
        return "PASS" if ok else "FAIL"
    k0_ok = sel == "K0" and k_ok
    if t in ("a", "d_high"):
        return "PASS" if k0_ok else "FAIL"
    if t == "f":
        return "no spurious k_dry" if k0_ok else "spurious k_dry"
    return "k not corrupted" if k0_ok else "k corrupted"        # (g)


def stops(name, v):
    return (name in STOP_ON_FAIL and v == "FAIL") or (name == "f" and v.startswith("spurious"))


def separable(c_vals, e_vals):
    """B3': separable if every (e) value is above the largest (c) value; a missing (e) value counts as not above.
    Returns (separable, threshold = max (c))."""
    c, e = np.asarray(c_vals, float), np.asarray(e_vals, float)
    if len(c) < B["b3_replicates"] or len(e) < B["b3_replicates"] or not np.isfinite(c).any():
        return False, np.nan
    thr = float(np.nanmax(c))
    return bool(np.all(np.isfinite(e) & (e > thr))), thr


def recovery_text(r):
    d = drying(r)
    if d.true_rung in ("K1", "K2"):
        return f"{100 * d.share_within_25pct:.0f} % / {100 * d.share_alt:.0f} % of the domain weight within 25 %"
    if d.true_rung == "none":
        return f"1 / k_dry {d.tau_fit_min:.1f} / {d.tau_fit_alt_min:.1f} min (fit; no single truth)"
    return f"1 / k_dry {d.tau_fit_min:.1f} / {d.tau_fit_alt_min:.1f} min, truth {d.tau_true_min:.0f}"


# =================================================================================================== keeping the machine awake, watching memory
ES_CONTINUOUS, ES_SYSTEM_REQUIRED = 0x80000000, 0x00000001


def keep_awake():
    """Ask Windows not to sleep while this (main) thread lives; cleared at exit. Returns True if Windows accepted."""
    if sys.platform != "win32":
        return False
    k = ctypes.windll.kernel32
    k.SetThreadExecutionState.restype = ctypes.c_uint
    ok = k.SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED) != 0
    atexit.register(lambda: k.SetThreadExecutionState(ES_CONTINUOUS))
    return ok


def programs():
    """Private memory (GB) per program name, summed over its processes."""
    out = {}
    for p in psutil.process_iter(["name"]):
        try:
            mi = p.memory_info()
            out[p.info["name"]] = out.get(p.info["name"], 0.0) + getattr(mi, "private", mi.rss) / GB
        except (psutil.Error, OSError):
            pass
    return out


def watch_memory(state):
    """Every minute: free memory and the gap since the last look (a gap far above a minute = the machine slept). Every
    `memory_log_minutes` a line in the log; below `min_free_gb` the programs that grew most since the start."""
    base, last, last_line, last_low = programs(), time.time(), 0.0, 0.0
    while True:
        time.sleep(60)
        now = time.time()
        free = psutil.virtual_memory().available / GB
        state["min_free"] = min(state["min_free"], free)
        if now - last > 180:
            state["gaps"].append((now - last) / 60)
            log(f"=== memory watch: {(now - last) / 60:.0f} min since the last look (the machine was asleep or suspended) ===")
        last = now
        if now - last_line >= 60 * B["memory_log_minutes"]:
            log(f"=== memory watch {time.strftime('%H:%M')}: free {free:.1f} GB ===")
            last_line = now
        if free < B["min_free_gb"] and now - last_low >= 600:
            cur = programs()
            grew = sorted(((cur[n] - base.get(n, 0.0), n, base.get(n, 0.0), cur[n]) for n in cur), reverse=True)[:5]
            log(f"=== LOW MEMORY {time.strftime('%H:%M')}: free {free:.2f} GB; grew most since the start: " +
                "; ".join(f"{n} {b:.2f} -> {c:.2f} GB ({g:+.2f})" for g, n, b, c in grew) + " ===")
            state["low"].append((time.strftime("%H:%M"), free, grew))
            last_low = now


# =================================================================================================== report
def report(R, det, state, stopped):
    rows, flags, cvs = [], [], []
    for n, r in R.items():
        d = drying(r)
        g, lo, hi = s3.k1_gain(r["cv"])
        rows.append({"run": n, "truth": r["truth"], "seed": r["seed"], "drying rung": d.selected_rung, "true rung": d.true_rung,
                     "K1 gain on L-dryhigh (dB)": round(g, 5), "interval": f"{lo:+.5f} to {hi:+.5f}", "recovery (k_wet gauge / radar)": recovery_text(r),
                     "events flagged": f"{r['events_flagged']} of {r['events_checked']}", "flagged with unseen water": r["flagged_with_unseen_water"],
                     "false removals": r["false_removals"], "scored values": r["scored_values"], "verdict": verdict(n, r),
                     "B3' (K1 fit, dB)": round(r["b3_prime"]["matched_high_minus_low_db"], 4), "B3 (selected fit, dB)": round(r["after_rain_selected"]["matched_high_minus_low_db"], 4),
                     "minutes": round(r["minutes"], 1), "cpu minutes": round(r["cpu_minutes"], 1), "peak memory (GB)": round(r.get("peak_gb", np.nan), 2)})
        flags.append(r["flags"].assign(run=n))
        cvs.append(r["cv"].assign(run=n))
    runs = pd.DataFrame(rows)
    save(runs, P.S3B_RUNS, index=False)
    if flags:
        save(pd.concat(flags, ignore_index=True), P.S3B_FLAGS, index=False)
        save(pd.concat(cvs, ignore_index=True), P.S3B_CV, index=False)
    m = s3.margin()
    L = ["# Step 3b: the tail flag in the procedure (Deviation 12), B3' calibrated (Deviation 13), (g), detectability\n",
         "Written by `research/waa_ksurface/s3b_recovery.py` (libraries `s3_recovery.py`, `a6_recovery.py`, `s1_separability.py`, `fake6.py`, `ladder.py`). "
         "Fake attenuation on all TRAIN link-events of the branch, as in step 3 (the radar of the events without a gauge ratio reads "
         f"{S1['radar_scale_added']} times the true rain), fresh seeds (base {B['seed']}). Drying K0 against K1 on L-dryhigh, margin {m:.5f} dB with a paired 90 % "
         "interval that excludes zero; wetting K0, not judged (R3). **New (Deviation 12): the tail flag is computed on the fake data before any fit, and the "
         "link-events of the flagged events are left out of both stages and of the scored set.** Where two numbers are given as `x / y`, the first is with "
         "k_wet at its gauge-input estimate, the second at its radar-input estimate. A returned drying K1 is named: *" + K1_NAME + "* (Deviation 13). "
         "No real attenuation is used as a target. VAL and TEST were not opened.\n"]
    main_rows = runs[runs.run.isin(TRUTHS + ["g"])]
    L += ["## Every truth with the flag (Deviation 12)\n",
          md_table(main_rows[["run", "seed", "true rung", "drying rung", "K1 gain on L-dryhigh (dB)", "interval", "recovery (k_wet gauge / radar)",
                              "events flagged", "flagged with unseen water", "false removals", "scored values", "verdict"]]), ""]
    if stopped:
        L.append(f"\n**STOP by the rule: {stopped}.** Nothing after it was run.\n")
    # B3'
    c_runs = [n for n in rep_names("c") if n in R]
    e_runs = [n for n in rep_names("e") if n in R]
    b3 = pd.DataFrame([{"truth": truth_of(n), "run": n, "seed": R[n]["seed"], "drying rung": drying(R[n]).selected_rung,
                        "B3' statistic, K1 fit (dB)": R[n]["b3_prime"]["matched_high_minus_low_db"], "matched bins": R[n]["b3_prime"]["matched_bins"],
                        "link-events": R[n]["b3_prime"]["link_events"], "plain high - low, K1 fit (dB)": R[n]["b3_prime"]["high_minus_low_db"],
                        "old B3, selected fit (dB)": R[n]["after_rain_selected"]["matched_high_minus_low_db"]} for n in c_runs + e_runs])
    if len(b3):
        save(b3, P.S3B_B3, index=False)
        L += ["\n## B3': the history test, calibrated (Deviation 13)\n",
              "Statistic: the after-rain residual (30 min after the last wet bin, short links of the unflagged learning events) at matched W (bins of 0.3 dB of the "
              "fitted W at the end of the rain), highest minus lowest third of preceding rain, from the drying K1 fit. Null = (c) (one state, k_dry rising with W); "
              "signal = (e) (two stores). Separable if every (e) value is above the largest (c) value.\n",
              md_table(b3.round(5)), ""]
        if len(c_runs) == len(e_runs) == B["b3_replicates"]:
            cv = b3[b3.truth == "c"]["B3' statistic, K1 fit (dB)"].to_numpy()
            ev = b3[b3.truth == "e"]["B3' statistic, K1 fit (dB)"].to_numpy()
            sep, thr = separable(cv, ev)
            L.append(f"(c): {np.nanmin(cv):+.4f} to {np.nanmax(cv):+.4f} dB; (e): {np.nanmin(ev):+.4f} to {np.nanmax(ev):+.4f} dB; "
                     f"(e) above max (c) in {int(np.sum(np.isfinite(ev) & (ev > thr)))} of {len(ev)}.\n")
            L.append(f"**Separable: B3' replaces B3 in Block B, threshold = max (c) = {thr:+.4f} dB.**\n" if sep else
                     "**Not separable: finding \"W-dependent rate vs second store not identifiable with this data\"; P5 becomes restriction R5.** "
                     "B3's 0.10 dB threshold is retired.\n")
        else:
            L.append(f"Not complete: {len(c_runs)} (c) and {len(e_runs)} (e) replicates.\n")
    if "g" in R:
        v = verdict("g", R["g"])
        L += ["\n## (g): a larger spread between links\n",
              f"Drying rung {drying(R['g']).selected_rung}, {recovery_text(R['g'])}: **{v}**; per-link question reopened: {'yes' if v == 'k corrupted' else 'no'}.\n"]
    if det is not None and len(det):
        save(det, P.S3B_DETECT, index=False)
        L += ["\n## Detectability of k_dry(W) under the updated procedure\n",
              "k_dry = (1/15)·(1 + a·W) / (1 + a), k_wet = 1/5; ratio k(W = 2 dB) / k(W = 0) = 1 + 2a. Detected = drying K1 returned in at least 4 of 5 replicates "
              "(the replicates of a ratio stop when decided; the sweep ends at the first ratio not detected). Ratio 3 is truth (c): its replicates are those of B3'.\n",
              md_table(det), ""]
        d_ok = det[det.detected]
        L.append(f"**Smallest ratio detected in at least 4 of 5: {d_ok.ratio.min():g}.**\n" if len(d_ok) else "**No ratio detected in at least 4 of 5.**\n")
    L += ["\n## Execution\n",
          f"Runs: {len(R)}; wall time of the runs {runs.minutes.sum() / 60:.1f} h, processor time {runs['cpu minutes'].sum() / 60:.1f} h; largest peak memory "
          f"{runs['peak memory (GB)'].max():.2f} GB.\n" if len(runs) else "No run.\n",
          md_table(runs[["run", "minutes", "cpu minutes", "peak memory (GB)", "B3' (K1 fit, dB)", "B3 (selected fit, dB)"]]) if len(runs) else "", ""]
    if state:
        L.append(f"Keep-awake request accepted by Windows: {state.get('awake')}; gaps of the memory watch above 3 min: {len(state['gaps'])} "
                 f"({', '.join(f'{g:.0f} min' for g in state['gaps']) or 'none'}); lowest free memory seen {state['min_free']:.2f} GB; low-memory events: {len(state['low'])}.\n")
        if state.get("incomplete"):
            L.append(f"\n**NOT COMPLETE: {state['incomplete']}. Run again to go on from the cache.**\n")
    save("\n".join(L) + "\n", P.S3B_REPORT)
    log(f"report: {rel(P.S3B_REPORT)}")


# =================================================================================================== main
def main():
    cache = P.CACHE / "s3b_recovery"
    cache.mkdir(parents=True, exist_ok=True)
    if "--one" in sys.argv:                                     # one run in this process: save it and its peak memory
        name = sys.argv[sys.argv.index("--one") + 1]
        r = run_one(name, s3.setting())
        r["peak_gb"] = psutil.Process().memory_info().peak_wset / GB
        (cache / f"{name}.pkl").write_bytes(pickle.dumps(r))
        log(f"=== ({name}): {verdict(name, r)}; drying {drying(r).selected_rung}; flagged {r['events_flagged']} of {r['events_checked']} "
            f"(false {r['false_removals']}); B3' {r['b3_prime']['matched_high_minus_low_db']:+.4f} dB; {r['minutes']:.1f} min wall, "
            f"{r['cpu_minutes']:.1f} min processor; peak memory {r['peak_gb']:.2f} GB ===")
        return
    only_report = "--report" in sys.argv
    state = {"min_free": np.inf, "gaps": [], "low": [], "peak": 0.0, "incomplete": ""}
    if not only_report:
        free = psutil.virtual_memory().available / GB
        log(f"=== step 3b started {time.strftime('%Y-%m-%d %H:%M:%S')}; free memory {free:.1f} GB ===")
        if free < B["min_free_start_gb"]:
            log(f"=== not started: free memory {free:.1f} GB is below {B['min_free_start_gb']} GB ===")
            sys.exit(3)
        state["awake"] = keep_awake()
        log(f"=== keep-awake (ES_CONTINUOUS | ES_SYSTEM_REQUIRED): {'on' if state['awake'] else 'NOT accepted'} ===")
        state["min_free"] = free
        threading.Thread(target=watch_memory, args=(state,), daemon=True).start()
    R = {}

    def get(name):
        """The result of a run: from the cache, or made in its own process. None if it cannot be made now."""
        f = cache / f"{name}.pkl"
        if not f.exists():
            if only_report or state["incomplete"]:
                return None
            free = psutil.virtual_memory().available / GB
            need = max(B["min_free_gb"], state["peak"])
            if free < need:
                state["incomplete"] = f"free memory {free:.1f} GB is below the {need:.1f} GB a run needs: stopped before ({name})"
                log("=== " + state["incomplete"] + " ===")
                return None
            log(f"=== ({name}) seed {spec_seed(name)[1]} (free memory {free:.1f} GB) {time.strftime('%H:%M:%S')} ===")
            subprocess.run([sys.executable, "-u", "-m", "research.waa_ksurface.s3b_recovery", "--one", name], check=True)
        r = pickle.loads(f.read_bytes())
        state["peak"] = max(state["peak"], r.get("peak_gb", 0.0))
        R[name] = r
        return r

    stopped = ""
    for n in TRUTHS:                                            # 1. Deviation 12
        r = get(n)
        if r is None:
            break
        if stops(n, verdict(n, r)):
            stopped = f"({n}) {verdict(n, r)}"
            break
    det = None
    if not stopped and all(n in R for n in TRUTHS):
        for t in ("c", "e"):                                    # 2. B3' replicates
            for n in rep_names(t)[1:]:
                if get(n) is None:
                    break
        if all(n in R for t in ("c", "e") for n in rep_names(t)) and get("g") is not None:     # 3. (g), then detectability
            rows, reps = [], B["b3_replicates"]
            for a in B["detect_sizes"]:
                names = rep_names("c") if a == 1.0 else [det_name(a, k) for k in range(reps)]
                hits, n_run, missing = 0, 0, False
                for nm in names:
                    r = get(nm)
                    if r is None:
                        missing = True
                        break
                    n_run += 1
                    hits += int(drying(r).selected_rung == "K1")
                    if a != 1.0 and (hits >= reps - 1 or n_run - hits >= 2):
                        break                                   # decided
                rows.append({"size_a": a, "ratio": 1 + 2 * a, "replicates_run": n_run, "hits (K1 returned)": hits, "detected": hits >= reps - 1,
                             "runs": ", ".join(names[:n_run])})
                if missing or hits < reps - 1:
                    break
            det = pd.DataFrame(rows)
    report(R, det, None if only_report else state, stopped)
    log(f"=== step 3b ended {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
    if stopped or state["incomplete"]:
        sys.exit(1 if stopped else 2)


if __name__ == "__main__":
    main()
