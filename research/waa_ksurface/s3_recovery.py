"""Step 3 after S1: the fake-data truths again, with the drying rung selected on the adopted loss (Deviation 9).

    python -m research.waa_ksurface.s3_recovery                 every truth in order, then the detectability study
    python -m research.waa_ksurface.s3_recovery --report        write the report from the cache
    python -m research.waa_ksurface.s3_recovery --one c         one truth or replicate in this process (used by the main run)

What differs from `a6_recovery.py` (decisions.md, "S1" and "Step 3 under Deviation 9", written before this code):

    fake world   all TRAIN events of the branch; the radar of the events without a gauge ratio reads 2.5 times the true rain
    stage 1      wetting K0 only, on the gauge-input set of the 14 learning events (Deviation 11: no wetting ladder;
                 k_wet is not identified and is carried as a nuisance at its gauge-input and its radar-input estimate)
    stage 2      the learning set plus the L-dry values of the other events; every drying rung is scored on L-dryhigh
                 (the dry 120 min after the rain of the link-events whose first dry 5 min value is at least 1 dB);
                 the drying ladder is K0 against K1 only (Deviation 11, `select_k0_k1`): K1 needs a gain above the
                 margin measured in S1 (the largest gain on a constant truth) with an interval that excludes zero
    judging      on the learning set, as in A6

Order: (c), (d), (d-high), (e), (f), (g), then the detectability study of k_dry(W). In (d) and (d-high) the true k depends
on R, which the reduced ladder cannot return (restriction R4): they pass if drying stays K0 and the constant is recovered.
A FAIL on (c), (d) or (d-high), or a spurious k_dry in (f), is a STOP. The cached (c) and (d) were fitted with the full
ladders; their K0 and K1 fits are those of the reduced procedure and are judged again, not refitted. One process per truth, one at a time; peak memory is logged; a clean stop if the machine is short of memory.
Fake data only. No real attenuation is used as a target.
"""
import pickle
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import psutil

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a6_recovery as a6, data, fake6, ladder, s1_separability as s1
from research.waa_ksurface import paths as P
from research.wet_antenna.system_id import PADM

A6, S1, S3 = a6.A6, s1.S1, data.branch_config()["s3"]
ORDER = ["c", "d", "d_high", "e", "f", "g"]
WET_RUNGS, DRY_RUNGS = ["K0"], ["K0", "K1"]                     # Deviation 11: no wetting ladder; drying K0 against K1
FORMS = {"cw": ("drying",)}                                     # detectability: k_dry(W) only (R4: no k(R) is claimed)
STOP_ON_FAIL = ["c", "d", "d_high"]
MIN_FREE_GB = 2.0
log = a6.log


class Files:
    """Where `a6_recovery.report` writes for this step (the same roles as the A6 files)."""
    A6_VERDICTS, A6_CV, A6_CELLS, A6_EXTENSION = P.S3_VERDICTS, P.S3_CV, P.S3_CELLS, P.S3_EXTENSION
    A6_FLAGS, A6_LAG, A6_DETECT_RUNS, A6_DETECT, A6_FIG, A6_REPORT = P.S3_FLAGS, None, P.S3_DETECT_RUNS, P.S3_DETECT, P.S3_FIG, P.S3_REPORT


def margin():
    """The margin of Deviation 9: the largest gain of K1 over K0 on a constant truth, on the adopted loss and variant (S1)."""
    t = pd.read_csv(P.S1_TABLE)
    r = t[(t.loss == S3["loss"]) & (t.variant == S3["variant"])].iloc[0]
    assert bool(r.separable)
    return max(float(r.a_gain_max_db), 0.0)


UNRESOLVED = "dependence detected, attribution W vs R not resolved"


def select_nonnested(compare, margin):
    """The drying selection of Deviation 10. compare(richer, simpler) = (gain in dB, lower end, upper end of its paired
    interval); a rung "beats" another if its gain is above the margin and the interval excludes zero.

        neither K1 nor K2 beats K0             -> K0, "no dependence" (K3 is not tried)
        one or both beat K0                    -> the one with the larger gain: K1 "W" or K2 "R"
        both beat K0, gains closer than margin -> the better one, "dependence detected, attribution W vs R not resolved"
        K3 beats the rung selected so far      -> K3, "W and R"

    Returns (rung, attribution, steps); steps[rung] = (compared with, gain, lower, upper, selected) for the table."""
    def beats(c):
        return bool(c[0] > margin and not (c[1] <= 0.0 <= c[2]))

    c = {"K1": compare("K1", "K0"), "K2": compare("K2", "K0")}
    won = [r for r in ("K1", "K2") if beats(c[r])]
    steps = {r: ("K0", *c[r], False) for r in ("K1", "K2")}
    if not won:
        return "K0", "no dependence", steps
    rung = max(won, key=lambda r: c[r][0])                      # the larger gain (K1 if they are equal)
    resolved = len(won) == 1 or abs(c["K1"][0] - c["K2"][0]) >= margin
    c3 = compare("K3", rung)
    steps["K3"] = (rung, *c3, beats(c3))
    if beats(c3):
        return "K3", "W and R", steps
    steps[rung] = ("K0", *c[rung], True)
    return rung, ({"K1": "W", "K2": "R"}[rung] if resolved else UNRESOLVED), steps


def select_k0_k1(compare, margin):
    """The drying selection of Deviation 11: K1 against K0 and nothing else. K1 is returned if its gain is above the margin
    and the paired interval excludes zero. Returns (rung, attribution, steps) like `select_nonnested`."""
    g, lo, hi = compare("K1", "K0")
    take = bool(g > margin and not (lo <= 0.0 <= hi))
    return ("K1", "W", {"K1": ("K0", g, lo, hi, True)}) if take else ("K0", "no dependence", {"K1": ("K0", g, lo, hi, False)})


def verdict(name, r):
    """The verdict of one truth under Deviation 11. Only drying is judged (wetting is K0 by construction, k_wet is not
    identified). (c): K1 and the surface recovered. (d), (d-high): K0 (no spurious W-dependence) and k_dry recovered."""
    if name == "e":
        return "report"
    d = r["judge"].set_index("regime").loc["drying"]
    sel = d.selected_rung
    if name == "c":
        return "PASS" if (sel == "K1" and bool(d.k_ok)) else "FAIL"
    if name == "d":                                             # the true k depends on R (outside the ladder): a robustness check
        ok = sel == "K0" and d.share_within_25pct >= A6["judge"]["min_share"] and (not np.isfinite(d.share_alt) or d.share_alt >= A6["judge"]["min_share"])
        return "PASS" if ok else "FAIL"
    k0_ok = sel == "K0" and bool(d.k_ok)                        # the remaining truths have a constant k_dry (inside F's measured range)
    if name == "d_high":
        return "PASS" if k0_ok else "FAIL"
    if name == "f":
        d2 = r["judge_excluded"].set_index("regime").loc["drying"]
        return ("no spurious k_dry" if k0_ok else "spurious k_dry") + (", restored by the exclusion" if (d2.selected_rung == "K0" and bool(d2.k_ok)) else ", NOT restored by the exclusion")
    if name == "g":
        return "k not corrupted" if k0_ok else "k corrupted"
    return "detectability"


def shown(r):
    """A result as the report shows it: the drying row only, and for (d) the note that its true rung is outside the ladder."""
    r = dict(r)
    for key in ("judge", "judge_excluded"):
        if key in r:
            j = r[key][r[key].regime == "drying"].copy()
            if r["name"] == "d":
                j["true_rung"] = "K2 (outside the ladder, R4)"
            r[key] = j
    return r


def setting():
    """The S1 setting (all TRAIN events) with what the judging of A6 reads: short links and the gauge-input set, both of
    the learning events only, and the frozen A4 domain."""
    S = s1.setting()
    dom = pd.read_csv(P.A4_DOMAIN)
    S["dom"] = dom[(dom.band == "pooled") & dom.in_domain][["regime", "rcol", "wbin", "R_class", "W_lo"]]
    S["short"] = (S["le"].L_km < a6.LT["short_km"]).to_numpy() & S["learn"]
    S["gauge"] = S["gauge"].astype(float)
    return S


def seed_of(name):
    if name in ORDER:
        return S3["seed"] + ORDER.index(name)
    form, a, rep = name.split("_")
    return S3["seed"] + 100 + 1000 * list(a6.FORMS).index(form) + 10 * A6["detect"]["sizes"].index(float(a)) + int(rep[1:])


def run_truth(name, S):
    """One truth: fake data on all events, the two stages with the drying selection of Deviation 9, judged as in A6."""
    spec, seed = a6.spec_of(name)[0], seed_of(name)
    D, le, learn = S["D"], S["le"], S["learn"]
    t0 = time.time()
    sc = np.repeat(np.where(learn, 1.0, S1["radar_scale_added"]), D["le_n"] + 2 * PADM)
    Dt = dict(D)
    Dt["rain"] = np.asarray(D["rain"], float) / sc              # the true rain of an added event = its radar / 2.5
    Dk, ex = fake6.make(Dt, le, S["n_cells"], spec, S["curves"], S["iF"], S["iP"], seed, S["step"], S["dist"], ghost_pool=np.sort(le.event_id[learn].unique()))
    Dk["rain"] = Dk["rain"] * sc                                # the fit sees the radar value
    bi = ex["band_idx"]
    counted, dry, high, _ = s1.score_sets(S, Dk)
    w5 = np.where(learn[S["bin_le"]], counted, dry)             # stage 2: of an added link-event only its L-dry values count
    stage2 = {"w5": w5, "score": high, "margin_db": margin(), "selector": select_k0_k1, "rungs": DRY_RUNGS}
    E = A6["extension"]
    start_F = [np.full(len(i), E["f_start_inc_db"]) for i in S["iF"]]
    start_P = [E["p_start_share"] * i for i in S["iP"]]
    tag = f"({name}) "
    out = {"name": name, "spec": spec, "seed": seed, "scored_values": int(high.sum()), "scored_link_events": int(len(np.unique(S["bin_le"][high])))}
    res = ladder.two_stage(Dk, le, bi, S["curves"], S["gauge"], inc_F=start_F, inc_P=start_P, hold_ext=False, log=log, tag=tag, stage2=stage2, wet_rungs=WET_RUNGS)
    dens = a6.state_density(S, Dk, ex)
    out["cv"] = pd.concat([res["wet"]["table"], res["dry"]["table"]], ignore_index=True)
    out["tau_wet"] = res.get("tau_wet")
    out["attribution"] = res["dry"]["attribution"]
    out["after_rain"] = a6.after_rain(S, Dk, res)
    out["judge"], out["cells"] = a6.judge(S, spec, res, dens, True)
    out["extension"] = a6.extension_table(S, res["model"], res["x"], start_F, start_P)
    if "x_alt" in res:
        out["extension_alt"] = a6.extension_table(S, res["model"], res["x_alt"], start_F, start_P)
    if "ghost" in spec:                                         # exclude the flagged events and fit again
        fl = a6.tail_flags(S, Dk)
        fl["changed_in_truth"] = fl.event_id.isin(ex["ghost_events"])
        out["flags"] = fl
        use = (~le.event_id.isin(set(fl.event_id[fl.flagged]))).to_numpy().astype(float)
        used = (use > 0)[S["bin_le"]]
        res2 = ladder.two_stage(Dk, le, bi, S["curves"], S["gauge"], inc_F=start_F, inc_P=start_P, hold_ext=False, use=use, log=log, tag=tag + "excl. ",
                                stage2={**stage2, "score": high & used}, wet_rungs=WET_RUNGS)
        out["judge_excluded"], out["cells_excluded"] = a6.judge(S, spec, res2, dens, True)
        out["cv_excluded"] = pd.concat([res2["wet"]["table"], res2["dry"]["table"]], ignore_index=True)
        out["attribution_excluded"] = res2["dry"]["attribution"]
    out["minutes"] = (time.time() - t0) / 60
    return out


def stops(name, v):
    """The STOP rules of this step: a FAIL on (c), (d) or (d-high), or a spurious k_dry in (f)."""
    return (name in STOP_ON_FAIL and v == "FAIL") or (name == "f" and v.startswith("spurious"))


def k1_gain(cv):
    """The gain of drying K1 over K0 on the adopted loss (dB) and its interval."""
    g = cv[(cv["select"] == "drying") & (cv.rung == "K1") & cv.best_lambda].iloc[0]
    return float(g.gain_db), float(g.gain_lo), float(g.gain_hi)


def main():
    only_report = "--report" in sys.argv
    cache = P.CACHE / "s3_recovery"
    cache.mkdir(parents=True, exist_ok=True)
    if "--one" in sys.argv:                                     # one truth in this process: save it and its peak memory
        name = sys.argv[sys.argv.index("--one") + 1]
        r = run_truth(name, setting())
        r["peak_gb"] = psutil.Process().memory_info().peak_wset / 2 ** 30
        (cache / f"{name}.pkl").write_bytes(pickle.dumps(r))
        log(f"=== ({name}): {verdict(name, r) if name in ORDER else 'replicate'}; {r['minutes']:.1f} min; peak memory {r['peak_gb']:.2f} GB ===")
        return
    R, state = {}, {"stopped": "", "peak": 0.0}

    def get(name):
        """The result of a truth: from the cache, or made in its own process. None if it cannot be made now."""
        f = cache / f"{name}.pkl"
        if not f.exists():
            if only_report or state["stopped"]:
                return None
            free = psutil.virtual_memory().available / 2 ** 30
            need = max(MIN_FREE_GB, state["peak"])
            if free < need:
                state["stopped"] = f"free memory {free:.1f} GB is below the {need:.1f} GB a run needs: stopped before ({name})"
                log("=== " + state["stopped"] + " ===")
                return None
            log(f"=== ({name}) (free memory {free:.1f} GB) ===")
            subprocess.run([sys.executable, "-m", "research.waa_ksurface.s3_recovery", "--one", name], check=True)
        r = pickle.loads(f.read_bytes())
        state["peak"] = max(state["peak"], r.get("peak_gb", 0.0))
        return r

    failed = False
    for n in ORDER:
        r = get(n)
        if r is None:
            break
        R[n] = r
        v = verdict(n, r)
        if (n in STOP_ON_FAIL and v == "FAIL") or stops(n, v):
            failed = True
            break                                               # the rule: a STOP
    det = None
    if not failed and all(n in R for n in ORDER):               # the detectability study, with the same rules as in A6 final
        Dt, rows = A6["detect"], []
        for form, regimes in FORMS.items():
            for a in Dt["sizes"]:
                hits, n_run, missing = {r: 0 for r in regimes}, 0, False
                for rep in range(Dt["replicates"]):
                    name = a6.det_name(form, a, rep)
                    r = get(name)
                    if r is None:
                        missing = True
                        break
                    R[name] = r
                    n_run += 1
                    j = r["judge"].set_index("regime")
                    for reg in regimes:
                        hits[reg] += int(j.rung_ok[reg])
                        rows.append({"form": form, "size_a": a, "replicate": rep, "regime": reg, "true_rung": j.true_rung[reg], "rung_returned": j.selected_rung[reg],
                                     "hit": bool(j.rung_ok[reg]), "gain_db_of_true_rung": a6.gain_of(r["cv"], reg, j.true_rung[reg]), "share_within_25pct": j.share_within_25pct[reg]})
                    if all(hits[reg] >= Dt["replicates"] - 1 or (n_run - hits[reg]) >= 2 for reg in regimes):
                        break                                   # decided for every regime
                if missing or not any(hits[reg] >= Dt["replicates"] - 1 for reg in regimes):
                    break                                       # not finished, or nothing detected at this size
        det = pd.DataFrame(rows)

    m = margin()
    intro = ("Written by `research/waa_ksurface/s3_recovery.py` (libraries `a6_recovery.py`, `s1_separability.py`, `fake6.py`, `ladder.py`, `ksurf.py`). Fake attenuation on all "
             f"TRAIN link-events of the branch (real links, real rain; the radar of the events without a gauge ratio reads {S1['radar_scale_added']} times the true rain), the noise and "
             "the input errors of A6, extension free (Deviation 7), fresh seeds. **Drying selection (Deviation 9): every drying rung is scored on L-dryhigh (held-out 5 min values of "
             "the dry 120 min after the rain, link-events whose first dry value is at least 1 dB), learning set = the 14 events with a gauge ratio plus the dry-path values of the "
             f"other events; a richer rung needs a gain above {m:.5f} dB (the largest gain on a constant truth in S1) with a paired 90 % interval that excludes zero. The drying ladder is K0 against K1 only, and there is no wetting "
             "ladder: wetting is K0 and is not judged, k_wet is not identified (Deviation 11, restrictions R3 and R4).** In (d) and (d-high) the true k depends on R, which this "
             "ladder cannot return: they PASS if drying stays K0 and the constant is recovered. (c) and (d) were fitted before Deviation 11 with the full ladders; their K0 and K1 "
             "fits are those of the reduced procedure, so they are judged again and not refitted, and their tables still list the K2 and K3 rows, which take no part. In the tables "
             "`cv_db` of a drying rung is its error on L-dryhigh. Wetting: Deviation 5 (0.005 dB on every counted value of the gauge-input set). Where two numbers are given as "
             "`x / y`, the first is with k_wet at its gauge-input estimate, the second at its radar-input estimate. Judged on the learning set as in A6. "
             "No real attenuation is used as a target. VAL and TEST were not opened.\n")
    stop = a6.report({n: shown(R[n]) for n in ORDER if n in R}, det, P=Files, title="# Step 3: the fake-data truths with the reduced ladders (Deviations 9 and 11)\n", intro=intro,
                     stop_on=stops, verdict=verdict)
    extra = []
    if "f" in R:
        g, lo, hi = k1_gain(R["f"]["cv"])
        g2, lo2, hi2 = k1_gain(R["f"]["cv_excluded"])
        extra += ["\n## (f): the spurious gain on the adopted loss against the margin\n",
                  f"Gain of drying K1 over K0 on L-dryhigh: {g:+.5f} dB (interval {lo:+.5f} to {hi:+.5f}); with the flagged events excluded {g2:+.5f} dB "
                  f"({lo2:+.5f} to {hi2:+.5f}); margin {m:.5f} dB. "
                  f"Verdict of (f): {verdict('f', R['f'])}.\n"]
    rows = [{"truth": n, "drying rung": R[n]["judge"].set_index("regime").selected_rung["drying"], "attribution (Deviation 10)": R[n]["attribution"], "K1 gain on L-dryhigh (dB)": round(k1_gain(R[n]["cv"])[0], 5), "margin (dB)": round(m, 5), "scored values": R[n]["scored_values"],
             "scored link-events": R[n]["scored_link_events"], "peak memory (GB)": round(R[n].get("peak_gb", np.nan), 2), "minutes": round(R[n]["minutes"], 1)} for n in R]
    extra += ["\n## Every run: the drying gain against the margin, the size of the scored set, peak memory\n", md_table(pd.DataFrame(rows)), ""]
    if state["stopped"]:
        extra.append(f"\n**NOT COMPLETE: {state['stopped']}. Run again to go on from the cache.**\n")
    with open(P.S3_REPORT, "a", encoding="utf-8") as f:
        f.write("\n".join(extra) + "\n")
    if stop or state["stopped"]:
        sys.exit(1 if stop else 2)


if __name__ == "__main__":
    main()
