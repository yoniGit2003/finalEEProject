"""Step 3 after S1: the fake-data truths again, with the drying rung selected on the adopted loss (Deviation 9).

    python -m research.waa_ksurface.s3_recovery                 every truth in order, then the detectability study
    python -m research.waa_ksurface.s3_recovery --report        write the report from the cache
    python -m research.waa_ksurface.s3_recovery --one c         one truth or replicate in this process (used by the main run)

What differs from `a6_recovery.py` (decisions.md, "S1" and "Step 3 under Deviation 9", written before this code):

    fake world   all TRAIN events of the branch; the radar of the events without a gauge ratio reads 2.5 times the true rain
    stage 1      unchanged: wetting ladder on the gauge-input set of the 14 learning events, margin 0.005 dB (Deviation 5)
    stage 2      the learning set plus the L-dry values of the other events; every drying rung is scored on L-dryhigh
                 (the dry 120 min after the rain of the link-events whose first dry 5 min value is at least 1 dB);
                 a richer rung needs a gain above the margin measured in S1 (the largest gain on a constant truth)
    judging      on the learning set, as in A6

Order: (c), (d), (d-high), (e), (f), (g), then the detectability study. A FAIL on (c) or (d), or a spurious k_dry in (f),
is a STOP. One process per truth, one at a time; peak memory is logged; a clean stop if the machine is short of memory.
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
    stage2 = {"w5": w5, "score": high, "margin_db": margin()}
    E = A6["extension"]
    start_F = [np.full(len(i), E["f_start_inc_db"]) for i in S["iF"]]
    start_P = [E["p_start_share"] * i for i in S["iP"]]
    tag = f"({name}) "
    out = {"name": name, "spec": spec, "seed": seed, "scored_values": int(high.sum()), "scored_link_events": int(len(np.unique(S["bin_le"][high])))}
    res = ladder.two_stage(Dk, le, bi, S["curves"], S["gauge"], inc_F=start_F, inc_P=start_P, hold_ext=False, log=log, tag=tag, stage2=stage2)
    dens = a6.state_density(S, Dk, ex)
    out["cv"] = pd.concat([res["wet"]["table"], res["dry"]["table"]], ignore_index=True)
    out["tau_wet"] = res.get("tau_wet")
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
                                stage2={**stage2, "score": high & used})
        out["judge_excluded"], out["cells_excluded"] = a6.judge(S, spec, res2, dens, True)
        out["cv_excluded"] = pd.concat([res2["wet"]["table"], res2["dry"]["table"]], ignore_index=True)
    out["minutes"] = (time.time() - t0) / 60
    return out


def stops(name, v):
    """The extra STOP rule of this step: a spurious k_dry in (f)."""
    return name == "f" and v.startswith("spurious")


def best_richer(cv, regime="drying"):
    """The largest gain of a richer rung of a regime over the rung it was compared with (dB), and that rung."""
    g = cv[(cv["select"] == regime) & cv.best_lambda & cv.gain_db.notna()]
    i = g.gain_db.astype(float).idxmax()
    return g.rung[i], float(g.gain_db[i]), float(g.gain_lo[i]), float(g.gain_hi[i])


def main():
    only_report = "--report" in sys.argv
    cache = P.CACHE / "s3_recovery"
    cache.mkdir(parents=True, exist_ok=True)
    if "--one" in sys.argv:                                     # one truth in this process: save it and its peak memory
        name = sys.argv[sys.argv.index("--one") + 1]
        r = run_truth(name, setting())
        r["peak_gb"] = psutil.Process().memory_info().peak_wset / 2 ** 30
        (cache / f"{name}.pkl").write_bytes(pickle.dumps(r))
        log(f"=== ({name}): {a6.verdict(name, r) if name in ORDER else 'replicate'}; {r['minutes']:.1f} min; peak memory {r['peak_gb']:.2f} GB ===")
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
        v = a6.verdict(n, r)
        if (n in a6.GATES and v == "FAIL") or stops(n, v):
            failed = True
            break                                               # the rule: a STOP
    det = None
    if not failed and all(n in R for n in ORDER):               # the detectability study, with the same rules as in A6 final
        Dt, rows = A6["detect"], []
        for form, regimes in a6.FORMS.items():
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
             f"other events; a richer rung needs a gain above {m:.5f} dB (the largest gain on a constant truth in S1) with a paired 90 % interval that excludes zero.** In the tables "
             "`cv_db` of a drying rung is its error on L-dryhigh. Wetting: Deviation 5 (0.005 dB on every counted value of the gauge-input set). Where two numbers are given as "
             "`x / y`, the first is with k_wet at its gauge-input estimate, the second at its radar-input estimate. Judged on the learning set as in A6. "
             "No real attenuation is used as a target. VAL and TEST were not opened.\n")
    stop = a6.report({n: R[n] for n in ORDER if n in R}, det, P=Files, title="# Step 3: the fake-data truths with the drying selection of Deviation 9\n", intro=intro, stop_on=stops)
    extra = []
    if "f" in R:
        rung, g, lo, hi = best_richer(R["f"]["cv"])
        extra += ["\n## (f): the spurious gain on the adopted loss against the margin\n",
                  f"Largest gain of a richer drying rung on L-dryhigh: {rung}, {g:+.5f} dB (interval {lo:+.5f} to {hi:+.5f}); margin {m:.5f} dB. "
                  f"Verdict of (f): {a6.verdict('f', R['f'])}.\n"]
    rows = [{"truth": n, "drying rung": R[n]["judge"].set_index("regime").selected_rung["drying"], "best richer drying rung": best_richer(R[n]["cv"])[0],
             "its gain on L-dryhigh (dB)": round(best_richer(R[n]["cv"])[1], 5), "margin (dB)": round(m, 5), "scored values": R[n]["scored_values"],
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
