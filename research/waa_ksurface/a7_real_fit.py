"""A7: the real fit, exactly as pre-registered (decisions.md, "A7: pre-registration of the real fit", commit bbb5d2f).

    python -m research.waa_ksurface.a7_real_fit --flag         step 1: the tail flag on the real learning events (no fit)
    python -m research.waa_ksurface.a7_real_fit                step 2: the fit (k_wet estimated, Deviation 15), then the bootstrap at both estimates
    python -m research.waa_ksurface.a7_real_fit --report       write the report from the cache
    python -m research.waa_ksurface.a7_real_fit --one fit      one part in this process (used by the main run: fit, boot_gauge, boot_radar)

What it does:
    1. the flag      the tail flag of Deviation 12 on the REAL attenuation of the 14 learning events; the events it removes
                     take no part in anything below. Its count is the first number of A7, saved before any fit.
    2. the fit       the step 3c procedure itself (`s3c_a7_variant.procedure`, Deviation 15) on the real TRAIN data of the
                     unflagged learning events (rain = radar / g_e, the data layer of this branch): stage 1 (gauge input,
                     wetting K0, k_wet estimated), stage 2 drying K0 against K1 on L-dryhigh (margin 0.00291 dB), selected
                     once; the selected rung refitted with k_wet at its radar-input estimate. F and P: the measured curves
                     (E4) fixed in their measured ranges, free monotone extensions above (started as in step 3c).
    3. bootstrap     at each k_wet estimate, the selected rung refitted on 200 event-bootstrap draws (seed 3900, the same
                     draws at both), 1 / k_dry at W = 0, 0.5, 1, 2 dB; interval 5-95 %; the reported interval is the union
                     of the two (R3 revised).
    4. readings      only the pre-registered sentences. VAL and TEST are not read.
"""
import pickle
import sys
import time

import numpy as np
import pandas as pd
import psutil

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a6_recovery as a6, data, fp, ksurf, s1_separability as s1, s3_recovery as s3, s3b_recovery as s3b
from research.waa_ksurface import s3c_a7_variant as s3c, s3d_kwet_fixed as s3d
from research.waa_ksurface import paths as P

A6, A7 = a6.A6, data.branch_config()["a7"]
GB = s3b.GB
log = a6.log
CACHE = P.CACHE / "a7_real_fit"
K0_TEXT = ("no dependence of the drying rate on W larger than a ratio k(W = 2 dB) / k(W = 0) of 2 is present (2 = the smallest ratio "
           "detected for this variant in step 3c; 1.5 was detected in 2 of 5 and is not excluded)")
K1_TEXT = ("non-exponential drying: the drying rate falls as the antenna dries; the mechanism (W-dependent rate or a second, slow store) is "
           "not claimed (R5); a spread of the baseline between links is excluded as its cause at the measured size (0.088 dB, step 3c (a) PASS, "
           "also at 0.15 dB)")


# =================================================================================================== the real setting
def real():
    """The branch TRAIN setting (all events; the fit uses the learning events only), the measured curves and the band index."""
    S = s3.setting()
    M = fp.measured()
    bands = list(M)
    band_idx = S["le"].band_name.map({b: i for i, b in enumerate(bands)}).to_numpy()
    return S, [M[b] for b in bands], bands, band_idx


def flag_step(S):
    """Step 1: the tail flag on the real attenuation of the learning events (short links). No fit."""
    fl = a6.tail_flags(S, S["D"])
    return fl


def use_of(S, fl):
    removed = set(fl.event_id[fl.flagged])
    return (S["learn"] & ~S["le"].event_id.isin(removed).to_numpy()).astype(float)


# =================================================================================================== the fit and the bootstrap
ESTIMATES = ("gauge", "radar")                                  # k_wet at its gauge-input and its radar-input estimate (Deviation 15)


def fit_part(S, curves, band_idx, use):
    """The step 3c procedure itself (`s3c_a7_variant.procedure`) on the real attenuation: the drying rung is selected once
    (k_wet at the gauge-input estimate) and refitted at the radar-input estimate (x_alt)."""
    D = S["D"]
    res, (high, used) = s3c.procedure(S, D, band_idx, curves, use, "(A7) ")
    counted = s1.score_sets(S, D)[0]
    st = res["dry"]
    g, lo, hi = s3.k1_gain(st["table"])
    m = res["model"]
    xs = {"gauge": res["x"], "radar": res["x_alt"]}
    return {"rung": st["rung"], "lam": st["lam"], "gain": g, "lo": lo, "hi": hi, "table": pd.concat([res["wet"]["table"], st["table"]], ignore_index=True),
            "tau_wet": res["tau_wet"], "x": xs,
            "tau_dry": {e: dict(zip(A7["w_report_db"], s3d.k_dry_at(m, xs[e], A7["w_report_db"]))) for e in ESTIMATES},
            "extension": {e: {f"{k[0]} band {k[1]}": v for k, v in m.extension(xs[e]).items()} for e in ESTIMATES},
            "scored_values": int((high & used).sum()), "scored_link_events": int(len(np.unique(S["bin_le"][high & used]))),
            "stage2_w5": counted & used}


def boot_part(S, curves, band_idx, use, fit, est):
    """The selected rung refitted on event-bootstrap draws with k_wet held at one estimate (extension free)."""
    D, le = S["D"], S["le"]
    prob = ksurf.prepare(D, band_idx, curves, D["rain"], 0)
    prob["w5"] = np.asarray(fit["stage2_w5"], float)
    m = ksurf.Model(prob, "K0", fit["rung"], fit["lam"])
    m.hold("gw", "hw")
    x0 = fit["x"][est]
    ev_all = le.event_id.to_numpy()
    events = np.sort(np.unique(ev_all[use > 0]))
    rng = np.random.default_rng(A7["boot_seed"])                # the same draws at both estimates
    rows = []
    t0 = time.time()
    for b in range(A7["boot_draws"]):
        pick = rng.choice(events, len(events), replace=True)
        cnt = pd.Series(pick).value_counts()
        w = use * pd.Series(ev_all).map(cnt).fillna(0).to_numpy()
        xb = m.fit(x0.copy(), w)
        rows.append({"draw": b, **{f"1/k_dry at W = {W:g} dB": v for W, v in zip(A7["w_report_db"], s3d.k_dry_at(m, xb, A7["w_report_db"]))},
                     "converged": m.last["converged"]})
        if (b + 1) % 20 == 0:
            log(f"   (A7) [k_wet {est} estimate] bootstrap {b + 1} of {A7['boot_draws']} ({(time.time() - t0) / 60:.1f} min)")
    return pd.DataFrame(rows)


# =================================================================================================== readings and report
def reading(fit):
    """(rung, sentence): the pre-registered reading. One rung is selected in this path (Deviation 15), valid at both estimates."""
    rung = fit["rung"]
    return rung, K0_TEXT if rung == "K0" else K1_TEXT


def report(fl, fit, boots):
    L = ["# A7: the real fit (pre-registered, decisions.md, commit bbb5d2f; Deviation 15)\n",
         "Written by `research/waa_ksurface/a7_real_fit.py`. Real TRAIN data of this branch (28 and 38 GHz), the 14 learning events (rain = radar / g_e). "
         "VAL and TEST were not read. The procedure is that of step 3c itself (`s3c_a7_variant.procedure`, Deviation 15): the tail flag first; wetting K0 "
         "with k_wet ESTIMATED (a nuisance, R3 revised); drying K0 against K1 on L-dryhigh, margin "
         f"{s3c.margin():.5f} dB with a paired 90 % interval that excludes zero, selected once with k_wet at its gauge-input estimate and refitted at "
         "its radar-input estimate; F and P measured parts fixed, extensions free (R1).\n"]
    checked, removed = len(fl), fl[fl.flagged]
    L += ["## 1. The flag on the real learning events (before any fit)\n",
          f"Events checked: {checked} (E20150831_0455 is not checkable). **Removed: {len(removed)}**"
          + (f" ({', '.join(removed.event_id)})" if len(removed) else "") + ".\n", md_table(fl), ""]
    if fit is None:
        save("\n".join(L) + "\n", P.A7_REPORT)
        return
    tw = fit["tau_wet"]
    L += ["## 2. k_wet (nuisance estimates, no claim; R3 revised)\n",
          f"1 / k_wet: gauge input {tw['gauge']:.2f} min, radar input {tw['radar']:.2f} min.\n",
          "## 3. The drying rung and k_dry\n",
          f"Scored values (L-dryhigh): {fit['scored_values']} in {fit['scored_link_events']} link-events. **Drying rung: {fit['rung']}** "
          f"(lambda {fit['lam']:g}); K1 gain over K0 {fit['gain']:+.5f} dB, paired 90 % interval {fit['lo']:+.5f} to {fit['hi']:+.5f}; margin "
          f"{s3c.margin():.5f} dB.\n"]
    rows = []
    for e in ESTIMATES:
        b = boots.get(e)
        r = {"k_wet estimate": e, "1/k_wet (min)": round(tw[e], 2)}
        for W in A7["w_report_db"]:
            col = f"1/k_dry at W = {W:g} dB"
            iv = "" if b is None else f" ({np.percentile(b[col], 5):.1f} to {np.percentile(b[col], 95):.1f})"
            r[col + " (min)"] = f"{fit['tau_dry'][e][W]:.1f}{iv}"
        r["bootstrap draws converged"] = "" if b is None else f"{int(b.converged.sum())} of {len(b)}"
        rows.append(r)
    tab = pd.DataFrame(rows)
    save(tab, P.A7_FIT, index=False)
    L += [f"1 / k_dry with its 5-95 % event-bootstrap interval ({A7['boot_draws']} draws, seed {A7['boot_seed']}, the selected rung refitted with "
          "k_wet held at that estimate) in brackets; for K0 the value does not depend on W.\n", md_table(tab), ""]
    if len(boots) == len(ESTIMATES):
        un = []
        for W in A7["w_report_db"]:
            col = f"1/k_dry at W = {W:g} dB"
            lo = min(np.percentile(boots[e][col], 5) for e in ESTIMATES)
            hi = max(np.percentile(boots[e][col], 95) for e in ESTIMATES)
            un.append({"W (dB)": W, "1/k_dry, union of the two intervals (min)": f"{lo:.1f} to {hi:.1f}",
                       "point values at gauge / radar estimate": " / ".join(f"{fit['tau_dry'][e][W]:.1f}" for e in ESTIMATES)})
        L += ["### k_dry over the two k_wet estimates (R3 revised)\n", md_table(pd.DataFrame(un)), ""]
    save(fit["table"], P.A7_CV, index=False)
    rung, text = reading(fit)
    L += ["## 4. Result under the pre-registered readings\n", f"Selected drying rung: {rung} (one selection, valid at both k_wet estimates).\n",
          f"**{text}.**\n", "Stated with every A7 result: the flag cannot check E20150831_0455, so unseen water in it, if any, would go undetected.\n"]
    ext = []
    for e in ESTIMATES:
        for k, v in fit["extension"][e].items():
            ext.append({"k_wet estimate": e, "curve and band index": k, "values at the knots (EFFECTIVE, R1)": ", ".join(f"{a:.3f}" for a in np.atleast_1d(v))})
    save(pd.DataFrame(ext), P.A7_EXT, index=False)
    L += ["\n## 5. Extension knots above the measured ranges (labeled EFFECTIVE, R1)\n",
          f"Band index: {', '.join(f'{i} = {b}' for i, b in enumerate(fp.measured()))}. These are effective quantities for this rain input, "
          f"not physical claims (`{rel(P.A7_EXT)}`).\n",
          md_table(pd.DataFrame(ext)), ""]
    save("\n".join(L) + "\n", P.A7_REPORT)
    log(f"report: {rel(P.A7_REPORT)}")


# =================================================================================================== main
def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    flag_file, fit_file = CACHE / "flag.pkl", CACHE / "fit.pkl"
    if "--flag" in sys.argv:                                    # step 1, before any fit
        S, _, _, _ = real()
        fl = flag_step(S)
        flag_file.write_bytes(pickle.dumps(fl))
        save(fl, P.A7_FLAG, index=False)
        log(f"=== A7 flag: {len(fl)} learning events checked; removed {int(fl.flagged.sum())}: {', '.join(fl.event_id[fl.flagged]) or 'none'} ===")
        report(fl, None, {})
        return
    assert flag_file.exists(), "run --flag first (the flag count is the first number of A7)"
    fl = pickle.loads(flag_file.read_bytes())
    if "--one" in sys.argv:
        part = sys.argv[sys.argv.index("--one") + 1]
        S, curves, _, band_idx = real()
        use = use_of(S, fl)
        t0 = time.time()
        if part == "fit":
            out = fit_part(S, curves, band_idx, use)
            fit_file.write_bytes(pickle.dumps(out))
            log(f"=== (A7) drying {out['rung']}; K1 gain {out['gain']:+.5f} dB ({out['lo']:+.5f} to {out['hi']:+.5f}); 1/k_wet gauge "
                f"{out['tau_wet']['gauge']:.2f} / radar {out['tau_wet']['radar']:.2f} min (nuisance) ===")
            for e in ESTIMATES:
                log(f"=== (A7) 1/k_dry at the {e} estimate: " + ", ".join(f"W {W:g}: {v:.1f}" for W, v in out["tau_dry"][e].items()) + " min ===")
        else:
            b = boot_part(S, curves, band_idx, use, pickle.loads(fit_file.read_bytes()), part.split("_")[1])
            (CACHE / f"{part}.pkl").write_bytes(pickle.dumps(b))
        log(f"=== (A7 {part}) done: {(time.time() - t0) / 60:.1f} min; peak memory {psutil.Process().memory_info().peak_wset / GB:.2f} GB ===")
        return
    parts = ["fit"] + [f"boot_{e}" for e in ESTIMATES]
    if "--report" not in sys.argv:
        s3d.start_state("A7")                                   # keep-awake, memory watch, refuses to start below the floor
        for part in parts:                                      # the fit writes fit.pkl; each bootstrap writes boot_<tau>.pkl
            f = fit_file if part == "fit" else CACHE / f"{part}.pkl"
            if f.exists():
                continue
            free = psutil.virtual_memory().available / GB
            if free < s3d.D3["min_free_gb"]:
                log(f"=== free memory {free:.1f} GB is below {s3d.D3['min_free_gb']} GB: stopped before ({part}); run again to go on ===")
                break
            log(f"=== (A7 {part}) (free memory {free:.1f} GB) {time.strftime('%H:%M:%S')} ===")
            before = s3b.programs()
            rc = s3d.subprocess.run([sys.executable, "-u", "-m", "research.waa_ksurface.a7_real_fit", "--one", part]).returncode
            if rc != 0:
                cur = s3b.programs()
                grew = sorted(((cur[n] - before.get(n, 0.0), n) for n in cur), reverse=True)[:5]
                log(f"=== (A7 {part}) died with exit code {rc}; grew most during it: " + "; ".join(f"{n} {g:+.2f} GB" for g, n in grew) + "; run again to go on ===")
                break
    fit = pickle.loads(fit_file.read_bytes()) if fit_file.exists() else None
    boots = {e: pickle.loads((CACHE / f"boot_{e}.pkl").read_bytes()) for e in ESTIMATES if (CACHE / f"boot_{e}.pkl").exists()}
    report(fl, fit, boots)
    log(f"=== A7 ended {time.strftime('%Y-%m-%d %H:%M:%S')} ===")


if __name__ == "__main__":
    main()
