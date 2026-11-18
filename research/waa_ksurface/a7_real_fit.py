"""A7: the real fit, exactly as pre-registered (decisions.md, "A7: pre-registration of the real fit", commit bbb5d2f).

    python -m research.waa_ksurface.a7_real_fit --flag         step 1: the tail flag on the real learning events (no fit)
    python -m research.waa_ksurface.a7_real_fit                step 2: the fit at 1 / k_wet = 2.5, 5, 10 min, then the bootstrap
    python -m research.waa_ksurface.a7_real_fit --report       write the report from the cache
    python -m research.waa_ksurface.a7_real_fit --one fit      one part in this process (used by the main run: fit, boot_2.5, ...)

What it does:
    1. the flag      the tail flag of Deviation 12 on the REAL attenuation of the 14 learning events; the events it removes
                     take no part in anything below. Its count is the first number of A7, saved before any fit.
    2. the fit       the procedure of step 3c / 3d (`s3d_kwet_fixed.fit_fixed_kwet`) on the real TRAIN data of the unflagged
                     learning events (rain = radar / g_e, the data layer of this branch): stage 1 (gauge input, wetting K0)
                     once, then drying K0 against K1 on L-dryhigh (margin 0.00291 dB) with 1 / k_wet held at 2.5, 5, 10 min.
                     F and P: the measured curves (E4) fixed in their measured ranges, free monotone extensions above
                     (started at the model defaults of A5).
    3. bootstrap     per value of 1 / k_wet, the selected rung refitted on 200 event-bootstrap draws (seed 3900, the same
                     draws at every value), 1 / k_dry at W = 0, 0.5, 1, 2 dB; interval 5-95 %; the reported interval is
                     the union over the three values (R3).
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
TAUS = s3d.D3["tau_wet_min"]
GB = s3b.GB
log = a6.log
CACHE = P.CACHE / "a7_real_fit"
K0_TEXT = ("no dependence of the drying rate on W larger than a ratio k(W = 2 dB) / k(W = 0) of 2 is present (2 = the smallest ratio "
           "detected for this variant in step 3c; 1.5 was detected in 2 of 5 and is not excluded)")
K1_TEXT = ("non-exponential drying: the drying rate falls as the antenna dries; the mechanism (W-dependent rate or a second, slow store) is "
           "not claimed (R5); a spread of the baseline between links is excluded as its cause at the measured size (0.088 dB, step 3c (a) PASS, "
           "also at 0.15 dB)")
ND_TEXT = "rung not determined over the k_wet range: no rung reading applies"


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
def fit_part(S, curves, band_idx, use):
    D = S["D"]
    used = (use > 0)[S["bin_le"]]
    counted, dry, high, _ = s1.score_sets(S, D)
    stage2 = {"w5": counted & used, "score": high & used, "margin_db": s3c.margin(), "selector": s3.select_k0_k1, "rungs": s3.DRY_RUNGS}
    fits = s3d.fit_fixed_kwet(D, S["le"], band_idx, curves, S["gauge"].astype(float), use, stage2, None, None, tag="(A7) ")
    out = {}
    for tau, res in fits.items():
        st = res["dry"]
        g, lo, hi = s3.k1_gain(st["table"])
        m, x = res["model"], res["x"]
        out[tau] = {"rung": st["rung"], "lam": st["lam"], "x": x, "table": st["table"], "gain": g, "lo": lo, "hi": hi,
                    "tau_dry": dict(zip(A7["w_report_db"], s3d.k_dry_at(m, x, A7["w_report_db"]))),
                    "extension": {f"{k[0]} band {k[1]}": v for k, v in m.extension(x).items()}}
    out["scored_values"] = int((high & used).sum())
    out["scored_link_events"] = int(len(np.unique(S["bin_le"][high & used])))
    out["stage2_w5"] = counted & used
    return out


def boot_part(S, curves, band_idx, use, fit, tau):
    """The selected rung at this 1 / k_wet refitted on event-bootstrap draws (k_wet held, extension free)."""
    D, le = S["D"], S["le"]
    f = fit[tau]
    prob = ksurf.prepare(D, band_idx, curves, D["rain"], 0)
    prob["w5"] = np.asarray(fit["stage2_w5"], float)
    m = ksurf.Model(prob, "K0", f["rung"], f["lam"])
    m.hold("gw", "hw")
    ev_all = le.event_id.to_numpy()
    events = np.sort(np.unique(ev_all[use > 0]))
    rng = np.random.default_rng(A7["boot_seed"])
    rows = []
    t0 = time.time()
    for b in range(A7["boot_draws"]):
        pick = rng.choice(events, len(events), replace=True)
        cnt = pd.Series(pick).value_counts()
        w = use * pd.Series(ev_all).map(cnt).fillna(0).to_numpy()
        xb = m.fit(f["x"].copy(), w)
        rows.append({"draw": b, **{f"1/k_dry at W = {W:g} dB": v for W, v in zip(A7["w_report_db"], s3d.k_dry_at(m, xb, A7["w_report_db"]))},
                     "converged": m.last["converged"]})
        if (b + 1) % 20 == 0:
            log(f"   (A7) [1/k_wet {tau:g}] bootstrap {b + 1} of {A7['boot_draws']} ({(time.time() - t0) / 60:.1f} min)")
    return pd.DataFrame(rows)


# =================================================================================================== readings and report
def reading(fit):
    """(rung, sentence): the pre-registered reading of the three selections."""
    rungs = {fit[t]["rung"] for t in TAUS}
    if len(rungs) > 1:
        return "not determined", ND_TEXT
    rung = rungs.pop()
    return rung, K0_TEXT if rung == "K0" else K1_TEXT


def report(fl, fit, boots):
    L = ["# A7: the real fit (pre-registered, decisions.md; commit bbb5d2f)\n",
         "Written by `research/waa_ksurface/a7_real_fit.py`. Real TRAIN data of this branch (28 and 38 GHz), the 14 learning events (rain = radar / g_e). "
         "VAL and TEST were not read. The procedure is that of steps 3c and 3d: the tail flag first; drying K0 against K1 on L-dryhigh, margin "
         f"{s3c.margin():.5f} dB with a paired 90 % interval that excludes zero; 1 / k_wet held at 2.5, 5, 10 min (R3); F and P measured parts fixed, "
         "extensions free (R1).\n"]
    checked, removed = len(fl), fl[fl.flagged]
    L += ["## 1. The flag on the real learning events (before any fit)\n",
          f"Events checked: {checked} (E20150831_0455 is not checkable). **Removed: {len(removed)}**"
          + (f" ({', '.join(removed.event_id)})" if len(removed) else "") + ".\n", md_table(fl), ""]
    if fit is None:
        save("\n".join(L) + "\n", P.A7_REPORT)
        return
    rows = []
    for t in TAUS:
        f = fit[t]
        b = boots.get(t)
        r = {"1/k_wet (min)": t, "drying rung": f["rung"], "lambda": f["lam"], "K1 gain (dB)": round(f["gain"], 5),
             "interval": f"{f['lo']:+.5f} to {f['hi']:+.5f}", "margin (dB)": round(s3c.margin(), 5)}
        for W in A7["w_report_db"]:
            col = f"1/k_dry at W = {W:g} dB"
            iv = "" if b is None else f" ({np.percentile(b[col], 5):.1f} to {np.percentile(b[col], 95):.1f})"
            r[col + " (min)"] = f"{f['tau_dry'][W]:.1f}{iv}"
        r["bootstrap draws converged"] = "" if b is None else f"{int(b.converged.sum())} of {len(b)}"
        rows.append(r)
    tab = pd.DataFrame(rows)
    save(tab, P.A7_FIT, index=False)
    L += ["## 2. The fit at the three values of 1 / k_wet\n",
          f"Scored values (L-dryhigh): {fit['scored_values']} in {fit['scored_link_events']} link-events. 1 / k_dry with its 5-95 % event-bootstrap "
          f"interval ({A7['boot_draws']} draws, seed {A7['boot_seed']}) in brackets; for K0 the value does not depend on W.\n", md_table(tab), ""]
    if len(boots) == len(TAUS):
        un = []
        for W in A7["w_report_db"]:
            col = f"1/k_dry at W = {W:g} dB"
            lo = min(np.percentile(boots[t][col], 5) for t in TAUS)
            hi = max(np.percentile(boots[t][col], 95) for t in TAUS)
            un.append({"W (dB)": W, "1/k_dry, union of the three intervals (min)": f"{lo:.1f} to {hi:.1f}",
                       "point values at 2.5 / 5 / 10 min": " / ".join(f"{fit[t]['tau_dry'][W]:.1f}" for t in TAUS)})
        L += ["### k_dry with its spread across 1 / k_wet = 2.5-10 min (R3)\n", md_table(pd.DataFrame(un)), ""]
    rung, text = reading(fit)
    L += ["## 3. Result under the pre-registered readings\n", f"Selected drying rung at 2.5 / 5 / 10 min: {' / '.join(fit[t]['rung'] for t in TAUS)}.\n",
          f"**{text}.**\n", "Stated with every A7 result: the flag cannot check E20150831_0455, so unseen water in it, if any, would go undetected.\n"]
    ext = []
    for t in TAUS:
        for k, v in fit[t]["extension"].items():
            ext.append({"1/k_wet (min)": t, "curve and band index": k, "values at the knots (EFFECTIVE, R1)": ", ".join(f"{a:.3f}" for a in np.atleast_1d(v))})
    save(pd.DataFrame(ext), P.A7_EXT, index=False)
    L += ["\n## 4. Extension knots above the measured ranges (labeled EFFECTIVE, R1)\n",
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
            for t in TAUS:
                log(f"=== (A7) 1/k_wet {t:g} min: drying {out[t]['rung']}; K1 gain {out[t]['gain']:+.5f} dB ({out[t]['lo']:+.5f} to {out[t]['hi']:+.5f}); "
                    "1/k_dry " + ", ".join(f"W {W:g}: {v:.1f}" for W, v in out[t]["tau_dry"].items()) + " min ===")
        else:
            tau = float(part.split("_")[1])
            b = boot_part(S, curves, band_idx, use, pickle.loads(fit_file.read_bytes()), tau)
            (CACHE / f"{part}.pkl").write_bytes(pickle.dumps(b))
        log(f"=== (A7 {part}) done: {(time.time() - t0) / 60:.1f} min; peak memory {psutil.Process().memory_info().peak_wset / GB:.2f} GB ===")
        return
    parts = ["fit"] + [f"boot_{t:g}" for t in TAUS]
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
    boots = {t: pickle.loads((CACHE / f"boot_{t:g}.pkl").read_bytes()) for t in TAUS if (CACHE / f"boot_{t:g}.pkl").exists()}
    report(fl, fit, boots)
    log(f"=== A7 ended {time.strftime('%Y-%m-%d %H:%M:%S')} ===")


if __name__ == "__main__":
    main()
