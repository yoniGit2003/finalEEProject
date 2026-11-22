"""Offset check (decisions.md, "Offset check", logged before the run): can a plain exponential drying, seen through the
measured baseline offsets and the 0.3 dB grid, make the step 8 switch (drying stops below W_half ~ 0.3 dB)?

    python -m research.waa_ksurface.a8o_offset              every seed in order (a finished seed is read from the cache)
    python -m research.waa_ksurface.a8o_offset --report     write the report from the cache
    python -m research.waa_ksurface.a8o_offset --one r0     one seed in this process (used by the main run)

Truth: K0 / K0, 1 / k_dry 9.7 min, 1 / k_wet 10.92 min, the step 3c fake world, with the normal offset error replaced
by offsets drawn from the MEASURED per-link-event offsets of the learning events. Per seed: the full A7 procedure (step
3c), then the step 8 form plugged back at the radar-input estimate of k_wet. Fake data only; VAL and TEST not opened.
"""
import pickle
import sys
import time

import numpy as np
import pandas as pd
import psutil

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a6_recovery as a6, a8_forms as F8, data, fake6, ksurf, s1_separability as s1
from research.waa_ksurface import s3_recovery as s3, s3b_recovery as s3b, s3c_a7_variant as s3c, s3d_kwet_fixed as s3d
from research.waa_ksurface import paths as P
from research.wet_antenna.system_id import PADM

CO = data.branch_config()["a8o"]
GB = s3b.GB
log = a6.log
CACHE = P.CACHE / "a8o_offset"
NAMES = [f"r{i}" for i in range(CO["replicates"])]
SPEC = {"kind": "const", "k_wet": 1.0 / CO["tau_wet_min"], "k_dry": 1.0 / CO["tau_dry_min"]}


# =================================================================================================== the measured offsets
def measured_offsets(S):
    """The E3 offset of every learning link-event: mean of A - c over the counted dry 5 min bins of the 30 min before the
    first wet bin (at least 3 such bins). Pre-rain windows only."""
    D, le = S["D"], S["le"]
    out = []
    for e in np.nonzero(S["learn"])[0]:
        o5, first = D["le_boff5"][e], D["le_first"][e]
        sl = slice(o5 + max(first - 6, 0), o5 + first)
        dry = (~D["wet5"][sl]) & (D["w5"][sl] > 0)
        if dry.sum() >= 3:
            out.append(float(np.mean(D["y5"][sl][dry]) - D["le_c"][e]))
    return np.array(out)


def add_offsets(Dk, le_n, off):
    """Add one offset per link-event to every minute of it (after the grid: equal in law to adding it before, since the
    grid phase is uniform and independent; see a8o_test)."""
    Dk["y1"] = Dk["y1"] + np.repeat(off, le_n)
    Dk["y5"] = Dk["y1"].reshape(-1, 5).mean(axis=1)
    return Dk


def seed_of(name):
    return CO["seed"] + NAMES.index(name)


# =================================================================================================== one seed
def fake(S, seed, pool):
    D, le, learn = S["D"], S["le"], S["learn"]
    sc = np.repeat(np.where(learn, 1.0, s1.S1["radar_scale_added"]), D["le_n"] + 2 * PADM)     # the fake world as in step 3c
    Dt = dict(D)
    Dt["rain"] = np.asarray(D["rain"], float) / sc
    keep = fake6.A6["noise"]
    fake6.A6["noise"] = {**keep, "offset_sd_db": 0.0}           # the normal offset off; the same random numbers are drawn
    try:
        Dk, ex = fake6.make(Dt, le, S["n_cells"], SPEC, S["curves"], S["iF"], S["iP"], seed, S["step"], S["dist"],
                            ghost_pool=np.sort(le.event_id[learn].unique()))
    finally:
        fake6.A6["noise"] = keep
    Dk["rain"] = Dk["rain"] * sc
    off = np.random.default_rng([seed, 1]).choice(pool, len(le), replace=True)
    return add_offsets(Dk, D["le_n"], off), ex, off


def w_half(theta):
    c0, c1, c2 = theta
    if c1 < 0 and c2 < 0:
        return float(np.log(-c1) / c0)
    if c1 > 0 and c2 > 0:
        return float(np.log(c1) / c0)
    return np.nan


def run_one(name, S):
    seed = seed_of(name)
    t0 = time.time()
    pool = measured_offsets(S)
    Dk, ex, off = fake(S, seed, pool)
    bi = ex["band_idx"]
    fl = a6.tail_flags(S, Dk)                                   # the flag, before any fit (as step 3c)
    flagged = set(fl.event_id[fl.flagged])
    use = (S["learn"] & ~S["le"].event_id.isin(flagged).to_numpy()).astype(float)
    res, (high, used) = s3c.procedure(S, Dk, bi, S["curves"], use, f"({name}) ")
    dens = a6.state_density(S, Dk, ex)
    out = {"name": name, "seed": seed, "spec": SPEC, "flags": fl, "events_flagged": int(fl.flagged.sum()), "events_checked": int(len(fl)),
           "offsets": {"n_pool": len(pool), "drawn_mean": float(off.mean()), "drawn_sd": float(off.std(ddof=1))}}
    out["cv"] = pd.concat([res["wet"]["table"], res["dry"]["table"]], ignore_index=True)
    out["tau_wet"] = res.get("tau_wet")
    out["judge"], _ = a6.judge(S, SPEC, res, dens, True)
    out["procedure_minutes"] = (time.time() - t0) / 60
    log(f"   ({name}) procedure: drying {s3b.drying(out).selected_rung}; {s3b.recovery_text(out)} ({out['procedure_minutes']:.1f} min)")

    # the switch fit at the radar-input estimate (step 8's plug-back, unchanged)
    t1 = time.time()
    counted, _, _, _ = s1.score_sets(S, Dk)
    prob = ksurf.prepare(Dk, bi, S["curves"], Dk["rain"], 0)
    prob["w5"] = (counted & used).astype(float)
    tau_r = res["tau_wet"]["radar"]
    rk = pickle.loads((P.CACHE / "a8_name" / "rank.pkl").read_bytes())     # step 8's chosen form and its curve fit (read only)
    fit8 = rk["fits"][CO["form_name"]]
    form = F8.from_sympy(fit8["spec"]["key"], fit8["spec"]["expr"])
    fm = F8.FormModel(prob, form, (np.array([np.log(1.0 / tau_r)]), np.zeros(0)))
    m2, xa = res["dry"]["model"], res["x_alt"]
    z0 = np.r_[np.asarray(fit8["theta"]["radar"], float), xa[m2.sl["F"]], xa[m2.sl["P"]]]
    z = fm.fit(z0, use)
    th = z[:form.p]
    out["switch"] = {"theta": th, "w_half_db": w_half(th), "tau_inf_min": float(th[1] / th[2]) if th[2] != 0 else np.nan,
                     "converged": fm.last["converged"], "iterations": fm.last["nit"], "loss": fm.last["loss"], "tau_wet_radar": tau_r,
                     "minutes": (time.time() - t1) / 60}
    log(f"   ({name}) switch fit: theta {np.round(th, 5).tolist()}; W_half {out['switch']['w_half_db']:.3f} dB; 1/k_inf "
        f"{out['switch']['tau_inf_min']:.2f} min; converged {fm.last['converged']} ({out['switch']['minutes']:.1f} min)")
    out["minutes"] = (time.time() - t0) / 60
    return out


# =================================================================================================== reading and report
def reading(R):
    rungs = [s3b.drying(r).selected_rung for r in R.values()]
    n1, n0 = rungs.count("K1"), rungs.count("K0")
    wh = [r["switch"]["w_half_db"] for r in R.values() if s3b.drying(r).selected_rung == "K1"]
    med = float(np.nanmedian(wh)) if len(wh) and np.isfinite(wh).any() else np.nan
    lo, hi = CO["w_half_range_db"]
    if len(R) < len(NAMES):
        return "incomplete", n1, n0, med
    if n1 >= CO["k1_min"] and lo <= med <= hi:
        return "consistent with an instrument floor (R7)", n1, n0, med
    if n0 >= CO["k0_min"]:
        return "the non-exponential claim stands", n1, n0, med
    return "undecided", n1, n0, med


def report(R, state):
    rows = []
    for n, r in R.items():
        d = s3b.drying(r)
        g, lo, hi = s3.k1_gain(r["cv"])
        sw = r["switch"]
        rows.append({"seed name": n, "seed": r["seed"], "drying rung": d.selected_rung, "K1 gain on L-dryhigh (dB)": round(g, 5),
                     "interval": f"{lo:+.5f} to {hi:+.5f}", "1/k_dry gauge / radar (min)": f"{d.tau_fit_min:.2f} / {d.tau_fit_alt_min:.2f}",
                     "1/k_wet gauge / radar (min)": f"{r['tau_wet']['gauge']:.2f} / {r['tau_wet']['radar']:.2f}",
                     "events flagged": f"{r['events_flagged']} of {r['events_checked']}", "offsets drawn mean / sd (dB)":
                     f"{r['offsets']['drawn_mean']:+.3f} / {r['offsets']['drawn_sd']:.3f}",
                     "switch c0, c1, c2": ", ".join(f"{v:.4g}" for v in sw["theta"]), "W_half (dB)": round(sw["w_half_db"], 3),
                     "1/k_inf (min)": round(sw["tau_inf_min"], 2), "switch converged": sw["converged"], "minutes": round(r["minutes"], 1)})
    t = pd.DataFrame(rows)
    save(t, P.A8O_RUNS, index=False)
    verdict, n1, n0, med = reading(R)
    pool = measured_offsets(s3.setting()) if R else np.zeros(1)
    L = ["# Offset check: plain exponential drying, measured offsets, the 0.3 dB grid\n",
         "Written by `research/waa_ksurface/a8o_offset.py` (decisions.md, \"Offset check\", logged before the run). Fake data only; VAL and TEST were "
         "not opened.\n",
         f"Truth K0 / K0: 1 / k_dry {CO['tau_dry_min']} min, 1 / k_wet {CO['tau_wet_min']} min; the step 3c fake world with the normal offset replaced by "
         f"offsets drawn from the {len(pool)} measured learning-set offsets (sd {pool.std(ddof=1):.3f} dB, mean {pool.mean():+.3f} dB, 10-90 % "
         f"{np.percentile(pool, 10):+.3f} to {np.percentile(pool, 90):+.3f} dB, 1-99 % {np.percentile(pool, 1):+.3f} to {np.percentile(pool, 99):+.3f} dB). "
         "Full A7 procedure (step 3c); then the step 8 form k = c2/(c1 - exp(c0 W)) plugged back at the radar-input estimate of k_wet.\n",
         md_table(t), "",
         f"**Reading: {verdict}.** K1 in {n1} of {len(R)}, K0 in {n0} of {len(R)}; median W_half over the K1 seeds: "
         f"{'n/a' if not np.isfinite(med) else f'{med:.3f} dB'} (floor if K1 >= {CO['k1_min']} and median W_half in {CO['w_half_range_db']}; "
         f"stands if K0 >= {CO['k0_min']}).\n"]
    if state:
        L.append(f"Lowest free memory {state['min_free']:.2f} GB; low-memory events {len(state['low'])}." +
                 (f" **NOT COMPLETE: {state['incomplete']}**" if state.get("incomplete") else "") + "\n")
    save("\n".join(L) + "\n", P.A8O_REPORT)
    log(f"report: {rel(P.A8O_REPORT)}")


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    if "--one" in sys.argv:
        name = sys.argv[sys.argv.index("--one") + 1]
        r = run_one(name, s3.setting())
        r["peak_gb"] = psutil.Process().memory_info().peak_wset / GB
        (CACHE / f"{name}.pkl").write_bytes(pickle.dumps(r))
        log(f"=== (a8o {name}) done: {r['minutes']:.1f} min; peak memory {r['peak_gb']:.2f} GB ===")
        return
    if "--report" in sys.argv:
        report({n: pickle.loads((CACHE / f"{n}.pkl").read_bytes()) for n in NAMES if (CACHE / f"{n}.pkl").exists()}, None)
        return
    state = s3d.start_state("offset check")
    R, stop = s3d.run_children("research.waa_ksurface.a8o_offset", CACHE, NAMES, state, lambda n: f"seed {seed_of(n)}")
    report(R, state)
    log(f"=== offset check ended {time.strftime('%Y-%m-%d %H:%M:%S')} ===")


if __name__ == "__main__":
    main()
