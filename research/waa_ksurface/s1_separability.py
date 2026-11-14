"""Step S1 of the k-surface branch: on which held-out loss does a real k_dry(W) gain more than noise does?

    python -m research.waa_ksurface.s1_separability             run what is missing (finished replicates are read from the cache)
    python -m research.waa_ksurface.s1_separability --report    write the report from the cache

Fake data only (the fake world of A6, library `fake6.py`), fresh seeds. Truth (a) has a constant k_dry, truth (c) a k_dry
that triples from W = 0 to 2 dB. On five replicates of each, drying K0 and drying K1 are fitted as in stage 2 of
`ladder.py`, and the held-out gain of K1 over K0 is scored on three sets of 5 min values:

    L-all       every counted held-out value (the current rule)
    L-dry       from the bin after the last wet bin of a link-event to 120 min after it, where the radar input is below 0.1 mm/h
    L-dryhigh   L-dry, only the link-events whose measured A - c in the first dry bin is at least 1 dB

and in two variants of the learning set:

    current     the 14 TRAIN events with a gauge ratio
    extended    plus the L-dry values of the other TRAIN events, whose radar reads 2.5 times the true rain in the fake world

A (loss, variant) is separable if the smallest gain on (c) is above the largest gain on (a). If none is, truth (c) is run
at k ratios 5 and 10: the detection limit. Definitions and the adoption rule: decisions.md, "S1" (written before this code).
"""
import pickle
import sys
import time

import numpy as np
import pandas as pd

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import data, fake6, fp, ksurf, ladder
from research.waa_ksurface import paths as P
from research.wet_antenna.system_id import PADM

A5, A6, S1 = ksurf.config(), data.branch_config()["a6"], data.branch_config()["s1"]
LOSSES = ["L-all", "L-dry", "L-dryhigh"]
VARIANTS = ["current", "extended"]


def log(msg):
    print(msg, flush=True)


def huber(r, delta):
    a = np.abs(r)
    return np.where(a <= delta, 0.5 * r * r, delta * (a - 0.5 * delta))


# =================================================================================================== the setting
def setting():
    """All TRAIN link-events of the branch (28 and 38 GHz). learn = the link-events of the events with a gauge ratio (the
    learning set); the others are the added events, whose radar is not corrected (g_e = 1)."""
    d = data.open_release()
    _, D, le, _ = data.branch_train(d)
    g = data.g_table()
    learn = le.event_id.isin(set(g[g.path_used].event_id)).to_numpy()
    gz = d.table("gauges").set_index("gauge_id")
    curves, iF, iP = fake6.fake_curves(fp.measured())
    n5 = D["le_n"] // 5
    bin_le = np.repeat(np.arange(len(le)), n5)
    return {"D": D, "le": le, "learn": learn, "curves": curves, "iF": iF, "iP": iP,
            "n_cells": np.maximum(le.cml_id.map(d.links().set_index("cml_id").n_cells).to_numpy(float), 1.0),
            "step": le.gauge_id.map(gz.quantization_mm).to_numpy(float), "dist": le.gauge_distance_m.to_numpy(float) / 1000.0,
            "gauge": (le.gauge_id != "").to_numpy() & learn, "bin_le": bin_le, "pos": np.arange(len(bin_le)) - D["le_boff5"][bin_le]}


def spec_of(truth):
    """The truth of a name: "a", "c", or "c_r5" = the k_dry(W) family at the ratio k(W = 2 dB) / k(W = 0) = 5."""
    if truth in ("a", "c"):
        return A6["truths"][truth]
    a = (float(truth.split("_r")[1]) - 1) / 2                   # ratio = 1 + 2 a
    return {"kind": "dry_w", "k_wet": 0.2, "k_dry0": (1 / 15) / (1 + a), "w_scale_db": 1 / a}


def fake(S, truth, rep):
    """One fake dataset on all TRAIN events. The true rain of an added link-event is its radar value / 2.5; the fit sees
    the radar value (times the T1 error), so the radar of the added events reads 2.5 times too much."""
    D, le = S["D"], S["le"]
    sc = np.repeat(np.where(S["learn"], 1.0, S1["radar_scale_added"]), D["le_n"] + 2 * PADM)
    Dt = dict(D)
    Dt["rain"] = np.asarray(D["rain"], float) / sc
    seed = S1["seed"] + 10 * S1["truths"].index(truth) + rep
    Dk, ex = fake6.make(Dt, le, S["n_cells"], spec_of(truth), S["curves"], S["iF"], S["iP"], seed, S["step"], S["dist"])
    Dk["rain"] = Dk["rain"] * sc
    return Dk, ex, seed


def score_sets(S, Dk):
    """The 5 min values of L-dry and L-dryhigh (boolean per bin), and the level (dB) of every link-event at the end of its rain."""
    D = S["D"]
    n5 = D["le_n"] // 5
    last = D["le_last"]
    counted = D["w5"] > 0
    rain5 = fake6.window(D, Dk["rain"]).reshape(-1, 5).mean(axis=1)
    lb = last[S["bin_le"]]
    dry = counted & (S["pos"] > lb) & (S["pos"] <= lb + S1["dry_bins"]) & (rain5 < S1["dry_below_mmh"])
    ok = (last + 1 >= 0) & (last + 1 < n5)
    first = np.where(ok, D["le_boff5"] + last + 1, 0)
    ok &= dry[first]
    level = np.where(ok, Dk["y5"][first], np.nan)
    high = dry & (np.nan_to_num(level, nan=-1.0) >= S1["high_db"])[S["bin_le"]]
    return counted, dry, high, level


# =================================================================================================== the fits
def held_out(m, x, use, fold, target, delta):
    """The Huber loss of every 5 min value when its event is held out (NaN for the values of link-events that are not used)."""
    H = np.full(len(target), np.nan)
    le_of = m.p["bin_le"]
    for f in range(fold.max() + 1):
        xf = m.fit(x, use * (fold != f), maxiter=A5["maxiter_fold"])
        out5 = m.simulate(xf)[1]
        b = ((fold == f) & (use > 0))[le_of]
        H[b] = huber(out5[b] - target[b], delta)
    return H


def one(S, truth, rep):
    """One replicate: k_wet from the gauge-input set, then drying K0 and K1 with the 5-fold CV in both variants.
    Returns the rows of the gains (one per variant and loss)."""
    t0 = time.time()
    D, le, learn = S["D"], S["le"], S["learn"]
    Dk, ex, seed = fake(S, truth, rep)
    bi = ex["band_idx"]
    E = A6["extension"]
    start_F = [np.full(len(i), E["f_start_inc_db"]) for i in S["iF"]]
    start_P = [E["p_start_share"] * i for i in S["iP"]]
    lam0 = A5["lambdas"][0]
    # ---- k_wet: the fit the wetting ladder returned in A6 final (wetting K0, drying K1 as a nuisance, gauge input, learning set)
    rain_g, ok5 = ladder.gauge_input(Dk)
    prob_g = ksurf.prepare(Dk, bi, S["curves"], rain_g, 0)
    prob_g["w5"] = prob_g["w5"] * ok5
    m1 = ksurf.Model(prob_g, "K0", "K1", lam0)
    x1 = m1.fit(m1.start(start_F, start_P), S["gauge"].astype(float))
    ext = m1.parts(x1)[4:]
    tau_wet = float(1 / m1.k(x1, "wetting", 0, 0))
    counted, dry, high, level = score_sets(S, Dk)
    rows = []
    for variant in VARIANTS:
        use = learn.astype(float) if variant == "current" else np.ones(len(le))
        used = (use > 0)[S["bin_le"]]
        w5 = np.where(learn[S["bin_le"]], counted, dry) & used      # of an added link-event only its L-dry values count
        prob = ksurf.prepare(Dk, bi, S["curves"], Dk["rain"], 0)
        prob["w5"] = w5.astype(float)
        prob["bin_le"] = S["bin_le"]
        fold = ladder.folds(le, use)
        sets = {"L-all": w5, "L-dry": dry & used, "L-dryhigh": high & used}
        H, info, prev = {}, {}, None
        for rung, lams in (("K0", [lam0]), ("K1", A5["lambdas"])):
            for lam in lams:
                m = ksurf.Model(prob, "K0", rung, lam)
                x0 = m.start([e.copy() for e in ext[0]], [e.copy() for e in ext[1]]) if prev is None else m.carry(*prev)
                x0[m.sl["gw"]] = x1[m1.sl["gw"]]
                m.hold("gw", "hw")
                x = m.fit(x0, use)
                H[(rung, lam)] = held_out(m, x, use, fold, prob["target"], A5["huber_delta_db"])
                info[(rung, lam)] = (m, x)
            if rung == "K0":
                prev = info[("K0", lam0)]
        m0, x0_ = info[("K0", lam0)]
        for loss, mask in sets.items():
            e0 = float(np.sqrt(2 * H[("K0", lam0)][mask].mean()))
            e1 = {lam: float(np.sqrt(2 * H[("K1", lam)][mask].mean())) for lam in A5["lambdas"]}
            lam = min(e1, key=e1.get)
            mk, xk = info[("K1", lam)]
            seg = np.unique(S["bin_le"][mask])
            rows.append({"truth": truth, "replicate": rep, "seed": seed, "variant": variant, "loss": loss, "values": int(mask.sum()), "link_events": len(seg),
                         "events": int(le.event_id.iloc[seg].nunique()), "error_K0_db": e0, "error_K1_db": e1[lam], "lambda_K1": lam, "gain_db": e0 - e1[lam],
                         "tau_wet_min": tau_wet, "tau_dry_K0_min": float(1 / m0.k(x0_, "drying", 0, 0)),
                         **{f"tau_dry_K1_W{w:g}_min": float(1 / mk.k(xk, "drying", 0, w)) for w in (0, 1, 2)}})
        log(f"   {truth} r{rep} {variant}: " + ", ".join(f"{r['loss']} {r['gain_db']:+.5f}" for r in rows[-3:]) + f"; 1/k_dry K0 {rows[-1]['tau_dry_K0_min']:.1f} min")
    log(f"   {truth} r{rep}: {(time.time() - t0) / 60:.1f} min")
    return pd.DataFrame(rows)


# =================================================================================================== the rule
def summary(runs, truth_c):
    """Per loss and variant: the gains on truth_c against the gains on (a), the ratio, separable or not."""
    rows = []
    for variant in VARIANTS:
        for loss in LOSSES:
            c = runs[(runs.truth == truth_c) & (runs.variant == variant) & (runs.loss == loss)].gain_db
            a = runs[(runs.truth == "a") & (runs.variant == variant) & (runs.loss == loss)].gain_db
            if len(c) < S1["replicates"] or len(a) < S1["replicates"]:
                continue
            ratio = np.inf if (a.max() <= 0 < c.min()) else (c.min() / a.max() if a.max() > 0 else 0.0)
            rows.append({"truth": truth_c, "loss": loss, "variant": variant, "c_gain_min_db": c.min(), "c_gain_max_db": c.max(), "a_gain_min_db": a.min(), "a_gain_max_db": a.max(),
                         "ratio": ratio, "separable": bool(c.min() > a.max()), "c_replicates_above_max_a": int((c > a.max()).sum())})
    return pd.DataFrame(rows)


def extension_check(runs):
    """The fake check of the extended variant: the K0 drying constant on truth (a), per replicate and variant."""
    t = runs[(runs.truth == "a") & (runs.loss == "L-all")].pivot(index="replicate", columns="variant", values="tau_dry_K0_min").reset_index()
    t["extended_within_25pct"] = (15.0 / t.extended - 1).abs() <= A6["judge"]["within"]
    return t


def adopt(summ, eligible_extended):
    """The eligible combination with the largest ratio, if that ratio is above 1 (ties at infinity: the larger smallest (c) gain)."""
    s = summ[summ.separable & ((summ.variant == "current") | eligible_extended)]
    if not len(s):
        return None
    return s.sort_values(["ratio", "c_gain_min_db"], ascending=False).iloc[0]


def main():
    only_report = "--report" in sys.argv
    cache = P.CACHE / "s1_separability"
    cache.mkdir(parents=True, exist_ok=True)
    S = None
    got = {}

    def run(truths):
        nonlocal S
        for truth in truths:
            for rep in range(S1["replicates"]):
                f = cache / f"{truth}_r{rep}.pkl"
                if not f.exists():
                    if only_report:
                        continue
                    S = S or setting()
                    log(f"=== truth ({truth}), replicate {rep} ===")
                    f.write_bytes(pickle.dumps(one(S, truth, rep)))
                got[(truth, rep)] = pickle.loads(f.read_bytes())
        return pd.concat(got.values(), ignore_index=True) if got else pd.DataFrame()

    runs = run(["a", "c"])
    summ = summary(runs, "c")
    chk = extension_check(runs)
    eligible = bool(chk.extended_within_25pct.all()) and len(chk) == S1["replicates"]
    best = adopt(summ, eligible) if len(summ) else None
    limit = None
    if best is None and len(summ):                              # nothing is separable: the detection limit
        runs = run(["a", "c", "c_r5", "c_r10"])
        limit = pd.concat([summary(runs, t) for t in ("c_r5", "c_r10")], ignore_index=True)
    save(runs.round(6), P.S1_RUNS, index=False)
    save(summ.round(6), P.S1_TABLE, index=False)
    save(chk.round(4), P.S1_CHECK, index=False)

    cov = runs[(runs.truth == "c") & (runs.loss == "L-dryhigh")].groupby("variant")[["values", "link_events", "events"]].mean().round(0).reset_index()
    cov_dry = runs[(runs.truth == "c") & (runs.loss == "L-dry")].groupby("variant")[["values", "link_events", "events"]].mean().round(0).reset_index()
    show = summ.drop(columns="truth").copy()
    for c in ("c_gain_min_db", "c_gain_max_db", "a_gain_min_db", "a_gain_max_db"):
        show[c] = show[c].round(5)
    show["ratio"] = show.ratio.round(2)
    a = ["# S1: signal against spurious gain under targeted losses\n",
         "Written by `research/waa_ksurface/s1_separability.py`. Fake data only (the fake world of A6 on all TRAIN events of the branch, fresh seeds from "
         f"{S1['seed']}), {S1['replicates']} replicates per truth. Truth (a): constant k_dry (1 / k = 15 min). Truth (c): k_dry = (1 + W) / 30 per min (k triples from W = 0 to 2 dB). "
         "Gain = held-out error of drying K0 minus that of drying K1, error = sqrt(2 x mean Huber loss) in dB, K1 at the lambda with the lowest held-out error on the same set. "
         "The fits are trained on the loss of A5 in every case; only the set of held-out values that is scored differs. In the extended variant the radar of the added events "
         f"reads {S1['radar_scale_added']} times the true rain. No real attenuation is used. VAL and TEST were not opened.\n",
         "## The table\n", md_table(show),
         f"\nSeparable = the smallest gain on (c) is above the largest gain on (a). Ratio = smallest (c) gain / largest (a) gain. (`{rel(P.S1_TABLE)}`; every replicate: `{rel(P.S1_RUNS)}`)\n",
         "## The fake check of the extended variant (radar of the added events 2.5 times too high)\n",
         "1 / k_dry of the K0 fit on truth (a), minutes (truth 15; within 25 % = 12.0 to 20.0):\n", md_table(chk.round(2)),
         f"\nThe extended variant is **{'eligible' if eligible else 'NOT eligible'}** for adoption (rule: within 25 % in 5 of 5). (`{rel(P.S1_CHECK)}`)\n",
         "Coverage (fake world, mean of the (c) replicates). L-dryhigh (level at the end of the rain at least 1 dB):\n", md_table(cov),
         "\nL-dry:\n", md_table(cov_dry), ""]
    a.append("## Verdict\n")
    if best is not None:
        a.append(f"**Adopted: {best.loss}, {best.variant} variant** (ratio {best.ratio:.2f}). Margin of Deviation 9 = the largest (a) gain on it = "
                 f"**{max(best.a_gain_max_db, 0.0):.5f} dB**.\n")
    else:
        a.append("**Nothing is separable: no loss or variant is adopted.** Deviation 9 does not come into force.\n")
    if limit is not None and len(limit):
        save(limit.round(6), P.S1_LIMIT, index=False)
        lim = limit.copy()
        lim["k_ratio"] = lim.truth.map({"c_r5": 5, "c_r10": 10})
        for c in ("c_gain_min_db", "c_gain_max_db", "a_gain_min_db", "a_gain_max_db"):
            lim[c] = lim[c].round(5)
        lim["ratio"] = lim.ratio.round(2)
        a += ["## The detection limit: truth (c) at larger sizes\n",
              "k_dry = (1/15)(1 + a W) / (1 + a), so 1 / k_dry = 15 min at W = 1 dB; `k_ratio` = k(W = 2 dB) / k(W = 0). Against the same (a) gains.\n",
              md_table(lim[["k_ratio", "loss", "variant", "c_gain_min_db", "c_gain_max_db", "a_gain_max_db", "ratio", "separable", "c_replicates_above_max_a"]]),
              f"\n(`{rel(P.S1_LIMIT)}`)\n"]
    a += ["## Every replicate\n", md_table(runs[["truth", "replicate", "variant", "loss", "values", "link_events", "error_K0_db", "error_K1_db", "lambda_K1", "gain_db",
                                                 "tau_wet_min", "tau_dry_K0_min", "tau_dry_K1_W0_min", "tau_dry_K1_W1_min", "tau_dry_K1_W2_min"]].round(5)), ""]
    save("\n".join(a) + "\n", P.S1_REPORT)
    log(f"report: {rel(P.S1_REPORT)}")


if __name__ == "__main__":
    main()
