"""wet_antenna step S8b: a POST-FREEZE EVALUATION OF A PRE-FREEZE MODEL - can F be written as the closed form
F_band(R) = A_band (R / 5)^b?

TEST was read once in S8 with the frozen S6 candidate (table F). This step refits nothing on TEST and does not change the S8
verdict. It runs ONE other model, fitted on TRAIN in S6 before TEST was read: the S6 plug-back fit "F and P" (F = A_band
(R / 5)^b, P = c_band R, everything else as S6; the "closed F" of S7).

  Part 1  TRAIN only: the closed-F model's cold event bootstrap (40 draws) and the S6 identifiability check (3 fake
          replicates); the closed F against the table F per band.
  Part 2  TEST: the closed-F model forward on the S8 TEST link-events (the S8 cache, V1 with the fixed g_e), P1-P5 exactly as
          S8, the paired RMSE difference to the table F, the dry-path errors; the decision rule written before this run.

Settings and the rule: config `wet_antenna.test_s8b` (written before the run). Caches: research/.cache/wet_antenna/test_s8b/.
Run from the repository root:
    python -m research.wet_antenna.test_s8b
"""
import hashlib
import json
import os
import pickle
import sys
import time
import zlib
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cmlpipe.api import CMLData
from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import length_test as LTL
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.event_factors import rescaled
from research.wet_antenna.final_model import CFG, DATA, FIX, FK, M, WA, Curves, load_sets, plan_start, predict, s6_models, through
from research.wet_antenna.tails import itu_p838
from research.wet_antenna.validation import LB, boot_items, ci, comp_setup, dry_windows, ev_draws, need, radar_max, w_minutes

S, S8 = WA["test_s8b"], WA["test_s8"]
CACHE = P.S8B_CACHE
K5 = int(np.nonzero(np.asarray(M["knots_mmh"], float)[1:] == 5.0)[0][0]) + 1
KN = [0.5, 1, 2, 5, 10, 20, 40]
STATUS = ("**Status.** TEST was read once in S8 with the frozen S6 candidate (the table F). This step does NOT refit anything on TEST and does "
          "not change the S8 verdict. It runs ONE other model that was fitted on TRAIN in S6, before TEST was read: the S6 plug-back fit "
          "\"F and P\" (F = A_band (R / 5)^b, P = c_band R, everything else as S6; the \"closed F\" of S7). Everything here is a **post-freeze "
          "evaluation of a pre-freeze model**.")


# =================================================================================================== TRAIN jobs (cached)
_W = {}


def closed_model(bands):
    if "C" not in _W:
        _W["Dt1"] = sid.load_data(P.S5A_CACHE / "data_train_v1.npz", bands)
        _W["let"] = pickle.loads((P.SID3_CACHE / "le_train_all.pkl").read_bytes())[0].reset_index(drop=True)
        mods, _ = s6_models(_W["Dt1"], _W["let"], bands)
        _W["C"], _W["spec"] = mods["closed F"]["C"], mods["closed F"]["spec"]
    return _W["Dt1"], _W["let"], _W["C"], _W["spec"]


def job(jb):
    file = CACHE / f"{jb['name']}.pkl"
    if file.exists():
        return jb["name"], pickle.loads(file.read_bytes())
    t0 = time.time()
    D, le, C, spec = closed_model(jb["bands"])
    if jb["kind"] == "fake":
        D0 = dict(D)
        D0["le_c"] = np.zeros(len(le))
        tr = sid.Problem2(D0, "M4b", C.M, agg=1).predict(np.asarray(jb["x_truth"]), C.d)
        Df = sid.fake_from_truth(D, le, tr, jb["seed"], sid.noise_table(), FK, DATA["offset_min_minutes"])
        Cv = Curves(Df, spec, C.x0, C.d, Mx=C.M)
    else:
        n_ev = int(D["le_event"].max()) + 1
        rng = np.random.default_rng(jb["seed"])
        mults = [np.bincount(rng.integers(0, n_ev, n_ev), minlength=n_ev).astype(float) for _ in range(jb["i"] + 1)]
        Cv = Curves(D, spec, C.x0, C.d, Mx=C.M, wle=mults[jb["i"]][D["le_event"]])
    z = np.asarray(jb["z0"])
    nit = 0
    for cap in (1000, 1000):
        f = Cv.fit(z, cap)
        z, nit = f["z"], nit + f["nit"]
    out = {"z": z, "loss": f["loss"], "converged": f["converged"], "nit": nit, "seconds": time.time() - t0}
    file.write_bytes(pickle.dumps(out))
    return jb["name"], out


def run_pool(jobs, workers, clock, results):
    jobs = [j for j in jobs if j["name"] not in results]
    if not jobs:
        return results
    with ProcessPoolExecutor(workers) as ex:
        running = {ex.submit(job, j): j for j in jobs[:workers]}
        queue = jobs[workers:]
        while running:
            done_, _ = wait(running, return_when=FIRST_COMPLETED)
            for fu in done_:
                running.pop(fu)
                name, out = fu.result()
                results[name] = out
                print(f"   {name} done ({(time.time() - clock) / 60:.1f} min)", flush=True)
                if queue:
                    j = queue.pop(0)
                    running[ex.submit(job, j)] = j
    return results


def readout(C, z, bands):
    nb = len(bands)
    _, _, a = C.values(z)
    out = {f"A[{b}]": float(z[i]) for i, b in enumerate(bands)}
    out["b"] = float(z[nb])
    out.update({f"c[{b}]": float(z[nb + 1 + i] / 5) for i, b in enumerate(bands)})
    la = np.log(a)
    out.update({"a_p10": float(np.exp(np.percentile(la, 10))), "a_p90": float(np.exp(np.percentile(la, 90)))})
    return out


def closed_curve(C, z, bands, R):
    g = C.spec["F"][1]["g"]
    nb = len(bands)
    return np.array([z[i] * np.asarray(g(np.asarray(R, float), z[nb]), float) for i in range(nb)])


# =================================================================================================== TEST helpers
def b_models(Dt1, let, bands):
    """The S8 comparison fits B0-B2 (S8 cache, read only) as fz."""
    ids = LTL.link_ids(Dt1, let)
    out = {}
    for b in ("B0", "B1", "B2"):
        spec, Mx = comp_setup(b, len(bands))
        Pm = sid.Problem2(Dt1, "M1", Mx)
        Cv = Curves(Dt1, spec, np.zeros(len(Pm.names)), S8["comparison"]["lag_min"], Mx=Mx, Pm=Pm)
        z = np.asarray(pickle.loads((P.S8_CACHE / f"fit_{b}.pkl").read_bytes())["z"])
        _, _, a = Cv.values(z)
        out[b] = {"x": Cv.build(z), "d": Cv.d, "a": dict(zip(ids, a)), "model": "M1", "M": Mx}
    return out


def pooled(D, pred, mask=None):
    w = D["w5"] > 0
    m = w if mask is None else (w & mask)
    r2 = (D["y5"] - pred) ** 2
    return float(np.sqrt(r2[m].sum() / m.sum()))


def ev_sums(D, pred, mask):
    w = D["w5"] > 0
    m = w & mask
    r2 = np.where(m, (D["y5"] - pred) ** 2, 0.0)
    ev = D["le_event"][D["bin_le"]]
    E = int(D["le_event"].max()) + 1
    return np.bincount(ev[m], weights=r2[m], minlength=E), np.bincount(ev[m], minlength=E).astype(float)


def rules_for(name, f, ctx):
    """P1-P5 of one model on TEST, exactly as S8 (the same draws)."""
    D1, le_tr, D1n, le_te, pr, prn, Dv1, lev, pv, cells, win, curves = (ctx[k] for k in ("D1", "le_tr", "D1n", "le_te", "pr", "prn", "Dv1", "lev", "pv",
                                                                                          "cells", "win", "curves"))
    PR = S8["pass_rules"]
    BR = S8["boot_reps"]
    out = {}
    e_ = cells[cells.eligible & (cells.fit == name)]
    n_el, n_ok = len(e_), int((e_.I_inside & e_.s_inside).sum())
    out["P1"] = (bool(n_el and n_ok >= need(n_el)), f"{n_ok} of {n_el} eligible cells")
    Lb = le_tr.L_km.to_numpy()[D1["bin_le"]]
    E = int(D1["le_event"].max()) + 1
    Mlt = ev_draws(E, BR, S8["seed"])
    parts, ok2 = [], True
    for sub, mask in (("all", np.ones(len(Lb), bool)), ("L < 1 km", Lb < 1)):
        sm = ev_sums(D1, pr[name], mask)
        bm = np.sqrt((Mlt @ sm[0]) / np.maximum(Mlt @ sm[1], 1))
        for b in ("B0", "B1", "B2"):
            sb = ev_sums(D1, pr[b], mask)
            bb = np.sqrt((Mlt @ sb[0]) / np.maximum(Mlt @ sb[1], 1))
            dpt = float(np.sqrt(sb[0].sum() / sb[1].sum()) - np.sqrt(sm[0].sum() / sm[1].sum()))
            lo, hi = ci(bb - bm)
            ok2 &= bool(dpt > 0 and lo > 0)
            parts.append(f"{sub} {b} {dpt:+.4f} [{lo:+.4f}, {hi:+.4f}]")
    out["P2"] = (ok2, "; ".join(parts))
    obs = np.array([c["obs"] for c in curves])
    mod = np.array([c[name] for c in curves])
    nmin = win.dry_minutes.to_numpy()
    draws = boot_items(pd.factorize(win.event_id)[0], BR, S8["seed"] + 1)
    t = 10
    ok = (nmin > t) & np.isfinite(obs[:, t])
    s_lo, s_hi = ok & (win.length.to_numpy() == LB[0][1]), ok & (win.length.to_numpy() == LB[2][1])
    err = obs[:, t] - mod[:, t]
    dpt = float(np.median(err[s_hi]) - np.median(err[s_lo]))
    lo, hi = ci([np.median(err[i[s_hi[i]]]) - np.median(err[i[s_lo[i]]]) for i in draws if s_hi[i].any() and s_lo[i].any()])
    out["P3"] = (bool(lo <= 0 <= hi), f"L > 3 minus L < 1 at 10 min {dpt:+.3f} [{lo:+.3f}, {hi:+.3f}] dB")
    Lv = lev.L_km.to_numpy()[Dv1["bin_le"]]
    rat = []
    for nm, mt, mv in [("all", np.ones(len(Lb), bool), np.ones(len(Lv), bool))] + [(nm, (Lb >= lo_) & (Lb < hi_), (Lv >= lo_) & (Lv < hi_)) for (lo_, hi_), nm in LB]:
        rat.append((nm, pooled(D1, pr[name], mt) / pooled(Dv1, pv[name], mv)))
    allr = rat[0][1]
    out["P4"] = (all(abs(r / allr - 1) <= PR["p4_tol"] for _, r in rat[1:]), "; ".join(f"{nm} {r:.3f} ({r / allr:.2f} x)" for nm, r in rat))
    s_tr, n_tr = ev_sums(D1, pr[name], np.ones(len(Lb), bool))
    s_te, n_te = ev_sums(D1n, prn[name], np.ones(len(D1n["bin_le"]), bool))
    r5 = float(np.sqrt(s_te.sum() / n_te.sum()) / np.sqrt(s_tr.sum() / n_tr.sum()))
    out["P5"] = (bool(r5 <= PR["p5_max_ratio"]), f"TEST links / TRAIN links {r5:.3f}")
    return out


# =================================================================================================== figures
def fig_f(Cc, zc, draws_z, table, bands):
    R = np.geomspace(0.3, 40, 200)
    fig, axes = plt.subplots(1, len(bands), figsize=(4.2 * len(bands), 3.8), sharey=True)
    band_draws = np.array([closed_curve(Cc, z, bands, R) for z in draws_z]) if len(draws_z) else None
    for i, (ax, b) in enumerate(zip(axes, bands)):
        t = table[table.band == b]
        ax.errorbar(t.R, t.table_F, yerr=[t.table_F - t.table_lo, t.table_hi - t.table_F], fmt="o", color=INK, capsize=3, label="table F (S6 knots, 90 %)")
        ax.plot(R, closed_curve(Cc, zc, bands, R)[i], "-", color=INK, lw=1.4, label="closed F = A (R / 5)^b")
        if band_draws is not None:
            lo, hi = np.percentile(band_draws[:, i, :], [5, 95], axis=0)
            ax.plot(R, lo, ":", color=INK2, lw=1)
            ax.plot(R, hi, ":", color=INK2, lw=1, label="closed F, 90 % of 40 cold draws")
        ax.set_xscale("log")
        ax.set_title(b)
        ax.set_xlabel("rain R / g_e (mm/h)")
        style_axes(ax)
    axes[0].set_ylabel("F (dB, a_link = 1)")
    axes[0].legend(frameon=False, fontsize=7)
    fig.suptitle("S8b (TRAIN): the table F of the S6 candidate against the closed form", fontsize=9)
    fig.tight_layout()
    return fig


def fig_rmse(rt):
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    for k, (m, mk) in enumerate((("table F", "o"), ("closed F", "s"))):
        g = rt[rt.model == m]
        x = np.arange(len(g)) + (k - 0.5) * 0.15
        ax.errorbar(x, g.rmse, yerr=[g.rmse - g.lo, g.hi - g.rmse], fmt=mk, mfc="white" if k else INK, color=INK, capsize=3, label=m)
        ax.set_xticks(np.arange(len(g)), g.subset)
    ax.set_ylabel("TEST RMSE (dB), 90 % event bootstrap")
    ax.set_title("S8b: post-freeze evaluation of a pre-freeze model (TEST, V1)", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    style_axes(ax)
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    clock = time.time()
    workers = max(min((os.cpu_count() or 2) - 1, 3), 1)
    fit_file = P.S6_CACHE / "plug_shared_FP_F_and_P.pkl"
    sha = hashlib.sha256(fit_file.read_bytes()).hexdigest()
    mtime = pd.Timestamp(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(fit_file))))   # local machine time, as the freeze
    freeze = json.loads(P.S8_FREEZE.read_text())
    if sha != S["fit_sha256"]:
        print(f"STOP: the closed-F fit changed (sha256 {sha} against {S['fit_sha256']}).")
        return
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    bands_cfg = sid.merged_bands(d, DATA)[0]
    BANDS = list(bands_cfg)
    Dt1, let, Cc, spec = closed_model(BANDS)
    mods, (Ck, zk) = s6_models(Dt1, let, BANDS)
    zc = np.asarray(mods["closed F"]["z"])
    f0, fc = mods["candidate"]["fz"], mods["closed F"]["fz"]

    # ------------------------------------------------------------------------------------------- Part 1: TRAIN
    Pg = sid.Problem2(Dt1, "M4b", M)
    sid.fix_params(Pg, FIX)
    freq = let.groupby("band_name").frequency_ghz.median()
    kn = np.asarray(M["knots_mmh"], float)
    itu_inc = np.array([np.diff(itu_p838(float(freq[b]), "v")[0] * kn ** itu_p838(float(freq[b]), "v")[1]) for b in BANDS])   # START only
    xg = plan_start(Pg, "g4+Phalf", itu_inc)
    nb = len(BANDS)
    F5g = xg[Pg.names == Pg.fname].reshape(nb, -1)[:, :K5].sum(axis=1)
    P5g = xg[Pg.names == Pg.pname].reshape(nb, -1)[:, :K5].sum(axis=1)
    b0 = float(spec["F"][1]["starts"][0][0])
    z_gen = np.r_[F5g, b0, P5g, np.zeros(Cc.nl)]
    TR = S["train"]
    base = {"bands": BANDS, "z0": z_gen.tolist()}
    jobs = [{**base, "name": f"boot_{i:02d}", "kind": "boot", "i": i, "seed": [TR["bootstrap"]["seed"], zlib.crc32(b"closed F")]} for i in range(TR["bootstrap"]["reps"])]
    jobs += [{**base, "name": f"fake_{r}", "kind": "fake", "x_truth": Cc.build(zc).tolist(), "seed": [TR["fake"]["seed"], r]} for r in range(TR["fake"]["replicates"])]
    res = {}
    run_pool(jobs, workers, clock, res)
    est = readout(Cc, zc, BANDS)
    bz = [np.asarray(res[f"boot_{i:02d}"]["z"]) for i in range(TR["bootstrap"]["reps"]) if res[f"boot_{i:02d}"]["converged"]]
    dr = pd.DataFrame([readout(Cc, z, BANDS) for z in bz])
    fk = pd.DataFrame([readout(Cc, np.asarray(res[f"fake_{r}"]["z"]), BANDS) for r in range(TR["fake"]["replicates"])])
    ib = Cc.blocks[0][1] + nb                                     # index of b in z
    rows = []
    for k, v in est.items():
        lo, hi = ci(dr[k])
        within = all((abs(f - v) <= 0.01) if abs(v) < 0.05 else (abs(f / v - 1) <= 0.20) for f in fk[k])
        if k.startswith("a_"):
            at_b = bool(lo <= np.exp(-3) + 1e-9 or hi >= np.exp(3) - 1e-9)
        elif k == "b":
            at_b = bool(lo <= Cc.lo[ib] * 1.01 or hi >= Cc.hi[ib] * 0.99)
        else:
            at_b = bool(lo <= 0.001)
        rows.append({"quantity": k, "value": v, "boot_lo": lo, "boot_hi": hi, "draws": len(dr), "fake": ", ".join(f"{f:.4g}" for f in fk[k]),
                     "fake_within_3_of_3": within, "interval_at_bound": at_b, "identified": bool(within and not at_b)})
    trt = pd.DataFrame(rows)
    # closed F against table F (S6 knots and their bootstrap intervals)
    s6b = pd.read_csv(P.S6_BOOT).set_index("quantity")
    Rd = np.geomspace(0.5, 40, 400)
    cf_d = closed_curve(Cc, zc, BANDS, Rd)
    cf_k = closed_curve(Cc, zc, BANDS, np.array(KN, float))
    cmp_rows = []
    for i, b in enumerate(BANDS):
        tv = np.array([s6b.value[f"F[{b}]({k:g})"] for k in KN])
        td = np.interp(Rd, np.r_[0.0, KN], np.r_[0.0, tv])
        j = int(np.argmax(np.abs(cf_d[i] - td)))
        for kk, k in enumerate(KN):
            lo, hi = s6b.boot_lo[f"F[{b}]({k:g})"], s6b.boot_hi[f"F[{b}]({k:g})"]
            cmp_rows.append({"band": b, "R": k, "table_F": tv[kk], "table_lo": lo, "table_hi": hi, "closed_F": float(cf_k[i, kk]),
                             "difference": float(cf_k[i, kk] - tv[kk]), "closed_inside_table_interval": bool(lo <= cf_k[i, kk] <= hi),
                             "band_max_abs_difference_0.5_40": float(abs(cf_d[i, j] - td[j])), "band_max_at_R": float(Rd[j])})
    cmpt = pd.DataFrame(cmp_rows)
    print(f"Part 1 done ({(time.time() - clock) / 60:.1f} min)", flush=True)

    # ------------------------------------------------------------------------------------------- Part 2: TEST (the S8 cache)
    le_tr, _ = pickle.loads((P.S8_CACHE / "le_test_train.pkl").read_bytes())
    le_te, _ = pickle.loads((P.S8_CACHE / "le_test_test.pkl").read_bytes())
    le_tr, le_te = le_tr.reset_index(drop=True), le_te.reset_index(drop=True)
    Dt, Dn = sid.load_data(P.S8_CACHE / "data_test_train.npz", BANDS), sid.load_data(P.S8_CACHE / "data_test_test.npz", BANDS)
    gt = pd.read_csv(P.S8_G_E)
    gm = {e: 1.0 for e in set(le_tr.event_id) | set(le_te.event_id)}
    gm.update({r.event_id: (r.g_path if r.path_used else 1.0) for r in gt.itertuples()})
    D1, D1n = rescaled(Dt, le_tr, gm), rescaled(Dn, le_te, gm)
    g4 = pd.read_csv(P.COND4C_G)
    gmap = {s: dict(zip(g4[g4.set == s].event_id, g4[g4.set == s].g_path)) for s in ("TRAIN", "VAL")}
    lev = pickle.loads((P.SID3_CACHE / "le_val_all.pkl").read_bytes())[0].reset_index(drop=True)
    _, Dv = load_sets(BANDS, let, lev, gmap)
    Dv1 = Dv["V1"]
    allm = {"table F": f0, "closed F": fc, **b_models(Dt1, let, BANDS)}
    pr = {m: predict(D1, le_tr, f, f["M"])[0] for m, f in allm.items()}
    prn = {m: predict(D1n, le_te, allm[m], allm[m]["M"])[0] for m in ("table F", "closed F")}
    pv = {m: predict(Dv1, lev, allm[m], allm[m]["M"])[0] for m in ("table F", "closed F")}
    cells = through(D1, le_tr, {"table F": f0, "closed F": fc}, "TEST", "V1")
    Rm, T5, col = radar_max(d, {"TEST": (D1, le_tr)}, blank_test=False)
    W1 = {"TEST": {"table F": w_minutes(D1, le_tr, f0), "closed F": w_minutes(D1, le_tr, fc)}}
    win, curves = dry_windows({"TEST": (D1, le_tr)}, {"TEST": gm}, Rm, T5, col, W1)
    ctx = {"D1": D1, "le_tr": le_tr, "D1n": D1n, "le_te": le_te, "pr": pr, "prn": prn, "Dv1": Dv1, "lev": lev, "pv": pv, "cells": cells, "win": win, "curves": curves}
    rul = {m: rules_for(m, allm[m], ctx) for m in ("table F", "closed F")}
    s8r = pd.read_csv(P.S8_RULES)
    s8v = {r.rule.split()[0]: r.verdict for r in s8r.itertuples()}
    side = pd.DataFrame([{"rule": k, "S8 verdict (table F)": s8v.get(k, ""), "table F (this run)": "pass" if rul["table F"][k][0] else "FAIL",
                          "table F numbers": rul["table F"][k][1], "closed F": "pass" if rul["closed F"][k][0] else "FAIL", "closed F numbers": rul["closed F"][k][1]}
                         for k in ("P1", "P2", "P3", "P4", "P5")])
    reproduces = bool((side["S8 verdict (table F)"] == side["table F (this run)"]).all())
    same_verdicts = bool((side["table F (this run)"] == side["closed F"]).all())
    # the paired difference closed minus table
    Lb = le_tr.L_km.to_numpy()[D1["bin_le"]]
    E = int(D1["le_event"].max()) + 1
    Mp = ev_draws(E, S["boot_reps"], S["seed"])
    pd_rows, rm_rows = [], []
    for sub, mask in (("all", np.ones(len(Lb), bool)), ("L < 1 km", Lb < 1)):
        bt = {}
        for m in ("table F", "closed F"):
            s_, n_ = ev_sums(D1, pr[m], mask)
            bt[m] = np.sqrt((Mp @ s_) / np.maximum(Mp @ n_, 1))
            lo, hi = ci(bt[m])
            rm_rows.append({"subset": sub, "model": m, "rmse": float(np.sqrt(s_.sum() / n_.sum())), "lo": lo, "hi": hi})
        pt = rm_rows[-1]["rmse"] - rm_rows[-2]["rmse"]
        lo, hi = ci(bt["closed F"] - bt["table F"])
        pd_rows.append({"subset": sub, "closed_minus_table_db": pt, "lo": lo, "hi": hi, "within_band": bool(-S["rule_band_db"] <= lo and hi <= S["rule_band_db"])})
    pdt, rmt = pd.DataFrame(pd_rows), pd.DataFrame(rm_rows)
    obs = np.array([c["obs"] for c in curves])
    nmin = win.dry_minutes.to_numpy()
    dry_rows = []
    for t in (5, 10, 20):
        sel = (nmin > t) & np.isfinite(obs[:, t])
        r = {"t_min": t, "windows": int(sel.sum())}
        for m in ("table F", "closed F"):
            r[f"err {m}"] = float(np.median(obs[sel, t] - np.array([c[m] for c in curves])[sel, t]))
        dry_rows.append(r)
    dryt = pd.DataFrame(dry_rows)
    adopted = bool(same_verdicts and pdt.within_band.all())
    print(f"Part 2 done ({(time.time() - clock) / 60:.1f} min): same verdicts {same_verdicts}, paired within band {bool(pdt.within_band.all())} -> adopted {adopted}",
          flush=True)

    # save, figures, report
    testt = pd.concat([side.assign(part="P1-P5"), pdt.assign(part="paired difference"), rmt.assign(part="RMSE"), dryt.assign(part="dry path")], ignore_index=True)
    for tab, path in ((trt.round(5), P.S8B_TRAIN), (cmpt.round(4), P.S8B_COMPARE_F), (testt.round(4), P.S8B_TEST)):
        save(tab, path, index=False)
    for fig, path in ((fig_f(Cc, zc, bz, cmpt, BANDS), P.S8B_FIG_F), (fig_rmse(rmt), P.S8B_FIG_RMSE)):
        save(fig, path, dpi=100)
        plt.close(fig)
    q = dict(locals())
    report = write_report(q)
    save(report, P.S8B_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {(time.time() - clock) / 60:.1f} min")


def write_report(q):
    L = []
    a = L.append
    trt, cmpt, est = q["trt"], q["cmpt"], q["est"]
    t = trt.set_index("quantity")
    bands = q["BANDS"]
    a("# wet_antenna step S8b: post-freeze evaluation of a pre-freeze model (the closed form of F)\n")
    a("Written by `research/wet_antenna/test_s8b.py` (libraries `final_model.py`, `validation.py`, `system_id.py`, `length_test.py`). The rule: "
      "`wet_antenna.test_s8b` in `research/config.yaml` and `docs/ROADMAP.md` (S8b), written on 2026-10-19 after TEST was read, BEFORE this evaluation.\n")
    a("## 1. Status and the file times\n")
    a(STATUS + "\n")
    a(f"- The closed-F fit `{rel(P.S6_CACHE / 'plug_shared_FP_F_and_P.pkl')}`: sha256 `{q['sha']}`, modified {q['mtime']:%Y-%m-%d %H:%M:%S} (local time {time.strftime('%z')}).")
    a(f"- The S8 freeze `{rel(P.S8_FREEZE)}`: written {q['freeze']['written']} (commit bb30b5c). The fit is **{'BEFORE' if q['mtime'] < pd.Timestamp(q['freeze']['written']) else 'NOT before'}** the freeze.")
    a("- Nothing is refitted on TEST; the S8 verdict (P1 fails, P2-P5 pass, with the table F) is unchanged.\n")
    a("## 2. The closed-F model on TRAIN: numbers and identifiability\n")
    a(f"Cold event bootstrap ({int(t.draws.iloc[0])} converged draws of 40) and the S6 identifiability check (3 fake replicates; `{rel(P.S8B_TRAIN)}`):\n")
    a(md_table(trt.round(4)) + "\n")
    nid = trt[~trt.identified].quantity.tolist()
    a(f"**READING.** b = {est['b']:.3f} [{t.boot_lo['b']:.3f}, {t.boot_hi['b']:.3f}]; A = " + ", ".join(f"{est[f'A[{b}]']:.2f} [{t.boot_lo[f'A[{b}]']:.2f}, {t.boot_hi[f'A[{b}]']:.2f}]" for b in bands)
      + " dB at 5 mm/h; c = " + ", ".join(f"{est[f'c[{b}]']:.3f} [{t.boot_lo[f'c[{b}]']:.3f}, {t.boot_hi[f'c[{b}]']:.3f}]" for b in bands)
      + f" dB/km per mm/h. Identified by the S6 rule: {int(trt.identified.sum())} of {len(trt)}" + (f"; NOT identified: {', '.join(nid)}." if nid else ".") + "\n")
    a("## 3. The closed F against the table F (TRAIN)\n")
    a(f"Per band and knot: the table F of the S6 candidate with its S6 bootstrap interval, the closed F, and whether it lies inside (`{rel(P.S8B_COMPARE_F)}`):\n")
    a(md_table(cmpt[["band", "R", "table_F", "table_lo", "table_hi", "closed_F", "difference", "closed_inside_table_interval"]].round(3)) + "\n")
    mx = cmpt.groupby("band")[["band_max_abs_difference_0.5_40", "band_max_at_R"]].first()
    a("Largest |closed F - table F| over 0.5-40 mm/h: " + "; ".join(f"{b} {r['band_max_abs_difference_0.5_40']:.2f} dB at {r['band_max_at_R']:.1f} mm/h" for b, r in mx.iterrows())
      + f". Inside the table's interval at {int(cmpt.closed_inside_table_interval.sum())} of {len(cmpt)} knots.\n")
    a("![F](../figures/test_s8b_F.png)\n")
    a("## 4. TEST: P1-P5 side by side, and the paired difference\n")
    a(f"The table F recomputed here reproduces the S8 verdicts: **{q['reproduces']}**.\n")
    a(md_table(q["side"]) + "\n")
    a("TEST RMSE (pooled 5 min bins, 90 % event bootstrap) and the paired difference closed minus table (rule: the 90 % interval within [-0.01, +0.01] dB):\n")
    a(md_table(q["rmt"].round(4)) + "\n")
    a(md_table(q["pdt"].round(4)) + "\n")
    a("Dry-path median errors (W_obs - W_model, dB):\n")
    a(md_table(q["dryt"].round(3)) + "\n")
    a("![rmse](../figures/test_s8b_rmse.png)\n")
    a("## 5. The decision by the rule\n")
    a(f"(1) The P1-P5 verdicts of the closed F are the same as the table F's: **{q['same_verdicts']}**. (2) The paired 90 % interval lies within [-0.01, +0.01] dB "
      f"for all link-events and L < 1 km: **{bool(q['pdt'].within_band.all())}** ("
      + "; ".join(f"{r.subset} {r.closed_minus_table_db:+.4f} [{r.lo:+.4f}, {r.hi:+.4f}]" for r in q["pdt"].itertuples()) + ").\n")
    a(("**The closed form is ADOPTED as the written form of F.**" if q["adopted"] else
       "**The table F stays the equation; the closed form is reported as an approximate description.**") + "\n")
    a("## 6. The equation of F to use in S9\n")
    if q["adopted"]:
        a("```\nF_band(R) = A_band (R / 5)^b   (closed form; post-freeze evaluation of a pre-freeze model, adopted by the S8b rule)\n"
          f"b = {est['b']:.3f} [{t.boot_lo['b']:.3f}, {t.boot_hi['b']:.3f}]\n"
          + "\n".join(f"A[{b}] = {est[f'A[{b}]']:.2f} [{t.boot_lo[f'A[{b}]']:.2f}, {t.boot_hi[f'A[{b}]']:.2f}] dB" for b in bands)
          + "\nwith P_band = c_band R: " + ", ".join(f"c[{b}] = {est[f'c[{b}]']:.3f}" for b in bands) + " dB/km per mm/h (the closed-F fit)\n```\n")
        a("**Range.** The closed form was fitted on 0-40 mm/h, and F above 10 mm/h rests on 0.4 % of the TRAIN minutes (S6). It has no saturation: at 40 "
          "mm/h it is " + ", ".join(f"{r.difference:+.1f} dB" for r in cmpt[cmpt.R == 40].itertuples())
          + " above the table F (inside the table's wide interval). It must not be extrapolated above 40 mm/h. b of this plug-back refit "
          f"({est['b']:.2f}) differs from the knot fit of S6 (0.47): the knot fit described the table F, this is the model refitted with the form.\n")
    else:
        a("The table F of the S6 candidate (`tighten_s6_candidate.json`): F_band(R) = A_band G(R), linear between the knots; the closed form A_band (R / 5)^b "
          f"(b = {est['b']:.3f}) is an approximate description (section 3).\n")
    return "\n".join(L)


if __name__ == "__main__":
    main()
