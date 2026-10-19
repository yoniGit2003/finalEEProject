"""wet_antenna step S8: freeze, then read TEST ONCE.

  Part 1  comparison models fitted on TRAIN V1 only: B0 no WAA, B1 a constant WAA while it rains, B2 WAA without memory
          (B3 skipped: no literature WAA model is implemented in the repo); reported on TRAIN and VAL (used, optimistic).
  Part 2  the freeze: products/test_s8_freeze.json with the sha256 of the candidate, of each comparison fit, of the config
          block wet_antenna.test_s8 and the git commit of the code. The first run writes it and STOPS; it is committed; the
          second run checks every hash (and that the freeze file is committed) before TEST is read.
  Part 3  TEST: the baselines of the TEST events (the frozen baseline code), the TEST link-events (the S2 rules, TRAIN links
          and TEST links), g_e from the TEST events' City gauges, every frozen model run forward; P1-P5 and the information.

The candidate is M0 = the S6 candidate (decision 2026-10-19). STANDING RULE: WAA is never estimated by subtracting a
power-law rain term; ITU-R P.838 only as a size reference; the rain correction comes from the City gauges only (g_e, for TEST
events from the TEST events' City gauges); nothing of the candidate is refitted on TEST; nothing is changed after TEST is read.
Every rule is in config `wet_antenna.test_s8` (written before TEST is read). Caches: research/.cache/wet_antenna/test_s8/.
Run from the repository root (twice: the freeze, then TEST):
    python -m research.wet_antenna.test_s8
"""
import copy
import hashlib
import json
import os
import pickle
import subprocess
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cmlpipe.api import CMLData
from research.baseline import baseline as BL
from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import length_test as LTL
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.event_factors import radar_gauge_primary, rescaled
from research.wet_antenna.final_model import CFG, DATA, J, M, SHAPES, WA, Curves, load_sets, predict, s6_models, through
from research.wet_antenna.validation import LB, boot_items, ci, dry_windows, error_budget, ev_draws, need, radar_max, w_minutes
from research.wet_antenna.wind_drying import ProblemW

S = WA["test_s8"]
CACHE = P.S8_CACHE
VL = "used (optimistic)"
BS = ("B0", "B1", "B2")
KB = tuple(WA["return_s6b"]["refinements"]["k_bounds"])
DRY = P.REPO / "research" / "explorations" / "s8_checks" / "dryrun"      # the VAL dry run of the evaluation code
STANDING = ("**Standing rule.** WAA is never estimated by subtracting a power-law rain term; ITU-R P.838 only as a size reference; the rain "
            "correction comes from the City gauges only (g_e, for TEST events from the TEST events' City gauges); nothing of the candidate is "
            "refitted on TEST; nothing is changed after TEST is read.")


# =================================================================================================== the comparison models
def step_g(R):
    """B1: 0 below 0.1 mm/h, 1 from 0.1001 mm/h on (a near-step on the knots)."""
    return np.clip((np.asarray(R, float) - 0.1) / 1e-4, 0.0, 1.0)


def comp_setup(name, nb):
    """(spec, M of the model, the static problem name) of a comparison model."""
    lin = {**SHAPES["linear"]}
    Mx = copy.deepcopy(M)
    if name == "B1":
        Mx["knots_mmh"] = [0.0, 0.1, 0.1001] + [k for k in M["knots_mmh"] if k > 0.1001]
        F = ("form", {"g": step_g, "lo": [], "hi": [], "starts": [[]], "text": "w0 (R > 0.1 mm/h)", "names": []}, True)
    elif name == "B0":
        F = ("fixed", np.zeros((nb, len(Mx["knots_mmh"]))))
    else:
        F = "shared"
    return {"F": F, "P": ("form", lin, True)}, Mx


_W = {}


def job(jb):
    file = CACHE / f"{jb['name']}.pkl"
    if file.exists():
        return jb["name"], pickle.loads(file.read_bytes())
    t0 = time.time()
    if "Dt1" not in _W:
        _W["Dt1"] = sid.load_data(P.S5A_CACHE / "data_train_v1.npz", jb["bands"])
    D = _W["Dt1"]
    spec, Mx = comp_setup(jb["model"], len(jb["bands"]))
    Pm = sid.Problem2(D, "M1", Mx)
    Cv = Curves(D, spec, np.zeros(len(Pm.names)), S["comparison"]["lag_min"], Mx=Mx, Pm=Pm)
    z = np.asarray(jb["z0"])
    nit = 0
    for cap in jb["caps"]:
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


def comp_models(Dt1, let, BANDS, mods, Ck, zk, clock, workers):
    """Fit B0-B2 on TRAIN (cached) and return {name: fz}."""
    ids = LTL.link_ids(Dt1, let)
    Cm0, zc = mods["candidate"]["C"], mods["candidate"]["z"]
    _, Pc, _ = Cm0.values(zc)
    c5 = Pc[:, -1] / Cm0.kn[-1] * 5
    la0 = zc[Cm0.nc:Cm0.nc + Cm0.nl]
    nF = Ck.blocks[0][2]
    z0 = {"B0": np.r_[c5, np.zeros(len(la0))], "B1": np.r_[np.ones(len(BANDS)), c5, np.zeros(len(la0))], "B2": np.r_[zk[:nF], c5, la0]}
    caps = [S["comparison"]["caps"]["start"], S["comparison"]["caps"]["polish"]]
    res = {}
    run_pool([{"name": f"fit_{b}", "model": b, "bands": BANDS, "z0": z0[b].tolist(), "caps": caps} for b in BS], workers, clock, res)
    fz = {}
    for b in BS:
        spec, Mx = comp_setup(b, len(BANDS))
        Pm = sid.Problem2(Dt1, "M1", Mx)
        Cv = Curves(Dt1, spec, np.zeros(len(Pm.names)), S["comparison"]["lag_min"], Mx=Mx, Pm=Pm)
        z = np.asarray(res[f"fit_{b}"]["z"])
        F, Pv, a = Cv.values(z)
        fz[b] = {"x": Cv.build(z), "d": Cv.d, "a": dict(zip(ids, a)), "model": "M1", "M": Mx, "loss": res[f"fit_{b}"]["loss"],
                 "converged": res[f"fit_{b}"]["converged"], "c": Pv[:, -1] / Cv.kn[-1], "F5": F[:, int(np.argmin(np.abs(Cv.kn - 5)))]}
    return fz


# =================================================================================================== RMSE helpers
def bin_err(D, pred):
    w = D["w5"] > 0
    return np.where(w, (D["y5"] - pred) ** 2, 0.0), w


def pooled(D, pred, mask=None):
    r2, w = bin_err(D, pred)
    m = w if mask is None else (w & mask)
    return float(np.sqrt(r2[m].sum() / m.sum())) if m.any() else np.nan


def ev_sums(D, pred, mask):
    r2, w = bin_err(D, pred)
    ev = D["le_event"][D["bin_le"]]
    E = int(D["le_event"].max()) + 1
    m = w & mask
    return np.bincount(ev[m], weights=r2[m], minlength=E), np.bincount(ev[m], minlength=E).astype(float)


def file_sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def config_sha():
    block = {k: v for k, v in S.items() if k != "results"}
    return hashlib.sha256(yaml.safe_dump(block, sort_keys=True, allow_unicode=True).encode("utf-8")).hexdigest()


def git(*args):
    return subprocess.run(["git", *args], cwd=P.REPO, capture_output=True, text=True).stdout.strip()


# =================================================================================================== TEST data (Part 3)
def test_baselines(d):
    f = P.S8_BASELINES
    if f.exists():
        return pd.read_csv(f)
    s, tl, state, dew, amp = BL.table_inputs(d, ("train", "test"))
    unc = pd.read_csv(P.BASELINE_UNCERTAINTY)
    t = BL.event_table(d, s, tl, state, dew, amp, unc, CFG["baseline"], ("test",))
    t["use_in_waa_main"] = ~(t.missing | t.excursion | t.daily_cycle_large | t.degraded | t.toggling)
    save(t.round(4), f, index=False)
    return pd.read_csv(f)


def test_set(d, bands_cfg, link_split, tb):
    f, fl = CACHE / f"data_test_{link_split}.npz", CACHE / f"le_test_{link_split}.pkl"
    if f.exists() and fl.exists():
        le, drops = pickle.loads(fl.read_bytes())
        return sid.load_data(f, list(bands_cfg)), le.reset_index(drop=True), drops
    D, le, drops = sid.build_dataset(d, {**DATA, "bands_ghz": bands_cfg, "max_length_km": WA["system_id_s3"]["data"]["max_length_km_all"]},
                                     event_split="test", link_split=link_split, baselines=tb)
    sid.save_data(D, f)
    fl.write_bytes(pickle.dumps((le, drops)))
    return sid.load_data(f, list(bands_cfg)), le.reset_index(drop=True), drops


def g_test(D, le):
    rg = radar_gauge_primary(D, le)
    mn = WA["conditions_s4c"]["g_e"]["min_gauge_mm"]
    rg["path_used"] = rg.gauge_mm_per_link >= mn
    rg["g_path"] = np.where(rg.path_used, rg.radar_over_gauge, 1.0)
    u = rg[rg.path_used]
    G = float(u.radar_mm_sum.sum() / u.gauge_mm_sum.sum()) if len(u) else 1.0
    rg["G_set"] = G
    return rg, dict(zip(rg.event_id, rg.g_path)), G


def m1_fz(Dt1, let, Ck, zk):
    """M1 of S6b (box B): the S6b cache fit, tau_f_dry free, k = 0 (no covariate)."""
    spec = {"F": "shared", "P": ("form", {**SHAPES["linear"]}, True)}
    x_ref = np.r_[Ck.build(zk), 0.0]
    Cv = Curves(Dt1, spec, x_ref, Ck.d, Pm=ProblemW(Dt1, M, np.zeros(len(Dt1["rain"])), KB), free=["tau_f_dry"])
    z = np.asarray(pickle.loads((P.S6B_CACHE / "fit_M1.pkl").read_bytes())["z"])
    ids = LTL.link_ids(Dt1, let)
    return {"x": Cv.build(z), "d": Cv.d, "a": dict(zip(ids, np.exp(z[Cv.nc:Cv.nc + Cv.nl]))), "model": "M4b", "M": M}


def m1_minutes(D, le, fz):
    D0 = dict(D)
    D0["le_c"] = np.zeros(len(le))
    Pm = ProblemW(D0, fz["M"], np.zeros(len(D["rain"])), KB, agg=1)
    x, _ = LTL.x_for(Pm, fz, LTL.link_ids(D, le))
    x[Pm.names == Pm.pname] = 0.0
    return Pm.predict(x, fz["d"])


def funnel(d, drops_tr, drops_te, let_, lte_, Dv1):
    rows = []
    for lab, drops in (("TEST events, TRAIN links", drops_tr), ("TEST events, TEST links", drops_te)):
        items = list(drops.items())
        rows.append({"stage": f"{lab}: (event, sub-link) pairs with a baseline row", "count": items[0][1]})
        for k, v in items[1:-1]:
            rows.append({"stage": f"{lab}:   minus {k}", "count": v})
        rows.append({"stage": f"{lab}: link-events kept", "count": items[-1][1]})
    ev = d.events()
    rows.insert(0, {"stage": "TEST events", "count": int((ev.split_event == "test").sum())})
    ft = pd.DataFrame(rows)
    lc = []
    for (lo, hi), nm in LB:
        r = {"length": nm}
        for lab, le in (("TRAIN links", let_), ("TEST links", lte_)):
            m = (le.L_km >= lo) & (le.L_km < hi)
            r[f"{lab} link-events"], r[f"{lab} links"] = int(m.sum()), int(le[m].cml_id.nunique())
        lc.append(r)
    lct = pd.DataFrame(lc)
    lct.loc[len(lct)] = {"length": "all", **{c: lct[c].sum() for c in lct.columns if c != "length"}}
    return ft, lct


# =================================================================================================== figures
def fig_rmse(rm):
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    for ax, sub in zip(axes, ("all", "L < 1 km")):
        g = rm[(rm.subset == sub) & rm.model.isin(["candidate", *BS])]
        x = np.arange(len(g))
        ax.errorbar(x, g.rmse, yerr=[g.rmse - g.lo, g.hi - g.rmse], fmt="o", color=INK, capsize=3)
        ax.set_xticks(x, g.model)
        ax.set_title(f"TEST, {sub} (V1)")
        style_axes(ax)
    axes[0].set_ylabel("RMSE of 5 min A (dB), 90 % event bootstrap")
    fig.suptitle("P2: the candidate against no WAA (B0), constant WAA (B1), WAA without memory (B2)", fontsize=9)
    fig.tight_layout()
    return fig


def fig_cells(cells):
    e = cells[cells.eligible]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for ax, (q, lo, hi, mod, lab) in zip(axes, (("I_obs", "I_lo", "I_hi", "I_model", "intercept I (dB)"), ("s_obs", "s_lo", "s_hi", "s_model", "slope s (dB/km)"))):
        x = np.arange(len(e))
        ax.errorbar(x, e[q], yerr=[e[q] - e[lo], e[hi] - e[q]], fmt="o", color=INK, capsize=3, label="TEST observed, 90 %")
        ax.plot(x, e[mod], "s", mfc="white", color=INK, label="candidate through the test")
        ax.set_xticks(x, [f"{b}\n{c}" for b, c in zip(e.band, e["class"])], fontsize=7)
        ax.set_ylabel(lab)
        style_axes(ax)
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("P1: the TEST length test cells (eligible: >= 8 links under 1 km)", fontsize=9)
    fig.tight_layout()
    return fig


def fig_dry(win, curves):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), sharey=True)
    obs = np.array([c["obs"] for c in curves]) if len(curves) else np.zeros((0, 60))
    nmin = win.dry_minutes.to_numpy() if len(win) else np.zeros(0)
    for ax, (_, nm) in zip(axes, LB):
        sel = (win.length.to_numpy() == nm) if len(win) else np.zeros(0, bool)
        for lab, ls in (("measured", "-"), ("M0", "--"), ("M1", ":")):
            arr = obs if lab == "measured" else np.array([c[lab] for c in curves])
            med = [np.nanmedian(arr[sel & (nmin > t), t]) if (sel & (nmin > t)).sum() >= 3 else np.nan for t in range(60)]
            ax.plot(range(60), med, ls, color=INK, lw=1.4, label=lab if lab == "measured" else f"model {lab}{' (candidate)' if lab == 'M0' else ' (box B)'}")
        ax.set_title(f"{nm}: {int(sel.sum())} windows")
        ax.set_xlabel("minutes after the path became dry")
        style_axes(ax)
    axes[0].set_ylabel("median W (dB)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("P3: TEST dry-path W, measured against the candidate (M0) and M1 (box B)", fontsize=9)
    fig.tight_layout()
    return fig


def fig_ratio(rt):
    fig, ax = plt.subplots(figsize=(6.5, 3.6))
    g = rt[rt.subset != "all"]
    allr = float(rt[rt.subset == "all"].ratio.iloc[0])
    ax.plot(range(len(g)), g.ratio, "o", color=INK)
    ax.axhline(allr, color=INK2, ls="--", lw=1)
    ax.axhspan(allr * (1 - S["pass_rules"]["p4_tol"]), allr * (1 + S["pass_rules"]["p4_tol"]), facecolor="none", edgecolor=INK2, hatch="//", lw=0)
    ax.text(len(g) - 0.6, allr, " all classes", va="bottom", fontsize=8)
    ax.set_xticks(range(len(g)), g.subset)
    ax.set_ylabel("TEST / VAL RMSE")
    ax.set_title("P4: TEST / VAL per length class (hatched: +-15 % of the all-class ratio)", fontsize=9)
    style_axes(ax)
    fig.tight_layout()
    return fig


def fig_new(nl):
    fig, ax = plt.subplots(figsize=(6, 3.6))
    x = np.arange(len(nl))
    ax.errorbar(x, nl.rmse, yerr=[nl.rmse - nl.lo, nl.hi - nl.rmse], fmt="o", color=INK, capsize=3)
    ax.set_xticks(x, nl.links)
    ax.set_ylabel("TEST RMSE (dB)")
    ax.set_title("P5: TEST links (a_link = 1, never fitted) against TRAIN links, TEST events", fontsize=9)
    style_axes(ax)
    fig.tight_layout()
    return fig


def fig_budget(bud):
    b = bud[bud.part != "all"]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    hatches = ["", "//", "..", "xx", "\\\\"]
    for i, r in enumerate(b.itertuples()):
        ax.barh(i, r.share_of_squared_error, color="white", edgecolor=INK, hatch=hatches[i % 5])
        ax.text(r.share_of_squared_error + 0.01, i, f"{r.share_of_squared_error:.0%} (RMSE {r.rmse_db:.2f} dB)", va="center", fontsize=8)
    ax.set_yticks(range(len(b)), b.part, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.15)
    ax.set_xlabel("share of the squared error on TEST (candidate, V1)")
    ax.set_title("The TEST error budget", fontsize=9)
    style_axes(ax)
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    clock = time.time()
    workers = max(min((os.cpu_count() or 2) - 1, 3), 1)
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    bands_cfg = sid.merged_bands(d, DATA)[0]
    BANDS = list(bands_cfg)
    let = pickle.loads((P.SID3_CACHE / "le_train_all.pkl").read_bytes())[0].reset_index(drop=True)
    lev = pickle.loads((P.SID3_CACHE / "le_val_all.pkl").read_bytes())[0].reset_index(drop=True)
    g4 = pd.read_csv(P.COND4C_G)
    gmap = {s: dict(zip(g4[g4.set == s].event_id, g4[g4.set == s].g_path)) for s in ("TRAIN", "VAL")}
    Dt1 = sid.load_data(P.S5A_CACHE / "data_train_v1.npz", BANDS)
    _, Dv = load_sets(BANDS, let, lev, gmap)
    Dv1 = Dv["V1"]
    mods, (Ck, zk) = s6_models(Dt1, let, BANDS)
    f0 = mods["candidate"]["fz"]

    # ------------------------------------------------------------------------------------------- Part 1: the comparison models
    fzB = comp_models(Dt1, let, BANDS, mods, Ck, zk, clock, workers)
    allm = {"candidate": f0, **fzB}
    cmp_rows = []
    for s_, (D, le) in (("TRAIN", (Dt1, let)), ("VAL", (Dv1, lev))):
        for m, f in allm.items():
            p = predict(D, le, f, f["M"])[0]
            r = {"set": s_ if s_ == "TRAIN" else f"VAL ({VL})", "model": m, "rmse_pooled": pooled(D, p),
                 "rmse_pooled_L<1km": pooled(D, p, (le.L_km.to_numpy()[D["bin_le"]] < 1)), "rmse_judge": LTL.val_metrics(D, p, J)["rmse"]}
            if m in fzB:
                r.update({f"c[{b}]": float(fzB[m]["c"][i]) for i, b in enumerate(BANDS)})
                r.update({f"W5[{b}]": float(fzB[m]["F5"][i]) for i, b in enumerate(BANDS)})
                r["train_loss"], r["converged"] = fzB[m]["loss"], fzB[m]["converged"]
            cmp_rows.append(r)
    cmp_t = pd.DataFrame(cmp_rows)
    save(cmp_t.round(4), P.S8_COMPARE, index=False)
    print(f"Part 1 done ({(time.time() - clock) / 60:.1f} min)\n{cmp_t.round(4).to_string()}", flush=True)

    if os.environ.get("S8_DRYRUN") == "VAL":                     # check the evaluation code on VAL (no TEST read; outputs to explorations)
        DRY.mkdir(parents=True, exist_ok=True)
        half = (lev.cml_id % 2 == 0).to_numpy()
        D1, D1n = sid.subset(Dv1, half), sid.subset(Dv1, ~half)
        le_a, le_b = lev[half].reset_index(drop=True), lev[~half].reset_index(drop=True)
        D0a = sid.subset(Dv["V0"], half)
        rg, g_v, G = g_test(D0a, le_a)
        drops = pickle.loads((P.SID3_CACHE / "le_val_all.pkl").read_bytes())[1]
        evaluate({"d": d, "D1": D1, "D1n": D1n, "le_tr": le_a, "le_te": le_b, "T": {"V0": (D0a, sid.subset(Dv["V0"], ~half)), "V1": (D1, D1n), "V3": (D1, D1n)},
                  "rg": rg, "g_te": g_v, "G": G, "drops_tr": drops, "drops_te": drops, "allm": allm, "f0": f0, "fm1": m1_fz(Dt1, let, Ck, zk), "Dv1": Dv1,
                  "lev": lev, "cmp_t": cmp_t, "freeze": {"written": "dry run", "candidate": {"file": "-", "sha256": "-"}, "comparison_fits": {},
                  "m1_box_b_fit": {"file": "-", "sha256": "-"}, "config_block_test_s8": "-", "code_commit": "-"}, "code_diff": [], "clock": clock, "dry": True})
        return
    # ------------------------------------------------------------------------------------------- Part 2: the freeze
    hashes = {"candidate": {"file": rel(P.S6_CANDIDATE), "sha256": file_sha(P.S6_CANDIDATE)},
              "comparison_fits": {b: {"file": rel(CACHE / f"fit_{b}.pkl"), "sha256": file_sha(CACHE / f"fit_{b}.pkl")} for b in BS},
              "m1_box_b_fit": {"file": rel(P.S6B_CACHE / "fit_M1.pkl"), "sha256": file_sha(P.S6B_CACHE / "fit_M1.pkl")},
              "config_block_test_s8": config_sha()}
    if not P.S8_FREEZE.exists():
        freeze = {"step": "S8", "written": time.strftime("%Y-%m-%d %H:%M"), **hashes, "code_commit": git("rev-parse", "HEAD"),
                  "code_dirty_files_at_freeze": [l for l in git("status", "--porcelain", "--", "research/wet_antenna/*.py", "research/baseline/*.py",
                                                                    "research/config.yaml").splitlines()],
                  "rule": "TEST is read only after this file is committed and every hash still matches"}
        save(json.dumps(freeze, indent=1), P.S8_FREEZE)
        print(f"\nFREEZE WRITTEN: {rel(P.S8_FREEZE)}. Commit it, then run this step again to read TEST.")
        return
    freeze = json.loads(P.S8_FREEZE.read_text())
    bad = [k for k in ("candidate", "comparison_fits", "m1_box_b_fit", "config_block_test_s8") if freeze[k] != hashes[k]]
    committed = git("ls-files", "--error-unmatch", str(P.S8_FREEZE.relative_to(P.REPO)).replace("\\", "/")) != "" and \
        git("status", "--porcelain", "--", str(P.S8_FREEZE.relative_to(P.REPO)).replace("\\", "/")) == ""
    if bad or not committed:
        print(f"STOP: the freeze does not hold (hashes differ: {bad}; freeze file committed and unchanged: {committed}). TEST NOT READ.")
        return
    code_diff = git("diff", "--name-only", freeze["code_commit"], "--", "research/wet_antenna/*.py", "research/baseline/*.py").splitlines()
    print(f"freeze verified (commit {freeze['code_commit'][:10]}); code changed since the freeze commit: {code_diff or 'none'}", flush=True)

    # ------------------------------------------------------------------------------------------- Part 3: TEST, read once
    tb = test_baselines(d)
    Dt_tr, le_tr, drops_tr = test_set(d, bands_cfg, "train", tb)
    Dt_te, le_te, drops_te = test_set(d, bands_cfg, "test", tb)
    rg, g_te, G = g_test(Dt_tr, le_tr)
    T = {"V0": (Dt_tr, Dt_te), "V1": (rescaled(Dt_tr, le_tr, g_te), rescaled(Dt_te, le_te, g_te)),
         "V3": (rescaled(Dt_tr, le_tr, {e: G for e in le_tr.event_id.unique()}), rescaled(Dt_te, le_te, {e: G for e in le_te.event_id.unique()}))}
    D1, D1n = T["V1"]
    print(f"TEST: {len(le_tr)} link-events on TRAIN links, {len(le_te)} on TEST links; g_e used for {int(rg.path_used.sum())} of {len(rg)} events "
          f"({(time.time() - clock) / 60:.1f} min)", flush=True)
    fm1 = m1_fz(Dt1, let, Ck, zk)
    evaluate({"d": d, "D1": D1, "D1n": D1n, "le_tr": le_tr, "le_te": le_te, "T": T, "rg": rg, "g_te": g_te, "G": G, "drops_tr": drops_tr,
              "drops_te": drops_te, "allm": allm, "f0": f0, "fm1": fm1, "Dv1": Dv1, "lev": lev, "cmp_t": cmp_t, "freeze": freeze,
              "code_diff": code_diff, "clock": clock, "dry": False})


def evaluate(c):
    """Part 3 on the data in c (TEST; or VAL in the dry run): P1-P5, the information items, tables, figures and the report."""
    d, D1, D1n, le_tr, le_te, T, rg, g_te, G = (c[k] for k in ("d", "D1", "D1n", "le_tr", "le_te", "T", "rg", "g_te", "G"))
    drops_tr, drops_te, allm, f0, fm1, Dv1, lev = (c[k] for k in ("drops_tr", "drops_te", "allm", "f0", "fm1", "Dv1", "lev"))
    cmp_t, freeze, code_diff, clock, dry = (c[k] for k in ("cmp_t", "freeze", "code_diff", "clock", "dry"))
    pr = {m: predict(D1, le_tr, f, f["M"])[0] for m, f in allm.items()}
    Lb = le_tr.L_km.to_numpy()[D1["bin_le"]]
    BR = S["boot_reps"]

    # P1 the length test
    cells = through(D1, le_tr, {"candidate": f0}, "TEST", "V1")
    e_ = cells[cells.eligible]
    n_el, n_ok = len(e_), int((e_.I_inside & e_.s_inside).sum())
    p1 = bool(n_el and n_ok >= need(n_el))

    # P2 the candidate against B0-B2, all and L < 1 km
    E = int(D1["le_event"].max()) + 1
    Mlt = ev_draws(E, BR, S["seed"])
    rm_rows, p2_parts = [], []
    for sub, mask in (("all", np.ones(len(Lb), bool)), ("L < 1 km", Lb < 1)):
        sums = {m: ev_sums(D1, pr[m], mask) for m in pr}
        bc = {m: np.sqrt((Mlt @ s_) / np.maximum(Mlt @ n_, 1)) for m, (s_, n_) in sums.items()}
        for m in pr:
            lo, hi = ci(bc[m])
            rm_rows.append({"set": "TEST", "subset": sub, "model": m, "rmse": float(np.sqrt(sums[m][0].sum() / sums[m][1].sum())), "lo": lo, "hi": hi})
        for b in BS:
            dpt = float(np.sqrt(sums[b][0].sum() / sums[b][1].sum()) - np.sqrt(sums["candidate"][0].sum() / sums["candidate"][1].sum()))
            lo, hi = ci(bc[b] - bc["candidate"])
            p2_parts.append({"subset": sub, "model": b, "diff_db": dpt, "lo": lo, "hi": hi, "ok": bool(dpt > 0 and lo > 0)})
    p2t = pd.DataFrame(p2_parts)
    p2 = bool(p2t.ok.all())

    # P3 the dry path on TEST (M0; M1 for information)
    Rm, T5, col = radar_max(d, {"TEST": (D1, le_tr)}, blank_test=False)
    W1 = {"TEST": {"M0": w_minutes(D1, le_tr, f0), "M1": m1_minutes(D1, le_tr, fm1)}}
    win, curves = dry_windows({"TEST": (D1, le_tr)}, {"TEST": g_te}, Rm, T5, col, W1)
    dry_rows, p3, p3_txt = [], None, "no windows"
    if len(win):
        obs = np.array([c["obs"] for c in curves])
        mods_ = {m: np.array([c[m] for c in curves]) for m in ("M0", "M1")}
        nmin = win.dry_minutes.to_numpy()
        draws = boot_items(pd.factorize(win.event_id)[0], BR, S["seed"] + 1)
        for t in (5, 10, 20):
            for _, nm in LB + [((0, 100), "all")]:
                sel = (nmin > t) & np.isfinite(obs[:, t]) & ((win.length.to_numpy() == nm) if nm != "all" else True)
                r = {"t_min": t, "length": nm, "windows": int(sel.sum()), "events": int(win[sel].event_id.nunique())}
                for m in ("M0", "M1"):
                    err = obs[:, t] - mods_[m][:, t]
                    r[f"err_{m}"] = float(np.median(err[sel])) if sel.any() else np.nan
                    r[f"err_{m}_lo"], r[f"err_{m}_hi"] = ci([np.median(err[i[sel[i]]]) for i in draws if sel[i].any()])
                dry_rows.append(r)
        t = 10
        ok = (nmin > t) & np.isfinite(obs[:, t])
        s_lo, s_hi = ok & (win.length.to_numpy() == LB[0][1]), ok & (win.length.to_numpy() == LB[2][1])
        if s_lo.sum() >= 10 and s_hi.sum() >= 10:
            err = obs[:, t] - mods_["M0"][:, t]
            dpt = float(np.median(err[s_hi]) - np.median(err[s_lo]))
            lo, hi = ci([np.median(err[i[s_hi[i]]]) - np.median(err[i[s_lo[i]]]) for i in draws if s_hi[i].any() and s_lo[i].any()])
            p3 = bool(lo <= 0 <= hi)
            p3_txt = f"L > 3 minus L < 1 at 10 min {dpt:+.3f} [{lo:+.3f}, {hi:+.3f}] dB ({int(s_hi.sum())} / {int(s_lo.sum())} windows)"
        else:
            p3_txt = f"not testable ({int(s_hi.sum())} windows L > 3 km, {int(s_lo.sum())} L < 1 km; rule >= 10 each)"
    dryt = pd.DataFrame(dry_rows)

    # P4 TEST / VAL per length class
    pv = predict(Dv1, lev, f0, f0["M"])[0]
    Lv = lev.L_km.to_numpy()[Dv1["bin_le"]]
    rt = []
    for nm, mt, mv in [("all", np.ones(len(Lb), bool), np.ones(len(Lv), bool))] + [(nm, (Lb >= lo) & (Lb < hi), (Lv >= lo) & (Lv < hi)) for (lo, hi), nm in LB]:
        a_, b_ = pooled(D1, pr["candidate"], mt), pooled(Dv1, pv, mv)
        rt.append({"subset": nm, "rmse_TEST": a_, "rmse_VAL": b_, "ratio": a_ / b_})
    rt = pd.DataFrame(rt)
    allr = float(rt.ratio.iloc[0])
    rt["rel_to_all"] = rt.ratio / allr
    rt["ok"] = (rt.rel_to_all - 1).abs() <= S["pass_rules"]["p4_tol"]
    p4 = bool(rt[rt.subset != "all"].ok.all())

    # P5 new links
    pn = predict(D1n, le_te, f0, f0["M"])[0]
    En = int(D1n["le_event"].max()) + 1
    s_tr, n_tr = ev_sums(D1, pr["candidate"], np.ones(len(Lb), bool))
    s_te, n_te = ev_sums(D1n, pn, np.ones(len(D1n["bin_le"]), bool))
    ev_tr = pd.factorize(le_tr.event_id)[1]
    ev_te = pd.factorize(le_te.event_id)[1]
    all_ev = sorted(set(ev_tr) | set(ev_te))
    a_tr = pd.Series(0.0, index=all_ev); a_tr[list(ev_tr)] += s_tr
    b_tr = pd.Series(0.0, index=all_ev); b_tr[list(ev_tr)] += n_tr
    a_te = pd.Series(0.0, index=all_ev); a_te[list(ev_te)] += s_te
    b_te = pd.Series(0.0, index=all_ev); b_te[list(ev_te)] += n_te
    Mn = ev_draws(len(all_ev), BR, S["seed"] + 2)
    r_tr, r_te = float(np.sqrt(a_tr.sum() / b_tr.sum())), float(np.sqrt(a_te.sum() / b_te.sum()))
    bt_tr = np.sqrt((Mn @ a_tr.to_numpy()) / np.maximum(Mn @ b_tr.to_numpy(), 1))
    bt_te = np.sqrt((Mn @ a_te.to_numpy()) / np.maximum(Mn @ b_te.to_numpy(), 1))
    ratio5 = r_te / r_tr
    lo5, hi5 = ci(bt_te / bt_tr)
    p5 = bool(ratio5 <= S["pass_rules"]["p5_max_ratio"])
    nl = pd.DataFrame([{"links": "TRAIN links", "n_links": le_tr.cml_id.nunique(), "rmse": r_tr, "lo": ci(bt_tr)[0], "hi": ci(bt_tr)[1]},
                       {"links": "TEST links (a_link = 1)", "n_links": le_te.cml_id.nunique(), "rmse": r_te, "lo": ci(bt_te)[0], "hi": ci(bt_te)[1]}])
    for (lo, hi), nm in LB:                                       # information: by length class
        mt = (Lb >= lo) & (Lb < hi)
        mn = (le_te.L_km.to_numpy()[D1n["bin_le"]] >= lo) & (le_te.L_km.to_numpy()[D1n["bin_le"]] < hi)
        nl.loc[len(nl)] = {"links": f"{nm}: TRAIN links / TEST links", "n_links": np.nan, "rmse": pooled(D1n, pn, mn) / pooled(D1, pr["candidate"], mt), "lo": np.nan, "hi": np.nan}
    print(f"P1-P5 done ({(time.time() - clock) / 60:.1f} min)", flush=True)

    rules = pd.DataFrame([
        {"rule": "P1 length test through the model", "numbers": f"{n_ok} of {n_el} eligible cells (needed {need(n_el)})", "verdict": "pass" if p1 else "FAIL",
         "claim": "W independent of L"},
        {"rule": "P2 the candidate beats B0, B1, B2", "numbers": "; ".join(f"{r.subset} {r.model} {r.diff_db:+.4f} [{r.lo:+.4f}, {r.hi:+.4f}] dB" for r in p2t.itertuples()),
         "verdict": "pass" if p2 else "FAIL", "claim": "the value of the memory and of the antenna term"},
        {"rule": "P3 dry-path error equal across lengths", "numbers": p3_txt, "verdict": "not testable" if p3 is None else ("pass" if p3 else "FAIL"),
         "claim": "W independent of L"},
        {"rule": "P4 TEST / VAL per length class within +-15 % of all", "numbers": "; ".join(f"{r.subset} {r.ratio:.3f} ({r.rel_to_all:.2f} x)" for r in rt.itertuples()),
         "verdict": "pass" if p4 else "FAIL", "claim": "the length behaviour"},
        {"rule": "P5 TEST links at most 20 % above TRAIN links", "numbers": f"{r_te:.3f} / {r_tr:.3f} = {ratio5:.3f} [{lo5:.3f}, {hi5:.3f}]", "verdict": "pass" if p5 else "FAIL",
         "claim": "the link multiplier"}])
    holds = bool(p1 and p2 and (p3 is not False) and p4 and p5)

    # information items
    info = []
    le_b = D1["bin_le"]
    w5 = D1["w5"] > 0
    r2c = np.where(w5, (D1["y5"] - pr["candidate"]) ** 2, np.nan)
    r2b = np.where(w5, (D1["y5"] - pr["B0"]) ** 2, np.nan)
    per_c = pd.Series(r2c).groupby(le_b).mean()
    per_b = pd.Series(r2b).groupby(le_b).mean()
    info.append({"item": "share of TEST link-events where the candidate beats B0", "value": float((per_c < per_b).mean())})
    res_ = np.where(w5, D1["y5"] - pr["candidate"], np.nan)
    lk = le_tr.cml_id.to_numpy()[le_b]
    mb = pd.Series(res_).groupby(lk).mean()
    info.append({"item": "share of TRAIN links with |mean residual| < 0.3 dB on TEST (all minutes)", "value": float((mb.abs() < 0.3).mean())})
    rain_b = np.nan_to_num(D1["r5"], nan=-1) >= 0.1
    mbr = pd.Series(np.where(rain_b, res_, np.nan)).groupby(lk).mean().dropna()
    info.append({"item": f"share of TRAIN links with |mean residual| < 0.3 dB on TEST (rain minutes only, {len(mbr)} links)", "value": float((mbr.abs() < 0.3).mean())})
    for inp in ("V0", "V1", "V3"):
        Dx = T[inp][0]
        info.append({"item": f"candidate RMSE with input {inp}", "value": pooled(Dx, predict(Dx, le_tr, f0, f0["M"])[0])})
    hasg = (le_tr.gauge_id.astype(str) != "").to_numpy()
    if hasg.any():
        Dg = sid.subset(D1, hasg)
        leg = le_tr[hasg].reset_index(drop=True)
        Dgg = dict(Dg)
        Dgg["rain"] = np.where(np.isfinite(Dg["rain_gauge"]), Dg["rain_gauge"], 0.0)
        info.append({"item": f"gauge links ({hasg.sum()} link-events): candidate RMSE with V1", "value": pooled(Dg, predict(Dg, leg, f0, f0["M"])[0])})
        info.append({"item": "gauge links: candidate RMSE with the gauge rain as input", "value": pooled(Dgg, predict(Dgg, leg, f0, f0["M"])[0])})
    if len(dryt):
        q = dryt[(dryt.length == "all")].set_index("t_min")
        for t in (5, 10):
            info.append({"item": f"dry-path median error at {t} min, M0 (candidate) / M1 (box B)", "value": f"{q.err_M0[t]:+.3f} / {q.err_M1[t]:+.3f} dB ({int(q.windows[t])} windows)"})
    dr = pd.read_csv(P.DEW_RISK_HOURLY, parse_dates=["time_utc"]).set_index("time_utc").dew_risk
    fwb = T5.searchsorted(pd.to_datetime(le_tr.first_wet).to_numpy())
    dew_flag = []
    for e, r in enumerate(le_tr.itertuples()):
        fw = pd.Timestamp(r.first_wet)
        hrs = (fw - pd.Timedelta("6h")).floor("h") + pd.to_timedelta(np.arange(6), "h")
        rmx = Rm[max(fwb[e] - 72, 0):fwb[e], col[r.cml_id]]
        dew_flag.append(bool(int(dr.reindex(hrs).fillna(False).astype(bool).sum()) >= 4 and len(rmx) == 72 and np.isfinite(rmx).all() and (rmx <= 0.05).all()))
    dew_flag = np.array(dew_flag)
    info.append({"item": "dew-before-rain link-events on TEST (S7.4 rule)", "value": f"{int(dew_flag.sum())} of {len(dew_flag)} ({le_tr[dew_flag].event_id.nunique()} events)"})
    if dew_flag.sum():
        bias = pd.Series(res_).groupby(le_b).mean()
        info.append({"item": "mean residual with / without dew before the rain", "value": f"{bias[dew_flag].mean():+.3f} / {bias[~dew_flag].mean():+.3f} dB"})
    infot = pd.DataFrame(info)
    infot["value"] = [f"{v:.3f}" if isinstance(v, float) else v for v in infot.value]
    bud = error_budget(D1, le_tr, pr["candidate"])
    ft, lct = funnel(d, drops_tr, drops_te, le_tr, le_te, Dv1)

    # save, figures, report
    rmt = pd.DataFrame(rm_rows)
    for tab, path in ((ft, P.S8_FUNNEL), (lct, P.S8_LENGTHS), (rg.round(4), P.S8_G_E), (rules, P.S8_RULES),
                      (pd.concat([rmt, p2t.assign(set="P2 differences"), rt.assign(set="P4 TEST / VAL"), nl.assign(set="P5")], ignore_index=True).round(4), P.S8_RMSE),
                      (cells.round(4), P.S8_CELLS), (pd.concat([dryt, win.round(4).assign(kind="window")], ignore_index=True), P.S8_DRY), (infot, P.S8_INFO),
                      (bud.round(4), P.S8_BUDGET)):
        save(tab, path if not dry else DRY / path.name, index=False)
    figs = [(fig_rmse(rmt), P.S8_FIG_RMSE), (fig_cells(cells), P.S8_FIG_CELLS), (fig_dry(win, curves), P.S8_FIG_DRY), (fig_ratio(rt), P.S8_FIG_RATIO),
            (fig_new(nl.iloc[:2]), P.S8_FIG_NEW), (fig_budget(bud), P.S8_FIG_BUDGET)]
    for fig, path in figs:
        save(fig, path if not dry else DRY / path.name, dpi=100)
        plt.close(fig)
    q = dict(locals())
    report = write_report(q)
    save(report, P.S8_REPORT if not dry else DRY / P.S8_REPORT.name)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {(time.time() - clock) / 60:.1f} min")


def readings(q):
    """The verdict and readings, generated from the numbers (written before TEST was read)."""
    R = {}
    c = q["cmp_t"]
    v = c[c.set.str.startswith("VAL")].set_index("model").rmse_pooled
    t = c[c.set == "TRAIN"].set_index("model").rmse_pooled
    R["READ2"] = ("**READING.** Pooled 5 min RMSE on TRAIN (fitted) and VAL (used, optimistic): "
                  + "; ".join(f"{m} {t[m]:.3f} / {v[m]:.3f} dB" for m in t.index)
                  + f". The candidate is {'the best' if v.idxmin() == 'candidate' else 'NOT the best (' + v.idxmin() + ' is)'} on VAL. These numbers are a sanity "
                  "check only; nothing was tuned after them.")
    rules = q["rules"]
    failed = rules[rules.verdict == "FAIL"]
    rmt = q["rmt"].set_index(["subset", "model"]).rmse
    head = ("**The box-A equation holds on TEST**: P1-P5 pass" + (" (P3 not testable)" if (rules.verdict == "not testable").any() else "") + "."
            if q["holds"] else "**The box-A equation does NOT fully hold on TEST**: " + "; ".join(f"{r.rule} failed ({r.numbers}) - this touches: {r.claim}"
                                                                                                for r in failed.itertuples()) + ".")
    R["READ6"] = (head + f" On TEST the candidate's pooled RMSE is {rmt[('all', 'candidate')]:.3f} dB against B0 {rmt[('all', 'B0')]:.3f}, B1 "
                  f"{rmt[('all', 'B1')]:.3f}, B2 {rmt[('all', 'B2')]:.3f} dB (L < 1 km: {rmt[('L < 1 km', 'candidate')]:.3f} against "
                  f"{rmt[('L < 1 km', 'B0')]:.3f} / {rmt[('L < 1 km', 'B1')]:.3f} / {rmt[('L < 1 km', 'B2')]:.3f} dB). The length test: {q['n_ok']} of "
                  f"{q['n_el']} cells. The TEST / VAL ratio is {q['allr']:.2f} overall. TEST links (never fitted) against TRAIN links: {q['ratio5']:.2f} x. "
                  "Nothing was refitted or changed after TEST was read.")
    R["READ7"] = ("1. The equation (box A) with its numbers and intervals (S6), and what it does NOT claim (boxes B and C, S6b).\n"
                  "2. The evidence that W does not depend on L: S7 (four checks) and TEST (P1, P3).\n"
                  "3. The value of the antenna term and of its memory: TEST against B0, B1, B2 (P2), all links and short links.\n"
                  "4. The data funnel (S6 TRAIN / VAL and S8 TEST) and the error budget (VAL and TEST).\n"
                  "5. The decisions and deviations, dated (M4c -> M4b-B, the 2 / 3 threshold, M1 -> box B, the undone cuts), and the lessons.\n"
                  "6. The rain reference: V1 = radar / g_e from the City gauges, and its limits (the radar is gauge-adjusted; g_e for few events).")
    R["READ8"] = ("- **TEST baselines were computed in this step** with the frozen baseline code; the TEST links have no fake-event uncertainty (the "
                  "offset_std rule does not apply to them) and no daily-cycle amplitude (that flag is not raised).\n"
                  f"- **g_e** is used for {int(q['rg'].path_used.sum())} of {len(q['rg'])} TEST events; the others keep the raw radar (g_e = 1).\n"
                  "- **P5 compares TEST links with TRAIN links on the same events**, but their length and band mix differ (the per-class ratios are in "
                  "the P5 table).\n"
                  "- **P2 uses the pooled RMSE**, dominated by the rain minutes (S6b budget: during rain 68 % of the squared error on VAL).\n"
                  "- **VAL was used for decisions**; TEST was read once, after the freeze.")
    return R


def write_report(q):
    q = {**q, **readings(q)}
    L = []
    a = L.append
    a("# wet_antenna step S8: freeze, then TEST once\n")
    a("Written by `research/wet_antenna/test_s8.py` (libraries `final_model.py`, `validation.py`, `system_id.py`, `length_test.py`, `event_factors.py`, "
      "`wind_drying.py`, `baseline/baseline.py`). Every rule: `wet_antenna.test_s8` in `research/config.yaml`, written on 2026-10-19 BEFORE TEST was read.\n")
    a(STANDING + "\n")
    a("The candidate is **M0 = the S6 candidate** (decision 2026-10-19 in `docs/ROADMAP.md`; M1 of S6b is box B and is shown for information only).\n")
    fr = q["freeze"]
    a("## 1. The freeze\n")
    a(f"`{rel(P.S8_FREEZE)}` (written {fr['written']}, committed before TEST was read; checked by this run):\n")
    a(f"- candidate `{fr['candidate']['file']}`: sha256 `{fr['candidate']['sha256']}`")
    for b, v in fr["comparison_fits"].items():
        a(f"- comparison fit {b} `{v['file']}`: sha256 `{v['sha256']}`")
    a(f"- M1 (box B, information) `{fr['m1_box_b_fit']['file']}`: sha256 `{fr['m1_box_b_fit']['sha256']}`")
    a(f"- config block `wet_antenna.test_s8`: sha256 `{fr['config_block_test_s8']}`")
    a(f"- code commit `{fr['code_commit']}`; code changed since then: {', '.join(q['code_diff']) or 'none'}\n")
    a("## 2. The comparison models (fitted on TRAIN only)\n")
    a(f"B3 skipped: {S['comparison']['B3']}.\n")
    a(md_table(q["cmp_t"].round(4)) + "\n")
    a(q.get("READ2", "") + "\n")
    a("## 3. The TEST funnel\n")
    a(md_table(q["ft"]) + "\n")
    a(md_table(q["lct"]) + "\n")
    rg = q["rg"]
    a(f"g_e of the TEST events (S4.3 path-mean rule): used for {int(rg.path_used.sum())} of {len(rg)} events; pooled V3 ratio {q['G']:.3f} (`{rel(P.S8_G_E)}`).\n")
    a(md_table(rg[["event_id", "gauge_links", "gauge_mm_per_link", "radar_over_gauge", "path_used", "g_path"]].round(3)) + "\n")
    a("## 4. P1-P5 on TEST\n")
    a(md_table(q["rules"]) + "\n")
    a("RMSE (pooled 5 min bins, 90 % event bootstrap):\n")
    a(md_table(q["rmt"].round(4)) + "\n")
    a("P2 differences (Bx minus candidate):\n")
    a(md_table(q["p2t"].round(4)) + "\n")
    a("![rmse](../figures/test_s8_rmse.png)\n")
    a("P1 cells:\n")
    a(md_table(q["cells"][["band", "class", "R_cell", "links_under_1km", "eligible", "I_obs", "I_lo", "I_hi", "I_model", "s_obs", "s_lo", "s_hi", "s_model",
                           "I_inside", "s_inside"]].round(3)) + "\n")
    a("![cells](../figures/test_s8_length_cells.png)\n")
    a("P3 dry-path errors (W_obs - W_model, median, 90 % event bootstrap):\n")
    a(md_table(q["dryt"].round(3)) + "\n" if len(q["dryt"]) else "No TEST dry-path windows.\n")
    a("![dry path](../figures/test_s8_dry_path.png)\n")
    a("P4 TEST / VAL:\n")
    a(md_table(q["rt"].round(4)) + "\n")
    a("![ratio](../figures/test_s8_test_over_val.png)\n")
    a("P5 new links:\n")
    a(md_table(q["nl"].round(4)) + "\n")
    a("![new links](../figures/test_s8_new_links.png)\n")
    a("## 5. Information items\n")
    a(md_table(q["infot"]) + "\n")
    a("The TEST error budget (candidate, V1; the S6b split):\n")
    a(md_table(q["bud"].round(4)) + "\n")
    a("![budget](../figures/test_s8_error_budget.png)\n")
    a("## 6. The verdict\n")
    a(q.get("READ6", "") + "\n")
    a("## 7. What S9 should present\n")
    a(q.get("READ7", "") + "\n")
    a("## 8. Unclear / possibly wrong\n")
    a(q.get("READ8", "") + "\n")
    return "\n".join(L)


if __name__ == "__main__":
    main()
