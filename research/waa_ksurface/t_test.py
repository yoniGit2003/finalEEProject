"""The final TEST read (docs/waa_ksurface/test_prereg.md): the frozen equation against the existing models on TEST.

    python -m research.waa_ksurface.t_test --dryrun-train   the whole path on TRAIN (scratch cache only; not a result)
    python -m research.waa_ksurface.t_test --freeze         the freeze record t_test_freeze.json (before TEST)
    python -m research.waa_ksurface.t_test --test           check every hash, read TEST ONCE, predict, all metrics
    python -m research.waa_ksurface.t_test --report         the report from the saved results (no data read)

Sets: known links (TEST events, TRAIN links), new links (TEST events, the never-fitted TEST links; the Part 1 models
with a_link = 1), gauge-ratio events (TEST events with a gauge ratio, TRAIN links). Models, metrics and readings as Block B
(b_val.py, whose functions are used unchanged).
"""
import json
import pickle
import sys
import time

import numpy as np
import pandas as pd

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a7_review as a7r, b_val as B, data, fake6, ksurf, s1_separability as s1
from research.waa_ksurface import paths as P
from research.waa_ksurface.e2_candidates import sizes_and_rule
from research.wet_antenna import paths as WP, system_id as sid
from research.wet_antenna.event_factors import rescaled

CB = data.branch_config()["b_val"]                            # the Block B settings: seed 4100, 1000 draws, interval, M2 lags
CACHE = P.CACHE / "t_test"
log = B.log
MODELS = B.MODELS
SETS = ("known links", "new links", "gauge-ratio events")
FREEZE = P.PRODUCTS / "t_test_freeze.json"
PREREG = P.REPO / "docs" / "waa_ksurface" / "test_prereg.md"
OUT = {"m1": P.TABLES / "t_test_m1_rmse.csv", "m1p": P.TABLES / "t_test_m1_paired.csv", "m2": P.TABLES / "t_test_m2_after_rain.csv",
       "m3": P.TABLES / "t_test_m3_length_test.csv", "report": P.REPORTS / "t_test.md"}


# =================================================================================================== the data
def load(split):
    """{"known": (D, le), "new": (D, le)}, drops, gmap, ratio events. split = "test" (TEST read) or "train" (the dry run:
    TRAIN events; "new links" played by 20 % of the TRAIN links, seed 0)."""
    if split == "test":
        d = data.open_release()
        bands_cfg = sid.merged_bands(d, data.DATA)[0]
        tb = pd.read_csv(WP.S8_BASELINES)
        cfg = {**data.DATA, "bands_ghz": bands_cfg, "max_length_km": data.S3["data"]["max_length_km_all"]}
        raw, drops = {}, {}
        for name, ls in (("known", "train"), ("new", "test")):
            D, le, dr = sid.build_dataset(d, cfg, event_split="test", link_split=ls, baselines=tb)
            raw[name], drops[name] = (D, le.reset_index(drop=True)), dr
        g = pd.read_csv(WP.S8_G_E)
        gmap = {r.event_id: (r.g_path if r.path_used else 1.0) for r in g.itertuples()}
        ratio = set(g[g.path_used].event_id)
    else:
        D, _, le, dr, _ = data.train()
        le = le.reset_index(drop=True)
        links = np.sort(le.cml_id.unique())
        new_links = set(np.random.default_rng(0).choice(links, len(links) // 5, replace=False))
        is_new = le.cml_id.isin(new_links).to_numpy()
        raw = {"known": (sid.subset(D, ~is_new), le[~is_new].reset_index(drop=True)), "new": (sid.subset(D, is_new), le[is_new].reset_index(drop=True))}
        drops = {"known": dr, "new": dr}
        g = data.g_table()
        gmap = {r.event_id: (r.g_path if r.path_used else 1.0) for r in g.itertuples()}
        ratio = set(g[g.path_used].event_id)
    out, info = {}, {}
    for name, (D, le) in raw.items():
        keep = le.band_name.isin(data.BRANCH_BANDS).to_numpy()
        n_band = int(keep.sum())
        if name == "known":
            keep = keep & B.predictable(le)                         # Block B drop rule: a link without a stored a_link
        le_k = le[keep].reset_index(drop=True)
        out[name] = (rescaled(sid.subset(D, keep), le_k, {e: gmap.get(e, 1.0) for e in le_k.event_id.unique()}), le_k)
        info[name] = {"built": len(le), "at_28_38": n_band, "dropped_no_a_link": n_band - int(keep.sum()), "kept": int(keep.sum()),
                      "links": int(le_k.cml_id.nunique()), "events": int(le_k.event_id.nunique())}
    return out, drops, info, ratio


def predict(D, le, new):
    """Every model on one dataset (A - c per 5 min). new = the Part 1 models get a_link = 1 (no stored factor)."""
    curves, band_idx = B.band_setting(le)
    rec = json.loads(P.B_RECIPES.read_text(encoding="utf-8"))
    pr = {}
    prob = ksurf.prepare(D, band_idx, curves, D["rain"], 0)
    m, x = B.ours_model(prob)
    pr["ours"] = m.simulate(x)[1]
    del prob, m
    da = 1.0 if new else None
    for n, f in B.lit_fits().items():
        if new:
            f = {**f, "a_link_per_link": {}}                    # a_link = 1 on a never-fitted link
        pr[n] = B.lit_predict(D, le, f, default_a=da)
    pr["Part 1 equation"] = B.part1_predict(D, le, default_a=da) if not new else part1_new(D, le)
    for k in B.RECIPES:
        pr[k] = B.recipe_predict(D, le, k, rec["recipes"][k]["increments"], curves, band_idx)
    return pr


def part1_new(D, le):
    """The Part 1 equation with a_link = 1 on every link (none of the new links has a stored factor)."""
    from research.wet_antenna.waa_equation import WAAEquation
    eq = WAAEquation()
    rain = np.nan_to_num(np.asarray(D["rain"], float))
    out = np.empty(int(D["le_n"].sum()))
    for e in range(len(le)):
        o, n, o1 = D["le_roff"][e], D["le_n"][e], D["le_boff1"][e]
        r = eq.run(rain[o + B.PADM:o + B.PADM + n], float(D["le_L"][e]), le.band_name.iloc[e], a_link=1.0, rain_before=rain[o + B.PADM - eq.lag:o + B.PADM])
        out[o1:o1 + n] = r["A_minus_c"]
    return out.reshape(-1, 5).mean(axis=1)


# =================================================================================================== metrics
def metrics(parts, ratio, t0):
    """parts[name] = dict(D, le, preds). Sets: known, new, gauge-ratio (known links in ratio events)."""
    all_events = sorted(set().union(*[set(p["le"].event_id) for p in parts.values()]))
    code = {e: i for i, e in enumerate(all_events)}
    E = len(all_events)
    rng = np.random.default_rng(CB["seed"])
    mult = np.array([np.bincount(rng.integers(0, E, E), minlength=E) for _ in range(CB["boot_reps"])], float)
    lo_q, hi_q = CB["interval"]
    for p in parts.values():
        D, le = p["D"], p["le"]
        p["bin_le"] = np.repeat(np.arange(len(le)), D["le_n"] // 5)
        p["target"] = D["y5"] - D["le_c"][p["bin_le"]]
        p["scored"] = (D["w5"] > 0) & np.isfinite(p["target"]) & np.logical_and.reduce([np.isfinite(v) for v in p["preds"].values()])
        p["ev_bin"] = le.event_id.map(code).to_numpy()[p["bin_le"]]
        p["L_bin"] = le.L_km.to_numpy()[p["bin_le"]]
        p["ratio_bin"] = le.event_id.isin(ratio).to_numpy()[p["bin_le"]]
    sets = {"known links": ("known", None), "new links": ("new", None), "gauge-ratio events": ("known", "ratio_bin")}
    m1, m1p, m2 = [], [], []
    for sn, (pn, extra) in sets.items():
        p = parts[pn]
        base = p["scored"] & (p[extra] if extra else True)
        for gn, gm in (("all links", np.ones(len(base), bool)), ("links < 1 km", p["L_bin"] < 1.0)):
            msk = base & gm
            bt = {}
            for mn in MODELS:
                r2 = (p["preds"][mn] - p["target"]) ** 2
                S_ = np.bincount(p["ev_bin"][msk], weights=r2[msk], minlength=E)
                N_ = np.bincount(p["ev_bin"][msk], minlength=E).astype(float)
                bt[mn] = B.boot_rmse(S_, N_, mult)
                m1.append({"set": sn, "links": gn, "model": mn, "RMSE (dB)": float(np.sqrt(S_.sum() / max(N_.sum(), 1))), "values": int(N_.sum()),
                           "events": int((N_ > 0).sum()), "link-events": int(len(np.unique(p["bin_le"][msk]))),
                           "90 % (dB)": f"{np.percentile(bt[mn], lo_q):.4f} to {np.percentile(bt[mn], hi_q):.4f}"})
            r_o = m1[-len(MODELS)]["RMSE (dB)"]
            for k, mn in enumerate(MODELS[1:], 1):
                lo, hi = np.percentile(bt[mn] - bt["ours"], [lo_q, hi_q])
                m1p.append({"set": sn, "links": gn, "model": mn, "RMSE(model) - RMSE(ours) (dB)": m1[-len(MODELS) + k]["RMSE (dB)"] - r_o,
                            "90 % lo": lo, "90 % hi": hi,
                            "reading": f"{mn} better than ours" if hi < 0 else (f"ours better than {mn}" if lo > 0 else "tie")})
        # M2 (a check)
        lags = CB["m2_lags_min"]
        D = p["D"]
        idx, _ = a7r.dry_values(D, fake6.window(D, np.asarray(D["rain"], float)), lags, CB["m2_rain_before_min"], s1.S1["dry_below_mmh"], np.ones(len(p["le"])))
        for i, lag in enumerate(lags):
            b = idx[:, i][idx[:, i] >= 0]
            b = b[base[b]]
            eb = p["ev_bin"][b]
            row = {"set": sn, "lag (min)": lag, "values": int(len(b)), "events": int(len(np.unique(eb))),
                   "measured median (dB)": float(np.median(p["target"][b])) if len(b) else np.nan}
            for mn in MODELS:
                if not len(b):
                    row[mn] = "n/a"
                    continue
                e0 = np.median(p["preds"][mn][b]) - np.median(p["target"][b])
                dr = [np.median(p["preds"][mn][rw]) - np.median(p["target"][rw]) for w in mult for rw in [np.repeat(b, w[eb].astype(int))] if len(rw)]
                lo, hi = np.percentile(dr, [lo_q, hi_q])
                row[mn] = f"{e0:+.3f} ({lo:+.3f} to {hi:+.3f})"
            m2.append(row)
    log(f"=== M1, M2 done ({(time.time() - t0) / 60:.1f} min) ===")
    # M3
    _, rule, rule_name, n_short, _ = sizes_and_rule()
    e4 = pd.read_csv(P.E4_CELLS)
    m3 = []
    for sn, (pn, extra) in sets.items():
        p = parts[pn]
        ev_keep = np.ones(len(p["le"]), bool) if extra is None else p["le"].event_id.isin(ratio).to_numpy()
        Ds = sid.subset(p["D"], ev_keep)
        les = p["le"][ev_keep].reset_index(drop=True)
        bk = ev_keep[p["bin_le"]]
        ps = {mn: np.where(p["scored"][bk], p["preds"][mn][bk], np.nan) for mn in MODELS}
        Ds["le_event"] = pd.factorize(les.event_id)[0]
        tab, ncell = B.m3_cells(Ds, les, ps, rule, n_short, e4[e4.eligible])
        m3.append(tab.assign(set=sn) if len(tab) else pd.DataFrame([{"set": sn, "model": "(none)", "inside": False, "cells": 0}]).iloc[0:0])
        log(f"   (M3) {sn}: {ncell} eligible cells")
    m3 = pd.concat(m3, ignore_index=True) if m3 else pd.DataFrame()
    log(f"=== M3 done ({(time.time() - t0) / 60:.1f} min) ===")
    return {"m1": pd.DataFrame(m1), "m1p": pd.DataFrame(m1p), "m2": pd.DataFrame(m2), "m3": m3, "events": E, "rule": rule_name}


def run(split):
    t0 = time.time()
    if split == "test":
        bad = B.check_freeze(FREEZE)
        if bad:
            log("=== STOP: a frozen file changed: " + "; ".join(bad) + " ===")
            return None
        log("=== every hash of t_test_freeze.json checked: unchanged ===")
    sets, drops, info, ratio = load(split)                     # ===== TEST is read here, once (split "test") =====
    log(f"=== {split.upper()} read: " + "; ".join(f"{k}: {v}" for k, v in info.items()) + f" ({(time.time() - t0) / 60:.1f} min) ===")
    parts = {}
    for name, (D, le) in sets.items():
        parts[name] = {"D": D, "le": le, "preds": predict(D, le, new=(name == "new"))}
    pre = {k: {"le": v["le"], "preds": v["preds"], "y5": v["D"]["y5"], "le_c": v["D"]["le_c"], "w5": v["D"]["w5"]} for k, v in parts.items()}
    (CACHE / f"predictions_{split}.pkl").write_bytes(pickle.dumps(pre))    # saved before any metric
    log(f"=== every model predicted and saved ({(time.time() - t0) / 60:.1f} min) ===")
    r = metrics(parts, ratio, t0)
    r.update({"info": info, "drops": drops, "ratio_events": sorted(ratio & set(parts["known"]["le"].event_id)), "minutes": (time.time() - t0) / 60, "split": split})
    return r


# =================================================================================================== report
def report(r, out):
    for k in ("m1", "m1p", "m2", "m3"):
        save(r[k].round(5) if k != "m2" else r[k], out[k], index=False)
    info = pd.DataFrame([{"set": k, **v} for k, v in r["info"].items()])
    m3s = (r["m3"].groupby(["set", "model"]).inside.agg(["sum", "size"]).reset_index().rename(columns={"sum": "cells inside", "size": "cells"})
           if len(r["m3"]) else pd.DataFrame([{"result": "no cell is both an E4 cell and eligible on the measurements of any set"}]))
    L = [f"# The final TEST read: the frozen equation against the existing models{' (DRY RUN ON TRAIN: not a result)' if r['split'] != 'test' else ''}\n",
         "Written by `research/waa_ksurface/t_test.py` under `docs/waa_ksurface/test_prereg.md` (hashes in `t_test_freeze.json`). TEST read once by this "
         "branch (the second read of TEST overall; the first was Part 1's, for its own model). Nothing refitted.\n",
         "## Data\n", md_table(info), "", f"Events: {r['events']}; with a gauge ratio (known links): {len(r['ratio_events'])}. Run time {r['minutes']:.1f} min.\n",
         "## M1: RMSE of the 5 min A - c (paired event bootstrap, 1000 draws, seed 4100)\n", md_table(r["m1"].round(4)), "",
         "### Paired differences and the pre-registered readings (the only claims)\n", md_table(r["m1p"].round(4)), "",
         "## M2: after-rain error (a CHECK; dry-shape caveat; no claim)\n", md_table(r["m2"]), "",
         f"## M3: the length test (rule {r['rule']}) on every model's predictions (no claim)\n", md_table(m3s), ""]
    save("\n".join(L) + "\n", out["report"])
    log(f"report: {rel(out['report'])}")


def freeze():
    files = {"TEST rules": PREREG, "TEST code": P.REPO / "research/waa_ksurface/t_test.py", "Block B code": P.REPO / "research/waa_ksurface/b_val.py",
             "frozen equation": P.B_EQUATION, "recipe fits": P.B_RECIPES,
             "Block B record 1": P.B_FREEZE, "Block B record 2": P.B_FREEZE2, "Block B record 3": P.B_FREEZE3, "Block B record 4": P.PRODUCTS / "a8_block_b_freeze4.json",
             "Part 1 equation": P.P1_EQUATION, "Part 1 literature fits": P.P1_LITERATURE_FITS, "Part 1 TEST baselines": WP.S8_BASELINES,
             "Part 1 TEST g_e": WP.S8_G_E, "ksurf.py": P.REPO / "research/waa_ksurface/ksurf.py", "a8_forms.py": P.REPO / "research/waa_ksurface/a8_forms.py",
             "fp.py": P.REPO / "research/waa_ksurface/fp.py", "E4 cells": P.E4_CELLS}
    rec = {"date": time.strftime("%Y-%m-%d %H:%M:%S"), "note": "TEST freeze record (test_prereg.md); written and committed before TEST is read",
           "sha256": {k: {"path": rel(v), "sha256": B.sha(v)} for k, v in files.items()}}
    save(json.dumps(rec, indent=1), FREEZE)
    for k, v in rec["sha256"].items():
        log(f"   {v['sha256'][:16]}  {v['path']}  ({k})")


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    if "--freeze" in sys.argv:
        freeze()
        return
    if "--dryrun-train" in sys.argv:
        d = CACHE / "dryrun"
        d.mkdir(exist_ok=True)
        r = run("train")
        report(r, {k: d / v.name for k, v in OUT.items()})
        log(f"=== dry run on TRAIN ended ({r['minutes']:.1f} min); written to {d} ===")
        return
    if "--test" in sys.argv:
        f = CACHE / "test.pkl"
        if f.exists() or (CACHE / "predictions_test.pkl").exists():
            log("=== TEST was already read: not read again (use --report) ===")
            sys.exit(2)
        r = run("test")
        if r is None:
            sys.exit(1)
        f.write_bytes(pickle.dumps(r))                          # every result saved before the report
        report(r, OUT)
        log(f"=== TEST ended {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
        return
    if "--report" in sys.argv:
        report(pickle.loads((CACHE / "test.pkl").read_bytes()), OUT)


if __name__ == "__main__":
    main()
