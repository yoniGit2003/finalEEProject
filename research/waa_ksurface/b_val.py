"""Block B: the frozen equation against the existing models on VAL (docs/waa_ksurface/block_b_prereg.md, frozen in 37470fa).

    python -m research.waa_ksurface.b_val --recipes     (before VAL) code checks on TRAIN, then fit the two recipes on TRAIN
    python -m research.waa_ksurface.b_val --freeze      (before VAL) hash the recipe fits and this code (second freeze record)
    python -m research.waa_ksurface.b_val --val         check every hash, then read VAL ONCE, predict every model, all metrics
    python -m research.waa_ksurface.b_val --report      the report from the saved VAL results (no data read)

Models: ours (the frozen equation file, the unchanged ksurf kernel), Schleiss (b) and Pastorek (b) (Part 1's stored TRAIN
fits, Part 1's forward code, read only), the Part 1 equation (its frozen file, Part 1's WAAEquation, read only), and two
recipes of this branch fitted on TRAIN: no WAA (A - c = L P_full(R)) and no memory (A - c = F_full(R) + L P_full(R)).
Rain input R = radar path mean / g_e for every model. Metrics M1 (RMSE, paired event bootstrap), M2 (after-rain check),
M3 (length test on every model's predictions). TEST is not opened.
"""
import hashlib
import json
import pickle
import sys
import time

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a7_review as a7r, a8_name as a8, data, fake6, fakeworld as fw, fp, ksurf, ltest, s1_separability as s1
from research.waa_ksurface import paths as P
from research.waa_ksurface.e2_candidates import sizes_and_rule
from research.wet_antenna import literature_s8c as lit, system_id as sid
from research.wet_antenna.event_factors import rescaled
from research.wet_antenna.waa_equation import WAAEquation

CB = data.branch_config()["b_val"]
CACHE = P.CACHE / "b_val"
log = a8.log
PADM = sid.PADM
MODELS = ["ours", "Schleiss b", "Pastorek b", "Part 1 equation", "no WAA", "no memory"]
RECIPES = ("no WAA", "no memory")
LTC = ltest.config()


# =================================================================================================== hashes
def sha(path):
    return a8.sha256(path)                                      # line endings normalized to LF, as the step 8 freeze


def check_freeze(record):
    """Every file of a freeze record must have its recorded hash. Returns a list of mismatches."""
    rec = json.loads(record.read_text(encoding="utf-8"))
    bad = []
    for name, v in rec["sha256"].items():
        h = sha(P.REPO / v["path"])
        if h != v["sha256"]:
            bad.append(f"{name}: {v['path']} {h[:16]} != {v['sha256'][:16]}")
    return bad


# =================================================================================================== our frozen equation
def ours_model(prob):
    eq = json.loads(P.B_EQUATION.read_text(encoding="utf-8"))
    c = dict(eq["ksurf_config_a5"])
    c["w_knots_dry_db"] = [float(v) for v in eq["w_knots_dry_db"]]
    m = ksurf.Model(prob, "K0", "K1", 0.0, cfg=c)
    x = np.asarray(eq["ksurf_parameter_vector"], float)
    assert len(x) == m.n, f"frozen vector has {len(x)} values, the model {m.n}"
    return m, x


def band_setting(le):
    M = fp.measured()
    bands = list(M)
    return [M[b] for b in bands], le.band_name.map({b: i for i, b in enumerate(bands)}).to_numpy()


# =================================================================================================== the recipes (no WAA, no memory)
def recipe_columns(D, le, curves, band_idx, kind, cfg=None):
    """5 min model output = base + C @ inc. The model is instantaneous and affine in the extension increments, so the 5 min
    mean of every column is exact. Parameter order: [F incs per band (no memory only)] + [P incs per band]."""
    cfg = cfg or ksurf.config()
    n = D["le_n"]
    idx = np.concatenate([o + PADM + np.arange(k) for o, k in zip(D["le_roff"], n)])
    rep = np.repeat(np.arange(len(n)), n)
    R = np.nan_to_num(np.asarray(D["rain"], float)[idx])       # as ksurf.prepare: the rain that drives the model
    L = D["le_L"][rep]
    b1 = band_idx[rep]
    parts = (["F"] if kind == "no memory" else []) + ["P"]
    nF, nP = len(cfg["f_knots_mmh"]), len(cfg["p_knots_mmh"])
    base = np.zeros(len(R))
    cols, names = [], []
    for q in parts:
        nk = nF if q == "F" else nP
        for b, cv in enumerate(curves):
            m = b1 == b
            mult = 1.0 if q == "F" else L[m]
            z = np.zeros(nk)
            base[m] += mult * ksurf.full_curve(cv, q, z, R[m], cfg)
            for j in range(nk):
                e = z.copy()
                e[j] = 1.0
                col = np.zeros(len(R))
                col[m] = mult * (ksurf.full_curve(cv, q, e, R[m], cfg) - ksurf.full_curve(cv, q, z, R[m], cfg))
                cols.append(col.reshape(-1, 5).mean(axis=1))
                names.append(f"{q} band {b} inc {j}")
    return base.reshape(-1, 5).mean(axis=1), np.column_stack(cols), names


def recipe_bounds(names, cfg=None):
    cfg = cfg or ksurf.config()
    return [(0.0, cfg["f_inc_max_db"]) if n.startswith("F") else (0.0, cfg["p_inc_max_db_km"]) for n in names]


def huber_fit(base, C, y, w, x0, bounds, delta):
    """Mean Huber loss of (base + C x - y) over the values with w > 0 (L-BFGS-B, exact gradient)."""
    m = w > 0
    b, Cm, ym = base[m], C[m], y[m]

    def f(x):
        r = b + Cm @ x - ym
        a = np.abs(r)
        loss = np.where(a <= delta, 0.5 * r * r, delta * (a - 0.5 * delta)).mean()
        return float(loss), Cm.T @ np.clip(r, -delta, delta) / len(r)
    r = minimize(f, np.clip(x0, [lo for lo, _ in bounds], [hi for _, hi in bounds]), jac=True, method="L-BFGS-B", bounds=bounds,
                 options={"maxiter": CB["recipe_maxiter"], "ftol": 1e-12, "gtol": 1e-9})
    return r.x, float(r.fun), bool(r.success), int(r.nit)


def recipe_predict(D, le, kind, inc, curves, band_idx):
    base, C, _ = recipe_columns(D, le, curves, band_idx, kind)
    return base + C @ np.asarray(inc, float)


# =================================================================================================== Part 1 competitors (read-only code)
def lit_fits():
    J = json.loads(P.P1_LITERATURE_FITS.read_text(encoding="utf-8"))["fits"]
    return {n: J[n] for n in ("Schleiss b", "Pastorek b")}


def lit_predict(D, le, fit, default_a=None):
    """Part 1's literature_s8c.Set.predict (its lag of 1 min, its W_max / C per band, its a_link per link); returns A - c.
    default_a = None: every link must have a stored a_link (the drop rule removes the others before)."""
    T = lit.Set(D, le)
    a = {str(k): v for k, v in fit["a_link_per_link"].items()}
    if default_a is None:
        assert all(str(i) in a for i in T.ids), "a link without a stored a_link reached the prediction"
    f = {"kind": fit["kind"], "shape": fit["shape"], "amp": fit["amp"], "c_band": fit["c_band"], "a": {i: a.get(str(i), default_a or 1.0) for i in T.ids}}
    pred5, _ = T.predict(f)
    return pred5 - T.c5


def part1_predict(D, le, default_a=None):
    """The Part 1 equation (WAAEquation, its file, its lag of 1 min, a_link per link) per link-event; A - c per 5 min."""
    eq = WAAEquation()
    a = {str(k): v for k, v in eq.eq["a_link"]["per_link"].items()}
    rain = np.nan_to_num(np.asarray(D["rain"], float))
    out = np.empty(int(D["le_n"].sum()))
    for e in range(len(le)):
        o, n, o1 = D["le_roff"][e], D["le_n"][e], D["le_boff1"][e]
        cid = str(le.cml_id.iloc[e])
        al = a[cid] if default_a is None else a.get(cid, default_a)
        r = eq.run(rain[o + PADM:o + PADM + n], float(D["le_L"][e]), le.band_name.iloc[e], a_link=al, rain_before=rain[o + PADM - eq.lag:o + PADM])
        out[o1:o1 + n] = r["A_minus_c"]
    return out.reshape(-1, 5).mean(axis=1)


def predictable(le):
    """A link-event every model can predict: its link has a stored a_link in Schleiss b, Pastorek b and the Part 1 equation."""
    keys = [set(map(str, f["a_link_per_link"])) for f in lit_fits().values()]
    keys.append(set(map(str, WAAEquation().eq["a_link"]["per_link"])))
    cid = le.cml_id.astype(str)
    return np.logical_and.reduce([cid.isin(k).to_numpy() for k in keys])


# =================================================================================================== before VAL: code checks and the recipes
def rmse_sub(D, pred_A, mask):
    """Part 1's pooled RMSE (test_s8b.pooled): counted bins, y5 against the prediction of A (= A - c + c)."""
    w = (D["w5"] > 0) & mask
    c5 = D["le_c"][D["bin_le"]]
    r2 = (D["y5"] - (pred_A + c5)) ** 2
    return float(np.sqrt(r2[w].sum() / w.sum())), int(w.sum())


def recipes_part():
    out = {"checks": {}}
    # (1) our frozen equation reproduces its step 8 training loss on TRAIN (same values, same kernel)
    S, prob, use, score, fold, fit = a8.setting()
    m, x = ours_model(prob)
    val = m.loss(x, use, penalty=False)[0]
    ref = pickle.loads((P.CACHE / "a8_name" / "plug_0.pkl").read_bytes())["est"]["radar"]["train_loss"]
    dev = abs(val / ref - 1)
    out["checks"]["ours TRAIN loss"] = {"value": val, "reference": ref, "rel_dev": dev, "ok": dev <= CB["reproduce_rtol"]}
    log(f"   (check) ours on TRAIN: loss {val:.10f}, step 8 plug-back {ref:.10f}, relative difference {dev:.1e}")

    # (2) the recipes on the same values (stage 2 of A7: unflagged learning events, stage2_w5, radar input)
    D = S["D"]
    curves, band_idx = band_setting(S["le"])
    sF, sP = a7r.starts(S)
    delta = ksurf.config()["huber_delta_db"]
    out["recipes"] = {}
    for kind in RECIPES:
        base, C, names = recipe_columns(D, S["le"], curves, band_idx, kind)
        x0 = np.r_[np.concatenate(sF) if kind == "no memory" else [], np.concatenate(sP)]
        inc, loss, conv, nit = huber_fit(base, C, prob["target"], prob["w5"], x0, recipe_bounds(names), delta)
        out["recipes"][kind] = {"increments": inc.tolist(), "names": names, "train_loss": loss, "converged": conv, "iterations": nit}
        log(f"   (recipe) {kind}: mean Huber {loss:.6f} on TRAIN stage 2; converged {conv} ({nit} it); increments {np.round(inc, 4).tolist()}")
    del S, prob, D

    # (3) the Part 1 wrappers reproduce Part 1's own TRAIN numbers (literature_s8c_rmse.csv, TRAIN, all / L < 1 km)
    D0, D1, le, _, _ = data.train()
    tab = pd.read_csv(P.P1_LIT_RMSE)
    tab = tab[tab.set == "TRAIN"].set_index(["subset", "model"])
    Lb = le.L_km.to_numpy()[D1["bin_le"]]
    preds = {"Schleiss b": lit_predict(D1, le, lit_fits()["Schleiss b"], default_a=1.0),
             "Pastorek b": lit_predict(D1, le, lit_fits()["Pastorek b"], default_a=1.0),
             "ours": part1_predict(D1, le, default_a=1.0)}             # Part 1 calls its equation "ours" in that table
    for name, pr in preds.items():
        for sub, mask in (("all", np.ones(len(Lb), bool)), ("L < 1 km", Lb < 1)):
            v, nb = rmse_sub(D1, pr, mask)
            r = tab.loc[(sub, name)]
            ok = bool(abs(v - r.rmse) <= CB["p1_rmse_tol_db"] and nb == r.bins)
            out["checks"][f"Part 1 {name} TRAIN {sub}"] = {"value": v, "reference": float(r.rmse), "bins": nb, "bins_ref": int(r.bins), "ok": ok}
            log(f"   (check) Part 1 {name} on TRAIN {sub}: RMSE {v:.4f} (table {r.rmse:.4f}); bins {nb} (table {r.bins}): {'ok' if ok else 'NOT REPRODUCED'}")
    out["ok"] = all(c["ok"] for c in out["checks"].values())
    return out


def write_recipes(r):
    doc = {"written_by": "research/waa_ksurface/b_val.py --recipes", "date": time.strftime("%Y-%m-%d %H:%M:%S"),
           "fitted_on": "TRAIN, stage 2 of A7 (unflagged learning events, stage2_w5, radar input / g_e), mean Huber 0.3 dB, L-BFGS-B from the A6 extension starts",
           "recipes": {k: {**v, "form": ("A - c = L P_full(R(t)), W = 0" if k == "no WAA" else "A - c = F_full(R(t)) + L P_full(R(t)) (W = its target every minute)")}
                       for k, v in r["recipes"].items()},
           "code_checks": r["checks"]}
    save(json.dumps(doc, indent=1, default=float), P.B_RECIPES)


def freeze3():
    """Deviation 16: b_val.py changed (one-line fix); hash it again with the two earlier records and the recipe fits."""
    files = {"recipe fits": P.B_RECIPES, "Block B code (Deviation 16 fix)": P.REPO / "research" / "waa_ksurface" / "b_val.py",
             "first freeze record": P.B_FREEZE, "second freeze record": P.B_FREEZE2}
    rec = {"date": time.strftime("%Y-%m-%d %H:%M:%S"), "note": "third freeze record of Block B (Deviation 16); written before VAL is read again",
           "sha256": {k: {"path": rel(v), "sha256": sha(v)} for k, v in files.items()}}
    save(json.dumps(rec, indent=1), P.B_FREEZE3)
    for k, v in rec["sha256"].items():
        log(f"   {v['sha256'][:16]}  {v['path']}  ({k})")


def freeze_recipes():
    files = {"recipe fits": P.B_RECIPES, "Block B code": P.REPO / "research" / "waa_ksurface" / "b_val.py",
             "first freeze record": P.B_FREEZE}
    rec = {"date": time.strftime("%Y-%m-%d %H:%M:%S"), "note": "second freeze record of Block B (block_b_prereg.md, section 3); written before VAL is read",
           "sha256": {k: {"path": rel(v), "sha256": sha(v)} for k, v in files.items()}}
    save(json.dumps(rec, indent=1), P.B_FREEZE2)
    for k, v in rec["sha256"].items():
        log(f"   {v['sha256'][:16]}  {v['path']}  ({k})")


# =================================================================================================== VAL
def val_data(split="val"):
    """VAL read once by the Part 1 builder (W1): event_split val, the TRAIN-link rules, the data-layer offset c; 28 and 38 GHz.
    split = "train": the dry run of Deviation 16 (the TRAIN link-events from the data-layer cache, TRAIN g_e); not a result."""
    if split == "train":
        D, _, le, drops, _ = data.train()
    else:
        d = data.open_release()
        bands_cfg = sid.merged_bands(d, data.DATA)[0]
        D, le, drops = sid.build_dataset(d, {**data.DATA, "bands_ghz": bands_cfg, "max_length_km": data.S3["data"]["max_length_km_all"]}, event_split="val")
    le = le.reset_index(drop=True)
    keep = le.band_name.isin(data.BRANCH_BANDS).to_numpy()
    g = pd.read_csv(P.P1_G_E)
    g = g[g.set == split.upper()]
    gmap = {r.event_id: (r.g_path if r.path_used else 1.0) for r in g.itertuples()}
    ratio_events = set(g[g.path_used].event_id)
    return D, le, drops, keep, gmap, ratio_events


def boot_rmse(S_, N_, mult):
    return np.sqrt(mult @ S_ / np.maximum(mult @ N_, 1))


def m3_cells(D, le, preds, rule, n_short, e4_cells):
    """The adopted length test (E4 candidate (ii)) on the measurements and on every model's predictions (same bins)."""
    z = fw.area_mean(D["r5"].astype(float), le)
    t_meas = ltest.steady_values(D, le, LTC, rule, sel=z, r_stat="mean")
    mult = ltest.multiplicities(int(D["le_event"].max()) + 1, LTC["boot_reps"], LTC["seed"])
    meas, _ = ltest.cell_fits(t_meas, LTC, n_short, mult, LTC["interval"], r_stat="mean")
    keys = set(zip(e4_cells.band, e4_cells.cls))
    ok = meas.eligible.to_numpy(bool) & np.array([k in keys for k in zip(meas.band, meas.cls)], bool)     # Deviation 16: numpy masks
    meas = meas[ok].reset_index(drop=True)
    c5 = D["le_c"][D["bin_le"]]
    rows = []
    for name, pr in preds.items():
        Dp = dict(D)
        Dp["y5"] = np.where(np.isfinite(D["y5"]), pr + c5, np.nan)
        cp, _ = ltest.cell_fits(ltest.steady_values(Dp, le, LTC, rule, sel=z, r_stat="mean"), LTC, n_short, None, LTC["interval"], r_stat="mean")
        cp = cp.set_index(["band", "cls"])
        for r in meas.itertuples():
            q = cp.loc[(r.band, r.cls)] if (r.band, r.cls) in cp.index else None
            I, s_ = (float(q.I), float(q.s)) if q is not None else (np.nan, np.nan)
            rows.append({"model": name, "band": r.band, "class": meas.loc[r.Index, "class"], "I measured": r.I,
                         "I 90 %": f"{r.I_lo:.3f} to {r.I_hi:.3f}", "I model": I, "s measured": r.s, "s 90 %": f"{r.s_lo:.4f} to {r.s_hi:.4f}", "s model": s_,
                         "inside": bool(r.I_lo <= I <= r.I_hi and r.s_lo <= s_ <= r.s_hi)})
    return pd.DataFrame(rows), len(meas)


def val_part(split="val"):
    t0 = time.time()
    bad = [] if split == "train" else check_freeze(P.B_FREEZE) + check_freeze(P.B_FREEZE3)     # Deviation 16: the first and the third record
    if bad:
        log("=== STOP: a frozen file changed: " + "; ".join(bad) + " ===")
        return {"stop": bad}
    log("=== every frozen hash checked (both freeze records): unchanged ===")
    rec = json.loads(P.B_RECIPES.read_text(encoding="utf-8"))
    assert all(c["ok"] for c in rec["code_checks"].values()), "the pre-VAL code checks did not pass"

    D0, le0, drops, keep, gmap, ratio_events = val_data(split)       # ===== VAL is read here, once =====
    log(f"=== VAL read: {len(le0)} link-events ({int(keep.sum())} at 28 / 38 GHz) of {le0.event_id.nunique()} events ({(time.time() - t0) / 60:.1f} min) ===")
    can = predictable(le0)
    keep_all = keep & can
    n_drop = int((keep & ~can).sum())
    le = le0[keep_all].reset_index(drop=True)
    D = rescaled(sid.subset(D0, keep_all), le, {e: gmap.get(e, 1.0) for e in le.event_id.unique()})
    del D0
    curves, band_idx = band_setting(le)
    n5 = D["le_n"] // 5
    bin_le = np.repeat(np.arange(len(le)), n5)
    c5 = D["le_c"][bin_le]
    target = D["y5"] - c5

    preds = {}
    prob = ksurf.prepare(D, band_idx, curves, D["rain"], 0)
    m, x = ours_model(prob)
    preds["ours"] = m.simulate(x)[1]
    del prob, m
    for n, f in lit_fits().items():
        preds[n] = lit_predict(D, le, f)
    preds["Part 1 equation"] = part1_predict(D, le)
    for k in RECIPES:
        preds[k] = recipe_predict(D, le, k, rec["recipes"][k]["increments"], curves, band_idx)
    log(f"=== every model predicted ({(time.time() - t0) / 60:.1f} min) ===")

    scored = (D["w5"] > 0) & np.isfinite(target) & np.logical_and.reduce([np.isfinite(p) for p in preds.values()])
    ev_code, ev_names = pd.factorize(le.event_id)
    ev_bin = ev_code[bin_le]
    E = len(ev_names)
    rng = np.random.default_rng(CB["seed"])
    mult = np.array([np.bincount(rng.integers(0, E, E), minlength=E) for _ in range(CB["boot_reps"])], float)
    L_bin = le.L_km.to_numpy()[bin_le]
    is_ratio = le.event_id.isin(ratio_events).to_numpy()[bin_le]
    sets = {"S-all": np.ones(len(bin_le), bool), "S-ratio": is_ratio}
    groups = {"all links": np.ones(len(bin_le), bool), "links < 1 km": L_bin < 1.0}
    lo_q, hi_q = CB["interval"]

    # ---- M1
    m1, m1p = [], []
    for sn, sm in sets.items():
        for gn, gm in groups.items():
            msk = scored & sm & gm
            bt = {}
            for mn in MODELS:
                r2 = (preds[mn] - target) ** 2
                S_ = np.bincount(ev_bin[msk], weights=r2[msk], minlength=E)
                N_ = np.bincount(ev_bin[msk], minlength=E).astype(float)
                bt[mn] = boot_rmse(S_, N_, mult)
                m1.append({"set": sn, "links": gn, "model": mn, "RMSE (dB)": float(np.sqrt(S_.sum() / N_.sum())), "values": int(N_.sum()),
                           "events": int((N_ > 0).sum()), "90 % (dB)": f"{np.percentile(bt[mn], lo_q):.4f} to {np.percentile(bt[mn], hi_q):.4f}"})
            r_o = m1[-len(MODELS)]["RMSE (dB)"]
            for k, mn in enumerate(MODELS[1:], 1):
                dlt = bt[mn] - bt["ours"]
                lo, hi = np.percentile(dlt, [lo_q, hi_q])
                verdict = f"{mn} better than ours" if hi < 0 else ("ours better than " + mn if lo > 0 else "tie")
                m1p.append({"set": sn, "links": gn, "model": mn, "RMSE(model) - RMSE(ours) (dB)": m1[-len(MODELS) + k]["RMSE (dB)"] - r_o,
                            "90 % lo": lo, "90 % hi": hi, "reading": verdict})
    log(f"=== M1 done ({(time.time() - t0) / 60:.1f} min) ===")

    # ---- M2 (a check)
    lags = CB["m2_lags_min"]
    use_all = np.ones(len(le))
    idx, _ = a7r.dry_values(D, fake6.window(D, np.asarray(D["rain"], float)), lags, CB["m2_rain_before_min"], s1.S1["dry_below_mmh"], use_all)
    m2 = []
    for sn, sm in sets.items():
        for i, lag in enumerate(lags):
            b = idx[:, i][idx[:, i] >= 0]
            b = b[sm[b] & scored[b]]
            eb = ev_bin[b]
            row = {"set": sn, "lag (min)": lag, "values": int(len(b)), "events": int(len(np.unique(eb))),
                   "measured median (dB)": float(np.median(target[b])) if len(b) else np.nan}
            for mn in MODELS:
                if not len(b):
                    row[mn] = "n/a"
                    continue
                e0 = np.median(preds[mn][b]) - np.median(target[b])
                draws = []
                for w in mult:
                    rw = np.repeat(b, w[eb].astype(int))
                    if len(rw):
                        draws.append(np.median(preds[mn][rw]) - np.median(target[rw]))
                lo, hi = np.percentile(draws, [lo_q, hi_q])
                row[mn] = f"{e0:+.3f} ({lo:+.3f} to {hi:+.3f})"
            m2.append(row)
    log(f"=== M2 done ({(time.time() - t0) / 60:.1f} min) ===")

    # ---- M3
    _, rule, rule_name, n_short, _ = sizes_and_rule()
    e4 = pd.read_csv(P.E4_CELLS)
    m3 = []
    for sn in sets:
        ev_keep = np.ones(len(le), bool) if sn == "S-all" else le.event_id.isin(ratio_events).to_numpy()
        Ds = sid.subset(D, ev_keep)
        les = le[ev_keep].reset_index(drop=True)
        bins_keep = ev_keep[bin_le]
        ps = {mn: np.where(scored[bins_keep], preds[mn][bins_keep], np.nan) for mn in MODELS}
        Ds["le_event"] = pd.factorize(les.event_id)[0]
        tab, ncell = m3_cells(Ds, les, ps, rule, n_short, e4[e4.eligible] if "eligible" in e4 else e4)
        m3.append(tab.assign(set=sn))
        log(f"   (M3) {sn}: {ncell} eligible cells")
    m3 = pd.concat(m3, ignore_index=True)
    log(f"=== M3 done ({(time.time() - t0) / 60:.1f} min) ===")
    return {"m1": pd.DataFrame(m1), "m1p": pd.DataFrame(m1p), "m2": pd.DataFrame(m2), "m3": m3, "dropped": n_drop, "drops": drops,
            "n_le": len(le), "events": E, "ratio_events": int(len(set(le.event_id) & ratio_events)), "minutes": (time.time() - t0) / 60,
            "rule": rule_name}


# =================================================================================================== report
def report(r):
    save(r["m1"].round(5), P.B_M1, index=False)
    save(r["m1p"].round(5), P.B_M1P, index=False)
    save(r["m2"], P.B_M2, index=False)
    save(r["m3"].round(5), P.B_M3, index=False)
    m3s = r["m3"].groupby(["set", "model"]).inside.agg(["sum", "size"]).reset_index().rename(columns={"sum": "cells inside", "size": "cells"})
    L = ["# Block B: the frozen equation against the existing models on VAL\n",
         "Written by `research/waa_ksurface/b_val.py` under `docs/waa_ksurface/block_b_prereg.md` (frozen in 37470fa; recipe fits frozen in the second "
         "record `a8_block_b_freeze2.json`). VAL read once; TEST not opened.\n",
         f"VAL link-events at 28 / 38 GHz scored: {r['n_le']} of {r['events']} events ({r['ratio_events']} with a gauge ratio); dropped because a model "
         f"cannot predict them (no stored a_link): {r['dropped']}. Run time {r['minutes']:.1f} min.\n",
         "## M1: RMSE of the 5 min A - c (paired event bootstrap, 1000 draws, seed 4100)\n", md_table(r["m1"].round(4)), "",
         "### Paired differences and the pre-registered readings (the only claims)\n", md_table(r["m1p"].round(4)), "",
         "## M2: after-rain error (a CHECK, not a measurement of the drying)\n",
         "median(model) - median(measured) at the bin ending 5 / 10 / 20 min after the path went dry (review check (5) windows), 90 % from the same "
         "event draws. Caveat (dry-shape control): the data-layer offsets and the selection of dry-path windows already distort a plain exponential "
         "drying (ratio 0.30 / 0.44 for a K0 truth), so these numbers are not a measurement of the drying dynamics. No claim.\n",
         md_table(r["m2"]), "",
         f"## M3: the length test (E4 candidate (ii), rule {r['rule']}) on every model's predictions\n",
         "Cells: eligible in E4 and on the VAL measurements. A model's cell is inside if both its intercept and slope lie in the measured 90 % "
         "event-bootstrap interval. No claim.\n", md_table(m3s), "", md_table(r["m3"].round(4)), ""]
    save("\n".join(L) + "\n", P.B_REPORT)
    log(f"report: {rel(P.B_REPORT)}")


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    if "--recipes" in sys.argv:
        r = recipes_part()
        (CACHE / "recipes.pkl").write_bytes(pickle.dumps(r))
        write_recipes(r)
        log(f"=== recipes and code checks: {'ALL OK' if r['ok'] else 'A CHECK FAILED: STOP'} ===")
        sys.exit(0 if r["ok"] else 1)
    if "--freeze" in sys.argv:
        freeze_recipes()
        return
    if "--freeze3" in sys.argv:
        freeze3()
        return
    if "--dryrun-train" in sys.argv:                            # Deviation 16: the whole --val path on TRAIN, to scratch; not a result
        out = CACHE / "dryrun"
        out.mkdir(parents=True, exist_ok=True)
        for k in ("B_M1", "B_M1P", "B_M2", "B_M3", "B_REPORT"):
            setattr(P, k, out / getattr(P, k).name)
        r = val_part("train")
        report(r)
        log(f"=== dry run on TRAIN ended ({r['minutes']:.1f} min); written to {out} ===")
        return
    if "--val" in sys.argv:
        f = CACHE / "val.pkl"
        if f.exists():
            log("=== VAL results exist: not read again (use --report) ===")
            sys.exit(2)
        r = val_part()
        if "stop" in r:
            sys.exit(1)
        f.write_bytes(pickle.dumps(r))
        report(r)
        log(f"=== Block B ended {time.strftime('%Y-%m-%d %H:%M:%S')} ===")
        return
    if "--report" in sys.argv:
        report(pickle.loads((CACHE / "val.pkl").read_bytes()))


if __name__ == "__main__":
    main()
