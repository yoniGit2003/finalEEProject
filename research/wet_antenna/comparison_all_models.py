"""wet_antenna: the physical consistency of every WAA model in the comparison, on TEST (a POST-FREEZE SUMMARY, nothing refitted).

TEST was read in S8 / S8b / S8c. This script runs the frozen fits of those steps forward once more and puts, side by side, what
each model says about the physics against what the data say without a model:

  1. the length test through each model (exactly the S8 P1 procedure, `final_model.through`): per eligible cell the model's own
     slope s (path loss per km) and intercept I (the WAA level at L = 0) against the model-free interval of the same cell;
  2. the dry-path error at 5, 10 and 20 min after the path became dry (measured pure WAA minus the model's W; the S7 / S8 windows);
  3. the wetting and drying times each model uses, as 1/e times, next to the gauge timing of our equation;
  4. the number of fitted parameters and how many of them are identified (fake recovery, interval not at a bound), or that the
     model has no such check;
  5. (added after the first run, dated in the config) the Spearman correlation of each model's link factor a_link with the link
     length over the TRAIN links: a WAA that does not depend on L should give a_link independent of L.

Models: B0 no WAA, B1 constant WAA while it rains, B2 WAA without memory (S8), Schleiss a / b and Pastorek a / b (S8c), ours (the
S8b closed-F model = the equation file) and the S8 table-F candidate (reference). Our equation is unchanged.

Settings: config `wet_antenna.comparison_all_models`. Run from the repository root:
    python -m research.wet_antenna.comparison_all_models
"""
import hashlib
import json
import pickle
import sys
import time

import numpy as np
import pandas as pd

from cmlpipe.api import CMLData
from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import length_test as LTL
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.event_factors import rescaled, spearman_boot
from research.wet_antenna.final_model import CFG, DATA, LT, WA, predict, s6_models
from research.wet_antenna.literature_s8c import LIT, Set
from research.wet_antenna.test_s8b import b_models
from research.wet_antenna.validation import LB, boot_items, ci, dry_windows, need, radar_max, w_minutes

S = WA["comparison_all_models"]
S8C = WA["literature_s8c"]
ORDER = ["B0", "B1", "B2", "Schleiss a", "Schleiss b", "Pastorek a", "Pastorek b", "ours", "table F"]
LABEL = {"B0": "no WAA (B0)", "B1": "constant WAA while it rains (B1)", "B2": "WAA without memory (B2)", "Schleiss a": "Schleiss (a)",
         "Schleiss b": "Schleiss (b)", "Pastorek a": "Pastorek (a)", "Pastorek b": "Pastorek (b)", "ours": "ours (final equation)",
         "table F": "ours, table F (S8 candidate)"}


def w_minutes_any(D, le, fz):
    """W alone at 1 min of a frozen model of our family (path term 0, c = 0); `validation.w_minutes` for M4b, the same for M1."""
    if fz["model"] == "M4b":
        return w_minutes(D, le, fz)
    D0 = dict(D)
    D0["le_c"] = np.zeros(len(le))
    Pm = sid.Problem2(D0, fz["model"], fz["M"], agg=1)
    x, _ = LTL.x_for(Pm, fz, LTL.link_ids(D, le))
    x[Pm.names == Pm.pname] = 0.0
    return Pm.predict(x, fz["d"])


def efold(tau_nominal):
    """The 1/e time (min) of the discrete Schleiss step W <- W + (target - W) 3 dt / tau with dt = 1 min."""
    return float(-1.0 / np.log(1.0 - 3.0 / tau_nominal))


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    clock = time.time()
    log = lambda s: print(s, flush=True)                                                    # noqa: E731
    fit_file = P.S6_CACHE / "plug_shared_FP_F_and_P.pkl"
    sha = hashlib.sha256(fit_file.read_bytes()).hexdigest()
    if sha != S8C["ours"]["fit_sha256"]:
        log(f"STOP: the fit behind our equation changed (sha256 {sha}).")
        return
    eq = json.loads(P.S9_EQUATION.read_text(encoding="utf-8"))
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    BANDS = list(sid.merged_bands(d, DATA)[0])
    let = pickle.loads((P.SID3_CACHE / "le_train_all.pkl").read_bytes())[0].reset_index(drop=True)
    Dt1 = sid.load_data(P.S5A_CACHE / "data_train_v1.npz", BANDS)
    mods, _ = s6_models(Dt1, let, BANDS)
    ours = {"ours": mods["closed F"]["fz"], "table F": mods["candidate"]["fz"], **b_models(Dt1, let, BANDS)}
    lit = {}
    for n in LIT:
        lit[n] = pickle.loads((P.S8C_CACHE / f"fit_{n.replace(' ', '_')}.pkl").read_bytes())
    nl_train = int(Dt1["nl"])
    Lk = let.groupby("cml_id").L_km.first()
    rng = np.random.default_rng(S["a_link_vs_L"]["seed"])
    rho = {}
    for m in ORDER:
        a_ = (ours[m]["a"] if m in ours else lit[m]["a"]) if m not in ("B0", "Schleiss a", "Pastorek a") else None
        if a_ is None:
            rho[m] = (np.nan, np.nan, np.nan, 0)
            continue
        ids = [k for k in a_ if k in Lk.index]
        r_, _, lo, hi, n = spearman_boot(Lk.loc[ids].to_numpy(float), np.array([a_[k] for k in ids], float), S["a_link_vs_L"]["boot_reps"], rng)
        rho[m] = (r_, lo, hi, n)
    del Dt1
    log(f"frozen fits loaded ({(time.time() - clock) / 60:.1f} min)")

    # ------------------------------------------------------------------------------------------- TEST (TRAIN links), V1 as S8b / S8c
    le_tr = pickle.loads((P.S8_CACHE / "le_test_train.pkl").read_bytes())[0].reset_index(drop=True)
    gt = pd.read_csv(P.S8_G_E)
    gm = {e: 1.0 for e in set(le_tr.event_id)}
    gm.update({r.event_id: (r.g_path if r.path_used else 1.0) for r in gt.itertuples()})
    D = rescaled(sid.load_data(P.S8_CACHE / "data_test_train.npz", BANDS), le_tr, gm)
    pred = {m: predict(D, le_tr, f, f["M"])[0] for m, f in ours.items()}
    Ts = Set(D, le_tr)
    W1 = {}
    for n, f in lit.items():
        pred[n], W1[n] = Ts.predict(f)
    del Ts
    log(f"TEST predicted ({(time.time() - clock) / 60:.1f} min)")

    # ------------------------------------------------------------------------------------------- 1. the length test through every model (S8 P1)
    t = LTL.steady_table(D, le_tr, LT, {m: pred[m] for m in ORDER})
    n_ev = int(D["le_event"].max()) + 1
    obs = LTL.length_fits(t, n_ev, "TEST", LT)
    cell_rows = []
    for m in ORDER:
        mod = LTL.length_fits(t, n_ev, "TEST", LT, ycol=f"pred_{m}", boot=False).set_index(["band", "class"])
        for r in obs.dropna(subset=["I"]).to_dict("records"):
            key = (r["band"], r["class"])
            if key not in mod.index or not np.isfinite(mod.loc[key, "I"]):
                continue
            q = mod.loc[key]
            cell_rows.append({"model": m, "band": r["band"], "class": r["class"], "R_cell": r["R_cell"], "links_under_1km": r["links_under_1km"],
                              "eligible": bool(r["links_under_1km"] >= LT["verdict"]["min_short_links"]),
                              "I_obs": r["I"], "I_lo": r["I_lo"], "I_hi": r["I_hi"], "I_model": q.I,
                              "s_obs": r["s"], "s_lo": r["s_lo"], "s_hi": r["s_hi"], "s_model": q.s,
                              "I_inside": bool(r["I_lo"] <= q.I <= r["I_hi"]), "s_inside": bool(r["s_lo"] <= q.s <= r["s_hi"])})
    cells = pd.DataFrame(cell_rows)
    e8 = pd.read_csv(P.S8_CELLS)
    e8 = e8[e8.eligible]
    mine = cells[(cells.model == "table F") & cells.eligible]
    repro = bool(len(mine) == len(e8) and np.allclose(np.sort(mine.s_model), np.sort(e8.s_model), atol=1e-3)
                 and np.allclose(np.sort(mine.I_model), np.sort(e8.I_model), atol=1e-3))   # the S8 table is rounded to 4 decimals
    log(f"length test done: the S8 P1 cells of table F reproduced: {repro} ({(time.time() - clock) / 60:.1f} min)")

    # ------------------------------------------------------------------------------------------- 2. the dry path on TEST (S7 / S8 windows)
    W1t = {m: w_minutes_any(D, le_tr, ours[m]) for m in ("ours", "table F", "B0", "B1", "B2")}
    W1t.update(W1)
    Rm, T5, col = radar_max(d, {"TEST": (D, le_tr)}, blank_test=False)
    win, curves = dry_windows({"TEST": (D, le_tr)}, {"TEST": gm}, Rm, T5, col, {"TEST": W1t})
    ob = np.array([c["obs"] for c in curves])
    nmin = win.dry_minutes.to_numpy()
    draws = boot_items(pd.factorize(win.event_id)[0], S["boot_reps"], S["seed"])
    arr = {m: np.array([c[m] for c in curves]) for m in ORDER}
    dry_rows = []
    for tm in S["dry_t_min"]:
        sel = (nmin > tm) & np.isfinite(ob[:, tm])
        for m in ORDER:
            err = ob[:, tm] - arr[m][:, tm]
            lo, hi = ci([np.median(err[i[sel[i]]]) for i in draws if sel[i].any()])
            dry_rows.append({"model": m, "t_min": tm, "windows": int(sel.sum()), "measured_W": float(np.median(ob[sel, tm])),
                             "model_W": float(np.median(arr[m][sel, tm])), "err": float(np.median(err[sel])), "err_lo": lo, "err_hi": hi})
    dry = pd.DataFrame(dry_rows)
    cv_rows = []                                                                            # median curves over all windows (figure)
    for tm in range(60):
        ok = nmin > tm
        if ok.sum() < 3:
            continue
        cv_rows.append({"t_min": tm, "windows": int(ok.sum()), "measured": float(np.nanmedian(ob[ok, tm])),
                        **{m: float(np.nanmedian(arr[m][ok, tm])) for m in ORDER}})
    save(pd.DataFrame(cv_rows).round(4), P.CMP_DRY_CURVES, index=False)
    s8c = pd.read_csv(P.S8C_DRY)
    s8c = s8c[s8c.length == "all"].set_index("t_min")
    dry_repro = bool(all(abs(dry[(dry.model == "ours") & (dry.t_min == tm)].err.iloc[0] - s8c.loc[tm, "err ours"]) < 1e-3 for tm in S["dry_t_min"]))
    log(f"dry path done: {len(win)} windows; S8c reproduced for ours: {dry_repro} ({(time.time() - clock) / 60:.1f} min)")

    # ------------------------------------------------------------------------------------------- 3. timing, 4. parameters
    dy = eq["dynamics_fixed"]
    timing = {
        "B0": ("no antenna term", np.nan, np.nan),
        "B1": ("none: a step with the radar wet flag (V1 > 0.1 mm/h)", 0.0, 0.0),
        "B2": ("none: follows the rain input", 0.0, 0.0),
        "Pastorek a": ("none: follows the rain input", 0.0, 0.0),
        "Pastorek b": ("none: follows the rain input", 0.0, 0.0),
        "ours": (f"fast store (share {dy['phi']:.2f}) plus a slow store (tau_s {dy['tau_s']:.0f} min, not identified)", dy["tau_f_wet"], dy["tau_f_dry"]),
        "table F": ("as ours", dy["tau_f_wet"], dy["tau_f_dry"]),
    }
    for n in ("Schleiss a", "Schleiss b"):
        tw, td = lit[n]["shape"]
        timing[n] = (f"nominal tau_wet {tw:.1f} min, tau_dry {td:.1f} min in the step 3 dt / tau; 1/e times given", efold(tw), efold(td))
    gauge = S["gauge_timing_min"]

    bnd = {"amp": S8C["fit"]["amp_bounds_db"], "c": S8C["fit"]["c_band_bounds"], "loga": S8C["fit"]["log_a_bounds"],
           "tau": S8C["schleiss"]["tau_bounds_min"], "zeta": S8C["pastorek"]["zeta_bounds"]}

    def at_bound(v, b, rel_tol=1e-3):
        return bool(min(abs(v - b[0]), abs(v - b[1])) <= rel_tol * max(1.0, abs(b[1] - b[0])))

    par = {}
    for n, f in lit.items():
        b_ = n.endswith("b")
        hits = []
        if b_:
            hits += [f"amplitude {v:.2f}" for v in f["amp"] if at_bound(v, bnd["amp"])]
            hits += [f"shape {v:.3f}" for v in f["shape"] if at_bound(v, bnd["tau"] if "Schleiss" in n else bnd["zeta"])]
        hits += [f"c_band {v:.3f}" for v in f["c_band"] if at_bound(v, bnd["c"])]
        n_a = sum(at_bound(np.log(v), bnd["loga"]) for v in f["a"].values()) if b_ else 0
        k = 3 + (3 + len(f["shape"]) if b_ else 0)
        par[n] = {"global_fitted": k, "a_link": "fitted per link (ridge)" if b_ else "1 (fixed)",
                  "identified": "no check: one TRAIN fit, no fake-data recovery, no intervals",
                  "at_bound": (", ".join(hits) if hits else "none") + (f"; {n_a} a_link at a bound" if b_ else "")}
    for b in ("B0", "B1", "B2"):
        z = pickle.loads((P.S8_CACHE / f"fit_{b}.pkl").read_bytes())["z"]
        k = len(z) - nl_train
        par[b] = {"global_fitted": k, "a_link": "in the vector, no effect (no antenna term)" if b == "B0" else "fitted per link (ridge)",
                  "identified": "no check: one TRAIN fit, no fake-data recovery, no intervals", "at_bound": "not checked"}
    par["ours"] = {"global_fitted": 7, "a_link": "fitted per link (ridge), p10-p90 0.77-1.24",
                   "identified": f"{eq['checks']['fake_data']}; 90 % intervals from 40 event-bootstrap refits, none at a bound",
                   "at_bound": "fixed: tau_wet, tau_dry, phi (gauge fit), tau_s at its search bound (not identified)"}
    par["table F"] = {"global_fitted": "free curves (S6 table)", "a_link": "as ours", "identified": "S6 fake recovery (62 numbers within 20 %)",
                      "at_bound": "as ours"}

    # ------------------------------------------------------------------------------------------- the one table
    rows = []
    for m in ORDER:
        e = cells[(cells.model == m) & cells.eligible]
        r = {"model": m, "label": LABEL[m], "cells": len(e), "s_inside": int(e.s_inside.sum()), "I_inside": int(e.I_inside.sum()),
             "both_inside": int((e.s_inside & e.I_inside).sum()), "P1_pass": bool(len(e) and (e.s_inside & e.I_inside).sum() >= need(len(e))),
             "s_ratio_model_over_obs_median": float(np.median(e.s_model / e.s_obs)), "I_ratio_model_over_obs_median": float(np.median(e.I_model / e.I_obs))}
        for tm in S["dry_t_min"]:
            q = dry[(dry.model == m) & (dry.t_min == tm)].iloc[0]
            r[f"dry_err_{tm}min_db"], r[f"dry_err_{tm}min_lo"], r[f"dry_err_{tm}min_hi"] = q.err, q.err_lo, q.err_hi
        tx, tw, td = timing[m]
        r.update({"a_link_vs_L_rho": rho[m][0], "a_link_vs_L_lo": rho[m][1], "a_link_vs_L_hi": rho[m][2], "a_link_links": rho[m][3]})
        r.update({"wet_time_min": tw, "dry_time_min": td, "timing": tx, "gauge_wet_train_min": gauge["train"][0], "gauge_dry_train_min": gauge["train"][1],
                  "gauge_wet_val_min": gauge["val"][0], "gauge_dry_val_min": gauge["val"][1], **{f"params_{k}": v for k, v in par[m].items()}})
        rows.append(r)
    tab = pd.DataFrame(rows)
    measured = dry[dry.model == "ours"].set_index("t_min").measured_W
    save(tab.round(4), P.CMP_TABLE, index=False)
    save(cells.round(4), P.CMP_CELLS, index=False)
    save(dry.round(4), P.CMP_DRY, index=False)
    save(write_report(tab, cells, dry, measured, repro, dry_repro, len(win), sha), P.CMP_REPORT)
    log(f"done in {(time.time() - clock) / 60:.1f} min")


def write_report(tab, cells, dry, measured, repro, dry_repro, nwin, sha):
    L = []
    a = L.append
    a("# wet_antenna: the physical consistency of every WAA model in the comparison (TEST)\n")
    a("Written by `research/wet_antenna/comparison_all_models.py`. Settings: `wet_antenna.comparison_all_models` in `research/config.yaml`.\n")
    a("**Status.** A post-freeze summary. TEST was read in S8 / S8b / S8c. The frozen fits of those steps are run forward again; nothing "
      f"is refitted and our equation is unchanged (fit sha256 `{sha}`, checked).\n")
    a(f"- The S8 P1 cells of the table-F candidate are reproduced exactly: **{repro}**. The S8c dry-path errors of our equation are reproduced: "
      f"**{dry_repro}** ({nwin} dry-path windows).")
    a("- Length test: exactly the S8 P1 procedure. The model-free intercept I and slope s of each cell (Huber line of A - c against L, event "
      "bootstrap 90 %) against the same line fitted through the model's own predictions on the same bins. Eligible cell: >= 8 links under 1 km.")
    a("- Dry path: median over the windows where the radar shows no rain on any pixel of the path (S7 / S8 rules) of measured W minus model W; "
      "positive = the model holds too little water. 90 % event bootstrap.\n")
    a("## The table\n")
    show = tab[["label", "cells", "s_inside", "I_inside", "both_inside", "s_ratio_model_over_obs_median", "I_ratio_model_over_obs_median",
                "dry_err_5min_db", "dry_err_10min_db", "dry_err_20min_db", "a_link_vs_L_rho", "a_link_vs_L_lo", "a_link_vs_L_hi",
                "wet_time_min", "dry_time_min", "params_global_fitted"]].copy()
    a(md_table(show.round(3)) + "\n")
    a(f"Measured dry-path WAA (median): {', '.join(f'{v:.2f} dB at {t} min' for t, v in measured.items())}. Gauge timing of our equation: "
      f"wetting {tab.gauge_wet_train_min.iloc[0]} min, drying {tab.gauge_dry_train_min.iloc[0]} min on TRAIN (VAL {tab.gauge_wet_val_min.iloc[0]} and "
      f"{tab.gauge_dry_val_min.iloc[0]} min).\n")
    a("Timing and parameters:\n")
    a(md_table(tab[["label", "timing", "params_global_fitted", "params_a_link", "params_identified", "params_at_bound"]]) + "\n")
    a("Dry-path errors with intervals:\n")
    dd = dry.copy()
    dd["err (90 %)"] = [f"{r.err:+.3f} [{r.err_lo:+.3f}, {r.err_hi:+.3f}]" for r in dd.itertuples()]
    a(md_table(dd.pivot_table(index="model", columns="t_min", values="err (90 %)", aggfunc="first").loc[ORDER].reset_index()) + "\n")
    a("Length-test cells (eligible), the model's slope and intercept against the model-free interval:\n")
    e = cells[cells.eligible].copy()
    e["s"] = [f"{r.s_model:.3f} ({'in' if r.s_inside else 'out'})" for r in e.itertuples()]
    e["I"] = [f"{r.I_model:.2f} ({'in' if r.I_inside else 'out'})" for r in e.itertuples()]
    obs = e[e.model == "ours"][["band", "class", "R_cell", "s_obs", "s_lo", "s_hi", "I_obs", "I_lo", "I_hi"]].round(3)
    a(md_table(obs) + "\n")
    a(md_table(e.pivot_table(index=["band", "class"], columns="model", values="s", aggfunc="first")[ORDER].reset_index()) + "\n")
    a(md_table(e.pivot_table(index=["band", "class"], columns="model", values="I", aggfunc="first")[ORDER].reset_index()) + "\n")
    a("## Reading\n")
    o = tab.set_index("model")
    a(f"- Path loss: our model's slope lies inside the measured interval in {o.loc['ours', 's_inside']} of {o.loc['ours', 'cells']} cells "
      f"(median model / measured {o.loc['ours', 's_ratio_model_over_obs_median']:.2f}); Schleiss (b) {o.loc['Schleiss b', 's_inside']} "
      f"({o.loc['Schleiss b', 's_ratio_model_over_obs_median']:.2f}), Pastorek (b) {o.loc['Pastorek b', 's_inside']} "
      f"({o.loc['Pastorek b', 's_ratio_model_over_obs_median']:.2f}), no WAA {o.loc['B0', 's_inside']} ({o.loc['B0', 's_ratio_model_over_obs_median']:.2f}).")
    a(f"- WAA level at L = 0: ours inside in {o.loc['ours', 'I_inside']} of {o.loc['ours', 'cells']} (the light-rain cells of S8 P1 fail); "
      f"Schleiss (b) {o.loc['Schleiss b', 'I_inside']}, Pastorek (b) {o.loc['Pastorek b', 'I_inside']}, no WAA {o.loc['B0', 'I_inside']}.")
    a(f"- After rain (10 min): ours {o.loc['ours', 'dry_err_10min_db']:+.2f} dB, Schleiss (b) {o.loc['Schleiss b', 'dry_err_10min_db']:+.2f}, "
      f"Pastorek (b) {o.loc['Pastorek b', 'dry_err_10min_db']:+.2f}, no memory (B2) {o.loc['B2', 'dry_err_10min_db']:+.2f}.")
    a("- a_link against L (Spearman, 90 % link bootstrap): " + "; ".join(f"{LABEL[m]} {o.loc[m, 'a_link_vs_L_rho']:+.2f} "
      f"[{o.loc[m, 'a_link_vs_L_lo']:+.2f}, {o.loc[m, 'a_link_vs_L_hi']:+.2f}]" for m in ORDER if np.isfinite(o.loc[m, 'a_link_vs_L_rho'])) + ".")
    a("- These are numbers, not verdicts: no pass rule was written for this summary. P1 (both inside in >= 6 of 8 cells) is shown for each "
      "model only as the S8 rule would read it.\n")
    a("## Unclear / possibly wrong\n")
    a("- Post-freeze: the literature fits were made after TEST had been read (S8c); nothing here was tuned on TEST.")
    a("- B0-B2 and the literature models have no fake-data recovery and no intervals; their parameter counts are the fitted global numbers "
      "(a_link per link counted apart).")
    a("- The Schleiss times are the 1/e times of its discrete step (W += (target - W) 3 dt / tau, dt = 1 min), not its nominal tau.")
    return "\n".join(L)


if __name__ == "__main__":
    main()
