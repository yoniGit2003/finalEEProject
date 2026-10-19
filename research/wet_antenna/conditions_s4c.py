"""wet_antenna step S4.3: is the VAL gap the radar? (1) the gauge-rescaled length test, (2) per-event factors on heavy events.

S3b: every event-independent fit fails the VAL length test (VAL attenuates ~0.6-0.7 x per mm/h of radar rain). S4.1: VAL is
a rare event set, and on its large events the radar reads ~1.7 x the City gauges against ~0.9 x on TRAIN. Here:

  Part 1  the radar rain of every event is divided by its radar / gauge ratio g_e, taken from the City GAUGES only (the S4.1
          definitions: radar path mean of the gauge links, or the radar pixel at the gauges), and the S3b length test is
          redone through the FROZEN M4b-B (M3-B as a check; nothing refitted). The observed side is rescaled too: the
          steady-rain bins and their rain classes are recomputed from the rescaled R. Variants V0 (none), V1 (path mean),
          V2 (pixel), V3 (one pooled ratio per set).
  Part 2  a new fake check of the per-event factors on the HEAVY TRAIN events only (3 replicates); only if it passes, the
          S4.1 factors of the heavy TRAIN and VAL events against the conditions and against g_e (descriptive).

STANDING RULE: WAA is never estimated by subtracting a power-law rain term. The path part is separated only through its
dependence on L at matched rain. ITU-R P.838 only as a size reference, never to correct data. The rain correction comes
from the City gauges only, never from the links.
Every rule is in config `wet_antenna.conditions_s4c` (written before the run). Results kept in
research/.cache/wet_antenna/conditions_s4c/ so that the report can be re-made without recomputing.
Run from the repository root (about 10 min the first time):
    python -m research.wet_antenna.conditions_s4c
"""
import pickle
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cmlpipe.api import CMLData
from research.common.config import load_config
from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import length_test as LTL
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.event_factors import (AG, KG, M, EventParts, bh, estimate, fake_obs, frozen, radar_gauge_pixel,
                                               radar_gauge_primary, rescaled, spearman_boot)

CFG = load_config()
WA = CFG["wet_antenna"]
S2, S3, S4, C = WA["system_id_s2"], WA["system_id_s3"], WA["conditions_s4"], WA["conditions_s4c"]
LT = S3["length_test"]
MAIN, CHECK = S4["frozen"]["main"], S4["frozen"]["check"]
FITS = [MAIN, CHECK]
SETS = ("TRAIN", "VAL")
VARS = list(C["variants"])                                     # V0 .. V3
VSTYLE = {"V0": ("-", "x"), "V1": ("--", "o"), "V2": (":", "s"), "V3": ("-.", "^")}
MK = {"TRAIN": "o", "VAL": "s"}
RULES = C["rules"]
HV = C["heavy"]


# =================================================================================================== Part 1
def g_table(Ds, les, d, ev_tab):
    """g_e per event, both versions (S4.1 definitions), and the pooled ratio G per set."""
    rows = []
    for s in SETS:
        prim = radar_gauge_primary(Ds[s], les[s]).set_index("event_id")
        ids = list(pd.unique(les[s].event_id))
        pix = radar_gauge_pixel(d, ev_tab.loc[ids]).set_index("event_id")
        for e in ids:
            r = prim.loc[e] if e in prim.index else None
            q = pix.loc[e]
            use_p = r is not None and r.gauge_mm_per_link >= C["g_e"]["min_gauge_mm"]
            use_x = q.pixel_gauge_mm >= C["g_e"]["min_gauge_mm"]
            rows.append({"set": s, "event_id": e, "intensity_class": ev_tab.loc[e, "intensity_class"],
                         "gauge_links": int(r.gauge_links) if r is not None else 0, "gauge_mm_per_link": r.gauge_mm_per_link if r is not None else np.nan,
                         "radar_mm_sum": r.radar_mm_sum if r is not None else np.nan, "gauge_mm_sum": r.gauge_mm_sum if r is not None else np.nan,
                         "g_path_raw": r.radar_over_gauge if r is not None else np.nan, "path_used": bool(use_p),
                         "g_path": float(r.radar_over_gauge) if use_p else 1.0,
                         "pixel_gauge_mm": q.pixel_gauge_mm, "g_pixel_raw": q.pixel_radar_over_gauge, "pixel_used": bool(use_x),
                         "g_pixel": float(q.pixel_radar_over_gauge) if use_x else 1.0})
    g = pd.DataFrame(rows)
    G = {s: float(g[(g.set == s) & g.path_used].radar_mm_sum.sum() / g[(g.set == s) & g.path_used].gauge_mm_sum.sum()) for s in SETS}
    g["G_set"] = g.set.map(G)
    return g, G


def length_cells(t, n_ev, set_, variant):
    """Observed (event bootstrap) and every fit's own line on the same bins; one row per cell and fit."""
    obs = LTL.length_fits(t, n_ev, set_, LT)
    rows = []
    for f in FITS:
        mod = LTL.length_fits(t, n_ev, set_, LT, ycol=f"pred_{f}", boot=False).set_index(["band", "class"])
        for r in obs.dropna(subset=["I"]).itertuples():
            if (r.band, r._3) not in mod.index:
                continue
            q = mod.loc[(r.band, r._3)]
            if not np.isfinite(q.I):
                continue
            rows.append({"set": set_, "variant": variant, "fit": f, "band": r.band, "class": r._3, "R_cell": r.R_cell, "link_events": r.link_events,
                         "links_under_1km": r.links_under_1km, "I_obs": r.I, "I_lo": r.I_lo, "I_hi": r.I_hi, "I_model": q.I,
                         "s_obs": r.s, "s_lo": r.s_lo, "s_hi": r.s_hi, "s_model": q.s,
                         "eligible": bool(r.links_under_1km >= LT["verdict"]["min_short_links"]),
                         "I_inside": bool(r.I_lo <= q.I <= r.I_hi), "s_inside": bool(r.s_lo <= q.s <= r.s_hi)})
    return obs.assign(set=set_, variant=variant), pd.DataFrame(rows)


def verdicts(cells):
    lo, hi = RULES["ratio_range"]
    rows = []
    for (s, v, f), g in cells[cells.eligible].groupby(["set", "variant", "fit"], sort=False):
        rI, rs = float(np.median(g.I_obs / g.I_model)), float(np.median(g.s_obs / g.s_model))
        both = float((g.I_inside & g.s_inside).mean())
        ok = bool(lo <= rI <= hi and lo <= rs <= hi and both >= RULES["both_inside_share"])
        rows.append({"set": s, "variant": v, "fit": f, "eligible_cells": len(g), "I_ratio_median": rI, "s_ratio_median": rs,
                     "I_inside": int(g.I_inside.sum()), "s_inside": int(g.s_inside.sum()), "both_inside": int((g.I_inside & g.s_inside).sum()),
                     "both_share": both, "all_three_criteria": ok})
    return pd.DataFrame(rows)


# =================================================================================================== Part 2
def heavy_fake(D, le, f, ev_tab, noise):
    """The S4.1 fake design on the heavy TRAIN events, HV['replicates'] replicates."""
    ids = list(pd.unique(le.event_id))
    heavy = [e for e in ids if ev_tab.loc[e, "intensity_class"] == "heavy"]
    pairs = [(k, a) for k in S4["fake"]["k_values"] for a in S4["fake"]["alpha_values"]]
    rows = []
    for rep in range(HV["replicates"]):
        rng = np.random.default_rng([HV["seed"], rep])
        order = np.concatenate([rng.permutation(len(pairs)) for _ in range(int(np.ceil(len(heavy) / len(pairs))))])[:len(heavy)]
        rows.append([(e, *pairs[o]) for e, o in zip(heavy, order)])
    out = []
    for j, e in enumerate(heavy):
        idx = np.nonzero(D["le_event"] == ids.index(e))[0]
        ev = EventParts(D, le, f, idx)
        for rep in range(HV["replicates"]):
            _, kt, at = rows[rep][j]
            y5, w5, cf = fake_obs(ev, ev.truth_1min(kt, at), np.random.default_rng([HV["seed"], rep, j]), noise)
            fe = estimate(ev.grid(y5, w5, cf), S4["factors"]["ridge_per_link"])
            out.append({"event_id": e, "replicate": rep, "links": ev.nl, "k_true": kt, "alpha_true": at, "k_est": KG[fe["ik"]], "alpha_est": AG[fe["ia"]]})
        print(f"   heavy fake {j + 1}/{len(heavy)} {e}", flush=True)
        del ev
    fk = pd.DataFrame(out)
    fk["k_ok"] = (fk.k_est / fk.k_true - 1).abs() <= HV["tol"]
    fk["alpha_ok"] = (fk.alpha_est / fk.alpha_true - 1).abs() <= HV["tol"]
    fk["both_ok"] = fk.k_ok & fk.alpha_ok
    return fk


def heavy_tests(g):
    """If the heavy fake check passed: the S4.1 factors of the heavy events against the conditions and against g_e."""
    fac = pd.read_csv(P.COND4_FACTORS)
    fac = fac[fac.fit == MAIN]
    ev = pd.read_csv(P.COND4_EVENTS, parse_dates=["start", "end"])
    h = fac.merge(ev, on=["set", "event_id"], suffixes=("", "_ev")).merge(g[["set", "event_id", "g_path", "path_used", "g_path_raw"]], on=["set", "event_id"])
    h = h[h.intensity_class == "heavy"].sort_values("start").reset_index(drop=True)
    conds = {"temperature": "mean temperature (C)", "wind_speed": "mean wind speed (m/s)", "humidity": "mean humidity (%)",
             "dew_point_depression": "dew-point depression (C)", "dry_before_h": "dry hours before", "peak_over_mean": "peak / mean rain"}
    rng = np.random.default_rng(HV["boot_seed"])
    rows = []
    for col, lab in conds.items():
        rho, p, lo, hi, n = spearman_boot(np.log(h.alpha.to_numpy()), h[col].to_numpy(float), HV["boot_reps"], rng)
        rows.append({"test": f"ln alpha_e vs {lab}", "column": col, "rho": rho, "lo": lo, "hi": hi, "p": p, "events": n})
    t = pd.DataFrame(rows)
    t["bh_discovery"] = bh(t.p.fillna(1).to_numpy(), HV["bh_q"])
    u = h[h.path_used]
    rho, p, lo, hi, n = spearman_boot(np.log(u.k.to_numpy()), np.log(u.g_path.to_numpy()), HV["boot_reps"], rng)
    t = pd.concat([t, pd.DataFrame([{"test": "ln k_e vs ln g_e (path mean, events with g_e)", "column": "g_path", "rho": rho, "lo": lo, "hi": hi, "p": p,
                                     "events": n, "bh_discovery": np.nan}])], ignore_index=True)
    return h, t, conds


# =================================================================================================== compute
def compute():
    clock = time.time()
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    BANDS = list(sid.merged_bands(d, S2["data"])[0])
    ev_tab = pd.concat([d.events(split="train"), d.events(split="val")]).set_index("event_id")      # TEST events never read
    Ds = {"TRAIN": sid.load_data(P.SID3_CACHE / "data_train_all.npz", BANDS), "VAL": sid.load_data(P.SID3_CACHE / "data_val_all.npz", BANDS)}
    les = {"TRAIN": pickle.loads((P.SID3_CACHE / "le_train_all.pkl").read_bytes())[0].reset_index(drop=True),
           "VAL": pickle.loads((P.SID3_CACHE / "le_val_all.pkl").read_bytes())[0].reset_index(drop=True)}
    fz = {k: {**frozen(k, Ds["TRAIN"], les["TRAIN"]), "name": k} for k in FITS}
    g, G = g_table(Ds, les, d, ev_tab)
    print(f"g_e done: path-mean used for {int(g.path_used.sum())} events, pixel for {int(g.pixel_used.sum())}; G {G}", flush=True)

    scales = {"V0": lambda s: dict.fromkeys(g[g.set == s].event_id, 1.0),
              "V1": lambda s: dict(zip(g[g.set == s].event_id, g[g.set == s].g_path)),
              "V2": lambda s: dict(zip(g[g.set == s].event_id, g[g.set == s].g_pixel)),
              "V3": lambda s: dict.fromkeys(g[g.set == s].event_id, G[s])}
    obs_all, cells = [], []
    for v in VARS:
        for s in SETS:
            D2 = rescaled(Ds[s], les[s], scales[v](s))
            preds = {}
            for k, f in fz.items():
                Pm = sid.Problem2(D2, f["model"], M)
                x, _ = LTL.x_for(Pm, f, LTL.link_ids(D2, les[s]))
                preds[k] = Pm.predict(x, f["d"])
            t = LTL.steady_table(D2, les[s], LT, preds)
            o, c = length_cells(t, int(D2["le_event"].max()) + 1, s, v)
            obs_all.append(o)
            cells.append(c)
            print(f"{v} {s} done ({(time.time() - clock) / 60:.1f} min)", flush=True)
    obs_all, cells = pd.concat(obs_all, ignore_index=True), pd.concat(cells, ignore_index=True)

    fake = heavy_fake(Ds["TRAIN"], les["TRAIN"], fz[MAIN], ev_tab, sid.noise_table())
    print(f"heavy fake done ({(time.time() - clock) / 60:.1f} min)", flush=True)
    return {"g": g, "G": G, "obs": obs_all, "cells": cells, "fake": fake, "minutes": (time.time() - clock) / 60}


def model_free(obs):
    """Observed I, s of VAL against TRAIN per variant, on the cells eligible on VAL in that variant."""
    rows = []
    for v in VARS:
        o = obs[obs.variant == v].set_index(["set", "band", "class"])
        val = obs[(obs.variant == v) & (obs.set == "VAL") & (obs.links_under_1km >= LT["verdict"]["min_short_links"]) & obs.I.notna()]
        for r in val.itertuples():
            key = ("TRAIN", r.band, r._3)
            if key not in o.index or not np.isfinite(o.loc[key, "I"]):
                continue
            tr = o.loc[key]
            rows.append({"variant": v, "band": r.band, "class": r._3, "R_cell_VAL": r.R_cell, "R_cell_TRAIN": tr.R_cell, "I_VAL": r.I, "I_VAL_lo": r.I_lo,
                         "I_VAL_hi": r.I_hi, "I_TRAIN": tr.I, "I_TRAIN_lo": tr.I_lo, "I_TRAIN_hi": tr.I_hi, "s_VAL": r.s, "s_VAL_lo": r.s_lo, "s_VAL_hi": r.s_hi,
                         "s_TRAIN": tr.s, "s_TRAIN_lo": tr.s_lo, "s_TRAIN_hi": tr.s_hi, "I_VAL_over_TRAIN": r.I / tr.I, "s_VAL_over_TRAIN": r.s / tr.s,
                         "I_intervals_overlap": bool(r.I_lo <= tr.I_hi and tr.I_lo <= r.I_hi), "s_intervals_overlap": bool(r.s_lo <= tr.s_hi and tr.s_lo <= r.s_hi)})
    return pd.DataFrame(rows)


# =================================================================================================== figures
def fig_ratios(cells):
    c = cells[cells.eligible & (cells.fit == MAIN)].copy()
    c["I_ratio"], c["s_ratio"] = c.I_obs / c.I_model, c.s_obs / c.s_model
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.4))
    lo, hi = RULES["ratio_range"]
    for i, s in enumerate(SETS):
        for j, q in enumerate(("I", "s")):
            ax = axes[i][j]
            g = c[c.set == s]
            labels = list(dict.fromkeys(f"{b}\n{k}" for b, k in zip(g.band, g["class"])))
            for n, v in enumerate(VARS):
                h = g[g.variant == v]
                x = [labels.index(f"{b}\n{k}") + (n - 1.5) * 0.12 for b, k in zip(h.band, h["class"])]
                ls, mk = VSTYLE[v]
                ax.plot(x, h[f"{q}_ratio"], ls=ls, marker=mk, mfc="none", color=INK, ms=6, lw=1.0, label=f"{v}: median {h[f'{q}_ratio'].median():.2f}")
            for b in (lo, hi):
                ax.axhline(b, color=INK2, lw=0.8, ls=":")
            ax.axhline(1, color=INK2, lw=0.8)
            ax.set_xticks(range(len(labels))), ax.set_xticklabels(labels, fontsize=7)
            ax.set(ylabel=f"observed / {MAIN} {q}", title=f"{s}: {q} per eligible cell (dotted = {lo}-{hi})")
            ax.legend(frameon=False, fontsize=7.5)
            style_axes(ax)
    fig.suptitle("Step S4.3: the length test through the frozen M4b-B with the radar rescaled by the gauges. V0 none (solid, x), V1 path-mean g_e\n"
                 "(dashed, circles), V2 pixel g_e (dotted, squares), V3 one pooled ratio per set (dash-dot, triangles). Cells are recomputed per variant.", fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig


def fig_curves(obs, bands):
    fig, axes = plt.subplots(2, len(bands), figsize=(5.4 * len(bands), 8.2), squeeze=False)
    style = {("TRAIN", "V0"): ("-", "o"), ("VAL", "V0"): ("-", "s"), ("TRAIN", "V1"): ("--", "o"), ("VAL", "V1"): ("--", "s")}
    for j, band in enumerate(bands):
        for i, q in enumerate(("I", "s")):
            ax = axes[i][j]
            for (s, v), (ls, mk) in style.items():
                g = obs[(obs.set == s) & (obs.variant == v) & (obs.band == band)].dropna(subset=[q]).sort_values("R_cell")
                if not len(g):
                    continue
                y = g[q].to_numpy()
                ax.errorbar(g.R_cell * (1.03 if s == "VAL" else 1.0), y, yerr=[np.clip(y - g[f"{q}_lo"], 0, None), np.clip(g[f"{q}_hi"] - y, 0, None)], ls=ls, marker=mk,
                            mfc="none" if v == "V0" else INK2, color=INK if v == "V0" else INK2, ms=5, capsize=2, lw=1.1)
                ax.annotate(f"{s} {v}", (g.R_cell.iloc[-1], y[-1]), xytext=(4, 0), textcoords="offset points", fontsize=7, va="center")
            ax.set_xscale("log")
            ax.set(xlabel="rain R of the class (mm/h; V1 = R / g_e)", ylabel=f"{q} ({'dB' if q == 'I' else 'dB/km'})", title=f"{band}: observed {q}, model-free")
            style_axes(ax)
    fig.suptitle("Step S4.3: observed I(R) and s(R) of TRAIN (circles) and VAL (squares), before (V0: solid, hollow) and after the path-mean gauge rescaling\n"
                 "(V1: dashed, filled grey); 90 % event-bootstrap intervals; every curve is labelled at its end", fontsize=9.5)
    fig.tight_layout(rect=(0, 0, 0.97, 0.93))
    return fig


def fig_fake(fk):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.6))
    mk = {0.6: "v", 0.7: "v", 1.0: "o", 1.4: "^"}
    rj = np.random.default_rng(0)
    for ax, q in zip(axes, ("k", "alpha")):
        for t in sorted(fk[f"{q}_true"].unique()):
            g = fk[fk[f"{q}_true"] == t]
            ax.plot(t * np.exp(rj.uniform(-0.04, 0.04, len(g))), g[f"{q}_est"], mk[t], mfc="none", color=INK, ms=6, ls="none", label=f"true {q} = {t:g}")
            for b in (1 - HV["tol"], 1 + HV["tol"]):
                ax.plot([t * 0.9, t * 1.1], [t * b] * 2, "-", color=INK2, lw=0.8)
        ax.plot([0.3, 3], [0.3, 3], ":", color=INK2, lw=0.8)
        ax.set(xscale="log", yscale="log", xlabel=f"true {q}", ylabel=f"estimated {q}", title=f"{q} (short bars = +-{HV['tol'] * 100:.0f} %)")
        ax.legend(frameon=False, fontsize=8)
        style_axes(ax)
    fig.suptitle(f"Step S4.3: fake check on the heavy TRAIN events ({fk.event_id.nunique()} events x {HV['replicates']} replicates), true against estimated", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def fig_heavy(h, tests, conds):
    cols = list(conds) + ["g_path"]
    fig, axes = plt.subplots(2, 4, figsize=(14, 6.6))
    tt = tests.set_index("column")
    for ax, col in zip(axes.ravel(), cols):
        y = "k" if col == "g_path" else "alpha"
        for s in SETS:
            g = h[(h.set == s) & (h.path_used if col == "g_path" else True)]
            ax.plot(g[col], g[y], MK[s], mfc="none", color=INK, ms=6, ls="none", label=s)
        r = tt.loc[col]
        ax.set_yscale("log")
        if col == "g_path":
            ax.set_xscale("log")
        ax.set(xlabel=conds.get(col, "g_e (radar / gauge, path mean)"), ylabel=f"{y}_e", title=f"rho {r.rho:+.2f} [{r.lo:+.2f}, {r.hi:+.2f}]", )
        ax.legend(frameon=False, fontsize=7)
        style_axes(ax)
    for ax in axes.ravel()[len(cols):]:
        ax.set_visible(False)
    fig.suptitle("Step S4.3: heavy events, alpha_e against the conditions and k_e against g_e (S4.1 fits of M4b-B; TRAIN hollow circles, VAL hollow squares)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    clock = time.time()
    P.COND4C_CACHE.mkdir(parents=True, exist_ok=True)
    cache = P.COND4C_CACHE / "results.pkl"
    if cache.exists():
        R = pickle.loads(cache.read_bytes())
        print(f"results read from the cache ({rel(cache)}); delete that file to recompute", flush=True)
    else:
        R = compute()
        cache.write_bytes(pickle.dumps(R))
    cells = R["cells"]
    ver = verdicts(cells)
    mf = model_free(R["obs"])
    fk = R["fake"]
    share = float(fk.both_ok.mean())
    passed = share >= HV["pass_share"]
    save(R["g"].round(4), P.COND4C_G, index=False)
    save(cells.round(4), P.COND4C_CELLS, index=False)
    save(ver.round(4), P.COND4C_VERDICT, index=False)
    save(mf.round(4), P.COND4C_MODEL_FREE, index=False)
    save(fk.round(4), P.COND4C_FAKE, index=False)
    bands = [b for b in ("28 GHz", "38 GHz") if b in set(R["obs"].band)]
    figs = [(fig_ratios(cells), P.COND4C_FIG_RATIOS), (fig_curves(R["obs"], bands), P.COND4C_FIG_CURVES), (fig_fake(fk), P.COND4C_FIG_FAKE)]
    H = None
    if passed:
        h, tests, conds = heavy_tests(R["g"])
        save(h.drop(columns=["start", "end"]).round(4), P.COND4C_HEAVY, index=False)
        save(tests.round(4), P.COND4C_HEAVY_TESTS, index=False)
        figs.append((fig_heavy(h, tests, conds), P.COND4C_FIG_HEAVY))
        H = (h, tests, conds)
    for fig, path in figs:
        save(fig, path, dpi=100)
        plt.close(fig)
    report = write_report(R, ver, mf, share, passed, H, clock)
    save(report, P.COND4C_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {(time.time() - clock) / 60:.1f} min")


def readings(R, ver, mf, share, passed, H):
    """The READING paragraphs, written after looking at the numbers of the run (every number is taken from its tables)."""
    vt = ver.set_index(["set", "variant", "fit"])
    cells = R["cells"]
    elig = cells[cells.eligible & (cells.fit == MAIN)]
    n_val = {v: int(((elig.set == "VAL") & (elig.variant == v)).sum()) for v in C["variants"]}
    lost = sorted(set(elig[(elig.set == "VAL") & (elig.variant == "V0")]["class"]) - set(elig[(elig.set == "VAL") & (elig.variant == "V1")]["class"]))
    ms = mf.groupby("variant").agg(I=("I_VAL_over_TRAIN", "median"), s=("s_VAL_over_TRAIN", "median"))
    worst = mf[mf.variant == "V1"].sort_values("I_VAL_over_TRAIN").iloc[0]

    def v(set_, var, fit=MAIN):
        r = vt.loc[(set_, var, fit)]
        return f"I {r.I_ratio_median:.2f}, s {r.s_ratio_median:.2f}, both inside {int(r.both_inside)} / {int(r.eligible_cells)}"
    m3 = vt.loc[("VAL", "V1", CHECK)]
    q = {}
    q["verdicts"] = f"""**READING.** On VAL without rescaling (V0) the model over-predicts the length dependence: {v('VAL', 'V0')}.
Dividing the radar rain by the event's radar / gauge ratio closes the gap: V1 {v('VAL', 'V1')}; V2 (pixel) {v('VAL', 'V2')};
V3 (one number per set) {v('VAL', 'V3')}. TRAIN stays inside the rules but loses some sharpness: V0 {v('TRAIN', 'V0')}, V1 {v('TRAIN', 'V1')}.
So by the rules written before the run the VAL gap is **explained by the rain reference**, and TRAIN is not broken. Three caveats:
(a) the number of eligible VAL cells falls from {n_val['V0']} (V0) to {n_val['V1']} (V1): the class(es) {', '.join(lost) or '-'} disappear
because VAL rain is divided by about 1.7, and those were among the worst V0 cells, so part of the improvement is that the hardest
cells leave the test; (b) the check fit {CHECK} misses the share rule in V1 on VAL by rounding ({int(m3.both_inside)} / {int(m3.eligible_cells)} =
{m3.both_share:.5f} against the written {RULES['both_inside_share']}) and on TRAIN ({vt.loc[('TRAIN', 'V1', CHECK)].both_share:.2f}); its medians move the same way;
(c) the per-set V3 does as well as the per-event V1 on VAL, so the data cannot say whether the error is per event or one bias of the
VAL period."""
    q["curves"] = f"""**READING.** The model-free curves tell the same story without any model: the median VAL / TRAIN ratio of I goes from
{ms.loc['V0', 'I']:.2f} (V0) to {ms.loc['V1', 'I']:.2f} (V1), of s from {ms.loc['V0', 's']:.2f} to {ms.loc['V1', 's']:.2f}. The gap closes most of the way but
not all: the {worst.band} {worst['class']} mm/h cell still has VAL I at {worst.I_VAL_over_TRAIN:.2f} x TRAIN with non-overlapping intervals (V1), and V3 over-shoots s
(median {ms.loc['V3', 's']:.2f}). At light rain (0.5-2 mm/h) VAL and TRAIN agree after rescaling within their intervals."""
    if H is not None:
        h, tests, _ = H
        kt = tests[tests.column == "g_path"].iloc[0]
        ct = tests[tests.column != "g_path"].sort_values("p")
        top = ct.iloc[0]
        edge = h[(h.k <= h.k.min() + 1e-9) | (h.alpha_hi >= 2.99)]
        q["fake"] = f"""**READING.** On heavy events the factors are identifiable: {share * 100:.0f} % of the event-replicates recover both k and alpha
within +-{HV['tol'] * 100:.0f} % (rule {HV['pass_share'] * 100:.0f} %), as S4.1 already hinted (16 / 16). Part 2 therefore goes on."""
        q["heavy"] = f"""**READING.** (Descriptive only, {len(h)} events; not a law.) **k_e follows g_e**: Spearman {kt.rho:.2f}
(90 % [{kt.lo:.2f}, {kt.hi:.2f}], {int(kt.events)} events), above the expectation {HV['expect_k_g_spearman']} written before the run. The rain-scale factor that the
links ask for is the same scale the gauges see in the radar (links against radar, gauges against radar; the radar is shared, see
section 6): a second sign that the radar, not the antenna, makes events differ. **No condition is a BH discovery** for alpha_e (6 tests, q = {HV['bh_q']}); the strongest
is {top.test} (rho {top.rho:.2f}, p {top.p:.3f}), then dew-point depression and humidity (wetter air, larger alpha_e), all with
intervals that are wide for {len(h)} events. Events at the edge of the fit grid ({', '.join(edge.event_id)}: k at its lower edge or alpha
interval to the upper edge) weigh on these ranks and should not be read alone."""
    else:
        q["fake"] = "**READING.** The heavy-event factors are not identifiable either; Part 2 stops."
    ok_val = bool(vt.loc[("VAL", "V1", MAIN), "all_three_criteria"])
    ok_tr = bool(vt.loc[("TRAIN", "V1", MAIN), "all_three_criteria"])
    q["meaning"] = f"""- **WAA v3: keep {MAIN} unchanged.** With the rain reference corrected by the gauges, the frozen model passes the VAL length
  test ({'explained' if ok_val else 'not explained'}; TRAIN {'not broken' if ok_tr else 'broken'}), so the S3b gap does not call for a new antenna term or event-dependent
  WAA parameters. Nothing was refitted here; whether a refit on gauge-rescaled rain moves the parameters is a question for S5
  (it would change the science, so only as a written S5 sub-step).
- **S5 evaluation, rain reference:** report every error twice, on the raw radar (V0, as before) and on the radar divided by the
  event's path-mean g_e (V1, events with >= {C['g_e']['min_gauge_mm']} mm per gauge link; g_e = 1 otherwise), with V1 as the main reference and the
  pooled V3 as the sensitivity check. The gauge links alone (City gauge within 1 km) stay the independent check, since they do not use
  g_e from a different place. The per-event factors k_e are descriptive and are not used to correct anything.
- The weather question (alpha_e against conditions) stays open: nothing passes BH on {len(H[0]) if H is not None else 0} heavy events."""
    q["unclear"] = f"""- **Circularity.** The radar is already gauge-adjusted, and g_e uses the same kind of gauges, so V1 partly measures how well the
  adjustment worked in this period; it is not an independent test. g_e also compares point gauges with path means (or 2 x 2 km pixels).
- **Few events.** The VAL correction rests on {int((R['g'].path_used & (R['g'].set == 'VAL')).sum())} VAL events with a path-mean g_e; the large VAL events dominate the steady-rain bins.
- **Cells leave the test.** Rescaling lowers R, so the 5-10 mm/h VAL cells drop below 8 links; the improvement is partly fewer, easier cells.
- **Rule rounding.** "2/3" was written as {RULES['both_inside_share']}, so 4 / 6 counts as a miss ({CHECK} on VAL). Kept as written, not changed after the run.
- **The 2-5 mm/h cells** of VAL stay low against TRAIN even after rescaling (model-free); a radar error that depends on rain rate,
  or a rain type in VAL, would do that; not tested.
- **Part 2** factors come from the S4.1 fits (grid-bounded; some events at the grid edge); the k_e vs g_e agreement uses the same radar
  on both sides of the ratio, so part of it may be shared radar noise."""
    return q


def write_report(R, ver, mf, share, passed, H, clock):
    q = readings(R, ver, mf, share, passed, H)
    g, G, fk = R["g"], R["G"], R["fake"]
    vt = ver.set_index(["set", "variant", "fit"])

    def rule(set_):
        k = (set_, "V1", MAIN)
        return ("YES" if vt.loc[k, "all_three_criteria"] else "NO") if k in vt.index else "no eligible cell"
    standing = ("**Standing rule.** WAA is never estimated by subtracting a power-law rain term. The path part is separated only through its "
                "dependence on L at matched rain. ITU-R P.838 appears only as a size reference, never to correct data. In this step the rain "
                "correction comes from the City gauges only, never from the links.")
    used = g[g.path_used | g.pixel_used][["set", "event_id", "intensity_class", "gauge_links", "gauge_mm_per_link", "g_path_raw", "path_used", "pixel_gauge_mm",
                                           "g_pixel_raw", "pixel_used"]].round(3)
    not_used = g[~(g.path_used | g.pixel_used)]
    ver_show = ver.round(3)
    cells_show = R["cells"][R["cells"].eligible & (R["cells"].fit == MAIN)][["set", "variant", "band", "class", "R_cell", "I_obs", "I_lo", "I_hi", "I_model", "s_obs",
                                                                              "s_lo", "s_hi", "s_model", "I_inside", "s_inside"]].round(3)
    mf_sum = mf.groupby("variant").agg(cells=("band", "size"), I_VAL_over_TRAIN_median=("I_VAL_over_TRAIN", "median"),
                                      s_VAL_over_TRAIN_median=("s_VAL_over_TRAIN", "median"), I_overlap=("I_intervals_overlap", "sum"),
                                      s_overlap=("s_intervals_overlap", "sum")).reset_index().round(3)
    fk_by = fk.groupby(["k_true", "alpha_true"]).agg(event_replicates=("event_id", "size"), recovered=("both_ok", "sum")).reset_index()
    fk_rep = fk.groupby("replicate").agg(event_replicates=("event_id", "size"), recovered=("both_ok", "sum"), k_ok=("k_ok", "sum"), alpha_ok=("alpha_ok", "sum")).reset_index()
    part2 = ""
    if H is not None:
        h, tests, _ = H
        part2 = f"""Real factors of the heavy events ({len(h)} events: {int((h.set == 'TRAIN').sum())} TRAIN, {int((h.set == 'VAL').sum())} VAL; the S4.1 fits of {MAIN},
`{rel(P.COND4_FACTORS)}`, no rerun; `{rel(P.COND4C_HEAVY)}`):

{md_table(h[['set', 'event_id', 'links', 'k', 'k_lo', 'k_hi', 'alpha', 'alpha_lo', 'alpha_hi', 'temperature', 'wind_speed', 'humidity', 'dew_point_depression', 'dry_before_h', 'peak_over_mean', 'g_path_raw', 'path_used']].round(3))}

Spearman tests (90 % event-bootstrap intervals; BH at q = {HV['bh_q']} over the {int(tests.bh_discovery.notna().sum())} condition tests; the k_e vs g_e test has its own
expectation, >= {HV['expect_k_g_spearman']}; `{rel(P.COND4C_HEAVY_TESTS)}`):

{md_table(tests.round(3))}

![heavy events](../figures/{P.COND4C_FIG_HEAVY.name})

{q.get('heavy', '')}"""
    else:
        part2 = "**Part 2 stops here (fake check FAIL), as the rule says.**"

    return f"""# wet_antenna step S4.3: is the VAL gap the radar? Gauge-rescaled length test, and per-event factors on heavy events

Written by `research/wet_antenna/conditions_s4c.py` (libraries `event_factors.py`, `system_id.py`, `length_test.py`). Settings and every
rule: `wet_antenna.conditions_s4c` in `research/config.yaml`, written on 2026-10-15 BEFORE the run. Frozen model: **{MAIN}** of S3b (WAA v3,
provisional), {CHECK} as a check; NOTHING is refitted. Plan: `docs/ROADMAP.md`, S4.3.

{standing}

## 1. What was run

- **g_e** = radar rain / City-gauge rain per event, the two S4.1 definitions (code moved unchanged to `event_factors.py`; the S4.1 tables
  are reproduced byte for byte after the move): **path mean** = the radar path mean of the gauge links (City gauge within 1 km) against
  their gauges, used if >= {C['g_e']['min_gauge_mm']} mm per gauge link; **pixel** = the radar pixel at the 10 City gauges against the gauges, used if the
  gauges' SUM is >= {C['g_e']['min_gauge_mm']} mm (the S4.1 code rule). Other events: g_e = 1. Pooled per set (path mean): G_TRAIN = {G['TRAIN']:.3f},
  G_VAL = {G['VAL']:.3f}. (`{rel(P.COND4C_G)}`.)
- Path-mean g_e used for {int(g.path_used.sum())} events ({int((g.path_used & (g.set == 'TRAIN')).sum())} TRAIN, {int((g.path_used & (g.set == 'VAL')).sum())} VAL); pixel g_e for {int(g.pixel_used.sum())}
  ({int((g.pixel_used & (g.set == 'TRAIN')).sum())} TRAIN, {int((g.pixel_used & (g.set == 'VAL')).sum())} VAL). The events that got a g_e:

{md_table(used)}

  g_e = 1 (too little gauge rain) for the other {len(not_used)} events: {', '.join(f"{r.event_id} ({r.set})" for r in not_used.itertuples())}.
- Variants (R -> R / g in the model input AND in the observed bins; the steady-rain bins and classes are recomputed per variant):
  {'; '.join(f'**{k}** {v}' for k, v in C['variants'].items())}.
- The length test exactly as S3b check 1 (Huber line; observed with the event bootstrap, 1000 draws; eligible cell = >= 8 links < 1 km).
- Part 2: the S4.1 fake design on the heavy TRAIN events, {HV['replicates']} replicates.
- **The circularity, plainly:** the radar product is gauge-adjusted (`pipeline/docs/methods.md`, section 9: "it is gauge-adjusted, so it is
  not independent of the gauges"), so dividing it by a gauge ratio partly undoes or redoes an adjustment that already used gauges; and
  g_e compares point gauges with path means (or a 2 x 2 km pixel), so it holds sampling error as well as radar error. Only the City
  gauges are used here, never the links.

## 2. Part 1: the verdicts

Rules (written before the run): **explained** on VAL = in V1 ({MAIN}) the median observed / model ratio of I AND of s within
{RULES['ratio_range'][0]}-{RULES['ratio_range'][1]} AND both inside in >= 2/3 of the eligible cells; **TRAIN not broken** = the same three on TRAIN in V1.

- **Explained on VAL: {rule('VAL')}.**
- **TRAIN not broken: {rule('TRAIN')}.**

Per set, variant and fit (`{rel(P.COND4C_VERDICT)}`):

{md_table(ver_show)}

{q.get('verdicts', '')}

Every eligible cell, {MAIN} (`{rel(P.COND4C_CELLS)}`, both fits):

{md_table(cells_show)}

![ratios](../figures/{P.COND4C_FIG_RATIOS.name})

## 3. The observed curves before and after

Model-free S3.1 I and s of VAL against TRAIN per variant, on the cells eligible on VAL (`{rel(P.COND4C_MODEL_FREE)}`):

{md_table(mf_sum)}

{md_table(mf.round(3))}

![observed curves](../figures/{P.COND4C_FIG_CURVES.name})

{q.get('curves', '')}

## 4. Part 2: heavy events only

Fake check: the S4.1 design (frozen {MAIN}; k in {S4['fake']['k_values']}, alpha in {S4['fake']['alpha_values']}; S2 noise; ridge {S4['factors']['ridge_per_link']} per link) on
the {fk.event_id.nunique()} heavy TRAIN events x {HV['replicates']} replicates. Rule: PASS if >= {HV['pass_share'] * 100:.0f} % of the event-replicates recover both within
+-{HV['tol'] * 100:.0f} %. **Result: {'PASS' if passed else 'FAIL'}** ({int(fk.both_ok.sum())} of {len(fk)} = {share * 100:.0f} %; k alone {int(fk.k_ok.sum())}, alpha alone {int(fk.alpha_ok.sum())}).
(`{rel(P.COND4C_FAKE)}`.)

{md_table(fk_rep)}

{md_table(fk_by)}

![heavy fake](../figures/{P.COND4C_FIG_FAKE.name})

{q.get('fake', '')}

{part2}

## 5. What this means for WAA v3 and for the S5 evaluation

{q.get('meaning', '')}

## 6. Unclear / possibly wrong

{q.get('unclear', '')}

Run: computation {R['minutes']:.1f} min; this report {(time.time() - clock) / 60:.1f} min (budget {C['time_budget_h']} h).
"""


if __name__ == "__main__":
    main()
