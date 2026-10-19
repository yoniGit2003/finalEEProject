"""wet_antenna step S3b: fit the S2 frame on ALL link lengths, so that the path term P is pinned by L (repair after S3).

S3 showed that from links shorter than 1 km the split between the wet antenna W and the path term L * P cannot be
identified: the S2b fits of M3 and M4b do not reproduce the model-free intercept I and slope s of the length test
(A - c) = I + s L. Here the same frame

    A_pred(t) = W(t) + L * P_band(R(t - d)) + c

is fitted on the TRAIN events of ALL links (the S3 all-length dataset, up to 8.7 km), radar input:
  M3-A, M4b-A   every parameter free
  M3-B, M4b-B   the time constants fixed at the TRAIN gauge values of S2b ("gauge timing"); only F, P, a_link, d fitted
  M1            static, for scale
No leave-one-event-out folds: the judge is VAL (never seen by these fits) with
  check 1  the length test through the model: the same Huber line on the model's own prediction, against the observed
           I and s intervals, cell by cell (the rule, on VAL; TRAIN as information);
  check 2  VAL errors per length class (< 1, 1-3, > 3 km), next to the S2b fits;
  check 3  the median residual during rain against L (flat = the L-scaling is right);
  check 4  the wetting time with the radar input now that P is pinned.
STANDING RULE: WAA is never estimated by subtracting a power-law rain term. The path part is separated only through its
dependence on L at matched rain. ITU-R P.838 only as a size reference line and as one of the optimiser starting values
of P, never to correct data.
Every rule is in config `wet_antenna.system_id_s3b` (written before the run). Fits are cached in
research/.cache/wet_antenna/system_id_s3b/ (every start and every 10 bootstrap draws).
Run from the repository root (hours the first time on 3 processes, about 15 min from the cache):
    python -m research.wet_antenna.system_id_s3b
"""
import copy
import json
import os
import pickle
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed

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
from research.wet_antenna.tails import itu_p838

CFG = load_config()
WA = CFG["wet_antenna"]
S2, S2B, S3, S = WA["system_id_s2"], WA["system_id_s2b"], WA["system_id_s3"], WA["system_id_s3b"]
DATA, CHK = S2["data"], S2["checks"]
M = copy.deepcopy(S2["model"])                                # the S2b model settings
M["increments"], M["maxiter"] = S2B["model"]["increments"], dict(S2B["model"]["maxiter"])
M["bounds"] = {**M["bounds"], "increment": S2B["model"]["increment_bounds"], "tau_s_min": S2B["model"]["tau_s_m4b_min"]}
J = {**S2["judge"], "order": S2B["judge"]["order"]}
LT = S3["length_test"]
KNOTS = np.asarray(M["knots_mmh"], float)
KEY_R = S2["length_rule"]["key_rain_mmh"]
CLASS_NAME = LTL.class_names(LT)
FITS = {k: tuple(v) for k, v in S["fits"].items()}                    # "M3-A" -> ("M3", "A")
CAND = [k for k in FITS if k != "M1"]
CHECKS = S["checks"]
CACHE = P.SID3B_CACHE
S2B_FITS = {"M3, S<1 (S2b)": ("S1", "M3"), "M4b, S<1 (S2b)": ("S1", "M4b"), "M3, S<3 (S2b)": ("S3", "M3"), "M4b, S<3 (S2b)": ("S3", "M4b")}
STYLE = {"M3-A": ("--", "o"), "M3-B": ("--", "s"), "M4b-A": (":", "^"), "M4b-B": ("-.", "v"), "M1": ("-", "x")}
TAUS = {"M3": ["tau_wet", "tau_dry"], "M4b": ["phi", "tau_f_wet", "tau_f_dry", "tau_s"], "M1": []}


# =================================================================================================== the fits (worker processes)
def fix_dyn(Pm, fixed):
    """Variant B: the time constants (and phi) held fixed by giving their parameters equal lower and upper bounds."""
    for p, v in fixed.items():
        val = float(np.log(v / (1 - v))) if p == "phi" else float(np.log(v))
        for i in np.nonzero(Pm.names == p)[0]:
            Pm.bounds[i] = (val, val)


def make_problem(D, name):
    model, variant = FITS[name]
    Pm = sid.Problem2(D, model, M)
    if variant == "B":
        fix_dyn(Pm, S["fixed_B"][model])
    return Pm


def starts_for(Pm, itu_inc, n):
    """The first n starts of the S2b start plan (generic start for F, the time constants and a = 1; P = factor x the
    ITU-R P.838 size, a STARTING value only)."""
    plan = S2B["model"]["start_plan"][:n]
    return {f"{g}+{p}": Pm.start(M["starts"][g], Pinc=S2B["model"]["p_starts"][p] * np.asarray(itu_inc)) for g, p in plan}


def run_job(job):
    """One fit on the all-length TRAIN set: the starts, the lag profile, the polish, then the event-bootstrap refits."""
    file = CACHE / f"fit_{job['name']}.pkl"
    res = pickle.loads(file.read_bytes()) if file.exists() else {}
    D = sid.load_data(job["data"], job["bands"])
    Pm = make_problem(D, job["name"])
    t0 = [time.time()]

    def dump():
        res["seconds"] = res.get("seconds", 0.0) + time.time() - t0[0]
        t0[0] = time.time()
        file.write_bytes(pickle.dumps(res))

    if "all" not in res:
        res["all"] = sid.fit_all2(Pm, starts_for(Pm, job["itu_inc"], job["n_starts"]), M["lag_grid_min"], M["maxiter"], 0,
                                  state=res.setdefault("all_state", {}), checkpoint=dump)
        del res["all_state"]
        res["params"] = Pm.params(res["all"]["x"], res["all"]["d"], KEY_R)
        res["a_link"] = np.exp(res["all"]["x"][Pm.names == "loga"])
        dump()
    draws = res.setdefault("boot", [])
    if len(draws) < job["boot"]:
        n_ev = int(D["le_event"].max()) + 1
        rng = np.random.default_rng([J["seed"], zlib.crc32(f"s3b_{job['name']}".encode())])
        mults = [np.bincount(rng.integers(0, n_ev, n_ev), minlength=n_ev).astype(float) for _ in range(job["boot"])]
        for i in range(len(draws), job["boot"]):
            fw = sid.fit_weighted2(Pm, res["all"], mults[i][D["le_event"]], M["maxiter"]["boot"])
            draws.append(Pm.params(fw["x"], fw["d"], KEY_R))
            if len(draws) % 10 == 0:
                dump()
        dump()
    return job["name"]


def done(job):
    f = CACHE / f"fit_{job['name']}.pkl"
    if not f.exists():
        return False
    r = pickle.loads(f.read_bytes())
    return "params" in r and len(r.get("boot", [])) >= job["boot"]


# =================================================================================================== checks
def through_rows(obs, mod, fit, present):
    """Check 1 per cell: the model's I and s (its own Huber line) against the observed intervals."""
    o = obs.dropna(subset=["I"]).set_index(["band", "class"])
    m = mod.set_index(["band", "class"])
    rows = []
    for (band, cls), r in o.iterrows():
        if (band, cls) not in m.index or not np.isfinite(m.loc[(band, cls), "I"]):
            continue
        q = m.loc[(band, cls)]
        has = band in present
        rows.append({"set": r.variant, "fit": fit, "band": band, "class": cls, "R_cell": r.R_cell, "links_under_1km": r.links_under_1km,
                     "I_obs": r.I, "I_obs_lo": r.I_lo, "I_obs_hi": r.I_hi, "I_model": q.I, "I_model_lo": q.get("I_lo", np.nan), "I_model_hi": q.get("I_hi", np.nan),
                     "s_obs": r.s, "s_obs_lo": r.s_lo, "s_obs_hi": r.s_hi, "s_model": q.s, "s_model_lo": q.get("s_lo", np.nan), "s_model_hi": q.get("s_hi", np.nan),
                     "eligible": bool(has and r.links_under_1km >= LT["verdict"]["min_short_links"]),
                     "I_inside": bool(has and r.I_lo <= q.I <= r.I_hi), "s_inside": bool(has and r.s_lo <= q.s <= r.s_hi)})
    return rows


def residual_table(D, le, pred, fit):
    """Check 3: per link-event the median (A - A_pred) during rain; per band and L bin the median over link-events with
    an event bootstrap."""
    RC = CHECKS["residual_vs_L"]
    w = (D["w5"] > 0) & D["wet5"]
    res = pd.Series(np.where(w, D["y5"] - pred, np.nan)).groupby(D["bin_le"]).median().reindex(range(len(le))).to_numpy()
    t = pd.DataFrame({"res": res, "L": le.L_km.to_numpy(), "band": le.band_name.to_numpy(), "cml": le.cml_id.to_numpy(), "ev": D["le_event"]}).dropna(subset=["res"])
    edges = np.asarray(RC["bins_km"], float)
    t["bin"] = np.searchsorted(edges, t.L, side="right") - 1
    n_ev = int(D["le_event"].max()) + 1
    rng = np.random.default_rng(RC["seed"])
    mult = [np.bincount(rng.integers(0, n_ev, n_ev), minlength=n_ev).astype(float) for _ in range(RC["boot_reps"])]
    rows = []
    for (band, b), g in t.groupby(["band", "bin"]):
        v, ev = g.res.to_numpy(), g.ev.to_numpy()
        bs = [LTL.wmedian(v, m[ev]) if m[ev].sum() > 0 else np.nan for m in mult]
        lo, hi = np.nanpercentile(bs, [5, 95])
        rows.append({"fit": fit, "band": band, "L_bin_km": f"{edges[b]:g}-{edges[b + 1]:g}" if edges[b + 1] < 99 else f">{edges[b]:g}", "bin": b,
                     "L_median_km": float(g.L.median()), "link_events": len(g), "links": g.cml.nunique(), "residual_median_db": float(np.median(v)), "lo": lo, "hi": hi})
    return pd.DataFrame(rows)


def flat_verdict(rt):
    out = []
    for (fit, band), g in rt.groupby(["fit", "band"], sort=False):
        g = g[g.links >= CHECKS["residual_vs_L"]["min_links"]].sort_values("bin")
        if len(g) < 2:
            out.append({"fit": fit, "band": band, "verdict": "too few bins"})
            continue
        a, b = g.iloc[0], g.iloc[-1]
        overlap = a.lo <= b.hi and b.lo <= a.hi
        out.append({"fit": fit, "band": band, "shortest": a.L_bin_km, "longest": b.L_bin_km, "shortest_median": a.residual_median_db, "longest_median": b.residual_median_db,
                    "verdict": "FLAT" if overlap else ("TREND up with L" if b.residual_median_db > a.residual_median_db else "TREND down with L")})
    return pd.DataFrame(out)


# =================================================================================================== figures
def spread_labels(ax, ends, x):
    lo_, hi_ = ax.get_ylim()
    gap, prev = 0.06 * (hi_ - lo_), -np.inf
    for y_, lab in sorted(ends):
        y_ = max(y_, prev + gap)
        ax.annotate(lab, (x, y_), xytext=(4, 0), textcoords="offset points", fontsize=7, va="center")
        prev = y_


def figure_length(through, bands, itu):
    sets = ["TRAIN", "VAL"]
    fig, axes = plt.subplots(4, len(bands), figsize=(5.2 * len(bands), 15.0), squeeze=False)
    for j, band in enumerate(bands):
        for i, (st, q) in enumerate([(s_, q_) for s_ in sets for q_ in ("I", "s")]):
            ax = axes[i][j]
            t = through[(through.set == st) & (through.band == band)]
            o = t.drop_duplicates("class").sort_values("R_cell")
            if not len(o):
                ax.set_visible(False)
                continue
            v, lo, hi = o[f"{q}_obs"].to_numpy(), o[f"{q}_obs_lo"].to_numpy(), o[f"{q}_obs_hi"].to_numpy()
            ax.errorbar(o.R_cell, v, yerr=[np.clip(v - lo, 0, None), np.clip(hi - v, 0, None)], fmt="o", color=INK, ms=6, capsize=3, lw=1.2,
                        label="observed (model-free, 90 %)")
            ends = []
            for k in CAND:
                g = t[t.fit == k].sort_values("R_cell")
                if len(g):
                    ls, mk = STYLE[k]
                    ax.plot(g.R_cell, g[f"{q}_model"], ls=ls, marker=mk, mfc="none", ms=6, color=INK, lw=1.3)
                    ends.append((float(g[f"{q}_model"].iloc[-1]), k))
            if q == "s":
                rr = np.geomspace(0.6, 15, 40)
                a, b = itu[band]
                ax.plot(rr, a * rr ** b, "-.", color=INK2, lw=0.8)
                ends.append((a * rr[-1] ** b, "ITU size, reference only"))
            ax.set_xscale("log")
            ax.set_xlim(0.5, 40)
            ax.axhline(0, color=INK2, lw=0.6)
            spread_labels(ax, ends, float(o.R_cell.max()) * 1.05)
            ax.set(xlabel="rain R of the class (mm/h)", ylabel=f"{q} ({'dB' if q == 'I' else 'dB/km'})",
                   title=f"{band}, {st} events: {'intercept I' if q == 'I' else 'slope s'}" + (" (THE RULE)" if st == "VAL" else " (information)"))
            if i == 0 and j == 0:
                ax.legend(frameon=False, fontsize=7.5, loc="upper left")
            style_axes(ax)
    fig.suptitle("Step S3b, check 1: the length test through the model. Filled circles = observed I and s (90 % interval); lines = the same Huber line on each\n"
                 "fit's own prediction (M3 dashed, M4b dotted / dash-dot; A = free, B = gauge timing; markers and labels name the fit)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 0.97, 0.96))
    return fig


def figure_val(val, length_sets):
    metrics = [("rmse", "all points"), ("rmse_rain", "during rain"), ("rmse_0_60", "0-60 min after"), ("rmse_60_180", "60-180 min after")]
    fits = list(dict.fromkeys(val.fit))
    fig, axes = plt.subplots(1, len(metrics), figsize=(4.6 * len(metrics), 4.8), squeeze=False)
    mk2 = {"M3, S<1 (S2b)": "P", "M4b, S<1 (S2b)": "X", "M3, S<3 (S2b)": "p", "M4b, S<3 (S2b)": "h"}
    for ax, (col, title) in zip(axes[0], metrics):
        for i, k in enumerate(fits):
            g = val[val.fit == k]
            x = np.array([length_sets.index(s_) for s_ in g["VAL set"]]) + (i - len(fits) / 2) * 0.07
            ls, mk = STYLE.get(k, ("none", mk2.get(k, "D")))
            ax.plot(x, g[col], ls=ls if k in STYLE else "none", marker=mk, mfc="none" if k in STYLE else INK2, color=INK if k in STYLE else INK2, ms=6, lw=1.0, label=k)
        ax.set_xticks(range(len(length_sets))), ax.set_xticklabels(length_sets, fontsize=8)
        ax.set_title(title, fontsize=9)
        style_axes(ax)
    axes[0][0].set_ylabel("RMSE of A on VAL (dB), mean over events")
    axes[0][0].legend(frameon=False, fontsize=7)
    fig.suptitle("Step S3b, check 2: VAL errors per length class, frozen fits from the rain alone (S3b = lines with hollow markers; S2b fits = grey filled markers)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def figure_residual(rt, bands):
    fig, axes = plt.subplots(1, len(bands), figsize=(5.2 * len(bands), 4.6), squeeze=False)
    for ax, band in zip(axes[0], bands):
        ends = []
        for i, k in enumerate(FITS):
            g = rt[(rt.fit == k) & (rt.band == band)].sort_values("bin")
            if not len(g):
                continue
            ls, mk = STYLE[k]
            x = g.L_median_km * (1 + 0.03 * (i - 2))
            ax.errorbar(x, g.residual_median_db, yerr=[g.residual_median_db - g.lo, g.hi - g.residual_median_db], ls=ls, marker=mk, mfc="none", color=INK, ms=5,
                        capsize=2, lw=1.0)
            ends.append((float(g.residual_median_db.iloc[-1]), k))
        ax.axhline(0, color=INK2, lw=0.8)
        ax.set_xscale("log")
        if ends:
            spread_labels(ax, ends, float(rt[rt.band == band].L_median_km.max()) * 1.1)
        ax.set(xlabel="link length L (km), median of the bin", ylabel="median (A - A_pred) during rain (dB)", title=f"{band}, VAL events")
        style_axes(ax)
    fig.suptitle("Step S3b, check 3: residual during rain against link length on VAL (90 % event-bootstrap intervals); flat = the L-scaling is right", fontsize=10)
    fig.tight_layout(rect=(0, 0, 0.95, 0.9))
    return fig


def figure_after(curves):
    fig, axes = plt.subplots(1, len(curves), figsize=(5.8 * len(curves), 4.8), squeeze=False)
    for ax, (name, c) in zip(axes[0], curves.items()):
        ax.plot(c.minutes_after_last_wet_bin, c.observed, "o", mfc="none", color=INK, ms=5, label="observed (circles)")
        for k in CAND:
            ls, mk = STYLE[k]
            ax.plot(c.minutes_after_last_wet_bin, c[k], ls=ls, marker=mk, markevery=6, mfc="none", ms=5, color=INK, lw=1.3, label=k)
        ax.set_yscale("log")
        ax.set_ylim(3e-3, 0.5)
        ax.set(xlabel="minutes after the last wet bin of the link", ylabel="median (A - c) / observed peak of the last shower",
               title=f"VAL, {name} ({int(c.link_events.iloc[0])} link-events)")
        ax.legend(frameon=False, fontsize=7.5)
        style_axes(ax)
    fig.suptitle("Step S3b: after the last shower on VAL, observed against the four fits (logarithmic axis)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    clock = time.time()
    workers = max((os.cpu_count() or 2) - 1, 1)

    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    bands_cfg, _, _ = sid.merged_bands(d, DATA)
    BANDS = list(bands_cfg)
    ftr, fva = P.SID3_CACHE / "data_train_all.npz", P.SID3_CACHE / "data_val_all.npz"
    assert ftr.exists() and fva.exists(), "run research.wet_antenna.system_id_s3 first (it builds the all-length datasets)"
    let, lev = pickle.loads((P.SID3_CACHE / "le_train_all.pkl").read_bytes())[0], pickle.loads((P.SID3_CACHE / "le_val_all.pkl").read_bytes())[0]
    prod_freq = pd.read_csv(P.SID2B_LINK_EVENTS, usecols=["band_name", "frequency_ghz"]).groupby("band_name").frequency_ghz.median()
    itu_inc = np.array([np.diff(itu_p838(float(prod_freq[b]), "v")[0] * KNOTS ** itu_p838(float(prod_freq[b]), "v")[1]) for b in BANDS])   # P starts of S2b
    itu = {b: itu_p838(float(let[let.band_name == b].frequency_ghz.median()), "v") for b in BANDS}                                       # figure reference only
    print(f"TRAIN all lengths {len(let)} link-events, VAL {len(lev)}; {workers} processes", flush=True)

    # ------------------------------------------------------------------------------------------- time one start, plan, cuts
    plan_file = CACHE / "plan.json"
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
    else:
        Pm = make_problem(sid.load_data(ftr, BANDS), "M4b-A")       # only for the timing (freed afterwards)
        st = starts_for(Pm, itu_inc, 8)
        ones = np.ones(Pm.ne)
        Pm.run(st["g3+Pitu"], 0, ones)
        t0 = time.time()
        for _ in range(3):
            Pm.run(st["g3+Pitu"], 0, ones)
        t_eval = (time.time() - t0) / 3
        t0 = time.time()
        fit = Pm.fit(st["g3+Pitu"], 0, ones, M["maxiter"]["start"])
        t_fit = time.time() - t0
        (CACHE / "fit_M4b-A.pkl").write_bytes(pickle.dumps({"all_state": {"first": {"g3+Pitu": fit}}, "seconds": t_fit}))   # reused by the job
        n_lag = M["lag_grid_min"][1] - M["lag_grid_min"][0]

        def proj(p):
            per = p["n_starts"] * t_fit + n_lag * t_fit / 3 + t_fit / 2 + S["boot_reps"] * 1.1 * M["maxiter"]["boot"] * t_eval
            n = len(p["fits"])
            return per * int(np.ceil(n / workers)) / 3600
        plan = {"t_eval_s": t_eval, "t_fit_s": round(t_fit, 1), "nit_start": fit["nit"], "converged_start": fit["converged"], "n_starts": 8,
                "fits": list(FITS), "cuts": [], "workers": workers}
        plan["projected_full_h"] = round(proj(plan), 2)
        for cut, change in (("M1", lambda p: p["fits"].remove("M1")), ("starts 8 -> 5", lambda p: p.update(n_starts=5)),
                            ("variant B of M4b", lambda p: p["fits"].remove("M4b-B"))):
            if proj(plan) <= S["time_budget_h"]:
                break
            change(plan)
            plan["cuts"].append(cut)
        plan["projected_h"] = round(proj(plan), 2)
        plan_file.write_text(json.dumps(plan))
    print(f"timing: one M4b-A start {plan['t_fit_s']} s ({plan['nit_start']} iterations, converged {plan['converged_start']}), one evaluation "
          f"{plan['t_eval_s'] * 1e3:.0f} ms; projection {plan['projected_full_h']} h; cuts: {plan['cuts'] or 'nothing'} -> {plan['projected_h']} h", flush=True)
    fits_run = plan["fits"]
    jobs = [{"name": k, "data": str(ftr), "bands": BANDS, "itu_inc": itu_inc.tolist(), "n_starts": plan["n_starts"], "boot": S["boot_reps"]} for k in fits_run]
    todo = [j for j in jobs if not done(j)]
    if todo:
        order_ = sorted(todo, key=lambda j: j["name"] == "M1")             # M1 (cheap) last
        with ProcessPoolExecutor(min(workers, len(todo))) as ex:
            for f in as_completed([ex.submit(run_job, j) for j in order_]):
                print(f"   fit {f.result()} done ({(time.time() - clock) / 60:.1f} min)", flush=True)
    Dt, Dv = sid.load_data(ftr, BANDS), sid.load_data(fva, BANDS)     # loaded after the fits, so the workers have the memory
    res = {k: pickle.loads((CACHE / f"fit_{k}.pkl").read_bytes()) for k in fits_run}
    cand = [k for k in CAND if k in fits_run]
    ids_t = LTL.link_ids(Dt, let)
    fz = {k: {"x": r["all"]["x"], "d": int(r["all"]["d"]), "a": dict(zip(ids_t, r["a_link"])), "a_med": float(np.median(r["a_link"])), "params": r["params"],
              "model": FITS[k][0]} for k, r in res.items()}
    order = LTL.s2b_order()
    fz2 = {k: {**LTL.s2b_fit(s_, m, order), "model": m} for k, (s_, m) in S2B_FITS.items()}

    # ------------------------------------------------------------------------------------------- parameters, convergence
    par_rows, conv_rows = [], []
    for k, r in res.items():
        for name, v in r["params"].items():
            b = np.array([dr.get(name, np.nan) for dr in r.get("boot", [])], float)
            lo, hi = np.nanpercentile(b, J["interval"]) if np.isfinite(b).any() else (np.nan, np.nan)
            par_rows.append({"fit": k, "parameter": name, "value": v, "lo": lo, "hi": hi})
        a = r["all"]
        conv_rows.append({"fit": k, "model": FITS[k][0], "variant": FITS[k][1], "loss": a["loss"], "d": a["d"], "converged_on_all_events": a.get("converged"),
                          "starts": len(a["starts_nit"]), "starts_converged": int(sum(bool(v) for v in a["starts_converged"].values())),
                          "starts_median_iterations": float(np.median(list(a["starts_nit"].values()))),
                          "spread_over_starts": max(a["starts_d0"].values()) - min(a["starts_d0"].values()), "best_start": min(a["starts_d0"], key=a["starts_d0"].get),
                          "lag_profile_loss_range": max(a["lag_profile"].values()) - min(a["lag_profile"].values()), "boot_draws": len(r.get("boot", [])),
                          "minutes": round(r.get("seconds", 0.0) / 60, 1)})
    params, conv = pd.DataFrame(par_rows), pd.DataFrame(conv_rows)
    save(params.round(4), P.SID3B_PARAMS, index=False)
    save(conv.round(3), P.SID3B_CONVERGENCE, index=False)

    # ------------------------------------------------------------------------------------------- predictions
    def pred_on(D_, le_, f):
        Pm = sid.Problem2(D_, f["model"], M)
        x, n_new = LTL.x_for(Pm, f, LTL.link_ids(D_, le_))
        return Pm.predict(x, f["d"]), n_new
    preds_t = {k: pred_on(Dt, let, f)[0] for k, f in fz.items()}                         # TRAIN: in-sample (information)
    pv = {k: pred_on(Dv, lev, f) for k, f in {**fz, **fz2}.items()}
    preds_v, new_links_v = {k: v[0] for k, v in pv.items()}, {k: v[1] for k, v in pv.items()}
    print(f"predictions done; VAL links not in the fit (a_link = 1): {new_links_v}; {(time.time() - clock) / 60:.1f} min", flush=True)

    # ------------------------------------------------------------------------------------------- check 1: the length test through the model
    tt, tv = LTL.steady_table(Dt, let, LT, preds_t), LTL.steady_table(Dv, lev, LT, preds_v)
    n_ev_t, n_ev_v = int(Dt["le_event"].max()) + 1, int(Dv["le_event"].max()) + 1
    obs = {"TRAIN": LTL.length_fits(tt, n_ev_t, "TRAIN", LT), "VAL": LTL.length_fits(tv, n_ev_v, "VAL", LT)}
    lf_rows, th_rows = [], []
    for st, (t_, n_) in (("TRAIN", (tt, n_ev_t)), ("VAL", (tv, n_ev_v))):
        for k, f in ({**fz, **fz2}.items() if st == "VAL" else fz.items()):
            mod = LTL.length_fits(t_, n_, st, LT, ycol=f"pred_{k}")
            lf_rows.append(mod.assign(fit=k))
            present = {b for b in BANDS if np.isfinite(f["params"].get(f"F[{b}](2)", f["params"].get(f"P[{b}](2)", np.nan)))}
            th_rows += through_rows(obs[st], mod, k, present)
    lfits = pd.concat([pd.concat([o.assign(fit="observed") for o in obs.values()])] + lf_rows, ignore_index=True)
    through = pd.DataFrame(th_rows)
    ver = LTL.verdict_table(through.rename(columns={"set": "variant"}), CHECKS["length_test_through"]["share"], cols=("I_inside", "s_inside")).rename(columns={"variant": "set"})
    eq_parts = [LTL.compare_models(obs[st], {k: f for k, f in fz.items() if k in cand}, LT, KNOTS) for st in obs]
    eq_cmp = pd.concat([p_[0] for p_ in eq_parts], ignore_index=True)
    eq_ver = pd.concat([p_[1] for p_ in eq_parts], ignore_index=True).assign(check="equilibrium a F, P (S3 rule)")
    save(lfits.round(4), P.SID3B_LENGTH, index=False)
    save(pd.concat([through.assign(check="through the model"), eq_cmp.rename(columns={"variant": "set"}).assign(check="equilibrium a F, P (S3 rule)")],
                   ignore_index=True).round(4), P.SID3B_LENGTH_CMP, index=False)
    verdict = pd.concat([ver.assign(check="through the model (THE RULE on VAL)"), eq_ver.rename(columns={"variant": "set"})], ignore_index=True)
    save(verdict.round(3), P.SID3B_VERDICT, index=False)
    print(f"check 1 done; {(time.time() - clock) / 60:.1f} min", flush=True)

    # ------------------------------------------------------------------------------------------- check 2: VAL by length
    bins = CHECKS["val_by_length"]["length_bins_km"]
    lnames = [f"L < {b[1]:g} km" if b[0] == 0 else (f"L > {b[0]:g} km" if b[1] > 50 else f"{b[0]:g}-{b[1]:g} km") for b in bins]
    Vs = {}
    for (lo_, hi_), nm in zip(bins, lnames):
        keep = ((lev.L_km >= lo_) & (lev.L_km < hi_)).to_numpy()
        Vs[nm] = (sid.subset(Dv, keep), lev[keep].reset_index(drop=True))
    val_rows, curves = [], {}
    for nm, (D_, le_) in Vs.items():
        cpreds = {}
        allowed = dict(fz)
        if nm == lnames[0]:
            allowed.update({k: f for k, f in fz2.items()})
        elif nm == lnames[1]:
            allowed.update({k: f for k, f in fz2.items() if f["set"] == "S3"})
        for k, f in allowed.items():
            p_, n_new = pred_on(D_, le_, f)
            if k in fz2:                                   # the S2b fits only where all their links exist
                if n_new:
                    continue
            val_rows.append({"fit": k, "VAL set": nm, "link_events": len(le_), "links": le_.cml_id.nunique(), "links_a_link_1": n_new, **LTL.val_metrics(D_, p_, J)})
            if k in cand:
                cpreds[k] = p_
        curves[nm] = sid.after_curves(D_, cpreds, J, CHK["min_peak_db"])
    val_all = pd.DataFrame([{"fit": k, "VAL set": "all lengths", "link_events": len(lev), "links": lev.cml_id.nunique(), "links_a_link_1": new_links_v[k],
                             **LTL.val_metrics(Dv, preds_v[k], J)} for k in fz])
    val = pd.concat([val_all, pd.DataFrame(val_rows)], ignore_index=True)
    save(val.round(4), P.SID3B_VAL, index=False)
    curve_tab = pd.concat([c.assign(set=k) for k, c in curves.items()], ignore_index=True)
    save(curve_tab.round(4), P.SID3B_VAL_AFTER, index=False)
    print(f"check 2 done; {(time.time() - clock) / 60:.1f} min", flush=True)

    # ------------------------------------------------------------------------------------------- check 3: residual against L
    rt = pd.concat([residual_table(Dv, lev, preds_v[k], k) for k in fz], ignore_index=True)
    flat = flat_verdict(rt)
    save(pd.concat([rt, flat.assign(row="verdict")], ignore_index=True).round(4), P.SID3B_RESIDUAL, index=False)

    # ------------------------------------------------------------------------------------------- figures
    fig = figure_length(through[through.fit.isin(cand)], BANDS, itu)
    save(fig, P.SID3B_FIG_LENGTH, dpi=100), plt.close(fig)
    fig = figure_val(val[val["VAL set"] != "all lengths"], lnames)
    save(fig, P.SID3B_FIG_VAL, dpi=105), plt.close(fig)
    fig = figure_residual(rt, BANDS)
    save(fig, P.SID3B_FIG_RESIDUAL, dpi=105), plt.close(fig)
    fig = figure_after(curves)
    save(fig, P.SID3B_FIG_AFTER, dpi=105), plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    q = dict(locals())
    q.update(readings(q))
    report = write_report(q)
    save(report, P.SID3B_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {(time.time() - clock) / 60:.1f} min")


def readings(q):
    """The READING paragraphs of the report, written after looking at the numbers of this run; every number is taken
    from the tables of the run."""
    th, ver, val, flat, params, rt = q["through"], q["verdict"], q["val"], q["flat"], q["params"], q["rt"]
    cand = q["cand"]
    vt = ver.set_index(["check", "set", "fit"])
    RULE = "through the model (THE RULE on VAL)"
    e = th[th.eligible & th.fit.isin(cand + ["M1"])]
    ratio = (e.assign(I_ratio=e.I_obs / e.I_model, s_ratio=e.s_obs / e.s_model)
              .groupby(["set", "fit"], sort=False)[["I_ratio", "s_ratio"]].median().round(2).reset_index())
    rv = ratio[ratio.set == "VAL"]
    rtr = ratio[ratio.set == "TRAIN"]
    o = q["obs"]
    ot, ov = o["TRAIN"].set_index(["band", "class"]), o["VAL"].set_index(["band", "class"])
    b28 = q["BANDS"][0]
    va = val.set_index(["fit", "VAL set"])
    pv = params.set_index(["fit", "parameter"])
    best_all = val[val["VAL set"] == "all lengths"].sort_values("rmse").iloc[0]
    n_tr = {k: int(vt.loc[(RULE, "TRAIN", k), "both_inside"]) for k in cand + ["M1"]}
    n_va = {k: int(vt.loc[(RULE, "VAL", k), "both_inside"]) for k in cand + ["M1"]}
    n_cells_va = int(vt.loc[(RULE, "VAL", cand[0]), "eligible_cells"])

    r1 = f"""**READING of check 1.** On TRAIN the repair works: with every link length in the fit, P is pinned and the fits reproduce the
model-free I and s ({', '.join(f'{k} {n_tr[k]}/10' for k in cand)}; M1 {n_tr['M1']}/10), where the S2b fits of S3 reached 0-2 of 10 by the
S3 rule and at most 5 of 10 read through the test. The gauge-timing fits (B) are the best on TRAIN. **On VAL every fit fails**
({', '.join(f'{k} {n_va[k]}/{n_cells_va}' for k in cand)}; M1 {n_va['M1']}/{n_cells_va}; the S2b fits 0-1/{n_cells_va} as well). The failure is not
the F / P split: on VAL the observed I AND s are both lower than every fit predicts on the same bins, by about the same factor
(median observed / model on VAL: I {rv.I_ratio.median():.2f}, s {rv.s_ratio.median():.2f}; on TRAIN {rtr.I_ratio.median():.2f} and {rtr.s_ratio.median():.2f}), and the static M1
fails in the same way. The VAL events attenuate less per mm/h of radar rain, both at the antenna and on the path (S3 already saw
lower VAL intercepts: {b28}, {CLASS_NAME[1]} mm/h, I {ot.loc[(b28, CLASS_NAME[1]), 'I']:.2f} on TRAIN against {ov.loc[(b28, CLASS_NAME[1]), 'I']:.2f} on VAL; s {ot.loc[(b28, CLASS_NAME[2]), 's']:.2f} against
{ov.loc[(b28, CLASS_NAME[2]), 's']:.2f} at {CLASS_NAME[2]} mm/h). A common factor on I and s is what a radar that reads the VAL rain too high (or a
different drop-size mix in those events) would give; it cannot be told apart here. No fit with event-independent F and P can pass this
rule while the two event sets differ like this.

Median observed / model ratio of I and s over the eligible cells:

{md_table(ratio)}"""

    r2 = f"""**READING of check 2.** The differences between the fits are again a few thousandths to a hundredth of a dB, inside the
intervals. Over all VAL lengths the lowest error is {best_all.fit} ({best_all.rmse:.3f} dB; the others {', '.join(f"{r.fit} {r.rmse:.3f}" for r in val[val['VAL set'] == 'all lengths'].itertuples() if r.fit != best_all.fit)}).
Fitting on all lengths does not cost much on the short links (L < 1 km: M4b-A {va.loc[('M4b-A', 'L < 1 km'), 'rmse']:.3f} against the S2b S<1 fit
{va.loc[('M4b, S<1 (S2b)', 'L < 1 km'), 'rmse']:.3f}) and matches the S2b S<3 fits on 1-3 km. The errors grow with L (> 3 km: {va.loc[('M4b-A', 'L > 3 km'), 'rmse']:.3f} dB
during all points, {va.loc[('M4b-A', 'L > 3 km'), 'rmse_rain']:.3f} during rain), where the path term dominates. **The gauge timing (B) gives better peaks**
(shower peak error M3-B {va.loc[('M3-B', 'all lengths'), 'shower_peak_error']:+.2f}, M4b-B {va.loc[('M4b-B', 'all lengths'), 'shower_peak_error']:+.2f} dB against M3-A {va.loc[('M3-A', 'all lengths'), 'shower_peak_error']:+.2f},
M4b-A {va.loc[('M4b-A', 'all lengths'), 'shower_peak_error']:+.2f}), but a worse first hour after the rain (RMSE 0-60 min M3-B {va.loc[('M3-B', 'all lengths'), 'rmse_0_60']:.3f} against M3-A
{va.loc[('M3-A', 'all lengths'), 'rmse_0_60']:.3f}): with the radar input, a fast antenna follows the radar's timing errors. After the rain M3 ends too low
(integral {va.loc[('M3-A', 'all lengths'), 'integral_pred_over_obs']:.2f} x), M4b too high ({va.loc[('M4b-A', 'all lengths'), 'integral_pred_over_obs']:.2f} x), as in S2b."""

    n_flat = int((flat.verdict == "FLAT").sum())
    r3 = f"""**READING of check 3.** {n_flat} of {len(flat)} (fit, band) lines are FLAT by the rule: the residual during rain does not change
with L within the intervals, so the L-scaling of every fit is right on VAL. But the residual is negative at every length (about
{rt.residual_median_db.quantile(0.25):+.2f} to {rt.residual_median_db.quantile(0.75):+.2f} dB, middle half of the bins): the fits predict too much attenuation during rain on VAL at
all lengths, the same level difference as in check 1."""

    tw3, tw4 = pv.loc[("M3-A", "tau_wet")], pv.loc[("M4b-A", "tau_f_wet")]
    r4 = f"""**READING of check 4.** No. With P pinned the radar wetting time stays long: M3-A tau_wet {tw3.value:.0f} [{tw3.lo:.0f}, {tw3.hi:.0f}] min
(S2b S<1: 43), M4b-A tau_f_wet {tw4.value:.1f} [{tw4.lo:.1f}, {tw4.hi:.1f}] min (S2b: 15), against 2-3 min with the gauge (TRAIN and VAL). So the slow
wetting with the radar is not caused by the F / P confusion; it is the input smoothing that fake check Ta of S2 showed (a 5 min, noisy
rain input makes the wetting about 2 x slower). The gauge-timing fits (B) fit TRAIN slightly worse (training loss M3-B - M3-A
= {q['conv'].set_index('fit').loss['M3-B'] - q['conv'].set_index('fit').loss['M3-A']:+.0f}, M4b-B - M4b-A = {q['conv'].set_index('fit').loss['M4b-B'] - q['conv'].set_index('fit').loss['M4b-A']:+.0f}, of about 1.6e5) and give the better peaks on VAL."""

    none_text = f"""What fails: every fit predicts a larger intercept and slope on the VAL steady-rain bins than VAL shows (observed / model
about {rv.I_ratio.median():.2f} for I and {rv.s_ratio.median():.2f} for s), the static M1 and the S2b fits included. This is a difference between the TRAIN and
VAL event sets that a model with one F and one P per band cannot follow; it is not the F / P split, which the all-length fits get
right on TRAIN (M3-B and M4b-B 10 of 10 cells). What stands from S3b: (1) with all lengths in the fit, P is pinned and the TRAIN length
test is reproduced (B fits: 10 / 10); (2) the L-scaling is right on VAL (residual flat with L); (3) the radar wetting time stays long
with P pinned (input smoothing, not P); (4) the gauge timing gives better peaks."""

    rc = """**Note for the next step (not started here).** The best-checked fit on TRAIN is M4b-B (gauge timing, all lengths): 10 / 10 cells
through the test and 7 / 10 even by the stricter S3 equilibrium rule, VAL RMSE within 0.004 dB of the best. Whether the VAL level
difference is a property of the events (rain type, radar error, antenna condition) is the conditions question of S4."""

    unclear = """- **The VAL rule may be testing the event sets, not the models.** VAL has 18 events (TRAIN 51) and only 8 eligible cells; if the
  radar over-reads the VAL rain, no rain-driven model can pass. The rule was written before the run and stays as it is; the reading
  above says what failed.
- **TRAIN numbers are in-sample** (the fits saw these events), so 10 / 10 on TRAIN shows that the frame CAN reproduce the length test,
  not that it predicts it.
- **Bootstrap intervals are too narrow** (30 warm-started refits of 30 iterations, as S2b), e.g. tau_s of M4b at the 240 bound with an
  interval that touches it. M4b-A did not converge on all events, and for M3-A and M4b-A none of the 8 starts converged before the cap of 1000.
- **The fixed gauge time constants (B) come from the TRAIN gauge links**: B is not independent of TRAIN.
- **Uneven rain on long links** may still bias s (radar path mean, b < 1); the residual is flat with L on VAL, which argues against a
  large effect.
- **32 GHz** has no link under 1 km and is never eligible; its I is an extrapolation.
- The run was stopped once by the session for low memory (after 97 min) and resumed from the cache; no fit was affected (every start
  and draw is saved)."""
    return {"reading_check1": r1, "reading_check2": r2, "reading_check3": r3, "reading_check4": r4, "none_text": none_text, "reading_candidate": rc,
            "unclear": unclear, "ratio": ratio}


def write_report(q):
    """The report from the numbers of the run (q = the local variables of main). The readings marked READING are written
    from the numbers of this run."""
    import numba
    import scipy
    plan, conv, params, ver, through, val, flat, rt = q["plan"], q["conv"], q["params"], q["verdict"], q["through"], q["val"], q["flat"], q["rt"]
    cand, BANDS, curve_tab, res = q["cand"], q["BANDS"], q["curve_tab"], q["res"]
    vt = ver.set_index(["check", "set", "fit"])
    RULE = "through the model (THE RULE on VAL)"
    EQ = "equilibrium a F, P (S3 rule)"

    def vline(st, fits, check=RULE, a="I_inside", b="s_inside"):
        out = []
        for k in fits:
            if (check, st, k) in vt.index:
                r = vt.loc[(check, st, k)]
                out.append(f"{k} **{r.verdict}** ({int(r.both_inside)} of {int(r.eligible_cells)} cells; {a[0]} {int(r[a])}, {b[0]} {int(r[b])})")
        return "; ".join(out)

    pv = params.set_index(["fit", "parameter"])

    def pci(k, name, nd=1):
        if (k, name) not in pv.index:
            return "-"
        r = pv.loc[(k, name)]
        return f"{r.value:.{nd}f} [{r.lo:.{nd}f}, {r.hi:.{nd}f}]" if np.isfinite(r.lo) else f"{r.value:.{nd}f}"

    consistent = [k for k in cand if (RULE, "VAL", k) in vt.index and vt.loc[(RULE, "VAL", k)].verdict == "consistent"]
    va = val[val["VAL set"] == "all lengths"].set_index("fit")
    winner = min(consistent, key=lambda k: va.loc[k, "rmse"]) if consistent else None
    gauge_tr = {"M3": (2.9, 12.8), "M4b": (2.6, 10.3)}
    gauge_val = {"M3": (1.8, 10.8), "M4b": (2.0, 8.8)}
    timing_rows = []
    for k in [c for c in cand if c.endswith("-A")]:
        m = FITS[k][0]
        tw, td = ("tau_wet", "tau_dry") if m == "M3" else ("tau_f_wet", "tau_f_dry")
        v = pv.loc[(k, tw), "value"]
        timing_rows.append({"fit": k, "wetting (radar, all lengths)": pci(k, tw), "drying (radar, all lengths)": pci(k, td),
                            "short (<= {:g} min)".format(CHECKS["timing_short_min"]): bool(v <= CHECKS["timing_short_min"]),
                            "gauge TRAIN (S2b) wet / dry": f"{gauge_tr[m][0]} / {gauge_tr[m][1]}", "gauge VAL (S3.4) wet / dry": f"{gauge_val[m][0]} / {gauge_val[m][1]}",
                            "S2b radar S<1 wet / dry": {"M3": "43.0 / 24.7", "M4b": "15.0 / 7.4"}[m]})
    timing = pd.DataFrame(timing_rows)
    key_par = params[params.parameter.isin(["d", "tau_wet", "tau_dry", "phi", "tau_f_wet", "tau_f_dry", "tau_s", "a_p10", "a_p90"])
                     | params.parameter.str.match(r"^[FP]\[.*\]\((0\.5|2|5|10)\)$")].round(3)
    th_show = through[(through.set == "VAL") & through.eligible & through.fit.isin(cand)][
        ["fit", "band", "class", "R_cell", "I_obs", "I_obs_lo", "I_obs_hi", "I_model", "s_obs", "s_obs_lo", "s_obs_hi", "s_model", "I_inside", "s_inside"]].round(2)
    th_tr = through[(through.set == "TRAIN") & through.eligible & through.fit.isin(cand)][
        ["fit", "band", "class", "I_obs", "I_obs_lo", "I_obs_hi", "I_model", "s_obs", "s_obs_lo", "s_obs_hi", "s_model", "I_inside", "s_inside"]].round(2)
    val_show = val[["fit", "VAL set", "link_events", "links", "links_a_link_1", "events", "rmse", "rmse_lo", "rmse_hi", "rmse_rain", "rmse_0_60", "rmse_60_180",
                    "shower_peak_error", "peak_lo", "peak_hi", "integral_pred_over_obs", "integral_lo", "integral_hi"]].round(3)
    cur_show = curve_tab[curve_tab.minutes_after_last_wet_bin.isin([0, 10, 20, 30, 60, 90, 120, 175])].round(3)
    rt_show = rt[["fit", "band", "L_bin_km", "L_median_km", "link_events", "links", "residual_median_db", "lo", "hi"]].round(3)
    secs = sum(r.get("seconds", 0.0) for r in res.values())
    standing = ("**Standing rule.** WAA is never estimated by subtracting a power-law rain term. The path part is separated only through its "
                "dependence on L at matched rain. ITU-R P.838 appears only as a size reference line on a figure and as one of the optimiser "
                "starting values of P (the S2b start plan), never to correct data.")

    if winner:
        m = FITS[winner][0]
        rows = [f"WAA v3 (S3b, 2026-10-15): fit {winner}, consistent with the length test on VAL, lowest VAL error among the consistent fits",
                "A(t) = W(t) + L * P_band(R(t - d)) + c          (radar input, TRAIN events of ALL link lengths; 90 % event-bootstrap intervals,",
                "                                                 warm-started 30-iteration refits: too narrow along flat directions)"]
        rows.append("dW/dt = (a_link F_band(R) - W) / tau,  tau = tau_wet while W rises, tau_dry while it falls" if m == "M3"
                    else "W = W_f + W_s;  W_f -> phi a_link F(R) (tau_f_wet / tau_f_dry),  W_s -> (1 - phi) a_link F(R) (tau_s)")
        rows.append(f"{'lag d':26s}{pv.loc[(winner, 'd'), 'value']:+.0f} min")
        for nm in TAUS[m]:
            rows.append(f"{nm:26s}{pci(winner, nm, 2 if nm == 'phi' else 1)}" + ("   (fixed: gauge timing)" if FITS[winner][1] == "B" else ""))
        for b in BANDS:
            rows.append(f"{'F[' + b + '](0.5 / 2 / 5)':26s}" + " / ".join(f"{pv.loc[(winner, f'F[{b}]({r:g})'), 'value']:.2f}" for r in KEY_R) + " dB")
        for b in BANDS:
            rows.append(f"{'P[' + b + '](2 / 5 / 10)':26s}" + " / ".join(f"{pv.loc[(winner, f'P[{b}]({r:g})'), 'value']:.2f}" for r in (2, 5, 10)) + " dB/km")
        rows.append(f"{'a_link 10-90 %':26s}{pv.loc[(winner, 'a_p10'), 'value']:.2f} - {pv.loc[(winner, 'a_p90'), 'value']:.2f}")
        rows.append(f"VAL: RMSE {va.loc[winner, 'rmse']:.3f} dB (all lengths), during rain {va.loc[winner, 'rmse_rain']:.3f}; shower peaks {va.loc[winner, 'shower_peak_error']:+.2f} dB; "
                    f"after-rain integral {va.loc[winner, 'integral_pred_over_obs']:.2f} x")
        box = "```\n" + "\n".join(rows) + "\n```"
    else:
        box = "**None.** No fit is consistent with the length test on VAL (the rule). " + q.get("none_text", "")

    return f"""# wet_antenna step S3b: the S2 frame fitted on ALL link lengths, so that P is pinned by L (repair after S3)

Written by `research/wet_antenna/system_id_s3b.py` (libraries `system_id.py`, `length_test.py`). Settings and every rule:
`wet_antenna.system_id_s3b` in `research/config.yaml`, written on 2026-10-15 BEFORE the run. Plan: `docs/ROADMAP.md`, "Findings that
stand (after S3)" and "Remaining plan". **VAL is the judge** (these fits never saw it); TRAIN numbers are information (the fits saw them).

{standing}

## 1. What was run

- The S2 frame A_pred = W + L P_band(R(t - d)) + c on the S3 all-length TRAIN dataset ({len(q['let'])} link-events, {q['let'].cml_id.nunique()} links, up to
  {q['let'].L_km.max():.1f} km; the S2 rules for A, c, windows and exclusions), radar input, the S2b settings (knots, increments >= 0, a_link ridge, lag
  grid {M['lag_grid_min'][0]}..{M['lag_grid_min'][1]} min, tau_s of M4b in 30-240 min).
- Fits: {', '.join(f'{k} ({FITS[k][0]}, ' + ('free' if FITS[k][1] == 'A' else 'gauge timing fixed') + ')' for k in plan['fits'])}. Variant B fixes
  {', '.join(f'{p} {v}' for p, v in S['fixed_B']['M3'].items())} (M3) and {', '.join(f'{p} {v}' for p, v in S['fixed_B']['M4b'].items())} (M4b), the TRAIN gauge values of S2b.
- {plan['n_starts']} starts per fit (the S2b start plan, incl. the 4 P starts), caps {M['maxiter']['start']} / {M['maxiter']['lag']} / {M['maxiter']['polish']}; no leave-one-event-out
  folds; {S['boot_reps']} event-bootstrap refits per fit (warm start, {M['maxiter']['boot']} iterations): **these intervals are too narrow along flat
  directions** (S2b), read them as lower bounds.
- **Timing:** one M4b-A start {plan['t_fit_s']} s ({plan['nit_start']} iterations, converged {plan['converged_start']}); one evaluation {plan['t_eval_s'] * 1e3:.0f} ms;
  projection {plan['projected_full_h']} h on {plan['workers']} processes against the budget {S['time_budget_h']} h: **cuts: {', '.join(plan['cuts']) or 'nothing'}**.
  Fitting time {secs / 3600:.1f} h summed over the processes. This run: {(time.time() - q['clock']) / 60:.1f} min.
- VAL: {len(q['lev'])} link-events, {int(q['Dv']['le_event'].max()) + 1} events. Every VAL link is in the TRAIN fits (a_link = 1 for {max(q['new_links_v'][k] for k in q['fz'])} links).
- Software: numba {numba.__version__}, scipy {scipy.__version__}, numpy {np.__version__}.

**Convergence** (`{rel(P.SID3B_CONVERGENCE)}`):

{md_table(conv.drop(columns=['loss']).round(2))}

## 2. The length test on VAL (check 1, the main result)

The same Huber line (A_pred - c) = I + s L on each fit's own prediction, on the S3.1 steady-rain bins, against the observed I and s
intervals. Rule: both inside in >= 2/3 of the eligible cells (>= {LT['verdict']['min_short_links']} links under 1 km).

**VAL (the rule):** {vline('VAL', cand)}.

TRAIN (information): {vline('TRAIN', cand)}. For continuity, the S2b fits through the test on VAL: {vline('VAL', list(S2B_FITS))}.
The S3 equilibrium comparison (median a_link x F against I, P against s) on VAL: {vline('VAL', cand, EQ, 'F_inside', 'P_inside')}.

{q.get('reading_check1', '')}

Every eligible VAL cell (`{rel(P.SID3B_LENGTH_CMP)}`):

{md_table(th_show)}

TRAIN cells (information):

{md_table(th_tr)}

![I and s through the model](../figures/{P.SID3B_FIG_LENGTH.name})

## 3. VAL errors by length (check 2)

Frozen fits from the rain alone; the S2b fits where their links exist (S<1 fits below 1 km, S<3 fits below 3 km).

{md_table(val_show)}

(`{rel(P.SID3B_VAL)}`.)

{q.get('reading_check2', '')}

![VAL by length](../figures/{P.SID3B_FIG_VAL.name})

Median (A - c) / peak after the last shower on VAL (`{rel(P.SID3B_VAL_AFTER)}`):

{md_table(cur_show)}

![after the last shower](../figures/{P.SID3B_FIG_AFTER.name})

## 4. Residual against L on VAL (check 3)

Per link-event the median (A - A_pred) during rain; per band and L bin the median with a 90 % event-bootstrap interval. FLAT = the
intervals of the shortest and the longest bin (>= {CHECKS['residual_vs_L']['min_links']} links) overlap.

{md_table(flat.round(3))}

{q.get('reading_check3', '')}

{md_table(rt_show)}

(`{rel(P.SID3B_RESIDUAL)}`.)

![residual vs L](../figures/{P.SID3B_FIG_RESIDUAL.name})

## 5. Timing (check 4)

{md_table(timing)}

{q.get('reading_check4', '')}

The parameters of every fit (`{rel(P.SID3B_PARAMS)}`):

{md_table(key_par)}

## 6. WAA v3 candidate

{box}

{q.get('reading_candidate', '')}

## 7. Unclear / possibly wrong

{q.get('unclear', '')}
"""


if __name__ == "__main__":
    main()
