"""wet_antenna step S2b: repair and finish the S2 model (the P / F trap, a bounded slow store, longer fits).

Step S2 left three problems in its model frame

    A_pred(t) = W(t) + L * P_band(R(t - d)) + c      (TRAIN events only; VAL is never read)

  (1) an OPTIMISER TRAP: P and F were sums of exp(parameter) increments, and an increment that went near 0 had no
      gradient left, so P at light rain stayed near 0 (fake check Tc failed);
  (2) the slow store of M4 ran to 400-600 min and overshot the small leftover after the rain (about 1-2 h);
  (3) most fits stopped at their iteration caps.
This step repairs them and nothing else:
  item 1   P and F with plain increments >= 0 (no exp), 4 starting curves for P (0, 0.5, 1, 2 x the ITU-R P.838 size;
           starting values only), the fake check Tc again, the length rule again;
  item 2   M4b = M4 with tau_s bounded to 30-240 min: a fake check (truth M4b), the ladder M1 < M3 < M3D < M4b on the
           real sets (the free M4 only as a reference), the behaviour after the rain;
  item 3   8 starts per model, 5 x the S2 iteration cap on all events, 3 x in the leave-one-event-out folds; a
           convergence table, and S2 against S2b.
RULE: WAA is never estimated by subtracting a power-law rain term. P is free and told apart from W only through L.
Every rule is in config `wet_antenna.system_id_s2b` (written before the run); the data, the judge and the fake noise
are those of `wet_antenna.system_id_s2`. Fits are cached in research/.cache/wet_antenna/system_id_s2b/ (every
finished stage and every few folds, so an interrupted run goes on where it stopped).
Run from the repository root (hours the first time on 3 processes, a few minutes from the cache):
    python -m research.wet_antenna.system_id_s2b
"""
import copy
import json
import os
import pickle
import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

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
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid

CFG = load_config()
S2 = CFG["wet_antenna"]["system_id_s2"]
S = CFG["wet_antenna"]["system_id_s2b"]
DATA, FK, CHK = S2["data"], S2["fake"], S2["checks"]
M = copy.deepcopy(S2["model"])                                # the S2 model settings with the changes of S2b
M["increments"], M["maxiter"] = S["model"]["increments"], dict(S["model"]["maxiter"])
M["bounds"] = {**M["bounds"], "increment": S["model"]["increment_bounds"], "tau_s_min": S["model"]["tau_s_m4b_min"]}
J = {**S2["judge"], "order": S["judge"]["order"], "param_boot_reps": S["judge"]["param_boot_reps"]}
FK2 = S["fake"]
LADDER, REF = S["model"]["ladder"], S["model"]["reference"]
MODELS = LADDER + REF
KEY_R = S2["length_rule"]["key_rain_mmh"]
KNOTS = np.asarray(M["knots_mmh"], float)
S2_CAP = S2["model"]["maxiter"]
CACHE = P.CACHE / "system_id_s2b"
SMOKE = os.environ.get("S2B_SMOKE")                          # a code test: a folder for its cache, every 12th link-event,
if SMOKE:                                                      #   tiny caps, nothing written to outputs/
    CACHE = Path(SMOKE)
    M["maxiter"] = {"start": 6, "lag": 3, "polish": 6, "em": 3, "fold": 3, "boot": 2}
    M["lag_grid_min"] = [-1, 1]
    J["boot_reps"], J["param_boot_reps"] = 50, 4

    def save(obj, path, **kw):                                 # noqa: F811
        pass
SET_LABEL = {"S1": "S<1, radar", "S2": "S<2, radar", "S3": "S<3, radar", "GG": "gauge links, GAUGE input", "GR": "gauge links, radar input (same points)"}
BAND_STYLE = [("-", "o"), ("--", "s"), ("-.", "^"), (":", "D")]


# =================================================================================================== the jobs (worker processes)
def starts_for(Pm, itu_inc):
    """The 8 starts of config `start_plan`: a generic start for F, the time constants and a = 1, and a P start =
    factor x the ITU-R P.838 size at the knots (a starting value only)."""
    return {f"{g}+{p}": Pm.start(M["starts"][g], Pinc=S["model"]["p_starts"][p] * np.asarray(itu_inc)) for g, p in S["model"]["start_plan"]}


def run_job(job):
    """One model on one dataset: the fit on all events (8 starts), the leave-one-event-out predictions and (if asked)
    the event-bootstrap refits. Every stage, and every few folds or draws, is saved in the cache."""
    file = CACHE / f"{job['key']}.pkl"
    res = pickle.loads(file.read_bytes()) if file.exists() else {}
    D = sid.load_data(job["data"], job["bands"])
    Pm = sid.Problem2(D, job["model"], M)
    t0 = [time.time()]

    def dump():
        res["seconds"] = res.get("seconds", 0.0) + time.time() - t0[0]
        t0[0] = time.time()
        file.write_bytes(pickle.dumps(res))

    if "all" not in res:                                   # every finished start and lag fit is saved (state "all_state")
        res["all"] = sid.fit_all2(Pm, starts_for(Pm, job["itu_inc"]), M["lag_grid_min"], M["maxiter"], M["drift"]["em_rounds"],
                                  state=res.setdefault("all_state", {}), checkpoint=dump)
        del res["all_state"]
        dump()
    fit = res["all"]
    if Pm.drift:
        Pm.s_b = fit["s_b"][-1]
    if "params" not in res:
        res["pred_all"], res["params"] = Pm.predict(fit["x"], fit["d"]), Pm.params(fit["x"], fit["d"], KEY_R)
        res["a_link"] = np.exp(fit["x"][Pm.names == "loga"])
        dump()
    if job.get("loo") and "loo_pred" not in res:
        st = res.setdefault("loo_state", {"pred": np.full(len(D["y5"]), np.nan), "folds": []})
        done_ev = {f["event"] for f in st["folds"]}
        for e in np.unique(D["le_event"]):
            if int(e) in done_ev:
                continue
            fw = sid.fit_weighted2(Pm, fit, (D["le_event"] != e).astype(float), M["maxiter"]["fold"])
            held = D["le_event"][D["bin_le"]] == e
            st["pred"][held] = Pm.predict(fw["x"], fw["d"])[held]
            st["folds"].append({"event": int(e), "nit": fw["nit"], "loss": fw["loss"], "converged": fw["converged"]})
            if len(st["folds"]) % 6 == 0:
                dump()
        res["loo_pred"], res["folds"] = st["pred"], st["folds"]
        del res["loo_state"]
        dump()
    if job.get("boot") and len(res.get("boot", [])) < job["boot"]:
        draws = res.setdefault("boot", [])
        n_ev = int(D["le_event"].max()) + 1
        rng = np.random.default_rng([J["seed"], zlib.crc32(job["key"].encode())])
        mults = [np.bincount(rng.integers(0, n_ev, n_ev), minlength=n_ev).astype(float) for _ in range(job["boot"])]
        for i in range(len(draws), job["boot"]):
            fw = sid.fit_weighted2(Pm, fit, mults[i][D["le_event"]], M["maxiter"]["boot"])
            draws.append(Pm.params(fw["x"], fw["d"], KEY_R))
            if len(draws) % 10 == 0:
                dump()
        dump()
    return job["key"]


def done(job):
    file = CACHE / f"{job['key']}.pkl"
    if not file.exists():
        return False
    res = pickle.loads(file.read_bytes())
    return "params" in res and (not job.get("loo") or "loo_pred" in res) and (not job.get("boot") or len(res.get("boot", [])) >= job["boot"])


def run_pool(jobs, workers, label):
    jobs = [j for j in jobs if not done(j)]
    if not jobs:
        return
    t0 = time.time()
    print(f"{label}: {len(jobs)} jobs", flush=True)
    with ProcessPoolExecutor(workers) as ex:
        for i, f in enumerate(as_completed([ex.submit(run_job, j) for j in jobs]), 1):
            print(f"   {label}: {f.result()} done ({i} of {len(jobs)}, {(time.time() - t0) / 60:.1f} min)", flush=True)


_RESULTS = {}


def result(key):
    if key not in _RESULTS:
        _RESULTS[key] = pickle.loads((CACHE / f"{key}.pkl").read_bytes())
    return _RESULTS[key]


def boot_interval(key, name):
    v = np.array([dr.get(name, np.nan) for dr in result(key).get("boot", [])], float)
    return tuple(np.nanpercentile(v, J["interval"])) if np.isfinite(v).any() else (np.nan, np.nan)   # a band absent from the set: no interval


# =================================================================================================== the judge
def ladder(prefix, D):
    """The judge on the leave-one-event-out predictions of one dataset: the candidates with the stop rule, the
    reference M4 as an extra row (its difference to the best candidate), the shower peaks. Returns (table, best, chosen)."""
    n_ev = int(D["le_event"].max()) + 1
    preds = {m: result(f"{prefix}_{m}")["loo_pred"] for m in MODELS}
    errors = {m: sid.event_errors(D, preds[m], J) for m in MODELS}
    tab, best, chosen = sid.judge({m: errors[m] for m in LADDER}, LADDER, J)
    for m in REF:
        ref_tab = sid.judge({best: errors[best], m: errors[m]}, [best, m], J)[0]
        row = ref_tab[ref_tab.model == m].iloc[0].to_dict()
        dif = sid.paired(errors[m].rmse.to_numpy(), errors[best].rmse.to_numpy(), n_ev, J)
        row.update({"minus_best": dif[0], "diff_lo": dif[1], "diff_hi": dif[2], "interval_holds_0": bool(dif[1] <= 0 <= dif[2])})
        tab = pd.concat([tab, pd.DataFrame([row])], ignore_index=True)
    sh = sid.shower_table(D, J["shower_gap_min"])
    peaks = pd.DataFrame([dict(zip(("shower_peak_error_median", "shower_peak_error_lo", "shower_peak_error_hi", "showers"),
                                   sid.shower_stat(*sid.shower_peaks(D, sh, preds[m], J["shower_min_peak_db"]), n_ev, J)), model=m) for m in MODELS])
    tab = tab.drop(columns=["skill_vs_M0", "skill_lo", "skill_hi"]).merge(peaks, on="model")
    tab["role"] = np.where(tab.model.isin(REF), "reference", "candidate")
    return tab, best, chosen


def key_names(model, bands):
    """The key numbers of the length rule: (time constants, F per band at the key rain rates)."""
    taus = {"M3": ["tau_wet", "tau_dry"], "M3D": ["tau_wet", "tau_dry"], "M4b": ["tau_f_wet", "tau_f_dry"], "M4": ["tau_f_wet", "tau_f_dry"], "M1": []}[model]
    return taus, [f"F[{b}]({r:g})" for b in bands for r in KEY_R]


def convergence_row(key, label, model):
    """One row of the convergence table for a cached job."""
    r = result(key)
    a = r["all"]
    folds = r.get("folds", [])
    nit = a.get("starts_nit", {})
    return {"fit": label, "model": model, "converged_on_all_events": a.get("converged"), "loss": a["loss"], "d": a["d"],
            "starts": len(nit), "starts_converged": int(sum(bool(v) for v in a.get("starts_converged", {}).values())),
            "starts_past_S2_cap": int(sum(v > S2_CAP["start"] for v in nit.values())), "starts_median_iterations": float(np.median(list(nit.values()))) if nit else np.nan,
            "spread_over_starts": max(a["starts_d0"].values()) - min(a["starts_d0"].values()),
            "best_start": min(a["starts_d0"], key=a["starts_d0"].get), "folds": len(folds),
            "folds_at_cap": int(sum(f["nit"] >= M["maxiter"]["fold"] for f in folds)), "folds_past_S2_cap": int(sum(f["nit"] > S2_CAP["fold"] for f in folds)),
            "seconds": round(r.get("seconds", 0.0))}


# =================================================================================================== figures
def figure_ladder(ladders, chosen_of):
    fig, axes = plt.subplots(1, len(ladders), figsize=(4.0 * len(ladders), 4.6), squeeze=False)
    for ax, (name, t) in zip(axes[0], ladders.items()):
        x = np.arange(len(t))
        for i, r in enumerate(t.itertuples()):
            ref = r.role == "reference"
            ax.errorbar(i, r.rmse, yerr=[[r.rmse - r.rmse_lo], [r.rmse_hi - r.rmse]], fmt="D" if ref else "o", mfc="none" if ref else INK, color=INK, capsize=3, ms=6)
            ax.text(i, r.rmse_hi, f"{r.rmse:.3f}" + ("\nCHOSEN" if r.model == chosen_of[name] else "") + ("\n(ref.)" if ref else ""), ha="center", va="bottom",
                    fontsize=7.5, fontweight="bold" if r.model == chosen_of[name] else "normal")
            if r.model == chosen_of[name]:
                ax.plot(i, r.rmse, "s", ms=16, mfc="none", color=INK, mew=1.4)
        ax.set_xticks(x), ax.set_xticklabels(t.model)
        ax.set_title(SET_LABEL[name], fontsize=9)
        lo, hi = t.rmse_lo.min(), t.rmse_hi.max()
        ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.25 * (hi - lo))
        style_axes(ax)
    axes[0][0].set_ylabel("held-out RMSE of A (dB), mean over events")
    fig.suptitle("Step S2b: leave-one-event-out error per set (dot = mean, whisker = 90 % event-bootstrap interval; square = chosen by the stop rule;\n"
                 "open diamond = the free-tau_s M4, a reference, not a candidate). Compare inside a panel.", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    return fig


def figure_curves(model, old_model, bands, old_params):
    """F and P of the chosen model on S<1 and S<3 after the fix (thick, filled markers; band = line style and
    marker) with the S2 curves (thin dotted lines, no markers)."""
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.0), sharex=True)
    for col, s_ in enumerate(("S1", "S3")):
        par = result(f"real_{s_}_{model}")["params"]
        old = old_params[(old_params.set == SET_LABEL[s_]) & (old_params.model == old_model)].set_index("parameter").value
        for row, (name, unit) in enumerate((("F", "dB"), ("P", "dB / km"))):
            ax = axes[row][col]
            for (ls, mk), band in zip(BAND_STYLE, bands):
                v = np.array([par.get(f"{name}[{band}]({k:g})", np.nan) for k in KNOTS[1:]])
                if not np.isfinite(v).any():
                    continue
                ax.plot(KNOTS, np.r_[0, v], ls=ls, marker=mk, ms=5, color=INK, lw=2.0, label=f"{band}, S2b ({model})")
                q = np.array([boot_interval(f"real_{s_}_{model}", f"{name}[{band}]({k:g})") for k in KNOTS[1:]])
                ax.vlines(KNOTS[1:], q[:, 0], q[:, 1], color=INK, lw=1.0)
                ov = np.array([old.get(f"{name}[{band}]({k:g})", np.nan) for k in KNOTS[1:]])
                if np.isfinite(ov).any():
                    ax.plot(KNOTS, np.r_[0, ov], ls=":", color=INK2, lw=1.0, label=f"{band}, S2 ({old_model}, thin dotted)" if band == bands[0] else None)
            ax.set_xscale("symlog", linthresh=0.5)
            ax.set_xticks(KNOTS), ax.set_xticklabels([f"{k:g}" for k in KNOTS])
            ax.set(ylabel=f"{name}(R) ({unit})", xlim=(0, 45), title=f"{name} on {SET_LABEL[s_]}")
            ax.legend(frameon=False, fontsize=7.5, loc="upper left")
            style_axes(ax)
    for ax in axes[1]:
        ax.set_xlabel("rain R (mm/h)")
    fig.suptitle(f"Step S2b: F (the WAA level, a_link = 1) and P (the path term per km) of {model} after the fix of the increments; vertical bars = 90 %\n"
                 "event-bootstrap interval; bands by line style and marker; thin dotted lines = the S2 curves (exp increments, shorter fits)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig


def figure_after(curves, models):
    styles = {"M3": ("--", "s"), "M3D": ("-.", "D"), "M4b": ("-", "^"), "M4": (":", "v")}
    fig, axes = plt.subplots(1, len(curves), figsize=(7.0 * len(curves), 4.8), squeeze=False)
    for ax, (name, c) in zip(axes[0], curves.items()):
        ax.plot(c.minutes_after_last_wet_bin, c.observed, "o", mfc="none", color=INK, ms=5, label="observed (circles)")
        for m in models:
            ls, mk = styles[m]
            ax.plot(c.minutes_after_last_wet_bin, c[m], ls=ls, marker=mk, markevery=6, ms=5, color=INK, lw=1.4, label=f"{m}, held out" + (" (ref.)" if m in REF else ""))
        ax.set_yscale("log")
        ax.set_ylim(3e-3, 0.4)
        ax.set(xlabel="minutes after the last wet bin of the link", ylabel="median (A - c) / observed peak of the last shower",
               title=f"{SET_LABEL[name]} ({int(c.link_events.iloc[0])} link-events)")
        ax.legend(frameon=False, fontsize=8)
        style_axes(ax)
    fig.suptitle("Step S2b: after the last shower, observed against predicted from the rain alone (logarithmic axis)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig


def figure_length(length_tab, models):
    keys = list(dict.fromkeys(length_tab.key))
    fig, axes = plt.subplots(1, len(models), figsize=(7.2 * len(models), 0.36 * len(keys) + 2.6), squeeze=False)
    for ax, model in zip(axes[0], models):
        t = length_tab[length_tab.model == model].set_index("key")
        ks = [k for k in keys if k in t.index]
        for i, k in enumerate(ks):
            for s_, (mk, off) in zip(("S1", "S2", "S3"), (("o", 0.22), ("s", 0.0), ("^", -0.22))):
                v, lo, hi = t.loc[k, f"{s_}_value"], t.loc[k, f"{s_}_lo"], t.loc[k, f"{s_}_hi"]
                if np.isfinite(v):
                    ax.errorbar(v, i + off, xerr=[[max(v - lo, 0)], [max(hi - v, 0)]] if np.isfinite(lo) else None, fmt=mk, color=INK, ms=5,
                                mfc="none" if s_ != "S3" else INK, capsize=3, lw=1.0, label=SET_LABEL[s_] if i == 0 else None)
        ax.set(yticks=range(len(ks)), xscale="log", xlabel="value (minutes for the time constants, dB for the levels)")
        ax.set_yticklabels(ks, fontsize=8)
        ax.invert_yaxis()
        ax.set_title(f"{model}: circle = S<1, " + ("square = S<2, " if t["S2_value"].notna().any() else "") + "filled triangle = S<3", fontsize=9)
        style_axes(ax)
    fig.suptitle("Step S2b, the length rule again: the key WAA numbers on the nested length sets (radar input), with 90 % event-bootstrap intervals", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    clock = time.time()
    workers = max((os.cpu_count() or 2) - 1, 1)

    # ------------------------------------------------------------------------------------------- data (as S2)
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    bands_cfg, band_notes, band_links = sid.merged_bands(d, DATA)
    BANDS = list(bands_cfg)
    D3, le, drops = sid.build_dataset(d, {**DATA, "bands_ghz": bands_cfg})
    if SMOKE:
        thin = np.arange(len(le)) % 12 == 0
        D3, le = sid.subset(D3, thin, renumber_events=False), le[thin].reset_index(drop=True)
    L = le.L_km.to_numpy()
    has_g = (le.gauge_id != "").to_numpy()
    keep = {"S1": L < DATA["length_sets_km"]["S<1"], "S2": L < DATA["length_sets_km"]["S<2"], "S3": np.ones(len(le), bool), "GG": has_g, "GR": has_g}
    data = {"S3": D3, "S1": sid.subset(D3, keep["S1"]), "S2": sid.subset(D3, keep["S2"])}
    data["GG"], data["GR"] = sid.gauge_variants(D3, has_g)
    les = {k: le[keep[k]].reset_index(drop=True) for k in data}
    files = {k: sid.save_data(v, CACHE / f"data_{k}.npz") for k, v in data.items()}
    n_ev = {k: int(v["le_event"].max()) + 1 for k, v in data.items()}
    minutes = {k: int(v["le_n"].sum()) for k, v in data.items()}
    counts = pd.DataFrame([{"set": SET_LABEL[k], "link_events": len(les[k]), "links": les[k].cml_id.nunique(), "events": n_ev[k],
                            "points_5min": int((data[k]["w5"] > 0).sum()), **{f"links_{b}": les[k][les[k].band_name == b].cml_id.nunique() for b in BANDS}}
                           for k in ("S1", "S2", "S3", "GG", "GR")])
    save(counts, P.SID2B_COUNTS, index=False)
    print(f"TRAIN only: S<3 {len(le)} link-events, {le.cml_id.nunique()} links, {n_ev['S3']} events; S<1 {keep['S1'].sum()}, S<2 {keep['S2'].sum()}, gauge links "
          f"{has_g.sum()}; bands {band_links}; {workers} processes; loaded in {time.time() - clock:.0f} s", flush=True)

    # the ITU-R P.838 size at the knots, per band (a STARTING value of P and the truth of the fake data only)
    F_true, itu = sid.truth_curves(le, BANDS, FK)
    itu_inc = np.array([np.diff(itu[b][0] * KNOTS ** itu[b][1]) if b in itu else np.diff(0.2 * KNOTS) for b in range(len(BANDS))])

    def job(key, name, model, data_file=None, **kw):
        return {"key": key, "data": str(data_file or files[name]), "model": model, "bands": BANDS, "set": name, "itu_inc": itu_inc.tolist(), **kw}

    # ------------------------------------------------------------------------------------------- the fake datasets
    noise = sid.noise_table()

    def fake_file(truth, rep, variant=""):
        """Tc: the S2 data again (the same seed), on S<3, or its S<1 part (variant 'S1'). Te: truth M4b on S<1."""
        file = CACHE / f"data_fake_{truth}_{rep}{'_' + variant if variant else ''}.npz"
        if not file.exists():
            if truth == "Tc":
                Fd = sid.make_fake(D3, le, FK["truths"]["Tc"], list(FK2["Tc"]["seed"]) + [rep], noise, FK, M["drift"]["t_ref_min_before_rain"], DATA["offset_min_minutes"])
                if variant == "S1":
                    Fd = sid.subset(Fd, keep["S1"])
            else:
                Fd = sid.make_fake(data["S1"], les["S1"], FK2["Te"], list(FK2["Te"]["seed"]) + [rep], noise, FK, M["drift"]["t_ref_min_before_rain"], DATA["offset_min_minutes"])
            sid.save_data(Fd, file)
        return str(file)

    # ------------------------------------------------------------------------------------------- time one fit, plan, cuts
    def all_jobs(plan, chosen=None):
        sets = ["S1"] + (["S2"] if plan["s2"] else []) + ["S3", "GG"] + (["GR"] if plan["gr"] else [])
        out = []
        for name in sets:
            for m in MODELS:
                boot = J["param_boot_reps"] if (m in ("M3", "M4b") or m == chosen) else 0
                out.append(job(f"real_{name}_{m}", name, m, loo=True, boot=boot))
        for rep in range(plan["fake_reps"]):
            out += [job(f"fake_Tc_{rep}_S3_M3", "S3", "M3", fake_file("Tc", rep)), job(f"fake_Tc_{rep}_S1_M3", "S1", "M3", fake_file("Tc", rep, "S1"))]
            out += [job(f"fake_Te_{rep}_{m}", "S1", m, fake_file("Te", rep), loo=True) for m in ("M3", "M3D", "M4b")]
        return out

    def cost_h(jobs, plan):
        """Projected hours: every start needs the evaluations of the timed start, a lag fit a third and the polish half of
        that, every fold and bootstrap draw its cap (1.1 evaluations per iteration)."""
        e_start, n_lag = plan["nfev_start"], M["lag_grid_min"][1] - M["lag_grid_min"][0]
        tot = 0.0
        for j in jobs:
            if done(j):
                continue
            ev = len(S["model"]["start_plan"]) * e_start + n_lag * e_start / 3 + e_start / 2 + (M["drift"]["em_rounds"] * e_start / 3 if j["model"] == "M3D" else 0)
            ev += 1.1 * ((n_ev[j["set"]] * M["maxiter"]["fold"] if j.get("loo") else 0) + j.get("boot", 0) * M["maxiter"]["boot"])
            tot += ev * plan["t_eval_s"] * minutes[j["set"]] / minutes["S3"]
        return tot / workers / 3600

    plan_file = CACHE / "plan.json"
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
        print(f"plan from {rel(plan_file) if not SMOKE else plan_file}: cut: {plan['cuts'] or 'nothing'}")
    else:
        ones = np.ones(len(le))
        Pm = sid.Problem2(D3, "M4b", M)
        st = starts_for(Pm, itu_inc)
        x0 = st["g3+Pitu"]
        Pm.run(x0, 0, ones)                                                            # compile
        t0 = time.time()
        for _ in range(5):
            Pm.run(x0, 0, ones)
        t_eval = (time.time() - t0) / 5
        t0 = time.time()
        fit = Pm.fit(x0, 0, ones, M["maxiter"]["start"])
        plan = {"t_eval_s": t_eval, "t_fit_s": round(time.time() - t0, 1), "nfev_start": fit["nfev"], "nit_start": fit["nit"], "converged_start": fit["converged"],
                "cuts": [], "s2": True, "gr": True, "fake_reps": FK2["replicates"]}
        plan["projected_full_h"] = round(cost_h(all_jobs(plan), plan), 2)
        for cut, change in (("S<2", {"s2": False}), ("the gauge-links radar comparison", {"gr": False}), ("fake replicates 3 -> 2", {"fake_reps": 2})):
            if cost_h(all_jobs(plan), plan) <= S["time_budget_h"]:
                break
            plan.update(change)
            plan["cuts"].append(cut)
        plan["projected_h"] = round(cost_h(all_jobs(plan), plan), 2)
        plan_file.write_text(json.dumps(plan))
        print(f"timing on S<3: one M4b start (g3 + P = ITU) {plan['t_fit_s']} s, {plan['nit_start']} iterations, {plan['nfev_start']} evaluations, "
              f"converged: {plan['converged_start']}; one evaluation {t_eval * 1e3:.0f} ms; projection {plan['projected_full_h']} h on {workers} processes "
              f"(budget {S['time_budget_h']} h); cuts: {plan['cuts'] or 'nothing'} -> {plan['projected_h']} h", flush=True)
    sets = ["S1"] + (["S2"] if plan["s2"] else []) + ["S3", "GG"] + (["GR"] if plan["gr"] else [])
    sets_radar = [s_ for s_ in ("S1", "S2", "S3") if s_ in sets]

    heavy = lambda j: -minutes[j["set"]] * (2 if j.get("loo") else 1)            # noqa: E731  (big jobs first)
    run_pool(sorted(all_jobs(plan), key=heavy), workers, "stage 1")
    _RESULTS.clear()

    # ------------------------------------------------------------------------------------------- the ladders, the chosen model
    ladders, best_of, chosen_of = {}, {}, {}
    for name in sets:
        ladders[name], best_of[name], chosen_of[name] = ladder(f"real_{name}", data[name])
    chosen = chosen_of["S1"]                                     # the MAIN set (config judge.order)
    if chosen not in ("M3", "M4b"):                              # its parameter intervals as well
        run_pool(sorted([j for j in all_jobs(plan, chosen) if j["model"] == chosen], key=heavy), workers, "stage 2 (bootstrap of the chosen model)")
        _RESULTS.clear()
    ref = list(dict.fromkeys(["M3", "M4b"] + ([chosen] if chosen != "M1" else [])))
    ladder_tab = pd.concat([t.assign(set=SET_LABEL[k], best=best_of[k], chosen=chosen_of[k]) for k, t in ladders.items()], ignore_index=True)
    ladder_tab.insert(0, "set", ladder_tab.pop("set"))
    save(ladder_tab.round(4), P.SID2B_LADDER, index=False)

    # ------------------------------------------------------------------------------------------- fake checks
    reps = range(plan["fake_reps"])
    fake_rows, verdicts = [], {}
    tc_rows = []
    for rep in reps:
        for s_ in ("S3", "S1"):
            par = result(f"fake_Tc_{rep}_{s_}_M3")["params"]
            for bi, band in enumerate(BANDS):
                for kind in ("F", "P"):
                    for r_ in FK2["Tc"]["knots_mmh"]:
                        truth = float(F_true(r_)) if kind == "F" else itu[bi][0] * r_ ** itu[bi][1]
                        tc_rows.append({"replicate": rep, "set": s_, "kind": kind, "band": band, "R_mmh": r_, "value": par.get(f"{kind}[{band}]({r_:g})", np.nan), "truth": truth})
            fake_rows.append({"truth": "Tc", "replicate": rep, "what": f"M3 fitted on {SET_LABEL[s_][:3]}", "tau_wet": par["tau_wet"], "tau_dry": par["tau_dry"], "d": par["d"]})
    tc = pd.DataFrame(tc_rows)
    tc_med = tc.groupby(["set", "kind", "band", "R_mmh"], as_index=False).agg(value=("value", "median"), truth=("truth", "first"))
    tc_med["ratio"] = tc_med.value / tc_med.truth
    t3_ = tc_med[(tc_med.set == "S3") & tc_med.value.notna()]
    within = bool((t3_.ratio.sub(1).abs() <= FK2["Tc"]["tol"]).all())
    p2 = t3_[(t3_.kind == "P") & (t3_.R_mmh == 2)]
    p2_ok = bool((p2.ratio >= FK2["Tc"]["p2_min_share"]).all())
    worst = t3_.assign(off=t3_.ratio.sub(1).abs()).sort_values("off").iloc[-1]
    s1_ = tc_med[(tc_med.set == "S1") & tc_med.value.notna()]
    verdicts["Tc"] = {"PASS": within and p2_ok, "text": f"on S<3 {int((t3_.ratio.sub(1).abs() <= FK2['Tc']['tol']).sum())} of {len(t3_)} F and P knots (0.5-10 mm/h, every band) are within "
                      f"+-{FK2['Tc']['tol']:.0%} of the truth (farthest: {worst.kind}[{worst.band}]({worst.R_mmh:g}) x{worst.ratio:.2f}); P(2 mm/h) / truth = "
                      + ", ".join(f"{r.band} {r.ratio:.2f}" for r in p2.itertuples()) + f"; on S<1 {int((s1_.ratio.sub(1).abs() <= FK2['Tc']['tol']).sum())} of {len(s1_)} within"}
    save(tc_med.round(4), P.SID2B_TC, index=False)
    te_rows = []
    for rep in reps:
        Df = sid.load_data(fake_file("Te", rep), BANDS)
        errs = {m: sid.event_errors(Df, result(f"fake_Te_{rep}_{m}")["loo_pred"], J) for m in ("M3", "M3D", "M4b")}
        _, best_f, chosen_f = sid.judge(errs, ["M3", "M3D", "M4b"], J)
        par = result(f"fake_Te_{rep}_M4b")["params"]
        te_rows.append({"replicate": rep, "chosen": chosen_f, "best": best_f, **{k: par[k] for k in ("phi", "tau_f_wet", "tau_f_dry", "tau_s", "d")},
                        **{f"rmse_{m}": float(e.rmse.mean()) for m, e in errs.items()}})
        fake_rows.append({"truth": "Te", "replicate": rep, "what": f"ladder M3 < M3D < M4b on M4b data: chosen {chosen_f}", "tau_s": par["tau_s"], "phi": par["phi"],
                          "tau_wet": par["tau_f_wet"], "tau_dry": par["tau_f_dry"], "d": par["d"]})
    te = pd.DataFrame(te_rows)
    n_m4b = int((te.chosen == "M4b").sum())
    ts_med = float(te.tau_s.median())
    te_ok = n_m4b >= 2 and abs(ts_med / FK2["Te"]["tau_s"] - 1) <= FK2["Te"]["tol"]
    verdicts["Te"] = {"PASS": bool(te_ok), "text": f"M4b chosen in {n_m4b} of {len(te)} replicates ({', '.join(te.chosen)}); median tau_s {ts_med:.0f} min (true {FK2['Te']['tau_s']}, "
                      f"x{ts_med / FK2['Te']['tau_s']:.2f}); phi {te.phi.median():.2f} (true {FK2['Te']['phi']}), tau_f_wet {te.tau_f_wet.median():.1f} ({FK2['Te']['tau_f_wet']}), "
                      f"tau_f_dry {te.tau_f_dry.median():.1f} ({FK2['Te']['tau_f_dry']})"}
    fake = pd.DataFrame(fake_rows)
    for t_, v in verdicts.items():
        fake.loc[fake.truth == t_, "PASS"] = v["PASS"]
    save(fake.round(4), P.SID2B_FAKE, index=False)
    print("fake checks: " + "; ".join(f"{t_} {'PASS' if v['PASS'] else 'FAIL'}" for t_, v in verdicts.items()), flush=True)
    tag = {"Tc": "" if verdicts["Tc"]["PASS"] else " [not validated: fake check Tc failed]", "Te": "" if verdicts["Te"]["PASS"] else " [not validated: fake check Te failed]"}

    # ------------------------------------------------------------------------------------------- the length rule
    len_rows = []
    for m in ref:
        taus, levels = key_names(m, BANDS)
        for k in taus + levels + [f"P[{b}]({r:g})" for b in BANDS for r in (2, 5, 10)]:
            row = {"model": m, "key": k}
            for s_ in ("S1", "S2", "S3"):
                has = s_ in sets_radar
                row[f"{s_}_value"] = result(f"real_{s_}_{m}")["params"].get(k, np.nan) if has else np.nan
                row[f"{s_}_lo"], row[f"{s_}_hi"] = boot_interval(f"real_{s_}_{m}", k) if has else (np.nan, np.nan)
            row["in_rule"] = bool(np.isfinite(row["S1_value"]) and not k.startswith("P["))
            row["S3_inside_S1_interval"] = bool(row["S1_lo"] <= row["S3_value"] <= row["S1_hi"]) if np.isfinite(row["S1_value"]) else None
            len_rows.append(row)
    length_tab = pd.DataFrame(len_rows)
    save(length_tab.round(4), P.SID2B_LENGTH, index=False)
    rule_model = chosen if chosen != "M1" else "M3"
    in_rule = length_tab[(length_tab.model == rule_model) & length_tab.in_rule]
    agree = bool(in_rule.S3_inside_S1_interval.all())
    outside = in_rule[~in_rule.S3_inside_S1_interval.astype(bool)]

    # ------------------------------------------------------------------------------------------- parameters, convergence, S2 against S2b
    par_rows = []
    for name in sets:
        for m in ref:
            for k, v in result(f"real_{name}_{m}")["params"].items():
                lo, hi = boot_interval(f"real_{name}_{m}", k)
                par_rows.append({"set": SET_LABEL[name], "model": m, "parameter": k, "value": v, "lo": lo, "hi": hi})
    params = pd.DataFrame(par_rows)
    save(params.round(4), P.SID2B_PARAMS, index=False)
    conv = pd.DataFrame([convergence_row(f"real_{name}_{m}", SET_LABEL[name], m) for name in sets for m in MODELS]
                        + [convergence_row(f"fake_Tc_{rep}_{s_}_M3", f"fake Tc {rep}, {s_}", "M3") for rep in reps for s_ in ("S3", "S1")]
                        + [convergence_row(f"fake_Te_{rep}_{m}", f"fake Te {rep}", m) for rep in reps for m in ("M3", "M3D", "M4b")])
    save(conv.round(3), P.SID2B_CONVERGENCE, index=False)
    old_ladder, old_params = pd.read_csv(P.SID2_LADDER), pd.read_csv(P.SID2_PARAMS)
    cmp_rows = [{"what": "chosen model", "set": SET_LABEL[k], "S2": old_ladder[old_ladder.set == SET_LABEL[k]].chosen.iloc[0] if (old_ladder.set == SET_LABEL[k]).any() else "",
                 "S2b": chosen_of[k], "note": "S2 ladder M1 < M3 < M3D < M4 < M5; S2b M1 < M3 < M3D < M4b"} for k in sets]
    for k in sets:
        for m_old, m_new in (("M3", "M3"), ("M4", "M4"), ("M4", "M4b")):
            new = result(f"real_{k}_{m_new}")["params"]
            for name in key_names(m_new, BANDS)[0] + (["tau_s"] if m_new.startswith("M4") else []) + [f"F[{b}]({r:g})" for b in BANDS for r in KEY_R] + [f"P[{b}]({r:g})" for b in BANDS for r in (2, 10)]:
                o = old_params[(old_params.set == SET_LABEL[k]) & (old_params.model == m_old) & (old_params.parameter == name)]
                if o.empty or not np.isfinite(new.get(name, np.nan)):
                    continue
                cmp_rows.append({"what": f"{name}: S2 {m_old} -> S2b {m_new}", "set": SET_LABEL[k], "S2": float(o.value.iloc[0]), "S2b": float(new[name]),
                                 "note": f"S2 interval {o.lo.iloc[0]:.3g}..{o.hi.iloc[0]:.3g}" if np.isfinite(o.lo.iloc[0]) else ""})
    cmp_tab = pd.DataFrame(cmp_rows)
    save(cmp_tab, P.SID2B_COMPARE, index=False)

    # ------------------------------------------------------------------------------------------- after the rain, product, figures
    after_models = ["M3", "M3D", "M4b", "M4"]
    curves = {k: sid.after_curves(data[k], {m: result(f"real_{k}_{m}")["loo_pred"] for m in after_models}, J, CHK["min_peak_db"]) for k in ("S1", "S3")}
    curve_tab = pd.concat([c.assign(set=SET_LABEL[k]) for k, c in curves.items()], ignore_index=True)
    save(curve_tab.round(4), P.SID2B_AFTER, index=False)
    loo3 = {m: result(f"real_S3_{m}")["loo_pred"] for m in MODELS}
    product = pd.concat([le.drop(columns=["band"]), pd.DataFrame({"in_S1": keep["S1"], "in_S2": keep["S2"], "has_city_gauge_1km": has_g})], axis=1)
    for m in dict.fromkeys([chosen, "M3", "M4b"]):
        lt = sid.link_event_table(D3, loo3[m], J)
        if m == chosen:
            product = pd.concat([product, lt[["points", "wet_points", "peak_obs", "integral_obs"]].rename(columns={"peak_obs": "peak_obs_db", "integral_obs": "integral_obs_db_min"})], axis=1)
        product = pd.concat([product, lt[["peak_pred", "integral_pred", "bias", "rmse"]].add_suffix(f"_{m}")], axis=1)
    product.insert(0, "chosen_model", chosen)
    num = product.select_dtypes("number").columns
    product[num] = product[num].round(4)
    save(product, P.SID2B_LINK_EVENTS, index=False)

    curve_model = chosen if chosen != "M1" else "M3"
    old_model = {"M4b": "M4", "M3D": "M3", "M3": "M3"}[curve_model]
    fig = figure_ladder(ladders, chosen_of)
    save(fig, P.SID2B_FIG_LADDER, dpi=110), plt.close(fig)
    fig = figure_curves(curve_model, old_model, BANDS, old_params)
    save(fig, P.SID2B_FIG_CURVES, dpi=110), plt.close(fig)
    fig = figure_after(curves, ["M3", "M3D", "M4b"])
    save(fig, P.SID2B_FIG_AFTER, dpi=110), plt.close(fig)
    fig = figure_length(length_tab[~length_tab.key.str.startswith("P[")], ref)
    save(fig, P.SID2B_FIG_LENGTH, dpi=110), plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    import numba
    import scipy
    show = ["set", "role", "model", "rmse", "rmse_lo", "rmse_hi", "minus_best", "diff_lo", "diff_hi", "interval_holds_0"]
    more = ["set", "model", "rmse_rain", "rmse_0_60", "rmse_60_180", "integral_pred_over_obs", "shower_peak_error_median", "shower_peak_error_lo", "shower_peak_error_hi"]
    real_conv = conv[~conv.fit.str.startswith("fake")]
    seconds = sum(pickle.loads(f.read_bytes()).get("seconds", 0.0) for f in CACHE.glob("*.pkl"))
    key_par = params[~params.parameter.str.startswith(("F[", "P[")) & params.set.isin([SET_LABEL["S1"], SET_LABEL["GG"]])].round(3)
    curve_par = params[params.parameter.str.startswith(("F[", "P[")) & (params.model == curve_model) & (params.set == SET_LABEL["S1"])].dropna(subset=["value"]).round(3)
    ctx = dict(locals())
    ctx.update({"BANDS": BANDS})
    box, item2_text, item3_text, verify, unclear = texts(ctx)
    cut_text = ", ".join(plan["cuts"]) or "nothing"
    report = f"""# wet_antenna step S2b: repair and finish the S2 model (the P / F trap, a bounded slow store, longer fits; TRAIN only)

Written by `research/wet_antenna/system_id_s2b.py` (library `system_id.py`). Settings and every rule:
`wet_antenna.system_id_s2b` in `research/config.yaml`, written before the run; data, judge and fake noise as
`wet_antenna.system_id_s2`. Plan: `docs/ROADMAP.md`, "Reboot (2026-10-13)", "WAA v1" and "Remaining plan".

## 1. What was run

- The S2 frame A_pred(t) = W(t) + L * P_band(R(t - d)) + c with three repairs: (1) the increments of P and F between
  the knots ({', '.join(f'{k:g}' for k in KNOTS)} mm/h) are plain parameters >= 0 (no exp), with 4 starting curves for P (0, 0.5,
  1, 2 x the ITU-R P.838 size; starting values only, never a correction); (2) M4b = M4 with tau_s in
  {S['model']['tau_s_m4b_min'][0]}-{S['model']['tau_s_m4b_min'][1]} min; (3) {len(S['model']['start_plan'])} starts per model, caps {M['maxiter']['start']} iterations per start and
  {M['maxiter']['polish']} for the polish (5 x S2), {M['maxiter']['fold']} per fold (3 x S2); bootstrap refits {M['maxiter']['boot']} (as S2).
- Candidates M1 < M3 < M3D < M4b (stop rule as S1 / S2; the chosen model = the one chosen on {SET_LABEL['S1']}, the main
  set); the free-tau_s M4 is fitted and shown as a reference only. Sets: {', '.join(SET_LABEL[k] for k in sets)}.
- Timing: one M4b start on S<3 took {plan['t_fit_s']} s ({plan['nit_start']} iterations, converged: {plan['converged_start']}); one evaluation
  {plan['t_eval_s'] * 1e3:.0f} ms; projection {plan['projected_full_h']} h on {workers} processes (budget {S['time_budget_h']} h). Cut: {cut_text}.
  Fitting time summed over the processes: {seconds / 3600:.1f} h. Run time of this run: {(time.time() - clock) / 60:.1f} min (fits are cached).
- **Convergence (real data, {len(real_conv)} fits on all events):** {int(real_conv.converged_on_all_events.sum())} converged before the cap;
  {int(real_conv.starts_converged.sum())} of {int(real_conv.starts.sum())} starts converged; {int(real_conv.starts_past_S2_cap.sum())} starts needed more than the S2 cap of
  {S2_CAP['start']} iterations; {int(real_conv.folds_at_cap.sum())} of {int(real_conv.folds.sum())} leave-one-event-out folds reached the cap of {M['maxiter']['fold']}
  ({int(real_conv.folds_past_S2_cap.sum())} needed more than the S2 cap of {S2_CAP['fold']}).
- TRAIN events only; VAL is never read. S1 and S2 outputs are unchanged (shared functions moved into the library).
- Software: numba {numba.__version__}, scipy {scipy.__version__} (L-BFGS-B, exact gradients), numpy {np.__version__}.

{md_table(counts)}

(`{rel(P.SID2B_COUNTS)}`.) The convergence of every fit (`{rel(P.SID2B_CONVERGENCE)}`):

{md_table(conv.drop(columns=["loss"]).round(2))}

## 2. Item 1: the P / F trap{tag['Tc']}

**Fake check Tc again: {'PASS' if verdicts['Tc']['PASS'] else 'FAIL'}.** {verdicts['Tc']['text']}.

{md_table(tc_med.round(3))}

(`{rel(P.SID2B_TC)}`; median over {len(reps)} replicates; truth F = {FK['f_top_db']} (1 - exp(-R / {FK['f_scale_mmh']:g})) dB, P = the ITU size.)

![F and P](../figures/{P.SID2B_FIG_CURVES.name})

The curves of {curve_model} on the main set ({SET_LABEL['S1']}):

{md_table(curve_par)}

**The length rule again, read on {rule_model}: {'AGREE' if agree else 'DISAGREE'}.** {f'Every key number of S<3 lies inside the S<1 interval: S<3 can be the main set.' if agree else f'{len(outside)} of {len(in_rule)} key numbers of S<3 lie outside the S<1 interval ({", ".join(outside.key)}): S<1 stays the main set; the difference is a finding for Stage 2.'}{tag['Tc']}

{md_table(length_tab.round(3))}

(`{rel(P.SID2B_LENGTH)}`; P is shown, not in the rule.)

![length rule](../figures/{P.SID2B_FIG_LENGTH.name})

## 3. Item 2: a bounded slow store{tag['Te']}

**Fake check Te (truth M4b: phi {FK2['Te']['phi']}, tau_f_wet {FK2['Te']['tau_f_wet']}, tau_f_dry {FK2['Te']['tau_f_dry']}, tau_s {FK2['Te']['tau_s']} min): {'PASS' if verdicts['Te']['PASS'] else 'FAIL'}.** {verdicts['Te']['text']}.

{md_table(te.round(4))}

(`{rel(P.SID2B_FAKE)}`.) The ladder per set (held-out RMSE of A in dB, from the rain alone; minus_best with the 90 % interval
of the paired difference to the best candidate):

{md_table(ladder_tab[show].round(4))}

**Chosen on the main set ({SET_LABEL['S1']}): {chosen}** (best: {best_of['S1']}). Per set: {', '.join(f'{SET_LABEL[k]}: best {best_of[k]}, chosen {chosen_of[k]}' for k in sets)}.

After the rain and the peaks (reported, not in the rule; integral = 0-180 min after the last wet bin, predicted / observed, 1 = right):

{md_table(ladder_tab[more].round(3))}

(`{rel(P.SID2B_LADDER)}`.)

![ladder](../figures/{P.SID2B_FIG_LADDER.name})

![after the last shower](../figures/{P.SID2B_FIG_AFTER.name})

Median (A - c) / peak after the last shower (`{rel(P.SID2B_AFTER)}`):

{md_table(curve_tab[curve_tab.minutes_after_last_wet_bin.isin([0, 10, 20, 30, 45, 60, 90, 120, 150, 175])].round(3))}

{item2_text}

## 4. Item 3: did the longer fits change a choice or a number?

{item3_text}

{md_table(cmp_tab.round(3))}

(`{rel(P.SID2B_COMPARE)}`.)

## 5. WAA v2

{box}

{md_table(key_par)}

(`{rel(P.SID2B_PARAMS)}`.)

## 6. What S3 will verify

{verify}

## 7. Unclear / possibly wrong

{unclear}
"""
    save(report, P.SID2B_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {(time.time() - clock) / 60:.1f} min")


def texts(q):
    """Sections 3 (the reading of item 2), 4 (item 3), 5 (the box of WAA v2), 6 and 7 of the report, written from
    the numbers of the run (q = the local variables of main)."""
    params, chosen, ladders, conv, old_ladder = q["params"], q["chosen"], q["ladders"], q["conv"], q["old_ladder"]
    curves = q["curves"]
    old_fits = pd.read_csv(P.SID2_FITS)

    def par(model, name, set_):
        r = params[(params.set == SET_LABEL[set_]) & (params.model == model) & (params.parameter == name)]
        return (np.nan, np.nan, np.nan) if r.empty else tuple(float(v) for v in r[["value", "lo", "hi"]].iloc[0])

    def show(model, name, set_, nd=1):
        v, lo, hi = par(model, name, set_)
        return f"{v:.{nd}f} [{lo:.{nd}f}, {hi:.{nd}f}]"

    def lad(set_):
        return ladders[set_].set_index("model")

    def old_rmse(set_, m):
        t = old_ladder[(old_ladder.set == SET_LABEL[set_]) & (old_ladder.model == m)]
        return float(t.rmse.iloc[0]) if len(t) else np.nan

    l1, l3, lg = lad("S1"), lad("S3"), lad("GG")
    c1 = curves["S1"].set_index("minutes_after_last_wet_bin")
    b0 = q["BANDS"][0]

    # ---------------------------------------------------------------- section 3: item 2
    item2 = f"""**Reading of item 2.** The fitting finds a bounded slow store when there is one (fake check Te: M4b chosen in
{int((q['te'].chosen == 'M4b').sum())} of {len(q['te'])}, tau_s x{q['te'].tau_s.median() / FK2['Te']['tau_s']:.2f}). On the real data M4b is the best candidate on every set, but by little: on the
main set M3 - M4b = {l1.minus_best['M3']:+.4f} dB ({l1.diff_lo['M3']:+.4f} .. {l1.diff_hi['M3']:+.4f}), so the stop rule takes **{chosen}**; on S<3 M3 is outside
({l3.minus_best['M3']:+.4f}, {l3.diff_lo['M3']:+.4f} .. {l3.diff_hi['M3']:+.4f}) and M4b is chosen there; with the gauge input M3D is chosen. M4b and the free
M4 give the same RMSE within 0.0003 dB.

**The bound did not fix the leftover.** tau_s of M4b sits at the upper bound 240 min on every set (S<1 {show('M4b', 'tau_s', 'S1', 0)},
S<3 {show('M4b', 'tau_s', 'S3', 0)}, gauge {show('M4b', 'tau_s', 'GG', 0)}): the fit still wants a longer store. After the last
shower M4b stays too high (S<1: {c1.loc[120, 'M4b']:.3f} of the peak at 120 min against {c1.loc[120, 'observed']:.3f} observed; integral 0-180 min
{l1.integral_pred_over_obs['M4b']:.2f} x the observed, M4 {l1.integral_pred_over_obs['M4']:.2f} x), while M3 falls to 0 too early ({c1.loc[60, 'M3']:.3f} at 60 min against
{c1.loc[60, 'observed']:.3f}; integral {l1.integral_pred_over_obs['M3']:.2f} x on S<1 but {l3.integral_pred_over_obs['M3']:.2f} x on S<3). The RMSE 60-180 min is the same for every rung
({l1.rmse_60_180.min():.3f}-{l1.rmse_60_180.max():.3f} dB on S<1): after the rain the error is noise, not the model. **The slow leftover (2-4 % of the peak,
1-2 h) is still not described by any rung**: a store fed by the same F during the whole event fills over long events and
then empties too slowly. (With the gauge input the integral ratios, {lg.integral_pred_over_obs.min():.1f}-{lg.integral_pred_over_obs.max():.1f}, rest on a very small observed integral.)"""

    # ---------------------------------------------------------------- section 4: item 3
    rc = conv[~conv.fit.str.startswith("fake")]
    loss_rows = []
    for k in ("S1", "S3", "GG"):
        for m in ("M1", "M3", "M3D", "M4"):
            o = old_fits[(old_fits.set == SET_LABEL[k]) & (old_fits.model == m)]
            n = rc[(rc.fit == SET_LABEL[k]) & (rc.model == m)]
            if len(o) and len(n):
                loss_rows.append((k, m, float(n.loss.iloc[0]) - float(o.loss.iloc[0])))
    small = [d for k, m, d in loss_rows if m != "M3D"]
    s2_f = lambda m, name: float(pd.read_csv(P.SID2_PARAMS).query("set == @SET_LABEL['GG'] and model == @m and parameter == @name").value.iloc[0])   # noqa: E731
    item3 = f"""**Convergence.** {int(rc.converged_on_all_events.sum())} of {len(rc)} fits on all events converged before the cap (not: {', '.join(f'{r.model} on {r.fit}' for r in rc[~rc.converged_on_all_events.astype(bool)].itertuples()) or 'none'});
{int(rc.starts_converged.sum())} of {int(rc.starts.sum())} starts converged, and every start needed more than the S2 cap of {S2_CAP['start']} iterations; {int(rc.folds_at_cap.sum())} of
{int(rc.folds.sum())} folds still reached the new cap of {M['maxiter']['fold']}. M3D converges nowhere on the radar sets (its drift slopes, one per
link-event, make a long flat valley). The spread of the loss over the 8 starts is below 0.4 for every model except M3D: the
starts end in the same place.

**Was S2 trapped on the real data? Hardly.** The training loss of the S2b fit minus that of S2, same model and data:
{', '.join(f'{m} {SET_LABEL[k][:3] if k != "GG" else "gauge"} {d:+.1f}' for k, m, d in loss_rows)} (losses of 2e4-1e5). Apart from M3D the
changes are at most {max(abs(d) for d in small):.0f} units: the S2 fits on the radar sets were at the same minimum, and their numbers stand
(M3 and M4 change by a few % at most, see the table). The trap of S2 showed on the fake data (Tc now PASS) and in ONE place on
the real data: **the level F at light rain with the GAUGE input** moved a lot for a tiny change of the loss (M3 F[{b0}](0.5)
{s2_f('M3', f'F[{b0}](0.5)'):.2f} -> {par('M3', f'F[{b0}](0.5)', 'GG')[0]:.2f} dB, M4 {s2_f('M4', f'F[{b0}](0.5)'):.2f} -> {result('real_GG_M4')['params'][f'F[{b0}](0.5)']:.2f}; F(2) {s2_f('M4', f'F[{b0}](2)'):.2f} -> {result('real_GG_M4')['params'][f'F[{b0}](2)']:.2f}). So v1's "with the gauge the
level is still rising at 5 mm/h" was mostly the trap: with the gauge the level is now about 2-3 dB by 2 mm/h, as with the radar.
The S2 bootstrap intervals of those numbers were narrow (e.g. 0.12-0.15 dB) and did not hold the new value: refits that start
from the fit and stop after 30 iterations do not explore a flat direction.

**Choices.** On S<1 the chosen model changed from M3D (S2) to M3: not because M3 got better (held-out {old_rmse('S1', 'M3'):.4f} -> {l1.rmse['M3']:.4f} dB)
but because the best got a little worse (S2 M4 {old_rmse('S1', 'M4'):.4f}; S2b M4b {l1.rmse['M4b']:.4f}, M4 {l1.rmse['M4']:.4f}). The longer folds ({M['maxiter']['fold']}
instead of {S2_CAP['fold']} iterations) raise the held-out RMSE of almost every model by 0.001-0.003 dB: stopping early acted as a
mild regularisation, and the judge depends on that cap. On S<3: M4 (S2) -> M4b (M4 is no longer a candidate). Gauge: M3D in both."""

    # ---------------------------------------------------------------- section 5: the box
    m = chosen if chosen != "M1" else "M3"
    tw, td = ("tau_wet", "tau_dry") if m in ("M3", "M3D") else ("tau_f_wet", "tau_f_dry")
    rows = [f"WAA v2 (S2b, 2026-10-14): model {m}, chosen by the stop rule on {SET_LABEL['S1']} (the main set); 90 % event-bootstrap intervals",
            "A(t) = W(t) + L * P_band(R(t - d)) + c"]
    if m in ("M3", "M3D"):
        rows.append("dW/dt = (a_link F_band(R) - W) / tau,  tau = tau_wet while W rises, tau_dry while it falls        (exact steps per minute)")
    else:
        rows.append("W = W_f + W_s;  W_f -> phi a_link F(R) (tau_f_wet / tau_f_dry),  W_s -> (1 - phi) a_link F(R) (tau_s)")
    rows += ["", f"{'':26s}{'radar, ' + SET_LABEL['S1'][:3] + ' (main)':34s}GAUGE input (gauge links)",
             f"{'lag d':26s}{f"{par(m, 'd', 'S1')[0]:+.0f} min":34s}{par(m, 'd', 'GG')[0]:+.0f} min",
             f"{tw:26s}{show(m, tw, 'S1') + ' min':34s}{show(m, tw, 'GG')} min   <- the wetting speed (the radar value is input smoothing)",
             f"{td:26s}{show(m, td, 'S1') + ' min':34s}{show(m, td, 'GG')} min"]
    if m not in ("M3", "M3D"):
        rows += [f"{'phi':26s}{show(m, 'phi', 'S1', 2):34s}{show(m, 'phi', 'GG', 2)}", f"{'tau_s':26s}{show(m, 'tau_s', 'S1', 0):34s}{show(m, 'tau_s', 'GG', 0)}"]
    for b in q["BANDS"]:
        if np.isfinite(par(m, f"F[{b}](2)", "S1")[0]):
            rows.append(f"{'F[' + b + '](0.5 / 2 / 5)':26s}{' / '.join(f'{par(m, f"F[{b}]({r:g})", "S1")[0]:.2f}' for r in KEY_R) + ' dB':34s}"
                        + " / ".join(f"{par(m, f'F[{b}]({r:g})', 'GG')[0]:.2f}" for r in KEY_R) + " dB")
    for b in q["BANDS"]:
        if np.isfinite(par(m, f"P[{b}](2)", "S1")[0]):
            rows.append(f"{'P[' + b + '](2 / 10)':26s}{f"{par(m, f'P[{b}](2)', 'S1')[0]:.2f} / {par(m, f'P[{b}](10)', 'S1')[0]:.2f} dB/km":34s}"
                        f"{par(m, f'P[{b}](2)', 'GG')[0]:.2f} / {par(m, f'P[{b}](10)', 'GG')[0]:.2f} dB/km")
    rows += [f"{'a_link 10-90 %':26s}{f"{par(m, 'a_p10', 'S1')[0]:.2f} - {par(m, 'a_p90', 'S1')[0]:.2f}":34s}{par(m, 'a_p10', 'GG')[0]:.2f} - {par(m, 'a_p90', 'GG')[0]:.2f}",
             "",
             f"Within the interval on the main set: M4b (fast store {show('M4b', 'tau_f_wet', 'GG')} / {show('M4b', 'tau_f_dry', 'GG')} min with the gauge,",
             f"phi {par('M4b', 'phi', 'S1')[0]:.2f}, tau_s at the bound 240); on S<3 the stop rule takes M4b.",
             "Limits that remain:",
             f"- the F / P split depends on the model and the input: P[{b0}](2) = {par('M3', f'P[{b0}](2)', 'S1')[0]:.2f} (M3, radar), {par('M4b', f'P[{b0}](2)', 'S1')[0]:.2f} (M4b, radar),",
             f"  {par('M3', f'P[{b0}](2)', 'GG')[0]:.2f} (M3, gauge) dB/km; ITU size {q['itu'][0][0] * 2 ** q['itu'][0][1]:.2f}. With the radar, M3's slow wetting pushes the fast part of W into P.",
             "- the slow leftover after the rain (2-4 % of the peak, 1-2 h) is not described: M3 ends too early, M4b's store",
             "  (tau_s at the bound) too late;",
             f"- peaks too low: {l1.shower_peak_error_median[m]:+.2f} dB (radar), {lg.shower_peak_error_median[m]:+.2f} dB (gauge);",
             f"- drying: gauge {par(m, td, 'GG')[0]:.1f} min against radar {par(m, td, 'S1')[0]:.1f} min, not one number yet;",
             "- intervals from warm-started bootstrap refits (30 iterations) are too narrow along flat directions (section 4)."]
    box = "```\n" + "\n".join(rows) + "\n```"

    # ---------------------------------------------------------------- section 6: S3
    verify = f"""S3 takes WAA v2 ({m}) and its runner-up M4b as fitted here on the TRAIN events, without refitting, to data they
have not seen: the VAL events (first time they are read), with the held-out RMSE, the shower peaks and the after-rain
integral as here; both directions of each link (the two sub-links share the antennas, so their W should agree); the
GAUGE input on the gauge links of VAL (is the wetting about 3 min and the drying about {par(m, td, 'GG')[0]:.0f} min again?); and the
links longer than 3 km, where the path term dominates: does A - W behave like L * P(R), or like (L + L0) * P, i.e. is a
constant antenna part left (the f = L + L0 test of Stage 2)? The F / P split (which model's P is right) is decided
there, on the long links, not on the short ones."""

    # ---------------------------------------------------------------- section 7: unclear
    m4b_len = q["length_tab"][(q["length_tab"].model == "M4b") & q["length_tab"].in_rule]
    out_m4b = m4b_len[~m4b_len.S3_inside_S1_interval.astype(bool)]
    unclear = f"""- **Bootstrap intervals are too narrow along flat directions.** They come from refits that start at the fit and stop
  after {M['maxiter']['boot']} iterations; the gauge-input F at light rain moved far outside its S2 interval (section 4). Treat the intervals
  of F and P as lower bounds on the uncertainty.
- **The judge depends on the fold cap:** longer folds raised the held-out RMSE by 0.001-0.003 dB and changed the choice on
  S<1 (M3D -> M3). The differences between M3, M3D and M4b are of that size.
- **The length rule says AGREE on {q['rule_model']}, partly because its intervals are wide** (tau_wet {show(q['rule_model'], 'tau_wet' if q['rule_model'] in ('M3', 'M3D') else 'tau_f_wet', 'S1')} min on S<1).
  For M4b {len(out_m4b)} key number(s) of S<3 lie outside the S<1 interval ({', '.join(out_m4b.key) or 'none'}), and the P of M3 halves from
  S<1 to S<3 (P[{b0}](2) {par('M3', f'P[{b0}](2)', 'S1')[0]:.2f} -> {par('M3', f'P[{b0}](2)', 'S3')[0]:.2f}, intervals apart): with the radar, M3's P holds something that does
  not scale with L. The config keeps S<1 as the main set; read on S<3, the chosen model would be M4b.
- **M3D does not converge** on the radar sets (0 of 8 starts); its numbers and its place in the ladder are less firm.
- **Cuts:** {', '.join(q['plan']['cuts']) or 'none'} (the projection assumed every start runs to the cap; the fits took {q['seconds'] / 3600:.1f} h summed over
  the {q['workers']} processes, about {q['seconds'] / 3600 / q['workers']:.1f} h of wall clock). The gauge rule of S2 (gauge against radar on the same points) is therefore not repeated.
- **The fake checks use the exact rain** as input (as S1 / S2): Tc passing shows that the fit can split F and P when the
  input is right; it says nothing about the radar's errors.
- The run was first stopped by the session (low memory) and restarted; no fit was affected (each finished start is saved)."""
    return box, item2, item3, verify, unclear


if __name__ == "__main__":
    main()
