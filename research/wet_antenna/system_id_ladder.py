"""wet_antenna step S1: system identification of the wet-antenna attenuation on short links, the model ladder M0-M4.

The reboot of 2026-10-13 (docs/ROADMAP.md): rain R(t) in, attenuation A(t) out, on the links below 1 km, on whole
events (all showers, no tail cutting). Every model predicts

    A_pred(t) = W(t) + L * P_band(R(t - d)) + c

and the ladder adds one thing per rung: M0 no WAA, M1 static W = a F(R), M2 one memory time constant, M3 separate
wetting and drying speeds, M4 two stores (fast + slow). The models and their fit are in the library system_id.py.
RULE: WAA is never estimated by subtracting a power-law rain term. The path term P is a free curve that is told
apart from W only by the path length L. ITU-R P.838 is used here only to BUILD the fake data and as a size reference.
What the script does (TRAIN events only; VAL is never read):
  data     the link-events of the short links (window from 2 h before the rain on the link to 3 h after it), A as
           5 min means, the radar rain on the path as input, the offset c from the dry minutes before the rain;
  fake     the fake check first: truths M3 and M4 driven by the real rain and windows, with noise; the same ladder
           and judge must choose the true rung and recover its numbers;
  ladder   every model fitted on all events and judged by leave-one-event-out prediction; the stop rule chooses the
           simplest model that is not worse than the best;
  after    parameters of the chosen model with event-bootstrap intervals, five sensitivity runs, three checks
           against earlier facts, figures, the product per link-event, the report.
Every rule is in config `wet_antenna.system_id`, written before the run. One fit is timed first and the projection
decides the cuts. Fits are cached in research/.cache/wet_antenna/system_id_ladder/ (delete the folder to fit again).
Run from the repository root (about 2.5 h the first time on 3 processes, under 1 min from the cache):
    python -m research.wet_antenna.system_id_ladder
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "1")            # one thread per process: the processes run in parallel

from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import pickle
import sys
import time
import warnings

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
from research.wet_antenna.simulator import ar1, quantize
from research.wet_antenna.tails import itu_p838

CFG = load_config()
S = CFG["wet_antenna"]["system_id"]
DATA, M, J, FK, SENS, CHK = S["data"], S["model"], S["judge"], S["fake"], S["sensitivity"], S["checks"]
CACHE = P.CACHE / "system_id_ladder"
BANDS = list(DATA["bands_ghz"])
KNOTS = np.asarray(M["knots_mmh"], float)
PADM = sid.PADM
LINE = {"M0": ":", "M1": "--", "M2": "-.", "M3": "-", "M4": (0, (5, 1, 1, 1, 1, 1))}


# =================================================================================================== the data
# (the functions are in the library system_id.py, shared with step S2; here they get this step's settings)
def build_dataset(d):
    return sid.build_dataset(d, DATA)


def finish(D):
    return sid.finish(D, BANDS)


def subset(D, keep):
    return sid.subset(D, keep)


def save_data(D, name):
    return sid.save_data(D, CACHE / f"data_{name}.npz")


def load_data(file):
    return sid.load_data(file, BANDS)


# =================================================================================================== the fake data
def make_fake(D, le, truth, rep, noise):
    """A fake dataset: the truth `truth` (M3 or M4, config `fake`) driven by the real radar rain, windows, lengths and
    missing minutes, observed with noise, quantization and an unknown offset; c by the same rule as the real data."""
    rng = np.random.default_rng([FK["seed_base"], sid.MODELS.index(truth), rep])
    t = FK["truths"][truth]
    a = np.exp(rng.normal(0, FK["a_sd"], D["nl"]))
    freq = le.groupby("band").frequency_ghz.median()
    itu = {b: itu_p838(freq[b], "v") for b in freq.index}
    y1, c = np.full(len(D["y1"]), np.nan), np.zeros(len(le))
    f_var = np.sqrt(FK["noise_variance_factor"])
    for e, r in enumerate(le.itertuples()):
        n, o1 = D["le_n"][e], D["le_boff1"][e]
        rain = D["rain"][D["le_roff"][e] + PADM - FK["lag_min"]:D["le_roff"][e] + PADM - FK["lag_min"] + n]
        u = a[D["le_link"][e]] * FK["f_top_db"] * (1 - np.exp(-rain / FK["f_scale_mmh"]))
        w = sid.simulate_stores(u, t.get("phi", 1.0), t.get("tau_f_wet", t.get("tau_wet")), t.get("tau_f_dry", t.get("tau_dry")), t.get("tau_s", 1.0))
        truth_a = w + r.L_km * itu[r.band][0] * np.power(rain, itu[r.band][1])
        sw, ss, tn = noise.get(r.sublink_id, noise["median"])
        off = rng.standard_normal() * (r.offset_std if np.isfinite(r.offset_std) else 0.0)
        slow = np.repeat(ar1(rng.standard_normal(n), ss * f_var, tn, rng.standard_normal() * ss * f_var), 6) if ss > 0 else 0.0
        x = np.repeat(truth_a, 6) + off + slow + sw * f_var * rng.standard_normal(n * 6)
        obs = quantize(x, rng.uniform(0, FK["rsl_step_db"]), FK["rsl_step_db"]).reshape(n, 6).mean(axis=1)
        obs[D["w1"][o1:o1 + n] == 0] = np.nan                      # the real missing minutes
        y1[o1:o1 + n] = obs
        pre = D["pre1"][o1:o1 + n] & np.isfinite(obs)
        c[e] = float(np.median(obs[pre])) if pre.sum() >= DATA["offset_min_minutes"] else 0.0
    F = {k: v for k, v in D.items()}
    F["y1"], F["le_c"] = y1, c
    m5 = y1.reshape(-1, 5)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        F["y5"] = np.where(np.isfinite(m5).sum(axis=1) >= 3, np.nanmean(m5, axis=1), np.nan)
    return finish(F), {"a_log_sd": float(np.log(a).std())}


# =================================================================================================== the jobs (worker processes)
def run_job(job):
    """One model on one dataset: the fit on all events, the leave-one-event-out predictions, and (if asked) the
    event-bootstrap refits. Every stage is saved in the cache as soon as it is finished."""
    file = CACHE / f"{job['key']}.pkl"
    res = pickle.loads(file.read_bytes()) if file.exists() else {}
    D = load_data(job["data"])
    opt = {"no_path": job.get("no_path", False), "tau_per_band": job.get("tau_per_band", False)}
    Pm = sid.Problem(D, job["model"], M, agg=job.get("agg", 5), **opt)
    P5 = Pm if Pm.agg == 5 else sid.Problem(D, job["model"], M, agg=5, **opt)      # the judge always works on 5 min means
    t0 = time.time()

    def dump():
        res["seconds"] = res.get("seconds", 0.0) + time.time() - t0
        file.write_bytes(pickle.dumps(res))

    if "all" not in res:
        res["all"] = sid.fit_all(Pm, M["starts"], M["lag_grid_min"], M["maxiter"])
        res["pred_all"], res["params"] = P5.predict(res["all"]["x"], res["all"]["d"]), Pm.params(res["all"]["x"], res["all"]["d"])
        dump()
    if job.get("loo") and "loo_pred" not in res:
        pred, folds = np.full(len(D["y5"]), np.nan), []
        for e in np.unique(D["le_event"]):
            fit = sid.fit_weighted(Pm, res["all"], (D["le_event"] != e).astype(float), M["lag_candidates"], M["maxiter"]["fold"])
            held = D["le_event"][D["bin_le"]] == e
            pred[held] = P5.predict(fit["x"], fit["d"])[held]
            folds.append({"event": int(e), "d": fit["d"], "nit": fit["nit"], **{k: v for k, v in Pm.params(fit["x"], fit["d"]).items() if k[0] in "tp"}})
        res["loo_pred"], res["folds"] = pred, folds
        dump()
    if job.get("boot") and "boot" not in res:
        rng, n_ev, draws = np.random.default_rng([J["seed"], sid.MODELS.index(job["model"])]), int(D["le_event"].max()) + 1, []
        for _ in range(job["boot"]):
            mult = np.bincount(rng.integers(0, n_ev, n_ev), minlength=n_ev).astype(float)
            fit = sid.fit_weighted(Pm, res["all"], mult[D["le_event"]], M["lag_candidates"], M["maxiter"]["boot"])
            draws.append(Pm.params(fit["x"], fit["d"]))
        res["boot"] = draws
        dump()
    return job["key"]


def run_pool(jobs, workers, label):
    jobs = [j for j in jobs if not done(j)]
    if not jobs:
        return
    t0 = time.time()
    with ProcessPoolExecutor(workers) as ex:
        for i, f in enumerate(as_completed([ex.submit(run_job, j) for j in jobs]), 1):
            print(f"   {label}: {f.result()} done ({i} of {len(jobs)}, {(time.time() - t0) / 60:.1f} min)", flush=True)


def done(job):
    file = CACHE / f"{job['key']}.pkl"
    if not file.exists():
        return False
    res = pickle.loads(file.read_bytes())
    return "all" in res and (not job.get("loo") or "loo_pred" in res) and (not job.get("boot") or "boot" in res)


def result(key):
    return pickle.loads((CACHE / f"{key}.pkl").read_bytes())


# =================================================================================================== the judge
def event_errors(D, pred, w=None):
    return sid.event_errors(D, pred, J, w)


def link_event_table(D, pred, w=None):
    return sid.link_event_table(D, pred, J, w)


def judge(errors, order):
    return sid.judge(errors, order, J)


def ladder(name, D, models):
    """The judge on the cached leave-one-event-out predictions of one dataset."""
    errors = {m: event_errors(D, result(f"{name}_{m}")["loo_pred"]) for m in models}
    return (*judge(errors, [m for m in J["order"] if m in models]), errors)


# =================================================================================================== the checks
def showers(D, e):
    return sid.showers(D, e, CHK["shower_gap_min"])


def checks(D, le, preds):
    """Three descriptive checks on held-out predictions (preds = {label: predicted A per 5 min bin})."""
    n_after = J["after_windows_min"][-1][1] // 5
    curves = {k: [] for k in ["observed"] + list(preds)}
    shower_rows, drop_rows = [], []
    f0, f1 = CHK["fast_window_min"][0] // 5, CHK["fast_window_min"][1] // 5
    for e in range(len(le)):
        o, n, c = D["le_boff5"][e], D["le_n"][e] // 5, D["le_c"][e]
        w = D["w5"][o:o + n] > 0
        obs = np.where(w, D["y5"][o:o + n] - c, np.nan)
        pr = {k: np.where(w, p[o:o + n] - c, np.nan) for k, p in preds.items()}
        sh = showers(D, e)
        for i, (s0, s1) in enumerate(sh):
            if np.isfinite(obs[s0:s1 + 1]).any():
                shower_rows.append({"le": e, "shower": i, "n_showers": len(sh), "minutes": (s1 - s0 + 1) * 5, "peak_R_mmh": float(np.nanmax(np.r_[D["r5"][o + s0:o + s1 + 1], 0.0])),
                                    "peak_obs": float(np.nanmax(obs[s0:s1 + 1])), **{f"peak_{k}": float(np.nanmax(p[s0:s1 + 1])) for k, p in pr.items()}})
        last = D["le_last"][e]
        if sh and np.isfinite(obs[sh[-1][0]:last + 1]).any():
            peak = float(np.nanmax(obs[sh[-1][0]:last + 1]))
            if peak >= CHK["min_peak_db"]:
                pad = lambda v: np.r_[v[last + 1:last + 1 + n_after], np.full(n_after, np.nan)][:n_after] / peak   # noqa: E731
                curves["observed"].append(pad(obs))
                for k, p in pr.items():
                    curves[k].append(pad(p))
        a = obs[last + 1:last + 1 + f1]
        if len(a) >= f1 and np.isfinite(a[0]) and np.isfinite(a[f0:f1]).any() and a[0] >= CHK["fast_min_start_db"]:
            drop_rows.append({"le": e, "L_km": le.L_km[e], "start_obs": a[0], "drop_obs": a[0] - np.nanmean(a[f0:f1]),
                              **{f"drop_{k}": p[last + 1] - np.nanmean(p[last + 1 + f0:last + 1 + f1]) for k, p in pr.items()}})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        curve = pd.DataFrame({"minutes_after_last_wet_bin": np.arange(n_after) * 5, "link_events": len(curves["observed"]),
                              **{k: np.nanmedian(np.array(v), axis=0) for k, v in curves.items()}})
    return curve, pd.DataFrame(shower_rows), pd.DataFrame(drop_rows)


# =================================================================================================== figures
def figure_ladder(tab, chosen, best):
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    x = np.arange(len(tab))
    ax.errorbar(x, tab.rmse, yerr=[tab.rmse - tab.rmse_lo, tab.rmse_hi - tab.rmse], fmt="o", color=INK, ms=7, capsize=5, lw=1.5)
    for i, r in enumerate(tab.itertuples()):
        note = "  CHOSEN" if r.model == chosen else ""
        note += "  (best)" if r.model == best else ""
        ax.text(i + 0.08, r.rmse, f"{r.rmse:.3f}{note}", fontsize=8.5, va="center", fontweight="bold" if r.model == chosen else "normal")
        if r.model == chosen:
            ax.plot(i, r.rmse, "s", ms=16, mfc="none", mec=INK, mew=1.6)
    ax.set(xticks=x, xticklabels=tab.model, xlim=(-0.5, len(tab) - 0.1), xlabel="model (simple to complex)", ylabel="held-out RMSE of A (dB), mean over events")
    style_axes(ax)
    ax.set_title("Step S1: leave-one-event-out error of the model ladder (dot = mean over the held-out events,\nwhisker = 90 % event-bootstrap interval; "
                 "square = the model chosen by the stop rule)", fontsize=9.5)
    fig.tight_layout()
    return fig


def figure_gallery(D, le, picks, preds, chosen):
    fig = plt.figure(figsize=(17, 9.6))
    gs = fig.add_gridspec(4, 4, height_ratios=[3, 1.1, 3, 1.1], hspace=0.55, wspace=0.22, top=0.91, bottom=0.06, left=0.045, right=0.985)
    for j, e in enumerate(picks):
        ax, bx = fig.add_subplot(gs[2 * (j // 4), j % 4]), fig.add_subplot(gs[2 * (j // 4) + 1, j % 4])
        o, n, r = D["le_boff5"][e], D["le_n"][e] // 5, le.iloc[e]
        x = (np.arange(n) - r.last_bin - 1) * 5
        w = D["w5"][o:o + n] > 0
        ax.plot(x[w], D["y5"][o:o + n][w], "o", ms=2.6, color="#8a8a8a")
        ax.plot(x, preds[chosen][o:o + n], "-", color=INK, lw=1.5)
        ax.plot(x, preds["M1"][o:o + n], "--", color=INK, lw=1.0)
        ax.axhline(r.offset_c_db, color=INK2, lw=0.6)
        bx.step(x, np.nan_to_num(D["r5"][o:o + n]), where="post", color=INK, lw=1.2)
        for a_ in (ax, bx):
            a_.axvline(0, color=INK, lw=0.8)
            a_.set_xlim(x[0], x[-1])
            style_axes(a_)
        ax.set_title(f"{r.sublink_id}, event {r.event_id[1:]}\nL = {r.L_km:.2f} km, {r.band_name}", fontsize=8)
        if j % 4 == 0:
            ax.set_ylabel("A (dB)"), bx.set_ylabel("mm/h")
        if j // 4 == 1:
            bx.set_xlabel("minutes since the last wet bin of the link")
    fig.suptitle(f"Step S1: 8 held-out link-events (seed {S['gallery']['seed']}, observed peak >= {CHK['min_peak_db']} dB). Top: observed A (dots, 5 min means), the chosen model "
                 f"{chosen} (solid) and the static model M1 (dashed),\nboth predicted from the rain alone with the event left out of the fit; thin line = the offset c. "
                 "Bottom: the radar rain on the path. Vertical line: the end of the rain on the link.", fontsize=9.5)
    return fig


def figure_curves(params, boot, model):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.9))
    for ax, name, unit in zip(axes, ("F", "P"), ("dB", "dB / km")):
        for b, (band, ls, mk) in enumerate(zip(BANDS, ("-", "--"), ("o", "s"))):
            keys = [f"{name}[{band}]({k:g})" for k in KNOTS[1:]]
            if keys[0] not in params:
                continue
            v = np.r_[0.0, [params[k] for k in keys]]
            ax.plot(KNOTS, v, ls=ls, marker=mk, ms=5, mfc="none", color=INK, lw=1.7, label=f"{band} ({'solid, circles' if b == 0 else 'dashed, squares'})")
            if boot:
                q = np.percentile(np.array([[dr[k] for k in keys] for dr in boot]), J["interval"], axis=0)
                ax.vlines(KNOTS[1:] * (1.0 if b == 0 else 1.06), q[0], q[1], color=INK, lw=1.0)
            if name == "P":
                a_, b_ = itu_p838({"28 GHz": 28.7, "38 GHz": 38.0}[band], "v")
                ax.plot(KNOTS, a_ * KNOTS ** b_, ls=":", color=INK2, lw=1.0)
                ax.text(KNOTS[-1], a_ * KNOTS[-1] ** b_, f" ITU-R P.838, {band}\n (size reference only)", fontsize=7, va="center")
        ax.set(xscale="symlog", xlabel="rain R (mm/h)", ylabel=f"{name}(R) ({unit})", xlim=(0, 75 if name == "P" else 45))
        ax.set_xscale("symlog", linthresh=0.5)
        ax.set_xticks(KNOTS), ax.set_xticklabels([f"{k:g}" for k in KNOTS])
        ax.legend(frameon=False, fontsize=8.5, loc="upper left")
        style_axes(ax)
    axes[0].set_title("F: the WAA level the rain drives toward (for a_link = 1)", fontsize=9.5)
    axes[1].set_title("P: the path term per km, learned only through the path length", fontsize=9.5)
    fig.suptitle(f"Step S1: the two learned curves of the chosen model {model}, per band; vertical bars = 90 % event-bootstrap interval at each knot", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def figure_after(curve, chosen):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    for ax, logy in zip(axes, (False, True)):
        ax.plot(curve.minutes_after_last_wet_bin, curve.observed, "o", ms=4.5, color=INK, mfc="none", label="observed (circles)")
        ax.plot(curve.minutes_after_last_wet_bin, curve[chosen], "-", color=INK, lw=1.7, label=f"{chosen}, held out (solid)")
        ax.plot(curve.minutes_after_last_wet_bin, curve["M1"], "--", color=INK, lw=1.2, label="M1 static, held out (dashed)")
        ax.set(xlabel="minutes after the last wet bin of the link", title="the same, logarithmic" if logy else "median (A - c) / observed peak of the last shower")
        if logy:
            ax.set_yscale("log")
            ax.set_ylim(3e-3, 0.5)                         # values at 0 (M1 after the rain) fall below the axis
        ax.legend(frameon=False, fontsize=8.5)
        style_axes(ax)
    fig.suptitle(f"Step S1: after the last shower ({int(curve.link_events.iloc[0])} link-events with an observed peak >= {CHK['min_peak_db']} dB): observed against predicted", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    clock = time.time()
    workers = max((os.cpu_count() or 2) - 1, 1)

    # ------------------------------------------------------------------------------------------- data
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    D, le, drops = build_dataset(d)
    real_file = save_data(D, "real")
    n_ev = int(D["le_event"].max()) + 1
    w = D["w5"] > 0
    counts = pd.DataFrame([{"group": name, "links": g.cml_id.nunique(), "sub_links": g.sublink_id.nunique(), "link_events": len(g), "events": g.event_id.nunique(),
                            "points_5min": int(w[np.isin(D["bin_le"], g.index)].sum()), "wet_points": int((w & D["wet5"])[np.isin(D["bin_le"], g.index)].sum()),
                            "dry_points": int((w & ~D["wet5"])[np.isin(D["bin_le"], g.index)].sum()), "links_with_city_gauge_1km": g.cml_id[g.gauge_id != ""].nunique(),
                            "offset_flagged": int(g.offset_flag.sum()), "cut_at_tsl_step": int(g.tsl_cut.sum()), "median_window_h": round(g.n_bins.median() / 12, 1)}
                           for name, g in [("all", le)] + [(b, le[le.band_name == b]) for b in BANDS]
                           + [(f"< {SENS['a_max_length_km']} km", le[le.L_km < SENS["a_max_length_km"]])]])
    save(counts, P.SID_COUNTS, index=False)
    drops_tab = pd.DataFrame({"rule": list(drops), "pairs": list(drops.values())})
    print(f"TRAIN only: {len(le)} link-events, {le.cml_id.nunique()} links, {le.sublink_id.nunique()} sub-links, {n_ev} events, {int(w.sum())} points of 5 min; "
          f"{workers} processes; loaded in {time.time() - clock:.0f} s", flush=True)

    # ------------------------------------------------------------------------------------------- time one fit, plan
    plan_file = CACHE / "plan.json"
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
        print(f"plan from {rel(plan_file)}: cut: {plan['cuts'] or 'nothing'}")
    else:
        Pm = sid.Problem(D, "M3", M)
        ones = np.ones(len(le))
        Pm.run(Pm.start(M["starts"]["g1"]), 0, ones)                                   # compile
        t0 = time.time()
        fit = Pm.fit(Pm.start(M["starts"]["g1"]), 0, ones, M["maxiter"]["start"])
        t_fit = time.time() - t0
        t_eval = t_fit / fit["nfev"]
        # evaluations of one model on one dataset: the starts twice, the lag profile, the folds (the caps are usually reached)
        n_lag = M["lag_grid_min"][1] - M["lag_grid_min"][0]
        ev_job = 1.1 * (2 * len(M["starts"]) * M["maxiter"]["start"] + n_lag * M["maxiter"]["lag"] + n_ev * M["lag_candidates"] * M["maxiter"]["fold"])
        plan = {"t_eval_ms": round(t_eval * 1e3, 1), "t_fit_s": round(t_fit, 1), "nfev": fit["nfev"], "cuts": [], "fake_reps": FK["replicates"], "sens": list("abcde")}

        def project():
            ladders = (1 + len(FK["truths"]) * plan["fake_reps"]) * 4.2            # M0 and M1 are cheaper than M2-M4
            sens = sum({"a": 0.8, "b": 2, "c": 1, "d": 2.5, "e": 2}[k] for k in plan["sens"])
            boots = 2 * J["param_boot_reps"] * M["lag_candidates"] * M["maxiter"]["boot"] * 1.1
            return ((ladders + sens) * ev_job + boots) * t_eval / workers

        proj = project()
        plan["projected_full_h"] = round(proj / 3600, 2)
        print(f"timing: one M3 fit from a generic start {t_fit:.1f} s ({fit['nfev']} evaluations, {t_eval * 1e3:.1f} ms each); "
              f"projection {proj / 3600:.2f} h on {workers} processes (budget {S['time_budget_h']} h)", flush=True)
        for cut, change in (("sensitivity (d)", {"sens": list("abce")}), ("sensitivity (e)", {"sens": list("abc")}), ("fake replicates 3 -> 2", {"fake_reps": 2})):
            if proj <= S["time_budget_h"] * 3600:
                break
            plan.update(change), plan["cuts"].append(cut)
            proj = project()
            print(f"   cut: {cut} -> projected {proj / 3600:.2f} h", flush=True)
        plan_file.write_text(json.dumps(plan))

    # ------------------------------------------------------------------------------------------- fake check + real ladder
    nz = pd.read_csv(P.SIM_NOISE).set_index("sublink_id")
    noise = {s_: (r.sigma_white_db, r.sigma_slow_db, r.tau_noise_min) for s_, r in nz.iterrows()}
    noise["median"] = (float(nz.sigma_white_db.median()), float(nz.sigma_slow_db.median()), float(nz.tau_noise_min.median()))
    fakes, jobs = {}, []
    for truth in FK["truths"]:
        for rep in range(plan["fake_reps"]):
            name = f"fake_{truth}_{rep}"
            file = CACHE / f"data_{name}.npz"
            if not file.exists():
                save_data(make_fake(D, le, truth, rep, noise)[0], name)
            fakes[name] = (truth, rep, str(file))
            jobs += [{"key": f"{name}_{m}", "data": str(file), "model": m, "loo": True} for m in sid.MODELS]
    jobs += [{"key": f"real_{m}", "data": real_file, "model": m, "loo": True} for m in sid.MODELS]
    run_pool(sorted(jobs, key=lambda j: -sid.MODELS.index(j["model"])), workers, "ladder")

    fake_rows = []
    for name, (truth, rep, file) in fakes.items():
        tab, best, chosen, _ = ladder(name, load_data(file), sid.MODELS)
        par = result(f"{name}_{truth}")["params"]
        row = {"truth": truth, "replicate": rep, "best": best, "chosen": chosen, "chosen_is_truth": chosen == truth,
               **{f"rmse_{m}": float(tab.set_index("model").rmse[m]) for m in sid.MODELS}, "d": par["d"], "d_true": FK["lag_min"]}
        for k, v in FK["truths"][truth].items():
            row[k], row[k + "_true"] = par[k], v
        fake_rows.append(row)
    fake = pd.DataFrame(fake_rows)
    verdicts, confused = {}, []
    for truth, g in fake.groupby("truth", sort=False):
        ratios = {k: float(g[k].median() / v) for k, v in FK["truths"][truth].items()}
        ok_par = all(abs(r_ - 1) <= FK["pass"]["tol"] for r_ in ratios.values()) and abs(g.d.median() - FK["lag_min"]) <= FK["pass"]["lag_tol_min"]
        ok_cho = int(g.chosen_is_truth.sum()) >= FK["pass"]["min_chosen"]
        verdicts[truth] = {"chosen_true": int(g.chosen_is_truth.sum()), "of": len(g), "chosen": ", ".join(g.chosen), "parameters_ok": bool(ok_par),
                           "ratios": ratios, "d_median": float(g.d.median()), "PASS": bool(ok_par and ok_cho)}
        if not ok_cho:
            confused.append(f"truth {truth} is read as {', '.join(sorted(set(g.chosen[~g.chosen_is_truth])))}")
    fake_pass = all(v["PASS"] for v in verdicts.values())
    label = "" if fake_pass else " [the method did not pass on fake data]"
    save(fake.round(4), P.SID_FAKE, index=False)
    print("fake check: " + "; ".join(f"{t} {'PASS' if v['PASS'] else 'FAIL'} (chosen {v['chosen']})" for t, v in verdicts.items()), flush=True)

    tab, best, chosen, errors = ladder("real", D, sid.MODELS)
    save(tab.round(4), P.SID_LADDER, index=False)
    res = {m: result(f"real_{m}") for m in sid.MODELS}
    start_rows = [{"model": m, "d": res[m]["all"]["d"], "loss": res[m]["all"]["loss"], **{f"loss_{s}": v for s, v in res[m]["all"]["starts_dbest"].items()},
                   "spread_over_starts": max(res[m]["all"]["starts_dbest"].values()) - min(res[m]["all"]["starts_dbest"].values()),
                   "fold_d_values": ", ".join(f"{k}: {v}" for k, v in sorted(pd.Series([f["d"] for f in res[m]["folds"]]).value_counts().items())),
                   "folds_at_iteration_cap": int(sum(f["nit"] >= M["maxiter"]["fold"] for f in res[m]["folds"])), "seconds": round(res[m]["seconds"])} for m in sid.MODELS]
    starts = pd.DataFrame(start_rows)
    print(f"real ladder: best {best}, chosen by the stop rule {chosen}{label}", flush=True)

    # ------------------------------------------------------------------------------------------- parameters with intervals, sensitivity
    ref = sorted({chosen, "M1"}, key=sid.MODELS.index)
    short = le.L_km.to_numpy() < SENS["a_max_length_km"]
    has_g = (le.gauge_id != "").to_numpy()
    variants = {"a": {"data": save_data(subset(D, short), "real_a"), "keep": short, "text": f"only links < {SENS['a_max_length_km']} km"}}
    Db = dict(D)
    Db["ok5"] = D["ok5"] & (np.nan_to_num(D["r5"]) <= SENS["b_max_rain_mmh"])
    Db["y5"], Db["y1"] = np.where(D["w5"] > 0, D["y5"], np.nan), np.where(D["w1"] > 0, D["y1"], np.nan)
    variants["b"] = {"data": save_data(finish(Db), "real_b"), "keep": np.ones(len(le), bool), "no_path": True,
                     "text": f"only 5 min points with R <= {SENS['b_max_rain_mmh']:g} mm/h, no path term"}
    if le.cml_id[has_g].nunique() >= SENS["c_min_links"]:
        Dc = subset(D, has_g)
        g_nan = np.isnan(Dc["rain_gauge"])
        bad5 = np.concatenate([g_nan[o + PADM:o + PADM + n].reshape(-1, 5).any(axis=1) for o, n in zip(Dc["le_roff"], Dc["le_n"])])
        Dc["rain"] = np.nan_to_num(Dc["rain_gauge"])
        Dc["ok5"] = Dc["ok5"] & ~bad5
        Dc["y5"], Dc["y1"] = np.where(Dc["w5"] > 0, Dc["y5"], np.nan), np.where(Dc["w1"] > 0, Dc["y1"], np.nan)
        variants["c"] = {"data": save_data(finish(Dc), "real_c"), "keep": has_g, "text": "City gauge (1 min) as input, links with a gauge within 1 km"}
    variants["d"] = {"data": real_file, "keep": np.ones(len(le), bool), "agg": 1, "text": "the loss on 1 min A"}
    variants["e"] = {"data": real_file, "keep": np.ones(len(le), bool), "tau_per_band": True, "text": "separate time constants per band"}
    variants = {k: v for k, v in variants.items() if k in plan["sens"]}
    jobs = [{"key": f"real_{m}", "data": real_file, "model": m, "loo": True, "boot": J["param_boot_reps"]} for m in ref]
    for k, v in variants.items():
        for m in ref:
            if k == "e" and not sid.DYN[m]:
                continue                                   # a model without a time constant has nothing to separate
            jobs.append({"key": f"sens_{k}_{m}", "data": v["data"], "model": m, "loo": True, **{o: v[o] for o in ("agg", "no_path", "tau_per_band") if o in v}})
    run_pool(jobs, workers, "parameters and sensitivity")
    res = {m: result(f"real_{m}") for m in sid.MODELS}

    par_rows = []
    for m in ref:
        boot = pd.DataFrame(res[m]["boot"])
        for k, v in res[m]["params"].items():
            lo, hi = np.percentile(boot[k], J["interval"])
            par_rows.append({"model": m, "parameter": k, "value": v, "lo": lo, "hi": hi})
    params = pd.DataFrame(par_rows)
    save(params.round(4), P.SID_PARAMS, index=False)

    loo = {m: res[m]["loo_pred"] for m in sid.MODELS}
    sens_rows = []
    for k, v in variants.items():
        Dv = load_data(v["data"])
        main_keep = np.isin(D["bin_le"], np.nonzero(v["keep"])[0])
        for m in ref:
            key = f"sens_{k}_{m}"
            if not (CACHE / f"{key}.pkl").exists():
                continue
            r = result(key)
            wv = Dv["w5"] > 0
            err = event_errors(Dv, r["loo_pred"])
            # the main run on the same link-events and points (for (b): the same low-rain points)
            w_main = (D["w5"] > 0) & main_keep & ((np.nan_to_num(D["r5"]) <= SENS["b_max_rain_mmh"]) if k == "b" else True)
            if k == "c":
                w_main = np.zeros(len(D["w5"]), bool)
                w_main[main_keep] = wv
            main_err = event_errors(D, loo[m], w_main)
            row = {"variant": f"({k}) {v['text']}", "model": m, "link_events": len(Dv["le_n"]), "points": int(wv.sum()),
                   "rmse_variant": float(err.rmse.mean()), "rmse_main_run_same_points": float(main_err.rmse.mean()), "d": r["params"]["d"]}
            row.update({p: v_ for p, v_ in r["params"].items() if p.split("[")[0] in sid.DYN[m]})
            for b in BANDS:
                row[f"F[{b}](2)"], row[f"F[{b}](10)"], row[f"P[{b}](10)"] = r["params"].get(f"F[{b}](2)"), r["params"].get(f"F[{b}](10)"), r["params"].get(f"P[{b}](10)")
            sens_rows.append(row)
    for m in ref:                                          # the main run as the reference row
        row = {"variant": "main run", "model": m, "link_events": len(le), "points": int(w.sum()), "rmse_variant": float(errors[m].rmse.mean()),
               "rmse_main_run_same_points": float(errors[m].rmse.mean()), "d": res[m]["params"]["d"]}
        row.update({p: v_ for p, v_ in res[m]["params"].items() if p in sid.DYN[m]})
        for b in BANDS:
            row[f"F[{b}](2)"], row[f"F[{b}](10)"], row[f"P[{b}](10)"] = res[m]["params"].get(f"F[{b}](2)"), res[m]["params"].get(f"F[{b}](10)"), res[m]["params"][f"P[{b}](10)"]
        sens_rows.insert(ref.index(m), row)
    sens = pd.DataFrame(sens_rows)
    save(sens.round(4), P.SID_SENS, index=False)

    # ------------------------------------------------------------------------------------------- checks, product, figures
    curve, shower_tab, drop_tab = checks(D, le, {m: loo[m] for m in ref})
    save(curve.round(4), P.SID_AFTER, index=False)
    multi = shower_tab[shower_tab.n_showers >= 2].assign(which=lambda t_: np.where(t_.shower == 0, "first shower", "later showers"))
    sh_sum = multi.groupby("which").agg(showers=("le", "size"), link_events=("le", "nunique"), peak_R_median_mmh=("peak_R_mmh", "median"),
                                        peak_obs_median_db=("peak_obs", "median"), **{f"peak_{m}_median_db": (f"peak_{m}", "median") for m in ref}).reset_index()
    drop_tab["length"] = np.where(drop_tab.L_km < CHK["length_split_km"], f"< {CHK['length_split_km']} km", f"{CHK['length_split_km']}-{DATA['max_length_km']:g} km")
    dr_sum = drop_tab.groupby("length").agg(link_events=("le", "size"), median_L_km=("L_km", "median"), start_obs_median_db=("start_obs", "median"),
                                            drop_obs_median_db=("drop_obs", "median"), **{f"drop_{m}_median_db": (f"drop_{m}", "median") for m in ref}).reset_index()
    check_tab = pd.concat([sh_sum.assign(check="per-shower peaks, link-events with >= 2 showers").rename(columns={"which": "group"}),
                           dr_sum.assign(check="fast drop after the last wet bin").rename(columns={"length": "group"})], ignore_index=True)
    save(check_tab.round(3), P.SID_CHECKS, index=False)
    lt = link_event_table(D, loo[chosen])
    lt1 = link_event_table(D, loo["M1"])
    product = pd.concat([le.drop(columns=["band"]), lt.add_suffix(f"_{chosen}").rename(columns={f"points_{chosen}": "points", f"wet_points_{chosen}": "wet_points",
                         f"peak_obs_{chosen}": "peak_obs_db", f"integral_obs_{chosen}": "integral_obs_db_min"})], axis=1)
    if chosen != "M1":
        product = pd.concat([product, lt1[["peak_pred", "integral_pred", "bias", "rmse"]].add_suffix("_M1")], axis=1)
    product.insert(0, "chosen_model", chosen)
    num = product.select_dtypes("number").columns
    product[num] = product[num].round(4)
    save(product, P.SID_LINK_EVENTS, index=False)

    fig = figure_ladder(tab, chosen, best)
    save(fig, P.SID_FIG_LADDER, dpi=110), plt.close(fig)
    cand = np.nonzero((lt.peak_obs - le.offset_c_db >= CHK["min_peak_db"]).to_numpy())[0]
    picks = np.sort(np.random.default_rng(S["gallery"]["seed"]).choice(cand, S["gallery"]["link_events"], replace=False))
    fig = figure_gallery(D, le, picks, loo, chosen)
    save(fig, P.SID_FIG_GALLERY, dpi=100), plt.close(fig)
    fig = figure_curves(res[chosen]["params"], res[chosen]["boot"], chosen)
    save(fig, P.SID_FIG_CURVES, dpi=110), plt.close(fig)
    fig = figure_after(curve, chosen)
    save(fig, P.SID_FIG_AFTER, dpi=110), plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    import numba
    import scipy
    show = ["model", "rmse", "rmse_lo", "rmse_hi", "minus_best", "diff_lo", "diff_hi", "interval_holds_0", "skill_vs_M0", "skill_lo", "skill_hi"]
    a0, a1 = J["after_windows_min"]
    more = ["model", "rmse_rain", f"rmse_{a0[0]}_{a0[1]}", f"rmse_{a1[0]}_{a1[1]}", "peak_error", "peak_error_lo", "peak_error_hi", "integral_obs", "integral_pred", "integral_pred_over_obs"]
    fake_show = fake[["truth", "replicate", "best", "chosen", "chosen_is_truth"] + [f"rmse_{m}" for m in sid.MODELS] + ["d"]
                     + [c for k in ("tau_wet", "tau_dry", "phi", "tau_f_wet", "tau_f_dry", "tau_s") if k in fake for c in (k,)]].round(3)
    key_par = params[~params.parameter.str.startswith(("F[", "P["))].round(3)
    curve_par = params[params.parameter.str.startswith(("F[", "P[")) & (params.model == chosen)].round(3)
    ctx = {k: v for k, v in locals().items()}
    plain, unclear, equation = texts(ctx)
    report = f"""# wet_antenna step S1: system identification of WAA on short links, the model ladder M0-M4 (TRAIN only){label}

Written by `research/wet_antenna/system_id_ladder.py` (library `system_id.py`). Settings and every rule:
`wet_antenna.system_id` in `research/config.yaml`, written before the run. Plan: `docs/ROADMAP.md`, "Reboot (2026-10-13)".
{'' if fake_pass else chr(10) + '**THE METHOD DID NOT PASS ON FAKE DATA (section 3): ' + '; '.join(confused or ['the parameters are not recovered']) + '. Every real result below carries that label.**' + chr(10)}
## 1. What was run

- Every model predicts A_pred(t) = W(t) + L * P_band(R(t - d)) + c at 1 min and is compared with the 5 min means of A.
  M0 no WAA; M1 static W = a_link F_band(R); M2 one time constant; M3 separate wetting / drying time constants; M4 two
  stores. P and F are free piecewise-linear curves on the knots {', '.join(f'{k:g}' for k in KNOTS)} mm/h (0 at R = 0, never
  falling). **No power-law rain term is subtracted anywhere**; P is told apart from W only through L. ITU-R P.838 is
  used only to build the fake data and as a dotted size reference in one figure.
- Loss: Huber (delta {M['huber_delta_db']} dB) on the 5 min residuals, minimised with L-BFGS-B and the exact gradient; {len(M['starts'])} generic starts
  per model; the lag d from the grid {M['lag_grid_min'][0]}..{M['lag_grid_min'][1]} min by a lag profile (config `lag_candidates`).
- Judge: leave-one-event-out over the {n_ev} TRAIN events ({J['boot_reps']} event-bootstrap draws for the intervals).
- TRAIN events only; VAL is never read. Nothing earlier is changed.
- Timing: one M3 fit from a generic start {plan['t_fit_s']} s ({plan['nfev']} evaluations of {plan['t_eval_ms']} ms); projection {plan['projected_full_h']} h on {workers}
  processes against a budget of {S['time_budget_h']} h. Cut: {', '.join(plan['cuts']) or 'nothing'}. Fitting time summed over the processes:
  {sum(result(j_.stem)['seconds'] for j_ in CACHE.glob('*.pkl')) / 3600:.2f} h. Run time of this run: {(time.time() - clock) / 60:.1f} min (fits are cached).
- Software: numba {numba.__version__} (the model and its gradient), scipy {scipy.__version__} (L-BFGS-B), numpy {np.__version__}. No pysindy, no PySR in this step.

## 2. Data

{md_table(counts)}

(`{rel(P.SID_COUNTS)}`.) A link-event = one TRAIN event on one sub-link. How the pairs were reduced (rules decided
without any fit):

{md_table(drops_tab)}

Path lengths {le.L_km.min():.2f} to {le.L_km.max():.2f} km (median {le.L_km.median():.2f}). Radar R_mean on wet points: median {np.nanmedian(D['r5'][w & D['wet5']]):.2f}, 90th percentile
{np.nanpercentile(D['r5'][w & D['wet5']], 90):.1f}, largest {np.nanmax(D['r5'][w]):.0f} mm/h; points above 10 mm/h: {int((D['r5'][w] > 10).sum())}.

## 3. The fake check (first)

Truths M3 (tau_wet {FK['truths']['M3']['tau_wet']}, tau_dry {FK['truths']['M3']['tau_dry']} min) and M4 (phi {FK['truths']['M4']['phi']}, tau_f_wet {FK['truths']['M4']['tau_f_wet']}, tau_f_dry {FK['truths']['M4']['tau_f_dry']}, tau_s {FK['truths']['M4']['tau_s']} min), F = {FK['f_top_db']} (1 - exp(-R / {FK['f_scale_mmh']:g})) dB,
P = ITU-R P.838 size, a_link log-normal (sd {FK['a_sd']}), lag {FK['lag_min']} min, driven by the real radar rain, windows and lengths; noise
variance x {FK['noise_variance_factor']:g}, {FK['rsl_step_db']} dB quantization, unknown offsets. {plan['fake_reps']} replicates per truth, the same ladder and judge.

{md_table(fake_show)}

(`{rel(P.SID_FAKE)}`; rmse = held-out RMSE; the parameters are those of the fit of the TRUE rung on all events.)

{chr(10).join(f"- **Truth {t}: {'PASS' if v_['PASS'] else 'FAIL'}.** Chosen: {v_['chosen']} (true rung in {v_['chosen_true']} of {v_['of']}; rule >= {FK['pass']['min_chosen']}). Fitted / true: " + ', '.join(f'{k} x{r_:.2f}' for k, r_ in v_['ratios'].items()) + f"; d {v_['d_median']:.0f} (true {FK['lag_min']}). Parameters within +-{FK['pass']['tol']:.0%}: {'yes' if v_['parameters_ok'] else 'no'}." for t, v_ in verdicts.items())}

**Fake check: {'PASS' if fake_pass else 'FAIL'}.**{'' if fake_pass else ' Rungs confused: ' + ('; '.join(confused) or 'none (parameters only)') + '.'}

## 4. The ladder on the real data{label}

{md_table(tab[show].round(4))}

(`{rel(P.SID_LADDER)}`; rmse = mean over the held-out events of the RMSE of A (dB); minus_best with the 90 % interval of
the paired difference; skill = 1 - RMSE / RMSE(M0).) More metrics (integrals in dB min, mean per link-event):

{md_table(tab[more].round(3))}

Fits on all events, and how the starts and folds behaved:

{md_table(starts.round(2))}

**Best: {best}. Chosen by the stop rule (the simplest model whose paired difference to the best has an interval that
holds 0): {chosen}.**{label}

![ladder](../figures/{P.SID_FIG_LADDER.name})

![gallery](../figures/{P.SID_FIG_GALLERY.name})

## 5. The chosen model{label}

{equation}

{md_table(key_par)}

F and P of {chosen} at the knots (value, 90 % interval from {J['param_boot_reps']} event-bootstrap refits):

{md_table(curve_par)}

(`{rel(P.SID_PARAMS)}`.)

![curves](../figures/{P.SID_FIG_CURVES.name})

## 6. Sensitivity (information, not part of the rule)

{md_table(sens.round(3))}

(`{rel(P.SID_SENS)}`. rmse_variant = held-out RMSE of the variant; rmse_main_run_same_points = the main run judged on
the same link-events and points.)

## 7. Checks against earlier facts (held-out predictions)

{md_table(check_tab.round(3))}

(`{rel(P.SID_CHECKS)}`; the curve after the last shower: `{rel(P.SID_AFTER)}`.)

![after the last shower](../figures/{P.SID_FIG_AFTER.name})

## 8. In plain words

{plain}

## 9. Unclear / possibly wrong

{unclear}
"""
    save(report, P.SID_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {(time.time() - clock) / 60:.1f} min")


def texts(q):
    """Sections 5 (the equation), 8 and 9 of the report. The numbers are taken from the run; the sentences that READ
    them were written after the run of 2026-10-13 (chosen M4) and say so where they depend on that outcome."""
    chosen, tab, params, sens, curve, le, D, res = q["chosen"], q["tab"].set_index("model"), q["params"], q["sens"], q["curve"], q["le"], q["D"], q["res"]
    par = params[params.model == chosen].set_index("parameter")
    pv = lambda k, n=1: f"{par.value[k]:.{n}f} ({par.lo[k]:.{n}f}..{par.hi[k]:.{n}f})"        # noqa: E731
    sv = lambda variant, model, col: float(sens[sens.variant.str.startswith(variant) & (sens.model == model)][col].iloc[0])   # noqa: E731
    has = lambda variant: bool(sens.variant.str.startswith(variant).any())                     # noqa: E731
    dyn = sid.DYN[chosen]
    W = {"M0": "W = 0", "M1": "W = a_link F_band(R(t - d))", "M2": "dW/dt = (a_link F_band(R(t - d)) - W) / tau",
         "M3": "dW/dt = (a_link F_band(R(t - d)) - W) / tau,   tau = tau_wet while a_link F > W, tau_dry otherwise",
         "M4": "W = W_f + W_s\n    dW_f/dt = (phi a_link F_band(R(t - d)) - W_f) / tau_f,   tau_f = tau_f_wet while phi a_link F > W_f, tau_f_dry otherwise\n"
               "    dW_s/dt = ((1 - phi) a_link F_band(R(t - d)) - W_s) / tau_s"}[chosen]
    at_bound = [k for k in dyn if k.startswith("tau") and par.hi[k] >= 0.99 * M["bounds"]["tau_min"][1]]
    dyn_lines = "\n".join(f"- {k} = {pv(k, 2 if k == 'phi' else 1)}"
                          + (f"  **(the interval reaches the upper bound of {M['bounds']['tau_min'][1]} min: not identified, only 'long')**" if k in at_bound else "") for k in dyn)
    equation = f"""```
A_pred(t) = W(t) + L * P_band(R(t - d)) + c
    {W}
```

Fit on all {q['n_ev']} TRAIN events; 90 % intervals from {J['param_boot_reps']} event-bootstrap refits (time constants and d in minutes;
positive d = A follows the radar):

- radar lag d = {pv('d', 0)} min;
{dyn_lines}
- a_link: sd of log a_link {pv('a_log_sd', 2) if 'a_log_sd' in par.index else 'n/a'}; 10th to 90th percentile {par.value.get('a_p10', np.nan):.2f} to {par.value.get('a_p90', np.nan):.2f}."""
    if chosen == "M0":
        return "The stop rule chose M0: no WAA term is needed to predict A on these links. See the tables.", "See the tables.", equation

    # ---- numbers for the plain-words part
    Lmed = float(le.L_km.median())
    b0, b1 = BANDS
    F = lambda b, r: par.value[f"F[{b}]({r:g})"]            # noqa: E731
    Pk = lambda b, r: par.value[f"P[{b}]({r:g})"]           # noqa: E731
    itu = {b: itu_p838(float(le[le.band_name == b].frequency_ghz.median()), "v") for b in BANDS}
    itu10 = {b: itu[b][0] * 10 ** itu[b][1] for b in BANDS}
    p_all = ", ".join(f"{m} {res[m]['params'][f'P[{b0}](10)']:.1f} / {res[m]['params'][f'P[{b1}](10)']:.1f}" for m in sid.MODELS)
    cur = curve.set_index("minutes_after_last_wet_bin")
    at = lambda col, t: float(cur[col][t])                  # noqa: E731
    chk = q["check_tab"].set_index("group")
    w = q["w"]
    n_hi = int((D["r5"][w] > 10).sum())
    fake_rmse = float(q["fake"][f"rmse_{chosen}"].mean())
    dry_name = {"M4": "tau_f_dry", "M3": "tau_dry"}.get(chosen, "tau")
    wet_name = {"M4": "tau_f_wet", "M3": "tau_wet"}.get(chosen, "tau")
    dry_sens = ", ".join(f"{sv(f'({k})', chosen, dry_name):.1f}" for k in "abcd" if has(f"({k})")) if dyn else ""
    gauge = has("(c)")
    top = chosen == sid.MODELS[-1]

    plain = f"""**Static or dynamic?** {'Static by the stop rule.' if chosen == 'M1' else 'Dynamic.'} The static model M1 is {tab.minus_best['M1']:.3f} dB worse than the best held out
({tab.diff_lo['M1']:.3f}..{tab.diff_hi['M1']:.3f}; the interval {'holds' if tab.interval_holds_0['M1'] else 'does not hold'} 0), and after the last shower it falls to 0 at once, while the observed A / peak is still
{at('observed', 0):.2f} at 0-5 min and {at('observed', 15):.2f} at 15-20 min (M1: {at('M1', 0):.2f} and {at('M1', 15):.2f}). **But the gain is small in dB**: the held-out RMSE goes from
{tab.rmse['M0']:.3f} (M0) over {tab.rmse['M1']:.3f} (M1) to {tab.rmse[chosen]:.3f} dB ({chosen}), a skill of {tab.skill_vs_M0[chosen]:.1%} ({tab.skill_lo[chosen]:.1%}..{tab.skill_hi[chosen]:.1%}) over "no WAA at all". Most of the error is
not removed by any rung: {tab.rmse_60_180[chosen]:.2f} dB remain 1-3 h after the rain (all rungs alike: noise and baseline), and {tab.rmse_rain[chosen]:.2f} dB during rain.
The stop rule chose {'the TOP rung, so the ladder did not stop: nothing here says that M4 is enough (that is S2)' if top else 'a rung below the top'}.
"""
    if dyn:
        plain += f"""
**How fast does it dry?** The one number that is stable: {dry_name} = {pv(dry_name)} min. In the sensitivity runs (a)-(d) it is
{dry_sens} min, with the radar or with the gauge as input. It agrees with the fast first phase of the tails (T17 / T18: fast
for 15-20 min).

**How fast does it wet?** NOT determined. With the radar {wet_name} = {pv(wet_name)} min, slower than drying, the opposite of what
[1] describes"""
        if gauge:
            plain += f"""; with the City gauge (1 min, a point on the ground) as input on the {int(sv('(c)', chosen, 'link_events'))} link-events near a gauge it is
{sv('(c)', chosen, wet_name):.1f} min, and the held-out RMSE falls from {sv('(c)', chosen, 'rmse_main_run_same_points'):.3f} (radar, the same points) to {sv('(c)', chosen, 'rmse_variant'):.3f} dB. So the slow wetting is most
likely the radar (5 min steps, 2 km cells, measured aloft) and not the antenna: a smoothed input makes the model rise slowly"""
        plain += "."
    if chosen == "M4":
        plain += f"""

**The slow store is not identified.** tau_s = {pv('tau_s', 0)} min: the interval ends at the bound. A store that slow is almost an
integrator of the rain: it predicts A / peak = {at(chosen, 60):.3f} at 60 min and {at(chosen, 120):.3f} at 120 min after the last shower, where the
observed values are {at('observed', 60):.3f} and {at('observed', 120):.3f}, and an integral of A over 0-180 min after the rain {tab.integral_pred_over_obs[chosen]:.1f} times the observed one
({tab.integral_pred[chosen]:.1f} against {tab.integral_obs[chosen]:.1f} dB min per link-event; M3: {tab.integral_pred['M3']:.1f}). It helps DURING long events (a level that builds up over several
showers) and is wrong AFTER them. It is not a drying time of {par.value['tau_s'] / 60:.0f} h; it may as well be a slow drift of the level
that the offset rule does not see. As a "false rain after the event" proxy M4 is WORSE than M2 / M3."""
    plain += f"""

**How big is the path part on short links?** For {chosen}, P at 10 mm/h is {pv(f'P[{b0}](10)')} dB/km at {b0} and {pv(f'P[{b1}](10)')} at {b1}
(ITU-R P.838, size reference only: {itu10[b0]:.1f} and {itu10[b1]:.1f}), so on the median link ({Lmed:.2f} km) about {Pk(b0, 10) * Lmed:.1f} and {Pk(b1, 10) * Lmed:.1f} dB, against a WAA
level F(10) of {F(b0, 10):.1f} and {F(b1, 10):.1f} dB (for a_link = 1). At 2 mm/h (the median wet point has {np.nanmedian(D['r5'][w & D['wet5']]):.2f} mm/h) the path part is {Pk(b0, 2) * Lmed:.2f} and {Pk(b1, 2) * Lmed:.2f} dB:
most of A on these links, most of the time, is the antenna term. **The split is not firm**: P(10) changes with the rung
({p_all} dB/km)"""
    if has("(a)"):
        plain += f" and with the length range (below {SENS['a_max_length_km']} km only: {sv('(a)', chosen, f'P[{b0}](10)'):.1f} / {sv('(a)', chosen, f'P[{b1}](10)'):.1f})"
    gain_c = (f"the gauge run gains {sv('(c)', chosen, 'rmse_main_run_same_points') - sv('(c)', chosen, 'rmse_variant'):.3f} dB on its link-events, the whole ladder M0 -> {chosen} gains "
              f"{tab.rmse['M0'] - tab.rmse[chosen]:.3f} dB") if gauge else "see (c)"
    plain += f""". With L between
{le.L_km.min():.2f} and {le.L_km.max():.2f} km the lengths do not separate P from F well; what a rung puts into W it takes out of P.

**What changes for S2** (proposals; the learners decide):
1. The input limits the result more than the model form: {gain_c}.
   S2 should be run (or at least checked) with the gauge input on the links near a gauge, and the wetting speed read there.
2. The peaks are too low for every rung (mean peak error {tab.peak_error[chosen]:.2f} dB for {chosen}, {tab.peak_error['M1']:.2f} for M1): a flexible G(R, W) cannot
   repair what the radar does not show. S2 needs a metric on the peaks, not only the RMSE.
3. The slow store needs another form or another explanation before it is "named": S2 should leave the slow part
   free (or test it against a level drift), and judge it on the 60-180 min window and the integral after the rain.
4. F above 10 mm/h rests on {n_hi} points of 5 min ({n_hi / w.sum():.2%}); its rise there (F(40) = {pv(f'F[{b0}](40)')} dB at {b0}) is not a result.
5. P and F are only separated by L: the long links (planned for S3) are needed to pin P; until then "the WAA size" of
   this step is F plus or minus an unknown share of P."""

    sh = lambda g, c: float(chk[c][g])                      # noqa: E731
    groups = sorted(g for g in chk.index if "km" in g)
    low = 1 - sh("later showers", f"peak_{chosen}_median_db") / sh("later showers", "peak_obs_median_db")
    unclear = f"""1. **The fake check passed, but it is easier than the real data.** The fake A is driven by the SAME radar series that the
   fit gets, so the input has no error; the real input has (sensitivity (c)). The fake held-out RMSE of {chosen} is {fake_rmse:.3f} dB,
   the real one {tab.rmse[chosen]:.3f} dB: the real data hold structure that the fake data do not have. PASS means "the code and the
   judge find the rung and the numbers when the model class is right and the rain is known", not more.
2. **The fake truths had fast wetting (3 min)**; the real fit gives slow wetting with the radar. A fake check with a
   smoothed or displaced rain input was not run, so "the slow wetting is the radar" is a reading of (c), not a test.
   (c) also differs from the main run in its links ({int(sv('(c)', chosen, 'link_events')) if gauge else 0} of {len(le)} link-events).
3. **Per-shower peaks** (link-events with >= 2 showers). Observed median peak {sh('first shower', 'peak_obs_median_db'):.2f} dB (first shower) and {sh('later showers', 'peak_obs_median_db'):.2f} dB (later
   showers) for a median radar peak of {sh('first shower', 'peak_R_median_mmh'):.2f} and {sh('later showers', 'peak_R_median_mmh'):.2f} mm/h. {chosen}: {sh('first shower', f'peak_{chosen}_median_db'):.2f} and {sh('later showers', f'peak_{chosen}_median_db'):.2f} (the right direction, {low:.0%} too low);
   M1: {sh('first shower', 'peak_M1_median_db'):.2f} and {sh('later showers', 'peak_M1_median_db'):.2f} (the wrong direction). Only a rung with memory makes the later showers higher.
4. **Fast drop after the last wet bin** (A at 0-5 min minus A at {CHK['fast_window_min'][0]}-{CHK['fast_window_min'][1]} min): observed {', '.join(f"{sh(g, 'drop_obs_median_db'):.2f} dB ({g})" for g in groups)}; {chosen}
   {', '.join(f"{sh(g, f'drop_{chosen}_median_db'):.2f}" for g in groups)}. It does not grow with L in the data (T18) nor in the model, but the model gives less than half of
   it: the predicted A is already too low when the rain ends (A / peak {at(chosen, 0):.2f} against {at('observed', 0):.2f}).
5. **The radar lag** d = {pv('d', 0)} min for {chosen}, but {', '.join(f"{m} {res[m]['params']['d']}" for m in sid.MODELS if m != chosen)} for the other rungs: d and the wetting time constant
   trade against each other (a slow wetting does the work of a lag). d is not a measured delay.
6. **Optimiser.** Spread of the loss over the {len(M['starts'])} starts: {', '.join(f"{r.model} {r.spread_over_starts:.1f}" for r in q['starts'].itertuples())} (the best loss is about {q['starts'].loss.min():.0f}).
   The leave-one-event-out and bootstrap refits are warm starts with an iteration cap (config `maxiter`; column
   folds_at_iteration_cap in section 4): the held-out errors are those of nearly, not fully, converged fits, for every rung alike.
7. **The offset c.** {int(le.offset_flag.sum())} of {len(le)} link-events had fewer than {DATA['offset_min_minutes']} dry minutes before the rain and got c = 0; they are in the
   fit and in the judge. c is measured 30-120 min before the rain and held for a window of {le.n_bins.median() / 12:.0f} h (median): a baseline
   that moves inside the window goes into the residual, or into the slow store.
8. **a_link of a link that is not in the training events is 1**; with leave-one-EVENT-out nearly every link is seen, so
   the judge does not test new links (that is S3).
9. **Huber loss** (delta {M['huber_delta_db']} dB) gives the large residuals at the peaks a linear, not a squared, weight: part of the low
   peaks is the loss, not only the radar. The judge's RMSE is squared, the fit is not.
10. The intervals are over events ({q['n_ev']}); the events differ a lot in size, so the mean over events is not the mean over points."""
    return plain, unclear, equation


if __name__ == "__main__":
    main()
