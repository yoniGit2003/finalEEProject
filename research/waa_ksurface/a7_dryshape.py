"""Diagnostic of the dry-path conflict (decisions.md, "Diagnostic of the dry-path conflict", logged before the run).

The A7 review found that the CV selection picks K1 while the dry-path medians favour K0. This script asks why, with
no refitting: K1 is the A7 fit (a7_real_fit/fit.pkl), K0 the refit of the review (a7_review/sel_mid.pkl).

    python -m research.waa_ksurface.a7_dryshape

    (1) per dry-path window: the decay rate of ln(A - c) in the early half (5-15 min) and the late half (15-30 min), and
        their ratio late / early (< 1: drying slows within a window, the K1 shape; about 1: exponential); the same on the
        K1 and K0 predictions in the same windows; medians with event-bootstrap intervals
    (2) every quantity for all windows and per start-level class (first dry bin >= 1 dB / < 1 dB)
    (3) level error at the first dry bin and decay-rate error, per class and model
    (4) 1 / k_dry of K0 at each k_wet estimate

Real TRAIN data of the 14 learning events. VAL and TEST are not opened.
"""
import pickle

import numpy as np
import pandas as pd

from research.common.report import md_table
from research.common.save import save
from research.waa_ksurface import a7_review as rv, data, fake6, ksurf
from research.waa_ksurface import paths as P

C = data.branch_config()["a7d"]
log = rv.log


# =================================================================================================== the windows and rates
def windows(D, rain1, use):
    """(link-event, bin indices of the lags) of every window kept: dry path, counted, finite, measured A - c > min_level."""
    idx, target = rv.dry_values(D, rain1, C["lags_min"], rv.R7["rain_before_min"], rv.S1["dry_below_mmh"], use)
    keep = np.all(idx >= 0, axis=1)
    keep[keep] = np.all(target[idx[keep]] > C["min_level_db"], axis=1)
    return np.nonzero(keep)[0], idx[keep], target


def rate(t, y):
    """-(OLS slope of ln y against t) per row (1/min); NaN where any y <= 0."""
    with np.errstate(divide="ignore", invalid="ignore"):
        ly = np.where(y > 0, np.log(np.where(y > 0, y, 1.0)), np.nan)
    tc = t - t.mean()
    return -(ly - ly.mean(axis=1, keepdims=True)) @ tc / (tc @ tc)


def ratios(Y):
    """Y = values at the lags (rows: windows) -> early rate, late rate, ratio late / early (NaN if early <= 0)."""
    lags = np.asarray(C["lags_min"], float)
    ie = [C["lags_min"].index(m) for m in C["early_min"]]
    il = [C["lags_min"].index(m) for m in C["late_min"]]
    re, rl = rate(lags[ie], Y[:, ie]), rate(lags[il], Y[:, il])
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(re > 0, rl / re, np.nan)
    return re, rl, r


# =================================================================================================== the bootstrap
def boot(ev, fns, reps, seed):
    """90 % event-bootstrap intervals of each statistic fn(rows) (rows = window indices; a drawn event's windows once per copy)."""
    u = np.unique(ev)
    by = {e: np.nonzero(ev == e)[0] for e in u}
    rng = np.random.default_rng(seed)
    out = np.full((reps, len(fns)), np.nan)
    for b in range(reps):
        rows = np.concatenate([by[e] for e in rng.choice(u, len(u), replace=True)])
        out[b] = [f(rows) for f in fns]
    return np.nanpercentile(out, 5, axis=0), np.nanpercentile(out, 95, axis=0)


def summary(T, est):
    """Per class: n, median ratios (measured, K1, K0), level and rate errors, each with its interval."""
    t = T[T["k_wet estimate"] == est]
    t = t[np.isfinite(t.ratio_meas) & np.isfinite(t.ratio_K1) & np.isfinite(t.ratio_K0)]
    rows = []
    for cls, g in [("all", t), (f">= {C['class_db']:g} dB", t[t.start_meas >= C["class_db"]]), (f"< {C['class_db']:g} dB", t[t.start_meas < C["class_db"]])]:
        if not len(g):
            rows.append({"k_wet estimate": est, "class": cls, "windows": 0, "events": 0})
            continue
        a = {k: g[k].to_numpy() for k in ("ratio_meas", "ratio_K1", "ratio_K0", "start_meas", "start_K1", "start_K0")}
        stats = {"ratio measured": lambda r: np.median(a["ratio_meas"][r]),
                 "ratio K1": lambda r: np.median(a["ratio_K1"][r]),
                 "ratio K0": lambda r: np.median(a["ratio_K0"][r]),
                 "level error K1 (dB)": lambda r: np.median(a["start_K1"][r]) - np.median(a["start_meas"][r]),
                 "level error K0 (dB)": lambda r: np.median(a["start_K0"][r]) - np.median(a["start_meas"][r]),
                 "rate error K1": lambda r: np.median(a["ratio_K1"][r]) - np.median(a["ratio_meas"][r]),
                 "rate error K0": lambda r: np.median(a["ratio_K0"][r]) - np.median(a["ratio_meas"][r])}
        allr = np.arange(len(g))
        lo, hi = boot(g.event.to_numpy(), list(stats.values()), C["boot_reps"], C["boot_seed"])
        row = {"k_wet estimate": est, "class": cls, "windows": len(g), "events": int(g.event.nunique()),
               "median start A - c measured (dB)": float(np.median(a["start_meas"]))}
        for (k, f), l, h in zip(stats.items(), lo, hi):
            row[k] = float(f(allr))
            row[f"{k} 5 %"], row[f"{k} 95 %"] = float(l), float(h)
        rows.append(row)
    return rows


# =================================================================================================== main
def main():
    S, curves, band_idx, use = rv.setting("mid")
    D = S["D"]
    sel = pickle.loads((rv.CACHE / "sel_mid.pkl").read_bytes())
    fit = pickle.loads((rv.A7_CACHE / "fit.pkl").read_bytes())
    prob = ksurf.prepare(D, band_idx, curves, D["rain"], 0)
    le_w, idx, target = windows(D, fake6.window(D, np.asarray(D["rain"], float)), use)
    ev = S["le"].event_id.to_numpy()[le_w]
    log(f"   (dryshape) windows kept: {len(le_w)} link-events, {len(np.unique(ev))} events")
    Ym = target[idx]
    re_m, rl_m, r_m = ratios(Ym)
    tabs, rows, kd = [], [], []
    for e in rv.ESTIMATES:
        k0_lam, k0_x = sel["est"][e]["fits"]["K0"]
        m1 = ksurf.Model(prob, "K0", fit["rung"], fit["lam"])
        m0 = ksurf.Model(prob, "K0", "K0", k0_lam)
        out = {"K1": m1.simulate(fit["x"][e])[1][idx], "K0": m0.simulate(k0_x)[1][idx]}
        T = pd.DataFrame({"k_wet estimate": e, "link_event": le_w, "event": ev, "start_meas": Ym[:, 0],
                          "rate_early_meas": re_m, "rate_late_meas": rl_m, "ratio_meas": r_m})
        for k, Y in out.items():
            a, b, r = ratios(Y)
            T[f"start_{k}"], T[f"rate_early_{k}"], T[f"rate_late_{k}"], T[f"ratio_{k}"] = Y[:, 0], a, b, r
        tabs.append(T)
        ok = np.isfinite(T.ratio_meas) & np.isfinite(T.ratio_K1) & np.isfinite(T.ratio_K0)
        log(f"   (dryshape) [{e}] windows with all three ratios: {int(ok.sum())} of {len(T)}")
        rows += summary(T, e)
        kd.append({"k_wet estimate": e, "K0: 1/k_dry (min)": round(float(1 / m0.k(k0_x, "drying", 0, 0)), 2),
                   **{f"K1 (A7): 1/k_dry at W = {W:g} dB (min)": round(float(1 / m1.k(fit["x"][e], "drying", 0, W)), 1) for W in (0, 0.5, 1, 2)}})
    T, Sm, K = pd.concat(tabs, ignore_index=True), pd.DataFrame(rows), pd.DataFrame(kd)
    save(T, P.A7D_WIN, index=False)
    save(Sm, P.A7D_SUM, index=False)
    report(T, Sm, K)


def ci(row, k, fmt="{:.3f}"):
    return f"{fmt.format(row[k])} ({fmt.format(row[k + ' 5 %'])} to {fmt.format(row[k + ' 95 %'])})"


def report(T, Sm, K):
    L = ["# Diagnostic of the dry-path conflict\n",
         "Written by `research/waa_ksurface/a7_dryshape.py` (decisions.md, \"Diagnostic of the dry-path conflict\", logged before the run). No "
         "refitting: K1 = the A7 fit, K0 = the refit of the A7 review. Real TRAIN data of the 14 learning events; VAL and TEST not opened.\n",
         f"Windows: the dry-path windows of review check (5) with the bins ending 5-30 min after the path went dry all present and measured "
         f"A - c > {C['min_level_db']} dB in each. Rate = -(slope of ln(A - c)), early on {C['early_min']} min, late on {C['late_min']} min; "
         f"ratio = late / early. Intervals: 90 % event bootstrap ({C['boot_reps']} draws, seed {C['boot_seed']}).\n"]
    for e in rv.ESTIMATES:
        t = T[T["k_wet estimate"] == e]
        ok = np.isfinite(t.ratio_meas) & np.isfinite(t.ratio_K1) & np.isfinite(t.ratio_K0)
        L.append(f"- k_wet {e}: {len(t)} windows kept, {int(ok.sum())} with all three ratios (dropped: measured {int((~np.isfinite(t.ratio_meas)).sum())}, "
                 f"K1 {int((~np.isfinite(t.ratio_K1)).sum())}, K0 {int((~np.isfinite(t.ratio_K0)).sum())}).")
    L.append("")
    a = []
    for _, r in Sm.iterrows():
        if not r.windows:
            a.append({"k_wet estimate": r["k_wet estimate"], "class": r["class"], "windows / events": "0 / 0"})
            continue
        a.append({"k_wet estimate": r["k_wet estimate"], "class": r["class"], "windows / events": f"{int(r.windows)} / {int(r.events)}",
                  "ratio measured": ci(r, "ratio measured"), "ratio K1": ci(r, "ratio K1"), "ratio K0": ci(r, "ratio K0")})
    L += ["## (1), (2) Ratio rate(late) / rate(early), median (90 %)\n", md_table(pd.DataFrame(a)), ""]
    b = [{"k_wet estimate": r["k_wet estimate"], "class": r["class"], "median start measured (dB)": round(r["median start A - c measured (dB)"], 3),
          "level error K1 (dB)": ci(r, "level error K1 (dB)"), "level error K0 (dB)": ci(r, "level error K0 (dB)"),
          "rate error K1": ci(r, "rate error K1"), "rate error K0": ci(r, "rate error K0")} for _, r in Sm.iterrows() if r.windows]
    L += ["## (3) Errors: level at the first dry bin (model - measured, medians) and decay-rate ratio (model - measured)\n", md_table(pd.DataFrame(b)), "",
          "## (4) The drying rates used\n", md_table(K), ""]
    save("\n".join(L) + "\n", P.A7D_REPORT)
    log(f"report: {P.A7D_REPORT}")


if __name__ == "__main__":
    main()
