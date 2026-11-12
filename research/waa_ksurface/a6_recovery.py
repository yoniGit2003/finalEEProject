"""Step A6 of the k-surface branch: does the fit find a KNOWN k on fake data? The gates for the real fit (A7).

    python -m research.waa_ksurface.a6_recovery                 every truth in order (a finished truth is read from the cache)
    python -m research.waa_ksurface.a6_recovery c d             only these truths
    python -m research.waa_ksurface.a6_recovery --force         run again, ignoring the cache
    python -m research.waa_ksurface.a6_recovery --report        write the report from the cache

Fake attenuation on the real links and the real rain (library `fake6.py`), fitted with the procedure of A5 (library
`ladder.py`): stage 1 learns k_wet from the gauge-input set, stage 2 learns k_dry from all links. Truths, PASS rules and
what a failure means: decisions.md, step A6 (written before any run). A FAIL on (a), (c) or (d) is a STOP: the script
writes the report and exits with code 1. No real attenuation is used as a target.
"""
import pickle
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import data, fake6, fp, ksurf, ladder, ltest
from research.waa_ksurface import paths as P

A5, A6 = ksurf.config(), data.branch_config()["a6"]
LT = ltest.config()
EDGES = np.asarray(LT["class_edges_mmh"], float)
ORDER = ["a", "a_lag", "b", "c", "d", "e", "f", "g"]
GATES = ["a", "c", "d"]                                         # a FAIL here is a STOP
J = A6["judge"]


def log(msg):
    print(msg, flush=True)


# =================================================================================================== the setting
def setting():
    d = data.open_release()
    D, le, n_cells = data.gauge_ratio_train(d)
    gz = d.table("gauges").set_index("gauge_id")
    step = le.gauge_id.map(gz.quantization_mm).to_numpy(float)              # NaN for a link without a gauge
    dist = le.gauge_distance_m.to_numpy(float) / 1000.0
    curves, iF, iP = fake6.fake_curves(fp.measured())
    dom = pd.read_csv(P.A4_DOMAIN)
    dom = dom[(dom.band == "pooled") & dom.in_domain][["regime", "rcol", "wbin", "R_class", "W_lo"]]
    return {"D": D, "le": le, "n_cells": n_cells, "step": step, "dist": dist, "curves": curves, "iF": iF, "iP": iP, "dom": dom,
            "gauge": (le.gauge_id != "").to_numpy().astype(float), "short": (le.L_km < LT["short_km"]).to_numpy()}


# =================================================================================================== judging k
def state_density(S, Dk, ex):
    """Minutes the fake world's TRUE state spends in every cell of the frozen A4 domain (short links; R as the fit sees it)."""
    D = S["D"]
    rep = np.repeat(np.arange(len(S["le"])), D["le_n"])
    sh = S["short"][rep]
    r = fake6.window(D, Dk["rain"])[sh]
    W, v = ex["W_true"][sh], ex["v_true"][sh]
    t = pd.DataFrame({"regime": np.where(v > W, "wetting", "drying"),
                      "rcol": np.where(r < EDGES[0], 0, np.searchsorted(EDGES, r, side="right")), "wbin": np.floor((W + 0.3) / 0.3).astype(int)})
    cnt = t.groupby(["regime", "rcol", "wbin"]).size().rename("minutes").reset_index()
    return S["dom"].merge(cnt, on=["regime", "rcol", "wbin"], how="left").fillna({"minutes": 0})


def judge(S, spec, res, dens):
    """Per regime: the selected rung against the true one, and the learned k against the true k on the domain."""
    m, x = res["model"], res["x"]
    rows, cells = [], []
    for regime, key in (("wetting", "wet"), ("drying", "dry")):
        g = dens[dens.regime == regime].copy()
        g["R"] = np.where(g.rcol == 0, 0.0, np.sqrt(np.r_[EDGES[0], EDGES][g.rcol.to_numpy()] * np.r_[EDGES[0], EDGES][np.minimum(g.rcol.to_numpy() + 1, len(EDGES))]))
        g["W"] = np.maximum(g.W_lo + 0.15, 0.0)
        g["k_fit"] = m.k(x, regime, g.R.to_numpy(), g.W.to_numpy())
        g["k_true"] = fake6.k_true(spec, regime, g.R.to_numpy(), g.W.to_numpy())
        g["ratio"] = g.k_fit / g.k_true
        g["inside"] = (g.ratio - 1).abs() <= J["within"]
        visited = g[g.minutes > 0]
        share = float(visited.minutes[visited.inside].sum() / visited.minutes.sum()) if len(visited) and np.isfinite(visited.ratio).all() else np.nan
        tr = fake6.true_rung(spec, regime)
        k0_fit, k0_true = float(m.k(x, regime, 0.0, 0.0)), float(fake6.k_true(spec, regime, 0.0, 0.0))
        const_ok = bool(abs(k0_fit / k0_true - 1) <= J["within"]) if tr == "K0" else None
        rows.append({"regime": regime, "true_rung": tr, "selected_rung": res[key]["rung"], "lowest_cv_rung": res[key]["lowest_cv"], "lambda": res[key]["lam"],
                     "rung_ok": res[key]["rung"] == tr, "tau_fit_at_R0_W0_min": 1 / k0_fit, "tau_true_at_R0_W0_min": 1 / k0_true if np.isfinite(k0_true) else np.nan,
                     "constant_within_25pct": const_ok, "share_within_25pct": share, "domain_cells_visited": int((g.minutes > 0).sum()),
                     "k_ok": (const_ok if tr == "K0" else bool(share >= J["min_share"])) if tr != "none" else None})
        cells.append(g)
    return pd.DataFrame(rows), pd.concat(cells, ignore_index=True)


# =================================================================================================== (e): the residual after rain
def after_rain(S, Dk, res):
    """Short links: the mean residual (fake A - c minus model) in the 30 min after the last wet bin, per link-event; the
    difference between the highest and the lowest third of preceding rain, plain and at matched model W (bins of 0.3 dB)."""
    D, le = S["D"], S["le"]
    AR = A6["after_rain"]
    W1, out5 = res["model"].simulate(res["x"])
    rain5 = fake6.window(D, Dk["rain"]).reshape(-1, 5).mean(axis=1)
    rows = []
    for e in np.nonzero(S["short"])[0]:
        o5, n5, last = D["le_boff5"][e], D["le_n"][e] // 5, D["le_last"][e]
        a, b = last + 1, min(last + 1 + AR["minutes"] // 5, n5)
        w = Dk["w5"][o5 + a:o5 + b] > 0
        if w.sum() < 3:
            continue
        rows.append({"le": e, "resid": float((Dk["y5"][o5 + a:o5 + b] - out5[o5 + a:o5 + b])[w].mean()), "rain_mm": float(rain5[o5:o5 + a].sum() / 12.0),
                     "W_end": float(W1[D["le_boff1"][e] + a * 5 - 1])})
    t = pd.DataFrame(rows)
    q1, q2 = t.rain_mm.quantile([1 / 3, 2 / 3])
    t["third"] = np.where(t.rain_mm <= q1, "low", np.where(t.rain_mm > q2, "high", "mid"))
    plain = float(t.resid[t.third == "high"].mean() - t.resid[t.third == "low"].mean())
    t["wbin"] = np.floor(t.W_end / (2 * AR["match_db"])).astype(int)
    num, den, nb = 0.0, 0.0, 0
    for _, g in t.groupby("wbin"):
        hi, lo = g.resid[g.third == "high"], g.resid[g.third == "low"]
        if len(hi) and len(lo):
            w = min(len(hi), len(lo))
            num, den, nb = num + w * (hi.mean() - lo.mean()), den + w, nb + 1
    return {"link_events": len(t), "high_minus_low_db": plain, "matched_high_minus_low_db": num / den if den else np.nan, "matched_bins": nb,
            "mean_resid_db": float(t.resid.mean())}


# =================================================================================================== (f): the flag
def tail_flags(S, Dk):
    """Per event: how many short-link link-events qualify, how many of them wait or rise, flagged or not (decisions.md, A6)."""
    D, le, F = S["D"], S["le"], A6["flag"]
    rows = []
    for e in np.nonzero(S["short"])[0]:
        o, n = D["le_boff1"][e], D["le_n"][e]
        end = (D["le_last"][e] + 1) * 5
        if end + 20 > n or end < 5:
            continue
        y = Dk["y1"][o:o + n]
        a0, a1 = y[end - 5:end + 5].mean(), y[end + 10:end + 20].mean()
        rows.append({"event_id": le.event_id.iloc[e], "qualifies": a0 >= F["a0_min_db"], "waits": (a0 >= F["a0_min_db"]) and (a1 >= a0 - F["fall_db"])})
    t = pd.DataFrame(rows).groupby("event_id").agg(link_events=("qualifies", "size"), qualifying=("qualifies", "sum"), waiting_or_rising=("waits", "sum")).reset_index()
    t["flagged"] = (t.qualifying >= F["min_link_events"]) & (t.waiting_or_rising >= F["share"] * t.qualifying)
    return t


# =================================================================================================== one truth
def extension_table(S, res, start_F, start_P):
    m, x = res["model"], res["x"]
    iF_fit, iP_fit = m.parts(x)[4:]
    E = A6["extension"]
    rows = []
    for b, c in enumerate(S["curves"]):
        for kind, knots, inc_true, inc_fit, inc_start, cap in (("F", A5["f_knots_mmh"], S["iF"][b], iF_fit[b], start_F[b], A5["f_inc_max_db"]),
                                                               ("P", A5["p_knots_mmh"], S["iP"][b], iP_fit[b], start_P[b], A5["p_inc_max_db_km"])):
            vt, vf, vs = (ksurf.full_curve(c, kind, i, knots) for i in (inc_true, inc_fit, inc_start))
            for j, r in enumerate(knots):
                at_bound = bool(inc_fit[j] >= cap - 1e-6 or (inc_fit[j] <= 1e-6 and inc_true[j] > E["at_zero_truth_above"]))
                rows.append({"band": c.band, "quantity": kind, "knot_mmh": r, "true": vt[j], "start": vs[j], "fitted": vf[j], "fitted_over_true": vf[j] / vt[j],
                             "within_25pct": bool(abs(vf[j] / vt[j] - 1) <= J["within"]), "increment_true": inc_true[j], "increment_fitted": inc_fit[j], "at_bound": at_bound})
    return pd.DataFrame(rows)


def slim(res):
    """What is kept of a two-stage fit (no arrays of the dataset)."""
    return {"cv": pd.concat([res["wet"]["table"], res["dry"]["table"]], ignore_index=True)}


def run_truth(name, S):
    spec = A6["truths"][name]
    seed = A6["seed"] + ORDER.index(name)
    D, le = S["D"], S["le"]
    t0 = time.time()
    Dk, ex = fake6.make(D, le, S["n_cells"], spec, S["curves"], S["iF"], S["iP"], seed, S["step"], S["dist"])
    bi = ex["band_idx"]
    out = {"name": name, "spec": spec, "seed": seed}
    free = bool(spec.get("free_extension", False))
    E = A6["extension"]
    start_F = [np.full(len(i), E["f_start_inc_db"]) for i in S["iF"]] if free else S["iF"]
    start_P = [E["p_start_share"] * i for i in S["iP"]] if free else S["iP"]
    tag = f"({name}) "
    if name == "a_lag":                                         # the rule for d on a truth with d = 2
        d_r, tab_r = ladder.lag_rule(Dk, le, bi, S["curves"], False, np.ones(len(le)), S["iF"], S["iP"], True, log, tag + "radar ")
        d_g, tab_g = ladder.lag_rule(Dk, le, bi, S["curves"], True, S["gauge"], S["iF"], S["iP"], True, log, tag + "gauge ")
        out.update({"lag": pd.concat([tab_r.assign(input="radar"), tab_g.assign(input="gauge")], ignore_index=True), "d_radar": d_r, "d_gauge": d_g})
        out["minutes"] = (time.time() - t0) / 60
        return out
    res = ladder.two_stage(Dk, le, bi, S["curves"], S["gauge"], inc_F=start_F, inc_P=start_P, hold_ext=not free, log=log, tag=tag)
    dens = state_density(S, Dk, ex)
    out.update(slim(res))
    out["after_rain"] = after_rain(S, Dk, res)
    if spec["kind"] != "two_store":
        out["judge"], out["cells"] = judge(S, spec, res, dens)
    else:
        m, x = res["model"], res["x"]
        out["judge"] = pd.DataFrame([{"regime": r, "true_rung": "none (two stores)", "selected_rung": res[k]["rung"], "lowest_cv_rung": res[k]["lowest_cv"], "lambda": res[k]["lam"],
                                      "tau_fit_at_R0_W0_min": float(1 / m.k(x, r, 0.0, 0.0))} for r, k in (("wetting", "wet"), ("drying", "dry"))])
        g = dens.copy()
        g["R"] = np.where(g.rcol == 0, 0.0, np.sqrt(np.r_[EDGES[0], EDGES][g.rcol.to_numpy()] * np.r_[EDGES[0], EDGES][np.minimum(g.rcol.to_numpy() + 1, len(EDGES))]))
        g["W"] = np.maximum(g.W_lo + 0.15, 0.0)
        g["k_fit"] = [float(m.k(x, r.regime, r.R, r.W)) for r in g.itertuples()]
        out["cells"] = g
    if free:
        out["extension"] = extension_table(S, res, start_F, start_P)
    if "ghost" in spec:                                         # exclude the flagged events and fit again
        fl = tail_flags(S, Dk)
        fl["changed_in_truth"] = fl.event_id.isin(ex["ghost_events"])
        out["flags"] = fl
        use = (~le.event_id.isin(set(fl.event_id[fl.flagged]))).to_numpy().astype(float)
        res2 = ladder.two_stage(Dk, le, bi, S["curves"], S["gauge"], inc_F=start_F, inc_P=start_P, hold_ext=True, use=use, log=log, tag=tag + "excl. ")
        out["judge_excluded"], out["cells_excluded"] = judge(S, spec, res2, dens)
        out["cv_excluded"] = slim(res2)["cv"]
    out["minutes"] = (time.time() - t0) / 60
    return out


def verdict(name, r):
    """(verdict, text) of one truth by the rules of decisions.md, step A6."""
    if name == "a_lag":
        return "report", f"d by the rule: radar {r['d_radar']}, gauge {r['d_gauge']} (truth 2)"
    j = r["judge"]
    if name == "e":
        return "report", "one-state ladder on a two-store truth"
    ok = bool(j.rung_ok.all() and j.k_ok.all())
    if name in GATES:
        return ("PASS" if ok else "FAIL"), ""
    if name == "b":
        e = r["extension"]
        f = e[e.quantity == "F"]
        ident = bool(f.within_25pct.all() and not f.at_bound.any() and j.k_ok.all())
        return ("identified" if ident else "NOT identified"), ""
    if name == "f":
        j2 = r["judge_excluded"]
        d1, d2 = j[j.regime == "drying"].iloc[0], j2[j2.regime == "drying"].iloc[0]
        spur = not (d1.rung_ok and d1.k_ok)
        rest = bool(d2.rung_ok and d2.k_ok)
        return ("spurious k_dry" if spur else "no spurious k_dry") + (", restored by the exclusion" if rest else ", NOT restored by the exclusion"), ""
    if name == "g":
        return ("k corrupted" if not ok else "k not corrupted"), ""
    return "", ""


# =================================================================================================== report
def fig_surfaces(R):
    names = [n for n in ("c", "d") if n in R]
    if not names:
        return None
    fig, axes = plt.subplots(len(names), 2, figsize=(11, 4 * len(names)), squeeze=False)
    for i, n in enumerate(names):
        c = R[n]["cells"]
        for j, regime in enumerate(("wetting", "drying")):
            ax = axes[i, j]
            g = c[(c.regime == regime) & (c.minutes > 0)].sort_values(["R", "W"])
            xs, lab = (g.W, "W (dB)") if n == "c" else (g.R, "R (mm/h)")
            ax.plot(xs, 1 / g.k_true, "o", color=INK, mfc="white", ms=5, label="truth (1 / k)")
            ax.plot(xs, 1 / g.k_fit, "x", color=INK2, ms=6, label="fitted")
            ax.set_yscale("log")
            if n == "d":
                ax.set_xscale("symlog", linthresh=0.1)
            ax.set_xlabel(f"{lab} of the domain cell")
            ax.set_ylabel("1 / k (min)")
            ax.set_title(f"truth ({n}), {regime}: selected {R[n]['judge'].set_index('regime').selected_rung[regime]}", fontsize=9)
            ax.legend(fontsize=7)
            style_axes(ax)
    fig.suptitle("A6: the learned rate against the true one, at the centres of the visited domain cells (circles = truth, crosses = fit)", fontsize=10)
    fig.tight_layout()
    return fig


def report(R):
    rows, stop = [], None
    for n in ORDER:
        if n not in R:
            continue
        v, txt = verdict(n, R[n])
        if n == "a_lag":
            rows.append({"truth": "(a-lag)", "regime": "", "true rung": "K0", "rung returned": "", "k within 25 %": txt, "verdict": v})
            continue
        for r in R[n]["judge"].itertuples():
            share = getattr(r, "share_within_25pct", np.nan)
            const = getattr(r, "constant_within_25pct", None)
            rows.append({"truth": f"({n})", "regime": r.regime, "true rung": r.true_rung, "rung returned": r.selected_rung,
                         "k within 25 %": (f"1/k {r.tau_fit_at_R0_W0_min:.1f} min (truth {r.tau_true_at_R0_W0_min:.1f}): {'yes' if const else 'no'}" if const is not None and r.true_rung == "K0"
                                           else (f"{100 * share:.0f} % of the domain weight" if np.isfinite(share) else f"1/k {r.tau_fit_at_R0_W0_min:.1f} min at R = 0, W = 0")), "verdict": v})
        if n == "f":
            for r in R[n]["judge_excluded"].itertuples():
                rows.append({"truth": "(f) flagged events excluded", "regime": r.regime, "true rung": r.true_rung, "rung returned": r.selected_rung,
                             "k within 25 %": f"1/k {r.tau_fit_at_R0_W0_min:.1f} min (truth {r.tau_true_at_R0_W0_min:.1f}): {'yes' if r.constant_within_25pct else 'no'}", "verdict": v})
        if n in GATES and v == "FAIL" and stop is None:
            stop = n
    table = pd.DataFrame(rows)
    save(table, P.A6_VERDICTS, index=False)
    cv = pd.concat([R[n]["cv"].assign(truth=n) for n in ORDER if n in R and "cv" in R[n]] + [R["f"]["cv_excluded"].assign(truth="f_excluded")] if "f" in R else
                   [R[n]["cv"].assign(truth=n) for n in ORDER if n in R and "cv" in R[n]], ignore_index=True)
    save(cv.round(6), P.A6_CV, index=False)
    cells = pd.concat([R[n]["cells"].assign(truth=n) for n in ORDER if n in R and "cells" in R[n]], ignore_index=True)
    save(cells.round(5), P.A6_CELLS, index=False)
    fig = fig_surfaces(R)
    if fig is not None:
        save(fig, P.A6_FIG, dpi=120)
        plt.close(fig)

    a = ["# A6: fake-data recovery (the gates for the real fit)\n",
         "Written by `research/waa_ksurface/a6_recovery.py` (libraries `fake6.py`, `ladder.py`, `ksurf.py`). Fake attenuation on the 4 563 TRAIN link-events of the learning set "
         "(real links, real rain), 0.3 dB quantization, the measured noise, the T1 rain error at its high end for the radar input and the T2 gauge model for the gauge input. "
         "One replicate per truth. No real attenuation is used as a target. VAL and TEST were not opened. Rules: `docs/waa_ksurface/decisions.md`, step A6.\n",
         "## The table\n", md_table(table),
         f"\n(`{rel(P.A6_VERDICTS)}`; every rung and λ with its held-out loss: `{rel(P.A6_CV)}`; every domain cell: `{rel(P.A6_CELLS)}`"
         + (f"; figure `{rel(P.A6_FIG)}`" if fig is not None else "") + ")\n"]
    if stop:
        a.append(f"**STOP: truth ({stop}) FAILS its rule. The truths after it were not run. Nothing was patched.**\n")
    for n in ORDER:
        if n not in R or n == "a_lag":
            continue
        a.append(f"## Truth ({n}): held-out loss per rung (λ of the lowest held-out loss)\n")
        c = R[n]["cv"]
        a.append(md_table(c[c.best_lambda][["select", "rung", "lambda", "parameters_free", "train_loss", "cv_loss", "minus_best_lo", "minus_best_hi", "lowest_cv", "selected"]].round(6)))
        a.append("\nPer regime:\n")
        a.append(md_table(R[n]["judge"].round(3)))
        ar = R[n]["after_rain"]
        a.append(f"\nResidual in the 30 min after rain on the short links ({ar['link_events']} link-events): mean {ar['mean_resid_db']:+.3f} dB; highest minus lowest third of preceding rain "
                 f"{ar['high_minus_low_db']:+.3f} dB, at matched model W {ar['matched_high_minus_low_db']:+.3f} dB ({ar['matched_bins']} W bins).\n")
        if "extension" in R[n]:
            save(R[n]["extension"].round(5), P.A6_EXTENSION, index=False)
            a.append("The extension (fitted with k; started from F held and P at half its true increments):\n")
            a.append(md_table(R[n]["extension"].round(3)))
            a.append(f"\n(`{rel(P.A6_EXTENSION)}`)\n")
        if "flags" in R[n]:
            save(R[n]["flags"], P.A6_FLAGS, index=False)
            a.append("The tail flag per event (`changed_in_truth` = the events that were given unseen rain):\n")
            a.append(md_table(R[n]["flags"]))
            a.append("\nWith the flagged events excluded:\n")
            a.append(md_table(R[n]["judge_excluded"].round(3)))
            a.append(f"\n(`{rel(P.A6_FLAGS)}`)\n")
    if "a_lag" in R:
        save(R["a_lag"]["lag"].round(6), P.A6_LAG, index=False)
        a.append("## Truth (a-lag): the rule for d on a truth with d = 2 min\n")
        a.append(md_table(R["a_lag"]["lag"].round(5)))
        a.append(f"\nChosen: radar d = {R['a_lag']['d_radar']}, gauge d = {R['a_lag']['d_gauge']}. (`{rel(P.A6_LAG)}`)\n")
    a.append("## Run times\n")
    a.append(md_table(pd.DataFrame([{"truth": n, "minutes": round(R[n]["minutes"], 1)} for n in ORDER if n in R])))
    save("\n".join(a) + "\n", P.A6_REPORT)
    log(f"report: {rel(P.A6_REPORT)}")
    return stop


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    force, only_report = "--force" in sys.argv, "--report" in sys.argv
    cache = P.CACHE / "a6_recovery"
    cache.mkdir(parents=True, exist_ok=True)
    R, S = {}, None
    for n in ORDER:
        f = cache / f"{n}.pkl"
        wanted = (not args) or (n in args)
        if f.exists() and not (force and wanted):
            R[n] = pickle.loads(f.read_bytes())
        elif wanted and not only_report:
            S = S or setting()
            log(f"=== truth ({n}) ===")
            R[n] = run_truth(n, S)
            f.write_bytes(pickle.dumps(R[n]))
            log(f"=== truth ({n}): {verdict(n, R[n])[0]} ({R[n]['minutes']:.1f} min) ===")
        if n in R and n in GATES and verdict(n, R[n])[0] == "FAIL":
            break                                               # the rule: a FAIL on a gate is a STOP
    stop = report(R)
    if stop:
        sys.exit(1)


if __name__ == "__main__":
    main()
