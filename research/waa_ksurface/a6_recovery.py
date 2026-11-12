"""Step A6 of the k-surface branch: does the fit find a KNOWN k on fake data? The gates for the real fit (A7).

    python -m research.waa_ksurface.a6_recovery                 every truth in order (a finished truth is read from the cache)
    python -m research.waa_ksurface.a6_recovery c d             only these truths
    python -m research.waa_ksurface.a6_recovery --force         run again, ignoring the cache
    python -m research.waa_ksurface.a6_recovery --report        write the report from the cache

Fake attenuation on the real links and the real rain (library `fake6.py`), fitted with the procedure of A5 (library
`ladder.py`): stage 1 learns k_wet from the gauge-input set, stage 2 learns k_dry from all links. Truths, PASS rules and
what a failure means: decisions.md, step A6 and "A6 again under Deviations 5 and 6" (written before any run).
Order: the k_wet sensitivity gate (Deviation 6), truth (a), then (a-lag), (b)-(g), then the detectability sweep. A failed
gate, or a FAIL on (a), (c) or (d), is a STOP: the script writes the report and exits with code 1.
No real attenuation is used as a target.
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
DET = [f"{t}_{a:g}" for t in ("c", "d") for a in A6["detect"]["grid"]]      # the detectability sweep: truth form and size a
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
    g = S["dom"].merge(cnt, on=["regime", "rcol", "wbin"], how="left").fillna({"minutes": 0})
    e = np.r_[EDGES[0], EDGES]
    g["R"] = np.where(g.rcol == 0, 0.0, np.sqrt(e[g.rcol.to_numpy()] * e[np.minimum(g.rcol.to_numpy() + 1, len(EDGES))]))
    g["W"] = np.maximum(g.W_lo + 0.15, 0.0)
    return g


def compare(m, x, spec, regime, dens):
    """The learned k of one regime against the truth: (cells, share of the domain weight within 25 %, 1 / k at R = 0, W = 0)."""
    g = dens[dens.regime == regime].copy()
    g["k_fit"] = m.k(x, regime, g.R.to_numpy(), g.W.to_numpy())
    g["k_true"] = fake6.k_true(spec, regime, g.R.to_numpy(), g.W.to_numpy())
    g["ratio"] = g.k_fit / g.k_true
    v = g[g.minutes > 0]
    share = float(v.minutes[(v.ratio - 1).abs() <= J["within"]].sum() / v.minutes.sum()) if len(v) and np.isfinite(v.ratio).all() else np.nan
    return g, share, float(1 / m.k(x, regime, 0.0, 0.0))


def judge(S, spec, res, dens, gate_ok):
    """Per regime: the rung returned against the true one, and the learned k against the true k (decisions.md, A6 and
    Deviation 6). Drying is judged with k_wet at its gauge-input estimate and, if k_wet is a nuisance, also at its
    radar-input estimate: both must hold."""
    m, x = res["model"], res["x"]
    rows, cells = [], []
    for regime, key in (("wetting", "wet"), ("drying", "dry")):
        tr = fake6.true_rung(spec, regime)
        g, share, tau = compare(m, x, spec, regime, dens)
        tau_true = float(1 / fake6.k_true(spec, regime, 0.0, 0.0)) if tr != "none" else np.nan
        row = {"regime": regime, "true_rung": tr, "selected_rung": res[key]["rung"], "lowest_cv_rung": res[key]["lowest_cv"], "lambda": res[key]["lam"],
               "rung_ok": res[key]["rung"] == tr, "tau_fit_min": tau, "tau_true_min": tau_true, "share_within_25pct": share,
               "domain_cells_visited": int((g.minutes > 0).sum()), "tau_fit_alt_min": np.nan, "share_alt": np.nan, "k_note": ""}
        ok = (abs(tau_true / tau - 1) <= J["within"]) if tr == "K0" else (share >= J["min_share"] if tr != "none" else None)
        if regime == "wetting" and "tau_wet" in res:            # the wetting ladder returned K0: k_wet is the nuisance interval
            row["tau_fit_alt_min"] = res["tau_wet"]["radar"]
            if tr == "K0" and not ok and gate_ok:
                ok, row["k_note"] = True, "governed by the gate"
        if regime == "drying" and "x_alt" in res:
            g2, share2, tau2 = compare(m, res["x_alt"], spec, regime, dens)
            row["tau_fit_alt_min"], row["share_alt"] = tau2, share2
            g["k_fit_alt"] = g2.k_fit
            if tr != "none":
                ok = bool(ok and ((abs(tau_true / tau2 - 1) <= J["within"]) if tr == "K0" else share2 >= J["min_share"]))
        row["k_ok"] = None if ok is None else bool(ok)
        rows.append(row)
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


# =================================================================================================== the k_wet gate (Deviation 6)
def wet_gate(S):
    """Truth (a), all links, radar input, K0 / K0 with 1 / k_wet HELD at half the truth, the truth, twice the truth; k_dry and
    the extension free. Per value: the grouped-CV error (dB) and 1 / k_dry."""
    spec, G = A6["truths"]["a"], A6["gate"]
    le = S["le"]
    Dk, ex = fake6.make(S["D"], le, S["n_cells"], spec, S["curves"], S["iF"], S["iP"], A6["seed"] + ORDER.index("a"), S["step"], S["dist"])
    prob = ksurf.prepare(Dk, ex["band_idx"], S["curves"], Dk["rain"], 0)
    E = A6["extension"]
    start_F, start_P = [np.full(len(i), E["f_start_inc_db"]) for i in S["iF"]], [E["p_start_share"] * i for i in S["iP"]]
    use = np.ones(len(le))
    fold = ladder.folds(le, use)
    ev, n_ev = prob["event"], int(prob["event"].max()) + 1
    rows = []
    for tau in G["tau_wet_min"]:
        m = ksurf.Model(prob, "K0", "K0")
        x0 = m.start(start_F, start_P)
        x0[m.sl["gw"]] = np.log(1.0 / tau)
        m.hold("gw")
        x = m.fit(x0, use)
        Ssum, N = np.zeros(n_ev), np.zeros(n_ev)
        for f in range(fold.max() + 1):
            xf = m.fit(x, use * (fold != f), maxiter=A5["maxiter_fold"])
            ls, cnt = m.per_link_event(xf)
            Ssum += np.bincount(ev[fold == f], weights=ls[fold == f], minlength=n_ev)
            N += np.bincount(ev[fold == f], weights=cnt[fold == f], minlength=n_ev)
        tau_d = float(1 / m.k(x, "drying", 0, 0))
        rows.append({"tau_wet_held_min": tau, "cv_error_db": ladder.cv_db(Ssum, N), "train_loss": m.last["loss"], "tau_dry_fit_min": tau_d,
                     "tau_dry_true_min": 1 / spec["k_dry"], "k_dry_within_25pct": bool(abs((1 / spec["k_dry"]) / tau_d - 1) <= J["within"])})
        log(f"   gate: 1/k_wet held at {tau:g} min: cv {rows[-1]['cv_error_db']:.5f} dB, 1/k_dry {tau_d:.2f} min")
    t = pd.DataFrame(rows)
    change = float(t.cv_error_db.max() - t.cv_error_db.min())
    ok = bool(change < G["margin_db"] and t.k_dry_within_25pct.all())
    return {"table": t, "cv_change_db": change, "pass": ok}


# =================================================================================================== one truth
def extension_table(S, m, x, start_F, start_P):
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


def spec_of(name):
    """The truth of a name: one of the written truths, or a point of the detectability sweep ("c_0.25": form (c), size a)."""
    if name in A6["truths"]:
        return A6["truths"][name], A6["seed"] + ORDER.index(name)
    form, a = name.split("_")
    a = float(a)
    seed = A6["seed"] + A6["detect"]["seed_offset"] + DET.index(name)
    if form == "c":                                             # k_dry = (1/15) (1 + a W) / (1 + a): 15 min at W = 1 dB for every a
        return {"kind": "dry_w", "k_wet": 0.2, "k_dry0": (1 / 15) / (1 + a), "w_scale_db": 1 / a, "size_a": a}, seed
    return {"kind": "both_r", "k0": 0.1, "r_scale_mmh": 2.0 / a, "size_a": a}, seed


def run_truth(name, S, gate_ok):
    spec, seed = spec_of(name)
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
    out["cv"] = pd.concat([res["wet"]["table"], res["dry"]["table"]], ignore_index=True)
    out["tau_wet"] = res.get("tau_wet")
    out["after_rain"] = after_rain(S, Dk, res)
    out["judge"], out["cells"] = judge(S, spec, res, dens, gate_ok)
    if free:
        out["extension"] = extension_table(S, res["model"], res["x"], start_F, start_P)
        if "x_alt" in res:
            out["extension_alt"] = extension_table(S, res["model"], res["x_alt"], start_F, start_P)
    if "ghost" in spec:                                         # exclude the flagged events and fit again
        fl = tail_flags(S, Dk)
        fl["changed_in_truth"] = fl.event_id.isin(ex["ghost_events"])
        out["flags"] = fl
        use = (~le.event_id.isin(set(fl.event_id[fl.flagged]))).to_numpy().astype(float)
        res2 = ladder.two_stage(Dk, le, bi, S["curves"], S["gauge"], inc_F=start_F, inc_P=start_P, hold_ext=True, use=use, log=log, tag=tag + "excl. ")
        out["judge_excluded"], out["cells_excluded"] = judge(S, spec, res2, dens, gate_ok)
        out["cv_excluded"] = pd.concat([res2["wet"]["table"], res2["dry"]["table"]], ignore_index=True)
    out["minutes"] = (time.time() - t0) / 60
    return out


def verdict(name, r):
    """The verdict of one truth by the rules of decisions.md (step A6, and A6 again under Deviations 5 and 6)."""
    if name in ("a_lag", "e"):
        return "report"
    j = r["judge"]
    ok = bool(j.rung_ok.all() and j.k_ok.all())
    if name in GATES:
        return "PASS" if ok else "FAIL"
    if name == "b":
        ident = True
        for key in ("extension", "extension_alt"):
            if key in r:
                f = r[key][r[key].quantity == "F"]
                ident &= bool(f.within_25pct.all() and not f.at_bound.any())
        return "identified" if (ident and bool(j.k_ok.all())) else "NOT identified"
    if name == "f":
        d1 = j[j.regime == "drying"].iloc[0]
        d2 = r["judge_excluded"][r["judge_excluded"].regime == "drying"].iloc[0]
        spur, rest = not (d1.rung_ok and d1.k_ok), bool(d2.rung_ok and d2.k_ok)
        return ("spurious k_dry" if spur else "no spurious k_dry") + (", restored by the exclusion" if rest else ", NOT restored by the exclusion")
    if name == "g":
        return "k not corrupted" if ok else "k corrupted"
    return "detectability"


# =================================================================================================== report
def k_text(r):
    if r.true_rung == "none":
        return f"1/k {r.tau_fit_min:.1f} min at R = 0, W = 0"
    alt = f" / {r.tau_fit_alt_min:.1f}" if np.isfinite(r.tau_fit_alt_min) else ""
    if r.true_rung == "K0":
        return f"1/k {r.tau_fit_min:.1f}{alt} min (truth {r.tau_true_min:.1f})" + (f": {r.k_note}" if r.k_note else (": yes" if r.k_ok else ": no"))
    sh = f"{100 * r.share_within_25pct:.0f} %" + (f" / {100 * r.share_alt:.0f} %" if np.isfinite(r.share_alt) else "")
    return f"{sh} of the domain weight" + (": yes" if r.k_ok else ": no")


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


def cv_show(c):
    cols = ["select", "rung", "lambda", "parameters_free", "cv_db", "compared_with", "gain_db", "gain_lo", "gain_hi", "replaces", "lowest_cv", "selected"]
    t = c[c.best_lambda][cols].copy()
    for k in ("cv_db", "gain_db", "gain_lo", "gain_hi"):
        t[k] = pd.to_numeric(t[k]).round(5)
    return md_table(t)


def gain_of(cv, regime, rung):
    g = cv[(cv.select == regime) & (cv.rung == rung) & cv.best_lambda]
    return round(float(g.gain_db.iloc[0]), 5) if len(g) and g.gain_db.iloc[0] is not None else np.nan


def detect_table(R):
    rows = []
    for n in DET + ["c", "d"]:
        if n not in R:
            continue
        form = n.split("_")[0]
        a = R[n]["spec"].get("size_a", 1.0)
        j = R[n]["judge"].set_index("regime")
        for regime in (("drying",) if form == "c" else ("wetting", "drying")):
            ratio = (1 + 2 * a) if form == "c" else ((1 + 2.5 * a) if regime == "wetting" else 1 / (1 + 2.5 * a))
            rows.append({"form": "k_dry(W), truth (c)" if form == "c" else f"k(R), truth (d), {regime}", "size_a": a,
                         "k_ratio_across_domain": round(ratio, 3), "true_rung": j.true_rung[regime], "rung_returned": j.selected_rung[regime],
                         "detected": bool(j.rung_ok[regime]), "gain_db_of_true_rung": gain_of(R[n]["cv"], regime, j.true_rung[regime])})
    return pd.DataFrame(rows).sort_values(["form", "size_a"]).reset_index(drop=True) if rows else pd.DataFrame()


def report(R, gate):
    rows, stop = [], None
    if gate is not None and not gate["pass"]:
        stop = "the k_wet gate"
    for n in ORDER:
        if n not in R:
            continue
        v = verdict(n, R[n])
        if n == "a_lag":
            rows.append({"truth": "(a-lag)", "regime": "", "true rung": "K0", "rung returned": "", "recovery": f"d by the rule: radar {R[n]['d_radar']}, gauge {R[n]['d_gauge']} (truth 2)", "verdict": v})
            continue
        for r in R[n]["judge"].itertuples():
            rows.append({"truth": f"({n})", "regime": r.regime, "true rung": r.true_rung, "rung returned": r.selected_rung, "recovery": k_text(r), "verdict": v})
        if n == "f":
            for r in R[n]["judge_excluded"].itertuples():
                rows.append({"truth": "(f), flagged events excluded", "regime": r.regime, "true rung": r.true_rung, "rung returned": r.selected_rung, "recovery": k_text(r), "verdict": v})
        if n in GATES and v == "FAIL" and stop is None:
            stop = f"truth ({n})"
    table = pd.DataFrame(rows)
    save(table, P.A6_VERDICTS, index=False)
    parts = [R[n]["cv"].assign(truth=n) for n in ORDER + DET if n in R and "cv" in R[n]] + ([R["f"]["cv_excluded"].assign(truth="f_excluded")] if "f" in R else [])
    if parts:
        save(pd.concat(parts, ignore_index=True), P.A6_CV, index=False)
        save(pd.concat([R[n]["cells"].assign(truth=n) for n in ORDER + DET if n in R and "cells" in R[n]], ignore_index=True).round(5), P.A6_CELLS, index=False)
    fig = fig_surfaces(R)
    if fig is not None:
        save(fig, P.A6_FIG, dpi=120)
        plt.close(fig)
    a = ["# A6: fake-data recovery (the gates for the real fit), under Deviations 5 and 6\n",
         "Written by `research/waa_ksurface/a6_recovery.py` (libraries `fake6.py`, `ladder.py`, `ksurf.py`). Fake attenuation on the 4 563 TRAIN link-events of the learning set "
         "(real links, real rain), 0.3 dB quantization, the measured noise, the T1 rain error at its high end for the radar input and the T2 gauge model for the gauge input. "
         "One replicate per truth. Selection: a richer rung must lower the CV error by more than 0.005 dB with a paired 90 % interval that excludes zero (Deviation 5). "
         "k_wet: a bounded nuisance if the gate passes (Deviation 6); where two numbers are given as `x / y`, the first is with k_wet at its gauge-input estimate, the second "
         "at its radar-input estimate. No real attenuation is used as a target. VAL and TEST were not opened.\n"]
    if gate is not None:
        save(gate["table"].round(6), P.A6_GATE, index=False)
        a += ["## The k_wet sensitivity gate (Deviation 6)\n",
              "Truth (a): 1 / k_wet = 5 min, 1 / k_dry = 15 min. All links, radar input, K0 in both regimes, 1 / k_wet HELD; k_dry and the F / P extension free.\n",
              md_table(gate["table"].round(5)),
              f"\nCV error changes by **{gate['cv_change_db']:.5f} dB** across the three values (margin {A6['gate']['margin_db']} dB). "
              f"**Gate: {'PASS: k_wet is a bounded nuisance' if gate['pass'] else 'FAIL: k_wet matters. STOP'}.** (`{rel(P.A6_GATE)}`)\n"]
    a += ["## The table\n", md_table(table) if len(table) else "no truth was run",
          f"\n(`{rel(P.A6_VERDICTS)}`; every rung and lambda: `{rel(P.A6_CV)}`; every domain cell: `{rel(P.A6_CELLS)}`" + (f"; figure `{rel(P.A6_FIG)}`" if fig is not None else "") + ")\n"]
    if stop:
        a.append(f"**STOP: {stop} FAILS its rule. What comes after it was not run. Nothing was patched.**\n")
    det = detect_table(R)
    if len(det):
        save(det, P.A6_DETECT, index=False)
        a += ["## Detectability: the smallest dependence the ladder returns at the 0.005 dB margin\n",
              "The dependence of truths (c) and (d) scaled by a (a = 1 is the written truth). `k_ratio_across_domain` = k(W = 2 dB) / k(W = 0) for k_dry(W), "
              "k(R = 5 mm/h) / k(R = 0) for k(R). `gain_db_of_true_rung` = how much the true rung lowered the CV error against the rung selected before it.\n",
              md_table(det), f"\n(`{rel(P.A6_DETECT)}`)\n"]
    for n in ORDER:
        if n not in R or n == "a_lag":
            continue
        a.append(f"## Truth ({n})\n")
        a.append(cv_show(R[n]["cv"]))
        a.append("\nPer regime:\n")
        a.append(md_table(R[n]["judge"].round(3)))
        if R[n].get("tau_wet"):
            a.append(f"\n1 / k_wet: gauge-input estimate {R[n]['tau_wet']['gauge']:.2f} min, radar-input estimate {R[n]['tau_wet']['radar']:.2f} min.")
        ar = R[n]["after_rain"]
        a.append(f"\nResidual in the 30 min after rain on the short links ({ar['link_events']} link-events): mean {ar['mean_resid_db']:+.3f} dB; highest minus lowest third of preceding rain "
                 f"{ar['high_minus_low_db']:+.3f} dB, at matched model W {ar['matched_high_minus_low_db']:+.3f} dB ({ar['matched_bins']} W bins).\n")
        if "extension" in R[n]:
            ext = R[n]["extension"].assign(k_wet_at="gauge estimate")
            if "extension_alt" in R[n]:
                ext = pd.concat([ext, R[n]["extension_alt"].assign(k_wet_at="radar estimate")], ignore_index=True)
            save(ext.round(5), P.A6_EXTENSION, index=False)
            a.append("The extension (fitted with k; started from F held and P at half its true increments):\n")
            a.append(md_table(ext.round(3)))
            a.append(f"\n(`{rel(P.A6_EXTENSION)}`)\n")
        if "flags" in R[n]:
            save(R[n]["flags"], P.A6_FLAGS, index=False)
            a.append("The tail flag per event (`changed_in_truth` = the events that were given unseen rain):\n")
            a.append(md_table(R[n]["flags"]))
            a.append("\nWith the flagged events excluded:\n")
            a.append(cv_show(R[n]["cv_excluded"]))
            a.append("")
            a.append(md_table(R[n]["judge_excluded"].round(3)))
            a.append(f"\n(`{rel(P.A6_FLAGS)}`)\n")
    if "a_lag" in R:
        save(R["a_lag"]["lag"].round(6), P.A6_LAG, index=False)
        a.append("## Truth (a-lag): the rule for d on a truth with d = 2 min\n")
        a.append(md_table(R["a_lag"]["lag"].round(5)))
        a.append(f"\nChosen: radar d = {R['a_lag']['d_radar']}, gauge d = {R['a_lag']['d_gauge']}. (`{rel(P.A6_LAG)}`)\n")
    a.append("## Run times\n")
    a.append(md_table(pd.DataFrame([{"truth": n, "minutes": round(R[n]["minutes"], 1)} for n in ORDER + DET if n in R])))
    save("\n".join(a) + "\n", P.A6_REPORT)
    log(f"report: {rel(P.A6_REPORT)}")
    return stop


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    force, only_report = "--force" in sys.argv, "--report" in sys.argv
    cache = P.CACHE / "a6_recovery_dev56"
    cache.mkdir(parents=True, exist_ok=True)
    R, S, gate = {}, None, None
    fg = cache / "gate.pkl"
    if fg.exists() and not (force and not args):
        gate = pickle.loads(fg.read_bytes())
    elif not only_report:
        S = setting()
        log("=== the k_wet gate ===")
        gate = wet_gate(S)
        fg.write_bytes(pickle.dumps(gate))
        log(f"=== gate: {'PASS' if gate['pass'] else 'FAIL'} (cv change {gate['cv_change_db']:.5f} dB) ===")
    if gate is not None and gate["pass"]:
        for n in ORDER + DET:
            f = cache / f"{n}.pkl"
            wanted = (not args) or (n in args)
            if f.exists() and not (force and wanted):
                R[n] = pickle.loads(f.read_bytes())
            elif wanted and not only_report:
                S = S or setting()
                log(f"=== truth ({n}) ===")
                R[n] = run_truth(n, S, gate["pass"])
                f.write_bytes(pickle.dumps(R[n]))
                log(f"=== truth ({n}): {verdict(n, R[n])} ({R[n]['minutes']:.1f} min) ===")
            if n in R and n in GATES and verdict(n, R[n]) == "FAIL":
                break                                           # the rule: a FAIL on a gate is a STOP
    stop = report(R, gate)
    if stop:
        sys.exit(1)


if __name__ == "__main__":
    main()
