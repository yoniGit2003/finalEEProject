"""Step A6 of the k-surface branch: does the fit find a KNOWN k on fake data? The gates for the real fit (A7).

    python -m research.waa_ksurface.a6_recovery                 every truth in order (a finished truth is read from the cache)
    python -m research.waa_ksurface.a6_recovery c d             only these truths
    python -m research.waa_ksurface.a6_recovery --force         run again, ignoring the cache
    python -m research.waa_ksurface.a6_recovery --report        write the report from the cache

Fake attenuation on the real links and the real rain (library `fake6.py`), fitted with the procedure of A5 (library
`ladder.py`): stage 1 learns k_wet from the gauge-input set, stage 2 learns k_dry from all links. Truths, PASS rules and
what a failure means: decisions.md, step A6 and "A6 again under Deviations 5 and 6" (written before any run).
Final form (Deviation 7): every truth frees the F / P extension, as the real fit will, and is drawn with a fresh seed.
Order: truth (a), then (a-lag), (b), (c), (d), (d-high), (e), (f), (g), then the detectability study (5 replicates per size).
A FAIL on (a), (c) or (d) is a STOP: the script writes the report and exits with code 1.
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
ORDER = ["a", "a_lag", "b", "c", "d", "d_high", "e", "f", "g"]
FORMS = {"cw": ("drying",), "kr": ("wetting", "drying")}        # detectability forms: k_dry(W); k(R) inside the measured range
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


R2_LAST_COL = 5                                                 # the R column of the class 2-5 mm/h: the last one inside F's measured range (R2)


def _share(v):
    return float(v.minutes[(v.ratio - 1).abs() <= J["within"]].sum() / v.minutes.sum()) if len(v) and np.isfinite(v.ratio).all() else np.nan


def compare(m, x, spec, regime, dens):
    """The learned k of one regime against the truth: (cells, share of the domain weight within 25 %, 1 / k at R = 0, W = 0,
    the same share on the whole domain). For a truth in which k depends on R the share is taken inside F's measured range
    (restriction R2); for the others the two shares are the same."""
    g = dens[dens.regime == regime].copy()
    g["k_fit"] = m.k(x, regime, g.R.to_numpy(), g.W.to_numpy())
    g["k_true"] = fake6.k_true(spec, regime, g.R.to_numpy(), g.W.to_numpy())
    g["ratio"] = g.k_fit / g.k_true
    v = g[g.minutes > 0]
    whole = _share(v)
    share = _share(v[v.rcol <= R2_LAST_COL]) if spec["kind"] == "both_r" else whole
    return g, share, float(1 / m.k(x, regime, 0.0, 0.0)), whole


def judge(S, spec, res, dens, gate_ok):
    """Per regime: the rung returned against the true one, and the learned k against the true k (decisions.md, A6 and
    Deviation 6). Drying is judged with k_wet at its gauge-input estimate and, if k_wet is a nuisance, also at its
    radar-input estimate: both must hold."""
    m, x = res["model"], res["x"]
    rows, cells = [], []
    for regime, key in (("wetting", "wet"), ("drying", "dry")):
        tr = fake6.true_rung(spec, regime)
        g, share, tau, whole = compare(m, x, spec, regime, dens)
        tau_true = float(1 / fake6.k_true(spec, regime, 0.0, 0.0)) if tr != "none" else np.nan
        row = {"regime": regime, "true_rung": tr, "selected_rung": res[key]["rung"], "lowest_cv_rung": res[key]["lowest_cv"], "lambda": res[key]["lam"],
               "rung_ok": res[key]["rung"] == tr, "tau_fit_min": tau, "tau_true_min": tau_true, "share_within_25pct": share, "share_whole_domain": whole,
               "domain_cells_visited": int((g.minutes > 0).sum()), "tau_fit_alt_min": np.nan, "share_alt": np.nan, "k_note": ""}
        ok = (abs(tau_true / tau - 1) <= J["within"]) if tr == "K0" else (share >= J["min_share"] if tr != "none" else None)
        if regime == "wetting" and "tau_wet" in res:            # the wetting ladder returned K0: k_wet is the nuisance interval
            row["tau_fit_alt_min"] = res["tau_wet"]["radar"]
            if tr == "K0" and not ok and gate_ok:
                ok, row["k_note"] = True, "governed by the gate"
        if regime == "drying" and "x_alt" in res:
            g2, share2, tau2, _ = compare(m, res["x_alt"], spec, regime, dens)
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


def det_name(form, a, rep):
    return f"{form}_{a:g}_r{rep}"


def spec_of(name):
    """(truth, seed) of a name: a written truth, or a replicate of the detectability study ("cw_0.5_r2": form, size a, replicate)."""
    base = A6["seed_final"]
    if name in A6["truths"]:
        return A6["truths"][name], base + ORDER.index(name)
    form, a, rep = name.split("_")
    a, rep = float(a), int(rep[1:])
    seed = base + 100 + 1000 * list(FORMS).index(form) + 10 * A6["detect"]["sizes"].index(a) + rep
    if form == "cw":                                            # k_dry = (1/15) (1 + a W) / (1 + a): 15 min at W = 1 dB for every a
        return {"kind": "dry_w", "k_wet": 0.2, "k_dry0": (1 / 15) / (1 + a), "w_scale_db": 1 / a, "size_a": a}, seed
    return {"kind": "both_r", "k0": 0.1, "r_scale_mmh": 2.0 / a, "r_cap_mmh": A6["detect"]["r_cap_mmh"], "size_a": a}, seed


def run_truth(name, S, gate_ok=True):
    """One truth, fitted exactly as A7 will fit: extension free (Deviation 7), started from F held and P at half its increments."""
    spec, seed = spec_of(name)
    D, le = S["D"], S["le"]
    t0 = time.time()
    Dk, ex = fake6.make(D, le, S["n_cells"], spec, S["curves"], S["iF"], S["iP"], seed, S["step"], S["dist"])
    bi = ex["band_idx"]
    out = {"name": name, "spec": spec, "seed": seed}
    E = A6["extension"]
    start_F = [np.full(len(i), E["f_start_inc_db"]) for i in S["iF"]]
    start_P = [E["p_start_share"] * i for i in S["iP"]]
    tag = f"({name}) "
    if name == "a_lag":                                         # the rule for d on a truth with d = 2
        d_r, tab_r = ladder.lag_rule(Dk, le, bi, S["curves"], False, np.ones(len(le)), start_F, start_P, False, log, tag + "radar ")
        d_g, tab_g = ladder.lag_rule(Dk, le, bi, S["curves"], True, S["gauge"], start_F, start_P, False, log, tag + "gauge ")
        out.update({"lag": pd.concat([tab_r.assign(input="radar"), tab_g.assign(input="gauge")], ignore_index=True), "d_radar": d_r, "d_gauge": d_g})
        out["minutes"] = (time.time() - t0) / 60
        return out
    res = ladder.two_stage(Dk, le, bi, S["curves"], S["gauge"], inc_F=start_F, inc_P=start_P, hold_ext=False, log=log, tag=tag)
    dens = state_density(S, Dk, ex)
    out["cv"] = pd.concat([res["wet"]["table"], res["dry"]["table"]], ignore_index=True)
    out["tau_wet"] = res.get("tau_wet")
    out["after_rain"] = after_rain(S, Dk, res)
    out["judge"], out["cells"] = judge(S, spec, res, dens, gate_ok)
    out["extension"] = extension_table(S, res["model"], res["x"], start_F, start_P)
    if "x_alt" in res:
        out["extension_alt"] = extension_table(S, res["model"], res["x_alt"], start_F, start_P)
    if "ghost" in spec:                                         # exclude the flagged events and fit again
        fl = tail_flags(S, Dk)
        fl["changed_in_truth"] = fl.event_id.isin(ex["ghost_events"])
        out["flags"] = fl
        use = (~le.event_id.isin(set(fl.event_id[fl.flagged]))).to_numpy().astype(float)
        res2 = ladder.two_stage(Dk, le, bi, S["curves"], S["gauge"], inc_F=start_F, inc_P=start_P, hold_ext=False, use=use, log=log, tag=tag + "excl. ")
        out["judge_excluded"], out["cells_excluded"] = judge(S, spec, res2, dens, gate_ok)
        out["cv_excluded"] = pd.concat([res2["wet"]["table"], res2["dry"]["table"]], ignore_index=True)
    out["minutes"] = (time.time() - t0) / 60
    return out


def verdict(name, r):
    """The verdict of one truth by the rules of decisions.md ("A6 final: the truths, amended")."""
    if name in ("a_lag", "e"):
        return "report"
    j = r["judge"]
    ok = bool(j.rung_ok.all() and j.k_ok.all())
    if name in GATES:
        return "PASS" if ok else "FAIL"
    if name == "b":
        return "k_dry recovered" if ok else "k_dry NOT recovered"
    if name == "d_high":
        w, d = (j[j.regime == r_].iloc[0].selected_rung for r_ in ("wetting", "drying"))
        return f"dependence above 3.2 mm/h: wetting {'detected' if w != 'K0' else 'not detected'}, drying {'detected' if d != 'K0' else 'not detected'}"
    if name == "f":
        d1 = j[j.regime == "drying"].iloc[0]
        d2 = r["judge_excluded"][r["judge_excluded"].regime == "drying"].iloc[0]
        spur, rest = not (d1.rung_ok and d1.k_ok), bool(d2.rung_ok and d2.k_ok)
        return ("spurious k_dry" if spur else "no spurious k_dry") + (", restored by the exclusion" if rest else ", NOT restored by the exclusion")
    if name == "g":
        return "k not corrupted" if ok else "k corrupted"
    return "detectability"


def detect_study(S, cache, R, only_report):
    """The detectability study: per form and size (from the largest) up to 5 replicates; detected in a regime = the true rung
    is returned in at least 4 of 5. Replicates stop when every regime is decided; a form stops at the first size at which
    no regime is detected. Returns the table of every replicate run."""
    Dt = A6["detect"]
    rows = []
    for form, regimes in FORMS.items():
        for a in Dt["sizes"]:
            hits = {r: 0 for r in regimes}
            n = 0
            for rep in range(Dt["replicates"]):
                name = det_name(form, a, rep)
                f = cache / f"{name}.pkl"
                if f.exists():
                    R[name] = pickle.loads(f.read_bytes())
                elif only_report:
                    break
                else:
                    log(f"=== detectability {name} ===")
                    R[name] = run_truth(name, S)
                    f.write_bytes(pickle.dumps(R[name]))
                n += 1
                j = R[name]["judge"].set_index("regime")
                for r in regimes:
                    hits[r] += int(j.rung_ok[r])
                    rows.append({"form": form, "size_a": a, "replicate": rep, "regime": r, "true_rung": j.true_rung[r], "rung_returned": j.selected_rung[r],
                                 "hit": bool(j.rung_ok[r]), "gain_db_of_true_rung": gain_of(R[name]["cv"], r, j.true_rung[r]), "share_within_25pct": j.share_within_25pct[r]})
                need = Dt["replicates"] - 1
                if all(hits[r] >= need or (n - hits[r]) >= 2 for r in regimes):
                    break                                       # decided for every regime
            if not any(hits[r] >= Dt["replicates"] - 1 for r in regimes):
                break                                           # nothing detected at this size: smaller sizes are not tried
    return pd.DataFrame(rows)


def detect_summary(t):
    """Per form, regime and size: hits / replicates run, detected (at least 4 hits), the ratio of k at the domain ends."""
    if not len(t):
        return t
    Dt = A6["detect"]
    g = t.groupby(["form", "regime", "size_a"]).agg(replicates_run=("hit", "size"), hits=("hit", "sum"), median_gain_db=("gain_db_of_true_rung", "median")).reset_index()
    g["detected"] = g.hits >= Dt["replicates"] - 1
    g["k_ratio_at_domain_ends"] = [round(1 + 2 * a, 3) if f == "cw" else round((1 + Dt["r_cap_mmh"] / 2 * a) ** (1 if r == "wetting" else -1), 3)
                                   for f, r, a in zip(g.form, g.regime, g.size_a)]
    g["form"] = g.form.map({"cw": "k_dry(W)", "kr": "k(R) inside the measured range"})
    return g.sort_values(["form", "regime", "size_a"], ascending=[True, True, False]).reset_index(drop=True)


# =================================================================================================== report
def k_text(r):
    if r.true_rung == "none":
        return f"1/k {r.tau_fit_min:.1f} min at R = 0, W = 0"
    alt = f" / {r.tau_fit_alt_min:.1f}" if np.isfinite(r.tau_fit_alt_min) else ""
    if r.true_rung == "K0":
        return f"1/k {r.tau_fit_min:.1f}{alt} min (truth {r.tau_true_min:.1f})" + (f": {r.k_note}" if r.k_note else (": yes" if r.k_ok else ": no"))
    sh = f"{100 * r.share_within_25pct:.0f} %" + (f" / {100 * r.share_alt:.0f} %" if np.isfinite(r.share_alt) else "")
    whole = f" (whole domain {100 * r.share_whole_domain:.0f} %)" if abs(r.share_whole_domain - r.share_within_25pct) > 1e-9 else ""
    return f"{sh} of the domain weight{whole}" + (": yes" if r.k_ok else ": no")


def fig_surfaces(R):
    names = [n for n in ("c", "d") if n in R]
    if not names:
        return None
    regimes = [r for r in ("wetting", "drying") if r in set(R[names[0]]["judge"].regime)]      # the regimes that are judged (both in A6)
    fig, axes = plt.subplots(len(names), len(regimes), figsize=(5.5 * len(regimes), 4 * len(names)), squeeze=False)
    for i, n in enumerate(names):
        c = R[n]["cells"]
        for j, regime in enumerate(regimes):
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
    g = cv[(cv["select"] == regime) & (cv.rung == rung) & cv.best_lambda]
    return round(float(g.gain_db.iloc[0]), 5) if len(g) and g.gain_db.iloc[0] is not None else np.nan


A6_INTRO = ("Written by `research/waa_ksurface/a6_recovery.py` (libraries `fake6.py`, `ladder.py`, `ksurf.py`). Fake attenuation on the 4 563 TRAIN link-events of the learning set "
         "(real links, real rain), 0.3 dB quantization, the measured noise, the T1 rain error at its high end for the radar input and the T2 gauge model for the gauge input. "
         "**Every truth frees the F / P extension, as the real fit will (Deviation 7), and is drawn with a fresh seed.** One replicate per truth (five per size in the detectability study). "
         "Selection: a richer rung must lower the CV error by more than 0.005 dB with a paired 90 % interval that excludes zero (Deviation 5). k_wet is a bounded nuisance "
         "(Deviation 6): where two numbers are given as `x / y`, the first is with k_wet at its gauge-input estimate, the second at its radar-input estimate. "
         "For the truths in which k depends on R the share is taken inside F's measured range (restriction R2). No real attenuation is used as a target. VAL and TEST were not opened.\n")


def report(R, det, P=P, title=None, intro=None, stop_on=None, verdict=None):
    """Tables, figure and report of the truths in R. P = where they are written (default: the A6 files); title, intro =
    the head of the report; stop_on(name, verdict) = an extra STOP rule; verdict = another verdict function."""
    verdict = verdict or globals()["verdict"]
    rows, stop = [], None
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
        if ((n in GATES and v == "FAIL") or (stop_on is not None and stop_on(n, v))) and stop is None:
            stop = f"truth ({n})"
    table = pd.DataFrame(rows)
    save(table, P.A6_VERDICTS, index=False)
    names = [n for n in R if "cv" in R[n]]
    if names:
        save(pd.concat([R[n]["cv"].assign(truth=n) for n in names] + ([R["f"]["cv_excluded"].assign(truth="f_excluded")] if "f" in R else []), ignore_index=True), P.A6_CV, index=False)
        save(pd.concat([R[n]["cells"].assign(truth=n) for n in names if n in ORDER], ignore_index=True).round(5), P.A6_CELLS, index=False)
        ext = pd.concat([R[n][k].assign(truth=n, k_wet_at=w) for n in names if n in ORDER for k, w in (("extension", "gauge estimate"), ("extension_alt", "radar estimate")) if k in R[n]], ignore_index=True)
        save(ext.round(5), P.A6_EXTENSION, index=False)
    fig = fig_surfaces(R)
    if fig is not None:
        save(fig, P.A6_FIG, dpi=120)
        plt.close(fig)
    a = [title or "# A6, final form: fake-data recovery with the procedure of the real fit\n", intro or A6_INTRO,
         "## The table\n", md_table(table) if len(table) else "no truth was run",
         f"\n(`{rel(P.A6_VERDICTS)}`; every rung and lambda: `{rel(P.A6_CV)}`; every domain cell: `{rel(P.A6_CELLS)}`; every extension knot: `{rel(P.A6_EXTENSION)}`"
         + (f"; figure `{rel(P.A6_FIG)}`" if fig is not None else "") + ")\n"]
    if stop:
        a.append(f"**STOP: {stop} FAILS its rule. What comes after it was not run. Nothing was patched.**\n")
    if det is not None and len(det):
        save(det, P.A6_DETECT_RUNS, index=False)
        summ = detect_summary(det)
        save(summ, P.A6_DETECT, index=False)
        a += ["## Detectability: the smallest dependence the ladder returns at the margin of its selection\n",
              "Detected = the true rung is returned in at least 4 of 5 replicates. Sizes are tried from the largest; the replicates of a size stop once the outcome is decided, "
              "and a form stops at the first size at which nothing is detected. `k_ratio_at_domain_ends` = k(W = 2 dB) / k(W = 0) for k_dry(W); k(3.2 mm/h) / k(0) for k(R).\n",
              md_table(summ.round(5)), f"\n(`{rel(P.A6_DETECT)}`; every replicate: `{rel(P.A6_DETECT_RUNS)}`)\n"]
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
        if n in ("b", "d_high"):
            e = R[n]["extension"].assign(k_wet_at="gauge estimate")
            if "extension_alt" in R[n]:
                e = pd.concat([e, R[n]["extension_alt"].assign(k_wet_at="radar estimate")], ignore_index=True)
            a.append("The extension knots (fitted with k; started from F held and P at half its true increments):\n")
            a.append(md_table(e[["band", "quantity", "knot_mmh", "true", "start", "fitted", "fitted_over_true", "within_25pct", "at_bound", "k_wet_at"]].round(3)))
            a.append("")
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
    a.append(md_table(pd.DataFrame([{"truth": n, "minutes": round(R[n]["minutes"], 1)} for n in R])))
    save("\n".join(a) + "\n", P.A6_REPORT)
    log(f"report: {rel(P.A6_REPORT)}")
    return stop


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    force, only_report = "--force" in sys.argv, "--report" in sys.argv
    cache = P.CACHE / "a6_recovery_final"
    cache.mkdir(parents=True, exist_ok=True)
    R, S, det, failed = {}, None, None, False
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
            log(f"=== truth ({n}): {verdict(n, R[n])} ({R[n]['minutes']:.1f} min) ===")
        if n in R and n in GATES and verdict(n, R[n]) == "FAIL":
            failed = True
            break                                               # the rule: a FAIL on a gate is a STOP
    if not failed and all(n in R for n in ORDER) and not args:
        S = S or (None if only_report else setting())
        det = detect_study(S, cache, R, only_report)
        R = {n: R[n] for n in R}
    stop = report({n: R[n] for n in ORDER if n in R}, det)
    if stop:
        sys.exit(1)


if __name__ == "__main__":
    main()
