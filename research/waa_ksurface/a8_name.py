"""Step 8: naming the drying law k_dry(W) (decisions.md, "Step 8: naming the drying law", logged before this code was run).

    python -m research.waa_ksurface.a8_name                 every part in order (a finished part is read from the cache)
    python -m research.waa_ksurface.a8_name --report        write the report from the cache
    python -m research.waa_ksurface.a8_name --one refs      one part in this process (used by the main run)
    python -m research.waa_ksurface.a8_name --freeze        Block B: write the frozen equation and the hashes (no fit, no data)

The parts, in order (one child process each):
    points          the A4 drying points (W bins, density weights) and ln k_dry of the A7 K1 fit at both k_wet estimates
    pysr_<e>_<s>    one PySR search (estimate e, seed s) on those points; its equations
    rank            the physics forms and the PySR Pareto fronts refitted at both estimates; AICc; the chosen form, the trial order
    refs            the free K1 and K0 of the review (sel_mid) at both estimates: held-out sums; code check against the review table
    plug_<j>        trial j: the form put back into the A7 model, refitted on the learning events, its 5-fold CV on L-dryhigh at
                    both estimates; accepted if within 0.00291 dB of the free K1 and better than K0 by more than 0.00291 dB at both.
                    The trials stop at the first accepted form.
    boot_<e>_<b>    block b of the 200 A7 bootstrap draws at estimate e: the K1 refit (the band; code check against the stored
                    A7 draws) and the accepted form refitted (its parameter intervals)

Real TRAIN data of the 14 learning events (the A7 setting, the A7 flag: 0 removed). VAL and TEST are not opened.
"""
import hashlib
import json
import pickle
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import psutil
import sympy as sp

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a6_recovery as a6, a7_real_fit as a7, a7_review as a7r, a8_forms as F8, ksurf, ladder
from research.waa_ksurface import s1_separability as s1, s3b_recovery as s3b, s3c_a7_variant as s3c, s3d_kwet_fixed as s3d
from research.waa_ksurface import paths as P

C8 = F8.C8
GB = s3b.GB
log = a6.log
CACHE = P.CACHE / "a8_name"
A7_CACHE = P.CACHE / "a7_real_fit"
R_CACHE = P.CACHE / "a7_review"
EST = ("gauge", "radar")
W_REPORT = [0.0, 0.5, 1.0, 2.0]                                 # the W of the stored A7 draws (the band's code check)
PYSR_LOSS = "myloss(x, y, w) = x > 0 ? w * (log(x) - log(y))^2 : w * (100 + abs(x))"


def load(name):
    return pickle.loads((CACHE / f"{name}.pkl").read_bytes())


def a7fit():
    return pickle.loads((A7_CACHE / "fit.pkl").read_bytes())


# =================================================================================================== points
def points_part():
    fit = a7fit()
    Wp, w, minutes = F8.points()
    return {"W": Wp, "w": w, "minutes": minutes, "y": {e: F8.k1_lnk(fit["x"][e], Wp) for e in EST}}


# =================================================================================================== PySR
def pysr_part(est, seed):
    from pysr import PySRRegressor                              # Julia starts in this child only
    pts = load("points")
    Y = C8["pysr"]
    m = PySRRegressor(niterations=Y["niterations"], populations=Y["populations"], maxsize=Y["maxsize"], binary_operators=["+", "-", "*", "/", "^"],
                      unary_operators=["exp"], constraints={"^": (-1, 1)}, nested_constraints={"exp": {"exp": 0}}, elementwise_loss=PYSR_LOSS,
                      precision=64, random_state=seed, deterministic=True, parallelism="serial", progress=False, verbosity=0,
                      timeout_in_seconds=Y["timeout_s"], output_directory=str(CACHE / "pysr_runs"), run_id=f"{est}_{seed}")
    t0 = time.time()
    m.fit(pts["W"][:, None], np.exp(pts["y"][est]), weights=pts["w"], variable_names=["W"])
    eq = m.equations_
    rows = [{"estimate": est, "seed": seed, "complexity": int(r.complexity), "loss": float(r.loss), "equation": str(r.equation),
             "sympy": str(r.sympy_format)} for r in eq.itertuples()]
    log(f"   (pysr {est} {seed}) {len(rows)} equations in {(time.time() - t0) / 60:.1f} min")
    return pd.DataFrame(rows)


def pysr_names():
    return [f"pysr_{e}_{s}" for e in EST for s in C8["pysr"]["seeds"]]


def pysr_front(tab):
    """Per estimate: per complexity the lowest loss over the seeds; kept if below every lower complexity's loss."""
    keep = []
    for e, g in tab.groupby("estimate"):
        best = g.loc[g.groupby("complexity").loss.idxmin()].sort_values("complexity")
        run = np.inf
        for i, r in best.iterrows():
            if r.loss < run:
                keep.append(i)
                run = r.loss
    t = tab.copy()
    t["front"] = t.index.isin(keep)
    return t


def structure(expr_text):
    ex = sp.sympify(expr_text, locals={"W": F8.W})
    return str(ex.xreplace({a: sp.Symbol("c") for a in ex.atoms(sp.Float)}))


# =================================================================================================== rank
def make_form(spec):
    return F8.physics()[spec["key"]] if spec["source"] == "physics" else F8.from_sympy(spec["key"], spec["expr"])


def rank_part():
    pts = load("points")
    tab = pd.concat([load(n) for n in pysr_names()], ignore_index=True)
    tab = pysr_front(tab)
    specs = [{"key": k, "source": "physics"} for k in F8.physics()]
    seen = {}
    for i, r in tab[tab.front].iterrows():
        specs.append({"key": f"pysr_{r.estimate}_c{r.complexity}", "source": "PySR", "expr": r.sympy, "struct": structure(r.sympy),
                      "pysr_complexity": int(r.complexity), "pysr_estimate": r.estimate})
    rows, fits = [], {}
    for s in specs:
        f = make_form(s)
        row = {"name": s["key"], "source": s["source"], "p": f.p, "form": f.text, "struct": s.get("struct", ""), "pysr_complexity": s.get("pysr_complexity")}
        th = {}
        for e in EST:
            th[e], wrss = F8.fit_curve(f, pts["W"], pts["y"][e], pts["w"])
            row[f"wrss_{e}"] = wrss
            row[f"aicc_{e}"] = F8.aicc(wrss, len(pts["W"]), f.p)
            row[f"formula_{e}"] = f.formula(th[e])
        rows.append(row)
        fits[s["key"]] = {"spec": s, "theta": th}
    t = pd.DataFrame(rows)
    # the same PySR expression found twice (two estimates or complexities): keep the one with the lower summed AICc
    t["aicc_sum"] = t.aicc_gauge + t.aicc_radar
    dup = t[t.source == "PySR"].sort_values("aicc_sum").duplicated("struct")
    t = t.drop(index=dup[dup].index).drop(columns="aicc_sum").reset_index(drop=True)
    t, order = F8.rank(t)
    return {"pysr": tab, "table": t, "order": [t.name[i] for i in order], "fits": fits}


# =================================================================================================== the real setting
def setting():
    S, curves, band_idx, use = a7r.setting("mid")
    D = S["D"]
    fit = a7fit()
    prob = ksurf.prepare(D, band_idx, curves, D["rain"], 0)
    prob["w5"] = np.asarray(fit["stage2_w5"], float)
    used = (use > 0)[S["bin_le"]]
    _, _, high, _ = s1.score_sets(S, D)
    return S, prob, use, high & used, ladder.folds(S["le"], use), fit


def refs_part():
    S, prob, use, score, fold, fit = setting()
    sel = pickle.loads((R_CACHE / "sel_mid.pkl").read_bytes())
    assert sel.get("reproduces")
    tab = pd.read_csv(P.A7R_CV)
    tab = tab[tab.f_curve == "mid"]
    out, worst = {}, 0.0
    for e in EST:
        for rung in ("K0", "K1"):
            lam, x = sel["est"][e]["fits"][rung]
            m = ksurf.Model(prob, "K0", rung, lam)
            m.hold("gw", "hw")
            S_, N_ = F8.held_out(m, x, use, fold, score, prob)
            ref = tab[(tab.k_wet_estimate == e) & (tab.rung == rung) & np.isclose(tab["lambda"], lam)].cv_loss.item()
            dev = abs((S_.sum() / N_.sum()) / ref - 1)
            worst = max(worst, dev)
            out[(e, rung)] = {"lam": lam, "x": x, "S": S_, "N": N_, "cv_db": ladder.cv_db(S_, N_), "dev": dev}
            log(f"   (refs) {e} {rung} (lambda {lam:g}): CV {ladder.cv_db(S_, N_):.6f} dB; relative difference to the review table {dev:.1e}")
    ok = worst <= C8["reproduce_rtol"]
    log(f"=== (refs) code check against a7_review_cv.csv: {'REPRODUCED' if ok else 'NOT REPRODUCED'} (largest {worst:.1e}) ===")
    return {"refs": out, "ok": ok, "worst": worst, "wet": {e: (sel["est"][e]["fits"]["K1"][1][ksurf.Model(prob, 'K0', 'K1').sl["gw"]].copy(),
                                                             np.zeros(0)) for e in EST}}


# =================================================================================================== plug-back
def plug_part(j):
    rk, rf = load("rank"), load("refs")
    assert rf["ok"]
    name = rk["order"][j]
    spec, theta = rk["fits"][name]["spec"], rk["fits"][name]["theta"]
    f = make_form(spec)
    S, prob, use, score, fold, fit = setting()
    margin = s3c.margin()
    Sel = ladder.C["select"]
    res = {"name": name, "j": j, "spec": spec, "est": {}}
    for e in EST:
        t0 = time.time()
        fm = F8.FormModel(prob, f, rf["wet"][e])
        z = fm.fit(fm.start(theta[e], fit["x"][e]), use)
        info = dict(fm.last)
        S_, N_ = F8.held_out(fm, z, use, fold, score, prob)
        k1, k0 = rf["refs"][(e, "K1")], rf["refs"][(e, "K0")]
        assert np.array_equal(N_, k1["N"]) and np.array_equal(N_, k0["N"])
        cv = ladder.cv_db(S_, N_)
        d_k1, g_k0 = cv - k1["cv_db"], k0["cv_db"] - cv
        i_k1 = ladder.paired_interval(S_, k1["S"], N_, Sel["boot_reps"], Sel["seed"], Sel["interval"])
        i_k0 = ladder.paired_interval(k0["S"], S_, N_, Sel["boot_reps"], Sel["seed"], Sel["interval"])
        ok = bool(d_k1 <= margin and g_k0 > margin)
        res["est"][e] = {"z": z, "theta": z[:f.p], "wet": rf["wet"][e], "train_loss": info["loss"], "iterations": info["nit"], "converged": info["converged"],
                         "cv_db": cv, "cv_K1_db": k1["cv_db"], "cv_K0_db": k0["cv_db"], "form_minus_K1": d_k1, "form_minus_K1_90": tuple(i_k1),
                         "K0_minus_form": g_k0, "K0_minus_form_90": tuple(i_k0), "pass": ok, "formula": f.formula(z[:f.p]),
                         "tau_dry": dict(zip(W_REPORT, 1.0 / f.k(np.array(W_REPORT), z[:f.p]))), "minutes": (time.time() - t0) / 60}
        log(f"   (plug {j} {name} {e}) CV {cv:.6f} dB; form - K1 {d_k1:+.5f} ({i_k1[0]:+.5f} to {i_k1[1]:+.5f}); K0 - form {g_k0:+.5f} "
            f"({i_k0[0]:+.5f} to {i_k0[1]:+.5f}); margin {margin:.5f}: {'pass' if ok else 'fail'} ({res['est'][e]['minutes']:.1f} min)")
    res["accepted"] = all(res["est"][e]["pass"] for e in EST)
    return res


def plug_names():
    rk = load("rank")
    return [f"plug_{j}" for j in range(len(rk["order"]))]


def accepted():
    """The accepted plug-back trial (or None), from the cache."""
    for n in plug_names():
        f = CACHE / f"{n}.pkl"
        if f.exists() and load(n)["accepted"]:
            return load(n)
    return None


# =================================================================================================== bootstrap
def draws():
    S, _, _, use = a7r.setting("mid")
    ev_all = S["le"].event_id.to_numpy()
    events = np.sort(np.unique(ev_all[use > 0]))
    rng = np.random.default_rng(a7.A7["boot_seed"])            # the A7 sequence
    return [rng.choice(events, len(events), replace=True) for _ in range(C8["boot_draws"])]


def boot_part(est, block):
    S, prob, use, score, fold, fit = setting()
    stored = pickle.loads((A7_CACHE / f"boot_{est}.pkl").read_bytes())
    m = ksurf.Model(prob, "K0", fit["rung"], fit["lam"])
    m.hold("gw", "hw")
    acc = accepted()
    fm = F8.FormModel(prob, make_form(acc["spec"]), acc["est"][est]["wet"]) if acc else None
    ev_all = S["le"].event_id.to_numpy()
    picks = draws()
    lo, hi = block * C8["boot_block"], min((block + 1) * C8["boot_block"], C8["boot_draws"])
    rows, worst = [], 0.0
    t0 = time.time()
    for b in range(lo, hi):
        cnt = pd.Series(picks[b]).value_counts()
        w = use * pd.Series(ev_all).map(cnt).fillna(0).to_numpy()
        xb = m.fit(fit["x"][est].copy(), w)
        vals = s3d.k_dry_at(m, xb, W_REPORT)
        ref = stored.loc[stored.draw == b, [f"1/k_dry at W = {W:g} dB" for W in W_REPORT]].to_numpy()[0]
        dev = float(np.max(np.abs(vals / ref - 1)))
        worst = max(worst, dev)
        r = {"draw": b, "x_K1": xb, "K1_converged": m.last["converged"], "dev": dev}
        if fm is not None:
            zb = fm.fit(acc["est"][est]["z"].copy(), w)
            r.update({"theta": zb[:fm.form.p], "form_converged": fm.last["converged"]})
        rows.append(r)
        log(f"   (boot {est} draw {b}) K1 check {dev:.1e}" + (f"; form {'converged' if r['form_converged'] else 'NOT converged'}" if fm else "")
            + f" ({(time.time() - t0) / 60:.1f} min in this block)")
    return {"est": est, "block": block, "rows": rows, "worst": worst, "ok": worst <= C8["reproduce_rtol"], "form": acc["name"] if acc else None}


def boot_names():
    nb = -(-C8["boot_draws"] // C8["boot_block"])
    return [f"boot_{e}_{b}" for e in EST for b in range(nb)]


# =================================================================================================== report
def band(boots, est, Wv):
    """5 / 50 / 95 % of 1 / k_dry of the K1 refits of the draws, at Wv."""
    m = F8.k1_model()
    tau = np.array([1.0 / m.k(r["x_K1"], "drying", 0.0, Wv) for b in boots if b["est"] == est for r in b["rows"]])
    return np.percentile(tau, [5, 50, 95], axis=0), len(tau)


def figure(pts, rk, acc, boots):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    Wf = np.geomspace(0.1, 9.0, 300)
    styles = {"linear": ":", "sat": "--", "power": "-.", "sat_floor": (0, (5, 1, 1, 1, 1, 1))}
    fit = a7fit()
    for ax, e in zip(axes, EST):
        if boots:
            (lo, _, hi), n = band(boots, e, Wf)
            ax.fill_between(Wf, lo, hi, facecolor="none", edgecolor="0.55", hatch="////", lw=0, label=f"K1 band 5-95 % ({n} draws, hatched)")
        ax.plot(Wf, 1.0 / np.exp(F8.k1_lnk(fit["x"][e], Wf)), "-", color="k", lw=2.4, label="free K1 (A7, thick solid)")
        ax.scatter(pts["W"], 1.0 / np.exp(pts["y"][e]), s=8 + 6 * pts["w"], color="k", zorder=5, label="A4 points (size = weight)")
        for key, ls in styles.items():
            f = F8.physics()[key]
            th = rk["fits"][key]["theta"][e]
            ax.plot(Wf, 1.0 / f.k(Wf, th), linestyle=ls, color="0.3", lw=1.3, label=f"{f.text} ({key})")
        if acc is not None:
            f = make_form(acc["spec"])
            ax.plot(Wf, 1.0 / f.k(Wf, acc["est"][e]["theta"]), "-", color="0.15", lw=1.0, marker="o", markevery=25, ms=4,
                    label=f"accepted, plugged back ({acc['name']}, circles)")
        ax.set(xscale="log", yscale="log", xlabel="W (dB)", title=f"k_wet at its {e}-input estimate")
        ax.grid(alpha=0.25, which="both")
        for sp_ in ("top", "right"):
            ax.spines[sp_].set_visible(False)
    axes[0].set_ylabel("1 / k_dry (min)")
    axes[1].legend(frameon=False, fontsize=7, loc="upper right")
    fig.suptitle("Step 8: the A7 drying curve and its closed-form candidates (fits on ln k, A4 density weights)", fontsize=10)
    fig.tight_layout()
    return fig


def report(R, state):
    pts, rk, rf = R.get("points"), R.get("rank"), R.get("refs")
    m = s3c.margin()
    L = ["# Step 8: naming the drying law\n",
         "Written by `research/waa_ksurface/a8_name.py` (decisions.md, \"Step 8: naming the drying law\", logged before the run). Real TRAIN data of "
         "the 14 learning events (the A7 setting; A7 flag: 0 removed). VAL and TEST were not opened.\n"]
    if state and state.get("incomplete"):
        L += [f"**Run incomplete: {state['incomplete']}**\n"]
    if pts is None:
        save("\n".join(L) + "\n", P.A8_REPORT)
        return
    pt = pd.DataFrame({"W (dB)": pts["W"], "A4 minutes": pts["minutes"].astype(int), "weight": pts["w"],
                       **{f"1/k_dry K1 {e} (min)": 1.0 / np.exp(pts["y"][e]) for e in EST}})
    boots = [R[n] for n in boot_names() if n in R]
    band_ok = bool(boots) and all(b["ok"] for b in boots)
    if band_ok:
        for e in EST:
            (lo, _, hi), _ = band(boots, e, pts["W"])
            pt[f"band {e} 5-95 % (min)"] = [f"{a:.1f} to {b:.1f}" for a, b in zip(lo, hi)]
    save(pt, P.A8_POINTS, index=False)
    L += ["## 1. The input\n", f"{len(pts['W'])} points (the A4 drying W bins, pooled), weight = A4 minutes scaled to mean 1; y = ln k_dry of the A7 K1 fit. "
          f"The 0-0.3 dB bin carries {pts['w'][0]:.1f} of {len(pts['W'])}.\n"]
    if boots:
        L += [f"Band: {sum(len(b['rows']) for b in boots)} bootstrap K1 refits; code check against the stored A7 draws: "
              f"{'REPRODUCED' if band_ok else 'NOT REPRODUCED'} (largest relative difference {max(b['worst'] for b in boots):.1e}).\n"]
    L += [md_table(pt.round(4)), ""]
    if rk is None:
        save("\n".join(L) + "\n", P.A8_REPORT)
        return
    t = rk["table"]
    save(rk["pysr"], P.A8_PYSR, index=False)
    save(t, P.A8_CAND, index=False)
    show = t[["name", "source", "p", "aicc_gauge", "daicc_gauge", "aicc_radar", "daicc_radar", "dmax", "front", "best", "chosen", "trial", "formula_radar"]]
    L += ["## 2. Candidates (AICc on the same weighted points)\n",
          f"PySR: {len(rk['pysr'])} equations from {len(pysr_names())} searches; {int(rk['pysr'].front.sum())} on the per-estimate fronts; "
          f"{int((t.source == 'PySR').sum())} distinct expressions kept as candidates. Rule: the fewest parameters within {C8['aicc_within']:g} of the best "
          f"Delta_max wins; trial order = the chosen form, then the (p, Delta_max) front by Delta_max (Delta_max <= {C8['dmax_try']:g}, at most "
          f"{C8['max_trials']}).\n", md_table(show.round(3)), "",
          f"Chosen form: **{rk['order'][0]}**. Trial order: {', '.join(rk['order'])}.\n"]
    if rf is not None:
        L += [f"References (review sel_mid, held-out sums recomputed): code check {'REPRODUCED' if rf['ok'] else 'NOT REPRODUCED'} "
              f"(largest relative difference {rf['worst']:.1e}). CV (dB): " + "; ".join(
                  f"{e}: K1 {rf['refs'][(e, 'K1')]['cv_db']:.5f}, K0 {rf['refs'][(e, 'K0')]['cv_db']:.5f}" for e in EST) + ".\n"]
    plugs = [R[n] for n in plug_names() if n in R]
    acc = next((p for p in plugs if p["accepted"]), None)
    if plugs:
        rows = []
        for p in plugs:
            for e in EST:
                q = p["est"][e]
                rows.append({"trial": p["j"], "form": p["name"], "k_wet estimate": e, "CV form (dB)": q["cv_db"], "CV free K1 (dB)": q["cv_K1_db"],
                             "CV K0 (dB)": q["cv_K0_db"], "form - K1 (dB)": q["form_minus_K1"], "90 % (form - K1)": "{:+.5f} to {:+.5f}".format(*q["form_minus_K1_90"]),
                             "K0 - form (dB)": q["K0_minus_form"], "90 % (K0 - form)": "{:+.5f} to {:+.5f}".format(*q["K0_minus_form_90"]),
                             "pass": q["pass"], "converged": q["converged"], "formula": q["formula"]})
        pb = pd.DataFrame(rows)
        save(pb, P.A8_PLUG, index=False)
        L += ["## 3. Plug-back\n", f"Accept: at BOTH estimates, CV(form) - CV(free K1) <= {m:.5f} dB and CV(K0) - CV(form) > {m:.5f} dB (intervals reported, "
              "not part of the rule).\n", md_table(pb.drop(columns="formula").round(5)), "",
              f"**Accepted: {acc['name'] if acc else 'none'}**" + ("" if acc else " (the free K1 curve remains the stated law; no closed form claimed)") + ".\n"]
    if acc is not None and band_ok:
        f = make_form(acc["spec"])
        prow = []
        for e in EST:
            th = np.array([r["theta"] for b in boots if b["est"] == e for r in b["rows"]])
            conv = sum(r["form_converged"] for b in boots if b["est"] == e for r in b["rows"])
            nat = pd.DataFrame([f.natural(x) for x in th])
            full = f.natural(acc["est"][e]["theta"])
            for k, v in full.items():
                prow.append({"k_wet estimate": e, "parameter": k, "value": v, "5 %": np.percentile(nat[k], 5), "95 %": np.percentile(nat[k], 95),
                             "draws": len(th), "converged": conv})
        pr = pd.DataFrame(prow)
        save(pr, P.A8_PARAMS, index=False)
        L += ["## 4. The accepted formula, parameters with 90 % event-bootstrap intervals\n", f"`{f.text}`\n", md_table(pr.round(5)), ""]
    if pts is not None and rk is not None:
        fig = figure(pts, rk, acc, boots if band_ok else [])
        save(fig, P.A8_FIG)
        plt.close(fig)
        L += [f"![curves](../figures/{P.A8_FIG.name})\n"]
    save("\n".join(L) + "\n", P.A8_REPORT)
    log(f"report: {rel(P.A8_REPORT)}")


# =================================================================================================== Block B freeze (no data, no fit)
def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def frozen_equation():
    """The equation frozen for Block B: k_wet at the RADAR estimate; drying = the accepted form (its plug-back refit at the radar
    estimate) or, if none was accepted, the free K1 of A7 at the radar estimate (x_alt); F, P measured (E4) plus the extensions
    of that fit (EFFECTIVE, R1). Every number needed to run it, in the ksurf parameter layout."""
    fit = a7fit()
    acc = accepted()
    cfg = ksurf.config()
    bands = list(a7.fp.measured())
    if acc is not None:
        f = make_form(acc["spec"])
        q = acc["est"]["radar"]
        m = ksurf.Model({"nb": len(bands)}, "K0", "K1", 0.0, cfg={**cfg, "w_knots_dry_db": [float(v) for v in F8.GRID]})
        p = f.p
        z = np.asarray(q["z"], float)
        L = np.clip(f.lnk(F8.GRID_EVAL, z[:p]), *m.lnk)
        x = np.zeros(m.n)
        x[m.sl["gw"]] = q["wet"][0]
        x[m.sl["gd"]], x[m.sl["hd"]] = L[0], L[1:] - L[0]
        nF = m.sl["F"].stop - m.sl["F"].start
        x[m.sl["F"]], x[m.sl["P"]] = z[p:p + nF], z[p + nF:]
        drying = {"kind": "closed form (step 8, accepted)", "name": acc["name"], "formula": f.text, "sympy": sp.sstr(f.expr),
                  "theta": [float(v) for v in z[:p]], "theta_symbols": [str(s) for s in f.params],
                  "natural": {k: float(v) for k, v in f.natural(z[:p]).items()},
                  "kernel": "ln k tabulated on w_knots_dry_db (value at W = 0 = the form at 1e-6 dB), clipped to the ln k bounds, read as K1 knots"}
        knots = [float(v) for v in F8.GRID]
    else:
        x = np.asarray(fit["x"]["radar"], float)
        m = ksurf.Model({"nb": len(bands)}, "K0", "K1", 0.0)
        drying = {"kind": "free K1 of A7 (no closed form accepted in step 8)", "lambda": fit["lam"]}
        knots = [float(v) for v in cfg["w_knots_dry_db"]]
    ext = {}
    for b, band_ in enumerate(bands):
        iF = x[m.sl["F"]].reshape(len(bands), -1)[b]
        iP = x[m.sl["P"]].reshape(len(bands), -1)[b]
        c = a7.fp.measured()[band_]
        ext[band_] = {"F_knots_mmh": cfg["f_knots_mmh"], "F_at_knots_db": [float(v) for v in ksurf.full_curve(c, "F", iF, cfg["f_knots_mmh"], cfg)],
                      "P_knots_mmh": cfg["p_knots_mmh"], "P_at_knots_db_km": [float(v) for v in ksurf.full_curve(c, "P", iP, cfg["p_knots_mmh"], cfg)],
                      "F_measured_max_mmh": float(c.F_max), "P_measured_max_mmh": float(c.P_max)}
    return {"written_by": "research/waa_ksurface/a8_name.py --freeze", "date": time.strftime("%Y-%m-%d"),
            "equation": ["v = F_full(R(t)) (the steady WAA at the rain of this minute, per band)",
                         "wetting (v > W): k = k_wet (constant, K0); drying (v <= W): k = k_dry(W)",
                         "W <- v + (W - v) * exp(-k * 1 min), W = 0 at the start of the window",
                         "A - c = W + L * P_full(R(t)) (lag d = 0)",
                         "after rain (R = 0, F = 0): dW/dt = -k_dry(W) * W"],
            "rain_input": "R = radar path mean / g_e, 1 min (data layer W1); k_wet at its RADAR-input estimate",
            "bands": bands, "tau_wet_min": float(1.0 / np.exp(x[m.sl["gw"]][0])), "drying": drying,
            "ksurf_parameter_vector": [float(v) for v in x], "w_knots_dry_db": knots, "ksurf_config_a5": cfg,
            "F_P_measured": {"file": rel(P.E4_CELLS), "builder": "research/waa_ksurface/fp.py: measured()"}, "extensions_EFFECTIVE_R1": ext,
            "claims": "R1, R3 (revised), R4, R5, R6 and the E20150831_0455 caveat (decisions.md, A7 pre-registration and Deviation 15)"}


def freeze():
    eq = frozen_equation()
    save(json.dumps(eq, indent=1) + "\n", P.B_EQUATION)
    assert P.B_PREREG.exists(), "write docs/waa_ksurface/block_b_prereg.md first"
    files = {"frozen equation": P.B_EQUATION, "Block B rules": P.B_PREREG,
             "competitor: Schleiss (b), Pastorek (b) stored TRAIN fits (W3)": P.P1_LITERATURE_FITS,
             "competitor: the Part 1 equation, stored frozen fit (W3)": P.WP.S9_EQUATION,
             "competitor: the Part 1 freeze record (read only)": P.WP.S8_FREEZE,
             "F and P measured (E4 cells)": P.E4_CELLS,
             "model code: ksurf.py": P.MODULE / "ksurf.py", "model code: a8_forms.py": P.MODULE / "a8_forms.py", "model code: fp.py": P.MODULE / "fp.py"}
    out = {"date": time.strftime("%Y-%m-%d %H:%M:%S"), "sha256": {k: {"path": rel(v), "sha256": sha256(v)} for k, v in files.items()},
           "VAL": "not opened", "TEST": "not opened"}
    save(json.dumps(out, indent=1) + "\n", P.B_FREEZE)
    for k, v in out["sha256"].items():
        log(f"   {v['sha256'][:16]}  {v['path']}  ({k})")


# =================================================================================================== main
def run_part(name):
    if name == "points":
        return points_part()
    if name.startswith("pysr_"):
        _, e, s = name.split("_")
        return pysr_part(e, int(s))
    if name == "rank":
        return rank_part()
    if name == "refs":
        return refs_part()
    if name.startswith("plug_"):
        return plug_part(int(name.split("_")[1]))
    _, e, b = name.split("_")
    return boot_part(e, int(b))


def cached():
    names = ["points"] + pysr_names() + ["rank", "refs"]
    R = {n: load(n) for n in names if (CACHE / f"{n}.pkl").exists()}
    if "rank" in R:
        R.update({n: load(n) for n in plug_names() + boot_names() if (CACHE / f"{n}.pkl").exists()})
    return R


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    if "--one" in sys.argv:
        name = sys.argv[sys.argv.index("--one") + 1]
        t0 = time.time()
        r = run_part(name)
        (CACHE / f"{name}.pkl").write_bytes(pickle.dumps(r))
        log(f"=== (a8 {name}) done: {(time.time() - t0) / 60:.1f} min; peak memory {psutil.Process().memory_info().peak_wset / GB:.2f} GB ===")
        return
    if "--report" in sys.argv:
        report(cached(), None)
        return
    if "--freeze" in sys.argv:
        freeze()
        return
    state = s3d.start_state("step 8")
    mod = "research.waa_ksurface.a8_name"
    lab = lambda n: ""                                          # noqa: E731
    R, stop = s3d.run_children(mod, CACHE, ["points"] + pysr_names() + ["rank"], state, lab)
    if not stop and not state["incomplete"]:
        R2, stop = s3d.run_children(mod, CACHE, ["refs"], state, lab, lambda n, r: "" if r["ok"] else "the references did not reproduce the review table")
        R.update(R2)
    if not stop and not state["incomplete"]:
        for n in plug_names():
            R3, stop = s3d.run_children(mod, CACHE, [n], state, lab)
            R.update(R3)
            if stop or state["incomplete"] or R3[n]["accepted"]:
                break
    if not stop and not state["incomplete"]:
        R4, stop = s3d.run_children(mod, CACHE, boot_names(), state, lab, lambda n, r: "" if r["ok"] else f"({n}) the K1 refits did not reproduce the stored A7 draws")
        R.update(R4)
    if stop:
        state["incomplete"] = f"STOP: {stop}"
        log(f"=== STOP: {stop} ===")
    report(R, state)
    log(f"=== step 8 ended {time.strftime('%Y-%m-%d %H:%M:%S')} ===")


if __name__ == "__main__":
    main()
