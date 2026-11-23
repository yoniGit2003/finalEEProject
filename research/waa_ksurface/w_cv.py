"""WAA-only comparison (decisions.md, "WAA-only comparison", logged before any code was run): every model gets the SAME
fixed path law; only its WAA part is refitted, leave-one-event-out on the 14 TRAIN learning events. VAL and TEST closed.

    python -m research.waa_ksurface.w_cv            every fit (cached per path, model, fold), then metrics, figures, report
    python -m research.waa_ksurface.w_cv --smoke    fold 0 of every model with a few iterations, to a scratch cache (not a result)
    python -m research.waa_ksurface.w_cv --report   metrics, figures and report from the cached fits

Each model fits W to (A - c) - path with the mean Huber 0.3 dB of the 5 min values; its total = W + path. Path A = L P_full(R)
(measured P + the frozen extension), path B = L c_band R (Part 1's c_band). Lag 0 everywhere.
"""
import json
import pickle
import sys
import time

import numpy as np
import pandas as pd
import psutil
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numba import njit
from scipy.optimize import minimize, minimize_scalar

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a7_review as a7r, a8_forms as F8, b_val as B, data, fake6, ksurf, s1_separability as s1, s3_recovery as s3
from research.waa_ksurface import paths as P
from research.wet_antenna import literature_s8c as lit, system_id as sid

CW = data.branch_config()["w_cv"]
CACHE = P.CACHE / "w_cv"
GB = 1024 ** 3
PADM = sid.PADM
MODELS = ["ours", "Part 1 equation", "Schleiss", "Pastorek", "no memory", "no WAA"]
WITH_A = ("Part 1 equation", "Schleiss", "Pastorek")
PATHS = ("A", "B")
DELTA = 0.3
COLOR = {"ours": "#0072B2", "Part 1 equation": "#D55E00", "Schleiss": "#009E73", "Pastorek": "#CC79A7", "no memory": "#E69F00", "no WAA": "#56B4E9"}
GREY, SPAN = "#8c8c8c", "#e3e3e3"


def log(s):
    print(f"{time.strftime('%H:%M:%S')} {s}", flush=True)


def huber(r):
    a = np.abs(r)
    return np.where(a <= DELTA, 0.5 * r * r, DELTA * (a - 0.5 * DELTA))


# =================================================================================================== the setting
def setting():
    S = s3.setting()
    learn = S["learn"]
    D = sid.subset(S["D"], learn)
    le = S["le"][learn].reset_index(drop=True)
    curves, band_idx = B.band_setting(le)
    n = D["le_n"]
    idx = np.concatenate([o + PADM + np.arange(k) for o, k in zip(D["le_roff"], n)])
    rep = np.repeat(np.arange(len(le)), n)
    R1 = np.nan_to_num(np.asarray(D["rain"], float)[idx])
    bin_le = np.repeat(np.arange(len(le)), n // 5)
    target = D["y5"] - D["le_c"][bin_le]
    scored = (D["w5"] > 0) & np.isfinite(target)
    events = sorted(le.event_id.unique())
    ev_le = le.event_id.map({e: i for i, e in enumerate(events)}).to_numpy()
    eq = json.loads(P.B_EQUATION.read_text(encoding="utf-8"))
    return {"D": D, "le": le, "curves": curves, "band": np.asarray(band_idx, np.int64), "R1": R1, "rep": rep, "bin_le": bin_le,
            "target": target, "scored": scored, "events": events, "ev_le": ev_le, "ev_bin": ev_le[bin_le], "eq": eq,
            "starts": D["le_boff1"].astype(np.int64), "lens": n.astype(np.int64), "L": D["le_L"].astype(float)}


def frozen_parts(st):
    """(theta, F increments per band, P increments per band) of the frozen equation."""
    prob = ksurf.prepare(st["D"], st["band"], st["curves"], st["D"]["rain"], 0)
    m, x = B.ours_model(prob)
    nb = len(st["curves"])
    return np.asarray(st["eq"]["drying"]["theta"], float), x[m.sl["F"]].reshape(nb, -1), x[m.sl["P"]].reshape(nb, -1)


def to5(a):
    return a.reshape(-1, 5).mean(axis=1)


def path5(st, which, incP):
    bm = st["band"][st["rep"]]
    p = np.empty(len(st["R1"]))
    for b, cv in enumerate(st["curves"]):
        m = bm == b
        if which == "A":
            p[m] = ksurf.full_curve(cv, "P", incP[b], st["R1"][m])
        else:
            p[m] = CW["path_b_c_band"][cv.band] * st["R1"][m]
    return to5(st["L"][st["rep"]] * p)


# =================================================================================================== the WAA laws
@njit(cache=True)
def two_store(R, starts, lens, band, A, b, kfw, kfd, ks, phi):
    """Part 1's two stores (exact per-minute steps of WAAEquation, the output after the update), u = A_band (R / 5)^b."""
    W = np.zeros(len(R))
    for e in range(len(starts)):
        wf, ws = 0.0, 0.0
        for t in range(lens[e]):
            i = starts[e] + t
            r = R[i]
            u = A[band[e]] * (r / 5.0) ** b if r > 0 else 0.0
            uf, us = phi * u, (1.0 - phi) * u
            wf = uf + (wf - uf) * (kfw if uf > wf else kfd)
            ws = us + (ws - us) * ks
            W[i] = wf + ws
    return W


def waa_param_model(st, name, p):
    """W per 5 min of the three parametric laws for a parameter vector p (log-parametrized where positive)."""
    bm = st["band"]
    if name == "Part 1 equation":
        A = np.exp(p[:2])
        W1 = two_store(st["R1"], st["starts"], st["lens"], bm, A, float(p[2]), np.exp(-1 / np.exp(p[3])), np.exp(-1 / np.exp(p[4])),
                       np.exp(-1 / CW["part1"]["tau_s_min"]), CW["part1"]["phi"])
    elif name == "Schleiss":
        u = lit.schleiss_unit(st["R1"], st["starts"], st["lens"], float(np.exp(p[2])), float(np.exp(p[3])), CW["schleiss"]["wet_mmh"], 1.0)
        W1 = np.exp(p[:2])[bm[st["rep"]]] * u
    else:                                                       # Pastorek
        W1 = np.exp(p[:2])[bm[st["rep"]]] * lit.pastorek_unit(st["R1"], float(p[2]), CW["pastorek"]["d"])
    return to5(W1)


def start_bounds(name):
    lit_f = B.lit_fits()
    if name == "Part 1 equation":
        e = json.loads(P.P1_EQUATION.read_text(encoding="utf-8"))
        x0 = [np.log(e["A_band_db"]["28 GHz"]["value"]), np.log(e["A_band_db"]["38 GHz"]["value"]), e["b"]["value"],
              np.log(e["dynamics_fixed"]["tau_f_wet"]), np.log(e["dynamics_fixed"]["tau_f_dry"])]
        bnd = [(np.log(0.01), np.log(50))] * 2 + [(0.05, 2.0)] + [(np.log(0.5), np.log(600))] * 2
    elif name == "Schleiss":
        f = lit_f["Schleiss b"]
        x0 = [np.log(max(f["amp"][0], 0.02)), np.log(max(f["amp"][2], 0.02)), np.log(f["shape"][0]), np.log(f["shape"][1])]
        bnd = [(np.log(0.01), np.log(50))] * 2 + [(np.log(3.0), np.log(600))] * 2
    else:
        f = lit_f["Pastorek b"]
        x0 = [np.log(max(f["amp"][0], 0.02)), np.log(max(f["amp"][2], 0.02)), f["shape"][0]]
        bnd = [(np.log(0.01), np.log(50))] * 2 + [(0.1, 2.0)]
    return np.clip(np.asarray(x0, float), [b[0] for b in bnd], [b[1] for b in bnd]), bnd


def f_columns(st):
    """No memory: W5 = base + C @ incF (affine in the F extension increments)."""
    bm = st["band"][st["rep"]]
    nk = len(ksurf.config()["f_knots_mmh"])
    base = np.zeros(len(st["R1"]))
    cols = []
    for b, cv in enumerate(st["curves"]):
        m = bm == b
        z = np.zeros(nk)
        f0 = ksurf.full_curve(cv, "F", z, st["R1"][m])
        base[m] = f0
        for j in range(nk):
            e = z.copy()
            e[j] = 1.0
            c = np.zeros(len(st["R1"]))
            c[m] = ksurf.full_curve(cv, "F", e, st["R1"][m]) - f0
            cols.append(to5(c))
    return to5(base), np.column_stack(cols)


# =================================================================================================== one fit
def fit_fold(st, ctx, name, path, k, maxiter):
    """Fit model `name` on every event but k; return its parameters and W5 for ALL link-events of the setting."""
    res = ctx["res"][path]
    train = st["scored"] & (st["ev_bin"] != k)
    t0 = time.time()
    if name == "no WAA":
        return {"params": [], "W5": np.zeros(len(res)), "converged": True, "loss": float(huber(res[train]).mean()), "minutes": 0.0}
    if name == "no memory":
        base, C = ctx["fcols"]
        x, loss, conv, nit = B.huber_fit(base, C, res, train.astype(float), ctx["incF0"].ravel(), [(0.0, ksurf.config()["f_inc_max_db"])] * C.shape[1], DELTA)
        return {"params": x.tolist(), "W5": base + C @ x, "converged": conv, "loss": loss, "minutes": (time.time() - t0) / 60}
    if name == "ours":
        prob = ctx["prob"]
        prob["target"] = res.astype(float)
        fm = F8.FormModel(prob, ctx["form"], (np.array([np.log(1.0 / CW["ours"]["tau_wet_min"])]), np.zeros(0)))
        z0 = np.r_[ctx["theta0"], ctx["incF0"].ravel(), np.zeros(fm.nP)]
        use = (st["ev_le"] != k).astype(float)
        z = fm.fit(z0, use, maxiter=maxiter)
        return {"params": z[:fm.form.p + fm.nF].tolist(), "W5": fm.simulate(z)[1], "converged": fm.last["converged"], "loss": fm.last["loss"],
                "minutes": (time.time() - t0) / 60}
    x0, bnd = start_bounds(name)

    def f(p):
        return float(huber(waa_param_model(st, name, p)[train] - res[train]).mean())
    r = minimize(f, x0, method="L-BFGS-B", bounds=bnd, options={"maxiter": maxiter or CW["maxiter"]})
    return {"params": r.x.tolist(), "W5": waa_param_model(st, name, r.x), "converged": bool(r.success), "loss": float(r.fun),
            "minutes": (time.time() - t0) / 60}


def a_link_fit(st, W5, res, k):
    """Per-link factor on W (Part 1's ridge), fitted on the training events of fold k, one link at a time."""
    pen = 0.5 * (DELTA / CW["a_prior_sd"]) ** 2
    train = st["scored"] & (st["ev_bin"] != k)
    link_bin = st["le"].cml_id.to_numpy()[st["bin_le"]]
    a = {}
    for c in np.unique(link_bin[train]):
        m = train & (link_bin == c)
        w, r = W5[m], res[m]
        o = minimize_scalar(lambda la: float(huber(np.exp(la) * w - r).sum() + pen * la * la), bounds=tuple(CW["log_a_bounds"]), method="bounded")
        a[c] = float(np.exp(o.x))
    return a


# =================================================================================================== the run
def context(st):
    theta0, incF0, incP0 = frozen_parts(st)
    ctx = {"theta0": theta0, "incF0": incF0, "incP0": incP0, "fcols": f_columns(st)}
    ctx["path"] = {p: path5(st, p, incP0) for p in PATHS}
    ctx["res"] = {p: st["target"] - ctx["path"][p] for p in PATHS}
    prob = ksurf.prepare(st["D"], st["band"], st["curves"], st["D"]["rain"], 0)
    prob["Pfix"][:] = 0.0                                       # the kernel's path term is off: the fixed path is in the target
    prob["jP"][:] = -1
    prob["w5"] = st["scored"].astype(float)
    ctx["prob"] = prob
    rk = pickle.loads((P.CACHE / "a8_name" / "rank.pkl").read_bytes())     # step 8's form (read only)
    spec = rk["fits"][CW["ours"]["form_name"]]["spec"]
    ctx["form"] = F8.from_sympy(spec["key"], spec["expr"])
    return ctx


def run_all(st, ctx, cache, maxiter=None, folds=None):
    out = {}
    for path in PATHS:
        for name in MODELS:
            for k in (folds if folds is not None else range(len(st["events"]))):
                f = cache / f"{path}_{name.replace(' ', '_')}_{k:02d}.pkl"
                if f.exists():
                    out[(path, name, k)] = pickle.loads(f.read_bytes())
                    continue
                free = psutil.virtual_memory().available / GB
                if free < CW["min_free_gb"]:
                    raise SystemExit(f"free memory {free:.1f} GB below {CW['min_free_gb']} GB: stopped before {path} {name} fold {k}")
                r = fit_fold(st, ctx, name, path, k, maxiter)
                held = st["ev_bin"] == k
                keep = {"params": r["params"], "converged": r["converged"], "loss": r["loss"], "minutes": r["minutes"], "W5_held": r["W5"][held]}
                if path == "A" and name in WITH_A:
                    a = a_link_fit(st, r["W5"], ctx["res"][path], k)
                    lb = st["le"].cml_id.to_numpy()[st["bin_le"]][held]
                    keep["a_held"] = np.array([a.get(c, 1.0) for c in lb])
                    keep["a_links"] = len(a)
                f.write_bytes(pickle.dumps(keep))
                out[(path, name, k)] = keep
                log(f"   ({path} {name} fold {k:2d} = {st['events'][k]}) loss {r['loss']:.5f}; converged {r['converged']}; {r['minutes']:.1f} min; "
                    f"free {free:.1f} GB")
    return out


def assemble(st, ctx, fits):
    """Held-out W5 and totals per path and model (each bin is held out in exactly one fold); plus the a_link variant."""
    pred, waa = {}, {}
    for path in PATHS:
        for name in MODELS:
            W = np.full(len(st["target"]), np.nan)
            for k in range(len(st["events"])):
                W[st["ev_bin"] == k] = fits[(path, name, k)]["W5_held"]
            waa[(path, name)] = W
            pred[(path, name)] = W + ctx["path"][path]
            if path == "A" and name in WITH_A:
                Wa = np.full(len(st["target"]), np.nan)
                for k in range(len(st["events"])):
                    h = st["ev_bin"] == k
                    Wa[h] = fits[(path, name, k)]["a_held"] * fits[(path, name, k)]["W5_held"]
                waa[("A+a", name)] = Wa
                pred[("A+a", name)] = Wa + ctx["path"]["A"]
    for name in MODELS:                                         # the a_link variant: the others unchanged
        if ("A+a", name) not in pred:
            pred[("A+a", name)], waa[("A+a", name)] = pred[("A", name)], waa[("A", name)]
    return pred, waa


# =================================================================================================== metrics
def metrics(st, pred):
    E = len(st["events"])
    rng = np.random.default_rng(CW["seed"])
    mult = np.array([np.bincount(rng.integers(0, E, E), minlength=E) for _ in range(CW["boot_reps"])], float)
    lo_q, hi_q = CW["interval"]
    sc = st["scored"] & np.logical_and.reduce([np.isfinite(v) for v in pred.values()])
    Lb = st["L"][st["bin_le"]]
    rm, pr = [], []
    for path in ("A", "B", "A+a"):
        for gn, gm in (("all links", np.ones(len(sc), bool)), ("links < 1 km", Lb < 1.0)):
            msk = sc & gm
            bt = {}
            for mn in MODELS:
                r2 = (pred[(path, mn)] - st["target"]) ** 2
                S_ = np.bincount(st["ev_bin"][msk], weights=r2[msk], minlength=E)
                N_ = np.bincount(st["ev_bin"][msk], minlength=E).astype(float)
                bt[mn] = np.sqrt(mult @ S_ / np.maximum(mult @ N_, 1))
                rm.append({"path": path, "links": gn, "model": mn, "RMSE (dB)": float(np.sqrt(S_.sum() / N_.sum())), "values": int(N_.sum()),
                           "90 % (dB)": f"{np.percentile(bt[mn], lo_q):.4f} to {np.percentile(bt[mn], hi_q):.4f}"})
            ro = rm[-len(MODELS)]["RMSE (dB)"]
            for j, mn in enumerate(MODELS[1:], 1):
                lo, hi = np.percentile(bt[mn] - bt["ours"], [lo_q, hi_q])
                pr.append({"path": path, "links": gn, "model": mn, "RMSE(model) - RMSE(ours) (dB)": rm[-len(MODELS) + j]["RMSE (dB)"] - ro,
                           "90 % lo": lo, "90 % hi": hi, "reading": f"{mn} better than ours" if hi < 0 else (f"ours better than {mn}" if lo > 0 else "tie")})
    rm, pr = pd.DataFrame(rm), pd.DataFrame(pr)

    def paired_abs(rows, name, vals, ev_of, key):
        """vals[model] = per-item errors; median and paired |median| comparison over event draws."""
        rr = {"item": name, **key}
        drw = {}
        for mn in MODELS:
            e0 = float(np.median(vals[mn])) if len(vals[mn]) else np.nan
            d = []
            for w in mult:
                rep = np.repeat(np.arange(len(ev_of)), w[ev_of].astype(int))
                d.append(np.median(vals[mn][rep]) if len(rep) else np.nan)
            drw[mn] = np.array(d)
            lo, hi = np.nanpercentile(drw[mn], [lo_q, hi_q]) if len(vals[mn]) else (np.nan, np.nan)
            rr[mn] = f"{e0:+.3f} ({lo:+.3f} to {hi:+.3f})"
        for mn in MODELS[1:]:
            lo, hi = np.nanpercentile(np.abs(drw[mn]) - np.abs(drw["ours"]), [lo_q, hi_q]) if len(ev_of) else (np.nan, np.nan)
            rr[f"{mn}: |e| - |e ours|"] = f"{lo:+.3f} to {hi:+.3f} ({'closer than ours' if hi < 0 else ('ours closer' if lo > 0 else 'tie')})"
        rows.append(rr)

    # after-rain (path A)
    D = st["D"]
    lags = CW["after_lags_min"]
    idx, _ = a7r.dry_values(D, fake6.window(D, np.asarray(D["rain"], float)), lags, CW["rain_before_min"], s1.S1["dry_below_mmh"], np.ones(len(st["le"])))
    ar = []
    for i, lag in enumerate(lags):
        b = idx[:, i][idx[:, i] >= 0]
        b = b[sc[b]]
        vals = {mn: pred[("A", mn)][b] - st["target"][b] for mn in MODELS}
        paired_abs(ar, "after-rain error", vals, st["ev_bin"][b], {"lag (min)": lag, "windows": int(len(b)),
                                                                   "measured median (dB)": float(np.median(st["target"][b])) if len(b) else np.nan})
    # peak error (path A)
    pk = []
    pk_bins = []
    o5 = D["le_boff5"]
    n5 = D["le_n"] // 5
    for e in range(len(st["le"])):
        bb = o5[e] + np.nonzero(sc[o5[e]:o5[e] + n5[e]])[0]
        if len(bb) and np.max(st["target"][bb]) >= CW["peak_min_db"]:
            pk_bins.append(bb[np.argmax(st["target"][bb])])
    pk_bins = np.array(pk_bins, int)
    for cname, cm in (("links < 1 km", Lb[pk_bins] < 1.0), ("links > 3 km", Lb[pk_bins] > 3.0)):
        b = pk_bins[cm]
        vals = {mn: pred[("A", mn)][b] - st["target"][b] for mn in MODELS}
        paired_abs(pk, "peak error", vals, st["ev_bin"][b], {"links": cname, "link-events": int(len(b)),
                                                             "measured peak median (dB)": float(np.median(st["target"][b])) if len(b) else np.nan})
    # ranking
    order = {p: list(rm[(rm.path == p) & (rm.links == "all links")].sort_values("RMSE (dB)").model) for p in ("A", "B")}
    return {"rmse": rm, "paired": pr, "after": pd.DataFrame(ar), "peak": pd.DataFrame(pk), "order": order, "sc": sc}


# =================================================================================================== figures
def select(st, pred, sc, Lmask, rule):
    le = st["le"]
    n5 = st["D"]["le_n"] // 5
    dry_after = n5 - st["D"]["le_last"] - 1
    o5 = st["D"]["le_boff5"]
    rm = np.full(len(le), np.nan)
    for e in range(len(le)):
        m = sc[o5[e]:o5[e] + n5[e]]
        if m.any():
            rm[e] = np.sqrt(np.mean((pred[("A", "ours")][o5[e]:o5[e] + n5[e]][m] - st["target"][o5[e]:o5[e] + n5[e]][m]) ** 2))
    cand = np.nonzero(Lmask & (dry_after >= CW["dry_min_bins"]) & np.isfinite(rm))[0]
    med = float(np.median(rm[cand]))
    pick = int(cand[np.argmin(np.abs(rm[cand] - med))])
    rank = int(np.sum(rm[cand] < rm[pick])) + 1
    row = {"rule": rule, "candidates": len(cand), "median held-out RMSE ours (dB)": med, "event": le.event_id[pick], "sublink": le.sublink_id[pick],
           "L (km)": float(le.L_km[pick]), "band": le.band_name[pick], "held-out RMSE ours (dB)": float(rm[pick]), "rank": f"{rank} of {len(cand)}"}
    log(f"   ({rule}) {len(cand)} candidates; median {med:.4f} dB; selected {row['event']} {row['sublink']} (L {row['L (km)']:.2f} km): "
        f"RMSE {row['held-out RMSE ours (dB)']:.4f}, rank {row['rank']}")
    return pick, row


def figure(st, pred, waa, sc, e, rule):
    le, D = st["le"], st["D"]
    o5, n = int(D["le_boff5"][e]), int(D["le_n"][e] // 5)
    sl = slice(o5, o5 + n)
    t = np.arange(n) * 5 + 2.5
    R = to5(st["R1"][int(D["le_boff1"][e]):int(D["le_boff1"][e]) + 5 * n])
    wet = R > CW["wet_mmh"]
    cnt = (D["w5"][sl] > 0) & np.isfinite(st["target"][sl])
    meas = np.where(cnt, st["target"][sl], np.nan)
    tot = {m: pred[("A", m)][sl] for m in MODELS}
    w = {m: waa[("A", m)][sl] for m in MODELS}
    v = np.concatenate([a[np.isfinite(a)] for a in [meas] + list(tot.values()) + list(w.values())])
    pad = 0.06 * (v.max() - v.min() if v.max() > v.min() else 1.0)
    yl = (v.min() - pad, v.max() + pad)
    fig, axes = plt.subplots(7, 1, figsize=(8.5, 13), sharex=True)
    for i, ax in enumerate(axes):
        if i > 0:
            for kk in np.nonzero(wet)[0]:
                ax.axvspan(kk * 5, kk * 5 + 5, color=SPAN, lw=0, zorder=0)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(axis="y", color="#ececec", lw=0.6)
        ax.tick_params(labelsize=8)
        ax.set_ylim(*yl)
        ax.set_ylabel("dB", fontsize=8)
    ax2 = axes[0].twinx()
    for kk in np.nonzero(wet)[0]:
        ax2.axvspan(kk * 5, kk * 5 + 5, color=SPAN, lw=0, zorder=0)
    ax2.bar(t, R, width=4.6, color="#6baed6", alpha=0.8, lw=0, zorder=1, label="radar rain R = radar / g_e (mm/h)")
    ax2.set_ylabel("R (mm/h)", fontsize=8, color="#4a7fa8")
    ax2.tick_params(labelsize=8, colors="#4a7fa8")
    ax2.spines["top"].set_visible(False)
    ax2.set_ylim(0, max(1.0, 1.15 * float(R.max())))
    axes[0].set_zorder(ax2.get_zorder() + 1)
    axes[0].patch.set_visible(False)
    axes[0].plot(t, meas, color="black", lw=1.6, label="measured A - c (5 min)")
    axes[0].set_title("measured A - c, with the radar rain input", fontsize=9, loc="left")
    h1, l1 = axes[0].get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    axes[0].legend(h1 + h2, l1 + l2, fontsize=7, loc="upper right", frameon=False)
    for ax, m in zip(axes[1:], MODELS):
        ax.plot(t, meas, color=GREY, lw=1.0, ls="--", label="measured")
        ax.plot(t, tot[m], color=COLOR[m], lw=1.8, label="held-out total (WAA + path)")
        ax.plot(t, w[m], color=COLOR[m], lw=1.0, ls=":", label="held-out WAA part")
        s_ = sc[sl]
        r = np.sqrt(np.mean((tot[m][s_] - meas[s_]) ** 2)) if s_.any() else np.nan
        ax.set_title(f"{m}  -  held-out RMSE {r:.3f} dB on this link-event", fontsize=9, loc="left")
    axes[1].legend(fontsize=7, loc="upper right", frameon=False)
    axes[-1].set_xlabel(f"minutes from the start of the link-event window (grey: radar rain input R > {CW['wet_mmh']} mm/h)", fontsize=9)
    fig.suptitle(f"TRAIN cross-validation, same path law (path A)\n{le.event_id[e]}   sub-link {le.sublink_id[e]}   L = {le.L_km[e]:.2f} km   "
                 f"{le.band_name[e]}   rule {rule}   (held-out fold: this event)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    stem = P.FIGURES / f"w_cv_{rule}_{le.event_id[e]}_{le.sublink_id[e]}"
    save(fig, stem.with_suffix(".png"), dpi=300)
    save(fig, stem.with_suffix(".pdf"))
    plt.close(fig)
    return stem


# =================================================================================================== report
def finish(st, ctx, fits, t0):
    pred, waa = assemble(st, ctx, fits)
    np.savez_compressed(CACHE / "preds.npz", **{f"total|{p}|{m}": v for (p, m), v in pred.items()}, **{f"waa|{p}|{m}": v for (p, m), v in waa.items()},
                        target=st["target"], bin_le=st["bin_le"])
    M = metrics(st, pred)
    params = pd.DataFrame([{"path": p, "model": m, "fold": k, "held-out event": st["events"][k], "loss": v["loss"], "converged": v["converged"],
                            "minutes": v["minutes"], "params": json.dumps(np.round(v["params"], 6).tolist()), "a_links": v.get("a_links", "")}
                           for (p, m, k), v in fits.items()])
    Lk = st["L"]
    rows, figs = [], []
    for rule, Lmask in (("short", Lk < 1.0), ("long", Lk > 3.0)):
        e, row = select(st, pred, M["sc"], Lmask, rule)
        stem = figure(st, pred, waa, M["sc"], e, rule)
        row["figure"] = rel(stem.with_suffix(".png"))
        rows.append(row)
    sel = pd.DataFrame(rows)
    for tab, f in ((params, P.WCV_PARAMS), (M["rmse"].round(5), P.WCV_RMSE), (M["paired"].round(5), P.WCV_PAIRED), (M["after"], P.WCV_AFTER),
                   (M["peak"], P.WCV_PEAK), (sel.round(4), P.WCV_SELECTION)):
        save(tab, f, index=False)
    same = M["order"]["A"] == M["order"]["B"]
    pa, pb = M["paired"][M["paired"].path == "A"], M["paired"][M["paired"].path == "B"]
    changed = [f"{r.model} ({r.links}): A '{r.reading}' vs B '{q.reading}'" for (_, r), (_, q) in zip(pa.iterrows(), pb.iterrows()) if r.reading != q.reading]
    L = ["# WAA-only comparison: the same path law for every model, leave-one-event-out on TRAIN\n",
         "Written by `research/waa_ksurface/w_cv.py` (decisions.md, \"WAA-only comparison\", logged before any code was run). 14 TRAIN learning events, "
         f"{len(st['le'])} link-events (28 / 38 GHz); held-out predictions only. VAL and TEST not opened. Disclosure: the forms were chosen on these "
         "events (ours) or on TRAIN + VAL (Part 1); only parameters are cross-validated.\n",
         "## RMSE of the 5 min A - c (held-out; paired event bootstrap, 1000 draws, seed 4400)\n", md_table(M["rmse"].round(4)), "",
         "### Paired RMSE(model) - RMSE(ours) and readings (path A = main; B = Part 1's straight path; A+a = competitors with a_link)\n",
         md_table(M["paired"].round(4)), "",
         f"**Ranking by RMSE (all links):** path A: {', '.join(M['order']['A'])}; path B: {', '.join(M['order']['B'])}. "
         f"**{'The same' if same else 'NOT the same'} under A and B.**" + (" Readings that change: " + "; ".join(changed) + "." if changed else " No reading changes.") + "\n",
         "## After-rain error (path A; model - measured, median with 90 %; paired |e| comparison)\n", md_table(M["after"]), "",
         f"## Peak error (path A; link-events with a measured peak >= {CW['peak_min_db']} dB)\n", md_table(M["peak"]), "",
         "## Figures (path A, held-out)\n", md_table(sel.round(4)), "",
         f"Saved predictions (total and WAA part per model, path and held-out bin): `{rel(CACHE / 'preds.npz')}` (cache, not committed). "
         f"Fitted parameters per fold: `{rel(P.WCV_PARAMS)}`. Run time {(time.time() - t0) / 60:.1f} min (this pass).\n"]
    save("\n".join(L) + "\n", P.WCV_REPORT)
    log(f"report: {rel(P.WCV_REPORT)}")


def main():
    t0 = time.time()
    st = setting()
    ctx = context(st)
    log(f"setting: {len(st['le'])} link-events, {len(st['events'])} events, {int(st['scored'].sum())} scored values ({(time.time() - t0) / 60:.1f} min)")
    if "--smoke" in sys.argv:
        c = CACHE / "smoke"
        c.mkdir(parents=True, exist_ok=True)
        run_all(st, ctx, c, maxiter=3, folds=[0])
        log("=== smoke ended (fold 0, 3 iterations): every model fits and predicts ===")
        return
    CACHE.mkdir(parents=True, exist_ok=True)
    fits = run_all(st, ctx, CACHE) if "--report" not in sys.argv else {
        (p, m, k): pickle.loads((CACHE / f"{p}_{m.replace(' ', '_')}_{k:02d}.pkl").read_bytes()) for p in PATHS for m in MODELS for k in range(len(st["events"]))}
    log(f"=== every fit done ({(time.time() - t0) / 60:.1f} min) ===")
    finish(st, ctx, fits, t0)
    log(f"=== WAA-only comparison ended ({(time.time() - t0) / 60:.1f} min) ===")


if __name__ == "__main__":
    main()
