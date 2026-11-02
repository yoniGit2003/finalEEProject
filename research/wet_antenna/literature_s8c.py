"""wet_antenna step S8c: our final WAA equation against two literature WAA models (a POST-FREEZE COMPARISON).

TEST was read in S8 / S8b. Here the literature models are fitted on TRAIN only, our equation is not refitted, and no choice
about our equation depends on this step.

  The task (as S8): every model gets only the rain input V1 (radar / g_e), the link length and band, and predicts the 5 min
  attenuation A - c = W + L c_band R(t - d). The same link-events, baseline, offset c and lag d (1 min) as S8.

  Pastorek et al. 2022 (KR-alt):  W = C (1 - exp(-d R^zeta)), d = 0.1. No memory, an upper limit C.
  Schleiss et al. 2013:           a model with memory. While the radar says wet (V1 > 0.1 mm/h) W rises toward W_max
                                  (pycomlink's step); when dry it decays (an exponential step that replaces pycomlink's use of
                                  the measured attenuation, which would hand the model the answer).
  Two versions of each: (a) the published / pycomlink parameters, a_link = 1, only c_band fitted on TRAIN (Schleiss a = "pycomlink
  defaults plus an assumed drying time equal to the wetting time (15 min), not taken from the paper"); (b) everything fitted on
  TRAIN (amplitude per band, shared shape, c_band, a_link with the ridge of our model). The verdict that counts is version (b).
  Schleiss b: a grid over (tau_wet, tau_dry), then Nelder-Mead from its best point. Every outer step is cached at once.

  The fit: for given shape parameters the model is linear in the amplitude and in c_band, so an inner L-BFGS-B fits
  (amplitude per band, c_band, log a_link) with the Huber loss and the a_link ridge of our model, and an outer search moves
  the shape parameters (zeta; tau_wet, tau_dry).

Settings: config `wet_antenna.literature_s8c` (written before the run). Caches: research/.cache/wet_antenna/literature_s8c/.
Run from the repository root:
    python -m research.wet_antenna.literature_s8c
"""
import hashlib
import json
import pickle
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numba import njit
from scipy.optimize import minimize, minimize_scalar

from cmlpipe.api import CMLData
from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import length_test as LTL
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.event_factors import rescaled
from research.wet_antenna.final_model import CFG, DATA, WA, predict, s6_models
from research.wet_antenna.test_s8b import b_models, ev_sums, pooled
from research.wet_antenna.validation import LB, boot_items, ci, dry_windows, ev_draws, radar_max, w_minutes

S = WA["literature_s8c"]
F_ = S["fit"]
CACHE = P.S8C_CACHE
PADM = sid.PADM
LAG = 1                                                           # d of our equation (checked against the fit in main)
VL = "VAL (used, optimistic)"
LIT = ("Schleiss a", "Schleiss b", "Pastorek a", "Pastorek b")
ORDER = ["ours", "table F", "B0", "B1", "B2", *LIT]
STATUS = ("**Status.** TEST was read in S8 / S8b. This is a **post-freeze comparison**. The literature models are fitted on TRAIN only, our "
          "equation is not refitted, and no choice about our equation depends on this step.")

# The pycomlink 0.6.0 code these models follow (processing/wet_antenna.py), quoted in the report:
PYCOMLINK_PASTOREK = '''def waa_pastorek_2021(R, A_max=14, zeta=0.55, d=0.1):
    R = np.asanyarray(R)
    waa = A_max * (1 - np.exp(-d * (R**zeta)))
    # Assure that numeric inaccuracy does not lead to waa > 0 for R == 0
    waa[R == 0] = 0
    return waa'''
PYCOMLINK_SCHLEISS = '''def _numba_waa_schleiss_2013(rsl, baseline, wet, waa_max, delta_t, tau):
    waa = np.zeros_like(rsl, dtype=np.float64)
    A = rsl - baseline
    for i in range(1, len(rsl)):
        if wet[i] == True:
            waa[i] = min(A[i], waa_max, waa[i - 1] + (waa_max - waa[i - 1]) * 3 * delta_t / tau)
        else:
            waa[i] = min(A[i], waa_max)
    return waa'''
OURS_SCHLEISS = '''for every link-event, W = 0 before its first minute; for every minute i:
    if V1(i - d) > 0.1 mm/h:  W[i] = min(W_max, W[i-1] + (W_max - W[i-1]) * 3 * dt / tau_wet)   # pycomlink's step, without A[i]
    else:                     W[i] = W[i-1] - W[i-1] * 3 * dt / tau_dry                          # replaces min(A[i], W_max)
with dt = 1 min; the prediction is a_link W + L c_band V1(i - d) + c'''


# =================================================================================================== the two literature models (unit amplitude)
def pastorek_unit(R, zeta, d=S["pastorek"]["d"]):
    """KR-alt with C = 1: the pycomlink expression A_max (1 - exp(-d R^zeta)), 0 where R == 0."""
    R = np.asarray(R, float)
    w = 1.0 * (1 - np.exp(-d * (np.maximum(R, 0.0) ** zeta)))
    w[R <= 0] = 0.0
    return w


@njit(cache=True)
def schleiss_unit(R1, starts, lens, tau_wet, tau_dry, thr, dt):
    """Schleiss with W_max = 1 (W scales exactly with W_max), driven by the radar wet flag R1 > thr."""
    g = np.zeros(len(R1))
    for e in range(len(starts)):
        w = 0.0
        for t in range(lens[e]):
            i = starts[e] + t
            if R1[i] > thr:
                w = min(1.0, w + (1.0 - w) * 3.0 * dt / tau_wet)
            else:
                w = w - w * 3.0 * dt / tau_dry
            g[i] = w
    return g


class Set:
    """One dataset prepared for the literature models: the lagged rain per minute, the 5 min path input, the bins that count."""

    def __init__(self, D, le, train_ids=None):
        self.D, self.le = D, le
        n, roff = D["le_n"], D["le_roff"]
        self.starts, self.lens = D["le_boff1"].astype(np.int64), n.astype(np.int64)
        idx = np.concatenate([roff[e] + PADM - LAG + np.arange(n[e]) for e in range(len(n))])
        self.R1 = D["rain"][idx].astype(float)                    # V1(t - d) per minute, aligned with y1
        self.p5 = D["le_L"][D["bin_le"]] * self.R1.reshape(-1, 5).mean(axis=1)
        self.m = D["w5"] > 0
        self.band5 = D["le_band"][D["bin_le"]]
        self.c5 = D["le_c"][D["bin_le"]]
        ids = LTL.link_ids(D, le)
        self.link5 = D["le_link"][D["bin_le"]]
        self.ids = ids
        self.nl, self.nb = int(D["nl"]), int(D["nb"])
        m = self.m                                                      # the counted bins, once (used by every inner fit)
        self.pm, self.ym, self.cm, self.bm, self.lm = self.p5[m], D["y5"][m], self.c5[m], self.band5[m], self.link5[m]

    def unit(self, kind, shape):
        """Unit-amplitude W per minute of a model."""
        if kind == "pastorek":
            return pastorek_unit(self.R1, shape[0])
        sc = S["schleiss"]
        return schleiss_unit(self.R1, self.starts, self.lens, float(shape[0]), float(shape[1]), sc["wet_mmh"], sc["dt_min"])

    def a_le(self, a_by_id):
        """a_link per link index of this dataset (1 for a link not in the fit)."""
        return np.array([a_by_id.get(c, 1.0) for c in self.ids])

    def predict(self, fit):
        g1 = self.unit(fit["kind"], fit["shape"])
        a = self.a_le(fit["a"])
        amp = np.asarray(fit["amp"], float)
        W1 = a[self.D["le_link"]].repeat(self.D["le_n"]) * amp[self.D["le_band"]].repeat(self.D["le_n"]) * g1
        pred5 = W1.reshape(-1, 5).mean(axis=1) + np.asarray(fit["c_band"])[self.band5] * self.p5 + self.c5
        return pred5, W1


def inner(T, g5, amp, c, la, free_amp, free_a):
    """Fit (amplitude per band, c_band, log a_link) for one unit curve g5 (5 min means): Huber loss + the a_link ridge of our model."""
    g, p, y, cc, band, link = g5[T.m], T.pm, T.ym, T.cm, T.bm, T.lm
    nb, nl, dl = T.nb, T.nl, F_["huber_delta_db"]
    pen = 0.5 * (dl / F_["a_prior_sd"]) ** 2

    def split(x):
        k = 0
        A = x[k:k + nb] if free_amp else np.asarray(amp, float)
        k += nb if free_amp else 0
        C = x[k:k + nb]
        k += nb
        L = x[k:k + nl] if free_a else np.zeros(nl)
        return A, C, L

    def fun(x):
        A, C, L = split(x)
        ea = np.exp(L)[link]
        Wb = ea * A[band] * g
        r = Wb + C[band] * p + cc - y
        ar = np.abs(r)
        loss = float(np.where(ar <= dl, 0.5 * r * r, dl * (ar - 0.5 * dl)).sum() + pen * L @ L)
        psi = np.clip(r, -dl, dl)
        gr = []
        if free_amp:
            gr.append(np.bincount(band, psi * ea * g, nb))
        gr.append(np.bincount(band, psi * p, nb))
        if free_a:
            gr.append(np.bincount(link, psi * Wb, nl) + 2 * pen * L)
        return loss, np.concatenate(gr)

    x0 = np.r_[np.asarray(amp, float) if free_amp else [], c, la if free_a else []]
    bnd = ([tuple(F_["amp_bounds_db"])] * nb if free_amp else []) + [tuple(F_["c_band_bounds"])] * nb + ([tuple(F_["log_a_bounds"])] * nl if free_a else [])
    r = minimize(fun, np.clip(x0, [b[0] for b in bnd], [b[1] for b in bnd]), jac=True, method="L-BFGS-B", bounds=bnd,
                 options={"maxiter": F_["inner"]["maxiter"], "ftol": 1e-12, "gtol": 1e-7})
    A, C, L = split(r.x)
    return {"amp": A, "c_band": C, "la": L, "loss": float(r.fun), "converged": bool(r.success), "nit": int(r.nit)}


def fit_model(T, kind, version, c_start, log):
    """One literature model version on TRAIN (cached)."""
    name = f"{'Schleiss' if kind == 'schleiss' else 'Pastorek'} {version}"
    file = CACHE / f"fit_{name.replace(' ', '_')}.pkl"
    if file.exists():
        return pickle.loads(file.read_bytes())
    t0 = time.time()
    df = S[kind]["default"]
    amp0 = np.full(T.nb, df["C_db"] if kind == "pastorek" else df["W_max_db"])
    shape0 = [df["zeta"]] if kind == "pastorek" else [df["tau_wet_min"], df["tau_dry_min"]]
    if version == "a":
        g5 = T.unit(kind, shape0).reshape(-1, 5).mean(axis=1)
        r = inner(T, g5, amp0, c_start, np.zeros(T.nl), False, False)
        shape, evals = shape0, 1
    else:
        steps_file = CACHE / f"steps_{name.replace(' ', '_')}.pkl"          # every outer evaluation, saved at once
        done = pickle.loads(steps_file.read_bytes()) if steps_file.exists() else {}
        state = {"amp": amp0, "c": np.asarray(c_start, float), "la": np.zeros(T.nl), "best": None, "n": 0}
        for key, r in done.items():                                        # resume: the best saved step is the warm start
            if state["best"] is None or r["loss"] < state["best"][1]["loss"]:
                state["best"] = (list(key), r)
                state["amp"], state["c"], state["la"] = r["amp"], r["c_band"], r["la"]

        def outer(shape):
            key = tuple(round(float(s), 6) for s in shape)
            if key in done:
                r = done[key]
            else:
                g5 = T.unit(kind, shape).reshape(-1, 5).mean(axis=1)
                r = inner(T, g5, state["amp"], state["c"], state["la"], True, True)
                done[key] = r
                steps_file.write_bytes(pickle.dumps(done))
                log(f"   {name} step {len(done)}: shape {np.round(shape, 4).tolist()} loss {r['loss']:.2f} ({(time.time() - t0) / 60:.1f} min)")
            state["n"] += 1
            if state["best"] is None or r["loss"] < state["best"][1]["loss"]:
                state["best"] = (list(key), r)
                state["amp"], state["c"], state["la"] = r["amp"], r["c_band"], r["la"]      # warm start of the next inner fit
            return r["loss"]

        if kind == "pastorek":
            lo, hi = S["pastorek"]["zeta_bounds"]
            minimize_scalar(lambda z: outer([z]), bounds=(lo, hi), method="bounded", options={"xatol": 1e-3})
        else:
            lo, hi = np.log(S["schleiss"]["tau_bounds_min"])
            o = F_["outer_schleiss"]
            grid = [(tw, td) for tw in o["grid"]["tau_wet_min"] for td in o["grid"]["tau_dry_min"]]
            gl = {g: outer(list(map(float, g))) for g in grid}
            q0 = np.log(min(gl, key=gl.get))
            st = o["refine_step_log"]
            simplex = np.array([q0, q0 + [st, 0.0], q0 + [0.0, st]])
            simplex = np.where(simplex > hi, simplex - 2 * st, simplex)    # keep the first simplex inside the bounds
            minimize(lambda q: outer(np.exp(np.clip(q, lo, hi))), q0, method="Nelder-Mead",
                     options={"initial_simplex": simplex, "xatol": o["xatol"], "fatol": o["fatol"], "maxiter": o["maxiter"]})
        shape, r = state["best"]
        evals = len(done)
    out = {"name": name, "kind": kind, "version": version, "shape": [float(s) for s in shape], "amp": [float(v) for v in r["amp"]],
           "c_band": [float(v) for v in r["c_band"]], "a": dict(zip(T.ids, np.exp(r["la"]).tolist())), "loss": r["loss"],
           "inner_converged": r["converged"], "outer_evals": evals, "seconds": time.time() - t0}
    file.write_bytes(pickle.dumps(out))
    return out


# =================================================================================================== evaluation helpers
def verdict(lo, hi):
    return "ours better" if lo > 0 else ("ours worse" if hi < 0 else "no difference")


def fig_rmse(rm, BANDS):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    marks = {"ours": "o", "table F": "D", "B0": "x", "B1": "+", "B2": "1", "Schleiss a": "v", "Schleiss b": "^", "Pastorek a": "<", "Pastorek b": ">"}
    for ax, sub in zip(axes, ("all", "L < 1 km")):
        g = rm[(rm.set == "TEST") & (rm.subset == sub)].set_index("model").loc[ORDER].reset_index()
        for i, r in g.iterrows():
            ax.errorbar(i, r.rmse, yerr=[[r.rmse - r.lo], [r.hi - r.rmse]], fmt=marks[r.model], color=INK, mfc="white" if r.model != "ours" else INK,
                        capsize=3, ms=7)
        ours = float(g[g.model == "ours"].rmse.iloc[0])
        ax.axhline(ours, color=INK2, ls=":", lw=1)
        ax.text(len(g) - 0.5, ours, " ours", va="bottom", ha="right", fontsize=8)
        ax.set_xticks(range(len(g)), g.model, rotation=40, ha="right", fontsize=8)
        ax.set_title(f"TEST, {sub} (V1)", fontsize=9)
        style_axes(ax)
    axes[0].set_ylabel("RMSE of 5 min A (dB), 90 % event bootstrap")
    fig.suptitle("S8c (post-freeze): ours against literature WAA models. a: pycomlink defaults (Schleiss a: plus an assumed drying time = "
                 "wetting time, 15 min, not from the paper); b: fitted on TRAIN", fontsize=8)
    fig.tight_layout()
    return fig


def fig_dry(cv):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True)
    styles = (("measured", "-", 1.8), ("ours", "--", 1.4), ("Schleiss b", "-.", 1.2), ("Pastorek b", ":", 1.6))
    for ax, (_, nm) in zip(axes, LB):
        g = cv[cv.length == nm]
        for lab, ls, lw in styles:
            ax.plot(g.t_min, g[lab], ls, color=INK, lw=lw, label=lab)
        ax.set_title(f"{nm}: {int(g.windows.iloc[0]) if len(g) else 0} windows", fontsize=9)
        ax.set_xlabel("minutes after the path became dry")
        style_axes(ax)
    axes[0].set_ylabel("median W (dB)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("S8c: TEST dry-path W, measured against our equation and the fitted literature models (b)", fontsize=9)
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    clock = time.time()
    log = lambda s: print(s, flush=True)                                                    # noqa: E731
    fit_file = P.S6_CACHE / "plug_shared_FP_F_and_P.pkl"
    sha = hashlib.sha256(fit_file.read_bytes()).hexdigest()
    eq_sha = hashlib.sha256(P.S9_EQUATION.read_bytes()).hexdigest()
    if sha != S["ours"]["fit_sha256"]:
        log(f"STOP: the fit behind our equation changed (sha256 {sha} against {S['ours']['fit_sha256']}).")
        return
    eq = json.loads(P.S9_EQUATION.read_text(encoding="utf-8"))
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    BANDS = list(sid.merged_bands(d, DATA)[0])
    let = pickle.loads((P.SID3_CACHE / "le_train_all.pkl").read_bytes())[0].reset_index(drop=True)
    Dt1 = sid.load_data(P.S5A_CACHE / "data_train_v1.npz", BANDS)
    mods, _ = s6_models(Dt1, let, BANDS)
    fc, f0 = mods["closed F"]["fz"], mods["candidate"]["fz"]
    assert int(fc["d"]) == LAG == int(eq["lag_min"]), "the lag of our equation is not 1 min"
    log(f"TRAIN loaded ({(time.time() - clock) / 60:.1f} min)")

    # ------------------------------------------------------------------------------------------- the literature models on TRAIN
    T = Set(Dt1, let)
    c_ours = [eq["c_band_db_km_per_mmh"][b]["value"] for b in BANDS]
    fits = {}
    for kind in ("schleiss", "pastorek"):
        for v in ("a", "b"):
            f = fit_model(T, kind, v, c_ours, log)
            fits[f["name"]] = f
            log(f"{f['name']}: shape {np.round(f['shape'], 3).tolist()}, amp {np.round(f['amp'], 3).tolist()}, c_band {np.round(f['c_band'], 4).tolist()}, "
                f"loss {f['loss']:.1f} ({(time.time() - clock) / 60:.1f} min)")
    del T
    fits_out = {n: {k: v for k, v in f.items() if k != "a"} | {"a_link_p10_p90": [float(np.percentile(list(f["a"].values()), q)) for q in (10, 90)],
                                                             "a_link_per_link": {str(k): float(v) for k, v in f["a"].items()}} for n, f in fits.items()}
    save(json.dumps({"bands": BANDS, "shape_names": {"schleiss": ["tau_wet_min", "tau_dry_min"], "pastorek": ["zeta"]},
                     "amp_name": {"schleiss": "W_max_db per band", "pastorek": "C_db per band"}, "fits": fits_out}, indent=1), P.S8C_FITS)

    # ------------------------------------------------------------------------------------------- the sets, loaded one at a time (memory)
    g4 = pd.read_csv(P.COND4C_G)
    gval = dict(zip(g4[g4.set == "VAL"].event_id, g4[g4.set == "VAL"].g_path))
    le_tr = pickle.loads((P.S8_CACHE / "le_test_train.pkl").read_bytes())[0].reset_index(drop=True)
    le_te = pickle.loads((P.S8_CACHE / "le_test_test.pkl").read_bytes())[0].reset_index(drop=True)
    gt = pd.read_csv(P.S8_G_E)
    gm = {e: 1.0 for e in set(le_tr.event_id) | set(le_te.event_id)}           # TEST: V1 with the S8 g_e after its fix (as S8b)
    gm.update({r.event_id: (r.g_path if r.path_used else 1.0) for r in gt.itertuples()})

    def load(s_):
        if s_ == "TRAIN":
            return Dt1, let
        if s_ == VL:
            lev = pickle.loads((P.SID3_CACHE / "le_val_all.pkl").read_bytes())[0].reset_index(drop=True)
            return rescaled(sid.load_data(P.SID3_CACHE / "data_val_all.npz", BANDS), lev, gval), lev
        if s_ == "TEST":
            return rescaled(sid.load_data(P.S8_CACHE / "data_test_train.npz", BANDS), le_tr, gm), le_tr
        return rescaled(sid.load_data(P.S8_CACHE / "data_test_test.npz", BANDS), le_te, gm), le_te

    ours_models = {"ours": fc, "table F": f0, **b_models(Dt1, let, BANDS)}
    BR = S["evaluation"]["boot_reps"]
    rm_rows, pd_rows = [], []
    for s_ in ("TRAIN", VL, "TEST", "TEST links (a_link = 1)"):
        D, le = load(s_)
        pred = {m: predict(D, le, f, f["M"])[0] for m, f in ours_models.items()}
        Ts = Set(D, le)
        W1 = {}
        for n, f in fits.items():
            pred[n], W1[n] = Ts.predict(f)
        del Ts
        Lb = le.L_km.to_numpy()[D["bin_le"]]
        subs = [("all", np.ones(len(Lb), bool)), ("L < 1 km", Lb < 1)]
        if s_ == "TEST":
            subs += [(nm, (Lb >= lo) & (Lb < hi)) for (lo, hi), nm in LB if lo > 0]
        E = int(D["le_event"].max()) + 1
        Mlt = ev_draws(E, BR, S["evaluation"]["seed"])
        for sub, mask in subs:
            bt = {}
            for m in ORDER:
                s__, n__ = ev_sums(D, pred[m], mask)
                bt[m] = np.sqrt((Mlt @ s__) / np.maximum(Mlt @ n__, 1))
                lo, hi = ci(bt[m])
                rm_rows.append({"set": s_, "subset": sub, "model": m, "rmse": float(np.sqrt(s__.sum() / n__.sum())), "lo": lo, "hi": hi,
                                "bins": int(n__.sum())})
            if s_ in ("TEST", "TEST links (a_link = 1)"):
                r_o = rm_rows[-len(ORDER)]["rmse"]
                for k, m in enumerate(ORDER[1:], 1):
                    lo, hi = ci(bt[m] - bt["ours"])
                    pd_rows.append({"set": s_, "subset": sub, "model": m, "minus_ours_db": rm_rows[-len(ORDER) + k]["rmse"] - r_o, "lo": lo, "hi": hi,
                                    "verdict": verdict(lo, hi), "rule_applies": m in LIT and sub in ("all", "L < 1 km") and s_ == "TEST"})
        if s_ == "TEST":
            D1, W1test = D, W1                                                      # kept for the dry path
        del pred, W1, D
        log(f"{s_} scored ({(time.time() - clock) / 60:.1f} min)")
    rm, pdt = pd.DataFrame(rm_rows), pd.DataFrame(pd_rows)
    ro = rm[(rm.set == "TEST") & (rm.model == "ours")].set_index("subset").rmse
    s8b = pd.read_csv(P.S8B_TEST)
    s8b_ours = s8b[(s8b.part == "RMSE") & (s8b.model == "closed F")].set_index("subset").rmse
    repro = bool(all(abs(ro[k] - s8b_ours[k]) < 1e-3 for k in ("all", "L < 1 km")))
    log(f"RMSE done: ours on TEST {ro['all']:.4f} / {ro['L < 1 km']:.4f} dB (S8b {s8b_ours['all']:.4f} / {s8b_ours['L < 1 km']:.4f}; reproduced {repro})")

    # ------------------------------------------------------------------------------------------- the dry path on TEST
    W1t = {"ours": w_minutes(D1, le_tr, fc), "table F": w_minutes(D1, le_tr, f0), **W1test}
    Rm, T5, col = radar_max(d, {"TEST": (D1, le_tr)}, blank_test=False)
    win, curves = dry_windows({"TEST": (D1, le_tr)}, {"TEST": gm}, Rm, T5, col, {"TEST": W1t})
    obs = np.array([c["obs"] for c in curves])
    nmin = win.dry_minutes.to_numpy()
    draws = boot_items(pd.factorize(win.event_id)[0], BR, S["evaluation"]["seed"] + 1)
    dm = ["ours", "table F", *LIT]
    arr = {m: np.array([c[m] for c in curves]) for m in dm}
    dry_rows = []
    for t in S["evaluation"]["dry_t_min"]:
        for _, nm in LB + [((0, 100), "all")]:
            sel = (nmin > t) & np.isfinite(obs[:, t]) & ((win.length.to_numpy() == nm) if nm != "all" else True)
            r = {"t_min": t, "length": nm, "windows": int(sel.sum()), "events": int(win[sel].event_id.nunique()),
                 "measured_W": float(np.median(obs[sel, t])) if sel.any() else np.nan}
            for m in dm:
                err = obs[:, t] - arr[m][:, t]
                r[f"err {m}"] = float(np.median(err[sel])) if sel.any() else np.nan
                lo, hi = ci([np.median(err[i[sel[i]]]) for i in draws if sel[i].any()])
                r[f"err {m} lo"], r[f"err {m} hi"] = lo, hi
            dry_rows.append(r)
    dryt = pd.DataFrame(dry_rows)
    cv_rows = []
    for _, nm in LB:
        sel_c = win.length.to_numpy() == nm
        for t in range(60):
            ok = sel_c & (nmin > t)
            if ok.sum() < 3:
                continue
            row = {"length": nm, "windows": int(sel_c.sum()), "t_min": t, "n_at_t": int(ok.sum()), "measured": float(np.nanmedian(obs[ok, t]))}
            row.update({m: float(np.nanmedian(arr[m][ok, t])) for m in dm})
            cv_rows.append(row)
    cv = pd.DataFrame(cv_rows)
    log(f"dry path done: {len(win)} windows ({(time.time() - clock) / 60:.1f} min)")

    # ------------------------------------------------------------------------------------------- save, figures, report
    for tab, path in ((rm.round(4), P.S8C_RMSE), (pdt.round(4), P.S8C_PAIRED), (dryt.round(4), P.S8C_DRY), (cv.round(4), P.S8C_DRY_CURVES)):
        save(tab, path, index=False)
    for fig, path in ((fig_rmse(rm, BANDS), P.S8C_FIG_RMSE), (fig_dry(cv), P.S8C_FIG_DRY)):
        save(fig, path, dpi=100)
        plt.close(fig)
    report = write_report({"fits": fits, "rm": rm, "pdt": pdt, "dryt": dryt, "BANDS": BANDS, "sha": sha, "eq_sha": eq_sha, "repro": repro,
                           "s8b_ours": s8b_ours, "ro": ro, "nwin": len(win), "eq": eq})
    save(report, P.S8C_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    log(f"done in {(time.time() - clock) / 60:.1f} min")


def write_report(q):
    fits, rm, pdt, dryt, BANDS = q["fits"], q["rm"], q["pdt"], q["dryt"], q["BANDS"]
    L = []
    a = L.append
    a("# wet_antenna step S8c: our final WAA equation against two literature WAA models\n")
    a("Written by `research/wet_antenna/literature_s8c.py` (libraries `final_model.py`, `validation.py`, `system_id.py`, `test_s8b.py`). Every "
      "setting: `wet_antenna.literature_s8c` in `research/config.yaml`, written on 2026-10-28 BEFORE the run.\n")
    a("## 1. Status\n")
    a(STATUS + "\n")
    a(f"- Our equation: `{rel(P.S9_EQUATION)}` (sha256 `{q['eq_sha']}`), run forward as the S8b closed-F model (fit `"
      f"{rel(P.S6_CACHE / 'plug_shared_FP_F_and_P.pkl')}`, sha256 `{q['sha']}`, checked). Nothing of it is refitted.")
    a(f"- The forward run reproduces S8b on TEST: **{q['repro']}** (ours {q['ro']['all']:.4f} / {q['ro']['L < 1 km']:.4f} dB, S8b "
      f"{q['s8b_ours']['all']:.4f} / {q['s8b_ours']['L < 1 km']:.4f} dB, all / L < 1 km).")
    a("- The task is the S8 task: every model gets only V1 (radar / g_e, the S8 g_e after its fix), L and the band, and predicts A - c at 5 min. "
      "Same link-events, baseline, c and lag d = 1 min. The S6 table F and the S8 comparison models B0-B2 are shown as references.\n")
    a("## 2. The models as implemented\n")
    a("pycomlink version **0.6.0** (wheel from PyPI, `pycomlink/processing/wet_antenna.py`). The code used, verbatim (docstrings removed):\n")
    a("```python\n" + PYCOMLINK_PASTOREK + "\n```\n")
    a("```python\n" + PYCOMLINK_SCHLEISS + "\n```\n")
    a("**Pastorek et al. 2022 (KR-alt; preprint 2021).** Used as in pycomlink, with R = V1(t - d): W = C (1 - exp(-d R^zeta)), d = 0.1. No change. The same form is "
      "equation (5) / (A.9) of the paper (the open preprint on engrXiv), with d fixed to 0.1 by the authors. The pycomlink function "
      "`waa_pastorek_2021_from_A_obs`, which derives R from the link's own attenuation, is NOT used.\n")
    a("**Schleiss et al. 2013.** pycomlink's loop uses the measured attenuation A = rsl - baseline twice: as a cap in the wet step and as the "
      "whole value in the dry step (W = min(A, W_max), i.e. after the rain W simply follows the measurement). Both would hand the model the "
      "answer. What was changed:\n")
    a("1. the wet flag is the radar: V1(t - d) > 0.1 mm/h (not the link's own wet / dry classification);")
    a("2. the cap at A[i] is removed from the wet step (the step itself, W + (W_max - W) 3 dt / tau, and the cap at W_max are kept);")
    a("3. the dry step min(A[i], W_max) is replaced by an exponential decay W - W 3 dt / tau_dry (the paper's abstract: the WAA \"decreases "
      "exponentially as soon as the rain stops\"); without the measured A, pycomlink's dry step would hold W at W_max for ever;")
    a("4. W = 0 before each link-event window and the update starts at its first minute (pycomlink keeps waa[0] = 0 and updates from the "
      "second step; the windows start dry, so this does not matter); dt = 1 min.\n")
    a("```\n" + OURS_SCHLEISS + "\n```\n")
    a("`research/explorations/s8c_checks/check_pycomlink.py` checks the coded forms against pycomlink 0.6.0 itself (Pastorek: identical; "
      "Schleiss: identical on wet minutes with the A cap inactive, and different on dry minutes by design).\n")
    og = F_["outer_schleiss"]
    a("**Versions.** (a) the published / pycomlink parameters, a_link = 1, only c_band fitted on TRAIN. **Schleiss a = pycomlink defaults plus an "
      "assumed drying time equal to the wetting time (15 min), not taken from the paper.** Pastorek a = the pycomlink defaults C = 14 dB, "
      "zeta = 0.55, d = 0.1. (b) everything fitted on TRAIN: the amplitude per band (as our A_band), the shape shared over bands, c_band, and "
      "a_link with the ridge of our model (Huber 0.3 dB, a_prior_sd 0.5). **The verdict that counts is version (b).** For given shape "
      "parameters both models are linear in the amplitude, so the amplitudes, c_band and log a_link are fitted by L-BFGS-B, and the shape by "
      f"an outer search: Pastorek b a bounded scalar search over zeta, Schleiss b a grid (tau_wet {og['grid']['tau_wet_min']} min x tau_dry "
      f"{og['grid']['tau_dry_min']} min) and then Nelder-Mead from its best point (first steps x3 in each direction, maxiter {og['maxiter']}). "
      "The grid replaced the first plan (Nelder-Mead from the defaults) on 2026-10-28, before any result: the first run was stopped by low "
      "memory. Every outer step is cached at once.\n")
    rows = []
    for n, f in fits.items():
        sh = (f"zeta {f['shape'][0]:.3f}" if f["kind"] == "pastorek" else f"tau_wet {f['shape'][0]:.1f} min, tau_dry {f['shape'][1]:.1f} min")
        av = list(f["a"].values())
        rows.append({"model": n, "shape": sh, "amplitude per band (dB; C or W_max)": ", ".join(f"{b} {v:.2f}" for b, v in zip(BANDS, f["amp"])),
                     "c_band (dB/km per mm/h)": ", ".join(f"{v:.3f}" for v in f["c_band"]),
                     "a_link p10-p90": f"{np.percentile(av, 10):.2f}-{np.percentile(av, 90):.2f}", "TRAIN loss": round(f["loss"], 1),
                     "outer evals": f["outer_evals"], "inner converged": f["inner_converged"]})
    pt = pd.DataFrame(rows)
    a(md_table(pt) + "\n")
    eq, dy = q["eq"], q["eq"]["dynamics_fixed"]
    a(f"Our equation for comparison (not refitted): b = {eq['b']['value']:.3f} (u = a_link A_band (R / 5)^b), A_band "
      + " / ".join(f"{eq['A_band_db'][b]['value']:.2f}" for b in BANDS) + " dB, c_band "
      + " / ".join(f"{eq['c_band_db_km_per_mmh'][b]['value']:.3f}" for b in BANDS)
      + f", tau_f_wet {dy['tau_f_wet']:.1f}, tau_f_dry {dy['tau_f_dry']:.1f}, phi {dy['phi']:.2f}, tau_s {dy['tau_s']:.0f} min. Fitted parameters "
      "of every version: `" + rel(P.S8C_FITS) + "`.\n")
    a("**Parameters checked against the papers.**\n")
    a("- Pastorek: the form matches equation (5) / (A.9) of the preprint, d = 0.1 is the paper's choice. The paper's Table II lists NO literature "
      "values for KR-alt, so the pycomlink defaults C = 14 dB, zeta = 0.55 could not be checked against the paper text (they may be the "
      "calibrated values of its Fig. 5, which cannot be read from the text). The paper assumes the WAA of both antennas, as our W.")
    a("- Schleiss: W_max = 2.3 dB matches the paper's abstract (\"bounded by a maximum value of about 2.3 dB\"); W_max = 2.3 dB and tau = 15 min "
      "are the values Pastorek et al. (Table II) take from Schleiss et al. The full text of Schleiss et al. was not available, so the exact "
      "equations, the meaning of the factor 3 in pycomlink's step and the drying time constant could not be checked. tau_dry = tau_wet = 15 min "
      "in version (a) is an ASSUMPTION, so Schleiss a is labelled \"pycomlink defaults plus an assumed drying time\". The pycomlink example notebook uses W_max = 2.2 dB, tau = 15, delta_t = 1.")
    a("- Pastorek's form is W = C (1 - exp(-d R^zeta)). Its zeta acts inside a saturating function, so it is NOT comparable to our b (the "
      "exponent of a plain power law that sets the target of two stores with memory). No agreement between the exponents is claimed, "
      "whatever the numbers.\n")
    a("## 3. Results and verdicts\n")
    a("RMSE of the 5 min A (dB, pooled bins, 90 % event bootstrap; VAL was used for decisions in S3-S6b, so its numbers are optimistic):\n")
    w = rm[rm.subset.isin(["all", "L < 1 km"])].copy()
    w["value"] = [f"{r.rmse:.3f} [{r.lo:.3f}, {r.hi:.3f}]" for r in w.itertuples()]
    tab = w.pivot_table(index="model", columns=["set", "subset"], values="value", aggfunc="first").loc[ORDER]
    tab.columns = [f"{s} {u}" for s, u in tab.columns]
    a(md_table(tab.reset_index()) + "\n")
    a("TEST per length class:\n")
    w2 = rm[(rm.set == "TEST")].copy()
    w2["value"] = [f"{r.rmse:.3f} [{r.lo:.3f}, {r.hi:.3f}]" for r in w2.itertuples()]
    t2 = w2.pivot_table(index="model", columns="subset", values="value", aggfunc="first").loc[ORDER]
    a(md_table(t2.reset_index()) + "\n")
    a("![rmse](../figures/literature_s8c_rmse.png)\n")
    a("Paired differences on TEST (model RMSE minus ours, the same 1000 event draws; verdict by the rule written before the run for the "
      "literature models, all link-events and L < 1 km; the other rows are information):\n")
    pv = pdt[(pdt.set == "TEST")].copy()
    pv["difference"] = [f"{r.minus_ours_db:+.4f} [{r.lo:+.4f}, {r.hi:+.4f}]" for r in pv.itertuples()]
    a(md_table(pv[["subset", "model", "difference", "verdict", "rule_applies"]]) + "\n")
    vr = pv[pv.rule_applies]
    a("**Verdicts (as they came out).** " + "; ".join(f"{r.model}, {r.subset}: {r.verdict} ({r.minus_ours_db:+.3f} dB [{r.lo:+.3f}, {r.hi:+.3f}])"
                                                    for r in vr.itertuples()) + ".\n")
    pn = pdt[pdt.set == "TEST links (a_link = 1)"].copy()
    pn["difference"] = [f"{r.minus_ours_db:+.4f} [{r.lo:+.4f}, {r.hi:+.4f}]" for r in pn.itertuples()]
    a("Information: TEST links never fitted (a_link = 1 for every model):\n")
    a(md_table(pn[pn.subset.isin(["all", "L < 1 km"])][["subset", "model", "difference", "verdict"]]) + "\n")
    a(f"Dry-path error on TEST (W_obs - W_model, median over the S7 dry-path windows, {q['nwin']} windows, 90 % event bootstrap; positive = the "
      "model holds too little water):\n")
    dcols = ["t_min", "length", "windows", "measured_W"] + [f"err {m}" for m in ("ours", "table F", *LIT)]
    dd = dryt.copy()
    for m in ("ours", "table F", *LIT):
        dd[f"err {m}"] = [f"{v:+.3f} [{lo:+.3f}, {hi:+.3f}]" for v, lo, hi in zip(dryt[f"err {m}"], dryt[f"err {m} lo"], dryt[f"err {m} hi"])]
    a(md_table(dd[dcols].round(3)) + "\n")
    a("![dry path](../figures/literature_s8c_dry_path.png)\n")
    a("## 4. What changed in the documents\n")
    a("Written by hand after this run, from the tables of this step (`docs/WAA_research_overview.tex`, the one WAA document; the files "
      "`docs/overview/waa_overview.tex` and `waa_part1_full.tex` no longer exist, they were merged into it on 2026-10-19):\n")
    a("- Section 6 (the final test): a paragraph \"Comparison with published WAA models\" (what each model assumes, memory or not, the "
      "references; the Schleiss model driven by the radar wet flag with an exponential drying in place of the measured attenuation; added "
      "after the freeze; Schleiss a labelled as pycomlink defaults with an assumed drying time; Pastorek's zeta not comparable to b), "
      "Table 7 (TEST RMSE, the difference to ours and the verdicts), Figure 7 (`docs/WAA_research_overview_figures/literature.png`, made by "
      "`figs.py` from `literature_s8c_rmse.csv`) and a short reading.")
    a("- References: Schleiss et al. 2013 (IEEE GRSL 10(5), 1195-1199, DOI 10.1109/LGRS.2012.2236074) and Pastorek et al. 2022 (IEEE TGRS "
      "60, 1-9, DOI 10.1109/TGRS.2021.3110004), checked against Crossref; the DOI and title of Schleiss et al. also against the EPFL record.")
    a("- `docs/ROADMAP.md`: S8c in the step list and the table.\n")
    a("## 5. Unclear / possibly wrong\n")
    a("- **Schleiss a is pycomlink defaults plus an assumed drying time equal to the wetting time (15 min), not taken from the paper.** "
      "Its drying is ours, not pycomlink's: pycomlink has no drying law without the measured signal. The exponential drying step "
      "follows the paper's abstract, but its exact form and the published drying constant could not be checked. Version (b) fits tau_dry, so "
      "the fair comparison does not depend on that assumption, version (a) does.")
    a("- **The pycomlink defaults of KR-alt** (C = 14 dB, zeta = 0.55) are not in the paper's text, they could not be checked.")
    a("- **The literature models are used outside their setting.** Both papers fit the WAA so that the corrected link rain matches gauges, "
      "here the models predict A from radar rain. Version (a) uses parameters calibrated on other links (Schleiss: one 38 GHz link of 1.85 km, "
      "Pastorek: 16 links near Prague), so a poor (a) result says little about the papers.")
    a("- **The radar wet flag** (V1 > 0.1 mm/h) is a choice made here, not the paper's. The radar is coarse (5 min, 2 km), so short showers can "
      "be missed or smeared.")
    a("- **The comparison is post-freeze.** TEST had been read in S8 / S8b. The literature models were fitted on TRAIN only and nothing was "
      "tuned on TEST, but the outer search ran once and its result was not checked on fake data.")
    a("- **The parameters of the literature models have no intervals** (one TRAIN fit each, no bootstrap).")
    return "\n".join(L)


if __name__ == "__main__":
    main()
