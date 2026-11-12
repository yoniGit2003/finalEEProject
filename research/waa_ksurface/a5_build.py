"""Step A5 of the k-surface branch: the model is built and its code is checked. NO real attenuation is fitted here.

    python -m research.waa_ksurface.a5_build

The model (library `ksurf.py`): dW/dt = k(R, W) * (F_full(R) - W), A - c = W + L * P_full(R(t - d)), with the measured F and P
inside their valid ranges and a free monotone extension above. Three unit tests (pass rules in decisions.md, step A5):

    1. the exact gradient of the loss against finite differences (rung K3 in both regimes, extension free);
    2. a constant k reproduces the analytic exponential response to a step of rain;
    3. with its knots set to the values of the measured curve, the extension reduces to the measured curve.

The script exits with code 1 if a test fails.
"""
import sys

import numpy as np
import pandas as pd

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import data, fp, ksurf
from research.waa_ksurface import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.system_id import PADM

C = ksurf.config()
T = C["tests"]


def toy_dataset(rain_windows, L_km, y5=None):
    """A dataset in the data-layer layout from a list of 1 min rain series (one link-event each, a multiple of 5 minutes).
    The rain before and after each window (the padding) is zero."""
    n = np.array([len(r) for r in rain_windows])
    D = {"le_n": n, "le_boff1": np.r_[0, np.cumsum(n)[:-1]], "le_boff5": np.r_[0, np.cumsum(n // 5)[:-1]], "le_roff": np.r_[0, np.cumsum(n + 2 * PADM)[:-1]],
         "le_L": np.asarray(L_km, float), "le_c": np.zeros(len(n)), "le_event": np.arange(len(n)),
         "y5": np.zeros(int((n // 5).sum())) if y5 is None else y5, "w5": np.ones(int((n // 5).sum()))}
    rain = np.concatenate([np.r_[np.zeros(PADM), r, np.zeros(PADM)] for r in rain_windows])
    return D, rain


def test_gradient(M):
    """Exact gradient against central finite differences on real link-events with a fake target."""
    D, le, _ = data.gauge_ratio_train()
    rng = np.random.default_rng(T["seed"])
    peak = np.array([np.nanmax(D["r5"][o:o + n // 5]) for o, n in zip(D["le_boff5"], D["le_n"])])
    pick = np.zeros(len(le), bool)
    pick[rng.choice(np.nonzero(peak > C["f_knots_mmh"][0])[0], T["n_link_events"], replace=False)] = True       # link-events that reach the extension
    Ds = sid.subset(D, pick)
    les = le[pick].reset_index(drop=True)
    bands = list(M)
    band_idx = les.band_name.map({b: i for i, b in enumerate(bands)}).to_numpy()
    prob = ksurf.prepare(Ds, band_idx, [M[b] for b in bands], Ds["rain"], d=2)
    others = []                                                 # beyond the written test: the same check for the other rungs
    for rw, rd in (("K0", "K0"), ("K1", "K2"), ("K2", "K1")):
        others.append({"wetting": rw, "drying": rd, "relative_difference": _gradient_error(ksurf.Model(dict(prob), rw, rd, lam=1e-3), rng)[0]})
    m = ksurf.Model(prob, "K3", "K3", lam=1e-3)
    err, val, g, fd = _gradient_error(m, rng)
    names = np.empty(m.n, object)
    for k, s in m.sl.items():
        names[s] = k
    tab = pd.DataFrame({"block": names, "exact": g, "finite_difference": fd, "difference": g - fd})
    return {"link_events": int(pick.sum()), "minutes": int(Ds["le_n"].sum()), "parameters": m.n, "loss": val, "largest_gradient": float(np.abs(g).max()),
            "largest_abs_difference": float(np.abs(g - fd).max()), "relative_difference": err, "pass": err < T["grad_tol"], "others": pd.DataFrame(others)}, tab


def _gradient_error(m, rng):
    """(largest |exact - finite difference| relative to the largest gradient entry, loss, exact, finite difference)."""
    prob = m.p
    x_true = m.start() + rng.normal(0, 0.3, m.n) * (~np.isin(np.arange(m.n), np.r_[np.arange(m.n)[m.sl["F"]], np.arange(m.n)[m.sl["P"]]]))
    x_true[m.sl["F"]] = rng.uniform(0.1, 0.5, m.nb * m.nF)
    x_true[m.sl["P"]] = rng.uniform(0.1, 0.8, m.nb * m.nP)
    prob["target"] = m.simulate(x_true)[1] + rng.normal(0, 0.2, len(prob["target"]))        # a fake target: the model itself plus noise
    x = x_true + rng.normal(0, 0.1, m.n)
    x[m.sl["F"]], x[m.sl["P"]] = np.abs(x[m.sl["F"]]), np.abs(x[m.sl["P"]])
    val, g = m.loss(x)
    eps, fd = 1e-6, np.empty(m.n)
    for i in range(m.n):
        e = np.zeros(m.n)
        e[i] = eps
        fd[i] = (m.loss(x + e)[0] - m.loss(x - e)[0]) / (2 * eps)
    return float(np.abs(g - fd).max() / np.abs(g).max()), val, g, fd


def test_step(M):
    """Constant k: rain 0 for 30 min, 4 mm/h for 120 min, 0 for 180 min. W against the analytic exponentials."""
    band = list(M)[0]
    c = M[band]
    r0, tau_w, tau_d, t_on, t_len = 4.0, 5.0, 20.0, 30, 120                      # round numbers
    rain1 = np.r_[np.zeros(t_on), np.full(t_len, r0), np.zeros(180)]
    D, rain = toy_dataset([rain1], [0.5])
    prob = ksurf.prepare(D, np.array([0]), [c], rain, d=0)
    m = ksurf.Model(prob, "K0", "K0")
    x = m.start()
    x[m.sl["gw"]], x[m.sl["gd"]] = np.log(1 / tau_w), np.log(1 / tau_d)
    W = m.simulate(x)[0]                                                         # W after each minute
    F = float(ksurf.full_curve(c, "F", x[m.sl["F"]][:m.nF], [r0])[0])
    t = np.arange(1, len(rain1) + 1)
    exact = np.where(t <= t_on, 0.0, np.where(t <= t_on + t_len, F * (1 - np.exp(-(t - t_on) / tau_w)),
                                              F * (1 - np.exp(-t_len / tau_w)) * np.exp(-(t - t_on - t_len) / tau_d)))
    err = float(np.abs(W - exact).max())
    return {"rain_step_mmh": r0, "F_of_step_db": F, "tau_wet_min": tau_w, "tau_dry_min": tau_d, "minutes": len(rain1), "W_peak_db": float(W.max()),
            "largest_difference_db": err, "pass": err < T["step_tol"]}


def test_extension(M):
    """Knots set to the values of the measured E4 curve (F held, P with its last slope): F_full and P_full against that curve."""
    grid = np.r_[0.0, np.geomspace(0.01, 30.0, 2000)]
    rows = []
    for band, c in M.items():
        for kind, knots in (("F", C["f_knots_mmh"]), ("P", C["p_knots_mmh"])):
            curve = c.F if kind == "F" else c.P
            inc = ksurf.increments_of(c, kind, curve(np.asarray(knots, float)))
            full = ksurf.full_curve(c, kind, inc, grid)
            rows.append({"band": band, "quantity": kind, "knots_mmh": " ".join(f"{k:g}" for k in knots), "values_at_knots": " ".join(f"{v:.4f}" for v in curve(np.asarray(knots, float))),
                         "largest_difference": float(np.abs(full - curve(grid)).max())})
    tab = pd.DataFrame(rows)
    return {"largest_difference": float(tab.largest_difference.max()), "pass": bool(tab.largest_difference.max() < T["ext_tol"])}, tab


def main():
    M = fp.measured()
    g, g_tab = test_gradient(M)
    s = test_step(M)
    e, e_tab = test_extension(M)
    tests = pd.DataFrame([{"test": "1 gradient vs finite differences", "number": g["relative_difference"], "tolerance": T["grad_tol"], "pass": g["pass"]},
                          {"test": "2 constant-k step response (dB)", "number": s["largest_difference_db"], "tolerance": T["step_tol"], "pass": s["pass"]},
                          {"test": "3 extension = measured curve at measured knot values", "number": e["largest_difference"], "tolerance": T["ext_tol"], "pass": e["pass"]}])
    save(tests, P.A5_TESTS, index=False)
    save(g_tab, P.A5_GRADIENT, index=False)
    print(tests.to_string(), flush=True)
    ok = bool(tests["pass"].all())
    m0 = ksurf.Model({"nb": len(M), "curves": list(M.values())}, "K3", "K3")
    sizes = pd.DataFrame([{"wetting": w, "drying": d_, "parameters_k": ksurf.Model({"nb": 2, "curves": list(M.values())}, w, d_).sl["F"].start,
                           "parameters_extension": m0.nb * (m0.nF + m0.nP)} for w, d_ in (("K0", "K0"), ("K1", "K1"), ("K2", "K2"), ("K3", "K3"))])
    a = ["# A5: the model, built and checked (no real attenuation fitted)\n",
         "Written by `research/waa_ksurface/a5_build.py`; the model is in `research/waa_ksurface/ksurf.py`. Definitions: `docs/waa_ksurface/decisions.md`, step A5.\n",
         "## The model\n",
         "Per link-event, one step per minute, W = 0 at the start of the window:\n",
         "    v = F_full(R(t));  wetting if v > W, else drying;  k = exp(g(R(t)) + h(W));  W <- v + (W - v) exp(-k);  A - c = W + L P_full(R(t - d))\n",
         f"- F_full, P_full: the measured E4 curve inside its valid range; above it a monotone, continuous extension, linear between the knots "
         f"**F: {C['f_knots_mmh']} mm/h, P: {C['p_knots_mmh']} mm/h** (per band; parameters = increments >= 0; above 20 mm/h F is held, P keeps its last slope).",
         f"- k: ln k piecewise linear, R knots {C['r_knots_mmh']} mm/h (in ln(1 + R)), W knots {C['w_knots_wet_db']} dB (wetting) and {C['w_knots_dry_db']} dB (drying); "
         f"1 / k between {C['tau_min_min'] * 60:.0f} s and {C['tau_max_min'] / 60:.0f} h.",
         f"- Loss: mean Huber ({C['huber_delta_db']} dB) of the 5 min residuals over whole events + λ·(second differences)², λ from {C['lambdas']} by CV "
         f"({C['folds']} folds grouped by event). No per-link factor. Exact gradient, L-BFGS-B.\n",
         md_table(sizes),
         "\n## The rule for d (decisions.md, 2026-11-12)\n",
         "> d is a whole number of minutes from 0 to 5: the path term reads the rain of d minutes earlier (5 min is one radar bin; a longer lag cannot be placed "
         "with this input). d is chosen once per rain input (radar; gauge), with the K0 model (constant k_wet and k_dry, extension free) on the TRAIN learning set: "
         "the d with the smallest loss. If the paired 90 % event-bootstrap interval of loss(d = 0) − loss(d_best) contains zero, the lag is not resolved and d = 0. "
         "d is then held for every rung, fold and refit. It is never chosen per link or per event.\n",
         "## Unit tests\n", md_table(tests),
         f"\n1. **Gradient.** {g['link_events']} real link-events ({g['minutes']} minutes), a fake target, rung K3 in both regimes, extension free, d = 2: {g['parameters']} parameters. "
         f"Largest |exact − finite difference| = {g['largest_abs_difference']:.2e} against a largest gradient entry of {g['largest_gradient']:.2e}: relative "
         f"**{g['relative_difference']:.1e}** (tolerance {T['grad_tol']:.0e}). Every parameter: `{rel(P.A5_GRADIENT)}`. Beyond the written test, the same check for the "
         "other rungs (wetting / drying): " + "; ".join(f"{r.wetting} / {r.drying}: {r.relative_difference:.1e}" for r in g["others"].itertuples()) + ".",
         f"2. **Step response.** Rain 0 → {s['rain_step_mmh']:g} mm/h for 120 min → 0; 1 / k = {s['tau_wet_min']:g} min wetting, {s['tau_dry_min']:g} min drying; F = {s['F_of_step_db']:.3f} dB. "
         f"Largest |W − analytic| over {s['minutes']} minutes = **{s['largest_difference_db']:.1e} dB** (tolerance {T['step_tol']:.0e}).",
         f"3. **Extension.** Knots set to the values of the measured curve: largest |F_full − F| and |P_full − P| on 2 000 values of R up to 30 mm/h = "
         f"**{e['largest_difference']:.1e}** (tolerance {T['ext_tol']:.0e}).\n", md_table(e_tab),
         f"\n**Verdict: {'PASS' if ok else 'FAIL'}.** (`{rel(P.A5_TESTS)}`)\n"]
    save("\n".join(a) + "\n", P.A5_REPORT)
    print(f"report: {rel(P.A5_REPORT)}", flush=True)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
