"""WAA_estm step 4: the drying law on the REAL tails (frozen tail set v1, TRAIN tails only).

The first step that fits real data. Only tools that passed on fake data are used, with their settings unchanged
(read from the config sections of the fake-data steps; `waa.real4` holds what is new here):
  (a) menu fits F1, F3, F4, F5, F6 with the pooled likelihood of waa_pooled.py (offset prior, ARMA(1,1) noise, speed
      spread marginalised), 4 starts each (3 generic + the F1 solution; F1: 4 generic), BIC on the whole dataset;
  (b) the free curve of step 3b-3b: dW/dt = -f_i W k(W), 5 knots from 0.35 dB to W_hi, read ONLY for W >= 0.35 dB,
      no step tool; band = event-bootstrap refits (each from 2 starts);
  (c) jump 2 on that curve, h(W) = W k(W), library W, W^2, W^3, max(W - 0.5, 0), max(W - 0.75, 0): E-SINDy (bagged
      backward elimination of waa_sindy.py; main), PySR (cross-check), best-subset + BIC (list);
  (d) literature: F1's tau against drying time constants quoted in the project documents (none are invented).
Checks: 1 the four sensitivity sets; 2 leaving out each of the 3 big events; 3 menu fits on the first 6 h;
4 residuals vs the dry-weather noise and a gallery of 12 tails; 5 the real numbers against the fake truths.

The 123 VAL tails are never read (every read filters split_event == 'train'); nothing in outputs/products changes.
One fit is timed first and the projection decides the cuts (order: check 3, check 2 to the biggest event, bootstraps
20 -> 10). Every fit is cached in .cache/fits4/, so a restart resumes; delete that folder to redo.
Run from the repository root (about 1 h the first time on 3 processes, under a minute when cached):
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step4_real_drying.py
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step4_real_drying.py --pysr   # only the PySR part (run by the line above)
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "1")            # one thread per process: the processes run in parallel

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import json
import subprocess
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.signal import lfilter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MODULE = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE))
import WAA_estm_paths as P
from research.common.save import save
import waa_pooled as wp
import waa_curve as wc
import waa_sindy as sd
import step3b1c_stress_test as st                          # load_noise, f1_start, shape
import step3b3_curve_then_equation as s3                   # load_tails, start_theta, md, rel

CFG = s3.CFG
W = CFG["waa"]
S, S1C, FIT, S3B, A3C, D3C = W["real4"], W["stress3b1c"], W["curve3b3"], W["curve3b3b"], W["curve3b3c"]["a"], W["curve3b3c"]["d"]
W_LO = S3B["w_lo_db"]
JTERMS = wc.curve_terms(S3B["j2_switch_grid_db"])          # W, W^2, W^3, sw0.5, sw0.75
WIN, WIN_LONG = S["window_min"], S["window_long_min"]
FAMILY = {"F1": "exponential", "F3": "nonlinear in W (F3 / F5)", "F5": "nonlinear in W (F3 / F5)",
          "F4": "two layers", "F6": "state switch"}
INK, INK2 = "#1a1a1a", "#52514e"
CACHE = P.CACHE / "fits4"
md, rel = s3.md, s3.rel


# =================================================================================================== data
class _Noise(dict):
    """Noise parameters per sub-link; a sub-link without control tails (alternative set only) gets the medians."""
    def __missing__(self, key):
        return self.default


def load_noise():
    nz, noise, offset_sd = st.load_noise()
    out = _Noise(noise)
    out.default = tuple(float(nz[c].median()) for c in ("sigma_white_db", "sigma_slow_db", "tau_noise_min"))
    return nz, out, offset_sd


def data_file(name):
    return CACHE / f"data_{name}.parquet"


_INPUTS = {}


def pack_set(name, window):
    """Tails of one dataset and its packed arrays for the first `window` minutes; N of the BIC = (link, event) pairs."""
    if not _INPUTS:
        _, _INPUTS["noise"], _INPUTS["offset_sd"] = load_noise()
    tails = s3.load_tails(data_file(name))
    D = wp.pack(tails, _INPUTS["noise"], _INPUTS["offset_sd"], window, FIT["min_points"])
    used = set(D["tail_id"])
    n_ind = len({(tl["cml_id"], tl["event_copy"]) for tl in tails if tl["tail_id"] in used})
    return tails, D, n_ind


def _save(f_, res):
    tmp = f_.with_suffix(".tmp")
    tmp.write_text(json.dumps(res))
    tmp.replace(f_)


# =================================================================================================== workers
def menu_starts(cand):
    return list(S1C["t2_starts"][cand]) + ([] if cand == "F1" else ["f1sol"])


def menu_file(name, window, cand, label):
    return CACHE / f"menu_{name}_{window}_{cand}_{label}.json"


def run_menu(job):
    """(a) one candidate on one dataset and window, from every start (generic starts of step 3b-1c + the start built
    from the F1 solution). Each fit is one cached json file."""
    name, window, cand = job["name"], job["window"], job["cand"]
    data, out = {}, []
    for label in menu_starts(cand):
        f_ = menu_file(name, window, cand, label)
        if f_.exists():
            out.append(json.loads(f_.read_text()))
            continue
        if not data:
            _, data["D"], data["n"] = pack_set(name, window)
        if label == "f1sol":                               # the best F1 fit of this dataset (finished in phase 1)
            f1 = max((json.loads(menu_file(name, window, "F1", lb).read_text()) for lb in menu_starts("F1")),
                     key=lambda r: r["ll"])
            theta0, s0 = st.f1_start(cand, f1["theta"][0]), f1["s"]
        else:
            theta0, s0 = S1C["t2_starts"][cand][label]
        t0 = time.time()
        r = wp.fit(cand, data["D"], theta0, s0, S1C["maxfev_generic"])
        r.update({"dataset": name, "window": window, "candidate": cand, "start": label,
                  "theta0": [float(v) for v in theta0], "s0": float(s0), "n_tails": len(data["D"]["tail_id"]),
                  "n_independent": data["n"], "n_points": int(len(data["D"]["T"])), "k": len(wp.NAMES[cand]) + 1,
                  "seconds": round(time.time() - t0, 1)})
        _save(f_, r)
        out.append(r)
    return out


def run_curve(job):
    """(b) the free-curve pooled fit of step 3b-3b on one dataset (first 60 min): 4 starts, the best continued, the
    other penalties from that solution, penalty chosen by 2-fold cross-validation over EVENTS. Stages are cached."""
    f_ = CACHE / f"curve_{job['name']}.json"
    res = json.loads(f_.read_text()) if f_.exists() else {}
    if res.get("done"):
        return res
    t0 = time.time()
    tails, D, n_ind = pack_set(job["name"], WIN)
    cv = wc.Curve(job["knots"], prior_sd=FIT["prior_sd_logk"])
    lam0, mf = FIT["lambda_start"], FIT["maxfev"]
    res.update({"dataset": job["name"], "knots": job["knots"], "w_hi": job["w_hi"], "n_tails": len(D["tail_id"]),
                "n_independent": n_ind, "n_points": int(len(D["T"]))})
    res.setdefault("starts", {})
    for name, spec in FIT["starts"].items():
        if name not in res["starts"]:
            res["starts"][name] = cv.fit(D, np.r_[s3.start_theta(spec, cv.knots), np.log(FIT["s_start"])], lam0, mf["start"])
            _save(f_, res)
    res.setdefault("lam", {})
    if str(lam0) not in res["lam"]:
        best = min(res["starts"], key=lambda n: res["starts"][n]["cost"])
        res["lam"][str(lam0)] = {**cv.fit(D, res["starts"][best]["x"], lam0, mf["finish"]), "from": best}
        _save(f_, res)
    for lam in FIT["lambda_grid"]:
        if str(lam) not in res["lam"]:
            res["lam"][str(lam)] = cv.fit(D, res["lam"][str(lam0)]["x"], lam, mf["other_lambda"])
            _save(f_, res)
    res.setdefault("cv", {})
    evs = np.random.default_rng(S["cv_seed"]).permutation(sorted({tl["event_copy"] for tl in tails}))
    half, Dh = [set(evs[::2].tolist()), set(evs[1::2].tolist())], None
    for lam in FIT["lambda_grid"]:
        if str(lam) not in res["cv"]:
            if Dh is None:
                Dh = [wp.pack([tl for tl in tails if tl["event_copy"] in h], _INPUTS["noise"], _INPUTS["offset_sd"], WIN,
                              FIT["min_points"]) for h in half]
            res["cv"][str(lam)] = float(sum(cv.loglik(np.asarray(cv.fit(Dh[a], res["lam"][str(lam)]["x"], lam, mf["cv"])["x"]),
                                                      Dh[b]) for a, b in ((0, 1), (1, 0))))
            _save(f_, res)
    lam_star = max(FIT["lambda_grid"], key=lambda v: res["cv"][str(v)])
    fit = res["lam"][str(lam_star)]
    res.update({"lam_star": lam_star, "x": fit["x"], "ll": fit["ll"], "s": float(np.exp(fit["x"][-1])),
                "bic": -2 * fit["ll"] + (cv.K + 1) * np.log(n_ind), "seconds": round(res.get("seconds", 0.0) + time.time() - t0, 1),
                "done": True})
    _save(f_, res)
    return res


def run_boot(job):
    """Band: the free curve refitted on ONE event bootstrap of the main set (events drawn with replacement, same knots
    and penalty), from 2 starts: the fit of the dataset and the flat start. The lower cost is kept."""
    f_ = CACHE / f"boot_{job['b']}.json"
    if f_.exists():
        return json.loads(f_.read_text())
    t0 = time.time()
    tails, _, _ = pack_set("main", WIN)
    by_ev = {}
    for tl in tails:
        by_ev.setdefault(tl["event_copy"], []).append(tl)
    evs = sorted(by_ev)
    rng = np.random.default_rng([S["bootstrap_seed"], job["b"]])
    drawn = [int(e) for e in rng.choice(evs, len(evs), replace=True)]
    Db = wp.pack([tl for e in drawn for tl in by_ev[e]], _INPUTS["noise"], _INPUTS["offset_sd"], WIN, FIT["min_points"])
    cv = wc.Curve(job["knots"], prior_sd=FIT["prior_sd_logk"])
    x0 = {"fit": np.asarray(job["x"]),
          "flat": np.r_[s3.start_theta(FIT["starts"]["flat"], cv.knots), np.log(FIT["s_start"])]}
    fits = {k: cv.fit(Db, x0[k], job["lam"], D3C["maxfev"]) for k in S["bootstrap_starts"]}
    best = min(fits, key=lambda k: fits[k]["cost"])
    res = {"b": job["b"], "events": drawn, "x": fits[best]["x"], "from": best,
           "cost": {k: v["cost"] for k, v in fits.items()}, "seconds": round(time.time() - t0, 1)}
    _save(f_, res)
    return res


def run_any(job):
    return {"menu": run_menu, "curve": run_curve, "boot": run_boot}[job["part"]](job)


# =================================================================================================== PySR (own process)
def run_pysr(names):
    """(c) cross-check: PySR on h(W) = W k(W) of the fitted curve of each dataset, settings of step 3b-3c part A.
    Stores PySR's own 'best' formula and its whole Pareto front (complexity, loss, formula)."""
    from pysr import PySRRegressor
    Y = A3C["pysr"]
    f_ = CACHE / "pysr.json"
    done = json.loads(f_.read_text()) if f_.exists() else {}
    for name in names:
        if name in done:
            continue
        r = json.loads((CACHE / f"curve_{name}.json").read_text())
        w = np.geomspace(W_LO, r["w_hi"], FIT["grid_points"])
        h = w * wc.Curve(r["knots"]).k(np.asarray(r["x"])[:-1], w)
        t0 = time.time()
        model = PySRRegressor(niterations=Y["niterations"], maxsize=Y["maxsize"], timeout_in_seconds=Y["timeout_s"],
                              binary_operators=["+", "-", "*", "/", "^"], unary_operators=["exp"],
                              constraints={"^": (-1, 1)}, parallelism="serial", progress=False, verbosity=0,
                              output_directory=str(P.CACHE / "pysr"))
        for attempt in range(3):                           # Windows sometimes locks PySR's checkpoint file: try again
            try:
                model.fit(w[:, None], h, weights=1 / h ** 2, variable_names=["W"])
                break
            except PermissionError:
                if attempt == 2:
                    raise
                time.sleep(2)
        eq = model.equations_
        done[name] = {"best": str(model.sympy()), "seconds": round(time.time() - t0, 1),
                      "front": [{"complexity": int(c), "loss": float(l), "equation": str(e)}
                                for c, l, e in zip(eq.complexity, eq.loss, eq.equation)]}
        _save(f_, done)
        print(f"   PySR {name}: {done[name]['best']} ({done[name]['seconds']} s)", flush=True)


# =================================================================================================== reading helpers
def menu_k(cand, theta, w):
    """k(W) = -(dW/dt) / W of a menu law at speed factor 1; None for F4 (two states: not a function of W)."""
    w = np.asarray(w, float)
    if cand == "F1":
        return np.full_like(w, 1 / theta[0])
    if cand == "F3":
        return theta[0] * w ** (theta[1] - 1)
    if cand == "F5":
        return theta[0] + theta[1] * w
    if cand == "F6":
        return np.where(w > theta[2], theta[1] + theta[0] * (1 - theta[2] / w), theta[1])
    return None


def read_menu(best):
    """The 3b-1c reading rules on the best fit of every candidate of one dataset (rows: candidate, bic).
    F3 and F5 are one family; F4 and F6 closer than f4_f6_bic are not separated. Returns winner, family, margins."""
    b = best.set_index("candidate").bic.sort_values()
    win = b.index[0]
    fam, members = FAMILY[win], {c for c in b.index if FAMILY[c] == FAMILY[win]}
    if win in ("F4", "F6"):
        other = "F6" if win == "F4" else "F4"
        if b[other] - b[win] < S["f4_f6_bic"]:
            fam, members = "two layers or state switch (F4 / F6 not separated)", {"F4", "F6"}
    rest = b[[c for c in b.index if c not in members]]
    return {"winner": win, "family": fam, "runner_up": b.index[1], "margin": float(b.iloc[1] - b.iloc[0]),
            "next_family": rest.index[0], "margin_family": float(rest.iloc[0] - b.iloc[0]),
            "f1_minus_winner": float(b["F1"] - b.iloc[0])}


def read_slope(slope):
    R = S["reading"]
    if abs(slope) < R["flat_abs_slope"]:
        return "exponential in this range"
    if slope > R["quadratic_slope"]:
        return "quadratic family in this range"
    if slope > 0:
        return "mildly rising, law not identified by the curve"
    return "falling (no fake truth looked like this)"


def formula_family(terms):
    t = set(terms)
    if t == {"W"}:
        return "exponential"
    if t in ({"W^2"}, {"W", "W^2"}):
        return "quadratic family"
    return "multi-term (not a named law)"


def formula_text(coef):
    """{'W': 0.03, 'sw0.5': 0.01} -> '0.03 W + 0.01 max(W - 0.5, 0)'."""
    name = lambda t: f"max(W - {t[2:]}, 0)" if t.startswith("sw") else t   # noqa: E731
    return " + ".join(f"{v:.4g} {name(t)}" for t, v in coef.items()).replace("+ -", "- ") or "(nothing)"


def agree_class(text):
    """One of three classes for the agreement test: exponential / quadratic / other."""
    return "exponential" if text.startswith("exponential") else ("quadratic" if text.startswith(("quadratic", "nonlinear")) else "other")


def pysr_misfit(equation, w, h):
    """Relative rms misfit of a PySR formula for h(W) (its own `loss` is weighted and not comparable)."""
    import sympy as sp
    try:
        Wsym = sp.Symbol("W", positive=True)
        val = np.broadcast_to(np.asarray(sp.lambdify(Wsym, sp.sympify(equation.replace("^", "**"), locals={"W": Wsym}), "numpy")(w),
                                         float), w.shape)
        return float(np.sqrt(np.mean((val / h - 1) ** 2)))
    except Exception:                                      # noqa: BLE001  (a formula sympy cannot read)
        return np.nan


def raw_rate(d, edges, smooth):
    """Model-free look (information): mean drying rate per level of A. A is smoothed per tail (centred mean of
    `smooth` minutes), dA/dt = central difference; per bin of smoothed A: k = -mean(dA/dt) / median(A)."""
    d = d.sort_values(["tail_id", "minutes_since_t0_trace"])
    sm = d.groupby("tail_id").A_db.transform(lambda a: a.rolling(smooth, center=True, min_periods=smooth - 2).mean())
    dd = sm.groupby(d.tail_id).transform(lambda v: (v.shift(-1) - v.shift(1)) / 2)
    x = pd.DataFrame({"tail_id": d.tail_id, "A": sm, "dadt": dd}).dropna()
    g = x.groupby(pd.cut(x.A, edges), observed=True)
    out = g.agg(tail_minutes=("dadt", "size"), tails=("tail_id", "nunique"), median_A_db=("A", "median"),
                mean_dadt_db_per_min=("dadt", "mean"), median_dadt_db_per_min=("dadt", "median")).reset_index(names="A_bin_db")
    out["k_from_mean_per_min"] = -out.mean_dadt_db_per_min / out.median_A_db
    out["A_bin_db"] = out.A_bin_db.astype(str)
    return out.round(4)


def jump2(w, k, rng):
    """The three jump-2 tools on one curve (PySR is added from its cache by the caller)."""
    h = w * k
    E = A3C["esindy"]
    inc, med = wc.bagged_select(w, h, JTERMS, E["threshold"], E["bootstraps"], rng)
    kept = [t for t, p_ in zip(JTERMS, inc) if p_ >= E["inclusion"]]
    # the bagged medians are taken term by term and need not fit together when many terms are kept (real curve: all
    # five), so the printed coefficients are one least-squares fit (relative residuals) on the kept terms
    G = wc.term_matrix(w, kept) if kept else np.zeros((len(w), 0))
    es = dict(zip(kept, map(float, np.linalg.lstsq(G / h[:, None], np.ones(len(w)), rcond=None)[0]))) if kept else {}
    subs = wc.best_subset(w, h, JTERMS, A3C["max_terms"])
    by_size = {n: next(s_ for s_ in subs if s_["size"] == n) for n in range(1, A3C["max_terms"] + 1)}
    return {"esindy": es, "esindy_inclusion": dict(zip(JTERMS, map(float, inc))),
            "esindy_bagged_median": {t: float(c) for t, c in zip(JTERMS, med) if t in kept},
            "esindy_rms": float(np.sqrt(np.mean((G @ np.array(list(es.values())) / h - 1) ** 2))) if es else np.nan,
            "subset_best": subs[0], "subset_by_size": by_size}


# =================================================================================================== check 4 helpers
def tail_fits(cand, theta, s, D):
    """Per tail of D: the fitted curve of a menu law at the tail's speed factor (posterior mode, as
    waa_pooled.tail_speed_factors), with its profiled start W0 and offset c. Returns a list of dicts in D's order."""
    code, lin, th = wp.CAND_CODE[cand], cand in wp.LINEAR, np.asarray(theta, float)
    _, gstar, kap = wp._tails_laplace(code, lin, th, float(s), D["T"], D["FY"], D["F1"], D["start"], D["length"], D["rho"],
                                      D["ma"], D["s2"], D["g0"], D["o2"], D["a0"], wp.G_LATTICE)
    f = np.exp(s * kap * gstar / (1 + kap))
    out = []
    for i, (a, n) in enumerate(zip(D["start"], D["length"])):
        t = D["T"][a:a + n]
        _, w0, c = wp._profile(code, lin, th, t, D["FY"][a:a + n], D["F1"][a:a + n], D["rho"][i], D["ma"][i], D["s2"][i],
                               D["g0"][i], D["o2"][i], D["a0"][i], f[i])
        w = np.empty(n)
        wp._curve(code, th, t, w0, f[i], w)
        out.append({"tail_id": D["tail_id"][i], "t": t, "y": D["Y"][a:a + n], "model": w + c, "f": float(f[i]),
                    "w0": float(w0), "c": float(c)})
    return out


def series_stats(series, lags):
    """series: list of (t minutes, values). Median within-tail std, and the pooled autocorrelation at each lag (each
    tail's own mean removed; only pairs of minutes exactly `lag` apart)."""
    std = float(np.median([np.std(v) for _, v in series]))
    num, cnt, var, nv = {k: 0.0 for k in lags}, {k: 0 for k in lags}, 0.0, 0
    for t, v in series:
        x = np.full(int(t.max()) + 1, np.nan)
        x[t.astype(int)] = v - np.mean(v)
        var, nv = var + float(np.nansum(x * x)), nv + int(np.isfinite(x).sum())
        for k in lags:
            if k < len(x):
                p = x[:-k] * x[k:]
                num[k], cnt[k] = num[k] + float(np.nansum(p)), cnt[k] + int(np.isfinite(p).sum())
    return std, {k: (num[k] / cnt[k]) / (var / nv) if cnt[k] > 50 else np.nan for k in lags}


def quality_rows(label, window, fits, D, sub_of, noise, rng):
    """Residual statistics of one fitted law in one window, next to the same statistics of pure dry-weather noise
    drawn from the sub-link's noise model (same tails, same minutes)."""
    lags = S["acf_lags_min"]
    res = [(r["t"], r["y"] - r["model"]) for r in fits]
    std, acf = series_stats(res, lags)
    q = n = 0.0                                            # whitened residuals: 1 if they are exactly the dry noise
    for i, (t, r) in enumerate(res):
        eta = np.empty(len(r))
        wp._whiten(r, D["rho"][i], D["ma"][i], eta)
        q += eta[0] ** 2 / D["g0"][i] + float(np.sum(eta[1:] ** 2)) / D["s2"][i]
        n += len(r)
    sim = []
    for _ in range(S["noise_draws"]):
        ser = []
        for r in fits:
            sw, ss, tau = noise[sub_of[r["tail_id"]]]
            rho = np.exp(-1 / tau)
            m = len(r["t"])
            ar = lfilter([1.0], [1.0, -rho], np.r_[rng.normal(0, ss), rng.normal(0, ss * np.sqrt(1 - rho ** 2), m - 1)])
            ser.append((r["t"], ar + rng.normal(0, np.sqrt((sw ** 2 + W["sim"]["rsl_step_db"] ** 2 / 12) / W["sim"]["samples_per_min"]), m)))
        sim.append(series_stats(ser, lags))
    row = {"law": label, "window_min": window, "tails": len(fits), "within_tail_std_db": round(std, 4),
           "noise_model_std_db": round(float(np.median([s_[0] for s_ in sim])), 4),
           "whitened_variance_ratio": round(q / n, 3)}
    for k in lags:
        row[f"acf_{k}"] = round(acf[k], 3) if np.isfinite(acf[k]) else np.nan
        v = [s_[1][k] for s_ in sim]
        row[f"noise_model_acf_{k}"] = round(float(np.nanmedian(v)), 3) if np.isfinite(v).any() else np.nan
    return row


# =================================================================================================== figures
def _axes_style(ax):
    ax.grid(alpha=0.25, which="both")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def k_figure(cur, band, boot_n, f1, win, slope):
    """The real k(W): free curve (thick solid), event-bootstrap band (hatched), F1 (dashed) and the menu winner
    (dash-dotted) on top. Text labels and line styles, no colour coding."""
    fig, ax = plt.subplots(figsize=(8.6, 5.6))
    w = band.w_db.to_numpy()
    wf = np.geomspace(0.1, w[-1], 160)
    cv = wc.Curve(cur["knots"])
    kf = cv.k(np.asarray(cur["x"])[:-1], wf)
    ax.axvspan(w[0], w[-1], color="#ececec", lw=0, label=f"range that is read: {w[0]:.2f} to {w[-1]:.2f} dB (shaded)")
    ax.fill_between(w, band.k_lo, band.k_hi, facecolor="none", edgecolor="#6f6f6f", hatch="///", lw=0,
                    label=f"band of {boot_n} event-bootstrap refits, 5-95 % (hatched)")
    lo = wf < w[0]
    ax.plot(wf[lo], kf[lo], ":", color=INK2, lw=1.4)
    ax.plot(wf[~lo], kf[~lo], "-", color=INK, lw=2.6, label="free curve k(W) of the real tails (thick solid)")
    ax.text(wf[-1] * 1.03, kf[-1], "free curve", fontsize=9, va="center")
    ax.axhline(1 / f1["theta"][0], ls="--", color=INK, lw=1.6, label=f"F1 exponential, tau = {f1['theta'][0]:.1f} min (dashed, flat)")
    ax.text(wf[0] * 1.05, 1.04 / f1["theta"][0], "F1 (exponential)", fontsize=9, va="bottom")
    if win["candidate"] == "F4":                           # two states: the rate depends on where the tail started
        tf, ts, phi = win["theta"]
        for a0, ls, nm in ((1.0, "-.", "dash-dotted"), (3.0, (0, (5, 1, 1, 1, 1, 1)), "dash-dot-dotted")):
            t = np.linspace(0, WIN, 2000)
            wt = a0 * ((1 - phi) * np.exp(-t / tf) + phi * np.exp(-t / ts))
            ax.plot(wt[1:], -np.diff(wt) / np.diff(t) / wt[1:], ls=ls, color=INK, lw=1.6,
                    label=f"menu winner F4 (two layers): rate along a tail that starts at {a0:g} dB ({nm})")
            ax.text(wt[1] * 1.03, (-np.diff(wt) / np.diff(t) / wt[1:])[0], f"F4 from {a0:g} dB", fontsize=8.5, va="bottom")
    elif win["candidate"] != "F1":
        kw = menu_k(win["candidate"], win["theta"], wf)
        ax.plot(wf, kw, "-.", color=INK, lw=1.6, label=f"menu winner {win['candidate']} ({FAMILY[win['candidate']]}) (dash-dotted)")
        ax.text(wf[len(wf) // 3], kw[len(wf) // 3] * 0.93, win["candidate"], fontsize=9, va="top")
    ax.axvline(w[0], color=INK, lw=1.0)
    ax.set(xscale="log", yscale="log", xlim=(wf[0], wf[-1] * 1.45), xlabel="W (dB of water on the antennas)",
           ylabel="k(W) = drying rate per dB of water (1 / min)")
    ax.text(w[0] * 0.96, ax.get_ylim()[0] * 1.15, "below 0.35 dB: not read (dotted)  ", fontsize=8.5, ha="right", va="bottom")
    _axes_style(ax)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.set_title(f"Real drying tails (TRAIN, v1), first {WIN} min: slope of log k vs log W = {slope:+.2f} over the read range\n"
                 f"On fake data this band was about 3 times too narrow: read the level as +-{S['quote']['level_pct']} %, "
                 f"the slope as +-{S['quote']['slope']}", fontsize=9.5)
    fig.tight_layout()
    return fig


def sets_figure(curves, labels):
    """k(W) of every dataset on its own read range: one line style per group, a text label at the right end."""
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    styles = {"main": ("-", 2.6), "no_": ("--", 1.5), "alt": ("-.", 1.5), "loo": (":", 1.7)}
    for name, g in curves.groupby("dataset", sort=False):
        ls, lw = next(v for k, v in styles.items() if name.startswith(k))
        ax.plot(g.w_db, g.k, ls=ls, color=INK, lw=lw)
        ax.text(g.w_db.iloc[-1] * 1.02, g.k.iloc[-1], labels[name], fontsize=7.5, va="center")
    ax.set(xscale="log", yscale="log", xlabel="W (dB), read range of each set", ylabel="k(W) (1 / min)",
           xlim=(W_LO * 0.9, curves.w_db.max() * 2.3))
    _axes_style(ax)
    ax.set_title("Step 4 checks 1 and 2: free curve of every set (main thick solid; sensitivity dashed; alternative\n"
                 "dash-dotted; leave-one-event-out dotted)", fontsize=9.5)
    fig.tight_layout()
    return fig


def gallery_figure(picks, fw, f1, win_name):
    fig, axes = plt.subplots(3, 4, figsize=(14, 8.2), sharex=True)
    for ax, tid in zip(axes.ravel(), picks):
        a, b = fw[tid], f1[tid]
        ax.plot(a["t"], a["y"], "o", ms=2.6, color="#8a8a8a", label="measured A, 1 min (dots)")
        ax.plot(a["t"], a["model"], "-", color=INK, lw=2.0, label=f"fitted winner {win_name} (solid)")
        ax.plot(b["t"], b["model"], "--", color=INK, lw=1.3, label="fitted F1 exponential (dashed)")
        ax.axhline(0, color=INK2, lw=0.6)
        ax.set_title(f"tail {tid}: speed factor {a['f']:.2f}, offset {a['c']:+.2f} dB", fontsize=8.5)
        _axes_style(ax)
    axes[0, 0].legend(frameon=False, fontsize=7.5)
    for ax in axes[-1]:
        ax.set_xlabel("minutes since t0_trace")
    for ax in axes[:, 0]:
        ax.set_ylabel("A (dB)")
    fig.suptitle(f"Step 4 check 4: {len(picks)} random real tails (seed {S['gallery_seed']}), first {WIN} min, with the pooled fits",
                 fontsize=10)
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    if "--pysr" in sys.argv:
        return run_pysr(sys.argv[sys.argv.index("--pysr") + 1:])
    clock = time.time()
    workers = max((os.cpu_count() or 2) - 1, 1)

    # ------------------------------------------------------------------------------------------- datasets (TRAIN only)
    cols = ["tail_id", "sublink_id", "cml_id", "event_id", "minutes_since_t0_trace", "A_db"]
    ix = pd.read_csv(P.TAILS_MODEL_INDEX).query("split_event == 'train'").sort_values("tail_id").reset_index(drop=True)
    alt_ix = pd.read_csv(P.TAILS_MODEL_ALT_INDEX).query("split_event == 'train'").sort_values("tail_id").reset_index(drop=True)
    rows = pq.read_table(P.TAILS_MODEL, filters=[("split_event", "=", "train")], columns=cols).to_pandas()
    alt_rows = pq.read_table(P.TAILS_MODEL_ALT, filters=[("split_event", "=", "train")], columns=cols).to_pandas()
    assert set(rows.tail_id) == set(ix.tail_id) and set(alt_rows.tail_id) == set(alt_ix.tail_id)   # no VAL tail
    ev_code = {e: i for i, e in enumerate(sorted(set(ix.event_id) | set(alt_ix.event_id)))}
    big = list(ix.event_id.value_counts().index[:S["big_events"]])
    sets = {"main": ("v1 main", ix, rows),
            "no_rise_heavy": ("without event_rise_heavy", ix[~ix.event_rise_heavy], rows),
            "no_shift_undershoot": ("without level_shift / undershoot", ix[~ix.level_shift & ~ix.undershoot], rows),
            "alt": ("alternative set (A_start >= 0.6 dB)", alt_ix, alt_rows)}
    for e in big:
        sets[f"loo_{e}"] = (f"without event {e}", ix[ix.event_id != e], rows)
    SENS = ["main", "no_rise_heavy", "no_shift_undershoot", "alt"]
    meta = {}
    for name, (label, sel, src) in sets.items():
        d = src[src.tail_id.isin(sel.tail_id) & (src.minutes_since_t0_trace < WIN_LONG)].merge(
            sel[["tail_id", "A_start_trace_db"]], on="tail_id")
        d["event_copy"] = d.event_id.map(ev_code)
        if not data_file(name).exists():
            d.drop(columns="event_id").to_parquet(data_file(name))
        a = d.A_db[(d.minutes_since_t0_trace < WIN) & d.A_db.notna()].to_numpy()
        w_hi = float(np.sort(a)[-FIT["support_min_minutes"]])           # the 200 tail-minutes rule of step 3b-3
        inside = np.log(a[(a >= W_LO) & (a <= w_hi)])
        q = np.quantile(inside, np.linspace(0, 1, S3B["knots"])[1:-1])
        meta[name] = {"label": label, "tails": len(sel), "events": int(sel.event_id.nunique()), "links": int(sel.cml_id.nunique()),
                      "sublinks": int(sel.sublink_id.nunique()), "w_hi": w_hi, "points": int(len(a)),
                      "minutes_in_range": int(len(inside)), "knots": [W_LO] + [float(v) for v in np.exp(q)] + [w_hi],
                      "big3_tails": int(sel.event_id.isin(big).sum())}
    M = meta["main"]
    nz, noise, _ = load_noise()
    n_default = len(set(alt_ix.sublink_id) - set(nz.sublink_id))
    print(f"TRAIN tails only: main {M['tails']} tails, {M['events']} events, {M['links']} links; W_hi = {M['w_hi']:.2f} dB; "
          f"alternative set {meta['alt']['tails']} tails ({n_default} sub-links with median noise parameters); {workers} processes")

    # ------------------------------------------------------------------------------------------- time one fit, plan
    plan_file = CACHE / "plan.json"
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
        print(f"plan from {rel(plan_file)}: cut: {plan['cuts'] or 'nothing'}")
    else:
        _, D, _ = pack_set("main", WIN)
        cv = wc.Curve(M["knots"], prior_sd=FIT["prior_sd_logk"])
        x_flat = np.r_[s3.start_theta(FIT["starts"]["flat"], cv.knots), 0.0]
        wp.pooled_loglik("F1", [30.0], 1.0, D), cv.loglik(x_flat, D)        # compile
        t0 = time.time()
        for _ in range(5):
            wp.pooled_loglik("F5", [0.03, 0.01], 1.0, D)
        t_menu = (time.time() - t0) / 5
        t0 = time.time()
        for _ in range(5):
            cv.loglik(x_flat, D)
        t_curve = (time.time() - t0) / 5
        t0 = time.time()
        one = run_menu({"name": "main", "window": WIN, "cand": "F1"})       # the 4 generic F1 fits of the main set
        t_fit = (time.time() - t0) / len(one)
        nf = float(np.mean([r["nfev"] for r in one]))
        print(f"timing on the main set: one likelihood evaluation {t_menu * 1e3:.0f} ms (menu law), {t_curve * 1e3:.0f} ms "
              f"(free curve); one F1 fit {t_fit:.1f} s ({nf:.0f} evaluations)", flush=True)
        ev_curve = (len(FIT["starts"]) * FIT["maxfev"]["start"] + FIT["maxfev"]["finish"]
                    + (len(FIT["lambda_grid"]) - 1) * FIT["maxfev"]["other_lambda"] + len(FIT["lambda_grid"]) * FIT["maxfev"]["cv"])
        n_menu = sum(len(menu_starts(c)) for c in S["candidates"])
        pts = {n: m["points"] / M["points"] for n, m in meta.items()}

        def project(loo, long_window, boots):
            sec = sum(pts[n] * (n_menu * S["typical_menu_nfev"] * t_menu + ev_curve * t_curve) for n in SENS + loo)
            sec += long_window * n_menu * S["typical_menu_nfev"] * t_menu * WIN_LONG / WIN
            return (sec + boots * len(S["bootstrap_starts"]) * D3C["maxfev"] * t_curve) / workers

        plan = {"loo": [f"loo_{e}" for e in big], "long_window": True, "bootstraps": S["bootstraps"], "cuts": [],
                "t_menu_ms": round(t_menu * 1e3, 1), "t_curve_ms": round(t_curve * 1e3, 1), "t_f1_fit_s": round(t_fit, 1),
                "wall_s": 0.0}
        proj = project(plan["loo"], True, plan["bootstraps"])
        plan["projected_full_s"] = round(proj)
        print(f"projection: {proj / 3600:.2f} h on {workers} processes (budget {S['time_budget_h']} h)", flush=True)
        for cut, change in (("check 3 (6 h window)", {"long_window": False}),
                            ("check 2 reduced to the biggest event", {"loo": [f"loo_{big[0]}"]}),
                            (f"bootstraps {S['bootstraps']} -> {S['bootstraps_cut']}", {"bootstraps": S["bootstraps_cut"]})):
            if proj <= S["time_budget_h"] * 3600:
                break
            plan.update(change)
            plan["cuts"].append(cut)
            proj = project(plan["loo"], plan["long_window"], plan["bootstraps"])
            print(f"   cut: {cut} -> projected {proj / 3600:.2f} h", flush=True)
        plan["projected_s"] = round(proj)
        plan_file.write_text(json.dumps(plan))
    names = SENS + plan["loo"]
    cells = [(n, WIN) for n in names] + ([("main", WIN_LONG)] if plan["long_window"] else [])

    def run(jobs, label):
        t_run, out = time.time(), []
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(run_any, j): j for j in jobs}
            for k, fu in enumerate(as_completed(futs)):
                out.append((futs[fu], fu.result()))
                j = futs[fu]
                print(f"   {label} {k + 1}/{len(jobs)}: {j['part']} {j.get('name', j.get('b'))} {j.get('cand', '')} "
                      f"({time.time() - t_run:.0f} s)", flush=True)
        return out, time.time() - t_run

    size = lambda j: -meta[j.get("name", "main")]["points"] * (j.get("window", WIN) / WIN)   # noqa: E731  long jobs first
    # phase 1: F1 of every cell (its solution seeds the other candidates) and the free curves
    jobs1 = [{"part": "menu", "name": n, "window": w_, "cand": "F1"} for n, w_ in cells] + [
        {"part": "curve", "name": n, "knots": meta[n]["knots"], "w_hi": meta[n]["w_hi"]} for n in names]
    out1, wall1 = run(sorted(jobs1, key=size), "phase 1")
    curve = {r["dataset"]: r for j, r in out1 if j["part"] == "curve"}
    # phase 2: the other candidates, and the bootstrap refits of the main curve
    jobs2 = [{"part": "menu", "name": n, "window": w_, "cand": c} for n, w_ in cells for c in S["candidates"] if c != "F1"] + [
        {"part": "boot", "b": b, "knots": M["knots"], "x": curve["main"]["x"], "lam": curve["main"]["lam_star"]}
        for b in range(plan["bootstraps"])]
    out2, wall2 = run(sorted(jobs2, key=size), "phase 2")
    boots = sorted((r for j, r in out2 if j["part"] == "boot"), key=lambda r: r["b"])
    fits = pd.DataFrame([r for j, rs in out1 + out2 if j["part"] == "menu" for r in rs])
    fresh = sum(r["seconds"] for r in boots) + fits.seconds.sum() + sum(c["seconds"] for c in curve.values())
    if wall1 + wall2 > 60:                                 # a run that only read the cache does not change the run time
        plan["wall_s"] = round(plan["wall_s"] + wall1 + wall2, 1)
        plan_file.write_text(json.dumps(plan))

    # ------------------------------------------------------------------------------------------- (c) PySR, own process
    pysr_names = [n for n in SENS if n in curve]
    pysr_file = CACHE / "pysr.json"
    have = json.loads(pysr_file.read_text()) if pysr_file.exists() else {}
    if any(n not in have for n in pysr_names):
        print("PySR (separate process) ...", flush=True)
        try:
            subprocess.run([sys.executable, str(Path(__file__).resolve()), "--pysr", *pysr_names], check=True,
                           timeout=60 * 30)
        except Exception as e:                             # noqa: BLE001  (PySR / Julia missing: the cross-check is left out)
            print(f"   PySR failed: {type(e).__name__}: {e}")
        have = json.loads(pysr_file.read_text()) if pysr_file.exists() else {}

    # ------------------------------------------------------------------------------------------- (a) menu tables
    fits["bic"] = -2 * fits.ll + fits.k * np.log(fits.n_independent)
    menu_rows, reads = [], {}
    for (name, window), g in fits.groupby(["dataset", "window"], sort=False):
        best = g.loc[g.groupby("candidate").ll.idxmax()].copy()
        best["dbic"] = best.bic - best.bic.min()
        reads[(name, window)] = read_menu(best)
        for r in best.sort_values("bic").itertuples():
            gs = g[g.candidate == r.candidate]
            menu_rows.append({"dataset": name, "window_min": window, "candidate": r.candidate, "family": FAMILY[r.candidate],
                              "n_tails": r.n_tails, "n_independent": r.n_independent, "k": r.k, "ll": round(r.ll, 2),
                              "bic": round(r.bic, 2), "dbic": round(r.dbic, 2),
                              "parameters": ", ".join(f"{n_} {v:.4g}" for n_, v in zip(wp.NAMES[r.candidate], r.theta)),
                              "s": round(r.s, 3), "best_start": r.start, "starts": len(gs),
                              "bic_spread_starts": round(float(2 * (gs.ll.max() - gs.ll.min())), 1),
                              "theta": json.dumps([round(float(v), 6) for v in r.theta])})
    menu = pd.DataFrame(menu_rows)
    # Like-for-like reading (added after the first run, see the report). F1 and F4 are linear in the start W0, and
    # waa_pooled profiles it without a bound: a tail that RISES gets a negative W0 (negative water). F3, F5, F6 and the
    # free curve keep W0 >= 0.001 dB. "F1+" = the exponential with W0 >= 0: the F3 code at n = 1 with F1's fitted tau
    # and s (one likelihood evaluation; F3's own fit is used instead if it ended at n = 1 with a higher likelihood).
    plus_rows, neg_start = [], {}
    for (name, window), g in menu.groupby(["dataset", "window_min"], sort=False):
        f1_, f3_ = g[g.candidate == "F1"].iloc[0], g[g.candidate == "F3"].iloc[0]
        _, Dn, _ = pack_set(name, window)
        tau = json.loads(f1_.theta)[0]
        ll = wp.pooled_loglik("F3", [1 / tau, 1.0], f1_.s, Dn)
        if json.loads(f3_.theta)[1] < 1.001:
            ll = max(ll, f3_.ll)
        w0 = np.array([r["w0"] for r in tail_fits("F1", [tau], f1_.s, Dn)])
        neg_start[(name, window)] = int((w0 < 0).sum())
        plus_rows.append({**f1_.to_dict(), "candidate": "F1+", "family": "exponential", "ll": round(ll, 2),
                          "bic": round(-2 * ll + 2 * np.log(f1_.n_independent), 2), "best_start": "F1 solution", "starts": 1,
                          "bic_spread_starts": np.nan, "dbic": np.nan, "parameters": f"tau {tau:.4g} (start W0 >= 0)"})
    menu = pd.concat([menu, pd.DataFrame(plus_rows)], ignore_index=True)
    menu["dbic_like_for_like"] = np.nan
    fair = menu.candidate.isin(["F1+", "F3", "F5", "F6"])
    menu.loc[fair, "dbic_like_for_like"] = (menu[fair].bic - menu[fair].groupby(["dataset", "window_min"]).bic.transform("min")).round(2)
    fair_read = {}
    for (name, window), g in menu[fair].groupby(["dataset", "window_min"], sort=False):
        b = g.set_index("candidate").bic.sort_values()
        fam = lambda c: "exponential" if c == "F1+" else FAMILY[c]   # noqa: E731
        rest = b[[c for c in b.index if fam(c) != fam(b.index[0])]]
        fair_read[(name, window)] = {"winner": b.index[0], "family": fam(b.index[0]), "next_family": rest.index[0],
                                     "margin_family": float(rest.iloc[0] - b.iloc[0]),
                                     "f1_minus_winner": float(b["F1+"] - b.iloc[0]),
                                     "f3_at_bound": bool(json.loads(g[g.candidate == "F3"].theta.iloc[0])[1] < 1.001),
                                     "f5_b": float(json.loads(g[g.candidate == "F5"].theta.iloc[0])[1])}
    save(menu, P.STEP4_MENU, index=False)
    best_of = lambda name, window, cand: menu[(menu.dataset == name) & (menu.window_min == window) & (menu.candidate == cand)].iloc[0]   # noqa: E731
    as_fit = lambda r: {"candidate": r.candidate, "theta": json.loads(r.theta), "s": float(r.s)}   # noqa: E731

    # ------------------------------------------------------------------------------------------- (b) curves, band, slopes
    w_common = np.geomspace(W_LO, min(meta[n]["w_hi"] for n in names), FIT["grid_points"])
    cur_rows, slopes, j2, wj, hj = [], {}, {}, {}, {}
    for name in names:
        r = curve[name]
        cv = wc.Curve(r["knots"])
        w = np.geomspace(W_LO, r["w_hi"], FIT["grid_points"])
        k = cv.k(np.asarray(r["x"])[:-1], w)
        slopes[name] = (float(np.log(k) @ wc.slope_weights(np.log(w))),
                        float(np.log(cv.k(np.asarray(r["x"])[:-1], w_common)) @ wc.slope_weights(np.log(w_common))))
        j2[name] = jump2(w, k, np.random.default_rng([S["bootstrap_seed"], names.index(name)]))
        wj[name], hj[name] = w, w * k
        cur_rows += [{"dataset": name, "w_db": w[i], "k": k[i]} for i in range(len(w))]
    cur = pd.DataFrame(cur_rows)
    wm = cur[cur.dataset == "main"].w_db.to_numpy()
    aw = wc.slope_weights(np.log(wm))
    cvm = wc.Curve(M["knots"])
    KB = np.array([cvm.k(np.asarray(b["x"])[:-1], wm) for b in boots])
    b_slopes = np.log(KB) @ aw
    lo, hi = np.percentile(KB, D3C["band"], axis=0)
    for col, v in (("k_lo", lo), ("k_hi", hi), ("k_min", KB.min(axis=0)), ("k_max", KB.max(axis=0))):
        cur.loc[cur.dataset == "main", col] = v
    save(cur.round(6), P.STEP4_CURVES,
                index=False)
    band_half = float(np.median(np.sqrt(hi / lo) - 1))
    slope, slope_c = slopes["main"]
    k_main = cur[cur.dataset == "main"].k.to_numpy()
    # jump 2 on every bootstrap curve (stability of the E-SINDy formula; information)
    boot_formulas = pd.Series([" + ".join(t for t, c in zip(JTERMS, sd.select(
        (wc.term_matrix(wm, JTERMS) / (wm * kb)[:, None]) / np.linalg.norm(wc.term_matrix(wm, JTERMS) / (wm * kb)[:, None], axis=0),
        np.ones(len(wm)), A3C["esindy"]["threshold"])) if c != 0) or "(nothing)" for kb in KB]).value_counts()

    # ------------------------------------------------------------------------------------------- (c) jump 2 table
    j2_rows = []
    for name in names:
        J = j2[name]
        j2_rows.append({"dataset": name, "tool": "E-SINDy (waa_sindy elimination, bagged)", "formula": formula_text(J["esindy"]),
                        "terms": " + ".join(J["esindy"]), "family": formula_family(J["esindy"]),
                        "relative_rms_misfit": round(J["esindy_rms"], 4),
                        "note": "inclusion: " + ", ".join(f"{t} {p_:.2f}" for t, p_ in J["esindy_inclusion"].items())
                                + "; coefficients = least squares on the kept terms"})
        bic1 = J["subset_by_size"][1]["bic"]
        for n_, sub in J["subset_by_size"].items():
            j2_rows.append({"dataset": name, "tool": f"best-subset, {n_} term{'s' if n_ > 1 else ''}",
                            "formula": formula_text(dict(zip(sub["terms"], sub["coef"]))), "terms": " + ".join(sub["terms"]),
                            "family": formula_family(sub["terms"]),
                            "relative_rms_misfit": round(float(np.sqrt(sub["rss"] / len(wm))), 4),
                            "note": f"BIC {sub['bic'] - bic1:+.0f} vs the 1-term formula"
                                    + ("; BIC pick" if sub["terms"] == J["subset_best"]["terms"] else "")})
        if name in have:
            front = sorted(have[name]["front"], key=lambda e: e["complexity"])
            for i, e in enumerate(front[:3]):
                j2_rows.append({"dataset": name, "tool": f"PySR, simplest {i + 1} of its Pareto front", "formula": e["equation"],
                                "terms": "", "family": "", "relative_rms_misfit": round(pysr_misfit(e["equation"], wj[name], hj[name]), 4),
                                "note": f"complexity {e['complexity']}"})
            j2_rows.append({"dataset": name, "tool": "PySR, its own 'best' pick", "formula": have[name]["best"], "terms": "",
                            "family": "", "relative_rms_misfit": round(pysr_misfit(have[name]["best"], wj[name], hj[name]), 4),
                            "note": f"Pareto front of {len(front)} formulas"})
    jt = pd.DataFrame(j2_rows)
    save(jt, P.STEP4_JUMP2, index=False)

    # ------------------------------------------------------------------------------------------- summary per dataset
    sum_rows = []
    for name in names:
        rd, c_ = reads[(name, WIN)], curve[name]
        cr, jf = read_slope(slopes[name][0]), formula_family(j2[name]["esindy"])
        classes = [agree_class(rd["family"]), agree_class(cr), agree_class(jf)]
        if cr.startswith("falling"):                       # no fake truth gave a falling curve: never counted as agreement
            classes[1] = "falling"
        sum_rows.append({"dataset": name, "what": meta[name]["label"], "tails": c_["n_tails"], "events": meta[name]["events"],
                         "links": meta[name]["links"], "tails_big3": meta[name]["big3_tails"],
                         "w_hi_db": round(meta[name]["w_hi"], 2), "minutes_in_range": meta[name]["minutes_in_range"],
                         "menu_winner": rd["winner"], "menu_family": rd["family"], "runner_up": rd["runner_up"],
                         "margin_bic": round(rd["margin"], 1), "next_family": rd["next_family"],
                         "margin_family_bic": round(rd["margin_family"], 1), "f1_minus_winner_bic": round(rd["f1_minus_winner"], 1),
                         "f1_negative_start_tails": neg_start[(name, WIN)],
                         "like_for_like_winner": fair_read[(name, WIN)]["winner"],
                         "like_for_like_margin_family_bic": round(fair_read[(name, WIN)]["margin_family"], 1),
                         "like_for_like_f1_minus_winner_bic": round(fair_read[(name, WIN)]["f1_minus_winner"], 1),
                         "f3_n_at_lower_bound": fair_read[(name, WIN)]["f3_at_bound"],
                         "f1_tau_min": round(json.loads(best_of(name, WIN, "F1").theta)[0], 1),
                         "s_winner": best_of(name, WIN, rd["winner"]).s, "curve_slope": round(slopes[name][0], 3),
                         "curve_slope_common_range": round(slopes[name][1], 3), "curve_reading": cr,
                         "curve_s": round(c_["s"], 3), "curve_penalty": c_["lam_star"],
                         "curve_minus_winner_bic": round(c_["bic"] - best_of(name, WIN, rd["winner"]).bic, 1),
                         "jump2_esindy": formula_text(j2[name]["esindy"]), "jump2_family": jf,
                         "agree": len(set(classes)) == 1})
    summ = pd.DataFrame(sum_rows)
    save(summ, P.STEP4_SUMMARY, index=False)
    R0, S0 = reads[("main", WIN)], summ.iloc[0]
    winner = best_of("main", WIN, R0["winner"])
    f1 = best_of("main", WIN, "F1")

    # ------------------------------------------------------------------------------------------- check 4: residuals, gallery
    tails, D, _ = pack_set("main", WIN)
    sub_of = {tl["tail_id"]: tl["sublink_id"] for tl in tails}
    rng = np.random.default_rng(S["noise_seed"])
    fit_w, fit_1 = tail_fits(winner.candidate, json.loads(winner.theta), winner.s, D), tail_fits("F1", json.loads(f1.theta), f1.s, D)
    q_rows = [quality_rows(f"{winner.candidate} (winner)", WIN, fit_w, D, sub_of, noise, rng)]
    if winner.candidate != "F1":
        q_rows.append(quality_rows("F1", WIN, fit_1, D, sub_of, noise, rng))
    if plan["long_window"]:
        _, DL, _ = pack_set("main", WIN_LONG)
        for c in dict.fromkeys([winner.candidate, "F1"]):
            lab = c + (" (winner of the first hour)" if c == winner.candidate else "")
            b = best_of("main", WIN_LONG, c)
            q_rows.append(quality_rows(lab, WIN_LONG, tail_fits(c, json.loads(b.theta), b.s, DL), DL, sub_of, noise, rng))
    _, Dq, _ = pack_set("no_rise_heavy", WIN)               # the same statistics on the set without the rise-heavy events
    for c in dict.fromkeys([reads[("no_rise_heavy", WIN)]["winner"], "F1"]):
        b = best_of("no_rise_heavy", WIN, c)
        q_rows.append(quality_rows(f"{c}, set without event_rise_heavy", WIN, tail_fits(c, json.loads(b.theta), b.s, Dq), Dq,
                                   sub_of, noise, rng))
    used = nz[nz.sublink_id.isin(ix.sublink_id)]
    q_rows.append({"law": "real control tails (dry weather, median over the sub-links of the main set)", "window_min": int(used.median_control_len_min.median()),
                   "tails": int(used.n_control_tails.sum()), "within_tail_std_db": round(float(used.within_std_db.median()), 4),
                   **{f"acf_{k}": round(float(used[f"acf_{k}"].median()), 3) for k in S["acf_lags_min"]}})
    quality = pd.DataFrame(q_rows)
    save(quality, P.STEP4_QUALITY, index=False)
    picks = sorted(np.random.default_rng(S["gallery_seed"]).choice(D["tail_id"], S["gallery_n"], replace=False).tolist())
    fig = gallery_figure(picks, {r["tail_id"]: r for r in fit_w}, {r["tail_id"]: r for r in fit_1}, winner.candidate)
    save(fig, P.STEP4_FIG_GALLERY, dpi=110)
    plt.close(fig)

    # ------------------------------------------------------------------------------------------- check 5: vs the fake truths
    mins = S["shape_minutes"]
    main_rows = rows[rows.tail_id.isin(ix.tail_id)]
    a0 = ix.set_index("tail_id").A_start_trace_db
    per_tail = pd.DataFrame({k: main_rows[(main_rows.minutes_since_t0_trace >= k) & (main_rows.minutes_since_t0_trace < k + 5)]
                             .groupby("tail_id").A_db.mean() / a0 for k in mins})
    per_tail["event"] = ix.set_index("tail_id").event_id
    real_shape = per_tail[mins].median()
    per_tail["rise_heavy"] = ix.set_index("tail_id").event_rise_heavy
    stay = pd.DataFrame([{"tails": lab, "n": len(g_), **{f"share_A_ge_{int(th * 100)}pct_of_start_at_{k}min": round(float((g_[k] >= th).mean()), 3)
                                                       for k in mins[:2] for th in (0.9, 1.0)}}
                         for lab, g_ in (("all", per_tail), ("events flagged event_rise_heavy", per_tail[per_tail.rise_heavy]),
                                         ("other events", per_tail[~per_tail.rise_heavy]))])
    raw = raw_rate(main_rows[main_rows.minutes_since_t0_trace < WIN], S["raw_rate"]["edges_db"], S["raw_rate"]["smooth_min"])
    rb = np.random.default_rng(S["noise_seed"] + 1)
    groups = [g[mins] for _, g in per_tail.groupby("event")]
    bs = np.array([pd.concat([groups[i] for i in rb.integers(0, len(groups), len(groups))]).median().to_numpy()
                   for _ in range(S["ci_reps"])])
    shape_ci = np.nanpercentile(bs, [5, 95], axis=0)
    f_rows = [{"source": "REAL tails (TRAIN v1)", "slope": round(slope, 3),
               "slope_range": f"{np.percentile(b_slopes, 5):.2f}..{np.percentile(b_slopes, 95):.2f} (bootstrap 5-95 %)",
               **{f"ratio_{k}min": round(float(real_shape[k]), 3) for k in mins},
               **{f"range_{k}min": f"{shape_ci[0, i]:.3f}..{shape_ci[1, i]:.3f}" for i, k in enumerate(mins)}}]
    cs3 = pd.read_csv(P.STEP3B3B_CURVE_SHAPE).query("kind == 'main'").set_index("truth")
    sim_dir = P.CACHE / "fits3b3b"
    for truth in cs3.index:
        row = {"source": f"fake truth {truth} (3b-3b datasets, 5 replicates)", "slope": cs3.loc[truth, "slope_median"],
               "slope_range": f"{cs3.loc[truth, 'slope_min']:.2f}..{cs3.loc[truth, 'slope_max']:.2f}",
               "slope_inside": bool(cs3.loc[truth, "slope_min"] - S["quote"]["slope"] <= slope <= cs3.loc[truth, "slope_max"] + S["quote"]["slope"])}
        files = [sim_dir / f"sim_main_{truth}_{r}.parquet" for r in range(S3B["replicates"])]
        if all(f_.exists() for f_ in files):
            sh = []
            for f_ in files:
                d = pd.read_parquet(f_, columns=["tail_id", "minutes_since_t0_trace", "A_db", "A_start_trace_db"])
                sh.append(st.shape(d, d.groupby("tail_id").A_start_trace_db.first().to_dict(), mins))
            for k in mins:
                v = [s_[k] for s_ in sh]
                row[f"ratio_{k}min"] = round(float(np.median(v)), 3)
                row[f"range_{k}min"] = f"{min(v):.3f}..{max(v):.3f}"
                row[f"real_inside_{k}min"] = bool(min(v) <= real_shape[k] <= max(v))
        f_rows.append(row)
    rec = pd.read_csv(P.STEP3B1C_RECAL).query("test == 'T3'")
    for r in rec.itertuples():
        f_rows.append({"source": f"fake truth {r.truth} recalibrated to the real early decay (3b-1c T3)",
                       **{f"ratio_{k}min": getattr(r, f"new_sim_{k}min") for k in mins},
                       **{f"range_{k}min": f"{getattr(r, f'new_sim_min_{k}min'):.3f}..{getattr(r, f'new_sim_max_{k}min'):.3f}" for k in mins},
                       **{f"real_inside_{k}min": bool(getattr(r, f"new_sim_min_{k}min") <= real_shape[k] <= getattr(r, f"new_sim_max_{k}min"))
                          for k in mins}})
    fake = pd.DataFrame(f_rows)
    save(fake, P.STEP4_FAKE, index=False)

    # ------------------------------------------------------------------------------------------- figures
    fig = k_figure(curve["main"], cur[cur.dataset == "main"], len(boots), as_fit(f1), as_fit(winner), slope)
    save(fig, P.STEP4_FIG_K, dpi=120)
    plt.close(fig)
    fig = sets_figure(cur, {n: f"{n} ({slopes[n][0]:+.2f})" for n in names})
    save(fig, P.STEP4_FIG_SETS, dpi=110)
    plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    def menu_md(name, window):
        g = menu[(menu.dataset == name) & (menu.window_min == window)]
        return md(g[["candidate", "family", "k", "ll", "bic", "dbic", "dbic_like_for_like", "parameters", "s", "best_start",
                     "bic_spread_starts"]])

    def read_md(rd):
        return (f"winner **{rd['winner']}** ({rd['family']}); runner-up {rd['runner_up']} at +{rd['margin']:.1f} BIC; best "
                f"candidate outside the winner's family: {rd['next_family']} at +{rd['margin_family']:.1f}; F1 at "
                f"+{rd['f1_minus_winner']:.1f}")

    q = S["quote"]
    rate = (f"is constant (slope {slope:+.2f} +- {q['slope']})" if abs(slope) < S["reading"]["flat_abs_slope"] else
            f"{'rises' if slope > 0 else 'falls'} with slope {slope:.2f} +- {q['slope']}")
    F0 = fair_read[("main", WIN)]
    n_terms = len(j2["main"]["esindy"])
    sentence = (f"Between {W_LO} and {M['w_hi']:.2f} dB, in the first hour after rain, the drying rate per dB of water {rate} "
                f"on the main set (bootstrap refits {np.percentile(b_slopes, 5):+.2f}..{np.percentile(b_slopes, 95):+.2f}; "
                f"{slopes['no_rise_heavy'][0]:+.2f} without the rise-heavy events), a shape no fake truth had; the data prefer "
                f"family \"{R0['family']}\" by BIC (margin {R0['margin_family']:.1f} to {R0['next_family']}), but that "
                f"ranking against F3 / F5 / F6 rests on {neg_start[('main', WIN)]} rising tails fitted with a negative start, "
                f"and like for like the winner is {F0['winner']} ({F0['family']}) by {F0['margin_family']:.1f}, not decisive "
                f"with residuals {q_rows[0]['whitened_variance_ratio']:.1f} times the dry-noise variance; the formula found is "
                f"not a law (E-SINDy keeps {n_terms} of {len(JTERMS)} library terms: h(W) = {formula_text(j2['main']['esindy'])}); "
                f"below {W_LO} dB the method cannot tell.")
    ev_tab = ix.event_id.value_counts().head(S["big_events"]).rename_axis("event").reset_index(name="tails")
    cinfo = pd.DataFrame([{"dataset": n, "knots_db": " ".join(f"{v:.3g}" for v in curve[n]["knots"]),
                           "k_at_knots_per_min": " ".join(f"{np.exp(v):.4f}" for v in curve[n]["x"][:-1]),
                           "penalty": curve[n]["lam_star"], "s": round(curve[n]["s"], 3), "bic": round(curve[n]["bic"], 1),
                           "best_start": curve[n]["lam"][str(FIT["lambda_start"])]["from"],
                           "start_cost_spread": round(max(v["cost"] for v in curve[n]["starts"].values())
                                                      - min(v["cost"] for v in curve[n]["starts"].values()), 1)} for n in names])
    grid_tab = pd.DataFrame({"W_db": wm, "k_per_min": k_main, "band_lo": lo, "band_hi": hi}).iloc[::5].round(4)
    sens_cols = ["dataset", "what", "tails", "events", "tails_big3", "w_hi_db", "menu_winner", "margin_family_bic",
                 "f1_negative_start_tails", "like_for_like_winner", "like_for_like_margin_family_bic",
                 "like_for_like_f1_minus_winner_bic", "f3_n_at_lower_bound", "f1_tau_min", "curve_slope",
                 "curve_slope_common_range", "curve_reading", "jump2_family", "agree"]
    j2_sets = jt[jt.tool.str.startswith(("E-SINDy", "PySR, its own"))][["dataset", "tool", "formula", "relative_rms_misfit"]]
    long_txt = "Cut for time." if not plan["long_window"] else (
        f"{read_md(reads[('main', WIN_LONG)])}.\n\n{menu_md('main', WIN_LONG)}")
    report = f"""# WAA_estm step 4: the drying law on the REAL tails

Written by `research/WAA_estm/step4_real_drying.py` (libraries `waa_pooled.py`, `waa_curve.py`, `waa_sindy.py`; settings
`waa.real4`, the fit settings are read unchanged from `waa.stress3b1c`, `waa.curve3b3`, `waa.curve3b3b`, `waa.curve3b3c`).
Frozen tail set v1, TRAIN tails only; the {int((pd.read_csv(P.TAILS_MODEL_INDEX).split_event == 'val').sum())} VAL tails were not read.

**Result.** {sentence}

## 1. What was run, cut, run time

- **Recipe** (only tools that passed on fake data): pooled likelihood of 3b-3b (offset prior, ARMA(1,1) noise, speed spread
  marginalised), first {WIN} min of each tail. (a) Menu fits {', '.join(S['candidates'])}: 4 starts each (3 generic + the
  F1 solution; F1: 4 generic; {S1C['maxfev_generic']} evaluations), BIC with N = distinct (link, event) pairs. (b) Free
  curve: {S3B['knots']} knots from {W_LO} dB to W_hi, {len(FIT['starts'])} starts, penalty by 2-fold cross-validation over
  events, read only for W >= {W_LO} dB, no step tool; band = {len(boots)} event-bootstrap refits, each from
  {len(S['bootstrap_starts'])} starts ({', '.join(S['bootstrap_starts'])}). (c) Jump 2 on that curve. (d) Literature.
- **Fits**: {len(fits)} menu fits ({len(cells)} dataset-windows), {len(curve)} free-curve fits, {len(boots)} bootstrap refits
  (best start: {', '.join(f'{k} x{v}' for k, v in pd.Series([b['from'] for b in boots]).value_counts().items())}).
- **Timing**: one likelihood evaluation on the main set {plan['t_menu_ms']:.0f} ms (menu law) / {plan['t_curve_ms']:.0f} ms
  (free curve); one F1 fit {plan['t_f1_fit_s']} s. Projection {plan['projected_full_s'] / 3600:.2f} h on {workers} processes
  (budget {S['time_budget_h']} h). **Cut: {'; '.join(plan['cuts']) if plan['cuts'] else 'nothing'}.**
- **Run time**: {plan['wall_s'] / 60:.0f} min wall ({fresh / 3600:.1f} process-hours; cached in `.cache/fits4/`). This run:
  {time.time() - clock:.0f} s.
- **Software**: numpy / scipy / numba (own code); PySR {('used for ' + str(len(have)) + ' curves') if have else 'NOT run'}.
- **Added after the first run, from the cached fits (no fit was rerun, no setting changed)**: the like-for-like menu
  reading (row `F1+`, section 3), the model-free rate table and the share of tails that do not decay (section 4), the
  residual rows of the set without the rise-heavy events (check 4).

## 2. Real W range, tails, events, links

Main set: **{M['tails']} tails, {M['events']} events, {M['links']} links** ({M['sublinks']} sub-links);
{curve['main']['n_tails']} tails have >= {FIT['min_points']} valid minutes in the window; N of the BIC =
{curve['main']['n_independent']} (link, event) pairs. **W_hi = {M['w_hi']:.2f} dB** (the level with
{FIT['support_min_minutes']} tail-minutes of A above it in the first {WIN} min); {M['minutes_in_range']} of {M['points']}
tail-minutes lie inside [{W_LO}, W_hi]. Knots: {', '.join(f'{v:.2f}' for v in M['knots'])} dB. The 3 big events:

{md(ev_tab)}

## 3. Menu fits, main set, first {WIN} min (`{rel(P.STEP4_MENU)}`)

{menu_md('main', WIN)}

**Reading (3b-1c rules: F3 / F5 one family; F4 vs F6 by < {S['f4_f6_bic']:g} BIC not decisive): {read_md(R0)}.**

**Found after the run: this ranking is not like for like.** F1 and F4 are linear in the start W0, and `waa_pooled.py`
profiles it without a bound; F3, F5, F6 and the free curve keep W0 >= 0.001 dB. {neg_start[('main', WIN)]} of the
{curve['main']['n_tails']} tails (rising tails) get a NEGATIVE start under F1, and they alone are worth the whole gap
between F1 and F3 at n = 1, which is the same law (log-likelihood {best_of('main', WIN, 'F1').ll:.0f} vs
{best_of('main', WIN, 'F1+').ll:.0f}). Row `F1+` = the exponential with W0 >= 0; column `dbic_like_for_like` compares
F1+, F3, F5, F6 only. **Like for like: winner {F0['winner']} ({F0['family']}); next family {F0['next_family']} at
+{F0['margin_family']:.1f} BIC; F1+ at +{F0['f1_minus_winner']:.1f}.** F3 ends at its lower bound n = 1 and F5 at
b = {F0['f5_b']:.1g}: the data do not want a rate that RISES with W. On fake data no tail rises, so this did not show.

- `dbic` = BIC minus the best BIC. `bic_spread_starts` = best minus worst start of that candidate.
- F6's A_c is a fitted number, **not a physical level** (3b-3c: the likelihood prefers a feature near 0.2-0.4 dB
  whatever the truth). F5's linear rate a is not interpretable from the first hour (3b-1c).
- Fitted speed spread s: {', '.join(f"{r.candidate} {r.s}" for r in menu[(menu.dataset == 'main') & (menu.window_min == WIN)].itertuples())}
  (the fake data were made with s = {S1C['s_true']}).

### (d) Literature

{LITERATURE.format(tau=json.loads(f1.theta)[0])}

## 4. Free curve (`{rel(P.STEP4_CURVES)}`)

![real k(W)](../figures/{P.STEP4_FIG_K.name})

*Band (hatched): {len(boots)} event-bootstrap refits, 5-95 %. On fake data this band was about 3 times too narrow.*

- **Slope of log k vs log W over [{W_LO}, {M['w_hi']:.2f}] dB: {slope:+.2f}** (quoted +-{q['slope']}; bootstrap refits
  {np.percentile(b_slopes, 5):+.2f}..{np.percentile(b_slopes, 95):+.2f}, min-max {b_slopes.min():+.2f}..{b_slopes.max():+.2f}).
  Reading by the rule written before the run: **{S0.curve_reading}**.
- Level: k = {k_main[0]:.4f} / min at {W_LO} dB, {np.exp(np.interp(0.0, np.log(wm), np.log(k_main))):.4f} at 1 dB,
  {k_main[-1]:.4f} at {wm[-1]:.2f} dB (1 / k = {1 / k_main[0]:.0f}, {1 / np.exp(np.interp(0.0, np.log(wm), np.log(k_main))):.0f},
  {1 / k_main[-1]:.0f} min); quoted +-{q['level_pct']} %. Bootstrap band half-width (median over the range): +-{band_half:.0%}.
- Information only: BIC(free curve) - BIC(menu winner) = {S0.curve_minus_winner_bic:+.1f} (the free curve had a lower BIC
  than the TRUE law on fake data too, 3b-3; not a test).

{md(grid_tab)}

{md(cinfo)}

`start_cost_spread` = worst minus best of the 4 starts after 150 evaluations (fake data, 3b-3b: median 8).

**Model-free look (information, not a tool of the recipe).** Main set, first {WIN} min: A smoothed over
{S['raw_rate']['smooth_min']} min per tail, dA/dt by central difference, k = -mean(dA/dt) / median(A) per level bin.
No offset, no speed factor, no noise model: only a check that the falling trend is in the data and not made by the fit.

{md(raw)}

Tails that do not decay (A / A_start, mean of A over 5 min):

{md(stay)}

## 5. Jump 2 on the curve above {W_LO} dB (`{rel(P.STEP4_JUMP2)}`)

Library {', '.join(JTERMS)} (sw x = max(W - x dB, 0)); h(W) = W k(W) on {FIT['grid_points']} log-spaced points.
`relative_rms_misfit` = rms of formula / curve - 1. Main set:

{md(jt[jt.dataset == 'main'].drop(columns=['dataset', 'terms']))}

- **E-SINDy (main tool): h(W) = {formula_text(j2['main']['esindy'])}** ({S0.jump2_family}).
- The same elimination on each of the {len(boots)} bootstrap curves: {', '.join(f'"{k}" x{v}' for k, v in boot_formulas.items())}.
- Best-subset is read as "true terms + at most one extra" (3b-3c): the curve is smooth and has no noise of its own,
  so BIC always prefers more terms; the misfit column says how much each extra term buys.

## 6. Agreement between menu, curve and jump 2

| tool | says | class |
|---|---|---|
| menu fits (BIC), as run | {R0['winner']}: {R0['family']} | {agree_class(R0['family'])} |
| menu fits, like for like (W0 >= 0) | {F0['winner']}: {F0['family']}; F3 at n = 1, F5 at b = 0 | {agree_class(F0['family'])} |
| free curve (slope {slope:+.2f}) | {S0.curve_reading} | {agree_class(S0.curve_reading)} |
| jump 2 (E-SINDy) | {formula_text(j2['main']['esindy'])}: {S0.jump2_family} | {agree_class(S0.jump2_family)} |

**{'They agree' if S0.agree else 'They do NOT agree'}** (classes: exponential / quadratic family / other; "mildly
rising" + a multi-term formula + F4 or F6 counts as agreement on "other", which is what a switch or two layers looked
like on fake data).

## 7. Checks

### Checks 1 and 2: sensitivity sets and leaving out each big event (`{rel(P.STEP4_SUMMARY)}`)

{md(summ[sens_cols])}

Jump 2 per sensitivity set (E-SINDy formula; PySR's own pick):

{md(j2_sets[j2_sets.dataset.isin(SENS)])}

`curve_slope` is over each set's own [{W_LO}, W_hi]; `curve_slope_common_range` over [{W_LO}, {w_common[-1]:.2f}] dB.
`tails_big3` = tails from the 3 big events.

![sets](../figures/{P.STEP4_FIG_SETS.name})

### Check 3: menu fits on the first {WIN_LONG // 60} h (information only)

{long_txt}

### Check 4: fit quality (`{rel(P.STEP4_QUALITY)}`)

Residual = A minus the fitted curve of the tail (pooled law at the tail's speed factor, profiled start and offset).
`noise_model_*` = the same statistic of pure dry-weather noise drawn from each sub-link's noise model, on the same
minutes. `whitened_variance_ratio` = variance of the residuals after the noise filter of the likelihood, divided by
what the dry noise would give (1 = the residuals are exactly dry noise). Statistics inside a {WIN} min window are not
comparable with the control-tail row (tails about 12 h long): compare each row with its own `noise_model_*`.

{md(quality)}

![gallery](../figures/{P.STEP4_FIG_GALLERY.name})

### Check 5: real numbers against the fake data (`{rel(P.STEP4_FAKE)}`)

Slope over the read range, and median A / A_start at 10 / 30 / 60 min (mean of A over 5 min). `slope_inside` allows
+-{q['slope']}. Real ranges: bootstrap over events, 5-95 %.

{md(fake)}

## 8. The result in one sentence

**{sentence}**

## 9. In plain words: what this means for the WAA equation

{PLAIN}

## 10. Unclear / possibly wrong

{UNCLEAR}
"""
    save(report, P.STEP4_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {time.time() - clock:.0f} s")


# Written after the run (the tables are computed on every run; these texts are the reading of that run and quote its
# numbers: update them if a rerun changes the tables).
LITERATURE = """F1's fitted time constant is **tau = {tau:.1f} min** in the first hour (the median tail; the speed spread s = 1.40
moves single tails around it). Over the sets it is 13.9 min (alternative set), 15.3 min (without the rise-heavy events)
and 25-28 min (the others); on the first 6 h of the main set 18.7 min.

**No literature value can be compared: none is quoted in the repository, and none is invented here.** Checked:
`docs/proposal.md` names the dynamic model of Schleiss et al. [1] and the corrections of Pastorek et al. [9] without any
number; `pipeline/docs/openmrg_summary.md` only says that the RSL recovery lags the rain in the Torslanda case;
`pipeline/docs/reference/cml_libraries_howto.md` lists pycomlink's `waa_schleiss_2013(..., tau)` without a value;
`roadmap.md` names "Minda & Nakamura exponential drying", which is not in the proposal's reference list and has no
value either. The only paper stored in the repository is the OpenMRG paper [13]. To do before Step 8a: read the time
constants from [1] and [9] themselves and add them, with page numbers, to `roadmap.md`."""
PLAIN = """**The tools were built to answer "exponential, quadratic, or something in between". The real tails gave an answer
that was not on the list, and they also exposed a fault in the menu comparison. Step 4 does NOT deliver a drying law.**

1. **The curve falls.** In k(W) = "share of the remaining water that leaves per minute", fake truths were flat
   (exponential) or rising (quadratic, switch, two layers). The real curve goes the other way: 1 / k is about 6 min
   at 0.35 dB, 33 min at 1 dB and 71 min at 2.8 dB (slope -0.87). In words: while A is high it hardly moves, and once
   it is lower it drops fast.
2. **This is in the data, not made by the fit.** The model-free table shows the same trend (k from 0.052 / min at
   0.4 dB down to 0.016 above 1.4 dB), and the gallery shows why: some tails stay flat or RISE after t0_trace before
   they drop. 21 % of the tails are still at >= 90 % of their start level after 10 min, and 13 % are at or above it.
3. **It comes mostly from the events flagged event_rise_heavy** (164 of the 429 tails; 38 % of them are still at
   >= 90 % after 10 min, against 11 % of the others). Without those events the slope is -0.05, which the pre-set rule
   reads as "exponential in this range" (F1 tau 15 min). In all the other sets it is between -0.46 and -1.60.
4. **The curve itself is badly determined.** The four starts end 800-1700 apart in cost (fake data: about 8), 5 of
   the 20 bootstrap refits found a better fit from the flat start, the band is +-94 % wide, and the dip near 0.9 dB
   moves with the third knot. Only the overall trend (not rising) is read; the bumps are not.
5. **No formula.** E-SINDy keeps all 5 library terms (that is a description of a bumpy curve, not a law). The best
   one-term formula, h = 0.032 W (exponential), misses the curve by 61 % rms; PySR's simplest formula is a constant,
   h = 0.044 dB / min (49 %).
6. **The menu cannot decide either.** As run, F1 wins, but only because F1 and F4 may give a rising tail a negative
   start (3 tails, worth 2630 BIC). Like for like, F6 leads the exponential by 21 BIC on the main set (over the sets: from a tie, F1+ ahead by 0.9,
   to F6 ahead by 216), which is not decisive: the residuals are 4.1 times larger (in variance) than the dry-weather noise
   the likelihood assumes, so every BIC margin is overstated by about that factor, and F6 depends on its start (192
   BIC between starts). **One thing is stable in all 7 sets: F3 ends at n = 1 and F5 at b = 0.** The data never ask
   for a rate that rises with W, so the quadratic family (m3 / m5) is not supported in the first hour.
7. **What it means for the WAA equation.** dW/dt = -h(W) with one curve h and a random speed per tail does not
   describe the first hour of the TRAIN tails. A tail that stays flat or rises is not "drying slowly": something is
   still adding attenuation after the last rain trace (rain or drizzle below what radar and gauges see, re-wetting,
   condensation, or attenuation along the path). That is an input, i.e. the g(R, W) side of
   dW/dt = g(R, W) - h(W), or a wrong start time, not a property of h. Until this is separated, the exponential
   h(W) = W / tau with tau of 15-28 min is a placeholder, not a discovered law, and nothing can be said about
   switch vs two layers.
8. **Decisions for the learners before Step 6** (not taken here): (a) which tails define "drying": all of v1, or
   without the rise-heavy events (a data choice, made without looking at a model fit, so not circular), or tails
   re-started where the decay begins (a new frozen set); (b) repair the menu comparison (the same bound on W0 for every
   candidate) and re-run the fake-data tests with some rising tails; (c) the noise during drying is about twice the
   dry-weather noise (std 0.14 vs 0.09 dB in a 60 min window): the likelihood needs a larger or fitted noise level
   before BIC margins mean anything on real tails."""
UNCLEAR = """1. **The menu comparison has a fault that only real data showed.** `waa_pooled.py` profiles the start W0 of F1 and F4
   without a bound and keeps W0 >= 0.001 dB for F3, F5, F6 and the free curve. I report both readings. `F1+` is ONE
   likelihood evaluation of the F3 code at n = 1 with F1's tau and s (or F3's own fit where it ended at n = 1), not a
   re-optimised fit; F4 has no like-for-like version. The fake-data steps (3b-1b, 3b-1c, 3b-3c) used the same code;
   their tails never rise, so I expect no effect there, but I did not recheck. I did not change the library.
2. **"Falls with slope -0.87 +- 0.1" is the requested form, but +-0.1 is too optimistic here.** It comes from fake
   data whose curves were well determined. The honest spread on real data is the range over the sets (-1.60 to
   -0.05) and the bootstrap refits (-1.43 to -0.82). The same holds for "+-20 %" on the level.
3. **The free-curve fit may not be at its optimum.** Start spread 800-1700 in cost; the chosen penalty differs between
   the sets (0.01, 0.1, 1.0); the budgets are those of 3b-3b, unchanged, as required. A longer search could move the
   bumps and the slope.
4. **The dip of k near 0.9 dB is suspicious.** In every set the knot values go high, high, low, high, low; the low one
   is the third knot (0.91 dB on the main set, 0.68 dB on the alternative set), just below each set's A_start
   threshold (1.0 and 0.6 dB). That looks like an effect of the knots or of selecting tails by their start level,
   not physics. Not tested.
5. **The noise model does not fit drying tails**: whitened variance ratio 4.1 (first hour), 2.1 (6 h), 3.5 without
   the rise-heavy events. So the likelihood is misspecified, BIC margins are overstated, and the "free curve beats
   every menu law by about 10000 BIC" line mostly says that tails differ from each other in ways no single law
   covers. The `noise_model_*` columns ignore that a fit absorbs a start, an offset and a speed per tail, so the true
   gap is a little larger than shown.
6. **"Something still adds attenuation after t0_trace" is my interpretation, not a result.** The rise-heavy flag was
   set in step 2c from the rise rule, before any model fit, so dropping those events is not circular, but it is a
   choice about the data that the learners must make. It also removes 9 of the 24 events.
7. **The model-free table is rough**: it ignores the per-tail offset (it uses A, not W), smooths over 7 min, and bins
   by the same noisy A that enters the derivative. It supports the sign of the trend, not the numbers. I added it
   after seeing the curve; it is not one of the tools tested on fake data.
8. **Jump 2**: the E-SINDy terms come from the tested tool, but I refit their coefficients by least squares (the
   bagged medians, term by term, gave a formula 112 % off the curve when all five terms are kept). The library
   cannot write a falling k except with negative coefficients, and PySR cannot write max(). PySR's Pareto front
   starts with trivial formulas; its three simplest are listed as asked.
9. **Alternative set**: 119 of its sub-links have no control tails, so no noise parameters; they got the medians of
   the table (my choice). Its "F6 by 216 BIC" should be read with that in mind.
10. **Check 2**: without E20150826_0520 the slope is -1.60 and the menu is a tie (F1+ by 0.9); without E20150628_0530
    F6 leads by 106. So single events move both the curve and the menu a lot: 3 events hold 246 of 429 tails.
11. **Check 3 (6 h)**: F4 beats F1 by 5.5 BIC with tau_s = 16000 min and phi = 0.015, i.e. F1 plus a small constant
    floor. That is the offset trade-off seen in 3b-1, not a second layer. The lag-60 autocorrelation of check 4 can
    only come from this window.
12. **Check 5**: the three "recalibrated" truths of 3b-1c T3 were tuned to exactly these real numbers, so "inside" for
    m1 there mostly repeats 3a-fix (an exponential with spread reproduces the median curve). The real slope is outside
    every fake truth, which means the fake data did not contain the real behaviour: the tools were not tested on
    tails that plateau or rise.
13. **Cross-validation halves**: events are split 12 / 12 at random; with 3 big events the halves are unequal in
    tails. Seed fixed in `waa.real4.cv_seed`; not varied.
14. **Key figure**: `step4_k_real.png` is saved as diag (the module has 6 key outputs). If it is promoted, the output to
    demote is `outputs/tables/step3b1b_recovery.csv` (as proposed in 3b-2 and 3b-3b). My advice is to wait: the curve
    is not yet a result to show as "the law"."""


if __name__ == "__main__":
    main()
