"""Library of steps A6 and A7: how k is learned with the model of `ksurf.py` (the procedure of decisions.md, step A5).

    stage 1   the gauge-input set, rain = the 1 min gauge rate: the WETTING rung is selected (drying at K1, a nuisance)
    stage 2   all links, rain = the radar: k_wet held, the DRYING rung is selected

In each stage every rung K0..K3 is fitted with every smoothness weight lambda; a 5-fold cross-validation grouped by event
gives the held-out loss; lambda per rung = the smallest held-out loss; the selected rung = the simplest one whose held-out
loss differs from the best rung's by a paired 90 % event-bootstrap interval that contains zero (rule B1).

    f = folds(le, use)                               the fold of every link-event (events dealt round-robin by start time)
    s = stage(prob, le, use, "wetting", ...)         one stage: table of all rungs, the selected rung, its fit
    r = two_stage(D, le, curves, ...)                both stages
    d = lag_rule(D, le, curves, rain, use, ...)      the rule for d
"""
import numpy as np
import pandas as pd

from research.waa_ksurface import ksurf

C = ksurf.config()


def folds(le, use, n_folds=None):
    """The fold number of every link-event: the events that have a used link-event, sorted by start, dealt round-robin."""
    n_folds = n_folds or C["folds"]
    ev = le[np.asarray(use) > 0].groupby("event_id").window_start.min().sort_values().index
    fold_of = {e: i % n_folds for i, e in enumerate(ev)}
    return le.event_id.map(fold_of).fillna(-1).astype(int).to_numpy()


def _event_sums(ls, cnt, ev, n_ev, mask):
    return np.bincount(ev[mask], weights=ls[mask], minlength=n_ev), np.bincount(ev[mask], weights=cnt[mask], minlength=n_ev)


def paired_interval(Sa, Sb, N, reps, seed, interval):
    """The event-bootstrap interval of (mean loss a - mean loss b) from per-event loss sums Sa, Sb and counts N."""
    rng = np.random.default_rng(seed)
    idx = np.nonzero(N > 0)[0]
    draws = rng.integers(0, len(idx), (reps, len(idx)))
    num = (Sa[idx] - Sb[idx])[draws].sum(axis=1)
    return np.percentile(num / N[idx][draws].sum(axis=1), interval)


def stage(prob, le, use, select, other_rung, other_x=None, inc_F=None, inc_P=None, hold_ext=False, lambdas=None, log=print, tag=""):
    """One stage. select = the regime whose rung is selected ("wetting" or "drying"); the other regime has the rung
    `other_rung`: free if other_x is None, else held at the values (g, h) = other_x. inc_F, inc_P = the starting extension
    (held if hold_ext). Returns a dict: table (every rung and lambda), rung, lam, x, model, per-rung held-out sums."""
    lambdas = lambdas or C["lambdas"]
    use = np.asarray(use, float)
    fold = folds(le, use)
    ev = prob["event"]
    n_ev = int(ev.max()) + 1
    rows, best = [], {}
    prev = None
    for rung in ksurf.ORDER:
        for lam in lambdas:
            m = ksurf.Model(prob, rung if select == "wetting" else other_rung, other_rung if select == "wetting" else rung, lam)
            if not m.pen and lam != lambdas[0]:
                continue                                        # no knot values to smooth: lambda does nothing
            x0 = m.start(inc_F, inc_P) if prev is None else m.carry(*prev)
            if other_x is not None:
                g, h = ("gd", "hd") if select == "wetting" else ("gw", "hw")
                x0[m.sl[g]], x0[m.sl[h]] = other_x
                m.hold(g, h)
            if hold_ext:
                x0[m.sl["F"]], x0[m.sl["P"]] = np.concatenate(inc_F), np.concatenate(inc_P)
                m.hold("F", "P")
            x = m.fit(x0, use)
            info = dict(m.last)
            S, N = np.zeros(n_ev), np.zeros(n_ev)
            for f in range(fold.max() + 1):
                out = (fold == f) & (use > 0)
                xf = m.fit(x, use * (fold != f), maxiter=C["maxiter_fold"])
                ls, cnt = m.per_link_event(xf)
                s_, n_ = _event_sums(ls * use, cnt * use, ev, n_ev, out)
                S += s_
                N += n_
            cv = float(S.sum() / N.sum())
            rows.append({"select": select, "rung": rung, "lambda": lam, "parameters_free": int((~m.fixed).sum()), "train_loss": info["loss"], "cv_loss": cv,
                         "iterations": info["nit"], "converged": info["converged"]})
            log(f"   {tag}{select:8s} {rung} lambda {lam:g}: train {info['loss']:.5f}  cv {cv:.5f}  ({info['nit']} it)")
            if rung not in best or cv < best[rung]["cv"]:
                best[rung] = {"cv": cv, "lam": lam, "x": x, "model": m, "S": S, "N": N}
        prev = (best[rung]["model"], best[rung]["x"])
    tab = pd.DataFrame(rows)
    top = min(best, key=lambda r: best[r]["cv"])
    S_ = C["select"]
    chosen, iv = None, {}
    for rung in ksurf.ORDER:
        iv[rung] = paired_interval(best[rung]["S"], best[top]["S"], best[top]["N"], S_["boot_reps"], S_["seed"], S_["interval"]) if rung != top else (0.0, 0.0)
        if chosen is None and iv[rung][0] <= 0.0 <= iv[rung][1]:
            chosen = rung
    tab["best_lambda"] = [best[r]["lam"] == l for r, l in zip(tab.rung, tab["lambda"])]
    tab["minus_best_lo"] = [iv[r][0] if b else np.nan for r, b in zip(tab.rung, tab.best_lambda)]
    tab["minus_best_hi"] = [iv[r][1] if b else np.nan for r, b in zip(tab.rung, tab.best_lambda)]
    tab["lowest_cv"] = tab.rung == top
    tab["selected"] = (tab.rung == chosen) & tab.best_lambda
    b = best[chosen]
    log(f"   {tag}{select}: lowest cv {top}, selected {chosen} (lambda {b['lam']:g})")
    return {"table": tab, "rung": chosen, "lam": b["lam"], "x": b["x"], "model": b["model"], "lowest_cv": top}


def two_stage(D, le, band_idx, curves, gauge_mask, d_radar=0, d_gauge=0, inc_F=None, inc_P=None, hold_ext=False, use=None, log=print, tag=""):
    """The whole procedure. gauge_mask = the link-events of the gauge-input set; use = link-events to leave out get 0.
    Returns {"wet": stage 1, "dry": stage 2, "x": the final parameters, "model": the final model}."""
    use = np.ones(len(le)) if use is None else np.asarray(use, float)
    rain_g, ok5 = gauge_input(D)
    prob_g = ksurf.prepare(D, band_idx, curves, rain_g, d_gauge)
    prob_g["w5"] = prob_g["w5"] * ok5
    s1 = stage(prob_g, le, use * gauge_mask, "wetting", "K1", None, inc_F, inc_P, hold_ext, log=log, tag=tag)
    m1, x1 = s1["model"], s1["x"]
    wet = (x1[m1.sl["gw"]].copy(), x1[m1.sl["hw"]].copy())
    ext = m1.parts(x1)[4:]
    prob_r = ksurf.prepare(D, band_idx, curves, D["rain"], d_radar)
    s2 = stage(prob_r, le, use, "drying", s1["rung"], wet, [e.copy() for e in ext[0]], [e.copy() for e in ext[1]], hold_ext, log=log, tag=tag)
    return {"wet": s1, "dry": s2, "x": s2["x"], "model": s2["model"]}


def gauge_input(D):
    """The rain input of stage 1 and which 5 min values count: the 1 min gauge rate; a missing gauge minute gets the radar
    value, and a 5 min bin with a missing gauge minute is left out of the loss."""
    g = np.asarray(D["rain_gauge"], float)
    rain = np.where(np.isfinite(g), g, D["rain"])
    ok1 = np.concatenate([np.isfinite(g[o + ksurf.PADM:o + ksurf.PADM + n]) for o, n in zip(D["le_roff"], D["le_n"])])
    return rain, ok1.reshape(-1, 5).all(axis=1).astype(float)


def lag_rule(D, le, band_idx, curves, gauge, use, inc_F=None, inc_P=None, hold_ext=False, log=print, tag=""):
    """The rule for d (decisions.md, step A5): K0 in both regimes at every d of the grid; the d with the smallest loss, unless
    the paired 90 % event-bootstrap interval of loss(0) - loss(d_best) contains zero (then d = 0). gauge = True: the gauge
    input, else the radar. Returns (d, table)."""
    use = np.asarray(use, float)
    rows, per = [], {}
    rain, ok5 = gauge_input(D) if gauge else (D["rain"], 1.0)
    for d in C["lag_grid_min"]:
        prob = ksurf.prepare(D, band_idx, curves, rain, d)
        prob["w5"] = prob["w5"] * ok5
        m = ksurf.Model(prob, "K0", "K0")
        x0 = m.start(inc_F, inc_P)
        if hold_ext:
            m.hold("F", "P")
        x = m.fit(x0, use)
        ls, cnt = m.per_link_event(x)
        n_ev = int(prob["event"].max()) + 1
        per[d] = _event_sums(ls * use, cnt * use, prob["event"], n_ev, use > 0)
        rows.append({"d_min": d, "loss": float(per[d][0].sum() / per[d][1].sum()), "tau_wet_min": float(1 / m.k(x, "wetting", 0, 0)), "tau_dry_min": float(1 / m.k(x, "drying", 0, 0))})
        log(f"   {tag}d = {d}: loss {rows[-1]['loss']:.6f}")
    tab = pd.DataFrame(rows)
    d_best = int(tab.d_min[tab.loss.idxmin()])
    S_ = C["select"]
    lo, hi = paired_interval(per[0][0], per[d_best][0], per[0][1], S_["boot_reps"], S_["seed"], S_["interval"]) if d_best != 0 else (0.0, 0.0)
    d = 0 if lo <= 0.0 <= hi else d_best
    tab["lowest_loss"] = tab.d_min == d_best
    tab["loss0_minus_best_lo"], tab["loss0_minus_best_hi"] = lo, hi
    tab["chosen"] = tab.d_min == d
    return d, tab
