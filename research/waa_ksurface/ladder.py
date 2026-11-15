"""Library of steps A6 and A7: how k is learned with the model of `ksurf.py` (the procedure of decisions.md, step A5).

    stage 1   the gauge-input set, rain = the 1 min gauge rate: the WETTING rung is selected (drying at K1, a nuisance)
    stage 2   all links, rain = the radar: k_wet held, the DRYING rung is selected

In each stage every rung K0..K3 is fitted with every smoothness weight lambda; a 5-fold cross-validation grouped by event
gives the held-out loss; lambda per rung = the smallest held-out loss. Selection (rule B1 with Deviation 5): the ladder is
walked from K0, and a richer rung replaces the current one only if its CV error is lower by MORE than 0.005 dB and the
paired 90 % event-bootstrap interval of the difference excludes zero. CV error (dB) = sqrt(2 x mean held-out Huber loss).
k_wet (Deviation 6): if the wetting ladder returns K0, k_wet is a bounded nuisance: the drying fit is made with k_wet at
its gauge-input estimate and again at its radar-input estimate.

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


def cv_db(S, N):
    """The CV error in dB: sqrt(2 x mean held-out Huber loss) (the RMS residual where the residuals are below the threshold)."""
    return float(np.sqrt(2.0 * S.sum() / N.sum()))


def paired_interval(Sa, Sb, N, reps, seed, interval):
    """The event-bootstrap interval of (CV error a - CV error b), in dB, from per-event loss sums Sa, Sb and counts N."""
    rng = np.random.default_rng(seed)
    idx = np.nonzero(N > 0)[0]
    draws = rng.integers(0, len(idx), (reps, len(idx)))
    n = N[idx][draws].sum(axis=1)
    return np.percentile(np.sqrt(2.0 * Sa[idx][draws].sum(axis=1) / n) - np.sqrt(2.0 * Sb[idx][draws].sum(axis=1) / n), interval)


def _huber(r, delta):
    a = np.abs(r)
    return np.where(a <= delta, 0.5 * r * r, delta * (a - 0.5 * delta))


def stage(prob, le, use, select, other_rung, other_x=None, inc_F=None, inc_P=None, hold_ext=False, lambdas=None, log=print, tag="", score=None, margin_db=None, selector=None, rungs=None):
    """One stage. select = the regime whose rung is selected ("wetting" or "drying"); the other regime has the rung
    `other_rung`: free if other_x is None, else held at the values (g, h) = other_x. inc_F, inc_P = the starting extension
    (held if hold_ext). Returns a dict: table (every rung and lambda), rung, lam, x, model, per-rung held-out sums.
    Deviation 9: score = a boolean per 5 min value; if given, the held-out error is taken on these values only (the fits
    are still trained on every counted value). margin_db = the margin of the selection (default: config, Deviation 5).
    Deviation 10: selector(compare, margin) -> (rung, attribution, steps) replaces the walk from K0; compare(richer, simpler)
    = (gain in dB, lower and upper end of its paired interval). The result then also holds "attribution".
    Deviation 11: rungs = the rungs that are fitted and compared, from the simplest (default: K0, K1, K2, K3)."""
    rungs = list(rungs or ksurf.ORDER)
    lambdas = lambdas or C["lambdas"]
    bin_le = np.repeat(np.arange(len(prob["le_n"])), prob["le_n"] // 5)
    use = np.asarray(use, float)
    fold = folds(le, use)
    ev = prob["event"]
    n_ev = int(ev.max()) + 1
    rows, best = [], {}
    prev = None
    for rung in rungs:
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
                if score is None:
                    ls, cnt = m.per_link_event(xf)
                    s_, n_ = _event_sums(ls * use, cnt * use, ev, n_ev, out)
                else:                                           # only the scored values of the held-out link-events
                    b = out[bin_le] & score
                    h = _huber(m.simulate(xf)[1][b] - prob["target"][b], m.cfg["huber_delta_db"])
                    s_, n_ = np.bincount(ev[bin_le[b]], weights=h, minlength=n_ev), np.bincount(ev[bin_le[b]], minlength=n_ev).astype(float)
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
    tab["cv_db"] = np.sqrt(2.0 * tab.cv_loss)
    top = min(best, key=lambda r: best[r]["cv"])
    S_ = C["select"]
    margin = S_["margin_db"] if margin_db is None else margin_db
    chosen, steps, attribution = "K0", {}, None

    def compare(richer, simpler):
        gain = cv_db(best[simpler]["S"], best[simpler]["N"]) - cv_db(best[richer]["S"], best[richer]["N"])
        lo, hi = paired_interval(best[simpler]["S"], best[richer]["S"], best[richer]["N"], S_["boot_reps"], S_["seed"], S_["interval"])
        return gain, lo, hi

    if selector is not None:                                    # Deviation 10: K1 and K2 are each compared with K0
        chosen, attribution, steps = selector(compare, margin)
    else:
        for rung in rungs[1:]:                                  # walk the ladder: a richer rung must beat the current one by more than the margin
            gain, lo, hi = compare(rung, chosen)
            take = bool(gain > margin and not (lo <= 0.0 <= hi))
            steps[rung] = (chosen, gain, lo, hi, take)
            if take:
                chosen = rung
    tab["best_lambda"] = [best[r]["lam"] == l for r, l in zip(tab.rung, tab["lambda"])]
    for k, name in enumerate(("compared_with", "gain_db", "gain_lo", "gain_hi", "replaces")):
        tab[name] = [steps[r][k] if (b and r in steps) else None for r, b in zip(tab.rung, tab.best_lambda)]
    tab["lowest_cv"] = tab.rung == top
    tab["selected"] = (tab.rung == chosen) & tab.best_lambda
    b = best[chosen]
    log(f"   {tag}{select}: lowest cv {top}, selected {chosen} (lambda {b['lam']:g})")
    return {"table": tab, "rung": chosen, "lam": b["lam"], "x": b["x"], "model": b["model"], "lowest_cv": top, "S": b["S"], "N": b["N"], "attribution": attribution}


def plain_k0(prob, use, inc_F, inc_P, hold_ext):
    """A K0 / K0 fit with both regimes free: (1 / k_wet, 1 / k_dry, model, x)."""
    m = ksurf.Model(prob, "K0", "K0")
    if hold_ext:
        m.hold("F", "P")
    x = m.fit(m.start(inc_F, inc_P), use)
    return float(1 / m.k(x, "wetting", 0, 0)), float(1 / m.k(x, "drying", 0, 0)), m, x


def two_stage(D, le, band_idx, curves, gauge_mask, d_radar=0, d_gauge=0, inc_F=None, inc_P=None, hold_ext=False, use=None, log=print, tag="", stage2=None, wet_rungs=None):
    """The whole procedure. gauge_mask = the link-events of the gauge-input set; use = link-events to leave out get 0.
    Returns {"wet": stage 1, "dry": stage 2, "x": the final parameters, "model": the final model}; if the wetting ladder
    returned K0 (Deviation 6) also "tau_wet" = {"gauge", "radar"} (1 / k_wet from each input, min) and "x_alt" = the selected
    drying rung refitted with k_wet at the radar-input estimate.
    stage2 (Deviation 9) = {"w5": the 5 min values that count in stage 2, "score": the values the drying rungs are scored
    on, "margin_db": the margin of the drying selection, "selector": the selection rule (Deviation 10; optional)};
    None = the rule of Deviation 5 on every counted value. Deviation 11: wet_rungs = the wetting rungs of stage 1 and
    stage2["rungs"] = the drying rungs of stage 2 (default: the whole ladder)."""
    use = np.ones(len(le)) if use is None else np.asarray(use, float)
    rain_g, ok5 = gauge_input(D)
    prob_g = ksurf.prepare(D, band_idx, curves, rain_g, d_gauge)
    prob_g["w5"] = prob_g["w5"] * ok5
    s1 = stage(prob_g, le, use * gauge_mask, "wetting", "K1", None, inc_F, inc_P, hold_ext, log=log, tag=tag, rungs=wet_rungs)
    m1, x1 = s1["model"], s1["x"]
    wet = (x1[m1.sl["gw"]].copy(), x1[m1.sl["hw"]].copy())
    ext = m1.parts(x1)[4:]
    prob_r = ksurf.prepare(D, band_idx, curves, D["rain"], d_radar)
    extra = {}
    if stage2 is not None:
        prob_r["w5"] = np.asarray(stage2["w5"], float)
        extra = {"score": np.asarray(stage2["score"], bool), "margin_db": stage2["margin_db"], "selector": stage2.get("selector"), "rungs": stage2.get("rungs")}
    s2 = stage(prob_r, le, use, "drying", s1["rung"], wet, [e.copy() for e in ext[0]], [e.copy() for e in ext[1]], hold_ext, log=log, tag=tag, **extra)
    out = {"wet": s1, "dry": s2, "x": s2["x"], "model": s2["model"]}
    if s1["rung"] == "K0":                                      # Deviation 6: k_wet is a bounded nuisance, carried at both its estimates
        tau_r = plain_k0(prob_r, use, inc_F, inc_P, hold_ext)[0]
        m2 = s2["model"]
        x0 = s2["x"].copy()
        x0[m2.sl["gw"]] = np.log(1.0 / tau_r)
        out["x_alt"] = m2.fit(x0, use)                          # gw is held in this model, so it stays at the radar-input estimate
        out["tau_wet"] = {"gauge": float(1 / m1.k(x1, "wetting", 0, 0)), "radar": tau_r}
        log(f"   {tag}1 / k_wet: gauge input {out['tau_wet']['gauge']:.2f} min, radar input {tau_r:.2f} min")
    return out


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
