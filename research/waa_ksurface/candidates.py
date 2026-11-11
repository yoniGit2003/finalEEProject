"""Library of step E2: three ways to run the length test when the rain input has an error.

Sorting 5 min values into rain classes by a noisy rain mixes the classes: a bin sorted into a high class had, on average,
less rain than its radar value says, and the other way round in low classes. The measured F(R) then comes out too flat.

    plain(...)        (i)   the plain test: bins selected by each link's own radar rain
    instrument(...)   (ii)  bins selected by Z = the mean radar rain over all links of the event at that time (the same for
                            every link, so the selection error does not depend on the link); the link's own L stays in the line
    corrected(...)    (iii) the plain test, then a forward-bias correction: simulate the test on fake data built from the
                            current curves and the measured rain error, and move the cell values until the simulated test
                            returns what was observed

All three return the cell table of `ltest.cell_fits`. Definitions and settings: decisions.md (steps E3, E2), config `ltest`, `fake`.
"""
import numpy as np

from research.waa_ksurface import fakeworld as fw, ltest

LT, FK = ltest.config(), fw.config()


def plain(D, le, rule, n_short, mult=None):
    """(i) The plain length test. Returns (cells, draws)."""
    return ltest.cell_fits(ltest.steady_values(D, le, LT, rule), LT, n_short, mult, LT["interval"])


def instrument(D, le, rule, n_short):
    """(ii) Instrument binning: class and steadiness from the event's area-mean radar rain Z; the cell sits at the mean of the
    links' own rain."""
    z = fw.area_mean(D["r5"].astype(float), le)
    t = ltest.steady_values(D, le, LT, rule, sel=z, r_stat="mean")
    return ltest.cell_fits(t, LT, n_short, r_stat="mean")[0]


def proxy_rain(D, le, n_cells, sigma, scale_sd):
    """A stand-in for the TRUE rain, made from the observed rain: ln R is shrunk about the mean of its link-event, and the
    link-event means about the mean of their event, so that stand-in + rain error has the observed variances.

        x' = M + b (m - M) + a (x - m) + v / 2        x = ln R of a wet bin, m = its link-event mean, M = the event mean of m
        a^2 = 1 - (per-bin error variance) / (variance of x inside link-events)
        b^2 = 1 - (error variance of m) / (variance of m inside events)          both floored at shrink_min^2
        v = the total error variance (the errors have mean one, so observed logs sit v / 2 below the true ones)
    """
    R = FK["rain_error"]
    q, rho, lo = R["shared_share"], R["rho"], FK["correction"]["shrink_min"] ** 2
    n5 = D["le_n"] // 5
    rep = np.repeat(np.arange(len(le)), n5)
    r = D["r5"].astype(float)
    wet = np.isfinite(r) & (r > 0)
    x = np.log(np.where(wet, r, 1.0))
    cnt = np.bincount(rep, weights=wet.astype(float), minlength=len(le))
    use = cnt >= 3                                              # link-events with enough wet bins to have a mean and a spread
    m = np.bincount(rep, weights=np.where(wet, x, 0.0), minlength=len(le)) / np.maximum(cnt, 1.0)
    ok = wet & use[rep]
    v_within = float(((x - m[rep])[ok] ** 2).sum() / max(ok.sum() - use.sum(), 1))
    ev = fw.time_index(le)[1]
    n_use = np.bincount(ev, weights=use.astype(float))
    M = (np.bincount(ev, weights=np.where(use, m, 0.0)) / np.maximum(n_use, 1.0))[ev]
    v_between = float((((m - M)[use]) ** 2).sum() / max(use.sum() - (n_use > 0).sum(), 1))
    v_own = sigma ** 2 * (q + (1 - q) / n_cells)                # per-bin error variance of every link-event
    v_mean = scale_sd ** 2 + (1 - q) * sigma ** 2 / n_cells * (1 + rho) / ((1 - rho) * np.maximum(cnt, 1.0))
    a = np.sqrt(np.maximum(1.0 - v_own / v_within, lo))
    b = np.sqrt(np.maximum(1.0 - v_mean / v_between, lo))
    xp = M[rep] + b[rep] * (m - M)[rep] + a[rep] * (x - m[rep]) + 0.5 * (scale_sd ** 2 + v_own)[rep]
    return np.where(ok, np.exp(xp), r)


def corrected(D, le, n_cells, sigma, rule, n_short, seed, scale_sd=None):
    """(iii) The plain cells, then the forward-bias correction for a rain error of per-bin sd `sigma` and scale part
    `scale_sd`. Returns the cell table with the corrected I, s (the plain values in I_plain, s_plain)."""
    scale_sd = FK["rain_error"]["scale_sd"] if scale_sd is None else scale_sd
    cells, _ = plain(D, le, rule, n_short)
    cells = cells.assign(I_plain=cells.I, s_plain=cells.s)
    if sigma <= 0 and scale_sd <= 0:
        return cells
    K = FK["correction"]
    rng = np.random.default_rng(seed)
    proxy = proxy_rain(D, le, n_cells, sigma, scale_sd)
    n5 = D["le_n"] // 5
    rep = np.repeat(np.arange(len(le)), n5)
    band_bin = le.band_name.to_numpy()[rep]
    L_bin = le.L_km.to_numpy(float)[rep]
    el = cells[cells.eligible].sort_values(["band", "R_cell"])
    obs = {(r.band, r.cls): np.array([r.I, r.s]) for r in el.itertuples()}
    cur = {k: v.copy() for k, v in obs.items()}
    base = dict(D)
    base["le_c"] = np.zeros(len(le))
    base["w5"] = ((D["w5"] > 0) & np.isfinite(proxy)).astype(float)
    for _ in range(K["iterations"]):
        y = np.zeros(len(proxy))
        for band, g in el.groupby("band"):
            if len(g) < 2:
                continue
            rc = g.R_cell.to_numpy()
            f = ltest.monotone_curve(rc, [cur[(band, c)][0] for c in g.cls], np.ones(len(g)), "F")[0]
            p = ltest.monotone_curve(rc, [cur[(band, c)][1] for c in g.cls], np.ones(len(g)), "P")[0]
            mb = band_bin == band
            pr = np.nan_to_num(proxy[mb])
            y[mb] = f(pr) + L_bin[mb] * p(pr)
        acc = {k: [] for k in obs}
        for _s in range(K["sims"]):
            sim = dict(base)
            sim["y5"] = y
            sim["r5"] = proxy * np.exp(fw.rain_error(le, n_cells, sigma, FK, rng, scale_sd=scale_sd))
            tab = plain(sim, le, rule, n_short)[0].set_index(["band", "cls"])
            for k in obs:
                if k in tab.index and np.isfinite(tab.loc[k, "I"]):
                    acc[k].append([tab.loc[k, "I"], tab.loc[k, "s"]])
        for k in obs:
            if acc[k]:
                cur[k] = cur[k] + (obs[k] - np.mean(acc[k], axis=0))    # move the cell until the simulated test returns the observed value
    for k, v in cur.items():
        m = (cells.band == k[0]) & (cells.cls == k[1])
        cells.loc[m, "I"], cells.loc[m, "s"] = v[0], v[1]
    return cells
