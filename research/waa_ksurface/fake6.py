"""Library of step A6: the fake world in which the k-surface fit must find a KNOWN k before it is used on real data.

What is real: the links, the storms, the 5 min rain of every link-event. What is made up, with round numbers:

    true rain per minute   the 5 min rain joined by straight lines, with a lognormal fluctuation inside the 5 minutes
    F_true, P_true         3 R / (R + 2) dB and 0.1 R / 0.2 R dB/km, handed to the model as "measured" cells plus knots
    W                      one of the truths of the A6 table (constant k, k_dry(W), k(R), two stores, ...)
    A - c                  W + L P_true(R) + offset error + white noise, rounded to the 0.3 dB grid
    radar input            the true 5 min rain x the T1 error (per event, shared, the link's own; slow)
    gauge input            the true 1 min rain x (own error + representativeness), in whole gauge steps

Every size is this branch's own measurement or a stated round number: decisions.md, step A6; config `a6`.
"""
import numpy as np
from numba import njit
from scipy.signal import lfilter

from research.waa_ksurface import data, fakeworld as fw, fp, ksurf
from research.wet_antenna.system_id import PADM

A6 = data.branch_config()["a6"]
A5 = ksurf.config()
FK = fw.config()
T3 = data.branch_config()["t3"]
KIND = {"const": 0, "dry_w": 1, "both_r": 2, "two_store": 3}


# =================================================================================================== F and P of the fake world
def fake_curves(M):
    """([Curves per band], [true F increments], [true P increments]): the fake truth in the layout of the real measurement:
    cells at the R of the real valid cells, and the extension knots at the true values."""
    curves, iF, iP = [], [], []
    for band, c in M.items():
        k = fp.from_values(band, c.F_r, fw.f_true(c.F_r, FK), c.P_r, fw.p_true(band, c.P_r, FK))
        curves.append(k)
        iF.append(ksurf.increments_of(k, "F", fw.f_true(np.asarray(A5["f_knots_mmh"], float), FK)))
        iP.append(ksurf.increments_of(k, "P", fw.p_true(band, np.asarray(A5["p_knots_mmh"], float), FK)))
    return curves, iF, iP


# =================================================================================================== the true k
def k_true(spec, regime, R, W):
    """The true rate (1 / min) of a truth at rain R and level W (not defined for the two-store truth)."""
    R, W = np.broadcast_arrays(np.asarray(R, float), np.asarray(W, float))
    kind = spec["kind"]
    if kind == "const":
        return np.full(R.shape, spec["k_wet"] if regime == "wetting" else spec["k_dry"])
    if kind == "dry_w":
        return np.full(R.shape, spec["k_wet"]) if regime == "wetting" else spec["k_dry0"] * (1 + W / spec["w_scale_db"])
    if kind == "both_r":                                        # the dependence acts between r0 and r_cap (default: everywhere)
        Re = np.clip(R, spec.get("r0_mmh", 0.0), spec.get("r_cap_mmh", 1e9)) - spec.get("r0_mmh", 0.0)
        return spec["k0"] * (1 + Re / spec["r_scale_mmh"]) if regime == "wetting" else spec["k0"] / (1 + Re / spec["r_scale_mmh"])
    return np.full(R.shape, np.nan)


def true_rung(spec, regime):
    if spec.get("true_rung"):                                   # a truth may name its rung itself (d-high: constant inside the measured range)
        return spec["true_rung"]
    return {"const": "K0", "dry_w": "K0" if regime == "wetting" else "K1", "both_r": "K2", "two_store": "none"}[spec["kind"]]


@njit(cache=True)
def _truth(kind, par, v, r, le_n, boff1):
    """W after every minute for the true dynamics (the same one-step-per-minute rule as the model)."""
    out = np.empty(len(v))
    for e in range(len(le_n)):
        w, ws = 0.0, 0.0                                        # (fast) store, slow store
        for t in range(le_n[e]):
            q = boff1[e] + t
            if kind == 3:                                       # two stores, each a share of F
                vf, vsl = par[0] * v[q], (1.0 - par[0]) * v[q]
                w = vf + (w - vf) * np.exp(-(par[1] if vf > w else par[2]))
                ws = vsl + (ws - vsl) * np.exp(-par[3])
                out[q] = w + ws
                continue
            wet = v[q] > w
            if kind == 0:
                k = par[0] if wet else par[1]
            elif kind == 1:
                k = par[0] if wet else par[1] * (1.0 + w / par[2])
            else:
                re = min(max(r[q], par[2]), par[3]) - par[2]
                k = par[0] * (1.0 + re / par[1]) if wet else par[0] / (1.0 + re / par[1])
            w = v[q] + (w - v[q]) * np.exp(-k)
            out[q] = w
    return out


def _params(spec):
    k = spec["kind"]
    if k == "const":
        return np.array([spec["k_wet"], spec["k_dry"]])
    if k == "dry_w":
        return np.array([spec["k_wet"], spec["k_dry0"], spec["w_scale_db"]])
    if k == "both_r":
        return np.array([spec["k0"], spec["r_scale_mmh"], spec.get("r0_mmh", 0.0), spec.get("r_cap_mmh", 1e9)])
    return np.array([spec["phi"], spec["k_fast_wet"], spec["k_fast_dry"], spec["k_slow"]])


# =================================================================================================== the rain
def rain_per_minute(D, rng):
    """The true 1 min rain of every link-event over its padded span (layout of D["rain"]): the 5 min values joined by
    straight lines between bin centres, times a lognormal AR(1) fluctuation, every bin rescaled to its 5 min value."""
    R = A6["rain_1min"]
    base = np.asarray(D["rain"], float)
    z = lfilter([np.sqrt(1 - R["rho"] ** 2)], [1.0, -R["rho"]], rng.standard_normal(len(base))) * R["sd"]
    out = np.empty(len(base))
    for o, n in zip(D["le_roff"], D["le_n"]):
        m = n + 2 * PADM
        c5 = base[o:o + m:5]
        lin = np.interp(np.arange(m) + 0.5, 2.5 + 5.0 * np.arange(len(c5)), c5)
        r = (lin * np.exp(z[o:o + m])).reshape(-1, 5)
        mean = r.mean(axis=1)
        out[o:o + m] = (r * np.where(mean > 0, c5 / np.where(mean > 0, mean, 1.0), 0.0)[:, None]).ravel()
    return out


def window(D, x):
    """The values of a padded 1 min array inside the windows (the layout of D["y1"])."""
    return np.concatenate([x[o + PADM:o + PADM + n] for o, n in zip(D["le_roff"], D["le_n"])])


def radar_seen(D, le, n_cells, rng):
    """The radar input of the fit (padded, held over 5 min): the true 5 min rain x the slow T1 error."""
    err = fw.rain_error_slow(le, n_cells, A6["sigma_c"], T3["slow_error"], rng)
    base = np.asarray(D["rain"], float)
    out = np.empty(len(base))
    for o, n, b in zip(D["le_roff"], D["le_n"], D["le_boff5"]):
        e = err[b:b + n // 5]
        e = np.r_[np.full(PADM // 5, e[0]), e, np.full(PADM // 5, e[-1])]          # the padding takes the error of the nearest bin
        out[o:o + n + 2 * PADM] = base[o:o + n + 2 * PADM] * np.repeat(np.exp(e), 5)
    return out


def gauge_seen(D, le, r1, step_mm, dist_km, rng):
    """The gauge input of the fit (padded, 1 min): the true 1 min rain x exp(own error + representativeness), per 5 min bin,
    then whole gauge steps per minute. Link-events without a gauge get NaN."""
    G = T3["gauge"]
    out = np.full(len(r1), np.nan)
    sd_own = G["own_sd_15min"] * np.sqrt(3.0)
    for i, (o, n) in enumerate(zip(D["le_roff"], D["le_n"])):
        if not np.isfinite(step_mm[i]):
            continue
        m = n + 2 * PADM
        nb = m // 5
        sd_rep = np.sqrt(G["variogram_slope_per_km"] * (dist_km[i] + D["le_L"][i] / 3.0) / G["f15"])
        e_rep = lfilter([np.sqrt(1 - G["rho"] ** 2)], [1.0, -G["rho"]], rng.standard_normal(nb)) * sd_rep
        err = sd_own * rng.standard_normal(nb) + e_rep - 0.5 * (sd_own ** 2 + sd_rep ** 2)
        rate = r1[o:o + m] * np.repeat(np.exp(err), 5)
        start = rng.uniform(0, 1)
        cum = start * step_mm[i] + np.cumsum(rate / 60.0)       # mm
        out[o:o + m] = np.diff(np.r_[np.floor(start), np.floor(cum / step_mm[i])]) * step_mm[i] * 60.0
    return out


# =================================================================================================== a fake dataset
def make(D, le, n_cells, spec, curves, inc_F, inc_P, seed, step_mm, dist_km, noise=True, exact_gauge=False, exact_radar=False):
    """A fake dataset for one truth. Returns (Dk, extra): Dk = the dataset as the fit reads it (y5 = fake A - c, rain = the
    radar input, rain_gauge = the gauge input); extra = the truth behind it (true rain, true W, changed events).
    The three switches are for the diagnostic only: noise = False removes the offset error, the link shift, the white noise
    and the rounding; exact_gauge / exact_radar hand the fit the true rain instead of the input with its error. The random
    numbers are drawn in the same order whatever the switches, so the truth is the same in every variant."""
    N = A6["noise"]
    rng = np.random.default_rng(seed)
    n, rep = D["le_n"], np.repeat(np.arange(len(le)), D["le_n"])
    r1_pad = rain_per_minute(D, rng)
    ghost_events = []
    if "ghost" in spec:                                         # rain that goes on after the last wet bin, unseen by the input
        g = spec["ghost"]
        ev_names = np.sort(le.event_id.unique())
        ghost_events = list(rng.choice(ev_names, g["events"], replace=False))
        for i in np.nonzero(le.event_id.isin(ghost_events).to_numpy())[0]:
            a = D["le_roff"][i] + PADM + (D["le_last"][i] + 1) * 5
            b = min(a + g["minutes"], D["le_roff"][i] + n[i] + 2 * PADM)
            r1_pad[a:b] = np.maximum(r1_pad[a:b], g["rain_mmh"])
    idx = np.concatenate([o + PADM + np.arange(k) for o, k in zip(D["le_roff"], n)])
    r1 = r1_pad[idx]
    band = le.band_name.map({c.band: i for i, c in enumerate(curves)}).to_numpy()
    v, p = np.empty(len(r1)), np.empty(len(r1))
    d_true = int(spec.get("d_true", 0))
    for b, c in enumerate(curves):
        m = band[rep] == b
        v[m] = ksurf.full_curve(c, "F", inc_F[b], r1[m])
        p[m] = ksurf.full_curve(c, "P", inc_P[b], r1_pad[idx[m] - d_true])
    link_shift = rng.normal(0, spec.get("link_sd_db", N["link_sd_db"]), D["nl"])[D["le_link"]] * float(noise)
    v = np.where(r1 > 0, np.maximum(v + link_shift[rep], 0.0), 0.0)        # the level of a link while it rains
    W = _truth(KIND[spec["kind"]], _params(spec), v, r1, n.astype(np.int64), D["le_boff1"].astype(np.int64))
    y = W + D["le_L"][rep] * p + float(noise) * (rng.normal(0, N["offset_sd_db"], len(le))[rep] + N["white_1min_sd_db"] * rng.standard_normal(len(W)))
    q, phase = N["quant_step_db"], rng.uniform(0, N["quant_step_db"], len(le))[rep]
    if noise:
        y = q * np.round((y + phase) / q) - phase               # the 0.3 dB grid of the instrument
    Dk = dict(D)
    Dk["y1"], Dk["y5"] = y, y.reshape(-1, 5).mean(axis=1)
    Dk["le_c"] = np.zeros(len(le))
    Dk["w5"] = (D["w5"] > 0).astype(float)
    Dk["rain"] = radar_seen(D, le, n_cells, rng)
    Dk["rain_gauge"] = gauge_seen(D, le, r1_pad, step_mm, dist_km, rng)
    if exact_radar:
        Dk["rain"] = np.asarray(D["rain"], float).copy()        # the true 5 min rain
    if exact_gauge:
        Dk["rain_gauge"] = np.where(np.isfinite(Dk["rain_gauge"]), r1_pad, np.nan)      # the true 1 min rain
    return Dk, {"r1_true": r1, "W_true": W, "v_true": v, "ghost_events": ghost_events, "band_idx": band}
