"""Tail simulator (library): drying laws for the WAA state W(t) and an observation model that mimics tail set v1.

Used by calibrate_simulator.py (calibration + checks) and by step 3b (which laws can the data tell apart?).
Nothing here reads files; the caller passes the skeleton (one row per tail) and the settings (`wet_antenna.simulator`).

Truth laws (t = minutes since t0_trace, W(0) = A_start of the tail; closed-form solutions of the ODEs):
    m1  dW/dt = -W / tau                          W = A0 exp(-t / tau)
    m3  dW/dt = -c W^n  (n = 2)                   W = (A0^(1-n) + (n - 1) c t)^(1 / (1 - n))
    m4  W = Wf + Ws, dWf/dt = -Wf / tau_f, dWs/dt = -Ws / tau_s, Ws(0) = phi A0
    m5  dW/dt = -a W - b W^2                      W = a A0 e^(-a t) / (a + b A0 (1 - e^(-a t)))
    m6  state switch (step 3b-1): dW/dt = -k_run (W - A_c) - k_evap W for W > A_c, else -k_evap W
        (run-off above a film level A_c, evaporation always). Above A_c: W = W* + (A0 - W*) e^(-(k_run + k_evap) t),
        W* = k_run A_c / (k_run + k_evap) < A_c, until W = A_c at t1; then W = A_c e^(-k_evap (t - t1)).
        A0 <= A_c: W = A0 e^(-k_evap t). A_c is in dB (absolute), so the shape depends on A_start.
Tail-to-tail spread (3a-fix, every law): the law's rate parameter(s) are scaled per tail (m1 tau, m3 c, m4 tau_f and
tau_s together, m5 a and b together, m6 k_run and k_evap together). For every law this is the same as stretching time: W_i(t) = W(t / f_i), with
the time-scale factor f_i = exp(s z_i), z_i = sqrt(1/2) z_link + sqrt(1/2) z_tail (half the variance per link).
(m2 of step 3a, "exponential with tau per tail", is m1 with spread.)

Observation model, per tail (each effect has a switch, so it can be turned off):
    10 s value  x = W(t / f) + offset + slow(minute) + sigma_white * eps
                offset ~ N(0, offset_std of the sub-link); eps ~ N(0, 1) per 10 s sample;
                slow = AR(1) at 1 min (3a-fix): slow_m = rho slow_(m-1) + sqrt(1 - rho^2) sigma_slow e_m,
                rho = exp(-1 / tau_noise), stationary start; constant within the minute
    quantize    q = phase + step * round((x - phase) / step)      step = 0.3 dB (RSL), phase ~ U(0, step) per tail
    1 min A     mean of the 6 quantized 10 s values (as the pipeline's 1 min tl_mean); NaN where the real tail is NaN
All random numbers are drawn whatever the switches are, so turning one effect off leaves the others identical.
Second slow component (step 3b-1c, misspecified noise): if a sub-link's noise tuple has 5 entries (sigma_white,
sigma_slow, tau_noise, sigma_slow2, tau_noise2), a second AR(1) wander is added to the first. It has its own random
stream, so every other random number stays the same as without it.
"""
import numpy as np
import pandas as pd
from scipy.signal import lfilter

LAWS = {"m1": "exponential", "m3": "nonlinear dW/dt = -c W^n", "m4": "two layers (fast + slow exponential)",
        "m5": "linear + quadratic dW/dt = -a W - b W^2", "m6": "state switch (run-off above A_c, evaporation)"}
PARAMS = {"m1": ["tau"], "m3": ["c", "n"], "m4": ["tau_f", "tau_s", "phi"], "m5": ["a", "b"],
          "m6": ["k_run", "k_evap", "A_c"]}


def law_curve(law, p, t, a0):
    """W(t) (dB) for the base law (no spread). t and a0 broadcast (arrays allowed)."""
    t = np.asarray(t, float)
    if law == "m1":
        return a0 * np.exp(-t / p["tau"])
    if law == "m3":
        n = p["n"]
        return (a0 ** (1 - n) + (n - 1) * p["c"] * t) ** (1 / (1 - n))
    if law == "m4":
        return a0 * ((1 - p["phi"]) * np.exp(-t / p["tau_f"]) + p["phi"] * np.exp(-t / p["tau_s"]))
    if law == "m5":
        e = np.exp(-p["a"] * t)
        return p["a"] * a0 * e / (p["a"] + p["b"] * a0 * (1 - e))
    if law == "m6":
        kr, ke, ac = p["k_run"], p["k_evap"], p["A_c"]
        kt = kr + ke
        ws_ = kr * ac / kt                                         # fixed point of the run-off phase (< A_c)
        a0 = np.asarray(a0, float)
        above = a0 > ac
        with np.errstate(divide="ignore", invalid="ignore"):
            t1 = np.where(above, np.log(np.maximum(a0 - ws_, 1e-12) / (ac - ws_)) / kt, 0.0)
        first = ws_ + (a0 - ws_) * np.exp(-kt * t)
        second = np.minimum(a0, ac) * np.exp(-ke * np.maximum(t - t1, 0))
        return np.where(above & (t < t1), first, second)
    raise ValueError(f"unknown law {law!r}")


def spread_z(skel, seed):
    """z_i = sqrt(1/2) z_link + sqrt(1/2) z_tail, fixed by the seed (the same z for every law and every s)."""
    rng = np.random.default_rng(seed)
    links = np.sort(skel.cml_id.unique())
    z_link = dict(zip(links, rng.standard_normal(len(links))))
    z_tail = rng.standard_normal(len(skel))
    return np.sqrt(0.5) * skel.cml_id.map(z_link).to_numpy() + np.sqrt(0.5) * z_tail


def truth_matrix(skel, law, p, f, horizon):
    """W / A0 for every tail (rows) and minute 0 .. horizon - 1; NaN after the tail's usable end. f: time-scale factors."""
    t = np.arange(int(horizon))[None, :]
    a0 = skel.A_start_trace_db.to_numpy()[:, None]
    w = law_curve(law, p, t / np.asarray(f)[:, None], a0) / a0
    return np.where(t < skel.usable_min.to_numpy()[:, None], w, np.nan)


def quantize(x, phase, step):
    return phase + step * np.round((x - phase) / step)


def ar1(e, sigma, tau, x_prev):
    """AR(1) with stationary std `sigma`, correlation time `tau` (steps), unit innovations e, previous value x_prev."""
    rho = np.exp(-1.0 / tau)
    return lfilter([np.sqrt(1 - rho ** 2) * sigma], [1, -rho], e, zi=[rho * x_prev])[0]


def simulate(skel, law, p, noise_by_sub, switches, seed, f, step=0.3, spm=6, nan_masks=None):
    """Simulate one tail set. Returns (rows, info).
    rows: one row per tail minute (ids, minutes_since_t0_trace, W_db = truth (mean of the 6 10 s values), A_db =
          observed 1 min A); info: one row per tail (f, offset, noise parameters, phase).
    skel: tail_id, sublink_id, cml_id, event_id, usable_min, A_start_trace_db, offset_std.
    noise_by_sub: {sublink_id: (sigma_white, sigma_slow, tau_noise)}, or with 5 entries (+ sigma_slow2, tau_noise2: a
    second AR(1) wander). f: time-scale factor per tail (ones = no spread).
    nan_masks: {tail_id: bool array of NaN minutes} (used if switches['real_nan_mask'])."""
    rng = np.random.default_rng(seed)
    rng2 = np.random.default_rng([seed, 2])                         # second slow component only (own stream)
    f = np.asarray(f) if switches.get("spread", True) else np.ones(len(skel))
    rows, info = [], []
    for i, r in enumerate(skel.itertuples()):
        L = int(r.usable_min)
        nz_ = noise_by_sub.get(r.sublink_id, (np.nan, 0.0, 1.0))
        sw, ss, tn = nz_[:3]
        t10 = np.arange(L * spm) / spm
        w10 = law_curve(law, p, t10 / f[i], r.A_start_trace_db)
        # draw everything first (same numbers whatever the switches are)
        off = rng.standard_normal() * (r.offset_std if np.isfinite(r.offset_std) else 0.0)
        eps = rng.standard_normal(L * spm)
        e_slow, x_prev = rng.standard_normal(L), rng.standard_normal() * ss
        phase = rng.uniform(0, step)
        slow = np.repeat(ar1(e_slow, ss, tn, x_prev), spm) if ss > 0 else np.zeros(L * spm)
        if len(nz_) > 3 and nz_[3] > 0:                             # second AR(1) wander (stationary start)
            slow = slow + np.repeat(ar1(rng2.standard_normal(L), nz_[3], nz_[4], rng2.standard_normal() * nz_[3]), spm)
        x = (w10 + (off if switches["offset"] else 0.0) + (slow if switches.get("slow_noise", True) else 0.0)
             + (sw * eps if switches["noise"] and np.isfinite(sw) else 0.0))
        q = quantize(x, phase, step) if switches["quantization"] else x
        a1 = q.reshape(L, spm).mean(axis=1)
        if switches.get("real_nan_mask") and nan_masks is not None and r.tail_id in nan_masks:
            a1 = np.where(nan_masks[r.tail_id][:L], np.nan, a1)
        rows.append(pd.DataFrame({"law": law, "tail_id": r.tail_id, "sublink_id": r.sublink_id, "cml_id": r.cml_id,
                                  "event_id": r.event_id, "minutes_since_t0_trace": np.arange(L, dtype="int16"),
                                  "W_db": w10.reshape(L, spm).mean(axis=1).astype("float32"),
                                  "A_db": a1.astype("float32")}))
        info.append({"law": law, "tail_id": r.tail_id, "f_time_scale": f[i], "offset_db": off,
                     "sigma_white_db": sw, "sigma_slow_db": ss, "tau_noise_min": tn, "phase_db": phase})
    return pd.concat(rows, ignore_index=True), pd.DataFrame(info)


def dry_stats(a, q=None, lags=()):
    """Targets of a dry series. a: 1 min A (blocks x minutes, NaN allowed); q: the quantized 10 s values
    (blocks x minutes x 6) or None. Returns within-block std, on-level share (all 6 samples equal), ACF at lags."""
    d = a - np.nanmean(a, axis=1, keepdims=True)
    n = np.isfinite(d).sum()
    ss = np.nansum(d ** 2)
    out = {"within_std_db": float(np.sqrt(ss / max(n - a.shape[0], 1)))}
    if q is not None:
        out["on_level_share"] = float(np.mean(np.all(q == q[:, :, :1], axis=2)))
    for k in lags:
        out[f"acf_{k}"] = float(np.nansum(d[:, :-k] * d[:, k:]) / ss) if a.shape[1] > k else np.nan
    return out


def dry_lookup(grid, blocks, block_min, step, spm, lags, seed):
    """Dry series (W = 0) for every grid point (sigma_white, sigma_slow, tau_noise): the targets of dry_stats.
    Common random numbers: the same unit white noise, unit AR(1) path (per tau) and phases at every grid point."""
    rng = np.random.default_rng(seed)
    eps = rng.standard_normal((blocks, block_min, spm))
    e_slow = rng.standard_normal((blocks, block_min))
    x_prev = rng.standard_normal(blocks)
    ph = rng.uniform(0, step, (blocks, 1, 1))
    out = []
    for tn in grid["tau_noise_min"]:
        u = np.stack([ar1(e_slow[b], 1.0, tn, x_prev[b]) for b in range(blocks)])[:, :, None]   # unit-std AR(1)
        for ss in grid["sigma_slow_db"]:
            for sw in grid["sigma_white_db"]:
                q = quantize(ss * u + sw * eps, ph, step)
                out.append({"sigma_white_db": sw, "sigma_slow_db": ss, "tau_noise_min": tn,
                            **dry_stats(q.mean(axis=2), q, lags)})
    return pd.DataFrame(out)
