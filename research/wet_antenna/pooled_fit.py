"""Pooled fit of a drying law to a whole dataset of tails (library, step 3b-1b).

Model for tail i (t = minutes since t0_trace, 1 min A):
    A_i(t) = W(t / f_i; theta, W0_i) + c_i + e_i(t)
    theta  shared law parameters of the candidate (below); W0_i per-tail start (free, profiled)
    f_i    per-tail speed factor, f_i = exp(s (u_link + v_tail)), u, v ~ N(0, 1/2) (half per link, half per tail);
           s is fitted; u and v are integrated out (marginalised)
    c_i    per-tail offset with the prior c_i ~ N(0, offset_std of the sub-link ^ 2) (MAP, profiled in closed form)
    e_i    noise: white + AR(1) at 1 min with the sub-link's parameters from the control tails (sim_noise_params);
           white + AR(1) is exactly an ARMA(1,1) process, so the residuals are whitened with that filter
           (Box-Jenkins conditional likelihood; the first minute gets its stationary variance). NaN minutes are
           dropped and the remaining minutes treated as contiguous (approximation; v1 tails have <= 10 % NaN).
Marginalisation of the speed factors: for each tail, its profile log-likelihood l_i(g) (g = u + v) is computed on a
lattice of g, the peak is refined and a LAPLACE (Gaussian) approximation in g is used, which makes the integral over
the per-tail part v analytic; the per-link part u is integrated by quadrature on a lattice of N(0, 1/2).

Candidates (theta):  F1 exponential (tau) | F4 two layers (tau_f, tau_s, phi) | F3 nonlinear dW/dt = -k W^n
(k, n in [1, 3]) | F5 linear + quadratic (a, b) | F6 state switch (k_run, k_evap, A_c), as in simulator.py.
FC (step 3b-3) is not a menu law but a FREE curve: dW/dt = -W k(W), with k(W) given as a table (free_curve.py builds
the table from a spline and calls pooled_loglik("FC", table, s, D)). Table layout in theta: [n, x_0, dx, T_0 .. T_(n-1),
k_0 .. k_(n-1)], x = log W on a regular grid, T(x) = integral from x_0 to x of dx / k (the time the law needs to dry
from x to x_0). A tail that starts at x(0) is then x(tau) = T^-1(T(x(0)) - tau). Outside the grid k is constant.

Two options added in step 4b (both OFF by default, so every earlier step gives the same numbers as before):
    bound_w0    the start of F1 and F4 is kept inside [0.001, 30] dB, as it always was for F3, F5, F6 and FC. Without
                it a tail that RISES gets a negative start (negative water) under F1 / F4, which the other laws cannot
                do, so their BIC values were not comparable on real tails (found in step 4).
    noise_mult  one number m > 0 that multiplies the noise VARIANCE of every tail of the dataset (the dry-weather
                noise model keeps its shape; m is its level during drying). The likelihood then needs the term
                -1/2 n log m per tail (it was a constant before). fit(..., fit_noise=True) fits m with theta and s.
"""
import numpy as np
from numba import njit, prange

CAND_CODE = {"F1": 1, "F3": 3, "F4": 4, "F5": 5, "F6": 6, "FC": 7}
LINEAR = {"F1", "F4"}                                               # W is proportional to W0: W0 profiled exactly
G_LATTICE = np.linspace(-4.0, 4.0, 17)                              # g = u + v, in units of z (sd 1)
U_NODES = np.linspace(-2.1, 2.1, 13)                                # u ~ N(0, 1/2): +-3 sd
U_LOGW = -U_NODES ** 2 / (2 * 0.5)
U_LOGW = U_LOGW - np.log(np.sum(np.exp(U_LOGW)))


def arma_params(sigma_white_10s, sigma_slow, tau_noise, step=0.3, spm=6):
    """1 min noise = AR(1) (sigma_slow, rho = exp(-1/tau)) + white (10 s white and quantization averaged over 6
    samples: (sigma_white^2 + step^2 / 12) / 6). Returns (rho, theta_ma, var_innovation, var_stationary)."""
    se2 = (sigma_white_10s ** 2 + step ** 2 / 12) / spm
    ss2 = sigma_slow ** 2
    rho = np.exp(-1.0 / tau_noise)
    a = ss2 * (1 - rho ** 2) + se2 * (1 + rho ** 2)
    b = se2 * rho
    if b > 1e-14:
        th = (a - np.sqrt(max(a * a - 4 * b * b, 0.0))) / (2 * b)
        s2 = b / th
    else:
        th, s2 = 0.0, a
    return rho, th, s2, ss2 + se2


@njit(cache=True)
def _curve(code, th, t, w0, f, out):
    """W(t / f) for one tail into `out`."""
    n = t.shape[0]
    if code == 7:                                                   # free curve: table look-up (see module docstring)
        m = int(th[0])
        x0, dx = th[1], th[2]
        T = th[3:3 + m]
        klo, khi, tmax, xmax = th[3 + m], th[3 + 2 * m - 1], th[3 + m - 1], th[1] + (m - 1) * th[2]
        xw = np.log(w0)
        j = 0
        if xw <= x0:
            ts = (xw - x0) / klo
        elif xw >= xmax:
            ts = tmax + (xw - xmax) / khi
            j = m - 2
        else:
            u = (xw - x0) / dx
            j = min(int(u), m - 2)
            ts = T[j] + (u - j) * (T[j + 1] - T[j])
        for k in range(n):
            tt = ts - t[k] / f
            if tt <= 0.0:
                out[k] = np.exp(x0 + klo * tt)
            elif tt >= tmax:
                out[k] = np.exp(xmax + khi * (tt - tmax))
            else:
                while j > 0 and T[j] > tt:                          # time runs forward, so the index only walks down
                    j -= 1
                while j < m - 2 and T[j + 1] < tt:
                    j += 1
                out[k] = np.exp(x0 + (j + (tt - T[j]) / (T[j + 1] - T[j])) * dx)
        return
    p3 = w0 ** (1 - th[1]) if code == 3 and abs(th[1] - 1.0) >= 1e-6 else 0.0   # F3: constant part, once per tail
    for k in range(n):
        x = t[k] / f
        if code == 1:
            out[k] = w0 * np.exp(-x / th[0])
        elif code == 4:
            out[k] = w0 * ((1 - th[2]) * np.exp(-x / th[0]) + th[2] * np.exp(-x / th[1]))
        elif code == 3:
            kk, nn = th[0], th[1]
            if abs(nn - 1.0) < 1e-6:
                out[k] = w0 * np.exp(-kk * x)
            else:
                out[k] = (p3 + (nn - 1) * kk * x) ** (1 / (1 - nn))
        elif code == 5:
            a, b = th[0], th[1]
            e = np.exp(-a * x)
            out[k] = a * w0 * e / (a + b * w0 * (1 - e))
        else:                                                       # 6: state switch
            kr, ke, ac = th[0], th[1], th[2]
            kt = kr + ke
            ws = kr * ac / kt
            if w0 > ac:
                t1 = np.log(max(w0 - ws, 1e-12) / (ac - ws)) / kt
                if x < t1:
                    out[k] = ws + (w0 - ws) * np.exp(-kt * x)
                else:
                    out[k] = ac * np.exp(-ke * (x - t1))
            else:
                out[k] = w0 * np.exp(-ke * x)


@njit(cache=True)
def _rate(code, th, w):
    """dW/dt of the base law (no speed factor) at state w (used for dW/dW0 = rate(W(t)) / rate(W0))."""
    if code == 3:
        return -th[0] * w ** th[1]
    if code == 5:
        return -th[0] * w - th[1] * w * w
    if code == 7:                                                   # free curve: -W k(W), k from the table
        m = int(th[0])
        u = (np.log(w) - th[1]) / th[2]
        if u <= 0.0:
            return -w * th[3 + m]
        if u >= m - 1:
            return -w * th[3 + 2 * m - 1]
        j = int(u)
        return -w * (th[3 + m + j] + (u - j) * (th[3 + m + j + 1] - th[3 + m + j]))
    if w > th[2]:                                                   # 6: state switch
        return -th[0] * (w - th[2]) - th[1] * w
    return -th[1] * w


@njit(cache=True)
def _whiten(x, rho, ma, out):
    """ARMA(1,1) innovations: eta_t = x_t - rho x_(t-1) + ma eta_(t-1); eta_0 = x_0."""
    out[0] = x[0]
    for k in range(1, x.shape[0]):
        out[k] = x[k] - rho * x[k - 1] + ma * out[k - 1]


@njit(cache=True)
def _profile(code, linear, th, t, fy, f1, rho, ma, s2, g0, o2, a0, f, bound=False, nm=1.0):
    """Profile log-likelihood of one tail at speed factor f: max over W0 and over the offset c (with its prior).
    fy, f1: the whitened data and the whitened constant 1 (they do not depend on the parameters).
    bound: keep the start of a linear law inside [0.001, 30] dB; nm: multiplier of the noise variance."""
    n = t.shape[0]
    w = np.empty(n)
    fw = np.empty(n)
    wt = np.full(n, 1.0 / (nm * s2))
    wt[0] = 1.0 / (nm * g0)
    lognm = 0.5 * n * np.log(nm)                                    # normalisation of the noise level (0 for nm = 1)
    s11 = np.sum(wt * f1 * f1) + 1.0 / o2
    if linear:
        _curve(code, th, t, 1.0, f, w)
        _whiten(w, rho, ma, fw)
        s12 = np.sum(wt * f1 * fw)
        s22 = np.sum(wt * fw * fw)
        r1 = np.sum(wt * f1 * fy)
        r2 = np.sum(wt * fw * fy)
        det = s11 * s22 - s12 * s12
        if det <= 1e-14:
            return -1e18, 0.0, 0.0
        c = (s22 * r1 - s12 * r2) / det
        w0 = (s11 * r2 - s12 * r1) / det
        if bound and (w0 < 1e-3 or w0 > 30.0):                      # at the bound: the best offset for that start
            w0 = min(max(w0, 1e-3), 30.0)
            c = (r1 - w0 * s12) / s11
        q = 0.0
        for k in range(n):
            e = fy[k] - c * f1[k] - w0 * fw[k]
            q += wt[k] * e * e
        return -0.5 * (q + c * c / o2) - lognm, w0, c
    # nonlinear in W0: two Gauss-Newton steps from W0 = A_start (sensitivity rate(W) / rate(W0)), then c
    w0 = max(a0, 0.05)
    wd = np.empty(n)
    fd = np.empty(n)
    for it in range(2):
        _curve(code, th, t, w0, f, w)
        r0 = _rate(code, th, w0)
        if r0 == 0.0:
            break
        for k in range(n):
            wd[k] = _rate(code, th, w[k]) / r0
        _whiten(w, rho, ma, fw)
        _whiten(wd, rho, ma, fd)
        s12 = np.sum(wt * f1 * fd)
        s22 = np.sum(wt * fd * fd)
        r1 = np.sum(wt * f1 * (fy - fw))
        r2 = np.sum(wt * fd * (fy - fw))
        det = s11 * s22 - s12 * s12
        if det <= 1e-14:
            break
        w0 = min(max(w0 + (s11 * r2 - s12 * r1) / det, 1e-3), 30.0)
    _curve(code, th, t, w0, f, w)
    _whiten(w, rho, ma, fw)
    c = np.sum(wt * f1 * (fy - fw)) / s11
    q = 0.0
    for k in range(n):
        e = fy[k] - fw[k] - c * f1[k]
        q += wt[k] * e * e
    return -0.5 * (q + c * c / o2) - lognm, w0, c


@njit(parallel=True, cache=True)
def _tails_laplace(code, linear, th, s, T, FY, F1, start, length, rho, ma, s2, g0, o2, a0, glat, bound=False, nm=1.0):
    """Per tail: Laplace summary (l_star, g_star, kappa) of l_i(g) over g = u + v."""
    nt = start.shape[0]
    lstar = np.empty(nt)
    gstar = np.empty(nt)
    kap = np.empty(nt)
    ng = glat.shape[0]
    for i in prange(nt):
        a, b = start[i], start[i] + length[i]
        t, fy, f1 = T[a:b], FY[a:b], F1[a:b]
        ll = np.empty(ng)
        for j in range(ng):
            ll[j] = _profile(code, linear, th, t, fy, f1, rho[i], ma[i], s2[i], g0[i], o2[i], a0[i], np.exp(s * glat[j]), bound, nm)[0]
        jm = np.argmax(ll)
        j0 = min(max(jm, 1), ng - 2)
        y0, y1, y2 = ll[j0 - 1], ll[j0], ll[j0 + 1]
        dg = glat[1] - glat[0]
        den = y0 - 2 * y1 + y2
        g = glat[j0] if den >= 0 else glat[j0] + 0.5 * dg * (y0 - y2) / den
        g = min(max(g, glat[0]), glat[ng - 1])
        hh = 0.05
        lm = _profile(code, linear, th, t, fy, f1, rho[i], ma[i], s2[i], g0[i], o2[i], a0[i], np.exp(s * (g - hh)), bound, nm)[0]
        lc = _profile(code, linear, th, t, fy, f1, rho[i], ma[i], s2[i], g0[i], o2[i], a0[i], np.exp(s * g), bound, nm)[0]
        lp = _profile(code, linear, th, t, fy, f1, rho[i], ma[i], s2[i], g0[i], o2[i], a0[i], np.exp(s * (g + hh)), bound, nm)[0]
        if lc < ll[jm]:                                             # refinement worse than the lattice: keep lattice
            g, lc = glat[jm], ll[jm]
            k2 = -den / (dg * dg)
        else:
            k2 = -(lp - 2 * lc + lm) / (hh * hh)
        lstar[i], gstar[i], kap[i] = lc, g, max(k2, 1e-3)
    return lstar, gstar, kap


def pooled_loglik(cand, theta, s, D, bound_w0=False, noise_mult=1.0):
    """Marginal log-likelihood of a dataset D (dict from pack()) for candidate `cand`, parameters theta, spread s.
    bound_w0, noise_mult: the two options of step 4b (module docstring); the defaults are the earlier behaviour."""
    code = CAND_CODE[cand]
    lstar, gstar, kap = _tails_laplace(code, cand in LINEAR, np.asarray(theta, float), float(s), D["T"], D["FY"], D["F1"],
                                       D["start"], D["length"], D["rho"], D["ma"], D["s2"], D["g0"], D["o2"],
                                       D["a0"], G_LATTICE, bool(bound_w0), float(noise_mult))
    # per tail and link node u: log h_i(u) = l* - 1/2 log(1 + k/2) - (k/2)(u - g*)^2 / (1 + k/2)
    k2 = kap[:, None] / 2
    logh = lstar[:, None] - 0.5 * np.log1p(k2) - k2 * (U_NODES[None, :] - gstar[:, None]) ** 2 / (1 + k2)
    per_link = np.zeros((D["n_links"], U_NODES.size))
    np.add.at(per_link, D["link"], logh)
    x = per_link + U_LOGW[None, :]
    m = x.max(axis=1, keepdims=True)
    return float(np.sum(m[:, 0] + np.log(np.exp(x - m).sum(axis=1))))


def pack(tails, noise, offset_sd, window=None, min_points=10):
    """Arrays for pooled_loglik. tails: list of dicts with tail_id, sublink_id, cml_id, A_start, t (minutes), a (1 min
    A, NaN allowed). noise: {sublink_id: (sigma_white, sigma_slow, tau_noise)}; offset_sd: {sublink_id: sd (dB)}.
    window: keep t < window minutes (None = whole tail). Links are numbered 0 .. n_links - 1."""
    T, Y, start, length, rho, ma, s2, g0, o2, a0, link, ids = [], [], [], [], [], [], [], [], [], [], [], []
    link_ix, pos = {}, 0
    for tl in tails:
        t, a = np.asarray(tl["t"], float), np.asarray(tl["a"], float)
        keep = np.isfinite(a) & ((t < window) if window else True)
        if keep.sum() < min_points:
            continue
        r_, m_, s_, g_ = arma_params(*noise[tl["sublink_id"]])
        T.append(t[keep]), Y.append(a[keep])
        start.append(pos), length.append(int(keep.sum()))
        pos += int(keep.sum())
        rho.append(r_), ma.append(m_), s2.append(s_), g0.append(g_)
        sd = offset_sd.get(tl["sublink_id"], np.nan)
        o2.append(max(sd if np.isfinite(sd) else 0.1, 0.01) ** 2)
        a0.append(max(float(tl["A_start"]), 0.05))
        link.append(link_ix.setdefault(tl["cml_id"], len(link_ix)))
        ids.append(tl["tail_id"])
    arr = lambda v, d=float: np.asarray(v, d)                         # noqa: E731
    FY, F1 = [], []
    for y_, r_, m_ in zip(Y, rho, ma):                              # whitened data and constant: fixed per dataset
        fy, f1 = np.empty(len(y_)), np.empty(len(y_))
        _whiten(y_, r_, m_, fy)
        _whiten(np.ones(len(y_)), r_, m_, f1)
        FY.append(fy), F1.append(f1)
    return {"T": np.concatenate(T), "Y": np.concatenate(Y), "FY": np.concatenate(FY), "F1": np.concatenate(F1),
            "start": arr(start, np.int64),
            "length": arr(length, np.int64), "rho": arr(rho), "ma": arr(ma), "s2": arr(s2), "g0": arr(g0),
            "o2": arr(o2), "a0": arr(a0), "link": arr(link, np.int64), "n_links": len(link_ix), "tail_id": ids}


# ------------------------------------------------------------------------------------------------- parameters
def to_theta(cand, x):
    """Unconstrained vector x -> (theta, s). Last element of x is log s."""
    s = float(np.exp(x[-1]))
    if cand == "F1":
        return [np.exp(x[0])], s
    if cand == "F4":
        tf = np.exp(x[0])
        return [tf, tf + np.exp(x[1]), 1 / (1 + np.exp(-x[2]))], s
    if cand == "F3":
        return [np.exp(x[0]), 1 + 2 / (1 + np.exp(-x[1]))], s
    if cand == "F5":
        return [np.exp(x[0]), np.exp(x[1])], s
    return [np.exp(x[0]), np.exp(x[1]), np.exp(x[2])], s


def from_theta(cand, theta, s):
    """(theta, s) -> unconstrained vector (inverse of to_theta)."""
    lg = lambda p: np.log(p / (1 - p))                                # noqa: E731
    if cand == "F1":
        x = [np.log(theta[0])]
    elif cand == "F4":
        x = [np.log(theta[0]), np.log(theta[1] - theta[0]), lg(theta[2])]
    elif cand == "F3":
        x = [np.log(theta[0]), lg((theta[1] - 1) / 2)]
    elif cand == "F5":
        x = [np.log(theta[0]), np.log(theta[1])]
    else:
        x = [np.log(theta[0]), np.log(theta[1]), np.log(theta[2])]
    return np.asarray(x + [np.log(s)], float)


NAMES = {"F1": ["tau"], "F4": ["tau_f", "tau_s", "phi"], "F3": ["k", "n"], "F5": ["a", "b"],
         "F6": ["k_run", "k_evap", "A_c"]}


def fit(cand, D, theta0, s0, maxfev=400, step=0.4, fatol=0.1, bound_w0=False, fit_noise=False):
    """Maximise the pooled marginal log-likelihood over (theta, s) with Nelder-Mead in the unconstrained space,
    starting from an explicit simplex (`step` in every coordinate; the default simplex is too small at log s = 0).
    fit_noise: also fit the noise-variance multiplier (one more coordinate, log m, started at 0)."""
    from scipy.optimize import minimize
    x0 = from_theta(cand, theta0, s0)
    if fit_noise:
        x0 = np.r_[x0, 0.0]
    simplex = np.vstack([x0] + [x0 + step * np.eye(len(x0))[k] for k in range(len(x0))])

    def cost(x):
        if fit_noise:
            ll = pooled_loglik(cand, *to_theta(cand, x[:-1]), D, bound_w0, float(np.exp(x[-1])))
        else:
            ll = pooled_loglik(cand, *to_theta(cand, x), D, bound_w0)
        return -ll if np.isfinite(ll) else 1e18                     # a start far from the optimum can overflow

    res = minimize(cost, x0, method="Nelder-Mead",
                   options={"maxfev": maxfev, "initial_simplex": simplex, "xatol": 1e-3, "fatol": fatol})   # fatol: log-likelihood units
    theta, s = to_theta(cand, res.x[:-1] if fit_noise else res.x)
    out = {"ll": -float(res.fun), "theta": [float(v) for v in theta], "s": s, "nfev": int(res.nfev),
           "success": bool(res.success)}
    if fit_noise:
        out["noise_mult"] = float(np.exp(res.x[-1]))
    return out


def tail_speed_factors(cand, theta, s, D):
    """Estimated speed factor f_i of every tail of D at the fitted (theta, s): {tail_id: f}. Per tail the profile
    likelihood in g is summarised by its peak g* and curvature kappa (as in pooled_loglik); with the prior g ~ N(0, 1)
    the posterior mode is kappa g* / (1 + kappa), so a tail that says little about its speed is pulled to f = 1.
    Used by step 3b-2 to rescale the time axis of each tail (t / f_i)."""
    _, gstar, kap = _tails_laplace(CAND_CODE[cand], cand in LINEAR, np.asarray(theta, float), float(s), D["T"], D["FY"],
                                   D["F1"], D["start"], D["length"], D["rho"], D["ma"], D["s2"], D["g0"], D["o2"],
                                   D["a0"], G_LATTICE)
    return dict(zip(D["tail_id"], np.exp(s * kap * gstar / (1 + kap))))


def f1_start(cand, tau1):
    """A start for `cand` that is close to the F1 (exponential) solution with time constant tau1 (minutes)."""
    return {"F3": [1 / tau1, 1.05],                        # n near 1 = exponential with rate k = 1 / tau
            "F4": [tau1, 5 * tau1, 0.05],                  # nearly all weight in one layer
            "F5": [1 / tau1, 0.1 / tau1],                  # small quadratic term
            "F6": [0.1 / tau1, 1 / tau1, 1.0]}[cand]       # weak run-off, evaporation rate 1 / tau
