"""Library of step S6b: the two-store model M4b with a drying time of the fast store that may follow the wind.

    dW_f/dt = (phi u - W_f) / tau_f,   tau_f = tau_f_wet while W_f rises,
                                       tau_f = tau_f_dry * exp(-k_wind * u_w(t)) while it falls
    (the slow store and everything else exactly as M4b; u_w = the hourly SMHI wind, standardised, held over the hour)

k_wind = 0 gives M4b exactly (the kernel is checked against system_id.kernel2 on the S6 candidate). The gradient is exact
(the same backward pass as kernel2, with the per-minute drying factor). ProblemW has the interface of system_id.Problem2
plus one parameter "k_wind"; it is used through final_model.Curves(..., Pm=ProblemW(...), free=[...]).
Settings: config `wet_antenna.return_s6b`.
"""
import numpy as np
import pandas as pd
from numba import njit

from research.wet_antenna import system_id as sid

PADM = sid.PADM


@njit(cache=True)
def kernel_w(d, Finc, Pinc, a, tau, phi, kw, uw, wle, le_roff, le_n, le_boff, le_link, le_band, le_L, le_c, rain, y, w, agg, delta, knots,
             pred, le_loss, gF, gP, ga, gtau, gphi, gk, want_grad):
    """M4b (mode 2, no drift, one slow time constant) with the drying factor of the fast store per minute. gtau[b, 1] is
    dL / d tau_f_dry (as kernel2), gk[0] = dL / d k_wind. uw has the layout of rain (index le_roff + PADM + t, no lag)."""
    nb, nk = Pinc.shape
    cF, cP = np.zeros((nb, nk)), np.zeros((nb, nk))
    for b in range(nb):
        for j in range(1, nk):
            cF[b, j] = cF[b, j - 1] + Finc[b, j - 1]
            cP[b, j] = cP[b, j - 1] + Pinc[b, j - 1]
    maxn = le_n.max()
    Wf, Ws, fr, kdt, tet = np.empty(maxn), np.empty(maxn), np.empty(maxn), np.empty(maxn), np.empty(maxn)
    jj, br = np.empty(maxn, np.int64), np.empty(maxn, np.int64)
    psi = np.empty(maxn)
    G0F, G1F, G0P, G1P = np.zeros((nb, nk)), np.zeros((nb, nk)), np.zeros((nb, nk)), np.zeros((nb, nk))
    total = 0.0
    for e in range(len(le_n)):
        n, b, link, L, c = le_n[e], le_band[e], le_link[e], le_L[e], le_c[e]
        al, ph = a[link], phi[b]
        ro, uo, bo = le_roff[e] + PADM - d, le_roff[e] + PADM, le_boff[e]
        kfw, ks = np.exp(-1.0 / tau[b, 0]), np.exp(-1.0 / tau[b, 2])
        wf, ws, acc = 0.0, 0.0, 0.0
        for t in range(n):
            r = rain[ro + t]
            j, f_ = nk - 1, 1.0
            if r <= 0.0:
                j, f_ = 0, 0.0
            else:
                for q in range(nk):
                    if r < knots[q + 1]:
                        j, f_ = q, (r - knots[q]) / (knots[q + 1] - knots[q])
                        break
            jj[t], fr[t] = j, f_
            u = al * (cF[b, j] + f_ * Finc[b, j])
            uf, us = ph * u, (1.0 - ph) * u
            Wf[t], Ws[t] = wf, ws
            te = tau[b, 1] * np.exp(-kw * uw[uo + t])
            kd = np.exp(-1.0 / te)
            kdt[t], tet[t] = kd, te
            if uf > wf:
                br[t] = 0
                wf = uf + (wf - uf) * kfw
            else:
                br[t] = 1
                wf = uf + (wf - uf) * kd
            ws = us + (ws - us) * ks
            acc += wf + ws + L * (cP[b, j] + f_ * Pinc[b, j]) + c
            if (t + 1) % agg == 0:
                pred[bo + t // agg] = acc / agg
                acc = 0.0
        loss = 0.0
        for q in range(n // agg):
            psi[q] = 0.0
            if w[bo + q] > 0:
                res = pred[bo + q] - y[bo + q]
                if abs(res) <= delta:
                    loss += 0.5 * res * res
                    psi[q] = res
                else:
                    loss += delta * (abs(res) - 0.5 * delta)
                    psi[q] = delta if res > 0 else -delta
        le_loss[e] = loss
        total += wle[e] * loss
        if want_grad and wle[e] != 0.0:
            lf, ls, gkfw, gks, glt, gkk = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
            for t in range(n - 1, -1, -1):
                gp = wle[e] * psi[t // agg] / agg
                j, f_ = jj[t], fr[t]
                G0P[b, j] += gp * L
                G1P[b, j] += gp * L * f_
                fval = cF[b, j] + f_ * Finc[b, j]
                u = al * fval
                uf, us = ph * u, (1.0 - ph) * u
                lf += gp
                ls += gp
                if br[t] == 0:
                    k = kfw
                    gkfw += lf * (Wf[t] - uf)
                else:
                    k = kdt[t]
                    gl = lf * (Wf[t] - uf)                 # dL / d(drying factor of minute t)
                    glt += gl * k / tet[t]                 # ... -> d / d log tau_f_dry
                    gkk -= gl * k * uw[uo + t] / tet[t]    # ... -> d / d k_wind
                guf = lf * (1.0 - k)
                lf *= k
                gus = ls * (1.0 - ks)
                gks += ls * (Ws[t] - us)
                ls *= ks
                gu = ph * guf + (1.0 - ph) * gus
                gphi[b] += u * (guf - gus)
                ga[link] += gu * fval
                G0F[b, j] += gu * al
                G1F[b, j] += gu * al * f_
            gtau[b, 0] += gkfw * kfw / tau[b, 0] ** 2
            gtau[b, 1] += glt / tau[b, 1]
            gtau[b, 2] += gks * ks / tau[b, 2] ** 2
            gk[0] += gkk
    if want_grad:
        for b in range(nb):
            sF, sP = 0.0, 0.0
            for j in range(nk - 1, -1, -1):
                gF[b, j] = sF + G1F[b, j]
                gP[b, j] = sP + G1P[b, j]
                sF += G0F[b, j]
                sP += G0P[b, j]
    return total


class ProblemW(sid.Problem2):
    """M4b (linear increments) with the extra parameter k_wind (bounds k_bounds); uw = the standardised wind per minute in the
    layout of D["rain"]."""

    def __init__(self, D, M, uw, k_bounds, agg=5):
        super().__init__(D, "M4b", M, agg)
        assert self.lin and self.mode == 2 and not self.drift
        self.names = np.append(self.names, "k_wind")
        self.bounds = list(self.bounds) + [tuple(k_bounds)]
        self.uw = np.ascontiguousarray(np.asarray(uw, float))
        assert len(self.uw) == len(D["rain"])

    def run(self, x, d, wle, want_grad=True):
        D = self.D
        Finc, Pinc, a, tau, phi, _, _, _, _ = self.unpack(x)
        kw = float(x[self.names == "k_wind"][0])
        gF, gP, ga = np.zeros_like(Finc), np.zeros_like(Pinc), np.zeros(self.nl)
        gtau, gphi, gk = np.zeros((self.nb, 4)), np.zeros(self.nb), np.zeros(1)
        loss = kernel_w(int(d), Finc, Pinc, a, tau, phi, kw, self.uw, wle, D["le_roff"], D["le_n"], self.boff, D["le_link"], D["le_band"], D["le_L"],
                        D["le_c"], D["rain"], self.y, self.w, self.agg, self.M["huber_delta_db"], self.knots, self.pred, self.le_loss, gF, gP, ga, gtau,
                        gphi, gk, want_grad)
        g = np.zeros(len(x))
        la = x[self.names == "loga"]
        loss += self.pen * float(la @ la)
        if not want_grad:
            return loss, g
        g[self.names == self.pname] = gP.ravel()
        g[self.names == self.fname] = gF.ravel()
        g[self.names == "loga"] = ga * a + 2 * self.pen * la
        for i, nm in enumerate(("tau_f_wet", "tau_f_dry", "tau_s")):
            g[self.names == nm] = (gtau[:, i] * tau[:, i]).sum()
        g[self.names == "phi"] = (gphi * phi * (1 - phi)).sum()
        g[self.names == "k_wind"] = gk[0]
        return loss, g


def wind_hourly(path):
    """The hourly SMHI GbgA wind (m/s) by UTC hour (baseline product dew_risk_hourly.csv, made from the release met_1h)."""
    t = pd.read_csv(path, parse_dates=["time_utc"]).set_index("time_utc")
    return t.wind_speed.astype(float)


def wind_minutes(D, le, wind, mean, sd):
    """The standardised wind (wind - mean) / sd per minute in the layout of D["rain"] (window start - PADM ... end + PADM),
    the value of an hour held over [hour, hour + 1 h); a missing hour = 0 (the mean)."""
    out = np.zeros(len(D["rain"]))
    for e, r in enumerate(le.itertuples()):
        n, o = D["le_n"][e], D["le_roff"][e]
        t0 = pd.Timestamp(r.window_start) - pd.Timedelta(minutes=PADM)
        tm = t0 + pd.to_timedelta(np.arange(n + 2 * PADM), "min")
        v = wind.reindex(tm.floor("h")).to_numpy(float)
        out[o:o + n + 2 * PADM] = np.where(np.isfinite(v), (v - mean) / sd, 0.0)
    return out
