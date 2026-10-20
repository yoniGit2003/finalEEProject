"""The WAA equation of Part 1 (closed form adopted in S8b): rain series, link length, band and a_link in; W(t) and the path term out.

    A(t) - c = W_f(t) + W_s(t) + L c_band R(t - d)
    dW_f/dt = (phi u - W_f) / tau_f,  tau_f = tau_f_wet while W_f rises, tau_f_dry while it falls
    dW_s/dt = ((1 - phi) u - W_s) / tau_s
    u = a_link A_band (R(t - d) / 5)^b

R is the radar path-mean rain divided by the event's radar / gauge ratio g_e (mm/h, one value per minute), d the lag (1 min).
Every number comes from outputs/products/waa_equation_final.json (written by formulate_s9.py). The stores are updated once per
minute exactly as in the fit (system_id.kernel2, mode 2): the decay factors exp(-1 / tau) per minute, the output after the
update. The fit evaluated F on 61 knots (0 and 60 log-spaced knots from 0.05 to 40 mm/h), linear between them and constant
above 40 mm/h; `exact=False` (default) uses that grid and reproduces the fit, `exact=True` the pure formula. Valid 0-40 mm/h.

Example:
    from research.wet_antenna.waa_equation import WAAEquation
    eq = WAAEquation()
    out = eq.run(rain_mm_h, L_km=1.2, band="28 GHz", a_link=1.0)     # out["W"], out["path"], out["A_minus_c"] per minute
"""
import json

import numpy as np

from research.wet_antenna import paths as P


def load(path=None):
    """The equation file (outputs/products/waa_equation_final.json)."""
    return json.loads((path or P.S9_EQUATION).read_text(encoding="utf-8"))


class WAAEquation:
    def __init__(self, eq=None, exact=False):
        self.eq = eq or load()
        self.exact = exact
        self.knots = np.asarray(self.eq["tabulated_knots_mmh"], float)
        dyn = self.eq["dynamics_fixed"]
        self.phi, self.tau_wet, self.tau_dry, self.tau_s = dyn["phi"], dyn["tau_f_wet"], dyn["tau_f_dry"], dyn["tau_s"]
        self.lag = int(self.eq["lag_min"])
        self.b = self.eq["b"]["value"]

    def band_of(self, frequency_ghz):
        """The band name of a sub-link frequency (the bands of the equation file)."""
        for name, (lo, hi) in self.eq["bands_ghz"].items():
            if lo <= frequency_ghz < hi:
                return name
        raise ValueError(f"{frequency_ghz} GHz is outside the bands of the equation")

    def F(self, R, band):
        """The antenna level (dB) in steady rain R for a link with a_link = 1."""
        A = self.eq["A_band_db"][band]["value"]
        R = np.asarray(R, float)
        f = lambda r: A * np.power(np.maximum(r, 0.0) / 5.0, self.b)                       # noqa: E731
        if self.exact:
            return np.where(R > 0, f(R), 0.0)
        return np.interp(R, self.knots, np.r_[0.0, f(self.knots[1:])])

    def P(self, R, band):
        """The path loss per km (dB/km) in rain R."""
        c = self.eq["c_band_db_km_per_mmh"][band]["value"]
        R = np.asarray(R, float)
        if self.exact:
            return c * np.maximum(R, 0.0)
        return np.interp(R, self.knots, c * self.knots)

    def run(self, rain_mm_h, L_km, band, a_link=1.0, c=0.0, rain_before=None):
        """The equation along a rain series (one value per minute). Minute t responds to rain[t - lag]; for t < lag the rain
        before the series is `rain_before` (lag values, oldest first) or 0. The stores are empty before minute 0.
        Returns a dict of arrays: W_fast, W_slow, W, path (= L P), A_minus_c (= W + path), A (= A_minus_c + c)."""
        r = np.asarray(rain_mm_h, float)
        lead = np.zeros(self.lag) if rain_before is None else np.asarray(rain_before, float)
        rl = np.r_[lead, r][:len(r)]                                                      # rain[t - lag]
        u = a_link * self.F(rl, band)
        path = L_km * self.P(rl, band)
        kfw, kfd, ks = np.exp(-1.0 / self.tau_wet), np.exp(-1.0 / self.tau_dry), np.exp(-1.0 / self.tau_s)
        wf, ws = 0.0, 0.0
        Wf, Ws = np.empty(len(r)), np.empty(len(r))
        for t in range(len(r)):
            uf, us = self.phi * u[t], (1.0 - self.phi) * u[t]
            wf = uf + (wf - uf) * (kfw if uf > wf else kfd)
            ws = us + (ws - us) * ks
            Wf[t], Ws[t] = wf, ws
        W = Wf + Ws
        return {"W_fast": Wf, "W_slow": Ws, "W": W, "path": path, "A_minus_c": W + path, "A": W + path + c}
