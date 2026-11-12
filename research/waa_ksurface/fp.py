"""Library: the MEASURED curves F(R) (steady WAA) and P(R) (path loss per km) of step E4, as functions.

F comes from the F-valid cells and P from the P-valid cells of `outputs/tables/e4_remeasure_cells.csv` (Deviation 4).
Each curve is the monotone curve of step E4: through (0, 0) and the cell values, made non-decreasing first.

    m = measured()            {band: Curves}
    m["28 GHz"].F(r)          the measured F; above its valid range F is HELD (an extrapolation, not a measurement)
    m["28 GHz"].P(r)          the measured P; above its valid range it continues with the slope of its last segment (E4 rule)
    m["28 GHz"].P_prop(r)     the same P, but above its valid range the line through the origin and the last valid cell
    m["28 GHz"].F_max, .P_max the upper end of the valid range of each quantity (mm/h)

`from_values` builds the same object from any cell values: the fake worlds of step A6 use it with a known truth.
"""
from dataclasses import dataclass
from typing import Callable

import numpy as np
import pandas as pd

from research.waa_ksurface import ltest
from research.waa_ksurface import paths as P


@dataclass
class Curves:
    band: str
    F: Callable
    P: Callable
    P_prop: Callable
    F_r: np.ndarray            # R of the F-valid cells (mm/h)
    F_v: np.ndarray            # their values after the monotone step (dB)
    P_r: np.ndarray
    P_v: np.ndarray            # (dB/km)

    @property
    def F_min(self):
        return float(self.F_r[0])

    @property
    def F_max(self):
        return float(self.F_r[-1])

    @property
    def P_min(self):
        return float(self.P_r[0])

    @property
    def P_max(self):
        return float(self.P_r[-1])


def from_values(band, f_r, f_v, p_r, p_v, f_w=None, p_w=None):
    """The curves through the given cells. f_w, p_w = weights of the monotone step (default: equal)."""
    f_r, f_v, p_r, p_v = (np.asarray(x, float) for x in (f_r, f_v, p_r, p_v))
    F, fv = ltest.monotone_curve(f_r, f_v, np.ones(len(f_r)) if f_w is None else f_w, "F")
    Pc, pv = ltest.monotone_curve(p_r, p_v, np.ones(len(p_r)) if p_w is None else p_w, "P")
    p_top, r_top = float(pv[-1]), float(p_r[-1])

    def P_prop(r):
        r = np.asarray(r, float)
        return np.where(r <= r_top, Pc(np.minimum(r, r_top)), p_top * r / r_top)
    return Curves(band, F, Pc, P_prop, f_r, fv, p_r, pv)


def measured():
    """{band: Curves} from the E4 cells (the weights of the monotone step are those of E4: 1 / bootstrap variance)."""
    c = pd.read_csv(P.E4_CELLS)
    c = c[c.eligible]
    out = {}
    for band, g in c.groupby("band"):
        f = g[g.F_valid].sort_values("R_cell")
        p = g[g.P_valid].sort_values("R_cell")
        out[band] = from_values(band, f.R_cell, f.I, p.R_cell, p.s, 1.0 / np.maximum(f.I_se.to_numpy(), 1e-6) ** 2, 1.0 / np.maximum(p.s_se.to_numpy(), 1e-6) ** 2)
    return out
