"""Unit tests of the WAA-only comparison (w_cv.py) on synthetic numbers.

    python -m pytest research/waa_ksurface/w_cv_test.py
"""
import numpy as np

from research.waa_ksurface import w_cv as W
from research.wet_antenna.waa_equation import WAAEquation


def test_two_store_equals_part1_equation():
    """Our numba two-store kernel = WAAEquation (exact formula) with the rain fed at its 1 min lag."""
    rng = np.random.default_rng(0)
    r = np.clip(rng.gamma(0.6, 6.0, 300) - 1.0, 0, None)
    eq = WAAEquation(exact=True)
    ref = eq.run(r, 1.0, "28 GHz", a_link=1.0, rain_before=np.array([0.0]))["W"]
    rl = np.r_[0.0, r[:-1]]
    A = np.array([eq.eq["A_band_db"]["28 GHz"]["value"], 1.0])
    out = W.two_store(rl, np.array([0]), np.array([len(r)]), np.array([0]), A, eq.b, np.exp(-1 / eq.tau_wet), np.exp(-1 / eq.tau_dry),
                      np.exp(-1 / eq.tau_s), eq.phi)
    assert np.allclose(out, ref)


def test_huber():
    r = np.array([0.0, 0.1, -0.5])
    assert np.allclose(W.huber(r), [0.0, 0.005, 0.3 * (0.5 - 0.15)])
