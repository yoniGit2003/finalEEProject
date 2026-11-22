"""Unit tests of the offset check (a8o_offset.py). Synthetic numbers only.

    python -m pytest research/waa_ksurface/a8o_test.py
"""
import numpy as np

from research.waa_ksurface import a8o_offset as O


def grid(y, phase, q=0.3):
    return q * np.round((y + phase) / q) - phase


def test_offset_after_grid_equals_before_in_law():
    """Q_phi(y + c) = Q_{phi + c}(y) + c exactly, and the grid is periodic in the phase: adding c after the rounding is
    rounding y + c with the phase phi - c, which is again uniform."""
    rng = np.random.default_rng(1)
    y, c = rng.normal(0, 1, 1000), 0.137
    for phi in rng.uniform(0, 0.3, 20):
        assert np.allclose(grid(y + c, phi), grid(y, phi + c) + c)
        assert np.allclose(grid(y, phi), grid(y, phi + 0.3))


def test_add_offsets_per_link_event():
    Dk = {"y1": np.zeros(20)}
    out = O.add_offsets(Dk, np.array([10, 10]), np.array([0.1, -0.2]))
    assert np.allclose(out["y1"][:10], 0.1) and np.allclose(out["y1"][10:], -0.2)
    assert np.allclose(out["y5"], [0.1, 0.1, -0.2, -0.2])


def test_w_half():
    assert np.isclose(O.w_half([-18.13, -0.004, -0.0006]), np.log(0.004) / -18.13)
    assert np.isclose(O.w_half([-18.0, 0.004, 0.0006]), np.log(0.004) / -18.0)
    assert np.isnan(O.w_half([-18.0, -0.004, 0.0006]))


def test_reading_rules():
    class R(dict):
        pass
    import research.waa_ksurface.s3b_recovery as s3b
    import pandas as pd

    def run(rung, wh):
        return {"judge": pd.DataFrame([{"regime": "drying", "selected_rung": rung}]), "switch": {"w_half_db": wh}}
    k1 = {f"r{i}": run("K1", 0.3) for i in range(3)}
    assert O.reading({**k1, "r3": run("K0", 0.1), "r4": run("K0", 0.1)})[0].startswith("consistent")
    assert O.reading({**{f"r{i}": run("K1", 0.6) for i in range(3)}, "r3": run("K0", 0.1), "r4": run("K0", 0.1)})[0] == "undecided"
    assert O.reading({f"r{i}": run("K0", 0.3) for i in range(4)} | {"r4": run("K1", 0.3)})[0] == "the non-exponential claim stands"
