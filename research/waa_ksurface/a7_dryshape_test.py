"""Unit tests of the dry-path diagnostic (toy data only)."""
import numpy as np

from research.waa_ksurface import a7_dryshape as ds


def test_exponential_ratio_is_one():
    t = np.array([5, 10, 15, 20, 25, 30], float)
    Y = np.vstack([2.0 * np.exp(-t / 8.0), 0.5 * np.exp(-t / 20.0)])
    re, rl, r = ds.ratios(Y)
    assert np.allclose(re, [1 / 8, 1 / 20]) and np.allclose(r, 1.0)


def test_slowing_drying_ratio_below_one_and_nonpositive_dropped():
    t = np.array([5, 10, 15, 20, 25, 30], float)
    slow = 1.0 / (1.0 + 0.2 * t)                            # rate falls with time
    bad = slow.copy(); bad[3] = -0.1
    _, _, r = ds.ratios(np.vstack([slow, bad]))
    assert r[0] < 1 and np.isnan(r[1])


def test_boot_counts_copies():
    ev = np.array([0, 0, 1])
    lo, hi = ds.boot(ev, [lambda rows: float(len(rows))], 200, 1)
    assert lo[0] >= 2 and hi[0] <= 4


def test_control_offsets_off_keeps_random_stream():
    from research.waa_ksurface import a7_dryshape_control as cc, fake6
    keep = fake6.A6["noise"]
    assert keep["offset_sd_db"] > 0 and set(cc.CC["truths"]) == {"G", "R"}
    rng1, rng2 = np.random.default_rng(5), np.random.default_rng(5)
    a, b = rng1.normal(0, 0.136, 4), rng2.normal(0, 0.0, 4)        # sd 0 still draws: the later numbers stay the same
    assert np.allclose(b, 0) and rng1.random() == rng2.random()
