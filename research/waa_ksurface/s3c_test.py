"""Unit test of step 3c (`s3c_a7_variant.py`, Deviation 14). No data and no fit are used.

    python -m research.waa_ksurface.s3c_test        (also collected by pytest: the functions test_*)

    1  seeds: every run has its own seed, none of them used before (step 3b used 3600-3723)
    2  the truths: (g) at 0.088 / 0.15 dB with constant k; detectability k_dry(W) at the ratios 1.5, 2, 3, five each
    3  the margin is that of the current variant in S1 (0.00291 dB), not that of the extended one
    4  the verdicts and the STOP rule on toy judge rows
"""
import numpy as np
import pandas as pd

from research.waa_ksurface import s3b_recovery as s3b, s3c_a7_variant as s3c


def all_names():
    return list(s3c.SPREAD) + [n for a in s3c.C["detect_sizes"] for n in s3c.det_names(a)]


def test_seeds_fresh_and_unique():
    names = all_names()
    seeds = [s3c.spec_seed(n)[1] for n in names]
    assert len(names) == 17 and len(seeds) == len(set(seeds))
    old = [s3b.spec_seed(n)[1] for n in s3b.RUNS + [f"{t}_r{k}" for t in ("c", "e") for k in range(1, 5)]
           + [s3b.det_name(a, k) for a in (0.5, 0.25) for k in range(5)]]
    assert not (set(seeds) & set(old)) and min(seeds) == 3800 and max(seeds) == 3834


def test_truths():
    for n, sd in s3c.SPREAD.items():
        spec = s3c.spec_seed(n)[0]
        assert spec["kind"] == "const" and spec["link_sd_db"] == sd and abs(1 / spec["k_dry"] - 15) < 1e-6
    for a in (0.25, 0.5, 1.0):
        spec = s3c.spec_seed(s3c.det_names(a)[0])[0]
        k = lambda w: spec["k_dry0"] * (1 + w / spec["w_scale_db"])
        assert abs(k(2) / k(0) - (1 + 2 * a)) < 1e-9 and abs(1 / k(1) - 15) < 1e-6     # 15 min at W = 1 dB for every a
    assert sorted(1 + 2 * a for a in s3c.C["detect_sizes"]) == [1.5, 2.0, 3.0] and s3c.C["replicates"] == 5


def test_margin_current_variant():
    assert abs(s3c.margin() - 0.00291) < 1e-5


def toy(rung, k_ok):
    j = pd.DataFrame([{"regime": "wetting", "selected_rung": "K0", "k_ok": True},
                      {"regime": "drying", "selected_rung": rung, "k_ok": k_ok}])
    return {"judge": j}


def test_verdicts():
    assert s3c.verdict("g_088", toy("K0", True)) == "PASS"
    assert s3c.verdict("g_150", toy("K0", False)) == "FAIL"
    assert s3c.verdict("g_088", toy("K1", True)) == "FAIL"
    assert s3c.verdict("cw_0.5_r1", toy("K1", False)) == "K1 returned"
    assert s3c.STOP_RUN == "g_088" and np.isclose(s3c.SPREAD[s3c.STOP_RUN], 0.088)


if __name__ == "__main__":
    for f in (test_seeds_fresh_and_unique, test_truths, test_margin_current_variant, test_verdicts):
        f()
        print("ok", f.__name__)
