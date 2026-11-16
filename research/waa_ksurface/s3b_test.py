"""Unit test of step 3b (`s3b_recovery.py`, Deviations 12 and 13). No data and no fit are used.

    python -m research.waa_ksurface.s3b_test        (also collected by pytest: the functions test_*)

    1  seeds: every run has its own seed, none of them used before (1200+, 2200+, 3200-3214, 3300-3305, 3400-3434)
    2  ratio 3 of the detectability study is truth (c) itself (same k_dry(W))
    3  the verdicts of Deviation 12 on toy judge rows, and the STOP rules
    4  the B3' rule: separable only if every (e) value is above the largest (c) value
    5  the fit returns every rung's best fit (B3' reads the K1 fit), and the drying ladder is K0 against K1 only
"""
import inspect

import numpy as np
import pandas as pd

from research.waa_ksurface import fake6, ladder, s3b_recovery as s3b


def all_names():
    det = [s3b.det_name(a, k) for a in s3b.B["detect_sizes"] if a != 1.0 for k in range(s3b.B["b3_replicates"])]
    return s3b.RUNS + [n for t in ("c", "e") for n in s3b.rep_names(t)[1:]] + det


def test_seeds_fresh_and_unique():
    seeds = [s3b.spec_seed(n)[1] for n in all_names()]
    assert len(seeds) == len(set(seeds))
    used = set(range(1200, 1300)) | set(range(2200, 2500)) | set(range(3200, 3215)) | set(range(3300, 3306)) | set(range(3400, 3435))
    assert not (set(seeds) & used)
    assert s3b.spec_seed("a")[1] == 3600 and s3b.spec_seed("g")[1] == 3606 and s3b.spec_seed("f_added")[1] == 3607
    assert [s3b.spec_seed(n)[1] for n in s3b.rep_names("c")[1:]] == [3620, 3621, 3622, 3623]
    assert [s3b.spec_seed(n)[1] for n in s3b.rep_names("e")[1:]] == [3630, 3631, 3632, 3633]


def test_ratio3_is_truth_c():
    c = s3b.spec_seed("c")[0]
    cw = s3b.a6.spec_of("cw_1_r0")[0]
    W = np.linspace(0, 3, 31)
    assert np.allclose(fake6.k_true(c, "drying", 0.0, W), fake6.k_true(cw, "drying", 0.0, W), rtol=1e-8)
    assert s3b.truth_of("c_r3") == "c" and s3b.truth_of("cw_0.5_r1") == "cw" and s3b.truth_of("d_high") == "d_high"


def toy(rung, true, k_ok, share=0.9, share_alt=0.9):
    j = pd.DataFrame([{"regime": "wetting", "selected_rung": "K0", "true_rung": "K0", "k_ok": True, "share_within_25pct": np.nan, "share_alt": np.nan},
                      {"regime": "drying", "selected_rung": rung, "true_rung": true, "k_ok": k_ok, "share_within_25pct": share, "share_alt": share_alt}])
    return {"judge": j}


def test_verdicts_and_stops():
    v = s3b.verdict
    assert v("a", toy("K0", "K0", True)) == "PASS" and v("a", toy("K1", "K0", True)) == "FAIL"
    assert v("c", toy("K1", "K1", True)) == "PASS" and v("c", toy("K0", "K1", False)) == "FAIL"
    assert v("c_r2", toy("K1", "K1", True)) == "PASS"
    assert v("d", toy("K0", "K2", False, 0.85, 0.82)) == "PASS" and v("d", toy("K0", "K2", False, 0.85, 0.7)) == "FAIL"
    assert v("d_high", toy("K0", "K0", True)) == "PASS"
    assert v("e", toy("K1", "none", None)) == "report"
    assert v("f", toy("K1", "K0", True)) == "spurious k_dry" and v("f", toy("K0", "K0", True)) == "no spurious k_dry"
    assert v("g", toy("K0", "K0", False)) == "k corrupted"
    caught, missed = {"flagged_with_unseen_water": 7, "unseen_water_events": 7}, {"flagged_with_unseen_water": 5, "unseen_water_events": 7}
    assert v("f_added", {**toy("K0", "K0", True), **caught}) == "PASS"
    assert v("f_added", {**toy("K0", "K0", True), **missed}) == "FAIL"
    assert v("f_added", {**toy("K1", "K0", True), **caught}) == "FAIL, material"
    assert s3b.stops("a", "FAIL") and s3b.stops("f", "spurious k_dry") and not s3b.stops("e", "report") and not s3b.stops("c_r1", "FAIL")


def test_b3_prime_rule():
    assert s3b.separable([0.0, 0.01, -0.01, 0.02, 0.005], [0.03, 0.04, 0.021, 0.05, 0.06]) == (True, 0.02)
    assert s3b.separable([0.0, 0.01, -0.01, 0.02, 0.005], [0.03, 0.04, 0.02, 0.05, 0.06])[0] is False       # equal is not above
    assert s3b.separable([0.0, 0.01, -0.01, 0.02, 0.005], [0.03, np.nan, 0.03, 0.05, 0.06])[0] is False     # missing counts as not above
    assert s3b.separable([0.0, 0.01], [0.03, 0.04])[0] is False                                             # fewer than 5 replicates


def test_k1_fit_returned_and_ladder():
    assert '"fits": {r: (best[r]["model"], best[r]["x"]) for r in best}' in inspect.getsource(ladder.stage)
    assert s3b.s3.DRY_RUNGS == ["K0", "K1"] and s3b.s3.WET_RUNGS == ["K0"]
    src = inspect.getsource(s3b.run_one)
    assert '"score": high & used' in src and "use=use" in src and 'res["dry"]["fits"]["K1"]' in src


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
            print(f"{name}: pass")
