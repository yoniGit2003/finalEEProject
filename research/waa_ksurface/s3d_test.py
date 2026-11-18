"""Unit test of step 3d (`s3d_kwet_fixed.py`). No data and no fit are used.

    python -m research.waa_ksurface.s3d_test        (also collected by pytest: the functions test_*)

    1  seeds 3840-3843, fresh (step 3c used 3800-3834); truths (a) and (c); the three values of 1 / k_wet
    2  the verdicts: a rung flip is a STOP, (a) K0 outside 25 % is a FAIL, everything ok is a PASS
    3  the fit holds the wetting rate: the stage-2 call passes the held values (other_x) to ladder.stage
"""
import inspect

import numpy as np
import pandas as pd

from research.waa_ksurface import s3d_kwet_fixed as s3d


def test_seeds_truths_values():
    seeds = {n: s3d.spec_seed(n)[1] for n in s3d.names()}
    assert seeds == {"a_r0": 3840, "a_r1": 3841, "c_r0": 3842, "c_r1": 3843}
    assert s3d.spec_seed("a_r0")[0]["kind"] == "const" and s3d.spec_seed("c_r1")[0]["kind"] == "dry_w"
    assert s3d.D3["tau_wet_min"] == [2.5, 5.0, 10.0]


def rows(a_rungs, c_rungs, a_ok=True):
    r = [{"truth": "a", "drying rung": g, "k_ok": a_ok} for g in a_rungs] + [{"truth": "c", "drying rung": g, "k_ok": True} for g in c_rungs]
    return pd.DataFrame(r)


def test_verdicts():
    assert s3d.verdict(rows(["K0"] * 3, ["K1"] * 3)) == "PASS"
    assert s3d.verdict(rows(["K0", "K1", "K0"], ["K1"] * 3)) == "STOP (rung flip)"
    assert s3d.verdict(rows(["K0"] * 3, ["K1", "K0", "K1"])) == "STOP (rung flip)"
    assert s3d.verdict(rows(["K0"] * 3, ["K1"] * 3, a_ok=False)) == "FAIL"


def test_wetting_held():
    src = inspect.getsource(s3d.fit_fixed_kwet)
    assert "np.log(1.0 / tau)" in src and 'st1["rung"], wet,' in src
    assert np.isclose(np.exp(np.log(1.0 / 2.5)), 0.4)


if __name__ == "__main__":
    for f in (test_seeds_truths_values, test_verdicts, test_wetting_held):
        f()
        print("ok", f.__name__)
