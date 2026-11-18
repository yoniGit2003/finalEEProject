"""Unit test of A7 (`a7_real_fit.py`, Deviation 15): A7 runs the step 3c code path exactly. No data and no fit are used.

    python -m research.waa_ksurface.a7_test        (also collected by pytest: the functions test_*)

    1  step 3c's run_one and A7's fit_part both go through `s3c_a7_variant.procedure` (one code path, not a copy)
    2  procedure calls ladder.two_stage with the arguments of step 3c: the use weights, the gauge set, the A6 extension
       starts, extension free, wetting K0, stage 2 = (counted & used, high & used, margin 0.00291, K0 vs K1)
    3  A7 hands procedure the real attenuation, the band index, the measured curves and the use weights unchanged, and
       reads k_dry at the gauge estimate (x) and at the radar estimate (x_alt); no held 1 / k_wet is left in A7
"""
import inspect

import numpy as np
import pandas as pd

from research.waa_ksurface import a6_recovery as a6, a7_real_fit as a7, ladder, s1_separability as s1, s3_recovery as s3, s3c_a7_variant as s3c


def toy_setting():
    le = pd.DataFrame({"event_id": ["E1", "E1", "E2"]})
    return {"D": {"tag": "real"}, "le": le, "bin_le": np.array([0, 0, 1, 2, 2]), "gauge": np.array([1.0, 0.0, 1.0]),
            "iF": [np.zeros(3), np.zeros(2)], "iP": [np.ones(2)], "learn": np.array([True, True, True])}


def test_one_code_path():
    run_src, fit_src = inspect.getsource(s3c.run_one), inspect.getsource(a7.fit_part)
    assert "procedure(S, Dk, bi, S[\"curves\"], use," in run_src and "two_stage" not in run_src
    assert "s3c.procedure(S, D, band_idx, curves, use," in fit_src
    a7_src = inspect.getsource(a7)
    assert "fit_fixed_kwet" not in a7_src and "tau_wet_min" not in a7_src and "TAUS" not in a7_src


def test_procedure_arguments(monkeypatch):
    S, seen = toy_setting(), {}
    counted = np.array([True, True, False, True, True])
    high = np.array([True, False, False, True, True])
    monkeypatch.setattr(s1, "score_sets", lambda S_, Dk: (counted, None, high, None))
    monkeypatch.setattr(ladder, "two_stage", lambda *a, **k: seen.update(args=a, kw=k) or {"ok": True})
    use = np.array([1.0, 1.0, 0.0])                             # E2 flagged
    res, (h, used) = s3c.procedure(S, S["D"], "BI", "CURVES", use, "(t) ")
    a, k = seen["args"], seen["kw"]
    assert a[0] is S["D"] and a[1] is S["le"] and a[2] == "BI" and a[3] == "CURVES" and a[4] is S["gauge"]
    assert k["use"] is use and k["hold_ext"] is False and k["wet_rungs"] == s3.WET_RUNGS
    E = a6.A6["extension"]
    assert all(np.allclose(f, E["f_start_inc_db"]) for f in k["inc_F"]) and all(np.allclose(p, E["p_start_share"]) for p in k["inc_P"])
    st2 = k["stage2"]
    exp_used = np.array([True, True, True, False, False])
    assert np.array_equal(used, exp_used)
    assert np.array_equal(st2["w5"], counted & exp_used) and np.array_equal(st2["score"], high & exp_used)
    assert abs(st2["margin_db"] - 0.00291) < 1e-5 and st2["selector"] is s3.select_k0_k1 and st2["rungs"] == s3.DRY_RUNGS


class ToyModel:
    def k(self, x, regime, R, W):
        return np.full(np.shape(W), 1.0 / x[0])

    def extension(self, x):
        return {("F", 0): np.array([x[0]])}


def test_a7_reads_both_estimates(monkeypatch):
    S, seen = toy_setting(), {}
    tab = pd.DataFrame({"rung": ["K0", "K1"], "cv": [1.0, 1.0]})

    def fake_procedure(S_, D, band_idx, curves, use, tag):
        seen.update(D=D, band_idx=band_idx, curves=curves, use=use)
        res = {"wet": {"table": tab}, "dry": {"rung": "K0", "lam": 0.0, "table": tab}, "model": ToyModel(),
               "x": np.array([15.0]), "x_alt": np.array([17.0]), "tau_wet": {"gauge": 3.0, "radar": 6.0}}
        return res, (np.ones(5, bool), np.ones(5, bool))

    monkeypatch.setattr(s3c, "procedure", fake_procedure)
    monkeypatch.setattr(s1, "score_sets", lambda S_, Dk: (np.ones(5, bool), None, np.ones(5, bool), None))
    monkeypatch.setattr(s3, "k1_gain", lambda t: (0.0, -1.0, 1.0))
    use = np.array([1.0, 0.0, 1.0])
    out = a7.fit_part(S, "CURVES", "BI", use)
    assert seen["D"] is S["D"] and seen["band_idx"] == "BI" and seen["curves"] == "CURVES" and seen["use"] is use
    assert np.isclose(out["tau_dry"]["gauge"][1], 15.0) and np.isclose(out["tau_dry"]["radar"][1], 17.0)
    assert out["tau_wet"] == {"gauge": 3.0, "radar": 6.0} and a7.reading(out)[1] == a7.K0_TEXT
    assert a7.ESTIMATES == ("gauge", "radar")


def flag_toy(y_rows, c):
    """Six short link-events of one event, 40 minutes each, rain ending at minute 10 (le_last = 1)."""
    n, k = 40, len(y_rows)
    D = {"le_boff1": np.arange(k) * n, "le_n": np.full(k, n), "le_last": np.ones(k, int)}
    S = {"D": D, "le": pd.DataFrame({"event_id": ["E1"] * k}), "short": np.ones(k, bool)}
    return S, {"y1": np.concatenate(y_rows), "le_c": np.asarray(c, float)}


def old_flag(S, Dk):
    F, D, rows = a6.A6["flag"], S["D"], []
    for e in range(len(S["le"])):
        o, n, end = D["le_boff1"][e], D["le_n"][e], (D["le_last"][e] + 1) * 5
        y = Dk["y1"][o:o + n]
        a0, a1 = y[end - 5:end + 5].mean(), y[end + 10:end + 20].mean()
        rows.append((a0 >= F["a0_min_db"], (a0 >= F["a0_min_db"]) and (a1 >= a0 - F["fall_db"])))
    q, wt = sum(r[0] for r in rows), sum(r[1] for r in rows)
    return q, wt, bool(q >= F["min_link_events"] and wt >= F["share"] * q)


def test_flag_fake_unchanged_and_real_reading():
    rng = np.random.default_rng(0)
    for _ in range(20):                                         # fake data: c = 0, no NaN -> the old result
        ys = [rng.uniform(0, 1.5) + rng.normal(0, 0.1, 40) for _ in range(6)]
        S, Dk = flag_toy(ys, np.zeros(6))
        t = a6.tail_flags(S, Dk).iloc[0]
        assert (int(t.qualifying), int(t.waiting_or_rising), bool(t.flagged)) == old_flag(S, Dk)
    ys = [np.full(40, 1.0) for _ in range(6)]                   # A = 1 dB flat (waits), c = 0.5 -> A - c = 0.5 < 0.6: none qualifies
    S, Dk = flag_toy(ys, np.full(6, 0.5))
    assert int(a6.tail_flags(S, Dk).iloc[0].qualifying) == 0
    ys[0][:] = np.nan                                           # a link-event with no valid minute is not counted
    ys[1][5:15] = np.nan
    S, Dk = flag_toy(ys, np.zeros(6))
    t = a6.tail_flags(S, Dk).iloc[0]
    assert int(t.link_events) == 4 and int(t.qualifying) == 4


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-q"]))
