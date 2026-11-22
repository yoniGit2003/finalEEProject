"""Unit tests of Block B (b_val.py) on synthetic data (no VAL, no TEST).

    python -m pytest research/waa_ksurface/b_val_test.py
"""
import json

import numpy as np
import pandas as pd

from research.waa_ksurface import b_val as B, fp, ksurf
from research.wet_antenna.system_id import PADM
from research.wet_antenna.waa_equation import WAAEquation


def tiny(n_le=3, n=60, seed=0):
    rng = np.random.default_rng(seed)
    le_n = np.full(n_le, n)
    roff = np.arange(n_le) * (n + 2 * PADM)
    rain = np.clip(rng.gamma(0.5, 8.0, n_le * (n + 2 * PADM)) - 1.0, 0, None)
    D = {"le_n": le_n, "le_roff": roff, "rain": rain, "le_L": np.array([0.5, 1.5, 3.0])[:n_le], "le_boff1": np.arange(n_le) * n}
    le = pd.DataFrame({"cml_id": ["10001", "10003", "99999"][:n_le], "band_name": ["28 GHz", "38 GHz", "28 GHz"][:n_le]})
    M = fp.measured()
    curves = [M[b] for b in M]
    band_idx = le.band_name.map({b: i for i, b in enumerate(M)}).to_numpy()
    return D, le, curves, band_idx


def direct(D, le, curves, band_idx, kind, inc):
    cfg = ksurf.config()
    nF, nP = len(cfg["f_knots_mmh"]), len(cfg["p_knots_mmh"])
    out = []
    k = 0
    incF = []
    if kind == "no memory":
        incF = [inc[k + b * nF:k + (b + 1) * nF] for b in range(len(curves))]
        k += nF * len(curves)
    incP = [inc[k + b * nP:k + (b + 1) * nP] for b in range(len(curves))]
    for e in range(len(le)):
        o, n = D["le_roff"][e], D["le_n"][e]
        R = D["rain"][o + PADM:o + PADM + n]
        b = band_idx[e]
        y = D["le_L"][e] * ksurf.full_curve(curves[b], "P", incP[b], R)
        if kind == "no memory":
            y = y + ksurf.full_curve(curves[b], "F", incF[b], R)
        out.append(y.reshape(-1, 5).mean(axis=1))
    return np.concatenate(out)


def test_recipe_columns_are_exact():
    D, le, curves, bi = tiny()
    rng = np.random.default_rng(3)
    for kind in B.RECIPES:
        base, C, names = B.recipe_columns(D, le, curves, bi, kind)
        inc = rng.uniform(0, 1, C.shape[1])
        assert np.allclose(base + C @ inc, direct(D, le, curves, bi, kind, inc), atol=1e-12)
        assert len(B.recipe_bounds(names)) == C.shape[1]


def test_huber_fit_recovers_increments():
    D, le, curves, bi = tiny(n=600)
    base, C, names = B.recipe_columns(D, le, curves, bi, "no memory")
    used = np.abs(C).sum(axis=0) > 0                            # only increments the rain reaches are identifiable
    true = np.where(used, 0.4, 0.0)
    y = base + C @ true
    x, loss, conv, _ = B.huber_fit(base, C, y, np.ones(len(y)), np.full(C.shape[1], 0.1), B.recipe_bounds(names), 0.3)
    assert loss < 1e-8 and np.allclose(base + C @ x, y, atol=1e-3)     # the output is recovered (some increments touch few bins)
    assert np.median(np.abs(x[used] - 0.4)) < 1e-3


def test_part1_predict_equals_equation():
    D, le, curves, bi = tiny(n_le=2)
    pr = B.part1_predict(D, le)
    eq = WAAEquation()
    a = eq.eq["a_link"]["per_link"]
    for e in range(2):
        o, n = D["le_roff"][e], D["le_n"][e]
        r = eq.run(D["rain"][o + PADM:o + PADM + n], D["le_L"][e], le.band_name[e], a_link=a[le.cml_id[e]], rain_before=D["rain"][o + PADM - 1:o + PADM])
        assert np.allclose(pr[e * n // 5:(e + 1) * n // 5], r["A_minus_c"].reshape(-1, 5).mean(axis=1))


def test_predictable_drops_unknown_link():
    _, le, _, _ = tiny()
    ok = B.predictable(le)
    assert not ok[2]


def test_check_freeze_detects_change(tmp_path, monkeypatch):
    f = tmp_path / "x.txt"
    f.write_bytes(b"a\nb\n")
    monkeypatch.setattr(B.P, "REPO", tmp_path)
    rec = tmp_path / "rec.json"
    rec.write_text(json.dumps({"sha256": {"x": {"path": "x.txt", "sha256": B.sha(f)}}}))
    assert B.check_freeze(rec) == []
    f.write_bytes(b"a\nc\n")
    assert len(B.check_freeze(rec)) == 1
    f.write_bytes(b"a\r\nb\r\n")                                # line endings alone do not count as a change
    assert B.check_freeze(rec) == []
