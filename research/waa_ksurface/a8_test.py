"""Unit tests of step 8 (`a8_forms.py`). Toy data only (no real attenuation, no real fit).

    python -m research.waa_ksurface.a8_test        (also collected by pytest: the functions test_*)

    1  the tabulated form is within 0.002 in ln k of each physics form between the grid points (W 0.01-9 dB)
    2  the FormModel loss equals the plain ksurf.Model loss with the same knots (one kernel, not a copy)
    3  the FormModel gradient matches central finite differences
    4  the AICc rule: the fewest parameters within 2 of the best wins; the Pareto front and the trial order
    5  a PySR expression becomes a form with its constants as parameters; k1_lnk reads Model.k
    6  a curve fit recovers the parameters of a noise-free form
"""
import types

import numpy as np
import pandas as pd

from research.waa_ksurface import a8_forms as F8, ksurf
from research.waa_ksurface.a5_build import toy_dataset

LIKE_A7 = {"linear": [np.log(0.02), 0.05], "sat": [np.log(0.2), np.log(0.5)], "power": [np.log(0.1), 1.0],
           "sat_floor": [np.log(0.01), np.log(0.2), np.log(0.5)]}


def toy_prob():
    rng = np.random.default_rng(8)
    rains = []
    for k in range(6):                                          # rain, then a long dry tail
        r = np.zeros(180)
        r[10:40 + 10 * k] = rng.uniform(0.5, 8.0, 30 + 10 * k)
        rains.append(r)
    D, rain = toy_dataset(rains, [0.4, 0.8, 1.5, 3.0, 0.6, 2.0])
    cur = types.SimpleNamespace(F=lambda r: 3.0 * r / (r + 2.0), P=lambda r: 0.1 * r, F_max=3.2, P_max=1.4)
    prob = ksurf.prepare(D, np.zeros(6, int), [cur], rain, 0)
    prob["target"] = rng.normal(0.5, 0.4, len(prob["target"]))
    return prob


def test_table_accuracy():
    Wf = np.geomspace(0.01, 9.0, 3000)
    fm_shell = types.SimpleNamespace(m=types.SimpleNamespace(lnk=(np.log(1 / 600), np.log(3.0))))
    for name, f in F8.physics().items():
        th = LIKE_A7[name]
        lo, hi = fm_shell.m.lnk
        tab = np.clip(f.lnk(F8.GRID_EVAL, th), lo, hi)
        i = np.clip(np.searchsorted(F8.GRID, Wf) - 1, 0, len(F8.GRID) - 2)
        inside = (tab[i] > lo) & (tab[i] < hi) & (tab[i + 1] > lo) & (tab[i + 1] < hi)   # the kernel's own clip makes a kink: not a table error
        err = np.abs(np.interp(Wf, F8.GRID, tab) - f.lnk(Wf, th))[inside].max()
        assert inside.mean() > 0.9 and err < 0.002, (name, err, inside.mean())


def _pair(name="sat_floor"):
    prob = toy_prob()
    f = F8.physics()[name]
    fm = F8.FormModel(prob, f, (np.array([np.log(1 / 3.0)]), np.zeros(0)))
    shell = F8.k1_model(1)
    x_ext = np.zeros(shell.n)
    x_ext[shell.sl["F"]] = 0.2
    x_ext[shell.sl["P"]] = 0.3
    return prob, fm, fm.start(LIKE_A7[name], x_ext)


def test_same_loss_as_ksurf():
    prob, fm, z = _pair()
    plain = ksurf.Model(prob, "K0", "K1", 0.0, cfg=fm.cfg)
    a = fm.loss(z)[0]
    b = plain.loss(fm.xfull(z), penalty=False)[0]
    assert abs(a - b) <= 1e-12 * abs(b)
    assert len(plain.wk["drying"]) == len(F8.GRID)


def test_gradient():
    for name in ("sat_floor", "power"):
        prob, fm, z = _pair(name)
        _, g = fm.loss(z)
        fd = np.zeros_like(z)
        for j in range(len(z)):
            h = 1e-5 * max(1.0, abs(z[j]))
            a, b = z.copy(), z.copy()
            a[j] += h
            b[j] -= h
            fd[j] = (fm.loss(a)[0] - fm.loss(b)[0]) / (2 * h)
        assert np.abs(g - fd).max() <= 1e-4 * np.abs(g).max(), (name, np.abs(g - fd).max(), np.abs(g).max())


def test_rank_rule():
    t = pd.DataFrame({"name": ["A", "B", "C", "D"], "p": [3, 2, 1, 4],
                      "aicc_gauge": [0.0, 1.5, 9.0, -0.5], "aicc_radar": [0.0, 1.0, 12.0, 0.4]})
    r, order = F8.rank(t)
    # Delta gauge: A 0.5, B 2, C 9.5, D 0; Delta radar: A 0, B 1, C 12, D 0.4 -> Delta_max A 0.5, B 2, C 12, D 0.4
    # best D; within 2 of 0.4: A, B, D -> fewest parameters: B. Front (nothing with <= p and a lower Delta_max): all four.
    r = r.set_index("name")
    assert np.allclose(r.dmax[["A", "B", "C", "D"]], [0.5, 2.0, 12.0, 0.4])
    assert r.index[r.chosen].tolist() == ["B"] and r.index[r.best].tolist() == ["D"]
    assert r.front.all()
    assert [t.name[i] for i in order] == ["B", "D", "A"]        # then by ascending Delta_max; C has Delta_max 12 > 10: not tried
    assert r.loc["C", "trial"] == -1
    t2 = t.copy()
    t2.loc[0, "aicc_gauge"] = 3.0                               # A now has Delta_max 3.5 and is dominated by B (p 2, 2.0)
    assert not F8.rank(t2)[0].set_index("name").loc["A", "front"]


def test_from_sympy_and_k1():
    f = F8.from_sympy("pysr_x", "0.2*W/(W + 0.5)")
    assert f.p == 2 and f.source == "PySR"
    W = np.array([0.1, 1.0, 5.0])
    assert np.allclose(f.k(W, f.starts[0]), 0.2 * W / (W + 0.5))
    m = F8.k1_model(2)
    x = np.zeros(m.n)
    x[m.sl["gd"]] = np.log(1 / 60.0)
    x[m.sl["hd"]] = [1.0, 2.0, 1.5, 1.0, 0.5]
    assert np.allclose(F8.k1_lnk(x, W), np.log(m.k(x, "drying", 0.0, W)))


def test_curve_fit_recovers():
    Wp, w, _ = F8.points()
    for name, f in F8.physics().items():
        y = f.lnk(Wp, LIKE_A7[name])
        th, wrss = F8.fit_curve(f, Wp, y, w)
        assert wrss < 1e-10, (name, wrss)
        assert np.allclose(f.k(Wp, th), f.k(Wp, LIKE_A7[name]), rtol=1e-4), name


def test_pysr_front_and_structure():
    from research.waa_ksurface import a8_name as A8
    t = pd.DataFrame({"estimate": ["gauge"] * 5 + ["radar"] * 2, "seed": [0, 0, 1, 1, 0, 0, 1], "complexity": [1, 3, 3, 5, 7, 1, 3],
                      "loss": [0.9, 0.5, 0.4, 0.45, 0.2, 0.8, 0.85], "equation": [""] * 7, "sympy": [""] * 7})
    f = A8.pysr_front(t)
    # gauge: c1 0.9, c3 best 0.4 (seed 1), c5 0.45 > 0.4 not kept, c7 0.2 kept; radar: c1 0.8 kept, c3 0.85 not kept
    assert f.index[f.front].tolist() == [0, 2, 4, 5]
    assert A8.structure("0.2*W/(W + 0.5)") == A8.structure("0.3*W/(W + 1.5)") != A8.structure("0.2*W**0.5")


if __name__ == "__main__":
    for k, v in list(globals().items()):
        if k.startswith("test_"):
            v()
            print(f"{k}: ok")
