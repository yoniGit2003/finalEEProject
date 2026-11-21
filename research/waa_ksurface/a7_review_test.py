"""Unit test of the A7 review (`a7_review.py`): the helpers that are new. No data and no fit are used; the stage-2 function
itself is checked on the real data by its first part (sel_mid must reproduce the committed A7 stage table, else STOP).

    python -m research.waa_ksurface.a7_review_test        (also collected by pytest: the functions test_*)
"""
import numpy as np

from research.waa_ksurface import a7_review as rv, ladder


def test_unit_counts_change_nothing():
    S, N = np.array([0.0, 2.0, 3.0, 0.0, 5.0]), np.array([0.0, 4.0, 6.0, 0.0, 8.0])
    cnt = np.array([0.0, 1.0, 1.0, 1.0, 1.0])                  # event 3 used but scored nothing (N = 0)
    Se, Ne = rv.expand(S, N, cnt)
    assert ladder.cv_db(Se, Ne) == ladder.cv_db(S, N)
    a = ladder.paired_interval(S, S * 0.9, N, 200, 1, [5, 95])
    b = ladder.paired_interval(Se, Se * 0.9, Ne, 200, 1, [5, 95])
    assert np.allclose(a, b)


def test_drawn_events_count_as_copies():
    ev = np.array([0, 0, 1, 2])
    w = np.array([2.0, 2.0, 0.0, 1.0])
    cnt = rv.event_counts(ev, w, 3)
    assert list(cnt) == [2.0, 0.0, 1.0]
    Se, Ne = rv.expand(np.array([1.0, 9.0, 3.0]), np.array([2.0, 9.0, 4.0]), cnt)
    assert list(Se) == [1.0, 1.0, 3.0] and list(Ne) == [2.0, 2.0, 4.0]


def toy_dry():
    """One link-event, 60 min: rain 0-14 min (1 mm/h), dry after; le_last = 2 (bins 0, 1, 2 wet)."""
    n = 60
    rain1 = np.where(np.arange(n) < 15, 1.0, 0.0)
    y5 = np.arange(n // 5, dtype=float)
    D = {"le_n": np.array([n]), "le_boff1": np.array([0]), "le_boff5": np.array([0]), "le_last": np.array([2]), "le_c": np.array([0.5]),
         "y5": y5, "w5": np.ones(n // 5)}
    return D, rain1


def test_dry_bins():
    D, rain1 = toy_dry()
    idx, target = rv.dry_values(D, rain1, [5, 10, 20, 30], 10, 0.1, np.array([1.0]))
    assert list(idx[0]) == [3, 4, 6, 8]                         # the bins that end 5, 10, 20, 30 min after minute 15
    assert target[3] == 2.5
    idx, _ = rv.dry_values(D, rain1, [5], 16, 0.1, np.array([1.0]))
    assert idx[0, 0] == -1                                      # only 15 wet minutes before: not used
    r2 = rain1.copy()
    r2[27:30] = 1.0                                             # bin 5 wet again: the 20 and 30 min values are not used
    idx, _ = rv.dry_values(D, r2, [5, 10, 20, 30], 10, 0.1, np.array([1.0]))
    assert list(idx[0]) == [3, 4, -1, -1]
    D["w5"][4] = 0.0
    idx, _ = rv.dry_values(D, rain1, [10], 10, 0.1, np.array([1.0]))
    assert idx[0, 0] == -1                                      # not counted
    idx, _ = rv.dry_values(D, rain1, [5], 10, 0.1, np.array([0.0]))
    assert idx[0, 0] == -1                                      # not used


def test_parts_order():
    p = rv.parts()
    assert p[:4] == ["sel_mid", "dry", "sel_flo", "sel_fhi"] and len(p) == 4 + 2 * 5


if __name__ == "__main__":
    for name, f in list(globals().items()):
        if name.startswith("test_"):
            f()
    print("ok")
