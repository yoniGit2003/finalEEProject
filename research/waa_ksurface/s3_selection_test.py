"""Unit test of the drying selection of Deviation 10 (`s3_recovery.select_nonnested`): one toy case per outcome.

    python -m research.waa_ksurface.s3_selection_test           (also collected by pytest: the functions test_*)

A toy case is a table of held-out errors (dB) of the four rungs. The gain of a richer rung over a simpler one is the
difference of their errors; its interval is the gain -/+ a half-width. No data and no fit are used. Margin: 0.002 dB.

    1  neither K1 nor K2 beats K0                         -> K0, "no dependence", K3 not tried
    2  K1 clearly ahead of K2 (the numbers of truth (c))  -> K1, "W";   and the mirror case -> K2, "R"
    3  K1 and K2 both beat K0, closer than the margin     -> the better one, attribution not resolved
    4  K3 beats the selected rung                         -> K3, "W and R"

Further cases guard the edges of the rule: a gain above the margin whose interval holds zero does not count; only one of
K1, K2 beating K0 is "resolved" however close they are; K3 is compared with the selected rung, not with K0.
"""
import pandas as pd

from research.common.report import rel
from research.common.save import save
from research.waa_ksurface import paths as P
from research.waa_ksurface.s3_recovery import UNRESOLVED, select_nonnested

MARGIN = 0.002


def toy(errors, half_width=0.0005, wide=()):
    """compare(richer, simpler) for a table of errors; the pairs in `wide` get an interval that holds zero."""
    def compare(richer, simpler):
        g = errors[simpler] - errors[richer]
        h = abs(g) + 0.001 if (richer, simpler) in wide else half_width
        return g, g - h, g + h
    return compare


CASES = [  # name, errors of K0..K3 (dB), pairs with an interval that holds zero, expected rung, expected attribution, is K3 tried
    ("1 neither beats K0", {"K0": 0.2000, "K1": 0.1990, "K2": 0.1995, "K3": 0.1985}, (), "K0", "no dependence", False),
    ("2 K1 ahead of K2 (truth (c): gains 0.0135, 0.0059, K3 0.0136)", {"K0": 0.2134, "K1": 0.1999, "K2": 0.2075, "K3": 0.1998}, (), "K1", "W", True),
    ("2 mirror: K2 ahead of K1", {"K0": 0.2134, "K1": 0.2075, "K2": 0.1999, "K3": 0.1998}, (), "K2", "R", True),
    ("3 both beat K0, closer than the margin", {"K0": 0.2100, "K1": 0.2010, "K2": 0.2020, "K3": 0.2005}, (), "K1", UNRESOLVED, True),
    ("3 mirror: K2 slightly better", {"K0": 0.2100, "K1": 0.2020, "K2": 0.2010, "K3": 0.2005}, (), "K2", UNRESOLVED, True),
    ("4 K3 beats the selected rung", {"K0": 0.2100, "K1": 0.2000, "K2": 0.2050, "K3": 0.1950}, (), "K3", "W and R", True),
    ("4 K3 after an unresolved pair", {"K0": 0.2100, "K1": 0.2010, "K2": 0.2020, "K3": 0.1950}, (), "K3", "W and R", True),
    ("edge: K1 above the margin but its interval holds zero", {"K0": 0.2100, "K1": 0.2050, "K2": 0.2095, "K3": 0.2040}, (("K1", "K0"),), "K0", "no dependence", False),
    ("edge: only K2 beats K0 (K1's interval holds zero), close gains", {"K0": 0.2100, "K1": 0.2055, "K2": 0.2050, "K3": 0.2045}, (("K1", "K0"),), "K2", "R", True),
    ("edge: K3 far below K0 but not below the selected K1", {"K0": 0.2100, "K1": 0.1900, "K2": 0.2090, "K3": 0.1890}, (), "K1", "W", True),
]


def run(case):
    name, errors, wide, rung, attribution, k3_tried = case
    got_rung, got_attr, steps = select_nonnested(toy(errors, wide=wide), MARGIN)
    ok = (got_rung == rung) and (got_attr == attribution) and (("K3" in steps) == k3_tried)
    if "K3" in steps:                                           # K3 is compared with the rung selected before it, never with K0
        ok = ok and steps["K3"][0] in ("K1", "K2")
    ok = ok and sum(bool(s[4]) for s in steps.values()) == (0 if got_rung == "K0" else 1)      # one rung is marked as selected
    return {"case": name, "expected rung": rung, "returned rung": got_rung, "expected attribution": attribution, "returned attribution": got_attr,
            "K3 tried": "K3" in steps, "pass": bool(ok)}


def test_outcome_1_k0():
    assert all(run(c)["pass"] for c in CASES if c[0].startswith("1"))


def test_outcome_2_one_rung():
    assert all(run(c)["pass"] for c in CASES if c[0].startswith("2"))


def test_outcome_3_unresolved():
    assert all(run(c)["pass"] for c in CASES if c[0].startswith("3"))


def test_outcome_4_k3():
    assert all(run(c)["pass"] for c in CASES if c[0].startswith("4"))


def test_edges():
    assert all(run(c)["pass"] for c in CASES if c[0].startswith("edge"))


def main():
    t = pd.DataFrame([run(c) for c in CASES])
    save(t, P.S3_SELECTION_TESTS, index=False)
    print(t.to_string(index=False))
    print(f"{int(t['pass'].sum())} of {len(t)} cases pass; table: {rel(P.S3_SELECTION_TESTS)}")
    assert t["pass"].all(), "the selection of Deviation 10 does not do what is written"


if __name__ == "__main__":
    main()
