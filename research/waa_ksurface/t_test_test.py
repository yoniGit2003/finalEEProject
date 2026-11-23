"""Unit tests of the TEST script (t_test.py) on synthetic data (no TEST read).

    python -m pytest research/waa_ksurface/t_test_test.py
"""
import numpy as np

from research.waa_ksurface import b_val as B, t_test as T
from research.waa_ksurface.b_val_test import tiny


def test_new_links_use_a_link_one():
    D, le, _, _ = tiny(n_le=3)
    le = le.assign(cml_id=["10001", "10003", "10005"])
    assert np.allclose(T.part1_new(D, le), B.part1_predict(D, le.assign(cml_id=["x1", "x2", "x3"]), default_a=1.0))

