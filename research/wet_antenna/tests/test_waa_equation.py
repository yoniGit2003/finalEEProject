"""Unit test of research/wet_antenna/waa_equation.py: it reproduces the S8b model's 1 min W and path term of 3 TEST link-events
(outputs/products/waa_equation_reference.json, written by formulate_s9.py from the fit's own kernel) within 1e-6 dB.
Run from the repository root:  python -m pytest research/wet_antenna/tests
"""
import json

import numpy as np

from research.wet_antenna import paths as P
from research.wet_antenna.waa_equation import WAAEquation


def test_reproduces_the_fit():
    ref = json.loads(P.S9_REFERENCE.read_text(encoding="utf-8"))
    eq = WAAEquation()
    assert len(ref["link_events"]) == 3
    for le in ref["link_events"]:
        out = eq.run(le["rain_mm_h"], le["L_km"], le["band"], le["a_link"], rain_before=le["rain_before"])
        assert np.max(np.abs(out["W"] - np.asarray(le["W"]))) <= ref["tolerance_db"]
        assert np.max(np.abs(out["path"] - np.asarray(le["path"]))) <= ref["tolerance_db"]
