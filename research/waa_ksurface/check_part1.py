r"""Step 0 of the k-surface branch: confirm the Part 1 files the branch reuses exist and load. Reads only, writes nothing.

Run from the repository root:  C:\Users\yonic\miniconda3\envs\cmlpipe\python.exe research/waa_ksurface/check_part1.py
"""
import importlib
import json
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]       # the repository root
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "pipeline" / "src"))

from research.baseline import paths as BP  # noqa: E402
from research.wet_antenna import paths as WP  # noqa: E402

TAB = REPO / "data" / "work" / "openmrg" / "tables"
REL = REPO / "data" / "release" / "openmrg_v1"


def try_(name, path, fn):
    try:
        info = fn(path)
        print(f"OK   | {name} | {path.relative_to(REPO) if isinstance(path, Path) else path} | {info}")
    except Exception as e:  # report and go on
        print(f"FAIL | {name} | {path} | {type(e).__name__}: {e}")


csv = lambda p: f"{len(pd.read_csv(p))} rows"
js = lambda p: f"{len(json.loads(p.read_text(encoding='utf-8')))} keys"

try_("baseline v2 (frozen)", BP.BASELINES, csv)
try_("baseline uncertainty (offset_std)", BP.BASELINE_UNCERTAINTY, csv)
try_("baseline code", BP.MODULE / "baseline.py", lambda p: importlib.import_module("research.baseline.baseline") and "imports")
try_("TEST baselines (S8)", WP.S8_BASELINES, csv)
try_("dataset builder + offset c (system_id.build_dataset)", WP.MODULE / "system_id.py",
     lambda p: hasattr(importlib.import_module("research.wet_antenna.system_id"), "build_dataset") and "imports, build_dataset found")
try_("g_e TRAIN/VAL (S4.3)", WP.COND4C_G, csv)
try_("g_e TEST (S8)", WP.S8_G_E, csv)
try_("g_e library (event_factors.py)", WP.MODULE / "event_factors.py",
     lambda p: hasattr(importlib.import_module("research.wet_antenna.event_factors"), "rescaled") and "imports, rescaled found")
try_("length test library", WP.MODULE / "length_test.py",
     lambda p: hasattr(importlib.import_module("research.wet_antenna.length_test"), "length_fits") and "imports, length_fits found")
try_("length test S5a table", WP.S5A_LENGTH, csv)
try_("S6 data funnel (link selection)", WP.S6_FUNNEL, csv)
try_("S6 lengths per class", WP.S6_LENGTHS, csv)
try_("frozen equation json", WP.S9_EQUATION, js)
try_("equation reference json", WP.S9_REFERENCE, js)
try_("S8 freeze hashes", WP.S8_FREEZE, js)
try_("equation code (WAAEquation)", WP.MODULE / "waa_equation.py",
     lambda p: (lambda eq: f"runs, {len(eq.run([5.0] * 60, L_km=1.0, band='28 GHz')['W'])} min out")(
         importlib.import_module("research.wet_antenna.waa_equation").WAAEquation()))
try_("equation unit test", WP.MODULE / "tests" / "test_waa_equation.py", lambda p: p.exists() and "exists")
try_("config.yaml", WP.CONFIG_FILE, lambda p: "loads, sections: " + ", ".join(
     importlib.import_module("research.common.config").load_config()["wet_antenna"].keys().__iter__().__next__() for _ in [0]))
for f in ("events.csv", "split_events.csv", "split_links.csv", "event_links.csv"):
    try_(f"pipeline table {f}", TAB / f, csv)
try_("fit cache folder", WP.CACHE, lambda p: p.exists() and "exists" or (_ for _ in ()).throw(FileNotFoundError("missing")))


def release(p):
    from cmlpipe.api import CMLData
    d = CMLData.open(str(p))
    ev, ln = d.events(), d.links()
    stores = [s for s in ("radar_along_links_5min", "link_reference_state_5min", "cml_1min", "gauges_city_1min") if d.store(s) is not None]
    return f"opens; events by split {ev.split_event.value_counts().to_dict()}; links by split {ln.split_link.value_counts().to_dict()}; stores {stores}"


try_("release openmrg_v1 (CMLData)", REL, release)
