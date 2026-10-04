"""Per-link reference state rules, release key normalisation, metadata coverage, release overwrite guard (synthetic)."""

import json

import numpy as np
import pandas as pd
import pytest

from cmlpipe.events import link_state
from cmlpipe.release import build, metadata

LCFG = {"radar_wet_mmh": 0.1, "gauge_wet_mm": 0.1, "gauge_max_dist_m": 1000, "margin_min": 30}
STEP = pd.Timedelta("5min")


def test_window_any_truncates_at_edges():
    m = np.zeros(20, bool)
    m[0] = True
    w = link_state.window_any(m[:, None], 6)[:, 0]
    assert w[:7].all() and not w[7:].any()


def test_link_state_rules():
    n = 60
    r = np.zeros((n, 3))
    r[20, 0] = 0.5                  # link 0: radar wet at bin 20
    r[40, 1] = 0.05                 # link 1: light radar (not wet) at bin 40 -> margin around it
    r[10, 2] = np.nan               # link 2: radar missing at bin 10
    g = np.zeros((n, 1))
    g[50, 0] = 0.2                  # the gauge near link 0 is wet at bin 50
    st, rs = link_state.classify_links(r, g, [[0], [], []], LCFG, STEP)
    W, D, U = 1, 0, -1
    assert st[20, 0] == W and rs[20, 0] & 1 and st[50, 0] == W and rs[50, 0] & 2
    assert st[14, 0] == U and st[13, 0] == D            # +-30 min margin around the wet bin is not dry
    assert st[40, 1] == U and rs[40, 1] & 8 and st[33, 1] == D
    assert st[10, 2] == U and rs[10, 2] & 4 and st[4, 2] == U and st[3, 2] == D
    assert st[0, 1] == D                                 # a link without a near gauge can be dry


def test_link_state_missing_gauge_blocks_dry():
    r = np.zeros((30, 1))
    g = np.zeros((30, 1))
    g[15, 0] = np.nan
    st, _ = link_state.classify_links(r, g, [[0]], LCFG, STEP)
    assert (st[9:22, 0] == -1).all() and st[8, 0] == 0


def test_normalise_keys():
    df = pd.DataFrame({"cml_id": [10227, 10227], "sublink_id": ["sublink_1", "sublink_2"], "x": [1, 2]})
    out = build.normalise("sublinks", df)
    assert list(out.sublink_id) == ["10227_1", "10227_2"] and list(out.sublink) == [1, 2]
    g = build.normalise("gauges", pd.DataFrame({"id": ["Torsl"]}))
    assert list(g.columns) == ["gauge_id"]


@pytest.mark.data
def test_every_table_column_described():
    from conftest import REAL_CONFIG
    from cmlpipe.config import load_config
    cfg = load_config(REAL_CONFIG)
    missing = []
    for name, (rel, _, _) in metadata.TABLES.items():
        p = cfg.output_root / rel
        if not p.exists():
            pytest.skip("work tables not built on this machine")
        for c in build.normalise(name, pd.read_csv(p, nrows=2)).columns:
            if not metadata.column_description(name, c):
                missing.append((name, c))
    assert not missing, missing


def test_release_refuses_overwrite_without_force(tmp_path, monkeypatch):
    class Cfg:
        repo_root = tmp_path

        def repo_path(self, rel):
            return tmp_path / rel

        def __getitem__(self, k):
            return {"paths": {"release_root": "rel"}}[k]
    out = tmp_path / "rel" / "openmrg_v1"
    out.mkdir(parents=True)
    (out / "catalog.json").write_text(json.dumps({"input_signature": "old"}), encoding="utf-8")
    monkeypatch.setattr(build, "input_signature", lambda cfg: "new")

    class Log:
        def info(self, *a, **k):
            pass
    with pytest.raises(SystemExit, match="read-only"):
        build.run_release(Cfg(), False, Log(), [], version="v1")
    rec = []
    monkeypatch.setattr(build, "input_signature", lambda cfg: "old")
    build.run_release(Cfg(), False, Log(), rec, version="v1")
    assert rec[-1]["status"] == "skipped_same_config"
