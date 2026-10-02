"""Step 7: significance rule, wet runs, completeness, TSL classes, quality status, link filters, quicklook link
selection, and that the step 7 stages leave the frozen tables byte-identical. Synthetic data only."""

import hashlib

import numpy as np
import pandas as pd
import pytest

from cmlpipe import events
from cmlpipe.events import stage as events_runner
from cmlpipe.api import explore
from cmlpipe import inventory
from cmlpipe.inventory import stage as inventory_runner
from cmlpipe.viz import quicklooks
from test_events import _Log, ev_world  # noqa: F401  (fixture re-used)

SIG = {"city_total_mean_mm": 1.0, "radar_domain_total_mm": 1.0, "city_max_5min_rate_mmh": 5.0, "heavy_total_mm": 10.0,
       "heavy_rate_mmh": 30.0, "run_total_mm": 0.5}


def ev_table(rows):
    return pd.DataFrame(rows, columns=["event_id", "city_total_mean_mm", "radar_domain_total_mm", "city_max_5min_rate_mmh"])


def test_significance_rule_and_classes():
    et = ev_table([("A", 0.2, 0.3, 1.2), ("B", 1.0, 0.0, 0.0), ("C", 0.0, 1.5, 0.0), ("D", 0.0, 0.0, 5.0),
                   ("E", 12.0, 3.0, 10.0), ("F", 2.0, 2.0, 30.0), ("G", np.nan, np.nan, np.nan)])
    s = inventory.significance(et, SIG).set_index("event_id")
    assert s.significant.to_dict() == {"A": False, "B": True, "C": True, "D": True, "E": True, "F": True, "G": False}
    assert s.rain_class.to_dict() == {"A": "none_drizzle", "B": "significant", "C": "significant", "D": "significant",
                                      "E": "heavy", "F": "heavy", "G": "none_drizzle"}
    assert s.loc["B", "sig_reason"] == "city_total" and s.loc["A", "sig_reason"] == "below all"


def _runs_fixture():
    W, D = events.STATE_CODES["wet"], events.STATE_CODES["dry"]
    n = 60
    t = pd.date_range("2015-07-28", periods=n, freq="5min")
    st = np.full(n, D, "int8")
    st[2:4] = W          # run 1 (in event E1)
    st[10:14] = W        # run 2 (in event E1), 6 bins = 30 min after run 1
    st[40:42] = W        # a short run that was dropped (no event id)
    eid = np.full(n, "", object)
    eid[2:14] = "E1"
    dm = np.zeros(n)
    dm[10:14] = 12.0     # 4 bins x 12 mm/h x 5/60 h = 4 mm domain total in run 2
    amount = np.zeros((n, 2))
    amount[2:4, 0] = 0.1
    return t, st, eid, dm, amount


def test_run_extraction_and_wet_core():
    t, st, eid, dm, amount = _runs_fixture()
    runs = inventory.event_runs(t, st, eid, dm, np.zeros(len(t)), np.zeros(len(t)), amount, pd.Timedelta("5min"), SIG)
    assert list(runs.run_id) == ["E1_R01", "E1_R02"]  # the dropped short run is not in any event
    r1, r2 = runs.iloc[0], runs.iloc[1]
    assert r1.duration_min == 10 and r2.duration_min == 20
    assert r2.gap_to_previous_run_h == pytest.approx(0.5) and np.isnan(r1.gap_to_previous_run_h)
    assert r2.radar_domain_total_mm == pytest.approx(4.0) and r1.city_total_mean_mm == pytest.approx(0.1)
    assert not r1.significant_run and r2.significant_run
    core = inventory.wet_core(runs).iloc[0]
    assert core.n_runs == 2 and core.wet_core_start == t[10] and core.wet_core_hours == pytest.approx(20 / 60)


def test_completeness_from_daily_summary():
    d = pd.DataFrame({"day": ["2015-06-01", "2015-07-01"], "cml_id": [1, 1], "sublink_id": ["sublink_1"] * 2,
                      "n_samples": [100, 100], "n_valid_rsl": [90, 100], "n_valid_tl": [80, 100], "n_missing_rsl": [10, 0],
                      "n_missing_tsl": [15, 0], "n_out_of_range": [0, 0], "n_rx_floor": [5, 0], "n_spike": [0, 2],
                      "n_tsl_step": [0, 0], "n_flatline": [0, 0], "n_network_gap": [0, 0], "n_blocking": [25, 2]})
    c = inventory.completeness(d).iloc[0]
    assert c.completeness_rsl == pytest.approx(0.95) and c.completeness_tsl == pytest.approx(0.925)
    assert c.completeness_tl == pytest.approx(0.9) and c.share_not_blocked == pytest.approx(173 / 200)
    assert c.flag_share_blocking_present == pytest.approx((180 - 173) / 200)  # missing data not counted twice
    assert c.completeness_tl_jun == pytest.approx(0.8) and c.completeness_tl_jul == pytest.approx(1.0) and np.isnan(c.completeness_tl_aug)


def test_tsl_behaviour_classes():
    tcfg = {"toggling_step_fraction": 0.01, "constant_top_share": 0.99, "constant_max_steps_per_day": 0.1, "rare_max_steps_per_day": 2.0}
    n = 8640 * 90
    cls = inventory.tsl_behaviour([0, n, n, n, n, n], [0, 0, 5, 100, 1000, n // 20], [np.nan, 1.0, 0.999, 0.999, 0.6, 0.5], tcfg)
    assert list(cls) == ["no_tsl", "constant", "constant", "rare_steps", "frequent_steps", "toggling"]


def _inv(n=40):
    rng = np.random.default_rng(0)
    df = pd.DataFrame({"completeness_tl": np.ones(n), "share_not_blocked": np.ones(n), "flag_share_blocking_present": np.zeros(n),
                       "flag_share_rx_floor": np.zeros(n), "flag_share_spike": np.zeros(n), "flag_share_flatline": np.zeros(n),
                       "flag_share_tsl_step": np.zeros(n), "completeness_tl_jun": np.ones(n), "completeness_tl_jul": np.ones(n),
                       "completeness_tl_aug": np.ones(n), "dry_tl_std": rng.uniform(0.1, 0.3, n),
                       "dry_baseline_range_db": rng.uniform(0.1, 0.5, n)})
    return df


STATUS = {"unusable_completeness_tl_pct": 50, "unusable_flagged_pct": 20, "unusable_rx_floor_pct": 20,
          "degraded_month_completeness_tl_pct": 80, "degraded_rx_floor_pct": 5,
          "degraded_flag_caps_pct": {"spike": 0.1, "flatline": 10, "tsl_step": 1}, "degraded_dry_noise_quantile": 0.95,
          "degraded_baseline_range_quantile": 0.95}


def test_status_classification():
    inv = _inv()
    inv.loc[0, "completeness_tl"] = 0.3                       # unusable
    inv.loc[1, ["completeness_tl", "share_not_blocked"]] = 0.6  # missing data only: NOT unusable via the flag rule
    inv.loc[1, "completeness_tl_jul"] = 0.5                   # ... but degraded (a month < 80 %)
    inv.loc[2, "flag_share_tsl_step"] = 0.05                  # degraded by the tsl_step cap
    inv.loc[3, "dry_tl_std"] = 5.0                            # degraded: above the network p95
    inv.loc[4, "flag_share_blocking_present"] = 0.25          # unusable: present samples flagged
    s = inventory.quality_status(inv, STATUS)
    assert list(s.status[:6]) == ["unusable", "degraded", "degraded", "degraded", "unusable", "good"]
    assert "completeness" in s.reason[0] and "a month" in s.reason[1] and "tsl_step" in s.reason[2] and "noise" in s.reason[3]
    assert "blocking-flagged" in s.reason[4] and s.reason[5] == "ok"


def _links():
    return pd.DataFrame({"cml_id": [1, 2, 3, 4, 5], "length_m": [500, 1500, 2500, 5000, 1200],
                         "frequency_band": [">=25GHz", ">=25GHz", "<25GHz", "<25GHz", ">=25GHz"],
                         "length_class": ["<1km", "1-3km", "1-3km", ">3km", "1-3km"],
                         "status_best": ["good", "good", "degraded", "good", "unusable"],
                         "status_sublink_1": ["good", "good", "degraded", "good", "unusable"],
                         "status_sublink_2": ["good", "degraded", "good", "good", "unusable"],
                         "split_link": ["train", "test", "test", "train", "excluded"],
                         "n_gauges_within_250m": [0, 1, 0, 0, 0], "n_gauges_within_1000m": [1, 1, 0, 0, 1],
                         "has_gauge_within_1km": [True, True, False, False, True],
                         "nearest_gauge_distance_m": [800, 100, 3000, np.nan, 500],
                         "tsl_behaviour_sublink_1": ["constant"] * 4 + ["toggling"], "tsl_behaviour_sublink_2": ["constant"] * 5})


def test_find_links_filters():
    lk = _links()
    f = lambda **k: list(explore.find_links(None, links=lk, **k).cml_id)  # noqa: E731
    assert f(length_km=(1, 3)) == [2, 3, 5]
    assert f(length_km=(1, 3), frequency_band=">=25GHz", status="good") == [2]
    assert f(split="test", has_gauge_within_m=1000) == [2]
    assert f(has_gauge_within_m=250) == [2] and f(tsl_behaviour="toggling") == [5]
    with pytest.raises(ValueError):
        f(has_gauge_within_m=300)


def test_quicklook_link_selection_rule():
    lk = _links()
    # cells in sorted order: 1-3km x <25GHz [3], 1-3km x >=25GHz [2] (5 is excluded), <1km x >=25GHz [1], >3km x <25GHz [4]
    assert quicklooks.select_links(lk, {"n_links": 4, "per_cell": 1}) == [3, 2, 1, 4]
    lk2 = pd.concat([lk, lk.iloc[[0]].assign(cml_id=6, split_link="test")], ignore_index=True)  # a 2nd <1km link, test split
    # the test-split link 6 now comes first in its cell; with per_cell 1 and n_links 5 the fill adds link 1
    assert quicklooks.select_links(lk2, {"n_links": 5, "per_cell": 1}) == [3, 2, 6, 4, 1]


def test_quantization_from_histograms():
    hr = np.zeros((1, 50), int)
    hr[0, [10, 13, 16, 20]] = 5    # steps of 0.3 and 0.4 dB
    ht = np.zeros((1, 20), int)
    ht[0, 5] = 99
    ht[0, 15] = 1
    q = inventory.quantization(hr, ht, 0.1).iloc[0]
    assert q.rsl_min_step_db == pytest.approx(0.3) and q.rsl_n_distinct == 4
    assert q.tsl_n_distinct_values == 2 and q.tsl_top_share == pytest.approx(0.99)


def test_significance_stage_keeps_frozen_tables(ev_world):  # noqa: F811
    cfg = ev_world
    events_runner.run_events(cfg, False, _Log(), [])
    events_runner.run_split(cfg, False, _Log(), [])
    events_runner.run_event_links(cfg, False, _Log(), [])
    frozen = [events_runner.table_path(cfg, k) for k in ("events", "split_events", "split_links", "event_links", "dry_periods")]
    before = [hashlib.sha256(p.read_bytes()).hexdigest() for p in frozen]
    rec = []
    inventory_runner.run_significance(cfg, False, _Log(), rec)
    assert rec[-1]["status"] == "written"
    sg = pd.read_csv(inventory_runner.out_paths(cfg)["event_significance"])
    et = pd.read_csv(events_runner.table_path(cfg, "events"))
    assert set(sg.event_id) == set(et.event_id) and (sg.set_index("event_id").n_runs == et.set_index("event_id").merged_from).all()
    inventory_runner.run_significance(cfg, False, _Log(), rec)
    assert rec[-1]["status"] == "skipped_same_config"
    assert [hashlib.sha256(p.read_bytes()).hexdigest() for p in frozen] == before
