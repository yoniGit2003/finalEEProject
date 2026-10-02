"""Events, dry periods, event x link table and the frozen split: synthetic data only.

Part 1: pure functions on hand-made state series and tables (exact answers).
Part 2: the tiny 2-day world of test_aggregate (radar L1, gauges, CML L1/L2, geometry) run through the
aggregate stages and then the events / split / event_links stages.
"""

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from cmlpipe.aggregate import stage as aggregate_runner
from cmlpipe import events
from cmlpipe.events import stage as events_runner
from test_aggregate import DAY

STEP = pd.Timedelta("5min")
W, D, U = events.STATE_CODES["wet"], events.STATE_CODES["dry"], events.STATE_CODES["unknown"]
ECFG = {"radar_min_valid_cell_fraction": 0.5, "min_wet_area_fraction": 0.01, "min_gauges_wet": 1,
        "dry_max_wet_area_fraction": 0.0}


def states(*parts):
    """states((D, 10), (W, 3), ...) -> array"""
    return np.concatenate([np.full(n, s, dtype="int8") for s, n in parts])


def find(st, gap_h=2, min_min=15):
    return events.find_events(st, STEP, pd.Timedelta(hours=gap_h), pd.Timedelta(minutes=min_min))


# ---------------------------------------------------------------------------
# Part 1: events and dry periods
# ---------------------------------------------------------------------------
def test_merge_90_min_apart_with_2h_gap():
    st = states((D, 50), (W, 6), (D, 18), (W, 6), (D, 50))  # 18 bins = 90 min apart
    ev, dropped = find(st, gap_h=2)
    assert len(ev) == 1 and ev[0]["merged_from"] == 2 and (ev[0]["i0"], ev[0]["i1"]) == (50, 80)


def test_not_merged_3h_apart():
    st = states((D, 50), (W, 6), (D, 36), (W, 6), (D, 50))  # 36 bins = 3 h
    ev, _ = find(st, gap_h=2)
    assert len(ev) == 2


def test_unknown_counts_as_not_wet_for_merging():
    st = states((D, 50), (W, 6), (U, 36), (W, 6), (D, 50))
    assert len(find(st, gap_h=2)[0]) == 2  # 3 h of unknown separates like 3 h of dry


def test_short_run_dropped_and_counted():
    st = states((D, 50), (W, 2), (D, 50), (W, 3), (D, 50))  # 10 min and 15 min
    ev, dropped = find(st, min_min=15)
    assert dropped == 1 and len(ev) == 1 and ev[0]["i1"] - ev[0]["i0"] == 3


def test_event_at_period_start_kept_and_flagged():
    st = states((W, 6), (D, 100), (W, 6), (D, 100))
    ev, _ = find(st)
    assert ev[0]["i0"] == 0
    gap = int(pd.Timedelta("2h") / STEP)
    assert events.touches_edge(ev[0]["i0"], ev[0]["i1"], st.size, gap)
    assert not events.touches_edge(ev[1]["i0"], ev[1]["i1"], st.size, gap)


def test_event_id_stable_format():
    assert events.event_id(pd.Timestamp("2015-07-28 14:55")) == "E20150728_1455"


def test_unknown_never_dry():
    st = states((D, 20), (U, 5), (D, 20))
    dry = events.certified_dry(st, margin_bins=0)
    assert not dry[20:25].any() and dry[:20].all()
    assert events.runs(dry) == [(0, 20), (25, 45)]  # the unknown bins split the dry period


def test_dry_margin_removes_30_min():
    st = states((D, 40), (W, 4), (D, 40))
    dry = events.certified_dry(st, margin_bins=6)
    assert dry[:34].all() and not dry[34:50].any() and dry[50:].all()


def test_gauge_only_and_radar_only_events():
    n = 300
    ref = pd.DataFrame({"radar_valid_cell_fraction": np.ones(n), "radar_wet_area_fraction": np.zeros(n),
                        "n_city_gauges_wet": np.zeros(n, int)}, index=pd.date_range(DAY, periods=n, freq="5min"))
    ref.iloc[50:56, ref.columns.get_loc("n_city_gauges_wet")] = 1      # gauge-only rain, radar dry
    ref.iloc[200:206, ref.columns.get_loc("radar_wet_area_fraction")] = 0.2  # radar-only rain, gauges dry
    c = events.classify(ref, ECFG)
    ev, _ = find(c.state.to_numpy())
    assert [(e["i0"], e["i1"]) for e in ev] == [(50, 56), (200, 206)]
    assert c.gauge_wet.iloc[50] and not c.radar_wet.iloc[50] and c.radar_wet.iloc[200] and not c.gauge_wet.iloc[200]


def test_radar_missing_is_unknown_and_small_area_not_dry():
    ref = pd.DataFrame({"radar_valid_cell_fraction": [1.0, 0.0, 1.0, 1.0], "radar_wet_area_fraction": [0.0, np.nan, 0.005, 0.0],
                        "n_city_gauges_wet": [0, 0, 0, 1]})
    st = events.classify(ref, ECFG).state.tolist()
    assert st == [D, U, U, W]


def test_context_counts_dry_not_unknown():
    st = states((D, 30), (U, 6), (W, 6), (D, 10))
    ev, _ = find(st)
    cx = events.context(st, ev, 0, cap_bins=1000)
    assert (cx["gap_before"], cx["dry_before"], cx["gap_after"], cx["dry_after"]) == (36, 30, 10, 10)


def test_dry_stats_only_dry_bins():
    t = pd.date_range(DAY, periods=288 * 6, freq="5min")
    rng = np.random.default_rng(0)
    tl = 50 + rng.normal(0, 0.2, (t.size, 2))
    dry = np.ones(t.size, bool)
    dry[100:110] = False
    tl[105, 0] = 1e6  # a huge value in a wet bin
    core = {"tl_mean": tl, "tl_std": np.full_like(tl, 0.1), "rsl_mean": -tl, "rsl_std": np.full_like(tl, 0.1)}
    s = events.dry_stats(t, dry, core, min_days=5)
    assert s["tl_p99"][0] < 51 and s["tl_mean"][0] < 51
    assert s["n_valid_tl"][0] == t.size - 10 and s["n_dry_days"][0] == 6
    assert np.isfinite(s["tl_diurnal_range"][0]) and np.isfinite(s["tl_daily_median_drift_std"][0])
    s3 = events.dry_stats(t[:288 * 3], dry[:288 * 3], {k: v[:288 * 3] for k, v in core.items()}, min_days=5)
    assert np.isnan(s3["tl_diurnal_range"][0])  # fewer than 5 dry days


# ---------------------------------------------------------------------------
# Part 1: split
# ---------------------------------------------------------------------------
SCFG = {"seed": 42, "event_fractions": {"train": 0.6, "val": 0.2, "test": 0.2}, "link_fractions": {"train": 0.8, "test": 0.2},
        "length_classes_km": [1, 3], "frequency_split_ghz": 25}


def fake_events(n=30):
    return pd.DataFrame({"event_id": [f"E201507{d:02d}_0000" for d in range(1, n + 1)],
                         "intensity_class": ["light", "moderate", "heavy"] * (n // 3)})


def fake_links():
    rows = []
    for k in range(40):
        rows.append({"cml_id": 1000 + k, "site_0_id": 2 * k, "site_1_id": 2 * k + 1,
                     "length_m": [500, 2000, 6000][k % 3], "frequency_ghz_sublink_1": [18.0, 38.0][k % 2]})
    rows.append({"cml_id": 2000, "site_0_id": 1, "site_1_id": 0, "length_m": 500, "frequency_ghz_sublink_1": 18.0})  # parallel to 1000
    lk = pd.DataFrame(rows)
    st = pd.DataFrame([{"cml_id": c, "sublink_id": s, "status": "good"} for c in lk.cml_id for s in ("sublink_1", "sublink_2")])
    st.loc[st.cml_id == 1005, "status"] = "unusable"
    return lk, st


def test_split_events_sizes_stratified_deterministic():
    ev = fake_events(30)
    a, b = events.split_events(ev, SCFG["event_fractions"], 42), events.split_events(ev, SCFG["event_fractions"], 42)
    pd.testing.assert_frame_equal(a, b)
    for cls, g in a.groupby("intensity_class"):
        assert g.split_event.value_counts().to_dict() == {"train": 6, "val": 2, "test": 2}, cls
    assert set(a.event_id) == set(ev.event_id)
    assert not a.equals(events.split_events(ev, SCFG["event_fractions"], 7))  # the seed matters


def test_split_links_parallel_together_and_excluded():
    lk, st = fake_links()
    s = events.split_links(lk, st, SCFG).set_index("cml_id")
    assert s.loc[1000, "split_link"] == s.loc[2000, "split_link"]  # same hop
    assert s.loc[1005, "split_link"] == "excluded"
    used = s[s.split_link != "excluded"]
    share = (used.split_link == "test").mean()
    assert 0.15 <= share <= 0.25
    pd.testing.assert_frame_equal(events.split_links(lk, st, SCFG), events.split_links(lk, st, SCFG))
    for stratum, g in used.groupby("stratum"):
        assert (g.split_link == "train").any(), stratum  # every stratum is represented in train


# ---------------------------------------------------------------------------
# Part 2: stages on the tiny world
# ---------------------------------------------------------------------------
@pytest.fixture
def ev_world(aggregated_world, tmp_path, synthetic_raw_dir):
    """A private copy of the aggregated world (copying the small output folder is much faster than rebuilding it)."""
    import shutil
    from conftest import make_config
    from cmlpipe.config import load_config
    src = aggregated_world[0]
    shutil.copytree(src.output_root, tmp_path / "out")
    cfg = load_config(make_config(tmp_path, synthetic_raw_dir / "cml.nc", synthetic_raw_dir / "cml_metadata.csv",
                                  tmp_path / "out"), repo_root=tmp_path)
    cfg.data["period"] = dict(src.data["period"])  # the world's two-day period (set by make_world)
    pd.DataFrame(columns=["start", "end", "duration_min", "max_fraction_nan", "n_alive"]).to_csv(
        cfg.tables_dir / "network_gaps.csv", index=False)
    t = pd.date_range("2015-07-28 01:00", "2015-07-30 00:00", freq="1h")
    met = xr.Dataset({v: (("time", "id"), np.full((t.size, 1), 15.0)) for v in ("temperature", "humidity", "wind_speed", "wind_gust")},
                     coords={"time": t, "id": ["GbgA"]})
    cfg.l1_met_dir.mkdir(parents=True, exist_ok=True)
    met.attrs["stage_hash"] = "met"
    met.to_netcdf(cfg.l1_met_dir / "smhi_met_1h.nc")
    cfg.data["smhi_gauge"]["id"] = "G2"  # the world's SMHI gauge is called G2
    return cfg


class _Log:
    def info(self, *a, **k):
        pass


def test_definition_reads_no_smhi_and_no_cml(ev_world, monkeypatch):
    cfg = ev_world

    def boom(*a, **k):
        raise AssertionError("must not be read by the event definition")
    monkeypatch.setattr(events_runner, "read_smhi", boom)
    monkeypatch.setattr(events_runner, "read_core", boom)
    monkeypatch.setattr(events_runner, "build_core", boom)
    ref, info = events_runner.build_reference(cfg)
    ref, ev, dropped, dry = events_runner.define_events(cfg, ref)
    assert len(ref) == 2 * 288 and info["n_domain_cells"] > 0
    assert len(ev) > 0  # the world's City gauge rains often


def test_stages_end_to_end(ev_world, capsys):
    cfg = ev_world
    rec = []
    events_runner.run_events(cfg, False, _Log(), rec)
    events_runner.run_split(cfg, False, _Log(), rec)
    events_runner.run_event_links(cfg, False, _Log(), rec)
    et = pd.read_csv(events_runner.table_path(cfg, "events"), parse_dates=["start", "end"])
    events.validate_events(et, pd.read_csv(events_runner.table_path(cfg, "split_events")))
    assert "split_event" in et and et.split_event.notna().all()

    # event x link radar sums equal the radar product sums
    el = pd.read_csv(events_runner.table_path(cfg, "event_links"))
    with xr.open_dataset(aggregate_runner.radar_path(cfg, "along_links", "5min")) as r:
        R = r.R_mean.load()
    for row in el.sample(min(20, len(el)), random_state=0).itertuples():
        e = et.set_index("event_id").loc[row.event_id]
        v = R.sel(cml_id=row.cml_id, time=slice(e.start, e.end - STEP)).values.astype("float64")
        assert row.radar_total_mm == pytest.approx(np.nansum(v) * 5 / 60, abs=1e-6)
    assert el.split_link.notna().all() and set(el.event_id) == set(et.event_id)

    # second run: everything skipped, nothing rewritten
    before = {k: events_runner.table_path(cfg, k).read_bytes() for k in events_runner.TABLES}
    rec2 = []
    events_runner.run_events(cfg, False, _Log(), rec2)
    events_runner.run_split(cfg, False, _Log(), rec2)
    events_runner.run_event_links(cfg, False, _Log(), rec2)
    assert [r["status"] for r in rec2] == ["skipped_same_config", "frozen_kept", "skipped_same_config"]
    assert {k: events_runner.table_path(cfg, k).read_bytes() for k in events_runner.TABLES} == before

    # forced rebuild of the events gives identical tables (ids stable)
    events_runner.run_events(cfg, True, _Log(), [])
    events_runner.run_split(cfg, False, _Log(), [])
    for k in ("events", "dry_periods", "calibration", "dry_stats"):
        assert events_runner.table_path(cfg, k).read_bytes() == before[k], k


def test_split_frozen(ev_world, capsys):
    cfg = ev_world
    events_runner.run_events(cfg, False, _Log(), [])
    events_runner.run_split(cfg, False, _Log(), [])
    paths = [events_runner.table_path(cfg, k) for k in ("split_events", "split_links")]
    first = [p.read_bytes() for p in paths]
    cfg.data["split"]["seed"] = 7  # a new seed must NOT redraw a frozen split
    rec = []
    events_runner.run_split(cfg, False, _Log(), rec)
    assert rec[-1]["status"] == "frozen_kept" and "FROZEN" in capsys.readouterr().out
    assert [p.read_bytes() for p in paths] == first
    events_runner.run_split(cfg, True, _Log(), rec)  # --force-split redraws
    assert rec[-1]["status"] == "written"


def test_event_without_split_is_an_error(ev_world):
    cfg = ev_world
    events_runner.run_events(cfg, False, _Log(), [])
    events_runner.run_split(cfg, False, _Log(), [])
    sp = pd.read_csv(events_runner.table_path(cfg, "split_events"))
    sp.iloc[1:].to_csv(events_runner.table_path(cfg, "split_events"), index=False)  # drop one event's row
    with pytest.raises(events.EventsValidationError, match="no split row"):
        events_runner.run_split(cfg, False, _Log(), [])
