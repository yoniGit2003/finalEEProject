"""CMLData against the real release and work products (marker `data`: python -m pytest -m data).

Small windows only; skipped when the release has not been built on this machine."""

import netCDF4
import numpy as np
import pandas as pd
import pytest
import xarray as xr

from conftest import REAL_CONFIG

from cmlpipe.config import load_config

pytestmark = pytest.mark.data
DAY = pd.Timestamp("2015-07-28")


@pytest.fixture(scope="module")
def CFG():
    return load_config(REAL_CONFIG)


@pytest.fixture(scope="module")
def T(CFG):
    return CFG.tables_dir


@pytest.fixture(scope="module")
def d(CFG):
    from cmlpipe.api import CMLData
    rel = CFG.repo_path(CFG["paths"]["release_root"]) / "openmrg_v1"
    if not (rel / "catalog.json").exists():
        pytest.skip("release openmrg_v1 not built on this machine")
    return CMLData.open(rel)


def test_tables_match_work(d, T):
    assert len(d.links()) == len(pd.read_csv(T / "inventory/links.csv")) == 364
    s = d.sublinks()
    assert len(s) == 728 and s.sublink_id.str.match(r"^\d+_[12]$").all()
    assert len(d.status()) == 728 and len(d.gauges()) == 11 and "gauge_id" in d.gauges()
    assert d.dry_periods().shape[0] == len(pd.read_csv(T / "dry_periods.csv"))
    runs = d.event_runs("E20150727_2140")
    assert list(runs.run_id) == [f"E20150727_2140_R0{k}" for k in range(1, 7)]


def test_link_and_event_filters(d, T):
    w = pd.read_csv(T / "inventory/links.csv")
    exp = w[w.length_m.between(1000, 3000) & (w.frequency_band == ">=25GHz") & (w.status_best == "good")
            & (w.split_link == "train") & (w.n_gauges_within_1000m > 0)]
    got = d.links(length_km=(1, 3), frequency_band=">=25GHz", status="good", split="train", has_gauge_within_m=1000)
    assert sorted(got.cml_id) == sorted(exp.cml_id)
    sg = pd.read_csv(T / "event_significance.csv")
    ev = d.events(significant=True, split="test")
    assert sorted(ev.event_id) == sorted(sg[sg.significant & (sg.split_event == "test")].event_id)
    assert set(d.events(month=7).month) == {7}
    with pytest.raises(ValueError):
        d.links(no_such_filter=1)


def test_cml_l3_equals_work_file(d, CFG):
    ds = d.cml("10227_2", DAY, DAY + pd.Timedelta("1D"), "5min")
    with xr.open_dataset(CFG.output_root / "L3_products/cml/5min/cml_5min_2015-07-28.nc") as w:
        ref = w.sel(cml_id=10227, sublink_id="sublink_2")
        for v in ("tl_mean", "rsl_median", "n_valid_tl"):
            np.testing.assert_array_equal(ds[v].sel(sublink_id="10227_2").values, ref[v].values)


def test_cml_10s_raw_mask_flags(d, CFG):
    t0, t1 = DAY + pd.Timedelta("15h"), DAY + pd.Timedelta("17h")
    raw = d.cml(10227, t0, t1, "10s", qc="raw")
    fl = d.cml(10227, t0, t1, "10s", qc="flags")
    mk = d.cml(10227, t0, t1, "10s", qc="mask")
    with xr.open_dataset(CFG.l1_cml_day_path(DAY)) as l1:
        ref = l1.rsl.sel(cml_id=10227, sublink_id="sublink_1").sel(time=slice(t0, t1 - pd.Timedelta("1s"))).values
    np.testing.assert_allclose(raw.rsl.sel(sublink_id="10227_1").values, ref, atol=1e-6)
    assert "flags" in fl and "flags" not in raw and "flags" not in mk
    from cmlpipe.qc.flags import mask_of
    blocked = (fl.flags.values & mask_of(d.catalog["qc"]["blocking_flags_rsl"])) != 0
    assert np.isnan(mk.rsl.values[blocked]).all()
    assert np.array_equal(np.isnan(mk.rsl.values[~blocked]), np.isnan(raw.rsl.values[~blocked]))


def test_radar_gauge_met_states_equal_work(d, CFG):
    r = d.radar_along([10227], DAY, DAY + pd.Timedelta("1D"), "15min")
    with xr.open_dataset(CFG.output_root / "L3_products/radar/radar_along_links_15min.nc") as w:
        np.testing.assert_array_equal(r.R_mean.sel(cml_id=10227).values,
                                      w.R_mean.sel(cml_id=10227, time=slice(DAY, DAY + pd.Timedelta("23h45min"))).values)
    g = d.gauge(["Torsl", "GbgA"], DAY, DAY + pd.Timedelta("1D"), "15min")
    assert set(g.gauge_id.values) == {"Torsl", "GbgA"} and g.sizes["time"] == 96
    assert d.met(DAY, DAY + pd.Timedelta("1D")).sizes["time"] == 24
    ls = d.link_state([10227], DAY, DAY + pd.Timedelta("1D"))
    with xr.open_dataset(CFG.output_root / "L3_products/events/link_reference_state_5min.nc") as w:
        np.testing.assert_array_equal(ls.state.sel(cml_id=10227).values,
                                      w.state.sel(cml_id=10227, time=slice(DAY, DAY + pd.Timedelta("23h55min"))).values)
    assert d.network_state(DAY, DAY + pd.Timedelta("1D")).sizes["time"] == 288


def test_frame_aligned(d):
    f = d.frame(10227, "2015-07-28 12:00", "2015-07-28 20:00", "5min")
    assert len(f) == 96 and f.index.is_monotonic_increasing and (np.diff(f.index.values) == np.timedelta64(5, "m")).all()
    for col in ("tl_mean_s1", "tl_mean_s2", "flags_or_s1", "radar_R_mean", "gauge_rate_mm_h", "met_temperature",
                "met_filled", "link_state", "network_state"):
        assert col in f, col
    assert f.attrs["gauge"]["gauge_id"] == "Torsl"
    assert (~f.met_filled[f.index.minute == 0]).all() and f.met_filled[f.index.minute != 0].all()
    one = d.frame(10227, "2015-07-28 12:00", "2015-07-28 13:00", "1min")
    assert not any(c.startswith("radar_") for c in one) and any("radar" in n for n in one.attrs["notes"])


def test_event_window_markers(d):
    w = d.event_window("E20150727_2140", 10227, pad="1h")
    ev = d.events().set_index("event_id").loc["E20150727_2140"]
    assert w.in_event.sum() == int((ev.end - ev.start) / pd.Timedelta("5min"))
    assert set(w.run_id) - {""} == set(d.event_runs("E20150727_2140").run_id)
    net = d.event_window("E20150727_2140", pad="1h")
    assert "state" in net and net.in_event.sum() == w.in_event.sum()


def test_describe_and_errors(d):
    from cmlpipe.api import CMLDataError
    desc = d.describe("tl_mean")
    assert (desc.kind == "timeseries").any() and desc.description.str.len().gt(0).all()
    with pytest.raises(CMLDataError):
        d.describe("no_such_variable")
    with pytest.raises(CMLDataError, match="unknown cml_id"):
        d.cml(99999, DAY, DAY + pd.Timedelta("1h"))
    with pytest.raises(CMLDataError, match="large"):
        d.cml(10227, "2015-07-01", "2015-07-20", "10s")
    assert d.cml("10227_1", "2015-07-01", "2015-07-09", "10s", variables=["tl"], allow_large=True).sizes["time"] == 8 * 8640
    with pytest.raises(CMLDataError, match="not available"):
        d.radar_along(10227, DAY, DAY + pd.Timedelta("1h"), "1min")
    with pytest.raises(CMLDataError, match="outside the period"):
        d.cml(10227, "2015-05-01", "2015-06-02")


def test_times_utc_and_interval_start(d):
    a = d.cml(10227, pd.Timestamp("2015-07-28 14:00", tz="Europe/Stockholm"), pd.Timestamp("2015-07-28 15:00", tz="Europe/Stockholm"), "5min")
    assert pd.Timestamp(a.time.values[0]) == pd.Timestamp("2015-07-28 12:00")  # CEST = UTC+2
    assert a.sizes["time"] == 12  # [start, end): the bin labelled at the end is not included


def test_event_pairs_real(d):
    raw = d.event_pairs(usable_only=False, max_flatline_share=None, min_valid_fraction=None)
    assert len(raw) == 35 * 364  # significant events x all links, before any filter
    kept, dropped = d.event_pairs(return_dropped=True)
    assert len(kept) + len(dropped) == len(raw)
    flat = dropped[dropped.drop_reason.str.contains("flatline")]
    assert sorted(zip(flat.event_id, flat.cml_id)) == sorted(
        [("E20150828_2145", c) for c in (10004, 10015, 10077, 10115, 10126, 10330)]
        + [("E20150603_1500", 10329), ("E20150829_1445", 10076)])
    assert not kept.merge(flat[["event_id", "cml_id"]]).shape[0]
