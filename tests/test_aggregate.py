"""Aggregation (L3) and access: synthetic data only (no real files needed).

Part 1 tests the pure reductions on hand-made arrays (exact numbers).
Part 2 builds a tiny 2-day world (period 2015-07-28 .. 2015-07-29): synthetic L1/L2 CML files, gauge L1,
a radar L1 file and the geometry stage, then runs the three aggregate stages and the access function.
"""

import hashlib

import numpy as np
import pandas as pd
import pyproj
import pytest
import xarray as xr

from cmlpipe import aggregate
from cmlpipe.aggregate import stage as aggregate_runner
from cmlpipe import geometry
from cmlpipe.geometry import stage as geometry_runner
from cmlpipe.io import schema
from cmlpipe.api.access import load_link_window, window_edges
from cmlpipe.parallel import run_days

DAY = pd.Timestamp("2015-07-28")
NATIVE = "10s"
N_DAY = 8640
BIT = {k: np.uint16(1 << b) for k, b in schema.FLAGS.items()}
ACFG = {
    "min_valid_fraction": 0.8,
    "blocking_flags": {"rsl": ["missing_rsl", "out_of_range", "rx_floor", "spike", "network_gap"],
                       "tsl": ["missing_tsl", "out_of_range", "network_gap"]},
    "stats": {"default": {"tl": ["mean", "median", "min", "max", "std"], "rsl": ["mean", "median", "min", "max", "std"],
                          "tsl": ["mean", "min", "max"]}},
}


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def cml_day(rsl, tsl, flags, day=DAY, cml_ids=(1,)):
    """L1 + L2 datasets of one day from (time, sublink, link) arrays."""
    t = pd.date_range(day, periods=N_DAY, freq=NATIVE)
    dims = schema.DIMS
    coords = {"time": t, "sublink_id": list(schema.SUBLINK_IDS), "cml_id": list(cml_ids),
              "length": ("cml_id", np.full(len(cml_ids), 1.5)),
              "frequency": (("sublink_id", "cml_id"), np.full((2, len(cml_ids)), 28.0))}
    l1 = xr.Dataset({"rsl": (dims, rsl), "tsl": (dims, tsl), "tl": (dims, tsl - rsl)}, coords=coords)
    l2 = xr.Dataset({schema.FLAGS_VAR: (dims, flags.astype("uint16"), schema.flag_attrs())},
                    coords={k: coords[k] for k in ("time", "sublink_id", "cml_id")})
    return l1, l2


def base_arrays(n_links=1, seed=0):
    rng = np.random.default_rng(seed)
    rsl = np.round(-45 + rng.normal(0, 0.5, (N_DAY, 2, n_links)), 1)
    tsl = np.full((N_DAY, 2, n_links), 10.0)
    return rsl, tsl, np.zeros((N_DAY, 2, n_links), dtype="uint16")


def bin_of(ds, ts):
    return ds.sel(time=pd.Timestamp(ts))


# ---------------------------------------------------------------------------
# Part 1: CML reductions
# ---------------------------------------------------------------------------
def test_exact_bin_statistics():
    rsl, tsl, fl = base_arrays()
    l1, l2 = cml_day(rsl, tsl, fl)
    ds = aggregate.aggregate_cml_day(l1, l2, "15min", ACFG, NATIVE)
    v = rsl[90:180, 0, 0]  # the second 15 min bin of sublink_1
    b = ds.isel(time=1, sublink_id=0, cml_id=0)
    assert float(b.rsl_mean) == pytest.approx(v.mean(), abs=1e-5)
    assert float(b.rsl_median) == pytest.approx(np.median(v), abs=1e-5)
    assert float(b.rsl_min) == pytest.approx(v.min()) and float(b.rsl_max) == pytest.approx(v.max())
    assert float(b.rsl_std) == pytest.approx(v.std(ddof=0), abs=1e-5)  # population std
    assert int(b.n_valid_rsl) == 90
    assert float(b.tl_mean) == pytest.approx((10.0 - v).mean(), abs=1e-5)
    assert ds.sizes["time"] == 96 and ds.attrs["n_expected"] == 90 and ds.attrs["std_ddof"] == 0


def test_blocking_flags_excluded_informational_not():
    rsl, tsl, fl = base_arrays()
    rsl[0:90, 0, 0] = -45.0
    rsl[5, 0, 0] = -90.0
    fl[5, 0, 0] = BIT["spike"]          # blocking for rsl -> excluded
    rsl[6, 0, 0] = -60.0
    fl[6, 0, 0] = BIT["tsl_step"]       # informational -> included
    rsl[7, 0, 0] = -60.0
    fl[7, 0, 0] = BIT["flatline"]       # informational -> included
    l1, l2 = cml_day(rsl, tsl, fl)
    b = aggregate.aggregate_cml_day(l1, l2, "15min", ACFG, NATIVE).isel(time=0, sublink_id=0, cml_id=0)
    assert int(b.n_valid_rsl) == 89
    assert float(b.rsl_min) == pytest.approx(-60.0)
    assert float(b.rsl_mean) == pytest.approx((87 * -45.0 + 2 * -60.0) / 89, abs=1e-4)
    assert int(b.n_tsl_step) == 1
    assert int(b.flags_or) == int(BIT["spike"] | BIT["tsl_step"] | BIT["flatline"])


@pytest.mark.parametrize("n_missing, expect_value", [(19, False), (18, True)])
def test_valid_fraction_rule(n_missing, expect_value):
    # 15 min = 90 samples: 71 valid = 78.9 % -> NaN; 72 valid = 80.0 % -> value
    rsl, tsl, fl = base_arrays()
    rsl[:n_missing, 0, 0] = np.nan
    fl[:n_missing, 0, 0] = BIT["missing_rsl"]
    l1, l2 = cml_day(rsl, tsl, fl)
    b = aggregate.aggregate_cml_day(l1, l2, "15min", ACFG, NATIVE).isel(time=0, sublink_id=0, cml_id=0)
    assert int(b.n_valid_rsl) == 90 - n_missing  # n_valid always kept
    assert np.isfinite(float(b.rsl_mean)) == expect_value
    assert np.isfinite(float(b.rsl_median)) == expect_value
    assert np.isfinite(float(b.tsl_mean))  # tsl untouched


def test_unflagged_spike_moves_mean_not_median():
    rsl, tsl, fl = base_arrays()
    rsl[0:90, 0, 0] = -45.0
    rsl[10, 0, 0] = -75.0  # 30 dB dip, not flagged
    l1, l2 = cml_day(rsl, tsl, fl)
    b = aggregate.aggregate_cml_day(l1, l2, "15min", ACFG, NATIVE).isel(time=0, sublink_id=0, cml_id=0)
    assert float(b.rsl_mean) == pytest.approx(-45.0 - 30 / 90, abs=1e-4)
    assert float(b.rsl_median) == pytest.approx(-45.0)


def test_tl_needs_rsl_and_tsl_valid():
    rsl, tsl, fl = base_arrays()
    fl[0:30, 0, 0] = BIT["spike"]      # rsl blocked, tsl fine
    fl[30:50, 0, 0] = BIT["missing_tsl"]
    tsl[30:50, 0, 0] = np.nan
    l1, l2 = cml_day(rsl, tsl, fl)
    b = aggregate.aggregate_cml_day(l1, l2, "15min", ACFG, NATIVE).isel(time=0, sublink_id=0, cml_id=0)
    assert int(b.n_valid_rsl) == 60 and int(b.n_valid_tsl) == 70 and int(b.n_valid_tl) == 40
    assert np.isnan(float(b.tl_mean))  # 40 / 90 < 0.8
    assert np.isnan(float(b.rsl_mean)) and np.isnan(float(b.tsl_mean))  # 60/90 and 70/90 < 0.8


def test_flags_or_bitwise():
    rsl, tsl, fl = base_arrays()
    fl[1, 1, 0] = BIT["rx_floor"]
    fl[2, 1, 0] = BIT["network_gap"] | BIT["missing_rsl"]
    rsl[2, 1, 0] = np.nan
    l1, l2 = cml_day(rsl, tsl, fl)
    ds = aggregate.aggregate_cml_day(l1, l2, "1min", ACFG, NATIVE)
    assert int(ds.flags_or.isel(time=0, sublink_id=1, cml_id=0)) == int(BIT["rx_floor"] | BIT["network_gap"] | BIT["missing_rsl"])
    assert int(ds.flags_or.isel(time=1, sublink_id=1, cml_id=0)) == 0


def test_days_separately_equal_concatenated():
    """No bin crosses midnight: per-day aggregation == aggregation of the 2-day series at bin level."""
    a, b = base_arrays(seed=1), base_arrays(seed=2)
    out = []
    for (rsl, tsl, fl), day in ((a, DAY), (b, DAY + pd.Timedelta("1D"))):
        l1, l2 = cml_day(rsl, tsl, fl, day=day)
        out.append(aggregate.aggregate_cml_day(l1, l2, "1h", ACFG, NATIVE))
    daily = xr.concat(out, "time")
    both = np.concatenate([a[0], b[0]])
    stats, n = aggregate.bin_statistics(both, ~np.isnan(both), 360, ["mean", "median", "max"], 288)
    assert np.allclose(daily.rsl_mean.values, stats["mean"], atol=1e-4)
    assert np.allclose(daily.rsl_median.values, stats["median"], atol=1e-4)
    assert np.array_equal(daily.n_valid_rsl.values, n)


def test_min_count_integer_exact():
    assert aggregate.min_count(0.8, 90) == 72 and aggregate.min_count(0.6, 3) == 2 and aggregate.min_count(0.6, 12) == 8
    assert aggregate.min_count(1.0, 5) == 5


# ---------------------------------------------------------------------------
# Part 1: gauges and radar reductions
# ---------------------------------------------------------------------------
def gauge_l1(amount, start=DAY, native="1min"):
    t = pd.date_range(start, periods=amount.shape[0], freq=native)
    return xr.Dataset({"rainfall_amount": (("time", "id"), amount)},
                      coords={"time": t, "id": ["A", "B"], "quantization_mm": ("id", [0.1, 0.2])})


def test_gauge_sums_and_strict_rule():
    rng = np.random.default_rng(3)
    amount = np.round(rng.choice([0, 0, 0, 0.1, 0.2, 0.5], (1440, 2)), 1)
    amount[602, 1] = np.nan  # gauge B misses 10:02
    l1 = gauge_l1(amount)
    for res in ("1min", "5min", "15min", "1h"):
        ds = aggregate.aggregate_gauges(l1, res, "1min", DAY, DAY + pd.Timedelta("1D"), 1.0, peak_rate=True)
        assert float(ds.amount_mm.sel(id="A").sum()) == pytest.approx(amount[:, 0].sum(), abs=1e-3)
        dt_h = pd.Timedelta(res).total_seconds() / 3600
        assert np.allclose(ds.rate_mm_h.values, ds.amount_mm.values / dt_h, equal_nan=True)
        schema.validate_l3(ds, "gauges", res, start=DAY, n_bins=int(pd.Timedelta("1D") / pd.Timedelta(res)))
    ds5 = aggregate.aggregate_gauges(l1, "5min", "1min", DAY, DAY + pd.Timedelta("1D"), 1.0, peak_rate=True)
    b = ds5.sel(time="2015-07-28 10:00", id="B")
    assert np.isnan(float(b.amount_mm)) and int(b.n_valid) == 4  # strict: not filled
    a = ds5.sel(time="2015-07-28 10:00", id="A")
    assert float(a.peak_1min_rate_mm_h) == pytest.approx(amount[600:605, 0].max() * 60)


def test_smhi_not_disaggregated():
    with pytest.raises(ValueError):
        aggregate.samples_per_bin("5min", "15min")
    with pytest.raises(ValueError):
        aggregate.samples_per_bin("1min", "5min")  # radar at 1 min


def test_floor_sets_zero_not_phantom():
    dbz = np.array([-30.0, -5.0, 0.0, 20.0, np.nan])
    rain = schema.zr_rain_rate(dbz, 200.0, 1.5)
    assert rain[0] == pytest.approx(2.92e-4, rel=1e-2)  # (1e-3 / 200)^(1/1.5): the phantom rain the floor removes
    floored, n_zero, n_fin = aggregate.floor_rain(dbz, rain, 0.0)
    assert floored[0] == 0.0 and floored[1] == 0.0 and floored[2] > 0 and np.isnan(floored[4])
    assert (n_zero, n_fin) == (2, 4)


def test_stamp_alignment_and_15min_mean():
    stamps = pd.date_range(DAY, periods=289, freq="5min")  # includes the closing instant of the day
    vals = np.zeros((289, 1))
    vals[13, 0] = 12.0  # spike at 01:05
    grid, dropped = aggregate.stamps_to_grid(stamps.values, vals, DAY, DAY + pd.Timedelta("1D"), "5min")
    assert dropped == 1 and grid.shape == (288, 1)
    assert grid[13, 0] == 12.0  # bin 01:05, not 01:00
    out, n = aggregate.aggregate_stamps({"R": grid}, ~np.isnan(grid), 3, 0.6, {"R": "mean"})
    assert out["R"][4, 0] == pytest.approx(4.0)  # 15 min bin 01:00 = mean of 01:00, 01:05, 01:10 = 12/3
    assert out["R"][3, 0] == 0.0 and out["R"][5, 0] == 0.0
    grid[12, 0] = grid[14, 0] = np.nan  # only 1 of 3 stamps valid -> NaN (needs 2)
    out, n = aggregate.aggregate_stamps({"R": grid}, ~np.isnan(grid), 3, 0.6, {"R": "mean"})
    assert np.isnan(out["R"][4, 0]) and n[4, 0] == 1


# ---------------------------------------------------------------------------
# Part 2: the tiny 2-day world
# ---------------------------------------------------------------------------
PROJ = "+proj=stere +lat_ts=60 +ellps=bessel +lon_0=14 +lat_0=90"
IDS = [100, 200, 300, 400]  # conftest metadata links
DAYS = [DAY, DAY + pd.Timedelta("1D")]
SPIKE = pd.Timestamp("2015-07-28 06:10")  # uniform 10 mm/h over the whole radar grid at this stamp


def _write(ds, path, stage_hash, **attrs):
    path.parent.mkdir(parents=True, exist_ok=True)
    ds.attrs.update({"stage_hash": stage_hash, **attrs})
    ds.to_netcdf(path)


def make_world(cfg):
    """Write the tiny two-day world (L1/L2 CML, gauges, radar, geometry, QC status) into cfg's output folder."""
    cfg.data["period"]["start"], cfg.data["period"]["end"] = "2015-07-28", "2015-07-29"
    rng = np.random.default_rng(7)
    # --- CML L1 + L2 ---
    for k, day in enumerate(DAYS):
        rsl, tsl, fl = base_arrays(len(IDS), seed=10 + k)
        fl[100:103, 0, 2] = BIT["spike"]
        l1, l2 = cml_day(rsl, tsl, fl, day=day, cml_ids=IDS)
        _write(l1, cfg.l1_cml_day_path(day), "l1hash")
        from cmlpipe.qc import stage as qc_runner
        _write(l2, qc_runner.l2_day_path(cfg, day), "l2hash", stats_stage_hash="statshash")
    # --- gauges: City 'G1' on link 300, SMHI 'G2' ~500 m from it ---
    cfg.tables_dir.mkdir(parents=True, exist_ok=True)
    gtab = pd.DataFrame({"id": ["G1", "G2"], "network": ["city", "smhi"], "gauge_type": ["weighing"] * 2,
                         "lon": [11.85, 11.858], "lat": [57.75, 57.75], "quantization_mm": [0.1, 0.1]})
    gtab.to_csv(cfg.tables_dir / "gauges.csv", index=False)
    city = rng.choice([0.0, 0.0, 0.1, 0.3], (2880, 1))
    city[700, 0] = np.nan
    cds = xr.Dataset({"rainfall_amount": (("time", "id"), city)},
                     coords={"time": pd.date_range(DAY, periods=2880, freq="1min"), "id": ["G1"],
                             "network": ("id", ["city"]), "quantization_mm": ("id", [0.1])})
    _write(cds, cfg.l1_gauges_dir / "city_1min.nc", "ghash")
    smhi = rng.choice([0.0, 0.2], (192, 1))
    sds = xr.Dataset({"rainfall_amount": (("time", "id"), smhi), "qc_code": (("time", "id"), np.ones((192, 1), "uint8"))},
                     coords={"time": pd.date_range(DAY, periods=192, freq="15min"), "id": ["G2"],
                             "network": ("id", ["smhi"]), "quantization_mm": ("id", [0.1])})
    _write(sds, cfg.l1_gauges_dir / "smhi_15min.nc", "ghash")
    # --- radar L1: -30 dBZ everywhere (phantom 0.0003 mm/h) except a uniform 10 mm/h stamp ---
    rx, ry = geometry.project_to_radar([11.3, 12.5], [56.8, 57.9], PROJ)
    gx = np.arange(np.floor(rx.min() / 2000) * 2000 + 1000, rx.max() + 2000, 2000.0)
    gy = np.arange(np.ceil(ry.max() / 2000) * 2000 - 1000, ry.min() - 2000, -2000.0)
    crs = pyproj.CRS.from_proj4(PROJ)
    xx, yy = np.meshgrid(gx, gy)
    lon, lat = pyproj.Transformer.from_crs(crs, crs.geodetic_crs, always_xy=True).transform(xx, yy)
    x3006, y3006 = geometry.project_points(lon, lat, "EPSG:3006")
    stamps = pd.date_range(DAY, DAY + pd.Timedelta("2D"), freq="5min")  # incl. the closing instant
    dbz = np.full((len(stamps), gy.size, gx.size), -30.0, dtype="float32")
    r10 = 10 * np.log10(200 * 10.0 ** 1.5)
    dbz[stamps.get_loc(SPIKE)] = r10
    dbz[stamps.get_loc(pd.Timestamp("2015-07-29 12:00"))] = np.nan  # a fully missing stamp
    rds = xr.Dataset({"dbz": (("time", "y", "x"), dbz),
                      "rain_rate": (("time", "y", "x"), schema.zr_rain_rate(dbz, 200.0, 1.5).astype("float32"))},
                     coords={"time": stamps, "y": gy, "x": gx, "lon": (("y", "x"), lon), "lat": (("y", "x"), lat),
                             "x_3006": (("y", "x"), x3006), "y_3006": (("y", "x"), y3006)})
    _write(rds, cfg.l1_radar_dir / "radar_5min_2015-07.nc", "rhash", crs_proj4=PROJ, zr_a=200.0, zr_b=1.5)
    # --- geometry and QC status tables ---
    geometry_runner.run_geometry(cfg, True, _Log(), [])
    rows = [{"cml_id": c, "sublink_id": s, "status": "good", "reason": "ok"} for c in IDS for s in schema.SUBLINK_IDS]
    pd.DataFrame(rows).to_csv(cfg.tables_dir / "sublink_qc_status.csv", index=False)
    return cfg


@pytest.fixture
def world(synthetic_cfg):
    """A fresh world per test, for tests that run stages again or change files."""
    return make_world(synthetic_cfg)


class _Log:
    def info(self, *a, **k):
        pass


def run_all(cfg, force=False, workers=1):
    rec = []
    aggregate_runner.run_aggregate_gauges(cfg, force, _Log(), rec)
    aggregate_runner.run_aggregate_radar(cfg, force, _Log(), rec)
    days = list(run_days(aggregate_runner.aggregate_cml_day, cfg, DAYS, workers, force=force))
    assert all(d["status"] != "failed" for d in days), [d.get("traceback") for d in days]
    return rec, days


def digest(paths):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def data_of(path):
    with xr.open_dataset(path) as ds:
        return ds.load()


def all_products(cfg):
    out = []
    for res in cfg["aggregate"]["resolutions"]:
        out += [aggregate_runner.cml_path(cfg, res, d) for d in DAYS]
    out += sorted(aggregate_runner.l3_dir(cfg, "l3_gauges_dir").glob("*.nc"))
    out += sorted(aggregate_runner.l3_dir(cfg, "l3_radar_dir").glob("*.nc"))
    return out


def test_world_products_axes_and_inputs_unchanged(aggregated_world):
    cfg, before, after = aggregated_world
    assert after == before  # L1 and L2 are only read

    for res in cfg["aggregate"]["resolutions"]:
        axes = [xr.concat([data_of(aggregate_runner.cml_path(cfg, res, d)) for d in DAYS], "time").time.values,
                data_of(aggregate_runner.gauge_path(cfg, "city", res)).time.values]
        if res in cfg["aggregate"]["gauges"]["smhi_resolutions"]:
            axes.append(data_of(aggregate_runner.gauge_path(cfg, "smhi", res)).time.values)
        if res in cfg["aggregate"]["radar"]["resolutions"]:
            axes += [data_of(aggregate_runner.radar_path(cfg, w, res)).time.values for w in ("along_links", "at_gauges")]
        else:
            assert not aggregate_runner.radar_path(cfg, "along_links", res).exists()  # no radar at 1 min
        assert all(np.array_equal(a, axes[0]) for a in axes), res
        assert len(axes[0]) == 2 * pd.Timedelta("1D") / pd.Timedelta(res)
    assert not aggregate_runner.gauge_path(cfg, "smhi", "5min").exists()


def test_world_radar(aggregated_world):
    cfg = aggregated_world[0]
    r5 = data_of(aggregate_runner.radar_path(cfg, "along_links", "5min"))
    assert r5.attrs["n_stamps_outside_period_dropped"] == 1
    assert float(r5.attrs["floor_share_of_cell_steps_set_to_zero"]) > 0.99
    spike = r5.R_mean.sel(time=SPIKE)
    assert np.allclose(spike, 10.0, rtol=1e-4)                       # uniform field -> R_mean = value
    assert np.allclose(r5.R_mean.sel(time=SPIKE - pd.Timedelta("5min")), 0.0)  # not in the bin before
    assert np.nanmax(r5.R_mean.drop_sel(time=SPIKE).values) == 0.0  # floor: 0, not 0.0003 phantom rain
    assert np.allclose(r5.wet_fraction.sel(time=SPIKE), 1.0)
    r15 = data_of(aggregate_runner.radar_path(cfg, "along_links", "15min"))
    assert np.allclose(r15.R_mean.sel(time="2015-07-28 06:00"), 10.0 / 3, rtol=1e-4)  # mean of 3 stamps
    missing = r5.sel(time="2015-07-29 12:00")
    assert np.isnan(missing.R_mean).all() and (missing.valid_weight_fraction == 0).all()
    assert (r15.sel(time="2015-07-29 12:00").n_valid_stamps == 2).all()  # 2 of 3 still enough
    g5 = data_of(aggregate_runner.radar_path(cfg, "at_gauges", "5min"))
    assert np.allclose(g5.R_pixel.sel(time=SPIKE), 10.0, rtol=1e-4) and np.allclose(g5.R_3x3.sel(time=SPIKE), 10.0, rtol=1e-4)


def test_world_gauges(aggregated_world):
    cfg = aggregated_world[0]
    with xr.open_dataset(cfg.l1_gauges_dir / "city_1min.nc") as l1:
        raw = l1.rainfall_amount.values[:, 0]
    for res in ("1min", "5min", "15min", "1h"):
        g = data_of(aggregate_runner.gauge_path(cfg, "city", res))
        total = float(g.amount_mm.sum())
        lost = float(np.nansum(raw)) - total
        # the bin holding the missing minute is NaN (strict); everything else is summed exactly
        assert lost >= -1e-4
        if res == "1min":
            assert total == pytest.approx(np.nansum(raw), abs=1e-3)
    s = data_of(aggregate_runner.gauge_path(cfg, "smhi", "1h"))
    assert "qc_worst" in s and int(s.qc_worst.max()) == 1


def test_world_skip_force_identical(world):
    cfg = world
    run_all(cfg)
    paths = all_products(cfg)
    first = {p: data_of(p) for p in paths}
    rec, days = run_all(cfg)  # second run: everything up to date
    assert [r["status"] for r in rec if r.get("status") == "skipped_same_config"] and \
        all(d["status"] == "skipped_same_config" for d in days)
    run_all(cfg, force=True)
    for p in paths:
        xr.testing.assert_identical(first[p].drop_attrs(), data_of(p).drop_attrs())


def test_world_workers_identical(world):
    cfg = world
    run_all(cfg, workers=1)
    paths = [aggregate_runner.cml_path(cfg, r, d) for r in cfg["aggregate"]["resolutions"] for d in DAYS]
    one = {p: data_of(p) for p in paths}
    days = list(run_days(aggregate_runner.aggregate_cml_day, cfg, DAYS, 4, force=True))
    assert all(d["status"] == "written" for d in days), [d.get("traceback") for d in days]
    for p in paths:
        xr.testing.assert_identical(one[p].drop_attrs(), data_of(p).drop_attrs())


def test_world_validate_files(aggregated_world):
    cfg = aggregated_world[0]
    for res in cfg["aggregate"]["resolutions"]:
        for d in DAYS:
            schema.validate_l3(data_of(aggregate_runner.cml_path(cfg, res, d)), "cml", res, start=d)
    for res in cfg["aggregate"]["radar"]["resolutions"]:
        schema.validate_l3(data_of(aggregate_runner.radar_path(cfg, "along_links", res)), "radar_links", res, start=DAY)
        schema.validate_l3(data_of(aggregate_runner.radar_path(cfg, "at_gauges", res)), "radar_gauges", res, start=DAY)


# ---------------------------------------------------------------------------
# access
# ---------------------------------------------------------------------------
def test_window_edges():
    assert window_edges("2015-07-28 10:07", "2015-07-28 11:00", "1h", "30min", "15min") == \
        (pd.Timestamp("2015-07-28 09:00"), pd.Timestamp("2015-07-28 11:30"))
    assert window_edges("2015-07-28 10:00", "2015-07-28 11:00") == \
        (pd.Timestamp("2015-07-28 10:00"), pd.Timestamp("2015-07-28 11:00"))
    with pytest.raises(ValueError):
        window_edges("2015-07-28 10:00", "2015-07-28 10:00")
    with pytest.raises(ValueError):
        window_edges("2015-07-28 10:00", "2015-07-28 11:00", pad_before="-1h")


def test_access(aggregated_world):
    cfg = aggregated_world[0]
    ds = load_link_window(cfg, 300, "5min", "2015-07-28 06:00", "2015-07-28 07:00", pad_before="1h",
                          pad_after="23h", max_gauge_dist_m=2000)
    assert ds.sizes["time"] == (1 + 1 + 23) * 12  # [05:00, 2015-07-29 06:00) crosses midnight
    assert all(v.dims[0] == "time" for v in ds.data_vars.values())
    assert list(ds.gauge.values) == ["G1"] and "smhi" in ds.attrs["sources_absent"]  # SMHI is 15 min native
    assert float(ds.radar_R_mean.sel(time=SPIKE)) == pytest.approx(10.0, rel=1e-4)
    assert float(ds.distance_to_path_m.sel(gauge="G1")) < 50
    q = load_link_window(cfg, 300, "15min", "2015-07-28 06:00", "2015-07-28 07:00", max_gauge_dist_m=2000)
    assert set(q.gauge.values) == {"G1", "G2"} and q.attrs["sources_absent"] == "none"
    assert q.sizes["time"] == 4 and np.isnan(q.gauge_qc_worst.sel(gauge="G1")).all()
    one = load_link_window(cfg, 300, "1min", "2015-07-28 06:00", "2015-07-28 07:00")
    assert not any(n.startswith("radar_") for n in one.data_vars) and "radar" in one.attrs["sources_absent"]
    assert "smhi" in one.attrs["sources_absent"] and list(one.gauge.values) == ["G1"]
    assert one.sizes["time"] == 60
    with pytest.raises(ValueError):
        load_link_window(cfg, 300, "2min", "2015-07-28 06:00", "2015-07-28 07:00")
    with pytest.raises(KeyError):
        load_link_window(cfg, 999, "5min", "2015-07-28 06:00", "2015-07-28 07:00")
    with pytest.raises(ValueError):
        load_link_window(cfg, 300, "5min", "2015-07-27 06:00", "2015-07-28 07:00")


# ---------------------------------------------------------------------------
# real data (skipped when the pilot outputs are absent)
# ---------------------------------------------------------------------------
@pytest.mark.data
def test_real_pilot_day_matches_plain_numpy():
    from cmlpipe.config import load_config
    from conftest import REAL_CONFIG
    from cmlpipe.qc import stage as qc_runner
    cfg = load_config(REAL_CONFIG)
    path = aggregate_runner.cml_path(cfg, "15min", DAY)
    if not (path.exists() and cfg.l1_cml_day_path(DAY).exists()):
        pytest.skip("real L3 pilot day not available")
    acfg = cfg["aggregate"]
    with xr.open_dataset(cfg.l1_cml_day_path(DAY)) as l1, xr.open_dataset(qc_runner.l2_day_path(cfg, DAY)) as l2, \
            xr.open_dataset(path) as l3:
        j = list(l1.cml_id.values).index(10227)
        rsl = l1.rsl.values[:, 0, j]
        fl = l2.flags.values[:, 0, j]
        mask = np.uint16(sum(1 << schema.FLAGS[n] for n in acfg["blocking_flags"]["rsl"]))
        v = np.where(np.isnan(rsl) | ((fl & mask) != 0), np.nan, rsl).reshape(96, 90)
        n = np.isfinite(v).sum(axis=1)
        with np.errstate(all="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                mean = np.where(n >= 72, np.nanmean(v, axis=1), np.nan)
        got = l3.rsl_mean.values[:, 0, j]
        assert np.array_equal(l3.n_valid_rsl.values[:, 0, j], n)
        assert np.allclose(got, mean, atol=1e-4, equal_nan=True)
