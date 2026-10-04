"""Step 2 tests: reference data (gauges, weather, radar).

Part A - invariants on tiny synthetic raw files that copy the real layouts
(see docs/methods.md). No real data needed. They contain on purpose:
  City : end-of-interval labels, an "NA", a duplicated timestamp with a different value,
         a missing minute, a row outside the period, a weighing and a tipping-bucket gauge
  SMHI : end-of-interval Pvol labels, QC letters G/Y and one missing code, hourly weather
         values only on the full hour, Latin-1 metadata with two text blocks
  Radar: int32 pseudo-dBZ with 255 = missing, scale_factor/add_offset, zr_a/zr_b as strings,
         a stamp exactly at the period end (instant: must be kept) and one after it (dropped)
  The synthetic period is one day, 2015-07-28, so the period end is 2015-07-29 00:00.

(The acceptance against the published OpenMRG numbers was a development check; it is archived.)
"""

import json

import netCDF4
import numpy as np
import pandas as pd
import pytest
import xarray as xr
import yaml

from cmlpipe.io import schema
from cmlpipe.io import openmrg_reference as ref
from cmlpipe.config import load_config
from conftest import REAL_CONFIG

DAY = "2015-07-28"
T0 = pd.Timestamp(DAY)

# ---------------------------------------------------------------------------
# Synthetic raw files
# ---------------------------------------------------------------------------
# City rows: raw END-of-interval labels (minutes after T0) and values for stations A (weighing), B (tipping).
CITY_ROWS = [
    (1, "0", "0"),
    (2, "0.3", "0.2"),
    (3, "NA", "0.4"),     # NA -> NaN
    (3, "0.9", "0.6"),    # duplicate of minute 3 with different values: dropped, counted
    (5, "1.1", "0"),      # minute 4 is missing -> inserted as NaN
    (1441, "0.5", "0.2"),  # 2015-07-29 00:01 end label -> starts 00:00 next day: outside the period
]
CITY_META = """Name,Location,Longitude_DecDeg,Latitude_DecDeg,Type,Unit,Resolution (mm)
A,Place A,11.9,57.7,Weighing,millimetre,0.1
B,Place B,12.0,57.6,Tipping-bucket,millimetre,0.2
"""
SMHI_META = (
    "Station information\r\nStation name,Göteborg A\r\nStation ID (klimatnummer),71420\r\n"
    "Station height above sea level (m),3.038\r\nLatitude (decimal degrees),57.7156\r\n"
    "Longitude (decimal degrees),11.9924\r\nQuality Codes (QC),G = approved. Y = suspicious\r\n\r\n"
    "Parameter information\r\n"
    "Short name,English name,Swedish name,Unit,Measurement height above ground (m),Description\r\n"
    "Pvol,Precipitation volume,Nederbörd,millimetre,2,Precipitation sum over preceeding 15min period\r\n"
    "Ta,Air temperature,Lufttemperatur,degree celsius,2,Momentaneous value once per hour\r\n"
    "RH,Relative humidity,Fukt,percent,2,Momentaneous value once per hour\r\n"
    "Pmsl,Air pressure,Tryck,hectopascal,0,Momentaneous value once per hour at sea level\r\n"
    "Vmean,Wind speed,Vind,metre per second,10,Average value during the preceeding 10 minutes\r\n"
    "Vdir,Wind direction,Vindriktning,degree true,10,Average value during the preceeding 10 minutes\r\n"
    "Vg,Gust wind speed,Byvind,metre per second,10,Maximum value during the preceeding hour\r\n"
)
# minutes after T0 (end labels); 1440 = the period end 2015-07-29 00:00
SMHI_PVOL = [(15, "0.4", "G"), (30, "1.2", "Y"), (45, "0", "G"), (60, "0.1", "NA"), (1440, "0.2", "G")]
SMHI_MET_MINUTES = (60, 1440)  # weather values only on these (full-hour) rows
RADAR_RAW = np.array([[[0, 255], [80, 105]], [[255, 255], [255, 255]], [[222, 50], [75, 0]],
                      [[10, 20], [30, 40]], [[1, 1], [1, 1]]], dtype="int32")
# 00:15 missing -> inserted; 1440 = period end (kept, instant); 1445 = after the period (dropped)
RADAR_MINUTES = [5, 10, 20, 1440, 1445]


def _utc(minutes: int) -> str:
    return (T0 + pd.Timedelta(minutes=minutes)).strftime("%Y-%m-%dT%H:%M:%SZ")


def write_reference_raw(folder):
    folder.mkdir(parents=True, exist_ok=True)
    lines = ["Time_UTC,A,B"] + [f"{_utc(m)},{a},{b}" for m, a, b in CITY_ROWS]
    (folder / "city.csv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (folder / "city_meta.csv").write_text(CITY_META, encoding="utf-8-sig")

    header = ("Time_UTC,Pvol_mm,Pvol_QC,Ta_C,Ta_QC,RH_%,RH_QC,Pmsl_hPa,Pmsl_QC,"
              "Vmean_m/s,Vmean_QC,Vdir_deg,Vdir_QC,Vg_m/s,Vg_QC")
    rows = []
    for m, p, qc in SMHI_PVOL:
        met = "15.2,G,80,Y,1012.3,G,3.1,G,200,G,7.5,G" if m in SMHI_MET_MINUTES else ",".join(["NA"] * 12)
        rows.append(f"{_utc(m)},{p},{qc},{met}")
    (folder / "smhi.csv").write_text("\n".join([header] + rows) + "\n", encoding="utf-8")
    (folder / "smhi_meta.csv").write_bytes(SMHI_META.encode("latin-1"))

    with netCDF4.Dataset(folder / "radar.nc", "w", format="NETCDF4") as nc:
        nc.createDimension("time", len(RADAR_MINUTES))
        nc.createDimension("y", 2)
        nc.createDimension("x", 2)
        tv = nc.createVariable("time", "i4", ("time",))
        tv.setncatts({"units": "seconds since 1970-01-01 00:00:00", "calendar": "gregorian"})
        tv[:] = [int((T0 + pd.Timedelta(minutes=m)).value // 10**9) for m in RADAR_MINUTES]
        nc.createVariable("x", "f8", ("x",))[:] = [-120000.0, -118000.0]
        nc.createVariable("y", "f8", ("y",))[:] = [-3450000.0, -3452000.0]
        nc.createVariable("lat", "f4", ("y", "x"))[:] = [[57.70, 57.70], [57.68, 57.68]]
        nc.createVariable("lon", "f4", ("y", "x"))[:] = [[11.90, 11.93], [11.90, 11.93]]
        crs = nc.createVariable("crs", "i4")
        crs.setncatts({"grid_mapping_name": "polar_stereographic", "spatial_ref": "test-wkt"})
        d = nc.createVariable("data", "i4", ("time", "y", "x"), fill_value=np.int32(255))
        d.setncatts({"missing_value": np.int32(255), "scale_factor": np.float32(0.4),
                     "add_offset": np.float32(-30.0), "zr_a": "200", "zr_b": "1.5", "grid_mapping": "crs"})
        d.set_auto_maskandscale(False)  # write the raw integers, like the real file
        d[:] = RADAR_RAW
        nc.setncatts({"proj_string": "+proj=stere +lat_ts=60 +ellps=bessel +lon_0=14 +lat_0=90",
                      "doi": "test-doi", "version": "test"})


@pytest.fixture
def ref_cfg(tmp_path):
    """Real config with raw paths pointing to the synthetic files and a one-day period."""
    raw = tmp_path / "raw"
    write_reference_raw(raw)
    data = yaml.safe_load(REAL_CONFIG.read_text(encoding="utf-8"))
    data["paths"].update({
        "raw_city_gauges": str(raw / "city.csv"), "raw_city_gauges_metadata": str(raw / "city_meta.csv"),
        "raw_smhi": str(raw / "smhi.csv"), "raw_smhi_metadata": str(raw / "smhi_meta.csv"),
        "raw_radar": str(raw / "radar.nc"), "output_root": str(tmp_path / "out"),
    })
    data["period"].update({"start": DAY, "end": DAY})
    path = tmp_path / "ref_config.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return load_config(path, repo_root=tmp_path)


# ---------------------------------------------------------------------------
# Part A: invariants
# ---------------------------------------------------------------------------
def test_city_time_shift_and_sum(ref_cfg):
    ds = ref.load_city_gauges(ref_cfg)
    a = ds.rainfall_amount.sel(id="A").to_series()
    # raw end label 00:02 -> interval start 00:01, etc. (shift -60 s)
    assert ds.attrs["time_shift_applied_s"] == -60 and ds.attrs["time_convention"] == "interval_start"
    assert a[T0 + pd.Timedelta(minutes=1)] == pytest.approx(0.3)
    assert a[T0 + pd.Timedelta(minutes=4)] == pytest.approx(1.1)
    # grid only between the first and last raw stamp: 00:00 .. 00:04 (raw end labels 00:01 .. 00:05)
    assert pd.Timestamp(ds.time.values[0]) == T0 and len(ds.time) == 5
    # sum after the shift = raw sum of the rows that stay (inside the period, first copy of the duplicate)
    raw_kept = [(0, 0), (0.3, 0.2), (np.nan, 0.4), (1.1, 0)]
    assert float(ds.rainfall_amount.sum()) == pytest.approx(np.nansum(raw_kept), abs=1e-5)


def test_city_na_duplicate_and_inserted_counted(ref_cfg):
    ds = ref.load_city_gauges(ref_cfg)
    assert ds.attrs["n_na_converted"] == 1
    assert np.isnan(float(ds.rainfall_amount.sel(id="A", time=T0 + pd.Timedelta(minutes=2))))  # NA, first copy kept
    assert ds.attrs["n_duplicates_removed"] == 1 and ds.attrs["n_duplicates_with_different_values"] == 1
    assert ds.attrs["n_outside_period_dropped"] == 1
    assert ds.attrs["n_steps_inserted"] == 1  # minute 00:03, between real stamps
    assert list(ds.gauge_type.values) == ["weighing", "tipping_bucket"]
    assert list(ds.quantization_mm.values) == [0.1, 0.2]
    schema.validate_gauge_l1(ds, interval="1min")


def test_smhi_shift_sum_and_qc(ref_cfg):
    ds = ref.load_smhi_gauge(ref_cfg)
    s = ds.rainfall_amount.sel(id="GbgA").to_series()
    assert ds.attrs["time_shift_applied_s"] == -900
    # raw end labels 00:15, 00:30, 00:45, 01:00 -> starts 00:00, 00:15, 00:30, 00:45
    assert s.iloc[:4].tolist() == pytest.approx([0.4, 1.2, 0.0, 0.1])
    assert pd.Timestamp(ds.time.values[-1]) == T0 + pd.Timedelta(minutes=1425)  # end label 24:00 -> start 23:45
    assert float(s.sum()) == pytest.approx(sum(float(p) for _, p, _ in SMHI_PVOL), abs=1e-5)
    assert ds.qc_code.values[:4, 0].tolist() == [1, 2, 1, 0]
    assert float(ds.height_m.values[0]) == 2.0
    schema.validate_gauge_l1(ds, interval="15min")


def test_weather_kept_as_in_file(ref_cfg):
    ds = ref.load_smhi_weather(ref_cfg)
    assert ds.attrs["time_convention"] == "as_in_file" and ds.attrs["time_shift_applied_s"] == 0
    assert float(ds.temperature.sel(id="GbgA", time=T0 + pd.Timedelta(hours=1))) == pytest.approx(15.2)
    assert int(ds.qc_humidity.sel(id="GbgA", time=T0 + pd.Timedelta(hours=1))) == 2  # "Y"
    assert "10 min" in ds.wind_speed.attrs["time_meaning_original"]
    schema.validate_met_l1(ds, interval="1h")


def test_radar_decoded_once_and_zr(ref_cfg):
    ds = ref.load_radar_month(ref_cfg, 2015, 7)
    assert ds.dbz.dims == ("time", "y", "x")
    t = pd.Timestamp(T0 + pd.Timedelta(minutes=5))
    dbz = ds.dbz.sel(time=t).values
    assert dbz[0, 0] == pytest.approx(-30.0) and np.isnan(dbz[0, 1]) and dbz[1, 0] == pytest.approx(2.0)
    # the same numbers as xarray's own automatic decoding (so scaling was not applied twice)
    with xr.open_dataset(ref_cfg.raw_path("raw_radar")) as auto:
        common = pd.DatetimeIndex(ds.time.values).intersection(pd.DatetimeIndex(auto.time.values))
        np.testing.assert_array_equal(ds.dbz.sel(time=common).values,
                                      auto["data"].sel(time=common).values.astype("float32"))
    expected = ((10 ** (2.0 / 10)) / 200.0) ** (1 / 1.5)
    assert float(ds.rain_rate.sel(time=t).values[1, 0]) == pytest.approx(expected, rel=1e-6)
    assert ds.attrs["zr_a"] == 200.0 and ds.attrs["zr_b"] == 1.5
    assert ds.attrs["n_missing_cells_raw"] == 5 and ds.attrs["n_steps_all_missing_raw"] == 1
    assert np.isnan(ds.dbz.sel(time=T0 + pd.Timedelta(minutes=15))).all()  # inserted step (between real stamps)
    assert ds.attrs["n_steps_inserted"] == 288 - 4
    assert np.isfinite(ds.x_3006.values).all()
    schema.validate_radar_l1(ds, interval="5min")


def _raw_stamps(minutes) -> pd.DatetimeIndex:
    return pd.DatetimeIndex([T0 + pd.Timedelta(minutes=m) for m in minutes])


@pytest.mark.parametrize("source", ["radar", "weather"])
def test_time_labels_instant_sources(ref_cfg, source):
    """Every raw stamp inside the window [T0, T0 + 1 day] (closed: instant stamps) is in L1, the closing
    stamp included, and L1 has no stamp before the first or after the last raw stamp."""
    if source == "radar":
        ds = ref.load_radar_month(ref_cfg, 2015, 7)
        raw = _raw_stamps(RADAR_MINUTES)
    else:
        ds = ref.load_smhi_weather(ref_cfg)
        raw = _raw_stamps(SMHI_MET_MINUTES)
    l1 = pd.DatetimeIndex(ds.time.values)
    window_end = T0 + pd.Timedelta("1D")
    inside = raw[(raw >= T0) & (raw <= window_end)]
    assert inside.difference(l1).empty                       # nothing real is lost ...
    assert window_end in l1                                   # ... including the period's closing instant
    assert l1.min() == inside.min() and l1.max() == inside.max()  # nothing invented at the edges
    if source == "radar":
        assert ref.radar_steps_outside_months(ref_cfg) == 1   # the 00:05 stamp after the period
        assert "last month file" in ds.attrs["period_end_stamp"]
        assert ds.attrs["time_window"].endswith("]")


def test_validator_rejects_tz_aware_time(ref_cfg):
    ds = ref.load_smhi_gauge(ref_cfg)
    aware = pd.DatetimeIndex(ds.time.values).tz_localize("UTC")
    bad = ds.assign_coords(time=np.array(aware.to_pydatetime(), dtype=object))
    with pytest.raises(schema.CMLValidationError, match="tz-naive"):
        schema.validate_gauge_l1(bad, interval="15min")


def test_validator_rejects_wrong_units(ref_cfg):
    ds = ref.load_city_gauges(ref_cfg)
    ds["rainfall_amount"].attrs["units"] = "mm h-1"
    with pytest.raises(schema.CMLValidationError, match="units"):
        schema.validate_gauge_l1(ds, interval="1min")


def test_validator_rejects_radar_dim_order(ref_cfg):
    ds = ref.load_radar_month(ref_cfg, 2015, 7).transpose("time", "x", "y")
    with pytest.raises(schema.CMLValidationError, match="dims"):
        schema.validate_radar_l1(ds, interval="5min")


def test_cli_reference_steps(ref_cfg):
    """run.py writes all reference files, validates them, skips them on the second run."""
    from cmlpipe import cli as run

    args = ["--config", str(ref_cfg.config_file), "--steps", "ingest_gauges,ingest_met,ingest_radar", "--all"]
    assert run.main(args) == 0
    assert run.main(args) == 0
    runs = json.loads(ref_cfg.manifest_path.read_text(encoding="utf-8"))
    first = [r for r in runs[0]["files"] if r["status"] == "written"]
    assert {r["path"].split("/")[-1] for r in first} >= {"city_1min.nc", "smhi_15min.nc", "smhi_met_1h.nc",
                                                          "radar_5min_2015-07.nc", "gauges.csv"}
    city = next(r for r in first if r["path"].endswith("city_1min.nc"))
    assert city["shape"] == {"time": 5, "id": 2} and city["n_duplicates_removed"] == 1
    assert all(r["status"] in ("skipped_same_config", "note") for r in runs[1]["files"])
    table = pd.read_csv(ref_cfg.tables_dir / "gauges.csv")
    assert list(table["id"]) == ["A", "B", "GbgA"]
    assert (ref_cfg.tables_dir / "gauges.csv.meta.json").exists()
    with xr.open_dataset(ref_cfg.l1_gauges_dir / "smhi_15min.nc") as ds:
        assert ds.qc_code.dtype == np.uint8 and ds.attrs["config_hash"] == ref_cfg.hash
