"""OpenMRG reference adapter: City gauges, SMHI gauge, SMHI weather and radar -> L1 files.

Raw facts behind every choice here: docs/data_notes/openmrg_reference_inspection.md.

Allowed transforms (and nothing else):
  - NA token / 255 fill value -> NaN (counted)
  - timestamp convention shift to "interval start" where the raw meaning is documented
    (City: -60 s, SMHI Pvol: -900 s; weather and radar are kept as in the file)
  - duplicates dropped (first kept), off-grid or out-of-period steps dropped, missing
    steps inserted as NaN; every one of these is counted and stored in the attributes
  - reshaping, unit bookkeeping, radar pseudo-dBZ -> dBZ with the file's own
    scale_factor/add_offset, and rain_rate from dBZ with the file's own Z-R (a, b)
No smoothing, filling, correction or thresholding of measured values. Negative gauge
values would be kept and only counted.
"""

from __future__ import annotations

import csv
import io
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pyproj
import xarray as xr

from cmlpipe.io import schema
from cmlpipe.config import Config
from cmlpipe.provenance import provenance_attrs, write_sidecar

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def _to_numeric(tokens: pd.Series, na_tokens: list[str], what: str) -> tuple[np.ndarray, int]:
    """Text column -> float64 array. Tokens in `na_tokens` become NaN (and are counted).

    Any other token that is not a number raises, so a new spelling of "missing" in a
    future file version is noticed instead of silently becoming NaN.
    """
    tokens = tokens.astype(str).str.strip()
    is_na = tokens.isin(na_tokens)
    values = pd.to_numeric(tokens.where(~is_na), errors="coerce")
    bad = tokens[values.isna() & ~is_na]
    if len(bad):
        raise ValueError(f"{what}: non-numeric tokens {bad.value_counts().head(5).to_dict()}; "
                         f"add them to na_tokens in the config if they mean 'missing'")
    return values.to_numpy(dtype="float64"), int(is_na.sum())


def _parse_utc(text: pd.Series) -> pd.DatetimeIndex:
    """ISO 8601 strings with 'Z' -> tz-naive UTC DatetimeIndex in nanoseconds."""
    return pd.DatetimeIndex(pd.to_datetime(text, utc=True)).tz_convert(None).as_unit("ns")


def period_window(cfg: Config, stamp_kind: str) -> tuple[pd.Timestamp, pd.Timestamp, bool]:
    """(start, end, end_inclusive) of the period for a kind of time stamp (config period.window_end_inclusive).

    interval_start stamps: [start, end)  - a value stamped `end` would describe time after the period.
    instant stamps:        [start, end]  - a measurement taken exactly at `end` belongs to the period.
    """
    inclusive = bool(cfg["period"]["window_end_inclusive"][stamp_kind])
    return cfg.period_start, cfg.period_end, inclusive


def window_text(start: pd.Timestamp, end: pd.Timestamp, end_inclusive: bool) -> str:
    """Text form of a period window, e.g. '[start, end)'."""
    return f"[{start:%Y-%m-%d %H:%M}, {end:%Y-%m-%d %H:%M}{']' if end_inclusive else ')'}"


def _regularize(times: pd.DatetimeIndex, arrays: dict[str, np.ndarray], start: pd.Timestamp,
                end: pd.Timestamp, interval: pd.Timedelta, fill: dict | None = None,
                end_inclusive: bool = False) -> tuple[pd.DatetimeIndex, dict[str, np.ndarray], dict]:
    """Put rows (already shifted to the pipeline convention) onto a regular grid inside the window.

    The window is [start, end) or, with end_inclusive, [start, end].
    Steps, each counted:
      1. duplicated timestamps: keep the first; also count duplicates whose values differ
      2. rows outside the window: dropped
      3. rows not on the grid (start + k * interval): dropped (none expected in OpenMRG)
      4. the grid runs from the FIRST to the LAST kept raw stamp. Steps between them without a
         row are inserted with NaN (or fill[name] for integer arrays). No step is ever created
         before the first or after the last real stamp: such a step would be invented, not missing.

    arrays: name -> array with time as the first axis.
    Returns (grid, gridded arrays, counts).
    """
    fill = fill or {}
    times = pd.DatetimeIndex(times).as_unit("ns")
    counts = {"n_raw_steps": int(len(times))}

    is_dup = times.duplicated(keep="first")
    n_differ = 0
    if is_dup.any():
        first_pos = pd.Series(np.arange(len(times))).groupby(times).transform("first").to_numpy()
        for a in arrays.values():
            x, y = a[is_dup], a[first_pos[is_dup]]
            same = (x == y) | (pd.isna(x) & pd.isna(y))
            n_differ += int((~same.reshape(len(x), -1).all(axis=1)).sum())
    counts["n_duplicates_removed"] = int(is_dup.sum())
    counts["n_duplicates_with_different_values"] = n_differ
    keep = ~is_dup

    inside = (times >= start) & ((times <= end) if end_inclusive else (times < end))
    counts["n_outside_period_dropped"] = int((keep & ~inside).sum())
    keep &= inside

    offset = (times - start).asi8  # ns since start
    on_grid = offset % interval.value == 0
    counts["n_offgrid_dropped"] = int((keep & ~on_grid).sum())
    keep &= on_grid

    if keep.any():
        first, last = times[keep].min(), times[keep].max()
        grid = pd.date_range(first, last, freq=interval).as_unit("ns")  # both ends included
        pos = (times[keep] - first).asi8 // interval.value
    else:
        grid = pd.DatetimeIndex([], dtype="datetime64[ns]")
        pos = np.array([], dtype="int64")
    out = {}
    for name, a in arrays.items():
        full = np.full((len(grid),) + a.shape[1:], fill.get(name, np.nan), dtype=a.dtype)
        full[pos] = a[keep]
        out[name] = full
    counts["n_steps_inserted"] = int(len(grid) - keep.sum())
    return grid, out, counts


def lonlat_to_metric(cfg: Config, lon, lat) -> tuple[np.ndarray, np.ndarray]:
    """WGS84 lon/lat (degrees) -> x, y in metres of the metric CRS (EPSG:3006 by default)."""
    tr = pyproj.Transformer.from_crs(cfg["reference"]["lonlat_crs"], cfg["reference"]["metric_crs"], always_xy=True)
    x, y = tr.transform(np.asarray(lon, dtype="float64"), np.asarray(lat, dtype="float64"))
    return np.asarray(x), np.asarray(y)


def _time_attrs(convention: str, interval: str, meaning: str, shift_s: int, window: str) -> dict:
    return {
        "time_window": window,
        "time_convention": convention,
        "interval_s": int(pd.Timedelta(interval).total_seconds()),
        "time_meaning_original": " ".join(str(meaning).split()),
        "time_shift_applied_s": int(shift_s),
        "time_zone": "UTC",
    }


def _set_attrs(ds: xr.Dataset, data_vars: dict, coords: dict) -> None:
    """Units / long_name from the schema tables onto the variables of `ds`."""
    for name, spec in data_vars.items():
        if name in ds:
            ds[name].attrs.update({"units": spec[0], "long_name": spec[1]})
    for name, (_, units, long_name) in coords.items():
        if name in ds.coords:
            ds[name].attrs.update({"long_name": long_name} if units is None else {"units": units, "long_name": long_name})
    ds["time"].attrs.update({"long_name": "time (UTC)", "standard_name": "time"})


def _base_attrs(cfg: Config, title: str, raw_files: list[str]) -> dict:
    return {
        "title": f"{cfg['dataset']['title']} - {title}",
        "dataset": cfg.dataset_name,
        "processing_level": "L1_standard",
        "source": cfg["dataset"]["source"],
        "source_doi": cfg["dataset"]["source_doi"],
        "licence": cfg["dataset"]["licence"],
        "raw_files": ", ".join(Path(cfg["paths"][k]).as_posix() for k in raw_files),
        "period_utc": f"[{cfg.period_start:%Y-%m-%d %H:%M}, {cfg.period_end:%Y-%m-%d %H:%M})",
    }


# ---------------------------------------------------------------------------
# City gauges
# ---------------------------------------------------------------------------
def load_city_gauges(cfg: Config) -> xr.Dataset:
    """City network (10 gauges, 1 min) in the gauge L1 schema.

    Raw Time_UTC is the END of each 1 min interval, so time is shifted by
    city_gauges.time_shift_s (-60 s) to label each value by its interval start.
    Station metadata (lon/lat, type, resolution) come from CityGauges-metadata.csv;
    the height (1.5 m, "typically") comes from the config because no file states it.
    """
    c = cfg["city_gauges"]
    interval = pd.Timedelta(c["interval"])
    raw = pd.read_csv(cfg.raw_path("raw_city_gauges"), dtype=str, keep_default_na=False)
    meta = pd.read_csv(cfg.raw_path("raw_city_gauges_metadata"), encoding="utf-8-sig")  # file has a BOM

    stations = [s for s in raw.columns if s != "Time_UTC"]
    not_in_meta = sorted(set(stations) - set(meta["Name"]))
    if not_in_meta:
        raise ValueError(f"City stations without metadata: {not_in_meta}")
    only_in_meta = sorted(set(meta["Name"]) - set(stations))
    if only_in_meta:
        log.warning("City stations in the metadata but not in the values file (ignored): %s", only_in_meta)

    values = np.empty((len(raw), len(stations)), dtype="float64")
    n_na = 0
    for j, s in enumerate(stations):
        values[:, j], n = _to_numeric(raw[s], c["na_tokens"], f"City {s}")
        n_na += n
    raw_sum = float(np.nansum(values))

    times = _parse_utc(raw["Time_UTC"]) + pd.Timedelta(seconds=int(c["time_shift_s"]))
    start, end, incl = period_window(cfg, c["stamp_kind"])
    grid, arr, counts = _regularize(times, {"a": values}, start, end, interval, end_inclusive=incl)
    amount = arr["a"].astype("float32")

    m = meta.set_index("Name").loc[stations]
    type_key = m["Type"].astype(str).str.strip().str.lower()
    unknown = sorted(set(type_key) - set(c["gauge_type_map"]))
    if unknown:
        raise ValueError(f"Unknown City gauge types {unknown}; add them to city_gauges.gauge_type_map")
    lon = m["Longitude_DecDeg"].to_numpy(float)
    lat = m["Latitude_DecDeg"].to_numpy(float)
    x, y = lonlat_to_metric(cfg, lon, lat)

    ds = xr.Dataset(
        {"rainfall_amount": (schema.GAUGE_DIMS, amount)},
        coords={
            "time": grid.to_numpy(),
            "id": np.array(stations, dtype=str),
            "lon": ("id", lon), "lat": ("id", lat), "x": ("id", x), "y": ("id", y),
            "network": ("id", np.array(["city"] * len(stations))),
            "gauge_type": ("id", type_key.map(c["gauge_type_map"]).to_numpy(str)),
            "quantization_mm": ("id", m["Resolution (mm)"].to_numpy(float)),
            "height_m": ("id", np.full(len(stations), float(c["height_m"]))),
            "location": ("id", m["Location"].to_numpy(str)),
        },
    )
    _set_attrs(ds, schema.GAUGE_DATA_VARS, schema.GAUGE_COORDS)
    ds["location"].attrs["long_name"] = "station location name (raw metadata 'Location')"
    ds["height_m"].attrs["source"] = " ".join(c["height_source"].split())
    ds.attrs = {
        **_base_attrs(cfg, "L1 City gauges, 1 min", ["raw_city_gauges", "raw_city_gauges_metadata"]),
        **_time_attrs("interval_start", c["interval"], c["time_meaning_original"], c["time_shift_s"],
                      window_text(start, end, incl)),
        "network": "city",
        "units_original": "millimetre per interval (metadata column Unit)",
        "transforms": "NA -> NaN; time shifted to interval start; reindexed to regular grid (inserted steps NaN).",
        "n_na_converted": n_na,
        "n_negative": int(np.sum(amount < 0)),
        "raw_sum_mm": round(raw_sum, 6),
        **counts,
    }
    return ds


# ---------------------------------------------------------------------------
# SMHI: metadata, gauge, weather
# ---------------------------------------------------------------------------
def read_smhi_metadata(cfg: Config) -> dict:
    """Parse GbgA-71420-metadata.csv (Latin-1, two blocks of text, not one table).

    Returns {"station": {key: value}, "parameters": DataFrame indexed by Short name}.
    """
    text = cfg.raw_path("raw_smhi_metadata").read_text(encoding="latin-1")
    station_block, _, param_block = text.partition("Parameter information")
    station = {}
    for row in csv.reader(io.StringIO(station_block)):
        if len(row) >= 2 and row[0] != "Station information":
            station[row[0].strip()] = ",".join(row[1:]).strip()
    params = pd.read_csv(io.StringIO(param_block.strip())).set_index("Short name")
    return {"station": station, "parameters": params}


def _smhi_station_coords(cfg: Config, meta: dict) -> dict:
    st = meta["station"]
    lat = float(st["Latitude (decimal degrees)"])
    lon = float(st["Longitude (decimal degrees)"])
    x, y = lonlat_to_metric(cfg, [lon], [lat])
    return {"id": np.array([cfg["smhi_gauge"]["id"]]), "lon": ("id", [lon]), "lat": ("id", [lat]),
            "x": ("id", x), "y": ("id", y)}


def _qc_codes(tokens: pd.Series, code_map: dict, what: str) -> np.ndarray:
    """SMHI QC letters -> uint8 (0 = no code). An unknown letter raises."""
    t = tokens.astype(str).str.strip()
    t = t.where(~t.isin(["NA", ""]), "")
    unknown = sorted(set(t) - set(code_map) - {""})
    if unknown:
        raise ValueError(f"{what}: unknown QC codes {unknown}; add them to the config")
    return t.map({**code_map, "": 0}).to_numpy(dtype="uint8")


def _qc_attrs(cfg_section: dict) -> dict:
    codes = cfg_section["qc_codes"]
    return {
        "long_name": "data-provider (SMHI) quality code",
        "flag_values": np.array([0] + list(codes.values()), dtype="uint8"),
        "flag_meanings": " ".join(["no_code"] + [f"{k}" for k in codes]),
        "comment": "; ".join(f"{v} = {k}: {cfg_section['qc_meanings'][k]}" for k, v in codes.items())
                   + "; 0 = no code in the raw file",
    }


def load_smhi_gauge(cfg: Config) -> xr.Dataset:
    """SMHI GbgA precipitation (Pvol, 15 min) in the gauge L1 schema, with the SMHI QC code kept.

    Raw Time_UTC is the END of each 15 min interval, so time is shifted by -900 s.
    Gauge type and resolution are not in the SMHI metadata file; they come from the config
    (Andersson et al. 2022, Sect. 2.2). The height is the metadata's "Measurement height
    above ground" for Pvol (NaN if absent).
    """
    c = cfg["smhi_gauge"]
    interval = pd.Timedelta(c["interval"])
    raw = pd.read_csv(cfg.raw_path("raw_smhi"), dtype=str, keep_default_na=False)
    meta = read_smhi_metadata(cfg)

    amount, n_na = _to_numeric(raw["Pvol_mm"], c["na_tokens"], "SMHI Pvol_mm")
    qc = _qc_codes(raw["Pvol_QC"], c["qc_codes"], "SMHI Pvol_QC")
    raw_sum = float(np.nansum(amount))
    times = _parse_utc(raw["Time_UTC"]) + pd.Timedelta(seconds=int(c["time_shift_s"]))
    start, end, incl = period_window(cfg, c["stamp_kind"])
    grid, arr, counts = _regularize(times, {"a": amount[:, None], "qc": qc[:, None]},
                                    start, end, interval, fill={"qc": 0}, end_inclusive=incl)
    amount = arr["a"].astype("float32")

    height = np.nan
    params = meta["parameters"]
    if "Pvol" in params.index and "Measurement height above ground (m)" in params.columns:
        height = float(params.loc["Pvol", "Measurement height above ground (m)"])

    ds = xr.Dataset(
        {"rainfall_amount": (schema.GAUGE_DIMS, amount),
         "qc_code": (schema.GAUGE_DIMS, arr["qc"].astype("uint8"))},
        coords={
            "time": grid.to_numpy(),
            **_smhi_station_coords(cfg, meta),
            "network": ("id", ["smhi"]),
            "gauge_type": ("id", [c["gauge_type"]]),
            "quantization_mm": ("id", [float(c["quantization_mm"])]),
            "height_m": ("id", [height]),
        },
    )
    _set_attrs(ds, schema.GAUGE_DATA_VARS, schema.GAUGE_COORDS)
    ds["qc_code"].attrs = _qc_attrs(c)
    ds["height_m"].attrs["source"] = "GbgA-71420-metadata.csv, Pvol 'Measurement height above ground (m)'"
    st = meta["station"]
    ds.attrs = {
        **_base_attrs(cfg, "L1 SMHI gauge GbgA, 15 min", ["raw_smhi", "raw_smhi_metadata"]),
        **_time_attrs("interval_start", c["interval"], c["time_meaning_original"], c["time_shift_s"],
                      window_text(start, end, incl)),
        "network": "smhi",
        "station_name": st.get("Station name", ""),
        "station_number": st.get("Station ID (klimatnummer)", ""),
        "station_height_asl_m": st.get("Station height above sea level (m)", ""),
        "units_original": str(params.loc["Pvol", "Unit"]) if "Pvol" in params.index else "",
        "transforms": "NA -> NaN; QC letter -> uint8; time shifted to interval start; reindexed to regular grid.",
        "n_na_converted": n_na,
        "n_negative": int(np.sum(amount < 0)),
        "raw_sum_mm": round(raw_sum, 6),
        **counts,
    }
    return ds


def load_smhi_weather(cfg: Config) -> xr.Dataset:
    """SMHI GbgA hourly weather (temperature, humidity, pressure, wind) in the weather L1 schema.

    Values exist only on the full hour; the 15 min rows between them are NA in the raw file
    and are dropped as off-grid rows (after checking that they hold no value). Time is kept
    as in the file (time_convention "as_in_file") because the variables have different
    meanings: instantaneous (Ta, RH, Pmsl), 10 min mean before the stamp (wind), 1 h maximum
    before the stamp (gust). Each variable carries its own time_meaning_original. The SMHI
    QC letters are kept as qc_<name> (uint8).
    """
    c = cfg["smhi_met"]
    interval = pd.Timedelta(c["interval"])
    raw = pd.read_csv(cfg.raw_path("raw_smhi"), dtype=str, keep_default_na=False)
    meta = read_smhi_metadata(cfg)
    na_tokens = cfg["smhi_gauge"]["na_tokens"]
    times = _parse_utc(raw["Time_UTC"]) + pd.Timedelta(seconds=int(c["time_shift_s"]))

    arrays, n_na = {}, 0
    for col, spec in c["variables"].items():
        v, n = _to_numeric(raw[col], na_tokens, f"SMHI {col}")
        arrays[spec["name"]] = v[:, None]
        arrays["qc_" + spec["name"]] = _qc_codes(raw[spec["qc"]], cfg["smhi_gauge"]["qc_codes"], f"SMHI {spec['qc']}")[:, None]

    # The non-hourly rows must be empty; if one held a value, dropping it would lose data.
    offgrid = ((times - cfg.period_start).asi8 % interval.value) != 0
    for name, a in arrays.items():
        held = np.isfinite(a[offgrid]) if a.dtype.kind == "f" else a[offgrid] != 0
        if held.any():
            raise ValueError(f"SMHI weather: {int(held.sum())} {name} values at non-hourly times")
    hourly = ~offgrid
    for name, a in arrays.items():
        if a.dtype.kind == "f":
            n_na += int(np.isnan(a[hourly]).sum())
    fill = {k: 0 for k in arrays if k.startswith("qc_")}
    start, end, incl = period_window(cfg, c["stamp_kind"])
    grid, arr, counts = _regularize(times[hourly], {k: a[hourly] for k, a in arrays.items()},
                                    start, end, interval, fill=fill, end_inclusive=incl)
    counts["n_nonhourly_rows_skipped"] = int(offgrid.sum())

    data_vars = {}
    for col, spec in c["variables"].items():
        name = spec["name"]
        data_vars[name] = (schema.MET_DIMS, arr[name].astype("float32"))
        data_vars["qc_" + name] = (schema.MET_DIMS, arr["qc_" + name].astype("uint8"))
    ds = xr.Dataset(data_vars, coords={"time": grid.to_numpy(), **_smhi_station_coords(cfg, meta)})
    _set_attrs(ds, schema.MET_DATA_VARS, schema.MET_COORDS)
    params = meta["parameters"]
    for col, spec in c["variables"].items():
        name, short = spec["name"], col.split("_")[0]
        ds[name].attrs["time_meaning_original"] = spec["meaning"]
        ds[name].attrs["raw_column"] = col
        if short in params.index:
            ds[name].attrs["units_original"] = str(params.loc[short, "Unit"])
            ds[name].attrs["description_original"] = str(params.loc[short, "Description"])
            ds[name].attrs["height_above_ground_m"] = float(params.loc[short, "Measurement height above ground (m)"])
        ds["qc_" + name].attrs = _qc_attrs(cfg["smhi_gauge"])
    ds.attrs = {
        **_base_attrs(cfg, "L1 SMHI weather GbgA, hourly", ["raw_smhi", "raw_smhi_metadata"]),
        **_time_attrs("as_in_file", c["interval"],
                      "per variable, see each variable's time_meaning_original (GbgA-71420-metadata.csv)",
                      c["time_shift_s"], window_text(start, end, incl)),
        "station_name": meta["station"].get("Station name", ""),
        "station_number": meta["station"].get("Station ID (klimatnummer)", ""),
        "transforms": "NA -> NaN; QC letter -> uint8; non-hourly (empty) rows skipped; reindexed to the hourly grid between the first and last raw stamp.",
        "n_na_converted": n_na,
        **counts,
    }
    return ds


# ---------------------------------------------------------------------------
# Radar
# ---------------------------------------------------------------------------
def month_window(cfg: Config, year: int, month: int) -> tuple[pd.Timestamp, pd.Timestamp, bool]:
    """(start, end, end_inclusive) of a calendar month, clipped to the radar's period window.

    Months are [first of month, first of next month). The stamp at the first of the next month
    belongs to the next month's file, EXCEPT after the last month: radar stamps are instants, so
    the period window is [start, end] (period.window_end_inclusive) and the closing stamp
    2015-09-01 00:00 goes into the last month's file.
    """
    p_start, p_end, p_incl = period_window(cfg, cfg["radar"]["stamp_kind"])
    start = pd.Timestamp(year=year, month=month, day=1)
    end = start + pd.offsets.MonthBegin(1)
    if end >= p_end:
        return max(start, p_start), p_end, p_incl
    return max(start, p_start), end, False


def load_radar_month(cfg: Config, year: int, month: int) -> xr.Dataset:
    """One calendar month of radar data in the radar L1 schema.

    The raw pseudo-dBZ (int32, 255 = missing) is read WITHOUT automatic decoding and decoded
    here once, dBZ = value * scale_factor + add_offset, using the file's own attributes, so the
    scaling can never be applied twice. rain_rate = ((10^(dBZ/10)) / a)^(1/b) with a, b from the
    file's zr_a / zr_b attributes (config fallback with a warning). No 0 dBZ threshold and no
    clutter removal. Time is kept as in the file ("instantaneous", per the readme).
    """
    c = cfg["radar"]
    interval = pd.Timedelta(c["interval"])
    start, end, incl = month_window(cfg, year, month)
    if start >= end:
        raise ValueError(f"{year}-{month:02d} is outside the configured period")

    with xr.open_dataset(cfg.raw_path("raw_radar"), mask_and_scale=False, decode_times=True) as raw:
        var = raw["data"]
        attrs = dict(var.attrs)
        times = pd.DatetimeIndex(raw["time"].values).as_unit("ns") + pd.Timedelta(seconds=int(c["time_shift_s"]))
        sel = np.flatnonzero((times >= start) & ((times <= end) if incl else (times < end)))
        block = slice(int(sel.min()), int(sel.max()) + 1) if sel.size else slice(0, 0)
        values = var.isel(time=block).transpose("time", "y", "x").values
        t_block = times[block]
        x = raw["x"].values.astype("float64")
        y = raw["y"].values.astype("float64")
        lon = raw["lon"].transpose("y", "x").values
        lat = raw["lat"].transpose("y", "x").values
        gattrs = dict(raw.attrs)
        crs_attrs = dict(raw["crs"].attrs) if "crs" in raw else {}

    fill = attrs.get("_FillValue", attrs.get("missing_value"))
    scale = attrs.get("scale_factor", 1.0)
    offset = attrs.get("add_offset", 0.0)
    is_missing = values == fill
    dbz = np.where(is_missing, np.nan, values.astype("float64") * float(scale) + float(offset)).astype("float32")

    try:
        a, b = float(attrs["zr_a"]), float(attrs["zr_b"])  # stored as strings in OpenMRG
        zr_source = "radar.nc variable 'data' attributes zr_a, zr_b"
    except (KeyError, ValueError):
        a, b = float(c["zr_fallback"]["a"]), float(c["zr_fallback"]["b"])
        zr_source = "config radar.zr_fallback (file had no usable zr_a / zr_b)"
        log.warning("radar.nc has no usable zr_a/zr_b; using config fallback a=%s, b=%s", a, b)
    rain = schema.zr_rain_rate(dbz, a, b).astype("float32")

    grid, arr, counts = _regularize(t_block, {"dbz": dbz, "rain_rate": rain}, start, end, interval,
                                    end_inclusive=incl)
    x3006, y3006 = lonlat_to_metric(cfg, lon, lat)

    ds = xr.Dataset(
        {"dbz": (schema.RADAR_DIMS, arr["dbz"]), "rain_rate": (schema.RADAR_DIMS, arr["rain_rate"])},
        coords={"time": grid.to_numpy(), "y": y, "x": x, "lon": (("y", "x"), lon), "lat": (("y", "x"), lat),
                "x_3006": (("y", "x"), x3006), "y_3006": (("y", "x"), y3006)},
    )
    _set_attrs(ds, schema.RADAR_DATA_VARS, schema.RADAR_COORDS)
    ds["x"].attrs["standard_name"] = "projection_x_coordinate"
    ds["y"].attrs["standard_name"] = "projection_y_coordinate"
    n_missing = int(is_missing.sum())
    ds.attrs = {
        **_base_attrs(cfg, f"L1 radar, 5 min, {year}-{month:02d}", ["raw_radar"]),
        **_time_attrs("instantaneous", c["interval"], c["time_meaning_original"], c["time_shift_s"],
                      window_text(start, end, incl)),
        "month": f"{year}-{month:02d}",
        "period_end_stamp": (f"This is the last month file: it also holds the stamp {end:%Y-%m-%d %H:%M}, "
                             "the instant that closes the period (radar stamps are instants, so the period "
                             "window is closed at its end)." if incl else "none (not the last month file)"),
        "raw_file_doi": str(gattrs.get("doi", "")),
        "raw_file_version": str(gattrs.get("version", "")),
        "crs_proj4": str(gattrs.get("proj_string", "")),
        "crs_wkt": str(crs_attrs.get("spatial_ref", "")),
        "crs_grid_mapping_name": str(crs_attrs.get("grid_mapping_name", "")),
        "lonlat_assumed_crs": cfg["reference"]["lonlat_crs"],
        "metric_crs": cfg["reference"]["metric_crs"],
        "zr_a": a, "zr_b": b, "zr_source": zr_source,
        "decoding": f"dBZ = raw * {float(scale)} + {float(offset)}; raw == {fill} -> NaN (decoded once, here)",
        "transforms": "fill -> NaN; pseudo-dBZ -> dBZ; rain_rate from Z-R; reindexed to a regular grid between the first "
                      "and last raw stamp. No threshold, no clutter removal.",
        "n_missing_cells_raw": n_missing,
        "n_cells_raw": int(is_missing.size),
        "n_steps_all_missing_raw": int(is_missing.all(axis=(1, 2)).sum()) if is_missing.size else 0,
        **counts,
    }
    return ds


def radar_months(cfg: Config) -> list[tuple[int, int]]:
    """(year, month) of every calendar month that overlaps the period."""
    months = pd.period_range(cfg.period_start, cfg.period_end - pd.Timedelta("1ns"), freq="M")
    return [(p.year, p.month) for p in months]


def radar_steps_outside_months(cfg: Config) -> int:
    """Raw radar steps outside the radar's period window, i.e. in no month file."""
    start, end, incl = period_window(cfg, cfg["radar"]["stamp_kind"])
    with xr.open_dataset(cfg.raw_path("raw_radar"), mask_and_scale=False) as raw:
        t = pd.DatetimeIndex(raw["time"].values) + pd.Timedelta(seconds=int(cfg["radar"]["time_shift_s"]))
    after = (t > end) if incl else (t >= end)
    return int(((t < start) | after).sum())


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------
def _encoding(ds: xr.Dataset) -> dict:
    """zlib level 4; float variables float32 with NaN fill; integer QC uint8 without fill."""
    enc = {"time": {"units": "seconds since 1970-01-01 00:00:00", "calendar": "proleptic_gregorian", "dtype": "int64"}}
    for name, v in ds.data_vars.items():
        e = {"zlib": True, "complevel": 4}
        if v.dtype.kind == "f":
            e.update({"dtype": "float32", "_FillValue": np.float32(np.nan)})
        else:
            e.update({"dtype": str(v.dtype), "_FillValue": None})
        if v.dims == schema.RADAR_DIMS:
            e["chunksizes"] = (min(288, v.shape[0]), v.shape[1], v.shape[2])  # one day per chunk
        enc[name] = e
    return enc


def check_roundtrip_exact(ds_memory: xr.Dataset, ds_file: xr.Dataset) -> None:
    """Every data variable and coordinate must be identical after writing (NaN in the same places)."""
    for name in list(ds_memory.data_vars) + list(ds_memory.coords):
        a, b = ds_memory[name].values, ds_file[name].values
        if a.shape != b.shape:
            raise ValueError(f"round trip: {name} shape {a.shape} became {b.shape}")
        if a.dtype.kind in "UOS":
            same = np.array_equal(a.astype(str), b.astype(str))
        elif a.dtype.kind == "f":
            same = np.array_equal(a, b.astype(a.dtype), equal_nan=True)
        else:
            same = np.array_equal(a, b)
        if not same:
            raise ValueError(f"round trip: {name} changed when written")


def write_l1_reference(ds: xr.Dataset, path: Path, cfg: Config, stage: str) -> dict:
    """Write an L1 reference dataset to NetCDF4 (zlib 4, float32), re-open it and check it is identical.

    Written as <name>.tmp and renamed only after the check, so a file under its final name is
    always complete and verified. Returns {"path", "file_size_bytes"}.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.unlink(missing_ok=True)
    ds = ds.copy()
    ds.attrs.update(provenance_attrs(cfg, stage))
    ds.to_netcdf(tmp, format="NETCDF4", engine="netcdf4", encoding=_encoding(ds))
    try:
        with xr.open_dataset(tmp) as written:
            written.load()
            check_roundtrip_exact(ds, written)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    os.replace(tmp, path)
    return {"path": path, "file_size_bytes": path.stat().st_size}


# ---------------------------------------------------------------------------
# Gauge table
# ---------------------------------------------------------------------------
GAUGE_TABLE_COLUMNS = ["id", "network", "gauge_type", "lon", "lat", "x", "y", "quantization_mm", "height_m",
                       "interval_s", "n_valid", "completeness_pct", "total_mm", "n_negative"]


def gauge_table(datasets: list[xr.Dataset]) -> pd.DataFrame:
    """One row per gauge of every given gauge L1 dataset (columns GAUGE_TABLE_COLUMNS)."""
    rows = []
    for ds in datasets:
        amount = ds.rainfall_amount.values.astype("float64")
        n_time = amount.shape[0]
        for j, gid in enumerate(ds.id.values.astype(str)):
            col = amount[:, j]
            valid = np.isfinite(col)
            rows.append({
                "id": gid, "network": str(ds.network.values[j]), "gauge_type": str(ds.gauge_type.values[j]),
                "lon": float(ds.lon.values[j]), "lat": float(ds.lat.values[j]),
                "x": round(float(ds.x.values[j]), 2), "y": round(float(ds.y.values[j]), 2),
                "quantization_mm": float(ds.quantization_mm.values[j]), "height_m": float(ds.height_m.values[j]),
                "interval_s": int(ds.attrs["interval_s"]), "n_valid": int(valid.sum()),
                "completeness_pct": round(100.0 * valid.sum() / n_time, 4) if n_time else float("nan"),
                "total_mm": round(float(np.nansum(col)), 3), "n_negative": int(np.sum(col < 0)),
            })
    return pd.DataFrame(rows, columns=GAUGE_TABLE_COLUMNS)


def write_gauge_table(table: pd.DataFrame, cfg: Config) -> Path:
    """tables/gauges.csv plus a provenance sidecar, like sublink_id_map.csv."""
    cfg.tables_dir.mkdir(parents=True, exist_ok=True)
    path = cfg.tables_dir / "gauges.csv"
    table.to_csv(path, index=False)
    write_sidecar(path, cfg, "ingest_gauges", {
        "description": "one row per rain gauge (City + SMHI) with position (lon/lat WGS84, x/y EPSG:3006) and "
                       "L1 statistics over the period",
        "columns": GAUGE_TABLE_COLUMNS,
    })
    return path
