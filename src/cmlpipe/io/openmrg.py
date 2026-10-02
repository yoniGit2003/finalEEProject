"""OpenMRG adapter: raw cml.nc + cml_metadata.csv -> L1 CML files, one per UTC day.

Raw layout (checked in docs/data_notes/openmrg_cml_inspection.md):
  cml.nc            rsl, tsl  float32 (time, sublink), missing value 1e10,
                    time = int seconds since 1970 (UTC, start of sampling window)
  cml_metadata.csv  one row per sub-link: Link, Direction (A/B), Sublink,
                    Near/Far lat/lon, Frequency_GHz, Polarization, Length_km,
                    Near/FarAntennaDiameter_m

Mapping to L1 (see cmlpipe.schema):
  Link -> cml_id; Direction A -> sublink_1, B -> sublink_2 (from the config).
  In the raw metadata Near = transmitter and Far = receiver OF EACH SUB-LINK,
  so A.Near is the same place as B.Far. We fix the sites per link:
      site_0 = transmitter of sublink_1 (A.Near = B.Far)
      site_1 = receiver    of sublink_1 (A.Far  = B.Near)
  so sublink_1 goes site_0 -> site_1 and sublink_2 goes site_1 -> site_0
  (stored in the tx_site / rx_site coordinates).

Values are never changed: the only transforms are sentinel -> NaN, dropping
duplicate timestamps, reindexing onto the regular 10 s grid (new steps = NaN),
reshaping, and tl = tsl - rsl.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from cmlpipe.io import schema
from cmlpipe.config import Config
from cmlpipe.provenance import provenance_attrs, write_sidecar

log = logging.getLogger(__name__)

# Which site transmits on each sub-link. This follows from the site definition
# above (site_0 = transmitter of sublink_1), so it is fixed, not a config choice.
TX_SITE = {"sublink_1": 0, "sublink_2": 1}

# Columns of cml_metadata.csv that we read.
RAW_COLUMNS = [
    "Link", "Direction", "Sublink",
    "NearLatitude_DecDeg", "NearLongitude_DecDeg", "FarLatitude_DecDeg", "FarLongitude_DecDeg",
    "Frequency_GHz", "Polarization", "Length_km", "NearAntennaDiameter_m", "FarAntennaDiameter_m",
]
# Link-level columns of the tidy table (the same for both sub-links of a link).
LINK_COLUMNS = [
    "length", "site_0_lat", "site_0_lon", "site_1_lat", "site_1_lon",
    "site_0_antenna_diameter", "site_1_antenna_diameter",
]


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------
def _normalize_polarization(raw: pd.Series, pol_map: dict) -> pd.Series:
    """Map raw spellings ('Vertical', 'V', 'vertical ') to 'v' / 'h'. Unknown spellings raise."""
    key = raw.astype(str).str.strip().str.lower()
    unknown = sorted(set(key) - set(pol_map))
    if unknown:
        raise ValueError(f"Unknown polarization spellings {unknown}; add them to cml.polarization_map in the config")
    return key.map(pol_map)


def load_cml_metadata(cfg: Config) -> pd.DataFrame:
    """Read cml_metadata.csv and return a tidy table with one row per (cml_id, sublink_id).

    Columns: cml_id, sublink_id, raw_sublink_id, direction, tx_site, rx_site,
    frequency [GHz], polarization ('v'/'h'), length [km], site_0_lat, site_0_lon,
    site_1_lat, site_1_lon [deg], site_0_antenna_diameter, site_1_antenna_diameter [m].

    A link that lacks one of its two directions is kept: the missing sub-link gets
    raw_sublink_id = -1, empty direction/polarization and NaN frequency (its data
    will be all NaN). Such links are listed in the log.
    """
    raw = pd.read_csv(cfg.raw_cml_metadata)
    missing_cols = [c for c in RAW_COLUMNS if c not in raw.columns]
    if missing_cols:
        raise ValueError(f"{cfg.raw_cml_metadata} lacks columns {missing_cols}")

    dir_map = cfg["cml"]["direction_to_sublink"]
    direction = raw["Direction"].astype(str).str.strip()
    unknown = sorted(set(direction) - set(dir_map))
    if unknown:
        raise ValueError(f"Unknown Direction codes {unknown}; expected {sorted(dir_map)}")
    if raw.duplicated(["Link", "Direction"]).any():
        raise ValueError("cml_metadata.csv has more than one row for the same (Link, Direction)")
    if raw["Sublink"].duplicated().any():
        raise ValueError("cml_metadata.csv has duplicated Sublink numbers")

    sublink_id = direction.map(dir_map)
    tx_site = sublink_id.map(TX_SITE)
    # Near is the transmitter of THIS sub-link. If this sub-link transmits from
    # site_0 then Near = site_0, otherwise Near = site_1.
    near_is_site0 = (tx_site == 0).to_numpy()

    def pick(near_col: str, far_col: str, site: int) -> np.ndarray:
        near, far = raw[near_col].to_numpy(float), raw[far_col].to_numpy(float)
        return np.where(near_is_site0, near, far) if site == 0 else np.where(near_is_site0, far, near)

    tidy = pd.DataFrame({
        "cml_id": raw["Link"].astype("int64"),
        "sublink_id": sublink_id.astype(object),
        "raw_sublink_id": raw["Sublink"].astype("int64"),
        "direction": direction.astype(object),
        "tx_site": tx_site.astype("int8"),
        "rx_site": (1 - tx_site).astype("int8"),
        "frequency": raw["Frequency_GHz"].astype(float),
        "polarization": _normalize_polarization(raw["Polarization"], cfg["cml"]["polarization_map"]).astype(object),
        "length": raw["Length_km"].astype(float),
        "site_0_lat": pick("NearLatitude_DecDeg", "FarLatitude_DecDeg", 0),
        "site_0_lon": pick("NearLongitude_DecDeg", "FarLongitude_DecDeg", 0),
        "site_1_lat": pick("NearLatitude_DecDeg", "FarLatitude_DecDeg", 1),
        "site_1_lon": pick("NearLongitude_DecDeg", "FarLongitude_DecDeg", 1),
        "site_0_antenna_diameter": pick("NearAntennaDiameter_m", "FarAntennaDiameter_m", 0),
        "site_1_antenna_diameter": pick("NearAntennaDiameter_m", "FarAntennaDiameter_m", 1),
    })

    # Both sub-links of a link describe the same two sites. Check that they agree
    # after the mapping; if not, keep sublink_1's values and warn.
    link_values = tidy.groupby("cml_id")[LINK_COLUMNS]
    n_disagree = int((link_values.nunique(dropna=False) > 1).any(axis=1).sum())
    if n_disagree:
        log.warning("%d links have A/B metadata that disagree after the site mapping; using sublink_1", n_disagree)
    per_link = tidy.sort_values("sublink_id").groupby("cml_id")[LINK_COLUMNS].first()

    # Add a placeholder row for every missing direction, so every link has both sub-links.
    all_pairs = pd.MultiIndex.from_product([sorted(tidy["cml_id"].unique()), schema.SUBLINK_IDS],
                                           names=["cml_id", "sublink_id"])
    present = pd.MultiIndex.from_frame(tidy[["cml_id", "sublink_id"]])
    absent = all_pairs.difference(present)
    if len(absent):
        log.warning("%d links have only one direction; the other sub-link is filled with NaN: %s",
                    len(absent), sorted(set(absent.get_level_values("cml_id"))))
        filler = pd.DataFrame(absent.to_list(), columns=["cml_id", "sublink_id"])
        filler["raw_sublink_id"] = schema.MISSING_SUBLINK
        filler["direction"] = ""
        filler["tx_site"] = filler["sublink_id"].map(TX_SITE).astype("int8")
        filler["rx_site"] = (1 - filler["tx_site"]).astype("int8")
        filler["frequency"] = np.nan
        filler["polarization"] = ""
        tidy = pd.concat([tidy, filler], ignore_index=True)

    # Link-level columns always come from per_link, so the two rows of a link agree.
    tidy = tidy.drop(columns=LINK_COLUMNS).join(per_link, on="cml_id")
    tidy = tidy.sort_values(["cml_id", "sublink_id"]).reset_index(drop=True)
    log.info("Metadata: %d links, %d sub-links present, %d absent",
             tidy["cml_id"].nunique(), int((tidy.raw_sublink_id != schema.MISSING_SUBLINK).sum()), len(absent))
    return tidy


def write_sublink_id_map(meta: pd.DataFrame, cfg: Config) -> Path:
    """Write tables/sublink_id_map.csv (cml_id, sublink_id, raw_sublink_id, direction) plus a provenance sidecar."""
    cfg.tables_dir.mkdir(parents=True, exist_ok=True)
    path = cfg.tables_dir / "sublink_id_map.csv"
    meta[["cml_id", "sublink_id", "raw_sublink_id", "direction"]].to_csv(path, index=False)
    write_sidecar(path, cfg, "ingest_cml", {"description": "L1 (cml_id, sublink_id) <-> raw OpenMRG Sublink / Direction",
                              "missing_sublink_value": schema.MISSING_SUBLINK})
    return path


# ---------------------------------------------------------------------------
# Data: one day
# ---------------------------------------------------------------------------
def _read_raw_day(cfg: Config, day: pd.Timestamp) -> tuple[np.ndarray, dict, np.ndarray, dict]:
    """Read the raw rows whose timestamp falls in [day, day + 1 day).

    Only that block of rows is read from disk (the raw file is stored row by row,
    so a block of consecutive time steps is one contiguous read).
    Returns (times, {"rsl": array, "tsl": array}, raw sublink ids, raw global attrs).
    Arrays are (time, sublink), float32, with the sentinel still in place.
    """
    # mask_and_scale=False: keep 1e10 as a number so that WE convert it (and count it).
    with xr.open_dataset(cfg.raw_cml, mask_and_scale=False) as raw:
        times = raw["time"].values
        in_day = (times >= day.to_datetime64()) & (times < (day + pd.Timedelta("1D")).to_datetime64())
        rows = np.flatnonzero(in_day)
        sublinks = raw["sublink"].values.astype("int64")
        attrs = dict(raw.attrs)
        if rows.size == 0:
            empty = np.empty((0, sublinks.size), dtype="float32")
            return times[:0], {"rsl": empty, "tsl": empty.copy()}, sublinks, attrs
        block = slice(int(rows.min()), int(rows.max()) + 1)
        keep = in_day[block]  # the block may contain rows of other days if times are unsorted
        values = {}
        for name in ("rsl", "tsl"):
            # transpose makes the order explicit even if a future file stores (sublink, time)
            values[name] = raw[name].isel(time=block).transpose("time", "sublink").values[keep]
        return times[block][keep], values, sublinks, attrs


def load_cml_day(cfg: Config, date, meta: pd.DataFrame | None = None) -> xr.Dataset:
    """Return one UTC day [00:00, 24:00) of OpenMRG CML data in the L1 schema.

    Steps (each one is counted and stored in the dataset attributes):
      1. read only the rows of this day from cml.nc
      2. sentinel (1e10) -> NaN
      3. drop duplicated timestamps, keeping the first occurrence
      4. drop timestamps that are not on the 10 s grid (none in OpenMRG; warned)
      5. reindex onto the full regular grid; inserted steps are NaN
      6. reshape (time, raw sublink) -> (time, sublink_id, cml_id) and attach metadata
      7. tl = tsl - rsl

    Parameters
    ----------
    cfg : the loaded config.
    date : the UTC day, e.g. "2015-07-28".
    meta : output of load_cml_metadata (read from disk if not given).
    """
    day = pd.Timestamp(date).normalize()
    interval = cfg.sampling_interval
    if meta is None:
        meta = load_cml_metadata(cfg)
    sentinel = np.float32(cfg["cml"]["raw_missing_value"])

    times, values, raw_sublinks, raw_attrs = _read_raw_day(cfg, day)

    # 2. sentinel -> NaN
    n_sentinel = {}
    for name in ("rsl", "tsl"):
        is_sentinel = values[name] == sentinel
        n_sentinel[name] = int(is_sentinel.sum())
        values[name] = np.where(is_sentinel, np.float32(np.nan), values[name]).astype("float32")

    # 3. duplicated timestamps: keep the first. Also count how many duplicates
    #    carried DIFFERENT values, since that would mean real information is lost.
    t_index = pd.DatetimeIndex(times)
    is_dup = t_index.duplicated(keep="first")
    n_dup = int(is_dup.sum())
    n_dup_differ = 0
    if n_dup:
        first_pos = pd.Series(np.arange(len(t_index))).groupby(t_index).transform("first").to_numpy()
        for name in ("rsl", "tsl"):
            a, b = values[name][is_dup], values[name][first_pos[is_dup]]
            same = (a == b) | (np.isnan(a) & np.isnan(b))
            n_dup_differ += int((~same.all(axis=1)).sum())
        log.warning("%s: %d duplicated timestamps dropped (%d with different values)", day.date(), n_dup, n_dup_differ)
    times = times[~is_dup]
    for name in values:
        values[name] = values[name][~is_dup]

    # 4. timestamps off the regular grid cannot be placed on it without changing
    #    the time; drop them and warn loudly (count is stored).
    # Offset from midnight in nanoseconds. as_unit("ns") because pandas 3 may use
    # other time units, and interval.value below is in nanoseconds.
    offset = (pd.DatetimeIndex(times).as_unit("ns") - day).as_unit("ns").asi8
    on_grid = offset % interval.value == 0
    n_offgrid = int((~on_grid).sum())
    if n_offgrid:
        log.warning("%s: %d timestamps not on the %s grid were dropped", day.date(), n_offgrid, interval)

    # 5. reindex onto the full regular grid of the day
    # as_unit("ns"): the same time resolution xarray uses when it decodes the raw file
    grid = pd.date_range(day, day + pd.Timedelta("1D"), freq=interval, inclusive="left").as_unit("ns")
    pos = offset[on_grid] // interval.value
    gridded = {}
    for name in ("rsl", "tsl"):
        full = np.full((len(grid), raw_sublinks.size), np.nan, dtype="float32")
        full[pos] = values[name][on_grid]
        gridded[name] = full
    n_inserted = len(grid) - int(on_grid.sum())

    # 6. reshape to (time, sublink_id, cml_id)
    cml_ids = np.sort(meta["cml_id"].unique())
    cml_pos = {c: i for i, c in enumerate(cml_ids)}
    sub_pos = {s: i for i, s in enumerate(schema.SUBLINK_IDS)}
    raw_col = {int(s): i for i, s in enumerate(raw_sublinks)}
    meta_raw_ids = set(meta["raw_sublink_id"]) - {schema.MISSING_SUBLINK}
    not_in_file = sorted(meta_raw_ids - set(raw_col))
    not_in_meta = sorted(set(raw_col) - meta_raw_ids)
    if not_in_file:
        log.warning("%d sub-links in the metadata are not in cml.nc (their data are NaN): %s", len(not_in_file), not_in_file)
    if not_in_meta:
        log.warning("%d sub-links in cml.nc have no metadata and are dropped: %s", len(not_in_meta), not_in_meta)

    shape3 = (len(grid), len(schema.SUBLINK_IDS), len(cml_ids))
    data3 = {name: np.full(shape3, np.nan, dtype="float32") for name in ("rsl", "tsl")}
    have = meta[meta["raw_sublink_id"].isin(list(raw_col))]
    s_idx = have["sublink_id"].map(sub_pos).to_numpy()
    c_idx = have["cml_id"].map(cml_pos).to_numpy()
    cols = have["raw_sublink_id"].map(raw_col).to_numpy()
    for name in ("rsl", "tsl"):
        data3[name][:, s_idx, c_idx] = gridded[name][:, cols]

    # Metadata tables in (sublink_id, cml_id) and (cml_id,) shape
    wide = meta.set_index(["sublink_id", "cml_id"]).reindex(
        pd.MultiIndex.from_product([schema.SUBLINK_IDS, cml_ids]))
    per_link = meta.groupby("cml_id")[LINK_COLUMNS].first().reindex(cml_ids)

    def sub2d(col, dtype):
        return wide[col].to_numpy().astype(dtype).reshape(len(schema.SUBLINK_IDS), len(cml_ids))

    dims3, dims2 = schema.DIMS, ("sublink_id", "cml_id")
    ds = xr.Dataset(
        data_vars={
            "tsl": (dims3, data3["tsl"]),
            "rsl": (dims3, data3["rsl"]),
            # 7. total loss. Both inputs are float32 on a 0.1 dB grid.
            "tl": (dims3, data3["tsl"] - data3["rsl"]),
        },
        coords={
            "time": grid.to_numpy(),
            "sublink_id": np.array(schema.SUBLINK_IDS),
            "cml_id": cml_ids.astype("int64"),
            "raw_sublink_id": (dims2, sub2d("raw_sublink_id", "int64")),
            "direction": (dims2, sub2d("direction", str)),
            "tx_site": (dims2, sub2d("tx_site", "int8")),
            "rx_site": (dims2, sub2d("rx_site", "int8")),
            "frequency": (dims2, sub2d("frequency", "float64")),
            "polarization": (dims2, sub2d("polarization", str)),
            **{col: (("cml_id",), per_link[col].to_numpy(float)) for col in LINK_COLUMNS},
        },
    )

    # Units and descriptions from the schema
    for name, (units, long_name) in schema.DATA_VARS.items():
        ds[name].attrs = {"units": units, "long_name": long_name}
    for name, (_, units, long_name) in schema.COORDS.items():
        ds[name].attrs = {"long_name": long_name} if units is None else {"units": units, "long_name": long_name}
    ds["time"].attrs = {"long_name": "time (UTC)", "standard_name": "time", "time_meaning": cfg["cml"]["time_meaning"]}
    ds["tx_site"].attrs["comment"] = "sublink_1 transmits site_0 -> site_1, sublink_2 transmits site_1 -> site_0"

    ds.attrs = {
        "title": f"{cfg['dataset']['title']} - L1 CML, {day:%Y-%m-%d}",
        "dataset": cfg.dataset_name,
        "processing_level": "L1_standard",
        "source": cfg["dataset"]["source"],
        "source_doi": cfg["dataset"]["source_doi"],
        "licence": cfg["dataset"]["licence"],
        "raw_file": Path(cfg["paths"]["raw_cml"]).as_posix(),
        "raw_file_doi": str(raw_attrs.get("doi", "")),
        "raw_file_version": str(raw_attrs.get("version", "")),
        "time_meaning": cfg["cml"]["time_meaning"],
        "time_zone": "UTC",
        "sampling_interval": cfg["cml"]["sampling_interval"],
        "day_utc": f"{day:%Y-%m-%d}",
        "site_convention": ("site_0 = transmitter of sublink_1 (raw Near of Direction A = Far of B); "
                            "site_1 = receiver of sublink_1. Raw Near/Far = transmitter/receiver of each sub-link."),
        "sublink_mapping": ", ".join(f"{k} -> {v}" for k, v in cfg["cml"]["direction_to_sublink"].items()),
        "transforms": ("sentinel -> NaN; duplicate timestamps dropped (first kept); reindexed to regular grid "
                       "(inserted steps NaN); reshaped; tl = tsl - rsl. No other change to values."),
        "n_raw_steps_in_day": int(len(t_index)),
        "n_duplicates_removed": n_dup,
        "n_duplicates_with_different_values": n_dup_differ,
        "n_offgrid_steps_dropped": n_offgrid,
        "n_steps_inserted": int(n_inserted),
        "n_sentinel_rsl": n_sentinel["rsl"],
        "n_sentinel_tsl": n_sentinel["tsl"],
    }
    return ds


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------
def _encoding(cfg: Config, ds: xr.Dataset) -> dict:
    """NetCDF encoding: compressed int16 x 0.1 dB (lossless for OpenMRG) or float32."""
    enc_cfg = cfg["encoding"]
    enc = {"time": {"units": "seconds since 1970-01-01 00:00:00", "calendar": "proleptic_gregorian", "dtype": "int64"}}
    for name in schema.DATA_VARS:
        e = {"zlib": bool(enc_cfg["zlib"]), "complevel": int(enc_cfg["complevel"])}
        if enc_cfg["signal_dtype"] == "int16":
            # Make sure every value fits into int16 after scaling, else int16 would wrap around.
            limit = 32767 * enc_cfg["scale_factor"]
            v = ds[name].values
            if np.nanmax(np.abs(v), initial=0.0) >= limit:
                raise ValueError(f"{name} has |values| >= {limit}: cannot be stored as int16 x {enc_cfg['scale_factor']}")
            e.update({"dtype": "int16", "scale_factor": float(enc_cfg["scale_factor"]),
                      "_FillValue": np.int16(enc_cfg["fill_value"])})
        elif enc_cfg["signal_dtype"] == "float32":
            e.update({"dtype": "float32", "_FillValue": np.float32(np.nan)})
        else:
            raise ValueError(f"encoding.signal_dtype must be int16 or float32, not {enc_cfg['signal_dtype']!r}")
        enc[name] = e
    return enc


def check_roundtrip(ds_memory: xr.Dataset, ds_file: xr.Dataset, atol: float) -> dict:
    """Compare what we meant to write with what was written.

    For each data variable: NaN must be in exactly the same places, and finite
    values may differ by at most `atol` (float32 representation error only).
    Coordinates must be identical. Raises ValueError on any difference.
    Returns {var: {"max_abs_diff": ..., "bit_identical": ...}}.
    """
    result = {}
    for name in schema.DATA_VARS:
        a = ds_memory[name].values
        b = ds_file[name].values
        if a.shape != b.shape:
            raise ValueError(f"round trip: {name} shape {a.shape} became {b.shape}")
        nan_a, nan_b = np.isnan(a), np.isnan(b)
        if not np.array_equal(nan_a, nan_b):
            raise ValueError(f"round trip: NaN pattern of {name} changed ({int((nan_a != nan_b).sum())} cells)")
        diff = np.abs(a[~nan_a].astype("float64") - b[~nan_b].astype("float64"))
        max_diff = float(diff.max()) if diff.size else 0.0
        if max_diff > atol:
            raise ValueError(f"round trip: {name} changed by up to {max_diff:.3g} dB (> {atol}); encoding is lossy")
        identical = bool(np.array_equal(a[~nan_a], b[~nan_b].astype(a.dtype)))
        result[name] = {"max_abs_diff": max_diff, "bit_identical": identical}
    for name in ds_memory.coords:
        a, b = ds_memory[name].values, ds_file[name].values
        same = np.array_equal(a.astype(str), b.astype(str)) if a.dtype.kind in "UO" else np.array_equal(a, b, equal_nan=True)
        if not same:
            raise ValueError(f"round trip: coordinate {name} changed")
    return result


def write_l1_day(ds: xr.Dataset, path: Path, cfg: Config) -> dict:
    """Write one L1 day to NetCDF4 (zlib), then re-open it and verify it is lossless.

    The file is first written as <name>.tmp and only renamed to its final name after
    the round-trip check passed. So a file with the final name is always complete
    and verified; a crash leaves at most a .tmp file behind.

    Returns {"path", "file_size_bytes", "roundtrip"}.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    if tmp.exists():
        tmp.unlink()

    ds = ds.copy()
    ds.attrs.update(provenance_attrs(cfg, "ingest_cml"))
    ds.attrs["signal_encoding"] = (f"int16, scale_factor {cfg['encoding']['scale_factor']}"
                                   if cfg["encoding"]["signal_dtype"] == "int16" else "float32")
    ds.to_netcdf(tmp, format="NETCDF4", engine="netcdf4", encoding=_encoding(cfg, ds))

    try:
        with xr.open_dataset(tmp) as written:
            written.load()
            roundtrip = check_roundtrip(ds, written, float(cfg["encoding"]["roundtrip_atol_db"]))
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    os.replace(tmp, path)
    return {"path": path, "file_size_bytes": path.stat().st_size, "roundtrip": roundtrip}
