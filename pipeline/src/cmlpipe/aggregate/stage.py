"""Aggregation stages: aggregate_cml (per day, parallel), aggregate_gauges, aggregate_radar.

Outputs (under data/work/openmrg/L3_products/):
  cml/<res>/cml_<res>_YYYY-MM-DD.nc          (time, sublink_id, cml_id)  one file per day and resolution
  gauges/gauges_city_<res>.nc                (time, id)                  1min 5min 15min 1h
  gauges/gauges_smhi_<res>.nc                (time, id)                  15min 1h
  radar/radar_along_links_<res>.nc           (time, cml_id)              5min 15min 1h
  radar/radar_at_gauges_<res>.nc             (time, id)                  5min 15min 1h
Every file: provenance attributes (stage, stage_hash, ...), the time attributes of aggregate.time_attrs and
the hashes of its inputs. A file is written to <name>.tmp, re-opened and compared with the in-memory data,
and only then renamed; a failed file never appears under its final name. L1, L2 and the geometry tables
are only read.
"""

from __future__ import annotations

import json
import os
import time
import warnings
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr

from cmlpipe import aggregate
from cmlpipe import geometry
from cmlpipe.geometry import stage as geometry_runner
from cmlpipe.qc import stage as qc_runner
from cmlpipe.io import schema
from cmlpipe.config import Config
from cmlpipe.provenance import display_path
from cmlpipe.io.netcdf import file_attrs
from cmlpipe.io.netcdf import existing_stage_hash
from cmlpipe.parallel import config_from_parts
from cmlpipe.provenance import provenance_attrs

TIME_ENCODING = {"units": "seconds since 1970-01-01 00:00:00", "calendar": "proleptic_gregorian", "dtype": "int64"}


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
def l3_dir(cfg: Config, key: str) -> Path:
    """An L3 product folder (key of the config `paths` section)."""
    return cfg.output_root / cfg["paths"][key]


def cml_path(cfg: Config, res: str, day) -> Path:
    """L3 CML file of one resolution and day."""
    return l3_dir(cfg, "l3_cml_dir") / res / f"cml_{res}_{pd.Timestamp(day):%Y-%m-%d}.nc"


def gauge_path(cfg: Config, network: str, res: str) -> Path:
    """L3 gauge file of a network and resolution."""
    return l3_dir(cfg, "l3_gauges_dir") / f"gauges_{network}_{res}.nc"


def radar_path(cfg: Config, what: str, res: str) -> Path:
    """what = 'along_links' or 'at_gauges'."""
    return l3_dir(cfg, "l3_radar_dir") / f"radar_{what}_{res}.nc"


# ---------------------------------------------------------------------------
# Writing with verification
# ---------------------------------------------------------------------------
def _encoding(ds: xr.Dataset, enc_cfg: dict, int16_vars: set[str]) -> dict:
    enc = {"time": dict(TIME_ENCODING)}
    z = {"zlib": bool(enc_cfg.get("zlib", True)), "complevel": int(enc_cfg.get("complevel", 4))}
    for name, v in ds.data_vars.items():
        if name in int16_vars:
            enc[name] = {"dtype": "int16", "scale_factor": float(enc_cfg["int16_scale_factor"]), "add_offset": 0.0,
                         "_FillValue": np.int16(-32768), **z}
        elif v.dtype.kind == "f":
            enc[name] = {"dtype": "float32", "_FillValue": np.float32(np.nan), **z}
        else:
            enc[name] = {"dtype": v.dtype.name, "_FillValue": None, **z}
    return enc


def write_verified(ds: xr.Dataset, path: Path, enc_cfg: dict, int16_vars: set[str] = frozenset()) -> dict:
    """Write ds to path via <path>.tmp, re-open, compare, rename. Returns size and the max int16 error.

    float32 variables must come back bit-identical (they are stored as float32); int16 variables within
    int16_max_error_db (0.005 dB = half the 0.01 scale step); integer variables and time exactly.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.unlink(missing_ok=True)
    max_err = 0.0
    try:
        ds.to_netcdf(tmp, format="NETCDF4", engine="netcdf4", encoding=_encoding(ds, enc_cfg, int16_vars))
        with xr.open_dataset(tmp) as back:
            back.load()
            if not np.array_equal(back.time.values, ds.time.values):
                raise ValueError("time axis changed in the round trip")
            for name, v in ds.data_vars.items():
                a, b = v.values, back[name].values
                if name in int16_vars:
                    if not np.array_equal(np.isnan(a), np.isnan(b)):
                        raise ValueError(f"{name}: NaN positions changed in the int16 round trip")
                    err = float(np.nanmax(np.abs(a.astype("float64") - b))) if np.isfinite(a).any() else 0.0
                    max_err = max(max_err, err)
                    # half a scale step, plus float32 representation error of the input (~4e-6 at 100 dB)
                    if err > float(enc_cfg["int16_max_error_db"]) + 1e-5:
                        raise ValueError(f"{name}: int16 round trip error {err} dB")
                elif v.dtype.kind == "f":
                    if not np.array_equal(a.astype("float32"), b.astype("float32"), equal_nan=True):
                        raise ValueError(f"{name}: float32 round trip differs")
                elif not np.array_equal(a, b):
                    raise ValueError(f"{name}: round trip differs")
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    os.replace(tmp, path)
    return {"file_size_bytes": path.stat().st_size, "int16_max_error_db": max_err}


# ---------------------------------------------------------------------------
# aggregate_cml (per day; a module-level task for cmlpipe.parallel)
# ---------------------------------------------------------------------------
def cml_input_hashes(cfg: Config, day) -> dict:
    """Stage hashes of the L1 and L2 files of a day (inputs of aggregate_cml)."""
    l2 = file_attrs(qc_runner.l2_day_path(cfg, day))
    return {"input_stage_hash_l1": existing_stage_hash(cfg.l1_cml_day_path(day)),
            "input_stage_hash_l2": l2.get("stage_hash"), "input_stats_stage_hash_l2": l2.get("stats_stage_hash")}


def db_vars(ds: xr.Dataset) -> set[str]:
    """The dB statistics (tl_*, rsl_*, tsl_*): the variables the int16 option applies to."""
    return {n for n in ds.data_vars if n.split("_")[0] in aggregate.CML_VARS and not n.startswith("n_")}


def aggregate_cml_day(cfg_data: dict, repo_root: str, config_file: str, day_text: str, force: bool = False) -> dict:
    """One day -> one file per resolution. Skipped only if every file has the current stage hash AND was made
    from the current L1 and L2 files (their stage hashes, and the L2 statistics hash)."""
    t0 = time.perf_counter()
    cfg = config_from_parts(cfg_data, repo_root, config_file)
    acfg = cfg["aggregate"]
    day = pd.Timestamp(day_text)
    l1_path, l2_path = cfg.l1_cml_day_path(day), qc_runner.l2_day_path(cfg, day)
    record = {"day": day_text, "messages": [], "files": {}}
    for p, what in ((l1_path, "ingest_cml"), (l2_path, "qc_flags")):
        if not p.exists():
            raise FileNotFoundError(f"{p} not found: run {what} for {day_text} first")
    stage_hash = cfg.stage_hash("aggregate_cml")
    inputs = cml_input_hashes(cfg, day)
    paths = {res: cml_path(cfg, res, day) for res in acfg["resolutions"]}
    if not force and all(
            (a := file_attrs(p)).get("stage_hash") == stage_hash and all(a.get(k) == v for k, v in inputs.items())
            for p in paths.values()):
        record["messages"].append(("info", f"{day_text}: aggregate_cml up to date (stage, L1 and L2 hashes), skipped"))
        return {**record, "status": "skipped_same_config", "runtime_s": round(time.perf_counter() - t0, 2)}

    native = cfg["cml"]["sampling_interval"]
    int16 = acfg["encoding"]["db_dtype"] == "int16"
    with xr.open_dataset(l1_path) as l1, xr.open_dataset(l2_path) as l2:
        l1, l2 = l1.load(), l2[[schema.FLAGS_VAR]].load()
    for res, path in paths.items():
        t1 = time.perf_counter()
        ds = aggregate.aggregate_cml_day(l1, l2, res, acfg, native)
        ds.attrs.update({**provenance_attrs(cfg, "aggregate_cml"), **inputs,
                         "title": f"L3 CML statistics, {res}, {day_text}", "processing_level": "L3_products",
                         "input_files": f"{display_path(l1_path, cfg.repo_root)} {display_path(l2_path, cfg.repo_root)}",
                         "db_encoding": ("int16, scale_factor %g dB, max round-trip error %g dB" %
                                         (acfg["encoding"]["int16_scale_factor"], acfg["encoding"]["int16_max_error_db"])
                                         if int16 else "float32")})
        schema.validate_l3(ds, "cml", res, start=day, n_bins=int(pd.Timedelta("1D") / pd.Timedelta(res)))
        written = write_verified(ds, path, acfg["encoding"], db_vars(ds) if int16 else set())
        with xr.open_dataset(path) as back:  # validate as on disk
            schema.validate_l3(back.load(), "cml", res, start=day, n_bins=ds.sizes["time"])
        record["files"][res] = {"path": display_path(path, cfg.repo_root), **written,
                                "runtime_s": round(time.perf_counter() - t1, 2)}
    record.update({"status": "written", "runtime_s": round(time.perf_counter() - t0, 2)})
    sizes = ", ".join(f"{r} {f['file_size_bytes'] / 1e6:.1f} MB" for r, f in record["files"].items())
    record["messages"].append(("info", f"{day_text}: aggregate_cml written ({record['runtime_s']:.1f} s): {sizes}"))
    return record


# ---------------------------------------------------------------------------
# aggregate_gauges
# ---------------------------------------------------------------------------
def _period_bins(cfg: Config, res: str) -> int:
    return int((cfg.period_end - cfg.period_start) / pd.Timedelta(res))


def _up_to_date(paths: list[Path], stage_hash: str, inputs: dict) -> bool:
    for p in paths:
        a = file_attrs(p)
        if a.get("stage_hash") != stage_hash or a.get("input_hashes") != json.dumps(inputs, sort_keys=True):
            return False
    return True


def run_aggregate_gauges(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """City and SMHI L1 -> sums per bin at every configured resolution (whole period)."""
    g = cfg["aggregate"]["gauges"]
    sources = {
        "city": (cfg.l1_gauges_dir / "city_1min.nc", cfg["city_gauges"]["interval"], g["city_resolutions"]),
        "smhi": (cfg.l1_gauges_dir / "smhi_15min.nc", cfg["smhi_gauge"]["interval"], g["smhi_resolutions"]),
    }
    stage_hash = cfg.stage_hash("aggregate_gauges")
    inputs = {net: existing_stage_hash(src[0]) for net, src in sources.items()}
    if None in inputs.values():
        raise FileNotFoundError("gauge L1 files missing: run ingest_gauges first")
    all_paths = [gauge_path(cfg, net, r) for net, (_, _, rs) in sources.items() for r in rs]
    if not force and _up_to_date(all_paths, stage_hash, inputs):
        log.info("aggregate_gauges: up to date, skipped (use --force to rebuild)")
        records.append({"step": "aggregate_gauges", "status": "skipped_same_config"})
        return
    for net, (src, native, resolutions) in sources.items():
        with xr.open_dataset(src) as l1:
            l1 = l1.load()
        for res in resolutions:
            t0 = time.perf_counter()
            ds = aggregate.aggregate_gauges(l1, res, native, cfg.period_start, cfg.period_end, g["min_valid_fraction"],
                                            peak_rate=(net == "city"), qc_worst=(net == "smhi"))
            ds.attrs.update({**provenance_attrs(cfg, "aggregate_gauges"), "processing_level": "L3_products",
                             "title": f"L3 {net} gauges, {res}", "network": net,
                             "input_file": display_path(src, cfg.repo_root),
                             "input_hashes": json.dumps(inputs, sort_keys=True),
                             "quantization_note": "native quantization (quantization_mm coord) is unchanged; "
                                                  "amounts are never disaggregated"})
            path = gauge_path(cfg, net, res)
            schema.validate_l3(ds, "gauges", res, start=cfg.period_start, n_bins=_period_bins(cfg, res))
            written = write_verified(ds, path, cfg["aggregate"]["encoding"])
            records.append({"step": "aggregate_gauges", "path": display_path(path, cfg.repo_root), "status": "written",
                            **written, "runtime_s": round(time.perf_counter() - t0, 2),
                            "total_mm": {str(i): round(float(np.nansum(ds.amount_mm.sel(id=i))), 3)
                                         for i in ds.id.values}})
            log.info("aggregate_gauges: wrote %s (%.2f MB)", path.name, written["file_size_bytes"] / 1e6)


# ---------------------------------------------------------------------------
# aggregate_radar
# ---------------------------------------------------------------------------
def radar_inputs(cfg: Config) -> dict:
    """Stage hashes / signatures of the radar L1 files and the geometry outputs."""
    files = sorted(cfg.l1_radar_dir.glob("radar_5min_*.nc"))
    if not files:
        raise FileNotFoundError("no radar L1 files: run ingest_radar first")
    wpath = geometry_runner.weights_path(cfg)
    if not wpath.exists():
        raise FileNotFoundError(f"{wpath} not found: run the geometry stage first")
    wa = file_attrs(wpath)
    side = json.loads(Path(str(cfg.tables_dir / "gauge_radar_pixel.csv") + ".meta.json").read_text(encoding="utf-8"))
    return {"radar_l1": {f.name: existing_stage_hash(f) for f in files},
            "geometry": {"stage_hash": wa.get("stage_hash"), "input_signature": wa.get("input_signature"),
                         "gauge_radar_pixel_signature": side.get("input_signature")}}


def radar_stamps(cfg: Config, rcfg: dict, weights: pd.DataFrame, grid: geometry.RadarGrid,
                 pix: pd.DataFrame) -> dict:
    """Per radar stamp: path statistics for every link and pixel / 3x3 values at every gauge (floored rates)."""
    files = sorted(cfg.l1_radar_dir.glob("radar_5min_*.nc"))
    parts = {k: [] for k in ("time", "R_mean", "R_std", "R_max", "wet_fraction", "valid_weight_fraction",
                             "R_pixel", "R_3x3")}
    n_zero = n_finite = 0
    ny, nx = grid.shape
    for f in files:
        with xr.open_dataset(f) as ds:
            dbz = ds.dbz.transpose("time", "y", "x").values
            rain = ds.rain_rate.transpose("time", "y", "x").values
            times, x, y = ds.time.values, ds.x.values, ds.y.values
            zr = f"Z = {float(ds.attrs['zr_a']):g} R^{float(ds.attrs['zr_b']):g} ({ds.attrs.get('zr_source', '')})"
        floored, nz, nf = aggregate.floor_rain(dbz, rain, rcfg["min_dbz"])
        n_zero, n_finite = n_zero + nz, n_finite + nf
        field = xr.DataArray(floored, dims=("time", "y", "x"), coords={"time": times, "y": y, "x": x})
        st = geometry.apply_weights(field, weights, threshold=rcfg["wet_threshold_mmh"],
                                    time_chunk=int(rcfg["time_chunk"]), grid=grid)
        ok = st["valid_weight_fraction"].values >= rcfg["min_valid_weight_fraction"]
        parts["time"].append(times)
        for src, dst in (("mean", "R_mean"), ("std", "R_std"), ("max", "R_max"), ("wet_fraction", "wet_fraction")):
            parts[dst].append(np.where(ok, st[src].values, np.nan))
        parts["valid_weight_fraction"].append(np.nan_to_num(st["valid_weight_fraction"].values, nan=0.0))
        iy, ix = pix.iy.to_numpy(), pix.ix.to_numpy()
        parts["R_pixel"].append(floored[:, iy, ix].astype("float64"))
        blocks = []
        for a, b in zip(iy, ix):
            blk = floored[:, max(a - 1, 0):min(a + 2, ny), max(b - 1, 0):min(b + 2, nx)].reshape(len(times), -1)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)  # all 9 cells NaN -> NaN
                blocks.append(np.nanmean(blk, axis=1))
        parts["R_3x3"].append(np.column_stack(blocks))
    out = {k: np.concatenate(v) for k, v in parts.items()}
    out["cml_id"] = st.cml_id.values
    out["n_zero"], out["n_finite"], out["zr"] = n_zero, n_finite, zr
    return out


def run_aggregate_radar(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """Radar L1 + path weights -> radar along links and at gauges, at 5 / 15 / 60 min (whole period)."""
    rcfg = cfg["aggregate"]["radar"]
    stage_hash = cfg.stage_hash("aggregate_radar")
    inputs = radar_inputs(cfg)
    paths = [radar_path(cfg, w, r) for w in ("along_links", "at_gauges") for r in rcfg["resolutions"]]
    if not force and _up_to_date(paths, stage_hash, inputs):
        log.info("aggregate_radar: up to date, skipped (use --force to rebuild)")
        records.append({"step": "aggregate_radar", "status": "skipped_same_config"})
        return
    t0 = time.perf_counter()
    weights, grid = geometry_runner.read_weights(geometry_runner.weights_path(cfg))
    pix = pd.read_csv(cfg.tables_dir / "gauge_radar_pixel.csv")
    if not pix.inside_radar_grid.all():
        raise ValueError("a gauge lies outside the radar grid; R_pixel is undefined for it")
    lg = pd.read_csv(cfg.tables_dir / "link_geometry.csv").set_index("cml_id")
    s = radar_stamps(cfg, rcfg, weights, grid, pix)
    log.info("aggregate_radar: per-stamp statistics done (%.1f s)", time.perf_counter() - t0)
    native = cfg["radar"]["interval"]
    zero_share = s["n_zero"] / max(s["n_finite"], 1)
    common = {
        **provenance_attrs(cfg, "aggregate_radar"), "processing_level": "L3_products",
        "input_hashes": json.dumps(inputs, sort_keys=True),
        "radar_stamp_meaning": "instantaneous sample (OpenMRG readme); the stamp at t is assigned to the bin "
                               "[t, t + dt); peak gauge-radar correlation at lag 0 supports this "
                               "(docs/data_notes/openmrg_reference_inspection.md)",
        "min_dbz_floor": rcfg["min_dbz"],
        "floor_rule": f"cells with dBZ < {rcfg['min_dbz']} are set to rain rate 0 before any averaging; NaN stays NaN",
        "floor_share_of_cell_steps_set_to_zero": zero_share,
        "floor_n_cell_steps_set_to_zero": s["n_zero"], "floor_n_finite_cell_steps": s["n_finite"],
        "min_valid_fraction": rcfg["min_valid_fraction"],
        "zr": s["zr"],
    }
    # stamps -> native bins of the period
    grid_vals = {}
    n_dropped = 0
    for k in ("R_mean", "R_std", "R_max", "wet_fraction", "valid_weight_fraction", "R_pixel", "R_3x3"):
        grid_vals[k], n_dropped = aggregate.stamps_to_grid(s["time"], s[k], cfg.period_start, cfg.period_end, native)
    common["n_stamps_outside_period_dropped"] = n_dropped
    link_valid = ~np.isnan(grid_vals["R_mean"])
    units = {"R_mean": "mm h-1", "R_std": "mm h-1", "R_max": "mm h-1", "wet_fraction": "1",
             "valid_weight_fraction": "1", "R_pixel": "mm h-1", "R_3x3": "mm h-1"}
    longs = {
        "R_mean": "path-weighted mean rain rate (floored), mean over the valid stamps of the bin",
        "R_std": "path-weighted std of the cell rates along the path, mean over the valid stamps",
        "R_max": "largest cell rate under the path, max over the valid stamps",
        "wet_fraction": f"share of the valid path weight in cells > {rcfg['wet_threshold_mmh']} mm/h, mean over valid stamps",
        "valid_weight_fraction": "share of the path weight with radar data, mean over ALL stamps of the bin "
                                 "(a missing stamp counts as 0)",
        "R_pixel": "rain rate (floored) of the radar cell containing the gauge, mean over valid stamps",
        "R_3x3": "mean rain rate (floored) of the 3 x 3 cells around the gauge cell (NaN cells excluded), "
                 "mean over valid stamps",
    }
    for res in rcfg["resolutions"]:
        t1 = time.perf_counter()
        spb = aggregate.samples_per_bin(res, native)
        times = aggregate.bin_starts(cfg.period_start, _period_bins(cfg, res), res)
        lv, n_valid = aggregate.aggregate_stamps(
            {k: grid_vals[k] for k in ("R_mean", "R_std", "R_max", "wet_fraction", "valid_weight_fraction")},
            link_valid, spb, rcfg["min_valid_fraction"],
            {"R_mean": "mean", "R_std": "mean", "R_max": "max", "wet_fraction": "mean", "valid_weight_fraction": "mean_all"})
        cids = s["cml_id"]
        links = xr.Dataset(
            {k: (("time", "cml_id"), v.astype("float32"), {"units": units[k], "long_name": longs[k]}) for k, v in lv.items()},
            coords={"time": times, "cml_id": cids,
                    "n_cells": ("cml_id", lg.loc[cids, "n_radar_cells"].to_numpy("int32"),
                                {"long_name": "radar cells crossed by the path"}),
                    "radar_coverage": ("cml_id", lg.loc[cids, "radar_coverage"].to_numpy("float64"),
                                       {"units": "1", "long_name": "fraction of the path inside the radar grid"})})
        links["n_valid_stamps"] = (("time", "cml_id"), n_valid.astype(aggregate.count_dtype(spb)),
                                   {"units": "1", "long_name": "stamps with valid_weight_fraction >= "
                                                               f"{rcfg['min_valid_weight_fraction']} in the bin"})
        links.attrs = {**aggregate.time_attrs(res, spb, "radar stamps are instants; stamp t -> bin [t, t + dt)"),
                       **common, "title": f"L3 radar along link paths, {res}",
                       "min_valid_weight_fraction": rcfg["min_valid_weight_fraction"],
                       "wet_threshold_mmh": rcfg["wet_threshold_mmh"],
                       "valid_stamp_rule": f"a stamp is valid for a link if valid_weight_fraction >= "
                                           f"{rcfg['min_valid_weight_fraction']}; a bin is NaN if fewer than "
                                           f"{aggregate.min_count(rcfg['min_valid_fraction'], spb)} of {spb} stamps are valid"}

        gv = {}
        nv = {}
        for k in ("R_pixel", "R_3x3"):
            out, n = aggregate.aggregate_stamps({k: grid_vals[k]}, ~np.isnan(grid_vals[k]), spb,
                                                rcfg["min_valid_fraction"], {k: "mean"})
            gv[k], nv[k] = out[k], n
        gauges = xr.Dataset(
            {k: (("time", "id"), v.astype("float32"), {"units": units[k], "long_name": longs[k]}) for k, v in gv.items()},
            coords={"time": times, "id": pix.gauge_id.to_numpy().astype(str),
                    "network": ("id", pix.network.to_numpy().astype(str)),
                    "iy": ("id", pix.iy.to_numpy("int32")), "ix": ("id", pix.ix.to_numpy("int32")),
                    "lon": ("id", pix.lon.to_numpy()), "lat": ("id", pix.lat.to_numpy())})
        for k in ("R_pixel", "R_3x3"):
            gauges[f"n_valid_stamps_{k.split('_')[1]}"] = (("time", "id"), nv[k].astype(aggregate.count_dtype(spb)),
                                                           {"units": "1", "long_name": f"valid stamps of {k} in the bin"})
        gauges.attrs = {**aggregate.time_attrs(res, spb, "radar stamps are instants; stamp t -> bin [t, t + dt)"),
                        **common, "title": f"L3 radar at gauges, {res}"}
        for kind, ds, what in (("radar_links", links, "along_links"), ("radar_gauges", gauges, "at_gauges")):
            schema.validate_l3(ds, kind, res, start=cfg.period_start, n_bins=_period_bins(cfg, res))
            path = radar_path(cfg, what, res)
            written = write_verified(ds, path, cfg["aggregate"]["encoding"])
            records.append({"step": "aggregate_radar", "path": display_path(path, cfg.repo_root), "status": "written",
                            **written, "runtime_s": round(time.perf_counter() - t1, 2)})
            log.info("aggregate_radar: wrote %s (%.1f MB)", path.name, written["file_size_bytes"] / 1e6)
    records.append({"step": "aggregate_radar", "status": "note", "floor_share_set_to_zero": zero_share,
                    "n_stamps_outside_period_dropped": n_dropped, "runtime_s": round(time.perf_counter() - t0, 2)})
