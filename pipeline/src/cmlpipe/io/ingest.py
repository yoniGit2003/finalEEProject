"""Ingest stages: raw OpenMRG files -> L1 (ingest_cml per day; ingest_gauges, ingest_met, ingest_radar).

ingest_cml_day is a module-level function so that worker processes (Windows "spawn") can import it: it rebuilds
the Config from its parts, does ONE day, and returns a manifest record (see cmlpipe.parallel). Log messages are
returned in record["messages"] as (level, text) and logged by the parent, because a worker process does not share
the parent's log handlers.
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from cmlpipe.config import Config
from cmlpipe.io import openmrg, schema
from cmlpipe.io import openmrg_reference as ref
from cmlpipe.io.netcdf import existing_stage_hash
from cmlpipe.parallel import config_from_parts, run_day_task
from cmlpipe.provenance import display_path


def day_statistics(ds: xr.Dataset) -> dict:
    """Numbers for the manifest, computed on an L1 CML day as re-read from disk."""
    present = ds.raw_sublink_id.values != schema.MISSING_SUBLINK
    stats = {"n_time": int(ds.sizes["time"])}
    for name in ("rsl", "tsl"):
        is_nan = np.isnan(ds[name].values)
        stats[f"pct_nan_{name}"] = round(100.0 * float(is_nan.mean()), 4)
        # count only sub-links that exist (absent ones are NaN by construction)
        stats[f"n_sublinks_all_nan_{name}"] = int((is_nan.all(axis=0) & present).sum())
    for key in ("n_duplicates_removed", "n_duplicates_with_different_values", "n_offgrid_steps_dropped",
                "n_steps_inserted", "n_sentinel_rsl", "n_sentinel_tsl"):
        stats[key] = int(ds.attrs[key])
    return stats


def ingest_cml_day(cfg_data: dict, repo_root: str, config_file: str, day_text: str, force: bool = False) -> dict:
    """Raw cml.nc -> one L1 day file, validated as read back from disk."""
    t0 = time.perf_counter()
    cfg = config_from_parts(cfg_data, repo_root, config_file)
    day = pd.Timestamp(day_text)
    path = cfg.l1_cml_day_path(day)
    record = {"day": day_text, "path": display_path(path, cfg.repo_root), "messages": []}
    if not force and existing_stage_hash(path) == cfg.stage_hash("ingest_cml"):
        record["messages"].append(("info", f"{day_text}: output exists with the same stage hash, skipped "
                                           "(use --force to rebuild)"))
        return {**record, "status": "skipped_same_config", "runtime_s": round(time.perf_counter() - t0, 2)}
    meta = openmrg.load_cml_metadata(cfg)
    ds = openmrg.load_cml_day(cfg, day, meta=meta)
    written = openmrg.write_l1_day(ds, path, cfg)
    with xr.open_dataset(path) as on_disk:  # validate the file as it is on disk, not the in-memory copy
        on_disk.load()
        report = schema.validate_cml_l1(on_disk, rules=cfg["validation"], interval=cfg.sampling_interval,
                                        raise_on_error=True)
        stats = day_statistics(on_disk)
    record["messages"] += [("warning", f"{day_text}: {w}") for w in report.warnings]
    record.update({"status": "written", **stats, "file_size_bytes": written["file_size_bytes"],
                   "roundtrip": written["roundtrip"], "runtime_s": round(time.perf_counter() - t0, 2),
                   "validator_warnings": report.warnings})
    record["messages"].append(("info", f"{day_text}: written {path.name} ({record['file_size_bytes'] / 1e6:.1f} MB, "
                                       f"{record['runtime_s']:.1f} s, NaN rsl {stats['pct_nan_rsl']:.2f}%, "
                                       f"tsl {stats['pct_nan_tsl']:.2f}%)"))
    return record


def ingest_cml(cfg: Config, days: list[pd.Timestamp], force: bool, log, records: list[dict], workers: int = 1) -> None:
    """Run the ingest_cml stage for the given days (one manifest record per day).

    tables/sublink_id_map.csv is rewritten only when its content changes, so an up-to-date run writes nothing.
    """
    if cfg.dataset_name != "openmrg":
        raise NotImplementedError(f"ingest_cml has no adapter for dataset {cfg.dataset_name!r}")
    meta = openmrg.load_cml_metadata(cfg)
    map_path = cfg.tables_dir / "sublink_id_map.csv"
    new = meta[["cml_id", "sublink_id", "raw_sublink_id", "direction"]]
    if not map_path.exists() or not pd.read_csv(map_path).astype(str).equals(new.astype(str).reset_index(drop=True)):
        log.info("Wrote %s", openmrg.write_sublink_id_map(meta, cfg))
    run_day_task(ingest_cml_day, cfg, days, workers, force, log, records)


def file_statistics(ds: xr.Dataset) -> dict:
    """Numbers for the manifest, computed on a reference file as re-read from disk."""
    stats = {
        "shape": {k: int(v) for k, v in ds.sizes.items()},
        "time_first": str(pd.Timestamp(ds.time.values[0])),
        "time_last": str(pd.Timestamp(ds.time.values[-1])),
        "pct_nan": {name: round(100.0 * float(np.isnan(v.values).mean()), 4)
                    for name, v in ds.data_vars.items() if v.dtype.kind == "f"},
    }
    # the adapter's counts (NA converted, duplicates, inserted steps, negatives, ...)
    stats.update({k: v.item() if hasattr(v, "item") else v for k, v in ds.attrs.items() if k.startswith("n_")})
    return stats


def write_reference_file(step: str, path: Path, build, validate, cfg: Config, force: bool, log,
                         records: list[dict]) -> bool:
    """Build one L1 reference file, write it, validate it as read back from disk, and record it.

    build    : function () -> xr.Dataset (only called if the file must be (re)built)
    validate : function (ds) -> ValidationReport (raises on schema errors)
    Returns True if the file was written, False if it was skipped (same stage hash).
    """
    t0 = time.perf_counter()
    record = {"step": step, "path": display_path(path, cfg.repo_root)}
    if not force and existing_stage_hash(path) == cfg.stage_hash(step):
        log.info("%s: %s exists with the same stage hash, skipped (use --force to rebuild)", step, path.name)
        records.append({**record, "status": "skipped_same_config"})
        return False
    try:
        ds = build()
        written = ref.write_l1_reference(ds, path, cfg, step)
        with xr.open_dataset(path) as on_disk:
            on_disk.load()
            report = validate(on_disk)
            stats = file_statistics(on_disk)
    except Exception as err:  # noqa: BLE001 - record the failure, then stop the run loudly
        log.exception("%s: %s FAILED", step, path.name)
        records.append({**record, "status": "failed", "error": str(err),
                        "runtime_s": round(time.perf_counter() - t0, 2)})
        raise
    for w in report.warnings:
        log.warning("%s: %s: %s", step, path.name, w)
    record.update({"status": "written", **stats, "file_size_bytes": written["file_size_bytes"],
                   "runtime_s": round(time.perf_counter() - t0, 2), "validator_warnings": report.warnings})
    log.info("%s: written %s (%.2f MB, %.1f s, NaN %s)", step, path.name, record["file_size_bytes"] / 1e6,
             record["runtime_s"], stats["pct_nan"])
    records.append(record)
    return True


def ingest_gauges(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """City (1 min) and SMHI (15 min) gauges -> two L1 files, then tables/gauges.csv from those files."""
    rules = cfg["reference_validation"]
    files = {
        "city": (cfg.l1_gauges_dir / "city_1min.nc", ref.load_city_gauges, cfg["city_gauges"]["interval"]),
        "smhi": (cfg.l1_gauges_dir / "smhi_15min.nc", ref.load_smhi_gauge, cfg["smhi_gauge"]["interval"]),
    }
    any_written = False
    for path, loader, interval in files.values():
        any_written |= write_reference_file(
            "ingest_gauges", path, lambda loader=loader: loader(cfg),
            lambda ds, interval=interval: schema.validate_gauge_l1(ds, rules, interval), cfg, force, log, records)

    table_path = cfg.tables_dir / "gauges.csv"
    if any_written or force or not table_path.exists():
        # Built from the files on disk, so the table always describes what was written.
        datasets = []
        for path, _, _ in files.values():
            with xr.open_dataset(path) as ds:
                datasets.append(ds.load())
        path = ref.write_gauge_table(ref.gauge_table(datasets), cfg)
        log.info("ingest_gauges: wrote %s", path)
        records.append({"step": "ingest_gauges", "path": display_path(path, cfg.repo_root), "status": "written",
                        "n_rows": sum(ds.sizes["id"] for ds in datasets)})


def ingest_met(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """SMHI hourly weather -> one L1 file."""
    rules, interval = cfg["reference_validation"], cfg["smhi_met"]["interval"]
    write_reference_file("ingest_met", cfg.l1_met_dir / "smhi_met_1h.nc", lambda: ref.load_smhi_weather(cfg),
                         lambda ds: schema.validate_met_l1(ds, rules, interval), cfg, force, log, records)


def ingest_radar(cfg: Config, months: list[tuple[int, int]], force: bool, log, records: list[dict]) -> None:
    """Radar -> one L1 file per calendar month."""
    rules, interval = cfg["reference_validation"], cfg["radar"]["interval"]
    n_outside = ref.radar_steps_outside_months(cfg)
    if n_outside:
        log.info("ingest_radar: %d raw radar steps fall outside the period and are in no month file", n_outside)
    for year, month in months:
        write_reference_file(
            "ingest_radar", cfg.l1_radar_dir / f"radar_5min_{year}-{month:02d}.nc",
            lambda year=year, month=month: ref.load_radar_month(cfg, year, month),
            lambda ds: schema.validate_radar_l1(ds, rules, interval), cfg, force, log, records)
    records.append({"step": "ingest_radar", "status": "note", "n_raw_steps_outside_period": n_outside})
