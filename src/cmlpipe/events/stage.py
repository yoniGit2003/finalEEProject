"""Stages `events`, `split` and `event_links` (pure logic in cmlpipe.events.core).

events      reference series (radar domain + City gauges) -> states -> events -> dry periods, the calibration
            grid, and the per-sublink dry-weather statistics. The event DEFINITION (build_reference,
            define_events) reads no CML and no SMHI data; only the dry statistics read the CML core file.
split       frozen train/val/test split of events and train/test split of links.
event_links one row per (event, link): radar along the path, matched gauge, CML data quality.

Outputs: L3_products/events/network_reference_5min.nc and tables/event_calibration.csv, events.csv,
dry_periods.csv, sublink_dry_stats.csv, split_events.csv, split_links.csv, event_links.csv (each with a
.meta.json sidecar). CML is read from ONE consolidated file per resolution
(L3_products/cml/consolidated/cml_<res>_core.nc, built here with netCDF4 from the daily L3 files) because
opening 92 daily files with xarray is slow on Windows and dask is not installed.
"""

from __future__ import annotations

import hashlib
import json
import time
import warnings
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr

from cmlpipe import aggregate
from cmlpipe.aggregate import stage as aggregate_runner
from cmlpipe import events
from cmlpipe.io import schema
from cmlpipe.config import Config
from cmlpipe.provenance import display_path
from cmlpipe.io.netcdf import read_core_file
from cmlpipe.provenance import provenance_attrs, write_sidecar

TABLES = {
    "calibration": "event_calibration.csv", "events": "events.csv", "dry_periods": "dry_periods.csv",
    "dry_stats": "sublink_dry_stats.csv", "split_events": "split_events.csv", "split_links": "split_links.csv",
    "event_links": "event_links.csv",
}


def table_path(cfg: Config, key: str) -> Path:
    """Path of a step 6 table by key (see TABLES)."""
    return cfg.tables_dir / TABLES[key]


def reference_path(cfg: Config) -> Path:
    """Path of network_reference_<res>.nc."""
    return cfg.output_root / cfg["paths"]["l3_events_dir"] / f"network_reference_{cfg['events']['resolution']}.nc"


def core_path(cfg: Config, res: str) -> Path:
    """Path of the consolidated CML core file of a resolution."""
    return cfg.output_root / cfg["paths"]["l3_cml_consolidated"] / f"cml_{res}_core.nc"


def _sidecar_meta(path: Path) -> dict:
    side = Path(str(path) + ".meta.json")
    return json.loads(side.read_text(encoding="utf-8")) if side.exists() else {}


def _write_table(cfg: Config, key: str, df: pd.DataFrame, stage: str, description: str, extra: dict | None = None) -> Path:
    path = table_path(cfg, key)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, index=False)
    tmp.replace(path)
    write_sidecar(path, cfg, stage, extra={"description": description, "n_rows": len(df), "columns": list(df.columns),
                                           **(extra or {})})
    return path


# ---------------------------------------------------------------------------
# Reference series (NO CML, NO SMHI)
# ---------------------------------------------------------------------------
def period_grid(cfg: Config, res: str) -> pd.DatetimeIndex:
    """Bin starts of the whole period at `res`."""
    n = int((cfg.period_end - cfg.period_start) / pd.Timedelta(res))
    return pd.DatetimeIndex(aggregate.bin_starts(cfg.period_start, n, res))


def build_reference(cfg: Config) -> tuple[pd.DataFrame, dict]:
    """Per 5 min bin of the period: radar over the link domain and the 10 City gauges (columns of the reference
    file). Reads radar L1, the City gauge L3 product and link_geometry.csv. Never CML, never SMHI."""
    ecfg = cfg["events"]
    res = ecfg["resolution"]
    lg = pd.read_csv(cfg.tables_dir / "link_geometry.csv")
    files = sorted(cfg.l1_radar_dir.glob("radar_5min_*.nc"))
    if not files:
        raise FileNotFoundError("no radar L1 files: run ingest_radar first")
    parts, times, mask, info = [], [], None, {}
    for f in files:
        with xr.open_dataset(f) as ds:
            if mask is None:
                mask = events.domain_mask(ds.x_3006.values, ds.y_3006.values, lg, ecfg["domain_buffer_m"])
                cell_area = abs(float(ds.x.values[1] - ds.x.values[0]) * float(ds.y.values[1] - ds.y.values[0]))
                info = {"n_domain_cells": int(mask.sum()), "n_grid_cells": int(mask.size),
                        "domain_area_km2_radar_projection": float(mask.sum() * cell_area / 1e6)}
            dbz, rain = ds.dbz.values, ds.rain_rate.values
            times.append(ds.time.values)
        floored, _, _ = aggregate.floor_rain(dbz, rain, cfg["aggregate"]["radar"]["min_dbz"])
        parts.append(events.radar_domain_series(floored, mask, ecfg["wet_threshold_mmh"]))
    stamps = np.concatenate(times)
    radar = {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}
    grid = period_grid(cfg, res)
    ref = pd.DataFrame(index=grid)
    for k, v in radar.items():
        g, dropped = aggregate.stamps_to_grid(stamps, v[:, None], cfg.period_start, cfg.period_end, res)
        ref[k] = g[:, 0]
    ref["radar_valid_cell_fraction"] = ref.radar_valid_cell_fraction.fillna(0.0)  # no stamp = no radar data
    info["n_radar_stamps_outside_period"] = dropped

    with xr.open_dataset(aggregate_runner.gauge_path(cfg, "city", res)) as g:
        g = g.load()
    if not np.array_equal(g.time.values, grid.values):
        raise ValueError("City gauge product and reference grid differ")
    gs = events.gauge_series(g.amount_mm.transpose("time", "id").values.astype("float64"),
                             g.rate_mm_h.transpose("time", "id").values.astype("float64"), ecfg["city_gauge_wet_amount_mm"])
    for k, v in gs.items():
        ref[k] = v
    ref.index.name = "time"
    return ref, {**info, "domain_mask": mask}


def define_events(cfg: Config, ref: pd.DataFrame) -> tuple[pd.DataFrame, list[dict], int, np.ndarray]:
    """Classify the reference bins and find the events. Returns (ref with state, events [i0, i1), n_dropped, dry)."""
    ecfg = cfg["events"]
    step = pd.Timedelta(ecfg["resolution"])
    ref = events.classify(ref, ecfg)
    st = ref.state.to_numpy()
    ev, dropped = events.find_events(st, step, pd.Timedelta(hours=ecfg["min_dry_gap_h"]),
                                     pd.Timedelta(minutes=ecfg["min_event_duration_min"]))
    dry = events.certified_dry(st, int(round(pd.Timedelta(minutes=ecfg["dry_margin_min"]) / step)))
    eid = np.full(len(ref), "", dtype=object)
    for e in ev:
        e["event_id"] = events.event_id(ref.index[e["i0"]])
        eid[e["i0"]:e["i1"]] = e["event_id"]
    ref["event_id"] = eid
    return ref, ev, dropped, dry


# ---------------------------------------------------------------------------
# CML core file (one per resolution)
# ---------------------------------------------------------------------------
def build_core(cfg: Config, res: str, force: bool = False, log=None, variables: list[str] | None = None) -> Path:
    """Consolidate a few CML variables from the 92 daily L3 files into ONE file per resolution (netCDF4 reads).

    variables: default events.cml_core_vars. Skipped when the file exists and was built from exactly the current
    daily files (stage hashes, creation times) and the same variable list.
    """
    vars_ = list(variables) if variables is not None else cfg["events"]["cml_core_vars"]
    days = cfg.period_days()
    paths = [aggregate_runner.cml_path(cfg, res, d) for d in days]
    missing = [p.name for p in paths if not p.exists()]
    if missing:
        raise FileNotFoundError(f"{len(missing)} L3 CML days missing at {res} (e.g. {missing[0]}): run aggregate_cml --all")
    sig = hashlib.sha256()
    for p in paths:
        a = aggregate_runner.file_attrs(p)
        sig.update(f"{p.name}:{a.get('stage_hash')}:{a.get('created_utc')};".encode())
    sig.update(json.dumps(vars_).encode())
    signature = sig.hexdigest()
    out = core_path(cfg, res)
    if not force and aggregate_runner.file_attrs(out).get("input_signature") == signature:
        return out
    t0 = time.perf_counter()
    arrays = {v: [] for v in vars_}
    with netCDF4.Dataset(paths[0]) as nc:
        cml_ids = nc["cml_id"][:].astype("int64")
        sub = [str(s) for s in nc["sublink_id"][:]]
        stage = nc.getncattr("stage_hash")
    for p in paths:
        with netCDF4.Dataset(p) as nc:
            nc.set_auto_mask(False)
            for v in vars_:
                arrays[v].append(nc[v][:])
    data = {v: np.concatenate(a) for v, a in arrays.items()}
    ds = xr.Dataset({v: (schema.DIMS, a) for v, a in data.items()},
                    coords={"time": period_grid(cfg, res).values, "sublink_id": sub, "cml_id": cml_ids})
    for v in vars_:
        with xr.open_dataset(paths[0]) as f:
            ds[v].attrs = {k: val for k, val in f[v].attrs.items() if k in ("units", "long_name", "flag_masks", "flag_meanings")}
    ds.attrs = {**aggregate.time_attrs(res, int(pd.Timedelta(res) / pd.Timedelta(cfg["cml"]["sampling_interval"]))),
                "title": f"Consolidated L3 CML core variables, {res}, whole period",
                "stage": "cml_core", "input_signature": signature, "input_aggregate_cml_stage_hash": stage,
                "comment": "a copy of variables of the daily L3 files (values unchanged), for fast reading",
                **{k: v for k, v in provenance_attrs(cfg, "events").items() if k not in ("stage", "stage_hash")}}
    aggregate_runner.write_verified(ds, out, cfg["aggregate"]["encoding"])
    if log:
        log.info("events: consolidated %d days of %s CML into %s (%.1f s)", len(paths), res, out.name, time.perf_counter() - t0)
    return out


def read_core(cfg: Config, res: str, variables: list[str]) -> tuple[pd.DatetimeIndex, np.ndarray, list[str], dict]:
    """(time, cml_ids, sublink_ids, {var: (time, sublink, cml) array}) from the consolidated file of `res`."""
    return read_core_file(core_path(cfg, res), variables)


# ---------------------------------------------------------------------------
# Event table
# ---------------------------------------------------------------------------
def read_smhi(cfg: Config, res: str = "15min") -> pd.Series:
    """SMHI GbgA amounts. Used ONLY as information in events.csv and in the acceptance check, never to define events."""
    with xr.open_dataset(aggregate_runner.gauge_path(cfg, "smhi", res)) as s:
        return s.amount_mm.sel(id=cfg["smhi_gauge"]["id"]).to_series().astype("float64")


def read_city(cfg: Config, res: str) -> xr.Dataset:
    """City gauge amounts and rates at `res` (L3)."""
    with xr.open_dataset(aggregate_runner.gauge_path(cfg, "city", res)) as g:
        return g[["amount_mm", "rate_mm_h"]].load()


def _window_sum_15(series: pd.Series, start, end) -> float:
    """Sum of the 15 min bins that overlap [start, end)."""
    s = series[(series.index >= pd.Timestamp(start).floor("15min")) & (series.index < pd.Timestamp(end))]
    return float(s.sum(min_count=1))


def _nanmax(a, axis=0):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmax(a, axis=axis)


def _nanmean(a, axis=0):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(a, axis=axis)


def overlap_minutes(periods: pd.DataFrame, start, end) -> float:
    """Total minutes of [start, end) covered by the periods (start/end columns, datetime)."""
    if not len(periods):
        return 0.0
    s = periods.start.where(periods.start > start, start)
    e = periods.end.where(periods.end < end, end)
    return float(((e - s) / pd.Timedelta("1min")).clip(lower=0).sum())


def event_table(cfg: Config, ref: pd.DataFrame, ev: list[dict]) -> pd.DataFrame:
    """One row per event: timing, dry context, radar and City gauge summaries, SMHI (info), gaps, weather."""
    ecfg = cfg["events"]
    step = pd.Timedelta(ecfg["resolution"])
    hours = step / pd.Timedelta("1h")
    st = ref.state.to_numpy()
    cap = int(pd.Timedelta(hours=ecfg["dry_context_cap_h"]) / step)
    gap_bins = int(np.ceil(pd.Timedelta(hours=ecfg["min_dry_gap_h"]) / step))
    city5, city1 = read_city(cfg, "5min"), read_city(cfg, "1min")
    gids = [str(i) for i in city5.id.values]
    a5 = city5.amount_mm.transpose("time", "id").values.astype("float64")
    r5 = city5.rate_mm_h.transpose("time", "id").values.astype("float64")
    t1 = pd.DatetimeIndex(city1.time.values)
    r1 = city1.rate_mm_h.transpose("time", "id").values.astype("float64")
    smhi = read_smhi(cfg)
    gaps = pd.read_csv(cfg.tables_dir / "network_gaps.csv")
    gaps["start"], gaps["end"] = pd.to_datetime(gaps.start), pd.to_datetime(gaps.end)  # also for an empty table
    with xr.open_dataset(cfg.l1_met_dir / "smhi_met_1h.nc") as m:
        met = m[["temperature", "humidity", "wind_speed", "wind_gust"]].isel(id=0).drop_vars("id").to_dataframe()
    dm_all = ref.domain_mean_rain_mmh.to_numpy()
    waf_all = ref.radar_wet_area_fraction.to_numpy()
    rv_all = ref.radar_valid.to_numpy()
    rows = []
    for k, e in enumerate(ev):
        i0, i1 = e["i0"], e["i1"]
        start, end = ref.index[i0], ref.index[i1 - 1] + step
        sl = slice(i0, i1)
        cx = events.context(st, ev, k, cap)
        dm = dm_all[sl]
        tot = np.where(np.isfinite(a5[sl]).any(axis=0), np.nansum(a5[sl], axis=0), np.nan)
        m5 = _nanmax(r5[sl])
        m1 = _nanmax(r1[(t1 >= start) & (t1 < end)])
        ov = overlap_minutes(gaps, start, end)
        # weather: hourly instants from the hour before the start to the hour after the end (never interpolated)
        wx = met[(met.index >= start.floor("h")) & (met.index <= end.ceil("h"))]
        row = {
            "event_id": e["event_id"], "start": start, "end": end, "duration_h": (end - start) / pd.Timedelta("1h"),
            "gap_before_h": cx["gap_before"] * hours, "dry_before_h": cx["dry_before"] * hours,
            "gap_after_h": cx["gap_after"] * hours, "dry_after_h": cx["dry_after"] * hours, "merged_from": e["merged_from"],
            "touches_period_edge": events.touches_edge(i0, i1, len(ref), gap_bins),
            "radar_peak_domain_rain_mmh": float(_nanmax(dm)), "radar_mean_domain_rain_mmh": float(_nanmean(dm)),
            "radar_max_wet_area_fraction": float(_nanmax(waf_all[sl])),
            "radar_domain_total_mm": float(np.nansum(dm) * hours),
        }
        row.update({f"city_total_{g}_mm": float(v) for g, v in zip(gids, tot)})
        row.update({
            "city_total_mean_mm": float(_nanmean(tot)), "city_total_max_mm": float(_nanmax(tot)),
            "city_max_5min_rate_mmh": float(_nanmax(m5)),
            "city_max_5min_gauge": gids[int(np.nanargmax(m5))] if np.isfinite(m5).any() else "",
            "city_max_1min_rate_mmh": float(_nanmax(m1)),
            "city_max_1min_gauge": gids[int(np.nanargmax(m1))] if np.isfinite(m1).any() else "",
            "smhi_total_mm_info": _window_sum_15(smhi, start, end),
            "n_unknown_bins": int((st[sl] == events.STATE_CODES["unknown"]).sum()),
            "radar_missing_fraction": float((~rv_all[sl]).mean()),
            "network_gap_minutes": float(ov),
        })
        for v in ("temperature", "humidity", "wind_speed", "wind_gust"):
            row.update({f"{v}_min": wx[v].min(), f"{v}_max": wx[v].max(), f"{v}_mean": wx[v].mean()})
        row["n_weather_hours"] = len(wx)
        row["month"] = start.month
        rows.append(row)
    df = pd.DataFrame(rows)
    if len(df):
        r = df.city_max_5min_rate_mmh.fillna(0.0)  # a radar-only event (no gauge rain) is "light"
        q = np.quantile(r, [1 / 3, 2 / 3])
        df["intensity_class"] = np.where(r <= q[0], "light", np.where(r <= q[1], "moderate", "heavy"))
        df.attrs["intensity_terciles_mmh"] = [float(x) for x in q]
    return df


def dry_periods_table(ref: pd.DataFrame, dry: np.ndarray, step: pd.Timedelta, min_h: float) -> pd.DataFrame:
    """Runs of certified-dry bins lasting at least min_h hours."""
    rows = []
    for a, b in events.runs(dry):
        if (b - a) * step >= pd.Timedelta(hours=min_h):
            rows.append({"start": ref.index[a], "end": ref.index[b - 1] + step, "duration_h": (b - a) * step / pd.Timedelta("1h"),
                         "radar_valid_fraction": float(ref.radar_valid_cell_fraction.iloc[a:b].mean())})
    return pd.DataFrame(rows, columns=["start", "end", "duration_h", "radar_valid_fraction"])


def in_periods(times: pd.DatetimeIndex, periods: pd.DataFrame) -> np.ndarray:
    """True for times inside any [start, end) of the periods table."""
    m = np.zeros(len(times), bool)
    for s, e in zip(pd.to_datetime(periods.start), pd.to_datetime(periods.end)):
        m |= (times >= s) & (times < e)
    return m


def dry_stats_table(cfg: Config, periods: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """sublink_dry_stats.csv from the consolidated 5 min CML core file, dry-period bins only."""
    res = cfg["events"]["resolution"]
    t0 = time.perf_counter()
    times, cids, subs, d = read_core(cfg, res, ["tl_mean", "tl_std", "rsl_mean", "rsl_std"])
    read_s = time.perf_counter() - t0
    dry = in_periods(times, periods)
    core = {k: v.reshape(v.shape[0], -1).astype("float64") for k, v in d.items()}  # columns: sublink-major
    stats = events.dry_stats(times, dry, core, min_days=cfg["events"]["dry_stats_min_days"])
    df = pd.DataFrame({"cml_id": np.tile(cids, len(subs)), "sublink_id": np.repeat(subs, len(cids)), **stats})
    return df.sort_values(["cml_id", "sublink_id"]).reset_index(drop=True), read_s


# ---------------------------------------------------------------------------
# Stage: events
# ---------------------------------------------------------------------------
def events_input_signature(cfg: Config) -> str:
    """sha256 over the inputs of the events stage."""
    h = hashlib.sha256()
    for f in sorted(cfg.l1_radar_dir.glob("radar_5min_*.nc")):
        h.update(str(aggregate_runner.file_attrs(f).get("stage_hash")).encode())
    for net, r in (("city", "5min"), ("city", "1min"), ("smhi", "15min")):
        a = aggregate_runner.file_attrs(aggregate_runner.gauge_path(cfg, net, r))
        h.update(f"{a.get('stage_hash')}{a.get('created_utc')}".encode())
    for t in ("link_geometry.csv", "network_gaps.csv"):
        h.update((cfg.tables_dir / t).read_bytes())
    h.update(str(aggregate_runner.file_attrs(cfg.l1_met_dir / "smhi_met_1h.nc").get("stage_hash")).encode())
    return h.hexdigest()


def write_reference(cfg: Config, ref: pd.DataFrame, info: dict) -> Path:
    """Write the network reference series NetCDF."""
    res = cfg["events"]["resolution"]
    floats = {"radar_wet_area_fraction": "1", "radar_valid_cell_fraction": "1", "domain_mean_rain_mmh": "mm h-1",
              "max_city_rate_mmh": "mm h-1"}
    ints = ["radar_n_wet_cells", "n_city_gauges_wet", "n_city_gauges_valid"]
    ds = xr.Dataset(coords={"time": ref.index.values})
    for c, u in floats.items():
        ds[c] = ("time", ref[c].to_numpy().astype("float32"), {"units": u})
    for c in ints:
        ds[c] = ("time", np.nan_to_num(ref[c].to_numpy()).astype("int16"), {"units": "1"})
    ds["state"] = ("time", ref.state.to_numpy().astype("int8"),
                   {"flag_values": np.array(list(events.STATE_CODES.values()), "int8"),
                    "flag_meanings": " ".join(events.STATE_CODES), "long_name": "reference state of the bin"})
    ds["event_id"] = ("time", ref.event_id.to_numpy().astype(str), {"long_name": "event id ('' outside events)"})
    ds.attrs = {**aggregate.time_attrs(res, 1, "radar stamp t -> bin [t, t + dt); City gauge sums of the bin"),
                **provenance_attrs(cfg, "events"),
                "title": "Network reference series (radar over the link domain + City gauges)",
                **{k: v for k, v in info.items() if k != "domain_mask"},
                "definition": "see cmlpipe.events.core module docstring; thresholds in config events.*"}
    path = reference_path(cfg)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    ds.to_netcdf(tmp, encoding={"time": aggregate_runner.TIME_ENCODING})
    tmp.replace(path)
    return path


def run_events(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """Stage events: reference series, events, dry periods, calibration grid, dry statistics."""
    t0 = time.perf_counter()
    ecfg = cfg["events"]
    res = ecfg["resolution"]
    step = pd.Timedelta(res)
    stage_hash = cfg.stage_hash("events")
    sig = events_input_signature(cfg)
    outs = [table_path(cfg, k) for k in ("calibration", "events", "dry_periods", "dry_stats")]
    core = build_core(cfg, res, log=log)  # cheap when up to date
    core_sig = aggregate_runner.file_attrs(core).get("input_signature")
    metas = [_sidecar_meta(p) for p in outs]
    if not force and reference_path(cfg).exists() and all(
            m.get("stage_hash") == stage_hash and m.get("input_signature") == sig and m.get("core_signature") == core_sig
            for m in metas):
        log.info("events: up to date, skipped (use --force to rebuild)")
        records.append({"step": "events", "status": "skipped_same_config"})
        return
    # --- definition: radar + City gauges only ---
    ref, info = build_reference(cfg)
    ref, ev, dropped, dry = define_events(cfg, ref)
    t_def = time.perf_counter() - t0
    extra = {"input_signature": sig, "core_signature": core_sig}
    cal = events.sensitivity(ref, ecfg, step)
    _write_table(cfg, "calibration", cal, "events", "events and event hours for a grid of min_wet_area_fraction x "
                 "min_dry_gap_h (all other thresholds as configured)", extra)
    et = event_table(cfg, ref, ev)
    events.validate_events(et)
    _write_table(cfg, "events", et, "events", "one row per rain event (defined from radar + City gauges only)",
                 {**extra, "n_dropped_short": dropped, "intensity_terciles_mmh": et.attrs.get("intensity_terciles_mmh")})
    annotate_events(cfg)  # add split_event if a frozen split exists (one code path, so the bytes are reproducible)
    dp = dry_periods_table(ref, dry, step, ecfg["dry_min_duration_h"])
    events.validate_periods(dp)
    _write_table(cfg, "dry_periods", dp, "events", "certified dry periods (radar valid, nothing wet, more than "
                 "dry_margin_min from wet bins) lasting >= dry_min_duration_h", extra)
    # --- descriptive dry-weather CML statistics (the only CML read of this stage) ---
    dstats, read_s = dry_stats_table(cfg, dp)
    events.validate_unique(dstats, ["cml_id", "sublink_id"], "sublink_dry_stats.csv")
    _write_table(cfg, "dry_stats", dstats, "events", "descriptive CML statistics per sublink over the dry-period 5 min "
                 "bins (no baseline model)", {**extra, "core_read_s": round(read_s, 2)})
    write_reference(cfg, ref, info)
    rec = {"step": "events", "status": "written", "n_events": len(et), "n_dropped_short": dropped,
           "n_dry_periods": len(dp), "dry_hours": float(dp.duration_h.sum()) if len(dp) else 0.0,
           "definition_runtime_s": round(t_def, 2), "core_read_s": round(read_s, 2),
           "runtime_s": round(time.perf_counter() - t0, 2)}
    log.info("events: %d events (%d short dropped), %d dry periods (%.0f h), CML core read %.1f s, total %.1f s",
             len(et), dropped, len(dp), rec["dry_hours"], read_s, rec["runtime_s"])
    records.append(rec)


# ---------------------------------------------------------------------------
# Stage: split (frozen)
# ---------------------------------------------------------------------------
def _read_split(cfg: Config, key: str) -> pd.DataFrame | None:
    p = table_path(cfg, key)
    return pd.read_csv(p) if p.exists() else None


def annotate_events(cfg: Config) -> None:
    """Put split_event into events.csv (the split files are the source of truth; nothing else changes)."""
    ep, sp = table_path(cfg, "events"), _read_split(cfg, "split_events")
    if not ep.exists() or sp is None:
        return
    meta = _sidecar_meta(ep)
    cur = pd.read_csv(ep)
    et = cur.drop(columns=["split_event"], errors="ignore")
    et = et.merge(sp[["event_id", "split_event"]], on="event_id", how="left")
    events.validate_events(et, sp)
    if "split_event" in cur and cur.split_event.equals(et.split_event):
        return  # already annotated: do not rewrite (a rewrite would change the bytes and invalidate event_links)
    tmp = ep.with_name(ep.name + ".tmp")
    et.to_csv(tmp, index=False)
    tmp.replace(ep)
    meta["columns"] = list(et.columns)
    Path(str(ep) + ".meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")


def run_split(cfg: Config, force_split: bool, log, records: list[dict]) -> None:
    """Stage split: write the frozen split (refuses to overwrite without --force-split)."""
    t0 = time.perf_counter()
    et = pd.read_csv(table_path(cfg, "events"))
    se, sl = table_path(cfg, "split_events"), table_path(cfg, "split_links")
    if (se.exists() or sl.exists()) and not force_split:
        events.validate_events(et, _read_split(cfg, "split_events"))  # every event must already have a split row
        msg = ("split: tables/split_events.csv and split_links.csv exist and are FROZEN; not overwritten. "
               "Use --force-split only if you really mean to redraw the train/test split.")
        log.info(msg)
        print(msg)
        annotate_events(cfg)
        records.append({"step": "split", "status": "frozen_kept"})
        return
    scfg = cfg["split"]
    sp_e = events.split_events(et, scfg["event_fractions"], scfg["seed"])
    lg = pd.read_csv(cfg.tables_dir / "link_geometry.csv")
    status = pd.read_csv(cfg.tables_dir / "sublink_qc_status.csv")
    sp_l = events.split_links(lg, status, scfg)
    sp_l = sp_l.merge(lg[["cml_id", "length_m", "frequency_ghz_sublink_1", "site_0_id", "site_1_id"]], on="cml_id")
    events.validate_unique(sp_e, ["event_id"], "split_events.csv")
    events.validate_unique(sp_l, ["cml_id"], "split_links.csv")
    _write_table(cfg, "split_events", sp_e, "split", "FROZEN event split (train/val/test), stratified by intensity class")
    _write_table(cfg, "split_links", sp_l, "split", "FROZEN link split (train/test/excluded), stratified by length x "
                 "frequency band; parallel links (same hop) share a split")
    annotate_events(cfg)
    rec = {"step": "split", "status": "written", "events": sp_e.split_event.value_counts().to_dict(),
           "links": sp_l.split_link.value_counts().to_dict(), "runtime_s": round(time.perf_counter() - t0, 2)}
    log.info("split: events %s, links %s", rec["events"], rec["links"])
    records.append(rec)


# ---------------------------------------------------------------------------
# Stage: event_links
# ---------------------------------------------------------------------------
FLAG_BITS = {k: np.uint16(1 << schema.FLAGS[k]) for k in ("rx_floor", "spike", "network_gap")}


def event_links_table(cfg: Config, et: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    """One row per (event, link). Radar from radar_along_links_5min, gauges from the 5 min City / 15 min SMHI
    products, CML quality from the consolidated 5 min core file."""
    lcfg, ecfg = cfg["event_links"], cfg["events"]
    res = ecfg["resolution"]
    step = pd.Timedelta(res)
    dt_h = step / pd.Timedelta("1h")
    lg = pd.read_csv(cfg.tables_dir / "link_geometry.csv").set_index("cml_id")
    t0 = time.perf_counter()
    times, cids, subs, core = read_core(cfg, res, ["tl_mean", "flags_or"])
    read_s = time.perf_counter() - t0
    with xr.open_dataset(aggregate_runner.radar_path(cfg, "along_links", res)) as r:
        r = r.sel(cml_id=cids).load()
    if not np.array_equal(r.time.values, times.values):
        raise ValueError("radar product and CML core have different time axes")
    R = {k: r[k].transpose("time", "cml_id").values.astype("float64")
         for k in ("R_mean", "wet_fraction", "R_max", "valid_weight_fraction")}
    city = read_city(cfg, res)
    smhi15 = read_smhi(cfg)
    status = pd.read_csv(cfg.tables_dir / "sublink_qc_status.csv").set_index(["cml_id", "sublink_id"]).status
    m = pd.read_csv(cfg.tables_dir / "link_gauge_matches.csv")
    near = (m[m.distance_to_path_m <= lcfg["gauge_match_max_m"]].sort_values(["cml_id", "distance_to_path_m"])
            .groupby("cml_id").first())
    n1k = m[m.distance_to_path_m <= lcfg["reference_gauge_max_m"]].groupby("cml_id").size()
    pad_b = int(pd.Timedelta(hours=lcfg["pad_before_h"]) / step)
    pad_a = int(pd.Timedelta(hours=lcfg["pad_after_h"]) / step)
    tl, fl = core["tl_mean"], core["flags_or"]  # (time, sublink, cml)
    idx = pd.Index(times)
    rows = []
    for e in et.itertuples():
        i0 = idx.get_loc(pd.Timestamp(e.start))
        i1 = idx.get_loc(pd.Timestamp(e.end) - step) + 1
        sl, pre, post = slice(i0, i1), slice(max(i0 - pad_b, 0), i0), slice(i1, min(i1 + pad_a, len(idx)))
        vf = {w: (np.isfinite(tl[s]).mean(axis=0) if s.stop > s.start else np.full(tl.shape[1:], np.nan))
              for w, s in (("event", sl), ("before", pre), ("after", post))}  # (sublink, cml)
        fr = {k: ((fl[sl] & b) != 0).mean(axis=0) for k, b in FLAG_BITS.items()}
        Rm = R["R_mean"][sl]
        radar = {"radar_total_mm": np.nansum(Rm, axis=0) * dt_h, "radar_peak_mmh": _nanmax(Rm),
                 "radar_wet_minutes": (np.nan_to_num(Rm) >= ecfg["wet_threshold_mmh"]).sum(axis=0) * step / pd.Timedelta("1min"),
                 "radar_mean_wet_fraction": _nanmean(R["wet_fraction"][sl]),
                 "radar_max_R_max_mmh": _nanmax(R["R_max"][sl]),
                 "radar_mean_valid_weight_fraction": _nanmean(R["valid_weight_fraction"][sl])}
        for j, c in enumerate(cids):
            g = lg.loc[c]
            row = {"event_id": e.event_id, "cml_id": int(c), "length_m": g.length_m,
                   "frequency_ghz_sublink_1": g.frequency_ghz_sublink_1, "frequency_ghz_sublink_2": g.frequency_ghz_sublink_2,
                   "polarization": g.polarization, "n_cells": int(g.n_radar_cells),
                   **{k: float(v[j]) for k, v in radar.items()}}
            if c in near.index:
                gm = near.loc[c]
                if gm.network == "city":
                    a = city.amount_mm.sel(id=gm.gauge_id).values[sl].astype("float64")
                    rr = city.rate_mm_h.sel(id=gm.gauge_id).values[sl].astype("float64")
                    gt = float(np.nansum(a)) if np.isfinite(a).any() else np.nan
                    gr, gi = float(_nanmax(rr)), "5min"
                else:  # SMHI: 15 min native, never disaggregated: overlapping 15 min bins, max 15 min rate
                    gt = _window_sum_15(smhi15, e.start, e.end)
                    w = smhi15[(smhi15.index >= pd.Timestamp(e.start).floor("15min")) & (smhi15.index < pd.Timestamp(e.end))]
                    gr, gi = (float(w.max() / 0.25) if len(w) else np.nan), "15min"
                row.update({"gauge_id": gm.gauge_id, "gauge_network": gm.network, "gauge_distance_to_path_m": gm.distance_to_path_m,
                            "gauge_position_along_path": gm.position_along_path, "gauge_total_mm": gt,
                            "gauge_max_rate_mmh": gr, "gauge_max_rate_interval": gi})
            else:
                row.update({"gauge_id": "", "gauge_network": "", "gauge_distance_to_path_m": np.nan,
                            "gauge_position_along_path": np.nan, "gauge_total_mm": np.nan, "gauge_max_rate_mmh": np.nan,
                            "gauge_max_rate_interval": ""})
            row["n_gauges_within_1km"] = int(n1k.get(c, 0))
            for w in ("event", "before", "after"):
                col = vf[w][:, j]
                for k, s in enumerate(subs):
                    row[f"valid_fraction_{w}_{s}"] = float(col[k])
                row[f"valid_fraction_{w}_min"] = float(np.min(col)) if np.isfinite(col).all() else np.nan
            for f, a in fr.items():
                for k, s in enumerate(subs):
                    row[f"frac_{f}_event_{s}"] = float(a[k, j])
                row[f"frac_{f}_event_max"] = float(a[:, j].max())
            for s in subs:
                row[f"status_{s}"] = status.get((c, s), "")
            # cml_ok_<window>: at least ONE sub-link has >= min_valid_fraction valid tl_mean bins in the window
            for w in ("event", "before", "after"):
                col = vf[w][:, j]
                row[f"cml_ok_{w}"] = bool(np.isfinite(col).any() and np.nanmax(col) >= lcfg["min_valid_fraction"])
            row["has_reference_gauge"] = bool(row["n_gauges_within_1km"] > 0)
            rows.append(row)
    return pd.DataFrame(rows), read_s


def run_event_links(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """Stage event_links: one row per (event, link)."""
    t0 = time.perf_counter()
    et = pd.read_csv(table_path(cfg, "events"))
    se, sl = _read_split(cfg, "split_events"), _read_split(cfg, "split_links")
    if se is None or sl is None:
        raise FileNotFoundError("run the split stage before event_links")
    events.validate_events(et, se)
    h = hashlib.sha256()
    for p in (table_path(cfg, "events"), table_path(cfg, "split_events"), table_path(cfg, "split_links"),
              cfg.tables_dir / "link_gauge_matches.csv", cfg.tables_dir / "sublink_qc_status.csv"):
        h.update(p.read_bytes())
    h.update(str(aggregate_runner.file_attrs(core_path(cfg, cfg["events"]["resolution"])).get("input_signature")).encode())
    h.update(str(aggregate_runner.file_attrs(aggregate_runner.radar_path(cfg, "along_links", "5min")).get("created_utc")).encode())
    sig = h.hexdigest()
    out = table_path(cfg, "event_links")
    meta = _sidecar_meta(out)
    if not force and meta.get("stage_hash") == cfg.stage_hash("event_links") and meta.get("input_signature") == sig:
        log.info("event_links: up to date, skipped (use --force to rebuild)")
        records.append({"step": "event_links", "status": "skipped_same_config"})
        return
    df, read_s = event_links_table(cfg, et)
    df = (df.merge(se[["event_id", "split_event"]], on="event_id", how="left")
            .merge(sl[["cml_id", "split_link"]], on="cml_id", how="left"))
    events.validate_unique(df, ["event_id", "cml_id"], "event_links.csv")
    if df.split_event.isna().any() or df.split_link.isna().any():
        raise events.EventsValidationError("event_links.csv: rows without a split")
    _write_table(cfg, "event_links", df, "event_links", "one row per (event, link): radar along the path, matched gauge, "
                 "CML data quality (CSV because pyarrow is not installed)", {"input_signature": sig, "core_read_s": round(read_s, 2)})
    rec = {"step": "event_links", "status": "written", "n_rows": len(df), "core_read_s": round(read_s, 2),
           "file_size_bytes": out.stat().st_size, "runtime_s": round(time.perf_counter() - t0, 2)}
    log.info("event_links: %d rows (%.1f MB), %.1f s", len(df), rec["file_size_bytes"] / 1e6, rec["runtime_s"])
    records.append(rec)
