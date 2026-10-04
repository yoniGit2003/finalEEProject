"""QC stages for the CML data: period statistics (qc_stats), daily flags (qc_flags), sub-link status (qc_status).

The flag logic itself is in qc.py (pure functions). This module reads L1 files, calls those
functions, and writes:
  tables/sublink_period_stats.csv (+ .meta.json, + sublink_period_hist.npz)   qc_stats
  L2_qc/cml_flags_10s_YYYY-MM-DD.nc                                        qc_flags
  tables/qc_daily_summary.csv, tables/network_gaps.csv                      qc_flags
  tables/sublink_qc_status.csv                                              qc_status
L1 files are only read, never modified.

Internally a day's signals are 2-D arrays (time, column) with one column per (sublink_id, cml_id)
pair, column = s * n_cml + c, i.e. the L1 array (time, sublink_id, cml_id) reshaped in C order.
"""

from __future__ import annotations

import json
import logging
import os
import time
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr

from cmlpipe import qc
from cmlpipe.io import schema
from cmlpipe.config import Config
from cmlpipe.io.netcdf import existing_stage_hash
from cmlpipe.parallel import config_from_parts, run_day_task
from cmlpipe.provenance import display_path
from cmlpipe.provenance import provenance_attrs, write_sidecar

log = logging.getLogger(__name__)

# Histogram grids for qc_stats. L1 values are exact multiples of 0.1 dB, so 0.1 dB bins give exact percentiles.
RSL_BINS = (-150.0, 50.0)
TSL_BINS = (-50.0, 50.0)
DRSL_MAX = 100.0
BIN = 0.1


# ---------------------------------------------------------------------------
# Paths and reading
# ---------------------------------------------------------------------------
def l2_dir(cfg: Config) -> Path:
    """Folder of the L2 flag files."""
    return cfg.output_root / cfg["paths"].get("l2_qc_dir", "L2_qc")


def l2_day_path(cfg: Config, day: pd.Timestamp) -> Path:
    """L2 flag file of one day."""
    return l2_dir(cfg) / f"cml_flags_{cfg['cml']['sampling_interval']}_{pd.Timestamp(day):%Y-%m-%d}.nc"


def read_l1(cfg: Config, day: pd.Timestamp, hours: float | None = None, head: bool = True) -> xr.Dataset | None:
    """rsl/tsl of one L1 day (or only its first/last `hours`), or None if the file does not exist."""
    path = cfg.l1_cml_day_path(day)
    if not path.exists():
        return None
    with xr.open_dataset(path) as ds:
        if hours is not None:
            n = int(round(pd.Timedelta(hours=hours) / cfg.sampling_interval))
            ds = ds.isel(time=slice(0, n) if head else slice(-n, None)) if n > 0 else ds.isel(time=slice(0, 0))
        return ds[["rsl", "tsl"]].load()


def to_2d(a: np.ndarray) -> np.ndarray:
    """(time, sublink_id, cml_id) -> (time, column)."""
    return a.reshape(a.shape[0], -1)


def column_table(ds: xr.Dataset) -> pd.DataFrame:
    """One row per column, in column order: cml_id, sublink_id, raw_sublink_id, length, frequency."""
    s_ids, c_ids = ds.sublink_id.values, ds.cml_id.values
    return pd.DataFrame({
        "cml_id": np.tile(c_ids, len(s_ids)),
        "sublink_id": np.repeat(s_ids, len(c_ids)),
        "raw_sublink_id": ds.raw_sublink_id.values.ravel(),
        "length_km": np.tile(ds.length.values, len(s_ids)),
        "frequency_ghz": ds.frequency.values.ravel(),
    })


# ---------------------------------------------------------------------------
# qc_stats: pass 1 over the whole period
# ---------------------------------------------------------------------------
def _hist_index(v: np.ndarray, lo: float, hi: float) -> tuple[np.ndarray, np.ndarray, int]:
    n_bins = int(round((hi - lo) / BIN)) + 1
    ok = np.isfinite(v)
    idx = np.clip(np.rint((np.where(ok, v, lo) - lo) / BIN), 0, n_bins - 1).astype("int64")
    return idx, ok, n_bins


def _accumulate_hist(hist: np.ndarray, v: np.ndarray, lo: float, hi: float) -> None:
    """hist[column, bin] += count of values of each column in each 0.1 dB bin (vectorized with bincount)."""
    idx, ok, n_bins = _hist_index(v, lo, hi)
    cols = np.broadcast_to(np.arange(v.shape[1]), v.shape)
    flat = (cols[ok] * n_bins + idx[ok])
    hist += np.bincount(flat, minlength=hist.size).reshape(hist.shape)


def _hist_quantile(hist: np.ndarray, q: float, lo: float) -> np.ndarray:
    """Value at quantile q (0..1) per row of a histogram (lower value: smallest v with CDF >= q)."""
    n = hist.sum(axis=1)
    cdf = np.cumsum(hist, axis=1)
    target = np.maximum(np.ceil(q * n), 1)[:, None]
    idx = (cdf >= target).argmax(axis=1)
    return np.where(n > 0, lo + idx * BIN, np.nan)


def floor_level(hist: np.ndarray, lo: float, min_samples: int, window_db: float) -> np.ndarray:
    """Robust floor level per row of an rsl histogram (0.1 dB bins starting at `lo`).

    The lowest OCCUPIED value L such that at least min_samples samples lie in [L, L + window_db].
    One isolated low sample cannot set it; a receiver that sits at its floor for minutes does.
    NaN for rows with too few samples anywhere.
    """
    w = int(round(window_db / BIN))
    csum = np.cumsum(np.pad(hist, ((0, 0), (1, 0))), axis=1)
    n_bins = hist.shape[1]
    upper = np.minimum(np.arange(n_bins) + w + 1, n_bins)
    in_window = csum[:, upper] - csum[:, :n_bins]
    ok = (in_window >= min_samples) & (hist > 0)
    idx = ok.argmax(axis=1)
    return np.where(ok.any(axis=1), lo + idx * BIN, np.nan)


def period_stats(cfg: Config, days: list[pd.Timestamp]) -> tuple[pd.DataFrame, dict]:
    """Statistics per sub-link over the given days (qc_stats, pass 1).

    Reads each L1 day once. Carries state across midnight (last valid tsl, last valid rsl and the
    start of the current constant-rsl run) so steps and runs that cross days are counted correctly.
    Returns (table, histograms) where histograms = {"rsl", "tsl", "drsl"} arrays (column, bin).
    """
    interval = cfg.sampling_interval
    info = None
    h_rsl = h_tsl = h_drsl = None
    n_samples = n_rsl = n_tsl = n_tl = tsl_steps = None
    last_rsl = last_tsl = run_value = run_start = longest = None
    n_valid_before = None  # valid rsl samples counted so far, per column (for run lengths in valid samples)
    for day in days:
        ds = xr.open_dataset(cfg.l1_cml_day_path(day)).load()
        rsl, tsl = to_2d(ds.rsl.values.astype("float64")), to_2d(ds.tsl.values.astype("float64"))
        if info is None:
            info = column_table(ds)
            n_col = rsl.shape[1]
            h_rsl = np.zeros((n_col, int(round((RSL_BINS[1] - RSL_BINS[0]) / BIN)) + 1), "int64")
            h_tsl = np.zeros((n_col, int(round((TSL_BINS[1] - TSL_BINS[0]) / BIN)) + 1), "int64")
            h_drsl = np.zeros((n_col, int(round(DRSL_MAX / BIN)) + 1), "int64")
            n_samples = np.zeros(n_col, "int64")
            n_rsl, n_tsl, n_tl, tsl_steps = (np.zeros(n_col, "int64") for _ in range(4))
            last_rsl = last_tsl = run_value = np.full(n_col, np.nan)
            run_start = np.full(n_col, np.nan)
            longest = np.zeros(n_col, "int64")
            n_valid_before = np.zeros(n_col, "float64")
        ds.close()

        n_samples += rsl.shape[0]
        n_rsl += np.isfinite(rsl).sum(axis=0)
        n_tsl += np.isfinite(tsl).sum(axis=0)
        n_tl += (np.isfinite(rsl) & np.isfinite(tsl)).sum(axis=0)
        _accumulate_hist(h_rsl, rsl, *RSL_BINS)
        _accumulate_hist(h_tsl, tsl, *TSL_BINS)

        # |delta rsl| between ADJACENT valid samples (10 s apart), including across midnight
        prev = np.vstack([last_rsl[None, :], rsl[:-1]])
        _accumulate_hist(h_drsl, np.abs(rsl - prev), 0.0, DRSL_MAX)
        last_rsl = rsl[-1].copy()

        # tsl steps: versus the previous VALID tsl sample (carried across days)
        ff_tsl = pd.DataFrame(np.vstack([last_tsl[None, :], tsl])).ffill().to_numpy()
        prev_t = ff_tsl[:-1]
        tsl_steps += (np.isfinite(tsl) & np.isfinite(prev_t) & (np.abs(tsl - prev_t) >= 1.0 - 1e-6)).sum(axis=0)
        last_tsl = ff_tsl[-1]

        # longest exactly-constant rsl run, counted in VALID samples (NaN neither breaks nor extends
        # a run), as the flatline flag; carried across days
        valid = np.isfinite(rsl)
        g = n_valid_before[None, :] + np.cumsum(valid, axis=0)
        ff = pd.DataFrame(np.vstack([run_value[None, :], rsl])).ffill().to_numpy()
        prev_valid = ff[:-1]
        start = valid & ~(prev_valid == rsl)
        seg_start = pd.DataFrame(np.vstack([run_start[None, :], np.where(start, g, np.nan)])).ffill().to_numpy()[1:]
        span = np.where(valid, g - seg_start + 1, 0)
        longest = np.maximum(longest, np.nan_to_num(span, nan=0).max(axis=0).astype("int64"))
        run_value = ff[-1]
        run_start = seg_start[-1]
        n_valid_before = g[-1].astype("float64")
        log.info("qc_stats: %s read", day.date())

    per_min = pd.Timedelta("1min") / interval
    table = info.copy()
    table["n_samples"] = n_samples
    table["n_valid_rsl"], table["n_valid_tsl"], table["n_valid_tl"] = n_rsl, n_tsl, n_tl
    with np.errstate(invalid="ignore", divide="ignore"):
        table["completeness_rsl_pct"] = np.round(100.0 * n_rsl / n_samples, 4)
        table["completeness_tl_pct"] = np.round(100.0 * n_tl / n_samples, 4)
    lo = RSL_BINS[0]
    for name, q in (("rsl_min", 0.0), ("rsl_p0_1", 0.001), ("rsl_p1", 0.01), ("rsl_median", 0.5),
                    ("rsl_p99", 0.99)):
        table[name] = np.round(_hist_quantile(h_rsl, q, lo), 1)
    has = h_rsl.sum(axis=1) > 0
    top = h_rsl.shape[1] - 1 - np.argmax(h_rsl[:, ::-1] > 0, axis=1)
    table["rsl_max"] = np.where(has, np.round(lo + top * BIN, 1), np.nan)
    table["tsl_min"] = np.round(_hist_quantile(h_tsl, 0.0, TSL_BINS[0]), 1)
    table["tsl_median"] = np.round(_hist_quantile(h_tsl, 0.5, TSL_BINS[0]), 1)
    has_t = h_tsl.sum(axis=1) > 0
    top_t = h_tsl.shape[1] - 1 - np.argmax(h_tsl[:, ::-1] > 0, axis=1)
    table["tsl_max"] = np.where(has_t, np.round(TSL_BINS[0] + top_t * BIN, 1), np.nan)
    table["tsl_n_distinct"] = (h_tsl > 0).sum(axis=1)
    table["tsl_n_steps"] = tsl_steps
    nz = h_drsl.copy()
    nz[:, 0] = 0  # non-zero deltas only
    order = np.argsort(-nz, axis=1, kind="stable")[:, :3]
    for k in range(3):
        cnt = np.take_along_axis(nz, order[:, k:k + 1], axis=1)[:, 0]
        table[f"drsl_mode{k + 1}"] = np.where(cnt > 0, np.round(order[:, k] * BIN, 1), np.nan)
    table["rsl_longest_constant_min"] = np.round(longest / per_min, 2)
    table["floor_level_dbm"] = np.round(floor_level(h_rsl, lo, int(cfg["qc"]["stats"]["floor_min_samples"]),
                                                    float(cfg["qc"]["stats"]["floor_window_db"])), 1)
    return table, {"rsl": h_rsl, "tsl": h_tsl, "drsl": h_drsl}


def write_period_stats(cfg: Config, table: pd.DataFrame, hists: dict, days: list[pd.Timestamp],
                       l1_hashes: set) -> Path:
    """Write sublink_period_stats.csv, its histograms and sidecar."""
    cfg.tables_dir.mkdir(parents=True, exist_ok=True)
    path = cfg.tables_dir / "sublink_period_stats.csv"
    table.to_csv(path, index=False)
    np.savez_compressed(cfg.tables_dir / "sublink_period_hist.npz", rsl=hists["rsl"], tsl=hists["tsl"],
                        drsl=hists["drsl"], rsl_lo=RSL_BINS[0], tsl_lo=TSL_BINS[0], bin=BIN,
                        cml_id=table["cml_id"].to_numpy(), sublink_id=table["sublink_id"].to_numpy().astype(str))
    write_sidecar(path, cfg, "qc_stats", {
        "description": "per sub-link statistics of the L1 CML data over the period (qc_stats, pass 1)",
        "days": [f"{d:%Y-%m-%d}" for d in days],
        "input_stage_hash": sorted(l1_hashes),
        "floor_level": (f"lowest occupied rsl value L with >= {cfg['qc']['stats']['floor_min_samples']} samples in "
                        f"[L, L + {cfg['qc']['stats']['floor_window_db']} dB] (qc.stats)"),
        "histograms": "sublink_period_hist.npz: 0.1 dB bins; rsl from rsl_lo, tsl from tsl_lo, |delta rsl| from 0",
    })
    return path


def read_period_stats(cfg: Config) -> tuple[pd.DataFrame, dict]:
    """The qc_stats table and its sidecar (raises with a clear message if qc_stats has not run)."""
    path = cfg.tables_dir / "sublink_period_stats.csv"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found: run the qc_stats step first")
    meta = json.loads(Path(str(path) + ".meta.json").read_text(encoding="utf-8"))
    return pd.read_csv(path), meta


# ---------------------------------------------------------------------------
# qc_flags: one day
# ---------------------------------------------------------------------------
def load_with_context(cfg: Config, day: pd.Timestamp) -> tuple[xr.Dataset, np.ndarray, np.ndarray, slice]:
    """The L1 day plus `qc.context_hours` of the neighbouring L1 days (if those files exist).

    Returns (the L1 day dataset, rsl block, tsl block (time, column), slice of the day in the block).
    At the edges of the period the missing context is simply absent.
    """
    hours = float(cfg["qc"]["context_hours"])
    day = pd.Timestamp(day).normalize()
    with xr.open_dataset(cfg.l1_cml_day_path(day)) as ds:
        ds = ds.load()
    parts = []
    before = read_l1(cfg, day - pd.Timedelta("1D"), hours, head=False) if hours > 0 else None
    after = read_l1(cfg, day + pd.Timedelta("1D"), hours, head=True) if hours > 0 else None
    for part in (before, ds[["rsl", "tsl"]], after):
        if part is not None and part.sizes["time"] > 0:
            parts.append(part)
    block = xr.concat(parts, dim="time") if len(parts) > 1 else parts[0]
    t = pd.DatetimeIndex(block.time.values)
    if not (np.diff(t.asi8) == cfg.sampling_interval.value).all():
        raise ValueError(f"{day.date()}: time axis with context is not regular")
    n_before = before.sizes["time"] if before is not None else 0
    day_slice = slice(n_before, n_before + ds.sizes["time"])
    return (ds, to_2d(block.rsl.values.astype("float64")), to_2d(block.tsl.values.astype("float64")), day_slice)


def compute_flags(cfg: Config, rsl: np.ndarray, tsl: np.ndarray, day_slice: slice, floor: np.ndarray,
                  median: np.ndarray) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
    """All flags for one day. Returns ({flag: bool mask on the day}, network NaN fraction, alive columns)."""
    f = cfg["qc"]["flags"]
    interval = cfg.sampling_interval
    rsl_day, tsl_day = rsl[day_slice], tsl[day_slice]
    alive = qc.alive_columns(rsl_day, float(cfg["qc"]["stats"]["alive_min_valid_fraction"]))
    frac = qc.network_nan_fraction(rsl_day, alive)
    masks = {
        "missing_rsl": qc.flag_missing(rsl_day),
        "missing_tsl": qc.flag_missing(tsl_day),
        "out_of_range": qc.flag_out_of_range(rsl_day, tsl_day, f["valid_range"]["rsl_dbm"], f["valid_range"]["tsl_dbm"]),
        # the run-based flags see the context, then are trimmed to the day
        "rx_floor": qc.flag_rx_floor(rsl, floor, median, float(f["rx_floor"]["margin_db"]), float(f["rx_floor"]["depth_db"]),
                                     qc.n_samples(f["rx_floor"]["min_duration"], interval),
                                     float(f["rx_floor"]["receiver_floor_max_dbm"]))[day_slice],
        "spike": qc.flag_spike(rsl, qc.n_samples(f["spike"]["window"], interval) + 1, float(f["spike"]["n_mad"]),
                               float(f["spike"]["min_abs_db"]), int(f["spike"]["max_len"]),
                               float(f["spike"]["neighbour_agree_db"]),
                               float(f["spike"]["neighbour_max_dev_db"]))[day_slice],
        "tsl_step": qc.flag_tsl_step(tsl, float(f["tsl_step"]["min_db"]))[day_slice],
        "flatline": qc.flag_flatline(rsl, qc.n_samples(f["flatline"]["min_duration"], interval))[day_slice],
        "network_gap": qc.flag_network_gap(rsl_day, frac, float(f["network_gap"]["fraction"]), alive),
    }
    return masks, frac, alive


def pair_correlation(rsl3: np.ndarray, window: int, min_common: int) -> np.ndarray:
    """Per link: correlation of the two sub-links' rsl after subtracting a rolling median (NaN if too few samples)."""
    out = np.full(rsl3.shape[2], np.nan)
    det = []
    for s in range(rsl3.shape[1]):
        df = pd.DataFrame(rsl3[:, s, :])
        det.append((df - df.rolling(window, center=True, min_periods=window // 2 + 1).median()).to_numpy())
    a, b = det[0], det[1]
    ok = np.isfinite(a) & np.isfinite(b)
    n = ok.sum(axis=0)
    a0, b0 = np.where(ok, a, 0.0), np.where(ok, b, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        ma, mb = a0.sum(0) / n, b0.sum(0) / n
        cov = (np.where(ok, (a - ma) * (b - mb), 0.0)).sum(0)
        va = (np.where(ok, (a - ma) ** 2, 0.0)).sum(0)
        vb = (np.where(ok, (b - mb) ** 2, 0.0)).sum(0)
        r = cov / np.sqrt(va * vb)
    out[:] = np.where((n >= min_common) & (va > 0) & (vb > 0), r, np.nan)
    return out


def flag_day(cfg: Config, day, stats: pd.DataFrame) -> tuple[xr.Dataset, pd.DataFrame]:
    """Flags of one UTC day as an L2 dataset, plus the day's summary rows (one per sub-link)."""
    day = pd.Timestamp(day).normalize()
    ds, rsl, tsl, day_slice = load_with_context(cfg, day)
    cols = column_table(ds)
    s = stats.set_index(["cml_id", "sublink_id"]).reindex(pd.MultiIndex.from_frame(cols[["cml_id", "sublink_id"]]))
    masks, frac, alive = compute_flags(cfg, rsl, tsl, day_slice, s["floor_level_dbm"].to_numpy(float),
                                       s["rsl_median"].to_numpy(float))
    flags2d = qc.combine(masks)
    shape3 = ds.rsl.shape
    flags = flags2d.reshape(shape3)

    out = xr.Dataset(
        {schema.FLAGS_VAR: (schema.DIMS, flags),
         "network_nan_fraction": (("time",), frac.astype("float32"))},
        coords={"time": ds.time.values, "sublink_id": ds.sublink_id.values, "cml_id": ds.cml_id.values,
                "raw_sublink_id": (("sublink_id", "cml_id"), ds.raw_sublink_id.values)},
    )
    out[schema.FLAGS_VAR].attrs = schema.flag_attrs()
    out["network_nan_fraction"].attrs = {
        "long_name": "fraction of alive sub-links whose rsl is NaN", "units": "1",
        "comment": f"alive = >= {cfg['qc']['stats']['alive_min_valid_fraction']} valid rsl samples on this day"}
    out["raw_sublink_id"].attrs = {"long_name": "raw 'Sublink' number (-1 = absent)"}
    counts = {name: int(m.sum()) for name, m in masks.items()}
    out.attrs = {
        "title": f"{cfg['dataset']['title']} - L2 CML QC flags, {day:%Y-%m-%d}",
        "processing_level": "L2_qc",
        "day_utc": f"{day:%Y-%m-%d}",
        "time_meaning": ds.attrs.get("time_meaning", ""),
        "input_file": cfg.l1_cml_day_path(day).relative_to(cfg.repo_root).as_posix()
        if cfg.l1_cml_day_path(day).is_relative_to(cfg.repo_root) else str(cfg.l1_cml_day_path(day)),
        "input_stage_hash": ds.attrs.get("stage_hash", ""),
        "context_hours": float(cfg["qc"]["context_hours"]),
        "context_samples_before": int(day_slice.start),
        "context_samples_after": int(rsl.shape[0] - day_slice.stop),
        "blocking_flags": " ".join(cfg["qc"]["blocking_flags"]),
        "qc_config": json.dumps(cfg["qc"], sort_keys=True),
        "n_alive_sublinks": int(alive.sum()),
        "principle": "flag, don't fix: L1 values are unchanged; aggregation drops blocking flags later",
        **{f"n_{k}": v for k, v in counts.items()},
    }

    # summary rows
    min_common = qc.n_samples("1h", cfg.sampling_interval)
    win = qc.n_samples(cfg["qc"]["summary"]["pair_corr_detrend_window"], cfg.sampling_interval) + 1
    rsl3 = ds.rsl.values.astype("float64")
    corr = pair_correlation(rsl3, win, min_common)
    blocking = qc.blocking_mask(flags2d, cfg["qc"]["blocking_flags"])
    rsl_day, tsl_day = rsl[day_slice], tsl[day_slice]
    has = np.isfinite(rsl_day).any(axis=0)
    med = np.full(rsl_day.shape[1], np.nan)
    if has.any():
        med[has] = np.nanmedian(rsl_day[:, has], axis=0)  # only columns with data: no all-NaN warning
    summary = cols[["cml_id", "sublink_id", "raw_sublink_id"]].copy()
    summary.insert(0, "day", f"{day:%Y-%m-%d}")
    summary["n_samples"] = rsl_day.shape[0]
    summary["n_valid_rsl"] = np.isfinite(rsl_day).sum(axis=0)
    summary["n_valid_tl"] = (np.isfinite(rsl_day) & np.isfinite(tsl_day)).sum(axis=0)
    for name, m in masks.items():
        summary[f"n_{name}"] = m.sum(axis=0)
    summary["n_blocking"] = blocking.sum(axis=0)
    summary["rsl_median"] = np.round(med, 2)
    summary["pair_corr"] = np.round(np.tile(corr, rsl3.shape[1]), 4)
    return out, summary


def existing_l2_hashes(path: Path) -> tuple[str | None, str | None, str | None]:
    """(stage_hash, input_stage_hash, stats_stage_hash) of an existing L2 file, or Nones."""
    if not path.exists():
        return None, None, None
    try:
        with netCDF4.Dataset(path) as nc:
            a = {k: nc.getncattr(k) for k in nc.ncattrs()}
        return a.get("stage_hash"), a.get("input_stage_hash"), a.get("stats_stage_hash")
    except OSError:
        return None, None, None


def write_l2(ds: xr.Dataset, path: Path, cfg: Config, stats_stage_hash: str) -> dict:
    """Write an L2 flags file (uint16, zlib 4) via .tmp, re-open and check it is identical, then rename."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.unlink(missing_ok=True)
    ds = ds.copy()
    ds.attrs.update(provenance_attrs(cfg, "qc_flags"))
    ds.attrs["stats_stage_hash"] = stats_stage_hash
    enc = {
        "time": {"units": "seconds since 1970-01-01 00:00:00", "calendar": "proleptic_gregorian", "dtype": "int64"},
        schema.FLAGS_VAR: {"dtype": "uint16", "zlib": True, "complevel": 4, "_FillValue": None},
        "network_nan_fraction": {"dtype": "float32", "zlib": True, "complevel": 4, "_FillValue": np.float32(np.nan)},
    }
    ds.to_netcdf(tmp, format="NETCDF4", engine="netcdf4", encoding=enc)
    try:
        with xr.open_dataset(tmp) as back:
            back.load()
            if not (np.array_equal(back[schema.FLAGS_VAR].values, ds[schema.FLAGS_VAR].values)
                    and back[schema.FLAGS_VAR].dtype == np.uint16
                    and np.array_equal(back.time.values, ds.time.values)
                    and np.array_equal(back.network_nan_fraction.values, ds.network_nan_fraction.values, equal_nan=True)):
                raise ValueError(f"round trip of {path.name} failed")
    except Exception:
        tmp.unlink(missing_ok=True)
        raise
    os.replace(tmp, path)
    return {"path": path, "file_size_bytes": path.stat().st_size}


# ---------------------------------------------------------------------------
# Tables built from the L2 files
# ---------------------------------------------------------------------------
def upsert_daily_summary(cfg: Config, rows: pd.DataFrame) -> Path:
    """Replace the rows of the given days in tables/qc_daily_summary.csv (other days are kept)."""
    path = cfg.tables_dir / "qc_daily_summary.csv"
    cfg.tables_dir.mkdir(parents=True, exist_ok=True)
    if path.exists():
        old = pd.read_csv(path)
        rows = pd.concat([old[~old["day"].isin(rows["day"].unique())], rows], ignore_index=True)
    rows = rows.sort_values(["day", "cml_id", "sublink_id"]).reset_index(drop=True)
    rows.to_csv(path, index=False)
    write_sidecar(path, cfg, "qc_flags", {
        "description": "one row per day x sub-link: valid samples, count per flag, blocking count, rsl daily median, "
                       "pair_corr (correlation of the two sub-links' rsl after a rolling-median detrend)",
        "days": sorted(rows["day"].unique().tolist()),
    })
    return path


def network_gaps_table(cfg: Config) -> pd.DataFrame:
    """Periods where the network NaN fraction exceeds qc.flags.network_gap.fraction, from all L2 files."""
    thr = float(cfg["qc"]["flags"]["network_gap"]["fraction"])
    frames = []
    for path in sorted(l2_dir(cfg).glob("cml_flags_*.nc")):
        with xr.open_dataset(path) as ds:
            frames.append(pd.DataFrame({"frac": ds.network_nan_fraction.values.astype(float),
                                        "n_alive": int(ds.attrs.get("n_alive_sublinks", -1))},
                                       index=pd.DatetimeIndex(ds.time.values)))
    cols = ["start", "end", "duration_min", "max_fraction_nan", "n_alive"]
    if not frames:
        return pd.DataFrame(columns=cols)
    s = pd.concat(frames).sort_index()
    step = cfg.sampling_interval
    gap = s["frac"].to_numpy() > thr
    t = s.index
    # a new period starts at a gap sample whose previous sample is not a gap or not the previous time step
    new = gap & ~(np.r_[False, gap[:-1]] & np.r_[False, np.diff(t.asi8) == step.value])
    gid = np.cumsum(new)
    rows = []
    for g in np.unique(gid[gap]):
        sel = gap & (gid == g)
        tt = t[sel]
        rows.append({"start": tt[0], "end": tt[-1] + step, "duration_min": round(sel.sum() * step / pd.Timedelta("1min"), 2),
                     "max_fraction_nan": round(float(s["frac"].to_numpy()[sel].max()), 4),
                     "n_alive": int(s["n_alive"].to_numpy()[sel][0])})
    return pd.DataFrame(rows, columns=cols)


def write_network_gaps(cfg: Config) -> Path:
    """Write network_gaps.csv from the L2 files."""
    path = cfg.tables_dir / "network_gaps.csv"
    table = network_gaps_table(cfg)
    table.to_csv(path, index=False)
    write_sidecar(path, cfg, "qc_flags", {
        "description": "periods where more than qc.flags.network_gap.fraction of the alive sub-links are NaN; "
                       "end is exclusive (last gap sample + one step); built from all L2 files present",
        "n_l2_files": len(list(l2_dir(cfg).glob("cml_flags_*.nc"))),
    })
    return path


# ---------------------------------------------------------------------------
# qc_status
# ---------------------------------------------------------------------------
def sublink_status(cfg: Config, summary: pd.DataFrame) -> pd.DataFrame:
    """Per sub-link status over all days in the daily summary: good / degraded / unusable + reasons.

    flagged_pct = share of VALID rsl samples carrying a blocking flag other than the missing flags
    (out_of_range, rx_floor, spike): missing data is covered by the completeness instead.
    """
    st = cfg["qc"]["status"]
    g = summary.groupby(["cml_id", "sublink_id"], sort=True)
    agg = g[[c for c in summary.columns if c.startswith("n_")]].sum()
    agg["raw_sublink_id"] = g["raw_sublink_id"].first()
    agg["n_days"] = g["day"].nunique()
    out = agg[["raw_sublink_id", "n_days"]].copy()
    out["has_rsl"] = agg["n_valid_rsl"] > 0
    out["has_tl"] = agg["n_valid_tl"] > 0
    with np.errstate(invalid="ignore", divide="ignore"):
        out["completeness_rsl_pct"] = (100.0 * agg["n_valid_rsl"] / agg["n_samples"]).round(3)
        out["completeness_tl_pct"] = (100.0 * agg["n_valid_tl"] / agg["n_samples"]).round(3)
        for name in schema.FLAGS:
            out[f"pct_{name}"] = (100.0 * agg[f"n_{name}"] / agg["n_samples"]).round(4)
        non_missing = [n for n in cfg["qc"]["blocking_flags"] if n not in ("missing_rsl", "missing_tsl", "network_gap")]
        flagged = sum(agg[f"n_{n}"] for n in non_missing)
        out["flagged_pct"] = (100.0 * flagged / agg["n_valid_rsl"]).where(agg["n_valid_rsl"] > 0).round(4)
        rx = (100.0 * agg["n_rx_floor"] / agg["n_valid_rsl"]).where(agg["n_valid_rsl"] > 0)

    status, reason = [], []
    for i in range(len(out)):
        r_unusable, r_degraded = [], []
        if not out["has_rsl"].iloc[i]:
            r_unusable.append("no_rsl")
        elif not out["has_tl"].iloc[i]:
            r_unusable.append("no_tl")
        if out["has_rsl"].iloc[i]:
            if out["completeness_rsl_pct"].iloc[i] < st["unusable_completeness_pct"]:
                r_unusable.append(f"completeness_lt_{st['unusable_completeness_pct']:g}pct")
            if rx.iloc[i] > st["unusable_rx_floor_pct"]:
                r_unusable.append(f"rx_floor_gt_{st['unusable_rx_floor_pct']:g}pct")
            if out["flagged_pct"].iloc[i] > st["unusable_flagged_pct"]:
                r_unusable.append(f"flagged_gt_{st['unusable_flagged_pct']:g}pct")
            elif out["flagged_pct"].iloc[i] > st["degraded_flagged_pct"]:
                r_degraded.append(f"flagged_gt_{st['degraded_flagged_pct']:g}pct")
        if r_unusable:
            status.append("unusable"); reason.append(";".join(r_unusable + r_degraded))
        elif r_degraded:
            status.append("degraded"); reason.append(";".join(r_degraded))
        else:
            status.append("good"); reason.append("ok")
    out["status"], out["reason"] = status, reason
    return out.reset_index()


def write_status(cfg: Config, table: pd.DataFrame, days: list[str]) -> Path:
    """Write sublink_qc_status.csv and its sidecar."""
    path = cfg.tables_dir / "sublink_qc_status.csv"
    table.to_csv(path, index=False)
    write_sidecar(path, cfg, "qc_status", {
        "description": "per sub-link QC status over the period: good / degraded / unusable with machine-readable reasons",
        "days": days, "thresholds": cfg["qc"]["status"], "blocking_flags": cfg["qc"]["blocking_flags"],
    })
    return path


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------
def qc_flags_day(cfg_data: dict, repo_root: str, config_file: str, day_text: str, force: bool = False) -> dict:
    """L1 day -> L2 flags day. The daily-summary rows travel back in record['_summary'] (not in the manifest)."""
    t0 = time.perf_counter()
    cfg = config_from_parts(cfg_data, repo_root, config_file)
    day = pd.Timestamp(day_text)
    stats, stats_meta = read_period_stats(cfg)
    stats_hash = stats_meta.get("stage_hash", "")
    path = l2_day_path(cfg, day)
    l1_path = cfg.l1_cml_day_path(day)
    record = {"day": day_text, "path": display_path(path, cfg.repo_root), "messages": []}
    if not l1_path.exists():
        raise FileNotFoundError(f"{l1_path} not found: run ingest_cml for {day_text} first")
    # skip only if the flag settings, the L1 input and the period statistics are all unchanged
    want = (cfg.stage_hash("qc_flags"), existing_stage_hash(l1_path), stats_hash)
    if not force and existing_l2_hashes(path) == want:
        record["messages"].append(("info", f"{day_text}: L2 exists with the same stage, input and stats hashes, skipped"))
        return {**record, "status": "skipped_same_config", "runtime_s": round(time.perf_counter() - t0, 2)}
    flags, summary = flag_day(cfg, day, stats)
    written = write_l2(flags, path, cfg, stats_hash)
    with xr.open_dataset(path) as on_disk, xr.open_dataset(l1_path) as l1:
        on_disk.load()
        schema.validate_flags_l2(on_disk, l1)
        counts = {k: int(v) for k, v in on_disk.attrs.items() if k.startswith("n_")}
    record.update({"status": "written", "flag_counts": counts, "file_size_bytes": written["file_size_bytes"],
                   "runtime_s": round(time.perf_counter() - t0, 2), "_summary": summary})
    record["messages"].append(("info", f"{day_text}: flags written ({record['runtime_s']:.1f} s): "
                                       f"{ {k: v for k, v in counts.items() if k not in ('n_missing_rsl', 'n_missing_tsl')} }"))
    return record


def run_qc_stats(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """Pass 1 over ALL period days: per sub-link statistics (needed by qc_flags, e.g. the floor level)."""
    t0 = time.perf_counter()
    days = cfg.period_days()
    missing = [d for d in days if not cfg.l1_cml_day_path(d).exists()]
    if missing:
        raise FileNotFoundError(f"qc_stats needs every L1 CML day; missing {len(missing)}, e.g. {missing[0].date()}")
    l1_hashes = {existing_stage_hash(cfg.l1_cml_day_path(d)) for d in days}
    path = cfg.tables_dir / "sublink_period_stats.csv"
    record = {"step": "qc_stats", "path": display_path(path, cfg.repo_root)}
    if not force and path.exists():
        _, meta = read_period_stats(cfg)
        up_to_date = (meta.get("stage_hash") == cfg.stage_hash("qc_stats")
                      and set(meta.get("input_stage_hash", [])) == l1_hashes
                      and meta.get("days") == [f"{d:%Y-%m-%d}" for d in days])
        if up_to_date:
            log.info("qc_stats: %s is up to date, skipped (use --force to rebuild)", path.name)
            records.append({**record, "status": "skipped_same_config"})
            return
    table, hists = period_stats(cfg, days)
    write_period_stats(cfg, table, hists, days, l1_hashes)
    record.update({"status": "written", "n_rows": len(table), "n_days": len(days),
                   "runtime_s": round(time.perf_counter() - t0, 2)})
    log.info("qc_stats: wrote %s (%d sub-links, %.0f s)", path.name, len(table), record["runtime_s"])
    records.append(record)


def run_qc_flags(cfg: Config, days: list[pd.Timestamp], force: bool, log, records: list[dict],
                 workers: int = 1) -> None:
    """Flags per day -> L2 files; then the daily summary (rows of the written days) and the network-gap table.

    The network-gap table is derived from the L2 files; it is rebuilt only when a day was (re)written or it is
    missing, so an up-to-date run writes nothing.
    """
    read_period_stats(cfg)  # fail early (before starting workers) if pass 1 has not been run
    done = run_day_task(qc_flags_day, cfg, days, workers, force, log, records)
    summaries = [r["_summary"] for r in done if "_summary" in r]
    if summaries:
        upsert_daily_summary(cfg, pd.concat(summaries, ignore_index=True))
    if summaries or force or not (cfg.tables_dir / "network_gaps.csv").exists():
        write_network_gaps(cfg)


def run_qc_status(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """Per sub-link status from the daily summary. Refuses to run until every period day has been flagged.

    Skipped when the status table exists with the current stage hash and was built from the same daily summary.
    """
    import hashlib
    import json
    path = cfg.tables_dir / "qc_daily_summary.csv"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found: run qc_flags first")
    sig = hashlib.sha256(path.read_bytes()).hexdigest()
    out_path = cfg.tables_dir / "sublink_qc_status.csv"
    side = Path(str(out_path) + ".meta.json")
    if not force and side.exists():
        m = json.loads(side.read_text(encoding="utf-8"))
        if m.get("stage_hash") == cfg.stage_hash("qc_status") and m.get("input_signature", sig) == sig:
            log.info("qc_status: %s is up to date, skipped (use --force to rebuild)", out_path.name)
            records.append({"step": "qc_status", "status": "skipped_same_config"})
            return
    summary = pd.read_csv(path)
    days = [f"{d:%Y-%m-%d}" for d in cfg.period_days()]
    missing = sorted(set(days) - set(summary["day"]))
    if missing:
        raise SystemExit(f"qc_status needs qc_flags for every day of the period; {len(missing)} days missing "
                         f"(e.g. {missing[0]}). Run qc_flags --all first.")
    summary = summary[summary["day"].isin(days)]
    table = sublink_status(cfg, summary)
    out = write_status(cfg, table, days)
    m = json.loads(side.read_text(encoding="utf-8"))
    m["input_signature"] = sig
    side.write_text(json.dumps(m, indent=2), encoding="utf-8")
    counts = table["status"].value_counts().to_dict()
    log.info("qc_status: wrote %s: %s", out.name, counts)
    records.append({"step": "qc_status", "path": display_path(out, cfg.repo_root), "status": "written",
                    "status_counts": counts})
