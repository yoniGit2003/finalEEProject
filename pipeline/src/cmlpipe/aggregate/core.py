"""Aggregation to L3: pure, vectorised reductions (no file I/O; the stages are in cmlpipe.aggregate.stage).

Time convention of every product: UTC, tz-naive, a bin [t, t + dt) is labelled by its START t.
Every resolution divides a day, so a day can be aggregated on its own and no bin crosses midnight.

Vectorisation: a regular time axis of length n_bins * spb (spb = samples per bin) is reshaped to
(n_bins, spb, ...), and every statistic is a reduction over axis 1. No Python loop over samples.

CML (one L1 day + its L2 flags) -> per bin and variable v in {tl, rsl, tsl}:
  statistics of the VALID samples (not NaN and no blocking flag for v), n_valid_<v>, and the rule
  n_valid_<v> < min_count -> every statistic of v is NaN in that bin (n_valid is kept);
  flags_or = bitwise OR of all L2 flag bits in the bin; n_tsl_step = number of tsl_step samples.
Gauges -> SUM of the native amounts per bin (strict: NaN unless enough native steps are valid).
Radar  -> stamps are instantaneous samples; the stamp at t goes to the bin [t, t + dt) and rain RATES are
  averaged over the valid stamps of a bin.
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pandas as pd
import xarray as xr

from cmlpipe.io import schema

CML_VARS = ("tl", "rsl", "tsl")
REDUCERS = {
    "mean": np.nanmean,
    "median": np.nanmedian,
    "min": np.nanmin,
    "max": np.nanmax,
    "std": lambda a, axis: np.nanstd(a, axis=axis, ddof=0),  # population std, documented in the attributes
}
STAT_LONG = {"mean": "mean", "median": "median", "min": "minimum", "max": "maximum",
             "std": "standard deviation (ddof=0)"}
CML_UNITS = {"tl": "dB", "rsl": "dBm", "tsl": "dBm"}


# ---------------------------------------------------------------------------
# Time helpers
# ---------------------------------------------------------------------------
def res_timedelta(res: str) -> pd.Timedelta:
    """A resolution string as a Timedelta."""
    return pd.Timedelta(res)


def samples_per_bin(res: str, native: str) -> int:
    """How many native steps make one bin; the resolution must be a whole multiple of the native step."""
    r, n = pd.Timedelta(res), pd.Timedelta(native)
    if r < n or r % n != pd.Timedelta(0):
        raise ValueError(f"resolution {res} is not a whole multiple of the native step {native}")
    if pd.Timedelta("1D") % r != pd.Timedelta(0):
        raise ValueError(f"resolution {res} does not divide a day")
    return int(r // n)


def min_count(fraction: float, n_expected: int) -> int:
    """Smallest number of valid samples that satisfies n_valid / n_expected >= fraction.

    Computed in integers (with a tiny tolerance) so that e.g. 72 / 90 = 0.8 is not lost to float rounding.
    """
    return int(math.ceil(round(fraction * n_expected, 9)))


def bin_starts(start: pd.Timestamp, n_bins: int, res: str) -> np.ndarray:
    """Labels of n_bins consecutive bins starting at `start` (datetime64[ns])."""
    return pd.date_range(pd.Timestamp(start), periods=n_bins, freq=pd.Timedelta(res)).to_numpy().astype("datetime64[ns]")


def time_attrs(res: str, n_expected: int, extra: str = "") -> dict:
    """Attributes describing the time axis of an L3 product (written into every product)."""
    return {
        "time_convention": "interval_start",
        "time_zone": "UTC",
        "interval": res,
        "interval_s": int(pd.Timedelta(res).total_seconds()),
        "n_expected": int(n_expected),
        "time_meaning": f"bin [time, time + {res}), labelled by its start" + (f"; {extra}" if extra else ""),
    }


def count_dtype(n_expected: int):
    """Smallest unsigned integer type for counts up to n_expected."""
    return np.uint8 if n_expected <= np.iinfo(np.uint8).max else np.uint16


def _reshape(a: np.ndarray, spb: int) -> np.ndarray:
    if a.shape[0] % spb:
        raise ValueError(f"time length {a.shape[0]} is not a multiple of {spb}")
    return a.reshape((a.shape[0] // spb, spb) + a.shape[1:])


def bin_statistics(values: np.ndarray, valid: np.ndarray, spb: int, stats: list[str], min_n: int
                   ) -> tuple[dict, np.ndarray]:
    """Per-bin statistics of `values` over the samples where `valid` is True.

    values, valid : (time, ...) with time = n_bins * spb.
    Returns ({stat: (n_bins, ...) float array}, n_valid (n_bins, ...)). Where n_valid < min_n every statistic is NaN.
    """
    v = _reshape(np.where(valid, values, np.nan), spb)
    n_valid = _reshape(valid, spb).sum(axis=1)
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-NaN bins: NaN is the right answer
        for s in stats:
            out[s] = REDUCERS[s](v, axis=1)
    bad = n_valid < min_n
    for s in out:
        out[s] = np.where(bad, np.nan, out[s])
    return out, n_valid


def blocking_mask(flags: np.ndarray, names: list[str]) -> np.ndarray:
    """True where any of the named L2 flag bits is set."""
    unknown = [n for n in names if n not in schema.FLAGS]
    if unknown:
        raise ValueError(f"unknown flag names {unknown}; known: {list(schema.FLAGS)}")
    mask = np.uint16(sum(1 << schema.FLAGS[n] for n in names))
    return (flags & mask) != 0


def stats_for(acfg: dict, res: str) -> dict:
    """Statistics per variable at a resolution: the resolution's own entry, else 'default'."""
    return acfg["stats"].get(res, acfg["stats"]["default"])


# ---------------------------------------------------------------------------
# CML
# ---------------------------------------------------------------------------
def aggregate_cml_day(l1: xr.Dataset, l2: xr.Dataset, res: str, acfg: dict, native: str) -> xr.Dataset:
    """One L1 CML day + its L2 flags -> the L3 dataset at `res` (dims time, sublink_id, cml_id).

    l1 must hold the full day on the native grid (it does: L1 days are regular and complete, checked here).
    """
    day = pd.Timestamp(l1.time.values[0]).normalize()
    n_native = int(pd.Timedelta("1D") / pd.Timedelta(native))
    want = bin_starts(day, n_native, native)
    if l1.sizes["time"] != n_native or not np.array_equal(l1.time.values, want):
        raise ValueError(f"L1 day {day.date()} is not the full regular {native} grid of that day")
    if not np.array_equal(l2.time.values, l1.time.values):
        raise ValueError("L2 time axis differs from L1")
    order = schema.DIMS
    rsl = l1.rsl.transpose(*order).values.astype("float64")
    tsl = l1.tsl.transpose(*order).values.astype("float64")
    tl = l1.tl.transpose(*order).values.astype("float64")
    flags = l2[schema.FLAGS_VAR].transpose(*order).values

    spb = samples_per_bin(res, native)
    n_bins = n_native // spb
    min_n = min_count(acfg["min_valid_fraction"], spb)
    bf = acfg["blocking_flags"]
    valid = {"rsl": ~np.isnan(rsl) & ~blocking_mask(flags, bf["rsl"]),
             "tsl": ~np.isnan(tsl) & ~blocking_mask(flags, bf["tsl"])}
    valid["tl"] = valid["rsl"] & valid["tsl"] & ~np.isnan(tl)
    values = {"tl": tl, "rsl": rsl, "tsl": tsl}

    cdt = count_dtype(spb)
    data_vars = {}
    want_stats = stats_for(acfg, res)
    for v in CML_VARS:
        stats, n_valid = bin_statistics(values[v], valid[v], spb, want_stats[v], min_n)
        for s, a in stats.items():
            data_vars[f"{v}_{s}"] = (order, a.astype("float32"),
                                     {"units": CML_UNITS[v], "long_name": f"{STAT_LONG[s]} of valid {v} samples in the bin"})
        data_vars[f"n_valid_{v}"] = (order, n_valid.astype(cdt),
                                     {"long_name": f"number of valid {v} samples in the bin", "units": "1"})
    data_vars["flags_or"] = (order, np.bitwise_or.reduce(_reshape(flags, spb), axis=1).astype("uint16"),
                             {**schema.flag_attrs(), "long_name": "bitwise OR of all L2 flag bits of the samples in the bin"})
    tsl_step = (flags & np.uint16(1 << schema.FLAGS["tsl_step"])) != 0
    data_vars["n_tsl_step"] = (order, _reshape(tsl_step, spb).sum(axis=1).astype(cdt),
                               {"long_name": "number of samples flagged tsl_step in the bin", "units": "1"})

    coords = {"time": bin_starts(day, n_bins, res)}
    for name, c in l1.coords.items():
        if name != "time":
            coords[name] = c
    ds = xr.Dataset({k: (d, a, attrs) for k, (d, a, attrs) in data_vars.items()}, coords=coords)
    ds.attrs = {
        **time_attrs(res, spb, f"a {native} CML sample goes into the bin that contains its label (window start)"),
        "day_utc": f"{day:%Y-%m-%d}",
        "native_interval": native,
        "min_valid_fraction": acfg["min_valid_fraction"],
        "min_valid_count": min_n,
        "valid_rule": f"if n_valid_<v> < {min_n} (= ceil({acfg['min_valid_fraction']} * {spb})) all statistics of v "
                      "are NaN in that bin; n_valid_<v> is kept",
        "blocking_flags_rsl": " ".join(bf["rsl"]),
        "blocking_flags_tsl": " ".join(bf["tsl"]),
        "blocking_flags_tl": "union of the rsl and tsl lists (valid for tl = valid for rsl AND tsl)",
        "informational_flags": "tsl_step flatline (never exclude a sample)",
        "std_ddof": 0,
        "median_note": "the median is stored on purpose: an isolated unflagged dip barely moves it",
    }
    return ds


# ---------------------------------------------------------------------------
# Gauges
# ---------------------------------------------------------------------------
def aggregate_gauges(l1: xr.Dataset, res: str, native: str, period_start: pd.Timestamp, period_end: pd.Timestamp,
                     min_valid_fraction: float, peak_rate: bool = False, qc_worst: bool = False) -> xr.Dataset:
    """Gauge L1 (time, id) -> sums per bin at `res`, on the full period grid [period_start, period_end).

    amount_mm = SUM of the native amounts, NaN unless n_valid >= min_count (strict by default: a sum with a
    missing step would silently under-count); rate_mm_h = amount_mm / dt_hours;
    peak_1min_rate_mm_h (City, res >= 5 min) = max over the bin of native amount * 60, same NaN rule;
    qc_worst (SMHI) = max of the original qc code in the bin (0 = no code).
    """
    n_native = int((period_end - period_start) / pd.Timedelta(native))
    grid = bin_starts(period_start, n_native, native)
    l1 = l1.reindex(time=grid)  # L1 already is this grid (checked by its validator); reindex is a no-op guard
    amount = l1.rainfall_amount.transpose("time", "id").values.astype("float64")
    spb = samples_per_bin(res, native)
    n_bins = n_native // spb
    min_n = min_count(min_valid_fraction, spb)
    valid = ~np.isnan(amount)
    stats, n_valid = bin_statistics(amount, valid, spb, [], min_n)
    sums = _reshape(np.where(valid, amount, 0.0), spb).sum(axis=1)
    sums = np.where(n_valid >= min_n, sums, np.nan)
    dt_h = pd.Timedelta(res).total_seconds() / 3600.0
    cdt = count_dtype(spb)
    data_vars = {
        "amount_mm": (("time", "id"), sums.astype("float32"),
                      {"units": "mm", "long_name": "rainfall amount summed over the bin",
                       "comment": f"NaN unless n_valid >= {min_n} of {spb} native steps (strict: not filled)"}),
        "n_valid": (("time", "id"), n_valid.astype(cdt), {"units": "1", "long_name": "valid native steps in the bin"}),
        "rate_mm_h": (("time", "id"), (sums / dt_h).astype("float32"),
                      {"units": "mm h-1", "long_name": "mean rain rate over the bin, amount_mm / dt_hours"}),
    }
    if peak_rate and spb > 1:
        native_h = pd.Timedelta(native).total_seconds() / 3600.0
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            peak = np.nanmax(_reshape(amount, spb), axis=1) / native_h
        peak = np.where(n_valid >= min_n, peak, np.nan)
        data_vars["peak_1min_rate_mm_h"] = (("time", "id"), peak.astype("float32"),
                                            {"units": "mm h-1", "long_name": "largest native 1 min amount * 60 in the bin"})
    if qc_worst and "qc_code" in l1:
        q = l1.qc_code.transpose("time", "id").values
        data_vars["qc_worst"] = (("time", "id"), _reshape(q, spb).max(axis=1).astype("uint8"),
                                 {"long_name": "largest original SMHI qc code in the bin (0 = none)",
                                  **{k: v for k, v in l1.qc_code.attrs.items() if k.startswith("flag_")}})
    coords = {"time": bin_starts(period_start, n_bins, res)}
    for name, c in l1.coords.items():
        if name != "time":
            coords[name] = c
    ds = xr.Dataset({k: (d, a, at) for k, (d, a, at) in data_vars.items()}, coords=coords)
    ds.attrs = {**time_attrs(res, spb, f"sums of the {native} L1 amounts, whose labels are interval starts"),
                "native_interval": native, "min_valid_fraction": min_valid_fraction, "min_valid_count": min_n}
    return ds


# ---------------------------------------------------------------------------
# Radar
# ---------------------------------------------------------------------------
def floor_rain(dbz: np.ndarray, rain: np.ndarray, min_dbz: float) -> tuple[np.ndarray, int, int]:
    """Rain rate with dBZ < min_dbz set to 0 (NaN stays NaN). Returns (floored, n_set_to_zero, n_finite)."""
    finite = ~np.isnan(rain)
    below = finite & (dbz < min_dbz)
    return np.where(below, 0.0, rain).astype("float32"), int(below.sum()), int(finite.sum())


def stamps_to_grid(stamp_times: np.ndarray, values: np.ndarray, period_start: pd.Timestamp,
                   period_end: pd.Timestamp, native: str) -> tuple[np.ndarray, int]:
    """Place instantaneous stamps on the native bin grid of the period: stamp t -> bin [t, t + native).

    Bins without a stamp are NaN. Stamps outside [period_start, period_end) (the closing instant
    2015-09-01 00:00) belong to no bin of the period and are dropped; their number is returned.
    """
    n = int((period_end - period_start) / pd.Timedelta(native))
    grid = bin_starts(period_start, n, native)
    st = pd.DatetimeIndex(stamp_times).floor(native)  # stamps are on the grid already; floor = assignment rule
    inside = (st >= period_start) & (st < period_end)
    idx = ((st[inside] - period_start) // pd.Timedelta(native)).to_numpy().astype("int64")
    out = np.full((n,) + values.shape[1:], np.nan, dtype="float64")
    out[idx] = values[np.asarray(inside)]
    return out, int((~inside).sum())


def aggregate_stamps(values: dict, valid: np.ndarray, spb: int, min_valid_fraction: float,
                     how: dict) -> tuple[dict, np.ndarray]:
    """Stamps on the native grid (time, ...) -> bins: mean or max over the VALID stamps of each bin.

    values : {name: (time, ...) array}; valid : (time, ...) bool, one validity per stamp and column.
    how    : {name: "mean" | "max" | "mean_all"}; "mean_all" averages over all stamps (missing = 0), used for
             valid_weight_fraction, which must describe the whole bin.
    A bin with n_valid_stamps < min_count(min_valid_fraction, spb) is NaN for every "mean"/"max" variable.
    """
    min_n = min_count(min_valid_fraction, spb)
    n_valid = _reshape(valid, spb).sum(axis=1)
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for name, a in values.items():
            if how[name] == "mean_all":
                out[name] = _reshape(np.nan_to_num(a, nan=0.0), spb).mean(axis=1)
                continue
            r = _reshape(np.where(valid, a, np.nan), spb)
            red = np.nanmean(r, axis=1) if how[name] == "mean" else np.nanmax(r, axis=1)
            out[name] = np.where(n_valid >= min_n, red, np.nan)
    return out, n_valid
