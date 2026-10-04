"""QC flags for CML signals: pure functions, no I/O.

FLAG, DON'T FIX: every function returns a boolean mask (True = flagged) with the same shape as its
input; nothing is removed, smoothed or filled. The daily loop and file handling are in cmlpipe.qc.stage.

All signal arrays are 2-D (time, n): time on axis 0, one column per sub-link. Every function works
on all columns at once (vectorized); there are no per-sample Python loops. The time step is regular
(L1 guarantees it), so durations are converted to sample counts.

Uses only the CML signals and their metadata, never radar, gauges or weather, so the flags cannot
be biased by the references used later for validation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from cmlpipe.qc.flags import FLAGS, decode_flags

MAD_TO_SIGMA = 1.4826  # scaled MAD: for Gaussian noise 1.4826 * MAD estimates the standard deviation


# ---------------------------------------------------------------------------
# Run helpers
# ---------------------------------------------------------------------------
def _ffill(a: np.ndarray) -> np.ndarray:
    return pd.DataFrame(a).ffill().to_numpy()


def _bfill(a: np.ndarray) -> np.ndarray:
    return pd.DataFrame(a).bfill().to_numpy()


def run_lengths(mask: np.ndarray) -> np.ndarray:
    """For every True element, the length of the run of consecutive True values it belongs to (0 elsewhere).

    A run is broken by any False, so a NaN that was turned into False breaks runs.
    """
    mask = np.asarray(mask, dtype=bool)
    t = np.arange(mask.shape[0], dtype="float64")[:, None] * np.ones((1, mask.shape[1]))
    prev = np.vstack([np.zeros((1, mask.shape[1]), bool), mask[:-1]])
    nxt = np.vstack([mask[1:], np.zeros((1, mask.shape[1]), bool)])
    start = _ffill(np.where(mask & ~prev, t, np.nan))
    end = _bfill(np.where(mask & ~nxt, t, np.nan))
    return np.where(mask, end - start + 1, 0).astype("int64")


def run_bounds(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Index of the first and last element of the run each True element belongs to (-1 elsewhere)."""
    mask = np.asarray(mask, dtype=bool)
    t = np.arange(mask.shape[0], dtype="float64")[:, None] * np.ones((1, mask.shape[1]))
    prev = np.vstack([np.zeros((1, mask.shape[1]), bool), mask[:-1]])
    nxt = np.vstack([mask[1:], np.zeros((1, mask.shape[1]), bool)])
    start = _ffill(np.where(mask & ~prev, t, np.nan))
    end = _bfill(np.where(mask & ~nxt, t, np.nan))
    return (np.where(mask, start, -1).astype("int64"), np.where(mask, end, -1).astype("int64"))


def n_samples(duration: str | pd.Timedelta, interval: pd.Timedelta) -> int:
    """Number of samples covering `duration` at the given step, e.g. 10 min at 10 s -> 60."""
    return int(round(pd.Timedelta(duration) / interval))


# ---------------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------------
def flag_missing(x: np.ndarray) -> np.ndarray:
    """NaN samples (missing_rsl / missing_tsl)."""
    return np.isnan(x)


def flag_out_of_range(rsl: np.ndarray, tsl: np.ndarray, rsl_range, tsl_range) -> np.ndarray:
    """rsl or tsl finite but outside its physical range [lo, hi]."""
    def outside(v, lo, hi):
        return np.isfinite(v) & ((v < lo) | (v > hi))
    return outside(rsl, *rsl_range) | outside(tsl, *tsl_range)


def flag_rx_floor(rsl: np.ndarray, floor: np.ndarray, median: np.ndarray, margin_db: float, depth_db: float,
                  min_len: int, receiver_floor_max_dbm: float) -> np.ndarray:
    """Receiver clipping: rsl within margin_db of the sub-link's floor level for >= min_len consecutive samples.

    floor, median: per column (the sub-link's robust floor level and period median, from qc_stats).
    Only sub-links whose floor is a RECEIVER floor can be flagged: it must lie at least depth_db below
    their median and at or below receiver_floor_max_dbm (near the receiver sensitivity). A sub-link
    whose lowest dense level is, say, -70 dBm has no clipping; its lowest values are real fades.
    A NaN breaks a run, and a real fade that touches the floor for less than min_len samples is not
    flagged. A clipped value is a lower bound, not a measurement, which is why the flag is blocking.
    """
    floor = np.asarray(floor, dtype="float64")[None, :]
    median = np.asarray(median, dtype="float64")[None, :]
    deep_enough = (np.isfinite(floor) & np.isfinite(median) & (floor <= median - depth_db)
                   & (floor <= receiver_floor_max_dbm))
    at_floor = np.isfinite(rsl) & (rsl <= floor + margin_db) & deep_enough
    return run_lengths(at_floor) >= min_len


def spike_deviation(rsl: np.ndarray, window: int) -> tuple[np.ndarray, np.ndarray]:
    """Deviation from a centred rolling median and the rolling scaled MAD of that deviation.

    The MAD is the rolling median of |deviation| (a common, fast approximation of the per-window
    MAD). During a fade the MAD grows, which raises the threshold: that is what keeps smooth fades
    from being called spikes. Windows need at least half their samples valid, else NaN.
    """
    df = pd.DataFrame(rsl)
    min_periods = window // 2 + 1
    med = df.rolling(window, center=True, min_periods=min_periods).median().to_numpy()
    dev = rsl - med
    mad = pd.DataFrame(np.abs(dev)).rolling(window, center=True, min_periods=min_periods).median().to_numpy()
    return dev, MAD_TO_SIGMA * mad


def flag_spike(rsl: np.ndarray, window: int, n_mad: float, min_abs_db: float, max_len: int,
               neighbour_agree_db: float, neighbour_max_dev_db: float) -> np.ndarray:
    """Isolated short outliers in rsl.

    A sample is an outlier when |rsl - rolling median| > max(n_mad * scaled MAD, min_abs_db).
    A run of outliers is a spike only if
      - it has at most max_len samples,
      - the samples just before and just after the run exist, are not outliers, and agree with
        each other within neighbour_agree_db (the signal returns to where it was), and
      - those two neighbours are themselves CALM: within neighbour_max_dev_db of the rolling
        median. This is what separates a spike from the sharp bottom of a rain fade: at a V-shaped
        fade minimum the neighbours also agree with each other (they sit at the same depth on the
        two walls), but they lie far below the median. (Pilot gallery, first version of the rule.)
    """
    dev, scale = spike_deviation(rsl, window)
    threshold = np.fmax(n_mad * np.nan_to_num(scale, nan=0.0), min_abs_db)
    outlier = np.isfinite(dev) & (np.abs(dev) > threshold)
    length = run_lengths(outlier)
    start, end = run_bounds(outlier)
    n_t = rsl.shape[0]
    before_i = np.clip(start - 1, 0, n_t - 1)
    after_i = np.clip(end + 1, 0, n_t - 1)
    before = np.take_along_axis(rsl, before_i, axis=0)
    after = np.take_along_axis(rsl, after_i, axis=0)
    before_out = np.take_along_axis(outlier, before_i, axis=0)
    after_out = np.take_along_axis(outlier, after_i, axis=0)
    inside = (start >= 1) & (end <= n_t - 2)  # the run has a sample on both sides
    dev_before = np.take_along_axis(dev, before_i, axis=0)
    dev_after = np.take_along_axis(dev, after_i, axis=0)
    agree = inside & np.isfinite(before) & np.isfinite(after) & ~before_out & ~after_out & \
        (np.abs(before - after) <= neighbour_agree_db)
    calm = (np.abs(dev_before) <= neighbour_max_dev_db) & (np.abs(dev_after) <= neighbour_max_dev_db)
    return outlier & (length <= max_len) & agree & calm


def flag_flatline(rsl: np.ndarray, min_len: int) -> np.ndarray:
    """rsl exactly constant over >= min_len VALID samples.

    NaN samples do not break the run (the value is compared with the previous valid sample, so a
    stuck value with a few gaps is still one run), but they do not count towards its length either:
    a value seen before and after a long outage is not a long flatline. Only valid samples are flagged.
    """
    valid = np.isfinite(rsl)
    t = np.cumsum(valid, axis=0).astype("float64")  # index counted in valid samples only
    prev_valid = np.vstack([np.full((1, rsl.shape[1]), np.nan), _ffill(rsl)[:-1]])
    next_valid = np.vstack([_bfill(rsl)[1:], np.full((1, rsl.shape[1]), np.nan)])
    start = valid & ~(prev_valid == rsl)   # first sample of a constant segment
    end = valid & ~(next_valid == rsl)     # last sample of a constant segment
    first = _ffill(np.where(start, t, np.nan))
    last = _bfill(np.where(end, t, np.nan))
    span = last - first + 1
    return valid & (span >= min_len)


def flag_tsl_step(tsl: np.ndarray, min_db: float) -> np.ndarray:
    """The sample where tsl differs by >= min_db from the previous VALID tsl sample."""
    prev_valid = np.vstack([np.full((1, tsl.shape[1]), np.nan), _ffill(tsl)[:-1]])
    return np.isfinite(tsl) & np.isfinite(prev_valid) & (np.abs(tsl - prev_valid) >= min_db - 1e-6)


def alive_columns(rsl_day: np.ndarray, min_valid_fraction: float) -> np.ndarray:
    """Columns (sub-links) with at least min_valid_fraction valid rsl samples in this block (one day)."""
    return np.isfinite(rsl_day).mean(axis=0) >= min_valid_fraction


def network_nan_fraction(rsl: np.ndarray, alive: np.ndarray) -> np.ndarray:
    """Per time step: fraction of the alive sub-links whose rsl is NaN (NaN if none is alive)."""
    alive = np.asarray(alive, dtype=bool)
    if not alive.any():
        return np.full(rsl.shape[0], np.nan)
    return np.isnan(rsl[:, alive]).mean(axis=1)


def flag_network_gap(rsl: np.ndarray, nan_fraction: np.ndarray, threshold: float, alive: np.ndarray) -> np.ndarray:
    """NaN samples of ALIVE sub-links at time steps where more than `threshold` of the alive sub-links are NaN.

    It only explains WHY a sample is missing (the collector, not the link); it never flags a valid sample.
    A sub-link that is not alive that day (e.g. dead for the whole period) is missing for its own
    reason, so its NaN samples are not attributed to the network gap.
    """
    gap = np.nan_to_num(nan_fraction, nan=0.0) > threshold
    return np.isnan(rsl) & gap[:, None] & np.asarray(alive, dtype=bool)[None, :]


# ---------------------------------------------------------------------------
# Bitmask
# ---------------------------------------------------------------------------
def combine(masks: dict[str, np.ndarray]) -> np.ndarray:
    """Boolean masks by flag name -> one uint16 bitmask (bit numbers from schema.FLAGS)."""
    out = None
    for name, mask in masks.items():
        bit = np.uint16(1 << FLAGS[name])
        part = np.where(mask, bit, np.uint16(0)).astype("uint16")
        out = part if out is None else (out | part)
    return out


decode = decode_flags  # uint16 bitmask -> {flag name: boolean mask}; defined in schema


def blocking_mask(flags: np.ndarray, blocking: list[str]) -> np.ndarray:
    """True where any of the blocking flags is set."""
    bits = np.uint16(sum(1 << FLAGS[name] for name in blocking))
    return (np.asarray(flags) & bits) != 0
