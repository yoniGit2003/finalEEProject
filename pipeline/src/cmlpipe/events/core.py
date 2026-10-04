"""Rain events, dry periods and the frozen train/test split: pure functions (stages in cmlpipe.events.stage).

WHY events come from the references only
----------------------------------------
Events and dry periods are defined from the radar and the 10 City gauges, never from CML signals. A model
trained later on CML data must not have its training data selected by its own input (that would be
circular, like choosing links by how well the power law fits). The SMHI gauge GbgA is not used either:
it is kept as an independent check of the event definition.

States per 5 min bin (config events.*)
--------------------------------------
  radar valid : at least radar_min_valid_cell_fraction of the domain cells have data
  radar wet   : radar valid and wet-area fraction >= min_wet_area_fraction
                (wet cell = floored rain rate >= wet_threshold_mmh; floored = 0 below 0 dBZ)
  gauge wet   : >= min_gauges_wet City gauges with a 5 min amount >= city_gauge_wet_amount_mm
  wet         : radar wet OR gauge wet
  dry         : not wet, radar valid, and wet-area fraction <= dry_max_wet_area_fraction
  unknown     : everything else, i.e. radar missing (we cannot certify dry) or a few wet cells below the
                event threshold. unknown NEVER counts as dry; for event merging it counts as "not wet".
A City gauge whose 5 min amount is NaN (e.g. Barl's 48 h gap) is "gauge unknown", not "gauge dry".

Events: wet bins merged when separated by < min_dry_gap_h of non-wet time; events shorter than
min_event_duration_min are dropped (and counted). event_id = "E" + start, e.g. "E20150728_1455", so ids stay
stable when thresholds elsewhere change.

Evaluation scheme (frozen split, tables/split_events.csv and split_links.csv)
--------------------------------------------------------------------------
Events: train / val / test = 60 / 20 / 20 %, stratified by intensity class. Links: train / test = 80 / 20 %,
stratified by length class x frequency band; parallel links (same two sites) always share a split; links with
both sub-links "unusable" are "excluded".
  * fit models on      train events x train links (tune on val events x train links)
  * unseen events:     test events  x train links
  * unseen links:      train events x test links
  * both unseen:       test events  x test links
CAVEAT: neighbouring links see the same weather and many share sites, so the link split is NOT a spatial
hold-out. A spatial hold-out must be built later from site ids / tables/link_neighbours.csv if needed.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from cmlpipe import geometry

STATE_CODES = {"dry": 0, "wet": 1, "unknown": 2}
STATE_NAMES = {v: k for k, v in STATE_CODES.items()}


# ---------------------------------------------------------------------------
# Domain and reference series
# ---------------------------------------------------------------------------
def domain_mask(cell_x_3006: np.ndarray, cell_y_3006: np.ndarray, links: pd.DataFrame, buffer_m: float) -> np.ndarray:
    """(ny, nx) bool: radar cells whose CENTRE lies within buffer_m of any link path (distances in EPSG:3006)."""
    px, py = cell_x_3006.ravel()[:, None], cell_y_3006.ravel()[:, None]
    d, _ = geometry.point_to_path(px, py, links.site_0_x_3006_m.to_numpy()[None], links.site_0_y_3006_m.to_numpy()[None],
                                  links.site_1_x_3006_m.to_numpy()[None], links.site_1_y_3006_m.to_numpy()[None])
    return (d.min(axis=1) <= buffer_m).reshape(cell_x_3006.shape)


def radar_domain_series(floored: np.ndarray, mask: np.ndarray, wet_threshold: float) -> dict:
    """Per radar stamp over the domain cells: valid cell fraction, wet-area fraction (of VALID cells), mean rate.

    floored : (time, ny, nx) rain rate with the 0 dBZ floor applied (NaN = missing).
    """
    v = floored[:, mask]  # (time, n_domain)
    valid = ~np.isnan(v)
    n_valid = valid.sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        wet = (np.where(valid, v, 0.0) >= wet_threshold).sum(axis=1)
        frac_wet = np.where(n_valid > 0, wet / n_valid, np.nan)
        mean = np.where(n_valid > 0, np.nansum(v, axis=1) / n_valid, np.nan)
    return {"radar_valid_cell_fraction": n_valid / mask.sum(), "radar_wet_area_fraction": frac_wet,
            "domain_mean_rain_mmh": mean, "radar_n_wet_cells": wet.astype(float)}


def gauge_series(amount_5min: np.ndarray, rate_5min: np.ndarray, wet_amount_mm: float) -> dict:
    """City gauges, 5 min bins (time, gauge): number valid, number wet, max rate. NaN amount = gauge unknown."""
    valid = ~np.isnan(amount_5min)
    # 0.1 mm stored as float32 is 0.1000000015; compare with a tolerance far below the 0.1 mm quantum
    wet = valid & (np.where(valid, amount_5min, 0.0) >= wet_amount_mm - 1e-6)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        mx = np.nanmax(rate_5min, axis=1)
    return {"n_city_gauges_valid": valid.sum(axis=1), "n_city_gauges_wet": wet.sum(axis=1), "max_city_rate_mmh": mx}


def classify(ref: pd.DataFrame, ecfg: dict) -> pd.DataFrame:
    """Add radar_valid, radar_wet, gauge_wet, state (0 dry, 1 wet, 2 unknown) to the reference table."""
    ref = ref.copy()
    ref["radar_valid"] = ref.radar_valid_cell_fraction >= ecfg["radar_min_valid_cell_fraction"]
    waf = ref.radar_wet_area_fraction.fillna(0.0)
    ref["radar_wet"] = ref.radar_valid & (waf >= ecfg["min_wet_area_fraction"])
    ref["gauge_wet"] = ref.n_city_gauges_wet >= ecfg["min_gauges_wet"]
    wet = ref.radar_wet | ref.gauge_wet
    dry = ~wet & ref.radar_valid & (waf <= ecfg["dry_max_wet_area_fraction"])
    ref["state"] = np.where(wet, STATE_CODES["wet"], np.where(dry, STATE_CODES["dry"], STATE_CODES["unknown"])).astype("int8")
    return ref


# ---------------------------------------------------------------------------
# Runs, events, dry periods
# ---------------------------------------------------------------------------
def runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """[start, end) index pairs of the True runs of a 1-D bool array."""
    m = np.concatenate([[False], np.asarray(mask, bool), [False]])
    d = np.diff(m.astype(int))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def find_events(state: np.ndarray, step: pd.Timedelta, min_dry_gap: pd.Timedelta, min_duration: pd.Timedelta
                ) -> tuple[list[dict], int]:
    """Merge wet runs separated by less than min_dry_gap of non-wet time; drop events shorter than min_duration.

    Returns ([{i0, i1, merged_from}], n_dropped) with [i0, i1) the bin range from the first to the last wet bin.
    """
    wet_runs = runs(state == STATE_CODES["wet"])
    gap_bins = int(np.ceil(min_dry_gap / step))  # a gap of >= this many bins separates two events
    merged = []
    for a, b in wet_runs:
        if merged and a - merged[-1]["i1"] < gap_bins:
            merged[-1]["i1"] = b
            merged[-1]["merged_from"] += 1
        else:
            merged.append({"i0": a, "i1": b, "merged_from": 1})
    min_bins = int(np.ceil(min_duration / step))
    kept = [e for e in merged if e["i1"] - e["i0"] >= min_bins]
    return kept, len(merged) - len(kept)


def touches_edge(i0: int, i1: int, n: int, gap_bins: int) -> bool:
    """True when less than one min_dry_gap of non-wet time separates the event from the period start or end:
    the rain may have begun before the data start (or continued after it), so its start (end) is not known."""
    return bool(i0 < gap_bins or n - i1 < gap_bins)


def event_id(start: pd.Timestamp) -> str:
    """Stable event id 'E<start YYYYmmdd_HHMM>'."""
    return "E" + pd.Timestamp(start).strftime("%Y%m%d_%H%M")


def dilate(mask: np.ndarray, n: int) -> np.ndarray:
    """True within n bins of any True (both directions)."""
    if n <= 0:
        return np.asarray(mask, bool).copy()
    k = np.convolve(np.asarray(mask, int), np.ones(2 * n + 1, int), mode="same")
    return k > 0


def certified_dry(state: np.ndarray, margin_bins: int) -> np.ndarray:
    """Dry bins that are more than margin_bins away from any wet bin. unknown is never dry."""
    wet = state == STATE_CODES["wet"]
    return (state == STATE_CODES["dry"]) & ~dilate(wet, margin_bins)


def context(state: np.ndarray, ev: list[dict], k: int, cap_bins: int) -> dict:
    """Weather between event k and its neighbours, looking back from its start and forward from its end.

    gap_*  : non-wet bins back to the previous event's end / forward to the next event's start (or the period
             edge), capped at cap_bins.
    dry_*  : how many of those bins are certified dry (state dry). NOT required to be contiguous: an event start
             is usually preceded by a few "unknown" bins (rain on < min_wet_area_fraction of the domain as it moves
             in), so a contiguous count would be 0 for most events and say nothing.
    """
    i0, i1 = ev[k]["i0"], ev[k]["i1"]
    lo = max(ev[k - 1]["i1"] if k > 0 else 0, i0 - cap_bins)
    hi = min(ev[k + 1]["i0"] if k + 1 < len(ev) else state.size, i1 + cap_bins)
    dry = state == STATE_CODES["dry"]
    return {"gap_before": i0 - lo, "dry_before": int(dry[lo:i0].sum()),
            "gap_after": hi - i1, "dry_after": int(dry[i1:hi].sum())}


# ---------------------------------------------------------------------------
# Sensitivity grid (calibration)
# ---------------------------------------------------------------------------
def sensitivity(ref: pd.DataFrame, ecfg: dict, step: pd.Timedelta) -> pd.DataFrame:
    """Number of events and total event hours for a grid of min_wet_area_fraction x min_dry_gap_h."""
    rows = []
    g = ecfg["calibration_grid"]
    for f in g["min_wet_area_fraction"]:
        st = classify(ref, {**ecfg, "min_wet_area_fraction": f}).state.to_numpy()
        for gap in g["min_dry_gap_h"]:
            ev, dropped = find_events(st, step, pd.Timedelta(hours=gap), pd.Timedelta(minutes=ecfg["min_event_duration_min"]))
            hours = sum(e["i1"] - e["i0"] for e in ev) * step / pd.Timedelta("1h")
            rows.append({"min_wet_area_fraction": f, "min_dry_gap_h": gap, "n_events": len(ev), "event_hours": hours,
                         "n_dropped_short": dropped, "wet_bin_share": float((st == 1).mean())})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Robust statistics for the dry-weather table
# ---------------------------------------------------------------------------
def _nan(f, a, *args, **kw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return f(a, *args, **kw)


def scaled_mad(a: np.ndarray, axis=0) -> np.ndarray:
    """1.4826 * median absolute deviation: equals the std for Gaussian noise but ignores rare outliers."""
    med = _nan(np.nanmedian, a, axis=axis, keepdims=True)
    return 1.4826 * _nan(np.nanmedian, np.abs(a - med), axis=axis)


def dry_stats(times: pd.DatetimeIndex, dry: np.ndarray, core: dict, min_days: int = 5,
              min_bins_per_day: int = 12) -> dict:
    """Descriptive statistics per column (sublink) over the certified-dry bins only.

    core : {"tl_mean", "tl_std", "rsl_mean", "rsl_std"} arrays (time, column). Bins outside `dry` are ignored
    entirely (they are set to NaN first), so no wet value can leak in.
    A "dry day" is a UTC day with >= min_bins_per_day valid dry tl bins (1 h at 5 min).
    """
    sel = {k: np.where(dry[:, None], v, np.nan) for k, v in core.items()}
    tl, rsl = sel["tl_mean"], sel["rsl_mean"]
    n_dry = int(dry.sum())
    n_valid = np.isfinite(tl).sum(axis=0)
    out = {
        "n_dry_bins": np.full(tl.shape[1], n_dry), "n_valid_tl": n_valid,
        "valid_fraction": n_valid / max(n_dry, 1),
        "tl_median": _nan(np.nanmedian, tl, axis=0), "tl_mean": _nan(np.nanmean, tl, axis=0),
        "tl_std": _nan(np.nanstd, tl, axis=0), "tl_scaled_mad": scaled_mad(tl),
        "tl_iqr": _nan(np.nanpercentile, tl, 75, axis=0) - _nan(np.nanpercentile, tl, 25, axis=0),
        "tl_p1": _nan(np.nanpercentile, tl, 1, axis=0), "tl_p99": _nan(np.nanpercentile, tl, 99, axis=0),
        "tl_intra_bin_std_median": _nan(np.nanmedian, sel["tl_std"], axis=0),
        "rsl_median": _nan(np.nanmedian, rsl, axis=0), "rsl_std": _nan(np.nanstd, rsl, axis=0),
        "rsl_scaled_mad": scaled_mad(rsl), "rsl_intra_bin_std_median": _nan(np.nanmedian, sel["rsl_std"], axis=0),
    }
    day = times.normalize()
    days = day.unique()
    daily = np.full((len(days), tl.shape[1]), np.nan)
    counts = np.zeros((len(days), tl.shape[1]), int)
    for k, d in enumerate(days):
        rows = np.asarray(day == d)
        counts[k] = np.isfinite(tl[rows]).sum(axis=0)
        daily[k] = _nan(np.nanmedian, tl[rows], axis=0)
    daily = np.where(counts >= min_bins_per_day, daily, np.nan)
    n_days = np.isfinite(daily).sum(axis=0)
    hourly = np.vstack([_nan(np.nanmedian, tl[np.asarray(times.hour == h)], axis=0) for h in range(24)])
    out["n_dry_days"] = n_days
    out["tl_diurnal_range"] = np.where(n_days >= min_days, _nan(np.nanmax, hourly, axis=0) - _nan(np.nanmin, hourly, axis=0), np.nan)
    out["tl_daily_median_drift_std"] = np.where(n_days >= min_days, _nan(np.nanstd, daily, axis=0), np.nan)
    return out


# ---------------------------------------------------------------------------
# Split
# ---------------------------------------------------------------------------
def _allocate(n: int, fractions: list[float]) -> list[int]:
    """Integer counts for n items that sum to n (largest remainder), in the order of `fractions`."""
    raw = np.array(fractions) * n
    base = np.floor(raw).astype(int)
    for k in np.argsort(-(raw - base), kind="stable")[: n - base.sum()]:
        base[k] += 1
    return base.tolist()


def split_events(events: pd.DataFrame, fractions: dict, seed: int) -> pd.DataFrame:
    """event_id -> split, stratified by intensity_class; seeded permutation of the SORTED ids per class."""
    rng = np.random.default_rng(seed)
    names, fr = list(fractions), [fractions[k] for k in fractions]
    rows = []
    for cls in sorted(events.intensity_class.unique()):
        ids = sorted(events.loc[events.intensity_class == cls, "event_id"])
        perm = [ids[i] for i in rng.permutation(len(ids))]
        counts = _allocate(len(ids), fr)
        k = 0
        for name, c in zip(names, counts):
            rows += [{"event_id": e, "intensity_class": cls, "split_event": name} for e in perm[k:k + c]]
            k += c
    return pd.DataFrame(rows).sort_values("event_id").reset_index(drop=True)


def length_class(length_m: float, edges_km: list[float]) -> str:
    """Length class label (<1 km, 1-3 km, >3 km by default)."""
    km = length_m / 1000.0
    if km < edges_km[0]:
        return f"<{edges_km[0]:g}km"
    if km <= edges_km[1]:
        return f"{edges_km[0]:g}-{edges_km[1]:g}km"
    return f">{edges_km[1]:g}km"


def split_links(links: pd.DataFrame, status: pd.DataFrame, scfg: dict) -> pd.DataFrame:
    """cml_id -> split (train/test/excluded), stratified by length class x frequency band, hop groups kept together.

    links  : link_geometry rows (cml_id, site_0_id, site_1_id, length_m, frequency_ghz_sublink_1).
    status : sublink_qc_status rows (cml_id, sublink_id, status).
    A hop group = all links with the same unordered pair of sites (parallel links). The group's stratum is
    that of its smallest cml_id. Links with both sub-links "unusable" are "excluded" and leave their group.
    """
    rng = np.random.default_rng(scfg["seed"])
    lk = links.sort_values("cml_id").copy()
    lk["hop"] = [f"{min(a, b)}-{max(a, b)}" for a, b in zip(lk.site_0_id, lk.site_1_id)]
    lk["length_class"] = [length_class(v, scfg["length_classes_km"]) for v in lk.length_m]
    lk["frequency_band"] = np.where(lk.frequency_ghz_sublink_1 < scfg["frequency_split_ghz"],
                                    f"<{scfg['frequency_split_ghz']:g}GHz", f">={scfg['frequency_split_ghz']:g}GHz")
    unusable = status.assign(u=status.status == "unusable").groupby("cml_id").u.all()
    lk["excluded"] = lk.cml_id.map(unusable).fillna(False).astype(bool)
    lk["stratum"] = lk.length_class + " x " + lk.frequency_band
    use = lk[~lk.excluded]
    groups = use.groupby("hop").cml_id.min().reset_index().merge(use[["cml_id", "stratum"]], on="cml_id")
    assign = {}
    names = list(scfg["link_fractions"])
    for st in sorted(groups.stratum.unique()):
        hops = sorted(groups.loc[groups.stratum == st, "hop"])
        perm = [hops[i] for i in rng.permutation(len(hops))]
        counts = _allocate(len(hops), [scfg["link_fractions"][k] for k in names])
        k = 0
        for name, c in zip(names, counts):
            for h in perm[k:k + c]:
                assign[h] = name
            k += c
    lk["split_link"] = [("excluded" if ex else assign[h]) for ex, h in zip(lk.excluded, lk.hop)]
    return lk[["cml_id", "hop", "length_class", "frequency_band", "stratum", "excluded", "split_link"]].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Validators
# ---------------------------------------------------------------------------
class EventsValidationError(ValueError):
    """Raised when an events / periods / split table breaks a rule."""
    pass


def validate_events(ev: pd.DataFrame, split: pd.DataFrame | None = None) -> None:
    """Unique ids, start < end, no overlaps, and (if given) a split row for every event."""
    errs = []
    if ev.event_id.duplicated().any():
        errs.append("duplicate event ids")
    s, e = pd.to_datetime(ev.start), pd.to_datetime(ev.end)
    if not (s < e).all():
        errs.append("an event has start >= end")
    o = np.argsort(s.to_numpy())
    if (s.to_numpy()[o][1:] < e.to_numpy()[o][:-1]).any():
        errs.append("events overlap")
    if split is not None:
        missing = sorted(set(ev.event_id) - set(split.event_id))
        if missing:
            errs.append(f"{len(missing)} events have no split row (e.g. {missing[0]}); the split is frozen: "
                        "rerun `--steps split --force-split` only if you really mean to redraw it")
    if errs:
        raise EventsValidationError("events.csv: " + "; ".join(errs))


def validate_periods(dp: pd.DataFrame, what: str = "dry_periods.csv") -> None:
    """start < end and no overlapping periods."""
    s, e = pd.to_datetime(dp.start), pd.to_datetime(dp.end)
    if not (s < e).all():
        raise EventsValidationError(f"{what}: start >= end")
    o = np.argsort(s.to_numpy())
    if (s.to_numpy()[o][1:] < e.to_numpy()[o][:-1]).any():
        raise EventsValidationError(f"{what}: periods overlap")


def validate_unique(df: pd.DataFrame, keys: list[str], what: str) -> None:
    """No duplicate rows for the given key columns."""
    if df.duplicated(keys).any():
        raise EventsValidationError(f"{what}: duplicate rows for {keys}")
