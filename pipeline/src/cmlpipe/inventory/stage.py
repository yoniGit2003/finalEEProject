"""Stages `significance` and `inventory` (pure logic in cmlpipe.inventory.core).

significance -> tables/event_significance.csv, tables/event_runs.csv
inventory    -> tables/inventory/sublinks.csv, links.csv, gauges.csv, quality_status.csv

Additive only: events.csv, split_*.csv, event_links.csv and L1/L2/L3 are read, never written. Inputs are the
step 3-6 tables (qc_daily_summary, sublink_period_stats + histograms, sublink_dry_stats, split_links,
link_geometry, link_gauge_matches, link_neighbours, gauge_radar_pixel), the network reference series, the
consolidated 5 min CML core file and the per-day 5 min L3 files (only n_tsl_step, read with netCDF4).
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

from cmlpipe.aggregate import stage as aggregate_runner
from cmlpipe import events
from cmlpipe.events import stage as events_runner
from cmlpipe import geometry
from cmlpipe import inventory
from cmlpipe.config import Config
from cmlpipe.provenance import write_sidecar


def inv_dir(cfg: Config) -> Path:
    """Folder of the inventory tables."""
    return cfg.output_root / cfg["paths"]["inventory_dir"]


def out_paths(cfg: Config) -> dict[str, Path]:
    """All output paths of the step 7 stages."""
    d = inv_dir(cfg)
    return {"event_significance": cfg.tables_dir / "event_significance.csv", "event_runs": cfg.tables_dir / "event_runs.csv",
            "sublinks": d / "sublinks.csv", "links": d / "links.csv", "gauges": d / "gauges.csv",
            "quality_status": d / "quality_status.csv"}


def _meta(path: Path) -> dict:
    side = Path(str(path) + ".meta.json")
    return json.loads(side.read_text(encoding="utf-8")) if side.exists() else {}


def _write(cfg: Config, path: Path, df: pd.DataFrame, stage: str, description: str, extra: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, index=False)
    tmp.replace(path)
    write_sidecar(path, cfg, stage, extra={"description": description, "n_rows": len(df), "columns": list(df.columns), **extra})


def _signature(paths: list[Path], attrs_of: list[Path] = ()) -> str:
    h = hashlib.sha256()
    for p in paths:
        h.update(Path(p).read_bytes())
    for p in attrs_of:
        a = aggregate_runner.file_attrs(p)
        h.update(f"{p.name}{a.get('stage_hash')}{a.get('created_utc')}{a.get('input_signature')}".encode())
    return h.hexdigest()


def _up_to_date(cfg: Config, stage: str, paths: list[Path], sig: str) -> bool:
    return all(p.exists() and _meta(p).get("stage_hash") == cfg.stage_hash(stage) and _meta(p).get("input_signature") == sig
               for p in paths)


# ---------------------------------------------------------------------------
# Stage: significance (7.1, 7.2)
# ---------------------------------------------------------------------------
def run_significance(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """Stage significance: event_significance.csv and event_runs.csv."""
    t0 = time.perf_counter()
    scfg = cfg["event_significance"]
    P = out_paths(cfg)
    ref_path = events_runner.reference_path(cfg)
    sig = _signature([events_runner.table_path(cfg, "events"), events_runner.table_path(cfg, "split_events")],
                     [ref_path, aggregate_runner.gauge_path(cfg, "city", "5min")])
    targets = [P["event_significance"], P["event_runs"]]
    if not force and _up_to_date(cfg, "significance", targets, sig):
        log.info("significance: up to date, skipped (use --force to rebuild)")
        records.append({"step": "significance", "status": "skipped_same_config"})
        return
    et = pd.read_csv(events_runner.table_path(cfg, "events"))
    se = pd.read_csv(events_runner.table_path(cfg, "split_events"))
    with xr.open_dataset(ref_path) as r:
        r = r.load()
    times = pd.DatetimeIndex(r.time.values)
    city = events_runner.read_city(cfg, cfg["events"]["resolution"])
    if not np.array_equal(city.time.values, r.time.values):
        raise ValueError("City gauge product and reference series have different time axes")
    runs = inventory.event_runs(times, r.state.values, r.event_id.values.astype(str), r.domain_mean_rain_mmh.values.astype("float64"),
                                r.radar_wet_area_fraction.values.astype("float64"), r.max_city_rate_mmh.values.astype("float64"),
                                city.amount_mm.transpose("time", "id").values.astype("float64"),
                                pd.Timedelta(cfg["events"]["resolution"]), scfg)
    sg = inventory.significance(et, scfg)
    sg = sg.merge(inventory.wet_core(runs), on="event_id", how="left").merge(se[["event_id", "split_event"]], on="event_id", how="left")
    # consistency with events.csv: the run count must equal merged_from (same wet runs, same merging)
    mismatch = et.set_index("event_id").merged_from.ne(sg.set_index("event_id").n_runs)
    if mismatch.any():
        raise ValueError(f"event_runs disagree with events.csv merged_from for {list(mismatch[mismatch].index[:5])}")
    events.validate_unique(sg, ["event_id"], "event_significance.csv")
    extra = {"input_signature": sig, "rule": scfg,
             "radar_domain_total_mm_meaning": "domain-mean radar accumulation over the event: mean floored rate over the "
                                              "720 domain cells, summed over the 5 min bins x 5/60 h"}
    _write(cfg, P["event_significance"], sg, "significance", "significance and rain class per event (additive to events.csv; "
           "reference-only rule, not tuned on CML)", extra)
    _write(cfg, P["event_runs"], runs, "significance", "wet runs of the network reference inside events, before gap merging", extra)
    rec = {"step": "significance", "status": "written", "n_events": len(sg), "n_significant": int(sg.significant.sum()),
           "rain_class": sg.rain_class.value_counts().to_dict(), "n_runs": len(runs), "runtime_s": round(time.perf_counter() - t0, 2)}
    log.info("significance: %d of %d events significant, %s, %d runs", rec["n_significant"], len(sg), rec["rain_class"], len(runs))
    records.append(rec)


# ---------------------------------------------------------------------------
# Stage: inventory (7.3)
# ---------------------------------------------------------------------------
def tsl_steps_wet_dry(cfg: Config, cids: np.ndarray, subs: list[str]) -> tuple[np.ndarray, np.ndarray, float, float]:
    """tsl_step counts inside events and inside dry periods per column (sublink-major), from the 5 min L3 n_tsl_step
    (netCDF4, one variable per daily file). Returns (steps_wet, steps_dry, event_hours, dry_hours)."""
    res = cfg["events"]["resolution"]
    with xr.open_dataset(events_runner.reference_path(cfg)) as r:
        in_ev = r.event_id.values.astype(str) != ""
        times = pd.DatetimeIndex(r.time.values)
    dp = pd.read_csv(events_runner.table_path(cfg, "dry_periods"))
    in_dry = events_runner.in_periods(times, dp)
    parts = []
    for d in cfg.period_days():
        with netCDF4.Dataset(aggregate_runner.cml_path(cfg, res, d)) as nc:
            nc.set_auto_mask(False)
            if not np.array_equal(nc["cml_id"][:].astype("int64"), cids):
                raise ValueError(f"cml_id order differs in {d.date()}")
            parts.append(nc["n_tsl_step"][:].astype("int64"))
    n = np.concatenate(parts).reshape(len(times), -1)
    h = pd.Timedelta(res) / pd.Timedelta("1h")
    return n[in_ev].sum(0), n[in_dry].sum(0), in_ev.sum() * h, in_dry.sum() * h


def dry_period_medians(cfg: Config) -> tuple[pd.DatetimeIndex, np.ndarray, np.ndarray, list[str]]:
    """Median tl_mean of every dry period per column (sublink-major); NaN when < 50 % of its bins are valid.
    Returns (period mid-times, medians (period, column), cml_ids, sublink_ids)."""
    times, cids, subs, d = events_runner.read_core(cfg, cfg["events"]["resolution"], ["tl_mean"])
    tl = d["tl_mean"].reshape(len(times), -1).astype("float64")
    dp = pd.read_csv(events_runner.table_path(cfg, "dry_periods"), parse_dates=["start", "end"])
    idx = pd.Index(times)
    step = times[1] - times[0]
    meds = []
    for p in dp.itertuples():
        a, b = idx.get_loc(p.start), idx.get_loc(p.end - step) + 1
        blk = tl[a:b]
        ok = np.isfinite(blk).mean(0) >= 0.5
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)  # an all-NaN column -> NaN median
            meds.append(np.where(ok, np.nanmedian(blk, 0), np.nan))
    mids = pd.DatetimeIndex(dp.start + (dp.end - dp.start) / 2)
    return mids, np.array(meds), cids, subs


def dry_baseline_range(cfg: Config) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Per column: p90 - p10 of the dry-period medians of tl_mean (NaN if fewer than 5 usable periods)."""
    _, meds, cids, subs = dry_period_medians(cfg)
    n = np.isfinite(meds).sum(0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        rng = np.nanpercentile(meds, 90, axis=0) - np.nanpercentile(meds, 10, axis=0)
    return np.where(n >= 5, rng, np.nan), n, cids, subs


def build_sublinks(cfg: Config, log) -> pd.DataFrame:
    """The sub-link inventory table (metadata, completeness, flags, TSL, noise, status)."""
    T = cfg.tables_dir
    icfg = cfg["inventory"]
    # static metadata (from the L1 coords of one day: identical on every day)
    with xr.open_dataset(cfg.l1_cml_day_path(cfg.period_start)) as l1:
        meta = pd.DataFrame({
            "cml_id": np.tile(l1.cml_id.values, 2), "sublink_id": np.repeat(l1.sublink_id.values.astype(str), l1.sizes["cml_id"]),
            "raw_sublink_id": l1.raw_sublink_id.values.ravel(), "direction": l1.direction.values.astype(str).ravel(),
            "tx_site": l1.tx_site.values.ravel(), "rx_site": l1.rx_site.values.ravel(),
            "frequency_ghz": l1.frequency.values.ravel(), "polarization": l1.polarization.values.astype(str).ravel()})
    meta["sublink"] = meta.sublink_id.str[-1].astype(int)
    lg = pd.read_csv(T / "link_geometry.csv")
    sl = pd.read_csv(T / "split_links.csv")
    base = meta.merge(lg[["cml_id", "length_m", "site_0_id", "site_1_id", "site_0_antenna_diameter_m", "site_1_antenna_diameter_m",
                          "n_radar_cells"]], on="cml_id").merge(sl[["cml_id", "hop", "split_link", "length_class", "frequency_band"]], on="cml_id")
    base["tx_site_id"] = np.where(base.tx_site == 0, base.site_0_id, base.site_1_id)
    base["rx_site_id"] = np.where(base.tx_site == 0, base.site_1_id, base.site_0_id)
    # completeness + flags (qc_daily_summary)
    daily = pd.read_csv(T / "qc_daily_summary.csv")
    comp = inventory.completeness(daily)
    base = base.merge(comp, on=["cml_id", "sublink_id"], how="left")
    # first / last day with any valid tl sample: shows sub-links that stop (or start) reporting within the period
    act = daily[daily.n_valid_tl > 0].groupby(["cml_id", "sublink_id"]).day
    base = base.merge(act.min().rename("first_valid_day").reset_index(), on=["cml_id", "sublink_id"], how="left")                .merge(act.max().rename("last_valid_day").reset_index(), on=["cml_id", "sublink_id"], how="left")
    # quantization + TSL (period stats + histograms)
    ps = pd.read_csv(T / "sublink_period_stats.csv")
    z = np.load(T / "sublink_period_hist.npz")
    q = inventory.quantization(z["rsl"], z["tsl"], float(z["bin"]))
    q["cml_id"], q["sublink_id"] = z["cml_id"], z["sublink_id"].astype(str)
    base = base.merge(q, on=["cml_id", "sublink_id"]).merge(ps[["cml_id", "sublink_id", "n_valid_tsl", "tsl_n_steps"]],
                                                         on=["cml_id", "sublink_id"])
    base["tsl_step_fraction"] = base.tsl_n_steps / base.n_valid_tsl.replace(0, np.nan)
    base["tsl_steps_per_day"] = base.tsl_n_steps / (base.n_valid_tsl.replace(0, np.nan) / 8640)
    base["tsl_behaviour"] = inventory.tsl_behaviour(base.n_valid_tsl, base.tsl_n_steps, base.tsl_top_share, icfg["tsl"])
    # descriptive: do TSL steps follow rain (ATPC-like)?
    with netCDF4.Dataset(aggregate_runner.cml_path(cfg, "5min", cfg.period_start)) as nc:
        cids = nc["cml_id"][:].astype("int64")
        subs = [str(s) for s in nc["sublink_id"][:]]
    sw, sd, eh, dh = tsl_steps_wet_dry(cfg, cids, subs)
    st = pd.DataFrame({"cml_id": np.tile(cids, len(subs)), "sublink_id": np.repeat(subs, len(cids)),
                       "tsl_steps_in_events": sw, "tsl_steps_in_dry": sd})
    st["tsl_step_rate_events_per_h"] = st.tsl_steps_in_events / eh
    st["tsl_step_rate_dry_per_h"] = st.tsl_steps_in_dry / dh
    tc = icfg["tsl"]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = st.tsl_step_rate_events_per_h / st.tsl_step_rate_dry_per_h
    st["tsl_follows_fades_descriptive"] = (st.tsl_steps_in_events >= tc["follows_min_steps"]) & \
                                          ((ratio >= tc["follows_rate_ratio"]) | (st.tsl_steps_in_dry == 0))
    base = base.merge(st, on=["cml_id", "sublink_id"])
    # dry-weather noise and baseline stability
    ds = pd.read_csv(T / "sublink_dry_stats.csv")
    base = base.merge(ds[["cml_id", "sublink_id", "tl_median", "tl_std", "tl_scaled_mad", "tl_intra_bin_std_median", "rsl_std",
                          "rsl_intra_bin_std_median", "tl_diurnal_range", "tl_daily_median_drift_std", "n_dry_days"]]
                      .rename(columns={"tl_median": "dry_tl_median", "tl_std": "dry_tl_std", "tl_scaled_mad": "dry_tl_scaled_mad",
                                       "tl_intra_bin_std_median": "dry_tl_10s_noise_median", "rsl_std": "dry_rsl_std",
                                       "rsl_intra_bin_std_median": "dry_rsl_10s_noise_median",
                                       "tl_diurnal_range": "dry_tl_diurnal_range", "tl_daily_median_drift_std": "dry_tl_daily_drift_std"}),
                      on=["cml_id", "sublink_id"])
    rng, nper, c2, s2 = dry_baseline_range(cfg)
    br = pd.DataFrame({"cml_id": np.tile(c2, len(s2)), "sublink_id": np.repeat(s2, len(c2)),
                       "dry_baseline_range_db": rng, "n_dry_periods_used": nper})
    base = base.merge(br, on=["cml_id", "sublink_id"])
    stt = inventory.quality_status(base, icfg["status"])
    base["status"], base["status_reason"] = stt.status, stt.reason
    old = pd.read_csv(T / "sublink_qc_status.csv")[["cml_id", "sublink_id", "status"]].rename(columns={"status": "status_step3"})
    base = base.merge(old, on=["cml_id", "sublink_id"], how="left")
    base.attrs["thresholds"] = stt.attrs["thresholds"]
    return base.sort_values(["cml_id", "sublink_id"]).reset_index(drop=True)


def build_links(cfg: Config, sub: pd.DataFrame) -> pd.DataFrame:
    """The link inventory table (sub-links summarised, gauges, neighbours, radar cells)."""
    T = cfg.tables_dir
    radii = cfg["inventory"]["gauge_radii_m"]
    g = sub.groupby("cml_id")
    lk = pd.DataFrame({
        "length_m": g.length_m.first(), "frequency_ghz_sublink_1": g.apply(lambda d: d[d.sublink == 1].frequency_ghz.iloc[0], include_groups=False),
        "frequency_ghz_sublink_2": g.apply(lambda d: d[d.sublink == 2].frequency_ghz.iloc[0], include_groups=False),
        "polarization": g.polarization.first(), "length_class": g.length_class.first(), "frequency_band": g.frequency_band.first(),
        "split_link": g.split_link.first(), "hop": g.hop.first(), "site_0_id": g.site_0_id.first(), "site_1_id": g.site_1_id.first(),
        "n_cells": g.n_radar_cells.first(),
        "completeness_tl_min": g.completeness_tl.min(), "completeness_tl_max": g.completeness_tl.max(),
        "share_not_blocked_min": g.share_not_blocked.min(),
        "status_sublink_1": g.apply(lambda d: d[d.sublink == 1].status.iloc[0], include_groups=False),
        "status_sublink_2": g.apply(lambda d: d[d.sublink == 2].status.iloc[0], include_groups=False),
        "tsl_behaviour_sublink_1": g.apply(lambda d: d[d.sublink == 1].tsl_behaviour.iloc[0], include_groups=False),
        "tsl_behaviour_sublink_2": g.apply(lambda d: d[d.sublink == 2].tsl_behaviour.iloc[0], include_groups=False),
        "dry_tl_std_max": g.dry_tl_std.max(),
        "last_valid_day": g.last_valid_day.max(),   # last day on which ANY sub-link had valid tl
    }).reset_index()
    order = {"good": 0, "degraded": 1, "unusable": 2}
    lk["status_best"] = [min((a, b), key=order.get) for a, b in zip(lk.status_sublink_1, lk.status_sublink_2)]
    m = pd.read_csv(T / "link_gauge_matches.csv")
    for r in radii:
        lk[f"n_gauges_within_{int(r)}m"] = lk.cml_id.map(m[m.distance_to_path_m <= r].groupby("cml_id").size()).fillna(0).astype(int)
    # Nearest gauge of each network REGARDLESS of the matching radius (link_gauge_matches stops at 5 km): distance from
    # every gauge to the path segment in EPSG:3006, the same point-to-segment rule as the geometry step.
    lg = pd.read_csv(T / "link_geometry.csv").set_index("cml_id").loc[lk.cml_id]
    gauges = pd.read_csv(T / "gauges.csv")
    d, _ = geometry.point_to_path(gauges.x.to_numpy()[:, None], gauges.y.to_numpy()[:, None],
                                  lg.site_0_x_3006_m.to_numpy()[None], lg.site_0_y_3006_m.to_numpy()[None],
                                  lg.site_1_x_3006_m.to_numpy()[None], lg.site_1_y_3006_m.to_numpy()[None])  # (gauge, link)
    ids = gauges.id.to_numpy()
    k = np.argmin(d, axis=0)
    lk["nearest_gauge_id"] = ids[k]
    lk["nearest_gauge_network"] = gauges.network.to_numpy()[k]
    lk["nearest_gauge_distance_m"] = d[k, np.arange(d.shape[1])]
    for net in ("city", "smhi"):
        sel = (gauges.network == net).to_numpy()
        kk = np.argmin(d[sel], axis=0)
        lk[f"nearest_{net}_gauge_id"] = ids[sel][kk]
        lk[f"nearest_{net}_gauge_distance_m"] = d[sel][kk, np.arange(d.shape[1])]
    lk["has_gauge_within_1km"] = lk.n_gauges_within_1000m > 0
    nb = pd.read_csv(T / "link_neighbours.csv")
    both = pd.concat([nb[["cml_id_a", "min_path_distance_m"]].rename(columns={"cml_id_a": "c"}),
                      nb[["cml_id_b", "min_path_distance_m"]].rename(columns={"cml_id_b": "c"})])
    lk["n_neighbours_within_1km"] = lk.cml_id.map(both[both.min_path_distance_m <= 1000].groupby("c").size()).fillna(0).astype(int)
    lk["n_neighbours_within_5km"] = lk.cml_id.map(both.groupby("c").size()).fillna(0).astype(int)
    lk["single_radar_cell"] = lk.n_cells == 1
    return lk


def build_gauges(cfg: Config) -> pd.DataFrame:
    """The gauge inventory table (completeness, totals, radar pixel comparison)."""
    T = cfg.tables_dir
    g = pd.read_csv(T / "gauges.csv")
    rows = []
    for net, f in (("city", "city_1min.nc"), ("smhi", "smhi_15min.nc")):
        with xr.open_dataset(cfg.l1_gauges_dir / f) as ds:
            a = ds.rainfall_amount.load()
        for gid in a.id.values:
            s = a.sel(id=gid).to_series()
            q = float(g.loc[g.id == gid, "quantization_mm"].iloc[0])
            row = {"id": str(gid), "jja_total_mm": float(s.astype("float64").sum()),
                   "wet_fraction": float((s >= q - 1e-6).sum() / s.notna().sum())}
            for mo, name in ((6, "jun"), (7, "jul"), (8, "aug")):
                row[f"completeness_{name}"] = float(s[s.index.month == mo].notna().mean())
            rows.append(row)
    out = g.merge(pd.DataFrame(rows), on="id")
    pix = pd.read_csv(T / "gauge_radar_pixel.csv").rename(columns={"gauge_id": "id"})
    out = out.merge(pix[["id", "iy", "ix", "distance_to_cell_centre_m"]], on="id", how="left")
    with xr.open_dataset(aggregate_runner.radar_path(cfg, "at_gauges", "5min")) as r:
        tot = (r.R_pixel.astype("float64").sum("time") * 5 / 60).to_series()
    out["radar_pixel_jja_total_mm"] = out.id.map(tot)
    out["radar_over_gauge_ratio"] = out.radar_pixel_jja_total_mm / out.jja_total_mm
    return out


def run_inventory(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """Stage inventory: sublinks, links, gauges and quality_status tables."""
    t0 = time.perf_counter()
    P = out_paths(cfg)
    T = cfg.tables_dir
    inputs = [T / f for f in ("qc_daily_summary.csv", "sublink_period_stats.csv", "sublink_dry_stats.csv", "split_links.csv",
                              "link_geometry.csv", "link_gauge_matches.csv", "link_neighbours.csv", "gauge_radar_pixel.csv",
                              "gauges.csv", "sublink_qc_status.csv", "dry_periods.csv", "events.csv")]
    sig = _signature(inputs, [events_runner.core_path(cfg, "5min"), events_runner.reference_path(cfg),
                              aggregate_runner.radar_path(cfg, "at_gauges", "5min")])
    targets = [P[k] for k in ("sublinks", "links", "gauges", "quality_status")]
    if not force and _up_to_date(cfg, "inventory", targets, sig):
        log.info("inventory: up to date, skipped (use --force to rebuild)")
        records.append({"step": "inventory", "status": "skipped_same_config"})
        return
    sub = build_sublinks(cfg, log)
    thr = sub.attrs["thresholds"]
    lk = build_links(cfg, sub)
    gg = build_gauges(cfg)
    qs = sub[["cml_id", "sublink_id", "raw_sublink_id", "status", "status_reason", "status_step3", "completeness_tl",
              "completeness_tl_jun", "completeness_tl_jul", "completeness_tl_aug", "share_not_blocked", "flag_share_blocking_present",
              "flag_share_rx_floor",
              "flag_share_spike", "flag_share_flatline", "flag_share_tsl_step", "dry_tl_std", "dry_baseline_range_db"]]
    for df, what in ((sub, "sublinks"), (lk, "links"), (gg, "gauges"), (qs, "quality_status")):
        events.validate_unique(df, ["id"] if what == "gauges" else (["cml_id"] if what == "links" else ["cml_id", "sublink_id"]), what)
    extra = {"input_signature": sig, "status_thresholds_derived": thr, "status_rules": cfg["inventory"]["status"]}
    _write(cfg, P["sublinks"], sub, "inventory", "one row per sub-link: metadata, completeness, flags, quantization, TSL "
           "behaviour, dry-weather noise, status (tsl_follows_fades_descriptive is DESCRIPTIVE ONLY)", extra)
    _write(cfg, P["links"], lk, "inventory", "one row per link: the sub-link inventory summarised + gauges, neighbours, radar cells", extra)
    _write(cfg, P["gauges"], gg, "inventory", "gauges.csv plus monthly completeness, JJA total, wet fraction, radar pixel total", extra)
    _write(cfg, P["quality_status"], qs, "inventory", "step 7 quality status per sub-link (good / degraded / unusable + reason); "
           "the step 3 sublink_qc_status.csv is unchanged", extra)
    rec = {"step": "inventory", "status": "written", "n_sublinks": len(sub), "n_links": len(lk),
           "status_counts": sub.status.value_counts().to_dict(), "tsl_behaviour": sub.tsl_behaviour.value_counts().to_dict(),
           "thresholds": thr, "runtime_s": round(time.perf_counter() - t0, 2)}
    log.info("inventory: %d sub-links, %d links, status %s, tsl %s (%.1f s)", len(sub), len(lk), rec["status_counts"],
             rec["tsl_behaviour"], rec["runtime_s"])
    records.append(rec)


# ---------------------------------------------------------------------------
# Stage: event_links_extra (7b; does NOT edit event_links.csv)
# ---------------------------------------------------------------------------
def read_n_tsl_step(cfg: Config) -> tuple[np.ndarray, np.ndarray]:
    """n_tsl_step (time, sublink, cml) of the 5 min L3 files, whole period (netCDF4, one variable per daily file)."""
    parts, cids = [], None
    for d in cfg.period_days():
        with netCDF4.Dataset(aggregate_runner.cml_path(cfg, "5min", d)) as nc:
            nc.set_auto_mask(False)
            c = nc["cml_id"][:].astype("int64")
            if cids is not None and not np.array_equal(c, cids):
                raise ValueError(f"cml_id order differs in {d.date()}")
            cids = c
            parts.append(nc["n_tsl_step"][:].astype("int64"))
    return np.concatenate(parts), cids


def run_event_links_extra(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """tables/event_links_flags_extra.csv: per (event, link) the share of 5 min bins inside the event whose flags_or
    holds flatline / tsl_step, per sub-link and the max, and the number of 10 s tsl_step samples in the event.
    Information only (a link frozen during rain shows a high flatline share)."""
    t0 = time.perf_counter()
    out = cfg.tables_dir / "event_links_flags_extra.csv"
    sig = _signature([events_runner.table_path(cfg, "events"), events_runner.table_path(cfg, "split_links"),
                      events_runner.table_path(cfg, "split_events")], [events_runner.core_path(cfg, "5min")])
    if not force and _up_to_date(cfg, "event_links_extra", [out], sig):
        log.info("event_links_extra: up to date, skipped (use --force to rebuild)")
        records.append({"step": "event_links_extra", "status": "skipped_same_config"})
        return
    from cmlpipe.io import schema
    times, cids, subs, core = events_runner.read_core(cfg, "5min", ["flags_or"])
    nts, c2 = read_n_tsl_step(cfg)
    if not np.array_equal(c2, cids):
        raise ValueError("daily L3 files and the CML core list the links in different orders")
    fl = core["flags_or"]
    bits = {f: np.uint16(1 << schema.FLAGS[f]) for f in ("flatline", "tsl_step")}
    et = pd.read_csv(events_runner.table_path(cfg, "events"), parse_dates=["start", "end"])
    sl = pd.read_csv(events_runner.table_path(cfg, "split_links")).set_index("cml_id").split_link
    idx = pd.Index(times)
    step = times[1] - times[0]
    blocks = []
    for e in et.itertuples():
        a, b = idx.get_loc(e.start), idx.get_loc(e.end - step) + 1
        df = pd.DataFrame({"event_id": e.event_id, "cml_id": cids})
        for f, bit in bits.items():
            frac = ((fl[a:b] & bit) != 0).mean(axis=0)  # (sublink, cml)
            for k, s in enumerate(subs):
                df[f"frac_{f}_event_{s}"] = frac[k]
            df[f"frac_{f}_event_max"] = frac.max(axis=0)
        n = nts[a:b].sum(axis=0)
        for k, s in enumerate(subs):
            df[f"n_tsl_step_in_event_{s}"] = n[k]
        df["n_tsl_step_in_event"] = n.sum(axis=0)
        df["split_event"] = e.split_event
        blocks.append(df)
    res = pd.concat(blocks, ignore_index=True)
    res["split_link"] = res.cml_id.map(sl)
    events.validate_unique(res, ["event_id", "cml_id"], out.name)
    _write(cfg, out, res, "event_links_extra", "per (event, link): share of 5 min bins with flatline / tsl_step in the event "
           "and 10 s tsl_step counts (information only; event_links.csv is unchanged)", {"input_signature": sig})
    rec = {"step": "event_links_extra", "status": "written", "n_rows": len(res), "runtime_s": round(time.perf_counter() - t0, 2)}
    log.info("event_links_extra: %d rows (%.1f s)", len(res), rec["runtime_s"])
    records.append(rec)
