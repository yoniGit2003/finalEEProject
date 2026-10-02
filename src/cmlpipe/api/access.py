"""Read one link and everything that belongs to it, on ONE time axis, for a window you choose.

Usage (from the repository root, after the aggregate stages have run):

    from cmlpipe.config import load_config
    from cmlpipe.api.access import load_link_window

    cfg = load_config("configs/openmrg.yaml")
    ds = load_link_window(cfg, 10227, "5min", "2015-07-28 14:30", "2015-07-28 17:00",
                          pad_before="1h", pad_after="2h", max_gauge_dist_m=1000)
    ds.tl_mean.sel(sublink_id="sublink_1").plot()      # CML statistics, dims (time, sublink_id)
    ds.radar_R_mean.plot()                               # radar along the path, dims (time,)
    ds.gauge_rate_mm_h.plot.line(x="time")               # matched gauges, dims (time, gauge)
    print(ds.attrs["sources_absent"])                    # e.g. radar at 1 min: never interpolated

The window is [start - pad_before, end + pad_after), cut to whole bins of the resolution (window_edges).
Every source is opened lazily and sliced BEFORE it is loaded. A source that does not exist at the
requested resolution (radar at 1 min, the SMHI gauge below 15 min) is left out and named in the attributes.
CML days whose L3 file has not been produced yet are NaN and listed in attrs["cml_days_missing"].
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from cmlpipe.aggregate import stage as aggregate_runner
from cmlpipe.config import Config


def window_edges(start, end, pad_before="0h", pad_after="0h", res: str | None = None
                 ) -> tuple[pd.Timestamp, pd.Timestamp]:
    """[t0, t1) = [start - pad_before, end + pad_after), widened to whole bins of `res` if given.

    t0 is floored and t1 ceiled to the bin grid, so the window never cuts a bin and always contains the
    requested span. Raises if the window is empty or a padding is negative.
    """
    pb, pa = pd.Timedelta(pad_before), pd.Timedelta(pad_after)
    if pb < pd.Timedelta(0) or pa < pd.Timedelta(0):
        raise ValueError("pad_before and pad_after must be >= 0")
    t0, t1 = pd.Timestamp(start) - pb, pd.Timestamp(end) + pa
    if res is not None:
        t0, t1 = t0.floor(res), t1.ceil(res)
    if t1 <= t0:
        raise ValueError(f"empty window: [{t0}, {t1})")
    return t0, t1


def _sources_at(cfg: Config, res: str) -> dict:
    a = cfg["aggregate"]
    return {"cml": res in a["resolutions"], "radar": res in a["radar"]["resolutions"],
            "city": res in a["gauges"]["city_resolutions"], "smhi": res in a["gauges"]["smhi_resolutions"]}


def _slice(ds: xr.Dataset, t0, t1, res) -> xr.Dataset:
    """Time slice [t0, t1) (sel is inclusive at both ends, so stop one bin early)."""
    return ds.sel(time=slice(t0, t1 - pd.Timedelta(res)))


def load_link_window(cfg: Config, cml_id: int, res: str, start, end, pad_before="0h", pad_after="0h",
                     radar: bool = True, gauges: bool = True, max_gauge_dist_m: float = 2000) -> xr.Dataset:
    """The link's CML statistics, radar along its path and matched gauges on one time axis (see module doc).

    Variables: the CML statistics with their L3 names (tl_mean, rsl_median, n_valid_tl, flags_or, ...), dims
    (time, sublink_id); radar_<name> (radar_R_mean, ...), dims (time,); gauge_<name> (gauge_amount_mm,
    gauge_rate_mm_h, ...), dims (time, gauge) with coords network, distance_to_path_m, position_along_path.
    Static link metadata are attributes (length_m, azimuth_deg, polarization, n_radar_cells, radar_coverage)
    or sublink coords (frequency, qc_status, qc_reason).
    """
    acfg = cfg["aggregate"]
    if res not in acfg["resolutions"]:
        raise ValueError(f"unknown resolution {res!r}; available: {acfg['resolutions']}")
    lg = pd.read_csv(cfg.tables_dir / "link_geometry.csv").set_index("cml_id")
    if cml_id not in lg.index:
        raise KeyError(f"unknown cml_id {cml_id}")
    t0, t1 = window_edges(start, end, pad_before, pad_after, res)
    if t0 < cfg.period_start or t1 > cfg.period_end:
        raise ValueError(f"window [{t0}, {t1}) is outside the period [{cfg.period_start}, {cfg.period_end})")
    step = pd.Timedelta(res)
    times = pd.date_range(t0, t1 - step, freq=step)
    have = _sources_at(cfg, res)
    absent = []

    # --- CML: one file per day ---
    parts, missing_days = [], []
    for day in pd.date_range(t0.normalize(), (t1 - step).normalize(), freq="D"):
        p = aggregate_runner.cml_path(cfg, res, day)
        if not p.exists():
            missing_days.append(f"{day:%Y-%m-%d}")
            continue
        with xr.open_dataset(p) as f:
            parts.append(_slice(f.sel(cml_id=cml_id), t0, t1, res).load())
    if parts:
        cml = xr.concat(parts, "time") if len(parts) > 1 else parts[0]
        cml = cml.reindex(time=times)  # days not produced yet -> NaN; never interpolated
    else:
        raise FileNotFoundError(f"no L3 CML file at {res} for any day of [{t0}, {t1}): run aggregate_cml first")
    ds = cml.drop_vars([c for c in cml.coords if c not in ("time", "sublink_id") and cml[c].dims != ("sublink_id",)])
    ds = ds.drop_vars(["cml_id"], errors="ignore")

    # --- static link metadata ---
    row = lg.loc[cml_id]
    st = pd.read_csv(cfg.tables_dir / "sublink_qc_status.csv")
    st = st[st.cml_id == cml_id].set_index("sublink_id").reindex(ds.sublink_id.values)
    ds = ds.assign_coords(qc_status=("sublink_id", st["status"].fillna("").to_numpy(str)),
                          qc_reason=("sublink_id", st["reason"].fillna("").to_numpy(str)))
    ds.attrs = {
        "cml_id": int(cml_id), "interval": res, "time_convention": "interval_start",
        "window_start": str(t0), "window_end_exclusive": str(t1),
        "length_m": float(row.length_m), "azimuth_deg": float(row.azimuth_deg), "polarization": str(row.polarization),
        "n_radar_cells": int(row.n_radar_cells), "radar_coverage": float(row.radar_coverage),
        "frequency_ghz_sublink_1": float(row.frequency_ghz_sublink_1),
        "frequency_ghz_sublink_2": float(row.frequency_ghz_sublink_2),
        "cml_days_missing": " ".join(missing_days),
    }

    # --- radar along the path ---
    if radar:
        if have["radar"]:
            with xr.open_dataset(aggregate_runner.radar_path(cfg, "along_links", res)) as f:
                r = _slice(f.sel(cml_id=cml_id), t0, t1, res).load().reindex(time=times)
            for name, v in r.data_vars.items():
                ds[f"radar_{name}"] = v.drop_vars([c for c in v.coords if c != "time"])
            ds.attrs["radar_floor_min_dbz"] = f.attrs.get("min_dbz_floor")
        else:
            absent.append(f"radar: native 5 min, not available at {res} (never interpolated)")

    # --- gauges within max_gauge_dist_m of the path ---
    if gauges:
        m = pd.read_csv(cfg.tables_dir / "link_gauge_matches.csv")
        m = m[(m.cml_id == cml_id) & (m.distance_to_path_m <= max_gauge_dist_m)].sort_values("distance_to_path_m")
        blocks = []
        for net in ("city", "smhi"):
            ids = m[m.network == net].gauge_id.tolist()
            if not ids:
                continue
            if not have[net]:
                absent.append(f"{net} gauges {ids}: not available at {res} (never disaggregated)")
                continue
            with xr.open_dataset(aggregate_runner.gauge_path(cfg, net, res)) as f:
                g = _slice(f.sel(id=ids), t0, t1, res).load().reindex(time=times)
            g = g.drop_vars([c for c in g.coords if c not in ("time", "id")])
            blocks.append(g.rename({"id": "gauge"}))
        if blocks:
            # City has peak_1min_rate_mm_h, SMHI has qc_worst: give each block the other's variables as NaN
            names = sorted({n for b in blocks for n in b.data_vars})
            for k, b in enumerate(blocks):
                for n in names:
                    if n not in b:
                        b[n] = xr.full_like(b[next(iter(b.data_vars))], np.nan, dtype="float32")
                blocks[k] = b
            g = xr.concat(blocks, "gauge")
            mm = m.set_index("gauge_id").loc[g.gauge.values]
            g = g.assign_coords(network=("gauge", mm.network.to_numpy(str)),
                                distance_to_path_m=("gauge", mm.distance_to_path_m.to_numpy(float)),
                                position_along_path=("gauge", mm.position_along_path.to_numpy(float)))
            for name, v in g.data_vars.items():
                ds[f"gauge_{name}"] = v
        ds.attrs["max_gauge_dist_m"] = float(max_gauge_dist_m)
    ds.attrs["sources_absent"] = "; ".join(absent) if absent else "none"
    return ds
