"""Interactive access helpers: find links / events, and plot one link or one event.

    from cmlpipe.config import load_config
    from cmlpipe.api import explore
    cfg = load_config("configs/openmrg.yaml")

    explore.find_links(cfg, length_km=(1, 3), frequency_band=">=25GHz", status="good", split="test", has_gauge_within_m=1000)
    explore.find_events(cfg, significant=True, month=7, split="test")
    ds, fig = explore.show_link(cfg, 10227, "2015-07-28 12:00", "2015-07-28 20:00", res="5min")
    fig = explore.show_event(cfg, "E20150727_2140", cml_id=10227)   # one link; without cml_id: the network overview

WARNING: do not open the 92 daily L3 files (or L1 / L2 files) with xarray.open_mfdataset or in a loop with
xarray: on this Windows machine xarray inspects every variable of every file and it is very slow, and dask is not
installed. These helpers read the consolidated files (L3_products/cml/consolidated/cml_<res>_core.nc) or, for a
time window, only the daily files that cover it (cmlpipe.access.load_link_window), sliced before loading.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cmlpipe.api import access
from cmlpipe.events import stage as events_runner
from cmlpipe.inventory import stage as inventory_runner
from cmlpipe.viz import quicklooks
from cmlpipe.io import schema
from cmlpipe.config import Config
from cmlpipe.viz.style import C_AQUA, C_ORANGE, C_SURF, SUBLINK_C, label, style


def _links(cfg: Config) -> pd.DataFrame:
    return pd.read_csv(inventory_runner.out_paths(cfg)["links"])


def find_links(cfg: Config, length_km: tuple | None = None, frequency_band: str | None = None, status: str | None = None,
               split: str | None = None, has_gauge_within_m: float | None = None, tsl_behaviour: str | None = None,
               links: pd.DataFrame | None = None) -> pd.DataFrame:
    """Rows of tables/inventory/links.csv matching ALL given filters.

    length_km: (min, max) inclusive; frequency_band: "<25GHz" / ">=25GHz"; status: "good" / "degraded" / "unusable"
    (matches the BEST of the two sub-links, status_best); split: train / test / excluded; has_gauge_within_m: 250, 500
    or 1000 (a gauge within that distance of the path); tsl_behaviour: matches either sub-link.
    links: pass a table to filter instead of reading the file (used by the tests).
    """
    lk = _links(cfg) if links is None else links
    m = np.ones(len(lk), bool)
    if length_km is not None:
        m &= (lk.length_m >= 1000 * length_km[0]).to_numpy() & (lk.length_m <= 1000 * length_km[1]).to_numpy()
    if frequency_band is not None:
        m &= (lk.frequency_band == frequency_band).to_numpy()
    if status is not None:
        m &= (lk.status_best == status).to_numpy()
    if split is not None:
        m &= (lk.split_link == split).to_numpy()
    if has_gauge_within_m is not None:
        col = f"n_gauges_within_{int(has_gauge_within_m)}m"
        if col not in lk:
            raise ValueError(f"has_gauge_within_m must be one of the inventory radii: "
                             f"{[c.split('_')[-1] for c in lk.columns if c.startswith('n_gauges_within_')]}")
        m &= (lk[col] > 0).to_numpy()
    if tsl_behaviour is not None:
        m &= ((lk.tsl_behaviour_sublink_1 == tsl_behaviour) | (lk.tsl_behaviour_sublink_2 == tsl_behaviour)).to_numpy()
    return lk[m].reset_index(drop=True)


def find_events(cfg: Config, significant: bool | None = None, month: int | None = None, split: str | None = None,
                intensity: str | None = None, rain_class: str | None = None) -> pd.DataFrame:
    """Rows of events.csv joined with event_significance.csv, matching ALL given filters."""
    et = pd.read_csv(events_runner.table_path(cfg, "events"), parse_dates=["start", "end"])
    sg = pd.read_csv(inventory_runner.out_paths(cfg)["event_significance"])
    et = et.merge(sg[["event_id", "significant", "rain_class", "n_runs", "wet_core_start", "wet_core_end"]], on="event_id")
    m = np.ones(len(et), bool)
    if significant is not None:
        m &= (et.significant == significant).to_numpy()
    if month is not None:
        m &= (et.month == month).to_numpy()
    if split is not None:
        m &= (et.split_event == split).to_numpy()
    if intensity is not None:
        m &= (et.intensity_class == intensity).to_numpy()
    if rain_class is not None:
        m &= (et.rain_class == rain_class).to_numpy()
    return et[m].reset_index(drop=True)


def show_link(cfg: Config, cml_id: int, start, end, res: str = "5min", pad: str = "0h"):
    """(dataset, figure): tl_mean of both sub-links, radar along the path and the matched gauges in [start, end)."""
    ds = access.load_link_window(cfg, cml_id, res, start, end, pad_before=pad, pad_after=pad, max_gauge_dist_m=1000)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(11, 6), sharex=True, facecolor=C_SURF, layout="constrained")
    for a in (a1, a2):
        style(a)
    for s in schema.SUBLINK_IDS:
        a1.plot(ds.time, ds.tl_mean.sel(sublink_id=s), color=SUBLINK_C[s], lw=1.1, label=s)
    label(a1, "", "TL mean (dB)", f"link {cml_id}, {res}")
    a1.legend(fontsize=8, frameon=False)
    if "radar_R_mean" in ds:
        a2.step(ds.time, ds.radar_R_mean, where="post", color=C_AQUA, label="radar along path")
    if "gauge_rate_mm_h" in ds:
        for g in ds.gauge.values:
            a2.step(ds.time, ds.gauge_rate_mm_h.sel(gauge=g), where="post", color=C_ORANGE, lw=1,
                    label=f"gauge {g} ({float(ds.distance_to_path_m.sel(gauge=g)):.0f} m)")
    label(a2, "time (UTC)", "rain (mm/h)", ds.attrs.get("sources_absent", ""))
    a2.legend(fontsize=8, frameon=False)
    return ds, fig


def show_event(cfg: Config, event_id: str, cml_id: int | None = None, out_path=None):
    """The quicklook of (event, link), or the network overview of the event when cml_id is None. Returns the PNG path."""
    out = Path(out_path) if out_path else Path(tempfile.mkdtemp()) / f"{event_id}_{cml_id or 'overview'}.png"
    if cml_id is not None:
        return quicklooks.make_event_quicklook(cfg, event_id, cml_id, out)
    ev = pd.read_csv(events_runner.table_path(cfg, "events"), parse_dates=["start", "end"]).set_index("event_id").loc[event_id]
    pad = pd.Timedelta(hours=cfg["quicklooks"]["pad_h"])
    w = quicklooks.load_window(cfg, ev.start - pad, ev.end + pad)
    lengths = _links(cfg).set_index("cml_id").loc[w["cml_id"], "length_m"].to_numpy()
    return quicklooks.make_event_overview(cfg, event_id, w, ev.start, ev.end, lengths, out)
