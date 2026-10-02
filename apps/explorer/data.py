"""Data access for the explorer: everything is read through cmlpipe.api.CMLData (no file access here).

Plain functions (no Streamlit), so smoke.py and notebooks can call them; app.py wraps them in st.cache_data.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd

from cmlpipe.api import CMLData

RESOLUTIONS = ["10s", "1min", "5min", "15min", "1h"]
MAX_DAYS_10S = 3
GAUGE_RES = {"city": ["1min", "5min", "15min", "1h"], "smhi": ["15min", "1h"]}
RADAR_RES = ["5min", "15min", "1h"]


def open_release(path) -> CMLData:
    return CMLData.open(path)


def links_table(d: CMLData) -> pd.DataFrame:
    """links + site coordinates (lon_0 / lat_0 / lon_1 / lat_1) + polarization, for pickers and the map."""
    links = d.links()
    s = d.table("sites").set_index("site_id")
    for k in (0, 1):
        links[f"lon_{k}"] = links[f"site_{k}_id"].map(s.lon)
        links[f"lat_{k}"] = links[f"site_{k}_id"].map(s.lat)
    return links


def filter_links(links: pd.DataFrame, length_km=None, band=None, pol=None, status=None, split=None,
                 gauge_within_m=None, min_completeness=0.0) -> pd.DataFrame:
    m = pd.Series(True, index=links.index)
    if length_km:
        m &= links.length_m.between(1000 * length_km[0], 1000 * length_km[1])
    if band:
        m &= links.frequency_band.isin(band)
    if pol:
        m &= links.polarization.isin(pol)
    if status:
        m &= links.status_best.isin(status)
    if split:
        m &= links.split_link.isin(split)
    if gauge_within_m:
        m &= links.nearest_gauge_distance_m <= gauge_within_m
    if min_completeness:
        m &= links.completeness_tl_max.fillna(0) >= min_completeness
    return links[m]


def gauges_for(d: CMLData, cml_id: int, mode: str, radius_m: float) -> list[tuple[str, str, float]]:
    """(gauge_id, network, distance) for 'nearest City', 'nearest SMHI', 'all within radius' or 'none'."""
    m = d.table("link_gauge_matches")
    m = m[m.cml_id == cml_id].sort_values("distance_to_path_m")
    if mode == "nearest City":
        m = m[m.network == "city"].head(1)
    elif mode == "nearest SMHI":
        m = m[m.network == "smhi"].head(1)
    elif mode == "nearest City + SMHI":
        m = pd.concat([m[m.network == "city"].head(1), m[m.network == "smhi"].head(1)])
    elif mode == "all within radius":
        m = m[m.distance_to_path_m <= radius_m]
    else:
        m = m.iloc[0:0]
    return [(r.gauge_id, r.network, float(r.distance_to_path_m)) for r in m.itertuples()]


def _step_onto(s: pd.Series, index: pd.DatetimeIndex, native: str) -> pd.Series:
    """A coarser series shown on a finer grid as steps: each fine bin takes the value of the coarse bin holding it."""
    return s.reindex(index.floor(native)).set_axis(index)


def load_view(d: CMLData, cml_id: int, t0, t1, res: str, qc: str = "mask", gauge_mode: str = "nearest City",
              gauge_radius_m: float = 1000) -> tuple[pd.DataFrame, dict]:
    """The aligned frame of one link + everything the link view needs. Returns (frame, info)."""
    tic = time.perf_counter()
    t0, t1 = pd.Timestamp(t0), pd.Timestamp(t1)
    p = d.catalog["period"]
    t0, t1 = max(t0, pd.Timestamp(p["start"])), min(t1, pd.Timestamp(p["end_exclusive"]))
    if res == "10s" and (t1 - t0) > pd.Timedelta(days=MAX_DAYS_10S):
        raise ValueError(f"10 s is limited to {MAX_DAYS_10S} days in the explorer; choose 1 min or a shorter window")
    notes = []
    frame = d.frame(cml_id, t0, t1, res, include=("cml", "flags", "radar", "met", "state"), qc=qc)
    notes += [n for n in frame.attrs.get("notes", []) if not n.startswith("radar")]
    idx = frame.index
    # radar at 10 s / 1 min: the native 5 min values shown as steps (never interpolated)
    if res not in RADAR_RES:
        r = d.radar_along(cml_id, t0.floor("5min"), min(t1.ceil("5min"), pd.Timestamp(p["end_exclusive"])), "5min")
        r = r.sel(cml_id=cml_id)
        for v in r.data_vars:
            frame["radar_" + v] = _step_onto(r[v].to_series(), idx, "5min").to_numpy()
        notes.append("radar: native 5 min values drawn as steps")
    # gauges: rate on the finest native resolution >= res, plus cumulative amount from the window start
    used = []
    for gid, net, dist in gauges_for(d, cml_id, gauge_mode, gauge_radius_m):
        avail = GAUGE_RES[net]
        g_res = res if res in avail else next((r for r in avail if pd.Timedelta(r) >= pd.Timedelta(res)), avail[0])
        g = d.gauge(gid, t0.floor(g_res), min(t1.ceil(g_res), pd.Timestamp(p["end_exclusive"])), g_res).sel(gauge_id=gid)
        rate, amount = g["rate_mm_h"].to_series(), g["amount_mm"].to_series()
        if g_res != res:
            notes.append(f"gauge {gid}: native {g_res} drawn as steps")
            frame[f"gauge_{gid}_rate"] = _step_onto(rate, idx, g_res).to_numpy()
            frame[f"gauge_{gid}_cum"] = _step_onto(amount.cumsum(), idx, g_res).to_numpy()
        else:
            frame[f"gauge_{gid}_rate"] = rate.reindex(idx).to_numpy()
            frame[f"gauge_{gid}_cum"] = amount.reindex(idx).cumsum().to_numpy()
        used.append({"gauge_id": gid, "network": net, "distance_to_path_m": round(dist), "resolution": g_res})
    if gauge_mode != "none" and not used:
        notes.append("no gauge matches this choice for this link")
    subs = d.sublinks(cml_ids=[cml_id]).set_index("sublink")
    baselines = {f"s{k}": float(subs.dry_tl_median.get(k, np.nan)) for k in (1, 2)}
    ev = d.events()
    ev = ev[(ev.end > t0) & (ev.start < t1)]
    runs = d.table("event_runs")
    runs = runs[(runs.end > t0) & (runs.start < t1)]
    gaps = d.table("network_gaps")
    gaps = gaps[(gaps.end > t0) & (gaps.start < t1)]
    n_values = int(frame.select_dtypes("number").notna().sum().sum())
    info = {"baselines": baselines, "events": ev, "runs": runs, "gaps": gaps, "gauges": used, "notes": notes,
            "n_rows": len(frame), "n_values": n_values, "load_s": time.perf_counter() - tic,
            "cml_id": cml_id, "t0": t0, "t1": t1, "res": res, "qc": qc}
    return frame, info


def event_all_links(d: CMLData, event_id: str, pre_h: float = 2.0) -> tuple[pd.DataFrame, pd.Series]:
    """Descriptive event view of every link: max TL rise above the sub-link's dry median (max of the two sub-links)
    and the share of links >= 1 dB above it at each 5 min step (event start - pre_h to event end)."""
    ev = d.events().set_index("event_id").loc[event_id]
    p = d.catalog["period"]
    t0 = max(pd.Timestamp(ev.start) - pd.Timedelta(hours=pre_h), pd.Timestamp(p["start"]))
    t1 = pd.Timestamp(ev.end)
    links = d.links()
    subs = d.sublinks()
    ds = d.cml(links.cml_id.tolist(), t0, t1, "5min", variables=["tl_mean"])
    tl = ds.tl_mean.to_pandas()                                  # time x sublink_id
    base = subs.set_index("sublink_id").dry_tl_median.reindex(tl.columns)
    rise = tl - base.values
    by_link = rise.T.groupby(lambda s: int(s.split("_")[0])).max().T  # time x cml_id (max of the two sub-links)
    summary = links[["cml_id", "length_m", "split_link", "status_best"]].copy()
    summary["length_km"] = summary.length_m / 1000
    summary["max_rise_db"] = summary.cml_id.map(by_link.max())
    summary = summary.dropna(subset=["max_rise_db"])
    valid = by_link.notna().sum(axis=1).replace(0, np.nan)
    share = (by_link >= 1.0).sum(axis=1) / valid
    return summary, share


def availability(d: CMLData) -> pd.DataFrame:
    """Valid TL share per link and day (mean of the two sub-links), from the 1 h statistics."""
    p = d.catalog["period"]
    links = d.links().cml_id.tolist()
    n = d.cml(links, p["start"], p["end_exclusive"], "1h", variables=["n_valid_tl"]).n_valid_tl.to_pandas()
    day = n.groupby(n.index.floor("1D")).sum() / 8640.0
    return day.T.groupby(lambda s: int(s.split("_")[0])).mean().T.clip(0, 1)


def dry_baseline(d: CMLData, cml_id: int) -> tuple[pd.DataFrame, pd.Series, dict]:
    """Hourly TL median over the whole period, the dry-period mask and the dry medians of both sub-links."""
    p = d.catalog["period"]
    ds = d.cml(cml_id, p["start"], p["end_exclusive"], "1h", variables=["tl_median"])
    tl = ds.tl_median.to_pandas()
    tl.columns = [f"tl_median_s{c.split('_')[1]}" for c in tl.columns]
    mask = pd.Series(False, index=tl.index)
    for r in d.dry_periods().itertuples():
        mask |= (tl.index >= r.start) & (tl.index < r.end)
    subs = d.sublinks(cml_ids=[cml_id]).set_index("sublink")
    return tl, mask, {f"s{k}": float(subs.dry_tl_median.get(k, np.nan)) for k in (1, 2)}


def snippet(release: str, cml_id: int, t0, t1, res: str, qc: str, gauge_mode: str, radius: float) -> str:
    """Python that reproduces the data of the current view with CMLData only."""
    return f'''from cmlpipe.api import CMLData
d = CMLData.open(r"{release}")
df = d.frame({cml_id}, "{pd.Timestamp(t0)}", "{pd.Timestamp(t1)}", res="{res}", qc="{qc}",
             include=("cml", "flags", "radar", "met", "state"))
# gauges as in the explorer ("{gauge_mode}", radius {radius:g} m): rows of link_gauge_matches for this link
m = d.table("link_gauge_matches").query("cml_id == {cml_id}").sort_values("distance_to_path_m")
# e.g. g = d.gauge(m.gauge_id.iloc[0], "{pd.Timestamp(t0)}", "{pd.Timestamp(t1)}", res="15min")
'''


def timed(fn, *a, **k):
    t = time.perf_counter()
    out = fn(*a, **k)
    return out, time.perf_counter() - t
