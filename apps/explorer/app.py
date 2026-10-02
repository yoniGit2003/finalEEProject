"""CML data explorer (Streamlit). Start with `python -m cmlpipe explore` or start_explorer.bat.

Reads ONLY through cmlpipe.api.CMLData (via data.py); figures come from cmlpipe.viz.explorer. This file only wires
controls to those functions. Every selection lives in st.session_state, so it survives reruns and tab changes.
"""

from __future__ import annotations

import argparse
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
import data as X  # noqa: E402
from cmlpipe.viz import explorer as V  # noqa: E402

REPO = Path(__file__).resolve().parents[2]


def _release() -> str:
    p = argparse.ArgumentParser()
    p.add_argument("--release", default=str(REPO / "data" / "release" / "openmrg_v1"))
    return p.parse_known_args(sys.argv[1:])[0].release


RELEASE = _release()
PLOT_CONFIG = {"displaylogo": False, "toImageButtonOptions": {"format": "png", "scale": 2}}


# ------------------------------------------------------------------------------------------- cached data access
@st.cache_resource(show_spinner=False)
def get_d(path: str):
    return X.open_release(path)


@st.cache_data(show_spinner=False)
def get_links(path: str) -> pd.DataFrame:
    return X.links_table(get_d(path))


@st.cache_data(show_spinner=False)
def get_events(path: str) -> pd.DataFrame:
    return get_d(path).events()


@st.cache_data(show_spinner=False, max_entries=48)
def get_view(path, cml_id, t0, t1, res, qc, gauge_mode, radius):
    return X.load_view(get_d(path), cml_id, t0, t1, res, qc=qc, gauge_mode=gauge_mode, gauge_radius_m=radius)


@st.cache_data(show_spinner=False, max_entries=16)
def get_event_all(path, event_id):
    return X.event_all_links(get_d(path), event_id)


@st.cache_data(show_spinner=False)
def get_availability(path):
    return X.availability(get_d(path))


@st.cache_data(show_spinner=False, max_entries=32)
def get_baseline(path, cml_id):
    return X.dry_baseline(get_d(path), cml_id)


@st.cache_data(show_spinner=False)
def get_table(path, name):
    return get_d(path).table(name)


# ------------------------------------------------------------------------------------------- state helpers
S = st.session_state
PERIOD = None


def _init(key, value):
    if key not in S:
        S[key] = value


def _set_free_window(t0: pd.Timestamp, t1: pd.Timestamp):
    t0, t1 = max(t0, PERIOD[0]), min(t1, PERIOD[1])
    S.f_d0, S.f_tm0, S.f_d1, S.f_tm1 = t0.date(), t0.time(), t1.date(), t1.time()


def _free_window():
    return pd.Timestamp.combine(S.f_d0, S.f_tm0), pd.Timestamp.combine(S.f_d1, S.f_tm1)


def _shift(delta):
    t0, t1 = _free_window()
    _set_free_window(t0 + delta, t1 + delta)


def _whole_month():
    t0, _ = _free_window()
    m0 = t0.normalize().replace(day=1)
    _set_free_window(m0, m0 + pd.offsets.MonthBegin(1))


def _step_link(options, k):
    if S.cml_id in options:
        S.cml_id = options[(options.index(S.cml_id) + k) % len(options)]
    elif options:
        S.cml_id = options[0]


def _set_res(res):
    S.res = res


def _random_link(options):
    if options:
        S.cml_id = int(np.random.default_rng().choice(options))


def _event_from_table():
    sel = S.get("ev_table")
    rows = sel.selection.rows if sel is not None else []
    if rows:
        S.win_mode = "event"
        S.ev_filter = "All events"
        S.event_id = S.ev_table_df.event_id.iloc[rows[0]]


# ------------------------------------------------------------------------------------------- sidebar
def sidebar(links: pd.DataFrame, events: pd.DataFrame):
    sb = st.sidebar
    sb.title("CML explorer")
    sb.caption(f"Release: `{Path(RELEASE).name}` (read through CMLData)")

    # --- link picker
    with sb.expander("Filter links", expanded=False):
        lmin, lmax = float(links.length_m.min() / 1000), float(links.length_m.max() / 1000)
        length = st.slider("length (km)", 0.0, float(np.ceil(lmax)), (0.0, float(np.ceil(lmax))), 0.1, key="f_len")
        band = st.multiselect("frequency band", sorted(links.frequency_band.unique()), key="f_band")
        pol = st.multiselect("polarization", sorted(links.polarization.dropna().unique()), key="f_pol")
        status = st.multiselect("status (best sub-link)", ["good", "degraded", "unusable"], key="f_status")
        split = st.multiselect("split", sorted(links.split_link.unique()), key="f_split")
        gw = st.selectbox("has a gauge within", ["any distance", 250, 500, 1000, 2000, 5000], key="f_gauge")
        comp = st.slider("min TL completeness", 0.0, 1.0, 0.0, 0.05, key="f_comp")
    sub = X.filter_links(links, length if length != (0.0, float(np.ceil(lmax))) else None, band, pol, status, split,
                         None if gw == "any distance" else gw, comp)
    options = sub.cml_id.tolist() or links.cml_id.tolist()
    if len(sub) == 0:
        sb.warning("No link matches the filters; showing all links.")
    _init("cml_id", 10316 if 10316 in options else options[0])
    if S.cml_id not in options:
        S.cml_id = options[0]
    lab = links.set_index("cml_id")
    sb.selectbox(f"Link ({len(options)} of {len(links)})", options, key="cml_id",
                 format_func=lambda c: f"{c}  ({lab.length_m[c] / 1000:.2f} km, {lab.frequency_ghz_sublink_1[c]:.0f} GHz, "
                                       f"{lab.status_best[c]})")
    c1, c2, c3 = sb.columns(3)
    c1.button("◀ prev", on_click=_step_link, args=(options, -1), width="stretch")
    c2.button("next ▶", on_click=_step_link, args=(options, 1), width="stretch")
    c3.button("random", on_click=_random_link, args=(options,), width="stretch")

    # --- window
    sb.divider()
    _init("win_mode", "event")
    sb.radio("Window", ["event", "free"], key="win_mode", horizontal=True,
             format_func=lambda m: "event ± padding" if m == "event" else "no event, free window")
    if S.win_mode == "event":
        _init("ev_filter", "significant only")
        sb.radio("Events", ["All events", "significant only", "by split", "by month"], key="ev_filter", horizontal=True)
        ev = events
        if S.ev_filter == "significant only":
            ev = ev[ev.significant]
        elif S.ev_filter == "by split":
            sp = sb.selectbox("split", ["train", "val", "test"], key="ev_split")
            ev = ev[ev.split_event == sp]
        elif S.ev_filter == "by month":
            mo = sb.selectbox("month", [6, 7, 8], key="ev_month", format_func=lambda m: {6: "June", 7: "July", 8: "August"}[m])
            ev = ev[ev.month == mo]
        ids = ev.event_id.tolist() or events.event_id.tolist()
        _init("event_id", "E20150727_2140" if "E20150727_2140" in ids else ids[0])
        if S.event_id not in ids:
            S.event_id = ids[0]
        evi = events.set_index("event_id")
        sb.selectbox("event", ids, key="event_id",
                     format_func=lambda e: f"{e}  ({evi.duration_h[e]:.1f} h, {evi.rain_class[e]}, {evi.split_event[e]})")
        _init("pad_h", 6)
        sb.slider("padding (hours)", 0, 48, key="pad_h")
        e = evi.loc[S.event_id]
        t0, t1 = e.start - pd.Timedelta(hours=S.pad_h), e.end + pd.Timedelta(hours=S.pad_h)
        t0, t1 = max(t0, PERIOD[0]), min(t1, PERIOD[1])
    else:
        if "f_d0" not in S:
            _set_free_window(pd.Timestamp("2015-07-28"), pd.Timestamp("2015-07-29"))
        a, b = sb.columns(2)
        a.date_input("start (UTC)", key="f_d0", min_value=PERIOD[0].date(), max_value=PERIOD[1].date())
        b.time_input("start time", key="f_tm0", step=300)
        a.date_input("end (UTC)", key="f_d1", min_value=PERIOD[0].date(), max_value=PERIOD[1].date())
        b.time_input("end time", key="f_tm1", step=300)
        q = sb.columns(3)
        for i, (lab_, dt) in enumerate((("1 h", "1h"), ("6 h", "6h"), ("1 day", "1D"))):
            q[i].button(f"◀ {lab_}", key=f"bk{dt}", on_click=_shift, args=(-pd.Timedelta(dt),), width="stretch")
            q[i].button(f"{lab_} ▶", key=f"fw{dt}", on_click=_shift, args=(pd.Timedelta(dt),), width="stretch")
        w = sb.columns(2)
        w[0].button("whole month", on_click=_whole_month, width="stretch")
        w[1].button("whole period", on_click=_set_free_window, args=PERIOD, width="stretch")
        t0, t1 = _free_window()
        if t1 <= t0:
            sb.error("End must be after start.")
            st.stop()

    # --- resolution
    sb.divider()
    _init("res", "5min")
    sb.radio("Resolution", X.RESOLUTIONS, key="res", horizontal=True)
    if S.res == "10s" and (t1 - t0) > pd.Timedelta(days=X.MAX_DAYS_10S):
        sb.warning(f"10 s is limited to {X.MAX_DAYS_10S} days (this window is {(t1 - t0) / pd.Timedelta('1D'):.1f} days).")
        sb.button("switch to 1 min", on_click=_set_res, args=("1min",))
        st.stop()

    # --- layers
    with sb.expander("Layers", expanded=True):
        defaults = dict(tl=True, rsl=False, tsl=False, radar=True, gauge=True, gauge_cum=False, met=False, flags=True,
                        state=True, events=True, gaps=True, baseline=True, band=True)
        names = dict(tl="TL (both sub-links)", rsl="RSL", tsl="TSL", radar="radar along path (R_mean, R_max, wet fraction)",
                     gauge="gauge rate", gauge_cum="gauge cumulative (mm)", met="meteorology", flags="QC flags",
                     state="reference state shading (wet / unknown)", events="event bounds + wet runs",
                     gaps="network gaps", baseline="dry-weather baseline (TL)", band="TL min-max band + std")
        layers = {}
        for k, v in defaults.items():
            _init(f"ly_{k}", v)
            layers[k] = st.checkbox(names[k], key=f"ly_{k}")
        _init("g_mode", "nearest City + SMHI")
        st.selectbox("gauges", ["nearest City", "nearest SMHI", "nearest City + SMHI", "all within radius", "none"],
                     key="g_mode")
        _init("g_radius", 2000)
        st.slider("radius (m)", 250, 5000, step=250, key="g_radius", disabled=S.g_mode != "all within radius")

    with sb.expander("Display", expanded=False):
        _init("qc", "mask")
        st.radio("blocked samples", ["mask", "raw"], key="qc", horizontal=True,
                 format_func=lambda q: "masked (NaN)" if q == "mask" else "raw (as recorded)",
                 help="10 s only: samples with a blocking QC flag. Coarser statistics always exclude them.")
        _init("mode", "lines")
        st.radio("draw as", ["lines", "markers"], key="mode", horizontal=True)
        _init("same_db", False)
        st.checkbox("same dB span on TL / RSL / TSL panels", key="same_db")
        _init("relative", False)
        st.checkbox("TL relative to dry baseline", key="relative")
        _init("log_rain", False)
        st.checkbox("log scale for rain", key="log_rain")
        _init("palette", "default")
        st.radio("palette", list(V.PALETTES), key="palette", horizontal=True)
    return int(S.cml_id), pd.Timestamp(t0), pd.Timestamp(t1), S.res, layers


# ------------------------------------------------------------------------------------------- tabs
def tab_link(frame, info, layers):
    with st.sidebar.expander("Advanced: any column", expanded=False):
        cols = [c for c in frame.columns if pd.api.types.is_numeric_dtype(frame[c])]
        extra_cols = st.multiselect("columns", cols, key="adv_cols")
        target = st.selectbox("add to panel", ["new panel", "tl", "rsl", "tsl", "radar", "gauge", "met"], key="adv_panel")
    opts = dict(palette=S.palette, mode=S.mode, same_db_span=S.same_db, relative=S.relative, log_rain=S.log_rain,
                extra={("extra" if target == "new panel" else target): extra_cols} if extra_cols else {},
                uirevision=f"{info['cml_id']}-{info['t0']}-{info['t1']}-{info['res']}")
    fig, notes = V.link_view_figure(frame, info, layers, opts)
    st.plotly_chart(fig, width="stretch", theme="streamlit", config=PLOT_CONFIG)
    for n in notes + info["notes"]:
        st.caption("ℹ️ " + n)
    if info["gauges"]:
        st.caption("Gauges: " + "; ".join(f"{g['gauge_id']} ({g['network']}, {g['distance_to_path_m']} m from the path, "
                                         f"{g['resolution']})" for g in info["gauges"]))


def tab_map(links, cml_id, info):
    c = st.radio("colour links by", ["status", "split", "completeness"], horizontal=True, key="map_color")
    dark = getattr(getattr(st.context, "theme", None), "type", "light") == "dark"
    fig = V.map_figure(links, get_d(RELEASE).gauges(), c, cml_id, [g["gauge_id"] for g in info["gauges"]], dark=dark)
    st.plotly_chart(fig, width="stretch", config=PLOT_CONFIG)
    st.caption("Radar cells along the path are not part of the release (only their number, `n_cells`), so they are "
               "not drawn. Map tiles need an internet connection.")


def tab_events(events):
    st.markdown("Click a row to open that event in the Link view (window = event ± padding).")
    cols = ["event_id", "start", "end", "duration_h", "significant", "rain_class", "intensity_class", "split_event",
            "month", "city_total_mean_mm", "city_max_5min_rate_mmh", "radar_domain_total_mm", "n_runs", "wet_core_hours"]
    df = events[[c for c in cols if c in events]].reset_index(drop=True)
    S.ev_table_df = df
    st.dataframe(df, key="ev_table", on_select=_event_from_table, selection_mode="single-row", hide_index=True,
                 height=320, width="stretch")
    ev = S.get("event_id") or df.event_id.iloc[0]
    st.subheader(f"All links in {ev}")
    try:
        summary, share = get_event_all(RELEASE, ev)
        a, b = V.event_response_figures(summary, share, ev)
        st.plotly_chart(a, width="stretch", theme="streamlit", config=PLOT_CONFIG)
        st.plotly_chart(b, width="stretch", theme="streamlit", config=PLOT_CONFIG)
        st.caption("Descriptive only: TL rise = 5 min mean TL minus the sub-link's dry-period median (no baseline model, "
                   "no wet-antenna correction); per link the larger of the two sub-links.")
    except Exception as err:  # noqa: BLE001
        st.warning(f"Could not build the all-links view: {err}")


def tab_info(links, cml_id):
    d = get_d(RELEASE)
    row = links[links.cml_id == cml_id].iloc[0]
    a, b = st.columns(2)
    a.markdown("**Link metadata**")
    a.dataframe(row.drop(labels=[c for c in row.index if c.startswith(("lon_", "lat_"))]).astype(str).rename("value"),
                width="stretch", height=420)
    subs = d.sublinks(cml_ids=[cml_id]).set_index("sublink_id")
    keep = ["status", "status_reason", "direction", "frequency_ghz", "completeness_tl", "completeness_tl_jun",
            "completeness_tl_jul", "completeness_tl_aug", "flag_share_missing_rsl", "flag_share_rx_floor", "flag_share_spike",
            "flag_share_flatline", "flag_share_tsl_step", "flag_share_network_gap", "dry_tl_median", "dry_tl_std",
            "dry_tl_10s_noise_median", "dry_baseline_range_db", "dry_tl_daily_drift_std", "tsl_behaviour", "rsl_min_step_db"]
    b.markdown("**Both sub-links: status and inventory metrics**")
    b.dataframe(subs[[k for k in keep if k in subs]].T.astype(str), width="stretch", height=420)
    st.markdown("**Dry-period baseline** (hourly TL median; dots = certified dry periods; dashed = dry median)")
    tl, mask, base = get_baseline(RELEASE, cml_id)
    st.plotly_chart(V.dry_baseline_figure(tl, mask, base, S.palette), width="stretch", theme="streamlit",
                    config=PLOT_CONFIG)
    st.markdown("**Gauges within 1 km of the path**")
    m = get_table(RELEASE, "link_gauge_matches")
    m = m[(m.cml_id == cml_id) & (m.distance_to_path_m <= 1000)]
    if m.empty:
        st.caption("none")
    else:
        st.dataframe(m[["gauge_id", "network", "gauge_type", "distance_to_path_m", "position_along_path"]],
                     hide_index=True)
    st.caption("Link cards (PNG) are made in the work folder, not in the release, so they are not shown here.")


def tab_network(links):
    st.markdown("**Data availability** (valid TL share per link and day, mean of the two sub-links)")
    st.plotly_chart(V.availability_figure(get_availability(RELEASE)), width="stretch", theme="streamlit",
                    config=PLOT_CONFIG)
    a, b = st.columns(2)
    a.plotly_chart(V.length_frequency_figure(links), width="stretch", theme="streamlit", config=PLOT_CONFIG)
    b.plotly_chart(V.distributions_figure(links), width="stretch", theme="streamlit", config=PLOT_CONFIG)


def tab_export(frame, info):
    st.markdown(f"Current frame: link **{info['cml_id']}**, {info['t0']} to {info['t1']} (UTC, [start, end)), "
                f"**{info['res']}**, {len(frame)} rows x {frame.shape[1]} columns (never decimated).")
    stem = f"cml_{info['cml_id']}_{info['t0']:%Y%m%d%H%M}_{info['t1']:%Y%m%d%H%M}_{info['res']}"
    a, b = st.columns(2)
    a.download_button("Download CSV", frame.to_csv().encode("utf-8"), f"{stem}.csv", "text/csv")
    buf = io.BytesIO()
    frame.to_parquet(buf)
    b.download_button("Download Parquet", buf.getvalue(), f"{stem}.parquet", "application/octet-stream")
    st.markdown("**Reproduce this view in Python**")
    st.code(X.snippet(RELEASE, info["cml_id"], info["t0"], info["t1"], info["res"], info["qc"], S.g_mode,
                      S.g_radius), language="python")
    st.info("**About this data.** OpenMRG (Gothenburg, JJA 2015): 364 links / 728 sub-links at 10 s, radar along each "
            "path, City and SMHI gauges, weather. Times are UTC and label the START of each interval; TL = TSL - RSL; "
            "NaN = unknown. Full explanation: `DATA_GUIDE.md` in the release folder (also `docs/DATA_GUIDE.md`).")


# ------------------------------------------------------------------------------------------- main
def main():
    global PERIOD
    st.set_page_config(page_title="CML explorer", layout="wide", page_icon="📡")
    try:
        d = get_d(RELEASE)
    except Exception as err:  # noqa: BLE001
        st.error(f"Cannot open the release at `{RELEASE}`: {err}. Build it with `python -m cmlpipe run all`.")
        return
    p = d.catalog["period"]
    PERIOD = (pd.Timestamp(p["start"]), pd.Timestamp(p["end_exclusive"]))
    links, events = get_links(RELEASE), get_events(RELEASE)
    cml_id, t0, t1, res, layers = sidebar(links, events)
    try:
        with st.spinner("loading ..."):
            (frame, info) = get_view(RELEASE, cml_id, t0, t1, res, S.qc, S.g_mode, S.g_radius)
    except Exception as err:  # noqa: BLE001 - friendly message, never a traceback
        st.error(f"Could not load link {cml_id} for {t0} - {t1} at {res}: {err}")
        return
    st.sidebar.caption(f"Loaded {info['n_rows']:,} time steps, {info['n_values']:,} values in {info['load_s']:.2f} s "
                       "(cached afterwards).")
    st.markdown(f"#### Link {cml_id} · {t0:%Y-%m-%d %H:%M} → {t1:%Y-%m-%d %H:%M} UTC · {res}")
    tabs = st.tabs(["Link view", "Map", "Events", "Link info", "Network", "Export"])
    for tab, fn in zip(tabs, (lambda: tab_link(frame, info, layers), lambda: tab_map(links, cml_id, info),
                              lambda: tab_events(events), lambda: tab_info(links, cml_id), lambda: tab_network(links),
                              lambda: tab_export(frame, info))):
        with tab:
            try:
                fn()
            except Exception as err:  # noqa: BLE001
                st.error(f"This tab could not be drawn: {type(err).__name__}: {err}")


main()
