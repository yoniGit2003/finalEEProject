"""Plotly figure builders for the data explorer (apps/explorer). No Streamlit and no file access here.

Every function takes DataFrames that came from cmlpipe.api.CMLData (see apps/explorer/data.py) and returns a
plotly Figure, so the same figures can be made in a notebook or a script (apps/explorer/smoke.py).

The link view is ONE figure with stacked panels sharing a single time axis. Column names follow CMLData.frame():
<var>_s1 / <var>_s2 for the two sub-links (tl_mean_s1 at 1 min and coarser, tl_s1 at 10 s), radar_*, met_*,
flags_or_s1 / flags_s1, link_state; gauge columns are gauge_<id>_rate and gauge_<id>_cum (added by the app).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from cmlpipe.viz.style import (C_AQUA, C_BLUE, C_MAGENTA, C_MUTED, C_ORANGE, C_RED, C_YELLOW, FLAG_C, SPLIT_C,
                               STATUS_C)

# Two palettes: the project palette and Okabe-Ito (colour-blind safe). Index 0 / 1 = sub-link 1 / 2.
PALETTES = {
    "default": [C_BLUE, C_ORANGE, C_AQUA, C_MAGENTA, C_YELLOW, C_RED, "#4a3aa7", C_MUTED],
    "colour-blind safe": ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#F0E442", "#D55E00", "#56B4E9", "#000000"],
}
STATE_FILL = {1: "rgba(42,120,214,0.16)", -1: "rgba(137,135,129,0.14)", 0: "rgba(27,175,122,0.06)"}
STATE_NAME = {1: "wet", 0: "dry", -1: "unknown"}
FLAG_BITS = {"missing_rsl": 1, "missing_tsl": 2, "out_of_range": 4, "rx_floor": 8, "spike": 16, "tsl_step": 32,
             "flatline": 64, "network_gap": 128}
MAX_POINTS = 20_000      # traces longer than this are min-max decimated for display
MAX_SHADES = 400         # more state runs than this: no background shading (the state row in the flags panel stays)


# ---------------------------------------------------------------------------------------------------- helpers
def minmax_decimate(s: pd.Series, max_points: int = MAX_POINTS) -> pd.Series:
    """Keep the min and the max of each of max_points / 2 equal buckets (peaks survive, unlike plain thinning)."""
    s = s.dropna()
    if len(s) <= max_points:
        return s
    n_b = max_points // 2
    b = np.arange(len(s)) * n_b // len(s)
    v = s.to_numpy()
    df = pd.DataFrame({"b": b, "v": v, "i": np.arange(len(s))})
    keep = np.union1d(df.loc[df.groupby("b").v.idxmin(), "i"], df.loc[df.groupby("b").v.idxmax(), "i"])
    return s.iloc[keep]


def runs(series: pd.Series):
    """(value, start, end) for each run of equal values; end = start of the next run (interval-start labels)."""
    s = series.dropna()
    if s.empty:
        return []
    v = s.to_numpy()
    idx = s.index
    step = (idx[1] - idx[0]) if len(idx) > 1 else pd.Timedelta("5min")
    cut = np.flatnonzero(v[1:] != v[:-1]) + 1
    starts = np.r_[0, cut]
    ends = np.r_[cut, len(v)]
    return [(v[a], idx[a], idx[e - 1] + step) for a, e in zip(starts, ends)]


def tl_cols(frame: pd.DataFrame, var: str = "tl"):
    """Names of the central value of `var` (tl / rsl / tsl) for the two sub-links, at any resolution."""
    out = []
    for s in ("s1", "s2"):
        for c in (f"{var}_mean_{s}", f"{var}_{s}"):
            if c in frame:
                out.append(c)
                break
    return out


def baseline_shift(frame: pd.DataFrame, baselines: dict) -> dict:
    """Column -> dry baseline (dB) to subtract for 'TL relative to dry baseline'."""
    out = {}
    for c in tl_cols(frame, "tl"):
        b = baselines.get(c[-2:])
        if b is not None and np.isfinite(b):
            out[c] = b
    return out


# ---------------------------------------------------------------------------------------------------- link view
def _trace(x, y, name, color, mode, dash=None, width=1.3, showlegend=True, legendgroup=None, fill=None, opacity=1.0):
    big = len(x) > 5000
    cls = go.Scattergl if big and fill is None else go.Scatter
    kw = dict(x=x, y=y, name=name, legendgroup=legendgroup or name, showlegend=showlegend, opacity=opacity)
    if mode == "markers":
        kw.update(mode="markers", marker=dict(size=3, color=color))
    else:
        kw.update(mode="lines", line=dict(color=color, width=width, dash=dash), connectgaps=False)
    if fill:
        kw.update(fill=fill, fillcolor=color, line=dict(width=0, color=color), mode="lines", hoverinfo="skip")
    return cls(**kw)


def link_view_figure(frame: pd.DataFrame, info: dict, layers: dict, opts: dict) -> tuple[go.Figure, list[str]]:
    """The stacked link view. Returns (figure, notes shown under the figure).

    layers: booleans tl, rsl, tsl, radar, gauge, met, flags, state, events, gaps, baseline, band.
    opts: palette, mode ('lines' | 'markers'), same_db_span, relative, log_rain, extra {panel: [columns]},
          max_points. info: baselines {'s1': dB, 's2': dB}, events, runs, gaps (DataFrames), gauges (list of ids).
    """
    pal = PALETTES[opts.get("palette", "default")]
    mode = opts.get("mode", "lines")
    maxp = int(opts.get("max_points", MAX_POINTS))
    notes: list[str] = []
    shift = baseline_shift(frame, info.get("baselines", {})) if opts.get("relative") else {}

    panels = []  # (key, title, secondary)
    if layers.get("tl") and tl_cols(frame, "tl"):
        panels.append(("tl", "TL relative to dry baseline (dB)" if shift else "TL = TSL - RSL (dB)", False))
    if layers.get("rsl") and tl_cols(frame, "rsl"):
        panels.append(("rsl", "RSL (dBm)", False))
    if layers.get("tsl") and tl_cols(frame, "tsl"):
        panels.append(("tsl", "TSL (dBm)", False))
    if layers.get("radar") and "radar_R_mean" in frame:
        panels.append(("radar", "radar along path (mm/h)", True))
    gauge_cols = [c for c in frame if c.startswith("gauge_") and c.endswith("_rate")]
    if layers.get("gauge") and gauge_cols:
        panels.append(("gauge", "gauge rate (mm/h)", True))
    if layers.get("met") and any(c.startswith("met_") for c in frame):
        panels.append(("met", "weather (degC, m/s | %)", True))
    flag_cols = [c for c in frame if c.startswith("flags_")]
    if layers.get("flags") and (flag_cols or "link_state" in frame):
        panels.append(("flags", "QC flags / reference state", False))
    extra = {k: [c for c in v if c in frame] for k, v in (opts.get("extra") or {}).items()}
    extra = {k: v for k, v in extra.items() if v}
    existing = {p[0] for p in panels}
    for k in extra:
        if k not in existing:
            panels.append((k, k, False))
    if "tl" in existing_keys(panels) and all(frame[c].isna().all() for c in tl_cols(frame, "tl")):
        notes.append("This link has no valid TL in this window (missing or blocked samples).")
    if not panels:
        fig = go.Figure()
        fig.add_annotation(text="No layer with data in this window", showarrow=False, x=0.5, y=0.5,
                           xref="paper", yref="paper")
        return fig, ["Nothing to draw: turn on a layer or pick another window."]

    heights = {"flags": 0.6, "met": 0.7, "tsl": 0.6}
    fig = make_subplots(rows=len(panels), cols=1, shared_xaxes=True, vertical_spacing=0.025,
                        row_heights=[heights.get(p[0], 1.0) for p in panels],
                        specs=[[{"secondary_y": p[2]}] for p in panels], subplot_titles=[p[1] for p in panels])
    row = {p[0]: i + 1 for i, p in enumerate(panels)}
    decimated = []

    def add(key, col, name, color, secondary=False, dash=None, values=None, showlegend=True, group=None):
        s = (values if values is not None else frame[col]).astype(float)
        d = minmax_decimate(s, maxp)
        if len(d) < s.notna().sum():
            decimated.append(name)
        fig.add_trace(_trace(d.index, d.values, name, color, mode, dash=dash, showlegend=showlegend,
                             legendgroup=group), row=row[key], col=1, secondary_y=secondary)

    # CML panels
    for key in ("tl", "rsl", "tsl"):
        if key not in row:
            continue
        for i, c in enumerate(tl_cols(frame, key)):
            sub = c[-2:]
            s = frame[c] - shift.get(c, 0.0)
            if key == "tl" and layers.get("band"):
                lo, hi = f"tl_min_{sub}", f"tl_max_{sub}"
                if lo in frame and hi in frame:
                    off = shift.get(c, 0.0)
                    col = pal[i]
                    fig.add_trace(_trace(frame.index, frame[lo] - off, f"min-max {sub}", col, "lines", showlegend=False,
                                         legendgroup=f"band{sub}", width=0), row=row[key], col=1)
                    fig.add_trace(_trace(frame.index, frame[hi] - off, f"min-max {sub}", _alpha(col, 0.18), "lines",
                                         fill="tonexty", showlegend=False, legendgroup=f"band{sub}"), row=row[key], col=1)
                std = f"tl_std_{sub}"
                if std in frame:
                    add(key, None, f"tl std {sub}", pal[i], dash="dot", values=frame[std], group=f"std{sub}")
            add(key, c, f"{c}", pal[i], values=s, group=f"{key}{sub}")
            if key == "tl" and layers.get("baseline"):
                b = info.get("baselines", {}).get(sub)
                if b is not None and np.isfinite(b):
                    y = 0.0 if c in shift else b
                    fig.add_hline(y=y, line=dict(color=pal[i], width=1, dash="dash"), row=row[key], col=1,
                                  annotation_text=f"dry median {sub}" if not shift else None)
    if opts.get("same_db_span"):
        spans = []
        for key in ("tl", "rsl", "tsl"):
            if key in row:
                v = np.concatenate([(frame[c] - shift.get(c, 0.0)).to_numpy(float) for c in tl_cols(frame, key)])
                v = v[np.isfinite(v)]
                if v.size:
                    spans.append((key, v.min(), v.max()))
        if spans:
            width = max(hi - lo for _, lo, hi in spans) + 1.0
            for key, lo, hi in spans:
                mid = (lo + hi) / 2
                fig.update_yaxes(range=[mid - width / 2, mid + width / 2], row=row[key], col=1)

    # radar
    if "radar" in row:
        add("radar", "radar_R_mean", "radar R_mean", pal[2])
        if "radar_R_max" in frame:
            add("radar", "radar_R_max", "radar R_max", pal[3], dash="dot")
        if "radar_wet_fraction" in frame:
            add("radar", "radar_wet_fraction", "wet fraction (right)", C_MUTED, secondary=True, dash="dash")
            fig.update_yaxes(range=[0, 1.05], row=row["radar"], col=1, secondary_y=True, showgrid=False)
        if opts.get("log_rain"):
            fig.update_yaxes(type="log", row=row["radar"], col=1, secondary_y=False)

    # gauges
    if "gauge" in row:
        for i, c in enumerate(gauge_cols):
            gid = c[len("gauge_"):-len("_rate")]
            add("gauge", c, f"{gid} rate", pal[(4 + i) % len(pal)])
            cum = f"gauge_{gid}_cum"
            if layers.get("gauge_cum") and cum in frame:
                add("gauge", cum, f"{gid} cumulative mm (right)", pal[(4 + i) % len(pal)], secondary=True, dash="dot")
        if opts.get("log_rain"):
            fig.update_yaxes(type="log", row=row["gauge"], col=1, secondary_y=False)

    # weather
    if "met" in row:
        for i, (c, nm) in enumerate((("met_temperature", "temperature"), ("met_wind_speed", "wind speed"),
                                     ("met_wind_gust", "wind gust"))):
            if c in frame:
                add("met", c, nm, pal[(2 + i) % len(pal)], dash="dot" if nm == "wind gust" else None)
        if "met_humidity" in frame:
            add("met", "met_humidity", "humidity % (right)", C_MUTED, secondary=True, dash="dash")

    # flags timeline: one row per flag and sub-link, plus the reference state
    if "flags" in row:
        cats = []
        for c in flag_cols:
            sub = c[-2:]
            v = frame[c].fillna(0).astype(np.int64).to_numpy()
            for name, bit in FLAG_BITS.items():
                hit = (v & bit) != 0
                if hit.any():
                    lab = f"{name} {sub}"
                    cats.append(lab)
                    x = frame.index[hit]
                    fig.add_trace(go.Scattergl(x=x, y=[lab] * len(x), mode="markers", name=lab, showlegend=False,
                                               marker=dict(size=4, color=FLAG_C.get(name, C_MUTED), symbol="line-ns-open"),
                                               hovertemplate=f"{lab}<br>%{{x}}<extra></extra>"),
                                  row=row["flags"], col=1)
        if "link_state" in frame:
            st = frame["link_state"].astype("float")
            for code, colr in ((1, C_BLUE), (0, C_AQUA), (-1, C_MUTED)):
                hit = (st == code).to_numpy()
                if hit.any():
                    x = frame.index[hit]
                    fig.add_trace(go.Scattergl(x=x, y=["reference state"] * len(x), mode="markers",
                                               name=f"state {STATE_NAME[code]}", showlegend=False,
                                               marker=dict(size=5, color=colr, symbol="square"),
                                               hovertemplate=f"{STATE_NAME[code]}<br>%{{x}}<extra></extra>"),
                                  row=row["flags"], col=1)
        if not cats and "link_state" not in frame:
            notes.append("No QC flag is set in this window.")

    # extra columns
    for i, (k, cols) in enumerate(extra.items()):
        for j, c in enumerate(cols):
            add(k, c, c, pal[(i + j) % len(pal)])

    # background shapes spanning all panels
    shapes = []
    if layers.get("state") and "link_state" in frame:
        rr = [r for r in runs(frame["link_state"].astype(float)) if r[0] in (1, -1)]
        if len(rr) <= MAX_SHADES:
            shapes += [dict(type="rect", xref="x", yref="paper", x0=a, x1=b, y0=0, y1=1, line_width=0,
                            fillcolor=STATE_FILL[int(v)], layer="below") for v, a, b in rr]
        else:
            notes.append(f"{len(rr)} wet / unknown runs: background shading skipped (see the reference state row).")
    t0, t1 = frame.index.min(), frame.index.max()
    if layers.get("events"):
        for r in info.get("runs", pd.DataFrame()).itertuples():
            shapes.append(dict(type="rect", xref="x", yref="paper", x0=max(r.start, t0), x1=min(r.end, t1), y0=0.985,
                               y1=1.0, line_width=0, fillcolor=_alpha(C_ORANGE, 0.7)))
        for e in info.get("events", pd.DataFrame()).itertuples():
            for x in (e.start, e.end):
                if t0 <= x <= t1:
                    shapes.append(dict(type="line", xref="x", yref="paper", x0=x, x1=x, y0=0, y1=1,
                                       line=dict(color=C_ORANGE, width=1.2, dash="dash")))
    if layers.get("gaps"):
        for g in info.get("gaps", pd.DataFrame()).itertuples():
            shapes.append(dict(type="rect", xref="x", yref="paper", x0=g.start, x1=g.end, y0=0, y1=1, line_width=0,
                               fillcolor="rgba(11,11,11,0.18)", layer="below"))
    fig.update_layout(shapes=list(fig.layout.shapes) + shapes)

    fig.update_layout(height=max(380, 210 * len(panels)), hovermode="x unified", margin=dict(l=60, r=30, t=40, b=30),
                      legend=dict(orientation="h", yanchor="bottom", y=1.01, x=0, font=dict(size=10)),
                      uirevision=opts.get("uirevision", "keep"))
    fig.update_xaxes(showspikes=True, spikemode="across", spikethickness=1)
    fig.update_annotations(font_size=11, x=0, xanchor="left")
    if decimated:
        notes.append(f"Display decimated (min-max, <= {maxp:,} points per trace) for: {', '.join(sorted(set(decimated)))}. "
                     "Zoom in or use a shorter window to see every sample; exports are never decimated.")
    return fig, notes


def existing_keys(panels) -> set:
    return {p[0] for p in panels}


def _alpha(hex_color: str, a: float) -> str:
    h = hex_color.lstrip("#")
    return f"rgba({int(h[0:2], 16)},{int(h[2:4], 16)},{int(h[4:6], 16)},{a})"


# ---------------------------------------------------------------------------------------------------- map
def map_figure(links: pd.DataFrame, gauges: pd.DataFrame, color_by: str, selected: int | None = None,
               near_gauges: list[str] | None = None, dark: bool = False) -> go.Figure:
    """All links as lines (lon/lat of site_0 and site_1 in `links`), coloured by status_best / split_link /
    completeness; the selected link drawn thick; gauges as markers (nearest ones labelled)."""
    fig = go.Figure()
    if color_by == "completeness":
        bins = pd.cut(links.completeness_tl_min, [-0.01, 0.5, 0.8, 0.95, 0.99, 1.0],
                      labels=["< 50 %", "50-80 %", "80-95 %", "95-99 %", ">= 99 %"])
        cats, key = list(bins.cat.categories), bins.astype(str)
        colors = dict(zip(cats, [C_RED, C_ORANGE, C_YELLOW, C_AQUA, C_BLUE]))
    else:
        col = {"status": "status_best", "split": "split_link"}[color_by]
        key = links[col].astype(str)
        cats = sorted(key.unique())
        colors = STATUS_C if color_by == "status" else SPLIT_C
    for cat in cats:
        sel = links[key == cat]
        lon, lat, txt = [], [], []
        for r in sel.itertuples():
            lab = f"link {r.cml_id}<br>{r.length_m / 1000:.2f} km, {r.frequency_ghz_sublink_1:.1f} GHz<br>{color_by}: {cat}"
            lon += [r.lon_0, r.lon_1, None]
            lat += [r.lat_0, r.lat_1, None]
            txt += [lab, lab, None]
        fig.add_trace(go.Scattermap(lon=lon, lat=lat, mode="lines", line=dict(width=2.5, color=colors.get(cat, C_MUTED)),
                                    name=f"{cat} ({len(sel)})", text=txt, hoverinfo="text"))
    fig.add_trace(go.Scattermap(lon=gauges.lon, lat=gauges.lat, mode="markers", name="gauges",
                                marker=dict(size=9, color="#4a3aa7"), text=gauges.gauge_id + " (" + gauges.network + ")",
                                hoverinfo="text"))
    center = dict(lon=float(np.nanmean(links.lon_0)), lat=float(np.nanmean(links.lat_0)))
    zoom = 9.3
    if selected is not None and (links.cml_id == selected).any():
        r = links[links.cml_id == selected].iloc[0]
        fig.add_trace(go.Scattermap(lon=[r.lon_0, r.lon_1], lat=[r.lat_0, r.lat_1], mode="lines+markers",
                                    line=dict(width=7, color=C_MAGENTA), marker=dict(size=10, color=C_MAGENTA),
                                    name=f"selected {selected}", hoverinfo="name"))
        center, zoom = dict(lon=(r.lon_0 + r.lon_1) / 2, lat=(r.lat_0 + r.lat_1) / 2), 12
        if near_gauges:
            g = gauges[gauges.gauge_id.isin(near_gauges)]
            fig.add_trace(go.Scattermap(lon=g.lon, lat=g.lat, mode="markers+text", text=g.gauge_id,
                                        textposition="top right", name="nearest gauges",
                                        marker=dict(size=14, color=C_YELLOW)))
    fig.update_layout(map=dict(style="carto-darkmatter" if dark else "carto-positron", center=center, zoom=zoom),
                      height=640, margin=dict(l=0, r=0, t=0, b=0), legend=dict(x=0, y=1, bgcolor="rgba(0,0,0,0)"))
    return fig


# ---------------------------------------------------------------------------------------------------- events
def event_response_figures(summary: pd.DataFrame, share: pd.Series, event_id: str) -> tuple[go.Figure, go.Figure]:
    """Descriptive only. summary: one row per link with length_km, max_rise_db, split_link; share: share of links
    whose TL is >= 1 dB above their dry median at each time."""
    a = go.Figure()
    for sp, g in summary.groupby("split_link"):
        a.add_trace(go.Scatter(x=g.length_km, y=g.max_rise_db, mode="markers", name=sp,
                               marker=dict(size=7, color=SPLIT_C.get(sp, C_MUTED), opacity=0.8),
                               text=g.cml_id.astype(str), hovertemplate="link %{text}<br>%{x:.2f} km<br>%{y:.1f} dB"))
    a.update_layout(title=f"{event_id}: max TL rise above dry median vs path length (descriptive)",
                    xaxis_title="path length (km)", yaxis_title="max TL rise (dB)", height=420)
    b = go.Figure(go.Scatter(x=share.index, y=share.values, mode="lines", line=dict(color=C_BLUE), name="share"))
    b.update_layout(title="share of links with TL >= 1 dB above their dry median", yaxis_title="share of links",
                    yaxis_range=[0, 1], height=320)
    return a, b


# ---------------------------------------------------------------------------------------------------- link info
def dry_baseline_figure(tl_hourly: pd.DataFrame, dry_mask: pd.Series, baselines: dict, pal: str = "default") -> go.Figure:
    """Hourly TL median over the summer: all hours faint, dry-period hours solid, dry median as a dashed line."""
    p = PALETTES[pal]
    fig = go.Figure()
    for i, c in enumerate(tl_hourly.columns):
        fig.add_trace(go.Scatter(x=tl_hourly.index, y=tl_hourly[c], mode="lines", line=dict(color=p[i], width=0.8),
                                 opacity=0.35, name=f"{c} (all hours)"))
        d = tl_hourly[c].where(dry_mask)
        fig.add_trace(go.Scatter(x=d.index, y=d, mode="markers", marker=dict(size=3, color=p[i]), name=f"{c} dry periods"))
        b = baselines.get(c[-2:])
        if b is not None and np.isfinite(b):
            fig.add_hline(y=b, line=dict(color=p[i], dash="dash", width=1), annotation_text=f"dry median {c[-2:]} {b:.1f} dB")
    fig.update_layout(height=360, yaxis_title="TL hourly median (dB)", margin=dict(t=30))
    return fig


# ---------------------------------------------------------------------------------------------------- network
def availability_figure(avail: pd.DataFrame) -> go.Figure:
    """Heatmap links x days of the valid TL share (avail: index day, columns cml_id)."""
    fig = go.Figure(go.Heatmap(z=avail.T.values, x=avail.index, y=[str(c) for c in avail.columns], zmin=0, zmax=1,
                               colorscale="Viridis", colorbar=dict(title="valid share"),
                               hovertemplate="link %{y}<br>%{x|%Y-%m-%d}<br>%{z:.2f}<extra></extra>"))
    fig.update_layout(height=700, yaxis=dict(type="category", showticklabels=False, title="links (sorted by cml_id)"),
                      margin=dict(t=30))
    return fig


def length_frequency_figure(links: pd.DataFrame) -> go.Figure:
    fig = go.Figure()
    for sp, g in links.groupby("split_link"):
        fig.add_trace(go.Scatter(x=g.length_m / 1000, y=g.frequency_ghz_sublink_1, mode="markers", name=sp,
                                 marker=dict(size=7, color=SPLIT_C.get(sp, C_MUTED), opacity=0.8),
                                 text=g.cml_id.astype(str), hovertemplate="link %{text}<br>%{x:.2f} km, %{y:.1f} GHz"))
    fig.update_layout(xaxis_title="length (km)", xaxis_type="log", yaxis_title="frequency sub-link 1 (GHz)", height=420)
    return fig


def distributions_figure(links: pd.DataFrame) -> go.Figure:
    fig = make_subplots(rows=1, cols=3, subplot_titles=["length (km)", "frequency (GHz)", "TL completeness (min of 2)"])
    for i, (v, col) in enumerate(((links.length_m / 1000, C_BLUE), (links.frequency_ghz_sublink_1, C_ORANGE),
                                  (links.completeness_tl_min, C_AQUA))):
        fig.add_trace(go.Histogram(x=v, marker_color=col, nbinsx=30, showlegend=False), row=1, col=i + 1)
    fig.update_layout(height=320, margin=dict(t=40))
    return fig
