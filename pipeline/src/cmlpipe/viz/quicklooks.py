"""Event quicklooks (step 7.5): one PNG per (significant event, selected link) and one overview per event.

    from cmlpipe.viz.quicklooks import make_event_quicklook
    make_event_quicklook(cfg, "E20150727_2140", 10227, "ql.png")

Data come only from stored products, read with netCDF4 slices (never 92 files with xarray): the consolidated
5 min CML core file (tl_mean, flags_or), radar_along_links_5min.nc, the 5 min City / 15 min SMHI gauge
products and the step 6/7 tables. The window is the event padded by quicklooks.pad_h on both sides (shaded as
context). The overview is DESCRIPTIVE (CML response of all links), never used to define anything.

Link selection (fixed rule, config quicklooks): links not "excluded" in the split; within each
length class x frequency band cell, prefer test split, then both sub-links good, then a gauge within 1 km, then
the nearest gauge (any distance), then cml_id; take per_cell links from every cell, then fill up to n_links
round-robin. Shading: grey = padding; magenta = OTHER events that fall in the padding (with their id).
"""

from __future__ import annotations

import html
import json
import multiprocessing as mp
import time
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import netCDF4
import numpy as np
import pandas as pd

from cmlpipe.aggregate import stage as aggregate_runner
from cmlpipe.events import stage as events_runner
from cmlpipe.inventory import stage as inventory_runner
from cmlpipe.io import schema
from cmlpipe.config import Config
from cmlpipe.parallel import config_from_parts
from cmlpipe.viz.style import C_BLUE, C_INK, C_MUTED, C_ORANGE, C_AQUA, C_SURF, FLAG_C, SUBLINK_C, label, style

warnings.filterwarnings("ignore", category=RuntimeWarning)
STEP = pd.Timedelta("5min")


def ql_dir(cfg: Config) -> Path:
    """Folder of the quicklooks."""
    return cfg.output_root / cfg["paths"]["reports_dir"] / "quicklooks"


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------
def select_links(links: pd.DataFrame, qcfg: dict) -> list[int]:
    """The fixed set of quicklook links (see module docstring)."""
    lk = links[links.split_link != "excluded"].copy()
    lk["cell"] = lk.length_class + " x " + lk.frequency_band
    lk["k_split"] = (lk.split_link != "test").astype(int)
    lk["k_status"] = ~((lk.status_sublink_1 == "good") & (lk.status_sublink_2 == "good"))
    lk["k_gauge"] = ~lk.has_gauge_within_1km
    lk["k_dist"] = lk.nearest_gauge_distance_m.fillna(1e9)
    lk = lk.sort_values(["cell", "k_split", "k_status", "k_gauge", "k_dist", "cml_id"])
    queues = {c: list(g.cml_id) for c, g in lk.groupby("cell", sort=True)}
    chosen = []
    for c in queues:
        take = queues[c][: qcfg["per_cell"]]
        chosen += take
        queues[c] = queues[c][qcfg["per_cell"]:]
    while len(chosen) < qcfg["n_links"] and any(queues.values()):
        for c in queues:
            if queues[c] and len(chosen) < qcfg["n_links"]:
                chosen.append(queues[c].pop(0))
    return [int(c) for c in chosen[: qcfg["n_links"]]]


# ---------------------------------------------------------------------------
# Window loading (netCDF4 slices)
# ---------------------------------------------------------------------------
def load_window(cfg: Config, t0: pd.Timestamp, t1: pd.Timestamp) -> dict:
    """All links, 5 min bins of [t0, t1) clipped to the period: tl_mean, flags_or, radar, gauges."""
    t0, t1 = max(t0, cfg.period_start), min(t1, cfg.period_end)
    a, b = int((t0 - cfg.period_start) / STEP), int((t1 - cfg.period_start) / STEP)
    w = {"times": pd.date_range(t0, periods=b - a, freq=STEP)}
    with netCDF4.Dataset(events_runner.core_path(cfg, "5min")) as nc:
        nc.set_auto_mask(False)
        w["cml_id"] = nc["cml_id"][:].astype("int64")
        w["tl_mean"] = nc["tl_mean"][a:b].astype("float64")   # (t, sublink, cml)
        w["flags_or"] = nc["flags_or"][a:b]
    with netCDF4.Dataset(aggregate_runner.radar_path(cfg, "along_links", "5min")) as nc:
        nc.set_auto_mask(False)
        if not np.array_equal(nc["cml_id"][:].astype("int64"), w["cml_id"]):
            raise ValueError("radar product and CML core list the links in different orders")
        for v in ("R_mean", "R_max", "wet_fraction"):
            w[v] = nc[v][a:b].astype("float64")             # (t, cml)
    with netCDF4.Dataset(aggregate_runner.gauge_path(cfg, "city", "5min")) as nc:
        nc.set_auto_mask(False)
        w["city_ids"] = [str(x) for x in nc["id"][:]]
        w["city_rate"] = nc["rate_mm_h"][a:b].astype("float64")
    a15, b15 = int((t0.floor("15min") - cfg.period_start) / pd.Timedelta("15min")), int((t1.ceil("15min") - cfg.period_start) / pd.Timedelta("15min"))
    with netCDF4.Dataset(aggregate_runner.gauge_path(cfg, "smhi", "15min")) as nc:
        nc.set_auto_mask(False)
        w["smhi_times"] = pd.date_range(cfg.period_start + a15 * pd.Timedelta("15min"), periods=b15 - a15, freq="15min")
        w["smhi_rate"] = nc["rate_mm_h"][a15:b15, 0].astype("float64")
    return w


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
BLOCK_MARK = ("rx_floor", "spike", "out_of_range", "network_gap")


C_OTHER = "#e87ba4"  # other events inside the padding window (magenta, distinct from the grey padding)


def _context(ax, w, start, end, core, others=(), label_others=False):
    """Padding shaded grey, event bounds solid, wet core (first..last significant run) dotted, OTHER events that fall in
    the padding shaded magenta (labelled with their id on the top panel), so rain in the padding is explained."""
    ax.axvspan(w["times"][0], start, color=C_MUTED, alpha=0.08, lw=0)
    ax.axvspan(end, w["times"][-1] + STEP, color=C_MUTED, alpha=0.08, lw=0)
    for oid, os_, oe in others:
        ax.axvspan(max(pd.Timestamp(os_), w["times"][0]), min(pd.Timestamp(oe), w["times"][-1] + STEP), color=C_OTHER, alpha=0.15, lw=0)
        if label_others:
            ax.text(max(pd.Timestamp(os_), w["times"][0]), 0.98, oid, transform=ax.get_xaxis_transform(), fontsize=7,
                    color=C_OTHER, va="top")
    for t in (start, end):
        ax.axvline(t, color=C_INK, lw=0.9)
    if core is not None and pd.notna(core[0]):
        for t in core:
            ax.axvline(pd.Timestamp(t), color=C_INK, lw=0.8, ls=":")


def context_handles():
    """Legend entries for the shading and vertical lines drawn by _context."""
    from matplotlib.patches import Patch
    return [plt.Line2D([], [], color=C_INK, lw=0.9, label="event start / end"),
            plt.Line2D([], [], color=C_INK, lw=0.8, ls=":", label="wet core (first to last significant run)"),
            Patch(color=C_MUTED, alpha=0.25, label="padding (context)"),
            Patch(color=C_OTHER, alpha=0.4, label="other event in the padding")]


def other_events(cfg: Config, event_id: str, t0, t1) -> list:
    """(event_id, start, end) of the other events that overlap [t0, t1)."""
    et = pd.read_csv(events_runner.table_path(cfg, "events"), parse_dates=["start", "end"])
    o = et[(et.event_id != event_id) & (et.start < t1) & (et.end > t0)]
    return list(zip(o.event_id, o.start, o.end))


def make_event_quicklook(cfg: Config, event_id: str, cml_id: int, out_path, window: dict | None = None,
                         event: dict | None = None, link: dict | None = None, others: list | None = None) -> Path:
    """Four panels for one (event, link): tl_mean (both sub-links, blocking flags marked), radar along the path,
    nearest gauge rate, QC flag timeline. window / event / link may be passed in by the batch to avoid re-reading."""
    pad = pd.Timedelta(hours=cfg["quicklooks"]["pad_h"])
    if event is None:
        ev = pd.read_csv(events_runner.table_path(cfg, "events"), parse_dates=["start", "end"]).set_index("event_id").loc[event_id]
        sg = pd.read_csv(inventory_runner.out_paths(cfg)["event_significance"]).set_index("event_id").loc[event_id]
        event = {"start": ev.start, "end": ev.end, "split_event": sg.split_event, "rain_class": sg.rain_class,
                 "core": (sg.wet_core_start, sg.wet_core_end)}
    if link is None:
        lk = pd.read_csv(inventory_runner.out_paths(cfg)["links"]).set_index("cml_id").loc[cml_id]
        link = lk.to_dict()
    start, end = pd.Timestamp(event["start"]), pd.Timestamp(event["end"])
    if window is None:
        window = load_window(cfg, start - pad, end + pad)
    w = window
    if others is None:
        others = other_events(cfg, event_id, w["times"][0], w["times"][-1] + STEP)
    j = int(np.flatnonzero(w["cml_id"] == cml_id)[0])
    t = w["times"]
    fig, axes = plt.subplots(4, 1, figsize=(11, 9.5), sharex=True, facecolor=C_SURF, layout="constrained",
                             gridspec_kw={"height_ratios": [1.4, 1, 1, 0.9]})
    for k_ax, ax in enumerate(axes):
        style(ax)
        _context(ax, w, start, end, event.get("core"), others, label_others=(k_ax == 0))
    a1, a2, a3, a4 = axes
    for k, s in enumerate(schema.SUBLINK_IDS):
        a1.plot(t, w["tl_mean"][:, k, j], color=SUBLINK_C[s], lw=1.1, label=f"{s} tl_mean")
        bl = np.zeros(len(t), bool)
        for f in BLOCK_MARK:
            bl |= (w["flags_or"][:, k, j] & np.uint16(1 << schema.FLAGS[f])) != 0
        if bl.any():
            y = np.nanmin(w["tl_mean"][:, :, j]) if np.isfinite(w["tl_mean"][:, :, j]).any() else 0
            a1.scatter(t[bl], np.full(bl.sum(), y - 0.5 - k * 0.5), marker="x", s=14, color=SUBLINK_C[s],
                       label=f"{s} blocking flag in bin")
    label(a1, "", "TL (dB)", f"{event_id} x link {cml_id}  ({link['length_m'] / 1000:.2f} km, "
          f"{link['frequency_ghz_sublink_1']:.1f}/{link['frequency_ghz_sublink_2']:.1f} GHz, {link['polarization']}; "
          f"event split {event['split_event']}, link split {link['split_link']}, {event['rain_class']})")
    a1.legend(fontsize=7, frameon=False, loc="upper right", ncol=2)
    a2.plot(t, w["R_mean"][:, j], color=C_AQUA, lw=1.2, label="R_mean (path-weighted)")
    a2.plot(t, w["R_max"][:, j], color=C_AQUA, lw=0.8, ls="--", label="R_max (wettest cell)")
    trans = mtransforms.blended_transform_factory(a2.transData, a2.transAxes)
    wf = np.nan_to_num(w["wet_fraction"][:, j])
    for ti, v in zip(t, wf):
        if v > 0:
            a2.add_patch(plt.Rectangle((ti, 0.0), STEP, 0.06, transform=trans, color=C_AQUA, alpha=float(v), lw=0))
    label(a2, "", "rain (mm/h)", "radar along the path (bottom strip: wet_fraction as opacity)")
    a2.legend(fontsize=7, frameon=False, loc="upper right")
    gid, gdist = link.get("nearest_gauge_id", ""), link.get("nearest_gauge_distance_m", np.nan)
    if isinstance(gid, str) and gid in w["city_ids"]:
        a3.step(t, w["city_rate"][:, w["city_ids"].index(gid)], where="post", color=C_ORANGE, lw=1.1)
        title = f"nearest gauge {gid} (City, 5 min rate), {gdist:.0f} m from the path"
    elif isinstance(gid, str) and gid:
        a3.step(w["smhi_times"], w["smhi_rate"], where="post", color=C_ORANGE, lw=1.1)
        title = f"nearest gauge {gid} (SMHI, 15 min rate), {gdist:.0f} m from the path"
    else:
        title = "no gauge within 5 km"
    label(a3, "", "gauge (mm/h)", title)
    for fi, (f, bit) in enumerate(schema.FLAGS.items()):
        for k, s in enumerate(schema.SUBLINK_IDS):
            m = (w["flags_or"][:, k, j] & np.uint16(1 << bit)) != 0
            if m.any():
                a4.scatter(t[m], np.full(m.sum(), fi + (k - 0.5) * 0.3), marker="|", s=40, color=FLAG_C[f])
    a4.legend(handles=context_handles(), fontsize=7, frameon=False, loc="upper left", ncol=4, bbox_to_anchor=(0, -0.3))
    a4.set_yticks(range(len(schema.FLAGS)), list(schema.FLAGS), fontsize=7)
    a4.set_ylim(-0.7, len(schema.FLAGS) - 0.3)
    label(a4, "time (UTC)", "", "QC flags in the 5 min bins (upper tick sublink_2, lower sublink_1)")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=90, facecolor=C_SURF)
    plt.close(fig)
    return out_path


def make_event_overview(cfg: Config, event_id: str, w: dict, start, end, lengths_m: np.ndarray, out_path,
                        others: list | None = None) -> Path:
    """All links (DESCRIPTIVE): max tl_mean rise in the event over the 3 h pre-event median vs path length, and the
    share of sub-links more than 3 dB above their pre-event median at each 5 min bin."""
    t = w["times"]
    i0, i1 = int(np.searchsorted(t, start)), int(np.searchsorted(t, end))
    tl = w["tl_mean"]
    pre = tl[max(i0 - 36, 0):i0]
    base = np.nanmedian(pre, axis=0)                      # (sublink, cml)
    ok = np.isfinite(pre).mean(axis=0) >= 0.8
    rise = np.nanmax(tl[i0:i1], axis=0) - base
    above = (tl - base[None]) >= 3.0
    share = np.where(ok[None], above, False).sum(axis=(1, 2)) / max(ok.sum(), 1)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.6), facecolor=C_SURF, layout="constrained")
    for a in (a1, a2):
        style(a)
    for k, s in enumerate(schema.SUBLINK_IDS):
        m = ok[k]
        a1.scatter(lengths_m[m] / 1000, rise[k, m], s=12, color=SUBLINK_C[s], edgecolor=C_SURF, lw=0.4, label=s)
    if ok.any():
        a1.set_xscale("log")
        a1.legend(fontsize=8, frameon=False)
    else:  # e.g. an event at the period start: no 3 h pre-event window
        a1.text(0.5, 0.5, "no sub-link has >= 80 % valid bins in the 3 h before the event", ha="center",
                transform=a1.transAxes, fontsize=9, color=C_MUTED)
    label(a1, "path length (km, log)", "max tl_mean rise (dB)", f"{event_id}: event max over the 3 h pre-event median")
    if others is None:
        others = other_events(cfg, event_id, t[0], t[-1] + STEP)
    _context(a2, w, start, end, None, others, label_others=True)
    a2.plot(t, 100 * share, color=C_BLUE, lw=1.2)
    a2.legend(handles=[h for h in context_handles() if "wet core" not in h.get_label()], fontsize=7, frameon=False,
              loc="upper right")
    label(a2, "time (UTC)", "sub-links > 3 dB above pre-event median (%)", "Network response (descriptive)")
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=90, facecolor=C_SURF)
    plt.close(fig)
    return out_path


# ---------------------------------------------------------------------------
# Batch
# ---------------------------------------------------------------------------
def _event_task(args) -> dict:
    """One worker task: all quicklooks of one event + its overview (module level for Windows spawn)."""
    cfg_data, repo_root, config_file, ev, link_rows, lengths = args
    cfg = config_from_parts(cfg_data, repo_root, config_file)
    t0 = time.perf_counter()
    pad = pd.Timedelta(hours=cfg["quicklooks"]["pad_h"])
    start, end = pd.Timestamp(ev["start"]), pd.Timestamp(ev["end"])
    w = load_window(cfg, start - pad, end + pad)
    others = other_events(cfg, ev["event_id"], w["times"][0], w["times"][-1] + STEP)
    d = ql_dir(cfg) / ev["event_id"]
    n = 0
    for lk in link_rows:
        make_event_quicklook(cfg, ev["event_id"], lk["cml_id"], d / f"{lk['cml_id']}.png", window=w, event=ev, link=lk,
                             others=others)
        n += 1
    make_event_overview(cfg, ev["event_id"], w, start, end, np.asarray(lengths), d / "overview.png", others=others)
    return {"event_id": ev["event_id"], "n": n, "runtime_s": round(time.perf_counter() - t0, 2)}


def write_index(cfg: Config, evs: list[dict], links: list[dict]) -> Path:
    """index.html with the quicklooks of every event."""
    rows = []
    for ev in evs:
        cells = "".join(f'<td><a href="{ev["event_id"]}/{lk["cml_id"]}.png"><img src="{ev["event_id"]}/{lk["cml_id"]}.png" '
                        f'width="150" loading="lazy"></a><br>{lk["cml_id"]} ({lk["split_link"]})</td>' for lk in links)
        rows.append(f'<tr><td><b>{ev["event_id"]}</b><br>{ev["start"]:%Y-%m-%d %H:%M} – {ev["end"]:%m-%d %H:%M}<br>'
                    f'split {ev["split_event"]}<br>{ev["rain_class"]}<br><a href="{ev["event_id"]}/overview.png">overview</a></td>{cells}</tr>')
    page = ("<!doctype html><html><head><meta charset='utf-8'><title>Event quicklooks</title><style>"
            "body{font-family:sans-serif;background:#fcfcfb;color:#0b0b0b;margin:16px}td{vertical-align:top;font-size:12px;padding:4px}"
            "table{border-collapse:collapse}tr{border-bottom:1px solid #e6e5e1}</style></head><body>"
            f"<h1>Event quicklooks</h1><p>{len(evs)} significant events x {len(links)} links: "
            + html.escape(", ".join(str(lk["cml_id"]) for lk in links)) + "</p><table>" + "".join(rows) + "</table></body></html>")
    p = ql_dir(cfg) / "index.html"
    p.write_text(page, encoding="utf-8")
    return p


def run_quicklooks(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """Stage quicklooks: PNGs for every significant event x the selected links."""
    t0 = time.perf_counter()
    P = inventory_runner.out_paths(cfg)
    qcfg = cfg["quicklooks"]
    sig = inventory_runner._signature([P["event_significance"], P["links"], events_runner.table_path(cfg, "events")],
                                      [events_runner.core_path(cfg, "5min")])
    stamp = ql_dir(cfg) / "quicklooks.meta.json"
    if not force and stamp.exists():
        m = json.loads(stamp.read_text(encoding="utf-8"))
        if m.get("stage_hash") == cfg.stage_hash("quicklooks") and m.get("input_signature") == sig:
            log.info("quicklooks: up to date, skipped (use --force to rebuild)")
            records.append({"step": "quicklooks", "status": "skipped_same_config"})
            return
    links = pd.read_csv(P["links"])
    chosen = select_links(links, qcfg)
    link_rows = [links.set_index("cml_id").loc[c].to_dict() | {"cml_id": c} for c in chosen]
    sg = pd.read_csv(P["event_significance"])
    et = pd.read_csv(events_runner.table_path(cfg, "events"), parse_dates=["start", "end"]).merge(
        sg[["event_id", "significant", "rain_class", "wet_core_start", "wet_core_end"]], on="event_id")
    et = et[et.significant].sort_values("start")
    evs = [{"event_id": r.event_id, "start": r.start, "end": r.end, "split_event": r.split_event, "rain_class": r.rain_class,
            "core": (r.wet_core_start, r.wet_core_end)} for r in et.itertuples()]
    with netCDF4.Dataset(events_runner.core_path(cfg, "5min")) as nc:
        order = nc["cml_id"][:].astype("int64")
    lengths = links.set_index("cml_id").loc[order, "length_m"].to_numpy().tolist()
    args = [(cfg.data, str(cfg.repo_root), str(cfg.config_file), ev, link_rows, lengths) for ev in evs]
    if qcfg["workers"] > 1:
        with mp.get_context("spawn").Pool(qcfg["workers"]) as pool:
            res = pool.map(_event_task, args, chunksize=1)
    else:
        res = [_event_task(a) for a in args]
    idx = write_index(cfg, evs, link_rows)
    n = sum(r["n"] for r in res)
    stamp.write_text(json.dumps({"stage_hash": cfg.stage_hash("quicklooks"), "input_signature": sig, "links": chosen,
                                 "n_events": len(evs), "n_quicklooks": n, "n_overviews": len(res)}, indent=2), encoding="utf-8")
    rec = {"step": "quicklooks", "status": "written", "n_events": len(evs), "n_quicklooks": n, "n_overviews": len(res),
           "links": chosen, "index": str(idx), "runtime_s": round(time.perf_counter() - t0, 2)}
    log.info("quicklooks: %d PNGs + %d overviews for %d events x %d links (%.0f s)", n, len(res), len(evs), len(chosen), rec["runtime_s"])
    records.append(rec)
