"""Link cards (step 7.6): one page per link, reports/link_cards/<cml_id>.png, plus index.html.

Each card shows the metadata block (length, frequencies, polarization, sites, split, status of both sub-links,
gauges within 1 km, TSL behaviour), 92-day 1 h strips of tl_mean (with network-outage bands and blocking-flag
ticks), rsl_mean and tsl_mean for both sub-links, the dry-period medians of tl (baseline stability), and an inset
map of the path on the radar grid (crossed cells shaded by weight, nearest gauge). Only stored products are used:
the consolidated 1 h CML core file (built here once from the daily L3 files with netCDF4), the dry-period medians
from the 5 min core, path_grid_weights.nc and the inventory tables.

No combined PDF: pypdf is not installed and re-rendering 364 figures into PdfPages doubles the runtime; the
index.html lists every card instead.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import time
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd

from cmlpipe.events import stage as events_runner
from cmlpipe import geometry
from cmlpipe.geometry import stage as geometry_runner
from cmlpipe.inventory import stage as inventory_runner
from cmlpipe.io import schema
from cmlpipe.config import Config
from cmlpipe.parallel import config_from_parts
from cmlpipe.viz.style import C_INK, C_MUTED, C_ORANGE, C_SURF, FLAG_C, SPLIT_C, STATUS_C, SUBLINK_C, label, style

warnings.filterwarnings("ignore", category=RuntimeWarning)
BLOCK = ("rx_floor", "spike", "out_of_range", "network_gap")


def cards_dir(cfg: Config) -> Path:
    """Folder of the link cards."""
    return cfg.output_root / cfg["paths"]["reports_dir"] / "link_cards"


def _load_shared(cfg: Config) -> dict:
    """Everything a worker needs, read once per worker."""
    times, cids, subs, d = events_runner.read_core(cfg, "1h", cfg["inventory"]["core_vars_1h"])
    mids, meds, c5, s5 = inventory_runner.dry_period_medians(cfg)
    if not np.array_equal(c5, cids):
        raise ValueError("1 h and 5 min core files list the links in different orders")
    gaps = pd.read_csv(cfg.tables_dir / "network_gaps.csv", parse_dates=["start", "end"])
    weights, grid = geometry_runner.read_weights(geometry_runner.weights_path(cfg))
    gauges = pd.read_csv(cfg.tables_dir / "gauges.csv")
    gx, gy = geometry.project_to_radar(gauges.lon, gauges.lat, grid.proj_string)
    return {"times": times, "cids": cids, "subs": subs, "core": d, "mids": mids, "meds": meds.reshape(len(mids), len(subs), len(cids)),
            "outages": gaps[gaps.duration_min >= 30][["start", "end"]], "weights": weights, "grid": grid,
            "gauge_xy": dict(zip(gauges.id, zip(gx, gy)))}


def make_link_card(cfg: Config, link: dict, sub: pd.DataFrame, S: dict, out_path) -> Path:
    """Draw one card (see module docstring). link: a links.csv row as dict; sub: its two sublinks.csv rows."""
    c = int(link["cml_id"])
    j = int(np.flatnonzero(S["cids"] == c)[0])
    t = S["times"]
    fig = plt.figure(figsize=(13, 11), facecolor=C_SURF, layout="constrained")
    gs = fig.add_gridspec(5, 3, height_ratios=[1.3, 1.2, 1, 1, 1], width_ratios=[1.2, 1.2, 1])
    # --- metadata block ---
    ax = fig.add_subplot(gs[0, :2])
    ax.axis("off")
    s1, s2 = sub.set_index("sublink_id").loc["sublink_1"], sub.set_index("sublink_id").loc["sublink_2"]
    lines = [f"Link {c}   split: {link['split_link']}   hop {link['hop']}   length {link['length_m'] / 1000:.3f} km "
             f"({link['length_class']}, {link['frequency_band']})",
             f"frequency {link['frequency_ghz_sublink_1']:.3f} / {link['frequency_ghz_sublink_2']:.3f} GHz, polarization "
             f"{link['polarization']}, sites {link['site_0_id']} -> {link['site_1_id']}, radar cells {link['n_cells']}",
             f"gauges within 250 / 500 / 1000 m: {link['n_gauges_within_250m']} / {link['n_gauges_within_500m']} / "
             f"{link['n_gauges_within_1000m']}; nearest gauge {link['nearest_gauge_id']} ({link['nearest_gauge_network']}, "
             f"{link['nearest_gauge_distance_m']:.0f} m from the path; nearest City {link['nearest_city_gauge_id']} "
             f"{link['nearest_city_gauge_distance_m']:.0f} m, SMHI {link['nearest_smhi_gauge_distance_m']:.0f} m)",
             ""]
    for name, r in (("sublink_1", s1), ("sublink_2", s2)):
        lines.append(f"{name}: status {r.status} ({r.status_reason}); completeness tl {100 * r.completeness_tl:.1f} %; "
                     f"TSL {r.tsl_behaviour}; dry tl noise {r.dry_tl_std:.2f} dB; dry baseline range {r.dry_baseline_range_db:.2f} dB")
    ax.text(0, 1, "\n".join(lines), va="top", ha="left", fontsize=9, color=C_INK, family="monospace", wrap=True)
    ax.text(0, 0.02, f"status sublink_1 / sublink_2", fontsize=8, color=C_MUTED)
    for k, r in enumerate((s1, s2)):
        ax.add_patch(plt.Rectangle((0.22 + 0.05 * k, 0.0), 0.04, 0.08, transform=ax.transAxes, color=STATUS_C[r.status]))
    ax.add_patch(plt.Rectangle((0.34, 0.0), 0.04, 0.08, transform=ax.transAxes, color=SPLIT_C[link["split_link"]]))
    ax.text(0.39, 0.02, "split", fontsize=8, color=C_MUTED, transform=ax.transAxes)
    # --- inset map (radar projection: cells are exact 2 km squares there) ---
    am = fig.add_subplot(gs[0:2, 2])
    style(am)
    grid, w = S["grid"], S["weights"]
    wl = w[w.cml_id == c]
    xe, ye = grid.x_edges(), grid.y_edges()
    for r in wl.itertuples():
        am.add_patch(plt.Rectangle((min(xe[r.ix], xe[r.ix + 1]) / 1000, min(ye[r.iy], ye[r.iy + 1]) / 1000), 2, 2,
                                   color=SUBLINK_C["sublink_1"], alpha=0.15 + 0.6 * r.weight, lw=0))
    x0, y0 = geometry.project_to_radar([link["site_0_lon"]], [link["site_0_lat"]], grid.proj_string)
    x1, y1 = geometry.project_to_radar([link["site_1_lon"]], [link["site_1_lat"]], grid.proj_string)
    am.plot([x0[0] / 1000, x1[0] / 1000], [y0[0] / 1000, y1[0] / 1000], color=C_INK, lw=2, marker="o", ms=4)
    gid = link["nearest_gauge_id"]
    if isinstance(gid, str) and gid in S["gauge_xy"]:
        gx, gy = S["gauge_xy"][gid]
        am.scatter([gx / 1000], [gy / 1000], marker="^", s=60, color=C_ORANGE, zorder=5, label=f"gauge {gid}")
        am.legend(fontsize=7, frameon=False, loc="lower left")
    cx, cy = (x0[0] + x1[0]) / 2000, (y0[0] + y1[0]) / 2000
    half = max(4.0, abs(x1[0] - x0[0]) / 2000 + 3, abs(y1[0] - y0[0]) / 2000 + 3)
    am.set_xlim(cx - half, cx + half)
    am.set_ylim(cy - half, cy + half)
    am.set_xticks(np.arange(np.floor((cx - half) / 2) * 2 + (xe[0] / 1000) % 2, cx + half, 2))
    am.set_yticks(np.arange(np.floor((cy - half) / 2) * 2 + (ye[0] / 1000) % 2, cy + half, 2))
    am.tick_params(labelbottom=False, labelleft=False)
    am.set_aspect("equal")
    label(am, "radar grid (2 km cells)", "", "path on the radar grid")
    # --- 1 h strips ---
    core = S["core"]
    axs = [fig.add_subplot(gs[r, :2]) for r in (1, 2, 3)]
    for a in axs[1:]:
        a.sharex(axs[0])
    for a, v, unit in zip(axs, ("tl_mean", "rsl_mean", "tsl_mean"), ("TL (dB)", "RSL (dBm)", "TSL (dBm)")):
        style(a)
        for o in S["outages"].itertuples():
            a.axvspan(o.start, o.end, color=C_MUTED, alpha=0.25, lw=0)
        for k, s in enumerate(schema.SUBLINK_IDS):
            a.plot(t, core[v][:, k, j], color=SUBLINK_C[s], lw=0.7, label=s)
        label(a, "", unit, f"{v}, 1 h" + (" (grey: network outages >= 30 min; ticks: hours with a blocking flag)" if v == "tl_mean" else ""))
    tlv = core["tl_mean"][:, :, j]
    if np.isfinite(tlv).any():
        lo = np.nanmin(tlv)
        for k, s in enumerate(schema.SUBLINK_IDS):
            m = np.zeros(len(t), bool)
            for f in BLOCK:
                m |= (core["flags_or"][:, k, j] & np.uint16(1 << schema.FLAGS[f])) != 0
            if m.any():
                axs[0].scatter(t[m], np.full(m.sum(), lo - 1 - k), marker="|", s=25, color=SUBLINK_C[s])
    axs[0].legend(fontsize=7, frameon=False, loc="upper right")
    # --- dry baseline ---
    ab = fig.add_subplot(gs[4, :2], sharex=axs[0])
    style(ab)
    for k, s in enumerate(schema.SUBLINK_IDS):
        ab.plot(S["mids"], S["meds"][:, k, j], marker="o", ms=3, lw=0.8, color=SUBLINK_C[s], label=s)
    label(ab, "time (UTC, 2015)", "dry median TL (dB)", "dry-period medians of tl_mean (baseline stability)")
    # --- TSL class / small text ---
    at = fig.add_subplot(gs[2:, 2])
    at.axis("off")
    txt = ["TSL behaviour", f"  sublink_1: {s1.tsl_behaviour}", f"  sublink_2: {s2.tsl_behaviour}", "",
           "flag shares (JJA, %)"]
    for f in schema.FLAGS:
        txt.append(f"  {f:13s} {100 * s1[f'flag_share_{f}']:6.2f} {100 * s2[f'flag_share_{f}']:6.2f}")
    txt += ["", "completeness tl (%)  jun   jul   aug"]
    for name, r in (("sublink_1", s1), ("sublink_2", s2)):
        txt.append(f"  {name}        {100 * r.completeness_tl_jun:5.1f} {100 * r.completeness_tl_jul:5.1f} {100 * r.completeness_tl_aug:5.1f}")
    at.text(0, 1, "\n".join(txt), va="top", fontsize=8.5, family="monospace", color=C_INK)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=80, facecolor=C_SURF)
    plt.close(fig)
    return out_path


def _chunk_task(args) -> dict:
    cfg_data, repo_root, config_file, ids = args
    cfg = config_from_parts(cfg_data, repo_root, config_file)
    t0 = time.perf_counter()
    S = _load_shared(cfg)
    P = inventory_runner.out_paths(cfg)
    links = pd.read_csv(P["links"]).set_index("cml_id")
    lg = pd.read_csv(cfg.tables_dir / "link_geometry.csv").set_index("cml_id")
    sub = pd.read_csv(P["sublinks"])
    for c in ids:
        row = links.loc[c].to_dict() | lg.loc[c, ["site_0_lon", "site_0_lat", "site_1_lon", "site_1_lat"]].to_dict() | {"cml_id": c}
        make_link_card(cfg, row, sub[sub.cml_id == c], S, cards_dir(cfg) / f"{c}.png")
    return {"n": len(ids), "runtime_s": round(time.perf_counter() - t0, 2)}


def write_index(cfg: Config, links: pd.DataFrame) -> Path:
    """index.html listing every card."""
    rows = "".join(
        f'<tr><td><a href="{r.cml_id}.png">{r.cml_id}</a></td><td>{r.length_m / 1000:.2f}</td><td>{r.frequency_ghz_sublink_1:.1f}</td>'
        f'<td>{r.split_link}</td><td>{r.status_sublink_1} / {r.status_sublink_2}</td><td>{r.n_gauges_within_1000m}</td></tr>'
        for r in links.sort_values("cml_id").itertuples())
    page = ("<!doctype html><html><head><meta charset='utf-8'><title>Link cards</title><style>body{font-family:sans-serif;"
            "background:#fcfcfb;color:#0b0b0b;margin:16px}td,th{padding:3px 8px;font-size:13px;text-align:left}"
            "tr:nth-child(even){background:#f1f0ec}</style></head><body><h1>Link cards</h1>"
            f"<p>{len(links)} links. Click an id to open its card.</p><table><tr><th>cml_id</th><th>length km</th>"
            "<th>f sublink_1 GHz</th><th>split</th><th>status s1 / s2</th><th>gauges within 1 km</th></tr>"
            + rows + "</table></body></html>")
    p = cards_dir(cfg) / "index.html"
    p.write_text(page, encoding="utf-8")
    return p


def run_link_cards(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """Stage link_cards: one PNG per link (workers in parallel) plus index.html."""
    t0 = time.perf_counter()
    P = inventory_runner.out_paths(cfg)
    events_runner.build_core(cfg, "1h", log=log, variables=cfg["inventory"]["core_vars_1h"])
    sig = inventory_runner._signature([P["links"], P["sublinks"], events_runner.table_path(cfg, "dry_periods")],
                                      [events_runner.core_path(cfg, "1h"), events_runner.core_path(cfg, "5min")])
    stamp = cards_dir(cfg) / "link_cards.meta.json"
    if not force and stamp.exists():
        m = json.loads(stamp.read_text(encoding="utf-8"))
        if m.get("stage_hash") == cfg.stage_hash("link_cards") and m.get("input_signature") == sig:
            log.info("link_cards: up to date, skipped (use --force to rebuild)")
            records.append({"step": "link_cards", "status": "skipped_same_config"})
            return
    links = pd.read_csv(P["links"])
    ids = links.cml_id.tolist()
    nw = max(1, int(cfg["link_cards"]["workers"]))
    chunks = [ids[k::nw] for k in range(nw)]
    args = [(cfg.data, str(cfg.repo_root), str(cfg.config_file), ch) for ch in chunks if ch]
    if nw > 1:
        with mp.get_context("spawn").Pool(nw) as pool:
            res = pool.map(_chunk_task, args, chunksize=1)
    else:
        res = [_chunk_task(a) for a in args]
    idx = write_index(cfg, links)
    n = sum(r["n"] for r in res)
    stamp.write_text(json.dumps({"stage_hash": cfg.stage_hash("link_cards"), "input_signature": sig, "n_cards": n}, indent=2),
                     encoding="utf-8")
    rec = {"step": "link_cards", "status": "written", "n_cards": n, "index": str(idx), "workers": nw,
           "runtime_s": round(time.perf_counter() - t0, 2)}
    log.info("link_cards: %d cards (%.0f s, %d workers)", n, rec["runtime_s"], nw)
    records.append(rec)
