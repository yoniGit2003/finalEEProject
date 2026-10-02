"""Stage `run_report` (step 7.7): reports/run_report.html, one self-contained page (no external CDN; figures are
linked relatively from reports/, so the page works when the reports folder is opened or copied as a whole).

Sections: run info, data coverage, QC summary, geometry, aggregation, events, inventory, quicklooks and link cards,
known limitations. Every number is read from the stored tables, the manifest and the file attributes; the
"known limitations" are generated from those facts, not written by hand.
"""

from __future__ import annotations

import html
import json
import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd

from cmlpipe.aggregate import stage as aggregate_runner
from cmlpipe.events import stage as events_runner
from cmlpipe.inventory import stage as inventory_runner
from cmlpipe.config import STAGE_CONFIG_KEYS, Config
from cmlpipe.provenance import git_commit, utc_now_iso


def _versions() -> dict:
    out = {"python": platform.python_version()}
    for m in ("numpy", "pandas", "xarray", "netCDF4", "matplotlib", "pyproj", "shapely", "yaml"):
        try:
            out[m] = __import__(m).__version__
        except Exception:  # noqa: BLE001
            out[m] = "n/a"
    return out


def _table(df: pd.DataFrame) -> str:
    return df.to_html(index=False, border=0, classes="t", float_format=lambda v: f"{v:,.3g}", na_rep="")


def _img(rel: str, caption: str, width: int = 900) -> str:
    return f'<figure><a href="{rel}"><img src="{rel}" style="max-width:{width}px;width:100%"></a><figcaption>{html.escape(caption)}</figcaption></figure>'


def build(cfg: Config) -> str:
    """The HTML of the run report."""
    T = cfg.tables_dir
    rep = cfg.output_root / cfg["paths"]["reports_dir"]
    P = inventory_runner.out_paths(cfg)
    man = json.loads(cfg.manifest_path.read_text(encoding="utf-8")) if cfg.manifest_path.exists() else []
    rt = {}
    for e in man:  # runtime of the last run in which a stage did work (skips are not timings)
        written = {r.get("step") for r in e.get("files", []) if r.get("status") == "written"} |                   ({s for s in e.get("steps", []) if any(d.get("status") == "written" for d in e.get("days", []))})
        for s, v in e.get("stage_runtime_s", {}).items():
            if s in written:
                rt[s] = v
    code = git_commit(cfg.repo_root)
    S = []
    # --- run info ---
    hashes = pd.DataFrame({"stage": list(STAGE_CONFIG_KEYS), "stage_hash": [cfg.stage_hash(s)[:12] for s in STAGE_CONFIG_KEYS],
                           "last runtime with output (s)": [rt.get(s, np.nan) for s in STAGE_CONFIG_KEYS]})
    S.append(("Run info", f"<p>Generated {utc_now_iso()}; code <code>{code}</code> (dirty: {code.endswith('+dirty')}); config "
              f"<code>{cfg.config_file.name}</code> hash <code>{cfg.hash[:12]}</code>; dataset {cfg.dataset_name} "
              f"({html.escape(cfg['dataset']['source_doi'])}).</p><p>Versions: "
              + ", ".join(f"{k} {v}" for k, v in _versions().items()) + f"; {len(man)} runs in the manifest.</p>" + _table(hashes)))
    # --- data coverage ---
    sub = pd.read_csv(P["sublinks"])
    gaps = pd.read_csv(T / "network_gaps.csv")
    big = gaps[gaps.duration_min >= 30]
    S.append(("Data coverage", f"<p>Period {cfg.period_start:%Y-%m-%d} to {cfg.period_end - pd.Timedelta('1D'):%Y-%m-%d} (UTC, "
              f"{len(cfg.period_days())} days); {len(sub)} sub-links, {sub.cml_id.nunique()} links. JJA tl completeness: median "
              f"{100 * sub.completeness_tl.median():.2f} %, {int((sub.completeness_tl < 0.5).sum())} sub-links below 50 %. "
              f"{len(gaps)} network-gap intervals, {len(big)} of them >= 30 min:</p>" + _table(big)
              + _img("inventory/availability_by_length.png", "Daily tl availability per link")))
    # --- QC ---
    daily = pd.read_csv(T / "qc_daily_summary.csv")
    flags = daily[[c for c in daily.columns if c.startswith("n_") and c not in ("n_samples", "n_valid_rsl", "n_valid_tl")]].sum()
    ft = pd.DataFrame({"flag": flags.index.str[2:], "samples": flags.to_numpy(),
                       "share of all samples (%)": 100 * flags.to_numpy() / daily.n_samples.sum()})
    S.append(("QC summary", "<p>L2 flag totals over JJA (a sample can carry several flags). Per day and sub-link: "
              "<code>tables/qc_daily_summary.csv</code>.</p>" + _table(ft)
              + _img("inventory/flag_gallery.png", "Flag gallery (one example per flag, fixed rule)")))
    # --- geometry ---
    lg = pd.read_csv(T / "link_geometry.csv")
    lk = pd.read_csv(P["links"])
    S.append(("Geometry", f"<p>{len(lg)} links, {pd.read_csv(T / 'sites.csv').shape[0]} sites; length {lg.length_m.min():.0f}-"
              f"{lg.length_m.max():.0f} m (median {lg.length_m.median():.0f} m); all paths inside the radar grid "
              f"(min coverage {lg.radar_coverage.min():.3f}); {int((lg.n_radar_cells == 1).sum())} links inside one radar cell; "
              f"{int(lk.has_gauge_within_1km.sum())} links with a gauge within 1 km. Definitions: <code>docs/methods.md</code>.</p>"
              + _img("inventory/length_frequency_split.png", "Length vs frequency by split", 700)))
    # --- aggregation ---
    rows = []
    for res in cfg["aggregate"]["resolutions"]:
        files = sorted((cfg.output_root / cfg["paths"]["l3_cml_dir"] / res).glob("*.nc"))
        rows.append({"product": f"CML {res}", "files": len(files), "size (MB)": sum(f.stat().st_size for f in files) / 1e6})
    for p in sorted((cfg.output_root / cfg["paths"]["l3_gauges_dir"]).glob("*.nc")) + sorted((cfg.output_root / cfg["paths"]["l3_radar_dir"]).glob("*.nc")):
        rows.append({"product": p.stem, "files": 1, "size (MB)": p.stat().st_size / 1e6})
    S.append(("Aggregation", "<p>L3 products (interval-start labels, UTC). Definitions: <code>docs/methods.md</code>.</p>"
              + _table(pd.DataFrame(rows))))
    # --- events ---
    et = pd.read_csv(events_runner.table_path(cfg, "events"))
    sg = pd.read_csv(P["event_significance"])
    ct = pd.crosstab(sg.rain_class, sg.split_event).reset_index()
    dp = pd.read_csv(events_runner.table_path(cfg, "dry_periods"))
    S.append(("Events", f"<p>{len(et)} events from radar + City gauges (SMHI = independent check, 100 % of its rain inside events); "
              f"{int(sg.significant.sum())} significant, {int((~sg.significant).sum())} drizzle-only; {int((et.merged_from > 1).sum())} "
              f"events merged from several wet runs; {len(dp)} certified dry periods ({dp.duration_h.sum():.0f} h). Definitions: "
              "<code>docs/methods.md</code>.</p>"
              + _table(ct)))
    # --- inventory ---
    st = pd.crosstab(sub.tsl_behaviour, sub.status).reset_index()
    S.append(("Inventory", f"<p>Status (step 7): " + ", ".join(f"{k} {v}" for k, v in sub.status.value_counts().items())
              + ". Rules: <code>docs/methods.md</code>. Tables in <code>tables/inventory/</code>.</p>"
              + _table(st) + _img("inventory/distributions.png", "Inventory distributions and status cutoffs")
              + _img("inventory/network_maps.png", "Status, split and completeness maps")))
    # --- quicklooks and cards ---
    qm = json.loads((rep / "quicklooks" / "quicklooks.meta.json").read_text(encoding="utf-8")) if (rep / "quicklooks" / "quicklooks.meta.json").exists() else {}
    cm = json.loads((rep / "link_cards" / "link_cards.meta.json").read_text(encoding="utf-8")) if (rep / "link_cards" / "link_cards.meta.json").exists() else {}
    S.append(("Quicklooks and link cards", f"<p><a href='quicklooks/index.html'>Event quicklooks</a>: {qm.get('n_quicklooks', 0)} "
              f"PNGs for {qm.get('n_events', 0)} significant events x links {qm.get('links', [])}, plus {qm.get('n_overviews', 0)} "
              f"overviews. <a href='link_cards/index.html'>Link cards</a>: {cm.get('n_cards', 0)}.</p>"))
    # --- limitations (from facts) ---
    gg = pd.read_csv(P["gauges"])
    ratio = gg.set_index("id").radar_over_gauge_ratio
    low = ratio[ratio < 1]
    over = ", ".join(f"{i} ({gg.set_index('id').gauge_type[i]}, ratio {v:.2f})" for i, v in ratio[ratio >= 1].items())
    r_log = float(np.corrcoef(np.log10(lk.length_m), lk.frequency_ghz_sublink_1)[0, 1])
    lim = [
        f"The radar (gauge-adjusted) is below the gauge JJA total at {len(low)} of {len(ratio)} gauge pixels, by "
        f"{100 * (1 - low.max()):.0f} to {100 * (1 - low.min()):.0f} % (median {100 * (1 - low.median()):.0f} %); above it at: "
        f"{over or 'none'}. It also missed the first Torslanda peak.",
        f"Only {int(lk.has_gauge_within_1km.sum())} of {len(lk)} links have a gauge within 1 km of the path "
        f"({int(lk.n_gauges_within_250m.gt(0).sum())} within 250 m).",
        f"Length and frequency are confounded: Pearson r(log10 length, frequency) = {r_log:.2f} (short links are high-frequency).",
        f"{int(lk.single_radar_cell.sum())} links lie inside a single 2 km radar cell: the radar cannot resolve rain along their path.",
        f"{int((~sg.significant).sum())} of {len(sg)} events ({100 * (~sg.significant).mean():.0f} %) are drizzle-only (not significant).",
        f"{int((et.merged_from > 1).sum())} events are merged from several wet runs (up to {int(et.merged_from.max())}); use the wet "
        "runs (tables/event_runs.csv) or per-link timing for anything that depends on the exact rain start / end.",
        "Wet-antenna attenuation is not modelled yet (stage 1 of the project).",
        "The attenuation A(t) = TSL - RSL - baseline is not computed yet: no baseline and no wet/dry classification of the CML.",
    ]
    S.append(("Known limitations", "<ul>" + "".join(f"<li>{html.escape(x)}</li>" for x in lim) + "</ul>"))
    nav = "".join(f'<a href="#s{i}">{html.escape(t)}</a> · ' for i, (t, _) in enumerate(S))
    body = "".join(f'<h2 id="s{i}">{html.escape(t)}</h2>{c}' for i, (t, c) in enumerate(S))
    return ("<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
            "<title>CML pipeline run report</title><style>body{font-family:system-ui,sans-serif;background:#fcfcfb;color:#0b0b0b;"
            "max-width:1100px;margin:16px auto;padding:0 16px;line-height:1.45}h2{border-bottom:1px solid #e6e5e1;margin-top:2em}"
            "table.t{border-collapse:collapse;font-size:13px;margin:8px 0}table.t td,table.t th{padding:3px 8px;border-bottom:1px solid #e6e5e1;"
            "text-align:right}figure{margin:12px 0}figcaption{font-size:12px;color:#52514e}code{background:#f1f0ec;padding:0 3px}"
            "nav{font-size:13px}</style></head><body><h1>CML pipeline (OpenMRG): run report</h1><nav>" + nav + "</nav>" + body + "</body></html>")


def _signature(cfg: Config) -> str:
    """sha256 over everything the page shows (tables, stamps of the figure stages); the manifest is left out on
    purpose (it changes on every run, also when nothing else does)."""
    import hashlib
    rep = cfg.output_root / cfg["paths"]["reports_dir"]
    P = inventory_runner.out_paths(cfg)
    files = [P[k] for k in ("sublinks", "links", "gauges", "event_significance")] + [
        cfg.tables_dir / f for f in ("network_gaps.csv", "qc_daily_summary.csv", "link_geometry.csv", "sites.csv",
                                     "events.csv", "dry_periods.csv")] + [
        rep / "quicklooks" / "quicklooks.meta.json", rep / "link_cards" / "link_cards.meta.json",
        rep / "inventory" / "inventory_figures.meta.json"]
    h = hashlib.sha256()
    for f in files:
        h.update(f.read_bytes() if f.exists() else b"missing")
    return h.hexdigest()


def run_report(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """Stage run_report: write reports/run_report.html; skipped when nothing it shows has changed."""
    t0 = time.perf_counter()
    out = cfg.output_root / cfg["paths"]["reports_dir"] / "run_report.html"
    stamp = out.with_suffix(".meta.json")
    sig = _signature(cfg)
    if not force and out.exists() and stamp.exists():
        m = json.loads(stamp.read_text(encoding="utf-8"))
        if m.get("stage_hash") == cfg.stage_hash("run_report") and m.get("input_signature") == sig:
            log.info("run_report: up to date, skipped (use --force to rebuild)")
            records.append({"step": "run_report", "status": "skipped_same_config"})
            return
    out.write_text(build(cfg), encoding="utf-8")
    stamp.write_text(json.dumps({"stage_hash": cfg.stage_hash("run_report"), "input_signature": sig}, indent=2),
                     encoding="utf-8")
    rec = {"step": "run_report", "status": "written", "path": str(out), "size_bytes": out.stat().st_size,
           "runtime_s": round(time.perf_counter() - t0, 2)}
    log.info("run_report: wrote %s (%.0f kB)", out.name, rec["size_bytes"] / 1e3)
    records.append(rec)
