"""Command-line interface: one entry point for every stage.

    python -m cmlpipe run <stage>[,<stage>...] [--all | --days D1,D2] [--force] [--workers N] [--config PATH]
    python -m cmlpipe run all                     # every stage in pipeline order; up-to-date stages are skipped
    python -m cmlpipe stages                      # list the stages
    python -m cmlpipe explore [--release PATH]    # start the data explorer (Streamlit) in the browser
    python -m cmlpipe guide [--release PATH]      # rewrite DATA_GUIDE.md, LICENSE.txt, README.md of a release

Stages (pipeline order): see STAGES. Day stages (ingest_cml, qc_flags, aggregate_cml) take --days or --all
(`run all` implies --all). Every stage skips outputs that are up to date (same stage hash and inputs); --force
rebuilds them. The frozen split is never redrawn unless --force-split is given (never implied by --force).
The release stage refuses to overwrite an existing release version without --force (--version v2 for a new one).
Each run appends an entry to <output_root>/run_manifest.json.

The old form `python pipeline/run.py --config ... --steps a,b` still works (deprecated shim).
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

from cmlpipe.config import DEFAULT_CONFIG, Config, load_config
from cmlpipe.logging_utils import new_run_id, setup_logging
from cmlpipe.provenance import append_to_manifest, display_path, git_commit, utc_now_iso

STAGES = ["ingest_cml", "ingest_gauges", "ingest_met", "ingest_radar", "qc_stats", "qc_flags", "qc_status",
          "geometry", "aggregate_gauges", "aggregate_radar", "aggregate_cml", "events", "split", "event_links",
          "significance", "inventory", "event_links_extra", "link_state", "inventory_figures", "quicklooks",
          "link_cards", "run_report", "release"]
# Stages that need an explicit day selection (--days or --all); the others always use the whole period.
DAY_STAGES = ("ingest_cml", "qc_flags", "aggregate_cml")


def _stage_calls(cfg: Config, args, days: list[pd.Timestamp], run_entry: dict, log) -> dict:
    """stage name -> zero-argument callable. Imports are local so `python -m cmlpipe stages` stays fast."""
    from cmlpipe.aggregate import stage as aggregate_stage
    from cmlpipe.events import link_state
    from cmlpipe.events import stage as events_stage
    from cmlpipe.geometry import stage as geometry_stage
    from cmlpipe.inventory import stage as inventory_stage
    from cmlpipe.io import ingest
    from cmlpipe.io import openmrg_reference as ref
    from cmlpipe.parallel import run_day_task
    from cmlpipe.qc import stage as qc_stage
    from cmlpipe.release import build as release_build
    from cmlpipe.viz import inventory_figures, link_cards, quicklooks, run_report

    f, d, force, w = run_entry["files"], run_entry["days"], args.force, args.workers
    return {
        "ingest_cml": lambda: ingest.ingest_cml(cfg, days, force, log, d, w),
        "ingest_gauges": lambda: ingest.ingest_gauges(cfg, force, log, f),
        "ingest_met": lambda: ingest.ingest_met(cfg, force, log, f),
        "ingest_radar": lambda: ingest.ingest_radar(cfg, [m for m in ref.radar_months(cfg)
                                                          if m in {(x.year, x.month) for x in days}], force, log, f),
        "qc_stats": lambda: qc_stage.run_qc_stats(cfg, force, log, f),
        "qc_flags": lambda: qc_stage.run_qc_flags(cfg, days, force, log, d, w),
        "qc_status": lambda: qc_stage.run_qc_status(cfg, force, log, f),
        "geometry": lambda: geometry_stage.run_geometry(cfg, force, log, f),
        "aggregate_gauges": lambda: aggregate_stage.run_aggregate_gauges(cfg, force, log, f),
        "aggregate_radar": lambda: aggregate_stage.run_aggregate_radar(cfg, force, log, f),
        "aggregate_cml": lambda: run_day_task(aggregate_stage.aggregate_cml_day, cfg, days, w, force, log, d),
        "events": lambda: events_stage.run_events(cfg, force, log, f),
        "split": lambda: events_stage.run_split(cfg, args.force_split, log, f),
        "event_links": lambda: events_stage.run_event_links(cfg, force, log, f),
        "significance": lambda: inventory_stage.run_significance(cfg, force, log, f),
        "inventory": lambda: inventory_stage.run_inventory(cfg, force, log, f),
        "event_links_extra": lambda: inventory_stage.run_event_links_extra(cfg, force, log, f),
        "link_state": lambda: link_state.run_link_state(cfg, force, log, f),
        "inventory_figures": lambda: inventory_figures.run_inventory_figures(cfg, force, log, f),
        "quicklooks": lambda: quicklooks.run_quicklooks(cfg, force, log, f),
        "link_cards": lambda: link_cards.run_link_cards(cfg, force, log, f),
        "run_report": lambda: run_report.run_report(cfg, force, log, f),
        "release": lambda: release_build.run_release(cfg, force, log, f, version=args.version),
    }


def build_parser() -> argparse.ArgumentParser:
    """The argparse parser of `python -m cmlpipe` (sub-commands run and stages; legacy --steps)."""
    p = argparse.ArgumentParser(prog="python -m cmlpipe", description="CML research pipeline")
    sub = p.add_subparsers(dest="command")
    r = sub.add_parser("run", help="run one or more stages (or 'all')")
    r.add_argument("stages", help=f"'all' or comma-separated stages: {', '.join(STAGES)}")
    sub.add_parser("stages", help="list the stages in pipeline order")
    for name, text in (("explore", "start the data explorer in the browser (Streamlit)"),
                       ("guide", "rewrite the text files of a release (DATA_GUIDE.md, LICENSE.txt, README.md); data untouched")):
        e = sub.add_parser(name, help=text)
        e.add_argument("--release", default=None, help="release folder (default <release_root>/openmrg_v1)")
        e.add_argument("--config", default=None, help=f"YAML config (default {DEFAULT_CONFIG.name} in configs/)")
        if name == "explore":
            e.add_argument("--port", type=int, default=8501, help="local port (default 8501)")
            e.add_argument("--headless", action="store_true", help="do not open a browser")
    for q in (p, r):
        q.add_argument("--config", default=None, help=f"YAML config (default {DEFAULT_CONFIG.name} in configs/)")
        q.add_argument("--steps", default=None, help=argparse.SUPPRESS)  # legacy form of `run`
        which = q.add_mutually_exclusive_group(required=False)
        which.add_argument("--days", nargs="+", help="UTC days, e.g. 2015-07-28 2015-07-29 or 2015-07-28,2015-07-29")
        which.add_argument("--all", action="store_true", help="every day of the configured period")
        q.add_argument("--force", action="store_true", help="rebuild outputs even if they are up to date")
        q.add_argument("--workers", type=int, default=1,
                       help="processes for the per-day stages (ingest_cml, qc_flags, aggregate_cml); default 1")
        q.add_argument("--force-split", action="store_true",
                       help="allow the split stage to OVERWRITE the frozen train/test split (never implied by --force)")
        q.add_argument("--version", default="v1", help="release version folder (release stage), default v1")
    return p


def main(argv=None) -> int:
    """Parse the command line, run the requested stages in pipeline order, write the manifest entry."""
    args = build_parser().parse_args(argv)
    if args.command == "stages":
        print("\n".join(STAGES))
        return 0
    if args.command in ("explore", "guide"):
        cfg = load_config(args.config or DEFAULT_CONFIG)
        release = Path(args.release) if args.release else cfg.repo_path(cfg["paths"]["release_root"]) / "openmrg_v1"
        if not (release / "catalog.json").exists():
            raise SystemExit(f"{release} is not a release folder (no catalog.json); build it with "
                             "`python -m cmlpipe run all` or pass --release")
        if args.command == "guide":
            from cmlpipe.release import guide
            for p in guide.write_guide(cfg, release):
                print(f"wrote {p}")
            return 0
        app = cfg.repo_root / "apps" / "explorer" / "app.py"
        cmd = [sys.executable, "-m", "streamlit", "run", str(app), "--server.port", str(args.port),
               "--server.headless", "true" if args.headless else "false",
               "--server.address", "localhost",          # only this computer can open it
               "--browser.gatherUsageStats", "false",    # no telemetry
               "--client.showErrorDetails", "false",     # friendly errors only, never a traceback in the page
               "--", "--release", str(release)]
        print("Starting the explorer:", " ".join(cmd), "\n(stop it with Ctrl+C)")
        return subprocess.call(cmd)
    spec = args.stages if args.command == "run" else args.steps
    if not spec:
        build_parser().print_help()
        return 2
    cfg = load_config(args.config or DEFAULT_CONFIG)
    if spec.strip() == "all":
        steps, args.all = list(STAGES), True
    else:
        steps = [s.strip() for s in spec.split(",") if s.strip()]
    unknown = [s for s in steps if s not in STAGES]
    if unknown:
        raise SystemExit(f"Unknown stages {unknown}; available: {STAGES}")
    if args.days is None and not args.all and any(s in DAY_STAGES for s in steps):
        raise SystemExit(f"stages {[s for s in steps if s in DAY_STAGES]} need --days or --all")
    if args.workers < 1:
        raise SystemExit("--workers must be >= 1")
    days = ([pd.Timestamp(d).normalize() for arg in args.days for d in arg.split(",") if d.strip()]
            if args.days else cfg.period_days())

    run_id = new_run_id()
    log = setup_logging(cfg.logs_dir, run_id)
    log.info("Run %s: config %s (hash %s), stages %s, %d days, force=%s",
             run_id, cfg.config_file.name, cfg.hash[:12], steps, len(days), args.force)
    run_entry = {
        "run_id": run_id, "started_utc": utc_now_iso(),
        "config_file": display_path(cfg.config_file, cfg.repo_root), "config_hash": cfg.hash,
        "stage_hashes": {s: cfg.stage_hash(s) for s in steps}, "code_version": git_commit(cfg.repo_root),
        "steps": steps, "force": args.force, "command": " ".join(sys.argv), "workers": args.workers,
        "days": [], "files": [], "stage_runtime_s": {},
    }
    if ("ingest_gauges" in steps or "ingest_met" in steps) and not args.all:
        log.info("ingest_gauges / ingest_met always process the whole period; --days is ignored for them")
    calls = _stage_calls(cfg, args, days, run_entry, log)
    status = "ok"
    try:
        for stage in STAGES:  # pipeline order, independent of the order on the command line
            if stage in steps:
                t0 = time.perf_counter()
                calls[stage]()
                run_entry["stage_runtime_s"][stage] = round(time.perf_counter() - t0, 2)
                log.info("%s: stage finished in %.1f s", stage, run_entry["stage_runtime_s"][stage])
    except Exception:
        status = "failed"
        raise
    finally:
        run_entry["status"] = status
        run_entry["finished_utc"] = utc_now_iso()
        append_to_manifest(cfg.manifest_path, run_entry)
        log.info("Manifest updated: %s", cfg.manifest_path)
    return 0
