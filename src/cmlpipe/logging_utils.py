"""Logging to the console and to one log file per run.

Each run gets an id (UTC time it started). The log file
data/work/<dataset>/logs/<run_id>.log keeps everything that was printed, so
we can look back later at what a run did (links dropped, days skipped, warnings).
"""

from __future__ import annotations

import logging
from datetime import datetime, timezone
from pathlib import Path

LOGGER_NAME = "cmlpipe"


def new_run_id() -> str:
    """Return a run id such as '20260930T141502Z' (UTC, second resolution)."""
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def setup_logging(logs_dir: Path, run_id: str, level: int = logging.INFO) -> logging.Logger:
    """Configure the 'cmlpipe' logger to write to the console and to <logs_dir>/<run_id>.log.

    All modules use logging.getLogger(__name__); their names start with 'cmlpipe',
    so their messages go through the handlers set here.
    Calling this twice (e.g. in tests) replaces the handlers instead of adding
    duplicates, so each message is printed only once.
    """
    logs_dir = Path(logs_dir)
    logs_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(level)
    for handler in list(logger.handlers):
        handler.close()
        logger.removeHandler(handler)

    fmt = logging.Formatter("%(asctime)s %(levelname)-7s %(name)s: %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler()
    console.setFormatter(fmt)
    logfile = logging.FileHandler(logs_dir / f"{run_id}.log", encoding="utf-8")
    logfile.setFormatter(fmt)
    logger.addHandler(console)
    logger.addHandler(logfile)
    return logger
