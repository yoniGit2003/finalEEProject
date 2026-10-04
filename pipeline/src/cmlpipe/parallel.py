"""Run per-day work on several processes (python -m cmlpipe run ... --workers N).

Why processes and not threads: the per-day work is NumPy-heavy Python code; separate processes avoid the
GIL and keep memory per day separate. Windows starts worker processes with "spawn": each worker is a fresh
Python that imports this package, so
  * the task functions must be importable at module level (no lambdas, no nested functions),
  * everything passed to them must be picklable: we pass the config as a plain dict plus paths, never
    open files or the logger,
  * run.py keeps its `if __name__ == "__main__":` guard, otherwise every worker would re-run the CLI.
Each worker handles WHOLE days, and a day's result does not depend on which worker computed it, so
workers=1 and workers=N give identical files (tested).

A task function has the signature task(cfg_data, repo_root, config_file, day_text, **kwargs) -> dict and
returns a manifest record. It must not raise: an exception is caught here and returned as
status "failed" with the traceback, so the parent can write the manifest before stopping.
"""

from __future__ import annotations

import multiprocessing as mp
import time
import traceback
from functools import partial
from pathlib import Path

import pandas as pd

from cmlpipe.config import Config, config_hash


def config_from_parts(cfg_data: dict, repo_root: str, config_file: str) -> Config:
    """Rebuild a Config inside a worker (the dict travels, not the object)."""
    return Config(data=cfg_data, repo_root=Path(repo_root), config_file=Path(config_file), hash=config_hash(cfg_data))


def _guarded(task, cfg_data, repo_root, config_file, day_text, kwargs) -> dict:
    """Call task; turn an exception into a 'failed' record (an exception would lose the other days' records)."""
    t0 = time.perf_counter()
    try:
        return task(cfg_data, repo_root, config_file, day_text, **kwargs)
    except Exception as err:  # noqa: BLE001 - reported to the parent, which stops the run
        return {"day": day_text, "status": "failed", "error": f"{type(err).__name__}: {err}",
                "traceback": traceback.format_exc(), "runtime_s": round(time.perf_counter() - t0, 2)}


def run_days(task, cfg: Config, days: list[pd.Timestamp], workers: int = 1, **kwargs):
    """Yield one record per day, IN DAY ORDER, computed by `task` on `workers` processes.

    workers <= 1 runs in this process (easier to debug, same code path otherwise).
    """
    args = [(cfg.data, str(cfg.repo_root), str(cfg.config_file), f"{pd.Timestamp(d):%Y-%m-%d}", kwargs) for d in days]
    if workers <= 1 or len(days) <= 1:
        for a in args:
            yield _guarded(task, *a)
        return
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=min(workers, len(days))) as pool:
        # imap keeps the order of `days`; chunksize 1 because one day is already a large unit of work
        yield from pool.imap(partial(_star, task), args, chunksize=1)


def _star(task, a):
    return _guarded(task, *a)


def run_day_task(task, cfg: Config, days: list[pd.Timestamp], workers: int, force: bool, log,
                 records: list[dict]) -> list[dict]:
    """Run a per-day task on `workers` processes; log its messages and append its records in day order.

    Records are appended as they arrive, so days finished before a failure still reach the manifest; the
    first failed day stops the run (after its record is written). Keys starting with "_" (e.g. the qc daily
    summary rows) are returned but kept out of the manifest.
    """
    done = []
    for rec in run_days(task, cfg, days, workers, force=force):
        for level, text in rec.pop("messages", []):
            getattr(log, level)(text)
        extra = {k: rec.pop(k) for k in [k for k in rec if k.startswith("_")]}
        records.append(rec)
        if rec["status"] == "failed":
            log.error("%s: FAILED: %s %s", rec["day"], rec["error"], rec.get("traceback", ""))
            raise RuntimeError(f"{rec['day']}: {rec['error']}")
        done.append({**rec, **extra})
    return done
