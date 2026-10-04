"""Provenance: record where every output came from.

Every output file carries in its global attributes:
  stage, stage_hash - which pipeline stage wrote it, and the hash of ONLY the config keys that
                  stage depends on (config.STAGE_CONFIG_KEYS). The skip rule compares stage_hash.
  config_hash   - hash of the whole configuration (information only: it changes with any key)
  code_version  - git commit of the code ("nogit" if the folder is not a git repo)
  created_utc   - when it was written
  source_doi    - the DOI of the raw dataset
The run manifest (run_manifest.json) collects one entry per run with per-day
statistics, so a whole processing history can be read in one place.
"""

from __future__ import annotations

import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import cmlpipe
from cmlpipe.config import Config


def git_commit(repo_root: Path) -> str:
    """Return the current git commit hash, with '+dirty' if there are uncommitted changes.

    Returns "nogit" when git is not installed or the folder is not a repository.
    """
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=repo_root, capture_output=True, text=True, timeout=10
        )
        if commit.returncode != 0:
            return "nogit"
        status = subprocess.run(
            ["git", "status", "--porcelain"], cwd=repo_root, capture_output=True, text=True, timeout=10
        )
        dirty = "+dirty" if status.stdout.strip() else ""
        return commit.stdout.strip() + dirty
    except (OSError, subprocess.SubprocessError):
        return "nogit"


def display_path(path: Path, repo_root: Path) -> str:
    """Path relative to the repository root if possible (portable in the manifest), else absolute."""
    path = Path(path).resolve()
    try:
        return path.relative_to(Path(repo_root).resolve()).as_posix()
    except ValueError:
        return path.as_posix()


def utc_now_iso() -> str:
    """Current UTC time as ISO 8601 text, e.g. '2026-09-30T14:15:02Z'."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def provenance_attrs(cfg: Config, stage: str) -> dict:
    """Global attributes that every output file must carry."""
    return {
        "stage": stage,
        "stage_hash": cfg.stage_hash(stage),
        "config_hash": cfg.hash,
        "config_file": cfg.config_file.name,
        "code_version": git_commit(cfg.repo_root),
        "cmlpipe_version": cmlpipe.__version__,
        "created_utc": utc_now_iso(),
        "source_doi": cfg["dataset"]["source_doi"],
    }


def write_sidecar(path: Path, cfg: Config, stage: str, extra: dict | None = None) -> Path:
    """Write <path>.meta.json with the provenance of a file that has no attributes (e.g. CSV)."""
    meta = provenance_attrs(cfg, stage)
    if extra:
        meta.update(extra)
    sidecar = Path(str(path) + ".meta.json")
    sidecar.write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return sidecar


def append_to_manifest(manifest_path: Path, run_entry: dict) -> None:
    """Append one run entry to the JSON manifest (a list of runs).

    The file is written to a temporary name and then renamed, so a crash in the
    middle never leaves a half-written, unreadable manifest.
    """
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    runs = []
    if manifest_path.exists():
        runs = json.loads(manifest_path.read_text(encoding="utf-8"))
    runs.append(run_entry)
    tmp = manifest_path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(runs, indent=2), encoding="utf-8")
    os.replace(tmp, manifest_path)
