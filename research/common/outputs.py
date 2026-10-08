"""Save a research output with an importance level and record it in the module's outputs/MANIFEST.csv.

Every file a research script writes into research/<module>/outputs/ goes through save_output() (or, for a file that
was written another way, register_output()). The level says how important the file is; /tidy reads it and deletes
the unimportant ones. The rules are in ORGANIZATION.md, section "Outputs and importance" (and only there).

    import WAA_estm_paths as P                      # the module's paths file also makes `research.` importable
    from research.common.outputs import save_output
    save_output(fig, P.STEP1_FIG_MEDIAN_CURVES, "diag", "median tail curve per A_start bin", dpi=110)
    save_output(df, P.TAILS_INDEX, "product", "one row per drying tail", index=False)

Levels:
    key      the few results the final report shows (at most 6 per module)
    product  data that a script or another module reads (frozen deliverables, inputs of later steps)
    diag     step figures, step reports and checks (the default for a step's figures and reports)
    scratch  throw-away files

MANIFEST.csv (one row per output file, rewritten on every save):
    path     relative to the module folder, e.g. outputs/figures/step1_median_curves.png
    level    key / product / diag / scratch
    purpose  one line: what the file is for
    script   the file whose code wrote it, relative to the repository root (or "by hand")
    commit   git HEAD when it was written; "+dirty" if that script had uncommitted changes
    date     UTC time of the write

Files outside an outputs/ folder (e.g. .cache/dev/ of a DEV run) are written but not recorded.
The root .gitignore ignores everything in research/*/outputs/ except MANIFEST.csv, TIDY_LOG.md and a generated list of
key / product files (<= 10 MB). That list is rewritten here, so diag / scratch files are never committed.
"""
from __future__ import annotations

import csv
import datetime as dt
import inspect
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]                 # the repository root
LEVELS = ("key", "product", "diag", "scratch")
MAX_KEY = 6                                                # at most this many key outputs per module
MAX_COMMIT_MB = 10                                         # larger files are never committed (ORGANIZATION.md)
COLUMNS = ["path", "level", "purpose", "script", "commit", "date"]
GITIGNORE = REPO / ".gitignore"
BEGIN = "# BEGIN generated: key / product outputs that are committed (research/common/outputs.py; do not edit)"
END = "# END generated"


# --------------------------------------------------------------------------------------------------- public helpers
def save_output(obj, path, level, purpose, *, script=None, **kwargs):
    """Write `obj` to `path` and record it in the module's MANIFEST.csv. Returns the path.

    obj: a matplotlib figure (savefig), a pandas DataFrame / Series (.csv, .csv.gz or .parquet), a pyarrow Table
    (.parquet), or a str (text, UTF-8). kwargs go to the writer (dpi=..., index=False, ...).
    script: the script that rebuilds the file (default: the file whose code called save_output); give it when the
    caller is a helper module that cannot be run by itself.
    """
    path = Path(path)
    _check_level(level)
    path.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(obj, "savefig"):                                          # matplotlib figure
        obj.savefig(path, **kwargs)
    elif hasattr(obj, "to_csv"):                                         # pandas DataFrame / Series
        if path.suffix == ".parquet":
            obj.to_parquet(path, **kwargs)
        else:
            obj.to_csv(path, **kwargs)
    elif type(obj).__module__.startswith("pyarrow"):                     # pyarrow Table
        import pyarrow.parquet as pq
        pq.write_table(obj, path, **kwargs)
    elif isinstance(obj, str):
        path.write_text(obj, encoding=kwargs.get("encoding", "utf-8"))
    else:
        raise TypeError(f"save_output: cannot write a {type(obj).__name__} to {path}")
    _record(path, level, purpose, script or _caller())
    return path


def register_output(path, level, purpose, script=None):
    """Record a file that was written another way (e.g. a streamed Parquet file, or a report written by hand)."""
    path = Path(path)
    _check_level(level)
    if not path.exists():
        raise FileNotFoundError(path)
    _record(path, level, purpose, script or _caller())
    return path


def read_manifest(module_dir):
    """The rows of research/<module>/outputs/MANIFEST.csv as a dict path -> row (empty if there is none)."""
    f = Path(module_dir) / "outputs" / "MANIFEST.csv"
    if not f.exists():
        return {}
    with open(f, newline="", encoding="utf-8") as fh:
        return {r["path"]: r for r in csv.DictReader(fh)}


def write_manifest(module_dir, rows):
    """Write the MANIFEST rows (dict path -> row), sorted by path, through a temporary file."""
    f = Path(module_dir) / "outputs" / "MANIFEST.csv"
    f.parent.mkdir(parents=True, exist_ok=True)
    tmp = f.with_suffix(".csv.tmp")
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        for p in sorted(rows):
            w.writerow({c: rows[p].get(c, "") for c in COLUMNS})
    tmp.replace(f)


def module_of(path):
    """(module folder, path relative to it) if `path` is inside research/<module>/outputs/, else (None, None)."""
    path = Path(path).resolve()
    parts = path.parts
    for i in range(len(parts) - 2, 1, -1):                              # .../research/<module>/outputs/...
        if parts[i] == "outputs" and parts[i - 2] == "research":
            mod = Path(*parts[:i])
            return mod, path.relative_to(mod).as_posix()
    return None, None


def sync_gitignore():
    """Rewrite the generated block of the root .gitignore: the key / product files (<= 10 MB) to commit."""
    keep = []
    for man in sorted((REPO / "research").glob("*/outputs/MANIFEST.csv")):
        mod = man.parent.parent
        for p, r in read_manifest(mod).items():
            f = mod / p
            if r["level"] in ("key", "product") and f.exists() and f.stat().st_size <= MAX_COMMIT_MB * 1e6:
                keep.append("!" + f.relative_to(REPO).as_posix())
    text = GITIGNORE.read_text(encoding="utf-8")
    if BEGIN not in text or END not in text:
        raise RuntimeError(f"{GITIGNORE}: the generated block ({BEGIN!r} ... {END!r}) is missing")
    head, rest = text.split(BEGIN, 1)
    tail = rest.split(END, 1)[1]
    new = head + BEGIN + "\n" + "".join(k + "\n" for k in sorted(keep)) + END + tail
    if new != text:
        GITIGNORE.write_text(new, encoding="utf-8")


# --------------------------------------------------------------------------------------------------- internals
def _check_level(level):
    if level not in LEVELS:
        raise ValueError(f"level must be one of {LEVELS}, not {level!r}")


def _caller():
    """The file (relative to the repository root) whose code called save_output / register_output."""
    here = Path(__file__).resolve()
    for fr in inspect.stack()[2:]:
        f = Path(fr.filename).resolve()
        if f != here:
            try:
                return f.relative_to(REPO).as_posix()
            except ValueError:
                return f.as_posix()
    return "unknown"


_COMMITS = {}


def _commit(script):
    """git HEAD (short), plus "+dirty" if the script has uncommitted changes or is not in git yet."""
    if script not in _COMMITS:
        def git(*a):
            r = subprocess.run(["git", *a], cwd=REPO, capture_output=True, text=True)
            return r.stdout.strip() if r.returncode == 0 else ""
        head = git("rev-parse", "--short", "HEAD") or "no-git"
        dirty = bool(git("status", "--porcelain", "--", script)) if (REPO / script).exists() else False
        _COMMITS[script] = head + ("+dirty" if dirty else "")
    return _COMMITS[script]


def _record(path, level, purpose, script):
    mod, rel = module_of(path)
    if mod is None:                                                      # e.g. .cache/dev/: not an output
        return
    rows = read_manifest(mod)
    rows[rel] = {"path": rel, "level": level, "purpose": purpose, "script": script,
                 "commit": "" if script == "by hand" else _commit(script),
                 "date": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
    write_manifest(mod, rows)
    n_key = sum(r["level"] == "key" for r in rows.values())
    if n_key > MAX_KEY:
        print(f"WARNING save_output: {mod.name} has {n_key} key outputs (at most {MAX_KEY}; ORGANIZATION.md)")
    if level in ("key", "product"):
        sync_gitignore()
