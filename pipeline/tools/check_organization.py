"""Check the repository against the rules in ORGANIZATION.md (exit 1 if anything is wrong).

What it checks (only files tracked by git, plus a few untracked ones that must not be tracked):
  1. only one .gitignore (the root one)
  2. no caches or IDE files tracked (__pycache__, .cache/, _cache/, .pytest_cache, .idea/, *.pyc, checkpoints)
  3. protected paths: nothing under research/.vibe-wise/ or data/ is tracked
  4. no versioned names (_v2, _new, _final, _old) on active paths (archive/ and docs/organization/ are exempt)
  5. root holds only the allowed files
  6. every research module has README.md, CATALOG.md and a <module>_paths.py; no other README inside it
  7. in a module, outputs (csv, png, pdf, gz, parquet, npz, log, ...) live only in outputs/ or archive/
  8. every file in a module's outputs/ has an entry in that module's CATALOG.md
  9. every archive/ subfolder has a NOTE.md
 10. module code does not hard-code output folders (they come from <module>_paths.py) and never reads
     explorations/ or archive/
 11. relative links in Markdown files point to files that exist

Run from anywhere:  python pipeline/tools/check_organization.py   (standard library only, a few seconds)
Used by the /tidy command; it can also be added to the pre-commit hook (see ORGANIZATION.md).
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ROOT_FILES = {"README.md", "ORGANIZATION.md", "CATALOG.md", "CLAUDE.md", "LICENSE", ".gitignore", ".gitattributes"}
CACHE_PARTS = {"__pycache__", ".cache", "_cache", ".pytest_cache", ".idea", ".ipynb_checkpoints", ".vibe-wise"}
OUTPUT_SUFFIXES = {".csv", ".png", ".jpg", ".svg", ".pdf", ".gz", ".parquet", ".npz", ".npy", ".log", ".nc", ".zarr",
                   ".pkl", ".xlsx"}
VERSIONED = re.compile(r"_(v\d+|new|final|old)(?=[._/-]|$)", re.IGNORECASE)
EXEMPT_VERSIONED = ("docs/organization/", "_archive/")   # records of past names; any archive/ folder is exempt too
MODULE_KINDS = ("products", "tables", "figures", "reports")


def tracked() -> list[str]:
    out = subprocess.run(["git", "ls-files"], cwd=REPO, capture_output=True, text=True, check=True).stdout
    return [line for line in out.splitlines() if line]


def research_modules(files: list[str]) -> list[str]:
    """research/<name>/ folders that hold code (explorations/ is the playground, not a module)."""
    names = {f.split("/")[1] for f in files if f.startswith("research/") and f.count("/") >= 2 and f.endswith(".py")}
    return sorted(n for n in names if n != "explorations")


def catalog_patterns(text: str) -> list[re.Pattern]:
    """Every `backticked` token of a CATALOG as a pattern: <x> = any name part, {a,b} = one of, * = anything."""
    pats = []
    for tok in re.findall(r"`([^`\s]+)`", text):
        rx = re.escape(tok)
        rx = re.sub(r"<[^>]*>", "[^/]+", rx.replace(r"\<", "<").replace(r"\>", ">"))
        rx = re.sub(r"\\\{([^}]*)\\\}", lambda m: "(" + "|".join(m.group(1).split(",")) + ")", rx)
        rx = rx.replace(r"\.\.\.", "[^/]*").replace(r"\*", "[^/]*")
        try:
            pats.append(re.compile(rx + "$"))
        except re.error:
            pass
    return pats


def check(files: list[str], repo: Path = REPO) -> list[str]:
    """All rule violations among `files` (paths relative to `repo`), one line each."""
    bad = []

    # 1-3: .gitignore, caches, protected paths
    for f in files:
        parts = f.split("/")
        if parts[-1] == ".gitignore" and f != ".gitignore":
            bad.append(f"{f}: nested .gitignore (put the rule in the root .gitignore)")
        if CACHE_PARTS & set(parts[:-1]) or f.endswith(".pyc"):
            bad.append(f"{f}: cache / IDE / plugin file is tracked (git rm --cached it; the root .gitignore covers it)")
        if f.startswith("data/"):
            bad.append(f"{f}: data/ must not be tracked")

    # 4: versioned names
    for f in files:
        if f.startswith(EXEMPT_VERSIONED) or "/archive/" in f:
            continue
        if VERSIONED.search(f):
            bad.append(f"{f}: versioned name on an active path (one current version; old ones go to archive/)")

    # 5: root files
    for f in files:
        if "/" not in f and f not in ROOT_FILES:
            bad.append(f"{f}: not an allowed root file (allowed: {', '.join(sorted(ROOT_FILES))})")

    # 6-10: research modules
    for m in research_modules(files):
        base = f"research/{m}/"
        mod = [f[len(base):] for f in files if f.startswith(base)]
        for need in ("README.md", "CATALOG.md", f"{m}_paths.py"):
            if need not in mod:
                bad.append(f"{base}{need}: missing (every module has README.md, CATALOG.md and <module>_paths.py)")
        for f in mod:
            if f.endswith("README.md") and f != "README.md" and not f.startswith("archive/"):
                bad.append(f"{base}{f}: a module has only one README (put it in README.md)")
            if Path(f).suffix.lower() in OUTPUT_SUFFIXES and not f.startswith(("outputs/", "archive/")):
                bad.append(f"{base}{f}: output file outside outputs/ (or archive/)")
            if f.startswith("outputs/") and (f.count("/") < 2 or f.split("/")[1] not in MODULE_KINDS):
                bad.append(f"{base}{f}: outputs/ holds only {', '.join(MODULE_KINDS)}/")
        cat = repo / base / "CATALOG.md"
        pats = catalog_patterns(cat.read_text(encoding="utf-8")) if cat.exists() else []
        for f in mod:
            if not f.startswith("outputs/") or f.count("/") < 2:
                continue
            rel = f.split("/", 2)[2]                                       # path inside outputs/<kind>/
            cands = {f, rel, Path(f).name}
            if not any(p.match(c) for p in pats for c in cands):
                bad.append(f"{base}{f}: no entry in {base}CATALOG.md")
        for d in {f.split("/")[1] for f in mod if f.startswith("archive/") and f.count("/") >= 2}:
            if f"archive/{d}/NOTE.md" not in mod:
                bad.append(f"{base}archive/{d}/: missing NOTE.md (what it was, what replaced it, why, when)")
        for f in mod:
            if not f.endswith(".py") or f == f"{m}_paths.py":
                continue
            for i, line in enumerate((repo / base / f).read_text(encoding="utf-8").splitlines(), 1):
                code = line.split("#", 1)[0]
                if re.search(r"""(?<!"")["'](outputs|products|tables|figures|reports)["'/]""", code):
                    bad.append(f"{base}{f}:{i}: hard-coded output folder (take it from {m}_paths.py)")
                if re.search(r"""["'/](explorations|archive)["'/]""", code):
                    bad.append(f"{base}{f}:{i}: module code must not read explorations/ or archive/")

    # 11: relative Markdown links
    for f in files:
        if not f.endswith(".md") or f.startswith("docs/organization/"):
            continue
        text = (repo / f).read_text(encoding="utf-8", errors="replace")
        for target in re.findall(r"\]\(([^)\s#]+)(?:#[^)]*)?\)", text):
            if re.match(r"^[a-z]+:", target):                             # http:, https:, mailto:
                continue
            p = (repo / f).parent / target
            if not p.exists():
                bad.append(f"{f}: broken link -> {target}")
    return bad


def main() -> int:
    files = tracked()
    bad = check(files)
    for b in bad:
        print(b)
    print(f"check_organization: {'FAIL' if bad else 'OK'} ({len(bad)} problems; {len(files)} tracked files checked)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
