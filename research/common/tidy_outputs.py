"""Move the unimportant outputs of the research modules to their local legacy/ folder (the output part of /tidy).
Nothing is ever deleted. Rules: ORGANIZATION.md, section "Outputs and importance"; this script only applies them.

    python research/common/tidy_outputs.py                  # dry run: what would move, what stays, and why
    python research/common/tidy_outputs.py --apply          # move, log, update the MANIFESTs and .gitignore
    python research/common/tidy_outputs.py --all            # every module, not only those touched since their last tidy

For every research module touched since its last tidy (or every module with --all):
  move     to <module>/legacy/<kind>/<name> (git-ignored, so local only): outputs with level diag or scratch, outputs
           not in outputs/MANIFEST.csv, and product files that no script reads and no script writes any more (orphans).
           An older legacy copy with the same name is replaced by the newer file;
  never    a file that a script reads (its name, or the paths-file constant that holds it, appears in the code of a
           script other than its producer), a file named in the module's DECISIONS.md, or a file linked from a
           Markdown file that stays. Such a file stays in outputs/ and is reported as a conflict;
  keep     key files, and products that their producer still writes (an unread product, e.g. a frozen set waiting
           for the next module, is listed but kept).
--apply logs every move to <module>/outputs/TIDY_LOG.md (path, level, size, the script that rebuilds it), moves the
rows from outputs/MANIFEST.csv to legacy/MANIFEST.csv, un-tracks moved and kept diag / scratch files (git rm --cached:
the file stays on disk) and rewrites the .gitignore block.
Only files inside <module>/outputs/ are ever moved; code, configs, docs, data/ and archive/ are never touched.
"""
from __future__ import annotations

import argparse
import ast
import datetime as dt
import os
import re
import subprocess
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from research.common.outputs import REPO, read_manifest, write_manifest, sync_gitignore  # noqa: E402

RESEARCH = REPO / "research"
BOOKKEEPING = ("MANIFEST.csv", "TIDY_LOG.md")                      # files in outputs/ that are not outputs
LEGACY = "legacy"                                                  # <module>/legacy/: local only (git-ignored)
MARK = re.compile(r"<!-- last tidy: (\S+) commit=(\S+) -->")


# ------------------------------------------------------------------------------------------------------- helpers
def git(*args):
    r = subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else ""


def modules():
    """research/<name>/ folders with a <name>_paths.py."""
    return sorted(p for p in RESEARCH.iterdir() if p.is_dir() and (p / f"{p.name}_paths.py").exists())


def size_str(n):
    return f"{n / 1e6:.1f} MB" if n >= 1e5 else f"{n / 1e3:.0f} kB"


def code_tokens(path):
    """(string constants, identifiers) used by the code of a .py file; docstrings and comments do not count."""
    try:
        tree = ast.parse(path.read_text(encoding="utf-8"))
    except (SyntaxError, UnicodeDecodeError):
        return set(), set()
    doc = {id(n.value) for n in ast.walk(tree) if isinstance(n, ast.Expr) and isinstance(n.value, ast.Constant)}
    strings = {n.value for n in ast.walk(tree)
               if isinstance(n, ast.Constant) and isinstance(n.value, str) and id(n) not in doc}
    names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
    names |= {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
    return strings, names


def path_constants(mod):
    """basename -> names of the paths-file constants that hold it (e.g. "tails_index.csv" -> {"TAILS_INDEX"})."""
    out = {}
    for line in (mod / f"{mod.name}_paths.py").read_text(encoding="utf-8").splitlines():
        m = re.match(r"""^([A-Z][A-Z0-9_]*)\s*=.*["']([^"'{}]+\.[A-Za-z0-9]+)["']""", line.split("#", 1)[0])
        if m:
            out.setdefault(m.group(2), set()).add(m.group(1))
    return out


def corpus():
    """Every research script except research/common/: path relative to REPO -> (strings, identifiers)."""
    files = [f for f in RESEARCH.rglob("*.py")
             if "common" != f.relative_to(RESEARCH).parts[0] and not {"__pycache__", ".cache"} & set(f.parts)]
    return {f.relative_to(REPO).as_posix(): code_tokens(f) for f in files}


def uses(tokens, base, consts):
    strings, names = tokens
    return any(base in s for s in strings) or bool(consts & names)


def md_links(md):
    """Absolute paths of the relative links in a Markdown file."""
    text = md.read_text(encoding="utf-8", errors="replace")
    out = set()
    for target in re.findall(r"\]\(([^)\s#]+)(?:#[^)]*)?\)", text):
        if not re.match(r"^[a-z]+:", target):
            out.add((md.parent / target).resolve())
    return out


def last_tidy(mod):
    log = mod / "outputs" / "TIDY_LOG.md"
    m = MARK.search(log.read_text(encoding="utf-8")) if log.exists() else None
    return (m.group(1), m.group(2)) if m else (None, None)


def touched(mod, rows):
    """Has the module changed since its last tidy (code, outputs, or new MANIFEST rows)?"""
    when, commit = last_tidy(mod)
    if when is None:
        return True
    rel = mod.relative_to(REPO).as_posix()
    spec = [rel] + [f":(exclude){rel}/outputs/{b}" for b in BOOKKEEPING] + [f":(exclude){rel}/DECISIONS.md"]
    if git("diff", "--name-only", commit, "--", *spec).strip() or git("status", "--porcelain", "--", *spec).strip():
        return True
    return any(r["date"] > when for r in rows.values())


# ------------------------------------------------------------------------------------------------------- the plan
def plan(mod, code, tracked):
    """Decide every output file of one module. Returns a dict of lists."""
    rows = read_manifest(mod)
    out_dir = mod / "outputs"
    files = sorted(f for f in out_dir.rglob("*") if f.is_file() and f.relative_to(out_dir).as_posix() not in BOOKKEEPING)
    consts = path_constants(mod)
    paths_file = (mod / f"{mod.name}_paths.py").relative_to(REPO).as_posix()
    decisions = (mod / "DECISIONS.md").read_text(encoding="utf-8") if (mod / "DECISIONS.md").exists() else ""
    res = {"legacy": [], "conflict": [], "key": [], "product": [], "unread_product": [], "missing": [],
           "stale_rows": [], "untrack": []}

    def readers(f, producer):
        base = f.name
        c = consts.get(base, set())
        return [p for p, tok in code.items() if p not in (producer, paths_file) and uses(tok, base, c)]

    cand = []
    for f in files:
        rel = f.relative_to(mod).as_posix()
        r = rows.get(rel)
        level = r["level"] if r else "unlisted"
        producer = r["script"] if r else ""
        if level == "key":
            res["key"].append(rel)
            continue
        if level == "product":
            if readers(f, producer):
                res["product"].append(rel)
                continue
            prod = REPO / producer
            still = prod.exists() and uses(code.get(producer, code_tokens(prod)), f.name, consts.get(f.name, set()))
            if still:
                res["unread_product"].append(rel)
                continue
            why = "product that no script reads or writes any more"
        elif level in ("diag", "scratch"):
            why = f"level {level}"
        else:
            why = "not in MANIFEST.csv"
        cand.append((f, rel, level, producer, why))

    # protection (rule b): read by a script, named in DECISIONS.md, linked from a Markdown file that stays
    moving = {f.resolve() for f, *_ in cand}
    mds = [REPO / p for p in tracked if p.endswith(".md")] + [m for m in mod.rglob("*.md") if ".cache" not in m.parts]
    linked = set()
    for md in set(mds):
        if md.exists() and md.resolve() not in moving:
            linked |= md_links(md)
    for f, rel, level, producer, why in cand:
        reasons = []
        rd = readers(f, producer)
        if rd:
            reasons.append("read by " + ", ".join(rd))
        if f.name in decisions:
            reasons.append("named in DECISIONS.md")
        if f.resolve() in linked:
            reasons.append("linked from a Markdown file that stays")
        if reasons:
            res["conflict"].append((rel, level, "; ".join(reasons)))
            if level in ("diag", "scratch") and f"{mod.relative_to(REPO).as_posix()}/{rel}" in tracked:
                res["untrack"].append(rel)
        else:
            res["legacy"].append((rel, level, f.stat().st_size, producer or "?", why))

    for rel, r in rows.items():
        if not (mod / rel).exists():
            (res["stale_rows"] if r["level"] in ("diag", "scratch") else res["missing"]).append(rel)
    return res


def reports_not_in_decisions(mod, res, rows):
    """Step reports that move to legacy/ although DECISIONS.md does not name their script yet (a note, not a blocker:
    the report is kept; DECISIONS.md is the committed summary that the research log is written from)."""
    reports = [(rel, rows.get(rel, {}).get("script", "")) for rel, *_ in res["legacy"] if rel.startswith("outputs/reports/")]
    f = mod / "DECISIONS.md"
    text = f.read_text(encoding="utf-8") if f.exists() else ""
    return [f"{rel} (script {s or '?'})" for rel, s in reports if not s or Path(s).name not in text]


def legacy_path(rel):
    """outputs/<kind>/<name> -> legacy/<kind>/<name> (both relative to the module folder)."""
    return LEGACY + "/" + rel.split("/", 1)[1]


# ------------------------------------------------------------------------------------------------------- apply
def apply(mod, res, rows):
    rel_mod = mod.relative_to(REPO).as_posix()
    now = dt.datetime.now(dt.timezone.utc)
    head = git("rev-parse", "--short", "HEAD").strip() or "no-git"
    tracked = set(git("ls-files", "--", rel_mod).splitlines())
    legacy_rows = read_manifest(mod, LEGACY)
    for rel, level, *_ in res["legacy"]:
        p = f"{rel_mod}/{rel}"
        if p in tracked:
            git("rm", "-q", "--cached", "--", p)                    # out of git only: the file stays on disk
        new = legacy_path(rel)
        (mod / new).parent.mkdir(parents=True, exist_ok=True)
        os.replace(mod / rel, mod / new)                            # a move; replaces an older legacy copy
        row = rows.pop(rel, None) or {"level": level, "purpose": "", "script": "?", "commit": "", "date": ""}
        legacy_rows[new] = {**row, "path": new}
    for rel in res["untrack"]:
        git("rm", "-q", "--cached", "--", f"{rel_mod}/{rel}")
    for rel in res["stale_rows"]:
        rows.pop(rel, None)
    for d in sorted((mod / "outputs").rglob("*"), key=lambda p: -len(p.parts)):   # empty sub-folders (not the 4 kinds)
        if d.is_dir() and d.parent != mod / "outputs" and not any(d.iterdir()):
            d.rmdir()
    write_manifest(mod, rows)
    if legacy_rows:
        write_manifest(mod, legacy_rows, LEGACY)

    log = mod / "outputs" / "TIDY_LOG.md"
    old = log.read_text(encoding="utf-8") if log.exists() else (
        f"# TIDY_LOG: outputs of {mod.name} moved to legacy/ by /tidy\n\n"
        "Written by `python research/common/tidy_outputs.py --apply` (rules: ORGANIZATION.md, \"Outputs and "
        "importance\").\nNewest first. `legacy/` is local only (git-ignored); a moved file can also be rebuilt by "
        "running the script in the last column from the repository root.\n\n<!-- last tidy: x commit=x -->\n")
    total = sum(s for _, _, s, *_ in res["legacy"])
    sec = [f"## {now:%Y-%m-%d %H:%M} UTC (commit {head}): {len(res['legacy'])} files moved to legacy/, "
           f"{size_str(total)}", ""]
    if res["legacy"]:
        sec += ["| path (now under legacy/) | level | size | rebuilt by | why |", "|---|---|---|---|---|"]
        sec += [f"| `{r}` | {lv} | {size_str(s)} | `{p}` | {w} |" for r, lv, s, p, w in res["legacy"]]
    else:
        sec.append("Nothing to move.")
    if res["conflict"]:
        sec += ["", "Kept in outputs/ although the level says legacy (conflicts):", ""]
        sec += [f"- `{r}` ({lv}): {why}" for r, lv, why in res["conflict"]]
    stamp = f"<!-- last tidy: {now:%Y-%m-%dT%H:%M:%SZ} commit={head} -->"
    new = MARK.sub(stamp, old, count=1)
    new = new.replace(stamp + "\n", stamp + "\n\n" + "\n".join(sec) + "\n", 1)
    log.write_text(new, encoding="utf-8")


# ------------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--apply", action="store_true", help="move to legacy/ and log (default: dry run)")
    ap.add_argument("--all", action="store_true", help="every module, not only the touched ones")
    ap.add_argument("modules", nargs="*", help="module names (default: all touched modules)")
    a = ap.parse_args()

    code = corpus()
    tracked = git("ls-files").splitlines()
    mods = [m for m in modules() if not a.modules or m.name in a.modules]
    print(f"tidy_outputs: {'APPLY' if a.apply else 'dry run'}")
    for mod in mods:
        rows = read_manifest(mod)
        if not (a.all or a.modules or touched(mod, rows)):
            print(f"\n== {mod.name}: not touched since the last tidy ({last_tidy(mod)[0]}), skipped")
            continue
        res = plan(mod, code, tracked)
        no_decision = reports_not_in_decisions(mod, res, rows)
        total = sum(s for _, _, s, *_ in res["legacy"])
        print(f"\n== {mod.name}")
        print(f"move to legacy/: {len(res['legacy'])} files, {size_str(total)}")
        for r, lv, s, p, w in res["legacy"]:
            print(f"   - {r}  [{w}, {size_str(s)}; rebuilt by {p}]")
        print(f"keep, key ({len(res['key'])}): " + ", ".join(Path(r).name for r in res["key"]))
        print(f"keep, product read by a script ({len(res['product'])}): " + ", ".join(Path(r).name for r in res["product"]))
        if res["unread_product"]:
            print(f"keep, product nothing reads yet ({len(res['unread_product'])}): "
                  + ", ".join(Path(r).name for r in res["unread_product"]))
        for r, lv, why in res["conflict"]:
            print(f"CONFLICT (kept in outputs/): {r} [{lv}]: {why}")
        for r in res["missing"]:
            print(f"missing key / product (rebuild with its script): {r}")
        if no_decision:
            print("note: DECISIONS.md does not summarise these step reports yet (they are kept in legacy/):")
            for m in no_decision:
                print(f"   - {m}")
        if a.apply:
            apply(mod, res, rows)
            print(f"applied: log in {mod.name}/outputs/TIDY_LOG.md")
    if a.apply:
        sync_gitignore()
        print("\n.gitignore: key / product list rewritten")
    return 0


if __name__ == "__main__":
    sys.exit(main())
