"""Delete the unimportant outputs of the research modules (the output part of /tidy). Rules: ORGANIZATION.md,
section "Outputs and importance"; this script only applies them.

    python research/common/tidy_outputs.py                  # dry run: what would be deleted, kept, and why
    python research/common/tidy_outputs.py --apply          # delete, log, update MANIFEST.csv and .gitignore
    python research/common/tidy_outputs.py --all            # every module, not only those touched since their last tidy

For every research module touched since its last tidy (or every module with --all):
  delete   outputs with level diag or scratch, outputs not in outputs/MANIFEST.csv, and product files that no script
           reads and no script writes any more (orphans);
  never    a file that a script reads (its name, or the paths-file constant that holds it, appears in the code of a
           script other than its producer), a file named in the module's DECISIONS.md, or a file linked from a
           Markdown file that stays. Such a file is kept and reported as a conflict.
  keep     key files, and products that their producer still writes (an unread product, e.g. a frozen set waiting
           for the next module, is listed but kept).
Before --apply, the module's DECISIONS.md must exist and name the script of every report that is deleted (step c of
/tidy: the decisions in those reports are summarised there first).
--apply logs every deletion to <module>/outputs/TIDY_LOG.md (path, level, size, the script that rebuilds it), drops
the deleted rows from MANIFEST.csv, un-tracks kept diag / scratch files and rewrites the .gitignore block.
Only files inside <module>/outputs/ are ever deleted; code, configs, docs, data/ and archive/ are never touched.
"""
from __future__ import annotations

import argparse
import ast
import datetime as dt
import re
import subprocess
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))
from research.common.outputs import REPO, read_manifest, write_manifest, sync_gitignore  # noqa: E402

RESEARCH = REPO / "research"
BOOKKEEPING = ("MANIFEST.csv", "TIDY_LOG.md")                      # files in outputs/ that are not outputs
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
    res = {"delete": [], "conflict": [], "key": [], "product": [], "unread_product": [], "missing": [],
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
    doomed = {f.resolve() for f, *_ in cand}
    mds = [REPO / p for p in tracked if p.endswith(".md")] + [m for m in mod.rglob("*.md") if ".cache" not in m.parts]
    linked = set()
    for md in set(mds):
        if md.exists() and md.resolve() not in doomed:
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
            res["delete"].append((rel, level, f.stat().st_size, producer or "?", why))

    for rel, r in rows.items():
        if not (mod / rel).exists():
            (res["stale_rows"] if r["level"] in ("diag", "scratch") else res["missing"]).append(rel)
    return res


def check_decisions(mod, res, rows):
    """Step c of /tidy: DECISIONS.md must name the script of every report that is about to be deleted."""
    reports = [(rel, rows.get(rel, {}).get("script", "")) for rel, *_ in res["delete"] if rel.startswith("outputs/reports/")]
    if not reports:
        return []
    f = mod / "DECISIONS.md"
    text = f.read_text(encoding="utf-8") if f.exists() else ""
    return [f"{rel} (script {s or '?'})" for rel, s in reports if not s or Path(s).name not in text]


# ------------------------------------------------------------------------------------------------------- apply
def apply(mod, res, rows):
    rel_mod = mod.relative_to(REPO).as_posix()
    now = dt.datetime.now(dt.timezone.utc)
    head = git("rev-parse", "--short", "HEAD").strip() or "no-git"
    tracked = set(git("ls-files", "--", rel_mod).splitlines())
    for rel, *_ in res["delete"]:
        p = f"{rel_mod}/{rel}"
        if p in tracked:
            git("rm", "-q", "--", p)
        (mod / rel).unlink(missing_ok=True)
        rows.pop(rel, None)
    for rel in res["untrack"]:
        git("rm", "-q", "--cached", "--", f"{rel_mod}/{rel}")
    for rel in res["stale_rows"]:
        rows.pop(rel, None)
    for d in sorted((mod / "outputs").rglob("*"), key=lambda p: -len(p.parts)):   # empty sub-folders (not the 4 kinds)
        if d.is_dir() and d.parent != mod / "outputs" and not any(d.iterdir()):
            d.rmdir()
    write_manifest(mod, rows)

    log = mod / "outputs" / "TIDY_LOG.md"
    old = log.read_text(encoding="utf-8") if log.exists() else (
        f"# TIDY_LOG: outputs of {mod.name} deleted by /tidy\n\n"
        "Written by `python research/common/tidy_outputs.py --apply` (rules: ORGANIZATION.md, \"Outputs and "
        "importance\").\nNewest first. Every file can be rebuilt by running the script in the last column from the "
        "repository root.\n\n<!-- last tidy: x commit=x -->\n")
    total = sum(s for _, _, s, *_ in res["delete"])
    sec = [f"## {now:%Y-%m-%d %H:%M} UTC (commit {head}): {len(res['delete'])} files deleted, {size_str(total)}", ""]
    if res["delete"]:
        sec += ["| path | level | size | rebuilt by | why |", "|---|---|---|---|---|"]
        sec += [f"| `{r}` | {lv} | {size_str(s)} | `{p}` | {w} |" for r, lv, s, p, w in res["delete"]]
    else:
        sec.append("Nothing to delete.")
    if res["conflict"]:
        sec += ["", "Kept although the level says delete (conflicts):", ""]
        sec += [f"- `{r}` ({lv}): {why}" for r, lv, why in res["conflict"]]
    stamp = f"<!-- last tidy: {now:%Y-%m-%dT%H:%M:%SZ} commit={head} -->"
    new = MARK.sub(stamp, old, count=1)
    new = new.replace(stamp + "\n", stamp + "\n\n" + "\n".join(sec) + "\n", 1)
    log.write_text(new, encoding="utf-8")


# ------------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--apply", action="store_true", help="delete and log (default: dry run)")
    ap.add_argument("--all", action="store_true", help="every module, not only the touched ones")
    ap.add_argument("modules", nargs="*", help="module names (default: all touched modules)")
    a = ap.parse_args()

    code = corpus()
    tracked = git("ls-files").splitlines()
    mods = [m for m in modules() if not a.modules or m.name in a.modules]
    bad = 0
    print(f"tidy_outputs: {'APPLY' if a.apply else 'dry run'}")
    for mod in mods:
        rows = read_manifest(mod)
        if not (a.all or a.modules or touched(mod, rows)):
            print(f"\n== {mod.name}: not touched since the last tidy ({last_tidy(mod)[0]}), skipped")
            continue
        res = plan(mod, code, tracked)
        missing_dec = check_decisions(mod, res, rows)
        total = sum(s for _, _, s, *_ in res["delete"])
        print(f"\n== {mod.name}")
        print(f"delete: {len(res['delete'])} files, {size_str(total)}")
        for r, lv, s, p, w in res["delete"]:
            print(f"   - {r}  [{w}, {size_str(s)}; rebuilt by {p}]")
        print(f"keep, key ({len(res['key'])}): " + ", ".join(Path(r).name for r in res["key"]))
        print(f"keep, product read by a script ({len(res['product'])}): " + ", ".join(Path(r).name for r in res["product"]))
        if res["unread_product"]:
            print(f"keep, product nothing reads yet ({len(res['unread_product'])}): "
                  + ", ".join(Path(r).name for r in res["unread_product"]))
        for r, lv, why in res["conflict"]:
            print(f"CONFLICT (kept): {r} [{lv}]: {why}")
        for r in res["missing"]:
            print(f"missing key / product (rebuild with its script): {r}")
        if missing_dec:
            print("DECISIONS.md does not yet summarise these reports (step c of /tidy):")
            for m in missing_dec:
                print(f"   - {m}")
        if a.apply:
            if missing_dec:
                print(f"NOT applied to {mod.name}: update DECISIONS.md first")
                bad = 1
                continue
            apply(mod, res, rows)
            print(f"applied: log in {mod.name}/outputs/TIDY_LOG.md")
    if a.apply:
        sync_gitignore()
        print("\n.gitignore: key / product list rewritten")
    return bad


if __name__ == "__main__":
    sys.exit(main())
