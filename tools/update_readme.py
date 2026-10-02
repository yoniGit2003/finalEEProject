"""Refresh the generated parts of README.md.

README.md has hand-written text and generated blocks. A generated block sits between two markers:

    <!-- AUTO:name -->
    ... replaced by this script ...
    <!-- /AUTO:name -->

Only the text between the markers is replaced; everything else is left as it is.
The script runs on every commit through the git hook in .githooks/pre-commit, and can also be run by hand:

    python tools/update_readme.py              # also runs the tests (under a minute)
    python tools/update_readme.py --no-tests   # faster, the test line says "not run"

It uses only the Python standard library, so any Python 3.9+ can run it. Facts about the code are read
by parsing the source files (ast), not by importing them, so numpy/xarray are not needed here. Only the
test run needs the project's packages; it uses the same Python that runs this script.
"""

from __future__ import annotations

import argparse
import ast
import datetime as dt
import json
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
README = REPO / "README.md"
ROADMAP = REPO / "docs" / "roadmap.md"
PACKAGE = REPO / "src" / "cmlpipe"
RUN_PY = PACKAGE / "cli.py"
SCHEMA_PY = PACKAGE / "io" / "schema.py"
TESTS = REPO / "tests"
L1_CML_DIR = REPO / "data" / "work" / "openmrg" / "L1_standard" / "cml"
MANIFEST = REPO / "data" / "work" / "openmrg" / "run_manifest.json"


# ---------------------------------------------------------------- small helpers

def rel(path: Path) -> str:
    """Path relative to the repository root, with forward slashes (as shown on GitHub)."""
    return path.relative_to(REPO).as_posix()


def parse(path: Path) -> ast.Module:
    return ast.parse(path.read_text(encoding="utf-8"), filename=str(path))


def first_doc_line(path: Path) -> str:
    """First line of a module's docstring, or '' if it has none."""
    doc = ast.get_docstring(parse(path)) or ""
    return doc.strip().splitlines()[0] if doc.strip() else ""


def module_constant(path: Path, name: str):
    """Value of a top-level constant such as STEPS = [...] (literal values only)."""
    for node in parse(path).body:
        targets = node.targets if isinstance(node, ast.Assign) else (
            [node.target] if isinstance(node, ast.AnnAssign) else [])
        for t in targets:
            if isinstance(t, ast.Name) and t.id == name and node.value is not None:
                return ast.literal_eval(node.value)
    return None


# ---------------------------------------------------------------- generated blocks

def block_implemented(run_tests: bool) -> str:
    lines = ["**Package `src/cmlpipe/`** (one line per module, from its docstring):", ""]
    for path in sorted(PACKAGE.rglob("*.py")):
        if path.name == "__init__.py" and path.parent == PACKAGE:
            continue
        name = rel(path).removeprefix("pipeline/")
        if path.name == "__init__.py":  # a sub-package, e.g. adapters/
            name = rel(path.parent).removeprefix("pipeline/") + "/"
        lines.append(f"- `{name}`: {first_doc_line(path) or '(no docstring)'}")

    steps = module_constant(RUN_PY, "STAGES") or []
    lines += ["", f"**CLI `python -m cmlpipe run <stage>`**, stages: "
                  + (", ".join(f"`{s}`" for s in steps) or "none")]

    # Count test functions per file by parsing, so the list is right even when tests are not run.
    lines += ["", "**Tests `tests/`**:", ""]
    total = 0
    for path in sorted(TESTS.glob("test_*.py")):
        n = sum(isinstance(n, ast.FunctionDef) and n.name.startswith("test_") for n in parse(path).body)
        total += n
        lines.append(f"- `{path.name}`: {n} tests. {first_doc_line(path)}".rstrip())
    lines += ["", f"Test result: {test_result(run_tests)} ({total} test functions in total)."]
    return "\n".join(lines)


def test_result(run_tests: bool) -> str:
    """Run pytest and return a one-line summary such as '25 passed in 7.4s'."""
    if not run_tests:
        return "not run"
    try:
        out = subprocess.run([sys.executable, "-m", "pytest", str(TESTS), "-q", "-p", "no:cacheprovider"],
                             cwd=REPO, capture_output=True, text=True, timeout=600)
    except (OSError, subprocess.TimeoutExpired) as err:
        return f"could not run ({type(err).__name__})"
    # pytest's last summary line looks like "25 passed, 3 skipped in 7.36s"
    summary = [l for l in out.stdout.splitlines() if re.search(r"\d+ (passed|failed|error)", l)]
    if not summary:
        return f"could not run with `{Path(sys.executable).name}` (missing packages?)"
    # drop the run time ("in 7.36s") so the README does not change just because tests ran faster
    text = re.sub(r"\s+in [\d.]+s.*$", "", summary[-1].strip("= ").strip())
    mark = "PASS" if out.returncode == 0 else "FAIL"
    return f"**{mark}** {text}, on {dt.date.today():%Y-%m-%d}"


def block_data() -> str:
    lines = []
    dims = module_constant(SCHEMA_PY, "DIMS")
    data_vars = module_constant(SCHEMA_PY, "DATA_VARS") or {}
    if dims:
        lines.append(f"L1 CML layout (from `src/cmlpipe/io/schema.py`): dims `{tuple(dims)}`; data variables "
                     + ", ".join(f"`{k}` [{u}]" for k, (u, _) in data_vars.items()) + ".")
        lines.append("")

    # Processed data is not in git, so what follows describes the computer that made the last commit.
    files = sorted(L1_CML_DIR.glob("cml_10s_*.nc")) if L1_CML_DIR.exists() else []
    if files:
        days = [f.stem.removeprefix("cml_10s_") for f in files]
        size_gb = sum(f.stat().st_size for f in files) / 1e9
        lines.append(f"- L1 CML files: **{len(files)} days** ({days[0]} to {days[-1]}), {size_gb:.2f} GB, "
                     f"in `{rel(L1_CML_DIR)}/`")
    else:
        lines.append(f"- L1 CML files: none found in `{rel(L1_CML_DIR)}/`; run the ingest command above.")

    if MANIFEST.exists():
        try:
            last = json.loads(MANIFEST.read_text(encoding="utf-8"))[-1]
            lines.append(f"- Last pipeline run: {last.get('finished_utc', '?')}, steps "
                         f"{', '.join(last.get('steps', []))}, {len(last.get('days', []))} days, "
                         f"status **{last.get('status', '?')}**")
        except (ValueError, IndexError, KeyError):
            lines.append("- Last pipeline run: manifest could not be read")
    lines.append("")
    lines.append("*(Processed data is not in git. These lines describe the computer that made the last commit.)*")
    return "\n".join(lines)


def block_progress() -> str:
    if not ROADMAP.exists():
        return f"`{rel(ROADMAP)}` is missing."
    sections: list[tuple[str, list[tuple[bool, str]]]] = []
    for line in ROADMAP.read_text(encoding="utf-8").splitlines():
        if line.startswith("## "):
            sections.append((line[3:].strip(), []))
        elif (m := re.match(r"\s*- \[([ xX])\] (.+)", line)) and sections:
            sections[-1][1].append((m.group(1) != " ", m.group(2).strip()))

    done = sum(d for _, items in sections for d, _ in items)
    total = sum(len(items) for _, items in sections)
    lines = [f"Overall: **{done} / {total}** items done. Edit `{rel(ROADMAP)}` to tick items.", "",
             "| Part | Done |", "|---|---|"]
    for title, items in sections:
        n = sum(d for d, _ in items)
        lines.append(f"| {title} | {n} / {len(items)}{' ✔' if items and n == len(items) else ''} |")

    lines += ["", "**Next up** (first open items):", ""]
    todo = [(title, text) for title, items in sections for d, text in items if not d]
    lines += [f"- [ ] {text} *({title})*" for title, text in todo[:5]] or ["- nothing left!"]

    lines += ["", "<details><summary>Full checklist</summary>", ""]
    for title, items in sections:
        lines += [f"**{title}**", ""] + [f"- [{'x' if d else ' '}] {t}" for d, t in items] + [""]
    lines.append("</details>")
    return "\n".join(lines)


# ---------------------------------------------------------------- main

def replace_block(text: str, name: str, make_body) -> str:
    """Replace the block `name`; make_body() is only called when the README has that block (no wasted test runs)."""
    pattern = re.compile(rf"(<!-- AUTO:{name} -->).*?(<!-- /AUTO:{name} -->)", re.S)
    if not pattern.search(text):
        return text
    body = make_body()
    return pattern.sub(lambda m: f"{m.group(1)}\n{body}\n{m.group(2)}", text)


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--no-tests", action="store_true", help="do not run pytest")
    args = p.parse_args(argv)

    text = README.read_text(encoding="utf-8")
    new = text
    new = replace_block(new, "updated", lambda: f"*Status blocks last refreshed: {dt.date.today():%Y-%m-%d}*")
    new = replace_block(new, "implemented", lambda: block_implemented(run_tests=not args.no_tests))
    new = replace_block(new, "data", block_data)
    new = replace_block(new, "progress", block_progress)
    if new != text:
        README.write_text(new, encoding="utf-8", newline="\n")
        print("update_readme: README.md refreshed")
    else:
        print("update_readme: README.md already up to date")
    return 0


if __name__ == "__main__":
    sys.exit(main())
