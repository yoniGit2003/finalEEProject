"""Refresh the generated parts of the three READMEs (root, pipeline/, research/).

Each README has hand-written text and generated blocks. A generated block sits between two markers:

    <!-- AUTO:name -->
    ... replaced by this script ...
    <!-- /AUTO:name -->

Only the text between the markers is replaced; everything else is left as it is. Blocks (see TARGETS):
    README.md            AUTO:status    one line per part, from the two roadmaps
    pipeline/README.md   AUTO:data      processed data on this computer; AUTO:progress from pipeline/docs/roadmap.md
    research/README.md   AUTO:progress  from research/ROADMAP.md
Progress is edited by hand ONLY in the two roadmap files (tick boxes, and the "Status:" line under the title).
The script runs on every commit through the git hook pipeline/tools/githooks/pre-commit (enable once per clone:
git config core.hooksPath pipeline/tools/githooks), and can also be run by hand:

    python pipeline/tools/update_readme.py

It uses only the Python standard library, so any Python 3.9+ can run it. Facts about the code are read
by parsing the source files (ast), not by importing them, so numpy/xarray are not needed here. Only the
test run needs the project's packages; it uses the same Python that runs this script.
"""

from __future__ import annotations

import argparse
import ast
import datetime as dt
import json
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]          # pipeline/ (this README, the code, the tests)
ROOT = REPO.parent                                   # repository root (data/)
README = REPO / "README.md"
ROADMAP = REPO / "docs" / "roadmap.md"                # pipeline roadmap
RESEARCH = ROOT / "research"
RESEARCH_ROADMAP = RESEARCH / "ROADMAP.md"
PACKAGE = REPO / "src" / "cmlpipe"
RUN_PY = PACKAGE / "cli.py"
SCHEMA_PY = PACKAGE / "io" / "schema.py"
TESTS = REPO / "tests"
L1_CML_DIR = ROOT / "data" / "work" / "openmrg" / "L1_standard" / "cml"
MANIFEST = ROOT / "data" / "work" / "openmrg" / "run_manifest.json"


# ---------------------------------------------------------------- small helpers

def rel(path: Path) -> str:
    """Path relative to the repository root, with forward slashes (as shown on GitHub)."""
    return path.relative_to(ROOT).as_posix()


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


def read_roadmap(roadmap: Path) -> tuple[str, list[tuple[str, list[tuple[bool, str]]]]]:
    """(hand-written 'Status:' line, [(section title, [(done, text), ...]), ...]) of a roadmap file."""
    status, sections = "", []
    for line in roadmap.read_text(encoding="utf-8").splitlines():
        if line.startswith("Status:"):
            status = line[len("Status:"):].strip()
        elif line.startswith("## "):
            sections.append((line[3:].strip(), []))
        elif (m := re.match(r"\s*- \[([ xX])\] (.+)", line)) and sections:
            sections[-1][1].append((m.group(1) != " ", m.group(2).strip()))
    return status, sections


def link(target: Path, base: Path) -> str:
    """Relative link from a README folder `base` to `target` (forward slashes)."""
    return Path(os.path.relpath(target, base)).as_posix()


def block_progress(roadmap: Path, base: Path, n_next: int = 6) -> str:
    """Status / next steps of ONE part, from its own roadmap (links relative to the README folder `base`)."""
    if not roadmap.exists():
        return f"`{rel(roadmap)}` is missing."
    status, sections = read_roadmap(roadmap)
    done = sum(d for _, items in sections for d, _ in items)
    total = sum(len(items) for _, items in sections)
    r = link(roadmap, base)
    lines = [f"Status: **{status or 'not set'}**; {done} / {total} roadmap items done (tick items in [`{r}`]({r})).", ""]
    if len(sections) > 1:
        lines += [f"- {title}: {sum(d for d, _ in items)} / {len(items)}" for title, items in sections] + [""]
    todo = [text for _, items in sections for d, text in items if not d]
    lines.append("Next steps:" if todo else "Next steps: none left.")
    lines += [f"- [ ] {t}" for t in todo[:n_next]]
    if len(todo) > n_next:
        lines.append(f"- ... and {len(todo) - n_next} more in the roadmap")
    return "\n".join(lines)


def block_status() -> str:
    """Root README: one line per part (status from each roadmap's 'Status:' line, counts from its boxes)."""
    out = []
    for name, roadmap, readme in (("Pipeline (data layer)", ROADMAP, README),
                                  ("Research", RESEARCH_ROADMAP, RESEARCH / "README.md")):
        status, sections = read_roadmap(roadmap) if roadmap.exists() else ("roadmap missing", [])
        done = sum(d for _, items in sections for d, _ in items)
        total = sum(len(items) for _, items in sections)
        r = link(readme, ROOT)
        out.append(f"- **{name}**: {status}; {done} / {total} roadmap items. Details: [`{r}`]({r})")
    return "\n".join(out)


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
    p.add_argument("--no-tests", action="store_true", help="do not run pytest (only used by an AUTO:implemented block)")
    args = p.parse_args(argv)
    targets = {
        ROOT / "README.md": {"status": block_status},
        README: {"updated": lambda: f"*Status blocks last refreshed: {dt.date.today():%Y-%m-%d}*",
                 "implemented": lambda: block_implemented(run_tests=not args.no_tests),
                 "data": block_data, "progress": lambda: block_progress(ROADMAP, REPO)},
        RESEARCH / "README.md": {"progress": lambda: block_progress(RESEARCH_ROADMAP, RESEARCH)},
    }
    for path, blocks in targets.items():
        if not path.exists():
            print(f"update_readme: {rel(path)} missing, skipped")
            continue
        text = path.read_text(encoding="utf-8")
        new = text
        for name, make in blocks.items():
            new = replace_block(new, name, make)
        if new != text:
            path.write_text(new, encoding="utf-8", newline="\n")
            print(f"update_readme: {rel(path)} refreshed")
        else:
            print(f"update_readme: {rel(path)} already up to date")
    return 0


if __name__ == "__main__":
    sys.exit(main())
