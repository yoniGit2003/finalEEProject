"""Check the boundary between the data layer and the research code (exit 1 on a violation).

Rules:
  * a file under research/ may import from cmlpipe only cmlpipe.api and cmlpipe.viz (never other internals);
  * no file under pipeline/ (code, explorer, tests, tools) may import research: the dependency goes one way only.
Imports are read with ast (the files are never run).

    python pipeline/tools/check_boundaries.py
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]   # repository root (pipeline/tools/check_boundaries.py)
ALLOWED = ("cmlpipe.api", "cmlpipe.viz")


def imports(path: Path) -> list[tuple[int, str]]:
    """(line, module) of every import in a .py file. `import a.b` gives a.b; `from a import b` gives a.b (so
    `from cmlpipe import io` is seen as cmlpipe.io); a bare `import cmlpipe` gives cmlpipe."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out += [(node.lineno, a.name) for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            out += [(node.lineno, f"{node.module}.{a.name}") for a in node.names]
    return out


def _is(mod: str, pkg: str) -> bool:
    return mod == pkg or mod.startswith(pkg + ".")


def violations(repo: Path = REPO) -> list[str]:
    """Every boundary violation as 'file:line: message'."""
    bad = []
    for f in sorted((repo / "research").rglob("*.py")):
        for line, mod in imports(f):
            if _is(mod, "cmlpipe") and not any(_is(mod, a) for a in ALLOWED):
                bad.append(f"{f.relative_to(repo)}:{line}: research imports {mod} (allowed: {', '.join(ALLOWED)})")
    for f in sorted((repo / "pipeline").rglob("*.py")):
        for line, mod in imports(f):
            if _is(mod, "research"):
                bad.append(f"{f.relative_to(repo)}:{line}: pipeline imports {mod}")
    return bad


def main() -> int:
    bad = violations()
    for b in bad:
        print(b)
    n_r = len(list((REPO / "research").rglob("*.py")))
    n_s = len(list((REPO / "pipeline").rglob("*.py")))
    print(f"check_boundaries: {'FAIL' if bad else 'OK'} ({len(bad)} violations; {n_r} research files, "
          f"{n_s} pipeline files checked)")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
