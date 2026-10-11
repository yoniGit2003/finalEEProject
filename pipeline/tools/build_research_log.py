"""Build docs/research_log.pdf from docs/research_log.tex (the research report at the level of ideas).

    python pipeline/tools/build_research_log.py

Runs pdflatex twice (the second run fills the table of contents and the references). The auxiliary files (.aux, .log,
.toc, .out) go to a temporary folder, so only the PDF is written into docs/. Needs a LaTeX installation with pdflatex
(MiKTeX or TeX Live) and the packages tikz and pgfplots. Used by the /tidy command after it has updated the .tex.
Exit code 1 if pdflatex fails or reports an overfull box / an undefined reference (the lines are printed).
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DOCS = REPO / "docs"
NAME = "research_log"


def main() -> int:
    if shutil.which("pdflatex") is None:
        print("build_research_log: pdflatex not found (install MiKTeX or TeX Live)")
        return 1
    with tempfile.TemporaryDirectory() as tmp:
        for run in (1, 2):
            # cwd = docs/, so that the figure paths in the .tex (../research/<module>/outputs/...) resolve
            r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"-output-directory={tmp}",
                                f"{NAME}.tex"], cwd=DOCS, capture_output=True, text=True, errors="replace")
            pdf = Path(tmp) / f"{NAME}.pdf"
            if not pdf.exists():                               # pdflatex's exit code is not reliable on old Windows
                print(r.stdout[-3000:])
                print(f"build_research_log: pdflatex failed (run {run})")
                return 1
        log = (Path(tmp) / f"{NAME}.log").read_text(encoding="utf-8", errors="replace")
        problems = [line for line in log.splitlines() if re.search(r"Overfull \\hbox|undefined|multiply defined", line)]
        shutil.copyfile(pdf, DOCS / f"{NAME}.pdf")
    pages = re.search(r"Output written on .*?\((\d+) pages?", log, re.S)
    print(f"build_research_log: wrote docs/{NAME}.pdf ({pages.group(1) if pages else '?'} pages)")
    for p in problems:
        print("  " + p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
