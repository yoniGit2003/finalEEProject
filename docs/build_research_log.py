"""Build docs/research_log.pdf from docs/research_log.tex (the research report at the level of ideas).

    python docs/build_research_log.py                      # the research log
    python docs/build_research_log.py project_structure    # another document of docs/, here project_structure.tex
    python docs/build_research_log.py WAA_research_overview # the WAA overview (figures in WAA_research_overview_figures/)

Runs pdflatex three times (the second run fills the table of contents and the references; the contents pages
shift every later page, so the third run corrects the page numbers in them). The auxiliary files (.aux, .log,
.toc, .out) go to a temporary folder, so only the PDF is written into docs/. Needs a LaTeX installation with pdflatex
(MiKTeX or TeX Live) and the packages tikz and pgfplots. Used by the /updatedocs command after it has updated the .tex.
Exit code 1 if pdflatex fails or reports an overfull box / an undefined reference (the lines are printed).
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

DOCS = Path(__file__).resolve().parent
DOCUMENTS = ("research_log", "project_structure", "WAA_research_overview")   # the .tex files of docs/ that this script builds


def main() -> int:
    NAME = sys.argv[1] if len(sys.argv) > 1 else DOCUMENTS[0]
    if NAME not in DOCUMENTS:
        print(f"build_research_log: unknown document {NAME!r} (choose from: {', '.join(DOCUMENTS)})")
        return 1
    if shutil.which("pdflatex") is None:
        print("build_research_log: pdflatex not found (install MiKTeX or TeX Live)")
        return 1
    where = DOCS / Path(NAME).parent                           # the folder of the .tex (figure paths are relative to it)
    NAME = Path(NAME).name
    with tempfile.TemporaryDirectory() as tmp:
        for run in (1, 2, 3):
            # cwd = docs/, so that the figure paths in the .tex (../research/<module>/outputs/...) resolve
            r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"-output-directory={tmp}",
                                f"{NAME}.tex"], cwd=where, capture_output=True, text=True, errors="replace")
            pdf = Path(tmp) / f"{NAME}.pdf"
            if not pdf.exists():                               # pdflatex's exit code is not reliable on old Windows
                print(r.stdout[-3000:])
                print(f"build_research_log: pdflatex failed (run {run})")
                return 1
        log = (Path(tmp) / f"{NAME}.log").read_text(encoding="utf-8", errors="replace")
        problems = [line for line in log.splitlines() if re.search(r"Overfull \\hbox|undefined|multiply defined", line)]
        shutil.copyfile(pdf, where / f"{NAME}.pdf")
    pages = re.search(r"Output written on .*?\((\d+) pages?", log, re.S)
    print(f"build_research_log: wrote {(where / NAME).relative_to(DOCS.parent).as_posix()}.pdf ({pages.group(1) if pages else '?'} pages)")
    for p in problems:
        print("  " + p)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
