# docs/: project documents

| file | what it is |
|---|---|
| [`REPRODUCE.md`](REPRODUCE.md) | how to rerun everything, command by command, with run times |
| [`ROADMAP.md`](ROADMAP.md) | the one roadmap: project checklist on top, then the detailed wet-antenna plan |
| [`research_log.pdf`](research_log.pdf) | the research log: what was done, what it found, what comes next, at the level of ideas (source `research_log.tex`) |
| [`project_structure.pdf`](project_structure.pdf) | the project structure: which folders and files there are and what each one is for (source `project_structure.tex`) |
| [`pipeline_log.pdf`](pipeline_log.pdf) | technical report of the data layer (no LaTeX source in the repository) |
| [`proposal.md`](proposal.md), `proposal.docx`, `media/` | the project proposal (Hebrew) and its figure |
| [`sindy_intro_report_2026-03.pdf`](sindy_intro_report_2026-03.pdf) | our report "Introduction to SINDy" (March 2026) |
| `papers/` | the OpenMRG paper (Andersson et al. 2022) |
| `build_research_log.py` | builds `research_log.pdf` from `research_log.tex`; with the argument `project_structure` it builds `project_structure.pdf` |

**Updating the research log.** In Claude Code, `/updatedocs` reads the roadmap, the modules' `README.md` and
`METHOD.md`, the step reports in `research/*/outputs/reports/`, the run logs and the git history since the log was
last updated; it edits `research_log.tex` where the research changed and rebuilds the PDF. By hand: edit the `.tex`,
then `python docs/build_research_log.py` (needs `pdflatex`, e.g. MiKTeX or TeX Live).

**Updating the project structure.** `/updatedocs` also checks which files were added, deleted or renamed since
`project_structure.tex` was last updated, edits the rows that changed and rebuilds the PDF. By hand: edit the `.tex`,
then `python docs/build_research_log.py project_structure`.

Documentation of the data and the pipeline is in [`../pipeline/docs/`](../pipeline/docs/) (start with
`DATA_GUIDE.md`); each research module documents itself in its `README.md` and `METHOD.md`.
