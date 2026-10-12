---
description: Bring docs/research_log.pdf, docs/project_structure.pdf and the reproduce page up to date with the project
argument-hint: "[full] [commit]   (full = recheck every section; commit = commit the result)"
---

# /updatedocs

Arguments: `$ARGUMENTS`

Bring the project's documents up to date with what the research has actually done:

1. `docs/research_log.tex` -> `docs/research_log.pdf`: the report of what the research has implemented, what it
   found, where it stands and what comes next, at the level of ideas (equations and plain language, no code).
2. `docs/project_structure.tex` -> `docs/project_structure.pdf`: which folders and files the repository has and what
   each one is for (short, mostly tables).
3. `docs/REPRODUCE.md`: the step table, generated from `research/steps.yaml`.

Do not ask questions; work through the steps and report at the end. Never change code, settings or results here.

## Steps

1. **The reproduce page.** Run `C:\Users\yoni\miniconda3\python.exe -m research docs`. If a research script with a
   `main()` is not in `research/steps.yaml` (the test `pipeline/tests/test_research_steps.py` says which), add its
   entry (name, module, what, needs, time from its run log or report) and run the command again.

2. **What changed.** Read the first line of `docs/research_log.tex`: `% last updated: <date> commit=<hash>`.
   List what changed since then: `git log --oneline <hash>..HEAD -- research/ docs/ROADMAP.md`,
   `git diff --stat <hash> -- research/ docs/ROADMAP.md`, and `git status --porcelain -- research/ docs/`
   (uncommitted work counts). If nothing changed and the argument `full` was not given, say so and skip steps 3
   to 6 (step 7, the project structure, has its own check and is always done).

3. **Read the sources.** Every number and every statement in the log must come from these:
   - `docs/ROADMAP.md`: part A (checklist) and part B (the WAA plan: status, decisions, next steps);
   - each module's `README.md` (status table) and `METHOD.md` (rules, what each step decided, open issues);
   - the step reports `research/*/outputs/reports/*.md` of the steps that changed (they hold the tables and the
     verdicts), and the tables in `outputs/tables/` when a report is not enough;
   - `research/config.yaml` for the settings and decision rules a step used;
   - the run logs in `research/.cache/<module>/logs/` for run times and for what a run printed (they are local:
     if they are missing, use the reports);
   - the git log, for the order and the dates of the work.

4. **Edit `docs/research_log.tex`** so that it describes the project as it is now. Only touch what changed: the
   status figure and table (Section 1), the section of every step that was finished, changed or started, "Where we
   stand", "What comes next", the open questions, the abstract and the date on the title page. A step that only
   progressed gets its section updated in place; a new step or module gets a new section in the same style. With
   `full`, check every section against the sources.

   Writing rules for this report:
   - ideas, not code: what a step does, why it is done that way, what it found; no file, function or setting names
     in the text (the appendix lists where the details are; keep its file names right);
   - the proposal's notation (A(t), A_rain, δ(t), R̄, L, f(L, R̄, Z)), equations where they explain, plain language
     for a final-year student;
   - state plainly what is finished, in progress and not started; a result on simulated data is called that; a step
     that failed is reported as failed; never present a planned step or an expected result as done;
   - figures: draw diagrams in TikZ; an existing figure may be included from `research/<module>/outputs/figures/`
     (they are committed, so the report builds from a clean clone); no figure may rely on colour alone (line styles,
     hatching, text labels);
   - references are the proposal's, with its numbers; never invent a citation.

5. **Build.** `C:\Users\yoni\miniconda3\python.exe docs/build_research_log.py`. It fails on a LaTeX error, an
   overfull box or an undefined reference: fix the `.tex`, do not ignore it. Then Read the pages of the PDF that
   changed and check for overlapping labels and misplaced figures.

6. **Stamp.** Rewrite the first line of the `.tex`: `% last updated: <today> commit=<git rev-parse --short HEAD>`,
   and rebuild if that changed the file after the last build.

7. **The project structure.** `docs/project_structure.tex` has the same first line as the research log
   (`% last updated: <date> commit=<hash>`). List the files that were added, deleted or renamed since then:
   `git diff --name-status --diff-filter=ADR <hash> HEAD -- . ":(exclude)research/*/outputs"` and
   `git status --porcelain` (uncommitted files count); also look at the local folders that are not in git
   (`data/`, `research/.cache/`, `research/explorations/`) one level deep. If nothing was added, deleted, renamed
   or given a new job, and `full` was not given, leave the document as it is. Otherwise:
   - edit only the rows that changed: a new file gets a row in the table of its folder, a deleted file loses its
     row, a new folder or module gets a table (or a subsection) in the same style; with `full`, check every row
     against `git ls-files` and against the first lines (the docstring) of each script;
   - what a file does comes from its docstring, its module's `README.md` or `research/steps.yaml`; never guess;
   - keep it short: one line per file, simple words, tables where there are several files, no numbers that go out
     of date (file counts, run times) and no research results (those belong in the research log);
   - file and folder names are written with `\code{...}` (it removes spaces: a command with spaces is `\texttt{...}`);
   - update the date on the title page, build with
     `C:\Users\yoni\miniconda3\python.exe docs/build_research_log.py project_structure` (same failure rules as in
     step 5), Read the PDF, then rewrite the first line of the `.tex` as in step 6 and rebuild.

8. **Report in the chat**: which sections changed and why, which numbers were added (with their source file), what
   could not be confirmed from the sources, and whether the builds passed. Leave the changed files uncommitted for
   review, unless the argument `commit` was given: then commit `docs/research_log.tex`, `docs/research_log.pdf`,
   `docs/project_structure.tex`, `docs/project_structure.pdf`, `docs/REPRODUCE.md` and `research/steps.yaml` with a
   message that says what was updated.

**Reminder to give in the report:** if every other item of "Stage 1: WAA" in part A of `docs/ROADMAP.md` is ticked
and the item "publish the fit caches" is not, tell the user that it is time to publish them (the item says how).

`docs/pipeline_log.pdf` (the data-layer report) has no LaTeX source in the repository and is not rebuilt. If the
pipeline changed in a way that makes it wrong, say so in the report.
