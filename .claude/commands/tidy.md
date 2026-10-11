---
description: Tidy the repository (or one path) according to ORGANIZATION.md, autonomously, on a branch; never deletes; updates the research log
argument-hint: "[path]  (empty = the whole repository, e.g. research/baseline_estm)"
---

# /tidy

Scope: `$ARGUMENTS` (if empty: the whole repository).

You tidy this scope **autonomously**. Do not ask the user organization questions and do not use AskUserQuestion:
decide yourself with the priorities in `ORGANIZATION.md` section 9, log every decision with a one-line reason, and
keep going. Only questions that depend on the **science** (e.g. "is the old method still needed?") are not decided:
keep both versions, cleanly named, and list the question under "Open questions" with a recommended answer.

**A tidy reorganizes; it never deletes content.** Files are moved (to their right place, to `archive/`, to a module's
local `legacy/` folder, or to `_review/`), never removed. If the only way to fix something is to delete a file, do not
fix it: report it.

## Steps

1. **Read the rules**: re-read `ORGANIZATION.md` and `CLAUDE.md` (they may have changed since the last tidy), and
   the scope's `README.md` / `CATALOG.md` / `METHOD.md` if there are any.
2. **Run the checker**: `C:\Users\yoni\miniconda3\python.exe pipeline/tools/check_organization.py`. Keep only the
   problems inside the scope.
3. **Audit the scope** beyond the checker, from the code (grep for `to_csv`, `savefig`, `to_parquet`, `np.savez`,
   `write_text`, `open(..., "w")`, `to_netcdf`, zarr writes): misplaced files, outputs without a CATALOG entry or
   whose producer no longer writes that path (orphans), duplicate or stale docs, `_v2`-style names, superseded
   outputs, stray caches, hard-coded output paths, broken path references, nested `.gitignore` files, untracked
   files that look like junk.
4. **Branch**: `git switch -c tidy/<scope-or-repo>-<YYYY-MM-DD>` (never work on `main`, never merge).
5. **Fix**, in small commits with clear messages:
   - tracked files: `git mv`; untracked files: plain move;
   - update every reference (code, configs, docs, notebooks), then grep for each old path; none may remain
     outside `archive/`, `docs/organization/` and the decision log;
   - output folders come from the module's `<module>_paths.py`; this and path strings are the only code changes
     allowed; never change a computation, parameter or result;
   - superseded files that nothing reads: move to `archive/<what>/` with a `NOTE.md` (what it was, what replaced
     it, why, when); keep active (named by content) if code still reads them;
   - untracked files that look like junk: move to `_review/<original path>` (the user decides);
   - add or fix CATALOG entries and READMEs (template: `ORGANIZATION.md` section 6).
6. **Outputs: important vs legacy** (rules: `ORGANIZATION.md` section 5; read them there, they are not repeated
   here). Important outputs (`key`, `product`) stay in `<module>/outputs/` and are committed; the others (`diag`,
   `scratch`, unlisted) are moved to `<module>/legacy/`, which is git-ignored, so they stay on this machine and are
   not pushed. The script applies the rules:
   a. Dry run: `C:\Users\yoni\miniconda3\python.exe research/common/tidy_outputs.py` (only modules touched since
      their last tidy; `--all` for every module). Read what it would move, keep, and report as conflicts.
      A script that still writes without `save_output` is a finding: refactor it first (a level for each output;
      `diag` for step figures and reports, `product` for frozen data products, at most 6 `key`), or its outputs
      move to legacy as "not in MANIFEST.csv".
   b. The protections (read by a script, named in `DECISIONS.md`, linked from a Markdown file that stays) are applied
      by the script. A conflict keeps the file in `outputs/`: report it, never relabel a file just to make a conflict
      go away.
   c. For every step report in the move list that the module's `DECISIONS.md` does not cover yet (the dry run lists
      them), add 5-10 lines there: what was decided, the numbers behind it, what stays open, naming the script that
      wrote the report (not the report file: naming it would keep it in `outputs/`). `DECISIONS.md` is the committed
      summary of the steps and a source of the research log (step 7).
   d. Apply: `... research/common/tidy_outputs.py --apply`. It moves each file to `legacy/<kind>/<name>` (a tracked
      file is first taken out of git with `git rm --cached`, the file itself stays), logs every move to
      `<module>/outputs/TIDY_LOG.md`, moves the rows from `outputs/MANIFEST.csv` to `legacy/MANIFEST.csv` and
      rewrites the `.gitignore` block. Fix README / CATALOG text that presents a moved file as being in `outputs/`
      (a CATALOG entry may still describe what a script writes). Commit.
   e. `check_organization.py` must pass (it also checks MANIFEST coverage, <= 6 key per module, no version suffix in
      output names, no tracked diag / scratch file, nothing tracked under `legacy/`).
   f. Print a short summary: moved to legacy (count, size) per module, the kept key and product files, the conflicts.
7. **Update the research log** (`docs/research_log.tex` -> `docs/research_log.pdf`): the report of what the research
   has implemented so far, at the level of ideas. Do this on every tidy of the whole repository or of anything under
   `research/`, whenever the research changed since the report was last updated:
   a. Read the first line of `docs/research_log.tex`: `% last updated: <date> commit=<hash>`. List what changed since
      then: `git diff --stat <hash> -- research/` plus `git status --porcelain -- research/` (uncommitted work
      counts). If nothing under `research/` changed (code, `ROADMAP.md`, a module's README, `DECISIONS.md`,
      `roadmap.md`, key outputs), skip this step and say so in the report.
   b. Re-read the sources: `research/README.md`, `research/ROADMAP.md`, and for every module its README (status
      table, open issues), `METHOD.md`, `DECISIONS.md`, `roadmap.md` and the key reports in `outputs/reports/`
      (step reports that moved are in `legacy/reports/`). Every number in the report must come from these.
   c. Edit the `.tex` so that it describes the project as it is now: the status figure and table (Section 1), the
      section of every step that was finished, changed or started, "Where we stand", "What comes next", the open
      questions, the abstract, the date on the title page, and the `% last updated` line (today, `git rev-parse
      --short HEAD`). A step that only progressed gets its section updated in place; a new step or module gets a new
      section in the same style. Do not rewrite sections whose content did not change.
   d. Writing rules for this report:
      - ideas, not code: what a step does, why it is done that way, what it found; no file, function or setting
        names in the text (Appendix B lists where the details are);
      - the proposal's notation (A(t), A_rain, δ(t), R̄, L, f(L, R̄, Z)), equations where they explain, plain language
        for a final-year student;
      - state plainly what is finished, in progress and not started; a result on simulated data is called that;
        never present a planned step or an expected result as done;
      - figures: draw diagrams in TikZ; include an existing figure only if it is a `key` output (committed, so the
        report builds from a clean clone); no figure may rely on colour alone (line styles, hatching, text labels);
      - references are the proposal's, with its numbers; never invent a citation.
   e. Build: `C:\Users\yoni\miniconda3\python.exe pipeline/tools/build_research_log.py` (it fails on a LaTeX error,
      an overfull box or an undefined reference: fix the `.tex`, do not ignore it). Then look at the rendered pages
      that changed (Read the PDF) for overlapping labels and misplaced floats. Commit the `.tex` and the `.pdf`.
   `docs/pipeline_log.pdf` (the data-layer report) has no LaTeX source in the repository: it is not rebuilt. If the
   pipeline changed in a way that makes it wrong, say so under "Open questions".
8. **Verify** after each part:
   - `C:\Users\yoni\miniconda3\python.exe -m pytest -q` from `pipeline/` (note which tests need the data);
   - `python pipeline/tools/check_boundaries.py` and `python pipeline/tools/check_organization.py`;
   - re-run at least one cheap producer of each touched module and compare the new files with the committed ones
     (`git status` must show no content change; for `.gz` compare the decompressed content). If a re-run gives
     different numbers, restore the committed bytes and report it; do not commit changed results.
   - if a step breaks something: fix the paths or revert that commit, and log why.
9. **Report** (in the chat, and in `docs/organization/TIDY_<YYYY-MM-DD>.md` if anything changed): what moved, what
   went to `legacy/`, `archive/` and `_review/`, what changed in the research log, the decision log, verification
   results, open questions (plain words, with a recommended answer), suspicious things noticed but not changed.

## Never

- delete a file, tracked or untracked (move it instead; if in doubt, `_review/`);
- touch `data/`, `research/.vibe-wise/`, or the files inside `research/explorations/` (a README there is fine);
- change the frozen pipeline tables (`events.csv`, `split_*.csv`, `event_links.csv`, `dry_periods.csv`);
- force-push, rewrite history, merge into `main`, or skip the git hook.
