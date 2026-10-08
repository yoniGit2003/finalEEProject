---
description: Tidy the repository (or one path) according to ORGANIZATION.md, autonomously, on a branch
argument-hint: "[path]  (empty = the whole repository, e.g. research/baseline_estm)"
---

# /tidy

Scope: `$ARGUMENTS` (if empty: the whole repository).

You tidy this scope **autonomously**. Do not ask the user organization questions and do not use AskUserQuestion:
decide yourself with the priorities in `ORGANIZATION.md` section 9, log every decision with a one-line reason, and
keep going. Only questions that depend on the **science** (e.g. "is the old method still needed?") are not decided:
keep both versions, cleanly named, and list the question under "Open questions" with a recommended answer.
The one exception: a decision in a step report about to be deleted that you cannot state clearly (step 6c).

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
   - superseded outputs: delete if tracked (git keeps them) and nothing reads them; archive (`archive/<what>/`
     + `NOTE.md`) if untracked or still cited by a report; keep active (named by content) if code still reads them;
   - **never delete an untracked file**: move it to `_review/<original path>` (exception: step 6, outputs);
   - add or fix CATALOG entries and READMEs (template: `ORGANIZATION.md` section 6).
6. **Outputs and importance** (rules: `ORGANIZATION.md` section 5; read them there, they are not repeated here).
   You have the authority to delete outputs; the script applies the rules:
   a. Dry run: `C:\Users\yoni\miniconda3\python.exe research/common/tidy_outputs.py` (only modules touched since
      their last tidy; `--all` for every module). Read what it would delete, keep, and report as conflicts.
      A script that still writes without `save_output` is a finding: refactor it first (a level for each output;
      `diag` for step figures and reports, `product` for frozen data products, at most 6 `key`), or its outputs are
      deleted as "not in MANIFEST.csv".
   b. The protections (read by a script, named in `DECISIONS.md`, linked from a Markdown file that stays) are applied
      by the script. A conflict keeps the file: report it, never relabel a file just to make a conflict go away.
   c. Before deleting: for every step report in the delete list, add 5-10 lines to the module's `DECISIONS.md` (what
      was decided, the numbers behind it, what stays open), naming the script that wrote the report, not the report
      file (naming it would protect it). Ask the user only if a decision is unclear.
   d. Apply: `... research/common/tidy_outputs.py --apply`. It deletes (tracked: `git rm`; untracked: removed), logs
      every deletion to `<module>/outputs/TIDY_LOG.md` (path, level, size, the script that rebuilds it), drops the
      rows from `MANIFEST.csv`, un-tracks kept diag / scratch files and rewrites the `.gitignore` block. Fix README /
      CATALOG text that presents a deleted file as kept (a CATALOG entry may still describe what a script writes).
      Commit.
   e. `check_organization.py` must pass (it also checks MANIFEST coverage, <= 6 key per module, no version suffix in
      output names, no tracked diag / scratch file).
   f. Print a short summary: deleted (count, size) per module, the kept key and product files, the conflicts.
7. **Verify** after each part:
   - `C:\Users\yoni\miniconda3\python.exe -m pytest -q` from `pipeline/` (note which tests need the data);
   - `python pipeline/tools/check_boundaries.py` and `python pipeline/tools/check_organization.py`;
   - re-run at least one cheap producer of each touched module and compare the new files with the committed ones
     (`git status` must show no content change; for `.gz` compare the decompressed content). If a re-run gives
     different numbers, restore the committed bytes and report it; do not commit changed results.
   - if a step breaks something: fix the paths or revert that commit, and log why.
8. **Report** (in the chat, and in `docs/organization/TIDY_<YYYY-MM-DD>.md` if anything changed): what moved, what
   was deleted, what went to `_review/`, the decision log, verification results, open questions (plain words, with a
   recommended answer), suspicious things noticed but not changed.

## Never

- touch `data/`, `research/.vibe-wise/`, or the files inside `research/explorations/` (a README there is fine);
- change the frozen pipeline tables (`events.csv`, `split_*.csv`, `event_links.csv`, `dry_periods.csv`);
- delete untracked files outside step 6, force-push, rewrite history, merge into `main`, or skip the git hook.
