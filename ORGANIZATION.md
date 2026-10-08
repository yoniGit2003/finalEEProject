# ORGANIZATION: where things go (the rulebook)

Short on purpose. Checked by `python pipeline/tools/check_organization.py`; applied by the `/tidy` command.

## 1. Layout

```
README.md  ORGANIZATION.md  CATALOG.md  CLAUDE.md  .gitignore (the only one)
docs/              project documents: proposal, reports, thesis, their figures; docs/organization/ = reorg records
pipeline/          the data layer: src/cmlpipe (package), tests/ (pytest), apps/, configs/, docs/, tools/, examples/
research/          config.yaml (all research settings), ROADMAP.md, README.md
  <module>/        one folder per research step (baseline_estm, WAA_estm, later path_law, evaluation)
  common/          shared research helpers: outputs.py (save_output), tidy_outputs.py (the output part of /tidy)
  explorations/    the playground
data/              raw / work / release data: NOT in git, never modified by hand          (protected)
_review/           git-ignored: files /tidy would have deleted but git cannot restore (you decide)
```

**A research module** always looks like this:

```
<module>/
  README.md              the only README of the module (template in section 6)
  METHOD.md              the scientific method, if longer than a README section
  CATALOG.md             one entry per output file
  DECISIONS.md           the decisions taken in the step reports (kept when /tidy deletes those reports)
  <module>_paths.py      every folder location of the module (and files other scripts read)
  *.py                   library code and entry scripts
  validation/            scientific checks: scripts only, they write to outputs/
  outputs/products/      deliverables other modules read
  outputs/tables/        other tables
  outputs/figures/       figures
  outputs/reports/       written reports (.md)
  outputs/MANIFEST.csv   level + purpose of every output file (written by save_output; section 5)
  outputs/TIDY_LOG.md    what /tidy deleted, and which script rebuilds it
  archive/               superseded outputs nothing reads; each subfolder has a NOTE.md
  .cache/                git-ignored: loaded data, DEV runs, logs; safe to delete
```

**Kinds stay apart**: source code, scientific validation (`<module>/validation/`), code tests (`pipeline/tests/`),
configuration (`research/config.yaml`, `pipeline/configs/`), generated outputs (`outputs/`, `data/`), documentation
(`README`, `METHOD`, `CATALOG`, `docs/`), caches (`.cache/`).

## 2. Naming

- Lower case with underscores; say what it is: `dew_risk_hourly.csv`, not `table3.csv`.
- Every output file name starts with the script or check that writes it: `c1_...`, `r1b_...`, `sample_two_window_...`.
- **No `_v2`, `_new`, `_final`, `_old` on active paths.** One current version per output; the old one goes to
  `archive/` (if untracked) or is deleted (if tracked: git keeps it).
- Two **method** versions that are both in use are named by what differs (`sample_two_window_*` vs
  `sample_pre_event_*`), and README / CATALOG say which method id (v1, v2 in `config.yaml`) each one is.
- Validation checks: `c<n>_<topic>.py`; a re-run of a whole validation: `r<n><letter>_<topic>.py` (e.g. `r1b_...`).

## 3. Where new things go

| new thing | goes to |
|---|---|
| script of an existing module | `research/<module>/`; output folders from `<module>_paths.py`; settings from `research/config.yaml`; every output written with `save_output` (section 5) |
| scientific check | `research/<module>/validation/` (writes to the module's `outputs/`) |
| code test (does the code work?) | `pipeline/tests/` (pytest; real-data tests marked `-m data`) |
| table / figure / report | `outputs/tables/`, `outputs/figures/`, `outputs/reports/` of the module that makes it, plus a CATALOG entry |
| doc about the project | `docs/`; about the data layer: `pipeline/docs/`; about a module: its README / METHOD |
| research module | `research/<name>/` with the layout above, a `config.yaml` section and a row in the root `CATALOG.md` |
| helper shared by research modules | `research/common/` (not a module: no README / CATALOG / paths file) |
| repo tool or checker | `pipeline/tools/` |
| an idea to try | `research/explorations/` |

**Graduating an exploration**: copy the script into a module, make it read locations from the module's paths file,
write to `outputs/`, add CATALOG entries. Modules never read files from `explorations/`.

## 4. Git: what is committed

- **Committed**: code, configs, docs, `MANIFEST.csv`, `TIDY_LOG.md`, and the outputs with level **key** or
  **product** that are <= 10 MB (section 5).
- **Ignored** (root `.gitignore` only): `data/`, `.cache/`, `_cache/`, `__pycache__/`, `.pytest_cache/`, `.idea/`,
  `_archive/`, `_review/`, `research/.vibe-wise/`, DEV runs, every **diag** / **scratch** output, anything > 10 MB.
  The `.gitignore` ignores all of `research/*/outputs/` and re-includes the key / product files in a generated block
  that `save_output` rewrites from the MANIFESTs: never edit that block by hand.
- A nested `.gitignore` is not allowed (tools such as PyCharm and pytest write their own inside ignored folders: fine).
- PySR runs set `output_directory` to the module's `.cache/` (the root `/outputs/` rule only covers the root).

## 5. Outputs and importance

These rules live **only here**; `research/common/tidy_outputs.py` and `/tidy` apply them, `check_organization.py`
checks them. Importance is set **when an output is written**: every file a research script writes into
`<module>/outputs/` goes through

```python
from research.common.outputs import save_output      # importable once the module's <module>_paths.py is imported
save_output(obj_or_fig, path, level, purpose, **writer_kwargs)
```

which writes the file (figure, DataFrame, pyarrow Table or text) and upserts its row in `outputs/MANIFEST.csv`
(`path`, `level`, `purpose`, `script`, `commit`, `date`). A file written another way (a streamed Parquet file, a
report written by hand) is recorded with `register_output(path, level, purpose)` (`script = "by hand"` for the
latter).

| level | what | default for | committed | deleted by /tidy |
|---|---|---|---|---|
| `key` | what the final report shows; **at most 6 per module** | (chosen) | yes | never |
| `product` | data a script or another module reads; frozen data products | frozen products, tables other scripts read | yes if <= 10 MB | only if no script reads **and** no script writes it any more |
| `diag` | step figures, step reports, checks, intermediate tables | **every step figure and step report** | never | yes |
| `scratch` | throw-away | | never | yes |

**What /tidy deletes** (for every module touched since its last tidy): outputs with level `diag` or `scratch`,
outputs not in `MANIFEST.csv`, and `product` files that no script reads and no script writes any more. A product that
its producer still writes but nothing reads yet (e.g. a frozen set waiting for the next module) is kept and listed.

**What /tidy never deletes**: code, configs, `README`, `CATALOG`, `METHOD`, `DECISIONS.md`, `MANIFEST.csv`,
`TIDY_LOG.md`, anything under `data/` or `archive/`, a file that a script reads (its name, or the paths-file constant
holding it, appears in the code of a script other than its producer; docstrings and comments do not count), a file
named in the module's `DECISIONS.md`, and a file linked from a Markdown file that stays. If the level says "delete"
but one of these applies, the file is kept and reported as a **conflict**.

**Before deleting**, the decisions in the step reports that are about to go are written to the module's
`DECISIONS.md`: 5-10 lines per step, naming the script that wrote the report (`tidy_outputs.py --apply` refuses
otherwise). Every deletion is logged in `outputs/TIDY_LOG.md` (path, level, size, the script that rebuilds it).
Deleting is safe because diag / scratch outputs are rebuilt by rerunning their script.

**Checks** (`check_organization.py`): every file in `outputs/` has a `MANIFEST.csv` row; at most 6 key per module;
no version suffix (`_v2`, `_new`, `_final`, `_old`) in output names; no diag / scratch file is tracked by git.

## 6. Catalog and README template

Every output file (or group of same-kind files) has a CATALOG entry: path, plain description, columns with units
(tables) or what is plotted and what to look for (figures), the question it answers, producer + command, inputs,
status (current / archived). Module README template:

```markdown
# <module>: <one line: what it computes>
## Purpose        why this step exists, in the proposal's notation
## Status         what is done / frozen (date, commit)
## How to run     the commands in order, with run times and what each needs
## Code           table: file -> what it does
## Outputs        the four outputs/ folders and a link to CATALOG.md
## Method         link to METHOD.md (or the method itself, if short)
## Open issues    what is undecided
```

## 7. Archive

`archive/<what>/` holds superseded outputs that **no code reads**, with a `NOTE.md`: what it was, what replaced it,
why, and the date. If code still reads an old version, it is not archived: it stays active, named by content.

## 8. Protected paths (never moved, edited or deleted by a tidy)

`data/` (raw data and pipeline outputs), `research/.vibe-wise/` (made by a plugin; ignored, never committed),
the files inside `research/explorations/` (README only), frozen pipeline tables (`events.csv`, `split_*.csv`,
`event_links.csv`, `dry_periods.csv`).

## 9. How Claude tidies (`/tidy [path]`)

Claude decides organization questions itself, using these priorities: (1) nothing breaks, (2) reversible over
irreversible, (3) one source of truth, (4) kinds stay apart, (5) the same layout in every module, (6) names describe
content, (7) regenerable means ignored, (8) fewer, shallower folders, (9) keep existing conventions that work. It works
on a branch, uses `git mv`, commits in small steps, logs each decision with a one-line reason, and **never deletes an
untracked file** (those go to `_review/`), with one exception: outputs inside `<module>/outputs/` that section 5 says
to delete are deleted (tracked or not) by `research/common/tidy_outputs.py --apply`, which logs them to
`TIDY_LOG.md`; their script rebuilds them. It never changes scientific logic, parameters or results. Questions that
depend on the science (is an old method still needed?) are not decided: both versions are kept, cleanly named, and
the question goes to the report's "Open questions". Claude does not stop to ask, except when a decision in a step
report about to be deleted is unclear (section 5, DECISIONS.md).
