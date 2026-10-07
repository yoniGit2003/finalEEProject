# Reorganization plan

Based on [`AUDIT.md`](AUDIT.md). Branch `reorg/project-structure`; nothing is merged into `main`.

## Target tree

```
README.md  ORGANIZATION.md  CATALOG.md  CLAUDE.md  LICENSE  .gitignore  .gitattributes
.claude/commands/tidy.md        the /tidy command
docs/                           project documents (proposal, reports, figures)
  organization/                 AUDIT.md, PLAN.md, REPORT.md of this reorganization
pipeline/                       the data layer (unchanged layout: src/, tests/, apps/, configs/, docs/, tools/, examples/)
  tools/check_organization.py   NEW: rule checker used by /tidy
research/
  README.md  ROADMAP.md  config.yaml
  explorations/                 playground (files untouched; README.md added)
  baseline_estm/                the standard research-module layout:
    README.md                   the only README of the module (fixed template)
    METHOD.md                   the scientific method (v2, plus v1 kept for comparison)
    CATALOG.md                  one entry per output file
    baseline_estm_paths.py      ALL folder / file locations of the module
    baseline.py  run_sample.py  make_dew_products.py  make_baseline_table.py
    validation/                 scientific checks (scripts only)
    outputs/
      products/                 the module's deliverables, read by later modules
      tables/                   every other table, file name starts with its producer (c1_, r1b_, sample_...)
      figures/                  every figure, same prefix rule
      reports/                  written reports (R1, R1b, sample summary)
    .cache/                     git-ignored: loaded data, DEV runs, logs
_review/                        git-ignored: untracked files that would otherwise be deleted (empty now)
data/                           protected, not in git
```

## Move table (all with `git mv`)

| old path (in `research/baseline_estm/`) | new path | why | references updated |
|---|---|---|---|
| `products/baselines_v2.csv` | `outputs/products/baselines.csv` | deliverable; "_v2" dropped (only one current version) | `baseline.py`, `make_baseline_table.py`, METHOD, README |
| `products/{baseline_uncertainty,dew_risk_hourly,dry_jumps}.csv` | `outputs/products/` (same names) | deliverables | `baseline.py`, both `make_*.py`, METHOD |
| `results/baselines.csv` | `outputs/tables/sample_two_window_baselines.csv` | method v1 = "two-window" (windows before AND after the event); named by content | `run_sample.py`, `c5_flagged.py`, README, CATALOG |
| `results/sensitivity_X.csv` | `outputs/tables/sample_two_window_sensitivity_X.csv` | same | `run_sample.py` |
| `results/summary.md` | `outputs/reports/sample_two_window_summary.md` | a report | `run_sample.py` |
| `figures/overview.png` | `outputs/figures/sample_two_window_overview.png` | figure | `run_sample.py` |
| `figures/examples/*.png` | `outputs/figures/sample_two_window_examples/*.png` | figures | `run_sample.py` |
| `results_v2/baselines.csv` | `outputs/tables/sample_pre_event_baselines.csv` | method v2 = "pre-event" (24 h median before the event only) | `run_sample.py`, R1b report |
| `validation/out/c*.csv(.gz)`, `r1b_*.csv(.gz)` | `outputs/tables/` | tables | `common.py` (OUT), `make_*.py`, `run_sample.py` |
| `validation/out/pass_fail_{per_sublink,summary}.csv` | `outputs/tables/c1_pass_fail_*.csv` | written by C1: prefix by producer | `c1_fake_events.py`, R1 report |
| `validation/out/runtimes.csv` | `outputs/tables/validation_runtimes.csv` | prefix by producer | `common.py` |
| `validation/out/r1b_real_events_v2.csv.gz`, `r1b_c2_delta_v2.csv`, `r1b_v2_splits.csv` | `outputs/tables/r1b_real_events.csv.gz`, `r1b_c2_delta.csv`, `r1b_splits.csv` | R1b is about v2 only: "_v2" is redundant | `r1b_event_level.py`, R1b report |
| `validation/out/figures/*.png` | `outputs/figures/` (`r1b_c2_delta_v2.png` -> `r1b_c2_delta.png`) | figures | `common.py` (FIG), `r1b_event_level.py` |
| `validation/out/REPORT.md` | `outputs/reports/R1_validation.md` | name says what it is (not "v1") | METHOD, README, config.yaml comment |
| `validation/out/REPORT_v2.md` | `outputs/reports/R1b_validation.md` | it is R1b, not "version 2" of R1 | METHOD, README, config.yaml comment |
| `validation/_cache/` (untracked) | `.cache/` (plain `mv`) | one ignored cache place per module | `common.py` (CACHE) |
| (DEV runs) `validation/out_dev/` | `.cache/dev/{tables,figures}` | DEV runs are throw-away | `common.py` |
| `docs/Rain_project.pdf` | `docs/sindy_intro_report_2026-03.pdf` | name describes content | `docs/README.md` |

## Deletion table (tracked files only)

| file | why |
|---|---|
| `research/baseline_estm/.gitignore` | merged into the root `.gitignore` (git history keeps it) |

No outputs or code are deleted: every output has a current producer and a reader or a report that cites it, and every
source file is used (AUDIT sections 2-3).

## `_review/`

Empty: no untracked file needed deleting (the only untracked files are caches and IDE files, which stay ignored).

## Decision log

| # | decision | reason |
|---|---|---|
| D1 | Keep the folder name `baseline_estm` (do not rename to `baseline/`) | nothing breaks; avoids a folder `baseline/` holding a module `baseline.py`; `research/README.md` is updated instead |
| D2 | Outputs split by kind: `outputs/{products,tables,figures,reports}/` | the suggested layout; `products/` added because later modules read them (the module's "API") |
| D3 | Every output file name starts with its producer (`c1_`, `r1b_`, `sample_two_window_`, `validation_`) | "which script made this?" is answered by the name |
| D4 | Keep method v1 active, named `sample_two_window_*`; v2 sample named `sample_pre_event_*` | C5 still reads the v1 sample; names by what differs (principle 6). In code and config the method ids stay `v1` / `v2` (they are documented method names; changing them would be a code change). Open question Q1 |
| D5 | `REPORT.md` / `REPORT_v2.md` -> `R1_validation.md` / `R1b_validation.md` | they are two different validations, not two versions |
| D6 | Drop the redundant `_v2` from `baselines_v2.csv` and the `r1b_*_v2` files | one current version; the method is stated in METHOD and CATALOG |
| D7 | Keep all currently tracked outputs tracked (all < 10 MB); policy: commit tables / figures / reports <= 10 MB, ignore caches and anything bigger | nothing breaks (`make_baseline_table.py` needs `c1_fake_events.csv.gz` etc.), repo stays small |
| D8 | Paths: one file per module, `baseline_estm_paths.py` (not `paths.py`) | the future WAA module will import `baseline.py`; two modules both called `paths` on `sys.path` would shadow each other |
| D9 | Folder locations and the files read across scripts are constants in the paths file; file names that only the writing script uses stay in that script | one source of truth for locations without moving 40 file names away from the code that writes them |
| D10 | Caches in `<module>/.cache/` (ignored by a root rule `.cache/`) | cache lives with the module that makes it; one ignore rule |
| D11 | One root `.gitignore`; change the PySR rule `outputs/` to `/outputs/` and ask PySR runs to write to `.cache/` | the old rule would hide every module's `outputs/` |
| D12 | Leave `.idea/.gitignore` and `.pytest_cache/.gitignore` alone | written by PyCharm / pytest inside ignored folders, not ours |
| D13 | Catalog: one `CATALOG.md` per module + a root `CATALOG.md` index with the regeneration order | the entry sits next to the code that makes the file; the root file answers "how do I rebuild everything" |
| D14 | Pipeline outputs are cataloged by the existing `pipeline/docs/data_dictionary.md` (root CATALOG points to it) | no duplicate doc |
| D15 | `pipeline/` layout unchanged | it already follows the rules (src/, tests/, one config with all paths) |
| D16 | The checker goes to `pipeline/tools/check_organization.py`, not a new top-level `tools/` | the repo-wide tools (`check_boundaries.py`, `update_readme.py`) and the hook already live there; one tools folder |
| D17 | Rewrite `baseline_estm/README.md` to the template; move its v1 method description into METHOD.md ("Method v1"), its column descriptions into CATALOG.md | one doc per topic; README is no longer stale |
| D18 | Fix the doc error "mean of the 10 s TL" -> "5 min `tl_mean`" while moving the v1 text | the code (`load_sublinks`, `res="5min"`) and METHOD.md say 5 min; documentation only, no code change |
| D19 | Do not tick `research/ROADMAP.md` items | progress is a scientific judgment: open question Q3 |
| D20 | `explorations/`: add `README.md` (conventions) only | rule: files there are not moved or renamed |
| D21 | Rename `docs/Rain_project.pdf` | name describes content |
| D22 | Report bugs (Timer appends, silent missing C3 file, C5 old figures) but do not fix | rule: no logic changes |
