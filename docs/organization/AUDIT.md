# Audit of the repository (before the reorganization)

Date: 2026-10-07. Branch point: `main` at `2962c41`. Read-only audit, built from the code, not from file names.
Excluded: `data/` (protected), `research/.vibe-wise/` (protected), `__pycache__`, `.pytest_cache`, `*.egg-info`.

## 1. Git layout

- **One git root**: `D:/UniShit/finalEEProject` (no submodules, no nested repositories).
- **Tracked**: 196 files. Untracked and not ignored: none (the working tree was clean).
- **Ignored, present on disk**: `.idea/`, `pipeline/.idea/` (PyCharm), `.pytest_cache/`, `pipeline/.pytest_cache/`,
  all `__pycache__/`, `pipeline/src/cmlpipe.egg-info/`, `_archive/` (old checks and notes, local only), `data/`,
  `research/.vibe-wise/`, `research/baseline_estm/validation/_cache/` (92 MB: `tl_5min.parquet`, `c1_slots.npz`,
  `run_r1b.log`).

| `.gitignore` | what it does | verdict |
|---|---|---|
| `/.gitignore` | data, Python caches, `.idea/`, Jupyter checkpoints, PySR output, `*.nc` / `*.zip`, Windows junk, `.pytest_cache/`, `check_*.png`, `_archive/`, `*.egg-info/`, `research/.vibe-wise/` | keep (the one root file). **Trap:** the PySR rule `outputs/` matches *every* folder named `outputs`, so a module `outputs/` folder would be silently ignored |
| `research/baseline_estm/.gitignore` | `validation/_cache/`, `validation/out_dev/` | merge into the root file |
| `.idea/.gitignore` | written by PyCharm, inside an ignored folder | not ours, never committed; leave |
| `.pytest_cache/.gitignore`, `pipeline/.pytest_cache/.gitignore` | written by pytest, inside ignored folders | not ours; leave |

## 2. Inventory

### Repository root and `docs/`

| path | type | producer / consumer | status |
|---|---|---|---|
| `README.md` | doc (front door); `AUTO:status` block generated | `pipeline/tools/update_readme.py` (pre-commit hook) | current |
| `CLAUDE.md` | instructions for Claude; imports proposal + OpenMRG summary | hand | current (layout block slightly stale) |
| `LICENSE`, `.gitattributes`, `.gitignore` | config / licence | hand | current |
| `docs/README.md` | doc: what `docs/` is for | hand | current |
| `docs/proposal.md`, `docs/proposal.docx`, `docs/media/fig_cml_path.png` | the proposal (md is imported by CLAUDE.md; docx is the submitted original) | hand | current; md and docx are two formats of one document, not duplicates |
| `docs/Rain_project.pdf` (1.2 MB) | our earlier report "Introduction to SINDy" (March 2026) | hand; read by `docs/README.md` | current; name does not describe the content |

### `pipeline/` (the data layer)

All 37 source files in `src/cmlpipe/` are imported by at least one other file (checked with grep; `__main__.py` is
the `python -m cmlpipe` entry point), so there is **no legacy code**. All pipeline outputs go to `data/work/openmrg/`
and `data/release/openmrg_<version>/`, with paths taken from one place (`configs/openmrg.yaml` +
`cmlpipe.config.REPO_ROOT`): already "one source of truth". The pipeline writes nothing inside the repository except
caches.

| path | type | notes | status |
|---|---|---|---|
| `src/cmlpipe/**` | source (package) | io, qc, geometry, aggregate, events, inventory, release, api, viz, cli | current |
| `tests/*.py` (15) | code tests (pytest); `-m data` marks the real-data tests (`test_real_openmrg.py`) | synthetic data, < 1 min | current |
| `apps/explorer/{app,data,smoke}.py` | Streamlit GUI + its smoke check | `smoke.py` is a check of the app, correctly next to it | current |
| `configs/openmrg.yaml` | config: every pipeline setting and output path | | current |
| `tools/update_readme.py`, `tools/check_boundaries.py`, `tools/githooks/pre-commit` | repo tools | the hook runs `update_readme.py` (which runs pytest) | current |
| `README.md`, `docs/DATA_GUIDE.md`, `docs/methods.md`, `docs/data_dictionary.md`, `docs/howto_access.md` | docs | each one topic, no overlap found | current |
| `docs/roadmap.md` | pipeline progress (source of the README status blocks) | not a duplicate of `research/ROADMAP.md`: different parts | current |
| `docs/openmrg_summary.md` | summary of the OpenMRG paper (imported by CLAUDE.md) | **only one copy exists** | current |
| `docs/reference/cml_libraries_howto.md` | how-to for pycomlink / poligrain | **only one copy exists** | current |
| `docs/reference/example_create_8d_OpenMRG_data.ipynb` | the OpenSense example notebook (third-party reference) | | current (reference) |
| `docs/papers/Andersson2022_OpenMRG.pdf` | the OpenMRG paper | | current (reference) |
| `docs/pipeline_log.pdf` (1.5 MB) | technical report of the data layer (commit `2962c41`) | **no LaTeX source in the repo** | current; source location unknown |
| `examples/quickstart.ipynb` | example notebook for `CMLData` | | current |
| `pyproject.toml`, `environment.yml`, `LICENSE-DATA.md`, `start_explorer.bat` | config / licence / launcher | | current |

### `research/`

| path | type | notes | status |
|---|---|---|---|
| `README.md` | doc; `AUTO:progress` block generated | "Planned layout" (`baseline/`, `notebooks/`, `figures/`) does not match reality (`baseline_estm/`) | current, partly stale |
| `ROADMAP.md` | research progress | **0 / 16 ticked** although the baseline work (v2, R1, R1b, A(t)) exists, so the generated status says "not started in code yet" | stale (science call: see the open questions) |
| `config.yaml` | config: all research settings | the only research config | current |
| `explorations/` | playground: `dry_excursions.py` -> `dry_excursions.csv`; `plot_link_baseline.py` -> `link_10200_{days,event}.png` | each output written next to its script; no README | current (not reorganized, by rule) |

### `research/baseline_estm/`: code

| file | role | reads | writes |
|---|---|---|---|
| `baseline.py` | library: method v1 (`usable_dry`, `estimate_event`, ...) and v2 (`baseline_v2`, `flags_v2`, `dew_*`, `a_of_t`) | `products/dew_risk_hourly.csv`, `products/baselines_v2.csv` | nothing |
| `run_sample.py [v1\|v2]` | sample run (24 sub-links x 18 events); default version from `config.yaml` (v2) | release; v2 also `validation/out/c3_daily_cycle_per_sublink.csv` (silently empty if missing) | v1: `results/{baselines,sensitivity_X}.csv`, `results/summary.md`, `figures/overview.png`, `figures/examples/*.png`; v2: `results_v2/baselines.csv` |
| `validation/common.py` | shared helpers; defines `OUT`, `FIG`, `CACHE` | release | cache `validation/_cache/` |
| `validation/c1_fake_events.py` ... `c6_sensitivity.py`, `r1b_event_level.py` | the scientific checks C1-C6 (R1) and R1b | see dependency map | `validation/out/*.csv(.gz)`, `validation/out/figures/*.png` |
| `validation/run_all.py` | runs C1..C6, R1b in order (~3.5 min) | | `validation/out/runtimes.csv` (appended) |
| `make_dew_products.py` | dew products | release, `validation/out/c4_jumps.csv` | `products/dew_risk_hourly.csv`, `products/dry_jumps.csv` |
| `make_baseline_table.py` | the frozen v2 products | release, `products/dew_risk_hourly.csv`, `validation/out/{c1_fake_events.csv.gz, c3_daily_cycle_per_sublink.csv, r1b_fake_event_offsets.csv.gz}` | `products/{baselines_v2,baseline_uncertainty}.csv` |

All code is used (no legacy code). Output paths are hard-coded in five files (`run_sample.py`, `common.py`,
`c5_flagged.py`, `make_dew_products.py`, `make_baseline_table.py`, `baseline.py`).

### `research/baseline_estm/`: outputs (all tracked, all dated 2026-10-07, all < 10 MB)

| group | files | producer | consumers | status |
|---|---|---|---|---|
| `products/` | `baselines_v2.csv` (8 MB), `baseline_uncertainty.csv`, `dew_risk_hourly.csv`, `dry_jumps.csv` | `make_baseline_table.py`, `make_dew_products.py` | `baseline.py` (`a_of_t`, `dew_*`), `make_baseline_table.py`, METHOD.md; the future WAA module | current: the module's deliverables |
| `results/` | `baselines.csv`, `sensitivity_X.csv`, `summary.md` | `run_sample.py v1` | **`validation/c5_flagged.py` reads `results/baselines.csv`**; README | current (method v1, sample) |
| `figures/` | `overview.png`, `examples/*.png` (11) | `run_sample.py v1` (examples are wiped and redrawn each run) | README | current (method v1, sample) |
| `results_v2/` | `baselines.csv` | `run_sample.py` (v2) | `REPORT_v2.md` | current (method v2, sample) |
| `validation/out/c1*` + `pass_fail_*` | 13 tables + 2 figures | `c1_fake_events.py` | C2, C6, R1b, `make_baseline_table.py`, REPORT.md | current (R1) |
| `validation/out/c2_*`...`c6_*` | 14 tables + 5 figures + 13 `c5_*` figures | `c2_...c6_*.py` | C5 (c4_jumps), R1b (c3), `run_sample.py v2` (c3), `make_dew_products.py` (c4), `make_baseline_table.py` (c3), REPORT.md | current (R1) |
| `validation/out/r1b_*` | 10 tables + 2 figures | `r1b_event_level.py` | `make_baseline_table.py` (offsets), REPORT_v2.md | current (R1b) |
| `validation/out/runtimes.csv` | runtimes | `common.Timer` (appends) | REPORT.md | current |
| `validation/out/REPORT.md`, `REPORT_v2.md` | reports | **hand-written** (no script writes them) | METHOD.md, README, config.yaml comment | current, both |
| `validation/_cache/` | cache (untracked, ignored) | `common.load_panel`, C1 | validation | cache |

Every file found in an output folder is written by the current code at exactly that path. **No orphan outputs**,
with one possible exception: C5 never deletes old `c5_*.png` figures, so a figure from an earlier run would stay.
Re-running checks this (see `REPORT.md`).

## 3. Versions: they are not "old and new copies"

| pair | what it really is | which is current |
|---|---|---|
| `results/` vs `results_v2/` | the **sample run** with **method v1** (6 h windows before/after, mean, flat or a line) vs **method v2** (24 h median before the event, flat, flags). Same sample, same data. Both made by `run_sample.py` (argument `v1`/`v2`). | both: v1 is "kept for comparison" (METHOD.md) and is **still read by C5** |
| `results/baselines.csv` vs `products/baselines_v2.csv` | not a pair: the first is the 24-sub-link **sample** (v1); the second is the **full frozen product** (v2, every train/val event x train sub-link, with uncertainty columns) | both current, different things |
| `REPORT.md` vs `REPORT_v2.md` | not versions of one report: **R1** (slot-level validation of v1 and variants M0-M3) and **R1b** (event-level validation of v2, adoption decision). R1b builds on R1. | both current |
| `c1_*` vs `r1b_*` | not versions: C1 = fake events in dry weather, slot-level errors of variants (R1); R1b = event-level offsets of v2 on the **same** fake events (it reads `c1_fake_events.csv.gz`) | both current |
| `r1b_*_v2.*` (`r1b_real_events_v2.csv.gz`, `r1b_c2_delta_v2.csv/.png`, `r1b_v2_splits.csv`) | "_v2" means "of method v2"; R1b is about v2 only, so the suffix is redundant | current |

History (`git log --follow`): all research files arrived in `a378802` (v1, v2, R1, R1b), then `ebf122b` (R1c freeze),
`afcd1d4` (`use_in_waa_main` excludes excursions). v2 replaced v1 as the method because of R1/R1b (METHOD.md
"Changes v1 -> v2"), not because of a bug fix or a data change.

## 4. Docs

| topic | authoritative doc | others / problems |
|---|---|---|
| project question | `docs/proposal.md` | short restatements in the three READMEs (fine: summaries) |
| OpenMRG dataset | `pipeline/docs/openmrg_summary.md` | no copy elsewhere |
| CML libraries | `pipeline/docs/reference/cml_libraries_howto.md` | no copy elsewhere |
| data / release | `pipeline/docs/DATA_GUIDE.md` (+ `data_dictionary.md`, `howto_access.md`, `methods.md`) | no overlap |
| progress | `pipeline/docs/roadmap.md` (pipeline), `research/ROADMAP.md` (research) | not duplicates; research roadmap stale |
| baseline method | `research/baseline_estm/METHOD.md` (v2) | `baseline_estm/README.md` still describes **v1** as if current, says "Status: sample run", and says the v1 window value is the "mean of the 10 s TL" while the code uses the **5 min** `tl_mean` (METHOD.md is right) |
| baseline validation | `validation/out/REPORT.md` (R1), `REPORT_v2.md` (R1b) | names suggest versions |

Duplicates found: **none** (the suspected copies of `openmrg_summary.md` / `cml_libraries_howto.md` do not exist;
the two roadmaps cover different parts). Stale: `research/README.md` (planned layout), `research/ROADMAP.md`,
`research/baseline_estm/README.md`, and the layout block of `CLAUDE.md`.

## 5. Tests vs validation

Confirmed: `pipeline/tests/` = code-correctness tests (pytest, synthetic data; `-m data` for the real data) and
`research/baseline_estm/validation/` = scientific checks that write tables, figures and feed hand-written reports.
Nothing misplaced. `pipeline/apps/explorer/smoke.py` and `pipeline/tools/check_boundaries.py` are code checks kept
next to what they check (fine). There are no tests for research code (none needed yet).

## 6. Dependency map ("what updates what, when")

```
0  python -m cmlpipe run all --workers 4          raw data -> data/work/openmrg -> data/release/openmrg_v1
1  run_sample.py v1                                release -> results/{baselines,sensitivity_X}.csv, summary.md,
                                                   figures/overview.png, figures/examples/*.png
2  validation/run_all.py  (C1 -> C2 -> C3 -> C4 -> C5 -> C6 -> R1b)
     C1   release -> c1_* , pass_fail_*, cache c1_slots.npz
     C2   c1_fake_events.csv.gz -> c2_*
     C3   release -> c3_*
     C4   release -> c4_*
     C5   results/baselines.csv (step 1!) + c4_jumps.csv -> c5_*
     C6   c1_fake_events.csv.gz -> c6_sensitivity.csv
     R1b  c1_fake_events.csv.gz + c3_daily_cycle_per_sublink.csv -> r1b_*
3  run_sample.py (v2)                              release + c3_daily_cycle_per_sublink.csv -> results_v2/baselines.csv
4  make_dew_products.py                            release + c4_jumps.csv -> products/dew_risk_hourly.csv, dry_jumps.csv
5  make_baseline_table.py                          release + dew_risk_hourly + c1_fake_events + c3 + r1b offsets
                                                   -> products/baselines_v2.csv, baseline_uncertainty.csv
(hand) REPORT.md / REPORT_v2.md / METHOD.md are written by hand from the numbers above.
```

## 7. Things noticed (not changed)

- Root `.gitignore` rule `outputs/` (for PySR) would hide every `outputs/` folder in the repo.
- `common.Timer` **appends** to `runtimes.csv`: running a single check twice adds duplicate rows.
- `run_sample.py v2` silently uses empty daily-cycle amplitudes (all `daily_cycle_large = False`) if C3 has not been
  run yet.
- C5 depends on the v1 sample run (`results/baselines.csv`); if v1 were removed, C5 would break.
- C5 does not remove figures from earlier runs.
- `research/ROADMAP.md` is unticked, so the root README says research is "not started in code yet".
