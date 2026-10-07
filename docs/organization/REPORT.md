# Reorganization report (branch `reorg/project-structure`, 2026-10-07)

Audit: [`AUDIT.md`](AUDIT.md). Plan, move table and decision log (D1-D22): [`PLAN.md`](PLAN.md). Rules: `ORGANIZATION.md`.

## What changed
- `research/baseline_estm/` outputs -> `outputs/{products,tables,figures,reports}/`; cache -> `.cache/`; all locations
  in `baseline_estm_paths.py`. Renamed: `baselines_v2.csv` -> `baselines.csv`; `results/` -> `sample_two_window_*`
  (method v1); `results_v2/` -> `sample_pre_event_*` (method v2); `REPORT.md` / `REPORT_v2.md` -> `R1_validation.md` /
  `R1b_validation.md`; `pass_fail_*` -> `c1_pass_fail_*`; `runtimes.csv` -> `validation_runtimes.csv`; `r1b_*_v2` -> `r1b_*`.
- `docs/Rain_project.pdf` -> `docs/sindy_intro_report_2026-03.pdf`.
- One root `.gitignore` (PySR rule `outputs/` -> `/outputs/`; added `.cache/`, `_cache/`, `_review/`).
- New: `ORGANIZATION.md` (imported by `CLAUDE.md`), `CATALOG.md` + `research/baseline_estm/CATALOG.md`,
  `.claude/commands/tidy.md`, `pipeline/tools/check_organization.py` + test, `research/explorations/README.md`.
- READMEs updated (baseline_estm README rewritten from the template; v1 method moved into METHOD.md).
- Deleted: only `research/baseline_estm/.gitignore` (merged). `_review/`: empty. `pipeline/`: no moves (already fine).

## Verification
- pytest: 130 passed (synthetic) + 15 passed (`-m data`, real data). `check_boundaries`, `check_organization`: OK.
- Re-ran `run_sample.py v1`, `run_sample.py`, `validation/run_all.py`, `make_dew_products.py`, `make_baseline_table.py`
  from the new locations: every file byte-identical (gz: identical after decompression), except below.

## Open questions (your scientific judgment)
1. **Are method v1 and check C5 still needed?** C5 reads the v1 sample run. Recommended: keep while writing up R1; later move v1 outputs to `archive/` and drop C5 from `run_all.py`.
2. **`r1b_real_events.csv.gz` is stale**: written before R1c, it lacks `dew_risk_during_event` / `dew_share_during_event`. A re-run gives identical values plus these two columns. Recommended: re-run R1b and commit.
3. **`research/ROADMAP.md` has 0/16 ticked** although the baseline work is done, so the README says "not started". Recommended: tick items 2-4 of section 1 (item 1 mentions the min-attenuation method [8], not done).

## Noticed, not changed
- `run_sample.py` (v2) silently sets all `daily_cycle_large = False` if C3 has not been run.
- C5 never deletes its old figures (no stale ones exist now).
- A DEV run of C1 overwrites `.cache/c1_slots.npz` (cache only, harmless).
- `pipeline/docs/pipeline_log.pdf` has no LaTeX source in the repository.

## Using /tidy
`/tidy` (whole repo) or `/tidy research/baseline_estm`. It runs the checker, fixes on a new branch, verifies and reports.
To run the checker on every commit, add `"$PY" pipeline/tools/check_organization.py || echo "pre-commit: organization problems" >&2`
to `pipeline/tools/githooks/pre-commit` (not installed). The test `test_check_organization.py` already runs in the hook's pytest.
