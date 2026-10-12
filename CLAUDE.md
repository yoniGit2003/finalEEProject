# finalEEProject

Final-year EE project (CellEnMon lab), Yonatan Haber and Yehonatan Tsitseluk: **data-driven discovery of rain-induced
and wet-antenna attenuation laws in commercial microwave links (CML)**. Project documents are in Hebrew (RTL) with
English technical terms. The proposal and the dataset summary are imported here:

@docs/proposal.md
@pipeline/docs/openmrg_summary.md

## Physics and notation (identical in code and text)
- A(t) = TSL(t) - RSL(t) - baseline = A_rain(t) + δ(t); δ(t) = non-rain terms (WAA, baseline residual).
- A_rain = ∫₀ᴸ a·R(x,t)^b dx ≈ a·R̄^b·L (exact only for uniform rain); extended model A_rain = a·R̄^b·f(L, R̄, Z), uniform rain gives f = L.
- WAA is local, dynamic (wetting / drying), does not scale with L; a constant WAA looks like f = L + L0, so isolate WAA (drying tails) before learning f. Never select links by how well the power law fits (circular).
- Methods: weak-form SINDy for WAA (rain as input), PySR for f; real data only (no simulator). Baseline is research and comes first. Gray-box: readable equations; the RNN [11] is only a baseline. Radar along the path is the preferred rain reference.

## Layout
```
README.md, CLAUDE.md, LICENSE (MIT), .gitignore (the only one, written by hand), .gitattributes,
environment.yml (+ environment.lock.yml = exact versions), .github/workflows/check.yml (fast tests on every push)
.claude/commands/updatedocs.md   the /updatedocs command (research log + reproduce page)
docs/        REPRODUCE.md (how to rerun everything), ROADMAP.md (the ONLY roadmap: checklist + WAA plan),
             proposal, research_log.tex / .pdf (build: python docs/build_research_log.py), pipeline_log.pdf, papers/
pipeline/    the data layer: src/cmlpipe (io, qc, geometry, aggregate, events, inventory, release, api, viz, cli),
             apps/explorer (Streamlit GUI), configs/openmrg.yaml, tests/, tools/check_boundaries.py, docs/,
             examples/, pyproject.toml (single dependency list), get_data.py, get_caches.py, LICENSE-DATA.md
research/    config.yaml (all settings), steps.yaml (every step in run order), __main__.py (python -m research)
  common/        save.py (save), report.py (md_table, rel), plotting.py (style_axes), config.py (load_config)
  baseline/      README, METHOD, paths.py, baseline.py, make_*.py, run_sample.py, validation/, outputs/
  wet_antenna/   README, METHOD, paths.py, libraries (simulator, per_tail_fit, pooled_fit, sindy, free_curve, tails,
                 fit_jobs), one entry script per step, outputs/
  explorations/  playground, local only (git-ignored except its README)
  .cache/        saved fits and logs, git-ignored (research/.cache/<module>/<script>/, logs/)
data/        NOT in git: raw/openmrg (never modified), work/openmrg (pipeline outputs), release/openmrg_<version>
```
- Outputs: every research script writes to `research/<module>/outputs/{products,tables,figures,reports}`; file names
  start with the script name; one report per script. **Everything there is committed**, except the big files named
  in `.gitignore`. There are no importance levels, manifests or catalogs.
- A module's docs are exactly `README.md` (what, status, how to run, files) and `METHOD.md` (method, decisions of
  each step with numbers, open issues). The roadmap is `docs/ROADMAP.md`: patch it in place, never rewrite it.
- Config paths (`data/...`) resolve against the repository root (`cmlpipe.config.REPO_ROOT`); code paths against
  `PIPELINE_ROOT` (= `pipeline/`). Research file locations come from the module's `paths.py`, never from strings in
  scripts.

## How to run (Windows; use `C:\Users\yoni\miniconda3\python.exe`, not the Anaconda `python`)
- Everything is run from the repository root. Full sequence: `docs/REPRODUCE.md`.
- Install: `conda env create -f environment.yml` or `pip install -e "pipeline[dev,explorer,research]"`.
- Data: `python pipeline/get_data.py` (Zenodo -> data/raw, then the pipeline build); fit caches: `python pipeline/get_caches.py`.
- Pipeline: `python -m cmlpipe run all --workers 4` (up-to-date stages skip), `python -m cmlpipe stages`, `python -m cmlpipe guide`.
- Research: `python -m research` (menu), `python -m research run <step>|all`, `python -m research compare`; a script directly: `python -m research.wet_antenna.extract_tails`.
- Data access: `from cmlpipe.api import CMLData; d = CMLData.open("data/release/openmrg_v1")` (see pipeline/docs/howto_access.md).
- Explorer: `python -m cmlpipe explore` (localhost:8501) or `pipeline/start_explorer.bat`; check: `python pipeline/apps/explorer/smoke.py`.
- Tests (run from `pipeline/`): `python -m pytest` (synthetic, < 1 min); `-m data` for the real-data tests. Boundary: `python pipeline/tools/check_boundaries.py`.
- Raw `cml.nc` is ~4.6 GB: never load it whole (slice by time / sub-link; missing value 1e10 -> NaN). Never open many NetCDF files with xarray at once. No dask.

## Rules
- Do not change the science by accident: a refactor keeps parameters, seeds, methods and numbers. After changing a
  research script, rerun it and compare its tables with the committed ones (`python -m research compare`).
- The fit caches in `research/.cache/wet_antenna/` hold about 9 hours of fitting: never delete them without asking.
- A new research script: `main()` + `if __name__ == "__main__"`, locations in `paths.py`, settings in `config.yaml`,
  outputs with `research.common.save.save`, an entry in `research/steps.yaml`, then `python -m research docs`
  (a test fails if the entry or the reproduce page is missing).
- Protected: `data/` (never edit by hand; `data/raw/` never at all), `research/.vibe-wise/` (plugin notes: never
  commit, move or delete), the files in `research/explorations/`.
- Frozen: `events.csv`, `split_events.csv`, `split_links.csv`, `event_links.csv`, `dry_periods.csv` (data/work/openmrg/tables and the release) never change unless explicitly asked (`--force-split` / `--force`). Rebuild the release only with `--force` or as a new `--version`.
- `pipeline/docs/DATA_GUIDE.md` and `pipeline/docs/data_dictionary.md` are written by the release stage: do not edit them by hand.
- One-way boundary: research may import only `cmlpipe.api` / `cmlpipe.viz`; nothing in `pipeline/` imports `research`.
- Write to the students' level: explain what a step does and why; simple, commented, rerunnable scripts over one-off notebooks.
- Figures must not rely on colour alone (line styles, hatching, text labels).
- Keep the proposal's notation (A_rain, R̄, L, f(L, R̄, Z), δ(t)). Verify literature claims against the proposal's references; never invent citations.
- Answer in English unless asked for Hebrew; in Hebrew keep English terms in Latin script and check RTL ordering.
