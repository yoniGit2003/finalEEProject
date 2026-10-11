# finalEEProject

Final-year EE project (CellEnMon lab), Yonatan Haber and Yehonatan Tsitseluk: **data-driven discovery of rain-induced
and wet-antenna attenuation laws in commercial microwave links (CML)**. Project documents are in Hebrew (RTL) with
English technical terms. The proposal and the dataset summary are imported here:

@docs/proposal.md
@pipeline/docs/openmrg_summary.md
@ORGANIZATION.md

## Physics and notation (identical in code and text)
- A(t) = TSL(t) - RSL(t) - baseline = A_rain(t) + δ(t); δ(t) = non-rain terms (WAA, baseline residual).
- A_rain = ∫₀ᴸ a·R(x,t)^b dx ≈ a·R̄^b·L (exact only for uniform rain); extended model A_rain = a·R̄^b·f(L, R̄, Z), uniform rain gives f = L.
- WAA is local, dynamic (wetting / drying), does not scale with L; a constant WAA looks like f = L + L0, so isolate WAA (drying tails) before learning f. Never select links by how well the power law fits (circular).
- Methods: weak-form SINDy for WAA (rain as input), PySR for f; real data only (no simulator). Baseline is research and comes first. Gray-box: readable equations; the RNN [11] is only a baseline. Radar along the path is the preferred rain reference.

## Layout
Where things go and how files are named: `ORGANIZATION.md` (imported above; `/tidy` applies it). What every output
file means: `CATALOG.md` (root index) -> each module's `CATALOG.md`.
```
README.md, ORGANIZATION.md, CATALOG.md, LICENSE (MIT), CLAUDE.md, .gitignore (the only one), .gitattributes
.claude/commands/tidy.md   the /tidy command
docs/        project-level material only (proposal, reports / figures); the two project reports: pipeline_log.pdf
             (data layer) and research_log.pdf + .tex (research at the level of ideas; /tidy updates it, build:
             python pipeline/tools/build_research_log.py); docs/organization/ = records of the reorg
research/    research code; uses ONLY cmlpipe.api (+ cmlpipe.viz) and the release pinned in research/config.yaml;
             research/ROADMAP.md = research tasks (pipeline tasks: pipeline/docs/roadmap.md);
             baseline_estm/ = the first module (README, METHOD, CATALOG, baseline_estm_paths.py, code, validation/,
             outputs/{products,tables,figures,reports}/ = important (key / product, committed), legacy/ = diag / scratch
             moved there by /tidy (git-ignored, local only), git-ignored .cache/); explorations/ = playground
pipeline/    the data layer: src/cmlpipe (io, qc, geometry, aggregate, events, inventory, release, api, viz, cli),
             apps/explorer (Streamlit GUI), configs/openmrg.yaml, tests/, tools/ (update_readme, check_boundaries,
             check_organization, build_research_log, githooks), docs/ (DATA_GUIDE, methods, data_dictionary,
             howto_access, roadmap, openmrg_summary, papers, reference), examples/, pyproject.toml (single dependency
             list), environment.yml, LICENSE-DATA.md (CC BY-SA 4.0)
data/        at the repository root, NOT in git: raw/openmrg (never modified), work/openmrg (pipeline outputs),
             release/openmrg_<version> (read-only release)
_archive/    local only, gitignored (old checks and notes); history before step 8: branch dev-history
_review/     local only, gitignored: untracked files that look like junk, moved there by a tidy (the user decides)
```
Config paths (`data/...`) resolve against the repository root (`cmlpipe.config.REPO_ROOT`); code paths against
`PIPELINE_ROOT` (= `pipeline/`). Dataset paths are never part of a stage hash. Research output locations come from
each module's `<module>_paths.py`, never from strings scattered in scripts.

## How to run (Windows; use `C:\Users\yoni\miniconda3\python.exe`, not the Anaconda `python`)
- Install: `conda env create -f pipeline/environment.yml` or `pip install -e "pipeline[dev,explorer]"`; hook: `git config core.hooksPath pipeline/tools/githooks`.
- READMEs: root = front door, `pipeline/README.md` = data layer, `research/README.md` = research; their status blocks are generated from the two roadmaps by `python pipeline/tools/update_readme.py` (also the pre-commit hook). Edit progress only in the roadmaps.
- Pipeline: `python -m cmlpipe run all --workers 4` (up-to-date stages skip), `python -m cmlpipe stages`, `python -m cmlpipe guide`.
- Data access: `from cmlpipe.api import CMLData; d = CMLData.open("data/release/openmrg_v1")` (see pipeline/docs/howto_access.md).
- Explorer: `python -m cmlpipe explore` (localhost:8501) or `pipeline/start_explorer.bat`; check: `python pipeline/apps/explorer/smoke.py`.
- Tests (run from `pipeline/`): `python -m pytest` (synthetic, < 1 min); `-m data` for the real-data tests. Boundary: `python pipeline/tools/check_boundaries.py`. Organization: `python pipeline/tools/check_organization.py`.
- Raw `cml.nc` is ~4.6 GB: never load it whole (slice by time / sub-link; missing value 1e10 -> NaN). Never open many NetCDF files with xarray at once. No dask.

## Rules
- A tidy never deletes a file: it moves (archive/, a module's local legacy/, _review/). See ORGANIZATION.md.
- Protected (see ORGANIZATION.md): `data/`, `research/.vibe-wise/` (plugin notes: never commit, move or delete), the files in `research/explorations/`.
- Frozen: `events.csv`, `split_events.csv`, `split_links.csv`, `event_links.csv`, `dry_periods.csv` (data/work/openmrg/tables and the release) never change unless explicitly asked (`--force-split` / `--force`). Do not modify `data/raw/`; rebuild the release only with `--force` or as a new `--version`.
- One-way boundary: research may import only `cmlpipe.api` / `cmlpipe.viz`; nothing in `pipeline/` imports `research`.
- Write to the students' level: explain what a step does and why; simple, commented, rerunnable scripts over one-off notebooks.
- Keep the proposal's notation (A_rain, R̄, L, f(L, R̄, Z), δ(t)). Verify literature claims against the proposal's references; never invent citations.
- Answer in English unless asked for Hebrew; in Hebrew keep English terms in Latin script and check RTL ordering.
