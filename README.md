# Discovering the Laws of Microwave-Link Attenuation

Final-year EE project, CellEnMon lab (opportunistic environmental sensing).
Team: Yonatan Haber and Yehonatan Tsitseluk.

**To rerun everything we did, in order: [`docs/REPRODUCE.md`](docs/REPRODUCE.md).**
What the research has found so far: [`docs/research_log.pdf`](docs/research_log.pdf).

## The question

When it rains, the microwave links that connect cellular towers lose signal, so the phone network can be used as a
rain sensor. The usual model says the rain loss is A_rain = a·R̄^b·L: a power of the average rain rate R̄ along the
link, times the link length L. That is exact only when the rain is uniform along the path, and real links of
different lengths do not always follow it. A second effect, wet-antenna attenuation (WAA), comes from water on the
antenna covers: it does not grow with L and it stays for a while after the rain stops. We ask how the loss really
depends on the link: A_rain = a·R̄^b·f(L, R̄, Z), where uniform rain gives f = L and Z stands for path and link
features. We keep the known physics and learn only the missing parts as readable equations (a gray-box approach):
a dynamic WAA model first, then the correction f, and finally rain estimates from the combined model.

## Setup

Windows, Linux or macOS with conda, from the repository root:

```
conda env create -f environment.yml
conda activate cmlpipe
```

This installs Python 3.12, the `cmlpipe` package (editable, from `pipeline/`) and what the research scripts need.

## Get the data

Nothing big is in git. Two scripts fetch or rebuild it:

```
python pipeline/get_data.py      # raw OpenMRG from Zenodo (318 MB) -> data/raw/, then builds the release in data/release/
python pipeline/get_caches.py    # the saved fits of the long wet-antenna steps -> research/.cache/ (saves about 9 hours)
```

## Run

| command | what it does |
|---|---|
| `python -m research` | the menu: every research step in order, what it does, how long it takes |
| `python -m research run all` | run every research step in order (about 20 minutes with the caches) |
| `python -m research run <step>` | run one step |
| `python -m research compare` | compare your tables with the committed ones |
| `python -m cmlpipe run all --workers 4` | build or update the data release (finished stages are skipped) |
| `python -m cmlpipe explore` | open the data explorer in the browser |
| `python -m pytest` (in `pipeline/`) | the fast tests (about a minute, no data needed) |
| `python docs/build_research_log.py` | build the research log PDF (needs LaTeX) |

Details, run times and the full list of steps: [`docs/REPRODUCE.md`](docs/REPRODUCE.md).

## Layout

```
docs/        proposal, research log (PDF + LaTeX), pipeline report, roadmap, REPRODUCE.md, papers
pipeline/    the data layer: package cmlpipe (raw data -> versioned data release), explorer, tests, get_data.py
research/    the research code: config.yaml (all settings), steps.yaml (all steps), one folder per module
  baseline/      the dry-weather baseline of every sub-link, and A(t)
  wet_antenna/   drying tails and the wet-antenna (WAA) law
  common/        small shared helpers (saving, report tables, plot style)
  explorations/  playground (local only, not in git)
data/        raw data, pipeline work files and the release (not in git; made by pipeline/get_data.py)
```

Each top-level folder has its own README. Two rules keep the parts apart: research code reads data only through
`cmlpipe.api.CMLData` (the release pinned in `research/config.yaml`), and the pipeline never imports research
(`python pipeline/tools/check_boundaries.py` checks both).

What is committed: code, settings, documents, and everything the research scripts write to
`research/<module>/outputs/` (reports, tables, figures), except three big files that the scripts rebuild in minutes.
Data and caches are never committed.

## Licence and citation

Code: MIT ([`LICENSE`](LICENSE)). Data release: CC BY-SA 4.0 ([`pipeline/LICENSE-DATA.md`](pipeline/LICENSE-DATA.md)).
It is derived from OpenMRG: Andersson, J. C. M., Olsson, J., van de Beek, R. (C. Z.), and Hansryd, J. (2022),
*Earth Syst. Sci. Data*, 14, 5411-5426, https://doi.org/10.5194/essd-14-5411-2022 (data:
https://doi.org/10.5281/zenodo.7107689).
