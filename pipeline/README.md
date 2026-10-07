# pipeline/: the data layer (cmlpipe)

The data layer turns the raw OpenMRG data (Gothenburg, summer 2015: 364 links at 10 s, radar, rain gauges, weather)
into a versioned, read-only **data release**, with an access API (`CMLData`) and an explorer GUI.
It cleans, flags, aligns and documents the data, and it defines rain events, dry periods and a frozen train / test
split. It does **not** estimate a baseline, compute the attenuation A(t) or model anything: that is research
(see [`../research/`](../research/)).

Start with the plain-language data guide: [`docs/DATA_GUIDE.md`](docs/DATA_GUIDE.md).

## Setup

From the repository root (Windows, conda; no dask):

```powershell
conda env create -f pipeline/environment.yml; conda activate cmlpipe   # installs cmlpipe too (editable)
pip install -e "pipeline[dev,explorer]"                                # same, in an existing environment
git config core.hooksPath pipeline/tools/githooks                      # once per clone: README hook
```

[`pyproject.toml`](pyproject.toml) is the single list of dependencies.

## Get the data

The release must end up in `data/release/openmrg_v1/` at the repository root (`data/` is not in git).

- **Option A, rebuild** from the public raw data (OpenMRG v1.1, https://doi.org/10.5281/zenodo.7107689). Unpack it
  so that `data/raw/openmrg/cml/cml.nc`, `cml/cml_metadata.csv`, `gauges/city/CityGauges-2015JJA.csv`,
  `gauges/smhi/GbgA-71420-2015JJA.csv` (each with its `-metadata.csv`) and `radar/radar.nc` exist, then run
  `python -m cmlpipe run all --workers 4`. `data/raw/` is never modified.
- **Option B, download** `openmrg_v1_release.zip` (1.6 GB, CC BY-SA 4.0) from
  https://github.com/yoniGit2003/finalEEProject/releases/tag/data-openmrg-v1 (private repository: a GitHub login
  with access is needed). Check it with `(Get-FileHash openmrg_v1_release.zip -Algorithm SHA256).Hash`, which must be
  `17E9A9048AD0DAE51166C59BA0D3FA760C2E090B4757D70A0B2CDD033AC009F2`, then unzip it into `data/release/`.

## How it is built

`python -m cmlpipe run all` runs the stages in order (`python -m cmlpipe stages` lists them): ingest, QC flags,
geometry, aggregation (1 / 5 / 15 / 60 min), events, dry periods and split, inventory and quality status, per-link
reference state, figures, release. Up-to-date stages are skipped: each output records a stage hash of only the
config keys it uses. Work files go to `data/work/openmrg/`. Every definition is in [`docs/methods.md`](docs/methods.md);
every setting in [`configs/openmrg.yaml`](configs/openmrg.yaml). The frozen tables (events, split, event_links,
dry periods) change only with `--force-split` / `--force`, and the release is never overwritten without `--force`.

## Read the data

```python
from cmlpipe.api import CMLData
d = CMLData.open("data/release/openmrg_v1")                     # path from the repository root
links = d.links(length_km=(1, 3), status="good", split="train")
df = d.event_window("E20150727_2140", cml_id=10316, pad="3h")   # CML + radar + gauges + weather + states
pairs = d.event_pairs(split=("train", "train"))                  # event x link pairs for modelling
```

More: [`docs/howto_access.md`](docs/howto_access.md), [`examples/quickstart.ipynb`](examples/quickstart.ipynb), every
column and variable in [`docs/data_dictionary.md`](docs/data_dictionary.md).

**Explorer.** `python -m cmlpipe explore` (or double-click [`start_explorer.bat`](start_explorer.bat)) opens
http://localhost:8501 on this computer only. `python pipeline/apps/explorer/smoke.py` checks it without a browser.

**Tests.** From `pipeline/`: `python -m pytest` (synthetic data, under a minute); `python -m pytest -m data` runs the
real-data tests.

**Data on this computer** *(auto)*:

<!-- AUTO:data -->
L1 CML layout (from `src/cmlpipe/io/schema.py`): dims `('time', 'sublink_id', 'cml_id')`; data variables `tsl` [dBm], `rsl` [dBm], `tl` [dB].

- L1 CML files: **92 days** (2015-06-01 to 2015-08-31), 0.49 GB, in `data/work/openmrg/L1_standard/cml/`
- Last pipeline run: 2026-10-04T15:08:19Z, steps ingest_cml, ingest_gauges, ingest_met, ingest_radar, qc_stats, qc_flags, qc_status, geometry, aggregate_gauges, aggregate_radar, aggregate_cml, events, split, event_links, significance, inventory, event_links_extra, link_state, inventory_figures, quicklooks, link_cards, run_report, release, 276 days, status **ok**

*(Processed data is not in git. These lines describe the computer that made the last commit.)*
<!-- /AUTO:data -->

## Map of pipeline/

```
src/cmlpipe/    io, qc, geometry, aggregate, events, inventory, release, api (CMLData), viz, cli
apps/explorer/  the explorer GUI (app.py, data.py, smoke.py)
configs/        openmrg.yaml (all settings)
docs/           DATA_GUIDE, methods, data_dictionary, howto_access, roadmap, OpenMRG notes and paper
tests/          pytest suite        tools/   update_readme.py, check_boundaries.py, githooks/
examples/       quickstart notebook
```

## Status and next steps *(auto, from [`docs/roadmap.md`](docs/roadmap.md))*

<!-- AUTO:progress -->
Status: **done (release openmrg_v1 built and published; the open items are extensions)**; 21 / 26 roadmap items done (tick items in [`docs/roadmap.md`](docs/roadmap.md)).

Next steps:
- [ ] Wet/dry classification from the CML signal itself
- [ ] Noise and outlier handling beyond QC flags (level shifts, pair mismatch)
- [ ] Link selection for Stage 1 and Stage 2
- [ ] OpenRainER adapter (sub-links named `channel1` / `channel2`)
- [ ] Check whether OpenMesh has a rain reference before using it
<!-- /AUTO:progress -->

Licences: data CC BY-SA 4.0 ([`LICENSE-DATA.md`](LICENSE-DATA.md)); code MIT ([`../LICENSE`](../LICENSE)).
