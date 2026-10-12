# pipeline/: the data layer (cmlpipe)

The data layer turns the raw OpenMRG data (Gothenburg, summer 2015: 364 links at 10 s, radar, rain gauges, weather)
into a versioned, read-only **data release**, with an access API (`CMLData`) and an explorer GUI.
It cleans, flags, aligns and documents the data, and it defines rain events, dry periods and a frozen train / test
split. It does **not** estimate a baseline, compute the attenuation A(t) or model anything: that is research
(see [`../research/`](../research/)).

Start with the plain-language data guide: [`docs/DATA_GUIDE.md`](docs/DATA_GUIDE.md).

## Setup

From the repository root (no dask):

```
conda env create -f environment.yml          # then: conda activate cmlpipe   (installs cmlpipe too, editable)
pip install -e "pipeline[dev,explorer]"      # the same, in an existing environment (add ,research for research/)
```

[`pyproject.toml`](pyproject.toml) is the single list of dependencies.

## Get the data

`data/` is not in git. From the repository root:

```
python pipeline/get_data.py
```

It downloads the public raw data (OpenMRG v1.1, https://doi.org/10.5281/zenodo.7107689, 318 MB) into
`data/raw/openmrg/`, checks it, and runs `python -m cmlpipe run all --workers 4`, which builds the release in
`data/release/openmrg_v1/`. If the raw files are already there, only the build runs. `data/raw/` is never modified.

By hand instead: unpack the Zenodo zip so that `data/raw/openmrg/cml/cml.nc`, `cml/cml_metadata.csv`,
`gauges/city/CityGauges-2015JJA.csv`, `gauges/smhi/GbgA-71420-2015JJA.csv` (each with its `-metadata.csv`) and
`radar/radar.nc` exist, then run `python -m cmlpipe run all --workers 4`.

A built copy of the release (`openmrg_v1_release.zip`, 1.6 GB, CC BY-SA 4.0) is also attached to the GitHub release
`data-openmrg-v1` of this repository; its SHA256 is
`17E9A9048AD0DAE51166C59BA0D3FA760C2E090B4757D70A0B2CDD033AC009F2`. Unzip it into `data/release/` to skip the build.

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
real-data tests. The fast tests also run on GitHub on every push (`.github/workflows/check.yml`).

## Map of pipeline/

```
src/cmlpipe/    io, qc, geometry, aggregate, events, inventory, release, api (CMLData), viz, cli
apps/explorer/  the explorer GUI (app.py, data.py, smoke.py)
configs/        openmrg.yaml (all settings)
docs/           DATA_GUIDE.md and data_dictionary.md (both written by the release stage: do not edit by hand),
                methods.md, howto_access.md, openmrg_summary.md (the source dataset), reference/
tests/          pytest suite
tools/          check_boundaries.py (research uses only cmlpipe.api / cmlpipe.viz; the pipeline never imports research)
examples/       quickstart notebook
get_data.py     download the raw data and build the release
get_caches.py   download the saved fits of the research scripts (or pack new ones)
```

## Status

Done: release `openmrg_v1` is built. The open items are extensions (wet/dry classification from the CML signal,
OpenRainER adapter, ...): see [`../docs/ROADMAP.md`](../docs/ROADMAP.md), part A.

Licences: data CC BY-SA 4.0 ([`LICENSE-DATA.md`](LICENSE-DATA.md)); code MIT ([`../LICENSE`](../LICENSE)).
