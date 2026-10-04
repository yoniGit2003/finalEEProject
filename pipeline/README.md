# pipeline/: the data layer (cmlpipe)

How the OpenMRG data was prepared for the research: a Python package (`cmlpipe`) that turns the raw OpenMRG files
into a versioned, read-only **release** (`data/release/openmrg_v1/`), an access API (`CMLData`), and an explorer
GUI. Plain-language description of the data: [`docs/DATA_GUIDE.md`](docs/DATA_GUIDE.md).

> Blocks marked *(auto)* are rewritten on every commit by `tools/update_readme.py` (hook `tools/githooks/pre-commit`).

## Install

From the repository root (Windows, conda; no dask):

```powershell
conda env create -f pipeline/environment.yml; conda activate cmlpipe   # also installs cmlpipe (editable)
pip install -e "pipeline[dev,explorer]"                                # same, in an existing environment
git config core.hooksPath pipeline/tools/githooks                      # once per clone: README hook
```

`pyproject.toml` is the single list of dependencies (`[dev]`: tests and notebook, `[explorer]`: Streamlit).

## Build the release

Raw inputs (OpenMRG v1.1, https://doi.org/10.5281/zenodo.7107689), unpacked at the **repository root** under `data/`
(which is not in git; `data/raw/` is never modified):

```
data/raw/openmrg/cml/cml.nc, cml/cml_metadata.csv          CML (about 4.6 GB) and link metadata
data/raw/openmrg/gauges/city/CityGauges-2015JJA.csv (+ -metadata.csv)
data/raw/openmrg/gauges/smhi/GbgA-71420-2015JJA.csv (+ -metadata.csv)
data/raw/openmrg/radar/radar.nc
```

```powershell
python -m cmlpipe run all --workers 4      # every stage in order; up-to-date stages are skipped
python -m cmlpipe stages                   # list the stages; one stage: python -m cmlpipe run events
python -m cmlpipe guide                    # rewrite DATA_GUIDE.md / LICENSE.txt / README.md of the release
```

Stages: ingest (CML 10 s, gauges, weather, radar) → QC flags (L2, flag-don't-fix) → geometry (radar cells along each
path, gauge matches) → aggregation (L3: 1 / 5 / 15 / 60 min) → events, dry periods, frozen split → significance,
inventory, quality status → per-link reference state → figures → release (Parquet + Zarr, catalog, data dictionary,
guide, licence). Work files go to `data/work/openmrg/`. Every output records its **stage hash** (only the config keys
that stage uses, never a dataset path) and the git commit, so a stage reruns only when its settings or inputs change.
The frozen files (`events.csv`, `split_events.csv`, `split_links.csv`, `event_links.csv`, `dry_periods.csv`) never
change unless `--force-split` / `--force` is given on purpose; the release is never overwritten without `--force`.
All settings: [`configs/openmrg.yaml`](configs/openmrg.yaml). Definitions: [`docs/methods.md`](docs/methods.md).

## Use the data

```python
from cmlpipe.api import CMLData
d = CMLData.open("data/release/openmrg_v1")
df = d.event_window("E20150727_2140", cml_id=10316, pad="3h")   # CML + radar + gauges + weather + states
pairs = d.event_pairs(split=("train", "train"))                  # event x link pairs for modelling
```

Examples: [`docs/howto_access.md`](docs/howto_access.md), [`examples/quickstart.ipynb`](examples/quickstart.ipynb).
Every column and variable: [`docs/data_dictionary.md`](docs/data_dictionary.md).

**Explorer** (Streamlit + Plotly, reads only through `CMLData`): `python -m cmlpipe explore`, or double-click
`start_explorer.bat`; opens http://localhost:8501 (this computer only). `python apps/explorer/smoke.py` checks it
without a browser.

**Data status on this computer** *(auto)*:

<!-- AUTO:data -->
L1 CML layout (from `src/cmlpipe/io/schema.py`): dims `('time', 'sublink_id', 'cml_id')`; data variables `tsl` [dBm], `rsl` [dBm], `tl` [dB].

- L1 CML files: **92 days** (2015-06-01 to 2015-08-31), 0.49 GB, in `data/work/openmrg/L1_standard/cml/`
- Last pipeline run: 2026-10-04T15:08:19Z, steps ingest_cml, ingest_gauges, ingest_met, ingest_radar, qc_stats, qc_flags, qc_status, geometry, aggregate_gauges, aggregate_radar, aggregate_cml, events, split, event_links, significance, inventory, event_links_extra, link_state, inventory_figures, quicklooks, link_cards, run_report, release, 276 days, status **ok**

*(Processed data is not in git. These lines describe the computer that made the last commit.)*
<!-- /AUTO:data -->

## Layout

```
src/cmlpipe/     io (ingest), qc, geometry, aggregate, events, inventory, release, api (CMLData), viz, cli
apps/explorer/   the explorer GUI (app.py), its data access (data.py) and smoke check (smoke.py)
configs/         openmrg.yaml          docs/   DATA_GUIDE, methods, data dictionary, how-to, roadmap, OpenMRG notes
tests/           python -m pytest (synthetic data, < 1 min); python -m pytest -m data (real data and products)
tools/           update_readme.py, check_boundaries.py, githooks/      examples/   quickstart notebook
```

**Boundary.** `research/` (at the repository root) may use only `cmlpipe.api` and `cmlpipe.viz`; nothing in
`pipeline/` imports research. Checked by `python pipeline/tools/check_boundaries.py` and by the tests.

## Progress *(auto)*

<!-- AUTO:progress -->
Overall: **21 / 45** items done. Edit `pipeline/docs/roadmap.md` to tick items.

| Part | Done |
|---|---|
| 1. Pipeline (data) | 21 / 27 |
| 2. Simulator (known ground truth) | 0 / 4 |
| 3. Stage 1: WAA | 0 / 5 |
| 4. Stage 2: path law f(L, R̄, Z) | 0 / 3 |
| 5. Rain estimation and evaluation | 0 / 6 |

**Next up** (first open items):

- [ ] Baseline estimation [8] *(1. Pipeline (data))*
- [ ] Wet/dry classification from the CML signal itself *(1. Pipeline (data))*
- [ ] Noise and outlier handling beyond QC flags (level shifts, pair mismatch) *(1. Pipeline (data))*
- [ ] Link selection for Stage 1 and Stage 2 *(1. Pipeline (data))*
- [ ] OpenRainER adapter (sub-links named `channel1` / `channel2`) *(1. Pipeline (data))*

<details><summary>Full checklist</summary>

**1. Pipeline (data)**

- [x] Inspect the raw OpenMRG CML file (`_archive/docs/data_notes/openmrg_cml_inspection.md`, local only)
- [x] Config file + config hash, so every output records the settings that made it
- [x] L1 CML schema (OpenSense layout, GHz / km) and validator
- [x] OpenMRG adapter: raw `cml.nc` + metadata -> daily L1 files (missing value -> NaN, regular 10 s grid, `tl = tsl - rsl`)
- [x] Provenance: run manifest, logs, sub-link id map
- [x] Tests on synthetic data and on the real file
- [x] Ingest all 92 days of JJA 2015
- [x] Ingest gauges (SMHI GbgA 15 min, City network 1 min) into a common gauge format
- [x] Ingest SMHI hourly weather (temperature, humidity, pressure, wind)
- [x] Ingest radar (dBZ -> R with the file's Z = 200 R^1.5), one file per month
- [x] Acceptance of the reference data against the OpenMRG paper (`_archive/docs/acceptance/step2_acceptance.md`, local only)
- [x] Radar rain along each link path (geometry step; `docs/methods.md` sections 4-5)
- [x] Synchronize CML, gauges and radar on common time steps (L3 aggregates at 1 min - 1 h)
- [ ] Baseline estimation [8]
- [x] Rain events, dry periods and per-link wet/dry state from the references (radar + gauges)
- [ ] Wet/dry classification from the CML signal itself
- [x] Per-stage config hashes (stage_hash) for all outputs
- [x] QC calibration (`_archive/docs/data_notes/openmrg_qc_calibration.md`, local only) and L2 QC flags (qc_stats, qc_flags, qc_status)
- [x] QC pilot on five days, gallery reviewed (`_archive/docs/acceptance/step3_pilot.md`, local only)
- [x] QC flags on the whole period + sub-link status table (`quality_status.csv`)
- [x] Check TSL behaviour (toggling, ATPC-like TSL increase during rain) before baseline estimation (TSL classes in the inventory)
- [ ] Noise and outlier handling beyond QC flags (level shifts, pair mismatch)
- [ ] Link selection for Stage 1 and Stage 2
- [x] Train / val / test split of rain events and train / test split of links (frozen)
- [x] Versioned release `data/release/openmrg_v1` and access API (`cmlpipe.api.CMLData`)
- [ ] OpenRainER adapter (sub-links named `channel1` / `channel2`)
- [ ] Check whether OpenMesh has a rain reference before using it

**2. Simulator (known ground truth)**

- [ ] Rain fields with controlled cell size
- [ ] Exact path integral of a·R(x,t)^b along each link
- [ ] Known WAA model (wetting and drying)
- [ ] Noise and quantization (RSL 0.3 dB, TSL 1 dB)

**3. Stage 1: WAA**

- [ ] Find links with a drying tail after rain ends
- [ ] Fit literature WAA models [1], [9] per link
- [ ] Learn a dynamic WAA model with weak-form SINDy (rain as input), first in simulation
- [ ] Compare the learned models across links (same form, different coefficients?)
- [ ] Drying time from SINDy vs drying time measured directly from tails

**4. Stage 2: path law f(L, R̄, Z)**

- [ ] Simulation: exact integral vs the linear law a·R̄^b·L
- [ ] Learn f(L, R̄, Z) with PySR in simulation
- [ ] Real data with radar along the path as reference (WAA removed first)

**5. Rain estimation and evaluation**

- [ ] Invert the combined model (A_rain + WAA) to get rain
- [ ] Normalized bias and normalized RMSE vs gauges and radar, split by link length
- [ ] Define and compute the false-rain-after-event metric
- [ ] Compare with power law, literature WAA models and the RNN baseline [11]
- [ ] Bootstrap stability of learned terms across links and events
- [ ] Final report

</details>
<!-- /AUTO:progress -->

Data licence: CC BY-SA 4.0 ([`LICENSE-DATA.md`](LICENSE-DATA.md)); code: MIT (`../LICENSE`).
