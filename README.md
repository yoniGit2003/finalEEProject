# Discovering the Laws of Microwave-Link Attenuation

Final-year EE project, CellEnMon lab. Yonatan Haber and Yehonatan Tsitseluk.

## The project *(in progress)*

Rain along the path of a commercial microwave link (CML) attenuates its signal, so cellular backhaul networks can
measure rain. The measured attenuation is A(t) = TSL(t) − RSL(t) − baseline = A_rain(t) + δ(t), with
A_rain = ∫₀ᴸ a·R(x,t)^b dx ≈ a·R̄^b·L (exact only for uniform rain) and δ(t) the non-rain part, mainly wet-antenna
attenuation (WAA). We ask how A_rain really depends on link length, A_rain = a·R̄^b·f(L, R̄, Z), and how WAA wets and
dries, and we want readable equations (gray-box): weak-form SINDy for the WAA dynamics, PySR for f, both tested
first in a simulator with known truth, then inverted to estimate rain on held-out events.

> *[Placeholder: research question in our own words, current status, main results.]*
>
> *[Placeholder: supervisor, dates, links to the report / thesis.]*

Background: [`docs/proposal.md`](docs/proposal.md). Plan and progress: [`docs/roadmap.md`](docs/roadmap.md).

## Data layer *(done)*

**The release** `data/release/openmrg_v1/` is a read-only, versioned copy of the processed OpenMRG data
(Gothenburg, 1 June – 31 August 2015): 364 links / 728 sub-links at 10 s and aggregated to 1 / 5 / 15 / 60 min,
QC flags, radar along every link path, City and SMHI gauges, weather, rain events, dry periods, per-link reference
states, a quality status and a frozen train / test split. Parquet tables + Zarr time series (about 1.6 GB).
**Start with [`docs/DATA_GUIDE.md`](docs/DATA_GUIDE.md)** (plain-language guide, also in the release folder);
every column: [`docs/data_dictionary.md`](docs/data_dictionary.md); definitions: [`docs/methods.md`](docs/methods.md).

**Setup** (Windows, conda; no dask):

```powershell
conda env create -f environment.yml; conda activate cmlpipe; pip install -e ".[dev,explorer]"
```

**Build it.** Download OpenMRG v1.1 (https://doi.org/10.5281/zenodo.7107689) and unpack it so that these exist
(`data/` is not in git and `data/raw/` is never modified):

```
data/raw/openmrg/cml/cml.nc, cml/cml_metadata.csv          CML (about 4.6 GB) and link metadata
data/raw/openmrg/gauges/city/CityGauges-2015JJA.csv (+ -metadata.csv)
data/raw/openmrg/gauges/smhi/GbgA-71420-2015JJA.csv (+ -metadata.csv)
data/raw/openmrg/radar/radar.nc
```

```powershell
python -m cmlpipe run all --workers 4      # every stage; up-to-date stages are skipped; work files in data/work/openmrg
python -m cmlpipe stages                   # the stages in order (a single one: python -m cmlpipe run events)
```

The last stage writes the release; an existing release is never overwritten (use `--force`, or `--version v2`).

**Open it** from Python with [`CMLData`](src/cmlpipe/api/cmldata.py) (examples:
[`docs/howto_access.md`](docs/howto_access.md), [`examples/quickstart.ipynb`](examples/quickstart.ipynb)):

```python
from cmlpipe.api import CMLData
d = CMLData.open("data/release/openmrg_v1")
df = d.event_window("E20150727_2140", cml_id=10316, pad="3h")   # CML + radar + gauges + weather + states
```

**Look at it** in the browser with the explorer: double-click `start_explorer.bat`, or

```powershell
python -m cmlpipe explore                  # opens http://localhost:8501 ; stop with Ctrl+C
```

## Repository map

```
src/cmlpipe/     the data layer: pipeline stages, release builder, CMLData API, plotting (viz)
apps/explorer/   the explorer GUI (Streamlit + Plotly; reads only through CMLData); smoke.py checks it
research/        research code (baselines, WAA, models, notebooks, figures, thesis); pinned release in config.yaml
configs/         openmrg.yaml: every setting of the pipeline
docs/            proposal, data guide, data dictionary, methods, how-to, roadmap, papers
examples/        quickstart notebook        tools/   README refresher, boundary check
tests/           python -m pytest           data/    raw, work and release data (not in git)
```

**Boundary rule.** Research code uses the data only through `cmlpipe.api.CMLData` (and `cmlpipe.viz`) and a pinned
release version; it never imports cmlpipe internals, and cmlpipe and the explorer never import `research`.
`python tools/check_boundaries.py` checks it.

## Data citation and licence

**Licences:** code MIT ([`LICENSE`](LICENSE)); the data release is CC BY-SA 4.0 with attribution to the OpenMRG
authors ([`LICENSE-DATA.md`](LICENSE-DATA.md)).

Data: Andersson, J. C. M., Olsson, J., van de Beek, R. (C. Z.), and Hansryd, J. (2022): OpenMRG: Open data from
Microwave links, Radar, and Gauges for rainfall quantification in Gothenburg, Sweden, *Earth Syst. Sci. Data*, 14,
5411–5426, https://doi.org/10.5194/essd-14-5411-2022. Dataset: https://doi.org/10.5281/zenodo.7107689 (CC BY-SA 4.0).
