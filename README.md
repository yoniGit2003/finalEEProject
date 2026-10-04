# Discovering the Laws of Microwave-Link Attenuation

Data-driven discovery of rain-induced and wet-antenna attenuation laws in commercial microwave links (CML).
Cellular backhaul links lose signal when it rains; we learn readable equations for how that loss depends on rain,
link length and the wetting and drying of the antennas (weak-form SINDy, PySR), using open CML, radar and gauge data.
Final-year EE project, CellEnMon lab, by Yonatan Haber and Yehonatan Tsitseluk. **Research in progress.**

## Where to start

- [`research/`](research/): the research (models, experiments, figures). Uses the data only through a fixed API.
- [`pipeline/`](pipeline/): how the data was prepared (code, explorer GUI, tests). Start with the plain-language
  data guide [`pipeline/docs/DATA_GUIDE.md`](pipeline/docs/DATA_GUIDE.md).
- [`docs/`](docs/): project documents (proposal, report material).

## Get the data

The processed data release `openmrg_v1` (OpenMRG, Gothenburg, summer 2015: 364 links at 10 s, radar, gauges,
weather) is not in git. It must end up in `data/release/openmrg_v1/` (so that `data/release/openmrg_v1/catalog.json`
exists). Either:

- **Option A, rebuild it** from the public raw OpenMRG data (Zenodo, https://doi.org/10.5281/zenodo.7107689):
  unpack it into `data/raw/openmrg/` and run `python -m cmlpipe run all` (details: [`pipeline/README.md`](pipeline/README.md)).
- **Option B, download it**: `openmrg_v1_release.zip` (1.6 GB, CC BY-SA 4.0) from the Releases page,
  https://github.com/yoniGit2003/finalEEProject/releases/tag/data-openmrg-v1 (the repository is private: a GitHub
  login with access is needed). Check it, then unzip it into `data/release/`:

  ```powershell
  (Get-FileHash openmrg_v1_release.zip -Algorithm SHA256).Hash
  # must be 17E9A9048AD0DAE51166C59BA0D3FA760C2E090B4757D70A0B2CDD033AC009F2
  ```

## Quick start (Windows, conda)

```powershell
conda env create -f pipeline/environment.yml     # Python 3.12 + all packages + cmlpipe itself (editable)
conda activate cmlpipe                           # (or, in an existing env: pip install -e "pipeline[dev,explorer]")
python -m cmlpipe explore                        # data explorer in the browser (http://localhost:8501)
```
In Python: `from cmlpipe.api import CMLData; d = CMLData.open("data/release/openmrg_v1")`.

## Licences

Code: MIT ([`LICENSE`](LICENSE)). Data release: CC BY-SA 4.0, derived from OpenMRG (Andersson et al. 2022,
https://doi.org/10.5194/essd-14-5411-2022), with attribution to the original authors
([`pipeline/LICENSE-DATA.md`](pipeline/LICENSE-DATA.md)).
