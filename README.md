# Discovering the Laws of Microwave-Link Attenuation

Final-year EE project, CellEnMon lab (opportunistic environmental sensing).
Team: Yonatan Haber and Yehonatan Tsitseluk.
*[Placeholder: supervisor, project dates.]*

## The question

When it rains, the microwave links that connect cellular towers lose signal, so the phone network can be used as a
rain sensor. The usual model says the rain loss is A_rain = a·R̄^b·L: a power of the average rain rate R̄ along the
link, times the link length L. That is exact only when the rain is uniform along the path, and real links of
different lengths do not always follow it. A second effect, wet-antenna attenuation (WAA), comes from water on the
antenna covers: it does not grow with L and it stays for a while after the rain stops. We ask how the loss really
depends on the link: A_rain = a·R̄^b·f(L, R̄, Z), where uniform rain gives f = L and Z stands for path and link
features. We keep the known physics and learn only the missing parts as readable equations (a gray-box approach):
a dynamic WAA model first, then the correction f, and finally rain estimates from the combined model.

## How the repository is organized

- [`pipeline/`](pipeline/): the data layer. Turns the raw OpenMRG data into a versioned, read-only data release.
- [`research/`](research/): the research code (simulator, WAA, path law, evaluation).
- [`docs/`](docs/): project documents (proposal, report material).
- `data/`: raw data, work files and the release. Not in git; see [`pipeline/README.md`](pipeline/README.md).

## Status

<!-- AUTO:status -->
- **Pipeline (data layer)**: done (release openmrg_v1 built and published; the open items are extensions); 21 / 27 roadmap items. Details: [`pipeline/README.md`](pipeline/README.md)
- **Research**: in progress (not started in code yet); 0 / 18 roadmap items. Details: [`research/README.md`](research/README.md)
<!-- /AUTO:status -->

## Boundary rule

Research code reads data only through `cmlpipe.api.CMLData` and the release pinned in `research/config.yaml`, and the
pipeline never imports research; `python pipeline/tools/check_boundaries.py` checks this.

## Setup

Windows with conda: `conda env create -f pipeline/environment.yml`, then `conda activate cmlpipe`. This installs
Python 3.12, all packages and the `cmlpipe` package. How to get or build the data: [`pipeline/README.md`](pipeline/README.md).

## Licence and citation

Code: MIT ([`LICENSE`](LICENSE)). Data release: CC BY-SA 4.0 ([`pipeline/LICENSE-DATA.md`](pipeline/LICENSE-DATA.md)).
It is derived from OpenMRG: Andersson, J. C. M., Olsson, J., van de Beek, R. (C. Z.), and Hansryd, J. (2022),
*Earth Syst. Sci. Data*, 14, 5411-5426, https://doi.org/10.5194/essd-14-5411-2022 (data:
https://doi.org/10.5281/zenodo.7107689).
