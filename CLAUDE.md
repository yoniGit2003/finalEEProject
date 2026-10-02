# finalEEProject

Final-year EE project in the CellEnMon lab: **data-driven discovery of rain-induced and wet-antenna attenuation laws in commercial microwave links (CML)**.
Students: Yonatan Haber and Yehonatan Tsitseluk. Project documents are written in Hebrew (RTL) with English technical terms.
Working title: "Discovering the Laws of Microwave-Link Attenuation" (Hebrew: גילוי משוואות מנתונים להנחתת גשם ולהרטבת אנטנה בלינקים סלולריים).

## Project documents
These two files are imported, so their full text is loaded at the start of every session:

@docs/proposal.md
@docs/openmrg_summary.md

- `docs/proposal.md` is the project proposal, revised after the guide's comments (same content as `docs/proposal.docx`): goals, stages, deliverables, success criteria, references [1]-[13].
- `docs/openmrg_summary.md` covers the OpenMRG dataset: variables, units, quantization, missing values, gauges, radar.
- `docs/papers/Andersson2022_OpenMRG.pdf` is the OpenMRG data paper (ESSD 2022, 16 pages). Not imported; `docs/openmrg_summary.md` already holds its facts. Open the PDF only for figures or exact wording.
- `docs/Rain_project.pdf` is not imported (large PDF). Read it on demand: our report "Introduction to SINDy" (March 2026) with SINDy/PySR tutorials (RLC circuit, heat-equation PDE, Term Robustness Index), Section 5 on CML rain attenuation and goals, and Appendix I on continuous-time WAA modelling.

## Physics and notation (keep consistent in code and text)
- Measured attenuation: A(t) = TSL(t) - RSL(t) - baseline, and A(t) = A_rain(t) + δ(t), where δ(t) collects non-rain terms (WAA and baseline residual). WAA is the first component of δ(t) we model.
- Path integral: A_rain = ∫₀ᴸ a·R(x,t)^b dx  ≈  a·R̄^b·L (ITU power law; exact only for uniform rain). R̄ = path-averaged rain rate, L = link length in km, a and b depend on frequency and polarization (ITU-R P.838-3).
- Extended model: A_rain = a·R̄^b·f(L, R̄, Z). Uniform rain gives f = L. Z = path/link features (rain non-uniformity along the path, frequency, polarization, time resolution).
- A constant WAA looks like f = L + L0, so WAA must be isolated before learning f.
- WAA is local to the antennas, dynamic with memory (wetting and drying), and does not scale with L.

## Plan
1. **Pipeline**: extract, clean, synchronize and process raw CML data and metadata; baseline estimation; rain-event detection; noise/outlier handling; link selection.
2. **Simulator** with known ground truth: rain fields with controlled cell size, exact path integral per link, known WAA model, noise and quantization. Test every method here before real data.
3. **Stage 1 (WAA)**: identify links with a drying tail; the tail after rain ends has zero path attenuation, so what remains is WAA (do NOT select links by how well the power law fits, that is circular). Fit literature models [1], [9]; learn a dynamic model with **weak-form SINDy** (rain as input); compare across links.
4. **Stage 2 (path law)**: learn f(L, R̄, Z) with **PySR**; compare exact integral vs the linear law in simulation, then real data with radar along the path as reference.
5. **Rain estimation**: invert the combined model on events not used for training.
- Approach is gray-box: keep known physics, learn the missing relations. We want readable equations, not just a black-box predictor. A black-box RNN [11] is only a performance baseline.

## Success criteria
- Simulation: recover the correct terms and coefficients as a function of noise, quantization and link length.
- Real data: normalized bias and normalized RMSE of reconstructed rain vs gauges and radar, split by link length, compared with literature models; amount of false rain after events ends (metric to define).
- Physical consistency: drying time learned by SINDy vs drying time estimated directly from tails; stability of terms across links and events (bootstrap).
- Rain reference: radar along the path is preferred over a point gauge.

## Data (`data/`, not tracked in docs)
OpenMRG (Gothenburg, JJA 2015) is in `data/raw/openmrg/`: `cml/cml.nc` (~4.6 GB), `cml/cml_metadata.csv`, `gauges/city/`, `gauges/smhi/`, `radar/radar.nc`, `Torslanda-event-2015-07-28/`, `readme.txt`. Planned additional datasets: **OpenRainER** (Emilia-Romagna, Italy; https://zenodo.org/records/14731404) and **OpenMesh** (check that a rain reference exists before committing to it). Follow OpenSense conventions.
Key facts (details in `docs/openmrg_summary.md`):
- 728 sub-links (364 links), 10 s sampling, `rsl` and `tsl` in dBm, RSL quantized 0.3 dB, TSL 1 dB, ATPC off. Missing value in the NetCDF is 1e10, so convert to NaN. Time is UTC seconds since 1970.
- `cml.nc` is large: never load it whole. Slice by time or sublink (xarray with chunks, or netCDF4 slices). In the example script the array order is [time, sublink]; check with `ncdump -h` or `nc.variables['rsl'].dimensions` before assuming.
- `cml_metadata.csv` columns: Link, Direction, Sublink, Near/Far Latitude/Longitude_DecDeg, Frequency_GHz, Polarization, Length_km, Near/FarAntennaDiameter_m. Join to the data via `Sublink`.
- Radar: 5 min, 2 km grid, pseudo-dBZ decoded by netCDF4 (dBZ = 0.4*value - 30); Z = 200 R^1.5 (SMHI summer). Gauges: SMHI GbgA 15 min (`Pvol_mm`), City network 1 min (10 gauges).

## Code and environment
- Windows machine, PyCharm project at `D:\UniShit\finalEEProject`. Package `src/cmlpipe/` (io, qc, geometry, aggregate, events, inventory, viz, release, api), installed with `pip install -e .`; CLI `python -m cmlpipe run <stage>` (`pipeline/run.py` is a deprecated shim); config `configs/openmrg.yaml`; tests `tests/` (`python -m pytest`). Final definitions of every product: `docs/methods.md`; data dictionary: `docs/data_dictionary.md`; access: `docs/howto_access.md` and `cmlpipe.api.CMLData`. Development checks, acceptance reports and calibration notes are in `_archive/` (gitignored, local only); history before step 8 is on branch `dev-history` (tag `pre-step8`).
- Explorer GUI: `apps/explorer/` (Streamlit + Plotly), started with `python -m cmlpipe explore` or `start_explorer.bat`; it reads only through `CMLData`, figures live in `cmlpipe.viz.explorer`; `python apps/explorer/smoke.py` checks it. Plain data guide: `docs/DATA_GUIDE.md` (generated from the release catalog by `python -m cmlpipe guide`).
- `research/` holds research code. It may use only `cmlpipe.api` (and `cmlpipe.viz`) and the release pinned in `research/config.yaml`; cmlpipe and apps never import research (`python tools/check_boundaries.py`).
- Stack: Python, NumPy, pandas, xarray, NetCDF (netCDF4), pycomlink, poligrain, PySR (needs Julia), SINDy (pysindy). Add real dependencies to `environment.yml` and `pyproject.toml` as we use them. No dask.
- Hardware: personal computers and lab compute.
- The example script `data/raw/openmrg/radar/example_read_radar_nc.py` has a hard-coded Linux path; do not reuse it as is.

## Working rules
- Write to the student's level: explain what a step does and why. Keep code simple and commented; prefer scripts we can rerun over one-off notebooks.
- Do not modify anything in `data/raw/`. Pipeline outputs go to `data/work/openmrg/`; the versioned, read-only release to `data/release/openmrg_<version>/` (rebuild only with `--force`). `events.csv` and the frozen split (`split_events.csv`, `split_links.csv`) must not change unless explicitly asked (`--force-split`).
- Keep notation identical to the proposal (A_rain, R̄, L, f(L, R̄, Z), δ(t)).
- Verify literature claims against the references in the proposal before stating them; do not invent citations.
- Answer in English unless asked for Hebrew. When writing Hebrew documents, keep English terms in Latin script and check RTL ordering.
