# Pipeline roadmap

What is done and what is left in the DATA LAYER (pipeline/). The research tasks (simulator, WAA, path law, rain
estimation) are in `research/ROADMAP.md`. Tick an item by changing `- [ ]` to `- [x]`: the status block of
`pipeline/README.md` and the one-line summary in the root `README.md` are generated from this file on every commit
(`pipeline/tools/update_readme.py`), so this is the only place where pipeline progress is edited by hand.

Status: done (release openmrg_v1 built and published; the open items are extensions)

## 1. Pipeline (data)
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
