# Data guide: openmrg v1

Generated from `catalog.json` (code 211f942743ac, built 2026-10-02T17:42:28Z). Column-level
details: `DATA_DICTIONARY.md`. Methods: `docs/methods.md` in the cmlpipe repository.
**Licence:** CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/), derived from OpenMRG (Andersson et al. 2022, https://doi.org/10.5281/zenodo.7107689); see LICENSE.txt.

## 1. What this is

- OpenMRG: commercial microwave links (CML) in Gothenburg, Sweden, 2015-06-01 to 2015-09-01 (exclusive), UTC.
- 364 links = 728 sub-links (one per direction), received and transmitted signal level every 10 s.
- References: NORDRAD radar (5 min, 2 km) along every link path, 10 City gauges (1 min), SMHI gauge GbgA (15 min), hourly SMHI weather.
- Added by cmlpipe: QC flags, statistics at 1 / 5 / 15 / 60 min, rain events and dry periods (from the references
  only, never from CML), per-link wet / dry reference state, a quality status and a frozen train / test split.
- Source: Andersson, J. C. M., Olsson, J., van de Beek, R. (C. Z.), and Hansryd, J. (2022): OpenMRG: Open data from Microwave links, Radar, and Gauges for rainfall quantification in Gothenburg, Sweden, Earth Syst. Sci. Data, 14, 5411-5426, https://doi.org/10.5194/essd-14-5411-2022. Data: https://doi.org/10.5281/zenodo.7107689 (CC BY-SA 4.0).

## 2. Folder map

```
openmrg_v1/
  README.md, DATA_GUIDE.md      short intro; this guide
  LICENSE.txt                   data licence: CC BY-SA 4.0, attribution to the OpenMRG authors
  DATA_DICTIONARY.md            every column and variable with type, unit and description
  catalog.json                  the same, machine-readable (+ conventions, QC sets, stage hashes)
  verification.json             checks run when the release was built
  tables/<name>.parquet + .csv  (same content, two formats):
    links                     one row per link: radio parameters, split, status of both sub-links, completeness, TSL behaviour, nearest gauges, neighbours, radar cells (364 rows)
    sublinks                  one row per sub-link: metadata, completeness (JJA and per month), flag shares, quantization, TSL behaviour, dry-weather noise, quality status (728 rows)
    gauges                    City (10) and SMHI (1) gauges: position, type, completeness, JJA totals, radar pixel comparison (11 rows)
    sites                     physical sites (masts): link ends merged within 10 m (418 rows)
    events                    rain events defined from radar + City gauges (never CML) (89 rows) FROZEN
    event_significance        significance, rain class and wet core of every event (89 rows)
    event_runs                wet runs inside events (before gap merging) (368 rows)
    dry_periods               certified dry periods (>= 6 h) (62 rows) FROZEN
    split_events              FROZEN event split (train / val / test) (89 rows) FROZEN
    split_links               FROZEN link split (train / test / excluded) (364 rows) FROZEN
    event_links               one row per (event, link): radar along the path, matched gauge, CML data quality (32396 rows) FROZEN
    event_links_flags_extra   per (event, link): flatline / tsl_step shares and counts in the event (information only) (32396 rows)
    link_gauge_matches        every (link, gauge) pair within 5 km of the path (852 rows)
    network_gaps              intervals when many alive sub-links were missing (collector outages) (93 rows)
    quality_status            quality status per sub-link (good / degraded / unusable) and the reason (728 rows)
    sublink_dry_stats         descriptive CML statistics per sub-link over the certified dry periods (728 rows)
  timeseries/<name>.zarr        (Zarr v2, consolidated; open with xarray):
    cml_10s                   CML at 10 s: tsl, rsl, tl (L1, int16 x 0.1 dB, exact copies) and the L2 QC flag bitmask; dims time=794880, sublink_id=728
    cml_1min                  CML statistics per 1 min bin (L3); dims time=132480, sublink_id=728
    cml_5min                  CML statistics per 5 min bin (L3); dims time=26496, sublink_id=728
    cml_15min                 CML statistics per 15 min bin (L3); dims time=8832, sublink_id=728
    cml_1h                    CML statistics per 1 h bin (L3); dims time=2208, sublink_id=728
    radar_along_links_5min    radar rain rate along every link path (floored at 0 dBZ), 5 min; dims time=26496, cml_id=364
    radar_along_links_15min   radar rain rate along every link path, 15 min; dims time=8832, cml_id=364
    radar_along_links_1h      radar rain rate along every link path, 1 h; dims time=2208, cml_id=364
    radar_at_gauges_5min      radar rain rate at every gauge (pixel and 3 x 3 mean), 5 min; dims time=26496, gauge_id=11
    radar_at_gauges_15min     radar rain rate at every gauge, 15 min; dims time=8832, gauge_id=11
    radar_at_gauges_1h        radar rain rate at every gauge, 1 h; dims time=2208, gauge_id=11
    gauges_city_1min          City gauges, 1 min amounts and rates; dims time=132480, gauge_id=10
    gauges_city_5min          City gauges, 5 min sums (strict: NaN with a missing minute); dims time=26496, gauge_id=10
    gauges_city_15min         City gauges, 15 min sums; dims time=8832, gauge_id=10
    gauges_city_1h            City gauges, 1 h sums; dims time=2208, gauge_id=10
    gauges_smhi_15min         SMHI gauge GbgA, 15 min (native); dims time=8832, gauge_id=1
    gauges_smhi_1h            SMHI gauge GbgA, 1 h sums; dims time=2208, gauge_id=1
    met_1h                    SMHI weather at GbgA, hourly (instantaneous / as in the file; see each variable); dims time=2208, gauge_id=1
    network_reference_5min    network reference series: radar over the link domain + City gauges, state and event id per 5 min bin; dims time=26496
    link_reference_state_5min per-link reference state (wet / dry / unknown) and reasons, from references only; dims time=26496, cml_id=364
```

## 3. Conventions that cause mistakes

- **time** (catalog): UTC, tz-naive; every value is labelled by the START of its interval [t, t + dt) ('interval_start'); radar stamps are instants assigned to [t, t + 5 min)
- **keys** (catalog): cml_id (int, link); sublink_id (str '<cml_id>_<1|2>'; sub-link 1 transmits site_0 -> site_1); gauge_id (City station name or 'GbgA'); event_id ('E' + start, e.g. E20150728_1455); run_id ('<event_id>_Rnn')
- **units** (catalog): dB / dBm for signal levels, mm for amounts, mm/h for rates, m for distances (EPSG:3006), GHz
- **missing** (catalog): NaN (float), fill value -32768 (int16, decoded to NaN), 0 for counts / flags
- **splits** (catalog): split_events.csv and split_links.csv are frozen; train on (train events x train links)
- **Signals**: `tsl` transmitted and `rsl` received level (dBm); `tl = tsl - rsl` total loss (dB). Rain makes `tl`
  go UP. Attenuation needs a baseline: A(t) = tl(t) - baseline (no baseline is applied in this release).
- **Sub-links**: `sublink_id = "<cml_id>_1"` is the raw direction **A** and transmits site_0 -> site_1;
  `"<cml_id>_2"` is direction **B** (site_1 -> site_0). The two directions of a link differ in frequency (up to 1.5 GHz).
- **10 s values** (`cml_10s`) are stored as int16 with scale 0.1 dB and fill -32768; xarray decodes them to float dB
  with NaN. Nothing is edited: bad samples are only flagged.
- **QC flags** are bits of a uint16 (`flags` at 10 s, `flags_or` = OR over the bin): missing_rsl = 1, missing_tsl = 2, out_of_range = 4, rx_floor = 8, spike = 16, tsl_step = 32, flatline = 64, network_gap = 128.
  Blocking sets (a sample with one of these bits is not "valid"): rsl missing_rsl, out_of_range, rx_floor, spike, network_gap;
  tsl missing_tsl, out_of_range, network_gap; tl the union. tsl_step and flatline are information only.
- **Valid** sample = value present and no blocking flag. **The 80% rule**: a statistic of a bin (mean, median, min,
  max, std over the 10 s samples) is NaN unless at least 80% of the bin's 10 s samples are valid; `n_valid_<v>`
  always tells how many were. Radar bins need 60% valid stamps; gauge sums are NaN if any native step is missing.
- **Radar** rain rates below 0 dBZ are set to 0 (no phantom drizzle). `R_mean` = path-weighted mean, `R_max` = wettest cell.
- **NaN** always means "unknown", never zero rain or zero loss.

## 4. Open it without cmlpipe

```python
import pandas as pd, xarray as xr
R = "data/release/openmrg_v1"
links = pd.read_parquet(f"{R}/tables/links.parquet")                 # any table (or read_csv on the .csv)
cml = xr.open_dataset(f"{R}/timeseries/cml_5min.zarr", engine="zarr", chunks=None, consolidated=True)
tl = cml.tl_mean.sel(sublink_id="10316_1", time=slice("2015-07-27 20:00", "2015-07-28 06:00")).to_series()
radar = xr.open_dataset(f"{R}/timeseries/radar_along_links_5min.zarr", engine="zarr", chunks=None).R_mean.sel(cml_id=10316)
```
Select before you load (`.sel(...)` then `.values`): the stores are lazy and only the needed chunks are read.

## 5. Open it with cmlpipe (`pip install -e .` in the repository)

```python
from cmlpipe.api import CMLData
d = CMLData.open("data/release/openmrg_v1")
d.links(length_km=(1, 3), status="good", split="train")         # filter links
d.events(significant=True, split="train")                         # events + significance
df = d.frame(10316, "2015-07-27 18:00", "2015-07-28 12:00", res="5min")   # one link, everything aligned
d.event_window("E20150727_2140", cml_id=10316, pad="3h")           # frame + event / run / dry markers
d.event_pairs(split=("train", "train"))                            # event x link pairs for modelling
d.describe("tl_mean")                                              # unit + description of any variable
```

## 6. Definitions in short

- **Wet 5 min bin (network)**: >= 1% of the radar cells near the links have >= 0.1 mm/h, or a City gauge has >= 0.1 mm.
- **Event**: wet bins merged across gaps shorter than 2 h; events shorter than 15 min dropped. `event_runs` are the wet runs before merging.
- **Significant event**: mean City total >= 1.0 mm, or radar domain total >= 1.0 mm, or a 5 min City rate >= 5.0 mm/h; otherwise drizzle. Heavy: total >= 10.0 mm or rate >= 30.0 mm/h.
- **Dry period**: certified dry bins (radar valid, no wet cell, no wet gauge) more than 30 min from any wet bin, in runs >= 6 h. Use these for baselines.
- **Per-link state** (`link_reference_state_5min`): wet if `R_max` >= 0.1 mm/h or a gauge within 1000 m has >= 0.1 mm; dry if `R_max` = 0 and every gauge within 1000 m = 0 in the bin and +-30 min; else unknown.
- **Split (frozen)**: events {'train': 0.6, 'val': 0.2, 'test': 0.2} stratified by intensity; links {'train': 0.8, 'test': 0.2} stratified by length class x band (< / >= 25 GHz), parallel links together, seed 42. Train on train events x train links; the link split is not spatially independent.
- **Status per sub-link**: good / degraded / unusable (completeness, flag shares, dry noise); `links.status_best` is the better of the two.

## 7. Known limitations

- The radar is gauge-adjusted (not independent of the gauges) and below the gauge totals at most gauge pixels.
- Many links are short and lie in one 2 km radar cell; only some links have a gauge within 1 km.
- Link length and frequency are correlated; most events are drizzle.
- Not in this release: baseline, attenuation A(t), CML wet / dry classification, wet-antenna correction, rain retrieval.

## 8. Which file do I need?

| I want | Use |
|---|---|
| one link, all summer | `timeseries/cml_1h.zarr` or `cml_15min.zarr`, `.sel(sublink_id=["<id>_1", "<id>_2"])` |
| all links, one event | `tables/events` for start / end, then `cml_5min.zarr` `.sel(time=slice(start, end))` |
| QC flags | `cml_10s.zarr` variable `flags` (or `flags_or` in the aggregates); bits in section 3 |
| radar along a path | `radar_along_links_<res>.zarr`, dim `cml_id` |
| gauges near a link | `tables/link_gauge_matches` (distance to path), then `gauges_city_*` / `gauges_smhi_*` |
| train / test split | `tables/split_events`, `tables/split_links` (frozen) |
| dry periods for a baseline | `tables/dry_periods` (+ `sublink_dry_stats` for per-sub-link dry statistics) |
| event x link quality | `tables/event_links`, `event_links_flags_extra` (or `CMLData.event_pairs`) |
