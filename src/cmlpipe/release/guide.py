"""DATA_GUIDE.md: one plain-language page about a release, for a person or a program that knows nothing about cmlpipe.

Everything that is also in the data dictionary (conventions, products, descriptions, row counts, dims, flag bits,
blocking sets) is taken from catalog.json, the same source as DATA_DICTIONARY.md, so the two cannot disagree.
Rule thresholds come from the config that built the release. The guide is written next to the release
(<release>/DATA_GUIDE.md) and to docs/DATA_GUIDE.md by the release stage, or on its own with
`python -m cmlpipe guide` (adds only that one file; the data are not touched).
"""

from __future__ import annotations

import json
import os
import stat
from pathlib import Path

from cmlpipe.config import Config
from cmlpipe.release import metadata


def _rules(cfg: Config) -> dict:
    """The few config values the guide quotes."""
    ev, a, ls = cfg["events"], cfg["aggregate"], cfg["link_state"]
    return {"min_valid": a["min_valid_fraction"], "radar_min_valid": a["radar"]["min_valid_fraction"],
            "wet_mmh": ev["wet_threshold_mmh"], "wet_area": ev["min_wet_area_fraction"],
            "gauge_wet": ev["city_gauge_wet_amount_mm"], "gap_h": ev["min_dry_gap_h"],
            "dry_margin": ev["dry_margin_min"], "dry_h": ev["dry_min_duration_h"],
            "ls_radar": ls["radar_wet_mmh"], "ls_gauge": ls["gauge_wet_mm"], "ls_dist": ls["gauge_max_dist_m"],
            "ls_margin": ls["margin_min"], "seed": cfg["split"]["seed"],
            "ev_frac": cfg["split"]["event_fractions"], "link_frac": cfg["split"]["link_fractions"],
            "f_split": cfg["split"]["frequency_split_ghz"], "evl_valid": cfg["event_links"]["min_valid_fraction"],
            "min_ev": ev["min_event_duration_min"], "sig": cfg["event_significance"]}


def data_guide(catalog: dict, cfg: Config) -> str:
    """The DATA_GUIDE.md text of one release."""
    r = _rules(cfg)
    rel = f"{catalog['name']}_{catalog['version']}"
    t, ts, q = catalog["tables"], catalog["timeseries"], catalog["qc"]
    period = catalog["period"]
    bits = q["flag_bits"]
    flag_line = ", ".join(f"{n} = {v}" for n, v in bits.items())
    n_links = t["links"]["n_rows"]
    n_sub = t["sublinks"]["n_rows"]

    folder = [f"{rel}/",
              "  README.md, DATA_GUIDE.md      short intro; this guide",
              "  LICENSE.txt                   data licence: CC BY-SA 4.0, attribution to the OpenMRG authors",
              "  DATA_DICTIONARY.md            every column and variable with type, unit and description",
              "  catalog.json                  the same, machine-readable (+ conventions, QC sets, stage hashes)",
              "  verification.json             checks run when the release was built",
              "  tables/<name>.parquet + .csv  (same content, two formats):"]
    folder += [f"    {n:<26}{e['description']} ({e['n_rows']} rows)" + (" FROZEN" if e["frozen"] else "")
               for n, e in t.items()]
    folder.append("  timeseries/<name>.zarr        (Zarr v2, consolidated; open with xarray):")
    folder += [f"    {n:<26}{e['description']}; dims " + ", ".join(f"{d}={s}" for d, s in e["dims"].items())
               for n, e in ts.items()]

    conv = "\n".join(f"- **{k}** (catalog): {v}" for k, v in catalog["conventions"].items())
    pct = f"{r['min_valid']:.0%}"
    return f"""# Data guide: {catalog['name']} {catalog['version']}

Generated from `catalog.json` (code {catalog['code_version'][:12]}, built {catalog['created_utc']}). Column-level
details: `DATA_DICTIONARY.md`. Methods: `docs/methods.md` in the cmlpipe repository.
**Licence:** {metadata.DATA_LICENSE_SHORT}.

## 1. What this is

- OpenMRG: commercial microwave links (CML) in Gothenburg, Sweden, {period['start'][:10]} to {period['end_exclusive'][:10]} (exclusive), UTC.
- {n_links} links = {n_sub} sub-links (one per direction), received and transmitted signal level every 10 s.
- References: NORDRAD radar (5 min, 2 km) along every link path, 10 City gauges (1 min), SMHI gauge GbgA (15 min), hourly SMHI weather.
- Added by cmlpipe: QC flags, statistics at 1 / 5 / 15 / 60 min, rain events and dry periods (from the references
  only, never from CML), per-link wet / dry reference state, a quality status and a frozen train / test split.
- Source: {catalog['citation']}

## 2. Folder map

```
{chr(10).join(folder)}
```

## 3. Conventions that cause mistakes

{conv}
- **Signals**: `tsl` transmitted and `rsl` received level (dBm); `tl = tsl - rsl` total loss (dB). Rain makes `tl`
  go UP. Attenuation needs a baseline: A(t) = tl(t) - baseline (no baseline is applied in this release).
- **Sub-links**: `sublink_id = "<cml_id>_1"` is the raw direction **A** and transmits site_0 -> site_1;
  `"<cml_id>_2"` is direction **B** (site_1 -> site_0). The two directions of a link differ in frequency (up to 1.5 GHz).
- **10 s values** (`cml_10s`) are stored as int16 with scale 0.1 dB and fill -32768; xarray decodes them to float dB
  with NaN. Nothing is edited: bad samples are only flagged.
- **QC flags** are bits of a uint16 (`flags` at 10 s, `flags_or` = OR over the bin): {flag_line}.
  Blocking sets (a sample with one of these bits is not "valid"): rsl {', '.join(q['blocking_flags_rsl'])};
  tsl {', '.join(q['blocking_flags_tsl'])}; tl the union. tsl_step and flatline are information only.
- **Valid** sample = value present and no blocking flag. **The {pct} rule**: a statistic of a bin (mean, median, min,
  max, std over the 10 s samples) is NaN unless at least {pct} of the bin's 10 s samples are valid; `n_valid_<v>`
  always tells how many were. Radar bins need {r['radar_min_valid']:.0%} valid stamps; gauge sums are NaN if any native step is missing.
- **Radar** rain rates below 0 dBZ are set to 0 (no phantom drizzle). `R_mean` = path-weighted mean, `R_max` = wettest cell.
- **NaN** always means "unknown", never zero rain or zero loss.

## 4. Open it without cmlpipe

```python
import pandas as pd, xarray as xr
R = "data/release/{rel}"
links = pd.read_parquet(f"{{R}}/tables/links.parquet")                 # any table (or read_csv on the .csv)
cml = xr.open_dataset(f"{{R}}/timeseries/cml_5min.zarr", engine="zarr", chunks=None, consolidated=True)
tl = cml.tl_mean.sel(sublink_id="10316_1", time=slice("2015-07-27 20:00", "2015-07-28 06:00")).to_series()
radar = xr.open_dataset(f"{{R}}/timeseries/radar_along_links_5min.zarr", engine="zarr", chunks=None).R_mean.sel(cml_id=10316)
```
Select before you load (`.sel(...)` then `.values`): the stores are lazy and only the needed chunks are read.

## 5. Open it with cmlpipe (`pip install -e .` in the repository)

```python
from cmlpipe.api import CMLData
d = CMLData.open("data/release/{rel}")
d.links(length_km=(1, 3), status="good", split="train")         # filter links
d.events(significant=True, split="train")                         # events + significance
df = d.frame(10316, "2015-07-27 18:00", "2015-07-28 12:00", res="5min")   # one link, everything aligned
d.event_window("E20150727_2140", cml_id=10316, pad="3h")           # frame + event / run / dry markers
d.event_pairs(split=("train", "train"))                            # event x link pairs for modelling
d.describe("tl_mean")                                              # unit + description of any variable
```

## 6. Definitions in short

- **Wet 5 min bin (network)**: >= {r['wet_area']:.0%} of the radar cells near the links have >= {r['wet_mmh']} mm/h, or a City gauge has >= {r['gauge_wet']} mm.
- **Event**: wet bins merged across gaps shorter than {r['gap_h']} h; events shorter than {r['min_ev']} min dropped. `event_runs` are the wet runs before merging.
- **Significant event**: mean City total >= {r['sig']['city_total_mean_mm']} mm, or radar domain total >= {r['sig']['radar_domain_total_mm']} mm, or a 5 min City rate >= {r['sig']['city_max_5min_rate_mmh']} mm/h; otherwise drizzle. Heavy: total >= {r['sig']['heavy_total_mm']} mm or rate >= {r['sig']['heavy_rate_mmh']} mm/h.
- **Dry period**: certified dry bins (radar valid, no wet cell, no wet gauge) more than {r['dry_margin']} min from any wet bin, in runs >= {r['dry_h']} h. Use these for baselines.
- **Per-link state** (`link_reference_state_5min`): wet if `R_max` >= {r['ls_radar']} mm/h or a gauge within {r['ls_dist']} m has >= {r['ls_gauge']} mm; dry if `R_max` = 0 and every gauge within {r['ls_dist']} m = 0 in the bin and +-{r['ls_margin']} min; else unknown.
- **Split (frozen)**: events {r['ev_frac']} stratified by intensity; links {r['link_frac']} stratified by length class x band (< / >= {r['f_split']} GHz), parallel links together, seed {r['seed']}. Train on train events x train links; the link split is not spatially independent.
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
"""


def write_guide(cfg: Config, release_path: Path, log=None) -> list[Path]:
    """Write the text files of a release from its catalog: DATA_GUIDE.md (release + docs/), LICENSE.txt and README.md.
    The data files and catalog.json are not touched."""
    from cmlpipe.release import build
    release_path = Path(release_path)
    catalog = json.loads((release_path / "catalog.json").read_text(encoding="utf-8"))
    guide = data_guide(catalog, cfg)
    files = {release_path / "DATA_GUIDE.md": guide, cfg.repo_root / "docs" / "DATA_GUIDE.md": guide,
             release_path / "LICENSE.txt": metadata.data_license_text(catalog["name"], catalog["version"]),
             release_path / "README.md": build.readme(catalog)}
    out = []
    for p, text in files.items():
        if p.exists():
            os.chmod(p, stat.S_IWRITE | stat.S_IREAD)  # the release files are read-only
        p.write_text(text, encoding="utf-8")
        if p.parent == release_path:
            os.chmod(p, stat.S_IREAD)
        out.append(p)
    if log is not None:
        log.info("guide: wrote %s", ", ".join(str(p) for p in out))
    return out
