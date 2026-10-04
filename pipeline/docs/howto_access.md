# How to access the data

Two places hold the processed data:

- **The release** `data/release/openmrg_v1/` (read-only, versioned): tables as Parquet + CSV, time series as Zarr,
  `catalog.json`, `DATA_DICTIONARY.md`. Use it for analysis, through `cmlpipe.api.CMLData`.
- **The work folder** `data/work/openmrg/`: what the pipeline stages write (NetCDF per day, tables, reports). Use it
  only when developing the pipeline.

Setup once: `conda env create -f environment.yml`, `conda activate cmlpipe`, `pip install -e .` from the repository
root. On this machine the environment is `C:\Users\yoni\miniconda3` (the Anaconda `python` has no netCDF4).

**Do not** open many NetCDF or Zarr files with `xarray.open_mfdataset` (dask is not installed, and xarray inspects
every variable of every file, which is very slow on Windows). CMLData opens each Zarr store once and reads only the
chunks you slice.

Plain-language overview first: [`DATA_GUIDE.md`](DATA_GUIDE.md). To look at the data interactively:
`python -m cmlpipe explore` (or double-click `start_explorer.bat`).

## Python

```python
from cmlpipe.api import CMLData

d = CMLData.open("data/release/openmrg_v1")

# --- find links and events --------------------------------------------------------------
d.links(length_km=(1, 3), frequency_band=">=25GHz", status="good", split="train", has_gauge_within_m=1000)
d.sublinks(status="degraded")[["sublink_id", "status_reason"]]
d.events(significant=True, split="test")[["event_id", "start", "end", "rain_class"]]
d.event_runs("E20150727_2140")
d.dry_periods()

# --- time series (xarray; times are UTC interval starts, [start, end)) ---------------------
ds = d.cml([10227, "10001_2"], "2015-07-28", "2015-07-29", res="5min")      # dims (time, sublink_id)
raw = d.cml(10227, "2015-07-28 14:00", "2015-07-28 18:00", res="10s", qc="flags")   # values + QC bitmask
r = d.radar_along([10227], "2015-07-28", "2015-07-29", res="15min")         # dims (time, cml_id)
g = d.gauge(["Torsl", "GbgA"], "2015-07-28", "2015-07-29", res="15min")     # dims (time, gauge_id)
s = d.link_state([10227], "2015-07-28", "2015-07-29")                        # 0 dry, 1 wet, -1 unknown

# --- one link, everything aligned, as one DataFrame -----------------------------------------
df = d.frame(10227, "2015-07-28 12:00", "2015-07-28 20:00", res="5min")
df.attrs["gauge"], df.attrs["notes"]            # which gauge was used; sources left out at this resolution
w = d.event_window("E20150727_2140", cml_id=10227, pad="3h")   # + in_event, run_id, in_wet_core, in_dry_period

# --- documentation from the catalog ----------------------------------------------------------
d.describe("R_max")
```

QC options of `cml(..., res="10s")`: `qc="mask"` (default; samples with a blocking flag are NaN), `"raw"`, `"flags"`
(raw values plus the `flags` bitmask; decode with `cmlpipe.qc.flags.decode_flags`). At 1 min and coarser the
statistics already exclude blocked samples; `qc="flags"` adds `flags_or`. 10 s requests longer than 7 days (or more
than 50 sub-links for more than 2 days) need `allow_large=True`.

### Event x link pairs for modelling

`event_pairs()` gives one row per (event, link) with the columns of `event_links`, `event_links_flags_extra` and the
event significance, after a selection and three quality filters:

```python
pairs = d.event_pairs()                                  # significant events, all splits, default filters
pairs, dropped = d.event_pairs(split=("test", "train"), return_dropped=True)   # unseen events x train links
dropped.drop_reason.value_counts()                       # why each pair was left out
```

| argument | default | effect |
|---|---|---|
| `significant` | `True` | `True` / `False` / `None` (all events); a selection, not a drop |
| `split` | `None` | event split (`"train"`, `"val"`, `"test"`) or a tuple (event split, link split); a selection |
| `usable_only` | `True` | drop (`unusable`) when no sub-link has status good / degraded |
| `max_flatline_share` | `0.5` | drop (`flatline`) when `frac_flatline_event_max` is above it; `None` disables |
| `min_valid_fraction` | `0.8` | drop (`low_valid`) when the best usable sub-link has a smaller valid share of `tl_mean` bins in the event; `None` disables |

On openmrg_v1 with the defaults: 12740 significant pairs, 12289 kept, 451 dropped (280 `unusable; low_valid`,
163 `low_valid`, 8 `flatline`: six links in E20150828_2145, 10329 in E20150603_1500, 10076 in E20150829_1445).

### A bulk matrix: 5 min tl_mean of all train links in all train events

```python
import pandas as pd
links = d.links(split="train", status="good")
events = d.events(significant=True, split="train")
blocks = []
for e in events.itertuples():
    x = d.cml(links.cml_id, e.start, e.end, res="5min", variables=["tl_mean"]).tl_mean.to_pandas()
    blocks.append(x.assign(event_id=e.event_id))
matrix = pd.concat(blocks)                     # rows: (time, event), columns: sublink_id
```

### Without cmlpipe

```python
import pandas as pd, xarray as xr
links = pd.read_parquet("data/release/openmrg_v1/tables/links.parquet")
ds = xr.open_dataset("data/release/openmrg_v1/timeseries/cml_5min.zarr", engine="zarr", chunks=None, consolidated=True)
```

## PowerShell

```powershell
$py = "$env:USERPROFILE\miniconda3\python.exe"
& $py -m cmlpipe stages
& $py -m cmlpipe run all --workers 4          # everything up to date is skipped
& $py -c "from cmlpipe.api import CMLData; d = CMLData.open('data/release/openmrg_v1'); print(d.links(status='degraded')[['cml_id','status_sublink_1','status_sublink_2']])"
Start-Process data\work\openmrg\reports\run_report.html
```

## Work folder (development)

`cmlpipe.api.access.load_link_window(cfg, cml_id, res, start, end, pad_before, pad_after)` reads the daily L3 files of
one window; `cmlpipe.api.explore` has find_links / find_events / show_link / show_event on the work tables.
