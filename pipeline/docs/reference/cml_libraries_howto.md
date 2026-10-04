# How-to: CML data formats and libraries (xarray, pycomlink, poligrain)

Reference for *how to use* the tools. It is **not** our method: our processing choices are in `docs/proposal.md`.
Source: an example notebook by the OpenSense community (`docs/reference/example_create_8d_OpenMRG_data.ipynb`), plus function signatures checked against **pycomlink 0.6.0** and **poligrain 0.3.1** (September 2026). If installed versions differ, re-check signatures with `help(func)` before relying on this file.

## 1. The OpenSense data format (target format for all our datasets)
Convention paper: Fencl et al. 2023, Open Research Europe, https://open-research-europe.ec.europa.eu/articles/3-169 . poligrain and pycomlink functions expect these names.

**CML dataset** (`xr.Dataset`), as produced for OpenMRG:
```
Dimensions:   (time, sublink_id: 2, cml_id: 364)
Coordinates:  cml_id (int), sublink_id ('sublink_1','sublink_2'), time (datetime64)
              site_0_lat, site_0_lon, site_1_lat, site_1_lon   (cml_id)   WGS84 degrees
              frequency     (sublink_id, cml_id)  MHz   <- NOT GHz
              polarization  (sublink_id, cml_id)  'v' / 'h'
              length        (cml_id)              metres <- NOT km
Data vars:    tsl, rsl  (time, sublink_id, cml_id)  dBm
```
Unit conversions when calling pycomlink: `f_GHz = frequency / 1000`, `f_Hz = frequency * 1e6`, `L_km = length / 1000`.
OpenRainER uses the same format but names the sub-links `'channel1'`, `'channel2'`.

**Gauge dataset**: `rainfall_amount (id, time)` in mm per interval, coords `id`, `time`, `lon`, `lat` (OpenRainER uses `longitude`, `latitude`), optional `location`, `type`, `quantization`.

**Radar dataset**: variable on `(time, y, x)` with **2-D** coords `lon (y, x)`, `lat (y, x)`. poligrain needs the 2-D lon/lat grids.

## 2. Getting OpenMRG into that format
Our local raw data (`data/OpenMRG/`) is NOT in OpenSense format: `cml.nc` has dims (time, sublink) with 728 sub-links, and link metadata sit in `cml_metadata.csv` (Frequency_GHz, Length_km, Near/Far coordinates, `Link` + `Direction` A/B).
Conversion = group the two sub-links of each `Link` into `sublink_id`, rename `Link` → `cml_id`, Near/Far coords → site_0/site_1, convert GHz → MHz and km → m. Do this in our own code (`pipeline/`).
- The notebook used `opensense_data_downloader_and_transformer.transform_andersson_2022_OpenMRG(fn=..., path_to_extract_to=..., time_start_end=(start, end), restructure_data=True)`. That module lives in the OpenSense sandbox repo on GitHub, not on pip.
- Missing values in raw `cml.nc` are `1e10`; check whether xarray masks them, otherwise set to NaN.
- `cml.nc` is 4.6 GB: open with `xr.open_dataset(path, chunks={'time': ...})` and select a time slice before `.load()`.

**Quick-start shortcut (sample data only):** `poligrain.example_data.load_openmrg(data_dir='...', subset='8d')` downloads a ready-made OpenSense-format subset (22-29 July 2015, native 10 s) and returns `ds_rad, ds_cmls, ds_gauges_municp, ds_gauge_smhi`. `subset='5min_2h'` gives a smaller 5-min version. `load_openrainer(data_dir, subset='8d')` returns `ds_rad, ds_cmls, ds_gauges`.
Caution: that sample was made by the notebook. Its CML file already contains a rain-rate `R` derived with *their* choices (radar-based wet/dry, constant baseline, Pastorek WAA). Use only `tsl`/`rsl` from it; never treat its `R` as truth.

## 3. Gauges and radar (OpenMRG)
```python
df = pd.read_csv("data/OpenMRG/gauges/city/CityGauges-2015JJA.csv", index_col=0, parse_dates=True)
meta = pd.read_csv("data/OpenMRG/gauges/city/CityGauges-metadata.csv", index_col=0)
ds_g = xr.Dataset({"rainfall_amount": (["id", "time"], df.T.values)},
                  coords={"id": meta.index.values, "time": df.index.values,
                          "lon": ("id", meta.Longitude_DecDeg.values), "lat": ("id", meta.Latitude_DecDeg.values)})
# SMHI gauge: column Pvol_mm (15-min sum), station lon 11.9924, lat 57.7156
```
Radar: `xr.open_dataset("data/OpenMRG/radar/radar.nc")` decodes pseudo-dBZ to dBZ automatically. Variable is `data` (rename to `dBZ`). Rain rate with the SMHI Z-R: `R = (10**(dBZ/10) / 200) ** (1/1.5)` mm/h.
For 2-D lon/lat: either use the `lat`/`lon` variables in the file (see readme) or project x/y with
`pyproj.Transformer.from_crs("+proj=stere +lat_ts=60 +ellps=bessel +lon_0=14 +lat_0=90", "EPSG:4326", always_xy=True)` on `np.meshgrid(x, y)`.

## 4. pycomlink (processing CML signals)
```python
import pycomlink as pycml
tl = ds.tsl - ds.rsl                                   # total loss, dB
tl = tl.interpolate_na(dim="time", method="linear", max_gap="5min")
```
- Baseline: `pycml.processing.baseline.baseline_constant(trsl, wet, n_average_last_dry=1)` (holds last dry level during wet periods); also `baseline_linear(rsl, wet, ignore_nan=False)`. NaNs in `wet` propagate to NaNs in the baseline.
- Wet/dry classifiers: `pycml.processing.wet_dry.stft.stft_classification(...)`, `...wet_dry.nearby_wetdry.nearby_wetdry(...)`, `...wet_dry.mlp.mlp_wet_dry(trsl_channel_1, trsl_channel_2, model_sel='rad_online')`.
- WAA models (useful as literature baselines for Stage 1):
  - `wet_antenna.waa_schleiss_2013(rsl, baseline, wet, waa_max, delta_t, tau)`  [ref 1]
  - `wet_antenna.waa_pastorek_2021_from_A_obs(A_obs, f_Hz, pol, L_km, A_max=14, zeta=0.55, d=0.1)`  [ref 9]
  - `wet_antenna.waa_leijnse_2008_from_A_obs(A_obs, f_Hz, pol, L_km, ...)`  [ref 4]
- Power law: `k_R_relation.calc_R_from_A(A, L_km, f_GHz=None, pol=None, a=None, b=None, a_b_approximation='ITU_2005', R_min=0.1)`. Give either `f_GHz`+`pol` or `a`+`b`. `pol` accepts 'v'/'h'/'vertical'/... Handles xarray inputs.
- Typical order in the example: `A_obs = tl - baseline` (clip at 0) → `waa` → `A = tl - baseline - waa` (clip at 0) → `R = calc_R_from_A(...)`.

## 5. poligrain (sensors in space: lines, points, grids)
```python
import poligrain as plg
# radar (or any grid) averaged along each CML path -> (cml_id, time)
gal = plg.spatial.GridAtLines(da_gridded_data=ds_rad.R, ds_line_data=ds_cmls, grid_point_location="center")
radar_along_cml = gal(ds_rad.R)
```
- Lower-level version (loop over files, memory-friendly): `calc_sparse_intersect_weights_for_several_cmls(x1_line, y1_line, x2_line, y2_line, cml_id, x_grid, y_grid, grid_point_location='center')`, then `get_grid_time_series_at_intersections(grid_data, intersect_weights)`. Dimension order must be grid (time, y, x) and weights (cml_id, y, x).
- Projection to metres: `ds.coords["x"], ds.coords["y"] = plg.spatial.project_point_coordinates(ds.lon, ds.lat, target_projection="EPSG:...")`. The example used EPSG:25832 (UTM 32N, Italy). For Gothenburg choose a Swedish/UTM 33N system (e.g. EPSG:3006 SWEREF99 TM or EPSG:32633); our choice, not from the example.
- Nearest gauge to each CML: `get_closest_points_to_line(ds_cmls, ds_gauges, max_distance, n_closest)`. Needs projected `site_0_x/_y`, `site_1_x/_y` on the CMLs and `x`, `y` on the gauges, plus `length`. Unmatched pairs come back with `neighbor_id == None`.
- Also: `GridAtPoints`, `get_closest_points_to_point`, `calc_point_to_point_distances`.
- Plots: `plg.plot_map.plot_lines(da_per_cml, vmin, vmax)`, `plg.plot_map.plot_plg(da_grid=..., da_cmls=..., da_gauges=..., vmin, vmax)` (one time step or a time-sum).
- Metrics: `plg.validation.calculate_rainfall_metrics(reference, estimate, ref_thresh=0.0, est_thresh=0.0)`, `calculate_wet_dry_metrics(...)`, `print_metrics_table`, `plot_hexbin`, `plot_confusion_matrix_count/_sum`.

## 6. xarray idioms that recur
- Select: `ds.sel(time=slice(t0, t1))`, `ds.sel(cml_id=10227)`, `ds.isel(sublink_id=0)`.
- Resample/aggregate: `.resample(time="1h").mean()` / `.sum()` / `.min()`; cumulative: `.cumsum(dim="time")`.
- Rolling statistics (10 s data: 60 steps = 10 min): `tl.rolling(time=60, center=True).std()`.
- Filter sensors: `ds.where(mask, drop=True)`.
- Align files with different sensor sets: intersect ids across files, `.sel(cml_id=common_ids)`, then `xr.concat(list, dim="time")`.
- Fill time gaps: `ds.reindex(time=pd.date_range(t0, t1, freq="15min"))`.
- Month filter: `ds.where(ds.time.dt.month.isin([4, 5, 6, 7, 8, 9, 10]), drop=True)`.
- Save compactly: `ds.to_netcdf(fn, encoding={"rsl": {"zlib": True, "complevel": 4, "dtype": "int16", "scale_factor": 0.1, "_FillValue": -999}})`.

## 7. OpenRainER file layout (from the example)
Zenodo archive (https://zenodo.org/api/records/14731404/files-archive) → `.tar` files → monthly `.gz` NetCDFs:
`CML_YYYYMM*.nc` (CMLs), `AWS_YYYYMM*.nc` (gauges, `rainfall_amount`), `RADadj_YYYYMM*.nc` (gauge-adjusted radar, 15 min, `rainfall_amount`, 1-D `lon`/`lat` → build 2-D grids with `np.meshgrid`; drop `geo_dim`, `mesh_dim`, `mosaic`).
The set of CML/gauge ids changes between months, so intersect ids before concatenating. Liquid rain only in warm months (April-October).

## 8. Pitfalls seen in the example
- Frequency MHz vs GHz and length m vs km (see section 1).
- The example's quality-control flags (rolling-std thresholds 2 dB over 5 h, 0.8 dB over 10 min) and its radar-based wet/dry (`radar_along > 0.1` mm/h) are their choices, not ours.
- Radar-based wet/dry uses the reference to process the CML, which is circular if the same radar is later used to validate.
- File paths in the example are absolute paths on other people's machines; always use paths relative to the project.
