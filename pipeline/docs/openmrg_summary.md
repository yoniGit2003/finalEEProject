# OpenMRG dataset summary

Source: Andersson, Olsson, van de Beek, Hansryd (2022), "OpenMRG: Open data from Microwave links, Radar, and Gauges for rainfall quantification in Gothenburg, Sweden", Earth Syst. Sci. Data 14, 5411-5426, https://doi.org/10.5194/essd-14-5411-2022. Local copy of the paper: `docs/papers/Andersson2022_OpenMRG.pdf`. Data: https://doi.org/10.5281/zenodo.7107689 (v1.1, CC BY-SA 4.0, about 318 MB zipped).
Collected by SMHI, Ericsson and the operator Hi3G. Period: 1 June - 31 August 2015 (JJA 2015). All times UTC.

## CML network (data/OpenMRG/cml)
- 728 sub-links = 364 bi-directional links, 363 hops, 418 nodes. Ericsson MINI-LINK radios. Star-shaped topology, densest in central Gothenburg (>= 10 sub-links per km2).
- Distance monitored 1041 km; area monitored (convex hull) about 2800 km2.
- Frequencies 7-38 GHz (most above 25 GHz), mainly vertical polarization. Paired sub-links differ by up to 1.5 GHz.
- Lengths 0.1-15 km, median about 2 km (higher frequency means shorter links).
- Antenna diameter 0.2-1.2 m (gain 31-47 dBi), mounted about 30 m above ground.
- Sampled every 10 s via SNMP (150 nodes in parallel, 7 s sampling window, nodes sampled within about +-1 s). SMHI rounded time stamps to the nearest 10 s (this changed 67 % of steps).
- TSL quantization 1 dB, RSL quantization 0.3 dB. ATPC disabled, so TSL is constant for 68 % of sub-links (occasional +-1 dB jumps from temperature).
- Values set to missing by the operator: TSL <= -99, TSL = 255, RSL <= -99.9, RSL > -20. In the NetCDF the missing value is 1e10 (set to NaN before use).
- Stats: RSL min -93.7, max -25.3, mean -44.3, SD 3.2 dBm. TSL min -10, max 26, mean 8.5, SD 6.2 dBm. Below about -90 dBm RSL is only noise. Dynamic RSL range (median - min) about 30 dB on average (0 to 54.6 dB).
- Completeness: 99.994 % of time steps have data for at least one sub-link; 93 % of sub-links have valid RSL and TSL for >= 99 % of steps; 16 sub-links (8 links) lack RSL or TSL for the whole period.
- Known signal disturbances: wet-antenna attenuation, refractivity multipath fading, construction cranes, ferries, foliage, RSL jumps, diurnal cycles, delayed return to baseline after rain.

### Files
- `cml.nc` (about 4.6 GB, CF NetCDF): variables `rsl` and `tsl` in dBm, dimensions time (10 s steps, start of sampling window, seconds since 1970-01-01 UTC) and sublink. In the example script the array order is [time, sublink].
- `cml_metadata.csv` (728 rows): Link, Direction (A/B), Sublink, NearLatitude_DecDeg, NearLongitude_DecDeg, FarLatitude_DecDeg, FarLongitude_DecDeg, Frequency_GHz, Polarization, Length_km, NearAntennaDiameter_m, FarAntennaDiameter_m. `Sublink` matches the sublink ID in cml.nc; `Link` pairs the two directions.

## Gauges (data/OpenMRG/gauges)
- SMHI station GbgA (71420, Goteborg A, central Gothenburg, 57.7156 N, 11.9924 E): Geonor weighing gauge, 0.1 mm, 15 min, `Pvol_mm` with `Pvol_QC` (G = approved, Y = uncontrolled). Interval is the 15 min before Time_UTC. Also hourly temperature, humidity, pressure, wind speed / gust / direction. Files GbgA-71420-2015JJA.csv and GbgA-71420-metadata.csv. Total JJA rainfall 256 mm, wet fraction 5.5 % (15 min), max 7.1 mm / 15 min, max daily 25.5 mm (17 June).
- City network (Goteborg Stad): 10 gauges, 1 min. 7 Geonor weighing (0.1 mm: Jarn, Torp, Bergsj, Torsl, Chalm, Tole, Barl) and 3 tipping-bucket (0.2 mm: Drakeg, Lbom, Askim). Columns: Time_UTC plus one per station; interval ends at Time_UTC (starts at second 01 of the preceding minute); NA = known error removed; tipping-bucket volumes assumed to fall in the preceding minute. Weighing gauges mean 263 mm, tipping-bucket mean 219 mm (probable undercatch). The City data had weaker quality control than SMHI.
- Gauge height about 1.5 m for the City network (assumed).

## Radar (data/OpenMRG/radar/radar.nc)
- Subset of the Swedish NORDRAD composite at the lowest elevation (0.5 deg), gauge-adjusted. 5 min, 2 x 2 km grid, 37 x 48 cells (x, y), covering about 10 km beyond the CML network. Projection: `+proj=stere +lat_ts=60 +ellps=bessel +lon_0=14 +lat_0=90`.
- Main radar Vara (78 km from Gothenburg, 5.605 GHz, C-band); Vilebo and Karlskrona are backups covering only part of the area.
- Variable `data` in pseudo-dBZ (0-255, 255 = missing); netCDF4 applies scale_factor 0.4 and add_offset -30 automatically (dBZ = 0.4*value - 30). Variables `lat`, `lon` (cell centres), `crs`, `time` (UTC seconds since 1970).
- Z-R: Z = a R^b with a = 200, b = 1.5 (SMHI summer value; Marshall-Palmer would be 1.6). Below 0 dBZ is about zero rain (about 0.03 mm/h). Clutter peak around -15 dBZ. Max 58.8 dBZ. 99.6 % complete.

## Case study (data/OpenMRG/Torslanda-event-2015-07-28)
28 July 2015 event at the Torslanda gauge (Torsl): 2 h, peak 1.1 mm/min, total 18.5 mm. Link 10227, sub-links 115 and 118 (28.23 and 29.24 GHz, vertical, 1.8 km): RSL fell about 27 dB (to -73.6 dBm) at the peak. The radar pixel above the gauge missed the peak. The RSL recovery lags the rain, which is a wet-antenna effect. The folder holds an interactive time-series zip and a radar animation.

## Processing recipe (from the paper)
Attenuation = TSL - RSL; classify wet/dry; baseline from adjacent dry steps; specific attenuation k (dB/km) beyond baseline; invert k = c R^d (ITU-R P.838-3) to get R; compensate wet-antenna attenuation. Reference software: pycomlink.
