# Data dictionary: openmrg v1

Generated 2026-10-02T17:42:28Z from catalog.json (code 211f942743ac91f118c4bbcf3e65b5484d43f7ea). Conventions:

- **time**: UTC, tz-naive; every value is labelled by the START of its interval [t, t + dt) ('interval_start'); radar stamps are instants assigned to [t, t + 5 min)
- **keys**: cml_id (int, link); sublink_id (str '<cml_id>_<1|2>'; sub-link 1 transmits site_0 -> site_1); gauge_id (City station name or 'GbgA'); event_id ('E' + start, e.g. E20150728_1455); run_id ('<event_id>_Rnn')
- **units**: dB / dBm for signal levels, mm for amounts, mm/h for rates, m for distances (EPSG:3006), GHz
- **missing**: NaN (float), fill value -32768 (int16, decoded to NaN), 0 for counts / flags
- **splits**: split_events.csv and split_links.csv are frozen; train on (train events x train links)

## Tables (tables/<name>.parquet and .csv)

### links

one row per link: radio parameters, split, status of both sub-links, completeness, TSL behaviour, nearest gauges, neighbours, radar cells. Keys: cml_id. Rows: 364.

| column | type | unit | description |
|---|---|---|---|
| cml_id | int64 |  | link id (raw metadata 'Link') |
| length_m | float64 | m | path length (m, EPSG:3006) |
| frequency_ghz_sublink_1 | float64 |  | frequency of sub-link 1 (GHz) |
| frequency_ghz_sublink_2 | float64 |  | frequency of sub-link 2 (GHz) |
| polarization | str |  | 'v' or 'h' |
| length_class | str |  | <1 km, 1-3 km, >3 km |
| frequency_band | str |  | <25 GHz or >=25 GHz (sub-link 1) |
| split_link | str |  | frozen link split: train / test / excluded |
| hop | str |  | unordered pair of site ids; parallel links share it |
| site_0_id | int64 |  | site id of end 0 |
| site_1_id | int64 |  | site id of end 1 |
| n_cells | int64 |  | radar cells crossed by the path |
| completeness_tl_min | float64 |  | worse sub-link tl completeness |
| completeness_tl_max | float64 |  | better sub-link tl completeness |
| share_not_blocked_min | float64 |  | worse sub-link not-blocked share |
| status_sublink_1 | str |  | quality status of sublink_1 |
| status_sublink_2 | str |  | quality status of sublink_2 |
| tsl_behaviour_sublink_1 | str |  | TSL behaviour of sublink_1 |
| tsl_behaviour_sublink_2 | str |  | TSL behaviour of sublink_2 |
| dry_tl_std_max | float64 | dB | larger dry tl noise of the two sub-links |
| last_valid_day | str |  | last day with valid tl |
| status_best | str |  | better status of the two sub-links |
| n_gauges_within_250m | int64 |  | gauges within 250 m of the path |
| n_gauges_within_500m | int64 |  | gauges within 500 m of the path |
| n_gauges_within_1000m | int64 |  | gauges within 1000 m of the path |
| nearest_gauge_id | str |  | nearest gauge at any distance |
| nearest_gauge_network | str |  | its network |
| nearest_gauge_distance_m | float64 | m | its distance to the path (m) |
| nearest_city_gauge_id | str |  | nearest city gauge |
| nearest_city_gauge_distance_m | float64 | m | its distance (m) |
| nearest_smhi_gauge_id | str |  | nearest smhi gauge |
| nearest_smhi_gauge_distance_m | float64 | m | its distance (m) |
| has_gauge_within_1km | bool |  | a gauge within 1 km of the path |
| n_neighbours_within_1km | int64 |  | other links whose path is within 1 km |
| n_neighbours_within_5km | int64 |  | within 5 km |
| single_radar_cell | bool |  | the whole path lies in one radar cell |

### sublinks

one row per sub-link: metadata, completeness (JJA and per month), flag shares, quantization, TSL behaviour, dry-weather noise, quality status. Keys: sublink_id. Rows: 728.

| column | type | unit | description |
|---|---|---|---|
| cml_id | int64 |  | link id (raw metadata 'Link') |
| sublink_id | str |  | sub-link id '<cml_id>_<1|2>' |
| raw_sublink_id | int64 |  | raw 'Sublink' number |
| direction | str |  | raw direction code A / B |
| tx_site | int64 |  | transmitting end (0 = site_0, 1 = site_1) |
| rx_site | int64 |  | receiving end (0 / 1) |
| frequency_ghz | float64 | GHz | carrier frequency (GHz) |
| polarization | str |  | 'v' or 'h' |
| sublink | int64 |  | sub-link number (1: site_0 -> site_1, 2: site_1 -> site_0) |
| length_m | float64 | m | path length (m, EPSG:3006) |
| site_0_id | int64 |  | site id of end 0 |
| site_1_id | int64 |  | site id of end 1 |
| site_0_antenna_diameter_m | float64 | m | antenna diameter at site 0 (m) |
| site_1_antenna_diameter_m | float64 | m | antenna diameter at site 1 (m) |
| n_radar_cells | int64 |  | radar cells crossed by the path |
| hop | str |  | unordered pair of site ids; parallel links share it |
| split_link | str |  | frozen link split: train / test / excluded |
| length_class | str |  | <1 km, 1-3 km, >3 km |
| frequency_band | str |  | <25 GHz or >=25 GHz (sub-link 1) |
| tx_site_id | int64 |  | site id of the transmitter |
| rx_site_id | int64 |  | site id of the receiver |
| completeness_rsl | float64 |  | JJA share of samples with valid rsl |
| completeness_tsl | float64 |  | JJA share of samples with valid tsl |
| completeness_tl | float64 |  | JJA share of samples with valid tl |
| share_not_blocked | float64 |  | share of samples without a blocking flag |
| flag_share_blocking_present | float64 |  | share of samples present but blocking-flagged |
| flag_share_missing_rsl | float64 |  | JJA share of samples with the missing_rsl flag |
| flag_share_missing_tsl | float64 |  | JJA share of samples with the missing_tsl flag |
| flag_share_out_of_range | float64 |  | JJA share of samples with the out_of_range flag |
| flag_share_rx_floor | float64 |  | JJA share of samples with the rx_floor flag |
| flag_share_spike | float64 |  | JJA share of samples with the spike flag |
| flag_share_tsl_step | float64 |  | JJA share of samples with the tsl_step flag |
| flag_share_flatline | float64 |  | JJA share of samples with the flatline flag |
| flag_share_network_gap | float64 |  | JJA share of samples with the network_gap flag |
| completeness_tl_jun | float64 |  | share of samples with valid tl in jun |
| completeness_tl_jul | float64 |  | share of samples with valid tl in jul |
| completeness_tl_aug | float64 |  | share of samples with valid tl in aug |
| first_valid_day | str |  | first day with valid tl |
| last_valid_day | str |  | last day with valid tl |
| rsl_min_step_db | float64 | dB | smallest RSL step (dB) |
| rsl_n_distinct | int64 | dB | distinct RSL values |
| tsl_n_distinct_values | int64 |  | distinct TSL values |
| tsl_top_share | float64 |  | share of TSL at its most common value |
| n_valid_tsl | int64 |  | valid TSL samples |
| tsl_n_steps | int64 |  | TSL changes >= 1 dB |
| tsl_step_fraction | float64 |  | share of valid samples with a TSL change |
| tsl_steps_per_day | float64 | 1/day | TSL changes per day of data |
| tsl_behaviour | str |  | no_tsl / toggling / constant / rare_steps / frequent_steps |
| tsl_steps_in_events | int64 |  | TSL steps inside events |
| tsl_steps_in_dry | int64 |  | TSL steps in dry periods |
| tsl_step_rate_events_per_h | float64 | 1/h | TSL step rate inside events (per h) |
| tsl_step_rate_dry_per_h | float64 | 1/h | rate in dry periods |
| tsl_follows_fades_descriptive | bool |  | DESCRIPTIVE ONLY: TSL steps concentrate in rain |
| dry_tl_median | float64 | dB | median dry tl (dB) |
| dry_tl_std | float64 | dB | std of the dry 5 min tl means (dB) |
| dry_tl_scaled_mad | float64 | dB | 1.4826 x MAD of the dry tl means (dB) |
| dry_tl_10s_noise_median | float64 | dB | median 5 min tl_std in dry periods (10 s noise, dB) |
| dry_rsl_std | float64 | dB | std of the dry rsl means |
| dry_rsl_10s_noise_median | float64 | dB | median 5 min rsl_std in dry periods (dB) |
| dry_tl_diurnal_range | float64 | dB | max - min of the 24 hourly dry tl medians (dB) |
| dry_tl_daily_drift_std | float64 | dB | std of the daily dry tl medians (dB) |
| n_dry_days | int64 |  | days with >= 1 h of valid dry tl |
| dry_baseline_range_db | float64 | dB | p90 - p10 of the dry-period medians of tl (dB) |
| n_dry_periods_used | int64 |  | dry periods with >= 50 % valid bins |
| status | str |  | quality status good / degraded / unusable (step 7 rules) |
| status_reason | str |  | rules that fired |
| status_step3 | str |  | status of the step 3 QC (kept for reference) |

### gauges

City (10) and SMHI (1) gauges: position, type, completeness, JJA totals, radar pixel comparison. Keys: gauge_id. Rows: 11.

| column | type | unit | description |
|---|---|---|---|
| gauge_id | str |  | gauge id (City station name or 'GbgA') |
| network | str |  | 'city' or 'smhi' |
| gauge_type | str |  | weighing or tipping_bucket |
| lon | float64 | degrees_east | longitude (deg, WGS84) |
| lat | float64 | degrees_north | latitude (deg, WGS84) |
| x | float64 | m | x (m, EPSG:3006) |
| y | float64 | m | y (m, EPSG:3006) |
| quantization_mm | float64 | mm | gauge resolution (mm) |
| height_m | float64 | m | gauge height above ground (m) |
| interval_s | int64 |  | native interval (s) |
| n_valid | int64 |  | valid native steps |
| completeness_pct | float64 | % | valid share of native steps (%) |
| total_mm | float64 | mm | L1 total over JJA (mm) |
| n_negative | int64 |  | negative amounts in the raw file |
| jja_total_mm | float64 | mm | JJA total (mm) |
| wet_fraction | float64 |  | share of native steps with at least one quantum |
| completeness_jun | float64 |  | share of valid native steps in jun |
| completeness_jul | float64 |  | share of valid native steps in jul |
| completeness_aug | float64 |  | share of valid native steps in aug |
| iy | int64 |  | radar row index |
| ix | int64 |  | radar column index |
| distance_to_cell_centre_m | float64 | m | gauge distance to its radar cell centre (m) |
| radar_pixel_jja_total_mm | float64 | mm | JJA radar total at the gauge pixel (mm) |
| radar_over_gauge_ratio | float64 |  | radar / gauge JJA total |

### sites

physical sites (masts): link ends merged within 10 m. Keys: site_id. Rows: 418.

| column | type | unit | description |
|---|---|---|---|
| site_id | int64 |  | site (mast) id |
| lon | float64 | degrees_east | longitude (deg, WGS84) |
| lat | float64 | degrees_north | latitude (deg, WGS84) |
| x_3006_m | float64 | m | x (m, EPSG:3006) |
| y_3006_m | float64 | m | y (m, EPSG:3006) |
| n_links | int64 |  | links with an end at the site |
| n_ends | int64 |  | link ends merged into the site |
| max_member_offset_m | float64 | m | largest offset of a merged end from the site (m) |

### events

rain events defined from radar + City gauges (never CML). Keys: event_id. Rows: 89. **Frozen.**

| column | type | unit | description |
|---|---|---|---|
| event_id | str |  | event id 'E<start>' |
| start | datetime64[us] |  | start (UTC, interval start) |
| end | datetime64[us] |  | end (UTC, exclusive) |
| duration_h | float64 | h | duration (h) |
| gap_before_h | float64 | h | non-wet time back to the previous event, capped at 48 h |
| dry_before_h | float64 | h | certified-dry hours in that gap |
| gap_after_h | float64 | h | non-wet time to the next event, capped at 48 h |
| dry_after_h | float64 | h | certified-dry hours in that gap |
| merged_from | int64 |  | wet runs merged into the event |
| touches_period_edge | bool |  | less than one min_dry_gap from the period start / end |
| radar_peak_domain_rain_mmh | float64 | mm/h | largest domain-mean radar rate (mm/h) |
| radar_mean_domain_rain_mmh | float64 | mm/h | mean domain-mean radar rate (mm/h) |
| radar_max_wet_area_fraction | float64 |  | largest wet share of the domain cells |
| radar_domain_total_mm | float64 | mm | domain-mean radar accumulation (mm averaged over the 720 domain cells) |
| city_total_Jarn_mm | float64 | mm | total of City gauge Jarn in the event (mm) |
| city_total_Torp_mm | float64 | mm | total of City gauge Torp in the event (mm) |
| city_total_Bergsj_mm | float64 | mm | total of City gauge Bergsj in the event (mm) |
| city_total_Torsl_mm | float64 | mm | total of City gauge Torsl in the event (mm) |
| city_total_Chalm_mm | float64 | mm | total of City gauge Chalm in the event (mm) |
| city_total_Tole_mm | float64 | mm | total of City gauge Tole in the event (mm) |
| city_total_Barl_mm | float64 | mm | total of City gauge Barl in the event (mm) |
| city_total_Drakeg_mm | float64 | mm | total of City gauge Drakeg in the event (mm) |
| city_total_Lbom_mm | float64 | mm | total of City gauge Lbom in the event (mm) |
| city_total_Askim_mm | float64 | mm | total of City gauge Askim in the event (mm) |
| city_total_mean_mm | float64 | mm | mean of the City gauge totals (mm) |
| city_total_max_mm | float64 | mm | largest City gauge total (mm) |
| city_max_5min_rate_mmh | float64 | mm/h | largest 5 min City gauge rate (mm/h) |
| city_max_5min_gauge | str |  | gauge of that maximum |
| city_max_1min_rate_mmh | float64 | mm/h | largest 1 min City gauge rate (mm/h) |
| city_max_1min_gauge | str |  | gauge of that maximum |
| smhi_total_mm_info | float64 |  | SMHI GbgA total over the overlapping 15 min bins (information only) |
| n_unknown_bins | int64 |  | 5 min bins with state 'unknown' |
| radar_missing_fraction | float64 |  | share of bins without valid radar |
| network_gap_minutes | float64 | min | minutes of network gaps in the event |
| temperature_min | float64 | degC | min of the hourly SMHI temperature around the event |
| temperature_max | float64 | degC | max of the hourly SMHI temperature around the event |
| temperature_mean | float64 | degC | mean of the hourly SMHI temperature around the event |
| humidity_min | float64 | % | min of the hourly SMHI humidity around the event |
| humidity_max | float64 | % | max of the hourly SMHI humidity around the event |
| humidity_mean | float64 | % | mean of the hourly SMHI humidity around the event |
| wind_speed_min | float64 | m s-1 | min of the hourly SMHI wind_speed around the event |
| wind_speed_max | float64 | m s-1 | max of the hourly SMHI wind_speed around the event |
| wind_speed_mean | float64 | m s-1 | mean of the hourly SMHI wind_speed around the event |
| wind_gust_min | float64 | m s-1 | min of the hourly SMHI wind_gust around the event |
| wind_gust_max | float64 | m s-1 | max of the hourly SMHI wind_gust around the event |
| wind_gust_mean | float64 | m s-1 | mean of the hourly SMHI wind_gust around the event |
| n_weather_hours | int64 |  | hourly weather values used |
| month | int64 |  | month of the start |
| intensity_class | str |  | tercile of the max 5 min City rate over events |
| split_event | str |  | frozen event split: train / val / test |

### event_significance

significance, rain class and wet core of every event. Keys: event_id. Rows: 89.

| column | type | unit | description |
|---|---|---|---|
| event_id | str |  | event id 'E<start>' |
| significant | bool |  | event (or run) passes the significance rule |
| sig_reason | str |  | which rule made it significant |
| rain_class | str |  | none_drizzle / significant / heavy |
| city_total_mean_mm | float64 | mm | mean of the City gauge totals (mm) |
| radar_domain_total_mm | float64 | mm | domain-mean radar accumulation (mm averaged over the 720 domain cells) |
| city_max_5min_rate_mmh | float64 | mm/h | largest 5 min City gauge rate (mm/h) |
| n_runs | int64 |  | wet runs in the event |
| wet_core_start | datetime64[us] |  | start of the first significant run |
| wet_core_end | datetime64[us] |  | end of the last significant run |
| wet_core_hours | float64 | h | wet core duration (h) |
| split_event | str |  | frozen event split: train / val / test |

### event_runs

wet runs inside events (before gap merging). Keys: run_id. Rows: 368.

| column | type | unit | description |
|---|---|---|---|
| run_id | str |  | wet-run id '<event_id>_Rnn' |
| event_id | str |  | event id 'E<start>' |
| start | datetime64[us] |  | start (UTC, interval start) |
| end | datetime64[us] |  | end (UTC, exclusive) |
| duration_min | float64 | min | duration (min) |
| radar_peak_domain_rain_mmh | float64 | mm/h | largest domain-mean radar rate (mm/h) |
| radar_max_wet_area_fraction | float64 |  | largest wet share of the domain cells |
| radar_domain_total_mm | float64 | mm | domain-mean radar accumulation (mm averaged over the 720 domain cells) |
| city_total_mean_mm | float64 | mm | mean of the City gauge totals (mm) |
| city_max_5min_rate_mmh | float64 | mm/h | largest 5 min City gauge rate (mm/h) |
| gap_to_previous_run_h | float64 | h | time since the previous run of the same event (h) |
| significant_run | bool |  | run is significant |

### dry_periods

certified dry periods (>= 6 h). Keys: start. Rows: 62. **Frozen.**

| column | type | unit | description |
|---|---|---|---|
| start | datetime64[us] |  | start (UTC, interval start) |
| end | datetime64[us] |  | end (UTC, exclusive) |
| duration_h | float64 | h | duration (h) |
| radar_valid_fraction | float64 |  | mean share of domain cells with radar data |

### split_events

FROZEN event split (train / val / test). Keys: event_id. Rows: 89. **Frozen.**

| column | type | unit | description |
|---|---|---|---|
| event_id | str |  | event id 'E<start>' |
| intensity_class | str |  | tercile of the max 5 min City rate over events |
| split_event | str |  | frozen event split: train / val / test |

### split_links

FROZEN link split (train / test / excluded). Keys: cml_id. Rows: 364. **Frozen.**

| column | type | unit | description |
|---|---|---|---|
| cml_id | int64 |  | link id (raw metadata 'Link') |
| hop | str |  | unordered pair of site ids; parallel links share it |
| length_class | str |  | <1 km, 1-3 km, >3 km |
| frequency_band | str |  | <25 GHz or >=25 GHz (sub-link 1) |
| stratum | str |  | length class x frequency band |
| excluded | bool |  | link excluded from the split (both sub-links unusable at step 3) |
| split_link | str |  | frozen link split: train / test / excluded |
| length_m | float64 | m | path length (m, EPSG:3006) |
| frequency_ghz_sublink_1 | float64 |  | frequency of sub-link 1 (GHz) |
| site_0_id | int64 |  | site id of end 0 |
| site_1_id | int64 |  | site id of end 1 |

### event_links

one row per (event, link): radar along the path, matched gauge, CML data quality. Keys: event_id, cml_id. Rows: 32396. **Frozen.**

| column | type | unit | description |
|---|---|---|---|
| event_id | str |  | event id 'E<start>' |
| cml_id | int64 |  | link id (raw metadata 'Link') |
| length_m | float64 | m | path length (m, EPSG:3006) |
| frequency_ghz_sublink_1 | float64 |  | frequency of sub-link 1 (GHz) |
| frequency_ghz_sublink_2 | float64 |  | frequency of sub-link 2 (GHz) |
| polarization | str |  | 'v' or 'h' |
| n_cells | int64 |  | radar cells crossed by the path |
| radar_total_mm | float64 | mm | radar along the path, total in the event (mm) |
| radar_peak_mmh | float64 | mm/h | largest 5 min R_mean in the event (mm/h) |
| radar_wet_minutes | float64 | min | minutes with R_mean >= 0.1 mm/h |
| radar_mean_wet_fraction | float64 |  | mean wet share of the path |
| radar_max_R_max_mmh | float64 | mm/h | largest R_max (mm/h) |
| radar_mean_valid_weight_fraction | float64 |  | mean share of the path with radar data |
| gauge_id | str |  | gauge id (City station name or 'GbgA') |
| gauge_network | str |  | network of the matched gauge |
| gauge_distance_to_path_m | float64 | m | matched gauge distance to the path (m) |
| gauge_position_along_path | float64 |  | matched gauge position (0..1) |
| gauge_total_mm | float64 | mm | matched gauge total in the event (mm) |
| gauge_max_rate_mmh | float64 | mm/h | matched gauge max rate (mm/h) |
| gauge_max_rate_interval | str |  | interval of that rate (5min City, 15min SMHI) |
| n_gauges_within_1km | int64 |  | gauges within 1 km of the path |
| valid_fraction_event_sublink_1 | float64 |  | valid share of tl_mean bins event (sublink_1) |
| valid_fraction_event_sublink_2 | float64 |  | valid share of tl_mean bins event (sublink_2) |
| valid_fraction_event_min | float64 |  | valid share of tl_mean bins event (min) |
| valid_fraction_before_sublink_1 | float64 |  | valid share of tl_mean bins before (sublink_1) |
| valid_fraction_before_sublink_2 | float64 |  | valid share of tl_mean bins before (sublink_2) |
| valid_fraction_before_min | float64 |  | valid share of tl_mean bins before (min) |
| valid_fraction_after_sublink_1 | float64 |  | valid share of tl_mean bins after (sublink_1) |
| valid_fraction_after_sublink_2 | float64 |  | valid share of tl_mean bins after (sublink_2) |
| valid_fraction_after_min | float64 |  | valid share of tl_mean bins after (min) |
| frac_rx_floor_event_sublink_1 | float64 |  | share of 5 min bins in the event whose flags_or holds rx_floor (sublink_1) |
| frac_rx_floor_event_sublink_2 | float64 |  | share of 5 min bins in the event whose flags_or holds rx_floor (sublink_2) |
| frac_rx_floor_event_max | float64 |  | share of 5 min bins in the event whose flags_or holds rx_floor (max) |
| frac_spike_event_sublink_1 | float64 |  | share of 5 min bins in the event whose flags_or holds spike (sublink_1) |
| frac_spike_event_sublink_2 | float64 |  | share of 5 min bins in the event whose flags_or holds spike (sublink_2) |
| frac_spike_event_max | float64 |  | share of 5 min bins in the event whose flags_or holds spike (max) |
| frac_network_gap_event_sublink_1 | float64 |  | share of 5 min bins in the event whose flags_or holds network_gap (sublink_1) |
| frac_network_gap_event_sublink_2 | float64 |  | share of 5 min bins in the event whose flags_or holds network_gap (sublink_2) |
| frac_network_gap_event_max | float64 |  | share of 5 min bins in the event whose flags_or holds network_gap (max) |
| status_sublink_1 | str |  | quality status of sublink_1 |
| status_sublink_2 | str |  | quality status of sublink_2 |
| cml_ok_event | bool |  | a sub-link has >= 80 % valid tl_mean bins in the event |
| cml_ok_before | bool |  | ... in the 3 h before |
| cml_ok_after | bool |  | ... in the 3 h after |
| has_reference_gauge | bool |  | a gauge within 1 km of the path |
| split_event | str |  | frozen event split: train / val / test |
| split_link | str |  | frozen link split: train / test / excluded |

### event_links_flags_extra

per (event, link): flatline / tsl_step shares and counts in the event (information only). Keys: event_id, cml_id. Rows: 32396.

| column | type | unit | description |
|---|---|---|---|
| event_id | str |  | event id 'E<start>' |
| cml_id | int64 |  | link id (raw metadata 'Link') |
| frac_flatline_event_sublink_1 | float64 |  | share of 5 min bins in the event whose flags_or holds flatline (sublink_1) |
| frac_flatline_event_sublink_2 | float64 |  | share of 5 min bins in the event whose flags_or holds flatline (sublink_2) |
| frac_flatline_event_max | float64 |  | share of 5 min bins in the event whose flags_or holds flatline (max) |
| frac_tsl_step_event_sublink_1 | float64 |  | share of 5 min bins in the event whose flags_or holds tsl_step (sublink_1) |
| frac_tsl_step_event_sublink_2 | float64 |  | share of 5 min bins in the event whose flags_or holds tsl_step (sublink_2) |
| frac_tsl_step_event_max | float64 |  | share of 5 min bins in the event whose flags_or holds tsl_step (max) |
| n_tsl_step_in_event_sublink_1 | int64 |  | 10 s TSL steps in the event (sublink_1) |
| n_tsl_step_in_event_sublink_2 | int64 |  | 10 s TSL steps in the event (sublink_2) |
| n_tsl_step_in_event | int64 |  | 10 s TSL steps in the event (both sub-links) |
| split_event | str |  | frozen event split: train / val / test |
| split_link | str |  | frozen link split: train / test / excluded |

### link_gauge_matches

every (link, gauge) pair within 5 km of the path. Keys: cml_id, gauge_id. Rows: 852.

| column | type | unit | description |
|---|---|---|---|
| cml_id | int64 |  | link id (raw metadata 'Link') |
| gauge_id | str |  | gauge id (City station name or 'GbgA') |
| network | str |  | 'city' or 'smhi' |
| gauge_type | str |  | weighing or tipping_bucket |
| quantization_mm | float64 | mm | gauge resolution (mm) |
| distance_to_path_m | float64 | m | gauge distance to the path segment (m) |
| position_along_path | float64 |  | 0 = site_0 end, 1 = site_1 end |
| distance_to_midpoint_m | float64 | m | gauge distance to the path midpoint (m) |
| distance_to_nearest_end_m | float64 | m | to the nearest end (m) |
| within_500m | bool |  | gauge within 500 m |
| within_1000m | bool |  | gauge within 1 km |
| within_2000m | bool |  | gauge within 2 km |
| link_length_m | float64 | m | path length (m) |
| frequency_ghz | float64 | GHz | carrier frequency (GHz) |

### network_gaps

intervals when many alive sub-links were missing (collector outages). Keys: start. Rows: 93.

| column | type | unit | description |
|---|---|---|---|
| start | datetime64[us] |  | start (UTC, interval start) |
| end | datetime64[us] |  | end (UTC, exclusive) |
| duration_min | float64 | min | duration (min) |
| max_fraction_nan | float64 |  | largest NaN share of alive sub-links |
| n_alive | int64 |  | alive sub-links that day |

### quality_status

quality status per sub-link (good / degraded / unusable) and the reason. Keys: sublink_id. Rows: 728.

| column | type | unit | description |
|---|---|---|---|
| cml_id | int64 |  | link id (raw metadata 'Link') |
| sublink_id | str |  | sub-link id '<cml_id>_<1|2>' |
| sublink | int8 |  | sub-link number (1: site_0 -> site_1, 2: site_1 -> site_0) |
| raw_sublink_id | int64 |  | raw 'Sublink' number |
| status | str |  | quality status good / degraded / unusable (step 7 rules) |
| status_reason | str |  | rules that fired |
| status_step3 | str |  | status of the step 3 QC (kept for reference) |
| completeness_tl | float64 |  | JJA share of samples with valid tl |
| completeness_tl_jun | float64 |  | share of samples with valid tl in jun |
| completeness_tl_jul | float64 |  | share of samples with valid tl in jul |
| completeness_tl_aug | float64 |  | share of samples with valid tl in aug |
| share_not_blocked | float64 |  | share of samples without a blocking flag |
| flag_share_blocking_present | float64 |  | share of samples present but blocking-flagged |
| flag_share_rx_floor | float64 |  | JJA share of samples with the rx_floor flag |
| flag_share_spike | float64 |  | JJA share of samples with the spike flag |
| flag_share_flatline | float64 |  | JJA share of samples with the flatline flag |
| flag_share_tsl_step | float64 |  | JJA share of samples with the tsl_step flag |
| dry_tl_std | float64 | dB | std of the dry 5 min tl means (dB) |
| dry_baseline_range_db | float64 | dB | p90 - p10 of the dry-period medians of tl (dB) |

### sublink_dry_stats

descriptive CML statistics per sub-link over the certified dry periods. Keys: sublink_id. Rows: 728.

| column | type | unit | description |
|---|---|---|---|
| cml_id | int64 |  | link id (raw metadata 'Link') |
| sublink_id | str |  | sub-link id '<cml_id>_<1|2>' |
| sublink | int8 |  | sub-link number (1: site_0 -> site_1, 2: site_1 -> site_0) |
| n_dry_bins | int64 |  | 5 min bins in dry periods |
| n_valid_tl | int64 |  | valid tl_mean bins |
| valid_fraction | float64 |  | valid share of the dry bins |
| tl_median | float64 | dB | median of the dry 5 min tl means (dB) |
| tl_mean | float64 | dB | mean of the dry 5 min tl means (dB) |
| tl_std | float64 | dB | std of the dry 5 min tl means (dB) |
| tl_scaled_mad | float64 | dB | scaled_mad of the dry 5 min tl means (dB) |
| tl_iqr | float64 | dB | iqr of the dry 5 min tl means (dB) |
| tl_p1 | float64 | dB | p1 of the dry 5 min tl means (dB) |
| tl_p99 | float64 | dB | p99 of the dry 5 min tl means (dB) |
| tl_intra_bin_std_median | float64 | dB | median 5 min tl_std in dry periods (dB) |
| rsl_median | float64 | dB | median of the dry 5 min rsl means (dB) |
| rsl_std | float64 | dB | std of the dry 5 min rsl means (dB) |
| rsl_scaled_mad | float64 | dB | scaled_mad of the dry 5 min rsl means (dB) |
| rsl_intra_bin_std_median | float64 | dB | median 5 min rsl_std in dry periods (dB) |
| n_dry_days | int64 |  | days with >= 1 h of valid dry tl |
| tl_diurnal_range | float64 | dB | max - min of the 24 hourly dry tl medians (dB) |
| tl_daily_median_drift_std | float64 | dB | std of the daily dry tl medians (dB) |

## Time series (timeseries/<name>.zarr)

### cml_10s

CML at 10 s: tsl, rsl, tl (L1, int16 x 0.1 dB, exact copies) and the L2 QC flag bitmask. Resolution 10s, dims time=794880, sublink_id=728, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks [60480, 64].

| variable | dims | type | unit | description |
|---|---|---|---|---|
| flags | time, sublink_id | uint16 |  | L2 QC flags (bitmask) |
| rsl | time, sublink_id | float64 | dBm | received signal level |
| tl | time, sublink_id | float64 | dB | total loss, tsl - rsl |
| tsl | time, sublink_id | float64 | dBm | transmitted signal level |

### cml_1min

CML statistics per 1 min bin (L3). Resolution 1min, dims time=132480, sublink_id=728, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks [132480, 32].

| variable | dims | type | unit | description |
|---|---|---|---|---|
| flags_or | time, sublink_id | uint16 |  | bitwise OR of all L2 flag bits of the samples in the bin |
| n_tsl_step | time, sublink_id | uint8 | 1 | number of samples flagged tsl_step in the bin |
| n_valid_rsl | time, sublink_id | uint8 | 1 | number of valid rsl samples in the bin |
| n_valid_tl | time, sublink_id | uint8 | 1 | number of valid tl samples in the bin |
| n_valid_tsl | time, sublink_id | uint8 | 1 | number of valid tsl samples in the bin |
| rsl_max | time, sublink_id | float32 | dBm | maximum of valid rsl samples in the bin |
| rsl_mean | time, sublink_id | float32 | dBm | mean of valid rsl samples in the bin |
| rsl_median | time, sublink_id | float32 | dBm | median of valid rsl samples in the bin |
| rsl_min | time, sublink_id | float32 | dBm | minimum of valid rsl samples in the bin |
| tl_max | time, sublink_id | float32 | dB | maximum of valid tl samples in the bin |
| tl_mean | time, sublink_id | float32 | dB | mean of valid tl samples in the bin |
| tl_median | time, sublink_id | float32 | dB | median of valid tl samples in the bin |
| tl_min | time, sublink_id | float32 | dB | minimum of valid tl samples in the bin |
| tl_std | time, sublink_id | float32 | dB | standard deviation (ddof=0) of valid tl samples in the bin |
| tsl_max | time, sublink_id | float32 | dBm | maximum of valid tsl samples in the bin |
| tsl_mean | time, sublink_id | float32 | dBm | mean of valid tsl samples in the bin |
| tsl_min | time, sublink_id | float32 | dBm | minimum of valid tsl samples in the bin |

### cml_5min

CML statistics per 5 min bin (L3). Resolution 5min, dims time=26496, sublink_id=728, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks [8928, 64].

| variable | dims | type | unit | description |
|---|---|---|---|---|
| flags_or | time, sublink_id | uint16 |  | bitwise OR of all L2 flag bits of the samples in the bin |
| n_tsl_step | time, sublink_id | uint8 | 1 | number of samples flagged tsl_step in the bin |
| n_valid_rsl | time, sublink_id | uint8 | 1 | number of valid rsl samples in the bin |
| n_valid_tl | time, sublink_id | uint8 | 1 | number of valid tl samples in the bin |
| n_valid_tsl | time, sublink_id | uint8 | 1 | number of valid tsl samples in the bin |
| rsl_max | time, sublink_id | float32 | dBm | maximum of valid rsl samples in the bin |
| rsl_mean | time, sublink_id | float32 | dBm | mean of valid rsl samples in the bin |
| rsl_median | time, sublink_id | float32 | dBm | median of valid rsl samples in the bin |
| rsl_min | time, sublink_id | float32 | dBm | minimum of valid rsl samples in the bin |
| rsl_std | time, sublink_id | float32 | dBm | standard deviation (ddof=0) of valid rsl samples in the bin |
| tl_max | time, sublink_id | float32 | dB | maximum of valid tl samples in the bin |
| tl_mean | time, sublink_id | float32 | dB | mean of valid tl samples in the bin |
| tl_median | time, sublink_id | float32 | dB | median of valid tl samples in the bin |
| tl_min | time, sublink_id | float32 | dB | minimum of valid tl samples in the bin |
| tl_std | time, sublink_id | float32 | dB | standard deviation (ddof=0) of valid tl samples in the bin |
| tsl_max | time, sublink_id | float32 | dBm | maximum of valid tsl samples in the bin |
| tsl_mean | time, sublink_id | float32 | dBm | mean of valid tsl samples in the bin |
| tsl_min | time, sublink_id | float32 | dBm | minimum of valid tsl samples in the bin |

### cml_15min

CML statistics per 15 min bin (L3). Resolution 15min, dims time=8832, sublink_id=728, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks [8832, 64].

| variable | dims | type | unit | description |
|---|---|---|---|---|
| flags_or | time, sublink_id | uint16 |  | bitwise OR of all L2 flag bits of the samples in the bin |
| n_tsl_step | time, sublink_id | uint8 | 1 | number of samples flagged tsl_step in the bin |
| n_valid_rsl | time, sublink_id | uint8 | 1 | number of valid rsl samples in the bin |
| n_valid_tl | time, sublink_id | uint8 | 1 | number of valid tl samples in the bin |
| n_valid_tsl | time, sublink_id | uint8 | 1 | number of valid tsl samples in the bin |
| rsl_max | time, sublink_id | float32 | dBm | maximum of valid rsl samples in the bin |
| rsl_mean | time, sublink_id | float32 | dBm | mean of valid rsl samples in the bin |
| rsl_median | time, sublink_id | float32 | dBm | median of valid rsl samples in the bin |
| rsl_min | time, sublink_id | float32 | dBm | minimum of valid rsl samples in the bin |
| rsl_std | time, sublink_id | float32 | dBm | standard deviation (ddof=0) of valid rsl samples in the bin |
| tl_max | time, sublink_id | float32 | dB | maximum of valid tl samples in the bin |
| tl_mean | time, sublink_id | float32 | dB | mean of valid tl samples in the bin |
| tl_median | time, sublink_id | float32 | dB | median of valid tl samples in the bin |
| tl_min | time, sublink_id | float32 | dB | minimum of valid tl samples in the bin |
| tl_std | time, sublink_id | float32 | dB | standard deviation (ddof=0) of valid tl samples in the bin |
| tsl_max | time, sublink_id | float32 | dBm | maximum of valid tsl samples in the bin |
| tsl_mean | time, sublink_id | float32 | dBm | mean of valid tsl samples in the bin |
| tsl_min | time, sublink_id | float32 | dBm | minimum of valid tsl samples in the bin |

### cml_1h

CML statistics per 1 h bin (L3). Resolution 1h, dims time=2208, sublink_id=728, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks [2208, 728].

| variable | dims | type | unit | description |
|---|---|---|---|---|
| flags_or | time, sublink_id | uint16 |  | bitwise OR of all L2 flag bits of the samples in the bin |
| n_tsl_step | time, sublink_id | uint16 | 1 | number of samples flagged tsl_step in the bin |
| n_valid_rsl | time, sublink_id | uint16 | 1 | number of valid rsl samples in the bin |
| n_valid_tl | time, sublink_id | uint16 | 1 | number of valid tl samples in the bin |
| n_valid_tsl | time, sublink_id | uint16 | 1 | number of valid tsl samples in the bin |
| rsl_max | time, sublink_id | float32 | dBm | maximum of valid rsl samples in the bin |
| rsl_mean | time, sublink_id | float32 | dBm | mean of valid rsl samples in the bin |
| rsl_median | time, sublink_id | float32 | dBm | median of valid rsl samples in the bin |
| rsl_min | time, sublink_id | float32 | dBm | minimum of valid rsl samples in the bin |
| rsl_std | time, sublink_id | float32 | dBm | standard deviation (ddof=0) of valid rsl samples in the bin |
| tl_max | time, sublink_id | float32 | dB | maximum of valid tl samples in the bin |
| tl_mean | time, sublink_id | float32 | dB | mean of valid tl samples in the bin |
| tl_median | time, sublink_id | float32 | dB | median of valid tl samples in the bin |
| tl_min | time, sublink_id | float32 | dB | minimum of valid tl samples in the bin |
| tl_std | time, sublink_id | float32 | dB | standard deviation (ddof=0) of valid tl samples in the bin |
| tsl_max | time, sublink_id | float32 | dBm | maximum of valid tsl samples in the bin |
| tsl_mean | time, sublink_id | float32 | dBm | mean of valid tsl samples in the bin |
| tsl_min | time, sublink_id | float32 | dBm | minimum of valid tsl samples in the bin |

### radar_along_links_5min

radar rain rate along every link path (floored at 0 dBZ), 5 min. Resolution 5min, dims time=26496, cml_id=364, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks [8928, 64].

| variable | dims | type | unit | description |
|---|---|---|---|---|
| n_valid_stamps | time, cml_id | uint8 | 1 | stamps with valid_weight_fraction >= 0.5 in the bin |
| R_max | time, cml_id | float32 | mm h-1 | largest cell rate under the path, max over the valid stamps |
| R_mean | time, cml_id | float32 | mm h-1 | path-weighted mean rain rate (floored), mean over the valid stamps of the bin |
| R_std | time, cml_id | float32 | mm h-1 | path-weighted std of the cell rates along the path, mean over the valid stamps |
| valid_weight_fraction | time, cml_id | float32 | 1 | share of the path weight with radar data, mean over ALL stamps of the bin (a missing stamp counts as 0) |
| wet_fraction | time, cml_id | float32 | 1 | share of the valid path weight in cells > 0.1 mm/h, mean over valid stamps |

### radar_along_links_15min

radar rain rate along every link path, 15 min. Resolution 15min, dims time=8832, cml_id=364, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks [8832, 64].

| variable | dims | type | unit | description |
|---|---|---|---|---|
| n_valid_stamps | time, cml_id | uint8 | 1 | stamps with valid_weight_fraction >= 0.5 in the bin |
| R_max | time, cml_id | float32 | mm h-1 | largest cell rate under the path, max over the valid stamps |
| R_mean | time, cml_id | float32 | mm h-1 | path-weighted mean rain rate (floored), mean over the valid stamps of the bin |
| R_std | time, cml_id | float32 | mm h-1 | path-weighted std of the cell rates along the path, mean over the valid stamps |
| valid_weight_fraction | time, cml_id | float32 | 1 | share of the path weight with radar data, mean over ALL stamps of the bin (a missing stamp counts as 0) |
| wet_fraction | time, cml_id | float32 | 1 | share of the valid path weight in cells > 0.1 mm/h, mean over valid stamps |

### radar_along_links_1h

radar rain rate along every link path, 1 h. Resolution 1h, dims time=2208, cml_id=364, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks [2208, 364].

| variable | dims | type | unit | description |
|---|---|---|---|---|
| n_valid_stamps | time, cml_id | uint8 | 1 | stamps with valid_weight_fraction >= 0.5 in the bin |
| R_max | time, cml_id | float32 | mm h-1 | largest cell rate under the path, max over the valid stamps |
| R_mean | time, cml_id | float32 | mm h-1 | path-weighted mean rain rate (floored), mean over the valid stamps of the bin |
| R_std | time, cml_id | float32 | mm h-1 | path-weighted std of the cell rates along the path, mean over the valid stamps |
| valid_weight_fraction | time, cml_id | float32 | 1 | share of the path weight with radar data, mean over ALL stamps of the bin (a missing stamp counts as 0) |
| wet_fraction | time, cml_id | float32 | 1 | share of the valid path weight in cells > 0.1 mm/h, mean over valid stamps |

### radar_at_gauges_5min

radar rain rate at every gauge (pixel and 3 x 3 mean), 5 min. Resolution 5min, dims time=26496, gauge_id=11, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks one chunk per variable.

| variable | dims | type | unit | description |
|---|---|---|---|---|
| n_valid_stamps_3x3 | time, gauge_id | uint8 | 1 | valid stamps of R_3x3 in the bin |
| n_valid_stamps_pixel | time, gauge_id | uint8 | 1 | valid stamps of R_pixel in the bin |
| R_3x3 | time, gauge_id | float32 | mm h-1 | mean rain rate (floored) of the 3 x 3 cells around the gauge cell (NaN cells excluded), mean over valid stamps |
| R_pixel | time, gauge_id | float32 | mm h-1 | rain rate (floored) of the radar cell containing the gauge, mean over valid stamps |

### radar_at_gauges_15min

radar rain rate at every gauge, 15 min. Resolution 15min, dims time=8832, gauge_id=11, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks one chunk per variable.

| variable | dims | type | unit | description |
|---|---|---|---|---|
| n_valid_stamps_3x3 | time, gauge_id | uint8 | 1 | valid stamps of R_3x3 in the bin |
| n_valid_stamps_pixel | time, gauge_id | uint8 | 1 | valid stamps of R_pixel in the bin |
| R_3x3 | time, gauge_id | float32 | mm h-1 | mean rain rate (floored) of the 3 x 3 cells around the gauge cell (NaN cells excluded), mean over valid stamps |
| R_pixel | time, gauge_id | float32 | mm h-1 | rain rate (floored) of the radar cell containing the gauge, mean over valid stamps |

### radar_at_gauges_1h

radar rain rate at every gauge, 1 h. Resolution 1h, dims time=2208, gauge_id=11, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks one chunk per variable.

| variable | dims | type | unit | description |
|---|---|---|---|---|
| n_valid_stamps_3x3 | time, gauge_id | uint8 | 1 | valid stamps of R_3x3 in the bin |
| n_valid_stamps_pixel | time, gauge_id | uint8 | 1 | valid stamps of R_pixel in the bin |
| R_3x3 | time, gauge_id | float32 | mm h-1 | mean rain rate (floored) of the 3 x 3 cells around the gauge cell (NaN cells excluded), mean over valid stamps |
| R_pixel | time, gauge_id | float32 | mm h-1 | rain rate (floored) of the radar cell containing the gauge, mean over valid stamps |

### gauges_city_1min

City gauges, 1 min amounts and rates. Resolution 1min, dims time=132480, gauge_id=10, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks one chunk per variable.

| variable | dims | type | unit | description |
|---|---|---|---|---|
| amount_mm | time, gauge_id | float32 | mm | rainfall amount summed over the bin |
| n_valid | time, gauge_id | uint8 | 1 | valid native steps in the bin |
| rate_mm_h | time, gauge_id | float32 | mm h-1 | mean rain rate over the bin, amount_mm / dt_hours |

### gauges_city_5min

City gauges, 5 min sums (strict: NaN with a missing minute). Resolution 5min, dims time=26496, gauge_id=10, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks one chunk per variable.

| variable | dims | type | unit | description |
|---|---|---|---|---|
| amount_mm | time, gauge_id | float32 | mm | rainfall amount summed over the bin |
| n_valid | time, gauge_id | uint8 | 1 | valid native steps in the bin |
| peak_1min_rate_mm_h | time, gauge_id | float32 | mm h-1 | largest native 1 min amount * 60 in the bin |
| rate_mm_h | time, gauge_id | float32 | mm h-1 | mean rain rate over the bin, amount_mm / dt_hours |

### gauges_city_15min

City gauges, 15 min sums. Resolution 15min, dims time=8832, gauge_id=10, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks one chunk per variable.

| variable | dims | type | unit | description |
|---|---|---|---|---|
| amount_mm | time, gauge_id | float32 | mm | rainfall amount summed over the bin |
| n_valid | time, gauge_id | uint8 | 1 | valid native steps in the bin |
| peak_1min_rate_mm_h | time, gauge_id | float32 | mm h-1 | largest native 1 min amount * 60 in the bin |
| rate_mm_h | time, gauge_id | float32 | mm h-1 | mean rain rate over the bin, amount_mm / dt_hours |

### gauges_city_1h

City gauges, 1 h sums. Resolution 1h, dims time=2208, gauge_id=10, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks one chunk per variable.

| variable | dims | type | unit | description |
|---|---|---|---|---|
| amount_mm | time, gauge_id | float32 | mm | rainfall amount summed over the bin |
| n_valid | time, gauge_id | uint8 | 1 | valid native steps in the bin |
| peak_1min_rate_mm_h | time, gauge_id | float32 | mm h-1 | largest native 1 min amount * 60 in the bin |
| rate_mm_h | time, gauge_id | float32 | mm h-1 | mean rain rate over the bin, amount_mm / dt_hours |

### gauges_smhi_15min

SMHI gauge GbgA, 15 min (native). Resolution 15min, dims time=8832, gauge_id=1, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks one chunk per variable.

| variable | dims | type | unit | description |
|---|---|---|---|---|
| amount_mm | time, gauge_id | float32 | mm | rainfall amount summed over the bin |
| n_valid | time, gauge_id | uint8 | 1 | valid native steps in the bin |
| qc_worst | time, gauge_id | uint8 |  | largest original SMHI qc code in the bin (0 = none) |
| rate_mm_h | time, gauge_id | float32 | mm h-1 | mean rain rate over the bin, amount_mm / dt_hours |

### gauges_smhi_1h

SMHI gauge GbgA, 1 h sums. Resolution 1h, dims time=2208, gauge_id=1, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks one chunk per variable.

| variable | dims | type | unit | description |
|---|---|---|---|---|
| amount_mm | time, gauge_id | float32 | mm | rainfall amount summed over the bin |
| n_valid | time, gauge_id | uint8 | 1 | valid native steps in the bin |
| qc_worst | time, gauge_id | uint8 |  | largest original SMHI qc code in the bin (0 = none) |
| rate_mm_h | time, gauge_id | float32 | mm h-1 | mean rain rate over the bin, amount_mm / dt_hours |

### met_1h

SMHI weather at GbgA, hourly (instantaneous / as in the file; see each variable). Resolution 1h, dims time=2208, gauge_id=1, 2015-06-01 01:00:00 to 2015-09-01 01:00:00 (exclusive), chunks one chunk per variable.

| variable | dims | type | unit | description |
|---|---|---|---|---|
| humidity | time, gauge_id | float32 | % | relative humidity at 2 m |
| pressure | time, gauge_id | float32 | hPa | air pressure reduced to sea level |
| qc_humidity | time, gauge_id | uint8 |  | data-provider (SMHI) quality code |
| qc_pressure | time, gauge_id | uint8 |  | data-provider (SMHI) quality code |
| qc_temperature | time, gauge_id | uint8 |  | data-provider (SMHI) quality code |
| qc_wind_direction | time, gauge_id | uint8 |  | data-provider (SMHI) quality code |
| qc_wind_gust | time, gauge_id | uint8 |  | data-provider (SMHI) quality code |
| qc_wind_speed | time, gauge_id | uint8 |  | data-provider (SMHI) quality code |
| temperature | time, gauge_id | float32 | degC | air temperature at 2 m |
| wind_direction | time, gauge_id | float32 | degree | wind direction at 10 m (from, degrees true) |
| wind_gust | time, gauge_id | float32 | m s-1 | gust wind speed at 10 m |
| wind_speed | time, gauge_id | float32 | m s-1 | wind speed at 10 m |

### network_reference_5min

network reference series: radar over the link domain + City gauges, state and event id per 5 min bin. Resolution 5min, dims time=26496, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks one chunk per variable.

| variable | dims | type | unit | description |
|---|---|---|---|---|
| domain_mean_rain_mmh | time | float32 | mm h-1 | mean floored radar rate over the domain cells (mm/h) |
| event_id | time | <U14 |  | event id ('' outside events) |
| max_city_rate_mmh | time | float32 | mm h-1 | largest 5 min City gauge rate (mm/h) |
| n_city_gauges_valid | time | int16 | 1 | City gauges with a valid 5 min amount |
| n_city_gauges_wet | time | int16 | 1 | City gauges with >= 0.1 mm in the 5 min bin |
| radar_n_wet_cells | time | int16 | 1 | number of wet domain cells |
| radar_valid_cell_fraction | time | float32 | 1 | share of the 720 domain cells with radar data |
| radar_wet_area_fraction | time | float32 | 1 | share of the valid domain cells with floored radar rate >= 0.1 mm/h |
| state | time | int8 |  | reference state of the bin |

### link_reference_state_5min

per-link reference state (wet / dry / unknown) and reasons, from references only. Resolution 5min, dims time=26496, cml_id=364, 2015-06-01 00:00:00 to 2015-09-01 00:00:00 (exclusive), chunks [8928, 64].

| variable | dims | type | unit | description |
|---|---|---|---|---|
| reasons | time, cml_id | uint8 |  | why the state was given (bitmask) |
| state | time, cml_id | int8 |  | per-link reference state (references only, never CML) |

