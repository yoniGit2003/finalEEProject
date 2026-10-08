# CATALOG: every output file of WAA_estm

One entry per output file. Paths are relative to `research/WAA_estm/`. Commands run from the repository root with
`C:\Users\yoni\miniconda3\python.exe` (written `python` below). Folder locations are set in `WAA_estm_paths.py`.
Settings come from section `waa` of `research/config.yaml`. Status of every entry: **current**.

**Levels** (ORGANIZATION.md section 5): every file's level and purpose is in `outputs/MANIFEST.csv`. Only **key** and
**product** files are kept and committed; **diag** files (most step figures, tables and reports below) are
deleted by /tidy and listed in `outputs/TIDY_LOG.md`. An entry below describes what its script writes, so a
diag file is back after rerunning that script. Decisions from deleted reports: `DECISIONS.md`.

## Shared vocabulary

| term | meaning |
|---|---|
| A, `A_db` | A(t) = `tl_mean`(t) - `baseline_db` (dB), frozen baseline v2 (`research/baseline_estm/outputs/products/baselines.csv`) |
| tail | one (event, sub-link) pair from t0 to t_end; definition in README.md, "Tail definition" |
| t0 | `last_wet_end` of baselines.csv: the first non-wet link_state bin after the link's last wet bin in the event (UTC) |
| t_first_dry | the first link_state == dry bin at or after t0 |
| t_end | exclusive end: next wet bin of the link / start of a NaN run of A >= 15 min / t0 + 12 h / season end |
| `ended_by` | `next_wet`, `gap`, `cap`, `season_end` |
| link_state | 0 dry, 1 wet, -1 unknown (release `link_reference_state_5min`, per link, 5 min) |
| margin | the bins from t0 to t_first_dry (link_state unknown: rain within 30 min, or 0 < R_max < 0.1 mm/h) |
| wet phase | from the link's first wet bin in the event to t0 |
| A_start bins | `< 0.5`, `0.5-1.0`, `1.0-2.0`, `> 2.0` dB (`waa.report.a_start_bins_db`) |

## outputs/products/: deliverables (read by the later WAA steps)

### `outputs/products/tails_1min.parquet`
- **What**: every minute of every tail (long format, 12.4 M rows, zstd Parquet ~47 MB, **git-ignored**: rebuild with
  step 1). Tails of length 0 have no rows.
- **Answers**: "what does A do after the rain, minute by minute, and what do the references and flags say then?"
- **Columns**: `tail_id` (int, key to `tails_index.csv`), `sublink_id`, `cml_id`, `event_id`, `t` (UTC, minute start),
  `minutes_since_t0`, `A_db` (dB, float32; NaN = no valid 1 min `tl_mean`), `n_valid_tl` (valid 10 s samples of the
  6), `link_state` (5 min value on each of its minutes), `R_max_mmh`, `R_mean_mmh` (radar along the path, mm/h, 5 min
  broadcast), `gauge_rate_mmh` (nearest gauge if within 1 km: City 5 min / SMHI 15 min rate broadcast, mm/h; NaN
  otherwise), `flag_spike`, `flag_tsl_step`, `flag_network_gap`, `flag_missing` (missing_rsl or missing_tsl), decoded from
  the 1 min `flags_or`, `dew_risk` (hourly GbgA dew risk), `in_dry_jump` (inside a `dry_jumps.csv` interval of the
  sub-link).
- **Producer**: `python research/WAA_estm/step1_extract_tails.py`. **Inputs**: release (`cml_1min`,
  `link_reference_state_5min`, `radar_along_links_5min`, `gauges_city_5min`, `gauges_smhi_15min`, `links`, `sublinks`),
  baseline products `baselines.csv`, `dew_risk_hourly.csv`, `dry_jumps.csv`.

### `outputs/products/tails_index.csv`
- **What**: one row per tail (20 355 rows; floats rounded to 4 decimals).
- **Answers**: "which tails exist, how long, how did they end, how strong was the rain before, and which flags apply?"
- **Columns**:
  - ids: `tail_id`, `event_id`, `sublink_id`, `cml_id`, `split_event`;
  - link: `frequency_ghz`, `polarization`, `length_m`, `length_class`, `frequency_band`, `site_0_antenna_diameter_m`,
    `site_1_antenna_diameter_m`, `gauge_id` / `gauge_distance_m` (nearest gauge if within 1 km, else empty / NaN);
  - timing: `t0`, `t_first_dry`, `t_end` (UTC), `length_min`, `ended_by`, `first_dry_after_end` (no dry bin before
    t_end), `state_at_t0`;
  - wet phase: `wet_start`, `wet_span_min` (t0 - first wet bin), `wet_bins_min` (wet bins x 5), `wet_radar_R_mean_sum_mm`
    (path-mean radar rain, mm), `wet_radar_R_max_max_mmh`, `wet_gauge_total_mm`, `wet_A_max_db` (largest 1 min A);
  - start: `A_start_db` (mean A over [t0, t0 + 5 min)); margin: `margin_min`, `margin_R_max_max_mmh`,
    `margin_gauge_rate_max_mmh`;
  - in-tail shares of minutes: `frac_dew_risk`, `frac_dry_jump`, `frac_tsl_step`, `frac_spike`, `frac_nan`,
    `frac_unknown`;
  - weather at t0 (GbgA hour containing t0): `t0_temperature_c`, `t0_humidity_pct`, `t0_wind_speed_ms`,
    `t0_dew_point_depression_c`, `t0_sun_elevation_deg`, `t0_night`;
  - all columns of `baselines.csv` (`significant`, `run_count`, `event_start`, `event_end`, `link_has_wet_slot`,
    `last_wet_end`, `status`, `baseline_db`, `window_h`, `n_usable_slots`, `missing`, `summer_dry_median_db`, the flags
    and shares, `offset_std`, `offset_p90_abs`, `use_in_waa_main`).
- **Producer / inputs**: as `tails_1min.parquet`.

### `outputs/products/tails_selected.csv`
- **What**: `tails_index.csv` (all 20 355 tails, train AND val) plus the step 2 columns (~13 MB, **git-ignored**:
  rebuild with step 2). Val rows carry the same columns but no step 2 statistic uses them.
- **Answers**: "which tails are clean enough for the WAA model, over which time span, and what do they look like?"
- **New columns**: `exclusion_reason` (`;`-separated: `wet_t0`, `length_0`, `A_start_nan`, `not_use_in_waa_main`,
  `status_not_good`, `nan_share_gt_max`), `excluded`; `t_use_end` (UTC, exclusive), `truncated_by` (`tsl_step`,
  `dry_rise`, or the step 1 `ended_by`), `usable_min`; `margin_minutes` (t_first_dry - t0), `A_dry_db` (mean A over
  [t_first_dry, + 5 min)), `fast_drop_db` (A_start - A_dry), `fast_drop_share` (/ A_start), `margin_R_mean_sum_mm`
  (path-mean radar rain in the margin, mm), `R_max_after_dry_max_mmh` (largest R_max between t_first_dry and
  t_use_end); `A_<k>min_db` for k = 10, 30, 60, 180, 360 (mean A over [k, k + 5) min, usable part), `ratio_<k>min`
  (/ A_start), `t_half_min`, `t_e_min` (first minute the 5 min trailing mean of A < A_start / 2, / e; NaN = not
  reached), `plateau` (a 45 min window within 0.15 dB while A > 0.6 dB), `plateau_start_min`; `L_km`, `band` (the
  sub-link's own frequency < / >= 25 GHz), `antennas_m` (smaller / larger diameter), `margin_rain` (`R_max = 0` /
  `R_max > 0` in the margin); `in_main_set`, `in_alt_set`, `in_control`.
- **Step 2b columns** (added by `python research/WAA_estm/step2b_fix_tails.py`; run it after step 2):
  - start: `t0_trace` (UTC; end of the last rain-trace bin in the margin, or t0), `margin_trace_bins`,
    `A_start_trace_db` (mean A over [t0_trace, + 5 min));
  - ends: `t_use_end_lenient` / `truncated_by_lenient` / `usable_lenient_min` and `t_use_end_strict` /
    `truncated_by_strict` (adds `rain_trace`) / `usable_strict_min` (usable minutes from t0_trace);
  - level checks: `late_level_db` (median A over [t0 + 6 h, t0 + 12 h]), `late_valid_min`, `level_shift`,
    `undershoot_db` (min 15 min median within 3 h after t0_trace), `undershoot`, `dT_c` (T at t0 minus T 2 h before
    the event start, degC), `hour_local_t0_trace`;
  - noise: `offset_std_unc` (from `baseline_uncertainty.csv`, dB), `noisy_baseline` (> 0.3 dB);
  - sets: `exclusion_reason_2b` (step 2 reasons + `noisy_baseline`, `A_start_trace_nan`), `excluded_2b`,
    `in_main_strict`, `in_main_lenient`, `in_alt_strict`, `in_alt_lenient`, `in_control_2b`.
- **Producer**: `python research/WAA_estm/step2_select_tails.py`. **Inputs**: `tails_index.csv`, `tails_1min.parquet`.

### `outputs/products/tails_model.parquet`, `outputs/products/tails_model_index.csv` (FROZEN tail set v1)
- **What**: the tail set the WAA model is built on (rules: README.md, "Tail set v1"). TRAIN and VAL tails, same rules;
  VAL is not looked at. Small enough for git.
- **Answers**: "which minutes of which tails does the model see?"
- **`tails_model.parquet`** (1 min rows from t0_trace to the v1 end): `tail_set` (`v1`), `tail_id`, `sublink_id`,
  `cml_id`, `event_id`, `split_event`, `t` (UTC minute start), `minutes_since_t0_trace`, `A_db` (dB, float32; NaN = no
  valid 1 min tl_mean), `link_state` (5 min value: 0 dry, -1 unknown).
- **`tails_model_index.csv`** (one row per tail): `tail_set`, `tail_id`, `sublink_id`, `cml_id`, `event_id`,
  `split_event`, `t0`, `t0_trace`, `t_end` (UTC, exclusive), `usable_min`, `truncated_by` (`rain_trace`,
  `radar_nan`, `rise_v2`, `dry_rise`, `tsl_step`, `cap`, `next_wet`, `gap`), `A_start_trace_db`, flags `level_shift`,
  `undershoot`, `event_rise_heavy`, `in_main_set`; metadata `frequency_ghz`, `polarization`, `L_km`, `length_m`,
  `band`, `length_class`, `site_0_antenna_diameter_m`, `site_1_antenna_diameter_m`, `antennas_m`; `offset_std` (dB);
  weather at t0 (GbgA hour): `t0_temperature_c`, `t0_humidity_pct` (NaN in 102 June hours), `t0_wind_speed_ms`,
  `t0_dew_point_depression_c`, `t0_sun_elevation_deg`, `t0_night`.
- **Producer**: `python research/WAA_estm/step2c_freeze.py`. **Inputs**: `tails_selected.csv` (steps 2 + 2b),
  `tails_1min.parquet`, release radar / gauges / `link_gauge_matches`, baseline `dew_risk_hourly.csv`.

### `outputs/products/tails_model_alt.parquet`, `outputs/products/tails_model_alt_index.csv`
- **What**: the same as the frozen set, with A_start_trace >= 0.6 dB instead of 1.0 dB (`tail_set` = `v1_alt`;
  `in_main_set` marks the tails that are also in the main set). Same columns, same producer.

## outputs/tables/

| file | producer | what it contains (columns) | question |
|---|---|---|---|
| `step0_example.csv` | `step0_inventory.py` | one sub-link, one event, 5 min: `A_db`, `tl_mean_db`, `n_valid_tl`, `link_state`, `state_reasons`, `qc_flags_or` (decoded), `dew_risk_hour`, `in_dry_jump`, `network_event`, `radar_R_mean_mmh`, gauge rate | what do A, the state and the flags look like around one event? |
| `step1_a_check.csv` | `step1_extract_tails.py` | 3 random pairs: `n_minutes`, `n_nan`, `max_abs_diff_db`, `exactly_equal` (A built here vs `baseline.a_of_t`) | is A(t) built from the products identical to the baseline module's? |
| `step1_return_to_dry.csv` | `step1_extract_tails.py` | per `A_start_bin` (and `all`) and `minute` (2, 5, 10, 30, 60, 180): `n_tails`, `share_below` (A < max(0.3 dB, 2 x offset_std)) | how fast does A return to the dry level? |
| `step1_median_curves.csv` | `step1_extract_tails.py` | per `A_start_bin` and `minutes_since_t0`: `median` A (dB), `count` of tails | the typical tail shape per start level |
| `step2c_changes.csv` | `step2c_freeze.py` | `what`, `tails`, `minutes`: tails / tail time changed by the radar-NaN rule and rise cut v2; pre-rise, v1 and step 2b main-set sizes (TRAIN) | what did the last fixes cost? |
| `step2c_events.csv` | `step2c_freeze.py` | one row per TRAIN event with pre-rise main-set tails: local start / end / median t0_trace hour, `tails_before_rise_cut`, `tails_v1`, `links_v1`, `share_cut_by_rise_v2`, `share_level_flags_v1`, `radar_valid_fraction` ([start, end + 6 h], path R_max non-NaN share), mean T / RH / wind at t0, `event_rise_heavy` | which events dominate or behave oddly? |
| `step2c_median_curves.csv` | `step2c_freeze.py` | per `group` (strict v1, step 2b strict, control) and `rel` (minutes since each group's t0_trace): `median` A (dB), `count` (TRAIN) | how did the last fixes change the typical tail? |
| `step2b_tail_time.csv` | `step2b_fix_tails.py` | `quantity`, `minutes`, `share_of_step1_time`: tail time of the kept step 2 TRAIN tails before / after each step 2b change | how much tail time do the new start, the earlier rise cut and the strict end cost? |
| `step2b_noisy_sublinks.csv` | `step2b_fix_tails.py` | sub-links removed by offset_std > 0.3 dB: `tails` (in the lenient main set before the rule), `offset_std` (dB), `L_km`, `frequency_ghz`, `antennas_m` | which sub-links are too noisy in dry weather? |
| `step2b_counts.csv` | `step2b_fix_tails.py` | per `band` x `length_class` x `antennas_m`: `strict_tails/sublinks/links`, `lenient_tails/sublinks/links` (TRAIN) | what the step 2b main sets cover |
| `step2b_median_curves.csv` | `step2b_fix_tails.py` | per `group` (strict, lenient, strict no level flags, control) and `rel` (minutes since t0_trace): `median_A`, `median_A_norm` (A / A_start_trace), `count` | the typical tail after the fixes |
| `step2_truncation.csv` | `step2_select_tails.py` | per `truncated_by` (kept TRAIN tails): `tails`, `removed_min`, `median_removed_min`, `share_of_all_tail_time` | how much tail time does each truncation rule remove? |
| `step2_counts.csv` | `step2_select_tails.py` | TRAIN main set per `band` x `length_class` x `antennas_m`: `tails`, `sublinks`, `links` | which kinds of links does the main set cover? |
| `step2_correlations.csv` | `step2_select_tails.py` | per `test` (4a, 4c, 4d): `x`, `y`, `n`, `links`, Spearman `rho`, `ci_low`, `ci_high` (95 %, 1000 bootstrap resamples of links) | does fast_drop scale with L / margin rain, and do the two sub-links agree? |
| `step2_median_curves.csv` | `step2_select_tails.py` | per `group` (main set, control, main by margin rain) and `minutes_since_t0`: `median` A (dB), `count` (TRAIN, usable part) | the typical tail of the main set vs the control |
| `step1_dry_rises.csv` | `step1_extract_tails.py` | per tail with a dry-weather rise: `rise_start_min`, `rise_start_utc`, `hour_utc`, `hour_local`, `dew_risk`, `level_before_db`, `rise_db` (at the start), `rise_max_db` (largest later) | does A go up again in dry weather, and when? (definition in the step 1 report) |

## outputs/figures/

| file | producer | what is plotted / what to look for |
|---|---|---|
| `step0_example.png` | `step0_inventory.py` | top: A of one sub-link in one event with link_state shading (wet blue, unknown grey); bottom: radar R_mean and gauge rate. Look for the rain, the drop of A after the last wet bin, and dry-weather changes |
| `step2c_median_curves.png` | `step2c_freeze.py` | median A of strict v1 (solid), step 2b strict (dashed), control (dash-dot), 0-90 min and 0-6 h, linear and log y, labelled at the line ends |
| `step2c_event_<event_id>.png` | `step2c_freeze.py` | 6 random pre-rise main-set tails of each event in `waa.freeze.look_events` (seed 4), first 6 h: A with labelled line-style markers (t0, t0_trace, first dry, end v1) and hatched link_state; radar R_max (solid), R_mean (dashed), gauges within 1 km (dotted), radar NaN (x); hourly RH. Is the rise cut removing rain or a dry-weather rise? |
| `step2b_offset_std.png` | `step2b_fix_tails.py` | histogram of offset_std over the lenient main-set sub-links, dashed line at the 0.3 dB exclusion |
| `step2b_undershoot_vs_weather.png` | `step2b_fix_tails.py` | undershoot vs temperature change (flagged tails as x), flagged share by temperature change and by local hour (hatched bars, n on top) |
| `step2b_median_curves.png` | `step2b_fix_tails.py` | median A and A / A_start aligned at t0_trace, 0-90 min and 0-6 h, linear and log y; strict solid, lenient dashed, strict without level flags dotted, control dash-dot, labelled at the line ends; dotted vertical = median first dry |
| `step2b_examples_level_shift.png`, `step2b_examples_undershoot.png` | `step2b_fix_tails.py` | 6 random flagged TRAIN lenient main-set tails each (seed 3), 12 h / 6 h, same panel style as the gallery |
| `step2b_gallery.png` | `step2b_fix_tails.py` | 24 random TRAIN strict main-set tails (seed 2), first 6 h, 1 min A. Vertical lines differ by LINE STYLE and carry a text label (t0, t0_trace, first dry, end strict, end lenient); hatching /// unknown, xxx wet; flags in the panel title; radar R_max and rain-trace ticks below |
| `step2_fast_drop_vs_length.png` | `step2_select_tails.py` | per band: fast_drop (tails + binned medians) and the reference wet_A_max vs L on log-log with a slope-1 line; fast_drop / A_start vs L. Path rain should follow slope 1, WAA stay flat |
| `step2_margin_rain.png` | `step2_select_tails.py` | fast_drop histograms and 0-90 min median curves for margin R_max = 0 vs > 0; fast_drop vs margin radar rain (symlog) |
| `step2_sublink_pairs.png` | `step2_select_tails.py` | sub-link 1 vs 2 of the same link and event: A_start, fast_drop, A at t0 + 60 / 180 min, with the 1:1 line and Spearman CI |
| `step2_shape.png` | `step2_select_tails.py` | boxes of A(t0 + k) / A_start, histograms of t_half / t_e, plateau share vs L and by margin rain |
| `step2_median_curves.png` | `step2_select_tails.py` | median A of the main set vs control, 0-12 h and 0-90 min, linear and log y; dotted = median t_first_dry |
| `step2_gallery.png` | `step2_select_tails.py` | 24 random TRAIN main-set tails (seed 1), first 6 h, 1 min A, link_state shading, t_first_dry (green) and t_use_end (orange), radar R_max below: a manual check of the rules |
| `step1_distributions.png` | `step1_extract_tails.py` | histograms (log count) of `A_start_db` and `wet_A_max_db`: most tails start near 0 dB |
| `step1_median_curves.png` | `step1_extract_tails.py` | median A vs hours since t0 per A_start bin, linear and log y: the speed and shape of the decay; the far end mixes fewer tails |
| `step1_resolution_panels.png` | `step1_extract_tails.py` | 12 random tails with A_start >= 1 dB (seed 0), first 2 h: A at 10 s / 1 min / 5 min, link_state shading, radar R_max below. Does 1 min keep the decay that 10 s shows, and how much quantization noise does each resolution carry? |

## outputs/reports/

| file | producer | content |
|---|---|---|
| `step0_inventory.md` | `step0_inventory.py` (overwritten each run) | inventory of every input: paths, formats, producers, dims, spans, coverage, flags, splits, references, metadata, weather; one example; list of gaps |
| `step2c_freeze.md` | `step2c_freeze.py` (overwritten each run) | changes by the last fixes, event table and look figures, counts of the frozen set (train / val, with events), tails per event, median curves, unclear points |
| `step2b_selection.md` | `step2b_fix_tails.py` (overwritten each run) | new start, two ends, tail time, noise exclusion, level flags, main sets, median curves, gallery, unclear points (TRAIN only) |
| `step2_selection.md` | `step2_select_tails.py` (overwritten each run) | exclusions, truncation, main-set counts, fast-drop test with bootstrap CIs, shape descriptors, median curves, gallery, unclear points (TRAIN only) |
| `step1_tails.md` | `step1_extract_tails.py` (overwritten each run) | counts, distributions, return to the dry level, median curves, dry-weather rises, resolution check, list of unclear points |

## .cache/ (git-ignored, safe to delete)

`step0_run.log`, `step1_run.log`, `step2_run.log`, `step2b_run.log`, `step2c_run.log`, `check_new.py` (organization check incl. untracked files): console output of the last runs.
