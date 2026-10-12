# CATALOG: every output file of WAA_estm

One entry per output file. Paths are relative to `research/WAA_estm/`. Commands run from the repository root with
`C:\Users\yoni\miniconda3\python.exe` (written `python` below). Folder locations are set in `WAA_estm_paths.py`.
Settings come from section `waa` of `research/config.yaml`. Status of every entry: **current**.

**Levels** (ORGANIZATION.md section 5): every file's level and purpose is in `outputs/MANIFEST.csv`. Only **key** and
**product** files are kept and committed; **diag** files (most step figures, tables and reports below) are
moved by /tidy to the local `legacy/` folder (same sub-folder and name) and listed in `outputs/TIDY_LOG.md`. An entry below describes what its script writes, so a
diag file is back after rerunning that script. Decisions from the step reports: `DECISIONS.md`.

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

### `outputs/products/sim_tails_m1.parquet`, `sim_tails_m3.parquet`, `sim_tails_m4.parquet`, `sim_tails_m5.parquet`, `outputs/products/sim_tails_index.csv`
- **What**: simulated tail sets on the TRAIN v1 skeleton (same tail_id, sub-link, link, event, usable length,
  A_start_trace), one per drying law (m1 exponential, m3 nonlinear, m4 two layers, m5 linear + quadratic), each with
  its tail-to-tail time-scale spread and the full observation model (offset, 10 s white noise, 1 min AR(1) slow noise,
  0.3 dB quantization, real NaN minutes). Laws, parameters, spread: `outputs/tables/step3a_calibration.csv`.
- **Answers**: "what would tail set v1 look like if law X were true?" (input of the identifiability step 3b).
- **`sim_tails_<law>.parquet`** (1 min rows): `law`, `tail_id`, `sublink_id`, `cml_id`, `event_id`, `split_event`
  (`train`), `minutes_since_t0_trace`, `W_db` (truth, mean of the six 10 s values, dB), `A_db` (observed 1 min A, dB).
- **`sim_tails_index.csv`** (per law and tail): `law`, `tail_id`, `f_time_scale` (the tail's time-scale factor),
  `offset_db`, `sigma_white_db`, `sigma_slow_db`, `tau_noise_min`, `phase_db` (quantization grid phase), `sublink_id`,
  `cml_id`, `event_id`, `usable_min`, `A_start_trace_db`.
- **Producer**: `python research/WAA_estm/step3a_simulator.py` (library `waa_sim.py`). **Inputs**:
  `tails_model_index.csv`, `tails_model.parquet`, `sim_noise_params.csv` (same run).

### `outputs/products/sim_noise_params.csv`
- **What**: the noise parameters per sub-link of the TRAIN v1 skeleton, fitted on its real control tails (step 2b
  control, A_start < 0.5 dB) through a lookup of simulated dry series (blocks = the median control tail length).
- **Columns**: `sublink_id`, `n_control_tails`, `n_control_min`, `median_control_len_min`; real targets
  `within_std_db` (within-tail std of 1 min A), `on_level_share` (share of minutes with 1 min `tl_std` = 0 and 6 valid
  samples), `acf_1`, `acf_10`, `acf_60` (autocorrelation of 1 min A at lags in minutes, within-tail, pooled); fitted
  `sigma_white_db` (10 s), `sigma_slow_db`, `tau_noise_min` (AR(1) at 1 min); `sim_<target>` (the same targets of the
  chosen grid point); `source` (`fit` / `median_fallback`).
- **Producer / inputs**: `step3a_simulator.py`; `tails_selected.csv`, `tails_1min.parquet`, release `cml_1min`.

## outputs/tables/

| file | producer | what it contains (columns) | question |
|---|---|---|---|
| `step0_example.csv` | `step0_inventory.py` | one sub-link, one event, 5 min: `A_db`, `tl_mean_db`, `n_valid_tl`, `link_state`, `state_reasons`, `qc_flags_or` (decoded), `dew_risk_hour`, `in_dry_jump`, `network_event`, `radar_R_mean_mmh`, gauge rate | what do A, the state and the flags look like around one event? |
| `step1_a_check.csv` | `step1_extract_tails.py` | 3 random pairs: `n_minutes`, `n_nan`, `max_abs_diff_db`, `exactly_equal` (A built here vs `baseline.a_of_t`) | is A(t) built from the products identical to the baseline module's? |
| `step1_return_to_dry.csv` | `step1_extract_tails.py` | per `A_start_bin` (and `all`) and `minute` (2, 5, 10, 30, 60, 180): `n_tails`, `share_below` (A < max(0.3 dB, 2 x offset_std)) | how fast does A return to the dry level? |
| `step1_median_curves.csv` | `step1_extract_tails.py` | per `A_start_bin` and `minutes_since_t0`: `median` A (dB), `count` of tails | the typical tail shape per start level |
| `step3a_calibration.csv` | `step3a_simulator.py` (**key**) | per `law` (m1, m3, m4, m5): `description`, `n_free_params` (incl. the spread), the law's parameters (`tau`; `c`, `n`; `tau_f`, `tau_s`, `phi`; `a`, `b`; minutes, dB), `s` (spread), `f_10_90` (time-scale factor 10-90 %), `s_at_grid_edge`, `rmse_median_truth` (noise-free median vs real), `rmse_median_observed`, `rmse_band_observed` (observed simulation vs real median / 10-90 % band; units of A / A_start) | which laws can reproduce the real median curve, and with which parameters? |
| `step3b1b_fits.csv` | `step3b1b_pooled_fit.py` | one row per pooled fit: `law` (truth), `rep`, `window` (first_60min, first_6h, whole), `candidate`, `ll` (marginal log-likelihood), `k`, `n_independent` (distinct link x event draw pairs = BIC's N), `n_tails`, `bic`, `s` (fitted spread), `p1`-`p3` with `param_names`, `nfev`, `success`, `seconds` | the raw material of the confusion, recovery and decision tables |
| `step3b1b_recovery.csv` | `step3b1b_pooled_fit.py` (**key**) | per truth `law`, `window`, `parameter` of its own candidate (and `s`): `true`, `median_ratio`, `p10_ratio`, `p90_ratio` (fitted / true over the 10 event-bootstrap replicates), `within_tol` (+-30 %) | does the pooled fit recover the law parameters? |
| `step3b1b_decision.csv` | `step3b1b_pooled_fit.py` | per `law` and `window`: `in_rule`, `own_win_share`, `win_ok` (>= 80 %), `params_ok`, `params_off`, `main_confusion`, `met` | is the pre-set decision rule met? |
| `step3b1b_population.csv` | `step3b1b_pooled_fit.py` | per `statistic` (median t_half per A_start bin; median A / A_start at 10, 30, 60, 180 min): `real`, `real_ci_low`, `real_ci_high` (events resampled), `<law>_min`, `<law>_max` over the replicates, `<law>_real_inside`, `separates` (pairs of truths with non-overlapping ranges) | which population statistics tell which laws apart? |
| `step3b1c_win_shares.csv` | `step3b1c_stress_test.py` (product) | per `test` (T1 misspecified noise, T2 honest starts, T3 faster truths, T3s = T3 with s = 1.5), `truth` (m1, m4, m5, m6) and `candidate`: `wins` (replicates with the lowest pooled BIC), `replicates`, `share`, `own` (the truth's own candidate) | who wins when the pooled fit is stressed? |
| `step3b1c_decision.csv` | `step3b1c_stress_test.py` (product) | per `test` and `truth`: `own`, `own_wins`, `replicates`, `wins_instead`, `median_dbic_own_minus_best_other` (negative = own wins), `median_dbic_where_lost`, `win_ok` (>= 4 of 5), `params_off`, `params_ok` (key parameters within +-30 %), `met` | does each cell meet the rule fixed before the run? |
| `step3b1c_recovery.csv` | `step3b1c_stress_test.py` (product) | per `test`, `truth`, `parameter` of the own candidate (and `s`, not in the rule): `true`, `median_ratio`, `min_ratio`, `max_ratio` (fitted / true over the replicates), `in_rule`, `within_tol` | are the law parameters still recovered under stress? |
| `step3b1c_start_sensitivity.csv` | `step3b1c_stress_test.py` (product) | T2, per `truth` and `candidate`: `datasets`, `starts_differ` (best and worst start differ by > 2 BIC), `median_bic_range`, `max_bic_range`, `generic_worse_than_calib` / `generic_better_than_calib` (best generic start vs the fit started at the truth, > 2 BIC), `fits_stopped_at_maxfev`, `fits`, `best_start_counts` | does the result depend on where the optimiser starts? |
| `step3b1c_recalibration.csv` | `step3b1c_stress_test.py` (product) | T3, per `test` and `truth`: `s`, `old_parameters`, `new_parameters`, and A / A_start at 10, 30, 60 min: `real_*`, `old_sim_*` / `new_sim_*` (median over the replicates of the observed fake tails; `new_sim_min_*`, `new_sim_max_*`), `old_noise_free_*` / `new_noise_free_*` (calibration values) | which faster truths were used, and do they match the real early decay? |
| `step3b1c_noise_check.csv` | `step3b1c_stress_test.py` | T1, per `statistic` (within-tail std, on-level share, ACF at 1 / 10 / 60 min): `real_median`, `sim_3afix_median`, `sim_T1_median`, `abs_err_3afix`, `abs_err_T1` (median absolute error per sub-link) | does the two-time-scale noise still look like the real dry noise? |
| `step3b1c_fits.csv` | `step3b1c_stress_test.py` | one row per pooled fit: `test`, `truth`, `rep`, `candidate`, `start` (calib, g1-g4, f1sol), `s0`, `ll`, `k`, `n_independent`, `n_tails`, `bic`, `s`, `p1`-`p3` with `param_names`, `nfev`, `success`, `seconds` | the raw material of the stress-test tables |
| `step3b2_selection.csv` | `step3b2_sindy_fake.py` (product) | per `truth`, `spread`, `window`, `offset` (naive, handled, oracle), `rescale` (none, fit, true) and library `term` (1, A, A^2, A^3, sw0.1 ... sw1 = max(A - A_c, 0)): `role` (true / allowed / spurious; `-` for m4), `n_selected` (replicates with stability frequency >= 0.6), `replicates`, `mean_frequency` | which terms does weak-form SINDy select when the truth is known? |
| `step3b2_coefficients.csv` | `step3b2_sindy_fake.py` (product) | per configuration and `coefficient` of the true law (`switch` = sw0.1 + sw0.2): `true` (1 / min, base law), `median_ratio`, `min_ratio`, `max_ratio` (fitted / true; 0 where not selected), `replicates_nonzero`, `counted` (part of the rule) | are the coefficients of the true terms recovered? |
| `step3b2_decision.csv` | `step3b2_sindy_fake.py` (product) | `kind` = rule: per truth and configuration `true_selected`, `spurious`, `coef_ratios`, `true_ok`, `spurious_ok`, `coef_ok`, `met`, `why_not`, `main` (the configuration that counts: with spread, first 60 min, handled, not rescaled); `kind` = verdict: `question` (Q1-Q4), `verdict`, `basis` | is the rule fixed before the run met, and what are the verdicts? |
| `step3b2_agreement.csv` | `step3b2_sindy_fake.py` (product) | per spread dataset (`truth`, `rep`): `pooled_winner` and `pooled_family`, `pooled_dbic_to_second`, `sindy_terms` / `sindy_family` (plain) and `sindy_rescaled_*` (time rescaled by the pooled-fit speed), `pooled_F1_rate` (1 / tau), `sindy_A_coef`, `agree`, `agree_rescaled` | does SINDy name the same family of law as the pooled fit? |
| `step3b2_fits.csv` | `step3b2_sindy_fake.py` | one row per SINDy run (dataset x window x offset x rescale): `lam` (threshold, one-SE rule), `n_rows`, `n_tails_used`, `freq_<term>`, `coef_<term>`, `sel_<term>`, `selected`, `family`; post hoc diagnostic: `lam_min`, `cv_err_*`, `fmin_<term>`, `cmin_<term>` (threshold of lowest held-out error) | the raw material of the 3b-2 tables |
| `step3b3_curves.csv` | `step3b3_curve_then_equation.py` (product) | per `kind` (main = spread + 3a-fix noise, nospread, T1 = T1 noise), `truth`, `window` (first_60min, first_6h), `rep` and grid point `w_db` (40 log-spaced points of the supported range): `k` (learned drying rate per dB of water, 1 / min), `sd_log_k` (per-replicate sd from the fit's curvature; empty if not available), `k_true` (empty for m4) | what does the free-curve pooled fit draw, against the truth? |
| `step3b3_j1_curve.csv` | `step3b3_curve_then_equation.py` (product) | per `kind` (main = spread + 3a-fix noise, nospread, T1 = T1 noise), `truth`, `window` (first_60min, first_6h): `band_coverage` (share of the range where the true k is inside the min-max band of the replicates), `curve_ok` (>= 80 %), `median_k_ratio` (learned / true), `w_hi_db`, `median_slope` (log k vs log W), `with_interval`, `flat_reps`, `rising_reps`, `step_reps`, `median_step_ratio` (m6), `shape_wanted`, `shape_ok` | jump 1: is the drawn curve right in level and shape? |
| `step3b3_j2_selection.csv` | `step3b3_curve_then_equation.py` (product) | per `kind` (main = spread + 3a-fix noise, nospread, T1 = T1 noise), `truth`, `window` (first_60min, first_6h) and library `term` (W, W^2, W^3, sw0.05 ... sw0.5 = max(W - A_c, 0)): `role` (true / allowed / spurious; `-` for m4), `n_selected`, `replicates`, `threshold` | jump 2: which terms name the learned curve? |
| `step3b3_j2_coefficients.csv` | `step3b3_curve_then_equation.py` (product) | per `kind` (main = spread + 3a-fix noise, nospread, T1 = T1 noise), `truth`, `window` (first_60min, first_6h) and `coefficient` of the true law (`switch` = sw0.1 + sw0.15 + sw0.2): `true` (1 / min), `median_ratio`, `min_ratio`, `max_ratio` (fitted / true; 0 where not selected), `replicates_nonzero`, `counted` | jump 2: are the coefficients right? |
| `step3b3_menu_comparison.csv` | `step3b3_curve_then_equation.py` (product) | per main dataset (`truth`, `rep`): `menu_winner_60min`, `menu_dbic_to_second`, `curve_minus_best_menu_bic_60min`, `curve_minus_f4_bic_60min`, `curve_minus_f4_bic_6h` (m4 only), `both_intervals`, `windows_disagree_share`, `median_abs_log_ratio_6h_60min`, `warn_windows`, `warn_f4`, `warning` | does the free curve agree with the menu fit, and does m4 raise a warning sign? |
| `step3b3_decision.csv` | `step3b3_curve_then_equation.py` (product) | per `truth`: `window_used` (the better window), `curve_ok`, `shape_ok`, `equation_ok`, `verdict`, `basis`; last row `overall` | the pre-set rule and the verdicts |
| `step3b3_posthoc_upper_range.csv` | `step3b3_curve_then_equation.py` (product) | post hoc diagnostic (no verdict uses it), per `kind`, `truth`, `window`: the learned curves judged only above the second knot (`w_from_db` - `w_to_db`): `band_coverage`, `median_k_ratio`, `true_slope`, `slope_min` / `slope_median` / `slope_max` over the replicates, `equation` (jump 2 on that range), `coef_ratios` | is the curve right where it is well determined? |
| `step3b3_fits.csv` | `step3b3_curve_then_equation.py` | one row per free-curve fit (dataset x window): `n_tails`, `n_independent`, `w_hi_db`, `knots_db`, `lambda`, `s`, `ll`, `bic`, `best_start`, `start_cost_spread`, `starts_within_2`, `converged`, `nfev`, `seconds`, `has_interval`, `slope` with `slope_lo` / `slope_hi`, `k_knot1`-`k_knot5`, `step_ratio` (m6) | the raw material of the 3b-3 tables |
| `step3b3b_curves.csv` | `step3b3b_confirm.py` (product) | per `kind` (main = spread + 3a-fix noise, T1 = T1 noise), `truth`, `rep` and grid point `w_db` (40 log-spaced points of the judged range, 0.35 dB to W_hi): `k` (learned drying rate per dB of water, 1 / min), `k_true` (empty for m4) | what does the free-curve fit draw on fresh data, against the truth? |
| `step3b3b_curve_shape.csv` | `step3b3b_confirm.py` (product) | per `kind` (main = spread + 3a-fix noise, T1 = T1 noise), `truth`: `w_lo_db`, `w_hi_db` (judged range common to the replicates), `band_coverage` (share of the range where the true k is inside the min-max band of the replicates), `median_k_ratio` (learned / true), `true_slope`, `slope_min` / `slope_median` / `slope_max` and `slopes` (slope of log k vs log W per replicate), `flat_reps` (|slope| < 0.15), `rising_reps` (slope > 0.6), `shape_wanted`, `curve_ok`, `shape_ok` (m1, m3, m5 only) | is the drawn curve right in level and shape above 0.35 dB? |
| `step3b3b_j2_selection.csv` | `step3b3b_confirm.py` (product) | per `kind` (main = spread + 3a-fix noise, T1 = T1 noise), `truth` and library `term` (W, W^2, W^3, sw0.5, sw0.75 = max(W - x dB, 0)): `role` (true / allowed / spurious; `-` for m4, m6), `n_selected`, `replicates`, `threshold` | which terms name the learned curve? |
| `step3b3b_j2_coefficients.csv` | `step3b3b_confirm.py` (product) | per `kind` (main = spread + 3a-fix noise, T1 = T1 noise), `truth` and `coefficient` of the true law: `true` (1 / min), `median_ratio`, `min_ratio`, `max_ratio` (fitted / true; 0 where not selected), `replicates_nonzero`, `counted` | are the coefficients right? |
| `step3b3b_m4_bic.csv` | `step3b3b_confirm.py` (product) | per replicate of truth m4: `bic_free_curve`, `bic_f4` (best of 3 starts), `f4_best_start`, `curve_minus_f4_bic`, `warning` (> 10) | does a two-layer truth raise the warning sign? |
| `step3b3b_decision.csv` | `step3b3b_confirm.py` (product) | per `truth`: `curve_ok`, `shape_ok`, `equation_ok`, `verdict`, `basis` (m6: not judged); last row `overall` | the rule fixed before the run and the verdicts |
| `step3b3b_fits.csv` | `step3b3b_confirm.py` | one row per free-curve fit: `n_tails`, `n_independent`, `w_lo_db`, `w_hi_db`, `minutes_in_range`, `knots_db`, `lambda`, `s`, `ll`, `bic`, `best_start`, `start_cost_spread`, `nfev`, `seconds`, `slope`, `k_knot1`-`k_knot5` | the raw material of the 3b-3b tables |
| `step3b3c_a_jump2_alone.csv` | `step3b3c_method.py` (product) | part A, per `truth` (m1, m3, m5, m6), `noise` (relative, 0.02 / 0.05) and `method` (best-subset + BIC, E-SINDy with waa_sindy, E-SINDy with pysindy STLSQ): `exact` (draws with exactly the true terms), `true_included` (true terms present, extras allowed), `draws`, `most_common` answer, `coef_ratios` (median fitted / true over the exact draws), `ok` (rule: exact >= 16 of 20, coefficients within +-30 %) | does jump 2 name the true formula when it is given the true curve plus a little noise? |
| `step3b3c_a_draws.csv` | `step3b3c_method.py` (product) | part A, per `truth`, `noise`, `draw`, `method`: `returned` terms, `exact`, `true_included`, `ratio_<coefficient>` (fitted / true) | the raw material of the part A table |
| `step3b3c_b_step.csv` | `step3b3c_method.py` (product) | part B, per 3b-3b dataset (`truth`, `rep`): `bic_no_step`, `bic_step`, `bic_gain`, `step_kept` (gain > 6), `wc_db` (fitted step position), `height_ln` (step height, ln units), `k_low` (1 / min, the constant rate below W_c), `best_wc_start`, `seconds` | does the free curve with an optional step find the switch of m6 and no step for m1, m3, m5? |
| `step3b3c_c_candidates.csv` | `step3b3c_method.py` (product) | part C, per dataset and `candidate` (J2-1 / J2-2 / J2-3 = jump-2 formulas of 1-3 terms, step, F1, F3, F4, F5, F6, free curve as reference): `what` (terms or fitted parameters), `family`, `parameters`, `held_out_ll` (sum over the two halves of the event draws), `held_out_minus_best`, `ll` and `bic` on the whole dataset | how does every candidate law score on events it was not fitted to? |
| `step3b3c_c_winners.csv` | `step3b3c_method.py` (product) | part C, per dataset: `winner` (best held-out log-likelihood), `winner_what`, `family`, `in_true_family`, `margin_to_other_family`, `best_other_family`, `simplest_within_2`, `bic_winner`, `bic_winner_family`, `free_curve_minus_winner`, `step_kept`, `j2_bic_pick` | does jump 3 choose the true family? |
| `step3b3c_d_bootstrap_band.csv` | `step3b3c_method.py` (product) | part D, per `truth` (m1, m3, m5) and `dataset` (0-2): `w_hi_db`, `bootstraps`, `coverage_band` (share of [0.35 dB, W_hi] where the true k is inside the 5-95 % band of the event-bootstrap refits), `coverage_minmax`, `band_half_width_pct`, `minmax_half_width_pct`, `fit_k_ratio`, `slope_fit`, `slope_lo`, `slope_hi`, `true_slope`, `seconds`, `ok` | does a band made from ONE dataset cover the truth, and how wide is it? |
| `step3b3c_d_band_curves.csv` | `step3b3c_method.py` (product) | part D, per `truth`, `dataset` and grid point `w_db`: `k_fit`, `k_lo`, `k_hi` (5-95 %), `k_min`, `k_max`, `k_true` (1 / min) | the bands themselves |
| `step3b3c_e_f4_check.csv` | `step3b3c_method.py` (product) | part E, per m4 replicate: `bic_free_curve`, `bic_f4_3b3b` (3 generic starts, step 3b-3b), `bic_f4_generic_now`, `bic_f4_f1sol`, `best_start`, `curve_minus_f4_3b3b`, `curve_minus_f4_now`, `warning_now` (> 10) | was the failed m4 warning of 3b-3b caused by a badly started F4? |
| `step3b3c_decision.csv` | `step3b3c_method.py` (product) | per `part` (A-D) and `truth`: `verdict`, `basis`; then one row per part and `overall` | the rules fixed before the run and the verdicts |
| `step3b4_step.csv` | `step3b4_quantization.py` (product) | per `setting` (RSL quantization none / 0.3 dB / 0.1 dB) and `rep` (truth m1): `step_kept` (BIC gain > 6), `wc_db` (fitted step position), `bic_gain` = `bic_no_step` - `bic_step`, `height_ln` (step height, ln units), `best_wc_start` (dB), `lambda` (penalty chosen by cross-validation), `w_hi_db`, `no_step_k_ratio_below` / `no_step_k_ratio_above` (fitted / true k of the curve without the step, median over 0.1-0.35 dB and over 0.35 dB-W_hi), `no_step_slope_above` (true 0), `seconds` | is the false step of 3b-3c still drawn when the fake RSL is not quantized? |
| `step3b4_decision.csv` | `step3b4_quantization.py` (product) | per `setting`: `step_kept` (n of 5), `median_wc_kept_db`, `wc_kept_db`, `median_bic_gain`, `bic_gain_range`, the two k ratios; last row: `verdict` of the pre-set rule with its basis | the verdict: quantization is / is not the cause |
| `step4_menu_fits.csv` | `step4_real_drying.py` (product) | per `dataset` (main, no_rise_heavy, no_shift_undershoot, alt, loo_<event>), `window_min` (60; 360 for main) and `candidate` (F1, F3, F4, F5, F6, and `F1+` = F1 with a start W0 >= 0): `family`, `n_tails`, `n_independent` (N of the BIC), `k` (parameters), `ll`, `bic`, `dbic` (to the best), `dbic_like_for_like` (among F1+, F3, F5, F6), `parameters` (text; tau in min, rates in 1 / min, A_c in dB), `s` (speed spread), `best_start`, `starts`, `bic_spread_starts`, `theta` (json) | which menu law do the real tails prefer, and by how much? |
| `step4_summary.csv` | `step4_real_drying.py` (product) | per `dataset`: `tails`, `events`, `links`, `tails_big3`, `w_hi_db`, `minutes_in_range`, menu reading as run (`menu_winner`, `menu_family`, `runner_up`, `margin_bic`, `next_family`, `margin_family_bic`, `f1_minus_winner_bic`), like for like (`f1_negative_start_tails`, `like_for_like_winner`, `like_for_like_margin_family_bic`, `like_for_like_f1_minus_winner_bic`, `f3_n_at_lower_bound`), `f1_tau_min`, `s_winner`, `curve_slope` (log k vs log W over the set's own range), `curve_slope_common_range`, `curve_reading`, `curve_s`, `curve_penalty`, `curve_minus_winner_bic`, `jump2_esindy`, `jump2_family`, `agree` | do the result and the agreement of the three tools hold in every sensitivity set and without each big event? |
| `step4_curves.csv` | `step4_real_drying.py` (product) | per `dataset` and `w_db` (40 log-spaced points, 0.35 dB to W_hi): `k` (1 / min); main set also `k_lo`, `k_hi` (5-95 % of the 20 event-bootstrap refits), `k_min`, `k_max` | the real drying rate per dB of water, and how uncertain it is |
| `step4_jump2.csv` | `step4_real_drying.py` (product) | per `dataset` and `tool` (E-SINDy; best-subset with 1, 2, 3 terms; PySR's three simplest formulas and its own pick, sensitivity sets only): `formula` for h(W) = W k(W) in dB / min, `terms`, `family`, `relative_rms_misfit` (formula / curve - 1), `note` | is there a short formula for the real curve? |
| `step4_fit_quality.csv` | `step4_real_drying.py` (product) | per fitted `law` and `window_min`: `within_tail_std_db` and `acf_1`, `acf_10`, `acf_60` of the residuals, the same for dry noise drawn from the noise model (`noise_model_*`), `whitened_variance_ratio` (1 = residuals are dry noise); last row: the real control tails | are the residuals of the fit as small as dry-weather noise? |
| `step4_fake_comparison.csv` | `step4_real_drying.py` (product) | per `source` (real tails; each fake truth of the 3b-3b datasets; each recalibrated truth of 3b-1c T3): `slope` and `slope_range`, `ratio_<k>min` = median A / A_start at 10 / 30 / 60 min with `range_<k>min`, `slope_inside`, `real_inside_<k>min` | do the real tails look like any of the fake data the tools were tested on? |
| `step3b1_confusion.csv` | `step3b1_fitting_on_sim.py` | per `score` (bic, aic), `law` (truth), `variant` (full or an ablation), `candidate`: `per_tail_share` (mean over replicates), `per_tail_share_sd`, `per_dataset_share` (share of replicates won by the summed score) | which candidate does the method choose for each true law? |
| `step3b1_recovery.csv` | `step3b1_fitting_on_sim.py` (diag since 2026-10-10) | per `law` and `parameter` of the truth's own candidate: `measure` (fitted / true for scale parameters; fitted - true for phi, n; fitted - true offset for c; fitted - A_start for W0), `true_value`, `n`, `median`, `p10`, `p90`, `at_bound_share` | can the fit recover the per-tail parameters? |
| `step3b1_fits.parquet` | `step3b1_fitting_on_sim.py` | every per-tail fit: `tail_id`, `candidate`, `n`, `rho`, `n_eff`, `k`, `rss`, `bic`, `aic`, `converged`, the candidate's parameters (`W0`, `tau`, `tau_f`, `tau_s`, `phi`, `k`, `n`, `a`, `b`, `c`), the truth's `f_time_scale`, `offset_db`, `law`, `variant`, `rep`, `dataset` | the raw material of the confusion and recovery tables |
| `step3b1_m6_calibration.csv` | `step3b1_fitting_on_sim.py` | m6 calibration per spread `s`: `k_run`, `k_evap` (1/min), `A_c` (dB), `status`, `rmse_median_truth`, `rmse_median_observed`, `rmse_band_observed` | does the state-switch truth calibrate? |
| `step3b1_thalf_population.csv` | `step3b1_fitting_on_sim.py` | per `law` and A_start `bin`: `min`, `median`, `max` of the replicate medians of t_half (min); real `tails`, `real_median`, `real_ci_low`, `real_ci_high` (events resampled); `real_inside` | does t_half per A_start bin separate the laws? |
| `step3b1_dadt_population.csv` | `step3b1_fitting_on_sim.py` | per `law` and smoothed-A bin `A_bin` (`A_mid`, dB): `min`, `max` of the replicate median dA/dt (dB/min); `real_median_dadt`, `real_count`, `real_inside` (1 / 0 / NaN) | does the drying rate vs A separate the laws? |
| `step3b1_ablation.csv` | `step3b1_fitting_on_sim.py` | the BIC rows of `step3b1_confusion.csv` for truths m1, m4 and every variant (full, no_noise, no_quantization, no_offset, no_spread, all_off) | which observation effect causes the confusion? |
| `step3a_thalf_by_astart.csv` | `step3a_simulator.py` | per `source` (real, m1, m3, m4, m5) and `A_start_bin_db` (1-1.5, 1.5-2, 2-3, > 3): `tails`, `n_reached`, `median_t_half_min`; real also `ci_low`, `ci_high` (95 %, 1000 resamples of events) | does t_half depend on A_start (nonlinear) or not (exponential, two layers)? |
| `step3a_shape_checks.csv` | `step3a_simulator.py` | per `quantity` (`t_half_min`, `t_e_min`, `ratio_<k>min`, `last_hour_std_db`) and `source` (real, m1 ... m5): `n`, `p10`, `median`, `p90`, `not_reached` | does the simulator look like the data beyond the median? |
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
| `step3b1b_confusion.png` | `step3b1b_pooled_fit.py` (**key**) | one greyscale matrix per fitting window (first 60 min, first 6 h, whole): truth (rows m1, m3, m4, m5, m6) x chosen candidate (F1, F3, F4, F5, F6, lowest pooled BIC), share of the 10 event-bootstrap replicates in the cells, own candidate bold in a dashed box. A strong diagonal = the pooled fit tells the laws apart |
| `step3b1b_population.png` | `step3b1b_pooled_fit.py` | one panel per population statistic: replicate range per truth (bars) vs the real value (line) and its 95 % CI from resampling events (hatched); the title names the pairs of truths it separates |
| `step3b1c_win_share.png` | `step3b1c_stress_test.py` | one panel per test (T1, T2, T3, T3s): bars = replicates won by the truth's own candidate, per truth; dashed line = the rule (4 of 5); grey bar = cell passes, hatched white bar = fails, with a text label (who wins instead, parameters off) | look for hatched bars: where does the pooled fit break? |
| `step3b2_selection_frequency.png` | `step3b2_sindy_fake.py` | three greyscale grids (with spread; no spread; with spread and time rescaled), truth (rows) x library term (columns), first 60 min, handled offset; text = replicates in which the term is selected; dashed box, bold = a term of the true law; hatched = a term outside the truth selected in more than one replicate | look for dashed boxes with 0/5 (true term missed) and hatched cells (false term) |
| `step3b3_k_<truth>.png` | `step3b3_curve_then_equation.py` | one figure per truth (m1, m3, m4, m5, m6), one panel per window, log-log: learned k(W) (thick solid = median, thin solid = each replicate), band of the replicates (hatched), true k(W) (dashed, labelled; m4: the true rate along tails that start at 1 and 3 dB, dotted / dash-dotted) | flat k = exponential, rising k = nonlinear in W, a step = state switch: is the dashed line inside the hatched band? |
| `step3b3b_k_<truth>.png` | `step3b3b_confirm.py` | one figure per truth (m1, m3, m4, m5, m6), log-log, first 60 min: learned k(W) (thick solid = median, thin solid = each replicate), band of the replicates (hatched), true k(W) (dashed; m4: the true rate along tails that start at 1 and 3 dB, dotted / dash-dotted), the judged range shaded, a vertical line at 0.35 dB labelled "below: noise-dominated, not judged" | is the dashed line inside the hatched band over the shaded range? |
| `step3b3c_b_step_<truth>.png` | `step3b3c_method.py` | part B, one figure per truth (m1, m3, m5, m6), one panel per replicate, log-log: free curve without the step (thin solid), with the step (thick solid), true k (dashed), fitted W_c (dotted vertical line with a text label); the title says the BIC gain and whether the step is kept | is a step drawn only where the truth has one? |
| `step3b3c_d_band_<truth>.png` | `step3b3c_method.py` | part D, one figure per truth (m1, m3, m5), one panel per dataset: event-bootstrap band 5-95 % (hatched), bootstrap min / max (dotted), the fit (solid), true k (dashed) over the judged range | is the dashed line inside the hatched band? |
| `step3b4_m1_quantization.png` | `step3b4_quantization.py` | three panels side by side (RSL quantization none, 0.3 dB, 0.1 dB), log-log: the fit with the step of every replicate (solid, numbered at its left end; thick = step kept), true k (dashed, labelled), ticks at the bottom = fitted W_c per replicate, dotted vertical line at 0.35 dB | does the false step disappear without quantization? (it does not) |
| `step4_k_real.png` | `step4_real_drying.py` (diag; the key candidate) | log-log, main set: the free curve k(W) (thick solid; dotted below 0.35 dB = not read), the 5-95 % band of the 20 event-bootstrap refits (hatched), F1 (dashed, flat) and the menu winner if it is not F1 (dash-dotted), text labels on the lines; read range shaded. Look for: flat = exponential, rising = nonlinear; the real curve falls |
| `step4_k_sets.png` | `step4_real_drying.py` (diag) | log-log: the free curve of every set on its own read range (main thick solid, sensitivity dashed, alternative dash-dotted, leave-one-event-out dotted), labelled with the set name and slope. Look for: do the sets agree in trend and in the bumps? |
| `step4_gallery.png` | `step4_real_drying.py` (diag) | 12 random real tails (fixed seed), first 60 min: measured A (dots), fitted menu winner (solid), fitted F1 (dashed), with the tail's speed factor and offset. Look for: tails that stay flat or rise |
| `step3b1_confusion.png` | `step3b1_fitting_on_sim.py` (diag since 2026-10-10) | two greyscale matrices, truth (rows m1, m3, m4, m5, m6) x chosen candidate (F1, F3, F4, F5), numbers in the cells, the truth's own candidate bold in a dashed box: left per tail (share of tails), right per dataset (share of 10 replicates won). Look for a strong diagonal: only m1 has one |
| `step3b1_thalf_population.png` | `step3b1_fitting_on_sim.py` | one panel per truth: replicate min-max of the median t_half per A_start bin (hatched boxes) vs the real median with its 95 % CI (dots, error bars) |
| `step3b1_dadt_population.png` | `step3b1_fitting_on_sim.py` | one column per truth: drying rate -dA/dt vs smoothed A, replicate min-max (hatched) vs real (dots), linear (top) and log-log (bottom) |
| `step3b1_ablation.png` | `step3b1_fitting_on_sim.py` | truths m1, m4: share choosing each candidate per variant (rows), per tail and per dataset, own candidate bold |
| `step3a_real_vs_sim_median.png` | `step3a_simulator.py` (**key**) | rows = laws m1, m3, m4, m5 (with spread and slow noise), columns = 0-90 min / 0-6 h, linear / log: real median A / A_start (solid) with hatched 10-90 % band vs simulated median (dashed) and 10 / 90 % (dotted). Which law follows the median, which the spread? |
| `step3a_t_half_t_e.png` | `step3a_simulator.py` | histograms of the time to A_start / 2 and A_start / e: real (hatched bars) vs each law (dashed outline) |
| `step3a_noise_fit.png` | `step3a_simulator.py` | per sub-link, real vs simulated for the 5 noise targets (within-tail std, on-level share, ACF at 1 / 10 / 60 min; each its own marker), 1:1 dashed; title = median absolute error |
| `step3a_gallery.png` | `step3a_simulator.py` | 12 random TRAIN v1 tails (seed 5): real (solid) vs simulated m1 (dashed) and m4 (dotted) with the same tail_id, A_start, length, sub-link noise and NaN minutes |
| `step3a_thalf_by_astart.png` | `step3a_simulator.py` | median t_half per A_start bin: real (solid, filled circles, 95 % CI bars from events) vs each law (own line style and marker). Flat = exponential / two layers, falling = nonlinear |
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
| `step3b4_quantization.md` | `step3b4_quantization.py` (product; overwritten each run) | table setting x replicate (step kept, W_c, BIC gain); the verdict; what is quantized in the real data and whether the simulator matches; what it means for real tails below 0.35 dB; unclear points |
| `step4_real_drying.md` | `step4_real_drying.py` (product; overwritten each run) | what was run and cut; W range, tails, events, links; menu fits (as run and like for like) and the literature check; free curve with slope and band, model-free rate table; jump 2; agreement; checks 1-5; the result in one sentence; plain-words reading; unclear points |
| `step3b3c_method.md` | `step3b3c_method.py` (product; overwritten each run) | what was run, cut and installed; part A (jump 2 alone), B (step), C (jump 3), D (bootstrap band), E (F4 check) with tables and verdicts; the Step 4 recipe; what can be said about the real drying law and for which W; unclear points |
| `step3b3b_confirm.md` | `step3b3b_confirm.py` (product; overwritten each run) | what was run and cut, judged range per dataset, curve table, shape table, jump 2 terms and coefficients, m6 and m4, T1 robustness, verdict and the statement Step 4 can make (W >= 0.35 dB), unclear points |
| `step3b3_curve_then_equation.md` | `step3b3_curve_then_equation.py` (product; overwritten each run) | what was run and cut, spline / penalty / supported range / starts, jump 1 (band coverage, shape), jump 2 (terms, coefficients), m4 warning signs, robustness (no spread, T1 noise), comparison with the menu fit, verdicts and what they mean for step 4, unclear points |
| `step3b2_sindy_fake.md` | `step3b2_sindy_fake.py` (product; overwritten each run) | software and data, library / settings / offset handling, selection frequency per truth and term, the rule in every configuration, post hoc threshold check, coefficient recovery, m4 and the spread mixture, offset naive vs handled vs oracle, agreement with the pooled fit, verdicts Q1-Q4 and what they mean for step 5, unclear points |
| `step3b1c_stress_test.md` | `step3b1c_stress_test.py` (product; overwritten each run) | what was run and cut, T1 noise check, T3 recalibration, win shares and the pre-set rule per cell, parameter recovery, T2 start sensitivity, overall verdict and what it means for step 4, unclear points |
| `step3b1b_pooled_fit.md` | `step3b1b_pooled_fit.py` (overwritten each run) | plain-words summary incl. what the real tails can tell, method (marginalisation, BIC's N), confusion per window, recovery, population statistics, the pre-set decision rule, unclear points |
| `step3b1_fitting_on_sim.md` | `step3b1_fitting_on_sim.py` (overwritten each run) | plain-words summary, the fit and its n_eff, confusion (BIC and AIC), recovery, population tests, ablation, unclear points |
| `step3a_fix.md` | `step3a_simulator.py` (overwritten each run) | 3a-fix: noise targets and parameters, spread scan per law, calibration and misfit, checks, t_half per A_start bin, unclear points (replaces `step3a_simulator.md`, retired) |
| `step2c_freeze.md` | `step2c_freeze.py` (overwritten each run) | changes by the last fixes, event table and look figures, counts of the frozen set (train / val, with events), tails per event, median curves, unclear points |
| `step2b_selection.md` | `step2b_fix_tails.py` (overwritten each run) | new start, two ends, tail time, noise exclusion, level flags, main sets, median curves, gallery, unclear points (TRAIN only) |
| `step2_selection.md` | `step2_select_tails.py` (overwritten each run) | exclusions, truncation, main-set counts, fast-drop test with bootstrap CIs, shape descriptors, median curves, gallery, unclear points (TRAIN only) |
| `step1_tails.md` | `step1_extract_tails.py` (overwritten each run) | counts, distributions, return to the dry level, median curves, dry-weather rises, resolution check, list of unclear points |

## .cache/ (git-ignored, safe to delete)

`step0_run.log`, `step1_run.log`, `step2_run.log`, `step2b_run.log`, `step2c_run.log`, `step3a_run.log`, `step3b1_run.log`, `fits3b1/` (cached per-dataset fits of step 3b-1; delete to refit), `pooled3b1b/` (cached simulations and pooled fits of step 3b-1b), `step3b1b_run.log`, `fits3b1c/` (cached simulations, fits and the run plan of step 3b-1c; delete to redo), `step3b1c_run.log`, `sindy3b2/` (cached simulations, pooled fits, SINDy results and the run plan of step 3b-2; `first_run/` = the results before the post hoc diagnostic was added), `step3b2_run.log`, `fits3b3/` (cached simulations, every stage of every free-curve fit and the run plan of step 3b-3; delete to redo), `step3b3_run.log`, `fits3b3b/` (cached simulations and every stage of every fit of step 3b-3b; delete to redo), `step3b3b_run.log`, `fits3b3c/` (every stage of every fit of step 3b-3c; delete to redo), `step3b3c_run.log`, `fits3b4/` (simulations and fits of step 3b-4; delete to redo), `step3b4_run.log`, `check_new.py` (organization check incl. untracked files): console output of the last runs. `fits4/` (step 4: the real TRAIN datasets as parquet, every menu / curve / bootstrap fit, the PySR answers and the run plan; delete to refit), `step4_run.log`, `step4_run_cached.log`.
