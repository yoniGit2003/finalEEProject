# CATALOG: every output file of baseline_estm

One entry per output file (or group of similar files). Paths are relative to `research/baseline_estm/`. Commands run
from the repository root with `C:\Users\yoni\miniconda3\python.exe` (written `python` below). Folder locations are
set in `baseline_estm_paths.py`. Status of every entry: **current** (there is no `archive/` in this module yet).

## Shared vocabulary

| term | meaning |
|---|---|
| TL, `tl_mean` | total loss TSL - RSL (dB), 5 min mean from the release |
| baseline | the TL the sub-link would show without rain (dB) |
| A(t), e(t) | A = TL - baseline (dB). In dry weather the truth is A = 0, so e(t) = TL - baseline is the baseline error |
| usable dry slot | 5 min slot that is dry for the link (radar + gauges, never the CML) and >= `skip_h` = 1 h after its last wet slot |
| method v1, "two-window" | 6 h windows before and after the event, mean, flat (A) or a line (B). METHOD.md, last section |
| method v2, "pre-event" | median of usable slots in the 24 h before the event (48 h fallback), flat, with flags. The frozen method |
| variants (R1) | M0 = summer dry median of the sub-link (a constant); M1 = median of the 24 h before; M2 = mean of the 6 h before; M2b = median of the 6 h before; M3 = method v1 |
| fake event | a span of D hours (1, 3, 6, 12, 24, 40) in certified dry weather, treated as if it were a rain event (C1) |
| offset | mean of e(t) over one fake event = the baseline error of that event (R1b) |
| D, `D_h` | fake-event duration (h); "short" = D <= 6 h, "24h" / "long" = D = 24 h |
| k, `k_h`; `d1` ... `d24` | hours after the event end; `dk` = median of usable dry TL in [end + k, end + k + 1 h) minus the baseline (dB): what is left after the rain (drying) |
| groups | `main` = good train sub-links (513); `degraded`; `toggling` (TSL toggles; all also degraded) |
| statistic columns | `median`, `p10`, `p90` (percentiles), `median_abs` (median of \|e\|), `rms`, `n_*` (counts), `missing_share` (share of events without a baseline); all in dB unless a count or share |
| link columns | `cml_id`, `sublink` (1 / 2), `sublink_id` (`<cml_id>_<sublink>`), `status` (good / degraded), `tsl_behaviour`, `length_class`, `frequency_band`, `length_m`, `frequency_ghz`, `is_degraded`, `is_toggling` (from the release) |

## outputs/products/: the deliverables (read by later modules)

### `outputs/products/baselines.csv`
- **What**: the frozen v2 baseline, one row per (event, sub-link) for train and val events on train links, good and
  degraded sub-links (40 186 rows). Test events / links are never included.
- **Answers**: "what is the baseline of this sub-link for this event, and can I trust it?" Input of the WAA stage.
- **Columns**: `event_id`, `split_event` (train / val), `significant`, `run_count` (wet runs in the event),
  `event_start` / `event_end` (UTC), `link_has_wet_slot`, `last_wet_end` (end of the link's last wet slot),
  `sublink_id`, `cml_id`, `status`, `baseline_db` (dB; empty if missing), `window_h` (24, 48 or 0 = missing),
  `n_usable_slots`, `missing`, `summer_dry_median_db` (dB), flags `excursion`, `dew_risk_in_window`,
  `dew_risk_during_event`, `daily_cycle_large`, `degraded`, `toggling` with shares `dew_share_in_window`,
  `dew_share_during_event`; `offset_std`, `offset_p90_abs` (dB, from `baseline_uncertainty.csv`), `use_in_waa_main`
  (= not missing / excursion / daily_cycle_large / degraded / toggling). Full definitions: METHOD.md.
- **Producer**: `python research/baseline_estm/make_baseline_table.py`.
- **Inputs**: the release, `outputs/products/dew_risk_hourly.csv`, `outputs/tables/c1_fake_events.csv.gz`,
  `c3_daily_cycle_per_sublink.csv`, `r1b_fake_event_offsets.csv.gz` (consistency check).
- Before 2026-10-07 this file was `products/baselines_v2.csv`.

### `outputs/products/baseline_uncertainty.csv`
- **What**: per good or degraded train sub-link, the spread of its v2 baseline error over its C1 fake events.
- **Answers**: "how much should this sub-link count?" Weights for later stages (e.g. 1 / `offset_std`^2), not pass / fail.
- **Columns**: link columns; `n_fake_short`, `n_fake_24h`; `offset_median`, `offset_std`, `offset_p90_abs` (D <= 6 h),
  `offset_p90_abs_24h` (dB); the same with suffix `_no_humid_start` (fake events starting at RH >= 90 % left out).
- **Producer / inputs**: as `baselines.csv`.

### `outputs/products/dew_risk_hourly.csv`
- **What**: one row per hour of JJA 2015 at SMHI station GbgA (one station for the whole network).
- **Answers**: "could dew form on the antennas in this hour?"
- **Columns**: `time_utc` (hour start), `time_local` (UTC + 2 h), `temperature` (degC), `humidity` (RH %),
  `wind_speed` (m/s), `dew_point` (degC, Magnus), `dew_point_depression` (T - Td, degC), `dew_risk` (RH >= 90 % or
  depression < 2 degC), `sun_elevation_deg` (mid-hour), `night` (sun below the horizon).
- **Producer**: `python research/baseline_estm/make_dew_products.py`. **Inputs**: the release (weather).

### `outputs/products/dry_jumps.csv`
- **What**: the R1 C4 dry-day jumps for all train sub-links: runs of >= 3 h of dry slots >= 1 dB above (positive) or
  below (negative) the sub-link's summer dry median.
- **Answers**: "when does the dry baseline shift by itself (dew, multipath), and on how many links at once?"
- **Columns**: link columns, `group`, `sign`, `start_utc` / `end_utc` / `start_local` / `end_local`, `duration_h`,
  `median_excess_db`, `peak_excess_db`, `mean_rh`, `rh_max`, `mean_dew_point_depression`, `mean_temperature`,
  `mean_wind_speed`, `hours_since_wet`, `both_directions` (other direction jumps too), `neighbours_jumping` (other
  train links within 5 km jumping at the same time).
- **Producer**: `python research/baseline_estm/make_dew_products.py`. **Inputs**: `outputs/tables/c4_jumps.csv`.

## outputs/tables/: validation R1 (checks C1-C6)

All written by `python research/baseline_estm/validation/run_all.py` (or the single check script named in the
prefix, e.g. `validation/c1_fake_events.py`). Input: the release (5 min `tl_mean`, link state, weather). Read the
numbers in `outputs/reports/R1_validation.md`.

| file | producer | what it contains (columns) | question |
|---|---|---|---|
| `c1_fake_events.csv.gz` | C1 | one row per fake event: `event`, link ids, `D_h`, `t0` (start), `i0` (slot index), `n_candidates`, `m0`; per variant X: `base_X` (dB), `missing_X`, `mean_e_X` (mean error, dB); M3 details `m3_method`, `m3_flag`, `m3_before`, `m3_after` | the raw material of C1, C2, C6, R1b and the products |
| `c1_event_counts.csv` | C1 | per `group` and `D_h`: sub-links, events placed, sub-links with all / fewer / no events | could the fake events be placed? |
| `c1a_variants_by_D_{main,degraded,toggling}.csv` | C1 | per `variant` x `D_h`: error statistics | which variant has the smallest error, for each duration? |
| `c1b_M2b_by_class.csv`, `c1b_M2b_by_tsl.csv` | C1 | M2b error statistics per link class (`cls`) or `tsl_behaviour`, and `D_h` | which kinds of links are worse? |
| `c1c_error_vs_time.csv` | C1 | per `variant` and time-since-start bin (`bin_h`, `from_h`, `to_h`): error statistics | does the error grow during an event? |
| `c1d_M2b_per_sublink.csv` | C1 | per sub-link: `M*_short_median_abs`, short / long median, p10, p90, p90_abs, link columns | the per-sub-link error |
| `c1e_error_vs_hour_humidity.csv` | C1 | per `by` (hour or humidity) and `bin`: `n_events`, median error of M0 / M1 / M2b, M2b p10 / p90 | is the error larger at night or in humid air? |
| `c1_pass_fail_per_sublink.csv` | C1 | per sub-link: the C1 columns plus criteria `c1_median` (\|median e\| <= 0.1), `c2_spread` (p10 >= -0.3 and p90 <= 0.3), `c3_growth` (extra p90 at 24 h <= 0.2), `c4_beats` (M2b beats M0 and M1), `c4_tie`, `evaluable`, `pass_all`, `fail_reason` | does each sub-link pass the R1 criteria (fixed before the run)? |
| `c1_pass_fail_summary.csv` | C1 | per `group` / `split`: counts passing each criterion | how many pass? |
| `c2_real_cases.csv` | C2 | per real (event, sub-link): times, `case_status`, `before` (M2b, dB), M3 method / flag / `m3_after_minus_before`, `d1`...`d24`, link columns | is there drying in the window after real rain? |
| `c2_delta_vs_k.csv` | C2 | per `group` and `k_h`: `real_*` vs `ctrl_*` (fake-event control) n, median, p10, p90 of Delta_k | real events vs dry-weather control |
| `c2_case_status.csv`, `c2_m3_method_share.csv`, `c2_m3_line_cases_delta_vs_k.csv` | C2 | case counts; share of M3 methods A / B; Delta_k of the line (B) cases | do v1 line cases look like drying? |
| `c3_daily_cycle_per_sublink.csv` | C3 | per sub-link: `amplitude_db` (daily cycle of dry TL), `hour_max_utc`, `hour_min_utc`, `corr_temperature`, `corr_humidity`, `h00`...`h23` (hourly anomaly, dB), link columns | how big is each sub-link's daily cycle? (feeds the `daily_cycle_large` flag) |
| `c3_summary.csv`, `c3_amplitude_by_class.csv` | C3 | amplitude percentiles, shares above 0.2 / 0.4 dB, correlations; amplitude per link class | is the daily cycle a general problem? |
| `c4_jumps.csv` | C4 | one row per dry jump (see `dry_jumps.csv`; UTC / local hours, `rh_*`, `t_mean`, `dpd_mean`, `wind_mean`) | dry-weather baseline shifts |
| `c4_summary.csv` | C4 | per `group` and `sign`: counts, durations, shares in both directions / with neighbours / at RH >= 90 %, start hours | what do the jumps look like? |
| `c5_flagged_cases.csv` | C5 | the 13 flagged pairs of the v1 sample run: `flags`, before / after, M2b before and `d1`...`d24` for this (`self_`) and the other direction (`other_`), weather after the rain, `c4_jump_overlap`, `hint` | are the flagged cases drying or level changes? |
| `c6_sensitivity.csv` | C6 | M2b-type baseline per `window_h` x `skip_h` x `stat` (mean / median) x `D_h`: error statistics | does the window length / skip / statistic matter? |
| `validation_runtimes.csv` | `run_all.py` | `step`, `seconds` (one row per check, replaced on re-run) | how long does validation take? |

## outputs/tables/: validation R1b (method v2, event level)

Producer: `python research/baseline_estm/validation/r1b_event_level.py` (also in `run_all.py`). Inputs:
`c1_fake_events.csv.gz`, `c3_daily_cycle_per_sublink.csv`, the release. Numbers: `outputs/reports/R1b_validation.md`.

| file | what it contains (columns) | question |
|---|---|---|
| `r1b_fake_event_offsets.csv.gz` | per C1 fake event: `off_X` / `base_X` for X = v1, M0, M1, v2 with humidity exclusion (`v2_hx`), v2 without (`v2_nohx`); link columns, `daily_cycle_large`, `rh_t0`, `humid_start` | the event-level error of every method |
| `r1b_offset_by_variant_D.csv` | per `variant` x `D_h`: offset median, p10, p90, p90_abs, missing share, slot-level p10 / p90 | v2 vs v1 by duration |
| `r1b_splits.csv` | offset statistics split by `humid_start`, `split`, `daily_cycle_large`, class `cls` | where is v2 worse? |
| `r1b_humidity_exclusion_decision.csv` | p90 \|offset\| at humid starts with / without exclusion, gain, missing shares, the rule, `adopt` | is the humidity exclusion adopted? (no) |
| `r1b_pass_fail_per_sublink.csv`, `r1b_pass_fail_summary.csv` | per sub-link: `median_offset_short`, `p90_abs_short`, `p90_abs_24h`, criteria `c1_median`, `c2_p90`, `c3_growth`, `pass_all`, `fail_reason`; and counts per group / split | does v2 pass the R1b criteria? |
| `r1b_real_events.csv.gz` | per real (event, sub-link): v2 baseline, window, flags, `duration_h`, `d1`...`d24`, link columns. **Note:** the committed file lacks `dew_risk_during_event` / `dew_share_during_event` (written before R1c); a re-run adds them (open question in `docs/organization/REPORT.md`) | coverage and drying after real rain with v2 |
| `r1b_real_event_coverage.csv` | per group: pairs, with baseline, missing, 48 h fallback, flag counts | how many real pairs get a v2 baseline? |
| `r1b_c2_delta.csv` | as `c2_delta_vs_k.csv` but with the v2 baseline | does v2 absorb the drying? (no) |

## outputs/tables/: sample runs

| file | producer | what it contains | question |
|---|---|---|---|
| `sample_two_window_baselines.csv` | `python research/baseline_estm/run_sample.py v1` | method v1 on the sample (24 sub-links x 18 significant train events, X = 1 h): `before_db` / `after_db` (window means, dB), `before_hours` / `after_hours` (usable dry time), `*_std_db` (spread of 5 min TL in the window), `usual_level_db` (summer dry median), `closure_db` (after - before), `method` (A / B / missing), `flags` (`no_before`, `no_after`, `excursion_before`, `excursion_after`), `in_event_dry_hours`, `in_event_dry_A_mean_db` (out-of-sample check), `A_max_db` (largest A in the event), length / frequency class | how does v1 behave? **Read by C5.** |
| `sample_two_window_sensitivity_X.csv` | same | the same for X = 0.5, 1, 2, 3 h (column `X_h`) | does the skip X matter? |
| `sample_pre_event_baselines.csv` | `python research/baseline_estm/run_sample.py` (v2 is the default) | method v2 on the same sample: `baseline_db`, `window_h`, `missing`, `summer_dry_median_db`, the v2 flags and shares | quick look at v2 on the sample |

## outputs/figures/

| file | producer | what is plotted / what to look for |
|---|---|---|
| `sample_two_window_overview.png` | `run_sample.py v1` | the five v1 quality checks (coverage, precision, closure, out-of-sample A, sensitivity to X); see METHOD.md, v1 section |
| `sample_two_window_examples/<sublink>_<event>.png` | `run_sample.py v1` (redrawn each run) | top: TL (blue), usable dry samples (green dots), before / after windows (green bands), baseline (red), usual level (dotted), event start / end (dashed), background blue = wet, grey = unknown; bottom: A(t), which should return to ~0 after the rain and the tail |
| `c1c_error_vs_time.png` | C1 | baseline error vs time since the fake-event start, per variant: should stay flat |
| `c1e_error_vs_hour_humidity.png` | C1 | error vs hour of day and humidity at the start: humid / night starts are worse |
| `c2_delta_vs_k.png` | C2 | Delta_k after real rain vs the dry control: a gap = drying signal |
| `c3_daily_cycle.png` | C3 | daily cycle of dry TL for typical and largest sub-links |
| `c4_start_hour.png`, `c4_humidity_dewpoint.png` | C4 | when dry jumps start; humidity and dew-point depression during jumps (dew?) |
| `c5_<sublink>_<event>.png` | C5 | the 13 flagged v1 cases, with the weather after the rain: drying or a level change? |
| `r1b_offset_by_D.png` | R1b | v2 vs v1 / M0 / M1 offsets by duration |
| `r1b_c2_delta.png` | R1b | Delta_k after real rain with the v2 baseline vs control |

## outputs/reports/

| file | producer | content |
|---|---|---|
| `R1_validation.md` | written by hand from the R1 tables | validation of v1 and the variants at slot level, PASS / FAIL, conclusions (before 2026-10-07: `validation/out/REPORT.md`) |
| `R1b_validation.md` | written by hand from the R1b tables | v2 at event level, humidity-exclusion decision, real-event coverage (before: `validation/out/REPORT_v2.md`) |
| `sample_two_window_summary.md` | `run_sample.py v1` (overwritten each run) | the v1 sample quality numbers as tables |

## .cache/ (git-ignored, safe to delete)

`tl_5min.parquet` (loaded 5 min TL, by `validation/common.py`), `c1_slots.npz` (C1 per-slot errors), `run_r1b.log`,
`dev/` (outputs of `BV_SAMPLE=DEV` runs). Rebuilt automatically when missing (the first run is then slower).
