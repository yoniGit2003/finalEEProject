# Baseline method (v2)

A(t) = TL(t) - A_baseline(t), with TL = TSL - RSL (5 min `tl_mean`). The baseline is the TL a sub-link would show
without rain. It is estimated per (rain event, sub-link) from dry time BEFORE the event only, so it cannot follow the
drying tail (WAA) after the rain. Code: `baseline.py` (`baseline_v2`, `flags_v2`, `usable_mask`, `dew_risk`).
Settings: section `baseline` of `research/config.yaml`; `baseline.version` switches between v2 (current) and v1
(kept runnable for comparison: `python research/baseline_estm/run_sample.py v1`).

## Rules

1. **Usable dry slot.** A 5 min slot is usable if the per-link reference state is dry (radar along the path and
   gauges within 1 km show no rain in it or in the 30 min around it; the link's own signal is never used) and it is
   at least `skip_h` = 1 h after the end of the last wet slot of the link. Optional **humidity exclusion**: also drop
   slots in an hour at dew risk (below). **Not adopted** (R1b rule, see `validation/out/REPORT_v2.md`), so
   `baseline.v2.humidity_exclusion: false`.
2. **Baseline.** Flat over the event: the **median** of the usable 5 min `tl_mean` values in the **24 h** before the
   event start. If that window has less than `min_dry_h` = 1 h of usable slots, the **48 h** window is used.
3. **Missing.** If the 48 h window also has less than 1 h of usable slots, the pair gets no baseline (no guessing).
   On the significant train events this happens only for the first event of the season (no history).
4. **Flags** (per event and sub-link; the baseline is still given, later stages decide):

| flag | meaning |
|---|---|
| `excursion` | the baseline is more than `excursion_db` = 1 dB from the sub-link's summer dry median (as in v1) |
| `dew_risk_in_window` | more than `dew.window_share` = 25 % of the window's dry slots are in dew-risk hours |
| `daily_cycle_large` | the sub-link's dry-weather daily cycle (R1 C3 amplitude) is above `daily_cycle_large_db` = 0.4 dB |
| `degraded` | the sub-link's quality status is "degraded" |
| `toggling` | the sub-link's TSL class is "toggling" |

## Changes v1 -> v2

| item | v1 | v2 | why (R1 / R1b) |
|---|---|---|---|
| window | 6 h before (and 6 h after) the event | 24 h before; 48 h fallback | the 24 h median beat the 6 h window on most sub-links and for D >= 6 h (R1 C1, C6) |
| statistic | mean of 5 min `tl_mean` | median of 5 min `tl_mean` | median had the lower median \|e\| in every setting (R1 C6) |
| shape | flat (A), or line before -> after (B) if they differ > 0.5 dB | flat only | the after window carries little drying, but the B cases looked like level changes, not drying (R1 C2, C5) |
| humidity exclusion | none | optional, dew-risk hours dropped | tested in R1b: not adopted (no gain at humid starts) |
| new flags | excursion (before / after window) | excursion, dew_risk_in_window, daily_cycle_large, degraded, toggling | R1 C3 / C4: daily cycles and dew are the main dry-weather disturbances |
| skip after rain | 1 h | 1 h (unchanged) | the skip does not matter in dry weather (R1 C6) |
| missing rule | < 1 h of usable dry time in the 6 h before | < 1 h in 24 h and in 48 h before | longer window, fewer missing pairs |

## Dew products (`products/`, built by `make_dew_products.py`)

Used later by other research code (WAA stage: dry-weather wet antennas).

**`dew_risk_hourly.csv`**: one row per hour of JJA 2015, SMHI station GbgA (central Gothenburg), UTC hour start.
Columns: `time_utc`, `time_local` (UTC + 2 h), `temperature` (degC), `humidity` (RH %), `wind_speed` (m/s),
`dew_point` (Magnus formula, b = 17.62, c = 243.12 degC), `dew_point_depression` (T - Td), `dew_risk`
(RH >= `dew.rh_pct` = 90 % OR depression < `dew.dpd_c` = 2 degC), `sun_elevation_deg` (at mid-hour, NOAA formula,
Gothenburg 57.709 N 11.975 E), `night` (sun elevation < `dew.night_sun_elev_deg` = -0.833 deg, i.e. the sun below
the horizon). One station stands for the whole network.

**`dry_jumps.csv`**: the R1 C4 table of dry-day jumps for all train sub-links (good and degraded). A jump is a run of
>= 3 h of dry slots (gaps <= 15 min bridged) where `tl_mean` is >= 1 dB above (positive) or below (negative) the
sub-link's summer dry median, slot by slot and in the run median. Columns: `cml_id`, `sublink_id`, `group`,
`is_degraded`, `is_toggling`, `sign`, `start_utc`, `end_utc`, `start_local`, `end_local`, `duration_h`,
`median_excess_db`, `peak_excess_db`, `mean_rh`, `rh_max`, `mean_dew_point_depression`, `mean_temperature`,
`mean_wind_speed`, `hours_since_wet`, `both_directions` (the other direction of the link has a same-sign jump
overlapping >= 50 % of this one), `neighbours_jumping` (number of other train links with a path midpoint within 5 km
and such an overlapping jump), `length_class`, `frequency_band`, `tsl_behaviour`. Thresholds: `baseline_validation.c4`.

## Where the numbers are

R1 (v1 and variants, slot level): `validation/out/REPORT.md`. R1b (v2, event level, adoption decision, PASS/FAIL,
real-event coverage): `validation/out/REPORT_v2.md`.
