# Baseline method (v2) — FROZEN

**Status: frozen on 2026-10-07 (R1c), commit `ebf122b`.** No method
change after this point without a new version (v3) and a new validation. v1 stays runnable for comparison.

A(t) = TL(t) - A_baseline(t), with TL = TSL - RSL (5 min `tl_mean`). The baseline is the TL a sub-link would show
without rain. It is estimated per (rain event, sub-link) from dry time BEFORE the event only, so it cannot follow the
drying tail (WAA) after the rain. Code: `baseline.py` (`baseline_v2`, `flags_v2`, `usable_mask`, `dew_risk`,
`a_of_t`). Settings: section `baseline` of `research/config.yaml`; `baseline.version` switches between v2 (current)
and v1 (`python research/baseline_estm/run_sample.py v1`).

## Rules

1. **Usable dry slot.** A 5 min slot is usable if the per-link reference state is dry (radar along the path and
   gauges within 1 km show no rain in it or in the 30 min around it; the link's own signal is never used) and it is
   at least `skip_h` = 1 h after the end of the last wet slot of the link. Optional **humidity exclusion**: also drop
   slots in an hour at dew risk (below). **Not adopted** (R1b rule, see `validation/out/REPORT_v2.md`), so
   `baseline.v2.humidity_exclusion: false`.
2. **Baseline.** Flat over the event: the **median** of the usable 5 min `tl_mean` values in the **24 h** before the
   event start. If that window has less than `min_dry_h` = 1 h of usable slots, the **48 h** window is used.
3. **Missing.** If the 48 h window also has less than 1 h of usable slots, the pair gets no baseline (no guessing).
   In the train + val table this is 2.9 % of the pairs, almost all from the first two events of the season (no
   history before them).
4. **Flags** (per event and sub-link; the baseline is still given, later stages decide):

| flag | meaning |
|---|---|
| `excursion` | the baseline is more than `excursion_db` = 1 dB from the sub-link's summer dry median (as in v1) |
| `dew_risk_in_window` | more than `dew.window_share` = 25 % of the window's dry slots are in dew-risk hours |
| `dew_risk_during_event` | more than 25 % of the event's time is in dew-risk hours; the event's time runs from the event start to 1 h after the end of the link's last wet slot (the event end if the link had no wet slot); each 5 min slot carries its hour's `dew_risk` from `products/dew_risk_hourly.csv`. **Added in R1c.** Rain itself makes the air humid, so this flag is common (53 % of significant-event pairs): it marks where dew may add to the drying, it is not a rejection rule |
| `daily_cycle_large` | the sub-link's dry-weather daily cycle (R1 C3 amplitude) is above `daily_cycle_large_db` = 0.4 dB (amplitude known for good sub-links only; False for degraded ones) |
| `degraded` | the sub-link's quality status is "degraded" |
| `toggling` | the sub-link's TSL class is "toggling" |

`use_in_waa_main` = not (missing or excursion or daily_cycle_large or degraded or toggling): the pairs for the main WAA
analysis. `excursion` was added to this rule after R1c (the learners' earlier decision: skip pairs whose baseline was
taken during a dry jump, and test them later with the WAA model); the baseline itself is unchanged.

## Changes v1 -> v2

| item | v1 | v2 | why (R1 / R1b / R1c) |
|---|---|---|---|
| window | 6 h before (and 6 h after) the event | 24 h before; 48 h fallback | the 24 h median beat the 6 h window on most sub-links and for D >= 6 h (R1 C1, C6) |
| statistic | mean of 5 min `tl_mean` | median of 5 min `tl_mean` | median had the lower median \|e\| in every setting (R1 C6) |
| shape | flat (A), or line before -> after (B) if they differ > 0.5 dB | flat only | the after window carries little drying, but the B cases looked like level changes, not drying (R1 C2, C5) |
| humidity exclusion | none | optional, dew-risk hours dropped | tested in R1b: not adopted (no gain at humid starts) |
| new flags | excursion (before / after window) | excursion, dew_risk_in_window, daily_cycle_large, degraded, toggling | R1 C3 / C4: daily cycles and dew are the main dry-weather disturbances |
| flag dew_risk_during_event | none | added in R1c | dew during the event can add to the attenuation (R1b humid-start result) |
| skip after rain | 1 h | 1 h (unchanged) | the skip does not matter in dry weather (R1 C6) |
| missing rule | < 1 h of usable dry time in the 6 h before | < 1 h in 24 h and in 48 h before | longer window, fewer missing pairs |

## Known accuracy (R1b, event level, good train sub-links)

Measured on 61 560 fake events in certified dry weather (truth: A = 0); offset = mean of A(t) over one fake event of
duration D, i.e. the baseline error of that event.

| D | median offset | p10 / p90 | p90 \|offset\| |
|---|---|---|---|
| 1 h | -0.00 dB | -0.24 / +0.24 dB | 0.33 dB |
| 6 h | -0.00 dB | -0.21 / +0.22 dB | 0.31 dB |
| 24 h | +0.02 dB | -0.15 / +0.23 dB | 0.28 dB |

- **Humid starts** (RH >= 90 % at the fake-event start) have a larger error: p90 |offset| 0.80 (1 h), 0.41 (6 h),
  0.32 dB (24 h). Part of that error is probably dew on the antenna during the fake event, so the test's
  "truth = 0" does not hold there: those numbers overstate the baseline error and contain real (dew) attenuation.
- Sub-links with a large daily cycle (> 0.4 dB) have about twice the error (p90 |offset| ~0.6 dB at D <= 6 h).
- **Drying after rain is not absorbed by the baseline.** Level after the rain minus the v2 baseline, real events vs a
  dry-weather control (medians): k = 1 h +0.050 vs +0.010 dB; k = 6 h +0.025 vs +0.010 dB; k = 12 h +0.010 vs
  +0.010 dB. The small excess at k <= 6 h is the drying signal that remains in A(t), as intended.

## Products (`products/`)

**`baselines_v2.csv`** (`make_baseline_table.py`): one row per (event, sub-link) for TRAIN and VAL events
(significant or not) on train links, good and degraded sub-links (test events and test links are never included).
Columns: `event_id`, `split_event`, `significant`, `run_count` (wet runs of the event), `event_start`, `event_end`,
`link_has_wet_slot`, `last_wet_end`, `sublink_id`, `cml_id`, `status`, `baseline_db`, `window_h` (24, 48 or 0 =
missing), `n_usable_slots`, `missing`, `summer_dry_median_db`, the flags above with `dew_share_in_window` and
`dew_share_during_event`, `offset_std`, `offset_p90_abs` (from the uncertainty table) and `use_in_waa_main` (rule above).
**A(t) is not stored**: `baseline.a_of_t(cml_id, event_id, res="1min" | "5min", pad="6h")` returns A(t) for both
sub-links of the link through CMLData.

**`baseline_uncertainty.csv`** (`make_baseline_table.py`): one row per good or degraded train sub-link: the spread of
its v2 baseline error over its C1 fake events (20 per duration): `offset_median`, `offset_std` (standard deviation of
the offsets at D <= 6 h), `offset_p90_abs` (D <= 6 h), `offset_p90_abs_24h`, the numbers of fake events, and the same
columns with the suffix `_no_humid_start` (fake events starting at RH >= 90 % left out). **These are meant as
weights for later stages** (e.g. 1 / offset_std^2 when combining sub-links or events), **not as a pass / fail**: a
sub-link with a large spread is still usable, it just counts less.

**`dew_risk_hourly.csv`** (`make_dew_products.py`): one row per hour of JJA 2015, SMHI station GbgA (central
Gothenburg), UTC hour start. Columns: `time_utc`, `time_local` (UTC + 2 h), `temperature` (degC), `humidity` (RH %),
`wind_speed` (m/s), `dew_point` (Magnus formula, b = 17.62, c = 243.12 degC), `dew_point_depression` (T - Td),
`dew_risk` (RH >= `dew.rh_pct` = 90 % OR depression < `dew.dpd_c` = 2 degC), `sun_elevation_deg` (at mid-hour, NOAA
formula, Gothenburg 57.709 N 11.975 E), `night` (sun elevation < `dew.night_sun_elev_deg` = -0.833 deg, i.e. the sun
below the horizon). One station stands for the whole network.

**`dry_jumps.csv`** (`make_dew_products.py`): the R1 C4 table of dry-day jumps for all train sub-links (good and
degraded). A jump is a run of >= 3 h of dry slots (gaps <= 15 min bridged) where `tl_mean` is >= 1 dB above
(positive) or below (negative) the sub-link's summer dry median, slot by slot and in the run median. Columns:
`cml_id`, `sublink_id`, `group`, `is_degraded`, `is_toggling`, `sign`, `start_utc`, `end_utc`, `start_local`,
`end_local`, `duration_h`, `median_excess_db`, `peak_excess_db`, `mean_rh`, `rh_max`, `mean_dew_point_depression`,
`mean_temperature`, `mean_wind_speed`, `hours_since_wet`, `both_directions` (the other direction of the link has a
same-sign jump overlapping >= 50 % of this one), `neighbours_jumping` (number of other train links with a path
midpoint within 5 km and such an overlapping jump), `length_class`, `frequency_band`, `tsl_behaviour`. Thresholds:
`baseline_validation.c4`.

Build order: `validation/run_all.py` -> `make_dew_products.py` -> `make_baseline_table.py`.

## Where the numbers are

R1 (v1 and variants, slot level): `validation/out/REPORT.md`. R1b (v2, event level, adoption decision, PASS/FAIL,
real-event coverage): `validation/out/REPORT_v2.md`.

Frozen at commit: `ebf122b` (Research R1c: freeze baseline v2, baseline table and uncertainty products), 2026-10-07.
