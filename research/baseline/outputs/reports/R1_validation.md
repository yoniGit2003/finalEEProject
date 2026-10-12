# R1: validation of the baseline method

> Files: the tables named below are in [`../tables/`](../tables/), the figures in [`../figures/`](../figures/) (moved
> there from `validation/out/` in the reorganization of 2026-10-07; what each check decided: [`../../METHOD.md`](../../METHOD.md)).

Release `openmrg_v1`; train links and significant train events only; main numbers = 513 good train sub-links
(none of them TSL "toggling"). Degraded (53) and toggling (16, all also degraded) are reported separately.
5 min `tl_mean`, seed 42, skip S = 1 h unless stated. Settings: `research/config.yaml` (`baseline_validation`).
Code: `research/baseline/validation/` (`run_all.py`). Runtime of the full run: 148 s (table at the end).
Errors are pooled over all 5 min slots of all fake events unless stated ("per sub-link" = per sub-link first).

## Summary

| question | answer (numbers) |
|---|---|
| typical baseline error at D <= 6 h (M2b) | median 0.00 dB; p10 / p90 = -0.25 / +0.26 (D = 1 h), -0.28 / +0.28 (3 h), -0.30 / +0.30 (6 h); median \|e\| 0.09-0.13 dB |
| error growth with D | M2b median \|e\| 0.09 (1 h) -> 0.13 (6 h) -> 0.14 (24 h) -> 0.17 dB (40 h); p90 +0.26 -> +0.30 -> +0.33 -> +0.39 dB. Per sub-link, extra p90 \|e\| at 24 h: median +0.02 dB |
| best variant and by how much | D = 1 h: M2b best on median \|e\| (0.090 vs M1 0.105, M0 0.150); D >= 6 h: M1 (24 h median) best (0.115 vs M2b 0.130 at 6 h; 0.130 vs 0.140 at 24 h). Per sub-link at D <= 6 h, M2b beats M1 in only 33 % and M0 in 60 %. Window length matters more than mean/median or skip (C6) |
| drying signal in the after window? | barely: real vs control median Delta = +0.030 vs +0.010 dB at k = 1 h, +0.015 vs +0.015 dB at k = 6 h; p90 equal (0.40 vs 0.40). M3 line (B) cases do not look like drying (Delta rises from +0.23 at k = 1 h to +0.50 at 6 h) |
| daily cycle amplitude | p50 0.15 dB, p90 0.50 dB; 37 % of sub-links > 0.2 dB, 14 % > 0.4 dB. Matters for that minority (same size as the error budget); negligible for the typical sub-link |
| what the dry jumps look like | positive jumps (1532 on 329 / 513 sub-links): 66 % both directions, 63 % with >= 1 neighbour, 32 % at RH >= 90 % (vs 5 % of all dry slots), dew-point depression median 2.8 vs 5.4 degC, most start 22-01 UTC (00-03 local). Negative jumps (409): daytime (mode 10 UTC), dry air (RH >= 90 %: 1.5 %) |
| share of good train sub-links passing | all four criteria: 17 % (89 / 513). Per criterion: \|median e\| 98 %, p10/p90 within +-0.3 dB 57 %, growth at 24 h 89 %, M2b beats M0 and M1 29 % |
| anything surprising | (1) the +-0.3 dB p10/p90 criterion equals one RSL quantization step, so slot noise alone fails it on many sub-links; (2) the 24 h rolling median (M1) beats the 6 h window from D = 6 h on; (3) some sub-links dip ~2 dB every dry afternoon (e.g. 10011_1), which produced 4 of the 13 flagged cases; (4) none of the 13 flagged cases looks like drying |

## C1: fake events in dry weather (truth A = 0)

Fake events of D = 1-40 h placed at random in certified dry periods (every slot dry for the link); each variant
computed as if it were an event; e = tl_mean - baseline. 20 fake events per sub-link and D could always be placed
(513 x 6 x 20 = 61 560 for main; `c1_event_counts.csv`). Missing baselines: <= 0.3 % for every variant.

a) per variant x D, main (`c1a_variants_by_D_main.csv`; degraded / toggling: `c1a_variants_by_D_degraded.csv`,
`c1a_variants_by_D_toggling.csv`):

| variant | D | median e | p10 | p90 | median \|e\| | RMS |
|---|---|---|---|---|---|---|
| M0 | 1 | -0.01 | -0.32 | +0.30 | 0.150 | 0.32 |
| M1 | 1 | 0.00 | -0.27 | +0.27 | 0.105 | 0.27 |
| M2 | 1 | 0.00 | -0.24 | +0.23 | 0.101 | 0.26 |
| M2b | 1 | 0.00 | -0.25 | +0.26 | 0.090 | 0.27 |
| M3 | 1 | 0.00 | -0.23 | +0.22 | 0.099 | 0.25 |
| M0 | 6 | -0.01 | -0.32 | +0.30 | 0.150 | 0.34 |
| M1 | 6 | 0.00 | -0.28 | +0.28 | 0.115 | 0.30 |
| M2 | 6 | 0.00 | -0.30 | +0.29 | 0.132 | 0.34 |
| M2b | 6 | 0.00 | -0.30 | +0.30 | 0.130 | 0.35 |
| M3 | 6 | -0.01 | -0.30 | +0.27 | 0.129 | 0.32 |
| M0 | 24 | -0.01 | -0.32 | +0.31 | 0.150 | 0.36 |
| M1 | 24 | 0.00 | -0.28 | +0.31 | 0.130 | 0.34 |
| M2 | 24 | +0.01 | -0.30 | +0.33 | 0.142 | 0.36 |
| M2b | 24 | 0.00 | -0.30 | +0.33 | 0.140 | 0.37 |
| M3 | 24 | 0.00 | -0.31 | +0.32 | 0.143 | 0.36 |
| M2b | 40 | +0.01 | -0.31 | +0.39 | 0.170 | 0.41 |

Separate groups, M2b: degraded p10 / p90 = -0.55 / +0.64 (1 h), -0.97 / +0.99 dB (24 h), RMS 1.2-1.6 dB;
toggling -0.48 / +0.60 (1 h), -0.94 / +0.85 dB (24 h).

b) M2b by class (`c1b_M2b_by_class.csv`, `c1b_M2b_by_tsl.csv`), p10 / p90 at D = 6 h: 1-3 km >= 25 GHz
-0.30 / +0.29; < 1 km -0.30 / +0.30; > 3 km >= 25 GHz -0.30 / +0.30; **> 3 km < 25 GHz -0.37 / +0.45** (at 24 h
-0.40 / +0.58). TSL constant -0.30 / +0.30; rare_steps -0.30 / +0.30; frequent_steps (12 sub-links) -0.43 / +0.43.

c) e against time since t0 (`c1c_error_vs_time.csv`, figure `../figures/c1c_error_vs_time.png`), M2b p10 / p90:
0-1 h -0.25 / +0.27; 3-6 h -0.31 / +0.32; 6-12 h -0.33 / +0.35; 24-40 h -0.33 / +0.47 (median +0.03).

d) per sub-link M2b: `c1d_M2b_per_sublink.csv` (sub-link p90 at D <= 6 h: median +0.27, 90th pct +0.59 dB).

e) against hour of day and humidity at t0 (`c1e_error_vs_hour_humidity.csv`, `../figures/c1e_error_vs_hour_humidity.png`):
M2b median per-event e within -0.03 ... +0.01 dB at every hour; widest spread at 05-08 UTC (p10 -0.38). At RH >= 90 %
(323 events) M2b p10 / p90 = -0.51 / +0.37 dB vs about -0.24 / +0.24 otherwise; M0 is biased +0.09 dB there.

- In dry weather the baseline has no bias (median 0.00 dB); its spread is about +-0.25-0.30 dB, i.e. one RSL step.
- The error grows slowly with event length (+0.04 dB in median |e| from 6 h to 40 h), mostly in the upper tail.
- Long, low-frequency links (> 3 km, < 25 GHz), frequent TSL steps, degraded sub-links and humid starts are worse.

## C2: does the after window carry a drying signal?

Real: 8 554 event x sub-link cases (488 without a wet slot of the link, 192 without a before window). Delta_k =
median of usable dry slots in [end + k, end + k + 1 h) minus the M2b before. Control: C1 fake events matched by
duration (`c2_delta_vs_k.csv`, `c2_real_cases.csv`, `c2_case_status.csv`, figure `../figures/c2_delta_vs_k.png`).

| k (h) | real n | real median | real p10 / p90 | control median | control p10 / p90 |
|---|---|---|---|---|---|
| 1 | 7038 | +0.030 | -0.25 / +0.40 | +0.010 | -0.30 / +0.40 |
| 2 | 7366 | +0.030 | -0.26 / +0.39 | +0.013 | -0.30 / +0.40 |
| 3 | 7031 | +0.020 | -0.29 / +0.39 | +0.013 | -0.30 / +0.39 |
| 6 | 6270 | +0.015 | -0.28 / +0.38 | +0.015 | -0.30 / +0.38 |
| 12 | 4960 | 0.000 | -0.32 / +0.36 | +0.010 | -0.30 / +0.42 |
| 24 | 5525 | 0.000 | -0.32 / +0.33 | +0.015 | -0.31 / +0.44 |

Per class: > 3 km < 25 GHz real +0.083 vs control +0.048 at k = 1 h; > 3 km >= 25 GHz real +0.040 vs +0.053.
M3 on real events: A 94.2 %, B 5.8 % (`c2_m3_method_share.csv`). B cases (`c2_m3_line_cases_delta_vs_k.csv`):
median Delta +0.23 (k = 1 h), +0.37 (2 h), +0.45 (3 h), +0.51 (6 h), +0.13 (12 h), +0.12 dB (24 h).

- The real curve is above the control by only ~0.02 dB at k = 1-3 h and meets it at k = 6 h: a drying signal one hour
  after the rain is at most a few hundredths of a dB in the median.
- The real lower tail is slightly narrower than the control (p10 -0.25 vs -0.30), i.e. a small shift upwards.
- Where M3 draws a line, Delta grows for 6 h before falling: that is not a drying tail (which would start high and
  shrink); it looks like a level change after the rain, so the line is not following drying.

## C3: daily cycle in dry weather

Anomaly = usable dry slot minus the centred 24 h median of usable dry slots; median per UTC hour; amplitude = max - min
(`c3_daily_cycle_per_sublink.csv`, `c3_summary.csv`, `c3_amplitude_by_class.csv`, figure `../figures/c3_daily_cycle.png`).

| | p10 | p50 | p90 | > 0.2 dB | > 0.4 dB |
|---|---|---|---|---|---|
| amplitude (dB), 513 sub-links | 0.05 | 0.15 | 0.50 | 37 % | 14 % |

Correlation of the anomaly with temperature: median -0.12 (p10 -0.37, p90 +0.19); with RH: median +0.15 (p10 -0.24,
p90 +0.44). Largest amplitudes: 10019_1 2.40, 10019_2 2.24, 10252_2 2.02 dB, all a daytime swing (low in the morning,
07-10 UTC, high in the afternoon, 15-17 UTC).

- For the typical sub-link the daily cycle (0.15 dB) is smaller than the quantization step and can be ignored.
- About 1 in 7 sub-links has a cycle above 0.4 dB, which is larger than the baseline error budget.
- The cycle tends to go with humidity (+) and against temperature (-), but weakly and with both signs across sub-links.

## C4: what are the dry-day jumps?

Runs of >= 3 h of dry slots >= 1 dB above (or below) the sub-link's summer dry median, 15 min gaps bridged
(`c4_jumps.csv`, `c4_summary.csv`, figures `../figures/c4_start_hour.png`, `../figures/c4_humidity_dewpoint.png`).

| main sub-links | positive | negative |
|---|---|---|
| jumps (sub-links with one) | 1532 (329) | 409 (64) |
| median duration / median excess | 4.2 h / +1.55 dB | 5.0 h / -1.24 dB |
| both directions | 66 % | 42 % |
| >= 1 neighbour within 5 km | 63 % | 46 % |
| mean RH >= 90 % (all dry slots: 5 %) | 32 % | 1.5 % |
| dew-point depression, median (all dry slots: 5.4 degC) | 2.8 degC | 7.0 degC |
| wind, median | 1.9 m/s | 3.0 m/s |
| typical start (UTC) | 22-01 (mode 00); 36 % in 00-05 (all dry slots: 25 %) | mode 10; 17 % in 00-05 |
| > 24 h after the last rain | 43 % | 55 % |

Degraded: positive 1121 jumps, 55 % both directions, RH >= 90 % only 8 %; toggling: 264, 54 %, 9 %.

Reading: the positive jumps look like dew or condensation rather than equipment. They start mostly around local
midnight, a third of them happen at RH >= 90 % (six times the share of all dry time), the dew-point depression is
small (median 2.8 degC, a sharp peak at 1 degC), the wind is light, two thirds affect both directions of the link and
most have a neighbour jumping at the same time. Equipment faults would hit one direction and show no weather link.
About a third of positive jumps (one direction, no neighbour, 16 % neither) remain unexplained. The negative jumps are
different: daytime, dry air, large dew-point depression, more often one direction; they look like a daytime dip of
the dry level (see C3 and C5), not like dew.

## C5: the 13 flagged cases

Delta_k as in C2 for the flagged sub-link and the other direction; RH / dew-point depression in the 6 h after the
rain; overlap with a C4 jump (`c5_flagged_cases.csv`, figures `../figures/c5_<sublink>_<event>.png`).

| event | sub-link | flag | Delta 1 / 6 / 24 h (self) | Delta 1 / 6 h (other) | RH after | C4 jump | hint |
|---|---|---|---|---|---|---|---|
| E20150705_1910 | 10011_1 | before | +1.33 / +1.33 / +1.02 | +0.05 / +0.25 | 63 | negative | shift-like |
| E20150804_1830 | 10011_1 | before | +1.31 / +0.73 / +1.20 | +0.05 / 0.00 | 75 | negative | shift-like |
| E20150810_1545 | 10011_1 | before | +0.68 / -0.01 / +0.43 | +0.29 / 0.00 | 79 | none | shift-like |
| E20150815_1305 | 10011_1 | before | +1.39 / +0.80 / +1.72 | +0.43 / +0.13 | 79 | negative | shift-like |
| E20150705_1910 | 10155_1 | before | -1.06 / -1.06 / n.a. | -0.14 / -0.13 | 64 | positive | shift-like |
| E20150617_0705 | 10281_1 | after | n.a. (next event) | n.a. | n.a. | positive | unclear |
| E20150618_0530 | 10281_1 | before | -1.06 / -1.22 / -0.82 | -0.64 / -0.57 | 86 | positive | shift-like |
| E20150707_1235 | 10281_1 | before + after | -0.09 / -0.62 / -0.59 | -0.23 / -0.41 | 65 | positive | shift-like |
| E20150617_0705 | 10281_2 | after | n.a. (next event) | n.a. | n.a. | none | unclear |
| E20150707_1235 | 10281_2 | before | -0.23 / -0.41 / -0.79 | -0.09 / -0.62 | 65 | positive | shift-like |
| E20150705_1910 | 10338_1 | after | -1.26 / -1.03 / -1.09 | -1.26 / -1.18 | 66 | negative | shift-like |
| E20150725_0440 | 10338_1 | after | -0.85 / n.a. / +0.12 | -0.87 / n.a. | 77 | none | shift-like |
| E20150725_0440 | 10338_2 | after | -0.87 / n.a. / -0.02 | -0.85 / n.a. | 77 | none | shift-like |

- 11 shift-like, 2 unclear (a new event starts before any dry time), 0 drying-like; 9 of 13 overlap a dry-day jump.
- The four 10011_1 cases are the same story: its before window fell in a daytime dip of ~2 dB (negative jump), so the
  before is too low and Delta stays positive for 24 h; the other direction does not dip.
- 10338 (both directions about 1 dB LOWER after the rain) and 10155_1 / 10281 (before window during a positive jump)
  are level changes, not drying; drying would make the after level higher and shrinking.

## C6: sensitivity of M2b (C1 fake events)

`c6_sensitivity.csv` (all 72 settings). Median |e| (dB), skip S = 1 h, median vs mean:

| window | D = 1 h | D = 6 h | D = 24 h |
|---|---|---|---|
| 3 h | 0.070 / 0.080 | 0.115 / 0.117 | 0.140 / 0.141 |
| 6 h | 0.090 / 0.101 | 0.130 / 0.132 | 0.140 / 0.142 |
| 12 h | 0.115 / 0.124 | 0.140 / 0.142 | 0.140 / 0.143 |
| 24 h | 0.105 / 0.121 | 0.115 / 0.129 | 0.130 / 0.143 |

Skip S = 3 or 6 h changes median |e| by <= 0.003 dB but raises the missing share (3 h window: 0.3 % -> 2.3 % -> 6.7 %).
Median e is 0.00 for every median setting; the mean settings have small negative medians (down to -0.016 dB).

- A shorter window is better for short events (3 h: 0.070 vs 6 h: 0.090 dB at D = 1 h); the 24 h window is best for
  long events (0.130 vs 0.140 at D = 24 h). No single window wins everywhere.
- Median beats mean on median |e| in every setting; RMS is similar (mean slightly lower when outliers are rare).
- The skip S does not matter in dry weather (it only costs coverage), so it can be chosen for the wet case alone.

## PASS / FAIL (M2b on C1, criteria fixed before the run)

`c1_pass_fail_per_sublink.csv` (with `fail_reason`), `c1_pass_fail_summary.csv`. Shares among evaluable sub-links (all 513
were evaluable):

| group | n | \|median e\| <= 0.1 | p10 >= -0.3 and p90 <= 0.3 | extra p90 at 24 h <= 0.2 | beats M0 and M1 | all four |
|---|---|---|---|---|---|---|
| all main | 513 | 98.4 % | 57.1 % | 89.3 % | 28.8 % | **17.3 %** |
| 1-3 km, >= 25 GHz | 233 | 97.9 % | 62.7 % | 94.8 % | 24.0 % | 16.3 % |
| < 1 km, >= 25 GHz | 90 | 100 % | 56.7 % | 92.2 % | 20.0 % | 8.9 % |
| > 3 km, < 25 GHz | 48 | 100 % | 22.9 % | 68.8 % | 35.4 % | 10.4 % |
| > 3 km, >= 25 GHz | 142 | 97.9 % | 59.9 % | 85.2 % | 40.1 % | 26.8 % |
| TSL constant | 388 | 99.2 % | 57.0 % | 89.4 % | 25.8 % | 14.4 % |
| TSL rare_steps | 113 | 97.3 % | 59.3 % | 88.5 % | 39.8 % | 28.3 % |
| TSL frequent_steps | 12 | 83.3 % | 41.7 % | 91.7 % | 25.0 % | 8.3 % |
| degraded (separate) | 53 | 88.7 % | 24.5 % | 77.4 % | 43.4 % | 9.4 % |
| toggling (separate) | 16 | 93.8 % | 0 % | 68.8 % | 50.0 % | 0 % |

Fail reasons (main, a sub-link can have several): not better than M0 and M1: 365 (46 of them exact ties, quantized
values); p10 / p90 outside +-0.3 dB: 220; extra error at 24 h: 55; |median e| > 0.1 dB: 8.

- Bias is not a problem (98 % pass); growth with event length is mostly fine (89 %).
- The spread criterion fails on 43 % of sub-links, but +-0.3 dB is one RSL quantization step, so a sub-link sitting
  between two levels fails it even with a perfect baseline.
- The main failure is that the 6 h median is not better than a 24 h median (M1) on most sub-links.

DEV sample (24 sub-links) gave the same picture; only "all four" differed (8 % vs 17 %).

## Runtimes (`validation_runtimes.csv`, full run with a cold data cache)

| step | C1 | C2 | C3 | C4 | C5 | C6 | total |
|---|---|---|---|---|---|---|---|
| seconds | 70.9 | 24.7 | 8.6 | 6.9 | 5.8 | 30.9 | 148.3 |
