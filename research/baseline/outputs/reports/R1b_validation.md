# R1b: baseline v2, event-level validation and dew flag

> Files: the tables named below are in [`../tables/`](../tables/), the figures in [`../figures/`](../figures/) (moved
> there from `validation/out/` in the reorganization of 2026-10-07; file meanings: [`../../CATALOG.md`](../../CATALOG.md)).

Train links and significant train events only; main numbers = 513 good train sub-links. Method: `../../METHOD.md`.
Code: `validation/r1b_event_level.py` (in `run_all.py`), `make_dew_products.py`; settings: `research/config.yaml`
(`baseline`, `baseline_validation.r1b`). Event-level error: **offset** = mean of e(t) = tl_mean - baseline over one
C1 fake event (certified dry weather, truth 0), on the same 61 560 fake events as R1 (20 per sub-link and D).
"v2" below = v2 as adopted (no humidity exclusion).

## Summary

| question | answer (numbers) |
|---|---|
| event-level error of v2 (offset) | D = 1 h: median -0.00, p10 / p90 -0.24 / +0.24, p90 \|offset\| 0.33 dB; 6 h: -0.00, -0.21 / +0.22, 0.31 dB; 24 h: +0.02, -0.15 / +0.23, 0.28 dB |
| v2 vs v1 | p90 \|offset\| v1 -> v2: 0.29 -> 0.33 (1 h), 0.31 -> 0.32 (3 h), 0.34 -> 0.31 (6 h), 0.36 -> 0.28 (12 h), 0.30 -> 0.28 (24 h), 0.38 -> 0.29 dB (40 h): v2 is worse for 1-3 h events, better from 6 h on |
| humidity exclusion | **not adopted**. Humid starts (RH >= 90 % at t0, 570 fake events): p90 \|offset\| 0.392 without vs 0.408 dB with (gain -0.016 dB; rule needs >= +0.05). Missing share 0.0 % in both (rise 0.0 pp; rule allows <= 2 pp) |
| PASS share (v2) | all three: **42.7 %** (219 / 513). \|median offset\| <= 0.1: 98.4 %; p90 \|offset\| <= 0.25 dB: 44.4 %; growth at 24 h <= 0.15 dB: 93.2 % |
| daily_cycle_large and failures | the 72 sub-links with a daily cycle > 0.4 dB all fail (0 % pass, all on p90 \|offset\|); they are 24.5 % of the 294 failures. Without them: 49.7 % pass |
| real-event coverage (good) | 9 234 pairs: baseline 94.4 %, missing 5.6 % (all = the season's first event E20150601_0000, no history), 48 h fallback 0 %. Flags: excursion 0.7 %, dew_risk_in_window 2.2 %, daily_cycle_large 14.0 % |
| Delta after the rain (v2 before) | k = 1 h: real +0.050 vs control +0.010 dB (p90 0.35 vs 0.36); k = 6 h: +0.025 vs +0.010 dB; k = 24 h: 0.000 vs +0.010 dB |

## 1. Method

`METHOD.md`: the v2 rules, flags, the "Changes v1 -> v2" table and the two dew products. v1 stays runnable:
`baseline.version` in the config, or `python research/baseline_estm/run_sample.py v1` (v2: `... v2`, writes
`sample_pre_event_baselines.csv` for the 24-sub-link sample: missing 5.6 %, excursion 1.2 %, dew_risk_in_window 1.9 %,
daily_cycle_large 4.2 %).

## 2. Baseline v2

Median of usable dry 5 min slots in the 24 h before the event start (48 h fallback, else missing); flags excursion,
dew_risk_in_window, daily_cycle_large, degraded, toggling (`baseline.py`: `baseline_v2`, `flags_v2`).

- The 48 h fallback never triggered on fake or real events: every event after the first had >= 1 h of usable dry
  time in the 24 h before it. So v2 is numerically identical to M1 (24 h median) here.
- `daily_cycle_large` uses the R1 C3 amplitude, which exists for good sub-links only; it is False for degraded ones.
- Missing pairs are only the season's first event (the data start at its start time).

## 3. Dew products (`../products/`)

`dew_risk_hourly.csv`: 2 208 hours; dew risk 13.4 % of hours (16.8 % at night vs 12.0 % by day); night 29.7 %.
`dry_jumps.csv`: 3 934 dry-day jumps (2 653 positive) on 382 train sub-links (good and degraded). Columns and
thresholds: `METHOD.md`.

- Dew risk is more common at night, as expected, but also occurs by day (after rain, humid days).
- One SMHI station stands for the whole network; local dew can differ between sites.
- The jump table is the R1 C4 table with local times added; its reading is in `R1_validation.md`, section C4.

## 4. Event-level validation (fake events)

`r1b_offset_by_variant_D.csv`, `r1b_splits.csv`, `r1b_fake_event_offsets.csv.gz`,
`r1b_humidity_exclusion_decision.csv`; figure `../figures/r1b_offset_by_D.png`.

| variant | D | median offset | p10 | p90 | p90 \|offset\| | missing | slot p10 / p90 (reference) |
|---|---|---|---|---|---|---|---|
| v1 | 1 | -0.004 | -0.21 | +0.20 | 0.29 | 0.2 % | -0.23 / +0.22 |
| v1 | 6 | -0.005 | -0.25 | +0.23 | 0.34 | 0.3 % | -0.30 / +0.27 |
| v1 | 24 | +0.014 | -0.19 | +0.23 | 0.30 | 0 % | -0.31 / +0.32 |
| v1 | 40 | +0.034 | -0.25 | +0.28 | 0.38 | 0 % | -0.34 / +0.37 |
| M0 | 1 | -0.011 | -0.31 | +0.28 | 0.43 | 0 % | -0.32 / +0.30 |
| M0 | 6 | -0.017 | -0.29 | +0.25 | 0.40 | 0 % | -0.32 / +0.30 |
| M0 | 24 | -0.008 | -0.23 | +0.23 | 0.34 | 0 % | -0.32 / +0.31 |
| M1 = v2 | 1 | -0.001 | -0.24 | +0.24 | 0.33 | 0 % | -0.27 / +0.27 |
| M1 = v2 | 3 | -0.003 | -0.22 | +0.22 | 0.32 | 0 % | -0.27 / +0.27 |
| M1 = v2 | 6 | -0.002 | -0.21 | +0.22 | 0.31 | 0 % | -0.28 / +0.28 |
| M1 = v2 | 12 | -0.001 | -0.19 | +0.20 | 0.28 | 0 % | -0.28 / +0.29 |
| M1 = v2 | 24 | +0.024 | -0.15 | +0.23 | 0.28 | 0 % | -0.28 / +0.31 |
| M1 = v2 | 40 | +0.047 | -0.15 | +0.24 | 0.29 | 0 % | -0.29 / +0.33 |
| v2 with exclusion | 1 / 6 / 24 | 0.000 / 0.000 / +0.030 | | | 0.33 / 0.31 / 0.28 | 0 % | |

v2 split by humid start, daily cycle and class (p90 |offset|, dB):

| split | 1 h | 6 h | 24 h |
|---|---|---|---|
| humid start (RH >= 90 % at t0; n = 134 / 106 / 91) | 0.80 | 0.41 | 0.32 |
| not humid | 0.33 | 0.31 | 0.28 |
| daily_cycle_large | 0.60 | 0.59 | 0.36 |
| no large daily cycle | 0.29 | 0.27 | 0.26 |
| 1-3 km, >= 25 GHz | 0.31 | 0.29 | 0.24 |
| < 1 km, >= 25 GHz | 0.33 | 0.30 | 0.21 |
| > 3 km, < 25 GHz | 0.45 | 0.45 | 0.46 |
| > 3 km, >= 25 GHz | 0.33 | 0.30 | 0.32 |

Adoption rule (fixed before the run): adopt the exclusion if p90 |offset| at humid starts improves by >= 0.05 dB and
the overall missing share rises by <= 2 pp. Result: 0.392 -> 0.408 dB (worse by 0.016) and +0.0 pp -> **not adopted**.

- Event-level errors are unbiased (median within +-0.05 dB) and about +-0.2 dB (p10 / p90); averaging over an event
  removes part of the slot-level quantization noise.
- v2 wins for events of 6 h and longer (the real significant events last several hours), v1 for 1-3 h events.
- Humid starts and large daily cycles double the error; excluding dew-risk hours from the window does not help.
  A likely reading (not tested here): the error comes from what happens during the event (e.g. dew forming), which
  no choice of the before window can see.

## 5. PASS / FAIL (v2 as adopted)

`r1b_pass_fail_per_sublink.csv` (with `fail_reason`), `r1b_pass_fail_summary.csv`.

| group | n | \|median offset\| <= 0.1 | p90 \|offset\| <= 0.25 | growth <= 0.15 | all |
|---|---|---|---|---|---|
| all good | 513 | 98.4 % | 44.4 % | 93.2 % | **42.7 %** |
| 1-3 km, >= 25 GHz | 233 | 97.9 % | 53.6 % | 98.3 % | 52.4 % |
| < 1 km, >= 25 GHz | 90 | 97.8 % | 44.4 % | 96.7 % | 43.3 % |
| > 3 km, < 25 GHz | 48 | 100 % | 12.5 % | 83.3 % | 8.3 % |
| > 3 km, >= 25 GHz | 142 | 99.3 % | 40.1 % | 85.9 % | 38.0 % |
| TSL constant | 388 | 98.5 % | 44.6 % | 95.1 % | 43.0 % |
| TSL rare_steps | 113 | 98.2 % | 46.9 % | 87.6 % | 44.2 % |
| TSL frequent_steps | 12 | 100 % | 16.7 % | 83.3 % | 16.7 % |
| daily_cycle_large | 72 | 88.9 % | 0 % | 94.4 % | 0 % |
| no large daily cycle | 441 | 100 % | 51.7 % | 93.0 % | 49.7 % |

Fail reasons (294 failing sub-links; several possible): p90 |offset| > 0.25 dB: 285; growth at 24 h > 0.15 dB: 35;
|median offset| > 0.1 dB: 8.

- Almost all failures are the spread criterion; the median sub-link has p90 |offset| = 0.26 dB at D <= 6 h, just above
  the 0.25 dB limit, so the pass share is sensitive to that limit.
- daily_cycle_large explains a quarter of the failures (72 of 294) and those sub-links never pass.
- Long low-frequency links (> 3 km, < 25 GHz) are the weakest class (8 % pass).

## 6. Real events

`r1b_real_event_coverage.csv`, `r1b_real_events.csv.gz`, `r1b_c2_delta.csv`; figure `../figures/r1b_c2_delta.png`.

| group | pairs | baseline | missing | excursion | dew_risk_in_window | daily_cycle_large | degraded | toggling |
|---|---|---|---|---|---|---|---|---|
| good | 9 234 | 94.4 % | 5.6 % | 0.7 % | 2.2 % | 14.0 % | 0 % | 0 % |
| degraded (separate) | 954 | 93.9 % | 6.1 % | 19.4 % | 1.8 % | n.a. | 100 % | 30.2 % |
| toggling (separate) | 288 | 94.4 % | 5.6 % | 14.9 % | 1.4 % | n.a. | 100 % | 100 % |

| k (h) | real n | real median / p10 / p90 | control median / p10 / p90 |
|---|---|---|---|
| 1 | 7107 | +0.050 / -0.20 / +0.35 | +0.010 / -0.28 / +0.36 |
| 2 | 7435 | +0.040 / -0.22 / +0.34 | +0.010 / -0.28 / +0.36 |
| 3 | 7100 | +0.030 / -0.24 / +0.34 | +0.010 / -0.28 / +0.37 |
| 6 | 6333 | +0.025 / -0.26 / +0.35 | +0.010 / -0.29 / +0.37 |
| 12 | 5029 | +0.010 / -0.29 / +0.31 | +0.010 / -0.29 / +0.37 |
| 24 | 5588 | 0.000 / -0.28 / +0.30 | +0.010 / -0.30 / +0.39 |

- With the v2 baseline the after-rain level is +0.04 dB above the control at k = 1 h, shrinking to 0 by k = 12 h: a
  small drying tail is visible in the median and the v2 baseline does not absorb it.
- Flags are rare on good sub-links except daily_cycle_large (14 %); excursions are a degraded-link issue (19 %).
- The only missing baselines are the first event of the season.

## Runtimes (`validation_runtimes.csv`)

Full validation (C1-C6 + R1b): 209 s, of which R1b 61 s; dew products and the v2 sample: a few seconds each.
