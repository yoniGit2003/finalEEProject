# DECISIONS: baseline_estm

What each step decided, with the numbers behind it, so that the step reports can be deleted by /tidy and rebuilt on
demand (rerun the script named in each heading). The method itself: METHOD.md. The two validation reports were
written by hand and are key outputs: `R1_validation.md` (R1) and `R1b_validation.md` (R1b).

## Method v1 sample run (`run_sample.py v1`)

- Method v1 ("two-window"): mean of usable dry TL in 6 h windows before and after the event; flat (A) or a line (B)
  if they differ by > 0.5 dB. Sample: 24 sub-links x 18 significant train events = 432 pairs.
- Coverage: A 380, B 24, missing 28 (no usable before window). Precision: std of dry TL in the before window median
  0.09 dB (p90 0.22), below one RSL step (0.3 dB).
- Closure |after - before|: median 0.15 dB, p90 0.46 dB, 8 % above 0.5 dB; worst for > 3 km, < 25 GHz links (0.24).
- Out-of-sample check (mean A on dry bins inside the event, truth 0): median 0.05, p10 / p90 -0.18 / 0.42 dB.
- The skip X (0.5-3 h after the last wet bin) hardly matters (p90 change ≤ 0.03 dB).
- Decision: v1 stays runnable for comparison only; its line (B) cases do not look like drying (R1, C2 / C5).

## R1, slot-level validation (`validation/run_all.py`, checks C1-C6; report `R1_validation.md`)

- On fake events in certified dry weather, the 6 h median before the event (M2b) is not better than a 24 h median
  (M1) on most sub-links (beats both M0 and M1 on 29 %); M1 is best from D = 6 h on.
- Bias is not a problem (98 % pass |median e| ≤ 0.1 dB); the ±0.3 dB spread criterion equals one RSL step.
- Daily cycle of dry TL: median 0.15 dB, but 14 % of sub-links > 0.4 dB -> the `daily_cycle_large` flag.
- Positive dry jumps happen at night in humid air (dew?); none of the 13 flagged v1 cases looks like drying.
- Decision: method v2 = median of usable dry slots in the 24 h before the event (48 h fallback), flat, with flags.

## R1b, event-level validation of v2 (`validation/r1b_event_level.py`; report `R1b_validation.md`)

- v2 vs v1, p90 |offset|: worse for 1-3 h events (0.33 vs 0.29 dB at 1 h), better from 6 h on (0.28 vs 0.36 at 12 h).
- Humidity exclusion **not adopted**: p90 |offset| at humid starts 0.392 without vs 0.408 dB with (the rule, fixed
  before the run, needed a gain ≥ 0.05 dB).
- PASS (all three criteria): 42.7 % of 513 good sub-links; the 72 `daily_cycle_large` sub-links never pass.
- After real rain the level is +0.04 dB above the dry control at 1 h, 0 by 12 h: v2 does not absorb the drying tail.
- Decision: **v2 frozen on 2026-10-07** (commit `ebf122b`); products built by `make_dew_products.py` and
  `make_baseline_table.py`.
