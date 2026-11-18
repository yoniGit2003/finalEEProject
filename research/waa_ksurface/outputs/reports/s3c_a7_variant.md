# Step 3c: fake runs of the A7 procedure (Deviation 14: learning events only)

Written by `research/waa_ksurface/s3c_a7_variant.py` (libraries `s3b_recovery.py`, `s3_recovery.py`, `a6_recovery.py`, `s1_separability.py`, `fake6.py`, `ladder.py`). Fake attenuation as in step 3b, fresh seeds (base 3800). **The fit uses the 14 learning events only** (S1 current variant): the tail flag before any fit, the link-events of flagged events left out; the 37 added events take no part. Drying K0 against K1 on L-dryhigh, margin 0.00291 dB (current variant, S1) with a paired 90 % interval that excludes zero; wetting K0, not judged (R3). Where two numbers are given as `x / y`, the first is with k_wet at its gauge-input estimate, the second at its radar-input estimate. A returned drying K1 is named: *non-exponential drying: the drying rate falls as the antenna dries* (R5). No real attenuation is read. VAL and TEST were not opened.

## (a) The link spread (truth K0, 1 / k_dry = 15 min)

Every fake truth of A6, S1, step 3 and step 3b already carried a level shift per link of sd 0.088 dB, the measured spread (E4), as part of the default noise; so `g_088` is truth (a) on a fresh seed under the A7 procedure, and every truth that passed so far passed at the measured link spread. It is run first because its STOP rule is new.

| run | link sd (dB) | seed | drying rung | K1 gain on L-dryhigh (dB) | interval | recovery (k_wet gauge / radar) | events flagged | scored values | verdict |
|---|---|---|---|---|---|---|---|---|---|
| g_088 | 0.088 | 3800 | K0 | 0.00016 | -0.00024 to +0.00145 | 1 / k_dry 16.7 / 16.6 min, truth 15 | 0 of 13 | 3092 | PASS |
| g_150 | 0.15 | 3801 | K0 | 0.00158 | +0.00020 to +0.00410 | 1 / k_dry 13.3 / 13.9 min, truth 15 | 0 of 13 | 2380 | PASS |


## (b) Detectability of k_dry(W) for the A7 variant

k_dry = (1/15)·(1 + a·W) / (1 + a), k_wet = 1/5; ratio k(W = 2 dB) / k(W = 0) = 1 + 2a. Five replicates per ratio, no early stop. Detected = drying K1 returned in at least 4 of 5.

| size_a | ratio | replicates_run | hits (K1 returned) | detected | runs |
|---|---|---|---|---|---|
| 0.25 | 1.5 | 5 | 2 | False | cw_0.25_r0, cw_0.25_r1, cw_0.25_r2, cw_0.25_r3, cw_0.25_r4 |
| 0.5 | 2.0 | 5 | 5 | True | cw_0.5_r0, cw_0.5_r1, cw_0.5_r2, cw_0.5_r3, cw_0.5_r4 |
| 1.0 | 3.0 | 5 | 5 | True | cw_1_r0, cw_1_r1, cw_1_r2, cw_1_r3, cw_1_r4 |

**Smallest ratio detected in at least 4 of 5: 2.**

Why step 3b showed 4 of 4 and not 5: there the replicates of a ratio stopped as soon as the verdict was decided (4 hits, or 2 misses; Deviation 13), so after 4 hits the fifth replicate was never run. Here every ratio runs all 5.


## (c) The flag on the 14 learning events

Checkable = at least 5 short-link link-events with the 20 bins after the end of rain (event timing only; no attenuation read). **13 of 14 learning events are checkable.** Removed on real data: not counted here (it needs the real attenuation at the end of the rain, which is not read before A7; it is the first number of A7).

| event_id | short_link_events | with_tail | checkable |
|---|---|---|---|
| E20150617_0705 | 55 | 55 | True |
| E20150618_0530 | 67 | 67 | True |
| E20150622_0930 | 67 | 67 | True |
| E20150705_1910 | 67 | 67 | True |
| E20150707_1235 | 67 | 67 | True |
| E20150725_0440 | 67 | 67 | True |
| E20150727_2140 | 67 | 67 | True |
| E20150804_1830 | 67 | 67 | True |
| E20150810_1545 | 67 | 67 | True |
| E20150815_1305 | 67 | 67 | True |
| E20150825_0145 | 67 | 67 | True |
| E20150826_0520 | 67 | 67 | True |
| E20150826_1620 | 67 | 67 | True |
| E20150831_0455 | 67 | 0 | False |


## Execution

Runs: 17; wall time 1.3 h, processor time 4.8 h; largest peak memory 2.94 GB.

| run | minutes | cpu minutes | peak memory (GB) |
|---|---|---|---|
| g_088 | 5.0 | 16.6 | 2.37 |
| g_150 | 4.9 | 17.7 | 2.93 |
| cw_0.25_r0 | 5.0 | 16.9 | 2.93 |
| cw_0.25_r1 | 3.6 | 13.3 | 2.93 |
| cw_0.25_r2 | 4.1 | 15.1 | 2.93 |
| cw_0.25_r3 | 4.2 | 15.4 | 2.93 |
| cw_0.25_r4 | 4.1 | 15.8 | 2.93 |
| cw_0.5_r0 | 4.2 | 16.4 | 2.92 |
| cw_0.5_r1 | 3.9 | 14.7 | 2.93 |
| cw_0.5_r2 | 5.7 | 19.7 | 2.93 |
| cw_0.5_r3 | 5.6 | 20.0 | 2.93 |
| cw_0.5_r4 | 4.9 | 17.6 | 2.93 |
| cw_1_r0 | 5.4 | 19.5 | 2.94 |
| cw_1_r1 | 4.8 | 17.1 | 2.93 |
| cw_1_r2 | 5.2 | 18.0 | 2.93 |
| cw_1_r3 | 4.6 | 16.9 | 2.93 |
| cw_1_r4 | 5.6 | 20.0 | 2.93 |

Keep-awake request accepted by Windows: True; gaps of the memory watch above 3 min: 0 (none); lowest free memory seen 0.76 GB; low-memory events: 8.

