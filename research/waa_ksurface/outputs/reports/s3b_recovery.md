# Step 3b: the tail flag in the procedure (Deviation 12), B3' calibrated (Deviation 13), (g), detectability

Written by `research/waa_ksurface/s3b_recovery.py` (libraries `s3_recovery.py`, `a6_recovery.py`, `s1_separability.py`, `fake6.py`, `ladder.py`). Fake attenuation on all TRAIN link-events of the branch, as in step 3 (the radar of the events without a gauge ratio reads 2.5 times the true rain), fresh seeds (base 3600). Drying K0 against K1 on L-dryhigh, margin 0.00195 dB with a paired 90 % interval that excludes zero; wetting K0, not judged (R3). **New (Deviation 12): the tail flag is computed on the fake data before any fit, and the link-events of the flagged events are left out of both stages and of the scored set.** Where two numbers are given as `x / y`, the first is with k_wet at its gauge-input estimate, the second at its radar-input estimate. A returned drying K1 is named: *non-exponential drying: the drying rate falls as the antenna dries* (Deviation 13). No real attenuation is used as a target. VAL and TEST were not opened.

## Every truth with the flag (Deviation 12)

| run | seed | true rung | drying rung | K1 gain on L-dryhigh (dB) | interval | recovery (k_wet gauge / radar) | events flagged | flagged with unseen water | false removals | scored values | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| a | 3600 | K0 | K0 | 0.00074 | +0.00021 to +0.00212 | 1 / k_dry 15.3 / 15.7 min, truth 15 | 0 of 13 | 0 | 0 | 2405 | PASS |
| c | 3601 | K1 | K1 | 0.01575 | +0.01283 to +0.02048 | 100 % / 100 % of the domain weight within 25 % | 0 of 13 | 0 | 0 | 2062 | PASS |
| d | 3602 | K2 | K0 | -0.00205 | -0.00405 to -0.00121 | 95 % / 95 % of the domain weight within 25 % | 0 of 13 | 0 | 0 | 2233 | PASS |
| d_high | 3603 | K0 | K0 | 0.00103 | +0.00019 to +0.00231 | 1 / k_dry 9.0 / 9.6 min, truth 10 | 0 of 13 | 0 | 0 | 552 | PASS |
| e | 3604 | none | K1 | 0.02087 | +0.01203 to +0.02384 | 1 / k_dry 41.0 / 41.9 min (fit; no single truth) | 0 of 13 | 0 | 0 | 950 | report |
| f | 3605 | K0 | K0 | 0.00026 | -0.00056 to +0.00283 | 1 / k_dry 18.1 / 16.7 min, truth 15 | 2 of 13 | 2 | 0 | 3050 | no spurious k_dry |
| f_added | 3607 | K0 | K1 | 0.00347 | +0.00244 to +0.00440 | 1 / k_dry 16.5 / 17.9 min, truth 15 | 1 of 41 | 1 | 0 | 5440 | FAIL, material |
| g | 3606 | K0 | K1 | 0.00531 | +0.00383 to +0.00819 | 1 / k_dry 18.9 / 19.0 min, truth 15 | 0 of 13 | 0 | 0 | 4520 | k corrupted |


## B3': the history test, calibrated (Deviation 13)

Statistic: the after-rain residual (30 min after the last wet bin, short links of the unflagged learning events) at matched W (bins of 0.3 dB of the fitted W at the end of the rain), highest minus lowest third of preceding rain, from the drying K1 fit. Null = (c) (one state, k_dry rising with W); signal = (e) (two stores). Separable if every (e) value is above the largest (c) value.

| truth | run | seed | drying rung | B3' statistic, K1 fit (dB) | matched bins | link-events | plain high - low, K1 fit (dB) | old B3, selected fit (dB) |
|---|---|---|---|---|---|---|---|---|
| c | c | 3601 | K1 | -0.06556 | 6 | 859 | -0.06995 | -0.06556 |
| c | c_r1 | 3620 | K1 | -0.03708 | 6 | 859 | -0.03774 | -0.03708 |
| c | c_r2 | 3621 | K1 | -0.04991 | 6 | 859 | -0.0497 | -0.04991 |
| c | c_r3 | 3622 | K1 | -0.05921 | 5 | 859 | -0.06758 | -0.05921 |
| c | c_r4 | 3623 | K1 | -0.00129 | 4 | 859 | 0.00311 | -0.00129 |
| e | e | 3604 | K1 | -0.00284 | 5 | 859 | -0.00298 | -0.00284 |
| e | e_r1 | 3630 | K1 | 0.0237 | 5 | 859 | 0.03259 | 0.0237 |
| e | e_r2 | 3631 | K1 | 0.00087 | 5 | 859 | -0.02833 | 0.00087 |
| e | e_r3 | 3632 | K1 | 0.02466 | 5 | 859 | 0.04663 | 0.02466 |
| e | e_r4 | 3633 | K1 | 0.01633 | 5 | 859 | 0.02532 | 0.01633 |

(c): -0.0656 to -0.0013 dB; (e): -0.0028 to +0.0247 dB; (e) above max (c) in 4 of 5.

**Not separable: finding "W-dependent rate vs second store not identifiable with this data"; P5 becomes restriction R5.** B3's 0.10 dB threshold is retired.


## The flag's blind spot: unseen water in added events, flag extended to all 51 events

Unseen water in 7 of the added events (6 of them checkable: at least 5 checkable short link-events). Caught by the extended flag: 1 of 7 (of the checkable ones: 1 of 6); false removals 0. Drying rung K1, 1 / k_dry 16.5 / 17.9 min, truth 15. **FAIL, material**; material: yes.

The flag's table holds 41 of the 51 events: an event with no short-link link-event that has 20 bins after the end of its rain cannot be checked at all, so the other 10 are outside the flag whatever their tails do.


## Check-report line: the flag's blind spot

| flag blind-spot variant | flag extended? | catch rate | drying verdict | material? | restriction or none |
|---|---|---|---|---|---|
| f_added (seed 3607) | yes: all 51 events, 41 checkable at all | 1 of 7; false removals 0 | K1 (truth K0); 1 / k_dry 16.5 / 17.9 min, truth 15 | yes | added events enter A7's extended variant only where the flag covers them; otherwise A7 falls back to the learning events, with the S1 margin of that variant (L-dryhigh, the current variant: 0.00291 dB, s1_separability.md) |

## (g): a larger spread between links

Drying rung K1, 1 / k_dry 18.9 / 19.0 min, truth 15: **k corrupted**; per-link question reopened: yes.


## Detectability of k_dry(W) under the updated procedure

k_dry = (1/15)·(1 + a·W) / (1 + a), k_wet = 1/5; ratio k(W = 2 dB) / k(W = 0) = 1 + 2a. Detected = drying K1 returned in at least 4 of 5 replicates (the replicates of a ratio stop when decided; the sweep ends at the first ratio not detected). Ratio 3 is truth (c): its replicates are those of B3'.

| size_a | ratio | replicates_run | hits (K1 returned) | detected | runs |
|---|---|---|---|---|---|
| 1.0 | 3.0 | 5 | 5 | True | c, c_r1, c_r2, c_r3, c_r4 |
| 0.5 | 2.0 | 4 | 4 | True | cw_0.5_r0, cw_0.5_r1, cw_0.5_r2, cw_0.5_r3 |
| 0.25 | 1.5 | 4 | 4 | True | cw_0.25_r0, cw_0.25_r1, cw_0.25_r2, cw_0.25_r3 |

**Smallest ratio detected in at least 4 of 5: 1.5.**


## Execution

Runs: 24; wall time of the runs 1.4 h, processor time 6.1 h; largest peak memory 2.93 GB.

| run | minutes | cpu minutes | peak memory (GB) | B3' (K1 fit, dB) | B3 (selected fit, dB) |
|---|---|---|---|---|---|
| a | 4.9 | 19.0 | 2.93 | -0.0664 | -0.0663 |
| c | 3.2 | 15.0 | 2.93 | -0.0656 | -0.0656 |
| d | 3.5 | 16.4 | 2.85 | -0.0253 | -0.0259 |
| d_high | 3.9 | 17.1 | 2.93 | -0.0312 | -0.0323 |
| e | 3.2 | 13.9 | 2.93 | -0.0028 | -0.0028 |
| f | 3.4 | 14.0 | 2.93 | -0.0444 | -0.045 |
| f_added | 3.6 | 15.1 | 2.93 | -0.071 | -0.071 |
| c_r1 | 3.4 | 14.6 | 2.93 | -0.0371 | -0.0371 |
| c_r2 | 4.0 | 17.3 | 2.93 | -0.0499 | -0.0499 |
| c_r3 | 3.2 | 13.7 | 2.93 | -0.0592 | -0.0592 |
| c_r4 | 2.9 | 12.3 | 2.93 | -0.0013 | -0.0013 |
| e_r1 | 3.3 | 14.2 | 2.93 | 0.0237 | 0.0237 |
| e_r2 | 3.0 | 12.8 | 2.93 | 0.0009 | 0.0009 |
| e_r3 | 3.5 | 14.9 | 2.93 | 0.0247 | 0.0247 |
| e_r4 | 3.4 | 14.4 | 2.93 | 0.0163 | 0.0163 |
| g | 3.7 | 15.8 | 2.93 | -0.0985 | -0.0985 |
| cw_0.5_r0 | 4.0 | 17.4 | 2.93 | -0.0839 | -0.0839 |
| cw_0.5_r1 | 3.4 | 14.5 | 2.93 | -0.0114 | -0.0114 |
| cw_0.5_r2 | 3.9 | 17.0 | 2.93 | -0.0034 | -0.0034 |
| cw_0.5_r3 | 3.0 | 12.7 | 2.93 | -0.036 | -0.036 |
| cw_0.25_r0 | 3.5 | 14.9 | 2.93 | -0.0144 | -0.0144 |
| cw_0.25_r1 | 4.1 | 17.8 | 2.93 | -0.0464 | -0.0464 |
| cw_0.25_r2 | 3.3 | 14.0 | 2.93 | -0.0379 | -0.0379 |
| cw_0.25_r3 | 4.3 | 18.9 | 2.88 | -0.1304 | -0.1304 |


Keep-awake request accepted by Windows: True; gaps of the memory watch above 3 min: 0 (none); lowest free memory seen 0.48 GB; low-memory events: 1.
