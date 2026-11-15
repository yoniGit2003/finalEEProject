# Step 3: the fake-data truths with the reduced ladders (Deviations 9 and 11)

Written by `research/waa_ksurface/s3_recovery.py` (libraries `a6_recovery.py`, `s1_separability.py`, `fake6.py`, `ladder.py`, `ksurf.py`). Fake attenuation on all TRAIN link-events of the branch (real links, real rain; the radar of the events without a gauge ratio reads 2.5 times the true rain), the noise and the input errors of A6, extension free (Deviation 7), fresh seeds. **Drying selection (Deviation 9): every drying rung is scored on L-dryhigh (held-out 5 min values of the dry 120 min after the rain, link-events whose first dry value is at least 1 dB), learning set = the 14 events with a gauge ratio plus the dry-path values of the other events; a richer rung needs a gain above 0.00195 dB (the largest gain on a constant truth in S1) with a paired 90 % interval that excludes zero. The drying ladder is K0 against K1 only, and there is no wetting ladder: wetting is K0 and is not judged, k_wet is not identified (Deviation 11, restrictions R3 and R4).** In (d) and (d-high) the true k depends on R, which this ladder cannot return: they PASS if drying stays K0 and the constant is recovered. (c) and (d) were fitted before Deviation 11 with the full ladders; their K0 and K1 fits are those of the reduced procedure, so they are judged again and not refitted, and their tables still list the K2 and K3 rows, which take no part. In the tables `cv_db` of a drying rung is its error on L-dryhigh. Wetting: Deviation 5 (0.005 dB on every counted value of the gauge-input set). Where two numbers are given as `x / y`, the first is with k_wet at its gauge-input estimate, the second at its radar-input estimate. Judged on the learning set as in A6. No real attenuation is used as a target. VAL and TEST were not opened.

## The table

| truth | regime | true rung | rung returned | recovery | verdict |
|---|---|---|---|---|---|
| (c) | drying | K1 | K1 | 97 % / 97 % of the domain weight: yes | PASS |
| (d) | drying | K2 (outside the ladder, R4) | K0 | 93 % / 93 % of the domain weight (whole domain 92 %): yes | PASS |
| (d_high) | drying | K0 | K0 | 1/k 8.9 / 9.9 min (truth 10.0): yes | PASS |
| (e) | drying | none | K1 | 1/k 39.4 min at R = 0, W = 0 | report |
| (f) | drying | K0 | K1 | 1/k 16.6 / 17.3 min (truth 15.0): yes | spurious k_dry, restored by the exclusion |
| (f), flagged events excluded | drying | K0 | K0 | 1/k 12.4 / 13.5 min (truth 15.0): yes | spurious k_dry, restored by the exclusion |

(`research/waa_ksurface/outputs/tables/s3_recovery_verdicts.csv`; every rung and lambda: `research/waa_ksurface/outputs/tables/s3_recovery_cv.csv`; every domain cell: `research/waa_ksurface/outputs/tables/s3_recovery_cells.csv`; every extension knot: `research/waa_ksurface/outputs/tables/s3_recovery_extension.csv`; figure `research/waa_ksurface/outputs/figures/s3_recovery_surfaces.png`)

**STOP: truth (f) FAILS its rule. What comes after it was not run. Nothing was patched.**

## Truth (c)

| select | rung | lambda | parameters_free | cv_db | compared_with | gain_db | gain_lo | gain_hi | replaces | lowest_cv | selected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | 0.0001 | 21 | 0.23519 |  |  |  |  |  | False | True |
| wetting | K1 | 0.0001 | 24 | 0.23452 | K0 | 0.00067 | 0.00022 | 0.00107 | False | False | False |
| wetting | K2 | 0.0001 | 25 | 0.23463 | K0 | 0.00056 | 0.0003 | 0.00089 | False | False | False |
| wetting | K3 | 0.0001 | 28 | 0.23382 | K0 | 0.00137 | 0.00084 | 0.00187 | False | True | False |
| drying | K0 | 0.0001 | 15 | 0.21337 |  |  |  |  |  | False | False |
| drying | K1 | 0.0001 | 20 | 0.19988 | K0 | 0.01348 | 0.00894 | 0.02041 | True | False | True |
| drying | K2 | 0.01 | 19 | 0.20739 | K0 | 0.00597 | 0.00347 | 0.00983 | False | False | False |
| drying | K3 | 0.0001 | 24 | 0.19978 | K1 | 0.00011 | -0.00375 | 0.00511 | False | True | False |

Per regime:

| regime | true_rung | selected_rung | lowest_cv_rung | lambda | rung_ok | tau_fit_min | tau_true_min | share_within_25pct | share_whole_domain | domain_cells_visited | tau_fit_alt_min | share_alt | k_note | k_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| drying | K1 | K1 | K3 | 0.0 | True | 33.458 | 30.0 | 0.972 | 0.972 | 51 | 34.161 | 0.972 |  | True |

1 / k_wet: gauge-input estimate 3.55 min, radar-input estimate 4.74 min.

Residual in the 30 min after rain on the short links (859 link-events): mean +0.003 dB; highest minus lowest third of preceding rain -0.098 dB, at matched model W -0.088 dB (6 W bins).

## Truth (d)

| select | rung | lambda | parameters_free | cv_db | compared_with | gain_db | gain_lo | gain_hi | replaces | lowest_cv | selected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | 0.0001 | 21 | 0.24261 |  |  |  |  |  | False | True |
| wetting | K1 | 0.0001 | 24 | 0.24119 | K0 | 0.00141 | 0.00091 | 0.00186 | False | False | False |
| wetting | K2 | 0.0001 | 25 | 0.24104 | K0 | 0.00157 | 0.00104 | 0.00221 | False | False | False |
| wetting | K3 | 0.0001 | 28 | 0.23937 | K0 | 0.00323 | 0.00254 | 0.00406 | False | True | False |
| drying | K0 | 0.0001 | 15 | 0.18335 |  |  |  |  |  | True | True |
| drying | K1 | 0.01 | 20 | 0.18667 | K0 | -0.00332 | -0.00402 | -0.0012 | False | False | False |
| drying | K2 | 0.0001 | 19 | 0.18666 | K0 | -0.00331 | -0.00417 | -0.00095 | False | False | False |
| drying | K3 | 0.01 | 24 | 0.18774 |  |  |  |  |  | False | False |

Per regime:

| regime | true_rung | selected_rung | lowest_cv_rung | lambda | rung_ok | tau_fit_min | tau_true_min | share_within_25pct | share_whole_domain | domain_cells_visited | tau_fit_alt_min | share_alt | k_note | k_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| drying | K2 (outside the ladder, R4) | K0 | K0 | 0.0 | False | 9.504 | 10.0 | 0.928 | 0.916 | 54 | 10.058 | 0.928 |  | True |

1 / k_wet: gauge-input estimate 3.27 min, radar-input estimate 4.74 min.

Residual in the 30 min after rain on the short links (859 link-events): mean -0.005 dB; highest minus lowest third of preceding rain -0.049 dB, at matched model W -0.046 dB (6 W bins).

## Truth (d_high)

| select | rung | lambda | parameters_free | cv_db | compared_with | gain_db | gain_lo | gain_hi | replaces | lowest_cv | selected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | 0.0001 | 21 | 0.23633 |  |  |  |  |  | True | True |
| drying | K0 | 0.0001 | 15 | 0.1944 |  |  |  |  |  | True | True |
| drying | K1 | 0.0001 | 20 | 0.19443 | K0 | -3e-05 | -0.00119 | 0.00034 | False | False | False |

Per regime:

| regime | true_rung | selected_rung | lowest_cv_rung | lambda | rung_ok | tau_fit_min | tau_true_min | share_within_25pct | share_whole_domain | domain_cells_visited | tau_fit_alt_min | share_alt | k_note | k_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| drying | K0 | K0 | K0 | 0.0 | True | 8.913 | 10.0 | 1.0 | 0.988 | 54 | 9.86 | 1.0 |  | True |

1 / k_wet: gauge-input estimate 4.17 min, radar-input estimate 7.00 min.

Residual in the 30 min after rain on the short links (859 link-events): mean -0.004 dB; highest minus lowest third of preceding rain -0.037 dB, at matched model W -0.035 dB (7 W bins).

The extension knots (fitted with k; started from F held and P at half its true increments):

| band | quantity | knot_mmh | true | start | fitted | fitted_over_true | within_25pct | at_bound | k_wet_at |
|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | F | 5 | 2.143 | 1.857 | 2.14 | 0.999 | True | False | gauge estimate |
| 28 GHz | F | 10 | 2.5 | 1.867 | 2.255 | 0.902 | True | False | gauge estimate |
| 28 GHz | F | 20 | 2.727 | 1.877 | 2.59 | 0.95 | True | False | gauge estimate |
| 28 GHz | P | 2 | 0.2 | 0.172 | 0.201 | 1.005 | True | False | gauge estimate |
| 28 GHz | P | 5 | 0.5 | 0.322 | 0.449 | 0.899 | True | False | gauge estimate |
| 28 GHz | P | 10 | 1.0 | 0.572 | 0.799 | 0.799 | True | False | gauge estimate |
| 28 GHz | P | 20 | 2.0 | 1.072 | 1.27 | 0.635 | False | False | gauge estimate |
| 38 GHz | F | 5 | 2.143 | 1.873 | 2.179 | 1.017 | True | False | gauge estimate |
| 38 GHz | F | 10 | 2.5 | 1.883 | 2.179 | 0.872 | True | True | gauge estimate |
| 38 GHz | F | 20 | 2.727 | 1.893 | 2.477 | 0.908 | True | False | gauge estimate |
| 38 GHz | P | 2 | 0.4 | 0.341 | 0.45 | 1.126 | True | False | gauge estimate |
| 38 GHz | P | 5 | 1.0 | 0.641 | 0.925 | 0.925 | True | False | gauge estimate |
| 38 GHz | P | 10 | 2.0 | 1.141 | 1.574 | 0.787 | True | False | gauge estimate |
| 38 GHz | P | 20 | 4.0 | 2.141 | 2.613 | 0.653 | False | False | gauge estimate |
| 28 GHz | F | 5 | 2.143 | 1.857 | 2.193 | 1.024 | True | False | radar estimate |
| 28 GHz | F | 10 | 2.5 | 1.867 | 2.208 | 0.883 | True | False | radar estimate |
| 28 GHz | F | 20 | 2.727 | 1.877 | 2.661 | 0.976 | True | False | radar estimate |
| 28 GHz | P | 2 | 0.2 | 0.172 | 0.23 | 1.149 | True | False | radar estimate |
| 28 GHz | P | 5 | 0.5 | 0.322 | 0.465 | 0.93 | True | False | radar estimate |
| 28 GHz | P | 10 | 1.0 | 0.572 | 0.833 | 0.833 | True | False | radar estimate |
| 28 GHz | P | 20 | 2.0 | 1.072 | 1.315 | 0.658 | False | False | radar estimate |
| 38 GHz | F | 5 | 2.143 | 1.873 | 2.199 | 1.026 | True | False | radar estimate |
| 38 GHz | F | 10 | 2.5 | 1.883 | 2.199 | 0.88 | True | True | radar estimate |
| 38 GHz | F | 20 | 2.727 | 1.893 | 2.46 | 0.902 | True | False | radar estimate |
| 38 GHz | P | 2 | 0.4 | 0.341 | 0.496 | 1.239 | True | False | radar estimate |
| 38 GHz | P | 5 | 1.0 | 0.641 | 0.963 | 0.963 | True | False | radar estimate |
| 38 GHz | P | 10 | 2.0 | 1.141 | 1.606 | 0.803 | True | False | radar estimate |
| 38 GHz | P | 20 | 4.0 | 2.141 | 2.698 | 0.674 | False | False | radar estimate |

## Truth (e)

| select | rung | lambda | parameters_free | cv_db | compared_with | gain_db | gain_lo | gain_hi | replaces | lowest_cv | selected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | 0.0001 | 21 | 0.22513 |  |  |  |  |  | True | True |
| drying | K0 | 0.0001 | 15 | 0.27929 |  |  |  |  |  | False | False |
| drying | K1 | 0.0001 | 20 | 0.26376 | K0 | 0.01554 | 0.01005 | 0.02081 | True | True | True |

Per regime:

| regime | true_rung | selected_rung | lowest_cv_rung | lambda | rung_ok | tau_fit_min | tau_true_min | share_within_25pct | share_whole_domain | domain_cells_visited | tau_fit_alt_min | share_alt | k_note | k_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| drying | none | K1 | K1 | 0.0 | False | 39.403 |  |  |  | 41 | 43.449 |  |  |  |

1 / k_wet: gauge-input estimate 7.00 min, radar-input estimate 15.27 min.

Residual in the 30 min after rain on the short links (859 link-events): mean -0.030 dB; highest minus lowest third of preceding rain +0.018 dB, at matched model W +0.014 dB (5 W bins).

## Truth (f)

| select | rung | lambda | parameters_free | cv_db | compared_with | gain_db | gain_lo | gain_hi | replaces | lowest_cv | selected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | 0.0001 | 21 | 0.23789 |  |  |  |  |  | True | True |
| drying | K0 | 0.0001 | 15 | 0.42639 |  |  |  |  |  | False | False |
| drying | K1 | 0.0001 | 20 | 0.4233 | K0 | 0.00309 | 0.0028 | 0.00345 | True | True | True |

Per regime:

| regime | true_rung | selected_rung | lowest_cv_rung | lambda | rung_ok | tau_fit_min | tau_true_min | share_within_25pct | share_whole_domain | domain_cells_visited | tau_fit_alt_min | share_alt | k_note | k_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| drying | K0 | K1 | K1 | 0.0 | False | 16.641 | 15.0 | 0.902 | 0.902 | 51 | 17.27 | 0.948 |  | True |

1 / k_wet: gauge-input estimate 3.42 min, radar-input estimate 6.00 min.

Residual in the 30 min after rain on the short links (859 link-events): mean +0.199 dB; highest minus lowest third of preceding rain +0.013 dB, at matched model W +0.074 dB (6 W bins).

The tail flag per event (`changed_in_truth` = the events that were given unseen rain):

| event_id | link_events | qualifying | waiting_or_rising | flagged | changed_in_truth |
|---|---|---|---|---|---|
| E20150617_0705 | 55 | 18 | 0 | False | False |
| E20150618_0530 | 67 | 37 | 0 | False | False |
| E20150622_0930 | 67 | 2 | 0 | False | False |
| E20150705_1910 | 67 | 8 | 8 | True | True |
| E20150707_1235 | 67 | 9 | 0 | False | False |
| E20150725_0440 | 67 | 14 | 0 | False | False |
| E20150727_2140 | 67 | 16 | 16 | True | True |
| E20150804_1830 | 67 | 8 | 0 | False | False |
| E20150810_1545 | 67 | 8 | 0 | False | False |
| E20150815_1305 | 67 | 3 | 0 | False | False |
| E20150825_0145 | 67 | 10 | 10 | True | True |
| E20150826_0520 | 67 | 47 | 0 | False | False |
| E20150826_1620 | 67 | 4 | 0 | False | False |

With the flagged events excluded:

| select | rung | lambda | parameters_free | cv_db | compared_with | gain_db | gain_lo | gain_hi | replaces | lowest_cv | selected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | 0.0001 | 21 | 0.23625 |  |  |  |  |  | True | True |
| drying | K0 | 0.0001 | 15 | 0.1869 |  |  |  |  |  | False | True |
| drying | K1 | 0.0001 | 20 | 0.18553 | K0 | 0.00137 | 0.00069 | 0.00364 | False | True | False |

| regime | true_rung | selected_rung | lowest_cv_rung | lambda | rung_ok | tau_fit_min | tau_true_min | share_within_25pct | share_whole_domain | domain_cells_visited | tau_fit_alt_min | share_alt | k_note | k_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| drying | K0 | K0 | K1 | 0.0 | True | 12.423 | 15.0 | 1.0 | 1.0 | 51 | 13.544 | 1.0 |  | True |

(`research/waa_ksurface/outputs/tables/s3_recovery_flags.csv`)

## Run times

| truth | minutes |
|---|---|
| c | 68.2 |
| d | 48.0 |
| d_high | 3.9 |
| e | 3.3 |
| f | 146.9 |

## (f): the spurious gain on the adopted loss against the margin

Gain of drying K1 over K0 on L-dryhigh: +0.00309 dB (interval +0.00280 to +0.00345); with the flagged events excluded +0.00137 dB (+0.00069 to +0.00364); margin 0.00195 dB. Verdict of (f): spurious k_dry, restored by the exclusion.


## Every run: the drying gain against the margin, the size of the scored set, peak memory

| truth | drying rung | attribution (Deviation 10) | K1 gain on L-dryhigh (dB) | margin (dB) | scored values | scored link-events | peak memory (GB) | minutes |
|---|---|---|---|---|---|---|---|---|
| c | K1 | W | 0.01348 | 0.00195 | 2409 | 109 | 2.84 | 68.2 |
| d | K0 | no dependence | -0.00332 | 0.00195 | 1728 | 82 | 2.86 | 48.0 |
| d_high | K0 | no dependence | -3e-05 | 0.00195 | 871 | 43 | 2.93 | 3.9 |
| e | K1 | W | 0.01554 | 0.00195 | 687 | 32 | 2.93 | 3.3 |
| f | K1 | W | 0.00309 | 0.00195 | 6733 | 297 | 3.16 | 146.9 |

