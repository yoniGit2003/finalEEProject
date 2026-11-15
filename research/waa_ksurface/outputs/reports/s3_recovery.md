# Step 3: the fake-data truths with the reduced ladders (Deviations 9 and 11)

Written by `research/waa_ksurface/s3_recovery.py` (libraries `a6_recovery.py`, `s1_separability.py`, `fake6.py`, `ladder.py`, `ksurf.py`). Fake attenuation on all TRAIN link-events of the branch (real links, real rain; the radar of the events without a gauge ratio reads 2.5 times the true rain), the noise and the input errors of A6, extension free (Deviation 7), fresh seeds. **Drying selection (Deviation 9): every drying rung is scored on L-dryhigh (held-out 5 min values of the dry 120 min after the rain, link-events whose first dry value is at least 1 dB), learning set = the 14 events with a gauge ratio plus the dry-path values of the other events; a richer rung needs a gain above 0.00195 dB (the largest gain on a constant truth in S1) with a paired 90 % interval that excludes zero. The drying ladder is K0 against K1 only, and there is no wetting ladder: wetting is K0 and is not judged, k_wet is not identified (Deviation 11, restrictions R3 and R4).** In (d) and (d-high) the true k depends on R, which this ladder cannot return: they PASS if drying stays K0 and the constant is recovered. (c) and (d) were fitted before Deviation 11 with the full ladders; their K0 and K1 fits are those of the reduced procedure, so they are judged again and not refitted, and their tables still list the K2 and K3 rows, which take no part. In the tables `cv_db` of a drying rung is its error on L-dryhigh. Wetting: Deviation 5 (0.005 dB on every counted value of the gauge-input set). Where two numbers are given as `x / y`, the first is with k_wet at its gauge-input estimate, the second at its radar-input estimate. Judged on the learning set as in A6. No real attenuation is used as a target. VAL and TEST were not opened.

## The table

| truth | regime | true rung | rung returned | recovery | verdict |
|---|---|---|---|---|---|
| (c) | drying | K1 | K1 | 97 % / 97 % of the domain weight: yes | PASS |
| (d) | drying | K2 (outside the ladder, R4) | K0 | 93 % / 93 % of the domain weight (whole domain 92 %): yes | PASS |

(`research/waa_ksurface/outputs/tables/s3_recovery_verdicts.csv`; every rung and lambda: `research/waa_ksurface/outputs/tables/s3_recovery_cv.csv`; every domain cell: `research/waa_ksurface/outputs/tables/s3_recovery_cells.csv`; every extension knot: `research/waa_ksurface/outputs/tables/s3_recovery_extension.csv`; figure `research/waa_ksurface/outputs/figures/s3_recovery_surfaces.png`)

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

## Run times

| truth | minutes |
|---|---|
| c | 68.2 |
| d | 48.0 |

## Every run: the drying gain against the margin, the size of the scored set, peak memory

| truth | drying rung | attribution (Deviation 10) | K1 gain on L-dryhigh (dB) | margin (dB) | scored values | scored link-events | peak memory (GB) | minutes |
|---|---|---|---|---|---|---|---|---|
| c | K1 | W | 0.01348 | 0.00195 | 2409 | 109 | 2.84 | 68.2 |
| d | K0 | no dependence | -0.00332 | 0.00195 | 1728 | 82 | 2.86 | 48.0 |

