# A6, final form: fake-data recovery with the procedure of the real fit

Written by `research/waa_ksurface/a6_recovery.py` (libraries `fake6.py`, `ladder.py`, `ksurf.py`). Fake attenuation on the 4 563 TRAIN link-events of the learning set (real links, real rain), 0.3 dB quantization, the measured noise, the T1 rain error at its high end for the radar input and the T2 gauge model for the gauge input. **Every truth frees the F / P extension, as the real fit will (Deviation 7), and is drawn with a fresh seed.** One replicate per truth (five per size in the detectability study). Selection: a richer rung must lower the CV error by more than 0.005 dB with a paired 90 % interval that excludes zero (Deviation 5). k_wet is a bounded nuisance (Deviation 6): where two numbers are given as `x / y`, the first is with k_wet at its gauge-input estimate, the second at its radar-input estimate. For the truths in which k depends on R the share is taken inside F's measured range (restriction R2). No real attenuation is used as a target. VAL and TEST were not opened.

## The table

| truth | regime | true rung | rung returned | recovery | verdict |
|---|---|---|---|---|---|
| (a) | wetting | K0 | K0 | 1/k 3.3 / 3.1 min (truth 5.0): governed by the gate | PASS |
| (a) | drying | K0 | K0 | 1/k 17.2 / 17.0 min (truth 15.0): yes | PASS |
| (a-lag) |  | K0 |  | d by the rule: radar 2, gauge 2 (truth 2) | report |
| (b) | wetting | K0 | K0 | 1/k 3.3 / 3.3 min (truth 5.0): governed by the gate | k_dry recovered |
| (b) | drying | K0 | K0 | 1/k 16.7 / 16.7 min (truth 15.0): yes | k_dry recovered |
| (c) | wetting | K0 | K0 | 1/k 3.5 / 3.7 min (truth 5.0): governed by the gate | FAIL |
| (c) | drying | K1 | K0 | 18 % / 18 % of the domain weight: no | FAIL |

(`research/waa_ksurface/outputs/tables/a6_recovery_verdicts.csv`; every rung and lambda: `research/waa_ksurface/outputs/tables/a6_recovery_cv.csv`; every domain cell: `research/waa_ksurface/outputs/tables/a6_recovery_cells.csv`; every extension knot: `research/waa_ksurface/outputs/tables/a6_recovery_extension.csv`; figure `research/waa_ksurface/outputs/figures/a6_recovery_surfaces.png`)

**STOP: truth (c) FAILS its rule. What comes after it was not run. Nothing was patched.**

## Truth (a)

| select | rung | lambda | parameters_free | cv_db | compared_with | gain_db | gain_lo | gain_hi | replaces | lowest_cv | selected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | 0.0001 | 21 | 0.23913 |  |  |  |  |  | False | True |
| wetting | K1 | 0.0001 | 24 | 0.23852 | K0 | 0.00061 | 0.00031 | 0.00088 | False | False | False |
| wetting | K2 | 0.0001 | 25 | 0.23839 | K0 | 0.00074 | 0.00043 | 0.00112 | False | False | False |
| wetting | K3 | 0.0001 | 28 | 0.23779 | K0 | 0.00134 | 0.00081 | 0.00182 | False | True | False |
| drying | K0 | 0.0001 | 15 | 0.23877 |  |  |  |  |  | False | True |
| drying | K1 | 0.01 | 20 | 0.23881 | K0 | -4e-05 | -0.00011 | 3e-05 | False | False | False |
| drying | K2 | 0.01 | 19 | 0.23821 | K0 | 0.00056 | -0.00211 | 0.00264 | False | False | False |
| drying | K3 | 0.01 | 24 | 0.23743 | K0 | 0.00133 | -0.00124 | 0.00346 | False | True | False |

Per regime:

| regime | true_rung | selected_rung | lowest_cv_rung | lambda | rung_ok | tau_fit_min | tau_true_min | share_within_25pct | share_whole_domain | domain_cells_visited | tau_fit_alt_min | share_alt | k_note | k_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | K0 | K3 | 0.0 | True | 3.331 | 5.0 | 0.0 | 0.0 | 37 | 3.119 |  | governed by the gate | True |
| drying | K0 | K0 | K3 | 0.0 | True | 17.158 | 15.0 | 1.0 | 1.0 | 54 | 17.004 | 1.0 |  | True |

1 / k_wet: gauge-input estimate 3.33 min, radar-input estimate 3.12 min.

Residual in the 30 min after rain on the short links (859 link-events): mean -0.021 dB; highest minus lowest third of preceding rain -0.038 dB, at matched model W -0.042 dB (6 W bins).

## Truth (b)

| select | rung | lambda | parameters_free | cv_db | compared_with | gain_db | gain_lo | gain_hi | replaces | lowest_cv | selected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | 0.0001 | 21 | 0.23646 |  |  |  |  |  | False | True |
| wetting | K1 | 0.0001 | 24 | 0.23589 | K0 | 0.00058 | 0.00025 | 0.00087 | False | False | False |
| wetting | K2 | 0.0001 | 25 | 0.2357 | K0 | 0.00077 | 0.0004 | 0.00121 | False | False | False |
| wetting | K3 | 0.0001 | 28 | 0.23501 | K0 | 0.00145 | 0.00089 | 0.002 | False | True | False |
| drying | K0 | 0.0001 | 15 | 0.25371 |  |  |  |  |  | True | True |
| drying | K1 | 0.0001 | 20 | 0.25372 | K0 | -1e-05 | -0.00011 | 9e-05 | False | False | False |
| drying | K2 | 0.01 | 19 | 0.25421 | K0 | -0.0005 | -0.00295 | 0.00122 | False | False | False |
| drying | K3 | 0.01 | 24 | 0.25376 | K0 | -5e-05 | -0.00233 | 0.00161 | False | False | False |

Per regime:

| regime | true_rung | selected_rung | lowest_cv_rung | lambda | rung_ok | tau_fit_min | tau_true_min | share_within_25pct | share_whole_domain | domain_cells_visited | tau_fit_alt_min | share_alt | k_note | k_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | K0 | K3 | 0.0 | True | 3.294 | 5.0 | 0.0 | 0.0 | 37 | 3.278 |  | governed by the gate | True |
| drying | K0 | K0 | K0 | 0.0 | True | 16.706 | 15.0 | 1.0 | 1.0 | 53 | 16.694 | 1.0 |  | True |

1 / k_wet: gauge-input estimate 3.29 min, radar-input estimate 3.28 min.

Residual in the 30 min after rain on the short links (859 link-events): mean -0.020 dB; highest minus lowest third of preceding rain +0.007 dB, at matched model W -0.004 dB (6 W bins).

The extension knots (fitted with k; started from F held and P at half its true increments):

| band | quantity | knot_mmh | true | start | fitted | fitted_over_true | within_25pct | at_bound | k_wet_at |
|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | F | 5 | 2.143 | 1.857 | 2.128 | 0.993 | True | False | gauge estimate |
| 28 GHz | F | 10 | 2.5 | 1.867 | 2.128 | 0.851 | True | True | gauge estimate |
| 28 GHz | F | 20 | 2.727 | 1.877 | 2.13 | 0.781 | True | False | gauge estimate |
| 28 GHz | P | 2 | 0.2 | 0.172 | 0.282 | 1.41 | False | False | gauge estimate |
| 28 GHz | P | 5 | 0.5 | 0.322 | 0.58 | 1.16 | True | False | gauge estimate |
| 28 GHz | P | 10 | 1.0 | 0.572 | 0.669 | 0.669 | False | False | gauge estimate |
| 28 GHz | P | 20 | 2.0 | 1.072 | 1.166 | 0.583 | False | False | gauge estimate |
| 38 GHz | F | 5 | 2.143 | 1.873 | 2.16 | 1.008 | True | False | gauge estimate |
| 38 GHz | F | 10 | 2.5 | 1.883 | 2.16 | 0.864 | True | True | gauge estimate |
| 38 GHz | F | 20 | 2.727 | 1.893 | 2.16 | 0.792 | True | False | gauge estimate |
| 38 GHz | P | 2 | 0.4 | 0.341 | 0.549 | 1.373 | False | False | gauge estimate |
| 38 GHz | P | 5 | 1.0 | 0.641 | 1.11 | 1.11 | True | False | gauge estimate |
| 38 GHz | P | 10 | 2.0 | 1.141 | 1.36 | 0.68 | False | False | gauge estimate |
| 38 GHz | P | 20 | 4.0 | 2.141 | 2.358 | 0.589 | False | False | gauge estimate |
| 28 GHz | F | 5 | 2.143 | 1.857 | 2.128 | 0.993 | True | False | radar estimate |
| 28 GHz | F | 10 | 2.5 | 1.867 | 2.128 | 0.851 | True | True | radar estimate |
| 28 GHz | F | 20 | 2.727 | 1.877 | 2.13 | 0.781 | True | False | radar estimate |
| 28 GHz | P | 2 | 0.2 | 0.172 | 0.282 | 1.409 | False | False | radar estimate |
| 28 GHz | P | 5 | 0.5 | 0.322 | 0.58 | 1.16 | True | False | radar estimate |
| 28 GHz | P | 10 | 1.0 | 0.572 | 0.669 | 0.669 | False | False | radar estimate |
| 28 GHz | P | 20 | 2.0 | 1.072 | 1.165 | 0.583 | False | False | radar estimate |
| 38 GHz | F | 5 | 2.143 | 1.873 | 2.16 | 1.008 | True | False | radar estimate |
| 38 GHz | F | 10 | 2.5 | 1.883 | 2.16 | 0.864 | True | True | radar estimate |
| 38 GHz | F | 20 | 2.727 | 1.893 | 2.16 | 0.792 | True | False | radar estimate |
| 38 GHz | P | 2 | 0.4 | 0.341 | 0.549 | 1.372 | False | False | radar estimate |
| 38 GHz | P | 5 | 1.0 | 0.641 | 1.109 | 1.109 | True | False | radar estimate |
| 38 GHz | P | 10 | 2.0 | 1.141 | 1.359 | 0.68 | False | False | radar estimate |
| 38 GHz | P | 20 | 4.0 | 2.141 | 2.358 | 0.589 | False | False | radar estimate |

## Truth (c)

| select | rung | lambda | parameters_free | cv_db | compared_with | gain_db | gain_lo | gain_hi | replaces | lowest_cv | selected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | 0.0001 | 21 | 0.2337 |  |  |  |  |  | False | True |
| wetting | K1 | 0.0001 | 24 | 0.23316 | K0 | 0.00055 | 0.0002 | 0.00082 | False | False | False |
| wetting | K2 | 0.0001 | 25 | 0.23293 | K0 | 0.00078 | 0.00044 | 0.00121 | False | False | False |
| wetting | K3 | 0.0001 | 28 | 0.23223 | K0 | 0.00148 | 0.00097 | 0.00194 | False | True | False |
| drying | K0 | 0.0001 | 15 | 0.24192 |  |  |  |  |  | False | True |
| drying | K1 | 0.0001 | 20 | 0.24001 | K0 | 0.00192 | 0.00023 | 0.00362 | False | True | False |
| drying | K2 | 0.01 | 19 | 0.24291 | K0 | -0.00099 | -0.003 | 0.00112 | False | False | False |
| drying | K3 | 0.0001 | 24 | 0.24072 | K0 | 0.0012 | -0.00064 | 0.00309 | False | False | False |

Per regime:

| regime | true_rung | selected_rung | lowest_cv_rung | lambda | rung_ok | tau_fit_min | tau_true_min | share_within_25pct | share_whole_domain | domain_cells_visited | tau_fit_alt_min | share_alt | k_note | k_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | K0 | K3 | 0.0 | True | 3.539 | 5.0 | 0.0 | 0.0 | 36 | 3.695 |  | governed by the gate | True |
| drying | K1 | K0 | K1 | 0.0 | False | 16.824 | 30.0 | 0.179 | 0.179 | 51 | 16.936 | 0.179 |  | False |

1 / k_wet: gauge-input estimate 3.54 min, radar-input estimate 3.69 min.

Residual in the 30 min after rain on the short links (859 link-events): mean +0.061 dB; highest minus lowest third of preceding rain -0.038 dB, at matched model W -0.048 dB (5 W bins).

## Truth (a-lag): the rule for d on a truth with d = 2 min

| d_min | loss | tau_wet_min | tau_dry_min | lowest_loss | loss0_minus_best_lo | loss0_minus_best_hi | chosen | input |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.03456 | 4.11024 | 14.47366 | False | 0.00144 | 0.00585 | False | radar |
| 1 | 0.03393 | 3.98757 | 14.61999 | False | 0.00144 | 0.00585 | False | radar |
| 2 | 0.03362 | 3.78192 | 14.53689 | True | 0.00144 | 0.00585 | True | radar |
| 3 | 0.03386 | 3.48075 | 14.19609 | False | 0.00144 | 0.00585 | False | radar |
| 4 | 0.03455 | 3.13879 | 13.64467 | False | 0.00144 | 0.00585 | False | radar |
| 5 | 0.03532 | 2.78491 | 12.9465 | False | 0.00144 | 0.00585 | False | radar |
| 0 | 0.02657 | 3.61563 | 17.93651 | False | 0.00019 | 0.00206 | False | gauge |
| 1 | 0.02638 | 3.58971 | 18.48096 | False | 0.00019 | 0.00206 | False | gauge |
| 2 | 0.0263 | 3.49971 | 18.70571 | True | 0.00019 | 0.00206 | True | gauge |
| 3 | 0.02663 | 3.38634 | 18.20779 | False | 0.00019 | 0.00206 | False | gauge |
| 4 | 0.0271 | 3.29735 | 17.35295 | False | 0.00019 | 0.00206 | False | gauge |
| 5 | 0.02749 | 3.23091 | 16.55679 | False | 0.00019 | 0.00206 | False | gauge |

Chosen: radar d = 2, gauge d = 2. (`research/waa_ksurface/outputs/tables/a6_recovery_lag.csv`)

## Run times

| truth | minutes |
|---|---|
| a | 90.9 |
| a_lag | 1.0 |
| b | 11.7 |
| c | 12.1 |
