# A6: fake-data recovery (the gates for the real fit), under Deviations 5 and 6

Written by `research/waa_ksurface/a6_recovery.py` (libraries `fake6.py`, `ladder.py`, `ksurf.py`). Fake attenuation on the 4 563 TRAIN link-events of the learning set (real links, real rain), 0.3 dB quantization, the measured noise, the T1 rain error at its high end for the radar input and the T2 gauge model for the gauge input. One replicate per truth. Selection: a richer rung must lower the CV error by more than 0.005 dB with a paired 90 % interval that excludes zero (Deviation 5). k_wet: a bounded nuisance if the gate passes (Deviation 6); where two numbers are given as `x / y`, the first is with k_wet at its gauge-input estimate, the second at its radar-input estimate. No real attenuation is used as a target. VAL and TEST were not opened.

## The k_wet sensitivity gate (Deviation 6)

Truth (a): 1 / k_wet = 5 min, 1 / k_dry = 15 min. All links, radar input, K0 in both regimes, 1 / k_wet HELD; k_dry and the F / P extension free.

| tau_wet_held_min | cv_error_db | train_loss | tau_dry_fit_min | tau_dry_true_min | k_dry_within_25pct |
|---|---|---|---|---|---|
| 2.5 | 0.24383 | 0.02811 | 12.54728 | 15.0 | True |
| 5.0 | 0.24002 | 0.0273 | 13.88221 | 15.0 | True |
| 10.0 | 0.24296 | 0.02803 | 16.24482 | 15.0 | True |

CV error changes by **0.00381 dB** across the three values (margin 0.005 dB). **Gate: PASS: k_wet is a bounded nuisance.** (`research/waa_ksurface/outputs/tables/a6_recovery_wet_gate.csv`)

## The table

| truth | regime | true rung | rung returned | recovery | verdict |
|---|---|---|---|---|---|
| (a) | wetting | K0 | K0 | 1/k 3.7 / 9.9 min (truth 5.0): governed by the gate | FAIL |
| (a) | drying | K0 | K0 | 1/k 11.6 / 14.5 min (truth 15.0): no | FAIL |

(`research/waa_ksurface/outputs/tables/a6_recovery_verdicts.csv`; every rung and lambda: `research/waa_ksurface/outputs/tables/a6_recovery_cv.csv`; every domain cell: `research/waa_ksurface/outputs/tables/a6_recovery_cells.csv`)

**STOP: truth (a) FAILS its rule. What comes after it was not run. Nothing was patched.**

## Truth (a)

| select | rung | lambda | parameters_free | cv_db | compared_with | gain_db | gain_lo | gain_hi | replaces | lowest_cv | selected |
|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | 0.0001 | 7 | 0.24141 |  |  |  |  |  | False | True |
| wetting | K1 | 0.0001 | 10 | 0.24041 | K0 | 0.00099 | 0.00057 | 0.00134 | False | False | False |
| wetting | K2 | 0.0001 | 11 | 0.24005 | K0 | 0.00136 | 0.00089 | 0.00196 | False | False | False |
| wetting | K3 | 0.0001 | 14 | 0.23885 | K0 | 0.00256 | 0.00194 | 0.0033 | False | True | False |
| drying | K0 | 0.0001 | 1 | 0.27326 |  |  |  |  |  | False | True |
| drying | K1 | 0.0001 | 6 | 0.27112 | K0 | 0.00214 | 0.00049 | 0.00363 | False | False | False |
| drying | K2 | 0.0001 | 5 | 0.27067 | K0 | 0.00258 | -0.00053 | 0.00516 | False | False | False |
| drying | K3 | 0.0001 | 10 | 0.2705 | K0 | 0.00276 | -0.00023 | 0.00522 | False | True | False |

Per regime:

| regime | true_rung | selected_rung | lowest_cv_rung | lambda | rung_ok | tau_fit_min | tau_true_min | share_within_25pct | domain_cells_visited | tau_fit_alt_min | share_alt | k_note | k_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | K0 | K3 | 0.0 | True | 3.728 | 5.0 | 0.0 | 38 | 9.908 |  | governed by the gate | True |
| drying | K0 | K0 | K3 | 0.0 | True | 11.601 | 15.0 | 0.0 | 51 | 14.458 | 1.0 |  | False |

1 / k_wet: gauge-input estimate 3.73 min, radar-input estimate 9.91 min.

Residual in the 30 min after rain on the short links (859 link-events): mean +0.035 dB; highest minus lowest third of preceding rain -0.055 dB, at matched model W -0.048 dB (5 W bins).

## Run times

| truth | minutes |
|---|---|
| a | 5.0 |
