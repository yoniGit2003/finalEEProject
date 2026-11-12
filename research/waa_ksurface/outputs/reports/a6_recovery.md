# A6: fake-data recovery (the gates for the real fit)

Written by `research/waa_ksurface/a6_recovery.py` (libraries `fake6.py`, `ladder.py`, `ksurf.py`). Fake attenuation on the 4 563 TRAIN link-events of the learning set (real links, real rain), 0.3 dB quantization, the measured noise, the T1 rain error at its high end for the radar input and the T2 gauge model for the gauge input. One replicate per truth. No real attenuation is used as a target. VAL and TEST were not opened. Rules: `docs/waa_ksurface/decisions.md`, step A6.

## The table

| truth | regime | true rung | rung returned | k within 25 % | verdict |
|---|---|---|---|---|---|
| (a) | wetting | K0 | K3 | 1/k 45.5 min (truth 5.0): no | FAIL |
| (a) | drying | K0 | K1 | 1/k 23.0 min (truth 15.0): no | FAIL |

(`research/waa_ksurface/outputs/tables/a6_recovery_verdicts.csv`; every rung and λ with its held-out loss: `research/waa_ksurface/outputs/tables/a6_recovery_cv.csv`; every domain cell: `research/waa_ksurface/outputs/tables/a6_recovery_cells.csv`)

**STOP: truth (a) FAILS its rule. The truths after it were not run. Nothing was patched.**

## Truth (a): held-out loss per rung (λ of the lowest held-out loss)

| select | rung | lambda | parameters_free | train_loss | cv_loss | minus_best_lo | minus_best_hi | lowest_cv | selected |
|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | 0.0001 | 7 | 0.029055 | 0.029139 | 0.000445 | 0.000833 | False | False |
| wetting | K1 | 0.0001 | 10 | 0.028884 | 0.0289 | 0.000254 | 0.000542 | False | False |
| wetting | K2 | 0.0001 | 11 | 0.028708 | 0.028813 | 0.00018 | 0.000388 | False | False |
| wetting | K3 | 0.0001 | 14 | 0.028503 | 0.028525 | 0.0 | 0.0 | True | True |
| drying | K0 | 0.0001 | 1 | 0.03751 | 0.037596 | 0.000243 | 0.001858 | False | False |
| drying | K1 | 0.0001 | 6 | 0.036655 | 0.036672 | -0.000205 | 0.000389 | False | True |
| drying | K2 | 0.0001 | 5 | 0.036628 | 0.036733 | 3.2e-05 | 0.000379 | False | False |
| drying | K3 | 0.0001 | 10 | 0.036444 | 0.036554 | 0.0 | 0.0 | True | False |

Per regime:

| regime | true_rung | selected_rung | lowest_cv_rung | lambda | rung_ok | tau_fit_at_R0_W0_min | tau_true_at_R0_W0_min | constant_within_25pct | share_within_25pct | domain_cells_visited | k_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|
| wetting | K0 | K3 | K3 | 0.0 | False | 45.463 | 5.0 | False | 0.139 | 38 | False |
| drying | K0 | K1 | K3 | 0.0 | False | 22.961 | 15.0 | False | 0.121 | 51 | False |

Residual in the 30 min after rain on the short links (859 link-events): mean +0.026 dB; highest minus lowest third of preceding rain -0.041 dB, at matched model W -0.029 dB (5 W bins).

## Run times

| truth | minutes |
|---|---|
| a | 4.0 |
