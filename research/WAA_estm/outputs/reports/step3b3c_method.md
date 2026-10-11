# WAA_estm step 3b-3c: the discovery method finished on fake tails (rules fixed before the run)

Written by `research/WAA_estm/step3b3c_method.py` (libraries `waa_curve.py`, `waa_pooled.py`, `waa_sindy.py`; settings
`waa.curve3b3c`). Fake data only: no real tail is fitted, the frozen tail set v1 is not changed.

**Overall: PARTIAL.** Per part: A PARTIAL, B PARTIAL, C PARTIAL, D FAIL; E is information only.

## 1. What was run

- **Data of B-E**: the cached fake datasets and free-curve fits of step 3b-3b (event bootstrap seeds
  4000 + replicate, 5 replicates, 3a-fix noise, spread s = 1.25, first 60 min,
  knots from 0.35 dB to W_hi). Nothing of 3b-3b was refitted.
- **A** 640 small regressions (seconds: 65 s). **B + C** 25 datasets: the step fit
  (3 starts), three jump-2 formulas, the step law if kept, and 5 menu laws from 4
  starts each, every candidate also on the two halves of the event draws. **D** 9 datasets x
  20 event bootstraps. **E** read from the F4 fits of C.
- **Timing and cut**: replicate 0 of every part was run first: 18.1 min wall on 3
  processes; projected total 1.29 h (budget 3.0 h). Cut:
  nothing.
- **Run time**: 93 min wall (3.2 process-hours B + C, 1.3 D; cached in
  `.cache/fits3b3c/`). This run: 78 s.
- **Software**: pysindy 2.1.0 installed (STLSQ used inside the bagging loop); PySR: 160 fits (15 min in a separate process). `pip install pysindy pysr sympy` worked (pysindy 2.1.0, PySR 2.7.0, scikit-learn 1.9.1, sympy 1.14.0 added
  to the miniconda environment). PySR's Julia backend (Julia 1.11.9, about 1 GB under `miniconda3/julia_env`) failed
  twice because drive C: was full; after C: was cleaned it installed in 11 min (inside the 15 min allowed). PySR
  stopped once on a Windows file lock after 36 fits; the loop now retries and resumed from its cache.

## 2. Part A: jump 2 alone (`research/WAA_estm/outputs/tables/step3b3c_a_jump2_alone.csv`) - PARTIAL

Input: the TRUE h(W) = W k(W) on 40 log-spaced points in [0.1, 2.0] dB,
times (1 + noise), 20 draws per noise level. Library: W, W^2, W^3, sw0.1, sw0.15, sw0.2, sw0.3, sw0.5 (sw x = max(W - x dB, 0)). All
regressions on relative residuals (every point divided by h).
- **best-subset + BIC**: every combination of up to 3 terms; BIC = n log(RSS / n) + terms x log n.
- **E-SINDy**: 100 bootstraps of the points; a term is kept if included in >=
  60% of them; coefficient = median over those. Sparse step: the thresholded backward
  elimination of `waa_sindy.py`, or pysindy's STLSQ (both with relative threshold 0.05).
Rule (judged on best-subset, the main tool; the others are cross-checks): `exact` = exactly the true terms, in >=
16 of 20 draws, and the median coefficient within +-30%. m6: W plus one or both of
sw0.15, sw0.2 (the true A_c = 0.16 dB lies between them; their sum is compared with k_run).

| truth | noise | method | exact | true_included | draws | most_common | coef_ratios | ok |
|---|---|---|---|---|---|---|---|---|
| m1 | 2% | best-subset + BIC | 16 | 20 | 20 | W x16 | W x1.00 | True |
| m1 | 2% | E-SINDy (waa_sindy, bagged) | 20 | 20 | 20 | W x20 | W x1.00 | True |
| m1 | 2% | E-SINDy (pysindy STLSQ, bagged) | 14 | 20 | 20 | W x14 | W x1.00 | False |
| m1 | 5% | best-subset + BIC | 14 | 20 | 20 | W x14 | W x1.00 | False |
| m1 | 5% | E-SINDy (waa_sindy, bagged) | 20 | 20 | 20 | W x20 | W x1.00 | True |
| m1 | 5% | E-SINDy (pysindy STLSQ, bagged) | 0 | 20 | 20 | W + W^2 + W^3 + sw0.1 + sw0.15 + sw0.2 + sw0.3 + sw0.5 x11 | W xnan | False |
| m3 | 2% | best-subset + BIC | 12 | 20 | 20 | W^2 x12 | W^2 x1.00 | False |
| m3 | 2% | E-SINDy (waa_sindy, bagged) | 20 | 20 | 20 | W^2 x20 | W^2 x1.00 | True |
| m3 | 2% | E-SINDy (pysindy STLSQ, bagged) | 20 | 20 | 20 | W^2 x20 | W^2 x1.00 | True |
| m3 | 5% | best-subset + BIC | 15 | 20 | 20 | W^2 x15 | W^2 x0.99 | False |
| m3 | 5% | E-SINDy (waa_sindy, bagged) | 20 | 20 | 20 | W^2 x20 | W^2 x1.00 | True |
| m3 | 5% | E-SINDy (pysindy STLSQ, bagged) | 8 | 20 | 20 | W^2 x8 | W^2 x0.99 | False |
| m5 | 2% | best-subset + BIC | 15 | 20 | 20 | W + W^2 x15 | W x1.00, W^2 x1.00 | False |
| m5 | 2% | E-SINDy (waa_sindy, bagged) | 20 | 20 | 20 | W + W^2 x20 | W x1.00, W^2 x1.00 | True |
| m5 | 2% | E-SINDy (pysindy STLSQ, bagged) | 20 | 20 | 20 | W + W^2 x20 | W x1.01, W^2 x1.00 | True |
| m5 | 5% | best-subset + BIC | 17 | 20 | 20 | W + W^2 x17 | W x1.02, W^2 x0.99 | True |
| m5 | 5% | E-SINDy (waa_sindy, bagged) | 20 | 20 | 20 | W + W^2 x20 | W x1.03, W^2 x0.99 | True |
| m5 | 5% | E-SINDy (pysindy STLSQ, bagged) | 5 | 20 | 20 | W + W^2 x5 | W x1.04, W^2 x1.00 | False |
| m6 | 2% | best-subset + BIC | 20 | 20 | 20 | W + sw0.15 + sw0.2 x20 | W x0.96, switch x1.02 | True |
| m6 | 2% | E-SINDy (waa_sindy, bagged) | 20 | 20 | 20 | W + sw0.15 x20 | W x0.96, switch x0.96 | True |
| m6 | 2% | E-SINDy (pysindy STLSQ, bagged) | 6 | 20 | 20 | W + sw0.15 + sw0.2 x6 | W x1.09, switch x1.12 | False |
| m6 | 5% | best-subset + BIC | 20 | 20 | 20 | W + sw0.15 + sw0.2 x20 | W x0.96, switch x1.02 | True |
| m6 | 5% | E-SINDy (waa_sindy, bagged) | 20 | 20 | 20 | W + sw0.15 x20 | W x0.97, switch x0.95 | True |
| m6 | 5% | E-SINDy (pysindy STLSQ, bagged) | 0 | 20 | 20 | W + W^2 + W^3 + sw0.1 + sw0.15 + sw0.2 + sw0.3 + sw0.5 x12 | W xnan, switch xnan | False |
| m1 | 2% | PySR | 20 | 20 | 20 | W x20 | W x1.00 | True |
| m1 | 5% | PySR | 20 | 20 | 20 | W x20 | W x1.00 | True |
| m3 | 2% | PySR | 20 | 20 | 20 | W^2 x20 | W^2 x1.00 | True |
| m3 | 5% | PySR | 20 | 20 | 20 | W^2 x20 | W^2 x1.00 | True |
| m5 | 2% | PySR | 18 | 18 | 20 | W + W^2 x18 | W x1.01, W^2 x1.00 | True |
| m5 | 5% | PySR | 18 | 18 | 20 | W + W^2 x18 | W x1.02, W^2 x0.99 | True |
| m6 | 2% | PySR | 0 | 0 | 20 | W**2*exp(-W) + W*exp(-W) x2 | W xnan, switch xnan | False |
| m6 | 5% | PySR | 0 | 0 | 20 | W*exp(-0.3041961/W) x1 | W xnan, switch xnan | False |

## 3. Part B: jump 1 with an optional step (`research/WAA_estm/outputs/tables/step3b3c_b_step.csv`) - PARTIAL

k(W) = k_low (a constant) below W_c and the 3b-3b spline above it, joined by a logistic of width 0.02 in
ln W; W_c free in [0.08, 0.6] dB, step height >= 0 (ln units). Fitted from the 3b-3b
solution of the dataset (same knots, same penalty) with W_c starting at [0.12, 0.2, 0.4] dB (250
evaluations each, the best continued for 400). The step costs 2 parameters; it is kept only if
BIC(no step) - BIC(step) > 6. Rule: m6 step kept in >= 4 of 5 and the median W_c of the kept steps
within +-30% of 0.16 dB; m1, m3, m5 no step in >= 4 of 5 and the true k inside the band
of the chosen curves on >= 80% of [0.35 dB, W_hi].

| truth | step_kept | median_bic_gain | coverage_above_0.35 | k_ratio_above_0.35 | coverage_above_0.1 | k_ratio_above_0.1 | ok | median_wc_kept | wc_within_tol |
|---|---|---|---|---|---|---|---|---|---|
| m1 | 5/5 | 108.5 | 0.65 | 1.042 | 0.525 | 0.993 | False |  |  |
| m3 | 5/5 | 37.0 | 1.0 | 0.973 | 0.975 | 0.995 | False |  |  |
| m5 | 5/5 | 43.3 | 0.85 | 1.062 | 0.9 | 1.059 | False |  |  |
| m6 | 5/5 | 116.5 | 1.0 | 0.953 | 0.925 | 1.001 | True | 0.209 | 3.0 |

| truth | rep | bic_no_step | bic_step | bic_gain | step_kept | wc_db | height_ln | k_low | best_wc_start | seconds |
|---|---|---|---|---|---|---|---|---|---|---|
| m1 | 0 | 32750.5 | 32642.0 | 108.5 | True | 0.39 | 0.46 | 0.03146 | 0.4 | 89.2 |
| m1 | 1 | 41266.8 | 41124.3 | 142.5 | True | 0.212 | 0.74 | 0.02562 | 0.2 | 130.1 |
| m1 | 2 | 37264.7 | 37128.5 | 136.1 | True | 0.406 | 0.44 | 0.03131 | 0.4 | 220.9 |
| m1 | 3 | 31449.6 | 31368.3 | 81.3 | True | 0.382 | 0.47 | 0.03258 | 0.4 | 96.4 |
| m1 | 4 | 23990.1 | 23890.7 | 99.4 | True | 0.209 | 1.09 | 0.02135 | 0.2 | 66.4 |
| m3 | 0 | 31940.4 | 31903.9 | 36.6 | True | 0.398 | 0.67 | 0.01942 | 0.4 | 93.6 |
| m3 | 1 | 40590.7 | 40508.4 | 82.3 | True | 0.213 | 1.09 | 0.00468 | 0.2 | 132.1 |
| m3 | 2 | 36899.4 | 36862.4 | 37.0 | True | 0.393 | 0.76 | 0.01379 | 0.4 | 152.7 |
| m3 | 3 | 30755.9 | 30715.2 | 40.7 | True | 0.203 | 1.15 | 0.00357 | 0.2 | 97.7 |
| m3 | 4 | 23505.2 | 23489.9 | 15.3 | True | 0.328 | 0.86 | 0.01335 | 0.4 | 68.0 |
| m5 | 0 | 32166.8 | 32116.4 | 50.4 | True | 0.37 | 0.61 | 0.01481 | 0.4 | 85.1 |
| m5 | 1 | 40483.4 | 40393.5 | 89.9 | True | 0.229 | 0.94 | 0.00611 | 0.2 | 199.7 |
| m5 | 2 | 36884.8 | 36852.8 | 32.1 | True | 0.381 | 0.68 | 0.01541 | 0.4 | 153.6 |
| m5 | 3 | 30805.1 | 30761.8 | 43.3 | True | 0.203 | 1.0 | 0.00635 | 0.2 | 90.3 |
| m5 | 4 | 23546.5 | 23509.9 | 36.6 | True | 0.215 | 1.11 | 0.00637 | 0.2 | 63.3 |
| m6 | 0 | 32679.1 | 32567.0 | 112.1 | True | 0.209 | 1.27 | 0.0083 | 0.2 | 92.4 |
| m6 | 1 | 40693.4 | 40556.6 | 136.8 | True | 0.232 | 1.34 | 0.00544 | 0.2 | 176.7 |
| m6 | 2 | 37055.4 | 36938.9 | 116.5 | True | 0.197 | 1.2 | 0.00634 | 0.2 | 148.8 |
| m6 | 3 | 31352.3 | 31179.9 | 172.5 | True | 0.25 | 2.39 | 0.00643 | 0.4 | 82.8 |
| m6 | 4 | 23507.2 | 23445.4 | 61.8 | True | 0.207 | 1.13 | 0.00621 | 0.2 | 63.8 |

![B m1](../figures/step3b3c_b_step_m1.png)

![B m3](../figures/step3b3c_b_step_m3.png)

![B m5](../figures/step3b3c_b_step_m5.png)

![B m6](../figures/step3b3c_b_step_m6.png)

## 4. Part C: jump 3, refit and choose (`research/WAA_estm/outputs/tables/step3b3c_c_winners.csv`, `research/WAA_estm/outputs/tables/step3b3c_c_candidates.csv`) - PARTIAL

Candidates per dataset: **J2-1, J2-2, J2-3** = the best subset of 1, 2 and 3 terms on the chosen free curve of that
dataset (no step: library W, W^2, W^3, sw0.5, sw0.75 on [0.35 dB, W_hi]; step kept: the library of part A on
[0.1 dB, W_hi]), refitted as a law h(W) of the pooled model; **step** = the step law of B if kept;
**F1, F3, F4, F5, F6** = the menu laws, best of 4 starts (3 generic + the F1 solution; F1: 4 generic). Every candidate
is fitted on one half of the event draws and scored by the log-likelihood of the other half, both ways (the halves
of 3b-3b); winner = highest sum. Families: F1 / "W" = exponential; F3 / F5 / "W^2" / "W + W^2" = quadratic; F6 / step /
"W + switch at <= 0.3 dB" = switch; F4 = two layers; other formulas = other. Rule: winner in the true family in >=
4 of 5.

| truth | winners | families | in_true_family | bic_in_true_family | median_margin | free_curve_minus_winner |
|---|---|---|---|---|---|---|
| m1 | step x5 | switch x5 | 0 | 0 | 44.17 | -51.9 |
| m3 | step x4, J2-3 x1 | switch x4, other x1 | 0 | 0 | 21.47 | -20.04 |
| m4 | F4 x3, J2-3 x1, F6 x1 | two layers x3, other x1, switch x1 | 3 | 3 | 4.22 | -6.21 |
| m5 | step x4, J2-3 x1 | switch x4, other x1 | 0 | 0 | 16.7 | -19.72 |
| m6 | step x5 | switch x5 | 5 | 5 | 60.63 | -61.57 |

`margin_to_other_family` = held-out log-likelihood of the winner minus the best candidate of another family;
`simplest_within_2` = the candidate with the fewest parameters within 2 of the winner; `bic_winner` = lowest BIC on
the whole dataset; `free_curve_minus_winner` = the free curve (reference, not a candidate) minus the winner;
`j2_bic_pick` = what BIC alone would pick on the curve.

| truth | rep | winner | winner_what | family | in_true_family | margin_to_other_family | best_other_family | simplest_within_2 | bic_winner | free_curve_minus_winner | step_kept | j2_bic_pick |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| m1 | 0 | step | free curve with a step at 0.39 dB | switch | False | 55.6 | F4 (two layers) | step | step | -61.06 | True | W + W^2 |
| m1 | 1 | step | free curve with a step at 0.21 dB | switch | False | 41.64 | F5 (quadratic) | step | step | -50.06 | True | W + sw0.15 + sw0.3 |
| m1 | 2 | step | free curve with a step at 0.41 dB | switch | False | 59.36 | F1 (exponential) | step | step | -77.98 | True | W + W^2 + sw0.5 |
| m1 | 3 | step | free curve with a step at 0.38 dB | switch | False | 44.17 | F1 (exponential) | step | step | -51.9 | True | W + sw0.2 |
| m1 | 4 | step | free curve with a step at 0.21 dB | switch | False | 17.19 | F1 (exponential) | step | step | -20.05 | True | W + sw0.2 + sw0.3 |
| m3 | 0 | step | free curve with a step at 0.40 dB | switch | False | 35.32 | J2-1 (quadratic) | step | step | -26.72 | True | W + W^2 + sw0.15 |
| m3 | 1 | step | free curve with a step at 0.21 dB | switch | False | 31.47 | J2-3 (other) | step | step | -29.72 | True | W + W^3 + sw0.2 |
| m3 | 2 | step | free curve with a step at 0.39 dB | switch | False | 6.77 | J2-3 (other) | step | step | -8.9 | True | W + W^3 + sw0.3 |
| m3 | 3 | J2-3 | W + W^3 + sw0.2 | other | False | 14.88 | step (switch) | J2-3 | step | -7.63 | True | W + W^3 + sw0.2 |
| m3 | 4 | step | free curve with a step at 0.33 dB | switch | False | 21.47 | J2-1 (quadratic) | step | step | -20.04 | True | W + W^3 + sw0.3 |
| m4 | 0 | F4 | tau_f 16.36, tau_s 155.7, phi 0.1736 | two layers | True | 4.22 | F6 (switch) | F4 | F4 | -6.21 | False | W + W^2 + W^3 |
| m4 | 1 | J2-3 | W^2 + W^3 + sw0.5 | other | False | 2.45 | F4 (two layers) | J2-3 | J2-3 | 2.57 | False | W^2 + W^3 + sw0.5 |
| m4 | 2 | F4 | tau_f 18.86, tau_s 165.6, phi 0.1971 | two layers | True | 18.21 | F6 (switch) | F4 | F4 | -34.01 | False | W^2 + W^3 + sw0.5 |
| m4 | 3 | F4 | tau_f 19.01, tau_s 212.5, phi 0.1762 | two layers | True | 24.92 | F6 (switch) | F4 | F4 | -21.42 | False | W^2 + W^3 + sw0.5 |
| m4 | 4 | F6 | k_run 0.04938, k_evap 0.004542, A_c 0.1941 | switch | False | 3.52 | F4 (two layers) | F6 | F6 | 2.29 | False | W^2 + W^3 + sw0.5 |
| m5 | 0 | step | free curve with a step at 0.37 dB | switch | False | 16.7 | J2-2 (quadratic) | step | step | -19.72 | True | W + W^2 + sw0.1 |
| m5 | 1 | step | free curve with a step at 0.23 dB | switch | False | 24.56 | J2-2 (other) | step | step | -23.01 | True | W + W^3 + sw0.2 |
| m5 | 2 | step | free curve with a step at 0.38 dB | switch | False | 21.09 | J2-1 (quadratic) | step | step | -36.89 | True | W + W^3 + sw0.3 |
| m5 | 3 | step | free curve with a step at 0.20 dB | switch | False | 13.03 | J2-3 (other) | step | step | -14.25 | True | W + W^3 + sw0.2 |
| m5 | 4 | J2-3 | W + W^3 + sw0.2 | other | False | 12.68 | J2-1 (quadratic) | J2-3 | step | -17.84 | True | W + W^3 + sw0.2 |
| m6 | 0 | step | free curve with a step at 0.21 dB | switch | True | 60.63 | F4 (two layers) | step | step | -68.89 | True | W + sw0.2 + sw0.3 |
| m6 | 1 | step | free curve with a step at 0.23 dB | switch | True | 102.69 | F4 (two layers) | step | step | -61.57 | True | W + sw0.2 |
| m6 | 2 | step | free curve with a step at 0.20 dB | switch | True | 47.16 | F4 (two layers) | step | step | -51.95 | True | W + sw0.15 + sw0.2 |
| m6 | 3 | step | free curve with a step at 0.25 dB | switch | True | 75.41 | F4 (two layers) | step | step | -75.85 | True | W + sw0.15 + sw0.2 |
| m6 | 4 | step | free curve with a step at 0.21 dB | switch | True | 41.66 | F4 (two layers) | step | J2-3 | -39.1 | True | W + sw0.2 + sw0.3 |

Held-out log-likelihood of every candidate minus the best one of its dataset:

| truth | rep | F1 | F3 | F4 | F5 | F6 | J2-1 | J2-2 | J2-3 | free curve (reference) | step |
|---|---|---|---|---|---|---|---|---|---|---|---|
| m1 | 0 | -55.6 | -55.8 | -55.6 | -55.7 | -54.4 | -56.3 | -56.4 | -57.0 | -61.1 | 0.0 |
| m1 | 1 | -47.4 | -44.4 | -43.5 | -41.6 | -35.2 | -43.6 | -43.9 | -42.3 | -50.1 | 0.0 |
| m1 | 2 | -59.4 | -61.5 | -61.0 | -63.0 | -65.9 | -60.8 | -63.1 | -70.0 | -78.0 | 0.0 |
| m1 | 3 | -44.2 | -46.6 | -45.5 | -45.1 | -45.3 | -44.6 | -41.6 | -52.0 | -51.9 | 0.0 |
| m1 | 4 | -17.2 | -19.1 | -19.4 | -19.2 | -17.0 | -17.2 | -21.1 | -30.3 | -20.0 | 0.0 |
| m3 | 0 | -703.7 | -55.7 | -121.7 | -49.9 | -190.8 | -35.3 | -160.3 | -61.4 | -26.7 | 0.0 |
| m3 | 1 | -811.5 | -67.9 | -139.6 | -65.2 | -218.0 | -61.8 | -64.8 | -31.5 | -29.7 | 0.0 |
| m3 | 2 | -748.2 | -40.4 | -82.2 | -26.0 | -118.9 | -15.7 | -124.8 | -6.8 | -8.9 | 0.0 |
| m3 | 3 | -630.2 | -34.7 | -93.9 | -32.3 | -190.3 | -17.1 | -15.1 | 0.0 | -7.6 | -14.9 |
| m3 | 4 | -455.7 | -27.8 | -81.4 | -28.6 | -130.2 | -21.5 | -32.8 | -49.4 | -20.0 | 0.0 |
| m4 | 0 | -130.7 | -28.1 | 0.0 | -68.7 | -4.2 | -126.3 | -62.7 | -44.4 | -6.2 |  |
| m4 | 1 | -215.8 | -60.1 | -2.4 | -98.2 | -26.8 | -216.1 | -95.4 | 0.0 | 2.6 |  |
| m4 | 2 | -190.3 | -50.7 | 0.0 | -108.0 | -18.2 | -195.1 | -108.2 | -69.1 | -34.0 |  |
| m4 | 3 | -141.3 | -56.3 | 0.0 | -88.0 | -24.9 | -141.4 | -85.0 | -492.5 | -21.4 |  |
| m4 | 4 | -110.3 | -43.4 | -3.5 | -77.7 | 0.0 | -108.9 | -75.2 | -45.5 | 2.3 |  |
| m5 | 0 | -587.9 | -25.6 | -86.9 | -22.1 | -160.6 | -31.8 | -16.7 | -18.1 | -19.7 | 0.0 |
| m5 | 1 | -719.0 | -51.6 | -116.5 | -50.2 | -191.2 | -32.3 | -24.6 | -33.3 | -23.0 | 0.0 |
| m5 | 2 | -695.2 | -49.7 | -87.4 | -44.4 | -145.7 | -21.1 | -21.3 | -41.9 | -36.9 | 0.0 |
| m5 | 3 | -582.1 | -41.1 | -90.7 | -40.6 | -149.7 | -37.8 | -28.5 | -13.0 | -14.2 | 0.0 |
| m5 | 4 | -424.3 | -32.2 | -84.3 | -33.8 | -101.8 | -12.7 | -32.3 | 0.0 | -17.8 | -15.7 |
| m6 | 0 | -125.3 | -91.4 | -60.6 | -112.8 | -46.8 | -61.7 | -52.9 | -47.9 | -68.9 | 0.0 |
| m6 | 1 | -204.8 | -130.1 | -102.7 | -176.4 | -33.4 | -395.7 | -34.8 | -26.1 | -61.6 | 0.0 |
| m6 | 2 | -169.2 | -103.3 | -47.2 | -146.9 | -24.6 | -86.8 | -24.2 | -28.7 | -52.0 | 0.0 |
| m6 | 3 | -174.4 | -148.7 | -75.4 | -169.6 | -68.3 | -492.4 | -71.2 | -68.2 | -75.8 | 0.0 |
| m6 | 4 | -112.0 | -93.1 | -41.7 | -109.5 | -20.5 | -307.0 | -19.3 | -11.8 | -39.1 | 0.0 |

**Post hoc, added after the first complete run (no refit; no verdict uses it): the menu laws only.** The step law
wins for every truth (see part B), and the jump-2 formulas were read from the curve with the step. The same held-out
scores, restricted to F1, F3, F4, F5, F6:

| truth | held_out_winners | in_true_family | median_margin | bic_winners | bic_in_true_family |
|---|---|---|---|---|---|
| m1 | F6 x3, F1 x2 | 2 | 1.2 | F1 x4, F6 x1 | 4 |
| m3 | F5 x4, F3 x1 | 5 | 61.5 | F5 x3, F3 x2 | 5 |
| m4 | F4 x4, F6 x1 | 4 | 18.2 | F4 x4, F6 x1 | 4 |
| m5 | F5 x4, F3 x1 | 5 | 52.0 | F5 x3, F3 x2 | 5 |
| m6 | F6 x5 | 5 | 21.2 | F6 x5 | 5 |

## 5. Part D: the band for real data (`research/WAA_estm/outputs/tables/step3b3c_d_bootstrap_band.csv`) - FAIL

One dataset = replicate 0, 1, 2 of 3b-3b. Its event draws are resampled with replacement
20 times (same noise, same tails per event), the free curve is refitted on each (same knots and penalty,
start = the fit of the dataset, 400 evaluations). Band = 5-95 % of the refits at every W.
Rule: true k inside that band on >= 80% of [0.35 dB, W_hi] in >= 2 of
3 datasets per truth. `band_half_width_pct` = median half-width of the band (+- %);
`slope_lo`-`slope_hi` = the same percentiles of the slope of log k vs log W.

| truth | dataset | w_hi_db | bootstraps | coverage_band | coverage_minmax | band_half_width_pct | minmax_half_width_pct | fit_k_ratio | slope_fit | slope_lo | slope_hi | true_slope | seconds | ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| m1 | 0 | 3.523 | 20 | 0.775 | 0.875 | 3.4 | 5.0 | 0.965 | 0.01 | -0.03 | 0.05 | 0.0 | 388.9 | False |
| m1 | 1 | 2.466 | 20 | 0.525 | 0.75 | 3.4 | 5.3 | 0.971 | 0.06 | 0.03 | 0.11 | 0.0 | 545.4 | False |
| m1 | 2 | 2.407 | 20 | 0.3 | 0.375 | 3.1 | 3.7 | 0.964 | 0.04 | 0.03 | 0.1 | 0.0 | 456.3 | False |
| m3 | 0 | 2.631 | 20 | 0.0 | 0.0 | 4.7 | 6.2 | 1.231 | 1.04 | 1.01 | 1.09 | 1.0 | 324.2 | False |
| m3 | 1 | 2.111 | 20 | 0.275 | 0.45 | 1.3 | 2.7 | 0.975 | 1.04 | 1.02 | 1.04 | 1.0 | 612.7 | False |
| m3 | 2 | 2.026 | 20 | 0.75 | 0.85 | 6.0 | 7.7 | 0.979 | 1.09 | 1.03 | 1.12 | 1.0 | 668.7 | False |
| m5 | 0 | 2.74 | 20 | 0.65 | 0.675 | 6.3 | 8.3 | 0.928 | 0.95 | 0.91 | 1.03 | 0.93 | 323.0 | False |
| m5 | 1 | 2.156 | 20 | 0.675 | 0.95 | 5.1 | 6.9 | 1.008 | 0.97 | 0.89 | 1.03 | 0.92 | 487.0 | False |
| m5 | 2 | 2.079 | 20 | 0.275 | 0.375 | 3.2 | 4.6 | 1.038 | 1.03 | 0.99 | 1.07 | 0.92 | 745.0 | False |

![D m1](../figures/step3b3c_d_band_m1.png)

![D m3](../figures/step3b3c_d_band_m3.png)

![D m5](../figures/step3b3c_d_band_m5.png)

## 6. Part E: F4 check for the m4 warning of 3b-3b (`research/WAA_estm/outputs/tables/step3b3c_e_f4_check.csv`) - information only

3b-3b fitted F4 from 3 generic starts and found BIC(free curve) - BIC(F4) > 10 in 2 of 5. Here F4 also starts from
the F1 solution of the dataset. The 3b-3b verdict stays as written.

| rep | bic_free_curve | bic_f4_3b3b | bic_f4_generic_now | bic_f4_f1sol | best_start | curve_minus_f4_3b3b | curve_minus_f4_now | warning_now |
|---|---|---|---|---|---|---|---|---|
| 0 | 32462.1 | 32459.4 | 32459.4 | 32460.0 | g3 | 2.7 | 2.7 | False |
| 1 | 40712.4 | 40708.2 | 40708.2 | 40709.3 | g3 | 4.2 | 4.2 | False |
| 2 | 36876.8 | 36838.1 | 36838.1 | 36836.9 | f1sol | 38.7 | 39.9 | True |
| 3 | 31403.2 | 31350.0 | 31350.0 | 31349.9 | f1sol | 53.2 | 53.3 | True |
| 4 | 23305.8 | 23302.5 | 23302.5 | 23303.4 | g3 | 3.3 | 3.3 | False |

With the extra start: > 10 in 2 of 5.

## 7. The Step 4 recipe that follows

| part | truth | verdict | basis |
|---|---|---|---|
| A | m1 | PARTIAL | 2%: 16/20 (W x1.00); 5%: 14/20 (W x1.00) |
| A | m3 | FAIL | 2%: 12/20 (W^2 x1.00); 5%: 15/20 (W^2 x0.99) |
| A | m5 | PARTIAL | 2%: 15/20 (W x1.00, W^2 x1.00); 5%: 17/20 (W x1.02, W^2 x0.99) |
| A | m6 | PASS | 2%: 20/20 (W x0.96, switch x1.02); 5%: 20/20 (W x0.96, switch x1.02) |
| B | m1 | FAIL | no step in 0 of 5; coverage above 0.35 dB 65% |
| B | m3 | FAIL | no step in 0 of 5; coverage above 0.35 dB 100% |
| B | m5 | FAIL | no step in 0 of 5; coverage above 0.35 dB 85% |
| B | m6 | PASS | step kept in 5 of 5; median W_c 0.209 dB (true 0.16; within +-30% in 3) |
| C | m1 | FAIL | held-out winner in the true family (exponential) in 0 of 5: step x5 |
| C | m3 | FAIL | held-out winner in the true family (quadratic) in 0 of 5: step x4, J2-3 x1 |
| C | m4 | FAIL | held-out winner in the true family (two layers) in 3 of 5: F4 x3, J2-3 x1, F6 x1 |
| C | m5 | FAIL | held-out winner in the true family (quadratic) in 0 of 5: step x4, J2-3 x1 |
| C | m6 | PASS | held-out winner in the true family (switch) in 5 of 5: step x5 |
| D | m1 | FAIL | true k inside the bootstrap band on >= 80 % in 0 of 3 datasets (78%, 52%, 30%) |
| D | m3 | FAIL | true k inside the bootstrap band on >= 80 % in 0 of 3 datasets (0%, 28%, 75%) |
| D | m5 | FAIL | true k inside the bootstrap band on >= 80 % in 0 of 3 datasets (65%, 68%, 28%) |

**Read the table with care: the PASS cells of m6 in B and C are not evidence.** The step is kept, and the step
law wins jump 3, for EVERY truth, also the exponential one. A tool that says "switch" for everything has not found
the switch of m6.

What each part showed:
- **A (jump 2 alone).** Given the true curve with 2-5 % noise, all three working tools find the true terms; they
  differ in how often they add one. Best-subset + BIC: true terms present in 20 of 20 everywhere, but one extra
  small term in 3 to 8 of 20 draws for m1, m3, m5 (so 5 of 8 cells miss the 16-of-20 rule); m6 exact in 20 of 20 with
  W x0.96 and switch x1.02. E-SINDy with the backward elimination of `waa_sindy.py`: exact in 20 of 20 in all 8
  cells, m6 included. PySR: exact in 20 of 20 for m1 and m3, 18 of 20 for m5, coefficients within 2 %; it cannot
  write m6 (no max in its operators) and returns a different smooth formula in every draw. pysindy's STLSQ inside
  the same bagging fails at 5 % noise (keeps all 8 terms): plain thresholding cannot handle the switch columns.
- **B (step).** The step is kept in 5 of 5 for m1, m3, m5 and m6, with BIC gains of 15 to 170. For m1 the fitted
  curve is flat below W_c, jumps up by 30-60 % in a narrow spike just above W_c (0.2-0.4 dB) and comes back to the
  true level. That is not in the truth. W_c sits near 0.2 or 0.4 dB for every truth (m6: median 0.21, true 0.16).
- **C (jump 3).** The step law has the best HELD-OUT log-likelihood in 18 of the 20 single-curve datasets, by 16 to
  65 units over the true menu law for m1, m3, m5. So this is not overfitting to one half of the events: the pooled
  likelihood itself prefers a feature near 0.2-0.4 dB that the truth does not have. m4: F4 wins in 3 of 5 (rule: 4).
  Post hoc, menu laws only: BIC on the whole dataset names the true family in 4, 5, 4, 5, 5 of 5 (m1, m3, m4,
  m5, m6), as in 3b-1c; the held-out score does the same except for m1 (F6 beats F1 in 3 of 5, by about 1 unit).
- **D (band from one dataset).** The 5-95 % event-bootstrap band is +-1 to 6 % wide and holds the true k on 0 to
  78 % of the range; no dataset reaches 80 %. The fits themselves are off by -7 to +23 % in level. The band is too
  narrow by a factor of about 3. The slope is much steadier: fitted slopes are within 0.01-0.11 of the truth.
- **E (F4 check).** The F1-solution start changes BIC(F4) by at most 1.2. F4 was at its optimum in 3b-3b; the m4
  warning stays at 2 of 5.

**The Step 4 recipe that follows (only what passed, here or in 3b-3b):**
1. **Jump 1**: the free-curve pooled fit of 3b-3b, unchanged: first 60 min, 5 knots from 0.35 dB to W_hi, read only
   for W >= 0.35 dB. **No step**: the step tool is not used on real tails.
2. **Shape**: the slope of log k vs log W over [0.35 dB, W_hi]. About 0 = exponential there, about 1 = quadratic,
   0.2-0.5 = "mildly rising, law not identified by the curve" (3b-3b).
3. **Band**: none passed. Show the event-bootstrap band, but say next to it that on fake data it was about 3
   times too narrow. For the level use the spread seen on fake replicates (about +-20 %); for the slope about
   +-0.1. Do not write "the law is inside the band".
4. **Jump 2** on the curve above 0.35 dB (library W, W^2, W^3, switch at 0.5 / 0.75 dB): E-SINDy with the
   `waa_sindy.py` elimination as the main tool, PySR as cross-check, best-subset + BIC as the list of "true terms
   plus at most one extra". Not pysindy's STLSQ.
5. **Jump 3 is not used to choose the law.** Held-out likelihood rewards whatever bends the curve near
   0.2-0.4 dB. For the comparison of named laws use the menu fits F1, F3, F4, F5, F6 with BIC on the whole
   dataset and the 3b-1c rules (several starts, F3 / F5 one family, F4 vs F6 under about 10 BIC not decisive).
6. **Before any statement about a switch or about W < 0.35 dB**: repair the likelihood at low W and repeat B, C
   and D on the cached fake datasets (minutes to an hour each). See unclear point 1.

## 8. In plain words: what we will be able to say about the real drying law

- **For W between 0.35 dB and W_hi (about 2 dB on fake data; to be computed on the real tails), first 60 min:** we can
  say how the drying rate per dB of water changes with W, as a slope with an uncertainty of about +-0.1, and
  whether the curve there is a one-term law. "Exponential in this range" (slope near 0, equation W) and
  "quadratic in this range" (slope near 1, equation W^2) would both be recovered; this was confirmed on fresh
  fake data (3b-3b) and by jump 2 alone here.
- **The level of the rate (how fast, not which shape)** is known to about +-20 %, not to the +-5 % the bootstrap
  band suggests.
- **A mild rise (slope 0.2-0.5)** can be stated as such, but not explained: a switch and two layers look the same
  above 0.35 dB.
- **For W below 0.35 dB we can say nothing yet**, and in particular nothing about a switch level. The tool built
  for that draws a step for every truth, including a pure exponential. A step found on real tails would
  therefore mean nothing.
- **Which named law is best** can be answered by the menu fits with BIC (true family in 4-5 of 5 on fake data,
  F3 / F5 as one family), not by held-out likelihood and not by the free curve against F4.
- So the honest real-data sentence will be of the form: "between 0.35 and about 2 dB, in the first hour after
  rain, the drying rate per dB of water [is constant / rises with slope x +- 0.1]; among our candidate laws the
  data prefer family Y by BIC; below 0.35 dB the method cannot yet tell".

## 9. Unclear / possibly wrong

1. **Why every truth gets a step: my reading is the 0.3 dB quantization, and I did not test it.** The simulator
   quantizes RSL in 0.3 dB steps; the likelihood treats the noise as Gaussian. The false feature sits at 0.2-0.4 dB,
   it has the shape "flat, then a fast passage" (what a signal crossing one quantization step looks like), it wins
   on held-out events, and 3b-3 already saw the free curve beat the true law at the low end. All of that fits, none
   of it is proof. Other suspects: the offset prior, the start-level floor (0.05 dB), the per-tail start profile. A
   direct test: simulate m1 without quantization and refit the step (about 15 min). Real RSL is quantized the
   same way, so the problem is not confined to fake data.
2. **If point 1 is right, the menu fits are touched too.** F6 has its switch at 0.16 dB and beats F1 on m1 data
   in 3 of 5 held-out comparisons (by about 1 unit; by 12 in one). With BIC the true family still wins in 4-5
   of 5, so the damage is small there, but F6's A_c on real tails should not be read as a physical level.
3. **The step model is harsher than it had to be**: a constant below W_c and a logistic only 0.02 wide in ln W
   (my reading of "width 0.02 dB in log W"), narrower than one cell of the internal table (0.027). The narrow
   spike above W_c is partly the spline reacting to that jump. A gentler step might behave differently; not tried.
4. **m6's W_c**: the median is 0.21 dB, inside +-30 % of 0.16 dB, but only 3 of the 5 replicates are inside it, and
   m1 gives 0.21 dB in 2 of 5. Judging the median of the kept steps was my choice; judged per replicate the cell
   would fail (3 of 5, rule 4). Either way the cell carries no information (see section 7).
5. **Jump 3 candidates were my construction**: instead of one jump-2 answer I sent the best subset of 1, 2 and 3
   terms, because BIC on a smooth fitted curve always takes the largest formula (column `j2_bic_pick`). And when
   the step was kept (always), the formulas were read from the curve WITH the step on [0.1 dB, W_hi], so they
   inherit the artefact. Jump 3 on formulas read above 0.35 dB only was therefore never tested.
6. **Half-fits start from the full-data solution** (menu laws and step: 250 evaluations; formulas: 200). The
   fitted half has thus "seen" the other half through its start. It affects all candidates alike, but the
   held-out scores are slightly optimistic.
7. **The bootstrap band may be too narrow for two different reasons**, and I cannot separate them: (a) resampling
   event draws of one dataset keeps the same noise and the same link speed factors, which is what a real
   bootstrap also does; (b) each refit starts at the dataset's solution and gets 400 evaluations, so it may not
   move far enough. (b) is a weakness of my setup; a refit from the four starts would cost 4 times more.
8. **The band was judged up to each dataset's own W_hi** (2.0-3.5 dB), where the top knot has only 200
   tail-minutes; 3b-3b judged up to the smallest W_hi of the replicates. That makes D somewhat harsher. The level
   errors (x0.93 to x1.23) are medians over the whole range, though.
9. **"+-20 % for the level, +-0.1 for the slope" in the recipe are read off 9 fits**, not derived. They are a
   placeholder until a band passes.
10. **Part A is an easy test**: independent noise on the true curve. A fitted curve has smooth, correlated
    errors, which is why best-subset + BIC cannot be used on it directly. E-SINDy's 20 of 20 used one fixed
    threshold (0.05) that I chose; no other value was tried.
11. **PySR**: judged by splitting its formula into terms with sympy (terms under 2 % of h dropped, exponents
    within 0.1 of an integer rounded). The m5 misses are formulas like W^1.78 exp(0.12 W), which fit as well
    but are not "W + W^2". 30 iterations, about 5 s per fit.
12. **The post hoc "menu laws only" table** was added after seeing the results; it reuses the stored scores and
    the decision table does not use it.
13. **Cut order**: E cost nothing extra (C already fits F4 with the F1-solution start), so the first real cut
    would have been the bootstraps of D (20 -> 12, my reading of "D-replicates 5 -> 3"). Nothing was cut:
    projected 1.3 h, run 1.6 h, plus 15 min of PySR in a fourth process.
14. **Installed software stays on the machine** (about 1.2 GB on C:, mostly Julia). A temporary folder of my smoke
    test is left in `.cache/` (`tmp*`, safe to delete).
15. **Key figures**: none proposed. The step and band figures show tools that failed.
