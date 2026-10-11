# WAA_estm step 3b-3: "draw the curve, then name it" on fake tails (known truth)

Written by `research/WAA_estm/step3b3_curve_then_equation.py` (libraries `waa_curve.py`, `waa_pooled.py`, `waa_sindy.py`,
`waa_sim.py`; settings `waa.curve3b3`). Fake data only: no real tail is fitted, the frozen tail set v1 is not changed.

**Overall verdict: PARTIAL.** Per truth: m1 FAIL, m3 PARTIAL, m4 PASS, m5 PARTIAL, m6 PARTIAL.

## 1. What was run

- **Data**: the fake datasets of steps 3b-1c / 3b-2 (event bootstrap of the 429-tail TRAIN v1 skeleton, seeds
  3000 + replicate), 5 replicates. Main test: m1, m3, m4, m5, m6 with spread s = 1.25 and 3a-fix
  noise. Extras: no spread, 3a-fix noise (m1), spread 1.25, T1 noise (m1, m3).
  Windows: first_60min, first_6h (the extras: first 60 min only).
- **Jobs**: 65 free-curve pooled fits (dataset x window), each = 4 starts + 3
  penalties + 2-fold cross-validation + curvature; 72553 likelihood evaluations in the main fits.
- **Cut for time**: 6 h window of the extras (m1 no spread, T1 noise): they keep the first 60 min. Timing: 2.94
  microseconds per tail-minute and evaluation on one replicate (m6), 2053 evaluations planned per job;
  projected total 2.55 h on 3 processes (budget 3.0 h).
- **Run time**: 127 min wall time (6.3 process-hours; cached in
  `.cache/fits3b3/`). This run: 11 s.

| kind | window | jobs | tails | s_median | converged | with_interval | start_spread_median | starts_within_2_median | minutes_median |
|---|---|---|---|---|---|---|---|---|---|
| T1 | first_60min | 10 | 478.0 | 1.1212499999999999 | 10 | 0 | 39.845 | 1.0 | 2.4 |
| main | first_60min | 25 | 478.0 | 1.2545 | 25 | 1 | 33.97 | 1.0 | 2.3 |
| main | first_6h | 25 | 478.0 | 1.241 | 25 | 1 | 41.35 | 1.0 | 10.8 |
| nospread | first_60min | 5 | 478.0 | 0.0377 | 5 | 1 | 155.58 | 1.0 | 2.5 |

## 2. Spline, penalty, supported range, starts

- **Law**: dW/dt = -f_i W k(W); log k = natural cubic spline in log W with 5 knots (parameters: log k at the
  knots), straight-line continuation outside the knots; f_i = exp(s g), g ~ N(0, 1) (median f = 1), marginalised as
  in 3b-1b; s fitted. Likelihood, offset prior and noise model: exactly those of the menu fits (`waa_pooled.py`;
  checked: a flat curve gives the same log-likelihood as F1 to 0.01).
- **Supported range of W**: from 0.05 dB (fixed: the smallest switch level of jump 2; below the 1 min noise
  sd of about 0.14 dB the curve is held by the offset prior and the early part of the tails, not seen directly) up to
  the level that has 200 tail-minutes of A above it in the first 60 min: 1.52-3.30 dB
  over the datasets. Bands, tests and jump 2 use the range common to the replicates of a truth, on 40
  log-spaced points.
- **Knots**: the two ends of that range and the 3 equal quantiles of log A in between (tail-minutes of the
  first 60 min); the same knots for both windows of a dataset. Example (m6, replicate 0): 0.05 0.304 0.591 1.03 2.31 dB.
- **Penalty**: lambda * integral of (d2 log k / d(log W)2)2, lambda from [0.01, 0.1, 1.0] by 2-fold cross-validation
  over EVENT DRAWS (fit on one half of the event draws, log-likelihood of the other half, both ways; highest sum).
  Plus a weak prior on log k at every knot (sd 3.0 around log(1/30 min)). Chosen penalties (jobs):

| kind | window | 0.01 | 0.1 | 1.0 |
|---|---|---|---|---|
| T1 | first_60min | 2 | 1 | 7 |
| main | first_60min | 11 | 4 | 10 |
| main | first_6h | 6 | 8 | 11 |
| nospread | first_60min | 3 | 1 | 1 |

- **Starts**: flat, rising, falling, step (flat = exponential), 150 evaluations each with lambda =
  0.1; the best is continued (450 more). `start_spread` = worst minus best penalised
  cost after the first stage; `starts_within_2` = starts within 2 of the best. Best start (jobs):

| kind | window | falling | flat | rising | step |
|---|---|---|---|---|---|
| T1 | first_60min | 6 | 2 | 2 | 0 |
| main | first_60min | 5 | 9 | 11 | 0 |
| main | first_6h | 8 | 6 | 10 | 1 |
| nospread | first_60min | 0 | 1 | 1 | 3 |

- **Uncertainty**: the band of a truth = min-max of its 5 replicates (each a new event bootstrap). Per replicate: a
  1.96 sd interval from the curvature of the penalised fit (finite differences); `with_interval` = fits where that
  curvature was positive definite.

## 3. Jump 1: is the curve right? (`research/WAA_estm/outputs/tables/step3b3_j1_curve.csv`)

`band_coverage` = share of the supported range where the true k(W) lies inside the replicate band (rule: >=
80%); `median_k_ratio` = learned / true k (median over range and replicates); shape: slope of log k
vs log W over the range with its interval (flat = interval contains 0; rising = interval above 0), m6: mean k above
2 A_c / mean k below A_c > 2.0 (rule: in >= 4 of 5).

| truth | window | band_coverage | curve_ok | median_k_ratio | median_slope | with_interval | flat_reps | rising_reps | step_reps | median_step_ratio | shape_wanted | shape_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| m1 | first_60min | 0.675 | False | 1.073 | -0.243 | 1 | 1 | 0 |  |  | flat | False |
| m1 | first_6h | 0.475 | False | 1.075 | -0.239 | 0 | 0 | 0 |  |  | flat | False |
| m3 | first_60min | 0.725 | False | 0.975 | 0.964 | 0 | 0 | 0 |  |  | rising | False |
| m3 | first_6h | 0.8 | True | 1.016 | 0.444 | 0 | 0 | 0 |  |  | rising | False |
| m4 | first_60min |  |  |  | 0.755 | 0 | 0 | 0 |  |  | - |  |
| m4 | first_6h |  |  |  | 0.816 | 1 | 0 | 1 |  |  | - |  |
| m5 | first_60min | 0.675 | False | 0.958 | 0.453 | 0 | 0 | 0 |  |  | rising | False |
| m5 | first_6h | 0.925 | True | 1.106 | 0.51 | 0 | 0 | 0 |  |  | rising | False |
| m6 | first_60min | 0.8 | True | 1.058 | 0.87 | 0 | 0 | 0 | 4.0 | 7.71 | step | True |
| m6 | first_6h | 0.6 | False | 1.113 | 0.881 | 0 | 0 | 0 | 4.0 | 7.93 | step | True |

![k m1](../figures/step3b3_k_m1.png)

![k m3](../figures/step3b3_k_m3.png)

![k m4](../figures/step3b3_k_m4.png)

![k m5](../figures/step3b3_k_m5.png)

![k m6](../figures/step3b3_k_m6.png)

**Post hoc diagnostic, added after the first complete run (it changes no verdict; `research/WAA_estm/outputs/tables/step3b3_posthoc_upper_range.csv`).**
Two things went wrong with the judging itself: the per-replicate interval exists for only
3 of 65 fits (the curvature matrix was not positive definite), so every slope test
fails for lack of an interval; and the first knot (0.05 dB) is badly determined while the next knot is at
0.30-0.40 dB. Slope of log k vs log W per replicate over the WHOLE range (min / median / max):

| kind | truth | window | min | median | max |
|---|---|---|---|---|---|
| T1 | m1 | first_60min | -0.36 | -0.32 | -0.19 |
| T1 | m3 | first_60min | 0.08 | 0.49 | 1.04 |
| main | m1 | first_60min | -0.34 | -0.24 | 0.0 |
| main | m1 | first_6h | -0.28 | -0.24 | -0.2 |
| main | m3 | first_60min | 0.48 | 0.96 | 0.99 |
| main | m3 | first_6h | 0.28 | 0.44 | 1.02 |
| main | m4 | first_60min | 0.26 | 0.76 | 0.91 |
| main | m4 | first_6h | 0.09 | 0.82 | 0.92 |
| main | m5 | first_60min | 0.07 | 0.45 | 0.9 |
| main | m5 | first_6h | 0.18 | 0.51 | 0.94 |
| main | m6 | first_60min | 0.28 | 0.87 | 0.94 |
| main | m6 | first_6h | 0.14 | 0.88 | 0.9 |
| nospread | m1 | first_60min | -0.12 | -0.03 | 0.01 |

The same fits judged only ABOVE THE SECOND KNOT (band coverage, level, slope over the replicates against the true
slope there, and the equation jump 2 gives on that range):

| kind | truth | window | w_from_db | w_to_db | band_coverage | median_k_ratio | true_slope | slope_min | slope_median | slope_max | equation | coef_ratios |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| main | m1 | first_60min | 0.397 | 1.779 | 1.0 | 1.033 | 0.0 | -0.04 | -0.02 | 0.1 | W x5 | W x1.03 |
| main | m1 | first_6h | 0.397 | 1.779 | 0.625 | 1.024 | 0.0 | -0.1 | -0.02 | -0.01 | W + W^2 + W^3 + sw0.05 x4, W^2 + sw0.05 x1 | W x-4.05 |
| main | m3 | first_60min | 0.348 | 1.525 | 0.675 | 0.945 | 1.0 | 0.96 | 1.08 | 1.19 | W^2 x5 | W^2 x0.94 |
| main | m3 | first_6h | 0.348 | 1.525 | 0.875 | 0.935 | 1.0 | 0.62 | 0.97 | 1.02 | W^2 x5 | W^2 x0.94 |
| main | m4 | first_60min | 0.362 | 1.725 |  |  |  | 0.3 | 0.37 | 0.42 | W + W^2 + W^3 x2, W + sw0.05 x1, W + W^2 + W^3 + sw0.05 x1, W + W^2 x1 |  |
| main | m4 | first_6h | 0.362 | 1.725 |  |  |  | 0.33 | 0.41 | 0.45 | W + W^2 x3, W + sw0.05 x1, W^2 + W^3 x1 |  |
| main | m5 | first_60min | 0.364 | 1.549 | 0.6 | 0.962 | 0.92 | 0.94 | 0.96 | 0.98 | W^2 x5 | W^2 x1.06 |
| main | m5 | first_6h | 0.364 | 1.549 | 1.0 | 1.047 | 0.92 | 0.75 | 0.89 | 0.95 | W^2 x5 | W^2 x1.14 |
| main | m6 | first_60min | 0.348 | 1.723 | 1.0 | 1.019 | 0.26 | 0.23 | 0.27 | 0.3 | W + W^2 x4, W + sw0.05 x1 | W x5.32, switch x0.00 |
| main | m6 | first_6h | 0.348 | 1.723 | 0.675 | 1.105 | 0.26 | 0.12 | 0.27 | 0.29 | W + W^2 x2, W^2 + W^3 x1, W + W^2 + sw0.05 x1, W + W^2 + W^3 x1 | W x3.11, switch x0.00 |
| nospread | m1 | first_60min | 0.357 | 1.743 | 0.875 | 1.016 | -0.0 | -0.07 | -0.06 | -0.02 | W x5 | W x1.02 |
| T1 | m1 | first_60min | 0.392 | 1.83 | 0.6 | 1.029 | -0.0 | -0.08 | -0.02 | -0.0 | W + W^2 + W^3 + sw0.05 x2, W x1, W^2 + sw0.05 x1, W + W^2 + sw0.05 x1 | W x-2.22 |
| T1 | m3 | first_60min | 0.348 | 1.58 | 0.575 | 0.954 | 1.0 | 0.84 | 0.99 | 1.08 | W^2 x5 | W^2 x0.96 |

## 4. Jump 2: is the equation right? (`research/WAA_estm/outputs/tables/step3b3_j2_selection.csv`, `research/WAA_estm/outputs/tables/step3b3_j2_coefficients.csv`)

h(W) = W k(W) of every replicate on the grid, weights 1 / variance over the replicates, library W, W^2, W^3, sw0.05, sw0.1, sw0.15, sw0.2, sw0.3, sw0.5
(sw x = max(W - x dB, 0)), thresholded backward elimination of `waa_sindy.py`; threshold by leave-one-replicate-out
cross-validation (one-SE rule). Replicates in which a term is selected (**bold** = a term of the true law):

| truth | window | W | W^2 | W^3 | sw0.05 | sw0.1 | sw0.15 | sw0.2 | sw0.3 | sw0.5 |
|---|---|---|---|---|---|---|---|---|---|---|
| m1 | first_60min | **5/5** | 4/5 | 4/5 | 4/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 |
| m1 | first_6h | **5/5** | 5/5 | 2/5 | 5/5 | 1/5 | 0/5 | 0/5 | 4/5 | 0/5 |
| m3 | first_60min | 0/5 | **5/5** | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 |
| m3 | first_6h | 3/5 | **5/5** | 1/5 | 3/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 |
| m4 | first_60min | 0/5 | 5/5 | 5/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 |
| m4 | first_6h | 1/5 | 5/5 | 4/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 |
| m5 | first_60min | 3/5 | **5/5** | 1/5 | 1/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 |
| m5 | first_6h | 4/5 | **5/5** | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 |
| m6 | first_60min | **1/5** | 5/5 | 4/5 | 0/5 | **0/5** | **0/5** | **0/5** | 0/5 | 0/5 |
| m6 | first_6h | **0/5** | 5/5 | 5/5 | 0/5 | **0/5** | **0/5** | **0/5** | 0/5 | 0/5 |

| truth | window | threshold | returned | true_selected | spurious | coef_ratios | equation_ok |
|---|---|---|---|---|---|---|---|
| m1 | first_60min | 0.03 | W + W^2 + W^3 + sw0.05 x3, W + sw0.05 x1, W + W^2 + W^3 x1 | W 5/5 | W^2 4/5, W^3 4/5, sw0.05 4/5 | W x2.10 | False |
| m1 | first_6h | 0.02 | W + W^2 + sw0.05 + sw0.3 x3, W + W^2 + W^3 + sw0.05 + sw0.1 + sw0.3 x1, W + W^2 + W^3 + sw0.05 x1 | W 5/5 | W^2 5/5, W^3 2/5, sw0.05 5/5, sw0.1 1/5, sw0.3 4/5 | W x2.66 | False |
| m3 | first_60min | 0.45 | W^2 x5 | W^2 5/5 |  | W^2 x0.94 | True |
| m3 | first_6h | 0.05 | W + W^2 + sw0.05 x2, W^2 x2, W + W^2 + W^3 + sw0.05 x1 | W^2 5/5 | W 3/5, W^3 1/5, sw0.05 3/5 | W^2 x1.01 | False |
| m4 | first_60min | 0.12 | W^2 + W^3 x5 |  |  |  |  |
| m4 | first_6h | 0.08 | W^2 + W^3 x4, W + W^2 x1 |  |  |  |  |
| m5 | first_60min | 0.05 | W + W^2 x2, W^2 x2, W + W^2 + W^3 + sw0.05 x1 | W^2 5/5 | W^3 1/5, sw0.05 1/5 | W^2 x1.07 | True |
| m5 | first_6h | 0.12 | W + W^2 x4, W^2 x1 | W^2 5/5 |  | W^2 x0.97 | True |
| m6 | first_60min | 0.2 | W^2 + W^3 x4, W + W^2 x1 | W 1/5; sw0.1/sw0.15/sw0.2 0/5 | W^2 5/5, W^3 4/5 | W x0.00, switch x0.00 | False |
| m6 | first_6h | 0.12 | W^2 + W^3 x5 | W 0/5; sw0.1/sw0.15/sw0.2 0/5 | W^2 5/5, W^3 5/5 | W x0.00, switch x0.00 | False |

Coefficients, fitted / true (0 where not selected):

| truth | window | coefficient | true | median_ratio | min_ratio | max_ratio | replicates_nonzero | counted |
|---|---|---|---|---|---|---|---|---|
| m1 | first_60min | W | 0.03361 | 2.104 | 1.203 | 3.319 | 5 | True |
| m1 | first_6h | W | 0.03361 | 2.659 | 2.226 | 3.641 | 5 | True |
| m3 | first_60min | W^2 | 0.05296 | 0.942 | 0.867 | 0.985 | 5 | True |
| m3 | first_6h | W^2 | 0.05296 | 1.009 | 0.965 | 1.833 | 5 | True |
| m5 | first_60min | W | 0.00299 | 0.788 | 0.0 | 20.126 | 3 | False |
| m5 | first_60min | W^2 | 0.04499 | 1.07 | 0.859 | 2.826 | 5 | True |
| m5 | first_6h | W | 0.00299 | 1.773 | 0.0 | 2.315 | 4 | False |
| m5 | first_6h | W^2 | 0.04499 | 0.968 | 0.846 | 1.193 | 5 | True |
| m6 | first_60min | W | 0.00554 | 0.0 | 0.0 | 1.898 | 1 | True |
| m6 | first_60min | switch | 0.04374 | 0.0 | 0.0 | 0.0 | 0 | True |
| m6 | first_6h | W | 0.00554 | 0.0 | 0.0 | 0.0 | 0 | True |
| m6 | first_6h | switch | 0.04374 | 0.0 | 0.0 | 0.0 | 0 | True |

## 5. m4 (two layers: not one curve): warning sign or not (`research/WAA_estm/outputs/tables/step3b3_menu_comparison.csv`)

Warning = the 60 min and 6 h curves of a replicate do not overlap (their 1.96 sd intervals) on more than
20% of the range, OR BIC(free curve) - BIC(F4) > 10 in either window.

| rep | menu_winner_60min | menu_dbic_to_second | curve_minus_best_menu_bic_60min | curve_minus_f4_bic_60min | both_intervals | windows_disagree_share | median_abs_log_ratio_6h_60min | warn_windows | warn_f4 | warning | curve_minus_f4_bic_6h |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | F4 | 51.1 | 32.1 | 32.1 | False | 0.0 | 0.094 | False | True | True | 8.2 |
| 1 | F6 | 3.8 | -2.4 | -6.2 | False | 0.0 | 0.38 | False | False | False | -31.3 |
| 2 | F4 | 65.7 | 19.9 | 19.9 | False | 0.0 | 0.027 | False | True | True | 68.9 |
| 3 | F4 | 69.5 | 31.2 | 31.2 | False | 0.0 | 0.108 | False | True | True | 108.0 |
| 4 | F4 | 57.1 | 36.1 | 36.1 | False | 0.0 | 0.056 | False | True | True | 60.8 |

The same signs for every truth (false alarms on the single-curve truths; F4 in the 6 h window was fitted for m4 only):

| truth | replicates | menu_winners_60min | curve_minus_best_menu_bic | curve_minus_f4_bic_60min | windows_disagree | warn_f4 | warning |
|---|---|---|---|---|---|---|---|
| m1 | 5 | F1 x4, F6 x1 | -20.7 | -31.6 | 0 | 0 | 0 |
| m3 | 5 | F3 x3, F5 x2 | -31.5 | -151.4 | 0 | 0 | 0 |
| m4 | 5 | F4 x4, F6 x1 | 31.2 | 31.2 | 0 | 4 | 4 |
| m5 | 5 | F5 x4, F3 x1 | -61.8 | -141.1 | 0 | 0 | 0 |
| m6 | 5 | F6 x5 | 17.2 | -23.9 | 0 | 0 | 0 |

## 6. Robustness: no spread, T1 noise

| kind | truth | window | band_coverage | curve_ok | median_k_ratio | median_slope | with_interval | flat_reps | rising_reps | step_reps | median_step_ratio | shape_wanted | shape_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| nospread | m1 | first_60min | 0.8 | True | 1.049 | -0.027 | 1 | 1 | 0 |  |  | flat | False |
| T1 | m1 | first_60min | 0.325 | False | 1.079 | -0.315 | 0 | 0 | 0 |  |  | flat | False |
| T1 | m3 | first_60min | 0.8 | True | 0.97 | 0.49 | 0 | 0 | 0 |  |  | rising | False |

| kind | truth | window | threshold | returned | true_selected | spurious | coef_ratios | equation_ok |
|---|---|---|---|---|---|---|---|---|
| nospread | m1 | first_60min | 0.02 | W + W^2 x2, W x2, W + sw0.05 x1 | W 5/5 | W^2 2/5, sw0.05 1/5 | W x1.09 | False |
| T1 | m1 | first_60min | 0.02 | W + W^2 + sw0.05 + sw0.3 x2, W + sw0.05 + sw0.3 x2, W + W^2 + sw0.05 + sw0.1 + sw0.3 x1 | W 5/5 | W^2 3/5, sw0.05 5/5, sw0.1 1/5, sw0.3 5/5 | W x3.08 | False |
| T1 | m3 | first_60min | 0.45 | W^2 x4, W x1 | W^2 4/5 | W 1/5 | W^2 x0.96 | True |

Fitted spread s (median): T1 m1 first_60min: 1.10, T1 m3 first_60min: 1.17, main m1 first_60min: 1.15, main m1 first_6h: 1.31, main m3 first_60min: 1.26, main m3 first_6h: 1.21, main m4 first_60min: 1.26, main m4 first_6h: 1.31, main m5 first_60min: 1.25, main m5 first_6h: 1.23, main m6 first_60min: 1.25, main m6 first_6h: 1.22, nospread m1 first_60min: 0.04.

## 7. Comparison with the menu fit (first 60 min; best of 4 starts, steps 3b-1c / 3b-2)

| truth | menu_winners | menu_dbic_to_second | curve_minus_best_menu_bic |
|---|---|---|---|
| m1 | F1 x4, F6 x1 | 5.0 | -20.7 |
| m3 | F3 x3, F5 x2 | 0.2 | -31.5 |
| m4 | F4 x4, F6 x1 | 57.1 | 31.2 |
| m5 | F5 x4, F3 x1 | 1.9 | -61.8 |
| m6 | F6 x5 | 16.9 | 17.2 |

BIC of the free curve counts 6 parameters (5 knots + s); positive = the menu law is preferred.

## 8. Verdict and what it means for Step 4

| truth | window_used | curve_ok | shape_ok | equation_ok | verdict | basis |
|---|---|---|---|---|---|---|
| m1 | first_60min | False | False | False | FAIL | first_60min: curve False, shape False, equation False; first_6h: curve False, shape False, equation False |
| m3 | first_60min | False | False | True | PARTIAL | first_60min: curve False, shape False, equation True; first_6h: curve True, shape False, equation False |
| m4 | both |  |  |  | PASS | warning sign in 4 of 5 replicates |
| m5 | first_6h | True | False | True | PARTIAL | first_60min: curve False, shape False, equation True; first_6h: curve True, shape False, equation True |
| m6 | first_60min | True | True | False | PARTIAL | first_60min: curve True, shape True, equation False; first_6h: curve False, shape True, equation False |

**Overall: PARTIAL.** PASS = all three criteria (curve, shape, equation) in the better window; FAIL = none; m4: PASS =
a warning sign in >= 4 of 5 replicates.

(to be written after the first complete run)

## 9. Unclear / possibly wrong

(to be written after the first complete run)
