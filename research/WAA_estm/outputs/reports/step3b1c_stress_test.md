# WAA_estm step 3b-1c: stress test of the pooled fit (fake data only)

Written by `research/WAA_estm/step3b1c_stress_test.py` (libraries `waa_pooled.py`, `waa_sim.py`; settings
`waa.stress3b1c`). No real tail is fitted and the frozen tail set v1 is not changed.

**Overall verdict: PARTIAL** (tests with a failing cell: T1, T3).
PASS = every cell meets the rule; PARTIAL = one or two of the three tests have a failing cell; FAIL = all three do.

## 1. What was run

- Window: the first 60 min of every tail. Truths m1, m4, m5, m6; candidates
  F1, F3, F4, F5, F6. 5 replicates per truth and test, each an event bootstrap of the 429-tail
  TRAIN v1 skeleton (seeds 3000 + replicate: new draws, not those of 3b-1b; the same draws for every test,
  so the tests differ only in what they change). Replicate sizes 135-753 tails,
  N of the BIC 105-502.
- Same pooled likelihood as 3b-1b (offset prior, ARMA(1,1) whitening with the sub-link's single AR(1), Laplace over
  the per-tail speed), same BIC (N = distinct (link, event draw) pairs).
- **T1** misspecified noise: truths and spread as in 3b-1b, noise with two slow time scales (section 2). Optimiser
  start = the 3a-fix calibration, which here is the truth (so T1 changes only the noise).
- **T2** honest starts: data made exactly as in 3b-1b. 4
  starts per candidate, the same for every truth; the best likelihood counts. The calibration (= truth) start is
  fitted too, only as a reference.
- **T3** faster truths: recalibrated truths (section 3), 3a-fix noise. Optimiser start = the OLD 3a-fix calibration,
  which is no longer the truth. **T3s**: the same with s = 1.5 for m1, m4.
- Nelder-Mead budget: 400 evaluations from the calibration start (as 3b-1b), 800 from
  the T2 starts. Fits stopped at the budget (test, start: stopped / fits): T1 calib: 1/100, T2 calib: 5/100, T2 f1sol: 0/80, T2 g1: 0/100, T2 g2: 0/100, T2 g3: 0/100, T2 g4: 0/20, T3 calib: 1/100, T3s calib: 0/50.
- **Cut for time**: nothing. One complete replicate (T2, m4: 5
  candidates x all starts) took 154 s; projected total 0.42 h on
  3 processes (budget 3.0 h).
- **Total run time**: 38 min wall time for 750 fits (1.9
  process-hours of fitting; cached in `.cache/fits3b1c/`). This run: 21 s.

## 2. T1 noise check (`research/WAA_estm/outputs/tables/step3b1c_noise_check.csv`)

T1 noise per sub-link = 10 s white (the sub-link's 3a-fix sigma_white) + AR(1) with tau 20 min +
AR(1) with tau 300 min + 0.3 dB quantization. The two slow parts together have the sub-link's 3a-fix
slow variance (median sigma_slow 0.150 dB, so the total slow variance is about 0.15^2 dB^2); the share of
the tau 20 min part is 0.3 for every sub-link (median sigmas
0.082 and 0.125 dB), the grid value that fits the real
statistics best. Dry series of 680 min (30 per sub-link, 249
sub-links) against the real control tails; medians over sub-links and the median absolute error per sub-link:

| statistic | real_median | sim_3afix_median | sim_T1_median | abs_err_3afix | abs_err_T1 |
|---|---|---|---|---|---|
| within_std_db | 0.1414 | 0.134 | 0.1273 | 0.0111 | 0.0182 |
| on_level_share | 0.2495 | 0.266 | 0.2575 | 0.0367 | 0.037 |
| acf_1 | 0.8012 | 0.8026 | 0.7774 | 0.032 | 0.044 |
| acf_10 | 0.572 | 0.629 | 0.5887 | 0.052 | 0.0568 |
| acf_60 | 0.2189 | 0.2298 | 0.2186 | 0.0447 | 0.0655 |

Score per share of the fast part (sum of squared scaled errors; lower is better):

| fast_share | abs_err_within_std_db | abs_err_acf_1 | abs_err_acf_10 | abs_err_acf_60 | score |
|---|---|---|---|---|---|
| 0.1 | 0.0243 | 0.053 | 0.0678 | 0.0957 | 2.3118 |
| 0.2 | 0.0213 | 0.0494 | 0.0545 | 0.0699 | 1.5334 |
| 0.3 | 0.0182 | 0.044 | 0.0568 | 0.0655 | 1.313 |
| 0.4 | 0.0166 | 0.0387 | 0.065 | 0.08 | 1.5202 |
| 0.5 | 0.0149 | 0.0364 | 0.0733 | 0.1047 | 2.0121 |
| 0.6 | 0.0132 | 0.034 | 0.0814 | 0.1377 | 2.8678 |
| 0.7 | 0.013 | 0.0334 | 0.091 | 0.1685 | 3.9668 |
| 0.8 | 0.0136 | 0.0329 | 0.101 | 0.1989 | 5.2922 |
| 0.9 | 0.0129 | 0.0346 | 0.113 | 0.2259 | 6.6842 |

## 3. T3 recalibration (`research/WAA_estm/outputs/tables/step3b1c_recalibration.csv`)

Each truth is fitted (least squares on the three relative errors) so that its noise-free median W / A_start, with
its spread, matches the real median A / A_start at 10, 30, 60 min (mean over [k, k + 5) min; the
calibration averages over 16 spread draws on the full skeleton). The slow parameters
(m4 tau_s, m6 k_evap) keep their 3a-fix value: three points inside
the first hour cannot fix a slow time constant (a first attempt with every parameter free gave m4 tau_s = 7e11 min
and m6 k_evap = 4e-20 / min, i.e. a constant floor; no fit was run on those truths).

| test | truth | s | old_parameters | new_parameters |
|---|---|---|---|---|
| T3 | m1 | 1.25 | tau 29.75 | tau 25.74 |
| T3 | m4 | 1.25 | tau_f 18.73, tau_s 193, phi 0.1714 | tau_f 14.29, tau_s 193, phi 0.1736 |
| T3 | m5 | 1.25 | a 0.00299, b 0.04499 | a 0.006573, b 0.05574 |
| T3 | m6 | 1.25 | k_run 0.04374, k_evap 0.00554, A_c 0.1615 | k_run 0.05749, k_evap 0.00554, A_c 0.1691 |
| T3s | m1 | 1.5 | tau 29.75 | tau 25.62 |
| T3s | m4 | 1.5 | tau_f 18.73, tau_s 193, phi 0.1714 | tau_f 14.24, tau_s 193, phi 0.173 |

Real vs simulated A / A_start (simulated = median over the 5 replicates of the OBSERVED fake tails, with noise and
event bootstrap; old = the 3b-1b truths):

| test | truth | real_10min | old_sim_10min | new_sim_10min | real_30min | old_sim_30min | new_sim_30min | real_60min | old_sim_60min | new_sim_60min |
|---|---|---|---|---|---|---|---|---|---|---|
| T3 | m1 | 0.568 | 0.613 | 0.571 | 0.229 | 0.324 | 0.284 | 0.142 | 0.165 | 0.143 |
| T3 | m4 | 0.568 | 0.557 | 0.488 | 0.229 | 0.297 | 0.257 | 0.142 | 0.184 | 0.164 |
| T3 | m5 | 0.568 | 0.526 | 0.465 | 0.229 | 0.301 | 0.251 | 0.142 | 0.196 | 0.149 |
| T3 | m6 | 0.568 | 0.557 | 0.496 | 0.229 | 0.294 | 0.257 | 0.142 | 0.183 | 0.159 |
| T3s | m1 | 0.568 | 0.613 | 0.58 | 0.229 | 0.324 | 0.287 | 0.142 | 0.165 | 0.144 |
| T3s | m4 | 0.568 | 0.557 | 0.479 | 0.229 | 0.297 | 0.258 | 0.142 | 0.184 | 0.172 |

Noise-free calibration values (full skeleton):

| test | truth | old_noise_free_10min | new_noise_free_10min | old_noise_free_30min | new_noise_free_30min | old_noise_free_60min | new_noise_free_60min |
|---|---|---|---|---|---|---|---|
| T3 | m1 | 0.673 | 0.633 | 0.348 | 0.295 | 0.13 | 0.095 |
| T3 | m4 | 0.603 | 0.527 | 0.301 | 0.24 | 0.158 | 0.139 |
| T3 | m5 | 0.56 | 0.493 | 0.316 | 0.253 | 0.185 | 0.136 |
| T3 | m6 | 0.606 | 0.535 | 0.296 | 0.238 | 0.155 | 0.139 |
| T3s | m1 | 0.673 | 0.632 | 0.348 | 0.295 | 0.13 | 0.095 |
| T3s | m4 | 0.603 | 0.526 | 0.301 | 0.24 | 0.158 | 0.139 |

## 4. Win shares and the decision rule (`research/WAA_estm/outputs/tables/step3b1c_win_shares.csv`, `research/WAA_estm/outputs/tables/step3b1c_decision.csv`)

Replicates (of 5) won by each candidate (lowest pooled BIC); `win_ok` = own candidate wins >=
80%; `params_ok` = every key parameter within +-30%.

| test | truth | F1 | F3 | F4 | F5 | F6 | own | win_ok | params_ok | verdict |
|---|---|---|---|---|---|---|---|---|---|---|
| T1 | m1 | 5 | 0 | 0 | 0 | 0 | F1 | True | True | pass |
| T1 | m4 | 0 | 0 | 5 | 0 | 0 | F4 | True | True | pass |
| T1 | m5 | 0 | 5 | 0 | 0 | 0 | F5 | False | True | FAIL |
| T1 | m6 | 0 | 0 | 0 | 0 | 5 | F6 | True | True | pass |
| T2 | m1 | 4 | 0 | 0 | 0 | 1 | F1 | True | True | pass |
| T2 | m4 | 0 | 0 | 4 | 0 | 1 | F4 | True | True | pass |
| T2 | m5 | 0 | 1 | 0 | 4 | 0 | F5 | True | True | pass |
| T2 | m6 | 0 | 0 | 0 | 0 | 5 | F6 | True | True | pass |
| T3 | m1 | 5 | 0 | 0 | 0 | 0 | F1 | True | True | pass |
| T3 | m4 | 0 | 0 | 5 | 0 | 0 | F4 | True | True | pass |
| T3 | m5 | 0 | 1 | 0 | 4 | 0 | F5 | True | False | FAIL |
| T3 | m6 | 0 | 0 | 0 | 0 | 5 | F6 | True | True | pass |
| T3s | m1 | 5 | 0 | 0 | 0 | 0 | F1 | True | True | pass |
| T3s | m4 | 0 | 0 | 4 | 0 | 1 | F4 | True | True | pass |

Margins (dBIC = BIC of the own candidate minus the best other candidate; negative = the own candidate wins):

| test | truth | own | own_wins | wins_instead | median_dbic_own_minus_best_other | median_dbic_where_lost |
|---|---|---|---|---|---|---|
| T1 | m1 | F1 | 5 |  | -5.7 |  |
| T1 | m4 | F4 | 5 |  | -29.3 |  |
| T1 | m5 | F5 | 0 | F3 x5 | 1.8 | 1.8 |
| T1 | m6 | F6 | 5 |  | -30.3 |  |
| T2 | m1 | F1 | 4 | F6 x1 | -5.0 | 1.5 |
| T2 | m4 | F4 | 4 | F6 x1 | -57.1 | 3.8 |
| T2 | m5 | F5 | 4 | F3 x1 | -1.9 | 0.2 |
| T2 | m6 | F6 | 5 |  | -16.9 |  |
| T3 | m1 | F1 | 5 |  | -4.7 |  |
| T3 | m4 | F4 | 5 |  | -69.2 |  |
| T3 | m5 | F5 | 4 | F3 x1 | -7.4 | 0.5 |
| T3 | m6 | F6 | 5 |  | -30.0 |  |
| T3s | m1 | F1 | 5 |  | -4.9 |  |
| T3s | m4 | F4 | 4 | F6 x1 | -43.8 | 8.4 |

Failing cells:

- **T1 m5**: own candidate F5 wins 0/5; instead: F3 x5 (median dBIC own - winner in those replicates +1.8).
- **T3 m5**: own candidate F5 wins 4/5; instead: F3 x1 (median dBIC own - winner in those replicates +0.5); parameters off: a x0.654.

![win share](../figures/step3b1c_win_share.png)

## 5. Parameter recovery (`research/WAA_estm/outputs/tables/step3b1c_recovery.csv`)

fitted / true of the own candidate: median, minimum and maximum over the 5 replicates. `s` is shown but is not in
the rule.

| test | truth | parameter | true | median_ratio | min_ratio | max_ratio | in_rule | within_tol |
|---|---|---|---|---|---|---|---|---|
| T1 | m1 | tau | 29.75393 | 1.0 | 0.681 | 1.082 | True | True |
| T1 | m1 | s | 1.25 | 0.976 | 0.877 | 0.994 | False | True |
| T1 | m4 | tau_f | 18.72823 | 0.995 | 0.91 | 1.12 | True | True |
| T1 | m4 | tau_s | 193.00825 | 1.137 | 1.026 | 1.61 | True | True |
| T1 | m4 | phi | 0.17142 | 0.957 | 0.85 | 0.979 | True | True |
| T1 | m4 | s | 1.25 | 0.972 | 0.867 | 1.003 | False | True |
| T1 | m5 | a | 0.00299 | 1.234 | 0.815 | 1.788 | True | True |
| T1 | m5 | b | 0.04499 | 1.021 | 0.917 | 1.251 | True | True |
| T1 | m5 | s | 1.25 | 1.031 | 0.874 | 1.08 | False | True |
| T1 | m6 | k_run | 0.04374 | 1.013 | 0.928 | 1.027 | True | True |
| T1 | m6 | k_evap | 0.00554 | 1.107 | 0.6 | 1.173 | True | True |
| T1 | m6 | A_c | 0.16148 | 0.972 | 0.764 | 1.023 | True | True |
| T1 | m6 | s | 1.25 | 0.982 | 0.83 | 1.064 | False | True |
| T2 | m1 | tau | 29.75393 | 0.972 | 0.712 | 1.107 | True | True |
| T2 | m1 | s | 1.25 | 0.948 | 0.868 | 0.99 | False | True |
| T2 | m4 | tau_f | 18.72823 | 0.94 | 0.792 | 1.113 | True | True |
| T2 | m4 | tau_s | 193.00825 | 0.794 | 0.504 | 1.286 | True | True |
| T2 | m4 | phi | 0.17142 | 1.041 | 0.983 | 1.34 | True | True |
| T2 | m4 | s | 1.25 | 0.966 | 0.865 | 1.012 | False | True |
| T2 | m5 | a | 0.00299 | 0.758 | 0.52 | 1.249 | True | True |
| T2 | m5 | b | 0.04499 | 1.052 | 0.883 | 1.261 | True | True |
| T2 | m5 | s | 1.25 | 1.015 | 0.882 | 1.095 | False | True |
| T2 | m6 | k_run | 0.04374 | 1.024 | 0.924 | 1.293 | True | True |
| T2 | m6 | k_evap | 0.00554 | 0.862 | 0.645 | 1.478 | True | True |
| T2 | m6 | A_c | 0.16148 | 0.984 | 0.783 | 1.192 | True | True |
| T2 | m6 | s | 1.25 | 0.928 | 0.862 | 1.054 | False | True |
| T3 | m1 | tau | 25.73736 | 1.056 | 0.727 | 1.129 | True | True |
| T3 | m1 | s | 1.25 | 0.959 | 0.872 | 1.034 | False | True |
| T3 | m4 | tau_f | 14.28566 | 1.046 | 0.843 | 1.169 | True | True |
| T3 | m4 | tau_s | 193.00825 | 1.009 | 0.623 | 1.171 | True | True |
| T3 | m4 | phi | 0.1736 | 1.053 | 1.029 | 1.221 | True | True |
| T3 | m4 | s | 1.25 | 0.974 | 0.882 | 1.026 | False | True |
| T3 | m5 | a | 0.00657 | 0.654 | 0.471 | 1.095 | True | False |
| T3 | m5 | b | 0.05574 | 1.015 | 0.912 | 1.444 | True | True |
| T3 | m5 | s | 1.25 | 1.01 | 0.912 | 1.09 | False | True |
| T3 | m6 | k_run | 0.05749 | 0.905 | 0.867 | 1.051 | True | True |
| T3 | m6 | k_evap | 0.00554 | 0.906 | 0.652 | 1.253 | True | True |
| T3 | m6 | A_c | 0.16915 | 0.963 | 0.87 | 1.177 | True | True |
| T3 | m6 | s | 1.25 | 0.979 | 0.874 | 1.059 | False | True |
| T3s | m1 | tau | 25.61663 | 0.841 | 0.675 | 1.118 | True | True |
| T3s | m1 | s | 1.5 | 0.989 | 0.832 | 1.049 | False | True |
| T3s | m4 | tau_f | 14.24137 | 1.083 | 0.916 | 1.183 | True | True |
| T3s | m4 | tau_s | 193.00825 | 1.062 | 0.64 | 1.178 | True | True |
| T3s | m4 | phi | 0.17303 | 1.054 | 1.032 | 1.237 | True | True |
| T3s | m4 | s | 1.5 | 0.946 | 0.805 | 0.984 | False | True |

## 6. Start sensitivity (T2; `research/WAA_estm/outputs/tables/step3b1c_start_sensitivity.csv`)

Starts:

- F1 (tau): g1 = [10], s 1.0; g2 = [30], s 1.0; g3 = [100], s 1.0; g4 = [300], s 0.5
- F3 (k, n): g1 = [0.01, 1.5], s 1.0; g2 = [0.1, 1.5], s 1.0; g3 = [0.03, 2.5], s 0.5
- F4 (tau_f, tau_s, phi): g1 = [10, 100, 0.3], s 1.0; g2 = [30, 300, 0.1], s 1.0; g3 = [3, 30, 0.5], s 0.5
- F5 (a, b): g1 = [0.1, 0.01], s 1.0; g2 = [0.01, 0.1], s 1.0; g3 = [0.003, 0.01], s 0.5
- F6 (k_run, k_evap, A_c): g1 = [0.1, 0.01, 0.5], s 1.0; g2 = [0.02, 0.002, 0.3], s 1.0; g3 = [0.3, 0.03, 1.0], s 0.5
- `f1sol` (every candidate but F1), from the best F1 fit of the same dataset (tau, s): F3 k = 1 / tau, n = 1.05; F4 tau_f = tau, tau_s = 5 tau, phi = 0.05; F5 a = 1 / tau, b = 0.1 / tau; F6 k_run = 0.1 / tau, k_evap = 1 / tau, A_c = 1 dB; s of the F1 fit.

`starts_differ` = datasets where the best and the worst start differ by more than 2 BIC;
`generic_worse_than_calib` = datasets where even the best generic start is more than 2 BIC worse than the fit
started at the truth (the optimum was missed); `generic_better_than_calib` = the opposite (the truth start stopped
in a worse place).

Per candidate (all truths):

| candidate | datasets | starts_differ | median_bic_range | max_bic_range | generic_worse_than_calib | generic_better_than_calib |
|---|---|---|---|---|---|---|
| F1 | 20.0 | 5.0 | 0.81 | 3.4 | 0.0 | 0.0 |
| F3 | 20.0 | 4.0 | 1.06 | 4.5 | 0.0 | 0.0 |
| F4 | 20.0 | 10.0 | 1.93 | 108.8 | 0.0 | 2.0 |
| F5 | 20.0 | 11.0 | 2.94 | 12.4 | 0.0 | 5.0 |
| F6 | 20.0 | 20.0 | 58.91 | 266.5 | 3.0 | 13.0 |

Per truth and candidate:

| truth | candidate | datasets | starts_differ | median_bic_range | max_bic_range | generic_worse_than_calib | generic_better_than_calib | fits_stopped_at_maxfev | fits | best_start_counts |
|---|---|---|---|---|---|---|---|---|---|---|
| m1 | F1 | 5 | 1 | 0.73 | 2.9 | 0 | 0 | 0 | 20 | g2 2, g4 2, g1 1 |
| m1 | F3 | 5 | 1 | 1.71 | 4.5 | 0 | 0 | 0 | 20 | g3 3, g2 1, g1 1 |
| m1 | F4 | 5 | 1 | 0.93 | 5.4 | 0 | 0 | 0 | 20 | g2 3, f1sol 1, g3 1 |
| m1 | F5 | 5 | 3 | 3.72 | 12.4 | 0 | 1 | 0 | 20 | g2 3, g1 2 |
| m1 | F6 | 5 | 5 | 8.99 | 13.1 | 0 | 4 | 0 | 20 | g1 2, g3 2, f1sol 1 |
| m4 | F1 | 5 | 0 | 0.29 | 1.8 | 0 | 0 | 0 | 20 | g4 3, g3 1, g2 1 |
| m4 | F3 | 5 | 1 | 0.94 | 3.1 | 0 | 0 | 0 | 20 | f1sol 3, g1 2 |
| m4 | F4 | 5 | 2 | 1.78 | 10.3 | 0 | 0 | 0 | 20 | g1 3, g3 1, g2 1 |
| m4 | F5 | 5 | 1 | 1.5 | 3.8 | 0 | 1 | 0 | 20 | f1sol 2, g3 2, g2 1 |
| m4 | F6 | 5 | 5 | 198.18 | 266.5 | 1 | 2 | 0 | 20 | g1 2, g2 2, f1sol 1 |
| m5 | F1 | 5 | 2 | 1.7 | 2.3 | 0 | 0 | 0 | 20 | g2 3, g3 1, g1 1 |
| m5 | F3 | 5 | 0 | 0.37 | 1.7 | 0 | 0 | 0 | 20 | g2 2, g1 1, f1sol 1, g3 1 |
| m5 | F4 | 5 | 3 | 2.09 | 8.7 | 0 | 1 | 0 | 20 | g3 2, f1sol 1, g1 1, g2 1 |
| m5 | F5 | 5 | 5 | 4.49 | 7.7 | 0 | 0 | 0 | 20 | g1 2, g2 1, g3 1, f1sol 1 |
| m5 | F6 | 5 | 5 | 38.71 | 85.7 | 0 | 5 | 0 | 20 | f1sol 2, g2 1, g1 1, g3 1 |
| m6 | F1 | 5 | 2 | 1.15 | 3.4 | 0 | 0 | 0 | 20 | g3 2, g1 1, g2 1, g4 1 |
| m6 | F3 | 5 | 2 | 1.17 | 3.3 | 0 | 0 | 0 | 20 | f1sol 2, g1 2, g3 1 |
| m6 | F4 | 5 | 4 | 4.1 | 108.8 | 0 | 1 | 0 | 20 | f1sol 2, g1 2, g2 1 |
| m6 | F5 | 5 | 2 | 1.33 | 4.5 | 0 | 3 | 0 | 20 | g2 2, f1sol 2, g1 1 |
| m6 | F6 | 5 | 5 | 178.66 | 232.2 | 2 | 2 | 0 | 20 | g2 3, g1 2 |

Winner with the generic starts vs winner with the truth start: the same in 19 of 20
datasets: different in m4 rep 1 (F4 -> F6).

## 7. Overall verdict and what it means for Step 4

**PARTIAL.** Two of the 14 cells fail, both for truth m5 (linear + quadratic); T2 passes in every cell.

- **What held.** Exponential (m1), two layers (m4) and the state switch (m6) keep winning as their own candidate in
  every test (4-5 of 5 replicates) and their key parameters stay within +-30 % (median), also m4's tau_s and m6's
  A_c: with noise the fit does not know (T1), with starts that are not the truth (T2), with the faster truths (T3)
  and with the wider spread (T3s). The margins of m4 and m6 are large (median dBIC -17 to -69); m1 wins by only
  about 5 BIC in every test.
- **What broke: F5 against F3.** With the unknown noise (T1) the m5 data are won by F3 (dW/dt = -k W^n) in 5 of 5
  replicates, by a median of 1.8 BIC. Both candidates have two parameters, so this is a near tie in likelihood, and
  the m5 truth is almost a pure quadratic (its linear part a = 0.003 / min is a 330 min time constant, invisible in
  60 min), which F3 with n near 2 copies. With the right noise the same data give F5 4 of 5, by 1.9 BIC (T2). In T3
  F5 wins 4 of 5 but its linear rate a is biased (x0.65, range 0.47-1.10). So the fit can say "the decay is
  nonlinear in W" but NOT whether it is a W^n law or linear + quadratic, and it cannot measure the linear rate of F5
  from the first hour.
- **Starts matter, mostly for F6.** The best and worst of the four starts differ by more than 2 BIC in 50 of 100
  (dataset, candidate) pairs: F1 5 / 20, F3 4 / 20, F4 10 / 20, F5 11 / 20, F6 20 / 20 (F6: median spread 59 BIC,
  up to 267). For the true m6, even the best of four starts was more than 2 BIC worse than the truth start in 2 of 5
  datasets (it still won). In one m4 dataset F6 wins from a generic start (by 3.8 BIC) although F4 won when every
  candidate started at its calibration: the single start of 3b-1b had left F6 in a poor optimum there, which flattered F4.

**For Step 4 (real tails).** (1) Go ahead with the pooled fit, but run every candidate from several starts (at least
the four used here; more for F6, e.g. a small grid of A_c) and report the spread over the starts. (2) Read F3 and F5
as ONE family ("nonlinear in W"): a win of one over the other by a few BIC means nothing. (3) A win of F4 over F6
or of F6 over F4 by less than about 10 BIC is not decisive (F6 beat a true m4 by 3.8 and 8.4 BIC in single
replicates); a win of F1 is always by a small margin. (4) The single-AR(1) noise model is good enough for the
exponential / two-layer / switch question; it is not for F3 against F5. (5) Parameters to trust from the first hour:
tau, tau_f, phi, tau_s (the range over replicates reaches 0.5-1.6), k_run, A_c, b; not F5's a.
All of this is on fake tails whose truth is one of the candidates; real tails may follow none of them.

## 8. Unclear / possibly wrong

1. **T3 recalibration is a choice of mine.** With every parameter free, matching the three real points drove m4's
   tau_s and m6's k_evap to a constant floor (tau_s = 7e11 min), which no 60 min window could recover. I stopped
   that run before any T3 fit and kept the slow parameters at their 3a-fix value (`waa.stress3b1c.t3.keep`). The
   learners may prefer another rule (e.g. an upper bound on tau_s). This was a change after seeing the calibration,
   not after seeing fit results.
2. **The recalibrated truths do not match the real points exactly.** Noise-free they give 0.49-0.63 / 0.24-0.30 /
   0.10-0.14 against the real 0.568 / 0.229 / 0.142: m1 has one parameter and cannot follow the shape (too slow at
   10 and 30 min, too fast at 60), the others overshoot the 10 min drop. The observed fake medians (5 bootstrap
   replicates) are 0.25-0.29 at 30 min: closer to the real 0.229 than the old truths (0.29-0.32), still above it.
   The real early decay is therefore still a little faster than every T3 truth.
3. **Noise-free and observed medians differ** (m1 at 60 min: 0.095 noise-free, 0.143 observed). I believe this is
   the event bootstrap (3 events hold 246 of 429 tails, and links share half the spread) plus the offset, but I
   did not separate the causes.
4. **T1 noise matches the real dry statistics a little worse than the 3a-fix noise**: median absolute error per
   sub-link 0.018 vs 0.011 dB (within-tail std), 0.044 vs 0.032 (ACF 1 min), 0.057 vs 0.052 (10 min), 0.066 vs 0.045
   (60 min); the network medians match as well or better. Reason: one fast share (0.3) for the whole network
   instead of a fit per sub-link, and the total slow variance was kept, not refitted, so the within-tail std is
   about 10 % low (0.127 vs 0.141 dB). I read "total slow variance about 0.15^2" as the sub-link's own 3a-fix value
   (median 0.15 dB), not 0.15 dB for every sub-link.
5. **T1 is one kind of misspecification** (two AR(1) instead of one, same total variance, Gaussian). Real noise may
   also have steps, a daily cycle or TSL jumps; none is tested.
6. **Only 5 replicates**: 4 of 5 and 5 of 5 are not distinguishable statistically, and a median of 5 ratios is
   rough (several ranges reach 0.5-1.6). A cell at exactly 4 of 5 (T2 m1, m4, m5; T3s m4) passes the rule but is
   one replicate from failing.
7. **Replicate sizes vary a lot** (135-753 tails, event bootstrap), more than in 3b-1b (289-832), so power differs
   between replicates.
8. **The reference start in T2 is the truth only for the own candidate.** For the other candidates "calib" is just
   one more start, so `generic_better_than_calib` is expected there; the alarming column is
   `generic_worse_than_calib` (3 of 100, all F6).
9. **Budget**: the generic starts got 800 evaluations, the calibration start 400 (as in 3b-1b); 5 of 100 T2
   calibration-start fits stopped at 400, no generic one stopped at 800. Part of "generic better than calib" may
   be this budget and not the start.
10. **The key figure limit is reached** (6 key outputs in this module), so the win-share figure is saved as diag.
    To make it key, another key output has to be demoted: the learners' choice.
11. **Seeds**: new event draws (3000 + replicate), not the 3b-1b draws, so T2 is not a rerun of the same datasets;
    "the winner with the truth start" is fitted here on the new draws.

