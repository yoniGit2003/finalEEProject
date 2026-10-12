# WAA_estm step 4: the drying law on the REAL tails

Written by `research/WAA_estm/step4_real_drying.py` (libraries `waa_pooled.py`, `waa_curve.py`, `waa_sindy.py`; settings
`waa.real4`, the fit settings are read unchanged from `waa.stress3b1c`, `waa.curve3b3`, `waa.curve3b3b`, `waa.curve3b3c`).
Frozen tail set v1, TRAIN tails only; the 123 VAL tails were not read.

**Result.** Between 0.35 and 2.83 dB, in the first hour after rain, the drying rate per dB of water falls with slope -0.87 +- 0.1 on the main set (bootstrap refits -1.43..-0.82; -0.05 without the rise-heavy events), a shape no fake truth had; the data prefer family "exponential" by BIC (margin 11.4 to F4), but that ranking against F3 / F5 / F6 rests on 3 rising tails fitted with a negative start, and like for like the winner is F6 (state switch) by 21.5, not decisive with residuals 4.1 times the dry-noise variance; the formula found is not a law (E-SINDy keeps 5 of 5 library terms: h(W) = 0.1391 W + 0.1454 W^2 - 0.03962 W^3 - 0.5928 max(W - 0.5, 0) + 0.3592 max(W - 0.75, 0)); below 0.35 dB the method cannot tell.

## 1. What was run, cut, run time

- **Recipe** (only tools that passed on fake data): pooled likelihood of 3b-3b (offset prior, ARMA(1,1) noise, speed spread
  marginalised), first 60 min of each tail. (a) Menu fits F1, F3, F4, F5, F6: 4 starts each (3 generic + the
  F1 solution; F1: 4 generic; 800 evaluations), BIC with N = distinct (link, event) pairs. (b) Free
  curve: 5 knots from 0.35 dB to W_hi, 4 starts, penalty by 2-fold cross-validation over
  events, read only for W >= 0.35 dB, no step tool; band = 20 event-bootstrap refits, each from
  2 starts (fit, flat). (c) Jump 2 on that curve. (d) Literature.
- **Fits**: 160 menu fits (8 dataset-windows), 7 free-curve fits, 20 bootstrap refits
  (best start: fit x15, flat x5).
- **Timing**: one likelihood evaluation on the main set 121 ms (menu law) / 97 ms
  (free curve); one F1 fit 1.3 s. Projection 1.19 h on 3 processes
  (budget 3.0 h). **Cut: nothing.**
- **Run time**: 36 min wall (1.8 process-hours; cached in `.cache/fits4/`). This run:
  24 s.
- **Software**: numpy / scipy / numba (own code); PySR used for 4 curves.
- **Added after the first run, from the cached fits (no fit was rerun, no setting changed)**: the like-for-like menu
  reading (row `F1+`, section 3), the model-free rate table and the share of tails that do not decay (section 4), the
  residual rows of the set without the rise-heavy events (check 4).

## 2. Real W range, tails, events, links

Main set: **429 tails, 24 events, 159 links** (249 sub-links);
429 tails have >= 10 valid minutes in the window; N of the BIC =
300 (link, event) pairs. **W_hi = 2.83 dB** (the level with
200 tail-minutes of A above it in the first 60 min); 14402 of 25713
tail-minutes lie inside [0.35, W_hi]. Knots: 0.35, 0.60, 0.91, 1.30, 2.83 dB. The 3 big events:

| event | tails |
|---|---|
| E20150804_1830 | 115 |
| E20150826_0520 | 75 |
| E20150628_0530 | 56 |

## 3. Menu fits, main set, first 60 min (`research/WAA_estm/outputs/tables/step4_menu_fits.csv`)

| candidate | family | k | ll | bic | dbic | dbic_like_for_like | parameters | s | best_start | bic_spread_starts |
|---|---|---|---|---|---|---|---|---|---|---|
| F1 | exponential | 2 | -54183.86 | 108379.14 | 0.0 |  | tau 25.23 | 1.398 | g2 | 27.9 |
| F4 | two layers | 4 | -54183.87 | 108390.55 | 11.41 |  | tau_f 25.2, tau_s 3.953e+08, phi 1.309e-08 | 1.398 | g2 | 4.5 |
| F6 | state switch | 4 | -55483.09 | 110989.0 | 2609.86 | 0.0 | k_run 0.04274, k_evap 0.003278, A_c 0.05456 | 1.368 | g2 | 192.2 |
| F5 | nonlinear in W (F3 / F5) | 3 | -55499.56 | 111016.22 | 2637.09 | 27.22 | a 0.03961, b 8.652e-16 | 1.398 | g2 | 52.4 |
| F3 | nonlinear in W (F3 / F5) | 3 | -55499.56 | 111016.22 | 2637.09 | 27.22 | k 0.03961, n 1 | 1.398 | g1 | 9.5 |
| F1+ | exponential | 2 | -55499.56 | 111010.53 |  | 21.53 | tau 25.23 (start W0 >= 0) | 1.398 | F1 solution |  |

**Reading (3b-1c rules: F3 / F5 one family; F4 vs F6 by < 10 BIC not decisive): winner **F1** (exponential); runner-up F4 at +11.4 BIC; best candidate outside the winner's family: F4 at +11.4; F1 at +0.0.**

**Found after the run: this ranking is not like for like.** F1 and F4 are linear in the start W0, and `waa_pooled.py`
profiles it without a bound; F3, F5, F6 and the free curve keep W0 >= 0.001 dB. 3 of the
429 tails (rising tails) get a NEGATIVE start under F1, and they alone are worth the whole gap
between F1 and F3 at n = 1, which is the same law (log-likelihood -54184 vs
-55500). Row `F1+` = the exponential with W0 >= 0; column `dbic_like_for_like` compares
F1+, F3, F5, F6 only. **Like for like: winner F6 (state switch); next family F1+ at
+21.5 BIC; F1+ at +21.5.** F3 ends at its lower bound n = 1 and F5 at
b = 0: the data do not want a rate that RISES with W. On fake data no tail rises, so this did not show.

- `dbic` = BIC minus the best BIC. `bic_spread_starts` = best minus worst start of that candidate.
- F6's A_c is a fitted number, **not a physical level** (3b-3c: the likelihood prefers a feature near 0.2-0.4 dB
  whatever the truth). F5's linear rate a is not interpretable from the first hour (3b-1c).
- Fitted speed spread s: F1 1.398, F4 1.398, F6 1.368, F5 1.398, F3 1.398, F1+ 1.398
  (the fake data were made with s = 1.25).

### (d) Literature

F1's fitted time constant is **tau = 25.2 min** in the first hour (the median tail; the speed spread s = 1.40
moves single tails around it). Over the sets it is 13.9 min (alternative set), 15.3 min (without the rise-heavy events)
and 25-28 min (the others); on the first 6 h of the main set 18.7 min.

**No literature value can be compared: none is quoted in the repository, and none is invented here.** Checked:
`docs/proposal.md` names the dynamic model of Schleiss et al. [1] and the corrections of Pastorek et al. [9] without any
number; `pipeline/docs/openmrg_summary.md` only says that the RSL recovery lags the rain in the Torslanda case;
`pipeline/docs/reference/cml_libraries_howto.md` lists pycomlink's `waa_schleiss_2013(..., tau)` without a value;
`roadmap.md` names "Minda & Nakamura exponential drying", which is not in the proposal's reference list and has no
value either. The only paper stored in the repository is the OpenMRG paper [13]. To do before Step 8a: read the time
constants from [1] and [9] themselves and add them, with page numbers, to `roadmap.md`.

## 4. Free curve (`research/WAA_estm/outputs/tables/step4_curves.csv`)

![real k(W)](../figures/step4_k_real.png)

*Band (hatched): 20 event-bootstrap refits, 5-95 %. On fake data this band was about 3 times too narrow.*

- **Slope of log k vs log W over [0.35, 2.83] dB: -0.87** (quoted +-0.1; bootstrap refits
  -1.43..-0.82, min-max -1.63..-0.71).
  Reading by the rule written before the run: **falling (no fake truth looked like this)**.
- Level: k = 0.1550 / min at 0.35 dB, 0.0300 at 1 dB,
  0.0140 at 2.83 dB (1 / k = 6, 33,
  71 min); quoted +-20 %. Bootstrap band half-width (median over the range): +-94%.
- Information only: BIC(free curve) - BIC(menu winner) = -9923.3 (the free curve had a lower BIC
  than the TRUE law on fake data too, 3b-3; not a test).

| W_db | k_per_min | band_lo | band_hi |
|---|---|---|---|
| 0.35 | 0.155 | 0.0435 | 0.1911 |
| 0.4576 | 0.1986 | 0.0521 | 0.2203 |
| 0.5982 | 0.1311 | 0.0662 | 0.1478 |
| 0.782 | 0.0382 | 0.0378 | 0.0877 |
| 1.0223 | 0.0315 | 0.0312 | 0.0582 |
| 1.3364 | 0.0766 | 0.0352 | 0.0873 |
| 1.7471 | 0.0779 | 0.0178 | 0.0798 |
| 2.284 | 0.0351 | 0.0078 | 0.0358 |

| dataset | knots_db | k_at_knots_per_min | penalty | s | bic | best_start | start_cost_spread |
|---|---|---|---|---|---|---|---|
| main | 0.35 0.6 0.907 1.3 2.83 | 0.1550 0.1297 0.0267 0.0717 0.0140 | 0.01 | 0.882 | 98455.8 | step | 1700.8 |
| no_rise_heavy | 0.35 0.53 0.79 1.19 2.4 | 0.0476 0.0751 0.0439 0.0921 0.0278 | 1.0 | 0.926 | 49533.2 | step | 328.9 |
| no_shift_undershoot | 0.35 0.58 0.9 1.3 2.75 | 0.1085 0.1381 0.0301 0.0761 0.0136 | 0.1 | 0.938 | 85644.7 | step | 1222.1 |
| alt | 0.35 0.51 0.68 0.95 2.83 | 0.1081 0.0748 0.0187 0.0441 0.0080 | 0.1 | 1.077 | 208147.9 | falling | 1197.5 |
| loo_E20150804_1830 | 0.35 0.61 0.943 1.31 2.69 | 0.0743 0.0609 0.0133 0.0444 0.0119 | 0.1 | 1.244 | 76084.5 | falling | 1195.5 |
| loo_E20150826_0520 | 0.35 0.61 0.937 1.3 2.77 | 0.1034 0.1525 0.0880 0.0601 0.0034 | 1.0 | 1.327 | 83278.0 | step | 1116.6 |
| loo_E20150628_0530 | 0.35 0.59 0.9 1.3 2.81 | 0.0623 0.0900 0.0540 0.0883 0.0148 | 1.0 | 1.086 | 81826.2 | step | 809.0 |

`start_cost_spread` = worst minus best of the 4 starts after 150 evaluations (fake data, 3b-3b: median 8).

**Model-free look (information, not a tool of the recipe).** Main set, first 60 min: A smoothed over
7 min per tail, dA/dt by central difference, k = -mean(dA/dt) / median(A) per level bin.
No offset, no speed factor, no noise model: only a check that the falling trend is in the data and not made by the fit.

| A_bin_db | tail_minutes | tails | median_A_db | mean_dadt_db_per_min | median_dadt_db_per_min | k_from_mean_per_min |
|---|---|---|---|---|---|---|
| (0.35, 0.5] | 2212 | 342 | 0.4167 | -0.0218 | -0.0143 | 0.0522 |
| (0.5, 0.7] | 2514 | 371 | 0.5937 | -0.0283 | -0.019 | 0.0477 |
| (0.7, 1.0] | 3142 | 392 | 0.8457 | -0.036 | -0.025 | 0.0425 |
| (1.0, 1.4] | 3078 | 380 | 1.1814 | -0.0333 | -0.0214 | 0.0282 |
| (1.4, 2.0] | 2006 | 183 | 1.5767 | -0.0232 | -0.0036 | 0.0147 |
| (2.0, 2.83] | 455 | 47 | 2.3686 | -0.0389 | -0.0071 | 0.0164 |
| (2.83, 6.0] | 143 | 21 | 3.0171 | -0.0496 | -0.0036 | 0.0164 |

Tails that do not decay (A / A_start, mean of A over 5 min):

| tails | n | share_A_ge_90pct_of_start_at_10min | share_A_ge_100pct_of_start_at_10min | share_A_ge_90pct_of_start_at_30min | share_A_ge_100pct_of_start_at_30min |
|---|---|---|---|---|---|
| all | 429 | 0.214 | 0.133 | 0.1 | 0.061 |
| events flagged event_rise_heavy | 164 | 0.378 | 0.256 | 0.201 | 0.152 |
| other events | 265 | 0.113 | 0.057 | 0.038 | 0.004 |

## 5. Jump 2 on the curve above 0.35 dB (`research/WAA_estm/outputs/tables/step4_jump2.csv`)

Library W, W^2, W^3, sw0.5, sw0.75 (sw x = max(W - x dB, 0)); h(W) = W k(W) on 40 log-spaced points.
`relative_rms_misfit` = rms of formula / curve - 1. Main set:

| tool | formula | family | relative_rms_misfit | note |
|---|---|---|---|---|
| E-SINDy (waa_sindy elimination, bagged) | 0.1391 W + 0.1454 W^2 - 0.03962 W^3 - 0.5928 max(W - 0.5, 0) + 0.3592 max(W - 0.75, 0) | multi-term (not a named law) | 0.2259 | inclusion: W 0.85, W^2 0.85, W^3 0.96, sw0.5 1.00, sw0.75 1.00; coefficients = least squares on the kept terms |
| best-subset, 1 term | 0.03231 W | exponential | 0.6053 | BIC +0 vs the 1-term formula |
| best-subset, 2 terms | 0.07869 W - 0.07107 max(W - 0.5, 0) | multi-term (not a named law) | 0.5015 | BIC -11 vs the 1-term formula |
| best-subset, 3 terms | 0.3646 W^2 - 0.07489 W^3 - 0.5128 max(W - 0.5, 0) | multi-term (not a named law) | 0.3442 | BIC -38 vs the 1-term formula; BIC pick |
| PySR, simplest 1 of its Pareto front | 0.0436851 |  | 0.4937 | complexity 1 |
| PySR, simplest 2 of its Pareto front | (1.1669207 ^ W) * 0.036843102 |  | 0.4856 | complexity 5 |
| PySR, simplest 3 of its Pareto front | (exp(W) * 0.013390888) / W |  | 0.4798 | complexity 6 |
| PySR, its own 'best' pick | (W - (-1)*1.0516075/W)*0.019157318 |  | 0.4499 | Pareto front of 8 formulas |

- **E-SINDy (main tool): h(W) = 0.1391 W + 0.1454 W^2 - 0.03962 W^3 - 0.5928 max(W - 0.5, 0) + 0.3592 max(W - 0.75, 0)** (multi-term (not a named law)).
- The same elimination on each of the 20 bootstrap curves: "W + W^2 + W^3 + sw0.5 + sw0.75" x15, "W + W^2 + sw0.5 + sw0.75" x1, "W^2 + W^3 + sw0.5 + sw0.75" x1, "W + W^2 + W^3" x1, "W + W^3 + sw0.5 + sw0.75" x1, "W + W^2 + W^3 + sw0.75" x1.
- Best-subset is read as "true terms + at most one extra" (3b-3c): the curve is smooth and has no noise of its own,
  so BIC always prefers more terms; the misfit column says how much each extra term buys.

## 6. Agreement between menu, curve and jump 2

| tool | says | class |
|---|---|---|
| menu fits (BIC), as run | F1: exponential | exponential |
| menu fits, like for like (W0 >= 0) | F6: state switch; F3 at n = 1, F5 at b = 0 | other |
| free curve (slope -0.87) | falling (no fake truth looked like this) | other |
| jump 2 (E-SINDy) | 0.1391 W + 0.1454 W^2 - 0.03962 W^3 - 0.5928 max(W - 0.5, 0) + 0.3592 max(W - 0.75, 0): multi-term (not a named law) | other |

**They do NOT agree** (classes: exponential / quadratic family / other; "mildly
rising" + a multi-term formula + F4 or F6 counts as agreement on "other", which is what a switch or two layers looked
like on fake data).

## 7. Checks

### Checks 1 and 2: sensitivity sets and leaving out each big event (`research/WAA_estm/outputs/tables/step4_summary.csv`)

| dataset | what | tails | events | tails_big3 | w_hi_db | menu_winner | margin_family_bic | f1_negative_start_tails | like_for_like_winner | like_for_like_margin_family_bic | like_for_like_f1_minus_winner_bic | f3_n_at_lower_bound | f1_tau_min | curve_slope | curve_slope_common_range | curve_reading | jump2_family | agree |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| main | v1 main | 429 | 24 | 246 | 2.83 | F1 | 11.4 | 3 | F6 | 21.5 | 21.5 | True | 25.2 | -0.872 | -0.745 | falling (no fake truth looked like this) | multi-term (not a named law) | False |
| no_rise_heavy | without event_rise_heavy | 265 | 15 | 190 | 2.4 | F6 | 41.9 | 0 | F6 | 51.8 | 51.8 | True | 15.3 | -0.054 | -0.054 | exponential in this range | multi-term (not a named law) | False |
| no_shift_undershoot | without level_shift / undershoot | 379 | 24 | 224 | 2.75 | F1 | 10.7 | 3 | F6 | 27.4 | 27.4 | True | 25.2 | -0.791 | -0.656 | falling (no fake truth looked like this) | multi-term (not a named law) | False |
| alt | alternative set (A_start >= 0.6 dB) | 1155 | 34 | 429 | 2.83 | F4 | 11.0 | 15 | F6 | 216.5 | 216.5 | True | 13.9 | -0.689 | -0.482 | falling (no fake truth looked like this) | multi-term (not a named law) | False |
| loo_E20150804_1830 | without event E20150804_1830 | 314 | 23 | 131 | 2.69 | F1 | 10.8 | 3 | F6 | 35.1 | 35.1 | True | 28.0 | -0.595 | -0.51 | falling (no fake truth looked like this) | multi-term (not a named law) | False |
| loo_E20150826_0520 | without event E20150826_0520 | 354 | 23 | 171 | 2.77 | F1 | 11.6 | 3 | F1+ | 0.9 | 0.0 | True | 27.3 | -1.6 | -1.368 | falling (no fake truth looked like this) | multi-term (not a named law) | False |
| loo_E20150628_0530 | without event E20150628_0530 | 373 | 23 | 190 | 2.81 | F1 | 2.0 | 3 | F6 | 106.2 | 106.2 | True | 25.9 | -0.464 | -0.261 | falling (no fake truth looked like this) | multi-term (not a named law) | False |

Jump 2 per sensitivity set (E-SINDy formula; PySR's own pick):

| dataset | tool | formula | relative_rms_misfit |
|---|---|---|---|
| main | E-SINDy (waa_sindy elimination, bagged) | 0.1391 W + 0.1454 W^2 - 0.03962 W^3 - 0.5928 max(W - 0.5, 0) + 0.3592 max(W - 0.75, 0) | 0.2259 |
| main | PySR, its own 'best' pick | (W - (-1)*1.0516075/W)*0.019157318 | 0.4499 |
| no_rise_heavy | E-SINDy (waa_sindy elimination, bagged) | 0.1824 W^2 - 0.06467 W^3 - 0.172 max(W - 0.5, 0) + 0.1321 max(W - 0.75, 0) | 0.125 |
| no_rise_heavy | PySR, its own 'best' pick | W*(W*(-0.015824063) - 1*(-0.039747372))*exp(W) | 0.2221 |
| no_shift_undershoot | E-SINDy (waa_sindy elimination, bagged) | 0.099 W + 0.1503 W^2 - 0.04202 W^3 - 0.4622 max(W - 0.5, 0) + 0.2645 max(W - 0.75, 0) | 0.2339 |
| no_shift_undershoot | PySR, its own 'best' pick | 0.17854437/(W*(-1*(-1.1377202) + 8.072217/exp(W**4.0856757))) | 0.3161 |
| alt | E-SINDy (waa_sindy elimination, bagged) | 0.2258 W - 0.3579 W^2 + 0.05568 W^3 + 0.0547 max(W - 0.5, 0) + 0.4142 max(W - 0.75, 0) | 0.1771 |
| alt | PySR, its own 'best' pick | 0.037836667/exp(W/exp(W**4.2890086)) | 0.4606 |

`curve_slope` is over each set's own [0.35, W_hi]; `curve_slope_common_range` over [0.35, 2.40] dB.
`tails_big3` = tails from the 3 big events.

![sets](../figures/step4_k_sets.png)

### Check 3: menu fits on the first 6 h (information only)

winner **F4** (two layers); runner-up F1 at +5.5 BIC; best candidate outside the winner's family: F1 at +5.5; F1 at +5.5.

| candidate | family | k | ll | bic | dbic | dbic_like_for_like | parameters | s | best_start | bic_spread_starts |
|---|---|---|---|---|---|---|---|---|---|---|
| F4 | two layers | 4 | -124540.67 | 249104.15 | 0.0 |  | tau_f 17.53, tau_s 1.605e+04, phi 0.01496 | 1.299 | g1 | 33.7 |
| F1 | exponential | 2 | -124549.13 | 249109.66 | 5.51 |  | tau 18.73 | 1.288 | g3 | 17.1 |
| F6 | state switch | 4 | -125760.79 | 251544.4 | 2440.25 | 0.0 | k_run 0.0523, k_evap 0.001926, A_c 0.09453 | 1.287 | g2 | 113.8 |
| F5 | nonlinear in W (F3 / F5) | 3 | -125816.51 | 251650.13 | 2545.98 | 105.73 | a 0.04206, b 5.266e-16 | 1.34 | g3 | 160.1 |
| F3 | nonlinear in W (F3 / F5) | 3 | -125816.86 | 251650.83 | 2546.67 | 106.43 | k 0.04201, n 1 | 1.342 | g3 | 3.0 |
| F1+ | exponential | 2 | -125816.86 | 251645.13 |  | 100.73 | tau 18.73 (start W0 >= 0) | 1.288 | F1 solution |  |

### Check 4: fit quality (`research/WAA_estm/outputs/tables/step4_fit_quality.csv`)

Residual = A minus the fitted curve of the tail (pooled law at the tail's speed factor, profiled start and offset).
`noise_model_*` = the same statistic of pure dry-weather noise drawn from each sub-link's noise model, on the same
minutes. `whitened_variance_ratio` = variance of the residuals after the noise filter of the likelihood, divided by
what the dry noise would give (1 = the residuals are exactly dry noise). Statistics inside a 60 min window are not
comparable with the control-tail row (tails about 12 h long): compare each row with its own `noise_model_*`.

| law | window_min | tails | within_tail_std_db | noise_model_std_db | whitened_variance_ratio | acf_1 | noise_model_acf_1 | acf_10 | noise_model_acf_10 | acf_60 | noise_model_acf_60 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| F1 (winner) | 60 | 429 | 0.1397 | 0.0879 | 4.13 | 0.802 | 0.603 | 0.185 | 0.141 |  |  |
| F1 (winner of the first hour) | 360 | 429 | 0.146 | 0.119 | 2.079 | 0.802 | 0.764 | 0.5 | 0.511 | 0.076 | 0.029 |
| F6, set without event_rise_heavy | 60 | 265 | 0.1236 | 0.0884 | 3.475 | 0.744 | 0.611 | 0.19 | 0.139 |  |  |
| F1, set without event_rise_heavy | 60 | 265 | 0.1273 | 0.0897 | 3.479 | 0.746 | 0.607 | 0.194 | 0.139 |  |  |
| real control tails (dry weather, median over the sub-links of the main set) | 677 | 5198 | 0.1414 |  |  | 0.801 |  | 0.572 |  | 0.219 |  |

![gallery](../figures/step4_gallery.png)

### Check 5: real numbers against the fake data (`research/WAA_estm/outputs/tables/step4_fake_comparison.csv`)

Slope over the read range, and median A / A_start at 10 / 30 / 60 min (mean of A over 5 min). `slope_inside` allows
+-0.1. Real ranges: bootstrap over events, 5-95 %.

| source | slope | slope_range | ratio_10min | ratio_30min | ratio_60min | range_10min | range_30min | range_60min | slope_inside | real_inside_10min | real_inside_30min | real_inside_60min |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| REAL tails (TRAIN v1) | -0.872 | -1.43..-0.82 (bootstrap 5-95 %) | 0.568 | 0.229 | 0.142 | 0.386..0.795 | 0.168..0.483 | 0.091..0.255 |  |  |  |  |
| fake truth m1 (3b-3b datasets, 5 replicates) | 0.012 | -0.01..0.05 | 0.635 | 0.331 | 0.164 | 0.579..0.636 | 0.292..0.353 | 0.152..0.190 | False | False | False | False |
| fake truth m3 (3b-3b datasets, 5 replicates) | 1.048 | 0.99..1.12 | 0.516 | 0.302 | 0.202 | 0.468..0.519 | 0.269..0.327 | 0.169..0.214 | False | False | False | False |
| fake truth m4 (3b-3b datasets, 5 replicates) | 0.418 | 0.36..0.47 | 0.573 | 0.318 | 0.192 | 0.523..0.581 | 0.281..0.330 | 0.168..0.202 | False | True | False | False |
| fake truth m5 (3b-3b datasets, 5 replicates) | 0.962 | 0.93..1.03 | 0.539 | 0.312 | 0.204 | 0.488..0.541 | 0.276..0.331 | 0.174..0.215 | False | False | False | False |
| fake truth m6 (3b-3b datasets, 5 replicates) | 0.296 | 0.20..0.32 | 0.577 | 0.316 | 0.186 | 0.516..0.581 | 0.278..0.326 | 0.158..0.196 | False | True | False | False |
| fake truth m1 recalibrated to the real early decay (3b-1c T3) |  |  | 0.571 | 0.284 | 0.143 | 0.550..0.624 | 0.228..0.325 | 0.127..0.151 |  | True | True | True |
| fake truth m4 recalibrated to the real early decay (3b-1c T3) |  |  | 0.488 | 0.257 | 0.164 | 0.457..0.535 | 0.236..0.285 | 0.156..0.185 |  | False | False | False |
| fake truth m5 recalibrated to the real early decay (3b-1c T3) |  |  | 0.465 | 0.251 | 0.149 | 0.434..0.505 | 0.223..0.276 | 0.142..0.172 |  | False | True | False |
| fake truth m6 recalibrated to the real early decay (3b-1c T3) |  |  | 0.496 | 0.257 | 0.159 | 0.466..0.545 | 0.234..0.279 | 0.155..0.174 |  | False | False | False |

## 8. The result in one sentence

**Between 0.35 and 2.83 dB, in the first hour after rain, the drying rate per dB of water falls with slope -0.87 +- 0.1 on the main set (bootstrap refits -1.43..-0.82; -0.05 without the rise-heavy events), a shape no fake truth had; the data prefer family "exponential" by BIC (margin 11.4 to F4), but that ranking against F3 / F5 / F6 rests on 3 rising tails fitted with a negative start, and like for like the winner is F6 (state switch) by 21.5, not decisive with residuals 4.1 times the dry-noise variance; the formula found is not a law (E-SINDy keeps 5 of 5 library terms: h(W) = 0.1391 W + 0.1454 W^2 - 0.03962 W^3 - 0.5928 max(W - 0.5, 0) + 0.3592 max(W - 0.75, 0)); below 0.35 dB the method cannot tell.**

## 9. In plain words: what this means for the WAA equation

**The tools were built to answer "exponential, quadratic, or something in between". The real tails gave an answer
that was not on the list, and they also exposed a fault in the menu comparison. Step 4 does NOT deliver a drying law.**

1. **The curve falls.** In k(W) = "share of the remaining water that leaves per minute", fake truths were flat
   (exponential) or rising (quadratic, switch, two layers). The real curve goes the other way: 1 / k is about 6 min
   at 0.35 dB, 33 min at 1 dB and 71 min at 2.8 dB (slope -0.87). In words: while A is high it hardly moves, and once
   it is lower it drops fast.
2. **This is in the data, not made by the fit.** The model-free table shows the same trend (k from 0.052 / min at
   0.4 dB down to 0.016 above 1.4 dB), and the gallery shows why: some tails stay flat or RISE after t0_trace before
   they drop. 21 % of the tails are still at >= 90 % of their start level after 10 min, and 13 % are at or above it.
3. **It comes mostly from the events flagged event_rise_heavy** (164 of the 429 tails; 38 % of them are still at
   >= 90 % after 10 min, against 11 % of the others). Without those events the slope is -0.05, which the pre-set rule
   reads as "exponential in this range" (F1 tau 15 min). In all the other sets it is between -0.46 and -1.60.
4. **The curve itself is badly determined.** The four starts end 800-1700 apart in cost (fake data: about 8), 5 of
   the 20 bootstrap refits found a better fit from the flat start, the band is +-94 % wide, and the dip near 0.9 dB
   moves with the third knot. Only the overall trend (not rising) is read; the bumps are not.
5. **No formula.** E-SINDy keeps all 5 library terms (that is a description of a bumpy curve, not a law). The best
   one-term formula, h = 0.032 W (exponential), misses the curve by 61 % rms; PySR's simplest formula is a constant,
   h = 0.044 dB / min (49 %).
6. **The menu cannot decide either.** As run, F1 wins, but only because F1 and F4 may give a rising tail a negative
   start (3 tails, worth 2630 BIC). Like for like, F6 leads the exponential by 21 BIC on the main set (over the sets: from a tie, F1+ ahead by 0.9,
   to F6 ahead by 216), which is not decisive: the residuals are 4.1 times larger (in variance) than the dry-weather noise
   the likelihood assumes, so every BIC margin is overstated by about that factor, and F6 depends on its start (192
   BIC between starts). **One thing is stable in all 7 sets: F3 ends at n = 1 and F5 at b = 0.** The data never ask
   for a rate that rises with W, so the quadratic family (m3 / m5) is not supported in the first hour.
7. **What it means for the WAA equation.** dW/dt = -h(W) with one curve h and a random speed per tail does not
   describe the first hour of the TRAIN tails. A tail that stays flat or rises is not "drying slowly": something is
   still adding attenuation after the last rain trace (rain or drizzle below what radar and gauges see, re-wetting,
   condensation, or attenuation along the path). That is an input, i.e. the g(R, W) side of
   dW/dt = g(R, W) - h(W), or a wrong start time, not a property of h. Until this is separated, the exponential
   h(W) = W / tau with tau of 15-28 min is a placeholder, not a discovered law, and nothing can be said about
   switch vs two layers.
8. **Decisions for the learners before Step 6** (not taken here): (a) which tails define "drying": all of v1, or
   without the rise-heavy events (a data choice, made without looking at a model fit, so not circular), or tails
   re-started where the decay begins (a new frozen set); (b) repair the menu comparison (the same bound on W0 for every
   candidate) and re-run the fake-data tests with some rising tails; (c) the noise during drying is about twice the
   dry-weather noise (std 0.14 vs 0.09 dB in a 60 min window): the likelihood needs a larger or fitted noise level
   before BIC margins mean anything on real tails.

## 10. Unclear / possibly wrong

1. **The menu comparison has a fault that only real data showed.** `waa_pooled.py` profiles the start W0 of F1 and F4
   without a bound and keeps W0 >= 0.001 dB for F3, F5, F6 and the free curve. I report both readings. `F1+` is ONE
   likelihood evaluation of the F3 code at n = 1 with F1's tau and s (or F3's own fit where it ended at n = 1), not a
   re-optimised fit; F4 has no like-for-like version. The fake-data steps (3b-1b, 3b-1c, 3b-3c) used the same code;
   their tails never rise, so I expect no effect there, but I did not recheck. I did not change the library.
2. **"Falls with slope -0.87 +- 0.1" is the requested form, but +-0.1 is too optimistic here.** It comes from fake
   data whose curves were well determined. The honest spread on real data is the range over the sets (-1.60 to
   -0.05) and the bootstrap refits (-1.43 to -0.82). The same holds for "+-20 %" on the level.
3. **The free-curve fit may not be at its optimum.** Start spread 800-1700 in cost; the chosen penalty differs between
   the sets (0.01, 0.1, 1.0); the budgets are those of 3b-3b, unchanged, as required. A longer search could move the
   bumps and the slope.
4. **The dip of k near 0.9 dB is suspicious.** In every set the knot values go high, high, low, high, low; the low one
   is the third knot (0.91 dB on the main set, 0.68 dB on the alternative set), just below each set's A_start
   threshold (1.0 and 0.6 dB). That looks like an effect of the knots or of selecting tails by their start level,
   not physics. Not tested.
5. **The noise model does not fit drying tails**: whitened variance ratio 4.1 (first hour), 2.1 (6 h), 3.5 without
   the rise-heavy events. So the likelihood is misspecified, BIC margins are overstated, and the "free curve beats
   every menu law by about 10000 BIC" line mostly says that tails differ from each other in ways no single law
   covers. The `noise_model_*` columns ignore that a fit absorbs a start, an offset and a speed per tail, so the true
   gap is a little larger than shown.
6. **"Something still adds attenuation after t0_trace" is my interpretation, not a result.** The rise-heavy flag was
   set in step 2c from the rise rule, before any model fit, so dropping those events is not circular, but it is a
   choice about the data that the learners must make. It also removes 9 of the 24 events.
7. **The model-free table is rough**: it ignores the per-tail offset (it uses A, not W), smooths over 7 min, and bins
   by the same noisy A that enters the derivative. It supports the sign of the trend, not the numbers. I added it
   after seeing the curve; it is not one of the tools tested on fake data.
8. **Jump 2**: the E-SINDy terms come from the tested tool, but I refit their coefficients by least squares (the
   bagged medians, term by term, gave a formula 112 % off the curve when all five terms are kept). The library
   cannot write a falling k except with negative coefficients, and PySR cannot write max(). PySR's Pareto front
   starts with trivial formulas; its three simplest are listed as asked.
9. **Alternative set**: 119 of its sub-links have no control tails, so no noise parameters; they got the medians of
   the table (my choice). Its "F6 by 216 BIC" should be read with that in mind.
10. **Check 2**: without E20150826_0520 the slope is -1.60 and the menu is a tie (F1+ by 0.9); without E20150628_0530
    F6 leads by 106. So single events move both the curve and the menu a lot: 3 events hold 246 of 429 tails.
11. **Check 3 (6 h)**: F4 beats F1 by 5.5 BIC with tau_s = 16000 min and phi = 0.015, i.e. F1 plus a small constant
    floor. That is the offset trade-off seen in 3b-1, not a second layer. The lag-60 autocorrelation of check 4 can
    only come from this window.
12. **Check 5**: the three "recalibrated" truths of 3b-1c T3 were tuned to exactly these real numbers, so "inside" for
    m1 there mostly repeats 3a-fix (an exponential with spread reproduces the median curve). The real slope is outside
    every fake truth, which means the fake data did not contain the real behaviour: the tools were not tested on
    tails that plateau or rise.
13. **Cross-validation halves**: events are split 12 / 12 at random; with 3 big events the halves are unequal in
    tails. Seed fixed in `waa.real4.cv_seed`; not varied.
14. **Key figure**: `step4_k_real.png` is saved as diag (the module has 6 key outputs). If it is promoted, the output to
    demote is `outputs/tables/step3b1b_recovery.csv` (as proposed in 3b-2 and 3b-3b). My advice is to wait: the curve
    is not yet a result to show as "the law".
