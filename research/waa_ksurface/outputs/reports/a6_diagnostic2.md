# A6 diagnostic 2: the drying constant against the held k_wet

Written by `research/waa_ksurface/a6_diagnostic2.py`. Truth (a): 1 / k_wet = 5 min, 1 / k_dry = 15 min; the fake data of the A6 run. All links, radar input, K0 in both regimes, 1 / k_wet held. The one diagnostic allowed after the failed check; nothing was changed on the basis of it. Fake data only.

| extension | tau_wet_held_min | tau_dry_fit_min | k_dry_over_truth | within_25pct | cv_error_db | train_loss | F_knots_fitted_over_true_min | F_knots_fitted_over_true_max |
|---|---|---|---|---|---|---|---|---|
| held at the truth | 2.5 | 10.9377 | 1.3714 | False | 0.2779 | 0.0372 |  |  |
| held at the truth | 3.73 | 11.602 | 1.2929 | False | 0.2733 | 0.0361 |  |  |
| held at the truth | 5.0 | 12.2809 | 1.2214 | True | 0.27 | 0.0354 |  |  |
| held at the truth | 7.5 | 13.4724 | 1.1134 | True | 0.2667 | 0.0348 |  |  |
| held at the truth | 9.91 | 14.4583 | 1.0375 | True | 0.2659 | 0.0348 |  |  |
| free | 2.5 | 12.5473 | 1.1955 | True | 0.2438 | 0.0281 | 0.7732 | 0.8947 |
| free | 3.73 | 13.1897 | 1.1373 | True | 0.2412 | 0.0275 | 0.7719 | 0.9179 |
| free | 5.0 | 13.8822 | 1.0805 | True | 0.24 | 0.0273 | 0.7767 | 0.9245 |
| free | 7.5 | 15.1391 | 0.9908 | True | 0.2406 | 0.0275 | 0.7533 | 0.9258 |
| free | 9.91 | 16.2073 | 0.9255 | True | 0.2429 | 0.028 | 0.7293 | 0.9206 |

`k_dry_over_truth` must lie in 0.75–1.25 to be within 25 % (1 / k_dry between 12.0 and 20.0 min). 3.73 and 9.91 min are the gauge-input and radar-input estimates of 1 / k_wet in the failed run. F knots of the free fits: `research/waa_ksurface/outputs/tables/a6_diagnostic2_knots.csv`. (`research/waa_ksurface/outputs/tables/a6_diagnostic2.csv`)

