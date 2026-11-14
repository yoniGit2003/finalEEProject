# A6 diagnostic 3: why a k_dry that triples is worth 0.002 dB

Written by `research/waa_ksurface/a6_diagnostic3.py`. Truth (c): 1 / k_wet = 5 min, k_dry = (1 + W) / 30 per min; the fake data of the A6 final run. All links, radar input, 1 / k_wet held at its gauge-input estimate (3.54 min), extension free, drying K0 and drying K1 (lambda = 0.0001) fitted on all events. The one diagnostic allowed after the failed check; nothing was changed on the basis of it. Fake data only.

## (i) The K1 fit against the truth

| W_db | tau_true_min | tau_K1_min | tau_K0_min | K1_k_over_true | K1_within_25pct | share_of_domain_minutes_from_this_knot |
|---|---|---|---|---|---|---|
| 0.0 | 30.0 | 34.344 | 16.824 | 0.874 | True | 0.863 |
| 0.5 | 20.0 | 21.53 | 16.824 | 0.929 | True | 0.043 |
| 1.0 | 15.0 | 14.272 | 16.824 | 1.051 | True | 0.08 |
| 2.0 | 10.0 | 7.791 | 16.824 | 1.284 | False | 0.014 |
| 4.0 | 6.0 | 4.235 | 16.824 | 1.417 | False | 0.0 |
| 9.0 | 3.0 | 2.301 | 16.824 | 1.304 | False | 0.0 |

Share of the domain weight within 25 % of the true k_dry: K0 18 % (the fit the ladder returned), K1 97 % (the fit it did not return); the rule needs 80 %. (`research/waa_ksurface/outputs/tables/a6_diagnostic3_curve.csv`)

## (ii) The in-sample error of the two fits on parts of the data

| set | values | share_of_all_values | error_K0_db | error_K1_db | gain_db | share_of_the_total_loss_drop |
|---|---|---|---|---|---|---|
| all 5 min values | 1354791 | 1.0 | 0.2344 | 0.2319 | 0.0025 | 1.0 |
| links below 1 km | 299593 | 0.2211 | 0.2111 | 0.2073 | 0.0038 | 0.3 |
| links below 1 km, 60 min after the last wet bin | 10287 | 0.0076 | 0.1679 | 0.1614 | 0.0065 | 0.0139 |
| all, rainy part left out | 132677 | 0.0979 | 0.1514 | 0.1491 | 0.0023 | 0.0586 |
| links below 1 km, rainy part left out | 30588 | 0.0226 | 0.1534 | 0.1508 | 0.0027 | 0.0156 |
| links below 1 km, 60 min after the last wet bin, rainy part left out | 10287 | 0.0076 | 0.1679 | 0.1614 | 0.0065 | 0.0139 |

Error = sqrt(2 x mean Huber loss), dB; gain = K0 minus K1. "Rainy part left out" = only the 5 min bins after the last wet bin of the link-event, so the last row is the third row again (by construction). The held-out gain of K1 in the A6 final run was 0.0019 dB; the margin is 0.005 dB. (`research/waa_ksurface/outputs/tables/a6_diagnostic3.csv`)

