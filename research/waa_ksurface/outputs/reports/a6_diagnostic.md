# A6 diagnostic: which ingredient makes the fit invent a k(R, W)?

Written by `research/waa_ksurface/a6_diagnostic.py`. Truth (a): constant k, 1 / k = 5 min (wetting) and 15 min (drying); same seed as the A6 run. The one diagnostic allowed after the failed check; nothing was changed on the basis of it. Fake data only.

## The rungs returned by the two stages

| run | gauge_input | radar_input | noise | wetting_rung | drying_rung | share_wet_within_25pct | share_dry_within_25pct |
|---|---|---|---|---|---|---|---|
| D0 | exact | exact | off | K0 | K1 | 1.0 | 1.0 |
| D1 | exact | exact | on | K3 | K1 | 0.688 | 1.0 |
| D2 | with error | exact | on | K3 | K3 | 0.096 | 0.181 |
| D3 | exact | with error | on | K3 | K0 | 0.702 | 1.0 |
| A6 run | with error | with error | on | K3 | K1 | 0.139 | 0.121 |

True rung: K0 in both regimes. The shares are the part of the domain weight on which the returned k is within 25 % of the truth.

## Held-out loss: K0 against the best rung

| run | wet_cv_K0 | wet_cv_best | dry_cv_K0 | dry_cv_best |
|---|---|---|---|---|
| D0 | 0.0 | 0.0 | 5e-05 | 4.9e-05 |
| D1 | 0.010924 | 0.010903 | 0.011421 | 0.011417 |
| D2 | 0.029139 | 0.028525 | 0.013508 | 0.013313 |
| D3 | 0.010924 | 0.010903 | 0.036204 | 0.03554 |
| A6 run | 0.029139 | 0.028525 | 0.037596 | 0.036554 |

## A plain K0 / K0 fit (no richer rung offered): 1 / k in minutes

| run | K0_gauge_tau_wet | K0_gauge_tau_dry | K0_radar_tau_wet | K0_radar_tau_dry |
|---|---|---|---|---|
| D0 | 5.0 | 15.0 | 5.06 | 14.79 |
| D1 | 5.47 | 14.4 | 5.17 | 14.76 |
| D2 | 3.79 | 20.82 | 5.17 | 14.76 |
| D3 | 5.47 | 14.4 | 9.91 | 14.46 |
| A6 run | 3.79 | 20.82 | 9.91 | 14.46 |

Truth: 5 (wetting) and 15 (drying). `K0_gauge_*` = the gauge-input set with the gauge input; `K0_radar_*` = all links with the radar input.

(`research/waa_ksurface/outputs/tables/a6_diagnostic.csv`)

