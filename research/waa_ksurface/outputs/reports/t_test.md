# The final TEST read: the frozen equation against the existing models

Written by `research/waa_ksurface/t_test.py` under `docs/waa_ksurface/test_prereg.md` (hashes in `t_test_freeze.json`). TEST read once by this branch (the second read of TEST overall; the first was Part 1's, for its own model). Nothing refitted.

## Data

| set | built | at_28_38 | dropped_no_a_link | kept | links | events |
|---|---|---|---|---|---|---|
| known | 3355 | 3050 | 0 | 3050 | 187 | 18 |
| new | 1010 | 923 | 0 | 923 | 53 | 17 |

Events: 18; with a gauge ratio (known links): 4. Run time 1.6 min.

## M1: RMSE of the 5 min A - c (paired event bootstrap, 1000 draws, seed 4100)

| set | links | model | RMSE (dB) | values | events | link-events | 90 % (dB) |
|---|---|---|---|---|---|---|---|
| known links | all links | ours | 0.7207 | 472142 | 18 | 3050 | 0.5120 to 0.8383 |
| known links | all links | Schleiss b | 0.7112 | 472142 | 18 | 3050 | 0.5250 to 0.8187 |
| known links | all links | Pastorek b | 0.7283 | 472142 | 18 | 3050 | 0.5288 to 0.8392 |
| known links | all links | Part 1 equation | 0.7023 | 472142 | 18 | 3050 | 0.5066 to 0.8121 |
| known links | all links | no WAA | 0.9239 | 472142 | 18 | 3050 | 0.6480 to 1.0805 |
| known links | all links | no memory | 0.7537 | 472142 | 18 | 3050 | 0.5501 to 0.8654 |
| known links | links < 1 km | ours | 0.5854 | 94583 | 14 | 563 | 0.4370 to 0.6664 |
| known links | links < 1 km | Schleiss b | 0.5755 | 94583 | 14 | 563 | 0.4530 to 0.6432 |
| known links | links < 1 km | Pastorek b | 0.5946 | 94583 | 14 | 563 | 0.4499 to 0.6674 |
| known links | links < 1 km | Part 1 equation | 0.5637 | 94583 | 14 | 563 | 0.4260 to 0.6351 |
| known links | links < 1 km | no WAA | 0.8562 | 94583 | 14 | 563 | 0.5879 to 1.0000 |
| known links | links < 1 km | no memory | 0.6175 | 94583 | 14 | 563 | 0.4733 to 0.6906 |
| new links | all links | ours | 0.764 | 146015 | 17 | 923 | 0.5179 to 0.9103 |
| new links | all links | Schleiss b | 0.7497 | 146015 | 17 | 923 | 0.5235 to 0.8786 |
| new links | all links | Pastorek b | 0.7916 | 146015 | 17 | 923 | 0.5364 to 0.9414 |
| new links | all links | Part 1 equation | 0.7559 | 146015 | 17 | 923 | 0.5159 to 0.8963 |
| new links | all links | no WAA | 0.9581 | 146015 | 17 | 923 | 0.6510 to 1.1438 |
| new links | all links | no memory | 0.7939 | 146015 | 17 | 923 | 0.5527 to 0.9333 |
| new links | links < 1 km | ours | 0.5684 | 32960 | 12 | 190 | 0.4527 to 0.6280 |
| new links | links < 1 km | Schleiss b | 0.5967 | 32960 | 12 | 190 | 0.4812 to 0.6634 |
| new links | links < 1 km | Pastorek b | 0.5775 | 32960 | 12 | 190 | 0.4753 to 0.6347 |
| new links | links < 1 km | Part 1 equation | 0.5637 | 32960 | 12 | 190 | 0.4536 to 0.6195 |
| new links | links < 1 km | no WAA | 0.8275 | 32960 | 12 | 190 | 0.6005 to 0.9487 |
| new links | links < 1 km | no memory | 0.595 | 32960 | 12 | 190 | 0.4896 to 0.6554 |
| gauge-ratio events | all links | ours | 0.8453 | 293322 | 4 | 1335 | 0.6314 to 0.9627 |
| gauge-ratio events | all links | Schleiss b | 0.8234 | 293322 | 4 | 1335 | 0.6215 to 0.9535 |
| gauge-ratio events | all links | Pastorek b | 0.8525 | 293322 | 4 | 1335 | 0.6557 to 0.9631 |
| gauge-ratio events | all links | Part 1 equation | 0.8186 | 293322 | 4 | 1335 | 0.6174 to 0.9298 |
| gauge-ratio events | all links | no WAA | 1.1065 | 293322 | 4 | 1335 | 0.9034 to 1.2304 |
| gauge-ratio events | all links | no memory | 0.8833 | 293322 | 4 | 1335 | 0.6921 to 0.9970 |
| gauge-ratio events | links < 1 km | ours | 0.6676 | 61944 | 4 | 268 | 0.5424 to 0.7380 |
| gauge-ratio events | links < 1 km | Schleiss b | 0.6457 | 61944 | 4 | 268 | 0.5140 to 0.7132 |
| gauge-ratio events | links < 1 km | Pastorek b | 0.6794 | 61944 | 4 | 268 | 0.5743 to 0.7376 |
| gauge-ratio events | links < 1 km | Part 1 equation | 0.639 | 61944 | 4 | 268 | 0.5193 to 0.6973 |
| gauge-ratio events | links < 1 km | no WAA | 1.0112 | 61944 | 4 | 268 | 0.8340 to 1.1159 |
| gauge-ratio events | links < 1 km | no memory | 0.7063 | 61944 | 4 | 268 | 0.6041 to 0.7650 |

### Paired differences and the pre-registered readings (the only claims)

| set | links | model | RMSE(model) - RMSE(ours) (dB) | 90 % lo | 90 % hi | reading |
|---|---|---|---|---|---|---|
| known links | all links | Schleiss b | -0.0095 | -0.0357 | 0.0338 | tie |
| known links | all links | Pastorek b | 0.0075 | -0.0058 | 0.0242 | tie |
| known links | all links | Part 1 equation | -0.0184 | -0.0286 | 0.0002 | tie |
| known links | all links | no WAA | 0.2032 | 0.1107 | 0.2483 | ours better than no WAA |
| known links | all links | no memory | 0.033 | 0.017 | 0.0533 | ours better than no memory |
| known links | links < 1 km | Schleiss b | -0.0099 | -0.0378 | 0.0398 | tie |
| known links | links < 1 km | Pastorek b | 0.0091 | -0.0054 | 0.0292 | tie |
| known links | links < 1 km | Part 1 equation | -0.0217 | -0.0344 | -0.0007 | Part 1 equation better than ours |
| known links | links < 1 km | no WAA | 0.2708 | 0.1283 | 0.3424 | ours better than no WAA |
| known links | links < 1 km | no memory | 0.0321 | 0.0116 | 0.0605 | ours better than no memory |
| new links | all links | Schleiss b | -0.0142 | -0.041 | 0.0295 | tie |
| new links | all links | Pastorek b | 0.0277 | 0.0172 | 0.0339 | ours better than Pastorek b |
| new links | all links | Part 1 equation | -0.0081 | -0.0151 | 0.0045 | tie |
| new links | all links | no WAA | 0.1941 | 0.1037 | 0.2381 | ours better than no WAA |
| new links | all links | no memory | 0.03 | 0.0155 | 0.0505 | ours better than no memory |
| new links | links < 1 km | Schleiss b | 0.0283 | 0.0066 | 0.05 | ours better than Schleiss b |
| new links | links < 1 km | Pastorek b | 0.0091 | -0.009 | 0.0354 | tie |
| new links | links < 1 km | Part 1 equation | -0.0047 | -0.0133 | 0.0089 | tie |
| new links | links < 1 km | no WAA | 0.2591 | 0.1265 | 0.3319 | ours better than no WAA |
| new links | links < 1 km | no memory | 0.0265 | 0.0054 | 0.0583 | ours better than no memory |
| gauge-ratio events | all links | Schleiss b | -0.0218 | -0.0558 | 0.0425 | tie |
| gauge-ratio events | all links | Pastorek b | 0.0072 | -0.0092 | 0.0328 | tie |
| gauge-ratio events | all links | Part 1 equation | -0.0266 | -0.0359 | -0.0057 | Part 1 equation better than ours |
| gauge-ratio events | all links | no WAA | 0.2613 | 0.1954 | 0.2902 | ours better than no WAA |
| gauge-ratio events | all links | no memory | 0.038 | 0.0119 | 0.0714 | ours better than no memory |
| gauge-ratio events | links < 1 km | Schleiss b | -0.0218 | -0.0558 | 0.0448 | tie |
| gauge-ratio events | links < 1 km | Pastorek b | 0.0118 | -0.0071 | 0.0446 | tie |
| gauge-ratio events | links < 1 km | Part 1 equation | -0.0286 | -0.0408 | 0.0 | tie |
| gauge-ratio events | links < 1 km | no WAA | 0.3437 | 0.2247 | 0.3929 | ours better than no WAA |
| gauge-ratio events | links < 1 km | no memory | 0.0388 | 0.0076 | 0.0882 | ours better than no memory |

## M2: after-rain error (a CHECK; dry-shape caveat; no claim)

| set | lag (min) | values | events | measured median (dB) | ours | Schleiss b | Pastorek b | Part 1 equation | no WAA | no memory |
|---|---|---|---|---|---|---|---|---|---|---|
| known links | 5 | 2478 | 15 | 0.37000045776367185 | -0.084 (-0.431 to +0.114) | -0.172 (-0.531 to +0.033) | -0.236 (-0.581 to -0.026) | -0.050 (-0.405 to +0.147) | -0.370 (-0.709 to -0.170) | -0.370 (-0.702 to -0.149) |
| known links | 10 | 2464 | 15 | 0.26999969482421876 | -0.041 (-0.282 to +0.117) | -0.125 (-0.371 to +0.037) | -0.257 (-0.498 to -0.102) | -0.038 (-0.276 to +0.121) | -0.270 (-0.510 to -0.120) | -0.270 (-0.510 to -0.120) |
| known links | 20 | 2446 | 15 | 0.1676666259765625 | +0.017 (-0.116 to +0.126) | -0.073 (-0.209 to +0.042) | -0.168 (-0.300 to -0.060) | -0.043 (-0.172 to +0.065) | -0.168 (-0.300 to -0.060) | -0.168 (-0.300 to -0.060) |
| new links | 5 | 742 | 12 | 0.39499931335449223 | -0.102 (-0.479 to +0.154) | -0.195 (-0.579 to +0.070) | -0.276 (-0.637 to -0.007) | -0.073 (-0.463 to +0.184) | -0.395 (-0.755 to -0.150) | -0.395 (-0.750 to -0.140) |
| new links | 10 | 736 | 12 | 0.2899993896484375 | -0.061 (-0.359 to +0.132) | -0.146 (-0.447 to +0.059) | -0.290 (-0.577 to -0.090) | -0.057 (-0.357 to +0.146) | -0.290 (-0.588 to -0.110) | -0.290 (-0.588 to -0.110) |
| new links | 20 | 726 | 12 | 0.18833389282226562 | -0.003 (-0.171 to +0.136) | -0.092 (-0.268 to +0.053) | -0.188 (-0.353 to -0.050) | -0.058 (-0.219 to +0.079) | -0.188 (-0.353 to -0.050) | -0.188 (-0.353 to -0.050) |
| gauge-ratio events | 5 | 1196 | 4 | 0.6300016784667953 | -0.348 (-0.882 to -0.012) | -0.441 (-0.994 to -0.113) | -0.510 (-1.037 to -0.178) | -0.272 (-0.811 to +0.059) | -0.630 (-1.156 to -0.280) | -0.630 (-1.145 to -0.280) |
| gauge-ratio events | 10 | 1180 | 4 | 0.4500007629394531 | -0.223 (-0.580 to +0.037) | -0.315 (-0.694 to -0.059) | -0.450 (-0.790 to -0.180) | -0.185 (-0.563 to +0.067) | -0.450 (-0.807 to -0.180) | -0.450 (-0.807 to -0.180) |
| gauge-ratio events | 20 | 1162 | 4 | 0.2844988250732413 | -0.101 (-0.318 to +0.073) | -0.198 (-0.431 to -0.043) | -0.284 (-0.500 to -0.105) | -0.130 (-0.350 to +0.010) | -0.284 (-0.500 to -0.105) | -0.284 (-0.500 to -0.105) |

## M3: the length test (rule R2) on every model's predictions (no claim)

| result |
|---|
| no cell is both an E4 cell and eligible on the measurements of any set |

