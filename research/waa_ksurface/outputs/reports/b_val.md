# Block B: the frozen equation against the existing models on VAL

Written by `research/waa_ksurface/b_val.py` under `docs/waa_ksurface/block_b_prereg.md` (frozen in 37470fa; recipe fits frozen in the second record `a8_block_b_freeze2.json`, re-hashed after the fix of Deviation 16 in `a8_block_b_freeze3.json`). **VAL was read twice** (Deviation 16: the first run stopped in M3 before any result was saved or seen; the results below are from the second run, computed once; this report is written from its saved results). TEST not opened.

VAL link-events at 28 / 38 GHz scored: 2944 of 18 events (5 with a gauge ratio); dropped because a model cannot predict them (no stored a_link): 0. Run time 1.2 min.

## M1: RMSE of the 5 min A - c (paired event bootstrap, 1000 draws, seed 4100)

| set | links | model | RMSE (dB) | values | events | 90 % (dB) |
|---|---|---|---|---|---|---|
| S-all | all links | ours | 0.6385 | 391024 | 18 | 0.5186 to 0.7170 |
| S-all | all links | Schleiss b | 0.6628 | 391024 | 18 | 0.5503 to 0.7391 |
| S-all | all links | Pastorek b | 0.6315 | 391024 | 18 | 0.5156 to 0.7088 |
| S-all | all links | Part 1 equation | 0.6273 | 391024 | 18 | 0.5110 to 0.7016 |
| S-all | all links | no WAA | 0.7434 | 391024 | 18 | 0.5954 to 0.8525 |
| S-all | all links | no memory | 0.6554 | 391024 | 18 | 0.5329 to 0.7372 |
| S-all | links < 1 km | ours | 0.4953 | 78881 | 14 | 0.4408 to 0.5321 |
| S-all | links < 1 km | Schleiss b | 0.546 | 78881 | 14 | 0.4800 to 0.5943 |
| S-all | links < 1 km | Pastorek b | 0.495 | 78881 | 14 | 0.4309 to 0.5340 |
| S-all | links < 1 km | Part 1 equation | 0.4839 | 78881 | 14 | 0.4276 to 0.5202 |
| S-all | links < 1 km | no WAA | 0.6354 | 78881 | 14 | 0.5289 to 0.7217 |
| S-all | links < 1 km | no memory | 0.5191 | 78881 | 14 | 0.4528 to 0.5631 |
| S-ratio | all links | ours | 0.7055 | 282587 | 5 | 0.5795 to 0.8027 |
| S-ratio | all links | Schleiss b | 0.7302 | 282587 | 5 | 0.6124 to 0.8203 |
| S-ratio | all links | Pastorek b | 0.698 | 282587 | 5 | 0.5790 to 0.7917 |
| S-ratio | all links | Part 1 equation | 0.6902 | 282587 | 5 | 0.5698 to 0.7833 |
| S-ratio | all links | no WAA | 0.8278 | 282587 | 5 | 0.6797 to 0.9521 |
| S-ratio | all links | no memory | 0.7273 | 282587 | 5 | 0.6079 to 0.8357 |
| S-ratio | links < 1 km | ours | 0.5282 | 61015 | 5 | 0.4845 to 0.5681 |
| S-ratio | links < 1 km | Schleiss b | 0.589 | 61015 | 5 | 0.5373 to 0.6458 |
| S-ratio | links < 1 km | Pastorek b | 0.5303 | 61015 | 5 | 0.4798 to 0.5679 |
| S-ratio | links < 1 km | Part 1 equation | 0.5159 | 61015 | 5 | 0.4713 to 0.5557 |
| S-ratio | links < 1 km | no WAA | 0.6927 | 61015 | 5 | 0.5919 to 0.7898 |
| S-ratio | links < 1 km | no memory | 0.5588 | 61015 | 5 | 0.5077 to 0.6049 |

### Paired differences and the pre-registered readings (the only claims)

| set | links | model | RMSE(model) - RMSE(ours) (dB) | 90 % lo | 90 % hi | reading |
|---|---|---|---|---|---|---|
| S-all | all links | Schleiss b | 0.0243 | 0.0164 | 0.0355 | ours better than Schleiss b |
| S-all | all links | Pastorek b | -0.007 | -0.0144 | 0.0029 | tie |
| S-all | all links | Part 1 equation | -0.0112 | -0.0206 | -0.0009 | Part 1 equation better than ours |
| S-all | all links | no WAA | 0.1049 | 0.0629 | 0.1434 | ours better than no WAA |
| S-all | all links | no memory | 0.0169 | 0.0056 | 0.0283 | ours better than no memory |
| S-all | links < 1 km | Schleiss b | 0.0506 | 0.0257 | 0.0721 | ours better than Schleiss b |
| S-all | links < 1 km | Pastorek b | -0.0004 | -0.0136 | 0.0141 | tie |
| S-all | links < 1 km | Part 1 equation | -0.0114 | -0.0215 | -0.0004 | Part 1 equation better than ours |
| S-all | links < 1 km | no WAA | 0.1401 | 0.0712 | 0.1992 | ours better than no WAA |
| S-all | links < 1 km | no memory | 0.0237 | 0.0062 | 0.0393 | ours better than no memory |
| S-ratio | all links | Schleiss b | 0.0247 | 0.0169 | 0.0368 | ours better than Schleiss b |
| S-ratio | all links | Pastorek b | -0.0075 | -0.0178 | 0.0042 | tie |
| S-ratio | all links | Part 1 equation | -0.0153 | -0.0277 | -0.0045 | Part 1 equation better than ours |
| S-ratio | all links | no WAA | 0.1223 | 0.0746 | 0.17 | ours better than no WAA |
| S-ratio | all links | no memory | 0.0218 | 0.011 | 0.0343 | ours better than no memory |
| S-ratio | links < 1 km | Schleiss b | 0.0608 | 0.0352 | 0.0852 | ours better than Schleiss b |
| S-ratio | links < 1 km | Pastorek b | 0.0021 | -0.0134 | 0.0195 | tie |
| S-ratio | links < 1 km | Part 1 equation | -0.0123 | -0.0257 | 0.0015 | tie |
| S-ratio | links < 1 km | no WAA | 0.1645 | 0.0888 | 0.2359 | ours better than no WAA |
| S-ratio | links < 1 km | no memory | 0.0306 | 0.0137 | 0.0481 | ours better than no memory |

## M2: after-rain error (a CHECK, not a measurement of the drying)

median(model) - median(measured) at the bin ending 5 / 10 / 20 min after the path went dry (review check (5) windows), 90 % from the same event draws. Caveat (dry-shape control): the data-layer offsets and the selection of dry-path windows already distort a plain exponential drying (ratio 0.30 / 0.44 for a K0 truth), so these numbers are not a measurement of the drying dynamics. No claim.

| set | lag (min) | values | events | measured median (dB) | ours | Schleiss b | Pastorek b | Part 1 equation | no WAA | no memory |
|---|---|---|---|---|---|---|---|---|---|---|
| S-all | 5 | 2387 | 14 | 0.21999969482421874 | +0.063 (-0.069 to +0.181) | -0.041 (-0.163 to +0.089) | -0.113 (-0.224 to +0.005) | +0.090 (-0.041 to +0.215) | -0.220 (-0.320 to -0.120) | -0.220 (-0.320 to -0.110) |
| S-all | 10 | 2382 | 14 | 0.13333511352539062 | +0.094 (-0.018 to +0.183) | -0.005 (-0.114 to +0.096) | -0.133 (-0.230 to -0.053) | +0.086 (-0.037 to +0.188) | -0.133 (-0.230 to -0.060) | -0.133 (-0.230 to -0.060) |
| S-all | 20 | 2377 | 14 | 0.059999084472656256 | +0.124 (+0.053 to +0.176) | +0.026 (-0.052 to +0.093) | -0.060 (-0.127 to -0.010) | +0.054 (-0.026 to +0.117) | -0.060 (-0.127 to -0.010) | -0.060 (-0.127 to -0.010) |
| S-ratio | 5 | 1531 | 5 | 0.31000289916992185 | -0.038 (-0.220 to +0.128) | -0.140 (-0.314 to +0.068) | -0.219 (-0.390 to -0.117) | +0.024 (-0.199 to +0.232) | -0.310 (-0.470 to -0.213) | -0.310 (-0.470 to -0.213) |
| S-ratio | 10 | 1526 | 5 | 0.2100006103515625 | +0.015 (-0.099 to +0.151) | -0.082 (-0.206 to +0.089) | -0.210 (-0.302 to -0.107) | +0.022 (-0.135 to +0.214) | -0.210 (-0.302 to -0.107) | -0.210 (-0.302 to -0.107) |
| S-ratio | 20 | 1525 | 5 | 0.1 | +0.083 (-0.027 to +0.174) | -0.014 (-0.136 to +0.108) | -0.100 (-0.200 to -0.020) | +0.026 (-0.110 to +0.148) | -0.100 (-0.200 to -0.020) | -0.100 (-0.200 to -0.020) |

## M3: the length test (E4 candidate (ii), rule R2) on every model's predictions

Cells: eligible in E4 and on the VAL measurements. A model's cell is inside if both its intercept and slope lie in the measured 90 % event-bootstrap interval. No claim.

| set | result |
|---|---|
| S-all and S-ratio | no cell is both an E4 cell and eligible on the VAL measurements: M3 has nothing to test |



