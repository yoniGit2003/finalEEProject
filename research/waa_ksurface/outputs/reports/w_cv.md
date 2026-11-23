# WAA-only comparison: the same path law for every model, leave-one-event-out on TRAIN

Written by `research/waa_ksurface/w_cv.py` (decisions.md, "WAA-only comparison", logged before any code was run). 14 TRAIN learning events, 4563 link-events (28 / 38 GHz); held-out predictions only. VAL and TEST not opened. Disclosure: the forms were chosen on these events (ours) or on TRAIN + VAL (Part 1); only parameters are cross-validated.

## RMSE of the 5 min A - c (held-out; paired event bootstrap, 1000 draws, seed 4400)

| path | links | model | RMSE (dB) | values | 90 % (dB) |
|---|---|---|---|---|---|
| A | all links | ours | 0.9276 | 1354791 | 0.7623 to 1.0737 |
| A | all links | Part 1 equation | 0.9232 | 1354791 | 0.7610 to 1.0674 |
| A | all links | Schleiss | 1.0007 | 1354791 | 0.8379 to 1.1361 |
| A | all links | Pastorek | 0.9765 | 1354791 | 0.8068 to 1.1246 |
| A | all links | no memory | 0.9791 | 1354791 | 0.8125 to 1.1260 |
| A | all links | no WAA | 1.2327 | 1354791 | 1.0545 to 1.3953 |
| A | links < 1 km | ours | 0.7167 | 299593 | 0.5958 to 0.8279 |
| A | links < 1 km | Part 1 equation | 0.7137 | 299593 | 0.5918 to 0.8250 |
| A | links < 1 km | Schleiss | 0.802 | 299593 | 0.6773 to 0.9121 |
| A | links < 1 km | Pastorek | 0.7676 | 299593 | 0.6405 to 0.8784 |
| A | links < 1 km | no memory | 0.7688 | 299593 | 0.6430 to 0.8798 |
| A | links < 1 km | no WAA | 1.0617 | 299593 | 0.9133 to 1.1990 |
| B | all links | ours | 0.9353 | 1354791 | 0.7658 to 1.0857 |
| B | all links | Part 1 equation | 0.925 | 1354791 | 0.7587 to 1.0717 |
| B | all links | Schleiss | 0.9908 | 1354791 | 0.8315 to 1.1265 |
| B | all links | Pastorek | 0.9816 | 1354791 | 0.8098 to 1.1323 |
| B | all links | no memory | 0.9821 | 1354791 | 0.8130 to 1.1301 |
| B | all links | no WAA | 1.227 | 1354791 | 1.0508 to 1.3901 |
| B | links < 1 km | ours | 0.718 | 299593 | 0.5967 to 0.8305 |
| B | links < 1 km | Part 1 equation | 0.7127 | 299593 | 0.5916 to 0.8246 |
| B | links < 1 km | Schleiss | 0.7936 | 299593 | 0.6711 to 0.9020 |
| B | links < 1 km | Pastorek | 0.7678 | 299593 | 0.6417 to 0.8796 |
| B | links < 1 km | no memory | 0.7685 | 299593 | 0.6446 to 0.8778 |
| B | links < 1 km | no WAA | 1.0565 | 299593 | 0.9087 to 1.1938 |
| A+a | all links | ours | 0.9276 | 1354791 | 0.7623 to 1.0737 |
| A+a | all links | Part 1 equation | 0.9195 | 1354791 | 0.7555 to 1.0662 |
| A+a | all links | Schleiss | 0.9955 | 1354791 | 0.8337 to 1.1299 |
| A+a | all links | Pastorek | 0.9725 | 1354791 | 0.8012 to 1.1233 |
| A+a | all links | no memory | 0.9791 | 1354791 | 0.8125 to 1.1260 |
| A+a | all links | no WAA | 1.2327 | 1354791 | 1.0545 to 1.3953 |
| A+a | links < 1 km | ours | 0.7167 | 299593 | 0.5958 to 0.8279 |
| A+a | links < 1 km | Part 1 equation | 0.7078 | 299593 | 0.5827 to 0.8184 |
| A+a | links < 1 km | Schleiss | 0.7967 | 299593 | 0.6727 to 0.9086 |
| A+a | links < 1 km | Pastorek | 0.7626 | 299593 | 0.6344 to 0.8749 |
| A+a | links < 1 km | no memory | 0.7688 | 299593 | 0.6430 to 0.8798 |
| A+a | links < 1 km | no WAA | 1.0617 | 299593 | 0.9133 to 1.1990 |

### Paired RMSE(model) - RMSE(ours) and readings (path A = main; B = Part 1's straight path; A+a = competitors with a_link)

| path | links | model | RMSE(model) - RMSE(ours) (dB) | 90 % lo | 90 % hi | reading |
|---|---|---|---|---|---|---|
| A | all links | Part 1 equation | -0.0044 | -0.0113 | 0.0018 | tie |
| A | all links | Schleiss | 0.0731 | 0.0377 | 0.1102 | ours better than Schleiss |
| A | all links | Pastorek | 0.0489 | 0.0365 | 0.0632 | ours better than Pastorek |
| A | all links | no memory | 0.0514 | 0.0386 | 0.0657 | ours better than no memory |
| A | all links | no WAA | 0.3051 | 0.2203 | 0.4085 | ours better than no WAA |
| A | links < 1 km | Part 1 equation | -0.003 | -0.0098 | 0.0026 | tie |
| A | links < 1 km | Schleiss | 0.0853 | 0.0643 | 0.1099 | ours better than Schleiss |
| A | links < 1 km | Pastorek | 0.0509 | 0.0368 | 0.0661 | ours better than Pastorek |
| A | links < 1 km | no memory | 0.0522 | 0.0378 | 0.0676 | ours better than no memory |
| A | links < 1 km | no WAA | 0.3451 | 0.2657 | 0.4407 | ours better than no WAA |
| B | all links | Part 1 equation | -0.0103 | -0.0181 | -0.0026 | Part 1 equation better than ours |
| B | all links | Schleiss | 0.0555 | 0.0242 | 0.0887 | ours better than Schleiss |
| B | all links | Pastorek | 0.0463 | 0.033 | 0.0602 | ours better than Pastorek |
| B | all links | no memory | 0.0468 | 0.0333 | 0.0605 | ours better than no memory |
| B | all links | no WAA | 0.2917 | 0.2078 | 0.3935 | ours better than no WAA |
| B | links < 1 km | Part 1 equation | -0.0053 | -0.0119 | 0.0009 | tie |
| B | links < 1 km | Schleiss | 0.0756 | 0.0567 | 0.0983 | ours better than Schleiss |
| B | links < 1 km | Pastorek | 0.0498 | 0.0358 | 0.065 | ours better than Pastorek |
| B | links < 1 km | no memory | 0.0504 | 0.0363 | 0.0664 | ours better than no memory |
| B | links < 1 km | no WAA | 0.3385 | 0.2617 | 0.4319 | ours better than no WAA |
| A+a | all links | Part 1 equation | -0.0081 | -0.0148 | -0.0024 | Part 1 equation better than ours |
| A+a | all links | Schleiss | 0.0679 | 0.0343 | 0.1028 | ours better than Schleiss |
| A+a | all links | Pastorek | 0.0449 | 0.0333 | 0.0573 | ours better than Pastorek |
| A+a | all links | no memory | 0.0514 | 0.0386 | 0.0657 | ours better than no memory |
| A+a | all links | no WAA | 0.3051 | 0.2203 | 0.4085 | ours better than no WAA |
| A+a | links < 1 km | Part 1 equation | -0.0088 | -0.0178 | -0.0007 | Part 1 equation better than ours |
| A+a | links < 1 km | Schleiss | 0.0801 | 0.0596 | 0.1019 | ours better than Schleiss |
| A+a | links < 1 km | Pastorek | 0.0459 | 0.0309 | 0.0621 | ours better than Pastorek |
| A+a | links < 1 km | no memory | 0.0522 | 0.0378 | 0.0676 | ours better than no memory |
| A+a | links < 1 km | no WAA | 0.3451 | 0.2657 | 0.4407 | ours better than no WAA |

**Ranking by RMSE (all links):** path A: Part 1 equation, ours, Pastorek, no memory, Schleiss, no WAA; path B: Part 1 equation, ours, Pastorek, no memory, Schleiss, no WAA. **The same under A and B.** Readings that change: Part 1 equation (all links): A 'tie' vs B 'Part 1 equation better than ours'.

## After-rain error (path A; model - measured, median with 90 %; paired |e| comparison)

| item | lag (min) | windows | measured median (dB) | ours | Part 1 equation | Schleiss | Pastorek | no memory | no WAA | Part 1 equation: |e| - |e ours| | Schleiss: |e| - |e ours| | Pastorek: |e| - |e ours| | no memory: |e| - |e ours| | no WAA: |e| - |e ours| |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| after-rain error | 5 | 3502 | 0.25999984741210935 | +0.028 (-0.060 to +0.104) | +0.065 (-0.002 to +0.117) | +0.052 (-0.030 to +0.135) | -0.170 (-0.280 to -0.080) | -0.210 (-0.300 to -0.120) | -0.250 (-0.330 to -0.168) | -0.038 to +0.058 (tie) | -0.055 to +0.101 (tie) | -0.024 to +0.234 (tie) | +0.014 to +0.259 (ours closer) | +0.063 to +0.288 (ours closer) |
| after-rain error | 10 | 3466 | 0.17333297729492114 | +0.060 (-0.013 to +0.112) | +0.051 (+0.012 to +0.083) | +0.050 (+0.000 to +0.114) | -0.144 (-0.220 to -0.090) | -0.153 (-0.234 to -0.100) | -0.170 (-0.250 to -0.110) | -0.045 to +0.029 (tie) | -0.062 to +0.069 (tie) | -0.019 to +0.207 (tie) | -0.012 to +0.220 (tie) | -0.000 to +0.238 (tie) |
| after-rain error | 20 | 3430 | 0.09999847412109375 | +0.082 (+0.049 to +0.106) | +0.030 (+0.003 to +0.056) | +0.027 (-0.000 to +0.060) | -0.090 (-0.133 to -0.060) | -0.090 (-0.137 to -0.060) | -0.100 (-0.140 to -0.060) | -0.076 to -0.019 (closer than ours) | -0.084 to -0.009 (closer than ours) | -0.048 to +0.084 (tie) | -0.047 to +0.087 (tie) | -0.044 to +0.091 (tie) |

## Peak error (path A; link-events with a measured peak >= 2.0 dB)

| item | links | link-events | measured peak median (dB) | ours | Part 1 equation | Schleiss | Pastorek | no memory | no WAA | Part 1 equation: |e| - |e ours| | Schleiss: |e| - |e ours| | Pastorek: |e| - |e ours| | no memory: |e| - |e ours| | no WAA: |e| - |e ours| |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| peak error | links < 1 km | 798 | 5.248333740234375 | -2.783 (-3.698 to -2.195) | -2.946 (-3.824 to -2.350) | -3.972 (-4.712 to -3.066) | -2.678 (-3.438 to -2.300) | -2.807 (-3.655 to -2.359) | -4.783 (-5.581 to -3.960) | +0.058 to +0.289 (ours closer) | +0.654 to +1.459 (ours closer) | -0.456 to +0.145 (tie) | -0.187 to +0.282 (tie) | +1.392 to +2.481 (ours closer) |
| peak error | links > 3 km | 1279 | 7.026666259765625 | -3.676 (-4.497 to -3.046) | -3.726 (-4.588 to -3.098) | -4.302 (-5.210 to -3.513) | -3.425 (-4.397 to -2.840) | -3.439 (-4.397 to -2.851) | -5.174 (-6.037 to -4.438) | -0.031 to +0.126 (tie) | +0.333 to +0.877 (ours closer) | -0.368 to -0.077 (closer than ours) | -0.321 to -0.044 (closer than ours) | +1.167 to +1.864 (ours closer) |

## Figures (path A, held-out)

| rule | candidates | median held-out RMSE ours (dB) | event | sublink | L (km) | band | held-out RMSE ours (dB) | rank | figure |
|---|---|---|---|---|---|---|---|---|---|
| short | 859 | 0.5765 | E20150826_1620 | 10026_1 | 0.885 | 38 GHz | 0.5765 | 430 of 859 | research/waa_ksurface/outputs/figures/w_cv_short_E20150826_1620_10026_1.png |
| long | 975 | 0.7737 | E20150826_0520 | 10146_2 | 3.0468 | 28 GHz | 0.7737 | 488 of 975 | research/waa_ksurface/outputs/figures/w_cv_long_E20150826_0520_10146_2.png |

Saved predictions (total and WAA part per model, path and held-out bin): `research/.cache/waa_ksurface/w_cv/preds.npz` (cache, not committed). Fitted parameters per fold: `research/waa_ksurface/outputs/tables/w_cv_params.csv`. Run time 73.8 min (this pass).

