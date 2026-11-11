# E3: the bin definitions of the length test, and their fake test

Written by `research/waa_ksurface/e3_bin_rule.py`. 4563 TRAIN link-events, 187 links, 14 events with a gauge ratio (28 and 38 GHz). The definitions and their reasons are in `docs/waa_ksurface/decisions.md` (step E3). VAL and TEST were not opened.

## Sizes measured on the real TRAIN data

| what | value | unit | from |
|---|---|---|---|
| white noise of a 5 min value (dry offset window) | 0.054 | dB | 4361 link-events |
| offset error of a link-event (dry 30 min before the rain) | 0.136 | dB | 4547 link-events |
| between-link sd of the short links' residuals | 0.241 | dB | 41 links, 539 link-storm values |
| within-link sd (storm to storm) | 0.592 | dB | same |

**n_short = max(5, ceil((0.241 / 0.15)²)) = 5** distinct links below 1 km, and as many of 1 km or more, for a cell to be eligible.

32 GHz has 0 links below 1 km, so by this rule it has no eligible cell: **32 GHz stays excluded from the branch for every model.**

## The fake test of the rules (memory and noise, no rain error)

Truth: F = 3·R / (R + 2) dB, P = 0.1·R / 0.2·R dB/km, one store with 5 min (wetting) and 15 min (drying). 5 replicates, medians. Pass (Deviation 3): intercept within 0.15 dB and slope within 0.15 dB / L_med = 0.0710 dB/km (L_med = 2.113 km, the median of the 187 links) in all judged cells but at most 1.

### Rule R1: 20 min within a factor 2 — 14 cells judged, 3 outside: FAIL

| band | class | R_cell | F_true | I | I_err | P_true | s | s_err_abs | links_short | I_ok | s_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.144 | 0.201 | 0.3 | 0.099 | 0.014 | 0.002 | -0.013 | 17.0 | True | True |
| 28 GHz | 0.2-0.5 | 0.331 | 0.426 | 0.521 | 0.096 | 0.033 | 0.017 | -0.016 | 17.0 | True | True |
| 28 GHz | 0.5-1 | 0.72 | 0.795 | 0.946 | 0.152 | 0.072 | 0.052 | -0.021 | 17.0 | False | True |
| 28 GHz | 1-2 | 1.434 | 1.253 | 1.389 | 0.136 | 0.143 | 0.133 | -0.01 | 17.0 | True | True |
| 28 GHz | 2-5 | 3.246 | 1.856 | 2.022 | 0.166 | 0.325 | 0.295 | -0.03 | 17.0 | False | True |
| 28 GHz | 5-10 | 6.17 | 2.266 | 2.312 | 0.047 | 0.617 | 0.613 | -0.004 | 17.0 | True | True |
| 28 GHz | 10-20 | 11.602 | 2.559 | 2.635 | 0.076 | 1.16 | 1.14 | -0.02 | 12.0 | True | True |
| 38 GHz | 0.1-0.2 | 0.14 | 0.197 | 0.296 | 0.099 | 0.028 | -0.001 | -0.029 | 24.0 | True | True |
| 38 GHz | 0.2-0.5 | 0.33 | 0.425 | 0.513 | 0.089 | 0.066 | 0.05 | -0.016 | 24.0 | True | True |
| 38 GHz | 0.5-1 | 0.721 | 0.795 | 0.926 | 0.131 | 0.144 | 0.122 | -0.022 | 24.0 | True | True |
| 38 GHz | 1-2 | 1.434 | 1.253 | 1.35 | 0.097 | 0.287 | 0.275 | -0.012 | 24.0 | True | True |
| 38 GHz | 2-5 | 3.171 | 1.84 | 1.953 | 0.113 | 0.634 | 0.616 | -0.019 | 24.0 | True | True |
| 38 GHz | 5-10 | 6.17 | 2.266 | 2.24 | -0.026 | 1.234 | 1.267 | 0.033 | 24.0 | True | True |
| 38 GHz | 10-20 | 11.615 | 2.559 | 2.724 | 0.164 | 2.323 | 2.227 | -0.096 | 6.0 | False | False |

### Rule R2: 30 min within a factor 1.5 — 12 cells judged, 0 outside: PASS

| band | class | R_cell | F_true | I | I_err | P_true | s | s_err_abs | links_short | I_ok | s_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.134 | 0.188 | 0.289 | 0.101 | 0.013 | 0.015 | 0.001 | 7.0 | True | True |
| 28 GHz | 0.2-0.5 | 0.365 | 0.463 | 0.564 | 0.101 | 0.036 | 0.008 | -0.029 | 17.0 | True | True |
| 28 GHz | 0.5-1 | 0.709 | 0.785 | 0.811 | 0.026 | 0.071 | 0.067 | -0.004 | 17.0 | True | True |
| 28 GHz | 1-2 | 1.537 | 1.303 | 1.367 | 0.063 | 0.154 | 0.148 | -0.006 | 16.0 | True | True |
| 28 GHz | 2-5 | 3.204 | 1.847 | 1.976 | 0.129 | 0.32 | 0.291 | -0.029 | 17.0 | True | True |
| 28 GHz | 5-10 | 6.307 | 2.278 | 2.337 | 0.06 | 0.631 | 0.621 | -0.01 | 14.0 | True | True |
| 38 GHz | 0.1-0.2 | 0.13 | 0.183 | 0.166 | -0.017 | 0.026 | 0.046 | 0.02 | 13.0 | True | True |
| 38 GHz | 0.2-0.5 | 0.404 | 0.505 | 0.582 | 0.078 | 0.081 | 0.038 | -0.043 | 24.0 | True | True |
| 38 GHz | 0.5-1 | 0.664 | 0.748 | 0.734 | -0.014 | 0.133 | 0.198 | 0.065 | 21.0 | True | True |
| 38 GHz | 1-2 | 1.586 | 1.327 | 1.341 | 0.014 | 0.317 | 0.286 | -0.031 | 21.0 | True | True |
| 38 GHz | 2-5 | 2.943 | 1.786 | 1.894 | 0.108 | 0.589 | 0.582 | -0.006 | 24.0 | True | True |
| 38 GHz | 5-10 | 6.317 | 2.279 | 2.369 | 0.091 | 1.263 | 1.228 | -0.035 | 23.0 | True | True |

**Adopted rule: R2.** (`research/waa_ksurface/outputs/tables/e3_bin_rule_fake.csv`, `research/waa_ksurface/outputs/tables/e3_bin_rule_sizes.csv`)

## The real data under rule R2 (events with g_e = 1 excluded)

Events with g_e = 1 do not enter the length test: E1 measured the radar at about 2.5 times the gauges in them. Left: 14 events, 4563 link-events, 187 links (41 below 1 km). Eligible cells: 12 of 14; short links per eligible cell: min 7, median 17; events per eligible cell: min 7, median 9. (`research/waa_ksurface/outputs/tables/e3_bin_rule_adopted_real_cells.csv`)

| band | class | R_cell | link_events | links | events | links_short | links_long | eligible |
|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.134 | 219 | 66 | 9 | 7 | 59 | True |
| 28 GHz | 0.2-0.5 | 0.365 | 332 | 100 | 11 | 17 | 83 | True |
| 28 GHz | 0.5-1 | 0.709 | 319 | 98 | 11 | 17 | 81 | True |
| 28 GHz | 1-2 | 1.537 | 483 | 108 | 9 | 16 | 92 | True |
| 28 GHz | 2-5 | 3.204 | 565 | 108 | 9 | 17 | 91 | True |
| 28 GHz | 5-10 | 6.307 | 301 | 100 | 7 | 14 | 86 | True |
| 28 GHz | 10-20 | 11.297 | 31 | 16 | 4 | 0 | 16 | False |
| 38 GHz | 0.1-0.2 | 0.13 | 103 | 41 | 8 | 13 | 28 | True |
| 38 GHz | 0.2-0.5 | 0.404 | 164 | 65 | 9 | 24 | 41 | True |
| 38 GHz | 0.5-1 | 0.664 | 174 | 62 | 9 | 21 | 41 | True |
| 38 GHz | 1-2 | 1.586 | 282 | 71 | 9 | 21 | 50 | True |
| 38 GHz | 2-5 | 2.943 | 325 | 71 | 8 | 24 | 47 | True |
| 38 GHz | 5-10 | 6.317 | 216 | 69 | 8 | 23 | 46 | True |
| 38 GHz | 10-20 | 10.655 | 14 | 8 | 3 | 0 | 8 | False |

Not eligible: 28 GHz 10-20 mm/h (4 events, 0 short links); 38 GHz 10-20 mm/h (3 events, 0 short links).

