# T2 / T3: candidate (iv), and the fake test at the conditional rain error

Written by `research/waa_ksurface/t3_rerun.py`. Fake data only (real data is read for the counts of candidate (iv)). Rule R2, n_short = 5, tolerances 0.15 dB and 0.0710 dB/km (Deviation 3). Rain error of (i)-(iii), from T1: per-event scale part sd 0.45; inside an event a slow error (0.9 correlated over 30 min) with the link's own part sd 0.18 (low end) and 0.26 (high end) and a part shared by all links 0.92 times that. 5 replicates; a cell is judged if eligible in at least 3. VAL and TEST were not opened.

## T2: candidate (iv) on real data

1983 TRAIN link-events of 41 links with a City gauge within 1 km (18 of them below 1 km), 43 events (all TRAIN events). (`research/waa_ksurface/outputs/tables/t3_rerun_gauge_real_cells.csv`)

| band | class | R_cell | link_events | links | events | links_short | links_long | fitted | eligible |
|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.171 | 517 | 13 | 31 | 4 | 9 | True | False |
| 28 GHz | 0.2-0.5 | 0.343 | 505 | 18 | 22 | 6 | 12 | True | True |
| 28 GHz | 0.5-1 | 0.686 | 470 | 18 | 16 | 6 | 12 | True | True |
| 28 GHz | 1-2 | 1.371 | 375 | 18 | 16 | 6 | 12 | True | True |
| 28 GHz | 2-5 | 2.914 | 265 | 18 | 15 | 6 | 12 | True | True |
| 28 GHz | 5-10 | 5.829 | 130 | 17 | 9 | 6 | 11 | True | True |
| 28 GHz | 10-20 | 18.514 | 12 | 7 | 1 | 5 | 2 | False | False |
| 38 GHz | 0.1-0.2 | 0.171 | 435 | 12 | 31 | 6 | 6 | True | True |
| 38 GHz | 0.2-0.5 | 0.343 | 607 | 23 | 21 | 12 | 11 | True | True |
| 38 GHz | 0.5-1 | 0.686 | 555 | 23 | 16 | 12 | 11 | True | True |
| 38 GHz | 1-2 | 1.371 | 456 | 23 | 16 | 12 | 11 | True | True |
| 38 GHz | 2-5 | 3.086 | 323 | 23 | 13 | 12 | 11 | True | True |
| 38 GHz | 5-10 | 5.829 | 184 | 23 | 9 | 12 | 11 | True | True |
| 38 GHz | 10-20 | 12.0 | 20 | 13 | 2 | 7 | 6 | False | False |

Cells that each floor of short links (and as many longer links) would leave. Report only: the floor is not lowered.

| floor_short_links | cells |
|---|---|
| 5 | 11 |
| 4 | 12 |
| 3 | 12 |
| 2 | 12 |
| 1 | 12 |

## T3: verdict per candidate

| candidate | cells_judged_at_both_ends | valid_cells | valid_below_2mmh | largest_I_err_db | qualifies | adopted |
|---|---|---|---|---|---|---|
| (i) gauge-ratio events, plain | 12 | 4 | 2 | 0.422 | False | False |
| (ii) instrument binning (area-mean radar) | 10 | 8 | 8 | 0.086 | True | True |
| (iii) forward-bias correction | 12 | 5 | 5 | 0.537 | False | False |
| (iv) gauge-binned | 12 | 3 | 3 | 1.677 | False | False |

Adoption rule: the most valid cells, with at least 8 valid and at least 3 of them below 2 mm/h. **Adopted: (ii) instrument binning (area-mean radar).**

## Every cell

### (i) gauge-ratio events, plain

| band | class | I_err_lo | s_err_lo | inside_lo | R_cell_lo | I_err_hi | s_err_hi | inside_hi | R_cell_hi | valid | I_err_normal | s_err_normal | inside_normal | R_cell_normal | I_err_doubled | s_err_doubled | inside_doubled | R_cell_doubled |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.132 | -0.008 | True | 0.14 | 0.131 | -0.014 | True | 0.14 | True |  |  |  |  |  |  |  |  |
| 28 GHz | 0.2-0.5 | 0.247 | -0.012 | False | 0.338 | 0.257 | -0.017 | False | 0.343 | False |  |  |  |  |  |  |  |  |
| 28 GHz | 0.5-1 | 0.422 | -0.005 | False | 0.75 | 0.332 | -0.002 | False | 0.748 | False |  |  |  |  |  |  |  |  |
| 28 GHz | 1-2 | 0.285 | 0.004 | False | 1.405 | 0.258 | -0.007 | False | 1.422 | False |  |  |  |  |  |  |  |  |
| 28 GHz | 2-5 | 0.126 | 0.032 | True | 3.164 | 0.107 | 0.043 | True | 3.194 | True |  |  |  |  |  |  |  |  |
| 28 GHz | 5-10 | -0.013 | 0.024 | True | 6.202 | -0.044 | 0.004 | True | 6.187 | True |  |  |  |  |  |  |  |  |
| 38 GHz | 0.1-0.2 | -0.005 | 0.041 | True | 0.132 | 0.088 | -0.004 | True | 0.136 | True |  |  |  |  |  |  |  |  |
| 38 GHz | 0.2-0.5 | 0.258 | 0.001 | False | 0.325 | 0.265 | -0.008 | False | 0.343 | False |  |  |  |  |  |  |  |  |
| 38 GHz | 0.5-1 | 0.29 | 0.016 | False | 0.764 | 0.266 | 0.038 | False | 0.751 | False |  |  |  |  |  |  |  |  |
| 38 GHz | 1-2 | 0.139 | 0.036 | True | 1.383 | 0.083 | 0.074 | False | 1.418 | False |  |  |  |  |  |  |  |  |
| 38 GHz | 2-5 | 0.197 | 0.01 | False | 3.396 | 0.108 | 0.049 | True | 3.319 | False |  |  |  |  |  |  |  |  |
| 38 GHz | 5-10 | 0.082 | -0.028 | True | 6.362 | 0.017 | -0.132 | False | 6.411 | False |  |  |  |  |  |  |  |  |

### (ii) instrument binning (area-mean radar)

| band | class | I_err_lo | s_err_lo | inside_lo | R_cell_lo | I_err_hi | s_err_hi | inside_hi | R_cell_hi | valid | I_err_normal | s_err_normal | inside_normal | R_cell_normal | I_err_doubled | s_err_doubled | inside_doubled | R_cell_doubled |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.006 | -0.001 | True | 0.154 | 0.004 | 0.003 | True | 0.157 | True |  |  |  |  |  |  |  |  |
| 28 GHz | 0.2-0.5 | -0.002 | 0.005 | True | 0.379 | 0.036 | -0.005 | True | 0.379 | True |  |  |  |  |  |  |  |  |
| 28 GHz | 0.5-1 | -0.083 | -0.026 | True | 0.73 | -0.01 | -0.026 | True | 0.735 | True |  |  |  |  |  |  |  |  |
| 28 GHz | 1-2 | 0.063 | -0.044 | True | 1.447 | 0.086 | -0.046 | True | 1.449 | True |  |  |  |  |  |  |  |  |
| 28 GHz | 2-5 | 0.079 | -0.071 | False | 3.062 | 0.026 | -0.053 | True | 3.049 | False |  |  |  |  |  |  |  |  |
| 38 GHz | 0.1-0.2 | 0.012 | 0.0 | True | 0.121 | 0.023 | -0.001 | True | 0.12 | True |  |  |  |  |  |  |  |  |
| 38 GHz | 0.2-0.5 | 0.021 | 0.004 | True | 0.342 | 0.042 | -0.005 | True | 0.348 | True |  |  |  |  |  |  |  |  |
| 38 GHz | 0.5-1 | -0.072 | -0.041 | True | 0.677 | -0.053 | -0.031 | True | 0.663 | True |  |  |  |  |  |  |  |  |
| 38 GHz | 1-2 | 0.043 | -0.052 | True | 1.403 | 0.073 | -0.065 | True | 1.443 | True |  |  |  |  |  |  |  |  |
| 38 GHz | 2-5 | 0.065 | -0.138 | False | 3.0 | -0.029 | -0.111 | False | 3.069 | False |  |  |  |  |  |  |  |  |

### (iii) forward-bias correction

| band | class | I_err_lo | s_err_lo | inside_lo | R_cell_lo | I_err_hi | s_err_hi | inside_hi | R_cell_hi | valid | I_err_normal | s_err_normal | inside_normal | R_cell_normal | I_err_doubled | s_err_doubled | inside_doubled | R_cell_doubled |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.095 | -0.011 | True | 0.14 | 0.038 | -0.013 | True | 0.14 | True |  |  |  |  |  |  |  |  |
| 28 GHz | 0.2-0.5 | -0.022 | -0.011 | True | 0.338 | 0.149 | -0.039 | True | 0.343 | True |  |  |  |  |  |  |  |  |
| 28 GHz | 0.5-1 | 0.537 | -0.045 | False | 0.75 | 0.216 | -0.023 | False | 0.748 | False |  |  |  |  |  |  |  |  |
| 28 GHz | 1-2 | -0.018 | -0.014 | True | 1.405 | 0.024 | -0.052 | True | 1.422 | True |  |  |  |  |  |  |  |  |
| 28 GHz | 2-5 | 0.068 | 0.06 | True | 3.164 | 0.038 | 0.083 | False | 3.194 | False |  |  |  |  |  |  |  |  |
| 28 GHz | 5-10 | 0.041 | 0.111 | False | 6.202 | -0.029 | 0.103 | False | 6.187 | False |  |  |  |  |  |  |  |  |
| 38 GHz | 0.1-0.2 | 0.003 | 0.076 | False | 0.132 | 0.007 | -0.004 | True | 0.136 | False |  |  |  |  |  |  |  |  |
| 38 GHz | 0.2-0.5 | -0.061 | -0.012 | True | 0.325 | -0.052 | 0.034 | True | 0.343 | True |  |  |  |  |  |  |  |  |
| 38 GHz | 0.5-1 | 0.03 | -0.057 | True | 0.764 | -0.025 | -0.038 | True | 0.751 | True |  |  |  |  |  |  |  |  |
| 38 GHz | 1-2 | -0.203 | 0.04 | False | 1.383 | -0.235 | 0.081 | False | 1.418 | False |  |  |  |  |  |  |  |  |
| 38 GHz | 2-5 | 0.036 | 0.091 | False | 3.396 | 0.02 | 0.09 | False | 3.319 | False |  |  |  |  |  |  |  |  |
| 38 GHz | 5-10 | 0.458 | 0.07 | False | 6.362 | 0.4 | -0.014 | False | 6.411 | False |  |  |  |  |  |  |  |  |

### (iv) gauge-binned

| band | class | I_err_lo | s_err_lo | inside_lo | R_cell_lo | I_err_hi | s_err_hi | inside_hi | R_cell_hi | valid | I_err_normal | s_err_normal | inside_normal | R_cell_normal | I_err_doubled | s_err_doubled | inside_doubled | R_cell_doubled |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.2-0.5 |  |  |  |  |  |  |  |  | True | 0.052 | -0.043 | True | 0.343 | 0.061 | -0.024 | True | 0.343 |
| 28 GHz | 0.5-1 |  |  |  |  |  |  |  |  | True | 0.111 | -0.038 | True | 0.686 | 0.109 | -0.034 | True | 0.686 |
| 28 GHz | 1-2 |  |  |  |  |  |  |  |  | False | 0.293 | -0.068 | False | 1.2 | 0.295 | -0.04 | False | 1.2 |
| 28 GHz | 2-5 |  |  |  |  |  |  |  |  | False | 0.179 | -0.037 | False | 2.743 | 0.25 | -0.061 | False | 2.743 |
| 28 GHz | 5-10 |  |  |  |  |  |  |  |  | False | 0.086 | -0.029 | True | 6.0 | 0.077 | -0.109 | False | 6.171 |
| 38 GHz | 0.1-0.2 |  |  |  |  |  |  |  |  | True | -0.024 | -0.065 | True | 0.171 | -0.014 | -0.06 | True | 0.171 |
| 38 GHz | 0.2-0.5 |  |  |  |  |  |  |  |  | False | -0.026 | -0.079 | False | 0.343 | -0.008 | -0.079 | False | 0.343 |
| 38 GHz | 0.5-1 |  |  |  |  |  |  |  |  | False | 0.032 | -0.077 | False | 0.686 | 0.041 | -0.075 | False | 0.686 |
| 38 GHz | 1-2 |  |  |  |  |  |  |  |  | False | 0.164 | -0.022 | False | 1.2 | 0.183 | 0.076 | False | 1.2 |
| 38 GHz | 2-5 |  |  |  |  |  |  |  |  | False | 0.086 | 0.015 | True | 2.743 | 0.115 | 0.135 | False | 2.743 |
| 38 GHz | 5-10 |  |  |  |  |  |  |  |  | False | 0.034 | -0.255 | False | 6.0 | -0.151 | -0.142 | False | 6.171 |
| 38 GHz | 10-20 |  |  |  |  |  |  |  |  | False | 1.677 | -1.911 | False | 11.314 | 0.672 | -1.282 | False | 11.314 |

(`research/waa_ksurface/outputs/tables/t3_rerun_valid_cells.csv`, `research/waa_ksurface/outputs/tables/t3_rerun_fake.csv`, `research/waa_ksurface/outputs/tables/t3_rerun_verdict.csv`, `research/waa_ksurface/outputs/figures/t3_rerun_fake.png`)

