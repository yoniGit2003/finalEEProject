# A3: F(R) and P(R) from the length test (TRAIN, 28 and 38 GHz)

Written by `research/waa_ksurface/a3_length_test.py`. 8836 link-events, 187 links, 51 events. Rain = radar / g_e. 32 GHz is excluded from the branch (rule A1e). VAL and TEST were not opened.

## (a) The cells

Line (A − c) = I + s·L per band and rain class, in steady rain (Part 1 rule), 200 event-bootstrap draws, 90 % intervals. Eligible = fitted and at least 8 links below 1 km. (`research/waa_ksurface/outputs/tables/a3_length_test_cells.csv`)

| band | class | R_cell (mm/h) | intercept I (dB) | slope s (dB/km) | link_events | links | events | links_under_1km | eligible |
|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.25 | 0.233 | 0.567 [0.357, 0.806] | 0.092 [0.004, 0.248] | 175 | 65 | 12 | 12 | True |
| 28 GHz | 0.25-0.5 | 0.399 | 0.730 [0.623, 0.834] | 0.095 [0.066, 0.141] | 1341 | 113 | 19 | 17 | True |
| 28 GHz | 0.5-1 | 0.752 | 1.010 [0.839, 1.232] | 0.104 [0.065, 0.143] | 2072 | 113 | 20 | 17 | True |
| 28 GHz | 1-2 | 1.445 | 1.370 [1.108, 1.553] | 0.185 [0.143, 0.228] | 2047 | 113 | 19 | 17 | True |
| 28 GHz | 2-5 | 3.215 | 1.965 [1.613, 2.212] | 0.339 [0.278, 0.398] | 1881 | 113 | 18 | 17 | True |
| 28 GHz | 5-10 | 6.231 | 2.548 [2.291, 2.829] | 0.741 [0.581, 0.861] | 1513 | 113 | 15 | 17 | True |
| 28 GHz | 10-20 | 11.874 | 3.086 [2.123, 3.851] | 1.110 [0.583, 1.595] | 499 | 106 | 9 | 17 | True |
| 38 GHz | 0.1-0.25 | 0.222 | 0.594 [0.221, 0.979] | 0.063 [-0.166, 0.327] | 70 | 33 | 11 | 9 | True |
| 38 GHz | 0.25-0.5 | 0.4 | 0.798 [0.663, 0.943] | 0.188 [0.113, 0.290] | 827 | 74 | 18 | 24 | True |
| 38 GHz | 0.5-1 | 0.747 | 1.013 [0.818, 1.237] | 0.237 [0.177, 0.318] | 1126 | 74 | 20 | 24 | True |
| 38 GHz | 1-2 | 1.437 | 1.601 [1.285, 1.812] | 0.317 [0.247, 0.384] | 1203 | 74 | 20 | 24 | True |
| 38 GHz | 2-5 | 3.134 | 2.156 [1.842, 2.380] | 0.701 [0.547, 0.847] | 1159 | 74 | 17 | 24 | True |
| 38 GHz | 5-10 | 6.17 | 2.313 [1.998, 2.726] | 1.565 [1.201, 1.808] | 920 | 74 | 13 | 24 | True |
| 38 GHz | 10-20 | 12.027 | 3.454 [2.079, 4.379] | 2.079 [0.612, 2.773] | 282 | 68 | 9 | 23 | True |

The two added classes (0.1-0.25, 0.25-0.5 mm/h): 28 GHz 0.1-0.25: eligible; 28 GHz 0.25-0.5: eligible; 38 GHz 0.1-0.25: eligible; 38 GHz 0.25-0.5: eligible. A bin of these classes needs the 20 min before it to be wet (R >= 0.3 mm/h), so their values come from rain that is falling off; their intercept can sit above the steady level of that rain.

## (b) Is the WAA steady? The intercept against the time in the rain class

| band | class | I_5 | I_15 | I_30 | I_60 | short_links_5 | short_links_15 | short_links_30 | short_links_60 | diff | diff_lo | diff_hi | excludes_zero |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.25 | 0.567 |  |  |  | 12 | 0 | 0 | 0 |  |  |  | False |
| 28 GHz | 0.25-0.5 | 0.73 | 0.729 | 0.93 |  | 17 | 17 | 9 | 0 | 0.2 | -0.029 | 0.748 | False |
| 28 GHz | 0.5-1 | 1.01 | 0.902 | 0.848 | -0.061 | 17 | 17 | 10 | 0 | -0.162 | -0.415 | 0.389 | False |
| 28 GHz | 1-2 | 1.37 | 1.332 | 1.292 | 2.107 | 17 | 17 | 17 | 0 | -0.078 | -0.336 | 0.369 | False |
| 28 GHz | 2-5 | 1.965 | 1.973 | 1.832 | 2.347 | 17 | 17 | 17 | 0 | -0.133 | -0.339 | 0.129 | False |
| 28 GHz | 5-10 | 2.548 | 2.789 | 2.597 | 2.255 | 17 | 17 | 14 | 1 | 0.049 | -0.522 | 0.417 | False |
| 28 GHz | 10-20 | 3.086 | 5.353 | 3.241 |  | 17 | 3 | 0 | 0 | 0.155 | -7.999 | 8.882 | False |
| 38 GHz | 0.1-0.25 | 0.594 |  |  |  | 9 | 0 | 0 | 0 |  |  |  | False |
| 38 GHz | 0.25-0.5 | 0.798 | 0.758 | 0.832 |  | 24 | 24 | 7 | 0 | 0.034 | -1.58 | 0.232 | False |
| 38 GHz | 0.5-1 | 1.013 | 0.898 | 1.102 |  | 24 | 24 | 8 | 0 | 0.089 | -0.509 | 0.897 | False |
| 38 GHz | 1-2 | 1.601 | 1.506 | 1.88 |  | 24 | 24 | 20 | 0 | 0.279 | -0.379 | 0.696 | False |
| 38 GHz | 2-5 | 2.156 | 2.208 | 1.808 | 3.026 | 24 | 24 | 24 | 3 | -0.348 | -0.551 | 0.079 | False |
| 38 GHz | 5-10 | 2.313 | 2.948 | 2.776 | 3.329 | 24 | 24 | 23 | 4 | 0.462 | -0.305 | 1.172 | False |
| 38 GHz | 10-20 | 3.454 | 4.305 |  |  | 23 | 13 | 0 | 0 |  |  |  | False |

`diff` = I(>= 30 min) − I(>= 5 min), paired over the same event draws, 90 % interval. It excludes zero in 0 of 14 eligible cells (share 0.00; the rule fires above 0.25). **The rule does not fire: F is taken from all steady values.** (`research/waa_ksurface/outputs/tables/a3_length_test_duration.csv`, `research/waa_ksurface/outputs/figures/a3_length_test_duration.png`)

## (c) The rain reference: events with a gauge ratio against events with g_e = 1 (report only)

| band | class | I_gauge | I_lo_gauge | I_hi_gauge | I_ge1 | I_lo_ge1 | I_hi_ge1 | I_gauge_minus_ge1 | s_gauge | s_ge1 | short_links_gauge | short_links_ge1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.25 | 0.661 | 0.528 | 0.878 | 0.201 | 0.201 | 0.201 | 0.459 | 0.083 | -0.017 | 9 | 3 |
| 28 GHz | 0.25-0.5 | 0.766 | 0.673 | 0.874 | 0.267 | -0.316 | 0.338 | 0.499 | 0.107 | 0.038 | 17 | 9 |
| 28 GHz | 0.5-1 | 1.042 | 0.87 | 1.256 | 0.492 | 0.264 | 0.683 | 0.55 | 0.121 | 0.028 | 17 | 10 |
| 28 GHz | 1-2 | 1.403 | 1.146 | 1.59 | 0.478 | 0.059 | 1.033 | 0.926 | 0.194 | 0.128 | 17 | 5 |
| 28 GHz | 2-5 | 2.036 | 1.716 | 2.237 | 0.364 | -2.607 | 0.743 | 1.673 | 0.336 | 0.439 | 17 | 8 |
| 28 GHz | 5-10 | 2.598 | 2.362 | 2.873 | 0.906 | -2.22 | 1.074 | 1.692 | 0.74 | 0.559 | 17 | 5 |
| 28 GHz | 10-20 | 3.086 | 2.123 | 3.851 |  |  |  |  | 1.11 |  | 17 | 0 |
| 38 GHz | 0.1-0.25 | 0.607 | 0.256 | 0.994 | 0.157 | 0.073 | 0.279 | 0.45 | 0.132 | -0.004 | 9 | 1 |
| 38 GHz | 0.25-0.5 | 0.865 | 0.733 | 1.009 | 0.307 | -2.166 | 0.71 | 0.558 | 0.203 | 0.078 | 24 | 17 |
| 38 GHz | 0.5-1 | 1.05 | 0.841 | 1.258 | 0.548 | 0.301 | 1.134 | 0.502 | 0.254 | 0.072 | 24 | 15 |
| 38 GHz | 1-2 | 1.632 | 1.339 | 1.853 | 0.326 | -0.364 | 1.253 | 1.305 | 0.347 | 0.271 | 24 | 9 |
| 38 GHz | 2-5 | 2.175 | 1.878 | 2.399 | 0.755 | 0.588 | 3.866 | 1.42 | 0.745 | 0.359 | 24 | 5 |
| 38 GHz | 5-10 | 2.331 | 2.035 | 2.768 | -0.538 | -0.538 | -0.538 | 2.869 | 1.57 | 2.051 | 24 | 2 |
| 38 GHz | 10-20 | 3.454 | 2.079 | 4.379 |  |  |  |  | 2.079 |  | 23 | 0 |

The F of the gauge-ratio events is kept as the column `F_gauge_events` of the curve table (third sensitivity curve of rule A3f). (`research/waa_ksurface/outputs/tables/a3_length_test_rain_input.csv`, `research/waa_ksurface/outputs/figures/a3_length_test_rain_input.png`)

## (d) Errors-in-variables fake check: STOP

Measured radar-vs-gauge scatter on the gauge-input set (5 min bins with both >= 0.5 mm/h, 14949 bins): σ_1 = 0.708 (1.4826·MAD of ln(radar / gauge); plain sd 0.747), lag-1 correlation ρ = 0.403. The fake rain error has the sd σ_1 / √max(1, L / 2 km), is AR(1) with ρ, and has mean one (E[exp(e)] = 1). Truth: F = A_band·(R / 5)^0.542, P = c_band·R (Part 1 closed forms). 20 replicates; the numbers are the medians. Pass: |I − F_true| <= 0.15 dB and |s − P_true| / P_true <= 20 % in every eligible cell.

| band | class | R_cell | F_true | I | I_err | P_true | s | s_err_rel | I_ok | s_ok |
|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.25 | 0.23 | 0.484 | 0.675 | 0.192 | 0.024 | 0.013 | -0.475 | False | False |
| 28 GHz | 0.25-0.5 | 0.407 | 0.66 | 0.905 | 0.244 | 0.043 | 0.017 | -0.597 | False | False |
| 28 GHz | 0.5-1 | 0.759 | 0.925 | 1.211 | 0.286 | 0.08 | 0.05 | -0.372 | False | False |
| 28 GHz | 1-2 | 1.446 | 1.312 | 1.622 | 0.311 | 0.152 | 0.114 | -0.247 | False | False |
| 28 GHz | 2-5 | 3.021 | 1.956 | 2.165 | 0.201 | 0.317 | 0.305 | -0.044 | False | True |
| 28 GHz | 5-10 | 6.456 | 2.952 | 2.602 | -0.354 | 0.678 | 0.602 | -0.118 | False | True |
| 28 GHz | 10-20 | 12.238 | 4.175 | 2.873 | -1.3 | 1.285 | 1.032 | -0.187 | False | True |
| 38 GHz | 0.1-0.25 | 0.23 | 0.522 | 0.75 | 0.23 | 0.045 | 0.045 | 0.017 | False | True |
| 38 GHz | 0.25-0.5 | 0.408 | 0.713 | 0.981 | 0.268 | 0.079 | 0.07 | -0.114 | False | True |
| 38 GHz | 0.5-1 | 0.762 | 1.0 | 1.331 | 0.328 | 0.148 | 0.138 | -0.062 | False | True |
| 38 GHz | 1-2 | 1.447 | 1.414 | 1.778 | 0.362 | 0.281 | 0.269 | -0.041 | False | True |
| 38 GHz | 2-5 | 3.037 | 2.114 | 2.333 | 0.219 | 0.589 | 0.588 | 0.003 | False | True |
| 38 GHz | 5-10 | 6.478 | 3.188 | 2.803 | -0.376 | 1.257 | 1.065 | -0.155 | False | True |
| 38 GHz | 10-20 | 12.305 | 4.513 | 3.142 | -1.368 | 2.387 | 1.845 | -0.226 | False | False |

Largest intercept error 1.368 dB; largest slope error 59.7 %. (`research/waa_ksurface/outputs/tables/a3_length_test_eiv.csv`, `research/waa_ksurface/outputs/figures/a3_length_test_eiv.png`)

## (e) The curves

Per band, a monotone cubic (PCHIP) through (0, 0) and the eligible cell values; intervals from the bootstrap draws. Product: `research/waa_ksurface/outputs/tables/a3_length_test_curves.csv` (read with `research.waa_ksurface.lt.Curves`). The cell values (`research/waa_ksurface/outputs/tables/a3_length_test_knots.csv`):

| band | quantity | class | R_cell | cell_value | curve_value | cell_lo | cell_hi | changed_by_monotone_step |
|---|---|---|---|---|---|---|---|---|
| 28 GHz | F | 0.1-0.25 | 0.2332 | 0.5665 | 0.5665 | 0.3565 | 0.8063 | False |
| 28 GHz | F | 0.25-0.5 | 0.3988 | 0.7302 | 0.7302 | 0.623 | 0.8338 | False |
| 28 GHz | F | 0.5-1 | 0.7522 | 1.0097 | 1.0097 | 0.8387 | 1.2317 | False |
| 28 GHz | F | 1-2 | 1.4449 | 1.3703 | 1.3703 | 1.108 | 1.5534 | False |
| 28 GHz | F | 2-5 | 3.2152 | 1.9651 | 1.9651 | 1.6135 | 2.2117 | False |
| 28 GHz | F | 5-10 | 6.2311 | 2.548 | 2.548 | 2.2909 | 2.8292 | False |
| 28 GHz | F | 10-20 | 11.8744 | 3.0856 | 3.0856 | 2.123 | 3.8512 | False |
| 28 GHz | P | 0.1-0.25 | 0.2332 | 0.0925 | 0.0925 | 0.0036 | 0.2482 | False |
| 28 GHz | P | 0.25-0.5 | 0.3988 | 0.0954 | 0.0954 | 0.0663 | 0.1406 | False |
| 28 GHz | P | 0.5-1 | 0.7522 | 0.1043 | 0.1043 | 0.0653 | 0.1434 | False |
| 28 GHz | P | 1-2 | 1.4449 | 0.1851 | 0.1851 | 0.1433 | 0.2284 | False |
| 28 GHz | P | 2-5 | 3.2152 | 0.3387 | 0.3387 | 0.2784 | 0.3981 | False |
| 28 GHz | P | 5-10 | 6.2311 | 0.7411 | 0.7411 | 0.5809 | 0.8613 | False |
| 28 GHz | P | 10-20 | 11.8744 | 1.1099 | 1.1099 | 0.5831 | 1.5952 | False |
| 38 GHz | F | 0.1-0.25 | 0.2218 | 0.5942 | 0.5942 | 0.2212 | 0.9792 | False |
| 38 GHz | F | 0.25-0.5 | 0.4004 | 0.7982 | 0.7982 | 0.6634 | 0.9435 | False |
| 38 GHz | F | 0.5-1 | 0.7467 | 1.0132 | 1.0132 | 0.8176 | 1.2371 | False |
| 38 GHz | F | 1-2 | 1.4374 | 1.6007 | 1.6007 | 1.2849 | 1.8123 | False |
| 38 GHz | F | 2-5 | 3.1338 | 2.1565 | 2.1565 | 1.8419 | 2.38 | False |
| 38 GHz | F | 5-10 | 6.1695 | 2.3134 | 2.3134 | 1.9984 | 2.7258 | False |
| 38 GHz | F | 10-20 | 12.0274 | 3.4543 | 3.4543 | 2.079 | 4.379 | False |
| 38 GHz | P | 0.1-0.25 | 0.2218 | 0.0635 | 0.0635 | -0.166 | 0.3271 | False |
| 38 GHz | P | 0.25-0.5 | 0.4004 | 0.1884 | 0.1884 | 0.113 | 0.2899 | False |
| 38 GHz | P | 0.5-1 | 0.7467 | 0.237 | 0.237 | 0.1766 | 0.3176 | False |
| 38 GHz | P | 1-2 | 1.4374 | 0.3167 | 0.3167 | 0.2467 | 0.3838 | False |
| 38 GHz | P | 2-5 | 3.1338 | 0.7007 | 0.7007 | 0.5472 | 0.8466 | False |
| 38 GHz | P | 5-10 | 6.1695 | 1.5651 | 1.5651 | 1.201 | 1.8078 | False |
| 38 GHz | P | 10-20 | 12.0274 | 2.0785 | 2.0785 | 0.6121 | 2.7729 | False |

- **28 GHz: valid range 0.24–11.69 mm/h** (lowest to highest eligible R_cell).

- **38 GHz: valid range 0.23–11.95 mm/h** (lowest to highest eligible R_cell).


**Below the lowest eligible class the curve is an interpolation to zero, not a measurement** (lowest measured class: 28 GHz: 0.1-0.25 mm/h at R_cell 0.23, 38 GHz: 0.1-0.25 mm/h at R_cell 0.22). Above the highest cell F is held constant and P continues with its last slope: extrapolation.

The monotone step changed 0 cell value(s).


## STOP

The errors-in-variables check failed its rule. Nothing after A3 is run.

