# wet_antenna: the physical consistency of every WAA model in the comparison (TEST)

Written by `research/wet_antenna/comparison_all_models.py`. Settings: `wet_antenna.comparison_all_models` in `research/config.yaml`.

**Status.** A post-freeze summary. TEST was read in S8 / S8b / S8c. The frozen fits of those steps are run forward again; nothing is refitted and our equation is unchanged (fit sha256 `1d8163cd8b84aafc85f59c6c223f44b33923278265e3e40819923d8076c8e0b9`, checked).

- The S8 P1 cells of the table-F candidate are reproduced exactly: **True**. The S8c dry-path errors of our equation are reproduced: **True** (1791 dry-path windows).
- Length test: exactly the S8 P1 procedure. The model-free intercept I and slope s of each cell (Huber line of A - c against L, event bootstrap 90 %) against the same line fitted through the model's own predictions on the same bins. Eligible cell: >= 8 links under 1 km.
- Dry path: median over the windows where the radar shows no rain on any pixel of the path (S7 / S8 rules) of measured W minus model W; positive = the model holds too little water. 90 % event bootstrap.

## The table

| label | cells | s_inside | I_inside | both_inside | s_ratio_model_over_obs_median | I_ratio_model_over_obs_median | dry_err_5min_db | dry_err_10min_db | dry_err_20min_db | a_link_vs_L_rho | a_link_vs_L_lo | a_link_vs_L_hi | wet_time_min | dry_time_min | params_global_fitted |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| no WAA (B0) | 8 | 0 | 0 | 0 | 2.358 | 0.008 | 0.35 | 0.267 | 0.15 |  |  |  |  |  | 3 |
| constant WAA while it rains (B1) | 8 | 4 | 0 | 0 | 1.315 | 0.544 | 0.35 | 0.267 | 0.15 | -0.597 | -0.662 | -0.515 | 0.0 | 0.0 | 6 |
| WAA without memory (B2) | 8 | 7 | 6 | 6 | 0.973 | 0.806 | 0.347 | 0.267 | 0.15 | 0.204 | 0.078 | 0.321 | 0.0 | 0.0 | 12 |
| Schleiss (a) | 8 | 5 | 2 | 1 | 0.649 | 1.384 | -0.052 | 0.103 | 0.125 |  |  |  | 4.481 | 4.481 | 3 |
| Schleiss (b) | 8 | 6 | 7 | 5 | 1.149 | 0.777 | -0.239 | -0.249 | -0.183 | -0.453 | -0.554 | -0.346 | 102.578 | 25.137 | 8 |
| Pastorek (a) | 8 | 4 | 7 | 4 | 0.584 | 1.05 | 0.333 | 0.25 | 0.15 |  |  |  | 0.0 | 0.0 | 3 |
| Pastorek (b) | 8 | 7 | 6 | 6 | 0.968 | 0.824 | 0.333 | 0.25 | 0.15 | 0.308 | 0.184 | 0.42 | 0.0 | 0.0 | 7 |
| ours (final equation) | 8 | 7 | 4 | 4 | 1.045 | 0.716 | -0.185 | -0.146 | -0.107 | 0.087 | -0.034 | 0.203 | 2.6 | 10.3 | 7 |
| ours, table F (S8 candidate) | 8 | 7 | 4 | 4 | 1.051 | 0.734 | -0.184 | -0.147 | -0.11 | 0.047 | -0.077 | 0.172 | 2.6 | 10.3 | free curves (S6 table) |

Measured dry-path WAA (median): 0.35 dB at 5 min, 0.27 dB at 10 min, 0.15 dB at 20 min. Gauge timing of our equation: wetting 2.6 min, drying 10.3 min on TRAIN (VAL 2.0 and 8.8 min).

Timing and parameters:

| label | timing | params_global_fitted | params_a_link | params_identified | params_at_bound |
|---|---|---|---|---|---|
| no WAA (B0) | no antenna term | 3 | in the vector, no effect (no antenna term) | no check: one TRAIN fit, no fake-data recovery, no intervals | not checked |
| constant WAA while it rains (B1) | none: a step with the radar wet flag (V1 > 0.1 mm/h) | 6 | fitted per link (ridge) | no check: one TRAIN fit, no fake-data recovery, no intervals | not checked |
| WAA without memory (B2) | none: follows the rain input | 12 | fitted per link (ridge) | no check: one TRAIN fit, no fake-data recovery, no intervals | not checked |
| Schleiss (a) | nominal tau_wet 15.0 min, tau_dry 15.0 min in the step 3 dt / tau; 1/e times given | 3 | 1 (fixed) | no check: one TRAIN fit, no fake-data recovery, no intervals | none |
| Schleiss (b) | nominal tau_wet 309.2 min, tau_dry 76.9 min in the step 3 dt / tau; 1/e times given | 8 | fitted per link (ridge) | no check: one TRAIN fit, no fake-data recovery, no intervals | none; 0 a_link at a bound |
| Pastorek (a) | none: follows the rain input | 3 | 1 (fixed) | no check: one TRAIN fit, no fake-data recovery, no intervals | none |
| Pastorek (b) | none: follows the rain input | 7 | fitted per link (ridge) | no check: one TRAIN fit, no fake-data recovery, no intervals | none; 0 a_link at a bound |
| ours (final equation) | fast store (share 0.72) plus a slow store (tau_s 240 min, not identified) | 7 | fitted per link (ridge), p10-p90 0.77-1.24 | 9 of 9 numbers identified (S8b); 90 % intervals from 40 event-bootstrap refits, none at a bound | fixed: tau_wet, tau_dry, phi (gauge fit), tau_s at its search bound (not identified) |
| ours, table F (S8 candidate) | as ours | free curves (S6 table) | as ours | S6 fake recovery (62 numbers within 20 %) | as ours |

Dry-path errors with intervals:

| model | 5 | 10 | 20 |
|---|---|---|---|
| B0 | +0.350 [+0.300, +0.387] | +0.267 [+0.232, +0.300] | +0.150 [+0.100, +0.200] |
| B1 | +0.350 [+0.300, +0.387] | +0.267 [+0.232, +0.300] | +0.150 [+0.100, +0.200] |
| B2 | +0.347 [+0.300, +0.383] | +0.267 [+0.217, +0.300] | +0.150 [+0.100, +0.200] |
| Schleiss a | -0.052 [-0.104, +0.004] | +0.103 [+0.053, +0.123] | +0.125 [+0.073, +0.175] |
| Schleiss b | -0.239 [-0.393, -0.113] | -0.249 [-0.399, -0.117] | -0.183 [-0.319, -0.063] |
| Pastorek a | +0.333 [+0.300, +0.350] | +0.250 [+0.200, +0.300] | +0.150 [+0.083, +0.200] |
| Pastorek b | +0.333 [+0.300, +0.351] | +0.250 [+0.200, +0.300] | +0.150 [+0.092, +0.200] |
| ours | -0.185 [-0.242, -0.153] | -0.146 [-0.196, -0.103] | -0.107 [-0.135, -0.066] |
| table F | -0.184 [-0.250, -0.138] | -0.147 [-0.203, -0.098] | -0.110 [-0.129, -0.061] |

Length-test cells (eligible), the model's slope and intercept against the model-free interval:

| band | class | R_cell | s_obs | s_lo | s_hi | I_obs | I_lo | I_hi |
|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.5-1 | 0.761 | 0.037 | 0.014 | 0.077 | 1.354 | 0.994 | 1.517 |
| 28 GHz | 1-2 | 1.392 | 0.164 | 0.086 | 0.234 | 1.433 | 1.145 | 1.705 |
| 28 GHz | 2-5 | 2.69 | 0.313 | 0.205 | 0.439 | 1.794 | 1.216 | 2.129 |
| 28 GHz | 5-10 | 5.923 | 0.603 | 0.513 | 0.863 | 2.603 | 1.114 | 3.466 |
| 38 GHz | 0.5-1 | 0.754 | 0.15 | 0.082 | 0.221 | 1.437 | 1.106 | 1.599 |
| 38 GHz | 1-2 | 1.341 | 0.332 | 0.203 | 0.455 | 1.554 | 1.2 | 1.791 |
| 38 GHz | 2-5 | 2.668 | 0.598 | 0.296 | 0.787 | 2.004 | 1.541 | 2.311 |
| 38 GHz | 5-10 | 5.669 | 0.85 | 0.121 | 1.57 | 3.421 | 1.013 | 4.901 |

| band | class | B0 | B1 | B2 | Schleiss a | Schleiss b | Pastorek a | Pastorek b | ours | table F |
|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.5-1 | 0.210 (out) | 0.046 (in) | 0.098 (out) | 0.046 (in) | 0.042 (in) | 0.037 (in) | 0.103 (out) | 0.096 (out) | 0.095 (out) |
| 28 GHz | 1-2 | 0.374 (out) | 0.179 (in) | 0.143 (in) | 0.082 (out) | 0.188 (in) | 0.067 (out) | 0.153 (in) | 0.168 (in) | 0.167 (in) |
| 28 GHz | 2-5 | 0.708 (out) | 0.453 (out) | 0.273 (in) | 0.154 (out) | 0.396 (in) | 0.115 (out) | 0.266 (in) | 0.295 (in) | 0.296 (in) |
| 28 GHz | 5-10 | 1.537 (out) | 1.083 (out) | 0.588 (in) | 0.335 (out) | 0.984 (out) | 0.270 (out) | 0.552 (in) | 0.643 (in) | 0.653 (in) |
| 38 GHz | 0.5-1 | 0.365 (out) | 0.107 (in) | 0.201 (in) | 0.115 (in) | 0.124 (in) | 0.093 (in) | 0.200 (in) | 0.178 (in) | 0.184 (in) |
| 38 GHz | 1-2 | 0.669 (out) | 0.357 (in) | 0.322 (in) | 0.208 (in) | 0.324 (in) | 0.194 (out) | 0.333 (in) | 0.311 (in) | 0.317 (in) |
| 38 GHz | 2-5 | 1.289 (out) | 0.826 (out) | 0.513 (in) | 0.401 (in) | 0.689 (in) | 0.348 (in) | 0.503 (in) | 0.522 (in) | 0.527 (in) |
| 38 GHz | 5-10 | 2.839 (out) | 2.045 (out) | 1.129 (in) | 0.884 (in) | 1.721 (out) | 0.797 (in) | 1.085 (in) | 1.187 (in) | 1.213 (in) |

| band | class | B0 | B1 | B2 | Schleiss a | Schleiss b | Pastorek a | Pastorek b | ours | table F |
|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.5-1 | 0.00 (out) | 0.90 (out) | 0.88 (out) | 2.30 (out) | 1.14 (in) | 1.16 (in) | 0.87 (out) | 0.80 (out) | 0.83 (out) |
| 28 GHz | 1-2 | 0.02 (out) | 0.90 (out) | 1.20 (in) | 2.30 (out) | 1.11 (out) | 1.59 (in) | 1.22 (in) | 1.05 (out) | 1.08 (out) |
| 28 GHz | 2-5 | 0.07 (out) | 0.90 (out) | 1.82 (in) | 2.32 (out) | 1.38 (in) | 2.27 (out) | 1.71 (in) | 1.50 (in) | 1.59 (in) |
| 28 GHz | 5-10 | 0.13 (out) | 0.99 (out) | 2.44 (in) | 2.33 (in) | 1.46 (in) | 3.29 (in) | 2.47 (in) | 2.11 (in) | 2.12 (in) |
| 38 GHz | 0.5-1 | 0.01 (out) | 0.93 (out) | 0.93 (out) | 2.30 (out) | 1.18 (in) | 1.16 (in) | 0.93 (out) | 0.87 (out) | 0.89 (out) |
| 38 GHz | 1-2 | -0.01 (out) | 0.91 (out) | 1.21 (in) | 2.30 (out) | 1.23 (in) | 1.54 (in) | 1.24 (in) | 1.09 (out) | 1.11 (out) |
| 38 GHz | 2-5 | 0.04 (out) | 0.95 (out) | 2.00 (in) | 2.31 (out) | 1.57 (in) | 2.23 (in) | 1.86 (in) | 1.66 (in) | 1.76 (in) |
| 38 GHz | 5-10 | -0.01 (out) | 0.89 (out) | 2.64 (in) | 2.29 (in) | 1.68 (in) | 3.17 (in) | 2.60 (in) | 2.28 (in) | 2.29 (in) |

## Reading

- Path loss: our model's slope lies inside the measured interval in 7 of 8 cells (median model / measured 1.04); Schleiss (b) 6 (1.15), Pastorek (b) 7 (0.97), no WAA 0 (2.36).
- WAA level at L = 0: ours inside in 4 of 8 (the light-rain cells of S8 P1 fail); Schleiss (b) 7, Pastorek (b) 6, no WAA 0.
- After rain (10 min): ours -0.15 dB, Schleiss (b) -0.25, Pastorek (b) +0.25, no memory (B2) +0.27.
- a_link against L (Spearman, 90 % link bootstrap): constant WAA while it rains (B1) -0.60 [-0.66, -0.52]; WAA without memory (B2) +0.20 [+0.08, +0.32]; Schleiss (b) -0.45 [-0.55, -0.35]; Pastorek (b) +0.31 [+0.18, +0.42]; ours (final equation) +0.09 [-0.03, +0.20]; ours, table F (S8 candidate) +0.05 [-0.08, +0.17].
- These are numbers, not verdicts: no pass rule was written for this summary. P1 (both inside in >= 6 of 8 cells) is shown for each model only as the S8 rule would read it.

## Unclear / possibly wrong

- Post-freeze: the literature fits were made after TEST had been read (S8c); nothing here was tuned on TEST.
- B0-B2 and the literature models have no fake-data recovery and no intervals; their parameter counts are the fitted global numbers (a_link per link counted apart).
- The Schleiss times are the 1/e times of its discrete step (W += (target - W) 3 dt / tau, dt = 1 min), not its nominal tau.