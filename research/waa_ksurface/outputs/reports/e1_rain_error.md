# E1: the radar's own rain error (rain data only, TRAIN events)

Written by `research/waa_ksurface/e1_rain_error.py`. 53 TRAIN events, 6112 five-minute bins. No attenuation data. Errors are multiplicative: sd of ln(rate), on bins where every sensor of the comparison reads >= 0.5 mm/h. 90 % intervals from 200 event draws.

## (a) Triple collocation (15 min)

Error sd of ln(rate) of each sensor, with no sensor taken as the truth. `event_centred = True` is the per-bin part (event means removed).

| triplet | sensors | event_centred | sensor | bins | events | error_sd | sd_lo | sd_hi |
|---|---|---|---|---|---|---|---|---|
| T1 | radar@GbgA, Drakeg, GbgA | False | radar | 144 | 13 | 0.683 | 0.553 | 0.772 |
| T1 | radar@GbgA, Drakeg, GbgA | False | Drakeg | 144 | 13 | 0.211 | 0.0 | 0.367 |
| T1 | radar@GbgA, Drakeg, GbgA | False | GbgA | 144 | 13 | 0.407 | 0.277 | 0.524 |
| T1 | radar@GbgA, Drakeg, GbgA | True | radar | 144 | 13 | 0.561 | 0.456 | 0.645 |
| T1 | radar@GbgA, Drakeg, GbgA | True | Drakeg | 144 | 13 | 0.281 | 0.0 | 0.399 |
| T1 | radar@GbgA, Drakeg, GbgA | True | GbgA | 144 | 13 | 0.375 | 0.187 | 0.493 |
| T2 | radar@GbgA, Lbom, GbgA | False | radar | 153 | 13 | 0.708 | 0.57 | 0.792 |
| T2 | radar@GbgA, Lbom, GbgA | False | Lbom | 153 | 13 | 0.387 | 0.055 | 0.565 |
| T2 | radar@GbgA, Lbom, GbgA | False | GbgA | 153 | 13 | 0.437 | 0.333 | 0.517 |
| T2 | radar@GbgA, Lbom, GbgA | True | radar | 153 | 13 | 0.618 | 0.457 | 0.702 |
| T2 | radar@GbgA, Lbom, GbgA | True | Lbom | 153 | 13 | 0.433 | 0.143 | 0.586 |
| T2 | radar@GbgA, Lbom, GbgA | True | GbgA | 153 | 13 | 0.402 | 0.199 | 0.522 |
| T3 | radar@Drakeg+Lbom, Drakeg, Lbom | False | radar | 154 | 15 | 0.563 | 0.434 | 0.637 |
| T3 | radar@Drakeg+Lbom, Drakeg, Lbom | False | Drakeg | 154 | 15 | 0.251 | 0.0 | 0.348 |
| T3 | radar@Drakeg+Lbom, Drakeg, Lbom | False | Lbom | 154 | 15 | 0.412 | 0.317 | 0.494 |
| T3 | radar@Drakeg+Lbom, Drakeg, Lbom | True | radar | 154 | 15 | 0.451 | 0.367 | 0.502 |
| T3 | radar@Drakeg+Lbom, Drakeg, Lbom | True | Drakeg | 154 | 15 | 0.212 | 0.0 | 0.346 |
| T3 | radar@Drakeg+Lbom, Drakeg, Lbom | True | Lbom | 154 | 15 | 0.436 | 0.339 | 0.513 |
| T4 | radar@Jarn+Askim, Jarn, Askim | False | radar | 127 | 13 | 0.617 | 0.492 | 0.68 |
| T4 | radar@Jarn+Askim, Jarn, Askim | False | Jarn | 127 | 13 | 0.404 | 0.333 | 0.467 |
| T4 | radar@Jarn+Askim, Jarn, Askim | False | Askim | 127 | 13 | 0.22 | 0.0 | 0.322 |
| T4 | radar@Jarn+Askim, Jarn, Askim | True | radar | 127 | 13 | 0.53 | 0.409 | 0.607 |
| T4 | radar@Jarn+Askim, Jarn, Askim | True | Jarn | 127 | 13 | 0.37 | 0.301 | 0.435 |
| T4 | radar@Jarn+Askim, Jarn, Askim | True | Askim | 127 | 13 | 0.278 | 0.113 | 0.358 |

Linear version (error sd in mm/h, bins where any sensor reads >= 0.1 mm/h):

| triplet | sensor | bins | error_sd | sd_lo | sd_hi |
|---|---|---|---|---|---|
| T1 | radar | 541 | 1.926 | 1.296 | 2.278 |
| T1 | Drakeg | 541 | 0.322 | 0.0 | 0.842 |
| T1 | GbgA | 541 | 1.969 | 1.162 | 2.648 |
| T2 | radar | 540 | 2.199 | 1.404 | 2.699 |
| T2 | Lbom | 540 | 1.203 | 0.35 | 1.579 |
| T2 | GbgA | 540 | 1.145 | 0.0 | 1.832 |
| T3 | radar | 540 | 1.423 | 1.015 | 1.71 |
| T3 | Drakeg | 540 | 0.551 | 0.0 | 0.985 |
| T3 | Lbom | 540 | 1.633 | 1.096 | 1.927 |
| T4 | radar | 528 | 1.592 | 0.986 | 1.938 |
| T4 | Jarn | 528 | 0.862 | 0.7 | 0.984 |
| T4 | Askim | 528 | 1.519 | 0.222 | 2.223 |

(`research/waa_ksurface/outputs/tables/e1_rain_error_triple_collocation.csv`)

## (b) Gauge pairs against separation (15 min)

Line through the semivariance of ln(rate) of the 26 City pairs up to 8 km: γ(h) = 0.102 + 0.029·h. **Extrapolated to 0–2 km** (the closest pair is 1.5 km apart): γ(0) = 0.102 (gauge's own error sd 0.32), γ(1 km) = 0.131, γ(2 km) = 0.161. Point-vs-pixel representativeness: sd 0.17 (variance 0.030). (`research/waa_ksurface/outputs/tables/e1_rain_error_gauge_pairs.csv`, `research/waa_ksurface/outputs/figures/e1_rain_error_gauge_pairs.png`)

| gauge_a | gauge_b | separation_km | corr_rate | bins_both_rain | semivar_log |
|---|---|---|---|---|---|
| Jarn | Askim | 1.51 | 0.825 | 190 | 0.175 |
| Drakeg | Lbom | 1.828 | 0.75 | 203 | 0.128 |
| Barl | Lbom | 2.397 | 0.855 | 208 | 0.13 |
| Chalm | Drakeg | 2.641 | 0.815 | 193 | 0.169 |
| Torp | Drakeg | 2.895 | 0.736 | 185 | 0.173 |
| Chalm | Barl | 3.189 | 0.625 | 190 | 0.207 |
| Barl | Drakeg | 3.372 | 0.687 | 181 | 0.148 |
| Chalm | Lbom | 3.603 | 0.692 | 205 | 0.188 |
| Tole | Lbom | 3.764 | 0.653 | 192 | 0.305 |
| Tole | Barl | 3.775 | 0.662 | 183 | 0.198 |
| Torp | Lbom | 4.104 | 0.685 | 195 | 0.272 |
| Torp | Bergsj | 4.262 | 0.738 | 191 | 0.232 |

## (c) The binning reversal (weighing City gauges)

| resolution_min | class_mmh | n_radar_in_class | mean_radar | mean_gauge_given_radar | gauge_over_radar_given_radar | n_gauge_in_class | mean_gauge | mean_radar_given_gauge | radar_over_gauge_given_gauge |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 0.1-0.2 | 2025 | 0.15 | 0.315 | 2.095 | 0 |  |  |  |
| 5 | 0.2-0.5 | 2058 | 0.315 | 0.437 | 1.39 | 0 |  |  |  |
| 5 | 0.5-1 | 1374 | 0.745 | 1.046 | 1.404 | 0 |  |  |  |
| 5 | 1-2 | 1581 | 1.376 | 1.68 | 1.222 | 2120 | 1.2 | 1.03 | 0.858 |
| 5 | 2-5 | 1107 | 3.148 | 3.087 | 0.981 | 1627 | 3.256 | 2.495 | 0.766 |
| 5 | 5-10 | 487 | 6.726 | 5.47 | 0.813 | 544 | 6.889 | 4.11 | 0.597 |
| 5 | 10-20 | 102 | 14.288 | 8.024 | 0.562 | 157 | 13.827 | 4.91 | 0.355 |
| 15 | 0.1-0.2 | 706 | 0.145 | 0.336 | 2.309 | 0 |  |  |  |
| 15 | 0.2-0.5 | 739 | 0.324 | 0.424 | 1.307 | 857 | 0.4 | 0.536 | 1.339 |
| 15 | 0.5-1 | 567 | 0.731 | 0.927 | 1.269 | 321 | 0.8 | 0.943 | 1.178 |
| 15 | 1-2 | 546 | 1.377 | 1.804 | 1.31 | 311 | 1.376 | 1.3 | 0.945 |
| 15 | 2-5 | 418 | 3.16 | 3.011 | 0.953 | 559 | 3.161 | 2.277 | 0.72 |
| 15 | 5-10 | 140 | 6.831 | 6.177 | 0.904 | 202 | 6.762 | 4.354 | 0.644 |
| 15 | 10-20 | 37 | 13.544 | 7.286 | 0.538 | 48 | 13.275 | 8.025 | 0.605 |

The unit-gain split of the event-centred log variances (bins with both >= 0.5 mm/h):

| resolution_min | bins | events | radar_error_var | radar_error_var_lo | radar_error_var_hi | gauge_error_var | ratio_var | shared_var | rho_lag1 |
|---|---|---|---|---|---|---|---|---|---|
| 5 | 2997 | 16 | 0.2972 | 0.231 | 0.3701 | 0.2871 | 0.5843 | 0.2519 | 0.3493 |
| 15 | 1109 | 15 | 0.2152 | 0.1573 | 0.2764 | 0.2934 | 0.5086 | 0.3079 | 0.2016 |

(`research/waa_ksurface/outputs/tables/e1_rain_error_reversal.csv`, `research/waa_ksurface/outputs/tables/e1_rain_error_split.csv`, `research/waa_ksurface/outputs/figures/e1_rain_error_reversal.png`)

## (d) The result

| quantity | value | lo | hi |
|---|---|---|---|
| rho_lag1 of the 5 min ln(radar / gauge), event-centred | 0.3493 |  |  |
| factor var(15 min mean) / var(5 min bin) | 0.5157 |  |  |
| est_TC: per-bin radar error sd at 5 min (from 15 min triple collocation) | 0.7604 | 0.6036 | 0.8718 |
| est_rev: per-bin radar error sd at 5 min (unit-gain split) | 0.5452 | 0.4806 | 0.6084 |
| sigma_lo (per-bin part, 5 min, one pixel) | 0.4806 |  |  |
| sigma_hi (per-bin part, 5 min, one pixel) | 0.8718 |  |  |
| all-to-radar limit: sd of ln(radar / gauge) at 5 min, event-centred | 0.7644 |  |  |
| gauge error sd at 5 min (unit-gain split) | 0.5358 |  |  |
| sigma_s: scale part, sd over events of ln(sum radar / sum gauge) | 0.6965 |  |  |
| median over events of ln(sum radar / sum gauge) | -0.1016 |  |  |
| events in sigma_s | 14.0 |  |  |
| sigma_s_res: scale part left after g_e (gauge's own event ratio minus the pooled one) | 0.4284 |  |  |
| shared-in-space variance of the per-bin part (15 min, mean pair covariance) | 0.1377 |  |  |
| variance of ln(radar / gauge) at 15 min on the same bins | 0.4089 |  |  |
| gauge own error sd at 15 min (pair variogram at h = 0) | 0.3194 |  |  |
| point-vs-pixel representativeness sd at 15 min (pair variogram) | 0.1746 |  |  |
| g_e = 1 events: number | 39.0 |  |  |
| g_e = 1 events: pooled radar mm per gauge | 5.8472 |  |  |
| g_e = 1 events: pooled gauge mm per gauge | 2.29 |  |  |

(`research/waa_ksurface/outputs/tables/e1_rain_error_synthesis.csv`, `research/waa_ksurface/outputs/tables/e1_rain_error_event_scale.csv`, `research/waa_ksurface/outputs/tables/e1_rain_error_shared_in_space.csv`)

