# Step 4b: why do some real tails stall or rise? Contamination or antenna physics?

Script: `research/wet_antenna/diagnose_stalled_tails.py`. Settings and every rule: `research/config.yaml`, section
`wet_antenna.diagnose_stalled_tails` (written before the run; nothing was re-run with changed rules).
Data: the 429 TRAIN tails of the frozen set v1 (24 events, 159 links), first 60 min after
`t0_trace`. The 123 VAL tails were never read; no file of an earlier step in `outputs/products/` was changed.

**The question.** Step 4 found that about one tail in five does not dry in its first minutes and that the drying rate
per dB falls with the level. Either something still wets the antenna after the rain reference says "dry" (H1,
contamination: it belongs to the weather), or the antenna itself dries more slowly per dB when there is more water on
it (H2, antenna physics: it belongs to the drying law). The two explanations predict different things, and this step
looks for those differences.

## 1. What was run, what was cut, run time

- Classes of all 429 tails; part A (outside evidence); part B (B1 to B5); part C (two code fixes, the check on
  6 fake datasets, the menu fits and the free curve again on 4 real sets).
- Timing first: one likelihood evaluation 70 ms (menu law) and 96 ms (free curve);
  one F1 fit with both fixes 2.7 s. Projection 0.58 h on 3 processes, budget
  2.5 h. **Cut: nothing.**
- Fits: 200 menu fits and 4 free-curve fits, 31 min of wall time for the fits; everything else
  takes about a minute. Cached in `research/.cache/wet_antenna/diagnose_stalled_tails/`.
- The spec named files and tools of the old layout (`roadmap.md`, `DECISIONS.md`, `save_output` with importance
  levels, `check_organization.py`, `.cache/fits4b/`). The repository was reorganised since: this step uses
  `docs/ROADMAP.md`, `METHOD.md`, `save`, the cache folder above, and the checks that exist now. There are no
  importance levels any more, so no output is "key", "product" or "diag"; the per-tail table is in `outputs/products/`
  because later steps will read it.

## 2. The classes

Rules (fixed first). A is smoothed with a centred 5 min median; A_start = mean of minutes 0 to 4.
The classes are tried in this order, and a tail gets the first one it meets:

1. **decays**: A falls below 0.5 x A_start within 30 min and never climbs more than 0.3 dB above its lowest level so far.
2. **plateau-then-decays**: A stays at or above 0.8 x A_start for at least 10 min, then falls below 0.5 x A_start before 60 min.
3. **rises**: A exceeds A_start + 0.3 dB at some time.
4. **flat / slow**: everything else.

The order is my reading of the spec (it did not say what happens when two rules hold). Tails meeting two rules:
decays and plateau: 28, plateau and rises: 19, decays and rises: 35. The columns `meets_decays`, `meets_plateau`, `meets_rises` of the per-tail table
keep all three answers.

| group | tails | decays | plateau-then-decays | rises | flat / slow | share_not_decaying |
|---|---|---|---|---|---|---|
| all tails | 429 | 240 | 64 | 36 | 89 | 0.441 |
| other events | 265 | 182 | 24 | 11 | 48 | 0.313 |
| events flagged rise-heavy | 164 | 58 | 40 | 25 | 41 | 0.646 |

Per event (largest first):

| group | rise_heavy | tails | decays | plateau-then-decays | rises | flat / slow | share_not_decaying |
|---|---|---|---|---|---|---|---|
| E20150804_1830 | False | 115 | 87 | 10 | 0 | 18 | 0.243 |
| E20150826_0520 | False | 75 | 61 | 3 | 5 | 6 | 0.187 |
| E20150628_0530 | True | 56 | 28 | 17 | 0 | 11 | 0.5 |
| E20150617_0705 | True | 54 | 9 | 17 | 10 | 18 | 0.833 |
| E20150727_2140 | False | 21 | 10 | 1 | 3 | 7 | 0.524 |
| E20150705_1910 | True | 15 | 6 | 2 | 0 | 7 | 0.6 |
| E20150701_0030 | True | 15 | 0 | 0 | 12 | 3 | 1.0 |
| E20150712_2255 | True | 13 | 9 | 3 | 0 | 1 | 0.308 |
| E20150825_0145 | False | 9 | 4 | 2 | 0 | 3 | 0.556 |
| E20150803_0210 | False | 8 | 1 | 2 | 0 | 5 | 0.875 |
| E20150826_1620 | False | 7 | 0 | 2 | 1 | 4 | 1.0 |
| E20150607_0245 | False | 6 | 5 | 0 | 1 | 0 | 0.167 |
| E20150815_1305 | False | 5 | 4 | 0 | 0 | 1 | 0.2 |
| E20150810_1545 | False | 5 | 3 | 0 | 0 | 2 | 0.4 |
| E20150605_1745 | True | 5 | 2 | 0 | 2 | 1 | 0.6 |
| E20150802_0400 | True | 4 | 4 | 0 | 0 | 0 | 0.0 |
| E20150725_0440 | False | 4 | 3 | 0 | 1 | 0 | 0.25 |
| E20150707_1235 | False | 3 | 1 | 2 | 0 | 0 | 0.667 |
| E20150829_0545 | False | 3 | 3 | 0 | 0 | 0 | 0.0 |
| E20150813_0235 | False | 2 | 0 | 0 | 0 | 2 | 1.0 |
| E20150704_2125 | True | 1 | 0 | 1 | 0 | 0 | 1.0 |
| E20150731_1440 | True | 1 | 0 | 0 | 1 | 0 | 1.0 |
| E20150806_0330 | False | 1 | 0 | 1 | 0 | 0 | 1.0 |
| E20150811_2220 | False | 1 | 0 | 1 | 0 | 0 | 1.0 |

## 3. Part A: outside evidence of continued wetting (H1)

Per tail, from `t0_trace` to `t0_trace` + 60 min. An indicator is "true" when the reference shows any water at all
(no threshold). Odds ratio = odds in the non-decaying tails (classes 2 to 4, 189 tails) over the odds in the
decaying tails (240); the interval is the 90 % event bootstrap (1000 draws of whole events).

| indicator | tails_known | share_not_decaying | share_decays | share_plateau-then-decays | share_rises | share_flat / slow | odds_ratio | or_lo | or_hi | rule_met |
|---|---|---|---|---|---|---|---|---|---|---|
| radar R_max > 0 on the path of the link | 429 | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |  |  |  |  |
| radar R_max > 0 on a link path within 5 km | 429 | 0.73 | 0.879 | 0.703 | 0.778 | 0.73 | 0.38 | 0.21 | 0.66 |  |
| nearest City gauge > 0 mm | 429 | 0.201 | 0.412 | 0.188 | 0.083 | 0.258 | 0.36 | 0.2 | 0.57 |  |
| SMHI gauge > 0 mm | 429 | 0.243 | 0.433 | 0.203 | 0.278 | 0.258 | 0.42 | 0.28 | 0.98 |  |
| ANY of the four rain indicators | 429 | 0.778 | 0.925 | 0.703 | 0.833 | 0.809 | 0.29 | 0.1 | 0.58 |  |
| dew risk (hourly product of the baseline module) | 360 | 0.767 | 0.87 | 0.83 | 0.5 | 0.779 | 0.49 | 0.25 | 1.37 |  |
| ANY moisture evidence (rain indicators or dew risk) = THE RULE | 416 | 0.983 | 0.983 | 0.983 | 1.0 | 0.977 | 0.94 | 0.17 | 4.5 | False |
| night at t0_trace (NOT moisture evidence; information, added after the first look) | 429 | 0.291 | 0.075 | 0.188 | 0.5 | 0.281 | 4.96 | 1.64 | 14.53 |  |

**Read this table with one fact in mind.** The frozen tail set v1 starts every tail at the LAST rain trace
(`t0_trace`) and cuts it at the next one, where "trace" = radar R_max > 0 on the path, a radar gap, or a gauge within
1 km > 0. All 429 tails are longer than 60 min, so inside the window the radar on the path and the near gauges are
dry BY CONSTRUCTION (radar on the path > 0 in 0 tails; 0 tails with a missing radar bin).
The only rain references that CAN show something are the ones farther away: other link paths within 5 km,
the nearest City gauge when it is more than 1 km away (median distance 6.2 km), and the SMHI gauge
(median 10.3 km). Humidity exists for 84 % of the tails: it is missing in the
hourly SMHI record for the events E20150617_0705, E20150701_0030 (69 tails, of which
60 do not decay and 22 of the 36 "rises" tails). There the dew indicator is
"unknown" and those tails are left out of its row. The weather is ONE station for the whole city, hourly.

Weather and time of day per class (medians):

| class | tails | night_share | humidity_pct | temperature_c | dew_point_depression_c | wind_speed_ms | hour_local | sun_elevation_deg | radar_near_wet_share | city_gauge_km |
|---|---|---|---|---|---|---|---|---|---|---|
| decays | 240 | 0.08 | 95.0 | 17.0 | 0.8 | 2.22 | 9.29 | 29.52 | 0.06 | 6.66 |
| plateau-then-decays | 64 | 0.19 | 92.0 | 15.3 | 0.95 | 2.35 | 11.08 | 9.96 | 0.04 | 1.75 |
| rises | 36 | 0.5 | 87.0 | 13.1 | 1.49 | 3.4 | 12.5 | 0.64 | 0.05 | 12.37 |
| flat / slow | 89 | 0.28 | 92.5 | 15.8 | 0.95 | 2.65 | 11.08 | 16.07 | 0.03 | 5.91 |

Other parts of part A:

- **Same class in the event's other tails at the same time** (other links, start within 30 min): the share of
  "rises" or "flat / slow" neighbours is 0.39 for the non-decaying tails and 0.21 for the decaying ones,
  ratio **1.85** (event bootstrap 1.22..2.39; rule: > 2.0). 22 tails
  have no such neighbour.
- **The two directions of one link** (129 link-events with both sub-links): same answer to "decays or not" in
  81 % (expected 64 % if the tails of an event were independent); same class of four in 75 %.
  The two directions share both antennas and the path, so this shows the classes are not noise; it cannot separate H1 from H2.
- **Does the rise grow with the path length?** Rise = highest smoothed A minus A_start. Spearman with L: all tails
  -0.00 (-0.13..+0.13, n = 429); tails of class "rises" +0.08 (-0.12..+0.44, n = 36).

| links | tails | median_rise_db | share_rises | share_not_decaying |
|---|---|---|---|---|
| long, L >= 2 km | 250 | 0.13 | 0.092 | 0.44 |
| short, L < 2 km | 179 | 0.15 | 0.073 | 0.441 |

**H1 verdict by the rule: NOT SUPPORTED.** Odds-ratio part ("any moisture evidence" > 2.0 with the
interval excluding 1): not met. Coherence part (ratio > 2.0): not met.

## 4. Part B: tests of a law of the antenna (H2)

**B1, transition level.** 64 plateau-then-decays tails. Transition level: median 0.96 dB, IQR
**0.35 dB**; as a fraction of A_start: median 0.76, IQR 0.06. For comparison the IQR of
A_start itself is 0.52 dB in these tails and 0.48 dB in all tails. Links with two or more such tails:
15, median spread inside a link 0.20 dB.
Warning: by the definition in the rules the transition level is the first A below 0.8 x A_start, so it is
about 0.8 x A_start for every tail. Its spread is therefore the spread of A_start, and the fraction says nothing.

**B2, the same link in different events.** 79 links with tails in at least 2 events (220 link-events).
Statistic (mean product of "does not decay", event mean removed, over pairs of events of a link): +0.0038; the
1000 shuffles inside events give 0 +- 0.0096; **p = 0.350**.

**B3, plateau duration against wetness** (Spearman, 90 % event bootstrap):

| plateau duration against | plateau-then-decays tails only (the rule) | all tails, minutes until A < 0.8 A_start (60 if never; information) |
|---|---|---|
| A_start (dB) | +0.16 (-0.02..+0.40, n = 64) | +0.00 (-0.08..+0.12, n = 429) |
| rain amount of the event on the path (radar, mm) | -0.09 (-0.50..+0.11, n = 64) | -0.20 (-0.48..+0.11, n = 429) |
| rain duration on the link (min) | -0.14 (-0.50..+0.06, n = 64) | -0.19 (-0.47..+0.08, n = 429) |

**B4, master curve.** Built from the 182 "decays" tails of the events not flagged rise-heavy (182 of them cross
0.8 dB and gave the first alignment; 182 are in the final curve), one stretch per link,
3 rounds. The curve is defined from 2.56 dB down to 0.02 dB over 142 min.

| group | tails | placed | median_rms_db | p25_rms_db | p75_rms_db |
|---|---|---|---|---|---|
| 'decays', events not rise-heavy (they built the curve; own shift and stretch) | 182 | 182 | 0.124 | 0.091 | 0.171 |
| 'decays', rise-heavy events (NOT used to build the curve) | 58 | 58 | 0.126 | 0.088 | 0.156 |
| 'plateau-then-decays', decaying part only | 64 | 64 | 0.117 | 0.082 | 0.18 |

The plateau part of the plateau-then-decays tails, with the shift and stretch found on their decaying part:
64 of 64 tails could be placed, 61 have the master curve defined under their plateau.
Mean (A - master) over the plateau: median -0.27 dB; **on** the curve (within the median misfit
0.124 dB of the first row) 8, **below** it 45, **above** it 8. "Below" means: going back in time
the master curve is higher than the tail was, so the tail sat at its level instead of coming down from above; that
is a delayed start, not a slow part of one common curve.

**B5, the model-free drying rate per level.** k = -mean dA/dt / median A per bin of A (A smoothed over
7 min, central difference), bins with at least 200 tail-minutes and k > 0:

| subset | tails | bins used | slope of log k vs log A | event-bootstrap 90 % |
|---|---|---|---|---|
| decaying tails, events not rise-heavy (the rule) | 182 | 5 | 0.28 | +0.05..+0.77 |
| decaying tails, all events | 240 | 5 | 0.1 |  |
| all tails | 429 | 6 | -0.8 | -1.64..-0.40 |
| events not rise-heavy | 265 | 5 | -0.32 |  |
| rise-heavy events | 164 | 6 | -1.69 |  |
| not rise-heavy, A_start < 1.5 dB | 170 | 4 | -0.54 |  |
| not rise-heavy, A_start >= 1.5 dB | 95 | 5 | -0.12 |  |

k per bin (1 / min):

| A_bin_db | decaying tails, events not rise-heavy (the rule) | decaying tails, all events | all tails | events not rise-heavy | rise-heavy events | not rise-heavy, A_start < 1.5 dB | not rise-heavy, A_start >= 1.5 dB |
|---|---|---|---|---|---|---|---|
| (0.35, 0.5] | 0.0536 | 0.0531 | 0.0522 | 0.0513 | 0.0542 | 0.0487 | 0.0568 |
| (0.5, 0.7] | 0.0587 | 0.0567 | 0.0477 | 0.0535 | 0.0393 | 0.0514 | 0.0571 |
| (0.7, 1.0] | 0.0816 | 0.0727 | 0.0425 | 0.0559 | 0.0295 | 0.0495 | 0.0726 |
| (1.0, 1.4] | 0.0758 | 0.0623 | 0.0282 | 0.0397 | 0.0177 | 0.026 | 0.0746 |
| (1.4, 2.0] | 0.0763 | 0.059 | 0.0147 | 0.0344 | 0.0044 |  | 0.0389 |
| (2.0, 2.83] |  |  | 0.0164 |  | 0.0038 |  |  |

Added after the first look (information, no rule): the same k, split by WHEN in the tail the level is reached. If
only A decided the rate, the two columns of a subset would agree.

| A_bin_db | all tails, minutes 0-14 | all tails, minutes 15-59 | events not rise-heavy, minutes 0-14 | events not rise-heavy, minutes 15-59 | rise-heavy events, minutes 0-14 | rise-heavy events, minutes 15-59 |
|---|---|---|---|---|---|---|
| (0.35, 0.5] | 0.1534 | 0.0298 | 0.1534 | 0.0223 |  | 0.0443 |
| (0.5, 0.7] | 0.1212 | 0.026 | 0.1266 | 0.0214 |  | 0.0309 |
| (0.7, 1.0] | 0.0789 | 0.0212 | 0.0969 | 0.0196 | 0.0487 | 0.0224 |
| (1.0, 1.4] | 0.0389 | 0.0145 | 0.0576 | 0.0135 | 0.0181 | 0.0154 |
| (1.4, 2.0] | 0.0267 | 0.0063 | 0.0513 | 0.0124 | 0.0049 | 0.0043 |
| (2.0, 2.83] |  | 0.0128 |  |  |  | 0.0068 |

**The four rules:**

| test | what | value | rule | met |
|---|---|---|---|---|
| B1 | IQR of the transition level of the plateau-then-decays tails (dB) | 0.35 | < 0.3 | False |
| B2 | the same link in different events: permutation p | 0.35 | < 0.05 | False |
| B3 | Spearman of the plateau duration with A_start / with the rain amount | 0.162 | > 0.3, interval excluding 0 | False |
| B5 | slope of log k vs log A, decaying tails of the events not rise-heavy | 0.284 | < -0.3 | False |

**H2 verdict by the rule: NOT SUPPORTED** (0 of 4 met, 2 needed). B4 has no rule; it is reported above.

## 5. Part C: the two code fixes, the fake check, the fits again

**C1.** `pooled_fit.py`, option `bound_w0`: the start W0 of F1 and F4 is kept inside [0.001, 30] dB, as it always was
for F3, F5, F6 and the free curve. **C2.** option `noise_mult` / `fit(..., fit_noise=True)`: one fitted number m
multiplies the noise variance of every tail of the dataset; the likelihood gets the term -1/2 n log m. Both options
are OFF by default, so steps 3b-1b to 4 give exactly the numbers they gave before (checked: step 4 re-run from its
cache, `python -m research compare`).

**Check on fake data** (6 cached datasets of step 3b-3b; "before" = the cached menu fits of step 3b-3c part C).
Rule: same winner, the parameters of the truth's own candidate within +-30 %, the multiplier inside
0.8..1.25. **NOT PASSED: 4 of 6.**

| truth | replicate | winner_before | winner_after | margin_before_bic | margin_after_bic | own_candidate | parameters_before | parameters_after | largest_change_pct | s_before | s_after | noise_multiplier | multiplier_all_candidates | passed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| m1 | 0 | F1 | F1 | 5.8 | 6.4 | F1 | tau 25.57 | tau 25.59 | 0.1 | 1.207 | 1.209 | 1.05 | 1.05..1.05 | True |
| m1 | 1 | F6 | F5 | 10.1 | 1.0 | F1 | tau 27.26 | tau 26.98 | 1.0 | 1.205 | 1.204 | 1.075 | 1.07..1.08 | False |
| m3 | 0 | F5 | F5 | 0.7 | 0.4 | F3 | k 0.05534, n 1.987 | k 0.05742, n 1.989 | 3.8 | 1.299 | 1.228 | 1.034 | 1.03..1.07 | True |
| m3 | 1 | F3 | F3 | 2.2 | 2.2 | F3 | k 0.0534, n 2.041 | k 0.05346, n 2.053 | 0.6 | 1.278 | 1.275 | 1.062 | 1.06..1.09 | True |
| m6 | 0 | F6 | F6 | 40.3 | 25.8 | F6 | k_run 0.04246, k_evap 0.003074, A_c 0.1539 | k_run 0.04182, k_evap 0.01607, A_c 0.1528 | 422.8 | 1.159 | 1.285 | 1.049 | 1.05..1.05 | False |
| m6 | 1 | F6 | F6 | 68.0 | 63.5 | F6 | k_run 0.04787, k_evap 0.003449, A_c 0.1826 | k_run 0.04964, k_evap 0.003738, A_c 0.178 | 8.4 | 1.266 | 1.237 | 1.063 | 1.06..1.06 | True |

**The real sets with both fixes** (first 60 min; every candidate now has the same start bound, so the BIC column IS
like for like; k of the BIC includes the multiplier):

| dataset | what | tails | events | w_hi_db | menu_winner | menu_family | runner_up | margin_bic | next_family | margin_family_bic | f1_minus_winner_bic | noise_multiplier_winner |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| main | v1 main | 429 | 24 | 2.83 | F6 | state switch | F4 | 11.6 | F4 | 11.6 | 16.7 | 4.32 |
| no_rise_heavy | without event_rise_heavy | 265 | 15 | 2.4 | F6 | state switch | F4 | 21.2 | F4 | 21.2 | 28.1 | 3.54 |
| decays | class 'decays' only | 240 | 17 | 1.96 | F6 | state switch | F4 | 15.5 | F4 | 15.5 | 28.4 | 2.61 |
| decays_plus | 'decays' + 'plateau-then-decays' from the plateau's end | 304 | 21 | 2.01 | F6 | state switch | F4 | 23.0 | F4 | 23.0 | 37.2 | 2.69 |

| dataset | f1_tau_min | f3_n | f5_b | s_winner | curve_slope | curve_reading | curve_noise_multiplier | curve_minus_winner_bic | curve_start_cost_spread | k_at_knots_per_min | knots_db |
|---|---|---|---|---|---|---|---|---|---|---|---|
| main | 19.3 | 1.0 | 4.23e-07 | 1.296 | -0.469 | falling (no fake truth looked like this) | 3.82 | -2347.4 | 355.0 | 0.0347 0.0460 0.0265 0.0498 0.0069 | 0.35 0.6 0.907 1.3 2.83 |
| no_rise_heavy | 14.9 | 1.0 | 4.01e-06 | 1.15 | -0.506 | falling (no fake truth looked like this) | 3.05 | -1910.3 | 236.1 | 0.0385 0.0482 0.0328 0.0380 0.0121 | 0.35 0.53 0.79 1.19 2.4 |
| decays | 10.2 | 1.0 | 1.05e-06 | 0.769 | 0.005 | exponential in this range | 2.33 | -1546.0 | 54.6 | 0.0459 0.0557 0.0336 0.0613 0.0327 | 0.35 0.488 0.67 1.02 1.96 |
| decays_plus | 10.8 | 1.0 | 6.79e-06 | 0.817 | -0.145 | exponential in this range | 2.44 | -1556.7 | 71.0 | 0.0565 0.0704 0.0339 0.0586 0.0327 | 0.35 0.49 0.65 0.957 2.01 |

Against step 4 (old code):

| dataset | step4_winner_as_run | step4_like_for_like_winner | step4_like_for_like_margin_bic | menu_winner | margin_family_bic | step4_curve_slope | curve_slope | step4_f1_tau_min | f1_tau_min |
|---|---|---|---|---|---|---|---|---|---|
| main | F1 | F6 | 21.5 | F6 | 11.6 | -0.872 | -0.469 | 25.2 | 19.3 |
| no_rise_heavy | F6 | F6 | 51.8 | F6 | 21.2 | -0.054 | -0.506 | 15.3 | 14.9 |

All menu fits:

| dataset | candidate | family | n_tails | k | ll | bic | dbic | parameters | s | noise_mult | bic_spread_starts |
|---|---|---|---|---|---|---|---|---|---|---|---|
| main | F6 | state switch | 429 | 5 | -32231.1 | 64490.73 | 0.0 | k_run 0.04869, k_evap 0.005095, A_c 0.09425 | 1.296 | 4.32 | 28.2 |
| main | F4 | two layers | 429 | 5 | -32236.9 | 64502.32 | 11.59 | tau_f 18.74, tau_s 2203, phi 0.02943 | 1.302 | 4.317 | 31.9 |
| main | F1 | exponential | 429 | 3 | -32245.15 | 64507.42 | 16.69 | tau 19.29 | 1.301 | 4.316 | 4.1 |
| main | F3 | nonlinear in W (F3 / F5) | 429 | 4 | -32245.17 | 64513.15 | 22.42 | k 0.05227, n 1 | 1.304 | 4.317 | 28.9 |
| main | F5 | nonlinear in W (F3 / F5) | 429 | 4 | -32246.46 | 64515.74 | 25.01 | a 0.05472, b 4.227e-07 | 1.253 | 4.315 | 430.6 |
| no_rise_heavy | F6 | state switch | 265 | 5 | -18401.75 | 36829.47 | 0.0 | k_run 0.06673, k_evap 0.0005629, A_c 0.1187 | 1.15 | 3.537 | 76.0 |
| no_rise_heavy | F4 | two layers | 265 | 5 | -18412.35 | 36850.66 | 21.19 | tau_f 14.07, tau_s 473.1, phi 0.04217 | 1.107 | 3.554 | 23.4 |
| no_rise_heavy | F1 | exponential | 265 | 3 | -18420.97 | 36857.52 | 28.05 | tau 14.88 | 1.123 | 3.553 | 1.8 |
| no_rise_heavy | F3 | nonlinear in W (F3 / F5) | 265 | 4 | -18421.11 | 36862.99 | 33.52 | k 0.06844, n 1 | 1.129 | 3.554 | 254.0 |
| no_rise_heavy | F5 | nonlinear in W (F3 / F5) | 265 | 4 | -18423.96 | 36868.68 | 39.22 | a 0.08885, b 4.006e-06 | 1.17 | 3.546 | 1421.1 |
| decays | F6 | state switch | 240 | 5 | -14472.42 | 28970.67 | 0.0 | k_run 0.1054, k_evap 0.001365, A_c 0.09668 | 0.769 | 2.614 | 65.1 |
| decays | F4 | two layers | 240 | 5 | -14480.17 | 28986.17 | 15.5 | tau_f 9.766, tau_s 6330, phi 0.03199 | 0.798 | 2.611 | 23.5 |
| decays | F1 | exponential | 240 | 3 | -14491.8 | 28999.09 | 28.43 | tau 10.25 | 0.784 | 2.62 | 0.3 |
| decays | F3 | nonlinear in W (F3 / F5) | 240 | 4 | -14491.79 | 29004.24 | 33.57 | k 0.09674, n 1 | 0.789 | 2.614 | 1.8 |
| decays | F5 | nonlinear in W (F3 / F5) | 240 | 4 | -14491.8 | 29004.25 | 33.59 | a 0.09759, b 1.053e-06 | 0.784 | 2.614 | 1.7 |
| decays_plus | F6 | state switch | 304 | 5 | -18401.31 | 36829.57 | 0.0 | k_run 0.09524, k_evap 0.0007559, A_c 0.09718 | 0.817 | 2.688 | 48.6 |
| decays_plus | F4 | two layers | 304 | 5 | -18412.81 | 36852.57 | 23.0 | tau_f 10.47, tau_s 1986, phi 0.03278 | 0.838 | 2.687 | 28.5 |
| decays_plus | F1 | exponential | 304 | 3 | -18425.32 | 36866.8 | 37.24 | tau 10.81 | 0.848 | 2.689 | 0.7 |
| decays_plus | F3 | nonlinear in W (F3 / F5) | 304 | 4 | -18425.27 | 36872.1 | 42.53 | k 0.09057, n 1 | 0.855 | 2.688 | 3.8 |
| decays_plus | F5 | nonlinear in W (F3 / F5) | 304 | 4 | -18425.43 | 36872.41 | 42.85 | a 0.09078, b 6.792e-06 | 0.857 | 2.685 | 0.8 |

Reading rules as in step 4: |slope| < 0.15 exponential; > 0.6 quadratic family; in between mildly
rising; negative = falling (no fake truth). F3 and F5 are one family; F4 against F6 by less than 10 BIC is not decisive.
The two per-class sets are chosen by the SHAPE of the tails, not by a fit, but that shape is the thing the law
describes: "decays" keeps only tails that halve within 30 min. Their fits describe those tails and are biased
towards fast drying; they are not the drying law of all tails.

## 6. In plain words: which explanation the data support

**Neither H1 nor H2 passes its rule. What the data do show is a third picture: a tail stalls or rises because of
WHEN it happens (the event), not because of WHICH link it is, how long the link is, or how much water is on it.**

- **H2 as "the rate per dB falls with the level": not supported, and fairly sure.** In the tails that decay, the rate
  does not fall with the level: slope +0.28 (event bootstrap +0.05..+0.77) in the events not flagged rise-heavy, +0.10
  in all events. The falling curve of step 4 (-0.80 here, model-free) appears only when tails that wait are mixed in
  (rise-heavy events: -1.69). The same link does not repeat its class in another event (B2, p = 0.35), and the plateau
  is not longer when the tail starts higher or the event had more rain (B3, at most +0.16). After their plateau, the
  plateau tails follow the same curve as the decaying tails just as well (rms 0.117 dB against 0.124 dB), and their
  plateau lies BELOW that curve in 45 of 61 tails: the tail waited, then dried like the others. A wait is a delay, not
  a slow part of one common curve. Limit of this test: "decays" is defined by a fast fall, so a tail that is slow
  because it is high would not be in that class at all (see 9).
- **H1 as "the references show water after t0_trace": not supported, but the test is weak.** The tails that do not
  decay have LESS rain nearby than the decaying ones (any rain indicator: 78 % against 93 %, odds ratio 0.29, interval
  0.10..0.58), not more. But the radar on the path is dry in all 429 tails by the way the tail set was built, and the
  humidity is missing for the two events that hold 60 of the 189 non-decaying tails and 22 of the 36 rising tails.
  So drizzle seen by the references is unlikely; dew or condensation could not be tested where it matters most.
- **What is supported: the cause is at the scale of the event.** The share of non-decaying tails runs from 0 % to
  100 % between events (19 % and 24 % in the two biggest clean events, 83 % in E20150617_0705, 100 % in
  E20150701_0030). The other links of the same event at the same time show the same behaviour more often (ratio 1.85,
  interval 1.27..2.43; the rule asked for 2). The two directions of a link agree in 81 % (64 % expected by chance).
  The rise does not grow with the path length (Spearman 0.00), and short and long links stall equally often (44 %),
  so the extra attenuation is not rain along the path. Rising tails are night tails half of the time (decaying: 8 %).
- **For which tails.** 240 of 429 decay at once (56 %). 64 wait 10 min or more, then decay normally (15 %). 36 rise by
  more than 0.3 dB (8 %; 22 of them in two night events without humidity data). 89 are flat or slow (21 %).
- **How sure.** Sure: the stall is not a property of the level, the link or the path length. Not sure: what it is. A
  weather cause fits (something still wets the antennas without rain on the references, or nothing dries them: night,
  cool air), but with one weather station, hourly values and 24 events it is a description, not a proven cause.

## 7. What this means for the drying law

- **h(W) does not need a falling k.** Do not build "slower when wetter" into the drying law to explain the stalled tails.
- **Which tails should define h(W): the decaying phase.** Two candidate sets are ready in
  `outputs/products/tails_model_classes.csv`: class "decays" (240 tails), and "decays" plus the plateau tails from the end
  of their plateau (304 tails). On both the free curve above 0.35 dB is flat (slopes +0.01 and -0.15, "exponential in
  this range") and F1 gives tau = 10.2 and 10.8 min; on all tails tau is 19.3 min. Warning: these sets keep the tails
  that fall fast, so tau is biased short, and the fits were made with a noise option whose fake-data check did not
  pass (section 5).
- **The stalled and rising minutes need something that is not in h(W).** For the 36 rising tails that is an INPUT
  (the wetting side g), but the rain references show nothing there, so rain R cannot be that input; humidity, dew or
  time of day would have to be. For the tails that only wait, a start time of drying per tail (or a rate that depends
  on the weather) describes them as well as an input does; the data here cannot tell these apart. That question is
  Step 6 (what drying depends on), which should now include night / day, humidity and wind from the start.
- **A new doubt about the "decays" class.** These tails follow rain more closely (rain nearby in 93 %, the gauge wet
  in 41 % against 20 %), and several fall steeply in the first minutes (gallery). Part of that fall may be the last
  path rain fading, not the antenna drying (step 2 already could not prove the fast drop is WAA). If so, the short
  tau of 10 min is too short for the antenna.
- **The menu with both fixes says the same thing in all four sets**: F6 first, F4 second (11.6 to 23 BIC), F1 third (17
  to 37 BIC behind F6); F3 at n = 1 and F5 at b = 0. F6 places its switch at 0.09 to 0.12 dB and F4 gives 3 to 4 % of
  the water a very slow layer: both add a small floor far below 0.35 dB, the range where the fake tests found a false
  feature for every truth. On one of the two fake exponential datasets F6 beat F1 by about 10 BIC with the old code. So read it as
  "exponential, plus something near 0.1 dB that the method cannot name", not as a state switch.
- **The noise during drying is 2.6 to 4.3 times the dry-weather noise variance**, also in the cleanest tails (2.6 in
  "decays"). The dry-weather noise model is too optimistic for drying tails whatever set is chosen.

## 8. Decisions for the learners (not taken)

1. Which tails define h(W): all of v1; "decays" (240); "decays" + plateau tails from the plateau's end (304); or all
   tails with the start of drying as an unknown per tail. A new frozen set is needed for anything but the first.
2. Whether the model gets a per-tail start delay, a weather-dependent rate (Step 6), or an input term without rain.
3. C1 (start bound for F1 / F4): make it the default? It removes the fault of step 4; earlier fake-data steps would
   then have to be re-run to keep their numbers honest (their tails never rise, so little should change).
4. C2 (fitted noise multiplier): the fake check did NOT pass its rule (4 of 6). Accept the two failures (see 9), fix the
   rule, or keep the multiplier out?
5. Get humidity for E20150617_0705 and E20150701_0030 from another station (they hold most of the rising tails).
   Without it, dew cannot be tested.
6. Test whether the fast fall of the "decays" tails is the end of the path rain (compare by path length with enough
   power, or use the two directions and neighbours), before the 10 min tau is called a drying time.
7. Importance levels no longer exist in this repository, so there is no key figure to promote or demote. The figure
   that carries the result is `diagnose_stalled_tails_b5_rate_per_level.png`.

## 9. Unclear, or possibly wrong

- **The class order is mine.** The spec did not say which class wins when two rules hold; I tried them in the listed
  order. The overlap counts are in section 2.
- **B1 cannot work as defined**: the transition level is about 0.8 x A_start by construction, so its IQR (0.35 dB)
  measures the spread of A_start (0.52 dB in these tails). It fails its rule, but it would be uninformative either way.
- **B5 and B4 are biased against H2 by the classes.** "decays" requires half of A_start within 30 min. If the truth
  were "slow when high", high tails would be classed as non-decaying and removed from the very table that tests it.
  The split by A_start inside the clean events (-0.12 for A_start >= 1.5 dB, -0.54 below) does not settle it.
- **The early / late table (added after the first look) is not a clean test either.** k at the same level is 3 to 7
  times higher in minutes 0-14 than later, but a tail that is already low after a few minutes is a fast tail, so the
  spread of speeds between tails gives this pattern too.
- **The dew indicator was corrected during the work**: the hourly product says "no dew risk" where humidity is
  missing; the first version counted that as "no". Now it is "unknown" (69 tails left out). The rule and its limits
  were not changed. "Radar nearby within 5 km" is my addition to the indicators the spec listed (fixed before the run).
- **"Any moisture evidence" is true for 98 % of both groups**, so that rule could hardly have been met: with dew risk
  and far gauges included, almost every hour after rain counts as moist.
- **Coherence 1.85 with interval 1.27..2.43**: below the rule's 2, clearly above 1. I report "not met" as the rule
  says; a reader may fairly call it partial support for a shared, weather-scale cause.
- **Part C, fake check failed by the rule (4 of 6).** m1 replicate 1: the winner changed from F6 to F5, both wrong
  and both only about 10 BIC or less ahead of F1; the tau of F1 moved 1 %. m6 replicate 0: same winner, but k_evap moved
  from 0.0031 to 0.016 / min, while two starts with k_evap 0.0023 and 0.016 end 0.2 log-likelihood apart: the first
  hour does not determine k_evap (known since 3b-1c). The multiplier is 1.03 to 1.08 in all six, inside the range but
  always above 1. My reading: the fixes did not break the tool, but the rule says "not passed", so the C3 numbers are
  not validated.
- **The free curve is still unstable on real tails**: its starts end 55 to 355 apart in cost (fake data: about 8), the
  cross-validation picks a different penalty per set, and the slope without the rise-heavy events moved from -0.05
  (step 4) to -0.51 here with the noise multiplier. Only the sign pattern is worth reading: falling when stalled
  tails are in, flat when they are out. The curve beats every menu law by 1500 to 2300 BIC, and it does so below
  0.35 dB, where it drops to k = 0.002 to 0.010 / min at 0.05 to 0.1 dB; that range is not read.
- **The master curve ignores the per-tail baseline offset** and uses the 5 min median of A, so 0.1 dB of misfit is
  about the noise; "on the curve" was set to the median misfit (0.124 dB), my choice.
- **One weather station, hourly**, for a city 30 km wide: humidity, wind and "night" are the same for every link of
  an event, so they cannot explain differences inside an event.
- **Earlier steps are unchanged**: both options are off by default; step 4 re-run from its cache gives identical tables.


## Files

- Per-tail table (classes, outside evidence, B4 placement): `research/wet_antenna/outputs/products/tails_model_classes.csv`
- Tables: `research/wet_antenna/outputs/tables/diagnose_stalled_tails_class_counts.csv`, `research/wet_antenna/outputs/tables/diagnose_stalled_tails_a_outside_evidence.csv`, `research/wet_antenna/outputs/tables/diagnose_stalled_tails_b_antenna_tests.csv`, `research/wet_antenna/outputs/tables/diagnose_stalled_tails_b4_master_curve.csv`,
  `research/wet_antenna/outputs/tables/diagnose_stalled_tails_b5_rate_per_level.csv`, `research/wet_antenna/outputs/tables/diagnose_stalled_tails_c_fake_check.csv`, `research/wet_antenna/outputs/tables/diagnose_stalled_tails_c_menu_fits.csv`, `research/wet_antenna/outputs/tables/diagnose_stalled_tails_c_summary.csv`
- Figures: `research/wet_antenna/outputs/figures/diagnose_stalled_tails_gallery.png` (4 tails per class with radar and gauge under A),
  `research/wet_antenna/outputs/figures/diagnose_stalled_tails_b4_collapse.png` (B4), `research/wet_antenna/outputs/figures/diagnose_stalled_tails_b5_rate_per_level.png` (B5)
