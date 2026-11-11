# T1: the rain error in steady windows (the error the length test meets)

Written by `research/waa_ksurface/t1_conditional_error.py`. Steadiness rule: 30 min within a factor 1.5. Part (a) reads rain data only (radar pixel at the 10 City gauges divided by g_e, the 1 min gauges), inside the 14 TRAIN events with a gauge ratio. VAL and TEST were not opened.

## (a) Radar label against gauge rain in steady 35 min windows

d = ln(label / gauge); label = the median of the 7 radar values of the window, gauge = its mean over the same 35 min. Used: label >= 1 mm/h and gauge > 0. `sigma_ev` = between events (the scale part), `sigma_w` = within events, `sigma_cond` = √0.54·sigma_w (the link-independent part), `sigma_shared` = √0.46·sigma_w.

| set | windows_label_ge_min | dropped_gauge_zero | windows | events | sd_total | median_d | sigma_ev | sigma_w | sigma_cond | sigma_shared |
|---|---|---|---|---|---|---|---|---|---|---|
| steady windows | 175 | 3 | 172 | 8 | 0.501 | 0.2 | 0.449 | 0.301 | 0.221 | 0.204 |
| all windows (reference) | 4280 | 183 | 4097 | 13 | 0.786 | 0.101 | 0.359 | 0.71 | 0.522 | 0.482 |

**[σc_lo, σc_hi] = [0.182, 0.256]** (source: (a) steady windows at the gauges, event bootstrap). The gauge's own error and the point-vs-pixel difference are not subtracted. **Not measured below 1 mm/h**: there a 35 min window holds only a few gauge steps.

By label class (steady windows):

| label_class_mmh | windows | events | gauges | median_d | sd_d |
|---|---|---|---|---|---|
| 1-2 | 57 | 5 | 9 | -0.082 | 0.593 |
| 2-5 | 38 | 6 | 10 | 0.044 | 0.624 |
| 5-10 | 77 | 4 | 10 | 0.303 | 0.239 |

(`research/waa_ksurface/outputs/tables/t1_conditional_error_summary.csv`, `research/waa_ksurface/outputs/tables/t1_conditional_error_windows.csv`, `research/waa_ksurface/outputs/tables/t1_conditional_error_result.csv`)

## (b) The count of steady bins (TRAIN link-events, gauge-ratio events)

Real data: 8551 steady bins in 14 events. Fake replicates with the observed radar as the true rain and a per-bin error added (46 % shared, lag-1 correlation 0.35):

| added_per_bin_sd | steady_bins | share_of_real | events_with_fewer_than_half | median_event_share |
|---|---|---|---|---|
| 0.0 | 8551.0 | 1.0 | 0 | 1.0 |
| 0.05 | 8177.667 | 0.956 | 0 | 0.955 |
| 0.1 | 7410.667 | 0.867 | 0 | 0.873 |
| 0.2 | 5490.0 | 0.642 | 0 | 0.681 |
| 0.3 | 3659.667 | 0.428 | 9 | 0.456 |
| 0.48 | 1552.667 | 0.182 | 14 | 0.201 |

**Caveat: the σ that matches the real count is an EFFECTIVE bound, biased low because the fake truth (the observed radar) already carries error.** With the observed radar as the truth, the count at added error 0 is the real count by construction, so the match is at 0. The table shows how fast the count falls: the real radar cannot carry much MORE fast error than it already has.

Second reading: 2.5 % of the rainy bins with 30 min of history are steady in the real data (345226 bins). A rain that is perfectly constant in time gives this share with a per-bin error of sd **0.509** (lag-1 correlation 0.35). This is an upper bound on the fast per-bin error of a link's own radar value: real changes of the rain use up part of the allowed unsteadiness. (`research/waa_ksurface/outputs/tables/t1_conditional_error_counts.csv`)

