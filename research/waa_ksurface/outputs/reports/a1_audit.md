# A1: data audit and reproduction gate (TRAIN only)

Written by `research/waa_ksurface/a1_audit.py`. VAL and TEST attenuation was not opened.

## (a) The rebuilt dataset: PASS

206 links and 9780 TRAIN link-events (51 events); the gate needs 206 and 9780. Built with `system_id.build_dataset` called as in `system_id_s3.dataset`. The funnel (`research/waa_ksurface/outputs/tables/a1_audit_funnel.csv`):

| stage | count |
|---|---|
| pairs of the short TRAIN links in TRAIN events | 29998 |
| no baseline | 1148 |
| no wet bin of the link in the event | 13480 |
| sub-link not good | 1450 |
| use_in_waa_main false | 1971 |
| offset_std > 0.3 dB | 1139 |
| frequency outside the two bands | 862 |
| cut at a TSL step | 679 |
| TSL step before the rain | 16 |
| too few valid 5 min bins | 12 |
| too many 5 min bins without radar | 140 |
| link-events kept | 9780 |

## (b) The Part 1 length test reproduced: PASS

15 stored TRAIN cells (inputs: V1) of `research/wet_antenna/outputs/tables/final_s5a_length_test.csv`. Largest difference: I 0.00005 dB, s 0.00005 dB/km (allowed 0.01); interval ends 0.00005. Table: `research/waa_ksurface/outputs/tables/a1_audit_length_test_vs_part1.csv`.

| input | band | class | I_part1 | I_rebuilt | s_part1 | s_rebuilt | dI | ds | links_under_1km_rebuilt | eligible_part1 |
|---|---|---|---|---|---|---|---|---|---|---|
| V1 | 28 GHz | 0.5-1 | 1.0097 | 1.0097 | 0.1043 | 0.1043 | 0.0 | 0.0 | 17 | True |
| V1 | 28 GHz | 1-2 | 1.3703 | 1.3703 | 0.1851 | 0.1851 | 0.0 | 0.0 | 17 | True |
| V1 | 28 GHz | 2-5 | 1.9651 | 1.9651 | 0.3387 | 0.3387 | 0.0 | 0.0 | 17 | True |
| V1 | 28 GHz | 5-10 | 2.548 | 2.548 | 0.7411 | 0.7411 | 0.0 | 0.0 | 17 | True |
| V1 | 28 GHz | 10-20 | 3.0856 | 3.0856 | 1.1099 | 1.1099 | 0.0 | 0.0 | 17 | True |
| V1 | 32 GHz | 0.5-1 | 0.8701 | 0.8701 | 0.1837 | 0.1837 | 0.0 | 0.0 | 0 | False |
| V1 | 32 GHz | 1-2 | 1.2933 | 1.2933 | 0.2453 | 0.2453 | 0.0 | 0.0 | 0 | False |
| V1 | 32 GHz | 2-5 | 1.9679 | 1.9679 | 0.4181 | 0.4181 | 0.0 | 0.0 | 0 | False |
| V1 | 32 GHz | 5-10 | 2.7587 | 2.7587 | 0.6219 | 0.6219 | 0.0 | 0.0 | 0 | False |
| V1 | 32 GHz | 10-20 | 2.9306 | 2.9306 | 2.0244 | 2.0244 | 0.0 | 0.0 | 0 | False |
| V1 | 38 GHz | 0.5-1 | 1.0132 | 1.0132 | 0.237 | 0.237 | 0.0 | 0.0 | 24 | True |
| V1 | 38 GHz | 1-2 | 1.6007 | 1.6007 | 0.3167 | 0.3167 | 0.0 | 0.0 | 24 | True |
| V1 | 38 GHz | 2-5 | 2.1565 | 2.1565 | 0.7007 | 0.7007 | 0.0 | 0.0 | 24 | True |
| V1 | 38 GHz | 5-10 | 2.3134 | 2.3134 | 1.5651 | 1.5651 | 0.0 | 0.0 | 24 | True |
| V1 | 38 GHz | 10-20 | 3.4543 | 3.4543 | 2.0785 | 2.0785 | 0.0 | 0.0 | 23 | True |

**The Part 1 definition of a steady-rain 5 min value** (config `wet_antenna.system_id_s3.length_test`, code `length_test.steady_table`):
a 5 min bin of a link-event with a valid A (at least 3 of its 5 minutes) and a valid radar path-mean R, whose R lies in a rain
class (0.5-1, 1-2, 2-5, 5-10, 10-20 mm/h, [lo, hi)), where (1) the 20 min before it
(4 previous bins) are all wet (R >= 0.3 mm/h), and (2) the mean R of the 15 min before it
(3 previous bins, all with radar) is within 0.5-2.0 times the bin's R. Per link-event and class the kept bins are
reduced to one value: the median of (A - c) and the median of R.

**The Part 1 definition of an eligible cell** (a cell = band x rain class): the cell is fitted if it has at least
10 link-event values and 3 distinct link lengths (line (A - c) = I + s L, Huber k = 1.345); it is ELIGIBLE if it also has at least
8 distinct links shorter than 1 km.


## (c) The TRAIN data

Rain = radar / g_e (V1). A valid minute has A and radar; a rainy minute has R >= 0.3 mm/h. (`research/waa_ksurface/outputs/tables/a1_audit_counts.csv`)

| band_name | length_class | link_events | links | events | valid_minutes | rainy_minutes |
|---|---|---|---|---|---|---|
| 28 GHz | 1-3 km | 2273 | 45 | 50 | 2242650 | 334570 |
| 28 GHz | L < 1 km | 803 | 17 | 39 | 846040 | 127190 |
| 28 GHz | L > 3 km | 2594 | 51 | 51 | 2301580 | 329255 |
| 32 GHz | 1-3 km | 582 | 12 | 36 | 600875 | 89725 |
| 32 GHz | L > 3 km | 362 | 7 | 43 | 237250 | 27710 |
| 38 GHz | 1-3 km | 2080 | 46 | 50 | 2076520 | 312480 |
| 38 GHz | L < 1 km | 915 | 24 | 39 | 989910 | 150760 |
| 38 GHz | L > 3 km | 171 | 4 | 35 | 169725 | 26315 |

Valid minutes per rain class on the links below 1 km (`research/waa_ksurface/outputs/tables/a1_audit_short_minutes.csv`):

| rain_class_mmh | 28 GHz | 32 GHz | 38 GHz | all bands |
|---|---|---|---|---|
| 0-0.1 | 653960 | 0 | 761195 | 1415155 |
| 0.1-0.25 | 54450 | 0 | 65065 | 119515 |
| 0.25-0.5 | 32910 | 0 | 40865 | 73775 |
| 0.5-1 | 28275 | 0 | 31970 | 60245 |
| 1-2 | 29550 | 0 | 36590 | 66140 |
| 2-5 | 29515 | 0 | 34240 | 63755 |
| 5-10 | 14225 | 0 | 16585 | 30810 |
| 10-20 | 2510 | 0 | 2790 | 5300 |
| >= 20 | 645 | 0 | 610 | 1255 |

The gauge-input set and the g_e use (`research/waa_ksurface/outputs/tables/a1_audit_gauge_set.csv`):

| what | count |
|---|---|
| gauge-input set (City gauge within 1000 m of the path): links | 50 |
| gauge-input set: link-events | 2425 |
| gauge-input set: events | 44 |
| gauge-input set, L < 1 km: links | 18 |
| gauge-input set, L < 1 km: link-events | 749 |
| TRAIN events with link-events | 51 |
| events with a gauge ratio (g_e used) | 14 |
| events with g_e = 1 (fallback) | 37 |
| link-events in events with a gauge ratio | 5047 |
| link-events in events with g_e = 1 | 4733 |

## (d) The stored comparison fits (nothing refitted; `research/waa_ksurface/outputs/tables/a1_audit_stored_fits.csv`)

| model | path | loads | content |
|---|---|---|---|
| Schleiss b | research/wet_antenna/outputs/products/literature_s8c_fits.json | True | shape [309.236952, 76.92152], amp per band, c_band, a_link of 206 links |
| Pastorek b | research/wet_antenna/outputs/products/literature_s8c_fits.json | True | shape [0.566297], amp per band, c_band, a_link of 206 links |
| Part 1 equation | research/wet_antenna/outputs/products/waa_equation_final.json | True | every number of the equation, a_link of 206 links |
| no memory (B2 of S8) | research/.cache/wet_antenna/test_s8/fit_B2.pkl | False | MISSING: the Part 1 fit cache is not on this machine; only its sha256 and its error numbers are committed |
| no WAA (B0 of S8) | research/.cache/wet_antenna/test_s8/fit_B0.pkl | False | MISSING: the Part 1 fit cache is not on this machine; only its sha256 and its error numbers are committed |

## (e) The 32 GHz rule

Eligible cells per band (V1, Part 1 classes): 28 GHz: 5, 32 GHz: 0, 38 GHz: 5.

**Excluded from this branch for every model: 32 GHz** (no eligible length-test cell, so F and P cannot be measured).
