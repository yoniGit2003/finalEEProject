# E4: F(R) and P(R) re-measured (candidate (ii), instrument binning)

Written by `research/waa_ksurface/e4_remeasure.py`. TRAIN, 4563 link-events of 187 links in the 14 events with a gauge ratio (28 and 38 GHz). Bins selected by the event's area-mean radar rain, rule R2 (30 min within a factor 1.5); the line over L per band and class; the cell at the mean of the links' own rain. **Validity is per quantity (Deviation 4, 2026-11-12):** F from the 10 F-valid cells, P from the 8 P-valid cells. VAL and TEST were not opened.

## Validity per quantity (the committed T3 fake results, re-read; no new fit)

Fake-world error of candidate (ii) at the two ends of the measured rain error. F-valid: |intercept error| <= 0.15 dB at both ends. P-valid: |slope error| <= 0.0710 dB/km at both ends. (`research/waa_ksurface/outputs/tables/e4_remeasure_validity.csv`)

| band | class | I_err_lo | I_err_hi | s_err_lo | s_err_hi | valid_T3_both | F_valid | P_valid |
|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.0058 | 0.0036 | -0.0007 | 0.0029 | True | True | True |
| 28 GHz | 0.2-0.5 | -0.0017 | 0.0364 | 0.005 | -0.0053 | True | True | True |
| 28 GHz | 0.5-1 | -0.0835 | -0.0097 | -0.0258 | -0.026 | True | True | True |
| 28 GHz | 1-2 | 0.0627 | 0.0859 | -0.044 | -0.0458 | True | True | True |
| 28 GHz | 2-5 | 0.0792 | 0.0262 | -0.0714 | -0.0528 | False | True | False |
| 38 GHz | 0.1-0.2 | 0.0119 | 0.023 | 0.0004 | -0.0011 | True | True | True |
| 38 GHz | 0.2-0.5 | 0.0206 | 0.0421 | 0.0038 | -0.0046 | True | True | True |
| 38 GHz | 0.5-1 | -0.0724 | -0.053 | -0.0406 | -0.0309 | True | True | True |
| 38 GHz | 1-2 | 0.0435 | 0.0728 | -0.0523 | -0.0654 | True | True | True |
| 38 GHz | 2-5 | 0.0652 | -0.0291 | -0.1376 | -0.1106 | False | True | False |

## The cells

Value [bootstrap 90 %, 200 event draws] [error range: the value minus the fake-world error of the cell at the two ends of the measured rain error].

| band | class | R_cell (mm/h) | F = intercept (dB) | P = slope (dB/km) | links_short | links_long | events | link_events |
|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.142 | 0.375 [0.198, 0.574] [0.369, 0.375] | 0.045 [0.020, 0.078] [0.042, 0.046] | 17 | 96 | 11 | 2131 |
| 28 GHz | 0.2-0.5 | 0.367 | 0.525 [0.211, 0.928] [0.488, 0.526] | 0.062 [0.001, 0.094] [0.057, 0.067] | 17 | 96 | 7 | 1288 |
| 28 GHz | 0.5-1 | 0.639 | 0.725 [0.619, 0.826] [0.725, 0.808] | 0.093 [0.058, 0.140] [0.093, 0.119] | 17 | 96 | 9 | 1772 |
| 28 GHz | 1-2 | 1.436 | 1.406 [1.189, 1.624] [1.321, 1.406] | 0.098 [0.015, 0.172] [0.098, 0.144] | 17 | 96 | 8 | 1485 |
| 28 GHz | 2-5 | 3.204 | 1.962 [1.589, 2.288] [1.882, 1.962] | not P-valid | 17 | 96 | 8 | 1444 |
| 38 GHz | 0.1-0.2 | 0.133 | 0.350 [0.212, 0.501] [0.327, 0.350] | 0.075 [0.041, 0.126] [0.074, 0.076] | 24 | 50 | 11 | 1303 |
| 38 GHz | 0.2-0.5 | 0.286 | 0.619 [0.272, 1.196] [0.577, 0.619] | 0.096 [-0.007, 0.142] [0.092, 0.100] | 24 | 50 | 7 | 817 |
| 38 GHz | 0.5-1 | 0.653 | 0.863 [0.716, 0.966] [0.863, 0.936] | 0.122 [0.047, 0.222] [0.122, 0.163] | 24 | 50 | 9 | 1068 |
| 38 GHz | 1-2 | 1.405 | 1.660 [1.189, 2.101] [1.587, 1.660] | 0.179 [-0.019, 0.388] [0.179, 0.245] | 24 | 50 | 8 | 928 |
| 38 GHz | 2-5 | 3.277 | 2.029 [1.653, 2.407] [1.964, 2.059] | not P-valid | 24 | 50 | 8 | 894 |

Every F-valid or P-valid cell is eligible on real data.

- **28 GHz: F valid range 0.14–3.20 mm/h (5 cells); P valid range 0.14–1.44 mm/h (4 cells).**
- **38 GHz: F valid range 0.13–3.28 mm/h (5 cells); P valid range 0.13–1.41 mm/h (4 cells).**

**Below the lowest cell the curves are an interpolation to zero, not a measurement. Above the highest valid cell of a quantity nothing is measured for it** (F held, P with its last slope: extrapolation, marked). No 5–10 mm/h cell could be judged in the fake test. Curves: `research/waa_ksurface/outputs/tables/e4_remeasure_curves.csv` (columns `F_measured`, `P_measured`), figure `research/waa_ksurface/outputs/figures/e4_remeasure_F_P.png`.

## Against the superseded A3a values (information)

| band | class | R_cell | F_new | F_a3_superseded | shift_db |
|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.142 | 0.375 | 0.399 | -0.024 |
| 28 GHz | 0.2-0.5 | 0.367 | 0.525 | 0.703 | -0.178 |
| 28 GHz | 0.5-1 | 0.639 | 0.725 | 0.931 | -0.206 |
| 28 GHz | 1-2 | 1.436 | 1.406 | 1.367 | 0.04 |
| 28 GHz | 2-5 | 3.204 | 1.962 | 1.962 | -0.001 |
| 38 GHz | 0.1-0.2 | 0.133 | 0.35 | 0.411 | -0.061 |
| 38 GHz | 0.2-0.5 | 0.286 | 0.619 | 0.684 | -0.065 |
| 38 GHz | 0.5-1 | 0.653 | 0.863 | 0.955 | -0.091 |
| 38 GHz | 1-2 | 1.405 | 1.66 | 1.583 | 0.077 |
| 38 GHz | 2-5 | 3.277 | 2.029 | 2.169 | -0.14 |

Largest shift 0.21 dB. (`research/waa_ksurface/outputs/tables/e4_remeasure_vs_a3.csv`)

## The per-link residuals across storms

Short links (< 1 km), residual from the cell lines of the 8 cells valid in BOTH quantities (the criterion as written before Deviation 4), averaged per link and storm: 41 links, 533 link-storm values. Between-link sd **0.088 dB**; storm-to-storm sd inside a link 0.557 dB. Permutation test (link labels shuffled inside each band, 1000 times): p = 0.0979 (null median 0.016 dB, 95 % point 0.102 dB).

**Criterion (decisions.md, 2026-11-11): a per-link factor is needed if the between-link sd is >= 0.15 dB AND p < 0.05. Verdict: not needed.**

Report only: over 41 links, the mean residual in the two upper classes against the two lower classes has correlation 0.68 and slope 1.04 (a slope above 1 points to a factor on F, a slope near 1 to an offset). (`research/waa_ksurface/outputs/tables/e4_remeasure_link_residuals.csv`, `research/waa_ksurface/outputs/tables/e4_remeasure_link_test.csv`, `research/waa_ksurface/outputs/figures/e4_remeasure_link_residuals.png`)

