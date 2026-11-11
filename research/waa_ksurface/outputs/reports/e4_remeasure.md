# E4: F(R) and P(R) re-measured (candidate (ii), instrument binning)

Written by `research/waa_ksurface/e4_remeasure.py`. TRAIN, 4563 link-events of 187 links in the 14 events with a gauge ratio (28 and 38 GHz). Bins selected by the event's area-mean radar rain, rule R2 (30 min within a factor 1.5); the line over L per band and class; the cell at the mean of the links' own rain. Only the 8 cells valid in the fake test (T3) are used. VAL and TEST were not opened.

## The cells

Value [bootstrap 90 %, 200 event draws] [error range: the value minus the fake-world error of the cell at the two ends of the measured rain error].

| band | class | R_cell (mm/h) | F = intercept (dB) | P = slope (dB/km) | links_short | links_long | events | link_events |
|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.142 | 0.375 [0.198, 0.574] [0.369, 0.375] | 0.045 [0.020, 0.078] [0.042, 0.046] | 17 | 96 | 11 | 2131 |
| 28 GHz | 0.2-0.5 | 0.367 | 0.525 [0.211, 0.928] [0.488, 0.526] | 0.062 [0.001, 0.094] [0.057, 0.067] | 17 | 96 | 7 | 1288 |
| 28 GHz | 0.5-1 | 0.639 | 0.725 [0.619, 0.826] [0.725, 0.808] | 0.093 [0.058, 0.140] [0.093, 0.119] | 17 | 96 | 9 | 1772 |
| 28 GHz | 1-2 | 1.436 | 1.406 [1.189, 1.624] [1.321, 1.406] | 0.098 [0.015, 0.172] [0.098, 0.144] | 17 | 96 | 8 | 1485 |
| 38 GHz | 0.1-0.2 | 0.133 | 0.350 [0.212, 0.501] [0.327, 0.350] | 0.075 [0.041, 0.126] [0.074, 0.076] | 24 | 50 | 11 | 1303 |
| 38 GHz | 0.2-0.5 | 0.286 | 0.619 [0.272, 1.196] [0.577, 0.619] | 0.096 [-0.007, 0.142] [0.092, 0.100] | 24 | 50 | 7 | 817 |
| 38 GHz | 0.5-1 | 0.653 | 0.863 [0.716, 0.966] [0.863, 0.936] | 0.122 [0.047, 0.222] [0.122, 0.163] | 24 | 50 | 9 | 1068 |
| 38 GHz | 1-2 | 1.405 | 1.660 [1.189, 2.101] [1.587, 1.660] | 0.179 [-0.019, 0.388] [0.179, 0.245] | 24 | 50 | 8 | 928 |

All 8 valid cells are eligible on real data.

- **28 GHz: valid range 0.14–1.44 mm/h.**
- **38 GHz: valid range 0.13–1.41 mm/h.**

**Below the lowest cell the curves are an interpolation to zero, not a measurement. Above the highest cell (about 1.5 mm/h) nothing is measured:** the 2–5 mm/h cells failed the fake test and no 5–10 mm/h cell could be judged. Curves: `research/waa_ksurface/outputs/tables/e4_remeasure_curves.csv` (column `measured`), figure `research/waa_ksurface/outputs/figures/e4_remeasure_F_P.png`.

## Against the superseded A3a values (information)

| band | class | R_cell | F_new | F_a3_superseded | shift_db |
|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.142 | 0.375 | 0.399 | -0.024 |
| 28 GHz | 0.2-0.5 | 0.367 | 0.525 | 0.703 | -0.178 |
| 28 GHz | 0.5-1 | 0.639 | 0.725 | 0.931 | -0.206 |
| 28 GHz | 1-2 | 1.436 | 1.406 | 1.367 | 0.04 |
| 38 GHz | 0.1-0.2 | 0.133 | 0.35 | 0.411 | -0.061 |
| 38 GHz | 0.2-0.5 | 0.286 | 0.619 | 0.684 | -0.065 |
| 38 GHz | 0.5-1 | 0.653 | 0.863 | 0.955 | -0.091 |
| 38 GHz | 1-2 | 1.405 | 1.66 | 1.583 | 0.077 |

Largest shift 0.21 dB. (`research/waa_ksurface/outputs/tables/e4_remeasure_vs_a3.csv`)

## The per-link residuals across storms

Short links (< 1 km), residual from the cell lines of the valid cells, averaged per link and storm: 41 links, 533 link-storm values. Between-link sd **0.088 dB**; storm-to-storm sd inside a link 0.557 dB. Permutation test (link labels shuffled inside each band, 1000 times): p = 0.0979 (null median 0.016 dB, 95 % point 0.102 dB).

**Criterion (decisions.md, 2026-11-11): a per-link factor is needed if the between-link sd is >= 0.15 dB AND p < 0.05. Verdict: not needed.**

Report only: over 41 links, the mean residual in the two upper classes against the two lower classes has correlation 0.68 and slope 1.04 (a slope above 1 points to a factor on F, a slope near 1 to an offset). (`research/waa_ksurface/outputs/tables/e4_remeasure_link_residuals.csv`, `research/waa_ksurface/outputs/tables/e4_remeasure_link_test.csv`, `research/waa_ksurface/outputs/figures/e4_remeasure_link_residuals.png`)

