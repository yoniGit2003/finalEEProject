# A4: coverage of the (R, W) plane and the learning domain

Written by `research/waa_ksurface/a4_coverage.py`. TRAIN, links below 1 km (28 and 38 GHz), the 14 events with a gauge ratio. **Ŵ = (A − c) − L·P̂(R) is used for this map only, never as the target of a fit.** Above P's valid range P̂ is a linear continuation, not a measurement. VAL and TEST were not opened.

## The data

| set | events | short_links | link_events | valid_minutes | rainy_minutes |
|---|---|---|---|---|---|
| events with a gauge ratio (the learning set) | 14 | 41 | 926 | 1497827 | 377161 |
| g_e = 1 events (not used) | 28 | 41 | 792 | 337981 | 43570 |
| gauge-input set inside the learning set | 14 | 18 | 377 | 616121 | 154832 |

## The floors

- Events: at least **5** (with 14 events a cell fed by 5 events is empty in 0.2 % of the bootstrap draws, by 4 events in 0.9 %).
- Minutes: σ_1 = 0.072 dB (sd of one 1 min value on dry minutes, median of 870 link-events; lower limit 0.087 dB). n_min = 5·ceil((σ_1 / 0.0375)²) = **30 minutes**.
- A cell is in the domain only if it passes both floors under both continuations of P (last-segment slope, and proportional to R).

## The domain

| band | regime | domain_cells | minutes_in_domain_L2 | R_classes | R_edge_mmh | W_lo_db | W_hi_db | cells_with_R_above_F_range |
|---|---|---|---|---|---|---|---|---|
| 28 GHz | drying | 132 | 484766 | < 0.1 .. 10-20 | 20.0 | 0.0 | 7.5 | 35 |
| 28 GHz | wetting | 41 | 183480 | < 0.1 .. 10-20 | 20.0 | -0.3 | 2.1 | 13 |
| 38 GHz | drying | 150 | 559101 | < 0.1 .. 10-20 | 20.0 | 0.0 | 8.1 | 39 |
| 38 GHz | wetting | 46 | 228659 | < 0.1 .. 10-20 | 20.0 | -0.3 | 2.1 | 16 |
| pooled | drying | 174 | 1045175 | < 0.1 .. 10-20 | 20.0 | 0.0 | 9.0 | 47 |
| pooled | wetting | 46 | 412386 | < 0.1 .. 10-20 | 20.0 | -0.3 | 2.1 | 16 |

Pooled domain: R up to the class edge **20 mm/h**, W up to **9.0 dB**. Frozen table: `research/waa_ksurface/outputs/tables/a4_coverage_domain.csv` (column `in_domain`). Figure: `research/waa_ksurface/outputs/figures/a4_coverage_map.png`.

### Pooled, wetting: minutes per domain cell (continuation L2)

| W bin from (dB) | < 0.1 | 0.1-0.2 | 0.2-0.5 | 0.5-1 | 1-2 | 2-5 | 5-10 | 10-20 |
|---|---|---|---|---|---|---|---|---|
| 1.8 | 0 | 0 | 0 | 0 | 0 | 2857 | 1243 | 175 |
| 1.5 | 0 | 0 | 0 | 0 | 1705 | 5097 | 1503 | 197 |
| 1.2 | 0 | 0 | 0 | 32 | 5274 | 5130 | 1265 | 189 |
| 0.9 | 0 | 0 | 0 | 1437 | 7921 | 4459 | 1048 | 152 |
| 0.6 | 0 | 0 | 1328 | 6543 | 7973 | 3731 | 936 | 99 |
| 0.3 | 0 | 2414 | 10478 | 7359 | 5562 | 2828 | 732 | 72 |
| 0.0 | 18561 | 18377 | 16694 | 6295 | 4549 | 2000 | 561 | 113 |
| -0.3 | 220925 | 13966 | 12049 | 3532 | 2551 | 1492 | 794 | 188 |

### Pooled, drying: minutes per domain cell (continuation L2)

| W bin from (dB) | < 0.1 | 0.1-0.2 | 0.2-0.5 | 0.5-1 | 1-2 | 2-5 | 5-10 | 10-20 |
|---|---|---|---|---|---|---|---|---|
| 8.7 | 0 | 0 | 0 | 0 | 0 | 0 | 42 | 31 |
| 8.4 | 43 | 0 | 0 | 0 | 0 | 0 | 48 | 0 |
| 8.1 | 30 | 0 | 0 | 0 | 0 | 31 | 46 | 43 |
| 7.8 | 51 | 0 | 0 | 0 | 0 | 36 | 65 | 42 |
| 7.5 | 55 | 0 | 0 | 0 | 0 | 46 | 68 | 47 |
| 7.2 | 44 | 0 | 0 | 0 | 0 | 65 | 97 | 46 |
| 6.9 | 60 | 0 | 0 | 0 | 0 | 71 | 119 | 59 |
| 6.6 | 77 | 0 | 0 | 0 | 31 | 101 | 136 | 67 |
| 6.3 | 83 | 0 | 0 | 0 | 50 | 118 | 171 | 74 |
| 6.0 | 101 | 0 | 30 | 0 | 57 | 163 | 211 | 71 |
| 5.7 | 145 | 0 | 41 | 0 | 76 | 210 | 224 | 104 |
| 5.4 | 231 | 30 | 56 | 0 | 91 | 288 | 334 | 106 |
| 5.1 | 372 | 35 | 49 | 44 | 118 | 370 | 416 | 126 |
| 4.8 | 451 | 42 | 70 | 56 | 152 | 517 | 586 | 142 |
| 4.5 | 700 | 58 | 79 | 72 | 235 | 766 | 727 | 176 |
| 4.2 | 822 | 74 | 144 | 109 | 328 | 1147 | 969 | 192 |
| 3.9 | 1007 | 107 | 198 | 156 | 575 | 1615 | 1375 | 198 |
| 3.6 | 1343 | 149 | 283 | 260 | 763 | 2311 | 1805 | 227 |
| 3.3 | 1581 | 182 | 402 | 403 | 1232 | 3088 | 2111 | 232 |
| 3.0 | 2115 | 287 | 596 | 637 | 1812 | 3820 | 2487 | 281 |
| 2.7 | 2241 | 485 | 883 | 884 | 2560 | 4808 | 2681 | 306 |
| 2.4 | 2776 | 757 | 1263 | 1362 | 3474 | 5605 | 2817 | 308 |
| 2.1 | 3758 | 1134 | 1951 | 1981 | 4277 | 5871 | 2437 | 286 |
| 1.8 | 4967 | 1620 | 2736 | 2715 | 5211 | 3001 | 685 | 81 |
| 1.5 | 9605 | 2392 | 4112 | 4037 | 4137 | 88 | 0 | 0 |
| 1.2 | 15452 | 4034 | 6552 | 5770 | 1832 | 0 | 0 | 0 |
| 0.9 | 28478 | 6403 | 10304 | 6625 | 148 | 0 | 0 | 0 |
| 0.6 | 52303 | 9322 | 12386 | 2861 | 0 | 0 | 0 | 0 |
| 0.3 | 117316 | 9358 | 3610 | 0 | 0 | 0 | 0 | 0 |
| 0.0 | 603532 | 53 | 0 | 0 | 0 | 0 | 0 | 0 |

## Rain above the measured ranges

Share of the rainy minutes (R >= 0.1 mm/h) of the short links whose R is above the valid range of F, and of P:

| band | rainy_minutes | share_above_F_range | share_above_P_range | F_max_mmh | P_max_mmh |
|---|---|---|---|---|---|
| 28 GHz | 172583 | 0.184 | 0.338 | 3.204 | 1.436 |
| 38 GHz | 204578 | 0.171 | 0.345 | 3.277 | 1.405 |
| pooled | 377161 | 0.177 | 0.342 |  |  |

## The gauge-input set (short links with a City gauge within 1 km), same floors

| band | regime | domain_cells | minutes_in_domain_L2 | R_classes | R_edge_mmh | W_lo_db | W_hi_db | cells_with_R_above_F_range |
|---|---|---|---|---|---|---|---|---|
| 28 GHz | drying | 94 | 159217 | < 0.1 .. 10-20 | 20.0 | 0.0 | 5.7 | 20 |
| 28 GHz | wetting | 36 | 61613 | < 0.1 .. 5-10 | 10.0 | -0.3 | 2.1 | 8 |
| 38 GHz | drying | 124 | 270708 | < 0.1 .. 10-20 | 20.0 | 0.0 | 7.5 | 32 |
| 38 GHz | wetting | 40 | 106225 | < 0.1 .. 10-20 | 20.0 | -0.3 | 2.1 | 11 |
| pooled | drying | 135 | 430888 | < 0.1 .. 10-20 | 20.0 | 0.0 | 7.8 | 33 |
| pooled | wetting | 41 | 167976 | < 0.1 .. 10-20 | 20.0 | -0.3 | 2.1 | 12 |

(`research/waa_ksurface/outputs/tables/a4_coverage_domain_gauge_set.csv`, `research/waa_ksurface/outputs/tables/a4_coverage_summary.csv`, `research/waa_ksurface/outputs/tables/a4_coverage_shares.csv`, `research/waa_ksurface/outputs/tables/a4_coverage_edges.csv`)

