# Review of A7: the five checks

Written by `research/waa_ksurface/a7_review.py` (decisions.md, "Review of A7: the five checks", logged before the run). The A7 result is not changed; its cache is read only. Real TRAIN data of the 14 learning events (A7 flag: 0 removed). VAL and TEST were not opened.

## (1) Order

See decisions.md: the addendum (f58a0f5 19:51:30, a6e868a 19:52:02) was committed before the first real read of A7 (flag launcher 19:52:08, flag.pkl 19:52:20). No deviation.

## (2) The rung at each k_wet estimate

Code check: REPRODUCED (largest relative difference of the losses 2.2e-16, of the K1 gain 2.8e-15; rung K1).

| F curve | k_wet estimate | 1/k_wet (min) | drying rung | K1 gain (dB) | paired 90 % interval | margin (dB) | 1/k_dry at W = 0 dB (min) | 1/k_dry at W = 0.5 dB (min) | 1/k_dry at W = 1 dB (min) | 1/k_dry at W = 2 dB (min) |
|---|---|---|---|---|---|---|---|---|---|---|
| F as measured (A7) | gauge | 1.74 | K1 | 0.02013 | +0.01525 to +0.02419 | 0.00291 | 81.5 | 12.2 | 5.4 | 6.0 |
| F as measured (A7) | radar | 10.92 | K1 | 0.01837 | +0.01378 to +0.02293 | 0.00291 | 94.9 | 17.1 | 8.4 | 8.0 |

**Gauge estimate: K1; radar estimate: K1. Result (pre-registered rule): K1.** Caveat stated with it (Deviation 15): the selection with k_wet held at one value is the procedure of step 3d, which gave a false K1 on truth (a); it is run here because the review asks for it.

## (3) Event bootstrap of the selection

The first 50 draws of the A7 bootstrap (seed 3900) per estimate (user's decision; 200 with the full selection would take about 40 h), the full stage-2 selection in each draw. The 90 % intervals of 1 / k_dry are those of the 200-draw A7 bootstrap (refits of the selected rung), restated. Every draw: `research/waa_ksurface/outputs/tables/a7_review_boot.csv`.

| k_wet estimate | draws | select K1 | share K1 | K1 gain, median (5-95 %) (dB) | 1/k_dry W = 0.5 dB, A7 200 draws 5-95 % (min) | 1/k_dry W = 1 dB, A7 200 draws 5-95 % (min) | 1/k_dry W = 2 dB, A7 200 draws 5-95 % (min) |
|---|---|---|---|---|---|---|---|
| gauge | 50 | 50 | 1.0 | +0.0177 (+0.0107 to +0.0248) | 8.7 to 17.7 | 4.2 to 7.3 | 4.4 to 9.3 |
| radar | 50 | 50 | 1.0 | +0.0173 (+0.0101 to +0.0239) | 12.1 to 25.0 | 6.6 to 10.9 | 5.2 to 11.8 |

## (4) Sensitivity to F (rule A3f)

F through the 5-95 % bootstrap ends of the E4 cells (I_lo, I_hi); P unchanged; the whole of (2) at each end.

| F curve | k_wet estimate | 1/k_wet (min) | drying rung | K1 gain (dB) | paired 90 % interval | margin (dB) | 1/k_dry at W = 0 dB (min) | 1/k_dry at W = 0.5 dB (min) | 1/k_dry at W = 1 dB (min) | 1/k_dry at W = 2 dB (min) |
|---|---|---|---|---|---|---|---|---|---|---|
| F as measured (A7) | gauge | 1.74 | K1 | 0.02013 | +0.01525 to +0.02419 | 0.00291 | 81.5 | 12.2 | 5.4 | 6.0 |
| F as measured (A7) | radar | 10.92 | K1 | 0.01837 | +0.01378 to +0.02293 | 0.00291 | 94.9 | 17.1 | 8.4 | 8.0 |
| F at the lower ends (I_lo) | gauge | 1.74 | K1 | 0.02013 | +0.01465 to +0.02551 | 0.00291 | 144.7 | 20.7 | 9.0 | 8.6 |
| F at the lower ends (I_lo) | radar | 6.84 | K1 | 0.0186 | +0.01329 to +0.02429 | 0.00291 | 163.5 | 24.5 | 11.0 | 9.8 |
| F at the upper ends (I_hi) | gauge | 1.75 | K1 | 0.02144 | +0.01677 to +0.02504 | 0.00291 | 73.7 | 7.8 | 3.0 | 3.5 |
| F at the upper ends (I_hi) | radar | 17.39 | K1 | 0.01751 | +0.01288 to +0.02152 | 0.00291 | 77.6 | 13.8 | 6.6 | 5.6 |

## (5) Dry-path check (a check, not a selection)

Link-events of the learning events; the path went dry at the end of the last wet bin if the next bin is radar-dry (< 0.1 mm/h), at least 10 wet minutes before; the 5 min bin that ends at the lag, the path dry throughout. Error = median(model) - median(measured).

| k_wet estimate | lag (min) | values | measured A - c (dB) | K1 (A7 fit) (dB) | K0 refit (dB) | error K1 (dB) | error K0 (dB) |
|---|---|---|---|---|---|---|---|
| gauge | 5 | 3502 | 0.26 | 0.361 | 0.297 | 0.101 | 0.037 |
| gauge | 10 | 3466 | 0.173 | 0.289 | 0.173 | 0.116 | -0.0 |
| gauge | 20 | 3430 | 0.1 | 0.21 | 0.058 | 0.11 | -0.042 |
| gauge | 30 | 3404 | 0.07 | 0.163 | 0.019 | 0.093 | -0.051 |
| radar | 5 | 3502 | 0.26 | 0.324 | 0.257 | 0.064 | -0.003 |
| radar | 10 | 3466 | 0.173 | 0.274 | 0.178 | 0.101 | 0.004 |
| radar | 20 | 3430 | 0.1 | 0.213 | 0.085 | 0.113 | -0.015 |
| radar | 30 | 3404 | 0.07 | 0.174 | 0.04 | 0.104 | -0.03 |

Execution: keep-awake True; lowest free memory 0.51 GB; low-memory events 21.

