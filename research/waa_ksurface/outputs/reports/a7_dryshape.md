# Diagnostic of the dry-path conflict

Written by `research/waa_ksurface/a7_dryshape.py` (decisions.md, "Diagnostic of the dry-path conflict", logged before the run). No refitting: K1 = the A7 fit, K0 = the refit of the A7 review. Real TRAIN data of the 14 learning events; VAL and TEST not opened.

Windows: the dry-path windows of review check (5) with the bins ending 5-30 min after the path went dry all present and measured A - c > 0.15 dB in each. Rate = -(slope of ln(A - c)), early on [5, 10, 15] min, late on [15, 20, 25, 30] min; ratio = late / early. Intervals: 90 % event bootstrap (1000 draws, seed 3901).

- k_wet gauge: 1073 windows kept, 757 with all three ratios (dropped: measured 310, K1 8, K0 10).
- k_wet radar: 1073 windows kept, 752 with all three ratios (dropped: measured 310, K1 17, K0 11).

## (1), (2) Ratio rate(late) / rate(early), median (90 %)

| k_wet estimate | class | windows / events | ratio measured | ratio K1 | ratio K0 |
|---|---|---|---|---|---|
| gauge | all | 757 / 13 | 0.227 (0.108 to 0.339) | 0.653 (0.611 to 0.687) | 1.000 (1.000 to 1.000) |
| gauge | >= 1 dB | 199 / 13 | 0.347 (0.228 to 0.516) | 0.585 (0.494 to 0.631) | 1.000 (1.000 to 1.000) |
| gauge | < 1 dB | 558 / 13 | 0.156 (0.040 to 0.283) | 0.682 (0.644 to 0.697) | 1.000 (1.000 to 1.000) |
| radar | all | 752 / 13 | 0.227 (0.108 to 0.341) | 0.695 (0.654 to 0.729) | 1.000 (1.000 to 1.000) |
| radar | >= 1 dB | 197 / 13 | 0.347 (0.226 to 0.516) | 0.613 (0.548 to 0.676) | 1.000 (0.998 to 1.000) |
| radar | < 1 dB | 555 / 13 | 0.159 (0.046 to 0.289) | 0.721 (0.689 to 0.745) | 1.000 (1.000 to 1.000) |

## (3) Errors: level at the first dry bin (model - measured, medians) and decay-rate ratio (model - measured)

| k_wet estimate | class | median start measured (dB) | level error K1 (dB) | level error K0 (dB) | rate error K1 | rate error K0 |
|---|---|---|---|---|---|---|
| gauge | all | 0.65 | -0.278 (-0.400 to -0.207) | -0.335 (-0.460 to -0.253) | 0.426 (0.308 to 0.542) | 0.773 (0.661 to 0.892) |
| gauge | >= 1 dB | 1.42 | -0.978 (-1.092 to -0.812) | -1.019 (-1.131 to -0.792) | 0.238 (0.078 to 0.374) | 0.653 (0.482 to 0.772) |
| gauge | < 1 dB | 0.55 | -0.196 (-0.230 to -0.149) | -0.263 (-0.299 to -0.203) | 0.525 (0.386 to 0.631) | 0.844 (0.717 to 0.960) |
| radar | all | 0.652 | -0.292 (-0.416 to -0.221) | -0.346 (-0.479 to -0.275) | 0.468 (0.344 to 0.591) | 0.773 (0.659 to 0.892) |
| radar | >= 1 dB | 1.418 | -0.979 (-1.086 to -0.808) | -1.017 (-1.137 to -0.803) | 0.266 (0.111 to 0.412) | 0.653 (0.480 to 0.774) |
| radar | < 1 dB | 0.55 | -0.201 (-0.252 to -0.165) | -0.266 (-0.320 to -0.224) | 0.562 (0.413 to 0.688) | 0.841 (0.711 to 0.954) |

## (4) The drying rates used

| k_wet estimate | K0: 1/k_dry (min) | K1 (A7): 1/k_dry at W = 0 dB (min) | K1 (A7): 1/k_dry at W = 0.5 dB (min) | K1 (A7): 1/k_dry at W = 1 dB (min) | K1 (A7): 1/k_dry at W = 2 dB (min) |
|---|---|---|---|---|---|
| gauge | 8.88 | 81.5 | 12.2 | 5.4 | 6.0 |
| radar | 13.22 | 95.1 | 17.1 | 8.4 | 8.0 |

