# Control for the within-window ratio

Written by `research/waa_ksurface/a7_dryshape_control.py` (decisions.md, "Control for the within-window ratio", logged before the run). Fake data only (the step 3c / A7 fake world), no fit. Windows and ratio exactly as in `a7_dryshape`.

Truths (K0): G: 1/k_dry 8.88 min, 1/k_wet 1.74 min, seed 3920; R: 1/k_dry 13.22 min, 1/k_wet 10.92 min, seed 3921.

## Ratios

| truth | offsets | class | windows kept (A - c > 0.15) | windows with ratio | events | median ratio (90 %) |
|---|---|---|---|---|---|---|
| G | on | all | 480 | 413 | 14 | 0.296 (0.229 to 0.375) |
| G | on | >= 1 dB | 480 | 42 | 10 | 0.615 (0.528 to 0.685) |
| G | on | < 1 dB | 480 | 371 | 14 | 0.254 (0.200 to 0.313) |
| G | off | all | 20 | 18 | 8 | 0.098 (-0.134 to 0.407) |
| G | off | >= 1 dB | 20 | 3 | 2 | 0.545 (0.545 to 0.820) |
| G | off | < 1 dB | 20 | 15 | 7 | 0.000 (-0.498 to 0.260) |
| R | on | all | 562 | 435 | 13 | 0.436 (0.309 to 0.578) |
| R | on | >= 1 dB | 562 | 38 | 7 | 0.883 (0.674 to 0.996) |
| R | on | < 1 dB | 562 | 397 | 13 | 0.382 (0.284 to 0.450) |
| R | off | all | 76 | 72 | 10 | 0.721 (0.541 to 0.838) |
| R | off | >= 1 dB | 76 | 24 | 6 | 0.976 (0.798 to 1.007) |
| R | off | < 1 dB | 76 | 48 | 9 | 0.568 (0.395 to 0.681) |

## Reading (pre-registered; offsets on)

| truth | fake interval (all) contains 1 | real 0.347 below fake interval (>= 1 dB) | real 0.156 below fake interval (< 1 dB) | real 0.227 below fake interval (all) | reading | real 0.159 below fake interval (< 1 dB) |
|---|---|---|---|---|---|---|
| G | False | True | True | True | not explained by offset and selection |  |
| R | False | True |  | True | not explained by offset and selection | True |

