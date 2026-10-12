# Baseline estimation, sample run (openmrg_v1)

Settings: {'X_H': 1.0, 'WINDOW_H': 6.0, 'MIN_DRY_H': 1.0, 'DIFF_DB': 0.5, 'EXC_DB': 1.0}
Sample: 12 links / 24 sub-links x 18 significant train events = 432 sub-link x event pairs.

## 1. Coverage (how often we get a baseline at all)
| method | pairs |
| --- | --- |
| A | 380 |
| missing | 28 |
| B | 24 |

Flags (no_before = the missing pairs; the others are A pairs worth a look):
| flag | pairs |
| --- | --- |
| none | 363 |
| no_before | 28 |
| no_after | 28 |
| excursion_before | 7 |
| excursion_after | 5 |
| excursion_before;excursion_after | 1 |

## 2. Precision (spread of the dry TL inside the before window, dB)
median 0.09, 90th percentile 0.22. Compare with the RSL quantization step of ~0.3 dB.

## 3. Closure: after - before (how much the dry level moved over the event, dB)
|closure| median 0.15, 90th percentile 0.46; share above DIFF_DB = 0.5 dB: 8%. For method A this is the baseline error at the event end.

## 4. Out-of-sample: mean A on usable dry bins inside the event (should be ~0, dB)
pairs with such bins: 402; median 0.05, 10th / 90th percentile -0.18 / 0.42.

## 5. Sensitivity to X (|before(X) - before(1 h)|, dB)
- X = 0.5 h: median 0.00, 90th percentile 0.01, max 0.13
- X = 2 h: median 0.00, 90th percentile 0.01, max 0.22
- X = 3 h: median 0.00, 90th percentile 0.03, max 0.33

## 6. By link class (median |closure| and median in-event dry A, dB)
| length_class | frequency_band | pairs | abs_closure | in_event_dry_A |
| --- | --- | --- | --- | --- |
| 1-3km | >=25GHz | 100 | 0.16 | 0.03 |
| <1km | >=25GHz | 100 | 0.09 | 0.05 |
| >3km | <25GHz | 102 | 0.24 | 0.09 |
| >3km | >=25GHz | 102 | 0.13 | 0.05 |
