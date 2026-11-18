# A7: the real fit (pre-registered, decisions.md, commit bbb5d2f; Deviation 15)

Written by `research/waa_ksurface/a7_real_fit.py`. Real TRAIN data of this branch (28 and 38 GHz), the 14 learning events (rain = radar / g_e). VAL and TEST were not read. The procedure is that of step 3c itself (`s3c_a7_variant.procedure`, Deviation 15): the tail flag first; wetting K0 with k_wet ESTIMATED (a nuisance, R3 revised); drying K0 against K1 on L-dryhigh, margin 0.00291 dB with a paired 90 % interval that excludes zero, selected once with k_wet at its gauge-input estimate and refitted at its radar-input estimate; F and P measured parts fixed, extensions free (R1).

## 1. The flag on the real learning events (before any fit)

Events checked: 13 (E20150831_0455 is not checkable). **Removed: 0**.

| event_id | link_events | qualifying | waiting_or_rising | flagged |
|---|---|---|---|---|
| E20150617_0705 | 55 | 39 | 14 | False |
| E20150618_0530 | 67 | 1 | 1 | False |
| E20150622_0930 | 67 | 10 | 2 | False |
| E20150705_1910 | 67 | 15 | 4 | False |
| E20150707_1235 | 67 | 2 | 1 | False |
| E20150725_0440 | 67 | 3 | 2 | False |
| E20150727_2140 | 67 | 20 | 3 | False |
| E20150804_1830 | 67 | 41 | 1 | False |
| E20150810_1545 | 67 | 18 | 1 | False |
| E20150815_1305 | 67 | 14 | 3 | False |
| E20150825_0145 | 67 | 15 | 3 | False |
| E20150826_0520 | 67 | 51 | 0 | False |
| E20150826_1620 | 67 | 20 | 8 | False |

## 2. k_wet (nuisance estimates, no claim; R3 revised)

1 / k_wet: gauge input 1.74 min, radar input 10.92 min.

## 3. The drying rung and k_dry

Scored values (L-dryhigh): 9438 in 410 link-events. **Drying rung: K1** (lambda 0.0001); K1 gain over K0 +0.02013 dB, paired 90 % interval +0.01525 to +0.02419; margin 0.00291 dB.

1 / k_dry with its 5-95 % event-bootstrap interval (200 draws, seed 3900, the selected rung refitted with k_wet held at that estimate) in brackets; for K0 the value does not depend on W.

| k_wet estimate | 1/k_wet (min) | 1/k_dry at W = 0 dB (min) | 1/k_dry at W = 0.5 dB (min) | 1/k_dry at W = 1 dB (min) | 1/k_dry at W = 2 dB (min) | bootstrap draws converged |
|---|---|---|---|---|---|---|
| gauge | 1.74 | 81.5 (37.1 to 252.7) | 12.2 (8.7 to 17.7) | 5.4 (4.2 to 7.3) | 6.0 (4.4 to 9.3) | 200 of 200 |
| radar | 10.92 | 95.1 (45.3 to 278.9) | 17.1 (12.1 to 25.0) | 8.4 (6.6 to 10.9) | 8.0 (5.2 to 11.8) | 200 of 200 |

### k_dry over the two k_wet estimates (R3 revised)

| W (dB) | 1/k_dry, union of the two intervals (min) | point values at gauge / radar estimate |
|---|---|---|
| 0.0 | 37.1 to 278.9 | 81.5 / 95.1 |
| 0.5 | 8.7 to 25.0 | 12.2 / 17.1 |
| 1.0 | 4.2 to 10.9 | 5.4 / 8.4 |
| 2.0 | 4.4 to 11.8 | 6.0 / 8.0 |

## 4. Result under the pre-registered readings

Selected drying rung: K1 (one selection, valid at both k_wet estimates).

**non-exponential drying: the drying rate falls as the antenna dries; the mechanism (W-dependent rate or a second, slow store) is not claimed (R5); a spread of the baseline between links is excluded as its cause at the measured size (0.088 dB, step 3c (a) PASS, also at 0.15 dB).**

Stated with every A7 result: the flag cannot check E20150831_0455, so unseen water in it, if any, would go undetected.


## 5. Extension knots above the measured ranges (labeled EFFECTIVE, R1)

Band index: 0 = 28 GHz, 1 = 38 GHz. These are effective quantities for this rain input, not physical claims (`research/waa_ksurface/outputs/tables/a7_real_fit_extension.csv`).

| k_wet estimate | curve and band index | values at the knots (EFFECTIVE, R1) |
|---|---|---|
| gauge | F band 0 | 2.513, 2.889, 4.671 |
| gauge | P band 0 | 0.157, 0.474, 0.812, 1.309 |
| gauge | F band 1 | 2.694, 2.978, 5.393 |
| gauge | P band 1 | 0.358, 0.988, 1.402, 2.410 |
| radar | F band 0 | 2.912, 3.051, 6.051 |
| radar | P band 0 | 0.251, 0.512, 0.911, 1.659 |
| radar | F band 1 | 3.077, 3.810, 6.810 |
| radar | P band 1 | 0.497, 1.022, 1.458, 2.773 |

