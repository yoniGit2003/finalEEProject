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

