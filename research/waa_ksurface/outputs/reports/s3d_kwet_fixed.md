# Step 3d: the A7 procedure with 1 / k_wet fixed at 2.5, 5, 10 min (fake data)

Written by `research/waa_ksurface/s3d_kwet_fixed.py` (procedure of `s3c_a7_variant.py`). Fake attenuation as in step 3c, fresh seeds (base 3840). The tail flag before any fit; the 14 learning events only; stage 1 (gauge input, wetting K0) once, then stage 2 (drying K0 against K1 on L-dryhigh, margin 0.00291 dB, paired 90 % interval excluding zero) with the wetting rate HELD at 1 / k_wet = 2.5, 5, 10 min. Truth (a): constant k, 1 / k_dry = 15 min; truth (c): k_dry = (1 + W) / 30 (ratio 3). PASS = (a) K0 with 1 / k_dry within 25 % and (c) K1 at every value; a rung flip is a STOP. No real attenuation is read. VAL and TEST were not opened.

## Every fit

| replicate | seed | 1/k_wet (min) | true rung | drying rung | K1 gain (dB) | interval | 1/k_dry fit (min) | 1/k_dry truth (min) | share within 25 % | 1/k_dry at W = 0.5, 1, 2 dB | verdict |
|---|---|---|---|---|---|---|---|---|---|---|---|
| a_r0 | 3840 | 2.5 | K0 | K1 | 0.00328 | +0.00256 to +0.00509 | 18.2 | 14.9999999925 | 1.0 | 16.4, 14.7, 13.2 | FLIP |
| a_r0 | 3840 | 5.0 | K0 | K1 | 0.00296 | +0.00226 to +0.00478 | 19.3 | 14.9999999925 | 1.0 | 17.9, 16.6, 13.6 | FLIP |
| a_r0 | 3840 | 10.0 | K0 | K1 | 0.00566 | +0.00443 to +0.00870 | 23.4 | 14.9999999925 | 0.08 | 21.2, 19.1, 12.3 | FLIP |

**Verdict: STOP (rung flip)** (only 1 of 4 replicates).


**STOP by the rule: (a_r0) rung flip.** The 'rung not determined' reading would be expected on real data; the review decides before A7.


## Execution

| replicate | events flagged | minutes | cpu minutes | peak memory (GB) |
|---|---|---|---|---|
| a_r0 | 0 of 13 | 13.6 | 50.0 | 2.85 |

Keep-awake request accepted by Windows: True; gaps of the memory watch above 3 min: 0; lowest free memory seen 0.40 GB; low-memory events: 2.

