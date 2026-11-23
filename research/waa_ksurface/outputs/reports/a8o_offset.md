# Offset check: plain exponential drying, measured offsets, the 0.3 dB grid

Written by `research/waa_ksurface/a8o_offset.py` (decisions.md, "Offset check", logged before the run). Fake data only; VAL and TEST were not opened.

Truth K0 / K0: 1 / k_dry 9.7 min, 1 / k_wet 10.92 min; the step 3c fake world with the normal offset replaced by offsets drawn from the 4547 measured learning-set offsets (sd 0.136 dB, mean +0.015 dB, 10-90 % -0.115 to +0.149 dB, 1-99 % -0.313 to +0.442 dB). Full A7 procedure (step 3c); then the step 8 form k = c2/(c1 - exp(c0 W)) plugged back at the radar-input estimate of k_wet.

| seed name | seed | drying rung | K1 gain on L-dryhigh (dB) | interval | 1/k_dry gauge / radar (min) | 1/k_wet gauge / radar (min) | events flagged | offsets drawn mean / sd (dB) | switch c0, c1, c2 | W_half (dB) | 1/k_inf (min) | switch converged | minutes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| r0 | 4200 | K0 | -0.00145 | -0.00425 to +0.00071 | 13.28 / 13.69 | 4.02 / 4.68 | 0 of 13 | +0.013 / 0.130 | -52.15, -0.08313, -0.006123 | 0.048 | 13.58 | True | 6.5 |
| r1 | 4201 | K0 | 0.00134 | +0.00092 to +0.00223 | 8.71 / 10.35 | 4.14 / 8.58 | 0 of 13 | +0.013 / 0.131 | -53.89, -0.03895, -0.003822 | 0.06 | 10.19 | False | 6.3 |
| r2 | 4202 | K0 | -0.00048 | -0.00072 to -0.00019 | 10.69 / 11.66 | 4.11 / 6.12 | 0 of 13 | +0.014 / 0.135 | -31.37, -0.2277, -0.01975 | 0.047 | 11.53 | False | 6.8 |
| r3 | 4203 | K0 | 0.00096 | +0.00063 to +0.00174 | 10.23 / 11.62 | 4.33 / 7.34 | 0 of 13 | +0.019 / 0.134 | -80.21, -0.01439, -0.00125 | 0.053 | 11.51 | True | 10.1 |
| r4 | 4204 | K0 | -0.00038 | -0.00112 to +0.00150 | 10.68 / 11.28 | 4.11 / 5.30 | 0 of 13 | +0.016 / 0.138 | -71.96, -0.03347, -0.002992 | 0.047 | 11.18 | True | 784.9 |

**Reading: the non-exponential claim stands.** K1 in 0 of 5, K0 in 5 of 5; median W_half over the K1 seeds: n/a (floor if K1 >= 3 and median W_half in [0.2, 0.4]; stands if K0 >= 4).

## The naming cross-check (run 2026-11-22 11:53-12:19; report `reports/a8x_crosscheck.md`)

Sparse regression (best subset of up to 3 terms over {1, W, W^2, 1/W, sqrt(W), exp(-W), W/(1+W)}) and the Bayesian Machine Scientist (no fallback needed), scored by the step 8 objective and rule.

Chosen: **sparse: 1 + sqrt(W) + exp(-W)** (p = 3); AICc gauge -133.83, radar -163.17; k = c0*1 + c1*sqrt(W) + c2*exp(-W)  with  c0 = 0.4182, c1 = -0.1342, c2 = -0.4003 (gauge); k = c0*1 + c1*sqrt(W) + c2*exp(-W)  with  c0 = 0.2481, c1 = -0.06193, c2 = -0.2401 (radar).

MDL models: gauge `(W * (((_a1_ + W) ** _a0_) * W))`, radar `(_a3_ * (W ** (_a1_ ** W)))`.

## Agreement and the plug-back queue

Step 8's form: **pysr_gauge_c8** (accepted); AICc gauge -137.34, radar -207.05, p = 3.

| pair | same family | dAICc gauge | dAICc radar | agree |
|---|---|---|---|---|
| sparse vs step 8 | False | 3.51 | 43.88 | False |
| BMS vs step 8 | False | 11.14 | 31.41 | False |
| sparse vs BMS | False | -7.63 | 12.48 | False |

Queued for plug-back (p <= 3 and AICc lower by more than 2 at both estimates): **0**. They run after step 8 has finished and at least 4 GB is free.


