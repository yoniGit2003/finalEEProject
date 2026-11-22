# Naming cross-check: sparse regression and the Bayesian Machine Scientist

Written by `research/waa_ksurface/a8x_crosscheck.py` (decisions.md, "Naming cross-check", logged before the run; disclosure there: the step 8 ranking had been seen, its plug-back and bootstrap not). Input and scoring as step 8 (1)-(2). VAL and TEST were not opened.

## Tool 1: sparse regression (best subset, 63 subsets)

Chosen: **sparse: 1 + sqrt(W) + exp(-W)** (p = 3); AICc gauge -133.83, radar -163.17; k = c0*1 + c1*sqrt(W) + c2*exp(-W)  with  c0 = 0.4182, c1 = -0.1342, c2 = -0.4003 (gauge); k = c0*1 + c1*sqrt(W) + c2*exp(-W)  with  c0 = 0.2481, c1 = -0.06193, c2 = -0.2401 (radar).

The 10 lowest Delta_max:

| name | p | aicc_gauge | aicc_radar | daicc_gauge | daicc_radar | dmax | chosen | formula_radar |
|---|---|---|---|---|---|---|---|---|
| sparse: 1 + sqrt(W) + exp(-W) | 3 | -133.827 | -163.17 | 0.0 | 0.0 | 0.0 | True | k = c0*1 + c1*sqrt(W) + c2*exp(-W)  with  c0 = 0.2481, c1 = -0.06193, c2 = -0.2401 |
| sparse: 1 + W + exp(-W) | 3 | -128.487 | -157.103 | 5.34 | 6.067 | 6.067 | False | k = c0*1 + c1*W + c2*exp(-W)  with  c0 = 0.1758, c1 = -0.0124, c2 = -0.1817 |
| sparse: 1/W + sqrt(W) + W/(1+W) | 3 | -131.16 | -150.307 | 2.667 | 12.863 | 12.863 | False | k = c0*1/W + c1*sqrt(W) + c2*W/(1+W)  with  c0 = 0.0005512, c1 = -0.08835, c2 = 0.3687 |
| sparse: sqrt(W) + W/(1+W) | 2 | -118.843 | -149.227 | 14.984 | 13.943 | 14.984 | False | k = c0*sqrt(W) + c1*W/(1+W)  with  c0 = -0.05803, c1 = 0.3062 |
| sparse: 1 + sqrt(W) + W/(1+W) | 3 | -128.519 | -147.697 | 5.308 | 15.473 | 15.473 | False | k = c0*1 + c1*sqrt(W) + c2*W/(1+W)  with  c0 = 0.002388, c1 = -0.07301, c2 = 0.3327 |
| sparse: W^2 + sqrt(W) + W/(1+W) | 3 | -123.693 | -147.632 | 10.134 | 15.538 | 15.538 | False | k = c0*W^2 + c1*sqrt(W) + c2*W/(1+W)  with  c0 = -0.0003595, c1 = -0.0592, c2 = 0.3101 |
| sparse: sqrt(W) + exp(-W) + W/(1+W) | 3 | -128.242 | -147.535 | 5.585 | 15.635 | 15.635 | False | k = c0*sqrt(W) + c1*exp(-W) + c2*W/(1+W)  with  c0 = -0.07212, c1 = 0.002177, c2 = 0.334 |
| sparse: W + sqrt(W) + W/(1+W) | 3 | -126.934 | -147.302 | 6.894 | 15.867 | 15.867 | False | k = c0*W + c1*sqrt(W) + c2*W/(1+W)  with  c0 = -0.002022, c1 = -0.06047, c2 = 0.3161 |
| sparse: 1 + W^2 + exp(-W) | 3 | -116.731 | -152.672 | 17.096 | 10.498 | 17.096 | False | k = c0*1 + c1*W^2 + c2*exp(-W)  with  c0 = 0.1557, c1 = -0.00141, c2 = -0.1605 |
| sparse: W + exp(-W) + W/(1+W) | 3 | -124.309 | -145.95 | 9.518 | 17.22 | 17.22 | False | k = c0*W + c1*exp(-W) + c2*W/(1+W)  with  c0 = -0.01407, c1 = -0.01155, c2 = 0.2266 |

## Tool 2: the Bayesian Machine Scientist (bitbucket.org/rguimera/machine-scientist, commit f4cf1e8)

Prior `nv1.np5.dat` (nv = 1, np = 5), operators {'+': 2, '*': 2, '/': 2, '**': 2, 'exp': 1, '-': 1}, max_size 10, seed 0, 3000 steps after 1000 burn-in. Probability = share of the samples after burn-in (BMS's own unweighted likelihood on k). Each formula refitted and scored by the step 8 objective.

| bms_estimate | rank | probability | bms_expr | p | aicc_gauge | aicc_radar | dmax | chosen |
|---|---|---|---|---|---|---|---|---|
| gauge | 1 | 0.4637 | (W * (((_a1_ + W) ** _a0_) * W)) | 2 | -129.5851 | -107.9186 | 67.7273 | False |
| gauge | 2 | 0.135 | (W / (_a0_ * ((_a2_ + W) ** _a2_))) | 2 | -78.2545 | -104.7992 | 70.8467 | False |
| gauge | 3 | 0.1337 | (W / (_a0_ + ((_a3_ * W) ** _a3_))) | 2 | -108.6457 | -160.7169 | 20.9395 | False |
| gauge | 4 | 0.066 | (((_a1_ * W) ** _a0_) ** (_a3_ ** W)) | 3 | -68.0588 | -96.6885 | 78.9574 | False |
| gauge | 5 | 0.0537 | _a2_ | 1 | -3.9015 | -6.7959 | 168.8501 | False |
| radar | 1 | 0.1577 | (_a2_ * (W ** ((_a2_ * _a1_) ** W))) | 2 | -126.1969 | -175.646 | 3.3882 | False |
| radar | 2 | 0.1357 | (_a1_ * (W ** ((_a2_ * _a1_) ** W))) | 2 | -126.1969 | -175.646 | 3.3882 | False |
| radar | 3 | 0.1007 | (_a2_ * (W ** ((_a2_ + _a1_) ** W))) | 2 | -126.1969 | -175.646 | 3.3882 | False |
| radar | 4 | 0.0757 | (_a2_ * (W ** ((_a2_ / _a1_) ** W))) | 2 | -126.1969 | -175.646 | 3.3882 | True |
| radar | 5 | 0.0403 | (_a3_ * (W ** (_a1_ ** W))) | 2 | -126.1969 | -175.646 | 3.3882 | False |

MDL models: gauge `(W * (((_a1_ + W) ** _a0_) * W))`, radar `(_a3_ * (W ** (_a1_ ** W)))`.

## Agreement and the plug-back queue

Step 8's form: **pysr_gauge_c8** (accepted); AICc gauge -137.34, radar -207.05, p = 3.

| pair | same family | dAICc gauge | dAICc radar | agree |
|---|---|---|---|---|
| sparse vs step 8 | False | 3.51 | 43.88 | False |
| BMS vs step 8 | False | 11.14 | 31.41 | False |
| sparse vs BMS | False | -7.63 | 12.48 | False |

Queued for plug-back (p <= 3 and AICc lower by more than 2 at both estimates): **0**. They run after step 8 has finished and at least 4 GB is free.

