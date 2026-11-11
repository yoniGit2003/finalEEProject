# E3 diagnostic: which ingredient of the fake world breaks the bin rule

Written by `research/waa_ksurface/e3_diagnostic.py` after the STOP of E3. No rain error in any run. 5 replicates per row. Nothing is patched with this table. (`research/waa_ksurface/outputs/tables/e3_diagnostic_summary.csv`, `research/waa_ksurface/outputs/tables/e3_diagnostic.csv`)

| rule | memory | noise | cells | cells_outside | I_err_median | I_err_max_abs | I_outside | s_outside | s_err_abs_max | passes |
|---|---|---|---|---|---|---|---|---|---|---|
| R1 | True | True | 14 | 7 | 0.099 | 0.166 | 3 | 5 | 0.096 | False |
| R1 | True | False | 14 | 1 | 0.077 | 0.124 | 0 | 1 | 0.1 | True |
| R1 | False | True | 14 | 4 | 0.032 | 0.255 | 1 | 3 | 0.131 | False |
| R1 | False | False | 14 | 1 | 0.005 | 0.218 | 1 | 0 | 0.136 | True |
| R2 | True | True | 12 | 4 | 0.07 | 0.129 | 0 | 4 | 0.065 | False |
| R2 | True | False | 12 | 5 | 0.037 | 0.135 | 0 | 5 | 0.074 | False |
| R2 | False | True | 12 | 4 | 0.033 | 0.077 | 0 | 4 | 0.045 | False |
| R2 | False | False | 12 | 5 | 0.0 | 0.08 | 0 | 5 | 0.059 | False |

Intercept error (dB) per cell, rule R1:

| band | class | I_err memory=on noise=on | I_err memory=on noise=off | I_err memory=off noise=on | I_err memory=off noise=off |
|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.099 | 0.065 | 0.031 | 0.005 |
| 28 GHz | 0.2-0.5 | 0.096 | 0.072 | 0.02 | -0.003 |
| 28 GHz | 0.5-1 | 0.152 | 0.112 | 0.033 | 0.005 |
| 28 GHz | 1-2 | 0.136 | 0.091 | 0.035 | -0.008 |
| 28 GHz | 2-5 | 0.166 | 0.117 | 0.093 | 0.046 |
| 28 GHz | 5-10 | 0.047 | 0.011 | 0.044 | 0.01 |
| 28 GHz | 10-20 | 0.076 | 0.029 | 0.066 | 0.019 |
| 38 GHz | 0.1-0.2 | 0.099 | 0.06 | 0.031 | -0.001 |
| 38 GHz | 0.2-0.5 | 0.089 | 0.067 | -0.003 | -0.014 |
| 38 GHz | 0.5-1 | 0.131 | 0.118 | 0.022 | 0.006 |
| 38 GHz | 1-2 | 0.097 | 0.082 | 0.01 | -0.013 |
| 38 GHz | 2-5 | 0.113 | 0.108 | 0.045 | 0.039 |
| 38 GHz | 5-10 | -0.026 | -0.022 | -0.026 | -0.03 |
| 38 GHz | 10-20 | 0.164 | 0.124 | 0.255 | 0.218 |
