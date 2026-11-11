# A3d diagnostic: the fake check against the size of the rain error

Written by `research/waa_ksurface/a3_eiv_diagnostic.py` after the STOP of A3d. The fake check of A3d with the rain error scaled (σ = scale × 0.708 / √n_pix), 5 replicates per scale. Pass limits: 0.15 dB and 20 %. Nothing is patched with this table. (`research/waa_ksurface/outputs/tables/a3_eiv_diagnostic.csv`)

| scale | sigma_1 | largest_I_err_db | largest_s_err_pct | cells_I_ok | cells_s_ok | cells |
|---|---|---|---|---|---|---|
| 0.0 | 0.0 | 0.111 | 59.038 | 14 | 13 | 14 |
| 0.25 | 0.177 | 0.161 | 33.96 | 13 | 13 | 14 |
| 0.5 | 0.354 | 0.659 | 39.871 | 10 | 11 | 14 |
| 1.0 | 0.708 | 1.229 | 64.93 | 0 | 8 | 14 |

Intercept error (dB) per cell:

| band | class | I_err at x0 | I_err at x0.25 | I_err at x0.5 | I_err at x1 |
|---|---|---|---|---|---|
| 28 GHz | 0.1-0.25 | 0.071 | -0.011 | 0.067 | 0.202 |
| 28 GHz | 0.25-0.5 | -0.009 | 0.009 | 0.053 | 0.248 |
| 28 GHz | 0.5-1 | -0.0 | 0.019 | 0.083 | 0.27 |
| 28 GHz | 1-2 | 0.005 | 0.019 | 0.077 | 0.287 |
| 28 GHz | 2-5 | 0.03 | 0.067 | 0.108 | 0.201 |
| 28 GHz | 5-10 | 0.009 | -0.032 | -0.166 | -0.366 |
| 28 GHz | 10-20 | 0.111 | -0.161 | -0.619 | -1.215 |
| 38 GHz | 0.1-0.25 | 0.005 | 0.036 | 0.075 | 0.291 |
| 38 GHz | 0.25-0.5 | 0.007 | 0.035 | 0.07 | 0.266 |
| 38 GHz | 0.5-1 | 0.016 | 0.028 | 0.083 | 0.322 |
| 38 GHz | 1-2 | 0.043 | 0.059 | 0.132 | 0.317 |
| 38 GHz | 2-5 | 0.066 | 0.104 | 0.195 | 0.221 |
| 38 GHz | 5-10 | 0.096 | -0.021 | -0.112 | -0.438 |
| 38 GHz | 10-20 | 0.059 | -0.059 | -0.659 | -1.229 |
