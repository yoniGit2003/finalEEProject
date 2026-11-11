# E2: which length test survives the measured rain error (fake data)

Written by `research/waa_ksurface/e2_candidates.py`. Fake world on the 4563 TRAIN link-events of 187 links in the 14 events with a gauge ratio. Bin rule R2 (30 min within a factor 1.5), n_short = 5. Truth: F = 3·R / (R + 2) dB, P = 0.1·R / 0.2·R dB/km, one store (5 / 15 min), measured noise. Rain error (E1): scale part per link-event sd 0.43, per-bin part with 46 % shared in space, lag-1 correlation 0.35, per-bin sd 0.48 (low end) and 0.87 (high end). 5 replicates, medians. Pass (Deviation 3): intercept within 0.15 dB and slope within 0.0710 dB/km in all judged cells but at most 1, at both ends. No real attenuation is used as a target.

## Verdict

| candidate | config | cells_judged | cells_inside | I_outside | s_outside | largest_I_err_db | largest_s_err_db_km | passes |
|---|---|---|---|---|---|---|---|---|
| (i) gauge-ratio events, plain | sigma_lo | 11 | 2 | 8 | 2 | 0.517 | 0.231 | False |
| (ii) instrument binning (area-mean radar) | sigma_lo | 8 | 6 | 1 | 2 | 0.173 | 0.156 | False |
| (iii) forward-bias correction | sigma_lo | 11 | 3 | 7 | 6 | 0.8 | 0.451 | False |
| (i) gauge-ratio events, plain | sigma_hi | 0 | 0 | 0 | 0 |  |  | False |
| (ii) instrument binning (area-mean radar) | sigma_hi | 0 | 0 | 0 | 0 |  |  | False |
| (iii) forward-bias correction | sigma_hi | 0 | 0 | 0 | 0 |  |  | False |

**Adopted: none (STOP).**

## The largest rain error each candidate still passes

| candidate | largest rain error passed |
|---|---|
| (i) gauge-ratio events, plain | 0.25 x the low end (per-bin sd 0.12, scale part 0.11) |
| (ii) instrument binning (area-mean radar) | none (it fails even without rain error) |
| (iii) forward-bias correction | 0.5 x the low end (per-bin sd 0.24, scale part 0.21) |

The sweep (the whole rain error scaled by f; f = 1 is the low end):

| candidate | config | sigma | scale_sd | cells_judged | cells_inside | largest_I_err_db | largest_s_err_db_km | passes |
|---|---|---|---|---|---|---|---|---|
| (i) gauge-ratio events, plain | f=0 | 0.0 | 0.0 | 12 | 12 | 0.129 | 0.065 | True |
| (ii) instrument binning (area-mean radar) | f=0 | 0.0 | 0.0 | 10 | 6 | 0.106 | 0.105 | False |
| (iii) forward-bias correction | f=0 | 0.0 | 0.0 | 12 | 12 | 0.129 | 0.065 | True |
| (i) gauge-ratio events, plain | f=0.25 | 0.12 | 0.108 | 12 | 11 | 0.15 | 0.069 | True |
| (ii) instrument binning (area-mean radar) | f=0.25 | 0.12 | 0.108 | 10 | 8 | 0.118 | 0.134 | False |
| (iii) forward-bias correction | f=0.25 | 0.12 | 0.108 | 12 | 12 | 0.137 | 0.056 | True |
| (i) gauge-ratio events, plain | f=0.5 | 0.24 | 0.215 | 12 | 8 | 0.185 | 0.039 | False |
| (ii) instrument binning (area-mean radar) | f=0.5 | 0.24 | 0.215 | 10 | 7 | 0.098 | 0.173 | False |
| (iii) forward-bias correction | f=0.5 | 0.24 | 0.215 | 12 | 11 | 0.164 | 0.061 | True |
| (i) gauge-ratio events, plain | f=0.75 | 0.36 | 0.322 | 12 | 4 | 0.243 | 0.091 | False |
| (ii) instrument binning (area-mean radar) | f=0.75 | 0.36 | 0.322 | 10 | 8 | 0.106 | 0.148 | False |
| (iii) forward-bias correction | f=0.75 | 0.36 | 0.322 | 12 | 6 | 0.364 | 0.1 | False |

## Every cell at both ends of the range

### (i) gauge-ratio events, plain, sigma_lo

| band | class | R_cell | F_true | I | I_err | P_true | s | s_err_abs | links_short | I_ok | s_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.2-0.5 | 0.326 | 0.42 | 0.753 | 0.328 | 0.033 | -0.004 | -0.038 | 5.0 | False | True |
| 28 GHz | 0.5-1 | 0.743 | 0.813 | 1.17 | 0.357 | 0.074 | 0.092 | 0.018 | 8.0 | False | True |
| 28 GHz | 1-2 | 1.442 | 1.257 | 1.681 | 0.424 | 0.144 | 0.156 | 0.022 | 10.0 | False | True |
| 28 GHz | 2-5 | 3.071 | 1.817 | 2.083 | 0.276 | 0.307 | 0.32 | 0.013 | 13.0 | False | True |
| 28 GHz | 5-10 | 6.314 | 2.278 | 2.304 | 0.025 | 0.631 | 0.556 | -0.067 | 5.0 | True | True |
| 38 GHz | 0.1-0.2 | 0.128 | 0.181 | 0.264 | 0.074 | 0.026 | 0.036 | 0.009 | 6.0 | True | True |
| 38 GHz | 0.2-0.5 | 0.366 | 0.464 | 0.797 | 0.302 | 0.073 | 0.06 | -0.02 | 8.0 | False | True |
| 38 GHz | 0.5-1 | 0.769 | 0.833 | 1.278 | 0.416 | 0.154 | 0.128 | -0.022 | 8.0 | False | True |
| 38 GHz | 1-2 | 1.487 | 1.279 | 1.71 | 0.517 | 0.297 | 0.24 | -0.063 | 11.0 | False | True |
| 38 GHz | 2-5 | 2.982 | 1.796 | 1.944 | 0.168 | 0.596 | 0.699 | 0.116 | 15.0 | False | False |
| 38 GHz | 5-10 | 6.378 | 2.284 | 2.397 | 0.123 | 1.276 | 0.935 | -0.231 | 10.0 | True | False |

### (i) gauge-ratio events, plain, sigma_hi

no cell judged

### (ii) instrument binning (area-mean radar), sigma_lo

| band | class | R_cell | F_true | I | I_err | P_true | s | s_err_abs | links_short | I_ok | s_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.1-0.2 | 0.137 | 0.192 | 0.18 | -0.029 | 0.014 | 0.01 | -0.005 | 17.0 | True | True |
| 28 GHz | 0.5-1 | 0.763 | 0.828 | 0.704 | -0.137 | 0.076 | 0.063 | -0.011 | 17.0 | True | True |
| 28 GHz | 1-2 | 1.567 | 1.317 | 1.276 | -0.041 | 0.157 | 0.142 | -0.015 | 17.0 | True | True |
| 28 GHz | 2-5 | 2.892 | 1.774 | 1.931 | 0.123 | 0.289 | 0.163 | -0.12 | 17.0 | True | False |
| 38 GHz | 0.1-0.2 | 0.102 | 0.146 | 0.106 | -0.04 | 0.02 | 0.026 | 0.005 | 24.0 | True | True |
| 38 GHz | 0.5-1 | 0.674 | 0.756 | 0.661 | -0.082 | 0.135 | 0.076 | -0.057 | 24.0 | True | True |
| 38 GHz | 1-2 | 1.469 | 1.27 | 1.259 | -0.029 | 0.294 | 0.358 | 0.069 | 24.0 | True | True |
| 38 GHz | 2-5 | 3.116 | 1.827 | 2.0 | 0.173 | 0.623 | 0.467 | -0.156 | 24.0 | False | False |

### (ii) instrument binning (area-mean radar), sigma_hi

no cell judged

### (iii) forward-bias correction, sigma_lo

| band | class | R_cell | F_true | I | I_err | P_true | s | s_err_abs | links_short | I_ok | s_ok |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 28 GHz | 0.2-0.5 | 0.326 | 0.42 | 0.684 | 0.273 | 0.033 | -0.034 | -0.066 | 5.0 | False | True |
| 28 GHz | 0.5-1 | 0.743 | 0.813 | 0.788 | -0.025 | 0.074 | 0.149 | 0.075 | 8.0 | True | False |
| 28 GHz | 1-2 | 1.442 | 1.257 | 1.411 | 0.141 | 0.144 | 0.138 | 0.004 | 10.0 | True | True |
| 28 GHz | 2-5 | 3.071 | 1.817 | 2.254 | 0.447 | 0.307 | 0.286 | -0.02 | 13.0 | False | True |
| 28 GHz | 5-10 | 6.314 | 2.278 | 2.599 | 0.313 | 0.631 | 0.525 | -0.115 | 5.0 | False | False |
| 38 GHz | 0.1-0.2 | 0.128 | 0.181 | 0.058 | -0.106 | 0.026 | 0.019 | -0.008 | 6.0 | True | True |
| 38 GHz | 0.2-0.5 | 0.366 | 0.464 | 0.513 | 0.049 | 0.073 | 0.069 | -0.01 | 8.0 | True | True |
| 38 GHz | 0.5-1 | 0.769 | 0.833 | 1.602 | 0.8 | 0.154 | 0.048 | -0.101 | 8.0 | False | False |
| 38 GHz | 1-2 | 1.487 | 1.279 | 0.968 | -0.212 | 0.297 | 0.171 | -0.126 | 11.0 | False | False |
| 38 GHz | 2-5 | 2.982 | 1.796 | 2.059 | 0.264 | 0.596 | 0.741 | 0.145 | 15.0 | False | False |
| 38 GHz | 5-10 | 6.378 | 2.284 | 2.546 | 0.228 | 1.276 | 0.869 | -0.451 | 10.0 | False | False |

### (iii) forward-bias correction, sigma_hi

no cell judged

(`research/waa_ksurface/outputs/tables/e2_candidates_fake.csv`, `research/waa_ksurface/outputs/tables/e2_candidates_verdict.csv`, `research/waa_ksurface/outputs/figures/e2_candidates_fake.png`)

## STOP

No candidate passes at both ends of the measured rain error. E4 is not run.

