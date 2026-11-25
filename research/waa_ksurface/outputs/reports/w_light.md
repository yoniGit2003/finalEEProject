# Light-rain overestimate: a diagnostic

Written by `research/waa_ksurface/w_light.py` (decisions.md, "Light-rain overestimate: diagnostic", logged before any code was run). Saved held-out predictions of the WAA-only comparison only (path A), TRAIN (14 learning events, 4563 link-events); no model fitted; VAL and TEST not opened. Residual = model - measured (5 min A - c), scored bins. Intervals: event bootstrap, 1000 draws, seed 4500, 90 %.

## 1. Residual length test (per link-event median residual = a + b L, robust line)

| model | class (mm/h) | link-events | a [90 %] | b [90 %] | reading |
|---|---|---|---|---|---|
| ours | 0.1-0.2 | 4236 | -0.049 [-0.156, +0.039] | +0.0187 [+0.0026, +0.0358] | path |
| ours | 0.2-0.5 | 4354 | -0.038 [-0.146, +0.040] | -0.0009 [-0.0199, +0.0188] | neither |
| ours | 0.5-1 | 4036 | -0.080 [-0.185, +0.006] | -0.0224 [-0.0472, +0.0028] | neither |
| ours | 1-2 | 3569 | +0.104 [-0.088, +0.272] | -0.0925 [-0.1404, -0.0478] | neither |
| ours | 2-5 | 3000 | +0.036 [-0.155, +0.227] | -0.0050 [-0.0674, +0.0625] | neither |
| Part 1 equation | 0.1-0.2 | 4236 | +0.016 [-0.085, +0.103] | +0.0185 [+0.0041, +0.0341] | path |
| Part 1 equation | 0.2-0.5 | 4354 | +0.040 [-0.052, +0.117] | +0.0016 [-0.0158, +0.0196] | neither |
| Part 1 equation | 0.5-1 | 4036 | +0.004 [-0.102, +0.092] | -0.0151 [-0.0415, +0.0114] | neither |
| Part 1 equation | 1-2 | 3569 | +0.121 [-0.068, +0.283] | -0.0874 [-0.1334, -0.0432] | neither |
| Part 1 equation | 2-5 | 3000 | +0.039 [-0.144, +0.230] | -0.0054 [-0.0660, +0.0650] | neither |
| Schleiss | 0.1-0.2 | 4236 | +0.048 [-0.042, +0.120] | +0.0198 [+0.0049, +0.0370] | path |
| Schleiss | 0.2-0.5 | 4354 | +0.083 [-0.039, +0.182] | +0.0055 [-0.0116, +0.0266] | neither |
| Schleiss | 0.5-1 | 4036 | +0.015 [-0.110, +0.147] | -0.0114 [-0.0419, +0.0175] | neither |
| Schleiss | 1-2 | 3569 | -0.072 [-0.268, +0.097] | -0.0783 [-0.1231, -0.0370] | neither |
| Schleiss | 2-5 | 3000 | -0.628 [-0.801, -0.429] | +0.0208 [-0.0318, +0.0812] | neither |
| Pastorek | 0.1-0.2 | 4236 | -0.001 [-0.130, +0.108] | +0.0221 [+0.0034, +0.0422] | path |
| Pastorek | 0.2-0.5 | 4354 | +0.022 [-0.086, +0.102] | +0.0083 [-0.0137, +0.0306] | neither |
| Pastorek | 0.5-1 | 4036 | +0.004 [-0.137, +0.110] | -0.0135 [-0.0397, +0.0144] | neither |
| Pastorek | 1-2 | 3569 | +0.066 [-0.174, +0.319] | -0.0787 [-0.1362, -0.0296] | neither |
| Pastorek | 2-5 | 3000 | -0.038 [-0.329, +0.259] | +0.0080 [-0.0531, +0.0668] | neither |
| no memory | 0.1-0.2 | 4236 | -0.084 [-0.209, +0.023] | +0.0298 [+0.0112, +0.0495] | path |
| no memory | 0.2-0.5 | 4354 | -0.049 [-0.156, +0.034] | +0.0015 [-0.0205, +0.0237] | neither |
| no memory | 0.5-1 | 4036 | -0.072 [-0.208, +0.033] | -0.0172 [-0.0435, +0.0105] | neither |
| no memory | 1-2 | 3569 | +0.269 [+0.038, +0.513] | -0.0916 [-0.1492, -0.0424] | antenna |
| no memory | 2-5 | 3000 | +0.111 [-0.156, +0.389] | +0.0193 [-0.0410, +0.0765] | neither |
| no WAA | 0.1-0.2 | 4236 | -0.452 [-0.581, -0.345] | +0.0298 [+0.0107, +0.0501] | path |
| no WAA | 0.2-0.5 | 4354 | -0.634 [-0.741, -0.552] | +0.0188 [-0.0033, +0.0410] | neither |
| no WAA | 0.5-1 | 4036 | -0.959 [-1.100, -0.853] | +0.0022 [-0.0244, +0.0298] | neither |
| no WAA | 1-2 | 3569 | -1.262 [-1.492, -1.019] | -0.0558 [-0.1108, -0.0091] | neither |
| no WAA | 2-5 | 3000 | -1.955 [-2.225, -1.677] | +0.0458 [-0.0111, +0.0994] | neither |

## 2. Tail split (light rain, 0.1-1 mm/h)

| group | minutes | ours residual | Part 1 residual | our WAA part | minutes with W < 0.3 dB |
|---|---|---|---|---|---|
| (i) within 60 min after R >= 2 | 270680 | -0.061 [-0.167, +0.042] | +0.025 [-0.068, +0.119] | +0.891 [+0.858, +0.937] | 10490 (4 % [3, 5]) |
| (ii) no R >= 2 in the previous 60 min | 737275 | -0.003 [-0.094, +0.083] | +0.070 [-0.005, +0.149] | +0.403 [+0.377, +0.433] | 204785 (28 % [24, 31]) |

(i) - (ii): ours -0.058 [-0.187, +0.060] dB; Part 1 -0.044 [-0.160, +0.066] dB.

## 3. Median residual per model and class (* = interval above 0)

| model | 0.1-0.2 | 0.2-0.5 | 0.5-1 | 1-2 | 2-5 |
|---|---|---|---|---|---|
| ours | +0.022 [-0.069, +0.107] | -0.022 [-0.118, +0.065] | -0.079 [-0.166, -0.009] | -0.067 [-0.259, +0.139] | +0.032 [-0.259, +0.357] |
| Part 1 equation | +0.082 [+0.006, +0.158] * | +0.061 [-0.011, +0.137] | +0.021 [-0.067, +0.090] | -0.026 [-0.219, +0.195] | +0.043 [-0.255, +0.389] |
| Schleiss | +0.118 [+0.055, +0.179] * | +0.123 [+0.047, +0.213] * | +0.048 [-0.060, +0.169] | -0.171 [-0.337, +0.037] | -0.562 [-0.830, -0.223] |
| Pastorek | +0.078 [-0.047, +0.186] | +0.042 [-0.043, +0.136] | +0.012 [-0.073, +0.079] | -0.085 [-0.347, +0.156] | -0.058 [-0.455, +0.355] |
| no memory | +0.013 [-0.104, +0.118] | -0.047 [-0.127, +0.047] | -0.074 [-0.158, -0.006] | +0.088 [-0.153, +0.313] | +0.120 [-0.254, +0.504] |
| no WAA | -0.354 [-0.473, -0.251] | -0.592 [-0.675, -0.504] | -0.917 [-0.996, -0.853] | -1.364 [-1.622, -1.136] | -1.897 [-2.252, -1.506] |

Reading per class: 0.1-0.2: mixed; 0.2-0.5: mixed; 0.5-1: mixed; 1-2: mixed; 2-5: mixed.

## 4. Gauge check (links < 1 km with a City gauge within 1 km; light-rain bins with a valid gauge 5 min rate)

18 links, 377 link-events, 17457 bins; radar-only (gauge < 0.1 mm/h): 13993 bins = 80.2 % [77.5, 82.4].

| model | radar-only | others | difference |
|---|---|---|---|
| ours | +0.076 [+0.039, +0.116] | -0.516 [-0.731, -0.369] | +0.592 [+0.455, +0.808] |
| Part 1 equation | +0.138 [+0.099, +0.178] | -0.454 [-0.637, -0.292] | +0.592 [+0.434, +0.780] |
| Schleiss | +0.149 [+0.117, +0.203] | -0.372 [-0.607, -0.129] | +0.521 [+0.300, +0.749] |
| Pastorek | +0.170 [+0.120, +0.220] | -0.620 [-0.815, -0.428] | +0.790 [+0.604, +0.984] |
| no memory | +0.093 [+0.046, +0.139] | -0.682 [-0.886, -0.506] | +0.775 [+0.605, +0.979] |
| no WAA | -0.471 [-0.527, -0.398] | -1.373 [-1.563, -1.197] | +0.902 [+0.730, +1.098] |

## 5. Figure

`research/waa_ksurface/outputs/figures/w_light_classes.png`

## Readings (the logged rules)

- Check 1 (ours, light classes 0.1-0.2: path, 0.2-0.5: neither, 0.5-1: neither; heavier: 1-2: neither, 2-5: neither).
- Check 2: our (i) - (ii) = -0.058 [-0.187, +0.060] dB -> not a tail; group (ii) residual of ours -0.003 [-0.094, +0.083] dB -> no overestimate without heavy rain before.
- Check 3: 0.1-0.2: mixed, 0.2-0.5: mixed, 0.5-1: mixed.
- Check 4: radar-only share 80 % [77, 82]; our difference +0.592 [+0.455, +0.808] dB -> input (the overestimate sits where the gauge sees no rain).

Run time 5.5 min.
