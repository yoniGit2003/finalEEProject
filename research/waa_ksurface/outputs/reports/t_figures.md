# TEST comparison figures (a visualization only)

Written by `research/waa_ksurface/t_figures.py` (decisions.md, "TEST comparison figures replaced"). Measured A - c and predictions: only the saved outputs of the TEST run; no model rerun. Radar rain input R = radar / g_e of the two link-events reloaded from the TEST rain input (no attenuation read). No fitting, no metric, no decision follows.

## The two link-events (the rules of the first figures, unchanged)

| rule | candidates | median RMSE ours (dB) | event | sublink | L (km) | band | RMSE ours (dB) | rank (ascending) | min after last wet bin | g_e | peak R (mm/h) | figure |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| rule1_short | 563 | 0.428 | E20150828_0630 | 10031_1 | 0.876 | 38 GHz | 0.428 | 282 of 563 | 180 | 1.0 | 0.6249 | research/waa_ksurface/outputs/figures/waa_compare_rule1_short_E20150828_0630_10031_1.png |
| rule2_long | 897 | 0.427 | E20150628_2200 | 10210_2 | 3.8877 | 28 GHz | 0.427 | 449 of 897 | 180 | 0.4524 | 10.7803 | research/waa_ksurface/outputs/figures/waa_compare_rule2_long_E20150628_2200_10210_2.png |

## Not plotted

Not saved by the TEST run, so not shown: the WAA part of each model (only the total prediction was saved).

