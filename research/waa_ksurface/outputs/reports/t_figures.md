# TEST comparison figures (a visualization only)

Written by `research/waa_ksurface/t_figures.py` (decisions.md, "TEST comparison figures", logged before the figures). Data: only the saved outputs of the TEST run (`predictions_test.pkl`); no model rerun, TEST not read again. No fitting, no metric, no decision follows.

## Selected link-events (rules written before any curve was seen)

| rule | candidates | median RMSE ours (dB) | event | sublink | L (km) | band | RMSE ours (dB) | rank (ascending) | min after last wet bin | figure |
|---|---|---|---|---|---|---|---|---|---|---|
| rule1_short | 563 | 0.428 | E20150828_0630 | 10031_1 | 0.876 | 38 GHz | 0.428 | 282 of 563 | 180 | research/waa_ksurface/outputs/figures/waa_compare_rule1_short_E20150828_0630_10031_1.png |
| rule2_long | 897 | 0.427 | E20150628_2200 | 10210_2 | 3.8877 | 28 GHz | 0.427 | 449 of 897 | 180 | research/waa_ksurface/outputs/figures/waa_compare_rule2_long_E20150628_2200_10210_2.png |
| rule4_new | 190 | 0.418 | E20150613_1100 | 10107_2 | 0.5154 | 28 GHz | 0.4172 | 95 of 190 | 180 | research/waa_ksurface/outputs/figures/waa_compare_rule4_new_E20150613_1100_10107_2.png |

**Rule 3 (strong storm) not applied:** the peak rain rate needs the radar series, which the TEST run did not save.

## Summary figure

`research/waa_ksurface/outputs/figures/waa_compare_summary_short_links.png`: 753 link-events of 53 links < 1 km (known and new) with at least 60 min after the last wet bin, aligned at the end of that bin.

## Not plotted

Not saved by the TEST run, so not shown: radar rain (bars and wet-minute shading), the WAA part of each model. The grey span in the per-event figures marks the first to the last wet 5 min bin of the link (data layer); dry bins inside it are not visible.

