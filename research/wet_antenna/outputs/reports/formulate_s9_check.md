# wet_antenna step S9: the overview checked against the repo

Written by hand in step S9 (2026-10-19), with the rule of config `wet_antenna.formulate_s9.check`: every statement of the overview
`docs/WAA_research_overview.pdf` (13 pages, 19 October 2026) is checked against the files that cover it. The repo wins.

**The source.** The overview's LaTeX source and figure script were not in the repo; only the compiled PDF was. Its source was
therefore rebuilt from the PDF (same structure and wording, now `docs/WAA_research_overview.tex`) and its figures were redrawn from
saved tables (`docs/WAA_research_overview_figures/figs.py`). The corrected document replaces the first PDF.

Verdicts: **confirmed** (with its source), **corrected** (old, new, source), **unsupported** (no file backs it; then removed or
reworded). Paths: `R/` = `research/wet_antenna/outputs/reports/`, `T/` = `research/wet_antenna/outputs/tables/`.

## The table

| # | Section | Statement in the overview | Verdict | Source, and the correction |
|---|---|---|---|---|
| 1 | Front | authors, CellEnMon lab, 19 October 2026 | confirmed | `CLAUDE.md`, `docs/proposal.md` (the lab); the date of the document |
| 2 | 1.1 | A(t) = TSL - RSL, both in dBm, A in dB | confirmed | `docs/proposal.md` (A = TSL - RSL - baseline: the baseline is removed in Section 3 of the overview) |
| 3 | 1.1 | rain on the path makes A rise; networks can be used to measure rain | confirmed | `docs/proposal.md` (Messer et al. 2006, proposal ref. [5]) |
| 4 | 1.1 | A_rain = integral of a R^b dl approx a R-bar^b L; a, b depend on frequency | confirmed | `docs/proposal.md`, Eq. of section 2 ("frequency and polarization") |
| 5 | 1.2 | rain wets the radomes; WAA is local, does not grow with L, has memory | confirmed | `docs/proposal.md` section 2; `CLAUDE.md` (physics) |
| 6 | 1.2 | on a 0.5 km link WAA is most of what is measured, on 5 km a smaller share in heavy rain, all that is left after rain | confirmed | `T/formulate_s9_example_event.csv` (a 0.36 km and a 4.13 km link in one TRAIN storm); the illustration of Figure 1 is replaced by this real-data figure |
| 7 | 1.2 | if WAA is ignored, rain is overestimated, most strongly on short links | confirmed | `docs/proposal.md` section 2 (short links deviate; WAA mistaken for path loss); S8 B0 (no WAA) puts 2.5 x more loss on the path: `T/test_s8_comparison_train_val.csv` |
| 8 | 1.3 | Part 1 finds W(t); Part 2 the length dependence; WAA never computed by subtracting a power law | confirmed | `docs/proposal.md`; `docs/ROADMAP.md` (standing rule) |
| 9 | 2 | OpenMRG, Gothenburg, 1 June - 31 August 2015, ref. [1] | confirmed | `pipeline/docs/openmrg_summary.md`; Andersson et al. 2022 (proposal ref. [13]) |
| 10 | 2 | 728 sub-links, 364 links, both directions, 10 s | confirmed | `pipeline/docs/openmrg_summary.md` |
| 11 | 2 | RSL 0.3 dB, TSL 1 dB | confirmed | `pipeline/docs/openmrg_summary.md` |
| 12 | 2 | frequencies 7 to 39 GHz | confirmed | `T/tighten_s6_bands.csv` (7.46 to 38.68 GHz); the paper says 7-38 GHz |
| 13 | 2 | radar every 5 min, about 2 x 2 km, sees every path | confirmed | `pipeline/docs/openmrg_summary.md` |
| 14 | 2 | radar values off by a different factor in each storm | confirmed | `R/conditions_s4c.md` (g_e 0.9 on TRAIN, 1.7 on VAL) |
| 15 | 2 | ten City gauges, accurate at their point | confirmed | `pipeline/docs/openmrg_summary.md` (1 min, 0.1 / 0.2 mm) |
| 16 | 2 | data cleaned, synchronised and averaged to 5 min values in a pipeline | corrected | old: "averaged to 5 minute values"; new: "to 1 minute and 5 minute values" (the model runs on 1 min attenuation and is judged on 5 min means): `research/wet_antenna/system_id.py` build_dataset, `pipeline/docs/` |
| 17 | 2 | kept 28, 32, 38 GHz: 660 of 728 sub-links (91 %) | confirmed | `T/tighten_s6_funnel.csv` |
| 18 | 2 | the lower bands have too few links to study on their own | corrected | old: "too few links"; new: "the 68 sub-links below 25 GHz all belong to links of 4 km or more, so their antenna part cannot be separated by length": `research/wet_antenna/METHOD.md` (T18 inventory: all links below 25 GHz are >= 4 km), `T/tighten_s6_bands.csv` |
| 19 | 2 | whole events go to one of three sets, 20 % of links set aside | confirmed | `data/release` split tables: 71 of 364 links TEST (19.5 %), 8 excluded |
| 20 | 2 | fitted on TRAIN, VAL compares choices | confirmed | `docs/ROADMAP.md` integrity rule 6 |
| 21 | 2 | TEST has not been opened yet | corrected | old: "not opened yet, read once at the very end"; new: "read once, in S8, after the freeze": `R/test_s8.md` |
| 22 | 2 | Table 1: TRAIN 53 events, 9780 link-events | corrected | 53 TRAIN events exist, 51 of them hold link-events: `T/tighten_s6_funnel.csv`, `T/tighten_s6_minutes.csv`; the table now gives 53 (51 with link-events) |
| 23 | 2 | Table 1: VAL 18 events, 3250 link-events | confirmed | `T/tighten_s6_funnel.csv` |
| 24 | 2 | Table 1: TEST 18 events, not read | corrected | new: 18 events, 3355 link-events on known links and 1010 on new links: `T/test_s8_funnel.csv` |
| 25 | 2 | 12 % of the minutes >= 0.5 mm/h, 0.4 % >= 10 mm/h | confirmed | `R/tighten_s6.md` section 2 |
| 26 | 3 | baseline: one flat number per event and sub-link | confirmed | `research/baseline/METHOD.md` (v2) |
| 27 | 3.1 | the 24 h before the event, 48 h if too little | confirmed | `research/config.yaml` baseline.v2 (window 24 h, fallback 48 h if < 1 h of dry slots) |
| 28 | 3.1 | only slots marked dry by radar and gauges, >= 1 h after the last rain | confirmed | `research/config.yaml` baseline.skip_h = 1.0; the link reference state uses radar and gauges |
| 29 | 3.1 | the median of A over these slots | confirmed | `research/baseline/METHOD.md`; the median of TL (A = TL - baseline): wording kept |
| 30 | 3.1 | links with a strong daily swing or erratic behaviour are flagged and left out | confirmed | `research/baseline/make_baseline_table.py` use_in_waa_main (daily_cycle_large, excursion, degraded, toggling) |
| 31 | 3.2 | tested on fake events in dry periods where the truth is zero | confirmed | `research/baseline/outputs/reports/R1b_validation.md` |
| 32 | 3.2 | median error 0.00 dB, 80 % within about +-0.22 dB | confirmed | `R1b_validation.md`: median -0.00, p10 / p90 -0.21 / +0.22 dB (6 h) |
| 33 | 3.2 | each sub-link keeps its error, large errors not used | confirmed | `research/config.yaml` system_id_s2.data.max_offset_std_db = 0.3 |
| 34 | 3.2 | after real rain the level stays above a dry control and returns within hours | confirmed | `R1b_validation.md` (+0.05 dB at 1 h, 0 by 12 h against the control) |
| 35 | 3.3 | on humid nights A rises without rain (dew); such nights are flagged | confirmed | `research/baseline/METHOD.md` (dew_risk flags), `R/validate_s7.md` 7.4 |
| 36 | 3.3 | dew later used as an independent check | confirmed | `R/validate_s7.md` 7.4 |
| 37 | 4.1 | A - c = W + L P(R-bar), Eq. 3 | confirmed | `research/wet_antenna/system_id.py` kernel2, `final_model.py` |
| 38 | 4.1 | c is a small leftover offset measured just before the rain | corrected | made precise: c = the median of A over the dry minutes 120 to 30 min before the link's first wet bin (>= 30 minutes, else 0): `system_id.py` build_dataset, `research/config.yaml` offset_window_min [120, 30] |
| 39 | 4.1 | P is not a fixed number but an unknown curve learned from the data | confirmed | `R/tighten_s6.md` (free knots, then a closed form) |
| 40 | 4.1 | Table 2: the two parts add up, checked by the length test | confirmed | `R/system_id_s3.md`, `R/system_id_s3b.md` |
| 41 | 4.1 | Table 2: W does not depend on L, three checks | confirmed | `R/validate_s7.md` 7.1, 7.3a, 7.3c |
| 42 | 4.1 | Table 2: path part L times a loss per km, no bend found in A against L | corrected | old: "no bend found"; new: "a quadratic term differs from zero in 4 of 15 rain groups only": `R/system_id_s3.md` (the quadratic check) |
| 43 | 4.1 | Table 2: radar along the path corrected by the gauges closes the TRAIN / VAL gap | confirmed | `R/conditions_s4c.md` |
| 44 | 4.1 | Table 2: rain drives the wetting, dew outside the model | confirmed | `docs/ROADMAP.md` (three boxes: dew in B) |
| 45 | 4.1 | L P(R-bar) exact when P is a straight line | confirmed | algebra (integral of c R dl = c L R-bar) |
| 46 | 4.2 | the length test: the 3-link example | confirmed | an illustration (arithmetic) |
| 47 | 4.2 | pool all 5 min measurements of a rain group, fit a straight line | corrected | old: "all 5 minute measurements"; new: "for every link and event, the median of A - c over its 5 min values in steady rain of that group (rain in the group for 20 minutes, the previous 15 min mean within x0.5 to x2); a robust (Huber) line through these points, intervals by resampling events": `research/config.yaml` system_id_s3.length_test, `research/wet_antenna/length_test.py` |
| 48 | 4.2 | Figure 2 illustrated with points generated by the model | corrected | replaced by the real TRAIN data of S3: `T/system_id_s3_length_fits.csv`, `research/wet_antenna/outputs/products/system_id_s3_steady_link_events.csv` |
| 49 | 4.2 | the antenna part rises from about 1 dB in drizzle to about 2.5 to 3 dB, then slowly; path part roughly in proportion to R | confirmed | `R/system_id_s3.md` (I 1.02 to 2.82 dB at 28 GHz) |
| 50 | 4.2 | short links alone cannot separate; the final model is fitted on all lengths | confirmed | `R/system_id_s2b.md`, `R/system_id_s3b.md` |
| 51 | 4.3 | Eq. 5-7: the fast and slow stores, phi 0.72, 2.6 / 10.3 min, 240 min | confirmed | `research/config.yaml` final_s5a.fits M4b-B fixed; `waa_equation_final.json` |
| 52 | 4.3 | Eq. 8 and 63 % / 37 % | confirmed | algebra |
| 53 | 4.3 | the two times measured with the gauges, the same on TRAIN and VAL | corrected | old: "came out the same"; new: "came out similar on VAL (2.0 and 8.8 min against 2.6 and 10.3 min)": `docs/ROADMAP.md` (findings after S3), `R/system_id_s3.md` |
| 54 | 4.3 | Figure 3: 60 min of 5 mm/h on a typical 28 GHz link | confirmed | kept as an illustration from the equation (no real figure shows the stores separately): `T/formulate_s9_memory.csv` (`waa_equation.py`); the caption says it is an illustration |
| 55 | 4.4 | radar where and when, gauges how much; one number per event | confirmed | `R/conditions_s4c.md`; the ratio is used if the gauges caught >= 1 mm per gauge link, else 1 (`waa_equation_final.json` rain_input) |
| 56 | 4.4 | in the large VAL storms the radar read about 1.7 times too much | confirmed | `R/conditions_s4.md` (1.71 on VAL, 0.89 on TRAIN) |
| 57 | 4.4 | after the correction TRAIN and VAL agreed | confirmed | `R/conditions_s4c.md` (VAL length test 1 of 8 -> 5 of 6) |
| 58 | 4.4 | the radar error does not depend on the rain rate (a statistical illusion) | confirmed | `R/validate_s7.md` 7.6 (the two binnings move apart: regression to the mean) |
| 59 | 4.5 | F and P first fitted as free curves | confirmed | `R/final_s5a.md`, `R/tighten_s6.md` |
| 60 | 4.5 | the three bands share one shape of F and differ only in height | corrected | new: "one shape describes the bands as well on unseen events (VAL 0.3875 against 0.3878 dB), although the training events alone prefer separate shapes": `R/tighten_s6.md` section 3 |
| 61 | 4.5 | F rises quickly then slowly, P a straight line | confirmed | `R/tighten_s6.md` |
| 62 | 4.5 | physical forms and PySR forms judged by the same rule of fit against simplicity (AICc) | confirmed | `R/tighten_s6.md` section 6 |
| 63 | 4.5 | the simplest form is a power close to a square root, proposed by PySR, put back and refitted | confirmed | `R/tighten_s6.md` section 6, `R/test_s8b.md` |
| 64 | 4.5 | b = 0.54, 0.47 to 0.63, 0.5 inside, exponent 1 ruled out | confirmed | `R/test_s8b.md` (0.542 [0.472, 0.626]) |
| 65 | 4.5 | A_band is the level at 5 mm/h; 90 % intervals by refitting on resampled whole events | confirmed | `R/test_s8b.md` (40 cold draws) |
| 66 | 4.5 | no maximum; few data above 10 mm/h; not to be used above 40 mm/h | confirmed | `R/test_s8b.md` section 6 |
| 67 | 4.5 | ITU exponent roughly 0.85 to 1 at 28 to 38 GHz | corrected | new: "about 0.86 to 0.92" (computed from the ITU-R P.838 values in `T/validate_s7_p_size.csv`, vertical, the bands' median frequencies) |
| 68 | 4.5 | a power law on the points gave 0.97, a quadratic term zero | confirmed | `R/tighten_s6.md` section 6 |
| 69 | 4.5 | the size is about 0.6 of the ITU value | confirmed | `R/validate_s7.md` 7.3e (0.52 to 0.74) |
| 70 | 4.5 | Table 3: A, c, b with intervals | confirmed | `R/test_s8b.md`, `waa_equation_final.json` |
| 71 | 4.5 | Table 3: a_link 80 % between 0.76 and 1.23 | corrected | new: 0.77 to 1.24 (the closed-F model; 0.76 to 1.23 was the free-F model of S6): `R/test_s8b.md` |
| 72 | 4.5 | Table 3: all checked on fake data | confirmed | `R/test_s8b.md` (9 of 9) |
| 73 | 4.5 | Figure 4: free points with intervals and the formula | confirmed | redrawn: `T/tighten_s6_bootstrap.csv`, `T/formulate_s9_closed_f_band.csv`; F and P on linear axes |
| 74 | 5 | a result only if it passed checks written before the outcome | confirmed | `docs/ROADMAP.md` integrity rules 1-3 |
| 75 | 5.1 | all 62 numbers came back within 10 % | corrected | old: "within 10 %"; new: "within 20 % in all three fake sets, 61 of 62 within 10 % (largest 10.3 %); the 9 numbers of the final equation within 20 % (largest 13 %)": `T/tighten_s6_identified.csv`, `T/test_s8b_closed_f_train.csv` |
| 76 | 5.1 | eight numbers (path loss below 1 mm/h) too close to zero | confirmed | `R/tighten_s6.md` section 5 |
| 77 | 5.2 | every model choice judged by its error on left-out events, never by the fit to training data | corrected | new: "model structures were judged on events left out of the fit; the choice among closed formulas used a fit-against-simplicity rule on the learned points, and the chosen formula was then checked on unseen events": `R/system_id_ladder.md` (leave-one-event-out judge), `R/tighten_s6.md` (AICc), `R/test_s8b.md` |
| 78 | 5.2 | a more complex version kept only if better on unseen events and its numbers measurable | confirmed | `docs/ROADMAP.md` (M4c dropped, M1 to box B) |
| 79 | 5.3 | a_link vs L: +0.05 (-0.07 to +0.17) | confirmed | `R/validate_s7.md` 7.1 |
| 80 | 5.3 | length test with L < 3 km and with L > 1.5 km: 13 of 15 | confirmed | `R/validate_s7.md` 7.3c |
| 81 | 5.3 | 9081 dry-path minutes after rain; the same for short, medium, long links | confirmed | `R/validate_s7.md` 7.3a (windows of minutes, not single minutes: wording "9081 such moments") |
| 82 | 5.3 | Figure 5: the multiplier against length | confirmed | redrawn with the final equation's a_link: `T/validate_s7_links.csv` |
| 83 | 5.3 | Figure 6: dry-path curves | confirmed | redrawn with the final equation: `T/formulate_s9_dry_curves.csv` |
| 84 | 5.4 | wets in about 3 min, dries in about 10 min, a small part for hours | confirmed | Eq. 6-7 |
| 85 | 5.4 | models with memory beat models without | confirmed | `R/system_id_ladder.md`, `R/test_s8.md` P2 (B2) |
| 86 | 5.4 | about 1 dB in drizzle, about 2.5 dB at 5 mm/h | confirmed | `T/validate_s7_w_height.csv`, `R/test_s8b.md` (A 2.36 to 2.77) |
| 87 | 5.4 | no band difference detected, a difference up to about 25 % cannot be excluded | unsupported | no band-ratio interval was computed; reworded to "the intervals overlap widely (for example 2.00 to 3.21 against 2.25 to 3.40 dB)": `R/test_s8b.md` |
| 88 | 5.4 | 80 % of links within 0.76 to 1.23 | corrected | 0.77 to 1.24 (see 71) |
| 89 | 5.4 | the straight line still passing the length test | confirmed | `R/tighten_s6.md` (plug-back P only, 5 of 6) |
| 90 | 6 | equation, rules, three simpler models fixed and fingerprinted before TEST | confirmed | `R/test_s8.md` section 1, `research/wet_antenna/outputs/products/test_s8_freeze.json` |
| 91 | 6 | 18 storms, 3355 link-events on known links, 1010 on 59 new links | confirmed | `R/test_s8.md` section 3 |
| 92 | 6 | the simpler models fitted on TRAIN the same way; no WAA takes a larger path loss | confirmed | `T/test_s8_comparison_train_val.csv` (c 0.27 against 0.105 dB/km per mm/h) |
| 93 | 6 | Table 5: error 0.69 against 0.86, 0.76, 0.73 (short 0.56 against 0.84, 0.64, 0.60) | confirmed | `R/test_s8.md` P2 |
| 94 | 6 | WAA the same on short and long links: error after rain equal | confirmed | `R/test_s8.md` P3 |
| 95 | 6 | no length class out of line: every class within 5 % of the others | corrected | new: "each class within 5 % of the overall TEST / VAL ratio" (classes 1.04, 1.03, 0.95 times it): `R/test_s8.md` P4 |
| 96 | 6 | links never seen: only 7 % worse | confirmed | `R/test_s8.md` P5 (1.07) |
| 97 | 6 | length test on TEST: 4 of 8, 6 needed | confirmed | `R/test_s8.md` P1 |
| 98 | 6 | in light rain about 1.4 to 1.5 dB against 0.8 to 1.1 dB | corrected | precise: 1.35 to 1.55 dB against 0.83 to 1.11 dB: `T/test_s8_length_cells.csv` |
| 99 | 6 | the formula run after the test without refitting, under a rule written first; 0.0025 dB (-0.0065 to +0.0085) | confirmed | `R/test_s8b.md` |
| 100 | 6 | (missing) two evaluation errors fixed after the TEST read | added | the full version reports them; the short version names them in one sentence: `R/test_s8.md` section 1, config `test_s8.results` |
| 101 | 7 | three groups: in the equation, seen outside, not explained | confirmed | `docs/ROADMAP.md` (three boxes) |
| 102 | 7 | dew wets the same antennas, about 1.7 dB for hours | confirmed | `R/validate_s7.md` 7.4 (rho +0.26, median 1.7 dB, about 4 h) |
| 103 | 7 | not explained: the link multiplier, the light-rain level, rain estimates not tested, F above 10 mm/h, P below 1 mm/h | confirmed | `R/validate_s7.md`, `R/test_s8.md`, `docs/ROADMAP.md` |
| 104 | 7 | the model holds about 0.2 dB more water right after rain, 0.5 dB after heavy rain | confirmed | `R/validate_s7.md` 7.3a, 7.3d |
| 105 | 7 | the antenna dries at the model's speed; not radar timing; not too high in heavy rain; larger after heavier rain and in drier air | confirmed | `R/return_s6b.md` Part 1 |
| 106 | 7 | drops run off the radome after heavy rain, also seen in laboratory tests [2] | corrected | the lab paper reports runoff right after wetting and water as drops (`docs/ROADMAP.md`, recorded literature), not specifically after heavy rain; reworded to "runoff of water right after wetting has been seen in laboratory tests [2]" |
| 107 | 7 | windy storms seemed to dry faster, but not once the preceding rain is accounted for | confirmed | `R/return_s6b.md` (partial Spearman) |
| 108 | 7 | 3.7 min drying removes the misfit, not adopted, path loss halves, does not match the gauge timing | confirmed | `R/return_s6b.md`, `docs/ROADMAP.md` (decision 2026-10-19) |
| 109 | 7 | the misfit is small in the error budget | confirmed | `T/return_s6b_error_budget.csv` (2 % on VAL), `T/test_s8_error_budget.csv` (6 % on TEST) |
| 110 | 7 | most remaining error comes from the radar seeing rain more coarsely than a link | unsupported | an interpretation, not tested; reworded to "most remaining error is during rain, within storms (61 to 64 %), where we think the coarser radar (5 min, 2 km) is the main cause" |
| 111 | 7.1 | next steps | corrected | new order from `docs/ROADMAP.md` (2026-10-19): R1 (does the correction improve link rain against the gauges), then Part 2 |
| 112 | Refs | [1] Andersson et al. 2022, ESSD 14, 5411-5426 | confirmed | `docs/proposal.md` ref. [13], `docs/papers/Andersson2022_OpenMRG.pdf` |
| 113 | Refs | [2] J. Tiede et al., "Wetness-induced variations in the antenna reflection coefficient", Adv. Radio Sci. 23, 13-19, 2025 | corrected | the title is wrong: "Frequency-dependent variations of the antenna reflection coefficient due to different wetness conditions on the antenna radome", J. Tiede, C. Chwala, U. Siart and T. Eibert, Adv. Radio Sci. 23, 13-19, 2025, doi:10.5194/ars-23-13-2025 (`docs/ROADMAP.md`, recorded and checked on the publisher's page) |

## Every change, in short

1. Data averaged to 1 and 5 minute values (16).
2. The lower bands: all their links are 4 km or longer (18).
3. TEST read once in S8, 3355 + 1010 link-events; TRAIN 53 events, 51 with link-events (21, 22, 24).
4. c defined precisely (38).
5. "No bend" made precise (42).
6. The length test uses per-link medians in steady rain and a robust line (47); its figure now shows real data (48).
7. The VAL timing is similar, not the same (53).
8. One shape for the bands: chosen on unseen events, not preferred by the training events (60).
9. The ITU exponent 0.86 to 0.92 (67).
10. a_link 0.77 to 1.24 (71, 88).
11. The fake-data recovery precise (75).
12. How model choices were judged (77).
13. The band difference of 25 % removed as unsupported (87).
14. P4 compares each class with the overall ratio (95).
15. The light-rain numbers made precise (98).
16. The two evaluation errors named (100).
17. The runoff sentence matched to what the lab paper reports (106); the reference title corrected (113).
18. The cause of the remaining error marked as our interpretation (110).
19. The next steps: R1, then Part 2 (111).
20. Figures 1, 2, 4, 5, 6 redrawn from saved real data; Figure 3 kept as an illustration from the equation (6, 48, 54, 73, 82, 83).
