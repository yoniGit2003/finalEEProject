# wet_antenna: the WAA overview restructured as a two-column paper (v2)

Written by hand on 2026-11-02. The previous overview is kept as `docs/WAA_research_overview_v1.tex` (and its PDF); the new one is
`docs/WAA_research_overview.tex` / `.pdf` (REVTeX 4-2, reprint, 12 pages including the appendices), modelled on
`docs/reference/Pulsed_NMR_Haber_Tsitseluk.pdf`. Figures: `docs/WAA_research_overview_figures/v2_*.png`, drawn for one column
(8.6 cm, fonts 8 pt) by `figs_v2.py` from saved tables only, colour-free (line styles, markers, labels). Nothing about our
equation changed. One new computation: `research/wet_antenna/comparison_all_models.py` (report
`comparison_all_models.md`, table `comparison_all_models.csv`).

## 1. The new outline

| Section | Content | Where it came from (v1) |
|---|---|---|
| Abstract | problem, approach (length + dry paths, no power law), the equation in one line, the test against all models, the main limit | new |
| 1 Introduction | CMLs, WAA and short links, why no power law is subtracted, Schleiss / Pastorek / ATTRRA2 / runoff, aim, four contributions | v1 1.1-1.3, the literature paragraph of v1 6 |
| 2 Data and preprocessing | OpenMRG, bands, split + funnel (Table I), how dry, the rain input g_e, the baseline in one paragraph | v1 2, 3 (baseline shortened), 4.4 (g_e) |
| 3 Method | decomposition + assumptions (Table II), the length test (short), the two-store model, free curves then closed forms, the validation protocol stated once | v1 4.1-4.3, 4.5, 5.1-5.2 |
| 4.1 | the model-free length test | v1 4.2 results |
| 4.2 | the equation in a boxed figure* (Fig. 4, the old Fig. 8 moved here, parameter table inside), F / P figure, fake recovery | v1 4.5, 5.1, 8 |
| 4.3 | W independent of L (three checks) + the VAL gap = rain input | v1 5.3, 4.4 |
| 4.4 | ALL models in one section: Table III accuracy, Table IV physical consistency (new, Part 1), Fig. 8 accuracy, Fig. 9 after-rain curves, the precise conclusion | v1 6 (comparison models and literature paragraph), new |
| 4.5 | P1-P5 (Table V), the light-rain level, the S8b adoption of the closed F | v1 6 |
| 5 Discussion | equal accuracy is not equivalence; the light-rain level; the end-of-rain level and the faster drying (box B); dew; the error budget; P at 0.52-0.74 of ITU-R and Part 2; the three boxes (Table VI) | v1 7, 4.5 |
| 6 Conclusions | one paragraph and the next steps (R1, Part 2) | v1 7.1 |
| App. A | the baseline (method and validation) | v1 3 |
| App. B | the dated decisions (Table VII) | `docs/ROADMAP.md` part B (new in the document) |
| App. C | side roads, one paragraph each (tails, M5, per-event factors, M4c, wind, g(R), paired links) | `docs/research_log.tex` (new in the document) |
| App. D | uncertainty and identifiability methods | v1 5.1, the equation file, `tighten_s6.md` |
| App. E | full TEST tables: RMSE by length class (Table VIII), the eight length-test cells for every model (Table IX) | `literature_s8c.md`, `comparison_all_models.md` |

## 2. What was cut or shortened, and where it went

- The baseline (v1 Section 3, a page with a numbered list and a dew subsection): two sentences in 2.4, the details in Appendix A.
  Dew: one paragraph in the Discussion.
- The illustrative three-link example of the length test (v1 Table 3 and its paragraph): cut; 3.2 states the idea in one paragraph.
- "Evidence for the results" (v1 Section 5: fake data, unseen events, independence) and the summary table of results (v1
  Table 5): the protocol is stated once in 3.5, the fake recovery in 4.2, the independence in 4.3; the summary table is replaced
  by the three boxes (Table VI) and the Conclusions.
- The comparison models of the test (v1 Section 6) and the literature paragraph (v1 Section 6, "Comparison with published WAA
  models", Table 7, Figure 7): merged into one section, 4.4, with one protocol and one accuracy table.
- The parameter table of v1 Section 4.5 (Table 4) and the box of v1 Section 8: merged into the one boxed Figure 4 (the box keeps
  the parameter table; the separate Table 4 is gone). The box no longer repeats the TEST RMSE (it is in Table III).
- Implementation details (pycomlink code, the Schleiss changes): one sentence each in 4.4; the full description stays in
  `literature_s8c.md`.
- The step labels (S1-S9, M4b-B, B0-B2) are not used in the paper text; the models are named by what they are.
- Lists were turned into paragraphs, bold and semicolons removed.

## 3. Part 1: the consistency numbers (comparison_all_models)

Computed on TEST with the frozen fits of S8 / S8b / S8c, nothing refitted (the S8 P1 cells of the table F and the S8c dry-path
errors are reproduced exactly). Summary (full table: `research/wet_antenna/outputs/tables/comparison_all_models.csv`):

| Model | slope inside (of 8) | intercept inside | dry-path error 5/10/20 min (dB) | rho(a_link, L) | wet / dry 1/e time (min) | fitted global parameters |
|---|---|---|---|---|---|---|
| no WAA | 0 | 0 | +0.35 / +0.27 / +0.15 | - | - | 3 |
| constant WAA | 4 | 0 | +0.35 / +0.27 / +0.15 | -0.60 [-0.66, -0.52] | 0 / 0 | 6 |
| no memory | 7 | 6 | +0.35 / +0.27 / +0.15 | +0.20 [+0.08, +0.32] | 0 / 0 | 12 |
| Schleiss (a) | 5 | 2 | -0.05 / +0.10 / +0.13 | - | 4.5 / 4.5 | 3 |
| Schleiss (b) | 6 | 7 | -0.24 / -0.25 / -0.18 | -0.45 [-0.55, -0.35] | 103 / 25 | 8 |
| Pastorek (a) | 4 | 7 | +0.33 / +0.25 / +0.15 | - | 0 / 0 | 3 |
| Pastorek (b) | 7 | 6 | +0.33 / +0.25 / +0.15 | +0.31 [+0.18, +0.42] | 0 / 0 | 7 |
| ours | 7 | 4 | -0.19 / -0.15 / -0.11 | +0.09 [-0.03, +0.20] | 2.6 / 10.3 (gauge) | 7 (+ a_link) |

**The requested conclusion does not hold as worded, and the paper states what the numbers support.** The request was to state that
ours "is the only model whose path loss, after-rain water and timing agree with the independent measurements". The numbers show:

- Path loss: ours agrees in 7 of 8 cells, but so do the model without memory and the refitted Pastorek model (7 of 8).
- After-rain water: no model agrees within its interval. Ours holds 0.11-0.19 dB too much, Schleiss (b) 0.18-0.25 dB too much,
  the models without memory 0.15-0.35 dB too little. Schleiss (a) has a smaller error at 5 and 10 min.
- Timing: ours agrees with the gauges because it was fixed from the gauge fit; this is not independent evidence.
- Intercept (WAA level): ours is the weakest of the fitted WAA models (4 of 8; Schleiss (b) 7, Pastorek (b) 6, no memory 6).

What the paper states instead (Section 4.4): ours is better than or equal to every model in accuracy and ties only with the
refitted Schleiss model; of the models with a link factor it is the only one whose link factor does not depend on the length
(rho +0.09 [-0.03, +0.20]; Schleiss (b) -0.45, Pastorek (b) +0.31, constant -0.60, no memory +0.20); its slope agrees in 7 of 8
cells; its after-rain error is the smallest of all fitted models except the default Schleiss at 5 and 10 min; its timing agrees by
construction; its intercept is the weakest, the failed rule P1.

**Dated addition.** The a_link-against-L column was added on 2026-11-02 AFTER the first run of the script (its other numbers had
been seen), because Pastorek (b) has c_band 0.6-0.7 x ours but the measured slope. It is the S7 check 7.1 applied to every model,
written into the config (`comparison_all_models.a_link_vs_L`) with that date and reason, and said in the paper (4.4, Appendix B).
No pass rule was written for the consistency table; its numbers are reported as they came out.

## 4. The number check (S9 rule: the repo wins)

Every statement with a number in the new file was checked against the file that covers it. Numbers carried over from v1 were
checked in S9 (`formulate_s9_check.md`, row numbers in brackets) and are confirmed again from the same sources; new numbers and
changes are listed individually. Paths: `R/` = `research/wet_antenna/outputs/reports/`, `T/` = `research/wet_antenna/outputs/tables/`.

| Where | Statement | Verdict | Source |
|---|---|---|---|
| Abstract, 4.2, box | equation: 0.72, 2.6, 10.3, 240 min, b 0.54, R/5 | confirmed | `products/waa_equation_final.json` |
| Abstract, 4.5 | 18 held-out storms | confirmed | `T/test_s8_funnel.csv` |
| Abstract, 5 | light-rain level 0.4 to 0.6 dB above the final equation | new | `T/comparison_all_models_cells.csv`, ours at 0.5-2 mm/h: 1.354-0.80, 1.433-1.05, 1.437-0.87, 1.554-1.09 |
| Abstract, 4.4, 6 | only model with a link factor independent of L | new | `T/comparison_all_models.csv` (a_link_vs_L) |
| 1 | 0.36 km and 4.13 km link (Fig. 1) | confirmed (6) | `T/formulate_s9_example_event.csv` |
| 2.1 | 1 June-31 August 2015, 728 sub-links, 364 links, 10 s, 0.3 / 1 dB, 7.5-38.7 GHz, 5 min, 2 km, ten gauges, 1 min | confirmed (9-15) | `pipeline/docs/openmrg_summary.md`, `T/tighten_s6_bands.csv` |
| 2.1 | 660 of 728 (91 %), the other 68 on links of 4 km or more | confirmed (17, 18) | `T/tighten_s6_funnel.csv`, `research/wet_antenna/METHOD.md` |
| 2.2 | 71 of 364 links (20 %) held out | confirmed (19) | release split tables |
| 2.2 | offset spread at most 0.3 dB, the selection steps, largest loss = no rain on the link (13 480 of 29 998) | new | `T/tighten_s6_funnel.csv` |
| 2.2 | 12 % of TRAIN minutes >= 0.5 mm/h, 0.4 % >= 10 mm/h | confirmed | `T/tighten_s6_minutes.csv` (1 134 690 and 33 815 of 9 499 740) |
| Table I | 53 (51 with link-events) / 18 / 18 events; 9 780 / 3 250 / 3 355 + 1 010; 206 TRAIN links; 59 new links | confirmed (22-24) | `T/tighten_s6_funnel.csv`, `T/test_s8_funnel.csv`, `products/waa_equation_final.json` |
| 2.3 | g_e only if >= 1 mm per gauge link; radar 0.89 x gauges on TRAIN, 1.71 x on VAL | confirmed (v1: about 0.9 and 1.7) | `R/conditions_s4c.md` (G_TRAIN 0.894, G_VAL 1.707) |
| 2.4, App. A | 24 h (48 h), 1 h after rain, median 0.00, -0.21 / +0.22 dB, 0.3 dB exclusion, +0.04 dB after 1 h, zero within 12 h | confirmed | `research/baseline/outputs/reports` (as v1, S9 rows on Section 3) |
| 2.4, box | c over the dry minutes 120 to 30 min before the first rain | confirmed (38) | `products/waa_equation_final.json` |
| Table II | quadratic term in 4 of 15 classes | confirmed | as v1, S9 check |
| 3.2 | rain classes 0.5-1 ... 10-20 mm/h, 1000 draws, >= 8 links under 1 km | confirmed | config `system_id_s3.length_test`, `T/test_s8_length_cells.csv` |
| 3.3 | 5 min input makes the wetting 2.1 x slower | confirmed | as v1 |
| 3.4 | knots 0.05-40 mm/h, Huber 0.3 dB | confirmed | config `system_id_s2.model`, `literature_s8c.fit` |
| 4.1 | intercept about 1.0 dB at 0.5-1 mm/h, 2.4-2.8 dB at 5-20 mm/h; short-link path loss halved when links to 3 km were added | confirmed (47-50) | `T/system_id_s3_length_fits.csv`, `R/system_id_s2.md` |
| 4.2 | RMSE 0.3875 against 0.3878 dB | confirmed (60) | `R/tighten_s6.md` |
| 4.2, box | b 0.54 [0.47, 0.63], A 2.57 / 2.36 / 2.77 with intervals, c 0.105 / 0.141 / 0.194 with intervals, a_link 0.77-1.24, 206 links, 40 refits, 2.6 [2.4, 3.1], 10.3 [9.2, 11.8], VAL 2.0 [1.5, 3.2], 8.8 [5.8, 13.1], phi 0.72 [0.60, 0.81], tau_s 240, d 1 | confirmed | `products/waa_equation_final.json`, `R/final_s5a.md` |
| box | bands 25-30.5 / 30.5-35 / 35-45 GHz; validity 0-40 mm/h, 0.4 %, 0.14-8.7 km; 1e-15 dB | confirmed | `products/waa_equation_final.json`, `R/formulate_s9.md` (9e-16 dB) |
| 4.2 | 9 of 9 within 20 %; 62 numbers within 20 %, 61 within 10 %; 8 path numbers below 1 mm/h not identified | confirmed (75) | `R/test_s8b.md`, `R/tighten_s6.md` |
| 4.3 | rho +0.05 [-0.07, +0.17] (frozen model) | confirmed | `R/validate_s7.md` 7.1 |
| 4.3 | rho +0.09 [-0.03, +0.20] (final equation) | new | `T/comparison_all_models.csv` |
| 4.3 | 13 of 15 cells; 9 081 windows on 206 links in 23 TRAIN and 10 VAL events; -0.18 dB at 10 min in every class | confirmed; "23 TRAIN and 10 VAL" new | `R/validate_s7.md`, `T/validate_s7_dry_windows.csv` |
| 4.3 | VAL length test at most 1 of 8 (raw radar), 5 of 6 (with g_e) | confirmed | `R/final_s5a.md`, `docs/ROADMAP.md` |
| 4.4 | Pastorek d = 0.1; pycomlink 0.6.0; wet when R > 0.1 mm/h; Schleiss (a) 15 min | confirmed | `R/literature_s8c.md`, config `literature_s8c` |
| Table III | every RMSE and paired difference, all and L < 1 km; ours 0.695 [0.514, 0.803] and 0.564 [0.436, 0.633] | confirmed; the L < 1 km intervals are new in the document | `T/literature_s8c_rmse.csv`, `T/literature_s8c_paired.csv` |
| Table IV | every entry | new | `T/comparison_all_models.csv`, `R/comparison_all_models.md` |
| 4.4 | Pastorek (b) path loss 0.6-0.7 x ours | new | `products/literature_s8c_fits.json` (0.065 / 0.098 / 0.125 against 0.105 / 0.141 / 0.194: 0.62, 0.70, 0.64) |
| 4.4 | Schleiss (b) path loss 1.7-1.8 x ours | confirmed | `products/literature_s8c_fits.json` (1.82, 1.72, 1.73) |
| 4.4, 5 | Schleiss (b) wetting about 100 min, drying about 25 min | **corrected** | v1: "wetting time 309 min, dries in 77 min" (the nominal tau of the pycomlink step 3 dt / tau); new: the 1/e times of that step, 102.6 and 25.1 min (`T/comparison_all_models.csv`); the nominal values are in the Table IV caption's method |
| 4.4 | after-rain errors 0.11-0.19 (ours), 0.18-0.25 (Schleiss b), 0.15-0.35 (no memory) dB; 1791 windows; measured 0.35 / 0.27 / 0.15 dB | new | `T/comparison_all_models_dry_path.csv` |
| 4.5, Table V | 3 355 / 1 010 / 59; one ratio of zero and three missing; P1 4 of 8, slopes 7, intercepts 4; P2 numbers; P3 -0.04 [-0.15, +0.04]; P4 1.12, 5 %; P5 1.07, 1.20 | confirmed (95-100) | `R/test_s8.md` |
| 4.5 | 1.35-1.55 against 0.83-1.11 dB | confirmed (98) | `T/test_s8_length_cells.csv` |
| 4.5 | +0.0025 [-0.0065, +0.0085] dB, +-0.01 dB rule | confirmed | `R/test_s8b.md` |
| 5 | -0.25 dB at the end of rain, 0.5 dB after heavy rain, 0.87 [0.70, 1.06], 3.7 min, wind -0.76 over 17 storms | confirmed (104-108) | `R/return_s6b.md`, `R/validate_s7.md` |
| 5 | 125 link-nights, 1.7 dB, about 4 h, rho +0.26 | confirmed | `R/validate_s7.md` 7.4 |
| 5 | TEST budget: 61 % within storms, 8 % between, 6 % after rain, 25 % dry | confirmed; 8 % and 25 % new | `T/test_s8_error_budget.csv` |
| 5 | VAL budget: 65 % within storms, 2 % after rain | **corrected** | v1: "61 to 64 %" (TEST and VAL); `T/return_s6b_error_budget.csv` gives 0.6454, which rounds to 65 % |
| 5 | P exponent 0.97, quadratic 0; ITU-R exponent 0.86-0.92; size 0.52-0.74 x ITU-R | confirmed (67-69) | `R/final_s5a.md`, `T/validate_s7_p_size.csv`, `R/validate_s7.md` 7.3e |
| App. B | dates 7, 13, 17, 18, 19, 28 Oct, 2 Nov; 0.0005 dB; 4 of 6 against 0.6667 | confirmed | `docs/ROADMAP.md` part B, config `literature_s8c.fit.outer_schleiss`, `comparison_all_models.a_link_vs_L` |
| App. C | 20 355 tails, 429 TRAIN tails, one in five flat or rising, closed 13 Oct, 4 pairs | confirmed | `docs/research_log.tex` (T1-T18), `R/validate_s7.md` |
| App. D | 40 cold refits, profile 1-13 x narrower, three fake datasets within 20 %, 1000 TEST draws | confirmed | `products/waa_equation_final.json`, `R/tighten_s6.md`, config `literature_s8c.evaluation` |
| Table VIII | RMSE by length class | confirmed, new in the document | `T/literature_s8c_rmse.csv` |
| Table IX | every cell value | new | `T/comparison_all_models_cells.csv` |
| References | Messer 2006, Leijnse 2008, Schleiss 2013, Pastorek 2022, Tiede 2023, Andersson 2022, Cranmer 2023, Brunton 2016 as in the proposal; Tiede 2025 as v1 (checked in S9, row 113) | confirmed | `docs/proposal.md` |
| References | ITU-R P.838-3 | new citation | named in `pipeline/docs/openmrg_summary.md` (the processing recipe); title and year only, no further details claimed |

**Changes, in short.** (1) The Schleiss (b) times are now the 1/e times of its discrete step (about 100 and 25 min), not the
nominal 309 and 77 min. (2) The VAL error share within storms is 65 % (v1 "61 to 64 %" for both sets). (3) Our TEST RMSE interval
is given once, from the paired S8c draws (0.695 [0.514, 0.803]); v1 also gave [0.497, 0.795] from the S8b draws in its box. (4) New
numbers: the consistency table, the a_link-against-L correlations, the L < 1 km intervals, the per-class RMSE, the cell table, the
light-rain gap of the final equation (0.4-0.6 dB), the full TEST error budget. (5) The abstract and the conclusion of 4.4 are
worded to what the numbers support (Section 3 above).

## 5. Open items

- The conclusion asked for ("the only model whose path loss, after-rain water and timing agree with the independent
  measurements") is not supported (Section 3); the learners should read 4.4 and decide whether they accept the precise wording.
- The a_link-against-L column was added after seeing the other numbers (dated); it is information, not a rule.
- The literature models have no fake-data recovery and no intervals; their comparison is post-freeze.
- `docs/research_log.tex`, `docs/project_structure.tex` and `docs/REPRODUCE.md` (regenerated) are not otherwise updated for this
  change; run /updatedocs when this version is accepted.
- The old figures (`figs.py`, used by v1) and the new ones (`figs_v2.py`) live side by side in `docs/WAA_research_overview_figures/`.
