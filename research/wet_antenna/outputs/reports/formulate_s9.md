# wet_antenna step S9: formulate Part 1 (the WAA equation)

Step S9 (2026-10-19), settings `wet_antenna.formulate_s9` (written before the run). No new analysis, no refitting, no new reading
of TEST. **Part 1 is closed**: [FINAL.md](../../FINAL.md).

## 1. The official equation

`research/wet_antenna/outputs/products/waa_equation_final.json`, sha256 `b06ebc8b75867171368e8ea0540bb011014e9615d3dde74f8559b524a25be9cd`. It is the closed-F model of S8b (the S6 plug-back fit "F and P", fitted on TRAIN
before the S8 freeze, `research/.cache/wet_antenna/tighten_s6/plug_shared_FP_F_and_P.pkl`, sha256 `1d8163cd8b84aafc85f59c6c223f44b33923278265e3e40819923d8076c8e0b9`).

```
A(t) - c = W_f + W_s + L c_band R(t - 1 min)
dW_f/dt = (0.72 u - W_f) / tau_f,  tau_f = 2.6 min (wetting), 10.3 min (drying)
dW_s/dt = (0.28 u - W_s) / 240 min
u = a_link A_band (R / 5)^b
```

| number (90 % interval, 40 cold event-bootstrap draws) | 28 GHz | 32 GHz | 38 GHz |
|---|---|---|---|
| A_band (dB) | 2.57 [2, 3.21] | 2.36 [1.7, 3.06] | 2.77 [2.25, 3.4] |
| c_band (dB/km per mm/h) | 0.105 [0.0576, 0.138] | 0.141 [0.0474, 0.19] | 0.194 [0.111, 0.245] |

b = 0.542 [0.472, 0.626]. a_link 10 % / 90 %: 0.768 [0.716, 0.814] / 1.24 [1.29, 1.4] (per link in the file, 1 for a new link).
Valid 0-40 mm/h. Not identified: phi and tau_s (fixed, the slow store is not identified: S2b, S5a). the path term below 1 mm/h (S6). the level of F apart from a_link (the ridge prior sets the median a_link to about 1: S6).

## 2. The code and its unit test

`research/wet_antenna/waa_equation.py` (`WAAEquation`) takes a rain series (mm/h, one value per minute), L, the band and a_link,
and returns W_f, W_s, W and the path term. It steps the stores as the fit does (`system_id.kernel2`) and evaluates F on the fit's
61 knots by default. The reference `research/wet_antenna/outputs/products/waa_equation_reference.json` holds 3 TEST link-events of the S8 cache (already read once; the first
of each length class with >= 30 min of V1 >= 1 mm/h) with the fit's own 1 min W and path term:

| event | sub-link | L (km) | band | a_link | minutes |
|---|---|---|---|---|---|
| E20150601_1235 | 10040_1 | 3.65 | 28 GHz | 1.080 | 785 |
| E20150601_1235 | 10126_1 | 1.69 | 28 GHz | 0.970 | 770 |
| E20150602_0335 | 10001_1 | 0.69 | 28 GHz | 0.978 | 2220 |

`waa_equation.py` reproduces them within 8.9e-16 dB (W) and 4.4e-16 dB (path term), against the tolerance
1e-06 dB. Test: `python -m pytest research/wet_antenna/tests`.

## 3. The overview checked against the repo

`research/wet_antenna/outputs/reports/formulate_s9_check.md` (written by hand). The overview `docs/WAA_research_overview.pdf` had no LaTeX source in the repo, so
`docs/WAA_research_overview.tex` was rebuilt from the PDF. Verdicts of its 113 statements: confirmed 89, corrected 21, unsupported 2, added 1. The repo wins.
The changes, in short:

1. The pipeline averages to 1 and 5 minute values. The lower bands were dropped because all their links are 4 km or longer.
2. TEST is read (S8), with 3355 + 1010 link-events. TRAIN holds 53 events, 51 of them with link-events.
3. c is defined precisely. "No bend of A against L" became "a quadratic term in 4 of 15 rain groups".
4. The length test uses per-link medians in steady rain and a robust line, and its figure now shows the real TRAIN data.
5. The VAL timing is similar (2.0 / 8.8 min), not the same as TRAIN (2.6 / 10.3 min).
6. One shape for the three bands was chosen on unseen events. The training events alone prefer separate shapes.
7. The ITU exponent is 0.86-0.92. a_link is 0.77-1.24 for the final equation.
8. The fake recovery: 62 of 62 within 20 %, 61 within 10 %, and the 9 final numbers within 20 %.
9. How choices were judged: structure on left-out events, the closed form by AICc on the TRAIN points.
10. "A band difference up to 25 % cannot be excluded" was unsupported and was removed. P4 compares each class with the overall
    TEST / VAL ratio. The light-rain levels are 1.35-1.55 dB against 0.83-1.11 dB.
11. The two evaluation errors of S8 are named. The cause of the remaining error is marked as an interpretation.
12. The lab result on runoff now matches what the paper reports, and the title of reference [2] is corrected (Tiede, Chwala,
    Siart and Eibert 2025).
13. The next steps are R1, then Part 2.
14. Figures 1, 2, 4, 5 and 6 are real data. Figure 3 (the two stores) is an illustration from the equation and is captioned so.

## 4. The documents

- `docs/WAA_research_overview.tex` / `.pdf`: the one WAA document. It keeps the structure, wording and look of the first
  overview, with only the corrections above. S9 first wrote a short and a full (school) version in `docs/overview/`; the learners
  then merged them into this one document and deleted the folder (2026-10-19).
- A4, 11 pt, numbered sections and equations, no semicolons, no bold in the text. Built with
  `python docs/build_research_log.py WAA_research_overview` (no LaTeX error, no overfull box).
- Figures `docs/WAA_research_overview_figures/` (6), drawn by its `figs.py` from saved tables only (S3, S6, S7 and the tables of
  this step: `formulate_s9_example_event.csv`, `_memory.csv`, `_closed_f_band.csv`, `_dry_curves.csv` (9081 dry-path
  windows)). Black and grey only, with line styles, markers and text labels.

## 5. Open items

- R1, the next step: does the WAA correction improve link rain estimates against the gauges? It needs its own written plan and
  a decision about TEST, which was already read once.
- Then Part 2, the path law f(L, R̄, Z), with this equation removing W first.
- Box C stays open: the WAA level in light rain varies between storms (TEST P1), a_link is not explained by any link property,
  F above 10 mm/h, P below 1 mm/h, and most of the error during rain within storms.
- A small oddity, reported, not changed: the bootstrap interval of the 90th percentile of a_link
  (1.29-1.40) does not contain the fitted value (1.24).
  The cold draws spread the link multipliers more than the full fit (S8b table).
- The fit caches (`research/.cache/wet_antenna/`) are still to be published now that the WAA stage is finished (ROADMAP part A).
