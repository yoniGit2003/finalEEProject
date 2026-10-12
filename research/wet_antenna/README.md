# wet_antenna: drying tails and a dynamic equation for wet-antenna attenuation (WAA)

## Purpose

A(t) = TL(t) - baseline = A_rain(t) + δ(t), where δ(t) is the non-rain part (WAA + baseline residual). In the
**drying tail** after the rain on a link has stopped, A_rain = 0, so A(t) = δ(t): what is left there is
mostly WAA. This module isolates WAA from the tails **before** f(L, R̄, Z) is learned. Otherwise a constant WAA looks
like f = L + L0. The later steps learn a dynamic model of the wetting / drying with weak-form SINDy (rain as input) and
compare it with the literature models [1], [9]. Tails are never selected by how well a model fits them (that would be
circular).

## Status

| step | what | state |
|---|---|---|
| 0 | inventory of the inputs (`inventory.py`) | done 2026-10-08 |
| 1 | extract ALL drying tails with their features, descriptive report (`extract_tails.py`) | done 2026-10-08; no selection, no fitting |
| 2 | pre-defined exclusions / truncation, main set, fast-drop test, shape descriptors (`select_tails.py`) | done 2026-10-08; TRAIN only, no fitting |
| 2b | fixes after the step 2 gallery: t0_trace start, strict / lenient ends, earlier rise cut, noise exclusion, level flags (`refine_selection.py`) | done 2026-10-08; TRAIN only, no fitting |
| 2c | last fixes (radar NaN = unknown, rise cut v2), event diagnostics, **frozen tail set v1** (`freeze_tail_set.py`) | done 2026-10-08; frozen |
| 3a | tail simulator (drying laws + observation model) calibrated to v1 (`calibrate_simulator.py`, library `simulator.py`) | done 2026-10-08; 3a-fix the same day: slow AR(1) noise, spread for every law, m2 dropped |
| 3b-1 | the per-tail fitting method on simulated tails with a known truth (`check_per_tail_fit.py`, library `per_tail_fit.py`) | done 2026-10-09: only "exponential vs not" is identifiable, per dataset with BIC |
| 3b-1b | pooled fit (shared law, per-tail start and speed, offset prior, AR(1) noise likelihood, F6 added), re-tested on event-bootstrap simulated tails (`check_pooled_fit.py`, library `pooled_fit.py`) | done 2026-10-10: pre-set rule met in all three windows (own law 80-100 % for m1, m4, m5, m6; parameters within +-15 %); m3 vs m5 partly confused |
| 3b-1c | stress test of the pooled fit on fake tails: unknown noise, honest starts, faster truths (`stress_pooled_fit.py`) | done 2026-10-11: PARTIAL (12 of 14 cells; m1, m4, m6 pass every test; m5 fails: F3 / F5 not separable) |
| 3b-2 | weak-form SINDy on fake tails with a known truth (`check_sindy.py`, library `sindy.py`) | done 2026-10-11: Q1 FAIL, Q2 FAIL, Q3 PARTIAL, Q4 FAIL: it returns "A" for almost every truth; the pooled fit stays the working tool |
| 3b-3 | "draw the curve, then name it" on fake tails: free-curve pooled fit, then the equation of the curve (`check_free_curve.py`, library `free_curve.py`) | done 2026-10-11: PARTIAL by the pre-set rule (m1 FAIL, m3 / m5 / m6 PARTIAL, m4 PASS); the curve is right above about 0.35 dB and undetermined below |
| 3b-3b | confirmation run of the two jumps on fresh fake tails, judged only above 0.35 dB, rule fixed before the run (`confirm_free_curve.py`) | done 2026-10-11: PARTIAL by the pre-set rule (m1, m3, m5 PASS on curve, shape and equation above 0.35 dB; the m4 warning sign fails, 2 of 5; m6 not judged) |
| 3b-3c | the discovery method finished on fake tails: jump 2 alone, jump 1 with an optional step, jump 3 (refit and choose), bootstrap band, F4 check (`check_discovery_method.py`) | done 2026-10-11: PARTIAL (jump 2 alone works with E-SINDy / PySR; the step is drawn for every truth; jump 3 and the bootstrap band fail; F4 check: no change) |
| 3b-4 | is the 0.3 dB RSL quantization the cause of the false step of 3b-3c? (`check_quantization.py`) | done 2026-10-11: NOT the cause (without quantization the step is still kept in 5 of 5); cause open |
| 4 | the drying law on the REAL tails: menu fits, free curve above 0.35 dB with a bootstrap band, jump 2, five checks (`fit_real_drying.py`) | done 2026-10-12: NO law identified (the real curve falls with W, slope -0.87; -0.05 without the rise-heavy events; no formula; menu not decisive; a W0 fault found in the menu comparison) |
| 6+ | what drying depends on (6); full-event simulator (3c); wetting (7) | not started; Step 6 waits for the learners' decisions on Step 4 (METHOD.md) |

Inputs: release `openmrg_v1` (pinned in `research/config.yaml`) and the frozen baseline v2 products of
`research/baseline/` (`baselines.csv`, `baseline_uncertainty.csv`, `dew_risk_hourly.csv`, `dry_jumps.csv`).

The step numbers (0 ... 3b-4, 4) are the ones the research log ([`../../docs/research_log.pdf`](../../docs/research_log.pdf))
and [`METHOD.md`](METHOD.md) use; each row names the script that does the step.

## How to run

From the repository root: `python -m research` lists every step with its command and run time, and
`python -m research run <step>` runs one (full instructions: [`../../docs/REPRODUCE.md`](../../docs/REPRODUCE.md)).
Each script can also be run directly, e.g. `python -m research.wet_antenna.extract_tails`. Settings: section
`wet_antenna` of [`../config.yaml`](../config.yaml).

- Steps 0 to 3a read the data release and the baseline products. Steps 3b-1 to 4 read only files of this module.
- The steps that fit (3b-1 to 4) save every finished fit in `research/.cache/wet_antenna/<script name>/`, so a second
  run takes seconds to minutes. `python pipeline/get_caches.py` downloads our saved fits; delete a folder to fit again.
- Three big files are not in git and are rebuilt by steps 1, 2 and 2b: `outputs/products/tails_1min.parquet` (47 MB),
  `outputs/products/tails_selected.csv` (17 MB) and `outputs/tables/check_per_tail_fit_fits.parquet` (18 MB, step 3b-1).
- Step 2 rewrites `tails_selected.csv` without the step 2b columns, so always run step 2b after step 2.

## Code

| file | what it does |
|---|---|
| `paths.py` | every folder and file location of this module |
| `inventory.py` | measures every input (dims, columns, spans, coverage), prints one example event, writes `outputs/reports/inventory.md` |
| `extract_tails.py` | cuts every tail, writes the two products, the A(t) check, the descriptive tables / figures and `outputs/reports/extract_tails.md` |
| `select_tails.py` | applies the step 2 rules, writes `tails_selected.csv`, the fast-drop test, shape descriptors, curves, gallery and `outputs/reports/select_tails.md` |
| `refine_selection.py` | applies the step 2b fixes, adds their columns to `tails_selected.csv`, writes the step 2b tables / figures and `outputs/reports/refine_selection.md` |
| `freeze_tail_set.py` | applies the step 2c fixes, event diagnostics, writes the frozen set (`tails_model*`) and `outputs/reports/freeze_tail_set.md` |
| `simulator.py` | library: the drying laws m1, m3, m4, m5 (closed forms), the per-tail time-scale spread, the observation model with switches (offset, white noise, slow AR(1) noise, quantization, real NaN mask, spread), the dry-noise lookup |
| `per_tail_fit.py` | library: the candidate models F1, F3, F4, F5 (closed forms, free W0 and offset c), per-tail least squares, BIC / AIC with an effective n, winners, t_half and dA/dt population statistics |
| `pooled_fit.py` | library: the pooled marginal likelihood of a dataset for F1, F3, F4, F5, F6 (numba): ARMA(1,1) whitening of white + AR(1) noise, offset prior (MAP), per-tail W0 profiled, speed factors marginalised (Laplace per tail, quadrature per link), Nelder-Mead fit |
| `sindy.py` | library (numpy only; pysindy is not installed): weak-form rows of a tail, library 1, A, A^2, A^3, max(A - A_c, 0), thresholded backward elimination, threshold by cross-validation over events, stability selection; `python sindy.py` runs a noise-free self-test |
| `check_sindy.py` | weak-form SINDy on the fake datasets of 3b-1c (+ m3, + no-spread variants): two windows, three offset treatments, time-rescaling check, pre-set rule, agreement with the pooled fit; cached in `.cache/check_sindy/`; writes `outputs/reports/check_sindy.md` |
| `free_curve.py` | library: the free drying curve dW/dt = -f W k(W) (log k = natural cubic spline in log W) for the pooled likelihood of `pooled_fit.py` (law `FC`), penalised fit, curvature-based uncertainty; and the equation of a curve (sparse regression of h(W) = W k(W) on W, W^2, W^3, switch terms) |
| `check_free_curve.py` | "draw the curve, then name it" on the fake datasets of 3b-1c (+ m3, m1 without spread, T1 noise): free-curve pooled fit per dataset and window (4 starts, penalty by cross-validation over events; cached in `.cache/check_free_curve/`), band and shape tests, equation of the curve, m4 warning signs, comparison with the menu fit; writes `outputs/reports/check_free_curve.md` |
| `confirm_free_curve.py` | confirmation run of 3b-3 on fresh fake tails (seeds 4000 + replicate, first 60 min): the same free-curve pooled fit with 5 knots from 0.35 dB to W_hi, nothing judged below 0.35 dB, band = min-max of the replicates, slope per replicate, equation of the curve; rule fixed before the run (`wet_antenna.confirm_free_curve`); cached in `.cache/confirm_free_curve/`; writes `outputs/reports/confirm_free_curve.md` |
| `check_discovery_method.py` | the rest of the discovery method on the fake datasets of 3b-3b: A jump 2 alone on true curves (best-subset + BIC, E-SINDy), B free curve with an optional step (`free_curve.StepCurve`), C jump 3 = every candidate (jump-2 formulas via `waa_curve.FormulaLaw`, step law, menu laws) refitted in the pooled model and chosen by held-out log-likelihood over event draws, D band from an event bootstrap of one dataset, E F4 with the F1-solution start; rules fixed before the run (`wet_antenna.check_discovery_method`); cached in `.cache/check_discovery_method/`; writes `outputs/reports/check_discovery_method.md` |
| `check_quantization.py` | diagnosis only: truth m1 simulated with no / 0.3 dB / 0.1 dB RSL quantization (same random numbers), free curve of 3b-3b + optional step of 3b-3c part B on each; rule fixed before the run (`wet_antenna.check_quantization`); the 0.3 dB fits are read from `.cache/confirm_free_curve/` and `.cache/check_discovery_method/`; cached in `.cache/check_quantization/`; writes `outputs/reports/check_quantization.md` |
| `fit_real_drying.py` | the first real-data step: TRAIN tails of the frozen set v1, first 60 min; menu fits F1, F3, F4, F5, F6 from 4 starts, free curve of 3b-3b (5 knots from 0.35 dB to W_hi) with 20 event-bootstrap refits, jump 2 (E-SINDy, PySR in its own process, best-subset), four sensitivity sets, leave-one-big-event-out, 6 h window, residuals vs dry noise, gallery, comparison with the fake truths; settings `wet_antenna.fit_real_drying` (fit settings read from the fake-data steps); fits cached in `.cache/fit_real_drying/`; writes `outputs/reports/fit_real_drying.md` |
| `stress_pooled_fit.py` | stress test of the pooled fit on fake tails (first 60 min; truths m1, m4, m5, m6): T1 noise with two slow time scales, T2 generic optimiser starts, T3 truths recalibrated to the real early decay; parallel fits cached in `.cache/stress_pooled_fit/`; win shares, recovery, start sensitivity, pre-set rule per cell; writes `outputs/reports/stress_pooled_fit.md` |
| `check_pooled_fit.py` | event-bootstrap replicates per truth (m1, m3, m4, m5, m6), simulates them, pooled fit of every candidate in three windows (cached in `.cache/check_pooled_fit/`), confusion, recovery, pre-set decision rule, population statistics; writes `outputs/reports/check_pooled_fit.md` |
| `check_per_tail_fit.py` | simulates replicate datasets per truth (m1, m3, m4, m5, m6), fits every tail with every candidate (cached in `.cache/check_per_tail_fit/`), confusion, recovery, population tests, ablation; writes `outputs/reports/check_per_tail_fit.md` |
| `calibrate_simulator.py` | calibrates the noise (control tails: std, on-level share, ACF 1 / 10 / 60 min) and every law (median, then spread from the band), simulates every law on the TRAIN v1 skeleton, checks real vs simulated (incl. t_half per A_start bin), writes `outputs/reports/calibrate_simulator.md` |

## Outputs

`outputs/products/` (the tail sets and the simulated tails that later steps read: `tails_index.csv`,
`tails_model*` = the frozen set v1, `sim_tails_*`), `outputs/tables/`, `outputs/figures/` and `outputs/reports/`.
Every file starts with the name of the script that writes it, and each script writes one report,
`outputs/reports/<script>.md`, which lists its tables and figures and explains them. Everything is committed except
the three big files above. The file locations are all in [`paths.py`](paths.py).

## Method, decisions and open issues

[`METHOD.md`](METHOD.md): the rules that define a tail and the frozen tail set v1, what every step decided (with the
numbers behind it), and the open issues. The plan for the next steps: [`../../docs/ROADMAP.md`](../../docs/ROADMAP.md).
