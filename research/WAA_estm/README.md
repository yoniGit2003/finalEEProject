# WAA_estm: drying tails and a dynamic equation for wet-antenna attenuation (WAA)

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
| 0 | inventory of the inputs (`step0_inventory.py`) | done 2026-10-08 |
| 1 | extract ALL drying tails with their features, descriptive report (`step1_extract_tails.py`) | done 2026-10-08; no selection, no fitting |
| 2 | pre-defined exclusions / truncation, main set, fast-drop test, shape descriptors (`step2_select_tails.py`) | done 2026-10-08; TRAIN only, no fitting |
| 2b | fixes after the step 2 gallery: t0_trace start, strict / lenient ends, earlier rise cut, noise exclusion, level flags (`step2b_fix_tails.py`) | done 2026-10-08; TRAIN only, no fitting |
| 2c | last fixes (radar NaN = unknown, rise cut v2), event diagnostics, **frozen tail set v1** (`step2c_freeze.py`) | done 2026-10-08; frozen |
| 3a | tail simulator (drying laws + observation model) calibrated to v1 (`step3a_simulator.py`, library `waa_sim.py`) | done 2026-10-08; 3a-fix the same day: slow AR(1) noise, spread for every law, m2 dropped |
| 3b-1 | the per-tail fitting method on simulated tails with a known truth (`step3b1_fitting_on_sim.py`, library `waa_fit.py`) | done 2026-10-09: only "exponential vs not" is identifiable, per dataset with BIC |
| 3b-1b | pooled fit (shared law, per-tail start and speed, offset prior, AR(1) noise likelihood, F6 added), re-tested on event-bootstrap simulated tails (`step3b1b_pooled_fit.py`, library `waa_pooled.py`) | done 2026-10-10: pre-set rule met in all three windows (own law 80-100 % for m1, m4, m5, m6; parameters within +-15 %); m3 vs m5 partly confused |
| 3b-1c | stress test of the pooled fit on fake tails: unknown noise, honest starts, faster truths (`step3b1c_stress_test.py`) | done 2026-10-11: PARTIAL (12 of 14 cells; m1, m4, m6 pass every test; m5 fails: F3 / F5 not separable) |
| 3b-2 | weak-form SINDy on fake tails with a known truth (`step3b2_sindy_fake.py`, library `waa_sindy.py`) | done 2026-10-11: Q1 FAIL, Q2 FAIL, Q3 PARTIAL, Q4 FAIL: it returns "A" for almost every truth; the pooled fit stays the working tool |
| 3b-3 | "draw the curve, then name it" on fake tails: free-curve pooled fit, then the equation of the curve (`step3b3_curve_then_equation.py`, library `waa_curve.py`) | done 2026-10-11: PARTIAL by the pre-set rule (m1 FAIL, m3 / m5 / m6 PARTIAL, m4 PASS); the curve is right above about 0.35 dB and undetermined below |
| 3b-3b | confirmation run of the two jumps on fresh fake tails, judged only above 0.35 dB, rule fixed before the run (`step3b3b_confirm.py`) | done 2026-10-11: PARTIAL by the pre-set rule (m1, m3, m5 PASS on curve, shape and equation above 0.35 dB; the m4 warning sign fails, 2 of 5; m6 not judged) |
| 3b-3c | the discovery method finished on fake tails: jump 2 alone, jump 1 with an optional step, jump 3 (refit and choose), bootstrap band, F4 check (`step3b3c_method.py`) | done 2026-10-11: PARTIAL (jump 2 alone works with E-SINDy / PySR; the step is drawn for every truth; jump 3 and the bootstrap band fail; F4 check: no change) |
| 3b-4 | is the 0.3 dB RSL quantization the cause of the false step of 3b-3c? (`step3b4_quantization.py`) | done 2026-10-11: NOT the cause (without quantization the step is still kept in 5 of 5); cause open |
| 4 | the drying law on the REAL tails: menu fits, free curve above 0.35 dB with a bootstrap band, jump 2, five checks (`step4_real_drying.py`) | done 2026-10-12: NO law identified (the real curve falls with W, slope -0.87; -0.05 without the rise-heavy events; no formula; menu not decisive; a W0 fault found in the menu comparison) |
| 6+ | what drying depends on (6); full-event simulator (3c); wetting (7) | not started; Step 6 waits for the learners' decisions on Step 4 (DECISIONS.md) |

Inputs: release `openmrg_v1` (pinned in `research/config.yaml`) and the frozen baseline v2 products of
`research/baseline_estm/` (`baselines.csv`, `baseline_uncertainty.csv`, `dew_risk_hourly.csv`, `dry_jumps.csv`).

## How to run

Commands from the repository root with `C:\Users\yoni\miniconda3\python.exe` (written `python` below). Settings:
section `waa` of `research/config.yaml`.

```powershell
python research/WAA_estm/step0_inventory.py        # inventory + one example (prints the report)       ~ 10 s
python research/WAA_estm/step1_extract_tails.py    # all tails, products, report (prints the report)   ~ 2.5 min
python research/WAA_estm/step2_select_tails.py     # selection + fast-drop test (needs step 1)          ~ 1.5 min
python research/WAA_estm/step2b_fix_tails.py       # step 2b fixes, adds columns (run after step 2)    ~ 1 min
python research/WAA_estm/step2c_freeze.py          # last fixes + frozen tail set v1 (after 2b)          ~ 45 s
python research/WAA_estm/step3a_simulator.py       # simulator calibration + checks (needs v1, 2b set)   ~ 3.5 min
python research/WAA_estm/step3b1_fitting_on_sim.py # fits on simulated tails (needs 3a), 4 workers      ~ 60 min (cached: 2.5 min)
python research/WAA_estm/step3b1b_pooled_fit.py    # pooled fits on simulated tails (needs 3a, 3b-1)     ~ 3 h (cached: 1 min)
python research/WAA_estm/step3b1c_stress_test.py   # stress test of the pooled fit (needs 3a, 3b-1), cores - 1   ~ 45 min (cached: 20 s)
python research/WAA_estm/step3b2_sindy_fake.py     # weak-form SINDy on fake tails (needs 3a, 3b-1; reuses 3b-1c fits)  ~ 13 min (cached: 10 s)
python research/WAA_estm/step3b3_curve_then_equation.py   # free drying curve + its equation on fake tails (needs 3a, 3b-1; reads 3b-2's cached menu fits)  ~ 2.1 h (cached: 15 s)
python research/WAA_estm/step3b3b_confirm.py              # confirmation run on fresh fake tails, judged above 0.35 dB (needs 3a, 3b-1)  ~ 35 min (cached: 15 s)
python research/WAA_estm/step3b3c_method.py               # jump 2 alone, step, jump 3, bootstrap band, F4 check on fake tails (needs 3b-3b's cache)  ~ 1.6 h + 15 min PySR with `--pysr` (cached: 80 s)
python research/WAA_estm/step3b4_quantization.py         # is the 0.3 dB quantization the cause of the false step? m1, three settings (needs 3b-3b's and 3b-3c's caches)  ~ 17 min (cached: 20 s)
python research/WAA_estm/step4_real_drying.py            # the drying law on the REAL TRAIN tails: menu fits, free curve, jump 2, checks (needs steps 2c, 3a; reads 3b-3b's cache for check 5 if present; PySR optional)  ~ 40 min (cached: 25 s)
```

Both need the release and the baseline products. Step 1 reads the 1 min CML store once, chunk by chunk (32 sub-links).
`outputs/products/tails_1min.parquet` (~47 MB) and `outputs/products/tails_selected.csv` (~13 MB) are git-ignored
(over the 10 MB limit). Rebuild them with steps 1, 2 and 2b. Step 2 rewrites `tails_selected.csv` without the
step 2b columns, so always run 2b after 2.

## Code

| file | what it does |
|---|---|
| `WAA_estm_paths.py` | every folder and file location of this module |
| `step0_inventory.py` | measures every input (dims, columns, spans, coverage), prints one example event, writes `outputs/reports/step0_inventory.md` |
| `step1_extract_tails.py` | cuts every tail, writes the two products, the A(t) check, the descriptive tables / figures and `outputs/reports/step1_tails.md` |
| `step2_select_tails.py` | applies the step 2 rules, writes `tails_selected.csv`, the fast-drop test, shape descriptors, curves, gallery and `outputs/reports/step2_selection.md` |
| `step2b_fix_tails.py` | applies the step 2b fixes, adds their columns to `tails_selected.csv`, writes the step 2b tables / figures and `outputs/reports/step2b_selection.md` |
| `step2c_freeze.py` | applies the step 2c fixes, event diagnostics, writes the frozen set (`tails_model*`) and `outputs/reports/step2c_freeze.md` |
| `waa_sim.py` | library: the drying laws m1, m3, m4, m5 (closed forms), the per-tail time-scale spread, the observation model with switches (offset, white noise, slow AR(1) noise, quantization, real NaN mask, spread), the dry-noise lookup |
| `waa_fit.py` | library: the candidate models F1, F3, F4, F5 (closed forms, free W0 and offset c), per-tail least squares, BIC / AIC with an effective n, winners, t_half and dA/dt population statistics |
| `waa_pooled.py` | library: the pooled marginal likelihood of a dataset for F1, F3, F4, F5, F6 (numba): ARMA(1,1) whitening of white + AR(1) noise, offset prior (MAP), per-tail W0 profiled, speed factors marginalised (Laplace per tail, quadrature per link), Nelder-Mead fit |
| `waa_sindy.py` | library (numpy only; pysindy is not installed): weak-form rows of a tail, library 1, A, A^2, A^3, max(A - A_c, 0), thresholded backward elimination, threshold by cross-validation over events, stability selection; `python waa_sindy.py` runs a noise-free self-test |
| `step3b2_sindy_fake.py` | weak-form SINDy on the fake datasets of 3b-1c (+ m3, + no-spread variants): two windows, three offset treatments, time-rescaling check, pre-set rule, agreement with the pooled fit; cached in `.cache/sindy3b2/`; writes `outputs/reports/step3b2_sindy_fake.md` |
| `waa_curve.py` | library: the free drying curve dW/dt = -f W k(W) (log k = natural cubic spline in log W) for the pooled likelihood of `waa_pooled.py` (law `FC`), penalised fit, curvature-based uncertainty; and the equation of a curve (sparse regression of h(W) = W k(W) on W, W^2, W^3, switch terms) |
| `step3b3_curve_then_equation.py` | "draw the curve, then name it" on the fake datasets of 3b-1c (+ m3, m1 without spread, T1 noise): free-curve pooled fit per dataset and window (4 starts, penalty by cross-validation over events; cached in `.cache/fits3b3/`), band and shape tests, equation of the curve, m4 warning signs, comparison with the menu fit; writes `outputs/reports/step3b3_curve_then_equation.md` |
| `step3b3b_confirm.py` | confirmation run of 3b-3 on fresh fake tails (seeds 4000 + replicate, first 60 min): the same free-curve pooled fit with 5 knots from 0.35 dB to W_hi, nothing judged below 0.35 dB, band = min-max of the replicates, slope per replicate, equation of the curve; rule fixed before the run (`waa.curve3b3b`); cached in `.cache/fits3b3b/`; writes `outputs/reports/step3b3b_confirm.md` |
| `step3b3c_method.py` | the rest of the discovery method on the fake datasets of 3b-3b: A jump 2 alone on true curves (best-subset + BIC, E-SINDy), B free curve with an optional step (`waa_curve.StepCurve`), C jump 3 = every candidate (jump-2 formulas via `waa_curve.FormulaLaw`, step law, menu laws) refitted in the pooled model and chosen by held-out log-likelihood over event draws, D band from an event bootstrap of one dataset, E F4 with the F1-solution start; rules fixed before the run (`waa.curve3b3c`); cached in `.cache/fits3b3c/`; writes `outputs/reports/step3b3c_method.md` |
| `step3b4_quantization.py` | diagnosis only: truth m1 simulated with no / 0.3 dB / 0.1 dB RSL quantization (same random numbers), free curve of 3b-3b + optional step of 3b-3c part B on each; rule fixed before the run (`waa.quant3b4`); the 0.3 dB fits are read from `.cache/fits3b3b/` and `.cache/fits3b3c/`; cached in `.cache/fits3b4/`; writes `outputs/reports/step3b4_quantization.md` |
| `step4_real_drying.py` | the first real-data step: TRAIN tails of the frozen set v1, first 60 min; menu fits F1, F3, F4, F5, F6 from 4 starts, free curve of 3b-3b (5 knots from 0.35 dB to W_hi) with 20 event-bootstrap refits, jump 2 (E-SINDy, PySR in its own process, best-subset), four sensitivity sets, leave-one-big-event-out, 6 h window, residuals vs dry noise, gallery, comparison with the fake truths; settings `waa.real4` (fit settings read from the fake-data steps); fits cached in `.cache/fits4/`; writes `outputs/reports/step4_real_drying.md` |
| `step3b1c_stress_test.py` | stress test of the pooled fit on fake tails (first 60 min; truths m1, m4, m5, m6): T1 noise with two slow time scales, T2 generic optimiser starts, T3 truths recalibrated to the real early decay; parallel fits cached in `.cache/fits3b1c/`; win shares, recovery, start sensitivity, pre-set rule per cell; writes `outputs/reports/step3b1c_stress_test.md` |
| `step3b1b_pooled_fit.py` | event-bootstrap replicates per truth (m1, m3, m4, m5, m6), simulates them, pooled fit of every candidate in three windows (cached in `.cache/pooled3b1b/`), confusion, recovery, pre-set decision rule, population statistics; writes `outputs/reports/step3b1b_pooled_fit.md` |
| `step3b1_fitting_on_sim.py` | simulates replicate datasets per truth (m1, m3, m4, m5, m6), fits every tail with every candidate (cached in `.cache/fits3b1/`), confusion, recovery, population tests, ablation; writes `outputs/reports/step3b1_fitting_on_sim.md` |
| `step3a_simulator.py` | calibrates the noise (control tails: std, on-level share, ACF 1 / 10 / 60 min) and every law (median, then spread from the band), simulates every law on the TRAIN v1 skeleton, checks real vs simulated (incl. t_half per A_start bin), writes `outputs/reports/step3a_fix.md` |

## Tail definition (decided by the learners; re-stated here, not read from `explorations/`)

One tail per (event, sub-link) pair:

- **Pairs**: every row of `baselines.csv` with `missing == False` and `link_has_wet_slot == True`. These are train + val
  events on train links (test events and links have no baseline). All flag columns are kept, and `split_event` is
  recorded on every tail.
- **A(t)** = 1 min `tl_mean` (release) - `baseline_db` (`baselines.csv`). This is exactly `baseline.a_of_t` (checked on 3
  pairs: `outputs/tables/step1_a_check.csv`).
- **t0** = `last_wet_end`: the first non-wet link_state bin after the link's last wet bin in the network event. The
  ±30 min unknown margin of link_state is **inside** the tail. Short pauses inside the event are not tails.
- **t_first_dry** = the first link_state == dry bin at or after t0 (searched over the whole season).
- **t_end** = the earliest of three stops:
  - the next **wet** link_state bin of this link, searched over the whole season (not only inside the event);
  - the start of a NaN run of A ≥ `gap_min` = 15 min;
  - t0 + `cap_h` = 12 h.

  `ended_by` is next_wet / gap / cap, or season_end when the data stop first. The end never looks at the value of A,
  so a model is not forced to reach 0.
- **No selection** in step 1: no threshold on A, no shape filter, all flags kept as columns. Selection rules are a
  later, separate step.

Inherited from earlier learner decisions:
- the tail start is the pipeline's link_state as is (no extra gauges);
- tails identify only the **drying** (R = 0); the **wetting** needs in-event data where A_rain is mixed in;
- dry-weather excursions (dew?) are WAA without rain, to be used as a check, not as training data.

## Step 2 rules (fixed before the run; settings `waa.selection`)

TRAIN events only for every count, figure and test. VAL tails carry the same columns and are never looked at.

1. **Exclude** (reasons recorded in `exclusion_reason`): t0 is a wet bin; length 0; A_start NaN;
   `use_in_waa_main` False; sub-link status not good; more than 10 % of the tail's minutes without A.
2. **Truncate** the end (`t_use_end`, `truncated_by`) at the earliest of:
   - just before the first `tsl_step` minute;
   - the start of the first sustained dry-weather rise: the 15 min trailing median of A is >= 0.6 dB above its running
     minimum since t_first_dry, on dry minutes, for >= 30 consecutive minutes;
   - otherwise the step 1 `t_end`.
3. **Main set**: kept, A_start >= 1.0 dB, usable length >= 60 min (alternative A_start >= 0.6 dB: counts only).
   **Control**: kept, A_start < 0.5 dB.
4. Nothing is selected by shape. The fast-drop test and the shape descriptors are descriptive.

## Step 2b rules (after the step 2 gallery; settings `waa.selection_2b`)

TRAIN only, as in step 2. These rules replace the step 2 start and ends. Step 2 columns stay in the product for
comparison.

1. **Start t0_trace**: the end of the last 5 min bin with ANY rain trace (radar R_max > 0 on the path, or any gauge
   within 1 km > 0) between t0 and t_first_dry; t0 if the margin has no trace. The radar's presence counts even when
   its value is below the 0.1 mm/h wet threshold. `A_start_trace_db` = mean A over [t0_trace, + 5 min).
2. **Two ends**:
   - lenient = the step 2 ends with the earlier rise cut (3);
   - strict = lenient, or the first rain-trace bin after t0_trace if earlier.

   Usable length is counted from t0_trace.
3. **Dry-rise cut**: same detection as step 2. The cut goes back to the last minute before the rise run whose smoothed
   excess was < 0.15 dB, or to t_first_dry if there is none.
4. **Noise exclusion**: sub-links with `offset_std` > 0.3 dB (one RSL step; `baseline_uncertainty.csv`, measured in
   dry weather, never on tails).
5. **Flags, not exclusions**:
   - `level_shift`: |median A over [t0 + 6 h, t0 + 12 h]| > max(0.3 dB, 2 x offset_std);
   - `undershoot`: the 15 min median of A < -0.3 dB within 3 h after t0_trace.
6. **Main sets**: kept, A_start_trace >= 1.0 dB, usable length >= 60 min, for the strict and the lenient end.
   **Control**: kept, A_start_trace < 0.5 dB.

## Tail set v1 (FROZEN 2026-10-08; the input of the WAA model)

Files: `outputs/products/tails_model.parquet` (1 min rows) + `tails_model_index.csv` (one row per tail); the
alternative set `tails_model_alt.*`. The set id `v1` is in config `waa.freeze.tail_set` and in the `tail_set` column.
File names carry no version, because ORGANIZATION.md forbids `_v1` on active paths. A later set gets a new id and
the old files go to `archive/`. Every rule that defines v1, in order:

1. **Pairs**: `baselines.csv` (frozen baseline v2) rows with `missing == False` and `link_has_wet_slot == True`, i.e.
   train + val events on train links.
2. **A(t)** = 1 min `tl_mean` - `baseline_db`.
3. **t0** = `last_wet_end` (first non-wet link_state bin after the link's last wet bin in the event).
4. **Exclusions**:
   - t0 is a wet bin; length 0; A_start NaN;
   - `use_in_waa_main` False; sub-link status not good; NaN share of the step 1 tail > 10 %;
   - `offset_std` > 0.3 dB (`baseline_uncertainty.csv`);
   - A_start_trace NaN.
5. **Rain trace** (5 min bin): radar R_max > 0 on the path, OR any gauge within 1 km > 0, OR radar R_max NaN on the
   path (radar NaN = unknown).
6. **Start t0_trace** = end of the last rain-trace bin between t0 and t_first_dry (t0 if none). **A_start_trace** =
   mean A over [t0_trace, t0_trace + 5 min).
7. **End** = the earliest of:
   - the step 1 stops: next wet bin of the link (whole season), a NaN run of A ≥ 15 min, t0 + 12 h, season end;
   - just before the first `tsl_step` minute;
   - the step 2b dry-rise cut: 15 min median ≥ 0.6 dB above its running minimum since t_first_dry on dry minutes,
     for ≥ 30 min; cut at the last minute before it with excess < 0.15 dB;
   - the first rain-trace bin after t0_trace (rule 5);
   - **rise cut v2**: from t0_trace on ALL minutes, the 10 min trailing median of A ≥ 0.6 dB above its running
     minimum since t0_trace for ≥ 10 min; cut at the last minute before that run with excess < 0.15 dB (that minute
     is kept), or at t0_trace if there is none.
8. **Main set**: A_start_trace ≥ 1.0 dB and usable length (t0_trace .. end) ≥ 60 min. **Alternative set**:
   A_start_trace ≥ 0.6 dB, same length rule (it contains the main set; column `in_main_set`).
9. **Flags (kept as columns, no exclusion)**:
   - `level_shift`: |median A over [t0 + 6 h, t0 + 12 h]| > max(0.3 dB, 2 x offset_std), step 2b ends;
   - `undershoot`: 15 min median of A < -0.3 dB within 3 h after the step 2b t0_trace;
   - `event_rise_heavy`: more than 30 % of the event's pre-rise tails are cut by rise v2.
10. **Split**: TRAIN for development; VAL is built with the same rules and not looked at. Test events and links have
    no baseline yet, so there are no test tails.

Settings: `waa.tails`, `waa.selection`, `waa.selection_2b`, `waa.freeze` in `research/config.yaml`. Do not change them
for v1.

## Outputs

`outputs/products/` (`tails_1min.parquet`, `tails_index.csv`, `tails_selected.csv`, `tails_model*` = frozen set v1), `outputs/tables/`, `outputs/figures/`,
`outputs/reports/`. Every file is described in [`CATALOG.md`](CATALOG.md); its importance level (key / product / diag /
scratch, ORGANIZATION.md section 5) is in `outputs/MANIFEST.csv`.

- **Key** (kept for the final report, 6 = the limit): `outputs/reports/step2c_freeze.md` with its figure
  `step2c_median_curves.png`; the simulator calibration `outputs/tables/step3a_calibration.csv` and
  `outputs/figures/step3a_real_vs_sim_median.png`; the pooled fitting test on simulated tails
  `outputs/figures/step3b1b_confusion.png` and `outputs/tables/step3b1b_recovery.csv`. Demoted to diag: the two step
  2c event look figures (2026-10-08, learners' choice), the step 2 fast-drop figure and table (2026-10-09, learners'
  choice) and the step 3b-1 per-tail confusion figure and recovery table (2026-10-10, decision delegated to Claude:
  superseded by 3b-1b); their findings are in DECISIONS.md. The key report `step2c_freeze.md` still links the event figures, so /tidy keeps
  them and reports them as a conflict until that report is rerun or edited.
- **Product**: everything in `outputs/products/`; the step 3b-1c report `outputs/reports/step3b1c_stress_test.md`
  and its tables (win shares, decision, recovery, start sensitivity, recalibration). Its win-share figure is diag
  (the key limit is reached). The step 3b-2 report `outputs/reports/step3b2_sindy_fake.md` and its tables (selection,
  coefficients, decision, agreement); its selection figure is diag. The step 3b-3 report
  `outputs/reports/step3b3_curve_then_equation.md` and its tables (curves, jump 1, jump 2, menu comparison, decision,
  post hoc upper range); its five k(W) figures are diag. The step 3b-3b report
  `outputs/reports/step3b3b_confirm.md` and its tables (curves, curve and shape, jump 2, m4 BIC, decision); its five
  k(W) figures are diag. The step 3b-3c report `outputs/reports/step3b3c_method.md` and its tables (parts A-E,
  decision); its figures (step, bootstrap band) are diag. The step 3b-4 report
  `outputs/reports/step3b4_quantization.md` and its two tables; its figure is diag.
  The step 4 report `outputs/reports/step4_real_drying.md` and its six tables (menu fits, summary, curves, jump 2,
  fit quality, fake comparison); its three figures are diag, `step4_k_real.png` is the key candidate (not promoted:
  the key limit is reached).
- **Diag** (all other step figures, tables and reports): moved by /tidy to the local `legacy/` folder (git-ignored;
  log: `outputs/TIDY_LOG.md`). The files /tidy deleted before 2026-10-11 are not there: rerun the
  step's script to get them back. The decisions they held are in [`DECISIONS.md`](DECISIONS.md).

## Method

Step 1 is extraction and description only. The report definitions (return to the dry level, median curves,
dry-weather rises) are in the step 1 report (diag: rerun `step1_extract_tails.py`) and in the `waa.report` settings;
the decisions of every step are in [`DECISIONS.md`](DECISIONS.md). The WAA model method will
go to `METHOD.md` once it is decided.

## Open issues

- **Step 4 found no drying law on the real tails**: the free curve falls with W because about one tail in five
  stays flat or rises after t0_trace (mostly in the events flagged event_rise_heavy). To decide before Step 6:
  which tails define "drying"; and two repairs of the tools: `waa_pooled.py` lets F1 / F4 fit a rising tail with a
  negative start W0 while F3 / F5 / F6 and the free curve cannot (the menu BIC is not like for like), and the
  residuals of drying tails have 4 times the variance of the dry-weather noise model (BIC margins overstated).

- **The pooled likelihood prefers a false feature near 0.2-0.4 dB (3b-3c)**: a step is drawn for every fake
  truth, also the exponential, and wins on held-out events. The 0.3 dB quantization is NOT the cause (3b-4: the
  step is kept in 5 of 5 also without it, and W_c stays at its start value); the cause is open (suspects: the
  extra search budget of the step fit, the offset prior / start floor, the narrow step shape). Until found: no
  step tool, no jump 3 by held-out likelihood, nothing said below 0.35 dB. Also open: the event-bootstrap band of one dataset is about 3 times too narrow.

- **The free curve cannot tell a switch from two layers (3b-3b)**: above 0.35 dB both give a mildly rising curve
  (slope 0.2-0.5) with a multi-term equation, and "free curve worse than F4 by > 10 BIC" fired in only 2 of 5 on
  fresh data. Step 4 states its result for W >= 0.35 dB and the first 60 min only. Open: was the F4 fit of 3b-3b
  (3 generic starts) at its optimum?

- **The free drying curve is right only above about 0.35 dB (3b-3)**: below it the first spline knot is badly
  determined (and the likelihood without quantization may favour fast drying there); the per-fit interval from the
  curvature does not work (3 of 65). Open before step 4: repair interval, knots and search, or read the curve only
  above the second knot?

- **SINDy does not find the drying law on fake tails (3b-2)**: one pooled law dA/dt = h(A) returns "A" for the
  quadratic, switch and two-layer truths, and a false A^2 for exponential tails of mixed speed. Step 5 (SINDy on real
  tails) should not be run as planned; a redesign (speed per tail inside the regression, offset prior, noise in the
  library terms) has to pass the same fake-data test first. Open: keep step 5, redesign it, or drop it?

- **What fitting can show (3b-1, 3b-1b, 3b-1c)**: the per-tail fit (3b-1) separates only "exponential or not"; the pooled
  fit (3b-1b) separates all five forms on simulated data. The stress test (3b-1c) confirms this for exponential, two
  layers and the state switch under unknown noise, honest starts and faster truths, but F3 and F5 cannot be told
  apart (one family, "nonlinear in W"), F5's linear rate is not recoverable from the first hour, and F6 needs several
  optimiser starts. Open: is keeping m4's tau_s and m6's k_evap fixed in the T3 recalibration acceptable?

- **m1 with spread mimics m4 (3a-fix)**: a mix of exponentials gives the same long slow tail in the median as two
  layers. Whether the data can tell them apart is the question of 3b (per-tail shape, not the median).
- **Spread vs noise**: the spread s = 1.25 is set from the band after the noise is fixed; any noise the model still
  misses (several time scales, daily cycle) ends up in s.

- **Rise cut v2 vs night-time rises (v1)**: in E20150701_0030 A rises 0.5-1 dB after t0_trace with no radar or gauge
  trace, then decays; v1 cuts these. Missed drizzle (radar 78 km away) or condensation on wet radomes? RH at GbgA is
  missing for those hours (102 June hours missing in total), so this cannot be decided here.
- **Small v1 set**: 429 TRAIN tails from 24 events, but 3 events give 246 of them (E20150804_1830 115,
  E20150826_0520 75, E20150628_0530 56). Resample or weight by EVENT in the model step.

- Step 2b moved the start to t0_trace. 381 of the 959 step 2 main-set tails then fall below A_start = 1 dB: much of
  the step 2 "fast drop" happened while the radar or a gauge still showed a rain trace.
- The strict and lenient main sets hold the same 479 tails (the strict end shortens tails, it does not remove them);
  220 strict tails end at a rain trace.
- `undershoot` mostly marks a lasting level near -0.3 dB right after the rain (baseline offset for that event?), not a
  short dip.
- 104 tails have a **wet** t0 (link still wet when the network event ended; `last_wet_end` is searched only inside
  the event), and 566 tails of the last event were clipped to the season end. All have length 0.
- Dry-weather "rises" at 0.3 dB are frequent (56 % of tails) and partly RSL quantization. A firmer definition is needed
  before they are used.
- Light rain after t_first_dry: the step 2b strict end cuts at the first trace; the lenient end keeps it. Which
  end the WAA model uses (or rain as an input) is still open.
- The L test of the fast drop has little power: the path-rain reference `wet_A_max` itself barely grows with L
  (Spearman 0.13 for >= 25 GHz). See DECISIONS.md, step 2.
- `wet_A_max` spans the link's whole wet phase in the network event (several showers for multi-run events).
