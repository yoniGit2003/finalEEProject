# WAA k-surface branch: plan and standing rules

## THE CORE IDEA (re-read before every step; rewritten 2026-11-11 by the Detachment correction, Deviation 2)

1. The steady WAA F(R) and the path loss per km P(R) are MEASURED, model-free, from how the attenuation depends on link
   length.
2. The dynamics are dW/dt = k(R, W)·(F(R) − W), k > 0, with k learned from data separately for wetting (W below F) and
   drying (W above F).
3. At the end, the equation is compared with the existing models on held-out storms.

Only three questions exist:

- **Q1** Can F and P be trusted?
- **Q2** What is k?
- **Q3** Does the equation beat the existing models in the final comparison?

Before each step, one sentence in `decisions.md` says which question the step serves. An analysis that serves none of them
is not run: it is logged as "side question, not pursued". After a failed check, at most ONE diagnostic to locate the cause,
then stop.

## THE DETACHMENT RULE (2026-11-11, Deviation 2; it overrides every earlier rule that conflicts with it)

This branch is a FRESH research effort. It inherits nothing from Part 1 (the `research/wet_antenna` work: its models,
parameters, thresholds, definitions, closed forms or findings) except the whitelist below. Part 1 is treated like a published
paper by another group: it is a benchmark at the end, never a foundation. No number is copied from it and no definition is
adopted silently.

Whitelist:

- **W1. The data layer and its quality rules:** the OpenMRG release tables, the baseline and the offset c, the event and link
  splits, the link-event selection rules, and the rain input R = R_radar / g_e. Each use is logged as "inherited data layer".
- **W2. The 0.3 dB RSL quantization step,** as an instrument property.
- **W3. Part 1's stored frozen fits, ONLY as competitors in the final Block B comparison.** The no-WAA and no-memory
  baselines are not Part 1 property: they are defined by recipe in this branch and fitted on TRAIN at Block B time.
- **W4. Part 1's documented warnings, as things to TEST, not to assume:** the 5-min input distorts fast dynamics; some events
  hold unseen water after rain (waiting / rising tails); dew.

Everything else is out, in particular: a_link and its prior (a per-link factor may enter only if this branch's own
measurements demand it); Part 1's time constants, stores, closed forms, lags, rain classes, steadiness definition,
cell-eligibility rule, noise models and tolerances. Where a rule copied such a thing, it is redefined from this branch's
needs and justified from the physics, the instrument, or this branch's own data.

**Amendments to the frozen rules below (the frozen text is kept as written on 2026-11-11; these amendments override it):**

| rule | inherited item | amendment |
|---|---|---|
| core idea, A3f, A5, B6 | a_link scaling the target | removed. The model is dW/dt = k(R, W)·(F(R) − W). A per-link factor enters only if step E4's pre-registered test of the per-link residual spread demands it. |
| A1b | reproduction of the Part 1 length test as a gate | stands as a one-off data-loading sanity check, diagnostic only. It is not a foundation of any later number. |
| A3a, A1e | "the Part 1 classes", "the Part 1 eligibility rule", the Part 1 steady-rain rule | re-derived in step E3 (rain classes, steady rain, eligible cell), each with its reason, and tested on fake data. The 32 GHz exclusion is re-decided with the new eligibility rule. |
| A3d | fake truth = the Part 1 closed forms, a_link spread 0.19 | replaced in E2 / E3 by a generic truth with round numbers and by this branch's own measured noise and link spread. |
| A5 | "A − c = W + L·P(R(t − 1 min))": the 1 min lag | the lag is not inherited. It is measured from this branch's own data before A5 (or set to 0 if it cannot be resolved at 5 min). |
| A5 | "Loss as in Part 1: Huber 0.3 dB ..., weak prior on ln a_link" | Huber loss on the 5-min residuals with the threshold 0.3 dB justified by the instrument (W2: one RSL quantization step); no a_link, no prior. |
| A5 | CV "stratified by the Part 1 intensity classes" | stratified by the `intensity_class` column of the event table (inherited data layer, W1). |
| A5 | unit test against "the Part 1 time constants" and "the Part 1 reference implementation's fast store" | unit test against the closed-form solution of dW/dt = k·(F − W) for a constant k and a step in R, with round-number constants, within 1e-6 dB. |
| A6 | "Part 1 noise" | the noise is measured in this branch from the dry minutes of the TRAIN link-events. |
| A6 | truth (a) "constant k (2.6 and 10.3 min)"; truth (d) "the full Part 1 two-store equation"; truth (e) "as in research-log Section 8" | (a) a constant k with round-number time constants; (d) a generic two-store truth with round-number constants; (e) stays, as a W4 warning to test, with the unseen water defined in this branch. |
| B2 | "the events flagged for rising tails and the dew-risk link-events" | the rising-tail flag is defined by this branch's own criterion before use (W4: tested, not assumed); the dew-risk hours come from the baseline module (inherited data layer, W1). |
| B6 | "no memory and no WAA (their stored TRAIN fits)" | defined by recipe in this branch and fitted on TRAIN at Block B time (W3). Part 1's stored fits of its equation, Schleiss (b) and Pastorek (b) are competitors only. |

Branch `waa-ksurface`, created 2026-11-11 from `main` at `fde431c`. This branch is a side study of the WAA research. Part 1 (the
two-store WAA equation, `research/wet_antenna/FINAL.md`) is frozen, and its TEST set has already been read once. This branch
does not change Part 1.

Branch folders (the repo convention; changed in Step 0b, see Deviation 1 in `decisions.md`):

| folder | holds | in git |
|---|---|---|
| `research/waa_ksurface/` | the code of the branch | yes |
| `research/waa_ksurface/outputs/` | small tables and figures | yes |
| `research/.cache/waa_ksurface/` | large arrays and fits | no (git-ignored) |
| `docs/waa_ksurface/` | this plan and `decisions.md` | yes |

Python environment: `C:\Users\yonic\miniconda3\envs\cmlpipe\python.exe` (the `cmlpipe` conda env). Every script of the
branch is run with it, from the repository root. The base miniconda Python lacks `numba` and `xarray`.

## Goal

One local, interpretable, dynamic WAA equation for a single link, learned from all links pooled, with no assumed power law and
no assumed dynamic form. The only assumptions are:

1. W is one state, driven by the rain.
2. Wetting is separated from drying.

F(R) (steady WAA) and P(R) (path loss per km) come from the model-free length test, not from a fitted model.

## Planned steps

| step | name |
|---|---|
| 1 | theory note |
| 2 | F and P from the length test |
| 3 | coverage of the (R, W) plane |
| 4 | build the model |
| 5 | fake-data recovery |
| 6 | real fit |
| 7 | one-state test |
| 8 | naming the equation |
| 9 | physical consistency |
| 10 | freeze and test |

## Standing rules (they apply to every step)

### Isolation
- Part 1 code, outputs, the frozen equation (`research/wet_antenna/outputs/products/waa_equation_final.json`) and its hashes
  (`fit_sha256` in that file, `research/wet_antenna/outputs/products/test_s8_freeze.json`) are read-only.
- Reuse the existing baseline (v2), offset c, rain input R = R_radar / g_e, link selection and event split (TRAIN / VAL / TEST
  events, TRAIN / TEST links). This is whitelist item W1 of the Detachment rule ("inherited data layer"). The list of these files
  and their paths is in `decisions.md`, Step 0.
- If anything must change, copy it into a branch folder and record in `decisions.md` what was copied and why.

### Justification
- Before every step, write an entry in `decisions.md` with: what the step does, why it is needed (a physical or statistical
  reason, never "to get a better result"), and a dated pass rule.
- After the run, add the result and the verdict to the same entry.
- A choice made after seeing a result is marked **DEVIATION**, with its reason.

### Fake data first
- Every method must recover a known truth on fake data before it is used on real data.

### Stop rules
- If a step fails its pass rule, stop and report. Do not patch it until it is approved.
- Stop after each step and wait for approval before starting the next one.

### Check report
- Every step ends with a CHECK REPORT in this format (the fields are adapted to the step):

```
=== CHECK REPORT: STEP n ===
Branch name and the commit it was created from:
Folders created:
Files created or changed (full paths):
Files outside the branch folders that were changed (must be "none"):
<the step's own checks, each answered yes / no or with numbers>
Anything done beyond the step's tasks (must be "none"):
Open questions for us:
=== END ===
```

---

# Block A and Block B rules (frozen 2026-11-11, step A0, before any analysis)

**Read with the Detachment rule at the top of this file (2026-11-11, Deviation 2): where a rule below conflicts with it,
the amendment table there overrides the text below. The text is kept as it was frozen.**

General rules of Block A: it runs unattended. Where a rule below decides, it is followed. Where no rule decides, or a rule
says STOP, the run stops and the report is written. A failed check is never patched to continue. Part 1 files are read-only.
Code is in `research/waa_ksurface/`, outputs in `research/waa_ksurface/outputs/`, large arrays in
`research/.cache/waa_ksurface/`. Every step has a `decisions.md` entry (what, why, dated pass rule) BEFORE it runs, and the
result and verdict after. There is a commit after every step. VAL and TEST attenuation data are NOT opened in Block A: TRAIN
only.

## A0. Write the rules first
All rules of Block A (A1–A6) and Block B (B1–B6) are written here, dated, and committed BEFORE any analysis (commit message
"waa-ksurface: rules frozen").

## A1. Data audit and reproduction gate
a) Rebuild the Part 1 TRAIN dataset with the Part 1 code (the fit cache is missing). PASS only if it gives exactly 206 links
   and 9 780 TRAIN link-events. Else STOP.
b) Reproduce the Part 1 TRAIN length test and compare with `outputs/tables/final_s5a_length_test.csv`. PASS if every
   intercept and slope matches within 0.01 dB (or dB/km). Else STOP. Write down the exact Part 1 definition of a "steady rain"
   5-min value and of an eligible cell.
c) Report, for TRAIN: link-events and rainy minutes per band and length class; minutes per rain class on links below 1 km;
   number of links and link-events that have a 1-min City gauge within 1 km (the gauge-input set); events with a gauge ratio
   vs g_e = 1.
d) Locate the stored TRAIN fits of Schleiss (b), Pastorek (b), the no-memory model and the Part 1 equation (parameters and
   per-link factors). Report path and whether each loads. Do not refit them.
e) 32 GHz rule: if 32 GHz has no eligible length-test cell, it is excluded from this branch for every model, and this is
   stated in every later table.

## A2. Theory note (no code)
Write `docs/waa_ksurface/theory.md`: the one-state assumption dW/dt = G(R, W); steady state G(R, F) = 0 and the sign
conditions; the factorization k·(F − W) with the exact conditions it needs (G differentiable in W, one stable zero) and the
degenerate case ∂G/∂W = 0 at F; the wetting/drying split and why k may jump at W = target; where a_link enters; constant k
gives the Part 1 FAST store, and the full Part 1 two-store model is NOT a special case; a list of what is assumed and what is
not.

## A3. F(R) and P(R) from the length test (TRAIN)
a) Classes: the Part 1 classes, plus 0.1–0.25 and 0.25–0.5 mm/h if they are eligible by the Part 1 eligibility rule.
   Intervals: 200 event-bootstrap draws.
b) Steady-state check. Recompute the intercept using only values with at least 5, 15, 30, 60 min of rain in the same class.
   If the paired bootstrap interval of (30 min − 5 min) excludes zero in more than a quarter of the eligible cells, F is taken
   from the longest duration that still leaves 8 links below 1 km per cell; otherwise from all steady values. Report the curve
   of intercept against duration either way.
c) Rain-input diagnostic. Compute intercepts separately on events with a gauge ratio and on g_e = 1 events. Report only. F
   from the gauge-ratio events is kept as a third sensitivity curve for Block B.
d) Errors-in-variables fake check. Fake data on the real links and rain: a known W that does not depend on L, a known P, and a
   rain error that is larger on short links (one 2 km pixel) than on long links (several pixels), sized from the
   radar-vs-gauge scatter. PASS if the length test returns the true intercept within 0.15 dB and the true slope within 20 % in
   every eligible cell. Else STOP.
e) Output: F and P per band as monotone interpolations through the cell values, with F(0) = P(0) = 0, intervals, and the valid
   range of R. State openly that below the lowest eligible class the curve is an interpolation to zero, not a measurement.
f) In Block B, F and P are FIXED at these values (a_link scales F). Sensitivity refits use F at the lower and upper interval
   ends and the curve of (c). F is never free.

## A4. Coverage of the (R, W) plane (TRAIN, links below 1 km)
Compute Ŵ = A − c − L·P(R) for this map only, never for fitting. Density of minutes in the (R, W) plane, split into wetting
and drying. The domain for learning and judging k is the set of cells with at least 200 link-minutes from at least 5 events.
Freeze it. Report separately the coverage of the gauge-input set.

## A5. Build the model (no real fit)
- A − c = W + L·P(R(t − 1 min)), dW/dt = k·(a_link·F(R) − W), W = 0 before the event, one exact exponential step per minute:
  W ← v + (W − v)·exp(−k·Δt).
- Loss as in Part 1: Huber 0.3 dB on the 5-min residuals, weak prior on ln a_link. Output error on whole events; the measured
  signal is never differentiated.
- The ladder for k, for wetting and drying separately: K0 constant; K1 depends on W; K2 depends on R; K3 = (function of
  R)·(function of W); K4 free surface on the A4 domain. Free parts are positive (exp of a spline) with a smoothness penalty
  chosen by leave-events-out CV inside TRAIN (5 folds, grouped by event, stratified by the Part 1 intensity classes).
- Rain timing: k_wet is learned with the 1-min gauge input on the gauge-input set and then held fixed; k_dry is learned with
  the radar input on all links.
- Unit tests: gradients against finite differences; with K0 and the Part 1 time constants and no slow store, the output
  equals the Part 1 reference implementation's fast store within 1e-6 dB.

## A6. Fake-data recovery (no real attenuation used as a target)
Fake events on the real TRAIN links and rain, 0.3 dB quantization, Part 1 noise, the model's F and P from A3. The attenuation
is generated from 1-min rain; the fit sees the same input it will see on real data (5-min radar-like input, 1-min on the
gauge-input set). Truths: (a) constant k (2.6 and 10.3 min); (b) k_dry rising with W; (c) k depending on R; (d) the full Part 1
two-store equation; (e) truth (a) plus unseen water after the rain in 20 % of the events (tails that wait or rise, as in
research-log Section 8).

Rules:
- (a)–(c): PASS if the ladder selects the true rung and the learned k is within 25 % of the truth on at least 80 % of the
  domain, weighted by data density. For (a), selecting any rung above K0 is a FAIL.
- If a rung above K0 cannot be recovered for wetting on the gauge-input set, wetting is restricted to K0 in Block B (a
  restriction, not a stop). Same for any rung that fails for drying: the ladder for Block B contains only recovered rungs. If
  even K0 fails, STOP.
- (d): report only: what the one-state ladder returns, and the after-rain residual by preceding rain, as the reference for the
  Block B one-state test.
- (e): report whether the fit returns a k_dry that falls with W although the truth is constant, and whether excluding the
  affected events restores a constant within 25 %.

## Block B rules (NOT executed in Block A)
- **B1.** Selection inside TRAIN only: the simplest recovered rung whose CV error differs from the best by a paired 90 %
  interval containing zero.
- **B2.** Guard: a dependence of k_dry on W counts as a result only if it keeps its sign and size (within its interval) when
  the events flagged for rising tails and the dew-risk link-events are excluded, and is consistent with the WAA measured on dry
  paths.
- **B3.** Second state: added only if, at matched W (±0.15 dB), the residual in the 30 min after rain differs between the
  lowest and highest third of preceding rain amount by more than 0.10 dB with a paired interval excluding zero.
- **B4.** Closed form: accepted if its CV error is within +0.005 dB of the free rung and the length test on its predictions
  agrees in at least two thirds of the cells (inclusive).
- **B5.** "New dynamics" is claimed only if the selected rung is above K0. If K0 is selected, the result is stated as "constant
  k confirmed with a measured target".
- **B6.** VAL is read once, after the equation is frozen and hashed, for the comparison with the Part 1 equation, Schleiss (b),
  Pastorek (b), no memory and no WAA (their stored TRAIN fits, same link-events, 1000 paired event-bootstrap draws): RMSE on all
  links and below 1 km, and the physical-consistency table. TEST is not opened.

## Final report of Block A
At the end of Block A, or at the first STOP, a CHECK REPORT: BLOCK A is written with the fields given in the Block A prompt
(stop, rules commit, A1a–A1e, A2, A3a–A3e, A4, A5, A6, deviations, side questions, files changed outside the branch folders,
VAL / TEST opened, commits).

---

# The A3 rescue: rules of steps E1–E4 (written 2026-11-11, after the STOP at A3d, before any rescue analysis)

The A3d failure means the length test, as computed, is biased by rain-input error. The open question is the size of the real
error. E1 measures it, E2 redesigns the test to survive it, E3 re-derives the bin definitions, E4 re-measures F and P. Rules
and pass criteria are written in `decisions.md` before each run. A failed check is never patched to continue.

## E1. Partition the rain error (TRAIN period, no attenuation data needed)
a) Triple collocation: at locations where radar, a City gauge and the SMHI gauge (or two City gauges) overlap, estimate each
   sensor's error variance without assuming any of them is truth. Matched 15-min or event accumulations; the matching rules
   are logged.
b) Gauge-pair analysis: correlation and scatter of City gauge pairs against their separation, extrapolated to 0–2 km, to size
   the point-vs-area representativeness part.
c) The binning-reversal signature: bin radar by gauge and gauge by radar; the asymmetry constrains the split.
d) Output: an interval [σ_lo, σ_hi] for the radar's own error at the 5-min, per-link scale, split into a common-per-event scale
   part and an independent per-bin part, and one paragraph on what the g_e = 1 intercept collapse (A3c) says about the scale
   part. If the data cannot constrain the split, this is said and the conservative range is given.

## E2. A length test that survives the measured error
Pre-register, implement and compare at least these candidates on fake data:
 (i) Gauge-ratio events only, with the A3d check rerun at σ_lo and σ_hi.
 (ii) Instrument binning: the bin is selected by a rain estimate whose error is shared by all links in the cell (event-gauge
      rain, or radar averaged over a wide common area), so the selection error does not depend on L; each link's own L stays
      in the regression.
 (iii) Forward-bias correction: simulate the bias at the measured error sizes and correct the intercepts, propagating
      [σ_lo, σ_hi] into the F interval.
The fake data use the E1 error model (common scale part per event + independent part), not all-noise-to-radar. PASS rule per
candidate: intercept within 0.15 dB and slope within 20 % in at least 13 of 14 cells, at BOTH σ_lo and σ_hi. The simplest
passing candidate is adopted; if several pass, the one with the narrowest F intervals on real data; if none passes, STOP and
report which cells fail and why.

## E3. Re-derive the bin definitions from scratch, each with its reason
What "steady rain" must mean so that the intercept measures F (W near its steady value, and no systematic drying into the
bin); the minimum links and lengths for a cell to be fitted; both tested in the E2 fake world (a bin rule that passes with a
known W(t) that has realistic memory: a generic one-store truth with round-number time constants).

## E4. Re-measure F(R) and P(R)
With the adopted method and definitions, per band, with event-bootstrap intervals AND the E1 error range propagated; monotone
interpolation, F(0) = P(0) = 0; the valid R range is stated, and that below the lowest measured class the curve is
interpolation. The per-link residual spread across storms is reported (the fresh test of whether a per-link factor is needed,
with its pre-registered criterion). Then STOP for review before A4. A4–A6 are not started.
