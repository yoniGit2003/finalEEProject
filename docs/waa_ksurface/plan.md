# WAA k-surface branch: plan and standing rules

## THE CORE IDEA (re-read before every step; added 2026-11-11, Block A step A0)

1. The target F(R) (steady WAA) and the path loss P(R) are MEASURED, model-free, from how the attenuation depends on link
   length.
2. The dynamics are dW/dt = k(R, W)·(a_link·F(R) − W), k > 0. This is a rewrite of any one-state law with steady state F,
   not an assumed form. k is learned from data, separately for wetting (W below target) and drying (W above target).
3. The learned k is written as an equation and compared with Schleiss on held-out storms.

Only three questions exist:

- **Q1** Can F and P be trusted?
- **Q2** What is k?
- **Q3** Does the equation beat Schleiss?

Before each step, one sentence in `decisions.md` says which question the step serves. An analysis that serves none of them
is not run: it is logged as "side question, not pursued". After a failed check, at most ONE diagnostic to locate the cause,
then stop.

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
  events, TRAIN / TEST links). The list of these files and their paths is in `decisions.md`, Step 0.
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
