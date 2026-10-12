# WAA_estm roadmap

The one place for the WAA plan. Update it when a step finishes or the plan changes. What each finished step decided,
with numbers, is in `DECISIONS.md`; rules of the frozen tail set are in `README.md`. Status checked against the repo on 2026-10-11.

**Goal.** A dynamic WAA law, `dW/dt = g(R, W) - h(W)`.
- `h` = drying, learned from drying tails (R = 0).
- `g` = wetting, learned from very short links, where path rain is small.
- Compare with the literature: Minda & Nakamura exponential drying (`dW/dt = -W/tau`), Schleiss wetting (rise toward A_max).

## What the final result is

- One equation for the wet-antenna attenuation W (dB), per minute: `dW/dt = g(R, W) - h(W)`.
  h(W) = drying (no rain), learned from drying tails. g(R, W) = wetting, learned from very short links.
- With numbers: parameters with uncertainty, how much they vary between links, and what they depend on
  (frequency, antenna, weather).
- With plots of the drying and wetting curves, a comparison with the literature (Minda & Nakamura drying, Schleiss wetting),
  the prediction error on unseen events vs the literature, and the limits.
- Use: subtract the predicted WAA from a link's attenuation; what is left is rain attenuation (input to Stage 2).

## Method: how the equation is discovered from the data

The recipe is the "universal differential equation" approach of scientific machine learning (Rackauckas et al. 2020), as used for
population data with between-individual variability in Martensen et al. 2024 (DeepNLME + symbolic regression):
- Jump 1 (learn the shape): a flexible function inside the model where the law is unknown (spline for one input; small neural
  network for several inputs), fitted with our pooled likelihood (offset prior, AR(1) noise, per-tail speed = random effect).
- Jump 2 (name it): turn the learned function into a short formula. Candidates ordered from simple to complex.
- Jump 3 (refit and choose): put each candidate formula back into the full pooled model, refit on the tails, choose by BIC on
  held-out events. Only Jump 3 gives the final parameters and their uncertainty.

Every jump is tested on fake data with known truth before real data.

## Tools and where they are used

| Tool | What it does | Where |
|---|---|---|
| Free curve (spline) / small neural network | learns the shape from the data | Jump 1 in 3b-3, 3b-3c, 4 (spline); 6 and 7 (several inputs) |
| Best-subset fit + BIC | tries every combination of a few building blocks | Jump 2 for the drying curve (main) |
| PySR (symbolic regression) | searches any formula (powers, exp, fractions); simple-to-complex list | Jump 2 when the shape is not in our blocks; 7 |
| Ensemble SINDy (E-SINDy) | SINDy repeated over resampled events; keeps terms that keep coming back | cross-check in 4; main tool in 6 and 7 |
| Hybrid / switching idea (Mangan et al. 2019) | different equations in different W ranges | 3b-3c: Jump 1 may draw a step at a place it finds itself |
| SINDy with inputs (SINDYc) / parametric SINDy (SINDyCP) | equations driven by an input or with condition-dependent coefficients | 7 (rain as input); 6 (weather changes the drying speed) |
| Menu fits F1-F6 | tests our own candidate laws | comparison and literature check in 4 |

Note: pysindy and scikit-learn are not installed; pip install is allowed. PySR needs a Julia backend (installed by pip on first
use; may be slow). Say in each step which software was used.

**Standing rules**
- Uncertainty by resampling EVENTS; leave-one-event-out. 3 events give 246 of the 429 train tails.
- Sensitivity runs at every modelling step: v1 main set; without `event_rise_heavy`; without `level_shift` / `undershoot`; alternative set (A_start >= 0.6 dB).
- Tools are tested on simulated tails (known truth) before they are used on real tails.
- No tail is chosen by how well it fits a model (circular).
- Every output is saved with an importance level (key / product / diag); `/tidy` deletes diag and scratch.

## Status

| Part | Step | What | Status | Result |
|---|---|---|---|---|
| 1 Prepare data | 0 | Inventory of inputs | done | A(t) = TL - frozen baseline v2 |
| 1 | 1 | Extract all drying tails | done | 20,355 tails |
| 1 | 2 | Selection and fast-drop test | done | 959 main tails; fast drop not proven to be WAA |
| 1 | 2b | Fixes: t0_trace, strict/lenient ends, noise exclusion, level flags | done | |
| 1 | 2c | Frozen tail set v1 | done | 429 train / 123 val tails, 159 links, 24 events |
| 2 Practise on fake data | 3a, 3a-fix | Fake-tail simulator with realistic noise and spread | done | |
| 2 | 3b-1 | Fitting tool, one tail at a time | done | too weak: only exponential vs slower |
| 2 | 3b-1b | Pooled fitting tool (all tails together) | done | pre-set rule met |
| 2 | 3b-1c | Stress test of the pooled fit | done | PARTIAL: m1, m4, m6 hold; F3/F5 one family |
| 2 | 3b-2 | SINDy directly on raw tails | done | FAIL: returns "A" for almost every truth |
| 2 | 3b-3 | Two jumps: free-curve pooled fit -> drying curve, then SINDy -> equation | done | PARTIAL: curve and equation right above ~0.35 dB (mixing trap gone there), undetermined below; m1 fails the pre-set rule; the switch step and the two layers are not visible in the curve |
| 2 | 3b-3b | Confirmation run of the two jumps on fresh fake data, judged only above 0.35 dB, rule fixed before the run | done | PARTIAL: m1, m3, m5 PASS (curve, shape, equation); the m4 warning sign fails (2 of 5) |
| 2 | 3b-3c | Jump 2 on TRUE curves (seconds); Jump 1 that may draw a step at a fitted W_c; Jump 3 (refit and choose) | done | PARTIAL: Jump 2 alone works (E-SINDy, PySR); the step tool draws a step for EVERY truth; Jump 3 and the bootstrap band fail |
| 2 | 3b-4 | Diagnosis only: is the 0.3 dB RSL quantization the cause of the false step? m1 simulated with no / 0.3 / 0.1 dB quantization, step fit of 3b-3c part B | done | NOT the cause: without quantization the step is still kept in 5 of 5 (BIC gain 40-130); cause open |
| 3 Find the law on real data | 4 | The confirmed setup (3b-3b) + menu fits on the 429 real tails, four sensitivity sets | done 2026-10-12 | NO law identified: the real curve FALLS with W (slope -0.87; -0.05 without the rise-heavy events), no formula, menu not decisive (F1 as run; like for like F6 by 21 BIC); a W0 fault found in the menu comparison; F3 / F5 at their exponential limit in all 7 sets |
| 3 | 5 | SINDy on raw real tails | replaced | replaced by 3b-3 + 4 |
| 3 | 6 | What drying depends on (frequency, antenna, time of day, temperature, wind, humidity) | next (after the learners' decisions on Step 4) | |
| 3 | 3c | Fake full events with rain, to test the wetting tools | planned | moved here: only needed for Step 7 |
| 3 | 7 | Very short links: learn the wetting law g(R, W) with h(W) fixed | planned | the wetting law |
| 4 Prove it | 8a | Drying law predicts unseen tails (VAL, leave-one-event-out) vs literature exponential | planned | |
| 4 | 8b | Full law predicts unseen events (short links) vs Schleiss | planned | |
| 4 | 8c | Longer links: A minus predicted WAA behaves like pure rain (f = L + L0 test) | planned | |
| 4 | 8d | Consistency: curve vs fitted parameters, sub-link pairs, event resampling | planned | |
| 4 | 8e | Dew events dry by the same h(W) | planned | |
| 5 Final | 9 | Final equation, parameters, plots, limits, handoff to Stage 2 | planned | |

Order of work: 3b-3b -> 3b-3c -> 3b-4 -> 4 -> 6 -> 3c -> 7 -> 8a-8e -> 9. Plan changed on 2026-10-11 after 3b-2 (see Part 3).

## Plan change (2026-10-11, after 3b-2)

SINDy on raw tails failed the fake-data test (3b-2). Discovery is kept but split into two jumps: (1) the pooled fit with a
FREE smooth drying curve instead of a menu of laws draws the drying curve from all tails; (2) SINDy names the equation of
that clean curve. Old Step 5 (SINDy on raw real tails) is replaced by this; it becomes part of 3b-3 (fake) and Step 4 (real).
3c is moved to just before Step 7 (it is only needed for wetting). Step 7 keeps "learn g(R, W)"; the tool for it is decided
after 3c, by the same fake-data testing.

## Corrections found when checking the repo (DECISIONS.md vs what I said in chat)

1. **Fast drop is not proven to be WAA.** Step 2 found it does not grow with L, but the path-rain reference `wet_A_max`
   barely grows with L either (Spearman 0.13), so the L test has little power. Step 2b then showed that 381 of the 959 step 2 main tails
   fall below 1 dB at `t0_trace`: much of the step 2 fast drop happened while radar or gauge still showed rain. Only the
   sub-link agreement (rho 0.87) and the flat control curve support WAA; both are weaker evidence than I claimed.
   -> The very-short-link length check in Step 7 is therefore the real test of "WAA, not path rain", with enough contrast in L.
2. **Three events, not two, give most tails** (246 of 429).
3. **With spread, a plain exponential (m1) reproduces the slow tail as well as two layers (m4).** The median shape does not favour m4.

## Part 2: practise on fake data

**3a / 3a-fix (done).** Simulator on the v1 skeleton. Noise = 10 s white + 1 min AR(1) wander + 0.3 dB quantization + per-tail offset.
Laws: m1 exponential, m3 `-c W^2`, m4 two layers, m5 `-aW - bW^2`, each with a speed spread s = 1.25. m2 was dropped (= m1 with spread).

**3b-1 (done).** Fit F1 exponential, F4 two layers, F3 nonlinear, F5 linear+quadratic to each fake tail (free W0 and free offset c);
BIC with an effective sample size for slow noise. Truths m1, m3, m4, m5, m6 (m6 calibrated: k_run 0.044, k_evap 0.0055, A_c 0.16).
Result (10 replicates per truth):
- Per tail: no skill. F1 is chosen for 60-80% of tails whatever the truth.
- Per dataset (summed BIC): m1 -> F1 in 10/10; m3, m4, m5, m6 -> F4 in 8-10/10. Own candidate wins only for m1 (100%) and m4 (80%); never for m3, m5.
  So the method separates "plain exponential" from "slower than exponential", not which slow law. F4 winning does not mean two layers.
- Recovery: exponential tau unbiased (1.01); m4 tau_s x4.5 too large; m3 n at its bound for 40% of tails; m5 a, b not recoverable. Free offset c trades off with a long small slow layer.
- Ablation: no single effect (noise, quantization, offset, spread) is to blame; only "all off" recovers the truth.
- AIC picks F4 for every truth; use BIC.
- Population tests (t_half vs A_start, dA/dt vs A) do not separate the laws. t_half: real inside m1 range 1/4 bins, m3 3/4, m4 3/4, m5 4/4, m6 3/4
  (weak hint that wetter dries faster; replicate ranges hold noise only, not event-sampling uncertainty).
- Real dA/dt above ~1.5 dB is thin and zig-zags to ~0.

**3b-1b (done).** Raised the power of the fit (`waa_pooled.py`): (a) prior on the offset, c ~ N(0, offset_std^2), instead of a free c;
(b) AR(1) noise in the likelihood instead of the n_eff correction; (c) pooled fit per dataset: shared law parameters, per-tail
speed factor log-normal with a per-link part (marginalised); (d) F6 (switch fit) added; (e) replicates that resample events;
(f) first 60 min vs 6 h vs whole tail. Result: on 10 event-bootstrap replicates per truth the pre-set rule (own candidate >= 80 %
for m1, m4, m5, m6, key parameters within +-30 %) is met in all three windows; the first 60 min is the best window.

**3b-1c (done, PARTIAL).** Tried to break the pooled fit on fake data (first 60 min, 5 replicates): T1 noise with two
slow time scales, T2 generic optimiser starts, T3 truths recalibrated to the real early decay (+ s = 1.5). m1, m4, m6
pass every test. m5 fails twice: with the unknown noise F3 wins 5/5 (by 1.8 BIC), and with the faster truth F5's
linear rate a is x0.65. The failing part is understood (m5 is almost a pure quadratic, which F3 copies).
Carried to Step 4: several starts per candidate (F6 is start-sensitive in 20/20 datasets); F3 and F5 are one family;
F4 vs F6 by < ~10 BIC is not decisive.

**3b-2 (done; Q1 FAIL, Q2 FAIL, Q3 PARTIAL, Q4 FAIL).** Weak-form SINDy (own code, `waa_sindy.py`) on the 3b-1c fake
datasets. One pooled law dA/dt = h(A) returns "A" for the quadratic (m3, m5), switch (m6) and two-layer (m4) truths;
exponential tails of mixed speed get a false A^2 in 3 of 5 replicates (none without spread); it agrees with the pooled
fit in 1 of 20 datasets. With the true offset and no spread it does find A^2, so the offset and the speed spread are
what break it. Consequence: Step 5 cannot be run as planned (see the plan change above). Next: 3b-3. The original plan was:

**3b-2 (plan).** SINDy on fake tails. Library: 1, A, A^2, A^3, switch term max(A - A_c, 0) for several A_c.
Checks: does it return the true law for m1, m3, m5, m6? What does it return for the hidden-state m4? Does mixing tails with
different tau create a false nonlinear term?

**3b-3 (done, PARTIAL).** By the pre-set rule: m1 FAIL, m3 / m5 / m6 PARTIAL, m4 PASS (warning in 4 of 5). The failures sit
below about 0.35 dB, where the first knot (0.05 dB) is badly determined, and in a broken interval (available for 3 of 65
fits). Above the second knot, first 60 min, with spread (post hoc): m1 flat (slope -0.02, level x1.03, equation "W" 5/5),
m3 slope 1.08 with W^2 x0.94, m5 slope 0.96 with W^2 x1.06, m6 slope 0.27 (true 0.26): the mixing trap is gone there. The curve
cannot show m6's step at 0.16 dB nor m4's second state; the menu fits F6 / F4 are needed for those. Before Step 4: a working
interval (event bootstrap or profile), knots / range that do not leave 0.05-0.35 dB to one value, a longer search; and check
the suspicion that the likelihood without quantization favours fast drying at the low end. The plan was:

**3b-3 (plan).** "Draw the curve, then name it" on fake tails. Jump 1: the pooled fit of 3b-1b / 3b-1c with a free
smooth drying curve, dW/dt = -f_i W k(W), k(W) = exp(spline in log W), instead of a menu law. Jump 2: SINDy on that clean
curve h(W) = W k(W) to name its equation. Question: with a known truth, is the curve right and is the equation right?

**3b-3b (done, PARTIAL).** Fresh data, rule fixed before the run, judged range 0.35 dB to W_hi (1.9-2.3 dB), first 60 min.
m1, m3, m5 PASS on all three criteria: true k inside the replicate band on 100 % of the range, level within 3 %; slopes
of log k vs log W -0.01..0.05 (m1), 0.99..1.12 (m3), 0.93..1.03 (m5); one-term equation in 5 of 5 (W x0.99, W^2 x0.97,
W^2 x1.13); no false term for m1, so the mixing trap is gone. Same under T1 noise. m6 (switch below W_lo, not judged):
the curve is right (slope 0.20-0.32, true 0.23) but reads as "mildly rising", with a three-term equation. m4 FAIL: the
free curve is worse than F4 by > 10 BIC in only 2 of 5 (3b-3: 4 of 5), and its curve looks like m6's (slope 0.36-0.47).
So: the curve tells exponential from quadratic from "something else"; switch vs two layers stays with the menu fits, and
"free curve vs F4" is not a reliable check. The plan was:

**3b-3b (plan).** Confirmation run. The "right above 0.35 dB" result of 3b-3 was found after seeing the results, so it
is tested here on fresh fake data (seeds 4000 + replicate) with the rule written before the run: first 60 min only, 5 knots
from W_lo = 0.35 dB to W_hi, nothing judged or named below W_lo, band = min-max of the replicates (no curvature interval),
shape by the slope of log k vs log W in each replicate, jump 2 with switch levels 0.5 and 0.75 dB only.

**3b-3c (done, PARTIAL).** Rules fixed before the run; data = the 3b-3b datasets. A, Jump 2 alone on the true curve + 2-5 %
noise: E-SINDy (backward elimination of waa_sindy.py, bagged) exact in 20 of 20 for m1, m3, m5, m6; PySR 20 / 20 / 18 of 20
for m1 / m3 / m5 (cannot write the switch); best-subset + BIC always contains the true terms but adds one in 3-8 of 20;
pysindy's STLSQ fails at 5 % noise. B, step: kept in 5 of 5 for EVERY truth (BIC gain 15-170), also the exponential, at
0.2-0.4 dB: the tool cannot be used. C, Jump 3: the step law wins the held-out likelihood in 18 of 20 single-curve
datasets, by 16-65 over the true law, so the pooled likelihood itself prefers a false feature near 0.2-0.4 dB (suspected:
the 0.3 dB quantization is not in the likelihood; not tested). Menu laws alone with BIC: true family in 4-5 of 5.
D, band from an event bootstrap of one dataset: +-1 to 6 % wide, true k inside on 0-78 % (rule 80 %): about 3 times too
narrow; slopes within 0.01-0.11 of the truth. E: F4 with the F1-solution start changes nothing (m4 warning stays 2 of 5).
Software: pysindy and PySR (with Julia) are now installed. The plan was:

**3b-3c (plan).** (a) Jump 2 tested alone on the true k(W) of m1, m3, m5, m6 with a little noise added: best-subset + BIC,
PySR and SINDy must name the true formula (m6 included: is the switch named when the curve has the step?). (b) Jump 1 with
an optional step at an unknown W_c: must find W_c near 0.16 dB for m6 and "no step" for m1, m3, m5. (c) Jump 3: each candidate
put back into the pooled model and refitted; winner by BIC on held-out event draws. Rule fixed before the run, as always.

**3b-4 (done 2026-10-11, verdict "not the cause").** Diagnosis only, rule fixed before the run. Truth m1 with the 3b-3b
skeleton and seeds, simulated with no / 0.3 dB / 0.1 dB RSL quantization, free curve + optional step as in 3b-3c part B.
The step is kept in 5 of 5 in all three settings (BIC gain 40-130 without quantization), and the fitted W_c stays at its
start value (0.12, 0.2 or 0.4 dB) in 15 of 15 fits. So the 0.3 dB quantization is not the cause of the false step, and a
quantization-aware likelihood is not expected to repair it. Real data: RSL in 0.3 dB steps on 719 of 728 sub-links, 10 s
samples, 1 min means; the simulator matches. Cause still open (suspects: the extra search budget of the step fit, the
offset prior / start floor, the narrow step shape). Next = Step 4 with the reduced recipe of 3b-3c, unchanged.

## Part 3: find the law on real data

- **4 (done 2026-10-12: no drying law identified).** TRAIN v1, 429 tails, first 60 min, W_hi = 2.83 dB; 36 min, nothing
  cut. The free curve falls with W (slope -0.87; bootstrap refits -1.43..-0.82; sets -0.05 to -1.60), which no fake truth
  did: 21 % of the tails stay at >= 90 % of their start level for 10 min (38 % in the rise-heavy events) before they drop.
  Without the rise-heavy events the slope is -0.05 (exponential-like, tau 15 min). Jump 2 gives no formula (5 of 5 terms
  kept). Menu as run: F1 (tau 25.2 min), but F1 / F4 may fit a rising tail with a negative start (3 tails = the whole
  2600 BIC gap to F3 at n = 1); like for like F6 is ahead by 21 BIC, not decisive (residual variance 4.1 times the
  dry-noise model). Stable in all sets: F3 at n = 1, F5 at b = 0 (nothing rising with W). No literature value in the
  repo to compare tau with. Open before Step 6 (learners): which tails define "drying"; repair the W0 bound and re-test
  on fake data with rising tails; a noise level for drying tails. The recipe was: Only what passed on fake data is used: Jump 1 without a step, read for
  W >= 0.35 dB; shape = slope of log k vs log W (about +-0.1); Jump 2 by E-SINDy (waa_sindy elimination) with PySR as
  cross-check and best-subset as "true terms plus at most one extra"; NO Jump 3 by held-out likelihood; named laws compared
  by the menu fits with BIC; the event-bootstrap band is shown with the warning that it was about 3 times too narrow on fake
  data (level about +-20 %). Nothing is said about W < 0.35 dB or about a switch until the cause of the false step is found
  and B, C, D are repeated on the cached fake datasets (3b-4: the cause is NOT the 0.3 dB quantization). The earlier plan, kept for reference: The 429 real tails; every result is stated for W >= W_lo = 0.35 dB (from 3b-3b; up to W_hi)
  and the first 60 min only:
  (a) Jump 1: free-curve pooled fit (5 knots from 0.35 dB to W_hi, band from an event bootstrap) -> real drying curve and
  its slope of log k vs log W, judged only for W >= W_lo;
  (b) Jump 2: best-subset + BIC (main), PySR and E-SINDy as cross-checks (library W, W^2, W^3, switch at 0.5 / 0.75 dB) -> equation;
  (c) Jump 3: each candidate refitted in the pooled model, chosen by BIC on held-out events;
  (d) compare with the menu fits F1, F3, F4, F5, F6 under the 3b-1c rules (several starts, F3 / F5 one family, F4 vs F6
  < ~10 BIC not decisive) and with the literature exponential;
  (e) the four sensitivity sets.
- **5.** SINDy on raw real tails -> drying law h(W). **Replaced by 3b-3 + 4** (plan change 2026-10-11): as tested in 3b-2,
  SINDy on raw tails returns "A" whatever the true law.
- **6.** What drying depends on. Jump 1 learns drying as a function of W AND weather / link features (temperature, humidity,
  wind, time of day, frequency, antenna size); then E-SINDy / parametric SINDy picks which features matter, with stability over
  events and leave-one-event-out. First a fake-data test with a known weather effect and one with none (must not invent an
  effect). Caution: about 24 distinct events, 3 hold most tails; weather is hourly from one station; RH is missing for part of June.
- **3c (planned; moved here, only needed for Step 7).** Extend the simulator to full events: real radar rain series on short links, wetting law (Schleiss form), drying law,
  a small path term that grows with L, noise. Test the wetting tools of Step 7 (below) with rain as input.
- **7.** Full WAA law from very short links (< ~0.6 km), whole events.
  - Path rain must be small without assuming the power law: among short links, attenuation during rain must not grow with L; both directions agree; restrict to light/moderate rain if heavy rain shows an L effect.
  - Wetting as dW/dt = NN(R, W) - h(W), h fixed from Step 4 (universal differential equation), radar rain along the path as
    input; then PySR / E-SINDy on the learned NN -> g(R, W); Jump 3 refit. Fake-data test (3c) first; the tool is kept only
    if it passes there.
  - Compare with Schleiss (dynamic), Valtr (WAA = 2 k' R^alpha') and Kharadly-Ross (C (1 - exp(-d A))). Also test whether
    wetting needs time dependence at all (Pastorek et al. 2022 found the rain-rate models best and the dynamic Schleiss
    model adding little on their links).
  - Support: WAA at tail start vs event rain amount and duration (saturation, A_max).

## Part 4: verification (after discovery)

Final choice by BIC / prediction error on unseen events (VAL, leave-one-event-out), as in Martensen et al. 2024.

- **8a.** Drying on unseen data (VAL events, leave-one-event-out); our law vs the literature exponential.
- **8b.** Full law on unseen events (short links); vs Schleiss.
- **8c.** Transfer to longer links: A minus predicted WAA should grow with L with no constant offset (the f = L + L0 test). Handoff to Stage 2.
- **8d.** Consistency: SINDy tau vs directly fitted tau; sub-link pairs; stability under event resampling.
- **8e.** Dew events should dry by the same h(W).

## Part 5: final (9)

Full WAA equation (wetting + drying), parameters and their spread, comparison with the literature, limits, input to Stage 2.

## Naming: m and F

m = the law that MADE the fake data (the truth); F = the law the fitting tool tries as an answer. Same number, same law.
1 exponential (Minda & Nakamura) | 3 `-c W^n` | 4 two layers | 5 linear + quadratic | 6 state switch.
There is no 2 (m2 = m1 with spread, dropped in 3a-fix; not renumbered) and F2 never existed.
m3 was not a truth in 3b-1c; F3 still competes. F3 and F5 are read as one family ("nonlinear in W").

## Planned routine: /research-report

LaTeX technical report on the research part in the style of the pipeline report (pipeline_log.pdf): scope and sources,
overview, data and conventions, one section per step (Purpose / Method / Parameters / Results / Test and verification /
Issues found and corrections), claims-and-status table, limitations, chronology, references, glossary, an "In plain words"
box per chapter. Numbers injected from DECISIONS.md, outputs/tables and MANIFEST.csv by a script, never typed. Runs after
/tidy. Prompt prepared; run when the machine is free.

## Options kept open (not forced)

- m6 / F6 state-switching law: fast runoff above a threshold A_c, slow evaporation below. Tested in 3b-1 (as fake truth), 3b-2, 4, 5.
- Switching law (m6 / hybrid): tested in 3b-3c with a Jump 1 that may draw a step. Result: the step is drawn for every truth,
  so the tool says nothing yet. 3b-4: the 0.3 dB quantization is not the cause; the cause is open.

## Literature used for the method

- Rackauckas et al. 2020, Universal Differential Equations for Scientific Machine Learning, arXiv:2001.04385
- Martensen et al. 2024, Data-driven discovery of feedback mechanisms in AML (DeepNLME + symbolic regression), bioRxiv 10.1101/2024.06.17.599366
- Fasel et al. 2022, Ensemble-SINDy, arXiv:2111.10992
- Mangan et al. 2019, Model selection for hybrid dynamical systems via sparse regression, arXiv:1808.03251
- Cranmer 2023, Interpretable Machine Learning for Science with PySR and SymbolicRegression.jl, arXiv:2305.01582 (proposal ref. [7])
- Pastorek et al. 2022, Precipitation estimates from CMLs: practical approaches to wet-antenna correction (engrXiv; proposal ref. [9] is the IEEE TGRS version)
- Moroder / Chwala et al. 2023, ATTRRA2 field experiment: dynamics of WAA (drying after rain, dew before rain), IEEE / EGU23
  (the proposal cites the ATTRRA2 paper as ref. [2], Tiede, Chwala & Siart 2023: check the author list before citing)

Note: machine learning in CML work is mostly rain detection; discovering the WAA equation from data appears to be a gap
(not checked exhaustively).

## Open questions

- Night-time rises without any rain trace (E20150701_0030): drizzle under the radar beam, or condensation? RH is missing for those hours.
- Test events and test links have no baseline yet (needed for the final held-out check).
- `ORGANIZATION.md` / `check_organization.py`: confirm `roadmap.md` at the module root is allowed.
- Step 4: what keeps A flat or rising after the last rain trace in 21 % of the tails (rain / drizzle below the reference,
  re-wetting, condensation, path attenuation)? It decides whether those minutes belong to h(W) or to g(R, W).
- Step 4: is the dip of k(W) near the third knot (0.9 dB main set, 0.7 dB alternative set) made by the knots or by
  selecting tails by A_start?
- Step 4: drying time constants of the proposal's refs [1] and [9] are not in the repo; read them before Step 8a.
