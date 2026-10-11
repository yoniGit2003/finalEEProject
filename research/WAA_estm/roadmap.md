# WAA_estm roadmap

The one place for the WAA plan. Update it when a step finishes or the plan changes. What each finished step decided,
with numbers, is in `DECISIONS.md`; rules of the frozen tail set are in `README.md`. Status checked against the repo on 2026-10-09.

**Goal.** A dynamic WAA law, `dW/dt = g(R, W) - h(W)`.
- `h` = drying, learned from drying tails (R = 0).
- `g` = wetting, learned from very short links, where path rain is small.
- Compare with the literature: Minda & Nakamura exponential drying (`dW/dt = -W/tau`), Schleiss wetting (rise toward A_max).

**Standing rules**
- Uncertainty by resampling EVENTS; leave-one-event-out. 3 events give 246 of the 429 train tails.
- Sensitivity runs at every modelling step: v1 main set; without `event_rise_heavy`; without `level_shift` / `undershoot`; alternative set (A_start >= 0.6 dB).
- Tools are tested on simulated tails (known truth) before they are used on real tails.
- No tail is chosen by how well it fits a model (circular).
- Every output is saved with an importance level (key / product / diag); `/tidy` deletes diag and scratch.

## Status

| Step | What | Status |
|---|---|---|
| 0 | Inventory of inputs | done |
| 1 | Extract all drying tails (20,355) | done |
| 2 | Selection and fast-drop test | done, see correction below |
| 2b | Fixes: t0_trace, strict/lenient end, noise exclusion, level flags | done |
| 2c | Frozen tail set v1: 429 train / 123 val tails, 159 links, 24 events | done |
| 3a | Tail simulator calibrated to v1 | done |
| 3a-fix | Slow noise, spread for every law (s = 1.25) | done |
| 3b-1 | Fit candidate laws to fake tails, confusion matrix | done, see result below |
| 3b-1b | Pooled fit with offset prior; add F6 switch; event-bootstrap replicates | done: pre-set rule met in all 3 windows (own law 80-100 % for m1, m4, m5, m6; parameters within +-15 %); m3 vs m5 partly confused |
| 3b-1c | Stress test of the pooled fit: unknown noise, honest starts, faster truths | done 2026-10-11: PARTIAL (12 of 14 cells; m5 fails: F3 / F5 not separable, F5's a not recoverable) |
| 3b-2 | SINDy on fake tails | next |
| 3c | Simulator for full events (rain, wetting, path term) | not started |
| 4 | Fit candidate drying laws to real tails | not started |
| 5 | SINDy on real tails: drying law h(W) | not started |
| 6 | What drying depends on | not started |
| 7 | Full law from very short links (g and h together) | not started |
| 8a-8e | Verification | not started |
| 9 | Final result | not started |

## Corrections found when checking the repo (DECISIONS.md vs what I said in chat)

1. **Fast drop is not proven to be WAA.** Step 2 found it does not grow with L, but the path-rain reference `wet_A_max`
   barely grows with L either (Spearman 0.13), so the L test has little power. Step 2b then showed that 381 of the 959 step 2 main tails
   fall below 1 dB at `t0_trace`: much of the step 2 fast drop happened while radar or gauge still showed rain. Only the
   sub-link agreement (rho 0.87) and the flat control curve support WAA; both are weaker evidence than I claimed.
   -> The very-short-link length check in Step 7 is therefore the real test of "WAA, not path rain", with enough contrast in L.
2. **Three events, not two, give most tails** (246 of 429).
3. **With spread, a plain exponential (m1) reproduces the slow tail as well as two layers (m4).** The median shape does not favour m4.

## Part 2: test the tools on fake data

**3a / 3a-fix (done).** Simulator on the v1 skeleton. Noise = 10 s white + 1 min AR(1) wander + 0.3 dB quantization + per-tail offset.
Laws: m1 exponential, m3 `-c W^2`, m4 two layers, m5 `-aW - bW^2`, each with a speed spread s = 1.25. m2 was dropped (= m1 with spread).

**3b-1 (done).** Fit F1 exponential, F4 two layers, F3 nonlinear, F5 linear+quadratic to each fake tail; BIC with an
effective sample size for slow noise. Outputs: confusion matrix (true law x chosen model), parameter recovery,
dA/dt vs A, t_half vs A_start with replicate ranges, ablation (noise, quantization, offset, spread off).
Optional truth m6 (state switch) is calibrated. F6 (switch fit) was added in 3b-1b.

**3b-1b (done).** Pooled fit per dataset (`waa_pooled.py`): shared law, per-tail start and speed factor (marginalised),
offset prior, white + AR(1) noise in the likelihood, F6 added. On 10 event-bootstrap replicates per truth the pre-set
rule is met in all three windows (60 min, 6 h, whole); the first 60 min is the best window.

**3b-1c (done, PARTIAL).** Tried to break the pooled fit on fake data (first 60 min, 5 replicates): T1 noise with two
slow time scales, T2 generic optimiser starts, T3 truths recalibrated to the real early decay (+ s = 1.5). m1, m4, m6
pass every test. m5 fails twice: with the unknown noise F3 wins 5/5 (by 1.8 BIC), and with the faster truth F5's
linear rate a is x0.65. The failing part is understood (m5 is almost a pure quadratic, which F3 copies), so 3b-2 is next.
Carried to Step 4: several starts per candidate (F6 is start-sensitive in 20/20 datasets); F3 and F5 are one family;
F4 vs F6 by < ~10 BIC is not decisive.

**3b-2.** SINDy on fake tails. Library: 1, A, A^2, A^3, switch term max(A - A_c, 0) for several A_c.
Checks: does it return the true law for m1, m3, m5, m6? What does it return for the hidden-state m4? Does mixing tails with
different tau create a false nonlinear term?

**3c.** Extend the simulator to full events: real radar rain series on short links, wetting law (Schleiss form), drying law,
a small path term that grows with L, noise. Test SINDy with rain as input.

## Part 3: discover the law

- **4.** Fit candidate drying laws on real tails (literature exponential, two layers, nonlinear, switch). Real dA/dt vs A plot. Run the four sensitivity sets. From 3b-1c: several optimiser starts per candidate; F3 / F5 reported as one family.
- **5.** SINDy on real tails -> drying law h(W). Stability under event resampling.
- **6.** What drying depends on: tau vs frequency, antenna size, time of day, temperature, wind, humidity (RH is missing for part of June).
- **7.** Full WAA law from very short links (< ~0.6 km), whole events.
  - Path rain must be small without assuming the power law: among short links, attenuation during rain must not grow with L; both directions agree; restrict to light/moderate rain if heavy rain shows an L effect.
  - Keep h(W) fixed; learn g(R, W) with SINDy, radar rain along the path as input; also fit Schleiss.
  - Support: WAA at tail start vs event rain amount and duration (saturation, A_max).

## Part 4: verification (after discovery)

- **8a.** Drying on unseen data (VAL events, leave-one-event-out); our law vs the literature exponential.
- **8b.** Full law on unseen events (short links); vs Schleiss.
- **8c.** Transfer to longer links: A minus predicted WAA should grow with L with no constant offset (the f = L + L0 test). Handoff to Stage 2.
- **8d.** Consistency: SINDy tau vs directly fitted tau; sub-link pairs; stability under event resampling.
- **8e.** Dew events should dry by the same h(W).

## Part 5: final (9)

Full WAA equation (wetting + drying), parameters and their spread, comparison with the literature, limits, input to Stage 2.

## Options kept open (not forced)

- m6 / F6 state-switching law: fast runoff above a threshold A_c, slow evaporation below. Tested in 3b-1 (as fake truth), 3b-2, 4, 5.

## Open questions

- Night-time rises without any rain trace (E20150701_0030): drizzle under the radar beam, or condensation? RH is missing for those hours.
- Test events and test links have no baseline yet (needed for the final held-out check).
- `ORGANIZATION.md` / `check_organization.py`: confirm `roadmap.md` at the module root is allowed.
