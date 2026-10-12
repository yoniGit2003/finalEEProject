# DECISIONS: WAA_estm

What each step decided, with the numbers behind it, so that the step reports can stay local (/tidy moves them to `legacy/`, which is not in git) and be rebuilt on
demand (rerun the script named in each heading). Rules of the frozen tail set: README.md, "Tail set v1". Settings:
section `waa` of `research/config.yaml`. TRAIN events only for every number unless stated.

## Step 0, inventory (`step0_inventory.py`, 2026-10-08)

- A(t) is **not stored**: it is `tl_mean` (release) - `baseline_db` (`baselines.csv`), exactly `baseline.a_of_t`.
- `baselines.csv` covers train + val events on the 283 train links (566 sub-links); test has no baseline, so there are
  no test tails. 7375 significant train pairs have a wet bin and `use_in_waa_main`.
- link_state is per link, 5 min; the ±30 min margin is on the dry side only (unknown bins, not wet). So t0
  (`last_wet_end`) is normally an unknown margin bin and the first dry bin comes about 30 min later.
- QC flags that block a sample are already excluded from `tl_mean`; `tsl_step` and `flatline` are information only.
- Decision: build the tails on 1 min `tl_mean` with the frozen baseline v2; no new baseline in this module.

## Step 1, all drying tails (`step1_extract_tails.py`, 2026-10-08)

- Tail = (event, sub-link) from t0 to the earliest of next wet bin / NaN run ≥ 15 min / t0 + 12 h. **No selection.**
- 20 355 tails (19 659 non-empty) on 566 sub-links, 69 events; train 15 370, val 4 985. A(t) equals `a_of_t`
  exactly on 3 random pairs.
- Return to the dry level (A < max(0.3 dB, 2 offset_std)): A_start > 2 dB tails 25 % back at 10 min, 59 % at 30 min,
  79 % at 3 h; A_start < 0.5 dB tails 94 % at 10 min. Strong tails dry slowly: WAA-like memory exists.
- 670 tails start on a wet bin (`last_wet_end` is searched only inside the network event); 566 are clipped at the
  season end. Both kept as rows, excluded later.
- Dry-weather rises ≥ 0.3 dB later in the tail: 56 % of tails, partly RSL quantization (upper bound).
- Decision: `tails_1min.parquet` (47 MB) is git-ignored and rebuilt by this script.

## Step 2, selection and fast-drop test (`step2_select_tails.py`, 2026-10-08)

- Rules fixed before the run (`waa.selection`); no tail kept or dropped by its shape. VAL not looked at.
- Excluded 4023 of 15 370 (mostly `use_in_waa_main` False: 3421). Truncation: dry-weather rise (≥ 0.6 dB for ≥ 30
  min) cut 1249 tails (5.5 % of tail time), a TSL step 29.
- Main set (A_start ≥ 1 dB, ≥ 60 min usable): 959 tails, 373 sub-links, 33 events; control (< 0.5 dB) 8838.
- Fast drop (A_start - A at the first dry bin): median 1.11 dB = 83 % of A_start. It does **not** grow with L
  (Spearman ≥ 25 GHz -0.10 [-0.20, 0.00]; < 25 GHz 0.14 [-0.09, 0.31]), but the path-rain reference wet_A_max barely
  grows with L either (0.13), so the L test has little power. Sub-links 1 vs 2 agree (rho 0.87).
- Open after step 2: light rain after the first dry bin stays in 466 of 959 main tails -> fixed in step 2b.
- The two fast-drop outputs of step 2 (the L figure and the correlation table) were key until 2026-10-09, then
  demoted to diag (key limit; superseded by step 2b). Their numbers are above.

## Step 2b, fixes after the step 2 gallery (`step2b_fix_tails.py`, 2026-10-08)

- Start moved to **t0_trace** = end of the last 5 min bin with any rain trace (radar R_max > 0 on the path, or a gauge
  within 1 km > 0) between t0 and the first dry bin.
- 381 of the 959 step 2 main tails fall below A_start = 1 dB at t0_trace: much of the step 2 "fast drop" happened
  while the radar or a gauge still showed rain. So the fast drop is not evidence for WAA.
- Two ends: lenient (step 2 ends, rise cut moved back to the last quiet minute, excess < 0.15 dB) and strict (also the
  first rain trace after t0_trace). Same 479 tails in both; 220 strict tails end at a rain trace.
- Noise exclusion: sub-links with dry-weather `offset_std` > 0.3 dB (one RSL step): 39 sub-links, 1077 tails.
- `level_shift` (|A| at 6-12 h > max(0.3, 2 offset_std)) and `undershoot` (A < -0.3 dB within 3 h) are flags, not
  exclusions; undershoot looks like a lasting level near -0.3 dB (baseline offset for that event?), not a dip.

## Step 2c, frozen tail set v1 (`step2c_freeze.py`, 2026-10-08; its report is a key output)

- Radar NaN on the path counts as a rain trace (unknown); rise cut v2 on all minutes from t0_trace (10 min median ≥ 0.6
  dB above its running minimum for ≥ 10 min).
- **Frozen set v1**: main 429 train / 123 val tails; alternative (A_start ≥ 0.6 dB) 1155 / 458.
- 3 events give 246 of the 429 train tails: resample or weight by event in the model step.
- Open: night-time rises without any rain trace (E20150701_0030): drizzle under the radar beam or condensation? RH is
  missing for those hours, so this is not decided.

## Step 3a, tail simulator (`step3a_simulator.py` + library `waa_sim.py`, 2026-10-08)

- Skeleton = the 429 TRAIN v1 tails (same sub-links, links, events, usable lengths, A_start_trace); 1 min output.
- Laws calibrated by least squares on the real median of A / A_start over 0-6 h (noise-free truth, same survivors):
  m1 tau = 25.7 min; m2 median tau = 30.7 min, s = 1.5; m3 (n = 2) c = 0.064; m4 tau_f = 16.1, tau_s = 161 min,
  phi = 0.17; m5 a = 0.0036, b = 0.055. Median misfit (truth): m4 0.016 < m5 0.024 < m3 0.027 < m1 = m2 0.052.
- m2's spread s is not identifiable from the median (median of exp(-t/tau) = exp(-t/median tau)); s was taken from
  the real 10-90 % band (deviation from "median only", config `waa.sim.m2_s_from_band`). s = 1.5 lumps all
  tail-to-tail variation; only m2 matches the band (0.05 vs 0.15-0.18 for the deterministic laws).
- Noise: Gaussian 10 s sigma per sub-link from the step 2b control tails (median 0.09 dB). It matches the share of
  minutes on a 0.3 dB level but NOT the within-tail std (0.05 vs 0.14 dB real; last-hour std 0.057 vs 0.083 dB).
  Real dry noise has a slow part the model lacks, so 3b with this model overstates how well laws can be separated.
- Decision: products `sim_tails_m1..m5.parquet`, `sim_tails_index.csv` and a noise sigma table (both replaced in
  3a-fix); switches for offset /
  noise / quantization / real NaN mask in `waa.sim.switches` (all on here). Open: add a slow-noise term before 3b?

## Step 3a-fix, slow noise and a spread for every law (`step3a_simulator.py` + `waa_sim.py`, 2026-10-08)

- Noise = 10 s white + 1 min AR(1) wander + 0.3 dB quantization, per sub-link from the control tails (lookup blocks =
  the median control tail length, 680 min). All 5 targets now match (real / sim medians): within-tail std 0.141 /
  0.134 dB, on-level share 0.25 / 0.27, ACF 1 min 0.80 / 0.80, 10 min 0.57 / 0.63, 60 min 0.22 / 0.23. Fitted
  sigma_white 0.09, sigma_slow 0.15 dB, tau_noise 80 min (medians). Last-hour std now 0.083 real vs 0.088 sim.
- m2 dropped (= m1 with spread). Every law gets a time-scale spread f = exp(s z) (half link, half tail); base
  parameters from the median, s from the band: s = 1.25 for all four laws (f 10-90 %: 0.20-4.96).
  m1 tau = 29.8 min; m3 c = 0.053; m4 tau_f = 18.7, tau_s = 193 min, phi = 0.17; m5 a = 0.0030, b = 0.045.
- Observed misfit (median / band): m1 0.026 / 0.062, m4 0.033 / 0.062, m5 0.038 / 0.069, m3 0.044 / 0.072.
  **With spread, plain m1 reproduces the long slow tail as well as m4**: a mix of exponentials looks like two
  layers in the median, so the median shape alone does not favour m4.
- t_half per A_start bin (1-1.5 / 1.5-2 / 2-3 / > 3 dB): real 16 / 17 / 14 / 8 min; the real CIs overlap every law and
  the > 3 dB bin has 9 tails, so this check cannot separate exponential from nonlinear laws on one simulation.
- Retired (level scratch, moved to `legacy/` by the next /tidy): the m2 simulated set, the old noise sigma table, the step 3a
  report. One path constant per simulated set (`SIM_TAILS_M1` ...), so /tidy sees that they are still written.

## Step 3b-1, the fitting method on simulated tails (`step3b1_fitting_on_sim.py` + `waa_fit.py`, 2026-10-09)

- Truths m1, m3, m4, m5 (3a-fix) + m6 (state switch: run-off above A_c, evaporation; calibrated like the others:
  k_run, k_evap, A_c from the median, s = 1.25 from the band; converged). 10 replicates per truth; every tail fitted
  with F1, F3, F4, F5 (free W0 and offset c); BIC with one n_eff per tail (lag-1 rho of the F4 residuals; median
  n 210, rho 0.69, n_eff 64).
- **Per tail nothing is identifiable**: F1 wins 60-80 % of the tails for every truth.
- **Per dataset (summed BIC) only "exponential vs not" is identifiable**: m1 -> F1 in 10/10; m3, m4, m5, m6 -> F4 in
  80-100 %. m3 and m5 never get their own candidate. So F4 winning on real data would mean "not a single exponential",
  NOT "two layers".
- AIC picks F4 for every truth (m1 too): use BIC with n_eff.
- Recovery: m1 tau unbiased (fitted / true 1.01, 10-90 % 0.63-2.1); m4 tau_s not recoverable (x4.5, up to x266),
  it trades off with the offset c; m3 n and m5 a, b not recoverable.
- Ablation (m4): no single effect explains the confusion (own share 22-28 % with any one of noise / quantization /
  offset / spread off, 100 % with all off): noise alone or quantization alone hides the small slow layer.
- Population tests (t_half per A_start bin, dA/dt vs A) do not separate the laws (2 and 1 non-overlapping bins over
  all pairs). Real t_half falls outside m1's replicate range in 3 of 4 bins.

## Step 3b-1b, pooled fit on event-bootstrap simulated tails (`step3b1b_pooled_fit.py` + `waa_pooled.py`, 2026-10-10)

- Fit per DATASET, not per tail: shared law parameters, per-tail start W0 (profiled), per-tail speed factor
  exp(s (u_link + v_tail)) marginalised (Laplace over the tail part, quadrature over the link part), offset prior
  c ~ N(0, offset_std^2) (MAP), white + AR(1) noise per sub-link in the likelihood (ARMA(1,1) whitening). Candidates
  F1, F3, F4, F5 and the new F6 (state switch). BIC's N = distinct (link, event draw) pairs (median 285 vs 416 tails).
- Data: 10 event-bootstrap replicates (24 events with replacement; 289-832 tails), same draws for every truth.
- **Pre-set rule met in all three windows** (first 60 min, first 6 h, whole): own candidate wins 80-100 % for m1, m4,
  m5, m6; all key parameters within +-15 % (median), incl. m4's tau_s and phi and m6's A_c. Remaining confusion:
  m3 -> F5 in 20-40 %. Median dBIC to the runner-up 5-8.
- Caveat: the fit used the simulator's own noise model. Next: a misspecified-noise test (truth with two noise time
  scales, fit with one) before the real tails are fitted (step 4).
- Hint from the population statistics: the real median A / A_start at 30 and 60 min is below every calibrated truth's
  range, and at 180 min only m1's range holds it (real CIs wide).
- Decisions taken by Claude (delegated by the learners for this step): the step 3b-1 per-tail confusion figure and
  recovery table were demoted from key to diag (superseded; key limit 6); the m6 calibration table of step 3b-1 was
  promoted to product (read by this step).

## Step 3b-1c, stress test of the pooled fit (`step3b1c_stress_test.py`, 2026-10-11)

- Fake data only, first 60 min, truths m1, m4, m5, m6, candidates F1, F3, F4, F5, F6, 5 event-bootstrap replicates
  per truth and test (new draws; 135-753 tails). Rule fixed before the run: own candidate wins >= 4 of 5 and key
  parameters within +-30 % (median). 750 fits, 38 min on 3 processes; nothing was cut for time.
- **Verdict: PARTIAL.** 12 of 14 cells pass; both failing cells are truth m5.
- T1 (noise the fit does not know: white + AR(1) tau 20 min + AR(1) tau 300 min, same total slow variance, fast
  share 0.3; the fit keeps one AR(1)): m1, m4, m6 pass 5/5 (median dBIC -5.7, -29, -30). **m5 fails: F3 wins 5/5**, by
  1.8 BIC (a near tie; with the right noise F5 wins 4/5 by 1.9). T1 noise vs real dry statistics: within-tail std
  0.127 vs 0.141 dB, ACF 1 / 10 / 60 min 0.78 / 0.59 / 0.22 vs 0.80 / 0.57 / 0.22 (per sub-link a little worse than
  the 3a-fix noise).
- T2 (4 generic starts per candidate instead of the truth; best likelihood kept): all four truths pass (4/5, 4/5,
  4/5, 5/5). Start sensitivity (best vs worst start > 2 BIC): F1 5/20, F3 4/20, F4 10/20, F5 11/20, **F6 20/20**
  (median spread 59 BIC, max 267). Same winner as with the truth start in 19 of 20 datasets; in one m4 dataset F6
  wins from a generic start (by 3.8 BIC).
- T3 (truths recalibrated to the real median A / A_start at 10 / 30 / 60 min = 0.568 / 0.229 / 0.142; s = 1.25, and
  s = 1.5 for m1, m4): new m1 tau 25.7 min; m4 tau_f 14.3 min, phi 0.174; m5 a 0.0066, b 0.056; m6 k_run 0.057,
  A_c 0.17 dB. m1, m4, m6 pass (and m1, m4 with s = 1.5). **m5 fails on a parameter**: F5 wins 4/5 but a is x0.65.
- Decision by Claude, to be confirmed by the learners: in T3 the slow parameters (m4 tau_s = 193 min, m6 k_evap =
  0.0055 / min) were kept at their 3a-fix value, because with all parameters free they ran off to a constant floor
  (tau_s 7e11 min). Seen in the calibration printout, before any T3 fit.
- For step 4: fit every candidate from several starts (more for F6); read F3 and F5 as one family ("nonlinear in
  W"); F4 vs F6 by less than about 10 BIC is not decisive; do not interpret F5's linear rate a from the first hour.
- The win-share figure is diag, not key: the module already has 6 key outputs.

## Step 3b-2, weak-form SINDy on fake tails (`step3b2_sindy_fake.py` + `waa_sindy.py`, 2026-10-11)

- Fake data only: the T2 datasets of 3b-1c (seeds 3000+, 5 replicates; 20 of 20 bit-identical, their pooled fits
  reused) for m1, m4, m5, m6, plus m3, plus m1, m3, m5, m6 without spread. Own weak-form SINDy (pysindy not
  installed): library 1, A, A^2, A^3, max(A - A_c, 0) for A_c 0.1-1 dB; 21 min test windows; thresholded backward
  elimination (STLSQ failed a noise-free self-test on the near-duplicate switch terms); threshold by 5-fold CV over
  event draws (one-SE rule); stability selection (50 resamples, selected if >= 0.6). 190 SINDy runs, nothing cut.
- Offset: `handled` = a free constant per tail in dA/dt (exact for a linear law only); also `naive` and `oracle`.
- **Q1 FAIL (0 of 4 truths)**, first 60 min, handled, with spread: m3 and m5 return "A" in 5/5 (A^2 0/5); m6
  returns "A" with rate 0.05 / min = k_run + k_evap, no switch term; m1 keeps A (5/5) but gets a false A^2 in 3/5
  and a rate x1.79. Without spread only m1 is met (rate x0.94).
- Why: held-out error 0.36-0.58 for every threshold; the best model beats the one-term model by 0.00-0.06 (SE
  0.03-0.04), so the laws cannot be told apart. Post hoc check with the lowest-error threshold (added after the
  first run; decision table unchanged): unstable sets of 1-8 terms, rule met nowhere new.
- With the TRUE offset and no spread the first hour does give A^2 for m3 and m5 (5/5, coefficients x0.83, x0.77):
  the offset and the spread are what break it. Naive offset: a false constant term in nearly every result.
- **Q2 FAIL (misleading)**: m4 returns a clean single-variable law in 5/5 ("A" x3, "A + A^2" x2).
- **Q3 PARTIAL (mixture trap is real)**: m1 with spread gets A^2 in 3/5, 0/5 without spread. Time rescaling by the
  pooled-fit speed does not cure it (A^2 4/5, rate x0.48); with the true speed A^2 2/5, rate x0.28, so the
  rescaling as built is itself doubtful (fast tails enter with noise-only regressors times up to 20).
- **Q4 FAIL**: same family as the pooled-fit winner in 1 of 20 datasets (the pooled fit is right in 19 of 20).
- For step 5: do not run SINDy on the real tails as planned (one pooled law); the pooled fit is the working tool.
  A SINDy redesign (speed per tail inside the regression, offset prior, noise in the library terms, no switch
  terms below the noise level) must pass this fake-data test first. To be decided by the learners.
- The selection figure is diag (key limit 6 reached); if wanted as key, the proposed output to demote is
  `step3b1b_recovery.csv`.

### Plan change after 3b-2 (2026-10-11, decided by the learners)

- SINDy on raw tails failed the fake-data test (3b-2). Discovery is kept but split into two jumps: (1) the pooled
  fit with a FREE smooth drying curve instead of a menu of laws draws the drying curve from all tails; (2) SINDy
  names the equation of that clean curve.
- Old Step 5 (SINDy on raw real tails) is replaced by this: it becomes part of 3b-3 (fake) and Step 4 (real).
- New order: 3b-3 -> 4 -> 6 -> 3c -> 7 -> 8a-8e -> 9 (old step numbers kept). 3c moves to just before Step 7 (only
  needed for wetting); the tool for g(R, W) in Step 7 is decided after 3c, by the same fake-data testing.

## Step 3b-3, "draw the curve, then name it" on fake tails (`step3b3_curve_then_equation.py` + `waa_curve.py`, 2026-10-11)

- Jump 1: the pooled fit of 3b-1b / 3b-1c with a free curve, dW/dt = -f_i W k(W), log k = natural cubic spline in
  log W (5 knots: 0.05 dB, three quantiles of log A, top of the supported range 1.5-3.3 dB); 4 starts; penalty by
  2-fold CV over event draws. Jump 2: sparse regression of h(W) = W k(W) on W, W^2, W^3, switch terms. Data: the
  3b-1c / 3b-2 fake datasets (m1, m3, m4, m5, m6 with spread; extras m1 without spread, m1 and m3 with T1 noise),
  first 60 min and first 6 h. 65 fits, 2.1 h on 3 processes. Cut: the 6 h window of the extras.
- **Verdict by the pre-set rule: PARTIAL** (m1 FAIL, m3 PARTIAL, m5 PARTIAL, m6 PARTIAL, m4 PASS).
  m1: curve bends up at the lowest W (slope -0.24, true 0), coverage 68 % / 48 %, equation with false terms, W x2.1.
  m3, m5: equation right in the first 60 min (W^2 5/5, x0.94 and x1.07). m6: step test passes (ratio 7.7, 4/5) but
  the equation is W^2 + W^3, no switch term. m4: warning in 4/5 (free curve worse than F4 by 20-36 BIC).
- Broken in the judging: the curvature interval exists for 3 of 65 fits (every slope test fails for lack of an
  interval; the "60 min vs 6 h disagree" warning never fires), and the first knot (0.05 dB) is badly determined
  with the next knot at 0.30-0.40 dB (m1: k there 2-4 times too high in 9 of 10 fits).
- **Post hoc, above the second knot only (0.35-1.7 dB), first 60 min, with spread** (cached fits, verdicts
  unchanged): m1 slope -0.02 (true 0), level x1.03, equation "W" 5/5 x1.03: the mixing trap is gone there; m3
  slope 1.08 (true 1.0), W^2 x0.94; m5 slope 0.96 (true 0.92), W^2 x1.06; m6 slope 0.27 (true 0.26), level x1.02;
  m4 slope 0.37.
- The curve cannot show m6's step (A_c = 0.16 dB is in the undetermined stretch; above it the switch law is a
  smooth rising curve named "W + W^2") nor that m4 has two states (it draws a clean rising curve); only the menu
  fits F6 and F4 can.
- The free curve has a lower BIC than the TRUE menu law for m1, m3, m5 (median -21, -32, -62). Suspected cause:
  the likelihood has no quantization model, so the low end looks like faster drying. Not tested; it may also
  affect F6's A_c and F4's slow layer.
- For step 4: draw the curve in the first 60 min and read it only above the second knot (about 0.35 dB); keep the
  menu fits next to it (F6, F4; free curve worse than F4 by > 10 BIC = do not trust one curve). Repairs to test on
  the cached fake datasets first: a working interval (event bootstrap or profile), knots / range that do not leave
  0.05-0.35 dB to one value, a longer search. To be decided by the learners.
- The five k(W) figures are diag (key limit reached).

## Step 3b-3b, confirmation run of the free curve above 0.35 dB (`step3b3b_confirm.py`, 2026-10-11)

- Why: the "right above 0.35 dB" result of 3b-3 was found post hoc. Here it is tested on fresh fake data (event
  bootstrap seeds 4000 + replicate, 5 replicates, 3a-fix noise, s = 1.25) with range, tests and rule fixed in
  `waa.curve3b3b` before the run. Same fit as 3b-3, first 60 min; 5 knots from W_lo = 0.35 dB to W_hi (level with
  200 tail-minutes above; 1.9-3.5 dB per dataset), straight line in log-log below W_lo, nothing judged or named
  below W_lo. Band = min-max of the replicates; shape = slope of log k vs log W per replicate; jump 2 library W,
  W^2, W^3, switch at 0.5 and 0.75 dB. 35 fits, 32 min on 3 processes, nothing cut.
- **Verdict by the pre-set rule: PARTIAL** (m1 PASS, m3 PASS, m5 PASS, m4 FAIL; m6 not judged by design).
- m1, m3, m5: true k inside the band on 100 % of the judged range, level x0.99 / x0.97 / x1.03; slopes -0.01..0.05
  (true 0), 0.99..1.12 (true 1), 0.93..1.03 (true 0.92); equation "W" x0.99, "W^2" x0.97, "W^2" x1.13, one term in 5
  of 5, no false term. The mixing trap does not appear above 0.35 dB. T1 noise: same slopes and equations (only the
  m3 band coverage drops to 72 %).
- m6 (switch at 0.16 dB, below W_lo): curve right over the judged range (coverage 92 %, slopes 0.20-0.32, true
  0.23); jump 2 gives a three-term fit (W^2 + W^3 + switch at 0.5 dB), not a law. A user would conclude "mildly
  rising, neither exponential nor quadratic", which is true but does not say "switch".
- m4: BIC(free curve) - BIC(F4) > 10 in 2 of 5 (2.7, 4.2, 38.7, 53.2, 3.3): the warning sign is NOT confirmed. Its
  curve looks like m6's (slopes 0.36-0.47, three-term equation). Open: F4 had 3 generic starts here (4 in 3b-3);
  not rechecked because of the no-rerun rule.
- For step 4: run this setup on the real tails. The statement holds for W in [0.35 dB, W_hi] and the first 60 min:
  slope near 0 with "W" = exponential there; near 1 with "W^2" = quadratic; 0.2-0.5 with a multi-term equation =
  "mildly rising, law not identified by the curve" (switch or two layers: menu fits F6 / F4). Nothing is said
  about W < 0.35 dB. "Free curve vs F4" is reported as information, not used as a check. The band on real data
  needs an event bootstrap of the fit (not tested here as such).
- The five k(W) figures are diag (key limit reached).

## Step 3b-3c, the discovery method finished on fake tails (`step3b3c_method.py`, 2026-10-11)

- Five parts on the cached 3b-3b datasets (seeds 4000 + replicate, first 60 min), rules fixed in `waa.curve3b3c` before
  the run. 93 min on 3 processes + 15 min of PySR, nothing cut. pysindy 2.1.0 and PySR 2.7.0 (Julia 1.11.9) were
  installed with pip. **Overall PARTIAL: A PARTIAL, B PARTIAL, C PARTIAL, D FAIL; E information.** The PASS cells of
  m6 in B and C carry no information (see B).
- A, jump 2 alone on the TRUE h(W) + 2 % / 5 % noise, 20 draws: E-SINDy (bagged backward elimination of
  `waa_sindy.py`) exact in 20 of 20 in all 8 cells, m6 included; PySR exact 20 / 20 / 18 of 20 for m1 / m3 / m5,
  cannot write the switch; best-subset + BIC always holds the true terms but adds one in 3-8 of 20 draws for m1,
  m3, m5 (5 of 8 cells miss the 16-of-20 rule); pysindy's STLSQ fails at 5 % noise.
- B, free curve with an optional step (`waa_curve.StepCurve`, kept if BIC improves by > 6): kept in 5 of 5 for m1,
  m3, m5 AND m6 (gain 15-170), with W_c at 0.2-0.4 dB whatever the truth. The step tool is not usable.
- C, jump 3 (every candidate refitted in the pooled model, winner by held-out log-likelihood over event draws):
  the step law wins in 18 of the 20 single-curve datasets, by 16-65 over the true menu law; m4: F4 in 3 of 5. The
  pooled likelihood itself prefers a false feature near 0.2-0.4 dB. Suspected cause: the 0.3 dB quantization is
  not in the likelihood; NOT tested. Post hoc, menu laws only: BIC names the true family in 4, 5, 4, 5, 5 of 5.
- D, band from 20 event bootstraps of ONE dataset: 5-95 % band +-1 to 6 % wide, true k inside on 0-78 % of the
  range (rule 80 %) in 9 of 9 datasets: about 3 times too narrow; level of the fits x0.93-x1.23; slopes within
  0.01-0.11 of the truth.
- E: F4 with the F1-solution start changes BIC by at most 1.2; the m4 warning of 3b-3b stays at 2 of 5.
- Step 4 recipe (only what passed): free curve of 3b-3b without a step, read for W >= 0.35 dB, first 60 min;
  slope of log k vs log W (about +-0.1); jump 2 by E-SINDy with PySR as cross-check; no jump 3 by held-out
  likelihood; named laws by the menu fits with BIC (3b-1c rules); bootstrap band shown with the warning (level
  about +-20 %). Nothing about W < 0.35 dB or a switch until the low-W likelihood is repaired and B, C, D are
  repeated on the cached fake datasets. To be decided by the learners: repair first, or run Step 4 with the reduced
  recipe.
- The figures are diag; no key figure proposed.

## Step 3b-4, is the 0.3 dB quantization the cause of the false step? (`step3b4_quantization.py`, 2026-10-11)

- Diagnosis only (no new likelihood, no fix), fake data only, rule fixed in `waa.quant3b4` before the run. Truth m1,
  the skeleton, seeds (4000 + replicate), noise and spread of 3b-3b, 5 replicates, first 60 min; three simulator
  settings that differ only in the RSL quantization: none, 0.3 dB (the cached 3b-3b / 3b-3c fits), 0.1 dB. Fit: free
  curve as in 3b-3b, optional step as in 3b-3c part B. 17 min on 3 processes (projected 18), nothing cut.
- **Verdict: not the cause.** The step is kept in 5 of 5 in every setting; BIC gain 40-130 without quantization,
  81-143 at 0.3 dB, 49-136 at 0.1 dB (rule: > 6). The explanation suspected in 3b-3c is withdrawn.
- The fitted W_c is within 12 % of its start value (0.12, 0.2, 0.4 dB) in 15 of 15 fits: the position of the false
  step is set by the start list, not by the data.
- Real data (reading only): RSL in 0.3 dB steps on 719 of 728 sub-links, TSL in 1 dB steps (constant on 523; tails end
  before a TSL step), 10 s samples averaged to 1 min. The simulator matches this.
- Consequence: the Step 4 recipe of 3b-3c is unchanged (curve read for W >= 0.35 dB, no step tool, no jump 3 by
  held-out likelihood, nothing about a switch). Dropped: "add the quantization to the likelihood, then repeat B, C, D".
- Open: the cause of the false step. Untested suspects: the 1150 extra evaluations of the step fit (a no-step fit that
  was not fully converged), the offset prior / start floor / start profile at low W, the 0.02-wide step shape.
- The figure is diag; no key figure proposed.

## Step 4, the drying law on the REAL tails (`step4_real_drying.py`, 2026-10-12)

- First real-data step: frozen set v1, TRAIN only (429 tails, 24 events, 159 links; the 123 VAL tails were not read),
  first 60 min, recipe of 3b-3c with unchanged settings (`waa.real4`): menu fits F1, F3, F4, F5, F6 (4 starts each),
  free curve with 5 knots from 0.35 dB to W_hi = 2.83 dB, band from 20 event-bootstrap refits, jump 2 by E-SINDy /
  PySR / best-subset. 160 menu fits, 7 curve fits, 20 refits; 36 min on 3 processes (projected 1.2 h), nothing cut.
- **No drying law is identified.** The free curve FALLS with W: slope of log k vs log W = -0.87 on the main set
  (bootstrap refits -1.43..-0.82), a shape no fake truth had; 1 / k = 6 min at 0.35 dB, 33 min at 1 dB, 71 min at
  2.83 dB. A model-free look gives the same trend. Cause in the data: tails that stay flat or rise after t0_trace
  (21 % are still at >= 90 % of A_start after 10 min; 38 % in the events flagged event_rise_heavy, 11 % in the others).
- Sensitivity: without event_rise_heavy the slope is -0.05 ("exponential in this range", F1 tau 15 min); without
  level_shift / undershoot -0.79; alternative set -0.69; leaving out each big event -0.60, -1.60, -0.46.
- The curve is badly determined: the 4 starts end 800-1700 apart in cost (fake data: about 8), 5 of 20 refits were
  better from the flat start, band +-94 %. The dip near 0.9 dB moves with the third knot: not read.
- Jump 2: no formula. E-SINDy keeps all 5 library terms; the best one-term formula (0.032 W) misses the curve by 61 %
  rms; PySR's simplest is a constant h = 0.044 dB / min (49 %).
- **Fault found in the menu comparison**: `waa_pooled.py` profiles the start W0 of F1 and F4 without a bound, and keeps
  W0 >= 0.001 dB for F3, F5, F6 and the free curve. As run, F1 wins (tau 25.2 min, s 1.40; F4 +11.4 BIC, F6 +2610), but
  3 rising tails with a negative start are worth the whole gap to F3 at n = 1 (the same law). Like for like (row
  `F1+`): F6 ahead of the exponential by 21.5 BIC on the main set (over the sets: F1+ ahead by 0.9, to F6 ahead by
  216). Not decisive: the residuals are 4.1 times the dry-noise variance (2.1 on 6 h), so BIC margins are overstated,
  and F6 depends on its start (192 BIC). Stable in all 7 sets: F3 ends at n = 1 and F5 at b = 0, so a rate that rises
  with W (m3 / m5) is not supported in the first hour.
- Check 3 (first 6 h, information): F4 ahead of F1 by 5.5 BIC with tau_s = 16000 min, phi = 0.015 (a constant floor);
  F1 tau 18.7 min. Check 5: real A / A_start at 10 / 30 / 60 min = 0.568 / 0.229 / 0.142, inside the range of the
  recalibrated m1 only; the real slope is outside every fake truth.
- Literature: no drying time constant is quoted in the proposal, the OpenMRG summary or elsewhere in the repository;
  none was invented. To do before Step 8a: read the values from the proposal's refs [1] and [9].
- Reading: dW/dt = -h(W) with one curve and a random speed per tail does not describe the first hour of the TRAIN
  tails; a tail that stays flat or rises points to an input after the last rain trace (the g(R, W) side) or a wrong
  start time. The exponential with tau of 15-28 min is a placeholder, not a discovered law.
- To be decided by the learners before Step 6: (a) which tails define "drying" (all of v1; without the rise-heavy
  events; or tails re-started where the decay begins = a new frozen set); (b) repair the W0 bound and re-run the
  fake-data tests with some rising tails; (c) a larger or fitted noise level for drying tails.
- Claude's additions after the first run, from cached fits (no fit rerun, no setting changed): the like-for-like
  reading, the model-free rate table, the residual rows without the rise-heavy events; the E-SINDy coefficients are
  one least-squares fit on the kept terms (the bagged medians did not fit together).
- `step4_k_real.png` is the key candidate but saved as diag (6 key outputs); if promoted, demote
  `outputs/tables/step3b1b_recovery.csv`. Advice: wait.
