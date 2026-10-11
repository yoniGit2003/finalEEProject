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
