# wet_antenna: method, decisions and open issues

Three parts: (1) the rules that define a drying tail and the frozen tail set v1, (2) what each step decided, with the
numbers behind it, (3) the open issues. Settings: section `wet_antenna` of `research/config.yaml`. TRAIN events only
for every number unless stated. The step reports with all tables are in `outputs/reports/` (one per script).

# 1. The tail rules

### Tail definition (decided by the learners)

One tail per (event, sub-link) pair:

- **Pairs**: every row of `baselines.csv` with `missing == False` and `link_has_wet_slot == True`. These are train + val
  events on train links (test events and links have no baseline). All flag columns are kept, and `split_event` is
  recorded on every tail.
- **A(t)** = 1 min `tl_mean` (release) - `baseline_db` (`baselines.csv`). This is exactly `baseline.a_of_t` (checked on 3
  pairs: `outputs/tables/extract_tails_a_check.csv`).
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

### Step 2 rules (fixed before the run; settings `wet_antenna.selection`)

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

### Step 2b rules (after the step 2 gallery; settings `wet_antenna.refine_selection`)

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

### Tail set v1 (FROZEN 2026-10-08; the input of the WAA model)

Files: `outputs/products/tails_model.parquet` (1 min rows) + `tails_model_index.csv` (one row per tail); the
alternative set `tails_model_alt.*`. The set id `v1` is in config `wet_antenna.freeze.tail_set` and in the `tail_set` column.
File names carry no version: a later set gets a new id in the config and in the `tail_set` column. Every rule that defines v1, in order:

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

Settings: `wet_antenna.tails`, `wet_antenna.selection`, `wet_antenna.refine_selection`, `wet_antenna.freeze` in `research/config.yaml`. Do not change them
for v1.

# 2. What each step decided

### Step 0, inventory (`inventory.py`, 2026-10-08)

- A(t) is **not stored**: it is `tl_mean` (release) - `baseline_db` (`baselines.csv`), exactly `baseline.a_of_t`.
- `baselines.csv` covers train + val events on the 283 train links (566 sub-links); test has no baseline, so there are
  no test tails. 7375 significant train pairs have a wet bin and `use_in_waa_main`.
- link_state is per link, 5 min; the ±30 min margin is on the dry side only (unknown bins, not wet). So t0
  (`last_wet_end`) is normally an unknown margin bin and the first dry bin comes about 30 min later.
- QC flags that block a sample are already excluded from `tl_mean`; `tsl_step` and `flatline` are information only.
- Decision: build the tails on 1 min `tl_mean` with the frozen baseline v2; no new baseline in this module.

### Step 1, all drying tails (`extract_tails.py`, 2026-10-08)

- Tail = (event, sub-link) from t0 to the earliest of next wet bin / NaN run ≥ 15 min / t0 + 12 h. **No selection.**
- 20 355 tails (19 659 non-empty) on 566 sub-links, 69 events; train 15 370, val 4 985. A(t) equals `a_of_t`
  exactly on 3 random pairs.
- Return to the dry level (A < max(0.3 dB, 2 offset_std)): A_start > 2 dB tails 25 % back at 10 min, 59 % at 30 min,
  79 % at 3 h; A_start < 0.5 dB tails 94 % at 10 min. Strong tails dry slowly: WAA-like memory exists.
- 670 tails start on a wet bin (`last_wet_end` is searched only inside the network event); 566 are clipped at the
  season end. Both kept as rows, excluded later.
- Dry-weather rises ≥ 0.3 dB later in the tail: 56 % of tails, partly RSL quantization (upper bound).
- Decision: `tails_1min.parquet` (47 MB) is git-ignored and rebuilt by this script.

### Step 2, selection and fast-drop test (`select_tails.py`, 2026-10-08)

- Rules fixed before the run (`wet_antenna.selection`); no tail kept or dropped by its shape. VAL not looked at.
- Excluded 4023 of 15 370 (mostly `use_in_waa_main` False: 3421). Truncation: dry-weather rise (≥ 0.6 dB for ≥ 30
  min) cut 1249 tails (5.5 % of tail time), a TSL step 29.
- Main set (A_start ≥ 1 dB, ≥ 60 min usable): 959 tails, 373 sub-links, 33 events; control (< 0.5 dB) 8838.
- Fast drop (A_start - A at the first dry bin): median 1.11 dB = 83 % of A_start. It does **not** grow with L
  (Spearman ≥ 25 GHz -0.10 [-0.20, 0.00]; < 25 GHz 0.14 [-0.09, 0.31]), but the path-rain reference wet_A_max barely
  grows with L either (0.13), so the L test has little power. Sub-links 1 vs 2 agree (rho 0.87).
- Open after step 2: light rain after the first dry bin stays in 466 of 959 main tails -> fixed in step 2b.

### Step 2b, fixes after the step 2 gallery (`refine_selection.py`, 2026-10-08)

- Start moved to **t0_trace** = end of the last 5 min bin with any rain trace (radar R_max > 0 on the path, or a gauge
  within 1 km > 0) between t0 and the first dry bin.
- 381 of the 959 step 2 main tails fall below A_start = 1 dB at t0_trace: much of the step 2 "fast drop" happened
  while the radar or a gauge still showed rain. So the fast drop is not evidence for WAA.
- Two ends: lenient (step 2 ends, rise cut moved back to the last quiet minute, excess < 0.15 dB) and strict (also the
  first rain trace after t0_trace). Same 479 tails in both; 220 strict tails end at a rain trace.
- Noise exclusion: sub-links with dry-weather `offset_std` > 0.3 dB (one RSL step): 39 sub-links, 1077 tails.
- `level_shift` (|A| at 6-12 h > max(0.3, 2 offset_std)) and `undershoot` (A < -0.3 dB within 3 h) are flags, not
  exclusions; undershoot looks like a lasting level near -0.3 dB (baseline offset for that event?), not a dip.

### Step 2c, frozen tail set v1 (`freeze_tail_set.py`, 2026-10-08)

- Radar NaN on the path counts as a rain trace (unknown); rise cut v2 on all minutes from t0_trace (10 min median ≥ 0.6
  dB above its running minimum for ≥ 10 min).
- **Frozen set v1**: main 429 train / 123 val tails; alternative (A_start ≥ 0.6 dB) 1155 / 458.
- 3 events give 246 of the 429 train tails: resample or weight by event in the model step.
- Open: night-time rises without any rain trace (E20150701_0030): drizzle under the radar beam or condensation? RH is
  missing for those hours, so this is not decided.

### Step 3a, tail simulator (`calibrate_simulator.py` + library `simulator.py`, 2026-10-08)

- Skeleton = the 429 TRAIN v1 tails (same sub-links, links, events, usable lengths, A_start_trace); 1 min output.
- Laws calibrated by least squares on the real median of A / A_start over 0-6 h (noise-free truth, same survivors):
  m1 tau = 25.7 min; m2 median tau = 30.7 min, s = 1.5; m3 (n = 2) c = 0.064; m4 tau_f = 16.1, tau_s = 161 min,
  phi = 0.17; m5 a = 0.0036, b = 0.055. Median misfit (truth): m4 0.016 < m5 0.024 < m3 0.027 < m1 = m2 0.052.
- m2's spread s is not identifiable from the median (median of exp(-t/tau) = exp(-t/median tau)); s was taken from
  the real 10-90 % band (deviation from "median only", config `wet_antenna.simulator.m2_s_from_band`). s = 1.5 lumps all
  tail-to-tail variation; only m2 matches the band (0.05 vs 0.15-0.18 for the deterministic laws).
- Noise: Gaussian 10 s sigma per sub-link from the step 2b control tails (median 0.09 dB). It matches the share of
  minutes on a 0.3 dB level but NOT the within-tail std (0.05 vs 0.14 dB real; last-hour std 0.057 vs 0.083 dB).
  Real dry noise has a slow part the model lacks, so 3b with this model overstates how well laws can be separated.
- Decision: products `sim_tails_m1..m5.parquet`, `sim_tails_index.csv` and a noise sigma table (both replaced in
  3a-fix); switches for offset /
  noise / quantization / real NaN mask in `wet_antenna.simulator.switches` (all on here). Open: add a slow-noise term before 3b?

### Step 3a-fix, slow noise and a spread for every law (`calibrate_simulator.py` + `simulator.py`, 2026-10-08)

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
- Retired by the 3a-fix: the m2 simulated set, the old noise sigma table and the first step 3a report (no script
  writes them any more).

### Step 3b-1, the fitting method on simulated tails (`check_per_tail_fit.py` + `per_tail_fit.py`, 2026-10-09)

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

### Step 3b-1b, pooled fit on event-bootstrap simulated tails (`check_pooled_fit.py` + `pooled_fit.py`, 2026-10-10)

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
- The per-tail confusion figure and recovery table of step 3b-1 are superseded by this step; the m6 calibration
  table of step 3b-1 is read by this step.

### Step 3b-1c, stress test of the pooled fit (`stress_pooled_fit.py`, 2026-10-11)

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

### Step 3b-2, weak-form SINDy on fake tails (`check_sindy.py` + `sindy.py`, 2026-10-11)

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

#### Plan change after 3b-2 (2026-10-11, decided by the learners)

- SINDy on raw tails failed the fake-data test (3b-2). Discovery is kept but split into two jumps: (1) the pooled
  fit with a FREE smooth drying curve instead of a menu of laws draws the drying curve from all tails; (2) SINDy
  names the equation of that clean curve.
- Old Step 5 (SINDy on raw real tails) is replaced by this: it becomes part of 3b-3 (fake) and Step 4 (real).
- New order: 3b-3 -> 4 -> 6 -> 3c -> 7 -> 8a-8e -> 9 (old step numbers kept). 3c moves to just before Step 7 (only
  needed for wetting); the tool for g(R, W) in Step 7 is decided after 3c, by the same fake-data testing.

### Step 3b-3, "draw the curve, then name it" on fake tails (`check_free_curve.py` + `free_curve.py`, 2026-10-11)

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

### Step 3b-3b, confirmation run of the free curve above 0.35 dB (`confirm_free_curve.py`, 2026-10-11)

- Why: the "right above 0.35 dB" result of 3b-3 was found post hoc. Here it is tested on fresh fake data (event
  bootstrap seeds 4000 + replicate, 5 replicates, 3a-fix noise, s = 1.25) with range, tests and rule fixed in
  `wet_antenna.confirm_free_curve` before the run. Same fit as 3b-3, first 60 min; 5 knots from W_lo = 0.35 dB to W_hi (level with
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

### Step 3b-3c, the discovery method finished on fake tails (`check_discovery_method.py`, 2026-10-11)

- Five parts on the cached 3b-3b datasets (seeds 4000 + replicate, first 60 min), rules fixed in `wet_antenna.check_discovery_method` before
  the run. 93 min on 3 processes + 15 min of PySR, nothing cut. pysindy 2.1.0 and PySR 2.7.0 (Julia 1.11.9) were
  installed with pip. **Overall PARTIAL: A PARTIAL, B PARTIAL, C PARTIAL, D FAIL; E information.** The PASS cells of
  m6 in B and C carry no information (see B).
- A, jump 2 alone on the TRUE h(W) + 2 % / 5 % noise, 20 draws: E-SINDy (bagged backward elimination of
  `sindy.py`) exact in 20 of 20 in all 8 cells, m6 included; PySR exact 20 / 20 / 18 of 20 for m1 / m3 / m5,
  cannot write the switch; best-subset + BIC always holds the true terms but adds one in 3-8 of 20 draws for m1,
  m3, m5 (5 of 8 cells miss the 16-of-20 rule); pysindy's STLSQ fails at 5 % noise.
- B, free curve with an optional step (`free_curve.StepCurve`, kept if BIC improves by > 6): kept in 5 of 5 for m1,
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

### Step 3b-4, is the 0.3 dB quantization the cause of the false step? (`check_quantization.py`, 2026-10-11)

- Diagnosis only (no new likelihood, no fix), fake data only, rule fixed in `wet_antenna.check_quantization` before the run. Truth m1,
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

### Step 4, the drying law on the REAL tails (`fit_real_drying.py`, 2026-10-12)

- First real-data step: frozen set v1, TRAIN only (429 tails, 24 events, 159 links; the 123 VAL tails were not read),
  first 60 min, recipe of 3b-3c with unchanged settings (`wet_antenna.fit_real_drying`): menu fits F1, F3, F4, F5, F6 (4 starts each),
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
- **Fault found in the menu comparison**: `pooled_fit.py` profiles the start W0 of F1 and F4 without a bound, and keeps
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
- The main figure of this step is `outputs/figures/fit_real_drying_k_real.png`.

# 3. Open issues

- **Step 4 found no drying law on the real tails**: the free curve falls with W because about one tail in five
  stays flat or rises after t0_trace (mostly in the events flagged event_rise_heavy). To decide before Step 6:
  which tails define "drying"; and two repairs of the tools: `pooled_fit.py` lets F1 / F4 fit a rising tail with a
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
  (Spearman 0.13 for >= 25 GHz). See part 2, step 2.
- `wet_A_max` spans the link's whole wet phase in the network event (several showers for multi-run events).
