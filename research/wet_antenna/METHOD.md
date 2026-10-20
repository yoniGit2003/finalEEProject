# wet_antenna: method, decisions and open issues

**Final equation (Part 1 closed, S9, 2026-10-19):** A(t) - c = W_f + W_s + L c_band R(t - 1 min), with dW_f/dt = (0.72 u - W_f) /
tau_f (2.6 min wetting, 10.3 min drying), dW_s/dt = (0.28 u - W_s) / 240 min, u = a_link A_band (R / 5)^b, b = 0.54. The numbers,
files and documents: [FINAL.md](FINAL.md). How it was reached: "Part 2: S-steps" below, S9 last.

Parts: (1) the rules that define a drying tail and the frozen tail set v1, (2) what each step T0-T18 decided, with the
numbers behind it, then "Part 2: S-steps" (the reboot, S1 on), (3) the open issues. Settings: section `wet_antenna` of `research/config.yaml`. TRAIN events only
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

### T2 (was 2) rules (fixed before the run; settings `wet_antenna.selection`)

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

### T3 (was 2b) rules (after the step T2 (was 2) gallery; settings `wet_antenna.refine_selection`)

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

Step labels: T0-T18 since 2026-10-13 (the closed drying-tail exploration; plan: `docs/ROADMAP.md`). The headings
give the old number in brackets; the text under them, the reports and the script docstrings still use the old
numbers (there, "T1", "T2", "T3" are the three stress tests inside step T8, was 3b-1c, not step labels).
Old numbers: T0 = 0, T1 = 1, T2 = 2, T3 = 2b, T4 = 2c, T5 = 3a and 3a-fix, T6 = 3b-1, T7 = 3b-1b, T8 = 3b-1c, T9 = 3b-2, T10 = 3b-3, T11 = 3b-3b, T12 = 3b-3c, T13 = 3b-4, T14 = 4, T15 = 4b, T16 = 4c, T17 = 4c2, T18 = 7a.

### T0 (was 0), inventory (`inventory.py`, 2026-10-08)

- A(t) is **not stored**: it is `tl_mean` (release) - `baseline_db` (`baselines.csv`), exactly `baseline.a_of_t`.
- `baselines.csv` covers train + val events on the 283 train links (566 sub-links); test has no baseline, so there are
  no test tails. 7375 significant train pairs have a wet bin and `use_in_waa_main`.
- link_state is per link, 5 min; the ±30 min margin is on the dry side only (unknown bins, not wet). So t0
  (`last_wet_end`) is normally an unknown margin bin and the first dry bin comes about 30 min later.
- QC flags that block a sample are already excluded from `tl_mean`; `tsl_step` and `flatline` are information only.
- Decision: build the tails on 1 min `tl_mean` with the frozen baseline v2; no new baseline in this module.

### T1 (was 1), all drying tails (`extract_tails.py`, 2026-10-08)

- Tail = (event, sub-link) from t0 to the earliest of next wet bin / NaN run ≥ 15 min / t0 + 12 h. **No selection.**
- 20 355 tails (19 659 non-empty) on 566 sub-links, 69 events; train 15 370, val 4 985. A(t) equals `a_of_t`
  exactly on 3 random pairs.
- Return to the dry level (A < max(0.3 dB, 2 offset_std)): A_start > 2 dB tails 25 % back at 10 min, 59 % at 30 min,
  79 % at 3 h; A_start < 0.5 dB tails 94 % at 10 min. Strong tails dry slowly: WAA-like memory exists.
- 670 tails start on a wet bin (`last_wet_end` is searched only inside the network event); 566 are clipped at the
  season end. Both kept as rows, excluded later.
- Dry-weather rises ≥ 0.3 dB later in the tail: 56 % of tails, partly RSL quantization (upper bound).
- Decision: `tails_1min.parquet` (47 MB) is git-ignored and rebuilt by this script.

### T2 (was 2), selection and fast-drop test (`select_tails.py`, 2026-10-08)

- Rules fixed before the run (`wet_antenna.selection`); no tail kept or dropped by its shape. VAL not looked at.
- Excluded 4023 of 15 370 (mostly `use_in_waa_main` False: 3421). Truncation: dry-weather rise (≥ 0.6 dB for ≥ 30
  min) cut 1249 tails (5.5 % of tail time), a TSL step 29.
- Main set (A_start ≥ 1 dB, ≥ 60 min usable): 959 tails, 373 sub-links, 33 events; control (< 0.5 dB) 8838.
- Fast drop (A_start - A at the first dry bin): median 1.11 dB = 83 % of A_start. It does **not** grow with L
  (Spearman ≥ 25 GHz -0.10 [-0.20, 0.00]; < 25 GHz 0.14 [-0.09, 0.31]), but the path-rain reference wet_A_max barely
  grows with L either (0.13), so the L test has little power. Sub-links 1 vs 2 agree (rho 0.87).
- Open after step 2: light rain after the first dry bin stays in 466 of 959 main tails -> fixed in step 2b.

### T3 (was 2b), fixes after the step T2 (was 2) gallery (`refine_selection.py`, 2026-10-08)

- Start moved to **t0_trace** = end of the last 5 min bin with any rain trace (radar R_max > 0 on the path, or a gauge
  within 1 km > 0) between t0 and the first dry bin.
- 381 of the 959 step 2 main tails fall below A_start = 1 dB at t0_trace: much of the step 2 "fast drop" happened
  while the radar or a gauge still showed rain. So the fast drop is not evidence for WAA.
- Two ends: lenient (step 2 ends, rise cut moved back to the last quiet minute, excess < 0.15 dB) and strict (also the
  first rain trace after t0_trace). Same 479 tails in both; 220 strict tails end at a rain trace.
- Noise exclusion: sub-links with dry-weather `offset_std` > 0.3 dB (one RSL step): 39 sub-links, 1077 tails.
- `level_shift` (|A| at 6-12 h > max(0.3, 2 offset_std)) and `undershoot` (A < -0.3 dB within 3 h) are flags, not
  exclusions; undershoot looks like a lasting level near -0.3 dB (baseline offset for that event?), not a dip.

### T4 (was 2c), frozen tail set v1 (`freeze_tail_set.py`, 2026-10-08)

- Radar NaN on the path counts as a rain trace (unknown); rise cut v2 on all minutes from t0_trace (10 min median ≥ 0.6
  dB above its running minimum for ≥ 10 min).
- **Frozen set v1**: main 429 train / 123 val tails; alternative (A_start ≥ 0.6 dB) 1155 / 458.
- 3 events give 246 of the 429 train tails: resample or weight by event in the model step.
- Open: night-time rises without any rain trace (E20150701_0030): drizzle under the radar beam or condensation? RH is
  missing for those hours, so this is not decided.

### T5 (was 3a), tail simulator (`calibrate_simulator.py` + library `simulator.py`, 2026-10-08)

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

### T5 (was 3a-fix), slow noise and a spread for every law (`calibrate_simulator.py` + `simulator.py`, 2026-10-08)

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

### T6 (was 3b-1), the fitting method on simulated tails (`check_per_tail_fit.py` + `per_tail_fit.py`, 2026-10-09)

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

### T7 (was 3b-1b), pooled fit on event-bootstrap simulated tails (`check_pooled_fit.py` + `pooled_fit.py`, 2026-10-10)

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

### T8 (was 3b-1c), stress test of the pooled fit (`stress_pooled_fit.py`, 2026-10-11)

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

### T9 (was 3b-2), weak-form SINDy on fake tails (`check_sindy.py` + `sindy.py`, 2026-10-11)

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

#### Plan change after T9 (was 3b-2) (2026-10-11, decided by the learners)

- SINDy on raw tails failed the fake-data test (3b-2). Discovery is kept but split into two jumps: (1) the pooled
  fit with a FREE smooth drying curve instead of a menu of laws draws the drying curve from all tails; (2) SINDy
  names the equation of that clean curve.
- Old Step 5 (SINDy on raw real tails) is replaced by this: it becomes part of 3b-3 (fake) and Step 4 (real).
- New order: 3b-3 -> 4 -> 6 -> 3c -> 7 -> 8a-8e -> 9 (old step numbers kept). 3c moves to just before Step 7 (only
  needed for wetting); the tool for g(R, W) in Step 7 is decided after 3c, by the same fake-data testing.

### T10 (was 3b-3), "draw the curve, then name it" on fake tails (`check_free_curve.py` + `free_curve.py`, 2026-10-11)

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

### T11 (was 3b-3b), confirmation run of the free curve above 0.35 dB (`confirm_free_curve.py`, 2026-10-11)

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

### T12 (was 3b-3c), the discovery method finished on fake tails (`check_discovery_method.py`, 2026-10-11)

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

### T13 (was 3b-4), is the 0.3 dB quantization the cause of the false step? (`check_quantization.py`, 2026-10-11)

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

### T14 (was 4), the drying law on the REAL tails (`fit_real_drying.py`, 2026-10-12)

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

### T15 (was 4b), why do some real tails stall or rise? (`diagnose_stalled_tails.py`, 2026-10-12)

- Question: step 4's falling k(W) and the tails that do not dry at once. H1 contamination (something still wets the
  antenna after t0_trace) against H2 antenna physics (a law of A alone whose rate per dB falls with A). TRAIN tails of
  v1 only (429; VAL never read), first 60 min. All rules in `wet_antenna.diagnose_stalled_tails`, fixed before the run;
  200 menu fits and 4 curve fits, 31 min on 3 processes (projected 0.58 h), nothing cut.
- Classes (5 min median of A, A_start = mean of minutes 0-4; tried in this order): decays 240, plateau-then-decays
  64, rises 36, flat / slow 89. Not decaying: 44 % of all, 65 % in the rise-heavy events, 31 % in the others; per
  event from 0 % to 100 % (E20150804_1830 24 %, E20150826_0520 19 %, E20150617_0705 83 %, E20150701_0030 100 %).
  The per-tail table is `outputs/products/tails_model_classes.csv`.
- **H1 by its rule: not supported, but the test is weak.** The radar on the path is dry in all 429 tails by
  construction of v1 (tails start at the last trace and are cut at the next). Rain farther away (radar on link paths
  within 5 km, nearest City gauge, SMHI gauge) is LESS frequent in the non-decaying tails (78 % against 93 %, odds
  ratio 0.29, event bootstrap 0.10..0.58). Humidity is missing for E20150617_0705 and E20150701_0030 (69 tails, 60 of
  the 189 non-decaying, 22 of the 36 rising), so dew cannot be tested where it matters. "Any moisture" is true for
  98 % of both groups. Coherence: the other links of the event at the same time share the behaviour, ratio 1.85
  (1.27..2.43; rule 2). The two directions of a link agree in 81 % (64 % by chance). The rise does not grow with the
  path length (Spearman 0.00; 44 % not decaying for short and for long links). Rising tails: 50 % at night (decaying 8 %); odds ratio of night for not
  decaying 5.0 (1.6..14.5; information, added after the first look).
- **H2 by its rule: not supported (0 of 4).** B1 IQR of the transition level 0.35 dB (but that level is 0.8 x A_start by
  definition, so the test is empty); B2 the same link in another event: p = 0.35; B3 plateau duration against
  A_start +0.16 (-0.02..+0.40), against the rain amount -0.09; B5 slope of the model-free k(A) in the decaying tails
  of the clean events +0.28 (+0.05..+0.77); all tails -0.80, rise-heavy events -1.69. B4 (no rule): after its
  plateau a plateau tail follows the master curve of the decaying tails as well as they do (rms 0.117 against
  0.124 dB), and its plateau lies below that curve in 45 of 61: the tail waits, then dries like the others.
- Reading: the stall is a property of the EVENT (the weather at that time), not of the level, the link or the path
  length. The falling k(W) of step 4 comes from mixing waiting tails in; h(W) does not need a falling k. What the
  cause is stays open (wetting without rain on the references, or no drying at night).
- Caveats: "decays" is defined by a fast fall, which biases B4 / B5 against H2 and makes tau short; the decaying tails
  follow rain more closely (gauge wet in 41 % against 20 %), so part of their fast fall may be the last path rain.
- **Part C, two options in `pooled_fit.py`, off by default** (earlier steps unchanged; step 4 re-run from its cache
  gives identical tables): `bound_w0` (start of F1 / F4 inside [0.001, 30] dB like the other laws) and `noise_mult` /
  `fit_noise` (one fitted multiplier m on the noise variance of a dataset). Fake check on 6 cached datasets of 3b-3b
  against the cached fits of 3b-3c: **NOT passed by its rule (4 of 6)**: m1 replicate 1 changes its (wrong) winner F6 ->
  F5, tau of F1 within 1 %; m6 replicate 0 keeps F6 but k_evap moves from 0.0031 to 0.016 (two starts 0.2
  log-likelihood apart: not determined in the first hour); m = 1.03 to 1.08 in all six.
- Real sets with both options: F6 first in all four (F4 at +11.6 to +23, F1 at +17 to +37 BIC), with A_c = 0.09 to
  0.12 dB (F4: phi 3-4 %), i.e. a small floor far below 0.35 dB; F3 at n = 1, F5 at b = 0. F1 tau: main 19.3 min,
  without rise-heavy 14.9, "decays" 10.2, "decays" + plateau tails from the plateau's end 10.8. Free-curve slope
  -0.47, -0.51, +0.01, -0.15 (step 4: -0.87 and -0.05 for the first two; the curve is unstable, starts 55 to 355
  apart in cost). Noise multiplier 4.3, 3.5, 2.6, 2.7: drying tails are noisier than the dry-weather model in every set.
- Changed during the work (not the rules): the dew indicator is "unknown", not "no", where humidity is missing;
  added for information after the first look: the early / late split of k, the "night" row, the class-order choice.
- To be decided by the learners: which tails define h(W) (all; "decays"; "decays" + plateau ends; or a start of
  drying per tail); delay / weather-dependent rate / input without rain; make `bound_w0` the default; accept or
  repair the noise multiplier; humidity for the two events from another station; is the fast fall of the "decays"
  tails path rain.
- The main figure of this step is `outputs/figures/diagnose_stalled_tails_b5_rate_per_level.png`.

### T16 (was 4c), what is in a drying tail? (`characterise_tails.py`, 2026-10-12; Phase A of the revised plan)

- Why: steps 4 and 4b fitted laws to tails that wait, rise or stall. This step fits NO law and asks what the signal
  is (Q1-Q3 of `docs/ROADMAP.md`, plan change of 2026-10-12). TRAIN tails of v1 only (429 main, 1155 alternative; for
  the start rules the 10 270 TRAIN pairs that pass the v1 exclusions; VAL rows dropped on reading). All rules in
  `wet_antenna.characterise_tails`, written before the run; 3 min, nothing cut. Checks: A(t) equals the frozen set
  (3e-6 dB), the rise equals step 4b's, start rule S0 reproduces the frozen main set and the 4b classes exactly.
- **A1 (Q1), rise of a tail against the same sub-link at the same clock time on certified dry days** (link_state dry
  in the whole hour, >= 12 h after the last wet bin, not in a dry jump; median 46 control days per tail, every tail
  has a control). **Verdict by the rule: "rises are mostly water (not drift)"**: 24 of the 36 class-"rises" tails are
  above their dry-day 95th percentile (67 %; rule >= 60 %). Dry days move little in that hour: 95th percentile of the
  rise 0.30 dB by day, 0.36 dB at night; the level changes by -0.2..+0.3 dB over the hour (5-95 %), median 0.00. So no baseline v3
  for the daily cycle follows from this. Weaker than it looks: night-time rises alone 10 of 18 (56 %, "mixed" by the
  same rule), by day 14 of 18; "rises" is defined by rise > 0.3 dB, about the dry-day percentile (median rise 0.45
  against 0.37 dB); 21 % of ALL tails exceed (decays 17 %, plateau 34 %, flat / slow 3 %), 5 % expected by chance.
- **Dew: not decided.** Not "dew-like" by the rule (dew risk 76 % in the exceeding tails, 87 % in the decaying ones;
  the rule of "more than twice" cannot be met when the reference share is above 50 %). 22 of the 36 rising tails are
  in E20150617_0705 and E20150701_0030, where humidity is missing (69 tails). Daily-cycle amplitude of the sub-links:
  0.12-0.17 dB in every class. Context: 2653 positive dry jumps, 41 % start at night (night = 30 % of the hours),
  27 % in a dew-risk hour (14 % of the hours).
- **A2 (Q2), A(m) = alpha_event(m) + beta(m) L + gamma(m) [f >= 25 GHz] per minute, -30..+60 min, event bootstrap.**
  **Power check NOT passed** in both sets: beta is never above 0 at -30..-10 min (largest 0.08 dB/km), so by the rule
  A2 is not interpreted (N = 5 by the rule, not read). Why: the radar expects only 0.03-0.06 dB/km (0.08-0.13 dB for
  the mean tail) before t0_trace, against a mean A of 1.7 dB and an interval of about +-0.10 dB/km: the rain at the
  end of an event is light. Information: the interval of beta stays between about -0.09 and +0.07 dB/km at every minute after
  t0_trace (beta slightly negative around minute 0-5); alpha (>= 25 GHz, A at L = 0) falls 1.79 -> 0.98 dB in 10 min, 0.28 dB at 60 min.
- **A3 (Q3), three start rules with the v1 end rules and the main-set rule applied from the new start.**

  | | S0 = t0_trace | S1 = strict dry (no reference within 5 km, + 10 min) | S2 = signal onset (maximum within 60 min) |
  |---|---|---|---|
  | tails / events / links | 429 / 24 / 159 | 129 / 14 / 69 | 571 / 27 / 174 |
  | decays / plateau / rises / flat-slow | 56 / 15 / 8 / 21 % | 19 / 32 / 14 / 36 % | 61 / 11 / 3 / 26 % |
  | A / A_start at 10 / 30 / 60 min | 0.568 / 0.229 / 0.142 | 0.937 / 0.708 / 0.518 | 0.525 / 0.223 / 0.166 |
  | t_half | 14 min | 48 min | 12 min |
  | slope of the model-free k(A) | -0.80 (-1.55..-0.43) | -1.85 (-3.02..-1.05) | -0.56 (-0.92..-0.02) |

  S1 comes 25 min after t0_trace (median, quartiles 15-30) and 349 of the 429 v1 tails are below 1 dB by then: most
  of what v1 calls the tail happens while a reference within 5 km still shows rain; the S1 set is the slow tails (a
  selection on the outcome). S2 equals S0 for 70 % of the v1 tails (5 min or more later for 19 %) and adds 150 tails
  that reach 1 dB only after t0_trace. Stable: the falling k(A) (negative under every rule) and, between S0 and S2,
  the median shape and t_half.
- **Proposal for tail set v2 (not frozen; the learners decide):** keep S0 as the start and store S1 and S2 per tail
  (`outputs/products/tails_starts.csv`, one row per TRAIN pair), so that step 4d models the start of drying per tail
  instead of re-cutting the set. No law fit was used to choose.
- Added after the first run (information, no rule changed): the expected path attenuation in dB in the A2 table, the
  peak above the 5 min median AT t0_trace next to the "rise before onset" flag (26 % by the flag as defined, 14 %
  against the level at t0_trace), the night / day split of the class "rises".
- Changed outside this script: `classify` (the class rules of step 4b) moved to `tails.py` so that both steps use
  one copy; step 4b re-run, all its tables and its report identical. `tails.py` also gained `itu_p838` (ITU-R P.838-3
  k = a R^b, typed in and checked against the tabulated values at 10, 20, 30 GHz from memory, not against the document).
- To be decided by the learners: tail set v2 (the proposal, a re-cut, or an S1 set with a lower A_start rule);
  accept "not drift" (no baseline v3); humidity for the two June events from another station; accept that Q2 waits
  for Step 7 or design another test; whether the minutes with rain within 5 km (S0 to S1) belong to the drying law.
- The main figure of this step is `outputs/figures/characterise_tails_a2_length_per_minute.png`.

#### Decisions of the learners after T16 (was 4c) (2026-10-12)

1. Night / day rises are NOT baseline drift (A1: 67 % of the rising tails exceed their own dry-day percentile). No
   baseline v3. Dew or drizzle stays open (humidity missing for E20150617_0705 and E20150701_0030).
2. Tail set v2 = the v1 tails with the start S0 (t0_trace), plus S1 (strict dry: no rain within 5 km, + 10 min) and S2
   (signal onset) stored per tail (`outputs/products/tails_starts.csv`). Not re-cut. The start of drying is treated
   per tail.
3. New key finding: rain within 5 km stops about 25 min after t0_trace (median), and by then 81 % of the tails are
   below 1 dB. Most of what v1 calls the drying tail happens while it still rains nearby, so the fast early drop may
   be the antenna drying, light rain near or on the path ending, or both. Q2 (antenna or rain) is now the central
   open question; A2 (A against L) had no power on tails.
4. Order changed: 4c -> 4c2 -> 3c + 7 (WHOLE events on very short links with radar rain as input, so the model sees
   the rain end by itself) -> 4d (drying law on tails, informed by 7) -> 6 -> 8a-8e -> 9.

### T17 (was 4c2), does drying start when the rain nearby ends? (`check_drying_after_rain.py`, 2026-10-12)

- Question: once NO reference within 5 km shows rain (from S1 of step 4c on), how fast does A fall? Quick check,
  data only, no law fit. The 10 270 TRAIN pairs that pass the v1 exclusions (VAL never read), S0 and S1 from
  `tails_starts.csv`. Three starts per pair, each with the v1 end rules counted from it: S1 (also ends at the next
  rain within 5 km), a placebo start (S0 + the S1 - S0 delay of another pair, no "rain nearby" condition) and S0. A
  tail is kept if A at the start >= the threshold and it is >= 30 min long. t_half = Kaplan-Meier median of the time
  to half (a tail that ends first is censored). Rules in `wet_antenna.drying_after_rain`, written before the run;
  40 s, nothing cut.
- A1 / A2 (tails; t_half in min; A / A_start at 30 min):

  | threshold | S1 | placebo | S0 | S1 - placebo, t_half (90 %) | S1 - placebo, A / A_start at 30 min (90 %) |
  |---|---|---|---|---|---|
  | 0.25 dB | 1848; 41; 0.766 | 2211; 39; 0.734 | 3351; 20; 0.471 | +2 (-2..+5) | +0.031 (-0.006..+0.052) |
  | **0.4 dB** | 760; 40; 0.679 | 973; 37; 0.650 | 1915; 18; 0.390 | +3 (-2..+7) | +0.028 (-0.014..+0.064) |
  | 0.6 dB | 350; 40; 0.668 | 484; 34; 0.575 | 1176; 16; 0.340 | +6 (0..+20) | +0.093 (+0.046..+0.154) |

- **Verdict by the rule (0.4 dB): SAME speed.** S1 / S0 = 2.22 (the first half of SLOW is met), but S1 is within 8 %
  of the placebo and the interval of the difference holds 0. The slowness after S1 comes mostly from starting late,
  not from the rain within 5 km having ended.
- The contrast is diluted: 52 % of the placebo starts lie at or after the pair's own S1. Split (added after the
  first run, no rule), 0.4 dB: placebo starts while it still rains within 5 km t_half 33 min (465 tails, median delay
  15 min), at or after the own S1 43 min (508 tails, delay 20 min). A gap of 10 min in the direction of SLOW, mixed up
  with the delay. At 0.6 dB the S1 tails are slower than the placebo by the interval (A / A_start at 30 min).
- A3, k per level in the 1155 alternative-set tails, first 120 min: after S1 the rate per dB is 4 to 7 times lower
  than before it at the SAME level (k 0.004-0.012 against 0.023-0.069 per min); the placebo split gives the same
  ratios up to 0.7 dB (0.15, 0.17, 0.21) and higher ones above (0.31, 0.33, 0.29, 0.34 against 0.25, 0.23, 0.15,
  0.28 for S1). Caution: fast tails reach low levels early, so "before" holds the fast ones whatever the cause.
- Reading: a tail has two phases in time. About 15-20 min after t0_trace are fast; what is left (0.3-0.9 dB) halves
  in about 40 min at every threshold (tau about 58 min IF exponential), with no rain within 5 km on any reference.
  The end of the rain nearby is not a clear event in the signal, so the fast drop is mostly not "nearby rain
  ending". That does NOT make it the antenna: it starts at the last rain trace ON THE PATH. The reading written
  with the rule ("SAME = the fast drop is antenna drying") is stronger than the test. Q2 is narrowed, not closed.
- For Step 7: two time scales (about 10-15 min and about 40-60 min) must be possible; whether the fast one scales
  with L on very short links is the test.
- Observation (gallery): in E20150617_0705 (169 of the 760 S1 tails) several links rise by about 1 dB for 10-15 min
  at the same clock time while the radar is valid and shows nothing within 5 km.
- Changed outside this script: `reference_traces` (the rain traces per link) and `tail_from_start` (the v1 end
  rules from any start) moved from `characterise_tails.py` to `tails.py`; step 4c re-run, all its tables identical.
- The main figure of this step is `outputs/figures/check_drying_after_rain_median_curves.png`.

### T18 (was 7a), very short links: enough data, and the fast first phase? (`short_links.py`, 2026-10-13)

- Question (Q2): is the fast first phase of a tail (4c2: about 15-20 min) the antenna, or rain on the path that the
  radar no longer sees? Path rain scales with L, water on the antennas does not. TRAIN links and events only, no law
  fit; all rules in `wet_antenna.short_links`, written before the run (only counts per length class had been seen).
  21 s, nothing cut.
- **Part 1, inventory** (length classes < 0.6, 0.6-1, 1-2, 2-4, >= 4 km): TRAIN links 25 / 27 / 75 / 94 / 64; pairs
  that pass the v1 exclusions 685 / 978 / 2553 / 3929 / 2125 on 18 / 23 / 58 / 79 / 45 links; v1 main tails 30 / 35 /
  114 / 140 / 110; tails with A_start >= 0.4 dB and >= 60 min 113 / 163 / 451 / 689 / 483. Below 1 km there are only
  28 and 38 GHz sub-links (below 0.6 km: 12 and 38); all links below 25 GHz are >= 4 km. Rain on a link: 26 TRAIN
  events per link (median), radar R_mean in wet bins 0.5 mm/h median, 3.5 at the 90th percentile.
  **"Enough for Step 7" by the rule (>= 10 links, 20 events, 200 pairs): yes below 0.6 km** (18 links, 37 events, 685
  pairs); below 1 km 41 links, 41 events, 1663 pairs.
- **Part 2, the two phases by length.** 1899 tails from S0 (A_start >= 0.4 dB, >= 60 min, followed to 3 h). Fast
  drop = A_start - A at 15-20 min: 0.27 / 0.30 / 0.30 / 0.30 / 0.28 dB by length class, relative 0.44 / 0.40 / 0.48 /
  0.47 / 0.43; class shares of step 4b the same in every length class (decays 52-56 %).
  SHORT (< 1 km) against LONG (>= 2 km), matched on event and A_start bin (274 against 986 tails): relative fast drop
  0.407 against 0.435, **ratio 0.94 (0.84..1.15)**, difference -0.028 (-0.076..+0.059); below 0.6 km 1.01
  (0.78..1.16); at 28 GHz 0.84 (0.69..1.09), at 38 GHz 1.01 (0.77..1.41). Path rain would give about 0.2 (the ratio
  of the path lengths) to 0.26 (radar-expected attenuation).
- **Verdict by the rule: UNDECIDED.** The ratio passes the "antenna" limit (>= 0.7), but the power condition fails:
  the interval of the difference (0.136 wide) had to be below half the difference path rain alone would give, and
  that is only -0.086 (limit 0.043), because the radar expects 0.02 dB (SHORT) and 0.08 dB (LONG) of path rain in
  the 30 min before a tail: 3 % and 11 % of A_start. The rule was not changed.
- Reading (Claude's, beyond the rule): the fast drop does not scale with L over a factor of 10 in length, and it is 4
  (LONG) to 14 (SHORT) times the radar-expected path rain. This speaks for the antenna, or for anything else that
  does not scale with L (a rain cell smaller than every link; a level error shared inside an event).
- Slow phase (information): A at 50-60 min / A at 15-20 min is 0.53 (SHORT) against 0.63 (LONG), -0.10
  (-0.18..+0.05); at 38 GHz 0.40 against 0.58 (-0.17, -0.28..-0.02). A / A_start at 60 min 0.17 / 0.22 / 0.20 / 0.26 /
  0.34 by length. At 180 min 0.10-0.16, i.e. about 0.08 dB = the baseline uncertainty: after 1-2 h a tail shows the
  baseline more than water. From 50-60 to 170-180 min A only halves once more, so the slow phase is not one
  exponential (or it ends on a floor).
- **Part 3, whole events on the 41 links below 1 km** (1663 pairs; descriptive). The rain on a link in an event is
  several showers over 9 h (median), so the timing is taken on the last shower (changed after the first run; no rule
  involved). There A peaks 32 min before link_state ends the rain, is at half the peak 11 min after the peak, at
  19 % of the peak when the rain ends and at 4 % 30 min later (28 and 38 GHz alike). Peak A against peak radar
  rain: Spearman 0.85; median peak A 0.8 / 1.2 / 1.9 / 3.3 / 6.9 dB for 0.5-1 / 1-2 / 2-5 / 5-10 / > 10 mm/h, no
  saturation level; the expected path part at > 10 mm/h is 2.3 dB even on these links. A rising before the first
  radar trace: 4.2 % of the tested pairs (65 of 1531) against 1.4 % for a placebo an hour earlier.
- For Step 7: links below 0.6 km are enough by the rule, below 1 km as the larger set with L as an input; whole
  multi-shower events; two time scales, the fast one overlapping the fading rain; a path term also on short links
  (or light and moderate rain only); no wetting before rain as a general feature; frequency band as an input.
- To be decided by the learners: accept "the fast phase is not path rain" as a working conclusion although the
  pre-set rule says undecided (or define a new test first); the length cut of Step 7.
- Part 0: the paper Tiede, Chwala, Siart & Eibert 2025 (Adv. Radio Sci. 23, 13-19) is recorded in `docs/ROADMAP.md`;
  checked on the publisher's page only (the three frequencies and the weather finding are not confirmed from it).
- The main figure of this step is `outputs/figures/short_links_fast_drop_by_length.png`.

# Part 2: S-steps (the reboot of 2026-10-13)

System identification on short links: the rain R(t) is the input, the attenuation A(t) the output, on whole events
(no tail cutting). Plan: `docs/ROADMAP.md`, "Reboot (2026-10-13)". RULE: WAA is never estimated by subtracting a
power-law rain term; the path part is told apart from the antenna part only through the path length L. ITU-R P.838 is
used only to build fake data and as a size reference.

### S1, the model ladder M0-M4 on the links below 1 km (`system_id_ladder.py` + library `system_id.py`, 2026-10-13)

- **Data (TRAIN only).** 41 links, 67 sub-links, 42 events, 1718 link-events (one event on one sub-link; 803 at
  28 GHz, 915 at 38 GHz), 367 190 points of 5 min (93 081 wet, 274 109 dry); L 0.14 to 0.97 km (median 0.67); 18
  links have a City gauge within 1 km. Exclusions as in v1, decided without any fit (baseline present, sub-link
  good, use_in_waa_main, offset_std <= 0.3 dB). Window: 2 h before the first wet bin of the link to 3 h after its
  last one (median window 14.8 h). No window had a TSL step. Radar: median 0.48 mm/h on wet points, 1228 points
  above 10 mm/h.
- **Offset c, not fitted, the same for every model:** the median of A over the dry minutes from 120 to 30 min before
  the rain; 124 link-events had fewer than 30 such minutes and got c = 0 (flagged, kept).
- **Models.** A_pred(t) = W(t) + L P_band(R(t - d)) + c, simulated at 1 min, compared on 5 min means. P_band and
  F_band are free curves (piecewise linear on 0, 0.5, 1, 2, 5, 10, 20, 40 mm/h, 0 at R = 0, never falling); a_link
  is one amplitude per link with a ridge on log a_link (prior sd 0.5); d is one radar lag for all links (-10..10 min).
  M0: W = 0. M1: W = a F(R). M2: dW/dt = (a F - W) / tau. M3: tau_wet / tau_dry. M4: a fast store (share phi, tau_f_wet
  / tau_f_dry) and a slow store (tau_s). Loss: Huber (0.3 dB) on the 5 min residuals; 4 generic starts per model.
- **Judge and stop rule (fixed before the run).** Leave-one-event-out RMSE of A; 1000 event-bootstrap draws, paired
  between models. Chosen = the simplest model whose paired difference to the best has a 90 % interval that holds 0.
- **Fake check first: PASS.** Truths M3 (tau_wet 3, tau_dry 15) and M4 (phi 0.7, 3 / 10 / 60 min) driven by the real
  radar rain, windows and lengths, noise variance x 3, 0.3 dB quantization, unknown offsets, 3 replicates each. The
  true rung was chosen in 3 of 3 for both truths. Recovered / true: tau_wet x1.00, tau_dry x1.04; phi x0.98, tau_f_wet
  x0.99, tau_f_dry x1.03, tau_s x1.09; d = 3 (true 3) every time. **Limit of the check:** the fake A is driven by the
  same radar series that the fit gets (no input error), and its held-out RMSE (0.275 dB) is below the real one.
- **The ladder on the real data** (held-out RMSE, mean over events, 90 % interval; difference to the best):

  | model | RMSE (dB) | minus best | skill over M0 | RMSE in rain | integral after rain, predicted / observed |
  |---|---|---|---|---|---|
  | M0 | 0.417 (0.338..0.506) | +0.027 (0.019..0.036) | 0 | 0.626 | 0.28 |
  | M1 | 0.403 (0.328..0.489) | +0.013 (0.007..0.020) | 3.3 % | 0.601 | 0.32 |
  | M2 | 0.394 (0.317..0.482) | +0.003 (0.001..0.006) | 5.6 % | 0.568 | 1.19 |
  | M3 | 0.395 (0.318..0.483) | +0.004 (0.002..0.007) | 5.4 % | 0.569 | 1.04 |
  | M4 | 0.390 (0.315..0.477) | best | 6.4 % (4.7..8.2) | 0.563 | 2.39 |

  **Best M4; chosen by the stop rule: M4.** It is the top rung, so the ladder did not stop. All rungs leave 0.24 dB
  1-3 h after the rain (noise and baseline). The mean peak error is -0.67 dB (M1) to -0.87 dB (M3, M4).
- **M4 on all events** (90 % interval from 100 event-bootstrap refits): d = -1 min (-2..3); phi 0.71 (0.53..0.79);
  tau_f_wet 17.3 min (10.4..65.1); tau_f_dry 8.3 min (5.3..9.9); tau_s 513 min (210..600, the bound); sd of log a_link
  0.19. F (dB, a_link = 1) at 2 / 10 mm/h: 2.9 / 3.2 at 28 GHz, 3.3 / 3.4 at 38 GHz (flat between 2 and 10 mm/h; the
  rise above 10 mm/h rests on 0.33 % of the points). P at 10 mm/h: 1.5 (1.2..3.4) and 2.5 (2.0..5.0) dB/km (ITU-R
  P.838 size: 1.8 and 2.8).
- **Sensitivity (M4; information).** (a) links below 0.6 km: RMSE 0.386, tau_f_dry 9.6, tau_f_wet 42, tau_s 594,
  P(10) 4.7 / 6.3. (b) only points with R <= 2 mm/h, no path term: 0.361, tau_f_dry 9.3, tau_f_wet 21, tau_s at the
  bound. (c) City gauge (1 min) as input, 749 link-events: **0.339 against 0.402 for the radar run on the same
  points; tau_f_wet 2.4 min**, tau_f_dry 10.6, tau_s 488. (d) the loss on 1 min A: 0.390, tau_f_dry 7.4, tau_f_wet
  16. (e) time constants per band: 0.391, the same values in both bands (tau_f_dry 8.3 / 8.3).
- **Checks against earlier facts (held-out predictions).** After the last shower, median A / peak: observed 0.21,
  0.08, 0.03, 0.02, 0.01 at 0, 15, 30, 60, 120 min; M4 0.15, 0.07, 0.05, 0.04, 0.035 (the fast phase is there, the
  slow store is too slow); M1 0.11 then 0. Fast drop (0-5 min minus 15-20 min): observed 0.33 dB below 0.6 km and
  0.36 dB at 0.6-1 km; M4 0.14 and 0.14 (no growth with L, as T18, but less than half the size). Per-shower peaks:
  observed 0.82 dB (first) and 1.00 dB (later showers) for 0.86 and 0.81 mm/h; M4 0.58 and 0.62; M1 0.97 and 0.92.
- **What it means.**
  - WAA on these links is dynamic, not static, but the whole ladder removes only 6 % of the error of "no WAA".
  - The stable number is the fast drying time: about 8 min (7-11 min in every run), in line with the fast phase of
    T17 / T18.
  - The wetting time is NOT determined: 17 min with the radar, 2.4 min with the gauge. The slow wetting is most
    likely the radar (5 min, 2 km cells, aloft), not the antenna. This is a reading of run (c), not a test.
  - The slow store is NOT identified (tau_s at the bound; 2.4 x too much A after the rain). It may be a level drift
    that the offset rule does not see.
  - The split between the path part and the antenna part is not firm with L below 1 km: P(10) is 1.1-4.6 (28 GHz)
    and 2.2-6.1 (38 GHz) dB/km over the rungs.
- **Proposals for S2 (the learners decide):** check or run S2 with the gauge input on the links near a gauge; add a
  metric on the peaks; leave the slow part free and judge it on 60-180 min and on the integral after the rain; do
  not read F above 10 mm/h; P needs the long links (S3).
- **Unclear:** all leave-one-event-out refits stop at the iteration cap (warm starts, the same for every rung); the
  spread of the loss over the starts is 70 for M4 and 31 for M2 (one start each ends higher), below 1 for the
  others; d trades against the wetting time constant (M0-M3: d = 2-3, M4: d = -1); the Huber loss lowers the peaks.
- Report: `outputs/reports/system_id_ladder.md`; product for S3 / S4: `outputs/products/system_id_link_events.csv`.
  The main figures are `outputs/figures/system_id_ladder_rmse.png` and `system_id_ladder_after_last_shower.png`.
- Note added after S2: S1 uses the same parametrisation of P and F, so the optimiser trap found in S2 (P near 0 at
  R <= 2 mm/h) may be in the S1 numbers of P and F as well.

### S2, the WAA law on links up to 3 km, the slow part and the discovery rung M5 (`system_id_s2.py` + library `system_id.py`, 2026-10-14)

- **Question.** The five things S1 left open: P and F not separated by L below 1 km; the wetting speed depends on the
  rain input; the slow store too high after the rain (water or drift?); the peaks too low; no model learned freely.
  Every rule is in config `wet_antenna.system_id_s2`, written before the run. TRAIN only; VAL never read.
- **Data.** The TRAIN links below 3 km with the S1 rules: S<3 = 144 links, 257 sub-links, 51 events, 6653 link-events
  (1.35 M points of 5 min, 355 k wet); nested S<2 (4381) and S<1 (1718, the S1 set). Bands 28 GHz 75 links, 32 GHz 16,
  38 GHz 90 (none merged). Gauge sets: the 1998 link-events with a City gauge within 1 km, once with the 1 min gauge
  rain (GAUGE) and once with the radar on the same link-events and points.
- **Models.** M1 static; M3 one store (tau_wet / tau_dry); M3D = M3 + a linear level drift per link-event (slopes
  N(0, s_b^2), s_b by EM); M4 fast + slow store; M5 dW/dt = s_band a_link g(R, W) - h(W) with g free on 7 x 6 knots
  (g(0, W) = 0, rising with R) and h free and rising, smoothness lambda from {1, 100, 10 000} by the held-out RMSE on
  S<3 (picked: 1, the least smoothing). Fits cheaper than S1 (iteration caps, lag grid -4..6). Judge and stop rule as
  S1 (order M1 < M3 < M3D < M4 < M5), plus a peak metric per shower and the integral after the rain.
- **Fake checks (first, 3 replicates each).**
  - Ta (input explanation) **PASS**: with the exact 1 min rain M3 recovers tau_wet 3.0 / tau_dry 8.0; with a 5 min
    mean x log-normal(0.3) input the fitted tau_wet is 6.3 = 2.1 x the truth. The slow radar wetting of S1 is the input.
  - Tb (drift or slow store) **PASS**: truth M3 + drift read DRIFT-LIKE 3/3, truth M4 (tau_s 300) WATER-LIKE 3/3.
  - Tc (F and P on S<1 against S<3) **FAIL**: F within 16 %, but P(2 mm/h) comes out about 0 (x0.05 on S<3, x0.002 on
    S<1) and P(5), P(10) 13-26 % low, with the true rain as input. A check after the run (not in the config): a refit
    from the true P reaches a LOWER loss (13 741 against 13 792 on replicate 0) and recovers P(2) = 0.42 (truth 0.38).
    **An optimiser trap:** P and F are sums of exp(parameter) increments; an increment near exp(-8) has no gradient and
    stays. So the real P at R <= 2 mm/h and the F / P split are shaped by the optimiser.
  - Td (M5 on M3 data) **FAIL**: M5 predicts as well as M3 (0.280 against 0.279 dB) but its levels are 0.36-0.52 x
    the truth and t_wet x1.56: the M5 surfaces are not identified.
- **The ladder** (held-out RMSE, dB; the 90 % interval of the difference to the best):

  | set | M1 | M3 | M3D | M4 | M5 | best / chosen |
  |---|---|---|---|---|---|---|
  | S<1, radar | 0.403 | 0.394 | 0.392 | 0.390 | 0.390 | M4 / M3D |
  | S<2, radar | 0.440 | 0.429 | 0.428 | 0.425 | 0.425 | M5 / M3D |
  | **S<3, radar (decides)** | 0.458 | 0.448 | 0.448 | 0.444 | 0.442 | **M5 / M4** (M4 - M5 +0.0024, -0.0002..+0.0052) |
  | gauge links, GAUGE | 0.408 | 0.354 | 0.353 | 0.350 | 0.357 | M4 / M3D |
  | gauge links, radar | 0.428 | 0.419 | 0.416 | 0.414 | 0.412 | M5 / M3D |

  **Chosen for S2: M4.** Memory is needed (M1 -> M3 is 0.01-0.05 dB); every rung above M3 adds a few thousandths of a
  dB. Peaks too low in every rung (-0.4 to -1.0 dB, median per shower); the smallest peak error is M1's.
- **Length rule: DISAGREE** (read on M4): 1 of 8 key numbers of S<3 lies outside the S<1 interval, tau_f_dry 5.5 min
  (3.1..7.1) against 8.3 (5.7..9.8). All M3 numbers agree. So S<1 stays the main set. Not validated (Tc failed).
  Held-out residuals of M4 during rain, per length bin: -0.24 to -0.29 dB, no trend with L.
- **Gauge rule.** GAUGE input: M4 tau_f_wet 2.7 min (2.4..3.1), tau_f_dry 10.6 (9.5..12.2), peaks -0.67 dB; M3
  tau_wet 3.0, tau_dry 12.9. Radar on the same points: M4 17.4 / 7.1 min, peaks -0.89 dB. **Drying does NOT agree**
  (10.6 against 7.1 min, more than 30 %). With the gauge F keeps rising to 5 mm/h (28 GHz: F(0.5) 0.17, F(5) 3.8 dB);
  with the radar it is flat at about 3 dB from 2 mm/h.
- **The slow part: SLOW RISE, FORM UNDECIDED.** With a free drift allowed per held-out link-event M4 is not better
  than M3D (-0.0008, -0.0036..+0.0017); the M3D slopes average +0.006 dB/h (+0.002..+0.011), Spearman with the rain
  depth 0.15. After the last shower the observed median (A - c) / peak falls 0.044 -> 0.020 -> 0.010 at 30 / 90 / 175 min
  (e-folding 76, then 122 min); M3 and M3D reach 0 within an hour (integral 0.55 x, 0.27 x), M4 (tau_s 395 min on
  S<3) 1.24 x, M5 1.40 x. **A small slow part (2-4 % of the peak) with a time constant of about 1-2 h; store or ramp
  not decided.**
- **M5 and its name.** M5 (lambda 1) is the best on S<3 by a hair, but its surface is odd (g(R, 0) about 0.005 dB/min,
  t_wet 325 min), its starts disagree (spread 77 948 loss units) and Td failed. It qualified for naming: E-SINDy
  (200 bags, weights = visits of the 1 min states) found no good g (relative error >= 0.45 for every lambda); picked
  g = 0.0026 R + 0.0001 R^2 - 0.52 R/(1+R/0.5) + 0.35 R/(1+R/1) + 0.074 W, h = -0.53 W - 0.053 W^2 + 0.95 max(W - 0.2, 0).
  PySR (cross-check): g = 0.231 - exp(-R W (W - 0.384) - 1.475), h = exp(-0.198 - 1.094 / (W - 0.292)). **M5n (the
  E-SINDy formula refitted) could not start**: h stops rising above 3 dB while g grows with W, so W runs away on the
  links with a large a_link s_band (566 dB after 2 h of 2 mm/h); the refit stopped after 1 iteration; **no formula was
  named.** A naming that gives a usable model needs conditions (h rising, g not growing faster than h).
- **The current WAA equation (M4, S<1, radar; 90 % intervals).** A = W_f + W_s + L P_band(R(t + 1 min)) + c;
  u = a_link F_band(R); W_f -> phi u with tau_f_wet / tau_f_dry, W_s -> (1 - phi) u with tau_s. phi 0.71 (0.59..0.79);
  tau_f_wet 17.2 (8.9..35.7) with the radar, **2.7 (2.4..3.1) with the gauge (the one to use)**; tau_f_dry 8.3
  (5.7..9.8) radar, 10.6 (9.5..12.2) gauge; tau_s 512 (225..600, not identified). F(0.5 / 2 / 5 mm/h): 28 GHz 1.46 /
  2.82 / 3.18 dB, 38 GHz 1.35 / 3.18 / 3.22 dB; a_link 10-90 % 0.77-1.23. P(10): 1.5 (28 GHz), 2.5 (38 GHz) dB/km; P
  at R <= 2 mm/h not trusted.
- **Unclear.** The optimiser trap (above); 21 of 26 model-set fits reach the fold iteration cap (40) in every fold
  and the spread over starts is large for M3D, M4, M5 (the small differences between rungs may move with longer
  fits); tau_s not identified; the radar wetting time of M3 is only bounded from below (43 min, 15..174); the gauge is
  up to 1 km away (d = -3 min); the computer slept 11.6 h during the run (three fit times in the fits table include it).
- **Proposals for S3 (the learners decide):** fix the parametrisation of P and F (no vanishing gradient, or starts
  with P of the ITU size) and repeat Tc and the length rule; test a slow part with a time constant of 1-2 h (a bounded
  slow store, or a store plus a ramp), judged on the after-rain integral as well; take the wetting from the gauge input
  and find out why gauge and radar drying differ; look at the peaks; refit M5 longer and from more starts, and name it
  only with conditions that keep W bounded.
- Report: `outputs/reports/system_id_s2.md`; product: `outputs/products/system_id_s2_link_events.csv`; main figures
  `outputs/figures/system_id_s2_ladder.png`, `system_id_s2_after_last_shower.png`, `system_id_s2_m5_surfaces.png`.

### S2b, repair of the S2 model: the P / F trap, a bounded slow store, longer fits (`system_id_s2b.py` + library `system_id.py`, 2026-10-14)

- **Question.** S2 left an optimiser trap in P and F, a slow store that overshoots after the rain, and fits stopped at
  their caps. S2b repairs these three things only (no M5, no new data, no VAL). Every rule is in config
  `wet_antenna.system_id_s2b`, written before the run; data, judge and fake noise as S2. Recorded first in the ROADMAP:
  "WAA v1" (the fallback), the remaining plan and the parked open questions.
- **Changes.** (1) The increments of P and F between the knots are plain parameters >= 0 (L-BFGS-B bounds) instead of
  exp(parameter); 4 P starts = 0, 0.5, 1, 2 x the ITU-R P.838 size (starting values only). (2) M4b = M4 with tau_s in
  30-240 min; ladder M1 < M3 < M3D < M4b, the free M4 as a reference only. (3) 8 starts per model, caps 1000 per start
  and polish (5 x S2), 120 per fold (3 x S2); bootstrap refits 30 (as S2). Chosen model = the one chosen on S<1 (main).
- **Time and cuts.** One M4b start on S<3 ran 1000 iterations without converging (4.4 min), so the projection was
  7.6 h against 5 h and all three cuts of the rule were applied: S<2, the gauge-links radar comparison, fake replicates
  3 -> 2. 12.4 h of fitting on 3 processes (4.2 h wall clock; the first attempt was stopped by the session for low memory
  and restarted; every finished start is now saved).
- **Fake checks.** **Tc again: PASS** (the same data as S2): on S<3 all 30 F and P knots (0.5-10 mm/h, three bands)
  within 15 % of the truth, P(2 mm/h) / truth 1.01-1.04 (S2: 0.05). **Te (truth M4b, tau_s 90 min): PASS**: M4b chosen
  2 of 2, tau_s 96 min, phi 0.80, tau_f_wet 3.0, tau_f_dry 10.2.
- **Was S2 trapped on the real data? Hardly.** Training loss S2b - S2, same model and data: -0.1 to -10 for M1, M3, M4
  (losses 2e4-1e5); M3D -10 to -125. The radar-input numbers of S2 stand (changes of a few %). The trap showed on the
  real data only in the GAUGE-input F at light rain, a flat direction: M3 F[28](0.5) 0.14 -> 0.64 dB, M4 0.17 -> 0.98,
  M4 F(2) 0.76 -> 2.96 for a loss change of 2-10. So v1's "with the gauge the level still rises at 5 mm/h" was mostly
  the trap; the S2 bootstrap intervals (0.12-0.15) did not hold the new values.
- **The ladder** (held-out RMSE, dB): S<1 M1 0.403, M3 0.396, M3D 0.396, M4b 0.393 (M4 ref. 0.393): best M4b,
  **chosen M3** (M3 - M4b +0.0029, -0.0005..+0.0061). S<3: M1 0.458, M3 0.449, M3D 0.449, M4b 0.444: chosen M4b
  (M3 +0.0047, +0.0021..+0.0076). Gauge links, GAUGE input: chosen M3D (best M4b). The longer folds raised the held-out
  RMSE of almost every model by 0.001-0.003 dB (early stopping at 40 iterations had acted as a mild regularisation); this
  is why the S<1 choice moved from M3D (S2) to M3.
- **The bounded store did not fix the leftover.** tau_s of M4b is at the bound 240 min on every set. After the last
  shower (S<1) M4b stays too high (0.030 of the peak at 120 min against 0.011; integral 0-180 min 2.07 x), M3 ends too
  early (0.009 at 60 min against 0.020; integral 1.05 x on S<1, 0.55 x on S<3). RMSE 60-180 min is 0.239-0.241 dB for
  every rung (noise). The slow leftover (2-4 % of the peak, 1-2 h) is still not described.
- **Length rule again (read on M3): AGREE** (every key number of S<3 inside the S<1 interval; the intervals are wide).
  For M4b tau_f_dry is outside (4.9 against 5.6..9.7 min). M3's P halves from S<1 to S<3 (P[28](2) 1.31 -> 0.57 dB/km,
  intervals apart): with the radar, M3's P holds something that does not scale with L.
- **Convergence.** 13 of 15 real fits on all events converged (not: M3D on S<1 and S<3, 0 of 8 starts); 68 of 120 starts;
  every start needed more than the S2 cap of 200; 356 of 680 folds still reached the new cap of 120. Spread over the 8
  starts below 0.4 loss units except for M3D.
- **WAA v2 (M3, S<1 radar | GAUGE input; 90 % intervals).** A = W + L P_band(R(t - d)) + c, dW/dt = (a_link F(R) - W) /
  tau. d +2 | -3 min; tau_wet 43 [15, 274] | **2.9 [2.5, 3.6] min (the wetting speed)**; tau_dry 24.7 [14.2, 115] |
  12.8 [11.3, 15.5] min; F[28](0.5 / 2 / 5) 0.83 / 1.40 / 1.53 | 0.64 / 1.97 / 2.79 dB; F[38] 0.88 / 1.81 / 1.81 | 0.81 /
  2.48 / 3.47 dB; P[28](2 / 10) 1.31 / 3.11 | 0.39 / 0.54 dB/km (ITU size at 2 mm/h: 0.38); a_link 10-90 % 0.65-1.32 |
  0.80-1.24. Runner-up within the interval: M4b (gauge: tau_f_wet 2.6, tau_f_dry 10.3, phi 0.72, tau_s at the bound).
- **Limits that remain.** The F / P split depends on model and input (P[28](2): M3 radar 1.31, M4b radar 0.13, M3 gauge
  0.39 dB/km); the slow leftover; peaks too low (-0.83 dB radar, -0.70 gauge); drying gauge 12.8 against radar 24.7 min;
  bootstrap intervals (warm starts, 30 iterations) too narrow along flat directions; the judge depends on the fold cap.
- **For S3 (the learners decide):** apply v2 (M3) and M4b as fitted to VAL, to both directions of each link, to the
  gauge links of VAL, and to the links above 3 km (does A - W behave like L P(R) or (L + L0) P?); the F / P split is
  decided on the long links.
- Report: `outputs/reports/system_id_s2b.md`; product for S3: `outputs/products/system_id_s2b_link_events.csv`; figures
  `system_id_s2b_ladder.png`, `_F_and_P.png`, `_after_last_shower.png`, `_length_rule.png`.

### S3, verification of WAA v2: model-free length test, VAL, directions, gauge (`system_id_s3.py` + library `system_id.py`, 2026-10-15)

- **Question.** Do the S2b fits of M3 and M4b (TRAIN, radar input, S<1 and S<3) hold up on a test they were not fitted to, on the VAL events
  (read here for the first time), on both directions of a link, and with the gauge input on VAL? The S2b fits are FROZEN; the only refit is S3.4.
  Every rule is in config `wet_antenna.system_id_s3`, written before VAL was read. Standing rule: WAA is never estimated by subtracting a
  power-law rain term; the path part is separated only through L -> 0 at matched rain; ITU-R P.838 only as a size line on a figure.
- **Data.** The S2 rules for ALL links of the TRAIN link split (no length limit; bands 28 / 32 / 38 GHz): 9780 TRAIN link-events (51 events,
  206 links, up to 8.7 km), 3250 VAL (18 events). Cut at 3 km the TRAIN set is exactly the S<3 set of S2b. Library change: `build_dataset`
  takes `event_split` (default "train", so S1 / S2 / S2b are unchanged).
- **S3.1, the model-free length test (main).** Steady-rain 5 min bins (radar path-mean R in a class 0.5-1 / 1-2 / 2-5 / 5-10 / 10-20 mm/h, wet
  for the previous 20 min, previous 15 min mean within x0.5..x2); per link-event and class the median (A - c); per band and class a Huber line
  (A - c) = I + s L with an event bootstrap (1000). I = antenna part, s = path part per km. Verdict per fit: "consistent" if in >= 2/3 of the
  eligible cells (>= 8 links under 1 km; 10 cells, 28 and 38 GHz) a F lies inside the I interval AND P inside the s interval.
  **Result: no fit is consistent** (M3 S<1 0 / 10, M3 S<3 0 / 10, M4b S<1 2 / 10, M4b S<3 2 / 10; the same with R = 15 min mean, with bins
  >= 30 min into the rain, and on VAL). M3: the level a F agrees with I (6 / 10 and 9 / 10 cells) but P is 3-8 x (S<1) and 1.6-3 x (S<3) the
  slope s (28 GHz, 0.8 mm/h: P 0.64 against s 0.08 [0.05, 0.12] dB/km): M3 hides antenna wetting in P, as S2b warned. M4b: P is closer
  (inside in 4-6 cells) but near 0 at light rain, and a F at equilibrium is 1.2-1.5 x I; read through the test (the same fit on the M4b
  prediction itself, information only) its intercept matches in 10 / 10 cells, its slope in 5 / 10.
  **I and L0.** I rises and flattens (28 GHz 1.0, 1.3, 1.9, 2.5, 2.8 dB from 0.8 to 12 mm/h; 38 GHz 1.1 ... 2.4 dB) while s keeps growing
  (28 GHz 0.08 ... 1.44 dB/km), so L0 = I / s falls from about 12 km to 2 km (28 GHz; 38 GHz 7 to 0.9 km): a saturating antenna part of 2-3
  dB, NOT a constant extra length (f = L + L0 does not describe the data). Quadratic term small; links >= 1 km alone give nearly the same I.
- **S3.2, VAL.** Frozen fits from the rain alone: overall RMSE not worse than TRAIN leave-one-event-out (L < 1 km: M1 0.378, M3 0.365,
  M4b 0.360; L < 3 km: 0.374 / 0.366 / 0.363 dB), same order, gaps of thousandths; during rain worse on the short links (M3 0.76 against
  0.57 dB). Peaks too low again (-0.7 to -0.9 dB). After the rain M3 0.48 x, M4b 0.84 x the observed integral (L < 1 km). The S<1 M3 fit
  gets worse on 1-3 km links (0.450 dB), M4b does not (0.371): the S3.1 finding on unseen data. S3.1 on VAL: no fit consistent; the VAL
  intercepts and slopes are lower than TRAIN.
- **S3.3, directions.** a_link is one per link in S2b, so it cannot be compared (stated in config before the run). VAL pairs (965, 113
  links): peak of A - c Spearman 0.96 (the same by the rule >= 0.6); after-rain integral 0.34 (not the same by the rule; mostly noise, only
  380 pairs positive in both directions). Substitute amplitude (observed / predicted peak per sub-link, TRAIN): Spearman 0.78 (M3), 0.63 (M4b).
- **S3.4, gauge on VAL (the only refit).** 625 link-events, 42 links, 14 events; 8 starts, 100 bootstrap refits, nothing cut (8 min).
  M3 tau_wet 1.8 [1.5, 2.9] min (expected 2-5: just below, intervals overlap), tau_dry 10.8 [6.7, 20.6] (expected 8-18: yes); M4b
  tau_f_wet 2.0 [1.5, 3.2], tau_f_dry 8.8 [5.8, 13.1]; lag -3 min again; tau_s at 240 again. The M3 level F with the gauge is lower on VAL
  (28 GHz F(2) 1.38 against 1.97 dB).
- **WAA v3 candidate: neither.** What stands from v2: the gauge time constants (wetting a few minutes, drying 9-13 min, replicated) and the
  size of the antenna level (I about 1 dB at light rain to 2.5 dB above 5 mm/h). What fails: the F / P split of both models (P of M3 too large,
  P of M4b too small at light rain). S3.1 gives the model-free I(R) and s(R) per band that a v3 must reproduce.
- **Unclear.** The radar path mean on long links (uneven rain, b < 1) may bias s low and I up; I extrapolates the long links; the rule compares
  a F at equilibrium after only 20-30 min of rain (harsh on M4b); I also holds anything that does not grow with L (baseline drift, dew, radar
  bias); the antenna level differs between TRAIN and VAL events; 32 GHz has no link under 1 km.
- Report: `outputs/reports/system_id_s3.md`; products `system_id_s3_steady_link_events.csv`, `system_id_s3_direction_pairs.csv`; figures
  `system_id_s3_I_s_L0.png`, `system_id_s3_length_<band>.png`, `system_id_s3_val_errors.png`, `system_id_s3_val_after_last_shower.png`,
  `system_id_s3_directions.png`.

### S3b, the S2 frame on ALL link lengths so that P is pinned by L (`system_id_s3b.py` + libraries `system_id.py`, `length_test.py`, 2026-10-15)

- **Question.** S3 found that from links < 1 km the F / P split is not identifiable. Does fitting the same frame on the TRAIN events of
  ALL link lengths (9780 link-events, up to 8.7 km; radar input; the S2b settings) pin P, and do the fits then pass the length test on
  VAL? Rules in config `wet_antenna.system_id_s3b`, written before the run. VAL is the judge; TRAIN is information (in-sample).
- **Fits.** M3-A and M4b-A (free), M3-B and M4b-B (time constants fixed at the TRAIN gauge values of S2b: M3 2.9 / 12.8 min; M4b
  2.6 / 10.3 min, phi 0.72, tau_s 240), M1 for scale. 8 starts each, caps as S2b, no folds, 30 bootstrap refits (too narrow). One
  M4b-A start 389 s; projection 2.8 h against 5 h: nothing cut. Stopped once by the session for low memory and resumed from the cache.
  M4b-A did not converge on all events; for M3-A and M4b-A none of the 8 starts converged before the cap.
- **Library change.** The length test (S3.1) moved from `system_id_s3.py` into `length_test.py` (shared with S3b); S3 reproduces all its
  outputs byte for byte.
- **Check 1, the length test through the model.** TRAIN: M3-B 10 / 10, M4b-B 10 / 10, M4b-A 9 / 10, M1 9 / 10, M3-A 5 / 10 cells: with all
  lengths P is pinned (28 GHz P(2): M3-A 0.45, M4b-A 0.17 dB/km against 1.31 / 0.13 in S2b S<1). **VAL (the rule): every fit fails**
  (0-1 of 8 cells; M1 and the S2b fits too). On VAL the observed I and s are both about 0.6-0.7 x what every fit predicts on the same
  bins (TRAIN about 1.0): a difference between the event sets common to the antenna and the path part (e.g. a radar that reads the VAL
  rain too high), not the F / P split.
- **Check 2, VAL errors by length.** All fits within a hundredth of a dB (all lengths: M4b-A 0.398, M4b-B 0.402, M3-A 0.403, M3-B
  0.405, M1 0.407 dB); on L < 1 km the S2b S<1 fits stay 0.005 dB better; > 3 km about 0.48 dB. Gauge timing (B) gives better shower
  peaks (-0.47 / -0.51 against -0.87 / -0.75 dB) but a worse first hour after the rain.
- **Check 3, residual vs L on VAL.** FLAT for every fit and band (the L-scaling is right), but negative at every length (about -0.1 to
  -0.2 dB): too much attenuation predicted during rain on VAL, the same level difference.
- **Check 4, timing.** With P pinned the radar wetting stays slow: M3-A tau_wet 56 [17, 125] min, M4b-A tau_f_wet 13 [7, 72] min,
  against 2-3 min with the gauge: the slow radar wetting is input smoothing (fake check Ta of S2), not the F / P confusion.
- **WAA v3 candidate: none** (no fit consistent on VAL). What stands: P pinned by all lengths (TRAIN length test 10 / 10 for the B
  fits); the L-scaling flat on VAL; the radar wetting time is input smoothing; gauge timing gives better peaks. Best-checked on TRAIN:
  M4b-B. The VAL level difference goes to S4 (conditions).
- Report: `outputs/reports/system_id_s3b.md`; figures `system_id_s3b_I_s.png`, `_val_by_length.png`, `_residual_vs_L.png`,
  `_val_after_last_shower.png`.

### S4.1 + S4.2, why do events differ: per-event factors and the conditions (`conditions_s4.py` + libraries `system_id.py`, `length_test.py`, 2026-10-15)

- **Question.** Why does VAL attenuate ~0.6-0.7 x per radar mm/h? On the FROZEN M4b-B of S3b (M3-B as a check; nothing refitted).
  Rules in config `wet_antenna.conditions_s4`, written before the run; one change before the run after a smoke test of the code on 3
  events: the log-ridge 0.5 -> 0.05 per link (0.5 set the factor of a light event by itself).
- **Facts.** Split random within intensity terciles (seed 42; pipeline methods.md section 8). In the datasets 51 TRAIN / 18 VAL events
  (6 / 6 / 6 light / moderate / heavy in VAL), interleaved in time. 40 of 69 events have < 0.1 mm at the median link; the VAL length test
  rests on 5 large events.
- **Split halves (rule).** 200 random 18-event halves of TRAIN: 0.5 % (I) and 0 % (s) as low as VAL relative to all TRAIN (0.5 / 0 %
  relative to the model); half against the other 33 events (information) 1.5 / 1.0 %. **VAL is rare: more than ordinary event spread.**
- **Per-event factors.** A_pred = alpha_e W[a_link F(R / k_e)] + L P(R(t - d) / k_e) + c; W and L P simulated per k on a 122-point log
  grid (alpha only scales W), grid minimum of the Huber loss + log-ridge, 200 link-bootstrap draws inside the event. **Fake check
  (k in 0.7 / 1 / 1.4, alpha in 0.6 / 1 / 1.4, S2 noise): FAIL, 32 / 51 = 63 % (rule 80 %)**; heavy events 16 / 16, light 10 / 18, moderate
  6 / 17: with little rain at the links k and alpha trade against each other. **So the real k_e, alpha_e are reported, not interpreted**
  (M4b-B and M3-B give the same factors, Spearman 0.96; TRAIN vs VAL medians k 1.03 / 1.23, alpha 0.83 / 0.75, Mann-Whitney p 0.27 /
  0.42). Both factors per event pass the VAL length test in-sample (ratios 1.01 / 0.99): expected from two free numbers per event.
- **Radar vs gauge (factor-free).** On events with >= 1 mm per gauge link the radar reads 1.71 x the City gauges in VAL (pooled, 5
  events) against 0.89 x in TRAIN (14 events); pixel at the gauges 1.61 / 0.91. TRAIN / VAL = 0.52, the size of the VAL gap. Spearman
  k_e vs the ratio 0.65 (expected >= 0.5), not interpreted (rests on k_e).
- **Conditions (S4.2, descriptive).** Data in the repo: the baseline's `dew_risk_hourly.csv` (GbgA hourly), the events table, the
  datasets. 28 Spearman tests of k_e, alpha_e: 0 discoveries by BH at q = 0.1 (and not interpreted). Residual of the frozen M4b-B
  during rain (no factors): -0.08 dB at every band, day and night alike; warmest tertile (>= 17.3 C) -0.15 against -0.02 dB at < 15.8 C
  (between events, not shown to be temperature).
- **Candidate for WAA v3 (not final):** keep M4b-B as is; no event factor; no condition-dependent antenna level. The VAL gap most likely
  comes from the rain reference (radar over-reading the large VAL events), not from the antennas. Proposal for a later sub-step: a fake
  check restricted to heavy events, rule fixed before running.
- Report: `outputs/reports/conditions_s4.md`; figures `conditions_s4_timeline.png`, `_split_half.png`, `_fake.png`, `_factors.png`,
  `_close_gap.png`, `_radar_gauge.png`, `_conditions.png`, `_residual_check.png`.

### S4.3, is the VAL gap the radar? (`conditions_s4c.py` + libraries `event_factors.py`, `system_id.py`, `length_test.py`, 2026-10-17)

- Rules in config `wet_antenna.conditions_s4c`, written before the run (2026-10-15). Frozen M4b-B (main) and M3-B (check), no refit.
  Standing rule: WAA is never estimated by subtracting a power-law rain term; the rain correction comes from the City gauges only.
- **g_e** = radar / City-gauge rain per event (the S4.1 definitions, moved unchanged to `event_factors.py`; S4.1 tables identical after
  the move): path mean of the gauge links (used for 19 events: 14 TRAIN, 5 VAL; >= 1.0 mm per gauge link) and pixel at the 10 gauges
  (24 events); g_e = 1 otherwise. Pooled: G_TRAIN 0.894, G_VAL 1.707. Variants V0 (none), V1 (R / g_e path), V2 (R / g_e pixel),
  V3 (R / G per set), applied to the model input and to the observed steady-rain bins (classes recomputed).
- **Length test (S3b check 1).** VAL V0: I 0.72, s 0.61, both inside 1 / 8; V1: 0.99, 0.95, 5 / 6 -> **explained**; TRAIN V1: 1.11,
  0.92, 7 / 10 -> **not broken**. V2 and V3 agree. M3-B misses the 2/3 share rule (4 / 6 against the written 0.6667; TRAIN 0.6).
  Model-free VAL / TRAIN median: I 0.63 -> 0.84, s 0.67 -> 0.91; 28 GHz 2-5 mm/h stays at 0.53.
- **Heavy events.** Fake check (16 heavy TRAIN events x 3 replicates): 98 % recover k and alpha within +-25 % -> PASS. Real factors
  (22 events, S4.1 fits): Spearman(ln k_e, ln g_e) 0.64 [0.27, 0.89] (18 events; expectation >= 0.5 met); alpha_e against 6
  conditions: 0 BH discoveries (q = 0.1); strongest peak / mean rain rho -0.49 (p 0.02). Descriptive, not a law.
- **Decision:** WAA v3 = M4b-B unchanged. S5 evaluates against the raw radar and the g_e-rescaled radar (main), gauge links as the
  independent check. Open: the radar is gauge-adjusted (circular), only 5 VAL events carry g_e, the 5-10 mm/h VAL cells leave the test.
- Report: `outputs/reports/conditions_s4c.md`; figures `conditions_s4c_ratios.png`, `_observed_curves.png`, `_heavy_fake.png`,
  `_heavy_conditions.png`.

### S5a, WAA v3 refitted on gauge-corrected rain, the slow leftover, closed forms (`final_s5a.py` + libraries `system_id.py`, `length_test.py`, `event_factors.py`, 2026-10-17/18)

- Rules in config `wet_antenna.final_s5a`, written before the run. Standing rule as S4.3; closed forms are fitted to the learned curves,
  a power law only as one candidate. TEST never read. Library changes (old outputs byte-identical, checked by re-running S3b and S4.3):
  M4c in `system_id.kernel2` (a slow store with separate fill / empty times behind a flag), `fix_params`, `x_into`, `fake_from_truth`;
  `val_metrics` moved to `length_test.py`, `rescaled` to `event_factors.py`.
- **Part 1.** M4b-B refitted on R / g_e (TRAIN, all lengths, 8 starts; training loss 163 991 -> 158 513); 30 COLD event-bootstrap refits
  (much wider intervals than the warm S3b ones). Curves within the intervals of S3b except F(1) at 32 / 38 GHz. VAL: V1 better than
  V0 for every fit (M4b-B 0.388 against 0.397 dB); the VAL length test passes on V1 for every fit (5 / 6), on V0 for none; shower peaks
  worse on V1 (-1.05 against -0.53 dB); M4b-B's after-rain integral already inside 0.7-1.4 on V1 (0.79 / 1.21).
- **Part 2.** M4c fake check PASS (tau_s_wet 29.9 / 26.6, tau_s_dry 95 / 103 against 30 / 90). Real fit: tau_s_wet 600 (bound), tau_s_dry
  217 [10, 600], phi 0.43 [0.19, 0.84]: **not identified**. By the rule M4c replaces M4b-B (integral 0.89 / 1.29, RMSE -0.0005 dB, length
  test 5 / 6). The median after-rain curves stay too high from ~40 min for both models (M4c more).
- **Part 3.** Knots 0.5-40 mm/h (723 TRAIN link-events above 20 mm/h) with cold SE. Chosen by AICc: P = c R (28, 32 GHz), c R^b with
  b 0.88 [0.81, 0.95] (38 GHz), quadratic term 0 everywhere; F = F_max R / (R + R0) per band (F_max 5-7.5 dB, wide intervals), the capped
  power law pooled. PySR: P = 0.147 R (= the linear candidate), F = log(R^2 + 3.46) (not a candidate, no saturation). Plug-back holds
  (VAL V1 RMSE +0.002 dB, length test 5 / 6): **WAA v4 (closed form)**, box in the report.
- **Open:** M4c kept by a rule without an identifiability criterion; the peaks on V1; F weakly determined; the S5b model choice (M4c or
  M4b-B) must be settled before TEST is read.
- Report: `outputs/reports/final_s5a.md`; figures `final_s5a_F_and_P.png`, `_val.png`, `_after_rain.png`, `_m4c_fake.png`.

### S6, tighten the equation of M4b-B on gauge-corrected rain (`tighten_s6.py` + library `final_model.py`, 2026-10-18)

- Rules in config `wet_antenna.tighten_s6`, written before the run; integrity rules in `docs/ROADMAP.md`. Decision before S6: the
  final model is M4b-B (M4c dropped: not identified; decided after VAL, for identifiability). VAL numbers are "used (optimistic)".
  Shared code of S5a moved unchanged to `final_model.py` (S5a re-run: tables identical except 5th-digit noise in one unchosen candidate).
- **Funnel.** 728 sub-links -> 660 in the 28 / 32 / 38 GHz bands -> 518 on TRAIN links; 29 998 TRAIN (event, sub-link) pairs -> 9 780
  link-events (the largest loss: no wet bin of the link in the event, 13 480); VAL 3 250. 12 % of the TRAIN minutes have V1 rain
  >= 0.5 mm/h, 0.06 % above 20 mm/h. TEST: 64 links, 18 events (counted only).
- **Shared shapes** F = A_band G(R), P = B_band H(R) (G(5) = H(5) = 1): VAL V1 0.3875 against 0.3878 dB, length test 5 / 6: **kept**.
  TRAIN AICc prefers the band shapes (586 AICc; still ~80 after the autocorrelation factor): the band details do not generalise.
- **Error bars.** 60 cold event-bootstrap draws (all converged; the 60 -> 40 budget cut was undone, dated deviation). A = 2.74 [2.20, 3.37],
  2.56 [1.88, 3.39], 2.96 [2.49, 3.53] dB; B = 0.57 / 0.75 / 1.09 dB/km at 5 mm/h. The profile (threshold x kappa 7.4) is 5-13 x narrower
  for shapes and B (it ignores the correlation between the links of one event) but WIDER for A: A trades against the mean a_link, only the
  ridge prior pins it, so **the data fix a_link F, not F**. Timing: tau_f_dry +-14 % moves F by 4-11 % and P by 5-19 % (opposite signs).
- **Identifiability.** Fake from the kept model: 62 of 62 quantities within +-20 % in 3 of 3; P (and H) at 0.5 and 1 mm/h NOT identified
  (interval at 0).
- **Closed forms.** Physics and PySR compete by AICc (simpler within 2 wins). P: the line c R (c = 0.107 / 0.142 / 0.197 dB/km per mm/h);
  F: PySR's (R / 5)^0.47 (the capped power +0.19 with one more parameter, the hyperbola +11). Plug-back: P alone holds; F (alone or
  with P) lowers VAL RMSE but the length test drops to 4 / 6 = exactly 2 / 3, failing the written 0.6667: F stays free (shared) knots.
- **Open:** 0.6667 or 2 / 3; the F level and the a_link prior; F above 10 mm/h; P below 1 mm/h; the slow store (fixed). S7 checks.
- Report: `outputs/reports/tighten_s6.md`; product `outputs/products/tighten_s6_candidate.json`; figures `tighten_s6_funnel.png`,
  `_F.png`, `_P.png`, `_shared_shape.png`, `_fake.png`, `_intervals.png`.
- **Decision 2026-10-19 (before any TEST read):** 4 / 6 = 0.66667 against the written 0.6667 is a rounding of the written number (the S3b
  wording is "at least 2 / 3"); the candidate is NOT changed; every length-test threshold is now written "at least 4 of 6 cells".

### S7, validate the frozen S6 candidate without TEST (`validate_s7.py` + library `final_model.py`, 2026-10-19)

- Rules in config `wet_antenna.validate_s7`, written before the run. The frozen object is `tighten_s6_candidate.json` (sha256 in the
  report); nothing of it is refitted. Two F descriptions are carried: the candidate (free shared knots) and the closed F (R / 5)^0.47.
  The S6 parametrisations moved unchanged to `final_model.py` (S6 re-run: every table byte-identical).
- **7.1 W against L: pass.** a_link against L rho +0.05 [-0.07, +0.17]; W height a_link F at 5 mm/h 2.81 / 2.75 / 2.81 dB (L < 1, 1-3,
  > 3 km), intervals overlap. No metadata variable is a BH discovery (antenna diameter and orientation nominal, |rho| < 0.2).
- **7.2 short links on VAL.** Length test 5 / 6 (closed F 4 / 6); shared shapes as good as band shapes on short links. The written RMSE
  rule FAILS for L > 3 km (1.18 x the all-class VAL RMSE; 1.12 x on TRAIN; every VAL class below TRAIN): reported, not a return reason.
- **7.3a dry-path moments** (path-max radar <= 0.05 mm/h on every pixel, 9 081 windows): model error at 10 min -0.18 dB in every length
  class (difference L > 3 minus L < 1: 0.00 [-0.04, 0.03] dB): **pass**. The model is 0.2-0.3 dB high right after the rain, about
  0.5 dB high after heavy rain; long links keep more water late. **7.3b** paired links: 4 pairs, not testable. **7.3c** intercepts from
  short-only and long-only links overlap in 13 / 15 cells: **pass**. **7.3d** the two F descriptions cannot be told apart. **7.3e**
  c_band R is 0.52-0.74 x ITU-R P.838 (size reference only, no flag).
- **7.4 dew** (baseline `dry_jumps.csv` with the baseline dew rule; 125 link-nights): dew W 1.7 dB (peak 2.3) for about 4 h after
  midnight, larger than the rain W at 0.75 mm/h (1.0 dB); weakly tied to a_link (rho +0.26 [+0.06, +0.45]); dew before the rain raises
  c by 0.13 [0.05, 0.35] dB (8 events).
- **7.5 leftovers:** a_link: no relation found (best PySR form just below the null's 95 %). Per-event drying time (17 events): falls with
  the wind speed (Spearman -0.76; PySR exp(-0.32 u) form, leave-one-event-out skill 0.50 > null 95 % 0.34): **a relation (hint)**.
- **7.6 rain reference:** the radar / gauge ratio moves in opposite directions in the two binnings (gauge -57 %, radar +88 %):
  regression to the mean, **no rate dependence**; the parked g(R) is not justified.
- **Return decision:** only reason (d) is met (the wind relation) -> the one allowed return to S6, recorded 2026-10-19; not started here.
- Report: `outputs/reports/validate_s7.md`; figures `validate_s7_alink_vs_L.png`, `_w_height.png`, `_val_rmse.png`, `_dry_path.png`,
  `_decay.png`, `_pairs.png`, `_intercept.png`, `_dew.png`, `_rain_reference.png`.

### S6b, the one return to S6: drying, wind and the after-rain misfit (`return_s6b.py` + libraries `wind_drying.py`, `validation.py`, 2026-10-19)

- Rules in config `wet_antenna.return_s6b`, written before the run; the three boxes (A in the equation, B identified outside, C
  unexplained) and the success definition in `docs/ROADMAP.md`. One dated deviation: the projected cuts (bootstrap 40 -> 30, M3) were
  undone before Part 2 ran. S7 helpers moved unchanged to `validation.py`; `Curves` takes a problem and free dynamic numbers (S6 and S7
  re-run identical; the S7 pair table now in a fixed order).
- **Part 1, the after-rain misfit (M0 frozen, 9 081 dry-path windows):** (a) not the radar's timing (the radar ends a median 9 min after
  the gauge, but the error does not depend on the lag); (b) not a faster drying in general (decay 0-10 min measured / model 0.87
  [0.70, 1.06]); the error is already -0.25 [-0.35, -0.16] dB when the path becomes dry: a LEVEL at the end of the rain; (c) it grows
  with the preceding peak rain (partial rho -0.13 given the wind) and in dry air (+0.24); the wind does not survive the preceding rain held
  fixed; (d) F is not too high at heavy rain on short links (VAL even low, 4 events).
- **Part 2:** kernel `wind_drying.kernel_w` (exact M4b at k = 0; gradients checked). M1 tau_f_dry = 3.7 min; M2 tau_0 3.0 min, k +0.83
  per SD of the wind; M3 k -0.009 (none). Every extra number recovered on fake data (3 of 3).
- **Part 3, the box-A rule:** M1 passes all (VAL after-rain error -0.24 / -0.17 -> -0.07 / -0.05 dB at 5 / 10 min; RMSE -0.0024 dB;
  length test 4 of 6; 23 / 23 folds). M2 fails rules 3 (heavy) and 4 (length classes): box B. M3 wrong sign: box C. **But M1 halves
  c_band and its bootstrap reaches 0 in every band (8-11 of 40 draws): the path term is NOT identified with M1**, and 3.7 min is an
  effective time for the radar input, not the gauge-measured drying (10.3 / 8.8 min).
- **Error budget (VAL):** dry minutes 30 % (0.40 dB), during rain 68 % (1.22 dB; between events only 3 %), after rain 2 % (0.31 dB).
- **Candidate for S8:** M1 by the rules (`outputs/products/return_s6b_candidate.json`); a learners' decision is needed (M1, or the S6
  candidate with the faster drying in box B).
- Report: `outputs/reports/return_s6b.md`; figures `return_s6b_after_rain_causes.png`, `_fake.png`, `_dry_path.png`, `_boxes.png`,
  `_error_budget.png`.

### S8, freeze, then TEST once (`test_s8.py`, 2026-10-19)

- Decision 2026-10-19 (before TEST): the candidate is the S6 candidate; M1 of S6b goes to box B (it unbinds the path term). Lesson: a
  refinement rule must also require that the other numbers stay identified.
- Comparison models fitted on TRAIN V1: B0 no WAA, B1 a constant WAA while it rains, B2 WAA without memory (static, the S6 shapes); B3
  skipped (no literature WAA model is implemented in the repo). VAL (used): candidate 0.620, B0 0.717, B1 0.646, B2 0.630 dB.
- Freeze: `outputs/products/test_s8_freeze.json` (hashes of the candidate, the B fits, the rules; code commit 761c6f1), committed in
  bb30b5c before TEST was read. The evaluation code was first run on VAL (`S8_DRYRUN=VAL`).
- TEST: baselines for the TEST events computed with the frozen baseline code (TRAIN and TEST links); 3 355 link-events on TRAIN links,
  1 010 on TEST links; g_e from the TEST events' gauge links (4 events used). Two evaluation bugs fixed after the read (a zero radar /
  gauge ratio; events without a g_e entry), both giving undefined rain; recorded with the numbers before and after.
- **P2-P5 pass**: TEST RMSE 0.692 dB against 0.862 / 0.760 / 0.727 (B0 / B1 / B2), L < 1 km 0.562 against 0.838 / 0.644 / 0.601; dry-path
  error equal across lengths; TEST / VAL 1.12 in every length class within 5 %; TEST links 1.07 x TRAIN links.
- **P1 fails** (4 of 8 cells): the slopes fit in 7 of 8, the intercepts at 0.5-2 mm/h are higher on TEST (1.35-1.55 dB against
  0.83-1.11 dB in the model): the antenna level at light rain, not the dependence on L.
- Report: `outputs/reports/test_s8.md`; figures `test_s8_rmse.png`, `_length_cells.png`, `_dry_path.png`, `_test_over_val.png`,
  `_new_links.png`, `_error_budget.png`.

### S8b, post-freeze evaluation of a pre-freeze model: the closed form of F (`test_s8b.py`, 2026-10-19)

- The S6 plug-back fit "F and P" (F = A_band (R / 5)^b, P = c_band R; cache file modified 2026-10-18 12:45, before the S8 freeze). Rule
  (written after TEST, before this run): adopted as the written form of F iff the same P1-P5 verdicts as the table F and the paired TEST
  RMSE difference within [-0.01, +0.01] dB (all and L < 1 km).
- TRAIN: 40 cold draws, all 9 numbers identified (b 0.54 [0.47, 0.63]; A 2.57 / 2.36 / 2.77 dB; c 0.105 / 0.141 / 0.194 dB/km per mm/h);
  inside the table F's interval at 21 / 21 knots; +2.3-2.7 dB above it at 40 mm/h (no saturation).
- TEST: the same verdicts (P1 fails, P2-P5 pass); paired difference +0.0025 [-0.0065, +0.0085] dB (L < 1 km +0.0016 [-0.0082, +0.0078]).
  **The closed form is adopted** as the written form of F (valid 0-40 mm/h). S8's comparison-model setup moved to `validation.py`
  (S8 re-run identical).

### S9, formulate Part 1: the official equation and the documents (`formulate_s9.py`, library `waa_equation.py`, 2026-10-19)

No new analysis, no refitting, no new reading of TEST. Rules in config `wet_antenna.formulate_s9`, written before the run.

- **The final equation** (the closed-F model of S8b; FINAL.md): A(t) - c = W_f + W_s + L c_band R(t - 1 min), dW_f/dt = (0.72 u - W_f) /
  tau_f (2.6 min wetting, 10.3 min drying), dW_s/dt = (0.28 u - W_s) / 240 min, u = a_link A_band (R / 5)^b; b 0.542 [0.472, 0.626],
  A 2.57 / 2.36 / 2.77 dB, c 0.105 / 0.141 / 0.194 dB/km per mm/h (28 / 32 / 38 GHz), a_link 10-90 % 0.77-1.24; R = radar / g_e; valid
  0-40 mm/h. Stored in `outputs/products/waa_equation_final.json` (every number with its 90 % interval, the fit's sha256
  `1d8163cd...e0b9`, its own sha256 in the S9 report).
- **The code** `waa_equation.py` steps the stores once per minute exactly as `system_id.kernel2` and evaluates F on the fit's 61-knot
  grid by default (`exact=True` gives the pure formula; they differ by at most 0.05 dB, near R = 0 where R^0.54 is steep, and by
  < 0.004 dB from 0.5 to 40 mm/h). The unit test
  (`tests/test_waa_equation.py`) reproduces the S8b W and path term of 3 TEST link-events (one per length class, from the S8 cache)
  within 1e-6 dB (achieved 8.9e-16 dB).
- **The overview check** (`outputs/reports/formulate_s9_check.md`): the overview's LaTeX source was not in the repo, so
  `docs/WAA_research_overview.tex` was rebuilt from the PDF. 113 statements: 89 confirmed, 21 corrected, 2 unsupported (reworded), 1
  added. The main corrections: the lower bands were dropped because their links are all >= 4 km; the VAL timing is similar, not
  equal; the fake recovery is within 20 % (61 of 62 within 10 %); a_link 0.77-1.24; the ITU exponent 0.86-0.92; the reference [2] title.
- **The document** `docs/WAA_research_overview.pdf` (13 pages): the first overview's look and wording with the corrections, figures
  drawn by `docs/WAA_research_overview_figures/figs.py` from saved tables only (colour-free). Figures 1 and 2 of the first version
  (model illustrations) are now real data. The memory figure stays an illustration from the equation and is captioned so. (S9 first
  wrote a short and a full version in `docs/overview/`; the learners merged them into this one document on 2026-10-19.)

# 3. Open issues

**Superseded by the reboot of 2026-10-13 (`docs/ROADMAP.md`, "Reboot (2026-10-13): system identification on short
links").** The drying-tail exploration (T0-T18) is closed and the next steps are S1-S5; the issues below are kept
as history (they use the old step numbers). The open issues of the S-steps are in "Part 2: S-steps" above (S1:
the wetting time, the slow store, and the split between P and F are not identified).

- **The three questions that come before any drying law (plan change of 2026-10-12, Phase A = step 4c).** Full text,
  the later questions Q4-Q5 and the phases: `docs/ROADMAP.md`, section "Plan change (2026-10-12, after 4b): the signal
  before the law". Decided by the learners after 4c: the rises are not drift (no baseline v3); tail set v2 = v1 with
  S0 as start and S1, S2 stored per tail; **Q2 is the central open question**; the order is now 4c2 -> 3c + 7 (whole
  events on very short links) -> 4d. After step 4c:
  - **Q1. Is a night-time rise water at all?** Baseline drift with the daily cycle, dew / condensation on the radome
    (WAA without rain), or drizzle under the radar. Test: A1, the tail against the same sub-link at the same clock
    time on certified dry days. **4c: not drift for most rising tails (24 of 36 above their dry days), weaker for the
    night-time ones (10 of 18); dew or drizzle NOT decided (humidity missing for 22 of the 36).**
  - **Q2. Is the early fast drop after t0_trace the antenna drying, or path rain fading?** Test: A2, A against L per
    minute around t0_trace. **4c: OPEN. The test has no power (it does not see path rain before t0_trace either).
    4c2: narrowed. From S1 (no rain within 5 km) A halves in 40 min, but a placebo start just as late gives 37 min
    (+3, -2..+7): "SAME speed" by the rule. The fast first 15-20 min are mostly not "nearby rain ending"; whether
    they are the antenna or rain on the path that the radar no longer sees is for Step 7. New: two phases in time
    (t_half 18 min from t0_trace, about 40 min for what is left), and at the same level k is 4-7 times lower later.
    7a: UNDECIDED by its pre-set rule (no power: the radar expects almost no path rain before a tail), but the fast
    drop is the same on links of 0.5 km and of 5 km (matched SHORT / LONG ratio 0.94, 0.84..1.15; path rain would give
    about 0.2). The learners decide whether that closes Q2.**
  - **Q3. When does drying start?** The radar is 5 min / 2 km and misses drizzle, so t0_trace may be too early or too
    late. Test: A3, three start rules defined without any law fit. **4c: the references within 5 km are dry 25 min
    after t0_trace (median), when 81 % of the tails are already below 1 dB; no start rule removes the falling k(A).
    Proposal: keep t0_trace and model the start per tail (learners decide).**

- **Step 4b: why tails stall is not known, only where it does not come from.** Not the level (in decaying tails k does
  not fall with A), not the link (p = 0.35), not the path length; it is shared inside an event. The rain references
  show nothing (the path radar is dry by construction of v1), humidity is missing for the two events with most
  rising tails, and the weather is one hourly station. Also open: is the fast fall of the "decays" tails the end of
  the path rain; and the noise multiplier of `pooled_fit.py` did not pass its fake check (4 of 6).

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
