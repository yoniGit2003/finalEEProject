"""Where the wet_antenna module reads and writes its files: the ONLY place with file locations of this module.

    outputs/products/  data that later steps read (the tail sets, the simulated tails)
    outputs/tables/    all other tables; the file name starts with the script that writes it
    outputs/figures/   all figures, same naming
    outputs/reports/   the report of each script (<script>.md)
Everything in outputs/ is committed, except the few big files named in the root .gitignore.
The cache (fits that take hours, logs) is research/.cache/wet_antenna/: not in git, see pipeline/get_caches.py.
"""
from pathlib import Path

MODULE = Path(__file__).resolve().parent                   # research/wet_antenna
REPO = MODULE.parents[1]                                   # the repository root
CONFIG_FILE = REPO / "research" / "config.yaml"            # all research settings (section `wet_antenna`)

OUTPUTS = MODULE / "outputs"
PRODUCTS = OUTPUTS / "products"
TABLES = OUTPUTS / "tables"
FIGURES = OUTPUTS / "figures"
REPORTS = OUTPUTS / "reports"
CACHE = REPO / "research" / ".cache" / "wet_antenna"

# ---- step 0: inventory (inventory.py) -----------------------------------------------------------------------
INVENTORY_EXAMPLE_PNG = FIGURES / "inventory_example.png"
INVENTORY_EXAMPLE_CSV = TABLES / "inventory_example.csv"
INVENTORY_REPORT = REPORTS / "inventory.md"

# ---- step 1: all drying tails (extract_tails.py) ------------------------------------------------------------
TAILS_1MIN = PRODUCTS / "tails_1min.parquet"               # one row per tail minute (large: git-ignored)
TAILS_INDEX = PRODUCTS / "tails_index.csv"                 # one row per tail
EXTRACT_A_CHECK = TABLES / "extract_tails_a_check.csv"               # A(t) built here vs baseline.a_of_t on 3 pairs
EXTRACT_RETURN = TABLES / "extract_tails_return_to_dry.csv"          # share of tails back below the threshold after k min
EXTRACT_MEDIAN_CURVES = TABLES / "extract_tails_median_curves.csv"   # median A per minute and A_start bin
EXTRACT_RISES = TABLES / "extract_tails_dry_rises.csv"               # one row per tail with a later dry-weather rise
EXTRACT_FIG_DISTRIBUTIONS = FIGURES / "extract_tails_distributions.png"
EXTRACT_FIG_MEDIAN_CURVES = FIGURES / "extract_tails_median_curves.png"
EXTRACT_FIG_RESOLUTION = FIGURES / "extract_tails_resolution_panels.png"
EXTRACT_REPORT = REPORTS / "extract_tails.md"

# ---- step 2: selection and fast-drop test (select_tails.py) ------------------------------------------------
TAILS_SELECTED = PRODUCTS / "tails_selected.csv"           # index + exclusion, truncation, main set, shape descriptors
SELECT_COUNTS = TABLES / "select_tails_counts.csv"                 # main-set counts per band x length class x antenna
SELECT_CORRELATIONS = TABLES / "select_tails_correlations.csv"     # Spearman rho with link-bootstrap 95 % CI
SELECT_TRUNCATION = TABLES / "select_tails_truncation.csv"         # tail time removed by each truncation rule
SELECT_MEDIAN_CURVES = TABLES / "select_tails_median_curves.csv"   # median A per minute: main set, control, margin-rain split
SELECT_FIG_LENGTH = FIGURES / "select_tails_fast_drop_vs_length.png"
SELECT_FIG_MARGIN = FIGURES / "select_tails_margin_rain.png"
SELECT_FIG_PAIRS = FIGURES / "select_tails_sublink_pairs.png"
SELECT_FIG_SHAPE = FIGURES / "select_tails_shape.png"
SELECT_FIG_MEDIAN = FIGURES / "select_tails_median_curves.png"
SELECT_FIG_GALLERY = FIGURES / "select_tails_gallery.png"
SELECT_REPORT = REPORTS / "select_tails.md"

# ---- inputs from the baseline module (frozen v2 products; read only) ----------------------------------------------
BASELINE_MODULE = REPO / "research" / "baseline"

# ---- step 2b: fixes after the step 2 gallery (refine_selection.py; adds columns to TAILS_SELECTED) -----------------
REFINE_COUNTS = TABLES / "refine_selection_counts.csv"               # strict / lenient main-set counts per band x length x antenna
REFINE_NOISY = TABLES / "refine_selection_noisy_sublinks.csv"        # sub-links removed by the offset_std rule
REFINE_TIME = TABLES / "refine_selection_tail_time.csv"              # tail time kept / lost by each step 2b change
REFINE_MEDIAN_CURVES = TABLES / "refine_selection_median_curves.csv" # median A (and A / A_start) aligned at t0_trace
REFINE_FIG_NOISE = FIGURES / "refine_selection_offset_std.png"
REFINE_FIG_MEDIAN = FIGURES / "refine_selection_median_curves.png"
REFINE_FIG_UNDERSHOOT = FIGURES / "refine_selection_undershoot_vs_weather.png"
REFINE_FIG_EX_SHIFT = FIGURES / "refine_selection_examples_level_shift.png"
REFINE_FIG_EX_UNDER = FIGURES / "refine_selection_examples_undershoot.png"
REFINE_FIG_GALLERY = FIGURES / "refine_selection_gallery.png"
REFINE_REPORT = REPORTS / "refine_selection.md"
BASELINE_UNCERTAINTY = BASELINE_MODULE / "outputs" / "products" / "baseline_uncertainty.csv"
DEW_RISK_HOURLY = BASELINE_MODULE / "outputs" / "products" / "dew_risk_hourly.csv"
DRY_JUMPS = BASELINE_MODULE / "outputs" / "products" / "dry_jumps.csv"                     # dry-weather jumps (R1 C4)
DAILY_CYCLE = BASELINE_MODULE / "outputs" / "tables" / "c3_daily_cycle_per_sublink.csv"    # daily-cycle amplitude (R1 C3)

# ---- step 2c: last fixes and the frozen tail set (freeze_tail_set.py); the set id (v1) is in config wet_antenna.freeze -------
TAILS_MODEL = PRODUCTS / "tails_model.parquet"             # frozen main set, 1 min rows from t0_trace (train + val)
TAILS_MODEL_INDEX = PRODUCTS / "tails_model_index.csv"     # one row per frozen tail
TAILS_MODEL_ALT = PRODUCTS / "tails_model_alt.parquet"     # same rules, A_start_trace >= 0.6 dB
TAILS_MODEL_ALT_INDEX = PRODUCTS / "tails_model_alt_index.csv"
FREEZE_EVENTS = TABLES / "freeze_tail_set_events.csv"               # one row per TRAIN event of the frozen main set
FREEZE_CHANGES = TABLES / "freeze_tail_set_changes.csv"             # tails / time changed by each step 2c rule
FREEZE_MEDIAN_CURVES = TABLES / "freeze_tail_set_median_curves.csv"
FREEZE_FIG_MEDIAN = FIGURES / "freeze_tail_set_median_curves.png"
FREEZE_FIG_EVENT = FIGURES / "freeze_tail_set_event_{event_id}.png" # filled with the event id
FREEZE_REPORT = REPORTS / "freeze_tail_set.md"

# ---- step 3a (+ 3a-fix): tail simulator (simulator.py = library, calibrate_simulator.py = calibration + checks) ----------
SIM_TAILS_M1 = PRODUCTS / "sim_tails_m1.parquet"          # simulated tail set per law, 1 min rows (one file per law)
SIM_TAILS_M3 = PRODUCTS / "sim_tails_m3.parquet"          #   file)
SIM_TAILS_M4 = PRODUCTS / "sim_tails_m4.parquet"
SIM_TAILS_M5 = PRODUCTS / "sim_tails_m5.parquet"
SIM_INDEX = PRODUCTS / "sim_tails_index.csv"               # per law and tail: time-scale factor, offset, noise, phase
SIM_NOISE = PRODUCTS / "sim_noise_params.csv"              # noise parameters per sub-link (white, slow, tau_noise)
SIM_CALIBRATION = TABLES / "calibrate_simulator_calibration.csv"     # law parameters + median-curve misfit (key)
SIM_CHECKS = TABLES / "calibrate_simulator_shape_checks.csv"         # ratios, t_half, t_e, last-hour std: real vs each law
SIM_FIG_MEDIAN = FIGURES / "calibrate_simulator_real_vs_sim_median.png"   # key
SIM_FIG_CROSSING = FIGURES / "calibrate_simulator_t_half_t_e.png"
SIM_FIG_NOISE = FIGURES / "calibrate_simulator_noise_fit.png"
SIM_FIG_GALLERY = FIGURES / "calibrate_simulator_gallery.png"
SIM_THALF = TABLES / "calibrate_simulator_thalf_by_astart.csv"       # median t_half per A_start bin, real vs each law
SIM_FIG_THALF = FIGURES / "calibrate_simulator_thalf_by_astart.png"
SIM_REPORT = REPORTS / "calibrate_simulator.md"

# ---- step 3b-1: the fitting method on simulated tails (per_tail_fit.py = library, check_per_tail_fit.py) ----------------
PER_TAIL_FITS = TABLES / "check_per_tail_fit_fits.parquet"             # every per-tail fit: dataset, tail, candidate, params, BIC
PER_TAIL_CONFUSION = TABLES / "check_per_tail_fit_confusion.csv"       # true law x chosen candidate, per tail and per dataset
PER_TAIL_RECOVERY = TABLES / "check_per_tail_fit_recovery.csv"         # parameter recovery of the true law's own candidate (key)
PER_TAIL_M6 = TABLES / "check_per_tail_fit_m6_calibration.csv"         # m6 (state switch) calibration scan
PER_TAIL_THALF = TABLES / "check_per_tail_fit_thalf_population.csv"    # t_half per A_start bin: replicate range vs real
PER_TAIL_DADT = TABLES / "check_per_tail_fit_dadt_population.csv"      # dA/dt vs A: replicate range vs real
PER_TAIL_ABLATION = TABLES / "check_per_tail_fit_ablation.csv"         # confusion for m1 / m4 with effects switched off
PER_TAIL_FIG_CONFUSION = FIGURES / "check_per_tail_fit_confusion.png"  # key
PER_TAIL_FIG_THALF = FIGURES / "check_per_tail_fit_thalf_population.png"
PER_TAIL_FIG_DADT = FIGURES / "check_per_tail_fit_dadt_population.png"
PER_TAIL_FIG_ABLATION = FIGURES / "check_per_tail_fit_ablation.png"
PER_TAIL_REPORT = REPORTS / "check_per_tail_fit.md"

# ---- step 3b-1b: pooled fit on event-bootstrap simulated datasets (pooled_fit.py, check_pooled_fit.py) -------------
POOLED_FITS = TABLES / "check_pooled_fit_fits.csv"               # every pooled fit: truth, replicate, window, candidate, BIC
POOLED_RECOVERY = TABLES / "check_pooled_fit_recovery.csv"       # own-candidate recovery over replicates (key)
POOLED_DECISION = TABLES / "check_pooled_fit_decision.csv"       # the pre-set decision rule per truth and window
POOLED_POPULATION = TABLES / "check_pooled_fit_population.csv"   # t_half per bin and A / A_start at k min: ranges vs real
POOLED_FIG_CONFUSION = FIGURES / "check_pooled_fit_confusion.png"   # key
POOLED_FIG_POPULATION = FIGURES / "check_pooled_fit_population.png"
POOLED_REPORT = REPORTS / "check_pooled_fit.md"

# ---- step 3b-1c: stress test of the pooled fit (stress_pooled_fit.py; fake data only) -------------------------------
STRESS_FITS = TABLES / "stress_pooled_fit_fits.csv"               # every fit: test, truth, replicate, candidate, start, BIC
STRESS_WIN = TABLES / "stress_pooled_fit_win_shares.csv"          # test x truth x candidate: replicates won
STRESS_DECISION = TABLES / "stress_pooled_fit_decision.csv"       # the pre-set rule per test and truth (cell verdicts)
STRESS_RECOVERY = TABLES / "stress_pooled_fit_recovery.csv"       # own-candidate parameter ratios per test and truth
STRESS_STARTS = TABLES / "stress_pooled_fit_start_sensitivity.csv"   # T2: spread of the BIC over the starts
STRESS_RECAL = TABLES / "stress_pooled_fit_recalibration.csv"     # T3: old and new truth parameters, A / A_start at 10/30/60
STRESS_NOISE = TABLES / "stress_pooled_fit_noise_check.csv"       # T1: real vs simulated dry-noise statistics
STRESS_FIG_WIN = FIGURES / "stress_pooled_fit_win_share.png"      # own-candidate win share per test and truth
STRESS_REPORT = REPORTS / "stress_pooled_fit.md"

# ---- step 3b-2: weak-form SINDy on fake tails (sindy.py = library, check_sindy.py; fake data only) ----------
SINDY_FITS = TABLES / "check_sindy_fits.csv"                 # every SINDy run: dataset, window, offset, rescale, terms
SINDY_SELECTION = TABLES / "check_sindy_selection.csv"       # truth x term: replicates in which the term is selected
SINDY_COEFS = TABLES / "check_sindy_coefficients.csv"        # fitted / true coefficient of the true terms
SINDY_DECISION = TABLES / "check_sindy_decision.csv"         # the pre-set rule per truth and the verdicts Q1 to Q4
SINDY_AGREEMENT = TABLES / "check_sindy_agreement.csv"       # per dataset: pooled-fit winner vs SINDy family
SINDY_FIG_SELECTION = FIGURES / "check_sindy_selection_frequency.png"   # selection frequency of each term per truth
SINDY_REPORT = REPORTS / "check_sindy.md"

# ---- step 3b-3: free drying curve, then its equation (free_curve.py = library, check_free_curve.py) ----------
CURVE_FITS = TABLES / "check_free_curve_fits.csv"                 # one row per free-curve fit (dataset x window): knots, k, s, BIC
CURVE_CURVES = TABLES / "check_free_curve_curves.csv"             # k(W) on the grid: learned (per replicate, with sd) and true
CURVE_J1 = TABLES / "check_free_curve_j1_curve.csv"               # jump 1: band coverage and shape test per truth and window
CURVE_J2_SEL = TABLES / "check_free_curve_j2_selection.csv"       # jump 2: truth x term, replicates selected
CURVE_J2_COEF = TABLES / "check_free_curve_j2_coefficients.csv"   # jump 2: fitted / true coefficients
CURVE_MENU = TABLES / "check_free_curve_menu_comparison.csv"      # per dataset: menu-fit winner, free curve vs menu BIC, m4 warning
CURVE_DECISION = TABLES / "check_free_curve_decision.csv"         # the pre-set rule per truth and window, verdicts
CURVE_POSTHOC = TABLES / "check_free_curve_posthoc_upper_range.csv"   # post hoc: curves judged above the second knot
CURVE_FIG_K = FIGURES / "check_free_curve_k_{truth}.png"          # true vs learned k(W) with the band; filled with the truth
CURVE_REPORT = REPORTS / "check_free_curve.md"
CONFIRM_FITS = TABLES / "confirm_free_curve_fits.csv"               # one row per free-curve fit: judged range, knots, k, s, BIC, slope
CONFIRM_CURVES = TABLES / "confirm_free_curve_curves.csv"           # k(W) on the judged grid: learned per replicate, and true
CONFIRM_CURVE_SHAPE = TABLES / "confirm_free_curve_curve_shape.csv" # per dataset kind and truth: band coverage, level, slopes
CONFIRM_J2_SEL = TABLES / "confirm_free_curve_j2_selection.csv"     # jump 2: truth x term, replicates selected
CONFIRM_J2_COEF = TABLES / "confirm_free_curve_j2_coefficients.csv" # jump 2: fitted / true coefficients
CONFIRM_M4 = TABLES / "confirm_free_curve_m4_bic.csv"               # m4 per replicate: BIC(free curve) - BIC(F4)
CONFIRM_DECISION = TABLES / "confirm_free_curve_decision.csv"       # the pre-set rule per truth, verdicts
CONFIRM_FIG_K = FIGURES / "confirm_free_curve_k_{truth}.png"        # true vs learned k(W), band, judged range; filled with the truth
CONFIRM_REPORT = REPORTS / "confirm_free_curve.md"
METHOD_A = TABLES / "check_discovery_method_a_jump2_alone.csv"         # A: per truth, noise, method: draws with exactly the true terms, coefficients
METHOD_A_DRAWS = TABLES / "check_discovery_method_a_draws.csv"         # A: what every method returned on every noise draw
METHOD_B = TABLES / "check_discovery_method_b_step.csv"                # B: per dataset: BIC with and without the step, W_c, height, kept
METHOD_C = TABLES / "check_discovery_method_c_candidates.csv"          # C: per dataset and candidate: held-out log-likelihood, BIC
METHOD_C_WIN = TABLES / "check_discovery_method_c_winners.csv"         # C: per dataset: winner, its family, margin
METHOD_D = TABLES / "check_discovery_method_d_bootstrap_band.csv"      # D: per dataset: coverage of the true k by the bootstrap band, width
METHOD_D_CURVES = TABLES / "check_discovery_method_d_band_curves.csv"  # D: the bands on the grid
METHOD_E = TABLES / "check_discovery_method_e_f4_check.csv"            # E: m4, BIC(free curve) - BIC(F4) with and without the F1-solution start
METHOD_DECISION = TABLES / "check_discovery_method_decision.csv"       # the pre-set rules per part and truth, verdicts
METHOD_FIG_B = FIGURES / "check_discovery_method_b_step_{truth}.png"   # B: curve with and without the step, per replicate
METHOD_FIG_D = FIGURES / "check_discovery_method_d_band_{truth}.png"   # D: bootstrap band vs true k, per dataset
METHOD_REPORT = REPORTS / "check_discovery_method.md"
QUANT_TABLE = TABLES / "check_quantization_step.csv"                # per quantization setting and replicate: step kept, W_c, BIC gain
QUANT_DECISION = TABLES / "check_quantization_decision.csv"         # per setting: steps kept, median W_c; the verdict
QUANT_FIG = FIGURES / "check_quantization_m1_quantization.png"      # m1 fits for none / 0.3 / 0.1 dB side by side
QUANT_REPORT = REPORTS / "check_quantization.md"

# ---- step 4: the drying law on the REAL tails (fit_real_drying.py; TRAIN tails of the frozen set v1 only) ------------
REAL_MENU = TABLES / "fit_real_drying_menu_fits.csv"                # per dataset and candidate: BIC, parameters, s, start spread
REAL_SUMMARY = TABLES / "fit_real_drying_summary.csv"               # per dataset: tails, W_hi, menu winner + margin, slope, jump 2
REAL_CURVES = TABLES / "fit_real_drying_curves.csv"                 # k(W) on the read grid per dataset; band of the main set
REAL_JUMP2 = TABLES / "fit_real_drying_jump2.csv"                   # the formulas of the three jump-2 tools per dataset
REAL_QUALITY = TABLES / "fit_real_drying_fit_quality.csv"           # residuals of the winner and of F1 vs the dry-weather noise
REAL_FAKE = TABLES / "fit_real_drying_fake_comparison.csv"          # real slope and A / A_start vs the ranges of the fake truths
REAL_FIG_K = FIGURES / "fit_real_drying_k_real.png"                 # real k(W) above 0.35 dB, band, F1 line, menu winner
REAL_FIG_SETS = FIGURES / "fit_real_drying_k_sets.png"              # k(W) of the sensitivity and leave-one-event-out sets
REAL_FIG_GALLERY = FIGURES / "fit_real_drying_gallery.png"          # 12 real tails with the fitted winner and the fitted F1
REAL_REPORT = REPORTS / "fit_real_drying.md"

# ---- step 4b: why do some real tails stall or rise? (diagnose_stalled_tails.py; TRAIN tails of the frozen set v1 only) ----
TAIL_CLASSES = PRODUCTS / "tails_model_classes.csv"                  # one row per TRAIN tail: class of the first 60 min + the outside evidence (later steps read it)
STALL_COUNTS = TABLES / "diagnose_stalled_tails_class_counts.csv"    # tails per class: overall, per event, per event flag
STALL_PART_A = TABLES / "diagnose_stalled_tails_a_outside_evidence.csv"   # part A: each indicator, non-decaying vs decaying, odds ratio
STALL_PART_B = TABLES / "diagnose_stalled_tails_b_antenna_tests.csv"      # part B: B1 to B5, value, rule, met or not
STALL_MASTER = TABLES / "diagnose_stalled_tails_b4_master_curve.csv"      # B4: the master curve and how well each group follows it
STALL_RATE = TABLES / "diagnose_stalled_tails_b5_rate_per_level.csv"      # B5: model-free k(A) per subset
STALL_FAKE = TABLES / "diagnose_stalled_tails_c_fake_check.csv"           # part C: menu fits of 6 fake datasets before / after the two fixes
STALL_MENU = TABLES / "diagnose_stalled_tails_c_menu_fits.csv"            # part C: menu fits of the real sets with the two fixes
STALL_SUMMARY = TABLES / "diagnose_stalled_tails_c_summary.csv"           # part C: per set: winner, margins, slope, noise multiplier
STALL_FIG_GALLERY = FIGURES / "diagnose_stalled_tails_gallery.png"        # 4 tails per class with radar and gauge under A
STALL_FIG_MASTER = FIGURES / "diagnose_stalled_tails_b4_collapse.png"     # B4: the tails shifted onto one curve
STALL_FIG_RATE = FIGURES / "diagnose_stalled_tails_b5_rate_per_level.png" # B5: k(A) per subset
STALL_REPORT = REPORTS / "diagnose_stalled_tails.md"

# ---- step 4c: what is in a drying tail? (characterise_tails.py; TRAIN tails of the frozen set v1 only, no law fit) ----
TAIL_STARTS = PRODUCTS / "tails_starts.csv"                          # one row per TRAIN pair: the starts S0, S1, S2 with A_start, end and kept flag per rule (later steps read it)
CHAR_A1_TAILS = TABLES / "characterise_tails_a1_per_tail.csv"        # A1: per tail its rise, its control days, the control percentile
CHAR_A1 = TABLES / "characterise_tails_a1_summary.csv"               # A1: share of tails above their control, per class, night / day, dew risk
CHAR_A2 = TABLES / "characterise_tails_a2_per_minute.csv"            # A2: beta, gamma, alpha per minute with intervals, expected path-rain beta
CHAR_A3 = TABLES / "characterise_tails_a3_start_rules.csv"           # A3: one row per start rule: counts, classes, shape, t_half, k(A) slope, delay
CHAR_A3_RATE = TABLES / "characterise_tails_a3_rate_per_level.csv"   # A3: the model-free k(A) per level and start rule
CHAR_A3_CURVES = TABLES / "characterise_tails_a3_median_curves.csv"  # A3: median A / A_start per minute and start rule
CHAR_FIG_A1 = FIGURES / "characterise_tails_a1_rise_vs_control.png"
CHAR_FIG_A2 = FIGURES / "characterise_tails_a2_length_per_minute.png"
CHAR_FIG_A3 = FIGURES / "characterise_tails_a3_median_curves.png"
CHAR_FIG_GALLERY = FIGURES / "characterise_tails_a3_gallery.png"
CHAR_REPORT = REPORTS / "characterise_tails.md"

# ---- step 4c2: how fast does A fall once the rain nearby has ended? (check_drying_after_rain.py; TRAIN pairs, no law fit) ----
AFTER_RAIN_SETS = TABLES / "check_drying_after_rain_a1_a2_sets.csv"        # A1 + A2: one row per threshold and start (S1, placebo, S0, the same pairs from S0)
AFTER_RAIN_DIFF = TABLES / "check_drying_after_rain_a2_difference.csv"     # A2: S1 - placebo with the event-bootstrap interval, per threshold
AFTER_RAIN_PLACEBO_SPLIT = TABLES / "check_drying_after_rain_a2_placebo_split.csv"   # A2, information: placebo tails by where their start falls (before / after the own S1)
AFTER_RAIN_CURVES = TABLES / "check_drying_after_rain_median_curves.csv"   # median A / A_start per minute, threshold and start
AFTER_RAIN_RATE = TABLES / "check_drying_after_rain_rate_per_level.csv"    # k(A) per level, threshold and start
AFTER_RAIN_SPLIT = TABLES / "check_drying_after_rain_a3_before_after.csv"  # A3: k per level before and after S1, and for the placebo split
AFTER_RAIN_PAIRS = TABLES / "check_drying_after_rain_pairs.csv"            # per pair: the three starts, ends, usable lengths, A_start, time to half
AFTER_RAIN_FIG_CURVES = FIGURES / "check_drying_after_rain_median_curves.png"
AFTER_RAIN_FIG_SPLIT = FIGURES / "check_drying_after_rain_a3_before_after.png"
AFTER_RAIN_FIG_GALLERY = FIGURES / "check_drying_after_rain_gallery.png"
AFTER_RAIN_REPORT = REPORTS / "check_drying_after_rain.md"

# ---- step 7a: very short links: enough data, and the fast first phase? (short_links.py; TRAIN only, no law fit) ----
SHORT_INVENTORY = TABLES / "short_links_p1_inventory.csv"              # part 1: one row per length class: links, rain, pairs, tails
SHORT_BY_BAND = TABLES / "short_links_p1_length_by_band.csv"           # part 1: tails per length class x frequency band, per tail set
SHORT_ENOUGH = TABLES / "short_links_p1_enough.csv"                    # part 1: the "enough for Step 7" rule per length cut
SHORT_TAILS = TABLES / "short_links_p2_per_tail.csv"                   # part 2: per tail: drops, ratios, class, expected path rain
SHORT_BY_LENGTH = TABLES / "short_links_p2_by_length.csv"              # part 2: one row per length class: shape, drops, ratios, classes
SHORT_CURVES = TABLES / "short_links_p2_median_curves.csv"             # part 2: median A / A_start per minute and length class, to 3 h
SHORT_MATCHED = TABLES / "short_links_p2_matched.csv"                  # part 2: SHORT against LONG, every matching, band and quantity
SHORT_VERDICT = TABLES / "short_links_p2_verdict.csv"                  # part 2: the rule, its numbers and the verdict
SHORT_EVENTS = TABLES / "short_links_p3_per_event.csv"                 # part 3: per (event, sub-link) pair of a short link: peaks, times
SHORT_PEAKS = TABLES / "short_links_p3_peak_vs_rain.csv"               # part 3: peak A per bin of the peak radar rain
SHORT_P3_BANDS = TABLES / "short_links_p3_by_band.csv"                 # part 3: timing and drying per frequency band
SHORT_RISE = TABLES / "short_links_p3_rise_before_rain.csv"            # part 3: A rising before the radar shows rain, with its placebo
SHORT_FIG_CURVES = FIGURES / "short_links_median_curves.png"
SHORT_FIG_DROP = FIGURES / "short_links_fast_drop_by_length.png"
SHORT_FIG_GALLERY = FIGURES / "short_links_gallery.png"
SHORT_REPORT = REPORTS / "short_links.md"

# ---- step S1: system identification on short links, the model ladder M0-M4 (system_id_ladder.py, library system_id.py) ----
SID_LINK_EVENTS = PRODUCTS / "system_id_link_events.csv"               # one row per link-event: ids, L, band, window, offset, observed and held-out predicted summaries (S3 / S4 read it)
SID_COUNTS = TABLES / "system_id_ladder_counts.csv"                    # the data: links, sub-links, link-events, events, points, per band
SID_FAKE = TABLES / "system_id_ladder_fake_check.csv"                  # fake check: per truth and replicate the chosen rung and the recovered parameters
SID_LADDER = TABLES / "system_id_ladder_ladder.csv"                    # the ladder: every metric per model with intervals, the paired difference to the best
SID_PARAMS = TABLES / "system_id_ladder_parameters.csv"                # the chosen model and M1: every parameter with its event-bootstrap interval
SID_SENS = TABLES / "system_id_ladder_sensitivity.csv"                 # the sensitivity runs (a)-(e)
SID_CHECKS = TABLES / "system_id_ladder_checks.csv"                    # per-shower peaks and the fast drop by length, observed against predicted
SID_AFTER = TABLES / "system_id_ladder_after_last_shower.csv"          # median (A - c) / peak after the last shower, observed against predicted
SID_FIG_LADDER = FIGURES / "system_id_ladder_rmse.png"
SID_FIG_GALLERY = FIGURES / "system_id_ladder_gallery.png"
SID_FIG_CURVES = FIGURES / "system_id_ladder_F_and_P.png"
SID_FIG_AFTER = FIGURES / "system_id_ladder_after_last_shower.png"
SID_REPORT = REPORTS / "system_id_ladder.md"

# ---- step S2: the WAA law on links up to 3 km, the slow part, the discovery rung M5 (system_id_s2.py, library system_id.py) ----
SID2_LINK_EVENTS = PRODUCTS / "system_id_s2_link_events.csv"           # one row per link-event of S<3: ids, L, band, sets, gauge, held-out summaries of the chosen model and of M3 (S3 / S4 read it)
SID2_COUNTS = TABLES / "system_id_s2_counts.csv"                       # links, link-events, events, points per length set, band and input
SID2_FAKE = TABLES / "system_id_s2_fake_checks.csv"                    # the four fake checks: per truth and replicate what was recovered, and PASS / FAIL
SID2_LADDER = TABLES / "system_id_s2_ladder.csv"                       # the ladder per set and input: every metric per model with intervals (judged from R alone)
SID2_LADDER_DRIFT = TABLES / "system_id_s2_ladder_drift_allowed.csv"   # the same after one linear drift per held-out link-event is allowed (the slow part)
SID2_SLOW = TABLES / "system_id_s2_slow_part.csv"                      # slow water or drift: the numbers of the rule, real data and fake data
SID2_LENGTH = TABLES / "system_id_s2_length_rule.csv"                  # the key WAA numbers on S<1, S<2, S<3 with intervals, and whether S<3 lies inside S<1
SID2_RESID = TABLES / "system_id_s2_residuals_by_length.csv"           # held-out residual of the chosen model during rain per length bin
SID2_GAUGE = TABLES / "system_id_s2_gauge_rule.csv"                    # gauge input against radar input on the same link-events
SID2_PARAMS = TABLES / "system_id_s2_parameters.csv"                   # M3 and the chosen model on every set and input: every parameter with its interval
SID2_FITS = TABLES / "system_id_s2_fits.csv"                           # every fit: lag, loss, spread over the starts, s_b, seconds
SID2_AFTER = TABLES / "system_id_s2_after_last_shower.csv"             # median (A - c) / peak after the last shower, observed against M3, M3D, M4, M5
SID2_NAMING = TABLES / "system_id_s2_m5_naming.csv"                    # the M5 surfaces g and h on their grid, the visits per knot and (if M5 qualifies) M5n
SID2_NAMING_TERMS = TABLES / "system_id_s2_m5_naming_terms.csv"        # E-SINDy: every term per lambda, its bag frequency and coefficient
SID2_FIG_LADDER = FIGURES / "system_id_s2_ladder.png"
SID2_FIG_LENGTH = FIGURES / "system_id_s2_length_rule.png"
SID2_FIG_CURVES = FIGURES / "system_id_s2_F_and_P.png"
SID2_FIG_GALLERY = FIGURES / "system_id_s2_gallery_gauge.png"
SID2_FIG_AFTER = FIGURES / "system_id_s2_after_last_shower.png"
SID2_FIG_M5 = FIGURES / "system_id_s2_m5_surfaces.png"
SID2_REPORT = REPORTS / "system_id_s2.md"

# step S2b (system_id_s2b.py): repair of the S2 model (P / F increments, bounded slow store M4b, longer fits)
SID2B_LINK_EVENTS = PRODUCTS / "system_id_s2b_link_events.csv"         # one row per link-event of S<3: held-out predictions of the chosen model, M3 and M4b (for S3)
SID2B_COUNTS = TABLES / "system_id_s2b_counts.csv"                      # link-events, links, events, points per set and band
SID2B_FAKE = TABLES / "system_id_s2b_fake_checks.csv"                   # fake checks Tc (again) and Te (M4b), per replicate, PASS / FAIL
SID2B_TC = TABLES / "system_id_s2b_fake_tc_knots.csv"                   # Tc: F and P at the knots against the truth (median over replicates)
SID2B_LADDER = TABLES / "system_id_s2b_ladder.csv"                      # the ladder per set and input (held-out RMSE, differences, after-rain metrics, peaks)
SID2B_LENGTH = TABLES / "system_id_s2b_length_rule.csv"                 # key numbers on S<1 / S<2 / S<3 with intervals
SID2B_PARAMS = TABLES / "system_id_s2b_parameters.csv"                  # parameters with 90 % event-bootstrap intervals
SID2B_CONVERGENCE = TABLES / "system_id_s2b_convergence.csv"            # every fit: converged or not, iterations, starts past the S2 caps, folds at the cap
SID2B_COMPARE = TABLES / "system_id_s2b_compare_s2.csv"                 # S2 against S2b: chosen models and key numbers
SID2B_AFTER = TABLES / "system_id_s2b_after_last_shower.csv"            # median (A - c) / peak after the last shower, observed and predicted
SID2B_FIG_LADDER = FIGURES / "system_id_s2b_ladder.png"
SID2B_FIG_CURVES = FIGURES / "system_id_s2b_F_and_P.png"
SID2B_FIG_AFTER = FIGURES / "system_id_s2b_after_last_shower.png"
SID2B_FIG_LENGTH = FIGURES / "system_id_s2b_length_rule.png"
SID2B_REPORT = REPORTS / "system_id_s2b.md"

# step S3 (system_id_s3.py): verification of WAA v2 (frozen S2b fits; VAL read for the first time)
SID3_CACHE = CACHE / "system_id_s3"                                    # the S3 datasets and the S3.4 gauge refits
SID2B_CACHE = CACHE / "system_id_s2b"                                  # the frozen S2b fits (read only)
SID3_STEADY = PRODUCTS / "system_id_s3_steady_link_events.csv"          # S3.1: per link-event and rain class, median (A - c) and R (all variants)
SID3_PAIRS = PRODUCTS / "system_id_s3_direction_pairs.csv"              # S3.3: one row per link-event pair (both sub-links of a link)
SID3_COUNTS = TABLES / "system_id_s3_counts.csv"                       # link-events, links, events per dataset and band
SID3_BIN_COUNTS = TABLES / "system_id_s3_bin_counts.csv"               # S3.1: bins, link-events, links per variant, band, rain class
SID3_LENGTH = TABLES / "system_id_s3_length_fits.csv"                  # S3.1: I, s, q, L0 = I / s with 90 % intervals per variant, band, class
SID3_LENGTH_CMP = TABLES / "system_id_s3_length_vs_models.csv"         # S3.1: a F and P of every frozen fit against I and s, per cell
SID3_VERDICT = TABLES / "system_id_s3_length_verdict.csv"              # S3.1: the verdict per model and fit set (TRAIN, robustness, VAL)
SID3_VAL = TABLES / "system_id_s3_val_errors.csv"                      # S3.2: errors of the frozen fits on VAL, with the TRAIN leave-one-event-out numbers
SID3_VAL_AFTER = TABLES / "system_id_s3_val_after_last_shower.csv"     # S3.2: median (A - c) / peak after the last shower on VAL
SID3_DIRECTIONS = TABLES / "system_id_s3_directions.csv"               # S3.3: Spearman and median ratio between the directions
SID3_GAUGE = TABLES / "system_id_s3_gauge_val.csv"                     # S3.4: the VAL gauge refits next to the TRAIN gauge numbers
SID3_FIG_LENGTH = {b: FIGURES / f"system_id_s3_length_{b.replace(' ', '')}.png" for b in ("28 GHz", "32 GHz", "38 GHz")}   # S3.1 per band
SID3_FIG_SUMMARY = FIGURES / "system_id_s3_I_s_L0.png"                  # S3.1: I(R), s(R), L0(R) with the model curves
SID3_FIG_VAL = FIGURES / "system_id_s3_val_errors.png"                  # S3.2: VAL against TRAIN leave-one-event-out
SID3_FIG_AFTER = FIGURES / "system_id_s3_val_after_last_shower.png"     # S3.2
SID3_FIG_DIRECTIONS = FIGURES / "system_id_s3_directions.png"           # S3.3
SID3_REPORT = REPORTS / "system_id_s3.md"

# step S3b (system_id_s3b.py): the S2 frame fitted on ALL link lengths (P pinned by L); judged on VAL
SID3B_CACHE = CACHE / "system_id_s3b"                                  # the S3b fits (every start and bootstrap draw saved)
SID3B_PARAMS = TABLES / "system_id_s3b_parameters.csv"                 # parameters of every fit with 90 % event-bootstrap intervals
SID3B_CONVERGENCE = TABLES / "system_id_s3b_convergence.csv"           # every fit: converged, iterations, spread over the starts, lag profile
SID3B_LENGTH = TABLES / "system_id_s3b_length_through.csv"             # check 1: I and s of every fit's prediction, per set (TRAIN / VAL), band, class
SID3B_LENGTH_CMP = TABLES / "system_id_s3b_length_vs_observed.csv"     # check 1: model I, s against the observed intervals, and the equilibrium comparison
SID3B_VERDICT = TABLES / "system_id_s3b_length_verdict.csv"            # check 1: the verdict per fit and set
SID3B_VAL = TABLES / "system_id_s3b_val_by_length.csv"                 # check 2: VAL errors per length class, S3b and S2b fits
SID3B_VAL_AFTER = TABLES / "system_id_s3b_val_after_last_shower.csv"   # check 2: median (A - c) / peak after the last shower on VAL
SID3B_RESIDUAL = TABLES / "system_id_s3b_residual_vs_L.csv"            # check 3: median (A - A_pred) during rain per band and L bin on VAL
SID3B_FIG_LENGTH = FIGURES / "system_id_s3b_I_s.png"                    # check 1: I(R), s(R) observed against the four fits, TRAIN and VAL
SID3B_FIG_VAL = FIGURES / "system_id_s3b_val_by_length.png"             # check 2
SID3B_FIG_RESIDUAL = FIGURES / "system_id_s3b_residual_vs_L.png"        # check 3
SID3B_FIG_AFTER = FIGURES / "system_id_s3b_val_after_last_shower.png"   # check 2
SID3B_REPORT = REPORTS / "system_id_s3b.md"

# ---- step S4.1 + S4.2: why do events differ, and the conditions (conditions_s4.py) -------------------------
COND4_CACHE = CACHE / "conditions_s4"                                   # the per-event results of a run (re-reading only)
COND4_EVENTS = TABLES / "conditions_s4_events.csv"                     # S4.1 a: per event dates, split, rain facts, links; S4.2 conditions
COND4_SPLIT_HALF = TABLES / "conditions_s4_split_half.csv"             # S4.1 b: the 200 random 18-event halves of TRAIN and VAL
COND4_FAKE = TABLES / "conditions_s4_fake.csv"                         # S4.1 c: fake check, true and estimated k_e, alpha_e per event
COND4_FACTORS = TABLES / "conditions_s4_factors.csv"                   # S4.1 c: per event and fit k_e, alpha_e with link-bootstrap intervals
COND4_FACTOR_SUMMARY = TABLES / "conditions_s4_factor_summary.csv"     # S4.1 c: TRAIN vs VAL, separability
COND4_CLOSE_GAP = TABLES / "conditions_s4_close_gap.csv"               # S4.1 d: the length test through the model with the event factors
COND4_RADAR_GAUGE = TABLES / "conditions_s4_radar_gauge.csv"           # S4.1 e: radar / gauge per event and its Spearman with k_e
COND4_CONDITIONS = TABLES / "conditions_s4_condition_tests.csv"        # S4.2: Spearman of k_e, alpha_e with each condition (BH)
COND4_RESIDUAL = TABLES / "conditions_s4_residual_check.csv"           # S4.2: frozen M4b-B residual during rain per band, night / day, temperature
COND4_FIG_TIMELINE = FIGURES / "conditions_s4_timeline.png"
COND4_FIG_SPLIT_HALF = FIGURES / "conditions_s4_split_half.png"
COND4_FIG_FAKE = FIGURES / "conditions_s4_fake.png"
COND4_FIG_FACTORS = FIGURES / "conditions_s4_factors.png"
COND4_FIG_CLOSE_GAP = FIGURES / "conditions_s4_close_gap.png"
COND4_FIG_RADAR_GAUGE = FIGURES / "conditions_s4_radar_gauge.png"
COND4_FIG_CONDITIONS = FIGURES / "conditions_s4_conditions.png"
COND4_FIG_RESIDUAL = FIGURES / "conditions_s4_residual_check.png"
COND4_REPORT = REPORTS / "conditions_s4.md"

# ---- step S4.3: is the VAL gap the radar? gauge-rescaled length test, heavy-event factors (conditions_s4c.py) ----
COND4C_CACHE = CACHE / "conditions_s4c"                                 # the results of a run (re-reading only)
COND4C_G = TABLES / "conditions_s4c_g_e.csv"                           # g_e per event (both versions), which events got one
COND4C_CELLS = TABLES / "conditions_s4c_length_cells.csv"              # per set, variant, fit and cell: observed and model I, s, inside
COND4C_VERDICT = TABLES / "conditions_s4c_length_verdict.csv"          # per set, variant, fit: median ratios, both-inside share, rules
COND4C_MODEL_FREE = TABLES / "conditions_s4c_model_free.csv"           # observed I, s of TRAIN and VAL per variant and cell, VAL / TRAIN
COND4C_FAKE = TABLES / "conditions_s4c_heavy_fake.csv"                 # heavy-event fake check, per event and replicate
COND4C_HEAVY = TABLES / "conditions_s4c_heavy_factors.csv"             # heavy events: k_e, alpha_e (S4.1 fits), conditions, g_e (if PASS)
COND4C_HEAVY_TESTS = TABLES / "conditions_s4c_heavy_tests.csv"         # heavy events: Spearman tests (if PASS)
COND4C_FIG_RATIOS = FIGURES / "conditions_s4c_ratios.png"
COND4C_FIG_CURVES = FIGURES / "conditions_s4c_observed_curves.png"
COND4C_FIG_FAKE = FIGURES / "conditions_s4c_heavy_fake.png"
COND4C_FIG_HEAVY = FIGURES / "conditions_s4c_heavy_conditions.png"
COND4C_REPORT = REPORTS / "conditions_s4c.md"

# ---- step S5a: WAA v3 refitted on gauge-corrected rain, the slow leftover (M4c), closed forms of F and P (final_s5a.py) ----
S5A_CACHE = CACHE / "final_s5a"                                         # the fits (every start and every bootstrap draw)
S5A_PARAMS = TABLES / "final_s5a_parameters.csv"                        # per fit and parameter: value, cold-bootstrap interval, S3b value
S5A_CONVERGENCE = TABLES / "final_s5a_convergence.csv"                  # per fit: loss, lag, starts, convergence, minutes
S5A_VAL = TABLES / "final_s5a_val.csv"                                  # VAL errors per input (V0 / V1), fit and length class
S5A_LENGTH = TABLES / "final_s5a_length_test.csv"                       # the length test through the model per input, set, fit and cell
S5A_LENGTH_VERDICT = TABLES / "final_s5a_length_verdict.csv"            # per input, set and fit: cells, both inside, the S3b verdict
S5A_FAKE = TABLES / "final_s5a_m4c_fake.csv"                            # the M4c fake check per replicate
S5A_M4C_RULE = TABLES / "final_s5a_m4c_rule.csv"                        # the three M4c criteria on VAL V1
S5A_KNOTS = TABLES / "final_s5a_knots.csv"                              # F and P of the kept model at the knots, with the bootstrap SE
S5A_FORMS = TABLES / "final_s5a_closed_forms.csv"                       # every closed-form candidate per band and pooled: parameters, AICc
S5A_PYSR = TABLES / "final_s5a_pysr.csv"                                # the PySR Pareto fronts (if PySR ran)
S5A_PLUG = TABLES / "final_s5a_plug_back.csv"                           # the plug-back check on VAL V1
S5A_FIG_CURVES = FIGURES / "final_s5a_F_and_P.png"                      # F and P: S3b, S5a refit (with intervals), closed forms
S5A_FIG_VAL = FIGURES / "final_s5a_val.png"                             # VAL errors by length class, V0 and V1
S5A_FIG_AFTER = FIGURES / "final_s5a_after_rain.png"                    # after the last shower on VAL V1: M4b-B against M4c
S5A_FIG_FAKE = FIGURES / "final_s5a_m4c_fake.png"
S5A_REPORT = REPORTS / "final_s5a.md"

# step S6 (tighten_s6.py, library final_model.py): tighten the equation of M4b-B on gauge-corrected rain
S6_CACHE = CACHE / "tighten_s6"                                         # every fit, bootstrap draw, fake replicate, profile, plug-back
S6_CANDIDATE = PRODUCTS / "tighten_s6_candidate.json"                   # the WAA candidate (S6): form, every number, not identified, assumptions
S6_FUNNEL = TABLES / "tighten_s6_funnel.csv"                            # Part 1: every stage of the data funnel with its count and share
S6_BANDS = TABLES / "tighten_s6_bands.csv"                              # Part 1: sub-links per nominal frequency band, kept or not
S6_LENGTHS = TABLES / "tighten_s6_lengths.csv"                          # Part 1: link-events (TRAIN, VAL) and TEST links per length class
S6_MINUTES = TABLES / "tighten_s6_minutes.csv"                          # Part 1: minutes, V1 rain mm and minutes per rain class (TRAIN, VAL)
S6_KNOTS = TABLES / "tighten_s6_knots.csv"                              # Part 2: F, P at the knots, both intervals, W_30, TRAIN minutes per knot
S6_SHARED = TABLES / "tighten_s6_shared_shape.csv"                      # Part 2: the shared-shape tests (VAL rule, AICc / BIC)
S6_VAL = TABLES / "tighten_s6_val.csv"                                  # Part 2: VAL V1 errors and length test of the fits (used, optimistic)
S6_BOOT = TABLES / "tighten_s6_bootstrap.csv"                           # Part 3: estimate and cold-bootstrap interval of every quantity
S6_BOOT_DRAWS = TABLES / "tighten_s6_bootstrap_draws.csv"               # Part 3: every draw
S6_PROFILE = TABLES / "tighten_s6_profile.csv"                          # Part 3: profile against bootstrap intervals of the headline numbers
S6_PROFILE_POINTS = TABLES / "tighten_s6_profile_points.csv"            # Part 3: every profile step
S6_TIMING = TABLES / "tighten_s6_timing.csv"                            # Part 3: timing sensitivity
S6_FAKE = TABLES / "tighten_s6_fake.csv"                                # Part 4: truth and estimate per quantity and replicate
S6_IDENT = TABLES / "tighten_s6_identified.csv"                         # Part 4: identified or not, per quantity
S6_FORMS = TABLES / "tighten_s6_closed_forms.csv"                       # Part 5: every candidate (physics and PySR), AICc, chosen
S6_PYSR = TABLES / "tighten_s6_pysr.csv"                                # Part 5: the PySR Pareto fronts
S6_PLUG = TABLES / "tighten_s6_plug_back.csv"                           # Part 5: the plug-back tries and the rule
S6_FIG_FUNNEL = FIGURES / "tighten_s6_funnel.png"
S6_FIG_F = FIGURES / "tighten_s6_F.png"                                 # F per band: knots, both intervals, closed form, W_30
S6_FIG_P = FIGURES / "tighten_s6_P.png"                                 # P per band likewise
S6_FIG_SHARED = FIGURES / "tighten_s6_shared_shape.png"                 # the band shapes against one shared shape
S6_FIG_FAKE = FIGURES / "tighten_s6_fake.png"                           # fake recovery, truth against estimate
S6_FIG_INTERVALS = FIGURES / "tighten_s6_intervals.png"                 # bootstrap against profile intervals
S6_REPORT = REPORTS / "tighten_s6.md"

# step S7 (validate_s7.py, library final_model.py): validate the frozen S6 candidate without TEST
S7_CACHE = CACHE / "validate_s7"                                        # the PySR runs (real and shuffled) and the dry-path windows
S7_SUMMARY = TABLES / "validate_s7_summary.csv"                         # every check: its rule and its result (pass / fail / not testable)
S7_ALINK = TABLES / "validate_s7_alink_tests.csv"                       # 7.1: a_link against L and the metadata (Spearman, BH)
S7_LINKS = TABLES / "validate_s7_links.csv"                             # 7.1 / 7.5: one row per link (a_link, metadata, W height)
S7_WHEIGHT = TABLES / "validate_s7_w_height.csv"                        # 7.1: W height a_link F at 2 and 5 mm/h per length class (and band)
S7_RMSE = TABLES / "validate_s7_rmse.csv"                               # 7.2: RMSE per length class, band and rain class, VAL and TRAIN
S7_LENGTH = TABLES / "validate_s7_length_test.csv"                      # 7.2: the VAL length test through the model (cells)
S7_SHAPE = TABLES / "validate_s7_shape_short.csv"                       # 7.2: shared against band shapes on short VAL links
S7_AFTER = TABLES / "validate_s7_after_rain.csv"                        # 7.2: median residual against minutes after the rain
S7_DRY_WINDOWS = TABLES / "validate_s7_dry_windows.csv"                 # 7.3a: every dry-path window (counts, W at 5 / 10 / 20 / 40 min)
S7_DRY = TABLES / "validate_s7_dry_path.csv"                            # 7.3a: measured against model W per length class and band
S7_DECAY = TABLES / "validate_s7_decay.csv"                             # 7.3a: decay time and share left at 40-60 min per length class
S7_PAIRS = TABLES / "validate_s7_pairs.csv"                             # 7.3b: the link pairs and the estimator noise
S7_INTERCEPT = TABLES / "validate_s7_intercept.csv"                     # 7.3c: I(R) from short-only and long-only links
S7_PSIZE = TABLES / "validate_s7_p_size.csv"                            # 7.3e: c_band R against ITU-R P.838 (size reference only)
S7_DEW = TABLES / "validate_s7_dew.csv"                                 # 7.4: the dew-suspect jumps used
S7_DEW_EVENTS = TABLES / "validate_s7_dew_before_rain.csv"              # 7.4: link-events with and without dew before the rain
S7_LEFTOVER = TABLES / "validate_s7_leftovers.csv"                      # 7.5: screening, PySR candidates, leave-one-out skill, null
S7_DRYING = TABLES / "validate_s7_event_drying.csv"                     # 7.5: per-event decay time and weather
S7_RAIN_REF = TABLES / "validate_s7_rain_reference.csv"                 # 7.6: radar / gauge per rate class, two binnings
S7_FIG_ALINK = FIGURES / "validate_s7_alink_vs_L.png"
S7_FIG_WHEIGHT = FIGURES / "validate_s7_w_height.png"
S7_FIG_RMSE = FIGURES / "validate_s7_val_rmse.png"
S7_FIG_DRY = FIGURES / "validate_s7_dry_path.png"
S7_FIG_DECAY = FIGURES / "validate_s7_decay.png"
S7_FIG_PAIRS = FIGURES / "validate_s7_pairs.png"
S7_FIG_INTERCEPT = FIGURES / "validate_s7_intercept.png"
S7_FIG_DEW = FIGURES / "validate_s7_dew.png"
S7_FIG_RAIN_REF = FIGURES / "validate_s7_rain_reference.png"
S7_REPORT = REPORTS / "validate_s7.md"

# step S6b (return_s6b.py, libraries wind_drying.py, validation.py): the one return to S6 (drying and wind; the after-rain misfit)
S6B_CACHE = CACHE / "return_s6b"                                        # every fit, fake replicate, fold and bootstrap draw
S6B_CANDIDATE = PRODUCTS / "return_s6b_candidate.json"                  # the candidate for S8 (with the three boxes)
S6B_DIAG = TABLES / "return_s6b_diagnosis.csv"                          # Part 1: every cause with its numbers and verdict
S6B_DIAG_WINDOWS = TABLES / "return_s6b_windows.csv"                    # Part 1: per dry-path window (errors, lag, weather, preceding peak rain)
S6B_FITS = TABLES / "return_s6b_fits.csv"                               # Part 2: M0r, M1, M2 (, M3): loss, extra numbers, VAL numbers
S6B_FAKE = TABLES / "return_s6b_fake.csv"                               # Part 2: fake recovery of the extra numbers
S6B_RULE = TABLES / "return_s6b_box_a_rule.csv"                         # Part 3: every rule per refinement
S6B_FOLDS = TABLES / "return_s6b_folds.csv"                             # Part 3: leave-one-event-out and the 3 windiest dropped
S6B_BOOT = TABLES / "return_s6b_bootstrap.csv"                          # Part 3: the new candidate's cold bootstrap (if one passes)
S6B_BOXES = TABLES / "return_s6b_boxes.csv"                             # Part 4: every effect seen in S2-S7 with its box
S6B_BUDGET = TABLES / "return_s6b_error_budget.csv"                     # Part 4: the VAL error budget
S6B_FIG_DIAG = FIGURES / "return_s6b_after_rain_causes.png"
S6B_FIG_FAKE = FIGURES / "return_s6b_fake.png"
S6B_FIG_DRY = FIGURES / "return_s6b_dry_path.png"
S6B_FIG_BOXES = FIGURES / "return_s6b_boxes.png"
S6B_FIG_BUDGET = FIGURES / "return_s6b_error_budget.png"
S6B_REPORT = REPORTS / "return_s6b.md"

# step S8 (test_s8.py): freeze, then read TEST once
S8_CACHE = CACHE / "test_s8"                                            # the comparison fits (TRAIN) and the TEST datasets (built once)
S8_FREEZE = PRODUCTS / "test_s8_freeze.json"                            # the hashes, committed before TEST is read
S8_BASELINES = PRODUCTS / "test_s8_baselines.csv"                       # the baseline v2 of every (TEST event, sub-link), TRAIN and TEST links
S8_COMPARE = TABLES / "test_s8_comparison_train_val.csv"                # Part 1: the candidate and B0-B2 on TRAIN and VAL
S8_FUNNEL = TABLES / "test_s8_funnel.csv"                               # Part 3: the TEST funnel (S6 format)
S8_LENGTHS = TABLES / "test_s8_lengths.csv"                             # Part 3: TEST link-events per length class (TRAIN and TEST links)
S8_G_E = TABLES / "test_s8_g_e.csv"                                     # Part 3: g_e of the TEST events
S8_RULES = TABLES / "test_s8_rules.csv"                                 # P1-P5: rule, numbers, verdict
S8_RMSE = TABLES / "test_s8_rmse.csv"                                   # P2 / P4 / P5: RMSE of every model and set with intervals
S8_CELLS = TABLES / "test_s8_length_cells.csv"                          # P1: the TEST length test cells
S8_DRY = TABLES / "test_s8_dry_path.csv"                                # P3: the TEST dry-path windows and errors
S8_INFO = TABLES / "test_s8_information.csv"                            # the information items
S8_BUDGET = TABLES / "test_s8_error_budget.csv"                         # the TEST error budget
S8_FIG_RMSE = FIGURES / "test_s8_rmse.png"
S8_FIG_CELLS = FIGURES / "test_s8_length_cells.png"
S8_FIG_DRY = FIGURES / "test_s8_dry_path.png"
S8_FIG_RATIO = FIGURES / "test_s8_test_over_val.png"
S8_FIG_NEW = FIGURES / "test_s8_new_links.png"
S8_FIG_BUDGET = FIGURES / "test_s8_error_budget.png"
S8_REPORT = REPORTS / "test_s8.md"

# step S8b (test_s8b.py): a post-freeze evaluation of a pre-freeze model (the closed F of S6)
S8B_CACHE = CACHE / "test_s8b"                                          # the TRAIN bootstrap draws and fake replicates of the closed-F model
S8B_TRAIN = TABLES / "test_s8b_closed_f_train.csv"                      # Part 1: estimate, bootstrap interval, identifiability
S8B_COMPARE_F = TABLES / "test_s8b_closed_vs_table_f.csv"               # Part 1: closed F against table F per band and knot
S8B_TEST = TABLES / "test_s8b_test.csv"                                 # Part 2: P1-P5 side by side, the paired difference, dry-path errors
S8B_FIG_F = FIGURES / "test_s8b_F.png"
S8B_FIG_RMSE = FIGURES / "test_s8b_rmse.png"
S8B_REPORT = REPORTS / "test_s8b.md"

# step S9 (formulate_s9.py, library waa_equation.py): formulate Part 1
S9_EQUATION = PRODUCTS / "waa_equation_final.json"                      # the official WAA equation (closed F of S8b), every number
S9_REFERENCE = PRODUCTS / "waa_equation_reference.json"                 # the unit-test reference: 3 TEST link-events, inputs and outputs
S9_EXAMPLE = TABLES / "formulate_s9_example_event.csv"                  # figure: one event on a short and a long link
S9_MEMORY = TABLES / "formulate_s9_memory.csv"                          # figure: the equation's response to 60 min of 5 mm/h
S9_CLOSED_BAND = TABLES / "formulate_s9_closed_f_band.csv"              # figure: the closed F and its 90 % band per R
S9_DRY_CURVES = TABLES / "formulate_s9_dry_curves.csv"                  # figure: median dry-path curves, measured and model
S9_BOXES = TABLES / "formulate_s9_boxes.csv"                            # figure: the three boxes after the S8 decision, with sizes
S9_CHECK = REPORTS / "formulate_s9_check.md"                            # the overview checked against the repo (written by hand)
S9_REPORT = REPORTS / "formulate_s9.md"                                 # what S9 wrote, the corrections, open items
FINAL = MODULE / "FINAL.md"                                             # the one-page pointer to the Part 1 result
