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
SIM_TAILS_M3 = PRODUCTS / "sim_tails_m3.parquet"          #   file, so /tidy can see that calibrate_simulator.py writes it)
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
