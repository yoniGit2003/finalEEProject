"""Where the WAA_estm module reads and writes its files: the ONLY place with folder locations of this module.

Layout (rules: ORGANIZATION.md at the repository root; every output file is described in CATALOG.md):
    outputs/products/  deliverables that later steps read (tails_1min.parquet, tails_index.csv)
    outputs/tables/    all other tables; the file name starts with the script that writes it (step0_, step1_)
    outputs/figures/   all figures, same prefix rule
    outputs/reports/   written reports (.md)
    .cache/            git-ignored: logs, scratch (safe to delete)
"""
from pathlib import Path
import sys

MODULE = Path(__file__).resolve().parent                   # research/WAA_estm
REPO = MODULE.parents[1]                                   # the repository root
CONFIG_FILE = REPO / "research" / "config.yaml"            # all research settings (section `waa`)

OUTPUTS = MODULE / "outputs"
PRODUCTS = OUTPUTS / "products"
TABLES = OUTPUTS / "tables"
FIGURES = OUTPUTS / "figures"
REPORTS = OUTPUTS / "reports"
CACHE = MODULE / ".cache"
MANIFEST = OUTPUTS / "MANIFEST.csv"                        # level + purpose of every output file

# every output is written through research/common/outputs.py (save_output); this makes `research.` importable
if str(REPO) not in sys.path:
    sys.path.append(str(REPO))

# ---- step 0: inventory (step0_inventory.py) -----------------------------------------------------------------------
STEP0_EXAMPLE_PNG = FIGURES / "step0_example.png"
STEP0_EXAMPLE_CSV = TABLES / "step0_example.csv"
STEP0_REPORT = REPORTS / "step0_inventory.md"

# ---- step 1: all drying tails (step1_extract_tails.py) ------------------------------------------------------------
TAILS_1MIN = PRODUCTS / "tails_1min.parquet"               # one row per tail minute (large: git-ignored)
TAILS_INDEX = PRODUCTS / "tails_index.csv"                 # one row per tail
STEP1_A_CHECK = TABLES / "step1_a_check.csv"               # A(t) built here vs baseline.a_of_t on 3 pairs
STEP1_RETURN = TABLES / "step1_return_to_dry.csv"          # share of tails back below the threshold after k min
STEP1_MEDIAN_CURVES = TABLES / "step1_median_curves.csv"   # median A per minute and A_start bin
STEP1_RISES = TABLES / "step1_dry_rises.csv"               # one row per tail with a later dry-weather rise
STEP1_FIG_DISTRIBUTIONS = FIGURES / "step1_distributions.png"
STEP1_FIG_MEDIAN_CURVES = FIGURES / "step1_median_curves.png"
STEP1_FIG_RESOLUTION = FIGURES / "step1_resolution_panels.png"
STEP1_REPORT = REPORTS / "step1_tails.md"

# ---- step 2: selection and fast-drop test (step2_select_tails.py) ------------------------------------------------
TAILS_SELECTED = PRODUCTS / "tails_selected.csv"           # index + exclusion, truncation, main set, shape descriptors
STEP2_COUNTS = TABLES / "step2_counts.csv"                 # main-set counts per band x length class x antenna
STEP2_CORRELATIONS = TABLES / "step2_correlations.csv"     # Spearman rho with link-bootstrap 95 % CI
STEP2_TRUNCATION = TABLES / "step2_truncation.csv"         # tail time removed by each truncation rule
STEP2_MEDIAN_CURVES = TABLES / "step2_median_curves.csv"   # median A per minute: main set, control, margin-rain split
STEP2_FIG_LENGTH = FIGURES / "step2_fast_drop_vs_length.png"
STEP2_FIG_MARGIN = FIGURES / "step2_margin_rain.png"
STEP2_FIG_PAIRS = FIGURES / "step2_sublink_pairs.png"
STEP2_FIG_SHAPE = FIGURES / "step2_shape.png"
STEP2_FIG_MEDIAN = FIGURES / "step2_median_curves.png"
STEP2_FIG_GALLERY = FIGURES / "step2_gallery.png"
STEP2_REPORT = REPORTS / "step2_selection.md"

# ---- inputs from the baseline module (frozen v2 products; read only) ----------------------------------------------
BASELINE_MODULE = REPO / "research" / "baseline_estm"

# ---- step 2b: fixes after the step 2 gallery (step2b_fix_tails.py; adds columns to TAILS_SELECTED) -----------------
STEP2B_COUNTS = TABLES / "step2b_counts.csv"               # strict / lenient main-set counts per band x length x antenna
STEP2B_NOISY = TABLES / "step2b_noisy_sublinks.csv"        # sub-links removed by the offset_std rule
STEP2B_TIME = TABLES / "step2b_tail_time.csv"              # tail time kept / lost by each step 2b change
STEP2B_MEDIAN_CURVES = TABLES / "step2b_median_curves.csv" # median A (and A / A_start) aligned at t0_trace
STEP2B_FIG_NOISE = FIGURES / "step2b_offset_std.png"
STEP2B_FIG_MEDIAN = FIGURES / "step2b_median_curves.png"
STEP2B_FIG_UNDERSHOOT = FIGURES / "step2b_undershoot_vs_weather.png"
STEP2B_FIG_EX_SHIFT = FIGURES / "step2b_examples_level_shift.png"
STEP2B_FIG_EX_UNDER = FIGURES / "step2b_examples_undershoot.png"
STEP2B_FIG_GALLERY = FIGURES / "step2b_gallery.png"
STEP2B_REPORT = REPORTS / "step2b_selection.md"
BASELINE_UNCERTAINTY = BASELINE_MODULE / "outputs" / "products" / "baseline_uncertainty.csv"
DEW_RISK_HOURLY = BASELINE_MODULE / "outputs" / "products" / "dew_risk_hourly.csv"

# ---- step 2c: last fixes and the frozen tail set (step2c_freeze.py); the set id (v1) is in config waa.freeze -------
TAILS_MODEL = PRODUCTS / "tails_model.parquet"             # frozen main set, 1 min rows from t0_trace (train + val)
TAILS_MODEL_INDEX = PRODUCTS / "tails_model_index.csv"     # one row per frozen tail
TAILS_MODEL_ALT = PRODUCTS / "tails_model_alt.parquet"     # same rules, A_start_trace >= 0.6 dB
TAILS_MODEL_ALT_INDEX = PRODUCTS / "tails_model_alt_index.csv"
STEP2C_EVENTS = TABLES / "step2c_events.csv"               # one row per TRAIN event of the frozen main set
STEP2C_CHANGES = TABLES / "step2c_changes.csv"             # tails / time changed by each step 2c rule
STEP2C_MEDIAN_CURVES = TABLES / "step2c_median_curves.csv"
STEP2C_FIG_MEDIAN = FIGURES / "step2c_median_curves.png"
STEP2C_FIG_EVENT = FIGURES / "step2c_event_{event_id}.png" # filled with the event id
STEP2C_REPORT = REPORTS / "step2c_freeze.md"

# ---- step 3a (+ 3a-fix): tail simulator (waa_sim.py = library, step3a_simulator.py = calibration + checks) ----------
SIM_TAILS_M1 = PRODUCTS / "sim_tails_m1.parquet"          # simulated tail set per law, 1 min rows (one constant per
SIM_TAILS_M3 = PRODUCTS / "sim_tails_m3.parquet"          #   file, so /tidy can see that step3a_simulator.py writes it)
SIM_TAILS_M4 = PRODUCTS / "sim_tails_m4.parquet"
SIM_TAILS_M5 = PRODUCTS / "sim_tails_m5.parquet"
SIM_INDEX = PRODUCTS / "sim_tails_index.csv"               # per law and tail: time-scale factor, offset, noise, phase
SIM_NOISE = PRODUCTS / "sim_noise_params.csv"              # noise parameters per sub-link (white, slow, tau_noise)
STEP3A_CALIBRATION = TABLES / "step3a_calibration.csv"     # law parameters + median-curve misfit (key)
STEP3A_CHECKS = TABLES / "step3a_shape_checks.csv"         # ratios, t_half, t_e, last-hour std: real vs each law
STEP3A_FIG_MEDIAN = FIGURES / "step3a_real_vs_sim_median.png"   # key
STEP3A_FIG_CROSSING = FIGURES / "step3a_t_half_t_e.png"
STEP3A_FIG_NOISE = FIGURES / "step3a_noise_fit.png"
STEP3A_FIG_GALLERY = FIGURES / "step3a_gallery.png"
STEP3A_THALF = TABLES / "step3a_thalf_by_astart.csv"       # median t_half per A_start bin, real vs each law
STEP3A_FIG_THALF = FIGURES / "step3a_thalf_by_astart.png"
STEP3A_REPORT = REPORTS / "step3a_fix.md"                 # the 3a-fix report (replaces step3a_simulator.md)

# ---- step 3b-1: the fitting method on simulated tails (waa_fit.py = library, step3b1_fitting_on_sim.py) ----------------
STEP3B1_FITS = TABLES / "step3b1_fits.parquet"             # every per-tail fit: dataset, tail, candidate, params, BIC
STEP3B1_CONFUSION = TABLES / "step3b1_confusion.csv"       # true law x chosen candidate, per tail and per dataset
STEP3B1_RECOVERY = TABLES / "step3b1_recovery.csv"         # parameter recovery of the true law's own candidate (key)
STEP3B1_M6 = TABLES / "step3b1_m6_calibration.csv"         # m6 (state switch) calibration scan
STEP3B1_THALF = TABLES / "step3b1_thalf_population.csv"    # t_half per A_start bin: replicate range vs real
STEP3B1_DADT = TABLES / "step3b1_dadt_population.csv"      # dA/dt vs A: replicate range vs real
STEP3B1_ABLATION = TABLES / "step3b1_ablation.csv"         # confusion for m1 / m4 with effects switched off
STEP3B1_FIG_CONFUSION = FIGURES / "step3b1_confusion.png"  # key
STEP3B1_FIG_THALF = FIGURES / "step3b1_thalf_population.png"
STEP3B1_FIG_DADT = FIGURES / "step3b1_dadt_population.png"
STEP3B1_FIG_ABLATION = FIGURES / "step3b1_ablation.png"
STEP3B1_REPORT = REPORTS / "step3b1_fitting_on_sim.md"

# ---- step 3b-1b: pooled fit on event-bootstrap simulated datasets (waa_pooled.py, step3b1b_pooled_fit.py) -------------
STEP3B1B_FITS = TABLES / "step3b1b_fits.csv"               # every pooled fit: truth, replicate, window, candidate, BIC
STEP3B1B_RECOVERY = TABLES / "step3b1b_recovery.csv"       # own-candidate recovery over replicates (key)
STEP3B1B_DECISION = TABLES / "step3b1b_decision.csv"       # the pre-set decision rule per truth and window
STEP3B1B_POPULATION = TABLES / "step3b1b_population.csv"   # t_half per bin and A / A_start at k min: ranges vs real
STEP3B1B_FIG_CONFUSION = FIGURES / "step3b1b_confusion.png"   # key
STEP3B1B_FIG_POPULATION = FIGURES / "step3b1b_population.png"
STEP3B1B_REPORT = REPORTS / "step3b1b_pooled_fit.md"

# ---- step 3b-1c: stress test of the pooled fit (step3b1c_stress_test.py; fake data only) -------------------------------
STEP3B1C_FITS = TABLES / "step3b1c_fits.csv"               # every fit: test, truth, replicate, candidate, start, BIC
STEP3B1C_WIN = TABLES / "step3b1c_win_shares.csv"          # test x truth x candidate: replicates won
STEP3B1C_DECISION = TABLES / "step3b1c_decision.csv"       # the pre-set rule per test and truth (cell verdicts)
STEP3B1C_RECOVERY = TABLES / "step3b1c_recovery.csv"       # own-candidate parameter ratios per test and truth
STEP3B1C_STARTS = TABLES / "step3b1c_start_sensitivity.csv"   # T2: spread of the BIC over the starts
STEP3B1C_RECAL = TABLES / "step3b1c_recalibration.csv"     # T3: old and new truth parameters, A / A_start at 10/30/60
STEP3B1C_NOISE = TABLES / "step3b1c_noise_check.csv"       # T1: real vs simulated dry-noise statistics
STEP3B1C_FIG_WIN = FIGURES / "step3b1c_win_share.png"      # own-candidate win share per test and truth
STEP3B1C_REPORT = REPORTS / "step3b1c_stress_test.md"

# ---- step 3b-2: weak-form SINDy on fake tails (waa_sindy.py = library, step3b2_sindy_fake.py; fake data only) ----------
STEP3B2_FITS = TABLES / "step3b2_fits.csv"                 # every SINDy run: dataset, window, offset, rescale, terms
STEP3B2_SELECTION = TABLES / "step3b2_selection.csv"       # truth x term: replicates in which the term is selected
STEP3B2_COEFS = TABLES / "step3b2_coefficients.csv"        # fitted / true coefficient of the true terms
STEP3B2_DECISION = TABLES / "step3b2_decision.csv"         # the pre-set rule per truth and the verdicts Q1 to Q4
STEP3B2_AGREEMENT = TABLES / "step3b2_agreement.csv"       # per dataset: pooled-fit winner vs SINDy family
STEP3B2_FIG_SELECTION = FIGURES / "step3b2_selection_frequency.png"   # selection frequency of each term per truth
STEP3B2_REPORT = REPORTS / "step3b2_sindy_fake.md"

# ---- step 3b-3: free drying curve, then its equation (waa_curve.py = library, step3b3_curve_then_equation.py) ----------
STEP3B3_FITS = TABLES / "step3b3_fits.csv"                 # one row per free-curve fit (dataset x window): knots, k, s, BIC
STEP3B3_CURVES = TABLES / "step3b3_curves.csv"             # k(W) on the grid: learned (per replicate, with sd) and true
STEP3B3_J1 = TABLES / "step3b3_j1_curve.csv"               # jump 1: band coverage and shape test per truth and window
STEP3B3_J2_SEL = TABLES / "step3b3_j2_selection.csv"       # jump 2: truth x term, replicates selected
STEP3B3_J2_COEF = TABLES / "step3b3_j2_coefficients.csv"   # jump 2: fitted / true coefficients
STEP3B3_MENU = TABLES / "step3b3_menu_comparison.csv"      # per dataset: menu-fit winner, free curve vs menu BIC, m4 warning
STEP3B3_DECISION = TABLES / "step3b3_decision.csv"         # the pre-set rule per truth and window, verdicts
STEP3B3_POSTHOC = TABLES / "step3b3_posthoc_upper_range.csv"   # post hoc: curves judged above the second knot
STEP3B3_FIG_K = FIGURES / "step3b3_k_{truth}.png"          # true vs learned k(W) with the band; filled with the truth
STEP3B3_REPORT = REPORTS / "step3b3_curve_then_equation.md"
STEP3B3B_FITS = TABLES / "step3b3b_fits.csv"               # one row per free-curve fit: judged range, knots, k, s, BIC, slope
STEP3B3B_CURVES = TABLES / "step3b3b_curves.csv"           # k(W) on the judged grid: learned per replicate, and true
STEP3B3B_CURVE_SHAPE = TABLES / "step3b3b_curve_shape.csv" # per dataset kind and truth: band coverage, level, slopes
STEP3B3B_J2_SEL = TABLES / "step3b3b_j2_selection.csv"     # jump 2: truth x term, replicates selected
STEP3B3B_J2_COEF = TABLES / "step3b3b_j2_coefficients.csv" # jump 2: fitted / true coefficients
STEP3B3B_M4 = TABLES / "step3b3b_m4_bic.csv"               # m4 per replicate: BIC(free curve) - BIC(F4)
STEP3B3B_DECISION = TABLES / "step3b3b_decision.csv"       # the pre-set rule per truth, verdicts
STEP3B3B_FIG_K = FIGURES / "step3b3b_k_{truth}.png"        # true vs learned k(W), band, judged range; filled with the truth
STEP3B3B_REPORT = REPORTS / "step3b3b_confirm.md"
STEP3B3C_A = TABLES / "step3b3c_a_jump2_alone.csv"         # A: per truth, noise, method: draws with exactly the true terms, coefficients
STEP3B3C_A_DRAWS = TABLES / "step3b3c_a_draws.csv"         # A: what every method returned on every noise draw
STEP3B3C_B = TABLES / "step3b3c_b_step.csv"                # B: per dataset: BIC with and without the step, W_c, height, kept
STEP3B3C_C = TABLES / "step3b3c_c_candidates.csv"          # C: per dataset and candidate: held-out log-likelihood, BIC
STEP3B3C_C_WIN = TABLES / "step3b3c_c_winners.csv"         # C: per dataset: winner, its family, margin
STEP3B3C_D = TABLES / "step3b3c_d_bootstrap_band.csv"      # D: per dataset: coverage of the true k by the bootstrap band, width
STEP3B3C_D_CURVES = TABLES / "step3b3c_d_band_curves.csv"  # D: the bands on the grid
STEP3B3C_E = TABLES / "step3b3c_e_f4_check.csv"            # E: m4, BIC(free curve) - BIC(F4) with and without the F1-solution start
STEP3B3C_DECISION = TABLES / "step3b3c_decision.csv"       # the pre-set rules per part and truth, verdicts
STEP3B3C_FIG_B = FIGURES / "step3b3c_b_step_{truth}.png"   # B: curve with and without the step, per replicate
STEP3B3C_FIG_D = FIGURES / "step3b3c_d_band_{truth}.png"   # D: bootstrap band vs true k, per dataset
STEP3B3C_REPORT = REPORTS / "step3b3c_method.md"
STEP3B4_TABLE = TABLES / "step3b4_step.csv"                # per quantization setting and replicate: step kept, W_c, BIC gain
STEP3B4_DECISION = TABLES / "step3b4_decision.csv"         # per setting: steps kept, median W_c; the verdict
STEP3B4_FIG = FIGURES / "step3b4_m1_quantization.png"      # m1 fits for none / 0.3 / 0.1 dB side by side
STEP3B4_REPORT = REPORTS / "step3b4_quantization.md"
