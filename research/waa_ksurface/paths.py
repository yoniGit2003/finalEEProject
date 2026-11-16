"""Where the waa_ksurface branch reads and writes its files: the ONLY place with file locations of this branch.

    outputs/tables/    small tables; the file name starts with the script that writes it
    outputs/figures/   figures, same naming
    outputs/reports/   one report per script (<script>.md)
    research/.cache/waa_ksurface/   large arrays and fits (not in git)
Part 1 files (research/wet_antenna, research/baseline) are READ ONLY for this branch.
"""
from pathlib import Path

from research.wet_antenna import paths as WP          # Part 1 locations (read only)

MODULE = Path(__file__).resolve().parent              # research/waa_ksurface
REPO = MODULE.parents[1]
CONFIG_FILE = MODULE / "config.yaml"                  # the settings of this branch (research/config.yaml is Part 1, read only)

OUTPUTS = MODULE / "outputs"
TABLES = OUTPUTS / "tables"
FIGURES = OUTPUTS / "figures"
REPORTS = OUTPUTS / "reports"
CACHE = REPO / "research" / ".cache" / "waa_ksurface"

# ---- the rebuilt Part 1 TRAIN dataset (data.py) ----------------------------------------------------------------
DATA_TRAIN = CACHE / "data_train_all.npz"             # rain input V0 (radar)
LE_TRAIN = CACHE / "le_train_all.pkl"                 # (link-event table, drops)

# ---- Part 1 inputs (read only) ---------------------------------------------------------------------------------
P1_G_E = WP.COND4C_G                                  # g_e per event (TRAIN / VAL), column g_path
P1_LENGTH_TEST = WP.S5A_LENGTH                        # the stored length test (observed I, s per input, set, cell)
P1_LITERATURE_FITS = WP.S8C_FITS                      # Schleiss a / b, Pastorek a / b: parameters and a_link per link
P1_EQUATION = WP.S9_EQUATION                          # the Part 1 equation, with a_link per link
P1_S8_FREEZE = WP.S8_FREEZE                           # names the cache files of B0 (no WAA), B1 (constant), B2 (no memory)

# ---- A1: audit and reproduction gate (a1_audit.py) -------------------------------------------------------------
A1_FUNNEL = TABLES / "a1_audit_funnel.csv"            # the rebuilt funnel next to the Part 1 funnel
A1_LENGTH_CMP = TABLES / "a1_audit_length_test_vs_part1.csv"   # every cell: I, s rebuilt against stored
A1_COUNTS = TABLES / "a1_audit_counts.csv"            # link-events and rainy minutes per band and length class
A1_SHORT_MINUTES = TABLES / "a1_audit_short_minutes.csv"       # minutes per rain class on links < 1 km
A1_GAUGE = TABLES / "a1_audit_gauge_set.csv"          # the gauge-input set and the g_e use
A1_FITS = TABLES / "a1_audit_stored_fits.csv"         # the stored comparison fits: path, loads
A1_REPORT = REPORTS / "a1_audit.md"

# ---- A3: F(R) and P(R) from the length test (a3_length_test.py, library lt.py) ----------------------------------
A3_CELLS = TABLES / "a3_length_test_cells.csv"        # (a) per band and class: I, s with intervals, counts, eligible
A3_DURATION = TABLES / "a3_length_test_duration.csv"  # (b) I per cell and minimum duration, the paired difference
A3_INPUT = TABLES / "a3_length_test_rain_input.csv"   # (c) I and s: gauge-ratio events against g_e = 1 events
A3_EIV = TABLES / "a3_length_test_eiv.csv"            # (d) the errors-in-variables fake check per cell
A3_EIV_NOISE = TABLES / "a3_length_test_eiv_noise.csv"  # (d) the measured radar-vs-gauge scatter (sigma_1, rho)
A3_KNOTS = TABLES / "a3_length_test_knots.csv"        # (e) the cell values the curves go through (after the monotone step)
A3_CURVES = TABLES / "a3_length_test_curves.csv"      # (e) F and P per band on a dense R grid, with intervals (THE product)
A3_FIG_CURVES = FIGURES / "a3_length_test_F_P.png"
A3_FIG_DURATION = FIGURES / "a3_length_test_duration.png"
A3_FIG_INPUT = FIGURES / "a3_length_test_rain_input.png"
A3_FIG_EIV = FIGURES / "a3_length_test_eiv.png"
A3_REPORT = REPORTS / "a3_length_test.md"

# ---- E1: the radar rain error (e1_rain_error.py; rain data only) ------------------------------------------------
E1_TC = TABLES / "e1_rain_error_triple_collocation.csv"   # (a) error variance per sensor and triplet
E1_PAIRS = TABLES / "e1_rain_error_gauge_pairs.csv"       # (b) every City pair: separation, correlation, semivariance
E1_PAIR_FIT = TABLES / "e1_rain_error_gauge_pair_fit.csv" # (b) the line through the semivariance, extrapolated to 0-2 km
E1_REVERSAL = TABLES / "e1_rain_error_reversal.csv"       # (c) conditional means both ways
E1_SPLIT = TABLES / "e1_rain_error_split.csv"             # (c) the unit-gain split of the log variances
E1_SCALE = TABLES / "e1_rain_error_event_scale.csv"       # (d) per event: radar and gauge mm, log ratio
E1_SHARED = TABLES / "e1_rain_error_shared_in_space.csv"  # (d) covariance of ln(radar / gauge) between sites
E1_SYNTHESIS = TABLES / "e1_rain_error_synthesis.csv"     # (d) THE result: sigma_lo, sigma_hi, scale part, shared part
E1_FIG_PAIRS = FIGURES / "e1_rain_error_gauge_pairs.png"
E1_FIG_REVERSAL = FIGURES / "e1_rain_error_reversal.png"
E1_REPORT = REPORTS / "e1_rain_error.md"

# ---- E3: the bin definitions (e3_bin_rule.py; libraries ltest.py, fakeworld.py) ---------------------------------
E3_SIZES = TABLES / "e3_bin_rule_sizes.csv"               # measured noise, link spread, n_short, the adopted rule (later steps read it)
E3_REAL_CELLS = TABLES / "e3_bin_rule_real_cells.csv"     # the plain cells on real data per rule (counts; used only for sigma_link)
E3_FAKE = TABLES / "e3_bin_rule_fake.csv"                 # the fake test of each rule, per cell
E3_ADOPTED_CELLS = TABLES / "e3_bin_rule_adopted_real_cells.csv"   # real data under the adopted rule: counts per cell, eligible
E3_REPORT = REPORTS / "e3_bin_rule.md"

# ---- E2: the candidates (e2_candidates.py) ----------------------------------------------------------------------
E2_FAKE = TABLES / "e2_candidates_fake.csv"               # every candidate, sigma and cell: truth, estimate, errors
E2_VERDICT = TABLES / "e2_candidates_verdict.csv"         # per candidate and sigma: cells inside, verdict
E2_FIG = FIGURES / "e2_candidates_fake.png"
E2_REPORT = REPORTS / "e2_candidates.md"

# ---- E4: F(R) and P(R) re-measured (e4_remeasure.py) -------------------------------------------------------------
E4_CELLS = TABLES / "e4_remeasure_cells.csv"              # per band and class: value, bootstrap interval, error-range interval
E4_CURVES = TABLES / "e4_remeasure_curves.csv"            # F and P per band on a dense R grid (THE product, after review)
E4_LINKS = TABLES / "e4_remeasure_link_residuals.csv"     # per link: mean residual, number of storms
E4_LINK_TEST = TABLES / "e4_remeasure_link_test.csv"      # the per-link factor criterion
E4_VALIDITY = TABLES / "e4_remeasure_validity.csv"      # Deviation 4: fake-world errors per cell, F-valid and P-valid
E4_VS_A3 = TABLES / "e4_remeasure_vs_a3.csv"              # the new F against the superseded A3a values
E4_FIG = FIGURES / "e4_remeasure_F_P.png"
E4_FIG_LINKS = FIGURES / "e4_remeasure_link_residuals.png"
E4_REPORT = REPORTS / "e4_remeasure.md"

# ---- T1: the conditional rain error (t1_conditional_error.py) ---------------------------------------------------
T1_SUMMARY = TABLES / "t1_conditional_error_summary.csv"       # (a) sd of ln(label / gauge) and its split, steady and all windows
T1_BY_CLASS = TABLES / "t1_conditional_error_by_class.csv"     # (a) per label class
T1_WINDOWS = TABLES / "t1_conditional_error_windows.csv"       # (a) every used steady window
T1_COUNTS = TABLES / "t1_conditional_error_counts.csv"         # (b) steady bins against the added per-bin error
T1_RESULT = TABLES / "t1_conditional_error_result.csv"         # THE result: sigma_c_lo, sigma_c_hi and the parts (T3 reads it)
T1_REPORT = REPORTS / "t1_conditional_error.md"

# ---- T2 / T3: candidate (iv) and the rerun at the conditional error (t3_rerun.py) -------------------------------
T2_REAL_CELLS = TABLES / "t3_rerun_gauge_real_cells.csv"       # candidate (iv) on real data: counts per cell, eligible
T3_FAKE = TABLES / "t3_rerun_fake.csv"                         # every candidate, end and cell: truth, estimate, errors
T3_CELLS = TABLES / "t3_rerun_valid_cells.csv"                 # per candidate and cell: inside at both ends -> valid
T3_VERDICT = TABLES / "t3_rerun_verdict.csv"                   # per candidate: valid cells, of them below 2 mm/h, adopted
T3_FIG = FIGURES / "t3_rerun_fake.png"
T3_REPORT = REPORTS / "t3_rerun.md"

# ---- A4: coverage of the (R, W) plane (a4_coverage.py; library fp.py = the measured curves) ----------------------
A4_DOMAIN = TABLES / "a4_coverage_domain.csv"                  # THE frozen learning domain: per band / pooled, regime, cell: minutes, events, in_domain
A4_DOMAIN_GAUGE = TABLES / "a4_coverage_domain_gauge_set.csv"  # the same for the gauge-input set
A4_SUMMARY = TABLES / "a4_coverage_summary.csv"                # per band and regime: domain cells, R and W range
A4_SHARES = TABLES / "a4_coverage_shares.csv"                  # share of the rainy minutes above the measured range of F and of P
A4_INFO = TABLES / "a4_coverage_data.csv"                      # the learning set, the g_e = 1 events, the gauge-input set: counts
A4_EDGES = TABLES / "a4_coverage_edges.csv"                    # sigma_1, the floors, the R and W edge of the domain (A5 reads it)
A4_FIG = FIGURES / "a4_coverage_map.png"
A4_REPORT = REPORTS / "a4_coverage.md"

# ---- A5: the model (a5_build.py; library ksurf.py) ---------------------------------------------------------------
A5_TESTS = TABLES / "a5_build_unit_tests.csv"                  # the three unit tests: number, tolerance, pass
A5_GRADIENT = TABLES / "a5_build_gradient.csv"                 # test 1: exact gradient and finite difference per parameter
A5_REPORT = REPORTS / "a5_build.md"

# ---- A6: fake-data recovery (a6_recovery.py; libraries fake6.py, ladder.py) --------------------------------------
A6_VERDICTS = TABLES / "a6_recovery_verdicts.csv"              # truth, regime, true rung, rung returned, k within 25 %, verdict
A6_CV = TABLES / "a6_recovery_cv.csv"                          # every truth, rung and lambda: train loss, held-out loss, selected
A6_CELLS = TABLES / "a6_recovery_cells.csv"                    # every truth and domain cell: true k, fitted k, weight
A6_EXTENSION = TABLES / "a6_recovery_extension.csv"            # (b): the extension knots, true / start / fitted
A6_FLAGS = TABLES / "a6_recovery_flags.csv"                    # (f): the tail flag per event
A6_LAG = TABLES / "a6_recovery_lag.csv"                        # (a-lag): the loss per d and the d the rule returns
A6_GATE = TABLES / "a6_recovery_wet_gate.csv"                  # Deviation 6: CV error and k_dry with 1 / k_wet held at half / the truth / twice the truth
A6_DETECT_RUNS = TABLES / "a6_recovery_detectability_runs.csv"  # every replicate of the detectability study
A6_DETECT = TABLES / "a6_recovery_detectability.csv"           # the smallest dependence the ladder returns at the margin
A6_FIG = FIGURES / "a6_recovery_surfaces.png"
A6_REPORT = REPORTS / "a6_recovery.md"
A6_DIAG = TABLES / "a6_diagnostic.csv"                         # the one diagnostic: truth (a) with the ingredients switched
A6_DIAG_REPORT = REPORTS / "a6_diagnostic.md"
A6_DIAG2 = TABLES / "a6_diagnostic2.csv"                       # the diagnostic after the rerun: 1 / k_dry against the held k_wet, extension held / free
A6_DIAG2_KNOTS = TABLES / "a6_diagnostic2_knots.csv"           # the extension knots of the free fits
A6_DIAG2_REPORT = REPORTS / "a6_diagnostic2.md"
A6_DIAG3 = TABLES / "a6_diagnostic3.csv"                       # the diagnostic after truth (c) failed: the error of drying K0 and K1 on parts of the data
A6_DIAG3_CURVE = TABLES / "a6_diagnostic3_curve.csv"           # the fitted K1 curve against the true k_dry(W) at the W knots
A6_DIAG3_REPORT = REPORTS / "a6_diagnostic3.md"
S1_RUNS = TABLES / "s1_separability_runs.csv"                  # S1: every replicate, variant and loss: the held-out gain of drying K1 over K0
S1_TABLE = TABLES / "s1_separability.csv"                      # per loss and variant: (c) gains against (a) gains, ratio, separable
S1_CHECK = TABLES / "s1_separability_extension_check.csv"      # the fake check of the extended variant: 1 / k_dry on truth (a)
S1_LIMIT = TABLES / "s1_separability_limit.csv"                # the detection limit: truth (c) at k ratios 5 and 10
S1_REPORT = REPORTS / "s1_separability.md"
S3_VERDICTS = TABLES / "s3_recovery_verdicts.csv"              # step 3 (Deviation 9): the same tables as A6, for the truths rerun with the adopted drying selection
S3_CV = TABLES / "s3_recovery_cv.csv"
S3_CELLS = TABLES / "s3_recovery_cells.csv"
S3_EXTENSION = TABLES / "s3_recovery_extension.csv"
S3_FLAGS = TABLES / "s3_recovery_flags.csv"
S3_DETECT_RUNS = TABLES / "s3_recovery_detectability_runs.csv"
S3_DETECT = TABLES / "s3_recovery_detectability.csv"
S3_FIG = FIGURES / "s3_recovery_surfaces.png"
S3_REPORT = REPORTS / "s3_recovery.md"
S3_SELECTION_TESTS = TABLES / "s3_selection_test.csv"            # Deviation 10: the toy cases of the drying selection, expected against returned

# ---- step 3b: Deviations 12 (tail flag) and 13 (B3'), s3b_recovery.py --------------------------------------------
S3B_RUNS = TABLES / "s3b_recovery_runs.csv"                    # every run: truth, seed, drying rung, K1 gain, recovery, events flagged, false removals, verdict
S3B_CV = TABLES / "s3b_recovery_cv.csv"                        # every run, rung and lambda: train loss, held-out error
S3B_FLAGS = TABLES / "s3b_recovery_flags.csv"                  # every run and event: the tail flag
S3B_B3 = TABLES / "s3b_recovery_b3prime.csv"                   # B3': the statistic of every (c) and (e) replicate, from the K1 fit and from the selected fit
S3B_DETECT = TABLES / "s3b_recovery_detectability.csv"         # per ratio: replicates run, hits, detected
S3B_REPORT = REPORTS / "s3b_recovery.md"
S3B_LOG = REPORTS / "s3b_terminal.log"                         # the terminal output of the detached run
