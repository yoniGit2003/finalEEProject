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
