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
