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
