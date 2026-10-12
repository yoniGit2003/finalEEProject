"""Where the baseline module reads and writes its files: the ONLY place with file locations of this module.

    outputs/products/  data that later modules read (baselines.csv, ...)
    outputs/tables/    all other tables; the file name starts with the script that writes it (c1_, r1b_, sample_...)
    outputs/figures/   all figures, same naming
    outputs/reports/   written reports (R1, R1b, sample summary)
Everything in outputs/ is committed. The cache (loaded data, DEV runs, logs) is research/.cache/baseline/: not in
git, safe to delete, rebuilt on the next run.
"""
from pathlib import Path

MODULE = Path(__file__).resolve().parent                   # research/baseline
REPO = MODULE.parents[1]                                   # the repository root
CONFIG_FILE = REPO / "research" / "config.yaml"            # all research settings

OUTPUTS = MODULE / "outputs"
PRODUCTS = OUTPUTS / "products"
TABLES = OUTPUTS / "tables"
FIGURES = OUTPUTS / "figures"
REPORTS = OUTPUTS / "reports"
CACHE = REPO / "research" / ".cache" / "baseline"
DEV = CACHE / "dev"                                        # validation runs with BV_SAMPLE=DEV (throw-away)

# ---- products (make_dew_products.py, make_baseline_table.py) ----------------------------------------------------
BASELINES = PRODUCTS / "baselines.csv"                     # frozen baseline v2, one row per (event, sub-link)
BASELINE_UNCERTAINTY = PRODUCTS / "baseline_uncertainty.csv"
DEW_RISK_HOURLY = PRODUCTS / "dew_risk_hourly.csv"
DRY_JUMPS = PRODUCTS / "dry_jumps.csv"

# ---- sample run (run_sample.py): method v1 = "two-window", method v2 = "pre-event" --------------------------------
SAMPLE_TWO_WINDOW_BASELINES = TABLES / "sample_two_window_baselines.csv"
SAMPLE_TWO_WINDOW_SENSITIVITY = TABLES / "sample_two_window_sensitivity_X.csv"
SAMPLE_TWO_WINDOW_SUMMARY = REPORTS / "sample_two_window_summary.md"
SAMPLE_TWO_WINDOW_OVERVIEW = FIGURES / "sample_two_window_overview.png"
SAMPLE_TWO_WINDOW_EXAMPLES = FIGURES / "sample_two_window_examples"
SAMPLE_PRE_EVENT_BASELINES = TABLES / "sample_pre_event_baselines.csv"

# ---- validation tables that other scripts read (written by validation/, FULL sample) ------------------------------
C1_FAKE_EVENTS = TABLES / "c1_fake_events.csv.gz"
C3_DAILY_CYCLE = TABLES / "c3_daily_cycle_per_sublink.csv"
C4_JUMPS = TABLES / "c4_jumps.csv"
R1B_FAKE_EVENT_OFFSETS = TABLES / "r1b_fake_event_offsets.csv.gz"


def validation_dirs(sample):
    """(tables, figures) folders of a validation run: the real outputs for FULL, the cache for DEV."""
    if sample == "FULL":
        return TABLES, FIGURES
    return DEV / "tables", DEV / "figures"
