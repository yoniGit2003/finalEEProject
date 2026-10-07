# CATALOG: where every output file is described

Every generated file has exactly one catalog entry: what it contains (columns and units), what question it answers,
which script writes it and the command that rebuilds it. The entries live next to the code that writes the files:

| part | outputs live in | catalog |
|---|---|---|
| data layer (pipeline) | `data/work/openmrg/`, `data/release/openmrg_v1/` (not in git) | [`pipeline/docs/data_dictionary.md`](pipeline/docs/data_dictionary.md) (every table, column and variable) |
| research: baseline | `research/baseline_estm/outputs/` | [`research/baseline_estm/CATALOG.md`](research/baseline_estm/CATALOG.md) |
| research: explorations | next to each script (playground) | not cataloged; see [`research/explorations/README.md`](research/explorations/README.md) |
| documents | `docs/` | [`docs/README.md`](docs/README.md) |

A new research module adds its own `CATALOG.md` and a row here.

## How to regenerate everything (in this order)

Commands from the repository root, with `C:\Users\yoni\miniconda3\python.exe` as `python`.

| step | command | writes | needs | time |
|---|---|---|---|---|
| 0 | `python -m cmlpipe run all --workers 4` | `data/work/openmrg/`, `data/release/openmrg_v1/` (up-to-date stages are skipped; frozen tables never change without `--force-split` / `--force`) | `data/raw/openmrg/` | long (first run) |
| 1 | `python research/baseline_estm/run_sample.py v1` | `outputs/tables/sample_two_window_*`, `outputs/figures/sample_two_window_*`, `outputs/reports/sample_two_window_summary.md` | release | ~20 s |
| 2 | `python research/baseline_estm/validation/run_all.py` | `outputs/tables/c1_*` ... `c6_*`, `r1b_*`, `validation_runtimes.csv`, the matching figures | release; C5 needs step 1 | ~4 min |
| 3 | `python research/baseline_estm/run_sample.py` | `outputs/tables/sample_pre_event_baselines.csv` | release, `c3_daily_cycle_per_sublink.csv` (step 2) | ~5 s |
| 4 | `python research/baseline_estm/make_dew_products.py` | `outputs/products/dew_risk_hourly.csv`, `dry_jumps.csv` | release, `c4_jumps.csv` (step 2) | seconds |
| 5 | `python research/baseline_estm/make_baseline_table.py` | `outputs/products/baselines.csv`, `baseline_uncertainty.csv` | release, step 4, `c1_fake_events.csv.gz`, `c3_*`, `r1b_fake_event_offsets.csv.gz` (step 2) | ~1-2 min |
| - | by hand | `outputs/reports/R1_validation.md`, `R1b_validation.md`, `METHOD.md` | the numbers of steps 1-5 | |

(`research/baseline_estm/` is left out of the paths in steps 1-5.) Re-running a step only rewrites its own files.
