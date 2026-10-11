# CATALOG: where every output file is described

Every generated file has exactly one catalog entry: what it contains (columns and units), what question it answers,
which script writes it and the command that rebuilds it. The entries live next to the code that writes the files:

| part | outputs live in | catalog |
|---|---|---|
| data layer (pipeline) | `data/work/openmrg/`, `data/release/openmrg_v1/` (not in git) | [`pipeline/docs/data_dictionary.md`](pipeline/docs/data_dictionary.md) (every table, column and variable) |
| research: baseline | `research/baseline_estm/outputs/` | [`research/baseline_estm/CATALOG.md`](research/baseline_estm/CATALOG.md) |
| research: WAA (drying tails) | `research/WAA_estm/outputs/` | [`research/WAA_estm/CATALOG.md`](research/WAA_estm/CATALOG.md) |
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
| 4 | `python research/baseline_estm/make_dew_products.py` | `outputs/products/dew_risk_hourly.csv`, `dry_jumps.csv` | release, `c4_jumps.csv` (step 2) | ~3 s |
| 5 | `python research/baseline_estm/make_baseline_table.py` | `outputs/products/baselines.csv`, `baseline_uncertainty.csv` | release, step 4, `c1_fake_events.csv.gz`, `c3_*`, `r1b_fake_event_offsets.csv.gz` (step 2) | ~40 s |
| 6 | `python research/WAA_estm/step0_inventory.py` | `WAA_estm/outputs/{figures,tables,reports}/step0_*` | release, step 5 | ~10 s |
| 7 | `python research/WAA_estm/step1_extract_tails.py` | `WAA_estm/outputs/products/tails_1min.parquet`, `tails_index.csv`, `step1_*` | release, steps 4-5 | ~2.5 min |
| 8 | `python research/WAA_estm/step2_select_tails.py` | `WAA_estm/outputs/products/tails_selected.csv`, `step2_*` | step 7 | ~1.5 min |
| 9 | `python research/WAA_estm/step2b_fix_tails.py` | new columns in `WAA_estm/outputs/products/tails_selected.csv`, `step2b_*` | step 8 | ~1 min |
| 10 | `python research/WAA_estm/step2c_freeze.py` | `WAA_estm/outputs/products/tails_model*` (frozen tail set v1), `step2c_*` | step 9 | ~45 s |
| 11 | `python research/WAA_estm/step3a_simulator.py` | `WAA_estm/outputs/products/sim_tails_*`, `sim_noise_params.csv`, `step3a_*` | step 10 | ~3.5 min |
| 12 | `python research/WAA_estm/step3b1_fitting_on_sim.py` | `WAA_estm/outputs/tables/step3b1_*`, `figures/step3b1_*`, `reports/step3b1_fitting_on_sim.md` | step 11 | ~60 min (cached fits: ~2.5 min) |
| 13 | `python research/WAA_estm/step3b1b_pooled_fit.py` | `WAA_estm/outputs/tables/step3b1b_*`, `figures/step3b1b_*`, `reports/step3b1b_pooled_fit.md` | steps 11-12 | ~3 h (cached: ~1 min) |
| 14 | `python research/WAA_estm/step3b1c_stress_test.py` | `WAA_estm/outputs/tables/step3b1c_*`, `figures/step3b1c_win_share.png`, `reports/step3b1c_stress_test.md` | steps 11-12 | ~45 min (cached: ~20 s) |
| - | by hand | `outputs/reports/R1_validation.md`, `R1b_validation.md`, `METHOD.md` | the numbers of steps 1-5 | |

(`research/baseline_estm/` is left out of the paths in steps 1-5.) Re-running a step only rewrites its own files.
