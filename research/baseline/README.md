# baseline: the dry-weather baseline of every sub-link, and A(t)

## Purpose

A(t) = TL(t) - A_baseline(t), with TL = TSL - RSL. The baseline is the TL a sub-link would show without rain. Every
later step (WAA, the path law) uses A(t), so the baseline comes first, and it must NOT follow the drying tail after
the rain (otherwise the WAA disappears into it).

## Status

**Method v2 frozen** (2026-10-07, commit `ebf122b`): 24 h median of usable dry time before the event, flat, with
flags. Validated by R1 (slot level) and R1b (event level). The products in `outputs/products/` are what later modules
use. Method v1 ("two-window") stays runnable for comparison.

## How to run

From the repository root (full instructions: [`../../docs/REPRODUCE.md`](../../docs/REPRODUCE.md)); settings:
sections `baseline` and `baseline_validation` of [`../config.yaml`](../config.yaml).

```
python -m research.baseline.run_sample v1          # v1 sample run (check C5 reads it)             ~ 20 s
python -m research.baseline.validation.run_all     # validation: checks C1-C6 (R1) and R1b         ~ 4 min
python -m research.baseline.run_sample             # v2 sample run (needs C3)                      ~ 5 s
python -m research.baseline.make_dew_products      # dew products (needs C4)                       ~ 3 s
python -m research.baseline.make_baseline_table    # the frozen v2 products (needs C1, C3, R1b)    ~ 40 s
```

`python -m research run all` runs them in this order before the wet-antenna steps. Setting the environment variable
`BV_SAMPLE=DEV` before `run_all` gives a quick validation run on 24 sub-links; it writes to
`research/.cache/baseline/dev/`, not to `outputs/`. Using the products from another module:
`baseline.a_of_t(cml_id, event_id, res="1min", pad="6h")` returns A(t) of both sub-links of a link (A(t) itself is
never stored).

## Code

| file | what it does |
|---|---|
| `baseline.py` | the method (v1 and v2 functions, dew helpers, `a_of_t`) |
| `paths.py` | every folder and shared file location of this module |
| `run_sample.py [v1\|v2]` | the sample run (12 links / 24 sub-links x 18 significant train events) |
| `make_dew_products.py`, `make_baseline_table.py` | build `outputs/products/` |
| `validation/` | the scientific checks C1-C6 (R1) and R1b; `common.py` holds the shared pieces |

## Outputs

`outputs/products/` (used by later modules: `baselines.csv`, `baseline_uncertainty.csv`, `dew_risk_hourly.csv`,
`dry_jumps.csv`), `outputs/tables/`, `outputs/figures/`, `outputs/reports/`. Every file starts with the name of the
check or script that writes it (`c1_`, `r1b_`, `sample_...`). Everything is committed. The two validation reports
were written by hand from the tables: [`outputs/reports/R1_validation.md`](outputs/reports/R1_validation.md) and
[`outputs/reports/R1b_validation.md`](outputs/reports/R1b_validation.md).

## Method and decisions

[`METHOD.md`](METHOD.md): the v2 rules and flags, the changes v1 -> v2, the known accuracy, the products with their
columns, the v1 method ("two-window"), and what the sample run and the two validations decided.

## Open issues

- Is method v1 (and check C5, which reads the v1 sample) still needed now that v2 is frozen?
