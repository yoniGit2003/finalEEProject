# baseline_estm: the dry-weather baseline of every sub-link, and A(t)

## Purpose

A(t) = TL(t) - A_baseline(t), with TL = TSL - RSL. The baseline is the TL a sub-link would show without rain. Every
later step (WAA, the path law) uses A(t), so the baseline comes first, and it must NOT follow the drying tail after
the rain (otherwise the WAA disappears into it).

## Status

**Method v2 frozen** (2026-10-07, commit `ebf122b`): 24 h median of usable dry time before the event, flat, with
flags. Validated by R1 (slot level) and R1b (event level). The products in `outputs/products/` are what later modules
use. Method v1 ("two-window") stays runnable for comparison.

## How to run (in this order; settings: sections `baseline` and `baseline_validation` of `research/config.yaml`)

```powershell
python research/baseline_estm/run_sample.py v1          # v1 sample run (C5 reads it)                ~ 20 s
python research/baseline_estm/validation/run_all.py     # R1 checks C1-C6 and R1b                     ~ 4 min
python research/baseline_estm/run_sample.py             # v2 sample run (needs C3)                    ~ 5 s
python research/baseline_estm/make_dew_products.py      # dew products (needs C4)                     ~ 3 s
python research/baseline_estm/make_baseline_table.py    # the frozen v2 products (needs C1, C3, R1b)  ~ 40 s
```

`BV_SAMPLE=DEV` before `run_all.py` gives a quick validation run on 24 sub-links; it writes to `.cache/dev/`, not to
`outputs/`. Using the products from another module: `baseline.a_of_t(cml_id, event_id, res="1min", pad="6h")`
returns A(t) of both sub-links of a link (A(t) itself is never stored).

## Code

| file | what it does |
|---|---|
| `baseline.py` | the method (v1 and v2 functions, dew helpers, `a_of_t`) |
| `baseline_estm_paths.py` | every folder and shared file location of this module |
| `run_sample.py [v1\|v2]` | the sample run (12 links / 24 sub-links x 18 significant train events) |
| `make_dew_products.py`, `make_baseline_table.py` | build `outputs/products/` |
| `validation/` | the scientific checks C1-C6 (R1) and R1b; `common.py` holds the shared pieces |

## Outputs

`outputs/products/` (used by later modules), `outputs/tables/`, `outputs/figures/`, `outputs/reports/`. Every file
is described in [`CATALOG.md`](CATALOG.md): columns, units, producer and the command that rebuilds it. The two
validation reports: [`outputs/reports/R1_validation.md`](outputs/reports/R1_validation.md) and
[`outputs/reports/R1b_validation.md`](outputs/reports/R1b_validation.md).

Importance levels (ORGANIZATION.md section 5) are in `outputs/MANIFEST.csv`:
- **key**: the two validation reports (written by hand) and the figures `r1b_offset_by_D.png`, `r1b_c2_delta.png`,
  `c3_daily_cycle.png`, `c4_humidity_dewpoint.png`;
- **product**: the four files in `outputs/products/` and the tables other scripts read (`c1_fake_events.csv.gz`,
  `c3_daily_cycle_per_sublink.csv`, `c4_jumps.csv`, `r1b_fake_event_offsets.csv.gz`,
  `sample_two_window_baselines.csv`);
- **diag**: everything else; moved by /tidy to the local `legacy/` folder (git-ignored; log:
  `outputs/TIDY_LOG.md`). The files /tidy deleted before 2026-10-11 are rebuilt by rerunning the script. The
  decisions they held are in [`DECISIONS.md`](DECISIONS.md).

## Method

[`METHOD.md`](METHOD.md): the v2 rules and flags, the changes v1 -> v2, the known accuracy, the products, and the v1
method ("two-window").

## Open issues

- Is method v1 (and C5, which reads the v1 sample) still needed now that v2 is frozen? See
  `docs/organization/REPORT.md`, open questions.
