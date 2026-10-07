# baseline_estm: baseline estimation per sub-link and rain event

A(t) = TL(t) - A_baseline(t), with TL = TSL - RSL. The baseline is the TL a sub-link would show without rain. Every
later step (WAA, the path law) uses A(t), so the baseline comes first, and it must NOT follow the drying tail.

Status: **sample run** (12 train links / 24 sub-links x 18 significant train events). All sub-links come after the
sample has been analysed.

## The method (decided in the research discussion)

| Rule | Value | Why |
|---|---|---|
| Use dry samples only | link_state == dry (radar + gauges, never CML) | the baseline is the no-rain TL |
| ... at least X after the last wet bin | X = 1 h | keeps the start of a drying tail out of the baseline |
| Before / after windows | 6 h before the event start, 6 h after its end | dry TL there is the baseline we can see |
| Window value | mean of the 10 s TL (not median) | the mean sits between the 0.3 dB quantization levels |
| Method A (default) | constant = before | cannot follow anything that happens in the event |
| Method B | line from before to after, if \|after - before\| > 0.5 dB | follows a slow drift over a long event |
| No B if a window is an excursion | window mean > 1 dB from the sub-link's usual dry level | an extreme dry shift must not tilt the line |
| Missing | < 1 h of usable dry samples in the before window | no guessing |

Not decided by the learners (filled in by Claude, open for discussion): an excursion in the BEFORE window still uses
A, but the pair is flagged `excursion_before`.

> **Current method: v2** (24 h median before the event, flat, flags), see [`METHOD.md`](METHOD.md). The text below
> describes v1, the first version, which stays runnable for comparison (`run_sample.py v1`). Dew products for later
> research code: `products/` (built by `make_dew_products.py`).

## Validation (R1, R1b)

`validation/` measures how good the baseline is (fake events in dry weather, the after window on real events, daily
cycle, dry-day jumps, the flagged cases, sensitivity) without changing the method. Settings: section
`baseline_validation` of `research/config.yaml`. Run `python research/baseline_estm/validation/run_all.py` (~2.5 min;
`BV_SAMPLE=DEV` for a quick run on 24 sub-links). Results and the reports: `validation/out/REPORT.md` (R1, v1 and
variants) and `validation/out/REPORT_v2.md` (R1b, v2 at event level).

## Files

```
baseline.py            the method (functions, commented)
run_sample.py          runs the sample, writes everything below; settings: section `baseline` of research/config.yaml
results/baselines.csv  one row per sub-link x event (main run, X = 1 h)
results/sensitivity_X.csv   the same for X = 0.5, 1, 2, 3 h (column X_h)
results/summary.md     the quality numbers (rerun overwrites)
figures/overview.png   the quality numbers as plots
figures/examples/      a few sub-link x event pictures, two per (method, flag) kind
```

Rerun: `python research/baseline_estm/run_sample.py` (about a minute).

## How to read the results

There is no "true" baseline to compare with, so a result is judged by checks that should hold if the baseline is right.

**baselines.csv columns.** `before_db` / `after_db`: window means (dB); `*_hours`: how much usable dry data each
window had; `*_std_db`: spread of the 5 min TL in the window; `usual_level_db`: the sub-link's median dry TL over the
season; `closure_db` = after - before; `method` A / B / missing; `flags` (below); `in_event_dry_*`: the out-of-sample
check (below); `A_max_db`: the largest A(t) in the event (how strong the rain signal was, for context).

**Flags.** `no_before`: no baseline (method missing). `no_after`: A used, nothing to compare the end with.
`excursion_before` / `excursion_after`: a window is more than 1 dB from the usual level. That can be a dry-time
baseline shift OR a drying tail longer than X (both look the same here), so look at the example figure.

**The five checks (overview.png, summary.md):**

1. **Coverage.** Share of pairs with a baseline. "missing" pairs cannot be used later.
2. **Precision.** Std of the dry TL in the before window. Small compared to the 0.3 dB RSL step = the window is calm,
   and its mean is a stable number.
3. **Closure (after - before).** How much the dry level moved during the event. For method A it IS the baseline error
   at the event end. Mostly within +-0.5 dB = a constant baseline is reasonable; the tails are where B or flags apply.
4. **Out-of-sample A.** Some dry bins inside the event (gaps between showers, >= X after rain) were NOT used to make
   the baseline, so A there should be about 0. Two readings when it is not:
   positive = baseline too low OR a drying tail still present (this is WAA, what we want to study next);
   negative = baseline too high (A cannot really be negative for long).
5. **Sensitivity to X.** How much `before_db` changes when X is 0.5, 2 or 3 h instead of 1 h. Near 0 = the choice of
   X does not matter for the baseline.

**Example figures.** Top: TL (blue line), usable dry samples (green dots), before / after windows (green bands), the
baseline (red), usual dry level (dotted), event start / end (dashed); background blue = wet, grey = unknown.
Bottom: A(t) = TL - baseline, which should return to ~0 when the rain and the tail are over.
