# WAA_estm: drying tails and a dynamic equation for wet-antenna attenuation (WAA)

## Purpose

A(t) = TL(t) - baseline = A_rain(t) + δ(t), where δ(t) is the non-rain part (WAA + baseline residual). In the
**drying tail** after the rain on a link has stopped, A_rain = 0, so A(t) = δ(t): what is left there is
mostly WAA. This module isolates WAA from the tails **before** f(L, R̄, Z) is learned. Otherwise a constant WAA looks
like f = L + L0. The later steps learn a dynamic model of the wetting / drying with weak-form SINDy (rain as input) and
compare it with the literature models [1], [9]. Tails are never selected by how well a model fits them (that would be
circular).

## Status

| step | what | state |
|---|---|---|
| 0 | inventory of the inputs (`step0_inventory.py`) | done 2026-10-08 |
| 1 | extract ALL drying tails with their features, descriptive report (`step1_extract_tails.py`) | done 2026-10-08; no selection, no fitting |
| 2 | pre-defined exclusions / truncation, main set, fast-drop test, shape descriptors (`step2_select_tails.py`) | done 2026-10-08; TRAIN only, no fitting |
| 2b | fixes after the step 2 gallery: t0_trace start, strict / lenient ends, earlier rise cut, noise exclusion, level flags (`step2b_fix_tails.py`) | done 2026-10-08; TRAIN only, no fitting |
| 2c | last fixes (radar NaN = unknown, rise cut v2), event diagnostics, **frozen tail set v1** (`step2c_freeze.py`) | done 2026-10-08; frozen |
| 3+ | WAA model (weak-form SINDy), comparison with [1], [9] | not started |

Inputs: release `openmrg_v1` (pinned in `research/config.yaml`) and the frozen baseline v2 products of
`research/baseline_estm/` (`baselines.csv`, `baseline_uncertainty.csv`, `dew_risk_hourly.csv`, `dry_jumps.csv`).

## How to run

Commands from the repository root with `C:\Users\yoni\miniconda3\python.exe` (written `python` below). Settings:
section `waa` of `research/config.yaml`.

```powershell
python research/WAA_estm/step0_inventory.py        # inventory + one example (prints the report)       ~ 10 s
python research/WAA_estm/step1_extract_tails.py    # all tails, products, report (prints the report)   ~ 2.5 min
python research/WAA_estm/step2_select_tails.py     # selection + fast-drop test (needs step 1)          ~ 1.5 min
python research/WAA_estm/step2b_fix_tails.py       # step 2b fixes, adds columns (run after step 2)    ~ 1 min
python research/WAA_estm/step2c_freeze.py          # last fixes + frozen tail set v1 (after 2b)          ~ 45 s
```

Both need the release and the baseline products. Step 1 reads the 1 min CML store once, chunk by chunk (32 sub-links).
`outputs/products/tails_1min.parquet` (~47 MB) and `outputs/products/tails_selected.csv` (~13 MB) are git-ignored
(over the 10 MB limit). Rebuild them with steps 1, 2 and 2b. Step 2 rewrites `tails_selected.csv` without the
step 2b columns, so always run 2b after 2.

## Code

| file | what it does |
|---|---|
| `WAA_estm_paths.py` | every folder and file location of this module |
| `step0_inventory.py` | measures every input (dims, columns, spans, coverage), prints one example event, writes `outputs/reports/step0_inventory.md` |
| `step1_extract_tails.py` | cuts every tail, writes the two products, the A(t) check, the descriptive tables / figures and `outputs/reports/step1_tails.md` |
| `step2_select_tails.py` | applies the step 2 rules, writes `tails_selected.csv`, the fast-drop test, shape descriptors, curves, gallery and `outputs/reports/step2_selection.md` |
| `step2b_fix_tails.py` | applies the step 2b fixes, adds their columns to `tails_selected.csv`, writes the step 2b tables / figures and `outputs/reports/step2b_selection.md` |
| `step2c_freeze.py` | applies the step 2c fixes, event diagnostics, writes the frozen set (`tails_model*`) and `outputs/reports/step2c_freeze.md` |

## Tail definition (decided by the learners; re-stated here, not read from `explorations/`)

One tail per (event, sub-link) pair:

- **Pairs**: every row of `baselines.csv` with `missing == False` and `link_has_wet_slot == True`. These are train + val
  events on train links (test events and links have no baseline). All flag columns are kept, and `split_event` is
  recorded on every tail.
- **A(t)** = 1 min `tl_mean` (release) - `baseline_db` (`baselines.csv`). This is exactly `baseline.a_of_t` (checked on 3
  pairs: `outputs/tables/step1_a_check.csv`).
- **t0** = `last_wet_end`: the first non-wet link_state bin after the link's last wet bin in the network event. The
  ±30 min unknown margin of link_state is **inside** the tail. Short pauses inside the event are not tails.
- **t_first_dry** = the first link_state == dry bin at or after t0 (searched over the whole season).
- **t_end** = the earliest of three stops:
  - the next **wet** link_state bin of this link, searched over the whole season (not only inside the event);
  - the start of a NaN run of A ≥ `gap_min` = 15 min;
  - t0 + `cap_h` = 12 h.

  `ended_by` is next_wet / gap / cap, or season_end when the data stop first. The end never looks at the value of A,
  so a model is not forced to reach 0.
- **No selection** in step 1: no threshold on A, no shape filter, all flags kept as columns. Selection rules are a
  later, separate step.

Inherited from earlier learner decisions:
- the tail start is the pipeline's link_state as is (no extra gauges);
- tails identify only the **drying** (R = 0); the **wetting** needs in-event data where A_rain is mixed in;
- dry-weather excursions (dew?) are WAA without rain, to be used as a check, not as training data.

## Step 2 rules (fixed before the run; settings `waa.selection`)

TRAIN events only for every count, figure and test. VAL tails carry the same columns and are never looked at.

1. **Exclude** (reasons recorded in `exclusion_reason`): t0 is a wet bin; length 0; A_start NaN;
   `use_in_waa_main` False; sub-link status not good; more than 10 % of the tail's minutes without A.
2. **Truncate** the end (`t_use_end`, `truncated_by`) at the earliest of:
   - just before the first `tsl_step` minute;
   - the start of the first sustained dry-weather rise: the 15 min trailing median of A is >= 0.6 dB above its running
     minimum since t_first_dry, on dry minutes, for >= 30 consecutive minutes;
   - otherwise the step 1 `t_end`.
3. **Main set**: kept, A_start >= 1.0 dB, usable length >= 60 min (alternative A_start >= 0.6 dB: counts only).
   **Control**: kept, A_start < 0.5 dB.
4. Nothing is selected by shape. The fast-drop test and the shape descriptors are descriptive.

## Step 2b rules (after the step 2 gallery; settings `waa.selection_2b`)

TRAIN only, as in step 2. These rules replace the step 2 start and ends. Step 2 columns stay in the product for
comparison.

1. **Start t0_trace**: the end of the last 5 min bin with ANY rain trace (radar R_max > 0 on the path, or any gauge
   within 1 km > 0) between t0 and t_first_dry; t0 if the margin has no trace. The radar's presence counts even when
   its value is below the 0.1 mm/h wet threshold. `A_start_trace_db` = mean A over [t0_trace, + 5 min).
2. **Two ends**:
   - lenient = the step 2 ends with the earlier rise cut (3);
   - strict = lenient, or the first rain-trace bin after t0_trace if earlier.

   Usable length is counted from t0_trace.
3. **Dry-rise cut**: same detection as step 2. The cut goes back to the last minute before the rise run whose smoothed
   excess was < 0.15 dB, or to t_first_dry if there is none.
4. **Noise exclusion**: sub-links with `offset_std` > 0.3 dB (one RSL step; `baseline_uncertainty.csv`, measured in
   dry weather, never on tails).
5. **Flags, not exclusions**:
   - `level_shift`: |median A over [t0 + 6 h, t0 + 12 h]| > max(0.3 dB, 2 x offset_std);
   - `undershoot`: the 15 min median of A < -0.3 dB within 3 h after t0_trace.
6. **Main sets**: kept, A_start_trace >= 1.0 dB, usable length >= 60 min, for the strict and the lenient end.
   **Control**: kept, A_start_trace < 0.5 dB.

## Tail set v1 (FROZEN 2026-10-08; the input of the WAA model)

Files: `outputs/products/tails_model.parquet` (1 min rows) + `tails_model_index.csv` (one row per tail); the
alternative set `tails_model_alt.*`. The set id `v1` is in config `waa.freeze.tail_set` and in the `tail_set` column.
File names carry no version, because ORGANIZATION.md forbids `_v1` on active paths. A later set gets a new id and
the old files go to `archive/`. Every rule that defines v1, in order:

1. **Pairs**: `baselines.csv` (frozen baseline v2) rows with `missing == False` and `link_has_wet_slot == True`, i.e.
   train + val events on train links.
2. **A(t)** = 1 min `tl_mean` - `baseline_db`.
3. **t0** = `last_wet_end` (first non-wet link_state bin after the link's last wet bin in the event).
4. **Exclusions**:
   - t0 is a wet bin; length 0; A_start NaN;
   - `use_in_waa_main` False; sub-link status not good; NaN share of the step 1 tail > 10 %;
   - `offset_std` > 0.3 dB (`baseline_uncertainty.csv`);
   - A_start_trace NaN.
5. **Rain trace** (5 min bin): radar R_max > 0 on the path, OR any gauge within 1 km > 0, OR radar R_max NaN on the
   path (radar NaN = unknown).
6. **Start t0_trace** = end of the last rain-trace bin between t0 and t_first_dry (t0 if none). **A_start_trace** =
   mean A over [t0_trace, t0_trace + 5 min).
7. **End** = the earliest of:
   - the step 1 stops: next wet bin of the link (whole season), a NaN run of A ≥ 15 min, t0 + 12 h, season end;
   - just before the first `tsl_step` minute;
   - the step 2b dry-rise cut: 15 min median ≥ 0.6 dB above its running minimum since t_first_dry on dry minutes,
     for ≥ 30 min; cut at the last minute before it with excess < 0.15 dB;
   - the first rain-trace bin after t0_trace (rule 5);
   - **rise cut v2**: from t0_trace on ALL minutes, the 10 min trailing median of A ≥ 0.6 dB above its running
     minimum since t0_trace for ≥ 10 min; cut at the last minute before that run with excess < 0.15 dB (that minute
     is kept), or at t0_trace if there is none.
8. **Main set**: A_start_trace ≥ 1.0 dB and usable length (t0_trace .. end) ≥ 60 min. **Alternative set**:
   A_start_trace ≥ 0.6 dB, same length rule (it contains the main set; column `in_main_set`).
9. **Flags (kept as columns, no exclusion)**:
   - `level_shift`: |median A over [t0 + 6 h, t0 + 12 h]| > max(0.3 dB, 2 x offset_std), step 2b ends;
   - `undershoot`: 15 min median of A < -0.3 dB within 3 h after the step 2b t0_trace;
   - `event_rise_heavy`: more than 30 % of the event's pre-rise tails are cut by rise v2.
10. **Split**: TRAIN for development; VAL is built with the same rules and not looked at. Test events and links have
    no baseline yet, so there are no test tails.

Settings: `waa.tails`, `waa.selection`, `waa.selection_2b`, `waa.freeze` in `research/config.yaml`. Do not change them
for v1.

## Outputs

`outputs/products/` (`tails_1min.parquet`, `tails_index.csv`, `tails_selected.csv`, `tails_model*` = frozen set v1), `outputs/tables/`, `outputs/figures/`,
`outputs/reports/`. Every file is described in [`CATALOG.md`](CATALOG.md); its importance level (key / product / diag /
scratch, ORGANIZATION.md section 5) is in `outputs/MANIFEST.csv`.

- **Key** (kept for the final report): `outputs/reports/step2c_freeze.md` with its figures `step2c_median_curves.png`,
  `step2c_event_E20150617_0705.png`, `step2c_event_E20150701_0030.png`; the fast-drop test
  `outputs/figures/step2_fast_drop_vs_length.png` and `outputs/tables/step2_correlations.csv`.
- **Product**: everything in `outputs/products/`.
- **Diag** (all other step figures, tables and reports): deleted by /tidy (log: `outputs/TIDY_LOG.md`); rerun the
  step's script to get them back. The decisions they held are in [`DECISIONS.md`](DECISIONS.md).

## Method

Step 1 is extraction and description only. The report definitions (return to the dry level, median curves,
dry-weather rises) are in the step 1 report (diag: rerun `step1_extract_tails.py`) and in the `waa.report` settings;
the decisions of every step are in [`DECISIONS.md`](DECISIONS.md). The WAA model method will
go to `METHOD.md` once it is decided.

## Open issues

- **Rise cut v2 vs night-time rises (v1)**: in E20150701_0030 A rises 0.5-1 dB after t0_trace with no radar or gauge
  trace, then decays; v1 cuts these. Missed drizzle (radar 78 km away) or condensation on wet radomes? RH at GbgA is
  missing for those hours (102 June hours missing in total), so this cannot be decided here.
- **Small v1 set**: 429 TRAIN tails from 24 events, but 3 events give 246 of them (E20150804_1830 115,
  E20150826_0520 75, E20150628_0530 56). Resample or weight by EVENT in the model step.

- Step 2b moved the start to t0_trace. 381 of the 959 step 2 main-set tails then fall below A_start = 1 dB: much of
  the step 2 "fast drop" happened while the radar or a gauge still showed a rain trace.
- The strict and lenient main sets hold the same 479 tails (the strict end shortens tails, it does not remove them);
  220 strict tails end at a rain trace.
- `undershoot` mostly marks a lasting level near -0.3 dB right after the rain (baseline offset for that event?), not a
  short dip.
- 104 tails have a **wet** t0 (link still wet when the network event ended; `last_wet_end` is searched only inside
  the event), and 566 tails of the last event were clipped to the season end. All have length 0.
- Dry-weather "rises" at 0.3 dB are frequent (56 % of tails) and partly RSL quantization. A firmer definition is needed
  before they are used.
- Light rain after t_first_dry: the step 2b strict end cuts at the first trace; the lenient end keeps it. Which
  end the WAA model uses (or rain as an input) is still open.
- The L test of the fast drop has little power: the path-rain reference `wet_A_max` itself barely grows with L
  (Spearman 0.13 for >= 25 GHz). See DECISIONS.md, step 2.
- `wet_A_max` spans the link's whole wet phase in the network event (several showers for multi-run events).
