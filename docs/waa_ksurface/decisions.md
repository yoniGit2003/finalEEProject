# WAA k-surface branch: decisions log

One entry per step, written **before** the step runs (what, why, dated pass rule) and completed **after** it (result,
verdict). A choice made after seeing a result is marked **DEVIATION** with its reason. The rules are in `plan.md`.

---

## Step 0: set up the branch and find the Part 1 inputs (2026-11-11)

**What.** Create the branch `waa-ksurface` from `main` (`fde431c`), the three branch folders, `plan.md` and this log. List the
Part 1 files the branch will reuse and confirm that each exists and loads. No analysis.

**Why.** The branch must reuse the Part 1 baseline, offset, rain input, selection and split unchanged, so that any difference
in its result comes from the model and not from different data. That needs a fixed list of what is reused, checked before
anything is built on it.

**Pass rule (2026-11-11).** Every file on the list exists and loads (read-only, in the `cmlpipe` conda env); no file outside the
branch folders is changed; nothing is run beyond the load check.

**Result.** Load check run with `C:\Users\yonic\miniconda3\envs\cmlpipe\python.exe` (the base miniconda Python has no `numba`
or `xarray`, so the research libraries do not import there). Every file below loaded, except the fit cache, which is not on
this machine.

| what | path | loads |
|---|---|---|
| baseline v2 (frozen), one row per (event, sub-link) | `research/baseline/outputs/products/baselines.csv` | yes (40 186 rows) |
| baseline uncertainty (`offset_std`, used by the selection) | `research/baseline/outputs/products/baseline_uncertainty.csv` | yes (566 rows) |
| baseline code | `research/baseline/baseline.py` | yes |
| TEST baselines (written in S8 with the frozen code) | `research/wet_antenna/outputs/products/test_s8_baselines.csv` | yes (12 744 rows) |
| link-event builder: selection rules + offset c (`build_dataset`) | `research/wet_antenna/system_id.py` | yes |
| selection and offset settings (`data` blocks of each step) | `research/config.yaml` (section `wet_antenna`) | yes |
| selection funnel of the frozen fit (S6) | `research/wet_antenna/outputs/tables/tighten_s6_funnel.csv` | yes |
| link-events per length class (S6) | `research/wet_antenna/outputs/tables/tighten_s6_lengths.csv` | yes |
| g_e per event, TRAIN / VAL (S4.3, column `g_path`) | `research/wet_antenna/outputs/tables/conditions_s4c_g_e.csv` | yes (51 TRAIN, 18 VAL events) |
| g_e per event, TEST (S8) | `research/wet_antenna/outputs/tables/test_s8_g_e.csv` | yes (15 events) |
| R = R_radar / g_e (`rescaled`) | `research/wet_antenna/event_factors.py` | yes |
| length test library (`steady_table`, `length_fits`) | `research/wet_antenna/length_test.py` | yes |
| length test through the S5a model (example output) | `research/wet_antenna/outputs/tables/final_s5a_length_test.csv` | yes |
| frozen equation (numbers, intervals, `fit_sha256`) | `research/wet_antenna/outputs/products/waa_equation_final.json` | yes |
| frozen equation code (`WAAEquation`) | `research/wet_antenna/waa_equation.py` | yes (runs) |
| reference link-events for the unit test | `research/wet_antenna/outputs/products/waa_equation_reference.json` | yes |
| unit test of the frozen equation | `research/wet_antenna/tests/test_waa_equation.py` | exists (not run) |
| S8 freeze hashes | `research/wet_antenna/outputs/products/test_s8_freeze.json` | yes |
| event table, event split, link split, event-links | `data/work/openmrg/tables/{events,split_events,split_links,event_links}.csv` | yes |
| release (radar along links, link state, 1 min CML, City gauges) | `data/release/openmrg_v1` (via `cmlpipe.api.CMLData`) | yes |
| fit cache (frozen fit `tighten_s6/plug_shared_FP_F_and_P.pkl`, S3 datasets `data_*_all.npz`) | `research/.cache/wet_antenna/` | **no: folder missing** |

Event split: 53 TRAIN, 18 VAL, 18 TEST events. Link split: 285 TRAIN links, 71 TEST links (unseen), 8 excluded.

**g_e coverage (checked 2026-11-11 against the files).**
- TRAIN: `conditions_s4c_g_e.csv` has 51 of the 53 TRAIN events, so the statement "51/53 TRAIN events, 2 events have no usable
  link-event" matches. The two missing events are E20150601_0000 and E20150601_0755. Every one of their 1 132 (event,
  sub-link) pairs in `baselines.csv` has `missing = True`: they are the first events of the record (1 June 2015), with no 24 h
  of dry data before them. Of the 51 events, 14 use the gauge ratio and 37 have g_e = 1 (too little gauge rain).
- VAL: 18 of 18 events; 5 use the gauge ratio, 13 have g_e = 1.
- TEST: `test_s8_g_e.csv` has 15 of the 18 TEST events. The 3 missing events (E20150616_1020, E20150722_1200,
  E20150723_1355, all light) are the "three storms without a ratio" of the overview (`docs/WAA_research_overview.tex`,
  section "The test rules"). They had no gauge links and get the fallback g_e = 1. This matches the statement, but
  "15/18 with 3 fallbacks" undercounts the fallback: of the 15 events in the table, only 4 use the gauge ratio. The
  other 11 also have g_e = 1. One of these, E20150603_1500, has a ratio of 0 (radar sum 0), which was fixed to 1 after
  the TEST read. So 14 of 18 TEST events run on g_e = 1.

**Verdict.** Pass, with one open item: the fit cache is missing (see Step 0b). It is not needed to read the frozen
equation, which is fully in `waa_equation_final.json`.

---

## Step 0b: fixes to the setup (2026-11-11)

**What.**
1. Commit the two pending docs edits on `main`, then rebase `waa-ksurface` onto that commit.
2. Move the branch layout to the repo convention.
3. Record the Python environment in `plan.md`.
4. Fetch the fit cache.
5. Run the Part 1 unit test.
6. Copy the load check into `research/waa_ksurface/check_part1.py`.
7. Check the g_e coverage (above).

**Why.** To keep the branch history clean of unrelated edits, to follow the repo layout and its import boundary, and to
confirm that the frozen equation and its fit can be read before anything is built on them.

**Pass rule (2026-11-11).** `main` gets one commit with only the two docs files. The old folders are gone. The unit test
passes. The cache files load. `check_part1.py` runs from the repo. Only `docs/waa_ksurface/` and `research/waa_ksurface/` are
committed on the branch.

### Deviation 1: branch layout (2026-11-11)

The Step 0 folders `pipeline/waa_ksurface/` and `data_processed/waa_ksurface/` are replaced:

| what | new folder |
|---|---|
| code | `research/waa_ksurface/` |
| small tables and figures (tracked) | `research/waa_ksurface/outputs/` |
| large arrays (git-ignored) | `research/.cache/waa_ksurface/` |
| docs (unchanged) | `docs/waa_ksurface/` |

**Reason:** `pipeline/` is the data layer and must not import `research`, but the branch code reuses `research.wet_antenna`
and `research.baseline`. Also, `data_processed/` was not git-ignored, so its large arrays would have been committed.

**Result.**
- `main`: commit `e2cf9a0` holds only the two docs files; the branch is rebased onto it.
- Layout: moved as above.
- Fit cache: **not fetched**. `get_caches.py` exits with "The fit caches are not published yet". Its `TAG` and `SHA256` are
  empty. Its list of folders also leaves out `tighten_s6` and `system_id_s3`, which hold the files the branch needs. Those
  are `tighten_s6/plug_shared_FP_F_and_P.pkl` and `system_id_s3/data_{train,val}_all.npz`.
- Unit test `research/wet_antenna/tests/test_waa_equation.py`: 1 passed. The reference implementation reproduces the
  stored W and path term of 3 TEST link-events within 1e-6 dB.
- `check_part1.py`: rerun from the repo. Every file loads except the fit cache folder.

**Verdict.** Fail on task 4 (the cache cannot be fetched). Per the stop rule this is reported and not patched. Every other
task passes.

---

# Block A (2026-11-11, unattended)

## Step A0: rules frozen (2026-11-11)

The core idea, the Block A rules A1–A6 and the Block B rules B1–B6 were written to `plan.md` and committed as `f8f076f`
(2026-11-11 12:49 +02:00), before any analysis.

## Step A1: data audit and reproduction gate (2026-11-11)

**Question served.** Q1 (can F and P be trusted?): F and P come from the length test on the Part 1 TRAIN dataset, so the
dataset and the test must first be shown to be the same as in Part 1.

**What.**
- (a) Rebuild the Part 1 TRAIN dataset with the Part 1 code (`system_id.build_dataset`, called exactly as
  `system_id_s3.dataset` calls it: bands from `merged_bands`, all link lengths, `event_split="train"`), saved in
  `research/.cache/waa_ksurface/`.
- (b) Run the Part 1 length test (`length_test.steady_table` + `length_fits`, config `system_id_s3.length_test`) on TRAIN, for
  both rain inputs (V0 = radar, V1 = radar / g_e), and compare I and s with the TRAIN rows of
  `research/wet_antenna/outputs/tables/final_s5a_length_test.csv`.
- (c) Count the TRAIN data. (d) Locate the stored comparison fits. (e) Apply the 32 GHz rule.

**Why.** The fit cache with the Part 1 datasets is missing, so the dataset is rebuilt. A rebuilt dataset that differs from
Part 1 would make every later number incomparable with Part 1 (statistical reason: same sample).

**Pass rule (2026-11-11).** (a) exactly 206 links and 9 780 TRAIN link-events, else STOP. (b) every I and s within 0.01 dB
(dB/km) of the stored value, else STOP. (c), (d): report. (e) 32 GHz is excluded from the branch if it has no eligible cell.

**Note on isolation.** The branch settings go in `research/waa_ksurface/config.yaml`, not in `research/config.yaml`, because
files outside the branch folders must not change. `build_dataset` loads the 1 min store of the kept sub-links for the whole
summer and blanks the minutes inside VAL and TEST events (Part 1 behaviour, unchanged); no VAL or TEST link-event is built.

**Result (2026-11-11).** Script `research/waa_ksurface/a1_audit.py`, report `research/waa_ksurface/outputs/reports/a1_audit.md`.
- (a) **PASS.** 206 links, 9 780 TRAIN link-events (51 events). Every row of the selection funnel equals the Part 1 funnel.
- (b) **PASS.** The stored table holds 15 TRAIN cells (input V1 = radar / g_e). Largest difference: I 0.00005 dB, s 0.00005
  dB/km, interval ends 0.00005 (rounding of the stored table).
- Steady-rain 5 min value (Part 1): a 5 min bin with valid A (3 of 5 minutes) and radar, R in a rain class, the 20 min before
  it all wet (R >= 0.3 mm/h), and the mean R of the 15 min before it within 0.5–2 times the bin's R. Per link-event and class:
  the median of (A − c) and of R.
- Eligible cell (Part 1): a band x rain-class cell with at least 10 link-event values and 3 distinct lengths (so it is
  fitted), and at least 8 distinct links shorter than 1 km.
- (c) Counts: in the report. Gauge-input set: 50 links, 2 425 link-events (18 links and 749 link-events below 1 km). Gauge
  ratio used in 14 events (5 047 link-events); g_e = 1 in 37 events (4 733 link-events).
- (d) Schleiss b, Pastorek b (`literature_s8c_fits.json`) and the Part 1 equation (`waa_equation_final.json`) load, each with
  the a_link of 206 links. **The no-memory fit (B2) and the no-WAA fit (B0) of S8 are missing**: they exist only in the Part 1
  fit cache (`research/.cache/wet_antenna/test_s8/fit_B{2,0}.pkl`), which is not on this machine. Not refitted (rule).
  This does not block Block A; it is an open item for B6.
- (e) 32 GHz has 0 eligible cells (no link below 1 km); 28 GHz and 38 GHz have 5 each. **32 GHz is excluded from this branch
  for every model.** The branch dataset is 187 links and 8 836 link-events (28 and 38 GHz).

**Verdict.** PASS (gates a and b). 32 GHz excluded by rule (e).

**Environment note.** Without the conda env activated, `numpy.linalg.lstsq` crashed natively (Windows exception 0xc06d007f, a
DLL that fails to load). Every script of the branch is therefore run with the env folders (`Library\bin` etc.) first on PATH,
which is what `conda activate cmlpipe` does. No code was changed for this.
