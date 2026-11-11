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

## Step A2: theory note (2026-11-11)

**Question served.** Q2 (what is k?): the note fixes what k means, when the rewrite dW/dt = k·(a_link·F − W) is exact, and
what the one-state model cannot represent, before any k is fitted.

**What.** Write `docs/waa_ksurface/theory.md` (no code).

**Why.** The claim "k·(F − W) is a rewrite, not an assumed form" holds only under stated conditions. They must be written
before the fits, so that a later result (for example a k that falls towards the target) is read against a condition fixed
in advance and not explained after the fact.

**Pass rule (2026-11-11).** The note covers every item of rule A2: the one-state assumption; the steady state and the sign
conditions; the factorization with its exact conditions and the degenerate case; the wetting / drying split and the jump of
k; where a_link enters; constant k = the Part 1 fast store, the two-store model not a special case; the list of what is
assumed and what is not.

**Result (2026-11-11).** `docs/waa_ksurface/theory.md` written: 10 sections, every item of rule A2 covered (sections 2, 3, 4,
5, 6, 8, 10). One choice is made explicit there: k is evaluated at the link's own W (dB), not at W / a_link, as rule A5
writes it; the two differ only when k depends on W.

**Verdict.** PASS.

## Step A3: F(R) and P(R) from the length test (2026-11-11)

**Question served.** Q1 (can F and P be trusted?).

**What.** On TRAIN, 28 and 38 GHz, rain input V1 (radar / g_e): (a) the length test per band and rain class with 200
event-bootstrap draws; (b) the intercept against the rain duration; (c) the intercept on gauge-ratio events against g_e = 1
events; (d) an errors-in-variables fake check of the length test; (e) the curves F(R) and P(R) per band.

**Why.** F and P are fixed inputs of every later fit (rule A3f). A biased F would be absorbed by k (a target that is too
high looks like slow wetting), so F must be shown to be steady (b), not driven by the rain reference (c), and not an artefact
of a rain error that depends on the link length (d), before k is learned.

**Choices fixed before the run (2026-11-11), where the rule leaves room.**
1. *Classes (a).* The seven classes 0.1–0.25, 0.25–0.5, 0.5–1, 1–2, 2–5, 5–10, 10–20 mm/h. The Part 1 steady rule is used
   unchanged for all of them (wet = R >= 0.3 mm/h in the 20 min before; previous 15 min mean within 0.5–2 times the bin's R).
   Known consequence, stated now: a bin of the class 0.1–0.25 needs wetter bins before it, so its values come from rain that
   is falling off. A new class enters the curves only if it is eligible by the Part 1 rule (fitted, and >= 8 links below
   1 km).
2. *Duration (b).* The duration of a steady bin = 5 min times the number of consecutive bins ending at it whose R is in the
   same class. "At least 5 / 15 / 30 / 60 min" = at least 1 / 3 / 6 / 12 such bins. Paired interval: the same 200 event
   draws for both durations, 90 % interval of I(30) − I(5). A cell whose 30 min fit is not possible counts in the
   denominator (eligible cells) and not in the numerator. If the rule fires, ONE common duration is used for all cells: the
   longest of 5 / 15 / 30 / 60 min at which every eligible cell is still fitted with >= 8 links below 1 km.
3. *Rain-input diagnostic (c).* Two separate length tests: the link-events of the 14 events with a gauge ratio, and those of
   the 37 events with g_e = 1; 200 event draws inside each group. Report only.
4. *Errors-in-variables check (d).*
   - True rain = the real V1 rain of each link-event (5 min). Known truth, independent of this step's results: the Part 1
     closed forms F_true = A_band·(R / 5)^0.542 (A_band 2.57 / 2.77 dB) and P_true = c_band·R (0.105 / 0.194 dB/km per mm/h).
   - Fake A − c = a_link·F_true(R_true) + L·P_true(R_true) + noise: ln a_link normal with sd 0.19 (the Part 1 spread, 80 % of
     links within 0.77–1.24), independent of L; white noise of sd 0.2 dB per 5 min value; an offset error per link-event with
     the sd of the sub-link's `offset_std`.
   - Observed rain = R_true·exp(e), a classical error. e is AR(1) in time per link-event with the sd σ_1 / √n_pix and the
     lag-1 correlation ρ, n_pix = max(1, L / 2 km) (the radar pixels along the path). σ_1 and ρ are measured on the
     gauge-input set: σ_1 = 1.4826·MAD of ln(radar / gauge) over the 5 min bins with both >= 0.5 mm/h, ρ = the lag-1
     correlation of that log ratio inside link-events. The whole scatter is given to the radar (conservative: part of it is
     the point gauge).
   - The length test runs on the observed rain exactly as on real data. Truth of a cell = F_true and P_true at the cell's
     R_cell (the median observed R), because that is where the cell value is placed on the curve. 20 replicates; the number
     judged is the median over the replicates. PASS if in every eligible cell |I − F_true| <= 0.15 dB and
     |s − P_true| / P_true <= 20 %. Else STOP.
5. *Curves (e).* Per band: the points (R_cell, I) and (R_cell, s) of the eligible cells plus (0, 0). If the cell values do
   not rise with R, they are first made monotone by pooling adjacent violators (weights 1 / bootstrap variance). Monotone
   cubic interpolation (PCHIP). Interval: the same construction on each bootstrap draw, 5–95 % per R. Valid range: from the
   lowest to the highest eligible R_cell. Above it F is held constant and P continues with its last slope; both are marked
   as extrapolation. If (b) fires, the cell values are those of the chosen duration.

**Pass rule (2026-11-11).** (d) as stated, else STOP. (a), (b), (c), (e): report; (b) decides the duration by its rule.

**Result (2026-11-11).** Script `research/waa_ksurface/a3_length_test.py`, report `outputs/reports/a3_length_test.md`.
- (a) All 14 cells (7 classes x 2 bands) are eligible, including the two added classes (28 GHz 0.1–0.25: 12 links below
  1 km; 38 GHz 0.1–0.25: 9). Intercepts rise from 0.57 / 0.59 dB (0.1–0.25 mm/h) to 3.09 / 3.45 dB (10–20 mm/h).
- (b) I(30 min) − I(5 min) excludes zero in 0 of 14 eligible cells. The rule does not fire: F from all steady values.
- (c) Reported in the table `a3_length_test_rain_input.csv`.
- (d) **FAIL.** Measured radar-vs-gauge scatter: σ_1 = 0.708, ρ = 0.403 (14 949 bins). With this rain error the length test
  does not return the truth: the intercept is outside ±0.15 dB in 14 of 14 cells (from +0.19 to +0.36 dB at 0.1–5 mm/h,
  −0.35 dB at 5–10 mm/h, −1.30 / −1.37 dB at 10–20 mm/h), and the slope is outside ±20 % in 5 of 14 cells (largest −60 %,
  28 GHz 0.25–0.5 mm/h). The recovered F is flatter than the true F (too high in light rain, too low in heavy rain).
- (e) The curves were written by the script before the stop (`a3_length_test_curves.csv`). **They are not approved**: by
  rule (d) they must not be used.
- Detail fixed in the code before the run and not written in the entry above: the fake rain error has mean one
  (E[exp(e)] = 1), so the observed rain is unbiased on average.

**Verdict.** STOP at A3d (rule: "PASS if the length test returns the true intercept within 0.15 dB and the true slope within
20 % in every eligible cell. Else STOP").

### The one diagnostic allowed after a failed check (2026-11-11, written before it is run)

**What.** The same fake check with the rain error scaled by 0, 0.25, 0.5 and 1 (σ = scale·σ_1 / √n_pix), 5 replicates each,
everything else unchanged (`a3_eiv_diagnostic.py`).

**Why.** To locate the cause: if the errors vanish at scale 0, the cause is the rain error itself (selecting bins by a noisy
R pulls the true R of a class towards the middle of the rain distribution), and not the noise, a_link, the steady rule or the
Huber line. The scales in between show how small the rain error must be for the test to pass. Nothing is changed or patched
on the basis of this diagnostic; Block A stops after it.

**Result of the diagnostic (2026-11-11).** `outputs/reports/a3_eiv_diagnostic.md`, 5 replicates per scale.

| rain error scale | σ_1 | largest intercept error (dB) | cells with intercept inside ±0.15 dB | largest slope error | cells with slope inside ±20 % |
|---|---|---|---|---|---|
| 0 | 0 | 0.11 | 14 of 14 | 59 % | 13 of 14 |
| 0.25 | 0.18 | 0.16 | 13 of 14 | 34 % | 13 of 14 |
| 0.5 | 0.35 | 0.66 | 10 of 14 | 40 % | 11 of 14 |
| 1 | 0.71 | 1.23 | 0 of 14 | 65 % | 8 of 14 |

- **The cause is the rain error.** Without it the length test returns every intercept within 0.11 dB. The one slope miss at
  scale 0 is 28 GHz 0.1–0.25 mm/h, where the true slope is 0.024 dB/km (an error of 0.014 dB/km).
- The error grows with the size of the rain error and has one shape: too high below about 5 mm/h, too low above. This is the
  flattening expected when bins are sorted into rain classes by a noisy R.
- The check would pass only with a rain error of about a quarter of the measured radar-vs-gauge scatter or less.
- What the diagnostic does not say: how much of σ_1 = 0.71 is radar error. The scatter was all given to the radar (choice 4,
  conservative); part of it is the point gauge against the path mean. No number in this step separates the two.

**Block A stops here (2026-11-11).** A3e curves are not approved; A4, A5 and A6 were not run. Nothing was patched.

**Side observation, not pursued (report only, rule A3c).** The intercepts on the g_e = 1 events are far below those on the
gauge-ratio events (for example 28 GHz 2–5 mm/h: 0.36 dB against 2.04 dB), with few short links in the g_e = 1 cells (0 to 17).
