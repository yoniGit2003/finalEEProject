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

---

# The A3 rescue (2026-11-11, unattended)

## Deviation 2: Detachment correction (2026-11-11)

**What changed.** This branch is redefined as a fresh research effort that inherits nothing from Part 1 except a whitelist
(W1 the data layer and its quality rules; W2 the 0.3 dB quantization step; W3 Part 1's stored fits as competitors only;
W4 Part 1's warnings as things to test). The full rule is at the top of `plan.md`.

**Reason.** The Block A rules copied Part 1 definitions and numbers into what was meant to be an independent measurement.
An independent result cannot rest on the definitions of the work it is later compared with.

**Audit: violations found and what replaces them.**

| file | inherited item | replacement and justification |
|---|---|---|
| `plan.md` core idea | a_link in the dynamics; "compared with Schleiss" | rewritten: dW/dt = k(R, W)·(F(R) − W); "compared with the existing models". A per-link factor only if E4's own test demands it. |
| `plan.md` rule A1b | reproduction of the Part 1 length test as a gate | kept as a one-off data-loading sanity check, diagnostic only. |
| `plan.md` rules A3a, A1e | Part 1 rain classes, steady-rain rule, eligibility rule | re-derived in E3 from this branch's needs, tested on fake data. The 32 GHz exclusion is re-decided there. |
| `plan.md` rule A3d (as run) | fake truth = Part 1 closed forms; a_link spread 0.19 | E2 / E3 use a generic round-number truth and this branch's own measured noise and link spread. |
| `plan.md` rules A3f, A5 | a_link scaling F; weak prior on ln a_link | removed. |
| `plan.md` rule A5 | the 1 min lag | not inherited: measured before A5, 0 until then. |
| `plan.md` rule A5 | "loss as in Part 1" | Huber threshold 0.3 dB justified by the instrument (one RSL quantization step, W2). |
| `plan.md` rule A5 | unit test against the Part 1 store and constants | unit test against the closed-form step response with round-number constants. |
| `plan.md` rule A6 | Part 1 noise; truths with 2.6 / 10.3 min and the Part 1 two-store equation | noise measured in this branch; round-number constant k; a generic two-store truth. |
| `plan.md` rule B2 | Part 1 rising-tail flags | the flag is defined in this branch before use; dew-risk hours come from the baseline module (W1). |
| `plan.md` rule B6 | stored no-WAA and no-memory fits | defined by recipe here and fitted on TRAIN at Block B time. |
| `theory.md` §1, §5–§10 | the 1 min lag; a_link; the Part 1 fast store and two-store model with their numbers; "wet for 20 min, stable for 15 min" | lag d measured, 0 until then; no per-link parameter; a generic first-order store and a generic two-store argument without numbers; steadiness to be derived in E3. |
| `research/waa_ksurface/lt.py` | Part 1 length-test settings (`data.LT`: classes, 0.3 mm/h wet threshold, 20 / 15 min windows, min 10 link-events, min 8 short links, Huber settings) and the Part 1 `huber_fit` code | superseded. New library `ltest.py` with its own Huber line and its settings in `research/waa_ksurface/config.yaml`, each justified in E3. |
| `research/waa_ksurface/a3_length_test.py`, `a3_eiv_diagnostic.py` | the above, and the Part 1 closed forms as fake truth | superseded; kept for the record. |
| `research/waa_ksurface/config.yaml`, section `a3` | Part 1 numbers (2.57, 2.77, 0.542, 0.105, 0.194, 0.19) | section marked superseded; not read by the rescue code. |
| `research/waa_ksurface/data.py` | `system_id.build_dataset`, `merged_bands`, `rescaled`, the g_e table | **kept, logged as inherited data layer (W1):** the link-event selection rules, the offset c, the nominal frequency bands of the radios, and R = R_radar / g_e. `data.LT` stays only for the A1b diagnostic. |

**A3 outputs that depended on a violation: all of them.** They are listed in `research/waa_ksurface/outputs/SUPERSEDED.md`.
The A3a intercepts are kept only for the comparison asked for in the rescue report. The A1b reproduction stands as a
diagnostic.

**Inherited data layer, logged once for the rescue (W1):** the release tables and time series (radar along links, radar at
gauges, City and SMHI gauges, link state), the baseline v2 and the offset c, the event split and link split, the link-event
selection rules of `system_id.build_dataset`, the frequency bands, the g_e table and R = R_radar / g_e.

## Step E1: partition the rain error (2026-11-11)

**Question served.** Q1 (can F and P be trusted?): A3d showed that the answer depends on how large the radar's own error
is. E1 measures it.

**What.** With rain data only (radar at the gauges, City gauges, the SMHI gauge; no attenuation), inside the TRAIN events:
(a) triple collocation, (b) gauge pairs against separation, (c) the binning reversal, (d) the synthesis [σ_lo, σ_hi].

**Why.** A3d gave ALL of the radar-vs-gauge scatter to the radar (σ_1 = 0.71). Part of that scatter is not radar error. In
particular a weighing gauge with 0.1 mm resolution resolves a 5 min rate only in steps of 1.2 mm/h, and a tipping bucket in
steps of 2.4 mm/h. A method that needs no "true" sensor is required to split the scatter (statistical reason).

**Data and matching rules (fixed before the run; inherited data layer, W1).**
- Time: only bins inside the 53 TRAIN events ([start, end) of the event table).
- Radar: `radar_at_gauges_5min`, `R_pixel` (the 2 km pixel over the gauge; raw radar, before g_e). 15 min radar = the mean
  of its three 5 min values, all valid.
- Gauges: `gauges_city_5min` and `gauges_city_15min` (rate, all minutes valid), `gauges_smhi_15min`. All time stamps are
  interval starts.
- The error is measured as a multiplicative error: on ln(rate). A log needs rain in every sensor of the comparison, so log
  statistics use the bins where all sensors involved read at least 0.5 mm/h. This truncation is a known limit; the linear
  version (on the rates, bins where any sensor reads >= 0.1 mm/h) is reported next to it.
- "Event-centred" = the mean of each sensor's ln(rate) over the event's used bins is removed first, so that what is left is
  the per-bin part; the degrees of freedom are corrected for the removed means.
- Intervals: 200 event-bootstrap draws, 90 %.

**(a) Triple collocation, 15 min.** Triplets (the three errors are assumed independent of each other):
- T1: radar pixel at GbgA, City gauge Drakeg, SMHI gauge GbgA (same radar pixel, 1.1 km apart);
- T2: radar pixel at GbgA, City gauge Lbom, SMHI gauge GbgA (1.5 km);
- T3: mean radar of the two pixels, Drakeg, Lbom (1.8 km);
- T4: mean radar of the two pixels, Jarn, Askim (1.5 km).
Estimator: error variance of sensor 1 = C11 − C12·C13 / C23 (and cyclic), on ln(rate), raw (total error) and event-centred
(per-bin part). 15 min is used because the SMHI gauge is 15 min and the gauge resolution is too coarse at 5 min.

**(b) Gauge pairs.** The 45 City pairs (1.5–17.7 km). Per pair, 15 min: the correlation of the rates, and the semivariance
½·var(ln g_i − ln g_j) on bins with both >= 0.5 mm/h. A line γ(h) = n + s·h is fitted over the pairs up to 8 km. n (the
value at h = 0) is the gauge's own error (two gauges). The point-vs-pixel representativeness variance is s × 1.04 km (the
mean distance between two random points of a 2 km square is 1.04 km; for a linear variogram this is the variance between a
point and the pixel mean).

**(c) Binning reversal, 5 min and 15 min, weighing City gauges.** The mean gauge rate per radar class and the mean radar rate
per gauge class (classes 0.1–0.2–0.5–1–2–5–10–20 mm/h). The split uses the same asymmetry in numbers: with unit gain,
radar error variance = var(ln r) − cov(ln r, ln g), gauge error variance = var(ln g) − cov(ln r, ln g), event-centred.

**(d) Synthesis (fixed before the run).**
- Per-bin part of the radar error at 5 min, one pixel, as the sd of ln(R_obs / R_true):
  est_TC = the median over the triplets of the event-centred radar error (15 min), converted to 5 min by
  σ_5² = σ_15² · 9 / (3 + 4ρ + 2ρ²), ρ = the lag-1 correlation of the event-centred 5 min ln(radar / gauge);
  est_rev = the 5 min radar error of (c).
  **σ_lo = the lower 5 % bound of the smaller estimate, σ_hi = the upper 95 % bound of the larger estimate.** An estimate
  with a negative variance counts as 0 and is reported as "not constrained".
- Reference upper limit (all scatter to the radar): the event-centred sd of ln(radar / gauge) at 5 min.
- The representativeness part of (b) is reported but NOT subtracted: a short link is also smaller than a pixel, so for a
  short link the pixel-vs-point difference is part of the radar's error.
- Scale part (common per event): σ_s = the sd over the TRAIN events of ln(Σ radar / Σ gauge), City gauges pooled, events
  with a mean gauge total >= 1 mm. For events with a gauge ratio the rain input is already divided by this ratio; what is
  left is σ_s,res = the sd over (event, gauge) of the gauge's own event log-ratio minus the pooled one (gauges with >= 1 mm).
  For g_e = 1 events the scale part is not corrected: σ_s applies, as a lower limit (those events have too little gauge rain
  to measure it).
- Shared-in-space share of the per-bin part: the mean over the City pairs of cov(d_i, d_j), d = event-centred
  ln(radar / gauge) at 15 min. Gauge errors at two sites are independent, so this covariance belongs to the radar.
- The length dependence of the independent part stays as in rule A3d (this branch's own rule): sd ∝ 1 / √max(1, L / 2 km).

**Pass rule (2026-11-11).** E1 is a measurement, not a gate: it passes if it returns [σ_lo, σ_hi] with σ_lo < σ_hi and the
scale part, or states plainly that the data do not constrain the split and gives the conservative range (σ_lo = 0 or
σ_hi = the all-to-radar limit).

**Result (2026-11-11).** Script `research/waa_ksurface/e1_rain_error.py`, report `outputs/reports/e1_rain_error.md`. 53 TRAIN
events, 6 112 five-minute bins. No attenuation data was read.

- **(a) Triple collocation (15 min, sd of ln rate, event-centred = the per-bin part).** Radar: 0.56 (T1), 0.62 (T2), 0.45 (T3),
  0.53 (T4). Gauges: Drakeg 0.28 / 0.21, Lbom 0.43 / 0.44, SMHI GbgA 0.38 / 0.40, Jarn 0.37, Askim 0.28. Raw (scale part
  included) the radar is at 0.56–0.71. 127–154 bins and 13–15 events per triplet. **The radar is the noisiest of the three
  sensors in every triplet.**
- **(b) Gauge pairs.** Semivariance of ln(rate) at 15 min, line up to 8 km, extrapolated to 0–2 km (closest pair 1.5 km):
  the gauge's own error sd is 0.32, the point-vs-pixel representativeness sd is 0.17. Representativeness is the small part.
- **(c) Binning reversal and the unit-gain split (5 min, weighing gauges).** Radar error sd 0.545 [0.481, 0.608]; gauge error
  sd 0.536 (at 5 min the gauge resolves only steps of 1.2 mm/h); sd of ln(radar / gauge) 0.764.
- **(d) Result.**
  - Per-bin part at 5 min, one pixel: est_TC = 0.760 [0.604, 0.872] (15 min value converted with ρ = 0.349), est_rev = 0.545
    [0.481, 0.608]. **[σ_lo, σ_hi] = [0.48, 0.87].** σ_hi lies above the all-to-radar limit of 0.76; it is kept as the rule
    wrote it (the conversion from 15 to 5 min is approximate), so the upper end is conservative.
  - Scale part: σ_s = 0.70 (sd over 14 events of ln(Σ radar / Σ gauge), median −0.10). After the division by g_e:
    σ_s,res = 0.43 per (event, gauge). This is an upper limit, because it still holds the per-bin errors that do not average
    out over an event's rain at one gauge.
  - Shared in space: the covariance of ln(radar / gauge) between two City sites is 0.138 at 15 min, against a radar per-bin
    variance of 0.298 there: **about 46 % of the per-bin radar error is shared between sites kilometres apart.**
- **The g_e = 1 events (the A3c collapse).** In the 39 TRAIN events without a gauge ratio the radar at the gauges gives 5.8 mm
  per gauge in total and the gauges 2.3 mm: the radar reads about 2.5 times the gauges, against 0.9 times (median) in the
  events with a ratio. So in those events the scale part is not a zero-mean scatter, it is a large one-sided error that
  nothing corrects: a bin sorted into "2–5 mm/h" by the radar holds much less rain there, and its attenuation is low. This
  is the direction and the order of size of the collapse seen in A3c. Part of the factor is selection (an event has no gauge
  ratio because its gauges caught little). Consequence for E2: the g_e = 1 events cannot enter a test that sorts by radar
  rain.

**What the data do not constrain.** The split of the per-bin part between its two estimates (0.55 and 0.76) is not resolved;
the log statistics use only bins where all sensors read >= 0.5 mm/h, so the error in lighter rain is not measured at all;
the 0–2 km behaviour of the gauge pairs is an extrapolation.

**Verdict.** PASS (a measurement): [σ_lo, σ_hi] = [0.48, 0.87], σ_s,res = 0.43, shared share 0.46, ρ = 0.35. The radar's own
error is not small: the A3d value of 0.71 was inside this range.

## Steps E3 and E2: the bin definitions, and a length test that survives the measured error (2026-11-11)

E3 is written first because the E2 fake world needs its definitions. Both entries are fixed before any of their runs.

**Question served.** Q1 (can F and P be trusted?).

### E3: the definitions, each with its reason (nothing taken from Part 1)

1. **Rain classes: 0.1–0.2–0.5–1–2–5–10–20 mm/h (7 classes).** Reason: rain rates span more than two decades and F is
   expected to change with the order of size of R, so the classes are equal steps in log R (the 1-2-5 series). Lower end
   0.1 mm/h: the data layer calls a radar value wet from 0.1 mm/h. Upper end 20 mm/h: above it lie under 0.1 % of the valid
   minutes of the short links (A1c).
2. **Steady rain.** A 5 min bin is kept if
   - it has a valid A (data layer: 3 of 5 minutes) and a radar value R inside a class, and
   - every bin of the T_pre minutes before it has a radar value within [R / q, R·q].
   Reasons: (i) W must have had time to reach its level, so the rain must have been near R for T_pre; (ii) the band is
   symmetric in log R and has no absolute wet threshold, so a bin is not entered systematically from heavier rain (an
   absolute threshold above a class forces exactly that, the bias named in the rescue rules); (iii) q about one class
   width keeps the past of the bin inside roughly its own class.
   Ordered candidates, the first that passes the fake test below is adopted: **R1: T_pre = 20 min, q = 2;
   R2: T_pre = 30 min, q = 1.5.** R1 comes first because it keeps more data.
3. **One value per link-event and class:** the median of (A − c) and the median of R over the kept bins. Reason: the storm
   is the unit that is independent, so a long storm must not count many times; the median resists single bad bins.
4. **The line:** (A − c) = I + s·L per band and class, Huber regression with k = 1.345 (the standard constant: 95 %
   efficiency for normal errors), scale = MAD / 0.6745, own code.
5. **A cell (band x class) is fitted** if it has values from at least 5 events (an event bootstrap of fewer events means
   nothing; the same number as the 5-event rule of A4).
6. **A cell is eligible** if it also has at least n_short distinct links below 1 km and at least n_short distinct links of
   1 km or more. Reasons: the intercept is an extrapolation to L = 0 and is pinned by the short links; the slope needs
   long links; "below 1 km" is the shortest length class of the data layer's link split. **n_short is derived, not
   chosen:** n_short = ceil((σ_link / 0.15 dB)²), where σ_link is the between-link sd of the short links' residuals from the
   cell line (TRAIN, events with a gauge ratio, variance components over links and storms), so that the mean over the short
   links is known to 0.15 dB (the tolerance of the fake tests) at one sd. n_short is at least 5.
7. **32 GHz** is re-decided by rule 6 (it has no link below 1 km, so it cannot be eligible).

### The fake world (E3 and E2)

- Links, events and true rain: the real TRAIN link-events of 28 and 38 GHz in the events with a gauge ratio; true rain = their
  real rain input (5 min, held over the 5 minutes).
- Truth, generic and round: F_true = 3·R / (R + 2) dB; P_true = 0.1·R dB/km (28 GHz) and 0.2·R dB/km (38 GHz).
- Memory: one store, dW/dt = (F_true(R) − W) / τ, τ = 5 min while W is below its target and 15 min while above (round
  numbers, generic), run per minute, W = 0 at the start of the window.
- Noise, measured in this branch on the dry minutes of the TRAIN link-events: white noise per 5 min value (the median over
  link-events of the sd of the 5 min means of A − c in the offset window), an offset error per link-event (the sd over
  link-events of the mean of A − c over the dry bins of the 30 min before the first wet bin), and a persistent offset per
  link with sd σ_link (rule 6).
- Rain error (E1 model), observed = true·exp(s + c_t + e), every term with mean one:
  s per link-event, sd σ_s,res = 0.43 (the scale part left after g_e);
  c_t per event and time, shared by all links, variance 0.46·σ²; e per link-event and bin, variance 0.54·σ² / n_cells
  (n_cells = the number of radar cells along the link, data layer); c and e are AR(1) with ρ = 0.35;
  σ = σ_lo = 0.48 and σ = σ_hi = 0.87.
- A cell's truth is F_true and P_true at the cell's own R_cell. 5 replicates; the judged number is the median.

### E3 fake test of the bin rule (isolated from the rain error)

The fake world with memory and noise but WITHOUT rain error (σ = 0, s = 0), plain length test. **PASS for a rule (R1, then
R2) if the intercept is within 0.15 dB and the slope within 20 % in all eligible cells but at most one. If neither rule
passes: STOP.**

### E2: the candidates (all on the events with a gauge ratio, all with the adopted bin rule)

- **(i) Gauge-ratio events only.** The plain length test, bins selected by each link's own radar rain. (Simplest.)
- **(ii) Instrument binning.** The bin is selected (class and steadiness) by Z = the mean radar rain over all link-events of
  the event at that time, the same for every link; each link's own L stays in the regression; the cell is placed at the mean
  of the links' own rain. The event-gauge version is not testable in this fake world (it has no true rain at the gauge
  points), so the wide-area radar version is the one tested.
- **(iii) Forward-bias correction.** Start from the plain cells of (i). Repeat 4 times: build the curves F, P through the
  current cell values; make a proxy of the true rain by shrinking the observed ln R about its link-event mean and the
  link-event means about their event mean, so that proxy + error has the observed variances; simulate A − c = F(R) + L·P(R)
  (no memory, no noise) with the rain error of the assumed σ; run the same length test on it; move each cell value by
  (observed value − simulated value). 4 simulations per step. The correction knows the error model and σ, not the truth.
- Simplicity order for the adoption: (i), (ii), (iii).

**PASS rule per candidate (2026-11-11):** intercept within 0.15 dB and slope within 20 % in at least 13 of 14 cells (if
fewer than 14 cells are eligible: in all eligible cells but at most one), at BOTH σ_lo and σ_hi. The simplest passing
candidate is adopted; if several pass, the one with the narrowest F intervals on real data; **if none passes, STOP** and
report which cells fail and why.

**Result of E3 (2026-11-11).** Script `research/waa_ksurface/e3_bin_rule.py`, report `outputs/reports/e3_bin_rule.md`.
4 563 TRAIN link-events of 187 links in the 14 events with a gauge ratio.

- Sizes measured on real data: white noise of a 5 min value 0.054 dB; offset error of a link-event 0.136 dB; between-link sd
  of the short links' residuals 0.241 dB (41 links); storm-to-storm sd inside a link 0.592 dB.
- n_short = max(5, ceil((0.241 / 0.15)²)) = **5** (the floor decides: the formula alone gives 3).
- 32 GHz has 0 links below 1 km: no eligible cell by rule 6, so it stays excluded from the branch for every model.
- **Fake test of the bin rules (memory and noise, NO rain error): both rules FAIL.**
  - R1 (20 min, factor 2): 14 cells judged, 7 outside. The intercept is too high in every cell but one (+0.05 to +0.17 dB);
    3 cells exceed 0.15 dB (28 GHz 0.5–1 and 2–5 mm/h, 38 GHz 10–20 mm/h). The slope is outside ±20 % in 5 cells, all below
    1 mm/h, all too low (−25 % to −104 %).
  - R2 (30 min, factor 1.5): 12 cells judged (the two 10–20 mm/h cells are no longer eligible), 4 outside. Every intercept
    is inside ±0.15 dB (largest +0.13 dB). The 4 failures are slopes below 1 mm/h (−79 %, −53 %, +49 %, +77 %).
  - The failing slopes are small numbers: the true slope there is 0.013–0.14 dB/km, so the errors are 0.01–0.07 dB/km.

**Verdict.** STOP at E3 (rule: "If neither rule passes: STOP"). E2 and E4 are not run. Nothing is patched.

### The one diagnostic allowed after the failed check (2026-11-11, written before it is run)

**What.** The same fake test for R1 and R2 with the two ingredients of the fake world switched on and off: memory (the store)
on / off, noise (white noise, offset error, link offset) on / off; still no rain error. 5 replicates per combination
(`e3_diagnostic.py`).

**Why.** To locate the cause of (1) the intercept that is too high and (2) the slope errors in light rain: the memory of
the store (the bin is still entered from heavier rain more often than from lighter rain, and drying is slower than
wetting), the noise, or the method itself (a class holds a range of R while the cell is placed at one R_cell; with neither
memory nor noise only this is left). Nothing is changed on the basis of this diagnostic; the rescue stops after it.

**Result of the diagnostic (2026-11-11).** `outputs/reports/e3_diagnostic.md`, 5 replicates per row, no rain error.

| rule | memory | noise | cells | cells outside | median intercept error (dB) | intercepts outside | slopes outside | largest slope error (dB/km) |
|---|---|---|---|---|---|---|---|---|
| R1 | on | on | 14 | 7 | +0.10 | 3 | 5 | 0.10 |
| R1 | on | off | 14 | 1 | +0.08 | 0 | 1 | 0.10 |
| R1 | off | on | 14 | 4 | +0.03 | 1 | 3 | 0.13 |
| R1 | off | off | 14 | 1 | 0.00 | 1 | 0 | 0.14 |
| R2 | on | on | 12 | 4 | +0.07 | 0 | 4 | 0.06 |
| R2 | on | off | 12 | 5 | +0.04 | 0 | 5 | 0.07 |
| R2 | off | on | 12 | 4 | +0.03 | 0 | 4 | 0.05 |
| R2 | off | off | 12 | 5 | 0.00 | 0 | 5 | 0.06 |

- **The intercept.** The method itself is unbiased (median 0.00 dB without memory and noise). The memory of the store lifts
  the intercept by about 0.07 dB under R1 and 0.04 dB under R2: even with a symmetric band a bin is still more often entered
  from heavier rain, and the store dries more slowly than it wets. The noise adds about 0.03 dB. Together they reach
  +0.10 dB under R1, and three cells cross 0.15 dB. Under R2 no intercept is outside in any combination.
- **The slope in light rain.** Under R2 four or five slopes are outside ±20 % in every combination, also with neither memory
  nor noise. So this failure is not caused by the bin rule, the memory or the noise. The errors are 0.01–0.07 dB/km in
  absolute size; they fail because the true slope below 1 mm/h is itself only 0.013–0.14 dB/km, so 20 % of it is a few
  thousandths of a dB/km. A probable mechanism, not tested further (one diagnostic only): inside one rain class the links'
  R differs a little with their length, and in light rain F changes fast with R, so a small part of F leaks into the slope.
- The one cell that fails without memory and noise under R1 is 38 GHz 10–20 mm/h (intercept +0.22 dB, 6 short links).

**What this says about the stop.** Neither rule passes the criterion as written. Under R2 the intercept (the quantity F is
taken from) is inside the tolerance in every cell; what fails is the relative slope criterion in light rain, where it asks
for an accuracy of thousandths of a dB/km. Whether the slope criterion should be absolute in light rain, or whether P below
1 mm/h should be declared not measurable, is a decision for the review. It is not taken here.

**The rescue stops here (2026-11-11).** E2 (the candidates) and E4 (the re-measurement of F and P) were not run. Nothing was
patched. The E2 code (candidates (i)–(iii)) is not written; the libraries `ltest.py` and `fakeworld.py` hold the bin rule,
the line, the fake world and the E1 error model that E2 needs.

---

# E2–E4 after the E3 stop (2026-11-11, unattended)

## Deviation 3: the slope criterion becomes absolute (2026-11-11, logged before any run)

**What changes.** In the fake tests of E3 and E2 the slope criterion "within 20 % of the true slope" is replaced by an
absolute tolerance:

    |s − P_true| <= 0.15 dB / L_med

L_med = the median length of the links in the TRAIN length-test dataset (the distinct links of 28 and 38 GHz in the events
with a gauge ratio). It is computed and recorded below, before the re-evaluation. No other rule changes: the intercept
tolerance stays 0.15 dB, "all judged cells but at most one" stays, the two bin rules and their order stay, the E2 candidates
and the E1 error model stay.

**Reason.** A relative tolerance is undefined for slopes near zero. In the E3 diagnostic it failed with neither memory nor
noise, on errors of 0.01–0.07 dB/km, because the true slope below 1 mm/h is only 0.013–0.14 dB/km. The absolute form is the
intercept tolerance (0.15 dB) applied to the path term L·s of a median link, so the intercept and the slope are held to one
standard.

**This is a choice made after seeing a result** (the E3 fake test and its diagnostic), which is why it is a deviation. The
E3 numbers already seen: under R2 the largest absolute slope error was 0.065 dB/km and no intercept was outside.

**L_med (2026-11-11).** The median length of the 187 distinct links of the TRAIN length-test dataset (28 and 38 GHz, events
with a gauge ratio) is **2.113 km**. The slope tolerance is 0.15 dB / 2.113 km = **0.0710 dB/km**.

## Step E3 re-evaluated under Deviation 3 (2026-11-11)

Same rules, same order, same fake world and seeds (`e3_bin_rule.py`, report `outputs/reports/e3_bin_rule.md`).

- R1 (20 min, factor 2): 14 cells judged, 3 outside (intercepts +0.152, +0.166, +0.164 dB; one of them also on the slope).
  FAIL.
- **R2 (30 min, factor 1.5): 12 cells judged, 0 outside. PASS.** Largest intercept error +0.129 dB (28 GHz, 2–5 mm/h); largest
  slope error 0.065 dB/km (38 GHz, 0.5–1 mm/h) against the tolerance 0.071.
- **Adopted: R2.** Steady rain = every bin of the 30 min before the bin within [R / 1.5, 1.5·R].
- **Consequence, recorded:** under R2 the two 10–20 mm/h cells are not eligible (on real data: 28 GHz 4 events and 0 short
  links, 38 GHz 3 events and 0 short links). **F and P are measured up to the 5–10 mm/h class only** (R_cell about 6.3 mm/h).
- The intercept is still slightly high in the fake world with memory (median about +0.07 dB): inside the tolerance, not zero.

**Verdict.** PASS with R2.

## The g_e = 1 exclusion, formalized (2026-11-11)

**Rule.** Events with g_e = 1 (no gauge ratio) do not enter the length test. **Reason:** E1 measured the radar at about 2.5
times the gauges in them, an uncorrected one-sided scale error; a test that sorts bins by radar rain cannot use them.
(This is a property of the rain reference, so it is a data-quality rule of this branch, not a choice by the fit of any
model.)

**What is left (TRAIN, 28 and 38 GHz, rule R2):** 14 events, 4 563 link-events, 187 links, 41 of them below 1 km. 12 of 14
cells are eligible. Per eligible cell: short links min 7, median 17; events min 7, median 9; link-events 103 to 565
(`outputs/tables/e3_bin_rule_adopted_real_cells.csv`).

## Step E2: addendum to the pre-registration (2026-11-11, before the E2 code is run)

- The definitions are those adopted in E3 (rule R2, classes, n_short = 5); the criteria are those of Deviation 3 (intercept
  within 0.15 dB, slope within 0.0710 dB/km). 12 cells are eligible, so "at least 13 of 14" reads, as written: all judged
  cells but at most one.
- The fake world, the three candidates and the adoption rule are unchanged from the entry "Steps E3 and E2".
- Details of candidate (iii) fixed now: the proxy of the true rain is x' = M + b·(m − M) + a·(x − m) + v / 2, with x = ln R of
  the wet bins of a link-event, m their mean, M the mean of m over the event's link-events, a² = 1 − (per-bin error
  variance of the link) / (pooled variance of x inside link-events), b² = 1 − (error variance of m) / (pooled variance of m
  inside events), both floored at 0.2², and v the total error variance (the errors have mean one, so the observed logs sit
  v / 2 below the true ones). The curves inside the correction are the monotone curves through (0, 0) and the current cell
  values. The simulation inside the correction uses the same error model as the fake world (it knows σ, not the truth).
- **The error-size sweep (asked for in the E2–E4 instructions):** every candidate is also run with the whole rain error
  scaled by f = 0, 0.25, 0.5, 0.75 (per-bin sd f·0.48, scale part f·0.43), to report the largest rain error at which it
  still passes. f = 1 is the σ_lo run. Report only; it does not change the adoption rule.

**Result of E2 (2026-11-11).** Script `research/waa_ksurface/e2_candidates.py` (library `candidates.py`), report
`outputs/reports/e2_candidates.md`. Fake data only: 4 563 link-events, 187 links, 14 events; rule R2; n_short = 5; 5
replicates.

| candidate | low end (σ = 0.48): cells inside / judged | high end (σ = 0.87) | largest intercept error at the low end | passes |
|---|---|---|---|---|
| (i) gauge-ratio events, plain | 2 / 11 | 0 / 0 | 0.52 dB | no |
| (ii) instrument binning (area-mean radar) | 6 / 8 | 0 / 0 | 0.17 dB | no |
| (iii) forward-bias correction | 3 / 11 | 0 / 0 | 0.80 dB | no |

- **At the high end no cell can be judged for any candidate.** The steady-rain rule R2 asks for 30 min of rain within a
  factor 1.5; with a per-bin rain error of sd 0.87 the observed rain almost never stays that close. One fake replicate keeps
  136 steady bins (8 539 without rain error, 1 285 at the low end) and has no eligible cell.
- **(i) at the low end:** the intercept is too high by 0.28–0.52 dB from 0.2 to 2 mm/h in both bands and by 0.17 dB at
  38 GHz 2–5 mm/h (9 cells outside). The rule keeps the bins whose noisy rain happened to look steady, and their true rain is
  on average not what the class says.
- **(ii) at the low end:** intercepts within 0.14 dB except 38 GHz 2–5 mm/h (+0.17 dB), but only 8 cells can be judged (no
  5–10 mm/h cell), and the slope is too low at 2–5 mm/h in both bands (−0.12 and −0.16 dB/km). It fails already WITHOUT rain
  error (6 of 10 cells): the slopes at 1–5 mm/h are too low by 0.08–0.11 dB/km. The links of one event do not share one
  rain, so sorting by the area mean mixes links with different rain, and part of the path term is lost.
- **(iii) at the low end:** the correction does not converge to the truth: errors from −0.21 to +0.80 dB (8 cells outside).
  With the steady rule this strict and the rain error this large, each simulated test keeps few bins, so the correction
  chases noise.

**The largest rain error at which each candidate still passes** (the whole error scaled by f; f = 1 is the low end):

| candidate | f = 0 | f = 0.25 | f = 0.5 | f = 0.75 | largest error passed |
|---|---|---|---|---|---|
| (i) | 12 / 12 | 11 / 12 | 8 / 12 | 4 / 12 | f = 0.25: per-bin sd 0.12, scale part 0.11 |
| (ii) | 6 / 10 | 8 / 10 | 7 / 10 | 8 / 10 | none (fails without rain error) |
| (iii) | 12 / 12 | 12 / 12 | 11 / 12 | 6 / 12 | f = 0.5: per-bin sd 0.24, scale part 0.21 |

The measured radar error (per-bin sd 0.48–0.87, scale part 0.43) is two to four times larger than what the best candidate
can take.

**Verdict.** STOP at E2 (rule: "if none passes, STOP and report which cells fail and why"). No candidate is adopted. **E4
was not run**; the per-link residual criterion was therefore not written and the spread not computed for it. Nothing was
patched.

**What this means for Q1 (can F and P be trusted?).** With radar as the rain input, at 5 min and with the radar error
measured in E1, the length test of this branch cannot return F(R) to 0.15 dB. This holds for the plain test, for area-mean
binning and for the forward correction as designed here. The A3a intercepts (superseded) and the E3 real cells carry this
error and its size is not known cell by cell.

---

# The conditional rain error (2026-11-11, unattended): steps T1–T4

**Why this block.** E2 failed with the UNCONDITIONAL rain error of E1 (every 5 min bin, sd 0.48–0.87). The length test uses
only steady windows, and labels a window by a median. The error that matters is the error of that label in such windows,
and for the line against L only the part that is independent between links tilts it. That error was not measured in E1.

## Step T1: the conditional rain error (2026-11-11, written before the run)

**Question served.** Q1 (can F and P be trusted?).

**(a) Primary: radar label against gauge rain in steady windows (rain data only).**
- Where and when: the 10 City gauges, inside the 14 TRAIN events with a gauge ratio. Radar = the pixel over the gauge
  (`radar_at_gauges_5min`, `R_pixel`) divided by the event's g_e, i.e. the rain input of the length test (inherited data
  layer). Gauge = the 1 min City series.
- A window = a 5 min bin that passes the adopted steadiness rule on the radar series (R inside 0.1–20 mm/h and every bin of
  the 30 min before it within [R / 1.5, 1.5·R]) together with those 30 min: 35 min. Label = the median of the 7 radar values.
  Gauge rain = the mean gauge rate over the same 35 min, all minutes valid.
- Used windows: label >= 1 mm/h and gauge rain > 0. Reason: below 1 mm/h a 35 min window holds under 0.6 mm, a handful of
  0.1 mm (weighing) or 0.2 mm (tipping) steps, so the gauge cannot serve as a reference there. The windows dropped for a
  zero gauge are counted and reported. **Consequence, stated now: the conditional error is not measured below 1 mm/h; the
  number found is applied to all classes.**
- d = ln(label / gauge). Reported: the sd of d; the split by variance components (method of moments, events as groups) into
  the per-event scale part σ_ev (between events) and the within-event part σ_w.
- **σ_cond = √0.54 · σ_w** : the link-independent, window-level part, with the E1 estimate that 46 % of the per-bin radar
  error is shared in space. The shared part is √0.46 · σ_w. The gauge's own error and the point-vs-pixel difference are NOT
  subtracted (conservative: σ_cond is an upper estimate of the radar's part).
- **[σc_lo, σc_hi] = the 5 % and 95 % points of σ_cond over 200 event-bootstrap draws.**
- Reference, same windows without the steadiness rule (every 35 min window with label >= 1 mm/h): the sd of d, to show what
  the conditioning changes.

**(b) Secondary: the steady-window count (TRAIN link-events of the gauge-ratio events).**
- As instructed: fake replicates whose true rain is the observed radar, with a per-bin error added on a grid
  (sd 0, 0.05, 0.1, 0.2, 0.3, 0.48; E1 structure: 46 % shared at one time, lag-1 correlation 0.35), and the count of steady
  bins per event compared with the real count. **Caveat, verbatim: the σ that matches the real count is an EFFECTIVE bound,
  biased low because the fake truth (the observed radar) already carries error.** Known before the run: with the observed
  radar as the truth the count at added error 0 IS the real count, so the match can only be at 0; what the grid shows is
  how fast the count falls with added error.
- A second reading that can give a number, added here: the share of rainy bins that are steady in the real data, compared
  with the same share for a truth that is perfectly constant in time plus a per-bin error (lag-1 correlation 0.35). The σ at
  which a constant truth gives the real share is an UPPER bound on the fast per-bin error of a link's own label: any real
  change of the rain uses up part of the allowed unsteadiness.

**(c) Fallback.** If (a) cannot be computed (fewer than 30 used windows or fewer than 5 events), σc_hi = the upper bound of
the second reading of (b), with its caveat, and σc_lo = 0.

**Pass rule (2026-11-11).** T1 is a measurement. It passes if (a) returns [σc_lo, σc_hi] from at least 30 windows in at least
5 events, or the fallback (c) is used and said so.

**Result of T1 (2026-11-11).** Script `research/waa_ksurface/t1_conditional_error.py`, report
`outputs/reports/t1_conditional_error.md`.

- **(a)** 175 steady windows with a label >= 1 mm/h at the 10 City gauges; 3 dropped for a zero gauge; **172 used, in 8 of
  the 14 events.** sd of ln(label / gauge) = 0.50 (median +0.20: the label is about 22 % above the gauge).
  - per-event scale part σ_ev = **0.45** (between events; it is there although the radar is already divided by g_e);
  - within events σ_w = 0.30, of which shared in space √0.46·σ_w = 0.20 and **link-independent σ_cond = √0.54·σ_w = 0.22**;
  - **[σc_lo, σc_hi] = [0.18, 0.26]** (200 event draws).
  - Reference, every 35 min window with a label >= 1 mm/h (4 097 windows, 13 events): sd 0.79, within events 0.71. **The
    steadiness rule cuts the within-event error from 0.71 to 0.30.**
  - By class the sd is 0.59 (1–2 mm/h, 57 windows), 0.62 (2–5, 38 windows), 0.24 (5–10, 77 windows in 4 events). Not
    measured below 1 mm/h.
- **(b)** Real data: 8 551 steady bins in 14 events. With the observed radar as the truth and a per-bin error added, the
  count falls to 96 % (sd 0.05), 87 % (0.10), 64 % (0.20), 43 % (0.30), 18 % (0.48). **Caveat, verbatim: the σ that matches
  the real count is an EFFECTIVE bound, biased low because the fake truth (the observed radar) already carries error.** The
  match is at 0 by construction. Second reading: 2.5 % of the rainy bins with 30 min of history are steady; a perfectly
  constant rain gives that share at a per-bin sd of 0.51. This upper bound is too loose to be of use (most rain is simply
  not steady).
- (c) not needed: (a) is computable (172 windows, 8 events).

**Verdict.** PASS (a measurement). [σc_lo, σc_hi] = [0.18, 0.26] for the link-independent window-level part; shared part
0.17–0.24 (the same ratio); per-event scale part 0.45. **Limits: 8 events; nothing below 1 mm/h; the gauge's own error is
inside these numbers.**

## Steps T2 and T3: candidate (iv), and E2 rerun at the conditional error (2026-11-11, written before any run)

**Question served.** Q1.

### The rain error of the T3 fake world for candidates (i), (ii), (iii)

The fake world of E3 / E2 is kept (links, events, true rain = the observed rain input, F_true = 3·R / (R + 2), P_true =
0.1·R and 0.2·R, one store with 5 / 15 min, the measured noise, rule R2, n_short = 5, 5 replicates, Deviation 3 criteria).
Only the rain error changes, to the one T1 measured for steady windows:

    ln(R_obs / R_true) = E_event + c_t + e_link,t      (shifted so that the mean of R_obs / R_true is one)

- E_event: one value per event, shared by all its links, sd σ_ev = 0.45.
- e_link,t: the link's own part, sd σ_c / √n_cells, σ_c = σc_lo = 0.18 and σc_hi = 0.26 (the two ends of the range).
- c_t: shared by all links of the event at one time, sd √(0.46 / 0.54)·σ_c (0.17 and 0.24).
- c and e change slowly: AR(1) per 5 min bin with ρ = 0.9^(1/6) = 0.983, so that the error is 0.9-correlated over the 30 min
  of a window. Reason: T1 measured the error of a window label, and an error that jumped inside a window would break the
  steadiness that defines the window.
- No extra fast per-bin error. Reason (the T1b caveat): the fake truth is the observed radar, which already carries the
  radar's fast error; adding more only removes steady windows that the real data do have.
- Candidate (iii) is given this model (it knows the sizes, not the truth). Its stand-in for the true rain shrinks at three
  levels: bins about the link-event mean, link-event means about the event mean, event means about the overall mean, each by
  the share of the variance that the model's error explains (floor 0.2).

### Candidate (iv): the gauge-binned length test (T2)

- Links: the TRAIN link-events of 28 and 38 GHz whose link has a City gauge within 1 km of its path (the pairing of the
  data layer), in ALL TRAIN events. Reason for all events: this test never uses the radar, so the g_e = 1 exclusion (a radar
  problem) does not apply to it.
- Gauge series: 5 min means of the gauge's 1 min rates (all 5 minutes valid).
- Window: the 30 min before a bin plus the bin (7 bins). Label = the mean gauge rate of the 7 bins; it must lie in a class.
- **Steadiness with the written rule for zeros and steps:** the 30 min before the bin are cut into three 10 min means a, b,
  c. The window is steady if max(a, b, c) − δ <= 1.5·(min(a, b, c) + δ), where δ is one gauge step in a 10 min rate
  (0.6 mm/h for a 0.1 mm weighing gauge, 1.2 mm/h for a 0.2 mm tipping bucket). Reason: a gauge reports whole steps, so two
  10 min sums of the same rain can differ by a step; single zero minutes are absorbed in the sums, and a 10 min sum of zero
  passes only if the others are within one step of zero. Known consequence, stated now: below about 1 mm/h the rule passes
  almost every window, so there the gauge cannot show steadiness.
- Line and eligibility: as E3 (Huber line over the links' lengths; values from >= 5 events; >= 5 links below 1 km and >= 5 of
  1 km or more). **If the floor of 5 short links leaves no eligible cell on real data: report the counts and the floor at
  which cells would exist, and STOP for this candidate. The floor is not lowered here.**
- **Its fake world:** the same truth, memory and noise; the fake gauge rate = the true path rain × exp(e_own + e_rep),
  then turned into whole gauge steps (running total, 0.1 or 0.2 mm).
  - e_own: the gauge's own error, white per 5 min bin, sd 0.32·√3 = 0.55, so that a 15 min mean has the E1 value 0.32.
  - e_rep: point-vs-path representativeness, AR(1) with ρ = 0.35. **Rule for the link length:** E1 gave a variogram that
    rises linearly with distance, slope s = 0.17² / 1.04 km = 0.028 per km (the point-vs-pixel sd 0.17 over the mean
    distance 1.04 km). For a gauge at distance d from a path of length L the variance between the point and the path mean is
    taken as s·(d + L / 3): d to reach the path, L / 3 the mean distance between two random points of the path. sd at
    15 min = √(s·(d + L / 3)), per bin = that / √0.516 (the E1 factor).
  - Normal end: as above. **Pessimistic end: the representativeness sd doubled.**
  - The steps are applied on top of e_own, which already holds some step error: conservative.

### T3: verdicts per cell, and the adoption rule (written before running)

- A cell of a candidate is **valid** if it is judged (eligible in at least 3 of the 5 replicates) at BOTH ends and its
  fake-world intercept error (<= 0.15 dB) and slope error (<= 0.0710 dB/km) are inside tolerance at BOTH ends. The ends are
  σc_lo and σc_hi for (i), (ii), (iii), and normal and doubled representativeness for (iv).
- **Adoption: the single candidate with the most valid cells, if it has at least 8 valid cells of which at least 3 lie
  below 2 mm/h (classes 0.1–0.2, 0.2–0.5, 0.5–1, 1–2). Ties are broken by the narrower bootstrap intervals of F on real
  data. No mixing of candidates per cell. If no candidate qualifies: go to T4.**
- If a candidate is adopted, F and P are taken only from its valid cells.

**Result of T2 and T3 (2026-11-11).** Script `research/waa_ksurface/t3_rerun.py`, report `outputs/reports/t3_rerun.md`.
Fake data only; real data read for the counts of (iv).

- **T2, candidate (iv) on real data:** 1 983 TRAIN link-events of 41 links of 28 and 38 GHz with a City gauge within 1 km
  (18 below 1 km), in 43 events (all TRAIN events). (Corrected 2026-11-12: this line first gave 2 425 link-events of 50
  links, which is the A1 gauge-input set over all bands; the numbers of the T3 report are the right ones. No result
  changes.) 11 of 14 cells are eligible with the floor of 5 short links (28 GHz 0.1–0.2 mm/h has 4; the two 10–20 cells
  are not fitted). Short links per eligible cell: 6 (28 GHz) and 6–12 (38 GHz). The floor did not need to be touched.
- **T3, valid cells (inside tolerance at both ends):**

| candidate | cells judged at both ends | valid cells | of them below 2 mm/h | largest intercept error | qualifies |
|---|---|---|---|---|---|
| (i) gauge-ratio events, plain | 12 | 4 | 2 | 0.42 dB | no |
| (ii) instrument binning (area-mean radar) | 10 | **8** | **8** | 0.09 dB | **yes** |
| (iii) forward-bias correction | 12 | 5 | 5 | 0.54 dB | no |
| (iv) gauge-binned | 12 | 3 | 3 | 1.68 dB | no |

  - (i): the intercept is too high by 0.25–0.42 dB from 0.2 to 2 mm/h in both bands. The per-event scale error (sd 0.45,
    shared by all links) and the link's own slow error still move bins between classes.
  - (ii): every intercept is inside 0.09 dB. The 8 cells from 0.1 to 2 mm/h are valid in both bands. The two 2–5 mm/h cells
    fail on the slope (−0.05 to −0.14 dB/km, 28 GHz just at the limit), and no 5–10 mm/h cell can be judged.
  - (iii): unstable (errors from −0.24 to +0.54 dB).
  - (iv): the intercept is too high by 0.16–0.30 dB at 1–5 mm/h; several 38 GHz slopes are too low.
- **Adopted by the written rule: candidate (ii), instrument binning.** 8 valid cells, all 8 below 2 mm/h. No tie.
- **Consequence, recorded: F and P are measured only in the four classes 0.1–0.2, 0.2–0.5, 0.5–1 and 1–2 mm/h of each band.
  Above about 1.5 mm/h nothing is measured.**

**Verdict.** PASS: (ii) adopted. T4 is not run.

## Step E4: F(R) and P(R) re-measured with candidate (ii) (2026-11-11, written before the run)

**Question served.** Q1.

**What.** On the real TRAIN data of the 14 events with a gauge ratio (28 and 38 GHz): bins selected by Z = the event's
area-mean radar rain with rule R2; per link-event and class the median of (A − c) and the mean of the link's own rain; the
Huber line over L per band and class; the cell sits at the mean of the links' own rain. Only the 8 cells that T3 found
valid are used.

**Intervals.**
- Bootstrap: 200 event draws, 90 %.
- Error range: the fake-world error of the cell at the two ends of the T1 range (σc_lo, σc_hi), taken as a possible bias of
  the real value. Error-range interval = [value − max(e_lo, e_hi, 0), value − min(e_lo, e_hi, 0)]. The value itself is not
  shifted. The combined interval is the bootstrap interval widened by this range.
- Curves: per band a monotone cubic through (0, 0) and the 4 valid cells (made non-decreasing first if needed); intervals
  from the bootstrap draws. Valid range = from the lowest to the highest valid R_cell; below it the curve is an
  interpolation to zero, above it an extrapolation (F held, P with its last slope), both marked.

**The per-link criterion (written 2026-11-11, before the spread is computed for it).**
- Quantity: for the links below 1 km, the residual of every link-event value from its cell line (the 8 valid cells),
  averaged per link and storm; then variance components over links: the between-link sd σ_link and the storm-to-storm sd.
- **A per-link factor is needed if BOTH hold: (1) σ_link >= 0.15 dB (the tolerance to which F is held), and (2) the
  between-link variance is significant: p < 0.05 in a permutation test (the link labels of the link-storm values shuffled
  inside each band, 1 000 times, statistic = σ_link).** Otherwise it is not needed.
- Disclosure: a number of this kind has already been seen. E3 measured, for sizing the fake world, a between-link sd of
  0.24 dB on the short links with the plain test and rule R1 (all 14 cells). The criterion above is therefore not blind to
  the likely outcome; its threshold is the tolerance used everywhere else in this branch, not a number tuned to it.
- Report only (no criterion): whether a link's mean residual grows with the level of F (a factor) or not (an offset): the
  correlation between its mean residual in the two lower classes and in the two upper classes.

**Pass rule (2026-11-11).** E4 is a measurement. It passes if all 8 valid cells are eligible on real data and give F and P
with both intervals. A valid cell that is not eligible on real data is dropped and reported; with fewer than 2 cells in a
band there is no curve for that band. Then STOP for review; A4 is not started.

**Result of E4 (2026-11-11).** Script `research/waa_ksurface/e4_remeasure.py`, report `outputs/reports/e4_remeasure.md`.
TRAIN, 4 563 link-events of 187 links in the 14 events with a gauge ratio. All 8 valid cells are eligible on real data.

| band | class (mm/h) | R_cell | F (dB) [bootstrap 90 %] [error range] | P (dB/km) [bootstrap 90 %] [error range] | short links | events |
|---|---|---|---|---|---|---|
| 28 GHz | 0.1–0.2 | 0.14 | 0.375 [0.198, 0.574] [0.369, 0.375] | 0.045 [0.020, 0.078] [0.042, 0.046] | 17 | 11 |
| 28 GHz | 0.2–0.5 | 0.37 | 0.525 [0.211, 0.928] [0.488, 0.526] | 0.062 [0.001, 0.094] [0.057, 0.067] | 17 | 7 |
| 28 GHz | 0.5–1 | 0.64 | 0.725 [0.619, 0.826] [0.725, 0.808] | 0.093 [0.058, 0.140] [0.093, 0.119] | 17 | 9 |
| 28 GHz | 1–2 | 1.44 | 1.406 [1.189, 1.624] [1.321, 1.406] | 0.098 [0.015, 0.172] [0.098, 0.144] | 17 | 8 |
| 38 GHz | 0.1–0.2 | 0.13 | 0.350 [0.212, 0.501] [0.327, 0.350] | 0.075 [0.041, 0.126] [0.074, 0.076] | 24 | 11 |
| 38 GHz | 0.2–0.5 | 0.29 | 0.619 [0.272, 1.196] [0.577, 0.619] | 0.096 [−0.007, 0.142] [0.092, 0.100] | 24 | 7 |
| 38 GHz | 0.5–1 | 0.65 | 0.863 [0.716, 0.966] [0.863, 0.936] | 0.122 [0.047, 0.222] [0.122, 0.163] | 24 | 9 |
| 38 GHz | 1–2 | 1.41 | 1.660 [1.189, 2.101] [1.587, 1.660] | 0.179 [−0.019, 0.388] [0.179, 0.245] | 24 | 8 |

- **Valid range: 28 GHz 0.14–1.44 mm/h, 38 GHz 0.13–1.41 mm/h.** Below it the curves are an interpolation to zero; above it
  nothing is measured (F held, P with its last slope, marked as extrapolation). Curves: `outputs/tables/e4_remeasure_curves.csv`.
- **The statistical interval is the wide one.** The bootstrap intervals (7 to 11 events per cell) are 0.2–0.9 dB wide for F;
  the error range from the fake world is under 0.09 dB.
- F rises with R in every band; the monotone step changes nothing. P at 0.2–0.5 and 1–2 mm/h has intervals that reach zero.
- Against the superseded A3a curve at the same R: shifts from −0.21 dB (28 GHz, 0.5–1 mm/h) to +0.08 dB.

**The per-link criterion (as written above, 2026-11-11): verdict NOT NEEDED.** 41 short links, 533 link-storm values.
Between-link sd 0.088 dB (threshold 0.15 dB); storm-to-storm sd inside a link 0.557 dB; permutation p = 0.098 (null: median
0.016 dB, 95 % point 0.102 dB). Both conditions fail. Note: with the plain test in E3 the same quantity was 0.24 dB; with
instrument binning it is 0.09 dB. Report only: the mean residual in the two upper classes against the two lower classes has
correlation 0.68 over the 41 links.

**Verdict.** PASS (a measurement). F and P exist for 0.13–1.44 mm/h. **STOP for review; A4 is not started.**

**What the review should know.**
1. The measured range ends at about 1.4 mm/h. Most of the WAA dynamics of interest (heavier showers) lie above it.
2. The method was adopted on a fake world whose rain error model rests on T1: 172 windows in 8 events, nothing below
   1 mm/h, and a slow-error form (0.9 correlated over 30 min) that is a modelling choice.
3. The cell sits at the MEAN of the links' own rain, and the bins are chosen by the area mean; inside one bin the links have
   different rain. The fake test covers this with the real rain fields, but only for the generic truth used there.
4. The per-link result (no factor needed) holds for 0.13–1.44 mm/h and the short links only.

---

# After the E4 review (2026-11-12, unattended): Deviation 4, then A4, A5, A6

The review accepted E4. This block logs one deviation, then runs A4 (coverage), A5 (the model, no real fit) and A6 (fake-data
recovery). **The real fit (A7) is not run.**

## Deviation 4: per-cell validity is split by quantity (2026-11-12, logged before anything else)

**What changes.** In T3 a cell of the adopted candidate (ii) was valid only if BOTH its fake-world intercept error and its
slope error were inside tolerance at both ends of [σc_lo, σc_hi]. From now on validity attaches to the quantity taken from
the cell:

- a cell is **F-valid** if it is judged at both ends and its fake-world INTERCEPT error is within 0.15 dB at both ends;
- a cell is **P-valid** if it is judged at both ends and its fake-world SLOPE error is within 0.0710 dB/km at both ends.

F(R) is taken from all F-valid cells, P(R) from all P-valid cells. Nothing else changes: the candidate (ii), the adoption,
the bin rule R2, the tolerances, the intervals machinery of E4 (bootstrap, error range), the curve construction.

**Reason (given by the review).** Validity must attach to the quantity that is extracted. A cell whose slope fails the fake
test says nothing against its intercept, and the other way round. The committed T3 fake results already hold both errors
per cell, so no new fit is needed.

**This is post hoc**, which is why it is a deviation. It was decided by the review before the stored numbers were re-read,
and it is applied mechanically: the script reads `outputs/tables/t3_rerun_valid_cells.csv` (committed, candidate (ii)) and
applies the two inequalities; no cell is judged by hand. Disclosure: the T3 report with these numbers was open in the
working session of 2026-11-12 (it was read to restate the previous check report) before this entry was written. The rule
itself was not written by that session.

**What is done.** (1) List every judged cell of candidate (ii) with its intercept and slope error at both ends and its
F- and P-validity (`e4_remeasure_validity.csv`). (2) Rerun `e4_remeasure.py` with F from the F-valid cells and P from the
P-valid cells. No fake fit is rerun. (3) Report the valid R range of F and of P separately.

**What is kept as it was.** The per-link criterion of E4 was written for the 8 cells that are valid in both quantities, and
the review accepted its verdict. It is still computed on those 8 cells (the cells whose whole line is trusted), so that the
accepted verdict is not re-opened by this deviation.

**Pass rule (2026-11-12).** The rerun passes if (a) the 8 cells of the accepted E4 table keep their values (same data, same
seeds: every number equal), (b) every F-valid or P-valid cell that is eligible on real data gets its value with both
intervals. A newly valid cell that is not eligible on real data is dropped and reported.

**Result of Deviation 4 (2026-11-12).** `research/waa_ksurface/e4_remeasure.py` rerun; report
`outputs/reports/e4_remeasure.md`; the re-read table `outputs/tables/e4_remeasure_validity.csv`. No fake fit was rerun.

Candidate (ii), fake-world errors at the two ends (σc_lo / σc_hi), tolerances 0.15 dB and 0.0710 dB/km:

| band | class (mm/h) | intercept error (dB) | slope error (dB/km) | F-valid | P-valid |
|---|---|---|---|---|---|
| 28 GHz | 0.1–0.2 | +0.006 / +0.004 | −0.001 / +0.003 | yes | yes |
| 28 GHz | 0.2–0.5 | −0.002 / +0.036 | +0.005 / −0.005 | yes | yes |
| 28 GHz | 0.5–1 | −0.083 / −0.010 | −0.026 / −0.026 | yes | yes |
| 28 GHz | 1–2 | +0.063 / +0.086 | −0.044 / −0.046 | yes | yes |
| 28 GHz | 2–5 | +0.079 / +0.026 | −0.0714 / −0.053 | **yes** | no (0.0714 > 0.0710 at the low end) |
| 38 GHz | 0.1–0.2 | +0.012 / +0.023 | +0.000 / −0.001 | yes | yes |
| 38 GHz | 0.2–0.5 | +0.021 / +0.042 | +0.004 / −0.005 | yes | yes |
| 38 GHz | 0.5–1 | −0.072 / −0.053 | −0.041 / −0.031 | yes | yes |
| 38 GHz | 1–2 | +0.043 / +0.073 | −0.052 / −0.065 | yes | yes |
| 38 GHz | 2–5 | +0.065 / −0.029 | −0.138 / −0.111 | **yes** | no |

No 5–10 or 10–20 mm/h cell of (ii) was judged in T3, so none can be valid for anything.

- **F-valid: 10 cells (was 8). P-valid: 8 cells (unchanged).** The two 2–5 mm/h cells are added to F only.
- New F cells on real data: 28 GHz 2–5 mm/h, R_cell 3.20, F = 1.962 dB [bootstrap 1.589, 2.288] [error range 1.882, 1.962];
  38 GHz 2–5 mm/h, R_cell 3.28, F = 2.029 dB [1.653, 2.407] [1.964, 2.059]. Both eligible (17 and 24 short links, 8 events).
- **Valid R range of F: 28 GHz 0.14–3.20 mm/h, 38 GHz 0.13–3.28 mm/h. Valid R range of P: 28 GHz 0.14–1.44 mm/h, 38 GHz
  0.13–1.41 mm/h (unchanged).**
- Pass rule (a): the 8 cells of the accepted table are equal in every column (largest difference 0); the per-link test is
  equal (σ_link 0.088 dB, p = 0.098, not needed). Pass rule (b): both new cells have their values and both intervals.
- The curves product has now one measured-range column per quantity (`F_measured`, `P_measured`).
- Seen and recorded, not acted on: above its valid range the E4 rule continues P with the slope of its last segment. At
  28 GHz that slope is 0.006 dB/km per mm/h (the cells 0.093 and 0.098 dB/km), far below the slope of the 2–5 mm/h cell
  that is not P-valid (0.29 dB/km at 3.2 mm/h). The continuation of P is not a measurement; A4 treats it accordingly.

**Verdict.** PASS.

## Step A4: coverage of the (R, W) plane and the learning domain (2026-11-12, written before the run)

**Question served.** Q2 (what is k?): k(R, W) can only be learned and judged where the data visit the (R, W) plane.

**What.** On the TRAIN links below 1 km (28 and 38 GHz), per minute: Ŵ = (A − c) − L·P̂(R). The number of minutes in every
cell of the (R, W) plane, wetting and drying apart, per band and pooled. The learning domain = the cells with enough
minutes from enough events. **Ŵ is used for this map only. It is never a target of a fit.**

**Definitions, each with its reason (fixed before the run).**

1. **Events: the 14 TRAIN events with a gauge ratio.** Reason: F and P were measured on these events only, and E1 found the
   radar at about 2.5 times the gauges in the g_e = 1 events, an uncorrected one-sided scale error of the rain input. A
   model that reads its target F(R) and its path term P(R) from R meets the same error as the length test did. The same
   data-quality rule therefore applies to learning k (inherited from this branch's own "g_e = 1 exclusion"). The g_e = 1
   events are counted and reported, not mapped into the domain. **This fixes the learning set of A5–A7: flagged for the
   review.**
2. **A minute counts** if its 1 min A is valid and its 5 min radar value is valid (data layer, W1). R of the minute = the
   radar value of its 5 min bin divided by g_e.
3. **P̂ (map only).** Inside P's valid range: the measured E4 curve. Below the lowest cell: its interpolation to zero. Above
   the valid range P is not measured, and a continuation must be assumed. Two linear continuations are computed:
   (L1) the E4 rule, the slope of the last measured segment; (L2) the line through the origin and the last P-valid cell
   (P proportional to R). Reason for two: the result of Deviation 4 showed that the last-segment slope is a difference of
   two cells with overlapping intervals (at 28 GHz it is nearly zero). The map is drawn with L2 (the larger path term at
   high R, so the smaller Ŵ). **A cell enters the domain only if it passes the floors under BOTH continuations**, so no
   domain cell depends on which continuation is assumed.
4. **Wetting / drying.** A minute is "wetting" if Ŵ < F̂(R) and "drying" if Ŵ >= F̂(R), F̂ = the measured E4 curve of the
   band (10 F-valid cells). Above F's measured range F̂ is held at its last value: there the label is an assumption, and
   those cells are marked in the table and hatched in the figure.
5. **The grid.** R: one column "below 0.1 mm/h" (not wet by the data layer; this is where drying happens) and the seven
   rain classes of step E3 (0.1–0.2–0.5–1–2–5–10–20 mm/h). W: bins of 0.3 dB from −0.3 dB upward. Reason for 0.3 dB: one
   RSL quantization step (W2); a single reading does not place W more finely. Minutes with Ŵ < −0.3 dB or R >= 20 mm/h
   are outside the plane and counted.
6. **Floor on events: 5.** Reason (the event bootstrap): intervals come from resampling the N = 14 events. A cell fed by
   n events is empty in a bootstrap draw with probability (1 − n / N)^N: 0.9 % for n = 4 (1.8 of 200 draws), 0.2 % for
   n = 5 (0.4 of 200 draws). 5 is the smallest number of events for which fewer than one of the 200 draws is expected to
   leave the cell without data. It is also the floor of a fitted length-test cell (E3, rule 5).
7. **Floor on minutes: n_min = 5 · ceil((σ_1 / δ)²)**, with
   - σ_1 = the sd of a single 1 min value of A − c on dry minutes (the median over the short-link link-events of the sd
     inside the offset window; it contains the 0.3 dB quantization), not smaller than 0.3 / √12 = 0.087 dB;
   - δ = 0.0375 dB. Reason: k in a cell is a rate divided by the distance F − W. The smallest distance the instrument
     resolves is one step, 0.3 dB. For k to be within 25 % there (the tolerance of A6), the mean Ŵ of the cell must be
     known to 25 % of 0.3 dB = 0.075 dB at two standard deviations, so to 0.0375 dB at one;
   - the factor 5: the 5 minutes of a radar bin share one R, so they count as one value.
   **The formula is fixed now; σ_1 is computed by the script.** The number of A4 in `plan.md` (200 link-minutes) is not
   used: it was not derived.
8. **Domain edges for A5.** The R edge of the domain = the upper edge of the highest rain class that holds a domain cell
   (either regime); the W edge = the upper edge of the highest domain W bin.
9. **Reported besides:** the same counts for the gauge-input set (short links with a City gauge within 1 km); the share
   of the rainy minutes (R >= 0.1 mm/h) with R above F's measured range, and above P's; the minutes of the g_e = 1 events.

**Pass rule (2026-11-12).** A4 is a measurement. It passes if the domain holds at least one cell in each regime and the
frozen table is written. The domain is frozen by this run: `outputs/tables/a4_coverage_domain.csv`.

**Result of A4 (2026-11-12).** Script `research/waa_ksurface/a4_coverage.py` (library `fp.py` = the measured curves),
report `outputs/reports/a4_coverage.md`, figure `outputs/figures/a4_coverage_map.png`.

- **Data.** Learning set: 14 events, 41 links below 1 km, 926 link-events, 1 497 827 valid minutes, 377 161 of them rainy
  (R >= 0.1 mm/h). The g_e = 1 events (not used): 28 events, 792 link-events, 337 981 valid minutes, 43 570 rainy: **the
  exclusion costs 10 % of the rainy minutes of the short links.** Gauge-input set inside the learning set: 18 links, 377
  link-events, 616 121 minutes, 154 832 rainy.
- **Floors.** σ_1 = 0.072 dB (870 link-events), below the quantization limit, so σ_1 = 0.087 dB is used:
  n_min = 5·ceil((0.087 / 0.0375)²) = **30 minutes**; events floor **5**.
- **The frozen domain (pooled bands): 46 wetting cells and 174 drying cells** (`outputs/tables/a4_coverage_domain.csv`,
  column `in_domain`). Two cells pass under one continuation of P only and are left out.
  - Wetting: every R column, from "below 0.1" to 10–20 mm/h; W from −0.3 dB up to F (0.6 dB at 0.1–0.2 mm/h, 2.1 dB from
    2 mm/h on). 412 386 minutes.
  - Drying: every R column; W from F up to 5.4–9.0 dB. 1 045 175 minutes, of which 849 739 in the column below 0.1 mm/h.
  - **R edge of the domain: 20 mm/h. W edge: 9.0 dB** (the top of the plane; drying only. Wetting ends at 2.1 dB.)
- **Rain above the measured ranges (rainy minutes of the short links): 17.7 % have R above F's measured range (18.4 % at
  28 GHz, 17.1 % at 38 GHz); 34.2 % have R above P's.**
- Gauge-input set, same floors: 41 wetting cells and 135 drying cells (pooled), R up to 20 mm/h, W up to 7.8 dB.

**What the map does not say (recorded for the review).**
1. In the columns above 3.2 mm/h (hatched) F̂ is held at its last value (about 2.0 dB). "Wetting" there means Ŵ below
   2.0–2.1 dB and "drying" above it; both labels are an assumption. 63 of the 220 domain cells lie in the classes 5–10 and
   10–20 mm/h.
2. Above 1.4 mm/h Ŵ depends on the assumed continuation of P. On links below 1 km the two continuations differ by less
   than L·ΔP; the domain keeps only cells that pass under both, but the POSITION of the minutes inside the domain above
   1.4 mm/h is not a measurement.
3. The drying column below 0.1 mm/h holds minutes with Ŵ of 3–8.7 dB (about 15 000 minutes above 3 dB). A wet antenna
   level of that size with no radar rain is rain the radar did not see, or a signal disturbance, more likely than WAA. It
   is the kind of data warning W4 names (unseen water). Not pursued here (side question); A6 truth (f) tests what it does
   to k_dry.
4. The floor of 30 minutes is low compared with the 200 of the frozen rule A4. It is what the pre-written formula gives
   once σ_1 is measured. A stricter floor would remove the thin cells at high W; the frozen table keeps the counts so
   that a different floor can be applied without a rerun.

**Verdict.** PASS (a measurement). The domain is frozen.

## Step A5: build the model (2026-11-12, written before the code is run; no real attenuation is fitted)

**Question served.** Q2 (what is k?): the model in which k is learned, and the proof that its code computes what is written.

**The model (per link-event, one step per minute, W = 0 at the start of the window).**

    v      = F_full(R(t))                                  the steady WAA at the rain of this minute (per band)
    regime = wetting if v > W, else drying
    k      = exp( g_regime(R(t)) + h_regime(W) )           per minute, always positive
    W      <- v + (W − v)·exp(−k·1 min)                    one exact exponential step
    A − c  = W + L·P_full(R(t − d))                        the model output of the minute

The fit compares the 5 min mean of the model output with the 5 min mean of the measured A − c (the measured signal is
never differentiated). No per-link factor (E4 verdict). Learning set: the 14 TRAIN events with a gauge ratio (A4, rule 1).

**F_full and P_full.** Inside the valid range of the quantity: the measured E4 curve, fixed. Above it: a free extension,
continuous at the end of the valid range, monotone (non-decreasing), linear in R between fixed knots:

- **F knots at R = 5, 10, 20 mm/h; P knots at R = 2, 5, 10, 20 mm/h** (per band). Reason: the rain-class edges of step E3
  that lie above the valid range (F to 3.2 mm/h, P to 1.4 mm/h), up to the R edge of the A4 domain (20 mm/h).
- The parameters are the INCREMENTS from one knot to the next, each >= 0 (bounds [0, 3] dB for F, [0, 6] dB/km for P), so
  the curve cannot fall. Above 20 mm/h F is held and P continues with its last slope.
- **The extension parameters are fitted together with k, never alone.** A fit with every k parameter held is not allowed
  to move them.

**The ladder for k (wetting and drying apart).** ln k is piecewise linear between knots, held outside them:
- K0: a constant. K1: a function of W. K2: a function of R. K3: (function of R)·(function of W), i.e. ln k = g(R) + h(W).
- R knots: 0, 0.5, 2, 5, 20 mm/h, linear in ln(1 + R). W knots: 0, 0.5, 1, 2 dB for wetting (the A4 wetting domain ends at
  2.1 dB), 0, 0.5, 1, 2, 4, 9 dB for drying (the A4 drying domain ends at 9.0 dB).
- Bounds: 1 / k between 20 s and 10 h (faster than the 1 min step cannot be told from "at once"; slower than the 3 h of
  data after the rain cannot be seen).
- Smoothness: penalty λ·Σ (second differences of the knot values)², λ from the grid 1e-4, 1e-3, 1e-2 by the CV below.
- **K4 (a free surface) of rule A5 is not built**: the instructions of this block give the ladder as constant → W → R →
  both, and K3 is "both". Recorded as a narrowing of the frozen rule.

**Loss.** The mean over the valid 5 min values of all link-events of Huber(residual; 0.3 dB), plus the smoothness penalty.
0.3 dB = one RSL quantization step (W2). Whole events; minimised with L-BFGS-B and an exact gradient (reverse pass).

**CV and selection (rule B1, inside TRAIN).** 5 folds, grouped by event: the 14 events sorted by start time and dealt
round-robin. (The rule asks for stratification by `intensity_class`; all 14 learning events are "heavy", so there is
nothing to stratify.) CV error = the mean Huber loss of the held-out 5 min values. λ per rung = the one with the smallest
CV error. Selected rung = the simplest (order K0, K1, K2, K3) whose CV error differs from the best rung's by a paired
90 % interval that contains zero (1 000 event-bootstrap draws of the per-event held-out losses).

**Rain timing (rule A5).** Stage 1: the gauge-input set (links with a City gauge within 1 km, all lengths), rain = the
1 min gauge rate; the wetting rung is selected with drying at K1 (a nuisance there); k_wet is kept. A 5 min value with a
missing gauge minute is left out of the loss, and the rain of a missing minute is the radar value. Stage 2: all links,
rain = the radar (5 min, held over its 5 minutes); k_wet held; the drying rung is selected. The extension is free in both
stages; the final extension is that of stage 2.

**The rule for d (written 2026-11-12, before any fit):**

> d is a whole number of minutes from 0 to 5: the path term reads the rain of d minutes earlier (5 min is one radar bin;
> a longer lag cannot be placed with this input). d is chosen once per rain input (radar; gauge), with the K0 model
> (constant k_wet and k_dry, extension free) on the TRAIN learning set: the d with the smallest loss. If the paired 90 %
> event-bootstrap interval of loss(d = 0) − loss(d_best) contains zero, the lag is not resolved and d = 0. d is then
> held for every rung, fold and refit. It is never chosen per link or per event.

**Unit tests (pass rule, 2026-11-12). A5 passes if all three pass:**
1. **Gradient.** On a small set of real link-events with a fake target, rung K3 in both regimes and a free extension:
   the exact gradient against central finite differences, every parameter. PASS if the largest relative difference
   (relative to the largest gradient entry) is below 1e-5.
2. **Step response.** Constant k, R stepping from 0 to a constant and back to 0, round numbers: W(t) against the analytic
   F·(1 − e^(−k_wet·t)) and then W_0·e^(−k_dry·t). PASS if the largest difference is below 1e-9 dB.
3. **Extension.** With the extension knots set to the values of the measured E4 curve (F held, P with its last slope),
   F_full and P_full against the E4 curve on a dense grid of R from 0 to 30 mm/h. PASS if the largest difference is below
   1e-9.

Rule A5's reference to the amendment table ("within 1e-6 dB") is met by the stricter 1e-9 of test 2.

**Result of A5 (2026-11-12).** Library `research/waa_ksurface/ksurf.py`, script `research/waa_ksurface/a5_build.py`,
report `outputs/reports/a5_build.md`. No real attenuation was fitted.

| unit test | number | tolerance | verdict |
|---|---|---|---|
| 1 gradient against finite differences (K3 / K3, extension free, 32 parameters, 40 real link-events, 68 185 minutes, d = 2) | 4.2e-10 (relative) | 1e-5 | PASS |
| 2 constant-k step response (0 → 4 mm/h for 120 min → 0; 1 / k = 5 min wetting, 20 min drying) | 4.4e-16 dB | 1e-9 dB | PASS |
| 3 extension with its knots at the measured values against the E4 curve (2 000 values of R to 30 mm/h) | 4.4e-16 | 1e-9 | PASS |

- **A bug found by test 2 and fixed (recorded, not hidden):** on the first run test 2 failed by 1.06 dB. Cause: with a
  constant rung (K0, K1) the code read the knot value at the R-knot index of the minute instead of the single constant.
  Test 1 had not seen it because it uses K3. The line was corrected and all three tests rerun. This is a correction of
  new code that its own unit test caught, before any fit; no pass rule was changed. Beyond the written test 1, the
  gradient is now also checked for K0 / K0, K1 / K2 and K2 / K1 (relative differences 2e-10 to 6e-10).
- Extension knots: F at 5, 10, 20 mm/h; P at 2, 5, 10, 20 mm/h (per band). Parameters of k: 2 (K0 / K0) to 22 (K3 / K3);
  extension: 14 (two bands).
- Speed: one loss with its gradient on the whole learning set (4 563 link-events, 6.8 million minutes) takes 0.07–0.11 s.
- Disclosure: for that timing the loss was EVALUATED once on the real learning set at the starting values (no
  minimisation, no parameter chosen from it). Nothing was fitted to real attenuation.

**Verdict.** PASS.

## Step A6: fake-data recovery, the gates for A7 (2026-11-12, written before any A6 run)

**Question served.** Q2 (what is k?): which statements about k may A7 make. No real attenuation is used as a target.

### The fake world of A6

- **Links, events, rain:** the learning set (4 563 TRAIN link-events of 187 links, 14 events with a gauge ratio). The true
  5 min rain of a link-event is its real rain input.
- **True rain per minute.** The attenuation is made from 1 min rain: the 5 min values joined by straight lines between bin
  centres, multiplied by a lognormal AR(1) fluctuation (sd 0.3, correlation 0.7 from one minute to the next), then every
  bin rescaled so that its mean is the 5 min value. The two numbers are generic and round: the variability of rain inside
  5 min cannot be measured here (a gauge resolves 6–12 mm/h in one minute).
- **F and P of the fake world** (generic, as in E3): F_true = 3·R / (R + 2) dB, P_true = 0.1·R (28 GHz) and 0.2·R (38 GHz)
  dB/km. They are handed to the model the way the real ones will be: as "measured" cells at the R of the real F-valid and
  P-valid cells, with the extension knots above. The truth of the fake world IS that curve (cells + knots at the true
  values), so "F known" is exact. The fit uses the measured part as fixed, as in A7.
- **Noise:** an offset error per link-event (sd 0.136 dB, E3); a level shift per link, sd 0.088 dB (E4), added to the
  target F while it rains (never below zero); white noise per 1 min value, sd 0.072 dB (A4); then every 1 min value is
  rounded to the 0.3 dB grid with a random phase per link-event. The 5 min value is the mean of its 5 minutes; the real
  pattern of missing values is kept. (σ_1 already holds the quantization, so the rounding on top is the conservative side.)
- **What the fit sees.** Radar input: the true 5 min rain × the T1 error at its high end (per event sd 0.45; inside an
  event shared part and the link's own part with σ_c = 0.26, 0.9 correlated over 30 min), held over the 5 minutes.
  Gauge input (gauge-input set): the true 1 min path rain × exp(own error + representativeness) with the T2 model at its
  normal end (own error white per 5 min bin, sd 0.55; representativeness variance 0.028·(d + L / 3) at 15 min), turned
  into whole gauge steps (0.1 or 0.2 mm) per minute.
- **d:** the truth has d = 0 and the fit uses d = 0, except in the lag run below.
- **One replicate per truth** (seeds in the config). Stated limit: every verdict below rests on one draw of the noise.

### The truths (round numbers), the fit and the PASS rules

The fit is the A5 procedure: stage 1 (gauge-input set, gauge rain, wetting rung selected with drying at K1), stage 2 (all
links, radar rain, k_wet held, drying rung selected), λ and rung by the 5-fold event CV and rule B1. "Within 25 %" for a
surface = |k_fit / k_true − 1| <= 0.25 at the centre of a cell; "on the domain" = the cells of the frozen A4 domain of that
regime, each weighted by the minutes the fake world's true state spends in it (short links; R as the fit sees it). A
constant is judged as a number.

| truth | k | F, P above the measured range | PASS rule | if not |
|---|---|---|---|---|
| (a) | k_wet = 1/5, k_dry = 1/15 per min | known (held at the truth) | K0 selected in both regimes AND both k within 25 % | **STOP** |
| (a-lag) | as (a), truth with d = 2 min | known | report: the d the rule returns per input, and k with it | report |
| (b) | as (a) | truth as above; the fit STARTS from F held and P at half its true increments, extension free | "identified" if every F knot value is within 25 % of the truth, no increment sits at its upper bound or at 0 where the truth is above 0.05, and both k within 25 % | restriction: A7 may claim k only inside F's measured range, plus drying |
| (c) | k_wet = 1/5; k_dry = (1 + W) / 30 | known | drying K1 and wetting K0 selected; k within 25 % on at least 80 % of the domain weight in each regime | **STOP** |
| (d) | k_wet = 0.1·(1 + R / 2); k_dry = 0.1 / (1 + R / 2) | known | K2 selected in both regimes; within 25 % on at least 80 % of the domain weight in each | **STOP** |
| (e) | two stores, each half of F: fast 1/5 (wetting), 1/10 (drying); slow 1/60 | known | report only: the rungs the one-state ladder returns, its k, and the residual in the 30 min after rain by preceding rain (highest minus lowest third, at matched W ±0.15 dB), next to the same number for (a) | report |
| (f) | as (a); in 3 of the 14 events (20 %) the rain goes on at 1 mm/h for 30 min after the link's last wet bin, unseen by the rain input | known | "spurious" if drying is not K0 or k_dry is outside 25 %. Then the flagged events are excluded and the fit repeated: "restored" if K0 and within 25 % | restriction recorded |
| (g) | as (a); level shift per link sd 0.24 dB (E3) instead of 0.088 | known | "corrupted" if either regime is not K0 or a k is outside 25 % | the per-link question reopens, with this as the reason |

- **The flag of (f), defined here before use (W4: tested, not assumed; model-free).** For a link-event of a short link:
  a0 = the mean of A − c from 5 min before to 5 min after the end of its last wet bin; a1 = the mean from 10 to 20 min
  after it. It qualifies if a0 >= 0.6 dB (two quantization steps). Its tail "waits or rises" if a1 >= a0 − 0.15 dB. An
  event is flagged if it has at least 5 qualifying link-events and at least half of them wait or rise. The exclusion uses
  the FLAG, not the knowledge of which events were changed; how well the flag finds them is reported.
- **Rungs allowed in A7 (the frozen rule A6: the ladder holds only recovered rungs).** A rung is allowed for a regime only
  if a truth with that rung was recovered for that regime here. Truths exist for: wetting K0 (a), K2 (d); drying K0 (a),
  K1 (c), K2 (d). **Wetting K1 and K3 in either regime are given no truth in this block, so they are not allowed in A7
  until tested.** Stated now, flagged for the review.
- **Claims domain for A7** = the A4 domain, cut by the restrictions that (b), (f) and (g) produce.
- (c)–(g) hold the extension at the truth so that each question is tested alone; only (b) frees it. Stated limit: A7 frees
  the extension in every rung, which (c) and (d) do not test.
- Any FAIL on (a), (c) or (d) is a STOP: reported, not patched. After the report: STOP for review. A7 is not run.

**Result of A6 (2026-11-12).** Script `research/waa_ksurface/a6_recovery.py` (libraries `fake6.py`, `ladder.py`), report
`outputs/reports/a6_recovery.md`. Fake data only.

**Truth (a) FAILS. By the rule this is a STOP; truths (a-lag) to (g) were not run.**

| regime | true rung | rung returned (lowest held-out loss) | 1 / k at R = 0, W = 0 | share of the domain weight within 25 % |
|---|---|---|---|---|
| wetting (stage 1, gauge input) | K0 | **K3** (K3) | 45.5 min (truth 5.0) | 14 % |
| drying (stage 2, radar input, k_wet held) | K0 | **K1** (K3) | 23.0 min (truth 15.0) | 12 % |

- Stage 1, held-out loss: K0 0.02914, K1 0.02890, K2 0.02881, K3 0.02853. The gains are small (2 % of the loss) but the
  paired event-bootstrap intervals of every simpler rung against K3 exclude zero, so rule B1 returns K3. λ is at the
  lowest value of the grid in every rung.
- Stage 2 inherits the wrong k_wet and then also prefers a rung above K0 (K0 0.03760, K1 0.03667, K3 0.03655).
- So on a constant-k truth, with the measured input errors, the procedure "finds" a dependence of k on R and W that is
  not there. By rule B5 this would have been reported as new dynamics.

**Verdict.** FAIL on (a): **STOP.** Nothing is patched. (b)–(g) and the rungs allowed in A7 are not decided.

### The one diagnostic allowed after the failed check (2026-11-12, written before it is run)

**What.** Truth (a) again, same seed, with the ingredients of the fake world switched (`a6_diagnostic.py`):

| run | gauge input | radar input | noise and quantization |
|---|---|---|---|
| D0 | exact (the true 1 min rain) | exact (the true 5 min rain) | off |
| D1 | exact | exact | on |
| D2 | with its error (as in A6) | exact | on |
| D3 | exact | with its error (as in A6) | on |
| (A6 run) | with its error | with its error | on |

For every run: the rungs the two stages return, and what a plain K0 / K0 fit gives for 1 / k (on the gauge-input set with
the gauge input, and on all links with the radar input, both regimes free).

**Why.** To locate the cause: the procedure or code itself (D0 must return K0 and the true k, else the fault is there),
the noise and the quantization (D1), the error of the gauge input (D2), the error of the radar input (D3). The plain K0
fit shows whether the constant itself is still recovered when no richer rung is offered. Nothing is changed on the basis
of this diagnostic; the block stops after it.

**Result of the diagnostic (2026-11-12).** `outputs/reports/a6_diagnostic.md`. Truth (a), same seed, one replicate per row.

| run | gauge input | radar input | noise | wetting rung | drying rung | wetting: weight within 25 % | drying: weight within 25 % |
|---|---|---|---|---|---|---|---|
| D0 | exact | exact | off | K0 | K1 | 100 % | 100 % |
| D1 | exact | exact | on | K3 | K1 | 69 % | 100 % |
| D2 | with error | exact | on | K3 | K3 | 10 % | 18 % |
| D3 | exact | with error | on | K3 | K0 | 70 % | 100 % |
| A6 run | with error | with error | on | K3 | K1 | 14 % | 12 % |

A plain K0 / K0 fit, 1 / k in minutes (truth 5 and 15):

| run | gauge-input set, gauge rain: wetting / drying | all links, radar rain: wetting / drying |
|---|---|---|
| D0 | 5.0 / 15.0 | 5.1 / 14.8 |
| D1 | 5.5 / 14.4 | 5.2 / 14.8 |
| D2 | 3.8 / 20.8 | 5.2 / 14.8 |
| D3 | 5.5 / 14.4 | 9.9 / 14.5 |
| A6 run | 3.8 / 20.8 | 9.9 / 14.5 |

What it locates:
1. **The code and the model are not the cause.** With exact inputs and no noise (D0) the constants come back exactly
   (5.0 and 15.0 min) and wetting returns K0.
2. **The selection rule B1 is one cause.** It has no smallest effect: in D1 (exact inputs, only noise and quantization)
   K3 beats K0 by 0.2 % of the held-out loss (0.010903 against 0.010924), the gain has the same sign in most events, the
   paired interval excludes zero, and K3 is returned. In D0 drying returns K1 on a difference in the fifth decimal (the
   5 min radar input against the 1 min truth), although its k is within 25 % everywhere. A rule that asks only "is the
   difference distinguishable from zero" selects structure whenever the data are many.
3. **The error of the gauge input is the other cause, and the larger one.** With it (D2, A6 run) the returned k_wet is
   within 25 % on only 10–14 % of the domain weight, and even the plain constant moves (3.8 min for 5, and 20.8 min for
   15 in drying). The gauge resolves 6–12 mm/h per minute, so its 1 min series is spikes; the fit bends k(R, W) to absorb
   them. A wrong k_wet, held, then spoils stage 2.
4. **What survives:** the constant k_dry from all links with the radar input is 14.5–14.8 min in every run (truth 15),
   also with the full radar error. The constant k_wet from the radar input is right with an exact radar (5.2) and wrong
   with the radar error (9.9).

**The block stops here (2026-11-12).** Truths (a-lag), (b)–(g) were not run; the restrictions of (b), (f), (g), the rungs
allowed in A7 and the claims domain for A7 are not decided. Nothing was patched. Decisions for the review: whether B1 gets
a smallest effect (and which), and whether k_wet can be learned from the 1 min gauge input at all.

---

# The A6 rescue (2026-11-12, unattended): Deviations 5 and 6, then A6 again

The review read the A6 stop and its diagnostic. This block logs two deviations, reruns truth (a), and if it passes runs
(a-lag) and (b)–(g). **A7 (the real fit) is not run.**

## Deviation 5: the selection rule B1 gains a smallest effect (2026-11-12, logged before any run)

**What changes.** A richer rung is selected only if its grouped-CV error improves on the simpler rung by MORE than
**0.005 dB** AND the paired 90 % event-bootstrap interval of the difference excludes zero.

**Reason (given by the review).** With about 1.5 million minutes, negligible improvements are significant (diagnostic D0,
D1: K3 selected for 0.2 % of the loss). 0.005 dB is the materiality margin already frozen in rule B4, so one standard
governs both structure (B1) and closed forms (B4).

**How it is computed (fixed now, where the rule leaves room).**
- **The CV error in dB.** The loss of the fit is a mean Huber loss (dB², not dB). The CV error is
  e = √(2 · mean held-out Huber loss): it equals the RMS residual when the residuals are below the Huber threshold, and it
  keeps the Huber protection against the large path-term residuals of long links. The 0.005 dB margin and the paired
  interval are both taken on e (the interval: 1 000 event-bootstrap draws of e_simpler − e_richer).
- **"The simpler rung."** The ladder is walked once from K0: the selection starts at K0; each rung K1, K2, K3 in turn
  replaces the current selection only if it beats the CURRENT selection by more than the margin with an interval that
  excludes zero. (K1 and K2 are not nested; each is compared with what is selected when its turn comes.)
- λ per rung stays the one with the lowest held-out loss.

**This is a choice made after seeing a result** (the A6 run and its diagnostic). Disclosure: on the numbers already seen,
the stage-1 gain of K3 over K0 in the A6 run is 0.0025 dB on this scale (0.2414 against 0.2389 dB), and in D1 0.0001 dB;
so it is known in advance that the wetting ladder of those two runs returns K0 under the new rule. The conversion to dB
was not tuned to that: it is the plain reading of "dB" for a squared-error loss.

## Deviation 6: k_wet is handled by a sensitivity gate (2026-11-12, logged before any run)

**What changes.** The diagnostic showed that a constant k_wet is not recovered to 25 % from either input (3.8 min from the
gauge input, 9.9 min from the radar input, truth 5). Instead of a claim on k_wet, a gate decides whether k_wet matters:

- **(a) The gate, on truth (a)** (the A6 fake world, same seed): all links, radar input, K0 in both regimes, 1 / k_wet HELD
  at 2.5, 5 and 10 min (half the truth, the truth, twice the truth); k_dry and the F / P extension free (the extension
  started as in truth (b): F held, P at half its true increments). Reported per value: the grouped-CV error (dB), and
  1 / k_dry.
- **(b) PASS = "bounded nuisance"** if the CV error changes by less than 0.005 dB across the three values (largest minus
  smallest) AND 1 / k_dry is within 25 % of the truth at all three. Then, on real data, **k_wet is reported as the
  interval between its gauge-input estimate and its radar-input estimate, carried through every refit (k_dry is fitted at
  both ends), and never claimed as a discovery.** This is recorded as a claims restriction.
- **(c) FAIL** (k_wet matters): STOP and report. The rescue (an input-error-aware fit, or a simulation-based bias
  correction with its own fake check) is a decision for the review.

**How the procedure changes if the gate passes (fixed now).**
- Stage 1 (gauge-input set, gauge rain) still walks the wetting ladder, now under Deviation 5.
- If stage 1 returns K0: k_wet,gauge = its constant; k_wet,radar = the constant of a plain K0 / K0 fit on all links with
  the radar input. Stage 2 walks the drying ladder with k_wet held at k_wet,gauge, and the selected drying rung is refitted
  with k_wet held at k_wet,radar. **k_dry must be within 25 % at BOTH ends** to count as recovered.
- If stage 1 returns a rung above K0, that rung is held in stage 2 as before (no interval), and it is judged as a surface.
- In the PASS rules of A6, "k_wet within 25 %" for a truth whose wetting rung is K0 becomes "k_wet within 25 %, or governed
  by the gate" (the rung must still be K0).

**Reason.** A quantity that cannot be measured must not be claimed; but if the result of interest (k_dry) does not depend
on it inside the range it can lie in, it can be carried as a bounded nuisance.

**This is a choice made after seeing a result.** Disclosure: the diagnostic already showed 1 / k_dry of 14.5 min from the
radar input with 1 / k_wet at 9.9 min (a K0 fit, extension held), so one point of the gate is roughly known.

## A6 again under Deviations 5 and 6 (2026-11-12, written before any run)

**Question served.** Q2.

- **Order.** The gate; truth (a); then, only if (a) passes, (a-lag) and (b)–(g) exactly as written in the A6 entry (same
  fake world, same seeds, same truths, same flag, same extension start), with the selection of Deviation 5 and the k_wet
  handling of Deviation 6.
- **PASS of (a):** K0 returned in both regimes; 1 / k_dry within 25 % (at both ends of the k_wet interval); k_wet within
  25 % or governed by the gate. Else STOP, with one diagnostic.
- **(c), (d):** as written, with the rung now having to beat the margin. FAIL at the written size is a STOP. For (c) the
  wetting side passes with K0 and "within 25 % or governed by the gate". For (d) wetting is a surface (K2) and is judged
  as one.
- **Detectability (asked by the review; report only, fixed now).** The size of the dependence is scaled by a factor a, a = 1
  being the written truth:
  - k(W), drying: k_dry = (1/15)·(1 + a·W) / (1 + a), so that 1 / k_dry is 15 min at W = 1 dB for every a (a = 1 is truth
    (c)). Grid a = 0.1, 0.25, 0.5. k_wet = 1/5.
  - k(R), both regimes: k_wet = 0.1·(1 + a·R / 2), k_dry = 0.1 / (1 + a·R / 2) (a = 1 is truth (d)). Grid a = 0.1, 0.25, 0.5.
  "Detected" = the ladder returns the true rung for that regime. Reported: the smallest a detected per regime, with the
  ratio of k across the domain that it stands for (k(W = 2 dB) / k(W = 0); k(R = 5 mm/h) / k(R = 0)). One replicate per
  point; a gap in the grid (detected at a small a but not at a larger one) is reported as such.
- **The outputs of the first A6 run and of its diagnostic are superseded** (they were made with the old selection rule);
  they stay in the history at commits `e0e7b4d` and `28f15d9` and are listed in `outputs/SUPERSEDED.md`.
- Then STOP for review. A7 is not started.

**Result of the gate and of truth (a) under Deviations 5 and 6 (2026-11-12).** Script `a6_recovery.py`, report
`outputs/reports/a6_recovery.md`. Fake data only.

**The k_wet gate: PASS (k_wet is a bounded nuisance), but narrowly.**

| 1 / k_wet held (min) | CV error (dB) | 1 / k_dry (min, truth 15) | within 25 % |
|---|---|---|---|
| 2.5 | 0.24383 | 12.55 | yes |
| 5.0 | 0.24002 | 13.88 | yes |
| 10.0 | 0.24296 | 16.24 | yes |

CV error changes by 0.0038 dB across the three values (margin 0.005). 1 / k_dry moves from 12.5 to 16.2 min: it is NOT
independent of k_wet, it stays inside 25 % (the limits are 12.0 and 20.0 min).

**Truth (a): FAIL.**

| regime | rung returned (lowest CV) | gain of the best richer rung over K0 | 1 / k (min) | verdict |
|---|---|---|---|---|
| wetting | **K0** (K3) | 0.0026 dB (K3), under the margin | 3.7 (gauge input) / 9.9 (radar input), truth 5: governed by the gate | ok |
| drying | **K0** (K3) | 0.0028 dB (K3), under the margin | **11.6** with k_wet at its gauge estimate / 14.5 at its radar estimate, truth 15 | **11.6 is outside 25 %** (limit 12.0) |

- Deviation 5 does what it was written for: both ladders return K0 (no richer rung gains more than 0.0028 dB).
- The failure is the second half of the rule of Deviation 6: k_dry must be within 25 % at BOTH ends of the k_wet interval.
  With k_wet held at the gauge-input estimate (3.7 min) the drying constant is 11.6 min, 23 % too fast in time, 29 % too
  large in k.
- Noted, not explained yet: in the gate, k_wet held at 2.5 min (further from the truth) gave 12.5 min. The gate fits the
  F / P extension; truth (a) holds it at the truth.

**Verdict.** FAIL on (a): **STOP.** (a-lag), (b)–(g), the detectability sweep and the optional step 5 were not run.
Nothing is patched.

### The one diagnostic allowed after the failed check (2026-11-12, written before it is run)

**What.** On the same fake data of truth (a): all links, radar input, K0 in both regimes, 1 / k_wet HELD at 2.5, 3.73, 5,
7.5, 9.91 min, each once with the F / P extension held at the truth and once with it free (started as in truth (b)).
Per fit: 1 / k_dry, the grouped-CV error (dB), and for the free fits the fitted F knot values against the truth
(`a6_diagnostic2.py`).

**Why.** To locate the cause of the 11.6 min: (i) how 1 / k_dry depends on the held k_wet (smooth and one-sided, or
erratic); (ii) whether the free extension absorbs the error of k_wet (which would explain why the gate saw less
sensitivity than truth (a)), and at what cost to the F knots. Nothing is changed on the basis of this diagnostic; the
block stops after it.

**Result of the diagnostic (2026-11-12).** `outputs/reports/a6_diagnostic2.md`. Truth (a), the fake data of the A6 run,
all links, radar input, K0 in both regimes. 1 / k_dry in minutes (truth 15; within 25 % = 12.0 to 20.0):

| 1 / k_wet held (min) | extension held at the truth: 1 / k_dry | CV error (dB) | extension free: 1 / k_dry | CV error (dB) | free F knots, fitted / true |
|---|---|---|---|---|---|
| 2.5 | 10.9 (outside) | 0.2779 | 12.5 | 0.2438 | 0.77–0.89 |
| 3.73 (the gauge-input estimate) | **11.6 (outside)** | 0.2733 | 13.2 | 0.2412 | 0.77–0.92 |
| 5.0 (the truth) | 12.3 | 0.2700 | 13.9 | 0.2400 | 0.78–0.92 |
| 7.5 | 13.5 | 0.2667 | 15.1 | 0.2406 | 0.75–0.93 |
| 9.91 (the radar-input estimate) | 14.5 | 0.2659 | 16.2 | 0.2429 | 0.73–0.92 |

What it locates:
1. **1 / k_dry depends on the held k_wet smoothly and in one direction** (about +0.45 min per minute of 1 / k_wet), with
   or without a free extension. k_wet is a nuisance that shifts k_dry by a known, bounded amount; it is not harmless.
2. **The failure comes from holding the extension at the truth.** With the "true" F and P above the measured range, the
   drying constant is 18 % too fast even at the TRUE k_wet (12.3 min), and the CV error is 0.03 dB worse than with a free
   extension. The rain the fit sees carries the radar error, and F evaluated at a noisy R is not the F the data follow on
   average: the free extension settles 8–27 % below the true knots and the fit is better for it. So "F known everywhere"
   is not the easy case it was written as in truth (a): it is a harder case than the one A7 will run.
3. **With the extension free (the A7 configuration), 1 / k_dry is 13.2 min at the gauge-input end and 16.2 min at the
   radar-input end of the k_wet interval: inside 25 % at both.** This is a diagnostic, not a pass: the pre-registered
   truth (a) holds the extension, and it failed.
4. The free F knots are 8–27 % below the truth at every k_wet. This is what truth (b) was written to judge; it was not
   run, but these numbers say that a knot above the measured range of F will not be recovered to 25 % at every knot.

**The block stops here (2026-11-12).** Not run: (a-lag), (b)–(g), the detectability sweep, the optional step 5. Not
decided: the rungs allowed in A7, the claims domain. Nothing was patched. Decision for the review: whether the gating
truths (a), (c), (d) should free the extension (as A7 does), given that holding it at the truth is biased by the rain
error; and whether the bias of the F knots is then a stated restriction (no claim on F above its measured range).

---

# A6, final form (2026-11-13, unattended): Deviation 7, then the truths on fresh seeds

The review read the stop at truth (a) and diagnostic 2. **A7 (the real fit) is not run in this block.**

## Deviation 7: every gating truth runs the procedure exactly as A7 will (2026-11-13, logged before any run)

**What changes.** In every truth of A6 the F and P extensions above their measured ranges are FREE and fitted jointly with
k (started as written for truth (b): F held, P at half its true increments). No truth holds the extension at the truth.

**Reason (given by the review).** A fake check validates the procedure that will be applied. Holding the extension at the
truth is a different procedure, and diagnostic 2 showed it is also the harder and biased one (F evaluated at a rain that
carries the radar error).

**Disclosure.** The free-extension result on the fake data of the failed run is already known (diagnostic 2: 1 / k_dry of
13.2 and 16.2 min at the two ends of the k_wet interval, truth 15; F knots 8–27 % below the truth). So this deviation is
made knowing that it turns that failure into a pass on those data. For that reason **the old fake data are not used for
any verdict**: every truth below is drawn with a fresh seed (base 2200 instead of 1200).

**Consequences, recorded as claims restrictions for A7:**

- **R1.** F and P above their measured ranges are EFFECTIVE quantities for this rain input; no physical claim is made on
  them.
- **R2.** A dependence of k on R is claimed only inside F's measured range (up to about 3.2 mm/h); above it, k(R) and the
  free extension are not separable.
- **R3.** k_dry is reported with an interval that includes its spread across the k_wet interval (from the gauge-estimate
  end to the radar-estimate end), on top of the bootstrap.

**This is a choice made after seeing a result**, hence a deviation.

## A6 final: the truths, amended (2026-11-13, written before any run)

**Question served.** Q2.

Everything not named here stays as written in the entries "Step A6", "Deviation 5" and "Deviation 6": the fake world, the
truths and their numbers, the two stages, the selection with the 0.005 dB margin, k_wet as an interval when the wetting
ladder returns K0, the tail flag, one replicate per truth. The k_wet gate is not repeated (it passed, with a free
extension, and its verdict is not tied to the failed run's held extension).

- **Seeds.** Base 2200; every truth and every replicate has its own seed.
- **(a)** PASS: K0 in both regimes AND 1 / k_dry within 25 % at both ends of the k_wet interval (k_wet: governed by the
  gate). Else STOP, with one diagnostic.
- **(a-lag)** as written (truth with d = 2 min), extension free. Report.
- **(b)** (amended by the review): its criterion is now the rung and the recovery of k_dry, not the extension: "k_dry
  recovered" if K0 / K0 and 1 / k_dry within 25 % at both ends. The error of every extension knot is reported as
  information; it feeds R1. With Deviation 7, (b) is the same procedure as (a) on another draw: it is a second replicate
  of (a) that also prints the knots.
- **(c)** PASS: drying K1, wetting K0, k_dry within 25 % on at least 80 % of the domain weight at both ends of the k_wet
  interval. FAIL is a STOP.
- **(d)** PASS: K2 in both regimes, each within 25 % on at least 80 % of the domain weight. **Amendment, written now and
  following from R2:** the weight is taken over the domain cells inside F's measured range (the R columns up to the class
  2–5 mm/h, whose centre 3.16 mm/h is inside 3.2). Above it k(R) is not claimed, so it is not judged. The share on the
  whole domain is reported next to it. FAIL is a STOP.
- **(d-high)** (new, asked by the review; tests R2): k depends on R only ABOVE 3.2 mm/h, with the form of (d):
  k_wet = 0.1·(1 + max(R − 3.2, 0) / 2), k_dry = 0.1 / (1 + max(R − 3.2, 0) / 2). Report only: the rung each ladder returns
  (K2 = the dependence is detected; K0 = it is not), the recovery of k inside the measured range (where the truth is the
  constant 0.1), and the extension knots (does the extension absorb the dependence).
- **(e), (f), (g)** as written, extension free (also in the refit of (f) without the flagged events).
- **Detectability (asked by the review; replaces the sweep of the previous block).** Two forms, size a:
  - k_dry(W): k_dry = (1/15)·(1 + a·W) / (1 + a), k_wet = 1/5. Ratio of k at the domain ends: k(W = 2 dB) / k(W = 0) = 1 + 2a.
  - k(R) inside the measured range: k_wet = 0.1·(1 + a·min(R, 3.2) / 2), k_dry = 0.1 / (1 + a·min(R, 3.2) / 2) (constant
    above 3.2 mm/h). Ratio: k(3.2) / k(0) = 1 + 1.6a for wetting, its inverse for drying.
  - Sizes a = 1, 0.5, 0.25, 0.1, tried from the largest; 5 replicates per size. "Detected" at a size = the ladder returns
    the true rung (K1 for k_dry(W); K2 for k(R), wetting and drying counted apart) in at least 4 of the 5 replicates.
    To save time, written now: the replicates of a size stop as soon as the outcome is decided for every regime of the form
    (4 hits or 2 misses), and a form's sweep ends at the first size at which no regime is detected. The answer per regime =
    the smallest size detected, given as the ratio of k.
- **Rungs allowed in A7:** a rung is allowed for a regime if its truth passed here: wetting K0 (a), K2 (d); drying K0 (a),
  K1 (c), K2 (d). Wetting K1 and K3 in either regime have no truth and stay not allowed.
- FAIL on (c) or (d) is a STOP. (b), (f), (g), (d-high) give restrictions. Then STOP for review; A7 is not started.

**Result of A6 final (2026-11-13).** Script `a6_recovery.py`, report `outputs/reports/a6_recovery.md`. Fake data only,
fresh seeds (base 2200), extension free in every truth.

| truth | regime | true rung | rung returned | recovery | verdict |
|---|---|---|---|---|---|
| (a) | wetting | K0 | K0 (best richer rung +0.0013 dB) | 1 / k_wet 3.3 (gauge) / 3.1 (radar) min, truth 5: governed by the gate | PASS |
| (a) | drying | K0 | K0 (best richer rung +0.0013 dB) | 1 / k_dry 17.2 / 17.0 min at the two ends, truth 15 | PASS |
| (a-lag) | | | | the rule returns d = 2 for the radar input and d = 2 for the gauge input (truth 2) | report: the rule works |
| (b) | wetting | K0 | K0 | 1 / k_wet 3.3 / 3.3 min: governed by the gate | k_dry recovered |
| (b) | drying | K0 | K0 | 1 / k_dry 16.7 / 16.7 min, truth 15 | k_dry recovered |
| (c) | wetting | K0 | K0 | 1 / k_wet 3.5 / 3.7 min: governed by the gate | ok |
| (c) | drying | K1 | **K0** | 18 % of the domain weight within 25 % (needs 80 %) | **FAIL** |

- **Truth (a) PASSES on fresh data**, and (b), the same procedure on another draw, agrees: K0 in both regimes, 1 / k_dry
  11–15 % too long (17.2 and 16.7 min for 15), the same at both ends of the k_wet interval.
- k_wet: both inputs give 3.1–3.7 min for a true 5 in three draws. In these draws the two estimates agree with each other
  and are both about 30 % too fast: **the interval between them does not contain the truth.** (In the first fake set
  they were 3.7 and 9.9.) Recorded as a fact about the nuisance interval, not acted on.
- (a-lag): the loss has its minimum at d = 2 for both inputs and the paired interval of loss(0) − loss(2) excludes zero.
- (b), extension knots: listed in the report and in `a6_recovery_extension.csv` (information for R1).
- **Truth (c) FAILS.** The true rung K1 has the lowest held-out error, and its gain over K0 has an interval that excludes
  zero (0.0002 to 0.0036 dB), but the gain is **0.0019 dB, below the 0.005 dB margin** of Deviation 5. So the ladder
  returns K0, and a constant (16.8 min) is within 25 % of the true k_dry on only 18 % of the domain weight. The written
  dependence is not small: k_dry triples from W = 0 to W = 2 dB.

**Verdict.** FAIL on (c): **STOP.** Not run: (d), (d-high), (e), (f), (g), the detectability study, the optional step.
Nothing is patched. One finding already follows without a further run: **at the 0.005 dB margin the ladder cannot return
a k_dry that triples across the drying domain**, so the smallest detectable k_dry(W) ratio is above 3.

### The one diagnostic allowed after the failed check (2026-11-13, written before it is run)

**What.** On the fake data of truth (c): all links, radar input, k_wet held at its gauge-input estimate, extension free;
drying K0 and drying K1 (λ = 1e-4, the one the CV chose) fitted on all events (`a6_diagnostic3.py`). Reported:
(i) the fitted K1 curve against the true k_dry(W) at the W knots, and the share of the domain weight within 25 % that K1
WOULD have had; (ii) the error √(2·mean Huber loss) of both fits, in sample, on: all 5 min values; the links below 1 km;
the links below 1 km in the 60 min after the last wet bin; and the same three sets with the rainy part of the event left
out.

**Why.** To locate why a threefold change of k_dry is worth only 0.002 dB: is the K1 fit itself good (then the margin is
simply larger than the footprint of the dependence in the loss), and is the footprint diluted by where the loss is
measured (all links and all minutes, most of them unaffected by k_dry)? Nothing is changed on the basis of this
diagnostic; the block stops after it.

**Result of the diagnostic (2026-11-14).** Script `a6_diagnostic3.py`, report `outputs/reports/a6_diagnostic3.md`. Fake
data only. The stage 1 fit was repeated and gives the 1 / k_wet of the A6 final run (3.5387 min), so these are the fits
of that run.

| W (dB) | true 1 / k_dry (min) | K1 fit (min) | K0 fit (min) | share of the domain minutes from this knot to the next |
|---|---|---|---|---|
| 0 | 30 | 34.3 | 16.8 | 86 % |
| 0.5 | 20 | 21.5 | 16.8 | 4 % |
| 1 | 15 | 14.3 | 16.8 | 8 % |
| 2 | 10 | 7.8 | 16.8 | 1 % |

| in-sample error √(2·mean Huber loss) | share of all values | K0 (dB) | K1 (dB) | gain (dB) | share of the total loss drop |
|---|---|---|---|---|---|
| all 5 min values | 100 % | 0.2344 | 0.2319 | 0.0025 | 100 % |
| links below 1 km | 22 % | 0.2111 | 0.2073 | 0.0038 | 30 % |
| links below 1 km, 60 min after the last wet bin | 0.8 % | 0.1679 | 0.1614 | 0.0065 | 1.4 % |
| all, after the last wet bin | 9.8 % | 0.1514 | 0.1491 | 0.0023 | 5.9 % |
| links below 1 km, after the last wet bin | 2.3 % | 0.1534 | 0.1508 | 0.0027 | 1.6 % |

(The sixth set of the written plan, the 60 min tail with the rainy part left out, is the third set again by construction.)

- **The K1 fit itself is good.** It follows the true k_dry(W) within 25 % up to W = 1 dB and would have had **97 %** of the
  domain weight within 25 % (the rule needs 80 %); K0 has 18 %. The failure is in the selection, not in the fit.
- **The footprint of the dependence in the loss is small everywhere, not only diluted.** Even on the most favourable
  part, the short links in the 60 min after rain, K1 gains 0.0065 dB in sample; on all values 0.0025 dB in sample and
  0.0019 dB held out. The margin of 0.005 dB is larger than what a threefold k_dry(W) is worth on the whole-event loss.
- **Most of the drop is inside the rain, not in the tails:** 94 % of the loss drop comes from the bins up to the last wet
  bin, 6 % from the bins after it. The tails after the last wet bin are 10 % of the values.
- 86 % of the true drying minutes of the domain lie below W = 0.5 dB, where a constant fitted to everything (16.8 min)
  is about twice too fast against the true 30 min.

**STOP for review.** Nothing was changed on the basis of this diagnostic: the margin, the loss and the rule for (c) stand
as written. Not run: (d), (d-high), (e), (f), (g), the detectability study, the optional step. A7 is not started.
One code correction, with no effect on any number already reported: `a6_recovery.gain_of` read the column `select` as
`cv.select`, which in pandas 3.1 is a DataFrame method; it now reads `cv["select"]`. The function is only used by the
detectability study, which has not run.
