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
