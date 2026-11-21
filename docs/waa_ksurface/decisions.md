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

---

# S1: signal against spurious gain under targeted losses (2026-11-14, unattended)

The review read the stop at truth (c) and diagnostic 3. **A7 (the real fit) is not run in this block.**

## Review note on the failure of truth (c) (2026-11-14, logged before any run)

The failure of (c) is accepted as diagnosed. The whole-series loss cannot separate a k_dry(W) that triples (held-out gain
0.0019 dB) from the gains that noise alone makes on a constant truth (0.0026–0.0028 dB in the first fake set, 0.0013 dB
in the fresh one). Nothing is changed retroactively: the verdicts of A6 final stand as written.

## Deviation 8: R3 revised, k_wet is not identified (2026-11-14, logged before any run)

**What changes.** In the fresh draws of A6 final both estimators of 1 / k_wet gave 3.1–3.7 min for a true 5: the interval
between the gauge-input and the radar-input estimate did not contain the truth. So that interval is not a bound.
**R3 (revised):** k_wet is NOT IDENTIFIED. k_dry is reported as a function of an ASSUMED 1 / k_wet over 2.5–10 min (the
range over which the gate of Deviation 6 showed the CV error to be insensitive), plus the bootstrap. No claim is made on
k_wet. **This is a choice made after seeing a result**, hence a deviation.

## S1: definitions and the adoption rule (2026-11-14, written before any code is run)

**Question served.** Q2 (what is k?): is there a held-out loss on which a real dependence of k_dry on W gains more than
noise does?

**Fake world.** That of A6 (fake6.py), fresh seeds (base 3200; every truth and replicate its own). Truth (a): k_wet = 1/5,
k_dry = 1/15. Truth (c): k_wet = 1/5, k_dry = (1 + W) / 30. Five replicates each. Extension free (Deviation 7).

**The fits (the same for every loss; only the held-out SCORE differs).** Per replicate:
- k_wet: one fit on the gauge-input set of the learning set with the gauge input, wetting K0, drying K1, λ = 1e-4,
  started as in A6 final (the fit the wetting ladder returned in every A6 final truth). 1 / k_wet is held at this estimate;
  the extension of this fit is the start of the next ones. The wetting ladder is not walked (it is not the question).
- Drying K0, and drying K1 at λ = 1e-4, 1e-3, 1e-2 (K1 started from the K0 fit): all links, radar input, trained on the
  loss of A5 (Huber 0.3 dB on the counted 5 min values), 5 folds grouped by event (events sorted by start, dealt
  round-robin), fold refits started from the fit on all events. Exactly the steps of stage 2 of `ladder.py`.
- Error on a set of held-out 5 min values = √(2 · mean Huber loss) in dB. **Gain = error(K0) − error(K1)**, K1 at the λ
  with the lowest held-out error ON THAT SAME SET (as the ladder chooses λ).

**The three losses (sets of held-out 5 min values).**
- **L-all:** every counted held-out 5 min value (the current rule).
- **L-dry:** the counted held-out 5 min values of a link-event from the bin after its last wet bin (`last_bin` of the
  data layer) to 120 min after it (24 bins), in which the path is dry by the rain reference: the radar input of the
  link-event, as the fit sees it, is below 0.1 mm/h in that bin (the wet threshold of the data layer, A4 rule 5). All
  links (with a dry path the length does not enter).
- **L-dryhigh:** L-dry restricted to the segments (one per link-event) whose level at the end of the rain is at least
  1 dB. W is not observed, so the level is defined model-free: the measured A − c of the first 5 min bin after the last
  wet bin (the path is dry there, so A − c = W + noise); that bin must be counted and radar-dry, else the segment is left
  out. Stated limit: the segments are chosen on a noisy reading that is also one of the scored values; both rungs are
  scored on the same values, and what this does to a constant truth is exactly what the (a) replicates measure.

**The two variants.**
- **current:** the learning set (the 14 TRAIN events with a gauge ratio).
- **extended:** the learning set plus the dry-path windows of the other TRAIN events (g_e = 1; the branch dataset holds 51
  of the 53 TRAIN events: 14 + 37). Of an added link-event only the L-dry values count in the training loss and in every
  score; the model is still run through its rain with the radar input, to arrive at a level at the end of the rain.
  Folds: all events, dealt the same way. Stage 1 (k_wet) does not use the added events.
  **In the fake world the radar of an added event reads 2.5 times the true rain** (E1; on top of the T1 error): the true
  rain of an added link-event is its radar value / 2.5, and the fit sees the radar value.
- **The fake check of the extension (its own rule):** on the 5 replicates of truth (a), extended variant, the drying
  constant of the K0 fit must be within 25 % of the truth (1 / k_dry in 12.0–20.0 min) in 5 of 5. If not, the extended
  variant is NOT ELIGIBLE for adoption, whatever its ratio. Reported besides: 1 / k_dry of the current variant on the same
  data, and the coverage gain = the number of L-dryhigh segments and values, extended against current (fake world, mean of
  the (c) replicates; real attenuation is not opened for this).

**Separability and adoption (the rule, written now).** For each of the 6 combinations (loss, variant):
min_c = the smallest of the 5 gains on (c); max_a = the largest of the 5 gains on (a). **Separable if min_c > max_a**
(so all 5 of 5 (c) gains exceed every (a) gain). Ratio = min_c / max_a (if max_a <= 0 and min_c > 0: ratio = infinite,
ties broken by the larger min_c; if min_c <= 0: not separable). **Adopted: the eligible combination with the largest
ratio, if that ratio is above 1.** Nothing else decides.

**Deviation 9 (comes into force ONLY if a combination is adopted).** The drying rung selection uses the adopted loss (and
the adopted variant as the learning set of stage 2): a richer drying rung replaces the current one only if its gain on
the adopted loss is above **margin = max_a on that loss** (not below 0; a measured noise floor, not a chosen number) and
the paired 90 % event-bootstrap interval of the gain excludes zero. The same margin serves every step of the drying
ladder (measured for K1 over K0 only: a stated limit). Wetting keeps the rule of Deviation 5.

**If nothing is separable.** Truth (c) is run at larger sizes in the family of the detectability study,
k_dry = (1/15)·(1 + a·W) / (1 + a): ratio k(W = 2 dB) / k(W = 0) = 5 (a = 2) and 10 (a = 4.5), 5 fresh replicates each,
against the same (a) gains. Reported: the smallest ratio that separates on the best loss (the combination with the
largest ratio at that size, eligible ones only). Then STOP: that is the detection limit and the answer to Q2.

**If adopted (step 3, fresh seeds, base 3300).** (c) again under Deviation 9 (PASS as written in A6 final: drying K1,
wetting K0, k_dry within 25 % on at least 80 % of the domain weight; with Deviation 8 the drying fit is judged with k_wet
at its gauge-input estimate and at its radar-input estimate as before, both must hold); then (d), (d-high), (e), (f), (g)
and the detectability study as written in A6 final, with the drying selection of Deviation 9. For (f) the spurious gain
on the adopted loss is reported against the margin. FAIL on (c), (d) or (f) (spurious k_dry) at the written sizes is a
STOP. Then STOP for review before A7.

**S1, interruption and resume (2026-11-14, written before the resume).** The first S1 run was stopped by the machine
(low memory) after truth (a), replicates 0–2 were cached; replicate 3 was cut off and is made again from scratch. The
resume changes only how the run is executed, no seed, setting, loss or rule: every replicate runs in its own process, one
at a time, so its arrays are freed when its numbers are saved; inside a replicate the arrays of stage 1 and of each
variant are released when no longer needed; the peak memory of every replicate is logged; before a replicate starts the
free memory is checked and, if it is below 2 GB or below the largest peak seen, the run stops cleanly with what is cached.
Check that no number moved: the "current" gains of (a) replicate 3 printed by the interrupted run (L-all +0.00004, L-dry
+0.00017, L-dryhigh +0.00291 dB; 1 / k_dry 13.9 min) must come back the same.

**Result of S1 (2026-11-14).** Script `s1_separability.py`, report `outputs/reports/s1_separability.md`. Fake data only.
The run was finished in a terminal by the students after two stops by the machine; checked afterwards: the code and the
config are identical to commit `8a6fa71`, 60 rows (2 truths × 5 replicates × 2 variants × 3 losses), seeds 3200–3204 for
(a) and 3210–3214 for (c). Peak memory per replicate 2.20 GB (logged from (a) replicate 3 on).

Held-out gain of drying K1 over K0 (dB), smallest to largest of 5 replicates:

| loss | variant | (c) gains | (a) gains | ratio | separable |
|---|---|---|---|---|---|
| L-all | current | −0.00015 to 0.00144 | −0.00011 to 0.00014 | below 0 | no (4 of 5) |
| L-dry | current | 0.00076 to 0.00311 | −0.00027 to 0.00017 | 4.56 | yes |
| L-dryhigh | current | 0.01327 to 0.01904 | 0.00050 to 0.00291 | 4.55 | yes |
| L-all | extended | 0.00056 to 0.00175 | −0.00012 to 0.00012 | 4.65 | yes |
| L-dry | extended | −0.00063 to 0.00079 | −0.00061 to 0.00000 | below 0 | no (3 of 5) |
| **L-dryhigh** | **extended** | **0.01337 to 0.01785** | **0.00046 to 0.00195** | **6.86** | **yes** |

- **Fake check of the extension (radar of the added events 2.5 times too high): PASS.** 1 / k_dry of the K0 fit on (a):
  15.4, 16.9, 14.5, 13.7, 15.0 min (truth 15), 5 of 5 within 25 %; 0.2–0.3 min below the current variant on the same data.
- **Coverage gain at a level of at least 1 dB is small:** L-dryhigh has 101 link-events of 10 events (2 207 values) in the
  current variant and 106 of 13 events (2 322 values) in the extended one (fake world, mean of the (c) replicates). L-dry
  itself doubles (88 611 → 189 039 values, 14 → 51 events). The added events hold little water at the end of their rain in
  the fake world, because their true rain is the radar / 2.5.
- On L-dryhigh the gains on (c) are about ten times those on L-all, but so is the noise floor: on a constant truth K1
  always gains there (0.0005–0.0029 dB), the effect of choosing the segments on a noisy reading that is also scored.
- Stated limits: the adopted combination is ahead of three others (ratios 4.55–4.65) by a margin that rests on the largest
  of five (a) gains; and L-dryhigh rests on about 100 link-events of 10–13 events.

**Verdict by the rule written before the run: ADOPTED, L-dryhigh in the extended variant. Deviation 9 is in force:
margin = 0.00195 dB.**

## Step 3 under Deviation 9: what is fixed before the code is written (2026-11-14)

**Question served.** Q2. Everything is as written in "A6 final" and in the S1 entry; the points below settle what those
texts leave open now that the learning set of stage 2 is the extended one.

- **Fake world.** All TRAIN events of the branch (51), as in S1: the radar of the 37 added events reads 2.5 times the
  true rain. Seeds: base 3300 (+ position of the truth; detectability as in A6 final with this base).
- **Stage 1 (wetting):** unchanged, gauge-input set of the 14 learning events, selection of Deviation 5 (0.005 dB on the
  whole-series loss).
- **Stage 2 (drying):** all links, radar input, the learning set plus the L-dry values of the added events (training
  loss as in S1). Every rung and λ gets its held-out error on L-dryhigh (5 folds over all events); λ per rung = the
  lowest of these; the ladder is walked from K0 with margin 0.00195 dB on L-dryhigh and the paired 90 % event-bootstrap
  interval excluding zero. k_wet as an interval as before (the selected rung is refitted with k_wet at its radar-input
  estimate, which is a plain K0 / K0 fit on the same stage-2 data).
- **Judging stays on the learning set:** the domain weight (minutes of the true state, short links), the after-rain
  residual of (e) and the tail flag of (f) use the link-events of the 14 learning events, as written.
- **(f):** the 3 events with unseen rain are drawn among the 14 learning events (20 %, as written). Reported: the gain of
  the best richer drying rung on L-dryhigh against the margin. STOP if the verdict is "spurious k_dry" (drying not K0, or
  1 / k_dry outside 25 %).
- **Execution:** one process per truth or replicate, cached, peak memory logged, clean stop on low memory (as in S1).

**Step 3, first attempt (2026-11-14).** The run was stopped by the machine (low memory) inside truth (c), before its
selection and judging; nothing of step 3 is cached and there is no verdict. Its log shows the drying ladder of (c) on
L-dryhigh up to the second λ of K3: gain over K0 of K1 0.0135 dB, K2 0.0059 dB, K3 0.0136 dB (margin 0.00195 dB); the
wetting ladder returned K0. Also recorded: the interpreter of the `cmlpipe` environment must be started with the
environment activated (its `Library\bin` on PATH), otherwise numpy's BLAS calls crash; this concerns execution only.

## Deviation 10: non-nested rungs (2026-11-15, logged before step 3 is run)

**What changes (text of the review).** "Non-nested rungs. K1 and K2 are each compared with K0 on the adopted loss. If
neither beats K0 by the margin, K0. If one or both do, take the one with the larger gain. If both beat K0 and differ from
each other by less than the margin, the result is 'dependence detected, attribution W vs R not resolved', and K3 is tried
against the better one. K3 is selected only if it beats the selected rung by the margin."

**Reason (given by the review).** The sequential walk tries K1 before K2, so a true R-dependence could be absorbed by K1
(in (c), K2 gained 3 times the margin on a truth without R-dependence).

**Disclosure.** This deviation is made after seeing the partial log of truth (c) of the stopped step 3 run (seed 3300):
K1 gain 0.0135 dB, K2 gain 0.0059 dB, K3 0.0136 dB on L-dryhigh. It does not change the outcome of (c) on those numbers:
K1 − K2 = 0.0076 dB, above the margin, and K3 gains 0.0001 dB over K1. Truth (c) of step 3 keeps its seed: it never
produced a verdict, and the rule is not tuned to make it pass. **A choice made after seeing a result**, hence a deviation.

**How it is computed (fixed now, where the text leaves room; the reading is this log's, flagged for the review).**
- It applies to the DRYING ladder of step 3 and A7 (the ladder selected on the adopted loss). Wetting keeps the walk of
  Deviation 5.
- "Beats by the margin" is the test of Deviation 9: gain on L-dryhigh above 0.00195 dB AND the paired 90 %
  event-bootstrap interval of the gain excludes zero.
- "Differ from each other by less than the margin": |gain(K1) − gain(K2)| < 0.00195 dB (no interval).
- K3 is tried against the selected rung whenever K1 or K2 was selected (resolved or not); it is not tried if K0 stays.
- The four outcomes, reported in a new column `attribution` next to the rung:
  1. neither beats K0 → **K0**, "no dependence";
  2. one is selected and the attribution is resolved → **K1** "W" or **K2** "R";
  3. both beat K0 and are closer than the margin → the better one is returned as the rung, with **"dependence detected,
     attribution W vs R not resolved"**. For the PASS rules of (c) and (d) the returned rung must still be the true one,
     and an unresolved attribution is reported next to the verdict;
  4. K3 beats the selected rung → **K3**, "W and R".
- Unit test, written with the code: one toy case per outcome (`s3_selection_test.py`).

**Result of step 3 so far (run 2026-11-15 from a terminal; logged 2026-11-15).** Script `s3_recovery.py`, report
`outputs/reports/s3_recovery.md`. Fake data only, seeds 3300 (c) and 3301 (d), code as committed in `9e5b274`
(Deviations 9 and 10). One process per truth; peak memory 2.84 and 2.86 GB; 68 and 48 min.

| truth | regime | true rung | rung returned | recovery | verdict |
|---|---|---|---|---|---|
| (c) | wetting | K0 | K0 (K3 +0.0014 dB, under 0.005) | 1 / k_wet 3.5 (gauge) / 4.7 (radar) min, truth 5: governed by the gate | PASS |
| (c) | drying | K1 | **K1**, attribution "W" | 97 % / 97 % of the domain weight within 25 % | PASS |
| (d) | wetting | K2 | **K0** | 12 % of the domain weight within 25 % | **FAIL** |
| (d) | drying | K2 | **K0**, "no dependence" | 93 % / 93 % of the domain weight (whole domain 92 %) | **FAIL** (the rung) |

- **(c) PASSES.** On L-dryhigh (2 409 values of 109 link-events) K1 gains 0.01348 dB over K0 (interval 0.0089–0.0204;
  margin 0.00195), K2 gains 0.00597 dB, K3 adds 0.00011 dB to K1. K1 − K2 = 0.0075 dB, above the margin: resolved, "W".
  The fitted K1 is within 25 % of the true k_dry(W) on 97 % of the domain weight at both ends of the k_wet interval.
- **(d) FAILS on both rungs.** Drying: on L-dryhigh (1 728 values of 82 link-events) K1 and K2 are both WORSE than K0
  (−0.0033 dB each), so K0 is returned; the constant (9.5 / 10.1 min for a true 10 at R = 0) is within 25 % of the true
  k_dry on 93 % of the domain weight. Wetting: the true K2 gains 0.0016 dB on the whole-series loss of the gauge-input
  set (K3 0.0032 dB), under the 0.005 dB margin of Deviation 5, and the constant K0 is within 25 % on 12 % of the weight.

**Verdict.** FAIL on (d): **STOP** (the script stopped by the rule; not run: (d-high), (e), (f), (g), detectability).
Nothing was patched.

## Deviation 11: scope narrowed after truth (d) (2026-11-15, logged before any further run)

**Text of the review.** "Scope narrowed after truth (d). (d) cannot pass by construction under the adopted rules:
(1) the adopted drying loss L-dryhigh scores only after-rain time (R < 0.1 mm/h), so an R-dependence of k_dry cannot
appear on it; (2) k_wet is NOT IDENTIFIED (R3), so a wetting ladder that must detect R-dependence contradicts R3.
Response: the claims are narrowed, no criterion is loosened.
- Wetting: K0 only, reported over 1/k_wet = 2.5–10 min (R3). The wetting ladder is removed.
- Drying: the ladder is K0 vs K1 only (k_dry as a function of W, after rain). K2 and K3 are removed.
- New restriction R4: no dependence of k on R is claimed. k_dry is the drying law after rain (R = 0). A true k(R) of the
  size in (d) goes undetected, and the drying constant is still recovered within 25 % on 93 % of the domain weight.
- R2 is superseded by R4.
- Truth (d) becomes a robustness check: PASS if the reduced procedure returns K0 for drying (no spurious W-dependence)
  and k_dry within 25 % on at least 80 % of the domain weight. Re-judge the cached (d) result under this rule. Do not
  refit it.
- (d-high) is judged the same way: PASS if no spurious K1 appears."

**Disclosure.** This deviation is made after seeing the results of (c) and (d) of step 3 (the entry above): (c) passed
with K1; in (d) both ladders returned K0, the drying constant was within 25 % on 93 % of the domain weight and K1 was
0.0033 dB WORSE than K0 on L-dryhigh. So it is known in advance that the cached (d) passes the new robustness rule.
**A choice made after seeing a result**, hence a deviation. The stop rule's one diagnostic after the failed (d) is not
used: the review located the cause by construction.

**How it is computed (fixed now, where the text leaves room; this log's reading, flagged for the review).**
- **Stage 1** fits wetting K0 only (drying K1 as a nuisance, λ by the CV as before): it gives the gauge-input estimate of
  k_wet and the start of the extension. No wetting rung is compared. In the fake checks the drying fit is still made with
  k_wet at its gauge-input estimate and refitted at its radar-input estimate, and both must hold (as in every truth so
  far); the report over an assumed 1 / k_wet of 2.5–10 min is the form of the claim in A7 (R3).
- **Stage 2** fits drying K0 and K1 only. K1 replaces K0 by the test of Deviation 9 (gain on L-dryhigh above 0.00195 dB
  and the paired 90 % interval excludes zero). Deviation 10 has nothing left to decide (K2 and K3 are gone); its code and
  test stay.
- **The cached fits of (c) and (d) are kept**, not refitted: in both the wetting ladder had returned K0, and the K0 and K1
  drying fits do not depend on the richer rungs that were also tried, so they ARE the fits of the reduced procedure.
  (c) is re-judged: drying K1 and at least 80 % of the domain weight within 25 %. (d) is re-judged by the rule above. Their
  tables keep the K2 and K3 rows of the run, which no longer take part.
- **Wetting is not judged** in any truth (K0 by construction; k_wet not identified, R3).
- **(d), (d-high):** PASS if drying returns K0 and k_dry is recovered. For (d): within 25 % on at least 80 % of the domain
  weight, at both ends of the k_wet interval, with the share as it was computed in the cached run (the cells inside F's
  measured range: 93 %; the whole domain, 92 %, is printed next to it). For (d-high) the truth is the constant 0.1 per min
  inside the measured range, so k_dry is judged as a number (1 / k_dry within 25 % of 10 min at both ends).
- **(e)** report only, as written. **(f)** "spurious" if drying is not K0 or 1 / k_dry is outside 25 %; "restored" if K0
  and within 25 % with the flagged events excluded. **(g)** "corrupted" if drying is not K0 or 1 / k_dry is outside 25 %
  (the wetting half of the old rule falls away with R3).
- **Detectability:** only k_dry(W) (the k(R) form falls away with R4): the smallest ratio k(W = 2 dB) / k(W = 0)
  detected (K1 returned) in at least 4 of 5 replicates, sizes and seeds as written.
- **STOP rules:** FAIL on (c), (d) or (d-high), or a spurious k_dry in (f). (g) and (e) give restrictions or reports.
- Unit test, written with the code: the drying ladder compares K1 with K0 and nothing else.
- **Execution of the rest of step 3** (asked by the review): one detached process on this machine, output in
  `outputs/reports/s3_terminal.log`, one process per truth inside it, from the cache.

**Review of Deviation 11's readings (2026-11-15).** The four readings fixed above where the review's text left room were
read and ACCEPTED by the review before the rest of step 3 was launched: k_wet in the fake checks (the drying fit judged
with k_wet at its gauge-input and its radar-input estimate), (d-high) judged on 1 / k_dry as a number, (g) judged on
drying only, and the STOP rules (FAIL on (c), (d) or (d-high), or a spurious k_dry in (f)).

**Result of step 3 under Deviation 11 (2026-11-15).** Script `s3_recovery.py`, report `outputs/reports/s3_recovery.md`,
run log `outputs/reports/s3_terminal.log`. Fake data only; seeds 3300–3304; code as committed in `0c8cd3f`. (c) and (d)
are the cached fits, judged again and not refitted; (d-high), (e), (f) were run as one detached process (one child process
per truth). Only drying is judged; margin 0.00195 dB on L-dryhigh.

| truth | true drying rung | returned | K1 gain on L-dryhigh (dB) | recovery of k_dry (k_wet at gauge / radar estimate) | verdict |
|---|---|---|---|---|---|
| (c) | K1 | K1 | +0.01348 | 97 % / 97 % of the domain weight within 25 % | PASS |
| (d) | K2, outside the ladder (R4) | K0 | −0.00332 | 93 % / 93 % of the domain weight (whole domain 92 %) | PASS (robustness) |
| (d-high) | K0 inside the measured range | K0 | −0.00003 | 1 / k_dry 8.9 / 9.9 min, truth 10 | PASS |
| (e) two stores | none | **K1** | +0.01554 | 1 / k_dry 39.4 / 43.4 min at W = 0 | report |
| (f) unseen rain in 3 events | K0 | **K1** | **+0.00309** (interval 0.00280–0.00345) | 1 / k_dry 16.6 / 17.3 min at W = 0, truth 15 | **spurious k_dry** |
| (f), flagged events excluded | K0 | K0 | +0.00137 (0.00069–0.00364) | 1 / k_dry 12.4 / 13.5 min, truth 15 | restored by the exclusion |

- **(d-high) PASSES:** no spurious K1 (871 scored values of 43 link-events; K1 gains nothing), the constant is recovered.
- **(e), report:** on a two-store truth the one-state ladder returns K1 with a gain of 0.0155 dB, as large as the gain of
  the true K1 in (c) (0.0135 dB). The residual in the 30 min after rain, highest minus lowest third of preceding rain at
  matched W, is +0.014 dB (rule B3 asks for more than 0.10 dB). **So a returned K1 does not tell a k_dry that depends on W
  from a second, slow store, and the after-rain residual does not tell them apart either at this size.** Reported, as the
  rule says; proposed to the review as a restriction on what a K1 result may be called.
- **(f): SPURIOUS.** With unseen rain of 1 mm/h for 30 min after the last wet bin in 3 of the 14 learning events, K1 gains
  0.00309 dB on L-dryhigh, above the margin of 0.00195 dB, with an interval that excludes zero: K1 is returned on a
  constant truth. The scored set swells to 6 733 values of 297 link-events (the unseen rain lifts the first dry reading
  above 1 dB) and its error doubles (0.43 dB against about 0.2 dB). The fitted K1 is still within 25 % of the true constant
  on 90 % / 95 % of the domain weight; the failure is the rung.
- **The flag works and the exclusion restores K0:** the tail flag (defined before use, model-free) flags exactly the 3
  changed events and no other (3 of 3, 0 of 10 false). Without them K1 gains 0.00137 dB, under the margin, K0 is
  returned, and 1 / k_dry is 12.4 / 13.5 min for a true 15 (inside 25 %, the lower one narrowly: the limit is 12.0).
- As the review expected, the targeted loss is where unseen water acts: the margin, measured on clean constant truths,
  does not cover it.

**Verdict.** Spurious k_dry in (f): **STOP** by the rule. Not run: (g), the detectability study. Nothing is patched.
**The one diagnostic allowed after a failed check is not used:** (f) carries its own (the refit without the flagged
events), and it locates the cause: the 3 events with unseen rain, found by the flag. Decisions for the review: whether
the tail flag and the exclusion of flagged events become part of the procedure (then (f) is judged on the refit), and
what a returned K1 may be called given (e).

Execution: the detached process started 19:56:42 and exited 22:31:17 (exit code 1 = STOP by a rule), no restart. (f)
took 147 min of wall time but about 13 min of processor time: the machine was asleep in between. Peak memory per truth
2.93–3.16 GB. Free memory fell to 0.5 GB near the end (Chrome grew from 1.6 to 2.7 GB); the run was not affected.

## Deviation 12: the tail flag becomes part of the procedure (2026-11-16, logged before any run of it)

**Text of the review.** "The tail flag becomes part of the procedure. From now on the drying fit excludes flagged events.
Because the procedure changed, every truth is rerun with the flag, on FRESH seeds: (a), (c), (d), (d-high), (e), (f),
same pass rules as written. For each, report how many events the flag removes and whether any were removed in truths
without unseen water (false removals). (c) must still return K1: the flag must not remove a real drying law."

**Disclosure.** This deviation is made after seeing the refit of (f) without the flagged events (step 3 under Deviation
11, entry above): the flag found exactly the 3 changed events (0 of 10 false), and without them K1 gained 0.00137 dB
(under the margin), K0 was returned and 1 / k_dry was 12.4 / 13.5 min for a true 15. So it is known in advance that the
exclusion restored K0 on that draw. **A choice made after seeing a result**, hence a deviation. For that reason no number
of the earlier runs is used for a verdict: every truth is drawn again with a fresh seed (base 3600; seeds used so far:
1200+, 2200+, 3200–3214, 3300–3305).

**The flag, exactly as implemented** (`a6_recovery.tail_flags`, settings `a6.flag` in `config.yaml`; defined in step A6
before its first use, unchanged since):

```python
def tail_flags(S, Dk):
    """Per event: how many short-link link-events qualify, how many of them wait or rise, flagged or not (decisions.md, A6)."""
    D, le, F = S["D"], S["le"], A6["flag"]
    rows = []
    for e in np.nonzero(S["short"])[0]:
        o, n = D["le_boff1"][e], D["le_n"][e]
        end = (D["le_last"][e] + 1) * 5
        if end + 20 > n or end < 5:
            continue
        y = Dk["y1"][o:o + n]
        a0, a1 = y[end - 5:end + 5].mean(), y[end + 10:end + 20].mean()
        rows.append({"event_id": le.event_id.iloc[e], "qualifies": a0 >= F["a0_min_db"], "waits": (a0 >= F["a0_min_db"]) and (a1 >= a0 - F["fall_db"])})
    t = pd.DataFrame(rows).groupby("event_id").agg(link_events=("qualifies", "size"), qualifying=("qualifies", "sum"), waiting_or_rising=("waits", "sum")).reset_index()
    t["flagged"] = (t.qualifying >= F["min_link_events"]) & (t.waiting_or_rising >= F["share"] * t.qualifying)
    return t
```

with `a0_min_db: 0.6, fall_db: 0.15, min_link_events: 5, share: 0.5`. In words: for every short-link (L < 1 km)
link-event of a learning event, `end` = the end of the last wet 5 min bin of the rain reference; a0 = mean attenuation
(1 min values) from 5 min before to 5 min after `end`, a1 = mean from 10 to 20 min after it. A link-event qualifies if
a0 ≥ 0.6 dB; it "waits or rises" if it qualifies and a1 ≥ a0 − 0.15 dB. An EVENT is flagged if at least 5 of its
short-link link-events qualify and at least half of the qualifying ones wait or rise. It is event-level, uses the rain
reference (for `end`), the other links of the same event, and the measured attenuation only (no model). It is computed
on the training data only (in step 3 `S["short"]` = short links of the 14 learning events); the final equation never uses
it (it selects training data; it is not a term of the model).

**The procedure from now on.** The flag is computed on the data before any fit. The link-events of the flagged events get
use = 0 in both stages of `ladder.two_stage`, and their values are dropped from the scored set L-dryhigh. This is exactly
the refit of (f) that was seen (`s3_recovery.run_truth`, branch `"ghost" in spec`); there is no longer a fit without the
exclusion. When nothing is flagged, the fit is identical to the step 3 fit.

**Readings, fixed now where the text leaves room (this log's, flagged for the review).**
1. **Both stages.** The text says "the drying fit"; the exclusion seen in (f) also left the flagged events out of stage 1
   (wetting K0 on the gauge-input set; they are learning events). It is kept as seen, so the disclosure is exact. Stage 1
   gives only the nuisance k_wet and the start of the extension (nothing of it is judged, R3).
2. **Scope of the flag.** As implemented it looks at the 14 learning events only; the 37 added events (which add L-dry
   values to stage 2 and are scored on L-dryhigh) are never flagged. In (f) the unseen water is put in learning events
   only, so the fake check does not test this gap. **Open point for A7:** on real data, unseen water in an added event
   would not be removed. Not changed here (the text asks for the flag as implemented).
3. **False removal** = an event flagged in a truth without unseen water ((a), (c), (d), (d-high), (e), (g) and every
   replicate), or in (f) a flagged event that is not one of its 3 changed events. Reported per run.
4. **Pass rules, as written, judged on the fit of the procedure** (Deviation 11's rules; only drying is judged):
   (a) PASS if drying K0 and 1 / k_dry within 25 % of 15 min at both ends of the k_wet interval (A6 final); (c) PASS if
   drying K1 and k_dry within 25 % on at least 80 % of the domain weight at both ends; (d) PASS if drying K0 and k_dry within
   25 % on at least 80 % of the domain weight (cells inside F's measured range) at both ends; (d-high) PASS if drying K0 and
   1 / k_dry within 25 % of 10 min at both ends; (e) report; (f) "no spurious k_dry" if drying K0 and 1 / k_dry within 25 %
   of 15 min at both ends, else "spurious k_dry".
5. **STOP rules:** FAIL on (a), (c), (d) or (d-high), or a spurious k_dry in (f). (a) was not in the STOP list of step 3
   because it was not rerun there; A6 final makes it one. A STOP ends the run: nothing of Deviation 13, (g) or the
   detectability study is run after it.
6. **Seeds:** the runs [(a), (c), (d), (d-high), (e), (f), (g)] get 3600 + their position (3600–3606).
7. **Execution:** one detached process (`run_s3b_detached.bat` in the cache folder), one child process per run, cached in
   `research/.cache/waa_ksurface/s3b_recovery/`, log `outputs/reports/s3b_terminal.log`. The parent asks Windows to keep the
   machine awake (`SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED)`, cleared at exit); it does not start if
   less than 4 GB is free; before each child it stops cleanly if free memory is below max(2 GB, the largest peak so far);
   every 10 min it logs the free memory, and below 2 GB it logs the programs that grew most since the start. Per run the
   wall time and the processor time are logged (a large gap = the machine slept).

Script: `s3b_recovery.py` (new; reuses the functions of `s3_recovery.py`, `a6_recovery.py`, `s1_separability.py`).

**Review of the readings (2026-11-16, before launch).** The five readings of Deviations 12 and 13 (both stages excluded;
the flag's scope on the learning events only, as an open point; "P5" read as the second-store claim, R5 if not
separable; replicate 0 of B3' = the Deviation 12 run, ratio 3 = truth (c); (a) a STOP gate, STOP rules on the Deviation 12
truths only) were read and ACCEPTED by the review.

**Addition to Deviation 12 (asked by the review, logged before launch): the flag's blind spot is tested.** One more (f)
variant, `f_added`, run after the six truths (before the B3' replicates):
- **Do the flag's inputs exist for the 37 added events? Yes:** the end of rain comes from the rain reference (the radar,
  which every link-event has), and the attenuation and the other short links of the same event exist as for a learning
  event. Measured on the branch data: the added events have 792 checkable short-link link-events (learning: 859), but only
  **21 of the 37 added events have at least 5 of them** (learning: 13 of 14); since the flag needs 5 qualifying link-events,
  the other 16 added events can never be flagged, whatever their tails do.
- So the flag is **extended to all 51 events** (the same code, `S["short"]` = short links of every event) in this variant
  only; the six truths and everything else run exactly as committed. **Correction after the run (2026-11-18):** flag
  extended to all 51 events; 41 checkable; 10 have no short-link link-event with rain (none with the 20 bins after the
  end of its rain that the flag needs), so they cannot be flagged at all.
- Unseen water as in (f) (1 mm/h for 30 min after the last wet bin, unseen by the input) in **7 of the 37 added events**
  (20 %, rounded), drawn among all 37 (not only the 21 checkable ones: that is the gap being tested); seed 3607.
- **PASS** if the extended flag catches every changed event (catch rate 7 of 7) AND drying returns K0 with 1 / k_dry within
  25 % of 15 min at both ends of the k_wet interval. Reported: the catch rate over all changed events and over the changed
  events that were checkable (at least 5 checkable short link-events), the false removals, the drying verdict, and whether
  the contamination is **material** (drying not K0, or 1 / k_dry outside 25 %).
- Not a STOP. **If it FAILS and the contamination is material**, the restriction proposed for A7 is the review's: the
  extended variant may use added events only if the flag covers them; otherwise A7 falls back to the learning events, and
  the S1 margin for that variant (the current variant, L-dryhigh: 0.00291 dB) applies. If it FAILS only on the catch rate
  (not material), it is reported with no restriction proposed, flagged for the review.

## Deviation 13: naming of a K1 result, and the history test B3' (2026-11-16, logged before any run of it)

**Text of the review.** "Naming. A returned K1 is reported as 'non-exponential drying: the drying rate falls as the
antenna dries'. It does not by itself distinguish a W-dependent rate from a second, slow store (shown by (e)). Then
calibrate the history test (B3'): statistic = after-rain residual difference at matched W (±0.15 dB) between the lowest
and highest third of preceding rain amount, from the K1 fit. Null: 5 fresh replicates of (c) (one state). Signal: 5 fresh
replicates of (e) (two stores). Separable if min (e) > max (c) in 5 of 5. If separable: B3' replaces B3 in Block B, with
threshold = max (c). If not: record the finding 'W-dependent rate vs second store not identifiable with this data', and P5
becomes restriction R5. B3's old 0.10 dB threshold is retired either way (it was not measured)."

**Disclosure.** Made after seeing (e) of step 3 (one draw, seed 3304): K1 returned with a gain of 0.0155 dB, and the
after-rain residual at matched W, highest minus lowest third, was +0.014 dB, from the selected fit (K1). The (c) runs of
step 3 stored the same statistic, but it was not read or reported for this purpose. Hence a deviation; fresh seeds.

**Naming (in force from now on).** A returned drying K1 is reported as **"non-exponential drying: the drying rate falls as
the antenna dries"**, never as "k_dry depends on W" alone. By itself it does not tell a W-dependent rate from a second,
slow store.

**B3', how it is computed (fixed now; this log's readings, flagged for the review).**
- **Statistic** = `matched_high_minus_low_db` of `a6_recovery.after_rain`, computed with the **drying K1 fit** (its λ of
  lowest held-out error; k_wet at its gauge-input estimate), whatever rung is selected: short-link link-events of the
  learning events that are NOT flagged (Deviation 12: they are not training data); per link-event the mean residual
  (fake A − c minus the model) over the 30 min after the last wet bin (5 min values that count, at least 3); preceding rain
  = the input rain up to then (mm); thirds of preceding rain over these link-events; matched W = bins of 0.3 dB (±0.15) of
  the fitted W at the end of the rain; statistic = Σ_bins min(n_high, n_low)·(mean_high − mean_low) / Σ_bins min(n_high,
  n_low), high minus low (signed, dB). The same quantity from the selected fit (old B3) is printed next to it.
- **Null:** 5 replicates of (c); **signal:** 5 replicates of (e). Replicate 0 of each is the (c), resp. (e), run of
  Deviation 12 (same truth, same procedure, fresh seed); replicates 1–4 get seeds 3620–3623 ((c)) and 3630–3633 ((e)).
- **Separable** if every one of the 5 (e) values is above the largest of the 5 (c) values (a missing value, no matched
  bin, counts as not above). Then B3' replaces B3 in Block B: "a second state is added only if the B3' statistic of the
  real K1 fit is above max (c)". Else the finding "W-dependent rate vs second store not identifiable with this data" is
  recorded, and P5 becomes restriction R5.
- **"P5"** is not defined in this branch's documents (plan.md, decisions.md). Reading: P5 = the claim that a Block B
  second state would support, that a slow second store can be told apart from a W-dependent rate. **R5 (if adopted):** no
  second store is claimed or excluded; a K1 result is reported only by the name above. Flagged for the review.
- **B3's 0.10 dB threshold is retired** (it was not measured). These replicates are calibration runs: the STOP rules of
  Deviation 12 do not apply to them, but each is judged and reported as its truth.

**Then (g) and the detectability study under the updated procedure** (as written in A6 final and Deviation 11: k_dry =
(1/15)·(1 + a·W) / (1 + a), k_wet = 1/5, ratio k(W = 2 dB) / k(W = 0) = 1 + 2a; detected = drying K1 returned in at least
4 of 5 replicates; the replicates of a ratio stop when decided (4 hits or 2 misses); the sweep ends at the first ratio not
detected). Ratios tried: **3, 2, 1.5** (a = 1, 0.5, 0.25). Reading: ratio 3 (a = 1) IS truth (c) (k_dry0 = 1/30,
w_scale = 1 dB), so its 5 replicates are the 5 (c) runs of B3' (no extra runs). Ratios 2 and 1.5: seeds 3600 + 100 +
10·(position of a in [1, 0.5, 0.25]) + replicate. (g): seed 3606, "k corrupted" if drying is not K0 or 1 / k_dry is
outside 25 % at either end; if corrupted, the per-link question is reopened. Then STOP for review before A7.

## Step 3b: result and disclosures (2026-11-18)

Run detached 2026-11-18 00:26–01:54, exit 0, 24 runs, 1.4 h wall time (`outputs/reports/s3b_recovery.md`). (a), (c),
(d), (d-high) PASS; (e) K1 (report); (f) no spurious k_dry (flag 2 of 13, both with unseen water); `f_added` FAIL,
material (caught 1 of 7, drying K1); B3' not separable ((e) above max (c) in 4 of 5): finding "W-dependent rate vs second
store not identifiable with this data", P5 becomes R5, B3's 0.10 dB threshold retired; (g) at link sd 0.24 dB: k corrupted
(K1); detectability (extended variant): ratios 3, 2, 1.5 detected (5/5, 4/4, 4/4).

**Disclosures.**
- (a) The report was regenerated from the saved results (`s3b_recovery.py --report`, no refit) after the check-report
  line and the 41-of-51 note were added to the script. The regeneration has no machine state, so the line "Keep-awake
  request accepted by Windows: True; gaps of the memory watch above 3 min: 0 (none); lowest free memory seen 0.48 GB;
  low-memory events: 1." was restored by hand, word for word as the run wrote it.
- (b) The terminal log `outputs/reports/s3b_terminal.log` has one extra final line, "=== step 3b ended 2026-11-18
  01:55:23 ===", written by that regeneration, after the real end of the run (01:54:35, exit 0). The log is not edited.
- The first launch (2026-11-16 14:27) did not start: free memory 3.8 GB, below the 4.0 GB rule (logged in the same file).
  During the run free memory fell to 0.48 GB once (other programs; they were closed); no run was affected.

## Deviation 14: A7 uses the learning events only (2026-11-18, logged before any run of it)

**Text of the review.** "From step 3b, the drying fit in A7 uses the 14 learning events only (S1 'current' variant,
margin 0.00291 dB on L-dryhigh); the extended variant is dropped because the flag cannot see the added events (catch
rate 1 of 7). Restrictions: R5 (a returned K1 does not separate a W-dependent rate from a second store; B3 retired), R6
(unseen water in the added events is unchecked, so they are not used)."

**Disclosure.** Made after seeing step 3b (`f_added`: unseen water in 7 added events, caught 1 of 7, drying K1 where the
truth is K0; (g) at 0.24 dB: K1). Hence a deviation. The restriction itself was written before the run (addition to
Deviation 12: "otherwise A7 falls back to the learning events"); this entry makes it the procedure.

**What changes.**
- **Variant:** S1 "current": the 14 TRAIN events with a gauge ratio. The 37 added events take no part in either stage
  or in the scored set (`use` = learning and not flagged; stage-2 values = the counted values of those link-events;
  scored = L-dryhigh of those link-events). Margin = the largest (a) gain on L-dryhigh, current variant, in S1:
  **0.00291 dB** (it replaces 0.00195 dB, which belonged to the extended variant).
- Everything else as in step 3b: tail flag before any fit (learning events), drying K0 vs K1 on L-dryhigh with a paired
  90 % interval that excludes zero, wetting K0 not judged, extension free, k_wet as in step 3b.
- **R5** (adopted, from step 3b B3'): a returned drying K1 does not separate a W-dependent rate from a second, slow
  store; it is named "non-exponential drying: the drying rate falls as the antenna dries"; no second store is claimed or
  excluded. B3 (and its 0.10 dB threshold) is retired.
- **R6** (new): unseen water in the 37 added events cannot be checked by the flag (10 of 51 events are not checkable at
  all, and on the added events it caught 1 of 7), so the added events are not used.
- **The detectability of step 3b does not carry over:** it was measured on the extended variant (51 events in stage 2,
  scored set about twice as large, margin 0.00195 dB); the current variant has fewer scored values and a larger margin
  (0.00291 dB), so its smallest detectable ratio must be measured again (step 3c b).

## Step 3c: fake runs of the A7 procedure (2026-11-18, written before the code)

Script `s3c_a7_variant.py` (new; reuses `s3b_recovery.py` for the flag, the fits and the judging; only the variant
changes). Fresh seeds from **3800**. Fake data only; no real attenuation is read. VAL and TEST are not opened.

- **(a) Link spread, K0 truth** (k_wet = 1/5, k_dry = 1/15, truth (g) with its level shift per link): `g_088` at the
  MEASURED between-link sd **0.088 dB** (E4; it is also the default noise, so this is truth (a) on a fresh seed), seed 3800,
  run first; `g_150` at **0.15 dB** (the E4 threshold), seed 3801. **PASS** if drying returns K0 and 1 / k_dry is within
  25 % of 15 min at both ends of the k_wet interval (gauge and radar estimates), at each spread. **STOP** if `g_088` returns
  a false K1: the per-link question is then reopened by measurement, and nothing after it is run.
- **(b) Detectability for the A7 variant:** k_dry = (1/15)·(1 + a·W) / (1 + a), k_wet = 1/5; ratio k(W = 2 dB) /
  k(W = 0) = 1 + 2a = **1.5, 2, 3** (a = 0.25, 0.5, 1); **five replicates each, no early stop**; seed 3810 + 10·(position
  of a in [0.25, 0.5, 1]) + replicate (3810–3814, 3820–3824, 3830–3834). Detected = drying K1 in at least 4 of 5. Reported:
  the smallest ratio detected. Run only if (a) does not STOP.
- **(c) The flag on the 14 learning events:** checkable = an event with at least 5 short-link link-events that have the
  20 bins after the end of rain the flag needs (from the event timing only; no attenuation). **Removed on real data:**
  the flag decides from the attenuation at the end of the rain (a0, a1), and real attenuation is not read before A7, so the
  removal count on real data is NOT made in step 3c; it is the first number of A7, reported before any fit. Flagged for
  the review.
- **Execution** as step 3b: one detached process, one child per run, cached in `research/.cache/waa_ksurface/s3c_a7_variant/`,
  log `outputs/reports/s3c_terminal.log`, keep-awake, the same memory rules.

**Step 3c, execution setting (2026-11-18, before launch; not a scientific deviation).** By the user's decision the
free-memory rules of step 3c are lowered: the run starts with at least **1.5 GB** free (was 4 GB), and a child run is
started with at least 1.5 GB free (was max(2 GB, the largest peak so far)). Reason: the runs peak at about 3 GB and a run
that dies resumes from the cache, so a low-memory failure costs at most one run (about 4 min). If a child dies, the parent
logs which programs grew during it and stops cleanly; the run is then resumed from the cache. The memory log every 10 min
is kept. Nothing of the procedure, the truths or the seeds changes. Also added to the report before launch: every fake
truth so far already carried the measured link spread (0.088 dB, the default noise), so `g_088` is truth (a) on a fresh
seed; it is still run first, with its STOP rule.

## Step 3c: result (2026-11-18)

Run detached 2026-11-18 11:02–12:26, exit 0, 17 runs, 1.3 h wall time, no child died (lowest free memory 0.76 GB, 8
low-memory log lines; `outputs/reports/s3c_a7_variant.md`). A7 procedure: learning events only, margin 0.00291 dB.
- **(a) Link spread:** `g_088` (0.088 dB, seed 3800) drying K0, 1 / k_dry 16.7 / 16.6 min (truth 15): **PASS**, no STOP;
  `g_150` (0.15 dB, seed 3801) drying K0, 1 / k_dry 13.3 / 13.9 min: **PASS**. The per-link question is not reopened at the
  measured spread (nor at the E4 threshold); the false K1 of step 3b's (g) appeared only at 0.24 dB.
- **(b) Detectability (A7 variant):** ratio 1.5: 2 of 5 (not detected); ratio 2: 5 of 5; ratio 3: 5 of 5. **Smallest ratio
  detected: 2** (step 3b's 1.5 was on the extended variant and does not carry over, as Deviation 14 said). Step 3b's 4 of 4
  at ratios 2 and 1.5: there the replicates stopped once decided (4 hits), so the fifth was never run.
- **(c) Flag coverage:** 13 of 14 learning events checkable; E20150831_0455 is not (its 67 short-link link-events have no
  20 bins after the end of rain). Removal count on real data: not made (first number of A7, before any fit).

## A7: pre-registration of the real fit (2026-11-18, written before any real fit; NOT run)

**Question served.** Q2 (what is k_dry?), on the real TRAIN data of the 14 learning events. Everything below is fixed
before any real attenuation is read for A7.

**Order.**
1. **The flag first.** The tail flag of Deviation 12 (`a6_recovery.tail_flags`, rule a0 ≥ 0.6 dB, fall < 0.15 dB, at
   least 5 qualifying short-link link-events, share 0.5) on the real attenuation of the 14 learning events. Its removal
   count (events checked: 13; events removed, with their ids) is the **first number A7 reports, before any fit**. The
   link-events of removed events take no part in either stage or in the scored set.
2. **The fit** (the procedure of step 3c, Deviation 14):
   - Data: the unflagged learning events (S1 "current" variant); the 37 added events are not used (R6).
   - Loss: whole-event output error, Huber 0.3 dB on the 5 min residuals, weak prior on ln a_link (A5); the measured
     signal is never differentiated.
   - F and P: measured parts fixed (F up to 3.2 mm/h, P up to 1.4 mm/h); free monotone extensions above (Deviation 7).
   - k_wet: constant (wetting K0), not identified (R3): the fit is repeated with 1 / k_wet held at **2.5, 5 and 10 min**.
   - k_dry: ladder **K0 vs K1 only**; K1 selected only if its gain over K0 on L-dryhigh (held-out, 5-fold leave-events-out
     CV) is above **0.00291 dB** with a paired 90 % interval that excludes zero.
   - The rung is selected at each of the three k_wet values. **Reading (written now, flagged for the review):** if the
     three selections disagree, the result is "rung not determined over the k_wet range" and no rung reading below applies.
3. **Report:**
   - the flag's removal count (step 1);
   - the selected rung at each k_wet value, the gain and its interval;
   - **k_dry(W)** (K0: one constant; K1: k_dry as a function of W over the A4 domain) with an **event-bootstrap interval**
     (200 draws over the unflagged learning events, seed 3900) **and** its spread across 1 / k_wet = 2.5–10 min (R3); the
     reported interval is the union of the three bootstrap intervals;
   - the extension knots of F and P above the measured ranges, labeled **effective** (R1).

**Pre-registered readings (the only claims allowed).**
- **K0** → "no dependence of the drying rate on W larger than a ratio k(W = 2 dB) / k(W = 0) of **2** is present"
  (2 = the smallest ratio detected for this variant in step 3c; 1.5 was detected in 2 of 5 and is not excluded).
- **K1** → "non-exponential drying: the drying rate falls as the antenna dries; the mechanism (W-dependent rate or a
  second, slow store) is not claimed (R5); a spread of the baseline between links is excluded as its cause at the measured
  size (0.088 dB, step 3c (a) PASS, also at 0.15 dB)".
- No other claim is allowed: none on k_wet, on k(R), on F or P above the measured ranges, on the added events, or on a
  second store.

**Claims list in force for A7.**
- **R1.** F and P above their measured ranges (3.2 / 1.4 mm/h) are EFFECTIVE quantities for this rain input; no physical
  claim on them.
- **R3 (revised).** k_wet is NOT IDENTIFIED; k_dry is reported as a function of an ASSUMED 1 / k_wet over 2.5–10 min,
  plus the bootstrap. No claim on k_wet.
- **R4.** No dependence of k on R is claimed; k_dry is the drying law after rain (R = 0). (R2 is superseded by R4.)
- **R5.** A returned K1 does not separate a W-dependent rate from a second, slow store (B3' not separable, step 3b); no
  second store is claimed or excluded; B3 retired.
- **R6.** Unseen water in the 37 added events cannot be checked by the flag (caught 1 of 7 in step 3b; 10 of 51 events
  not checkable at all), so they are not used.
- **Flag coverage.** The flag can check 13 of the 14 learning events; **E20150831_0455 is not checkable** (no short-link
  link-event with 20 bins after the end of rain), so unseen water in it, if any, would go undetected. It stays in the fit;
  this is stated with every A7 result.

**Not read in A7:** VAL and TEST. **A7 is not run until the review accepts this pre-registration.**

**Two points for the review.** (i) The fake runs judged k_dry with k_wet at its gauge- and radar-input estimates (step 3b
machinery); A7 holds 1 / k_wet at 2.5, 5, 10 min as R3 (revised) says. The fake runs bracket the truth (5 min) between
3.3 and 6.5 min, inside 2.5–10, but the three fixed values were not run on fake data. (ii) The rule for disagreeing rungs
across k_wet (step 2, last item) is this log's reading.

## Review of the A7 pre-registration (2026-11-18, logged before any run)

- **Accepted:** the rule "rung not determined over the k_wet range" if the selected drying rung differs across 1 / k_wet =
  2.5, 5, 10 min (reviewed and accepted). Whatever the outcome, the rung and 1 / k_dry are reported at all three values.
- **Asked first: step 3d, closing the k_wet gap on fake data** (the fake runs judged k_wet at its estimates, A7 holds it
  fixed; review point (i)).

## Step 3d: the A7 procedure with 1 / k_wet FIXED, on fake data (2026-11-18, written before the code)

Script `s3d_kwet_fixed.py` (new; the procedure of step 3c, `s3c_a7_variant.py`: flag before any fit, learning events
only, drying K0 vs K1 on L-dryhigh, margin 0.00291 dB, extension free). **What changes:** stage 1 (gauge input, wetting K0)
runs as before and gives the start of the extension; then stage 2 is run three times, with 1 / k_wet HELD at **2.5, 5 and
10 min** (the wetting rate parameter set to ln(1 / τ) and held), as A7 will do. Judged as before (k_dry with that fit;
no radar-estimate refit, since k_wet is no longer estimated).
- **Truths:** (a) constant k (k_wet = 1/5, k_dry = 1/15); (c) k_dry = (1 + W) / 30 (ratio 3). **2 replicates each**,
  fresh seeds: (a) 3840, 3841; (c) 3842, 3843. Each replicate is one fake dataset fitted at the three values.
- **PASS** if (a) returns drying K0 with 1 / k_dry within 25 % of 15 min, and (c) returns drying K1, at all three values
  in both replicates. **STOP** if the rung flips at any value (any (a) fit not K0 or any (c) fit not K1): the "rung not
  determined" reading would then be expected on real data, and the review decides before A7. A (a) fit that is K0 with
  1 / k_dry outside 25 % is a FAIL (reported, not a flip).
- Execution: detached, one child per replicate, cache `research/.cache/waa_ksurface/s3d_kwet_fixed/`, log
  `outputs/reports/s3d_terminal.log`, keep-awake, free-memory floor 1.5 GB (the execution setting of step 3c), memory log
  every 10 min. Fake data only; VAL and TEST not opened.

## Step 3d: result — STOP (2026-11-18)

Run detached 12:37–12:51, exit 1 (STOP by the rule after the first replicate; `outputs/reports/s3d_kwet_fixed.md`).
Truth (a) (constant k, 1 / k_dry = 15 min), replicate `a_r0`, seed 3840, flag 0 of 13:

| 1 / k_wet held | drying rung | K1 gain (dB), 90 % interval | 1 / k_dry (min) |
|---|---|---|---|
| 2.5 min | **K1** | +0.00328 (+0.00256 to +0.00509) | 18.2 |
| 5 min | **K1** | +0.00296 (+0.00226 to +0.00478) | 19.3 |
| 10 min | **K1** | +0.00566 (+0.00443 to +0.00870) | 23.4 (outside 25 %) |

**Verdict: STOP.** Truth (a) gives a false K1 at every held value of 1 / k_wet, even at the true 5 min. This is not a flip
ACROSS the k_wet values (the three agree), so the report's line "the 'rung not determined' reading would be expected" is
the wrong wording for this case: the expected failure on real data is a **false K1 at every value**, read as
"non-exponential drying". Replicates a_r1, c_r0, c_r1 were not run. **A7 is not run; the flag was not run on real data.**

**For the review (a reading, not a test):** the margin 0.00291 dB was calibrated in S1 with k_wet ESTIMATED (stage 1),
and step 3c (`g_088`, the same truth (a), k_wet estimated at 3.3 / 6.5 min) gave a gain of +0.00016 dB. With k_wet HELD
the gain on the same kind of truth is 0.0030–0.0057 dB, at or above the margin. A likely reason (not tested): the stage-1
estimate (3.3 min against a true 5) absorbs the timing error of the rain input; a held k_wet leaves that misfit, and the
W-dependent drying rung takes part of it. So the margin is not calibrated for the fixed-k_wet procedure. Options for the
review: (1) A7 with k_wet at its two estimates, the procedure that was calibrated (steps 3b, 3c); (2) recalibrate the
margin with k_wet held (S1 on truths (a) and (c) at 2.5, 5, 10 min); (3) both.

`a7_real_fit.py` was written before this result (flag first, then the fit at the three held values, then the bootstrap);
it is committed unrun.

## Deviation 15: A7 returns to the validated procedure, k_wet estimated (2026-11-18, logged before any run of it)

**Disclosure:** this deviation is written AFTER the step 3d result was seen (truth (a) gave a false K1 at every held value).

"The fixed-k_wet grid (2.5/5/10 min) is withdrawn as a procedure: on truth (a) it returns a false K1 at every value,
including the true one, because a held k_wet leaves the input timing error to the K1 rung (step 3d). A7 uses the
validated procedure of steps 3b/3c: k_wet estimated (gauge-input and radar-input), k_dry selected at each. R3 revised:
k_wet is a nuisance estimate, not a physical wetting time; k_dry is reported at both estimates, plus the event bootstrap.
If the rung differs between the two estimates, the result is 'rung not determined'. Finding kept for the write-up: fixing
k_wet wrongly creates a spurious W-dependent drying law."

**How the validated code path does it (stated now, before the run, for the review).** A7 calls the step 3c path itself
(`s3c_a7_variant.procedure`, the code of `run_one` moved into a function without change; `ladder.two_stage` with the same
arguments: stage 2 sets and margin 0.00291 dB, selector K0 vs K1, wetting K0, extension starts of A6, extension free).
In that path the drying rung is **selected once**, with k_wet at its gauge-input estimate (stage 1 carried into stage 2);
at the radar-input estimate the selected rung is **refitted** (`x_alt`), not selected again. A second selection at the
radar estimate would be the held-k_wet procedure of step 3d at one value, which is not validated, so A7 does not run it.
Consequences: (i) one drying rung is reported, valid at both estimates; the "rung not determined" rule cannot trigger in
this path and is kept only as written; (ii) k_dry (1 / k_dry at W = 0, 0.5, 1, 2 dB) is reported at the gauge-input and
at the radar-input estimate, each with a 200-draw event bootstrap (seed 3900, the same draws, the selected rung refitted
with that k_wet held, extension free), and the reported interval is the union of the two. The k_wet estimates are
reported labeled **nuisance**, no claim on them. Everything else of the pre-registration (bbb5d2f) stands: the flag
count first; learning events only; the readings K0 / K1; claims R1, R3 (revised as above), R4, R5, R6, E20150831_0455.

**The step 3d wording.** The wrong STOP line of `s3d_kwet_fixed.md` ("rung not determined would be expected") is
corrected in this log (step 3d result, above); the report file itself is left as written by the run.

**Addendum to Deviation 15: the flag on real data (logged before the flag is run).** The flag rule (Deviation 12: a0 ≥
0.6 dB, fall < 0.15 dB) is defined on the model's quantity A − c (the target of the fit, `ksurf.prepare`: y − c). It was
only ever run on fake data, where c = 0 and every minute is valid. On real data c is the measured offset of each
link-event (`offset_c_db`) and invalid minutes are NaN. So the flag code reads **y1 − c** and averages the **valid
minutes** of each window (a0: the 10 min around the end of rain; a1: minutes 10–20 after it). A link-event with no valid
minute in either window is not counted (like a link-event without the 20 bins). On fake data (c = 0, no NaN) the result
is unchanged (unit test). Seen before writing this: nothing of the real attenuation; only the code of the data layer.

## A7: result — STOP for review (2026-11-18)

Run detached 19:52–21:52, exit 0 (`outputs/reports/a7_real_fit.md`, log `outputs/reports/a7_terminal.log`), procedure of
Deviation 15. VAL and TEST not read.
1. **Flag (before any fit):** 13 learning events checked, **0 removed**; E20150831_0455 not checkable.
2. **k_wet (nuisance, no claim):** 1 / k_wet = 1.74 min (gauge input), 10.92 min (radar input).
3. **Drying rung: K1**, gain over K0 +0.02013 dB (paired 90 % +0.01525 to +0.02419) against the margin 0.00291 dB; one
   selection, valid at both estimates. 1 / k_dry (min), gauge / radar estimate, 5–95 % event bootstrap (200 of 200
   converged at each), union in brackets: W = 0.5 dB 12.2 / 17.1 (8.7–25.0); W = 1 dB 5.4 / 8.4 (4.2–10.9); W = 2 dB
   6.0 / 8.0 (4.4–11.8); W = 0 81.5 / 95.1 (37.1–278.9).
4. **Reading (pre-registered, K1):** "non-exponential drying: the drying rate falls as the antenna dries; the mechanism
   (W-dependent rate or a second, slow store) is not claimed (R5); a spread of the baseline between links is excluded as
   its cause at the measured size (0.088 dB, step 3c (a) PASS, also at 0.15 dB)." With it: E20150831_0455 not checkable.

**Points for the review (observations, not claims):** (i) both k_wet estimates (1.74, 10.92 min) lie outside the range the
fake runs gave (3.3–6.5 min), so the real data sit outside the k_wet range where the procedure was validated; (ii) the
fitted K1 curve is not monotone between W = 1 and 2 dB (5.4 → 6.0 and 8.4 → 8.0 min), within the bootstrap intervals;
(iii) 1 / k_dry at W = 0 is poorly constrained (interval 37–279 min). Memory fell to 0.72 GB during the bootstraps
(Discord and Chrome grew); no part died.

## Review of A7: the five checks (2026-11-21, logged before any run of them)

**Text of the review.** "Do not change the result. (1) Confirm the order of the Deviation 15 addendum against the first
real read. (2) The selected rung at BOTH k_wet estimates separately, with the K1 gain vs margin at each; if they differ,
the result is 'rung not determined' (pre-registered). (3) From the event bootstrap: share of draws that select K1; 90 %
intervals of 1 / k_dry at W = 0.5, 1, 2 dB at each estimate. (4) Sensitivity to F: refit at F's lower and upper interval
ends (rule A3f); the rung and 1 / k_dry at each. (5) Independent check, no rain input: on the dry-path windows of the
learning events, the median measured A − c at 5, 10, 20, 30 min after the path went dry, against the fitted K1 and a K0
refit; both errors. A check, not a selection." Script `a7_review.py` (new). The A7 result (`a7_real_fit.md`, the cache
`research/.cache/waa_ksurface/a7_real_fit/`) is not changed or rerun; every number below is an addition. Real TRAIN data
of the 14 learning events only (the A7 setting, the A7 flag: 0 removed). VAL and TEST are not opened.

**(1) Order (no run; from the commit times and the session record).** f58a0f5 (addendum, 19:51:30 local) and a6e868a
(flag code, 19:52:02) were committed before the flag launcher started (19:52:08, `a7_terminal.log`; `flag.pkl` written
19:52:20). The session record has no command between the A7 pre-registration and 19:52:06 that loads the real data (the
unit tests use toy data). The real TRAIN attenuation itself was read by earlier steps of this branch (A1–E4, from
2026-11-11, by design: F and P are measured from it); that is not an A7 read. **No deviation.**

**(2) The rung at each k_wet estimate.** Stage 1 exactly as in A7 (gauge input, wetting K0, k_wet estimated; it must give
1 / k_wet = 1.74 min again). Stage 2, drying K0 against K1 (λ = 1e-4, 1e-3, 1e-2 for K1), L-dryhigh, margin 0.00291 dB,
paired 90 % interval that excludes zero, is run as a SELECTION twice: with k_wet held at its gauge-input estimate (A7's own
selection) and with k_wet held at its radar-input estimate (1 / k_wet from the radar-input K0 / K0 fit, as A7; must give
10.92 min again). Both start from the stage-1 extension. **Code check first:** the stage-2 function of the review
(`stage_w`, event weights, see (3)) run at the gauge estimate with every weight 1 must reproduce the committed A7 stage
table (`a7_real_fit_cv.csv`, drying rows: train and held-out loss to 1e-9 relative); if not, **STOP**, nothing after it
is run. **Reading (pre-registered in Deviation 15):** the same rung at both → that rung; different → "rung not
determined", and no rung reading applies. Disclosed now: Deviation 15 called a selection at the radar estimate "the
held-k_wet procedure of step 3d at one value, which is not validated" (step 3d: a false K1 at every held value on truth
(a)). The review asks for it; it is run and reported, and that caveat goes with its result.

**(3) Event bootstrap of the selection.** By the user's decision (2026-11-21; 200 draws with the full selection would take
about 40 h): **50 draws per estimate**, the full stage-2 selection of (2) in every draw (K0, K1 at the three λ, 5 folds,
margin 0.00291 dB, paired interval), at both k_wet estimates. The draws are the FIRST 50 of the A7 bootstrap (seed 3900,
the same sequence); k_wet held at the full-data estimate, extension started from the full-data stage 1 (as A7's
bootstrap). A drawn event counts as often as it is drawn: in the fits (weights, as A7), in the held-out score, and as
separate events in the paired interval (each copy one event). Folds are dealt over the distinct drawn events as in
`ladder.folds`. Reported: the share of the 50 draws that select K1 at each estimate. **The 90 % intervals of 1 / k_dry at
W = 0.5, 1, 2 dB are those of the 200-draw A7 bootstrap** (5–95 %, already in `a7_real_fit.md`; restated, not recomputed);
the 1 / k_dry of the selected rung in the 50 draws is reported besides.

**(4) Sensitivity to F (rule A3f, plan.md: "Sensitivity refits use F at the lower and upper interval ends").** The
interval = the 5–95 % event-bootstrap interval of each F-valid E4 cell (`I_lo`, `I_hi` of `e4_remeasure_cells.csv`); the
curve F_low (F_high) is built exactly as `fp.measured` but through I_lo (I_hi), with the same monotone step and weights; P
unchanged. At each end the whole of (2) is run (stage 1, k_wet re-estimated at both inputs, the selection at both
estimates). Reported: 1 / k_wet, the rung at each estimate (gain vs margin), and 1 / k_dry at W = 0, 0.5, 1, 2 dB of the
selected rung's fit. No bootstrap at the F ends.

**(5) Dry-path check (no rain input in the measured side; a check, not a selection).** Link-events of the unflagged
learning events, all lengths. "The path went dry" = the end of the last wet 5 min bin of the data layer (`le_last`), if
the next bin is radar-dry (the reference rule of S1: the radar input as the fit sees it, 5 min mean, below 0.1 mm/h).
"At least 10 min of rain before" = at least 10 minutes of the link-event before that time with the 1 min radar input at
or above 0.1 mm/h. The value "at ℓ min" (ℓ = 5, 10, 20, 30) = the 5 min bin that ends ℓ min after the path went dry; it is
used if every bin from the first dry one to it is radar-dry and it is counted (`w5`) with a finite A − c. Per ℓ: n, the
median measured A − c, the median model output of the K1 fit (A7's fit at that estimate) and of a K0 refit (the stage-2
drying K0 fit of (2) at the same k_wet estimate, all used data), on the same values. **Error = median(model) − median
(measured)** (dB), for K1 and for K0, at each estimate. Nothing is selected or changed by it.

**Execution.** One detached process, one child per part (`sel_mid`, `dry`, `sel_flo`, `sel_fhi`, then the bootstrap in
blocks of 10 draws per estimate, each draw cached), cache `research/.cache/waa_ksurface/a7_review/`, log
`outputs/reports/a7_review_terminal.log`, keep-awake, free-memory floor 1.5 GB, memory log every 10 min. Order: `sel_mid`
(with its STOP rule), `dry`, `sel_flo`, `sel_fhi`, bootstrap gauge, bootstrap radar. Then STOP for the review.

### Result of the A7 review (run 2026-11-21, exit 0; report `reports/a7_review.md`)

The A7 result is unchanged; nothing below selects or alters it.
- (1) No deviation (addendum committed 19:51:30 / 19:52:02, first real read 19:52:08).
- (2) Code check REPRODUCED (2.2e-16). Rung K1 at the gauge estimate (gain +0.02013, 90 % +0.01525 to +0.02419) and K1 at the radar estimate (+0.01837, +0.01378 to +0.02293), margin 0.00291. Same rung, so "rung not determined" does not apply. The radar selection carries the Deviation 15 caveat.
- (3) 50 draws per estimate (seed 3900): K1 in 50 of 50 at both (share 1.0). 1/k_dry 90 % intervals (A7, 200 draws) at W = 0.5 / 1 / 2 dB: gauge 8.7-17.7 / 4.2-7.3 / 4.4-9.3 min; radar 12.1-25.0 / 6.6-10.9 / 5.2-11.8 min.
- (4) F at I_lo: K1 at both; 1/k_dry 20.7 / 9.0 / 8.6 (gauge), 24.5 / 11.0 / 9.8 (radar). F at I_hi: K1 at both; 7.8 / 3.0 / 3.5 (gauge), 13.8 / 6.6 / 5.6 (radar). The rung does not move; 1/k_dry moves by about a factor 2 across F.
- (5) Dry-path check (not a selection): K1 is above the measured median by +0.06 to +0.12 dB at 5-30 min; the K0 refit is within 0.05 dB (−0.051 to +0.037). The check favours K0 although the CV selection picks K1. This is reported, not acted on.
- Execution: lowest free memory 0.51 GB, 21 low-memory events; all parts completed and cached.

## Diagnostic of the dry-path conflict (2026-11-22, logged before any run of it)

Asked after the A7 review: why the dry-path check favours K0 while the CV selection picks K1. A diagnostic only: no refitting,
no selection, the A7 result and every fit are read from the caches (A7 `fit.pkl` for K1, the K0 refit of `a7_review/sel_mid.pkl`).
Real TRAIN data of the 14 learning events (A7 flag, 0 removed). VAL and TEST are not opened. Script `a7_dryshape.py`.

- **Windows.** The dry-path windows of review check (5), at the 5 min resolution of the model output: the path went dry at
  the end of bin `le_last`, at least 10 wet minutes before, and the bins ending 5, 10, 15, 20, 25 and 30 min after are all
  radar-dry, counted and finite. A window is kept if the measured A − c > 0.15 dB in all six bins. "The first dry minute" is
  the bin ending 5 min after (the model has no 1 min output).
- **(1) Decay rates.** Per window: rate = −(OLS slope of ln(A − c) against time), early half on the bins ending 5, 10, 15 min,
  late half on 15, 20, 25, 30 min; ratio = rate(late) / rate(early). The same on the K1 and K0 predictions in the same windows.
  A window whose early rate is ≤ 0, or whose prediction is ≤ 0 in any bin, has no ratio for that series. Each comparison uses
  the windows where the measured ratio and both models' ratios exist (per k_wet estimate); the number dropped is reported.
  Reported: the median ratio, with a 90 % event-bootstrap interval (events drawn with replacement, all windows of a drawn event
  count once per copy; 1000 draws; seed 3901). Below 1 means drying slows within a window (the K1 shape); about 1 means exponential.
- **(2) Classes.** By the measured A − c of the first dry bin: ≥ 1 dB and < 1 dB; every quantity is given for all windows and
  per class, with n (windows and distinct events).
- **(3) Errors**, per class, model and estimate, each with its bootstrap interval: level error = median(model) − median(measured)
  at the first dry bin; decay-rate error = median ratio(model) − median ratio(measured), the ratio of (1).
- **(4)** 1 / k_dry of the K0 refit (min) at each estimate (K0 has one drying rate), with the A7 K1 values beside it.
- Then STOP.

### Result of the dry-path diagnostic (run 2026-11-22, exit 0; report `reports/a7_dryshape.md`)

No refitting; nothing selected. 1073 windows (13 events) kept; 757 (gauge) / 752 (radar) with all three ratios. 310 dropped
because the measured early rate is ≤ 0 (A − c not falling over 5-15 min).
- (1) Median ratio rate(late)/rate(early): measured 0.227 (90 % 0.108-0.339); K1 0.653 (gauge) / 0.695 (radar); K0 1.000.
  The measured drying slows within a window much more than K1 does; K0 is exponential by construction.
- (2) ≥ 1 dB: 199 / 197 windows, measured ratio 0.347; < 1 dB: 558 / 555 windows, measured ratio 0.156-0.159.
- (3) Level error at the first dry bin: K1 −0.28 / −0.29 dB, K0 −0.34 / −0.35 dB (≥ 1 dB class about −1.0 dB for both).
  Rate error: K1 +0.43 / +0.47, K0 +0.77 (both classes the same order, K1 nearer). On these windows K1 is nearer in shape
  and level; both start too low.
- (4) K0 1/k_dry: 8.88 min (gauge), 13.22 min (radar). K1 (A7 fit.pkl): 81.5 / 12.2 / 5.4 / 6.0 min (gauge) and
  95.1 / 17.1 / 8.4 / 8.0 min (radar) at W = 0 / 0.5 / 1 / 2 dB. The 95.1 vs 94.9 min of review table (2) at the radar
  estimate at W = 0: A7's own radar refit vs the review's refit of the same rung.

## Control for the within-window ratio (2026-11-22, logged before any run of it)

Asked after the dry-path diagnostic: can baseline offsets plus the window selection alone make a plain exponential drying
look like the measured slowing (ratio 0.23)? Fake data only; no fit of any kind (real or fake). Script `a7_dryshape_control.py`.

- **Fake world.** The step 3c / A7 fake world (`fake6.make` as in `s3c_a7_variant.run_one`: the TRAIN links and storms, the
  fake F and P, the radar scaling of the added events, radar and gauge inputs with their errors). Truth K0 (kind `const`),
  two truths, each the A7 K0 pair it comes from: (G) 1/k_dry 8.88 min with 1/k_wet 1.74 min (gauge estimate), (R) 1/k_dry
  13.22 min with 1/k_wet 10.92 min (radar estimate). One dataset per truth, seeds 3920 (G) and 3921 (R).
- **Noise.** As in step 3c: an offset error per link-event, normal, sd 0.136 dB (measured in E3; median 0, 10-90 % −0.17 to
  +0.17 dB), the link shift (sd 0.088 dB), white noise (0.072 dB per minute), the 0.3 dB grid. "Offsets off" = the same
  dataset with the offset sd set to 0 (same random numbers drawn, so the same storms and truth); everything else unchanged.
- **Windows and ratio: exactly the diagnostic's** (`a7_dryshape.windows`, `ratios`, the same config `a7d`): the learning
  events (no flag: the A7 flag removed none), the dry-path rule of review check (5) on the radar input the fake gives,
  A − c > 0.15 dB in all six bins 5-30 min, a ratio only where the early rate > 0. Median ratio with the 90 % event
  bootstrap (1000 draws, seed 3901), all windows and per class (first dry bin ≥ 1 dB / < 1 dB), with n.
- **Reading (pre-registered, from the request).** Per truth, with offsets on: "not explained by offset and selection" if the
  fake K0 ratio's interval (all windows) contains 1, or if the real median lies below the fake interval in both classes
  (real medians of the same class: 0.347 for ≥ 1 dB, 0.156 / 0.159 for < 1 dB, gauge / radar; the overall 0.227 shown beside).
  Otherwise "explained": then we stop to rethink. The real slowing counts as not explained only if both truths say so;
  if they differ, that is reported as it is.
- Then STOP.

### Result of the control (run 2026-11-22, exit 0; report `reports/a7_dryshape_control.md`)

Fake K0 truths, no fit. Median ratio (90 %), offsets on: G (8.88 min) 0.296 (0.229-0.375), 413 windows; ≥ 1 dB 0.615
(0.528-0.685, 42 windows); < 1 dB 0.254 (0.200-0.313, 371). R (13.22 min) 0.436 (0.309-0.578), 435; ≥ 1 dB 0.883
(0.674-0.996, 38); < 1 dB 0.382 (0.284-0.450, 397). Offsets off: G 20 windows kept, 0.098 (−0.134 to 0.407); R 76 kept,
0.721 (0.541-0.838).
- Pre-registered reading: neither fake interval contains 1; the real medians lie below the fake intervals in both classes
  for both truths, so **"not explained by offset and selection"** at both truths.
- Stated with it: exponential drying alone, through this selection with offsets and the 0.3 dB grid, already gives a ratio
  far below 1 (0.30 / 0.44), and the offsets are what puts most windows above 0.15 dB (480 / 562 kept with them, 20 / 76
  without). Overall, the real 0.227 is just below the G interval (lower end 0.229); the clear separation is in the ≥ 1 dB
  class (real 0.347 vs fake 0.53-0.68 and 0.67-1.00), which has far fewer fake windows (42 / 38) than real (199).
