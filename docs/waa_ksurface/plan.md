# WAA k-surface branch: plan and standing rules

Branch `waa-ksurface`, created 2026-11-11 from `main` at `fde431c`. This branch is a side study of the WAA research. Part 1 (the
two-store WAA equation, `research/wet_antenna/FINAL.md`) is frozen, and its TEST set has already been read once. This branch
does not change Part 1.

Branch folders (the repo convention; changed in Step 0b, see Deviation 1 in `decisions.md`):

| folder | holds | in git |
|---|---|---|
| `research/waa_ksurface/` | the code of the branch | yes |
| `research/waa_ksurface/outputs/` | small tables and figures | yes |
| `research/.cache/waa_ksurface/` | large arrays and fits | no (git-ignored) |
| `docs/waa_ksurface/` | this plan and `decisions.md` | yes |

Python environment: `C:\Users\yonic\miniconda3\envs\cmlpipe\python.exe` (the `cmlpipe` conda env). Every script of the
branch is run with it, from the repository root. The base miniconda Python lacks `numba` and `xarray`.

## Goal

One local, interpretable, dynamic WAA equation for a single link, learned from all links pooled, with no assumed power law and
no assumed dynamic form. The only assumptions are:

1. W is one state, driven by the rain.
2. Wetting is separated from drying.

F(R) (steady WAA) and P(R) (path loss per km) come from the model-free length test, not from a fitted model.

## Planned steps

| step | name |
|---|---|
| 1 | theory note |
| 2 | F and P from the length test |
| 3 | coverage of the (R, W) plane |
| 4 | build the model |
| 5 | fake-data recovery |
| 6 | real fit |
| 7 | one-state test |
| 8 | naming the equation |
| 9 | physical consistency |
| 10 | freeze and test |

## Standing rules (they apply to every step)

### Isolation
- Part 1 code, outputs, the frozen equation (`research/wet_antenna/outputs/products/waa_equation_final.json`) and its hashes
  (`fit_sha256` in that file, `research/wet_antenna/outputs/products/test_s8_freeze.json`) are read-only.
- Reuse the existing baseline (v2), offset c, rain input R = R_radar / g_e, link selection and event split (TRAIN / VAL / TEST
  events, TRAIN / TEST links). The list of these files and their paths is in `decisions.md`, Step 0.
- If anything must change, copy it into a branch folder and record in `decisions.md` what was copied and why.

### Justification
- Before every step, write an entry in `decisions.md` with: what the step does, why it is needed (a physical or statistical
  reason, never "to get a better result"), and a dated pass rule.
- After the run, add the result and the verdict to the same entry.
- A choice made after seeing a result is marked **DEVIATION**, with its reason.

### Fake data first
- Every method must recover a known truth on fake data before it is used on real data.

### Stop rules
- If a step fails its pass rule, stop and report. Do not patch it until it is approved.
- Stop after each step and wait for approval before starting the next one.

### Check report
- Every step ends with a CHECK REPORT in this format (the fields are adapted to the step):

```
=== CHECK REPORT: STEP n ===
Branch name and the commit it was created from:
Folders created:
Files created or changed (full paths):
Files outside the branch folders that were changed (must be "none"):
<the step's own checks, each answered yes / no or with numbers>
Anything done beyond the step's tasks (must be "none"):
Open questions for us:
=== END ===
```
