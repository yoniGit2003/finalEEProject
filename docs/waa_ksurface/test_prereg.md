# TEST pre-registration (written 2026-11-23, after Block B on VAL; TEST NOT opened by this branch)

**Question served.** Q3 on held-out storms that no model of this branch or of Part 1 was developed on after its own
freeze: does the equation of this branch beat the existing models on TEST, on the links it was fitted on and on links
it has never seen?

## 1. Disclosure

- **TEST is read a second time.** The first read was Part 1's (S8, October 2026), for a different model (the Part 1
  equation). This branch knew Part 1's TEST verdicts (its reports: P1 failed, P2-P5 passed, the RMSE per model) but never
  read TEST attenuation, rain or link-events.
- **Our equation is frozen as in the VAL run**: `a8_equation_frozen.json` (first Block B record, 37470fa), unchanged after
  VAL; the two recipes as fitted on TRAIN (`b_val_recipes.json`); nothing is refitted on VAL or TEST.
- VAL was read twice (Deviation 16, accepted); the Block B results on VAL are known when these rules are written. They
  changed nothing in any model.

## 2. Models (as Block B)

Ours (the frozen equation through the unchanged `ksurf` kernel), the Part 1 equation (`WAAEquation`, its file), Schleiss
(b) and Pastorek (b) (Part 1's stored TRAIN fits, Part 1's `literature_s8c.Set`), no memory and no WAA (the recipes).
Rain input R = radar path mean / g_e for every model; g_e of a TEST event = `g_path` of Part 1's TEST table
(`test_s8_g_e.csv`, read only; 1 where `path_used` is false). Each model keeps its own lag (Part 1 models 1 min, ours and
the recipes 0). 28 and 38 GHz only (32 GHz excluded for every model, A1e).

## 3. Data: TEST read once

One run builds the TEST link-events with the Part 1 builder (`system_id.build_dataset`, event_split "test", the same
selection rules and offset c, the TEST baselines of Part 1, `test_s8_baselines.csv`, read only), once with the TRAIN
links (link_split "train") and once with the TEST links (link_split "test"). Three sets, each reported on all links and on
links < 1 km:

| set | link-events | per-link factor of the Part 1 models |
|---|---|---|
| **known links** | TEST events, TRAIN links | the stored a_link; a link without one is dropped for every model (as Block B) |
| **new links** | TEST events, TEST links (the 59 links never fitted by any model) | a_link = 1 (none is stored); ours and the recipes need none |
| **gauge-ratio events** | TEST events with `path_used`, TRAIN links | as known links |

Scored values: counted 5 min bins with a finite A − c where every model's prediction is finite (the same values for
every model within a set).

## 4. Metrics and readings (as Block B)

- **M1:** RMSE of the 5 min A − c per model, set and link group; paired event bootstrap, 1000 draws, seed 4100, events
  drawn with replacement from the 18 TEST events, the same draws for every model, set and group; 90 % interval (5-95 %) of
  RMSE(model) − RMSE(ours). **Readings (the only claims):** "X better than ours" if the interval lies entirely below zero;
  "ours better than X" if entirely above; otherwise "tie".
- **M2 (a check, no claim):** after-rain error at the bins ending 5, 10, 20 min after the path went dry (review check (5)
  windows), median(model) − median(measured), interval from the same draws, per set; with the dry-shape caveat.
- **M3 (no claim):** the E4 length test on every model's predictions per set, cells eligible in E4 and on the TEST
  measurements; number of cells where intercept and slope are inside the measured 90 % interval.

## 5. Freeze and execution

- Code: `research/waa_ksurface/t_test.py` (imports the Block B functions of `b_val.py` unchanged). Before TEST: unit tests
  and a dry run of the whole TEST path on TRAIN (the "new links" played by 20 % of the TRAIN links with a_link = 1),
  written to the scratch cache only.
- Freeze record `t_test_freeze.json`: SHA-256 (line endings normalized) of these rules, `t_test.py`, `b_val.py`, the
  frozen equation, the recipe fits, the four Block B records, the Part 1 files used (equation, literature fits,
  TEST baselines, TEST g_e), `ksurf.py`, `a8_forms.py`, `fp.py` and the E4 cells; committed before TEST is read.
- The run checks every hash, then reads TEST ONCE, saves every result before writing the report, and refuses to run a
  second time. **If the run crashes: STOP and report before any rerun.**
