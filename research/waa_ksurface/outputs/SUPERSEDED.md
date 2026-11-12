# Superseded outputs (Detachment correction, 2026-11-11)

These outputs of step A3 depended on definitions or numbers inherited from Part 1 (the steady-rain rule, the rain classes,
the cell-eligibility rule, and a fake truth made of the Part 1 closed forms). They are kept for the record and must not be
used. The rescue steps E1-E4 replace them (files starting with `e1_`, `e2_`, `e4_`).

| file | why superseded |
|---|---|
| `tables/a3_length_test_cells.csv` | Part 1 steady-rain rule, rain classes and eligibility rule |
| `tables/a3_length_test_duration.csv` | same bins |
| `tables/a3_length_test_rain_input.csv` | same bins (kept as the observation that motivated E1: the g_e = 1 intercept collapse) |
| `tables/a3_length_test_eiv.csv`, `tables/a3_length_test_eiv_noise.csv` | fake truth = Part 1 closed forms and Part 1 a_link spread; all scatter given to the radar; the check FAILED |
| `tables/a3_eiv_diagnostic.csv`, `reports/a3_eiv_diagnostic.md` | diagnostic of the failed check |
| `tables/a3_length_test_knots.csv`, `tables/a3_length_test_curves.csv` | curves from the superseded cells; never approved (STOP at A3d) |
| `figures/a3_length_test_*.png`, `reports/a3_length_test.md` | figures and report of the superseded step |

Diagnostic only (not superseded, not a foundation): `tables/a1_audit_length_test_vs_part1.csv` (the one-off check that the
rebuilt dataset loads like the Part 1 dataset).

## The first A6 run and its diagnostic (superseded 2026-11-12 by Deviations 5 and 6)

| output | why |
|---|---|
| `reports/a6_recovery.md`, `tables/a6_recovery_*.csv` as committed at `e0e7b4d` | made with the selection rule B1 without a smallest effect; truth (a) FAILED. The files of the same name now hold the rerun under Deviations 5 and 6. |
| `reports/a6_diagnostic.md`, `tables/a6_diagnostic.csv` (commit `28f15d9`) | the one diagnostic of that failed check, made with the old selection rule. Kept as the record that motivated Deviations 5 and 6; rerunning `a6_diagnostic.py` today would apply the new rule. |
