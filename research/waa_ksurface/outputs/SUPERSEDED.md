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
