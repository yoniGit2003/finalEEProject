# The wet antenna equation: Part 1, final (2026-10-19)

Part 1 of the project (the WAA equation) is closed. This page points to what it produced.

## The equation

```
A(t) - c = W_f(t) + W_s(t) + L c_band R(t - 1 min)
dW_f/dt = (0.72 u - W_f) / tau_f,   tau_f = 2.6 min while W_f rises (wetting), 10.3 min while it falls (drying)
dW_s/dt = (0.28 u - W_s) / 240 min
u = a_link A_band (R(t - 1 min) / 5)^b
```

- R = the radar path-mean rain divided by the event's radar / gauge ratio g_e (City gauges, >= 1 mm per gauge link, else 1), mm/h.
- c = the median of A over the dry minutes 120 to 30 min before the link's first wet bin (0 if fewer than 30 such minutes).

| number | 28 GHz | 32 GHz | 38 GHz |
|---|---|---|---|
| A_band (dB, level at 5 mm/h) | 2.57 [2.00, 3.21] | 2.36 [1.70, 3.06] | 2.77 [2.25, 3.40] |
| c_band (dB/km per mm/h) | 0.105 [0.058, 0.138] | 0.141 [0.047, 0.190] | 0.194 [0.111, 0.245] |

- b = 0.542 [0.472, 0.626] for all bands.
- a_link: 80 % of links within 0.77-1.24 (1 for a new link).
- Fixed: phi 0.72, tau_f 2.6 / 10.3 min, tau_s 240 min, lag 1 min.
- 90 % intervals from 40 event-bootstrap refits.
- Valid 0-40 mm/h. F above 10 mm/h rests on 0.4 % of the minutes.
- Not identified: P below 1 mm/h.
- TEST RMSE (5 min means) 0.695 dB, against 0.86 / 0.76 / 0.73 dB for no WAA / constant WAA / WAA without memory.

## The files

| what | where |
|---|---|
| the official equation (every number, interval, the fit's sha256, its own sha256) | [`outputs/products/waa_equation_final.json`](outputs/products/waa_equation_final.json) |
| the code: rain series, L, band, a_link in, W(t) and the path term out | [`waa_equation.py`](waa_equation.py) |
| the unit test (reproduces the fit on 3 TEST link-events within 1e-6 dB) | [`tests/test_waa_equation.py`](tests/test_waa_equation.py) |
| the WAA document (15 pages): the problem, data, baseline, model, checks, the test, what is left | [`docs/WAA_research_overview.pdf`](../../docs/WAA_research_overview.pdf) |
| the overview checked statement by statement against the repo | [`outputs/reports/formulate_s9_check.md`](outputs/reports/formulate_s9_check.md) |
| the S9 report | [`outputs/reports/formulate_s9.md`](outputs/reports/formulate_s9.md) |

Use:

```python
from research.wet_antenna.waa_equation import WAAEquation
eq = WAAEquation()                                   # reads waa_equation_final.json
out = eq.run(rain_mm_h, L_km=1.2, band="28 GHz")     # one value per minute: out["W"], out["path"], out["A_minus_c"]
```

## What it does not say

These are listed in the documents (three boxes):

- the WAA level in light rain varies between storms (TEST P1 failed, 4 of 8 cells);
- dew and the end-of-rain misfit (0.2-0.5 dB) are outside the equation;
- whether the correction improves link rain estimates is not tested yet.

The next step is R1, then Part 2 (`docs/ROADMAP.md`).
