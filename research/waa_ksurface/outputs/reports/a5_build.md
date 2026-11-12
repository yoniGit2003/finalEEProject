# A5: the model, built and checked (no real attenuation fitted)

Written by `research/waa_ksurface/a5_build.py`; the model is in `research/waa_ksurface/ksurf.py`. Definitions: `docs/waa_ksurface/decisions.md`, step A5.

## The model

Per link-event, one step per minute, W = 0 at the start of the window:

    v = F_full(R(t));  wetting if v > W, else drying;  k = exp(g(R(t)) + h(W));  W <- v + (W - v) exp(-k);  A - c = W + L P_full(R(t - d))

- F_full, P_full: the measured E4 curve inside its valid range; above it a monotone, continuous extension, linear between the knots **F: [5, 10, 20] mm/h, P: [2, 5, 10, 20] mm/h** (per band; parameters = increments >= 0; above 20 mm/h F is held, P keeps its last slope).
- k: ln k piecewise linear, R knots [0, 0.5, 2, 5, 20] mm/h (in ln(1 + R)), W knots [0, 0.5, 1, 2] dB (wetting) and [0, 0.5, 1, 2, 4, 9] dB (drying); 1 / k between 20 s and 10 h.
- Loss: mean Huber (0.3 dB) of the 5 min residuals over whole events + λ·(second differences)², λ from [0.0001, 0.001, 0.01] by CV (5 folds grouped by event). No per-link factor. Exact gradient, L-BFGS-B.

| wetting | drying | parameters_k | parameters_extension |
|---|---|---|---|
| K0 | K0 | 2 | 14 |
| K1 | K1 | 10 | 14 |
| K2 | K2 | 10 | 14 |
| K3 | K3 | 18 | 14 |

## The rule for d (decisions.md, 2026-11-12)

> d is a whole number of minutes from 0 to 5: the path term reads the rain of d minutes earlier (5 min is one radar bin; a longer lag cannot be placed with this input). d is chosen once per rain input (radar; gauge), with the K0 model (constant k_wet and k_dry, extension free) on the TRAIN learning set: the d with the smallest loss. If the paired 90 % event-bootstrap interval of loss(d = 0) − loss(d_best) contains zero, the lag is not resolved and d = 0. d is then held for every rung, fold and refit. It is never chosen per link or per event.

## Unit tests

| test | number | tolerance | pass |
|---|---|---|---|
| 1 gradient vs finite differences | 4.1938327698685946e-10 | 1e-05 | True |
| 2 constant-k step response (dB) | 4.440892098500626e-16 | 1e-09 | True |
| 3 extension = measured curve at measured knot values | 4.440892098500626e-16 | 1e-09 | True |

1. **Gradient.** 40 real link-events (68185 minutes), a fake target, rung K3 in both regimes, extension free, d = 2: 32 parameters. Largest |exact − finite difference| = 3.54e-12 against a largest gradient entry of 8.45e-03: relative **4.2e-10** (tolerance 1e-05). Every parameter: `research/waa_ksurface/outputs/tables/a5_build_gradient.csv`. Beyond the written test, the same check for the other rungs (wetting / drying): K0 / K0: 6.0e-10; K1 / K2: 1.8e-10; K2 / K1: 2.2e-10.
2. **Step response.** Rain 0 → 4 mm/h for 120 min → 0; 1 / k = 5 min wetting, 20 min drying; F = 1.966 dB. Largest |W − analytic| over 330 minutes = **4.4e-16 dB** (tolerance 1e-09).
3. **Extension.** Knots set to the values of the measured curve: largest |F_full − F| and |P_full − P| on 2 000 values of R up to 30 mm/h = **4.4e-16** (tolerance 1e-09).

| band | quantity | knots_mmh | values_at_knots | largest_difference |
|---|---|---|---|---|
| 28 GHz | F | 5 10 20 | 1.9617 1.9617 1.9617 | 0.0 |
| 28 GHz | P | 2 5 10 20 | 0.1020 0.1217 0.1546 0.2202 | 5.551115123125783e-17 |
| 38 GHz | F | 5 10 20 | 2.0295 2.0295 2.0295 | 0.0 |
| 38 GHz | P | 2 5 10 20 | 0.2246 0.4522 0.8317 1.5905 | 4.440892098500626e-16 |

**Verdict: PASS.** (`research/waa_ksurface/outputs/tables/a5_build_unit_tests.csv`)

