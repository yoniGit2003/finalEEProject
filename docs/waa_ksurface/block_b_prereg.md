# Block B pre-registration (written 2026-11-22 with step 8; NOT run; VAL stays closed)

**Question served.** Q3: does the equation of this branch beat the existing models on held-out storms?

**Status.** These rules are written before any Block B code exists and before VAL is opened. They are hashed together with
the frozen equation and the competitor files (`research/waa_ksurface/outputs/products/a8_block_b_freeze.json`) and
committed. Any later change to a hashed file, or any read of VAL attenuation before the Block B run, is a DEVIATION with a
disclosure. TEST is not opened in Block B.

## 1. The frozen equation (ours)

The file `research/waa_ksurface/outputs/products/a8_equation_frozen.json`, written by `a8_name.py --freeze` from the step 8
cache. Nothing in it is refitted on VAL.

- **Drying:** the closed form accepted in step 8 (its plug-back refit at the radar estimate), or, if step 8 accepted no
  form, the free K1 curve of A7 at the radar estimate (`x_alt`). After rain dW/dt = −k_dry(W)·W.
- **Wetting:** K0, k_wet at its **radar-input** estimate (1 / k_wet = 10.92 min; a nuisance value, R3 revised).
- **F and P:** the measured E4 curves (`fp.measured()` on `e4_remeasure_cells.csv`) inside their measured ranges, and the
  extensions of that same fit above them (EFFECTIVE, R1). No per-link factor. Lag d = 0.
- **The step:** one exact exponential step per minute, W = 0 at the start of every link-event window (the `ksurf` kernel).

## 2. Rain input, the same for every model

R = radar path mean / g_e, 1 min, from the data layer (W1). g_e of a VAL event = its gauge ratio (`g_path` of the Part 1
g_e table, read only); events without a gauge ratio have g_e = 1. Our equation and both recipes below use exactly this R;
the Part 1 competitors are run with this R as their rain input.

## 3. Competitors

| model | what | parameters |
|---|---|---|
| Schleiss (b) | Part 1's stored TRAIN fit (W3), forward model = Part 1's code (`research/wet_antenna/literature_s8c.py`, read only) | `literature_s8c_fits.json` → `fits["Schleiss b"]`, its stored a_link per link, unchanged |
| Pastorek (b) | as above | `fits["Pastorek b"]`, unchanged |
| Part 1 equation | Part 1's stored frozen fit (W3), forward model = Part 1's code (read only) | `waa_equation_final.json`, unchanged |
| no WAA | this branch's recipe: A − c = L·P_full(R(t)), W ≡ 0 | the P extension increments; P measured part fixed |
| no memory | this branch's recipe: A − c = F_full(R(t)) + L·P_full(R(t)); W equals its target every minute (no dynamics) | the F and P extension increments; measured parts fixed |

- **The two recipes are fitted on TRAIN at Block B time**, on the same link-events and 5 min values as our stage 2 (the
  unflagged learning events, `stage2_w5` of A7), radar input, the Huber 0.3 dB loss, L-BFGS-B from the A6 extension starts,
  no per-link factor. Their fitted parameters are written to a file that is hashed and committed (a second freeze record)
  **before** VAL is opened.
- The Part 1 competitors keep their own definitions (a_link, path term, closed forms) as stored; this branch does not
  refit or alter them.

## 4. Data

- **VAL is read once:** one load of the VAL link-events by the data layer (W1: the Part 1 builder with event_split = "val",
  the same link-selection rules, the data-layer offset c), **28 and 38 GHz only** (32 GHz excluded for every model, A1e).
  All models are predicted, the metrics computed and saved in that one run; nothing is refitted or re-read afterwards.
- **The same link-events for every model:** a link-event that any model cannot predict (for example a link without a
  stored a_link) is dropped for all models; the number dropped is reported.
- **Scored values:** the counted 5 min bins (w5 > 0) with a finite A − c; the model value is the 5 min mean of its
  1 min output.
- **Two sets, both reported for every metric:** (S-all) all VAL events; (S-ratio) the VAL events with a gauge ratio. The
  g_e = 1 events carry a radar rain error about 2.5 times larger (E1), so S-ratio is the cleaner comparison.

## 5. Metrics

- **M1, RMSE** of the 5 min A − c (model − measured) over the scored values, on all links and on links shorter than
  1 km. **Paired event bootstrap:** 1000 draws, seed 4100, events drawn with replacement, the **same draws for every
  model**; per draw the RMSE of every model; 90 % interval (5-95 %) of RMSE(model) − RMSE(ours).
- **M2, after-rain error (a check):** on the dry-path windows of review check (5) (the path went dry at the end of the last
  wet 5 min bin; the next bin radar-dry; at least 10 wet minutes before; the bin ending 5, 10 and 20 min after, every bin
  up to it radar-dry, counted and finite): per lag, median(model) − median(measured), with its 90 % interval from the same
  bootstrap draws. **Labeled a check**, with the caveat of the dry-shape control: the data-layer offsets and the selection
  of dry-path windows alone already distort a plain exponential drying (A7 control: within-window ratio 0.30 / 0.44 for a
  K0 truth), so M2 is not a measurement of the drying dynamics.
- **M3, the length test on each model's predictions:** the branch's adopted length test (candidate (ii) of E4: the same
  bins, steady-rain rule, cells and eligibility) run on each model's VAL predictions and on the VAL measurements. Per
  eligible cell: the prediction's intercept and slope against the measured ones with the measured 90 % event-bootstrap
  interval; reported as the number of cells where both the intercept and the slope lie inside the measured interval.

## 6. Pre-registered readings (the only claims allowed)

- From **M1 only**, per set (S-all, S-ratio) and per link group (all, < 1 km): "**X better than ours**" if the paired 90 %
  interval of RMSE(X) − RMSE(ours) lies entirely below zero; "**ours better than X**" if it lies entirely above zero;
  otherwise "**tie**".
- M2 and M3 are reported as numbers; no "better" claim is made from them.
- No claim on k_wet, on F or P above their measured ranges (R1), on a mechanism of the drying (R5), on the 37 added
  events (R6), or from TEST.

## 7. Freeze

SHA-256 of: the frozen equation, these rules, Schleiss (b) / Pastorek (b) (`literature_s8c_fits.json`), the Part 1
equation (`waa_equation_final.json`) and the Part 1 freeze record (`test_s8_freeze.json`), the E4 cells, and the model
code (`ksurf.py`, `a8_forms.py`, `fp.py`), written to `a8_block_b_freeze.json` by `a8_name.py --freeze` and committed.
Block B order: write the Block B code; unit tests on fake data; fit the two recipes on TRAIN; hash and commit their fits;
then read VAL once.
