# Roadmap

The one roadmap of the project. **Part A** is the short project checklist (tick an item by changing `- [ ]` to
`- [x]`). **Part B** is the detailed plan of the wet-antenna work; update it in place when a step finishes or the
plan changes. Part B starts with the current plan (the reboot of 2026-10-13, steps S1-S5) and the status table;
the closed drying-tail exploration (steps T0-T18) follows as history. What each finished step decided, with numbers: `research/<module>/METHOD.md`.

# Part A: project checklist

Status: in progress (baseline frozen; Stage 1 WAA: the drying-tail exploration T0-T18 is closed, with no drying law
identified; rebooted on 2026-10-13 as system identification on short links, steps S1-S5; S1, S2, S2b (WAA v2), S3 and S3b done; WAA v3 provisional = M4b-B of S3b; S4.1 / S4.2 / S4.3 and S5a done; final model M4b-B (decision 2026-10-18); S5b not run; S6 done (WAA candidate S6); S7 done 2026-10-19; S6b (the one return to S6) done 2026-10-19; decision 2026-10-19: the S8 candidate is the S6 candidate (M1 to box B); S8 done 2026-10-19 (TEST read once: P2-P5 pass, P1 fails on the antenna level at light rain); S8b done 2026-10-19 (closed form of F adopted); S9 done 2026-10-19: **Part 1 (the WAA equation) closed**, `research/wet_antenna/FINAL.md`; next R1, then the path law).

## 0. Data layer (pipeline/): done, release `openmrg_v1`
Open extensions:
- [ ] Wet/dry classification from the CML signal itself
- [ ] Noise and outlier handling beyond QC flags (level shifts, pair mismatch)
- [ ] Link selection for Stage 1 and Stage 2
- [ ] OpenRainER adapter (sub-links named `channel1` / `channel2`)
- [ ] Check whether OpenMesh has a rain reference before using it

## 1. Baseline and attenuation A(t)
- [ ] Baseline per sub-link from the certified dry periods (constant / rolling) and the dynamic minimum-attenuation method [8]
- [x] Check that the baseline does not follow the drying tail (otherwise WAA ends up in the baseline)
- [x] Baseline residual in dry periods: size and spread per sub-link (the noise floor of δ(t))
- [x] A(t) = TSL(t) - RSL(t) - baseline per sub-link, as the input of the next stages

## 2. Stage 1: WAA
- [x] Find links with a drying tail after rain ends
- [ ] Characterise what a drying tail contains (daily cycle, dew, path rain, start of drying) before fitting a law
- [ ] Fit literature WAA models [1], [9] per link
- [ ] Learn a dynamic WAA model with weak-form SINDy (rain as input)
- [x] System identification of WAA on short links (S1-S5), formulated as the final equation (S6-S9, closed 2026-10-19: `research/wet_antenna/FINAL.md`)
- [ ] Compare the learned models across links (same form, different coefficients?)
- [ ] Drying time from SINDy vs drying time measured directly from tails
- [ ] When this stage is finished: publish the fit caches (`python pipeline/get_caches.py --pack`, upload the zip
      as a GitHub release, put the tag and the sha256 into `pipeline/get_caches.py`, remove the "not published yet"
      notes in `docs/REPRODUCE.md` and `README.md`)

## 3. Stage 2: path law f(L, R̄, Z)
- [ ] Learn f(L, R̄, Z) with PySR, radar along the path as reference (WAA removed first)

## 4. Rain estimation and evaluation
- [ ] Invert the combined model (A_rain + WAA) to get rain
- [ ] Normalized bias and normalized RMSE vs gauges and radar, split by link length
- [ ] Define and compute the false-rain-after-event metric
- [ ] Compare with power law, literature WAA models and the RNN baseline [11]
- [ ] Bootstrap stability of learned terms across links and events
- [ ] Final report

# Part B: the wet-antenna (WAA) plan

## Reboot (2026-10-13): system identification on short links

**Why.** The drying tail was the wrong object.
- A tail starts after the LAST rain of an event that lasts 9 h (median) and has several showers. By then A is already
  at about 19 % of its peak (T18).
- There is one tail per event, and its start is uncertain.
- Wetting and drying were split artificially: most of the fall of A happens while the rain fades.
- There was no success criterion, so the checks did not converge.

**What we hold on to (from T14-T18).**
- On links < 1 km, A is mostly WAA in light / moderate rain: the fast drop does not scale with L from 0.46 to 4.9 km.
- WAA follows the rain within minutes (peak A vs peak radar rain: Spearman 0.85).
- There are two phases after rain.
- The rises are water, not drift.
- There is enough data: 41 links < 1 km, about 1660 TRAIN link-events, 28 and 38 GHz.

**The approach.** Rain R(t) in, attenuation A(t) out, on short links, on whole events (all showers, no tail cutting):

    A_pred(t) = W(t) + L * P(R(t)) + offset per link-event

- W(t) = the wet-antenna attenuation (dB), the state we want the law for.
- P = the path term per km. It is learned freely, with P(0) = 0, and is identified ONLY through the length differences
  between links (L -> 0).
- **RULE: WAA is never estimated by subtracting a power-law rain term.** The ITU power law is a size reference at most.

**The model ladder** (each rung adds one thing to the one before):

| Rung | Model |
|---|---|
| M0 | no WAA |
| M1 | static: W = F(R) |
| M2 | one memory time constant: dW/dt = (F(R) - W) / tau |
| M3 | separate wetting and drying speeds (Schleiss-like) |
| M4 | two stores (fast + slow) |
| M5 | flexible G(R, W), then named with E-SINDy / PySR |

**The judge.** Leave-one-event-out prediction error of A (5 min means), with event-bootstrap intervals.
**Stop rule: the simplest model within the uncertainty of the best is the answer.**

**New steps.**
- **S1** fake check + the ladder M0-M4 on the real short links.
- **S2** the discovery rung M5.
- **S2b** repair of the S2 model (P / F trap, bounded slow store, longer fits); added 2026-10-14.
- **S3** verification: VAL, both directions, long links -> Stage 2.
- **S4** conditions (night, humidity, temperature, wind, frequency) on the residuals.
- **S5** the final equation.

**Discipline.** No new side checks unless they block an S-step. Every S-step ends with its deliverable (status table).

## WAA v1 (frozen 2026-10-14, from S1 / S2)

The fallback answer. Later steps can only replace it with a better-checked v2.

```
A(t) = W(t) + L * P_band(R(t - d)) + c          P free (never falling, P(0) = 0), told apart from W only by L
M3:  dW/dt = (a_link F(R) - W) / tau,  tau = tau_wet while W rises, tau_dry while it falls
M4:  W = W_f + W_s;  W_f -> phi a_link F(R) (tau_f_wet / tau_f_dry),  W_s -> (1 - phi) a_link F(R) (tau_s)
wetting   tau_wet about 3 min (GAUGE input: M4 2.7 [2.4, 3.1], M3 3.0); the radar's 15-50 min is input
          smoothing (fake check Ta: a 5 min noisy input makes the wetting 2.1 x slower)
drying    fast part about 5-13 min (gauge 10.6 [9.5, 12.2]; radar 5.5-8.3)
level     F about 3 dB by 2 mm/h with the radar (a_link 0.8-1.2 between links); with the gauge still rising to
          about 3.8 dB at 5 mm/h
bands     28 / 32 / 38 GHz alike (within the intervals)
leftover  a slow part after the rain: 2-4 % of the peak, time constant about 1-2 h (no rung gets it: M3 too
          fast, M4 / M5 too slow)
limits    the F / P split at light rain (an optimiser trap, fake check Tc failed); peaks too low by 0.7-0.9 dB
          (rain input); M5 not identified (fake check Td failed), no formula named
```

## WAA v2 (S2b, 2026-10-14)

Copied from the S2b report (`research/wet_antenna/outputs/reports/system_id_s2b.md`, section 5). Replaces v1 until S3
says otherwise.

```
WAA v2 (S2b, 2026-10-14): model M3, chosen by the stop rule on S<1, radar (the main set); 90 % event-bootstrap intervals
A(t) = W(t) + L * P_band(R(t - d)) + c
dW/dt = (a_link F_band(R) - W) / tau,  tau = tau_wet while W rises, tau_dry while it falls        (exact steps per minute)

                          radar, S<1 (main)                 GAUGE input (gauge links)
lag d                     +2 min                            -3 min
tau_wet                   43.0 [14.9, 274.0] min            2.9 [2.5, 3.6] min   <- the wetting speed (the radar value is input smoothing)
tau_dry                   24.7 [14.2, 114.7] min            12.8 [11.3, 15.5] min
F[28 GHz](0.5 / 2 / 5)    0.83 / 1.40 / 1.53 dB             0.64 / 1.97 / 2.79 dB
F[38 GHz](0.5 / 2 / 5)    0.88 / 1.81 / 1.81 dB             0.81 / 2.48 / 3.47 dB
P[28 GHz](2 / 10)         1.31 / 3.11 dB/km                 0.39 / 0.54 dB/km
P[38 GHz](2 / 10)         1.65 / 4.22 dB/km                 0.59 / 0.85 dB/km
a_link 10-90 %            0.65 - 1.32                       0.80 - 1.24

Within the interval on the main set: M4b (fast store 2.6 [2.4, 3.1] / 10.3 [9.2, 11.8] min with the gauge,
phi 0.72, tau_s at the bound 240); on S<3 the stop rule takes M4b.
Limits that remain:
- the F / P split depends on the model and the input: P[28 GHz](2) = 1.31 (M3, radar), 0.13 (M4b, radar),
  0.39 (M3, gauge) dB/km; ITU size 0.38. With the radar, M3's slow wetting pushes the fast part of W into P.
- the slow leftover after the rain (2-4 % of the peak, 1-2 h) is not described: M3 ends too early, M4b's store
  (tau_s at the bound) too late;
- peaks too low: -0.83 dB (radar), -0.70 dB (gauge);
- drying: gauge 12.8 min against radar 24.7 min, not one number yet;
- intervals from warm-started bootstrap refits (30 iterations) are too narrow along flat directions (section 4).
```

Correction to v1: with the gauge the level is NOT still rising at 5 mm/h; that was the S2 P/F trap. Level ~2-3 dB by
2 mm/h with both inputs.

**Length-rule warning (S2b):** M3's P halves from S<1 to S<3 (28 GHz, 2 mm/h: 1.31 -> 0.57 dB/km, intervals apart), so
M3 hides antenna wetting in P with the radar input. M4b's P is stable. RMSE M3 vs M4b differs by 0.003 dB (the size of
the fold-cap effect). The choice between them is made by S3.1, not by RMSE.

## Remaining plan

S2b (repair and finish the S2 model: the P / F trap, a bounded slow store, longer fits) -> S3 verification (done) ->
S3b (fit on all link lengths so that P is pinned by L; done) -> S4.1 (why events differ) -> S4.2 (conditions) -> S4.3 (done; S4 done) -> S5a (done) -> S6 (done) -> S7 (validate without TEST, sub-steps 7.1-7.6 below) -> S8 (freeze, then TEST once) -> S9 (formulate) (S5b not run; plan from S6 on below). A step may get sub-steps when needed, but every sub-step is added to this list and to the
status table BEFORE it runs, with its reason. (S2b, added 2026-10-14: S2 found an optimiser trap in P and F and a slow
store that overshoots; both must be fixed before S3 verifies anything.)

- **S2b: done 2026-10-14** (WAA v2, above).
- **S3 sub-steps** (added 2026-10-15, before running; reason: the S2b length-rule warning means the F / P split must be
  decided without a fit, and VAL must be read only with frozen models). TRAIN-fitted M3 and M4b, NO refit except S3.4:
  - **S3.1** model-free length test: steady-rain 5 min bins, per (band, rain class) fit (A - c) = I + s L over ALL links
    of all lengths (TRAIN events, VAL repeated in S3.2); I = antenna part, s = path part per km; compared with a_link F
    and P of M3 / M4b; verdict "consistent" per model.
  - **S3.2** VAL events (first read): frozen M3 / M4b (S<1 and S<3 radar fits) on VAL links <= 3 km, errors as S2b
    against the TRAIN leave-one-event-out numbers; M1 for scale; S3.1 repeated on VAL.
  - **S3.3** both directions: the two sub-links of a link (same antennas): a_link, after-rain integral and peak per
    direction (Spearman >= 0.6 = the same wetting).
  - **S3.4** gauge replication on VAL: the only refit, M3 and M4b with the GAUGE input on the VAL gauge links;
    expected tau_wet 2-5 min, tau_dry 8-18 min.
- **S3: done 2026-10-15** (`research/wet_antenna/outputs/reports/system_id_s3.md`). Results:
  - S3.1: model-free I(R) saturates ~1 -> ~2.5 dB, s(R) grows ~with R, q ~ 0; no constant L0 (7-12 km at light rain, 1-2 km at 10-20 mm/h).
  - M3 and M4b both NOT consistent (M3: P 3-8 x s; M4b: P ~0 at light rain, a F 1.2-1.5 x I); the same with R = 15 min mean, late bins, VAL.
  - S3.2: VAL errors as TRAIN (L < 1 km: M4b 0.360, M3 0.365, M1 0.378 dB; worse during rain); the S<1 M3 fit breaks on 1-3 km links.
  - S3.3: directions: peaks agree (Spearman 0.96), the after-rain integral is too noisy (0.34); a_link is shared by construction.
  - S3.4: gauge on VAL: drying 9-11 min replicated (M3 10.8, M4b 8.8), wetting ~2 min (1.8 / 2.0). WAA v3 candidate: neither.
- **S3b** (added 2026-10-15, before running; reason: S3 showed that P cannot be identified from links < 1 km alone): fit the S2
  frame (M3, M4b; M1 for scale) on the TRAIN events of ALL link lengths so that P is pinned by L, free and with the gauge time
  constants fixed; judged on VAL by the length test through the model and by the errors per length class. Then S4.
- **S3b: done 2026-10-15** (`research/wet_antenna/outputs/reports/system_id_s3b.md`). Results:
  - The all-length fits reproduce the TRAIN length test (M3-B, M4b-B 10 / 10 cells; M4b-A 9 / 10, M3-A 5 / 10).
  - Every fit (incl. M1 and the S2b fits) fails on VAL (0-1 / 8): VAL attenuates less per radar mm/h (observed / model I 0.72, s 0.61).
  - Residual during rain flat with L on VAL, but about -0.1..-0.2 dB at every length.
  - Radar tau_wet stays long with P pinned (M3-A 56 min, M4b-A 13 min): input smoothing, settled.
  - Gauge timing (B) gives better peaks (-0.47 against -0.87 dB) and a worse first hour after the rain.
- **S4 sub-steps** (added 2026-10-15, before running; reason: S3b showed an event-set difference that no event-independent model
  can follow):
  - **S4.1** why do events differ: event facts and the split, a random-split test of the TRAIN / VAL gap, per-event rain-scale k_e
    and antenna factor alpha_e on the frozen M4b-B (fake check first), radar vs gauge per event.
  - **S4.2** conditions behind the factors (descriptive only): the condition data already in the repo, per event, against k_e, alpha_e.
- **S4.1 + S4.2: done 2026-10-15** (`research/wet_antenna/outputs/reports/conditions_s4.md`). Results:
  - The event split is random per intensity tercile (seed 42); VAL (6 light / 6 moderate / 6 heavy) is spread over the summer.
  - The VAL gap is not ordinary event spread: 0-0.5 % of 200 random 18-event TRAIN halves are as low (1-1.5 % half vs the other 33).
  - Per-event factors: **fake check FAILED** (63 %; heavy events 16 / 16, light / moderate not), so k_e, alpha_e are not interpreted.
  - Radar / gauge about 1.7 x on the large VAL events against about 0.9 x on TRAIN, the size of the gap (a lead: 5 VAL events).
  - No condition linked (0 of 28 BH discoveries, uninterpretable); residual equal across bands and day / night, more negative in the warmest events (a lead).
- **S4.3** (added 2026-10-15, before running; reason: S4.1 points at the radar and found the factors identifiable only on heavy events):
  the VAL length test with the radar rescaled by the event's radar / gauge ratio (from the City gauges only), and per-event factors on
  heavy events only after a new heavy-event fake check.
- **S4.3: done 2026-10-17** (`research/wet_antenna/outputs/reports/conditions_s4c.md`). Results:
  - Dividing the radar by the event's radar / gauge ratio g_e closes the VAL gap (M4b-B: I 0.72 -> 0.99, s 0.61 -> 0.95; both inside 5 / 6); TRAIN not broken (1.11, 0.92, 7 / 10).
  - The model-free VAL / TRAIN curves move together (I 0.63 -> 0.84, s 0.67 -> 0.91); 28 GHz 2-5 mm/h stays low.
  - Heavy-event factors identifiable (fake check 98 %); k_e follows g_e (Spearman 0.64, 18 events).
  - No condition linked to alpha_e (22 heavy events, 0 of 6 BH); hints only: peak / mean rain rho -0.49, humidity +0.37.
  - Caveats: only 5 VAL events carry g_e; the radar is gauge-adjusted (circular); the VAL 5-10 mm/h cells drop out after rescaling.
- **S4: done** (S4.1-S4.3).
- **S5a** (added 2026-10-17, before running): refit WAA v3 (M4b-B) on the gauge-corrected rain R / g_e, one last candidate for the slow leftover (M4c), closed forms for F(R) and P(R) and the plug-back check -> WAA v4 or "free curves kept".
- **S5a: done 2026-10-18** (`research/wet_antenna/outputs/reports/final_s5a.md`). Results:
  - Corrected rain (R / g_e) is what passes the VAL length test (5 / 6 on V1, at most 1 / 8 on V0, also for the S3b fit that never saw it); the refit moves WAA v3 little.
  - Most of the slow leftover was radar over-read (M4b-B after-rain integral 0.79 / 1.21 x on V1); the median tail is still high after ~40 min.
  - M4c passed its rule but its slow store is not identified (tau_s_wet at the 600 min bound, phi 0.19-0.84).
  - P ~ c R (28, 32 GHz), c R^0.88 (38 GHz), PySR agrees; F saturating, poorly determined; plug-back +0.002 dB.
  - V1 lowers the shower peaks to -1.05 dB (V0 -0.53 dB).
- **DECISION (2026-10-18, before any TEST read; a deviation from the S5a rule, which kept M4c):** the final model is M4b-B (gauge
  timing fixed); M4c is dropped because its slow-store parameters are not identified on real data and it gains only 0.0005 dB on
  VAL. Said openly: this was decided after seeing VAL; the reason is identifiability, not the VAL score.
- **S5b** (added 2026-10-17): NOT run (2026-10-18). Its parts move into S6 (closed forms, intervals), S8 (freeze, TEST once) and S9 (write-up).
- **Plan from S6 on** (2026-10-18; under the integrity rules below):
  - **S6** tighten the equation: final model M4b-B on V1, data funnel, well-determined parametrisation, honest error bars, identifiability, closed forms (physics candidates and PySR compete).
  - **S6: done 2026-10-18** (`research/wet_antenna/outputs/reports/tighten_s6.md`). Results:
    - Shared shapes F = A_band G(R), P = B_band H(R) kept by the VAL rule (0.3875 against 0.3878 dB, length test 5 / 6); TRAIN AICc prefers band shapes.
    - 60 cold draws (the budget cut to 40 undone, dated): A = 2.74 / 2.56 / 2.96 dB at 5 mm/h (+-0.4-0.8); the profile shows the data fix a_link F, not F (the F level is set by the a_link prior).
    - Fake: 62 / 62 within +-20 %; P (and H) below 1 mm/h NOT identified (interval at 0); phi, tau_s fixed.
    - P = c_band R (0.107 / 0.142 / 0.197 dB/km per mm/h; plug-back holds). F stays free shared knots: the AICc form (R / 5)^0.47 (PySR) failed the length test at 4 / 6 = exactly 2 / 3 against the written 0.6667 (to decide).
  - **S6 results in short:**
    - F level A = 2.7 / 2.6 / 3.0 dB at 5 mm/h (28 / 32 / 38 GHz) with one shared shape G(R).
    - F no longer clearly saturates (hyperbola +11 AICc against a capped power; R^0.5 as good); F above 10 mm/h rests on 0.4 % of the minutes, so no F_max is claimed.
    - P = c R, c = 0.107 / 0.142 / 0.197 dB/km per mm/h; P below 1 mm/h is not identified (8 numbers).
    - The data fix a_link F, not F alone: the F level assumes median a_link = 1.
    - The intervals are dominated by the event-to-event rain reference (bootstrap 1-13 x wider than the profile).
  - **DECISION (2026-10-19, before any TEST read):** the S6 F closed form failed the plug-back length test at 4 of 6 = 0.66667
    against the written 0.6667, which is a rounding of the written number; the S3b wording is "at least 2 / 3". The S6 candidate is
    NOT changed (free shared F, P linear). From now on every length-test threshold is written "at least 4 of 6 cells". The closed
    F is carried in S7 as a second description of F.
  - **S7** validate the frozen S6 candidate without TEST (added 2026-10-19; rules in config `wet_antenna.validate_s7`, written before the run; nothing refitted):
    - 7.1 W independent of L: a_link against L and the link metadata (Spearman, link bootstrap, BH q 0.10); W height a_link F at 2 and 5 mm/h per length class.
    - 7.2 short links on VAL (frozen): RMSE per length class, band and rain class; the VAL length test (at least 4 of 6 cells); shared against band shapes on short links; the residual after the rain.
    - 7.3 long links without a power law: (a) dry-path moments (A - c = W measured against the model W); (b) paired links W = (L2 A1 - L1 A2) / (L2 - L1), counted and noise-checked first; (c) the intercept I(R) from short-only and long-only links; (d) does any of them separate a saturating F from R^0.5; (e) c_band against ITU-R P.838 as a size reference only.
    - 7.4 dew: size and timing of the dew W, dew W against a_link, events with dew before the rain against events without, dew water inside c.
    - 7.5 leftovers: a_link against link metadata and the per-event drying time against the weather (Spearman + BH, then PySR with a leave-one-out and a permutation guard); a hit is a hint only.
    - 7.6 rain reference: radar (V0, V1) against the City gauges by rate class, binned by the gauge and by the radar; decision point for a gauge-based g(R).
    - Return rule: back to S6 (at most once, dated) only if 7.1 finds a dependence on L, 7.3a or 7.3c fails, 7.6 finds a rate dependence, or 7.5 finds a relation that beats the null.
  - **S7: done 2026-10-19** (`research/wet_antenna/outputs/reports/validate_s7.md`). Results:
    - W does not depend on L: a_link against L rho +0.05; W height 2.8 dB at 5 mm/h in all three length classes; dry-path W (A - c with no rain on any path pixel) error -0.18 dB at 10 min in every class; intercepts from short-only and long-only links overlap 13 / 15. Paired links not testable (4 pairs).
    - No rate-dependent radar error (the two binnings move apart: regression to the mean): the parked g(R) is not justified. P = c R is 0.5-0.7 x ITU-R P.838 (size only).
    - Failed / misfit, reported: the RMSE rule fails for L > 3 km (1.18 x; also 1.12 x on TRAIN); after heavy rain the model is about 0.5 dB high on dry paths (both F forms; they cannot be told apart). Dew W 1.7 dB for hours, weakly tied to a_link; it raises c by about 0.1 dB.
    - The per-event drying time falls with the wind (17 events, beats the permutation null): return reason (d).
  - **RETURN TO S6 (2026-10-19, the one allowed return, by the S7 return rule (d)):** test, under the S6 rules, whether a wind-dependent tau_f_dry is identified (fake check), improves the model out of sample and moves the other numbers; look at the heavy-rain dry-path misfit. Not started; the S6 candidate stays frozen (hash in the S7 report) until then.
  - **S7 in short:**
    - W independent of L confirmed three ways: a_link against L rho +0.05; the W height is equal in the length classes; the dry-path error is equal in the length classes; the intercept from short-only and from long-only links agrees in 13 / 15 cells.
    - The model is HIGH after the rain: about 0.2 dB, about 0.5 dB after heavy rain. Drying is faster with wind (rho -0.76, 17 events, beats the null).
    - Dew is a real WAA (125 link-nights, median 1.7 dB, lasting hours, rho with a_link +0.26), outside the rain equation.
    - No rate-dependent radar error (the two binnings move in opposite directions = regression to the mean): the parked g(R) is CLOSED (not justified).
    - Paired links not testable; P is about 0.6 x ITU-R P.838 (size only).
  - **RETURN (2026-10-19):** the one allowed return to S6, triggered by S7 reason (d), is run as **S6b**. No further return is allowed:
    after S6b the candidate goes to S8 whatever the outcome.
  - **MISTAKE recorded (2026-10-19):** the 7.2 RMSE rule compared absolute errors across lengths (long links carry a larger path term); it
    failed for L > 3 km (1.18 x) while TRAIN shows the same (1.12 x). S8 keeps its own rule (TEST against VAL per class), which compares
    like with like.
  - **S6b** the return to S6 (added 2026-10-19; rules in config `wet_antenna.return_s6b`, written before the run): Part 1 diagnose the
    after-rain misfit cause by cause (radar end of rain late; faster drying; drying with the weather; F too high at heavy rain); Part 2
    refinements on TRAIN (M1 tau_f_dry free; M2 tau_f_dry = tau_0 exp(-k u) with the hourly wind; M3 only if the preceding rain and not
    the wind explains the drying), each with a fake check; Part 3 the box-A rule; Part 4 the three boxes and the VAL error budget; Part 5
    the candidate for S8 (sha256).
  - **S6b: done 2026-10-19** (`research/wet_antenna/outputs/reports/return_s6b.md`; one dated deviation: the projected cuts were undone). Results:
    - The after-rain misfit is a LEVEL at the end of the rain (-0.25 dB when the path becomes dry), growing with the preceding rain and in dry air; not the radar's timing, not a faster drying in general, not F at heavy rain.
    - M1 (tau_f_dry free, 3.7 min [1.7, 6.6]) passes every box-A rule (VAL after-rain error -0.24 / -0.17 -> -0.07 / -0.05 dB; RMSE -0.0024 dB; length test 4 / 6; 23 / 23 folds). M2 (wind, k +0.83 per SD) fails rules 3-4: box B. M3 (preceding rain): no effect, box C.
    - BUT with M1 the path term moves into the antenna: c_band halves and its bootstrap reaches 0 in every band: **c_band NOT identified with M1** (identified in S6); 3.7 min is an effective time for the radar input, not the gauge drying (10.3 / 8.8 min).
    - Error budget (VAL): dry minutes 30 %, during rain 68 % (between events only 3 %), after rain 2 %.
    - **Candidate for S8 by the rules: M1, `research/wet_antenna/outputs/products/return_s6b_candidate.json`, sha256 `06ac89c0abdd1bd49312607d2decad443716e6430910d5cea34f9a99c5879e1d`.** A learners' decision is needed before S8 (dated): M1 as written, or the S6 candidate (sha256 `36135658...d82b`) with the faster effective drying in box B. No further return to S6.
  - **S6b in short:**
    - The after-rain misfit is a LEVEL at the end of the rain (decay speed as the model, 0.87 [0.70, 1.06]; error at t = 0 -0.25 dB), larger after heavier rain and in drier air; not the radar timing.
    - F is not too high at heavy rain (VAL even low, 4 events).
    - Wind does not hold once the preceding rain is held fixed (M2 box B); M3 box C.
    - M1 (tau_f_dry free, 3.7 min) passed the written refinement rules but unbinds the path term (c_band halves, its interval reaches 0 in every band).
    - Error budget on VAL: dry minutes 30 %, during rain within events 65 %, between events 3 %, after rain 2 %.
  - **DECISION (2026-10-19, before any TEST read; a deviation from the outcome of the written S6b rule):** the candidate for S8 is
    **M0 = the S6 candidate** (sha256 `36135658c4362397754247faf4a744dc2e6e61b6737e9bcc112c83937519d82b`). M1 goes to box B: "a faster
    effective drying (3.7 min) removes the end-of-rain misfit but unbinds the path term; 3.7 min is an effective time for the 5 min radar
    input, not the gauge-measured drying time (10.3 min, VAL 8.8 min)". Reason: integrity rule 3 and the success definition (box-A numbers
    must be identified). **Lesson recorded:** a refinement rule must also require that the OTHER numbers of the model stay identified (the
    same lesson as M4c -> M4b-B).
  - **Plan (2026-10-19):** S8 (this step: freeze, then TEST once) -> S9 formulate. Part 2 of the project (the path law) comes after S9.
  - **S8** freeze (hashed), then TEST once (added 2026-10-19; rules in config `wet_antenna.test_s8`, written before TEST is read):
    comparison models on TRAIN (B0 no WAA, B1 constant WAA while it rains, B2 WAA without memory; B3 skipped: no literature WAA model is
    implemented in the repo); a freeze file with hashes, committed before TEST is read; pass rules P1-P5 on TEST; information items.
  - **S8: done 2026-10-19** (`research/wet_antenna/outputs/reports/test_s8.md`). Freeze `test_s8_freeze.json` committed in `bb30b5c` (code
    `761c6f1`); TEST read once. Two evaluation-code bugs found after the read and fixed (a radar / gauge ratio of 0 and 3 events without
    a g_e entry gave undefined rain; both now get the written fallback g_e = 1); recorded with the numbers before and after in config
    `test_s8.results`; nothing of the model or the rules changed. Results:
    - P2 pass: TEST RMSE candidate 0.692 dB against B0 0.862, B1 0.760, B2 0.727 (L < 1 km 0.562 against 0.838 / 0.644 / 0.601).
    - P3 pass (dry-path error L > 3 minus L < 1 km -0.04 [-0.15, +0.04] dB); P4 pass (TEST / VAL 1.12, classes 0.95-1.04 x); P5 pass (TEST links 1.07 x).
    - **P1 FAIL** (4 of 8 cells): the slopes s fit in 7 of 8, the intercepts I in 4 of 8 - the antenna level at 0.5-2 mm/h is higher on TEST
      (1.35-1.55 dB observed against 0.83-1.11 dB): the size of W at light rain, not its dependence on L (P3 passes).
    - Verdict: the box-A equation does not fully hold on TEST (P1). V1 (0.692 dB) beats V0 (0.724); the after-rain misfit is as on VAL
      (-0.18 / -0.15 dB at 5 / 10 min; M1 of box B removes it).
  - **S8 in short:**
    - P2-P5 pass: the candidate beats no WAA, constant WAA and WAA without memory on TEST (all links and L < 1 km); the dry-path error is equal across lengths; TEST / VAL is even across lengths; new links 1.07 x.
    - P1 fails (4 of 8 cells): the slopes match; the intercept (the WAA level) at 0.5-2 mm/h is higher on TEST (about 1.4-1.5 dB against 0.8-1.1 dB in the model) -> **box C**: the WAA level in light rain varies between storms more than the equation captures.
    - The after-rain misfit is 6 % of the TEST error (M1, box B, removes it).
    - Two evaluation bugs (g_e = 0 or missing) were fixed after reading TEST with the written fallback g_e = 1.
  - **S8b** (added 2026-10-19, after TEST was read): a **post-freeze evaluation of a pre-freeze model** - the S6 plug-back fit "F and P"
    (F = A_band (R / 5)^b, P = c_band R, everything else as S6; the "closed F" of S7). Nothing is refitted on TEST and the S8 verdict does
    not change. The fit: `research/.cache/wet_antenna/tighten_s6/plug_shared_FP_F_and_P.pkl`, sha256
    `1d8163cd8b84aafc85f59c6c223f44b33923278265e3e40819923d8076c8e0b9`, modified 2026-10-18 12:45:04 (+0300): BEFORE the S8 freeze
    (`test_s8_freeze.json`, 2026-10-19 18:14).
  - **DECISION RULE (2026-10-19, after TEST was read, before this evaluation):** the closed form is adopted as THE written form of F if and
    only if, on TEST, (1) its verdicts on P1-P5 are the same as the table F's, and (2) the paired event-bootstrap 90 % interval of its RMSE
    minus the table F's RMSE lies within [-0.01, +0.01] dB (all link-events and L < 1 km). Otherwise the table F stays the equation and
    the closed form is reported as an approximate description. Either outcome is reported.
  - **S8b: done 2026-10-19** (`research/wet_antenna/outputs/reports/test_s8b.md`; post-freeze evaluation of a pre-freeze model):
    - The same P1-P5 verdicts as the table F on TEST; the paired RMSE difference +0.0025 [-0.0065, +0.0085] dB (L < 1 km +0.0016 [-0.0082, +0.0078]); on TRAIN all 9 numbers identified, inside the table F's interval at 21 / 21 knots.
    - **By the rule the closed form is ADOPTED as the written form of F**: F_band(R) = A_band (R / 5)^b, b = 0.54 [0.47, 0.63], A = 2.57 / 2.36 / 2.77 dB (28 / 32 / 38 GHz), with P = c_band R (0.105 / 0.141 / 0.194 dB/km per mm/h); valid 0-40 mm/h only (no saturation; +2.3-2.7 dB above the table F at 40 mm/h).
  - **S8b in short:** the closed form F = A_band (R / 5)^b was fitted on TRAIN in S6, before the S8 freeze. It was run on TEST after S8
    without refitting, under a rule written first: the same P1-P5 verdicts, RMSE difference +0.0025 [-0.0065, +0.0085] dB. Adopted as the
    written form of F: b = 0.542 [0.472, 0.626], A = 2.57 / 2.36 / 2.77 dB, c = 0.105 / 0.141 / 0.194 dB/km per mm/h; valid 0-40 mm/h.
  - **S9** formulate Part 1 (added 2026-10-19): the official equation (`outputs/products/waa_equation_final.json` and
    `research/wet_antenna/waa_equation.py` with a unit test), a sentence-by-sentence check of the overview against the repo, the short
    overview and the full (school) version in `docs/overview/`, and `research/wet_antenna/FINAL.md` (the two documents were merged into
    `docs/WAA_research_overview.pdf` afterwards, see below). No new analysis, no refitting, no
    new reading of TEST.
  - **Plan after S9 (2026-10-19):** **R1** "does the WAA correction improve link rain estimates against the gauges" (its own written plan,
    and a decision about TEST, which was already read once), then **Part 2** (the path law).
  - **S9** formulate: the story, the assumptions with their checks, the data funnel, the core figures, FINAL.md. The last step.
  - **S9: done 2026-10-19** (`research/wet_antenna/outputs/reports/formulate_s9.md`; no new analysis, no refitting, no new reading of TEST):
    - The official equation `research/wet_antenna/outputs/products/waa_equation_final.json` (sha256 `b06ebc8b75867171368e8ea0540bb011014e9615d3dde74f8559b524a25be9cd`; the fit sha256 `1d8163cd...e0b9`) and `research/wet_antenna/waa_equation.py`; the unit test reproduces 3 TEST link-events of the S8 cache within 9e-16 dB.
    - The overview checked statement by statement (`formulate_s9_check.md`; its LaTeX source was missing, so it was rebuilt from the PDF): 113 statements, 21 corrected, 2 unsupported and reworded (a 25 % band difference, the cause of the remaining error), reference [2] title corrected.
    - `docs/overview/waa_overview.pdf` (short, 12 pages) and `docs/overview/waa_part1_full.pdf` (full school version, 19 pages), 12 colour-free figures from saved tables (`docs/overview/figs.py`); `research/wet_antenna/FINAL.md`.
    - **Change (2026-10-19, decided by the learners after S9):** one WAA document only. The corrected short overview was merged into
      `docs/WAA_research_overview.pdf` (the first overview's look and wording, the S9 corrections, real-data figures; source
      `docs/WAA_research_overview.tex`, figures by `docs/WAA_research_overview_figures/figs.py`). `docs/overview/` and the full school
      version were deleted.
  - **S8c: done 2026-10-28** (`research/wet_antenna/outputs/reports/literature_s8c.md`; a post-freeze comparison, our equation not refitted): Schleiss et al. 2013 and Pastorek et al. 2022 (KR-alt), pycomlink 0.6.0, fitted on TRAIN only (Schleiss driven by the radar wet flag, an exponential drying in place of the measured attenuation).
    On TEST (verdict = version b, fitted on TRAIN): Schleiss b no difference (+0.009 [-0.010, +0.039] dB; L < 1 km +0.012 [-0.009, +0.044]), with a path term 1.7-1.8 x ours and too much water after rain; Pastorek b ours better (+0.026 [+0.015, +0.038]; L < 1 km +0.031 [+0.017, +0.046]). Version a (Schleiss a = pycomlink defaults plus an assumed drying time of 15 min): ours better. Added to `docs/WAA_research_overview.pdf` (Table 7, Figure 7).
  - **Overview v2: done 2026-11-02** (`research/wet_antenna/outputs/reports/overview_v2.md`): `docs/WAA_research_overview.pdf` rewritten as a two-column paper (v1 kept as `docs/WAA_research_overview_v1.tex`); all WAA models compared in one section, with a new post-freeze consistency table (`comparison_all_models.py`, nothing refitted): ours ties Schleiss (b) in accuracy and is the only model whose a_link does not depend on L (rho +0.09 [-0.03, +0.20]; column added after the first run, dated), but its light-rain intercept is the weakest (4 of 8 cells) and no model matches the after-rain WAA.
  - **PART 1 CLOSED (2026-10-19).** Next: **R1** (its own written plan, and a decision about TEST), then **Part 2** (the path law).

## Integrity rules (from S6 on)

1. Every rule and threshold is written in config before the run.
2. A fake-data check runs before real numbers are trusted.
3. A quantity is reported as a result only if it is identified (fake recovery within +-20 %, interval not at a bound).
4. Every deviation from a written rule is dated, with its reason.
5. Failed checks are reported, never dropped. S7 may send the work back to S6 at most once, and that is recorded.
6. VAL has been used for decisions, so VAL numbers are labelled "used (optimistic)".
7. TEST is read once, in S8, after a hashed freeze; nothing changes after it.
8. ITU-R P.838 never corrects data.

## What the model claims - three boxes (2026-10-19)

- **A - in the equation:** strong, physical, identified, holds out of sample.
- **B - identified, outside the equation:** its size, conditions and evidence are known, but it is not clean or general enough to model.
- **C - unexplained:** part of the error budget.

**Success** means: box A explains most of the behaviour on unseen events and links, its numbers are identified and bounded, it beats
simpler models, and it does not depend on L; boxes B and C are named and sized. **Forcing a B effect into the equation is a failure, not
a success.**

## Findings that stand (after S3)

- **Time constants with the gauge input** (TRAIN S2b + VAL S3.4): wetting 2-3 min, drying 9-13 min (M3 2.9 / 12.8 TRAIN,
  1.8 / 10.8 VAL; M4b fast store 2.6 / 10.3 TRAIN, 2.0 / 8.8 VAL).
- **The F / P split from links < 1 km alone is not identifiable** (S2b: P of M3 halves from S<1 to S<3; S3: neither fit
  matches the model-free slope).
- **The frame W + L P(R) fitted on all lengths reproduces the model-free I(R), s(R) on TRAIN** (S3b: M3-B, M4b-B 10 / 10 cells).
- **Slow radar wetting = input smoothing** (S2 fake check Ta; S3b: tau_wet stays 13-56 min with P pinned, 2-3 min with the gauge).
- **Fast (gauge) timing reduces the peak deficit** (S3b, VAL: -0.47 / -0.51 dB with gauge timing against -0.87 / -0.75 dB free).
- **The VAL gap is the rain reference, not the antennas** (S4.1-S4.3): the radar over-reads the large VAL events ~1.7 x against the
  City gauges (TRAIN ~0.9 x); dividing by that ratio closes the gap. WAA v3 = M4b-B unchanged.
- **Evaluation convention from now on:** every error is given on the raw radar (V0) AND on the radar / g_e (V1, the main reference),
  with the pooled ratio per set (V3) as a sensitivity; the gauge links with the gauge input are the independent check.
- **Model-free I(R) and s(R)** (S3.1, TRAIN events, all link lengths, steady rain; (A - c) = I + s L per band and rain class,
  Huber, 90 % event-bootstrap intervals; `research/wet_antenna/outputs/tables/system_id_s3_length_fits.csv`):

| band | rain class (mm/h) | R of the cell | I (dB) | s (dB/km) | L0 = I / s (km) | links (under 1 km) |
|---|---|---|---|---|---|---|
| 28 GHz | 0.5-1 | 0.77 | 1.02 [0.84, 1.26] | 0.08 [0.05, 0.12] | 12.0 [7.6, 23.8] | 113 (17) |
| 28 GHz | 1-2 | 1.39 | 1.33 [1.07, 1.58] | 0.18 [0.14, 0.22] | 7.5 [5.7, 9.8] | 113 (17) |
| 28 GHz | 2-5 | 2.96 | 1.86 [1.52, 2.12] | 0.34 [0.28, 0.42] | 5.5 [4.1, 6.9] | 113 (17) |
| 28 GHz | 5-10 | 6.11 | 2.45 [2.08, 2.78] | 0.82 [0.72, 0.93] | 3.0 [2.3, 3.7] | 113 (17) |
| 28 GHz | 10-20 | 12.44 | 2.82 [1.64, 5.38] | 1.44 [0.68, 2.13] | 2.0 [0.8, 6.3] | 97 (15) |
| 32 GHz | 0.5-1 | 0.79 | 0.90 [0.68, 1.18] | 0.12 [0.05, 0.18] | 7.2 [3.8, 21.1] | 19 (0) |
| 32 GHz | 1-2 | 1.35 | 1.13 [0.80, 1.48] | 0.25 [0.16, 0.36] | 4.5 [2.7, 8.1] | 19 (0) |
| 32 GHz | 2-5 | 2.89 | 1.86 [1.39, 2.27] | 0.39 [0.22, 0.59] | 4.8 [2.8, 10.2] | 19 (0) |
| 32 GHz | 5-10 | 5.83 | 2.41 [1.96, 3.22] | 0.86 [0.41, 1.12] | 2.8 [1.8, 7.2] | 19 (0) |
| 32 GHz | 10-20 | 11.78 | 3.56 [-5.49, 5.09] | 0.74 [0.35, 6.13] | 4.8 [-0.9, 14.6] | 16 (0) |
| 38 GHz | 0.5-1 | 0.79 | 1.10 [0.91, 1.36] | 0.16 [0.11, 0.20] | 7.0 [5.1, 10.9] | 74 (24) |
| 38 GHz | 1-2 | 1.38 | 1.47 [1.13, 1.84] | 0.36 [0.29, 0.42] | 4.1 [2.9, 5.9] | 74 (24) |
| 38 GHz | 2-5 | 3.04 | 2.08 [1.74, 2.39] | 0.64 [0.47, 0.78] | 3.3 [2.5, 4.4] | 74 (24) |
| 38 GHz | 5-10 | 6.07 | 2.35 [2.00, 2.81] | 1.56 [1.20, 1.70] | 1.5 [1.3, 2.2] | 74 (24) |
| 38 GHz | 10-20 | 12.00 | 2.44 [2.02, 3.99] | 2.66 [2.24, 3.09] | 0.9 [0.7, 1.7] | 71 (23) |

### WAA v3 (provisional, 2026-10-15) = M4b-B of S3b

```
A(t) = W(t) + L * P_band(R(t - d)) + c,   W = W_f + W_s     (radar input, TRAIN events of ALL link lengths)
W_f -> phi a_link F_band(R)        tau_f_wet 2.6 / tau_f_dry 10.3 min   (fixed: gauge timing of S2b)
W_s -> (1 - phi) a_link F_band(R)  tau_s 240 min, phi 0.72               (fixed)
lag d                     +1 min
F[28 GHz](0.5 / 2 / 5 / 10)   0.98 / 1.63 / 2.62 / 3.38 dB      P[28 GHz](0.5 / 2 / 5 / 10)   0.01 / 0.24 / 0.63 / 1.14 dB/km
F[32 GHz](0.5 / 2 / 5 / 10)   0.90 / 1.70 / 2.94 / 3.71 dB      P[32 GHz](0.5 / 2 / 5 / 10)   0.01 / 0.33 / 0.70 / 0.97 dB/km
F[38 GHz](0.5 / 2 / 5 / 10)   0.96 / 2.02 / 3.00 / 3.48 dB      P[38 GHz](0.5 / 2 / 5 / 10)   0.02 / 0.47 / 1.09 / 2.00 dB/km
a_link 10-90 %            0.73 - 1.21
(90 % intervals in system_id_s3b_parameters.csv: 30 warm-started bootstrap refits, too narrow)
Limits: event-to-event level differences (VAL attenuates 0.6-0.7 x per radar mm/h; fails the VAL length test);
        the slow store tau_s sits at the 240 min bound and the tail after the rain is too high (after-rain integral 1.2-1.4 x).
```

## Open questions (parked, for the final report)

- M5 / naming with a boundedness condition (h rising, g not growing with W faster than h).
- Why the gauge and the radar give different drying times (10.6 against 7.1 min on the same points). **Radar vs gauge wetting: SETTLED (S3b): input smoothing.**
- The peaks: every model is 0.7-0.9 dB too low. **Mostly timing: fixed with gauge timing (S3b, -0.47 dB).**
- Dew / night rises of the level.
- What drives the slow leftover after the rain (water or a level change).
- (S2b) The slow leftover is not described by any rung: tau_s of M4b sits at the 240 min bound; a store fed by F
  overfills on long events. Candidate for S5: a store fed by something else (e.g. the rain amount).
- (S2b) Peaks 0.7-0.8 dB too low with both inputs (radar and gauge). **Mostly timing: fixed with gauge timing (S3b).**
- (S2b) Drying with the radar 25 min against 13 min with the gauge (M3).
- (S2b) Bootstrap intervals too narrow (warm-started 30-iteration refits); the final numbers in S5 need full refits.
- (S2b) The judge depends on the fold cap (longer folds raise the held-out RMSE by 0.001-0.003 dB, the size of the
  differences between rungs).
- (S3) The antenna level is lower on the VAL events (I 0.8-1.1 dB against 1.0-2.45 dB on TRAIN; also the gauge-input F
  of M3): event dependence -> S4.
- (S3) The intercept I holds everything that does not grow with L (a radar bias common to all links, baseline drift, dew),
  not only the antenna.
- (S3) Uneven rain on long links may lower s (the radar path mean on a 2 x 2 km grid; b < 1).
- (S3) 32 GHz has no link < 1 km: its intercepts are extrapolations.
- (S3b) VAL events attenuate ~0.6-0.7 x per mm/h of radar rain at the antenna and on the path (all fits, M1 too): radar error,
  drop sizes or antenna state per event? -> S4.
- (S4.2) Warm events over-predicted more (residual -0.15 vs -0.02 dB): antenna, rain type or radar?
- (S4.3) Does alpha_e depend on rain burstiness / humidity? (hints only, 22 heavy events)
- (S4.3) The radar error may depend on the rain rate (28 GHz 2-5 mm/h stays low after V1).
- (S5a) Rain-rate dependent radar error (peaks under V1; 28 GHz 2-5 mm/h).
- (S5a) F shape: saturation vs slow log rise.
- g(R) rain-rate dependent correction from the gauges - decision point in S7.6. **CLOSED (S7.6, 2026-10-19): not justified** (no rate dependence; the two binnings move apart).

## Status

The wet-antenna steps have two series. **T0-T18** = the closed drying-tail exploration (renamed on 2026-10-13; the old
number is given in brackets, and scripts, outputs, reports and caches keep their names). **S1-S5** = the reboot.
Old number of each T step: T0 = 0, T1 = 1, T2 = 2, T3 = 2b, T4 = 2c, T5 = 3a and 3a-fix, T6 = 3b-1, T7 = 3b-1b, T8 = 3b-1c, T9 = 3b-2, T10 = 3b-3, T11 = 3b-3b, T12 = 3b-3c, T13 = 3b-4, T14 = 4, T15 = 4b, T16 = 4c, T17 = 4c2, T18 = 7a.
The old planned steps that were never run (5, 6, 3c, 4d, 7, 8a-8e, 9) are CANCELLED; each row says which S step takes
its place.

| Part | Step | What | Status | Result |
|---|---|---|---|---|
| Part 1 (closed): 1 Prepare data | T0 (was 0) | Inventory of inputs | done | A(t) = TL - frozen baseline v2 |
| Part 1 (closed): 1 | T1 (was 1) | Extract all drying tails | done | 20,355 tails |
| Part 1 (closed): 1 | T2 (was 2) | Selection and fast-drop test | done | 959 main tails; fast drop not proven to be WAA |
| Part 1 (closed): 1 | T3 (was 2b) | Fixes: t0_trace, strict/lenient ends, noise exclusion, level flags | done |  |
| Part 1 (closed): 1 | T4 (was 2c) | Frozen tail set v1 | done | 429 train / 123 val tails, 159 links, 24 events |
| Part 1 (closed): 2 Practise on fake data | T5 (was 3a, 3a-fix) | Fake-tail simulator with realistic noise and spread | done |  |
| Part 1 (closed): 2 | T6 (was 3b-1) | Fitting tool, one tail at a time | done | too weak: only exponential vs slower |
| Part 1 (closed): 2 | T7 (was 3b-1b) | Pooled fitting tool (all tails together) | done | pre-set rule met |
| Part 1 (closed): 2 | T8 (was 3b-1c) | Stress test of the pooled fit | done | PARTIAL: m1, m4, m6 hold; F3/F5 one family |
| Part 1 (closed): 2 | T9 (was 3b-2) | SINDy directly on raw tails | done | FAIL: returns "A" for almost every truth |
| Part 1 (closed): 2 | T10 (was 3b-3) | Two jumps: free-curve pooled fit -> drying curve, then SINDy -> equation | done | PARTIAL: curve and equation right above ~0.35 dB (mixing trap gone there), undetermined below; m1 fails the pre-set rule; the switch step and the two layers are not visible in the curve |
| Part 1 (closed): 2 | T11 (was 3b-3b) | Confirmation run of the two jumps on fresh fake data, judged only above 0.35 dB, rule fixed before the run | done | PARTIAL: m1, m3, m5 PASS (curve, shape, equation); the m4 warning sign fails (2 of 5) |
| Part 1 (closed): 2 | T12 (was 3b-3c) | Jump 2 on TRUE curves (seconds); Jump 1 that may draw a step at a fitted W_c; Jump 3 (refit and choose) | done | PARTIAL: Jump 2 alone works (E-SINDy, PySR); the step tool draws a step for EVERY truth; Jump 3 and the bootstrap band fail |
| Part 1 (closed): 2 | T13 (was 3b-4) | Diagnosis only: is the 0.3 dB RSL quantization the cause of the false step? m1 simulated with no / 0.3 / 0.1 dB quantization, step fit of T12 (was 3b-3c) part B | done | NOT the cause: without quantization the step is still kept in 5 of 5 (BIC gain 40-130); cause open |
| Part 1 (closed): 3 Find the law on real data | T14 (was 4) | The confirmed setup (T11 (was 3b-3b)) + menu fits on the 429 real tails, four sensitivity sets | done 2026-10-12 | NO law identified: the real curve FALLS with W (slope -0.87; -0.05 without the rise-heavy events), no formula, menu not decisive (F1 as run; like for like F6 by 21 BIC); a W0 fault found in the menu comparison; F3 / F5 at their exponential limit in all 7 sets |
| Part 1 (closed): 3 | T15 (was 4b) | Why do some real tails stall or rise: contamination (H1) or antenna physics (H2)? Classes, outside evidence, antenna tests, two fixes of the pooled fit | done 2026-10-12 (set to in progress and finished the same day) | NEITHER rule met. The stall is a property of the event, not of the level, the link or the path length: in decaying tails k does not fall with A (+0.28), so h(W) needs no falling k; rain references show less, not more, near stalled tails (path radar dry by construction; humidity missing for the two worst events). Classes: 240 decay, 64 wait then decay, 36 rise, 89 flat / slow. Fixes: W0 bound and noise multiplier (off by default; fake check 4 of 6, not passed) |
| Part 1 (closed): 3 | T16 (was 4c) | Phase A: what a drying tail contains. A1 daily cycle / dew vs water, A2 antenna vs path rain (A vs L per minute), A3 three start rules; proposes the rules of tail set v2 (no law fit) | done 2026-10-12 | A1: the rises are mostly NOT drift (24 of 36 above the same sub-link's dry days at that clock time; night-time rises 10 of 18), so no baseline v3; dew or drizzle not decided (humidity missing for 22 of the 36). A2: NO POWER (beta is not above 0 before t0_trace either; the radar expects only 0.08-0.13 dB of path rain there), antenna vs path stays open. A3: the references within 5 km are dry 25 min after t0_trace (median), when 349 of 429 tails are below 1 dB; k(A) falls under every start rule (-0.80, -1.85, -0.56). Proposal: keep t0_trace, store the strict-dry and onset times per tail, model the start per tail in 4d |
| Part 1 (closed): 3 | T17 (was 4c2) | Quick check, data only: once no reference within 5 km shows rain (from S1 on), how fast does A fall? Lower A_start thresholds, a placebo start with the same delays, the same level before and after S1 | done 2026-10-12 | SAME speed by the rule (threshold 0.4 dB): t_half 40 min from S1 (760 tails), 37 min from the placebo start, 18 min from S0; S1 - placebo +3 min (-2..+7). The slowness after S1 comes from starting late, not from the rain nearby having ended (diluted contrast: 52 % of the placebo starts are after the own S1; split 33 against 43 min). A tail has two phases in time: fast for 15-20 min after t0_trace, then t_half about 40 min; at the same level k is 4-7 times lower later. Whether the fast phase is the antenna or rain on the path is NOT decided: Step 7 |
| Part 1 (closed): 3 | 4d (old plan) | Phase B: B0 the simulator gets what the real tails showed (start delays, rises, daily cycle, noise 2.6-4.3 x) and the tools must pass on that; B1 the drying law on tail set v2, with a per-tail start delay as a model option | CANCELLED 2026-10-13 (was: planned (moved after Step 7 on 2026-10-12: informed by the whole-event model)) | see the Reboot: S1 (drying is one rung of the model ladder) |
| Part 1 (closed): 3 | 5 (old plan) | SINDy on raw real tails | CANCELLED 2026-10-13 (was: replaced) | see the Reboot: S2. Was: replaced by T10 (was 3b-3) + 4 |
| Part 1 (closed): 3 | 6 (old plan) | Phase C: what the start delay and the drying speed depend on, per event (night / sun, temperature, humidity, wind; also frequency, antenna) | CANCELLED 2026-10-13 (was: planned (after 4d)) | see the Reboot: S4. Was: expect weak power: about 24 events, one hourly station, RH missing for two events |
| Part 1 (closed): 3 | 3c (old plan) | Fake full events with rain, to test the wetting tools | CANCELLED 2026-10-13 (was: next after T17 (was 4c2) (Phase D moved forward on 2026-10-12)) | see the Reboot: S1 (its fake check). Was: needed for Step 7. Decided: WHOLE events on very short links, radar rain as input |
| Part 1 (closed): 3 | T18 (was 7a) | Very short links: is there enough data, and do they show the fast first phase? Inventory per length class; the two phases by length (SHORT against LONG, matched on A_start); a first look at whole events below 1 km (no law fit) | done 2026-10-13 | UNDECIDED by the rule (power condition fails), but the fast drop does not grow with L: 0.27 / 0.30 / 0.30 / 0.30 / 0.28 dB from < 0.6 km to >= 4 km; matched SHORT / LONG ratio 0.94 (0.84..1.15), below 0.6 km 1.01 (0.78..1.16), against 0.2-0.26 if it were path rain. The radar-expected path rain before a tail is only 0.02 / 0.08 dB, which is why the rule has no power. Enough data below 0.6 km by the rule (18 links, 37 events, 685 pairs; below 1 km 41 links, 1663 pairs); only 28 and 38 GHz below 1 km. Whole events: several showers over 9 h; A is at 19 % of the last shower's peak when link_state ends the rain; no saturation level seen; a rise before the radar rain is rare (4.2 % against 1.4 % placebo) |
| Part 1 (closed): 3 | 7 (old plan) | Very short links, WHOLE events with radar rain as input: the model sees the rain end by itself, so nobody has to guess where drying starts; wetting g(R, W) and drying h(W) together | CANCELLED 2026-10-13 (was: planned (after 3c; before 4d)) | see the Reboot: S1 and S2. Was: the test of Q2 (antenna or rain) with enough contrast in L |
| Part 1 (closed): 4 Prove it | 8a (old plan) | Drying law predicts unseen tails (VAL, leave-one-event-out) vs literature exponential | CANCELLED 2026-10-13 (was: planned) | see the Reboot: S3 |
| Part 1 (closed): 4 | 8b (old plan) | Full law predicts unseen events (short links) vs Schleiss | CANCELLED 2026-10-13 (was: planned) | see the Reboot: S3 |
| Part 1 (closed): 4 | 8c (old plan) | Longer links: A minus predicted WAA behaves like pure rain (f = L + L0 test) | CANCELLED 2026-10-13 (was: planned) | see the Reboot: S3 |
| Part 1 (closed): 4 | 8d (old plan) | Consistency: curve vs fitted parameters, sub-link pairs, event resampling | CANCELLED 2026-10-13 (was: planned) | see the Reboot: S3 |
| Part 1 (closed): 4 | 8e (old plan) | Dew events dry by the same h(W) | CANCELLED 2026-10-13 (was: planned) | see the Reboot: S3 |
| Part 1 (closed): 5 Final | 9 (old plan) | Final equation, parameters, plots, limits, handoff to Stage 2 | CANCELLED 2026-10-13 (was: planned) | see the Reboot: S5 |
| Part 2 (reboot) | S1 | Fake check of the fitting, then the model ladder M0-M4 on the real short links (whole events, radar rain in, A out) | done 2026-10-13 | fake check PASS (truth M3 chosen 3 of 3, truth M4 chosen 3 of 3; time constants, phi and d within 10 %, tau_s x1.09). Held-out RMSE of A, 42 events: M0 0.417, M1 0.403, M2 0.394, M3 0.395, M4 0.390 dB; best M4, **chosen by the stop rule: M4** (M3 - M4 = +0.004, 0.0015..0.0067), the TOP rung, so the ladder did not stop. Dynamic, not static (M1 - M4 = +0.013, 0.007..0.020), but the skill over M0 is only 6.4 % (4.7..8.2). M4: phi 0.71 (0.53..0.79), tau_f_dry 8.3 min (5.3..9.9; 7.4-10.6 in every sensitivity run), tau_f_wet 17 min (10..65) with the radar but 2.4 min with the City gauge as input (gauge run: RMSE 0.339 against 0.402 on the same points), tau_s 513 min (210..600 = the bound: not identified; the integral of A after the rain is 2.4 x the observed one). P(10 mm/h) 1.5 / 2.5 dB/km (28 / 38 GHz; ITU size 1.8 / 2.8) but 1.1-4.6 / 2.2-6.1 over the rungs: the split between P and F is not firm with L below 1 km. Peaks too low by 0.86 dB. Report: `research/wet_antenna/outputs/reports/system_id_ladder.md`. (Deliverable was: the leave-one-event-out error of M0-M4 with intervals, and the rung chosen by the stop rule) |
| Part 2 (reboot) | S2 | Discovery rung M5: a flexible G(R, W), then named with E-SINDy / PySR; also the links up to 3 km (length rule), the City-gauge input (gauge rule) and the slow part (water or drift) | done 2026-10-14 | fake checks: Ta (input) PASS, Tb (drift / slow store) PASS, **Tc (F and P on S<1 vs S<3) FAIL, Td (M5 on M3 data) FAIL**. Held-out RMSE on S<3 (51 events, 6653 link-events): M1 0.458, M3 0.448, M3D 0.448, M4 0.444, M5 0.442 dB; **chosen M4** (best M5; M4 - M5 +0.0024, -0.0002..+0.0052); every rung above M3 adds only a few thousandths of a dB. Length rule DISAGREE (M4 tau_f_dry 5.5 min on S<3 against 8.3 on S<1; not validated): S<1 stays the main set; no residual trend with L. Gauge input: wetting 2.7 min (M4; 3.0 for M3), drying 10.6 min against 7.1 with the radar (do not agree); the slow radar wetting of S1 is the input (Ta: a 5 min noisy input makes the wetting 2.1 x slower). Slow part: SLOW RISE, FORM UNDECIDED (observed after-rain e-folding 76-122 min; M3 too fast, M4 / M5 too slow). M5: best by a hair, but not identified (Td) and **no formula named** (E-SINDy error >= 0.45 on g; the E-SINDy formula makes W run away, M5n could not be refitted). Found from Tc: **an optimiser trap in P and F** (P near 0 at R <= 2 mm/h; a refit from the true P reaches a lower loss). Report: `research/wet_antenna/outputs/reports/system_id_s2.md`. (Deliverable was: a named equation, or "M5 is not better than the chosen rung") |
| Part 2 (reboot) | S2b | Repair and finish the S2 model: P and F with plain non-negative increments (the trap), the fake check Tc again and the length rule again; M4b = M4 with tau_s bounded to 30-240 min (fake check, ladder M1 < M3 < M3D < M4b, after-rain behaviour); longer fits (8 starts, 5 x the S2 cap, folds 3 x) | done 2026-10-14 | cuts by the time rule: S<2, the gauge-links radar comparison, fake replicates 3 -> 2. **Tc PASS** (all 30 F and P knots within 15 %, P(2) / truth 1.01-1.04), **Te PASS** (M4b chosen 2/2, tau_s 96 for 90). On the real radar data S2 was hardly trapped (training loss -0.1 to -10; numbers within a few %); the trap showed in the GAUGE-input F at light rain (F(0.5) 0.14 -> 0.64-0.98 dB), so v1's "still rising at 5 mm/h" was mostly the trap. Held-out RMSE S<1: M1 0.403, M3 0.396, M3D 0.396, M4b 0.393; **chosen M3 = WAA v2** (M4b on S<3, M3D with the gauge). tau_s of M4b at the bound 240 on every set; the slow leftover after the rain is still not described (M4b 2.07 x, M3 1.05 x the observed integral on S<1). Length rule AGREE (on M3). 13 of 15 real fits converged (not M3D). WAA v2 (M3): wetting 2.9 [2.5, 3.6] min (gauge), drying 12.8 (gauge) / 24.7 (radar) min, F about 2-3.5 dB by 5 mm/h; limits: the F / P split depends on model and input, the leftover, the peaks. Report: `research/wet_antenna/outputs/reports/system_id_s2b.md`. (Deliverable was: WAA v2, or v1 with the limits that remain) |
| Part 2 (reboot) | S3 | Verification: VAL events, both directions of a link, long links (handoff to Stage 2); sub-steps S3.1-S3.4 (see Remaining plan) | done 2026-10-15 | **Neither M3 nor M4b passes the model-free length test** (S3.1; M3 P 3-8 x the slope s, M4b P near 0 at light rain and level 1.2-1.5 x I); I saturates at 2-3 dB, L0 = I / s falls from about 12 to 2 km (no constant L0). VAL overall RMSE as TRAIN (M4b 0.360, M3 0.365 dB below 1 km), worse during rain. Directions: peaks agree (0.96), integrals noisy (0.34). Gauge on VAL: wetting 1.8-2.0, drying 8.8-10.8 min. WAA v3 candidate: neither. Report: `research/wet_antenna/outputs/reports/system_id_s3.md`. (Deliverable was: the error on unseen data and the f = L + L0 check on long links) |
| Part 2 (reboot) | S3b | Fit the S2 frame (M3, M4b; M1 for scale) on the TRAIN events of ALL link lengths so that P is pinned by L; variant A free, variant B with the gauge time constants fixed; judged on VAL (length test through the model, errors per length class, residual vs L) | done 2026-10-15 | TRAIN length test reproduced with all lengths (M3-B, M4b-B 10 / 10 cells); **on VAL every fit fails** (0-1 / 8): VAL attenuates ~0.6-0.7 x per mm/h at the antenna and on the path (event-set difference, M1 too). Residual flat with L; VAL RMSE 0.398-0.407 dB; gauge timing gives better peaks; radar wetting stays slow with P pinned (input smoothing). **WAA v3 by the S3b rule: none**; provisional WAA v3 = M4b-B (see the box under "Findings that stand"). Report: `research/wet_antenna/outputs/reports/system_id_s3b.md`. (Deliverable was: WAA v3, or what fails) |
| Part 2 (reboot) | S4.1 | Why do events differ: event facts and the split, random-split test of the TRAIN / VAL gap, per-event rain scale k_e and antenna factor alpha_e on the frozen M4b-B (fake check first), radar vs gauge per event | done 2026-10-15 | VAL is rare (0-1.5 % of random 18-event TRAIN halves as low); **fake check of the factors FAIL** (63 %; heavy events 16 / 16), factors not interpreted; factor-free: radar / gauge 1.71 on VAL against 0.89 on TRAIN, the size of the gap (rain reference most likely). Report: `research/wet_antenna/outputs/reports/conditions_s4.md`. (Deliverable was: whether the VAL gap is ordinary event spread, rain-input scale or antenna) |
| Part 2 (reboot) | S4.2 | Conditions behind the factors (descriptive): the condition data in the repo per event against k_e, alpha_e | done 2026-10-15 | no condition passes BH (28 tests; not interpreted after the fake FAIL); residual of M4b-B -0.08 dB at every band, day = night, lower in the warmest events. Candidate: keep M4b-B as is. (Deliverable was: which conditions go with the factors) |
| Part 2 (reboot) | S4.3 | Is the VAL gap the radar? The VAL length test with the radar rescaled by the event's radar / gauge ratio (City gauges only; frozen M4b-B, M3-B check, no refit); per-event factors on heavy events only after a new heavy-event fake check | done 2026-10-17 | **Explained** (V1: VAL I 0.99, s 0.95, 5 / 6 cells; TRAIN not broken); heavy fake PASS 98 %, k_e ~ g_e (rho 0.64), no condition. Keep M4b-B; S5 uses the g_e-rescaled radar as main reference. Report: `research/wet_antenna/outputs/reports/conditions_s4c.md`. (Deliverable was: whether the gauge rescaling explains the VAL gap without breaking TRAIN; the heavy-event factors if their fake check passes) |
| Part 2 (reboot) | S5a | Refit WAA v3 (M4b-B) on gauge-corrected rain (cold bootstrap), slow-leftover candidate M4c (fake check first), closed forms for F(R), P(R) with the plug-back check | done 2026-10-18 | Refit on V1: VAL length test passes (5 / 6), V0 fails; M4c replaces M4b-B by the rule but its slow store is not identified (tau_s_wet at bound); P = c R / c R^0.88, F = F_max R / (R + R0); plug-back holds: **WAA v4 (closed form)**. Report: `research/wet_antenna/outputs/reports/final_s5a.md`. (Deliverable was: WAA v4 (closed form) or the free curves kept; the S5b rules) |
| Part 2 (reboot) | S5b | TEST events and links, read once with the S5a model and rules; final write-up | NOT run (2026-10-18) | Its parts move into S6, S8, S9 |
| Part 2 (reboot) | S6 | Tighten the equation: M4b-B on V1, data funnel, well-determined parametrisation (shared shape), cold bootstrap + profile intervals, fake identifiability, closed forms (physics and PySR compete), plug-back | done 2026-10-18 | Shared shapes kept; P = c_band R; F free shared knots (closed form failed the length test at exactly 2 / 3); F level set by the a_link prior; P < 1 mm/h not identified. Report: `research/wet_antenna/outputs/reports/tighten_s6.md` |
| Part 2 (reboot) | S7 | Validate without TEST: W independent of L, short links on VAL, long links (dry-path moments, paired links, short-only / long-only intercepts), dew, leftovers vs metadata and weather (PySR), radar vs gauges by rate | done 2026-10-19 | W independent of L (a_link, W height, dry-path W, intercepts pass); no rate-dependent radar error; paired links not testable; RMSE rule fails for L > 3 km; model 0.5 dB high after heavy rain on dry paths; drying time falls with the wind -> RETURN TO S6 (2026-10-19, run as S6b). Report: `research/wet_antenna/outputs/reports/validate_s7.md` |
| Part 2 (reboot) | S6b | The one return to S6: the after-rain misfit cause by cause; tau_f_dry free (M1), with the wind (M2), with the preceding rain (M3); the box-A rule; the three boxes; the error budget | done 2026-10-19 | Misfit = a level at the end of rain; M1 passes box A but makes c_band NOT identified; M2 (wind) box B; M3 box C. Candidate by the rules: M1, sha256 06ac89c0...; a learners' decision (M1 or the S6 candidate) before S8. Report: `research/wet_antenna/outputs/reports/return_s6b.md` |
| Part 2 (reboot) | S8 | Freeze (hashes, committed), comparison models B0-B2 on TRAIN, then TEST once: pass rules P1-P5 and information items | done 2026-10-19 | P2-P5 pass (candidate 0.692 dB vs B0 0.862 / B1 0.760 / B2 0.727); P1 FAIL (4 of 8 cells: the antenna level at 0.5-2 mm/h is higher on TEST); two evaluation bugs fixed after the read (recorded). Report: `research/wet_antenna/outputs/reports/test_s8.md` |
| Part 2 (reboot) | S8b | Post-freeze evaluation of a pre-freeze model: the closed form F = A_band (R / 5)^b on TEST against the table F, by a rule written before the run | done 2026-10-19 | Same P1-P5 verdicts; paired RMSE difference within +-0.01 dB: the closed form is ADOPTED as the written form of F (b 0.54; valid 0-40 mm/h). Report: `research/wet_antenna/outputs/reports/test_s8b.md` |
| Part 2 (reboot) | S9 | Formulate: story, assumptions with checks, data funnel, core figures, FINAL.md | done 2026-10-19 | Equation file `waa_equation_final.json` (sha256 `b06ebc8b...25be9cd`), `waa_equation.py` with a unit test (1e-6 dB, achieved 9e-16), the overview checked (113 statements: 89 confirmed, 21 corrected, 2 unsupported, 1 added), one document `docs/WAA_research_overview.pdf` (after the learners merged the short and the full version), `FINAL.md`. **Part 1 closed.** Report: `research/wet_antenna/outputs/reports/formulate_s9.md` |
| Part 2 (reboot) | S8c | Post-freeze comparison with two literature WAA models (Schleiss 2013, Pastorek 2022 KR-alt; pycomlink 0.6.0), fitted on TRAIN only, our equation not refitted | done 2026-10-28 | TEST, version b: Schleiss no difference (+0.009 [-0.010, +0.039] dB), Pastorek ours better (+0.026 [+0.015, +0.038] dB); version a: ours better. Report: `research/wet_antenna/outputs/reports/literature_s8c.md` |
| Part 2 (reboot) | R1 | Does the WAA correction improve link rain estimates against the gauges? | next (added 2026-10-19) | Its own written plan first, and a decision about TEST (already read once) |

Order of work: S1 -> S2 -> S2b -> S3 -> S4 -> S5a -> S6 -> S7 -> S6b (the one return) -> S8 -> S9.

Order of Part 1 as it was run: T0 -> T1 -> T2 -> T3 -> T4 -> T5 -> T6 -> T7 -> T8 -> T9 -> T10 -> T11 -> T12 -> T13 -> T14 ->
T15 -> T16 -> T17 -> T18. Its plan was changed on 2026-10-11 after T9 (was 3b-2) and on 2026-10-12 after T15 (was 4b);
those sections are kept in Part 1 below.

## Part 1 (closed): drying-tail exploration T0-T18

Closed on 2026-10-13. The text below is the plan and the record of the tail exploration, kept as it was; only the step
numbers were replaced by the new labels (the old number in brackets). The steps it still calls planned (5, 6, 3c, 4d,
7, 8a-8e, 9) are CANCELLED: the Reboot above (S1-S5) is the plan. What each T step decided, with numbers:
`research/wet_antenna/METHOD.md`.
Do not confuse: "T1 noise", "T2 starts" and "T3 truths" in the text are the three stress tests INSIDE step T8 (was 3b-1c),
named before the renumbering; they are not the steps T1, T2, T3.

### The goal of the tail plan (as it was)

**Goal.** A dynamic WAA law, `dW/dt = g(R, W) - h(W)`.
- `h` = drying, learned from drying tails (R = 0).
- `g` = wetting, learned from very short links, where path rain is small.
- Compare with the literature: Minda & Nakamura exponential drying (`dW/dt = -W/tau`), Schleiss wetting (rise toward A_max).

### What the final result is

- One equation for the wet-antenna attenuation W (dB), per minute: `dW/dt = g(R, W) - h(W)`.
  h(W) = drying (no rain), learned from drying tails. g(R, W) = wetting, learned from very short links.
- With numbers: parameters with uncertainty, how much they vary between links, and what they depend on
  (frequency, antenna, weather).
- With plots of the drying and wetting curves, a comparison with the literature (Minda & Nakamura drying, Schleiss wetting),
  the prediction error on unseen events vs the literature, and the limits.
- Use: subtract the predicted WAA from a link's attenuation; what is left is rain attenuation (input to Stage 2).

### Method: how the equation is discovered from the data

The recipe is the "universal differential equation" approach of scientific machine learning (Rackauckas et al. 2020), as used for
population data with between-individual variability in Martensen et al. 2024 (DeepNLME + symbolic regression):
- Jump 1 (learn the shape): a flexible function inside the model where the law is unknown (spline for one input; small neural
  network for several inputs), fitted with our pooled likelihood (offset prior, AR(1) noise, per-tail speed = random effect).
- Jump 2 (name it): turn the learned function into a short formula. Candidates ordered from simple to complex.
- Jump 3 (refit and choose): put each candidate formula back into the full pooled model, refit on the tails, choose by BIC on
  held-out events. Only Jump 3 gives the final parameters and their uncertainty.

Every jump is tested on fake data with known truth before real data.

### Tools and where they are used

| Tool | What it does | Where |
|---|---|---|
| Free curve (spline) / small neural network | learns the shape from the data | Jump 1 in T10 (was 3b-3), T12 (was 3b-3c), T14 (was 4) (spline); 6 and 7 (several inputs) |
| Best-subset fit + BIC | tries every combination of a few building blocks | Jump 2 for the drying curve (main) |
| PySR (symbolic regression) | searches any formula (powers, exp, fractions); simple-to-complex list | Jump 2 when the shape is not in our blocks; 7 |
| Ensemble SINDy (E-SINDy) | SINDy repeated over resampled events; keeps terms that keep coming back | cross-check in T14 (was 4); main tool in 6 and 7 |
| Hybrid / switching idea (Mangan et al. 2019) | different equations in different W ranges | T12 (was 3b-3c): Jump 1 may draw a step at a place it finds itself |
| SINDy with inputs (SINDYc) / parametric SINDy (SINDyCP) | equations driven by an input or with condition-dependent coefficients | 7 (rain as input); 6 (weather changes the drying speed) |
| Menu fits F1-F6 | tests our own candidate laws | comparison and literature check in T14 (was 4) |

Note: pysindy and scikit-learn are not installed; pip install is allowed. PySR needs a Julia backend (installed by pip on first
use; may be slow). Say in each step which software was used.

**Standing rules**
- Uncertainty by resampling EVENTS; leave-one-event-out. 3 events give 246 of the 429 train tails.
- Sensitivity runs at every modelling step: v1 main set; without `event_rise_heavy`; without `level_shift` / `undershoot`; alternative set (A_start >= 0.6 dB).
- Tools are tested on simulated tails (known truth) before they are used on real tails.
- No tail is chosen by how well it fits a model (circular).
- Every script writes its report, tables and figures to the module's `outputs/`; all of it is committed.
- Each step must either produce a real-data result or answer a question that blocks one. (added 2026-10-12)
- Simulated data must contain every behaviour seen in the real data before a tool is trusted on real data. (added 2026-10-12)
- The tail definition is a result to be checked, not a given: report how conclusions change with the start rule. (added 2026-10-12)

### Plan change (2026-10-11, after T9 (was 3b-2))

SINDy on raw tails failed the fake-data test (T9 (was 3b-2)). Discovery is kept but split into two jumps: (1) the pooled fit with a
FREE smooth drying curve instead of a menu of laws draws the drying curve from all tails; (2) SINDy names the equation of
that clean curve. Old Step 5 (SINDy on raw real tails) is replaced by this; it becomes part of T10 (was 3b-3) (fake) and Step T14 (was 4) (real).
3c is moved to just before Step 7 (it is only needed for wetting). Step 7 keeps "learn g(R, W)"; the tool for it is decided
after 3c, by the same fake-data testing.

### Plan change (2026-10-12, after T15 (was 4b)): the signal before the law

**Why.** Steps 3 to T15 (was 4b) built and tested the fitting tools on simulated tails before the real tails were understood.
The real tails then showed behaviour the simulator never had: tails that wait, rise or stall, and noise 2.6 to 4.3
times the dry-weather model. The tail definition (start at t0_trace, A_start >= 1 dB, first 60 min) was frozen early
and shapes every result. So the order is turned around: first establish WHAT the signal in a tail is, then fit laws.

**Open questions, ranked, each with its test.**
- **Q1. Is a night-time rise water at all?** Candidates: baseline drift with the daily cycle (hardware / temperature);
  dew or condensation on the radome (this IS WAA, without rain); drizzle under the radar. Test: Phase A1.
- **Q2. Is the early fast drop after t0_trace the antenna drying, or path rain fading?** Test: Phase A2 (A vs L per minute).
- **Q3. When does drying start?** The radar is 5 min / 2 km and misses drizzle, so t0_trace may be too early or too
  late. Test: Phase A3 (several start rules).
- **Q4. Does drying (its start delay and its speed) depend on conditions** (night / sun, temperature, humidity, wind)?
  Phase C = Step 6. Expect weak power: about 24 events, one hourly station, RH missing for two events.
- **Q5. The form of h(W)** (exponential or not). Only after Q1-Q3.

**Phases and their decision points.**
- **Phase A (Step T16 (was 4c)).** A1 daily cycle / dew vs water; A2 antenna vs path; A3 start rules.
  - If A1 shows drift: correct the baseline for the daily cycle (baseline v3, new validation) and re-cut the tails.
  - If A1 shows dew: rises are WAA with a non-rain input (wetting side; check 8e).
  - If A2 shows path rain in the first N min: the drying law is fitted from minute N on.
  - A3 then proposes the rules of tail set v2 (chosen without any law fit); the learners decide and freeze v2.
- **Phase B (Step 4d).**
  - B0: the simulator gets what the real tails showed (start delays, rises, the daily cycle, noise 2.6-4.3 x), and
    the tools must pass on THAT before use.
  - B1: the drying law on tail set v2 (menu fits, free curve + E-SINDy), with a per-tail start delay as a model
    option instead of dropping waiting tails, and a realistic noise level.
- **Phase C (Step 6).** What the delay and the drying speed depend on (per event).
- **Phase D (Steps 3c, 7).** Wetting from very short links; consider modelling WHOLE events there (radar rain as
  input), which removes the start-of-drying problem.
- **Phase E (Steps 8, 9).** Verification on unseen events; final equation.

**The realistic form of the final drying statement.** "Once drying starts, the antenna part decays [roughly
exponentially] with tau = X +- Y min; drying starts later [under conditions] by D min; night rises are [dew / drift]."

#### Decisions of the learners after T16 (was 4c) (2026-10-12)

1. **Night / day rises are NOT baseline drift** (A1: 67 % of the rising tails exceed the 95th percentile of their own
   dry days). No baseline v3. Dew or drizzle stays open (humidity missing for E20150617_0705 and E20150701_0030).
2. **Tail set v2 = the v1 tails with the start S0 (t0_trace), plus two more times stored per tail**: S1 (strict dry:
   no rain within 5 km, + 10 min) and S2 (signal onset), in `research/wet_antenna/outputs/products/tails_starts.csv`.
   The set is not re-cut. The start of drying is treated per tail.
3. **New key finding.** Rain within 5 km stops about 25 min after t0_trace (median), and by then 81 % of the tails are
   below 1 dB. Most of what v1 calls the drying tail happens while it still rains nearby. So the fast early drop may
   be the antenna drying, light rain near or on the path ending, or both. **Q2 (antenna or rain) is now the central
   open question.** A2 (A against L) had no power on tails.
4. **Order changed.** T16 (was 4c) -> T17 (was 4c2) (quick check: how fast does A fall once the rain nearby has ended?) -> Phase D moved
   forward: 3c + 7 = WHOLE EVENTS on very short links with radar rain as input (the model sees the rain end by
   itself, so there is no guess of where drying starts) -> 4d (the drying law on tails, informed by 7) -> 6 -> 8a-8e
   -> 9.

### Corrections found when checking the repo (METHOD.md vs what I said in chat)

1. **Fast drop is not proven to be WAA.** Step T2 (was 2) found it does not grow with L, but the path-rain reference `wet_A_max`
   barely grows with L either (Spearman 0.13), so the L test has little power. Step T3 (was 2b) then showed that 381 of the 959 step T2 (was 2) main tails
   fall below 1 dB at `t0_trace`: much of the step T2 (was 2) fast drop happened while radar or gauge still showed rain. Only the
   sub-link agreement (rho 0.87) and the flat control curve support WAA; both are weaker evidence than I claimed.
   -> The very-short-link length check in Step 7 is therefore the real test of "WAA, not path rain", with enough contrast in L.
2. **Three events, not two, give most tails** (246 of 429).
3. **With spread, a plain exponential (m1) reproduces the slow tail as well as two layers (m4).** The median shape does not favour m4.

### Part 2: practise on fake data

**T5 (was 3a) / T5 (was 3a-fix) (done).** Simulator on the v1 skeleton. Noise = 10 s white + 1 min AR(1) wander + 0.3 dB quantization + per-tail offset.
Laws: m1 exponential, m3 `-c W^2`, m4 two layers, m5 `-aW - bW^2`, each with a speed spread s = 1.25. m2 was dropped (= m1 with spread).

**T6 (was 3b-1) (done).** Fit F1 exponential, F4 two layers, F3 nonlinear, F5 linear+quadratic to each fake tail (free W0 and free offset c);
BIC with an effective sample size for slow noise. Truths m1, m3, m4, m5, m6 (m6 calibrated: k_run 0.044, k_evap 0.0055, A_c 0.16).
Result (10 replicates per truth):
- Per tail: no skill. F1 is chosen for 60-80% of tails whatever the truth.
- Per dataset (summed BIC): m1 -> F1 in 10/10; m3, m4, m5, m6 -> F4 in 8-10/10. Own candidate wins only for m1 (100%) and m4 (80%); never for m3, m5.
  So the method separates "plain exponential" from "slower than exponential", not which slow law. F4 winning does not mean two layers.
- Recovery: exponential tau unbiased (1.01); m4 tau_s x4.5 too large; m3 n at its bound for 40% of tails; m5 a, b not recoverable. Free offset c trades off with a long small slow layer.
- Ablation: no single effect (noise, quantization, offset, spread) is to blame; only "all off" recovers the truth.
- AIC picks F4 for every truth; use BIC.
- Population tests (t_half vs A_start, dA/dt vs A) do not separate the laws. t_half: real inside m1 range 1/4 bins, m3 3/4, m4 3/4, m5 4/4, m6 3/4
  (weak hint that wetter dries faster; replicate ranges hold noise only, not event-sampling uncertainty).
- Real dA/dt above ~1.5 dB is thin and zig-zags to ~0.

**T7 (was 3b-1b) (done).** Raised the power of the fit (`pooled_fit.py`): (a) prior on the offset, c ~ N(0, offset_std^2), instead of a free c;
(b) AR(1) noise in the likelihood instead of the n_eff correction; (c) pooled fit per dataset: shared law parameters, per-tail
speed factor log-normal with a per-link part (marginalised); (d) F6 (switch fit) added; (e) replicates that resample events;
(f) first 60 min vs 6 h vs whole tail. Result: on 10 event-bootstrap replicates per truth the pre-set rule (own candidate >= 80 %
for m1, m4, m5, m6, key parameters within +-30 %) is met in all three windows; the first 60 min is the best window.

**T8 (was 3b-1c) (done, PARTIAL).** Tried to break the pooled fit on fake data (first 60 min, 5 replicates): T1 noise with two
slow time scales, T2 generic optimiser starts, T3 truths recalibrated to the real early decay (+ s = 1.5). m1, m4, m6
pass every test. m5 fails twice: with the unknown noise F3 wins 5/5 (by 1.8 BIC), and with the faster truth F5's
linear rate a is x0.65. The failing part is understood (m5 is almost a pure quadratic, which F3 copies).
Carried to Step T14 (was 4): several starts per candidate (F6 is start-sensitive in 20/20 datasets); F3 and F5 are one family;
F4 vs F6 by < ~10 BIC is not decisive.

**T9 (was 3b-2) (done; Q1 FAIL, Q2 FAIL, Q3 PARTIAL, Q4 FAIL).** Weak-form SINDy (own code, `sindy.py`) on the T8 (was 3b-1c) fake
datasets. One pooled law dA/dt = h(A) returns "A" for the quadratic (m3, m5), switch (m6) and two-layer (m4) truths;
exponential tails of mixed speed get a false A^2 in 3 of 5 replicates (none without spread); it agrees with the pooled
fit in 1 of 20 datasets. With the true offset and no spread it does find A^2, so the offset and the speed spread are
what break it. Consequence: Step 5 cannot be run as planned (see the plan change above). Next: T10 (was 3b-3). The original plan was:

**T9 (was 3b-2) (plan).** SINDy on fake tails. Library: 1, A, A^2, A^3, switch term max(A - A_c, 0) for several A_c.
Checks: does it return the true law for m1, m3, m5, m6? What does it return for the hidden-state m4? Does mixing tails with
different tau create a false nonlinear term?

**T10 (was 3b-3) (done, PARTIAL).** By the pre-set rule: m1 FAIL, m3 / m5 / m6 PARTIAL, m4 PASS (warning in 4 of 5). The failures sit
below about 0.35 dB, where the first knot (0.05 dB) is badly determined, and in a broken interval (available for 3 of 65
fits). Above the second knot, first 60 min, with spread (post hoc): m1 flat (slope -0.02, level x1.03, equation "W" 5/5),
m3 slope 1.08 with W^2 x0.94, m5 slope 0.96 with W^2 x1.06, m6 slope 0.27 (true 0.26): the mixing trap is gone there. The curve
cannot show m6's step at 0.16 dB nor m4's second state; the menu fits F6 / F4 are needed for those. Before Step T14 (was 4): a working
interval (event bootstrap or profile), knots / range that do not leave 0.05-0.35 dB to one value, a longer search; and check
the suspicion that the likelihood without quantization favours fast drying at the low end. The plan was:

**T10 (was 3b-3) (plan).** "Draw the curve, then name it" on fake tails. Jump 1: the pooled fit of T7 (was 3b-1b) / T8 (was 3b-1c) with a free
smooth drying curve, dW/dt = -f_i W k(W), k(W) = exp(spline in log W), instead of a menu law. Jump 2: SINDy on that clean
curve h(W) = W k(W) to name its equation. Question: with a known truth, is the curve right and is the equation right?

**T11 (was 3b-3b) (done, PARTIAL).** Fresh data, rule fixed before the run, judged range 0.35 dB to W_hi (1.9-2.3 dB), first 60 min.
m1, m3, m5 PASS on all three criteria: true k inside the replicate band on 100 % of the range, level within 3 %; slopes
of log k vs log W -0.01..0.05 (m1), 0.99..1.12 (m3), 0.93..1.03 (m5); one-term equation in 5 of 5 (W x0.99, W^2 x0.97,
W^2 x1.13); no false term for m1, so the mixing trap is gone. Same under T1 noise. m6 (switch below W_lo, not judged):
the curve is right (slope 0.20-0.32, true 0.23) but reads as "mildly rising", with a three-term equation. m4 FAIL: the
free curve is worse than F4 by > 10 BIC in only 2 of 5 (T10 (was 3b-3): 4 of 5), and its curve looks like m6's (slope 0.36-0.47).
So: the curve tells exponential from quadratic from "something else"; switch vs two layers stays with the menu fits, and
"free curve vs F4" is not a reliable check. The plan was:

**T11 (was 3b-3b) (plan).** Confirmation run. The "right above 0.35 dB" result of T10 (was 3b-3) was found after seeing the results, so it
is tested here on fresh fake data (seeds 4000 + replicate) with the rule written before the run: first 60 min only, 5 knots
from W_lo = 0.35 dB to W_hi, nothing judged or named below W_lo, band = min-max of the replicates (no curvature interval),
shape by the slope of log k vs log W in each replicate, jump 2 with switch levels 0.5 and 0.75 dB only.

**T12 (was 3b-3c) (done, PARTIAL).** Rules fixed before the run; data = the T11 (was 3b-3b) datasets. A, Jump 2 alone on the true curve + 2-5 %
noise: E-SINDy (backward elimination of sindy.py, bagged) exact in 20 of 20 for m1, m3, m5, m6; PySR 20 / 20 / 18 of 20
for m1 / m3 / m5 (cannot write the switch); best-subset + BIC always contains the true terms but adds one in 3-8 of 20;
pysindy's STLSQ fails at 5 % noise. B, step: kept in 5 of 5 for EVERY truth (BIC gain 15-170), also the exponential, at
0.2-0.4 dB: the tool cannot be used. C, Jump 3: the step law wins the held-out likelihood in 18 of 20 single-curve
datasets, by 16-65 over the true law, so the pooled likelihood itself prefers a false feature near 0.2-0.4 dB (suspected:
the 0.3 dB quantization is not in the likelihood; not tested). Menu laws alone with BIC: true family in 4-5 of 5.
D, band from an event bootstrap of one dataset: +-1 to 6 % wide, true k inside on 0-78 % (rule 80 %): about 3 times too
narrow; slopes within 0.01-0.11 of the truth. E: F4 with the F1-solution start changes nothing (m4 warning stays 2 of 5).
Software: pysindy and PySR (with Julia) are now installed. The plan was:

**T12 (was 3b-3c) (plan).** (a) Jump 2 tested alone on the true k(W) of m1, m3, m5, m6 with a little noise added: best-subset + BIC,
PySR and SINDy must name the true formula (m6 included: is the switch named when the curve has the step?). (b) Jump 1 with
an optional step at an unknown W_c: must find W_c near 0.16 dB for m6 and "no step" for m1, m3, m5. (c) Jump 3: each candidate
put back into the pooled model and refitted; winner by BIC on held-out event draws. Rule fixed before the run, as always.

**T13 (was 3b-4) (done 2026-10-11, verdict "not the cause").** Diagnosis only, rule fixed before the run. Truth m1 with the T11 (was 3b-3b)
skeleton and seeds, simulated with no / 0.3 dB / 0.1 dB RSL quantization, free curve + optional step as in T12 (was 3b-3c) part B.
The step is kept in 5 of 5 in all three settings (BIC gain 40-130 without quantization), and the fitted W_c stays at its
start value (0.12, 0.2 or 0.4 dB) in 15 of 15 fits. So the 0.3 dB quantization is not the cause of the false step, and a
quantization-aware likelihood is not expected to repair it. Real data: RSL in 0.3 dB steps on 719 of 728 sub-links, 10 s
samples, 1 min means; the simulator matches. Cause still open (suspects: the extra search budget of the step fit, the
offset prior / start floor, the narrow step shape). Next = Step T14 (was 4) with the reduced recipe of T12 (was 3b-3c), unchanged.

### Part 3: find the law on real data

- **T14 (was 4) (done 2026-10-12: no drying law identified).** TRAIN v1, 429 tails, first 60 min, W_hi = 2.83 dB; 36 min, nothing
  cut. The free curve falls with W (slope -0.87; bootstrap refits -1.43..-0.82; sets -0.05 to -1.60), which no fake truth
  did: 21 % of the tails stay at >= 90 % of their start level for 10 min (38 % in the rise-heavy events) before they drop.
  Without the rise-heavy events the slope is -0.05 (exponential-like, tau 15 min). Jump 2 gives no formula (5 of 5 terms
  kept). Menu as run: F1 (tau 25.2 min), but F1 / F4 may fit a rising tail with a negative start (3 tails = the whole
  2600 BIC gap to F3 at n = 1); like for like F6 is ahead by 21 BIC, not decisive (residual variance 4.1 times the
  dry-noise model). Stable in all sets: F3 at n = 1, F5 at b = 0 (nothing rising with W). No literature value in the
  repo to compare tau with. Open before Step 6 (learners): which tails define "drying"; repair the W0 bound and re-test
  on fake data with rising tails; a noise level for drying tails. The recipe was: Only what passed on fake data is used: Jump 1 without a step, read for
  W >= 0.35 dB; shape = slope of log k vs log W (about +-0.1); Jump 2 by E-SINDy (the elimination of `sindy.py`) with PySR as
  cross-check and best-subset as "true terms plus at most one extra"; NO Jump 3 by held-out likelihood; named laws compared
  by the menu fits with BIC; the event-bootstrap band is shown with the warning that it was about 3 times too narrow on fake
  data (level about +-20 %). Nothing is said about W < 0.35 dB or about a switch until the cause of the false step is found
  and B, C, D are repeated on the cached fake datasets (T13 (was 3b-4): the cause is NOT the 0.3 dB quantization). The earlier plan, kept for reference: The 429 real tails; every result is stated for W >= W_lo = 0.35 dB (from T11 (was 3b-3b); up to W_hi)
  and the first 60 min only:
  (a) Jump 1: free-curve pooled fit (5 knots from 0.35 dB to W_hi, band from an event bootstrap) -> real drying curve and
  its slope of log k vs log W, judged only for W >= W_lo;
  (b) Jump 2: best-subset + BIC (main), PySR and E-SINDy as cross-checks (library W, W^2, W^3, switch at 0.5 / 0.75 dB) -> equation;
  (c) Jump 3: each candidate refitted in the pooled model, chosen by BIC on held-out events;
  (d) compare with the menu fits F1, F3, F4, F5, F6 under the T8 (was 3b-1c) rules (several starts, F3 / F5 one family, F4 vs F6
  < ~10 BIC not decisive) and with the literature exponential;
  (e) the four sensitivity sets.
- **T15 (was 4b) (done 2026-10-12: neither explanation passes its rule; the stall belongs to the event).** Why 44 % of the TRAIN
  tails do not decay at once (64 wait 10 min or more and then decay, 36 rise, 89 are flat or slow). H1 (something
  still wets the antenna): rain references show LESS near these tails (odds ratio 0.29), but the path radar is dry by
  construction of v1 and humidity is missing for the two events with most rising tails. H2 (rate per dB falls with
  A): 0 of 4 rules; in decaying tails k does not fall with A (+0.28, all tails -0.80), the same link does not repeat
  (p = 0.35), and after its plateau a waiting tail follows the same curve as the others. Shared inside an event
  (coherence 1.85, interval 1.27..2.43), not related to path length. Consequences: h(W) needs no falling k; the tails
  that define h(W) are the decaying phase ("decays" 240 tails, or + plateau ends 304: free curve flat above 0.35 dB,
  F1 tau 10-11 min, biased short); the stalled minutes need a start delay, a weather-dependent rate or an input
  that is not rain: that is Step 6, which should start with night / day, humidity and wind. Tool fixes in
  `pooled_fit.py`, off by default: start bound for F1 / F4, fitted noise multiplier (2.6 to 4.3 on real tails; its
  fake check did not pass, 4 of 6). Open for the learners: which tails define h(W); delay, weather or input; make the
  fixes default; humidity from another station; is the fast fall of the decaying tails the end of the path rain.
- **T16 (was 4c) (done 2026-10-12, Phase A: drift ruled out for most rises; path rain and dew not decided).** What a drying tail
  contains, on the TRAIN tails of v1, no law fit. A1, the tail's rise against the same sub-link at the same clock
  time on certified dry days: "mostly water, not drift" by the rule (24 of 36 rising tails above their dry-day 95th
  percentile; dry days rise 0.30 dB by day and 0.36 dB at night at that percentile), weaker for the night-time rises
  (10 of 18) and partly built into the class; dew not decided (humidity missing for the two events with 22 of the 36
  rising tails). A2, A against L per minute from -30 to +60 min: power check not passed, nothing concluded (the rain
  before t0_trace is light: 0.03-0.06 dB/km expected). A3, three start rules: strict dry (no reference within 5 km)
  comes 25 min after t0_trace and leaves 129 slow tails; signal onset equals t0_trace for 70 % of the tails; k(A)
  falls under every rule. Proposal for tail set v2: keep t0_trace and store the two other starts per tail
  (`tails_starts.csv`); the learners decide and freeze. Decision points of the plan: no drift -> no baseline v3;
  "path rain in the first N min" not established -> no minute N; Q2 goes to Step 7.
- **T17 (was 4c2) (done 2026-10-12: SAME speed by the rule; two phases in time).** Quick check, data only, no law fit: tails
  started at S1 (no reference within 5 km shows rain) with lower thresholds (0.25 / 0.4 / 0.6 dB at S1), against a
  placebo start with the same delays but no "rain has ended" condition, and the model-free k at the same level before
  and after S1. Rule fixed first: SLOW (the fast early drop belongs to the end of the rain; tau about 10 min is not a
  drying time), SAME (the fast drop is the antenna drying) or unclear. Result at 0.4 dB: t_half 40 min from S1, 37
  from the placebo, 18 from S0; S1 - placebo +3 min (-2..+7), so SAME: starting late, not the end of the rain nearby,
  makes the S1 tails slow. The reading "SAME = antenna drying" is stronger than the test (the fast drop starts at
  the last rain trace on the path, and 52 % of the placebo starts are themselves after S1; at 0.6 dB the S1 tails are
  slower than the placebo). What Step 7 takes from it: two time scales (about 10-15 min, then t_half about 40 min with no
  rain within 5 km), and the question whether the fast one scales with L.
- **4d (planned, Phase B; moved after Step 7 on 2026-10-12).** CANCELLED 2026-10-13: see the Reboot (S1; drying is one rung of the ladder). B0 simulator with start delays, rises, the daily cycle and the real noise level; the
  tools must pass on it. B1 the drying law on tail set v2 with a per-tail start delay as a model option.
  Added 2026-10-13: look at tails up to 3-6 h, not only 60 min; the model must allow two phases (fast about
  10-15 min, slow about 40-60 min half-life, from T17 (was 4c2)).
- **5.** CANCELLED 2026-10-13 (see the Reboot, S2). SINDy on raw real tails -> drying law h(W). **Replaced by T10 (was 3b-3) + T14 (was 4)** (plan change 2026-10-11): as tested in T9 (was 3b-2),
  SINDy on raw tails returns "A" whatever the true law.
- **6.** CANCELLED 2026-10-13 (see the Reboot, S4). What drying depends on. Jump 1 learns drying as a function of W AND weather / link features (temperature, humidity,
  wind, time of day, frequency, antenna size); then E-SINDy / parametric SINDy picks which features matter, with stability over
  events and leave-one-event-out. First a fake-data test with a known weather effect and one with none (must not invent an
  effect). Caution: about 24 distinct events, 3 hold most tails; weather is hourly from one station; RH is missing for part of June.
  Added 2026-10-13: frequency band (23 / 28 / 38 GHz) explicitly (Tiede et al. 2025: the effect depends strongly on frequency).
- **3c (planned; moved here, only needed for Step 7).** CANCELLED 2026-10-13 (see the Reboot: the fake check of S1). Extend the simulator to full events: real radar rain series on short links, wetting law (Schleiss form), drying law,
  a small path term that grows with L, noise. Test the wetting tools of Step 7 (below) with rain as input.
- **7.** CANCELLED 2026-10-13 (see the Reboot, S1 and S2). Full WAA law from very short links (< ~0.6 km), whole events.
  - Path rain must be small without assuming the power law: among short links, attenuation during rain must not grow with L; both directions agree; restrict to light/moderate rain if heavy rain shows an L effect.
  - Wetting as dW/dt = NN(R, W) - h(W), h fixed from Step T14 (was 4) (universal differential equation), radar rain along the path as
    input; then PySR / E-SINDy on the learned NN -> g(R, W); Jump 3 refit. Fake-data test (3c) first; the tool is kept only
    if it passes there.
  - Compare with Schleiss (dynamic), Valtr (WAA = 2 k' R^alpha') and Kharadly-Ross (C (1 - exp(-d A))). Also test whether
    wetting needs time dependence at all (Pastorek et al. 2022 found the rain-rate models best and the dynamic Schleiss
    model adding little on their links).
  - Support: WAA at tail start vs event rain amount and duration (saturation, A_max).
  - Added 2026-10-13: look at tails up to 3-6 h, not only 60 min; the model must allow two phases (fast about
    10-15 min, slow about 40-60 min half-life, from T17 (was 4c2)).
  - **T18 (was 7a) (done 2026-10-13: undecided by the rule; the fast drop does not grow with L).** Very short links: is there
    enough data for Step 7, and do they show the fast first phase? Inventory per length class, the two phases by
    length (SHORT against LONG, matched on A_start and event, rule fixed before the run), a first look at whole events
    on links below 1 km. No law fit. Result: the fast drop is the same on every length (0.27-0.30 dB; matched ratio
    0.94, 0.84..1.15; path rain would give about 0.2), but the pre-set power condition fails, because the radar
    expects almost no path rain before a tail (0.02 / 0.08 dB). The slow phase is what differs, weakly (long links
    keep more). What Step 7 takes from it: links below 0.6 km are enough by the rule (below 1 km as the larger
    set); multi-shower events; most of the fall of A happens while the rain fades; no A_max visible, and the path
    term is not small in heavy rain even below 1 km; wetting before rain is rare; band (28 / 38 GHz) as an input.
    To decide (learners): accept the reading "fast phase = not path rain" although the rule says undecided; the
    length cut.

### Part 4: verification (after discovery)

CANCELLED 2026-10-13: 8a-8e are replaced by S3 of the Reboot.

Final choice by BIC / prediction error on unseen events (VAL, leave-one-event-out), as in Martensen et al. 2024.

- **8a.** Drying on unseen data (VAL events, leave-one-event-out); our law vs the literature exponential.
- **8b.** Full law on unseen events (short links); vs Schleiss.
- **8c.** Transfer to longer links: A minus predicted WAA should grow with L with no constant offset (the f = L + L0 test). Handoff to Stage 2.
- **8d.** Consistency: SINDy tau vs directly fitted tau; sub-link pairs; stability under event resampling.
- **8e.** Dew events should dry by the same h(W).

### Part 5: final (9)

CANCELLED 2026-10-13: replaced by S5 of the Reboot.

Full WAA equation (wetting + drying), parameters and their spread, comparison with the literature, limits, input to Stage 2.

### Naming: m and F

m = the law that MADE the fake data (the truth); F = the law the fitting tool tries as an answer. Same number, same law.
1 exponential (Minda & Nakamura) | 3 `-c W^n` | 4 two layers | 5 linear + quadratic | 6 state switch.
There is no 2 (m2 = m1 with spread, dropped in T5 (was 3a-fix); not renumbered) and F2 never existed.
m3 was not a truth in T8 (was 3b-1c); F3 still competes. F3 and F5 are read as one family ("nonlinear in W").

### Planned routine: /research-report

LaTeX technical report on the research part in the style of the pipeline report (pipeline_log.pdf): scope and sources,
overview, data and conventions, one section per step (Purpose / Method / Parameters / Results / Test and verification /
Issues found and corrections), claims-and-status table, limitations, chronology, references, glossary, an "In plain words"
box per chapter. Numbers injected from METHOD.md and outputs/tables by a script, never typed. Prompt prepared; run when the
machine is free. (The shorter research log is kept up to date by `/updatedocs`.)

### Options kept open (not forced)

- m6 / F6 state-switching law: fast runoff above a threshold A_c, slow evaporation below. Tested in T6 (was 3b-1) (as fake truth), T9 (was 3b-2), T14 (was 4), 5.
- Switching law (m6 / hybrid): tested in T12 (was 3b-3c) with a Jump 1 that may draw a step. Result: the step is drawn for every truth,
  so the tool says nothing yet. T13 (was 3b-4): the 0.3 dB quantization is not the cause; the cause is open.

### Literature used for the method

- Rackauckas et al. 2020, Universal Differential Equations for Scientific Machine Learning, arXiv:2001.04385
- Martensen et al. 2024, Data-driven discovery of feedback mechanisms in AML (DeepNLME + symbolic regression), bioRxiv 10.1101/2024.06.17.599366
- Fasel et al. 2022, Ensemble-SINDy, arXiv:2111.10992
- Mangan et al. 2019, Model selection for hybrid dynamical systems via sparse regression, arXiv:1808.03251
- Cranmer 2023, Interpretable Machine Learning for Science with PySR and SymbolicRegression.jl, arXiv:2305.01582 (proposal ref. [7])
- Pastorek et al. 2022, Precipitation estimates from CMLs: practical approaches to wet-antenna correction (engrXiv; proposal ref. [9] is the IEEE TGRS version)
- Moroder / Chwala et al. 2023, ATTRRA2 field experiment: dynamics of WAA (drying after rain, dew before rain), IEEE / EGU23
  (the proposal cites the ATTRRA2 paper as ref. [2], Tiede, Chwala & Siart 2023: check the author list before citing)
- Tiede, Chwala, Siart & Eibert 2025, Frequency-dependent variations of the antenna reflection coefficient due to
  different wetness conditions on the antenna radome, Adv. Radio Sci. 23, 13-19, doi:10.5194/ars-23-13-2025 (lab,
  Ericsson antennas 23 / 28 / 38 GHz). Recorded 2026-10-13 as given by the learners. Checked on the publisher's page
  only (citation right; drops, water amount, several hours and the spread between repeats are in the text; the page
  says "K and Ka band" and that temperature, humidity and air movement could not be held constant, so the
  three frequencies and the last finding below are not confirmed from it). Its relevant findings:
  - water sits as drops / runlets, not a film;
  - the total water amount matters more than the drop positions (supports one state W);
  - drying takes typically several hours;
  - the same antenna dries in very different times between repeats (drop formation, runoff right after wetting),
    which supports a real speed spread and two phases (runoff, then evaporation);
  - the effect depends strongly on frequency;
  - temperature, humidity and wind change the drying.

  We cannot use its method: the reflection coefficient needs a vector network analyser; OpenMRG has only TSL / RSL.

Note: machine learning in CML work is mostly rain detection; discovering the WAA equation from data appears to be a gap
(not checked exhaustively).

### Open questions

Q1-Q5 are in the section "Plan change (2026-10-12, after T15 (was 4b)): the signal before the law".

- Night-time rises without any rain trace (E20150701_0030), and what keeps A flat or rising after the last rain trace
  in 21 % of the tails (steps T4 (was 2c), T14 (was 4), T15 (was 4b)): now **Q1** (drift, dew or drizzle; Phase A1) and **Q3** (when drying starts;
  Phase A3). Still needed for Q1: humidity for E20150617_0705 and E20150701_0030 from another station.
- Is the stall "no drying" (night, cool, calm) rather than wetting (step T15 (was 4b))? Now **Q4** (Phase C = Step 6).
- Is the fast fall of the "decays" tails the antenna drying or the last path rain fading (step T15 (was 4b))? Now **Q2**
  (Phase A2). It decides whether tau = 10 min is a drying time.
- Test events and test links have no baseline yet (needed for the final held-out check).
- Step T14 (was 4): is the dip of k(W) near the third knot (0.9 dB main set, 0.7 dB alternative set) made by the knots or by
  selecting tails by A_start?
- Step T14 (was 4): drying time constants of the proposal's refs [1] and [9] are not in the repo; read them before Step 8a.
