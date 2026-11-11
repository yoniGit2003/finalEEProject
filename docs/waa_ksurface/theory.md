# Theory note of the k-surface branch (step A2, 2026-11-11; revised the same day by the Detachment correction)

Revision note: the first version took three things from Part 1: a per-link factor a_link, the 1 min lag, and a section that
used the Part 1 store and its numbers as the reference case. They are removed (Deviation 2 in `decisions.md`). The mathematics
of sections 2–5 is unchanged.

This note fixes, before any fit, what the model of this branch assumes, when its form is an exact rewrite, and what it
cannot represent. Notation as in the proposal: A(t) = TSL − RSL − baseline, c = the offset of the link-event, R = the rain
input (radar path mean / g_e, mm/h), L = the link length (km), W = the wet-antenna attenuation (dB).

## 1. The measurement equation

    A(t) − c = W(t) + L·P(R(t − d))

P(R) is the path loss per km (dB/km). d is a possible lag between the rain input and the link; it is not taken from
anywhere: it is measured from this branch's own data before the model is built, and it is 0 until then. W does not depend
on L; the path term grows with L. This is what the length test uses: at a fixed steady rain, the line of (A − c) against L
has the intercept W and the slope P.

## 2. The one-state assumption

W is ONE number per link, and its rate of change depends only on the rain now and on W now:

    dW/dt = G(R(t), W(t))                                                    (1)

This is the whole assumption about the dynamics. It says the antenna has no other memory: two moments with the same R and
the same W have the same dW/dt, whatever rain fell before. G itself is not given a form.

## 3. The steady state and the sign conditions

For a constant rain R, W settles where it stops changing. The steady WAA F(R) is defined by

    G(R, F(R)) = 0                                                           (2)

For F(R) to be the one value that W approaches from both sides, G must have the signs

    G(R, W) > 0  for  W < F(R)      (below the target: W rises, wetting)
    G(R, W) < 0  for  W > F(R)      (above the target: W falls, drying)      (3)

Condition (3) says: for each R there is exactly one zero of G in W, and it is stable. Two further conditions keep W
physical: F(0) = 0 (no rain, no steady water: a dry antenna stays dry), and F(R) ≥ 0. With (3) and F ≥ 0, a W that starts at
0 never becomes negative.

## 4. The factorization k·(F − W)

Define, wherever W ≠ F(R),

    k(R, W) = G(R, W) / (F(R) − W)                                           (4)

By (3), numerator and denominator have the same sign, so **k > 0**. Then (1) reads

    dW/dt = k(R, W)·(F(R) − W),   k > 0                                      (5)

This is a rewrite of (1), not an added assumption:

- Every G with one stable zero (3) gives exactly one pair (F, k): F is the zero, k is the ratio (4).
- Every pair (F, k > 0) gives a G = k·(F − W) that satisfies (3).

So (5) with a free k is as general as (1) + (3). What (5) buys: F can be measured separately (the length test, step A3), and
what is left to learn is one positive function k with a direct meaning: 1 / k is the local time constant, the time scale at
which W closes the gap to its target.

**At the target (W = F) the ratio (4) is 0 / 0.** What k is there depends on how smooth G is:

- **G differentiable in W at F, with ∂G/∂W < 0 there (the regular case).** Then k extends continuously to the target with
  k(R, F) = −∂G/∂W(R, F) > 0. Near the target W approaches it exponentially with the time constant 1 / k(R, F). If G is
  continuously differentiable in W, k is continuous everywhere:
  k(R, W) = −∫₀¹ ∂G/∂W(R, F + s·(W − F)) ds.
- **The degenerate case: ∂G/∂W = 0 at F.** Then k → 0 at the target. Example: G = −(W − F)³ gives k = (W − F)², and W
  approaches the target like 1 / √t, not exponentially: it never "arrives" at an exponential rate. (5) still holds away from
  the target, but k > 0 fails at the target itself. Our k is the exponential of a spline, which cannot be exactly 0. It can
  only become very small. So a fitted k that drops steeply towards the target is the fingerprint of this case, and it will
  be read as such, not as a new law.
- **G with a kink at F (only one-sided derivatives).** Then k has two different limits at the target: −∂G/∂W from below
  and −∂G/∂W from above. This is the wetting / drying jump of the next section.

## 5. Wetting and drying are fitted separately; k may jump at the target

Below the target the antenna gains water (rain arrives faster than it leaves). Above the target it loses water (runoff,
evaporation). These are different physical processes, so there is no reason for G to have the same slope on both sides of
its zero. We therefore learn two functions:

    k = k_wet(R, W)   if W < target      (wetting)
    k = k_dry(R, W)   if W > target      (drying)

A jump of k at the target does no harm to equation (5): at the target the factor (F − W) is zero, so dW/dt is continuous
there (it is 0 from both sides), and the solution is unique as long as k is bounded. A store that fills
faster than it empties is exactly such a case.

After the rain R = 0 and F = 0, so (5) becomes dW/dt = −k_dry(0, W)·W: a drying tail measures k_dry(0, W) directly. A
constant k gives an exponential tail; k ∝ Wⁿ⁻¹ gives the power-law drying dW/dt = −c·Wⁿ.

## 6. Differences between links

The model has no per-link parameter: every link of a band has the same F(R), P(R) and k.

    dW/dt = k(R, W)·(F(R) − W)                                               (6)

Links may differ in reality (radome, antenna size, age). Whether they differ enough to matter is a question for this branch's
own data: step E4 measures the spread of the per-link residuals of the length test across storms, with a criterion written
before the measurement. Only if that criterion is met does a per-link factor enter, and then its form (a factor on F, or
something else) is chosen from what the residuals show.

## 7. The numerical step

The model runs once per minute with R held constant over the minute, and k frozen at its value at the start of the minute:

    v = F(R),   W ← v + (W − v)·exp(−k·Δt),   Δt = 1 min                     (7)

For a constant k this is the exact solution of (6) over the minute. For a k that depends on W it is a first-order step
(k changes little within one minute if 1 / k is several minutes). The step never crosses the target, because exp(−k·Δt) lies
between 0 and 1. So W cannot switch between wetting and drying inside a step, and the step is stable for every k > 0. W = 0
before the event.

## 8. What one state can and cannot represent

**A constant k is the simplest member.** With k_wet and k_dry constant, (6) is a first-order store with one time constant
for wetting (1 / k_wet) and one for drying (1 / k_dry), and the step (7) is its exact solution over a minute. This is the
rung K0 of the ladder. Its unit test is the closed-form answer to a step in R: W(t) = F·(1 − exp(−k_wet·t)).

**A model with two stores is NOT a special case.** Suppose the antenna held water in two places with different time
constants, W = W_1 + W_2, each relaxing to its share of the target:

    dW/dt = (φ·F − W_1) / τ_1 + ((1 − φ)·F − W_2) / τ_2

The right-hand side depends on how W is split between the two stores, not on W alone. Two moments with the same R and the
same W can have different dW/dt: after a long rain most of W sits in the slow store and decays slowly; after a short shower
the same W sits in the fast store and decays fast. No function G(R, W) can do this, so no k(R, W) can.

What a one-state fit returns if the truth has two stores (this is why the fake truth (d) of A6 exists):

- its k_dry looks as if it **rises with W** (large W right after rain is mostly fast water; the small W left later is slow
  water), although neither store has a W-dependent rate;
- its residual after rain depends on the **preceding rain amount** at the same W. This is the signature that rule B3 tests
  before a second state is added.

The steady state is the same in both pictures: for a constant R, W_1 + W_2 → F. So F is comparable between them; only the
dynamics differ.

## 9. What "steady" means for the measured F

The length test measures the intercept of (A − c) against L in "steady rain". This
equals F(R) only if W has reached its target within that time: a fast process has, a process with a time constant of hours
has not. How long "steady" must last is therefore a definition this branch has to derive and test (step E3), and the
intercept against the rain duration is its check on real data. In the one-state
view, the F we use is "the WAA level reached in steady rain of the duration the data offers".

## 10. What is assumed and what is not

**Assumed**

1. The measurement equation: A − c = W + L·P(R(t − d)); W does not depend on L, the path term is linear in L.
2. One state: dW/dt depends on R now and W now only (1).
3. One stable steady state F(R) for each R (3), with F(0) = 0.
4. Wetting and drying may have different k; which one applies is decided by the sign of (target − W).
5. All links of a band share F, P and k (no per-link parameter), unless the pre-registered test of step E4 shows that the
   data demand one.
6. F(R) and P(R) are the same for all links of a band, rise with R (monotone interpolation between the measured classes),
   and are zero at R = 0.
7. R is constant within a minute, and k is frozen within a minute (7).
8. The rain input (radar path mean / g_e, or a 1 min gauge) is the rain the antenna sees, up to a lag d and up to the rain
   error measured in step E1.

**Not assumed**

1. No power law, for F, for P or for the drying (F and P are measured per rain class; k is free).
2. No form of the dynamics beyond one state: (5) is a rewrite of (1) + (3). k may depend on R, on W, on both, or on neither;
   the data choose (the ladder K0–K4).
3. No equality of the wetting and drying rates, and no continuity of k at the target.
4. No exponential approach to the target (it follows only if k is constant near it).
5. No second (slow) store. If the data need one, the one-state model is wrong, and rule B3 is the test for it.
6. No dependence of W on the link length.
7. Nothing about temperature, wind, humidity or dew: they are outside this model.
8. Nothing from Part 1: no per-link factor, no time constant, no store, no closed form, no lag, no rain class, no
   steadiness rule (the Detachment rule in `plan.md`).
