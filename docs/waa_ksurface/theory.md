# Theory note of the k-surface branch (step A2, 2026-11-11)

This note fixes, before any fit, what the model of this branch assumes, when its form is an exact rewrite, and what it
cannot represent. Notation as in the proposal: A(t) = TSL − RSL − baseline, c = the offset of the link-event, R = the rain
input (radar path mean / g_e, mm/h), L = the link length (km), W = the wet-antenna attenuation (dB).

## 1. The measurement equation

    A(t) − c = W(t) + L·P(R(t − d))

P(R) is the path loss per km (dB/km) and d = 1 min is the Part 1 lag between the rain input and the link. W does not depend
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
there (it is 0 from both sides), and the solution is unique as long as k is bounded. The Part 1 fast store is exactly such a
case (2.6 min while rising, 10.3 min while falling).

After the rain R = 0 and F = 0, so (5) becomes dW/dt = −k_dry(0, W)·W: a drying tail measures k_dry(0, W) directly. A
constant k gives an exponential tail; k ∝ Wⁿ⁻¹ gives the power-law drying dW/dt = −c·Wⁿ.

## 6. Where a_link enters

Links differ in their antennas (radome, size, age). As in Part 1, one positive factor per link scales the target:

    dW/dt = k(R, W)·(a_link·F(R) − W)                                        (6)

- F(R) is the steady WAA of the typical link (the length-test intercept, which pools the links); a_link·F(R) is the steady
  WAA of one link. The median of a_link is near 1; it has a weak prior on ln a_link (as in Part 1).
- a_link does **not** multiply P: the path loss per km does not depend on the antenna.
- k is evaluated at the link's own W in dB, not at W / a_link. This is a choice (rule A5), and it matters only when k depends
  on W: with K0 and K2 the two readings are identical. The alternative "the whole state scales with the link", which gives
  k(R, W / a_link), is not assumed and not tested in Block A.
- Wetting / drying is decided against the link's own target: wetting if W < a_link·F(R).

## 7. The numerical step

The model runs once per minute with R held constant over the minute, and k frozen at its value at the start of the minute:

    v = a_link·F(R),   W ← v + (W − v)·exp(−k·Δt),   Δt = 1 min              (7)

For a constant k this is the exact solution of (6) over the minute. For a k that depends on W it is a first-order step
(k changes little within one minute if 1 / k is several minutes). The step never crosses the target, because exp(−k·Δt) lies
between 0 and 1. So W cannot switch between wetting and drying inside a step, and the step is stable for every k > 0. W = 0
before the event.

## 8. Constant k and Part 1

**A constant k is the Part 1 FAST store.** Part 1 has

    dW_f/dt = (φ·u − W_f) / τ_f,   τ_f = 2.6 min while W_f rises, 10.3 min while it falls,   u = a_link·A_band·(R / 5)^b

"Rises" means W_f < φ·u, which is "below the target". So this is (6) with F(R) = φ·A_band·(R / 5)^b, k_wet = 1 / 2.6 per
min, k_dry = 1 / 10.3 per min, and the Part 1 code uses the same one-minute step (7). The ladder rung K0 with these numbers
and this F must reproduce the Part 1 fast store exactly (unit test of step A5).

**The full Part 1 two-store model is NOT a special case.** It has two states, W = W_f + W_s, with

    dW/dt = (φ·u − W_f) / τ_f + ((1 − φ)·u − W_s) / 240 min

The right-hand side depends on how W is split between the two stores, not on W alone. Two moments with the same R and the
same W can have different dW/dt: after a long rain most of W sits in the slow store and decays slowly; after a short shower
the same W sits in the fast store and decays fast. No function G(R, W) can do this, so no k(R, W) can.

What a one-state fit returns if the truth has two stores (this is why fake truth (d) of A6 exists):

- its k_dry looks as if it **rises with W** (large W right after rain is mostly fast water; the small W left later is slow
  water), although neither store has a W-dependent rate;
- its residual after rain depends on the **preceding rain amount** at the same W. This is the signature that rule B3 tests
  before a second state is added.

The steady state is the same in both models: for a constant R, W_f + W_s → u. So F is comparable between them; only the
dynamics differ.

## 9. What "steady" means for the measured F

The length test measures the intercept of (A − c) against L in "steady rain" (wet for 20 min, stable for 15 min). This
equals F(R) only if W has reached its target within that time. With time constants of a few minutes it has; a component
with a time constant of hours has not. Step A3b checks this with the intercept against the rain duration. In the one-state
view, the F we use is "the WAA level reached in steady rain of the duration the data offers".

## 10. What is assumed and what is not

**Assumed**

1. The measurement equation: A − c = W + L·P(R(t − 1 min)); W does not depend on L, the path term is linear in L.
2. One state: dW/dt depends on R now and W now only (1).
3. One stable steady state F(R) for each R (3), with F(0) = 0.
4. Wetting and drying may have different k; which one applies is decided by the sign of (target − W).
5. A link differs from the typical link only by a_link, which scales the target and nothing else; k is a function of the
   link's own W.
6. F(R) and P(R) are the same for all links of a band, rise with R (monotone interpolation between the measured classes),
   and are zero at R = 0.
7. R is constant within a minute, and k is frozen within a minute (7).
8. The rain input (radar path mean / g_e, or a 1 min gauge) is the rain the antenna sees, up to the lag of 1 min.

**Not assumed**

1. No power law, for F, for P or for the drying (F and P are measured per rain class; k is free).
2. No form of the dynamics beyond one state: (5) is a rewrite of (1) + (3). k may depend on R, on W, on both, or on neither;
   the data choose (the ladder K0–K4).
3. No equality of the wetting and drying rates, and no continuity of k at the target.
4. No exponential approach to the target (it follows only if k is constant near it).
5. No second (slow) store. If the data need one, the one-state model is wrong, and rule B3 is the test for it.
6. No dependence of W on the link length.
7. Nothing about temperature, wind, humidity or dew: they are outside this model.
