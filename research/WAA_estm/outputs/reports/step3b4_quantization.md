# WAA_estm step 3b-4: is the 0.3 dB quantization the cause of the false step? (rule fixed before the run)

Written by `research/WAA_estm/step3b4_quantization.py` (settings `waa.quant3b4`). Fake data only; diagnosis only: no new
likelihood, no fix. No real tail is fitted, the frozen tail set v1 is not changed.

**Verdict: not the cause.**

What was run: truth m1 (exponential; no step in the truth), the event-bootstrap skeleton, seeds
(4000 + replicate), noise and spread (s = 1.25) of step 3b-3b, 5 replicates, first
60 min. Three simulator settings that differ only in how the 10 s RSL values are rounded: not at all,
to 0.3 dB (the data of 3b-3b / 3b-3c), to 0.1 dB. The simulator draws the same random numbers in all three (re-simulated
0.3 dB data identical to the 3b-3b cache in 5 of 5 replicates). Fit: the free curve of 3b-3b
(same knot rule, 4 starts, penalty by cross-validation), then the optional step of 3b-3c part B (W_c starts
[0.12, 0.2, 0.4] dB, 250 + 400 evaluations, kept if BIC(no step) - BIC(step) >
6). The 0.3 dB rows are the cached fits of 3b-3b and 3b-3c, not refitted. The likelihood is unchanged
in all three settings.
Timing: replicate 0 took 4.8 min wall on 3 processes, projected total
18 min (budget 30 min); run time 17 min
(cached in `.cache/fits3b4/`; this run: 7 s). Software: numpy / scipy / numba only.

## 1. Setting x replicate (`research/WAA_estm/outputs/tables/step3b4_step.csv`)

`bic_gain` = BIC(no step) - BIC(step); `step_kept` = gain > 6; `wc_db` = fitted step position;
`height_ln` = step height in ln units; `no_step_k_ratio_below` / `_above` = fitted / true k of the curve WITHOUT the
step, median over [0.1, 0.35] dB and over [0.35 dB, W_hi] (information only);
`no_step_slope_above` = slope of log k vs log W above 0.35 dB (true: 0).

| setting | rep | step_kept | wc_db | bic_gain | height_ln | best_wc_start | lambda | w_hi_db | no_step_k_ratio_below | no_step_k_ratio_above | no_step_slope_above |
|---|---|---|---|---|---|---|---|---|---|---|---|
| none | 0 | True | 0.4 | 103.9 | 0.4 | 0.4 | 0.01 | 3.51 | 1.172 | 1.113 | -0.002 |
| none | 1 | True | 0.423 | 129.8 | 0.41 | 0.4 | 1.0 | 2.44 | 0.955 | 0.958 | 0.064 |
| none | 2 | True | 0.224 | 86.1 | 0.67 | 0.2 | 1.0 | 2.41 | 0.92 | 0.951 | 0.056 |
| none | 3 | True | 0.191 | 75.2 | 0.75 | 0.2 | 1.0 | 2.43 | 1.056 | 1.024 | -0.007 |
| none | 4 | True | 0.122 | 40.0 | 1.09 | 0.12 | 1.0 | 2.24 | 1.059 | 1.061 | -0.027 |
| 0.3 dB | 0 | True | 0.39 | 108.5 | 0.46 | 0.4 | 1.0 | 3.52 | 1.083 | 0.965 | 0.014 |
| 0.3 dB | 1 | True | 0.212 | 142.5 | 0.74 | 0.2 | 0.1 | 2.47 | 0.963 | 0.971 | 0.065 |
| 0.3 dB | 2 | True | 0.406 | 136.1 | 0.44 | 0.4 | 0.1 | 2.41 | 1.0 | 0.964 | 0.044 |
| 0.3 dB | 3 | True | 0.382 | 81.3 | 0.47 | 0.4 | 1.0 | 2.45 | 1.068 | 0.998 | 0.017 |
| 0.3 dB | 4 | True | 0.209 | 99.4 | 1.09 | 0.2 | 0.01 | 2.26 | 1.053 | 1.074 | 0.006 |
| 0.1 dB | 0 | True | 0.124 | 93.2 | 0.94 | 0.12 | 0.01 | 3.51 | 1.126 | 1.069 | 0.027 |
| 0.1 dB | 1 | True | 0.408 | 136.4 | 0.39 | 0.4 | 1.0 | 2.44 | 0.978 | 0.965 | 0.066 |
| 0.1 dB | 2 | True | 0.417 | 83.0 | 0.41 | 0.4 | 1.0 | 2.41 | 0.886 | 0.916 | 0.057 |
| 0.1 dB | 3 | True | 0.205 | 61.5 | 0.58 | 0.2 | 1.0 | 2.43 | 1.073 | 1.046 | 0.007 |
| 0.1 dB | 4 | True | 0.435 | 49.4 | 0.32 | 0.4 | 0.01 | 2.24 | 1.055 | 1.083 | -0.026 |

Per setting (`research/WAA_estm/outputs/tables/step3b4_decision.csv`):

| setting | step_kept | median_wc_kept_db | wc_kept_db | median_bic_gain | bic_gain_range | no_step_k_ratio_below | no_step_k_ratio_above |
|---|---|---|---|---|---|---|---|
| none | 5 of 5 | 0.224 | 0.40 0.42 0.22 0.19 0.12 | 86.1 | 40.0 to 129.8 | 1.056 | 1.024 |
| 0.3 dB | 5 of 5 | 0.382 | 0.39 0.21 0.41 0.38 0.21 | 108.5 | 81.3 to 142.5 | 1.053 | 0.971 |
| 0.1 dB | 5 of 5 | 0.408 | 0.12 0.41 0.42 0.20 0.43 | 83.0 | 49.4 to 136.4 | 1.055 | 1.046 |

![m1 fits per quantization setting](../figures/step3b4_m1_quantization.png)

## 2. Verdict

Rule (fixed before the run): "quantization is the cause" if without quantization the step is kept in <=
1 of 5 AND at 0.1 dB the median W_c is < 0.6 x the median W_c at 0.3 dB or the
step is kept in <= 1 of 5; "not the cause" if without quantization the step is still kept in >=
4 of 5; otherwise "partly / unclear".

**Not the cause: without quantization the step is kept in 5 of 5; at 0.1 dB in 5 of 5 with median W_c 0.408 dB = 1.07 x the 0.382 dB at 0.3 dB.**

## 3. What is quantized in the real data

Read from `pipeline/docs/methods.md`, `pipeline/docs/openmrg_summary.md` and the release table `sublinks`
(column `rsl_min_step_db`); nothing was fitted. **RSL** is reported in steps of 0.3 dB: the smallest RSL step is 0.3 dB
on 719 of the 728 sub-links (0.1 dB on one, no data on 8). **TSL** is reported in steps of 1 dB, but it is constant on
most sub-links (ATPC is off; `tsl_behaviour` = constant for 523 of 728), and a frozen tail ends just before its first
minute with a TSL step, so inside a tail only the RSL rounding acts. Both are sampled every **10 s**; the pipeline
stores them losslessly (L1) and the tails use the **1 min mean of the six 10 s values** of tl = tsl - rsl (L3), minus
the baseline. The simulator does the same: it rounds the 10 s value to a 0.3 dB grid and averages six values into one
minute. Its grid has a random position per tail; on real data the position is fixed per sub-link, but the baseline
is subtracted afterwards, so relative to W = 0 the position is also arbitrary. **So the fake simulator matches the real
data in what is quantized (RSL, 0.3 dB), where (10 s samples) and how it reaches the 1 min series.** Not in the
simulator: the 1 dB TSL steps (cut away by the tail rule) and the 7 s sampling window inside each 10 s.

## 4. In plain words: what this means for real tails below 0.35 dB

- **The 0.3 dB rounding of RSL is not what makes the false step.** With no rounding at all the step is still kept
  in 5 of 5 replicates, with BIC gains of 40 to 130 (rule: more than 6). With 0.1 dB rounding it is the same. My
  explanation in the 3b-3c report was wrong.
- **The fits look the same in all three settings** (figure): a flat rate below W_c (under the truth in 13 of 15 fits, by up to 60 %) and a
  narrow spike of 16 to 75 % over the truth just above W_c. Neither is in the truth.
- **Where the step sits is not information either.** In 15 of 15 fits the fitted W_c is within 12 % of the value
  it started from (0.12, 0.2 or 0.4 dB). The "feature at 0.2-0.4 dB" of 3b-3c is largely the list of start values.
- **For real tails nothing changes, and one plan is dropped.** The Step 4 recipe of 3b-3c stays as it is: read the
  curve only for W >= 0.35 dB, no step tool, no choice of the law by held-out likelihood, nothing said about a switch.
  What is dropped: "put the quantization into the likelihood, then repeat B, C, D". This test gives no reason to
  expect that it would remove the step.
- **A step on real tails would still mean nothing.** The tool draws one for an exponential truth whether or not the
  data are rounded, so the real rounding cannot be used as the excuse or as the cure.
- **The curve without the step is not broken by the rounding.** For m1 its level is within 10 % of the truth above
  0.35 dB in every setting but one fit (+11 %), and its slope within 0.07 of the true 0. The rounding does make the
  data noisier (the BIC of the same replicate is 55 to 75 % higher at 0.3 dB than without rounding), but it does not
  bend the curve.
- **The cause is not known.** What remains is listed in section 5; none of it was tested here.

## 5. Unclear / possibly wrong

1. **The cause of the false step is still open; this step only removes one suspect.** Not tested: (a) the step fit
   gets 1150 more evaluations than the curve without the step, and they also move the spline values, so part of the
   "gain" may be a no-step fit that had not fully converged (check: continue the no-step fit for the same number
   of evaluations and compare; minutes per dataset); (b) the offset prior, the start-level floor (0.05 dB) and the
   per-tail start profile, which all act at low W; (c) the speed spread, which is marginalised with one width for
   all links; (d) the step itself, a logistic only 0.02 wide in ln W, to which the spline answers with the spike.
2. **The likelihood was left unchanged on purpose, so it is slightly wrong for two settings.** Its white-noise
   variance contains the term for a 0.3 dB step also when the data are not rounded or rounded to 0.1 dB (up to
   0.035 dB of extra assumed noise per minute). That makes the fit a little too tolerant there. I expect this to
   lower the gains, not to create a step, but it was not checked.
3. **Knots and penalty were chosen per dataset by the 3b-3b rule**, not copied from the 0.3 dB dataset. W_hi
   agrees within 0.03 dB between the settings of a replicate; the cross-validated penalty differs in 4 of 5
   replicates (column `lambda`). "Differ ONLY in quantization" holds for the data, not exactly for these two
   fit settings.
4. **The W_c half of the rule did not matter and could not have.** The median W_c at 0.1 dB is 1.07 x that at
   0.3 dB, but since W_c stays at its start value the ratio measures which start won, not a position found in the
   data. The verdict rests on "5 of 5 without quantization".
5. **One truth, 5 replicates.** Only m1 was tested. The 0.3 dB rows are the cached fits of 3b-3b / 3b-3c, not
   refitted (the re-simulated data are identical in 5 of 5 replicates).
6. **The real-data paragraph is a reading, not a measurement.** The counts come from the release table
   `sublinks`, read today; that the grid position is "arbitrary relative to W = 0" on real data is my argument
   (baseline subtraction), not something checked on tails.
7. **This says nothing about other effects of the rounding** (the level error of +-20 %, the too narrow bootstrap
   band of 3b-3c part D). Only the step was tested.
8. **The figure is diag; no key figure proposed.**
