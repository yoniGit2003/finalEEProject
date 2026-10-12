"""wet_antenna step 3b-4: is the 0.3 dB RSL quantization the cause of the false step of 3b-3c? Diagnosis only.

Step 3b-3c found that the free curve with an optional step keeps a step at 0.2-0.4 dB for EVERY fake truth, also the
pure exponential m1. The suspected cause was the 0.3 dB quantization of RSL, which the simulator applies and the
likelihood does not model. This step tests that directly, on fake data, with the rule written before the run:
  - truth m1 (exponential, no step), the skeleton, seeds (4000 + replicate), spread and noise of step 3b-3b;
  - three simulator settings that differ ONLY in the RSL quantization: none, 0.3 dB (= the 3b-3b data), 0.1 dB;
  - the free curve as in 3b-3b (same knot rule, starts, penalty by cross-validation), then the optional step exactly
    as in 3b-3c part B (W_c starts, budgets, kept if BIC improves by > 6).
Nothing is fixed here: the likelihood is the one of 3b-3b / 3b-3c, unchanged (it still assumes a 0.3 dB step in its
white-noise variance for all three settings). The 0.3 dB fits are read from the caches of 3b-3b and 3b-3c.
Settings and the rule: `wet_antenna.check_quantization` in research/config.yaml. Real data is used only as in 3b-3b (NaN minutes of the
real tails, dry-noise statistics). The frozen tail set v1 is read, never changed.

Replicate 0 is run first and timed. Cached in .cache/check_quantization/ (delete to redo). Run from the repository root (needs
steps 3b-3b and 3b-3c run first; about 20 min the first time, seconds when cached):
    python -m research.wet_antenna.check_quantization
"""
import os
os.environ.setdefault("NUMBA_NUM_THREADS", "1")            # one thread per process: the processes run in parallel

from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import json
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from research.wet_antenna import paths as P
from research.common.save import save
from research.wet_antenna import simulator as ws
from research.wet_antenna import pooled_fit as wp
from research.wet_antenna import free_curve as wc
from research.wet_antenna import stress_pooled_fit as st                          # event bootstrap, noise tables
from research.wet_antenna import check_free_curve as s3                   # load_tails, true_k, md, rel
from research.wet_antenna import confirm_free_curve as s3b                             # run_job: the free-curve fit of 3b-3b

CFG = s3.CFG
SIM, S, S3B, FIT = (CFG["wet_antenna"][k] for k in ("simulator", "check_quantization", "confirm_free_curve", "check_free_curve"))
B = CFG["wet_antenna"]["check_discovery_method"]["b"]                         # the step fit of 3b-3c part B, unchanged
DEC, TRUTH, W_LO = S["decision"], S["truth"], S3B["w_lo_db"]
SETTINGS = S["settings"]                                   # name -> RSL step (dB); 0 = not quantized
REF, FINE = S["reference"], S["fine"]
SPM = SIM["samples_per_min"]
INK, INK2 = "#1a1a1a", "#52514e"
CACHE = P.CACHE / "check_quantization"
CACHE_3B, CACHE_3C = P.CACHE / "confirm_free_curve", P.CACHE / "check_discovery_method"
md, rel = s3.md, s3.rel
tag_of = lambda name, rep: f"q{SETTINGS[name]:g}_{TRUTH}_{rep}"   # noqa: E731


# =================================================================================================== worker
def run_job(job):
    """ONE new dataset: the free-curve fit of 3b-3b (its own function, cached here), then the step fit of 3b-3c part
    B from that solution (same knots, same penalty, one start per W_c, the best continued)."""
    s3b.CACHE = CACHE
    r3 = s3b.run_job(job)
    f_ = CACHE / f"step_{job['tag']}.json"
    res = json.loads(f_.read_text()) if f_.exists() else {}
    if res.get("done"):
        return r3, res

    def save():
        tmp = f_.with_suffix(".tmp")
        tmp.write_text(json.dumps(res))
        tmp.replace(f_)

    t0 = time.time()
    _, noise, offset_sd = st.load_noise()
    D = wp.pack(s3.load_tails(job["sim"]), noise, offset_sd, S3B["window_min"], FIT["min_points"])
    K, x3, lam = len(r3["knots"]), np.asarray(r3["x"]), r3["lam_star"]
    sc = wc.StepCurve(r3["knots"], B["wc_range_db"], B["width_logw"], prior_sd=FIT["prior_sd_logk"])
    res.setdefault("step_starts", {})
    for wc0 in B["wc_starts_db"]:
        if str(wc0) not in res["step_starts"]:
            p0 = np.r_[x3[:-1], sc.u_of(wc0, B["height_start"]), x3[-1]]
            res["step_starts"][str(wc0)] = sc.fit(D, p0, lam, B["maxfev"]["start"])
            save()
    best = min(res["step_starts"], key=lambda n: res["step_starts"][n]["cost"])
    r = sc.fit(D, res["step_starts"][best]["x"], lam, B["maxfev"]["finish"])
    wc_db, height = sc.step(r["x"][:-1])
    r.update({"from": best, "wc_db": wc_db, "height": height, "bic": -2 * r["ll"] + (K + 3) * np.log(r3["n_independent"])})
    res.update({"step": r, "seconds_b": round(time.time() - t0, 1), "done": True})
    save()
    return r3, res


# =================================================================================================== figure
def figure(fits, par, summ):
    """One panel per quantization setting: the fit with the step of every replicate (solid, numbered), the true k
    (dashed), the fitted W_c of every replicate (short tick with its number at the bottom)."""
    fig, ax = plt.subplots(1, len(SETTINGS), figsize=(5.0 * len(SETTINGS), 4.8), sharey=True)
    for a_, name in zip(ax, SETTINGS):
        w_hi = min(fits[(name, r)][0]["w_hi"] for r in range(S3B["replicates"]))
        w = np.geomspace(0.05, w_hi, 240)
        for rep in range(S3B["replicates"]):
            r3, res = fits[(name, rep)]
            kept = r3["bic"] - res["step"]["bic"] > B["bic_gain"]
            sc = wc.StepCurve(r3["knots"], B["wc_range_db"], B["width_logw"])
            k = sc.k(res["step"]["x"][:-1], w)
            a_.plot(w, k, "-", color=INK, lw=1.5 if kept else 0.8,
                    label="learned, with the step (solid; thick = step kept; number = replicate)" if rep == 0 else None)
            a_.text(w[0] * 0.97, k[0], str(rep), fontsize=8, ha="right", va="center")
            a_.plot([res["step"]["wc_db"]] * 2, [0, 0.05], "-", color=INK, lw=1.2, transform=a_.get_xaxis_transform())
            a_.text(res["step"]["wc_db"], 0.06 + 0.035 * rep, str(rep), fontsize=7.5, ha="center",
                    transform=a_.get_xaxis_transform())
        kt = s3.true_k(TRUTH, par, w)
        a_.plot(w, kt, "--", color=INK, lw=2.2, label="true k(W) (dashed)")
        a_.text(w[-1], kt[-1], " true", fontsize=9, va="center")
        a_.axvline(W_LO, color=INK2, lw=0.8, ls=":")
        a_.text(W_LO, 0.99, f" {W_LO} dB", transform=a_.get_xaxis_transform(), fontsize=8, va="top")
        q = summ[summ.setting == name].iloc[0]
        a_.set(xscale="log", yscale="log", xlim=(0.04, w[-1] * 1.5), xlabel="W (dB)",
               title=f"RSL quantization: {name}\nstep kept in {q.step_kept}; median W_c of the kept steps "
                     f"{q.median_wc_kept_db} dB")
        a_.title.set_fontsize(9.5)
        a_.grid(alpha=0.25, which="both")
        for sp in ("top", "right"):
            a_.spines[sp].set_visible(False)
    ax[0].set_ylabel("k(W) = drying rate per dB of water (1 / min)")
    ax[0].legend(frameon=False, fontsize=7.5, loc="upper left")
    fig.suptitle(f"Step 3b-4, truth {TRUTH} ({ws.LAWS[TRUTH]}, no step in the truth): free curve with the optional step "
                 f"of 3b-3c; ticks at the bottom = fitted W_c per replicate", fontsize=10)
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    clock = time.time()
    CACHE.mkdir(parents=True, exist_ok=True)
    workers = max((os.cpu_count() or 2) - 1, 1)
    R = S3B["replicates"]
    new = [n for n in SETTINGS if n != REF]

    # ------------------------------------------------------------------------------------------- simulate (cached)
    cal = pd.read_csv(P.SIM_CALIBRATION).set_index("law")
    par = {TRUTH: {k: float(cal.loc[TRUTH, k]) for k in ws.PARAMS[TRUTH]}}
    sim_file = lambda name, rep: (CACHE_3B / f"sim_main_{TRUTH}_{rep}.parquet" if name == REF   # noqa: E731
                                  else CACHE / f"sim_{tag_of(name, rep)}.parquet")
    check_file = CACHE / "reference_check.json"
    ref_check = json.loads(check_file.read_text()) if check_file.exists() else {}
    if any(not sim_file(n, r).exists() for n in new for r in range(R)) or len(ref_check) < R:
        skel = pd.read_csv(P.TAILS_MODEL_INDEX).query("split_event == 'train'").sort_values("tail_id").reset_index(drop=True)
        real = pq.read_table(P.TAILS_MODEL, filters=[("split_event", "=", "train")],
                             columns=["tail_id", "minutes_since_t0_trace", "A_db"]).to_pandas()
        nan_real = {t: g.sort_values("minutes_since_t0_trace").A_db.isna().to_numpy() for t, g in real.groupby("tail_id")}
        _, noise, _ = st.load_noise()
        for rep in range(R):
            seed = S3B["seed_base"] + rep
            sk = st.bootstrap_skeleton(skel, np.random.default_rng(seed))
            masks = {r.tail_id: nan_real[r.orig_tail_id] for r in sk.itertuples()}
            f = np.exp(S3B["s_true"] * ws.spread_z(sk, seed))
            for name, step in SETTINGS.items():
                # the simulator draws the same random numbers whatever the step is; only the rounding changes
                sw = {**SIM["switches"], "quantization": step > 0}
                rows, _ = ws.simulate(sk, TRUTH, par[TRUTH], noise, sw, seed, f, step if step > 0 else SIM["rsl_step_db"],
                                      SPM, masks)
                if name == REF:                            # must reproduce the 3b-3b dataset exactly
                    old = pd.read_parquet(sim_file(name, rep), columns=["A_db"]).A_db.to_numpy()
                    ref_check[str(rep)] = bool(len(old) == len(rows) and np.allclose(old, rows.A_db, equal_nan=True, atol=0))
                else:
                    rows.merge(sk[["tail_id", "event_copy", "A_start_trace_db"]], on="tail_id").to_parquet(sim_file(name, rep))
        check_file.write_text(json.dumps(ref_check))
    print(f"datasets ready ({time.time() - clock:.0f} s); re-simulated {REF} data identical to the 3b-3b cache in "
          f"{sum(ref_check.values())} of {R} replicates", flush=True)

    def job_of(name, rep):                                 # the knot rule of 3b-3b, applied to this dataset
        a = pd.read_parquet(sim_file(name, rep), columns=["minutes_since_t0_trace", "A_db"])
        a = a.A_db[(a.minutes_since_t0_trace < S3B["window_min"]) & a.A_db.notna()].to_numpy()
        w_hi = float(np.sort(a)[-FIT["support_min_minutes"]])
        inside = np.log(a[(a >= W_LO) & (a <= w_hi)])
        q = np.quantile(inside, np.linspace(0, 1, S3B["knots"])[1:-1])
        return {"tag": tag_of(name, rep), "setting": name, "truth": TRUTH, "kind": "main", "rep": rep,
                "sim": str(sim_file(name, rep)), "knots": [W_LO] + [float(v) for v in np.exp(q)] + [w_hi], "w_hi": w_hi,
                "points": int(len(a)), "minutes_in_range": int(len(inside)), "menu_f4": False}

    def run(jobs, label):
        t_run, out = time.time(), {}
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futs = {pool.submit(run_job, j): j for j in jobs}
            for k, fu in enumerate(as_completed(futs)):
                j = futs[fu]
                out[(j["setting"], j["rep"])] = fu.result()
                print(f"   {label} {k + 1}/{len(jobs)}: {j['tag']} ({time.time() - t_run:.0f} s)", flush=True)
        return out, time.time() - t_run

    # ------------------------------------------------------------------------------------------- time replicate 0, run
    plan_file = CACHE / "plan.json"
    if plan_file.exists():
        plan = json.loads(plan_file.read_text())
    else:
        print(f"\ntiming replicate 0 ({len(new)} datasets) ...", flush=True)
        _, wall0 = run([job_of(n, 0) for n in new], "replicate 0")
        per_job = wall0 * min(workers, len(new)) / len(new)                 # process-seconds per dataset
        proj = wall0 + per_job * len(new) * (R - 1) / workers
        plan = {"wall_rep0_s": round(wall0, 1), "projected_s": round(proj), "wall_s": round(wall0, 1)}
        print(f"   replicate 0: {wall0 / 60:.1f} min wall; projected total {proj / 60:.0f} min on {workers} processes "
              f"(budget {S['time_budget_h'] * 60:.0f} min)", flush=True)
        plan_file.write_text(json.dumps(plan))
    fits, wall = run([job_of(n, r) for n in new for r in range(R)], "fits")
    if wall > 30:
        plan["wall_s"] = round(plan["wall_s"] + wall, 1)
        plan_file.write_text(json.dumps(plan))
    for rep in range(R):                                   # the 0.3 dB setting: the fits of 3b-3b and 3b-3c, as they are
        fits[(REF, rep)] = (json.loads((CACHE_3B / f"job_main_{TRUTH}_{rep}.json").read_text()),
                            json.loads((CACHE_3C / f"bc_{TRUTH}_{rep}.json").read_text()))

    # ------------------------------------------------------------------------------------------- tables
    aw = lambda w: wc.slope_weights(np.log(w))             # noqa: E731
    rows = []
    for name in SETTINGS:
        for rep in range(R):
            r3, res = fits[(name, rep)]
            stp, gain = res["step"], r3["bic"] - res["step"]["bic"]
            plain, x = wc.Curve(r3["knots"]), np.asarray(r3["x"])[:-1]
            w_lo, w_hi = np.geomspace(*S["low_range_db"], 40), np.geomspace(W_LO, r3["w_hi"], FIT["grid_points"])
            rat = lambda w: float(np.exp(np.median(np.log(plain.k(x, w) / s3.true_k(TRUTH, par, w)))))   # noqa: E731
            rows.append({"setting": name, "rep": rep, "step_kept": bool(gain > B["bic_gain"]), "wc_db": round(stp["wc_db"], 3),
                         "bic_gain": round(gain, 1), "height_ln": round(stp["height"], 2),
                         "bic_no_step": round(r3["bic"], 1), "bic_step": round(stp["bic"], 1), "best_wc_start": stp["from"],
                         "lambda": r3["lam_star"], "w_hi_db": round(r3["w_hi"], 2),
                         "no_step_k_ratio_below": round(rat(w_lo), 3), "no_step_k_ratio_above": round(rat(w_hi), 3),
                         "no_step_slope_above": round(float(np.log(plain.k(x, w_hi)) @ aw(w_hi)), 3),
                         "seconds": round(r3.get("seconds", 0) + res.get("seconds_b", 0), 0)})
    tab = pd.DataFrame(rows)
    summ = []
    for name, g in tab.groupby("setting", sort=False):
        kept = g[g.step_kept]
        summ.append({"setting": name, "step_kept": f"{len(kept)} of {R}", "n_kept": len(kept),
                     "median_wc_kept_db": round(float(kept.wc_db.median()), 3) if len(kept) else np.nan,
                     "wc_kept_db": " ".join(f"{v:.2f}" for v in kept.wc_db), "median_bic_gain": round(float(g.bic_gain.median()), 1),
                     "bic_gain_range": f"{g.bic_gain.min():.1f} to {g.bic_gain.max():.1f}",
                     "no_step_k_ratio_below": round(float(g.no_step_k_ratio_below.median()), 3),
                     "no_step_k_ratio_above": round(float(g.no_step_k_ratio_above.median()), 3)})
    summ = pd.DataFrame(summ)
    q = summ.set_index("setting")
    n_none, n_fine = int(q.loc["none", "n_kept"]), int(q.loc[FINE, "n_kept"])
    wc_ratio = q.loc[FINE, "median_wc_kept_db"] / q.loc[REF, "median_wc_kept_db"]
    fine_ok = bool(n_fine <= DEC["cause_max_kept"] or wc_ratio < DEC["fine_wc_ratio"])
    if n_none <= DEC["cause_max_kept"] and fine_ok:
        verdict = "quantization is the cause"
    elif n_none >= DEC["not_cause_min_kept"]:
        verdict = "not the cause"
    else:
        verdict = "partly / unclear"
    basis = (f"without quantization the step is kept in {n_none} of {R}; at {FINE} in {n_fine} of {R} with median W_c "
             f"{q.loc[FINE, 'median_wc_kept_db']} dB = {wc_ratio:.2f} x the {q.loc[REF, 'median_wc_kept_db']} dB at {REF}")
    dec = pd.concat([summ.drop(columns="n_kept").assign(verdict=""),
                     pd.DataFrame([{"setting": "verdict", "verdict": verdict, "step_kept": basis}])], ignore_index=True)
    summ = summ.drop(columns="n_kept")
    print(f"\nVERDICT: {verdict} ({basis})")

    save(tab, P.QUANT_TABLE, index=False)
    save(dec, P.QUANT_DECISION, index=False)
    fig = figure(fits, par, summ)
    save(fig, P.QUANT_FIG, dpi=110)
    plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    show = tab[["setting", "rep", "step_kept", "wc_db", "bic_gain", "height_ln", "best_wc_start", "lambda", "w_hi_db",
                "no_step_k_ratio_below", "no_step_k_ratio_above", "no_step_slope_above"]]
    report = f"""# wet_antenna step 3b-4: is the 0.3 dB quantization the cause of the false step? (rule fixed before the run)

Written by `research/wet_antenna/check_quantization.py` (settings `wet_antenna.check_quantization`). Fake data only; diagnosis only: no new
likelihood, no fix. No real tail is fitted, the frozen tail set v1 is not changed.

**Verdict: {verdict}.**

What was run: truth {TRUTH} ({ws.LAWS[TRUTH]}; no step in the truth), the event-bootstrap skeleton, seeds
({S3B['seed_base']} + replicate), noise and spread (s = {S3B['s_true']}) of step 3b-3b, {R} replicates, first
{S3B['window_min']} min. Three simulator settings that differ only in how the 10 s RSL values are rounded: not at all,
to 0.3 dB (the data of 3b-3b / 3b-3c), to 0.1 dB. The simulator draws the same random numbers in all three (re-simulated
0.3 dB data identical to the 3b-3b cache in {sum(ref_check.values())} of {R} replicates). Fit: the free curve of 3b-3b
(same knot rule, 4 starts, penalty by cross-validation), then the optional step of 3b-3c part B (W_c starts
{B['wc_starts_db']} dB, {B['maxfev']['start']} + {B['maxfev']['finish']} evaluations, kept if BIC(no step) - BIC(step) >
{B['bic_gain']:g}). The 0.3 dB rows are the cached fits of 3b-3b and 3b-3c, not refitted. The likelihood is unchanged
in all three settings.
Timing: replicate 0 took {plan['wall_rep0_s'] / 60:.1f} min wall on {workers} processes, projected total
{plan['projected_s'] / 60:.0f} min (budget {S['time_budget_h'] * 60:.0f} min); run time {plan['wall_s'] / 60:.0f} min
(cached in `.cache/check_quantization/`; this run: {time.time() - clock:.0f} s). Software: numpy / scipy / numba only.

## 1. Setting x replicate (`{rel(P.QUANT_TABLE)}`)

`bic_gain` = BIC(no step) - BIC(step); `step_kept` = gain > {B['bic_gain']:g}; `wc_db` = fitted step position;
`height_ln` = step height in ln units; `no_step_k_ratio_below` / `_above` = fitted / true k of the curve WITHOUT the
step, median over [{S['low_range_db'][0]}, {S['low_range_db'][1]}] dB and over [{W_LO} dB, W_hi] (information only);
`no_step_slope_above` = slope of log k vs log W above {W_LO} dB (true: 0).

{md(show)}

Per setting (`{rel(P.QUANT_DECISION)}`):

{md(summ)}

![m1 fits per quantization setting](../figures/{P.QUANT_FIG.name})

## 2. Verdict

Rule (fixed before the run): "quantization is the cause" if without quantization the step is kept in <=
{DEC['cause_max_kept']} of {R} AND at 0.1 dB the median W_c is < {DEC['fine_wc_ratio']} x the median W_c at 0.3 dB or the
step is kept in <= {DEC['cause_max_kept']} of {R}; "not the cause" if without quantization the step is still kept in >=
{DEC['not_cause_min_kept']} of {R}; otherwise "partly / unclear".

**{verdict[0].upper() + verdict[1:]}: {basis}.**

## 3. What is quantized in the real data

{REAL}

## 4. In plain words: what this means for real tails below 0.35 dB

{PLAIN}

## 5. Unclear / possibly wrong

{UNCLEAR}
"""
    save(report, P.QUANT_REPORT)
    print(f"report: {rel(P.QUANT_REPORT)} ({time.time() - clock:.0f} s)")


REAL = """Read from `pipeline/docs/methods.md`, `pipeline/docs/openmrg_summary.md` and the release table `sublinks`
(column `rsl_min_step_db`); nothing was fitted. **RSL** is reported in steps of 0.3 dB: the smallest RSL step is 0.3 dB
on 719 of the 728 sub-links (0.1 dB on one, no data on 8). **TSL** is reported in steps of 1 dB, but it is constant on
most sub-links (ATPC is off; `tsl_behaviour` = constant for 523 of 728), and a frozen tail ends just before its first
minute with a TSL step, so inside a tail only the RSL rounding acts. Both are sampled every **10 s**; the pipeline
stores them losslessly (L1) and the tails use the **1 min mean of the six 10 s values** of tl = tsl - rsl (L3), minus
the baseline. The simulator does the same: it rounds the 10 s value to a 0.3 dB grid and averages six values into one
minute. Its grid has a random position per tail; on real data the position is fixed per sub-link, but the baseline
is subtracted afterwards, so relative to W = 0 the position is also arbitrary. **So the fake simulator matches the real
data in what is quantized (RSL, 0.3 dB), where (10 s samples) and how it reaches the 1 min series.** Not in the
simulator: the 1 dB TSL steps (cut away by the tail rule) and the 7 s sampling window inside each 10 s."""
PLAIN = """- **The 0.3 dB rounding of RSL is not what makes the false step.** With no rounding at all the step is still kept
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
- **The cause is not known.** What remains is listed in section 5; none of it was tested here."""
UNCLEAR = """1. **The cause of the false step is still open; this step only removes one suspect.** Not tested: (a) the step fit
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
8. **The figure is diag; no key figure proposed.**"""


if __name__ == "__main__":
    main()
