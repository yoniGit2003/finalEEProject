"""Review of A7: five checks, written AFTER the A7 result and logged before this code is run (decisions.md, "Review of A7:
the five checks"). The A7 result is not changed: its cache is only read.

    python -m research.waa_ksurface.a7_review                 every part in order (a finished part is read from the cache)
    python -m research.waa_ksurface.a7_review --report        write the report from the cache
    python -m research.waa_ksurface.a7_review --one sel_mid   one part in this process (used by the main run)

What it does (the numbers of the review):
    (2) sel_mid   stage 1 as in A7 (k_wet estimated), then the drying SELECTION (K0 vs K1, L-dryhigh, margin 0.00291 dB)
                  twice: k_wet held at its gauge-input estimate and at its radar-input estimate. First a code check: the
                  selection at the gauge estimate must reproduce the committed A7 stage table, else STOP.
    (5) dry       the dry-path windows of the learning events: median measured A - c at 5, 10, 20, 30 min after the path
                  went dry, against the A7 K1 fit and the K0 fit of sel_mid, at both estimates. A check, not a selection.
    (4) sel_flo, sel_fhi   as sel_mid, with F through the lower / upper bootstrap ends of the E4 cells (rule A3f).
    (3) boot_<est>_<block>  the selection of (2) on the first 50 A7 bootstrap draws (seed 3900) per estimate; a drawn
                  event counts as often as it is drawn (fits, held-out score, paired interval).

Real TRAIN data of the 14 learning events (the A7 setting and flag). VAL and TEST are not opened.
"""
import pickle
import sys
import time

import numpy as np
import pandas as pd
import psutil

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a6_recovery as a6, a7_real_fit as a7, data, fake6, fp, ksurf, ladder, s1_separability as s1
from research.waa_ksurface import s3_recovery as s3, s3b_recovery as s3b, s3c_a7_variant as s3c, s3d_kwet_fixed as s3d
from research.waa_ksurface import paths as P

A6, S1, A7C, R7 = a6.A6, s1.S1, data.branch_config()["a7"], data.branch_config()["a7r"]
GB = s3b.GB
log = a6.log
CACHE = P.CACHE / "a7_review"
A7_CACHE = P.CACHE / "a7_real_fit"
ESTIMATES = ("gauge", "radar")
CURVES = {"mid": "F as measured (A7)", "flo": "F at the lower ends (I_lo)", "fhi": "F at the upper ends (I_hi)"}
W_REPORT = A7C["w_report_db"]


# =================================================================================================== the setting
def f_curves(which):
    """The measured curves (fp.measured) with F through the cell values ("mid"), or through I_lo / I_hi (rule A3f); P unchanged."""
    if which == "mid":
        M = fp.measured()
        return [M[b] for b in M]
    col = R7["f_ends"][which]
    c = pd.read_csv(P.E4_CELLS)
    c = c[c.eligible]
    out = []
    for band, g in c.groupby("band"):
        f = g[g.F_valid].sort_values("R_cell")
        p = g[g.P_valid].sort_values("R_cell")
        out.append(fp.from_values(band, f.R_cell, f[col], p.R_cell, p.s, 1.0 / np.maximum(f.I_se.to_numpy(), 1e-6) ** 2,
                                  1.0 / np.maximum(p.s_se.to_numpy(), 1e-6) ** 2))
    return out


def setting(which="mid"):
    S, curves, bands, band_idx = a7.real()
    assert bands == list(fp.measured())
    if which != "mid":
        curves = f_curves(which)
    fl = pickle.loads((A7_CACHE / "flag.pkl").read_bytes())    # the A7 flag (0 removed); not run again
    use = a7.use_of(S, fl)
    return S, curves, band_idx, use


def starts(S):
    E = A6["extension"]                                         # the extension starts of step 3c / A7 (s3c.procedure)
    return [np.full(len(i), E["f_start_inc_db"]) for i in S["iF"]], [E["p_start_share"] * i for i in S["iP"]]


# =================================================================================================== weights of drawn events
def event_counts(ev, w, n_ev):
    """How often each event counts: the largest weight of its link-events (all link-events of a drawn event share it)."""
    c = np.zeros(n_ev)
    np.maximum.at(c, ev, np.asarray(w, float))
    return c


def expand(S, N, cnt):
    """Per-event sums of ONE copy -> one entry per copy (a drawn event counts as separate events in the paired interval).
    With every count 1 this is S, N restricted to the used events, in the same order."""
    idx = np.nonzero(cnt > 0)[0]
    k = np.round(cnt[idx]).astype(int)
    return np.repeat(S[idx], k), np.repeat(N[idx], k)


# =================================================================================================== stage 2 with event weights
def stage_w(prob, le, w, wet, inc_F, inc_P, score, margin_db, tag=""):
    """Stage 2 of `ladder.stage` (select = drying, wetting K0 held at `wet`, rungs K0, K1, selector K0 vs K1), with event
    weights w (= use x how often the event was drawn). With every weight 1 it is ladder.stage itself (checked in sel_mid)."""
    lambdas = ladder.C["lambdas"]
    bin_le = np.repeat(np.arange(len(prob["le_n"])), prob["le_n"] // 5)
    w = np.asarray(w, float)
    fold = ladder.folds(le, w)
    ev = prob["event"]
    n_ev = int(ev.max()) + 1
    cnt = event_counts(ev, w, n_ev)
    rows, best, prev = [], {}, None
    for rung in s3.DRY_RUNGS:
        for lam in lambdas:
            m = ksurf.Model(prob, "K0", rung, lam)
            if not m.pen and lam != lambdas[0]:
                continue
            x0 = m.start(inc_F, inc_P) if prev is None else m.carry(*prev)
            x0[m.sl["gw"]], x0[m.sl["hw"]] = wet
            m.hold("gw", "hw")
            x = m.fit(x0, w)
            info = dict(m.last)
            S1_, N1_ = np.zeros(n_ev), np.zeros(n_ev)
            for f in range(fold.max() + 1):
                out = (fold == f) & (w > 0)
                xf = m.fit(x, w * (fold != f), maxiter=ladder.C["maxiter_fold"])
                b = out[bin_le] & score
                h = ladder._huber(m.simulate(xf)[1][b] - prob["target"][b], m.cfg["huber_delta_db"])
                S1_ += np.bincount(ev[bin_le[b]], weights=h, minlength=n_ev)
                N1_ += np.bincount(ev[bin_le[b]], minlength=n_ev).astype(float)
            Se, Ne = expand(S1_, N1_, cnt)
            cv = float(Se.sum() / Ne.sum())
            rows.append({"rung": rung, "lambda": lam, "train_loss": info["loss"], "cv_loss": cv, "iterations": info["nit"], "converged": info["converged"]})
            log(f"   {tag}drying   {rung} lambda {lam:g}: train {info['loss']:.5f}  cv {cv:.5f}  ({info['nit']} it)")
            if rung not in best or cv < best[rung]["cv"]:
                best[rung] = {"cv": cv, "lam": lam, "x": x, "model": m, "S": Se, "N": Ne}
        prev = (best[rung]["model"], best[rung]["x"])
    Sel = ladder.C["select"]

    def compare(richer, simpler):
        gain = ladder.cv_db(best[simpler]["S"], best[simpler]["N"]) - ladder.cv_db(best[richer]["S"], best[richer]["N"])
        lo, hi = ladder.paired_interval(best[simpler]["S"], best[richer]["S"], best[richer]["N"], Sel["boot_reps"], Sel["seed"], Sel["interval"])
        return gain, lo, hi

    rung, _, steps = s3.select_k0_k1(compare, margin_db)
    _, g, lo, hi, _ = steps["K1"]
    tab = pd.DataFrame(rows)
    tab["best_lambda"] = [best[r]["lam"] == l for r, l in zip(tab.rung, tab["lambda"])]
    tab["selected"] = (tab.rung == rung) & tab.best_lambda
    b = best[rung]
    log(f"   {tag}drying: selected {rung} (lambda {b['lam']:g}); K1 gain {g:+.5f} dB ({lo:+.5f} to {hi:+.5f})")
    return {"table": tab, "rung": rung, "lam": b["lam"], "gain": g, "lo": lo, "hi": hi,
            "tau_dry": dict(zip(W_REPORT, s3d.k_dry_at(b["model"], b["x"], W_REPORT))),
            "tau_dry_K1": dict(zip(W_REPORT, s3d.k_dry_at(best["K1"]["model"], best["K1"]["x"], W_REPORT))),
            "fits": {r: (best[r]["lam"], best[r]["x"]) for r in best}}


# =================================================================================================== (2), (4): the selection at both estimates
def stage1_and_prob(S, curves, band_idx, use):
    """Stage 1 of A7 (ladder.two_stage, exactly) and the stage-2 problem; 1 / k_wet at the gauge and the radar input."""
    D, le = S["D"], S["le"]
    start_F, start_P = starts(S)
    rain_g, ok5 = ladder.gauge_input(D)
    prob_g = ksurf.prepare(D, band_idx, curves, rain_g, 0)
    prob_g["w5"] = prob_g["w5"] * ok5
    st1 = ladder.stage(prob_g, le, use * S["gauge"], "wetting", "K1", None, start_F, start_P, False, log=log, tag="(review) ", rungs=s3.WET_RUNGS)
    m1, x1 = st1["model"], st1["x"]
    assert st1["rung"] == "K0"
    gw, hw = x1[m1.sl["gw"]].copy(), x1[m1.sl["hw"]].copy()
    ext = m1.parts(x1)[4:]
    del prob_g, rain_g, ok5
    used = (use > 0)[S["bin_le"]]
    counted, _, high, _ = s1.score_sets(S, D)
    prob_r = ksurf.prepare(D, band_idx, curves, D["rain"], 0)
    prob_r["w5"] = np.asarray(counted & used, float)
    tau_r = ladder.plain_k0(prob_r, use, start_F, start_P, False)[0]
    tau_g = float(1 / m1.k(x1, "wetting", 0, 0))
    log(f"   (review) 1 / k_wet: gauge input {tau_g:.4f} min, radar input {tau_r:.4f} min")
    wet = {"gauge": (gw, hw), "radar": (np.full_like(gw, np.log(1.0 / tau_r)), hw)}
    return {"prob": prob_r, "score": high & used, "wet": wet, "ext": ext, "tau_wet": {"gauge": tau_g, "radar": tau_r}}


def sel_part(which):
    S, curves, band_idx, use = setting(which)
    st = stage1_and_prob(S, curves, band_idx, use)
    out = {"which": which, "tau_wet": st["tau_wet"], "est": {}}
    for e in ESTIMATES:
        r = stage_w(st["prob"], S["le"], use, st["wet"][e], [v.copy() for v in st["ext"][0]], [v.copy() for v in st["ext"][1]], st["score"],
                    s3c.margin(), tag=f"(review {which}, k_wet {e}) ")
        out["est"][e] = r
        if which == "mid" and e == "gauge":                     # the code check: A7's own selection must come back
            ok, text = reproduces(r)
            out["reproduces"], out["reproduce_text"] = ok, text
            log(f"=== (review) code check against the A7 stage table: {text} ===")
            if not ok:
                return out
    return out


def reproduces(r):
    """The selection at the gauge estimate, all weights 1, against the committed A7 drying table."""
    a = pd.read_csv(P.A7_CV)
    a = a[a["select"] == "drying"].reset_index(drop=True)
    t = r["table"]
    if len(a) != len(t) or list(a.rung) != list(t.rung) or not np.allclose(a["lambda"], t["lambda"]):
        return False, "the rung / lambda rows differ"
    dev = max(np.max(np.abs(t.train_loss.to_numpy() / a.train_loss.to_numpy() - 1)), np.max(np.abs(t.cv_loss.to_numpy() / a.cv_loss.to_numpy() - 1)))
    k1 = a[(a.rung == "K1") & a.best_lambda].iloc[0]
    dg = abs(r["gain"] - k1.gain_db) / abs(k1.gain_db)
    ok = bool(dev <= R7["reproduce_rtol"] and dg <= R7["reproduce_rtol"] and r["rung"] == a[a.selected].rung.iloc[0])
    return ok, f"{'REPRODUCED' if ok else 'NOT REPRODUCED'} (largest relative difference of the losses {dev:.1e}, of the K1 gain {dg:.1e}; rung {r['rung']})"


# =================================================================================================== (3): the selection bootstrap
def draws():
    """The first sel_draws draws of the A7 bootstrap (the same generator and sequence as a7_real_fit.boot_part)."""
    S, _, _, use = setting("mid")
    ev_all = S["le"].event_id.to_numpy()
    events = np.sort(np.unique(ev_all[use > 0]))
    rng = np.random.default_rng(A7C["boot_seed"])
    return [rng.choice(events, len(events), replace=True) for _ in range(R7["sel_draws"])]


def boot_part(est, block):
    S, curves, band_idx, use = setting("mid")
    sel = pickle.loads((CACHE / "sel_mid.pkl").read_bytes())
    assert sel.get("reproduces")
    st = stage1_and_prob(S, curves, band_idx, use)
    for e in ESTIMATES:                                         # stage 1 is deterministic: it must give sel_mid's k_wet again
        assert abs(st["tau_wet"][e] / sel["tau_wet"][e] - 1) < 1e-9
    ev_all = S["le"].event_id.to_numpy()
    picks = draws()
    lo, hi = block * R7["block"], min((block + 1) * R7["block"], R7["sel_draws"])
    t0 = time.time()
    for b in range(lo, hi):
        f = CACHE / f"draw_{est}_{b:03d}.pkl"
        if f.exists():
            continue
        cnt = pd.Series(picks[b]).value_counts()
        w = use * pd.Series(ev_all).map(cnt).fillna(0).to_numpy()
        r = stage_w(st["prob"], S["le"], w, st["wet"][est], [v.copy() for v in st["ext"][0]], [v.copy() for v in st["ext"][1]], st["score"],
                    s3c.margin(), tag=f"(draw {b} {est}) ")
        r.pop("fits")
        r["est"], r["draw"], r["events_distinct"] = est, b, int(len(cnt))
        f.write_bytes(pickle.dumps(r))
        log(f"   (review) [k_wet {est}] draw {b + 1} of {R7['sel_draws']}: {r['rung']} ({(time.time() - t0) / 60:.1f} min in this block)")
    return {"est": est, "block": block, "done": True}


# =================================================================================================== (5): the dry-path check
def dry_values(D, rain1, lags, rain_before_min, dry_below, use):
    """Per used link-event and lag: the index of the 5 min bin that ends `lag` min after the path went dry, or -1.
    rain1 = the 1 min radar input in the window layout (D["y1"]); the path went dry at the end of bin le_last if the next
    bin is radar-dry; at least rain_before_min wet minutes before; every bin up to the one used is radar-dry, counted, finite."""
    n5 = D["le_n"] // 5
    rain5 = rain1.reshape(-1, 5).mean(axis=1)
    target = D["y5"] - np.repeat(D["le_c"], n5)
    out = np.full((len(n5), len(lags)), -1, np.int64)
    for e in np.nonzero(use > 0)[0]:
        o1, o5, last = D["le_boff1"][e], D["le_boff5"][e], int(D["le_last"][e])
        t0 = (last + 1) * 5                                     # minute of the window at which the path went dry
        if last + 1 >= n5[e] or last + 1 < 0:
            continue
        if (rain1[o1:o1 + t0] >= dry_below).sum() < rain_before_min:
            continue
        for i, lag in enumerate(lags):
            j = last + lag // 5                                 # the bin that ends `lag` min after t0
            if j >= n5[e]:
                continue
            seg = slice(o5 + last + 1, o5 + j + 1)
            if np.all(rain5[seg] < dry_below) and D["w5"][o5 + j] > 0 and np.isfinite(target[o5 + j]):
                out[e, i] = o5 + j
    return out, target


def dry_part():
    S, curves, band_idx, use = setting("mid")
    D = S["D"]
    sel = pickle.loads((CACHE / "sel_mid.pkl").read_bytes())
    fit = pickle.loads((A7_CACHE / "fit.pkl").read_bytes())
    prob = ksurf.prepare(D, band_idx, curves, D["rain"], 0)
    lags = R7["dry_lags_min"]
    idx, target = dry_values(D, fake6.window(D, np.asarray(D["rain"], float)), lags, R7["rain_before_min"], S1["dry_below_mmh"], use)
    rows = []
    for e in ESTIMATES:
        k0_lam, k0_x = sel["est"][e]["fits"]["K0"]
        out_k1 = ksurf.Model(prob, "K0", fit["rung"], fit["lam"]).simulate(fit["x"][e])[1]          # A7's fit at this estimate
        out_k0 = ksurf.Model(prob, "K0", "K0", k0_lam).simulate(k0_x)[1]                           # the K0 refit (stage 2 of sel_mid)
        for i, lag in enumerate(lags):
            b = idx[:, i][idx[:, i] >= 0]
            meas, m1, m0 = np.median(target[b]), np.median(out_k1[b]), np.median(out_k0[b])
            rows.append({"k_wet estimate": e, "lag (min)": lag, "values": int(len(b)), "measured A - c (dB)": meas, f"{fit['rung']} (A7 fit) (dB)": m1,
                         "K0 refit (dB)": m0, f"error {fit['rung']} (dB)": m1 - meas, "error K0 (dB)": m0 - meas})
    return pd.DataFrame(rows)


# =================================================================================================== report
def f3(d):
    return ", ".join(f"W {W:g}: {v:.1f}" for W, v in d.items())


def report(R, state):
    fit = pickle.loads((A7_CACHE / "fit.pkl").read_bytes())
    m = s3c.margin()
    L = ["# Review of A7: the five checks\n",
         "Written by `research/waa_ksurface/a7_review.py` (decisions.md, \"Review of A7: the five checks\", logged before the run). The A7 result is not "
         "changed; its cache is read only. Real TRAIN data of the 14 learning events (A7 flag: 0 removed). VAL and TEST were not opened.\n",
         "## (1) Order\n", "See decisions.md: the addendum (f58a0f5 19:51:30, a6e868a 19:52:02) was committed before the first real read of A7 (flag "
         "launcher 19:52:08, flag.pkl 19:52:20). No deviation.\n"]
    rows, cvs = [], []
    for which in CURVES:
        r = R.get(f"sel_{which}")
        if r is None:
            continue
        for e, s in r["est"].items():
            rows.append({"F curve": CURVES[which], "k_wet estimate": e, "1/k_wet (min)": round(r["tau_wet"][e], 2), "drying rung": s["rung"],
                         "K1 gain (dB)": round(s["gain"], 5), "paired 90 % interval": f"{s['lo']:+.5f} to {s['hi']:+.5f}", "margin (dB)": round(m, 5),
                         **{f"1/k_dry at W = {W:g} dB (min)": round(v, 1) for W, v in s["tau_dry"].items()}})
            cvs.append(s["table"].assign(f_curve=which, k_wet_estimate=e))
    if rows:
        tab = pd.DataFrame(rows)
        save(tab, P.A7R_RUNGS, index=False)
        save(pd.concat(cvs, ignore_index=True), P.A7R_CV, index=False)
    mid = R.get("sel_mid")
    if mid is not None:
        L += ["## (2) The rung at each k_wet estimate\n", f"Code check: {mid.get('reproduce_text')}.\n"]
        if len(mid["est"]) == 2:
            rg = {e: mid["est"][e]["rung"] for e in ESTIMATES}
            verdict = rg["gauge"] if rg["gauge"] == rg["radar"] else "rung not determined"
            L += [md_table(tab[tab["F curve"] == CURVES["mid"]]), "",
                  f"**Gauge estimate: {rg['gauge']}; radar estimate: {rg['radar']}. Result (pre-registered rule): {verdict}.** Caveat stated with it "
                  "(Deviation 15): the selection with k_wet held at one value is the procedure of step 3d, which gave a false K1 on truth (a); it is "
                  "run here because the review asks for it.\n"]
    boot = sorted(CACHE.glob("draw_*.pkl"))
    if boot:
        B = pd.DataFrame([{"k_wet estimate": r["est"], "draw": r["draw"], "events distinct": r["events_distinct"],
                           "rung": r["rung"], "K1 gain (dB)": r["gain"], "lo": r["lo"], "hi": r["hi"],
                           **{f"1/k_dry at W = {W:g} dB (min)": v for W, v in r["tau_dry"].items()}}
                          for f in boot for r in [pickle.loads(f.read_bytes())]])
        save(B, P.A7R_BOOT, index=False)
        sh = []
        for e in ESTIMATES:
            b = B[B["k_wet estimate"] == e]
            if not len(b):
                continue
            a7b = pickle.loads((A7_CACHE / f"boot_{e}.pkl").read_bytes())
            sh.append({"k_wet estimate": e, "draws": len(b), "select K1": int((b.rung == "K1").sum()), "share K1": round(float((b.rung == "K1").mean()), 2),
                       "K1 gain, median (5-95 %) (dB)": f"{b['K1 gain (dB)'].median():+.4f} ({b['K1 gain (dB)'].quantile(.05):+.4f} to {b['K1 gain (dB)'].quantile(.95):+.4f})",
                       **{f"1/k_dry W = {W:g} dB, A7 200 draws 5-95 % (min)": f"{np.percentile(a7b[f'1/k_dry at W = {W:g} dB'], 5):.1f} to "
                          f"{np.percentile(a7b[f'1/k_dry at W = {W:g} dB'], 95):.1f}" for W in (0.5, 1.0, 2.0)}})
        L += ["## (3) Event bootstrap of the selection\n",
              f"The first {R7['sel_draws']} draws of the A7 bootstrap (seed {A7C['boot_seed']}) per estimate (user's decision; 200 with the full "
              "selection would take about 40 h), the full stage-2 selection in each draw. The 90 % intervals of 1 / k_dry are those of the 200-draw A7 "
              f"bootstrap (refits of the selected rung), restated. Every draw: `{rel(P.A7R_BOOT)}`.\n", md_table(pd.DataFrame(sh)), ""]
    if any(R.get(f"sel_{w}") is not None for w in ("flo", "fhi")):
        L += ["## (4) Sensitivity to F (rule A3f)\n", "F through the 5-95 % bootstrap ends of the E4 cells (I_lo, I_hi); P unchanged; the whole of (2) at "
              "each end.\n", md_table(tab), ""]
    if "dry" in R:
        dt = R["dry"]
        save(dt, P.A7R_DRY, index=False)
        L += ["## (5) Dry-path check (a check, not a selection)\n",
              f"Link-events of the learning events; the path went dry at the end of the last wet bin if the next bin is radar-dry (< "
              f"{S1['dry_below_mmh']} mm/h), at least {R7['rain_before_min']} wet minutes before; the 5 min bin that ends at the lag, the path dry "
              "throughout. Error = median(model) - median(measured).\n", md_table(dt.round(3)), ""]
    if state:
        L.append(f"Execution: keep-awake {state.get('awake')}; lowest free memory {state['min_free']:.2f} GB; low-memory events {len(state['low'])}.\n")
        if state.get("incomplete"):
            L.append(f"\n**NOT COMPLETE: {state['incomplete']}. Run again to go on from the cache.**\n")
    save("\n".join(L) + "\n", P.A7R_REPORT)
    log(f"report: {rel(P.A7R_REPORT)}")


# =================================================================================================== main
def parts():
    nb = -(-R7["sel_draws"] // R7["block"])
    return ["sel_mid", "dry", "sel_flo", "sel_fhi"] + [f"boot_{e}_{k}" for e in ESTIMATES for k in range(nb)]


def run_part(name):
    if name.startswith("sel_"):
        return sel_part(name[4:])
    if name == "dry":
        return dry_part()
    _, est, k = name.split("_")
    return boot_part(est, int(k))


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    if "--one" in sys.argv:
        name = sys.argv[sys.argv.index("--one") + 1]
        t0 = time.time()
        r = run_part(name)
        (CACHE / f"{name}.pkl").write_bytes(pickle.dumps(r))
        log(f"=== (review {name}) done: {(time.time() - t0) / 60:.1f} min; peak memory {psutil.Process().memory_info().peak_wset / GB:.2f} GB ===")
        return
    state = None
    R = {}
    if "--report" not in sys.argv:
        state = s3d.start_state("A7 review")

        def rule(name, r):
            return "the code check did not reproduce A7" if name == "sel_mid" and not r.get("reproduces") else ""
        R, stop = s3d.run_children("research.waa_ksurface.a7_review", CACHE, parts(), state, lambda n: "", rule)
        if stop:
            state["incomplete"] = f"STOP: {stop}"
            log(f"=== STOP: {stop} ===")
    else:
        R = {n: pickle.loads((CACHE / f"{n}.pkl").read_bytes()) for n in parts() if (CACHE / f"{n}.pkl").exists()}
    report(R, state)
    log(f"=== A7 review ended {time.strftime('%Y-%m-%d %H:%M:%S')} ===")


if __name__ == "__main__":
    main()
