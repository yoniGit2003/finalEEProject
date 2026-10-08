"""R1b: event-level validation of baseline v2, the humidity-exclusion decision, PASS/FAIL, and real-event coverage.

On the C1 fake events (certified dry weather, truth A = 0): offset_j = MEAN of e(t) = tl_mean - baseline over fake
event j, i.e. the baseline error of that event. Variants: v1 (A/B, 6 h mean), M0 (summer median), M1 (24 h median),
v2 with and without the humidity exclusion. The adoption rule (fixed before the run) decides the exclusion; PASS/FAIL
and the real-event part use v2 as adopted.

Run after c1_fake_events.py and c3_daily_cycle.py: python research/baseline_estm/validation/r1b_event_level.py
"""
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from common import (V, B, OUT, FIG, HERE, SPH, Timer, open_data, sublink_table, groups, load_panel, usable, window,
                    stat, m3_baseline)
from research.common.outputs import save_output
from c2_after_window import deltas, curve

sys.path.insert(0, str(HERE.parent))
import baseline as bl                                                          # noqa: E402  (the method)

R = V["r1b"]
VARS = ["v1", "M0", "M1", "v2_hx", "v2_nohx"]
LABEL = {"v1": "v1 (A/B, 6 h mean)", "M0": "M0 (summer median)",
         "M1": "M1 (24 h median; same as v2 without exclusion when the 48 h fallback is not needed)",
         "v2_hx": "v2 with humidity exclusion", "v2_nohx": "v2 without humidity exclusion"}


def offsets(sl, tl, state, met, dew):
    """Baseline and offset of every variant for every C1 fake event of the given sub-links."""
    ev = pd.read_csv(OUT / "c1_fake_events.csv.gz", parse_dates=["t0"])
    ev = ev[ev.sublink_id.isin(sl.sublink_id)]
    min_n, w1 = int(round(V["variants"]["min_dry_h"] * SPH)), int(V["variants"]["M1_window_h"] * SPH)
    rows, slots = [], {(v, D): [] for v in VARS for D in V["c1"]["durations_h"]}
    for s, g in ev.groupby("sublink_id"):
        c = int(s.split("_")[0])
        x, st = tl[s].to_numpy(), state[c].to_numpy()
        ok = usable(state[c], V["skip_h"])                                     # same as bl.usable_mask(st, 1 h)
        okhx = bl.usable_mask(st, B["skip_h"], dew)
        m0 = float(np.median(x[ok & ~np.isnan(x)]))
        for e in g.itertuples():
            n = int(e.D_h * SPH)
            i0, seg = e.i0, x[e.i0:e.i0 + n]
            base = {"v1": m3_baseline(x, ok, m0, i0, i0 + n)[0], "M0": np.full(n, m0),
                    "M1": np.full(n, stat(window(x, ok, i0 - w1, i0), "median", min_n)),
                    "v2_hx": np.full(n, bl.baseline_v2(x, okhx, i0, B)[0]),
                    "v2_nohx": np.full(n, bl.baseline_v2(x, ok, i0, B)[0])}
            r = {"sublink_id": s, "D_h": e.D_h, "t0": e.t0, "i0": i0}
            for v in VARS:
                err = seg - base[v]
                r[f"off_{v}"] = float(np.nanmean(err)) if np.isfinite(err).any() else np.nan
                r[f"base_{v}"] = base[v][0]
                slots[(v, e.D_h)].append(err)
            rows.append(r)
    o = pd.DataFrame(rows).merge(sl, on="sublink_id", how="left")
    o["rh_t0"] = met.humidity.reindex(o.t0).to_numpy()
    o["humid_start"] = o.rh_t0 >= R["humid_start_rh_pct"]
    return o, slots


def table(o, slots, v, by=None):
    """Event-level statistics of one variant per D (and per `by` value)."""
    rows = []
    keys = [None] if by is None else sorted(o[by].dropna().unique())
    for k in keys:
        ok_ = o if k is None else o[o[by] == k]
        for D, g in ok_.groupby("D_h"):
            x = g[f"off_{v}"]
            r = {"variant": v, "D_h": D, "n_events": len(g), "missing_share": x.isna().mean(),
                 "median_offset": x.median(), "p10": x.quantile(0.1), "p90": x.quantile(0.9),
                 "p90_abs": x.abs().quantile(0.9)}
            if by is None:
                e = np.concatenate(slots[(v, D)])
                r["slot_p10"], r["slot_p90"] = np.nanpercentile(e, 10), np.nanpercentile(e, 90)
            else:
                r = {by: k, **r}
            rows.append(r)
    return pd.DataFrame(rows)


def adoption(o):
    """The rule fixed before the run: adopt the humidity exclusion if p90 |offset| at humid starts improves by
    >= min_gain_db AND the overall missing share rises by <= max_missing_rise."""
    h = o[o.humid_start]
    p_no, p_hx = h.off_v2_nohx.abs().quantile(0.9), h.off_v2_hx.abs().quantile(0.9)
    m_no, m_hx = o.off_v2_nohx.isna().mean(), o.off_v2_hx.isna().mean()
    gain, rise = p_no - p_hx, m_hx - m_no
    adopt = bool(gain >= R["adopt_min_gain_db"] and rise <= R["adopt_max_missing_rise"])
    return pd.DataFrame([{"humid_start_events": len(h), "p90_abs_humid_without": p_no, "p90_abs_humid_with": p_hx,
                          "gain_db": gain, "missing_without": m_no, "missing_with": m_hx, "missing_rise": rise,
                          "rule_min_gain_db": R["adopt_min_gain_db"],
                          "rule_max_missing_rise": R["adopt_max_missing_rise"], "adopt": adopt}])


def pass_fail(o, v):
    P = R["pass_fail"]
    rows = []
    for s, g in o.groupby("sublink_id"):
        sh, lg = g[g.D_h <= V["c1"]["short_max_h"]][f"off_{v}"].dropna(), g[g.D_h == V["c1"]["long_h"]][f"off_{v}"].dropna()
        r = {"sublink_id": s, "n_short": len(sh), "n_24h": len(lg),
             "median_offset_short": sh.median() if len(sh) else np.nan,
             "p90_abs_short": sh.abs().quantile(0.9) if len(sh) else np.nan,
             "p90_abs_24h": lg.abs().quantile(0.9) if len(lg) else np.nan}
        r["c1_median"] = abs(r["median_offset_short"]) <= P["median_abs_db"] if len(sh) else np.nan
        r["c2_p90"] = r["p90_abs_short"] <= P["p90_abs_db"] if len(sh) else np.nan
        r["c3_growth"] = (r["p90_abs_24h"] - r["p90_abs_short"]) <= P["extra_p90_24h_db"] if len(sh) and len(lg) else np.nan
        rows.append(r)
    pf = pd.DataFrame(rows).merge(o.drop_duplicates("sublink_id")[["sublink_id", "length_class", "frequency_band",
                                                                    "tsl_behaviour", "daily_cycle_large"]], on="sublink_id")
    cs = ["c1_median", "c2_p90", "c3_growth"]
    pf["evaluable"] = pf[cs].notna().all(axis=1)
    pf["pass_all"] = pf.evaluable & (pf[cs] == True).all(axis=1)                     # noqa: E712
    names = {"c1_median": f"|median offset| > {P['median_abs_db']}", "c2_p90": f"p90 |offset| > {P['p90_abs_db']}",
             "c3_growth": f"p90 |offset| growth at 24 h > {P['extra_p90_24h_db']}"}
    pf["fail_reason"] = ["; ".join(names[c] for c in cs if r[c] is False or r[c] == 0) for _, r in pf.iterrows()]
    return pf


def pf_summary(pf):
    cs = ["c1_median", "c2_p90", "c3_growth"]
    out = []
    for split, by in (("all", None), ("length / frequency", ["length_class", "frequency_band"]),
                      ("TSL class", "tsl_behaviour"), ("daily_cycle_large", "daily_cycle_large")):
        for k, g in ([("all", pf)] if by is None else pf.groupby(by)):
            r = {"split": split, "group": " / ".join(map(str, k)) if isinstance(k, tuple) else str(k), "n": len(g)}
            for c in cs:
                r[f"pass_{c}"] = g[c].astype(float).mean()
            r["pass_all"] = g[g.evaluable].pass_all.mean()
            out.append(r)
    return pd.DataFrame(out)


def real_events(d, sl, tl, state, dew, amp, use_hx):
    """v2 (as adopted) on significant train events x train sub-links: baseline, flags, Delta_k."""
    ev = d.events(significant=True, split="train")
    T = tl.index
    rows = []
    for s in sl.sublink_id:
        c = int(s.split("_")[0])
        x, st = tl[s].to_numpy(), state[c].to_numpy()
        ok_plain = bl.usable_mask(st, B["skip_h"])
        ok = bl.usable_mask(st, B["skip_h"], dew) if use_hx else ok_plain
        summer = float(np.median(x[ok_plain & ~np.isnan(x)]))
        info = sl.set_index("sublink_id").loc[s]
        wet = st == 1
        for e in ev.itertuples():
            i0, i1 = T.searchsorted(pd.Timestamp(e.start)), T.searchsorted(pd.Timestamp(e.end))
            value, w, a = bl.baseline_v2(x, ok, i0, B)
            r = {"sublink_id": s, "event_id": e.event_id, "baseline_db": value, "window_h": w,
                 "missing": bool(np.isnan(value)),
                 **bl.flags_v2(value, a, i0, st, ok_plain, dew, summer, amp.get(s, np.nan), info.is_degraded,
                               info.is_toggling, B)}
            wp = np.nonzero(wet[i0:i1])[0]
            if len(wp) and not np.isnan(value):
                end_pos = i0 + wp[-1] + 1
                r.update({"duration_h": (end_pos - i0) / SPH, **deltas(x, ok_plain, value, end_pos)})
            rows.append(r)
    return pd.DataFrame(rows).merge(sl, on="sublink_id", how="left")


def plots(tab, o, cv):
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
    for v, col in zip(VARS, ["C7", "C2", "C1", "C0", "C3"]):
        t = tab[tab.variant == v]
        axes[0].plot(t.D_h, t.p90_abs, "o-", color=col, label=LABEL[v])
        axes[1].plot(t.D_h, t.median_offset, "o-", color=col, label=LABEL[v])
        axes[1].fill_between(t.D_h, t.p10, t.p90, color=col, alpha=0.07)
    for ax, lab in zip(axes, ("p90 of |offset| [dB]", "median offset (band p10-p90) [dB]")):
        ax.set_xscale("log")
        ax.set_xticks(V["c1"]["durations_h"], [str(D) for D in V["c1"]["durations_h"]])
        ax.set_xlabel("fake event duration D [h]")
        ax.set_ylabel(lab)
    axes[1].axhline(0, color="k", lw=0.8)
    axes[0].axhline(R["pass_fail"]["p90_abs_db"], color="k", ls=":", lw=1, label="PASS limit 0.25 dB")
    axes[0].legend(fontsize=8)
    fig.suptitle("R1b: event-level baseline error (offset = mean e over the fake event), main sub-links")
    fig.tight_layout()
    save_output(fig, FIG / "r1b_offset_by_D.png", "key", "v2 vs v1 / M0 / M1 event offsets by duration", dpi=110)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    k = cv.k_h
    ax.fill_between(k, cv.real_p10, cv.real_p90, color="C3", alpha=0.15)
    ax.plot(k, cv.real_median, "o-", color="C3", label="real events (v2 before)")
    ax.fill_between(k, cv.ctrl_p10, cv.ctrl_p90, color="C0", alpha=0.15)
    ax.plot(k, cv.ctrl_median, "s--", color="C0", label="fake control (v2 before)")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xscale("log")
    ax.set_xticks(V["c2"]["k_h"], [str(v) for v in V["c2"]["k_h"]])
    ax.set_xlabel("k = hours after the end of the last wet slot")
    ax.set_ylabel("Delta_k = level_k - v2 baseline [dB]")
    ax.legend(fontsize=8)
    ax.set_title("R1b: after-window level vs v2 baseline (median, band p10-p90)")
    fig.tight_layout()
    save_output(fig, FIG / "r1b_c2_delta.png", "key", "Delta_k after real rain with v2 vs control", dpi=110)
    plt.close(fig)


def main():
    with Timer("R1b"):
        d = open_data()
        sl = sublink_table(d)
        tl, state, cert, met = load_panel(d, sl)
        dew = bl.dew_risk(met.humidity, met.dpd, B["dew"]) & met.humidity.notna().to_numpy()
        amp = pd.read_csv(OUT / "c3_daily_cycle_per_sublink.csv").set_index("sublink_id").amplitude_db
        main_sl = sl[sl.group == "main"].assign(
            daily_cycle_large=lambda x: x.sublink_id.map(amp) > B["daily_cycle_large_db"])
        o, slots = offsets(main_sl, tl, state, met, dew)
        save_output(o, OUT / "r1b_fake_event_offsets.csv.gz", "product",
                    "R1b event offsets of every method (read by make_baseline_table)", index=False)
        tab = pd.concat([table(o, slots, v) for v in VARS])
        save_output(tab, OUT / "r1b_offset_by_variant_D.csv", "diag", "R1b offsets per variant and D", index=False)
        o["cls"] = o.length_class + " / " + o.frequency_band
        splits = pd.concat([table(o, slots, v, by).assign(split=by) for v in ("v2_hx", "v2_nohx")
                            for by in ("humid_start", "daily_cycle_large", "cls")])
        save_output(splits, OUT / "r1b_splits.csv", "diag",
                    "R1b offsets split by humidity / split / class", index=False)
        dec = adoption(o)
        save_output(dec, OUT / "r1b_humidity_exclusion_decision.csv", "diag",
                    "R1b humidity exclusion decision", index=False)
        adopted = "v2_hx" if dec.adopt.iloc[0] else "v2_nohx"
        if B["v2"]["humidity_exclusion"] != dec.adopt.iloc[0]:
            print(f"NOTE: research/config.yaml baseline.v2.humidity_exclusion = {B['v2']['humidity_exclusion']}, "
                  f"the rule says {dec.adopt.iloc[0]} -> update the config")
        pf = pass_fail(o, adopted)
        save_output(pf, OUT / "r1b_pass_fail_per_sublink.csv", "diag", "R1b PASS / FAIL per sub-link", index=False)
        pfs = pf_summary(pf)
        save_output(pfs, OUT / "r1b_pass_fail_summary.csv", "diag", "R1b PASS / FAIL counts", index=False)
        # real events: all train sub-links (good + degraded), numbers reported for good (main)
        re_ = real_events(d, sl, tl, state, dew, amp, adopted == "v2_hx")      # flags come from bl.flags_v2
        save_output(re_, OUT / "r1b_real_events.csv.gz", "diag",
                    "R1b real (event, sub-link) pairs with v2", index=False)
        cov = []
        for gname, g in (("good (main)", re_[re_.group == "main"]), ("degraded", re_[re_.is_degraded]),
                         ("toggling", re_[re_.is_toggling])):
            cov.append({"group": gname, "pairs": len(g), "with_baseline": 1 - g.missing.mean(),
                        "missing": g.missing.mean(), "fallback_48h": (g.window_h == 48).mean(),
                        **{f"flag_{f}": g[f].mean() for f in ("excursion", "dew_risk_in_window", "daily_cycle_large",
                                                               "degraded", "toggling")}})
        save_output(pd.DataFrame(cov), OUT / "r1b_real_event_coverage.csv", "diag",
                    "R1b real-event coverage", index=False)
        # Delta_k with the v2 before: real (main) vs duration-matched fake control
        rm = re_[(re_.group == "main") & re_.duration_h.notna()]
        ctrl = []
        for s, g in o.groupby("sublink_id"):
            c = int(s.split("_")[0])
            okp = usable(state[c], V["skip_h"])
            x = tl[s].to_numpy()
            for e in g.itertuples():
                b0 = getattr(e, f"base_{adopted}")
                if not np.isnan(b0):
                    ctrl.append({"D_h": e.D_h, **deltas(x, okp, b0, e.i0 + int(e.D_h * SPH))})
        cv = curve(rm, pd.DataFrame(ctrl), "all main, v2")
        save_output(cv, OUT / "r1b_c2_delta.csv", "diag", "R1b Delta_k with the v2 baseline", index=False)
        plots(tab, o, cv)
        print(dec.round(3).T.to_string(), "\n", tab.round(3).to_string(), "\n", pfs.round(3).to_string(), "\n",
              pd.DataFrame(cov).round(3).to_string(), "\n", cv.round(3).to_string())


if __name__ == "__main__":
    main()
