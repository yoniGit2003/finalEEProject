"""C1: fake events in certified dry weather, where the truth is A = 0, plus the PASS/FAIL criteria (they use C1).

For each sub-link and duration D, up to N random start times t0 are drawn such that every 5 min slot of
[t0, t0 + D) is dry for the link AND inside a network certified dry period. Each baseline variant is computed as if
that span were a rain event (with the real history before it), and e(t) = tl_mean(t) - baseline over the span.
Since nothing should be there, e is the baseline error.

Run: python -m research.baseline.validation.c1_fake_events   (settings: research/config.yaml)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from research.baseline.validation.common import (V, OUT, FIG, CACHE, SPH, VARIANTS, Timer, open_data, sublink_table, groups, load_panel, usable,
                    variants, summary)
from research.common.save import save

C = V["c1"]
PF = V["pass_fail"]


def place_and_compute(sl, tl, state, cert):
    """Draw the fake events and compute every variant. Returns (events table, slot arrays dict)."""
    rng = np.random.default_rng(V["seed"])
    T = tl.index
    rows, parts = [], {k: [] for k in ["ev", "off"] + VARIANTS}
    for s in sl.sublink_id:
        c = int(s.split("_")[0])
        x = tl[s].to_numpy()
        ok = usable(state[c], V["skip_h"])
        good = ok & ~np.isnan(x)
        m0 = float(np.median(x[good])) if good.any() else np.nan            # M0: whole summer (non-causal)
        place = (state[c].to_numpy() == 0) & cert.to_numpy()
        cs = np.concatenate([[0], np.cumsum(place)])
        for D in C["durations_h"]:
            n = int(D * SPH)
            starts = np.nonzero(cs[n:] - cs[:-n] == n)[0]                  # every slot of [i0, i0 + n) placeable
            k = min(C["max_events_per_sublink"], len(starts))
            for i0 in np.sort(rng.choice(starts, k, replace=False)) if k else []:
                base, info = variants(x, ok, m0, i0, i0 + n)
                seg = x[i0:i0 + n]
                ev = len(rows)
                row = {"event": ev, "sublink_id": s, "cml_id": c, "D_h": D, "t0": T[i0], "i0": int(i0),
                       "n_candidates": len(starts), "m0": m0, "m3_method": info["m3_method"],
                       "m3_flag": info["m3_flag"], "m3_before": info["before"], "m3_after": info["after"]}
                for v in VARIANTS:
                    e = seg - base[v]
                    row[f"base_{v}"] = base[v][0]
                    row[f"missing_{v}"] = bool(np.isnan(base[v]).all())
                    row[f"mean_e_{v}"] = float(np.nanmean(e)) if np.isfinite(e).any() else np.nan
                    parts[v].append(e.astype("float32"))
                parts["ev"].append(np.full(n, ev, "int32"))
                parts["off"].append(np.arange(n, dtype="int16"))
                rows.append(row)
    ev = pd.DataFrame(rows)
    slots = {k: np.concatenate(v) for k, v in parts.items()}
    return ev, slots


def table_a(ev, slots, ids):
    """Per variant x D: events, missing share, slot-level error statistics (slots pooled over all fake events)."""
    sel_ev = ev.sublink_id.isin(ids).to_numpy()
    in_slot = sel_ev[slots["ev"]]
    Ds = ev.D_h.to_numpy()[slots["ev"]]
    out = []
    for v in VARIANTS:
        for D in C["durations_h"]:
            m = in_slot & (Ds == D)
            e_ev = ev[sel_ev & (ev.D_h == D)]
            out.append({"variant": v, "D_h": D, "n_events": len(e_ev),
                        "missing_share": e_ev[f"missing_{v}"].mean() if len(e_ev) else np.nan,
                        **summary(slots[v][m])})
    return pd.DataFrame(out)


def by_class(ev, slots, ids, sl, col):
    """M2b slot statistics per D, split by a sub-link column (class or TSL behaviour)."""
    cls = sl.set_index("sublink_id")[col]
    ev_cls = ev.sublink_id.map(cls).to_numpy()
    sel = ev.sublink_id.isin(ids).to_numpy()
    s_cls, s_sel, s_D = ev_cls[slots["ev"]], sel[slots["ev"]], ev.D_h.to_numpy()[slots["ev"]]
    out = []
    for k in sorted(pd.unique(ev_cls[sel])):
        for D in C["durations_h"]:
            m = s_sel & (s_cls == k) & (s_D == D)
            out.append({col: k, "D_h": D, "n_sublinks": int(ev[sel & (ev_cls == k)].sublink_id.nunique()),
                        "n_events": int((sel & (ev_cls == k) & (ev.D_h == D)).sum()), **summary(slots["M2b"][m])})
    return pd.DataFrame(out)


def vs_time(ev, slots, ids):
    """Error against time since t0 (all D pooled), per variant."""
    sel = ev.sublink_id.isin(ids).to_numpy()[slots["ev"]]
    t = slots["off"] / SPH
    edges = C["time_bins_h"]
    out = []
    for v in VARIANTS:
        for a, b in zip(edges[:-1], edges[1:]):
            m = sel & (t >= a) & (t < b)
            out.append({"variant": v, "bin_h": f"{a}-{b}", "from_h": a, "to_h": b, **summary(slots[v][m])})
    return pd.DataFrame(out)


def per_sublink(ev, slots, sl):
    """Per sub-link M2b (and M0 / M1 median |e| for the PASS criterion) at D <= short_max_h and D = long_h."""
    s_id = ev.sublink_id.to_numpy()[slots["ev"]]
    s_D = ev.D_h.to_numpy()[slots["ev"]]
    short, long_ = s_D <= C["short_max_h"], s_D == C["long_h"]
    df = pd.DataFrame({"sid": s_id, "short": short, "long": long_,
                       **{v: slots[v] for v in ("M0", "M1", "M2b")}})
    rows = []
    for s, g in df.groupby("sid", sort=False):
        r = {"sublink_id": s}
        gs, gl = g[g.short], g[g.long]
        for v in ("M2b", "M0", "M1"):
            es = gs[v].dropna().to_numpy()
            r[f"{v}_short_median_abs"] = np.median(np.abs(es)) if len(es) else np.nan
        es, el = gs.M2b.dropna().to_numpy(), gl.M2b.dropna().to_numpy()
        r.update({"n_events_short": int(ev[(ev.sublink_id == s) & (ev.D_h <= C["short_max_h"])].shape[0]),
                  "n_events_24h": int(ev[(ev.sublink_id == s) & (ev.D_h == C["long_h"])].shape[0]),
                  "short_median": np.median(es) if len(es) else np.nan,
                  "short_p10": np.percentile(es, 10) if len(es) else np.nan,
                  "short_p90": np.percentile(es, 90) if len(es) else np.nan,
                  "short_p90_abs": np.percentile(np.abs(es), 90) if len(es) else np.nan,
                  "long_median": np.median(el) if len(el) else np.nan,
                  "long_p10": np.percentile(el, 10) if len(el) else np.nan,
                  "long_p90": np.percentile(el, 90) if len(el) else np.nan,
                  "long_p90_abs": np.percentile(np.abs(el), 90) if len(el) else np.nan})
        rows.append(r)
    out = pd.DataFrame(rows).merge(sl, on="sublink_id", how="left")
    return out


def pass_fail(ps):
    """The four fixed criteria per sub-link; NaN = not evaluable (no fake events of that duration)."""
    def crit(cond, *need):
        ok = np.all([ps[c].notna() for c in need], axis=0)
        return pd.Series(np.where(ok, cond, np.nan), index=ps.index)
    ps = ps.copy()
    ps["c1_median"] = crit(ps.short_median.abs() <= PF["median_abs_db"], "short_median")
    ps["c2_spread"] = crit((ps.short_p10 >= PF["p10_min_db"]) & (ps.short_p90 <= PF["p90_max_db"]), "short_p10")
    ps["c3_growth"] = crit((ps.long_p90_abs - ps.short_p90_abs) <= PF["extra_p90_24h_db"], "long_p90_abs",
                           "short_p90_abs")
    ps["c4_beats"] = crit((ps.M2b_short_median_abs < ps.M0_short_median_abs) &
                          (ps.M2b_short_median_abs < ps.M1_short_median_abs), "M2b_short_median_abs",
                          "M0_short_median_abs", "M1_short_median_abs")
    # information only (the criterion is unchanged): exact ties are common because TL is quantized
    ps["c4_tie"] = (ps.M2b_short_median_abs == ps.M0_short_median_abs) | \
                   (ps.M2b_short_median_abs == ps.M1_short_median_abs)
    cs = ["c1_median", "c2_spread", "c3_growth", "c4_beats"]
    ps["evaluable"] = ps[cs].notna().all(axis=1)
    ps["pass_all"] = ps.evaluable & (ps[cs] == 1).all(axis=1)
    names = {"c1_median": f"|median e| > {PF['median_abs_db']}",
             "c2_spread": f"p10 < {PF['p10_min_db']} or p90 > {PF['p90_max_db']}",
             "c3_growth": f"extra p90|e| at 24 h > {PF['extra_p90_24h_db']}",
             "c4_beats": "M2b not better than M0 and M1"}
    ps["fail_reason"] = [
        "; ".join([names[c] for c in cs if r[c] == 0] + [f"not evaluable ({c})" for c in cs if pd.isna(r[c])])
        for _, r in ps.iterrows()]
    return ps


def pf_summary(pf, by=None):
    cs = ["c1_median", "c2_spread", "c3_growth", "c4_beats"]
    g = pf.groupby(by) if by else [("all", pf)]
    rows = []
    for k, x in g:
        r = {"group": k if not isinstance(k, tuple) else " / ".join(map(str, k)), "n_sublinks": len(x),
             "n_evaluable_all": int(x.evaluable.sum())}
        for c in cs:
            r[f"pass_{c}"] = x[c].mean() if x[c].notna().any() else np.nan    # share among evaluable
        r["pass_all"] = x[x.evaluable].pass_all.mean() if x.evaluable.any() else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def plot_vs_time(vt):
    fig, axes = plt.subplots(1, len(VARIANTS), figsize=(16, 3.8), sharey=True)
    for ax, v in zip(axes, VARIANTS):
        g = vt[vt.variant == v]
        x = (g.from_h + g.to_h) / 2
        ax.fill_between(x, g.p10, g.p90, alpha=0.25, color="C0", label="p10-p90")
        ax.plot(x, g["median"], "o-", color="C0", label="median")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xscale("log")
        ax.set_xticks([0.5, 2, 4.5, 9, 18, 32], ["0-1", "1-3", "3-6", "6-12", "12-24", "24-40"], fontsize=8)
        ax.set_title(v)
        ax.set_xlabel("time since t0 [h]")
    axes[0].set_ylabel("e = tl_mean - baseline [dB]")
    axes[0].legend(fontsize=8)
    fig.suptitle("C1: baseline error in dry weather against time since the (fake) event start, main sub-links")
    fig.tight_layout()
    save(fig, FIG / "c1c_error_vs_time.png", dpi=110)
    plt.close(fig)


def hour_humidity(ev, ids, met):
    """Per-event mean e (D <= short_max_h) against the hour of day and the humidity at t0."""
    e = ev[ev.sublink_id.isin(ids) & (ev.D_h <= C["short_max_h"])].copy()
    e["hour_utc"] = pd.to_datetime(e.t0).dt.hour
    e["rh"] = met.humidity.reindex(pd.to_datetime(e.t0)).to_numpy()
    e["rh_bin"] = pd.cut(e.rh, [0, 60, 70, 80, 90, 101], right=False).astype(str)
    rows = []
    for key, col in (("hour_utc", "hour_utc"), ("rh_bin", "rh_bin")):
        for k, g in e.groupby(col):
            r = {"by": key, "bin": k, "n_events": len(g)}
            for v in ("M0", "M1", "M2b"):
                x = g[f"mean_e_{v}"].dropna()
                r[f"{v}_median"] = x.median() if len(x) else np.nan
                if v == "M2b":
                    r["M2b_p10"], r["M2b_p90"] = (x.quantile(0.1), x.quantile(0.9)) if len(x) else (np.nan, np.nan)
            rows.append(r)
    t = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4))
    for ax, key in zip(axes, ("hour_utc", "rh_bin")):
        g = t[t.by == key].reset_index(drop=True)
        x = np.arange(len(g))
        ax.fill_between(x, g.M2b_p10, g.M2b_p90, alpha=0.2, color="C0", label="M2b p10-p90")
        for v, c in (("M2b", "C0"), ("M1", "C1"), ("M0", "C2")):
            ax.plot(x, g[f"{v}_median"], "o-", color=c, ms=3, label=f"{v} median")
        ax.set_xticks(x, g["bin"], rotation=0 if key == "hour_utc" else 20, fontsize=8)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xlabel("hour of t0 (UTC; local = UTC + 2)" if key == "hour_utc" else "relative humidity at t0 [%]")
        ax.set_ylabel("mean e of the fake event [dB]")
    axes[0].legend(fontsize=8)
    fig.suptitle(f"C1: fake-event error (D <= {C['short_max_h']} h) against hour of day and humidity, main sub-links")
    fig.tight_layout()
    save(fig, FIG / "c1e_error_vs_hour_humidity.png", dpi=110)
    plt.close(fig)
    return t


def main():
    with Timer("C1"):
        d = open_data()
        sl = sublink_table(d)
        tl, state, cert, met = load_panel(d, sl)
        ev, slots = place_and_compute(sl, tl, state, cert)
        save(ev, OUT / "c1_fake_events.csv.gz", index=False)
        np.savez_compressed(CACHE / "c1_slots.npz", **slots)
        G = groups(sl)
        # counts: how many fake events could be placed
        cnt = (ev.groupby(["sublink_id", "D_h"]).size().unstack(fill_value=0)
               .reindex(sl.sublink_id, fill_value=0))
        crows = []
        for gname, ids in G.items():
            for D in C["durations_h"]:
                n = cnt.loc[ids, D]
                crows.append({"group": gname, "D_h": D, "n_sublinks": len(ids), "events": int(n.sum()),
                              "sublinks_full": int((n == C["max_events_per_sublink"]).sum()),
                              "sublinks_fewer": int(((n > 0) & (n < C["max_events_per_sublink"])).sum()),
                              "sublinks_none": int((n == 0).sum())})
        save(pd.DataFrame(crows), OUT / "c1_event_counts.csv", index=False)
        for gname, ids in G.items():
            save(table_a(ev, slots, ids), OUT / f"c1a_variants_by_D_{gname}.csv", index=False)
        tab_cls = by_class(ev, slots, G["main"], sl.assign(cls=sl.length_class + " / " + sl.frequency_band), "cls")
        save(tab_cls, OUT / "c1b_M2b_by_class.csv", index=False)
        save(by_class(ev, slots, G["main"], sl, "tsl_behaviour"), OUT / "c1b_M2b_by_tsl.csv", index=False)
        vt = vs_time(ev, slots, G["main"])
        save(vt, OUT / "c1c_error_vs_time.csv", index=False)
        plot_vs_time(vt)
        save(hour_humidity(ev, G["main"], met), OUT / "c1e_error_vs_hour_humidity.csv", index=False)
        ps = per_sublink(ev, slots, sl)
        save(ps, OUT / "c1d_M2b_per_sublink.csv", index=False)
        pf = pass_fail(ps)
        save(pf, OUT / "c1_pass_fail_per_sublink.csv", index=False)
        main_pf = pf[pf.group == "main"]
        summ = pd.concat([pf_summary(main_pf).assign(split="all main"),
                          pf_summary(main_pf, ["length_class", "frequency_band"]).assign(split="length / frequency"),
                          pf_summary(main_pf, "tsl_behaviour").assign(split="TSL class"),
                          pf_summary(pf[pf.is_degraded]).assign(split="degraded (separate)"),
                          pf_summary(pf[pf.is_toggling]).assign(split="toggling (separate)")])
        save(summ, OUT / "c1_pass_fail_summary.csv", index=False)
        print(summ.round(3).to_string())


if __name__ == "__main__":
    main()
