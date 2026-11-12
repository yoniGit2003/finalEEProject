"""Step A4 of the k-surface branch: where do the data visit the (R, W) plane? The learning domain for k(R, W).

    python -m research.waa_ksurface.a4_coverage

On the TRAIN links below 1 km, in the 14 events with a gauge ratio, every valid minute gets

    W_hat = (A - c) - L * P_hat(R)          R = the radar value of the minute's 5 min bin / g_e

P_hat is the measured P of step E4 inside its valid range and a LINEAR CONTINUATION above it. W_hat is used for this map
only; it is never the target of a fit. The minutes are counted in the cells of the (R, W) plane, wetting (W_hat below the
measured F) and drying (above) apart. A cell belongs to the learning domain if it has enough minutes from enough events
under both continuations of P. Every definition and both floors: decisions.md, step A4. Settings: config section `a4`.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from research.common.plotting import INK, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import data, fp, ltest
from research.waa_ksurface import paths as P

LT = ltest.config()
C = data.branch_config()["a4"]
EDGES = np.asarray(LT["class_edges_mmh"], float)
R_NAMES = [f"< {EDGES[0]:g}"] + ltest.class_names(LT)          # column 0 = not wet by the data layer
REGIMES = ["wetting", "drying"]


def minutes(D, le, keep):
    """One row per valid minute of the link-events `keep`: link-event, event, band, L, R, A - c."""
    n = D["le_n"]
    rep = np.repeat(np.arange(len(le)), n)
    r1 = np.repeat(D["r5"].astype(float), 5)
    ok = (D["w1"] > 0) & np.isfinite(r1) & keep[rep]
    t = pd.DataFrame({"lev": rep[ok], "R": r1[ok], "y": D["y1"][ok] - D["le_c"][rep[ok]]})
    t["ev"] = D["le_event"][t.lev.to_numpy()]
    t["band"] = le.band_name.to_numpy()[t.lev.to_numpy()]
    t["L"] = D["le_L"][t.lev.to_numpy()]
    return t


def one_minute_noise(D, le, keep):
    """sigma_1: the sd of a single 1 min value of A - c on dry minutes = the median over link-events of the sd inside the
    offset window (at least 30 valid minutes)."""
    sds = []
    for e in np.nonzero(keep)[0]:
        o, n = D["le_boff1"][e], D["le_n"][e]
        m = D["pre1"][o:o + n] & (D["w1"][o:o + n] > 0)
        if m.sum() >= 30:
            sds.append(np.std(D["y1"][o:o + n][m], ddof=1))
    return float(np.median(sds)), len(sds)


def place(t, M):
    """Add W_hat under both continuations of P, the regime, and the cell indices (R column, W bin)."""
    t = t.copy()
    for band, m in M.items():
        b = (t.band == band).to_numpy()
        r = t.R.to_numpy()[b]
        t.loc[b, "F"] = m.F(r)
        t.loc[b, "W_L1"] = t.y.to_numpy()[b] - t.L.to_numpy()[b] * m.P(r)             # E4 rule: last-segment slope
        t.loc[b, "W_L2"] = t.y.to_numpy()[b] - t.L.to_numpy()[b] * m.P_prop(r)        # proportional to R above the valid range
        t.loc[b, "above_F"] = r > m.F_max
        t.loc[b, "above_P"] = r > m.P_max
    t["rcol"] = np.where(t.R < EDGES[0], 0, np.searchsorted(EDGES, t.R.to_numpy(), side="right"))   # 1..7 = the classes, 8 = R >= 20
    return t


def count(t, wcol, n_w):
    """Minutes and events per (band, regime, R column, W bin) for one continuation of P."""
    w = t[wcol].to_numpy()
    wbin = np.floor((w - C["w_min_db"]) / C["w_step_db"]).astype(int)
    inside = (wbin >= 0) & (wbin < n_w) & (t.rcol.to_numpy() < len(R_NAMES))
    g = t[inside].assign(wbin=wbin[inside], regime=np.where(w[inside] < t.F.to_numpy()[inside], "wetting", "drying"))
    per_band = g.groupby(["band", "regime", "rcol", "wbin"]).agg(minutes=("ev", "size"), events=("ev", "nunique")).reset_index()
    pooled = g.groupby(["regime", "rcol", "wbin"]).agg(minutes=("ev", "size"), events=("ev", "nunique")).reset_index().assign(band="pooled")
    return pd.concat([per_band, pooled], ignore_index=True), int((~inside).sum())


def domain_table(t, n_min, n_w):
    """Per band (and pooled), regime and cell: the counts under both continuations and whether the cell is in the domain."""
    a, out1 = count(t, "W_L1", n_w)
    b, out2 = count(t, "W_L2", n_w)
    key = ["band", "regime", "rcol", "wbin"]
    d = a.merge(b, on=key, how="outer", suffixes=("_L1", "_L2")).fillna(0)
    for c in ("minutes_L1", "events_L1", "minutes_L2", "events_L2"):
        d[c] = d[c].astype(int)
    d["in_domain"] = (d[["minutes_L1", "minutes_L2"]].min(axis=1) >= n_min) & (d[["events_L1", "events_L2"]].min(axis=1) >= C["min_events"])
    d["R_class"] = [R_NAMES[i] for i in d.rcol]
    d["W_lo"] = C["w_min_db"] + d.wbin * C["w_step_db"]
    d["W_hi"] = d.W_lo + C["w_step_db"]
    return d.sort_values(key).reset_index(drop=True), out1, out2


def summary(d, M):
    """Per band and regime: the number of domain cells and the R and W range they span."""
    f_max = min(m.F_max for m in M.values())
    rows = []
    for (band, regime), g in d[d.in_domain].groupby(["band", "regime"]):
        top = int(g.rcol.max())
        rows.append({"band": band, "regime": regime, "domain_cells": len(g), "minutes_in_domain_L2": int(g.minutes_L2.sum()),
                     "R_classes": f"{R_NAMES[int(g.rcol.min())]} .. {R_NAMES[top]}", "R_edge_mmh": float(EDGES[top]) if top >= 1 else float(EDGES[0]),
                     "W_lo_db": float(g.W_lo.min()), "W_hi_db": float(g.W_hi.max()),
                     "cells_with_R_above_F_range": int((np.r_[0.0, EDGES][g.rcol.to_numpy()] >= f_max).sum())})
    return pd.DataFrame(rows)


def fig_maps(d, M, n_min, n_w):
    f_max = min(m.F_max for m in M.values())
    bands = list(M) + ["pooled"]
    fig, axes = plt.subplots(len(bands), 2, figsize=(11, 3.4 * len(bands)), sharex=True, sharey=True)
    for i, band in enumerate(bands):
        for j, regime in enumerate(REGIMES):
            ax = axes[i, j]
            g = d[(d.band == band) & (d.regime == regime)]
            Z = np.zeros((n_w, len(R_NAMES)))
            Z[g.wbin.to_numpy(), g.rcol.to_numpy()] = g.minutes_L2.to_numpy()
            im = ax.pcolormesh(np.arange(len(R_NAMES) + 1), C["w_min_db"] + C["w_step_db"] * np.arange(n_w + 1), np.log10(np.maximum(Z, 1.0)),
                               cmap="Greys", vmin=0, vmax=5.5)
            for r in g[g.in_domain].itertuples():                       # the domain: outlined cells
                ax.add_patch(Rectangle((r.rcol, r.W_lo), 1, C["w_step_db"], fill=False, ec=INK, lw=0.9))
            first_above = int(np.searchsorted(EDGES, f_max, side="right"))          # the class that holds F's last measured R
            ax.add_patch(Rectangle((first_above + 1, C["w_min_db"]), len(R_NAMES) - first_above - 1, n_w * C["w_step_db"], fill=False, hatch="//", ec="0.45", lw=0))
            ax.axvline(first_above + (f_max - EDGES[first_above - 1]) / (EDGES[first_above] - EDGES[first_above - 1]), color=INK, ls="--", lw=1)
            ax.set_title(f"{band}: {regime} (W_hat {'below' if regime == 'wetting' else 'at or above'} the measured F)", fontsize=9)
            ax.set_xticks(np.arange(len(R_NAMES)) + 0.5)
            ax.set_xticklabels(R_NAMES, fontsize=7)
            if i == len(bands) - 1:
                ax.set_xlabel("rain R (mm/h), radar / g_e")
            if j == 0:
                ax.set_ylabel("W_hat (dB)")
            style_axes(ax)
        fig.colorbar(im, ax=axes[i, :], label="log10(minutes)", shrink=0.8)
    fig.suptitle("A4 coverage map, TRAIN links below 1 km, events with a gauge ratio. MAP ONLY, never used for fitting.\n"
                 "W_hat = A - c - L*P_hat(R); above P's valid range P_hat is a LINEAR CONTINUATION (not measured).\n"
                 f"Outlined cells = the learning domain (>= {n_min} minutes from >= {C['min_events']} events, under both continuations). "
                 "Dashed line and hatching: R above F's measured range.", fontsize=9)
    return fig


def main():
    d = data.open_release()
    D, le, _ = data.gauge_ratio_train(d)
    short = (le.L_km < LT["short_km"]).to_numpy()
    M = fp.measured()
    sigma_raw, n_le = one_minute_noise(D, le, short)
    sigma_1 = max(sigma_raw, C["quant_step_db"] / np.sqrt(12.0))
    n_min = int(5 * np.ceil((sigma_1 / C["delta_db"]) ** 2))
    n_w = int(round((C["w_max_db"] - C["w_min_db"]) / C["w_step_db"]))
    print(f"sigma_1 = {sigma_raw:.4f} dB ({n_le} link-events) -> minutes floor {n_min}; events floor {C['min_events']}", flush=True)

    t = place(minutes(D, le, short), M)
    dom, out1, out2 = domain_table(t, n_min, n_w)
    summ = summary(dom, M)
    # the gauge-input set (short links with a City gauge within 1 km), same floors
    tg = place(minutes(D, le, short & (le.gauge_id != "").to_numpy()), M)
    dom_g = domain_table(tg, n_min, n_w)[0]
    summ_g = summary(dom_g, M)
    # shares of the rainy minutes above the measured ranges
    rainy = t[t.R >= EDGES[0]]
    shares = pd.DataFrame([{"band": b, "rainy_minutes": len(g), "share_above_F_range": float(g.above_F.mean()), "share_above_P_range": float(g.above_P.mean()),
                            "F_max_mmh": M[b].F_max, "P_max_mmh": M[b].P_max} for b, g in rainy.groupby("band")]
                          + [{"band": "pooled", "rainy_minutes": len(rainy), "share_above_F_range": float(rainy.above_F.mean()),
                              "share_above_P_range": float(rainy.above_P.mean()), "F_max_mmh": np.nan, "P_max_mmh": np.nan}])
    # the g_e = 1 events: counted, not mapped
    _, D1, le_all, _ = data.branch_train(d)
    ge1 = ~le_all.event_id.isin(set(le.event_id)).to_numpy() & (le_all.L_km < LT["short_km"]).to_numpy()
    t1 = minutes(D1, le_all, ge1)
    info = pd.DataFrame([{"set": "events with a gauge ratio (the learning set)", "events": le.event_id.nunique(), "short_links": le[short].cml_id.nunique(),
                          "link_events": int(short.sum()), "valid_minutes": len(t), "rainy_minutes": len(rainy)},
                         {"set": "g_e = 1 events (not used)", "events": le_all[ge1].event_id.nunique(), "short_links": le_all[ge1].cml_id.nunique(),
                          "link_events": int(ge1.sum()), "valid_minutes": len(t1), "rainy_minutes": int((t1.R >= EDGES[0]).sum())},
                         {"set": "gauge-input set inside the learning set", "events": tg.ev.nunique(), "short_links": le[short & (le.gauge_id != "").to_numpy()].cml_id.nunique(),
                          "link_events": int((short & (le.gauge_id != "").to_numpy()).sum()), "valid_minutes": len(tg), "rainy_minutes": int((tg.R >= EDGES[0]).sum())}])
    pooled = dom[(dom.band == "pooled") & dom.in_domain]
    edges = pd.DataFrame([{"quantity": "sigma_1_db", "value": sigma_raw}, {"quantity": "sigma_1_used_db", "value": sigma_1},
                          {"quantity": "sigma_1_link_events", "value": n_le}, {"quantity": "delta_db", "value": C["delta_db"]},
                          {"quantity": "min_minutes", "value": n_min}, {"quantity": "min_events", "value": C["min_events"]},
                          {"quantity": "domain_R_edge_mmh", "value": float(EDGES[int(pooled.rcol.max())]) if len(pooled) else np.nan},
                          {"quantity": "domain_W_edge_db", "value": float(pooled.W_hi.max()) if len(pooled) else np.nan},
                          {"quantity": "domain_cells_wetting", "value": int((pooled.regime == "wetting").sum())},
                          {"quantity": "domain_cells_drying", "value": int((pooled.regime == "drying").sum())},
                          {"quantity": "minutes_outside_plane_L1", "value": out1}, {"quantity": "minutes_outside_plane_L2", "value": out2}])
    for tab, f in ((dom, P.A4_DOMAIN), (summ, P.A4_SUMMARY), (dom_g, P.A4_DOMAIN_GAUGE), (shares.round(5), P.A4_SHARES), (info, P.A4_INFO), (edges, P.A4_EDGES)):
        save(tab, f, index=False)
    fig = fig_maps(dom, M, n_min, n_w)
    save(fig, P.A4_FIG, dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(summ.to_string(), flush=True)
    print(shares.round(4).to_string(), flush=True)
    print(info.to_string(), flush=True)

    def grid(band, regime, table):
        g = table[(table.band == band) & (table.regime == regime) & table.in_domain]
        piv = g.pivot(index="W_lo", columns="rcol", values="minutes_L2").sort_index(ascending=False)
        piv.columns = [R_NAMES[c] for c in piv.columns]
        return md_table(piv.fillna(0).astype(int).rename_axis("W bin from (dB)").reset_index().round(2))

    a = ["# A4: coverage of the (R, W) plane and the learning domain\n",
         f"Written by `research/waa_ksurface/a4_coverage.py`. TRAIN, links below {LT['short_km']:g} km (28 and 38 GHz), the 14 events with a gauge ratio. "
         "**Ŵ = (A − c) − L·P̂(R) is used for this map only, never as the target of a fit.** Above P's valid range P̂ is a linear continuation, not a measurement. "
         "VAL and TEST were not opened.\n",
         "## The data\n", md_table(info),
         "\n## The floors\n",
         f"- Events: at least **{C['min_events']}** (with 14 events a cell fed by 5 events is empty in 0.2 % of the bootstrap draws, by 4 events in 0.9 %).",
         f"- Minutes: σ_1 = {sigma_raw:.3f} dB (sd of one 1 min value on dry minutes, median of {n_le} link-events; lower limit 0.087 dB). "
         f"n_min = 5·ceil((σ_1 / {C['delta_db']})²) = **{n_min} minutes**.",
         "- A cell is in the domain only if it passes both floors under both continuations of P (last-segment slope, and proportional to R).\n",
         "## The domain\n", md_table(summ.round(3)),
         f"\nPooled domain: R up to the class edge **{edges.set_index('quantity').value['domain_R_edge_mmh']:g} mm/h**, W up to **{edges.set_index('quantity').value['domain_W_edge_db']:.1f} dB**. "
         f"Frozen table: `{rel(P.A4_DOMAIN)}` (column `in_domain`). Figure: `{rel(P.A4_FIG)}`.\n",
         "### Pooled, wetting: minutes per domain cell (continuation L2)\n", grid("pooled", "wetting", dom),
         "\n### Pooled, drying: minutes per domain cell (continuation L2)\n", grid("pooled", "drying", dom),
         "\n## Rain above the measured ranges\n",
         "Share of the rainy minutes (R >= 0.1 mm/h) of the short links whose R is above the valid range of F, and of P:\n", md_table(shares.round(3)),
         "\n## The gauge-input set (short links with a City gauge within 1 km), same floors\n", md_table(summ_g.round(3)) if len(summ_g) else "no domain cell",
         f"\n(`{rel(P.A4_DOMAIN_GAUGE)}`, `{rel(P.A4_SUMMARY)}`, `{rel(P.A4_SHARES)}`, `{rel(P.A4_EDGES)}`)\n"]
    save("\n".join(a) + "\n", P.A4_REPORT)
    print(f"report: {rel(P.A4_REPORT)}", flush=True)


if __name__ == "__main__":
    main()
