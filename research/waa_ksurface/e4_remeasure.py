"""Step E4 of the k-surface branch: F(R) (steady WAA) and P(R) (path loss per km) re-measured with the adopted method.

    python -m research.waa_ksurface.e4_remeasure

Method (candidate (ii), adopted in step T3): the bins are selected by Z = the mean radar rain over all links of the event at
that time (30 min within a factor 1.5), so that the selection does not depend on a single link's own rain error. Per band
and rain class the line (A - c) = I + s L over the links gives I = F and s = P; the cell sits at the mean of the links' own
rain. Validity is per quantity (Deviation 4, 2026-11-12): F is taken from the cells whose fake-world INTERCEPT error is
inside tolerance at both ends of the measured rain error, P from the cells whose SLOPE error is inside tolerance.
Also: the spread of the per-link residuals across storms, and the test of whether a per-link factor is needed (criterion
written in decisions.md before this run). TRAIN only, events with a gauge ratio.
"""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import data, fakeworld as fw, ltest
from research.waa_ksurface import paths as P
from research.waa_ksurface.e2_candidates import sizes_and_rule

LT, FK = ltest.config(), fw.config()
C = data.branch_config()["e4"]
MARK = {"28 GHz": ("o", "-"), "38 GHz": ("s", "--")}


def measure(D, le, rule, n_short):
    """The cells of candidate (ii) on real data with the event bootstrap, and the link-event values behind them."""
    z = fw.area_mean(D["r5"].astype(float), le)
    t = ltest.steady_values(D, le, LT, rule, sel=z, r_stat="mean")
    mult = ltest.multiplicities(int(D["le_event"].max()) + 1, LT["boot_reps"], LT["seed"])
    cells, draws = ltest.cell_fits(t, LT, n_short, mult, LT["interval"], r_stat="mean")
    return t, cells, draws


def validity(l_med):
    """Deviation 4: every judged cell of candidate (ii) with its fake-world errors (read from the committed T3 table, no new
    fit) and its validity per quantity. F-valid: |intercept error| <= tolerance at both ends. P-valid: |slope error| <=
    tolerance / L_med at both ends. A cell that was not judged at both ends is valid for nothing."""
    v = pd.read_csv(P.T3_CELLS)
    v = v[v.candidate == "ii"][["band", "class", "cls", "I_err_lo", "I_err_hi", "s_err_lo", "s_err_hi", "valid"]].reset_index(drop=True)
    tol_i, tol_s = FK["pass"]["intercept_db"], FK["pass"]["intercept_db"] / l_med
    both = v[["I_err_lo", "I_err_hi", "s_err_lo", "s_err_hi"]].notna().all(axis=1)
    v["F_valid"] = both & (v.I_err_lo.abs() <= tol_i) & (v.I_err_hi.abs() <= tol_i)
    v["P_valid"] = both & (v.s_err_lo.abs() <= tol_s) & (v.s_err_hi.abs() <= tol_s)
    assert ((v.F_valid & v.P_valid) == v.valid.astype(bool)).all(), "the split must reproduce the T3 verdict where both hold"
    return v.rename(columns={"valid": "valid_T3_both"}), tol_i, tol_s


def with_error_range(cells, val):
    """Keep the cells that are F-valid or P-valid and add the error-range interval from their fake-world errors at both ends."""
    v = val[val.F_valid | val.P_valid].set_index(["band", "cls"])
    rows = []
    for r in cells.itertuples():
        key = (r.band, r.cls)
        if key not in v.index:
            continue
        e = v.loc[key]
        row = r._asdict()
        row.pop("Index")
        row["class"] = row.pop("_2")
        for name, a, b in (("I", e.I_err_lo, e.I_err_hi), ("s", e.s_err_lo, e.s_err_hi)):
            row[f"{name}_fake_err_lo"], row[f"{name}_fake_err_hi"] = a, b
            row[f"{name}_range_lo"], row[f"{name}_range_hi"] = row[name] - max(a, b, 0.0), row[name] - min(a, b, 0.0)
            row[f"{name}_all_lo"] = row[f"{name}_lo"] - max(a, b, 0.0)
            row[f"{name}_all_hi"] = row[f"{name}_hi"] - min(a, b, 0.0)
        row["F_valid"], row["P_valid"] = bool(e.F_valid), bool(e.P_valid)
        rows.append(row)
    return pd.DataFrame(rows)


def curves(cells, draws):
    """F through the F-valid cells, P through the P-valid cells (each with its own measured range)."""
    grid = np.r_[0.0, np.geomspace(0.01, 20.0, 300)]
    out = []
    for band, gb in cells[cells.eligible].groupby("band"):
        col = {}
        for kind, name, k in (("F", "I", 0), ("P", "s", 1)):
            g = gb[gb[f"{kind}_valid"]].sort_values("R_cell")
            if len(g) < 2:
                break
            rc = g.R_cell.to_numpy()
            w = 1.0 / np.maximum(g[f"{name}_se"].to_numpy(), 1e-6) ** 2
            col[kind] = ltest.monotone_curve(rc, g[name].to_numpy(), w, kind)[0](grid)
            Dr = np.array([draws[(band, c)][:, k] for c in g["class"]])
            good = np.isfinite(Dr).all(axis=0)
            per = np.array([ltest.monotone_curve(rc, Dr[:, j], w, kind)[0](grid) for j in np.nonzero(good)[0]])
            col[f"{kind}_lo"], col[f"{kind}_hi"] = np.percentile(per, LT["interval"], axis=0)
            col[f"{kind}_all_lo"] = np.maximum(ltest.monotone_curve(rc, g[f"{name}_all_lo"].to_numpy(), w, kind)[0](grid), 0.0)
            col[f"{kind}_all_hi"] = ltest.monotone_curve(rc, g[f"{name}_all_hi"].to_numpy(), w, kind)[0](grid)
            col[f"{kind}_measured"] = (grid >= rc.min()) & (grid <= rc.max())
        else:
            out.append(pd.DataFrame({"band": band, "R": grid, **col}))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def link_test(t, cells):
    """The per-link residual spread across storms on the short links, and the criterion of decisions.md (step E4)."""
    cells = cells[cells.F_valid & cells.P_valid]                 # the 8 cells of the criterion as it was written (Deviation 4: kept)
    r = ltest.residuals(t, cells[cells.eligible].assign(cls=lambda x: x.cls))
    ok = set(zip(cells[cells.eligible].band, cells[cells.eligible].cls))
    in_valid = np.array([k in ok for k in zip(r.band_name, r.cls)])
    r = r[in_valid & (r.L_km < LT["short_km"]).to_numpy() & r.resid.notna().to_numpy()]
    per = r.groupby(["band_name", "cml_id", "event_id"]).resid.mean().reset_index()
    between, within, n_links, n_values = ltest.variance_components(per.resid.to_numpy(), per.cml_id.to_numpy())
    rng = np.random.default_rng(C["seed"])
    null = []
    for _ in range(C["permutations"]):
        lab = per.cml_id.to_numpy().copy()
        for b in per.band_name.unique():                        # shuffle the link labels inside each band
            m = (per.band_name == b).to_numpy()
            lab[m] = rng.permutation(lab[m])
        null.append(ltest.variance_components(per.resid.to_numpy(), lab)[0])
    p = float((np.sum(np.array(null) >= between) + 1) / (len(null) + 1))
    # report only: does a link's residual grow with the level of F (a factor) or stay (an offset)?
    lowhigh = r.assign(half=np.where(r.cls <= r.cls.min() + 1, "low", "high")).groupby(["cml_id", "half"]).resid.mean().unstack()
    lowhigh = lowhigh.dropna()
    corr = float(np.corrcoef(lowhigh.low, lowhigh.high)[0, 1]) if len(lowhigh) >= 5 else np.nan
    slope = float(np.polyfit(lowhigh.low, lowhigh.high, 1)[0]) if len(lowhigh) >= 5 else np.nan
    links = per.groupby(["band_name", "cml_id"]).agg(mean_resid=("resid", "mean"), storms=("resid", "size"), sd_over_storms=("resid", "std")).reset_index()
    links = links.merge(r.drop_duplicates("cml_id")[["cml_id", "L_km"]], on="cml_id")
    needed = bool(between >= C["link_sd_threshold_db"] and p < C["p_threshold"])
    test = {"sigma_link_db": between, "storm_sd_db": within, "short_links": n_links, "link_storm_values": n_values, "permutation_p": p,
            "null_median_db": float(np.median(null)), "null_95_db": float(np.percentile(null, 95)), "threshold_db": C["link_sd_threshold_db"],
            "per_link_factor_needed": needed, "links_low_high": len(lowhigh), "corr_low_high": corr, "slope_high_on_low": slope}
    return links, test


def versus_a3(cells):
    """The new F against the superseded A3a curve at the same R (information only)."""
    old = pd.read_csv(P.A3_CURVES)
    rows = []
    for r, cls_name in zip(cells.itertuples(), cells["class"]):
        g = old[old.band == r.band].sort_values("R")
        f_old = float(np.interp(r.R_cell, g.R, g.F))
        rows.append({"band": r.band, "class": cls_name, "R_cell": r.R_cell, "F_new": r.I, "F_a3_superseded": f_old, "shift_db": r.I - f_old})
    return pd.DataFrame(rows)


def fig_curves(cur, cells):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, kind, name, unit in ((axes[0], "F", "I", "dB"), (axes[1], "P", "s", "dB/km")):
        for band in cells.band.unique():
            mk, ls = MARK[band]
            g = cur[(cur.band == band) & (cur.R <= 6)]
            e = cells[(cells.band == band) & cells[f"{kind}_valid"]]
            ax.plot(g.R, g[kind], ls, color=INK, lw=1.4, label=f"{band}: curve")
            ax.fill_between(g.R, g[f"{kind}_all_lo"], g[f"{kind}_all_hi"], color=INK2, alpha=0.15, hatch="//" if band == "28 GHz" else "\\\\", lw=0)
            ax.errorbar(e.R_cell, e[name], yerr=[e[name] - e[f"{name}_lo"], e[f"{name}_hi"] - e[name]], fmt=mk, color=INK, mfc="white", capsize=3, lw=0.8,
                        label=f"{band}: cells (bootstrap 90 %)")
            ax.axvline(e.R_cell.min(), color=INK2, lw=0.6, ls=":")
            ax.axvline(e.R_cell.max(), color=INK2, lw=0.6, ls=":")
        ax.set_xscale("symlog", linthresh=0.1)
        ax.set_xlabel("rain R (mm/h)")
        ax.set_ylabel(f"{kind}(R) ({unit})")
        ax.set_title(("F: steady WAA" if kind == "F" else "P: path loss per km") + " (hatched: bootstrap + error range)\noutside the dotted lines: not measured", fontsize=9)
        ax.legend(fontsize=7)
        style_axes(ax)
    fig.tight_layout()
    return fig


def fig_links(links, test):
    fig, ax = plt.subplots(figsize=(7, 4))
    for band, (mk, _) in MARK.items():
        m = links[links.band_name == band]
        ax.errorbar(m.L_km, m.mean_resid, yerr=m.sd_over_storms / np.sqrt(m.storms.clip(lower=1)), fmt=mk, color=INK, mfc="white", capsize=2, lw=0.7, label=band)
    ax.axhline(0, color=INK2, lw=0.6)
    for lim in (C["link_sd_threshold_db"], -C["link_sd_threshold_db"]):
        ax.axhline(lim, color=INK2, ls=":", lw=1)
    ax.set_xlabel("link length (km)")
    ax.set_ylabel("mean residual of the link over its storms (dB)")
    ax.set_title(f"short links: residual from the cell lines (bar = standard error over storms)\nbetween-link sd {test['sigma_link_db']:.2f} dB, dotted = ±{C['link_sd_threshold_db']} dB", fontsize=9)
    ax.legend(fontsize=7)
    style_axes(ax)
    fig.tight_layout()
    return fig


def main():
    D, le, _ = data.gauge_ratio_train()
    _, rule, rule_name, n_short, l_med = sizes_and_rule()
    t, all_cells, draws = measure(D, le, rule, n_short)
    val, tol_i, tol_s = validity(l_med)
    save(val.round(5), P.E4_VALIDITY, index=False)
    print(val.round(4).to_string(), flush=True)
    cells = with_error_range(all_cells, val)
    dropped = cells[~cells.eligible]
    cells = cells[cells.eligible].reset_index(drop=True)
    cur = curves(cells, draws)
    links, test = link_test(t, cells)
    vs = versus_a3(cells)
    print(cells[["band", "class", "R_cell", "I", "I_lo", "I_hi", "I_range_lo", "I_range_hi", "s", "s_lo", "s_hi", "s_range_lo", "s_range_hi", "links_short", "events"]].round(3).to_string(), flush=True)
    print(test, flush=True)
    print(vs.round(3).to_string(), flush=True)
    for tab, f in ((cells, P.E4_CELLS), (cur, P.E4_CURVES), (links, P.E4_LINKS), (pd.DataFrame([test]), P.E4_LINK_TEST), (vs, P.E4_VS_A3)):
        save(tab.round(5), f, index=False)
    for fig, f in ((fig_curves(cur, cells), P.E4_FIG), (fig_links(links, test), P.E4_FIG_LINKS)):
        save(fig, f, dpi=120)
        plt.close(fig)

    show = cells.assign(F=[f"{r.I:.3f} [{r.I_lo:.3f}, {r.I_hi:.3f}] [{r.I_range_lo:.3f}, {r.I_range_hi:.3f}]" if r.F_valid else "not F-valid" for r in cells.itertuples()],
                        Pv=[f"{r.s:.3f} [{r.s_lo:.3f}, {r.s_hi:.3f}] [{r.s_range_lo:.3f}, {r.s_range_hi:.3f}]" if r.P_valid else "not P-valid" for r in cells.itertuples()])
    nF, nP = int(cells.F_valid.sum()), int(cells.P_valid.sum())
    a = ["# E4: F(R) and P(R) re-measured (candidate (ii), instrument binning)\n",
         f"Written by `research/waa_ksurface/e4_remeasure.py`. TRAIN, {len(le)} link-events of {le.cml_id.nunique()} links in the {le.event_id.nunique()} events with a gauge ratio "
         f"(28 and 38 GHz). Bins selected by the event's area-mean radar rain, rule {rule_name} ({rule['t_pre_min']} min within a factor {rule['q']:g}); the line over L per band and "
         "class; the cell at the mean of the links' own rain. **Validity is per quantity (Deviation 4, 2026-11-12):** F from the "
         f"{nF} F-valid cells, P from the {nP} P-valid cells. VAL and TEST were not opened.\n",
         "## Validity per quantity (the committed T3 fake results, re-read; no new fit)\n",
         f"Fake-world error of candidate (ii) at the two ends of the measured rain error. F-valid: |intercept error| <= {tol_i} dB at both ends. "
         f"P-valid: |slope error| <= {tol_s:.4f} dB/km at both ends. (`{rel(P.E4_VALIDITY)}`)\n",
         md_table(val.drop(columns=["cls"]).round(4)),
         "\n## The cells\n",
         f"Value [bootstrap 90 %, {LT['boot_reps']} event draws] [error range: the value minus the fake-world error of the cell at the two ends of the measured rain error].\n",
         md_table(show[["band", "class", "R_cell", "F", "Pv", "links_short", "links_long", "events", "link_events"]].round(3)
                  .rename(columns={"F": "F = intercept (dB)", "Pv": "P = slope (dB/km)", "R_cell": "R_cell (mm/h)"})),
         ("\nDropped (valid in the fake test, not eligible on real data): " + "; ".join(f"{r.band} {r._asdict().get('class', '')}" for r in dropped.itertuples()) + ".\n") if len(dropped) else "\nEvery F-valid or P-valid cell is eligible on real data.\n"]
    for band in cells.band.unique():
        g = cells[cells.band == band]
        gf, gp = g[g.F_valid], g[g.P_valid]
        a.append(f"- **{band}: F valid range {gf.R_cell.min():.2f}–{gf.R_cell.max():.2f} mm/h ({len(gf)} cells); P valid range {gp.R_cell.min():.2f}–{gp.R_cell.max():.2f} mm/h ({len(gp)} cells).**")
    a.append("\n**Below the lowest cell the curves are an interpolation to zero, not a measurement. Above the highest valid cell of a quantity nothing is measured for it** "
             "(F held, P with its last slope: extrapolation, marked). No 5–10 mm/h cell could be judged in the fake test. "
             f"Curves: `{rel(P.E4_CURVES)}` (columns `F_measured`, `P_measured`), figure `{rel(P.E4_FIG)}`.\n")
    a.append("## Against the superseded A3a values (information)\n")
    a.append(md_table(vs.round(3)))
    a.append(f"\nLargest shift {vs.shift_db.abs().max():.2f} dB. (`{rel(P.E4_VS_A3)}`)\n")
    a.append("## The per-link residuals across storms\n")
    a.append(f"Short links (< {LT['short_km']:g} km), residual from the cell lines of the 8 cells valid in BOTH quantities (the criterion as written before Deviation 4), averaged per link and storm: {test['short_links']} links, "
             f"{test['link_storm_values']} link-storm values. Between-link sd **{test['sigma_link_db']:.3f} dB**; storm-to-storm sd inside a link {test['storm_sd_db']:.3f} dB. "
             f"Permutation test (link labels shuffled inside each band, {C['permutations']} times): p = {test['permutation_p']:.4f} (null median {test['null_median_db']:.3f} dB, "
             f"95 % point {test['null_95_db']:.3f} dB).\n")
    a.append(f"**Criterion (decisions.md, 2026-11-11): a per-link factor is needed if the between-link sd is >= {C['link_sd_threshold_db']} dB AND p < {C['p_threshold']}. "
             f"Verdict: {'NEEDED' if test['per_link_factor_needed'] else 'not needed'}.**\n")
    a.append(f"Report only: over {test['links_low_high']} links, the mean residual in the two upper classes against the two lower classes has correlation "
             f"{test['corr_low_high']:.2f} and slope {test['slope_high_on_low']:.2f} (a slope above 1 points to a factor on F, a slope near 1 to an offset). "
             f"(`{rel(P.E4_LINKS)}`, `{rel(P.E4_LINK_TEST)}`, `{rel(P.E4_FIG_LINKS)}`)\n")
    save("\n".join(a) + "\n", P.E4_REPORT)
    print(f"report: {rel(P.E4_REPORT)}", flush=True)


if __name__ == "__main__":
    main()
