"""SUPERSEDED (Detachment correction, 2026-11-11): this step used Part 1 definitions (steady rain, classes, eligible cell) and a
fake truth made of Part 1 closed forms; its check (d) failed. Its outputs are superseded by the rescue steps E1-E4.

Step A3 of the k-surface branch: F(R) (steady WAA) and P(R) (path loss per km) from the model-free length test.

    python -m research.waa_ksurface.a3_length_test

TRAIN only, 28 and 38 GHz, rain = radar / g_e. In steady rain the line (A - c) = I + s L over the links of a band gives the
part that does not grow with the length (I = the WAA) and the part per km (s = the path loss).
(a) the cells (band x rain class) with 200 event-bootstrap draws;
(b) the intercept against how long the rain stayed in the class (is the WAA steady?);
(c) the intercept on the events with a gauge ratio against the events with g_e = 1 (does the rain reference matter?);
(d) a fake check: does the test return a known W and P when the rain error is larger on short links? (STOP if not);
(e) the curves F(R), P(R) per band, with intervals: the product that every later step reads.
"""
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import lfilter

from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import system_id as sid
from research.waa_ksurface import data, lt
from research.waa_ksurface import paths as P

C = data.branch_config()["a3"]
CLASSES = [tuple(c) for c in C["classes_mmh"]]
NAMES = lt.class_names(CLASSES)
BANDS = data.BRANCH_BANDS
MARK = {"28 GHz": ("o", "-"), "38 GHz": ("s", "--")}           # marker and line style per band (no reliance on colour)


# =================================================================================================== (b) the duration
def duration_check(bins, D, le, mult, cells):
    """I per cell for values that stayed at least 5 / 15 / 30 / 60 min in the class, and the paired difference of the rule."""
    fits = {}
    for dur in C["durations_min"]:
        fits[dur] = lt.cell_fits(lt.reduce(bins, D, le, CLASSES, min_run=dur // 5), mult, CLASSES, C["interval"])
    rows = []
    d_hi, d_lo = C["duration_test"]
    elig = cells[cells.eligible]
    for r in elig.itertuples():
        key = (r.band, r._2)
        row = {"band": r.band, "class": r._2}
        for dur, (tab, _) in fits.items():
            m = tab[(tab.band == r.band) & (tab["class"] == r._2)]
            ok = len(m) and np.isfinite(m.I.iloc[0])
            row[f"I_{dur}"] = float(m.I.iloc[0]) if ok else np.nan
            row[f"I_{dur}_lo"], row[f"I_{dur}_hi"] = (float(m.I_lo.iloc[0]), float(m.I_hi.iloc[0])) if ok else (np.nan, np.nan)
            row[f"s_{dur}"] = float(m.s.iloc[0]) if ok else np.nan
            row[f"short_links_{dur}"] = int(m.links_under_1km.iloc[0]) if len(m) else 0
            row[f"link_events_{dur}"] = int(m.link_events.iloc[0]) if len(m) else 0
        a, b = fits[d_hi][1].get(key), fits[d_lo][1].get(key)
        if a is not None and b is not None and np.isfinite(a[:, 0] - b[:, 0]).sum() >= 0.9 * len(mult):
            diff = a[:, 0] - b[:, 0]                            # the same event draws for both durations
            row["diff"], (row["diff_lo"], row["diff_hi"]) = row[f"I_{d_hi}"] - row[f"I_{d_lo}"], np.nanpercentile(diff, C["interval"])
            row["excludes_zero"] = bool(row["diff_lo"] > 0 or row["diff_hi"] < 0)
        else:
            row["diff"], row["diff_lo"], row["diff_hi"], row["excludes_zero"] = np.nan, np.nan, np.nan, False
        rows.append(row)
    tab = pd.DataFrame(rows)
    share = float(tab.excludes_zero.mean())
    fires = share > C["duration_share"]
    need = data.LT["verdict"]["min_short_links"]
    chosen = C["durations_min"][0]
    if fires:                                                   # the longest common duration at which every eligible cell keeps 8 short links
        for dur in C["durations_min"]:
            if (tab[f"short_links_{dur}"] >= need).all() and tab[f"I_{dur}"].notna().all():
                chosen = dur
    return tab, share, fires, chosen, fits


# =================================================================================================== (c) the rain input
def rain_input_check(bins, D, le, mult):
    """The length test on the events with a gauge ratio and on the events with g_e = 1, separately."""
    g = data.g_table()
    used = set(g[g.path_used].event_id)
    is_g = le.event_id.isin(used).to_numpy()
    out = {}
    for name, keep in (("gauge ratio", is_g), ("g_e = 1", ~is_g)):
        b = bins[keep[bins.le_idx.to_numpy()]]
        out[name] = lt.cell_fits(lt.reduce(b, D, le, CLASSES), mult, CLASSES, C["interval"])
    rows = []
    for (band, cls), _ in out["gauge ratio"][0].groupby(["band", "class"], sort=False):
        row = {"band": band, "class": cls}
        for name, short in (("gauge ratio", "gauge"), ("g_e = 1", "ge1")):
            m = out[name][0]
            m = m[(m.band == band) & (m["class"] == cls)]
            for k in ("I", "I_lo", "I_hi", "s", "s_lo", "s_hi", "R_cell"):
                row[f"{k}_{short}"] = float(m[k].iloc[0]) if len(m) and k in m and np.isfinite(m.I.iloc[0]) else np.nan
            row[f"short_links_{short}"] = int(m.links_under_1km.iloc[0]) if len(m) else 0
            row[f"eligible_{short}"] = bool(m.eligible.iloc[0]) if len(m) else False
        row["I_gauge_minus_ge1"] = row["I_gauge"] - row["I_ge1"]
        rows.append(row)
    return pd.DataFrame(rows), out


# =================================================================================================== (d) errors in variables
def radar_gauge_scatter(D, le):
    """sigma_1 and rho of ln(radar / gauge) over the 5 min bins of the gauge-input set where both are >= min_rain_mmh
    (one sub-link per link and event)."""
    E = C["eiv"]
    first = ~le.duplicated(["event_id", "cml_id"]) & (le.gauge_id != "")
    lr_all, pairs = [], []
    for e in np.nonzero(first.to_numpy())[0]:
        o5, n, ro = D["le_boff5"][e], D["le_n"][e], D["le_roff"][e]
        g1 = D["rain_gauge"][ro + sid.PADM:ro + sid.PADM + n].reshape(-1, 5)
        g5 = np.where(np.isfinite(g1).all(axis=1), g1.mean(axis=1), np.nan)
        r5 = D["r5"][o5:o5 + n // 5].astype(float)
        ok = np.isfinite(g5) & np.isfinite(r5) & (g5 >= E["min_rain_mmh"]) & (r5 >= E["min_rain_mmh"])
        lr = np.where(ok, np.log(np.where(ok, r5, 1.0) / np.where(ok, g5, 1.0)), np.nan)
        lr_all.append(lr[ok])
        both = ok[:-1] & ok[1:]
        pairs.append(np.c_[lr[:-1][both], lr[1:][both]])
    lr_all, pairs = np.concatenate(lr_all), np.concatenate(pairs)
    med = np.median(lr_all)
    sigma = 1.4826 * np.median(np.abs(lr_all - med))
    rho = float(np.corrcoef(pairs[:, 0], pairs[:, 1])[0, 1])
    return {"bins": len(lr_all), "link_events": int(first.sum()), "median_ln_ratio": float(med), "sigma_1": float(sigma), "sd_plain": float(lr_all.std()),
            "rho_lag1": rho, "pairs": len(pairs)}


def eiv_check(D, le, cells, noise):
    """Fake data with a known W (independent of L) and a known P, and a rain error that is larger on short links. The length
    test is run on the OBSERVED rain; its I and s are compared with the truth at the cell's R_cell."""
    E, T = C["eiv"], C["eiv"]["truth"]
    sigma1, rho = noise["sigma_1"], noise["rho_lag1"]
    band = le.band_name.to_numpy()
    A_b = np.array([T["A_band_db"][b] for b in band])
    c_b = np.array([T["c_band"][b] for b in band])
    L = le.L_km.to_numpy(float)
    sig_le = sigma1 / np.sqrt(np.maximum(1.0, L / E["pixel_km"]))  # one pixel on a short link, several on a long one
    off_sd = np.nan_to_num(le.offset_std.to_numpy(float))
    n5 = D["le_n"] // 5
    rep_le = np.repeat(np.arange(len(le)), n5)
    r_true = D["r5"].astype(float)
    elig = cells[cells.eligible]
    rows = []
    for rep in range(E["reps"]):
        rng = np.random.default_rng(E["seed"] + rep)
        a_link = np.exp(rng.normal(0, E["a_sd"], D["nl"]))[D["le_link"]]
        z = rng.standard_normal(len(r_true))
        err = np.empty(len(r_true))
        for e in range(len(le)):                                # AR(1) per link-event, stationary, sd sig_le, mean -sd^2/2 (E[exp] = 1)
            o, n = D["le_boff5"][e], n5[e]
            x = z[o:o + n] * np.sqrt(1 - rho ** 2)
            x[0] = z[o]
            err[o:o + n] = sig_le[e] * lfilter([1.0], [1.0, -rho], x) - 0.5 * sig_le[e] ** 2
        r_obs = r_true * np.exp(err)
        w = a_link[rep_le] * A_b[rep_le] * np.power(np.maximum(np.nan_to_num(r_true), 0.0) / 5.0, T["b"])
        y = w + L[rep_le] * c_b[rep_le] * np.nan_to_num(r_true) + E["noise_sd_db"] * rng.standard_normal(len(r_true)) + (off_sd * rng.standard_normal(len(le)))[rep_le]
        Dk = dict(D)
        Dk["r5"], Dk["y5"], Dk["le_c"] = r_obs, y, np.zeros(len(le))
        Dk["w5"] = ((D["w5"] > 0) & np.isfinite(r_true)).astype(float)
        tab, _ = lt.cell_fits(lt.reduce(lt.steady_bins(Dk, le, CLASSES), Dk, le, CLASSES), [], CLASSES)
        for r in elig.itertuples():
            m = tab[(tab.band == r.band) & (tab["class"] == r._2)]
            if not len(m) or not np.isfinite(m.I.iloc[0]):
                rows.append({"rep": rep, "band": r.band, "class": r._2})
                continue
            rc = float(m.R_cell.iloc[0])
            rows.append({"rep": rep, "band": r.band, "class": r._2, "R_cell": rc, "I": float(m.I.iloc[0]), "s": float(m.s.iloc[0]),
                         "F_true": T["A_band_db"][r.band] * (rc / 5.0) ** T["b"], "P_true": T["c_band"][r.band] * rc,
                         "short_links": int(m.links_under_1km.iloc[0])})
        print(f"   eiv replicate {rep + 1} / {E['reps']}", flush=True)
    reps = pd.DataFrame(rows)
    reps["I_err"], reps["s_err_rel"] = reps.I - reps.F_true, (reps.s - reps.P_true) / reps.P_true
    g = reps.groupby(["band", "class"], sort=False)
    out = g.agg(R_cell=("R_cell", "median"), F_true=("F_true", "median"), I=("I", "median"), I_err=("I_err", "median"), I_err_sd=("I_err", "std"),
                P_true=("P_true", "median"), s=("s", "median"), s_err_rel=("s_err_rel", "median"), s_err_rel_sd=("s_err_rel", "std"),
                reps_fitted=("I", "count")).reset_index()
    out["I_ok"] = out.I_err.abs() <= E["pass"]["intercept_db"]
    out["s_ok"] = out.s_err_rel.abs() <= E["pass"]["slope_rel"]
    return out, bool((out.I_ok & out.s_ok).all() and (out.reps_fitted == E["reps"]).all())


# =================================================================================================== (e) the curves
def build_curves(cells, draws, gauge_cells):
    """F and P per band on the dense R grid, with the 5-95 % band of the bootstrap draws, and the F of the gauge-ratio events."""
    G = C["grid"]
    grid = np.r_[0.0, np.geomspace(G["r_min"], G["r_max"], G["n"])]
    rows, knots = [], []
    for band in BANDS:
        e = cells[(cells.band == band) & cells.eligible].sort_values("R_cell")
        rc = e.R_cell.to_numpy()
        col = {}
        for kind, name in (("F", "I"), ("P", "s")):
            wts = 1.0 / np.maximum(e[f"{name}_se"].to_numpy(), 1e-6) ** 2
            f, v = lt.monotone_curve(rc, e[name].to_numpy(), wts, kind)
            col[kind] = f(grid)
            k = ("I", "s").index(name)
            D = np.array([draws[(band, c)][:, k] for c in e["class"]])           # cells x draws
            good = np.isfinite(D).all(axis=0)
            per_draw = np.array([lt.monotone_curve(rc, D[:, j], wts, kind)[0](grid) for j in np.nonzero(good)[0]])
            col[f"{kind}_lo"], col[f"{kind}_hi"] = np.percentile(per_draw, C["interval"], axis=0)
            for c_, r_, raw, mono, lo, hi in zip(e["class"], rc, e[name], v, e[f"{name}_lo"], e[f"{name}_hi"]):
                knots.append({"band": band, "quantity": kind, "class": c_, "R_cell": r_, "cell_value": raw, "curve_value": mono, "cell_lo": lo, "cell_hi": hi,
                              "changed_by_monotone_step": bool(abs(raw - mono) > 1e-9)})
        ge = gauge_cells[(gauge_cells.band == band) & gauge_cells.eligible].sort_values("R_cell")
        if len(ge) >= 2:
            col["F_gauge_events"] = lt.monotone_curve(ge.R_cell.to_numpy(), ge.I.to_numpy(), 1.0 / np.maximum(ge.I_se.to_numpy(), 1e-6) ** 2, "F")[0](grid)
        else:
            col["F_gauge_events"] = np.full(len(grid), np.nan)
        t = pd.DataFrame({"band": band, "R": grid, **col})
        t["measured"] = (t.R >= rc.min()) & (t.R <= rc.max())   # inside the range of the eligible cells
        rows.append(t)
    return pd.concat(rows, ignore_index=True), pd.DataFrame(knots)


# =================================================================================================== figures
def fig_curves(curves, cells):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, kind, name, unit in ((axes[0], "F", "I", "dB"), (axes[1], "P", "s", "dB/km")):
        for band in BANDS:
            mk, ls = MARK[band]
            g = curves[(curves.band == band) & (curves.R <= 20)]
            e = cells[(cells.band == band) & cells.eligible]
            ax.plot(g.R, g[kind], ls, color=INK, lw=1.4, label=f"{band}: curve")
            ax.fill_between(g.R, g[f"{kind}_lo"], g[f"{kind}_hi"], color=INK2, alpha=0.15, hatch="//" if band == BANDS[0] else "\\\\", lw=0)
            ax.errorbar(e.R_cell, e[name], yerr=[e[name] - e[f"{name}_lo"], e[f"{name}_hi"] - e[name]], fmt=mk, color=INK, mfc="white", capsize=3, lw=0.8,
                        label=f"{band}: cells (90 %)")
            lo = e.R_cell.min()
            ax.axvline(lo, color=INK2, lw=0.6, ls=":")
        ax.set_xscale("symlog", linthresh=0.1)
        ax.set_xlabel("rain R (mm/h)")
        ax.set_ylabel(f"{kind}(R) ({unit})")
        ax.set_title(("F: steady WAA (intercept of the length test)" if kind == "F" else "P: path loss per km (slope of the length test)")
                     + "\nleft of the dotted line: interpolation to zero, not measured", fontsize=9)
        ax.legend(fontsize=7)
        style_axes(ax)
    fig.tight_layout()
    return fig


def fig_duration(dur):
    fig, axes = plt.subplots(1, len(BANDS), figsize=(11, 4), sharey=True)
    styles = ["-", "--", ":", "-.", (0, (5, 1)), (0, (1, 1)), (0, (3, 1, 1, 1))]
    marks = ["o", "s", "^", "v", "D", "x", "P"]
    for ax, band in zip(axes, BANDS):
        for i, r in enumerate(dur[dur.band == band].itertuples()):
            d = C["durations_min"]
            v = [getattr(r, f"I_{k}") for k in d]
            ax.plot(d, v, marker=marks[i % 7], ls=styles[i % 7], color=INK, mfc="white", lw=1, label=f"{r._2} mm/h")
            ax.annotate(r._2, (d[0], v[0]), textcoords="offset points", xytext=(-4, 0), ha="right", fontsize=7)
        ax.set_title(f"{band}: intercept against the time in the rain class", fontsize=9)
        ax.set_xlabel("at least this many minutes in the class")
        ax.set_xticks(C["durations_min"])
        style_axes(ax)
    axes[0].set_ylabel("intercept I (dB)")
    axes[-1].legend(fontsize=7, title="rain class")
    fig.tight_layout()
    return fig


def fig_input(inp):
    fig, axes = plt.subplots(1, len(BANDS), figsize=(11, 4), sharey=True)
    for ax, band in zip(axes, BANDS):
        g = inp[inp.band == band]
        for short, mk, ls, lab in (("gauge", "o", "-", "events with a gauge ratio"), ("ge1", "s", "--", "events with g_e = 1")):
            m = g[g[f"I_{short}"].notna()]
            ax.errorbar(m[f"R_cell_{short}"], m[f"I_{short}"], yerr=[m[f"I_{short}"] - m[f"I_lo_{short}"], m[f"I_hi_{short}"] - m[f"I_{short}"]], fmt=mk, ls=ls,
                        color=INK, mfc="white", capsize=3, lw=0.9, label=lab)
        ax.set_xscale("log")
        ax.set_title(f"{band}: intercept by rain reference", fontsize=9)
        ax.set_xlabel("rain R (mm/h)")
        style_axes(ax)
    axes[0].set_ylabel("intercept I (dB)")
    axes[0].legend(fontsize=7)
    fig.tight_layout()
    return fig


def fig_eiv(eiv):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    E = C["eiv"]["pass"]
    for band in BANDS:
        mk, ls = MARK[band]
        g = eiv[eiv.band == band]
        axes[0].errorbar(g.R_cell, g.I_err, yerr=g.I_err_sd, fmt=mk, ls=ls, color=INK, mfc="white", capsize=3, lw=0.9, label=band)
        axes[1].errorbar(g.R_cell, 100 * g.s_err_rel, yerr=100 * g.s_err_rel_sd, fmt=mk, ls=ls, color=INK, mfc="white", capsize=3, lw=0.9, label=band)
    for ax, lim, lab, ttl in ((axes[0], E["intercept_db"], "intercept minus true F (dB)", "intercept error"), (axes[1], 100 * E["slope_rel"], "slope error (% of true P)", "slope error")):
        ax.axhline(lim, color=INK2, ls=":", lw=1)
        ax.axhline(-lim, color=INK2, ls=":", lw=1)
        ax.axhline(0, color=INK2, lw=0.6)
        ax.set_xscale("log")
        ax.set_xlabel("observed R of the cell (mm/h)")
        ax.set_ylabel(lab)
        ax.set_title(f"fake check, rain error larger on short links: {ttl}\n(median of the replicates, bar = their sd, dotted = the pass limit)", fontsize=9)
        ax.legend(fontsize=7)
        style_axes(ax)
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    D0, D1, le, _ = data.branch_train()
    n_ev = int(D1["le_event"].max()) + 1
    print(f"data: {len(le)} TRAIN link-events, {le.cml_id.nunique()} links, {n_ev} events (28 and 38 GHz)", flush=True)
    mult = lt.multiplicities(n_ev, C["boot_reps"], C["seed"])
    bins = lt.steady_bins(D1, le, CLASSES)

    # ---- (a) the cells
    cells, draws = lt.cell_fits(lt.reduce(bins, D1, le, CLASSES), mult, CLASSES, C["interval"])
    print("(a) cells:\n" + cells[["band", "class", "R_cell", "I", "I_lo", "I_hi", "s", "s_lo", "s_hi", "links_under_1km", "eligible"]].round(3).to_string(), flush=True)

    # ---- (b) the duration
    dur, share, fires, chosen, dfits = duration_check(bins, D1, le, mult, cells)
    print(f"(b) I({C['duration_test'][0]}) - I({C['duration_test'][1]}) excludes zero in {dur.excludes_zero.sum()} of {len(dur)} eligible cells "
          f"(share {share:.2f}); the rule {'FIRES' if fires else 'does not fire'}; duration used: >= {chosen} min", flush=True)
    used_cells, used_draws = (cells, draws) if chosen == C["durations_min"][0] else dfits[chosen]
    if chosen != C["durations_min"][0]:                         # eligibility stays that of (a); the values are those of the chosen duration
        used_cells = used_cells.copy()
        used_cells["eligible"] = [bool(cells[(cells.band == r.band) & (cells["class"] == r._2)].eligible.iloc[0]) for r in used_cells.itertuples()]

    # ---- (c) the rain input
    inp, inp_fits = rain_input_check(bins, D1, le, mult)

    # ---- (d) errors in variables
    noise = radar_gauge_scatter(D1, le)
    print(f"(d) radar-vs-gauge scatter: sigma_1 {noise['sigma_1']:.3f}, rho {noise['rho_lag1']:.3f} ({noise['bins']} bins)", flush=True)
    eiv, eiv_pass = eiv_check(D1, le, cells, noise)
    print("(d) " + eiv[["band", "class", "R_cell", "F_true", "I", "I_err", "P_true", "s", "s_err_rel", "I_ok", "s_ok"]].round(3).to_string(), flush=True)
    print(f"(d) verdict: {'PASS' if eiv_pass else 'STOP'}", flush=True)

    # ---- (e) the curves
    curves, knots = build_curves(used_cells, used_draws, inp_fits["gauge ratio"][0])

    for t, f in ((cells, P.A3_CELLS), (dur, P.A3_DURATION), (inp, P.A3_INPUT), (eiv, P.A3_EIV), (pd.DataFrame([noise]), P.A3_EIV_NOISE),
                 (knots, P.A3_KNOTS), (curves, P.A3_CURVES)):
        save(t.round(5), f, index=False)
    for fig, f in ((fig_curves(curves, used_cells), P.A3_FIG_CURVES), (fig_duration(dur), P.A3_FIG_DURATION), (fig_input(inp), P.A3_FIG_INPUT),
                   (fig_eiv(eiv), P.A3_FIG_EIV)):
        save(fig, f, dpi=120)
        plt.close(fig)

    write_report(le, n_ev, cells, dur, share, fires, chosen, inp, noise, eiv, eiv_pass, knots, curves)
    if not eiv_pass:
        sys.exit(1)


def write_report(le, n_ev, cells, dur, share, fires, chosen, inp, noise, eiv, eiv_pass, knots, curves):
    E = C["eiv"]
    a = ["# A3: F(R) and P(R) from the length test (TRAIN, 28 and 38 GHz)\n",
         f"Written by `research/waa_ksurface/a3_length_test.py`. {len(le)} link-events, {le.cml_id.nunique()} links, {n_ev} events. Rain = radar / g_e. "
         "32 GHz is excluded from the branch (rule A1e). VAL and TEST were not opened.\n",
         "## (a) The cells\n",
         f"Line (A − c) = I + s·L per band and rain class, in steady rain (Part 1 rule), {C['boot_reps']} event-bootstrap draws, 90 % intervals. "
         f"Eligible = fitted and at least {data.LT['verdict']['min_short_links']} links below 1 km. (`{rel(P.A3_CELLS)}`)\n"]
    t = cells.assign(I_=[f"{r.I:.3f} [{r.I_lo:.3f}, {r.I_hi:.3f}]" if np.isfinite(r.I) else "" for r in cells.itertuples()],
                     s_=[f"{r.s:.3f} [{r.s_lo:.3f}, {r.s_hi:.3f}]" if np.isfinite(r.s) else "" for r in cells.itertuples()])
    a.append(md_table(t[["band", "class", "R_cell", "I_", "s_", "link_events", "links", "events", "links_under_1km", "eligible"]].round(3)
                      .rename(columns={"I_": "intercept I (dB)", "s_": "slope s (dB/km)", "R_cell": "R_cell (mm/h)"})))
    new = cells[cells["class"].isin(NAMES[:2])]
    a.append(f"\nThe two added classes ({NAMES[0]}, {NAMES[1]} mm/h): " + "; ".join(f"{r.band} {r._2}: {'eligible' if r.eligible else 'NOT eligible'}" for r in new.itertuples()) + ". "
             "A bin of these classes needs the 20 min before it to be wet (R >= 0.3 mm/h), so their values come from rain that is falling off; "
             "their intercept can sit above the steady level of that rain.\n")
    a.append("## (b) Is the WAA steady? The intercept against the time in the rain class\n")
    cols = ["band", "class"] + [f"I_{d}" for d in C["durations_min"]] + [f"short_links_{d}" for d in C["durations_min"]] + ["diff", "diff_lo", "diff_hi", "excludes_zero"]
    a.append(md_table(dur[cols].round(3)))
    d_hi, d_lo = C["duration_test"]
    a.append(f"\n`diff` = I(>= {d_hi} min) − I(>= {d_lo} min), paired over the same event draws, 90 % interval. It excludes zero in {int(dur.excludes_zero.sum())} of "
             f"{len(dur)} eligible cells (share {share:.2f}; the rule fires above {C['duration_share']}). **The rule {'FIRES' if fires else 'does not fire'}: F is taken from "
             f"{'the values with at least ' + str(chosen) + ' min in the class' if fires else 'all steady values'}.** (`{rel(P.A3_DURATION)}`, `{rel(P.A3_FIG_DURATION)}`)\n")
    a.append("## (c) The rain reference: events with a gauge ratio against events with g_e = 1 (report only)\n")
    a.append(md_table(inp[["band", "class", "I_gauge", "I_lo_gauge", "I_hi_gauge", "I_ge1", "I_lo_ge1", "I_hi_ge1", "I_gauge_minus_ge1", "s_gauge", "s_ge1",
                           "short_links_gauge", "short_links_ge1"]].round(3)))
    a.append(f"\nThe F of the gauge-ratio events is kept as the column `F_gauge_events` of the curve table (third sensitivity curve of rule A3f). "
             f"(`{rel(P.A3_INPUT)}`, `{rel(P.A3_FIG_INPUT)}`)\n")
    a.append(f"## (d) Errors-in-variables fake check: {'PASS' if eiv_pass else 'STOP'}\n")
    a.append(f"Measured radar-vs-gauge scatter on the gauge-input set (5 min bins with both >= {E['min_rain_mmh']} mm/h, {noise['bins']} bins): σ_1 = {noise['sigma_1']:.3f} "
             f"(1.4826·MAD of ln(radar / gauge); plain sd {noise['sd_plain']:.3f}), lag-1 correlation ρ = {noise['rho_lag1']:.3f}. The fake rain error has the sd "
             f"σ_1 / √max(1, L / {E['pixel_km']:g} km), is AR(1) with ρ, and has mean one (E[exp(e)] = 1). Truth: F = A_band·(R / 5)^{E['truth']['b']}, P = c_band·R "
             f"(Part 1 closed forms). {E['reps']} replicates; the numbers are the medians. Pass: |I − F_true| <= {E['pass']['intercept_db']} dB and "
             f"|s − P_true| / P_true <= {100 * E['pass']['slope_rel']:.0f} % in every eligible cell.\n")
    a.append(md_table(eiv[["band", "class", "R_cell", "F_true", "I", "I_err", "P_true", "s", "s_err_rel", "I_ok", "s_ok"]].round(3)))
    a.append(f"\nLargest intercept error {eiv.I_err.abs().max():.3f} dB; largest slope error {100 * eiv.s_err_rel.abs().max():.1f} %. "
             f"(`{rel(P.A3_EIV)}`, `{rel(P.A3_FIG_EIV)}`)\n")
    a.append("## (e) The curves\n")
    a.append(f"Per band, a monotone cubic (PCHIP) through (0, 0) and the eligible cell values; intervals from the bootstrap draws. Product: `{rel(P.A3_CURVES)}` "
             f"(read with `research.waa_ksurface.lt.Curves`). The cell values (`{rel(P.A3_KNOTS)}`):\n")
    a.append(md_table(knots.round(4)))
    for band in BANDS:
        g = curves[(curves.band == band) & curves.measured]
        a.append(f"\n- **{band}: valid range {g.R.min():.2f}–{g.R.max():.2f} mm/h** (lowest to highest eligible R_cell).")
    lowest = knots.sort_values("R_cell").groupby("band").first()
    a.append("\n\n**Below the lowest eligible class the curve is an interpolation to zero, not a measurement** (lowest measured class: "
             + ", ".join(f"{b}: {r['class']} mm/h at R_cell {r.R_cell:.2f}" for b, r in lowest.iterrows()) + "). Above the highest cell F is held constant and P continues "
             "with its last slope: extrapolation.\n")
    ch = knots[knots.changed_by_monotone_step]
    a.append(f"The monotone step changed {len(ch)} cell value(s)" + (": " + "; ".join(f"{r.band} {r.quantity} {r._3}: {r.cell_value:.3f} -> {r.curve_value:.3f}" for r in ch.itertuples()) if len(ch) else "") + ".\n")
    if not eiv_pass:
        a.append("\n## STOP\n\nThe errors-in-variables check failed its rule. Nothing after A3 is run.\n")
    save("\n".join(a) + "\n", P.A3_REPORT)
    print(f"report: {rel(P.A3_REPORT)}", flush=True)


if __name__ == "__main__":
    main()
