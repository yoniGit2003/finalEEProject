"""Steps T2 and T3 of the k-surface branch: candidate (iv), and the E2 fake test rerun at the conditional rain error.

    python -m research.waa_ksurface.t3_rerun

E2 used the error of every single 5 min bin. The length test only meets the error of a window label in steady rain, which
step T1 measured: a per-event scale part (sd 0.45), and inside an event a slowly changing part of which the link's own
share has sd 0.18-0.26. Candidates (i), (ii), (iii) are rerun with that error. Candidate (iv) sorts the bins by the link's
own City gauge instead of the radar, and has its own fake world (gauge error, point-vs-path representativeness, gauge steps).
Verdicts are per cell: a cell is valid if its intercept and slope are inside tolerance at BOTH ends of the error range.
The candidate with the most valid cells is adopted if it has at least 8, of which at least 3 below 2 mm/h.
Fake data only; real data is read for the counts of candidate (iv).
"""
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import candidates as cand, data, fakeworld as fw, ltest
from research.waa_ksurface import paths as P
from research.waa_ksurface.e2_candidates import sizes_and_rule

LT, FK = ltest.config(), fw.config()
T3 = data.branch_config()["t3"]
NAMES = {"i": "(i) gauge-ratio events, plain", "ii": "(ii) instrument binning (area-mean radar)", "iii": "(iii) forward-bias correction",
         "iv": "(iv) gauge-binned"}
ORDER = ["i", "ii", "iii", "iv"]
ENDS = {"i": ("lo", "hi"), "ii": ("lo", "hi"), "iii": ("lo", "hi"), "iv": ("normal", "doubled")}
CLASS = ltest.class_names(LT)
BELOW2 = 4                                                      # classes 0.1-0.2, 0.2-0.5, 0.5-1, 1-2 lie below 2 mm/h


def radar_candidates(D, le, n_cells, sizes, rule, n_short, l_med):
    """(i), (ii), (iii) in the fake world with the slow error of T1, at sigma_c_lo and sigma_c_hi."""
    out = {}
    for end, sig in (("lo", T3["sigma_c_lo"]), ("hi", T3["sigma_c_hi"])):
        tabs = {c: [] for c in ("i", "ii", "iii")}
        for rep in range(FK["reps"]):
            Dk = fw.make(D, le, n_cells, sizes, 0.0, FK, FK["seed"] + rep)             # truth, memory, noise; no rain error yet
            err = fw.rain_error_slow(le, n_cells, sig, T3["slow_error"], np.random.default_rng(FK["seed"] + 500 + rep))
            Dk["r5"] = Dk["r_true"] * np.exp(err)
            tabs["i"].append(cand.plain(Dk, le, rule, n_short)[0])
            tabs["ii"].append(cand.instrument(Dk, le, rule, n_short))
            tabs["iii"].append(cand.corrected_slow(Dk, le, n_cells, sig, rule, n_short, FK["seed"] + 1000 + rep, T3))
        for c in tabs:
            out[(c, end)] = fw.judge(tabs[c], FK, min_reps=T3["min_reps_judged"], l_med=l_med)[0]
            print(f"   {NAMES[c]:45s} sigma_c {end}: {len(out[(c, end)])} cells judged", flush=True)
    return out


def gauge_candidate(sizes, n_short, l_med):
    """(iv): the real-data counts, and its fake world at normal and doubled representativeness."""
    Dg, leg, n_cells, g5, step, dist = data.gauge_links_train()
    real = cand.gauge_binned(Dg, leg, g5, step, n_short)[0]
    info = {"link_events": len(leg), "links": leg.cml_id.nunique(), "links_short": leg[leg.L_km < LT["short_km"]].cml_id.nunique(), "events": leg.event_id.nunique()}
    out = {}
    if not real.eligible.any():
        return real, info, out
    for end, factor in zip(ENDS["iv"], T3["gauge"]["rep_factors"]):
        tabs = []
        for rep in range(FK["reps"]):
            Dk = fw.make(Dg, leg, n_cells, sizes, 0.0, FK, FK["seed"] + rep)
            gf = fw.fake_gauge(Dk["r_true"], Dg, leg, step, dist, T3["gauge"], np.random.default_rng(FK["seed"] + 700 + rep), rep_factor=factor)
            tabs.append(cand.gauge_binned(Dk, leg, gf, step, n_short)[0])
        out[("iv", end)] = fw.judge(tabs, FK, min_reps=T3["min_reps_judged"], l_med=l_med)[0]
        print(f"   {NAMES['iv']:45s} representativeness {end}: {len(out[('iv', end)])} cells judged", flush=True)
    return real, info, out


def valid_cells(judged):
    """Per candidate and cell: judged and inside at both ends -> valid."""
    rows = []
    for c in ORDER:
        a, b = (judged.get((c, e), pd.DataFrame()) for e in ENDS[c])
        keys = sorted(set(map(tuple, a[["band", "cls"]].to_numpy())) | set(map(tuple, b[["band", "cls"]].to_numpy()))) if len(a) or len(b) else []
        for band, cls in keys:
            row = {"candidate": c, "band": band, "class": CLASS[int(cls)], "cls": int(cls)}
            ok = True
            for end, t in zip(ENDS[c], (a, b)):
                m = t[(t.band == band) & (t.cls == cls)] if len(t) else t
                if len(m):
                    r = m.iloc[0]
                    row.update({f"I_err_{end}": r.I_err, f"s_err_{end}": r.s_err_abs, f"inside_{end}": bool(r.I_ok and r.s_ok), f"R_cell_{end}": r.R_cell})
                    ok &= bool(r.I_ok and r.s_ok)
                else:
                    row.update({f"I_err_{end}": np.nan, f"s_err_{end}": np.nan, f"inside_{end}": False, f"R_cell_{end}": np.nan})
                    ok = False
            row["valid"] = ok
            rows.append(row)
    return rows


def fig_cells(judged):
    fig, axes = plt.subplots(2, 4, figsize=(16, 7), sharey=True, sharex=True)
    style = {"28 GHz": ("o", "-"), "38 GHz": ("s", "--")}
    for j, c in enumerate(ORDER):
        for i, end in enumerate(ENDS[c]):
            ax = axes[i, j]
            g = judged.get((c, end), pd.DataFrame())
            for band, (mk, ls) in style.items():
                m = g[g.band == band] if len(g) else g
                if len(m):
                    ax.plot(m.R_cell, m.I_err, marker=mk, ls=ls, color=INK, mfc="white", lw=1, label=band)
            for lim in (FK["pass"]["intercept_db"], -FK["pass"]["intercept_db"]):
                ax.axhline(lim, color=INK2, ls=":", lw=1)
            ax.axhline(0, color=INK2, lw=0.6)
            ax.set_xscale("log")
            ax.set_title(f"{NAMES[c]}\nend: {end}", fontsize=8)
            style_axes(ax)
            if i == 1:
                ax.set_xlabel("R of the cell (mm/h)")
            if j == 0:
                ax.set_ylabel("intercept minus true F (dB)")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("T3 fake test at the conditional rain error: intercept error per cell (dotted = ±0.15 dB)", fontsize=10)
    fig.tight_layout()
    return fig


def main():
    D, le, n_cells = data.gauge_ratio_train()
    sizes, rule, rule_name, n_short, l_med = sizes_and_rule()
    print(f"radar candidates: {len(le)} link-events, {le.cml_id.nunique()} links, {le.event_id.nunique()} events; rule {rule_name}; n_short {n_short}", flush=True)
    judged = radar_candidates(D, le, n_cells, sizes, rule, n_short, l_med)
    real_iv, info_iv, j_iv = gauge_candidate(sizes, n_short, l_med)
    judged.update(j_iv)
    save(real_iv[["band", "class", "R_cell", "link_events", "links", "events", "bins", "links_short", "links_long", "fitted", "eligible"]].round(4), P.T2_REAL_CELLS, index=False)
    floors = []
    for fl in range(n_short, 0, -1):                            # report only: how many cells of (iv) each floor would leave (the floor is NOT lowered)
        floors.append({"floor_short_links": fl, "cells": int((real_iv.fitted & (real_iv.links_short >= fl) & (real_iv.links_long >= fl)).sum())})
    floors = pd.DataFrame(floors)

    cells = pd.DataFrame(valid_cells(judged))
    ver = []
    for c in ORDER:
        g = cells[cells.candidate == c] if len(cells) else cells
        v = g[g.valid] if len(g) else g
        ver.append({"candidate": c, "cells_judged_at_both_ends": int((g[[f"inside_{e}" for e in ENDS[c]]].notna().all(axis=1) & g[f"R_cell_{ENDS[c][0]}"].notna()
                                                                      & g[f"R_cell_{ENDS[c][1]}"].notna()).sum()) if len(g) else 0,
                    "valid_cells": len(v), "valid_below_2mmh": int((v.cls < BELOW2).sum()) if len(v) else 0,
                    "largest_I_err_db": float(np.nanmax(np.abs(g[[f"I_err_{e}" for e in ENDS[c]]].to_numpy()))) if len(g) else np.nan})
    ver = pd.DataFrame(ver)
    A = T3["adopt"]
    ver["qualifies"] = (ver.valid_cells >= A["min_valid_cells"]) & (ver.valid_below_2mmh >= A["min_valid_below_2mmh"])
    q = ver[ver.qualifies]
    adopted, tie = None, False
    if len(q):
        best = q[q.valid_cells == q.valid_cells.max()]
        tie = len(best) > 1
        adopted = best.candidate.iloc[0] if not tie else None
    ver["adopted"] = ver.candidate == adopted
    parts = [t.assign(candidate=c, end=e) for (c, e), t in judged.items() if len(t)]
    fake = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    save(fake.round(5), P.T3_FAKE, index=False)
    save(cells.round(5), P.T3_CELLS, index=False)
    save(ver.round(5), P.T3_VERDICT, index=False)
    fig = fig_cells(judged)
    save(fig, P.T3_FIG, dpi=120)
    plt.close(fig)
    print(ver.to_string(), flush=True)

    S = T3["slow_error"]
    a = ["# T2 / T3: candidate (iv), and the fake test at the conditional rain error\n",
         f"Written by `research/waa_ksurface/t3_rerun.py`. Fake data only (real data is read for the counts of candidate (iv)). Rule {rule_name}, n_short = {n_short}, "
         f"tolerances {FK['pass']['intercept_db']} dB and {FK['pass']['intercept_db'] / l_med:.4f} dB/km (Deviation 3). Rain error of (i)-(iii), from T1: per-event scale part sd "
         f"{S['event_sd']}; inside an event a slow error ({S['window_corr']} correlated over 30 min) with the link's own part sd {T3['sigma_c_lo']} (low end) and {T3['sigma_c_hi']} "
         f"(high end) and a part shared by all links {np.sqrt(S['shared_share'] / (1 - S['shared_share'])):.2f} times that. {FK['reps']} replicates; a cell is judged if eligible "
         f"in at least {T3['min_reps_judged']}. VAL and TEST were not opened.\n",
         "## T2: candidate (iv) on real data\n",
         f"{info_iv['link_events']} TRAIN link-events of {info_iv['links']} links with a City gauge within 1 km ({info_iv['links_short']} of them below 1 km), "
         f"{info_iv['events']} events (all TRAIN events). (`{rel(P.T2_REAL_CELLS)}`)\n",
         md_table(real_iv[["band", "class", "R_cell", "link_events", "links", "events", "links_short", "links_long", "fitted", "eligible"]].round(3)),
         "\nCells that each floor of short links (and as many longer links) would leave. Report only: the floor is not lowered.\n", md_table(floors),
         ("\n**No cell of (iv) is eligible with the floor of 5: candidate (iv) stops here.**\n" if not real_iv.eligible.any() else ""),
         "## T3: verdict per candidate\n",
         md_table(ver.assign(candidate=ver.candidate.map(NAMES)).round(3)),
         f"\nAdoption rule: the most valid cells, with at least {A['min_valid_cells']} valid and at least {A['min_valid_below_2mmh']} of them below 2 mm/h. "
         f"**Adopted: {NAMES[adopted] if adopted else ('none (a tie, to be broken on real data)' if tie else 'none')}.**\n",
         "## Every cell\n"]
    for c in ORDER:
        g = cells[cells.candidate == c] if len(cells) else cells
        a.append(f"### {NAMES[c]}\n")
        a.append(md_table(g.drop(columns=["candidate", "cls"]).round(3)) if len(g) else "no cell judged")
        a.append("")
    a.append(f"(`{rel(P.T3_CELLS)}`, `{rel(P.T3_FAKE)}`, `{rel(P.T3_VERDICT)}`, `{rel(P.T3_FIG)}`)\n")
    save("\n".join(a) + "\n", P.T3_REPORT)
    print(f"adopted: {adopted} (tie: {tie}); report: {rel(P.T3_REPORT)}", flush=True)
    if adopted is None:
        sys.exit(2 if tie else 1)


if __name__ == "__main__":
    main()
