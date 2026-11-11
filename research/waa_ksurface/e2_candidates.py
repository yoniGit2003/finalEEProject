"""Step E2 of the k-surface branch: which length test survives the measured rain error? (fake data only)

    python -m research.waa_ksurface.e2_candidates

Three candidates (library `candidates.py`) are run in the fake world of step E3 (known F and P, a store with memory, the
measured noise) with the rain error measured in step E1: a scale part per link-event, a per-bin part of which 46 % is shared
by all links at one time, lag-1 correlation 0.35. The per-bin sd is set to both ends of the measured range (0.48 and 0.87).
A candidate passes if it returns the known F (within 0.15 dB) and P (within 0.15 dB / L_med) in all judged cells but at most
one, at BOTH ends. The simplest passing candidate is adopted. A sweep over smaller rain errors shows how much error each
candidate can take. If none passes the script exits with code 1 (the rule: STOP).
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

LT, FK = ltest.config(), fw.config()
NAMES = {"i": "(i) gauge-ratio events, plain", "ii": "(ii) instrument binning (area-mean radar)", "iii": "(iii) forward-bias correction"}
ORDER = ["i", "ii", "iii"]                                      # from the simplest
SWEEP = [0.0, 0.25, 0.5, 0.75]                                  # the whole rain error scaled by f (f = 1 is the sigma_lo run)


def sizes_and_rule():
    s = pd.read_csv(P.E3_SIZES).set_index("quantity").value
    sizes = {k: float(s[k]) for k in ("white_sd_db", "offset_sd_db", "link_sd_db")}
    return sizes, LT["rules"][s["adopted_rule"]], s["adopted_rule"], int(float(s["n_short"])), float(s["l_med_km"])


def run_config(D, le, n_cells, sizes, rule, n_short, l_med, sigma, scale_sd, label):
    """All candidates on the same fake replicates of one rain-error size. Returns (cell table, verdict rows)."""
    tabs = {c: [] for c in ORDER}
    for rep in range(FK["reps"]):
        Dk = fw.make(D, le, n_cells, sizes, sigma, FK, FK["seed"] + rep, scale_sd=scale_sd)
        tabs["i"].append(cand.plain(Dk, le, rule, n_short)[0])
        tabs["ii"].append(cand.instrument(Dk, le, rule, n_short))
        tabs["iii"].append(cand.corrected(Dk, le, n_cells, sigma, rule, n_short, FK["seed"] + 1000 + rep, scale_sd=scale_sd))
    parts, ver = [], []
    for c in ORDER:
        tab, n_fail, ok = fw.judge(tabs[c], FK, l_med=l_med)
        if len(tab):
            parts.append(tab.assign(candidate=c, config=label, sigma=sigma, scale_sd=scale_sd))
        ver.append({"candidate": c, "config": label, "sigma": sigma, "scale_sd": scale_sd, "cells_judged": len(tab), "cells_inside": len(tab) - n_fail,
                    "I_outside": int((~tab.I_ok).sum()) if len(tab) else 0, "s_outside": int((~tab.s_ok).sum()) if len(tab) else 0,
                    "largest_I_err_db": float(tab.I_err.abs().max()) if len(tab) else np.nan,
                    "largest_s_err_db_km": float(tab.s_err_abs.abs().max()) if len(tab) else np.nan, "passes": ok})
        print(f"   {label:9s} {NAMES[c]:45s} {ver[-1]['cells_inside']} / {ver[-1]['cells_judged']} inside, largest |I err| {ver[-1]['largest_I_err_db']:.3f} dB "
              f"-> {'PASS' if ok else 'fail'}", flush=True)
    return parts, ver


def fig_errors(cells):
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), sharey="row", sharex=True)
    style = {"28 GHz": ("o", "-"), "38 GHz": ("s", "--")}
    for j, c in enumerate(ORDER):
        for i, (cfg, ttl) in enumerate((("sigma_lo", "rain error at the low end"), ("sigma_hi", "rain error at the high end"))):
            ax = axes[i, j]
            g = cells[(cells.candidate == c) & (cells.config == cfg)]
            for band, (mk, ls) in style.items():
                m = g[g.band == band]
                ax.plot(m.R_cell, m.I_err, marker=mk, ls=ls, color=INK, mfc="white", lw=1, label=band)
            for lim in (FK["pass"]["intercept_db"], -FK["pass"]["intercept_db"]):
                ax.axhline(lim, color=INK2, ls=":", lw=1)
            ax.axhline(0, color=INK2, lw=0.6)
            ax.set_xscale("log")
            ax.set_title(f"{NAMES[c]}\n{ttl}", fontsize=8)
            style_axes(ax)
            if i == 1:
                ax.set_xlabel("R of the cell (mm/h)")
            if j == 0:
                ax.set_ylabel("intercept minus true F (dB)")
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("E2 fake test: intercept error per cell (dotted = the pass limit ±0.15 dB)", fontsize=10)
    fig.tight_layout()
    return fig


def main():
    D, le, n_cells = data.gauge_ratio_train()
    sizes, rule, rule_name, n_short, l_med = sizes_and_rule()
    R = FK["rain_error"]
    print(f"data: {len(le)} link-events, {le.cml_id.nunique()} links, {le.event_id.nunique()} events; rule {rule_name} {rule}; n_short {n_short}; "
          f"slope tolerance {FK['pass']['intercept_db'] / l_med:.4f} dB/km", flush=True)
    parts, ver = [], []
    configs = [("sigma_lo", R["sigma_lo"], R["scale_sd"]), ("sigma_hi", R["sigma_hi"], R["scale_sd"])] + \
              [(f"f={f:g}", f * R["sigma_lo"], f * R["scale_sd"]) for f in SWEEP]
    for label, sigma, scale in configs:
        p, v = run_config(D, le, n_cells, sizes, rule, n_short, l_med, sigma, scale, label)
        parts += p
        ver += v
    cells, ver = pd.concat(parts, ignore_index=True), pd.DataFrame(ver)
    both = {c: bool(ver[(ver.candidate == c) & ver.config.isin(["sigma_lo", "sigma_hi"])].passes.all()) for c in ORDER}
    adopted = next((c for c in ORDER if both[c]), None)
    ver["adopted"] = ver.candidate == adopted
    # the largest rain error each candidate still passes (the sweep, then sigma_lo = f 1, then sigma_hi)
    ladder = [(f"f={f:g}", f"{f:g} x the low end (per-bin sd {f * R['sigma_lo']:.2f}, scale part {f * R['scale_sd']:.2f})") for f in SWEEP] + \
             [("sigma_lo", f"the low end (per-bin sd {R['sigma_lo']}, scale part {R['scale_sd']})"), ("sigma_hi", f"the high end (per-bin sd {R['sigma_hi']}, scale part {R['scale_sd']})")]
    largest = {}
    for c in ORDER:
        ok = "none (it fails even without rain error)"
        for label, text in ladder:
            if bool(ver[(ver.candidate == c) & (ver.config == label)].passes.iloc[0]):
                ok = text
            else:
                break
        largest[c] = ok
    save(cells.round(5), P.E2_FAKE, index=False)
    save(ver.round(5), P.E2_VERDICT, index=False)
    fig = fig_errors(cells)
    save(fig, P.E2_FIG, dpi=120)
    plt.close(fig)

    a = ["# E2: which length test survives the measured rain error (fake data)\n",
         f"Written by `research/waa_ksurface/e2_candidates.py`. Fake world on the {len(le)} TRAIN link-events of {le.cml_id.nunique()} links in the {le.event_id.nunique()} events "
         f"with a gauge ratio. Bin rule {rule_name} ({rule['t_pre_min']} min within a factor {rule['q']:g}), n_short = {n_short}. Truth: F = 3·R / (R + 2) dB, P = 0.1·R / 0.2·R dB/km, "
         f"one store (5 / 15 min), measured noise. Rain error (E1): scale part per link-event sd {R['scale_sd']}, per-bin part with {100 * R['shared_share']:.0f} % shared in space, "
         f"lag-1 correlation {R['rho']}, per-bin sd {R['sigma_lo']} (low end) and {R['sigma_hi']} (high end). {FK['reps']} replicates, medians. Pass (Deviation 3): intercept within "
         f"{FK['pass']['intercept_db']} dB and slope within {FK['pass']['intercept_db'] / l_med:.4f} dB/km in all judged cells but at most {FK['pass']['cells_may_fail']}, at both ends. "
         "No real attenuation is used as a target.\n",
         "## Verdict\n",
         md_table(ver[ver.config.isin(["sigma_lo", "sigma_hi"])].assign(candidate=lambda t: t.candidate.map(NAMES))[
             ["candidate", "config", "cells_judged", "cells_inside", "I_outside", "s_outside", "largest_I_err_db", "largest_s_err_db_km", "passes"]].round(3)),
         f"\n**Adopted: {NAMES[adopted] if adopted else 'none (STOP)'}.**\n",
         "## The largest rain error each candidate still passes\n",
         md_table(pd.DataFrame([{"candidate": NAMES[c], "largest rain error passed": largest[c]} for c in ORDER])),
         "\nThe sweep (the whole rain error scaled by f; f = 1 is the low end):\n",
         md_table(ver[~ver.config.isin(["sigma_lo", "sigma_hi"])].assign(candidate=lambda t: t.candidate.map(NAMES))[
             ["candidate", "config", "sigma", "scale_sd", "cells_judged", "cells_inside", "largest_I_err_db", "largest_s_err_db_km", "passes"]].round(3)),
         "\n## Every cell at both ends of the range\n"]
    for c in ORDER:
        for cfg in ("sigma_lo", "sigma_hi"):
            g = cells[(cells.candidate == c) & (cells.config == cfg)]
            a.append(f"### {NAMES[c]}, {cfg}\n")
            a.append(md_table(g[["band", "class", "R_cell", "F_true", "I", "I_err", "P_true", "s", "s_err_abs", "links_short", "I_ok", "s_ok"]].round(3)) if len(g) else "no cell judged")
            a.append("")
    a.append(f"(`{rel(P.E2_FAKE)}`, `{rel(P.E2_VERDICT)}`, `{rel(P.E2_FIG)}`)\n")
    if adopted is None:
        a.append("## STOP\n\nNo candidate passes at both ends of the measured rain error. E4 is not run.\n")
    save("\n".join(a) + "\n", P.E2_REPORT)
    print(f"adopted: {adopted}; largest error passed: {largest}; report: {rel(P.E2_REPORT)}", flush=True)
    if adopted is None:
        sys.exit(1)


if __name__ == "__main__":
    main()
