"""The one diagnostic after the failed A6 check: which ingredient of the fake world makes the fit invent a k(R, W)?

    python -m research.waa_ksurface.a6_diagnostic

Truth (a) (constant k: 1 / k = 5 min wetting, 15 min drying), same seed as in A6, with the ingredients switched: the gauge
input exact or with its error, the radar input exact or with its error, the noise and the quantization on or off. For
every run: the rungs the two stages return, and what a plain K0 / K0 fit gives for 1 / k. Written in decisions.md before
it was run; nothing is changed on the basis of it. Fake data only.
"""
import numpy as np
import pandas as pd

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a6_recovery as a6, fake6, ksurf, ladder
from research.waa_ksurface import paths as P

RUNS = [("D0", True, True, False), ("D1", True, True, True), ("D2", False, True, True), ("D3", True, False, True), ("A6 run", False, False, True)]


def plain_k0(Dk, S, bi, gauge):
    """1 / k of a K0 / K0 fit (both regimes free, extension held at the truth) on one rain input."""
    if gauge:
        rain, ok5 = ladder.gauge_input(Dk)
        use = S["gauge"]
    else:
        rain, ok5, use = Dk["rain"], 1.0, np.ones(len(S["le"]))
    prob = ksurf.prepare(Dk, bi, S["curves"], rain, 0)
    prob["w5"] = prob["w5"] * ok5
    m = ksurf.Model(prob, "K0", "K0").hold("F", "P")
    x = m.fit(m.start(S["iF"], S["iP"]), use)
    return float(1 / m.k(x, "wetting", 0, 0)), float(1 / m.k(x, "drying", 0, 0)), m.last["loss"]


def main():
    S = a6.setting()
    spec = a6.A6["truths"]["a"]
    seed = a6.A6["seed"] + a6.ORDER.index("a")
    rows = []
    for name, exact_gauge, exact_radar, noise in RUNS:
        a6.log(f"=== {name} ===")
        Dk, ex = fake6.make(S["D"], S["le"], S["n_cells"], spec, S["curves"], S["iF"], S["iP"], seed, S["step"], S["dist"],
                            noise=noise, exact_gauge=exact_gauge, exact_radar=exact_radar)
        bi = ex["band_idx"]
        res = ladder.two_stage(Dk, S["le"], bi, S["curves"], S["gauge"], inc_F=S["iF"], inc_P=S["iP"], hold_ext=True, log=a6.log, tag=f"({name}) ")
        j = a6.judge(S, spec, res, a6.state_density(S, Dk, ex))[0].set_index("regime")
        gw, gd, gl = plain_k0(Dk, S, bi, True)
        rw, rd, rl = plain_k0(Dk, S, bi, False)
        cv = {s: res[s]["table"][res[s]["table"].best_lambda].set_index("rung").cv_loss for s in ("wet", "dry")}
        rows.append({"run": name, "gauge_input": "exact" if exact_gauge else "with error", "radar_input": "exact" if exact_radar else "with error",
                     "noise": "on" if noise else "off", "wetting_rung": res["wet"]["rung"], "drying_rung": res["dry"]["rung"],
                     "wet_cv_K0": cv["wet"]["K0"], "wet_cv_best": cv["wet"].min(), "dry_cv_K0": cv["dry"]["K0"], "dry_cv_best": cv["dry"].min(),
                     "share_wet_within_25pct": j.share_within_25pct["wetting"], "share_dry_within_25pct": j.share_within_25pct["drying"],
                     "K0_gauge_tau_wet": gw, "K0_gauge_tau_dry": gd, "K0_radar_tau_wet": rw, "K0_radar_tau_dry": rd})
        a6.log(str(rows[-1]))
        save(pd.DataFrame(rows).round(5), P.A6_DIAG, index=False)
    t = pd.DataFrame(rows)
    a = ["# A6 diagnostic: which ingredient makes the fit invent a k(R, W)?\n",
         "Written by `research/waa_ksurface/a6_diagnostic.py`. Truth (a): constant k, 1 / k = 5 min (wetting) and 15 min (drying); same seed as the A6 run. "
         "The one diagnostic allowed after the failed check; nothing was changed on the basis of it. Fake data only.\n",
         "## The rungs returned by the two stages\n",
         md_table(t[["run", "gauge_input", "radar_input", "noise", "wetting_rung", "drying_rung", "share_wet_within_25pct", "share_dry_within_25pct"]].round(3)),
         "\nTrue rung: K0 in both regimes. The shares are the part of the domain weight on which the returned k is within 25 % of the truth.\n",
         "## Held-out loss: K0 against the best rung\n",
         md_table(t[["run", "wet_cv_K0", "wet_cv_best", "dry_cv_K0", "dry_cv_best"]].round(6)),
         "\n## A plain K0 / K0 fit (no richer rung offered): 1 / k in minutes\n",
         md_table(t[["run", "K0_gauge_tau_wet", "K0_gauge_tau_dry", "K0_radar_tau_wet", "K0_radar_tau_dry"]].round(2)),
         "\nTruth: 5 (wetting) and 15 (drying). `K0_gauge_*` = the gauge-input set with the gauge input; `K0_radar_*` = all links with the radar input.\n",
         f"(`{rel(P.A6_DIAG)}`)\n"]
    save("\n".join(a) + "\n", P.A6_DIAG_REPORT)
    a6.log(f"report: {rel(P.A6_DIAG_REPORT)}")


if __name__ == "__main__":
    main()
