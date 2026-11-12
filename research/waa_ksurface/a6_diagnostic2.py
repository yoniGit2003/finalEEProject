"""The one diagnostic after truth (a) failed under Deviations 5 and 6: how does the drying constant depend on the held k_wet?

    python -m research.waa_ksurface.a6_diagnostic2

Truth (a) (1 / k_wet = 5 min, 1 / k_dry = 15 min), the same fake data as in A6. All links, radar input, K0 in both regimes,
1 / k_wet HELD on a grid; each value once with the F / P extension held at the truth and once with it free. Per fit:
1 / k_dry, the grouped-CV error (dB), and for the free fits the F knots against the truth. Written in decisions.md before
it was run; nothing is changed on the basis of it. Fake data only.
"""
import numpy as np
import pandas as pd

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a6_recovery as a6, fake6, ksurf, ladder
from research.waa_ksurface import paths as P

TAUS = [2.5, 3.73, 5.0, 7.5, 9.91]                              # 3.73 and 9.91 = the gauge-input and radar-input estimates of the failed run


def main():
    S = a6.setting()
    spec = a6.A6["truths"]["a"]
    le = S["le"]
    Dk, ex = fake6.make(S["D"], le, S["n_cells"], spec, S["curves"], S["iF"], S["iP"], a6.A6["seed"] + a6.ORDER.index("a"), S["step"], S["dist"])
    prob = ksurf.prepare(Dk, ex["band_idx"], S["curves"], Dk["rain"], 0)
    E = a6.A6["extension"]
    free_F, free_P = [np.full(len(i), E["f_start_inc_db"]) for i in S["iF"]], [E["p_start_share"] * i for i in S["iP"]]
    use = np.ones(len(le))
    fold = ladder.folds(le, use)
    ev, n_ev = prob["event"], int(prob["event"].max()) + 1
    rows, knots = [], []
    for held in (True, False):
        for tau in TAUS:
            m = ksurf.Model(prob, "K0", "K0")
            x0 = m.start(S["iF"], S["iP"]) if held else m.start(free_F, free_P)
            x0[m.sl["gw"]] = np.log(1.0 / tau)
            m.hold("gw")
            if held:
                m.hold("F", "P")
            x = m.fit(x0, use)
            Ssum, N = np.zeros(n_ev), np.zeros(n_ev)
            for f in range(fold.max() + 1):
                xf = m.fit(x, use * (fold != f), maxiter=a6.A5["maxiter_fold"])
                ls, cnt = m.per_link_event(xf)
                Ssum += np.bincount(ev[fold == f], weights=ls[fold == f], minlength=n_ev)
                N += np.bincount(ev[fold == f], weights=cnt[fold == f], minlength=n_ev)
            tau_d = float(1 / m.k(x, "drying", 0, 0))
            row = {"extension": "held at the truth" if held else "free", "tau_wet_held_min": tau, "tau_dry_fit_min": tau_d, "k_dry_over_truth": 15.0 / tau_d,
                   "within_25pct": bool(abs(15.0 / tau_d - 1) <= 0.25), "cv_error_db": ladder.cv_db(Ssum, N), "train_loss": m.last["loss"]}
            if not held:
                e = a6.extension_table(S, m, x, free_F, free_P)
                f_ = e[e.quantity == "F"]
                row["F_knots_fitted_over_true_min"], row["F_knots_fitted_over_true_max"] = float(f_.fitted_over_true.min()), float(f_.fitted_over_true.max())
                knots.append(e.assign(tau_wet_held_min=tau))
            rows.append(row)
            a6.log(str(row))
    t = pd.DataFrame(rows)
    save(t.round(5), P.A6_DIAG2, index=False)
    save(pd.concat(knots, ignore_index=True).round(5), P.A6_DIAG2_KNOTS, index=False)
    a = ["# A6 diagnostic 2: the drying constant against the held k_wet\n",
         "Written by `research/waa_ksurface/a6_diagnostic2.py`. Truth (a): 1 / k_wet = 5 min, 1 / k_dry = 15 min; the fake data of the A6 run. All links, radar input, K0 in both "
         "regimes, 1 / k_wet held. The one diagnostic allowed after the failed check; nothing was changed on the basis of it. Fake data only.\n",
         md_table(t.round(4)),
         "\n`k_dry_over_truth` must lie in 0.75–1.25 to be within 25 % (1 / k_dry between 12.0 and 20.0 min). 3.73 and 9.91 min are the gauge-input and radar-input "
         f"estimates of 1 / k_wet in the failed run. F knots of the free fits: `{rel(P.A6_DIAG2_KNOTS)}`. (`{rel(P.A6_DIAG2)}`)\n"]
    save("\n".join(a) + "\n", P.A6_DIAG2_REPORT)
    a6.log(f"report: {rel(P.A6_DIAG2_REPORT)}")


if __name__ == "__main__":
    main()
