"""The one diagnostic after truth (c) failed in A6 final: why is a k_dry that triples worth only 0.002 dB to the ladder?

    python -m research.waa_ksurface.a6_diagnostic3

Truth (c) (k_dry = (1 + W) / 30 per min: 30 min at W = 0, 10 min at W = 2 dB), the same fake data as in the A6 final run.
All links, radar input, k_wet held at its gauge-input estimate, extension free. Drying K0 and drying K1 (lambda = 1e-4,
the one the CV chose) are fitted on all events, by the same steps as stage 2 of the ladder. Reported:

    (i)  the fitted K1 curve against the true k_dry(W) at the W knots, and the share of the domain weight within 25 %
         that K1 WOULD have had if the ladder had returned it;
    (ii) the in-sample error sqrt(2 x mean Huber loss) of both fits on: all 5 min values; the links below 1 km; the links
         below 1 km in the 60 min after the last wet bin; and the same three sets with the rainy part of the event (the
         bins up to and including the last wet bin) left out.

Written in decisions.md before it was run; nothing is changed on the basis of it. Fake data only.
"""
import pickle

import numpy as np
import pandas as pd

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a6_recovery as a6, fake6, ksurf, ladder
from research.waa_ksurface import paths as P

LAM = ksurf.config()["lambdas"][0]                              # 1e-4: the lambda the CV chose for drying K1 on truth (c)
AFTER_BINS = 12                                                 # 60 min of 5 min bins


def huber(r, delta):
    a = np.abs(r)
    return np.where(a <= delta, 0.5 * r * r, delta * (a - 0.5 * delta))


def main():
    S = a6.setting()
    spec, seed = a6.spec_of("c")
    D, le = S["D"], S["le"]
    Dk, ex = fake6.make(D, le, S["n_cells"], spec, S["curves"], S["iF"], S["iP"], seed, S["step"], S["dist"])
    bi = ex["band_idx"]
    E = a6.A6["extension"]
    start_F = [np.full(len(i), E["f_start_inc_db"]) for i in S["iF"]]
    start_P = [E["p_start_share"] * i for i in S["iP"]]

    # ---- stage 1, only the fit that the ladder selected (wetting K0, drying K1 as a nuisance, gauge input): gives k_wet and the extension
    rain_g, ok5 = ladder.gauge_input(Dk)
    prob_g = ksurf.prepare(Dk, bi, S["curves"], rain_g, 0)
    prob_g["w5"] = prob_g["w5"] * ok5
    m1 = ksurf.Model(prob_g, "K0", "K1", LAM)
    x1 = m1.fit(m1.start(start_F, start_P), S["gauge"])
    tau_wet = float(1 / m1.k(x1, "wetting", 0, 0))
    cached = pickle.loads((P.CACHE / "a6_recovery_final" / "c.pkl").read_bytes())
    a6.log(f"1 / k_wet (gauge input): {tau_wet:.4f} min; in the A6 final run: {cached['tau_wet']['gauge']:.4f} min")
    assert abs(tau_wet / cached["tau_wet"]["gauge"] - 1) < 1e-3, "this is not the fit of the A6 final run"
    ext = m1.parts(x1)[4:]

    # ---- stage 2 (all links, radar input, k_wet held): drying K0, then drying K1 started from it
    prob = ksurf.prepare(Dk, bi, S["curves"], Dk["rain"], 0)
    use = np.ones(len(le))
    fits, prev = {}, None
    for rung in ("K0", "K1"):
        m = ksurf.Model(prob, "K0", rung, LAM)
        x0 = m.start([e.copy() for e in ext[0]], [e.copy() for e in ext[1]]) if prev is None else m.carry(*prev)
        x0[m.sl["gw"]] = x1[m1.sl["gw"]]
        m.hold("gw", "hw")
        x = m.fit(x0, use)
        fits[rung] = prev = (m, x)
        a6.log(f"drying {rung}: train loss {m.last['loss']:.6f} ({m.last['nit']} it)")

    # ---- (i) the K1 curve against the truth, and the share it would have had
    dens = a6.state_density(S, Dk, ex)
    m, x = fits["K1"]
    wk = m.wk["drying"]
    k0_fit = float(fits["K0"][0].k(fits["K0"][1], "drying", 0, 0))
    curve = pd.DataFrame({"W_db": wk, "tau_true_min": 1 / fake6.k_true(spec, "drying", 0.0, wk), "tau_K1_min": 1 / m.k(x, "drying", 0.0, wk), "tau_K0_min": 1 / k0_fit})
    curve["K1_k_over_true"] = curve.tau_true_min / curve.tau_K1_min
    curve["K1_within_25pct"] = (curve.K1_k_over_true - 1).abs() <= a6.J["within"]
    g = dens[dens.regime == "drying"]
    wb = np.digitize(g.W, wk[1:-1])                             # the minutes of the true state between the W knots
    curve["share_of_domain_minutes_from_this_knot"] = [float(g.minutes[wb == i].sum() / g.minutes.sum()) if i < len(wk) - 1 else 0.0 for i in range(len(wk))]
    shares = {r: a6.compare(fits[r][0], fits[r][1], spec, "drying", dens)[1] for r in fits}
    save(curve.round(4), P.A6_DIAG3_CURVE, index=False)

    # ---- (ii) where the two fits differ: the in-sample error on parts of the data
    n5 = D["le_n"] // 5
    bin_le = np.repeat(np.arange(len(le)), n5)
    pos = np.arange(len(bin_le)) - D["le_boff5"][bin_le]        # the position of a 5 min bin in its window
    last = D["le_last"][bin_le]
    counted = prob["w5"] > 0
    short = S["short"][bin_le]
    after = pos > last
    tail60 = after & (pos <= last + AFTER_BINS)
    sets = [("all 5 min values", counted), ("links below 1 km", counted & short), ("links below 1 km, 60 min after the last wet bin", counted & short & tail60),
            ("all, rainy part left out", counted & after), ("links below 1 km, rainy part left out", counted & short & after),
            ("links below 1 km, 60 min after the last wet bin, rainy part left out", counted & short & tail60 & after)]
    delta = m.cfg["huber_delta_db"]
    h = {r: huber(fits[r][0].simulate(fits[r][1])[1] - prob["target"], delta) for r in fits}
    rows = []
    for name, mask in sets:
        e0, e1 = (float(np.sqrt(2 * h[r][mask].mean())) for r in ("K0", "K1"))
        rows.append({"set": name, "values": int(mask.sum()), "share_of_all_values": float(mask.sum() / counted.sum()), "error_K0_db": e0, "error_K1_db": e1,
                     "gain_db": e0 - e1, "share_of_the_total_loss_drop": float((h["K0"][mask].sum() - h["K1"][mask].sum()) / (h["K0"][counted].sum() - h["K1"][counted].sum()))})
    t = pd.DataFrame(rows)
    save(t.round(5), P.A6_DIAG3, index=False)

    cv = cached["cv"]
    cvg = float(cv[(cv["select"] == "drying") & (cv.rung == "K1") & cv.best_lambda].gain_db.iloc[0])
    a = ["# A6 diagnostic 3: why a k_dry that triples is worth 0.002 dB\n",
         "Written by `research/waa_ksurface/a6_diagnostic3.py`. Truth (c): 1 / k_wet = 5 min, k_dry = (1 + W) / 30 per min; the fake data of the A6 final run. All links, radar input, "
         f"1 / k_wet held at its gauge-input estimate ({tau_wet:.2f} min), extension free, drying K0 and drying K1 (lambda = {LAM:g}) fitted on all events. "
         "The one diagnostic allowed after the failed check; nothing was changed on the basis of it. Fake data only.\n",
         "## (i) The K1 fit against the truth\n",
         md_table(curve.round(3)),
         f"\nShare of the domain weight within 25 % of the true k_dry: K0 {100 * shares['K0']:.0f} % (the fit the ladder returned), "
         f"K1 {100 * shares['K1']:.0f} % (the fit it did not return); the rule needs 80 %. (`{rel(P.A6_DIAG3_CURVE)}`)\n",
         "## (ii) The in-sample error of the two fits on parts of the data\n",
         md_table(t.round(4)),
         f"\nError = sqrt(2 x mean Huber loss), dB; gain = K0 minus K1. \"Rainy part left out\" = only the 5 min bins after the last wet bin of the link-event, so the last row is "
         f"the third row again (by construction). The held-out gain of K1 in the A6 final run was {cvg:.4f} dB; the margin is {ladder.C['select']['margin_db']} dB. (`{rel(P.A6_DIAG3)}`)\n"]
    save("\n".join(a) + "\n", P.A6_DIAG3_REPORT)
    a6.log(f"report: {rel(P.A6_DIAG3_REPORT)}")


if __name__ == "__main__":
    main()
