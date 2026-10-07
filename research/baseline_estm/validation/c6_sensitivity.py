"""C6: sensitivity of M2b on the C1 fake events (D in {1, 6, 24} h): window length, skip S, mean vs median.

For every setting the baseline at t0 is recomputed on the same fake events; e(t) = tl_mean(t) - baseline over the
fake event. Returns median e, p10, p90, median |e| and the missing share per setting and D.

Run after c1_fake_events.py: python research/baseline_estm/validation/c6_sensitivity.py
"""
import numpy as np
import pandas as pd
from common import V, OUT, SPH, Timer, open_data, sublink_table, groups, load_panel, usable, window, stat, summary

C = V["c6"]


def main():
    with Timer("C6"):
        d = open_data()
        sl = sublink_table(d)
        tl, state, cert, met = load_panel(d, sl)
        ids = groups(sl)["main"]
        ev = pd.read_csv(OUT / "c1_fake_events.csv.gz")
        ev = ev[ev.sublink_id.isin(ids) & ev.D_h.isin(C["durations_h"])]
        min_n = int(round(V["variants"]["min_dry_h"] * SPH))
        settings = [(w, s, h) for w in C["window_h"] for s in C["skip_h"] for h in C["stats"]]
        e = {(st, D): [] for st in settings for D in C["durations_h"]}
        miss = {(st, D): [] for st in settings for D in C["durations_h"]}
        for sid, g in ev.groupby("sublink_id"):
            c = int(sid.split("_")[0])
            x = tl[sid].to_numpy()
            oks = {S: usable(state[c], S) for S in C["skip_h"]}
            for r in g.itertuples():
                n = int(r.D_h * SPH)
                seg = x[r.i0:r.i0 + n]
                for (w, S, how) in settings:
                    b = stat(window(x, oks[S], r.i0 - int(w * SPH), r.i0), how, min_n)
                    e[((w, S, how), r.D_h)].append(seg - b)
                    miss[((w, S, how), r.D_h)].append(np.isnan(b))
        rows = []
        for (st, D), parts in e.items():
            w, S, how = st
            rows.append({"window_h": w, "skip_h": S, "stat": how, "D_h": D, "n_events": len(parts),
                         "missing_share": float(np.mean(miss[(st, D)])) if parts else np.nan,
                         **summary(np.concatenate(parts) if parts else [])})
        t = pd.DataFrame(rows).sort_values(["D_h", "window_h", "skip_h", "stat"])
        t.to_csv(OUT / "c6_sensitivity.csv", index=False)
        print(t.round(3).to_string())


if __name__ == "__main__":
    main()
