"""C2: on real events, does the "after" window carry a drying signal?

For every significant train event x main sub-link: before = M2b (median of usable dry slots in the 6 h before the
event start). end = end of the last wet slot of the link inside the event. For each k: level_k = median of usable
dry slots in [end + k, end + k + 1 h) (needs >= 30 min), Delta_k = level_k - before.
Control: the same Delta_k on the C1 fake events (end = t0 + D), weighted so that their durations match the real
events (each real case is matched to the nearest fake duration D). A real curve above the control at small k that
meets it at large k means the after window still contains the drying.
Also: how often the current method (M3) draws the line B on real events, and Delta_k for those cases.

Run after c1_fake_events.py: python research/baseline_estm/validation/c2_after_window.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from common import (V, OUT, FIG, SPH, Timer, open_data, sublink_table, groups, load_panel, usable, window, stat,
                    m3_baseline, wquantile)
from research.common.outputs import save_output

C = V["c2"]


def deltas(x, ok, before, end_pos):
    """Delta_k for every k: median of usable dry slots in [end + k, end + k + 1 h) minus before."""
    w, min_n = int(C["level_window_h"] * SPH), int(round(C["min_dry_level_h"] * SPH))
    out = {}
    for k in C["k_h"]:
        a = end_pos + int(k * SPH)
        out[f"d{k}"] = stat(window(x, ok, a, a + w), "median", min_n) - before
    return out


def real_cases(d, sl, ids, tl, state):
    ev = d.events(significant=True, split="train")
    T = tl.index
    vv = V["variants"]
    min_n, w2 = int(round(vv["min_dry_h"] * SPH)), int(vv["M2_window_h"] * SPH)
    rows = []
    for s in ids:
        c = int(s.split("_")[0])
        x = tl[s].to_numpy()
        ok = usable(state[c], V["skip_h"])
        wet = state[c].to_numpy() == 1
        good = ok & ~np.isnan(x)
        usual = float(np.median(x[good])) if good.any() else np.nan
        for e in ev.itertuples():
            i0, i1 = T.searchsorted(pd.Timestamp(e.start)), T.searchsorted(pd.Timestamp(e.end))
            r = {"sublink_id": s, "event_id": e.event_id, "event_start": e.start, "event_end": e.end}
            w_pos = np.nonzero(wet[i0:i1])[0]
            if not len(w_pos):
                r["case_status"] = "no_wet_slot"
                rows.append(r)
                continue
            end_pos = i0 + w_pos[-1] + 1                                   # end of the last wet slot
            before = stat(window(x, ok, i0 - w2, i0), "median", min_n)
            _, info = m3_baseline(x, ok, usual, i0, i1)
            r.update({"case_status": "ok" if not np.isnan(before) else "no_before", "last_wet_end": T[min(end_pos, len(T) - 1)],
                      "duration_h": (end_pos - i0) / SPH, "before": before, "m3_method": info["m3_method"],
                      "m3_flag": info["m3_flag"], "m3_after_minus_before": info["after"] - info["before"],
                      **deltas(x, ok, before, end_pos)})
            rows.append(r)
    return pd.DataFrame(rows).merge(sl, on="sublink_id", how="left")


def control_cases(sl, ids, tl, state):
    ev = pd.read_csv(OUT / "c1_fake_events.csv.gz", parse_dates=["t0"])
    ev = ev[ev.sublink_id.isin(ids)]
    rows = []
    for s, g in ev.groupby("sublink_id"):
        c = int(s.split("_")[0])
        x = tl[s].to_numpy()
        ok = usable(state[c], V["skip_h"])
        for e in g.itertuples():
            end_pos = e.i0 + int(e.D_h * SPH)
            rows.append({"sublink_id": s, "D_h": e.D_h, "before": e.base_M2b,
                         **deltas(x, ok, e.base_M2b, end_pos)})
    return pd.DataFrame(rows).merge(sl, on="sublink_id", how="left")


def curve(real, ctrl, label=""):
    """Median / p10 / p90 / count of Delta_k against k, real vs duration-matched control."""
    Ds = np.array(sorted(ctrl.D_h.unique()))
    real = real.assign(D_match=Ds[np.abs(real.duration_h.to_numpy()[:, None] - Ds[None, :]).argmin(1)])
    n_real_D = real.D_match.value_counts()
    n_ctrl_D = ctrl.D_h.value_counts()
    wts = ctrl.D_h.map(lambda D: n_real_D.get(D, 0) / n_ctrl_D[D]).to_numpy()
    rows = []
    for k in C["k_h"]:
        x = real[f"d{k}"].dropna()
        y = ctrl[f"d{k}"].to_numpy()
        rows.append({"group": label, "k_h": k, "real_n": len(x), "real_median": x.median(),
                     "real_p10": x.quantile(0.1), "real_p90": x.quantile(0.9),
                     "ctrl_n": int((~np.isnan(y) & (wts > 0)).sum()),
                     "ctrl_median": wquantile(y, wts, 0.5), "ctrl_p10": wquantile(y, wts, 0.1),
                     "ctrl_p90": wquantile(y, wts, 0.9)})
    return pd.DataFrame(rows)


def plot(cv, m3k):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
    ax = axes[0]
    k = cv.k_h
    ax.fill_between(k, cv.real_p10, cv.real_p90, color="C3", alpha=0.15)
    ax.plot(k, cv.real_median, "o-", color="C3", label="real events: median (band p10-p90)")
    ax.fill_between(k, cv.ctrl_p10, cv.ctrl_p90, color="C0", alpha=0.15)
    ax.plot(k, cv.ctrl_median, "s--", color="C0", label="fake control: median (band p10-p90)")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xscale("log")
    ax.set_xticks(C["k_h"], [str(v) for v in C["k_h"]])
    ax.set_xlabel("k = hours after the end of the last wet slot")
    ax.set_ylabel("Delta_k = level_k - before [dB]")
    ax.set_title("main sub-links, significant train events")
    ax.legend(fontsize=8)
    ax = axes[1]
    if len(m3k):
        ax.fill_between(m3k.k_h, m3k.p10, m3k.p90, color="C1", alpha=0.15)
        ax.plot(m3k.k_h, m3k["median"], "o-", color="C1", label="cases where M3 chose the line (B)")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xscale("log")
    ax.set_xticks(C["k_h"], [str(v) for v in C["k_h"]])
    ax.set_xlabel("k [h]")
    ax.set_title("M3 line (B) cases: does after - before shrink with k?")
    ax.legend(fontsize=8)
    fig.suptitle("C2: does the after window carry a drying signal?")
    fig.tight_layout()
    save_output(fig, FIG / "c2_delta_vs_k.png", "diag", "C2 Delta_k after real rain vs dry control", dpi=110)
    plt.close(fig)


def main():
    with Timer("C2"):
        d = open_data()
        sl = sublink_table(d)
        tl, state, cert, met = load_panel(d, sl)
        ids = groups(sl)["main"]
        real = real_cases(d, sl, ids, tl, state)
        save_output(real, OUT / "c2_real_cases.csv", "diag", "C2 real (event, sub-link) cases", index=False)
        ok = real[real.case_status == "ok"]
        ctrl = control_cases(sl, ids, tl, state)
        cv = curve(ok, ctrl, "all main")
        cls = []
        for k, g in ok.groupby(["length_class", "frequency_band"]):
            cc = ctrl[(ctrl.length_class == k[0]) & (ctrl.frequency_band == k[1])]
            cls.append(curve(g, cc, f"{k[0]} / {k[1]}"))
        save_output(pd.concat([cv] + cls), OUT / "c2_delta_vs_k.csv", "diag", "C2 Delta_k real vs control", index=False)
        status = real.case_status.value_counts().rename_axis("case_status").reset_index(name="cases")
        save_output(status, OUT / "c2_case_status.csv", "diag", "C2 case counts", index=False)
        # M3 behaviour on real events
        m3 = ok.m3_method.value_counts(normalize=True).rename_axis("m3_method").reset_index(name="share")
        save_output(m3, OUT / "c2_m3_method_share.csv", "diag", "C2 share of M3 methods A / B", index=False)
        bb = ok[ok.m3_method == "B"]
        m3k = pd.DataFrame([{"k_h": k, "n": bb[f"d{k}"].notna().sum(), "median": bb[f"d{k}"].median(),
                             "p10": bb[f"d{k}"].quantile(0.1), "p90": bb[f"d{k}"].quantile(0.9)} for k in C["k_h"]])
        save_output(m3k, OUT / "c2_m3_line_cases_delta_vs_k.csv", "diag",
                    "C2 Delta_k of the line (B) cases", index=False)
        plot(cv, m3k)
        print(status.to_string(), "\n", cv.round(3).to_string(), "\n", m3.round(3).to_string(), "\n",
              m3k.round(3).to_string())


if __name__ == "__main__":
    main()
