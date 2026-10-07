"""C5: the flagged cases of the first sample run (results/baselines.csv, flags excursion_before / excursion_after).

For each case: before / after of the sample run, Delta_k (C2 definition: level in [end + k, end + k + 1 h) minus
the M2b before) for this sub-link and the other direction, humidity and dew-point depression after the rain,
whether a C4 jump overlaps the case, and a one-line hint:
  drying-like : both directions high after the rain (Delta at k = 1 h >= high_delta_db) and Delta shrinks with k
                (the last Delta at k >= 6 h below shrink_share x Delta at k = 1 h);
  shift-like  : one direction only, or Delta does not shrink;
  unclear     : not enough after-rain data.
One figure per case: tl_mean of both sub-links with the before level, radar along the path, humidity.

Run after c4_jumps.py: python research/baseline_estm/validation/c5_flagged.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from common import (V, OUT, FIG, SPH, HERE, Timer, open_data, sublink_table, load_panel, usable, window, stat)
from c2_after_window import deltas

C = V["c5"]
K = V["c2"]["k_h"]


def case_numbers(s, x, ok, wet, T, i0, i1):
    min_n, w2 = int(round(V["variants"]["min_dry_h"] * SPH)), int(V["variants"]["M2_window_h"] * SPH)
    before = stat(window(x, ok, i0 - w2, i0), "median", min_n)
    w_pos = np.nonzero(wet[i0:i1])[0]
    end_pos = i0 + w_pos[-1] + 1 if len(w_pos) else i1
    return before, end_pos, deltas(x, ok, before, end_pos)


def hint(dk_self, dk_other):
    d1s, d1o = dk_self.get("d1", np.nan), dk_other.get("d1", np.nan)
    late = [dk_self[f"d{k}"] for k in K if k >= 6 and not np.isnan(dk_self[f"d{k}"])]
    if np.isnan(d1s) or not late:
        return "unclear (not enough dry data after the rain)"
    high_both = d1s >= C["high_delta_db"] and (not np.isnan(d1o)) and d1o >= C["high_delta_db"]
    shrinks = late[-1] < C["shrink_share"] * d1s
    if high_both and shrinks:
        return "drying-like (both directions high after rain, shrinking with k)"
    why = "not high in both directions" if not high_both else "Delta does not shrink"
    return f"shift-like ({why})"


def main():
    with Timer("C5"):
        d = open_data()
        sl = sublink_table(d)
        tl, state, cert, met = load_panel(d, sl)
        allsl = d.sublinks()
        T = tl.index
        res = pd.read_csv(HERE.parent / "results" / "baselines.csv", parse_dates=["event_start", "event_end"])
        cases = res[res["flags"].fillna("").str.contains("excursion")]
        jumps = pd.read_csv(OUT / "c4_jumps.csv", parse_dates=["start", "end"])
        rows = []
        for r in cases.itertuples():
            c = int(r.sublink_id.split("_")[0])
            other = [s for s in allsl[allsl.cml_id == c].sublink_id if s != r.sublink_id][0]
            i0, i1 = T.searchsorted(r.event_start), T.searchsorted(r.event_end)
            wet = state[c].to_numpy() == 1
            ok = usable(state[c], V["skip_h"])
            out = {"event_id": r.event_id, "sublink_id": r.sublink_id, "flags": r.flags,
                   "before_sample_run": r.before_db, "after_sample_run": r.after_db}
            dk = {}
            for tag, s in (("self", r.sublink_id), ("other", other)):
                x = tl[s].to_numpy() if s in tl else d.cml(s, T[0], T[-1] + pd.Timedelta("5min"), "5min",
                                                           variables=["tl_mean"]).tl_mean.to_series().to_numpy()
                before, end_pos, dd = case_numbers(s, x, ok, wet, T, i0, i1)
                dk[tag] = dd
                out[f"before_M2b_{tag}"] = before
                out.update({f"{tag}_d{k}": dd[f"d{k}"] for k in K})
            out["last_wet_end"] = T[min(end_pos, len(T) - 1)]
            m = met.iloc[end_pos:end_pos + int(C["after_met_h"] * SPH)]
            out.update({"rh_after_mean": m.humidity.mean(), "dpd_after_mean": m.dpd.mean()})
            span0, span1 = r.event_start - pd.Timedelta(hours=6), r.event_end + pd.Timedelta(hours=24)
            jj = jumps[(jumps.sublink_id == r.sublink_id) & (jumps.start < span1) & (jumps.end > span0)]
            out["c4_jump_overlap"] = ";".join(sorted(set(jj.sign))) if len(jj) else "none"
            out["hint"] = hint(dk["self"], dk["other"])
            rows.append(out)
            plot_case(d, r, other, c, tl, met, T, out)
        t = pd.DataFrame(rows)
        t.to_csv(OUT / "c5_flagged_cases.csv", index=False)
        print(t[["event_id", "sublink_id", "flags", "self_d1", "self_d6", "other_d1", "rh_after_mean",
                 "c4_jump_overlap", "hint"]].round(2).to_string())


def plot_case(d, r, other, c, tl, met, T, out):
    t0 = r.event_start - pd.Timedelta(hours=C["plot_before_h"])
    t1 = r.event_end + pd.Timedelta(hours=C["plot_after_h"])
    t0, t1 = max(t0, T[0]), min(t1, T[-1])
    sel = (T >= t0) & (T < t1)
    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(13, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1.2, 1.2]})
    for s, col, tag in ((r.sublink_id, "C0", "self"), (other, "C1", "other")):
        if s in tl:
            a1.plot(T[sel], tl[s][sel], color=col, lw=0.8, label=f"{s} tl_mean")
        a1.axhline(out[f"before_M2b_{tag}"], color=col, ls="--", lw=1.2, label=f"{s} before (M2b)")
    for a in (a1, a2, a3):
        a.axvspan(r.event_start, r.event_end, color="#4c8ed9", alpha=0.08)
        a.axvline(out["last_wet_end"], color="k", ls=":", lw=1)
    a1.set_ylabel("TL [dB]")
    a1.legend(fontsize=7, loc="upper left")
    a1.set_title(f"{r.sublink_id}  {r.event_id}  [{r.flags}]  ->  {out['hint']}\n"
                 f"blue band = event, dotted = end of last wet slot of the link")
    rad = d.radar_along(c, t0, t1 + pd.Timedelta("5min")).R_mean.sel(cml_id=c).to_series()
    a2.plot(rad.index, rad.values, "k", lw=0.8)
    a2.set_ylabel("radar R [mm/h]")
    a3.plot(T[sel], met.humidity[sel], "m", lw=0.8)
    a3.set_ylabel("RH [%]")
    a3.set_xlabel("time (UTC)")
    fig.tight_layout()
    fig.savefig(FIG / f"c5_{r.sublink_id}_{r.event_id}.png", dpi=100)
    plt.close(fig)


if __name__ == "__main__":
    main()
