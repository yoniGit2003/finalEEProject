"""C4: what are the dry-day jumps?

Jump = a run of >= 3 h of dry slots (link state 0; gaps of up to 15 min bridged) where tl_mean exceeds the sub-link's
whole-summer dry median (M0) by >= 1 dB slot by slot, and the median excess over the run is also >= 1 dB.
Negative jumps (<= -1 dB) are found the same way as a control. For each jump: timing, weather (SMHI hourly),
hours since the last wet slot, whether the other direction of the same link jumps at the same time
(both_directions, >= 50 % overlap) and how many other train links with a path midpoint within 5 km do
(neighbours_jumping).

Run: python research/baseline_estm/validation/c4_jumps.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from common import (V, OUT, FIG, SPH, Timer, open_data, sublink_table, groups, load_panel, usable, hours_since_wet,
                    link_midpoints)

C = V["c4"]


def find_runs(cond, max_gap):
    """Runs of True (inclusive [a, b]) with gaps of up to max_gap False slots bridged."""
    idx = np.nonzero(cond)[0]
    if not len(idx):
        return []
    breaks = np.nonzero(np.diff(idx) > max_gap + 1)[0]
    starts = np.concatenate([[idx[0]], idx[breaks + 1]])
    ends = np.concatenate([idx[breaks], [idx[-1]]])
    return list(zip(starts, ends))


def jumps_of(s, x, st, m0, met, T, since):
    dry = st == 0
    exc = x - m0
    rows = []
    for sign in (1, -1):
        cond = dry & ~np.isnan(x) & (sign * exc >= C["excess_db"])
        for a, b in find_runs(cond, int(C["max_gap_min"] / 5)):
            dur = (b - a + 1) / SPH
            if dur < C["min_duration_h"]:
                continue
            seg = slice(a, b + 1)
            v = exc[seg][dry[seg] & ~np.isnan(x[seg])]
            if sign * np.median(v) < C["excess_db"]:
                continue
            m = met.iloc[a:b + 1]
            rows.append({"sublink_id": s, "cml_id": int(s.split("_")[0]), "sign": "positive" if sign > 0 else "negative",
                         "start": T[a], "end": T[b] + pd.Timedelta("5min"), "duration_h": dur,
                         "median_excess_db": float(np.median(v)),
                         "peak_excess_db": float(v.max() if sign > 0 else v.min()),
                         "start_hour_utc": T[a].hour, "start_hour_local": (T[a].hour + 2) % 24,
                         "rh_mean": m.humidity.mean(), "rh_max": m.humidity.max(), "t_mean": m.temperature.mean(),
                         "dpd_mean": m.dpd.mean(), "wind_mean": m.wind_speed.mean(),
                         "hours_since_wet": float(since[a]), "a": a, "b": b})
    return rows


def overlaps(j):
    """both_directions and neighbours_jumping for every jump (same sign, >= overlap_share of THIS jump's time)."""
    mid = j.attrs["mid"]
    j = j.copy()
    j["both_directions"] = False
    j["neighbours_jumping"] = 0
    for sign, g in j.groupby("sign"):
        a, b = g.a.to_numpy(), g.b.to_numpy() + 1
        cml, sid = g.cml_id.to_numpy(), g.sublink_id.to_numpy()
        xy = mid.reindex(cml).to_numpy()
        for n, i in enumerate(g.index):
            ov = np.clip(np.minimum(b, b[n]) - np.maximum(a, a[n]), 0, None) / (b[n] - a[n])
            hit = ov >= C["overlap_share"]
            j.loc[i, "both_directions"] = bool((hit & (cml == cml[n]) & (sid != sid[n])).any())
            dist = np.hypot(*(xy - xy[n]).T) / 1000
            j.loc[i, "neighbours_jumping"] = len(set(cml[hit & (cml != cml[n]) & (dist <= C["neighbour_km"])]))
    return j


def summarise(j, label):
    out = []
    for sign, g in j.groupby("sign"):
        out.append({"group": label, "sign": sign, "count": len(g), "sublinks_with_jump": g.sublink_id.nunique(),
                    "duration_h_median": g.duration_h.median(), "median_excess_median": g.median_excess_db.median(),
                    "share_both_directions": g.both_directions.mean(),
                    "share_neighbour_jumping": (g.neighbours_jumping >= 1).mean(),
                    f"share_rh_mean_ge_{C['high_rh_pct']}": (g.rh_mean >= C["high_rh_pct"]).mean(),
                    "dpd_mean_median": g.dpd_mean.median(), "wind_mean_median": g.wind_mean.median(),
                    "start_hour_utc_mode": int(g.start_hour_utc.mode().iloc[0]),
                    "share_start_00_05_utc": g.start_hour_utc.between(0, 5).mean(),
                    "hours_since_wet_median": g.hours_since_wet.replace(np.inf, np.nan).median(),
                    "share_more_than_24h_after_rain": (g.hours_since_wet > 24).mean()})
    return pd.DataFrame(out)


def plots(j, dry_w, met):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4))
    for ax, sign in zip(axes, ("positive", "negative")):
        g = j[j.sign == sign]
        ax.bar(range(24), g.start_hour_utc.value_counts().reindex(range(24), fill_value=0), color="C3" if sign == "positive" else "C0")
        ax.set_xticks(range(0, 24, 3), [f"{h}\n({(h + 2) % 24})" for h in range(0, 24, 3)])
        ax.set_xlabel("start hour UTC (local)")
        ax.set_ylabel("jumps")
        ax.set_title(f"{sign} jumps (n = {len(g)})")
    fig.suptitle("C4: start hour of dry-day jumps, main sub-links")
    fig.tight_layout()
    fig.savefig(FIG / "c4_start_hour.png", dpi=110)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13, 4))
    for ax, var, bins, lab in ((axes[0], "humidity", np.arange(30, 102, 2), "relative humidity [%]"),
                               (axes[1], "dpd", np.arange(0, 20.5, 0.5), "dew-point depression T - Td [deg C]")):
        v = met[var].to_numpy()
        ok = ~np.isnan(v)
        ax.hist(v[ok], bins=bins, weights=dry_w["all"][ok], density=True, histtype="step", lw=2, color="k",
                label="all dry slots")
        for sign, col in (("positive", "C3"), ("negative", "C0")):
            ax.hist(v[ok], bins=bins, weights=dry_w[sign][ok], density=True, histtype="step", lw=1.5, color=col,
                    label=f"{sign} jumps")
        ax.set_xlabel(lab)
        ax.set_ylabel("density (sub-link slots)")
        ax.legend(fontsize=8)
    fig.suptitle("C4: weather during dry-day jumps vs all dry slots, main sub-links")
    fig.tight_layout()
    fig.savefig(FIG / "c4_humidity_dewpoint.png", dpi=110)
    plt.close(fig)


def main():
    with Timer("C4"):
        d = open_data()
        sl = sublink_table(d)
        tl, state, cert, met = load_panel(d, sl)
        T = tl.index
        rows = []
        for s in sl.sublink_id:                                   # all groups: needed for neighbour counts
            c = int(s.split("_")[0])
            x = tl[s].to_numpy()
            ok = usable(state[c], V["skip_h"])
            good = ok & ~np.isnan(x)
            if not good.any():
                continue
            rows += jumps_of(s, x, state[c].to_numpy(), float(np.median(x[good])), met, T,
                             hours_since_wet(state[c]))
        j = pd.DataFrame(rows)
        j.attrs["mid"] = link_midpoints(d)
        j = overlaps(j).merge(sl[["sublink_id", "group", "is_degraded", "is_toggling", "length_class",
                                  "frequency_band", "tsl_behaviour"]], on="sublink_id", how="left")
        j.drop(columns=["a", "b"]).to_csv(OUT / "c4_jumps.csv", index=False)
        G = groups(sl)
        summ = pd.concat([summarise(j[j.sublink_id.isin(ids)], g) for g, ids in G.items()])
        summ.to_csv(OUT / "c4_summary.csv", index=False)
        # slot weights for the weather distributions: how many main sub-links are dry / in a jump at each time
        main_ids = set(G["main"])
        jm = j[j.sublink_id.isin(main_ids)]
        w = {"all": np.zeros(len(T)), "positive": np.zeros(len(T)), "negative": np.zeros(len(T))}
        for s in G["main"]:
            c = int(s.split("_")[0])
            w["all"] += (state[c].to_numpy() == 0) & ~np.isnan(tl[s].to_numpy())
        for r in jm.itertuples():
            i0, i1 = T.searchsorted(r.start), T.searchsorted(r.end)
            w[r.sign][i0:i1] += 1
        plots(jm, w, met)
        print(summ.round(3).T.to_string())


if __name__ == "__main__":
    main()
