"""C3: daily cycle of TL in dry weather.

For each main sub-link: anomaly = usable dry slot tl_mean minus the median of usable dry slots in the centred 24 h
around it. Median anomaly per hour of day (UTC; local summer time = UTC + 2 h); amplitude = max - min of the 24
hourly medians. Also the correlation of the anomaly with temperature and with relative humidity (SMHI hourly).

Run: python research/baseline_estm/validation/c3_daily_cycle.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from common import V, OUT, FIG, SPH, Timer, open_data, sublink_table, groups, load_panel, usable

C = V["c3"]


def anomalies(x, ok):
    """Anomaly of every usable dry slot against the centred local median (NaN elsewhere)."""
    s = pd.Series(np.where(ok, x, np.nan))
    w = int(C["anomaly_window_h"] * SPH) + 1
    local = s.rolling(w, center=True, min_periods=C["anomaly_min_slots"]).median()
    return (s - local).to_numpy()


def main():
    with Timer("C3"):
        d = open_data()
        sl = sublink_table(d)
        tl, state, cert, met = load_panel(d, sl)
        ids = groups(sl)["main"]
        hour = tl.index.hour.to_numpy()
        T, RH = met.temperature.to_numpy(), met.humidity.to_numpy()
        rows, curves = [], {}
        for s in ids:
            c = int(s.split("_")[0])
            a = anomalies(tl[s].to_numpy(), usable(state[c], V["skip_h"]))
            m = ~np.isnan(a)
            hm = pd.Series(a[m]).groupby(hour[m]).median().reindex(range(24))
            curves[s] = hm.to_numpy()
            r = {"sublink_id": s, "n_slots": int(m.sum()), "amplitude_db": hm.max() - hm.min(),
                 "hour_max_utc": int(hm.idxmax()) if hm.notna().any() else -1,
                 "hour_min_utc": int(hm.idxmin()) if hm.notna().any() else -1}
            for name, v in (("temperature", T), ("humidity", RH)):
                mm = m & ~np.isnan(v)
                r[f"corr_{name}"] = np.corrcoef(a[mm], v[mm])[0, 1] if mm.sum() > 100 else np.nan
            rows.append(r)
        amp = pd.DataFrame(rows).merge(sl, on="sublink_id", how="left")
        cur = pd.DataFrame(curves, index=pd.Index(range(24), name="hour_utc")).T
        cur.columns = [f"h{h:02d}" for h in range(24)]
        amp.merge(cur, left_on="sublink_id", right_index=True).to_csv(OUT / "c3_daily_cycle_per_sublink.csv",
                                                                       index=False)
        a = amp.amplitude_db.dropna()
        summ = {"n_sublinks": len(a), "amp_p10": a.quantile(0.1), "amp_p50": a.median(), "amp_p90": a.quantile(0.9)}
        for th in C["amplitude_thresholds_db"]:
            summ[f"share_amp_gt_{th}"] = (a > th).mean()
        for name in ("temperature", "humidity"):
            x = amp[f"corr_{name}"].dropna()
            summ.update({f"corr_{name}_p10": x.quantile(0.1), f"corr_{name}_median": x.median(),
                         f"corr_{name}_p90": x.quantile(0.9)})
        summ = pd.DataFrame([summ])
        summ.to_csv(OUT / "c3_summary.csv", index=False)
        by_cls = amp.groupby(["length_class", "frequency_band"]).amplitude_db.describe(percentiles=[.1, .5, .9])
        by_cls.to_csv(OUT / "c3_amplitude_by_class.csv")
        # figure: amplitude distribution + 3 largest and 3 typical curves
        a_sorted = amp.dropna(subset=["amplitude_db"]).sort_values("amplitude_db")
        largest = a_sorted.tail(C["n_largest"]).sublink_id.tolist()
        mid = (a_sorted.amplitude_db - a.median()).abs().sort_values().index[:C["n_typical"]]
        typical = a_sorted.loc[mid].sublink_id.tolist()
        fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
        axes[0].hist(a, bins=40, color="C0")
        for th in C["amplitude_thresholds_db"]:
            axes[0].axvline(th, color="k", ls="--", lw=1)
        axes[0].set_xlabel("daily-cycle amplitude (max - min of hourly median anomaly) [dB]")
        axes[0].set_ylabel("sub-links")
        axes[0].set_title(f"{len(a)} main sub-links")
        for s, ls in [(s, "-") for s in largest] + [(s, "--") for s in typical]:
            axes[1].plot(range(24), cur.loc[s], ls, marker="o", ms=3,
                         label=f"{s} ({amp.set_index('sublink_id').amplitude_db[s]:.2f} dB)")
        axes[1].axhline(0, color="k", lw=0.8)
        axes[1].set_xticks(range(0, 24, 3), [f"{h}\n({(h + 2) % 24})" for h in range(0, 24, 3)])
        axes[1].set_xlabel("hour UTC (local)")
        axes[1].set_ylabel("median anomaly [dB]")
        axes[1].set_title("3 largest (solid) and 3 typical (dashed) amplitudes")
        axes[1].legend(fontsize=7)
        fig.suptitle("C3: daily cycle of dry-weather TL")
        fig.tight_layout()
        fig.savefig(FIG / "c3_daily_cycle.png", dpi=110)
        plt.close(fig)
        print(summ.round(3).T.to_string())


if __name__ == "__main__":
    main()
