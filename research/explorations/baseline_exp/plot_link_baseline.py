"""Exploration: TL of one link over ~12 days and around one event, with rain references and wet/dry labels.

Question: how much does the dry level (baseline) of a link drift, and what does TL look like around rain?
Run from anywhere: python research/explorations/plot_link_baseline.py  (figures are saved next to this script)
Change CML / DAYS / EVENT below to look at another example.
"""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
import yaml
from cmlpipe.api import CMLData

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
release = yaml.safe_load(open(REPO / "research" / "config.yaml"))["release"]   # pinned data version
d = CMLData.open(REPO / "data" / "release" / release)

CML = 10200                               # 6.5 km, 18.4 GHz, constant TSL on both sub-links
DAYS = ("2015-08-01", "2015-08-13")       # long view: [start, end)
EVENT = "E20150804_1830"                  # heavy event inside DAYS; zoom = event +- PAD
PAD = "12h"
STATE_COLORS = {1: "#4c8ed9", -1: "#bbbbbb"}  # link_state: wet blue, unknown grey, dry = no shading


def shade_state(ax, df):
    """Shade the background by the per-link reference state (radar + gauges, never CML)."""
    s = df.link_state.astype("float")
    for code, col in STATE_COLORS.items():
        ax.fill_between(df.index, 0, 1, where=(s == code).to_numpy(), color=col, alpha=0.25,
                        transform=ax.get_xaxis_transform(), lw=0)


def panel(df, title, path, event=None):
    """Three stacked panels on one time axis: TL per sub-link, radar along the path, nearest gauge.
    event: optional (start, end) drawn as dashed vertical lines."""
    step = df.index[1] - df.index[0]
    t_end = df.index[-1] + step                                          # the window is [first bin, last bin + step)
    title += f"\nwindow {df.index[0]:%Y-%m-%d %H:%M} to {t_end:%Y-%m-%d %H:%M} UTC"
    fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(14, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1.2, 1.2]})
    for s, c in (("s1", "C0"), ("s2", "C1")):
        a1.plot(df.index, df[f"tl_median_{s}"], c, lw=0.8, label=f"TL median 5 min, sub-link {s}")
    a1.set_ylabel("TL = TSL - RSL [dB]")
    a1.set_title(title)
    a1.legend(loc="upper left", fontsize=8)

    a2.plot(df.index, df.radar_R_mean, "k", lw=0.8, label="radar along path, R_mean (5 min)")
    miss = df.index[df.radar_R_mean.isna()]                              # missing radar bins (not zero rain!)
    a2.plot(miss, [0] * len(miss), "r|", ms=12, mew=1.5, label=f"radar missing ({len(miss)} bins)")
    a2.set_ylabel("radar R [mm/h]")
    a2.legend(loc="upper left", fontsize=8)

    g = df.attrs.get("gauge")
    if g is not None:
        a3.plot(df.index, df.gauge_rate_mm_h, "g", lw=0.8,
                label=f"gauge {g['gauge_id']} ({g['distance_to_path_m']:.0f} m from path), rate (5 min)")
        a3.legend(loc="upper left", fontsize=8)
    else:
        a3.text(0.5, 0.5, "no gauge within 1 km of the path", ha="center", transform=a3.transAxes)
    a3.set_ylabel("gauge R [mm/h]")

    for ax in (a1, a2, a3):
        shade_state(ax, df)
        if event is not None:
            for t in event:
                ax.axvline(t, color="k", ls="--", lw=1)
    if event is not None:
        a1.text(event[0], 1.0, " event start", transform=a1.get_xaxis_transform(), va="top", fontsize=8)
        a1.text(event[1], 0.95, " event end", transform=a1.get_xaxis_transform(), va="top", fontsize=8)  # lower: no overlap in long views
    # date + hour on every tick so the time is readable without guessing
    loc = mdates.AutoDateLocator(minticks=8, maxticks=16)
    a3.xaxis.set_major_locator(loc)
    a3.xaxis.set_major_formatter(mdates.DateFormatter("%d %b\n%H:%M"))
    a3.set_xlabel("time (UTC, bin start)")
    fig.text(0.01, 0.005, "shading from link_state: blue = wet, grey = unknown, white = dry", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    print("saved", path)


ev = d.events().set_index("event_id").loc[EVENT]
ev_span = (pd.Timestamp(ev.start), pd.Timestamp(ev.end))
panel(d.frame(CML, *DAYS), f"link {CML} (5 min)", HERE / f"link_{CML}_days.png", event=ev_span)
panel(d.event_window(EVENT, CML, pad=PAD), f"link {CML}, event {EVENT} +- {PAD} (5 min)",
      HERE / f"link_{CML}_event.png", event=ev_span)
