"""TEST comparison figures (decisions.md, "TEST comparison figures", logged before they were made). A VISUALIZATION ONLY:
no fitting, no metric, no decision follows from it.

    python -m research.waa_ksurface.t_figures

Measured A - c and the predictions: ONLY the saved outputs of the TEST run (research/.cache/waa_ksurface/t_test/
predictions_test.pkl). No model is rerun. The radar rain input R = radar path mean / g_e of the two plotted link-events is
reloaded from the TEST rain input (data-layer store radar_along_links_5min, R_mean; g_e of test_s8_g_e.csv); no attenuation
is read from it. Two figures: rule 1 (typical short link) and rule 2 (typical long link), decisions.md "TEST comparison
figures replaced". Not saved by the TEST run, so not shown: the WAA part of each model.
"""
import pickle

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import data, paths as P
from research.wet_antenna import paths as WP

SRC = P.CACHE / "t_test" / "predictions_test.pkl"
MODELS = ["ours", "Part 1 equation", "Schleiss b", "Pastorek b", "no memory", "no WAA"]
LABEL = {"ours": "ours (frozen equation)", "Part 1 equation": "Part 1 equation", "Schleiss b": "Schleiss (b)", "Pastorek b": "Pastorek (b)",
         "no memory": "no memory (recipe)", "no WAA": "no WAA (recipe)"}
COLOR = {"ours": "#0072B2", "Part 1 equation": "#D55E00", "Schleiss b": "#009E73", "Pastorek b": "#CC79A7",   # Okabe-Ito, validated
         "no memory": "#E69F00", "no WAA": "#56B4E9"}
GREY, SPAN = "#8c8c8c", "#e3e3e3"
DRY_MIN_BINS = 6                                               # 30 min after the last wet bin
NOT_SAVED = "Not saved by the TEST run, so not shown: the WAA part of each model (only the total prediction was saved)."
WET_MMH = 0.1                                                  # shading: radar rain input above the wet threshold of the data layer
EXPECTED = {"rule1_short": ("E20150828_0630", "10031_1"), "rule2_long": ("E20150628_2200", "10210_2")}   # the figures made before


def rain_input(row):
    """R = radar path mean / g_e of one link-event over its window's 5 min bins (the model input; missing = 0)."""
    ra = data.open_release().store("radar_along_links_5min").R_mean
    T5 = pd.DatetimeIndex(ra.time.values)
    i0 = T5.get_loc(pd.Timestamp(row.window_start))
    r = ra.sel(cml_id=int(row.cml_id)).isel(time=slice(i0, i0 + int(row.n_bins))).load().to_numpy().astype(float)
    g = pd.read_csv(WP.S8_G_E).set_index("event_id")
    ge = float(g.g_path[row.event_id]) if (row.event_id in g.index and bool(g.path_used[row.event_id])) else 1.0
    return np.nan_to_num(r) / ge, ge, int(np.isnan(r).sum())


def load():
    """Per set: le, target (A - c), counted, scored, predictions, bin offsets."""
    raw = pickle.loads(SRC.read_bytes())
    out = {}
    for name, v in raw.items():
        le = v["le"].reset_index(drop=True)
        n = le.n_bins.to_numpy(int)
        assert n.sum() == len(v["y5"]), "the table and the saved values do not line up"
        off = np.r_[0, np.cumsum(n)]
        bin_le = np.repeat(np.arange(len(le)), n)
        target = v["y5"] - v["le_c"][bin_le]
        counted = (v["w5"] > 0) & np.isfinite(target)
        scored = counted & np.logical_and.reduce([np.isfinite(p) for p in v["preds"].values()])
        out[name] = {"le": le, "off": off, "target": target, "counted": counted, "scored": scored, "preds": v["preds"]}
    return out


def rmse_per_le(s, model):
    """RMSE of a model on each link-event over its scored values (NaN if none)."""
    le, off = s["le"], s["off"]
    r = np.full(len(le), np.nan)
    for e in range(len(le)):
        sl = slice(off[e], off[e + 1])
        m = s["scored"][sl]
        if m.any():
            r[e] = np.sqrt(np.mean((s["preds"][model][sl][m] - s["target"][sl][m]) ** 2))
    return r


def select(s, rule, Lmask, rm):
    """The link-event whose RMSE (ours) is closest to the median among the candidates (tie: the first in table order)."""
    le = s["le"]
    dry_after = le.n_bins.to_numpy() - le.last_bin.to_numpy() - 1
    cand = np.nonzero(Lmask & (dry_after >= DRY_MIN_BINS) & np.isfinite(rm))[0]
    med = float(np.median(rm[cand]))
    pick = int(cand[np.argmin(np.abs(rm[cand] - med))])
    rank = int(np.sum(rm[cand] < rm[pick])) + 1
    row = {"rule": rule, "candidates": len(cand), "median RMSE ours (dB)": med, "event": le.event_id[pick], "sublink": le.sublink_id[pick],
           "L (km)": float(le.L_km[pick]), "band": le.band_name[pick], "RMSE ours (dB)": float(rm[pick]), "rank (ascending)": f"{rank} of {len(cand)}",
           "min after last wet bin": int(dry_after[pick] * 5)}
    print(f"   ({rule}) {len(cand)} candidates; median RMSE ours {med:.4f} dB; selected {row['event']} {row['sublink']} "
          f"(L {row['L (km)']:.2f} km, {row['band']}): RMSE {row['RMSE ours (dB)']:.4f} dB, rank {row['rank (ascending)']}", flush=True)
    return pick, row


def ylims(arrays):
    v = np.concatenate([a[np.isfinite(a)] for a in arrays])
    lo, hi = float(np.min(v)), float(np.max(v))
    pad = 0.06 * (hi - lo if hi > lo else 1.0)
    return lo - pad, hi + pad


def style(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.grid(axis="y", color="#ececec", lw=0.6)
    ax.tick_params(labelsize=8)


def event_figure(s, e, rule, known, R):
    le = s["le"]
    sl = slice(s["off"][e], s["off"][e + 1])
    n = int(le.n_bins[e])
    t = np.arange(n) * 5 + 2.5
    meas = np.where(s["counted"][sl], s["target"][sl], np.nan)
    pred = {m: np.where(np.isfinite(s["preds"][m][sl]), s["preds"][m][sl], np.nan) for m in MODELS}
    sc = s["scored"][sl]
    yl = ylims([meas] + list(pred.values()))
    wet = R > WET_MMH
    fig, axes = plt.subplots(7, 1, figsize=(8.5, 13), sharex=True)
    for i, ax in enumerate(axes):
        if i > 0:                                               # panel 1: the shading goes on the rain axis, behind the bars
            for k in np.nonzero(wet)[0]:
                ax.axvspan(k * 5, k * 5 + 5, color=SPAN, lw=0, zorder=0)
        style(ax)
        ax.set_ylim(*yl)
        ax.set_ylabel("dB", fontsize=8)
    ax2 = axes[0].twinx()                                       # requested: the rain input on a secondary axis, panel 1 only
    for k in np.nonzero(wet)[0]:
        ax2.axvspan(k * 5, k * 5 + 5, color=SPAN, lw=0, zorder=0)
    ax2.bar(t, R, width=4.6, color="#6baed6", alpha=0.8, lw=0, zorder=1, label="radar rain R = radar / g_e (mm/h)")
    ax2.set_ylabel("R (mm/h)", fontsize=8, color="#4a7fa8")
    ax2.tick_params(labelsize=8, colors="#4a7fa8")
    ax2.spines["top"].set_visible(False)
    ax2.set_ylim(0, max(1.0, 1.15 * float(np.max(R))))
    axes[0].set_zorder(ax2.get_zorder() + 1)
    axes[0].patch.set_visible(False)
    axes[0].plot(t, meas, color="black", lw=1.6, label="measured A - c (5 min)")
    axes[0].set_title("measured A - c, with the radar rain input", fontsize=9, loc="left")
    h1, l1 = axes[0].get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    axes[0].legend(h1 + h2, l1 + l2, fontsize=7, loc="upper right", frameon=False)
    for ax, m in zip(axes[1:], MODELS):
        ax.plot(t, meas, color=GREY, lw=1.0, ls="--", label="measured")
        ax.plot(t, pred[m], color=COLOR[m], lw=1.8, label=f"{LABEL[m]}: total (WAA + path)")
        r = np.sqrt(np.mean((pred[m][sc] - meas[sc]) ** 2)) if sc.any() else np.nan
        ax.set_title(f"{LABEL[m]}  -  RMSE {r:.3f} dB on this link-event", fontsize=9, loc="left")
    axes[1].legend(fontsize=7, loc="upper right", frameon=False)
    axes[-1].set_xlabel(f"minutes from the start of the link-event window (grey: radar rain input R > {WET_MMH} mm/h)", fontsize=9)
    fig.suptitle(f"{le.event_id[e]}   sub-link {le.sublink_id[e]}   L = {le.L_km[e]:.2f} km   {le.band_name[e]}   "
                 f"{'known link' if known else 'NEW link (never fitted; Part 1 models with a_link = 1)'}\n"
                 f"rule {rule}   (TEST; visualization only)", fontsize=10)
    fig.text(0.01, 0.005, NOT_SAVED, fontsize=7, color="#555555")
    fig.tight_layout(rect=(0, 0.015, 1, 0.97))
    stem = P.FIGURES / f"waa_compare_{rule}_{le.event_id[e]}_{le.sublink_id[e]}"
    save(fig, stem.with_suffix(".png"), dpi=300)
    save(fig, stem.with_suffix(".pdf"))
    plt.close(fig)
    return stem


def main():
    sets = load()
    rows = []
    kn = sets["known"]
    rk = rmse_per_le(kn, "ours")
    Lk = kn["le"].L_km.to_numpy()
    for rule, Lmask in (("rule1_short", Lk < 1.0), ("rule2_long", Lk > 3.0)):
        e, row = select(kn, rule, Lmask, rk)
        assert (row["event"], row["sublink"]) == EXPECTED[rule], "the rule no longer picks the link-event of the first figures"
        R, ge, n_nan = rain_input(kn["le"].iloc[e])
        print(f"   ({rule}) rain input: g_e {ge:.4f}; peak R {R.max():.2f} mm/h; {int((R > WET_MMH).sum())} bins above {WET_MMH} mm/h; "
              f"{n_nan} missing radar bins (0 in the input)", flush=True)
        stem = event_figure(kn, e, rule, True, R)
        row.update({"g_e": ge, "peak R (mm/h)": float(R.max()), "figure": rel(stem.with_suffix(".png"))})
        rows.append(row)
    sel = pd.DataFrame(rows)
    save(sel.round(4), P.TF_SELECTION, index=False)
    L = ["# TEST comparison figures (a visualization only)\n",
         "Written by `research/waa_ksurface/t_figures.py` (decisions.md, \"TEST comparison figures replaced\"). Measured A - c and predictions: only the "
         "saved outputs of the TEST run; no model rerun. Radar rain input R = radar / g_e of the two link-events reloaded from the TEST rain input (no "
         "attenuation read). No fitting, no metric, no decision follows.\n",
         "## The two link-events (the rules of the first figures, unchanged)\n", md_table(sel.round(4)), "",
         "## Not plotted\n", NOT_SAVED + "\n"]
    save("\n".join(L) + "\n", P.TF_REPORT)
    print(f"report: {rel(P.TF_REPORT)}", flush=True)


if __name__ == "__main__":
    main()
