"""TEST comparison figures (decisions.md, "TEST comparison figures", logged before they were made). A VISUALIZATION ONLY:
no fitting, no metric, no decision follows from it.

    python -m research.waa_ksurface.t_figures

Reads ONLY the saved outputs of the TEST run (research/.cache/waa_ksurface/t_test/predictions_test.pkl: the link-event
table, measured y5 and offset c, the counted flag w5, every model's total 5 min prediction). No model is rerun; TEST is not
read again. Not saved by the TEST run, so not shown: the radar rain series, the WAA part of each model, the peak rain rate
(rule 3 cannot be applied).
"""
import pickle

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import paths as P

SRC = P.CACHE / "t_test" / "predictions_test.pkl"
MODELS = ["ours", "Part 1 equation", "Schleiss b", "Pastorek b", "no memory", "no WAA"]
LABEL = {"ours": "ours (frozen equation)", "Part 1 equation": "Part 1 equation", "Schleiss b": "Schleiss (b)", "Pastorek b": "Pastorek (b)",
         "no memory": "no memory (recipe)", "no WAA": "no WAA (recipe)"}
COLOR = {"ours": "#0072B2", "Part 1 equation": "#D55E00", "Schleiss b": "#009E73", "Pastorek b": "#CC79A7",   # Okabe-Ito, validated
         "no memory": "#E69F00", "no WAA": "#56B4E9"}
GREY, SPAN = "#8c8c8c", "#e3e3e3"
DRY_MIN_BINS, SUMMARY_BINS = 6, 12                             # 30 min after the last wet bin; the summary covers -60..+60 min
NOT_SAVED = "Not saved by the TEST run, so not shown: radar rain (bars and wet-minute shading), the WAA part of each model."


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


def event_figure(s, e, rule, known):
    le = s["le"]
    sl = slice(s["off"][e], s["off"][e + 1])
    n = int(le.n_bins[e])
    t = np.arange(n) * 5 + 2.5
    meas = np.where(s["counted"][sl], s["target"][sl], np.nan)
    pred = {m: np.where(np.isfinite(s["preds"][m][sl]), s["preds"][m][sl], np.nan) for m in MODELS}
    sc = s["scored"][sl]
    yl = ylims([meas] + list(pred.values()))
    a, b = int(le.first_bin[e]) * 5, (int(le.last_bin[e]) + 1) * 5
    fig, axes = plt.subplots(7, 1, figsize=(8.5, 13), sharex=True)
    for i, ax in enumerate(axes):
        ax.axvspan(a, b, color=SPAN, lw=0, zorder=0)
        style(ax)
        ax.set_ylim(*yl)
        ax.set_ylabel("dB", fontsize=8)
    axes[0].plot(t, meas, color="black", lw=1.6, label="measured A - c (5 min)")
    axes[0].set_title("measured A - c", fontsize=9, loc="left")
    axes[0].text(0.99, 0.95, "grey span: first to last wet 5 min bin of the link (data layer)", transform=axes[0].transAxes,
                 ha="right", va="top", fontsize=7, color="#555555")
    for ax, m in zip(axes[1:], MODELS):
        ax.plot(t, meas, color=GREY, lw=1.0, ls="--", label="measured")
        ax.plot(t, pred[m], color=COLOR[m], lw=1.8, label=f"{LABEL[m]}: total (WAA + path)")
        r = np.sqrt(np.mean((pred[m][sc] - meas[sc]) ** 2)) if sc.any() else np.nan
        ax.set_title(f"{LABEL[m]}  -  RMSE {r:.3f} dB on this link-event", fontsize=9, loc="left")
    axes[1].legend(fontsize=7, loc="upper right", frameon=False)
    axes[-1].set_xlabel("minutes from the start of the link-event window", fontsize=9)
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


def summary_figure(sets):
    """All TEST link-events on links < 1 km with >= 60 min after the last wet bin, aligned at the end of that bin."""
    k = np.arange(-SUMMARY_BINS, SUMMARY_BINS)
    rows_m, rows_p, n_le, links = [], {m: [] for m in MODELS}, 0, set()
    for name, s in sets.items():
        le = s["le"]
        dry_after = le.n_bins.to_numpy() - le.last_bin.to_numpy() - 1
        for e in np.nonzero((le.L_km.to_numpy() < 1.0) & (dry_after >= SUMMARY_BINS))[0]:
            end = int(le.last_bin[e]) + 1                       # first bin after the last wet bin
            idx = end + k
            ok = (idx >= 0) & (idx < le.n_bins[e])
            g = s["off"][e] + np.clip(idx, 0, le.n_bins[e] - 1)
            rows_m.append(np.where(ok & s["counted"][g], s["target"][g], np.nan))
            for m in MODELS:
                rows_p[m].append(np.where(ok & s["scored"][g], s["preds"][m][g], np.nan))
            n_le += 1
            links.add((name, le.cml_id[e]))
    M = np.array(rows_m)
    t = k * 5 + 2.5
    with np.errstate(all="ignore"):
        q25, med, q75 = np.nanpercentile(M, [25, 50, 75], axis=0)
        pm = {m: np.nanmedian(np.array(rows_p[m]), axis=0) for m in MODELS}
    yl = ylims([q25, q75] + list(pm.values()))
    fig, axes = plt.subplots(7, 1, figsize=(8.5, 13), sharex=True)
    for ax in axes:
        ax.axvspan(t[0] - 2.5, 0, color=SPAN, lw=0, zorder=0)
        ax.axvline(0, color="#555555", lw=0.8)
        style(ax)
        ax.set_ylim(*yl)
        ax.set_ylabel("dB", fontsize=8)
    axes[0].fill_between(t, q25, q75, color="#bdbdbd", lw=0, label="measured 25-75 %")
    axes[0].plot(t, med, color="black", lw=1.6, label="measured median")
    axes[0].set_title("measured A - c: median and 25-75 % band", fontsize=9, loc="left")
    axes[0].legend(fontsize=7, loc="upper right", frameon=False)
    for ax, m in zip(axes[1:], MODELS):
        ax.plot(t, med, color=GREY, lw=1.0, ls="--", label="measured median")
        ax.plot(t, pm[m], color=COLOR[m], lw=1.8, label=f"{LABEL[m]}: median total")
        ax.set_title(f"{LABEL[m]}: median of the predicted total (WAA + path)", fontsize=9, loc="left")
    axes[1].legend(fontsize=7, loc="upper right", frameon=False)
    axes[-1].set_xlabel("minutes from the end of the last wet 5 min bin (grey: before it)", fontsize=9)
    fig.suptitle(f"TEST, links < 1 km: {n_le} link-events of {len(links)} links (known and new), at least 60 min after the last wet bin\n"
                 "aligned at the end of the last wet bin (typical behaviour; visualization only)", fontsize=10)
    fig.text(0.01, 0.005, NOT_SAVED, fontsize=7, color="#555555")
    fig.tight_layout(rect=(0, 0.015, 1, 0.97))
    stem = P.FIGURES / "waa_compare_summary_short_links"
    save(fig, stem.with_suffix(".png"), dpi=300)
    save(fig, stem.with_suffix(".pdf"))
    plt.close(fig)
    return stem, n_le, len(links)


def main():
    sets = load()
    rows, files = [], []
    kn, nw = sets["known"], sets["new"]
    rk, rn = rmse_per_le(kn, "ours"), rmse_per_le(nw, "ours")
    Lk, Ln = kn["le"].L_km.to_numpy(), nw["le"].L_km.to_numpy()
    for rule, s, Lmask, rm, known in (("rule1_short", kn, Lk < 1.0, rk, True), ("rule2_long", kn, Lk > 3.0, rk, True),
                                      ("rule4_new", nw, Ln < 1.0, rn, False)):
        e, row = select(s, rule, Lmask, rm)
        files.append(event_figure(s, e, rule, known))
        row["figure"] = rel(files[-1].with_suffix(".png"))
        rows.append(row)
    print("   (rule3_storm) not applied: the peak rain rate (radar series) was not saved by the TEST run", flush=True)
    stem, n_le, n_links = summary_figure(sets)
    sel = pd.DataFrame(rows)
    save(sel.round(4), P.TF_SELECTION, index=False)
    L = ["# TEST comparison figures (a visualization only)\n",
         "Written by `research/waa_ksurface/t_figures.py` (decisions.md, \"TEST comparison figures\", logged before the figures). Data: only the saved "
         "outputs of the TEST run (`predictions_test.pkl`); no model rerun, TEST not read again. No fitting, no metric, no decision follows.\n",
         "## Selected link-events (rules written before any curve was seen)\n", md_table(sel.round(4)), "",
         "**Rule 3 (strong storm) not applied:** the peak rain rate needs the radar series, which the TEST run did not save.\n",
         f"## Summary figure\n\n`{rel(stem.with_suffix('.png'))}`: {n_le} link-events of {n_links} links < 1 km (known and new) with at least 60 min "
         "after the last wet bin, aligned at the end of that bin.\n",
         "## Not plotted\n", NOT_SAVED + " The grey span in the per-event figures marks the first to the last wet 5 min bin of the link (data layer); "
         "dry bins inside it are not visible.\n"]
    save("\n".join(L) + "\n", P.TF_REPORT)
    print(f"report: {rel(P.TF_REPORT)}", flush=True)


if __name__ == "__main__":
    main()
