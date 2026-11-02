"""The figures of docs/WAA_research_overview.pdf (the two-column paper version), drawn for one column (8.6 cm) from saved
research tables only (nothing is fitted here).

The first version's figures (figs.py, used by WAA_research_overview_v1.tex) are left unchanged; this script writes v2_*.png
next to them. Black and grey only: line styles, markers and text labels carry the meaning, never colour. Fonts are 8 pt or
larger at the printed size. Run from the repository root:
    python docs/WAA_research_overview_figures/figs_v2.py
then build the PDF with: python docs/build_research_log.py WAA_research_overview
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
T = REPO / "research" / "wet_antenna" / "outputs" / "tables"
PR = REPO / "research" / "wet_antenna" / "outputs" / "products"
INK, GREY = "#1a1a1a", "#7a7a7a"
BANDS = ["28 GHz", "32 GHz", "38 GHz"]
MK = {"28 GHz": "o", "32 GHz": "s", "38 GHz": "^"}
LS = {"28 GHz": "-", "32 GHz": "--", "38 GHz": ":"}
CLASSES = ["L < 1 km", "1-3 km", "L > 3 km"]
W1 = 3.39                                          # one column, inches (8.6 cm)

plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 8, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.alpha": 0.25, "savefig.dpi": 300, "lines.linewidth": 1.0, "font.family": "serif"})


def save(fig, name):
    fig.savefig(OUT / f"v2_{name}.png", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def panel(ax, s):
    ax.set_title(s, loc="left", fontsize=8, pad=2)


def example_event():
    """One TRAIN storm on a short and a long link: measured, model total, W and path part."""
    e = pd.read_csv(T / "formulate_s9_example_event.csv")
    fig, axes = plt.subplots(3, 1, figsize=(W1, 3.9), sharex=True, gridspec_kw={"height_ratios": [0.55, 1, 1]})
    s = e[e.link == "short"]
    peak = int(s.minute[s.rain_mm_h.idxmax()])
    axes[0].fill_between(s.minute, s.rain_mm_h, color=GREY, step="mid", lw=0)
    axes[0].set_ylabel("$R$ (mm/h)")
    panel(axes[0], "(a)")
    for ax, lab, p in zip(axes[1:], ("short", "long"), ("(b)", "(c)")):
        g = e[e.link == lab]
        ax.plot(g.minute, g.A_minus_c, "-", color=GREY, lw=0.6, label="measured $A-c$")
        ax.plot(g.minute, g.W + g.path, "-", color=INK, lw=1.1, label="model, total")
        ax.plot(g.minute, g.W, "--", color=INK, lw=1.0, label="model, $W$")
        ax.plot(g.minute, g.path, ":", color=INK, lw=1.3, label="model, $L\\,P$")
        ax.set_ylabel("dB")
        panel(ax, f"{p} $L$ = {g.L_km.iloc[0]:.2f} km")
    axes[2].set_xlabel("time (min)")
    axes[2].set_xlim(peak - 240, peak + 300)
    fig.tight_layout(h_pad=0.3)
    h, lab = axes[1].get_legend_handles_labels()
    fig.legend(h, lab, frameon=False, loc="upper center", bbox_to_anchor=(0.55, 0.0), ncol=2, fontsize=8, handlelength=2.2)
    save(fig, "example_event")


def length_test():
    """The model-free length test on TRAIN: (a) lines per rain class at 28 GHz, (b) intercepts W(R), (c) slopes P(R)."""
    pts = pd.read_csv(PR / "system_id_s3_steady_link_events.csv")
    fits = pd.read_csv(T / "system_id_s3_length_fits.csv")
    fits = fits[fits.variant == "TRAIN"]
    pts = pts[(pts.variant == "TRAIN") & (pts.band_name == "28 GHz")]
    fig, axes = plt.subplots(3, 1, figsize=(W1, 5.4), gridspec_kw={"height_ratios": [1.2, 1, 1]})
    ax = axes[0]
    for cls, mk in (("0.5-1", "o"), ("2-5", "s"), ("5-10", "^"), ("10-20", "D")):
        g = pts[pts["class"] == cls]
        f = fits[(fits.band == "28 GHz") & (fits["class"] == cls)].iloc[0]
        ax.scatter(g.L_km, g.A_med, marker=mk, s=4, facecolors="none", edgecolors=GREY, linewidths=0.4)
        L = np.array([0, min(g.L_km.max(), 8.5)])
        ax.plot(L, f.I + f.s * L, "-", color=INK, lw=1.0)
        ax.plot([0], [f.I], mk, color=INK, ms=4)
        ax.text(L[1] + 0.1, f.I + f.s * L[1], f"{cls} mm/h", va="center", fontsize=8)
    ax.set_xlabel("link length $L$ (km)")
    ax.set_ylabel("$A-c$ (dB)")
    ax.set_xlim(-0.2, 12.0)
    ax.set_ylim(-2, 16)
    panel(ax, "(a) 28 GHz, one line per rain class")
    for ax, q, lab, p in ((axes[1], "I", "intercept $I$ (dB)", "(b) intercepts: $W(R)$"), (axes[2], "s", "slope $s$ (dB/km)", "(c) slopes: $P(R)$")):
        for k, b in enumerate(BANDS):
            g = fits[(fits.band == b) & np.isfinite(fits[q])].sort_values("R_cell")
            g = g[g.links_under_1km >= 8] if b != "32 GHz" else g[g["class"] != "10-20"]
            x = g.R_cell * (1 + (k - 1) * 0.04)
            ax.errorbar(x, g[q], yerr=[g[q] - g[f"{q}_lo"], g[f"{q}_hi"] - g[q]], fmt=MK[b] + LS[b], color=INK, mfc="white", ms=3.5,
                        capsize=1.5, lw=0.9, label=b)
        ax.set_xscale("log")
        ax.set_ylabel(lab)
        panel(ax, p)
    axes[1].legend(frameon=False, loc="upper left", fontsize=8)
    axes[2].set_xlabel("rain rate $R$ (mm/h)")
    fig.tight_layout(h_pad=0.4)
    save(fig, "length_test")


def memory():
    """Response of the equation to 60 min of 5 mm/h (computed from the equation, not measured)."""
    m = pd.read_csv(T / "formulate_s9_memory.csv")
    fig, (a, b) = plt.subplots(2, 1, figsize=(W1, 2.4), sharex=True, gridspec_kw={"height_ratios": [0.45, 1]})
    a.fill_between(m.minute, m.rain_mm_h, color=GREY, step="mid", lw=0)
    a.set_ylabel("$R$ (mm/h)")
    b.plot(m.minute, m.W, "-", color=INK, lw=1.3, label="$W = W_f + W_s$")
    b.plot(m.minute, m.W_fast, "--", color=INK, lw=1.0, label="fast store $W_f$")
    b.plot(m.minute, m.W_slow, ":", color=INK, lw=1.3, label="slow store $W_s$")
    b.set_xlabel("time (min)")
    b.set_ylabel("$W$ (dB)")
    b.legend(frameon=False, fontsize=8)
    fig.tight_layout(h_pad=0.2)
    save(fig, "memory")


def F_and_P():
    """(a) F: free points and A (R/5)^b per band; (b) P: free points and c R per band."""
    boot = pd.read_csv(T / "tighten_s6_bootstrap.csv").set_index("quantity")
    band = pd.read_csv(T / "formulate_s9_closed_f_band.csv")
    eq = json.loads((PR / "waa_equation_final.json").read_text())
    knots = np.array([0.5, 1, 2, 5, 10, 20, 40])
    fig, axes = plt.subplots(2, 1, figsize=(W1, 4.1))
    for ax, kind, p in zip(axes, ("F", "P"), ("(a)", "(b)")):
        for k, b in enumerate(BANDS):
            v = np.array([boot.value[f"{kind}[{b}]({x:g})"] for x in knots])
            lo = np.array([boot.boot_lo[f"{kind}[{b}]({x:g})"] for x in knots])
            hi = np.array([boot.boot_hi[f"{kind}[{b}]({x:g})"] for x in knots])
            x = knots + (k - 1) * 0.6
            ax.errorbar(x, v, yerr=[v - lo, hi - v], fmt=MK[b], color=INK, ecolor=GREY, mfc="white", ms=3.5, capsize=1.5, lw=0.8)
            g = band[band.band == b]
            y = g.F if kind == "F" else eq["c_band_db_km_per_mmh"][b]["value"] * g.R
            ax.plot(g.R, y, LS[b], color=INK, lw=1.0)
        ax.set_xlim(-1, 42)
        handles = [Line2D([], [], color=INK, ls=LS[b], marker=MK[b], mfc="white", ms=3.5, label=b) for b in BANDS]
        ax.legend(handles=handles, frameon=False, loc="upper left", fontsize=8)
        panel(ax, p)
    axes[0].axvspan(10, 42, color="0.92", zorder=0)
    axes[0].text(30, 0.4, "few data (0.4 % of minutes)", ha="center", fontsize=8)
    axes[0].set_ylabel("$F(R)$ (dB)")
    axes[1].set_ylabel("$P(R)$ (dB/km)")
    axes[1].set_xlabel("rain rate $R$ (mm/h)")
    fig.tight_layout(h_pad=0.4)
    save(fig, "F_P")


def a_link():
    """The link factor of the final equation against the link length."""
    lk = pd.read_csv(T / "validate_s7_links.csv")
    fig, ax = plt.subplots(figsize=(W1, 2.0))
    for b in BANDS:
        g = lk[lk.band == b]
        ax.scatter(g.L_km, g["a_link closed F"], marker=MK[b], s=9, facecolors="none", edgecolors=INK, linewidths=0.5, label=b)
    ax.set_xscale("log")
    ax.axhline(1, color=GREY, lw=0.8)
    ax.set_xlabel("link length $L$ (km)")
    ax.set_ylabel("$a_{\\mathrm{link}}$")
    ax.legend(frameon=False, fontsize=8, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    save(fig, "a_link")


def dry_path():
    """WAA measured on dry paths after rain against the model, per length class (TRAIN and VAL windows of S7)."""
    dc = pd.read_csv(T / "formulate_s9_dry_curves.csv")
    fig, ax = plt.subplots(figsize=(W1, 2.3))
    styles = {"L < 1 km": ("o", 0), "1-3 km": ("s", 2), "L > 3 km": ("^", 4)}
    for cls in CLASSES:
        g = dc[dc.length == cls]
        mk, off = styles[cls]
        ax.plot(g.minute, g.measured, "-", color=INK, lw=1.0, marker=mk, markevery=(off, 6), ms=3.5, mfc="white")
        ax.plot(g.minute, g.model, "--", color=GREY, lw=1.0, marker=mk, markevery=(off + 3, 6), ms=3.5, mfc=GREY)
    h = [Line2D([], [], color=INK, ls="-", label="measured (solid)"), Line2D([], [], color=GREY, ls="--", label="model $W$ (dashed)")]
    h += [Line2D([], [], color=INK, ls="none", marker=styles[c][0], mfc="white", ms=3.5, label=c) for c in CLASSES]
    ax.legend(handles=h, frameon=False, fontsize=8, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.0), columnspacing=0.8, handlelength=1.8)
    ax.set_xlabel("minutes after the path became dry")
    ax.set_ylabel("median $W$ (dB)")
    save(fig, "dry_path")


NAMES = {"ours": "ours", "B0": "no WAA", "B1": "constant WAA", "B2": "no memory", "Schleiss a": "Schleiss (a)",
         "Schleiss b": "Schleiss (b)", "Pastorek a": "Pastorek (a)", "Pastorek b": "Pastorek (b)"}


def accuracy():
    """TEST RMSE of every model, all links (filled) and links below 1 km (open), with 90 % event-bootstrap intervals."""
    rm = pd.read_csv(T / "literature_s8c_rmse.csv")
    order = ["ours", "Schleiss b", "Pastorek b", "B2", "Pastorek a", "B1", "Schleiss a", "B0"]
    fig, ax = plt.subplots(figsize=(W1, 2.6))
    y = np.arange(len(order))[::-1]
    for sub, mk, mf, dy in (("all", "o", INK, 0.14), ("L < 1 km", "s", "white", -0.14)):
        g = rm[(rm.set == "TEST") & (rm.subset == sub)].set_index("model").loc[order]
        ax.errorbar(g.rmse, y + dy, xerr=[g.rmse - g.lo, g.hi - g.rmse], fmt=mk, color=INK, mfc=mf, ms=3.5, capsize=1.5, lw=0.8)
        ax.axvline(g.loc["ours", "rmse"], color=GREY, ls=":" if sub == "all" else "--", lw=0.8)
    ax.set_yticks(y, [NAMES[m] for m in order])
    ax.set_xlabel("TEST RMSE of the 5 min attenuation (dB)")
    h = [Line2D([], [], color=INK, ls="none", marker="o", ms=3.5, label="all links (dotted: ours)"),
         Line2D([], [], color=INK, ls="none", marker="s", mfc="white", ms=3.5, label="$L < 1$ km (dashed: ours)")]
    ax.legend(handles=h, frameon=False, fontsize=8, loc="lower center", bbox_to_anchor=(0.4, 1.0), ncol=2, handletextpad=0.2)
    ax.grid(axis="y", alpha=0)
    save(fig, "accuracy")


def after_rain():
    """TEST dry-path WAA after rain: measured against the models (median over all windows)."""
    cv = pd.read_csv(T / "comparison_all_models_dry_curves.csv")
    cv = cv[cv.t_min <= 40]
    fig, ax = plt.subplots(figsize=(W1, 2.3))
    ax.plot(cv.t_min, cv.measured, "-", color=INK, lw=1.8, label="measured")
    for m, ls, mk, lab in (("ours", "--", None, "ours"), ("Schleiss b", "-.", None, "Schleiss (b)"),
                           ("Pastorek b", ":", None, "Pastorek (b)"), ("B2", "-", "x", "no memory")):
        ax.plot(cv.t_min, cv[m], ls, color=INK if mk is None else GREY, lw=1.1, marker=mk, markevery=5, ms=3.5, label=lab)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.set_xlim(0, 40)
    ax.set_xlabel("minutes after the path became dry")
    ax.set_ylabel("median $W$ (dB)")
    save(fig, "after_rain")


if __name__ == "__main__":
    for f in (example_event, length_test, memory, F_and_P, a_link, dry_path, accuracy):
        f()
    if (T / "comparison_all_models_dry_curves.csv").exists():
        after_rain()
    print("figures written to", OUT)
