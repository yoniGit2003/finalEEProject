"""The figures of docs/WAA_research_overview.pdf, made from saved research tables only (nothing is fitted here).

Every figure reads tables written by the wet_antenna steps (research/wet_antenna/outputs/tables and products). Black and
grey only: line styles, markers and text labels carry the meaning, never colour. Run from the repository root:
    python docs/WAA_research_overview_figures/figs.py
then build the PDF with: python docs/build_research_log.py WAA_research_overview
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
T = REPO / "research" / "wet_antenna" / "outputs" / "tables"
PR = REPO / "research" / "wet_antenna" / "outputs" / "products"
INK, GREY = "#1a1a1a", "#7a7a7a"
BANDS = ["28 GHz", "32 GHz", "38 GHz"]
MK = {"28 GHz": "o", "32 GHz": "s", "38 GHz": "^"}
LS = {"28 GHz": "-", "32 GHz": "--", "38 GHz": ":"}
CLASSES = ["L < 1 km", "1-3 km", "L > 3 km"]

plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9, "legend.fontsize": 8, "xtick.labelsize": 8,
                     "ytick.labelsize": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.alpha": 0.25, "savefig.dpi": 200, "savefig.bbox": "tight"})


def save(fig, name):
    fig.savefig(OUT / f"{name}.png")
    plt.close(fig)


def example_event():
    """Figure 1: one TRAIN storm on a short and a long link (measured, model total, W, path)."""
    e = pd.read_csv(T / "formulate_s9_example_event.csv")
    fig, axes = plt.subplots(3, 1, figsize=(8, 6.2), sharex=True, gridspec_kw={"height_ratios": [0.6, 1, 1]})
    s = e[e.link == "short"]
    peak = int(s.minute[s.rain_mm_h.idxmax()])                         # the main shower: 4 h before to 5 h after its peak
    axes[0].fill_between(s.minute, s.rain_mm_h, color=GREY, step="mid")
    axes[0].set_ylabel("rain (mm/h)")
    for ax, lab in zip(axes[1:], ("short", "long")):
        g = e[e.link == lab]
        ax.plot(g.minute, g.A_minus_c, "-", color=GREY, lw=0.8, label="measured attenuation")
        ax.plot(g.minute, g.W + g.path, "-", color=INK, lw=1.3, label="model, total")
        ax.plot(g.minute, g.W, "--", color=INK, lw=1.2, label="model, wet antenna part")
        ax.plot(g.minute, g.path, ":", color=INK, lw=1.4, label="model, rain on the path part")
        ax.set_ylabel("dB")
        ax.text(0.99, 0.9, f"link length {g.L_km.iloc[0]:.2f} km", transform=ax.transAxes, ha="right", fontsize=8)
    axes[1].legend(frameon=False, loc="upper left", fontsize=7)
    axes[2].set_xlabel("time (min)")
    axes[2].set_xlim(peak - 240, peak + 300)
    save(fig, "example_event")


def length_test():
    """Figure 2: the length test on the real TRAIN data."""
    pts = pd.read_csv(PR / "system_id_s3_steady_link_events.csv")
    fits = pd.read_csv(T / "system_id_s3_length_fits.csv")
    fits = fits[fits.variant == "TRAIN"]
    pts = pts[(pts.variant == "TRAIN") & (pts.band_name == "28 GHz")]
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))
    ax = axes[0]
    for cls, mk in (("0.5-1", "o"), ("2-5", "s"), ("5-10", "^"), ("10-20", "D")):
        g = pts[pts["class"] == cls]
        f = fits[(fits.band == "28 GHz") & (fits["class"] == cls)].iloc[0]
        ax.scatter(g.L_km, g.A_med, marker=mk, s=6, facecolors="none", edgecolors=GREY, linewidths=0.5)
        L = np.array([0, g.L_km.max()])
        ax.plot(L, f.I + f.s * L, "-", color=INK, lw=1.2)
        ax.plot([0], [f.I], mk, color=INK, ms=5)
        ax.text(L[1] + 0.1, f.I + f.s * L[1], f"{cls} mm/h", va="center", fontsize=7)
    ax.set_xlabel("link length L (km)")
    ax.set_ylabel("A minus c (dB)")
    ax.set_title("(a) one line per rain group (28 GHz)")
    ax.set_ylim(-2, 20)
    for ax, q, lab in ((axes[1], "I", "height at L = 0 (dB)"), (axes[2], "s", "slope (dB/km)")):
        for k, b in enumerate(BANDS):
            g = fits[(fits.band == b) & np.isfinite(fits[q])].sort_values("R_cell")
            g = g[g.links_under_1km >= 8] if b != "32 GHz" else g[g["class"] != "10-20"]
            x = g.R_cell * (1 + (k - 1) * 0.03)
            ax.errorbar(x, g[q], yerr=[g[q] - g[f"{q}_lo"], g[f"{q}_hi"] - g[q]], fmt=MK[b] + LS[b], color=INK, mfc="white", ms=4,
                        capsize=2, lw=1, label=b)
        ax.set_xlabel("rain rate (mm/h)")
        ax.set_ylabel(lab)
        ax.legend(frameon=False)
    axes[1].set_title("(b) heights give W(R)")
    axes[2].set_title("(c) slopes give P(R)")
    fig.tight_layout()
    save(fig, "length_test")


def memory():
    """Figure 3: an illustration computed from the equation (60 min of 5 mm/h)."""
    m = pd.read_csv(T / "formulate_s9_memory.csv")
    fig, (a, b) = plt.subplots(2, 1, figsize=(7, 4), sharex=True, gridspec_kw={"height_ratios": [0.5, 1]})
    a.fill_between(m.minute, m.rain_mm_h, color=GREY, step="mid")
    a.set_ylabel("rain (mm/h)")
    b.plot(m.minute, m.W, "-", color=INK, lw=1.5, label="total W")
    b.plot(m.minute, m.W_fast, "--", color=INK, lw=1.1, label="fast part (72 %)")
    b.plot(m.minute, m.W_slow, ":", color=INK, lw=1.4, label="slow part (28 %)")
    b.set_xlabel("time (min)")
    b.set_ylabel("W (dB)")
    b.legend(frameon=False)
    save(fig, "memory")


def F_and_P():
    """Figure 4: (a) F, free points and the formula per band, (b) P, free points and c R per band (linear axes)."""
    boot = pd.read_csv(T / "tighten_s6_bootstrap.csv").set_index("quantity")
    band = pd.read_csv(T / "formulate_s9_closed_f_band.csv")
    eq = json.loads((PR / "waa_equation_final.json").read_text())
    knots = np.array([0.5, 1, 2, 5, 10, 20, 40])
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    for ax, kind in zip(axes, ("F", "P")):
        for k, b in enumerate(BANDS):
            v = np.array([boot.value[f"{kind}[{b}]({x:g})"] for x in knots])
            lo = np.array([boot.boot_lo[f"{kind}[{b}]({x:g})"] for x in knots])
            hi = np.array([boot.boot_hi[f"{kind}[{b}]({x:g})"] for x in knots])
            x = knots + (k - 1) * 0.6
            ax.errorbar(x, v, yerr=[v - lo, hi - v], fmt=MK[b], color=INK, ecolor=GREY, mfc="white", ms=4, capsize=2, label=b)
            g = band[band.band == b]
            y = g.F if kind == "F" else eq["c_band_db_km_per_mmh"][b]["value"] * g.R
            ax.plot(g.R, y, LS[b], color=INK, lw=1.2)
        ax.set_xlabel("rain rate R (mm/h)")
        ax.set_xlim(-1, 42)
        handles = [matplotlib.lines.Line2D([], [], color=INK, ls=LS[b], marker=MK[b], mfc="white", label=b) for b in BANDS]
        ax.legend(handles=handles, frameon=False, loc="upper left", title="points: free curve, lines: formula", title_fontsize=7)
    axes[0].axvspan(10, 42, color="0.92", zorder=0)
    axes[0].text(25, 0.4, "few data", ha="center", fontsize=8)
    axes[0].set_ylabel("F(R) (dB)")
    axes[0].set_title("(a) steady-state WAA F(R), lines: A (R/5)^0.54")
    axes[1].set_ylabel("P(R) (dB/km)")
    axes[1].set_title("(b) path attenuation per km P(R), lines: c R")
    fig.tight_layout()
    save(fig, "F_P")


def a_link():
    """Figure 5: the link multiplier of the final equation against the link length."""
    lk = pd.read_csv(T / "validate_s7_links.csv")
    fig, ax = plt.subplots(figsize=(6, 3.2))
    for b in BANDS:
        g = lk[lk.band == b]
        ax.scatter(g.L_km, g["a_link closed F"], marker=MK[b], s=14, facecolors="none", edgecolors=INK, linewidths=0.6, label=b)
    ax.set_xscale("log")
    ax.axhline(1, color=GREY, lw=0.8)
    ax.set_xlabel("link length L (km)")
    ax.set_ylabel("link factor a_link")
    ax.legend(frameon=False, title="band (marker)", title_fontsize=7)
    save(fig, "a_link")


def dry_path():
    """Figure 6: the wet antenna measured on dry paths after rain against the model, per length class."""
    dc = pd.read_csv(T / "formulate_s9_dry_curves.csv")
    fig, axes = plt.subplots(1, 3, figsize=(10, 3), sharey=True)
    for ax, cls in zip(axes, CLASSES):
        g = dc[dc.length == cls]
        ax.plot(g.minute, g.measured, "-", color=INK, lw=1.3, label="measured A minus c")
        ax.plot(g.minute, g.model, "--", color=INK, lw=1.2, label="model W")
        ax.set_title(cls)
        ax.set_xlabel("minutes after the path became dry")
    axes[0].set_ylabel("median W (dB)")
    axes[0].legend(frameon=False)
    fig.tight_layout()
    save(fig, "dry_path")


if __name__ == "__main__":
    for f in (example_event, length_test, memory, F_and_P, a_link, dry_path):
        f()
    print("figures written to", OUT)
