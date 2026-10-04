"""Shared figure style for the report figures (same palette as the earlier steps: the dataviz reference palette)."""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")

C_BLUE, C_ORANGE, C_AQUA, C_YELLOW, C_MAGENTA, C_RED = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#e34948"
C_MUTED, C_INK, C_SURF, C_GRID, C_SPINE = "#898781", "#0b0b0b", "#fcfcfb", "#e6e5e1", "#d6d5d0"
SPLIT_C = {"train": C_BLUE, "val": C_YELLOW, "test": C_ORANGE, "excluded": C_MUTED}
STATUS_C = {"good": C_AQUA, "degraded": C_YELLOW, "unusable": C_RED}
SUBLINK_C = {"sublink_1": C_BLUE, "sublink_2": C_ORANGE}
FLAG_C = {"missing_rsl": C_MUTED, "missing_tsl": "#b5b3ac", "out_of_range": C_MAGENTA, "rx_floor": C_RED, "spike": C_ORANGE,
          "tsl_step": C_YELLOW, "flatline": "#4a3aa7", "network_gap": C_INK}


def style(ax):
    """Apply the shared axes style."""
    ax.set_facecolor(C_SURF)
    ax.grid(color=C_GRID, lw=0.6)
    ax.tick_params(labelsize=8, colors=C_MUTED)
    for sp in ax.spines.values():
        sp.set_color(C_SPINE)


def label(ax, xlabel, ylabel, title):
    """Set axis labels and a left-aligned title."""
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=9, color=C_INK)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=9, color=C_INK)
    if title:
        ax.set_title(title, fontsize=10, color=C_INK, loc="left")
