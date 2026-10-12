"""Plot style shared by the research figures.

Figures must be readable without colour: use line styles, hatching and text labels; colour only comes second.
"""
INK, INK2 = "#1a1a1a", "#52514e"                           # the two greys used for lines and text


def style_axes(ax):
    """Light grid, no top and right frame lines."""
    ax.grid(alpha=0.25)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
