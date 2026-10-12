"""Small helpers for the Markdown reports that the research scripts write."""
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]                 # the repository root


def md_table(df, index=False):
    """A DataFrame as a Markdown table (missing values as empty cells). index=True also writes the index."""
    df = df.reset_index() if index else df
    head = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns)
    return head + "\n" + "\n".join("| " + " | ".join("" if pd.isna(x) else str(x) for x in row) + " |"
                                   for row in df.itertuples(index=False))


def rel(path):
    """A path relative to the repository root, with forward slashes (the way reports name files)."""
    return Path(path).resolve().relative_to(REPO).as_posix()
