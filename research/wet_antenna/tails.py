"""Small helpers for working with the real drying tails (used by the selection steps 2, 2b, 2c)."""
import numpy as np
import pandas as pd


def first_run(cond, n):
    """Start of the first run of >= n consecutive True values (None if there is none)."""
    c = np.concatenate([[0], cond.astype(int), [0]])
    starts, ends = np.nonzero(np.diff(c) == 1)[0], np.nonzero(np.diff(c) == -1)[0]
    ok = np.nonzero(ends - starts >= n)[0]
    return int(starts[ok[0]]) if len(ok) else None


def window_mean(a, k, w, end):
    """Mean of a[k:k + w] if that window lies inside the usable part [0, end), else NaN."""
    if k + w > end:
        return np.nan
    v = a[k:k + w]
    return float(np.nanmean(v)) if np.isfinite(v).any() else np.nan


def minutes(t, t0):
    """Whole minutes from time t0 to time t."""
    return int(round((t - t0) / pd.Timedelta("1min")))


def set_counts(df):
    """How many tails, sub-links, links and events a set of tails holds."""
    return {"tails": len(df), "sub-links": df.sublink_id.nunique(), "links": df.cml_id.nunique(),
            "events": df.event_id.nunique()}


def spearman(x, y):
    """Spearman rank correlation of two Series."""
    return x.rank().corr(y.rank())
