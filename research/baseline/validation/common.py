"""Shared pieces of the baseline validation (R1): settings, data loading (with a local cache), masks, the baseline
variants and small statistics helpers. Every check script (c1_ ... c6_) imports this module.

Settings: section `baseline_validation` (and `baseline`) of research/config.yaml.
Data: only through cmlpipe.api.CMLData, on the release pinned in research/config.yaml.
Outputs: research/baseline/outputs/tables/ and outputs/figures/ (FULL sample) or .cache/dev/ (DEV sample);
locations in paths.py, every file described in CATALOG.md.
"""
import os
from pathlib import Path
import sys
import time
import numpy as np
import pandas as pd
import yaml
from cmlpipe.api import CMLData

HERE = Path(__file__).resolve().parent                     # research/baseline/validation
from research.baseline import paths as P                            # noqa: E402  (all file locations)
from research.common.save import save
REPO = P.REPO
CONFIG = yaml.safe_load(open(P.CONFIG_FILE, encoding="utf-8"))
V = CONFIG["baseline_validation"]                          # this step's settings
B = CONFIG["baseline"]                                     # the current method's settings (M3)
V["sample"] = os.environ.get("BV_SAMPLE", V["sample"])     # quick runs: set BV_SAMPLE=DEV (config stays FULL)
SPH = 12                                                   # 5 min slots per hour
VARIANTS = ["M0", "M1", "M2", "M2b", "M3"]

OUT, FIG = P.validation_dirs(V["sample"])                  # tables, figures (FULL: outputs/; DEV: .cache/dev/)
CACHE = P.CACHE                                            # local only (gitignored): loaded data
for p in (OUT, FIG, CACHE):
    p.mkdir(parents=True, exist_ok=True)


def open_data():
    return CMLData.open(REPO / "data" / "release" / CONFIG["release"])


class Timer:
    """with Timer("C1"): ...  -> prints the runtime and appends it to OUT/validation_runtimes.csv"""
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.t0 = time.perf_counter()
        return self

    def __exit__(self, *exc):
        dt = time.perf_counter() - self.t0
        print(f"[{self.name}] runtime {dt:.1f} s")
        f = OUT / "validation_runtimes.csv"
        old = pd.read_csv(f) if f.exists() else pd.DataFrame(columns=["step", "seconds"])
        old = old[old.step != self.name]
        new = pd.concat([old, pd.DataFrame([{"step": self.name, "seconds": round(dt, 1)}])])
        save(new, f, index=False)


# ------------------------------------------------------------------------------------------------- sub-links
def sublink_table(d):
    """Train sub-links that are good or degraded, with their group:
    main = status good (the main numbers); degraded = status degraded; toggling = TSL class toggling (reported
    separately; in this release every toggling train sub-link is also degraded). DEV keeps only the DEV list."""
    s = d.sublinks()
    s = s[s.split_link == "train"].copy()
    s["status"] = s.sublink_id.map(d.status().set_index("sublink_id").status)
    s = s[s.status.isin(["good", "degraded"])]
    s["group"] = np.where((s.status == "good") & (s.tsl_behaviour != "toggling"), "main", "secondary")
    s["is_degraded"] = s.status == "degraded"
    s["is_toggling"] = s.tsl_behaviour == "toggling"
    if V["sample"] == "DEV":
        s = s[s.sublink_id.isin(V["dev_sublinks"])]
    cols = ["sublink_id", "cml_id", "sublink", "status", "tsl_behaviour", "length_class", "frequency_band",
            "length_m", "frequency_ghz", "group", "is_degraded", "is_toggling"]
    return s[cols].reset_index(drop=True)


def groups(sl):
    """The reported groups: name -> list of sub-link ids. Main numbers use only 'main'."""
    g = {"main": sl[sl.group == "main"].sublink_id.tolist(),
         "degraded": sl[sl.is_degraded].sublink_id.tolist(),
         "toggling": sl[sl.is_toggling].sublink_id.tolist()}
    return {k: v for k, v in g.items() if v}


# ------------------------------------------------------------------------------------------------- data panel
def load_panel(d, sl):
    """All time series needed, on the common 5 min grid of the season:
    tl (time x sublink_id, 5 min tl_mean), state (time x cml_id: 0 dry, 1 wet, -1 unknown), cert (network
    certified dry period), met (temperature, humidity, wind_speed, dew point, dew-point depression; hourly SMHI
    values forward-filled within the hour). Cached in .cache/ because loading takes a while."""
    p = d.catalog["period"]
    t0, t1 = p["start"], p["end_exclusive"]
    f_tl = CACHE / "tl_5min.parquet"
    ids = sl.sublink_id.tolist()
    if f_tl.exists():
        tl = pd.read_parquet(f_tl)
    else:
        allsl = d.sublinks()
        allsl = allsl[allsl.split_link == "train"].sublink_id.tolist()   # cache every train sub-link once
        tl = d.cml(allsl, t0, t1, res="5min", variables=["tl_mean"]).tl_mean.to_pandas()
        tl.columns = [str(c) for c in tl.columns]
        tl.to_parquet(f_tl)
    tl = tl[ids].astype("float64")
    state = d.store("link_reference_state_5min").state.load().to_pandas().reindex(tl.index)
    cert = pd.Series(False, index=tl.index)
    for r in d.dry_periods().itertuples():
        cert[(tl.index >= r.start) & (tl.index < r.end)] = True
    met = d.store("met_1h").load()
    m = pd.DataFrame({v: met[v].isel(gauge_id=0).to_series() for v in ["temperature", "humidity", "wind_speed"]})
    m = m.reindex(tl.index.floor("1h")).set_axis(tl.index)          # hourly value carried through the hour
    m["dew_point"] = dew_point(m.temperature, m.humidity)
    m["dpd"] = m.temperature - m.dew_point                           # dew-point depression (0 = saturated)
    return tl, state, cert, m


def dew_point(t_c, rh_pct):
    """Dew point (deg C) from temperature and relative humidity, Magnus formula (b = 17.62, c = 243.12 deg C)."""
    b, c = 17.62, 243.12
    g = np.log(np.clip(rh_pct, 1, 100) / 100) + b * t_c / (c + t_c)
    return c * g / (b - g)


def hours_since_wet(state_col):
    """For every slot: hours from the END of the most recent wet slot (state == 1) to the slot start; inf if none."""
    wet = (state_col.to_numpy() == 1)
    pos = np.arange(len(wet))
    last = np.maximum.accumulate(np.where(wet, pos, -1))
    out = (pos - last - 1) / SPH
    out[last < 0] = np.inf
    out[wet] = -1.0                                                  # the wet slots themselves
    return out


def usable(state_col, skip_h):
    """Dry slot (link state 0) that is at least skip_h after the end of the last wet slot of the link."""
    return (state_col.to_numpy() == 0) & (hours_since_wet(state_col) >= skip_h)


# ------------------------------------------------------------------------------------------------- baselines
def window(tl, ok, i0, i1):
    """Values of tl in slots [i0, i1) that are usable (ok) and not NaN."""
    i0, i1 = max(i0, 0), min(i1, len(tl))
    if i1 <= i0:
        return np.empty(0)
    v, m = tl[i0:i1], ok[i0:i1]
    return v[m & ~np.isnan(v)]


def stat(v, how, min_n):
    if len(v) < min_n:
        return np.nan
    return float(np.mean(v)) if how == "mean" else float(np.median(v))


def variants(tl, ok, m0, i0, i1):
    """Baselines of every variant for an 'event' covering slots [i0, i1).
    Returns {variant: baseline over the event as an array of i1 - i0 values (NaN = missing)} plus M3 details."""
    n = i1 - i0
    vv = V["variants"]
    min_n = int(round(vv["min_dry_h"] * SPH))
    w1, w2 = int(vv["M1_window_h"] * SPH), int(vv["M2_window_h"] * SPH)
    out = {"M0": m0,
           "M1": stat(window(tl, ok, i0 - w1, i0), "median", min_n),
           "M2": stat(window(tl, ok, i0 - w2, i0), "mean", min_n),
           "M2b": stat(window(tl, ok, i0 - w2, i0), "median", min_n)}
    out = {k: np.full(n, v) for k, v in out.items()}
    m3, info = m3_baseline(tl, ok, m0, i0, i1)
    out["M3"] = m3
    return out, info


def m3_baseline(tl, ok, usual, i0, i1):
    """The current method A/B (research/baseline/baseline.py) on positions: mean of usable dry slots in the
    window before / after; line if |after - before| > diff_db unless a window is > excursion_db from the usual
    level (then flat + flag); missing if < min_dry_h before."""
    w, min_n = int(B["v1"]["window_h"] * SPH), int(round(B["min_dry_h"] * SPH))
    n = i1 - i0
    before = stat(window(tl, ok, i0 - w, i0), "mean", min_n)
    after = stat(window(tl, ok, i1, i1 + w), "mean", min_n)
    if np.isnan(before):
        return np.full(n, np.nan), {"m3_method": "missing", "m3_flag": "no_before", "before": before, "after": after}
    exc_b = abs(before - usual) > B["excursion_db"]
    exc_a = (not np.isnan(after)) and abs(after - usual) > B["excursion_db"]
    flag = ";".join([f for f, on in (("excursion_before", exc_b), ("excursion_after", exc_a)) if on])
    if np.isnan(after):
        meth, flag = "A", "no_after"
    elif exc_b or exc_a:
        meth = "A"
    elif abs(after - before) > B["v1"]["diff_db"]:
        meth = "B"
    else:
        meth = "A"
    base = np.full(n, before) if meth == "A" else before + (after - before) * np.arange(n) / n
    return base, {"m3_method": meth, "m3_flag": flag, "before": before, "after": after}


# ------------------------------------------------------------------------------------------------- statistics
def summary(e):
    """n, median, p10, p90, median |e|, RMS of the non-NaN values."""
    e = np.asarray(e, float)
    e = e[~np.isnan(e)]
    if not len(e):
        return {"n_slots": 0, "median": np.nan, "p10": np.nan, "p90": np.nan, "median_abs": np.nan, "rms": np.nan}
    return {"n_slots": len(e), "median": np.median(e), "p10": np.percentile(e, 10), "p90": np.percentile(e, 90),
            "median_abs": np.median(np.abs(e)), "rms": float(np.sqrt(np.mean(e ** 2)))}


def wquantile(x, w, q):
    """Weighted quantile (q in 0..1) of x with weights w (NaN x dropped)."""
    x, w = np.asarray(x, float), np.asarray(w, float)
    m = ~np.isnan(x)
    x, w = x[m], w[m]
    if not len(x):
        return np.nan
    o = np.argsort(x)
    x, w = x[o], w[o]
    c = (np.cumsum(w) - 0.5 * w) / w.sum()
    return float(np.interp(q, c, x))


def link_midpoints(d):
    """Path midpoint (x, y in m, EPSG:3006) per cml_id."""
    sites = d.table("sites").set_index("site_id")
    L = d.links().set_index("cml_id")
    x = (sites.x_3006_m.reindex(L.site_0_id).to_numpy() + sites.x_3006_m.reindex(L.site_1_id).to_numpy()) / 2
    y = (sites.y_3006_m.reindex(L.site_0_id).to_numpy() + sites.y_3006_m.reindex(L.site_1_id).to_numpy()) / 2
    return pd.DataFrame({"x": x, "y": y}, index=L.index)


def fmt_table(df, digits=2):
    """A DataFrame as a markdown table (no extra dependency)."""
    df = df.reset_index() if df.index.name or isinstance(df.index, pd.MultiIndex) else df
    rows = [" | ".join(map(str, df.columns)), " | ".join("---" for _ in df.columns)]
    for r in df.itertuples(index=False):
        rows.append(" | ".join(f"{v:.{digits}f}" if isinstance(v, (float, np.floating)) else str(v) for v in r))
    return "\n".join(f"| {r} |" for r in rows)
