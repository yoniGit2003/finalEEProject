"""Small helpers for working with the real drying tails (used by the selection steps 2, 2b, 2c and by steps 4, 4b, 4c, 4c2)."""
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


def raw_rate(d, edges, smooth):
    """Model-free look (information): mean drying rate per level of A. A is smoothed per tail (centred mean of
    `smooth` minutes), dA/dt = central difference; per bin of smoothed A: k = -mean(dA/dt) / median(A)."""
    d = d.sort_values(["tail_id", "minutes_since_t0_trace"])
    sm = d.groupby("tail_id").A_db.transform(lambda a: a.rolling(smooth, center=True, min_periods=smooth - 2).mean())
    dd = sm.groupby(d.tail_id).transform(lambda v: (v.shift(-1) - v.shift(1)) / 2)
    x = pd.DataFrame({"tail_id": d.tail_id, "A": sm, "dadt": dd}).dropna()
    g = x.groupby(pd.cut(x.A, edges), observed=True)
    out = g.agg(tail_minutes=("dadt", "size"), tails=("tail_id", "nunique"), median_A_db=("A", "median"),
                mean_dadt_db_per_min=("dadt", "mean"), median_dadt_db_per_min=("dadt", "median")).reset_index(names="A_bin_db")
    out["k_from_mean_per_min"] = -out.mean_dadt_db_per_min / out.median_A_db
    out["A_bin_db"] = out.A_bin_db.astype(str)
    return out.round(4)


def classify(a, CL):
    """Class of one tail from its first 60 minutes of A (1 min values, NaN allowed). CL = the rules, config
    `wet_antenna.diagnose_stalled_tails.classes`. Returns (the class and its numbers, the smoothed A)."""
    sm = pd.Series(a).rolling(CL["smooth_min"], center=True, min_periods=CL["smooth_min"] // 2 + 1).median().to_numpy()
    a_start = float(np.nanmean(a[:CL["start_minutes"]]))
    ok = np.isfinite(sm)
    run_min = pd.Series(sm).cummin().to_numpy()
    first = lambda cond: int(np.argmax(cond)) if cond.any() else np.nan   # noqa: E731  first minute a condition holds
    t_half = first(ok & (sm < CL["half"] * a_start))
    t_80 = first(ok & (sm < CL["plateau_level"] * a_start))
    back_up = float(np.nanmax(sm - run_min))               # largest climb above the lowest level reached so far
    rise = float(np.nanmax(sm) - a_start)
    is_decay = bool(t_half <= CL["half_within_min"] and back_up <= CL["rise_db"])
    is_plateau = bool(t_80 >= CL["plateau_min"] and np.isfinite(t_half))
    is_rise = bool(rise > CL["rise_db"])
    cls = "decays" if is_decay else "plateau-then-decays" if is_plateau else "rises" if is_rise else "flat / slow"
    return {"class": cls, "A_start_db": round(a_start, 3), "meets_decays": is_decay, "meets_plateau": is_plateau,
            "meets_rises": is_rise, "min_to_half": t_half, "min_to_80pct": t_80,
            "transition_level_db": round(float(sm[int(t_80)]), 3) if np.isfinite(t_80) else np.nan,
            "rise_db": round(rise, 3), "rise_at_min": int(np.nanargmax(sm)), "climb_above_running_min_db": round(back_up, 3)}, sm



def reference_traces(d, cml_ids, gauge_m, near_km, city_km, smhi_km):
    """What the rain references say per link, on the 5 min grid of the radar (d = the open data release).
    trace  the rain trace of tail set v1: radar R_max > 0 on the path, or a gauge within gauge_m of the path > 0
    nan    the radar on the path is NaN (v1 counts this as "unknown", like a trace)
    broad  rain NEAR the link: radar R_max > 0 on the path of another link whose midpoint is within near_km, or a
           City gauge within city_km of the path > 0, or the SMHI gauge > 0 if it is within smhi_km
    near_max / path_max  the largest radar R_max (mm/h) on the paths of those other links / on the own path
    Returns a dict: "T5" (the times) and, per link, the arrays above."""
    rmax5 = d.store("radar_along_links_5min").R_max.load().to_pandas()
    T5 = pd.DatetimeIndex(rmax5.index)
    r_col = {c: j for j, c in enumerate(rmax5.columns)}
    lgm, links, sites = d.table("link_gauge_matches"), d.links().set_index("cml_id"), d.table("sites").set_index("site_id")
    city, smhi = d.store("gauges_city_5min").amount_mm.load().to_pandas(), d.store("gauges_smhi_15min").amount_mm.load().to_pandas().reindex(T5, method="ffill", limit=2)
    g_amt = pd.concat([city, smhi], axis=1)
    near1 = lgm[lgm.distance_to_path_m <= gauge_m].groupby("cml_id").gauge_id.apply(list).to_dict()      # v1: gauges within 1 km
    far = lgm[(lgm.network == "city") & (lgm.distance_to_path_m <= city_km * 1000)].groupby("cml_id").gauge_id.apply(list).to_dict()
    xy = np.column_stack([(sites.x_3006_m.reindex(links.site_0_id).to_numpy() + sites.x_3006_m.reindex(links.site_1_id).to_numpy()) / 2,
                          (sites.y_3006_m.reindex(links.site_0_id).to_numpy() + sites.y_3006_m.reindex(links.site_1_id).to_numpy()) / 2])
    xy = pd.DataFrame(xy, index=links.index).reindex(rmax5.columns).to_numpy()
    near_link = np.hypot(xy[:, None, 0] - xy[None, :, 0], xy[:, None, 1] - xy[None, :, 1]) < near_km * 1000
    np.fill_diagonal(near_link, False)
    R = rmax5.to_numpy()
    wet5 = (rmax5.fillna(0) > 0).to_numpy()
    out = {"T5": T5, "trace": {}, "nan": {}, "broad": {}, "near_max": {}, "path_max": {}}
    for c in cml_ids:
        t = wet5[:, r_col[c]].copy()
        for g in near1.get(c, []):
            t |= (g_amt[g].fillna(0) > 0).to_numpy()
        out["trace"][c], out["nan"][c] = t, rmax5[c].isna().to_numpy()
        b = wet5[:, near_link[r_col[c]]].any(axis=1)                    # radar on the path of a link within near_km
        for g in far.get(c, []):
            b |= (city[g].fillna(0) > 0).to_numpy()
        if links.nearest_smhi_gauge_distance_m[c] <= smhi_km * 1000:
            b |= (smhi.iloc[:, 0].fillna(0) > 0).to_numpy()
        out["broad"][c] = b
        near = R[:, near_link[r_col[c]]]
        out["near_max"][c] = np.fmax.reduce(near, axis=1) if near.shape[1] else np.full(len(R), np.nan)
        out["path_max"][c] = R[:, r_col[c]]
    return out


def tail_from_start(a, sm10, k, unk, i05, k_len, len_label, F, a_window, extra=None):
    """The tail of one (event, sub-link) pair when it is started at minute k, with the END RULES of tail set v1
    (freeze_tail_set.py): the earliest of the step 2b end (k_len, named len_label), the next rain trace on the path,
    and the rise cut v2 counted from the start. Returns the end, the usable length and A_start.
    a      A per minute from t0 (the step 1 tail);  sm10 = its trailing median over F["rise_smooth_min"] minutes
    unk    rain trace or radar NaN on the path, per 5 min bin of the season;  i05 = the bin of t0
    F      config wet_antenna.freeze;  a_window = minutes of the A_start mean
    extra  optional {name: minute} of further ends (e.g. the next rain nearby); the earliest end wins."""
    L, m = len(a), np.arange(len(a))
    nxt = np.nonzero(unk[i05 + k // 5:])[0]
    k_s = (k // 5 + int(nxt[0])) * 5 if len(nxt) else 10 ** 9      # first rain trace (or radar NaN) after the start
    if k_s < k:                                            # only when the start lies inside a trace bin
        k_s = k
    mask = (m >= k) & np.isfinite(sm10)                             # rise cut v2, from the start
    excess = sm10 - np.minimum.accumulate(np.where(mask, sm10, np.inf))
    k_run, k_rise = first_run(mask & (excess >= F["rise_db"]), F["rise_min_duration_min"]), 10 ** 9
    if k_run is not None:
        quiet = np.nonzero(mask[:k_run] & (excess[:k_run] < F["rise_quiet_db"]))[0]
        k_rise = int(quiet[-1]) + 1 if len(quiet) else k
    cands = {len_label: k_len, "rain_trace": k_s, "rise_v2": k_rise, **(extra or {})}
    by = min(cands, key=lambda q: cands[q])
    end = max(min(cands[by], L), k)
    return {"k": k, "end": end, "usable_min": end - k, "truncated_by": by, "k_trace": k_s,
            "A_start_db": window_mean(a, k, a_window, L)}


# ITU-R P.838-3: specific rain attenuation k = a R^b (dB/km, R in mm/h). log10(a) and b are sums of Gaussians in
# log10(f) plus a line; one coefficient set per polarization. Each row: (a_j, b_j, c_j); then (m, c) of the line.
_P838 = {"aH": ([(-5.33980, -0.10008, 1.13098), (-0.35351, 1.26970, 0.45400), (-0.23789, 0.86036, 0.15354),
                 (-0.94158, 0.64552, 0.16817)], (-0.18961, 0.71147)),
         "aV": ([(-3.80595, 0.56934, 0.81061), (-3.44965, -0.22911, 0.51059), (-0.39902, 0.73042, 0.11899),
                 (0.50167, 1.07319, 0.27195)], (-0.16398, 0.63297)),
         "bH": ([(-0.14318, 1.82442, -0.55187), (0.29591, 0.77564, 0.19822), (0.32177, 0.63773, 0.13164),
                 (-5.37610, -0.96230, 1.47828), (16.1721, -3.29980, 3.43990)], (0.67849, -1.95537)),
         "bV": ([(-0.07771, 2.33840, -0.76284), (0.56727, 0.95545, 0.54039), (-0.20238, 1.14520, 0.26809),
                 (-48.2991, 0.791669, 0.116226), (48.5833, 0.791459, 0.116479)], (-0.053739, 0.83433))}


def itu_p838(f_ghz, polarization):
    """(a, b) of the power law k = a R^b (dB/km) of ITU-R P.838-3 for one frequency (GHz) and a polarization that
    starts with "v" or "h" (any case). Used only as a reference curve (the expected path-rain attenuation per km)."""
    x, p = np.log10(float(f_ghz)), "V" if str(polarization).strip().lower().startswith("v") else "H"
    val = {}
    for name in ("a", "b"):
        terms, (m, c) = _P838[name + p]
        val[name] = sum(aj * np.exp(-((x - bj) / cj) ** 2) for aj, bj, cj in terms) + m * x + c
    return float(10 ** val["a"]), float(val["b"])
