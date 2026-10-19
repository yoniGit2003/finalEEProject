"""Baseline estimation per sub-link and rain event (functions; scripts: run_sample.py, validation/).

A(t) = TL(t) - A_baseline(t), TL = TSL - RSL. The baseline is the TL the sub-link would show without rain.
It is estimated ONLY from dry samples, never from samples during or right after rain, so it cannot follow a
drying tail (WAA). Two versions (switch: `baseline.version` in research/config.yaml; full text: METHOD.md):

v2 (current): flat baseline = MEDIAN of the usable dry 5 min tl_mean in the 24 h before the event start, or in
the 48 h before if the 24 h window has < 1 h; else missing. Optional humidity exclusion of dew-risk hours.
Functions: baseline_v2, flags_v2, dew_risk, sun_elevation_deg (positions on the regular 5 min grid).

v1 (first version, kept for comparison; the rules below; parameters in `baseline` / `baseline.v1` of the config):
  usable dry sample : 5 min bin with link_state == dry (radar + gauges, never CML) AND at least X hours after the
                      last wet bin of this link (X = 1 h). The +-30 min "unknown" margins are not dry, so excluded.
  window means      : mean TL of the usable dry samples in the WINDOW_H hours before the event start ("before")
                      and after the event end ("after"); a window with < MIN_DRY_H hours of samples is missing.
  usual level       : median TL of all usable dry samples of the season (what the baseline "normally" is).
  excursion         : a window whose mean is more than EXC_DB from the usual level (an extreme dry-time shift).
  method A          : constant baseline = before.                                     (default)
  method B          : straight line from before (at event start) to after (at event end), used only when
                      |after - before| > DIFF_DB and neither window is an excursion.
  missing           : no usable before window -> no baseline (no guessing).
"""
import numpy as np
import pandas as pd

STEP = pd.Timedelta("5min")


def load_sublinks(d, sublink_ids, start, end):
    """5 min TL and its number of valid 10 s samples for the given sub-links, plus the link state of their links.

    Returns {sublink_id: DataFrame(index=time, columns=[tl, n, state])}. tl is the mean of the valid (QC-passed)
    10 s TL samples in the bin; state is 0 dry, 1 wet, -1 unknown.
    """
    ds = d.cml(sublink_ids, start, end, res="5min", variables=["tl_mean", "n_valid_tl"])
    cml_ids = sorted({int(s.split("_")[0]) for s in sublink_ids})
    state = d.link_state(cml_ids, start, end).state.to_pandas()          # (time, cml_id)
    out = {}
    for s in sublink_ids:
        df = pd.DataFrame({"tl": ds.tl_mean.sel(sublink_id=s).to_series(),
                           "n": ds.n_valid_tl.sel(sublink_id=s).to_series()})
        df["state"] = state[int(s.split("_")[0])].reindex(df.index)
        out[s] = df
    return out


def usable_dry(state, x_hours):
    """True where the bin is dry AND at least x_hours after the last wet bin (bins before any wet bin count)."""
    t = pd.Series(state.index, index=state.index)
    last_wet = t.where(state == 1).ffill()            # time of the most recent wet bin (NaT if none yet)
    since = t - last_wet - STEP                       # time from the END of that wet bin to this bin's start
    far = since.isna() | (since >= pd.Timedelta(hours=x_hours))
    return (state == 0) & far


def window_stats(df, mask, t0, t1):
    """Mean TL of the usable samples in [t0, t1), weighted by the number of 10 s samples (= mean of the 10 s TL).

    Returns (mean, hours of samples, std of the 5 min means). The mean, not the median, keeps the information of
    how long the quantized TL sits on each level, which places the estimate between the quantization levels.
    """
    w = df.loc[t0:t1 - pd.Timedelta("1ns")]
    m = mask.loc[w.index] & w.tl.notna() & (w.n > 0)
    if not m.any():
        return np.nan, 0.0, np.nan
    tl, n = w.tl[m], w.n[m]
    return float(np.average(tl, weights=n)), m.sum() * STEP / pd.Timedelta("1h"), float(tl.std())


def estimate_event(df, mask, usual, ev_start, ev_end, cfg):
    """Baseline of one sub-link for one event. Returns a dict (one row of the results table)."""
    win = pd.Timedelta(hours=cfg["WINDOW_H"])
    before, h_b, sd_b = window_stats(df, mask, ev_start - win, ev_start)
    after, h_a, sd_a = window_stats(df, mask, ev_end, ev_end + win)
    ok_b, ok_a = h_b >= cfg["MIN_DRY_H"], h_a >= cfg["MIN_DRY_H"]
    exc_b = ok_b and abs(before - usual) > cfg["EXC_DB"]
    exc_a = ok_a and abs(after - usual) > cfg["EXC_DB"]

    flags = []
    if not ok_b:
        method = "missing"
        flags.append("no_before")
    elif not ok_a:
        method = "A"
        flags.append("no_after")
    elif exc_b or exc_a:
        method = "A"                                   # extreme dry shift in a window: B is not allowed
        flags += ["excursion_before"] * exc_b + ["excursion_after"] * exc_a
    elif abs(after - before) > cfg["DIFF_DB"]:
        method = "B"
    else:
        method = "A"

    row = {"event_start": ev_start, "event_end": ev_end, "usual_level_db": usual,
           "before_db": before, "before_hours": h_b, "before_std_db": sd_b,
           "after_db": after, "after_hours": h_a, "after_std_db": sd_a,
           "closure_db": after - before if ok_a and ok_b else np.nan,   # how much the dry level moved over the event
           "method": method, "flags": ";".join(flags)}

    # out-of-sample check: usable dry bins INSIDE the event were not used to estimate the baseline, so A there
    # should be ~0. A clearly positive value can also be a drying tail longer than X (see README).
    if method != "missing":
        inside = df.loc[ev_start:ev_end - pd.Timedelta("1ns")]
        a = (inside.tl - baseline_series(row, inside.index))[mask.loc[inside.index]].dropna()
        row["in_event_dry_hours"] = len(a) * STEP / pd.Timedelta("1h")
        row["in_event_dry_A_mean_db"] = float(a.mean()) if len(a) else np.nan
        row["A_max_db"] = float((inside.tl - baseline_series(row, inside.index)).max())
    return row


def baseline_series(row, index):
    """The baseline of one event (row from estimate_event) on a time index (NaN outside the event)."""
    t0, t1 = pd.Timestamp(row["event_start"]), pd.Timestamp(row["event_end"])
    inside = (index >= t0) & (index < t1)
    out = pd.Series(np.nan, index=index)
    if row["method"] == "A":
        out[inside] = row["before_db"]
    elif row["method"] == "B":
        frac = (index[inside] - t0) / (t1 - t0)                       # 0 at the event start, 1 at the end
        out[inside] = row["before_db"] + frac * (row["after_db"] - row["before_db"])
    return out


# =================================================================================================== v2
SPH = 12                                                   # 5 min slots per hour


def dew_risk(rh_pct, dpd_c, dew_cfg):
    """True where dew can form: RH >= rh_pct OR dew-point depression (T - Td) < dpd_c."""
    rh, dpd = np.asarray(rh_pct, float), np.asarray(dpd_c, float)
    return (rh >= dew_cfg["rh_pct"]) | (dpd < dew_cfg["dpd_c"])


def dew_point(t_c, rh_pct):
    """Dew point (deg C) from temperature and relative humidity, Magnus formula (b = 17.62, c = 243.12 deg C)."""
    b, c = 17.62, 243.12
    g = np.log(np.clip(np.asarray(rh_pct, float), 1, 100) / 100) + b * np.asarray(t_c, float) / (c + np.asarray(t_c, float))
    return c * g / (b - g)


def sun_elevation_deg(times_utc, lat, lon):
    """Solar elevation angle (deg) at the given UTC times (NOAA approximation, accuracy ~0.5 deg)."""
    t = pd.DatetimeIndex(times_utc)
    hr = t.hour.to_numpy() + t.minute.to_numpy() / 60
    g = 2 * np.pi / 365 * (t.dayofyear.to_numpy() - 1 + (hr - 12) / 24)        # fractional year (rad)
    decl = (0.006918 - 0.399912 * np.cos(g) + 0.070257 * np.sin(g) - 0.006758 * np.cos(2 * g)
            + 0.000907 * np.sin(2 * g) - 0.002697 * np.cos(3 * g) + 0.00148 * np.sin(3 * g))
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(g) - 0.032077 * np.sin(g) - 0.014615 * np.cos(2 * g)
                    - 0.040849 * np.sin(2 * g))                                # equation of time (min)
    ha = np.deg2rad((hr * 60 + eqt + 4 * lon) / 4 - 180)                       # hour angle
    la = np.deg2rad(lat)
    cosz = np.sin(la) * np.sin(decl) + np.cos(la) * np.cos(decl) * np.cos(ha)
    return 90 - np.rad2deg(np.arccos(np.clip(cosz, -1, 1)))


def usable_mask(state, skip_h, dew=None):
    """Usable dry slots (numpy, 5 min grid): link state dry, at least skip_h after the end of the last wet slot,
    and, if `dew` (bool array of dew-risk slots) is given, not at dew risk (the v2 humidity exclusion)."""
    st = np.asarray(state)
    wet = st == 1
    pos = np.arange(len(st))
    last = np.maximum.accumulate(np.where(wet, pos, -1))
    since = np.where(last < 0, np.inf, (pos - last - 1) / SPH)
    ok = (st == 0) & (since >= skip_h)
    return ok & ~dew if dew is not None else ok


def _win(tl, ok, i0, n):
    a = max(i0 - n, 0)
    v, m = tl[a:i0], ok[a:i0]
    return v[m & ~np.isnan(v)], a


def baseline_v2(tl, ok, i0, B):
    """v2 baseline for an event starting at slot i0: (value, window hours used, first slot of the window).
    Median of the usable dry slots in the 24 h before i0; 48 h if fewer than min_dry_h; else (NaN, 0, i0)."""
    min_n = int(round(B["min_dry_h"] * SPH))
    for w in (B["v2"]["window_h"], B["v2"]["fallback_window_h"]):
        v, a = _win(tl, ok, i0, int(w * SPH))
        if len(v) >= min_n:
            return float(np.median(v)), w, a
    return np.nan, 0.0, i0


def flags_v2(value, window_first, i0, state, ok_plain, dew, summer_median, cycle_amp, is_degraded, is_toggling, B,
             event_end_pos=None):
    """The v2 flags of one (event, sub-link). ok_plain = usable mask WITHOUT humidity exclusion (the window's dry
    slots); dew = dew-risk slots (5 min, hourly values); cycle_amp = the sub-link's daily-cycle amplitude (R1 C3);
    event_end_pos = end of the link's last wet slot of the event (slot index; added in R1c): the event's time is
    [i0, event_end_pos + 1 h), and dew_risk_during_event is set if more than dew.window_share of it is dew risk."""
    dry = ok_plain[window_first:i0]
    share = float(dew[window_first:i0][dry].mean()) if dry.any() else np.nan
    ev_share = np.nan
    if event_end_pos is not None:
        span = dew[i0:min(event_end_pos + SPH, len(dew))]
        ev_share = float(span.mean()) if len(span) else np.nan
    return {"dew_risk_during_event": bool(ev_share > B["dew"]["window_share"]) if not np.isnan(ev_share) else False,
            "dew_share_during_event": ev_share,
            "excursion": bool(abs(value - summer_median) > B["excursion_db"]) if not np.isnan(value) else False,
            "dew_risk_in_window": bool(share > B["dew"]["window_share"]) if not np.isnan(share) else False,
            "dew_share_in_window": share,
            "daily_cycle_large": bool(cycle_amp > B["daily_cycle_large_db"]) if not np.isnan(cycle_amp) else False,
            "degraded": bool(is_degraded), "toggling": bool(is_toggling)}


# =================================================================================================== products
from research.baseline import paths as P                  # all file locations (this folder is on sys.path when baseline is imported)


def dew_hourly():
    """outputs/products/dew_risk_hourly.csv (built by make_dew_products.py), indexed by the UTC hour start."""
    return pd.read_csv(P.DEW_RISK_HOURLY, parse_dates=["time_utc", "time_local"]).set_index("time_utc")


def dew_mask_5min(index):
    """Dew-risk flag on a 5 min time index: each slot carries the dew_risk of its hour (False where unknown)."""
    h = dew_hourly().dew_risk
    return h.reindex(pd.DatetimeIndex(index).floor("1h")).fillna(False).to_numpy(bool)


def _open_release():
    import yaml
    from cmlpipe.api import CMLData
    cfg = yaml.safe_load(open(P.CONFIG_FILE, encoding="utf-8"))
    return CMLData.open(P.REPO / "data" / "release" / cfg["release"])


_TABLE = None


def a_of_t(cml_id, event_id, res="5min", pad="0h", d=None):
    """A(t) = tl_mean(t) - baseline for both sub-links of one link during one event, at any CML resolution
    ("1min", "5min", ...). The baseline comes from outputs/products/baselines.csv (A(t) itself is never stored).

    Returns a DataFrame (index: time UTC, interval start) with one column per sub-link (A in dB, NaN if the pair has no
    baseline) over [event start - pad, event end + pad). `d` = an open CMLData (opened from the config if None).
    Example:  A = a_of_t(10227, "E20150727_2140", res="1min", pad="6h")
    """
    global _TABLE
    if _TABLE is None:
        _TABLE = pd.read_csv(P.BASELINES, parse_dates=["event_start", "event_end"])
    rows = _TABLE[(_TABLE.cml_id == int(cml_id)) & (_TABLE.event_id == event_id)]
    if rows.empty:
        raise KeyError(f"no baseline rows for cml_id {cml_id} and event {event_id} in {P.BASELINES.name}")
    d = d or _open_release()
    t0 = rows.event_start.iloc[0] - pd.Timedelta(pad)
    t1 = rows.event_end.iloc[0] + pd.Timedelta(pad)
    p = d.catalog["period"]
    t0, t1 = max(t0, pd.Timestamp(p["start"])), min(t1, pd.Timestamp(p["end_exclusive"]))
    tl = d.cml(rows.sublink_id.tolist(), t0, t1, res=res, variables=["tl_mean"]).tl_mean.to_pandas()
    base = rows.set_index("sublink_id").baseline_db
    out = tl[base.index] - base                                   # broadcast: one baseline per sub-link column
    out.columns.name = "sublink_id"
    out.attrs.update({"cml_id": int(cml_id), "event_id": event_id, "res": res, "unit": "dB",
                      "definition": "A(t) = tl_mean(t) - baseline v2 (outputs/products/baselines.csv)"})
    return out


# =================================================================================================== the baseline table (moved from make_baseline_table.py)
def table_inputs(d, link_splits=("train",)):
    """Sub-links of the link splits (good or degraded), their 5 min tl_mean, the per-link state, dew-risk mask, daily-cycle amp.
    (Moved from make_baseline_table.load; step S8 also asks for the TEST links.)"""
    s = d.sublinks()
    s = s[s.split_link.isin(list(link_splits))].copy()
    s["status"] = s.sublink_id.map(d.status().set_index("sublink_id").status)
    s = s[s.status.isin(["good", "degraded"])].reset_index(drop=True)
    p = d.catalog["period"]
    tl = d.cml(s.sublink_id.tolist(), p["start"], p["end_exclusive"], res="5min",
               variables=["tl_mean"]).tl_mean.to_pandas()
    state = d.store("link_reference_state_5min").state.load().to_pandas().reindex(tl.index)
    dew = dew_mask_5min(tl.index)
    amp = pd.read_csv(P.C3_DAILY_CYCLE).set_index("sublink_id").amplitude_db
    return s, tl, state, dew, amp


def table_masks(state_col, dew, B):
    ok_plain = usable_mask(state_col.to_numpy(), B["skip_h"])
    ok = usable_mask(state_col.to_numpy(), B["skip_h"], dew) if B["v2"]["humidity_exclusion"] else ok_plain
    return ok_plain, ok


def event_table(d, s, tl, state, dew, amp, unc, B, splits=("train", "val")):
    """One row per (event, sub-link) of the event splits: the v2 baseline and its flags (moved from make_baseline_table.table;
    step S8 also asks for the TEST events)."""
    ev = d.events()
    ev = ev[ev.split_event.isin(list(splits))]
    T = tl.index
    u = unc.set_index("sublink_id")
    rows = []
    for r in s.itertuples():
        c = r.cml_id
        x, st = tl[r.sublink_id].to_numpy(), state[c].to_numpy()
        ok_plain, ok = table_masks(state[c], dew, B)
        summer = float(np.median(x[ok_plain & ~np.isnan(x)]))
        wet = st == 1
        for e in ev.itertuples():
            i0, i1 = T.searchsorted(pd.Timestamp(e.start)), T.searchsorted(pd.Timestamp(e.end))
            value, w, a = baseline_v2(x, ok, i0, B)
            n_use = int((ok[a:i0] & ~np.isnan(x[a:i0])).sum())
            wp = np.nonzero(wet[i0:i1])[0]
            end_pos = i0 + wp[-1] + 1 if len(wp) else i1          # no wet slot of the link: the event end
            fl = flags_v2(value, a, i0, st, ok_plain, dew, summer, amp.get(r.sublink_id, np.nan),
                             r.status == "degraded", r.tsl_behaviour == "toggling", B, event_end_pos=end_pos)
            rows.append({"event_id": e.event_id, "split_event": e.split_event, "significant": e.significant,
                         "run_count": e.n_runs, "event_start": e.start, "event_end": e.end,
                         "link_has_wet_slot": bool(len(wp)), "last_wet_end": T[min(end_pos, len(T) - 1)],
                         "sublink_id": r.sublink_id, "cml_id": c, "status": r.status,
                         "baseline_db": value, "window_h": w, "n_usable_slots": n_use,
                         "missing": bool(np.isnan(value)), "summer_dry_median_db": summer, **fl,
                         "offset_std": u.offset_std.get(r.sublink_id, np.nan),
                         "offset_p90_abs": u.offset_p90_abs.get(r.sublink_id, np.nan)})
    t = pd.DataFrame(rows)
    t["use_in_waa_main"] = ~(t.missing | t.excursion | t.daily_cycle_large | t.degraded | t.toggling)
    return t
