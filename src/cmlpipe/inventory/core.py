"""Step 7: event significance, wet runs, inventory metrics, TSL behaviour and quality status (pure functions).

Everything here is ADDITIVE: it reads the step 1-6 products and writes new tables; events.csv, the split
files and L1/L2/L3 are never changed. The stages are in cmlpipe.inventory.stage.

Significance (config event_significance, physics-based, NOT tuned on CML)
  significant = city_total_mean_mm >= 1 OR radar_domain_total_mm >= 1 OR city_max_5min_rate_mmh >= 5
  radar_domain_total_mm (events.csv) is the DOMAIN-MEAN radar accumulation of the event: the 5 min mean floored
  rain rate over the 720 domain cells, summed over the event bins x 5/60 h, i.e. the average rain depth over the
  domain, not the depth at any one place.
  rain_class: none_drizzle (not significant), significant, heavy (significant and total >= 10 mm or rate >= 30 mm/h).

Wet runs: the wet bins of the network reference (events stage, same wet definition) BEFORE gap merging. A run is
significant by the same rule with both totals lowered to run_total_mm; the event's wet core spans its first to its
last significant run.

Quality status (config inventory.status), checked in this order:
  unusable: JJA tl completeness < 50 %, or > 20 % of the samples present but blocking-flagged, or rx_floor > 20 %
  degraded: any month tl completeness < 80 %, rx_floor > 5 %, a spike / flatline / tsl_step share above its cap,
            dry tl noise (std of the dry 5 min means) above the network p95, or dry baseline range above the p95
  good:     otherwise. The reason string lists every rule that fired.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from cmlpipe import events
from cmlpipe.io import schema

RAIN_CLASSES = ("none_drizzle", "significant", "heavy")
TSL_CLASSES = ("no_tsl", "toggling", "constant", "rare_steps", "frequent_steps")


# ---------------------------------------------------------------------------
# Significance and runs
# ---------------------------------------------------------------------------
def significance(et: pd.DataFrame, scfg: dict) -> pd.DataFrame:
    """event_id, significant, sig_reason, rain_class and the three numbers used."""
    tot = et.city_total_mean_mm.fillna(0.0)
    rad = et.radar_domain_total_mm.fillna(0.0)
    rate = et.city_max_5min_rate_mmh.fillna(0.0)
    c1 = tot >= scfg["city_total_mean_mm"]
    c2 = rad >= scfg["radar_domain_total_mm"]
    c3 = rate >= scfg["city_max_5min_rate_mmh"]
    sig = c1 | c2 | c3
    reason = [", ".join(n for n, f in (("city_total", a), ("radar_total", b), ("city_rate", c)) if f) or "below all"
              for a, b, c in zip(c1, c2, c3)]
    heavy = sig & ((tot >= scfg["heavy_total_mm"]) | (rate >= scfg["heavy_rate_mmh"]))
    return pd.DataFrame({
        "event_id": et.event_id, "significant": sig.to_numpy(), "sig_reason": reason,
        "rain_class": np.where(heavy, "heavy", np.where(sig, "significant", "none_drizzle")),
        "city_total_mean_mm": et.city_total_mean_mm, "radar_domain_total_mm": et.radar_domain_total_mm,
        "city_max_5min_rate_mmh": et.city_max_5min_rate_mmh,
    })


def event_runs(times: pd.DatetimeIndex, state: np.ndarray, event_id: np.ndarray, domain_mean: np.ndarray,
               waf: np.ndarray, max_rate: np.ndarray, city_amount: np.ndarray, step: pd.Timedelta,
               scfg: dict) -> pd.DataFrame:
    """The wet runs (before gap merging) inside events, one row each.

    times, state, event_id, domain_mean, waf, max_rate : per bin of the reference series.
    city_amount : (time, gauge) 5 min City amounts (NaN = gauge unknown).
    """
    h = step / pd.Timedelta("1h")
    rows = []
    prev_end = {}
    for a, b in events.runs(state == events.STATE_CODES["wet"]):
        eid = event_id[a]
        if not eid:  # a short run that was dropped (no event)
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            tots = np.where(np.isfinite(city_amount[a:b]).any(axis=0), np.nansum(city_amount[a:b], axis=0), np.nan)
            city_mean = float(np.nanmean(tots))
            rad_tot = float(np.nansum(domain_mean[a:b]) * h)
            rate = float(np.nanmax(max_rate[a:b])) if np.isfinite(max_rate[a:b]).any() else 0.0
        start, end = times[a], times[b - 1] + step
        sig = (np.nan_to_num(city_mean) >= scfg["run_total_mm"]) or (rad_tot >= scfg["run_total_mm"]) or \
              (rate >= scfg["city_max_5min_rate_mmh"])
        rows.append({"event_id": eid, "start": start, "end": end, "duration_min": (end - start) / pd.Timedelta("1min"),
                     "radar_peak_domain_rain_mmh": float(np.nanmax(domain_mean[a:b])) if np.isfinite(domain_mean[a:b]).any() else np.nan,
                     "radar_max_wet_area_fraction": float(np.nanmax(waf[a:b])) if np.isfinite(waf[a:b]).any() else np.nan,
                     "radar_domain_total_mm": rad_tot, "city_total_mean_mm": city_mean, "city_max_5min_rate_mmh": rate,
                     "gap_to_previous_run_h": ((start - prev_end[eid]) / pd.Timedelta("1h")) if eid in prev_end else np.nan,
                     "significant_run": bool(sig)})
        prev_end[eid] = end
    df = pd.DataFrame(rows)
    if len(df):
        df.insert(0, "run_id", [f"{e}_R{k:02d}" for e, k in zip(df.event_id, df.groupby("event_id").cumcount() + 1)])
    return df


def wet_core(runs: pd.DataFrame) -> pd.DataFrame:
    """Per event: n_runs, wet_core_start / end (first / last significant run), wet_core_hours."""
    rows = []
    for eid, g in runs.groupby("event_id", sort=True):
        s = g[g.significant_run]
        rows.append({"event_id": eid, "n_runs": len(g),
                     "wet_core_start": s.start.min() if len(s) else pd.NaT, "wet_core_end": s.end.max() if len(s) else pd.NaT,
                     "wet_core_hours": ((s.end.max() - s.start.min()) / pd.Timedelta("1h")) if len(s) else np.nan})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Completeness, flags, quantization
# ---------------------------------------------------------------------------
FLAG_COLS = {name: f"n_{name}" for name in schema.FLAGS}


def completeness(daily: pd.DataFrame) -> pd.DataFrame:
    """Per (cml_id, sublink_id) from qc_daily_summary.csv: valid shares (JJA and per month), not-blocked share,
    per-flag shares. n_valid_tsl = n_samples - n_missing_tsl (tsl is missing exactly when the bit is set)."""
    d = daily.copy()
    d["n_valid_tsl"] = d.n_samples - d.n_missing_tsl
    d["n_not_blocked"] = d.n_samples - d.n_blocking
    d["month"] = pd.to_datetime(d.day).dt.month
    keys = ["cml_id", "sublink_id"]
    s = d.groupby(keys)[["n_samples", "n_valid_rsl", "n_valid_tsl", "n_valid_tl", "n_not_blocked", *FLAG_COLS.values()]].sum()
    out = pd.DataFrame(index=s.index)
    for v in ("rsl", "tsl", "tl"):
        out[f"completeness_{v}"] = s[f"n_valid_{v}"] / s.n_samples
    out["share_not_blocked"] = s.n_not_blocked / s.n_samples
    # Samples with valid rsl AND tsl that a blocking flag still removes (rx_floor, spike, out_of_range). Not-blocked
    # samples are a subset of valid-tl samples, so this is exact; missing data is NOT counted twice here.
    out["flag_share_blocking_present"] = (s.n_valid_tl - s.n_not_blocked).clip(lower=0) / s.n_samples
    for name, col in FLAG_COLS.items():
        out[f"flag_share_{name}"] = s[col] / s.n_samples
    m = d.groupby(keys + ["month"])[["n_samples", "n_valid_tl"]].sum()
    mt = (m.n_valid_tl / m.n_samples).unstack("month").reindex(columns=[6, 7, 8])  # a month without data -> NaN
    for mo, name in ((6, "jun"), (7, "jul"), (8, "aug")):
        out[f"completeness_tl_{name}"] = mt[mo]
    return out.reset_index()


def quantization(hist_rsl: np.ndarray, hist_tsl: np.ndarray, bin_db: float) -> pd.DataFrame:
    """From the period histograms (column, 0.1 dB bins): smallest rsl step, distinct rsl / tsl values, TSL top share."""
    rows = []
    for hr, ht in zip(hist_rsl, hist_tsl):
        occ = np.flatnonzero(hr)
        step = float(np.diff(occ).min() * bin_db) if occ.size > 1 else np.nan
        tot = ht.sum()
        rows.append({"rsl_min_step_db": round(step, 3) if np.isfinite(step) else np.nan, "rsl_n_distinct": int(occ.size),
                     "tsl_n_distinct_values": int((ht > 0).sum()), "tsl_top_share": float(ht.max() / tot) if tot else np.nan})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# TSL behaviour
# ---------------------------------------------------------------------------
def tsl_behaviour(n_valid_tsl, n_steps, top_share, tcfg: dict, samples_per_day: int = 8640) -> np.ndarray:
    """Class per sublink, rules applied in TSL_CLASSES order (first match wins, except no_tsl first)."""
    n_valid_tsl, n_steps, top_share = (np.asarray(x, dtype="float64") for x in (n_valid_tsl, n_steps, top_share))
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = n_steps / n_valid_tsl
        per_day = n_steps / (n_valid_tsl / samples_per_day)
    out = np.where(frac > tcfg["toggling_step_fraction"], "toggling",
          np.where((top_share >= tcfg["constant_top_share"]) & (per_day <= tcfg["constant_max_steps_per_day"]), "constant",
          np.where(per_day <= tcfg["rare_max_steps_per_day"], "rare_steps", "frequent_steps")))
    return np.where(n_valid_tsl > 0, out, "no_tsl")


# ---------------------------------------------------------------------------
# Quality status
# ---------------------------------------------------------------------------
def quality_status(inv: pd.DataFrame, scfg: dict) -> pd.DataFrame:
    """status (good / degraded / unusable) and reason for every row of the sublink inventory.

    Network quantile thresholds (dry noise, baseline range) are computed over the sublinks that are not unusable.
    Returns (status, reason) columns plus the thresholds used (in .attrs['thresholds']).
    """
    pct = 100.0
    reasons_u = [
        (inv.completeness_tl.fillna(0) * pct < scfg["unusable_completeness_tl_pct"], "JJA tl completeness < {:g} %".format(scfg["unusable_completeness_tl_pct"])),
        (inv.flag_share_blocking_present.fillna(0) * pct > scfg["unusable_flagged_pct"],
         "blocking-flagged present samples > {:g} %".format(scfg["unusable_flagged_pct"])),
        (inv.flag_share_rx_floor.fillna(0) * pct > scfg["unusable_rx_floor_pct"], "rx_floor > {:g} %".format(scfg["unusable_rx_floor_pct"])),
    ]
    unusable = np.zeros(len(inv), bool)
    for m, _ in reasons_u:
        unusable |= m.to_numpy()
    ok = ~unusable
    q_noise = float(np.nanquantile(inv.dry_tl_std[ok], scfg["degraded_dry_noise_quantile"]))
    q_base = float(np.nanquantile(inv.dry_baseline_range_db[ok], scfg["degraded_baseline_range_quantile"]))
    month_min = inv[["completeness_tl_jun", "completeness_tl_jul", "completeness_tl_aug"]].min(axis=1)
    reasons_d = [
        (month_min.fillna(0) * pct < scfg["degraded_month_completeness_tl_pct"],
         "a month with tl completeness < {:g} %".format(scfg["degraded_month_completeness_tl_pct"])),
        (inv.flag_share_rx_floor.fillna(0) * pct > scfg["degraded_rx_floor_pct"], "rx_floor > {:g} %".format(scfg["degraded_rx_floor_pct"])),
        *[(inv[f"flag_share_{f}"].fillna(0) * pct > cap, f"{f} > {cap:g} %") for f, cap in scfg["degraded_flag_caps_pct"].items()],
        (inv.dry_tl_std > q_noise, f"dry tl noise > network p{100 * scfg['degraded_dry_noise_quantile']:.0f} ({q_noise:.2f} dB)"),
        (inv.dry_baseline_range_db > q_base,
         f"dry baseline range > network p{100 * scfg['degraded_baseline_range_quantile']:.0f} ({q_base:.2f} dB)"),
    ]
    status, reason = [], []
    for i in range(len(inv)):
        ru = [txt for m, txt in reasons_u if bool(m.iloc[i])]
        if ru:
            status.append("unusable")
            reason.append("; ".join(ru))
            continue
        rd = [txt for m, txt in reasons_d if bool(m.iloc[i])]
        status.append("degraded" if rd else "good")
        reason.append("; ".join(rd) if rd else "ok")
    out = pd.DataFrame({"status": status, "reason": reason}, index=inv.index)
    out.attrs["thresholds"] = {"dry_noise_db": q_noise, "dry_baseline_range_db": q_base}
    return out


def length_class(length_m, edges_km):
    """Length class labels for an array of lengths (m)."""
    return [events.length_class(v, edges_km) for v in np.asarray(length_m, float)]
