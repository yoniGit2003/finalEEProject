"""WAA_estm step 1: extract ALL drying tails with their features. No selection, no fitting, no shape filter.

In a drying tail the path rain is over (A_rain = 0), so A(t) = δ(t) = WAA + baseline residual. This script cuts every
tail out of the data and describes them; what to do with them is decided later.

Tail (one per event x sub-link; README.md, "Tail definition"; settings: `waa.tails` in research/config.yaml):
  pairs        every row of baselines.csv (frozen baseline v2) with missing == False and link_has_wet_slot == True
               (train + val events, train links)
  A(t)         1 min tl_mean (release) - baseline_db (baselines.csv)
  t0           last_wet_end: the first non-wet link_state bin after the link's last wet bin in the event
  t_first_dry  the first link_state == dry bin at or after t0 (searched over the whole season)
  t_end        earliest of: the next WET link_state bin of the link (whole season), the start of a NaN run of A
               >= gap_min minutes, t0 + cap_h; ended_by = next_wet / gap / cap (season_end if the data stop first)

Writes (locations: WAA_estm_paths.py; every file: CATALOG.md):
  outputs/products/tails_1min.parquet   one row per tail minute
  outputs/products/tails_index.csv      one row per tail, with its features
  outputs/tables/step1_*.csv            check of A(t), return to the dry level, median curves, dry-weather rises
  outputs/figures/step1_*.png           distributions, median curves, resolution panels
  outputs/reports/step1_tails.md        the descriptive report (printed at the end)

Run from the repository root (a few minutes; reads the 1 min CML store once):
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step1_extract_tails.py
"""
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cmlpipe.api import CMLData

sys.stdout.reconfigure(encoding="utf-8")                            # the report has δ: Windows consoles are cp1252
MODULE = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE))
import WAA_estm_paths as P                                          # every location of this module
from research.common.save import save
sys.path.insert(0, str(P.BASELINE_MODULE))
import baseline as bl                                               # only for the a_of_t check (step 2 of the task)
import baseline_estm_paths as BP                                    # where the frozen baseline products are

CFG = yaml.safe_load(open(P.CONFIG_FILE, encoding="utf-8"))
TC, RC = CFG["waa"]["tails"], CFG["waa"]["report"]
CAP_MIN = int(round(TC["cap_h"] * 60))
GAP_MIN = int(TC["gap_min"])
WET, DRY, UNKNOWN = 1, 0, -1                                        # link_state codes
QC = {"spike": 16, "tsl_step": 32, "network_gap": 128, "missing": 1 | 2}   # flags_or bits (missing = rsl or tsl)
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]               # categorical slots 1-4 (fixed order)
INK, INK2 = "#0b0b0b", "#52514e"
LOCAL = pd.Timedelta(hours=2)                                       # Gothenburg summer time = UTC + 2 h
for p in (P.PRODUCTS, P.TABLES, P.FIGURES, P.REPORTS):
    p.mkdir(parents=True, exist_ok=True)
clock = time.time()

d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])


def next_true(mask):
    """For every position i: the first position j >= i where mask is True (len(mask) if there is none)."""
    n = len(mask)
    idx = np.where(mask, np.arange(n), n)
    return np.minimum.accumulate(idx[::-1])[::-1]


def frac(x):
    """Share of True in a boolean array (NaN for an empty one)."""
    return float(np.mean(x)) if len(x) else np.nan


# =================================================================================================== 1. the pairs
b = pd.read_csv(BP.BASELINES, parse_dates=["event_start", "event_end", "last_wet_end"])
pairs = b[~b.missing & b.link_has_wet_slot].sort_values(["sublink_id", "event_start"]).reset_index(drop=True)
pairs.insert(0, "tail_id", np.arange(len(pairs), dtype="int32"))
print(f"1. pairs: {len(pairs)} of {len(b)} baseline rows (missing == False and link_has_wet_slot == True)")

links = d.links().set_index("cml_id")
sub = d.sublinks().set_index("sublink_id")

# =================================================================================================== references
c1 = d.store("cml_1min")
T1 = pd.DatetimeIndex(c1.time.values)                               # 1 min grid, interval starts
ls = d.store("link_reference_state_5min")
T5 = pd.DatetimeIndex(ls.time.values)                               # 5 min grid
assert (T1[::5] == T5).all(), "1 min and 5 min grids do not line up"
n1, n5 = len(T1), len(T5)
state5 = ls.state.load().to_pandas()                                # (time, cml_id): 0 dry, 1 wet, -1 unknown
ra = d.store("radar_along_links_5min")
rmax5 = ra.R_max.load().to_pandas()                                 # mm/h, wettest cell under the path
rmean5 = ra.R_mean.load().to_pandas()                               # mm/h, path-weighted mean

# gauge: the link's nearest gauge, if within gauge_max_dist_m (5 min rate; SMHI 15 min rate on its three 5 min bins)
city = d.store("gauges_city_5min").rate_mm_h.load().to_pandas()
smhi = d.store("gauges_smhi_15min").rate_mm_h.load().to_pandas().reindex(T5, method="ffill", limit=2)
gauge_rate5 = pd.concat([city, smhi], axis=1)                       # (time, gauge_id), mm/h


def link_gauge(c):
    """(gauge_id, distance m) of the link's nearest gauge if it is within the limit, else (None, NaN)."""
    g, dist = links.nearest_gauge_id[c], links.nearest_gauge_distance_m[c]
    return (g, float(dist)) if dist <= TC["gauge_max_dist_m"] else (None, np.nan)


# dew risk and weather (hourly, baseline product): 1 min flag, and the row at t0
dew = pd.read_csv(BP.DEW_RISK_HOURLY, parse_dates=["time_utc"]).set_index("time_utc")
dew1 = dew.dew_risk.reindex(T1.floor("1h")).fillna(False).to_numpy(bool)
jumps = pd.read_csv(BP.DRY_JUMPS, parse_dates=["start_utc", "end_utc"])   # [start, end) intervals, 5 min edges

# =================================================================================================== 2-4. tails
SCHEMA = pa.schema([("tail_id", pa.int32()), ("sublink_id", pa.string()), ("cml_id", pa.int32()),
                    ("event_id", pa.string()), ("t", pa.timestamp("ns")), ("minutes_since_t0", pa.int16()),
                    ("A_db", pa.float32()), ("n_valid_tl", pa.uint8()), ("link_state", pa.int8()),
                    ("R_max_mmh", pa.float32()), ("R_mean_mmh", pa.float32()), ("gauge_rate_mmh", pa.float32()),
                    ("flag_spike", pa.bool_()), ("flag_tsl_step", pa.bool_()), ("flag_network_gap", pa.bool_()),
                    ("flag_missing", pa.bool_()), ("dew_risk", pa.bool_()), ("in_dry_jump", pa.bool_())])
writer = pq.ParquetWriter(P.TAILS_1MIN, SCHEMA, compression="zstd")
rng = np.random.default_rng(RC["seed"])
check_ids = set(rng.choice(pairs.tail_id.to_numpy(), 3, replace=False).tolist())   # pairs for the a_of_t check
check_A = {}                                                        # tail_id -> full-season A built here
index_rows, rise_rows = [], []
by_sub = {s: g for s, g in pairs.groupby("sublink_id")}
store_subs = c1.sublink_id.values.astype(str)

print("2. cutting the tails (1 min CML read in blocks of 32 sub-links = one store chunk) ...")
for a in range(0, len(store_subs), 32):                             # one chunk of the 1 min store at a time
    block = [s for s in store_subs[a:a + 32] if s in by_sub]
    if not block:
        continue
    ds = c1.isel(sublink_id=slice(a, a + 32))[["tl_mean", "n_valid_tl", "flags_or"]].load()
    parts = []
    for s in block:
        col = list(store_subs[a:a + 32]).index(s)
        tl = ds.tl_mean.transpose("time", "sublink_id").values[:, col].astype("float64")   # exact copy of float32
        nv = ds.n_valid_tl.transpose("time", "sublink_id").values[:, col]
        fl = ds.flags_or.transpose("time", "sublink_id").values[:, col].astype("int64")
        c = int(s.split("_")[0])
        st5 = state5[c].to_numpy()
        st1 = np.repeat(st5, 5)                                     # 5 min state on the 1 min grid
        next_wet5, next_dry5 = next_true(st5 == WET), next_true(st5 == DRY)
        isnan = np.isnan(tl)
        nan_fwd = next_true(~isnan) - np.arange(n1)                 # length of the NaN run starting at each minute
        next_gap = next_true(nan_fwd >= GAP_MIN)                    # first minute that starts a NaN run >= gap_min
        rmax, rmean = rmax5[c].to_numpy(), rmean5[c].to_numpy()
        g_id, g_dist = link_gauge(c)
        grate5 = gauge_rate5[g_id].to_numpy() if g_id else np.full(n5, np.nan)
        jump1 = np.zeros(n1, bool)
        for j in jumps[jumps.sublink_id == s].itertuples():
            jump1[T1.searchsorted(j.start_utc):T1.searchsorted(j.end_utc)] = True
        meta = sub.loc[s]

        for r in by_sub[s].itertuples():
            A = tl - r.baseline_db                                  # A(t) = tl_mean - baseline_db, float64
            if r.tail_id in check_ids:
                check_A[r.tail_id] = A.copy()
            i05 = T5.searchsorted(r.last_wet_end)                   # t0 on the 5 min grid ...
            assert T5[i05] == r.last_wet_end
            i0 = i05 * 5                                            # ... and on the 1 min grid
            # --- t_end: the earliest stop
            stops = {"next_wet": next_wet5[i05] * 5 if next_wet5[i05] < n5 else n1 + 1,
                     "gap": next_gap[i0] if next_gap[i0] < n1 else n1 + 1,
                     "cap": i0 + CAP_MIN}
            ended_by = min(stops, key=lambda k: stops[k])           # ties: next_wet before gap before cap
            i_end = stops[ended_by]
            if i_end > n1 or (i05 == n5 - 1 and st5[i05] == WET):  # the season ends first (or t0 was clipped to
                ended_by, i_end = "season_end", min(i_end, n1)     # the last bin of the season while still wet)
            fd5 = next_dry5[i05]
            t_first_dry = T5[fd5] if fd5 < n5 else pd.NaT
            seg = slice(i0, i_end)
            m = np.arange(i_end - i0)                               # minutes since t0
            st_seg = st1[seg]

            # --- wet phase of this link in the event: first wet bin of the link in the event .. t0
            e05 = T5.searchsorted(r.event_start)
            wet_pos = np.nonzero(st5[e05:i05] == WET)[0]
            w05 = e05 + wet_pos[0] if len(wet_pos) else i05
            wseg5 = slice(w05, i05)
            gtot = np.nansum(grate5[wseg5]) * 5 / 60 if g_id and np.isfinite(grate5[wseg5]).any() else np.nan
            a_wet = A[w05 * 5:i0]
            # --- margin: the 5 min bins from t0 to the first dry bin
            mseg5 = slice(i05, fd5 if fd5 < n5 else n5)
            wx = dew.loc[r.last_wet_end.floor("1h")] if r.last_wet_end.floor("1h") in dew.index else None
            a0 = A[i0:i0 + 5]

            index_rows.append({
                "tail_id": r.tail_id, "event_id": r.event_id, "sublink_id": s, "cml_id": c,
                "split_event": r.split_event, "frequency_ghz": meta.frequency_ghz, "polarization": meta.polarization,
                "length_m": meta.length_m, "length_class": meta.length_class, "frequency_band": meta.frequency_band,
                "site_0_antenna_diameter_m": meta.site_0_antenna_diameter_m,
                "site_1_antenna_diameter_m": meta.site_1_antenna_diameter_m,
                "gauge_id": g_id or "", "gauge_distance_m": g_dist,
                "t0": r.last_wet_end, "t_first_dry": t_first_dry, "t_end": T1[i_end - 1] + pd.Timedelta("1min")
                if i_end > i0 else r.last_wet_end, "length_min": i_end - i0, "ended_by": ended_by,
                "first_dry_after_end": bool(fd5 * 5 >= i_end), "state_at_t0": int(st5[i05]),
                "wet_start": T5[w05], "wet_span_min": (i05 - w05) * 5, "wet_bins_min": len(wet_pos) * 5,
                "wet_radar_R_mean_sum_mm": np.nansum(rmean[wseg5]) * 5 / 60,
                "wet_radar_R_max_max_mmh": np.nanmax(rmax[wseg5]) if np.isfinite(rmax[wseg5]).any() else np.nan,
                "wet_gauge_total_mm": gtot,
                "wet_A_max_db": np.nanmax(a_wet) if np.isfinite(a_wet).any() else np.nan,
                "A_start_db": np.nanmean(a0) if np.isfinite(a0).any() else np.nan,
                "margin_min": (mseg5.stop - i05) * 5,
                "margin_R_max_max_mmh": np.nanmax(rmax[mseg5]) if np.isfinite(rmax[mseg5]).any() else np.nan,
                "margin_gauge_rate_max_mmh": np.nanmax(grate5[mseg5]) if g_id and np.isfinite(grate5[mseg5]).any()
                else np.nan,
                "frac_dew_risk": frac(dew1[seg]), "frac_dry_jump": frac(jump1[seg]),
                "frac_tsl_step": frac((fl[seg] & QC["tsl_step"]) > 0), "frac_spike": frac((fl[seg] & QC["spike"]) > 0),
                "frac_nan": frac(isnan[seg]), "frac_unknown": frac(st_seg == UNKNOWN),
                "t0_temperature_c": wx.temperature if wx is not None else np.nan,
                "t0_humidity_pct": wx.humidity if wx is not None else np.nan,
                "t0_wind_speed_ms": wx.wind_speed if wx is not None else np.nan,
                "t0_dew_point_depression_c": wx.dew_point_depression if wx is not None else np.nan,
                "t0_sun_elevation_deg": wx.sun_elevation_deg if wx is not None else np.nan,
                "t0_night": wx.night if wx is not None else np.nan,
                **{k: getattr(r, k) for k in b.columns if k not in ("event_id", "sublink_id", "cml_id",
                                                                    "split_event")}})
            if i_end == i0:
                continue                                            # zero-length tail: index row only

            # --- dry-weather rise later in the tail (descriptive; README "Report definitions")
            sm = pd.Series(A[seg]).rolling(RC["rise_smooth_min"], min_periods=RC["rise_smooth_min"] * 2 // 3).median()
            sm = sm.to_numpy()
            ok = (m >= fd5 * 5 - i0) & (st_seg == DRY) & np.isfinite(sm)
            if ok.any():
                runmin = np.minimum.accumulate(np.where(ok, sm, np.inf))
                hit = np.nonzero(ok & (sm - runmin >= RC["rise_db"]))[0]
                if len(hit):
                    k = hit[0]
                    rise_rows.append({"tail_id": r.tail_id, "sublink_id": s, "event_id": r.event_id,
                                      "rise_start_min": int(k), "rise_start_utc": T1[i0 + k],
                                      "hour_utc": T1[i0 + k].hour, "hour_local": (T1[i0 + k] + LOCAL).hour,
                                      "dew_risk": bool(dew1[i0 + k]),
                                      "level_before_db": float(runmin[k]), "rise_db": float(sm[k] - runmin[k]),
                                      "rise_max_db": float(np.max((sm - runmin)[k:][ok[k:]]))})   # largest later

            parts.append(pd.DataFrame({
                "tail_id": np.full(len(m), r.tail_id, "int32"), "sublink_id": s, "cml_id": np.int32(c),
                "event_id": r.event_id, "t": T1[seg], "minutes_since_t0": m.astype("int16"),
                "A_db": A[seg].astype("float32"), "n_valid_tl": nv[seg].astype("uint8"),
                "link_state": st_seg.astype("int8"),
                "R_max_mmh": np.repeat(rmax, 5)[seg].astype("float32"),
                "R_mean_mmh": np.repeat(rmean, 5)[seg].astype("float32"),
                "gauge_rate_mmh": np.repeat(grate5, 5)[seg].astype("float32"),
                "flag_spike": (fl[seg] & QC["spike"]) > 0, "flag_tsl_step": (fl[seg] & QC["tsl_step"]) > 0,
                "flag_network_gap": (fl[seg] & QC["network_gap"]) > 0, "flag_missing": (fl[seg] & QC["missing"]) > 0,
                "dew_risk": dew1[seg], "in_dry_jump": jump1[seg]}))
    if parts:
        writer.write_table(pa.Table.from_pandas(pd.concat(parts, ignore_index=True), schema=SCHEMA,
                                                preserve_index=False))
    print(f"   sub-links {a + 1}-{a + 32} done ({time.time() - clock:.0f} s)")
writer.close()

idx = pd.DataFrame(index_rows).sort_values("tail_id").reset_index(drop=True)
save(idx.round(4), P.TAILS_INDEX,
            index=False)                                          # 4 decimals: keeps the file < 10 MB
rises = pd.DataFrame(rise_rows)
save(rises, P.STEP1_RISES, index=False)

# =================================================================================================== 2. A(t) check
rows = []
for tid in sorted(check_A):
    r = pairs.loc[pairs.tail_id == tid].iloc[0]
    ref = bl.a_of_t(r.cml_id, r.event_id, res="1min", pad="12h", d=d)[r.sublink_id]
    mine = pd.Series(check_A[tid], index=T1).reindex(ref.index)
    rows.append({"tail_id": tid, "event_id": r.event_id, "sublink_id": r.sublink_id, "n_minutes": len(ref),
                 "n_nan": int(ref.isna().sum()),
                 "max_abs_diff_db": float(np.nanmax(np.abs(ref.to_numpy() - mine.to_numpy()))),
                 "exactly_equal": bool(np.array_equal(ref.to_numpy(), mine.to_numpy(), equal_nan=True))})
acheck = pd.DataFrame(rows)
save(acheck, P.STEP1_A_CHECK, index=False)
print("3. A(t) check vs baseline.a_of_t:\n" + acheck.to_string(index=False))

# =================================================================================================== 5. report numbers
print("4. describing the tails ...")
long = pq.read_table(P.TAILS_1MIN, columns=["tail_id", "minutes_since_t0", "A_db", "link_state"]).to_pandas()
edges = RC["a_start_bins_db"]
labels = [f"< {edges[0]}"] + [f"{lo}-{hi}" for lo, hi in zip(edges[:-1], edges[1:])] + [f"> {edges[-1]}"]
idx["A_start_bin"] = pd.cut(idx.A_start_db, [-np.inf, *edges, np.inf], labels=labels, right=False)
idx["thr_db"] = np.maximum(RC["return_floor_db"], RC["return_std_factor"] * idx.offset_std.fillna(0))
long = long.merge(idx[["tail_id", "A_start_bin", "thr_db"]], on="tail_id")

# return to the dry level: A at minute k below the tail's threshold (tails still running at k with a value at k)
ret = []
for k in RC["return_minutes"]:
    at = long[long.minutes_since_t0 == k].dropna(subset=["A_db"])
    for name, g in [("all", at)] + [(str(lab), at[at.A_start_bin == lab]) for lab in labels]:
        ret.append({"A_start_bin": name, "minute": k, "n_tails": len(g),
                    "share_below": round(float((g.A_db < g.thr_db).mean()), 3) if len(g) else np.nan})
ret = pd.DataFrame(ret)
save(ret, P.STEP1_RETURN, index=False)
ret_wide = ret.pivot(index="A_start_bin", columns="minute", values="share_below").reindex(["all"] + labels)
ret_n = ret.pivot(index="A_start_bin", columns="minute", values="n_tails").reindex(["all"] + labels)

# median curves per A_start bin
curves = (long.dropna(subset=["A_db"]).groupby(["A_start_bin", "minutes_since_t0"], observed=True).A_db
          .agg(["median", "count"]).reset_index())
save(curves, P.STEP1_MEDIAN_CURVES, index=False)

# =================================================================================================== figures
# distributions of A_start and of the max A in the wet phase (one series per panel)
fig, ax = plt.subplots(1, 2, figsize=(12, 4))
for a_, col, title in ((ax[0], "A_start_db", "A_start = mean A over [t0, t0 + 5 min)"),
                       (ax[1], "wet_A_max_db", "max A during the wet phase")):
    v = idx[col].dropna().clip(-1, 15)
    a_.hist(v, bins=np.arange(-1, 15.25, 0.25), color=COLORS[0], edgecolor="white", linewidth=0.5)
    a_.set(title=title, xlabel="dB (clipped to [-1, 15])", ylabel="tails")
    a_.set_yscale("log")
for a_ in ax:
    a_.grid(alpha=0.25)
    for sp in ("top", "right"):
        a_.spines[sp].set_visible(False)
fig.tight_layout()
save(fig, P.STEP1_FIG_DISTRIBUTIONS, dpi=110)
plt.close(fig)

# median curves, linear and log y
fig, ax = plt.subplots(1, 2, figsize=(13, 4.5))
for i, lab in enumerate(labels):
    cv = curves[(curves.A_start_bin == lab) & (curves["count"] >= 20)]          # drawn only where >= 20 tails
    n0 = int(idx.A_start_bin.eq(lab).sum())
    for a_ in ax:
        y = cv["median"].to_numpy()
        a_.plot(cv.minutes_since_t0 / 60, np.where(y > 0, y, np.nan) if a_ is ax[1] else y, color=COLORS[i],
                lw=2, label=f"A_start {lab} dB (n = {n0})")
ax[0].axhline(0, color=INK2, lw=0.8)
ax[0].set(title="median A since t0 (linear)", ylabel="A (dB)")
ax[1].set(title="median A since t0 (log y from 0.01 dB; values <= 0 not drawn)", ylabel="A (dB), log")
ax[1].set_yscale("log")
ax[1].set_ylim(0.01, None)                                          # below 0.01 dB is float32 rounding of a ~0 median
for a_ in ax:
    a_.set_xlabel("hours since t0")
    a_.grid(alpha=0.25)
    for sp in ("top", "right"):
        a_.spines[sp].set_visible(False)
ax[0].legend(frameon=False, fontsize=8)
fig.suptitle("Median tail curve per A_start bin (minutes with >= 20 tails)", fontsize=10)
fig.tight_layout()
save(fig, P.STEP1_FIG_MEDIAN_CURVES, dpi=110)
plt.close(fig)

# resolution check: 10 s, 1 min, 5 min A for the first hours of random tails with a clear start
cand = idx[(idx.A_start_db >= RC["resolution_min_a_start_db"]) & (idx.length_min > 0)]
pick = cand.sample(min(RC["resolution_panels"], len(cand)), random_state=RC["seed"]).sort_values("tail_id")
H = pd.Timedelta(hours=RC["resolution_hours"])
fig = plt.figure(figsize=(16, 16))
outer = fig.add_gridspec(4, 3, hspace=0.35, wspace=0.25)          # one cell per tail: A on top, R_max below
one_min = pq.read_table(P.TAILS_1MIN, columns=["tail_id", "t", "A_db"],
                        filters=[("tail_id", "in", pick.tail_id.tolist())]).to_pandas()
for n, r in enumerate(pick.itertuples()):
    colp = n % 3
    inner = outer[n // 3, colp].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
    a_top = fig.add_subplot(inner[0])
    a_bot = fig.add_subplot(inner[1], sharex=a_top)
    t0, t1 = r.t0, r.t0 + H
    a10 = d.cml([r.sublink_id], t0, t1, res="10s", variables=["tl"]).tl.to_pandas().iloc[:, 0] - r.baseline_db
    a5 = d.cml([r.sublink_id], t0, t1, res="5min", variables=["tl_mean"]).tl_mean.to_pandas().iloc[:, 0] - r.baseline_db
    a1 = one_min[one_min.tail_id == r.tail_id].set_index("t").A_db.loc[t0:t1 - pd.Timedelta("1ns")]
    stt = d.link_state([r.cml_id], t0, t1).state.to_pandas().iloc[:, 0]
    rm = d.radar_along([r.cml_id], t0, t1).R_max.to_pandas().iloc[:, 0]
    h = lambda ix: (pd.DatetimeIndex(ix) - t0) / pd.Timedelta(hours=1)   # noqa: E731  hours since t0
    for a_ in (a_top, a_bot):
        for ti, sv in stt.items():
            if sv != DRY:
                a_.axvspan(h([ti])[0], h([ti])[0] + 5 / 60, color="#4c8ed9" if sv == WET else "#999999",
                           alpha=0.25 if sv == WET else 0.15, lw=0)
        a_.axvline((r.t_end - t0) / pd.Timedelta(hours=1), color=INK2, ls=":", lw=1) if r.t_end < t1 else None
    a_top.plot(h(a10.index), a10.to_numpy(), color=COLORS[0], lw=0.6, alpha=0.6, label="10 s")
    a_top.step(h(a1.index), a1.to_numpy(), where="post", color=COLORS[1], lw=1.2, label="1 min")
    a_top.step(h(a5.index), a5.to_numpy(), where="post", color=COLORS[2], lw=1.8, label="5 min")
    a_top.axhline(0, color=INK2, lw=0.6)
    a_top.set_title(f"{r.sublink_id} {r.event_id} ({r.length_m / 1000:.1f} km, {r.frequency_ghz:.0f} GHz), "
                    f"ended_by {r.ended_by}", fontsize=8)
    a_top.tick_params(labelsize=7, labelbottom=False)
    a_bot.step(h(rm.index), rm.to_numpy(), where="post", color=INK2, lw=1)
    a_bot.tick_params(labelsize=7)
    a_top.set_xlim(0, RC["resolution_hours"])
    a_bot.set_xlim(0, RC["resolution_hours"])
    if colp == 0:
        a_top.set_ylabel("A (dB)", fontsize=8)
        a_bot.set_ylabel("R_max\n(mm/h)", fontsize=7)
    if n == 0:
        a_top.legend(fontsize=7, frameon=False, ncol=3)
    for a_ in (a_top, a_bot):
        for sp in ("top", "right"):
            a_.spines[sp].set_visible(False)
fig.suptitle(f"First {RC['resolution_hours']:g} h of {len(pick)} random tails with A_start >= "
             f"{RC['resolution_min_a_start_db']} dB (seed {RC['seed']}). Top: A at 10 s / 1 min / 5 min; bottom: radar "
             "R_max on the path.\nShading: link_state wet (blue) / unknown (grey); dotted: t_end. x: hours since t0.",
             fontsize=10, y=0.98)
fig.subplots_adjust(top=0.95)
save(fig, P.STEP1_FIG_RESOLUTION, dpi=100, bbox_inches="tight")
plt.close(fig)

# =================================================================================================== 5. the report
def md(df, index=True):
    """DataFrame -> markdown table."""
    df = df.reset_index() if index else df
    head = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns)
    return head + "\n" + "\n".join("| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
                                   for row in df.itertuples(index=False))


def q(s):
    """Percentile summary of a column."""
    s = s.dropna()
    return {"n": len(s), **{k: round(float(s.quantile(v)), 2) for k, v in
                            (("p10", .1), ("p25", .25), ("median", .5), ("p75", .75), ("p90", .9), ("p99", .99))},
            "max": round(float(s.max()), 2)}


rel = lambda p: Path(p).resolve().relative_to(P.REPO).as_posix()    # noqa: E731
nonzero = idx[idx.length_min > 0]
per_sub = idx.groupby("sublink_id").size()
per_sub_hist = pd.cut(per_sub, np.arange(per_sub.min() // 5 * 5, per_sub.max() + 5, 5), right=False)
per_sub_hist = per_sub_hist.value_counts().sort_index()
split_tab = idx.groupby("split_event").agg(tails=("tail_id", "size"), events=("event_id", "nunique"),
                                           sublinks=("sublink_id", "nunique"))
ended_tab = idx.groupby("ended_by").agg(tails=("tail_id", "size"), median_length_min=("length_min", "median"))
band_tab = pd.crosstab(idx.length_class, idx.frequency_band, margins=True)
lens = q(idx.length_min)
dist = pd.DataFrame({"A_start_db": q(idx.A_start_db), "wet_A_max_db": q(idx.wet_A_max_db),
                     "length_min": lens}).T
n_rise = len(rises)
rise_hours = rises.hour_local.value_counts().sort_index() if n_rise else pd.Series(dtype=int)
size_mb = P.TAILS_1MIN.stat().st_size / 1e6
dup_t0 = int(idx.duplicated(["sublink_id", "t0"]).sum())
overlap = 0
for s, g in idx[idx.length_min > 0].sort_values("t0").groupby("sublink_id"):
    overlap += int((g.t0.to_numpy()[1:] < g.t_end.to_numpy()[:-1]).sum())   # a tail starts before the previous ends

report = f"""# WAA_estm step 1: all drying tails (descriptive, no selection)

Written by `research/WAA_estm/step1_extract_tails.py` (settings `waa` in `research/config.yaml`; tail definition in
the module README). Nothing here selects, fits or filters tails: every pair with a baseline and a wet bin is a tail.
In a tail A_rain = 0 by construction, so A(t) = δ(t) = WAA + baseline residual.

Run time {time.time() - clock:.0f} s. Products: `{rel(P.TAILS_1MIN)}` ({size_mb:.0f} MB, {len(long)} tail minutes),
`{rel(P.TAILS_INDEX)}` ({len(idx)} tails).

## A(t) check

A(t) is built here as 1 min `tl_mean` (float32 -> float64) - `baseline_db`. On {len(acheck)} random pairs (seed
{RC['seed']}) it is compared with `baseline.a_of_t(cml_id, event_id, res="1min", pad="12h")` over that function's whole
window (`{rel(P.STEP1_A_CHECK)}`):

{md(acheck, index=False)}

## Counts

Tails: **{len(idx)}** ({len(nonzero)} with at least one minute; {len(idx) - len(nonzero)} of length 0), on
{idx.sublink_id.nunique()} sub-links of {idx.cml_id.nunique()} links, {idx.event_id.nunique()} events.

By split (train + val events on train links only; test has no baseline):

{md(split_tab)}

By `ended_by`:

{md(ended_tab)}

Tails per sub-link (histogram: number of sub-links with that many tails; median {per_sub.median():.0f}, max
{per_sub.max()}):

{md(per_sub_hist.rename("sub-links").rename_axis("tails per sub-link").to_frame())}

By length class x frequency band (sub-link 1 band, as in the release):

{md(band_tab)}

## Distributions

Percentiles (dB; length in minutes). `A_start` = mean A over [t0, t0 + 5 min); `wet_A_max` = largest 1 min A between
the link's first wet bin in the event and t0. Figure: `{rel(P.STEP1_FIG_DISTRIBUTIONS)}`.

{md(dist)}

![distributions](../figures/{P.STEP1_FIG_DISTRIBUTIONS.name})

A_start bins: {idx.A_start_bin.value_counts().reindex(labels).to_dict()}; A_start NaN (no valid minute in the first
5): {int(idx.A_start_db.isna().sum())}.

## How fast A returns to the dry level

Share of tails with A(t0 + k) < max({RC['return_floor_db']} dB, {RC['return_std_factor']:g} x `offset_std`), among
the tails still running at minute k and with a valid A there (`{rel(P.STEP1_RETURN)}`). Rows: A_start bin.

{md(ret_wide)}

Number of tails behind each share:

{md(ret_n.astype("Int64"))}

## Median tail curve per A_start bin

`{rel(P.STEP1_MEDIAN_CURVES)}` (median and count per minute). Drawn where >= 20 tails remain; the curve mixes
fewer tails at later minutes (tails end by next_wet / gap / cap), so its far end is a different population.

![median curves](../figures/{P.STEP1_FIG_MEDIAN_CURVES.name})

## Dry-weather rises later in the tail

Definition (descriptive, settings `waa.report`): A smoothed by a trailing {RC['rise_smooth_min']} min rolling median;
from `t_first_dry` on, only minutes with link_state dry; a rise starts at the first such minute where the smoothed A is
>= {RC['rise_db']} dB above its running minimum since `t_first_dry`. Table: `{rel(P.STEP1_RISES)}`.

- Tails with a rise: **{n_rise}** of {len(nonzero)} ({n_rise / max(len(nonzero), 1):.0%}).
- Start, minutes since t0: {q(rises.rise_start_min) if n_rise else '-'}.
- Rise size at its start (dB): {q(rises.rise_db) if n_rise else '-'}.
- Largest rise after its start (dB): {q(rises.rise_max_db) if n_rise else '-'}; tails whose rise reaches >= 0.6 dB
  (two RSL steps): {int((rises.rise_max_db >= 0.6).sum()) if n_rise else 0}; >= 1 dB: {int((rises.rise_max_db >= 1).sum()) if n_rise else 0}.
- In a dew-risk hour at its start: {int(rises.dew_risk.sum()) if n_rise else 0} of {n_rise}.
- Start hour, local time (UTC + 2 h): {rise_hours.to_dict()}.

The threshold equals one RSL quantization step (0.3 dB) and a 15 min median of 1 min means can still flip between two
levels, so part of these "rises" may be quantization; the share is an upper bound until checked.

## Resolution check

{len(pick)} random tails with A_start >= {RC['resolution_min_a_start_db']} dB (seed {RC['seed']}): A at 10 s, 1 min
and 5 min over the first {RC['resolution_hours']:g} h, link_state shading, radar R_max below.

![resolution](../figures/{P.STEP1_FIG_RESOLUTION.name})

## Unclear or inconsistent

1. **Tails of length 0**: {len(idx) - len(nonzero)}. They stay in the index (no selection) with no rows in the 1 min
   product. Causes: t0 clipped to the season's last bin while still wet ({int(((idx.length_min == 0) & (idx.ended_by == 'season_end')).sum())}, all from the last event);
   **t0 itself is a WET bin** ({int(((idx.length_min == 0) & (idx.ended_by == 'next_wet')).sum())}: `last_wet_end` is searched only inside the network event, and the link
   was still wet at the network event end, so the "tail" stops at once); A NaN for >= {GAP_MIN} min from t0
   ({int(((idx.length_min == 0) & (idx.ended_by == 'gap')).sum())}). Tails with a wet t0 in total: {int((idx.state_at_t0 == WET).sum())}; state at t0:
   {idx.state_at_t0.map({0: 'dry', 1: 'wet', -1: 'unknown'}).value_counts().to_dict()}.
2. **`ended_by = season_end`**: {int((idx.ended_by == 'season_end').sum())} tails (the data stop first, or `last_wet_end`
   was clipped by `make_baseline_table.py` to the last bin of the season); not in the requested set
   {{next_wet, gap, cap}}.
3. **`t_first_dry` after `t_end`**: {int(idx.first_dry_after_end.sum())} tails never reach a dry bin before they end
   (all margin / unknown); {int(idx.t_first_dry.isna().sum())} have no dry bin left in the season.
4. **Overlapping tails on one sub-link**: {overlap}; same (sub-link, t0) twice: {dup_t0}. A tail can run past the
   network event bounds into the next event's time if the link itself stays non-wet.
5. **A_start NaN**: {int(idx.A_start_db.isna().sum())} tails have no valid A in their first 5 minutes.
6. **The margin is part of the tail.** t0 is normally an `unknown` bin (rain within 30 min, or 0 < R_max < 0.1 mm/h).
   Median `margin_min` {idx.margin_min.median():.0f} min; tails with R_max > 0 in the margin:
   {int((idx.margin_R_max_max_mmh > 0).sum())}; with a gauge rate > 0 there: {int((idx.margin_gauge_rate_max_mmh > 0).sum())}.
   So the first minutes of a tail can still hold light path rain.
7. **Gauge columns use the nearest gauge only** (if within {TC['gauge_max_dist_m']} m: {int((idx.gauge_id != '').sum())}
   tails), while link_state uses every gauge within 1 km. SMHI GbgA is 15 min (broadcast to 5 min, then 1 min).
8. **Broadcast resolutions.** link_state, radar and gauge are 5 min, dew risk 1 h, on a 1 min A: within a 5 min bin
   the reference does not change, so the exact minute rain stops is unknown to +-5 min.
9. **`wet_A_max` window** starts at the link's first wet bin in the network event, so for multi-run events it covers
   several showers and the pauses between them.
10. **The 1 min product is large** ({size_mb:.0f} MB): git-ignored, rebuilt by this script.
11. **Rises at 0.3 dB** may include quantization flips (see above).
"""
save(report, P.STEP1_REPORT)
print("\n" + "=" * 100 + "\n" + report)
print(f"done in {time.time() - clock:.0f} s")
