"""WAA_estm step 2c: last cleaning fixes, then FREEZE the tail set for modelling (tail set "v1", config waa.freeze).
TRAIN only for every number and figure; VAL gets the same rules and files but is only counted. No fitting.

Fixes on top of step 2b (README.md, "Tail set v1"):
  1. a radar NaN bin on the path counts as "unknown" (like a rain trace), not as "no trace": t0_trace and the strict end
     are recomputed with this rule;
  2. rise cut v2: from t0_trace, on ALL minutes, the first run of >= 10 min where the 10 min trailing median of A is
     >= 0.6 dB above its running minimum since t0_trace; cut at the last minute before that run with an excess
     < 0.15 dB (kept). Applied to the strict end.
  3. event diagnostics per TRAIN event; events where > 30 % of the tails are cut by the rise rule -> event_rise_heavy;
     plots of the two events named in config (waa.freeze.look_events).
  4. freeze: tails_model.parquet / tails_model_index.csv (main set) and tails_model_alt.* (A_start_trace >= 0.6 dB).

Reads tails_selected.csv (run steps 2 and 2b first), tails_1min.parquet, the release (radar, gauges) and the
baseline weather product. Run from the repository root (about 1 minute):
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step2c_freeze.py
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

sys.stdout.reconfigure(encoding="utf-8")
MODULE = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE))
import WAA_estm_paths as P

CFG = yaml.safe_load(open(P.CONFIG_FILE, encoding="utf-8"))
S2, B, F = CFG["waa"]["selection"], CFG["waa"]["selection_2b"], CFG["waa"]["freeze"]
TAIL_SET = F["tail_set"]
WET, DRY = 1, 0
INK, INK2 = "#1a1a1a", "#52514e"
LOCAL = pd.Timedelta(hours=2)
clock = time.time()
# median-curve groups: line style first, colour second (the reader may not separate colours)
CURVES = [(f"strict {TAIL_SET}", "-", "#2a78d6"), ("step 2b strict", "--", "#eb6834"), ("control", "-.", "#7a7a7a")]


def tidy(ax):
    ax.grid(alpha=0.25)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def first_run(cond, n):
    """Start of the first run of >= n consecutive True values (None if there is none)."""
    c = np.concatenate([[0], cond.astype(int), [0]])
    st, en = np.nonzero(np.diff(c) == 1)[0], np.nonzero(np.diff(c) == -1)[0]
    ok = np.nonzero(en - st >= n)[0]
    return int(st[ok[0]]) if len(ok) else None


def window_mean(a, k, w, end):
    if k + w > end:
        return np.nan
    v = a[k:k + w]
    return float(np.nanmean(v)) if np.isfinite(v).any() else np.nan


def minutes(t, t0):
    return int(round((t - t0) / pd.Timedelta("1min")))


# =================================================================================================== load
dates = ["t0", "t0_trace", "t_first_dry", "t_end", "event_start", "event_end", "t_use_end_lenient", "t_use_end_strict"]
sel = pd.read_csv(P.TAILS_SELECTED, parse_dates=dates)
assert "t0_trace" in sel.columns, "run step2b_fix_tails.py first (tails_selected.csv has no step 2b columns)"
long = pq.read_table(P.TAILS_1MIN, columns=["tail_id", "minutes_since_t0", "A_db", "link_state"]).to_pandas()
long = long.sort_values(["tail_id", "minutes_since_t0"], kind="stable").reset_index(drop=True)

d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
ra = d.store("radar_along_links_5min")
T5 = pd.DatetimeIndex(ra.time.values)
rmax5 = ra.R_max.load().to_pandas()
rmean5 = ra.R_mean.load().to_pandas()
city = d.store("gauges_city_5min")
smhi = d.store("gauges_smhi_15min")
g_amt = pd.concat([city.amount_mm.load().to_pandas(),
                   smhi.amount_mm.load().to_pandas().reindex(T5, method="ffill", limit=2)], axis=1)
g_rate = pd.concat([city.rate_mm_h.load().to_pandas(),
                    smhi.rate_mm_h.load().to_pandas().reindex(T5, method="ffill", limit=2)], axis=1)
lgm = d.table("link_gauge_matches")
near = lgm[lgm.distance_to_path_m <= B["gauge_max_dist_m"]].groupby("cml_id").gauge_id.apply(list).to_dict()
trace5, nan5 = {}, {}
for c in sel.cml_id.unique():
    t = (rmax5[c].fillna(0) > 0).to_numpy().copy()                  # radar trace (as in step 2b)
    for g in near.get(c, []):
        t |= (g_amt[g].fillna(0) > 0).to_numpy()                    # any gauge within 1 km
    trace5[c] = t
    nan5[c] = rmax5[c].isna().to_numpy()                            # radar NaN on the path: unknown
wx = pd.read_csv(P.DEW_RISK_HOURLY, parse_dates=["time_utc"]).set_index("time_utc")
print(f"loaded ({time.time() - clock:.0f} s)")

# =================================================================================================== 1 + 2. per tail
bounds = np.searchsorted(long.tail_id.to_numpy(), sel.tail_id.to_numpy(), side="left")
ends = np.searchsorted(long.tail_id.to_numpy(), sel.tail_id.to_numpy(), side="right")
A_all = long.A_db.to_numpy("float64")
rows = []
for i, r in enumerate(sel.itertuples()):
    lo, hi = bounds[i], ends[i]
    L = hi - lo
    if L == 0:
        rows.append({"tail_id": r.tail_id, "k_tr": 0, "end_v1": 0, "end_pre": 0, "usable_v1_min": 0,
                     "truncated_by_v1": r.ended_by})
        continue
    a = A_all[lo:hi]
    i05 = T5.searchsorted(r.t0)
    fd = int(r.margin_minutes) if pd.notna(r.margin_minutes) else L
    tr, nn = trace5[r.cml_id], nan5[r.cml_id]
    unk = tr | nn if F["radar_nan_is_unknown"] else tr

    # --- 1. t0_trace and strict trace end with radar NaN = unknown
    mb = np.nonzero(unk[i05:i05 + int(np.ceil(min(fd, L) / 5))])[0]
    k_tr = min(int((mb[-1] + 1) * 5) if len(mb) else 0, L)
    nxt = np.nonzero(unk[i05 + k_tr // 5:])[0]
    k_s = k_tr + int(nxt[0]) * 5 if len(nxt) else 10 ** 9
    s_label = ("rain_trace" if tr[i05 + k_s // 5] else "radar_nan") if k_s < 10 ** 9 else ""
    k_len = min(minutes(r.t_use_end_lenient, r.t0), L) if pd.notna(r.t_use_end_lenient) else L

    # --- 2. rise cut v2: from t0_trace, ALL minutes
    sm = pd.Series(a).rolling(F["rise_smooth_min"], min_periods=F["rise_smooth_min"] * 2 // 3).median().to_numpy()
    m = np.arange(L)
    mask = (m >= k_tr) & np.isfinite(sm)
    runmin = np.minimum.accumulate(np.where(mask, sm, np.inf))
    excess = sm - runmin
    k_run = first_run(mask & (excess >= F["rise_db"]), F["rise_min_duration_min"])
    k_rise = 10 ** 9
    if k_run is not None:
        quiet = np.nonzero(mask[:k_run] & (excess[:k_run] < F["rise_quiet_db"]))[0]
        k_rise = int(quiet[-1]) + 1 if len(quiet) else k_tr

    cands = {r.truncated_by_lenient: k_len, s_label or "rain_trace": k_s, "rise_v2": k_rise}
    by = min(cands, key=lambda k: cands[k])                          # ties: lenient reason, then trace, then rise
    end_v1 = max(min(cands[by], L), k_tr)
    end_pre = max(min(k_len, k_s, L), k_tr)                          # the strict end without the rise cut v2
    rows.append({"tail_id": r.tail_id, "k_tr": k_tr, "end_v1": end_v1, "end_pre": end_pre,
                 "usable_v1_min": end_v1 - k_tr, "usable_pre_min": end_pre - k_tr, "truncated_by_v1": by,
                 "rise_v2_found": k_run is not None, "A_start_v1_db": window_mean(a, k_tr, S2["a_dry_window_min"], L),
                 "k_tr_2b": minutes(r.t0_trace, r.t0),
                 "k_s_2b": minutes(r.t_use_end_strict, r.t0) if pd.notna(r.t_use_end_strict) else 0})
v = sel.merge(pd.DataFrame(rows), on="tail_id", how="left").copy()   # one block: no fragmentation warnings
v["t0_trace_v1"] = v.t0 + pd.to_timedelta(v.k_tr, unit="min")
v["t_end_v1"] = v.t0 + pd.to_timedelta(v.end_v1, unit="min")
print(f"per-tail pass done ({time.time() - clock:.0f} s)")

# --- sets: step 2 exclusions + noisy baseline (offset_std > 0.3) + A_start NaN, then thresholds
ok = (v.exclusion_reason.fillna("") == "") & ~(v.offset_std_unc > B["max_offset_std_db"]) & v.A_start_v1_db.notna()
v["pre_main"] = ok & (v.A_start_v1_db >= S2["main_a_start_db"]) & (v.usable_pre_min >= S2["min_usable_min"])
v["main_v1"] = ok & (v.A_start_v1_db >= S2["main_a_start_db"]) & (v.usable_v1_min >= S2["min_usable_min"])
v["alt_v1"] = ok & (v.A_start_v1_db >= S2["alt_a_start_db"]) & (v.usable_v1_min >= S2["min_usable_min"])
v["control_v1"] = ok & (v.A_start_v1_db < S2["control_max_a_start_db"]) & (v.usable_v1_min > 0)
v["cut_by_rise_v2"] = v.truncated_by_v1 == "rise_v2"

# =================================================================================================== 3. events
rv = rmax5.notna()


def radar_valid(c, t0, t1):
    s = rv[c].loc[t0:t1 - pd.Timedelta("1ns")]
    return float(s.mean()) if len(s) else np.nan


ev_rows = []
for (eid, split), g in v[v.pre_main].groupby(["event_id", "split_event"]):  # train AND val (val: flag only)
    e0, e1 = g.event_start.iloc[0], g.event_end.iloc[0]
    win1 = e1 + pd.Timedelta(hours=F["event_window_after_h"])
    gm = g[g.main_v1]
    ev_rows.append({"event_id": eid, "split_event": split, "event_start": e0, "event_end": e1,
                    "start_hour_local": (e0 + LOCAL).hour, "end_hour_local": (e1 + LOCAL).hour,
                    "median_t0_trace_hour_local": float((g.t0_trace_v1 + LOCAL).dt.hour.median()),
                    "tails_before_rise_cut": len(g), "tails_v1": len(gm), "links_v1": gm.cml_id.nunique(),
                    "share_cut_by_rise_v2": round(float(g.cut_by_rise_v2.mean()), 3),
                    "share_level_flags_v1": round(float((gm.level_shift | gm.undershoot).mean()), 3) if len(gm) else np.nan,
                    "radar_valid_fraction": round(float(np.nanmean([radar_valid(c, e0, win1) for c in g.cml_id.unique()])), 3),
                    "mean_T_t0_c": round(g.t0_temperature_c.mean(), 1), "mean_RH_t0_pct": round(g.t0_humidity_pct.mean(), 1),
                    "mean_wind_t0_ms": round(g.t0_wind_speed_ms.mean(), 1)})
events = pd.DataFrame(ev_rows)
events["event_rise_heavy"] = events.share_cut_by_rise_v2 > F["event_rise_heavy_share"]
v["event_rise_heavy"] = v.event_id.map(events.set_index("event_id").event_rise_heavy).fillna(False).astype(bool)
ev_tr = events[events.split_event == "train"].drop(columns="split_event").sort_values("event_start")
ev_tr.to_csv(P.STEP2C_EVENTS, index=False)

# =================================================================================================== 4. freeze
tr_ = v[v.split_event == "train"]
meta_cols = ["frequency_ghz", "polarization", "L_km", "length_m", "band", "length_class", "site_0_antenna_diameter_m",
             "site_1_antenna_diameter_m", "antennas_m"]
wx_cols = ["t0_temperature_c", "t0_humidity_pct", "t0_wind_speed_ms", "t0_dew_point_depression_c",
           "t0_sun_elevation_deg", "t0_night"]


def freeze(flag, idx_path, rows_path, label):
    """Write one frozen set: index (one row per tail) and 1 min rows from t0_trace to the v1 end."""
    s = v[v[flag]].copy()
    s.insert(0, "tail_set", f"{TAIL_SET}{'' if label == 'main' else '_alt'}")
    cols = (["tail_set", "tail_id", "sublink_id", "cml_id", "event_id", "split_event", "t0", "t0_trace_v1", "t_end_v1",
             "usable_v1_min", "truncated_by_v1", "A_start_v1_db", "level_shift", "undershoot", "event_rise_heavy",
             "main_v1"] + meta_cols + ["offset_std_unc"] + wx_cols)
    out = s[cols].rename(columns={"t0_trace_v1": "t0_trace", "t_end_v1": "t_end", "usable_v1_min": "usable_min",
                                  "truncated_by_v1": "truncated_by", "A_start_v1_db": "A_start_trace_db",
                                  "main_v1": "in_main_set", "offset_std_unc": "offset_std"})
    num = out.select_dtypes("number").columns
    out.assign(**out[num].round(4)).to_csv(idx_path, index=False)
    # 1 min rows: minutes [k_tr, end_v1) of each tail
    lr = long.merge(s[["tail_id", "sublink_id", "cml_id", "event_id", "split_event", "t0", "k_tr", "end_v1"]],
                    on="tail_id")
    lr = lr[(lr.minutes_since_t0 >= lr.k_tr) & (lr.minutes_since_t0 < lr.end_v1)]
    rows_ = pd.DataFrame({"tail_set": out.tail_set.iloc[0] if len(out) else TAIL_SET,
                          "tail_id": lr.tail_id.astype("int32"), "sublink_id": lr.sublink_id,
                          "cml_id": lr.cml_id.astype("int32"), "event_id": lr.event_id, "split_event": lr.split_event,
                          "t": lr.t0 + pd.to_timedelta(lr.minutes_since_t0, unit="min"),
                          "minutes_since_t0_trace": (lr.minutes_since_t0 - lr.k_tr).astype("int16"),
                          "A_db": lr.A_db.astype("float32"), "link_state": lr.link_state.astype("int8")})
    pq.write_table(pa.Table.from_pandas(rows_, preserve_index=False), rows_path, compression="zstd")
    return out, rows_


idx_main, rows_main = freeze("main_v1", P.TAILS_MODEL_INDEX, P.TAILS_MODEL, "main")
idx_alt, rows_alt = freeze("alt_v1", P.TAILS_MODEL_ALT_INDEX, P.TAILS_MODEL_ALT, "alt")

# =================================================================================================== changes (TRAIN)
k = tr_[(tr_.exclusion_reason.fillna("") == "") & (tr_.length_min > 0)]
changes = pd.DataFrame([
    ("kept step 2 TRAIN tails", len(k), np.nan),
    ("t0_trace moved by the radar-NaN rule", int((k.k_tr != k.k_tr_2b).sum()), float((k.k_tr - k.k_tr_2b).sum())),
    ("strict end moved earlier by the radar-NaN rule (before the rise cut)",
     int((k.end_pre < np.minimum(k.k_s_2b, k.length_min)).sum()),
     float((np.minimum(k.k_s_2b, k.length_min) - k.end_pre).clip(lower=0).sum())),
    ("rise v2 detected", int(k.rise_v2_found.sum()), np.nan),
    ("end set by rise v2 (cut_by_rise_v2)", int(k.cut_by_rise_v2.sum()), float((k.end_pre - k.end_v1).sum())),
    ("  of them in the pre-rise main set", int(k[k.pre_main].cut_by_rise_v2.sum()),
     float((k[k.pre_main].end_pre - k[k.pre_main].end_v1).sum())),
    ("pre-rise main set (after rule 1)", int(k.pre_main.sum()), float(k[k.pre_main].usable_pre_min.sum())),
    (f"main set {TAIL_SET} (after rules 1-2)", int(k.main_v1.sum()), float(k[k.main_v1].usable_v1_min.sum())),
    ("step 2b strict main set (reference)", int(k.in_main_strict.sum()), float(k[k.in_main_strict].usable_strict_min.sum())),
], columns=["what", "tails", "minutes"])
changes.to_csv(P.STEP2C_CHANGES, index=False)
by_v1 = k[k.main_v1].truncated_by_v1.value_counts()
lost_to_short = int((k.pre_main & ~k.main_v1).sum())

# =================================================================================================== median curves
use = long.merge(v[["tail_id", "split_event", "main_v1", "control_v1", "in_main_strict", "k_tr", "end_v1",
                    "k_tr_2b", "k_s_2b"]], on="tail_id")
use = use[use.split_event == "train"]
grp = {f"strict {TAIL_SET}": (use.main_v1, use.k_tr, use.end_v1),
       "step 2b strict": (use.in_main_strict, use.k_tr_2b, use.k_s_2b),
       "control": (use.control_v1, use.k_tr, use.end_v1)}
cv_list = []
for name, (m_, k0, k1) in grp.items():
    g = use[m_ & (use.minutes_since_t0 >= k0) & (use.minutes_since_t0 < k1)]
    rel_ = g.minutes_since_t0 - k0.loc[g.index]
    cv_list.append(g.A_db.groupby(rel_.to_numpy()).agg(["median", "count"]).rename_axis("rel").reset_index()
                   .assign(group=name))
curves = pd.concat(cv_list)
curves.to_csv(P.STEP2C_MEDIAN_CURVES, index=False)
fig, ax = plt.subplots(2, 2, figsize=(13, 8))
for row, log in enumerate((False, True)):
    for j, xmax in enumerate((90, 360)):
        a_ = ax[row, j]
        sc = 1 if xmax <= 90 else 60
        for name, ls, colr in CURVES:
            c_ = curves[(curves.group == name) & (curves["count"] >= 20) & (curves.rel <= xmax)]
            y, x = c_["median"].to_numpy(), c_.rel.to_numpy() / sc
            a_.plot(x, np.where(y > 0.01, y, np.nan) if log else y, ls=ls, color=colr, lw=2,
                    label=f"{name} (n = {int(c_['count'].max()) if len(c_) else 0})")
            fin = np.isfinite(y) & ((y > 0.01) if log else True)
            if fin.any():
                a_.text(x[fin][-1], y[fin][-1], f" {name}", fontsize=7, va="center")
        a_.set_yscale("log") if log else a_.axhline(0, color=INK2, lw=0.6)
        a_.set(xlabel="minutes since t0_trace" if sc == 1 else "hours since t0_trace",
               ylabel="median A (dB)" + (", log" if log else ""), xlim=(0, xmax / sc * 1.2))
        tidy(a_)
ax[0, 0].legend(frameon=False, fontsize=8)
fig.suptitle(f"Median tail curves (TRAIN), each aligned at its own t0_trace. Line styles: strict {TAIL_SET} solid, "
             "step 2b strict dashed, control dash-dot. Drawn where >= 20 tails; log y from 0.01 dB.", fontsize=10)
fig.tight_layout()
fig.savefig(P.STEP2C_FIG_MEDIAN, dpi=110)
plt.close(fig)

# =================================================================================================== event look plots
MARKS = [("t0", "t0", (0, (1, 0))), ("t0_trace", "t0_trace_v1", (0, (6, 3))), ("first dry", "t_first_dry", (0, (1, 2))),
         ("end v1", "t_end_v1", (0, (6, 2, 1, 2)))]
event_figs = {}
for n_ev, eid in enumerate(F["look_events"]):
    g = v[(v.event_id == eid) & v.pre_main & (v.split_event == "train")]
    pick = g.sample(min(F["look_tails_per_event"], len(g)), random_state=F["look_seed"]).sort_values("tail_id")
    path = Path(str(P.STEP2C_FIG_EVENT).format(event_id=eid))
    event_figs[eid] = (path, len(g), len(pick))
    if not len(pick):
        continue
    H = F["look_hours"]
    fig = plt.figure(figsize=(16, 4.6 * int(np.ceil(len(pick) / 3)) + 1))
    outer = fig.add_gridspec(int(np.ceil(len(pick) / 3)), 3, hspace=0.5, wspace=0.25)
    for n, r in enumerate(pick.itertuples()):
        inner = outer[n // 3, n % 3].subgridspec(3, 1, height_ratios=[3, 1.3, 1], hspace=0.1)
        top = fig.add_subplot(inner[0])
        mid = fig.add_subplot(inner[1], sharex=top)
        bot = fig.add_subplot(inner[2], sharex=top)
        lt = long[(long.tail_id == r.tail_id) & (long.minutes_since_t0 < H * 60)]
        h = lt.minutes_since_t0.to_numpy() / 60
        st5 = lt.link_state.to_numpy()[::5]
        for a_ in (top, mid, bot):
            for kk, s_ in enumerate(st5):
                if s_ == WET:
                    a_.axvspan(kk * 5 / 60, (kk + 1) * 5 / 60, facecolor="#cfcfcf", edgecolor="#555555", hatch="xxx", lw=0)
                elif s_ != DRY:
                    a_.axvspan(kk * 5 / 60, (kk + 1) * 5 / 60, facecolor="none", edgecolor="#8a8a8a", hatch="///", lw=0)
        top.plot(h, lt.A_db.to_numpy(), color=INK, lw=0.8)
        top.axhline(0, color=INK2, lw=0.5)
        for kk, (lab, col, ls) in enumerate(MARKS):
            x = (getattr(r, col) - r.t0) / pd.Timedelta(hours=1) if pd.notna(getattr(r, col)) else np.nan
            if np.isfinite(x) and x <= H:
                for a_ in (top, mid, bot):
                    a_.axvline(x, color=INK, ls=ls, lw=1.1)
                top.text(x, 0.98 - 0.12 * kk, f" {lab}", transform=top.get_xaxis_transform(), fontsize=6, va="top",
                         bbox=dict(facecolor="white", edgecolor="none", pad=0.5, alpha=0.8))
        i05 = T5.searchsorted(r.t0)
        hb = np.arange(H * 12) * 5 / 60
        mid.step(hb, rmax5[r.cml_id].to_numpy()[i05:i05 + H * 12], where="post", color=INK, lw=1.1, ls="-",
                 label="radar R_max")
        mid.step(hb, rmean5[r.cml_id].to_numpy()[i05:i05 + H * 12], where="post", color=INK, lw=1.1, ls="--",
                 label="radar R_mean")
        gl = near.get(r.cml_id, [])
        if gl:
            mid.step(hb, g_rate[gl].max(axis=1).to_numpy()[i05:i05 + H * 12], where="post", color=INK2, lw=1.3,
                     ls=":", label=f"gauge max ({', '.join(gl)})")
        nanb = nan5[r.cml_id][i05:i05 + H * 12]
        mid.plot(hb[nanb] + 2.5 / 60, np.zeros(nanb.sum()), "x", color=INK, ms=4, label="radar NaN")
        mid.legend(fontsize=5, frameon=False, loc="upper right")
        hrs = pd.date_range(r.t0.floor("1h"), r.t0 + pd.Timedelta(hours=H), freq="1h")
        bot.step((hrs - r.t0) / pd.Timedelta(hours=1), wx.humidity.reindex(hrs).to_numpy(), where="post", color=INK,
                 lw=1.1)
        bot.set_ylim(40, 102)
        flags = [f for f, val in (("level_shift", r.level_shift), ("undershoot", r.undershoot),
                                  ("cut by rise v2", r.cut_by_rise_v2)) if val]
        top.set_title(f"{r.sublink_id} | {r.L_km:.1f} km {r.frequency_ghz:.0f} GHz | A_start {r.A_start_v1_db:.1f} dB\n"
                      f"flags: {', '.join(flags) if flags else 'none'} | end v1: {r.truncated_by_v1}", fontsize=7)
        top.set_ylabel("A (dB)", fontsize=7)
        mid.set_ylabel("mm/h", fontsize=7)
        bot.set_ylabel("RH %", fontsize=7)
        for a_ in (top, mid):
            a_.tick_params(labelbottom=False, labelsize=6)
        bot.tick_params(labelsize=6)
        bot.set_xlim(0, H)
        bot.set_xlabel("hours since t0", fontsize=7)
        for a_ in (top, mid, bot):
            for sp in ("top", "right"):
                a_.spines[sp].set_visible(False)
    fig.text(0.5, 0.995, f"{eid}: {len(pick)} random tails of {len(g)} (pre-rise strict main set, seed {F['look_seed']}), "
             f"first {H} h.\nTop: A; lines differ by style and carry labels (t0 solid, t0_trace dashed, first dry dotted, "
             "end v1 dash-dot); hatching /// unknown, xxx wet. Middle: radar R_max solid, R_mean dashed, gauges "
             "within 1 km (max) dotted, x = radar NaN. Bottom: hourly RH at GbgA.", ha="center", va="top", fontsize=9)
    fig.subplots_adjust(top=0.9)
    fig.savefig(path, dpi=90, bbox_inches="tight")
    plt.close(fig)


# =================================================================================================== report
def md(df, index=True):
    df = df.reset_index() if index else df
    head = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns)
    return head + "\n" + "\n".join("| " + " | ".join("" if pd.isna(x) else str(x) for x in row) + " |"
                                   for row in df.itertuples(index=False))


def set_counts(df):
    return {"tails": len(df), "sub-links": df.sublink_id.nunique(), "links": df.cml_id.nunique(),
            "events": df.event_id.nunique()}


rel = lambda p: Path(p).resolve().relative_to(P.REPO).as_posix()    # noqa: E731
counts = pd.DataFrame({f"{name}, {sp}": set_counts(df[df.split_event == sp])
                       for name, df in (("main", idx_main), (f"alternative (>= {S2['alt_a_start_db']} dB)", idx_alt))
                       for sp in ("train", "val")}).T
per_ev = idx_main[idx_main.split_event == "train"].groupby("event_id").size()
heavy = ev_tr[ev_tr.event_rise_heavy]
look_txt = "\n\n".join(f"**{eid}** ({n_all} pre-rise tails, {n_p} shown):\n\n![{eid}](../figures/{p.name})"
                       if n_p else f"**{eid}**: no tail of this event in the pre-rise main set."
                       for eid, (p, n_all, n_p) in event_figs.items())
look_rows = ev_tr[ev_tr.event_id.isin(F["look_events"])]
ev_show = ev_tr.drop(columns=["event_start", "event_end"])

report = f"""# WAA_estm step 2c: last fixes and the frozen tail set {TAIL_SET} (TRAIN only)

Written by `research/WAA_estm/step2c_freeze.py` (settings `waa.freeze`; the rules of steps 2 and 2b stay). No fitting.
TRAIN events for every number and figure; VAL gets the same rules and files and is only counted. Run time
{time.time() - clock:.0f} s.

## 1-2. Radar NaN as unknown, rise cut v2

Changes on the {len(k)} kept step 2 TRAIN tails (`{rel(P.STEP2C_CHANGES)}`; minutes = tail time moved / removed):

{md(changes, index=False)}

Rise cut v2 (from t0_trace, all minutes, 10 min trailing median >= {F['rise_db']} dB above the running minimum for >= {F['rise_min_duration_min']} min;
cut at the last minute before that run with an excess < {F['rise_quiet_db']} dB): pre-rise main-set tails that drop below
{S2['min_usable_min']} min usable because of it: {lost_to_short}. End of the {TAIL_SET} main-set tails: {by_v1.to_dict()}.

## 3. Event diagnostics (TRAIN, `{rel(P.STEP2C_EVENTS)}`)

Per event: tails in the pre-rise main set and in {TAIL_SET}, links, share of the pre-rise tails whose end is set by the
rise cut v2, share of {TAIL_SET} tails with level_shift or undershoot, radar valid fraction on the paths over
[event start, event end + {F['event_window_after_h']} h], local hours (UTC + 2 h), mean weather at t0 (GbgA).
`event_rise_heavy` = share cut by rise > {F['event_rise_heavy_share']:.0%}: {len(heavy)} events ({', '.join(heavy.event_id) or 'none'}).

{md(ev_show, index=False)}

{look_txt}

## 4. Frozen set {TAIL_SET}

Files: `{rel(P.TAILS_MODEL)}` ({len(rows_main)} 1 min rows), `{rel(P.TAILS_MODEL_INDEX)}`; alternative:
`{rel(P.TAILS_MODEL_ALT)}` ({len(rows_alt)} rows), `{rel(P.TAILS_MODEL_ALT_INDEX)}`. Rules: README.md, "Tail set v1".

{md(counts)}

TRAIN main-set tails per event: median {per_ev.median():.0f}, max {per_ev.max()} ({per_ev.idxmax()}); events with 1-2
tails: {int((per_ev <= 2).sum())} of {len(per_ev)}. The largest events: {per_ev.sort_values(ascending=False).head(5).to_dict()}.
Flags in the TRAIN main set: level_shift {int(idx_main[idx_main.split_event == 'train'].level_shift.sum())}, undershoot
{int(idx_main[idx_main.split_event == 'train'].undershoot.sum())}, event_rise_heavy {int(idx_main[idx_main.split_event == 'train'].event_rise_heavy.sum())}.

## 5. Median curve

`{rel(P.STEP2C_MEDIAN_CURVES)}`.

![median curves](../figures/{P.STEP2C_FIG_MEDIAN.name})

## Unclear

1. **level_shift and undershoot are the step 2b values**: level_shift uses t0 and the lenient end, undershoot uses the
   step 2b t0_trace. They were not recomputed for the v1 start / end.
2. **The rise cut v2 runs on unknown minutes too**, so a rise caused by light rain the radar missed is cut (as
   intended), but also a real late WAA rise (dew?) is cut; the two cannot be told apart here.
3. **The rise cut can end a tail at t0_trace** (no quiet minute before the run): such tails get usable 0 and leave the
   set; this is a selection on the early shape of A, though by a fixed rule.
4. **Radar NaN** now blocks the start and ends the strict tail; links / hours with many NaN radar bins lose tails.
   Gauge NaN is still "no trace".
5. **Event shares are small-number statistics**: events with 1-3 pre-rise tails can reach 0 % or 100 % by chance; the
   {F['event_rise_heavy_share']:.0%} rule then flags a whole event.
6. **VAL events** get `event_rise_heavy` from the same rule, computed on their own tails (not inspected).
7. **Humidity is missing at GbgA** for {int(wx.humidity.isna().sum())} hours (all in June), including all of
   17-18 June and 30 June - 1 July 06 UTC, i.e. during both events looked at: {int(idx_main[idx_main.split_event == 'train'].t0_humidity_pct.isna().sum())}
   TRAIN main-set tails have no RH at t0, and the RH panels of the look figures are empty there.
8. **What the rise cut v2 removes differs by event** (look figures, seed {F['look_seed']}, observation only):
   in E20150617_0705 the rises come with new radar traces (rain returning), so the cut is right. In E20150701_0030 (night,
   ~02-03 local) A RISES by 0.5-1 dB right after t0_trace with no radar or gauge trace, stays up 1-2 h, then decays;
   the cut ends these tails within minutes. That can be drizzle below the radar beam (radar 78 km away) or
   condensation on wet radomes; with RH missing these hours the data here cannot decide.
9. **File names** carry no version (`ORGANIZATION.md` forbids `_v1` on active paths): the set id `{TAIL_SET}` is in
   config `waa.freeze.tail_set` and in the `tail_set` column of every frozen file.
"""
P.STEP2C_REPORT.write_text(report, encoding="utf-8")
print("\n" + "=" * 100 + "\n" + report)
print(f"done in {time.time() - clock:.0f} s")
