"""WAA_estm step 2b: fix the four issues seen in the step 2 gallery, then redo the main-set summaries.
TRAIN only for every number and figure (VAL rows get the same columns). No fitting, no shape-based selection.

Fixes (settings `waa.selection_2b` in research/config.yaml; README.md, "Step 2b rules"):
  1. tail start t0_trace = end of the last 5 min bin with ANY rain trace (radar R_max > 0 on the path, or any gauge
     within 1 km > 0) between t0 and t_first_dry; t0 if the margin has no trace. A_start_trace = mean A over
     [t0_trace, + 5 min).
  2. two tail ends: lenient = step 2 rules with the earlier dry-rise cut (3); strict = lenient, also cut at the first
     rain-trace bin after t0_trace.
  3. dry-rise cut: same detection as step 2, cut moved back to the last minute before the rise run whose smoothed
     excess was < 0.15 dB (the tail keeps that minute).
  4. noise exclusion: sub-links with offset_std (dry-weather baseline spread, baseline_uncertainty.csv) > 0.3 dB.
  5. level flags (descriptive, not exclusions): level_shift (late level far from 0), undershoot (A dips below -0.3 dB).
  6. main sets (strict / lenient) at t0_trace, counts, median curves; 7. gallery with line styles + text labels.

Writes: adds columns to outputs/products/tails_selected.csv (run step 2 first; rerunning step 2b replaces its own
columns), tables / figures step2b_*, outputs/reports/step2b_selection.md (printed at the end).
Run from the repository root (about 2 minutes):
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step2b_fix_tails.py
"""
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
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
S2, B = CFG["waa"]["selection"], CFG["waa"]["selection_2b"]
WET, DRY = 1, 0
INK, INK2 = "#1a1a1a", "#52514e"
LOCAL = pd.Timedelta(hours=2)
clock = time.time()
rng = np.random.default_rng(S2["seed"])

# Vertical markers: each has its own LINE STYLE and a text label (the reader may not separate colours).
MARKS = [("t0", "t0", (0, (1, 0))), ("t0_trace", "t0_trace", (0, (6, 3))), ("first dry", "first_dry", (0, (1, 2))),
         ("end strict", "end_strict", (0, (6, 2, 1, 2))), ("end lenient", "end_lenient", (0, (10, 2, 1, 2, 1, 2)))]
# Curves of the median plots: line style first, colour second.
CURVES = [("strict", "-", "#2a78d6"), ("lenient", "--", "#eb6834"), ("strict, no level flags", ":", "#1a1a1a"),
          ("control", "-.", "#7a7a7a")]

NEW_COLS = ["t0_trace", "margin_trace_bins", "A_start_trace_db", "t_use_end_lenient", "truncated_by_lenient",
            "usable_lenient_min", "t_use_end_strict", "truncated_by_strict", "usable_strict_min", "late_level_db",
            "late_valid_min", "level_shift", "undershoot_db", "undershoot", "dT_c", "hour_local_t0_trace",
            "offset_std_unc", "noisy_baseline", "exclusion_reason_2b", "excluded_2b", "in_main_strict",
            "in_main_lenient", "in_alt_strict", "in_alt_lenient", "in_control_2b"]


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


# =================================================================================================== load
sel = pd.read_csv(P.TAILS_SELECTED, parse_dates=["t0", "t_first_dry", "t_end", "event_start", "event_end", "t_use_end"])
sel = sel.drop(columns=[c for c in NEW_COLS if c in sel.columns])  # rerunnable: drop the old step 2b columns
long = pq.read_table(P.TAILS_1MIN, columns=["tail_id", "minutes_since_t0", "A_db", "link_state",
                                            "flag_tsl_step"]).to_pandas()
long = long.sort_values(["tail_id", "minutes_since_t0"], kind="stable").reset_index(drop=True)
unc = pd.read_csv(P.BASELINE_UNCERTAINTY).set_index("sublink_id")
wx = pd.read_csv(P.DEW_RISK_HOURLY, parse_dates=["time_utc"]).set_index("time_utc")

# --- rain trace per link and 5 min bin: radar R_max > 0 on the path OR any gauge within 1 km > 0 (references only)
d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
ra = d.store("radar_along_links_5min")
T5 = pd.DatetimeIndex(ra.time.values)
rmax5 = ra.R_max.load().to_pandas()                                 # (time, cml_id), NaN = no valid radar
city = d.store("gauges_city_5min").amount_mm.load().to_pandas()
smhi = d.store("gauges_smhi_15min").amount_mm.load().to_pandas().reindex(T5, method="ffill", limit=2)
gauge5 = pd.concat([city, smhi], axis=1)                            # 5 min amounts (SMHI: its 15 min amount x 3 bins)
lgm = d.table("link_gauge_matches")
near = lgm[lgm.distance_to_path_m <= B["gauge_max_dist_m"]].groupby("cml_id").gauge_id.apply(list).to_dict()
trace5 = {}
for c in sel.cml_id.unique():
    t = (rmax5[c].fillna(0) > 0).to_numpy().copy()                   # writable copy
    for g in near.get(c, []):
        t |= (gauge5[g].fillna(0) > 0).to_numpy()
    trace5[c] = t
print(f"loaded ({time.time() - clock:.0f} s); links with a gauge within {B['gauge_max_dist_m']} m: {len(near)}")

# =================================================================================================== per tail
bounds = np.searchsorted(long.tail_id.to_numpy(), sel.tail_id.to_numpy(), side="left")
ends = np.searchsorted(long.tail_id.to_numpy(), sel.tail_id.to_numpy(), side="right")
A_all, ST_all, TSL_all = long.A_db.to_numpy("float64"), long.link_state.to_numpy(), long.flag_tsl_step.to_numpy(bool)
rows = []
for i, r in enumerate(sel.itertuples()):
    lo, hi = bounds[i], ends[i]
    L = hi - lo
    if L == 0:
        rows.append({"tail_id": r.tail_id, "usable_lenient_min": 0, "usable_strict_min": 0, "t0_trace": r.t0,
                     "truncated_by_lenient": r.ended_by, "truncated_by_strict": r.ended_by})
        continue
    a, st, tsl = A_all[lo:hi], ST_all[lo:hi], TSL_all[lo:hi]
    i05 = T5.searchsorted(r.t0)
    fd = int(r.margin_minutes) if pd.notna(r.margin_minutes) else L          # minutes t0 -> first dry (or whole tail)
    tr = trace5[r.cml_id]

    # --- 1. t0_trace: end of the last trace bin in the margin [t0, t_first_dry)
    mb = np.nonzero(tr[i05:i05 + int(np.ceil(min(fd, L) / 5))])[0]
    k_tr = int((mb[-1] + 1) * 5) if len(mb) else 0                  # minutes since t0
    k_tr = min(k_tr, L)

    # --- 3 + 2. lenient end: tsl_step cut, earlier dry-rise cut, else t_end
    cut, by = L, r.ended_by
    k_tsl = np.nonzero(tsl)[0]
    if len(k_tsl) and k_tsl[0] < cut:
        cut, by = int(k_tsl[0]), "tsl_step"
    if pd.notna(r.margin_minutes) and fd < L:
        sm = pd.Series(a).rolling(S2["rise_smooth_min"], min_periods=S2["rise_smooth_min"] * 2 // 3).median().to_numpy()
        m = np.arange(L)
        mask = (m >= fd) & (st == DRY) & np.isfinite(sm)
        runmin = np.minimum.accumulate(np.where(mask, sm, np.inf))
        excess = sm - runmin
        k_rise = first_run(mask & (excess >= S2["rise_db"]), S2["rise_min_duration_min"])   # same detection as step 2
        if k_rise is not None:
            quiet = np.nonzero(mask[:k_rise] & (excess[:k_rise] < B["rise_quiet_db"]))[0]
            k_cut = int(quiet[-1]) + 1 if len(quiet) else fd       # keep the last quiet minute
            if k_cut < cut:
                cut, by = k_cut, "dry_rise"
    # strict end: also the first trace bin at or after t0_trace
    nxt = np.nonzero(tr[i05 + k_tr // 5:])[0]
    k_strict_trace = k_tr + int(nxt[0]) * 5 if len(nxt) else 10 ** 9
    cut_s, by_s = (k_strict_trace, "rain_trace") if k_strict_trace < cut else (cut, by)

    # --- 5. level checks
    w0, w1 = int(B["late_from_h"] * 60), min(int(B["late_to_h"] * 60) + 1, cut)
    late = a[w0:w1] if w1 > w0 else np.array([])
    n_late = int(np.isfinite(late).sum())
    late_level = float(np.nanmedian(late)) if n_late >= B["late_min_valid_min"] else np.nan
    sm15 = pd.Series(a).rolling(15, min_periods=10).median().to_numpy()
    u = sm15[k_tr:min(k_tr + int(B["undershoot_window_h"] * 60), cut)]
    undershoot = float(np.nanmin(u)) if np.isfinite(u).any() else np.nan

    rows.append({"tail_id": r.tail_id, "t0_trace": r.t0 + pd.Timedelta(minutes=k_tr), "margin_trace_bins": len(mb),
                 "A_start_trace_db": window_mean(a, k_tr, S2["a_dry_window_min"], L),
                 "t_use_end_lenient": r.t0 + pd.Timedelta(minutes=cut), "truncated_by_lenient": by,
                 "usable_lenient_min": max(cut - k_tr, 0),
                 "t_use_end_strict": r.t0 + pd.Timedelta(minutes=cut_s), "truncated_by_strict": by_s,
                 "usable_strict_min": max(cut_s - k_tr, 0),
                 "late_level_db": late_level, "late_valid_min": n_late, "undershoot_db": undershoot,
                 "_k_tr": k_tr, "_cut": cut, "_cut_s": cut_s})
print(f"per-tail pass done ({time.time() - clock:.0f} s)")

new = pd.DataFrame(rows)
sel = sel.merge(new, on="tail_id", how="left")
sel["offset_std_unc"] = sel.sublink_id.map(unc.offset_std)
thr = np.maximum(B["level_floor_db"], B["level_std_factor"] * sel.offset_std_unc.fillna(0))
sel["level_shift"] = sel.late_level_db.abs() > thr                 # NaN late level -> False (not evaluable)
sel["undershoot"] = sel.undershoot_db < B["undershoot_db"]
t_now = wx.temperature.reindex(sel.t0.dt.floor("1h")).to_numpy()
t_before = wx.temperature.reindex((sel.event_start - pd.Timedelta(hours=B["temp_before_h"])).dt.floor("1h")).to_numpy()
sel["dT_c"] = t_now - t_before
sel["hour_local_t0_trace"] = (pd.to_datetime(sel.t0_trace) + LOCAL).dt.hour

# --- 4. noise exclusion and the step 2b exclusion reason
sel["noisy_baseline"] = sel.offset_std_unc > B["max_offset_std_db"]
extra = np.where(sel.noisy_baseline, "noisy_baseline", "")
extra = np.where(sel.A_start_trace_db.isna() & (sel.exclusion_reason.fillna("") == ""),
                 np.char.add(extra.astype(str), np.where(extra != "", ";A_start_trace_nan", "A_start_trace_nan")),
                 extra)
base = sel.exclusion_reason.fillna("")
sel["exclusion_reason_2b"] = [";".join(x for x in (b0, e0) if x) for b0, e0 in zip(base, extra)]
sel["excluded_2b"] = sel.exclusion_reason_2b != ""
ok = ~sel.excluded_2b
A0 = sel.A_start_trace_db
sel["in_main_strict"] = ok & (A0 >= S2["main_a_start_db"]) & (sel.usable_strict_min >= S2["min_usable_min"])
sel["in_main_lenient"] = ok & (A0 >= S2["main_a_start_db"]) & (sel.usable_lenient_min >= S2["min_usable_min"])
sel["in_alt_strict"] = ok & (A0 >= S2["alt_a_start_db"]) & (sel.usable_strict_min >= S2["min_usable_min"])
sel["in_alt_lenient"] = ok & (A0 >= S2["alt_a_start_db"]) & (sel.usable_lenient_min >= S2["min_usable_min"])
sel["in_control_2b"] = ok & (A0 < S2["control_max_a_start_db"]) & (sel.usable_strict_min > 0)
helper = sel[["tail_id", "_k_tr", "_cut", "_cut_s"]].copy()
out = sel.drop(columns=["_k_tr", "_cut", "_cut_s"])
num = out.select_dtypes("number").columns
out.assign(**out[num].round(4)).to_csv(P.TAILS_SELECTED, index=False)

# ------------------------------------------------------------------------------------------------- TRAIN from here on
tr_ = sel[sel.split_event == "train"].copy()
ms, ml = tr_[tr_.in_main_strict], tr_[tr_.in_main_lenient]
ctrl = tr_[tr_.in_control_2b]
ms_clean = ms[~ms.level_shift & ~ms.undershoot]

# --- 4. noise: offset_std over the sub-links of the lenient main set BEFORE the noise rule
pre = tr_[(tr_.exclusion_reason.fillna("") == "") & (A0.loc[tr_.index] >= S2["main_a_start_db"])
          & (tr_.usable_lenient_min >= S2["min_usable_min"])]
pre_subs = pre.groupby("sublink_id").agg(tails=("tail_id", "size"), offset_std=("offset_std_unc", "first"))
noisy = pre_subs[pre_subs.offset_std > B["max_offset_std_db"]].sort_values("offset_std", ascending=False)
noisy = noisy.join(sel.drop_duplicates("sublink_id").set_index("sublink_id")[["L_km", "frequency_ghz", "antennas_m"]])
noisy.round(3).to_csv(P.STEP2B_NOISY)
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(pre_subs.offset_std.dropna(), bins=np.arange(0, pre_subs.offset_std.max() + 0.02, 0.02), color="#9a9a9a",
        edgecolor=INK, hatch="//", linewidth=0.5)
ax.axvline(B["max_offset_std_db"], color=INK, ls=(0, (6, 3)), lw=1.5)
ax.text(B["max_offset_std_db"], ax.get_ylim()[1] * 0.95, f"  exclusion > {B['max_offset_std_db']} dB", fontsize=8,
        va="top")
ax.set(xlabel="offset_std of the sub-link (dB, dry-weather baseline spread)", ylabel="sub-links",
       title=f"offset_std over the {len(pre_subs)} sub-links of the lenient main set (before the noise rule)")
tidy(ax)
fig.tight_layout()
fig.savefig(P.STEP2B_FIG_NOISE, dpi=110)
plt.close(fig)

# --- tail time kept / lost (kept TRAIN tails of step 2, minutes)
k2 = tr_[(tr_.exclusion_reason.fillna("") == "") & (tr_.length_min > 0)]   # sel already holds _k_tr, _cut, _cut_s
time_tab = pd.DataFrame([
    ("tail time t0 .. t_end (step 1)", k2.length_min.sum()),
    ("usable in step 2 (t0 .. t_use_end)", k2.usable_min.sum()),
    ("lenient end with the earlier rise cut (t0 .. end)", k2._cut.sum()),
    ("  removed by moving the rise cut earlier", (k2.usable_min - k2._cut).clip(lower=0).sum()),
    ("  removed by moving the start to t0_trace", np.minimum(k2._k_tr, k2._cut).sum()),
    ("usable lenient (t0_trace .. lenient end)", k2.usable_lenient_min.sum()),
    ("usable strict (t0_trace .. strict end)", k2.usable_strict_min.sum()),
    ("  lost strict vs lenient (rain trace after t0_trace)", (k2.usable_lenient_min - k2.usable_strict_min).sum()),
], columns=["quantity", "minutes"])
time_tab["share_of_step1_time"] = (time_tab.minutes / k2.length_min.sum()).round(3)
time_tab.to_csv(P.STEP2B_TIME, index=False)
trunc_l = k2.truncated_by_lenient.value_counts()
trunc_s = k2.truncated_by_strict.value_counts()


# --- 6. counts
def set_counts(df):
    return {"tails": len(df), "sub-links": df.sublink_id.nunique(), "links": df.cml_id.nunique(),
            "events": df.event_id.nunique()}


counts_sets = pd.DataFrame({
    "train tails": set_counts(tr_), "kept (step 2b exclusions)": set_counts(tr_[~tr_.excluded_2b]),
    "main strict": set_counts(ms), "main lenient": set_counts(ml),
    "main strict, no level flags": set_counts(ms_clean),
    f"alternative strict (>= {S2['alt_a_start_db']})": set_counts(tr_[tr_.in_alt_strict]),
    f"alternative lenient (>= {S2['alt_a_start_db']})": set_counts(tr_[tr_.in_alt_lenient]),
    "control": set_counts(ctrl), "step 2 main set (for reference)": set_counts(tr_[tr_.in_main_set])}).T
grp = pd.concat([g.groupby(["band", "length_class", "antennas_m"]).agg(
    tails=("tail_id", "size"), sublinks=("sublink_id", "nunique"), links=("cml_id", "nunique")).add_prefix(f"{k}_")
    for k, g in (("strict", ms), ("lenient", ml))], axis=1).fillna(0).astype(int).reset_index()
grp.to_csv(P.STEP2B_COUNTS, index=False)


# --- 5. level flags: shares, undershoot vs temperature change and time of day
def spearman(x, y):
    return x.rank().corr(y.rank())


def rho_ci(df, x, y):
    dd = df[[x, y, "cml_id"]].dropna()
    if len(dd) < 10:
        return np.nan, np.nan, np.nan, len(dd)
    groups = [g.index.to_numpy() for _, g in dd.groupby("cml_id")]   # resample LINKS, all their tails together
    reps = []
    for _ in range(S2["bootstrap_reps"]):
        ix = np.concatenate([groups[j] for j in rng.integers(0, len(groups), len(groups))])
        reps.append(spearman(dd.loc[ix, x], dd.loc[ix, y]))
    lo, hi = np.nanpercentile(reps, [2.5, 97.5])
    return round(spearman(dd[x], dd[y]), 3), round(lo, 3), round(hi, 3), len(dd)


flag_tab = pd.DataFrame({
    k: {"tails": len(g), "late level evaluable": int(g.late_level_db.notna().sum()),
        "level_shift": int(g.level_shift.sum()), "level_shift share": round(g.level_shift.mean(), 3),
        "level_shift share of evaluable": round(g.level_shift.sum() / max(g.late_level_db.notna().sum(), 1), 3),
        "undershoot": int(g.undershoot.sum()), "undershoot share": round(g.undershoot.mean(), 3),
        "both": int((g.level_shift & g.undershoot).sum())}
    for k, g in (("main strict", ms), ("main lenient", ml))}).T
rho_dT = rho_ci(ml, "dT_c", "undershoot_db")
dT_bins = pd.cut(ml.dT_c, [-np.inf, -4, -2, 0, 2, np.inf])
by_dT = ml.groupby(dT_bins, observed=False).agg(tails=("tail_id", "size"), undershoot_share=("undershoot", "mean"),
                                                median_undershoot_db=("undershoot_db", "median")).round(3)
by_dT.index = by_dT.index.astype(str)
h_bins = pd.cut(ml.hour_local_t0_trace, [0, 6, 12, 18, 24], right=False, labels=["00-06", "06-12", "12-18", "18-24"])
by_h = ml.groupby(h_bins, observed=False).agg(tails=("tail_id", "size"), undershoot_share=("undershoot", "mean"),
                                              median_undershoot_db=("undershoot_db", "median"),
                                              median_dT_c=("dT_c", "median")).round(3)
by_h.index = by_h.index.astype(str)

fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
for flag, mk, lab in ((False, "o", "not flagged"), (True, "x", f"undershoot < {B['undershoot_db']} dB")):
    g = ml[ml.undershoot == flag]
    ax[0].scatter(g.dT_c, g.undershoot_db, s=14, marker=mk, color=INK if flag else "#9a9a9a",
                  alpha=0.7 if flag else 0.4, lw=1 if flag else 0, label=lab)
ax[0].axhline(B["undershoot_db"], color=INK, ls=(0, (6, 3)), lw=1)
ax[0].set(xlabel=f"T(hour of t0) - T(event start - {B['temp_before_h']} h) (degC)", ylabel="undershoot (dB)",
          title=f"undershoot vs temperature change, rho {rho_dT[0]} [{rho_dT[1]}, {rho_dT[2]}]")
ax[0].legend(frameon=False, fontsize=8)
ax[1].bar(range(len(by_dT)), by_dT.undershoot_share, color="#bdbdbd", edgecolor=INK, hatch="//")
ax[1].set_xticks(range(len(by_dT)), by_dT.index, fontsize=8)
for j, n in enumerate(by_dT.tails):
    ax[1].text(j, by_dT.undershoot_share.iloc[j], f"n={n}", ha="center", va="bottom", fontsize=8)
ax[1].set(xlabel="temperature change (degC)", ylabel="share with undershoot", title="by temperature change")
ax[2].bar(range(len(by_h)), by_h.undershoot_share, color="#bdbdbd", edgecolor=INK, hatch="\\\\")
ax[2].set_xticks(range(len(by_h)), by_h.index)
for j, n in enumerate(by_h.tails):
    ax[2].text(j, by_h.undershoot_share.iloc[j], f"n={n}", ha="center", va="bottom", fontsize=8)
ax[2].set(xlabel="local hour of t0_trace (UTC + 2 h)", ylabel="share with undershoot", title="by time of day")
for a_ in ax:
    tidy(a_)
fig.suptitle("Undershoot (TRAIN lenient main set)", fontsize=10)
fig.tight_layout()
fig.savefig(P.STEP2B_FIG_UNDERSHOOT, dpi=110)
plt.close(fig)

# --- 6. median curves aligned at t0_trace (A and A / A_start_trace)
use = long.merge(sel[["tail_id", "split_event", "in_main_strict", "in_main_lenient", "in_control_2b", "level_shift",
                      "undershoot", "A_start_trace_db"]].merge(helper, on="tail_id"), on="tail_id")
use = use[use.split_event == "train"]
use["rel"] = use.minutes_since_t0 - use._k_tr
use = use[use.rel >= 0]
use["A_norm"] = use.A_db / use.A_start_trace_db
groups = {"strict": use[use.in_main_strict & (use.minutes_since_t0 < use._cut_s)],
          "lenient": use[use.in_main_lenient & (use.minutes_since_t0 < use._cut)],
          "strict, no level flags": use[use.in_main_strict & ~use.level_shift & ~use.undershoot
                                        & (use.minutes_since_t0 < use._cut_s)],
          "control": use[use.in_control_2b & (use.minutes_since_t0 < use._cut_s)]}
curves = pd.concat([g.groupby("rel").agg(median_A=("A_db", "median"), median_A_norm=("A_norm", "median"),
                                         count=("A_db", "count")).reset_index().assign(group=k)
                    for k, g in groups.items()])
curves.to_csv(P.STEP2B_MEDIAN_CURVES, index=False)
med_fd = (ms.t_first_dry - ms.t0_trace).dt.total_seconds().median() / 60     # median t_first_dry after t0_trace

fig, ax = plt.subplots(4, 2, figsize=(13, 15))
for row, (col, log, ylab) in enumerate((("median_A", False, "median A (dB)"), ("median_A", True, "median A (dB), log"),
                                        ("median_A_norm", False, "median A / A_start"),
                                        ("median_A_norm", True, "median A / A_start, log"))):
    for j, xmax in enumerate((90, 360)):
        a_ = ax[row, j]
        for name, ls, colr in CURVES:
            cv = curves[(curves.group == name) & (curves["count"] >= 20) & (curves.rel <= xmax)]
            if name == "control" and col == "median_A_norm":
                continue                                            # A_start ~ 0: the ratio is meaningless
            y = cv[col].to_numpy()
            x = cv.rel.to_numpy() / (1 if xmax <= 90 else 60)
            n0 = int(cv["count"].max()) if len(cv) else 0
            a_.plot(x, np.where(y > 0.01, y, np.nan) if log else y, ls=ls, color=colr, lw=2,
                    label=f"{name} (n = {n0})")
            if len(cv):                                             # direct label at the end of each curve
                yy = y[np.isfinite(y)][-1] if np.isfinite(y).any() else np.nan
                if np.isfinite(yy) and (not log or yy > 0.01):
                    a_.text(x[-1], yy, f" {name}", fontsize=7, va="center")
        if med_fd > 0:
            xv = med_fd / (1 if xmax <= 90 else 60)
            a_.axvline(xv, color=INK2, ls=(0, (1, 2)), lw=1)
            a_.text(xv, 1.0, " first dry (median)", transform=a_.get_xaxis_transform(), fontsize=7, va="top")
        a_.set_yscale("log") if log else a_.axhline(0, color=INK2, lw=0.6)
        a_.set(xlabel="minutes since t0_trace" if xmax <= 90 else "hours since t0_trace", ylabel=ylab,
               xlim=(0, (xmax if xmax <= 90 else xmax / 60) * 1.15))
        tidy(a_)
ax[0, 0].legend(frameon=False, fontsize=8)
fig.suptitle("Median tail curves aligned at t0_trace (TRAIN). Line styles: strict solid, lenient dashed, strict "
             "without level flags dotted, control dash-dot. Drawn where >= 20 tails; log y from 0.01.", fontsize=10)
fig.tight_layout()
fig.savefig(P.STEP2B_FIG_MEDIAN, dpi=100)
plt.close(fig)


# =================================================================================================== panels
def spans(states, x0, step):
    """Merge consecutive equal states into (start, end, state) spans (x in hours)."""
    out, s = [], 0
    for k in range(1, len(states) + 1):
        if k == len(states) or states[k] != states[s]:
            out.append((x0 + s * step, x0 + k * step, states[s]))
            s = k
    return out


def panel(fig, cell, r, hours, data):
    """One tail: A (top) with hatched link_state and labelled vertical markers; radar R_max + trace bins (bottom)."""
    inner = cell.subgridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
    top = fig.add_subplot(inner[0])
    bot = fig.add_subplot(inner[1], sharex=top)
    g = data[(data.tail_id == r.tail_id) & (data.minutes_since_t0 < hours * 60)]
    h = g.minutes_since_t0.to_numpy() / 60
    st5_ = g.link_state.to_numpy()[::5]
    for a_ in (top, bot):
        for x0, x1, s_ in spans(list(st5_), 0, 5 / 60):
            if s_ == WET:
                a_.axvspan(x0, x1, facecolor="#cfcfcf", edgecolor="#555555", hatch="xxx", lw=0)
            elif s_ != DRY:
                a_.axvspan(x0, x1, facecolor="none", edgecolor="#8a8a8a", hatch="///", lw=0)
    top.plot(h, g.A_db.to_numpy(), color=INK, lw=0.8)
    top.axhline(0, color=INK2, lw=0.5)
    i05 = T5.searchsorted(r.t0)
    nb = int(hours * 12)
    rm = rmax5[r.cml_id].to_numpy()[i05:i05 + nb]
    hb = np.arange(len(rm)) * 5 / 60
    bot.step(hb, rm, where="post", color=INK, lw=1)
    trb = trace5[r.cml_id][i05:i05 + nb]
    bot.plot(hb[trb] + 2.5 / 60, np.zeros(trb.sum()), "|", color=INK, ms=8, mew=1.2)
    times = {"t0": 0.0, "t0_trace": (r.t0_trace - r.t0) / pd.Timedelta(hours=1),
             "first_dry": (r.t_first_dry - r.t0) / pd.Timedelta(hours=1) if pd.notna(r.t_first_dry) else np.nan,
             "end_strict": (r.t_use_end_strict - r.t0) / pd.Timedelta(hours=1),
             "end_lenient": (r.t_use_end_lenient - r.t0) / pd.Timedelta(hours=1)}
    for k, (lab, key, ls) in enumerate(MARKS):
        x = times[key]
        if np.isfinite(x) and x <= hours:
            for a_ in (top, bot):
                a_.axvline(x, color=INK, ls=ls, lw=1.1)
            top.text(x, 0.98 - 0.11 * k, f" {lab}", transform=top.get_xaxis_transform(), fontsize=6, va="top",
                     bbox=dict(facecolor="white", edgecolor="none", pad=0.5, alpha=0.8))
    flags = [f for f, v in (("level_shift", r.level_shift), ("undershoot", r.undershoot),
                            ("noisy", r.noisy_baseline)) if v]
    top.set_title(f"{r.sublink_id} {r.event_id} | {r.L_km:.1f} km {r.frequency_ghz:.0f} GHz | A_start "
                  f"{r.A_start_trace_db:.1f} dB\nflags: {', '.join(flags) if flags else 'none'} | strict end: "
                  f"{r.truncated_by_strict}", fontsize=7)
    top.tick_params(labelsize=6, labelbottom=False)
    bot.tick_params(labelsize=6)
    bot.set_xlim(0, hours)
    bot.set_ylabel("R_max", fontsize=6)
    top.set_ylabel("A (dB)", fontsize=6)
    for a_ in (top, bot):
        for sp in ("top", "right"):
            a_.spines[sp].set_visible(False)


LEGEND = ("Vertical lines (each its own style, labelled at the top): t0 solid, t0_trace dashed, first dry dotted, "
          "end strict dash-dot, end lenient dash-dot-dot.\nHatching: /// unknown, xxx wet. Bottom: radar R_max (mm/h); "
          "| = rain-trace bin (radar > 0 or gauge within 1 km > 0). x: hours since t0.")


def panel_figure(df, n, seed, hours, ncol, path, title):
    pick = df.sample(min(n, len(df)), random_state=seed).sort_values("tail_id") if len(df) else df
    data = pq.read_table(P.TAILS_1MIN, columns=["tail_id", "minutes_since_t0", "A_db", "link_state"],
                         filters=[("tail_id", "in", pick.tail_id.tolist())]).to_pandas() if len(pick) else None
    nrow = int(np.ceil(max(len(pick), 1) / ncol))
    fig = plt.figure(figsize=(5 * ncol, 3.9 * nrow + 0.8))
    outer = fig.add_gridspec(nrow, ncol, hspace=0.55, wspace=0.22)
    for k, r in enumerate(pick.itertuples()):
        panel(fig, outer[k // ncol, k % ncol], r, hours, data)
    fig.subplots_adjust(top=1 - 1.6 / (3.9 * nrow + 0.8))           # room for the 3-line title above the panels
    fig.text(0.5, 1 - 0.2 / (3.9 * nrow + 0.8), title + "\n" + LEGEND, ha="center", va="top", fontsize=9)
    fig.savefig(path, dpi=90, bbox_inches="tight")
    plt.close(fig)
    return pick


ex_shift = panel_figure(ml[ml.level_shift], B["example_panels"], B["example_seed"], 12, 3, P.STEP2B_FIG_EX_SHIFT,
                        f"{B['example_panels']} random TRAIN lenient main-set tails flagged level_shift "
                        f"(seed {B['example_seed']}), first 12 h")
ex_under = panel_figure(ml[ml.undershoot], B["example_panels"], B["example_seed"], 6, 3, P.STEP2B_FIG_EX_UNDER,
                        f"{B['example_panels']} random TRAIN lenient main-set tails flagged undershoot "
                        f"(seed {B['example_seed']}), first 6 h")
gal = panel_figure(ms, B["gallery_n"], B["gallery_seed"], B["gallery_hours"], 4, P.STEP2B_FIG_GALLERY,
                   f"{B['gallery_n']} random TRAIN strict main-set tails (seed {B['gallery_seed']}), "
                   f"first {B['gallery_hours']} h, 1 min A")


# =================================================================================================== report
def md(df, index=True):
    df = df.reset_index() if index else df
    head = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns)
    return head + "\n" + "\n".join("| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
                                   for row in df.itertuples(index=False))


def q(s):
    s = s.replace([np.inf, -np.inf], np.nan).dropna()
    return {"n": len(s), **({k: round(float(s.quantile(v)), 2) for k, v in
                             (("p10", .1), ("p25", .25), ("median", .5), ("p75", .75), ("p90", .9))} if len(s) else {})}


rel = lambda p: Path(p).resolve().relative_to(P.REPO).as_posix()    # noqa: E731
kept2 = tr_[(tr_.exclusion_reason.fillna("") == "") & (tr_.length_min > 0)]
shift_tr = (kept2.t0_trace - kept2.t0).dt.total_seconds() / 60
cmp_start = pd.DataFrame({"A_start (t0, step 2)": q(ms.A_start_db), "A_start_trace (strict main)": q(ms.A_start_trace_db),
                          "t0_trace - t0 (min, kept tails)": q(shift_tr),
                          "t0_trace - t0 (min, strict main)": q((ms.t0_trace - ms.t0).dt.total_seconds() / 60),
                          "usable strict (min, strict main)": q(ms.usable_strict_min),
                          "usable lenient (min, lenient main)": q(ml.usable_lenient_min),
                          "late_level_db (lenient main)": q(ml.late_level_db),
                          "undershoot_db (lenient main)": q(ml.undershoot_db)}).T
# why the step 2 main set shrank: first rule that drops each step 2 main-set tail
m2 = tr_[tr_.in_main_set]
why = np.select([m2.noisy_baseline, m2.A_start_trace_db.isna(), m2.A_start_trace_db < S2["main_a_start_db"],
                 m2.usable_lenient_min < S2["min_usable_min"]],
                ["noisy_baseline", "A_start_trace NaN", f"A_start_trace < {S2['main_a_start_db']} dB",
                 f"usable lenient < {S2['min_usable_min']} min"], "still in the lenient main set")
why_tab = pd.Series(why).value_counts().rename("step 2 main-set tails").to_frame()
same_sets = bool((tr_.in_main_strict == tr_.in_main_lenient).all())
n_cut_trace = int((ms.truncated_by_strict == "rain_trace").sum())
noisy_tails_removed = int(((tr_.exclusion_reason.fillna("") == "") & tr_.noisy_baseline).sum())
noisy_main_removed = int(pre.noisy_baseline.sum())

report = f"""# WAA_estm step 2b: fixes after the step 2 gallery, main sets again (TRAIN only)

Written by `research/WAA_estm/step2b_fix_tails.py` (settings `waa.selection_2b`, plus `waa.selection` for the
step 2 rules it keeps). No fitting, no shape-based selection. TRAIN events only; VAL rows of
`{rel(P.TAILS_SELECTED)}` carry the same new columns. Run time {time.time() - clock:.0f} s.

## 1. New tail start t0_trace

t0_trace = end of the last 5 min bin with a rain trace (radar R_max > 0 on the path, or any gauge within
{B['gauge_max_dist_m']} m > 0; {len(near)} links have such a gauge) between t0 and t_first_dry; t0 if the margin has no trace.
A_start_trace = mean A over [t0_trace, t0_trace + 5 min). Kept step 2 TRAIN tails with a trace in the margin:
{int((kept2.margin_trace_bins > 0).sum())} of {len(kept2)}.

{md(cmp_start)}

## 2. + 3. Tail ends and tail time

Lenient end = earliest of the first `tsl_step`, the dry-rise cut (step 2 detection; cut moved back to the last minute
before the rise run with a smoothed excess < {B['rise_quiet_db']} dB), and t_end. Strict end = lenient end, or the first
rain-trace bin after t0_trace if earlier. Over the {len(k2)} kept step 2 TRAIN tails (`{rel(P.STEP2B_TIME)}`):

{md(time_tab, index=False)}

`truncated_by_lenient`: {trunc_l.to_dict()}. `truncated_by_strict`: {trunc_s.to_dict()}.

## 4. Noise exclusion (offset_std from dry weather)

offset_std (`baseline_uncertainty.csv`) over the {len(pre_subs)} sub-links of the lenient main set before this rule:
{q(pre_subs.offset_std)} dB.

![offset_std](../figures/{P.STEP2B_FIG_NOISE.name})

Excluded: offset_std > {B['max_offset_std_db']} dB. Sub-links removed: {len(noisy)}; their tails in that main set:
{noisy_main_removed}; kept step 2 TRAIN tails removed in total: {noisy_tails_removed}. List (`{rel(P.STEP2B_NOISY)}`):

{md(noisy.round(3))}

## 5. Level checks (flags only; nothing excluded)

late_level = median A over [t0 + {B['late_from_h']} h, t0 + {B['late_to_h']} h] inside the lenient tail, if it has >=
{B['late_min_valid_min']} valid minutes there; `level_shift` if |late_level| > max({B['level_floor_db']},
{B['level_std_factor']:g} x offset_std). undershoot = minimum of the 15 min trailing median of A within
{B['undershoot_window_h']} h after t0_trace (inside the lenient tail); flag if < {B['undershoot_db']} dB.

{md(flag_tab)}

Undershoot vs temperature change T(hour of t0) - T(event start - {B['temp_before_h']} h), lenient main set: Spearman
rho {rho_dT[0]} [{rho_dT[1]}, {rho_dT[2]}] (n = {rho_dT[3]}, bootstrap over links). By temperature change:

{md(by_dT)}

By local time of t0_trace:

{md(by_h)}

![undershoot](../figures/{P.STEP2B_FIG_UNDERSHOOT.name})

Examples, level_shift (first 12 h):

![level shift](../figures/{P.STEP2B_FIG_EX_SHIFT.name})

Examples, undershoot (first 6 h):

![undershoot examples](../figures/{P.STEP2B_FIG_EX_UNDER.name})

## 6. Main sets at t0_trace

Main = kept (step 2 exclusions + noise + A_start_trace NaN), A_start_trace >= {S2['main_a_start_db']} dB, usable length
(t0_trace .. end) >= {S2['min_usable_min']} min, for the strict and the lenient end. Control = kept, A_start_trace <
{S2['control_max_a_start_db']} dB, strict end.

{md(counts_sets)}

Strict and lenient main sets hold the same tails: {same_sets}. A dry bin needs no rain within +-30 min, so the first trace
after t_first_dry comes >= 35 min later and the strict tail almost always keeps >= {S2['min_usable_min']} min. The sets differ in
LENGTH: {n_cut_trace} of the {len(ms)} strict tails end at a rain trace (median usable strict {ms.usable_strict_min.median():.0f} min vs
lenient {ml.usable_lenient_min.median():.0f} min).

Why the step 2 main set ({len(m2)} tails) shrank (first rule that drops each tail):

{md(why_tab)}

Strict main: tails per sub-link median {ms.groupby('sublink_id').size().median():.0f} (max {ms.groupby('sublink_id').size().max()}),
per link median {ms.groupby('cml_id').size().median():.0f} (max {ms.groupby('cml_id').size().max()}). Per band (own frequency) x length
class x antenna diameters (`{rel(P.STEP2B_COUNTS)}`):

{md(grp, index=False)}

Median curves aligned at t0_trace (`{rel(P.STEP2B_MEDIAN_CURVES)}`), 0-90 min and 0-6 h, linear and log y, then the
same normalised by A_start_trace (control left out there: its A_start is ~0):

![median curves](../figures/{P.STEP2B_FIG_MEDIAN.name})

## 7. Gallery

![gallery](../figures/{P.STEP2B_FIG_GALLERY.name})

## Unclear

1. **Trace = presence, not amount.** Any R_max > 0 counts, also radar values at the 0 dBZ floor edge
   (~0.03 mm/h); radar NaN counts as no trace. Links without a gauge within 1 km rely on the radar alone.
2. **SMHI GbgA** is 15 min: one tip marks three 5 min bins as trace, which can move t0_trace up to 10 min too late.
3. **The rise cut** falls back to t_first_dry when no dry minute before the rise run had an excess < {B['rise_quiet_db']} dB.
4. **late_level** uses the lenient end and t0 (not t0_trace) as the clock, as specified; tails that end before
   t0 + {B['late_from_h']} h + {B['late_min_valid_min']} min are not evaluable (`late_level_db` NaN, `level_shift` False).
5. **The temperature change** uses one station (GbgA) and hourly values; the hour of t0 and of event start - 2 h are
   floored to the hour.
6. **offset_std** is missing for sub-links without fake events ({int(tr_.offset_std_unc.isna().sum())} TRAIN tails);
   they are not excluded by the noise rule.
7. **Strict = lenient in membership** (see section 6): the strict end changes the tail length, not the set. A strict
   set that differs would need a stricter start or a minimum strict length longer than 60 min.
8. **Undershoot looks like a level, not a dip** in the example panels: A falls to about -0.3 dB right after the rain and
   often stays there for hours (baseline too high for that event, or a level change after the rain). The flag cannot
   tell this apart from a short dip; `late_level` catches only the part that lasts beyond 6 h.
9. **A_start_trace** is measured right after the last trace; in the strict set the next trace can come soon after,
   so a short strict tail may still be dominated by the rain before.
"""
P.STEP2B_REPORT.write_text(report, encoding="utf-8")
print("\n" + "=" * 100 + "\n" + report)
print(f"done in {time.time() - clock:.0f} s")
