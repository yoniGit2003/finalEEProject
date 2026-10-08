"""WAA_estm step 2: clean and select tails with rules fixed in advance, and test whether the fast early drop of A after
the rain is WAA or missed path rain. No model fitting. No tail is selected or rejected by how much it looks like a decay.

TRAIN events only for every number, figure and test. VAL tails get the same columns in the product (exclusion,
truncation, descriptors), but they are never counted, plotted or used here.

Rules (settings: `waa.selection` in research/config.yaml; README.md, "Step 2 rules"):
  1. exclusions (reason recorded): wet t0, length 0, A_start NaN, use_in_waa_main False, status != good, NaN share > 10 %
  2. truncation: the first tsl_step minute (cut just before it), or the first sustained dry-weather rise (15 min
     trailing median of A >= 0.6 dB above its running minimum since t_first_dry, dry minutes, held >= 30 min; cut at
     the start of that run), else t_end. -> t_use_end, truncated_by
  3. main set: not excluded, A_start >= 1.0 dB, usable length >= 60 min (alternative: A_start >= 0.6 dB, counts only);
     control group: not excluded, A_start < 0.5 dB
  4. fast-drop test: fast_drop = A_start - A_dry (A_dry = mean A over [t_first_dry, + 5 min)) vs link length L, margin
     rain, and sub-link 1 vs 2; Spearman rho with a 95 % CI from a bootstrap over LINKS
  5. shape descriptors (descriptive only): A(t0 + k) / A_start, time to A_start / 2 and A_start / e, plateau
  6. median tail curves main set vs control; 7. gallery of 24 random main-set tails

Reads only this module's step 1 products (tails_index.csv, tails_1min.parquet).
Run from the repository root (about 1-2 minutes):
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step2_select_tails.py
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

sys.stdout.reconfigure(encoding="utf-8")
MODULE = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE))
import WAA_estm_paths as P

S = yaml.safe_load(open(P.CONFIG_FILE, encoding="utf-8"))["waa"]["selection"]
WET, DRY, UNKNOWN = 1, 0, -1
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]                    # categorical slots 1-4, fixed order
INK2 = "#52514e"
clock = time.time()
rng = np.random.default_rng(S["seed"])


def tidy(ax):
    """Recessive frame: no top / right spines, light grid."""
    ax.grid(alpha=0.25)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


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


# =================================================================================================== load
idx = pd.read_csv(P.TAILS_INDEX, parse_dates=["t0", "t_first_dry", "t_end", "wet_start", "event_start", "event_end",
                                                "last_wet_end"])
long = pq.read_table(P.TAILS_1MIN, columns=["tail_id", "minutes_since_t0", "A_db", "link_state", "flag_tsl_step",
                                            "R_mean_mmh", "R_max_mmh"]).to_pandas()
long = long.sort_values(["tail_id", "minutes_since_t0"], kind="stable").reset_index(drop=True)
print(f"loaded {len(idx)} tails, {len(long)} tail minutes ({time.time() - clock:.0f} s)")

# =================================================================================================== 1. exclusions
reasons = {
    "wet_t0": idx.state_at_t0 == WET,
    "length_0": idx.length_min == 0,
    "A_start_nan": idx.A_start_db.isna(),
    "not_use_in_waa_main": ~idx.use_in_waa_main.astype(bool),
    "status_not_good": idx.status != "good",
    "nan_share_gt_max": idx.frac_nan > S["max_nan_share"],
}
idx["exclusion_reason"] = [";".join(k for k, m in reasons.items() if m.iloc[i]) for i in range(len(idx))]
idx["excluded"] = idx.exclusion_reason != ""

# =================================================================================================== 2 + 4 + 5. per tail
# One pass over every tail with minutes (train AND val, so val carries the same columns): truncation, A_dry,
# margin radar, shape descriptors. Nothing computed here decides membership except the pre-defined rules.
bounds = np.searchsorted(long.tail_id.to_numpy(), idx.tail_id.to_numpy(), side="left")
ends = np.searchsorted(long.tail_id.to_numpy(), idx.tail_id.to_numpy(), side="right")
A_all, ST_all = long.A_db.to_numpy("float64"), long.link_state.to_numpy()
TSL_all, RM_all = long.flag_tsl_step.to_numpy(bool), long.R_mean_mmh.to_numpy("float64")
RX_all = long.R_max_mmh.to_numpy("float64")
rows = []
for i, r in enumerate(idx.itertuples()):
    lo, hi = bounds[i], ends[i]
    L = hi - lo                                                     # minutes in the product (= length_min)
    out = {"tail_id": r.tail_id}
    if L == 0:
        rows.append(out | {"t_use_end": r.t_end, "truncated_by": r.ended_by, "usable_min": 0})
        continue
    a, st, tsl = A_all[lo:hi], ST_all[lo:hi], TSL_all[lo:hi]
    fd = int((r.t_first_dry - r.t0) / pd.Timedelta("1min")) if pd.notna(r.t_first_dry) else None  # minutes t0 -> dry

    # --- truncation candidates (minutes since t0; the tail is used over [0, cut))
    cut, by = L, r.ended_by
    k_tsl = np.nonzero(tsl)[0]
    if len(k_tsl) and k_tsl[0] < cut:
        cut, by = int(k_tsl[0]), "tsl_step"
    if fd is not None and fd < L:
        sm = pd.Series(a).rolling(S["rise_smooth_min"], min_periods=S["rise_smooth_min"] * 2 // 3).median().to_numpy()
        m = np.arange(L)
        mask = (m >= fd) & (st == DRY) & np.isfinite(sm)
        runmin = np.minimum.accumulate(np.where(mask, sm, np.inf))
        k_rise = first_run(mask & (sm - runmin >= S["rise_db"]), S["rise_min_duration_min"])
        if k_rise is not None and k_rise < cut:
            cut, by = k_rise, "dry_rise"
    out.update({"t_use_end": r.t0 + pd.Timedelta(minutes=cut), "truncated_by": by, "usable_min": cut})

    # --- fast-drop quantities (A_dry inside the usable part only)
    out["margin_minutes"] = fd if fd is not None else np.nan
    out["A_dry_db"] = window_mean(a, fd, S["a_dry_window_min"], cut) if fd is not None else np.nan
    out["margin_R_mean_sum_mm"] = float(np.nansum(RM_all[lo:lo + min(fd if fd is not None else L, L)]) / 60)
    rx = RX_all[lo + fd:lo + cut] if fd is not None and fd < cut else np.array([])   # radar after t_first_dry
    out["R_max_after_dry_max_mmh"] = float(np.nanmax(rx)) if np.isfinite(rx).any() else np.nan   # light rain again?

    # --- shape descriptors (usable part only; descriptive)
    w = S["value_window_min"]
    for k in S["shape_minutes"]:
        out[f"A_{k}min_db"] = window_mean(a, k, w, cut)
    for k in S["pair_minutes"]:
        out.setdefault(f"A_{k}min_db", window_mean(a, k, w, cut))
    sm5 = pd.Series(a[:cut]).rolling(S["crossing_smooth_min"], min_periods=1).mean().to_numpy()
    for name, level in (("t_half_min", r.A_start_db / 2), ("t_e_min", r.A_start_db / np.e)):
        hit = np.nonzero(sm5 < level)[0] if np.isfinite(r.A_start_db) else []
        out[name] = int(hit[0]) if len(hit) else np.nan             # NaN = not reached before t_use_end (censored)
    n = S["plateau_min_duration_min"]
    s_ = pd.Series(a[:cut])
    mx, mn = s_.rolling(n, min_periods=n).max().to_numpy(), s_.rolling(n, min_periods=n).min().to_numpy()
    pl = np.nonzero((mx - mn <= S["plateau_band_db"]) & (mn > S["plateau_min_a_db"]))[0]
    out["plateau"] = bool(len(pl))
    out["plateau_start_min"] = int(pl[0] - n + 1) if len(pl) else np.nan
    rows.append(out)
print(f"per-tail pass done ({time.time() - clock:.0f} s)")

sel = idx.merge(pd.DataFrame(rows), on="tail_id", how="left")
sel["fast_drop_db"] = sel.A_start_db - sel.A_dry_db
sel["fast_drop_share"] = sel.fast_drop_db / sel.A_start_db
for k in S["shape_minutes"]:
    sel[f"ratio_{k}min"] = sel[f"A_{k}min_db"] / sel.A_start_db
sel["L_km"] = sel.length_m / 1000
sel["band"] = np.where(sel.frequency_ghz < 25, "<25GHz", ">=25GHz")   # the sub-link's OWN frequency
d0, d1 = sel.site_0_antenna_diameter_m, sel.site_1_antenna_diameter_m
sel["antennas_m"] = [f"{min(x, y):g}/{max(x, y):g}" for x, y in zip(d0, d1)]
sel["margin_rain"] = np.where(sel.margin_R_max_max_mmh > 0, "R_max > 0", "R_max = 0")
ok = ~sel.excluded
sel["in_main_set"] = ok & (sel.A_start_db >= S["main_a_start_db"]) & (sel.usable_min >= S["min_usable_min"])
sel["in_alt_set"] = ok & (sel.A_start_db >= S["alt_a_start_db"]) & (sel.usable_min >= S["min_usable_min"])
sel["in_control"] = ok & (sel.A_start_db < S["control_max_a_start_db"]) & (sel.usable_min > 0)
num = sel.select_dtypes('number').columns                          # round numbers only (keeps the file small)
sel.assign(**sel[num].round(4)).to_csv(P.TAILS_SELECTED, index=False)

# ------------------------------------------------------------------------------------------------- TRAIN from here on
tr = sel[sel.split_event == "train"].copy()
main = tr[tr.in_main_set].copy()
ctrl = tr[tr.in_control].copy()

# =================================================================================================== counts
excl_counts = pd.Series({k: int((m & (sel.split_event == "train")).sum()) for k, m in reasons.items()})
trunc = (tr[~tr.excluded & (tr.length_min > 0)]
         .assign(removed_min=lambda x: x.length_min - x.usable_min)
         .groupby("truncated_by").agg(tails=("tail_id", "size"), removed_min=("removed_min", "sum"),
                                      median_removed_min=("removed_min", "median")))
trunc["share_of_all_tail_time"] = (trunc.removed_min / tr[~tr.excluded].length_min.sum()).round(4)
trunc.to_csv(P.STEP2_TRUNCATION)


def set_counts(df):
    return {"tails": len(df), "sub-links": df.sublink_id.nunique(), "links": df.cml_id.nunique(),
            "events": df.event_id.nunique()}


counts_sets = pd.DataFrame({
    "train tails": set_counts(tr), "not excluded": set_counts(tr[~tr.excluded]),
    f"main (A_start >= {S['main_a_start_db']})": set_counts(main),
    f"alternative (A_start >= {S['alt_a_start_db']})": set_counts(tr[tr.in_alt_set]),
    f"control (A_start < {S['control_max_a_start_db']})": set_counts(ctrl)}).T
grp = (main.groupby(["band", "length_class", "antennas_m"])
       .agg(tails=("tail_id", "size"), sublinks=("sublink_id", "nunique"), links=("cml_id", "nunique"))
       .reset_index())
grp.to_csv(P.STEP2_COUNTS, index=False)
per_sub, per_link = main.groupby("sublink_id").size(), main.groupby("cml_id").size()


# =================================================================================================== Spearman + CI
def spearman(x, y):
    return x.rank().corr(y.rank())


def rho_ci(df, x, y, label):
    """Spearman rho of x vs y, 95 % CI from resampling LINKS (all their tails go together)."""
    d = df[[x, y, "cml_id"]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(d) < 10:
        return {"test": label, "x": x, "y": y, "n": len(d), "links": d.cml_id.nunique(), "rho": np.nan,
                "ci_low": np.nan, "ci_high": np.nan}
    groups = [g.index.to_numpy() for _, g in d.groupby("cml_id")]
    reps = []
    for _ in range(S["bootstrap_reps"]):
        pick = rng.integers(0, len(groups), len(groups))
        ix = np.concatenate([groups[j] for j in pick])
        reps.append(spearman(d.loc[ix, x], d.loc[ix, y]))
    lo, hi = np.nanpercentile(reps, [2.5, 97.5])
    return {"test": label, "x": x, "y": y, "n": len(d), "links": d.cml_id.nunique(),
            "rho": round(spearman(d[x], d[y]), 3), "ci_low": round(lo, 3), "ci_high": round(hi, 3)}


tests = []
for band in ("<25GHz", ">=25GHz"):
    m = main[main.band == band]
    tests += [rho_ci(m, "L_km", "fast_drop_db", f"4a fast_drop vs L, {band}"),
              rho_ci(m, "L_km", "fast_drop_share", f"4a fast_drop / A_start vs L, {band}"),
              rho_ci(m, "L_km", "wet_A_max_db", f"4a reference: wet_A_max vs L, {band}"),
              rho_ci(m.assign(k=m.wet_A_max_db / m.L_km), "L_km", "k", f"4a reference: wet_A_max / L vs L, {band}")]
tests += [rho_ci(main, "margin_R_mean_sum_mm", "fast_drop_db", "4c fast_drop vs margin_R_mean_sum"),
          rho_ci(main, "margin_R_max_max_mmh", "fast_drop_db", "4c fast_drop vs margin_R_max_max"),
          rho_ci(main, "margin_minutes", "fast_drop_db", "4c fast_drop vs margin_min")]

# 4d: sub-link 1 vs 2 of the same link and event, both in the main set
pv = main.assign(sl=main.sublink_id.str[-1]).pivot_table(
    index=["cml_id", "event_id"], columns="sl",
    values=["A_start_db", "fast_drop_db"] + [f"A_{k}min_db" for k in S["pair_minutes"]])
pair_vars = ["A_start_db", "fast_drop_db"] + [f"A_{k}min_db" for k in S["pair_minutes"]]
pairs = pd.DataFrame({"cml_id": pv.index.get_level_values(0)})
for v in pair_vars:
    if (v, "1") in pv.columns and (v, "2") in pv.columns:
        pairs[f"{v}_1"], pairs[f"{v}_2"] = pv[(v, "1")].to_numpy(), pv[(v, "2")].to_numpy()
pairs = pairs.dropna(subset=["A_start_db_1", "A_start_db_2"])
for v in pair_vars:
    tests.append(rho_ci(pairs, f"{v}_1", f"{v}_2", f"4d sub-link 1 vs 2: {v}"))
tests = pd.DataFrame(tests)
tests.to_csv(P.STEP2_CORRELATIONS, index=False)

# =================================================================================================== median curves
use = long.merge(sel[["tail_id", "usable_min", "in_main_set", "in_control", "split_event", "margin_rain"]],
                 on="tail_id")
use = use[(use.split_event == "train") & (use.minutes_since_t0 < use.usable_min)]
groups = {"main set": use[use.in_main_set], "control": use[use.in_control],
          "main, margin R_max = 0": use[use.in_main_set & (use.margin_rain == "R_max = 0")],
          "main, margin R_max > 0": use[use.in_main_set & (use.margin_rain == "R_max > 0")]}
curves = pd.concat([g.groupby("minutes_since_t0").A_db.agg(["median", "count"]).reset_index().assign(group=k)
                    for k, g in groups.items()])
curves.to_csv(P.STEP2_MEDIAN_CURVES, index=False)
med_fd = main.margin_minutes.median()


def draw_curves(ax, names, colors, xmax_h, log):
    for name, col in zip(names, colors):
        cv = curves[(curves.group == name) & (curves["count"] >= 20) & (curves.minutes_since_t0 <= xmax_h * 60)]
        y = cv["median"].to_numpy()
        ax.plot(cv.minutes_since_t0 / (60 if xmax_h > 2 else 1), np.where(y > 0.01, y, np.nan) if log else y,
                color=col, lw=2, label=f"{name} (n = {int(cv['count'].max()) if len(cv) else 0})")
    ax.axvline(med_fd / (60 if xmax_h > 2 else 1), color=INK2, ls=":", lw=1)
    if log:
        ax.set_yscale("log")
    else:
        ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xlabel("hours since t0" if xmax_h > 2 else "minutes since t0")
    ax.set_ylabel("median A (dB)" + (", log" if log else ""))
    tidy(ax)


fig, ax = plt.subplots(2, 2, figsize=(13, 8))
for j, (xh, title) in enumerate(((12, "0-12 h"), (1.5, "0-90 min"))):
    for i, log in enumerate((False, True)):
        draw_curves(ax[i, j], ["main set", "control"], C[:2], xh, log)
        ax[i, j].set_title(f"{title}, {'log' if log else 'linear'} y", fontsize=10)
ax[0, 0].legend(frameon=False, fontsize=8)
fig.suptitle(f"Median tail curve, TRAIN: main set vs control (A_start < {S['control_max_a_start_db']} dB). "
             f"Dotted: median t_first_dry ({med_fd:.0f} min). Drawn where >= 20 tails; log y from 0.01 dB.",
             fontsize=10)
fig.tight_layout()
fig.savefig(P.STEP2_FIG_MEDIAN, dpi=110)
plt.close(fig)

# =================================================================================================== 4a figure
edges = np.array(S["length_bins_km"])
mids = np.sqrt(edges[:-1] * edges[1:])


def binned(df, col):
    b = pd.cut(df.L_km, edges)
    g = df.groupby(b, observed=False)[col]
    return g.median().to_numpy(), g.size().to_numpy()


fig, ax = plt.subplots(2, 2, figsize=(13, 9))
for i, band in enumerate(("<25GHz", ">=25GHz")):
    m = main[main.band == band]
    a_ = ax[i, 0]
    pos = m[m.fast_drop_db > 0.01]                                  # log axis: values <= 0.01 dB not drawn
    a_.scatter(pos.L_km, pos.fast_drop_db, s=8, color=C[0], alpha=0.35, lw=0, label="fast_drop (tails, > 0.01 dB)")
    for col, cc, lab in (("fast_drop_db", C[0], "fast_drop, binned median"),
                         ("wet_A_max_db", C[1], "wet_A_max, binned median (reference)")):
        med, n = binned(m, col)
        keep = (n >= 5) & (med > 0)
        a_.plot(mids[keep], med[keep], "o-", color=cc, lw=2, ms=6, label=lab)
    xs = np.array([0.15, 15])
    ref, nref = binned(m, "wet_A_max_db")
    if ((nref >= 5) & (ref > 0)).any():
        j0 = np.nonzero((nref >= 5) & (ref > 0))[0][0]
        a_.plot(xs, ref[j0] * xs / mids[j0], color=INK2, ls="--", lw=1, label="slope 1 (A proportional to L)")
    a_.set(xscale="log", yscale="log", ylim=(0.01, 60), xlabel="link length L (km)", ylabel="dB",
           title=f"{band}: fast_drop and wet_A_max vs L (n = {len(m)}; fast_drop <= 0.01: {(m.fast_drop_db <= 0.01).sum()})")
    a_.legend(frameon=False, fontsize=7)
    tidy(a_)
    b_ = ax[i, 1]
    b_.scatter(m.L_km, m.fast_drop_share, s=8, color=C[0], alpha=0.35, lw=0)
    med, n = binned(m, "fast_drop_share")
    b_.plot(mids[n >= 5], med[n >= 5], "o-", color=C[0], lw=2, ms=6, label="binned median")
    b_.axhline(0, color=INK2, lw=0.8)
    b_.set(xscale="log", ylim=(-1, 1.2), xlabel="link length L (km)", ylabel="fast_drop / A_start",
           title=f"{band}: share of A_start lost in the margin")
    b_.legend(frameon=False, fontsize=7)
    tidy(b_)
fig.suptitle("4a. Expected: path rain grows with L (slope 1 on log-log), WAA does not (flat). TRAIN main set.", fontsize=10)
fig.tight_layout()
fig.savefig(P.STEP2_FIG_LENGTH, dpi=110)
plt.close(fig)

# =================================================================================================== 4b / 4c figure
fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
for k, (name, g) in enumerate(main.groupby("margin_rain")):
    ax[0].hist(g.fast_drop_db.clip(-2, 8), bins=np.arange(-2, 8.25, 0.25), color=C[k], alpha=0.6,
               label=f"{name} (n = {len(g)}, median {g.fast_drop_db.median():.2f} dB)")
ax[0].set(xlabel="fast_drop (dB, clipped to [-2, 8])", ylabel="tails", title="4b. fast_drop by margin radar")
ax[0].legend(frameon=False, fontsize=8)
draw_curves(ax[1], ["main, margin R_max = 0", "main, margin R_max > 0"], C[:2], 1.5, False)
ax[1].set_title("4b. median tail curve (0-90 min)")
ax[1].legend(frameon=False, fontsize=8)
ax[2].scatter(main.margin_R_mean_sum_mm, main.fast_drop_db, s=8, color=C[0], alpha=0.4, lw=0)
ax[2].set_xscale("symlog", linthresh=1e-3)
ax[2].set(xlabel="margin_R_mean_sum (mm, symlog)", ylabel="fast_drop (dB)", title="4c. fast_drop vs radar rain in the margin")
for a_ in (ax[0], ax[2]):
    tidy(a_)
fig.tight_layout()
fig.savefig(P.STEP2_FIG_MARGIN, dpi=110)
plt.close(fig)

# =================================================================================================== 4d figure
fig, ax = plt.subplots(1, len(pair_vars), figsize=(4.2 * len(pair_vars), 4.2))
for a_, v in zip(ax, pair_vars):
    x, y = pairs[f"{v}_1"], pairs[f"{v}_2"]
    a_.scatter(x, y, s=10, color=C[0], alpha=0.5, lw=0)
    lim = np.nanpercentile(np.r_[x, y], [1, 99]) if len(pairs) else (0, 1)
    a_.plot(lim, lim, color=INK2, ls="--", lw=1)
    t = tests[tests.test == f"4d sub-link 1 vs 2: {v}"].iloc[0]
    a_.set(xlabel=f"{v}, sub-link 1", ylabel=f"{v}, sub-link 2",
           title=f"rho {t.rho} [{t.ci_low}, {t.ci_high}], n = {t.n}", xlim=lim, ylim=lim)
    tidy(a_)
fig.suptitle("4d. Same link and event, both sub-links in the TRAIN main set. Dashed: 1:1.", fontsize=10)
fig.tight_layout()
fig.savefig(P.STEP2_FIG_PAIRS, dpi=110)
plt.close(fig)

# =================================================================================================== 5 figure
fig, ax = plt.subplots(1, 3, figsize=(16, 4.5))
ax[0].boxplot([main[f"ratio_{k}min"].dropna().clip(-1, 2) for k in S["shape_minutes"]],
              tick_labels=[f"{k}" for k in S["shape_minutes"]], showfliers=False)
ax[0].axhline(0, color=INK2, lw=0.8)
ax[0].set(xlabel="minutes since t0", ylabel="A(t0 + k) / A_start", title="5. decay ratios (box = p25-p75)")
for k, col in enumerate(("t_half_min", "t_e_min")):
    ax[1].hist(main[col].dropna().clip(0, 360), bins=np.arange(0, 365, 5), color=C[k], alpha=0.6,
               label=f"{col} (median {main[col].median():.0f} min; not reached {main[col].isna().sum()})")
ax[1].set(xlabel="minutes (clipped at 360)", ylabel="tails", title="5. time to A_start / 2 and A_start / e")
ax[1].legend(frameon=False, fontsize=8)
lb = pd.cut(main.L_km, edges)
pl_L = main.groupby(lb, observed=False).plateau.agg(["mean", "size"])
keep = pl_L["size"] >= 5
ax[2].plot(mids[keep.to_numpy()], pl_L["mean"][keep], "o-", color=C[0], lw=2, label="by link length")
pl_m = main.groupby("margin_rain").plateau.mean()
for k, (name, v) in enumerate(pl_m.items()):
    ax[2].axhline(v, color=C[1 + k], ls="--", lw=1.2, label=f"all L, margin {name}: {v:.0%}")
ax[2].set(xscale="log", xlabel="link length L (km)", ylabel="share of tails with a plateau",
          title="5. plateau share")
ax[2].legend(frameon=False, fontsize=8)
for a_ in ax:
    tidy(a_)
fig.tight_layout()
fig.savefig(P.STEP2_FIG_SHAPE, dpi=110)
plt.close(fig)

# =================================================================================================== 7 gallery
gal = main.sample(min(S["gallery_n"], len(main)), random_state=S["gallery_seed"]).sort_values("tail_id")
H = S["gallery_hours"] * 60
g_long = pq.read_table(P.TAILS_1MIN, columns=["tail_id", "minutes_since_t0", "A_db", "link_state", "R_max_mmh"],
                       filters=[("tail_id", "in", gal.tail_id.tolist())]).to_pandas()
fig = plt.figure(figsize=(20, 22))
outer = fig.add_gridspec(6, 4, hspace=0.45, wspace=0.22)
for n, r in enumerate(gal.itertuples()):
    inner = outer[n // 4, n % 4].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
    top = fig.add_subplot(inner[0])
    bot = fig.add_subplot(inner[1], sharex=top)
    g = g_long[(g_long.tail_id == r.tail_id) & (g_long.minutes_since_t0 < H)].sort_values("minutes_since_t0")
    h = g.minutes_since_t0.to_numpy() / 60
    for a_ in (top, bot):
        for x0, s_ in zip(h[::5], g.link_state.to_numpy()[::5]):
            if s_ != DRY:
                a_.axvspan(x0, x0 + 5 / 60, color="#4c8ed9" if s_ == WET else "#999999",
                           alpha=0.25 if s_ == WET else 0.15, lw=0)
        if pd.notna(r.margin_minutes):
            a_.axvline(r.margin_minutes / 60, color=C[2], ls="--", lw=1.2)
        if r.usable_min < H:
            a_.axvline(r.usable_min / 60, color=C[1], ls="-", lw=1.5)
    top.plot(h, g.A_db.to_numpy(), color="#1a1a1a", lw=0.8)
    top.axhline(0, color=INK2, lw=0.6)
    top.set_title(f"{r.sublink_id} {r.event_id} | {r.L_km:.1f} km {r.frequency_ghz:.0f} GHz | {r.truncated_by}",
                  fontsize=8)
    top.tick_params(labelsize=7, labelbottom=False)
    bot.step(h, g.R_max_mmh.to_numpy(), where="post", color=INK2, lw=1)
    bot.tick_params(labelsize=7)
    bot.set_xlim(0, S["gallery_hours"])
    if n % 4 == 0:
        top.set_ylabel("A (dB)", fontsize=8)
        bot.set_ylabel("R_max", fontsize=7)
    for a_ in (top, bot):
        for sp in ("top", "right"):
            a_.spines[sp].set_visible(False)
fig.suptitle(f"{len(gal)} random TRAIN main-set tails (seed {S['gallery_seed']}), first {S['gallery_hours']} h, 1 min A. "
             "Green dashed: t_first_dry; orange: t_use_end (if inside); shading: link_state wet (blue) / unknown "
             "(grey). Bottom: radar R_max (mm/h). x: hours since t0.", fontsize=11, y=0.995)
fig.subplots_adjust(top=0.97)
fig.savefig(P.STEP2_FIG_GALLERY, dpi=90, bbox_inches="tight")
plt.close(fig)


# =================================================================================================== report
def md(df, index=True):
    df = df.reset_index() if index else df
    head = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns)
    return head + "\n" + "\n".join("| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
                                   for row in df.itertuples(index=False))


def q(s):
    s = s.replace([np.inf, -np.inf], np.nan).dropna()
    if not len(s):
        return {"n": 0}
    return {"n": len(s), **{k: round(float(s.quantile(v)), 2) for k, v in
                            (("p10", .1), ("p25", .25), ("median", .5), ("p75", .75), ("p90", .9))}}


rel = lambda p: Path(p).resolve().relative_to(P.REPO).as_posix()    # noqa: E731
fd_by_margin = main.groupby("margin_rain").fast_drop_db.apply(lambda s: pd.Series(q(s))).unstack()
shape_tab = pd.DataFrame({**{f"ratio_{k}min": q(main[f"ratio_{k}min"]) for k in S["shape_minutes"]},
                          "t_half_min": q(main.t_half_min), "t_e_min": q(main.t_e_min),
                          "fast_drop_db": q(main.fast_drop_db), "fast_drop_share": q(main.fast_drop_share),
                          "margin_minutes": q(main.margin_minutes)}).T
pl_tab = pd.DataFrame({"plateau share": main.groupby(lb, observed=False).plateau.mean().round(3),
                       "tails": main.groupby(lb, observed=False).size()})
pl_tab.index = pl_tab.index.astype(str)
val_n = int((sel.split_event == "val").sum())

report = f"""# WAA_estm step 2: clean / select tails, fast-drop test (TRAIN only)

Written by `research/WAA_estm/step2_select_tails.py` (settings `waa.selection` in `research/config.yaml`, fixed
before the run). No fitting. No tail is kept or dropped by its shape. Every number below uses **TRAIN events only**.
The {val_n} VAL tails carry the same columns in `{rel(P.TAILS_SELECTED)}` but are not looked at here.
Run time {time.time() - clock:.0f} s.

## 1. Exclusions (pre-defined, non-shape)

TRAIN tails: {len(tr)}; excluded: {int(tr.excluded.sum())}; kept: {int((~tr.excluded).sum())}. A tail can have several
reasons (`exclusion_reason`, `;`-separated). Tails per reason:

{md(excl_counts.rename('train tails').to_frame())}

## 2. Truncation of the tail end

On the {int((~tr.excluded & (tr.length_min > 0)).sum())} kept TRAIN tails: `truncated_by` = the rule that set
`t_use_end` (`tsl_step`: just before the first TSL jump; `dry_rise`: start of the first run of >= {S['rise_min_duration_min']}
min where the {S['rise_smooth_min']} min trailing median of A is >= {S['rise_db']} dB above its running minimum since
t_first_dry, dry minutes only; otherwise the step 1 `ended_by`). Tail time removed (`{rel(P.STEP2_TRUNCATION)}`):

{md(trunc)}

## 3. Main set

Main set = kept, A_start >= {S['main_a_start_db']} dB, usable length >= {S['min_usable_min']} min.

{md(counts_sets)}

Main-set tails per sub-link: median {per_sub.median():.0f}, max {per_sub.max()} ({len(per_sub)} sub-links); per link:
median {per_link.median():.0f}, max {per_link.max()}. Per band (the sub-link's own frequency) x length class x
antenna diameters (smaller / larger end, m) (`{rel(P.STEP2_COUNTS)}`):

{md(grp, index=False)}

## 4. Fast-drop test: is the early drop WAA or missed path rain?

fast_drop = A_start - A_dry, A_dry = mean A over [t_first_dry, t_first_dry + {S['a_dry_window_min']} min) (NaN if that
window is not inside the usable tail). The margin (t0 .. t_first_dry) is the link_state `unknown` time: light rain
or rain within 30 min. **If the early drop is path rain**, fast_drop should grow with L (like wet_A_max, slope 1 on
log-log), with margin radar rain, and agree between the two sub-links only as much as rain does. **If it is WAA**, it should
not grow with L.

Distributions (main set):

{md(shape_tab.loc[['fast_drop_db', 'fast_drop_share', 'margin_minutes']])}

Spearman rho with 95 % CI from {S['bootstrap_reps']} bootstrap resamples of LINKS (`{rel(P.STEP2_CORRELATIONS)}`):

{md(tests, index=False)}

**4a.** fast_drop and the reference wet_A_max vs L, per band:

![4a](../figures/{P.STEP2_FIG_LENGTH.name})

**4b / 4c.** fast_drop by margin radar (R_max in the margin = 0 vs > 0):

{md(fd_by_margin)}

![4b 4c](../figures/{P.STEP2_FIG_MARGIN.name})

**4d.** Sub-link 1 vs 2 of the same link and event (both in the main set): {len(pairs)} pairs.

![4d](../figures/{P.STEP2_FIG_PAIRS.name})

## 5. Shape descriptors (descriptive only; not used for selection)

A(t0 + k) = mean over [k, k + {S['value_window_min']}) min, inside the usable tail (NaN after t_use_end). Crossing times
on a {S['crossing_smooth_min']} min trailing mean; NaN = not reached before t_use_end (censored, so the medians are lower
bounds). Plateau: a {S['plateau_min_duration_min']} min window with max - min <= {S['plateau_band_db']} dB and min >
{S['plateau_min_a_db']} dB.

{md(shape_tab.drop(['fast_drop_db', 'fast_drop_share', 'margin_minutes']))}

Plateau tails: {int(main.plateau.sum())} of {len(main)} ({main.plateau.mean():.0%}). By margin radar:
{ {k: f'{v:.0%}' for k, v in pl_m.items()} }. By link length:

{md(pl_tab)}

![shape](../figures/{P.STEP2_FIG_SHAPE.name})

## 6. Median tail curve: main set vs control

Control = kept TRAIN tails with A_start < {S['control_max_a_start_db']} dB ({len(ctrl)} tails). Usable part only.
Dotted line = median t_first_dry of the main set ({med_fd:.0f} min). `{rel(P.STEP2_MEDIAN_CURVES)}`.

![median curves](../figures/{P.STEP2_FIG_MEDIAN.name})

## 7. Gallery for a manual check

![gallery](../figures/{P.STEP2_FIG_GALLERY.name})

## Unclear

1. **Rise start.** The `dry_rise` cut is at the first minute of the >= {S['rise_min_duration_min']} min run in which
   the excess is >= {S['rise_db']} dB, so the cut comes after the rise began (by the time the median needs to climb
   {S['rise_db']} dB). An earlier cut (where the excess first left ~0) would need another threshold.
2. **The rise is searched only on dry minutes from t_first_dry on.** A rise inside the margin, or during `unknown`
   minutes later on, is not cut.
3. **A_dry needs a dry bin inside the usable tail**: main-set tails with A_dry NaN: {int(main.A_dry_db.isna().sum())}
   (t_first_dry after t_use_end or no valid A there). They have no fast_drop.
4. **Band** is the sub-link's own frequency here (< / >= 25 GHz), not the release's `frequency_band` (sub-link 1 of
   the link). Sub-links whose band differs from the release band: {int((tr.band != tr.frequency_band).sum())} train tails.
5. **wet_A_max is not a clean path-rain measure**: it spans the link's whole wet phase in the network event (several
   showers in multi-run events) and contains WAA itself. It is the reference curve only for the L scaling.
6. **Pairs (4d)** need both sub-links in the main set, so pairs where one direction has A_start < {S['main_a_start_db']} dB
   are missing (selection on A_start of each direction).
7. **Crossing times and plateau** use the 1 min A, which still carries the 0.3 dB RSL quantization; a {S['plateau_band_db']}
   dB band is half a step.
8. **Light rain after t_first_dry stays in the tail.** No rule cuts `unknown` bins or 0 < R_max < 0.1 mm/h after the
   first dry bin (only a WET bin ends a tail). Main-set tails with R_max > 0 somewhere between t_first_dry and
   t_use_end: {int((main.R_max_after_dry_max_mmh > 0).sum())} of {len(main)}. Seen in the gallery (e.g. bumps of A
   inside grey shading after the green line).
9. **Negative A in the tail**: main-set tails whose A at t0 + 60 min is below -0.3 dB: {int((main.A_60min_db < -0.3).sum())}
   (baseline too high for that event, or a level change); they stay in the set (no shape rule).
10. **Status** is the sub-link's quality status from `baselines.csv` (good / degraded); `use_in_waa_main` already
   drops degraded, so `status_not_good` adds nothing new here.
"""
P.STEP2_REPORT.write_text(report, encoding="utf-8")
print("\n" + "=" * 100 + "\n" + report)
print(f"done in {time.time() - clock:.0f} s")
