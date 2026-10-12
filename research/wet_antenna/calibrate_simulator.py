"""wet_antenna step 3a (with the 3a-fix): calibrate the tail simulator (simulator.py) to tail set v1 and check it.

Purpose: before fitting real tails, know which drying laws our data can tell apart, given noise, the 0.3 dB RSL
quantization, per-tail baseline offsets, tail lengths and the event / link structure.

Calibration (the only uses of real data):
  - noise, per sub-link, on its real CONTROL tails (step 2b control, TRAIN, A_start < 0.5 dB): 10 s white sigma +
    1 min AR(1) slow wander (sigma_slow, tau_noise), chosen from a lookup grid so that the simulated dry series match
    the within-tail std of 1 min A, the share of minutes exactly on a 0.3 dB level (1 min tl_std = 0) and the
    autocorrelation of 1 min A at lags 1, 10, 60 min;
  - laws m1, m3, m4, m5, each with a log-normal tail-to-tail spread s of its time scale (half per link, half per tail):
    for every s on a grid, the base parameters are fitted by least squares to the real TRAIN v1 median of A / A_start
    over 0-6 h (noise-free truth); s is then chosen so that the OBSERVED simulated 10-90 % band matches the real band.

Skeleton: the TRAIN rows of tails_model_index.csv (same tails, sub-links, links, events, usable lengths, A_start).
Writes: simulated sets + index + noise parameters -> product; calibration table and the real-vs-sim
figure -> key; everything else -> diag. Report: outputs/reports/calibrate_simulator.md.
Run from the repository root (about 5 minutes):
    python -m research.wet_antenna.calibrate_simulator
"""
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import yaml
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cmlpipe.api import CMLData

sys.stdout.reconfigure(encoding="utf-8")
from research.wet_antenna import paths as P
from research.common.save import save
from research.wet_antenna import simulator as ws

CFG = yaml.safe_load(open(P.CONFIG_FILE, encoding="utf-8"))
SIM = CFG["wet_antenna"]["simulator"]
SW = SIM["switches"]
STEP, SPM, H = SIM["rsl_step_db"], SIM["samples_per_min"], SIM["horizon_min"]
LAGS = SIM["noise_acf_lags_min"]
LAWS = SIM["laws"]
SIM_FILES = {"m1": P.SIM_TAILS_M1, "m3": P.SIM_TAILS_M3, "m4": P.SIM_TAILS_M4, "m5": P.SIM_TAILS_M5}   # per law
INK, INK2 = "#1a1a1a", "#52514e"
# per source: line style, marker, colour (style and marker first: the reader may not separate colours)
SRC_STYLE = {"real": ("-", "o", INK), "m1": ("--", "s", "#2a78d6"), "m3": (":", "^", "#eb6834"),
             "m4": ("-.", "D", "#1baf7a"), "m5": ((0, (8, 2, 1, 2, 1, 2)), "v", "#7a7a7a")}
clock = time.time()


def tidy(ax):
    ax.grid(alpha=0.25)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)


def rel(p):
    return Path(p).resolve().relative_to(P.REPO).as_posix()


def md(df, index=True):
    df = df.reset_index() if index else df
    head = "| " + " | ".join(map(str, df.columns)) + " |\n|" + "---|" * len(df.columns)
    return head + "\n" + "\n".join("| " + " | ".join("" if pd.isna(x) else str(x) for x in row) + " |"
                                   for row in df.itertuples(index=False))


def grid_of(g):
    return list(np.round(np.arange(g[0], g[1] + g[2] / 2, g[2]), 4))


# =================================================================================================== skeleton + real
idx = pd.read_csv(P.TAILS_MODEL_INDEX)
skel = idx[idx.split_event == "train"].sort_values("tail_id").reset_index(drop=True)
real = pq.read_table(P.TAILS_MODEL, filters=[("split_event", "=", "train")],
                     columns=["tail_id", "minutes_since_t0_trace", "A_db"]).to_pandas()
real = real.merge(skel[["tail_id", "A_start_trace_db"]], on="tail_id")
real["A_norm"] = real.A_db / real.A_start_trace_db
nan_masks = {t: g.sort_values("minutes_since_t0_trace").A_db.isna().to_numpy() for t, g in real.groupby("tail_id")}
z = ws.spread_z(skel, SIM["seed"])
print(f"skeleton: {len(skel)} TRAIN v1 tails, {skel.sublink_id.nunique()} sub-links ({time.time() - clock:.0f} s)")


def curve_stats(df, col="A_norm", horizon=H):
    d = df[(df.minutes_since_t0_trace <= horizon) & np.isfinite(df[col])]
    g = d.groupby("minutes_since_t0_trace")[col]
    return pd.DataFrame({"median": g.median(), "p10": g.quantile(0.1), "p90": g.quantile(0.9), "n": g.size()})


real_c = curve_stats(real)
ok_t = real_c.index[real_c.n >= SIM["min_tails_per_minute"]].to_numpy()
target = real_c.loc[ok_t, "median"].to_numpy()

# =================================================================================================== 1. noise
print("noise: real control targets ...")
sel = pd.read_csv(P.TAILS_SELECTED, usecols=["tail_id", "sublink_id", "split_event", "in_control_2b", "t0", "t0_trace",
                                             "t_use_end_strict"], parse_dates=["t0", "t0_trace", "t_use_end_strict"])
ctl = sel[(sel.split_event == "train") & sel.in_control_2b & sel.sublink_id.isin(skel.sublink_id)].copy()
ctl["k0"] = ((ctl.t0_trace - ctl.t0) / pd.Timedelta("1min")).round().astype(int)
ctl["k1"] = ((ctl.t_use_end_strict - ctl.t0) / pd.Timedelta("1min")).round().astype(int)
cl = pq.read_table(P.TAILS_1MIN, columns=["tail_id", "minutes_since_t0", "A_db"],
                   filters=[("tail_id", "in", ctl.tail_id.tolist())]).to_pandas()
cl = cl.merge(ctl[["tail_id", "sublink_id", "t0", "k0", "k1"]], on="tail_id")
cl = cl[(cl.minutes_since_t0 >= cl.k0) & (cl.minutes_since_t0 < cl.k1)].sort_values(["tail_id", "minutes_since_t0"])
cl["t"] = cl.t0 + pd.to_timedelta(cl.minutes_since_t0, unit="min")
# on-level from the release: 1 min tl_std == 0 with 6 valid 10 s samples
d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
c1 = d.store("cml_1min")
T1 = pd.DatetimeIndex(c1.time.values)
subs = c1.sublink_id.values.astype(str)
need = set(cl.sublink_id)
on_parts = []
for a in range(0, len(subs), 32):                                  # one store chunk (32 sub-links) at a time
    block = [s for s in subs[a:a + 32] if s in need]
    if not block:
        continue
    ds = c1.isel(sublink_id=slice(a, a + 32))[["tl_std", "n_valid_tl"]].load()
    for s in block:
        j = list(subs[a:a + 32]).index(s)
        g = cl[(cl.sublink_id == s) & cl.A_db.notna()]
        ti = T1.searchsorted(g.t)
        sd = ds.tl_std.transpose("time", "sublink_id").values[ti, j]
        nv = ds.n_valid_tl.transpose("time", "sublink_id").values[ti, j]
        on_parts.append(pd.Series((sd == 0) & (nv == SPM), index=g.index))
cl["on"] = pd.concat(on_parts).reindex(cl.index)

real_rows = []
for s, gs in cl.groupby("sublink_id"):
    n = ss = nt = 0
    lag_sum = {k: 0.0 for k in LAGS}
    for _, g in gs.groupby("tail_id"):
        a_ = g.A_db.to_numpy(float)
        if np.isfinite(a_).sum() < 2:
            continue
        dd = a_ - np.nanmean(a_)
        n, ss, nt = n + np.isfinite(dd).sum(), ss + np.nansum(dd ** 2), nt + 1
        for k in LAGS:
            if len(dd) > k:
                lag_sum[k] += np.nansum(dd[:-k] * dd[k:])
    real_rows.append({"sublink_id": s, "n_control_tails": nt, "n_control_min": int(n),
                      "within_std_db": np.sqrt(ss / max(n - nt, 1)), "on_level_share": float(gs.on.mean()),
                      **{f"acf_{k}": lag_sum[k] / ss if ss > 0 else np.nan for k in LAGS},
                      "median_control_len_min": float(gs.groupby("tail_id").size().median())})
real_noise = pd.DataFrame(real_rows).set_index("sublink_id")
print(f"   {len(real_noise)} sub-links with control data ({time.time() - clock:.0f} s); lookup ...")

G = SIM["noise_grid"]
grid = {"sigma_white_db": grid_of(G["sigma_white_db"]), "sigma_slow_db": grid_of(G["sigma_slow_db"]),
        "tau_noise_min": G["tau_noise_min"]}
BLOCK = SIM["noise_lookup_block_min"]
if BLOCK == "auto":                                                 # same length as the real control tails
    BLOCK = int(round(real_noise.median_control_len_min.median() / 10) * 10)
look = ws.dry_lookup(grid, SIM["noise_lookup_blocks"], BLOCK, STEP, SPM, LAGS, SIM["seed"])
W_ = SIM["noise_weights"]
TARGETS = ["within_std_db", "on_level_share"] + [f"acf_{k}" for k in LAGS]


def pick_noise(r):
    """Grid point with the smallest summed squared error over the targets (each scaled by its weight)."""
    e = ((look.within_std_db - r.within_std_db) / W_["std_db"]) ** 2 + \
        ((look.on_level_share - r.on_level_share) / W_["on_level"]) ** 2
    for k in LAGS:
        if np.isfinite(getattr(r, f"acf_{k}")):
            e = e + ((look[f"acf_{k}"] - getattr(r, f"acf_{k}")) / W_["acf"]) ** 2
    return look.loc[e.idxmin()]


nz = skel[["sublink_id"]].drop_duplicates().merge(real_noise.reset_index(), on="sublink_id", how="left")
fit_ok = (nz.n_control_min >= SIM["noise_min_control_min"]).to_numpy()
picked = pd.DataFrame([pick_noise(r) for r in nz[fit_ok].itertuples()], index=nz.index[fit_ok])
med = picked[["sigma_white_db", "sigma_slow_db", "tau_noise_min"]].median()
for c in ["sigma_white_db", "sigma_slow_db", "tau_noise_min"]:
    nz[c] = picked[c].reindex(nz.index).fillna(med[c])
for c in TARGETS:
    nz[f"sim_{c}"] = picked[c].reindex(nz.index)
nz["source"] = np.where(fit_ok, "fit", "median_fallback")
save(nz.round(4), P.SIM_NOISE,
            index=False)
noise_by_sub = {r.sublink_id: (r.sigma_white_db, r.sigma_slow_db, r.tau_noise_min) for r in nz.itertuples()}
print(f"   noise fitted: {int(fit_ok.sum())}, fallback {int((~fit_ok).sum())} ({time.time() - clock:.0f} s)")

# =================================================================================================== 2. laws + spread
fixed = {"m3": {"n": SIM["m3_n"]}}
x0 = {"m1": [np.log(30.0)], "m3": [np.log(0.05)], "m4": [np.log(10.0), np.log(140.0), np.log(0.25 / 0.75)],
      "m5": [np.log(0.01), np.log(0.05)]}


def to_p(law, x):
    """Free parameters (log scale; phi through a logistic; tau_s = tau_f + exp(x)) -> parameter dict."""
    if law == "m1":
        return {"tau": np.exp(x[0])}
    if law == "m3":
        return {"c": np.exp(x[0]), **fixed["m3"]}
    if law == "m4":
        return {"tau_f": np.exp(x[0]), "tau_s": np.exp(x[0]) + np.exp(x[1]), "phi": 1 / (1 + np.exp(-x[2]))}
    return {"a": np.exp(x[0]), "b": np.exp(x[1])}


def model_median(law, p, f):
    with np.errstate(all="ignore"):
        return np.nanmedian(ws.truth_matrix(skel, law, p, f, H + 1), axis=0)[ok_t]


def fit(law, s):
    f = np.exp(s * z)
    res = least_squares(lambda x: model_median(law, to_p(law, x), f) - target, x0[law], method="trf")
    return to_p(law, res.x), float(np.sqrt(np.mean(res.fun ** 2)))


def observed(law, p, s, switches=SW):
    rows, info = ws.simulate(skel, law, p, noise_by_sub, switches, SIM["seed"], np.exp(s * z), STEP, SPM, nan_masks)
    rows = rows.merge(skel[["tail_id", "A_start_trace_db"]], on="tail_id")
    rows["A_norm"] = rows.A_db / rows.A_start_trace_db
    return rows, info


def band_rmse(rows):
    c = curve_stats(rows).reindex(ok_t)
    r = real_c.loc[ok_t]
    return float(np.sqrt(np.nanmean((c.p10 - r.p10) ** 2 + (c.p90 - r.p90) ** 2) / 2))


print("laws: median fit per s, then s from the band ...")
scan, params, spread = [], {}, {}
for law in LAWS:
    for s in SIM["spread_grid"]:
        p, r_med = fit(law, s)
        rows, _ = observed(law, p, s)
        scan.append({"law": law, "s": s, **{k: round(float(v), 5) for k, v in p.items()}, "rmse_median_truth": r_med,
                     "rmse_band_observed": band_rmse(rows)})
    sc = pd.DataFrame([x for x in scan if x["law"] == law])
    b = sc.loc[sc.rmse_band_observed.idxmin()]
    spread[law] = float(b.s)
    params[law] = {k: float(b[k]) for k in ws.PARAMS[law]}
    print(f"   {law}: s = {b.s} ({time.time() - clock:.0f} s)")
scan = pd.DataFrame(scan)

# =================================================================================================== final simulation
sims, infos, cal_rows = {}, [], []
for law in LAWS:
    rows, info = observed(law, params[law], spread[law])
    sims[law] = rows
    infos.append(info)
    c = curve_stats(rows).reindex(ok_t)
    sc = scan[scan.law == law]
    cal_rows.append({"law": law, "description": ws.LAWS[law],
                     "n_free_params": {"m1": 1, "m3": 1, "m4": 3, "m5": 2}[law] + 1,   # + the spread s
                     **{k: round(v, 5) for k, v in params[law].items()}, "s": spread[law],
                     "f_10_90": f"{np.exp(-1.2816 * spread[law]):.2f}-{np.exp(1.2816 * spread[law]):.2f}",
                     "s_at_grid_edge": spread[law] in (min(SIM["spread_grid"]), max(SIM["spread_grid"])),
                     "rmse_median_truth": round(float(sc[sc.s == spread[law]].rmse_median_truth.iloc[0]), 4),
                     "rmse_median_observed": round(float(np.sqrt(np.nanmean((c["median"] - target) ** 2))), 4),
                     "rmse_band_observed": round(band_rmse(rows), 4)})
    tbl = pa.Table.from_pandas(rows.drop(columns=["A_start_trace_db", "A_norm"]).assign(split_event="train"),
                               preserve_index=False)
    save(tbl, SIM_FILES[law], compression="zstd")
cal = pd.DataFrame(cal_rows)
save(cal, P.SIM_CALIBRATION, index=False)
info_all = pd.concat(infos).merge(skel[["tail_id", "sublink_id", "cml_id", "event_id", "usable_min",
                                        "A_start_trace_db"]], on="tail_id")
save(info_all.round(5), P.SIM_INDEX,
            index=False)
print(f"simulated ({time.time() - clock:.0f} s)")


# =================================================================================================== 3. checks
def draw(ax, c_real, c_sim, xmax, log, sc_):
    """Real: solid median + hatched 10-90 % band; sim: dashed median + dotted 10 / 90 % lines."""
    rr, ss = c_real[c_real.index <= xmax], c_sim[c_sim.index <= xmax]
    x_r, x_s = rr.index / sc_, ss.index / sc_
    lo, hi = rr.p10.to_numpy(), rr.p90.to_numpy()
    if log:
        lo = np.where(lo > 0.01, lo, 0.01)
    ax.fill_between(x_r, lo, hi, facecolor="none", edgecolor="#9a9a9a", hatch="///", lw=0, label="real 10-90 %")
    f_ = (lambda y: np.where(y > 0.01, y, np.nan)) if log else (lambda y: y)
    ax.plot(x_r, f_(rr["median"].to_numpy()), "-", color=INK, lw=2, label="real median")
    ax.plot(x_s, f_(ss["median"].to_numpy()), "--", color="#2a78d6", lw=2, label="sim median")
    for col in ("p10", "p90"):
        ax.plot(x_s, f_(ss[col].to_numpy()), ":", color="#2a78d6", lw=1.4, label="sim 10 / 90 %" if col == "p10" else None)
    if log:
        ax.set_yscale("log")
        ax.set_ylim(0.01, 2)
    else:
        ax.axhline(0, color=INK2, lw=0.6)
    tidy(ax)


sim_c = {law: curve_stats(sims[law]) for law in LAWS}
fig, ax = plt.subplots(len(LAWS), 4, figsize=(17, 3.4 * len(LAWS)), sharey="col")
for i, law in enumerate(LAWS):
    r_ = cal.set_index("law").loc[law]
    for j, (xmax, log) in enumerate(((90, False), (90, True), (H, False), (H, True))):
        sc_ = 1 if xmax <= 90 else 60
        draw(ax[i, j], real_c, sim_c[law], xmax, log, sc_)
        ax[i, j].set_xlabel("minutes since t0_trace" if sc_ == 1 else "hours since t0_trace", fontsize=8)
    ax[i, 0].set_ylabel(f"{law}: A / A_start", fontsize=9)
    ax[i, 0].set_title(f"{law} {ws.LAWS[law]}, spread s = {r_.s}: RMSE median {r_.rmse_median_observed:.3f}, "
                       f"band {r_.rmse_band_observed:.3f}", fontsize=9, loc="left")
ax[0, 3].legend(frameon=False, fontsize=7, loc="upper right")
fig.suptitle(f"Real TRAIN v1 ({len(skel)} tails) vs simulated tails per law (with spread and slow noise), A / A_start."
             "\nReal: solid median, hatched 10-90 % band; simulated: dashed median, dotted 10 / 90 %. Columns: 0-90 min "
             "linear, log; 0-6 h linear, log.", fontsize=10)
fig.tight_layout(rect=(0.02, 0, 1, 0.965))
save(fig, P.SIM_FIG_MEDIAN, dpi=100)
plt.close(fig)

skel_a0 = skel.set_index("tail_id").A_start_trace_db.to_dict()


def crossings(df):
    """Per tail: t_half, t_e (first minute the 5 min trailing mean of A < A_start / 2, / e; NaN = not reached),
    ratios at k minutes, std of A in the last hour of long tails."""
    out = []
    for tid, g in df.groupby("tail_id"):
        a_ = g.sort_values("minutes_since_t0_trace").A_db.to_numpy()
        a0 = skel_a0[tid]
        sm = pd.Series(a_).rolling(5, min_periods=1).mean().to_numpy()
        r = {"tail_id": tid, "A_start_trace_db": a0}
        for name, lev in (("t_half_min", a0 / 2), ("t_e_min", a0 / np.e)):
            hit = np.nonzero(sm < lev)[0]
            r[name] = int(hit[0]) if len(hit) else np.nan
        for k in SIM["shape_minutes"]:
            v = a_[k:k + 5]
            r[f"ratio_{k}min"] = float(np.nanmean(v)) / a0 if len(v) == 5 and np.isfinite(v).any() else np.nan
        r["last_hour_std_db"] = float(np.nanstd(a_[-60:])) if len(a_) >= SIM["long_tail_min"] else np.nan
        out.append(r)
    return pd.DataFrame(out)


cr = {"real": crossings(real)}
for law in LAWS:
    cr[law] = crossings(sims[law])


def summ(s):
    s = s.dropna()
    return {"n": len(s), "p10": round(float(s.quantile(0.1)), 3), "median": round(float(s.median()), 3),
            "p90": round(float(s.quantile(0.9)), 3)} if len(s) else {"n": 0}


chk_rows = []
for col in ["t_half_min", "t_e_min"] + [f"ratio_{k}min" for k in SIM["shape_minutes"]] + ["last_hour_std_db"]:
    for src, df in cr.items():
        chk_rows.append({"quantity": col, "source": src, **summ(df[col]),
                         "not_reached": int(df[col].isna().sum()) if col.startswith("t_") else np.nan})
checks = pd.DataFrame(chk_rows)
save(checks, P.SIM_CHECKS, index=False)
chk_wide = checks.pivot_table(index="quantity", columns="source", values="median", sort=False)[["real"] + LAWS]
chk_p90 = checks[checks.quantity == "last_hour_std_db"].set_index("source").p90

fig, ax = plt.subplots(2, len(LAWS), figsize=(16, 6.5), sharex=True)
bins = np.arange(0, 245, 5)
for i, col in enumerate(("t_half_min", "t_e_min")):
    for j, law in enumerate(LAWS):
        a_ = ax[i, j]
        a_.hist(cr["real"][col].dropna().clip(upper=240), bins=bins, facecolor="#dddddd", edgecolor="#777777",
                hatch="///", label="real")
        a_.hist(cr[law][col].dropna().clip(upper=240), bins=bins, histtype="step", color=INK, lw=1.6, ls="--",
                label=f"sim {law}")
        a_.set_title(f"{col}, {law}: median real {cr['real'][col].median():.0f} / sim {cr[law][col].median():.0f} min",
                     fontsize=8)
        a_.legend(frameon=False, fontsize=7)
        tidy(a_)
for a_ in ax[1]:
    a_.set_xlabel("minutes (clipped at 240)")
fig.suptitle("Time to A_start / 2 and A_start / e: real (hatched bars) vs simulated (dashed outline)", fontsize=10)
fig.tight_layout()
save(fig, P.SIM_FIG_CROSSING, dpi=100)
plt.close(fig)

# noise fit: 1:1 per target (each target its own marker)
fz = nz[nz.source == "fit"]
fig, ax = plt.subplots(1, len(TARGETS), figsize=(4 * len(TARGETS), 4))
marks = ["o", "^", "s", "D", "v"]
for a_, c, mk in zip(ax, TARGETS, marks):
    a_.scatter(fz[c], fz[f"sim_{c}"], s=14, marker=mk, facecolor="none", edgecolor=INK)
    lo_, hi_ = np.nanmin(np.r_[fz[c], fz[f"sim_{c}"]]), np.nanmax(np.r_[fz[c], fz[f"sim_{c}"]])
    a_.plot([lo_, hi_], [lo_, hi_], ls="--", color=INK2, lw=1)
    a_.text(hi_, hi_, " 1:1", fontsize=8, va="bottom")
    err = (fz[f"sim_{c}"] - fz[c]).abs().median()
    a_.set(xlabel=f"real control: {c}", ylabel=f"simulated: {c}", title=f"{c}: median |error| {err:.3f}")
    tidy(a_)
fig.suptitle(f"Noise model per sub-link ({len(fz)} fitted): real control tails vs simulated dry series "
             "(white 10 s + AR(1) 1 min + 0.3 dB quantization)", fontsize=10)
fig.tight_layout()
save(fig, P.SIM_FIG_NOISE, dpi=100)
plt.close(fig)

# gallery: 12 real tails with their m1 and m4 twins
pick = skel.sample(SIM["gallery_n"], random_state=SIM["gallery_seed"]).sort_values("tail_id")
fig, ax = plt.subplots(3, 4, figsize=(18, 10))
for a_, r in zip(ax.ravel(), pick.itertuples()):
    for src, df, ls, lw, col in (("real", real, "-", 1.0, INK), ("sim m1", sims["m1"], "--", 1.1, "#2a78d6"),
                                 ("sim m4", sims["m4"], ":", 1.4, "#1baf7a")):
        g = df[(df.tail_id == r.tail_id) & (df.minutes_since_t0_trace < 360)]
        a_.plot(g.minutes_since_t0_trace / 60, g.A_db, ls=ls, lw=lw, color=col, label=src)
    a_.axhline(0, color=INK2, lw=0.5)
    sw, ss, tn = noise_by_sub[r.sublink_id]
    a_.set_title(f"tail {r.tail_id} {r.sublink_id} {r.event_id}\nwhite {sw:.2f}, slow {ss:.2f} dB, tau_noise {tn:g} min",
                 fontsize=7)
    a_.tick_params(labelsize=7)
    tidy(a_)
ax[0, 0].legend(frameon=False, fontsize=7)
for a_ in ax[-1]:
    a_.set_xlabel("hours since t0_trace", fontsize=8)
for a_ in ax[:, 0]:
    a_.set_ylabel("A (dB)", fontsize=8)
fig.suptitle(f"{SIM['gallery_n']} random TRAIN v1 tails (seed {SIM['gallery_seed']}): real (solid) vs simulated m1 "
             "(dashed) and m4 (dotted): same tail_id, A_start, length, sub-link noise and NaN minutes", fontsize=10)
fig.tight_layout()
save(fig, P.SIM_FIG_GALLERY, dpi=90)
plt.close(fig)

# =================================================================================================== 4. t_half by A_start
edges = SIM["thalf_bins_db"] + [np.inf]
blabels = [f"{lo:g}-{hi:g}" if np.isfinite(hi) else f"> {lo:g}" for lo, hi in zip(edges[:-1], edges[1:])]
rng_b = np.random.default_rng(SIM["seed"])
ev_of = skel.set_index("tail_id").event_id
th_rows = []
for src, df in cr.items():
    df = df.assign(bin=pd.cut(df.A_start_trace_db, edges, labels=blabels, right=False), event_id=df.tail_id.map(ev_of))
    for lab in blabels:
        g = df[df.bin == lab]
        v = g.t_half_min.dropna()
        row = {"source": src, "A_start_bin_db": lab, "tails": len(g), "n_reached": len(v),
               "median_t_half_min": float(v.median()) if len(v) else np.nan}
        if src == "real" and len(v) >= 5:                          # 95 % CI: resample EVENTS
            evs = g.event_id.unique()
            reps = []
            for _ in range(1000):
                ix = rng_b.choice(evs, len(evs))
                vv = pd.concat([g[g.event_id == e].t_half_min for e in ix]).dropna()
                reps.append(vv.median() if len(vv) else np.nan)
            row["ci_low"], row["ci_high"] = np.nanpercentile(reps, [2.5, 97.5])
        th_rows.append(row)
thalf = pd.DataFrame(th_rows)
save(thalf.round(2), P.SIM_THALF, index=False)
th_wide = thalf.pivot(index="A_start_bin_db", columns="source", values="median_t_half_min").reindex(blabels)[
    ["real"] + LAWS]
th_n = thalf[thalf.source == "real"].set_index("A_start_bin_db").reindex(blabels)

fig, ax = plt.subplots(figsize=(8, 4.8))
xs = np.arange(len(blabels))
for src in ["real"] + LAWS:
    ls, mk, col = SRC_STYLE[src]
    y = th_wide[src].to_numpy()
    ax.plot(xs, y, ls=ls, marker=mk, color=col, lw=2.4 if src == "real" else 1.5, ms=7,
            mfc="white" if src != "real" else col, label=src if src == "real" else f"sim {src}")
r_ci = th_n[["ci_low", "ci_high"]].to_numpy()
ax.vlines(xs, r_ci[:, 0], r_ci[:, 1], color=INK, lw=1.2)
ax.set_xticks(xs, [f"{b}\n(n real = {int(n)})" for b, n in zip(blabels, th_n.tails)])
ax.set(xlabel="A_start_trace bin (dB)", ylabel="median t_half (min)",
       title="t_half vs A_start: flat for exponential / two-layer laws, falling for nonlinear laws", xlim=(-0.3, len(xs) - 0.4))
ax.legend(frameon=False, fontsize=8)
tidy(ax)
fig.text(0.01, 0.01, "Vertical bars: 95 % CI of the real median (1000 resamples of events). Tails that never reach "
         "A_start / 2 are left out.", fontsize=7)
fig.tight_layout(rect=(0, 0.03, 1, 1))
save(fig, P.SIM_FIG_THALF, dpi=110)
plt.close(fig)

# =================================================================================================== report
par_txt = "\n".join(f"- **{law}** ({ws.LAWS[law]}): " + ", ".join(f"{k} = {v:.4g}" for k, v in params[law].items())
                    + f"; spread s = {spread[law]} (time-scale factor 10-90 %: {cal.set_index('law').loc[law, 'f_10_90']})"
                    for law in LAWS)
nf = nz[nz.source == "fit"]
noise_tab = pd.DataFrame({c: {"real median": round(float(nf[c].median()), 3),
                              "sim median": round(float(nf[f"sim_{c}"].median()), 3),
                              "median abs error": round(float((nf[f"sim_{c}"] - nf[c]).abs().median()), 3)}
                          for c in TARGETS}).T
par_tab = pd.DataFrame({c: {"p10": round(float(nf[c].quantile(.1)), 3), "median": round(float(nf[c].median()), 3),
                            "p90": round(float(nf[c].quantile(.9)), 3)}
                        for c in ["sigma_white_db", "sigma_slow_db", "tau_noise_min"]}).T
edge_tau = int((nf.tau_noise_min == max(grid["tau_noise_min"])).sum())
edge_ss = int((nf.sigma_slow_db == max(grid["sigma_slow_db"])).sum())
scan_show = scan.pivot(index="s", columns="law", values="rmse_band_observed").round(4)
report = f"""# wet_antenna step 3a-fix: slow noise and a spread for every law (TRAIN v1 skeleton)

Written by `research/wet_antenna/calibrate_simulator.py` (library `research/wet_antenna/simulator.py`, settings `wet_antenna.simulator`); this
replaces the step 3a calibration (m2 dropped: it is m1 with spread). Skeleton = the {len(skel)} TRAIN v1 tails.
Switches: {SW}. Run time {time.time() - clock:.0f} s.

## 1. Slow noise (per sub-link, real control tails)

Noise = 10 s white (sigma_white) + 1 min AR(1) wander (sigma_slow, correlation time tau_noise) + 0.3 dB quantization.
Chosen per sub-link from a lookup grid ({len(look)} points; dry series of {BLOCK} min = the median control tail length, {SIM['noise_lookup_blocks']} per point)
by the summed squared error over 5 targets (scales: std {W_['std_db']} dB, on-level {W_['on_level']}, ACF {W_['acf']}).
Fitted: {len(nf)} sub-links; median fallback: {int((nz.source != 'fit').sum())}.

{md(noise_tab)}

Fitted parameters (p10 / median / p90):

{md(par_tab)}

At the grid edge: tau_noise = {max(grid['tau_noise_min'])} min for {edge_tau} sub-links, sigma_slow = {max(grid['sigma_slow_db'])} dB for {edge_ss}.
Before the fix (step 3a, white noise only) the within-tail std was 0.051 dB simulated vs 0.141 dB real.

![noise](../figures/{P.SIM_FIG_NOISE.name})

## 2. Laws with spread

For every s on the grid, the base parameters are fitted to the real median of A / A_start over 0-{H} min (noise-free
truth with the spread, same survivors); s is the one whose OBSERVED 10-90 % band (full noise model) is closest to the
real band. Band misfit per s (units of A / A_start):

{md(scan_show)}

{par_txt}

## Calibration and misfit (`{rel(P.SIM_CALIBRATION)}`, key)

{md(cal.drop(columns='description'), index=False)}

![real vs sim](../figures/{P.SIM_FIG_MEDIAN.name})

## 3. Checks (medians; p10 / p90 in `{rel(P.SIM_CHECKS)}`)

{md(chk_wide.round(3))}

Last-hour std p90: {', '.join(f"{k} {v:.3f}" for k, v in chk_p90.items())} dB. Not reached (t_half / t_e):
{', '.join(f"{s}: {int(cr[s].t_half_min.isna().sum())} / {int(cr[s].t_e_min.isna().sum())}" for s in cr)}.

![t_half t_e](../figures/{P.SIM_FIG_CROSSING.name})

![gallery](../figures/{P.SIM_FIG_GALLERY.name})

## 4. t_half per A_start bin (`{rel(P.SIM_THALF)}`)

Exponential (m1) and two-layer (m4) laws give the same t_half for every A_start; nonlinear laws (m3, m5) give shorter
t_half for larger A_start. Median t_half (min); real with a 95 % CI from resampling events:

{md(th_wide.round(1))}

Real: {', '.join(f"{b}: n = {int(r.tails)}, CI [{r.ci_low:.0f}, {r.ci_high:.0f}]" if np.isfinite(r.get('ci_low', np.nan)) else f"{b}: n = {int(r.tails)}" for b, r in th_n.iterrows())}.

![t_half by A_start](../figures/{P.SIM_FIG_THALF.name})

## Unclear

1. **The lookup uses one block length** ({BLOCK} min, the median control tail length; per sub-link the median ranges
   {real_noise.median_control_len_min.min():.0f}-{real_noise.median_control_len_min.max():.0f} min). Within-tail std and ACF depend on the
   length through the per-tail mean, so the fit is approximate for sub-links with short control tails.
2. **The noise targets are from the step 2b control set** (not rebuilt with the v1 radar-NaN rule) and may contain a
   little late drying or light rain; this would push sigma_slow up.
3. **One spread for all rate parameters** (m4: tau_f and tau_s scaled together, phi fixed; m5: a and b together). Real
   tails may vary in phi or in the a / b balance; the band fit cannot see that.
4. **The spread is calibrated on the band, the base law on the median**, as two stages; they are not fitted jointly,
   so the median fit at the chosen s may be a little off from the best joint fit.
5. **t_half is censored**: tails that end before reaching A_start / 2 are left out of the median (counts in the table);
   censoring is stronger for slow tails, so all medians are biased low by about the same rule in real and simulated
   data (the simulator copies the real lengths).
6. **The quantization phase** is drawn uniformly per tail and TSL steps are not modelled.
7. **t_half per A_start bin cannot separate the laws yet**: the real 95 % CIs (resampling events) overlap every law in
   every bin, and the > 3 dB bin has {int(th_n.tails.iloc[-1])} tails. Even m1, whose truth t_half does not depend on A_start,
   drops in that bin ({th_wide.loc[blabels[-1], 'm1']:.0f} min): with so few tails the median follows the spread values these
   tails happened to draw. Step 3b should repeat the simulation with many seeds before reading this check.
8. **One AR(1) time scale per sub-link**: for some sub-links lags 1, 10 and 60 cannot all be matched (the horizontal
   stripes in the acf_60 panel: grid tau values); real dry noise likely has more than one time scale (minutes and the
   daily cycle).
"""
save(report, P.SIM_REPORT)
print("\n" + "=" * 100 + "\n" + report)
print(f"done in {time.time() - clock:.0f} s")
