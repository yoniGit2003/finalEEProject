"""wet_antenna step 4c (Phase A): what is in a drying tail? The signal before the law.

Steps 4 and 4b fitted drying laws to tails that wait, rise or stall, and found no law. Before any new fit this step
asks what the signal in a tail IS. Three questions (docs/ROADMAP.md, plan change of 2026-10-12), NO law is fitted:
  A1 (Q1)  Is the rise of a tail water at all? Each tail is compared with the SAME sub-link at the SAME clock time on
           certified dry days. A rise that dry days also show is drift (daily cycle); one they never show is water
           (dew, drizzle under the radar, or rain the references missed).
  A2 (Q2)  Is the early drop the antenna or the path? Path rain grows with the path length L, the antenna does not.
           Per minute m from -30 to +60 around t0_trace: A(m) = alpha_event(m) + beta(m) L + gamma(m) [f >= 25 GHz].
           beta(m) = dB per km = the path part; alpha(m) = what all links of an event share = the antenna part.
  A3 (Q3)  When does drying start? Three start rules (S0 = t0_trace of v1, S1 = strict dry, S2 = signal onset) and
           what each does to the counts, the classes, the shape and the model-free k(A). It ends with a PROPOSAL for
           tail set v2; nothing is frozen here.
TRAIN tails of the frozen set v1 only (VAL rows are dropped when a file is read and never used). Every rule is in
config `wet_antenna.characterise_tails`, written before the run. A1 is timed on 20 tails first.
Reads tails_selected.csv and tails_1min.parquet (steps 1, 2, 2b), the frozen set, the 4b classes, the baseline
products and the release. Results of A1 and of the start rules are cached in
research/.cache/wet_antenna/characterise_tails/ (delete the folder to compute again).
Run from the repository root (about 2 min):
    python -m research.wet_antenna.characterise_tails
"""
import json
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cmlpipe.api import CMLData
from research.common.config import load_config
from research.common.plotting import INK, INK2, style_axes
from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import datasets as dsets
from research.wet_antenna import free_curve as wc
from research.wet_antenna import paths as P
from research.wet_antenna.tails import classify, first_run, itu_p838, minutes, raw_rate, reference_traces, tail_from_start

CFG = load_config()
W = CFG["wet_antenna"]
S = W["characterise_tails"]
A1, A2, A3 = S["a1"], S["a2"], S["a3"]
S2SEL, B2B, FRZ = W["selection"], W["refine_selection"], W["freeze"]       # the rules of tail set v1 (unchanged)
CL, RR = W["diagnose_stalled_tails"]["classes"], W["fit_real_drying"]["raw_rate"]
W_LO = W["confirm_free_curve"]["w_lo_db"]
WIN, BEFORE = S["window_min"], S["before_min"]
CACHE = P.CACHE / "characterise_tails"
CLASSES = ["decays", "plateau-then-decays", "rises", "flat / slow"]
RULES = {"S0": "S0 = t0_trace (v1)", "S1": "S1 = strict dry", "S2": "S2 = signal onset"}
LINE = {"S0": "-", "S1": "--", "S2": ":"}                  # one line style per start rule, in every figure
MIN = pd.Timedelta("1min")


def rise_of(a):
    """The rise of 60 minutes of a signal, exactly as step 4b measures it for a tail: the maximum of the centred
    5 min median, minus the mean of the first 5 minutes. Also the change from the first to the last 5 minutes."""
    sm = pd.Series(a).rolling(CL["smooth_min"], center=True, min_periods=CL["smooth_min"] // 2 + 1).median().to_numpy()
    start = np.nanmean(a[:CL["start_minutes"]])
    return float(np.nanmax(sm) - start), float(np.nanmean(a[-CL["start_minutes"]:]) - start)


def event_boot(df, stat, rng, reps, interval):
    """Event bootstrap: draw whole events with replacement, compute `stat` on the tails of the draw (as step 4b)."""
    groups = [g for _, g in df.groupby("event_id")]
    vals = np.asarray([stat(pd.concat([groups[i].assign(draw=k) for k, i in enumerate(rng.integers(0, len(groups), len(groups)))]))
                       for _ in range(reps)], float)
    return (np.nan, np.nan) if np.isnan(vals).mean() > 0.5 else tuple(np.nanpercentile(vals, interval))


def share(x):
    """Share of True in a column that may hold True / False / NaN (NaN = unknown, left out). NaN if nothing is known."""
    x = pd.Series(x).dropna()
    return float(x.astype(bool).mean()) if len(x) else np.nan


# =================================================================================================== A1
def a1_controls(tails, TL, col, i_of, state, st_col, in_jump, n_days):
    """Per tail: its rise, and the rises of the same sub-link in the same clock window on certified dry days."""
    n5 = WIN // 5
    out = []
    for r in tails.itertuples():
        tl, i_tail = TL[:, col[r.sublink_id]], i_of(r.t0_trace)
        st, since = state["state"][:, st_col[r.cml_id]], state["since_wet_h"][:, st_col[r.cml_id]]
        rise, end_change = rise_of(tl[i_tail:i_tail + WIN])
        ctrl, ends = [], []
        for day in range(n_days):
            i0 = day * 1440 + i_tail % 1440
            if i0 + WIN > len(tl) or abs(i0 - i_tail) < 1440:
                continue
            b0 = i0 // 5
            if not (st[b0:b0 + n5] == 0).all() or since[b0] < A1["min_hours_since_wet"]:
                continue                                   # not dry in the whole window, or too soon after rain
            if in_jump[r.sublink_id][i0:i0 + WIN].any():
                continue                                   # inside a dry-weather jump of this sub-link (dry_jumps.csv)
            w = tl[i0:i0 + WIN]
            if np.isfinite(w).mean() < A1["min_valid_share"] or not np.isfinite(w[:CL["start_minutes"]]).any():
                continue
            c_rise, c_end = rise_of(w)
            ctrl.append(c_rise)
            ends.append(c_end)
        has = len(ctrl) >= A1["min_control_days"]
        p = float(np.percentile(ctrl, A1["control_percentile"])) if has else np.nan
        out.append({"tail_id": r.tail_id, "rise_db": round(rise, 3), "end_change_db": round(end_change, 3),
                    "control_days": len(ctrl), "has_control": has,
                    "control_p95_db": round(p, 3) if has else np.nan,
                    "control_median_db": round(float(np.median(ctrl)), 3) if has else np.nan,
                    "control_max_db": round(float(np.max(ctrl)), 3) if has else np.nan,
                    "control_end_change_median_db": round(float(np.median(ends)), 3) if has else np.nan,
                    "control_end_change_p5_db": round(float(np.percentile(ends, 5)), 3) if has else np.nan,
                    "control_end_change_p95_db": round(float(np.percentile(ends, 95)), 3) if has else np.nan,
                    "share_of_control_days_below": round(float(np.mean(np.asarray(ctrl) < rise)), 3) if has else np.nan,
                    "exceeds_control": bool(rise > p) if has else np.nan})
    return pd.DataFrame(out)


def a1_row(group, g):
    c = g[g.has_control]
    return {"group": group, "tails": len(g), "with_control": len(c), "exceed": int(c.exceeds_control.astype(bool).sum()),
            "share_exceed": round(float(c.exceeds_control.astype(bool).mean()), 3) if len(c) else np.nan,
            "median_rise_db": round(float(g.rise_db.median()), 3), "median_control_p95_db": round(float(c.control_p95_db.median()), 3) if len(c) else np.nan,
            "median_control_days": float(c.control_days.median()) if len(c) else np.nan,
            "median_control_end_change_db": round(float(c.control_end_change_median_db.median()), 3) if len(c) else np.nan}


# =================================================================================================== A2
def length_regression(Y, L, high, events, rng):
    """Per minute (column of Y): least squares of A on [L, high band] with one intercept per event (event fixed
    effects), and an event bootstrap of the coefficients. Works on sums per event, so a bootstrap draw is only a
    weighted sum. Returns per minute: beta, gamma, alpha (tail-weighted mean of the event intercepts) and the draws."""
    X = np.column_stack([L, high]).astype(float)
    ev = np.unique(events)
    M = Y.shape[1]
    Sxx, Sxy, Sx, Sy, N = np.zeros((len(ev), M, 2, 2)), np.zeros((len(ev), M, 2)), np.zeros((len(ev), M, 2)), np.zeros((len(ev), M)), np.zeros((len(ev), M))
    for k, e in enumerate(ev):
        rows = events == e
        ok, x = np.isfinite(Y[rows]).astype(float), X[rows]
        y = np.where(ok > 0, Y[rows], 0.0)
        n = ok.sum(axis=0)
        sx, sy = ok.T @ x, y.sum(axis=0)
        nn = np.maximum(n, 1)
        Sxx[k] = np.einsum("im,ij,ik->mjk", ok, x, x) - sx[:, :, None] * sx[:, None, :] / nn[:, None, None]   # centred inside the event
        Sxy[k] = np.einsum("im,ij->mj", y, x) - sx * (sy / nn)[:, None]
        Sx[k], Sy[k], N[k] = sx, sy, n

    def solve(w):
        coef = np.einsum("mjk,mk->mj", np.linalg.pinv(np.tensordot(w, Sxx, 1), hermitian=True), np.tensordot(w, Sxy, 1))
        alpha = (w @ Sy - np.einsum("mj,mj->m", coef, np.tensordot(w, Sx, 1))) / np.maximum(w @ N, 1)
        return coef[:, 0], coef[:, 1], alpha

    beta, gamma, alpha = solve(np.ones(len(ev)))
    draws = np.array([solve(w.astype(float)) for w in rng.multinomial(len(ev), np.full(len(ev), 1 / len(ev)), A2["boot_reps"])])
    return beta, gamma, alpha, draws, N.sum(axis=0)


def a2_set(name, ix, TL, col, i_of, base, R5, r_col, rng):
    """Part A2 on one tail set: the table per minute, the power check and N."""
    m = np.arange(-BEFORE, WIN + 1)
    Y = np.array([TL[i_of(r.t0_trace) - BEFORE:i_of(r.t0_trace) + WIN + 1, col[r.sublink_id]] - base[r.tail_id] for r in ix.itertuples()], float)
    L, high = ix.L_km.to_numpy(float), (ix.frequency_ghz >= A2["band_split_ghz"]).to_numpy(float)
    beta, gamma, alpha, draws, n = length_regression(Y, L, high, ix.event_id.to_numpy(), rng)
    lo, hi = np.percentile(draws[:, 0], A2["interval"], axis=0)
    a_lo, a_hi = np.percentile(draws[:, 2] + draws[:, 1], A2["interval"], axis=0)
    # the beta expected from path rain: mean over the tails of a R_mean^b (dB/km), ITU-R P.838 at the frequency and
    # polarization of the sub-link, R_mean = the radar along the path in the 5 min bin of the minute
    k_exp = np.full(Y.shape, np.nan)
    for j, r in enumerate(ix.itertuples()):
        a_, b_ = itu_p838(r.frequency_ghz, r.polarization)
        rr = R5[(i_of(r.t0_trace) + m) // 5, r_col[r.cml_id]]
        k_exp[j] = a_ * np.power(rr, b_)
    with np.errstate(all="ignore"):
        tab = pd.DataFrame({"set": name, "minute": m, "tails": n.astype(int), "mean_A_db": np.nanmean(Y, axis=0),
                            "beta_db_per_km": beta, "beta_lo": lo, "beta_hi": hi,
                            "expected_path_rain_db_per_km": np.nanmean(k_exp, axis=0),
                            "expected_path_rain_db": np.nanmean(k_exp * L[:, None], axis=0),
                            "gamma_high_band_db": gamma, "alpha_low_band_db": alpha, "alpha_high_band_db": alpha + gamma,
                            "alpha_high_lo": a_lo, "alpha_high_hi": a_hi,
                            "path_part_median_link_db": beta * float(np.median(L))}).round(4)
    pm = tab[tab.minute.between(*A2["power_minutes"])]
    power = bool((pm.beta_lo > 0).any())
    zero_in = ((tab.beta_lo <= 0) & (tab.beta_hi >= 0) & (tab.minute >= 0)).to_numpy()
    k = first_run(zero_in, A2["run_min"])
    n_min = int(tab.minute.to_numpy()[k]) if k is not None else None
    pre = tab[tab.minute < 0]
    info = {"set": name, "tails": len(ix), "events": int(ix.event_id.nunique()), "links": int(ix.cml_id.nunique()),
            "median_L_km": round(float(np.median(L)), 2), "L_p10_km": round(float(np.percentile(L, 10)), 2),
            "L_p90_km": round(float(np.percentile(L, 90)), 2), "share_high_band": round(float(high.mean()), 3),
            "power": power, "minutes_with_beta_above_0_in_power_range": int((pm.beta_lo > 0).sum()),
            "max_beta_in_power_range": round(float(pm.beta_db_per_km.max()), 3), "N": n_min,
            "minutes_ge0_with_beta_above_0": int(((tab.beta_lo > 0) & (tab.minute >= 0)).sum()),
            "minutes_ge0_with_beta_below_0": int(((tab.beta_hi < 0) & (tab.minute >= 0)).sum()),
            "expected_inside_interval_before_0": round(float(((pre.expected_path_rain_db_per_km >= pre.beta_lo) & (pre.expected_path_rain_db_per_km <= pre.beta_hi)).mean()), 3),
            "beta_over_expected_before_0": round(float(pre.beta_db_per_km.mean() / max(pre.expected_path_rain_db_per_km.mean(), 1e-9)), 2)}
    return tab, info


# =================================================================================================== A3
def start_rules(pool, long_a, bounds, T5, trace5, nan5, broad5):
    """The three starts and, for each, the end by the v1 end rules, A_start and the usable length. One row per pair.
    The code of the end rules follows freeze_tail_set.py (step 2c); only the start minute differs."""
    rows = []
    for i, r in enumerate(pool.itertuples()):
        a = long_a[bounds[i, 0]:bounds[i, 1]]
        L = len(a)
        i05 = T5.searchsorted(r.t0)
        fd = int(r.margin_minutes) if pd.notna(r.margin_minutes) else L
        tr, nn = trace5[r.cml_id], nan5[r.cml_id]
        unk = tr | nn if FRZ["radar_nan_is_unknown"] else tr
        n_margin = int(np.ceil(min(fd, L) / 5))

        def last_trace_end(mask):
            mb = np.nonzero(mask[i05:i05 + n_margin])[0]
            return min(int((mb[-1] + 1) * 5) if len(mb) else 0, L)

        k_len = min(minutes(r.t_use_end_lenient, r.t0), L) if pd.notna(r.t_use_end_lenient) else L
        sm10 = pd.Series(a).rolling(FRZ["rise_smooth_min"], min_periods=FRZ["rise_smooth_min"] * 2 // 3).median().to_numpy()
        sm5 = pd.Series(a).rolling(CL["smooth_min"], center=True, min_periods=CL["smooth_min"] // 2 + 1).median().to_numpy()

        def tail_from(k):
            """End, usable length and A_start of the tail that starts at minute k (v1 end rules; code in tails.py)."""
            return tail_from_start(a, sm10, k, unk, i05, k_len, r.truncated_by_lenient, FRZ, S2SEL["a_dry_window_min"])

        s0 = tail_from(last_trace_end(unk))
        s1 = tail_from(min(last_trace_end(unk | broad5[r.cml_id]) + A3["s1_buffer_min"], L))
        hi = min(s0["k"] + A3["s2_search_min"] + 1, s0["k_trace"], k_len, L)   # S2: never beyond a rain trace or the 2b end
        seg = sm5[s0["k"]:hi]
        k2 = s0["k"] + int(np.nanargmax(seg)) if len(seg) and np.isfinite(seg).any() else s0["k"]
        s2 = tail_from(k2)
        peak = float(np.nanmax(seg) - s0["A_start_db"]) if len(seg) and np.isfinite(seg).any() and np.isfinite(s0["A_start_db"]) else np.nan
        # added after the first run (information): the same peak above the 5 min median AT t0_trace. A_start is a mean
        # over 5 min, so a tail that falls fast has its first minute above A_start without any rise.
        peak0 = float(np.nanmax(seg) - seg[0]) if len(seg) and np.isfinite(seg[0]) else np.nan
        row = {"tail_id": r.tail_id, "s2_peak_above_A_start_db": round(peak, 3) if np.isfinite(peak) else np.nan,
               "s2_rise_before_onset": bool(peak >= A3["s2_rise_flag_db"]) if np.isfinite(peak) else False,
               "s2_peak_above_level_at_t0_trace_db": round(peak0, 3) if np.isfinite(peak0) else np.nan}
        for name, s in (("S0", s0), ("S1", s1), ("S2", s2)):
            row.update({f"{name}_k": s["k"], f"{name}_end_k": s["end"], f"{name}_usable_min": s["usable_min"],
                        f"{name}_truncated_by": s["truncated_by"],
                        f"{name}_A_start_db": round(s["A_start_db"], 4) if np.isfinite(s["A_start_db"]) else np.nan,
                        f"{name}_kept": bool(np.isfinite(s["A_start_db"]) and s["A_start_db"] >= S2SEL["main_a_start_db"]
                                             and s["usable_min"] >= S2SEL["min_usable_min"])})
        rows.append(row)
    return pd.DataFrame(rows)


def rate_slope(dd):
    """The model-free k(A) of step 4b (B5) and the slope of log k vs log A from 0.35 dB up."""
    if "draw" in dd:                                       # event bootstrap: a copy of an event is a new set of tails
        dd = dd.assign(tail_id=dd.tail_id * 100 + dd.draw)
    r_ = raw_rate(dd, RR["edges_db"], RR["smooth_min"])
    r_["used_for_slope"] = (r_.tail_minutes >= A3["rate_min_minutes"]) & (r_.k_from_mean_per_min > 0)
    u = r_[r_.used_for_slope]
    return r_, (float(np.log(u.k_from_mean_per_min) @ wc.slope_weights(np.log(u.median_A_db))) if len(u) >= 3 else np.nan)


def describe_rule(name, st, pool, long_a, bounds, night_of, rng, with_boot=True):
    """Everything the report says about one start rule, computed on the tails the rule keeps (first 60 min)."""
    keep = st[st[f"{name}_kept"]].merge(pool[["tail_id", "event_id", "cml_id", "sublink_id", "t0"]].assign(pos=np.arange(len(pool))), on="tail_id")
    parts, cls = [], []
    for r in keep.itertuples():
        k = getattr(r, f"{name}_k")
        a = long_a[bounds[r.pos, 0] + k:bounds[r.pos, 0] + min(k + WIN + S2SEL["value_window_min"], getattr(r, f"{name}_end_k"))]
        c, _ = classify(a[:WIN], CL)                       # the class from the first 60 min, as in step 4b
        cls.append({"tail_id": r.tail_id, "class": c["class"], "min_to_half": c["min_to_half"]})
        parts.append(pd.DataFrame({"tail_id": r.tail_id, "event_id": r.event_id, "minutes_since_t0_trace": np.arange(len(a)), "A_db": a}))
    rows, cls = pd.concat(parts, ignore_index=True), pd.DataFrame(cls)
    a0 = keep.set_index("tail_id")[f"{name}_A_start_db"]
    shape = dsets.shape(rows, a0, A3["shape_minutes"])     # A at k min = mean over [k, k + 5), so 60 needs minutes 60-64
    rows = rows[rows.minutes_since_t0_trace < WIN]
    rate, slope = rate_slope(rows)
    lo, hi = event_boot(rows, lambda q: rate_slope(q)[1], rng, A3["rate_boot_reps"], A2["interval"]) if with_boot else (np.nan, np.nan)
    delay = (keep[f"{name}_k"] - keep.S0_k).astype(float)
    start = keep.t0 + pd.to_timedelta(keep[f"{name}_k"], unit="min")
    night = start.dt.floor("1h").map(night_of)
    th = cls.min_to_half.fillna(np.inf)
    row = {"start_rule": RULES[name], "tails": len(keep), "events": int(keep.event_id.nunique()), "links": int(keep.cml_id.nunique()),
           **{f"share_{c}": round(float((cls["class"] == c).mean()), 3) for c in CLASSES},
           **{f"A_over_A_start_{k}min": round(v, 3) for k, v in shape.items()},
           "t_half_min": float(np.median(th)) if np.isfinite(np.median(th)) else np.nan,
           "share_half_within_60min": round(float(np.isfinite(th).mean()), 3),
           "k_slope": round(slope, 2), "k_slope_lo": round(lo, 2), "k_slope_hi": round(hi, 2), "bins_in_slope": int(rate.used_for_slope.sum()),
           "delay_median_min": float(delay.median()), "delay_q25_min": float(delay.quantile(0.25)), "delay_q75_min": float(delay.quantile(0.75)),
           "delay_median_night_min": float(delay[night == True].median()) if (night == True).any() else np.nan,     # noqa: E712
           "delay_median_day_min": float(delay[night == False].median()) if (night == False).any() else np.nan,     # noqa: E712
           "tails_night": int((night == True).sum()), "tails_day": int((night == False).sum())}                      # noqa: E712
    curve = rows.assign(ratio=rows.A_db / rows.tail_id.map(a0)).groupby("minutes_since_t0_trace").agg(
        median_A_over_A_start=("ratio", "median"), median_A_db=("A_db", "median"), tails=("A_db", "count")).reset_index()
    return row, rate.assign(start_rule=RULES[name], slope=round(slope, 3)), curve.assign(start_rule=RULES[name]).round(4), set(keep.tail_id), cls


# =================================================================================================== figures
def figure_a1(t):
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.6), sharex=True, sharey=True)
    top = 1.6
    for ax, c in zip(axes, CLASSES):
        g = t[(t["class"] == c) & t.has_control]
        for nightval, mk, fill, lab in ((True, "^", INK, "night (filled triangles)"), (False, "o", "none", "day (open circles)")):
            q = g[g.night == nightval]
            ax.plot(q.control_p95_db.clip(upper=top), q.rise_db.clip(upper=top), mk, ms=5, mfc=fill, mec=INK, mew=0.8, ls="none", label=lab)
        ax.plot([0, top], [0, top], "--", color=INK2, lw=1)
        ax.text(top * 0.97, top * 0.9, "rise = dry-day\n95th percentile\n(dashed)", fontsize=7.5, ha="right", va="top")
        n_ex = int(g.exceeds_control.astype(bool).sum())
        ax.set_title(f"{c.upper()}\n{n_ex} of {len(g)} above their dry days ({n_ex / max(len(g), 1):.0%})", fontsize=9)
        ax.set_xlabel("95th percentile of the rise on the\nsub-link's own dry days, same clock time (dB)")
        style_axes(ax)
    axes[0].set_ylabel("rise of the tail in its first 60 min (dB)")
    axes[0].legend(frameon=False, fontsize=8, loc="center right")
    fig.suptitle("Step 4c, A1: the rise of each tail (max of the 5 min median - A_start) against what the same sub-link does at the "
                 f"same clock time in dry weather. Above the dashed line = more than drift. Values above {top:g} dB are drawn at {top:g}.", fontsize=9.5)
    fig.tight_layout()
    return fig


def figure_a2(tabs, infos):
    fig, axes = plt.subplots(len(tabs), 2, figsize=(14, 4.6 * len(tabs)), squeeze=False)
    for row, (name, tab) in enumerate(tabs.items()):
        info, ax, bx = infos[name], axes[row, 0], axes[row, 1]
        ax.fill_between(tab.minute, tab.beta_lo, tab.beta_hi, facecolor="none", edgecolor="#8a8a8a", hatch="///", lw=0)
        ax.plot(tab.minute, tab.beta_db_per_km, "-", color=INK, lw=1.8)
        ax.plot(tab.minute, tab.expected_path_rain_db_per_km, "--", color=INK, lw=1.4)
        ax.axhline(0, color=INK2, lw=0.8)
        ax.text(0.98, 0.97, "solid: beta(m)\nhatched: its 90 % interval\ndashed: beta expected from\nthe radar on the path", transform=ax.transAxes,
                fontsize=8, ha="right", va="top", bbox=dict(facecolor="white", edgecolor="none", pad=1.5, alpha=0.85))
        ax.set(xlabel="minutes since t0_trace", ylabel="beta = dB per km of path")
        ax.set_title(f"{name}\n{info['tails']} tails, {info['events']} events. Power check: {'passed' if info['power'] else 'NOT passed'}; "
                     f"N = {info['N'] if info['N'] is not None else 'not reached'}", fontsize=9.5)
        bx.plot(tab.minute, tab.mean_A_db, "-", color=INK, lw=1.8)
        bx.plot(tab.minute, tab.alpha_high_band_db, "--", color=INK, lw=1.5)
        bx.plot(tab.minute, tab.alpha_low_band_db, ":", color=INK, lw=1.5)
        bx.plot(tab.minute, tab.path_part_median_link_db, "-.", color=INK, lw=1.3)
        bx.text(0.98, 0.97, "solid: mean A of the tails\ndashed: alpha for >= 25 GHz = A at L = 0\ndotted: alpha for < 25 GHz\n"
                "dash-dot: beta x median L = path part of a median link", transform=bx.transAxes, fontsize=8, ha="right", va="top",
                bbox=dict(facecolor="white", edgecolor="none", pad=1.5, alpha=0.85))
        bx.axhline(0, color=INK2, lw=0.8)
        bx.set(xlabel="minutes since t0_trace", ylabel="dB")
        bx.set_title(f"{name}\nthe part shared by the links of an event (alpha) and the path part", fontsize=9.5)
        for a_ in (ax, bx):
            a_.axvline(0, color=INK, lw=0.8)
            a_.set_xlim(-BEFORE, WIN + 1)
            style_axes(a_)
        if info["N"] is not None:
            ax.axvline(info["N"], color=INK, lw=1, ls=":")
            ax.text(info["N"] + 0.5, ax.get_ylim()[1], f" N = {info['N']}", fontsize=8, va="top")
    fig.suptitle("Step 4c, A2: does A grow with the path length, minute by minute? A(m) = alpha_event(m) + beta(m) L + gamma(m) [f >= 25 GHz]. "
                 "Vertical line: t0_trace.", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig


def figure_a3(curves):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    for name, key in RULES.items():
        c = curves[curves.start_rule == key]
        for ax, colname in zip(axes, ("median_A_over_A_start", "median_A_db")):
            ax.plot(c.minutes_since_t0_trace, c[colname], LINE[name], color=INK, lw=1.9, label=f"{key} ({int(c.tails.max())} tails)")
            ax.text(WIN, c[colname].iloc[-1], f" {name}", fontsize=8, va="center")
    axes[0].set(xlabel="minutes since the start of the rule", ylabel="median A / A_start")
    axes[1].set(xlabel="minutes since the start of the rule", ylabel="median A (dB)")
    axes[0].legend(frameon=False, fontsize=8.5)
    for ax in axes:
        ax.set_xlim(0, WIN + 5)
        style_axes(ax)
    fig.suptitle("Step 4c, A3: the median tail under the three start rules (S0 solid, S1 dashed, S2 dotted). "
                 "S2 starts every tail at its maximum, so its curve is biased towards decaying.", fontsize=10)
    fig.tight_layout()
    return fig


def figure_gallery(picks, st, TL, col, i_of, base):
    after = S["gallery"]["after_min"]
    fig, axes = plt.subplots(2, len(picks) // 2, figsize=(4.2 * (len(picks) // 2), 7.6), squeeze=False)
    for ax, r in zip(axes.ravel(), picks.itertuples()):
        i0 = i_of(r.t0_trace)
        x = np.arange(-BEFORE, after + 1)
        a = TL[i0 - BEFORE:i0 + after + 1, col[r.sublink_id]] - base[r.tail_id]
        ax.plot(x, a, "o", ms=2.2, color="#8a8a8a")
        ax.plot(x, pd.Series(a).rolling(5, center=True, min_periods=3).median(), "-", color=INK, lw=1.5)
        s = st.loc[r.tail_id]
        for j, name in enumerate(RULES):
            xk = s[f"{name}_k"] - s.S0_k
            ax.axvline(xk, color=INK, ls=LINE[name], lw=1.4)
            ax.text(xk, 0.99 - 0.09 * j, f" {name}", transform=ax.get_xaxis_transform(), fontsize=8, va="top",
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.4, alpha=0.85))
        ax.axhline(0, color=INK2, lw=0.6)
        ax.set_title(f"4b class: {r.class_4b.upper()}\ntail {r.tail_id}, event {r.event_id[1:9]}, L = {r.L_km:.1f} km", fontsize=8)
        ax.set_xlim(-BEFORE, after)
        style_axes(ax)
    for ax in axes[-1]:
        ax.set_xlabel("minutes since t0_trace")
    for ax in axes[:, 0]:
        ax.set_ylabel("A (dB)")
    fig.suptitle(f"Step 4c, A3: {len(picks)} TRAIN tails (two per 4b class, seed {S['gallery']['seed']}) with the three starts. Dots: A, 1 min; "
                 "line: 5 min median. Vertical lines: S0 solid, S1 dashed, S2 dotted (labelled).", fontsize=10)
    fig.tight_layout()
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    clock = time.time()

    # ------------------------------------------------------------------------------------------- data (TRAIN only)
    dates = ["t0", "t0_trace", "t_first_dry", "t_end", "t_use_end_lenient"]
    sel = pd.read_csv(P.TAILS_SELECTED, parse_dates=dates, low_memory=False).query("split_event == 'train'")   # VAL rows dropped here
    ok = (sel.exclusion_reason.fillna("") == "") & ~(sel.offset_std_unc > B2B["max_offset_std_db"]) & (sel.length_min > 0)
    pool = sel[ok].sort_values("tail_id").reset_index(drop=True)            # the TRAIN pairs that pass the v1 exclusions
    base = sel.set_index("tail_id").baseline_db.to_dict()
    main_ix = pd.read_csv(P.TAILS_MODEL_INDEX, parse_dates=["t0", "t0_trace", "t_end"]).query("split_event == 'train'").sort_values("tail_id").reset_index(drop=True)
    alt_ix = pd.read_csv(P.TAILS_MODEL_ALT_INDEX, parse_dates=["t0", "t0_trace", "t_end"]).query("split_event == 'train'").sort_values("tail_id").reset_index(drop=True)
    cls4b = pd.read_csv(P.TAIL_CLASSES)
    assert set(cls4b.tail_id) == set(main_ix.tail_id) and (main_ix.split_event == "train").all() and (alt_ix.split_event == "train").all()
    long = pq.read_table(P.TAILS_1MIN, columns=["tail_id", "minutes_since_t0", "A_db"], filters=[("tail_id", "in", pool.tail_id.tolist())]).to_pandas()
    long = long.sort_values(["tail_id", "minutes_since_t0"], kind="stable").reset_index(drop=True)
    long_a = long.A_db.to_numpy("float64")
    bounds = np.column_stack([np.searchsorted(long.tail_id.to_numpy(), pool.tail_id.to_numpy(), side=s_) for s_ in ("left", "right")])

    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    subs = sorted(alt_ix.sublink_id.unique())
    s1 = d.store("cml_1min")
    T1 = pd.DatetimeIndex(s1.time.values)
    TL = s1.tl_mean.sel(sublink_id=subs).load().to_numpy()                # (minutes of the season, sub-links): 1 min tl_mean
    col = {s_: j for j, s_ in enumerate(subs)}
    i_of = lambda t: int((t - T1[0]) / MIN)                               # noqa: E731  minute index of a time
    ra = d.store("radar_along_links_5min")
    T5 = pd.DatetimeIndex(ra.time.values)
    assert T5[0] == T1[0], "the 1 min and 5 min stores must start at the same time"
    rmean5 = ra.R_mean.load().to_pandas()
    r_col = {c: j for j, c in enumerate(rmean5.columns)}
    wx = pd.read_csv(P.DEW_RISK_HOURLY, parse_dates=["time_utc"]).set_index("time_utc")
    night_of = wx.night.astype(bool).to_dict()
    print(f"TRAIN only: {len(main_ix)} main tails, {len(alt_ix)} alternative tails, {len(pool)} pairs in the pool of A3; loaded in {time.time() - clock:.0f} s", flush=True)

    # A(t) built here from the 1 min store must be the A of the frozen set
    chk = pq.read_table(P.TAILS_MODEL, filters=[("split_event", "=", "train")], columns=["tail_id", "sublink_id", "t", "A_db"]).to_pandas().head(20000)
    mine = np.array([TL[i_of(t), col[s_]] - base[k] for k, s_, t in zip(chk.tail_id, chk.sublink_id, chk.t)])
    a_check = float(np.nanmax(np.abs(mine - chk.A_db.to_numpy())))
    assert a_check < 1e-3, f"A(t) of this script differs from the frozen set by {a_check} dB"

    # ------------------------------------------------------------------------------------------- A1
    T = main_ix.merge(cls4b[["tail_id", "class", "rise_db", "night", "ev_dew_risk", "humidity_pct"]].rename(columns={"rise_db": "rise_4b_db"}), on="tail_id")
    st5 = d.store("link_reference_state_5min").state.sel(cml_id=sorted(T.cml_id.unique())).load()
    assert pd.DatetimeIndex(st5.time.values)[0] == T1[0]
    st = st5.to_numpy()
    idx = np.where(st == 1, np.arange(len(st))[:, None], -1)
    state = {"state": st, "since_wet_h": (np.arange(len(st))[:, None] - np.maximum.accumulate(idx, axis=0)) * 5 / 60}   # no wet bin yet: since the season start
    st_col = {int(c): j for j, c in enumerate(st5.cml_id.values)}
    jumps = pd.read_csv(P.DRY_JUMPS, parse_dates=["start_utc", "end_utc"])
    in_jump = {s_: np.zeros(len(T1), bool) for s_ in T.sublink_id.unique()}
    for j in jumps[jumps.sublink_id.isin(in_jump)].itertuples():
        in_jump[j.sublink_id][max(i_of(j.start_utc), 0):i_of(j.end_utc) + 1] = True
    n_days = len(T1) // 1440
    plan_file, a1_file = CACHE / "plan.json", CACHE / "a1_per_tail.csv"
    if plan_file.exists() and a1_file.exists():
        plan, a1 = json.loads(plan_file.read_text()), pd.read_csv(a1_file)
        print(f"A1 from the cache ({rel(a1_file)}); plan: cut {plan['cuts'] or 'nothing'}")
    else:
        t0 = time.time()
        first = a1_controls(T.head(A1["timing_tails"]), TL, col, i_of, state, st_col, in_jump, n_days)
        per_tail = (time.time() - t0) / A1["timing_tails"]
        proj = per_tail * len(T) + 240                     # A2 and A3 have no fit: about 4 min together (bootstraps, figures)
        plan = {"a1_seconds_per_tail": round(per_tail, 3), "projected_s": round(proj), "cuts": []}
        print(f"timing: A1 takes {per_tail:.2f} s per tail on {A1['timing_tails']} tails -> the whole step is projected at "
              f"{proj / 60:.0f} min (budget {S['time_budget_h']} h)", flush=True)
        for cut in ("A2 on the alternative set", "the night / day split of A1 by dew risk", "S2"):
            if proj <= S["time_budget_h"] * 3600:
                break
            plan["cuts"].append(cut)
            proj *= 0.8
        a1 = pd.concat([first, a1_controls(T.iloc[A1["timing_tails"]:], TL, col, i_of, state, st_col, in_jump, n_days)], ignore_index=True)
        a1.to_csv(a1_file, index=False)
        plan_file.write_text(json.dumps(plan))
    T = T.merge(a1, on="tail_id")
    rise_check = float((T.rise_db - T.rise_4b_db).abs().max())             # the rise here must be the rise of step 4b
    assert rise_check < 2e-3, f"the rise of a tail differs from step 4b by {rise_check} dB"
    amp = pd.read_csv(P.DAILY_CYCLE).set_index("sublink_id").amplitude_db
    T["daily_cycle_amplitude_db"] = T.sublink_id.map(amp)
    T["dew"] = T.ev_dew_risk.map({True: "yes", False: "no", "True": "yes", "False": "no"}).fillna("unknown")
    T["night"] = T.night.astype(bool)
    save(T[["tail_id", "sublink_id", "cml_id", "event_id", "t0_trace", "L_km", "class", "night", "dew", "humidity_pct", "daily_cycle_amplitude_db"]
           + [c for c in a1.columns if c != "tail_id"]], P.CHAR_A1_TAILS, index=False)

    rows = [a1_row("all tails", T)] + [a1_row(f"class: {c}", T[T["class"] == c]) for c in CLASSES]
    rows += [a1_row("night", T[T.night]), a1_row("day", T[~T.night])]
    rows += [a1_row(f"dew risk: {v}", T[T.dew == v]) for v in ("yes", "no", "unknown")]
    rows += [a1_row(f"class: {c}, {'night' if n_ else 'day'}", T[(T["class"] == c) & (T.night == n_)]) for c in CLASSES for n_ in (True, False)]
    if "the night / day split of A1 by dew risk" not in plan["cuts"]:
        rows += [a1_row(f"{'night' if n_ else 'day'}, dew risk: {v}", T[(T.night == n_) & (T.dew == v)]) for n_ in (True, False) for v in ("yes", "no", "unknown")]
    a1_tab = pd.DataFrame(rows)
    save(a1_tab, P.CHAR_A1, index=False)
    g1 = a1_tab.set_index("group")
    rises = T[(T["class"] == "rises") & T.has_control]
    share_rises = float(rises.exceeds_control.astype(bool).mean()) if len(rises) else np.nan
    a1_verdict = ("no control days: not decided" if not len(rises) else "rises are mostly drift" if share_rises <= A1["drift_max_share"]
                  else "rises are mostly water (not drift)" if share_rises >= A1["water_min_share"] else "mixed")
    exc = T[T.has_control & (T.exceeds_control == True)]                    # noqa: E712
    rn, rd = g1.loc["class: rises, night"], g1.loc["class: rises, day"]     # Q1 asks about the NIGHT-time rises
    dew_exc, dew_dec = share(exc.ev_dew_risk), share(T[T["class"] == "decays"].ev_dew_risk)
    dew_exc_rises = share(exc[exc["class"] == "rises"].ev_dew_risk)
    dew_like = bool(dew_exc > A1["dew_ratio"] * dew_dec) if np.isfinite(dew_exc) and np.isfinite(dew_dec) else None
    amp_tab = T.groupby("class").daily_cycle_amplitude_db.agg(["count", "median", "max"]).reindex(CLASSES).round(3)
    ctrl_drift = T[T.has_control].groupby("night").agg(tails=("tail_id", "size"), control_end_change_median_db=("control_end_change_median_db", "median"),
                                                       control_end_change_p5_db=("control_end_change_p5_db", "median"),
                                                       control_end_change_p95_db=("control_end_change_p95_db", "median"),
                                                       control_p95_of_rise_db=("control_p95_db", "median")).round(3).rename(index={True: "night", False: "day"})
    pos = jumps[jumps.sign == "positive"].copy()
    hour = pos.start_utc.dt.floor("1h")
    pos["night"], pos["rh_known"], pos["dew"] = hour.map(night_of), hour.map(wx.humidity.notna().to_dict()), hour.map(wx.dew_risk.astype(bool).to_dict())
    neg = jumps[jumps.sign == "negative"]
    jump_ctx = {"positive": len(pos), "negative": len(neg), "sublinks": int(pos.sublink_id.nunique()),
                "night": float(pos.night.astype(bool).mean()), "night_all_hours": float(wx.night.astype(bool).mean()),
                "dew": float(pos[pos.rh_known == True].dew.astype(bool).mean()), "rh_unknown": int((pos.rh_known != True).sum()),   # noqa: E712
                "dew_all_hours": float(wx[wx.humidity.notna()].dew_risk.astype(bool).mean()),
                "neg_night": float(neg.start_utc.dt.floor("1h").map(night_of).astype(bool).mean())}
    gap = {e: (int((T.event_id == e).sum()), int(((T.event_id == e) & (T.dew == "unknown")).sum()),
               int(((T.event_id == e) & (T["class"] == "rises")).sum()),
               share(T[(T.event_id == e) & T.has_control].exceeds_control)) for e in FRZ["look_events"]}
    print(f"A1: {int(T.has_control.sum())} of {len(T)} tails have a control; class 'rises': {share_rises:.0%} above their dry days -> {a1_verdict}", flush=True)

    # ------------------------------------------------------------------------------------------- A2
    rng = np.random.default_rng(A2["seed"])
    R5 = rmean5.to_numpy()
    sets = {"main set (A_start >= 1 dB)": main_ix}
    if "A2 on the alternative set" not in plan["cuts"]:
        sets[f"alternative set (A_start >= {S2SEL['alt_a_start_db']} dB)"] = alt_ix
    a2_tabs, a2_info = {}, {}
    for name, ix_ in sets.items():
        a2_tabs[name], a2_info[name] = a2_set(name, ix_, TL, col, i_of, base, R5, r_col, rng)
        print(f"A2, {name}: power {a2_info[name]['power']}, N = {a2_info[name]['N']}", flush=True)
    save(pd.concat(a2_tabs.values(), ignore_index=True), P.CHAR_A2, index=False)

    # ------------------------------------------------------------------------------------------- A3
    st_file = CACHE / "a3_starts.csv"
    if st_file.exists():
        starts = pd.read_csv(st_file)
    else:
        ref = reference_traces(d, pool.cml_id.unique(), B2B["gauge_max_dist_m"], A3["s1_near_km"], A3["s1_city_gauge_km"], A3["s1_smhi_gauge_km"])
        starts = start_rules(pool, long_a, bounds, T5, ref["trace"], ref["nan"], ref["broad"])
        starts.to_csv(st_file, index=False)
    s0_same = set(starts.tail_id[starts.S0_kept]) == set(main_ix.tail_id)
    assert s0_same, "S0 with the v1 rules must give the frozen TRAIN main set"
    rng3 = np.random.default_rng(A3["seed"])
    use_rules = [n for n in RULES if not (n == "S2" and "S2" in plan["cuts"])]
    a3_rows, a3_rates, a3_curves, kept, cls_by = [], [], [], {}, {}
    for name in use_rules:
        row, rate, curve, kept[name], cls_by[name] = describe_rule(name, starts, pool, long_a, bounds, night_of, rng3)
        a3_rows.append(row), a3_rates.append(rate), a3_curves.append(curve)
        print(f"A3, {RULES[name]}: {row['tails']} tails, slope of k(A) {row['k_slope']} ({row['k_slope_lo']}..{row['k_slope_hi']})", flush=True)
    a3 = pd.DataFrame(a3_rows)
    save(a3, P.CHAR_A3, index=False)
    save(pd.concat(a3_rates, ignore_index=True), P.CHAR_A3_RATE, index=False)
    curves = pd.concat(a3_curves, ignore_index=True)
    save(curves, P.CHAR_A3_CURVES, index=False)
    cls_same = float((cls_by["S0"].set_index("tail_id")["class"].reindex(cls4b.tail_id).to_numpy() == cls4b["class"].to_numpy()).mean())
    overlap = pd.DataFrame([{"start_rule": RULES[a_], **{f"also kept by {b_}": len(kept[a_] & kept[b_]) for b_ in use_rules}} for a_ in use_rules])
    # the product: one row per TRAIN pair of the pool
    prod = pool[["tail_id", "sublink_id", "cml_id", "event_id", "t0", "L_km", "frequency_ghz"]].merge(starts, on="tail_id")
    for name in RULES:
        prod[f"{name}_start"] = prod.t0 + pd.to_timedelta(prod[f"{name}_k"], unit="min")
        prod[f"{name}_end"] = prod.t0 + pd.to_timedelta(prod[f"{name}_end_k"], unit="min")
    prod["in_v1_main"], prod["in_v1_alt"] = prod.tail_id.isin(main_ix.tail_id), prod.tail_id.isin(alt_ix.tail_id)
    cols = ["tail_id", "sublink_id", "cml_id", "event_id", "t0", "L_km", "frequency_ghz", "in_v1_main", "in_v1_alt"]
    for name in RULES:
        cols += [f"{name}_start", f"{name}_end", f"{name}_usable_min", f"{name}_truncated_by", f"{name}_A_start_db", f"{name}_kept"]
    save(prod[cols + ["s2_peak_above_A_start_db", "s2_rise_before_onset", "s2_peak_above_level_at_t0_trace_db"]], P.TAIL_STARTS, index=False)
    sm_ = starts.set_index("tail_id")
    in_main = sm_.loc[main_ix.tail_id]
    s1_delay, s2_delay = in_main.S1_k - in_main.S0_k, in_main.S2_k - in_main.S0_k
    s2_flag = float(in_main.s2_rise_before_onset.mean())
    s2_flag0 = float((in_main.s2_peak_above_level_at_t0_trace_db >= A3["s2_rise_flag_db"]).mean())
    s2_late = float((s2_delay >= 5).mean())
    lost = {n: (len(kept["S0"] - kept[n]), len(kept[n] - kept["S0"])) for n in use_rules if n != "S0"}
    why_lost = {n: starts[starts.tail_id.isin(kept["S0"] - kept[n])].pipe(
        lambda q: {"A_start below 1 dB at the new start": int((q[f"{n}_A_start_db"].fillna(-1) < S2SEL["main_a_start_db"]).sum()),
                   "shorter than 60 min": int(((q[f"{n}_A_start_db"].fillna(-1) >= S2SEL["main_a_start_db"]) & (q[f"{n}_usable_min"] < S2SEL["min_usable_min"])).sum())})
        for n in lost}

    # ------------------------------------------------------------------------------------------- figures
    fig = figure_a1(T)
    save(fig, P.CHAR_FIG_A1, dpi=110)
    plt.close(fig)
    fig = figure_a2(a2_tabs, a2_info)
    save(fig, P.CHAR_FIG_A2, dpi=110)
    plt.close(fig)
    fig = figure_a3(curves)
    save(fig, P.CHAR_FIG_A3, dpi=110)
    plt.close(fig)
    per = S["gallery"]["tails"] // len(CLASSES)
    picks = pd.concat([T[T["class"] == c].sample(per, random_state=S["gallery"]["seed"]) for c in CLASSES]).rename(columns={"class": "class_4b"})
    fig = figure_gallery(picks, sm_, TL, col, i_of, base)
    save(fig, P.CHAR_FIG_GALLERY, dpi=105)
    plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    main_name, alt_name = list(sets)[0], list(sets)[-1]
    im, ia = a2_info[main_name], a2_info[alt_name]
    tm = a2_tabs[main_name].set_index("minute")
    show_min = [-30, -20, -10, -5, 0, 2, 5, 10, 15, 20, 30, 45, 60]
    a2_show = pd.concat([t_[t_.minute.isin(show_min)] for t_ in a2_tabs.values()])[
        ["set", "minute", "tails", "mean_A_db", "beta_db_per_km", "beta_lo", "beta_hi", "expected_path_rain_db_per_km", "expected_path_rain_db", "alpha_high_band_db", "path_part_median_link_db"]].round(3)
    info_tab = pd.DataFrame(a2_info.values())
    r0 = a3.set_index("start_rule")
    row_of = lambda n: r0.loc[RULES[n]]                                    # noqa: E731
    a3_show = a3.set_index("start_rule").T.rename_axis("quantity").reset_index()
    n_txt = lambda i: "not reached within 60 min" if i["N"] is None else f"{i['N']}"   # noqa: E731
    power_txt = lambda i: ("POWER: yes" if i["power"] else "NO POWER")     # noqa: E731
    drop = lambda t_, a_, b_: float(t_.loc[a_, "alpha_high_band_db"] - t_.loc[b_, "alpha_high_band_db"])   # noqa: E731
    slopes = {n: (row_of(n).k_slope, row_of(n).k_slope_lo, row_of(n).k_slope_hi) for n in use_rules}
    pre_m = a2_tabs[main_name][a2_tabs[main_name].minute < 0]
    pre_exp_lo, pre_exp_hi = float(pre_m.expected_path_rain_db_per_km.min()), float(pre_m.expected_path_rain_db_per_km.max())
    pre_db_lo, pre_db_hi = float(pre_m.expected_path_rain_db.min()), float(pre_m.expected_path_rain_db.max())
    half_width = float(((pre_m.beta_hi - pre_m.beta_lo) / 2).median())
    post_m = a2_tabs[main_name][a2_tabs[main_name].minute >= 0]
    beta_hi_max = float(post_m.beta_hi.max())
    nd = {n: 1 - row_of(n)["share_decays"] for n in use_rules}
    rises_all = T[T["class"] == "rises"]
    rises_gap = int(rises_all.event_id.isin(FRZ["look_events"]).sum())
    med_rise, med_p95 = float(rises.rise_db.median()), float(rises.control_p95_db.median())
    all_share = float(T[T.has_control].exceeds_control.astype(bool).mean())
    below_s1 = why_lost["S1"]["A_start below 1 dB at the new start"]
    report = f"""# wet_antenna step 4c (Phase A): what is in a drying tail? (TRAIN tails of v1, no law fit)

Written by `research/wet_antenna/characterise_tails.py`. Settings and every rule: `wet_antenna.characterise_tails` in
`research/config.yaml`, written before the run. The plan this step belongs to: `docs/ROADMAP.md`, section "Plan change
(2026-10-12, after 4b): the signal before the law".

## 1. What was run

- Data: the TRAIN tails of the frozen set v1: {len(main_ix)} main tails ({main_ix.event_id.nunique()} events, {main_ix.cml_id.nunique()} links), the
  alternative set ({len(alt_ix)} tails, A_start >= {S2SEL['alt_a_start_db']} dB), and for the start rules the {len(pool)} TRAIN pairs that pass the v1
  exclusions. VAL rows are dropped when a file is read; no VAL tail is used. A(t) = 1 min `tl_mean` - `baseline_db`
  (frozen baseline v2), read from the release for 30 min before t0_trace as well; largest difference from the A of
  the frozen set: {a_check:.1e} dB.
- No drying law is fitted anywhere in this step. Tail set v1, the baseline and all earlier outputs are unchanged.
- Timing: A1 took {plan['a1_seconds_per_tail']:.2f} s per tail on the first {A1['timing_tails']} tails; projection {plan['projected_s'] / 60:.0f} min against a budget of
  {S['time_budget_h']} h. Cut: {', '.join(plan['cuts']) or 'nothing'}. Run time of this run: {time.time() - clock:.0f} s (A1 and the start rules are cached).
- Checks inside the script: the rise of a tail equals the rise of step 4b (largest difference {rise_check:.1e} dB); the start
  rule S0 with the v1 end rules gives exactly the frozen TRAIN main set ({'yes' if s0_same else 'NO'}); the classes recomputed from S0
  equal the 4b classes for {cls_same:.1%} of the tails.

## 2. Part 0: the revised plan

Recorded before the analysis: `docs/ROADMAP.md` (the new plan-change section with Q1-Q5 and Phases A-E, three new
standing rules, status rows 4c and 4d, the order line, the open questions pointing to Q1-Q3, one checklist item),
`research/wet_antenna/METHOD.md` (Q1-Q3 at the top of the open issues) and `research/wet_antenna/README.md` (status).

## 3. A1: is the rise of a tail more than the daily cycle? (Q1)

For each of the {len(T)} main tails: the same sub-link, the same clock window (60 min from the clock time of t0_trace), on
every other day of the season that is certified dry: link_state dry in the whole window, at least
{A1['min_hours_since_wet']} h since the last wet bin of the link, not inside a run of `dry_jumps.csv`, at least {A1['min_valid_share']:.0%} of the
minutes present. On each such day the same statistic as for the tail: rise = maximum of the 5 min median - mean of
the first 5 min. A tail "exceeds" when its rise is above the {A1['control_percentile']}th percentile of its own control days
(at least {A1['min_control_days']} days, otherwise "no control"). {int(T.has_control.sum())} of {len(T)} tails have a control (median
{T[T.has_control].control_days.median():.0f} control days). By chance alone about {100 - A1['control_percentile']} % would exceed.

{md_table(a1_tab)}

(`{rel(P.CHAR_A1)}`; per tail: `{rel(P.CHAR_A1_TAILS)}`.)

![A1](../figures/{P.CHAR_FIG_A1.name})

**What dry days do in the same clock window** (median over the tails of each tail's control statistics): how far the
level moves from the first to the last 5 min of the window, and the 95th percentile of the rise.

{md_table(ctrl_drift, index=True)}

**Verdict by the rule fixed before the run: {a1_verdict.upper()}.** {int(rises.exceeds_control.astype(bool).sum())} of the {len(rises)} class-"rises" tails with a
control are above their dry-day {A1['control_percentile']}th percentile ({share_rises:.0%}; rule: at most {A1['drift_max_share']:.0%} = drift, at least
{A1['water_min_share']:.0%} = water). The night-time rises alone: {int(rn.exceed)} of {int(rn.with_control)} ({rn.share_exceed:.0%}); by day: {int(rd.exceed)} of {int(rd.with_control)}
({rd.share_exceed:.0%}). On the night-time rises alone the same rule would say "mixed"; at night the dry days themselves move
more (table above).

**Dew-like?** Dew risk (where the humidity is known) of the tails that exceed: {dew_exc:.0%} (class "rises" only: {dew_exc_rises:.0%});
of the decaying tails: {dew_dec:.0%}. Rule (more than {A1['dew_ratio']:g} times): **{'dew-like' if dew_like else 'not dew-like by the rule' if dew_like is not None else 'not decided'}**.
The rule can hardly be met: rain leaves the air humid, so the hourly dew-risk flag is already true for most decaying
tails, and a share cannot be more than twice a share above 50 %.

**The humidity gap.** Humidity at GbgA is missing for {int(wx.humidity.isna().sum())} hours of the season, all in June.
{'; '.join(f'{e}: {n} main tails, {u} with unknown dew risk, {r_} of class "rises", {s_:.0%} of its tails with a control above their dry days' for e, (n, u, r_, s_) in gap.items())}.
So for the two events with most rising tails "dew or not" cannot be read from the humidity.

**Daily-cycle amplitude (R1 C3) of the sub-links, per class** (dB; v1 already excludes sub-links above
{CFG['baseline']['daily_cycle_large_db']} dB through `use_in_waa_main`):

{md_table(amp_tab, index=True)}

**Context from `dry_jumps.csv`** (dry-weather runs of >= 3 h at >= 1 dB from the dry median, all train sub-links):
{jump_ctx['positive']} positive jumps on {jump_ctx['sublinks']} sub-links ({jump_ctx['negative']} negative). {jump_ctx['night']:.0%} of the positive jumps
start at night (night is {jump_ctx['night_all_hours']:.0%} of all hours; negative jumps: {jump_ctx['neg_night']:.0%}), and {jump_ctx['dew']:.0%} start in an hour with
dew risk (dew risk holds in {jump_ctx['dew_all_hours']:.0%} of the hours with humidity; {jump_ctx['rh_unknown']} jumps start in an hour without humidity).
So these antennas do gain attenuation at night in humid air without any rain.

## 4. A2: antenna or path rain? (Q2)

Per minute m from -{BEFORE} to +{WIN} around t0_trace, across the tails of a set:
A(m) = alpha_event(m) + beta(m) L + gamma(m) [f >= {A2['band_split_ghz']:g} GHz] + error, with one intercept per event. beta(m) is the dB
per km of path at that minute (path rain grows with L; the antenna does not). 90 % interval of beta from
{A2['boot_reps']} event bootstraps. The dashed reference is the beta expected from the radar: the mean over the tails of
a R_mean^b (ITU-R P.838-3 at the frequency and polarization of the sub-link, R_mean = radar along the path).

{md_table(info_tab)}

- **Power check** (the interval of beta above 0 at some minute in {A2['power_minutes'][0]}..{A2['power_minutes'][1]}): main set **{power_txt(im)}**
  ({im['minutes_with_beta_above_0_in_power_range']} of {A2['power_minutes'][1] - A2['power_minutes'][0] + 1} minutes; largest beta {im['max_beta_in_power_range']} dB/km); alternative set **{power_txt(ia)}**
  ({ia['minutes_with_beta_above_0_in_power_range']} minutes; largest beta {ia['max_beta_in_power_range']} dB/km).
- **N by the rule** (first minute >= 0 from which the interval of beta contains 0 for {A2['run_min']} minutes in a row): main set
  **N = {n_txt(im)}**, alternative set **N = {n_txt(ia)}**. {'Where the power check fails, N is not interpreted (rule).' if not (im['power'] and ia['power']) else ''}
- Minutes >= 0 with the interval of beta above 0: main {im['minutes_ge0_with_beta_above_0']}, alternative {ia['minutes_ge0_with_beta_above_0']}; below 0:
  main {im['minutes_ge0_with_beta_below_0']}, alternative {ia['minutes_ge0_with_beta_below_0']}.
- **Against the expected path rain** (minutes before 0): the mean observed beta is {im['beta_over_expected_before_0']} times the expected
  one in the main set ({ia['beta_over_expected_before_0']} in the alternative set); the expected value lies inside the interval in
  {im['expected_inside_interval_before_0']:.0%} of those minutes (alternative {ia['expected_inside_interval_before_0']:.0%}). After t0_trace the expected beta is 0 by
  construction of v1 (no radar on the path).
- **alpha(m), the part shared by the links of an event** (>= 25 GHz, i.e. A at L = 0), main set: {tm.loc[-10, 'alpha_high_band_db']:.2f} dB at
  -10 min, {tm.loc[0, 'alpha_high_band_db']:.2f} at 0, {tm.loc[10, 'alpha_high_band_db']:.2f} at 10, {tm.loc[30, 'alpha_high_band_db']:.2f} at 30, {tm.loc[60, 'alpha_high_band_db']:.2f} at 60 min. It falls by
  {drop(tm, 0, 10):.2f} dB in the first 10 min and by {drop(tm, 10, 60):.2f} dB from 10 to 60 min; the path part of a median link
  (beta x {im['median_L_km']} km) is {tm.loc[0, 'path_part_median_link_db']:+.2f} dB at 0 and {tm.loc[10, 'path_part_median_link_db']:+.2f} dB at 10 min.
- **Why there is no power** (information): the path attenuation the radar expects in the last 30 min before t0_trace is
  small, {pre_exp_lo:.2f} to {pre_exp_hi:.2f} dB/km, i.e. {pre_db_lo:.2f} to {pre_db_hi:.2f} dB for the mean tail, against a mean A of {tm.loc[-10, 'mean_A_db']:.2f} dB;
  the 90 % interval of beta is about +-{half_width:.2f} dB/km wide there. The tails end with light rain, so the minutes before
  t0_trace are not a period of strong path rain, and the regression cannot see an effect of that size.

Selected minutes (`{rel(P.CHAR_A2)}` has every minute):

{md_table(a2_show)}

![A2](../figures/{P.CHAR_FIG_A2.name})

## 5. A3: when does drying start? (Q3)

Three start rules, none uses a law fit; for each the v1 end rules are applied again from the new start (next rain
trace on the path, rise cut v2 from the new start, the step 2b ends) and the main-set rule (A_start >= {S2SEL['main_a_start_db']} dB at the
new start, usable length >= {S2SEL['min_usable_min']} min).
- **S0** = t0_trace of v1.
- **S1** = strict dry: the end of the last 5 min bin, between t0 and the first dry bin, with the v1 trace on the
  path OR radar R_max > 0 on the path of a link whose midpoint is within {A3['s1_near_km']:g} km OR a City gauge within
  {A3['s1_city_gauge_km']:g} km of the path > 0 OR the SMHI gauge > 0 if within {A3['s1_smhi_gauge_km']:g} km; plus {A3['s1_buffer_min']} min. The release has no
  radar cells around a path, so "radar nearby" is the radar on the paths of the links within {A3['s1_near_km']:g} km, as in step 4b.
- **S2** = signal onset: the minute of the maximum of the 5 min median of A in [t0_trace, t0_trace + {A3['s2_search_min']} min]
  (never beyond the next rain trace on the path or the step 2b end). **Warning: S2 makes every tail start at a
  maximum, so the shape after the start is biased towards decaying. It is a sensitivity case, not a drying-law set.**

{md_table(a3_show)}

(`{rel(P.CHAR_A3)}`; k(A) per level: `{rel(P.CHAR_A3_RATE)}`; median curves: `{rel(P.CHAR_A3_CURVES)}`; the starts of every
pair: `{rel(P.TAIL_STARTS)}`.) The k(A) slope is the B5 method of step 4b: A smoothed over {RR['smooth_min']} min, central
difference, k = -mean dA/dt / median A per bin, bins with >= {A3['rate_min_minutes']} tail-minutes, slope of log k vs log A from
{W_LO} dB up; interval from {A3['rate_boot_reps']} event bootstraps. The delay is the time from t0_trace to the start of the rule.

How the kept sets overlap:

{md_table(overlap)}

- Among the {len(main_ix)} v1 main tails, S1 starts {s1_delay.median():.0f} min later (median; quartiles {s1_delay.quantile(0.25):.0f}-{s1_delay.quantile(0.75):.0f}) and S2
  {s2_delay.median():.0f} min later (quartiles {s2_delay.quantile(0.25):.0f}-{s2_delay.quantile(0.75):.0f}; S2 = S0 for {float((s2_delay == 0).mean()):.0%}). "Rise before onset" (peak
  >= {A3['s2_rise_flag_db']} dB above A_start at t0_trace): {s2_flag:.0%} of the v1 main tails. That flag is inflated: A_start is a mean
  over 5 min, so the first minute of a fast-falling tail is already above it. Against the 5 min median AT t0_trace
  (added after the first run, information) the peak is >= {A3['s2_rise_flag_db']} dB higher in {s2_flag0:.0%}; S2 is 5 min or more after S0 in {s2_late:.0%}.
- {'; '.join(f"{RULES[n]} loses {a_} of the v1 main tails ({', '.join(f'{k}: {v}' for k, v in why_lost[n].items())}) and adds {b_} new ones" for n, (a_, b_) in lost.items())}.

![A3 median curves](../figures/{P.CHAR_FIG_A3.name})

![A3 gallery](../figures/{P.CHAR_FIG_GALLERY.name})

### Proposal for tail set v2 (for the learners; nothing is frozen here)

No law fit was used for anything below.

- **Proposal: keep S0 (t0_trace) as the start of the tail, and store S1 and S2 as two more times per tail**
  (`{rel(P.TAIL_STARTS)}`), instead of cutting the tails again at one of them. Step 4d then treats the start of
  drying as a quantity per tail that lies between S2 (the earliest it can be read from the signal) and later, and
  uses S1 as the marker "from here no reference within {A3['s1_near_km']:g} km shows rain".
- **Why not S1 as the start.** S1 comes {s1_delay.median():.0f} min after t0_trace (median), and by then most of the signal is gone:
  {below_s1} of the {len(main_ix)} v1 tails are below {S2SEL['main_a_start_db']} dB at S1. What is left ({int(row_of('S1')['tails'])} tails, {int(row_of('S1')['events'])} events) is the tails
  that are STILL high that long after the rain, i.e. the slow ones: a selection on the outcome ({row_of('S1')['share_decays']:.0%} "decays",
  t_half {row_of('S1')['t_half_min']:.0f} min). S1 cannot be the start of a set that also has the {S2SEL['main_a_start_db']} dB rule.
- **Why not S2 as the start.** It starts every tail at its maximum, so decay after the start is partly built in
  (the warning above). It changes little for most tails anyway (S2 = S0 for {float((s2_delay == 0).mean()):.0%} of the v1 tails), and it adds
  {lost['S2'][1]} tails whose level reaches {S2SEL['main_a_start_db']} dB only after t0_trace, i.e. tails that RISE after the last rain trace.
- **How the conclusions of steps 4 and 4b change with the start rule.**
  - Step 4, "the rate per dB falls with the level": the slope of k(A) is {slopes['S0'][0]} ({slopes['S0'][1]}..{slopes['S0'][2]}) for S0,
    {slopes['S1'][0]} ({slopes['S1'][1]}..{slopes['S1'][2]}) for S1 and {slopes['S2'][0]} ({slopes['S2'][1]}..{slopes['S2'][2]}) for S2. It is negative under every rule, so the
    falling k is not made by the choice of the start; S2 only weakens it.
  - Step 4b, "44 % of the tails do not decay at once": {nd['S0']:.0%} for S0, {nd['S1']:.0%} for S1, {nd['S2']:.0%} for S2. Starting at the
    maximum removes most rising tails ({row_of('S0')['share_rises']:.0%} -> {row_of('S2')['share_rises']:.0%}) but not the waiting and the flat ones.
- **What stays stable** between S0 and S2: the median shape (A / A_start {row_of('S0')['A_over_A_start_10min']} / {row_of('S0')['A_over_A_start_30min']} / {row_of('S0')['A_over_A_start_60min']} at 10 / 30 / 60 min
  against {row_of('S2')['A_over_A_start_10min']} / {row_of('S2')['A_over_A_start_30min']} / {row_of('S2')['A_over_A_start_60min']}), t_half ({row_of('S0')['t_half_min']:.0f} and {row_of('S2')['t_half_min']:.0f} min) and the share that halves within
  60 min ({row_of('S0')['share_half_within_60min']:.0%}, {row_of('S2')['share_half_within_60min']:.0%}). S1 is different in everything, for the reason above.

## 6. In plain words: what a drying tail contains

- **Drift (the daily cycle)? No, not the rises.** On dry days the same sub-link at the same clock time moves little
  in an hour: the 95th percentile of the rise is {g1.loc['day', 'median_control_p95_db']:.2f} dB by day and {g1.loc['night', 'median_control_p95_db']:.2f} dB at night. {int(rises.exceeds_control.astype(bool).sum())} of the
  {len(rises)} rising tails rise more than that. So the first candidate of Q1, a baseline that drifts with the hour, does
  not explain most of them, and a baseline v3 with a daily-cycle correction is NOT called for by this result. For the
  night-time rises the answer is weaker: {int(rn.exceed)} of {int(rn.with_control)} are above their dry nights, so about half of them are no larger
  than what the sub-link does on a dry night. Two more cautions:
  a tail is in the class "rises" because it rises more than {CL['rise_db']} dB, which is about the dry-day percentile, so part
  of the {share_rises:.0%} is built in (median rise of these tails {med_rise:.2f} dB, median dry-day percentile {med_p95:.2f} dB); and {all_share:.0%} of
  ALL tails exceed their dry days ({g1.loc['class: decays', 'share_exceed']:.0%} of the decaying ones), against 5 % by chance: a link that has just been
  rained on is noisier than a dry one, whatever the tail does.
- **Dew? Not decided, and it cannot be decided with this humidity record.** {rises_gap} of the {len(rises_all)} rising tails are in the
  two events without humidity. Where humidity is known, the tails that exceed are not in more humid hours than the
  decaying ones ({dew_exc:.0%} against {dew_dec:.0%} with dew risk). The dry-jump table shows that these antennas do gain
  attenuation at night in humid air without rain ({jump_ctx['night']:.0%} of the positive jumps start at night), so dew stays a
  candidate.
- **Path rain in the first minutes? Not shown, and not excluded.** A does not grow with the path length at any minute
  from -{BEFORE} to +{WIN} (largest upper end of the interval after t0_trace: {beta_hi_max:.2f} dB/km), but the test did not pass its
  power check: it does not see path rain before t0_trace either, because the rain at the end of an event is light.
  By the rule nothing is concluded from N. Step 2 met the same wall.
- **From which minute is it the antenna?** Not answered. What can be said: the part shared by all links of an event
  falls from {tm.loc[0, 'alpha_high_band_db']:.2f} to {tm.loc[10, 'alpha_high_band_db']:.2f} dB in the first 10 min, and the radar expects only {pre_db_lo:.2f}-{pre_db_hi:.2f} dB of path
  attenuation before t0_trace, so by the radar the fast drop is not path rain. That rests on the radar, which is
  the reference whose blind spots (drizzle, 2 km cells, 5 min) are the question.
- **The start.** The references within {A3['s1_near_km']:g} km are dry only {s1_delay.median():.0f} min after t0_trace (median), and by then {below_s1 / len(main_ix):.0%} of
  the tails are below {S2SEL['main_a_start_db']} dB. So most of what v1 calls the drying tail happens while it still rains within {A3['s1_near_km']:g} km.
- **What changes for Step 4d.** (1) No baseline v3 for the daily cycle; the simulator needs rises of 0.3-1 dB that
  are not drift. (2) The start of drying is modelled per tail (between S2 and later) instead of re-cutting the set.
  (3) The simulator must contain nearby rain during the first 15-25 min, because that is when the real tails drop.
  (4) "tau = 10 min is the drying time" stays unproven: the length test cannot separate antenna from path here; the
  very short links of Step 7 remain the real test.

## 7. Decisions for the learners (not taken here)

1. Tail set v2: keep S0 and add S1 / S2 as columns (the proposal), or re-cut at S1 or S2, or lower the {S2SEL['main_a_start_db']} dB rule
   for an S1 set (not computed here).
2. Is "not drift" accepted for the rises, so that no baseline v3 is built?
3. Humidity for E20150617_0705 and E20150701_0030 from another station: without it "dew or drizzle" stays open.
4. A2 has no power on tails. Accept that Q2 is answered only in Step 7 (very short links, whole events), or design
   another test (for example the two directions of a link, or pairs of links that share one site)?
5. Should the minutes with rain within {A3['s1_near_km']:g} km (S0 to S1) be used for the drying law at all, be modelled with rain as
   an input (the wetting side), or be left out?

## 8. Unclear / possibly wrong

1. **The A1 verdict leans on a class defined by the same quantity.** "rises" = rise > {CL['rise_db']} dB; the dry-day 95th
   percentile is about {g1.loc['all tails', 'median_control_p95_db']:.2f} dB. The rule was fixed before the run, but a share of {share_rises:.0%} against a threshold of
   {A1['water_min_share']:.0%} is not a wide margin.
2. **Control days are dry by the references only** (radar on the path, gauges within 1 km): drizzle the radar misses
   can be in a control window too. That makes the control percentile larger, i.e. the test more cautious.
3. **A tail is at a higher level than a dry day.** With 0.3 dB RSL steps and a falling level, the 5 min median of a
   tail moves more than that of a flat dry signal. "More than dry days" is therefore not the same as "water".
4. **The first 12 h of the season** count as "no wet bin yet" from the season start, so no control window lies there.
5. **A2 is run on tails selected by A_start >= 1 dB (0.6 dB)**, a selection on the outcome at minute 0-4. If path
   rain made long links higher, short links would need more antenna water to pass, which pushes beta down. The
   negative beta around minute 0 (interval below 0 at {im['minutes_ge0_with_beta_below_0']} minutes) may come from this.
6. **alpha(m) is not "the antenna" by itself**: it is A at L = 0 in this regression. Rain that is the same for all
   links of an event and does not scale with L in the data would sit in alpha too.
7. **The ITU-R P.838-3 coefficients** are typed into `tails.py` and checked against the tabulated values at 10, 20
   and 30 GHz from memory, not against the document (it is not in `docs/papers/`). They are used only for the
   dashed reference curve.
8. **S1 uses the paths of the links within {A3['s1_near_km']:g} km as "radar nearby"** (the release has no radar cells around a
   path), and its end rules are the v1 ones (path trace only), so an S1 tail can still have rain within {A3['s1_near_km']:g} km later.
9. **S1 and S2 tails are not a frozen set**: the event flag `event_rise_heavy`, `level_shift` and `undershoot` were
   not recomputed for them.
10. **The humidity is one station** (GbgA) and hourly; "dew risk" is the same for the whole network.
"""
    save(report, P.CHAR_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {time.time() - clock:.0f} s")


if __name__ == "__main__":
    main()
