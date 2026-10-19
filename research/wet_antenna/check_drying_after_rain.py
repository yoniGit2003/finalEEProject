"""wet_antenna step 4c2: does drying really start when the rain nearby ends? A quick check, data only.

Step 4c found that rain within 5 km stops about 25 min after t0_trace, when most tails are already low. So the fast
early drop of a tail may be the antenna drying, light rain near or on the path ending, or both. This step asks: once
NO reference within 5 km shows rain (from the start S1 of step 4c on), how fast does A fall? The 129 "S1 tails" of
step 4c were slow (t_half 48 min), but that set kept only tails still above 1 dB 25 min after the rain: the slow ones.
Three ways around that bias, with NO law fit:
  A1  lower thresholds: tails from S1 with A at S1 >= 0.25 / 0.4 / 0.6 dB and at least 30 min long;
  A2  a placebo start: the same pairs started at S0 + a delay drawn from the S1 - S0 delays of OTHER pairs, with the
      same thresholds but without the condition "the rain nearby has ended". If the S1 tails are slower than the
      placebo tails, the slowness comes from the rain having ended, not from starting late;
  A3  the same level before and after S1: the model-free k(A) of the alternative-set tails on the minutes with rain
      within 5 km (S0 to S1) and on the minutes without (from S1), and the same split at the placebo start.
The rule that reads the result (SLOW / SAME speed / unclear) is in config `wet_antenna.drying_after_rain`, written
before the run. TRAIN pairs only (tails_starts.csv of step 4c holds only TRAIN pairs; VAL rows of the other files are
dropped on reading). The pass over the pairs is timed on 200 pairs first and cached in
research/.cache/wet_antenna/check_drying_after_rain/ (delete the folder to compute again).
Run from the repository root (about 1 min):
    python -m research.wet_antenna.check_drying_after_rain
"""
import json
import sys
import time
import warnings

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
from research.wet_antenna import free_curve as wc
from research.wet_antenna import paths as P
from research.wet_antenna.tails import minutes, raw_rate, reference_traces, tail_from_start

CFG = load_config()
W = CFG["wet_antenna"]
S = W["drying_after_rain"]
S2SEL, B2B, FRZ, A3C = W["selection"], W["refine_selection"], W["freeze"], W["characterise_tails"]["a3"]
RATE, RULE = S["rate"], S["rule"]
WIN, MAIN = S["window_min"], S["main_threshold_db"]
PAD = S2SEL["value_window_min"]                            # A at k min = mean over [k, k + PAD): the window needs WIN + PAD minutes
CACHE = P.CACHE / "check_drying_after_rain"
KINDS = {"S1": "S1 (no rain within 5 km)", "placebo": "placebo start (same delays)", "S0": "S0 = t0_trace"}
LINE = {"S1": "-", "placebo": "--", "S0": ":"}             # one line style per start, in every figure
EDGES = np.asarray(RATE["edges_db"], float)


# =================================================================================================== the pass over the pairs
def pair_pass(pool, long_a, bounds, ref, delay_p):
    """For every pair: the three starts (S1, placebo, S0), each with the v1 end rules counted from it; S1 also ends
    at the next rain within 5 km. Returns one row per pair and, per start, the first WIN + PAD minutes of A."""
    T5 = ref["T5"]
    win = {k: np.full((len(pool), WIN + PAD), np.nan, "float32") for k in KINDS}
    rows = []
    for i, r in enumerate(pool.itertuples()):
        a = long_a[bounds[i, 0]:bounds[i, 1]]
        L, i05 = len(a), T5.searchsorted(r.t0)
        unk = ref["trace"][r.cml_id] | ref["nan"][r.cml_id] if FRZ["radar_nan_is_unknown"] else ref["trace"][r.cml_id]
        near = unk | ref["broad"][r.cml_id]                # rain (or unknown) on the path or within 5 km
        k_len = min(minutes(r.t_use_end_lenient, r.t0), L) if pd.notna(r.t_use_end_lenient) else L
        sm10 = pd.Series(a).rolling(FRZ["rise_smooth_min"], min_periods=FRZ["rise_smooth_min"] * 2 // 3).median().to_numpy()
        k0, k1 = minutes(r.S0_start, r.t0), minutes(r.S1_start, r.t0)
        kp = min(k0 + int(delay_p[i]), L)
        nxt = np.nonzero(near[i05 + k1 // 5:])[0]
        k_near = max((k1 // 5 + int(nxt[0])) * 5, k1) if len(nxt) else 10 ** 9   # the next rain within 5 km after S1
        args = (unk, i05, k_len, r.truncated_by_lenient, FRZ, S2SEL["a_dry_window_min"])
        tails = {"S1": tail_from_start(a, sm10, k1, *args, extra={"rain_nearby": k_near}),
                 "placebo": tail_from_start(a, sm10, kp, *args), "S0": tail_from_start(a, sm10, k0, *args)}
        row = {"tail_id": r.tail_id, "S1_next_rain_nearby_k": min(k_near, L)}
        for name, t in tails.items():
            seg = a[t["k"]:min(t["end"], t["k"] + WIN + PAD)]
            win[name][i, :len(seg)] = seg
            row.update({f"{name}_k": t["k"], f"{name}_end_k": t["end"], f"{name}_usable_min": t["usable_min"],
                        f"{name}_truncated_by": t["truncated_by"], f"{name}_A_start_db": t["A_start_db"]})
        rows.append(row)
    return pd.DataFrame(rows), win


# =================================================================================================== statistics
def tail_stats(Wm, a_start, usable):
    """Per tail: the time to half (minute, and whether it was reached before the tail ended or 60 min), and
    A / A_start at the shape minutes."""
    sm = pd.DataFrame(Wm[:, :WIN].T).rolling(5, center=True, min_periods=3).median().to_numpy().T
    below = sm < S["half"] * a_start[:, None]
    reached = below.any(axis=1)
    t = np.where(reached, below.argmax(axis=1), np.minimum(usable, WIN)).astype(float)
    out = {"t_half": t, "halved": reached}
    with warnings.catch_warnings():                        # a tail shorter than k minutes has no value there (NaN)
        warnings.simplefilter("ignore", RuntimeWarning)
        for k in S["shape_minutes"]:
            out[f"ratio_{k}"] = np.nanmean(Wm[:, k:k + PAD], axis=1) / a_start
    return out


def km_half(t, halved):
    """Kaplan-Meier median of the time to half (a tail that ended before halving is censored at its end), capped at
    WIN; and the share that has halved by WIN according to the same curve."""
    if not len(t):
        return np.nan, np.nan
    order = np.lexsort((~halved, t))                       # at equal times the halvings come before the censorings
    t, e = t[order], halved[order]
    surv = np.cumprod(np.where(e, 1 - 1 / (len(t) - np.arange(len(t))), 1.0))
    hit = np.nonzero(surv <= 0.5)[0]
    return (float(t[hit[0]]) if len(hit) else float(WIN)), float(1 - surv[-1])


def rate_arrays(Wm):
    """A smoothed over 7 min per tail and its central difference (the B5 method of step 4b), as two matrices."""
    sm = pd.DataFrame(Wm.T).rolling(RATE["smooth_min"], center=True, min_periods=RATE["smooth_min"] - 2).mean().to_numpy().T
    d = np.full(sm.shape, np.nan)
    d[:, 1:-1] = (sm[:, 2:] - sm[:, :-2]) / 2
    return sm, d


def rate_table(sm, d, mask=None):
    """k per level of A: k = -mean dA/dt / median A in each bin; and the slope of log k vs log A over the bins with
    enough minutes. `mask` selects the minutes that count."""
    ok = np.isfinite(sm) & np.isfinite(d) & (True if mask is None else mask)
    a, dd = sm[ok], d[ok]
    b = np.digitize(a, EDGES, right=True) - 1
    rows = []
    for j in range(len(EDGES) - 1):
        q = b == j
        n = int(q.sum())
        med, mean = (float(np.median(a[q])), float(dd[q].mean())) if n else (np.nan, np.nan)
        rows.append({"A_bin_db": f"({EDGES[j]:g}, {EDGES[j + 1]:g}]", "tail_minutes": n, "median_A_db": med,
                     "mean_dadt_db_per_min": mean, "k_per_min": -mean / med if n else np.nan})
    tab = pd.DataFrame(rows)
    tab["used"] = (tab.tail_minutes >= RATE["min_minutes"]) & (tab.k_per_min > 0)
    u = tab[tab.used]
    return tab, (float(np.log(u.k_per_min) @ wc.slope_weights(np.log(u.median_A_db))) if len(u) >= 3 else np.nan)


def boot_rows(ev_rows, rng):
    """One event bootstrap draw: the row numbers of the tails of the drawn events (an event drawn twice counts twice)."""
    keys = list(ev_rows)
    return np.concatenate([ev_rows[keys[i]] for i in rng.integers(0, len(keys), len(keys))])


def describe(name, thr, mask, st, Wm, pool, sm, d, rng, with_boot=True):
    """One row of the A1 / A2 table: a set of tails (start `name`, threshold `thr`, rows `mask`)."""
    idx = np.nonzero(mask)[0]
    t_half, share = km_half(st["t_half"][idx], st["halved"][idx])
    tab, slope = rate_table(sm[idx], d[idx])
    lo = hi = np.nan
    if with_boot and len(idx):
        ev_rows = {e: g.to_numpy() for e, g in pd.Series(idx).groupby(pool.event_id.to_numpy()[idx])}
        vals = []
        for _ in range(RATE["boot_reps"]):
            rws = boot_rows(ev_rows, rng)
            vals.append(rate_table(sm[rws], d[rws])[1])
        lo, hi = (np.nanpercentile(vals, S["interval"]) if np.isfinite(vals).mean() > 0.5 else (np.nan, np.nan))
    row = {"threshold_db": thr, "start": name, "tails": len(idx), "events": int(pool.event_id.iloc[idx].nunique()),
           "links": int(pool.cml_id.iloc[idx].nunique()), "median_A_start_db": round(float(np.median(pool[f"{name.split()[0]}_A_start_db"].to_numpy()[idx])), 3) if len(idx) else np.nan,
           **{f"A_over_A_start_{k}min": round(float(np.nanmedian(st[f"ratio_{k}"][idx])), 3) for k in S["shape_minutes"]},
           "tails_with_60min": int(np.isfinite(st[f"ratio_{S['shape_minutes'][-1]}"][idx]).sum()),
           "t_half_min": t_half, "share_half_within_60min": round(share, 3),
           "share_censored": round(float((~st["halved"][idx] & (st["t_half"][idx] < WIN)).mean()), 3) if len(idx) else np.nan,
           "k_slope": round(slope, 2), "k_slope_lo": round(float(lo), 2), "k_slope_hi": round(float(hi), 2), "bins_in_slope": int(tab.used.sum())}
    return row, tab.assign(threshold_db=thr, start=name, slope=round(slope, 3))


# =================================================================================================== figures
def figure_curves(curves, sets):
    thr_list = sorted(curves.threshold_db.unique())
    fig, axes = plt.subplots(1, len(thr_list), figsize=(5.6 * len(thr_list), 4.8), sharey=True, squeeze=False)
    for ax, thr in zip(axes[0], thr_list):
        for name in KINDS:
            c = curves[(curves.threshold_db == thr) & (curves.start == name)]
            n = int(sets[(sets.threshold_db == thr) & (sets.start == name)].tails.iloc[0])
            ax.plot(c.minute, c.median_A_over_A_start, LINE[name], color=INK, lw=1.9, label=f"{KINDS[name]}: {n} tails")
            ax.text(WIN, c.median_A_over_A_start.iloc[-1], f" {name}", fontsize=8, va="center")
        ax.axhline(0.5, color=INK2, lw=0.7)
        ax.text(1, 0.51, "half", fontsize=7.5, va="bottom")
        ax.set(xlabel="minutes since the start", title=f"A at the start >= {thr:g} dB, at least {S['min_usable_min']} min long", xlim=(0, WIN + 9))
        ax.legend(frameon=False, fontsize=8, loc="upper right")
        style_axes(ax)
    axes[0, 0].set_ylabel("median A / A_start")
    fig.suptitle("Step 4c2: the median tail from S1 (solid: the rain within 5 km has ended), from a placebo start with the same delays "
                 "(dashed) and from S0 = t0_trace (dotted)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def figure_split(split):
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.9), sharey=True)
    for ax, kind, title in zip(axes, ("S1", "placebo"), ("split at S1 (the real end of the rain within 5 km)", "split at the placebo start")):
        for part, ls, mk in (("before", "-", "o"), ("after", "--", "s")):
            g = split[(split.split_at == kind) & (split.part == part) & split.used]
            ax.plot(g.median_A_db, g.k_per_min, ls=ls, marker=mk, ms=5, mfc="none", color=INK, lw=1.6,
                    label=f"{part} the split ({'solid, circles' if part == 'before' else 'dashed, squares'})")
        ax.set(xscale="log", yscale="log", xlabel="A (dB), median of the bin", title=title)
        ax.legend(frameon=False, fontsize=8.5)
        style_axes(ax)
    axes[0].set_ylabel("k = -mean dA/dt / median A (1 / min)")
    fig.suptitle("Step 4c2, A3: the drying rate per dB at the same level, before and after the split, alternative-set tails "
                 f"(first {S['a3_window_min']} min from S0)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def figure_gallery(picks, long_a, bounds, ref):
    G = S["gallery"]
    n = len(picks)
    fig = plt.figure(figsize=(4.3 * (n // 2), 9.6))
    gs = fig.add_gridspec(4, n // 2, height_ratios=[3, 1.1, 3, 1.1], hspace=0.5, wspace=0.24, top=0.92, bottom=0.06, left=0.05, right=0.98)
    for j, r in enumerate(picks.itertuples()):
        ax, bx = fig.add_subplot(gs[2 * (j // (n // 2)), j % (n // 2)]), fig.add_subplot(gs[2 * (j // (n // 2)) + 1, j % (n // 2)])
        a = long_a[bounds[r.pos, 0]:bounds[r.pos, 1]]
        lo, hi = max(r.S0_k - G["before_min"], 0), min(r.S0_k + G["after_min"], len(a))
        x = np.arange(lo, hi) - r.S0_k
        ax.plot(x, a[lo:hi], "o", ms=2.1, color="#8a8a8a")
        ax.plot(x, pd.Series(a[lo:hi]).rolling(5, center=True, min_periods=3).median(), "-", color=INK, lw=1.5)
        for i, name in enumerate(KINDS):
            xk = getattr(r, f"{name}_k") - r.S0_k
            for a_ in (ax, bx):
                a_.axvline(xk, color=INK, ls=LINE[name], lw=1.3)
            ax.text(xk, 0.99 - 0.1 * i, f" {name}", transform=ax.get_xaxis_transform(), fontsize=8, va="top",
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.4, alpha=0.85))
        ax.axhline(0, color=INK2, lw=0.6)
        ax.set_title(f"tail {r.tail_id}, event {r.event_id[1:9]}, L = {r.L_km:.1f} km\nA at S1 = {r.S1_A_start_db:.2f} dB", fontsize=8)
        i0 = ref["T5"].searchsorted(r.t0) + lo // 5
        nb = int(np.ceil((hi - lo) / 5)) + 1
        xb = (np.arange(nb) * 5 + (lo // 5) * 5) - r.S0_k
        near, path = np.nan_to_num(ref["near_max"][r.cml_id][i0:i0 + nb]), np.nan_to_num(ref["path_max"][r.cml_id][i0:i0 + nb])
        bx.step(xb[:len(near)], near, where="post", color=INK, lw=1.3, label="radar, largest R_max on the paths within 5 km (solid)")
        bx.step(xb[:len(path)], path, where="post", color=INK, lw=1.1, ls="--", label="radar R_max on the own path (dashed)")
        bx.set_ylim(0, max(1.0, float(np.max(np.r_[near, path, 0.0])) * 1.15))
        for a_ in (ax, bx):
            a_.set_xlim(-G["before_min"], G["after_min"])
            style_axes(a_)
        if j % (n // 2) == 0:
            ax.set_ylabel("A (dB)")
            bx.set_ylabel("mm/h")
        if j == 0:
            bx.legend(frameon=False, fontsize=6.5, loc="upper right")
        if j // (n // 2) == 1:
            bx.set_xlabel("minutes since S0 = t0_trace")
    fig.suptitle(f"Step 4c2: {n} pairs of the main S1 set (A at S1 >= {MAIN} dB; seed {G['seed']}). Dots: A, 1 min; line: 5 min median. Vertical "
                 "lines, labelled: S1 solid, placebo start dashed, S0 dotted. Below: radar nearby and on the path.", fontsize=9.5)
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    CACHE.mkdir(parents=True, exist_ok=True)
    clock = time.time()
    rng = np.random.default_rng(S["seed"])

    # ------------------------------------------------------------------------------------------- data (TRAIN only)
    sel = pd.read_csv(P.TAILS_SELECTED, parse_dates=["t0", "t_use_end_lenient"], low_memory=False,
                      usecols=["tail_id", "split_event", "t0", "t_use_end_lenient", "truncated_by_lenient"]).query("split_event == 'train'")   # VAL rows dropped here
    pool = pd.read_csv(P.TAIL_STARTS, parse_dates=["t0", "S0_start", "S1_start"])[
        ["tail_id", "sublink_id", "cml_id", "event_id", "t0", "L_km", "in_v1_main", "in_v1_alt", "S0_start", "S1_start"]]
    assert set(pool.tail_id) <= set(sel.tail_id), "tails_starts.csv must hold TRAIN pairs only"
    pool = pool.merge(sel[["tail_id", "t_use_end_lenient", "truncated_by_lenient"]], on="tail_id").sort_values("tail_id").reset_index(drop=True)
    long = pq.read_table(P.TAILS_1MIN, columns=["tail_id", "minutes_since_t0", "A_db"], filters=[("tail_id", "in", pool.tail_id.tolist())]).to_pandas()
    long = long.sort_values(["tail_id", "minutes_since_t0"], kind="stable").reset_index(drop=True)
    long_a = long.A_db.to_numpy("float64")
    bounds = np.column_stack([np.searchsorted(long.tail_id.to_numpy(), pool.tail_id.to_numpy(), side=s_) for s_ in ("left", "right")])
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    ref = reference_traces(d, pool.cml_id.unique(), B2B["gauge_max_dist_m"], A3C["s1_near_km"], A3C["s1_city_gauge_km"], A3C["s1_smhi_gauge_km"])
    # the placebo delay of a pair: the S1 - S0 delay of ANOTHER pair, drawn at random (seed fixed)
    delay = ((pool.S1_start - pool.S0_start) / pd.Timedelta("1min")).round().astype(int).to_numpy()
    other = np.random.default_rng(S["placebo_seed"]).integers(0, len(pool) - 1, len(pool))
    other = other + (other >= np.arange(len(pool)))        # never the pair itself
    delay_p = delay[other]
    print(f"TRAIN only: {len(pool)} pairs, {pool.event_id.nunique()} events, {pool.cml_id.nunique()} links; loaded in {time.time() - clock:.0f} s", flush=True)

    plan_file, rows_file, win_file = CACHE / "plan.json", CACHE / "pairs.csv", CACHE / "windows.npz"
    if plan_file.exists() and rows_file.exists() and win_file.exists():
        plan, rows = json.loads(plan_file.read_text()), pd.read_csv(rows_file)
        win = dict(np.load(win_file))
        print(f"pair pass from the cache ({rel(rows_file)}); plan: cut {plan['cuts'] or 'nothing'}")
    else:
        n_t = S["timing_pairs"]
        t0 = time.time()
        first_rows, first_win = pair_pass(pool.head(n_t), long_a, bounds[:n_t], ref, delay_p[:n_t])
        per = (time.time() - t0) / n_t
        proj = per * len(pool) + 120                       # the statistics and bootstraps: about 2 min
        plan = {"seconds_per_pair": round(per, 4), "projected_s": round(proj), "cuts": []}
        print(f"timing: {per * 1e3:.1f} ms per pair on {n_t} pairs -> the whole step is projected at {proj / 60:.1f} min "
              f"(budget {S['time_budget_min']} min)", flush=True)
        for cut in ("threshold 0.25", "A3"):
            if proj <= S["time_budget_min"] * 60:
                break
            plan["cuts"].append(cut)
            proj *= 0.7
        rest_rows, rest_win = pair_pass(pool.iloc[n_t:], long_a, bounds[n_t:], ref, delay_p[n_t:])
        rows = pd.concat([first_rows, rest_rows], ignore_index=True)
        win = {k: np.vstack([first_win[k], rest_win[k]]) for k in KINDS}
        rows.to_csv(rows_file, index=False)
        np.savez_compressed(win_file, **win)
        plan_file.write_text(json.dumps(plan))
    pool = pool.merge(rows, on="tail_id", validate="1:1")
    pool["placebo_delay_min"], pool["S1_delay_min"] = delay_p, delay
    assert (pool.S1_k - pool.S0_k == delay).all()
    thresholds = [t for t in S["thresholds_db"] if not (t == 0.25 and "threshold 0.25" in plan["cuts"])]

    # ------------------------------------------------------------------------------------------- A1 + A2
    stats, rate_m = {}, {}
    for name in KINDS:
        stats[name] = tail_stats(win[name], pool[f"{name}_A_start_db"].to_numpy(float), pool[f"{name}_usable_min"].to_numpy(float))
        rate_m[name] = rate_arrays(win[name][:, :WIN].astype(float))
        pool[f"{name}_t_half_min"], pool[f"{name}_halved"] = stats[name]["t_half"], stats[name]["halved"]
    keep = {(name, thr): ((pool[f"{name}_A_start_db"] >= thr) & (pool[f"{name}_usable_min"] >= S["min_usable_min"])).to_numpy()
            for name in KINDS for thr in thresholds}
    set_rows, rate_rows, curve_rows, diff_rows = [], [], [], []
    for thr in thresholds:
        for name in KINDS:
            row, tab = describe(name, thr, keep[name, thr], stats[name], win[name], pool, *rate_m[name], rng)
            row["median_delay_from_S0_min"] = float(np.median((pool[f"{name}_k"] - pool.S0_k).to_numpy()[keep[name, thr]]))
            row["median_usable_min"] = float(np.median(pool[f"{name}_usable_min"].to_numpy()[keep[name, thr]]))
            set_rows.append(row), rate_rows.append(tab)
            Wk, a0 = win[name][keep[name, thr]], pool[f"{name}_A_start_db"].to_numpy(float)[keep[name, thr]]
            with np.errstate(all="ignore"):
                curve_rows.append(pd.DataFrame({"threshold_db": thr, "start": name, "minute": np.arange(WIN),
                                                "median_A_over_A_start": np.nanmedian(Wk[:, :WIN] / a0[:, None], axis=0),
                                                "median_A_db": np.nanmedian(Wk[:, :WIN], axis=0), "tails": np.isfinite(Wk[:, :WIN]).sum(axis=0)}))
        # information: the SAME pairs as the S1 set, looked at from S0 (those with A at S0 >= the threshold)
        same = keep["S1", thr] & (pool.S0_A_start_db >= thr).to_numpy() & (pool.S0_usable_min >= S["min_usable_min"]).to_numpy()
        row, _ = describe("S0", thr, same, stats["S0"], win["S0"], pool, *rate_m["S0"], rng, with_boot=False)
        row["start"] = "S0, only the pairs of the S1 set"
        set_rows.append(row)
        # the difference S1 - placebo with an event bootstrap (whole events; both sets come from the same draw)
        i1, ip = np.nonzero(keep["S1", thr])[0], np.nonzero(keep["placebo", thr])[0]
        ev = pool.event_id.to_numpy()
        events = sorted(set(ev[i1]) | set(ev[ip]))
        by1, byp = {e: i1[ev[i1] == e] for e in events}, {e: ip[ev[ip] == e] for e in events}
        r30 = f"ratio_{30}"

        def both(r1, rp):
            return (km_half(stats["S1"]["t_half"][r1], stats["S1"]["halved"][r1])[0] - km_half(stats["placebo"]["t_half"][rp], stats["placebo"]["halved"][rp])[0],
                    float(np.nanmedian(stats["S1"][r30][r1]) - np.nanmedian(stats["placebo"][r30][rp])))

        obs = both(i1, ip)
        draws = []
        for _ in range(S["boot_reps"]):
            pick = rng.integers(0, len(events), len(events))
            r1, rp = np.concatenate([by1[events[q]] for q in pick]), np.concatenate([byp[events[q]] for q in pick])
            draws.append(both(r1, rp) if len(r1) and len(rp) else (np.nan, np.nan))
        lo, hi = np.nanpercentile(np.asarray(draws, float), S["interval"], axis=0)
        diff_rows.append({"threshold_db": thr, "t_half_S1_minus_placebo_min": obs[0], "t_half_lo": float(lo[0]), "t_half_hi": float(hi[0]),
                          "ratio30_S1_minus_placebo": round(obs[1], 3), "ratio30_lo": round(float(lo[1]), 3), "ratio30_hi": round(float(hi[1]), 3),
                          "pairs_in_both_sets": int((keep["S1", thr] & keep["placebo", thr]).sum())})
    sets = pd.DataFrame(set_rows)
    diff = pd.DataFrame(diff_rows)
    # added after the first run (information, no rule): many placebo starts fall after the pair's own S1, so the
    # placebo set mixes "rain still within 5 km" and "rain gone". The placebo tails split by that:
    info_rows = []
    for thr in thresholds:
        for lab, m_ in (("placebo start BEFORE the pair's S1 (rain still within 5 km)", (pool.placebo_k < pool.S1_k).to_numpy()),
                        ("placebo start at or after the pair's S1", (pool.placebo_k >= pool.S1_k).to_numpy())):
            row, _ = describe("placebo", thr, keep["placebo", thr] & m_, stats["placebo"], win["placebo"], pool, *rate_m["placebo"], rng, with_boot=False)
            info_rows.append({"threshold_db": thr, "placebo tails": lab, **{k: row[k] for k in ("tails", "events", "median_A_start_db",
                              "A_over_A_start_10min", "A_over_A_start_30min", "A_over_A_start_60min", "t_half_min", "share_half_within_60min")},
                              "median_delay_from_S0_min": float(np.median(pool.placebo_delay_min.to_numpy()[keep["placebo", thr] & m_]))})
    info = pd.DataFrame(info_rows)
    save(info, P.AFTER_RAIN_PLACEBO_SPLIT, index=False)
    curves = pd.concat(curve_rows, ignore_index=True).round(4)
    save(sets, P.AFTER_RAIN_SETS, index=False)
    save(diff, P.AFTER_RAIN_DIFF, index=False)
    save(curves, P.AFTER_RAIN_CURVES, index=False)
    save(pd.concat(rate_rows, ignore_index=True).round(5), P.AFTER_RAIN_RATE, index=False)
    out_cols = ["tail_id", "sublink_id", "cml_id", "event_id", "t0", "L_km", "in_v1_main", "in_v1_alt", "S1_delay_min", "placebo_delay_min", "S1_next_rain_nearby_k"]
    for name in KINDS:
        out_cols += [f"{name}_k", f"{name}_end_k", f"{name}_usable_min", f"{name}_truncated_by", f"{name}_A_start_db", f"{name}_t_half_min", f"{name}_halved"]
    save(pool[out_cols].round(4), P.AFTER_RAIN_PAIRS, index=False)

    # check: the fast k(A) of this script is the raw_rate of step 4b (same smoothing, same bins)
    i0 = np.nonzero(keep["S0", MAIN])[0][:400]
    chk = pd.DataFrame({"tail_id": np.repeat(pool.tail_id.to_numpy()[i0], WIN), "minutes_since_t0_trace": np.tile(np.arange(WIN), len(i0)),
                        "A_db": win["S0"][i0, :WIN].astype(float).ravel()})
    ref_rate = raw_rate(chk, RATE["edges_db"], RATE["smooth_min"]).set_index("A_bin_db").k_from_mean_per_min
    mine = rate_table(rate_m["S0"][0][i0], rate_m["S0"][1][i0])[0].set_index("A_bin_db").k_per_min
    rate_check = float(np.nanmax(np.abs(ref_rate.to_numpy() - mine.reindex(ref_rate.index).to_numpy())))
    assert rate_check < 2e-3, f"k(A) here differs from raw_rate by {rate_check}"

    # ------------------------------------------------------------------------------------------- the rule
    g = sets.set_index(["threshold_db", "start"])
    t1, tp, ts0 = g.loc[(MAIN, "S1"), "t_half_min"], g.loc[(MAIN, "placebo"), "t_half_min"], g.loc[(MAIN, "S0"), "t_half_min"]
    dm = diff.set_index("threshold_db").loc[MAIN]
    slow = bool(t1 >= RULE["slow_ratio"] * ts0 and dm.t_half_lo > 0)
    same_speed = bool(abs(t1 - tp) <= RULE["same_tol"] * tp and dm.t_half_lo <= 0 <= dm.t_half_hi)
    verdict = "drying after the rain ends is SLOW" if slow else "SAME speed" if same_speed else "unclear"
    print(f"main threshold {MAIN} dB: t_half S1 {t1:.0f}, placebo {tp:.0f}, S0 {ts0:.0f} min; S1 - placebo {dm.t_half_S1_minus_placebo_min:.0f} "
          f"({dm.t_half_lo:.0f}..{dm.t_half_hi:.0f}) -> {verdict}", flush=True)

    # ------------------------------------------------------------------------------------------- A3
    split = pd.DataFrame()
    if "A3" not in plan["cuts"]:
        alt = pool[pool.in_v1_alt].reset_index(drop=True)
        pos = pool.index[pool.in_v1_alt].to_numpy()
        n3 = S["a3_window_min"]
        M3 = np.full((len(alt), n3), np.nan)
        for j, r in enumerate(alt.itertuples()):
            seg = long_a[bounds[pos[j], 0] + r.S0_k:bounds[pos[j], 0] + min(r.S0_end_k, r.S0_k + n3)]
            M3[j, :len(seg)] = seg
        sm3, d3 = rate_arrays(M3)
        m = np.arange(n3)[None, :]
        d1, dpl, back = (alt.S1_k - alt.S0_k).to_numpy()[:, None], alt.placebo_delay_min.to_numpy()[:, None], (alt.S1_next_rain_nearby_k - alt.S0_k).to_numpy()[:, None]
        masks = {("S1", "before"): m < d1, ("S1", "after"): (m >= d1) & (m < back), ("placebo", "before"): m < dpl, ("placebo", "after"): m >= dpl}
        parts = []
        for (kind, part), mk in masks.items():
            tab, _ = rate_table(sm3, d3, mk)
            parts.append(tab.assign(split_at=kind, part=part, tails=int((np.isfinite(sm3) & mk).any(axis=1).sum())))
        split = pd.concat(parts, ignore_index=True)
        wide = split.pivot(index=["split_at", "A_bin_db"], columns="part", values=["k_per_min", "tail_minutes", "used"])
        ratio = (wide["k_per_min"]["after"] / wide["k_per_min"]["before"]).where(wide["used"]["after"].astype(bool) & wide["used"]["before"].astype(bool))
        split["k_after_over_before"] = [ratio.get((a_, b_), np.nan) for a_, b_ in zip(split.split_at, split.A_bin_db)]
        save(split.round(5), P.AFTER_RAIN_SPLIT, index=False)
        a3_show = pd.DataFrame({"A_bin_db": [f"({EDGES[j]:g}, {EDGES[j + 1]:g}]" for j in range(len(EDGES) - 1)]})
        for kind in ("S1", "placebo"):
            for part in ("before", "after"):
                q = split[(split.split_at == kind) & (split.part == part)].set_index("A_bin_db")
                a3_show[f"{kind}: minutes {part}"] = a3_show.A_bin_db.map(q.tail_minutes)
                a3_show[f"{kind}: k {part}"] = a3_show.A_bin_db.map(q.k_per_min.where(q.used)).round(4)
            a3_show[f"{kind}: k after / before"] = a3_show.A_bin_db.map(ratio[kind]).round(2)
        rr = a3_show.set_index("A_bin_db")
        both_ok = rr[["S1: k after / before", "placebo: k after / before"]].dropna()
        a3_n = len(alt)

    # ------------------------------------------------------------------------------------------- figures
    fig = figure_curves(curves, sets)
    save(fig, P.AFTER_RAIN_FIG_CURVES, dpi=110)
    plt.close(fig)
    if len(split):
        fig = figure_split(split)
        save(fig, P.AFTER_RAIN_FIG_SPLIT, dpi=110)
        plt.close(fig)
    picks = pool.assign(pos=np.arange(len(pool)))[keep["S1", MAIN]].sample(S["gallery"]["pairs"], random_state=S["gallery"]["seed"]).sort_values("tail_id")
    fig = figure_gallery(picks, long_a, bounds, ref)
    save(fig, P.AFTER_RAIN_FIG_GALLERY, dpi=105)
    plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    show_cols = ["threshold_db", "start", "tails", "events", "links", "median_A_start_db", "A_over_A_start_10min", "A_over_A_start_30min",
                 "A_over_A_start_60min", "tails_with_60min", "t_half_min", "share_half_within_60min", "share_censored", "k_slope", "k_slope_lo", "k_slope_hi",
                 "bins_in_slope", "median_delay_from_S0_min", "median_usable_min"]
    a1_show = sets[sets.start == "S1"][show_cols]
    a2_show = sets[sets.start != "S1"][show_cols]
    ratio_tab = pd.DataFrame([{"threshold_db": thr, "t_half S1": g.loc[(thr, "S1"), "t_half_min"], "t_half placebo": g.loc[(thr, "placebo"), "t_half_min"],
                               "t_half S0": g.loc[(thr, "S0"), "t_half_min"],
                               "S1 / S0": round(g.loc[(thr, "S1"), "t_half_min"] / g.loc[(thr, "S0"), "t_half_min"], 2),
                               "S1 / placebo": round(g.loc[(thr, "S1"), "t_half_min"] / g.loc[(thr, "placebo"), "t_half_min"], 2)} for thr in thresholds])
    gm = lambda name, colname: g.loc[(MAIN, name), colname]               # noqa: E731
    im_ = info[info.threshold_db == MAIN].reset_index(drop=True)
    ib, ia = im_.iloc[0], im_.iloc[1]                      # placebo starts before / at or after the pair's S1
    dh = diff.set_index("threshold_db").loc[thresholds[-1]]
    after_share = float((pool.placebo_k >= pool.S1_k).to_numpy()[keep["placebo", MAIN]].mean())
    report = f"""# wet_antenna step 4c2: does drying really start when the rain nearby ends? (TRAIN pairs, no law fit)

Written by `research/wet_antenna/check_drying_after_rain.py`. Settings and the rule: `wet_antenna.drying_after_rain` in
`research/config.yaml`, written before the run. Why this step: `docs/ROADMAP.md`, "Decisions of the learners after 4c".

## 1. What was run

- Data: the {len(pool)} TRAIN pairs that pass the v1 exclusions ({pool.event_id.nunique()} events, {pool.cml_id.nunique()} links), with S0 (t0_trace) and S1
  (no rain within {A3C['s1_near_km']:g} km, + {A3C['s1_buffer_min']} min) from `{rel(P.TAIL_STARTS)}` and the 1 min A(t) of v1. VAL is never read. No law
  is fitted. Nothing earlier is changed.
- Three starts per pair, each with the v1 end rules counted from that start (next rain trace on the path, rise cut
  v2, the step 2b ends): **S1**, which ALSO ends at the next rain within {A3C['s1_near_km']:g} km, so its window is truly "no rain
  nearby"; a **placebo** start = S0 + the S1 - S0 delay of another pair drawn at random (seed {S['placebo_seed']}), without the
  "rain nearby" condition; and **S0**. A tail is kept at a threshold if A at the start (mean of the first 5 min) is
  at least the threshold and the tail is at least {S['min_usable_min']} min long.
- t_half of a set is the Kaplan-Meier median of the time to half: a tail that ends before it has halved is censored
  at its end (S1 windows are shorter, because they also end when rain comes back nearby). Capped at {WIN} min.
- Timing: {plan['seconds_per_pair'] * 1e3:.1f} ms per pair on the first {S['timing_pairs']} pairs; projection {plan['projected_s'] / 60:.1f} min against a budget of
  {S['time_budget_min']} min. Cut: {', '.join(plan['cuts']) or 'nothing'}. Run time of this run: {time.time() - clock:.0f} s (the pair pass is cached).
- Check: the k(A) computed here equals `raw_rate` of step 4b on 400 tails (largest difference {rate_check:.1e} per min).

## 2. Part 0: the decisions

The four decisions of the learners after 4c were recorded before the analysis: `docs/ROADMAP.md` (section "Decisions
of the learners after 4c", status rows 4c2 / 3c / 7 / 4d, the order line), `research/wet_antenna/METHOD.md` (under
"Step 4c" and at the top of the open issues) and `research/wet_antenna/README.md` (status).

## 3. A1: tails from S1, lower thresholds

{md_table(a1_show)}

(`{rel(P.AFTER_RAIN_SETS)}`; k(A) per level: `{rel(P.AFTER_RAIN_RATE)}`; per pair: `{rel(P.AFTER_RAIN_PAIRS)}`.)
`share_censored` = tails that ended before halving and before {WIN} min. The k(A) slope is that of log k vs log A from
{EDGES[0]:g} dB up (bins with >= {RATE['min_minutes']} tail-minutes; interval from {RATE['boot_reps']} event bootstraps).

## 4. A2: the placebo start, S0, and the difference

{md_table(a2_show)}

"S0, only the pairs of the S1 set" (information): the pairs of the S1 set, looked at from S0, where A at S0 is also
above the threshold.

The difference S1 - placebo ({S['boot_reps']} event bootstraps, 90 % interval; `{rel(P.AFTER_RAIN_DIFF)}`):

{md_table(diff)}

{md_table(ratio_tab)}

Added after the first run (information, no rule; `{rel(P.AFTER_RAIN_PLACEBO_SPLIT)}`): the placebo set mixes starts
that fall while it still rains within {A3C['s1_near_km']:g} km and starts after the pair's own S1. Split by that:

{md_table(info)}

![median curves](../figures/{P.AFTER_RAIN_FIG_CURVES.name})

![gallery](../figures/{P.AFTER_RAIN_FIG_GALLERY.name})

## 5. A3: the same level before and after S1
"""
    if len(split):
        report += f"""
The {a3_n} alternative-set tails of v1 (A_start >= {S2SEL['alt_a_start_db']} dB at S0), first {S['a3_window_min']} min from S0. The minutes are split into
"before S1" (rain within {A3C['s1_near_km']:g} km) and "from S1 until rain comes back within {A3C['s1_near_km']:g} km"; k = -mean dA/dt / median A per
level in each part (bins with >= {RATE['min_minutes']} minutes). The same split at the placebo start.

{md_table(a3_show)}

(`{rel(P.AFTER_RAIN_SPLIT)}`.) **Caution:** tails that fall fast reach the low levels before S1, so "before" holds the
fast ones at each level and "after" the slow ones, whatever the cause. The placebo split has the same bias, so the
S1 ratio is to be read against the placebo ratio, not against 1. Bins with both ratios: {len(both_ok)}; S1 ratio
{', '.join(f'{v:.2f}' for v in both_ok.iloc[:, 0])}; placebo ratio {', '.join(f'{v:.2f}' for v in both_ok.iloc[:, 1])}.

![A3](../figures/{P.AFTER_RAIN_FIG_SPLIT.name})
"""
    else:
        report += "\nA3 was cut (time budget).\n"
    report += f"""
## 6. Verdict by the rule (threshold {MAIN} dB)

- t_half: S1 **{t1:.0f} min** ({int(gm('S1', 'tails'))} tails), placebo **{tp:.0f} min** ({int(gm('placebo', 'tails'))} tails), S0 **{ts0:.0f} min** ({int(gm('S0', 'tails'))} tails).
- S1 / S0 = {t1 / ts0:.2f} (rule for SLOW: >= {RULE['slow_ratio']:g}); S1 - placebo = {dm.t_half_S1_minus_placebo_min:+.0f} min, 90 % interval
  {dm.t_half_lo:+.0f}..{dm.t_half_hi:+.0f} (SLOW needs it above 0; SAME needs 0 inside and S1 within +-{RULE['same_tol']:.0%} of the placebo:
  S1 / placebo = {t1 / tp:.2f}).
- A / A_start at 30 min: S1 {gm('S1', 'A_over_A_start_30min')}, placebo {gm('placebo', 'A_over_A_start_30min')}, S0 {gm('S0', 'A_over_A_start_30min')}; S1 - placebo
  {dm.ratio30_S1_minus_placebo:+.3f} ({dm.ratio30_lo:+.3f}..{dm.ratio30_hi:+.3f}).

**Verdict: {verdict.upper()}.**

The reading written with the rule: SLOW = the fast early drop belongs to the end of the rain (path or nearby), not
to the antenna, and tau of about 10 min is not a drying time; SAME = the fast drop is antenna drying.

## 7. In plain words

- **From S1 on, A falls slowly: it halves in {t1:.0f} min, against {ts0:.0f} min from t0_trace.** That is the first half of the
  SLOW rule (ratio {t1 / ts0:.2f}). But a tail started just as late at a time that has nothing to do with the rain halves
  in {tp:.0f} min: the difference is {dm.t_half_S1_minus_placebo_min:+.0f} min ({dm.t_half_lo:+.0f}..{dm.t_half_hi:+.0f}). So most of the slowness comes from STARTING LATE, not from
  the rain within {A3C['s1_near_km']:g} km having ended. By the rule: SAME speed.
- **The verdict rests on a diluted contrast.** {after_share:.0%} of the placebo starts lie at or after the pair's own S1, so they
  are "no rain nearby" too. Split (added after the first run, no rule): placebo tails that start while it still
  rains within {A3C['s1_near_km']:g} km halve in {ib.t_half_min:.0f} min ({int(ib.tails)} tails, median delay {ib.median_delay_from_S0_min:.0f} min), those that start at or after their S1
  in {ia.t_half_min:.0f} min ({int(ia.tails)} tails, median delay {ia.median_delay_from_S0_min:.0f} min); the S1 tails in {t1:.0f} min. That gap of about {ia.t_half_min - ib.t_half_min:.0f} min goes in the
  direction of SLOW, but it is mixed up with the delay itself (a start before the pair's S1 is an earlier start).
  So: the end of the rain nearby accounts for at most a small part ({ib.t_half_min:.0f} -> {ia.t_half_min:.0f} min) of the step from {ts0:.0f} to {t1:.0f} min;
  most of it is the time since t0_trace.
- **What that does and does not show.** The end of the rain within {A3C['s1_near_km']:g} km is not a clear event in the signal, so
  the fast early drop is mostly not "nearby rain ending". It does NOT show that the fast drop is the antenna: the
  drop starts at t0_trace, which is the last rain trace ON THE PATH, and rain on the path that the radar no longer
  sees would end there too. The reading "SAME = antenna drying" written with the rule is stronger than this test;
  Q2 is narrowed (mostly not the rain nearby), not closed.
- **A tail has two phases in TIME.** About the first 15-20 min after t0_trace are fast; what is left after that
  (0.3-0.9 dB) halves in about {t1:.0f} min at every threshold ({', '.join(f"{g.loc[(t_, 'S1'), 't_half_min']:.0f}" for t_ in thresholds)} min), which would be tau of about {t1 / np.log(2):.0f} min IF it
  were exponential. At the same level of A the rate per dB is {1 / both_ok.iloc[:, 0].max():.0f} to {1 / both_ok.iloc[:, 0].min():.0f} times lower after the split than before
  it (A3), and nearly the same for the placebo split. Either the tails differ a lot in speed (the fast ones are gone
  early: the caution of A3), or the rate depends on the time since the rain and not on the level alone. Step 4b
  pointed the same way (the stall belongs to the event, not to the level).
- **A hint at the higher levels, not at the rule's threshold.** At {thresholds[-1]:g} dB the S1 tails are slower than the placebo
  tails (t_half {dh.t_half_S1_minus_placebo_min:+.0f} min, {dh.t_half_lo:+.0f}..{dh.t_half_hi:+.0f}; A / A_start at 30 min {dh.ratio30_S1_minus_placebo:+.3f}, {dh.ratio30_lo:+.3f}..{dh.ratio30_hi:+.3f}), and above
  0.7 dB the k ratio of the S1 split is below that of the placebo split ({', '.join(f"{a_:.2f} vs {b_:.2f}" for a_, b_ in both_ok.iloc[3:].itertuples(index=False))}). So for
  the wetter tails some of the early speed may go with rain that is still around.
- **For the drying time.** "tau about 10 min" describes the first phase only. It is not shown to be a drying time, and
  it is not shown to be rain. The second phase (t_half about {t1:.0f} min) is there with no rain within {A3C['s1_near_km']:g} km on any
  reference, so it is the better candidate for antenna water (or a baseline that returns slowly; step 2b's
  `undershoot` and `level_shift` flags).
- **For Step 7 (whole events on very short links).** The model needs room for two time scales (about 10-15 min and
  about 40-60 min) and must not be forced into one exponential. The radar rain as input lets it put the fast phase
  on the path if it scales with L; very short links are where the path part is smallest, so if the fast phase is
  still there at L below about 0.6 km it is the antenna. That is the test this step could not make.

## 8. Unclear / possibly wrong

1. **The placebo is a weak contrast.** {int(dm.pairs_in_both_sets)} pairs are in both the S1 and the placebo set at {MAIN} dB, and {after_share:.0%} of the
   placebo starts of that set lie at or after the pair's own S1 (the delays have the same distribution). The split
   in section 4 is the cleaner comparison, but it was added after the first run and has no rule.
2. **The sets are still selected on the outcome**: a tail is in the S1 (or placebo) set only if it is still above the
   threshold that late. That is the same for both, which is the point of the placebo, but it means neither t_half
   is "the" drying time of all antennas.
3. **Censoring.** S1 windows end when rain comes back within {A3C['s1_near_km']:g} km (median usable {gm('S1', 'median_usable_min'):.0f} min against {gm('placebo', 'median_usable_min'):.0f}); the
   Kaplan-Meier median assumes that this end is unrelated to how the tail would have gone on.
4. **The S0 sets here are not tail set v1**: {int(gm('S0', 'tails'))} tails at {MAIN} dB with a 30 min minimum, against 429 at 1 dB and 60 min; their
   t_half is {ts0:.0f} min (v1 main set: 14 min).
5. **"No rain within {A3C['s1_near_km']:g} km" is as good as the references**: radar on link paths (2 km, 5 min) and gauges. Drizzle that the
   radar misses is "no rain" here. The {A3C['s1_buffer_min']} min buffer is a choice.
6. **A3 mixes tails of different speed** (the caution in section 5); its ratios are descriptive.
7. **Low levels and the 0.3 dB RSL steps**: at 0.25-0.4 dB a tail is one quantization step above its baseline, so
   "half" is crossed by noise. The 0.25 dB row is the least reliable.
8. **t_half -> tau** uses tau = t_half / ln 2, true only for an exponential; nothing here shows the second phase is one.
9. **"No rain nearby" windows still hold bumps that look like rain** (gallery, observation only): three of the eight
   pairs shown, all in E20150617_0705, rise by about 1 dB for 10-15 min at the same clock time while the radar is
   valid and shows nothing on any path within {A3C['s1_near_km']:g} km. That event gives {int((pool.event_id.to_numpy()[keep['S1', MAIN]] == 'E20150617_0705').sum())} of the {int(gm('S1', 'tails'))} S1 tails at {MAIN} dB. The rise
   cut v2 ends such tails only from 0.6 dB for 10 min on.
"""
    save(report, P.AFTER_RAIN_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {time.time() - clock:.0f} s")


if __name__ == "__main__":
    main()
