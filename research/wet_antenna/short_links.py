"""wet_antenna step 7a: very short links. Is there enough data for Step 7, and do they show the fast first phase?

Step 4c2 found two phases in a drying tail: fast for 15-20 min after t0_trace, then slow (half-life about 40 min).
Open (Q2 of docs/ROADMAP.md): is the fast phase the antenna drying, or rain on the path that the radar no longer sees?
Rain on the path attenuates in proportion to the path length L; water on the antennas does not. On a very short link
there is almost no path rain, so if the fast phase is just as big there, it is the antenna. Three parts, NO law fit:
  part 1  inventory per length class (< 0.6, 0.6-1, 1-2, 2-4, >= 4 km): links, frequency bands, rain on the link in
          the TRAIN events, pairs and tails; and the rule "is there enough for Step 7 below 0.6 km, or below 1 km?";
  part 2  the two phases by length: per tail the fast-phase drop (A_start minus A at 15-20 min) and the slow-phase
          ratios; SHORT (< 1 km) against LONG (>= 2 km) tails matched on A_start (and on the event), with an event
          bootstrap, next to the difference that path rain alone would give (ITU-R P.838-3 on the radar rain before
          the tail). The rule that reads the result is in config `wet_antenna.short_links.part2.rule`;
  part 3  a first look at whole events on the links below 1 km: the peak of A against the peak of the rain, the
          time from the end of the rain to the peak and to half of it, A rising before the radar shows rain, the
          drying per frequency band; and a gallery.
Every rule was written in the config before the run. TRAIN links and TRAIN events only: tails_starts.csv (step 4c)
holds only TRAIN pairs, VAL rows of the other files are dropped on reading, and minutes that fall inside an event
that is not a TRAIN event are blanked in part 3. Nothing is cached (the whole step takes about a minute).
Run from the repository root:
    python -m research.wet_antenna.short_links
"""
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
from research.wet_antenna import paths as P
from research.wet_antenna.tails import classify, first_run, itu_p838, minutes, spearman

CFG = load_config()
W = CFG["wet_antenna"]
S = W["short_links"]
P1, P2, P3, RULE, RISE = S["part1"], S["part2"], S["part3"], S["part2"]["rule"], S["part3"]["rise"]
CL = W["diagnose_stalled_tails"]["classes"]                # the class rules of step 4b
CLASSES = ["decays", "plateau-then-decays", "rises", "flat / slow"]
EDGES = [0.0] + list(S["length_edges_km"]) + [np.inf]
LENGTHS = [f"< {EDGES[1]:g} km"] + [f"{a:g}-{b:g} km" for a, b in zip(EDGES[1:-2], EDGES[2:-1])] + [f">= {EDGES[-2]:g} km"]
BANDS = list(S["bands_ghz"])
FOLLOW, PAD = P2["follow_min"], P2["value_window_min"]
A_BINS = np.asarray(P2["a_start_bins_db"], float)
MATCHINGS = ["event and A_start bin", "A_start bin"]
QUANTITIES = {"rel_fast_drop": "relative fast drop", "fast_drop_db": "fast drop (dB)", "exp_path_rel": "expected path rain / A_start",
              "exp_path_db": "expected path rain (dB)", "slow_ratio": "slow ratio (50-60 / 15-20 min)", "late_ratio": "late ratio (170-180 / 15-20 min)"}
LINE = ["-", "--", "-.", ":", (0, (5, 1, 1, 1, 1, 1))]     # one line style per length class, in every figure
MARK = ["o", "s", "^", "D", "v"]
MIN = pd.Timedelta("1min")


def length_class(km):
    """The length class of each path length (km), as an ordered category."""
    return pd.cut(km, EDGES, right=False, labels=LENGTHS)


def band_of(f_ghz):
    """The frequency band of each sub-link frequency (GHz): config `bands_ghz`."""
    out = pd.Series("other", index=f_ghz.index)
    for name, (lo, hi) in S["bands_ghz"].items():
        out[(f_ghz >= lo) & (f_ghz < hi)] = name
    return out


def counts_text(s, order=None):
    """'name: n, name: n' of the values of a Series (for one table cell)."""
    c = s.value_counts()
    return ", ".join(f"{k}: {c[k]}" for k in (order or sorted(c.index)) if k in c.index)


# =================================================================================================== part 1: inventory
def inventory(d, pool, kept, state, T5, st_col, R):
    """One row per length class (and the two cuts < 0.6 km, < 1 km): the links, the rain on them in the TRAIN events,
    the pairs and the tails."""
    links = d.links()
    links = links[links.split_link == "train"].copy()
    subs = d.sublinks()
    subs = subs[subs.split_link == "train"].copy()
    base = pd.read_csv(P.BASELINE_MODULE / "outputs" / "products" / "baselines.csv",
                       usecols=["event_id", "split_event", "sublink_id", "cml_id", "missing", "link_has_wet_slot", "use_in_waa_main"]).query("split_event == 'train'")
    ev = d.events()
    ev = ev[ev.split_event == "train"]
    in_train = np.zeros(len(T5), bool)                     # the 5 min bins inside a TRAIN event
    wet_bins = np.zeros((len(ev), state.shape[1]), int)    # (TRAIN event, link): wet bins of the link in the event
    for i, e in enumerate(ev.itertuples()):
        i0, i1 = T5.searchsorted(e.start), T5.searchsorted(e.end)
        in_train[i0:i1] = True
        wet_bins[i] = (state[i0:i1] == 1).sum(axis=0)
    links["lc"], subs["lc"] = length_class(links.length_m / 1000), length_class(subs.length_m / 1000)
    subs["band"] = band_of(subs.frequency_ghz)
    use_any = base.groupby("sublink_id").use_in_waa_main.any()
    subs["use_in_waa_main"] = subs.sublink_id.map(use_any).fillna(False).astype(bool)
    base = base.merge(links[["cml_id", "lc"]], on="cml_id")
    pool = pool.assign(kept=kept)
    groups = [(name, links.lc == name) for name in LENGTHS] + [(f"CUT < {c:g} km", links.length_m / 1000 < c) for c in (P2["very_short_max_km"], P2["short_max_km"])]
    rows = []
    for name, lm in groups:
        ids = links.cml_id[lm].tolist()
        cols = [st_col[c] for c in ids]
        s_, b_, p_ = subs[subs.cml_id.isin(ids)], base[base.cml_id.isin(ids)], pool[pool.cml_id.isin(ids)]
        wb = wet_bins[:, cols]
        wet_r = R[:, cols][(state[:, cols] == 1) & in_train[:, None]]
        wet_r = wet_r[np.isfinite(wet_r)]
        rows.append({"length": name, "links": len(ids), "sub_links": len(s_), "sub_links_good": int((s_.status == "good").sum()),
                     "sub_links_use_in_waa_main": int(s_.use_in_waa_main.sum()),
                     "bands (sub-links)": counts_text(s_.band, BANDS), "polarization (sub-links)": counts_text(s_.polarization),
                     "train_events_with_rain": int((wb.sum(axis=1) > 0).sum()), "link_events_with_rain": int((wb > 0).sum()),
                     "events_with_rain_per_link_median": float(np.median((wb > 0).sum(axis=0))) if len(ids) else np.nan,
                     "wet_minutes_total": int(wb.sum() * 5), "wet_minutes_per_link_median": float(np.median(wb.sum(axis=0) * 5)) if len(ids) else np.nan,
                     "radar_R_mean_in_rain_median_mmh": round(float(np.median(wet_r)), 2) if len(wet_r) else np.nan,
                     "radar_R_mean_in_rain_p90_mmh": round(float(np.percentile(wet_r, 90)), 2) if len(wet_r) else np.nan,
                     "pairs_with_baseline_and_rain": int((~b_.missing & b_.link_has_wet_slot).sum()),
                     "pairs_pass_v1_exclusions": len(p_), "links_in_those_pairs": int(p_.cml_id.nunique()), "events_in_those_pairs": int(p_.event_id.nunique()),
                     "tails_v1_main": int(p_.in_v1_main.sum()), "tails_0.6dB_set_S0": int(p_.in_v1_alt.sum()),
                     f"tails_{P2['a_start_min_db']:g}dB_{P2['min_usable_min']}min": int(p_.kept.sum())})
    inv = pd.DataFrame(rows)
    by_band = []
    for label, m in (("v1 main set (1 dB)", pool.in_v1_main), ("0.6 dB set (S0)", pool.in_v1_alt),
                     (f"{P2['a_start_min_db']:g} dB, >= {P2['min_usable_min']} min (part 2)", pool.kept), ("all pairs that pass the v1 exclusions", pool.tail_id >= 0)):
        t = pd.crosstab(pool.lc[m], pool.band[m]).reindex(index=LENGTHS, columns=BANDS, fill_value=0)
        by_band.append(t.assign(tail_set=label).reset_index(names="length"))
    by_band = pd.concat(by_band, ignore_index=True)[["tail_set", "length"] + BANDS]
    E = P1["enough"]
    enough = inv[inv.length.str.startswith("CUT")][["length", "links_in_those_pairs", "events_in_those_pairs", "pairs_pass_v1_exclusions",
                                                     "link_events_with_rain", f"tails_{P2['a_start_min_db']:g}dB_{P2['min_usable_min']}min"]].copy()
    enough["enough"] = ((enough.links_in_those_pairs >= E["min_links"]) & (enough.events_in_those_pairs >= E["min_events"])
                        & (enough.pairs_pass_v1_exclusions >= E["min_pairs"])).to_numpy()
    return inv, by_band, enough


# =================================================================================================== part 2: per tail
def tail_pass(T, long_a, bounds, R5, T5, r_col):
    """Per kept tail: its first FOLLOW + PAD minutes of A from S0, the class of its first 60 min (rules of step 4b)
    and the expected path rain before S0 (dB). `bounds` = the rows of each tail in long_a."""
    Wm = np.full((len(T), FOLLOW + PAD), np.nan)
    cls, exp_db, r_before = [], np.full(len(T), np.nan), np.full(len(T), np.nan)
    nb = P2["path_before_min"] // 5
    for i, r in enumerate(T.itertuples()):
        a = long_a[bounds[i, 0]:bounds[i, 1]]
        seg = a[r.k:min(r.end_k, r.k + FOLLOW + PAD)]
        Wm[i, :len(seg)] = seg
        cls.append(classify(Wm[i, :60], CL)[0]["class"])
        i5 = T5.searchsorted(r.S0_start)
        rr = R5[max(i5 - nb, 0):i5, r_col[r.cml_id]]       # radar R_mean on the path, the 30 min before S0
        if np.isfinite(rr).any():
            a_, b_ = itu_p838(r.frequency_ghz, r.polarization)
            exp_db[i], r_before[i] = r.L_km * float(np.nanmean(a_ * np.power(rr, b_))), float(np.nanmean(rr))
    return Wm, cls, exp_db, r_before


def window_mean(Wm, lo, hi):
    """Mean of A over the minutes [lo, hi) per tail (NaN where the tail has no value there)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(Wm[:, lo:hi], axis=1)


def weighted_median(v, w):
    ok = np.isfinite(v) & (w > 0)
    if not ok.any():
        return np.nan
    v, w = v[ok], w[ok]
    o = np.argsort(v, kind="stable")
    c = np.cumsum(w[o])
    return float(v[o][np.searchsorted(c, 0.5 * c[-1])])


def matched(v, short, strata, mult):
    """The matched values of SHORT and LONG (config: MATCHING). v = the quantity per tail, short = True for a SHORT
    tail (False = LONG), strata = the stratum number of each tail, mult = how often each tail counts (1, or the
    number of times its event was drawn in a bootstrap). Returns (SHORT, LONG, SHORT tails used, LONG tails used)."""
    ok = np.isfinite(v) & (mult > 0)
    n = strata.max() + 1
    ns = np.bincount(strata[ok & short], weights=mult[ok & short], minlength=n)
    nl = np.bincount(strata[ok & ~short], weights=mult[ok & ~short], minlength=n)
    both = (ns > 0) & (nl > 0)
    w = np.where(ok & both[strata], mult, 0.0)
    w_long = w * np.divide(ns, nl, out=np.zeros(n), where=nl > 0)[strata]
    return (weighted_median(v[short], w[short]), weighted_median(v[~short], w_long[~short]),
            int((w[short] > 0).sum()), int((w_long[~short] > 0).sum()))


def compare(T, short, long_, group, band, rng):
    """SHORT against LONG for every matching and quantity, with the event bootstrap (the same draws for all)."""
    g = T[short | long_]
    is_short = short[short | long_].to_numpy()
    ev, ev_names = pd.factorize(g.event_id)
    mults = np.array([np.bincount(rng.integers(0, len(ev_names), len(ev_names)), minlength=len(ev_names)) for _ in range(P2["boot_reps"])], float)
    rows = []
    for kind in MATCHINGS:
        strata = ev * len(A_BINS) + g.a_bin.to_numpy() if kind == MATCHINGS[0] else g.a_bin.to_numpy()
        for q, text in QUANTITIES.items():
            v = g[q].to_numpy(float)
            rs, rl, n_s, n_l = matched(v, is_short, strata, np.ones(len(g)))
            draws = np.array([matched(v, is_short, strata, m[ev])[:2] for m in mults])
            dd, rr = draws[:, 0] - draws[:, 1], draws[:, 0] / np.where(draws[:, 1] > 0, draws[:, 1], np.nan)
            enough = np.isfinite(dd).mean() > 0.5
            lo, hi = np.nanpercentile(dd, P2["interval"]) if enough else (np.nan, np.nan)
            r_lo, r_hi = np.nanpercentile(rr, P2["interval"]) if np.isfinite(rr).mean() > 0.5 else (np.nan, np.nan)
            used = np.isfinite(v)
            rows.append({"group": group, "band": band, "matching": kind, "quantity": text, "short_tails": int((is_short & used).sum()),
                         "long_tails": int((~is_short & used).sum()), "short_used": n_s, "long_used": n_l,
                         "SHORT": rs, "LONG": rl, "difference": rs - rl, "diff_lo": lo, "diff_hi": hi,
                         "ratio": rs / rl if rl > 0 else np.nan, "ratio_lo": r_lo, "ratio_hi": r_hi})
    return rows


def judge(M, group, band):
    """The rule of config part2.rule on one comparison. The matching is chosen first (config main_matching_*)."""
    m = M[(M.group == group) & (M.band == band)].set_index(["matching", "quantity"])
    first = m.loc[(MATCHINGS[0], QUANTITIES["rel_fast_drop"])]
    kind = MATCHINGS[0] if (first.short_used >= P2["main_matching_min_share"] * first.short_tails
                            and first.short_used >= P2["main_matching_min_tails"]) else MATCHINGS[1]
    r, e = m.loc[(kind, QUANTITIES["rel_fast_drop"])], m.loc[(kind, QUANTITIES["exp_path_rel"])]
    width, need = r.diff_hi - r.diff_lo, RULE["power_fraction"] * abs(e.difference)
    power = bool(width < need)
    if not np.isfinite(r.ratio):
        verdict, reason = "undecided", "the LONG drop is not positive, so there is no ratio"
    elif r.ratio <= RULE["path_ratio"]:
        verdict, reason = "fast phase = path rain", f"the SHORT drop is {r.ratio:.2f} x the LONG one (rule: <= {RULE['path_ratio']})"
    elif r.ratio >= RULE["antenna_ratio"] and power:
        verdict, reason = "fast phase = antenna", f"the SHORT drop is {r.ratio:.2f} x the LONG one and the interval is narrow enough"
    elif r.ratio >= RULE["antenna_ratio"]:
        verdict = "undecided"
        reason = (f"the SHORT drop is {r.ratio:.2f} x the LONG one (>= {RULE['antenna_ratio']}), but the power condition fails: the interval of the "
                  f"difference is {width:.3f} wide, and it must be below {need:.3f} (half the expected path-rain difference of {e.difference:+.3f})")
    else:
        verdict, reason = "undecided", f"the SHORT drop is {r.ratio:.2f} x the LONG one, between {RULE['path_ratio']} and {RULE['antenna_ratio']}"
    return {"group": group, "band": band, "matching_used": kind, "short_used": int(r.short_used), "long_used": int(r.long_used),
            "r_SHORT": r.SHORT, "r_LONG": r.LONG, "ratio": r.ratio, "ratio_lo": r.ratio_lo, "ratio_hi": r.ratio_hi,
            "difference": r.difference, "diff_lo": r.diff_lo, "diff_hi": r.diff_hi, "interval_width": width,
            "expected_path_SHORT": e.SHORT, "expected_path_LONG": e.LONG, "expected_path_difference": e.difference,
            "power_needs_width_below": need, "power": power, "verdict": verdict, "reason": reason}


def boot_median(v, events, rng):
    """Median of v with its event-bootstrap interval."""
    ok = np.isfinite(v)
    v, (ev, names) = v[ok], pd.factorize(events[ok])
    if not len(v):
        return np.nan, np.nan, np.nan
    draws = [weighted_median(v, np.bincount(rng.integers(0, len(names), len(names)), minlength=len(names))[ev].astype(float))
             for _ in range(P2["boot_reps"])]
    return (float(np.median(v)), *np.nanpercentile(draws, P2["interval"]))


# =================================================================================================== part 3: whole events
def event_pass(E, TL, col, i_of, blank, R5mean, R5max, r_col, state):
    """Per (event, sub-link) pair of a short link: the peak of A during the rain and the rain's peak; for the LAST
    SHOWER (the wet bins of link_state at the end of the rain, gaps up to last_shower_gap_min bridged) the time of
    the peak of A and of half the peak, and A after the rain; and the test "A rises before the radar shows rain"
    with its placebo. Minute 0 of every pair is t0 = the end of the rain on the link (link_state)."""
    sm_w, hb = P3["smooth_min"], P3["before_min"]
    rows = []
    for r in E.itertuples():
        iw, i0 = i_of(r.wet_start), i_of(r.t0)
        lo, hi = iw - 3 * hb, i0 + P3["after_min"] + PAD   # 3 x before_min: the level and placebo windows of the rise test
        a = np.full(hi - lo, np.nan)
        s0, s1 = max(lo, 0), min(hi, len(TL))
        a[s0 - lo:s1 - lo] = np.where(blank[s0:s1], np.nan, TL[s0:s1, col[r.sublink_id]] - r.baseline_db)
        sm = pd.Series(a).rolling(sm_w, center=True, min_periods=sm_w // 2 + 1).median().to_numpy()
        w, z = iw - lo, i0 - lo                            # positions of the start and the end of the rain in `a`
        c = r_col[r.cml_id]
        rain_mean, rain_max = R5mean[iw // 5:i0 // 5, c], R5max[iw // 5:i0 // 5, c]
        row = {"tail_id": r.tail_id, "rain_minutes": i0 - iw,
               "R_peak_mmh": float(np.nanmax(rain_mean)) if np.isfinite(rain_mean).any() else np.nan,
               "R_max_peak_mmh": float(np.nanmax(rain_max)) if np.isfinite(rain_max).any() else np.nan,
               "radar_rain_mm": float(np.nansum(rain_mean) * 5 / 60) if np.isfinite(rain_mean).any() else np.nan,
               "A_peak_in_rain_db": float(np.nanmax(sm[w:z])) if np.isfinite(sm[w:z]).any() else np.nan}
        if np.isfinite(row["R_peak_mmh"]):
            a_, b_ = itu_p838(r.frequency_ghz, r.polarization)
            row["expected_path_at_R_peak_db"] = r.L_km * a_ * row["R_peak_mmh"] ** b_
        b, gap, first_bin = i0 // 5 - 1, 0, i0 // 5 - 1    # the last shower: walk back over the wet bins
        while b >= iw // 5 and gap <= P3["last_shower_gap_min"] // 5:
            if state[b, c] == 1:
                first_bin, gap = b, 0
            else:
                gap += 1
            b -= 1
        s = first_bin * 5 - lo                             # position of the start of the last shower in `a`
        last_mean = R5mean[first_bin:i0 // 5, c]
        row.update({"last_shower_minutes": i0 - first_bin * 5,
                    "last_shower_R_peak_mmh": float(np.nanmax(last_mean)) if np.isfinite(last_mean).any() else np.nan})
        seg = sm[s:z + P3["after_min"]]                    # the last shower and the 3 h after it
        if np.isfinite(seg).any():
            kp = int(np.nanargmax(seg))
            peak = float(seg[kp])
            row.update({"A_peak_db": peak, "peak_min_after_rain_end": kp + s - z})
            if peak >= P3["min_peak_db"]:
                below = np.nonzero(seg[kp:] < P3["half"] * peak)[0]
                row["half_min_after_rain_end"] = int(below[0]) + kp + s - z if len(below) else np.nan
                for k in [0] + list(P3["after_minutes"]):
                    v = a[z + k:z + k + PAD]
                    row[f"A_over_peak_{k}min"] = float(np.nanmean(v)) / peak if np.isfinite(v).any() else np.nan
        # A rising before the radar shows rain. t_r = first 5 min bin with R_max > 0 at or after wet_start - before_min
        b0 = (iw - hb) // 5
        trace = np.nonzero(np.nan_to_num(R5max[b0:i0 // 5, c]) > 0)[0]
        if len(trace):
            tr = (b0 + int(trace[0])) * 5                  # minute index of t_r
            row["radar_onset_min_before_rain_end"] = i0 - tr
            for name, shift in (("rise", 0), ("placebo", RISE["placebo_shift_min"])):
                t = tr - shift
                quiet = R5max[(t - RISE["level"][0]) // 5:tr // 5, c]      # the radar from the level window up to t_r
                if (t - RISE["level"][0] - lo) < 0 or not np.isfinite(quiet).all() or (quiet > 0).any():
                    continue                               # radar trace (or radar missing) in the hour(s) before: not tested
                lev = a[t - RISE["level"][0] - lo:t - RISE["level"][1] - lo]
                look = sm[t - RISE["search"][0] - lo:t - lo]
                if np.isfinite(lev).sum() < len(lev) // 2 or np.isfinite(look).sum() < len(look) // 2:
                    continue
                k = first_run(np.nan_to_num(look - np.nanmedian(lev), nan=-9.0) >= RISE["rise_db"], RISE["run_min"])
                row[f"{name}_tested"], row[f"{name}_early"] = True, k is not None
                row[f"{name}_lead_min"] = RISE["search"][0] - k if k is not None else np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    for c_ in ("rise_tested", "rise_early", "placebo_tested", "placebo_early"):
        out[c_] = out[c_].fillna(False).astype(bool) if c_ in out else False
    return out


# =================================================================================================== figures
def figure_curves(curves, by_len):
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
    for ax, logy in zip(axes, (False, True)):
        for name, ls in zip(LENGTHS, LINE):
            c = curves[(curves.length == name) & curves.median_A_over_A_start.notna()]
            n = int(by_len.set_index("length").tails[name])
            ax.plot(c.minute, c.median_A_over_A_start, ls=ls, color=INK, lw=1.8, label=f"{name}: {n} tails")
        for lo, hi in (P2["fast_window"], P2["slow_window"], P2["late_window"]):
            ax.axvspan(lo, hi, color="#d9d9d9", lw=0)
            ax.text((lo + hi) / 2, 0.985, f"{lo}-{hi}\nmin", transform=ax.get_xaxis_transform(), fontsize=7, ha="center", va="top")
        ax.set(xlabel="minutes since S0 = t0_trace", xlim=(0, FOLLOW + 5), title="the same, logarithmic" if logy else "median A / A_start")
        if logy:
            ax.set_yscale("log")
            ax.set_ylim(0.03, 1.2)
        style_axes(ax)
    axes[0].set_ylabel("median A / A_start")
    axes[0].legend(frameon=False, fontsize=8.5, loc="center", bbox_to_anchor=(0.6, 0.62), title="line style = length class")
    axes[1].legend(frameon=False, fontsize=7.5, loc="lower left", ncol=2)
    fig.suptitle(f"Step 7a: the median tail by path length, up to 3 h (TRAIN tails from S0, A_start >= {P2['a_start_min_db']} dB, >= {P2['min_usable_min']} min "
                 "long). Grey bands: the three windows of the phase numbers.", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig


def figure_drop(by_len, V):
    fig, ax = plt.subplots(figsize=(9.6, 5.4))
    x = np.arange(len(LENGTHS))
    b = by_len.set_index("length").loc[LENGTHS]
    ax.errorbar(x, b.rel_fast_drop_median, yerr=[b.rel_fast_drop_median - b.rel_fast_drop_lo, b.rel_fast_drop_hi - b.rel_fast_drop_median],
                fmt="o", color=INK, ms=7, capsize=5, lw=1.6, label="measured: median relative fast drop, 90 % event-bootstrap interval (circles)")
    ax.plot(x, b.exp_path_rel_median, "s--", color=INK, mfc="none", ms=7, lw=1.4,
            label="expected from path rain alone: median of a R^b L / A_start (squares, dashed)")
    for i, name in enumerate(LENGTHS):
        ax.text(i + 0.07, b.rel_fast_drop_median[name], f"{b.rel_fast_drop_median[name]:.2f}\n(n = {int(b.tails[name])})", fontsize=8, va="center")
        ax.text(i + 0.07, b.exp_path_rel_median[name] + 0.02, f"{b.exp_path_rel_median[name]:.3f}", fontsize=8, va="bottom")
    v = V.iloc[0]
    ax.text(0.02, 0.83, f"SHORT (< {P2['short_max_km']:g} km) against LONG (>= {P2['long_min_km']:g} km), matched on {v.matching_used}:\n"
            f"relative fast drop {v.r_SHORT:.3f} against {v.r_LONG:.3f}, ratio {v.ratio:.2f}; difference {v.difference:+.3f} ({v.diff_lo:+.3f}..{v.diff_hi:+.3f})\n"
            f"path rain alone would give a difference of {v.expected_path_difference:+.3f}\nverdict by the rule: {v.verdict}",
            transform=ax.transAxes, fontsize=8.5, va="top", bbox=dict(facecolor="white", edgecolor=INK2, lw=0.6))
    ax.set(xticks=x, xticklabels=LENGTHS, xlim=(-0.5, len(LENGTHS) - 0.2), ylim=(0, 1.2), xlabel="path length",
           ylabel="(A_start - A at 15-20 min) / A_start")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    style_axes(ax)
    ax.set_title("Step 7a: the fast first phase of a tail by path length, next to what path rain alone would give", fontsize=10)
    fig.tight_layout()
    return fig


def figure_gallery(picks, TL, col, i_of, blank, R5mean, R5max, r_col):
    n, ncol = len(picks), 4 if len(picks) > 6 else 3
    nrow = int(np.ceil(n / ncol))
    fig = plt.figure(figsize=(4.4 * ncol, 4.3 * nrow))
    gs = fig.add_gridspec(2 * nrow, ncol, height_ratios=[3, 1.2] * nrow, hspace=0.6, wspace=0.25, top=0.93, bottom=0.05, left=0.05, right=0.98)
    for j, r in enumerate(picks.itertuples()):
        ax, bx = fig.add_subplot(gs[2 * (j // ncol), j % ncol]), fig.add_subplot(gs[2 * (j // ncol) + 1, j % ncol])
        iw, i0 = i_of(r.wet_start), i_of(r.t0)
        lo, hi = max(iw - P3["before_min"], 0), min(i0 + P3["after_min"], len(TL))
        a = np.where(blank[lo:hi], np.nan, TL[lo:hi, col[r.sublink_id]] - r.baseline_db)
        x = np.arange(lo, hi) - i0
        ax.plot(x, a, "o", ms=1.8, color="#8a8a8a")
        ax.plot(x, pd.Series(a).rolling(P3["smooth_min"], center=True, min_periods=3).median(), "-", color=INK, lw=1.4)
        ax.axhline(0, color=INK2, lw=0.6)
        xb = np.arange(lo // 5, hi // 5 + 1) * 5 - i0
        rm, rx = np.nan_to_num(R5mean[lo // 5:hi // 5 + 1, r_col[r.cml_id]]), np.nan_to_num(R5max[lo // 5:hi // 5 + 1, r_col[r.cml_id]])
        bx.step(xb[:len(rm)], rm, where="post", color=INK, lw=1.4)
        bx.step(xb[:len(rx)], rx, where="post", color=INK, lw=1.0, ls="--")
        bx.set_ylim(0, max(1.0, float(np.max(np.r_[rx, 0.0])) * 1.15))
        for a_ in (ax, bx):
            a_.axvline(iw - i0, color=INK, ls=":", lw=1.2)
            a_.axvline(0, color=INK, ls="-", lw=1.0)
            a_.set_xlim(x[0], x[-1])
            style_axes(a_)
        ax.text(iw - i0, 0.98, " rain starts", transform=ax.get_xaxis_transform(), fontsize=7, va="top")
        ax.text(0, 0.88, " rain ends", transform=ax.get_xaxis_transform(), fontsize=7, va="top")
        ax.set_title(f"{r.sublink_id}, event {r.event_id[1:]}\nL = {r.L_km:.2f} km, {r.frequency_ghz:.0f} GHz; peak A {r.A_peak_db:.1f} dB, peak R {r.R_peak_mmh:.1f} mm/h", fontsize=8)
        if j % ncol == 0:
            ax.set_ylabel("A (dB)")
            bx.set_ylabel("mm/h")
        if j // ncol == nrow - 1:
            bx.set_xlabel("minutes since the rain ended on the link")
    fig.suptitle(f"Step 7a: {n} whole events on links below {P2['short_max_km']:g} km (seed {P3['gallery']['seed']}), same time axis in each pair of panels. "
                 "Top: A (dots = 1 min, line = 5 min median). Bottom: radar rain on the path (solid = R_mean, dashed = R_max).\n"
                 "Dotted vertical line: the rain on the link starts (link_state); solid vertical line: it ends.", fontsize=10)
    return fig


# =================================================================================================== main
def main():
    sys.stdout.reconfigure(encoding="utf-8")
    clock = time.time()
    rng = np.random.default_rng(P2["seed"])

    # ------------------------------------------------------------------------------------------- data (TRAIN only)
    ti = pd.read_csv(P.TAILS_INDEX, usecols=["tail_id", "split_event", "polarization", "wet_start", "baseline_db"],
                     parse_dates=["wet_start"]).query("split_event == 'train'")                     # VAL rows dropped here
    pool = pd.read_csv(P.TAIL_STARTS, parse_dates=["t0", "S0_start"])[
        ["tail_id", "sublink_id", "cml_id", "event_id", "t0", "L_km", "frequency_ghz", "in_v1_main", "in_v1_alt", "S0_start", "S0_usable_min", "S0_A_start_db"]]
    assert set(pool.tail_id) <= set(ti.tail_id), "tails_starts.csv must hold TRAIN pairs only"
    pool = pool.merge(ti.drop(columns="split_event"), on="tail_id", validate="1:1").sort_values("tail_id").reset_index(drop=True)
    pool["lc"], pool["band"] = length_class(pool.L_km), band_of(pool.frequency_ghz)
    pool["k"] = [minutes(s, t) for s, t in zip(pool.S0_start, pool.t0)]
    pool["end_k"] = pool.k + pool.S0_usable_min
    kept = ((pool.S0_A_start_db >= P2["a_start_min_db"]) & (pool.S0_usable_min >= P2["min_usable_min"])).to_numpy()
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    ev_all = d.events()
    assert set(pool.event_id) <= set(ev_all.event_id[ev_all.split_event == "train"]), "a pair of an event that is not TRAIN"
    ra = d.store("radar_along_links_5min")
    T5 = pd.DatetimeIndex(ra.time.values)
    rmean5, rmax5 = ra.R_mean.load().to_pandas(), ra.R_max.load().to_pandas()
    R5mean, R5max, r_col = rmean5.to_numpy(), rmax5.to_numpy(), {int(c): j for j, c in enumerate(rmean5.columns)}
    print(f"TRAIN only: {len(pool)} pairs that pass the v1 exclusions, {pool.event_id.nunique()} events, {pool.cml_id.nunique()} links; "
          f"{int(kept.sum())} tails with A_start >= {P2['a_start_min_db']} dB and >= {P2['min_usable_min']} min; loaded in {time.time() - clock:.0f} s", flush=True)

    # ------------------------------------------------------------------------------------------- part 1
    st = d.store("link_reference_state_5min").state.load()                # per link and 5 min bin: 1 wet, 0 dry, -1 unknown
    assert (pd.DatetimeIndex(st.time.values) == T5).all() and list(st.cml_id.values) == list(rmean5.columns)
    state, st_col = st.to_numpy(), r_col                   # same link order as the radar store
    inv, by_band, enough = inventory(d, pool, kept, state, T5, st_col, R5mean)
    save(inv, P.SHORT_INVENTORY, index=False)
    save(by_band, P.SHORT_BY_BAND, index=False)
    save(enough, P.SHORT_ENOUGH, index=False)
    en = enough.set_index("length").enough
    cut06, cut1 = f"CUT < {P2['very_short_max_km']:g} km", f"CUT < {P2['short_max_km']:g} km"
    enough_text = (f"By the rule there is enough below {P2['very_short_max_km']:g} km." if en[cut06] else
                   f"By the rule there is NOT enough below {P2['very_short_max_km']:g} km; Step 7 must use the links below {P2['short_max_km']:g} km." if en[cut1] else
                   f"By the rule there is not enough below {P2['very_short_max_km']:g} km and not below {P2['short_max_km']:g} km either.")
    print("part 1: " + enough_text, flush=True)

    # ------------------------------------------------------------------------------------------- part 2: the pass, timed first
    T = pool[kept].reset_index(drop=True)
    long = pq.read_table(P.TAILS_1MIN, columns=["tail_id", "minutes_since_t0", "A_db"], filters=[("tail_id", "in", T.tail_id.tolist())]).to_pandas()
    long = long.sort_values(["tail_id", "minutes_since_t0"], kind="stable").reset_index(drop=True)
    long_a = long.A_db.to_numpy("float64")
    bounds = np.column_stack([np.searchsorted(long.tail_id.to_numpy(), T.tail_id.to_numpy(), side=s_) for s_ in ("left", "right")])
    n_t = min(P2["timing_pairs"], len(T))
    t0 = time.time()
    first = tail_pass(T.head(n_t), long_a, bounds[:n_t], R5mean, T5, r_col)
    per = (time.time() - t0) / n_t
    n_short_pairs = int((pool.L_km < P2["short_max_km"]).sum())
    proj = per * len(T) + 3 * per * n_short_pairs + 90     # part 3 costs about 3 x per pair; bootstraps and figures about 90 s
    cuts = []
    print(f"timing: {per * 1e3:.1f} ms per pair on {n_t} pairs -> the whole step is projected at {proj / 60:.1f} min (budget {S['time_budget_min']} min)", flush=True)
    for cut in ("the Part 3 gallery to 6 events", "the < 0.6 km subclass", "the 3 h ratios"):
        if proj <= S["time_budget_min"] * 60:
            break
        cuts.append(cut)
        proj *= 0.8
    rest = tail_pass(T.iloc[n_t:], long_a, bounds[n_t:], R5mean, T5, r_col)
    Wm = np.vstack([first[0], rest[0]])
    T["class"], T["exp_path_db"], T["radar_R_mean_before_mmh"] = first[1] + rest[1], np.r_[first[2], rest[2]], np.r_[first[3], rest[3]]
    a_start = T.S0_A_start_db.to_numpy()
    a_check = float(np.nanmax(np.abs(window_mean(Wm, 0, PAD) - a_start)))
    assert a_check < 1e-3, f"A_start here differs from tails_starts.csv by {a_check} dB"
    fast, slow, late = (window_mean(Wm, *P2[w]) for w in ("fast_window", "slow_window", "late_window"))
    if "the 3 h ratios" in cuts:
        late[:] = np.nan
    ok = fast >= P2["ratio_floor_db"]
    T["A_fast_db"], T["A_slow_db"], T["A_late_db"] = fast, slow, late
    T["fast_drop_db"], T["rel_fast_drop"] = a_start - fast, (a_start - fast) / a_start
    T["slow_ratio"], T["late_ratio"] = np.where(ok, slow / np.where(ok, fast, 1), np.nan), np.where(ok, late / np.where(ok, fast, 1), np.nan)
    T["exp_path_rel"] = T.exp_path_db / a_start
    T["a_bin"] = np.clip(np.digitize(a_start, A_BINS) - 1, 0, len(A_BINS) - 1)
    for k in P2["shape_minutes"]:
        T[f"ratio_{k}"] = (window_mean(Wm, k - PAD, k) if k == FOLLOW else window_mean(Wm, k, k + PAD)) / a_start
    save(T.drop(columns=["k", "end_k", "t0", "wet_start", "baseline_db"]).round(4), P.SHORT_TAILS, index=False)

    # per length class
    len_rows, curve_rows = [], []
    for name in LENGTHS:
        m = (T.lc == name).to_numpy()
        g = T[m]
        med, lo, hi = boot_median(g.rel_fast_drop.to_numpy(), g.event_id.to_numpy(), rng)
        row = {"length": name, "tails": len(g), "events": g.event_id.nunique(), "links": g.cml_id.nunique(), "median_A_start_db": g.S0_A_start_db.median(),
               **{f"A_over_A_start_{k}min": g[f"ratio_{k}"].median() for k in P2["shape_minutes"]},
               "tails_at_120min": int(g.ratio_120.notna().sum()), "tails_at_180min": int(g.ratio_180.notna().sum()),
               "fast_drop_db_median": g.fast_drop_db.median(), "rel_fast_drop_median": med, "rel_fast_drop_lo": lo, "rel_fast_drop_hi": hi,
               "slow_ratio_median": g.slow_ratio.median(), "slow_ratio_tails": int(g.slow_ratio.notna().sum()),
               "slow_ratio_pooled": g.A_slow_db[g.A_slow_db.notna()].sum() / g.A_fast_db[g.A_slow_db.notna()].sum(),
               "late_ratio_median": g.late_ratio.median(), "late_ratio_tails": int(g.late_ratio.notna().sum()),
               "late_ratio_pooled": g.A_late_db[g.A_late_db.notna()].sum() / g.A_fast_db[g.A_late_db.notna()].sum() if g.A_late_db.notna().any() else np.nan,
               **{f"share_{c}": (g["class"] == c).mean() for c in CLASSES},
               "exp_path_db_median": g.exp_path_db.median(), "exp_path_rel_median": g.exp_path_rel.median(),
               "radar_R_mean_before_median_mmh": g.radar_R_mean_before_mmh.median(), "median_L_km": g.L_km.median()}
        len_rows.append(row)
        ratio = Wm[m, :FOLLOW] / a_start[m, None]
        n_min = np.isfinite(ratio).sum(axis=0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            med_curve = np.where(n_min >= 10, np.nanmedian(ratio, axis=0), np.nan)
        curve_rows.append(pd.DataFrame({"length": name, "minute": np.arange(FOLLOW), "median_A_over_A_start": med_curve, "tails": n_min}))
    by_len = pd.DataFrame(len_rows)
    curves = pd.concat(curve_rows, ignore_index=True)
    save(by_len.round(4), P.SHORT_BY_LENGTH, index=False)
    save(curves.round(4), P.SHORT_CURVES, index=False)

    # SHORT against LONG
    short, long_ = T.L_km < P2["short_max_km"], T.L_km >= P2["long_min_km"]
    vshort = T.L_km < P2["very_short_max_km"]
    groups = [(f"SHORT < {P2['short_max_km']:g} km", short)]
    vs_ok = bool(vshort.sum() >= P2["very_short_min_tails"] and T.event_id[vshort].nunique() >= P2["very_short_min_events"])
    if vs_ok and "the < 0.6 km subclass" not in cuts:
        groups.append((f"VERY SHORT < {P2['very_short_max_km']:g} km", vshort))
    band_share = pd.DataFrame({"SHORT": T.band[short].value_counts(normalize=True), "LONG": T.band[long_].value_counts(normalize=True)}).reindex(BANDS).fillna(0.0)
    uneven = bool((band_share.SHORT - band_share.LONG).abs().max() > P2["band_uneven_share"])
    m_rows, v_rows = [], []
    for gname, gm in groups:
        bands = ["all bands"] + [b for b in BANDS if (gm & (T.band == b)).sum() >= P2["band_min_tails"] and (long_ & (T.band == b)).sum() >= P2["band_min_tails"]]
        for b in bands:
            bm = (T.band == b) if b != "all bands" else pd.Series(True, index=T.index)
            m_rows += compare(T, gm & bm, long_ & bm, gname, b, rng)
    M = pd.DataFrame(m_rows)
    for (gname, b), _ in M.groupby(["group", "band"], sort=False):
        v_rows.append(judge(M, gname, b))
    V = pd.DataFrame(v_rows)
    save(M.round(4), P.SHORT_MATCHED, index=False)
    save(V.round(4), P.SHORT_VERDICT, index=False)
    v = V.iloc[0]                                          # the rule is read on SHORT < 1 km, all bands
    print(f"part 2: SHORT {v.r_SHORT:.3f}, LONG {v.r_LONG:.3f}, ratio {v.ratio:.2f}, difference {v.difference:+.3f} ({v.diff_lo:+.3f}..{v.diff_hi:+.3f}), "
          f"expected from path rain {v.expected_path_difference:+.3f} -> {v.verdict}", flush=True)

    # ------------------------------------------------------------------------------------------- part 3
    E = pool[pool.L_km < P2["short_max_km"]].reset_index(drop=True)
    subs = sorted(E.sublink_id.unique())
    s1 = d.store("cml_1min")
    T1 = pd.DatetimeIndex(s1.time.values)
    assert T1[0] == T5[0], "the 1 min and 5 min stores must start at the same time"
    TL = s1.tl_mean.sel(sublink_id=subs).load().to_numpy()                # (minutes of the season, short sub-links)
    col = {s_: j for j, s_ in enumerate(subs)}
    i_of = lambda t: int((t - T1[0]) / MIN)                               # noqa: E731  minute index of a time
    blank = np.zeros(len(T1), bool)                                       # minutes inside an event that is not TRAIN
    for e in ev_all[ev_all.split_event != "train"].itertuples():
        blank[max(i_of(e.start), 0):max(i_of(e.end), 0)] = True
    # A(t) built here from the 1 min store must be the A of the tails (three short tails of part 2)
    chk = T[T.L_km < P2["short_max_km"]].head(3)
    for r in chk.itertuples():
        i = int(np.nonzero(T.tail_id.to_numpy() == r.tail_id)[0][0])
        mine = TL[i_of(r.S0_start):i_of(r.S0_start) + 30, col[r.sublink_id]] - r.baseline_db
        assert np.nanmax(np.abs(mine - Wm[i, :30])) < 1e-3, "A(t) of part 3 differs from the tails"
    ev3 = E[["tail_id", "sublink_id", "cml_id", "event_id", "L_km", "frequency_ghz", "band", "lc", "wet_start", "t0"]].merge(
        event_pass(E, TL, col, i_of, blank, R5mean, R5max, r_col, state), on="tail_id")
    save(ev3.round(4), P.SHORT_EVENTS, index=False)
    with_r = ev3[ev3.R_peak_mmh.notna() & ev3.A_peak_in_rain_db.notna()].copy()
    r_edges = [0.0] + list(P3["r_peak_bins_mmh"]) + [np.inf]
    with_r["R_bin"] = pd.cut(with_r.R_peak_mmh, r_edges, right=False)
    with_r["A_minus_path"] = with_r.A_peak_in_rain_db - with_r.expected_path_at_R_peak_db
    peaks = with_r.groupby("R_bin", observed=True).agg(pairs=("tail_id", "size"), events=("event_id", "nunique"), links=("cml_id", "nunique"),
                                                    R_peak_median_mmh=("R_peak_mmh", "median"), A_peak_median_db=("A_peak_in_rain_db", "median"),
                                                    A_peak_p90_db=("A_peak_in_rain_db", lambda x: x.quantile(0.9)), A_peak_max_db=("A_peak_in_rain_db", "max"),
                                                    expected_path_median_db=("expected_path_at_R_peak_db", "median"),
                                                    A_peak_minus_path_median_db=("A_minus_path", "median"),
                                                    A_peak_minus_path_p90_db=("A_minus_path", lambda x: x.quantile(0.9))).reset_index()
    peaks["R_bin"] = peaks.R_bin.astype(str)
    save(peaks.round(3), P.SHORT_PEAKS, index=False)
    rho_peak = float(spearman(with_r.R_peak_mmh, with_r.A_peak_in_rain_db))
    rho_excess = float(spearman(with_r.R_peak_mmh, with_r.A_minus_path))
    strong = ev3[ev3.A_peak_db >= P3["min_peak_db"]].copy()

    def drying_rows(key, names):
        rows = []
        for name in names:
            g = strong[strong[key] == name]
            if not len(g):
                continue
            half = g.half_min_after_rain_end.fillna(np.inf)     # not halved within 3 h: later than everything
            med = lambda x: float(x.median()) if np.isfinite(x.median()) else np.nan   # noqa: E731
            rows.append({"split": key, "group": name, "pairs": len(g), "events": g.event_id.nunique(), "links": g.cml_id.nunique(),
                         "last_shower_minutes_median": g.last_shower_minutes.median(), "last_shower_R_peak_median_mmh": g.last_shower_R_peak_mmh.median(),
                         "A_peak_median_db": g.A_peak_db.median(), "peak_min_after_rain_end_median": g.peak_min_after_rain_end.median(),
                         "share_peak_after_rain_end": (g.peak_min_after_rain_end > 0).mean(),
                         "share_halved_within_3h": g.half_min_after_rain_end.notna().mean(),
                         "half_min_after_rain_end_median": med(half), "half_min_after_peak_median": med(half - g.peak_min_after_rain_end),
                         **{f"A_over_peak_{k}min_median": g[f"A_over_peak_{k}min"].median() for k in [0] + list(P3["after_minutes"])}})
        return rows

    strong["all"] = "all short links"
    p3_bands = pd.DataFrame(drying_rows("all", ["all short links"]) + drying_rows("band", BANDS) + drying_rows("lc", LENGTHS[:2]))
    save(p3_bands.round(3), P.SHORT_P3_BANDS, index=False)
    rise_rows = []
    for name, g in [("all short links", ev3)] + [(n_, ev3[ev3.lc == n_]) for n_ in LENGTHS[:2]] + [(b, ev3[ev3.band == b]) for b in BANDS if (ev3.band == b).any()]:
        rise_rows.append({"group": name, "pairs": len(g), "tested": int(g.rise_tested.sum()), "early_rise": int(g.rise_early.sum()),
                          "share_early_rise": g.rise_early[g.rise_tested].mean() if g.rise_tested.any() else np.nan,
                          "median_lead_min": g.rise_lead_min.median() if "rise_lead_min" in g else np.nan,
                          "placebo_tested": int(g.placebo_tested.sum()), "placebo_early": int(g.placebo_early.sum()),
                          "share_placebo_early": g.placebo_early[g.placebo_tested].mean() if g.placebo_tested.any() else np.nan,
                          "events_tested": int(g.event_id[g.rise_tested].nunique()), "events_with_early_rise": int(g.event_id[g.rise_early].nunique())})
    rise = pd.DataFrame(rise_rows)
    save(rise.round(3), P.SHORT_RISE, index=False)
    print(f"part 3: {len(ev3)} pairs on {ev3.cml_id.nunique()} short links, {len(strong)} with a peak >= {P3['min_peak_db']} dB", flush=True)

    # ------------------------------------------------------------------------------------------- figures
    fig = figure_curves(curves, by_len)
    save(fig, P.SHORT_FIG_CURVES, dpi=110)
    plt.close(fig)
    fig = figure_drop(by_len, V)
    save(fig, P.SHORT_FIG_DROP, dpi=110)
    plt.close(fig)
    n_gal = S["gallery_cut"] if "the Part 3 gallery to 6 events" in cuts else P3["gallery"]["events"]
    cand = strong[strong.R_peak_mmh > 0].merge(E[["tail_id", "baseline_db", "polarization"]], on="tail_id")
    cand = cand.sort_values("tail_id").drop_duplicates(["cml_id", "event_id"])       # one direction per link and event
    picks = cand.sample(min(n_gal, len(cand)), random_state=P3["gallery"]["seed"]).sort_values(["L_km", "t0"])
    fig = figure_gallery(picks, TL, col, i_of, blank, R5mean, R5max, r_col)
    save(fig, P.SHORT_FIG_GALLERY, dpi=100)
    plt.close(fig)

    # ------------------------------------------------------------------------------------------- report
    inv_show = inv.set_index("length").T.reset_index(names="")
    whole = [c_ for c_ in by_len if c_ in ("tails", "events", "links") or c_.startswith("tails_at") or c_.endswith("_tails")]
    len_show = by_len.round(3).astype({c_: int for c_ in whole}).astype(object).set_index("length").T.reset_index(names="")
    key = [QUANTITIES[q] for q in ("rel_fast_drop", "fast_drop_db", "exp_path_rel", "exp_path_db")]
    m_show = M.round(3).drop(columns=["ratio_lo", "ratio_hi"])
    main_rows = m_show[(m_show.band == "all bands") & m_show.quantity.isin(key)]
    band_rows = m_show[(m_show.band != "all bands") & (m_show.quantity == QUANTITIES["rel_fast_drop"])]
    slow_rows = m_show[m_show.quantity.isin([QUANTITIES["slow_ratio"], QUANTITIES["late_ratio"]]) & (m_show.matching == v.matching_used)]
    v_show = V.round(3)[["group", "band", "matching_used", "short_used", "long_used", "r_SHORT", "r_LONG", "ratio", "ratio_lo", "ratio_hi", "difference",
                         "diff_lo", "diff_hi", "interval_width", "expected_path_difference", "power_needs_width_below", "power", "verdict"]]
    bl = by_len.set_index("length")
    ra_ = rise.set_index("group").loc["all short links"]
    al = p3_bands.set_index("group").loc["all short links"]
    report = f"""# wet_antenna step 7a: very short links. Enough data, and do they show the fast first phase? (TRAIN only, no law fit)

Written by `research/wet_antenna/short_links.py`. Settings and every rule: `wet_antenna.short_links` in
`research/config.yaml`, written before the run. Why this step: `docs/ROADMAP.md` (Q2, Steps 7 and 7a).

## 1. What was run

- TRAIN links and TRAIN events only; VAL is never read (VAL rows are dropped when a file is read; in part 3 the minutes
  inside an event that is not TRAIN are blanked). No law is fitted. Nothing earlier is changed.
- Part 1: inventory per length class from the release (links, sub-links, link_state, radar on the path) and
  `baselines.csv`, with the pairs of `{rel(P.TAIL_STARTS)}` (the {len(pool)} TRAIN pairs that pass the v1 exclusions).
- Part 2: the {len(T)} tails of those pairs with start S0 (t0_trace), the v1 end rules, A_start (mean of the first 5 min)
  >= {P2['a_start_min_db']} dB and a usable length >= {P2['min_usable_min']} min ({T.event_id.nunique()} events, {T.cml_id.nunique()} links), followed for up to {FOLLOW} min.
- Part 3: the {len(ev3)} pairs on the {ev3.cml_id.nunique()} links below {P2['short_max_km']:g} km, from {P3['before_min']} min before the rain on the link to {P3['after_min']} min after it.
- Timing: {per * 1e3:.1f} ms per pair on the first {n_t} pairs; projection {proj / 60:.1f} min against a budget of {S['time_budget_min']} min.
  Cut: {', '.join(cuts) or 'nothing'}. Run time: {time.time() - clock:.0f} s.
- Checks: A_start computed here equals `tails_starts.csv` (largest difference {a_check:.1e} dB); A(t) built from the 1 min
  store in part 3 equals the A of the tails on three short tails.

## 2. Part 0: what was recorded

Before the analysis, in `docs/ROADMAP.md`: the paper Tiede, Chwala, Siart & Eibert 2025 (Adv. Radio Sci. 23, 13-19)
with its findings and why its method cannot be used here (section "Literature used for the method"); "tails up to
3-6 h, two phases" under Steps 7 and 4d; "frequency band explicitly" under Step 6; the status row 7a.

## 3. Part 1: inventory

Rows = quantities, columns = length classes; the last two columns are the two candidate cuts for Step 7.
"Rain on the link" = a wet bin of link_state inside a TRAIN event.

{md_table(inv_show)}

(`{rel(P.SHORT_INVENTORY)}`.) Tails per length class and frequency band, so that length and frequency are not confused
(`{rel(P.SHORT_BY_BAND)}`):

{md_table(by_band)}

**Enough for Step 7?** Rule (written first): at least {P1['enough']['min_links']} links, {P1['enough']['min_events']} events and {P1['enough']['min_pairs']} pairs that pass the v1
exclusions (`{rel(P.SHORT_ENOUGH)}`):

{md_table(enough)}

**{enough_text}**

## 4. Part 2: the two phases by length

Per tail: fast-phase drop = A_start - mean A over minutes {P2['fast_window'][0]}-{P2['fast_window'][1]}; relative fast drop = that / A_start; slow ratio =
mean A over {P2['slow_window'][0]}-{P2['slow_window'][1]} min / mean A over {P2['fast_window'][0]}-{P2['fast_window'][1]} min; late ratio = mean A over {P2['late_window'][0]}-{P2['late_window'][1]} min / the same (tails that
last). A ratio of one tail needs A at {P2['fast_window'][0]}-{P2['fast_window'][1]} min >= {P2['ratio_floor_db']} dB; the pooled ratio (sum / sum) uses every tail.
Expected path rain = L x mean of a R_mean^b over the {P2['path_before_min']} min before S0 (ITU-R P.838-3).

{md_table(len_show)}

(`{rel(P.SHORT_BY_LENGTH)}`; per tail: `{rel(P.SHORT_TAILS)}`; curves: `{rel(P.SHORT_CURVES)}`.)

![median curves](../figures/{P.SHORT_FIG_CURVES.name})

### SHORT against LONG, matched

SHORT = L < {P2['short_max_km']:g} km, LONG = L >= {P2['long_min_km']:g} km. Inside a stratum every LONG tail is weighted so that both groups have the
same distribution over the strata; the value of a group is its weighted median; {P2['boot_reps']} event bootstraps, 90 % interval
of the difference SHORT - LONG. The rule is read on "{MATCHINGS[0]}" if it keeps at least {P2['main_matching_min_share']:.0%} of the SHORT
tails (and {P2['main_matching_min_tails']}), else on "{MATCHINGS[1]}".

{md_table(main_rows)}

Bands of the two groups (share of the tails): {'; '.join(f'{b}: SHORT {band_share.SHORT[b]:.0%}, LONG {band_share.LONG[b]:.0%}' for b in BANDS)}.
Unevenly spread by the rule (a share differs by more than {P2['band_uneven_share']:.0%}): **{'yes' if uneven else 'no'}**. The relative fast drop inside each band
with at least {P2['band_min_tails']} tails in both groups:

{md_table(band_rows) if len(band_rows) else 'No band has enough tails in both groups.'}

(`{rel(P.SHORT_MATCHED)}`.)

### The rule

- "fast phase = antenna" if the matched relative fast drop of SHORT is >= {RULE['antenna_ratio']} x that of LONG AND the interval of the
  difference is narrower than {RULE['power_fraction']} x the difference path rain alone would give (power);
- "fast phase = path rain" if the SHORT drop is <= {RULE['path_ratio']} x the LONG one;
- otherwise, or if the power condition fails, "undecided".

{md_table(v_show)}

(`{rel(P.SHORT_VERDICT)}`.) The first row is the rule; the other rows are the same rule on the sub-groups, for information.

**Verdict by the rule: {v.verdict.upper()}** ({v.reason}).

![fast drop by length](../figures/{P.SHORT_FIG_DROP.name})

### The slow phase (information, no rule)

{md_table(slow_rows)}

## 5. Part 3: whole events on the links below {P2['short_max_km']:g} km (descriptive)

One row per pair: `{rel(P.SHORT_EVENTS)}`. Minute 0 = the end of the rain on the link (t0 of link_state). A is read on
its centred {P3['smooth_min']} min median.

**Peak of A during the rain against the peak of the radar rain** ({len(with_r)} pairs with radar rain; Spearman of the peak A with the
peak R {rho_peak:+.2f}; of the peak A minus the expected path attenuation {rho_excess:+.2f}). Expected path = a R_peak^b L (ITU-R
P.838-3), shown only as a reference:

{md_table(peaks.round(2))}

(`{rel(P.SHORT_PEAKS)}`.)

**Timing and drying, on the last shower** (pairs with a peak of A >= {P3['min_peak_db']} dB; `{rel(P.SHORT_P3_BANDS)}`). The rain on a link in
one event is often several showers over many hours (median {ev3.rain_minutes.median():.0f} min from the first to the last wet bin), so the times
are measured on the LAST shower = the wet bins of link_state at the end of the rain (gaps up to {P3['last_shower_gap_min']} min bridged).
The peak is searched from the start of that shower to {P3['after_min']} min after its end; "half" = the first minute after the
peak below {P3['half']} x the peak (not halved within 3 h counts as later than everything); negative = before the rain ended:

{md_table(p3_bands.round(2))}

**Does A rise before the radar shows rain?** t_r = the first radar trace (R_max > 0) on the path. Early rise = A at
least {RISE['rise_db']} dB above its level of {RISE['level'][0]}-{RISE['level'][1]} min before t_r, for {RISE['run_min']} min in a row, starting {RISE['search'][0]}-{RISE['search'][1]} min before t_r; only pairs
with no radar trace in the hour before t_r. Placebo = the same test {RISE['placebo_shift_min']} min earlier (its share is the false-alarm rate):

{md_table(rise.round(3))}

(`{rel(P.SHORT_RISE)}`.)

![gallery](../figures/{P.SHORT_FIG_GALLERY.name})

## 6. In plain words

PLAIN_WORDS

## 7. Unclear / possibly wrong

UNCLEAR
"""
    plain, unclear = texts(locals())
    report = report.replace("PLAIN_WORDS", plain).replace("UNCLEAR", unclear)
    save(report, P.SHORT_REPORT)
    print("\n" + "=" * 100 + "\n" + report)
    print(f"done in {time.time() - clock:.0f} s")


def texts(q):
    """Sections 6 and 7 of the report: the reading of the run, filled with its numbers."""
    bl, V, M, T, pool, inv, ev3, rise, peaks, p3, bs = (q[k] for k in ("bl", "V", "M", "T", "pool", "inv", "ev3", "rise", "peaks", "p3_bands", "band_share"))
    v = V.iloc[0]
    vv = V[V.group.str.startswith("VERY") & (V.band == "all bands")]
    vb = V[V.group.str.startswith("SHORT") & (V.band != "all bands")]
    short, long_ = T.L_km < P2["short_max_km"], T.L_km >= P2["long_min_km"]
    l_ratio = T.L_km[short].median() / T.L_km[long_].median()
    mm = M[(M.group == v.group) & (M.matching == v.matching_used)].set_index(["band", "quantity"])
    drop_db, exp_db, slow, late = (mm.loc[("all bands", QUANTITIES[k])] for k in ("fast_drop_db", "exp_path_db", "slow_ratio", "late_ratio"))
    slow_b = {b: mm.loc[(b, QUANTITIES["slow_ratio"])] for b in vb.band}
    fmt = lambda col, f="{:.2f}": " / ".join(f.format(bl[col][n_]) for n_ in LENGTHS)   # noqa: E731  one value per length class
    kept_share = (T.lc.value_counts() / pool.lc.value_counts()).reindex(LENGTHS)
    iv = inv.set_index("length")
    c06, c1 = f"CUT < {P2['very_short_max_km']:g} km", f"CUT < {P2['short_max_km']:g} km"
    tcol = f"tails_{P2['a_start_min_db']:g}dB_{P2['min_usable_min']}min"
    al, pb = p3.set_index("group").loc["all short links"], p3.set_index("group")
    ra = rise.set_index("group").loc["all short links"]
    pk = peaks.set_index("R_bin")
    top = pk.iloc[-1]
    floor_db = bl.A_over_A_start_180min.median() * bl.median_A_start_db.median()
    very = (f"Below {P2['very_short_max_km']:g} km alone the ratio is {vv.iloc[0].ratio:.2f} ({vv.iloc[0].ratio_lo:.2f}..{vv.iloc[0].ratio_hi:.2f}; {int(vv.iloc[0].short_used)} tails)."
            if len(vv) else f"The class below {P2['very_short_max_km']:g} km was not compared on its own.")
    bands = "; ".join(f"{r.band}: {r.ratio:.2f} ({r.ratio_lo:.2f}..{r.ratio_hi:.2f})" for r in vb.itertuples())
    slow_bands = "; ".join(f"{b}: {r.SHORT:.2f} against {r.LONG:.2f} ({r.difference:+.2f}, {r.diff_lo:+.2f}..{r.diff_hi:+.2f})" for b, r in slow_b.items())
    plain = f"""- **Is the fast phase the antenna? By the rule: {v.verdict}.** The first half of the "antenna" rule is met (the SHORT drop is
  {v.ratio:.2f} x the LONG one, rule >= {RULE['antenna_ratio']}); the power condition is not: the interval of the difference is {v.interval_width:.3f} wide and
  had to be below {v.power_needs_width_below:.3f}. So no verdict is claimed.
- **What the numbers show all the same: the fast drop does not grow with the path length.** Median fast drop by
  length class ({', '.join(LENGTHS)}): {fmt('fast_drop_db_median')} dB; relative: {fmt('rel_fast_drop_median')}; the median
  path length goes from {bl.median_L_km.iloc[0]:.2f} to {bl.median_L_km.iloc[-1]:.2f} km. Matched (same event, same A_start bin) SHORT against LONG:
  {v.r_SHORT:.3f} against {v.r_LONG:.3f}, ratio {v.ratio:.2f} ({v.ratio_lo:.2f}..{v.ratio_hi:.2f}), in dB {drop_db.SHORT:.2f} against {drop_db.LONG:.2f}. {very}
  A fast phase made of path rain would scale with L: the SHORT links are {l_ratio:.2f} x as long as the LONG ones (medians), and
  the radar-expected path attenuation is {exp_db.ratio:.2f} x. The whole interval of the ratio lies above the "path rain" limit of
  {RULE['path_ratio']}. This is the test with contrast in L that steps 2, 4c and 4c2 could not make, and it speaks for the antenna
  (or for anything else that does not scale with L). It is my reading; the pre-set rule says undecided.
- **Why the rule has no power.** The power condition measures the interval against the path rain the RADAR saw in the
  30 min before the tail: {exp_db.SHORT:.2f} dB on SHORT and {exp_db.LONG:.2f} dB on LONG tails ({v.expected_path_SHORT:.0%} and {v.expected_path_LONG:.0%} of A_start). That rain is so
  light that it predicts a difference of only {v.expected_path_difference:+.3f} in the relative fast drop; the data cannot resolve half of that.
  The same numbers say something else: the measured fast drop ({drop_db.LONG:.2f} dB) is {drop_db.LONG / exp_db.LONG:.0f} times the radar-expected path rain
  on LONG tails and {drop_db.SHORT / exp_db.SHORT:.0f} times on SHORT ones. So the fast phase is not the rain the radar saw. Rain the radar
  does not see would have to give the same dB on a {bl.median_L_km.iloc[0]:.1f} km and on a {bl.median_L_km.iloc[-1]:.1f} km path.
- **The slow phase is what differs with length, weakly.** Slow ratio (A at 50-60 min / A at 15-20 min), matched:
  SHORT {slow.SHORT:.2f} against LONG {slow.LONG:.2f} ({slow.difference:+.2f}, {slow.diff_lo:+.2f}..{slow.diff_hi:+.2f}; the interval holds 0). Inside the bands: {slow_bands}.
  Median A / A_start at 60 min by length: {fmt('A_over_A_start_60min')}. Long links keep more of their level after the first
  20 min. That could be light rain on or near a long path, the lower frequencies of the longest links, or chance.
- **Up to 3 h.** Median A / A_start at 120 and 180 min: {fmt('A_over_A_start_120min')} and {fmt('A_over_A_start_180min')}. From 15-20 min to
  50-60 min A falls to about half (about {37.5 / np.log2(1 / slow.LONG):.0f} min half-life on LONG, {37.5 / np.log2(1 / slow.SHORT):.0f} on SHORT, if exponential); from there to 170-180 min it only halves once
  more (late ratio {late.SHORT:.2f} and {late.LONG:.2f}). So the slow phase is not one exponential either: it slows down again, or it ends on a
  floor. The level at 3 h is about {floor_db:.2f} dB, the size of the baseline uncertainty of a sub-link, so the stretch
  after about 1-2 h says more about the baseline than about water.
- **Enough data for Step 7?** {q['enough_text']} Below {P2['very_short_max_km']:g} km: {int(iv.links_in_those_pairs[c06])} links, {int(iv.events_in_those_pairs[c06])} events, {int(iv.pairs_pass_v1_exclusions[c06])} pairs
  ({int(iv.link_events_with_rain[c06])} link-events with rain), but only {int(iv[tcol][c06])} tails above {P2['a_start_min_db']} dB and two bands ({iv['bands (sub-links)'][c06]} sub-links).
  Below {P2['short_max_km']:g} km: {int(iv.links_in_those_pairs[c1])} links, {int(iv.pairs_pass_v1_exclusions[c1])} pairs, {int(iv[tcol][c1])} tails. The two short classes behave alike in part 2.
  Suggestion: below {P2['very_short_max_km']:g} km as the main set and below {P2['short_max_km']:g} km as the larger sensitivity set, L as an input. For the learners.
- **What Step 7 needs, from the first look at whole events.**
  - *Events*: the rain on a link in one event is several showers over {ev3.rain_minutes.median() / 60:.0f} h (median). Between showers A often goes back
    to about 0 (gallery), also while the radar still shows light rain. The model must run through such events.
  - *Two time scales*: on the last shower A peaks {-al.peak_min_after_rain_end_median:.0f} min BEFORE link_state ends the rain and is at half the peak
    {al.half_min_after_peak_median:.0f} min after the peak; when the rain "ends", A is at {al.A_over_peak_0min_median:.0%} of the peak, 30 min later at {al.A_over_peak_30min_median:.0%}. Most of the
    fall happens while the rain fades, so it cannot be read off without the rain as input. The slow phase of part 2
    (tens of minutes to hours, a few tenths of a dB) comes after that.
  - *Saturation*: no level A_max is visible. The median peak A grows with the peak radar rain in every bin, up to
    {top.A_peak_median_db:.1f} dB at {top.R_peak_median_mmh:.0f} mm/h; minus the expected path attenuation still {top.A_peak_minus_path_median_db:.1f} dB. But on these links the expected path
    part in heavy rain is {top.expected_path_median_db:.1f} dB (not small), and the radar (2 km, 5 min) underestimates a local peak. Step 7 must
    either model the path term also on short links or stay with light and moderate rain, as the roadmap says.
  - *Wetting before rain*: rare. A rises before the first radar trace in {ra.share_early_rise:.1%} of the tested pairs ({int(ra.early_rise)} of {int(ra.tested)}),
    against {ra.share_placebo_early:.1%} for the placebo an hour earlier. It is not a general feature; single cases exist (gallery:
    10039_1 in E20150617_0705 climbs to 2-4 dB hours before the radar shows rain on the path).
  - *Frequency*: only two bands exist below {P2['short_max_km']:g} km (28 and 38 GHz; none at 23 or 32 GHz). After the last shower they look
    alike (half {pb.half_min_after_peak_median['28 GHz']:.0f} and {pb.half_min_after_peak_median['38 GHz']:.0f} min after the peak; {pb.A_over_peak_0min_median['28 GHz']:.0%} and {pb.A_over_peak_0min_median['38 GHz']:.0%} of the peak at the end of the rain). In the
    slow phase of part 2 the SHORT tails at 38 GHz fall further than at 28 GHz ({slow_b['38 GHz'].SHORT:.2f} against {slow_b['28 GHz'].SHORT:.2f}), so the band
    should be an input of Step 7."""
    unclear = f"""1. **The power condition of the rule.** I read "half the expected path-rain difference" as half of E = the matched
   (radar-expected path rain / A_start) of SHORT minus that of LONG. E is small ({v.expected_path_difference:+.3f}) because the rain before a tail is
   light, so the condition fails whatever the data say. A condition on the ratio (is the "path rain" limit {RULE['path_ratio']} outside
   the interval?) would have passed; it was not the rule, and the rule was not changed.
2. **Two readings of "or if the power condition fails".** I applied the power condition to the "antenna" verdict
   only. If it also guards "path rain", nothing changes here (the ratio is far from {RULE['path_ratio']}).
3. **"Does not scale with L" is not yet "antenna".** A rain cell smaller than every link, or a level error shared by
   all links of an event, would not scale with L either. The matching inside events removes the second only partly.
4. **Selection on A_start >= {P2['a_start_min_db']} dB.** The share of the pairs that reach it is {', '.join(f'{x:.0%}' for x in kept_share)} by length class: about
   the same, which fits a level that does not depend on L. But A_start itself may hold some path rain on long links,
   and the A_start bins are wide.
5. **Frequency and length are mixed.** SHORT tails are {bs.SHORT['38 GHz']:.0%} 38 GHz, LONG tails {bs.LONG['38 GHz']:.0%}; no link below {P2['short_max_km']:g} km is below 25 GHz.
   Inside the bands the ratio is {bands}: the same direction, wider intervals (28 GHz reaches down to its lower end).
6. **The expected path rain is an order of magnitude**, not a measurement: radar R_mean of 5 min and 2 km cells, and
   the ITU-R P.838-3 coefficients were typed in from memory in step 4c (checked at 10, 20, 30 GHz only).
7. **Few independent units.** The SHORT tails come from {T.event_id[short].nunique()} events and {T.cml_id[short].nunique()} links; three events hold most tails of v1. The
   bootstrap resamples events, not links.
8. **The ratios of one tail** leave out tails that are below {P2['ratio_floor_db']} dB at 15-20 min (they dried completely), which biases the
   median ratios upward; the pooled ratios in the table use every tail.
9. **Part 3 was changed after the first run** (it has no rule): the peak of the WHOLE rain lies hours before its end
   (median rain span {ev3.rain_minutes.median():.0f} min), so the timing numbers are now taken on the last shower (`last_shower_gap_min`).
   The peak A and the peak R of the first table are still those of the whole rain and need not be at the same time.
10. **"The rain ends" is link_state** (radar >= 0.1 mm/h on the path or a gauge within 1 km, 5 min bins). "A is at {al.A_over_peak_0min_median:.0%} of
    the peak when the rain ends" depends on that definition as much as on the antenna.
11. **The rise test**: {RISE['rise_db']} dB is one RSL step; the placebo is the same link an hour earlier, which is not the same
    weather (air gets more humid as rain comes closer); rain seen first by a gauge is not "radar rain" here.
12. **The paper of Part 0** (Tiede, Chwala, Siart & Eibert 2025) was recorded as the learners gave it. Checked on the
    publisher's page (2026-10-13): the DOI, title, authors, journal and pages are right, and the text supports drops /
    runlets, "the water amount matters more than the drop positions", "several hours" and "very different drying times
    between repeats". Not confirmed there: the three frequencies (the page says K and Ka band), and "temperature,
    humidity and wind change the drying" (the page says only that they could not be held constant). I did not read
    the full paper.
13. **Inventory details**: "use_in_waa_main" of a sub-link = true in at least one TRAIN event; the rain counts use all
    TRAIN events, the pair counts only those with a baseline."""
    return plain, unclear


if __name__ == "__main__":
    main()
