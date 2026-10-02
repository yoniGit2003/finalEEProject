"""Inventory figures (stage inventory_figures): reports/inventory/*.png.

availability_by_length / _by_completeness (daily tl availability per link), length_frequency_split (confounding and
split), distributions (completeness, worst month, dry noise, baseline range, TSL change rate and classes, with the
status cutoffs), flag_gallery (one 10 s example per flag, chosen by fixed rules), network_maps (status, split,
completeness). Reads the inventory tables and, for the gallery, single L1 / L2 day files with netCDF4.
"""

from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import netCDF4  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from cmlpipe import events  # noqa: E402
from cmlpipe.inventory import stage as inventory_runner  # noqa: E402
from cmlpipe.io import schema  # noqa: E402
from cmlpipe.qc import stage as qc_runner  # noqa: E402
from cmlpipe.viz.style import (C_BLUE, C_INK, C_MUTED, C_ORANGE, C_SURF, FLAG_C, SPLIT_C, STATUS_C,  # noqa: E402
                               C_YELLOW, label, style)

warnings.filterwarnings("ignore", category=RuntimeWarning)


def read_l1_column(cfg, day, cml_id, sub):
    """10 s rsl, tsl (NaN where missing) and L2 flags of one sub-link for one day (netCDF4)."""
    with netCDF4.Dataset(cfg.l1_cml_day_path(day)) as nc:
        j = int(np.flatnonzero(nc["cml_id"][:] == cml_id)[0])
        k = [str(s) for s in nc["sublink_id"][:]].index(sub)
        rsl = np.ma.filled(nc["rsl"][:, k, j].astype("float64"), np.nan)
        tsl = np.ma.filled(nc["tsl"][:, k, j].astype("float64"), np.nan)
    with netCDF4.Dataset(qc_runner.l2_day_path(cfg, day)) as nc:
        nc.set_auto_mask(False)
        fl = nc["flags"][:, k, j]
    return rsl, tsl, fl


def pick_spike(cfg, usable: pd.DataFrame, daily: pd.DataFrame):
    """The clearest spike, by a fixed rule: among usable sub-links whose dry tl noise is below the network median, the ISOLATED
    spike-flagged 10 s sample (its unflagged 6 + 6 neighbours span <= 1 dB) with the largest jump |rsl - their median| relative to the
    sub-link's dry tl noise (std of its dry 5 min means). Returns (sublink row, day, sample index, jump / noise)."""
    quiet = usable[usable.dry_tl_std < usable.dry_tl_std.median()].set_index(["cml_id", "sublink_id"])
    cand = daily[(daily.n_spike > 0)].set_index(["cml_id", "sublink_id"])
    cand = cand[cand.index.isin(quiet.index)].reset_index()
    bit = np.uint16(1 << schema.FLAGS["spike"])
    best = (-1.0, None, None, None)
    for day, g in cand.groupby("day"):
        d = pd.Timestamp(day)
        with netCDF4.Dataset(cfg.l1_cml_day_path(d)) as nc:
            ids = nc["cml_id"][:]
            subs = [str(x) for x in nc["sublink_id"][:]]
            rsl_all = np.ma.filled(nc["rsl"][:].astype("float64"), np.nan)
        with netCDF4.Dataset(qc_runner.l2_day_path(cfg, d)) as nc:
            nc.set_auto_mask(False)
            fl_all = nc["flags"][:]
        for r in g.itertuples():
            j, k = int(np.flatnonzero(ids == r.cml_id)[0]), subs.index(r.sublink_id)
            rsl, m = rsl_all[:, k, j], (fl_all[:, k, j] & bit) != 0
            noise = float(quiet.loc[(r.cml_id, r.sublink_id), "dry_tl_std"])
            for i in np.flatnonzero(m):
                nb = np.r_[np.arange(max(i - 6, 0), i), np.arange(i + 1, min(i + 7, rsl.size))]
                ctx = rsl[nb][np.isfinite(rsl[nb]) & ~m[nb]]  # unflagged neighbours only
                # isolated: the unflagged neighbours must be calm (span <= 1 dB), else it is part of a fade, not a spike
                if ctx.size < 6 or not np.isfinite(rsl[i]) or noise <= 0 or np.ptp(ctx) > 1.0:
                    continue
                ratio = abs(rsl[i] - np.median(ctx)) / noise
                if ratio > best[0]:
                    best = (ratio, usable[(usable.cml_id == r.cml_id) & (usable.sublink_id == r.sublink_id)].iloc[0], d, int(i))
    ratio, row, d, i = best
    return row, d, i, ratio


def outages(cfg) -> pd.DataFrame:
    """Network gaps lasting at least 30 min."""
    g = pd.read_csv(cfg.tables_dir / "network_gaps.csv", parse_dates=["start", "end"])
    return g[g.duration_min >= 30]


# =====================================================================================
# Figures + data note
# =====================================================================================
def make_inventory_figures(cfg) -> tuple[list[Path], list]:
    """All inventory figures (availability heatmaps, length vs frequency, distributions with the status cutoffs,
    flag gallery, network maps) into reports/inventory/. Returns (figure paths, gallery examples)."""
    t0 = time.perf_counter()
    P = inventory_runner.out_paths(cfg)
    sub, lk, gg = pd.read_csv(P["sublinks"]), pd.read_csv(P["links"]), pd.read_csv(P["gauges"])
    icfg = cfg["inventory"]
    out = cfg.output_root / cfg["paths"]["reports_dir"] / "inventory"
    out.mkdir(parents=True, exist_ok=True)
    figs = []
    meta = json.loads(Path(str(P["sublinks"]) + ".meta.json").read_text(encoding="utf-8"))  # derived thresholds
    thr = meta["status_thresholds_derived"]

    # --- 1. availability heatmaps ---
    daily = pd.read_csv(cfg.tables_dir / "qc_daily_summary.csv")
    daily["f"] = daily.n_valid_tl / daily.n_samples
    per_link = daily.groupby(["cml_id", "day"]).f.max().unstack("day")  # at least one usable sub-link
    days = pd.to_datetime(per_link.columns)
    out_days = sorted({d.normalize() for d in outages(cfg).start})
    for name, order in (("by_length", lk.sort_values("length_m").cml_id), ("by_completeness", lk.sort_values("completeness_tl_max").cml_id)):
        m = per_link.loc[order.to_numpy()].to_numpy()
        fig, ax = plt.subplots(figsize=(13, 9), facecolor=C_SURF, layout="constrained")
        im = ax.imshow(m, aspect="auto", cmap="Blues", vmin=0, vmax=1, interpolation="nearest",
                       extent=[mdates.date2num(days[0]), mdates.date2num(days[-1] + pd.Timedelta("1D")), len(m), 0])
        for d in out_days:
            ax.axvline(mdates.date2num(d + pd.Timedelta("12h")), color=C_ORANGE, lw=0.8, ls="--")
        ax.xaxis_date()
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
        cb = fig.colorbar(im, ax=ax, shrink=0.7)
        cb.set_label("daily share of valid tl samples (best sub-link)", fontsize=9)
        label(ax, "day (UTC, 2015)", f"links sorted {name.replace('_', ' ')} (top = {'shortest' if name == 'by_length' else 'least complete'})",
              "Data availability per link and day (orange dashes: days with a network outage >= 30 min)")
        figs.append(out / f"availability_{name}.png")
        fig.savefig(figs[-1], dpi=100, facecolor=C_SURF)
        plt.close(fig)

    # --- 2. length vs frequency ---
    r_log = float(np.corrcoef(np.log10(lk.length_m), lk.frequency_ghz_sublink_1)[0, 1])
    r_sp = float(lk.length_m.rank().corr(lk.frequency_ghz_sublink_1.rank()))
    fig = plt.figure(figsize=(10, 8), facecolor=C_SURF, layout="constrained")
    gs = fig.add_gridspec(2, 2, width_ratios=[4, 1], height_ratios=[1, 4])
    ax = fig.add_subplot(gs[1, 0])
    ah = fig.add_subplot(gs[0, 0], sharex=ax)
    av = fig.add_subplot(gs[1, 1], sharey=ax)
    for a in (ax, ah, av):
        style(a)
    for s in ("train", "test", "excluded"):
        d = lk[lk.split_link == s]
        ax.scatter(d.length_m / 1000, d.frequency_ghz_sublink_1, s=22, color=SPLIT_C[s], edgecolor=C_SURF, lw=0.6, label=f"{s} ({len(d)})")
    ax.set_xscale("log")
    for xe in cfg["split"]["length_classes_km"]:
        ax.axvline(xe, color=C_MUTED, lw=0.8, ls=":")
    ax.axhline(cfg["split"]["frequency_split_ghz"], color=C_MUTED, lw=0.8, ls=":")
    cells = lk.groupby(["length_class", "frequency_band"]).size()
    ax.text(0.01, 0.01, "links per cell:\n" + "\n".join(f"{a} x {b}: {n}" for (a, b), n in cells.items())
            + f"\nPearson r(log10 L, f) = {r_log:.2f}; Spearman = {r_sp:.2f}", transform=ax.transAxes, fontsize=8, color=C_INK, va="bottom")
    label(ax, "link length (km, log)", "frequency, sublink_1 (GHz)", "")
    ax.legend(fontsize=8, frameon=False, loc="upper right")
    ah.hist(lk.length_m / 1000, bins=np.logspace(np.log10(0.1), np.log10(16), 30), color=C_BLUE, edgecolor=C_SURF)
    av.hist(lk.frequency_ghz_sublink_1, bins=30, orientation="horizontal", color=C_BLUE, edgecolor=C_SURF)
    label(ah, "", "links", "Length vs frequency by split (confounded: short links are high-frequency)")
    label(av, "links", "", "")
    figs.append(out / "length_frequency_split.png")
    fig.savefig(figs[-1], dpi=100, facecolor=C_SURF)
    plt.close(fig)

    # --- 3. distributions with the status cutoffs (also the calibration view of the data note) ---
    sc = icfg["status"]
    fig, axs = plt.subplots(2, 3, figsize=(15, 8), facecolor=C_SURF, layout="constrained")
    for a in axs.ravel():
        style(a)
    a = axs[0, 0]
    a.hist(100 * sub.completeness_tl, bins=np.arange(0, 101, 2), color=C_BLUE, edgecolor=C_SURF)
    a.axvline(sc["unusable_completeness_tl_pct"], color=STATUS_C["unusable"], ls="--", label="unusable cutoff")
    a.set_yscale("log")
    a.legend(fontsize=8, frameon=False)
    label(a, "JJA tl completeness (%)", "sub-links (log)", "Completeness")
    a = axs[0, 1]
    mm = 100 * sub[["completeness_tl_jun", "completeness_tl_jul", "completeness_tl_aug"]].min(axis=1)
    a.hist(mm, bins=np.arange(0, 101, 2), color=C_BLUE, edgecolor=C_SURF)
    a.axvline(sc["degraded_month_completeness_tl_pct"], color=STATUS_C["degraded"], ls="--", label="degraded cutoff")
    a.set_yscale("log")
    a.legend(fontsize=8, frameon=False)
    label(a, "worst month tl completeness (%)", "sub-links (log)", "Worst month")
    a = axs[0, 2]
    a.hist(sub.dry_tl_std.dropna(), bins=60, color=C_BLUE, edgecolor=C_SURF)
    a.axvline(thr["dry_noise_db"], color=STATUS_C["degraded"], ls="--", label=f"network p95 = {thr['dry_noise_db']:.2f} dB")
    a.set_yscale("log")
    a.legend(fontsize=8, frameon=False)
    label(a, "std of the dry 5 min tl means (dB)", "sub-links (log)", "Dry-weather noise")
    a = axs[1, 0]
    a.hist(sub.dry_baseline_range_db.dropna(), bins=60, color=C_BLUE, edgecolor=C_SURF)
    a.axvline(thr["dry_baseline_range_db"], color=STATUS_C["degraded"], ls="--", label=f"network p95 = {thr['dry_baseline_range_db']:.2f} dB")
    a.set_yscale("log")
    a.legend(fontsize=8, frameon=False)
    label(a, "p90 - p10 of dry-period medians (dB)", "sub-links (log)", "Baseline stability")
    a = axs[1, 1]
    sf = sub.tsl_step_fraction[sub.tsl_step_fraction > 0]
    a.hist(np.log10(sf), bins=np.arange(-6, 0.01, 0.25), color=C_BLUE, edgecolor=C_SURF)
    a.axvline(np.log10(icfg["tsl"]["toggling_step_fraction"]), color=STATUS_C["degraded"], ls="--", label="toggling cutoff (1 %)")
    a.legend(fontsize=8, frameon=False)
    n_tog = int((sub.tsl_step_fraction > icfg["tsl"]["toggling_step_fraction"]).sum())
    label(a, "log10(share of valid samples with a TSL change)", "sub-links",
          f"TSL change rate: {len(sf)} sub-links plotted (>= 1 change), {n_tog} above the toggling cutoff;\n"
          f"not shown: {int((sub.tsl_step_fraction == 0).sum())} with no change, {int(sub.tsl_step_fraction.isna().sum())} without TSL")
    a = axs[1, 2]
    vc = sub.tsl_behaviour.value_counts().reindex(["no_tsl", "constant", "rare_steps", "frequent_steps", "toggling"]).fillna(0)
    a.bar(range(len(vc)), vc.to_numpy(), color=C_BLUE, edgecolor=C_SURF, lw=2)
    for i, v in enumerate(vc.to_numpy()):
        a.text(i, v + 5, f"{int(v)}", ha="center", fontsize=8, color=C_INK)
    a.set_xticks(range(len(vc)), vc.index, fontsize=8)
    label(a, "", "sub-links", "TSL behaviour classes")
    fig.suptitle("Inventory distributions and the status cutoffs", fontsize=12, color=C_INK)
    figs.append(out / "distributions.png")
    fig.savefig(figs[-1], dpi=100, facecolor=C_SURF)
    plt.close(fig)

    # --- 4. flag gallery: one example per flag, fixed rule ---
    # Rule: among sub-links that are not unusable, the one with the highest JJA share of the flag; on that sub-link
    # the day with most samples of that flag; the window is -1 h / +1 h around the START of the longest flagged run,
    # so it shows the normal signal turning into the flagged one (for missing data the densest window is empty).
    # spike: its own rule (pick_spike): the clearest spike, zoomed to +-5 min.
    usable = sub[sub.status != "unusable"]
    fig, axs = plt.subplots(4, 2, figsize=(15, 14), facecolor=C_SURF, layout="constrained")
    gal = []
    for ax, (f, bit) in zip(axs.ravel(), schema.FLAGS.items()):
        style(ax)
        best = usable.sort_values(f"flag_share_{f}", ascending=False).iloc[0]
        if best[f"flag_share_{f}"] <= 0:
            ax.text(0.5, 0.5, f"{f}: never set on a usable sub-link (unit-tested only)", ha="center", transform=ax.transAxes)
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        if f == "spike":
            best, day, c, ratio = pick_spike(cfg, usable, daily)
            rsl, tsl, fl = read_l1_column(cfg, day, best.cml_id, best.sublink_id)
            m = (fl & np.uint16(1 << bit)) != 0
            a, b = max(c - 30, 0), min(c + 31, len(m))  # +-5 min: single-sample spikes stay visible
        else:
            dd = daily[(daily.cml_id == best.cml_id) & (daily.sublink_id == best.sublink_id)]
            day = pd.Timestamp(dd.sort_values(f"n_{f}", ascending=False).day.iloc[0])
            rsl, tsl, fl = read_l1_column(cfg, day, best.cml_id, best.sublink_id)
            m = (fl & np.uint16(1 << bit)) != 0
            rr = events.runs(m)
            c = max(rr, key=lambda ab: ab[1] - ab[0])[0]
            a, b = max(c - 360, 0), min(c + 360, len(m))
        t = pd.date_range(day, periods=8640, freq="10s")[a:b]
        ax.plot(t, tsl[a:b] - rsl[a:b], color=C_BLUE, lw=0.8, label="tl")
        ax.plot(t, -rsl[a:b], color=C_MUTED, lw=0.6, label="-rsl")
        ax.plot(t, tsl[a:b] * 0 + np.nanmedian(tsl[a:b] - rsl[a:b]) - 10 + (tsl[a:b] - np.nanmedian(tsl[a:b])), color=C_YELLOW,
                lw=0.8, label="tsl (offset)")
        mm_ = m[a:b]
        tlw = tsl[a:b] - rsl[a:b]
        y = np.nanmin(tlw) - 2 if np.isfinite(tlw).any() else (np.nanmin(-rsl[a:b]) - 2 if np.isfinite(rsl[a:b]).any() else 0)
        ax.scatter(t[mm_], np.full(mm_.sum(), y), marker="|", s=50, color=FLAG_C[f], label=f"{f} flagged ({int(mm_.sum())})")
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        extra = (f"; jump {ratio:.0f} x dry noise, +-5 min" if f == "spike" else "")
        label(ax, f"UTC {day:%Y-%m-%d}", "dB", f"{f}: link {best.cml_id} {best.sublink_id} (JJA share {100 * best[f'flag_share_{f}']:.2f} %{extra})")
        if f == "spike":
            ax.plot(t, tsl[a:b] - rsl[a:b], color=C_BLUE, lw=0, marker=".", ms=4)
        ax.legend(fontsize=7, frameon=False, loc="upper right")
        gal.append((f, int(best.cml_id), best.sublink_id, f"{day:%Y-%m-%d}"))
    fig.suptitle("Flag gallery: 10 s data around the flagged samples (one example per flag, fixed rule)", fontsize=12, color=C_INK)
    figs.append(out / "flag_gallery.png")
    fig.savefig(figs[-1], dpi=90, facecolor=C_SURF)
    plt.close(fig)

    # --- 5. maps ---
    lg = pd.read_csv(cfg.tables_dir / "link_geometry.csv").merge(lk[["cml_id", "status_best", "split_link", "completeness_tl_min"]], on="cml_id")
    fig, axs = plt.subplots(1, 3, figsize=(18, 7), facecolor=C_SURF, layout="constrained")
    cmap = plt.get_cmap("Blues")
    for a, what in zip(axs, ("status", "split", "completeness")):
        style(a)
        for r in lg.sort_values("length_m", ascending=False).itertuples():
            col = (STATUS_C[r.status_best] if what == "status" else SPLIT_C[r.split_link] if what == "split"
                   else cmap(0.25 + 0.75 * max(min((r.completeness_tl_min - 0.5) / 0.5, 1), 0)))
            a.plot([r.site_0_x_3006_m / 1000, r.site_1_x_3006_m / 1000], [r.site_0_y_3006_m / 1000, r.site_1_y_3006_m / 1000], color=col, lw=1.4)
        city = gg.network == "city"
        a.scatter(gg.x[city] / 1000, gg.y[city] / 1000, marker="^", s=50, color=C_ORANGE, edgecolor=C_SURF, zorder=5, label="City gauge")
        a.scatter(gg.x[~city] / 1000, gg.y[~city] / 1000, marker="s", s=50, color=C_INK, edgecolor=C_SURF, zorder=5, label="SMHI gauge")
        a.set_aspect("equal")
        handles = ([plt.Line2D([], [], color=STATUS_C[k], lw=3, label=f"{k} ({(lk.status_best == k).sum()})") for k in STATUS_C] if what == "status"
                   else [plt.Line2D([], [], color=SPLIT_C[k], lw=3, label=f"{k} ({(lk.split_link == k).sum()})") for k in ("train", "test", "excluded")]
                   if what == "split" else [plt.Line2D([], [], color=cmap(v), lw=3, label=l) for v, l in ((0.25, "<= 50 %"), (0.6, "75 %"), (1.0, "100 %"))])
        a.legend(handles=handles + a.get_legend_handles_labels()[0], fontsize=8, frameon=False, loc="upper right")
        label(a, "x EPSG:3006 (km)", "y EPSG:3006 (km)" if what == "status" else "",
              {"status": "Best sub-link status", "split": "Link split", "completeness": "Worse sub-link tl completeness"}[what])
    figs.append(out / "network_maps.png")
    fig.savefig(figs[-1], dpi=90, facecolor=C_SURF)
    plt.close(fig)

    (out / "gallery_examples.json").write_text(json.dumps(gal, indent=1), encoding="utf-8")
    return figs, gal


def run_inventory_figures(cfg, force: bool, log, records: list[dict]) -> None:
    """Stage inventory_figures: reports/inventory/*.png. Skipped when the inventory tables and the inputs of the
    flag gallery are unchanged (stamp file reports/inventory/inventory_figures.meta.json)."""
    t0 = time.perf_counter()
    P = inventory_runner.out_paths(cfg)
    out = cfg.output_root / cfg["paths"]["reports_dir"] / "inventory"
    sig = inventory_runner._signature([P["sublinks"], P["links"], P["gauges"], cfg.tables_dir / "qc_daily_summary.csv",
                                       cfg.tables_dir / "network_gaps.csv", cfg.tables_dir / "link_geometry.csv"])
    stamp = out / "inventory_figures.meta.json"
    if not force and stamp.exists():
        m = json.loads(stamp.read_text(encoding="utf-8"))
        if m.get("stage_hash") == cfg.stage_hash("inventory_figures") and m.get("input_signature") == sig:
            log.info("inventory_figures: up to date, skipped (use --force to rebuild)")
            records.append({"step": "inventory_figures", "status": "skipped_same_config"})
            return
    figs, gal = make_inventory_figures(cfg)
    stamp.write_text(json.dumps({"stage_hash": cfg.stage_hash("inventory_figures"), "input_signature": sig,
                                 "figures": [f.name for f in figs], "gallery": gal}, indent=2), encoding="utf-8")
    rec = {"step": "inventory_figures", "status": "written", "n_figures": len(figs), "runtime_s": round(time.perf_counter() - t0, 2)}
    log.info("inventory_figures: %d figures (%.0f s)", len(figs), rec["runtime_s"])
    records.append(rec)
