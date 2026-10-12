"""Baseline estimation on a SAMPLE of sub-links and the significant TRAIN events (see README.md for how to read
the outputs). The method itself is in baseline.py.

Run from anywhere: python research/baseline_estm/run_sample.py
Writes (method v1 = "two-window"; locations in baseline_estm_paths.py, described in CATALOG.md):
         outputs/tables/sample_two_window_baselines.csv      one row per sub-link x event (main run, X = 1 h)
         outputs/tables/sample_two_window_sensitivity_X.csv  the same for X = 0.5, 1, 2, 3 h
         outputs/reports/sample_two_window_summary.md        the quality numbers
         outputs/figures/sample_two_window_overview.png      the quality numbers as plots
         outputs/figures/sample_two_window_examples/*.png    a few sub-link x event examples
         (method v2 = "pre-event": outputs/tables/sample_pre_event_baselines.csv only)
"""
from pathlib import Path
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from cmlpipe.api import CMLData

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import baseline as bl                                                          # noqa: E402  (this folder)
import baseline_estm_paths as P                                                # noqa: E402  (all file locations)
from research.common.save import save

REPO = P.REPO
CONFIG = yaml.safe_load(open(P.CONFIG_FILE))     # all settings: research/config.yaml
RELEASE = CONFIG["release"]
B = CONFIG["baseline"]
CFG = {"X_H": B["skip_h"],          # usable dry: at least this many hours after the last wet bin
       "WINDOW_H": B["v1"]["window_h"],   # length of the before / after windows
       "MIN_DRY_H": B["min_dry_h"], # a window needs at least this many hours of usable dry samples
       "DIFF_DB": B["v1"]["diff_db"],     # |after - before| above this -> method B (about 1.5 RSL quantization steps)
       "EXC_DB": B["excursion_db"]} # window mean further than this from the usual level -> excursion (no B)
X_VALUES = (0.5, 1.0, 2.0, 3.0)    # sensitivity check of X (this script only)
LINKS_PER_CLASS = B["sample"]["links_per_class"]   # train links per (length class x frequency band)
SEED = B["sample"]["seed"]

for folder in (P.TABLES, P.REPORTS, P.SAMPLE_TWO_WINDOW_EXAMPLES):
    folder.mkdir(parents=True, exist_ok=True)
d = CMLData.open(REPO / "data" / "release" / RELEASE)
period = d.catalog["period"]

# ---- sample: good train links, a few per class; significant train events --------------------------------------
links = d.links(status="good", split="train")
sample = (links.sample(frac=1, random_state=SEED)                     # shuffle with a fixed seed, then take
          .groupby(["length_class", "frequency_band"]).head(LINKS_PER_CLASS))   # the first few of each class
sublinks = [f"{c}_{k}" for c in sample.cml_id for k in (1, 2)]
good_subs = set(d.status().query("status == 'good'").sublink_id)
sublinks = [s for s in sublinks if s in good_subs]                       # a "good" link can have one weaker sub-link
events = d.events(significant=True, split="train")
print(f"sample: {len(sample)} links, {len(sublinks)} sub-links, {len(events)} events")

data = bl.load_sublinks(d, sublinks, period["start"], period["end_exclusive"])
meta = sample.set_index("cml_id")[["length_m", "length_class", "frequency_band"]]
VERSION = sys.argv[1] if len(sys.argv) > 1 else B["version"]     # v1 | v2 (default: research/config.yaml)


def run_v2():
    """v2 on the same sample: one flat baseline per sub-link x event, with its flags -> outputs/tables/sample_pre_event_baselines.csv."""
    T = data[sublinks[0]].index
    met = d.store("met_1h").load()
    m = pd.DataFrame({v: met[v].isel(gauge_id=0).to_series() for v in ("temperature", "humidity")})
    m = m.reindex(T.floor("1h"))
    dew = bl.dew_risk(m.humidity, m.temperature - bl.dew_point(m.temperature, m.humidity), B["dew"])
    f_amp = P.C3_DAILY_CYCLE                                                     # R1 C3 amplitudes
    amp = pd.read_csv(f_amp).set_index("sublink_id").amplitude_db if f_amp.exists() else pd.Series(dtype=float)
    st = d.sublinks().set_index("sublink_id").tsl_behaviour
    status = d.status().set_index("sublink_id").status
    rows = []
    for s in sublinks:
        x, state = data[s].tl.to_numpy(), data[s].state.to_numpy()
        ok_plain = bl.usable_mask(state, B["skip_h"])
        ok = bl.usable_mask(state, B["skip_h"], dew if B["v2"]["humidity_exclusion"] else None)
        summer = float(np.nanmedian(np.where(ok_plain, x, np.nan)))
        for ev in events.itertuples():
            i0 = T.searchsorted(pd.Timestamp(ev.start))
            value, w, a = bl.baseline_v2(x, ok, i0, B)
            rows.append({"sublink_id": s, "event_id": ev.event_id, "event_start": ev.start, "event_end": ev.end,
                         "baseline_db": value, "window_h": w, "missing": bool(np.isnan(value)),
                         "summer_dry_median_db": summer,
                         **bl.flags_v2(value, a, i0, state, ok_plain, dew, summer, amp.get(s, np.nan),
                                       status[s] == "degraded", st[s] == "toggling", B)})
    r = pd.DataFrame(rows)
    save(r, P.SAMPLE_PRE_EVENT_BASELINES, index=False)
    print(f"v2 sample: {len(r)} pairs, missing {r.missing.mean():.1%}, 48 h fallback {(r.window_h == 48).mean():.1%}, "
          f"excursion {r.excursion.mean():.1%}, dew_risk_in_window {r.dew_risk_in_window.mean():.1%}, "
          f"daily_cycle_large {r.daily_cycle_large.mean():.1%}")


if VERSION == "v2":
    run_v2()
    sys.exit(0)


def run(x_hours):
    """All sub-link x event baselines for one value of X."""
    cfg = dict(CFG, X_H=x_hours)
    rows = []
    for s in sublinks:
        df = data[s]
        mask = bl.usable_dry(df.state, x_hours)
        usual = float(df.tl[mask].median())
        for ev in events.itertuples():
            r = bl.estimate_event(df, mask, usual, pd.Timestamp(ev.start), pd.Timestamp(ev.end), cfg)
            rows.append({"sublink_id": s, "event_id": ev.event_id, "X_h": x_hours, **r})
    out = pd.DataFrame(rows)
    out["cml_id"] = out.sublink_id.str.split("_").str[0].astype(int)
    return out.merge(meta, left_on="cml_id", right_index=True, how="left")   # add length / frequency class


runs = {x: run(x) for x in X_VALUES}
main = runs[CFG["X_H"]]
save(main, P.SAMPLE_TWO_WINDOW_BASELINES, index=False)
sens = pd.concat(runs.values(), ignore_index=True)
save(sens, P.SAMPLE_TWO_WINDOW_SENSITIVITY, index=False)


def md(df):
    """A DataFrame as a small markdown table (avoids the optional 'tabulate' dependency of DataFrame.to_markdown)."""
    df = df.reset_index()
    rows = [" | ".join(map(str, df.columns)), " | ".join("---" for _ in df.columns)]
    rows += [" | ".join(f"{v:.2f}" if isinstance(v, float) else str(v) for v in r) for r in df.itertuples(index=False)]
    return "\n".join(f"| {r} |" for r in rows)


# ---- quality numbers ------------------------------------------------------------------------------------------
key = ["sublink_id", "event_id"]
ref = main.set_index(key).before_db
dx = {x: (r.set_index(key).before_db - ref).abs() for x, r in runs.items() if x != CFG["X_H"]}
has = main[main.method != "missing"]
q = lambda s, p: f"{np.nanpercentile(s, p):.2f}" if s.notna().any() else "n/a"   # noqa: E731

lines = [f"# Baseline estimation, sample run ({RELEASE})", "",
         f"Settings: {CFG}", f"Sample: {len(sample)} links / {len(sublinks)} sub-links x {len(events)} significant "
         f"train events = {len(main)} sub-link x event pairs.", "",
         "## 1. Coverage (how often we get a baseline at all)",
         md(main.method.value_counts().rename_axis("method").to_frame("pairs")), "",
         "Flags (no_before = the missing pairs; the others are A pairs worth a look):",
         md(main["flags"].replace("", "none").value_counts().rename_axis("flag").to_frame("pairs")), "",
         "## 2. Precision (spread of the dry TL inside the before window, dB)",
         f"median {q(has.before_std_db, 50)}, 90th percentile {q(has.before_std_db, 90)}. "
         "Compare with the RSL quantization step of ~0.3 dB.", "",
         "## 3. Closure: after - before (how much the dry level moved over the event, dB)",
         f"|closure| median {q(has.closure_db.abs(), 50)}, 90th percentile {q(has.closure_db.abs(), 90)}; "
         f"share above DIFF_DB = {CFG['DIFF_DB']} dB: {(has.closure_db.abs() > CFG['DIFF_DB']).mean():.0%}. "
         "For method A this is the baseline error at the event end.", "",
         "## 4. Out-of-sample: mean A on usable dry bins inside the event (should be ~0, dB)",
         f"pairs with such bins: {has.in_event_dry_A_mean_db.notna().sum()}; median "
         f"{q(has.in_event_dry_A_mean_db, 50)}, 10th / 90th percentile {q(has.in_event_dry_A_mean_db, 10)} / "
         f"{q(has.in_event_dry_A_mean_db, 90)}.", "",
         "## 5. Sensitivity to X (|before(X) - before(1 h)|, dB)"]
lines += [f"- X = {x:g} h: median {q(v, 50)}, 90th percentile {q(v, 90)}, max {q(v, 100)}" for x, v in dx.items()]
lines += ["", "## 6. By link class (median |closure| and median in-event dry A, dB)",
          md(has.groupby(["length_class", "frequency_band"]).agg(
              pairs=("method", "size"), abs_closure=("closure_db", lambda s: s.abs().median()),
              in_event_dry_A=("in_event_dry_A_mean_db", "median")).round(2))]
save("\n".join(lines) + "\n", P.SAMPLE_TWO_WINDOW_SUMMARY)
print("\n".join(lines))

# ---- overview figure ------------------------------------------------------------------------------------------
fig, ax = plt.subplots(2, 3, figsize=(15, 8))
mc = main.method.value_counts().reindex(["A", "B", "missing"], fill_value=0)
ax[0, 0].bar(mc.index, mc.values, color=["C0", "C1", "C7"])
ax[0, 0].set_title("1. coverage: method per sub-link x event")
ax[0, 1].hist(has.before_std_db.dropna(), bins=30, color="C0")
ax[0, 1].axvline(0.3, color="k", ls="--", lw=1, label="RSL step 0.3 dB")
ax[0, 1].set_title("2. precision: std of dry TL in before window")
ax[0, 1].set_xlabel("dB")
ax[0, 1].legend(fontsize=8)
ax[0, 2].hist(has.closure_db.dropna(), bins=40, color="C2")
for v in (-CFG["DIFF_DB"], CFG["DIFF_DB"]):
    ax[0, 2].axvline(v, color="k", ls="--", lw=1)
ax[0, 2].set_title("3. closure: after - before")
ax[0, 2].set_xlabel("dB  (dashed: +-DIFF_DB, beyond -> B)")
ax[1, 0].hist(has.in_event_dry_A_mean_db.dropna(), bins=40, color="C3")
ax[1, 0].axvline(0, color="k", lw=1)
ax[1, 0].set_title("4. out-of-sample: mean A on dry bins inside event")
ax[1, 0].set_xlabel("dB  (should be ~0)")
dur = (has.event_end - has.event_start) / pd.Timedelta("1h")
ax[1, 1].scatter(dur, has.closure_db.abs(), s=8, alpha=0.6)
ax[1, 1].set_title("|closure| vs event duration")
ax[1, 1].set_xlabel("event duration [h]")
ax[1, 1].set_ylabel("|after - before| [dB]")
ax[1, 2].boxplot([v.dropna() for v in dx.values()], showfliers=True)
ax[1, 2].set_xticks(range(1, len(dx) + 1), [f"X={x:g} h" for x in dx])
ax[1, 2].set_title(f"5. sensitivity: |before(X) - before({CFG['X_H']:g} h)|")
ax[1, 2].set_ylabel("dB")
fig.suptitle(f"Baseline estimation, sample: {len(sublinks)} sub-links x {len(events)} train events (X = {CFG['X_H']:g} h)")
fig.tight_layout()
save(fig, P.SAMPLE_TWO_WINDOW_OVERVIEW, dpi=110)
plt.close(fig)


# ---- example figures: a few pairs per method / flag ------------------------------------------------------------
def plot_example(r):
    df = data[r.sublink_id]
    pad = pd.Timedelta(hours=CFG["WINDOW_H"] + 2)
    t0, t1 = r.event_start - pad, r.event_end + pad
    w = df.loc[t0:t1]
    mask = bl.usable_dry(df.state, CFG["X_H"]).loc[w.index]
    base = bl.baseline_series(r._asdict(), w.index)
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(13, 6), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
    for code, col in ((1, "#4c8ed9"), (-1, "#bbbbbb")):
        for a in (a1, a2):
            a.fill_between(w.index, 0, 1, where=(w.state == code).to_numpy(), color=col, alpha=0.25,
                           transform=a.get_xaxis_transform(), lw=0)
    win = pd.Timedelta(hours=CFG["WINDOW_H"])
    for a, b in ((r.event_start - win, r.event_start), (r.event_end, r.event_end + win)):
        a1.axvspan(a, b, color="green", alpha=0.08)
    a1.plot(w.index, w.tl, color="C0", lw=0.8, label="TL (5 min mean)")
    a1.plot(w.index[mask], w.tl[mask], ".", color="green", ms=3, label=f"usable dry (>= {CFG['X_H']:g} h after wet)")
    a1.plot(w.index, base, color="red", lw=2, label=f"baseline, method {r.method}")
    a1.axhline(r.usual_level_db, color="k", ls=":", lw=1, label="usual dry level")
    for t in (r.event_start, r.event_end):
        a1.axvline(t, color="k", ls="--", lw=1)
        a2.axvline(t, color="k", ls="--", lw=1)
    a1.set_ylabel("TL [dB]")
    a1.legend(loc="upper left", fontsize=8)
    a1.set_title(f"{r.sublink_id}  {r.event_id}  method {r.method} {r.flags}\n"
                 f"before {r.before_db:.2f} dB, after {r.after_db:.2f} dB, closure {r.closure_db:+.2f} dB; "
                 f"green bands = before / after windows; blue = wet, grey = unknown")
    a2.plot(w.index, w.tl - base, color="C3", lw=0.8)
    a2.axhline(0, color="k", lw=0.8)
    a2.set_ylabel("A = TL - baseline [dB]")
    a2.set_xlabel("time (UTC)")
    fig.tight_layout()
    save(fig, P.SAMPLE_TWO_WINDOW_EXAMPLES / f"{r.sublink_id}_{r.event_id}.png", dpi=100)
    plt.close(fig)


for f in P.SAMPLE_TWO_WINDOW_EXAMPLES.glob("*.png"):
    f.unlink()                                                          # only the current run's examples
cat = has.assign(cat=has.method + "|" + has["flags"])
picks = cat.sample(frac=1, random_state=SEED).groupby("cat").head(2)
for r in picks.itertuples():
    plot_example(r)
print(f"saved {len(picks)} example figures, the overview figure, the tables and the summary (sample_two_window_* in outputs/)")
