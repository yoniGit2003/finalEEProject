"""WAA_estm step 0: INVENTORY of the data the WAA stage will use. No analysis, no fitting, no selection.

What it does (rerunnable, read-only on data/ and on the other modules):
  1. measures every input of the WAA stage (release tables and Zarr stores, frozen baseline products): dimensions,
     columns, time span, resolution, coverage in links / sub-links;
  2. loads ONE sub-link for ONE rain event and prints, 2 h before the event start to 12 h after the event end, per
     5 min bin: A (baseline-subtracted), link_state, and every flag;  plot -> outputs/figures/step0_example.png;
  3. writes reports/step0_inventory.md (numbers from step 1 and 2, plus the text descriptions) and prints it.

How the example is chosen (references only, never A, so nothing is chosen by how the tail looks): significant TRAIN
event that is one wet run of <= 6 h, train link, pair usable for WAA (`use_in_waa_main`), the link has a gauge
within 1 km, the largest radar rain
total along the path among the candidates whose link stays NOT wet for 12 h after the last wet bin.

Run from the repository root (about a minute; loads all link states and 5 min flags once):
    C:\\Users\\yoni\\miniconda3\\python.exe research/WAA_estm/step0_inventory.py
"""
from pathlib import Path
import inspect
import sys

import numpy as np
import pandas as pd
import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cmlpipe.api import CMLData

MODULE = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE))
import WAA_estm_paths as P                                          # every location of this module
from research.common.save import save
sys.path.insert(0, str(P.BASELINE_MODULE))
import baseline as bl                                               # frozen baseline v2: a_of_t, dew_mask_5min
import baseline_estm_paths as BP                                    # where the frozen baseline products are

PRE_H, POST_H = 2, 12                                               # example window: event start - 2 h ... end + 12 h
MAX_EVENT_H = 6                                                     # example: one short event (a single wet run)
QC_BITS = {1: "missing_rsl", 2: "missing_tsl", 4: "out_of_range", 8: "rx_floor", 16: "spike", 32: "tsl_step",
           64: "flatline", 128: "network_gap"}                     # pipeline/docs/methods.md, section 2
STATE_REASON_BITS = {1: "radar_wet", 2: "gauge_wet", 4: "radar_missing", 8: "margin"}   # methods.md, section 7
STATE_NAME = {0: "dry", 1: "wet", -1: "unknown"}

sys.stdout.reconfigure(encoding="utf-8")                            # the report has δ, R̄: Windows consoles are cp1252
for p in (P.FIGURES, P.TABLES, P.REPORTS):
    p.mkdir(parents=True, exist_ok=True)

cfg = yaml.safe_load(open(P.CONFIG_FILE, encoding="utf-8"))
RELEASE = P.REPO / "data" / "release" / cfg["release"]
d = CMLData.open(RELEASE)
N = {}                                                              # every number that goes into the report


def bits(value, names):
    """Bitmask -> 'name1|name2' ('' if no bit is set)."""
    v = int(value) if not pd.isna(value) else 0
    return "|".join(n for b, n in names.items() if v & b)


def span(ds):
    """'first .. last (n steps)' of a dataset's time axis."""
    t = pd.DatetimeIndex(ds.time.values)
    return f"{t[0]:%Y-%m-%d %H:%M} .. {t[-1]:%Y-%m-%d %H:%M} ({len(t)} steps)"


def md_table(df):
    """A small DataFrame as a markdown table (no extra dependency)."""
    cols = [str(c) for c in df.columns]
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    out += ["| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(out)


print("1. measuring the inputs ...")
# ---------------------------------------------------------------------------------------------- a. A(t) / baseline
b = pd.read_csv(BP.BASELINES, parse_dates=["event_start", "event_end", "last_wet_end"])
N["b_rows"], N["b_sub"], N["b_links"] = len(b), b.sublink_id.nunique(), b.cml_id.nunique()
N["b_events"] = b.event_id.nunique()
N["b_split"] = b.groupby("split_event").event_id.nunique().to_dict()
N["b_missing"] = int(b.missing.sum())
N["b_waa_main"] = int(b.use_in_waa_main.sum())
sig_tr = b[(b.split_event == "train") & b.significant]
N["b_waa_sig_train"] = int((sig_tr.use_in_waa_main & sig_tr.link_has_wet_slot).sum())
N["b_waa_sig_train_sub"] = sig_tr[sig_tr.use_in_waa_main & sig_tr.link_has_wet_slot].sublink_id.nunique()
N["b_flag_counts"] = {f: int(b[f].sum()) for f in ("missing", "excursion", "dew_risk_in_window",
                                                   "dew_risk_during_event", "daily_cycle_large", "degraded",
                                                   "toggling", "use_in_waa_main")}
N["b_status"] = b.drop_duplicates("sublink_id").status.value_counts().to_dict()
N["b_columns"] = list(b.columns)
src, first = inspect.getsourcelines(bl.a_of_t)                      # show THE line that subtracts the baseline
sub = [(first + i, s.rstrip()) for i, s in enumerate(src) if "tl[base.index] - base" in s]
N["sub_line"] = sub[0] if sub else (None, "NOT FOUND (baseline.py changed?)")
unc = pd.read_csv(BP.BASELINE_UNCERTAINTY)
N["unc_rows"] = len(unc)

# ---------------------------------------------------------------------------------------------- b. link_state
ls = d.store("link_reference_state_5min")
st = ls.state.load().to_pandas()                                    # (time, cml_id) int8
rs = ls.reasons.load().to_pandas()
N["ls_dims"], N["ls_span"] = dict(ls.sizes), span(ls)
vals = st.to_numpy().ravel()
N["ls_share"] = {STATE_NAME[k]: round(float((vals == k).mean()), 3) for k in (0, 1, -1)}
unk = st.to_numpy() == -1
rr = rs.to_numpy()
N["ls_unknown_margin_share"] = round(float(((rr[unk] & 8) > 0).mean()), 3)
N["ls_unknown_radar_missing_share"] = round(float(((rr[unk] & 4) > 0).mean()), 3)

# ---------------------------------------------------------------------------------------------- c. flags
c5 = d.store("cml_5min")
fl = c5.flags_or.load().to_numpy()                                  # (time, sublink_id) uint16, OR over the bin
N["c5_dims"], N["c5_span"] = dict(c5.sizes), span(c5)
N["flag_share"] = {n: round(float(((fl & bit) > 0).mean()), 4) for bit, n in QC_BITS.items()}
tl5 = c5.tl_mean.load().to_numpy()
N["tl_nan_share"] = round(float(np.isnan(tl5).mean()), 4)
N["sub_all_nan"] = int(np.isnan(tl5).all(axis=0).sum())
N["c10_dims"], N["c1_dims"] = dict(d.store("cml_10s").sizes), dict(d.store("cml_1min").sizes)
jumps = pd.read_csv(BP.DRY_JUMPS, parse_dates=["start_utc", "end_utc"])
N["jumps"], N["jumps_sub"] = len(jumps), jumps.sublink_id.nunique()
N["jumps_sign"] = jumps.sign.value_counts().to_dict()
dew = pd.read_csv(BP.DEW_RISK_HOURLY, parse_dates=["time_utc"])
N["dew_rows"], N["dew_share"] = len(dew), round(float(dew.dew_risk.mean()), 3)
N["dew_span"] = f"{dew.time_utc.min():%Y-%m-%d %H:%M} .. {dew.time_utc.max():%Y-%m-%d %H:%M}"
gaps = d.table("network_gaps")
N["gaps"] = len(gaps)
qs = d.status()
N["status_all"] = qs.status.value_counts().to_dict()

# ---------------------------------------------------------------------------------------------- d. splits
se, sl = d.table("split_events"), d.table("split_links")
ev = d.events()
N["split_events"] = se.split_event.value_counts().to_dict()
N["split_links"] = sl.split_link.value_counts().to_dict()
N["sig_by_split"] = ev[ev.significant].split_event.value_counts().to_dict()

# ---------------------------------------------------------------------------------------------- e. rain references
links = d.links()
ra = d.store("radar_along_links_5min")
N["radar_dims"], N["radar_span"] = dict(ra.sizes), span(ra)
N["radar_valid_share"] = round(float(ra.R_mean.notnull().mean()), 4)
N["gauge_1km"] = int(links.has_gauge_within_1km.sum())
N["gauge_1km_train"] = int((links.has_gauge_within_1km & (links.split_link == "train")).sum())
N["single_cell"] = int(links.single_radar_cell.sum())
gauges = d.gauges()
N["gauges"] = gauges[[c for c in ("gauge_id", "network", "gauge_type") if c in gauges.columns]].to_dict("records")
lgm = d.table("link_gauge_matches")
N["lgm_rows"] = len(lgm)

# ---------------------------------------------------------------------------------------------- f. metadata
sb = d.sublinks()
N["sub_rows"] = len(sb)
N["pairs_complete"] = int((sb.groupby("cml_id").size() == 2).sum())
N["dir_map"] = sb.groupby("sublink").direction.unique().apply(list).to_dict()
N["pol"] = sb.polarization.value_counts().to_dict()
N["freq_range"] = (round(sb.frequency_ghz.min(), 2), round(sb.frequency_ghz.max(), 2))
N["len_range_km"] = (round(sb.length_m.min() / 1000, 2), round(sb.length_m.median() / 1000, 2),
                     round(sb.length_m.max() / 1000, 2))
ant = sb[["site_0_antenna_diameter_m", "site_1_antenna_diameter_m"]]
N["ant_missing"] = int(ant.isna().any(axis=1).sum())
N["ant_values"] = sorted(pd.unique(ant.to_numpy().ravel()[~np.isnan(ant.to_numpy().ravel())]).round(2).tolist())
fd = links.frequency_ghz_sublink_2 - links.frequency_ghz_sublink_1
N["freq_diff"] = (round(float(fd.abs().min()), 3), round(float(fd.abs().max()), 3))
N["freq_same"] = int((fd.abs() < 1e-6).sum())                       # links whose two directions have one frequency

# ---------------------------------------------------------------------------------------------- g. weather
met = d.store("met_1h")
N["met_vars"], N["met_span"] = [v for v in met.data_vars if not v.startswith("qc_")], span(met)

# ---------------------------------------------------------------------------------------------- h. tail scripts

print("2. the example ...")
# ---------------------------------------------------------------------------------------------- example choice
el = d.table("event_links")[["event_id", "cml_id", "radar_total_mm", "gauge_id", "gauge_distance_to_path_m",
                             "has_reference_gauge"]]
pool = b[(b.split_event == "train") & b.significant & b.use_in_waa_main & b.link_has_wet_slot]
pool = pool.merge(el, on=["event_id", "cml_id"])
dur_h = (pool.event_end - pool.event_start) / pd.Timedelta(hours=1)
pool = pool[pool.has_reference_gauge & (pool.run_count == 1) & (dur_h <= MAX_EVENT_H)].sort_values(["radar_total_mm", "sublink_id"], ascending=[False, True])
ex = None
for r in pool.itertuples():
    after = st[r.cml_id].loc[r.last_wet_end:r.last_wet_end + pd.Timedelta(hours=POST_H)]
    if (after != 1).all():                                          # no wet bin for 12 h after the rain ended
        ex = r
        break
if ex is None:
    raise SystemExit("no example found with the rule above")
t0 = ex.event_start - pd.Timedelta(hours=PRE_H)
t1 = ex.event_end + pd.Timedelta(hours=POST_H)
pad = f"{max(PRE_H, POST_H)}h"

A = bl.a_of_t(ex.cml_id, ex.event_id, res="5min", pad=pad, d=d)[ex.sublink_id].loc[t0:t1 - pd.Timedelta("1ns")]
cm = d.cml([ex.sublink_id], t0, t1, res="5min", variables=["tl_mean", "n_valid_tl", "flags_or"], qc="flags")
lst = d.link_state([ex.cml_id], t0, t1)
rad = d.radar_along([ex.cml_id], t0, t1)
g_res = "5min" if ex.gauge_id != "GbgA" else "15min"
gau = d.gauge([ex.gauge_id], t0, t1, res=g_res).rate_mm_h.to_pandas().iloc[:, 0]
idx = A.index
jump_on = np.zeros(len(idx), bool)
for j in jumps[jumps.sublink_id == ex.sublink_id].itertuples():
    jump_on |= (idx >= j.start_utc) & (idx < j.end_utc)
net = d.network_state(t0, t1).event_id.to_pandas()

w = pd.DataFrame({
    "A_db": A.round(2),
    "tl_mean_db": cm.tl_mean.to_pandas().iloc[:, 0].round(2),
    "n_valid_tl": cm.n_valid_tl.to_pandas().iloc[:, 0],
    "link_state": lst.state.to_pandas().iloc[:, 0].map(STATE_NAME),
    "state_reasons": [bits(v, STATE_REASON_BITS) for v in lst.reasons.to_pandas().iloc[:, 0]],
    "qc_flags_or": [bits(v, QC_BITS) for v in cm.flags_or.to_pandas().iloc[:, 0]],
    "dew_risk_hour": bl.dew_mask_5min(idx),
    "in_dry_jump": jump_on,
    "network_event": net.reindex(idx).replace("", "-"),
    "radar_R_mean_mmh": rad.R_mean.to_pandas().iloc[:, 0].round(2),
    f"gauge_{ex.gauge_id}_mmh": gau.reindex(idx, method="ffill" if g_res == "15min" else None).round(2),
}, index=idx)
w.index.name = "time_utc"
save(w, P.STEP0_EXAMPLE_CSV)
pair_flags = {k: getattr(ex, k) for k in ("baseline_db", "window_h", "n_usable_slots", "missing", "excursion",
                                          "dew_risk_in_window", "dew_risk_during_event", "daily_cycle_large",
                                          "degraded", "toggling", "use_in_waa_main", "offset_std",
                                          "offset_p90_abs")}
N["ex"] = dict(event_id=ex.event_id, sublink_id=ex.sublink_id, cml_id=int(ex.cml_id),
               event_start=ex.event_start, event_end=ex.event_end, last_wet_end=ex.last_wet_end,
               radar_total_mm=round(ex.radar_total_mm, 1), gauge=ex.gauge_id,
               gauge_dist_m=round(ex.gauge_distance_to_path_m), n_candidates=len(pool),
               length_km=round(links.set_index("cml_id").length_m[ex.cml_id] / 1000, 2),
               freq=float(sb.set_index("sublink_id").frequency_ghz[ex.sublink_id]))

pd.set_option("display.width", 250, "display.max_rows", 1000, "display.max_columns", 30)
print(f"\nExample: {ex.sublink_id}, event {ex.event_id}; pair-level flags (one value per event and sub-link):")
print(pair_flags)
print(w.to_string())

# ---------------------------------------------------------------------------------------------- the figure
# Two panels on one time axis (no second y axis): A with the link_state as background, radar + gauge below.
fig, (a1, a2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1.3]})
stc = lst.state.to_pandas().iloc[:, 0]
for ax in (a1, a2):
    for ti, s in stc.items():
        if s == 1:
            ax.axvspan(ti, ti + pd.Timedelta("5min"), color="#4c8ed9", alpha=0.25, lw=0)
        elif s == -1:
            ax.axvspan(ti, ti + pd.Timedelta("5min"), color="#999999", alpha=0.15, lw=0)
    for t, ls_ in ((ex.event_start, "--"), (ex.event_end, "--"), (ex.last_wet_end, ":")):
        ax.axvline(t, color="#333333", ls=ls_, lw=1)
a1.step(A.index, A.to_numpy(), where="post", color="#1a1a1a", lw=1.2)
a1.axhline(0, color="#777777", lw=0.8)
a1.set_ylabel("A = tl_mean - baseline v2 (dB)")
a1.set_title(f"{ex.sublink_id} ({N['ex']['length_km']} km, {N['ex']['freq']:.1f} GHz), event {ex.event_id}. "
             "Background: link_state wet (blue) / unknown (grey) / dry (white). "
             "Dashed: event start / end; dotted: end of last wet bin", fontsize=9)
a2.step(rad.time.values, rad.R_mean.to_pandas().iloc[:, 0].to_numpy(), where="post", color="#1f5fa8", lw=1.2,
        label="radar R_mean along path")
a2.step(gau.index, gau.to_numpy(), where="post", color="#c0661b", lw=1.2,
        label=f"gauge {ex.gauge_id} ({N['ex']['gauge_dist_m']} m from path)")
a2.set_ylabel("rain (mm/h)")
a2.legend(fontsize=8, loc="upper right", frameon=False)
for ax in (a1, a2):
    ax.grid(alpha=0.25)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
a2.set_xlabel("time (UTC)")
fig.tight_layout()
save(fig, P.STEP0_EXAMPLE_PNG, dpi=120)
plt.close(fig)

print("3. writing the report ...")
REL = lambda p: Path(p).resolve().relative_to(P.REPO).as_posix()     # noqa: E731  repo-relative path for the report
x = N["ex"]
gauge_rows = pd.DataFrame(N["gauges"])
example_rows = w.reset_index()
example_rows["time_utc"] = example_rows.time_utc.dt.strftime("%m-%d %H:%M")

report = f"""# WAA_estm step 0: inventory of the inputs

Written by `research/WAA_estm/step0_inventory.py` (rerunnable; every number below is measured by that script
on release `{cfg['release']}` and the frozen baseline products). Inventory only: nothing is analysed, fitted or
selected here.

Notation: TL = TSL - RSL (`tl_mean`), A(t) = TL(t) - baseline = A_rain(t) + δ(t); in a drying tail A_rain = 0, so
A(t) = δ(t) (WAA + baseline residual).

## 2a. A(t) after baseline subtraction (method v2, frozen)

**There is no stored A(t) file.** A(t) is computed on demand from two stored parts:

| part | path | format | produced by | content |
|---|---|---|---|---|
| TL | `data/release/{cfg['release']}/timeseries/cml_<res>.zarr`, variable `tl_mean` | Zarr, dims (time, sublink_id) | pipeline stage `aggregate` (L3) | mean of the valid 10 s TL samples of the bin, dB; NaN if < 80 % of the bin's samples are valid |
| baseline | `{REL(BP.BASELINES)}` | CSV, {N['b_rows']} rows | `research/baseline_estm/make_baseline_table.py` | one flat baseline per (event, sub-link), `baseline_db` (dB) |

- Baseline table: {N['b_rows']} rows = {N['b_events']} events (train {N['b_split'].get('train')}, val
  {N['b_split'].get('val')}; **no test events**) x {N['b_sub']} sub-links on {N['b_links']} **train** links (status:
  {N['b_status']}). Missing baselines: {N['b_missing']}. `use_in_waa_main` = True: {N['b_waa_main']} rows.
  Significant train events with a wet bin on the link and `use_in_waa_main`: **{N['b_waa_sig_train']} pairs on
  {N['b_waa_sig_train_sub']} sub-links**.
- Resolution: the baseline is one number per event; A(t) can be had at any CML resolution (10 s via
  `cml_10s`, 1 min, 5 min, 15 min, 1 h), but the baseline itself was computed from 5 min `tl_mean`.
- Time span: JJA 2015 (1 Jun - 31 Aug), only around events that have a baseline row.
- **How it is subtracted** (`research/baseline_estm/baseline.py`, function `a_of_t`, line {N['sub_line'][0]}):

```python
{N['sub_line'][1].strip()}
```

  i.e. A(t) = tl_mean(t) - baseline_db, one constant per sub-link column, over [event start - pad, event end + pad).
  The release itself has no baseline ("no baseline is applied in this release", DATA_GUIDE section 4), so the
  only subtraction is this line. Columns of `baselines.csv`: {', '.join(f'`{c}`' for c in N['b_columns'])}.
- `{REL(BP.BASELINE_UNCERTAINTY)}`: {N['unc_rows']} sub-links, the spread of the v2 baseline error in dry weather
  (`offset_std`, `offset_p90_abs`, dB), meant as weights.

## 2b. link_state (wet / dry per bin per link)

- Path: `data/release/{cfg['release']}/timeseries/link_reference_state_5min.zarr`; API `d.link_state(cml_ids, start, end)`.
  Produced by `pipeline/src/cmlpipe/events/link_state.py` (stage `link_state`, config section `link_state`).
- Dims {N['ls_dims']}, **5 min**, **per link (cml_id), not per sub-link**. Span {N['ls_span']}.
- Values: `state` int8: 0 dry, 1 wet, -1 unknown; `reasons` uint8 bitmask: 1 radar_wet, 2 gauge_wet,
  4 radar_missing (in the +-30 min window), 8 margin. Shares over all link-bins: {N['ls_share']}.
- Rules (references only, never CML): **wet** if radar `R_max` along the path >= 0.1 mm/h or a gauge within 1 km has
  >= 0.1 mm in the bin; **dry** if not wet and `R_max` == 0 and every gauge within 1 km == 0 (and valid) in the bin
  **and in every bin within +-30 min**; **unknown** otherwise.
- **The 30 min margin is already included, on the dry side only**: wet bins are NOT widened; the bins within 30 min of
  any non-zero rain become `unknown` (reason bit 8), not wet. Of the unknown bins, {N['ls_unknown_margin_share']:.0%}
  carry the margin bit and {N['ls_unknown_radar_missing_share']:.0%} the radar-missing bit. So the first non-wet bin
  after the rain (`last_wet_end` in `baselines.csv`) is normally an `unknown` margin bin, and the first `dry` bin
  comes ~30 min later.

## 2c. Flags

**Per 10 s sample / per CML bin (pipeline QC, L2; CML signals only).** `flags` (10 s, `cml_10s`) and `flags_or`
(bitwise OR over the bin, every L3 resolution), uint16; made by `pipeline/src/cmlpipe/qc/`. Share of 5 min sub-link
bins with the bit set (all 728 sub-links, JJA):

| bit | flag | meaning | blocks the sample? | share of 5 min bins |
|---|---|---|---|---|
| 1 | missing_rsl | RSL missing | yes | {N['flag_share']['missing_rsl']} |
| 2 | missing_tsl | TSL missing | yes (tsl, tl) | {N['flag_share']['missing_tsl']} |
| 4 | out_of_range | RSL outside [-100, -20] dBm or TSL outside [-15, 30] dBm | yes | {N['flag_share']['out_of_range']} |
| 8 | rx_floor | RSL within 1 dB of a receiver floor for >= 10 min (attenuation only a lower bound) | yes (rsl, tl) | {N['flag_share']['rx_floor']} |
| 16 | spike | isolated RSL outlier (<= 3 samples) vs a 5 min rolling median | yes (rsl, tl) | {N['flag_share']['spike']} |
| 32 | tsl_step | TSL change >= 1 dB | no (information) | {N['flag_share']['tsl_step']} |
| 64 | flatline | RSL exactly constant >= 12 h | no (information) | {N['flag_share']['flatline']} |
| 128 | network_gap | missing while >= 5 % of the alive sub-links are missing (collector outage) | yes | {N['flag_share']['network_gap']} |

Blocked samples are excluded from the bin statistics; a bin statistic is NaN unless >= 80 % of its 10 s samples are
valid (`n_valid_tl` gives the count, 30 per 5 min bin). Share of NaN `tl_mean` 5 min bins: {N['tl_nan_share']};
sub-links with no valid 5 min `tl_mean` at all: {N['sub_all_nan']}. `network_gaps` table: {N['gaps']} outage
intervals.

**Per sub-link (pipeline).** `quality_status` (`d.status()`): good / degraded / unusable from completeness, flag
shares and dry noise (methods.md section 3); all 728 sub-links: {N['status_all']}. `tsl_behaviour`: no_tsl / toggling /
constant / rare_steps / frequent_steps.

**Per (event, link) (pipeline).** `event_links`: `frac_rx_floor_event_*`, `frac_spike_event_*`,
`frac_network_gap_event_*`, `valid_fraction_{{event,before,after}}_*`, `cml_ok_{{event,before,after}}` (>= 80 % valid
bins in the event / 3 h before / 3 h after). `event_links_flags_extra`: flatline / tsl_step shares (information).

**Per (event, sub-link) (frozen baseline v2, `baselines.csv`).** Counts over the {N['b_rows']} rows:
{N['b_flag_counts']}.

| flag | meaning | how | resolution |
|---|---|---|---|
| `missing` | no baseline | < 1 h usable dry time in the 24 h and in the 48 h before the event | per event |
| `excursion` | "flagged baseline window": baseline > 1 dB from the sub-link's summer dry median | value check | per event |
| `dew_risk_in_window` | > 25 % of the baseline window's dry slots in dew-risk hours | hourly dew risk on 5 min slots | per event (from hourly) |
| `dew_risk_during_event` | > 25 % of [event start, last wet end + 1 h) in dew-risk hours | same | per event (from hourly) |
| `daily_cycle_large` | dry daily-cycle amplitude > 0.4 dB (R1 C3; good sub-links only) | per sub-link | per sub-link |
| `degraded`, `toggling` | quality status degraded / TSL class toggling | from the release | per sub-link |
| `use_in_waa_main` | not (missing, excursion, daily_cycle_large, degraded, toggling) | rule | per event |

**Time-resolved, not in the pipeline (baseline module products).**
- `{REL(BP.DRY_JUMPS)}`: {N['jumps']} dry-day jumps ("baseline jumps") on {N['jumps_sub']} train sub-links
  (sign: {N['jumps_sign']}): runs >= 3 h of dry slots with TL >= 1 dB above / below the summer dry median; start / end
  times (5 min resolution), weather, `both_directions`, `neighbours_jumping`.
- `{REL(BP.DEW_RISK_HOURLY)}`: {N['dew_rows']} hours ({N['dew_span']}), `dew_risk` = RH >= 90 % or T - Td < 2 degC
  (true in {N['dew_share']:.0%} of hours); "possible dew" flag, **hourly, one station for the whole network**.
- `dry_excursions.csv` (in the `baseline_exp` exploration): an exploration (hourly, good sub-links, 1 / 2 / 3 dB
  above the dry level, episodes >= 3 h). Not a product; modules must not read it.

## 2d. Learn / test split

- Frozen pipeline split (`split_events`, `split_links`; `pipeline/src/cmlpipe/events/stage.py`, seed 42, never
  changed without `--force-split`): events {N['split_events']} (significant events by split: {N['sig_by_split']}),
  stratified by intensity tercile; links {N['split_links']}, stratified by length class x frequency band, parallel
  links kept together, links with both sub-links unusable excluded. Not a spatial hold-out.
- Intended use: fit on train events x train links, tune on val events; unseen events = test x train links; unseen
  links = train events x test links.
- `baselines.csv` already follows it: only train + val events on train links. **Test events and test links have no
  baseline yet**, so A(t) cannot be computed for them with `a_of_t`.
- No WAA-specific split exists.

## 2e. Rain references per link

| reference | path | resolution | coverage |
|---|---|---|---|
| radar along the path | `radar_along_links_{{5min,15min,1h}}.zarr`: `R_mean` (path-weighted), `R_max`, `R_std`, `wet_fraction`, `valid_weight_fraction` | **5 min** native (2 km grid), mm/h, floored at 0 dBZ | dims {N['radar_dims']}, {N['radar_span']}; non-NaN `R_mean` {N['radar_valid_share']:.1%}; {N['single_cell']} links lie in one radar cell |
| City gauges (10) | `gauges_city_{{1min,5min,15min,1h}}.zarr` (`amount_mm`, `rate_mm_h`) | **1 min** native | {N['gauge_1km']} of 364 links ({N['gauge_1km_train']} train links) have a gauge within 1 km (`links.has_gauge_within_1km`) |
| SMHI GbgA | `gauges_smhi_{{15min,1h}}.zarr` | **15 min** native | one station |
| radar at gauges | `radar_at_gauges_*.zarr` | 5 min | 11 gauges |
| matching | `link_gauge_matches` ({N['lgm_rows']} link-gauge pairs within 5 km); `event_links.gauge_id` = matched gauge per (event, link) | - | - |

Gauges: {', '.join(f"{g.get('gauge_id')} ({g.get('network')}, {g.get('gauge_type')})" for g in N['gauges'])}.

## 2f. Link metadata and A / B pairing

- `tables/sublinks` ({N['sub_rows']} rows) and `tables/links` (364 rows); the same values are coordinates of the CML
  Zarr stores. Columns: `frequency_ghz` ({N['freq_range'][0]}-{N['freq_range'][1]} GHz), `polarization`
  ({N['pol']}), `length_m` (min / median / max {N['len_range_km'][0]} / {N['len_range_km'][1]} /
  {N['len_range_km'][2]} km, EPSG:3006), `site_0_antenna_diameter_m`, `site_1_antenna_diameter_m` (values
  {N['ant_values']} m; rows with a missing diameter: {N['ant_missing']}), `tx_site`, `rx_site`, `hop`, `length_class`,
  `frequency_band`, `tsl_behaviour`, `status`.
- Pairing: `sublink_id = "<cml_id>_<1|2>"`; sub-link 1 = raw direction A (site_0 -> site_1), sub-link 2 = raw
  direction B (site_1 -> site_0); measured mapping sub-link -> direction: {N['dir_map']}. Links with both sub-links:
  {N['pairs_complete']} of 364. |f2 - f1| per link: {N['freq_diff'][0]}-{N['freq_diff'][1]} GHz
  ({N['freq_same']} links list the same frequency for both directions). Both directions share the path, the two sites
  (`site_0_antenna_diameter_m`, `site_1_antenna_diameter_m` are per site, not per direction) and the link_state. Whether
  the two directions use the same physical antenna at a site is not stated in the metadata.

## 2g. Weather

- `met_1h.zarr` (SMHI GbgA, central Gothenburg): {', '.join(N['met_vars'])} (+ SMHI qc codes); hourly, {N['met_span']}.
  Units: degC, %, hPa, m/s, degree. API `d.met(start, end)`.
- Loaded and derived in `{REL(BP.DEW_RISK_HOURLY)}` (temperature, humidity, wind speed, dew point, dew-point
  depression, dew_risk, sun elevation, night), and in `dry_jumps.csv` (mean weather per jump).

## 2h. Existing tail-related scripts (the `WAA_exp` exploration folder; not a module input)

- `tail_sample.py`: first look at 12 random drying tails (seed 1) from significant train events with
  `use_in_waa_main` and a wet bin on the link. Tail start = `last_wet_end` (first bin after the link's last wet
  link_state bin in the event; the unknown margin is inside the tail); tail end = next wet bin or 6 h. A = 5 min
  `tl_mean` - `baseline_db`. Writes `tail_sample.csv` (ids, length, frequency, tail start / end, `ended_by`, A at
  start, mean A in the last hour), `tail_sample_panels.png`, `tail_sample_overlay.png` (linear and log y).
- `gauge_radius.py`: how far a gauge can be from a path and still say "no rain on the path". In-event 5 min bins:
  (A) gauge vs gauge P(B wet | A dry) by distance; (B) P(path radar wet | gauge dry) for every link x City gauge;
  (C) links covered within a radius. Writes `gauge_radius_pairs.csv`, `gauge_radius_link_gauge.csv`,
  `gauge_radius_by_distance.csv`, `gauge_radius.png`. Result used as background: the tail start stays the
  pipeline's link_state (no extra gauges).

## 3. Example: {x['sublink_id']}, event {x['event_id']}

Chosen by references only (rule in the script docstring; one wet run, event <= {MAX_EVENT_H} h): {x['n_candidates']} candidate pairs, the largest radar
total along the path ({x['radar_total_mm']} mm) whose link has no wet bin in the 12 h after its rain.
Link {x['cml_id']}: {x['length_km']} km, {x['freq']:.2f} GHz; gauge {x['gauge']} at {x['gauge_dist_m']} m from the
path. Event {x['event_start']} .. {x['event_end']} (network bounds); last wet bin of the link ends
{x['last_wet_end']}. Window: event start - {PRE_H} h .. event end + {POST_H} h, 5 min.

Pair-level flags (one value for the whole event): {pair_flags}

![example](../figures/step0_example.png)

Columns: `A_db` = tl_mean - baseline_db; `n_valid_tl` (of 30); `state_reasons` and `qc_flags_or` are decoded
bitmasks; `dew_risk_hour` = hourly dew risk of GbgA; `in_dry_jump` = inside a `dry_jumps.csv` interval of this
sub-link; `network_event` = network event id. Full table also in `{REL(P.STEP0_EXAMPLE_CSV)}`.

{md_table(example_rows)}

## 4. Missing, unclear or inconsistent

1. **No stored A(t).** A(t) exists only through `baseline.a_of_t`; the WAA stage has to import
   `research/baseline_estm/baseline.py` (module-to-module code import) or rebuild A from `tl_mean` and `baselines.csv`.
2. **No baseline for test events / test links** (and none for events not in `baselines.csv`), so the held-out
   evaluation cannot use A(t) yet.
3. **Sub-link vs link resolution.** A(t) and the QC flags are per sub-link; link_state, radar, gauges and
   `event_links` are per link. Both directions of a link get the same wet / dry label.
4. **Resolution mismatches.** CML 10 s (aggregates 1 min - 1 h); link_state, radar and the baseline input 5 min; City
   gauges 1 min; SMHI gauge 15 min (its amount labels all three 5 min bins in link_state); weather and dew risk 1 h,
   one station for the whole network. Radar has no 1 min product, so a 1 min WAA model would have a 5 min rain input.
5. **Tail start sits in the margin.** `last_wet_end` is the first non-wet bin, which is normally `unknown` (30 min
   margin, reason 8): a tail that starts there still contains bins with light rain (0 < R_max < 0.1 mm/h) or rain
   within 30 min. The first `dry` bin is ~30 min later. Which one is "the tail start" is a WAA-stage decision.
6. **Radar floor.** Rates below 0 dBZ (~0.03 mm/h) are set to 0, so very light drizzle can be labelled dry.
7. **`last_wet_end` is searched only inside the network event** [start, end) (`make_baseline_table.py`), and is
   clipped at the last time step. Wet bins of the link outside the network event bounds are not seen.
8. **Overlapping flag concepts.** `excursion` (baseline level vs summer median, per event), `dry_jumps.csv` (time
   intervals, >= 3 h, >= 1 dB), `dry_excursions.csv` (exploration, hourly, 1 / 2 / 3 dB) and the C4 table measure
   related "dry level shift" ideas with different windows and resolutions; `dew_risk_during_event` covers the event
   plus 1 h after the last wet bin, i.e. it overlaps the start of every tail, and it is true for about half of the
   significant pairs (by design, rain makes the air humid).
9. **`daily_cycle_large` is known for good sub-links only** (False for degraded ones, which are excluded from
   `use_in_waa_main` anyway).
10. **Weather timestamps.** `met_1h` stamps are instants (first stamp 01:00); `dew_risk_hourly.csv` uses each stamp as
    the START of its hour. A possible 1 h ambiguity for dew risk.
11. **Sub-links with no data.** {N['sub_all_nan']} sub-links have no valid 5 min `tl_mean` at all; status over all 728:
    {N['status_all']}. Unusable sub-links have no baseline rows.
12. **Committed table lacks columns:** `r1b_real_events.csv.gz` has no `dew_risk_during_event` (already noted in the
    baseline CATALOG).
13. **Exploration results are not module inputs.** `tail_sample.py` / `gauge_radius.py` hold the learners' tail
    decisions (start = first non-wet bin; end = next wet bin or a maximum length; max 6 h in the sample, where 10 of
    12 tails reached 6 h, i.e. the 6 h cap was often binding). ORGANIZATION.md: modules never read from
    `explorations/`, so the decisions must be re-stated in this module.
14. **Folder layout.** Step 0 first used `scripts/`, `outputs/`, `reports/`; since step 1 the module follows
    ORGANIZATION.md (code at the root, `outputs/{{products,tables,figures,reports}}/`, README.md, CATALOG.md).
15. **Flag gap seen in the example.** In {x['sublink_id']} a dry-weather rise of A (link_state dry, no rain on radar /
    gauge, dew-risk hours) is below the 1 dB `dry_jumps` threshold, so no time-resolved flag marks it; only the
    hourly `dew_risk` is set. A tail that ran into such a period would carry it. (Observation, not analysed.)
"""
save(report, P.STEP0_REPORT)
print("\n" + "=" * 100 + "\n")
print(report)
print(f"\nwrote {REL(P.STEP0_REPORT)}, {REL(P.STEP0_EXAMPLE_PNG)}, {REL(P.STEP0_EXAMPLE_CSV)}")
