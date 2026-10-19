"""Frozen baseline v2 -> outputs/products/ (R1c). Method: METHOD.md; settings: research/config.yaml (section baseline).

baseline_uncertainty.csv : per good or degraded train sub-link, the spread of the v2 baseline error on the
                                    C1 fake events (certified dry weather): offset_std, offset_p90_abs (D <= 6 h),
                                    offset_p90_abs_24h; and the same without fake events that start humid.
                                    Meant as WEIGHTS for later stages, not as a pass / fail.
baselines.csv            : one row per (event, sub-link) for TRAIN and VAL events (significant or not) on
                                    train links (good and degraded sub-links): the baseline, its window, the number of
                                    usable slots, all flags, the sub-link's uncertainty and use_in_waa_main.
A(t) is not stored: baseline.a_of_t(cml_id, event_id, res) computes it through CMLData.

Needs: outputs/products/dew_risk_hourly.csv (make_dew_products.py), outputs/tables/c1_fake_events.csv.gz,
c3_daily_cycle_per_sublink.csv and r1b_fake_event_offsets.csv.gz (validation/run_all.py).
Locations: paths.py; columns: METHOD.md.
Run: python -m research.baseline.make_baseline_table
"""
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from cmlpipe.api import CMLData

HERE = Path(__file__).resolve().parent
from research.baseline import baseline as bl                                                          # noqa: E402
from research.baseline import paths as P                                                # noqa: E402  (all file locations)
from research.common.save import save

REPO = P.REPO
CONFIG = yaml.safe_load(open(P.CONFIG_FILE, encoding="utf-8"))
B = CONFIG["baseline"]
SPH = bl.SPH
HUMID_RH = CONFIG["baseline_validation"]["r1b"]["humid_start_rh_pct"]
SHORT_MAX_H, LONG_H = CONFIG["baseline_validation"]["c1"]["short_max_h"], CONFIG["baseline_validation"]["c1"]["long_h"]
assert B["version"] == "v2", "the products are made with baseline v2 (set baseline.version: v2)"


def load(d):
    return bl.table_inputs(d, ("train",))


def masks(state_col, dew):
    return bl.table_masks(state_col, dew, B)


def uncertainty(s, tl, state, dew):
    """v2 offsets (mean error over each C1 fake event) and their per-sub-link spread."""
    ev = pd.read_csv(P.C1_FAKE_EVENTS, parse_dates=["t0"])
    ev = ev[ev.sublink_id.isin(s.sublink_id)]
    rh = bl.dew_hourly().humidity
    ev["humid_start"] = rh.reindex(ev.t0.dt.floor("1h")).to_numpy() >= HUMID_RH
    off = np.full(len(ev), np.nan)
    for sid, g in ev.groupby("sublink_id"):
        c = int(sid.split("_")[0])
        x = tl[sid].to_numpy()
        _, ok = masks(state[c], dew)
        for k, e in zip(g.index, g.itertuples()):
            n = int(e.D_h * SPH)
            b = bl.baseline_v2(x, ok, e.i0, B)[0]
            seg = x[e.i0:e.i0 + n] - b
            off[ev.index.get_loc(k)] = np.nanmean(seg) if np.isfinite(seg).any() else np.nan
    ev["offset"] = off
    rows = []
    for sid, g in ev.groupby("sublink_id"):
        r = {"sublink_id": sid}
        for tag, gg in (("", g), ("_no_humid_start", g[~g.humid_start])):
            sh = gg[gg.D_h <= SHORT_MAX_H].offset.dropna()
            lg = gg[gg.D_h == LONG_H].offset.dropna()
            r.update({f"n_fake_short{tag}": len(sh), f"n_fake_24h{tag}": len(lg),
                      f"offset_median{tag}": sh.median() if len(sh) else np.nan,
                      f"offset_std{tag}": sh.std() if len(sh) > 1 else np.nan,
                      f"offset_p90_abs{tag}": sh.abs().quantile(0.9) if len(sh) else np.nan,
                      f"offset_p90_abs_24h{tag}": lg.abs().quantile(0.9) if len(lg) else np.nan})
        rows.append(r)
    u = pd.DataFrame(rows)
    return s[["sublink_id", "cml_id", "status", "tsl_behaviour", "length_class", "frequency_band"]].merge(
        u, on="sublink_id", how="left"), ev


def table(d, s, tl, state, dew, amp, unc):
    return bl.event_table(d, s, tl, state, dew, amp, unc, B, ("train", "val"))


def main():
    d = CMLData.open(REPO / "data" / "release" / CONFIG["release"])
    s, tl, state, dew, amp = load(d)
    unc, ev = uncertainty(s, tl, state, dew)
    # consistency with R1b (main sub-links, v2 as adopted): the offsets must be identical
    r1b = pd.read_csv(P.R1B_FAKE_EVENT_OFFSETS, parse_dates=["t0"])
    col = "off_v2_hx" if B["v2"]["humidity_exclusion"] else "off_v2_nohx"
    chk = ev.merge(r1b[["sublink_id", "D_h", "i0", col]], on=["sublink_id", "D_h", "i0"])
    print(f"check vs R1b: {len(chk)} fake events, max |difference| {np.nanmax(np.abs(chk.offset - chk[col])):.2e} dB")
    save(unc.round(4), P.BASELINE_UNCERTAINTY, index=False)
    t = table(d, s, tl, state, dew, amp, unc)
    save(t.round(4), P.BASELINES, index=False)
    print(f"baselines (v2): {len(t)} rows; by split:", t.split_event.value_counts().to_dict())


if __name__ == "__main__":
    main()
