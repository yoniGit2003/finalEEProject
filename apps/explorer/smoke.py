"""Smoke check of the explorer without a browser: load data and build every figure for 5 scenarios, print timings.

    python apps/explorer/smoke.py [--release data/release/openmrg_v1]
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import data as X  # noqa: E402
from cmlpipe.viz import explorer as V  # noqa: E402

LAYERS = dict(tl=True, rsl=True, tsl=True, radar=True, gauge=True, gauge_cum=True, met=True, flags=True, state=True,
              events=True, gaps=True, baseline=True, band=True)
OPTS = dict(palette="default", mode="lines", same_db_span=False, relative=False, log_rain=False, extra={})


def scenarios(d):
    ev = d.events().set_index("event_id")
    e = ev.loc["E20150727_2140"]
    pad = X.pd.Timedelta("3h")
    return [
        ("Torslanda event E20150727_2140, link 10316, 5 min", 10316, e.start - pad, e.end + pad, "5min"),
        ("same event, 10 s", 10316, e.start - pad, e.end + pad, "10s"),
        ("flatline sub-link: link 10116 in E20150729_1905, 1 min", 10116,
         ev.loc["E20150729_1905"].start - pad, ev.loc["E20150729_1905"].end + pad, "1min"),
        ("link 10041 (no TL, unusable), 15 min, July", 10041, "2015-07-01", "2015-08-01", "15min"),
        ("link 10227 whole period, 1 h", 10227, "2015-06-01", "2015-09-01", "1h"),
    ]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--release", default=str(Path(__file__).resolve().parents[2] / "data/release/openmrg_v1"))
    a = ap.parse_args()
    t = time.perf_counter()
    d = X.open_release(a.release)
    links = X.links_table(d)
    print(f"open release + links table: {time.perf_counter() - t:.2f} s")
    ok = True
    for name, cml, t0, t1, res in scenarios(d):
        try:
            (frame, info), t_load = X.timed(X.load_view, d, cml, t0, t1, res, gauge_mode="nearest City + SMHI")
            (fig, notes), t_fig = X.timed(V.link_view_figure, frame, info, LAYERS, OPTS)
            _, t_json = X.timed(fig.to_json)  # what the browser receives
            print(f"OK   {name}: {info['n_rows']} rows, {info['n_values']:,} values, {len(fig.data)} traces; "
                  f"load {t_load:.2f} s + figure {t_fig:.2f} s + serialise {t_json:.2f} s = "
                  f"{t_load + t_fig + t_json:.2f} s" + (f"  [{'; '.join(notes + info['notes'])}]" if notes or info["notes"] else ""))
        except Exception as err:  # noqa: BLE001 - report every scenario
            ok = False
            print(f"FAIL {name}: {type(err).__name__}: {err}")
    # the other tabs
    (m, tm) = X.timed(V.map_figure, links, d.gauges(), "status", 10316, [])
    (summ_share, te) = X.timed(X.event_all_links, d, "E20150727_2140")
    figs, tf = X.timed(V.event_response_figures, *summ_share, "E20150727_2140")
    (av, ta) = X.timed(X.availability, d)
    (_, tah) = X.timed(V.availability_figure, av)
    (bl, tb) = X.timed(X.dry_baseline, d, 10316)
    (_, tbf) = X.timed(V.dry_baseline_figure, *bl)
    print(f"OK   map {tm:.2f} s; event all-links {te + tf:.2f} s ({len(summ_share[0])} links); availability "
          f"{ta + tah:.2f} s ({av.shape[1]} links x {av.shape[0]} days); dry baseline {tb + tbf:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
