"""Exploration: how common are large dry-time TL excursions (like link 10200 on 7-8 Aug 2015)?

Question: is the dry TL of a sub-link usually confined to a few quantization levels, or do many sub-links have
episodes where TL rises well above its usual dry level with no rain at all?

Method (hourly, sub-links with quality status "good"):
  1. an hour is dry if all twelve 5 min link_state bins in it are dry (radar + gauges, never CML)
  2. the usual dry level of a sub-link = median of its hourly mean TL over all dry hours of JJA 2015
  3. excursion = dry hour with mean TL more than THRESH_DB above that level
  4. an episode = >= MIN_RUN_H consecutive excursion hours
Run from anywhere: python research/explorations/dry_excursions.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from cmlpipe.api import CMLData

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
release = yaml.safe_load(open(REPO / "research" / "config.yaml"))["release"]
d = CMLData.open(REPO / "data" / "release" / release)

THRESHOLDS_DB = (1.0, 2.0, 3.0)
MIN_RUN_H = 3

# 1. dry hours per link (all 5 min bins dry)
state = d.store("link_reference_state_5min").state.load().to_pandas()          # (time, cml_id)
dry_h = (state == 0).groupby(state.index.floor("1h")).all()

# hourly mean TL per good sub-link
good = d.status().query("status == 'good'").sublink_id.tolist()
tl = d.store("cml_1h").tl_mean.sel(sublink_id=good).load().to_pandas()       # (time, sublink_id)
cml_of = {s: int(s.split("_")[0]) for s in good}

rows = []
for s in good:
    dry = dry_h[cml_of[s]].reindex(tl.index, fill_value=False).to_numpy()
    x = tl[s].where(dry).dropna()
    if len(x) < 100:                       # too few dry hours to say anything
        continue
    level = x.median()                     # 2. usual dry level
    row = {"sublink_id": s, "dry_hours": len(x), "level_db": level}
    for th in THRESHOLDS_DB:
        exc = (tl[s] > level + th) & dry   # 3. excursion hours (dry only)
        runs = exc.groupby((exc != exc.shift()).cumsum()).sum()   # length of each True run
        row[f"share_dry_h_above_{th:g}dB"] = exc.sum() / dry.sum()
        row[f"episodes_{th:g}dB"] = int((runs >= MIN_RUN_H).sum())
    rows.append(row)

r = pd.DataFrame(rows)
r.to_csv(HERE / "dry_excursions.csv", index=False)
print(f"{len(r)} good sub-links with >= 100 dry hours")
for th in THRESHOLDS_DB:
    ep = r[f"episodes_{th:g}dB"]
    sh = r[f"share_dry_h_above_{th:g}dB"]
    print(f"> {th:g} dB above usual dry level: sub-links with >= 1 episode of >= {MIN_RUN_H} h: "
          f"{(ep > 0).sum()} ({(ep > 0).mean():.0%}); median share of dry hours: {sh.median():.2%}, "
          f"90th pct: {sh.quantile(0.9):.2%}")
print("10200:", r[r.sublink_id.str.startswith("10200_")].round(3).to_dict("records"))
