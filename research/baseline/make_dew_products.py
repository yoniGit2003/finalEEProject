"""Dew products for later research code (documented in METHOD.md; thresholds in research/config.yaml, baseline.dew).

outputs/products/dew_risk_hourly.csv : one row per hour of JJA 2015 (SMHI GbgA): temperature, RH, dew point (Magnus),
                               dew-point depression, wind, dew_risk (RH >= 90 % or depression < 2 degC), solar
                               elevation at mid-hour and night (sun below the horizon in Gothenburg).
outputs/products/dry_jumps.csv       : the R1 C4 dry-day jump table for all train sub-links (good and degraded), with
                               both_directions, neighbours_jumping, weather and start / end in local time (UTC + 2 h).
Needs outputs/tables/c4_jumps.csv (python -m research.baseline.validation.c4_jumps).

Run: python -m research.baseline.make_dew_products
"""
from pathlib import Path
import sys
import pandas as pd
import yaml
from cmlpipe.api import CMLData

HERE = Path(__file__).resolve().parent
from research.baseline import baseline as bl                                                          # noqa: E402
from research.baseline import paths as P                                                # noqa: E402  (all file locations)
from research.common.save import save

REPO = P.REPO
CONFIG = yaml.safe_load(open(P.CONFIG_FILE, encoding="utf-8"))
DEW = CONFIG["baseline"]["dew"]
P.PRODUCTS.mkdir(parents=True, exist_ok=True)
LOCAL = pd.Timedelta(hours=2)                                                  # CEST = UTC + 2 h all summer

d = CMLData.open(REPO / "data" / "release" / CONFIG["release"])
met = d.store("met_1h").load()
h = pd.DataFrame({v: met[v].isel(gauge_id=0).to_series()
                  for v in ("temperature", "humidity", "wind_speed")}).rename_axis("time_utc")
h["dew_point"] = bl.dew_point(h.temperature, h.humidity)
h["dew_point_depression"] = h.temperature - h.dew_point
h["dew_risk"] = bl.dew_risk(h.humidity, h.dew_point_depression, DEW) & h.humidity.notna()
h["sun_elevation_deg"] = bl.sun_elevation_deg(h.index + pd.Timedelta("30min"), DEW["site_lat"], DEW["site_lon"])
h["night"] = h.sun_elevation_deg < DEW["night_sun_elev_deg"]
h.insert(0, "time_local", h.index + LOCAL)
save(h.round(3), P.DEW_RISK_HOURLY)

j = pd.read_csv(P.C4_JUMPS, parse_dates=["start", "end"])
j = j.rename(columns={"start": "start_utc", "end": "end_utc", "rh_mean": "mean_rh", "dpd_mean": "mean_dew_point_depression",
                      "t_mean": "mean_temperature", "wind_mean": "mean_wind_speed"})
j["start_local"], j["end_local"] = j.start_utc + LOCAL, j.end_utc + LOCAL
cols = ["cml_id", "sublink_id", "group", "is_degraded", "is_toggling", "sign", "start_utc", "end_utc", "start_local",
        "end_local", "duration_h", "median_excess_db", "peak_excess_db", "mean_rh", "rh_max",
        "mean_dew_point_depression", "mean_temperature", "mean_wind_speed", "hours_since_wet", "both_directions",
        "neighbours_jumping", "length_class", "frequency_band", "tsl_behaviour"]
save(j[cols].round(3), P.DRY_JUMPS, index=False)
print(f"dew_risk_hourly: {len(h)} hours, dew risk {h.dew_risk.mean():.1%}, night {h.night.mean():.1%}, "
      f"dew risk at night {h[h.night].dew_risk.mean():.1%} vs day {h[~h.night].dew_risk.mean():.1%}")
print(f"dry_jumps: {len(j)} jumps ({(j.sign == 'positive').sum()} positive) on {j.sublink_id.nunique()} sub-links")
