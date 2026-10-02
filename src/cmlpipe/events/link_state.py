"""Per-link reference state at 5 min (stage link_state): wet / dry / unknown for every link and bin.

From the REFERENCES only (radar along the path, gauges near the path), never from CML, so it can label CML data
for training without the labels depending on the CML itself. Rules (config link_state):

  wet     : R_max >= radar_wet_mmh in the bin (wettest radar cell under the path, floored at 0 dBZ),
            or a gauge within gauge_max_dist_m of the path has >= gauge_wet_mm in the bin
  dry     : not wet, and R_max == 0 in the bin and in every bin within +-margin_min, and every gauge within
            gauge_max_dist_m is 0 (and not missing) over the same window (a link without such a gauge passes this part)
  unknown : everything else (radar missing in the window, light rain or rain in the margin, a gauge missing)

Stored as int8 state (0 dry, 1 wet, -1 unknown) and a uint8 reason mask:
  1 radar_wet   (R_max >= threshold in the bin)
  2 gauge_wet   (a near gauge >= threshold in the bin)
  4 radar_missing (R_max NaN in the bin or within the margin window)
  8 margin      (not wet in the bin, but R_max > 0 or a near gauge > 0 somewhere in the +-margin window, bin included)

Gauges: City gauges at their native 5 min (sums of 1 min). The SMHI gauge is 15 min native: its 15 min amount is
used for each of the three 5 min bins of that quarter hour (a state label, not a disaggregated amount). At the period
edges the window is truncated to the available bins.
Output: L3_products/events/link_reference_state_5min.nc, dims (time, cml_id).
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from cmlpipe import aggregate
from cmlpipe.aggregate import stage as aggregate_stage
from cmlpipe.config import Config
from cmlpipe.events import stage as events_stage
from cmlpipe.io.netcdf import file_attrs
from cmlpipe.provenance import provenance_attrs

STATE = {"dry": 0, "wet": 1, "unknown": -1}
REASON = {"radar_wet": 1, "gauge_wet": 2, "radar_missing": 4, "margin": 8}


def state_path(cfg: Config) -> Path:
    """Where the per-link reference state is written."""
    return cfg.output_root / cfg["paths"]["l3_events_dir"] / "link_reference_state_5min.nc"


def window_any(mask: np.ndarray, half: int) -> np.ndarray:
    """True where `mask` is True anywhere within +-half bins along axis 0 (truncated at the edges)."""
    c = np.concatenate([np.zeros((1,) + mask.shape[1:], int), np.cumsum(mask, axis=0)])
    n = mask.shape[0]
    hi = np.minimum(np.arange(n) + half + 1, n)
    lo = np.maximum(np.arange(n) - half, 0)
    return (c[hi] - c[lo]) > 0


def classify_links(r_max: np.ndarray, gauge_amounts: np.ndarray, gauge_index: list[list[int]], lcfg: dict,
                   step: pd.Timedelta) -> tuple[np.ndarray, np.ndarray]:
    """(state int8, reasons uint8), both (time, link).

    r_max         : (time, link) R_max along the path (NaN = no valid radar)
    gauge_amounts : (time, gauge) 5 min amounts (NaN = missing)
    gauge_index   : per link, the columns of gauge_amounts within gauge_max_dist_m of its path
    """
    half = int(round(pd.Timedelta(minutes=lcfg["margin_min"]) / step))
    t, n = r_max.shape
    radar_missing = np.isnan(r_max)
    rm = np.nan_to_num(r_max, nan=0.0)
    radar_wet = rm >= lcfg["radar_wet_mmh"]
    radar_nonzero = rm > 0
    g_wet = np.zeros((t, n), bool)
    g_nonzero = np.zeros((t, n), bool)
    g_missing = np.zeros((t, n), bool)
    for j, idx in enumerate(gauge_index):
        if idx:
            a = gauge_amounts[:, idx]
            g_wet[:, j] = (np.nan_to_num(a) >= lcfg["gauge_wet_mm"] - 1e-6).any(axis=1)
            g_nonzero[:, j] = (np.nan_to_num(a) > 0).any(axis=1)
            g_missing[:, j] = np.isnan(a).any(axis=1)
    wet = radar_wet | g_wet
    miss_win = window_any(radar_missing, half)
    rain_win = window_any(radar_nonzero | g_nonzero, half)
    gmiss_win = window_any(g_missing, half)
    dry = ~wet & ~miss_win & ~rain_win & ~gmiss_win
    state = np.where(wet, STATE["wet"], np.where(dry, STATE["dry"], STATE["unknown"])).astype("int8")
    reasons = (radar_wet * REASON["radar_wet"] + g_wet * REASON["gauge_wet"] + miss_win * REASON["radar_missing"]
               + (~wet & rain_win) * REASON["margin"]).astype("uint8")
    return state, reasons


def build(cfg: Config) -> xr.Dataset:
    """The per-link reference state dataset (reads radar_along_links_5min, the City 5 min and SMHI 15 min gauge
    products and link_gauge_matches.csv)."""
    lcfg, res = cfg["link_state"], cfg["events"]["resolution"]
    step = pd.Timedelta(res)
    with xr.open_dataset(aggregate_stage.radar_path(cfg, "along_links", res)) as r:
        r_max = r.R_max.transpose("time", "cml_id").values.astype("float64")
        times, cids = r.time.values, r.cml_id.values
    city = events_stage.read_city(cfg, res)
    if not np.array_equal(city.time.values, times):
        raise ValueError("City gauge product and radar product have different time axes")
    smhi15 = events_stage.read_smhi(cfg)
    smhi5 = smhi15.reindex(pd.DatetimeIndex(times).floor("15min")).to_numpy()  # each 5 min bin gets its quarter hour
    amounts = np.column_stack([city.amount_mm.transpose("time", "id").values.astype("float64"), smhi5])
    names = [str(i) for i in city.id.values] + [cfg["smhi_gauge"]["id"]]
    m = pd.read_csv(cfg.tables_dir / "link_gauge_matches.csv")
    near = m[m.distance_to_path_m <= lcfg["gauge_max_dist_m"]]
    gidx = [[names.index(g) for g in near[near.cml_id == c].gauge_id] for c in cids]
    state, reasons = classify_links(r_max, amounts, gidx, lcfg, step)
    ds = xr.Dataset(
        {"state": (("time", "cml_id"), state, {"long_name": "per-link reference state (references only, never CML)",
                                               "flag_values": np.array(list(STATE.values()), "int8"),
                                               "flag_meanings": " ".join(STATE)}),
         "reasons": (("time", "cml_id"), reasons, {"long_name": "why the state was given (bitmask)",
                                                   "flag_masks": np.array(list(REASON.values()), "uint8"),
                                                   "flag_meanings": " ".join(REASON)})},
        coords={"time": times, "cml_id": cids,
                "n_gauges_near": ("cml_id", np.array([len(g) for g in gidx], "int16"),
                                  {"long_name": f"gauges within {lcfg['gauge_max_dist_m']} m of the path"})})
    ds.attrs = {**aggregate.time_attrs(res, 1, "per 5 min bin of the radar / gauge products"),
                **provenance_attrs(cfg, "link_state"), "title": "Per-link reference state (wet / dry / unknown), 5 min",
                "rules": json.dumps(lcfg), "definition": "see cmlpipe.events.link_state module docstring"}
    return ds


def run_link_state(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """Stage link_state: write the per-link reference state; skipped when the stage hash and inputs are unchanged."""
    t0 = time.perf_counter()
    res = cfg["events"]["resolution"]
    h = hashlib.sha256()
    for p in (aggregate_stage.radar_path(cfg, "along_links", res), aggregate_stage.gauge_path(cfg, "city", res),
              aggregate_stage.gauge_path(cfg, "smhi", "15min")):
        a = file_attrs(p)
        h.update(f"{p.name}{a.get('stage_hash')}{a.get('created_utc')}".encode())
    h.update((cfg.tables_dir / "link_gauge_matches.csv").read_bytes())
    sig = h.hexdigest()
    out = state_path(cfg)
    a = file_attrs(out)
    if not force and a.get("stage_hash") == cfg.stage_hash("link_state") and a.get("input_signature") == sig:
        log.info("link_state: up to date, skipped (use --force to rebuild)")
        records.append({"step": "link_state", "status": "skipped_same_config"})
        return
    ds = build(cfg)
    ds.attrs["input_signature"] = sig
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    ds.to_netcdf(tmp, encoding={"time": aggregate_stage.TIME_ENCODING, "state": {"zlib": True, "complevel": 4},
                                "reasons": {"zlib": True, "complevel": 4}})
    tmp.replace(out)
    st = ds.state.values
    shares = {k: round(float((st == v).mean()), 4) for k, v in STATE.items()}
    rec = {"step": "link_state", "status": "written", "shares": shares, "runtime_s": round(time.perf_counter() - t0, 2)}
    log.info("link_state: wrote %s, shares %s (%.1f s)", out.name, shares, rec["runtime_s"])
    records.append(rec)
