"""Write the CML time series of the release as Zarr: cml_10s (L1 + L2) and cml_<res> (L3), dims (time, sublink_id).

sublink_id is "<cml_id>_<1|2>", ordered by cml_id then sub-link, so the two directions of a link are neighbours.
Values are copied, never recomputed: the 10 s tsl / rsl / tl are the L1 int16 values (scale 0.1 dB, fill -32768) and
the flags the L2 uint16 bitmask; the L3 statistics keep their float32 / uint8 / uint16 types. Inputs are read with
netCDF4 one day at a time.
"""

from __future__ import annotations

from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

from cmlpipe.aggregate import stage as aggregate_stage
from cmlpipe.config import Config
from cmlpipe.qc import stage as qc_stage
from cmlpipe.release import zarrstore

L1_VARS = ("tsl", "rsl", "tl")


def sublink_coords(cfg: Config) -> tuple[dict, np.ndarray]:
    """Coordinates on sublink_id (from the first L1 day + link_geometry) and the column order (index into the
    flattened (cml, sublink) axis). Returns (coords, order)."""
    with netCDF4.Dataset(cfg.l1_cml_day_path(cfg.period_start)) as nc:
        cids = nc["cml_id"][:].astype("int64")
        subs = [str(s) for s in nc["sublink_id"][:]]
        freq = nc["frequency"][:]                          # (sublink, cml)
        pol = np.array(nc["polarization"][:]).astype(str)  # (sublink, cml)
        raw = nc["raw_sublink_id"][:]
        direction = np.array(nc["direction"][:]).astype(str)
    lg = pd.read_csv(cfg.tables_dir / "link_geometry.csv").set_index("cml_id").loc[cids]
    n_s = len(subs)
    labels, cml, sub_no, f, p, r, d, length = [], [], [], [], [], [], [], []
    for j, c in enumerate(cids):
        for k in range(n_s):
            labels.append(f"{c}_{k + 1}")
            cml.append(c)
            sub_no.append(k + 1)
            f.append(float(freq[k, j]))
            p.append(pol[k, j])
            r.append(int(raw[k, j]))
            d.append(direction[k, j])
            length.append(float(lg.length_m.iloc[j]))
    coords = {"sublink_id": np.array(labels), "cml_id": ("sublink_id", np.array(cml)),
              "sublink": ("sublink_id", np.array(sub_no, "int8")), "frequency_ghz": ("sublink_id", np.array(f)),
              "polarization": ("sublink_id", np.array(p)), "raw_sublink_id": ("sublink_id", np.array(r)),
              "direction": ("sublink_id", np.array(d)), "length_m": ("sublink_id", np.array(length))}
    # flattened source axis is (sublink, cml) -> column k * n_cml + j; we want (cml, sublink)
    order = np.array([k * len(cids) + j for j in range(len(cids)) for k in range(n_s)])
    return coords, order


def _flat(a: np.ndarray, order: np.ndarray) -> np.ndarray:
    """(time, sublink, cml) -> (time, sublink_id) in release column order."""
    return a.reshape(a.shape[0], -1)[:, order]


def write_cml_10s(cfg: Config, path: Path, chunks: tuple[int, int], days=None, attrs: dict | None = None) -> None:
    """cml_10s: tsl, rsl, tl (int16, scale 0.1, raw L1 values) and flags (uint16, L2) for the given days."""
    days = list(days) if days is not None else cfg.period_days()
    coords, order = sublink_coords(cfg)
    n_day = int(pd.Timedelta("1D") / cfg.sampling_interval)
    times = pd.date_range(days[0], periods=n_day * len(days), freq=cfg.sampling_interval)
    if len(days) > 1 and not (pd.DatetimeIndex(days).to_series().diff().dropna() == pd.Timedelta("1D")).all():
        raise ValueError("days must be consecutive")
    zarrstore.write_skeleton(path, times, coords, attrs or {})
    shape = (len(times), len(order))
    with netCDF4.Dataset(cfg.l1_cml_day_path(days[0])) as nc:
        v_attrs = {v: {k: nc[v].getncattr(k) for k in nc[v].ncattrs()} for v in L1_VARS}
    arrs = {v: zarrstore.add_array(path, v, ("time", "sublink_id"), shape, chunks, "int16",
                                   {"units": v_attrs[v]["units"], "long_name": v_attrs[v]["long_name"],
                                    "scale_factor": 0.1}, fill_value=-32768) for v in L1_VARS}
    from cmlpipe.qc.flags import flag_attrs
    fa = flag_attrs()
    arrs["flags"] = zarrstore.add_array(path, "flags", ("time", "sublink_id"), shape, chunks, "uint16",
                                        {"long_name": "L2 QC flags (bitmask)", "flag_masks": fa["flag_masks"],
                                         "flag_meanings": fa["flag_meanings"]}, fill_value=None)  # 0 = no flag, not missing
    for i, d in enumerate(days):
        a, b = i * n_day, (i + 1) * n_day
        with netCDF4.Dataset(cfg.l1_cml_day_path(d)) as nc:
            nc.set_auto_maskandscale(False)
            for v in L1_VARS:
                arrs[v][a:b, :] = _flat(nc[v][:], order)
        with netCDF4.Dataset(qc_stage.l2_day_path(cfg, d)) as nc:
            nc.set_auto_maskandscale(False)
            arrs["flags"][a:b, :] = _flat(nc["flags"][:], order)
    zarrstore.finalize(path)


def write_cml_res(cfg: Config, path: Path, res: str, chunks: tuple[int, int], attrs: dict | None = None) -> None:
    """cml_<res>: every L3 variable of the daily files, one variable at a time over the whole period."""
    days = cfg.period_days()
    coords, order = sublink_coords(cfg)
    files = [aggregate_stage.cml_path(cfg, res, d) for d in days]
    with netCDF4.Dataset(files[0]) as nc:
        nb = nc.dimensions["time"].size
        names = [v for v in nc.variables if nc[v].dimensions == ("time", "sublink_id", "cml_id")]
        meta = {v: (nc[v].dtype, {k: nc[v].getncattr(k) for k in nc[v].ncattrs()}) for v in names}
    times = pd.date_range(days[0], periods=nb * len(days), freq=res)
    zarrstore.write_skeleton(path, times, coords, attrs or {})
    for v in names:
        dtype, va = meta[v]
        fill = np.nan if np.dtype(dtype).kind == "f" else None  # counts / bitmasks: NO fill value (0 is data, not missing)
        keep = {k: va[k] for k in ("units", "long_name", "flag_masks", "flag_meanings", "comment") if k in va}
        arr = zarrstore.add_array(path, v, ("time", "sublink_id"), (len(times), len(order)), chunks, dtype, keep, fill)
        buf = np.empty((len(times), len(order)), dtype=dtype)
        for i, f in enumerate(files):
            with netCDF4.Dataset(f) as nc:
                nc.set_auto_maskandscale(False)
                buf[i * nb:(i + 1) * nb] = _flat(nc[v][:], order)
        arr[:] = buf
    zarrstore.finalize(path)
