"""The one place that reads NetCDF metadata and the consolidated CML files with netCDF4.

Why netCDF4 and not xarray here: these helpers are called in loops over the 92 daily files (stage hashes, skip
checks) and on the 174 MB consolidated file (one or two variables). xarray opens every variable of every file it
touches, which is slow on this Windows machine, and dask is not installed. Bulk products are opened with xarray
elsewhere, one file at a time.
"""

from __future__ import annotations

from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd


def file_attrs(path: Path) -> dict:
    """Global attributes of a NetCDF file ({} if it is missing or unreadable)."""
    if not Path(path).exists():
        return {}
    try:
        with netCDF4.Dataset(path) as nc:
            return {k: nc.getncattr(k) for k in nc.ncattrs()}
    except OSError:
        return {}


def existing_stage_hash(path: Path) -> str | None:
    """stage_hash attribute of an existing output file, or None if the file is missing or unreadable."""
    return file_attrs(path).get("stage_hash")


def read_variables(path: Path, variables: list[str], time_slice: slice = slice(None)) -> tuple[pd.DatetimeIndex, dict, dict]:
    """(time, {var: array}, {coord: array}) of a (time, ...) NetCDF file, read with netCDF4 (no masking).

    Only the requested variables and the time slice are read. coords holds cml_id and sublink_id if present.
    """
    with netCDF4.Dataset(path) as nc:
        nc.set_auto_mask(False)
        tv = nc["time"]
        t = pd.DatetimeIndex(netCDF4.num2date(tv[time_slice], tv.units, only_use_cftime_datetimes=False,
                                              only_use_python_datetimes=True))
        data = {v: nc[v][time_slice] for v in variables}
        coords = {}
        if "cml_id" in nc.variables:
            coords["cml_id"] = nc["cml_id"][:].astype("int64")
        if "sublink_id" in nc.variables:
            coords["sublink_id"] = [str(s) for s in nc["sublink_id"][:]]
    return t, data, coords


def read_core_file(path: Path, variables: list[str]) -> tuple[pd.DatetimeIndex, np.ndarray, list[str], dict]:
    """(time, cml_ids, sublink_ids, {var: (time, sublink, cml) array}) of a consolidated CML core file."""
    t, data, coords = read_variables(path, variables)
    return t, coords["cml_id"], coords["sublink_id"], data
