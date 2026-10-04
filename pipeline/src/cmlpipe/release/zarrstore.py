"""Zarr writing / reading helpers for the release (no dask).

Layout: Zarr format 2 (read by xarray, zarr-python 2 and 3, and other tools) with consolidated metadata, Blosc-zstd
compression. A store is written in two steps: xarray writes the small skeleton (time + coordinates), then every data
variable is created with the zarr API and filled block by block (a day, or one variable at a time), so the full
10 s arrays (794 880 x 728 per variable) are never held in memory. CF attributes (_ARRAY_DIMENSIONS, units,
scale_factor, fill value) make xarray decode the arrays exactly as the work-folder NetCDF files.

Read with:  xr.open_dataset(path, engine="zarr", chunks=None, consolidated=True)  (lazy, loads only what is sliced)
"""

from __future__ import annotations

from pathlib import Path

import numcodecs
import numpy as np
import pandas as pd
import xarray as xr
import zarr

COMPRESSOR = numcodecs.Blosc(cname="zstd", clevel=3, shuffle=numcodecs.Blosc.SHUFFLE)


def write_skeleton(path: Path, times: pd.DatetimeIndex, coords: dict, attrs: dict) -> None:
    """Create the store with the time axis and the coordinates (dict name -> (dims, values) or values)."""
    ds = xr.Dataset(coords={"time": pd.DatetimeIndex(times), **coords})
    ds.attrs = {k: (v if isinstance(v, (str, int, float, list)) else str(v)) for k, v in attrs.items()}
    ds["time"].encoding = {"units": "seconds since 1970-01-01 00:00:00", "calendar": "proleptic_gregorian", "dtype": "int64"}
    ds.to_zarr(path, mode="w", zarr_format=2, consolidated=False)


def add_array(path: Path, name: str, dims: tuple, shape: tuple, chunks: tuple, dtype, attrs: dict, fill_value=None):
    """Create an empty data array in an existing store and return it (fill it with arr[a:b] = ...)."""
    g = zarr.open_group(path, mode="a", zarr_format=2)
    if fill_value is None:
        fill_value = np.nan if np.dtype(dtype).kind == "f" else None
    arr = g.create_array(name, shape=shape, chunks=chunks, dtype=dtype, fill_value=fill_value, compressors=COMPRESSOR,
                         overwrite=True)
    clean = {k: (v.tolist() if isinstance(v, np.ndarray) else (v.item() if isinstance(v, np.generic) else v))
             for k, v in attrs.items() if k not in ("_FillValue", "missing_value", "coordinates", "dtype")}
    arr.attrs.update({"_ARRAY_DIMENSIONS": list(dims), **clean})
    return arr


def finalize(path: Path) -> None:
    """Consolidate the metadata (one small read opens the whole store)."""
    zarr.consolidate_metadata(str(path), zarr_format=2)


def open_store(path: Path) -> xr.Dataset:
    """Lazy xarray view of a store (no dask; slicing loads only the chunks needed)."""
    return xr.open_dataset(path, engine="zarr", chunks=None, consolidated=True)


def store_size(path: Path) -> int:
    """Bytes on disk."""
    return sum(f.stat().st_size for f in Path(path).rglob("*") if f.is_file())
