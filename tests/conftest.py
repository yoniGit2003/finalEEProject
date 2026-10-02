"""Shared test fixtures: a tiny synthetic OpenMRG-like raw dataset.

The synthetic raw files copy the REAL raw layout (see
docs/data_notes/openmrg_cml_inspection.md): cml.nc with rsl/tsl float32 on
(time, sublink), missing value 1e10, time as int seconds since 1970; and a
metadata CSV with the same column names. It contains on purpose:
  - a 1e10 sentinel,
  - a missing timestamp (00:00:20),
  - a duplicated timestamp (00:00:10, second copy has a different value),
  - timestamps of the previous and next day (must be ignored),
  - Link 200 with only direction A, Link 400 with only direction B,
  - polarization spelled "Vertical", "V" and "horizontal ",
  - sub-link 9 in cml.nc without metadata (must be dropped),
  - Link 100 whose A sub-link has a higher number than its B sub-link.
"""

import copy
import os
import sys
import tempfile
from pathlib import Path

# netCDF-C 4.9 on Windows scans the HDF5 filter-plugin folder on EVERY Variable.filters() call, which xarray makes
# for each variable it opens: about 90 ms per variable, so a tiny test file costs seconds. The tests (and the
# pipeline) only use zlib, which is built into HDF5, so point the plugin path at a folder that does not exist.
# HDF5 reads this once when the library loads, so it must be set before netCDF4 is imported.
os.environ.setdefault("HDF5_PLUGIN_PATH", str(Path(tempfile.gettempdir()) / "cmlpipe-tests-no-hdf5-plugins"))

import netCDF4  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
import yaml  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
sys.path.insert(0, str(SRC_DIR))  # also works without `pip install -e .`

from cmlpipe.config import load_config  # noqa: E402

REAL_CONFIG = REPO_ROOT / "configs" / "openmrg.yaml"
DAY = pd.Timestamp("2015-07-28")
SENTINEL = np.float32(1e10)

# Raw sub-link numbers stored in the synthetic cml.nc (columns, in this order).
RAW_SUBLINKS = [1, 2, 3, 4, 5, 6, 9]
# Seconds after DAY 00:00 of each raw row. -10 and 86400 belong to other days.
RAW_OFFSETS = [-10, 0, 10, 10, 30, 40, 86400]

METADATA = pd.DataFrame(
    [
        # Link, Dir, Sublink, NearLat, NearLon, FarLat, FarLon, f_GHz, Pol, L_km, NearDiam, FarDiam
        (100, "A", 2, 57.0, 12.0, 57.1, 12.1, 28.0, "Vertical", 1.5, 0.3, 0.6),
        (100, "B", 1, 57.1, 12.1, 57.0, 12.0, 29.0, "Vertical", 1.5, 0.6, 0.3),
        (200, "A", 3, 57.5, 11.5, 57.6, 11.6, 38.0, "V", 0.8, 0.2, 0.2),
        (300, "A", 5, 57.7, 11.9, 57.8, 11.8, 18.0, "horizontal ", 5.0, 0.6, 0.6),
        (300, "B", 4, 57.8, 11.8, 57.7, 11.9, 19.0, "Horizontal", 5.0, 0.6, 0.6),
        (400, "B", 6, 57.3, 12.3, 57.2, 12.2, 23.0, "Vertical", 2.0, 0.3, 0.2),
    ],
    columns=["Link", "Direction", "Sublink", "NearLatitude_DecDeg", "NearLongitude_DecDeg",
             "FarLatitude_DecDeg", "FarLongitude_DecDeg", "Frequency_GHz", "Polarization",
             "Length_km", "NearAntennaDiameter_m", "FarAntennaDiameter_m"],
)


def raw_values():
    """rsl, tsl arrays (time, sublink) of the synthetic raw file, on a 0.1 dB grid."""
    n_t, n_s = len(RAW_OFFSETS), len(RAW_SUBLINKS)
    t = np.arange(n_t)[:, None]
    s = np.arange(n_s)[None, :]
    rsl = (-400 - 3 * t - 10 * s).astype("float64") / 10  # e.g. -40.0, -40.3, ...
    tsl = (100 + 10 * s + 0 * t).astype("float64") / 10
    rsl = rsl.astype("float32")
    tsl = tsl.astype("float32")
    rsl[1, 0] = SENTINEL          # 00:00:00, sub-link 1: missing value
    rsl[3, 1] = np.float32(-99.9)  # second copy of 00:00:10 differs from the first
    return rsl, tsl


def write_raw_files(folder: Path):
    """Write the synthetic cml.nc and cml_metadata.csv into `folder`."""
    folder.mkdir(parents=True, exist_ok=True)
    rsl, tsl = raw_values()
    t0 = int(DAY.value // 10**9)
    with netCDF4.Dataset(folder / "cml.nc", "w", format="NETCDF4") as nc:
        nc.createDimension("time", len(RAW_OFFSETS))
        nc.createDimension("sublink", len(RAW_SUBLINKS))
        tv = nc.createVariable("time", "i4", ("time",))
        tv.setncatts({"units": "seconds since 1970-01-01 00:00:00", "calendar": "gregorian",
                      "missing_value": np.int32(-99)})
        tv[:] = np.array(RAW_OFFSETS) + t0
        sv = nc.createVariable("sublink", "i4", ("sublink",))
        sv[:] = RAW_SUBLINKS
        for name, values in (("rsl", rsl), ("tsl", tsl)):
            v = nc.createVariable(name, "f4", ("time", "sublink"), fill_value=SENTINEL)
            v.setncatts({"missing_value": SENTINEL, "units": "dBm"})
            v.set_auto_mask(False)  # write the sentinel as a plain number, like the real file
            v[:] = values
        nc.setncatts({"title": "synthetic OpenMRG-like", "version": "test", "doi": "test-doi"})
    METADATA.to_csv(folder / "cml_metadata.csv", index=False)


_REAL_CONFIG_DATA: dict | None = None


def real_config_data() -> dict:
    """A fresh copy of the real config as a dict (the YAML is parsed once per test session)."""
    global _REAL_CONFIG_DATA
    if _REAL_CONFIG_DATA is None:
        _REAL_CONFIG_DATA = yaml.safe_load(REAL_CONFIG.read_text(encoding="utf-8"))
    return copy.deepcopy(_REAL_CONFIG_DATA)


def make_config(tmp_path: Path, raw_cml: Path, raw_meta: Path, output_root: Path) -> Path:
    """Copy the real config, point its paths to the given (absolute) locations, save it in tmp_path."""
    data = real_config_data()
    data["paths"]["raw_cml"] = str(raw_cml)
    data["paths"]["raw_cml_metadata"] = str(raw_meta)
    data["paths"]["output_root"] = str(output_root)
    path = tmp_path / "test_config.yaml"
    path.write_text(yaml.safe_dump(data), encoding="utf-8")
    return path


@pytest.fixture(scope="session")
def synthetic_raw_dir(tmp_path_factory):
    """The synthetic raw files, written once per session (tests only read them)."""
    raw_dir = tmp_path_factory.mktemp("synthetic_raw")
    write_raw_files(raw_dir)
    return raw_dir


@pytest.fixture
def synthetic_cfg(tmp_path, synthetic_raw_dir):
    """Config whose raw files are the synthetic ones and whose outputs go to tmp_path/out (a new one per test,
    because tests change the config and write outputs)."""
    cfg_file = make_config(tmp_path, synthetic_raw_dir / "cml.nc", synthetic_raw_dir / "cml_metadata.csv",
                           tmp_path / "out")
    return load_config(cfg_file, repo_root=tmp_path)


@pytest.fixture(scope="session")
def _synthetic_day_once(tmp_path_factory, synthetic_raw_dir):
    tmp = tmp_path_factory.mktemp("synthetic_day")
    cfg = load_config(make_config(tmp, synthetic_raw_dir / "cml.nc", synthetic_raw_dir / "cml_metadata.csv",
                                  tmp / "out"), repo_root=tmp)
    from cmlpipe.io import openmrg
    return openmrg.load_cml_day(cfg, DAY)


@pytest.fixture
def synthetic_day(_synthetic_day_once):
    """The synthetic day loaded through the adapter (in memory): loaded once, a deep copy per test."""
    return _synthetic_day_once.copy(deep=True)


@pytest.fixture(scope="session")
def aggregated_world(tmp_path_factory, synthetic_raw_dir):
    """ONE tiny world (see test_aggregate.make_world) after one aggregation run, shared by every test that only reads
    its products; tests that change files copy it first (test_events.ev_world). Defined here, not in
    test_aggregate, so that all modules share the same instance.
    Returns (cfg, sha256 of the L1 / L2 inputs before the run, the same after the run)."""
    from test_aggregate import DAYS, digest, make_world, run_all
    from cmlpipe.qc import stage as qc_runner
    tmp = tmp_path_factory.mktemp("aggregated_world")
    cfg = make_world(load_config(make_config(tmp, synthetic_raw_dir / "cml.nc", synthetic_raw_dir / "cml_metadata.csv",
                                             tmp / "out"), repo_root=tmp))
    inputs = [cfg.l1_cml_day_path(d) for d in DAYS] + [qc_runner.l2_day_path(cfg, d) for d in DAYS]
    before = digest(inputs)
    run_all(cfg)
    return cfg, before, digest(inputs)
