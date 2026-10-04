"""Integration test on the real OpenMRG data for 2015-07-28.

Skipped automatically when data/raw/openmrg is not on this machine. Outputs go to a
temporary folder, never to data/work/.
"""

import netCDF4
import numpy as np
import pandas as pd
import pytest
import xarray as xr

from cmlpipe.io import schema
from cmlpipe.io import openmrg
from cmlpipe.config import load_config
from conftest import REPO_ROOT, make_config

RAW_CML = REPO_ROOT / "data" / "raw" / "openmrg" / "cml" / "cml.nc"
RAW_META = REPO_ROOT / "data" / "raw" / "openmrg" / "cml" / "cml_metadata.csv"

pytestmark = [pytest.mark.data, pytest.mark.skipif(not (RAW_CML.exists() and RAW_META.exists()),
                                reason="real OpenMRG data not available")]


@pytest.fixture(scope="module")
def real(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("real")
    cfg = load_config(make_config(tmp, RAW_CML, RAW_META, tmp / "out"), repo_root=REPO_ROOT)
    ds = openmrg.load_cml_day(cfg, "2015-07-28")
    return cfg, ds


def test_shape_and_validation(real):
    cfg, ds = real
    assert ds.sizes == {"time": 8640, "sublink_id": 2, "cml_id": 364}
    schema.validate_cml_l1(ds, cfg["validation"], cfg.sampling_interval)
    # every raw sub-link is present exactly once
    ids = ds.raw_sublink_id.values.ravel()
    assert sorted(ids.tolist()) == list(range(1, 729))


def test_torslanda_link_matches_raw(real):
    """Link 10227 (the Torslanda case study): values equal a direct netCDF4 read of the raw file."""
    cfg, ds = real
    link = ds.sel(cml_id=10227)
    assert sorted(link.raw_sublink_id.values.tolist()) == [115, 118]
    with netCDF4.Dataset(RAW_CML) as nc:
        t = np.asarray(nc.variables["time"][:], dtype="int64")
        i0 = int(np.searchsorted(t, pd.Timestamp("2015-07-28").value // 10**9))
        for sub in ("sublink_1", "sublink_2"):
            raw_id = int(link.raw_sublink_id.sel(sublink_id=sub))
            raw = np.ma.filled(nc.variables["rsl"][i0:i0 + 8640, raw_id - 1].astype("float32"), np.nan)
            np.testing.assert_array_equal(link.rsl.sel(sublink_id=sub).values, raw)
    # The paper reports RSL fell to about -73.6 dBm at the peak of the event.
    assert float(link.rsl.min()) == pytest.approx(-73.6, abs=0.05)


def test_real_day_roundtrip(real):
    cfg, ds = real
    info = openmrg.write_l1_day(ds, cfg.l1_cml_day_path(pd.Timestamp("2015-07-28")), cfg)
    assert info["roundtrip"]["rsl"]["bit_identical"] and info["roundtrip"]["tsl"]["bit_identical"]
    with xr.open_dataset(info["path"]) as back:
        back.load()
        schema.validate_cml_l1(back, cfg["validation"], cfg.sampling_interval)
