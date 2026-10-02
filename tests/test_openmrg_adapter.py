"""Tests of the OpenMRG adapter on the tiny synthetic raw dataset (see conftest.py)."""

import json

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from cmlpipe.io import schema
from cmlpipe.io import openmrg
from conftest import DAY, raw_values

SEC = pd.Timedelta("1s")


def at(ds, seconds):
    """Select the time step `seconds` after midnight."""
    return ds.sel(time=DAY + seconds * SEC)


def test_metadata_mapping(synthetic_cfg):
    meta = openmrg.load_cml_metadata(synthetic_cfg).set_index(["cml_id", "sublink_id"])
    # A -> sublink_1, B -> sublink_2, raw Sublink kept for traceability
    assert meta.loc[(100, "sublink_1"), "raw_sublink_id"] == 2
    assert meta.loc[(100, "sublink_1"), "direction"] == "A"
    assert meta.loc[(100, "sublink_2"), "raw_sublink_id"] == 1
    assert meta.loc[(100, "sublink_2"), "direction"] == "B"
    # Link 200 has only A: sublink_2 is a placeholder
    assert meta.loc[(200, "sublink_2"), "raw_sublink_id"] == schema.MISSING_SUBLINK
    assert np.isnan(meta.loc[(200, "sublink_2"), "frequency"])
    # Link 400 has only B: sublink_1 is a placeholder
    assert meta.loc[(400, "sublink_1"), "raw_sublink_id"] == schema.MISSING_SUBLINK
    assert len(meta) == 8  # 4 links x 2 sub-links


def test_polarization_normalized(synthetic_day):
    pol = synthetic_day.polarization
    assert pol.sel(cml_id=100).values.tolist() == ["v", "v"]
    assert pol.sel(cml_id=200, sublink_id="sublink_1").item() == "v"   # raw "V"
    assert pol.sel(cml_id=300).values.tolist() == ["h", "h"]          # raw "horizontal " / "Horizontal"
    assert pol.sel(cml_id=200, sublink_id="sublink_2").item() == ""    # absent sub-link


def test_sites_follow_sublink_1_transmitter(synthetic_day):
    ds = synthetic_day
    # Link 100: site_0 = Near of A = Far of B
    l100 = ds.sel(cml_id=100)
    assert (l100.site_0_lat.item(), l100.site_0_lon.item()) == (57.0, 12.0)
    assert (l100.site_1_lat.item(), l100.site_1_lon.item()) == (57.1, 12.1)
    assert (l100.site_0_antenna_diameter.item(), l100.site_1_antenna_diameter.item()) == (0.3, 0.6)
    # Link 400 has only B: site_0 must be B's Far (= the transmitter of the missing A)
    l400 = ds.sel(cml_id=400)
    assert (l400.site_0_lat.item(), l400.site_0_lon.item()) == (57.2, 12.2)
    assert (l400.site_0_antenna_diameter.item(), l400.site_1_antenna_diameter.item()) == (0.2, 0.3)
    # tx/rx sites
    assert ds.tx_site.sel(sublink_id="sublink_1").values.tolist() == [0] * 4
    assert ds.tx_site.sel(sublink_id="sublink_2").values.tolist() == [1] * 4
    assert (ds.rx_site == 1 - ds.tx_site).all()


def test_shape_and_grid(synthetic_day):
    ds = synthetic_day
    assert ds.rsl.dims == schema.DIMS
    assert ds.sizes == {"time": 8640, "sublink_id": 2, "cml_id": 4}
    assert ds.cml_id.values.tolist() == [100, 200, 300, 400]
    assert ds.time.values[0] == np.datetime64("2015-07-28T00:00:00")
    assert ds.time.values[-1] == np.datetime64("2015-07-28T23:59:50")
    assert (np.diff(ds.time.values) == np.timedelta64(10, "s")).all()
    # In-day raw steps: 00:00:00, 00:00:10 (twice), 00:00:30, 00:00:40 -> 4 unique; 8636 inserted
    assert ds.attrs["n_raw_steps_in_day"] == 5
    assert ds.attrs["n_duplicates_removed"] == 1
    assert ds.attrs["n_duplicates_with_different_values"] == 1
    assert ds.attrs["n_steps_inserted"] == 8636
    assert ds.attrs["n_offgrid_steps_dropped"] == 0


def test_values_placed_correctly(synthetic_day):
    ds = synthetic_day
    rsl_raw, tsl_raw = raw_values()
    # raw row 4 = 00:00:30; Link 100 sublink_1 = raw Sublink 2 = raw column 1
    assert at(ds, 30).rsl.sel(cml_id=100, sublink_id="sublink_1").item() == rsl_raw[4, 1]
    # Link 300 sublink_2 = raw Sublink 4 = raw column 3
    assert at(ds, 40).tsl.sel(cml_id=300, sublink_id="sublink_2").item() == tsl_raw[5, 3]
    # duplicate 00:00:10: the FIRST copy (raw row 2) is kept, not -99.9 from row 3
    assert at(ds, 10).rsl.sel(cml_id=100, sublink_id="sublink_1").item() == rsl_raw[2, 1]


def test_sentinel_and_missing_become_nan(synthetic_day):
    ds = synthetic_day
    # raw Sublink 1 (Link 100, B) at 00:00:00 held 1e10
    assert np.isnan(at(ds, 0).rsl.sel(cml_id=100, sublink_id="sublink_2").item())
    assert ds.attrs["n_sentinel_rsl"] == 1
    assert np.nanmax(np.abs(ds.rsl.values)) < 1e3
    # 00:00:20 was missing in the raw file: inserted, all NaN
    assert at(ds, 20).rsl.isnull().all()
    # absent sub-links are all NaN
    assert ds.rsl.sel(cml_id=200, sublink_id="sublink_2").isnull().all()
    assert ds.rsl.sel(cml_id=400, sublink_id="sublink_1").isnull().all()
    # previous / next day rows are not included: only 4 time steps have any data
    assert int(ds.rsl.notnull().any(dim=["sublink_id", "cml_id"]).sum()) == 4


def test_tl(synthetic_day):
    ds = synthetic_day
    np.testing.assert_array_equal(ds.tl.values, ds.tsl.values - ds.rsl.values)


def test_validator_accepts_adapter_output(synthetic_day, synthetic_cfg):
    rep = schema.validate_cml_l1(synthetic_day, synthetic_cfg["validation"], synthetic_cfg.sampling_interval)
    assert rep.ok, rep.errors


def test_write_roundtrip_is_lossless(synthetic_day, synthetic_cfg):
    path = synthetic_cfg.l1_cml_day_path(DAY)
    info = openmrg.write_l1_day(synthetic_day, path, synthetic_cfg)
    assert path.exists() and not path.with_name(path.name + ".tmp").exists()
    for name in ("rsl", "tsl"):
        assert info["roundtrip"][name]["bit_identical"]
    with xr.open_dataset(path) as back:
        back.load()
        # stored as int16 x 0.1 on disk
        assert back.rsl.encoding["dtype"] == np.dtype("int16")
        assert back.rsl.encoding["scale_factor"] == pytest.approx(0.1)
        np.testing.assert_array_equal(np.isnan(back.rsl.values), np.isnan(synthetic_day.rsl.values))
        np.testing.assert_allclose(back.rsl.values, synthetic_day.rsl.values, rtol=0, atol=1e-5)
        assert back.raw_sublink_id.sel(cml_id=100).values.tolist() == [2, 1]
        assert back.attrs["config_hash"] == synthetic_cfg.hash
        assert back.attrs["source_doi"] == "10.5281/zenodo.7107689"
        for key in ("code_version", "created_utc", "time_meaning"):
            assert back.attrs[key]
        schema.validate_cml_l1(back, synthetic_cfg["validation"], synthetic_cfg.sampling_interval)


def test_write_refuses_lossy_encoding(synthetic_day, synthetic_cfg):
    """A value that is not a multiple of 0.1 dB must make the writer fail, and leave no file."""
    ds = synthetic_day.copy(deep=True)
    ds["rsl"][3, 0, 0] = np.float32(-42.25)
    path = synthetic_cfg.l1_cml_day_path(DAY)
    with pytest.raises(ValueError, match="lossy"):
        openmrg.write_l1_day(ds, path, synthetic_cfg)
    assert not path.exists()
    assert not path.with_name(path.name + ".tmp").exists()


def test_cli_run_skip_and_manifest(synthetic_cfg):
    """run.py writes the day, skips it on the second run (same config hash), rebuilds it with --force."""
    from cmlpipe import cli as run

    args = ["--config", str(synthetic_cfg.config_file), "--steps", "ingest_cml", "--days", "2015-07-28"]
    assert run.main(args) == 0
    assert run.main(args) == 0
    assert run.main(args + ["--force"]) == 0

    runs = json.loads(synthetic_cfg.manifest_path.read_text(encoding="utf-8"))
    assert [r["days"][0]["status"] for r in runs] == ["written", "skipped_same_config", "written"]
    day = runs[0]["days"][0]
    assert runs[0]["config_hash"] == synthetic_cfg.hash
    assert day["n_time"] == 8640 and day["n_duplicates_removed"] == 1 and day["n_steps_inserted"] == 8636
    assert day["n_sublinks_all_nan_rsl"] == 0  # absent sub-links are not counted
    assert (synthetic_cfg.tables_dir / "sublink_id_map.csv").exists()
    id_map = pd.read_csv(synthetic_cfg.tables_dir / "sublink_id_map.csv")
    assert len(id_map) == 8
    assert list(synthetic_cfg.logs_dir.glob("*.log"))
