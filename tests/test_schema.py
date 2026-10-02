"""Validator tests: a valid L1 dataset passes, deliberately broken ones fail."""

import numpy as np
import pytest

from cmlpipe.io import schema


@pytest.fixture
def valid(synthetic_day):
    """A small valid L1 dataset: the first 30 steps of the synthetic day."""
    return synthetic_day.isel(time=slice(0, 30))


def test_valid_dataset_passes(valid):
    rep = schema.validate_cml_l1(valid)
    assert rep.ok and rep.warnings == []


def test_wrong_dim_order_fails(valid):
    broken = valid.transpose("time", "cml_id", "sublink_id")
    with pytest.raises(schema.CMLValidationError, match="dims"):
        schema.validate_cml_l1(broken)


def test_missing_units_fails(valid):
    broken = valid.copy()
    broken["frequency"].attrs = {"long_name": "carrier frequency"}
    with pytest.raises(schema.CMLValidationError, match="frequency' units"):
        schema.validate_cml_l1(broken)


def test_non_monotonic_time_fails(valid):
    order = list(range(valid.sizes["time"]))
    order[1], order[2] = order[2], order[1]
    broken = valid.isel(time=order)
    with pytest.raises(schema.CMLValidationError, match="strictly increasing"):
        schema.validate_cml_l1(broken)


def test_irregular_time_fails(valid):
    broken = valid.isel(time=[0, 1, 3, 4])
    with pytest.raises(schema.CMLValidationError, match="not regular"):
        schema.validate_cml_l1(broken)


def test_frequency_in_mhz_fails(valid):
    freq = valid.frequency * 1000
    freq.attrs = valid.frequency.attrs  # keep units = "GHz" so only the range check can catch it
    broken = valid.assign_coords(frequency=freq)
    with pytest.raises(schema.CMLValidationError, match="MHz"):
        schema.validate_cml_l1(broken)


def test_length_in_metres_fails(valid):
    length = valid.length * 1000
    length.attrs = valid.length.attrs
    with pytest.raises(schema.CMLValidationError, match="length"):
        schema.validate_cml_l1(valid.assign_coords(length=length))


def test_bad_polarization_fails(valid):
    pol = valid.polarization.copy()
    pol[0, 0] = "Vertical"
    with pytest.raises(schema.CMLValidationError, match="polarization"):
        schema.validate_cml_l1(valid.assign_coords(polarization=pol))


def test_tl_mismatch_fails(valid):
    broken = valid.copy(deep=True)
    broken["tl"][4, 0, 0] = broken["tl"][4, 0, 0] + 1.0
    with pytest.raises(schema.CMLValidationError, match="tl != tsl - rsl"):
        schema.validate_cml_l1(broken)


def test_implausible_values_only_warn(valid):
    odd = valid.copy(deep=True)
    odd["rsl"][4, 0, 0] = -5.0
    odd["tl"][4, 0, 0] = odd["tsl"][4, 0, 0] - odd["rsl"][4, 0, 0]
    rep = schema.validate_cml_l1(odd)
    assert rep.ok
    assert len(rep.warnings) == 1 and "rsl" in rep.warnings[0]


def test_report_without_raising(valid):
    rep = schema.validate_cml_l1(valid.transpose("time", "cml_id", "sublink_id"), raise_on_error=False)
    assert not rep.ok and rep.errors
