"""Per-stage config hashes: a stage's hash changes only when a config key IT uses changes."""

import copy

import pytest

from cmlpipe.config import STAGE_CONFIG_KEYS, Config, config_hash, load_config
from conftest import REAL_CONFIG


def with_change(cfg: Config, *path, value) -> Config:
    """A copy of cfg with one key changed (path = nested keys)."""
    data = copy.deepcopy(cfg.data)
    node = data
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = value
    return Config(data=data, repo_root=cfg.repo_root, config_file=cfg.config_file, hash=config_hash(data))


@pytest.fixture
def cfg():
    return load_config(REAL_CONFIG)


def test_every_stage_has_a_hash(cfg):
    hashes = {s: cfg.stage_hash(s) for s in STAGE_CONFIG_KEYS}
    assert len(set(hashes.values())) == len(hashes)  # the stage name is hashed too


def test_unrelated_key_does_not_change_other_stages(cfg):
    changed = with_change(cfg, "radar", "zr_fallback", "a", value=300.0)
    assert changed.hash != cfg.hash  # the global hash does change ...
    assert changed.stage_hash("ingest_radar") != cfg.stage_hash("ingest_radar")  # ... and so does radar's
    for stage in ("ingest_cml", "ingest_gauges", "ingest_met", "qc_stats", "qc_flags", "qc_status", "geometry"):
        assert changed.stage_hash(stage) == cfg.stage_hash(stage), stage


def test_used_key_changes_its_stage(cfg):
    assert with_change(cfg, "cml", "raw_missing_value", value=-999.0).stage_hash("ingest_cml") != cfg.stage_hash("ingest_cml")
    assert with_change(cfg, "city_gauges", "height_m", value=2.0).stage_hash("ingest_gauges") != cfg.stage_hash("ingest_gauges")
    assert with_change(cfg, "qc", "flags", "spike", "min_abs_db", value=99.0).stage_hash("qc_flags") != cfg.stage_hash("qc_flags")
    # a qc_flags threshold does not touch the L1 stages
    assert with_change(cfg, "qc", "flags", "spike", "min_abs_db", value=99.0).stage_hash("ingest_cml") == cfg.stage_hash("ingest_cml")


def test_output_folder_is_not_part_of_the_hash(cfg):
    moved = with_change(cfg, "paths", "output_root", value="somewhere/else")
    assert all(moved.stage_hash(s) == cfg.stage_hash(s) for s in STAGE_CONFIG_KEYS)


def test_missing_key_raises(cfg):
    broken = copy.deepcopy(cfg.data)
    del broken["radar"]
    with pytest.raises(KeyError):
        Config(data=broken, repo_root=cfg.repo_root, config_file=cfg.config_file, hash="x").stage_hash("ingest_radar")
