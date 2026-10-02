"""Load the pipeline configuration (YAML) and compute its hash.

Why a hash: every output file stores the hash of the config that produced it.
When the config changes, the hash changes, so we can tell stale outputs from
current ones and re-run only what is needed.

All paths in the YAML are relative to the repository root, so the same config
works on every machine (Windows or Linux) as long as the folder layout is the same.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml

# src/cmlpipe/config.py -> parents[0] = cmlpipe, [1] = src, [2] = repository root.
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG = REPO_ROOT / "configs" / "openmrg.yaml"


def config_hash(data: dict) -> str:
    """Return the sha256 hex digest of the config as canonical JSON.

    Canonical means: keys sorted, no extra whitespace. Two configs with the same
    content therefore get the same hash even if the YAML is formatted differently.
    """
    canonical = json.dumps(data, sort_keys=True, separators=(",", ":"), ensure_ascii=True, default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


# Config keys each stage depends on (dotted paths into the YAML). A stage's hash covers ONLY these
# keys, so editing, say, the radar section does not make the CML files look stale.
# NO dataset path is part of any stage hash: moving the raw data or the output folder does not change what is
# inside a file. (Output locations were never hashed; raw-file paths were until step 8, see LEGACY_PATH_HASH below.)
STAGE_CONFIG_KEYS = {
    "ingest_cml": ["dataset", "cml", "encoding", "validation"],
    "ingest_gauges": ["dataset", "period", "reference", "city_gauges", "smhi_gauge", "reference_validation"],
    "ingest_met": ["dataset", "period", "reference", "smhi_gauge", "smhi_met", "reference_validation"],
    "ingest_radar": ["dataset", "period", "reference", "radar", "reference_validation"],
    "qc_stats": ["dataset", "qc.stats"],
    "qc_flags": ["dataset", "qc.context_hours", "qc.flags", "qc.blocking_flags", "qc.summary"],
    "qc_status": ["dataset", "qc.blocking_flags", "qc.status"],
    # metadata parsing keys too: they decide which end is site_0 and how polarization is spelled
    "geometry": ["dataset", "cml.direction_to_sublink", "cml.polarization_map",
                 "reference.lonlat_crs", "reference.metric_crs", "geometry"],
    "aggregate_cml": ["dataset", "cml.sampling_interval", "aggregate.time_convention", "aggregate.resolutions",
                      "aggregate.min_valid_fraction", "aggregate.blocking_flags", "aggregate.stats",
                      "aggregate.encoding"],
    "aggregate_gauges": ["dataset", "period", "aggregate.time_convention", "aggregate.gauges"],
    "aggregate_radar": ["dataset", "period", "aggregate.time_convention", "aggregate.radar"],
    "events": ["dataset", "period", "events", "aggregate.radar.min_dbz", "aggregate.time_convention"],
    "split": ["dataset", "split"],
    "event_links": ["dataset", "event_links", "events.resolution", "events.wet_threshold_mmh"],
    "significance": ["dataset", "event_significance", "events.resolution"],
    "inventory": ["dataset", "inventory", "aggregate.blocking_flags", "event_links.reference_gauge_max_m"],
    "quicklooks": ["dataset", "quicklooks", "event_significance"],
    "link_cards": ["dataset", "link_cards", "inventory.core_vars_1h"],
    "run_report": ["dataset"],
    "event_links_extra": ["dataset", "events.resolution"],
    "link_state": ["dataset", "link_state", "events.resolution"],
    "inventory_figures": ["dataset", "inventory.status", "inventory.tsl"],
    "release": ["dataset"],
}


# Before step 8 five stages also hashed the raw-file PATHS. Those outputs (all of L1, the geometry tables) are
# still valid after the data moved to data/raw/openmrg, so these stages hash the pre-move path strings as fixed
# constants: the live path never enters a hash, and the hashes of existing outputs stay the same (nothing is rebuilt
# because of a move). A new dataset simply has no entry here.
LEGACY_PATH_HASH = {
    "paths.raw_cml": "data/OpenMRG/cml/cml.nc",
    "paths.raw_cml_metadata": "data/OpenMRG/cml/cml_metadata.csv",
    "paths.raw_city_gauges": "data/OpenMRG/gauges/city/CityGauges-2015JJA.csv",
    "paths.raw_city_gauges_metadata": "data/OpenMRG/gauges/city/CityGauges-metadata.csv",
    "paths.raw_smhi": "data/OpenMRG/gauges/smhi/GbgA-71420-2015JJA.csv",
    "paths.raw_smhi_metadata": "data/OpenMRG/gauges/smhi/GbgA-71420-metadata.csv",
    "paths.raw_radar": "data/OpenMRG/radar/radar.nc",
}
STAGE_LEGACY_PATH_KEYS = {
    "ingest_cml": ["paths.raw_cml", "paths.raw_cml_metadata"],
    "ingest_gauges": ["paths.raw_city_gauges", "paths.raw_city_gauges_metadata", "paths.raw_smhi", "paths.raw_smhi_metadata"],
    "ingest_met": ["paths.raw_smhi", "paths.raw_smhi_metadata"],
    "ingest_radar": ["paths.raw_radar"],
    "geometry": ["paths.raw_cml_metadata"],
}


def _get_dotted(data: dict, key: str):
    """data["a"]["b"] for key "a.b". A missing key raises, so a typo in STAGE_CONFIG_KEYS is noticed."""
    node = data
    for part in key.split("."):
        if not isinstance(node, dict) or part not in node:
            raise KeyError(f"config key '{key}' not found (used in STAGE_CONFIG_KEYS)")
        node = node[part]
    return node


@dataclass
class Config:
    """The loaded configuration.

    Attributes
    ----------
    data : dict
        The YAML content, unchanged.
    repo_root : Path
        Folder that all relative paths are resolved against.
    config_file : Path
        The YAML file this config was read from.
    hash : str
        sha256 of the canonical JSON of `data`.
    """

    data: dict
    repo_root: Path
    config_file: Path
    hash: str

    # ---- short accessors for the config sections --------------------------
    def __getitem__(self, key):
        return self.data[key]

    def stage_hash(self, stage: str) -> str:
        """sha256 of the canonical JSON of only the config keys `stage` depends on (STAGE_CONFIG_KEYS).

        The stage name is part of the hashed content, so two stages that happen to use the same keys
        still get different hashes.
        """
        if stage not in STAGE_CONFIG_KEYS:
            raise KeyError(f"unknown stage {stage!r}; known: {sorted(STAGE_CONFIG_KEYS)}")
        subset = {key: _get_dotted(self.data, key) for key in STAGE_CONFIG_KEYS[stage]}
        if self.data.get("dataset", {}).get("name") == "openmrg":
            subset.update({k: LEGACY_PATH_HASH[k] for k in STAGE_LEGACY_PATH_KEYS.get(stage, [])})
        return config_hash({"stage": stage, "keys": subset})

    @property
    def dataset_name(self) -> str:
        """Dataset name from the config (e.g. 'openmrg')."""
        return self.data["dataset"]["name"]

    @property
    def sampling_interval(self) -> pd.Timedelta:
        """CML sampling step (10 s)."""
        return pd.Timedelta(self.data["cml"]["sampling_interval"])

    # ---- paths ------------------------------------------------------------
    def repo_path(self, relative: str) -> Path:
        """Resolve a path from the config relative to the repository root."""
        return (self.repo_root / relative).resolve()

    @property
    def raw_cml(self) -> Path:
        """Raw cml.nc path."""
        return self.repo_path(self.data["paths"]["raw_cml"])

    @property
    def raw_cml_metadata(self) -> Path:
        """Raw cml_metadata.csv path."""
        return self.repo_path(self.data["paths"]["raw_cml_metadata"])

    @property
    def output_root(self) -> Path:
        """Root of the work folder (data/work/<dataset>)."""
        return self.repo_path(self.data["paths"]["output_root"])

    @property
    def l1_cml_dir(self) -> Path:
        """Folder of the L1 CML day files."""
        return self.output_root / self.data["paths"]["l1_cml_dir"]

    def raw_path(self, key: str) -> Path:
        """A raw input file from the `paths` section, e.g. raw_path("raw_radar")."""
        return self.repo_path(self.data["paths"][key])

    @property
    def l1_gauges_dir(self) -> Path:
        """Folder of the L1 gauge files."""
        return self.output_root / self.data["paths"]["l1_gauges_dir"]

    @property
    def l1_met_dir(self) -> Path:
        """Folder of the L1 weather file."""
        return self.output_root / self.data["paths"]["l1_met_dir"]

    @property
    def l1_radar_dir(self) -> Path:
        """Folder of the L1 radar month files."""
        return self.output_root / self.data["paths"]["l1_radar_dir"]

    @property
    def period_start(self) -> pd.Timestamp:
        """First instant of the period (00:00 UTC of the first day)."""
        return pd.Timestamp(self.data["period"]["start"])

    @property
    def period_end(self) -> pd.Timestamp:
        """First instant AFTER the period (00:00 UTC of the day after the last day)."""
        return pd.Timestamp(self.data["period"]["end"]) + pd.Timedelta("1D")

    @property
    def tables_dir(self) -> Path:
        """Folder of the tables."""
        return self.output_root / self.data["paths"]["tables_dir"]

    @property
    def logs_dir(self) -> Path:
        """Folder of the run logs."""
        return self.output_root / self.data["paths"]["logs_dir"]

    @property
    def manifest_path(self) -> Path:
        """Path of run_manifest.json."""
        return self.output_root / self.data["paths"]["manifest"]

    def l1_cml_day_path(self, day: pd.Timestamp) -> Path:
        """Output file of one UTC day, e.g. .../cml_10s_2015-07-28.nc."""
        interval = self.data["cml"]["sampling_interval"]
        return self.l1_cml_dir / f"cml_{interval}_{pd.Timestamp(day):%Y-%m-%d}.nc"

    def period_days(self) -> list[pd.Timestamp]:
        """All UTC days of the configured period (both ends included)."""
        start = pd.Timestamp(self.data["period"]["start"])
        end = pd.Timestamp(self.data["period"]["end"])
        return list(pd.date_range(start, end, freq="D"))


def load_config(config_file: str | Path, repo_root: str | Path | None = None) -> Config:
    """Read a YAML config file.

    Parameters
    ----------
    config_file : path to the YAML file.
    repo_root : folder that relative paths are resolved against. Default: the
        repository that holds this code. Tests pass a temporary folder here.
    """
    config_file = Path(config_file).resolve()
    with open(config_file, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f)
    root = Path(repo_root).resolve() if repo_root is not None else REPO_ROOT
    return Config(data=data, repo_root=root, config_file=config_file, hash=config_hash(data))
