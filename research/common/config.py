"""Read research/config.yaml: every research setting, one section per module."""
from pathlib import Path

import yaml

CONFIG_FILE = Path(__file__).resolve().parents[1] / "config.yaml"


def load_config():
    """The whole config as a dict, e.g. load_config()["wet_antenna"]["freeze"]."""
    with open(CONFIG_FILE, encoding="utf-8") as fh:
        return yaml.safe_load(fh)
