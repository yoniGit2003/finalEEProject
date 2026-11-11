"""Library: the Part 1 TRAIN dataset, rebuilt with the Part 1 code (the Part 1 fit cache is not on this machine).

It calls `system_id.build_dataset` exactly as `system_id_s3.dataset` does (bands from `merged_bands`, links of every length of
the TRAIN link split, TRAIN events only) and saves the result in research/.cache/waa_ksurface/. Nothing of Part 1 is changed.
Only TRAIN link-events are built: VAL and TEST attenuation is not opened here.
"""
import pickle

import yaml

from cmlpipe.api import CMLData
from research.common.config import load_config
from research.wet_antenna import system_id as sid
from research.wet_antenna.event_factors import rescaled
from research.waa_ksurface import paths as P

import numpy as np
import pandas as pd

CFG1 = load_config()                                           # Part 1 settings (read only)
WA = CFG1["wet_antenna"]
DATA = WA["system_id_s2"]["data"]                              # the data rules of S2 (windows, exclusions, offset c)
S3 = WA["system_id_s3"]
LT = S3["length_test"]                                         # the Part 1 length test settings


def branch_config():
    return yaml.safe_load(P.CONFIG_FILE.read_text(encoding="utf-8"))


def open_release():
    return CMLData.open(P.REPO / "data" / "release" / CFG1["release"])


def train(d=None):
    """(D0, D1, le, drops, BANDS): the TRAIN dataset with the radar input (V0) and with radar / g_e (V1), the link-event table,
    the counts of the selection rules, the band names. Built once, then read from the cache."""
    d = d or open_release()
    bands_cfg = sid.merged_bands(d, DATA)[0]
    BANDS = list(bands_cfg)
    if not (P.DATA_TRAIN.exists() and P.LE_TRAIN.exists()):
        P.CACHE.mkdir(parents=True, exist_ok=True)
        D, le, drops = sid.build_dataset(d, {**DATA, "bands_ghz": bands_cfg, "max_length_km": S3["data"]["max_length_km_all"]}, event_split="train")
        sid.save_data(D, P.DATA_TRAIN)
        P.LE_TRAIN.write_bytes(pickle.dumps((le, drops)))
    le, drops = pickle.loads(P.LE_TRAIN.read_bytes())
    le = le.reset_index(drop=True)
    D0 = sid.load_data(P.DATA_TRAIN, BANDS)
    return D0, rescaled(D0, le, g_map()), le, drops, BANDS


def g_table():
    """g_e of the TRAIN events (Part 1, S4.3)."""
    g = pd.read_csv(P.P1_G_E)
    return g[g.set == "TRAIN"].reset_index(drop=True)


def g_map():
    g = g_table()
    return dict(zip(g.event_id, g.g_path))


BRANCH_BANDS = ["28 GHz", "38 GHz"]                            # 32 GHz is excluded by rule A1e (no eligible length-test cell)


def branch_train(d=None):
    """(D0, D1, le, BANDS): the TRAIN dataset of this branch = the Part 1 TRAIN dataset without the 32 GHz link-events.
    D0 = radar input, D1 = radar / g_e (V1). The band index of a link-event stays the Part 1 index (0 = 28, 2 = 38 GHz)."""
    D0, _, le, _, BANDS = train(d)
    keep = le.band_name.isin(BRANCH_BANDS).to_numpy()
    le_b = le[keep].reset_index(drop=True)
    D0b = sid.subset(D0, keep)
    return D0b, rescaled(D0b, le_b, g_map()), le_b, BANDS


def gauge_ratio_train(d=None):
    """(D, le, n_cells): the branch TRAIN link-events of the events WITH a gauge ratio (rain = radar / g_e), and the number of
    radar cells along the link of every link-event (data layer: links table, column n_cells)."""
    d = d or open_release()
    _, D1, le, _ = branch_train(d)
    g = g_table()
    keep = le.event_id.isin(set(g[g.path_used].event_id)).to_numpy()
    le_g = le[keep].reset_index(drop=True)
    n_cells = le_g.cml_id.map(d.links().set_index("cml_id").n_cells).to_numpy(float)
    return sid.subset(D1, keep), le_g, np.maximum(n_cells, 1.0)
