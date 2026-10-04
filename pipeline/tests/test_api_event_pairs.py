"""CMLData.event_pairs on a tiny synthetic release (three tables, two events x four links)."""

import json

import pandas as pd
import pytest

from cmlpipe.api import CMLData


@pytest.fixture(scope="module")
def rel(tmp_path_factory):
    """A release folder holding only the tables event_pairs reads (parquet + a minimal catalog.json)."""
    root = tmp_path_factory.mktemp("release")
    ev = ["E1"] * 4 + ["E2"] * 4
    links = [1, 2, 3, 4] * 2
    event_links = pd.DataFrame({
        "event_id": ev, "cml_id": links,
        # link 2: no usable sub-link; link 4: only sub-link 2 usable
        "status_sublink_1": ["good", "unusable", "degraded", "unusable"] * 2,
        "status_sublink_2": ["good", "unusable", "good", "good"] * 2,
        # E1 link 3 has too few valid bins; link 4's unusable sub-link 1 has none, its usable sub-link 2 is fine
        "valid_fraction_event_sublink_1": [1.0, 0.0, 0.5, 0.0, 1.0, 0.0, 1.0, 0.0],
        "valid_fraction_event_sublink_2": [1.0, 0.0, 0.5, 0.9, 1.0, 0.0, 1.0, 0.9],
        "split_event": ["train"] * 4 + ["test"] * 4,
        "split_link": ["train", "train", "test", "train"] * 2,
    })
    flags = pd.DataFrame({"event_id": ev, "cml_id": links,
                          "frac_flatline_event_max": [0.0, 0.0, 0.0, 0.0, 0.9, 0.0, 0.0, 0.0],  # E2 link 1 frozen
                          "split_event": event_links.split_event, "split_link": event_links.split_link})
    sig = pd.DataFrame({"event_id": ["E1", "E2"], "significant": [True, True], "sig_reason": ["city_total"] * 2,
                        "rain_class": ["significant", "heavy"], "wet_core_start": pd.to_datetime(["2015-07-01"] * 2),
                        "wet_core_end": pd.to_datetime(["2015-07-02"] * 2), "wet_core_hours": [24.0, 24.0],
                        "split_event": ["train", "test"]})
    tables = {}
    for name, df in (("event_links", event_links), ("event_links_flags_extra", flags), ("event_significance", sig)):
        df.to_parquet(root / f"{name}.parquet")
        tables[name] = {"paths": {"parquet": f"{name}.parquet"}}
    (root / "catalog.json").write_text(json.dumps({"name": "test", "version": "v0", "tables": tables,
                                                   "timeseries": {}}), encoding="utf-8")
    return CMLData.open(root)


def test_event_pairs_filters_and_reasons(rel):
    kept, dropped = rel.event_pairs(return_dropped=True)
    assert len(kept) + len(dropped) == 8
    assert sorted(zip(kept.event_id, kept.cml_id)) == [("E1", 1), ("E1", 4), ("E2", 3), ("E2", 4)]
    why = {(r.event_id, r.cml_id): r.drop_reason for r in dropped.itertuples()}
    assert why == {("E1", 2): "unusable; low_valid", ("E2", 2): "unusable; low_valid",
                   ("E1", 3): "low_valid", ("E2", 1): "flatline"}
    assert {"frac_flatline_event_max", "rain_class", "wet_core_hours"} <= set(kept.columns)

    # disabling every filter keeps all pairs; a split tuple selects (event split, link split)
    assert len(rel.event_pairs(usable_only=False, max_flatline_share=None, min_valid_fraction=None)) == 8
    test_train = rel.event_pairs(split=("test", "train"))
    assert sorted(test_train.cml_id) == [4] and set(test_train.event_id) == {"E2"}
