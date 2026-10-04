"""Step 3 tests: QC flags (L2) on synthetic L1 data with known, injected events.

Two synthetic L1 days (2015-07-28 and -29), 27 links = 54 sub-links, 10 s steps. The background is
quantization jitter on the 1/3 dB RSL grid (stored with 0.1 dB rounding: steps of 0.3 / 0.4 dB) around
-45 dBm, TSL constant at 10 dBm. Injected (column = sub-link):
  A  link 100 s1 : +8 dB 1-sample spike, +8 dB 3-sample spike (flagged); 4- and 10-sample excursions,
                   and a smooth 20 dB rain-like fade (NOT flagged)
  B  link 100 s2 : jitter only; NaN for 1 h on its own (missing_rsl only, no network_gap)
  C  link 200 s1 : 20 min at the floor -90 (rx_floor); a 5 min touch (not flagged); a 12 min floor
                   plateau crossing midnight, 6 min on each side (flagged on both days only thanks to context)
  D  link 200 s2 : rsl exactly constant for 14 h (flatline; the calibrated threshold is 12 h)
  Q  link 800 s1 : rsl exactly constant for 8 h (NOT flatline at 12 h)
  E  link 300 s1 : rsl constant 2 h (not flatline); tsl +1 dB step at one sample (tsl_step)
  A  day 2       : a sharp V-shaped fade (25 dB in 2 min, one sample 8 dB deeper at the bottom): NOT a spike
  G  link 2300 s1: rsl -46 for 1 h, NaN for 11 h, -46 again for 1 h: NOT a flatline (2 h of valid samples)
  F  link 300 s2 : jitter only
  links 400-600 and 11 of the 20 quiet links 800-2700: all their sub-links NaN for 1 h
                   (28 of the 53 alive sub-links, 53 %: network_gap)
  links 800-2700 : 20 quiet links (jitter only), so that ONE missing sub-link is < 5 % of the network
  link 700 s1    : all NaN (dead);  link 700 s2 : rsl but no tsl (has_tl False)
Every assertion checks exactly which samples are flagged, and that nothing else is.
"""

import hashlib
import json

import numpy as np
import pandas as pd
import pytest
import xarray as xr
import yaml

from cmlpipe import qc
from cmlpipe.qc import stage as qc_runner
from cmlpipe.io import schema
from cmlpipe.io import openmrg
from cmlpipe.config import load_config
from conftest import REAL_CONFIG

DAY1, DAY2 = pd.Timestamp("2015-07-28"), pd.Timestamp("2015-07-29")
N = 8640  # samples per day
QUIET = [800 + 100 * k for k in range(20)]
LINKS = [100, 200, 300, 400, 500, 600, 700] + QUIET
GAP_LINKS = [400, 500, 600] + QUIET[:11]  # 28 sub-links NaN together
S1, S2 = 0, 1
RNG = np.random.default_rng(42)


def q(level_db):
    """Put dBm values on the RSL grid: 1/3 dB steps stored rounded to 0.1 dB."""
    return np.round(np.round(np.asarray(level_db) * 3) / 3, 1)


def jitter(n, level=-45.0):
    return q(level + RNG.integers(-1, 2, n) / 3)


def build_days():
    """rsl, tsl arrays (2 days, time, sublink, link) with the injected events; plus expected positions."""
    rsl = np.stack([np.stack([np.stack([jitter(N) for _ in LINKS], axis=1) for _ in (S1, S2)], axis=1)
                    for _ in (DAY1, DAY2)])  # (day, time, sublink, link)
    tsl = np.full(rsl.shape, 10.0)
    L = {c: i for i, c in enumerate(LINKS)}

    # A: spikes and excursions on day 1
    a = rsl[0, :, S1, L[100]]
    a[1000] += 8
    a[2000:2003] += 8
    a[3000:3004] += 8
    a[4000:4010] += 8
    # smooth fade: 20 dB down over 10 min and back over 10 min (rain-like), on the jitter
    tri = np.r_[np.linspace(0, 1, 60), np.linspace(1, 0, 60)]
    a[5000:5120] = q(a[5000:5120] - 20 * tri)

    # A, day 2: sharp V-shaped fade, the bottom sample 8 dB below its neighbours
    v = np.r_[np.linspace(0, 25, 12), [33.0], np.linspace(25, 0, 12)]
    a2 = rsl[1, :, S1, L[100]]
    a2[3000:3025] = q(a2[3000:3025] - v)

    # G: a constant value on both sides of a 13 h gap
    g = rsl[0, :, S1, L[2300]]
    g[0:360] = -46.0
    g[360:4320] = np.nan  # 11 h: the sub-link stays alive (>= 50 % valid)
    g[4320:4680] = -46.0

    # B: one sub-link NaN for 1 h on its own
    rsl[0, 8000:8360, S2, L[100]] = np.nan
    tsl[0, 8000:8360, S2, L[100]] = np.nan

    # C: floor plateau 20 min, a 5 min touch, and a plateau across midnight (6 min each side)
    c1, c2 = rsl[0, :, S1, L[200]], rsl[1, :, S1, L[200]]
    c1[6000:6120] = -90.0
    c1[7000:7030] = -90.0
    c1[N - 36:] = -90.0
    c2[:36] = -90.0

    # D: constant 8 h; E: constant 2 h and a tsl step
    rsl[0, 0:5040, S2, L[200]] = -50.0   # 14 h
    rsl[0, 0:2880, S1, L[800]] = -45.0   # 8 h
    rsl[0, 1000:1720, S1, L[300]] = -47.0
    tsl[0, 4321:, S1, L[300]] = 11.0
    tsl[1, :, S1, L[300]] = 11.0

    # network gap: all sub-links of the GAP_LINKS NaN for 1 h
    for c in GAP_LINKS:
        rsl[0, 7200:7560, :, L[c]] = np.nan
        tsl[0, 7200:7560, :, L[c]] = np.nan

    # link 700: s1 dead, s2 without tsl
    rsl[:, :, S1, L[700]] = np.nan
    tsl[:, :, :, L[700]] = np.nan
    return rsl, tsl


def make_l1(day, rsl, tsl) -> xr.Dataset:
    """A schema-valid L1 dataset for the synthetic signals."""
    n_c = len(LINKS)
    dims2 = ("sublink_id", "cml_id")
    ds = xr.Dataset(
        {"tsl": (schema.DIMS, tsl.astype("float32")), "rsl": (schema.DIMS, rsl.astype("float32")),
         "tl": (schema.DIMS, (tsl - rsl).astype("float32"))},
        coords={
            "time": pd.date_range(day, periods=N, freq="10s").as_unit("ns").to_numpy(),
            "sublink_id": np.array(schema.SUBLINK_IDS), "cml_id": np.array(LINKS, dtype="int64"),
            "raw_sublink_id": (dims2, np.arange(2 * n_c).reshape(2, n_c) + 1),
            "direction": (dims2, np.array([["A"] * n_c, ["B"] * n_c])),
            "tx_site": (dims2, np.array([[0] * n_c, [1] * n_c], dtype="int8")),
            "rx_site": (dims2, np.array([[1] * n_c, [0] * n_c], dtype="int8")),
            "frequency": (dims2, np.full((2, n_c), 28.0)),
            "polarization": (dims2, np.full((2, n_c), "v")),
            "length": ("cml_id", np.full(n_c, 2.0)),
            "site_0_lat": ("cml_id", np.full(n_c, 57.7)), "site_0_lon": ("cml_id", np.full(n_c, 11.9)),
            "site_1_lat": ("cml_id", np.full(n_c, 57.71)), "site_1_lon": ("cml_id", np.full(n_c, 11.91)),
            "site_0_antenna_diameter": ("cml_id", np.full(n_c, 0.3)),
            "site_1_antenna_diameter": ("cml_id", np.full(n_c, 0.3)),
        })
    for name, (units, long_name) in schema.DATA_VARS.items():
        ds[name].attrs = {"units": units, "long_name": long_name}
    for name, (_, units, long_name) in schema.COORDS.items():
        ds[name].attrs = {"long_name": long_name} if units is None else {"units": units, "long_name": long_name}
    ds.attrs = {"time_meaning": "synthetic"}
    return ds


@pytest.fixture(scope="module")
def qc_setup(tmp_path_factory):
    """Config (period = the two synthetic days) with the L1 files written, qc_stats and qc_flags run via run.py."""
    from cmlpipe import cli as run

    tmp = tmp_path_factory.mktemp("qc")
    data = yaml.safe_load(REAL_CONFIG.read_text(encoding="utf-8"))
    data["paths"]["output_root"] = str(tmp / "out")
    data["period"].update({"start": f"{DAY1:%Y-%m-%d}", "end": f"{DAY2:%Y-%m-%d}"})
    cfg_file = tmp / "cfg.yaml"
    cfg_file.write_text(yaml.safe_dump(data), encoding="utf-8")
    cfg = load_config(cfg_file, repo_root=tmp)
    rsl, tsl = build_days()
    for i, day in enumerate((DAY1, DAY2)):
        openmrg.write_l1_day(make_l1(day, rsl[i], tsl[i]), cfg.l1_cml_day_path(day), cfg)
    l1_digest = {d: hashlib.sha256(cfg.l1_cml_day_path(d).read_bytes()).hexdigest() for d in (DAY1, DAY2)}
    base = ["--config", str(cfg_file)]
    assert run.main(base + ["--steps", "qc_stats", "--all"]) == 0
    assert run.main(base + ["--steps", "qc_flags", "--all"]) == 0
    assert run.main(base + ["--steps", "qc_status", "--all"]) == 0
    flags = {}
    for d in (DAY1, DAY2):
        with xr.open_dataset(qc_runner.l2_day_path(cfg, d)) as ds:  # close it: Windows locks open files
            flags[d] = ds.load()
    return cfg, flags, l1_digest, run


def positions(flags: xr.Dataset, name: str) -> set:
    """{(time index, sublink index, link id)} where flag `name` is set."""
    m = schema.decode_flags(flags[schema.FLAGS_VAR].values)[name]
    return {(int(t), int(s), LINKS[int(c)]) for t, s, c in zip(*np.nonzero(m))}


def span(t0, t1, s, link):
    return {(t, s, link) for t in range(t0, t1)}


def test_spike_flags_exactly_the_short_spikes(qc_setup):
    _, flags, _, _ = qc_setup
    assert positions(flags[DAY1], "spike") == span(1000, 1001, S1, 100) | span(2000, 2003, S1, 100)
    assert positions(flags[DAY2], "spike") == set()  # jitter only: nothing


def test_rx_floor_plateau_and_midnight(qc_setup):
    _, flags, _, _ = qc_setup
    assert positions(flags[DAY1], "rx_floor") == span(6000, 6120, S1, 200) | span(N - 36, N, S1, 200)
    assert positions(flags[DAY2], "rx_floor") == span(0, 36, S1, 200)  # needs the previous day as context


def test_flatline_14h_not_8h_or_2h(qc_setup):
    _, flags, _, _ = qc_setup
    assert positions(flags[DAY1], "flatline") == span(0, 5040, S2, 200)
    assert positions(flags[DAY2], "flatline") == set()


def test_tsl_step_single_sample(qc_setup):
    _, flags, _, _ = qc_setup
    assert positions(flags[DAY1], "tsl_step") == {(4321, S1, 300)}
    assert positions(flags[DAY2], "tsl_step") == set()  # tsl stays at 11 across midnight


def test_network_gap_only_for_the_shared_gap(qc_setup):
    _, flags, _, _ = qc_setup
    expected = set()
    for c in GAP_LINKS:
        expected |= span(7200, 7560, S1, c) | span(7200, 7560, S2, c)
    assert positions(flags[DAY1], "network_gap") == expected
    missing = positions(flags[DAY1], "missing_rsl")
    assert expected <= missing and span(8000, 8360, S2, 100) <= missing  # B: missing only
    assert not (span(8000, 8360, S2, 100) & positions(flags[DAY1], "network_gap"))
    frac = flags[DAY1].network_nan_fraction.values
    assert frac[7300] == pytest.approx(28 / 53)  # 53 alive sub-links (link 700 s1 is dead)
    assert frac[8100] == pytest.approx(1 / 53)   # sub-link B alone: 1.9 %, below the 5 % threshold


def test_nothing_else_is_flagged(qc_setup):
    _, flags, _, _ = qc_setup
    assert positions(flags[DAY1], "out_of_range") == set() and positions(flags[DAY2], "out_of_range") == set()
    # missing flags only where the data are NaN
    for day, f in flags.items():
        miss = positions(f, "missing_rsl")
        expected = span(0, N, S1, 700)
        if day == DAY1:
            expected |= span(8000, 8360, S2, 100) | span(360, 4320, S1, 2300)
            for c in GAP_LINKS:
                expected |= span(7200, 7560, S1, c) | span(7200, 7560, S2, c)
        assert miss == expected


def test_status_table(qc_setup):
    cfg, _, _, _ = qc_setup
    st = pd.read_csv(cfg.tables_dir / "sublink_qc_status.csv").set_index(["cml_id", "sublink_id"])
    assert not st.loc[(700, "sublink_1"), "has_rsl"] and st.loc[(700, "sublink_1"), "reason"].startswith("no_rsl")
    assert st.loc[(700, "sublink_2"), "has_rsl"] and not st.loc[(700, "sublink_2"), "has_tl"]
    assert st.loc[(700, "sublink_2"), "status"] == "unusable" and "no_tl" in st.loc[(700, "sublink_2"), "reason"]
    assert st.loc[(300, "sublink_2"), "status"] == "good"


def test_summary_and_gap_tables(qc_setup):
    cfg, _, _, _ = qc_setup
    summ = pd.read_csv(cfg.tables_dir / "qc_daily_summary.csv")
    assert len(summ) == 2 * 2 * len(LINKS)
    row = summ[(summ.day == "2015-07-28") & (summ.cml_id == 100) & (summ.sublink_id == "sublink_1")].iloc[0]
    assert row.n_spike == 4 and row.n_valid_rsl == N
    assert summ[summ.cml_id == 700].pair_corr.isna().all()  # one sub-link dead -> NaN
    gaps = pd.read_csv(cfg.tables_dir / "network_gaps.csv")
    assert len(gaps) == 1 and gaps.duration_min.iloc[0] == 60.0
    assert pd.Timestamp(gaps.start.iloc[0]) == DAY1 + pd.Timedelta(seconds=7200 * 10)


def test_l1_unchanged_and_rerun_identical(qc_setup):
    cfg, flags, l1_digest, run = qc_setup
    for d, digest in l1_digest.items():
        assert hashlib.sha256(cfg.l1_cml_day_path(d).read_bytes()).hexdigest() == digest
    base = ["--config", str(cfg.config_file), "--steps", "qc_flags", "--all"]
    assert run.main(base) == 0  # skipped: same stage and input hashes
    last = json.loads(cfg.manifest_path.read_text(encoding="utf-8"))[-1]
    assert all(r["status"] == "skipped_same_config" for r in last["days"])
    assert run.main(base + ["--force"]) == 0
    with xr.open_dataset(qc_runner.l2_day_path(cfg, DAY1)) as again:
        assert np.array_equal(again[schema.FLAGS_VAR].values, flags[DAY1][schema.FLAGS_VAR].values)


def test_validate_flags_l2(qc_setup):
    cfg, flags, _, _ = qc_setup
    with xr.open_dataset(cfg.l1_cml_day_path(DAY1)) as l1:
        l1 = l1.load()
    schema.validate_flags_l2(flags[DAY1], l1)
    with pytest.raises(schema.CMLValidationError, match="dims"):
        schema.validate_flags_l2(flags[DAY1].transpose("time", "cml_id", "sublink_id"), l1)
    bad = flags[DAY1].copy(deep=True)
    bad[schema.FLAGS_VAR].values[0, 0, 0] |= np.uint16(1 << 9)
    with pytest.raises(schema.CMLValidationError, match="not defined"):
        schema.validate_flags_l2(bad, l1)


def test_decode_roundtrip_and_blocking():
    rng = np.random.default_rng(0)
    masks = {name: rng.random((50, 3)) < 0.3 for name in schema.FLAGS}
    packed = qc.combine(masks)
    assert packed.dtype == np.uint16
    back = schema.decode_flags(packed)
    assert all(np.array_equal(back[k], masks[k]) for k in masks)
    blk = qc.blocking_mask(packed, ["spike", "rx_floor"])
    assert np.array_equal(blk, masks["spike"] | masks["rx_floor"])


def test_out_of_range_unit():
    rsl = np.array([[-45.0], [-15.0], [np.nan], [-101.0]])
    tsl = np.array([[10.0], [10.0], [40.0], [10.0]])
    assert qc.flag_out_of_range(rsl, tsl, [-100, -20], [-15, 30]).ravel().tolist() == [False, True, True, True]


def test_run_lengths_unit():
    m = np.array([[1, 1, 0, 1, 1, 1, 0]], dtype=bool).T
    assert qc.run_lengths(m).ravel().tolist() == [2, 2, 0, 3, 3, 3, 0]
