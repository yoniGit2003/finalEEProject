"""Geometry stage: synthetic geometry only (no real data needed).

Synthetic radar grid used below: 5 x 5 cells of 2 km, centres x = 1000, 3000, ..., 9000 and
y = 9000, 7000, ..., 1000 (y DEScending, like the OpenMRG radar), so edges are at 0, 2000, ..., 10000 m.
"""

import filecmp

import numpy as np
import pandas as pd
import pyproj
import pytest
import shapely
import xarray as xr

from cmlpipe import geometry
from cmlpipe.geometry import stage as geometry_runner
from cmlpipe.geometry import RadarGrid, apply_weights, merge_sites, path_grid_weights, point_to_path

PROJ = "+proj=stere +lat_ts=60 +ellps=bessel +lon_0=14 +lat_0=90"


@pytest.fixture
def grid():
    return RadarGrid(np.arange(1000.0, 10000.0, 2000.0), np.arange(9000.0, 0.0, -2000.0), PROJ)


def weights_dict(iy, ix, w):
    return {(int(a), int(b)): float(c) for a, b, c in zip(iy, ix, w)}


# ---------------------------------------------------------------------------
# point_to_path
# ---------------------------------------------------------------------------
def test_point_perpendicular_to_midpoint():
    d, pos = point_to_path(500.0, 300.0, 0.0, 0.0, 1000.0, 0.0)
    assert d == pytest.approx(300.0) and pos == pytest.approx(0.5)


def test_point_beyond_site_1():
    d, pos = point_to_path(1300.0, 400.0, 0.0, 0.0, 1000.0, 0.0)
    assert d == pytest.approx(500.0) and pos == pytest.approx(1.0)  # distance to the end (3-4-5)


def test_point_beyond_site_0():
    d, pos = point_to_path(-200.0, 0.0, 0.0, 0.0, 1000.0, 0.0)
    assert d == pytest.approx(200.0) and pos == pytest.approx(0.0)


def test_point_on_path():
    d, pos = point_to_path(250.0, 250.0, 0.0, 0.0, 1000.0, 1000.0)
    assert d == pytest.approx(0.0, abs=1e-9) and pos == pytest.approx(0.25)


# ---------------------------------------------------------------------------
# path_grid_weights
# ---------------------------------------------------------------------------
def test_horizontal_path_three_cells(grid):
    # y = 5000 is a cell centre (row 2); x 1000 -> 5000 crosses edges 2000 and 4000: pieces 1, 2, 1 km
    iy, ix, w = path_grid_weights(1000.0, 5000.0, 5000.0, 5000.0, grid)
    assert weights_dict(iy, ix, w) == pytest.approx({(2, 0): 0.25, (2, 1): 0.5, (2, 2): 0.25})


def test_diagonal_through_corners_no_double_count(grid):
    # passes exactly through the corners (2000, 2000) and (4000, 4000)
    iy, ix, w = path_grid_weights(1000.0, 1000.0, 5000.0, 5000.0, grid)
    assert w.sum() == pytest.approx(1.0, abs=1e-12)
    assert len(w) == 3  # only the 3 cells on the diagonal, none of the corner-touching neighbours
    assert weights_dict(iy, ix, w) == pytest.approx({(4, 0): 0.25, (3, 1): 0.5, (2, 2): 0.25})


def test_path_along_a_cell_edge(grid):
    # runs exactly on the edge x = 2000: each piece must go to ONE side only
    iy, ix, w = path_grid_weights(2000.0, 500.0, 2000.0, 3500.0, grid)
    assert w.sum() == pytest.approx(1.0)
    assert len(set(ix)) == 1


def test_path_starting_outside(grid):
    # 4 km long, the first 1 km is west of the grid (x < 0)
    iy, ix, w = path_grid_weights(-1000.0, 3000.0, 3000.0, 3000.0, grid)
    assert w.sum() == pytest.approx(0.75)
    assert weights_dict(iy, ix, w) == pytest.approx({(3, 0): 0.5, (3, 1): 0.25})


def test_path_inside_one_cell(grid):
    iy, ix, w = path_grid_weights(1200.0, 1200.0, 1800.0, 1800.0, grid)
    assert list(w) == [1.0] and (iy[0], ix[0]) == (4, 0)


def test_weights_match_polygon_intersection(grid):
    """For generic paths (not on edges) the result equals shapely's polygon-intersection lengths."""
    rng = np.random.default_rng(1)
    xe, ye = grid.x_edges(), grid.y_edges()
    for _ in range(50):
        x0, y0, x1, y1 = rng.uniform(-1500, 11500, 4)
        iy, ix, w = path_grid_weights(x0, y0, x1, y1, grid)
        line = shapely.LineString([(x0, y0), (x1, y1)])
        ref = {}
        for j in range(5):
            for i in range(5):
                cell = shapely.box(xe[i], ye[j + 1], xe[i + 1], ye[j])
                length = line.intersection(cell).length
                if length > 0:
                    ref[(j, i)] = length / line.length
        assert weights_dict(iy, ix, w) == pytest.approx(ref, abs=1e-9)


def test_weights_match_poligrain(grid):
    """Independent implementation (poligrain, polygon intersection) on the same descending-y grid."""
    spatial = pytest.importorskip("poligrain.spatial")
    xg, yg = np.meshgrid(grid.x, grid.y)
    rng = np.random.default_rng(2)
    for _ in range(10):
        x0, y0, x1, y1 = rng.uniform(500, 9500, 4)
        iy, ix, w = path_grid_weights(x0, y0, x1, y1, grid)
        ours = np.zeros(grid.shape)
        ours[iy, ix] = w
        theirs = spatial.calc_intersect_weights(x0, y0, x1, y1, xg, yg)
        assert np.abs(ours - theirs).max() < 1e-9


# ---------------------------------------------------------------------------
# apply_weights
# ---------------------------------------------------------------------------
def field(grid, values_yx, n_t=3):
    data = np.broadcast_to(values_yx, (n_t,) + grid.shape).astype("float32").copy()
    return xr.DataArray(data, dims=("time", "y", "x"),
                        coords={"time": pd.date_range("2015-07-28", periods=n_t, freq="5min"),
                                "y": grid.y, "x": grid.x})


def weights_df(grid, cml_id, *path):
    iy, ix, w = path_grid_weights(*path, grid)
    return pd.DataFrame({"cml_id": cml_id, "iy": iy, "ix": ix, "weight": w})


def test_apply_uniform_field(grid):
    w = weights_df(grid, 1, 1000.0, 5000.0, 5000.0, 5000.0)
    out = apply_weights(field(grid, np.full(grid.shape, 7.0)), w, threshold=1.0, time_chunk=2, grid=grid)
    assert np.allclose(out["mean"], 7.0) and np.allclose(out["std"], 0.0)
    assert np.allclose(out["wet_fraction"], 1.0) and np.allclose(out["valid_weight_fraction"], 1.0)


def test_apply_two_cells_50_50(grid):
    values = np.zeros(grid.shape)
    values[2, 1] = 10.0
    w = weights_df(grid, 1, 1000.0, 5000.0, 3000.0, 5000.0)  # 1 km in cell (2,0), 1 km in (2,1)
    out = apply_weights(field(grid, values), w, threshold=5.0)
    assert np.allclose(out["mean"], 5.0) and np.allclose(out["std"], 5.0)
    assert np.allclose(out["wet_fraction"], 0.5)


def test_apply_nan_cell_quarter(grid):
    values = np.full(grid.shape, 4.0)
    values[2, 2] = 8.0
    values[2, 0] = np.nan  # under the first quarter of the path
    w = weights_df(grid, 1, 1000.0, 5000.0, 5000.0, 5000.0)  # 0.25 / 0.5 / 0.25
    out = apply_weights(field(grid, values), w)
    assert np.allclose(out["valid_weight_fraction"], 0.75)
    assert np.allclose(out["mean"], (0.5 * 4 + 0.25 * 8) / 0.75)  # over the valid cells only
    assert np.allclose(out["max"], 8.0)


def test_apply_single_hot_cell_and_two_links(grid):
    values = np.ones(grid.shape)
    values[3, 3] = 50.0
    w = pd.concat([weights_df(grid, 2, 1000.0, 3000.0, 9000.0, 3000.0),   # row 3, crosses (3, 3)
                   weights_df(grid, 1, 1000.0, 9000.0, 9000.0, 9000.0)])  # row 0, does not
    out = apply_weights(field(grid, values), w)
    assert list(out.cml_id.values) == [1, 2]
    assert np.allclose(out["max"].sel(cml_id=2), 50.0) and np.allclose(out["max"].sel(cml_id=1), 1.0)


def test_apply_all_nan_step(grid):
    data = field(grid, np.ones(grid.shape))
    data[1] = np.nan
    out = apply_weights(data, weights_df(grid, 1, 1000.0, 5000.0, 5000.0, 5000.0))
    assert np.isnan(out["mean"][1]).all() and float(out["valid_weight_fraction"][1, 0]) == 0.0


# ---------------------------------------------------------------------------
# sites, units
# ---------------------------------------------------------------------------
def test_sites_5m_merge_20m_do_not():
    assert len(set(merge_sites([0.0, 5.0], [0.0, 0.0], 10.0))) == 1
    assert len(set(merge_sites([0.0, 20.0], [0.0, 0.0], 10.0))) == 2


def synthetic_inputs():
    """A small network near Gothenburg, metadata lengths = true geodesic lengths, and a radar grid around it."""
    geod = pyproj.Geod(ellps="WGS84")
    ends = [  # (lon0, lat0, lon1, lat1)
        (11.90, 57.70, 11.95, 57.71),
        (11.95, 57.71, 11.97, 57.68),          # shares a site with the first link
        (11.9501, 57.71, 12.00, 57.72),        # ~6 m east of that site: merged
        (12.10, 57.60, 12.20, 57.65),
    ]
    rows = []
    for k, (a, b, c, d) in enumerate(ends):
        rows.append({"cml_id": 100 + k, "site_0_lon": a, "site_0_lat": b, "site_1_lon": c, "site_1_lat": d,
                     "length_km_metadata": geod.line_length([a, c], [b, d]) / 1000.0,
                     "frequency_ghz_sublink_1": 28.0, "frequency_ghz_sublink_2": 29.0, "polarization": "v",
                     "site_0_antenna_diameter_m": 0.3, "site_1_antenna_diameter_m": 0.6})
    links = pd.DataFrame(rows)
    gauges = pd.DataFrame({"gauge_id": ["G1", "G2"], "network": ["city", "smhi"], "gauge_type": ["weighing"] * 2,
                           "lon": [11.93, 12.15], "lat": [57.705, 57.62], "quantization_mm": [0.1, 0.1]})
    rx, ry = geometry.project_to_radar([11.8, 12.3], [57.55, 57.8], PROJ)
    gx = np.arange(np.floor(rx.min() / 2000) * 2000 + 1000, rx.max() + 2000, 2000.0)
    gy = np.arange(np.ceil(ry.max() / 2000) * 2000 - 1000, ry.min() - 2000, -2000.0)
    grid = RadarGrid(gx, gy, PROJ)
    crs = pyproj.CRS.from_proj4(PROJ)
    xx, yy = np.meshgrid(gx, gy)
    lon, lat = pyproj.Transformer.from_crs(crs, crs.geodetic_crs, always_xy=True).transform(xx, yy)
    x3006, y3006 = geometry.project_points(lon, lat, "EPSG:3006")
    return links, gauges, grid, (x3006, y3006), lon, lat


GCFG = {"crs_metric": "EPSG:3006", "site_merge_tol_m": 10, "gauge_match_max_m": 5000,
        "gauge_flag_radii_m": [500, 1000, 2000], "neighbour_max_m": 5000,
        "length_mismatch_tol": {"relative": 0.02, "absolute_m": 30}}


def test_length_in_metres_matches_metadata_km():
    links, gauges, grid, cells, _, _ = synthetic_inputs()
    lg = geometry.build_geometry(links, gauges, grid, cells, GCFG)["link_geometry"]
    ratio = lg.length_m / (1000.0 * lg.length_km_metadata)
    assert np.all(np.abs(ratio - 1.0) < 0.02)  # a km/m mix-up would be off by 1000x
    assert not lg.length_mismatch.any()


def test_build_geometry_tables():
    links, gauges, grid, cells, _, _ = synthetic_inputs()
    t = geometry.build_geometry(links, gauges, grid, cells, GCFG)
    lg = t["link_geometry"].set_index("cml_id")
    assert len(t["sites"]) == 6  # 8 ends; 101 and 102 share one end with 100, which 102 shares within ~6 m
    assert lg.loc[100, "site_1_id"] == lg.loc[101, "site_0_id"] == lg.loc[102, "site_0_id"]
    assert np.allclose(lg.radar_coverage, 1.0)
    sums = t["path_weights"].groupby("cml_id").weight.sum()
    assert np.allclose(sums.reindex(lg.index), lg.radar_coverage)
    nb = t["link_neighbours"].set_index(["cml_id_a", "cml_id_b"])
    assert nb.loc[(100, 101), "shares_site"] and nb.loc[(100, 101), "min_path_distance_m"] == pytest.approx(0, abs=1e-6)
    assert (t["link_gauge_matches"].distance_to_path_m <= 5000).all()
    assert t["gauge_radar_pixel"].inside_radar_grid.all()


# ---------------------------------------------------------------------------
# stage: end to end on the synthetic config (raw metadata from conftest), idempotence and skip
# ---------------------------------------------------------------------------
@pytest.fixture
def geometry_cfg(synthetic_cfg):
    """synthetic_cfg + a gauges.csv and a one-file radar L1 grid that covers the conftest links."""
    cfg = synthetic_cfg
    _, gauges, _, _, _, _ = synthetic_inputs()
    cfg.tables_dir.mkdir(parents=True, exist_ok=True)
    gauges.rename(columns={"gauge_id": "id"}).to_csv(cfg.tables_dir / "gauges.csv", index=False)
    rx, ry = geometry.project_to_radar([11.3, 12.5], [56.8, 57.9], PROJ)  # conftest links: 11.5-12.3 E, 57.0-57.8 N
    gx = np.arange(np.floor(rx.min() / 2000) * 2000 + 1000, rx.max() + 2000, 2000.0)
    gy = np.arange(np.ceil(ry.max() / 2000) * 2000 - 1000, ry.min() - 2000, -2000.0)
    crs = pyproj.CRS.from_proj4(PROJ)
    xx, yy = np.meshgrid(gx, gy)
    lon, lat = pyproj.Transformer.from_crs(crs, crs.geodetic_crs, always_xy=True).transform(xx, yy)
    x3006, y3006 = geometry.project_points(lon, lat, "EPSG:3006")
    ds = xr.Dataset(coords={"y": gy, "x": gx, "lon": (("y", "x"), lon), "lat": (("y", "x"), lat),
                            "x_3006": (("y", "x"), x3006), "y_3006": (("y", "x"), y3006)},
                    attrs={"crs_proj4": PROJ})
    cfg.l1_radar_dir.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(cfg.l1_radar_dir / "radar_5min_2015-07.nc")
    return cfg


class _Log:
    def info(self, *a, **k):
        pass


def test_stage_idempotent_and_skips(geometry_cfg, tmp_path):
    cfg = geometry_cfg
    rec = []
    geometry_runner.run_geometry(cfg, False, _Log(), rec)
    paths = geometry_runner.output_paths(cfg)
    first = tmp_path / "first"
    first.mkdir()
    for name, p in paths.items():
        (first / p.name).write_bytes(p.read_bytes())
    weights_1 = xr.load_dataset(paths["path_weights"])

    geometry_runner.run_geometry(cfg, False, _Log(), rec)  # up to date -> skipped
    assert rec[-1]["status"] == "skipped_same_config"

    geometry_runner.run_geometry(cfg, True, _Log(), rec)   # forced rebuild -> identical content
    assert rec[-1]["status"] == "written"
    for name in geometry_runner.TABLES:
        assert filecmp.cmp(first / f"{name}.csv", paths[name], shallow=False), name
    weights_2 = xr.load_dataset(paths["path_weights"])
    for v in ("cml_id", "iy", "ix", "weight"):
        assert np.array_equal(weights_1[v].values, weights_2[v].values)

    lg = pd.read_csv(paths["link_geometry"])
    assert len(lg) == 4  # conftest metadata: links 100, 200, 300, 400
    w, grid = geometry_runner.read_weights(paths["path_weights"])
    assert np.allclose(w.groupby("cml_id").weight.sum().reindex(lg.cml_id).fillna(0), lg.radar_coverage)


def test_stage_rejects_other_proj_string(geometry_cfg):
    geometry_cfg.data["geometry"]["radar_proj_string"] = "+proj=stere +lat_ts=60 +ellps=WGS84 +lon_0=14 +lat_0=90"
    with pytest.raises(ValueError, match="proj string"):
        geometry_runner.run_geometry(geometry_cfg, True, _Log(), [])
