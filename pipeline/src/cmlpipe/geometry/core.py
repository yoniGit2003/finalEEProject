"""Geometry of the CML network: where links, gauges and radar pixels are, and how they relate.

Pure functions only (no file I/O, no config object); the stage that reads inputs and writes tables
is cmlpipe.geometry.stage. No time series are used here, except by apply_weights, which the
NEXT pipeline step calls to average radar fields along link paths.

Two coordinate systems, on purpose:
  * EPSG:3006 (SWEREF 99 TM, metres) for every LENGTH and DISTANCE (link length, gauge-to-path
    distance, neighbour distance, hull area). It is the national metric CRS of Sweden, so its scale
    error over Gothenburg is tiny (< 1e-4).
  * The radar's own polar stereographic projection for the path/cell INTERSECTION. Only there are the
    radar cells exact 2 km squares. The stereographic scale over Gothenburg differs from 1 (about
    1.4 %), so we keep only FRACTIONS of the path length from it. A fraction is a ratio of two lengths
    measured in the same projection, so a (locally) uniform scale factor cancels.
All lengths and distances carry their unit in the name (..._m). Nothing here works in km.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
import pyproj
import shapely
import xarray as xr

WGS84 = "EPSG:4326"

# Two path parameters t closer than this are the same crossing. Needed when a path passes exactly
# through a cell corner: the x-edge and the y-edge crossing then differ only by rounding (~1e-16), and
# without merging they would create a zero-length piece that could be assigned to a third cell.
T_MERGE_TOL = 1e-12


# ---------------------------------------------------------------------------
# Projections
# ---------------------------------------------------------------------------
def project_points(lon, lat, crs: str, src_crs: str = WGS84) -> tuple[np.ndarray, np.ndarray]:
    """lon/lat in degrees (src_crs) -> x, y in metres of `crs`. always_xy: input order is lon, lat."""
    tr = pyproj.Transformer.from_crs(src_crs, crs, always_xy=True)
    x, y = tr.transform(np.asarray(lon, dtype="float64"), np.asarray(lat, dtype="float64"))
    return np.asarray(x, dtype="float64"), np.asarray(y, dtype="float64")


def project_to_radar(lon, lat, proj_string: str) -> tuple[np.ndarray, np.ndarray]:
    """WGS84 lon/lat -> x, y of the radar projection, with NO datum shift.

    The radar proj string names only the Bessel ellipsoid, no datum. The radar file's lon/lat are the
    inverse projection of its x/y on that ellipsoid (checked to 0.2 m in the inspection note). So
    "lon/lat as written in the radar file" and "WGS84 lon/lat of the links" are treated as the same
    numbers: we project the link lon/lat with the radar's own geographic CRS as the source, which applies
    the map projection only. Whether the composite really sits on another datum is undocumented; the
    acceptance report measures what a shift of that size would change.
    """
    crs = pyproj.CRS.from_proj4(proj_string)
    return project_points(lon, lat, crs.to_wkt(), src_crs=crs.geodetic_crs.to_wkt())


def azimuth_deg(x0, y0, x1, y1) -> np.ndarray:
    """Grid azimuth from (x0, y0) to (x1, y1): 0 = north (+y), 90 = east (+x), clockwise, in [0, 360)."""
    az = np.degrees(np.arctan2(np.asarray(x1) - x0, np.asarray(y1) - y0))
    return np.mod(az, 360.0)


# ---------------------------------------------------------------------------
# Radar grid
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class RadarGrid:
    """A regular grid of square cells given by its 1-D cell CENTRES (in the radar projection, metres).

    The step may be negative (OpenMRG: x step +2000, y step -2000, i.e. y runs north to south).
    Cell i spans centre[i] - step/2 .. centre[i] + step/2.
    """

    x: np.ndarray
    y: np.ndarray
    proj_string: str = ""

    def __post_init__(self):
        for name in ("x", "y"):
            c = np.asarray(getattr(self, name), dtype="float64")
            if c.ndim != 1 or c.size < 2:
                raise ValueError(f"grid {name} must be 1-D with at least 2 cell centres")
            d = np.diff(c)
            if not np.allclose(d, d[0], rtol=0, atol=1e-6 * abs(d[0])):
                raise ValueError(f"grid {name} is not regular (steps {d.min()} .. {d.max()})")
            object.__setattr__(self, name, c)

    @property
    def dx(self) -> float:
        """Cell step along x (m, may be negative)."""
        return float(self.x[1] - self.x[0])

    @property
    def dy(self) -> float:
        """Cell step along y (m, negative for OpenMRG)."""
        return float(self.y[1] - self.y[0])

    @property
    def shape(self) -> tuple[int, int]:
        """(ny, nx), the order of the radar arrays (time, y, x)."""
        return self.y.size, self.x.size

    def x_edges(self) -> np.ndarray:
        """nx + 1 cell edges along x, in the same order as the centres."""
        return self.x[0] - self.dx / 2 + self.dx * np.arange(self.x.size + 1)

    def y_edges(self) -> np.ndarray:
        """ny + 1 cell edges along y, in the same order as the centres."""
        return self.y[0] - self.dy / 2 + self.dy * np.arange(self.y.size + 1)

    def cell_index(self, px, py) -> tuple[np.ndarray, np.ndarray]:
        """Indices (iy, ix) of the cell that contains each point; -1 where the point is outside the grid.

        Index = floor((p - first edge) / step). A point exactly on a shared edge gets ONE cell (the floor
        decides), so nothing is ever counted twice.
        """
        fx = (np.asarray(px, dtype="float64") - (self.x[0] - self.dx / 2)) / self.dx
        fy = (np.asarray(py, dtype="float64") - (self.y[0] - self.dy / 2)) / self.dy
        ix = np.floor(fx).astype("int64")
        iy = np.floor(fy).astype("int64")
        inside = (ix >= 0) & (ix < self.x.size) & (iy >= 0) & (iy < self.y.size)
        return np.where(inside, iy, -1), np.where(inside, ix, -1)

    def outline(self) -> shapely.Polygon:
        """The outer boundary of the grid as a polygon (radar projection)."""
        xe, ye = self.x_edges(), self.y_edges()
        return shapely.box(min(xe[0], xe[-1]), min(ye[0], ye[-1]), max(xe[0], xe[-1]), max(ye[0], ye[-1]))


def path_grid_weights(x0: float, y0: float, x1: float, y1: float, grid: RadarGrid
                      ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Cells crossed by the straight path (x0, y0) -> (x1, y1) and the FRACTION of the path in each.

    Method (exact for straight lines on a rectangular grid): parametrize the path as
    p(t) = p0 + t (p1 - p0), t in [0, 1]. Collect every t where the path crosses a vertical or horizontal
    cell edge. Between two consecutive crossings the path lies in a single cell, found from the piece's
    MIDPOINT; the piece's weight is its share of t (= share of the path length). Why midpoints and not
    polygon intersections: a path running exactly along an edge, or through a corner, touches two or four
    cells with zero area; polygon intersection would give the same piece to both neighbours (double count),
    while the midpoint falls in exactly one cell.

    Pieces outside the grid are dropped, so sum(weights) = radar_coverage (1 if the path is fully inside).
    A zero-length path gets weight 1 in the cell that contains it (if any).

    Returns iy, ix, weight (1-D arrays, in order along the path from site_0 to site_1).
    """
    dxp, dyp = x1 - x0, y1 - y0
    if dxp == 0 and dyp == 0:
        iy, ix = grid.cell_index([x0], [y0])
        keep = iy >= 0
        return iy[keep], ix[keep], np.ones(int(keep.sum()))

    ts = [np.array([0.0, 1.0])]
    if dxp != 0:
        ts.append((grid.x_edges() - x0) / dxp)
    if dyp != 0:
        ts.append((grid.y_edges() - y0) / dyp)
    t = np.concatenate(ts)
    t = np.sort(t[(t >= 0.0) & (t <= 1.0)])
    t = t[np.concatenate([[True], np.diff(t) > T_MERGE_TOL])]  # merge corner duplicates
    t[-1] = 1.0  # the merge may have kept a value a hair below 1 instead of 1 itself

    dt = np.diff(t)
    tm = 0.5 * (t[:-1] + t[1:])
    iy, ix = grid.cell_index(x0 + tm * dxp, y0 + tm * dyp)
    keep = (iy >= 0) & (dt > 0)
    iy, ix, dt = iy[keep], ix[keep], dt[keep]
    if iy.size == 0:
        return iy, ix, dt
    # A straight line enters a convex cell at most once, so consecutive pieces in the same cell can only
    # appear through rounding; merge them anyway so each (iy, ix) appears once.
    new = np.concatenate([[True], (np.diff(iy) != 0) | (np.diff(ix) != 0)])
    groups = np.cumsum(new) - 1
    w = np.bincount(groups, weights=dt)
    return iy[new], ix[new], w


def point_to_path(px, py, x0, y0, x1, y1) -> tuple[np.ndarray, np.ndarray]:
    """Distance from point(s) P to the path SEGMENT site_0 -> site_1, and where along it the closest point is.

    position = projection of P onto the line, as a fraction of the path (0 = site_0, 1 = site_1), clipped
    to [0, 1]; the distance is to that clipped point, i.e. to the nearest END when P lies beyond it.
    All arguments broadcast, so one point vs many links (or the reverse) works. Units: those of the input.
    """
    px, py, x0, y0, x1, y1 = (np.asarray(v, dtype="float64") for v in (px, py, x0, y0, x1, y1))
    vx, vy = x1 - x0, y1 - y0
    vv = vx * vx + vy * vy
    with np.errstate(invalid="ignore", divide="ignore"):
        t = ((px - x0) * vx + (py - y0) * vy) / vv
    t = np.where(vv > 0, np.clip(t, 0.0, 1.0), 0.0)  # zero-length path: the point itself
    cx, cy = x0 + t * vx, y0 + t * vy
    return np.hypot(px - cx, py - cy), t


# ---------------------------------------------------------------------------
# Sites
# ---------------------------------------------------------------------------
def merge_sites(x, y, tol_m: float) -> np.ndarray:
    """Group points that lie within tol_m of each other into one site; return a 0-based label per point.

    Single linkage: A and C end up together if A-B and B-C are both within tol_m (the acceptance report
    gives the largest member spread, so a chain longer than tol_m would be visible). Labels are numbered in
    order of first appearance, so the result does not depend on anything but the input order.
    """
    x, y = np.asarray(x, dtype="float64"), np.asarray(y, dtype="float64")
    n = x.size
    close = np.hypot(x[:, None] - x[None, :], y[:, None] - y[None, :]) <= tol_m
    parent = np.arange(n)

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in zip(*np.nonzero(np.triu(close, k=1))):
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)
    roots = np.array([find(i) for i in range(n)])
    _, first = np.unique(roots, return_index=True)
    order = {r: k for k, r in enumerate(roots[np.sort(first)])}
    return np.array([order[r] for r in roots], dtype="int64")


# ---------------------------------------------------------------------------
# Network geometry (pure: DataFrames in, DataFrames out)
# ---------------------------------------------------------------------------
LINK_INPUT_COLUMNS = ["cml_id", "site_0_lon", "site_0_lat", "site_1_lon", "site_1_lat", "length_km_metadata",
                      "frequency_ghz_sublink_1", "frequency_ghz_sublink_2", "polarization",
                      "site_0_antenna_diameter_m", "site_1_antenna_diameter_m"]
GAUGE_INPUT_COLUMNS = ["gauge_id", "network", "gauge_type", "lon", "lat", "quantization_mm"]


def build_geometry(links: pd.DataFrame, gauges: pd.DataFrame, grid: RadarGrid, grid_cells_3006: tuple,
                   gcfg: dict) -> dict:
    """All geometry tables from link and gauge positions and the radar grid.

    links  : one row per link, columns LINK_INPUT_COLUMNS (lon/lat WGS84 degrees).
    gauges : one row per gauge, columns GAUGE_INPUT_COLUMNS.
    grid   : the radar grid (cell centres in the radar projection).
    grid_cells_3006 : (x_3006, y_3006) 2-D arrays (ny, nx) of the radar cell centres in EPSG:3006, used only
             for the gauge-to-cell-centre distance (a distance, so it is measured in EPSG:3006).
    gcfg   : the config's `geometry` section.
    Returns a dict of DataFrames: sites, link_geometry, path_weights, link_gauge_matches,
    gauge_radar_pixel, link_neighbours.
    """
    crs = gcfg["crs_metric"]
    links = links.sort_values("cml_id").reset_index(drop=True)
    n = len(links)

    # --- metric coordinates of both ends ---
    x0, y0 = project_points(links.site_0_lon, links.site_0_lat, crs)
    x1, y1 = project_points(links.site_1_lon, links.site_1_lat, crs)

    # --- sites: both ends of every link, merged within site_merge_tol_m ---
    ex, ey = np.concatenate([x0, x1]), np.concatenate([y0, y1])
    elon = np.concatenate([links.site_0_lon, links.site_1_lon])
    elat = np.concatenate([links.site_0_lat, links.site_1_lat])
    # Order the ends as (link 0 site_0, link 0 site_1, link 1 site_0, ...) so site numbers follow cml_id.
    order = np.column_stack([np.arange(n), np.arange(n) + n]).ravel()
    labels_ordered = merge_sites(ex[order], ey[order], gcfg["site_merge_tol_m"])
    labels = np.empty(2 * n, dtype="int64")
    labels[order] = labels_ordered
    site_id = labels + 1  # 1-based ids are easier to read in tables
    s0, s1 = site_id[:n], site_id[n:]
    ends = pd.DataFrame({"site_id": site_id, "lon": elon, "lat": elat, "x": ex, "y": ey})
    ends["cml_id"] = np.concatenate([links.cml_id, links.cml_id])
    g = ends.groupby("site_id")
    sites = pd.DataFrame({
        "lon": g.lon.mean(), "lat": g.lat.mean(), "x_3006_m": g.x.mean(), "y_3006_m": g.y.mean(),
        "n_links": g.cml_id.nunique(),  # distinct links with an end at this site
        "n_ends": g.size(),             # differs from n_links only if both ends of a link merged here
        # largest distance of a merged end from the site centre: shows how much merging moved anything
        "max_member_offset_m": g.apply(lambda d: float(np.hypot(d.x - d.x.mean(), d.y - d.y.mean()).max()),
                                       include_groups=False),
    }).reset_index()

    # --- link geometry ---
    length_m = np.hypot(x1 - x0, y1 - y0)
    length_km_meta = links.length_km_metadata.to_numpy(float)
    length_diff_m = length_m - 1000.0 * length_km_meta
    with np.errstate(invalid="ignore", divide="ignore"):
        length_rel_diff = length_diff_m / (1000.0 * length_km_meta)
    tol = gcfg["length_mismatch_tol"]
    mismatch = (np.abs(length_rel_diff) > tol["relative"]) & (np.abs(length_diff_m) > tol["absolute_m"])

    # --- radar path weights (in the radar projection) ---
    rx0, ry0 = project_to_radar(links.site_0_lon, links.site_0_lat, grid.proj_string)
    rx1, ry1 = project_to_radar(links.site_1_lon, links.site_1_lat, grid.proj_string)
    rows, n_cells, coverage = [], np.zeros(n, dtype="int64"), np.zeros(n)
    for k in range(n):
        iy, ix, w = path_grid_weights(rx0[k], ry0[k], rx1[k], ry1[k], grid)
        n_cells[k], coverage[k] = iy.size, w.sum()
        rows.append(pd.DataFrame({"cml_id": links.cml_id[k], "iy": iy, "ix": ix, "weight": w}))
    weights = pd.concat(rows, ignore_index=True)
    weights[["iy", "ix"]] = weights[["iy", "ix"]].astype("int32")

    link_geometry = pd.DataFrame({
        "cml_id": links.cml_id.to_numpy(),
        "site_0_id": s0, "site_1_id": s1,
        "site_0_lon": links.site_0_lon, "site_0_lat": links.site_0_lat,
        "site_1_lon": links.site_1_lon, "site_1_lat": links.site_1_lat,
        "site_0_x_3006_m": x0, "site_0_y_3006_m": y0, "site_1_x_3006_m": x1, "site_1_y_3006_m": y1,
        "mid_x_3006_m": 0.5 * (x0 + x1), "mid_y_3006_m": 0.5 * (y0 + y1),
        "azimuth_deg": azimuth_deg(x0, y0, x1, y1),
        "length_m": length_m, "length_km_metadata": length_km_meta,
        "length_diff_m": length_diff_m, "length_rel_diff": length_rel_diff, "length_mismatch": mismatch,
        "frequency_ghz_sublink_1": links.frequency_ghz_sublink_1,
        "frequency_ghz_sublink_2": links.frequency_ghz_sublink_2,
        "polarization": links.polarization,
        "site_0_antenna_diameter_m": links.site_0_antenna_diameter_m,
        "site_1_antenna_diameter_m": links.site_1_antenna_diameter_m,
        "n_radar_cells": n_cells, "radar_coverage": coverage,
    })

    # --- gauges ---
    gauges = gauges.reset_index(drop=True)
    gx, gy = project_points(gauges.lon, gauges.lat, crs)
    # (n_gauges, n_links) matrices by broadcasting
    dist, pos = point_to_path(gx[:, None], gy[:, None], x0[None, :], y0[None, :], x1[None, :], y1[None, :])
    mx, my = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
    d_mid = np.hypot(gx[:, None] - mx[None, :], gy[:, None] - my[None, :])
    d_end = np.minimum(np.hypot(gx[:, None] - x0[None, :], gy[:, None] - y0[None, :]),
                       np.hypot(gx[:, None] - x1[None, :], gy[:, None] - y1[None, :]))
    gi, li = np.nonzero(dist <= gcfg["gauge_match_max_m"])
    matches = pd.DataFrame({
        "cml_id": links.cml_id.to_numpy()[li], "gauge_id": gauges.gauge_id.to_numpy()[gi],
        "network": gauges.network.to_numpy()[gi], "gauge_type": gauges.gauge_type.to_numpy()[gi],
        "quantization_mm": gauges.quantization_mm.to_numpy()[gi],
        "distance_to_path_m": dist[gi, li], "position_along_path": pos[gi, li],
        "distance_to_midpoint_m": d_mid[gi, li], "distance_to_nearest_end_m": d_end[gi, li],
    })
    for r in gcfg["gauge_flag_radii_m"]:
        matches[f"within_{int(r)}m"] = matches.distance_to_path_m <= r
    matches["link_length_m"] = length_m[li]
    matches["frequency_ghz"] = links.frequency_ghz_sublink_1.to_numpy()[li]
    matches = matches.sort_values(["cml_id", "distance_to_path_m", "gauge_id"]).reset_index(drop=True)

    # --- gauge radar pixel ---
    grx, gry = project_to_radar(gauges.lon, gauges.lat, grid.proj_string)
    giy, gix = grid.cell_index(grx, gry)
    inside = giy >= 0
    cx3006, cy3006 = grid_cells_3006
    dcell = np.full(len(gauges), np.nan)
    dcell[inside] = np.hypot(gx[inside] - cx3006[giy[inside], gix[inside]],
                             gy[inside] - cy3006[giy[inside], gix[inside]])
    gauge_pixel = pd.DataFrame({
        "gauge_id": gauges.gauge_id, "network": gauges.network, "lon": gauges.lon, "lat": gauges.lat,
        "x_3006_m": gx, "y_3006_m": gy, "iy": giy, "ix": gix,
        "distance_to_cell_centre_m": dcell, "inside_radar_grid": inside,
    })

    # --- neighbours ---
    lines = shapely.linestrings(np.stack([np.column_stack([x0, y0]), np.column_stack([x1, y1])], axis=1))
    ia, ib = np.triu_indices(n, k=1)
    dmin = shapely.distance(lines[ia], lines[ib])
    keep = dmin <= gcfg["neighbour_max_m"]
    ia, ib, dmin = ia[keep], ib[keep], dmin[keep]
    az = link_geometry.azimuth_deg.to_numpy()
    dang = np.mod(np.abs(az[ia] - az[ib]), 180.0)
    angle = np.minimum(dang, 180.0 - dang)  # undirected lines: 0..90
    shares = ((s0[ia] == s0[ib]) | (s0[ia] == s1[ib]) | (s1[ia] == s0[ib]) | (s1[ia] == s1[ib]))
    neighbours = pd.DataFrame({
        "cml_id_a": links.cml_id.to_numpy()[ia], "cml_id_b": links.cml_id.to_numpy()[ib],
        "min_path_distance_m": dmin, "midpoint_distance_m": np.hypot(mx[ia] - mx[ib], my[ia] - my[ib]),
        "shares_site": shares, "angle_between_deg": angle,
        "length_a_m": length_m[ia], "length_b_m": length_m[ib],
    })

    return {"sites": sites, "link_geometry": link_geometry, "path_weights": weights,
            "link_gauge_matches": matches, "gauge_radar_pixel": gauge_pixel, "link_neighbours": neighbours}


def convex_hull_area_m2(x, y) -> float:
    """Area of the convex hull of the points (m^2, in the CRS of x and y)."""
    return float(shapely.MultiPoint(np.column_stack([x, y])).convex_hull.area)


# ---------------------------------------------------------------------------
# Applying the weights to a radar field (used by the next step)
# ---------------------------------------------------------------------------
def apply_weights(field: xr.DataArray, weights: pd.DataFrame, threshold: float | None = None,
                  time_chunk: int = 288, grid: RadarGrid | None = None) -> xr.Dataset:
    """Path statistics of a (time, y, x) field for every link in `weights`.

    weights : DataFrame with cml_id, iy, ix, weight (as in path_grid_weights.nc).
    threshold : wet_fraction = share of the valid path weight in cells with value > threshold
              (None: wet_fraction is not computed). No threshold is applied to the field itself.
    time_chunk : number of time steps read at once. Only the bounding box of the crossed cells is read, so
              a lazily opened month file is never loaded whole.
    grid    : if given, the field's x / y must equal the grid the weights were computed on.

    Returns a Dataset on (time, cml_id) with:
      mean   weighted mean over the NON-NaN cells, sum(w v) / sum(w valid)
      std    weighted standard deviation over the same cells (population form)
      max    largest value of a crossed non-NaN cell
      wet_fraction           (if threshold is given) see above
      valid_weight_fraction  sum(w valid) / sum(w): share of the link's (in-grid) path weight that had data
    A time step where every crossed cell is NaN gives NaN statistics and valid_weight_fraction 0.
    """
    field = field.transpose(..., "time", "y", "x")
    if grid is not None:
        if not (np.allclose(field.x.values, grid.x) and np.allclose(field.y.values, grid.y)):
            raise ValueError("the field's x/y grid differs from the grid the weights were computed on")
    w_df = weights.sort_values("cml_id", kind="stable").reset_index(drop=True)
    cml_ids, starts = np.unique(w_df.cml_id.to_numpy(), return_index=True)
    link_of = np.searchsorted(cml_ids, w_df.cml_id.to_numpy())
    n_links = cml_ids.size
    w = w_df.weight.to_numpy(float)
    iy, ix = w_df.iy.to_numpy(), w_df.ix.to_numpy()
    # Read only the bounding box of all crossed cells.
    y0, y1, x0, x1 = iy.min(), iy.max() + 1, ix.min(), ix.max() + 1
    liy, lix = iy - y0, ix - x0
    onehot = np.zeros((w.size, n_links))  # entries x links, sums per link via one matrix product
    onehot[np.arange(w.size), link_of] = 1.0
    w_total = w @ onehot

    n_t = field.sizes["time"]
    out = {k: np.full((n_t, n_links), np.nan) for k in ("mean", "std", "max", "valid_weight_fraction")}
    if threshold is not None:
        out["wet_fraction"] = np.full((n_t, n_links), np.nan)
    for t0 in range(0, n_t, time_chunk):
        block = field.isel(time=slice(t0, t0 + time_chunk), y=slice(y0, y1), x=slice(x0, x1)).values
        v = block[:, liy, lix].astype("float64")  # (t, entries)
        valid = ~np.isnan(v)
        wv = w[None, :] * valid
        v0 = np.where(valid, v, 0.0)
        wsum = wv @ onehot
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = (wv * v0) @ onehot / wsum
            dev = np.where(valid, v0 - mean[:, link_of], 0.0)
            var = (wv * dev * dev) @ onehot / wsum
            sl = slice(t0, t0 + block.shape[0])
            out["mean"][sl] = mean
            out["std"][sl] = np.sqrt(np.maximum(var, 0.0))
            out["valid_weight_fraction"][sl] = wsum / w_total
            if threshold is not None:
                out["wet_fraction"][sl] = (wv * (v0 > threshold)) @ onehot / wsum
        # np.fmax ignores NaN; a group that is all NaN stays NaN
        out["max"][t0:t0 + block.shape[0]] = np.fmax.reduceat(v, starts, axis=1)
    ds = xr.Dataset({k: (("time", "cml_id"), a) for k, a in out.items()},
                    coords={"time": field.time.values, "cml_id": cml_ids})
    ds.attrs["comment"] = ("path statistics from cmlpipe.geometry.apply_weights; NaN cells excluded from mean/std/"
                           "max/wet_fraction and reported via valid_weight_fraction")
    if threshold is not None:
        ds["wet_fraction"].attrs["threshold"] = threshold
    return ds
