"""The `geometry` stage: read positions, call cmlpipe.geometry, write the tables.

Inputs (positions only, no time series):
  * link ends, lengths, frequencies, polarization, antenna diameters: raw cml_metadata.csv through the
    same adapter function the L1 ingest uses (openmrg.load_cml_metadata), so site_0 / site_1 are
    exactly the L1 convention.
  * gauges: tables/gauges.csv (written by ingest_gauges).
  * radar grid: x, y, x_3006, y_3006 and the proj string of the L1 radar month files.
Outputs (data/work/openmrg/...):
  tables/sites.csv, link_geometry.csv, link_gauge_matches.csv, gauge_radar_pixel.csv, link_neighbours.csv
  (each with a .meta.json sidecar) and L3_products/geometry/path_grid_weights.nc.

Skip rule: all outputs exist, carry the current geometry stage hash, AND were built from the same inputs
(input_signature = sha256 of the metadata CSV, gauges.csv and the radar grid arrays). --force rebuilds.
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr

from cmlpipe import geometry
from cmlpipe.io import openmrg
from cmlpipe.config import Config
from cmlpipe.provenance import provenance_attrs, write_sidecar

STAGE = "geometry"
TABLES = ["sites", "link_geometry", "link_gauge_matches", "gauge_radar_pixel", "link_neighbours"]

DESCRIPTIONS = {
    "sites": "unique physical sites: link ends merged within geometry.site_merge_tol_m (single linkage); "
             "lon/lat WGS84 and x/y EPSG:3006 are the mean of the merged ends",
    "link_geometry": "one row per link: sites, coordinates (WGS84, EPSG:3006), azimuth (0 = north, clockwise, "
                     "site_0 -> site_1), length computed in EPSG:3006 vs metadata, radio parameters, radar cells "
                     "crossed and radar_coverage (fraction of the path inside the radar grid)",
    "link_gauge_matches": "every (link, gauge) pair with gauge-to-PATH (segment) distance <= "
                          "geometry.gauge_match_max_m; position_along_path 0 = site_0, 1 = site_1",
    "gauge_radar_pixel": "radar cell (iy, ix into the radar L1 y / x) containing each gauge, located in the radar "
                         "projection; distance to the cell centre measured in EPSG:3006",
    "link_neighbours": "every pair of links whose PATHS come within geometry.neighbour_max_m (EPSG:3006); "
                       "angle_between_deg 0..90 for undirected lines",
}

# Rounding of the written CSVs: sub-millimetre digits are noise and make diffs hard to read.
ROUND = {"_m": 3, "lon": 7, "lat": 7, "azimuth_deg": 4, "angle_between_deg": 4, "length_rel_diff": 6,
         "position_along_path": 6}


def weights_path(cfg: Config) -> Path:
    """Path of path_grid_weights.nc."""
    return cfg.output_root / cfg["paths"]["l3_geometry_dir"] / "path_grid_weights.nc"


def output_paths(cfg: Config) -> dict[str, Path]:
    """All output paths of the stage."""
    paths = {name: cfg.tables_dir / f"{name}.csv" for name in TABLES}
    paths["path_weights"] = weights_path(cfg)
    return paths


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------
def link_inputs(cfg: Config) -> pd.DataFrame:
    """One row per link with the columns geometry.LINK_INPUT_COLUMNS, from the L1 metadata adapter."""
    meta = openmrg.load_cml_metadata(cfg)
    per = meta.pivot(index="cml_id", columns="sublink_id")
    pol = per["polarization"].replace("", np.nan)
    # polarization is a sub-link property in the raw file; both directions agree for OpenMRG, else say so
    differ = pol.notna().all(axis=1) & (pol["sublink_1"] != pol["sublink_2"])
    if differ.any():
        raise ValueError(f"links with different polarization per direction: {list(pol.index[differ])}")
    first = meta.drop_duplicates("cml_id").set_index("cml_id")
    out = pd.DataFrame({
        "cml_id": per.index.to_numpy(),
        "site_0_lon": first.site_0_lon, "site_0_lat": first.site_0_lat,
        "site_1_lon": first.site_1_lon, "site_1_lat": first.site_1_lat,
        "length_km_metadata": first.length,  # cml_metadata Length_km, the ONLY km value in this stage
        "frequency_ghz_sublink_1": per["frequency"]["sublink_1"],
        "frequency_ghz_sublink_2": per["frequency"]["sublink_2"],
        "polarization": pol["sublink_1"].fillna(pol["sublink_2"]),
        "site_0_antenna_diameter_m": first.site_0_antenna_diameter,
        "site_1_antenna_diameter_m": first.site_1_antenna_diameter,
    }).reset_index(drop=True)
    return out[geometry.LINK_INPUT_COLUMNS]


def gauge_inputs(cfg: Config) -> pd.DataFrame:
    """Gauges (gauge_id, network, type, lon, lat, quantization) from gauges.csv."""
    path = cfg.tables_dir / "gauges.csv"
    if not path.exists():
        raise FileNotFoundError(f"{path} not found: run ingest_gauges first")
    g = pd.read_csv(path).rename(columns={"id": "gauge_id"})
    return g[geometry.GAUGE_INPUT_COLUMNS]


def radar_grid_inputs(cfg: Config) -> tuple[geometry.RadarGrid, tuple[np.ndarray, np.ndarray], dict]:
    """Radar grid from the L1 month files (all must agree), checked against geometry.radar_proj_string."""
    files = sorted(cfg.l1_radar_dir.glob("radar_5min_*.nc"))
    if not files:
        raise FileNotFoundError(f"no radar L1 files in {cfg.l1_radar_dir}: run ingest_radar first")
    want = cfg["geometry"]["radar_proj_string"]
    ref = None
    for f in files:
        with xr.open_dataset(f) as ds:
            cur = {"x": ds.x.values, "y": ds.y.values, "x_3006": ds.x_3006.values, "y_3006": ds.y_3006.values,
                   "lon": ds.lon.values, "lat": ds.lat.values, "proj": ds.attrs.get("crs_proj4", "")}
        if cur["proj"].split() != want.split():
            raise ValueError(f"{f.name}: radar proj string {cur['proj']!r} differs from geometry.radar_proj_string")
        if ref is None:
            ref = cur
        elif not all(np.array_equal(ref[k], cur[k]) for k in ("x", "y", "x_3006", "y_3006")):
            raise ValueError(f"{f.name}: radar grid differs from {files[0].name}")
    grid = geometry.RadarGrid(ref["x"], ref["y"], want)
    return grid, (ref["x_3006"], ref["y_3006"]), {"files": [f.name for f in files], **ref}


def input_signature(cfg: Config, grid_info: dict) -> str:
    """sha256 over everything the stage reads: if it is unchanged, the outputs are too."""
    h = hashlib.sha256()
    h.update(cfg.raw_cml_metadata.read_bytes())
    h.update((cfg.tables_dir / "gauges.csv").read_bytes())
    for k in ("x", "y", "x_3006", "y_3006"):
        h.update(np.ascontiguousarray(grid_info[k], dtype="float64").tobytes())
    h.update(grid_info["proj"].encode())
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------
def _rounded(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for col in df.columns:
        if df[col].dtype.kind != "f":
            continue
        for suffix, nd in ROUND.items():
            if col.endswith(suffix):
                df[col] = df[col].round(nd)
                break
    return df


def write_tables(cfg: Config, tables: dict, grid: geometry.RadarGrid, signature: str) -> dict[str, Path]:
    """Write the five CSVs (+ sidecars) and the weights NetCDF. Content depends only on the inputs."""
    paths = output_paths(cfg)
    cfg.tables_dir.mkdir(parents=True, exist_ok=True)
    extra_common = {"input_signature": signature, "units": "all *_m columns in metres (EPSG:3006 for distances)"}
    for name in TABLES:
        df = _rounded(tables[name])
        df.to_csv(paths[name], index=False)
        write_sidecar(paths[name], cfg, STAGE, extra={**extra_common, "description": DESCRIPTIONS[name],
                                                       "n_rows": len(df), "columns": list(df.columns)})

    w = tables["path_weights"]
    ds = xr.Dataset(
        {"cml_id": ("entry", w.cml_id.to_numpy("int64")), "iy": ("entry", w.iy.to_numpy("int32")),
         "ix": ("entry", w.ix.to_numpy("int32")), "weight": ("entry", w.weight.to_numpy("float64"))})
    ds["iy"].attrs = {"long_name": "row index into the radar L1 'y' coordinate (0-based)"}
    ds["ix"].attrs = {"long_name": "column index into the radar L1 'x' coordinate (0-based)"}
    ds["weight"].attrs = {"long_name": "fraction of the link path length inside the cell", "units": "1",
                          "comment": "sum over a link = radar_coverage in link_geometry.csv"}
    ds.attrs = {
        **provenance_attrs(cfg, STAGE),
        "title": "Link path / radar cell intersection weights (sparse long format)",
        "method": "exact edge crossings of the straight path in the radar projection; each piece between two "
                  "crossings is assigned to the cell of its midpoint (no double counting on edges/corners)",
        "grid_proj_string": grid.proj_string,
        "grid_x_centres_m": grid.x, "grid_y_centres_m": grid.y,
        "grid_dx_m": grid.dx, "grid_dy_m": grid.dy,
        "input_signature": signature,
    }
    paths["path_weights"].parent.mkdir(parents=True, exist_ok=True)
    tmp = paths["path_weights"].with_suffix(".nc.tmp")
    ds.to_netcdf(tmp, format="NETCDF4")
    tmp.replace(paths["path_weights"])
    return paths


def read_weights(path: Path) -> tuple[pd.DataFrame, geometry.RadarGrid]:
    """path_grid_weights.nc -> (DataFrame cml_id, iy, ix, weight; the grid it refers to)."""
    with xr.open_dataset(path) as ds:
        df = pd.DataFrame({k: ds[k].values for k in ("cml_id", "iy", "ix", "weight")})
        grid = geometry.RadarGrid(np.atleast_1d(ds.attrs["grid_x_centres_m"]),
                                  np.atleast_1d(ds.attrs["grid_y_centres_m"]), ds.attrs["grid_proj_string"])
    return df, grid


def _existing(paths: dict[str, Path]) -> list[tuple[str, str]] | None:
    """(stage_hash, input_signature) of every output, or None if one is missing."""
    found = []
    for name, p in paths.items():
        if not p.exists():
            return None
        if p.suffix == ".nc":
            with netCDF4.Dataset(p) as nc:
                found.append((nc.getncattr("stage_hash"), nc.getncattr("input_signature")))
        else:
            side = Path(str(p) + ".meta.json")
            if not side.exists():
                return None
            m = json.loads(side.read_text(encoding="utf-8"))
            found.append((m.get("stage_hash"), m.get("input_signature")))
    return found


def run_geometry(cfg: Config, force: bool, log, records: list[dict]) -> None:
    """Run the stage (skips when outputs are up to date)."""
    t0 = time.perf_counter()
    gcfg = cfg["geometry"]
    if gcfg["crs_metric"] != cfg["reference"]["metric_crs"]:
        raise ValueError("geometry.crs_metric must equal reference.metric_crs (used for radar x_3006/y_3006)")
    grid, cells3006, info = radar_grid_inputs(cfg)
    signature = input_signature(cfg, info)
    paths = output_paths(cfg)
    existing = _existing(paths)
    if not force and existing and set(existing) == {(cfg.stage_hash(STAGE), signature)}:
        log.info("geometry: outputs up to date (same stage hash and inputs), skipped (use --force to rebuild)")
        records.append({"step": STAGE, "status": "skipped_same_config"})
        return
    tables = geometry.build_geometry(link_inputs(cfg), gauge_inputs(cfg), grid, cells3006, gcfg)
    written = write_tables(cfg, tables, grid, signature)
    lg = tables["link_geometry"]
    record = {"step": STAGE, "status": "written", "input_signature": signature,
              "n_links": len(lg), "n_sites": len(tables["sites"]), "n_weight_entries": len(tables["path_weights"]),
              "n_gauge_matches": len(tables["link_gauge_matches"]), "n_neighbour_pairs": len(tables["link_neighbours"]),
              "n_length_mismatch": int(lg.length_mismatch.sum()),
              "files": [str(p.relative_to(cfg.repo_root)) if p.is_relative_to(cfg.repo_root) else str(p)
                        for p in written.values()],
              "runtime_s": round(time.perf_counter() - t0, 2)}
    log.info("geometry: %d links, %d sites, %d weight entries, %d gauge matches, %d neighbour pairs (%.1f s)",
             record["n_links"], record["n_sites"], record["n_weight_entries"], record["n_gauge_matches"],
             record["n_neighbour_pairs"], record["runtime_s"])
    records.append(record)
