"""Stage `release`: build the read-only, versioned data release data/release/openmrg_<version>/.

    python -m cmlpipe run release                 # builds openmrg_v1 (skips if it is complete and up to date)
    python -m cmlpipe run release --version v2    # a new version next to v1
    python -m cmlpipe run release --force         # rebuild openmrg_v1 (it is read-only otherwise)

Contents: tables/ (Parquet + CSV, consistent keys), timeseries/ (Zarr, see release.metadata.TIMESERIES for dims and
chunks), catalog.json (every product, variable, unit, time span, stage hashes, sha256 of the tables), DATA_DICTIONARY.md
and README.md (generated from the same metadata), verification.json, DATA_GUIDE.md (release.guide), LICENSE.txt (CC BY-SA 4.0, release.metadata).

The release is built into <name>.building, verified (20 random (sub-link, day, resolution) samples equal the work
NetCDF values, table row counts equal the work tables, frozen tables unchanged), then renamed and made read-only.
Nothing in the work folder is changed.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import stat
import time
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import xarray as xr

import cmlpipe
from cmlpipe.aggregate import stage as aggregate_stage
from cmlpipe.config import STAGE_CONFIG_KEYS, Config
from cmlpipe.events import link_state as link_state_stage
from cmlpipe.events import stage as events_stage
from cmlpipe.io.netcdf import file_attrs
from cmlpipe.provenance import git_commit, utc_now_iso
from cmlpipe.qc import stage as qc_stage
from cmlpipe.qc.flags import QCFlag
from cmlpipe.release import cmlstores, metadata, zarrstore


def release_dir(cfg: Config, version: str) -> Path:
    """data/release/openmrg_<version> (paths.release_root, relative to the repository root)."""
    return cfg.repo_path(cfg["paths"]["release_root"]) / f"{metadata.RELEASE_NAME}_{version}"


def sha256(path: Path) -> str:
    """sha256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------
def timeseries_sources(cfg: Config) -> dict:
    """Release product -> work-folder source file(s) (NetCDF)."""
    R, G = aggregate_stage.radar_path, aggregate_stage.gauge_path
    src = {f"cml_{r}": [aggregate_stage.cml_path(cfg, r, d) for d in cfg.period_days()] for r in cfg["aggregate"]["resolutions"]}
    src["cml_10s"] = [cfg.l1_cml_day_path(d) for d in cfg.period_days()] + [qc_stage.l2_day_path(cfg, d) for d in cfg.period_days()]
    for r in cfg["aggregate"]["radar"]["resolutions"]:
        src[f"radar_along_links_{r}"] = [R(cfg, "along_links", r)]
        src[f"radar_at_gauges_{r}"] = [R(cfg, "at_gauges", r)]
    for r in cfg["aggregate"]["gauges"]["city_resolutions"]:
        src[f"gauges_city_{r}"] = [G(cfg, "city", r)]
    for r in cfg["aggregate"]["gauges"]["smhi_resolutions"]:
        src[f"gauges_smhi_{r}"] = [G(cfg, "smhi", r)]
    src["met_1h"] = [cfg.l1_met_dir / "smhi_met_1h.nc"]
    src["network_reference_5min"] = [events_stage.reference_path(cfg)]
    src["link_reference_state_5min"] = [link_state_stage.state_path(cfg)]
    return src


def input_signature(cfg: Config) -> str:
    """sha256 over the work tables (bytes) and the stage hash / creation time of every source NetCDF file."""
    h = hashlib.sha256()
    for name, (rel, _, _) in metadata.TABLES.items():
        h.update(name.encode())
        h.update((cfg.output_root / rel).read_bytes())
    for name, files in sorted(timeseries_sources(cfg).items()):
        for f in files:
            a = file_attrs(f)
            h.update(f"{name}{f.name}{a.get('stage_hash')}{a.get('created_utc')}".encode())
    h.update(json.dumps({k: list(v[1:]) for k, v in metadata.TIMESERIES.items()}, default=str).encode())
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------
def normalise(name: str, df: pd.DataFrame) -> pd.DataFrame:
    """Consistent keys: sublink_id '<cml_id>_<1|2>' (+ int column sublink), gauge_id instead of id."""
    df = df.copy()
    if "sublink_id" in df and df.sublink_id.astype(str).str.startswith("sublink_").any():
        n = df.sublink_id.astype(str).str[-1].astype(int)
        df["sublink_id"] = df.cml_id.astype(int).astype(str) + "_" + n.astype(str)
        if "sublink" not in df:
            df.insert(df.columns.get_loc("sublink_id") + 1, "sublink", n.astype("int8"))
    if name == "gauges":
        df = df.rename(columns={"id": "gauge_id"})
    for c in ("start", "end", "wet_core_start", "wet_core_end"):
        if c in df:
            df[c] = pd.to_datetime(df[c])
    return df


def write_tables(cfg: Config, out: Path, log) -> dict:
    """tables/<name>.parquet and .csv; returns catalog entries. Frozen tables: the CSV is a byte copy of the work file."""
    (out / "tables").mkdir(parents=True, exist_ok=True)
    cat = {}
    for name, (rel, keys, desc) in metadata.TABLES.items():
        src = cfg.output_root / rel
        df = normalise(name, pd.read_csv(src))
        pq, csv = out / "tables" / f"{name}.parquet", out / "tables" / f"{name}.csv"
        df.to_parquet(pq, index=False)
        if name in metadata.FROZEN_TABLES:
            shutil.copyfile(src, csv)   # byte-identical copy of the frozen table
        else:
            df.to_csv(csv, index=False, date_format="%Y-%m-%d %H:%M:%S")
        cols = [{"name": c, "dtype": str(df[c].dtype), "unit": metadata.column_unit(c),
                 "description": metadata.column_description(name, c)} for c in df.columns]
        cat[name] = {"kind": "table", "description": desc, "keys": [k if k != "id" else "gauge_id" for k in keys],
                     "n_rows": len(df), "paths": {"parquet": f"tables/{name}.parquet", "csv": f"tables/{name}.csv"},
                     "sha256": {"parquet": sha256(pq), "csv": sha256(csv)}, "columns": cols,
                     "source": rel, "source_stage_hash": json.loads(Path(str(src) + ".meta.json").read_text(
                         encoding="utf-8")).get("stage_hash") if Path(str(src) + ".meta.json").exists() else None,
                     "frozen": name in metadata.FROZEN_TABLES}
        log.info("release: table %s (%d rows)", name, len(df))
    return cat


# ---------------------------------------------------------------------------
# Time series
# ---------------------------------------------------------------------------
def _common_attrs(cfg: Config, name: str, src_files: list[Path]) -> dict:
    a0 = file_attrs(src_files[0])
    return {"title": f"{metadata.RELEASE_NAME} release: {name}", "description": metadata.TIMESERIES[name][0],
            "time_convention": "interval_start", "time_zone": "UTC", "interval": metadata.TIMESERIES[name][1],
            "source_stage": a0.get("stage", ""), "source_stage_hash": a0.get("stage_hash", ""),
            "code_version": git_commit(cfg.repo_root), "cmlpipe_version": cmlpipe.__version__,
            "dataset_doi": cfg["dataset"]["source_doi"], "created_utc": utc_now_iso()}


def write_netcdf_product(cfg: Config, name: str, src: Path, out: Path) -> None:
    """A product that is one work NetCDF file: copied to Zarr (dims renamed id -> gauge_id)."""
    desc, res, dims, chunks = metadata.TIMESERIES[name]
    with xr.open_dataset(src) as ds:
        ds = ds.load()
    if "id" in ds.dims:
        ds = ds.rename({"id": "gauge_id"})
    keep = ("units", "long_name", "flag_masks", "flag_meanings", "flag_values", "comment")
    attrs = {**_common_attrs(cfg, name, [src])}
    for k, v in ds.attrs.items():  # keep the source description (rules, thresholds) under source_*
        if k in ("time_meaning", "valid_rule", "floor_rule", "rules", "definition", "valid_stamp_rule", "min_dbz_floor",
                 "floor_share_of_cell_steps_set_to_zero", "wet_threshold_mmh", "min_valid_fraction", "n_expected"):
            attrs[f"source_{k}"] = v if isinstance(v, (str, int, float)) else str(v)
    enc = {}
    for v in ds.variables:
        ds[v].attrs = {k: (x.tolist() if isinstance(x, np.ndarray) else x) for k, x in ds[v].attrs.items() if k in keep}
        ds[v].encoding = {}  # drop NetCDF-specific encodings (zlib, contiguous, source); xarray re-encodes for Zarr
        if v in ds.data_vars and chunks is not None and ds[v].dims == dims:
            enc[v] = {"chunks": chunks}
    enc["time"] = {"units": "seconds since 1970-01-01 00:00:00", "dtype": "int64"}
    ds.attrs = attrs
    ds.to_zarr(out, mode="w", zarr_format=2, consolidated=True, encoding=enc)


def write_timeseries(cfg: Config, out: Path, log) -> dict:
    """timeseries/<name>.zarr for every product; returns catalog entries."""
    ts = out / "timeseries"
    ts.mkdir(parents=True, exist_ok=True)
    sources = timeseries_sources(cfg)
    cat = {}
    for name, (desc, res, dims, chunks) in metadata.TIMESERIES.items():
        t0 = time.perf_counter()
        p = ts / f"{name}.zarr"
        attrs = _common_attrs(cfg, name, sources[name])
        if name == "cml_10s":
            cmlstores.write_cml_10s(cfg, p, chunks, attrs=attrs)
        elif name.startswith("cml_"):
            cmlstores.write_cml_res(cfg, p, res, chunks, attrs=attrs)
        else:
            write_netcdf_product(cfg, name, sources[name][0], p)
        ds = zarrstore.open_store(p)
        t = pd.DatetimeIndex(ds.time.values)
        cat[name] = {
            "kind": "timeseries", "description": desc, "path": f"timeseries/{name}.zarr", "format": "zarr (format 2, "
            "consolidated, blosc-zstd)", "resolution": res, "dims": {d: int(ds.sizes[d]) for d in ds.dims},
            "chunks": list(chunks) if chunks else "one chunk per variable", "time_start": str(t[0]),
            "time_end_exclusive": str(t[-1] + pd.Timedelta(res)), "time_convention": "interval_start (UTC)",
            "variables": {v: {"dims": list(ds[v].dims), "dtype": str(ds[v].dtype), "units": ds[v].attrs.get("units", ""),
                              "description": metadata.variable_description(v, ds[v].attrs)} for v in ds.data_vars},
            "coords": {c: {"dims": list(ds[c].dims), "dtype": str(ds[c].dtype)} for c in ds.coords if c != "time"},
            "source_stage_hash": attrs["source_stage_hash"], "source_files": len(sources[name]),
            "size_bytes": zarrstore.store_size(p),
        }
        ds.close()
        log.info("release: %s (%.0f MB, %.0f s)", name, cat[name]["size_bytes"] / 1e6, time.perf_counter() - t0)
    return cat


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------
def verify(cfg: Config, out: Path, catalog: dict, frozen_before: dict, seed: int = 8) -> dict:
    """20 random (sub-link, day, resolution) samples: Zarr == work NetCDF (ints exactly, floats to 1e-6);
    table row counts; frozen work tables unchanged."""
    rng = np.random.default_rng(seed)
    days = cfg.period_days()
    res_list = ["10s"] + list(cfg["aggregate"]["resolutions"])
    samples, worst = [], 0.0
    stores = {r: zarrstore.open_store(out / "timeseries" / f"cml_{r}.zarr") for r in res_list}
    for _ in range(20):
        r = res_list[rng.integers(len(res_list))]
        d = days[rng.integers(len(days))]
        ds = stores[r]
        sid = str(ds.sublink_id.values[rng.integers(ds.sizes["sublink_id"])])
        cml, k = int(sid.split("_")[0]), int(sid.split("_")[1]) - 1
        z = ds.sel(sublink_id=sid, time=slice(d, d + pd.Timedelta("1D") - pd.Timedelta("1s")))
        files = [(cfg.l1_cml_day_path(d), ("tsl", "rsl", "tl")), (qc_stage.l2_day_path(cfg, d), ("flags",))] if r == "10s" \
            else [(aggregate_stage.cml_path(cfg, r, d), tuple(v for v in ds.data_vars))]
        ok = True
        for f, names in files:
            with netCDF4.Dataset(f) as nc:
                j = int(np.flatnonzero(nc["cml_id"][:] == cml)[0])
                for v in names:
                    a = np.ma.filled(nc[v][:, k, j].astype("float64"), np.nan)
                    b = z[v].values.astype("float64")
                    if np.issubdtype(nc[v].dtype, np.integer) and "scale_factor" not in nc[v].ncattrs():
                        same = np.array_equal(a, b)
                    else:
                        diff = np.nanmax(np.abs(a - b)) if np.isfinite(a).any() else 0.0
                        worst = max(worst, float(diff))
                        same = np.array_equal(np.isnan(a), np.isnan(b)) and diff <= 1e-6
                    ok &= bool(same)
        samples.append({"resolution": r, "day": f"{d:%Y-%m-%d}", "sublink_id": sid, "equal": ok})
    rows = {n: {"release": e["n_rows"], "work": len(pd.read_csv(cfg.output_root / metadata.TABLES[n][0]))}
            for n, e in catalog["tables"].items()}
    frozen_now = {n: sha256(cfg.output_root / metadata.TABLES[n][0]) for n in metadata.FROZEN_TABLES}
    frozen_copy = {n: catalog["tables"][n]["sha256"]["csv"] for n in metadata.FROZEN_TABLES}
    return {"samples": samples, "all_samples_equal": all(s["equal"] for s in samples), "max_float_diff": worst,
            "row_counts": rows, "row_counts_equal": all(v["release"] == v["work"] for v in rows.values()),
            "frozen_sha256_work": frozen_now, "frozen_unchanged": frozen_now == frozen_before,
            "frozen_release_csv_equal_work": frozen_copy == frozen_now}


# ---------------------------------------------------------------------------
# Docs generated from the catalog
# ---------------------------------------------------------------------------
def data_dictionary(catalog: dict) -> str:
    """DATA_DICTIONARY.md from the catalog (single source: release.metadata + the files' own attributes)."""
    L = [f"# Data dictionary: {catalog['name']} {catalog['version']}\n",
         f"Generated {catalog['created_utc']} from catalog.json (code {catalog['code_version']}). Conventions:\n"]
    L += [f"- **{k}**: {v}" for k, v in catalog["conventions"].items()]
    L.append("\n## Tables (tables/<name>.parquet and .csv)\n")
    for n, e in catalog["tables"].items():
        L.append(f"### {n}\n\n{e['description']}. Keys: {', '.join(e['keys'])}. Rows: {e['n_rows']}."
                 + (" **Frozen.**" if e["frozen"] else "") + "\n")
        L.append("| column | type | unit | description |\n|---|---|---|---|")
        L += [f"| {c['name']} | {c['dtype']} | {c['unit']} | {c['description']} |" for c in e["columns"]]
        L.append("")
    L.append("## Time series (timeseries/<name>.zarr)\n")
    for n, e in catalog["timeseries"].items():
        L.append(f"### {n}\n\n{e['description']}. Resolution {e['resolution']}, dims "
                 + ", ".join(f"{d}={s}" for d, s in e["dims"].items()) + f", {e['time_start']} to "
                 f"{e['time_end_exclusive']} (exclusive), chunks {e['chunks']}.\n")
        L.append("| variable | dims | type | unit | description |\n|---|---|---|---|---|")
        L += [f"| {v} | {', '.join(m['dims'])} | {m['dtype']} | {m['units']} | {m['description']} |" for v, m in e["variables"].items()]
        L.append("")
    return "\n".join(L) + "\n"


def readme(catalog: dict) -> str:
    """README.md of the release."""
    return f"""# {catalog['name']} {catalog['version']}: OpenMRG commercial microwave links, processed

Processed, quality-controlled and aggregated OpenMRG data (Gothenburg, 1 June - 31 August 2015): 364 links / 728
sub-links at 10 s, their statistics at 1 / 5 / 15 / 60 min, radar along every link path, City and SMHI gauges,
weather, rain events, dry periods, a frozen train / test split, per-link reference states and an inventory with a
quality status. Built by cmlpipe {catalog['cmlpipe_version']} (code {catalog['code_version']}) on {catalog['created_utc']}.
Read-only: rebuild only with `python -m cmlpipe run release --force` (or build `--version v2`).
**Licence:** {metadata.DATA_LICENSE_SHORT}.

## Open it

```python
from cmlpipe.api import CMLData
d = CMLData.open("data/release/{catalog['name']}_{catalog['version']}")
df = d.frame(10227, "2015-07-28 12:00", "2015-07-28 20:00", res="5min")
```

Without cmlpipe: tables are Parquet / CSV (`pandas.read_parquet`), time series are Zarr
(`xarray.open_dataset(path, engine="zarr", chunks=None, consolidated=True)`).

## Conventions

""" + "\n".join(f"- **{k}**: {v}" for k, v in catalog["conventions"].items()) + f"""

Every product, variable, unit and time span is listed in `catalog.json` and `DATA_DICTIONARY.md`; `verification.json`
holds the checks run when the release was built. Methods: docs/methods.md of the cmlpipe repository.

## Citation

Data: {metadata.CITATION}
Processing: Haber, Y. and Tsitseluk, Y., cmlpipe (CellEnMon lab), {catalog['name']} release {catalog['version']}.
"""


# ---------------------------------------------------------------------------
# Stage
# ---------------------------------------------------------------------------
def _make_writable(path: Path) -> None:
    for p in list(path.rglob("*")) + [path]:
        try:
            os.chmod(p, stat.S_IWRITE | stat.S_IREAD)
        except OSError:
            pass


def _make_readonly(path: Path) -> None:
    for p in path.rglob("*"):
        if p.is_file():
            os.chmod(p, stat.S_IREAD)


def run_release(cfg: Config, force: bool, log, records: list[dict], version: str = "v1") -> None:
    """Build (or skip) the release version `version` (see module docstring)."""
    t0 = time.perf_counter()
    out = release_dir(cfg, version)
    sig = input_signature(cfg)
    cat_path = out / "catalog.json"
    if out.exists() and not force:
        if cat_path.exists() and json.loads(cat_path.read_text(encoding="utf-8")).get("input_signature") == sig:
            log.info("release: %s is complete and up to date (read-only), skipped", out.name)
            records.append({"step": "release", "status": "skipped_same_config", "path": str(out)})
            return
        raise SystemExit(f"release: {out} exists and is read-only; its inputs changed. Rebuild with --force, or build a "
                         f"new version with --version v{int(version.lstrip('v')) + 1 if version.lstrip('v').isdigit() else 2}.")
    frozen_before = {n: sha256(cfg.output_root / metadata.TABLES[n][0]) for n in metadata.FROZEN_TABLES}
    tmp = out.with_name(out.name + ".building")
    if tmp.exists():
        _make_writable(tmp)
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    catalog = {"name": metadata.RELEASE_NAME, "version": version, "created_utc": utc_now_iso(),
               "code_version": git_commit(cfg.repo_root), "cmlpipe_version": cmlpipe.__version__,
               "dataset_doi": cfg["dataset"]["source_doi"], "citation": metadata.CITATION,
               "conventions": metadata.CONVENTIONS, "input_signature": sig,
               "qc": {"blocking_flags_rsl": cfg["aggregate"]["blocking_flags"]["rsl"],
                      "blocking_flags_tsl": cfg["aggregate"]["blocking_flags"]["tsl"],
                      "blocking_flags_l2": cfg["qc"]["blocking_flags"],
                      "flag_bits": {f.name: int(f) for f in QCFlag}},
               "stage_hashes": {s: cfg.stage_hash(s) for s in STAGE_CONFIG_KEYS},
               "period": {"start": str(cfg.period_start), "end_exclusive": str(cfg.period_end)}}
    catalog["tables"] = write_tables(cfg, tmp, log)
    catalog["timeseries"] = write_timeseries(cfg, tmp, log)
    check = verify(cfg, tmp, catalog, frozen_before)
    (tmp / "verification.json").write_text(json.dumps(check, indent=2, default=str), encoding="utf-8")
    if not (check["all_samples_equal"] and check["row_counts_equal"] and check["frozen_unchanged"]
            and check["frozen_release_csv_equal_work"]):
        raise RuntimeError(f"release verification failed, see {tmp / 'verification.json'} (release not published)")
    catalog["verification"] = {k: v for k, v in check.items() if k != "samples"}
    cat_path_tmp = tmp / "catalog.json"
    cat_path_tmp.write_text(json.dumps(catalog, indent=2, default=str), encoding="utf-8")
    dd = data_dictionary(catalog)
    (tmp / "DATA_DICTIONARY.md").write_text(dd, encoding="utf-8")
    (cfg.repo_root / "docs" / "data_dictionary.md").write_text(dd, encoding="utf-8")
    (tmp / "README.md").write_text(readme(catalog), encoding="utf-8")
    (tmp / "LICENSE.txt").write_text(metadata.data_license_text(catalog["name"], version), encoding="utf-8")
    if out.exists():  # only reached with --force
        _make_writable(out)
        shutil.rmtree(out)
    tmp.rename(out)
    _make_readonly(out)
    from cmlpipe.release import guide
    guide.write_guide(cfg, out, log)  # DATA_GUIDE.md (release root + docs/), generated from the catalog
    size = zarrstore.store_size(out)
    rec = {"step": "release", "status": "written", "path": str(out), "size_bytes": size,
           "runtime_s": round(time.perf_counter() - t0, 2), "verification": catalog["verification"]}
    log.info("release: wrote %s (%.2f GB, %.0f s); verification passed", out.name, size / 1e9, rec["runtime_s"])
    records.append(rec)
