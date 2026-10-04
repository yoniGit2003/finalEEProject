"""The common L1 formats (CML, gauges, weather, radar) and their validators.

The CML part is described below; the reference formats (gauges, weather, radar) and
their time convention are described where they are defined, further down.

L1 ("standard") is the first processed level: raw values, unchanged except for
  - the missing-value sentinel turned into NaN,
  - reshaping to (time, sublink_id, cml_id) and a regular time grid,
  - tl = tsl - rsl.
It follows the OpenSense naming (Fencl et al. 2023) so pycomlink and poligrain
can use it, with OUR unit choices: frequency in GHz and length in km. (The
OpenSense sample files use MHz and m; convert when a library expects those.)

Layout
------
dims       (time, sublink_id, cml_id), sublink_id = "sublink_1", "sublink_2"
data vars  tsl, rsl [dBm], tl [dB]                    on (time, sublink_id, cml_id)
coords     cml_id                                     link id (raw metadata "Link")
           raw_sublink_id, direction                  raw "Sublink" and "Direction" (traceability)
           tx_site, rx_site                           0 or 1: which site transmits / receives
           frequency [GHz], polarization ("v"/"h")    per sub-link
           length [km]                                per link
           site_0_lat/lon, site_1_lat/lon [deg]       per link (WGS84)
           site_0/1_antenna_diameter [m]              per link
A sub-link that is absent from the raw metadata has raw_sublink_id = -1, empty
direction and polarization, NaN frequency and NaN data.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import xarray as xr

DIMS = ("time", "sublink_id", "cml_id")
SUBLINK_IDS = ("sublink_1", "sublink_2")
MISSING_SUBLINK = -1  # raw_sublink_id value for a sub-link absent from the metadata

# name -> (units, long_name)
DATA_VARS = {
    "tsl": ("dBm", "transmitted signal level"),
    "rsl": ("dBm", "received signal level"),
    "tl": ("dB", "total loss, tsl - rsl"),
}

# name -> (dims, units, long_name). units = None marks identifiers and labels,
# which have no physical unit; for those we require a long_name instead.
COORDS = {
    "cml_id": (("cml_id",), None, "CML (link) id, = raw metadata 'Link'"),
    "sublink_id": (("sublink_id",), None, "sub-link (direction) label"),
    "raw_sublink_id": (("sublink_id", "cml_id"), None, "raw 'Sublink' number (-1 = absent)"),
    "direction": (("sublink_id", "cml_id"), None, "raw 'Direction' code A/B"),
    "tx_site": (("sublink_id", "cml_id"), None, "transmitting site (0 or 1)"),
    "rx_site": (("sublink_id", "cml_id"), None, "receiving site (0 or 1)"),
    "frequency": (("sublink_id", "cml_id"), "GHz", "carrier frequency"),
    "polarization": (("sublink_id", "cml_id"), None, "polarization, 'v' or 'h'"),
    "length": (("cml_id",), "km", "link length"),
    "site_0_lat": (("cml_id",), "degrees_north", "latitude of site 0 (WGS84)"),
    "site_0_lon": (("cml_id",), "degrees_east", "longitude of site 0 (WGS84)"),
    "site_1_lat": (("cml_id",), "degrees_north", "latitude of site 1 (WGS84)"),
    "site_1_lon": (("cml_id",), "degrees_east", "longitude of site 1 (WGS84)"),
    "site_0_antenna_diameter": (("cml_id",), "m", "antenna diameter at site 0"),
    "site_1_antenna_diameter": (("cml_id",), "m", "antenna diameter at site 1"),
}

# Used when the validator is called without rules from a config (same values as openmrg.yaml).
DEFAULT_RULES = {
    "frequency_ghz": [1.0, 100.0],
    "length_km": [0.01, 50.0],
    "lat_deg": [-90.0, 90.0],
    "lon_deg": [-180.0, 180.0],
    "rsl_plausible_dbm": [-110.0, -10.0],
    "tsl_plausible_dbm": [-20.0, 40.0],
    "tl_tolerance_db": 1.0e-3,
}


class CMLValidationError(ValueError):
    """Raised when a dataset does not follow the L1 CML schema."""


@dataclass
class ValidationReport:
    """Result of validate_cml_l1: hard errors and soft warnings (plain text)."""

    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """True when there are no errors."""
        return not self.errors

    def raise_if_errors(self, what: str = "L1 CML") -> None:
        """Raise CMLValidationError listing the errors, if any."""
        if self.errors:
            raise CMLValidationError(f"{what} validation failed:\n  - " + "\n  - ".join(self.errors))


def _out_of_range(values: np.ndarray, lo: float, hi: float) -> int:
    """Number of finite values outside [lo, hi]. NaN is not counted."""
    v = np.asarray(values, dtype="float64")
    return int(np.sum(np.isfinite(v) & ((v < lo) | (v > hi))))


def validate_cml_l1(
    ds: xr.Dataset,
    rules: dict | None = None,
    interval: str | pd.Timedelta | None = "10s",
    raise_on_error: bool = True,
) -> ValidationReport:
    """Check that `ds` follows the L1 CML schema.

    Parameters
    ----------
    ds : the dataset to check.
    rules : thresholds (the `validation` section of the config). Default: DEFAULT_RULES.
    interval : expected time step, e.g. "10s". None skips the regularity check.
    raise_on_error : raise CMLValidationError if there are errors (default True).

    Returns
    -------
    ValidationReport with `errors` (schema broken) and `warnings` (data suspicious,
    e.g. rsl outside the plausible range; those values are kept, not removed).
    """
    rules = {**DEFAULT_RULES, **(rules or {})}
    rep = ValidationReport()

    # ---- 1. dimensions ------------------------------------------------------
    for name, (units, _) in DATA_VARS.items():
        if name not in ds.data_vars:
            rep.errors.append(f"data variable '{name}' is missing")
            continue
        if ds[name].dims != DIMS:
            rep.errors.append(f"'{name}' has dims {ds[name].dims}, expected {DIMS}")
        if ds[name].attrs.get("units") != units:
            rep.errors.append(f"'{name}' units = {ds[name].attrs.get('units')!r}, expected {units!r}")
    if "sublink_id" in ds.coords and tuple(ds.sublink_id.values.tolist()) != SUBLINK_IDS:
        rep.errors.append(f"sublink_id values are {ds.sublink_id.values.tolist()}, expected {list(SUBLINK_IDS)}")

    # ---- 2. coordinates: present, right dims, units (or long_name) -----------
    for name, (dims, units, _) in COORDS.items():
        if name not in ds.coords:
            rep.errors.append(f"coordinate '{name}' is missing")
            continue
        if ds[name].dims != dims:
            rep.errors.append(f"coordinate '{name}' has dims {ds[name].dims}, expected {dims}")
        if units is not None and ds[name].attrs.get("units") != units:
            rep.errors.append(f"coordinate '{name}' units = {ds[name].attrs.get('units')!r}, expected {units!r}")
        if units is None and not ds[name].attrs.get("long_name"):
            rep.errors.append(f"coordinate '{name}' has no long_name")

    # ---- 3. time axis ---------------------------------------------------------
    if "time" not in ds.coords:
        rep.errors.append("coordinate 'time' is missing")
    elif ds.time.dtype.kind != "M":
        # tz-aware times end up as dtype object in xarray, strings as '<U'.
        rep.errors.append(f"time has dtype {ds.time.dtype}; expected tz-naive datetime64 (UTC)")
    else:
        t = ds.time.values.astype("datetime64[ns]").astype("int64")
        steps = np.diff(t)
        if len(np.unique(t)) != len(t):
            rep.errors.append(f"time has {len(t) - len(np.unique(t))} duplicate values")
        if np.any(steps <= 0):
            rep.errors.append(f"time is not strictly increasing ({int(np.sum(steps <= 0))} bad steps)")
        if interval is not None and len(steps) > 0:
            expected = pd.Timedelta(interval).value
            n_irregular = int(np.sum(steps != expected))
            if n_irregular:
                rep.errors.append(f"time is not regular at {interval}: {n_irregular} irregular steps")

    # Stop here if the structure is broken: value checks would only add noise.
    if rep.errors:
        if raise_on_error:
            rep.raise_if_errors()
        return rep

    # ---- 4. metadata ranges -----------------------------------------------------
    present = ds.raw_sublink_id.values != MISSING_SUBLINK  # sub-links that exist in the raw metadata
    freq = ds.frequency.values.astype("float64")
    lo, hi = rules["frequency_ghz"]
    if np.any(~np.isfinite(freq[present])):
        rep.errors.append("frequency is NaN for a sub-link that exists")
    n_bad = _out_of_range(freq, lo, hi)
    if n_bad:
        rep.errors.append(f"{n_bad} frequencies outside [{lo}, {hi}] GHz (wrong unit, e.g. MHz?)")

    lo, hi = rules["length_km"]
    length = ds.length.values.astype("float64")
    if np.any(~np.isfinite(length)) or _out_of_range(length, lo, hi):
        rep.errors.append(f"length NaN or outside [{lo}, {hi}] km (wrong unit, e.g. m?)")

    pol = ds.polarization.values.astype(str)
    bad_pol = sorted(set(pol[present].tolist()) - {"v", "h"})
    if bad_pol:
        rep.errors.append(f"polarization values {bad_pol} not in {{'v', 'h'}}")

    for name, key in [("site_0_lat", "lat_deg"), ("site_1_lat", "lat_deg"),
                      ("site_0_lon", "lon_deg"), ("site_1_lon", "lon_deg")]:
        lo, hi = rules[key]
        v = ds[name].values.astype("float64")
        if np.any(~np.isfinite(v)) or _out_of_range(v, lo, hi):
            rep.errors.append(f"{name} NaN or outside [{lo}, {hi}]")
    at_null_island = (ds.site_0_lat.values == 0) & (ds.site_0_lon.values == 0)
    if np.any(at_null_island):
        rep.errors.append(f"{int(at_null_island.sum())} links have site_0 at (0, 0): coordinates missing?")

    # ---- 5. data values -----------------------------------------------------------
    tsl = ds.tsl.values.astype("float64")
    rsl = ds.rsl.values.astype("float64")
    tl = ds.tl.values.astype("float64")
    both = np.isfinite(tsl) & np.isfinite(rsl)
    diff = np.abs(tl[both] - (tsl[both] - rsl[both]))
    if diff.size and np.max(diff) > rules["tl_tolerance_db"]:
        rep.errors.append(f"tl != tsl - rsl (max difference {np.max(diff):.3g} dB)")
    if np.any(np.isfinite(tl) & ~both):
        rep.errors.append("tl is finite where tsl or rsl is NaN")
    if np.any(np.isfinite(tsl[:, ~present]) | np.isfinite(rsl[:, ~present])):
        rep.errors.append("an absent sub-link (raw_sublink_id = -1) has data")

    for name, values, key in [("rsl", rsl, "rsl_plausible_dbm"), ("tsl", tsl, "tsl_plausible_dbm")]:
        lo, hi = rules[key]
        n_out = _out_of_range(values, lo, hi)
        if n_out:
            rep.warnings.append(f"{n_out} {name} values outside plausible range [{lo}, {hi}] dBm (kept)")

    if raise_on_error:
        rep.raise_if_errors()
    return rep


# ===========================================================================
# Reference data: gauges, weather, radar
# ===========================================================================
# Time convention of all L1 reference files (dataset attribute `time_convention`):
#   interval_start  value covers [time, time + interval_s)      (gauges)
#   instantaneous   value was sampled at `time`                   (radar, per the readme)
#   as_in_file      time kept exactly as in the raw file; each variable's meaning is in its
#                   `time_meaning_original` attribute            (SMHI weather: mixed meanings)
# Every file also carries interval_s, time_meaning_original and time_shift_applied_s.
TIME_CONVENTIONS = ("interval_start", "instantaneous", "as_in_file")
TIME_ATTRS = ("time_convention", "interval_s", "time_meaning_original", "time_shift_applied_s")

# ---- Gauge L1: one file per network, native interval -----------------------
GAUGE_DIMS = ("time", "id")
GAUGE_DATA_VARS = {
    # name -> (units, long_name, dtype)
    "rainfall_amount": ("mm", "rainfall amount accumulated over [time, time + interval_s)", "float32"),
}
GAUGE_QC_VAR = "qc_code"  # optional (SMHI only): data-provider quality code, uint8, 0 = no code
GAUGE_COORDS = {
    "id": (("id",), None, "station name as in the raw files"),
    "lon": (("id",), "degrees_east", "station longitude (WGS84)"),
    "lat": (("id",), "degrees_north", "station latitude (WGS84)"),
    "x": (("id",), "m", "station easting, EPSG:3006 (SWEREF 99 TM)"),
    "y": (("id",), "m", "station northing, EPSG:3006 (SWEREF 99 TM)"),
    "network": (("id",), None, "gauge network: city or smhi"),
    "gauge_type": (("id",), None, "gauge type: weighing or tipping_bucket"),
    "quantization_mm": (("id",), "mm", "measurement resolution"),
    "height_m": (("id",), "m", "gauge height above ground (NaN = unknown)"),
}
NETWORKS = ("city", "smhi")
GAUGE_TYPES = ("weighing", "tipping_bucket")

# ---- Weather L1: SMHI station, native hourly step ---------------------------
MET_DIMS = ("time", "id")
MET_DATA_VARS = {
    # name -> (units, long_name)
    "temperature": ("degC", "air temperature at 2 m"),
    "humidity": ("%", "relative humidity at 2 m"),
    "pressure": ("hPa", "air pressure reduced to sea level"),
    "wind_speed": ("m s-1", "wind speed at 10 m"),
    "wind_gust": ("m s-1", "gust wind speed at 10 m"),
    "wind_direction": ("degree", "wind direction at 10 m (from, degrees true)"),
}
MET_COORDS = {k: GAUGE_COORDS[k] for k in ("id", "lon", "lat", "x", "y")}

# ---- Radar L1: one file per month ------------------------------------------
RADAR_DIMS = ("time", "y", "x")
RADAR_DATA_VARS = {
    "dbz": ("dBZ", "radar reflectivity (decoded from pseudo-dBZ)", "float32"),
    "rain_rate": ("mm h-1", "rain rate from dbz with the file's Z-R relation Z = a R^b", "float32"),
}
RADAR_COORDS = {
    "x": (("x",), "m", "x of cell centre, radar projection"),
    "y": (("y",), "m", "y of cell centre, radar projection"),
    "lon": (("y", "x"), "degrees_east", "longitude of cell centre"),
    "lat": (("y", "x"), "degrees_north", "latitude of cell centre"),
    "x_3006": (("y", "x"), "m", "easting of cell centre, EPSG:3006"),
    "y_3006": (("y", "x"), "m", "northing of cell centre, EPSG:3006"),
}

# Used when a validator is called without rules (same values as reference_validation in openmrg.yaml).
REFERENCE_RULES = {
    "lat_deg": [50.0, 70.0],
    "lon_deg": [5.0, 25.0],
    "rain_rate_plausible_max_mm_h": 300.0,
    "dbz_plausible": [-35.0, 75.0],
    "zr_rel_tolerance": 1.0e-5,
}


def _check_structure(ds: xr.Dataset, rep: ValidationReport, dims: tuple, data_vars: dict, coords: dict,
                     interval, allowed_conventions: tuple) -> None:
    """Checks shared by the reference validators: variables, dims, units, dtypes, coords, time axis, time attrs."""
    for name, spec in data_vars.items():
        units, dtype = spec[0], (spec[2] if len(spec) > 2 else None)
        if name not in ds.data_vars:
            rep.errors.append(f"data variable '{name}' is missing")
            continue
        if ds[name].dims != dims:
            rep.errors.append(f"'{name}' has dims {ds[name].dims}, expected {dims}")
        if ds[name].attrs.get("units") != units:
            rep.errors.append(f"'{name}' units = {ds[name].attrs.get('units')!r}, expected {units!r}")
        if dtype is not None and ds[name].dtype != np.dtype(dtype):
            rep.errors.append(f"'{name}' has dtype {ds[name].dtype}, expected {dtype}")

    for name, (cdims, units, _) in coords.items():
        if name not in ds.coords:
            rep.errors.append(f"coordinate '{name}' is missing")
            continue
        if ds[name].dims != cdims:
            rep.errors.append(f"coordinate '{name}' has dims {ds[name].dims}, expected {cdims}")
        if units is not None and ds[name].attrs.get("units") != units:
            rep.errors.append(f"coordinate '{name}' units = {ds[name].attrs.get('units')!r}, expected {units!r}")
        if units is None and not ds[name].attrs.get("long_name"):
            rep.errors.append(f"coordinate '{name}' has no long_name")

    # Time axis: tz-naive datetime64 (tz-aware times end up as dtype object), unique, increasing, regular.
    if "time" not in ds.coords:
        rep.errors.append("coordinate 'time' is missing")
    elif ds.time.dtype.kind != "M":
        rep.errors.append(f"time has dtype {ds.time.dtype}; expected tz-naive datetime64 (UTC)")
    else:
        t = ds.time.values.astype("datetime64[ns]").astype("int64")
        steps = np.diff(t)
        if len(np.unique(t)) != len(t):
            rep.errors.append(f"time has {len(t) - len(np.unique(t))} duplicate values")
        if np.any(steps <= 0):
            rep.errors.append(f"time is not strictly increasing ({int(np.sum(steps <= 0))} bad steps)")
        if interval is not None and len(steps) > 0:
            n_irregular = int(np.sum(steps != pd.Timedelta(interval).value))
            if n_irregular:
                rep.errors.append(f"time is not regular at {interval}: {n_irregular} irregular steps")

    # Time convention attributes
    missing = [a for a in TIME_ATTRS if a not in ds.attrs]
    if missing:
        rep.errors.append(f"dataset attributes missing: {missing}")
    else:
        if ds.attrs["time_convention"] not in allowed_conventions:
            rep.errors.append(f"time_convention = {ds.attrs['time_convention']!r}, "
                              f"expected one of {allowed_conventions}")
        if interval is not None and int(ds.attrs["interval_s"]) != int(pd.Timedelta(interval).total_seconds()):
            rep.errors.append(f"interval_s = {ds.attrs['interval_s']}, "
                              f"expected {pd.Timedelta(interval).total_seconds():.0f}")


def _check_station_coords(ds: xr.Dataset, rep: ValidationReport, rules: dict) -> None:
    """lat/lon inside the plausible box, x/y finite, ids unique."""
    for name, key in (("lat", "lat_deg"), ("lon", "lon_deg")):
        lo, hi = rules[key]
        v = ds[name].values.astype("float64")
        if np.any(~np.isfinite(v)) or _out_of_range(v, lo, hi):
            rep.errors.append(f"{name} NaN or outside [{lo}, {hi}] degrees")
    for name in ("x", "y"):
        if np.any(~np.isfinite(ds[name].values.astype("float64"))):
            rep.errors.append(f"{name} (EPSG:3006) has NaN")
    if len(np.unique(ds["id"].values.astype(str))) != ds.sizes["id"]:
        rep.errors.append("station ids are not unique")


def validate_gauge_l1(ds: xr.Dataset, rules: dict | None = None, interval=None,
                      raise_on_error: bool = True) -> ValidationReport:
    """Check that `ds` follows the gauge L1 schema (GAUGE_* above).

    Errors: missing variables/coords, wrong dims/units/dtype, bad time axis or time attributes,
    time convention other than interval_start, unknown network or gauge type, bad qc_code.
    Warnings (values kept): negative amounts, amounts above the plausible rate, values that are
    not multiples of the gauge resolution.
    """
    rules = {**REFERENCE_RULES, **(rules or {})}
    rep = ValidationReport()
    _check_structure(ds, rep, GAUGE_DIMS, GAUGE_DATA_VARS, GAUGE_COORDS, interval, ("interval_start",))
    if rep.errors:
        if raise_on_error:
            rep.raise_if_errors("L1 gauge")
        return rep

    _check_station_coords(ds, rep, rules)
    bad = sorted(set(ds.network.values.astype(str)) - set(NETWORKS))
    if bad:
        rep.errors.append(f"network values {bad} not in {NETWORKS}")
    bad = sorted(set(ds.gauge_type.values.astype(str)) - set(GAUGE_TYPES))
    if bad:
        rep.errors.append(f"gauge_type values {bad} not in {GAUGE_TYPES}")
    q = ds.quantization_mm.values.astype("float64")
    if np.any(~np.isfinite(q) | (q <= 0)):
        rep.errors.append("quantization_mm must be finite and > 0")

    if GAUGE_QC_VAR in ds:
        qc = ds[GAUGE_QC_VAR]
        if qc.dims != GAUGE_DIMS or qc.dtype != np.uint8:
            rep.errors.append(f"{GAUGE_QC_VAR} must be uint8 on {GAUGE_DIMS}")
        elif "flag_values" not in qc.attrs or "flag_meanings" not in qc.attrs:
            rep.errors.append(f"{GAUGE_QC_VAR} lacks flag_values / flag_meanings")
        else:
            allowed = set(np.atleast_1d(qc.attrs["flag_values"]).tolist())
            unknown = sorted(set(np.unique(qc.values).tolist()) - allowed)
            if unknown:
                rep.errors.append(f"{GAUGE_QC_VAR} has values {unknown} not listed in flag_values")

    amount = ds.rainfall_amount.values.astype("float64")
    n_neg = int(np.sum(amount < 0))
    if n_neg:
        rep.warnings.append(f"{n_neg} negative rainfall amounts (kept)")
    limit = rules["rain_rate_plausible_max_mm_h"] * float(ds.attrs["interval_s"]) / 3600.0
    n_high = int(np.sum(amount > limit))
    if n_high:
        rep.warnings.append(f"{n_high} amounts above {limit:.3g} mm per interval "
                            f"({rules['rain_rate_plausible_max_mm_h']} mm/h) (kept)")
    # Amounts should be multiples of the gauge resolution (float rounding allowed).
    ratio = amount / q[None, :]
    off_grid = np.isfinite(ratio) & (np.abs(ratio - np.round(ratio)) > 1e-3)
    if np.any(off_grid):
        rep.warnings.append(f"{int(off_grid.sum())} amounts are not multiples of quantization_mm (kept)")

    if raise_on_error:
        rep.raise_if_errors("L1 gauge")
    return rep


def validate_met_l1(ds: xr.Dataset, rules: dict | None = None, interval=None,
                    raise_on_error: bool = True) -> ValidationReport:
    """Check that `ds` follows the weather L1 schema (MET_* above).

    Every variable must carry a `time_meaning_original` attribute, because the SMHI weather
    variables have different time meanings on one shared time axis.
    """
    rules = {**REFERENCE_RULES, **(rules or {})}
    rep = ValidationReport()
    _check_structure(ds, rep, MET_DIMS, MET_DATA_VARS, MET_COORDS, interval, TIME_CONVENTIONS)
    if rep.errors:
        if raise_on_error:
            rep.raise_if_errors("L1 weather")
        return rep
    _check_station_coords(ds, rep, rules)
    for name in MET_DATA_VARS:
        if not ds[name].attrs.get("time_meaning_original"):
            rep.errors.append(f"'{name}' has no time_meaning_original attribute")
    for name, lo, hi in [("humidity", 0, 100), ("wind_speed", 0, 100), ("wind_gust", 0, 100),
                         ("wind_direction", 0, 360), ("temperature", -60, 60), ("pressure", 850, 1100)]:
        n = _out_of_range(ds[name].values, lo, hi)
        if n:
            rep.warnings.append(f"{n} {name} values outside [{lo}, {hi}] (kept)")
    if raise_on_error:
        rep.raise_if_errors("L1 weather")
    return rep


def zr_rain_rate(dbz, a: float, b: float) -> np.ndarray:
    """Rain rate (mm/h) from reflectivity: Z = 10^(dBZ/10) in mm^6 m^-3, R = (Z / a)^(1/b).

    Computed in float64. NaN stays NaN. No threshold: -30 dBZ gives a tiny positive rate.
    """
    z = np.power(10.0, np.asarray(dbz, dtype="float64") / 10.0)
    return np.power(z / a, 1.0 / b)


def validate_radar_l1(ds: xr.Dataset, rules: dict | None = None, interval=None,
                      raise_on_error: bool = True) -> ValidationReport:
    """Check that `ds` follows the radar L1 schema (RADAR_* above).

    Besides the structure, rain_rate must be NaN exactly where dbz is NaN and must equal
    zr_rain_rate(dbz, zr_a, zr_b) with the a, b stored in the attributes.
    """
    rules = {**REFERENCE_RULES, **(rules or {})}
    rep = ValidationReport()
    _check_structure(ds, rep, RADAR_DIMS, RADAR_DATA_VARS, RADAR_COORDS, interval,
                     ("instantaneous", "interval_start"))
    for key in ("zr_a", "zr_b", "zr_source", "crs_proj4"):
        if key not in ds.attrs:
            rep.errors.append(f"dataset attribute '{key}' is missing")
    if rep.errors:
        if raise_on_error:
            rep.raise_if_errors("L1 radar")
        return rep

    dbz = ds.dbz.values
    rr = ds.rain_rate.values
    if not np.array_equal(np.isnan(dbz), np.isnan(rr)):
        rep.errors.append("rain_rate and dbz have different NaN patterns")
    else:
        ok = ~np.isnan(dbz)
        expected = zr_rain_rate(dbz[ok], float(ds.attrs["zr_a"]), float(ds.attrs["zr_b"]))
        rel = np.abs(rr[ok].astype("float64") - expected) / np.maximum(expected, 1e-30)
        if rel.size and rel.max() > rules["zr_rel_tolerance"]:
            rep.errors.append(f"rain_rate differs from Z-R of dbz (max relative difference {rel.max():.3g})")
    lo, hi = rules["dbz_plausible"]
    n = _out_of_range(dbz, lo, hi)
    if n:
        rep.warnings.append(f"{n} dbz values outside [{lo}, {hi}] (kept)")
    if raise_on_error:
        rep.raise_if_errors("L1 radar")
    return rep


# ===========================================================================
# L2 QC flags for the CML data
# ===========================================================================
# One uint16 bitmask per sample, dims (time, sublink_id, cml_id) = the L1 dims. The file documents
# itself with the CF attributes flag_masks / flag_meanings. Bits 8-15 are reserved: level_shift and
# pair_mismatch are deliberately NOT flags yet (not calibrated), so they must stay 0.
# The flag bits are defined once, in cmlpipe.qc.flags (QCFlag); re-exported here for the L2 validator.
from cmlpipe.qc.flags import FLAG_RESERVED_BITS, FLAGS, FLAGS_VAR, QCFlag, decode_flags, flag_attrs  # noqa: E402,F401


def validate_flags_l2(flags_ds: xr.Dataset, l1_ds: xr.Dataset, raise_on_error: bool = True) -> ValidationReport:
    """Check an L2 flags dataset against the L1 day it was computed from.

    Errors: flags variable missing; dims or sizes different from L1; dtype not uint16; bits set that
    are not defined in FLAGS; time axis or sub-link/link ids not equal to L1; flag_masks/flag_meanings
    missing or different from FLAGS.
    """
    rep = ValidationReport()
    if FLAGS_VAR not in flags_ds:
        rep.errors.append(f"variable '{FLAGS_VAR}' is missing")
    else:
        f = flags_ds[FLAGS_VAR]
        if f.dims != DIMS:
            rep.errors.append(f"'{FLAGS_VAR}' has dims {f.dims}, expected {DIMS}")
        elif any(flags_ds.sizes[d] != l1_ds.sizes[d] for d in DIMS):
            rep.errors.append(f"sizes {dict(flags_ds.sizes)} differ from L1 {dict(l1_ds.sizes)}")
        if f.dtype != np.uint16:
            rep.errors.append(f"'{FLAGS_VAR}' has dtype {f.dtype}, expected uint16")
        else:
            known = np.uint16(sum(1 << b for b in FLAGS.values()))
            unknown = int(np.count_nonzero(f.values & ~known))
            if unknown:
                rep.errors.append(f"{unknown} samples have bits set that are not defined in FLAGS")
        masks = f.attrs.get("flag_masks")
        if masks is None or list(np.atleast_1d(masks)) != [1 << b for b in FLAGS.values()] \
                or f.attrs.get("flag_meanings") != " ".join(FLAGS):
            rep.errors.append("flag_masks / flag_meanings missing or not equal to schema.FLAGS")
    for name in ("time", "sublink_id", "cml_id"):
        if name not in flags_ds.coords or flags_ds[name].shape != l1_ds[name].shape \
                or not np.array_equal(flags_ds[name].values, l1_ds[name].values):
            rep.errors.append(f"coordinate '{name}' is not equal to L1")
    if raise_on_error:
        rep.raise_if_errors("L2 flags")
    return rep


# ===========================================================================
# L3 products (aggregated): CML statistics, gauge sums, radar along links / at gauges
# ===========================================================================
L3_DIMS = {
    "cml": ("time", "sublink_id", "cml_id"),
    "gauges": ("time", "id"),
    "radar_links": ("time", "cml_id"),
    "radar_gauges": ("time", "id"),
}


def validate_l3(ds: xr.Dataset, kind: str, res: str, start=None, n_bins: int | None = None,
                raise_on_error: bool = True) -> ValidationReport:
    """Check an L3 product: dims, a regular START-labelled time axis at `res`, time attributes, units, counts.

    kind    : one of L3_DIMS ("cml", "gauges", "radar_links", "radar_gauges").
    start   : expected first label (a day for CML, the period start otherwise); n_bins: expected length.
    Errors: wrong dims; time not regular at res, not aligned to res, or not the expected start/length;
    time_convention / interval / interval_s / n_expected missing or wrong; a float variable without units;
    a count variable (n_valid*) larger than n_expected.
    """
    rep = ValidationReport()
    dims = L3_DIMS[kind]
    step = pd.Timedelta(res)
    for name, v in ds.data_vars.items():
        if v.dims != dims:
            rep.errors.append(f"'{name}' has dims {v.dims}, expected {dims}")
        if v.dtype.kind == "f" and not v.attrs.get("units"):
            rep.errors.append(f"float variable '{name}' has no units")
    t = pd.DatetimeIndex(ds.time.values)
    if len(t) == 0:
        rep.errors.append("empty time axis")
    else:
        if len(t) > 1 and not (np.diff(t.values) == step.to_timedelta64()).all():
            rep.errors.append(f"time axis is not regular at {res}")
        if not (t == t.floor(res)).all():
            rep.errors.append(f"time labels are not aligned to {res} bin starts")
        if start is not None and t[0] != pd.Timestamp(start):
            rep.errors.append(f"first label {t[0]} is not the expected {pd.Timestamp(start)}")
        if n_bins is not None and len(t) != n_bins:
            rep.errors.append(f"{len(t)} bins, expected {n_bins}")
    a = ds.attrs
    if a.get("time_convention") != "interval_start":
        rep.errors.append("attribute time_convention must be 'interval_start'")
    if a.get("interval") != res or int(a.get("interval_s", -1)) != int(step.total_seconds()):
        rep.errors.append(f"attributes interval/interval_s do not describe {res}")
    n_exp = a.get("n_expected")
    if n_exp is None:
        rep.errors.append("attribute n_expected missing")
    else:
        for name, v in ds.data_vars.items():
            if name.startswith("n_valid") and v.size and int(v.values.max()) > int(n_exp):
                rep.errors.append(f"'{name}' exceeds n_expected = {n_exp}")
    if raise_on_error:
        rep.raise_if_errors(f"L3 {kind} {res}")
    return rep
