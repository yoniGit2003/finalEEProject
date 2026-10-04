"""Release metadata: the ONE source for catalog.json, DATA_DICTIONARY.md and CMLData.describe().

PRODUCTS lists every table and time series of the release (name, kind, source, description, keys / dims, chunks);
column_description() describes every table column (explicit entries first, then name patterns). Time-series
variables are described by their own units / long_name attributes, copied from the work-folder NetCDF files.
"""

from __future__ import annotations

import re

RELEASE_NAME = "openmrg"
CITATION = ("Andersson, J. C. M., Olsson, J., van de Beek, R. (C. Z.), and Hansryd, J. (2022): OpenMRG: Open data from "
            "Microwave links, Radar, and Gauges for rainfall quantification in Gothenburg, Sweden, Earth Syst. Sci. Data, "
            "14, 5411-5426, https://doi.org/10.5194/essd-14-5411-2022. Data: https://doi.org/10.5281/zenodo.7107689 "
            "(CC BY-SA 4.0).")
DATA_LICENSE_SHORT = ("CC BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/), derived from OpenMRG "
                      "(Andersson et al. 2022, https://doi.org/10.5281/zenodo.7107689); see LICENSE.txt")


def data_license_text(name: str, version: str) -> str:
    """LICENSE.txt of a release (also the source of LICENSE-DATA.md in the repository)."""
    return f"""Data licence of the release {name}_{version}

This data release is distributed under the Creative Commons Attribution-ShareAlike 4.0 International licence
(CC BY-SA 4.0): https://creativecommons.org/licenses/by-sa/4.0/ (full legal code:
https://creativecommons.org/licenses/by-sa/4.0/legalcode).

It is a derived work (processed, quality-controlled and aggregated) of the OpenMRG dataset, which is itself
published under CC BY-SA 4.0. Any use must give attribution to the original authors of OpenMRG, and works built
on this release must be shared under the same licence.

Original data (cite this):
  {CITATION}

Processing: Haber, Y. and Tsitseluk, Y., cmlpipe (CellEnMon lab), {name} release {version}. Changes made to the
original data: format conversion, QC flags, aggregation, radar along link paths, events, dry periods, splits and
inventory tables (see DATA_GUIDE.md and DATA_DICTIONARY.md). The code that made it (cmlpipe) is MIT-licensed;
that licence does not apply to the data.
"""


CONVENTIONS = {
    "time": "UTC, tz-naive; every value is labelled by the START of its interval [t, t + dt) ('interval_start'); "
            "radar stamps are instants assigned to [t, t + 5 min)",
    "keys": "cml_id (int, link); sublink_id (str '<cml_id>_<1|2>'; sub-link 1 transmits site_0 -> site_1); gauge_id "
            "(City station name or 'GbgA'); event_id ('E' + start, e.g. E20150728_1455); run_id ('<event_id>_Rnn')",
    "units": "dB / dBm for signal levels, mm for amounts, mm/h for rates, m for distances (EPSG:3006), GHz",
    "missing": "NaN (float), fill value -32768 (int16, decoded to NaN), 0 for counts / flags",
    "splits": "split_events.csv and split_links.csv are frozen; train on (train events x train links)",
}

# ------------------------------------------------------------------ products
TABLES = {
    "links": ("tables/inventory/links.csv", ["cml_id"], "one row per link: radio parameters, split, status of both "
              "sub-links, completeness, TSL behaviour, nearest gauges, neighbours, radar cells"),
    "sublinks": ("tables/inventory/sublinks.csv", ["sublink_id"], "one row per sub-link: metadata, completeness (JJA "
                 "and per month), flag shares, quantization, TSL behaviour, dry-weather noise, quality status"),
    "gauges": ("tables/inventory/gauges.csv", ["gauge_id"], "City (10) and SMHI (1) gauges: position, type, "
               "completeness, JJA totals, radar pixel comparison"),
    "sites": ("tables/sites.csv", ["site_id"], "physical sites (masts): link ends merged within 10 m"),
    "events": ("tables/events.csv", ["event_id"], "rain events defined from radar + City gauges (never CML)"),
    "event_significance": ("tables/event_significance.csv", ["event_id"], "significance, rain class and wet core of "
                           "every event"),
    "event_runs": ("tables/event_runs.csv", ["run_id"], "wet runs inside events (before gap merging)"),
    "dry_periods": ("tables/dry_periods.csv", ["start"], "certified dry periods (>= 6 h)"),
    "split_events": ("tables/split_events.csv", ["event_id"], "FROZEN event split (train / val / test)"),
    "split_links": ("tables/split_links.csv", ["cml_id"], "FROZEN link split (train / test / excluded)"),
    "event_links": ("tables/event_links.csv", ["event_id", "cml_id"], "one row per (event, link): radar along the "
                    "path, matched gauge, CML data quality"),
    "event_links_flags_extra": ("tables/event_links_flags_extra.csv", ["event_id", "cml_id"], "per (event, link): "
                                "flatline / tsl_step shares and counts in the event (information only)"),
    "link_gauge_matches": ("tables/link_gauge_matches.csv", ["cml_id", "gauge_id"], "every (link, gauge) pair within "
                           "5 km of the path"),
    "network_gaps": ("tables/network_gaps.csv", ["start"], "intervals when many alive sub-links were missing "
                     "(collector outages)"),
    "quality_status": ("tables/inventory/quality_status.csv", ["sublink_id"], "quality status per sub-link "
                       "(good / degraded / unusable) and the reason"),
    "sublink_dry_stats": ("tables/sublink_dry_stats.csv", ["sublink_id"], "descriptive CML statistics per sub-link "
                          "over the certified dry periods"),
}
FROZEN_TABLES = ("events", "split_events", "split_links", "event_links", "dry_periods")

# name -> (description, resolution, dims, chunks)
TIMESERIES = {
    "cml_10s": ("CML at 10 s: tsl, rsl, tl (L1, int16 x 0.1 dB, exact copies) and the L2 QC flag bitmask", "10s",
                ("time", "sublink_id"), (60480, 64)),
    "cml_1min": ("CML statistics per 1 min bin (L3)", "1min", ("time", "sublink_id"), (132480, 32)),
    "cml_5min": ("CML statistics per 5 min bin (L3)", "5min", ("time", "sublink_id"), (8928, 64)),
    "cml_15min": ("CML statistics per 15 min bin (L3)", "15min", ("time", "sublink_id"), (8832, 64)),
    "cml_1h": ("CML statistics per 1 h bin (L3)", "1h", ("time", "sublink_id"), (2208, 728)),
    "radar_along_links_5min": ("radar rain rate along every link path (floored at 0 dBZ), 5 min", "5min",
                               ("time", "cml_id"), (8928, 64)),
    "radar_along_links_15min": ("radar rain rate along every link path, 15 min", "15min", ("time", "cml_id"), (8832, 64)),
    "radar_along_links_1h": ("radar rain rate along every link path, 1 h", "1h", ("time", "cml_id"), (2208, 364)),
    "radar_at_gauges_5min": ("radar rain rate at every gauge (pixel and 3 x 3 mean), 5 min", "5min", ("time", "id"), None),
    "radar_at_gauges_15min": ("radar rain rate at every gauge, 15 min", "15min", ("time", "id"), None),
    "radar_at_gauges_1h": ("radar rain rate at every gauge, 1 h", "1h", ("time", "id"), None),
    "gauges_city_1min": ("City gauges, 1 min amounts and rates", "1min", ("time", "id"), None),
    "gauges_city_5min": ("City gauges, 5 min sums (strict: NaN with a missing minute)", "5min", ("time", "id"), None),
    "gauges_city_15min": ("City gauges, 15 min sums", "15min", ("time", "id"), None),
    "gauges_city_1h": ("City gauges, 1 h sums", "1h", ("time", "id"), None),
    "gauges_smhi_15min": ("SMHI gauge GbgA, 15 min (native)", "15min", ("time", "id"), None),
    "gauges_smhi_1h": ("SMHI gauge GbgA, 1 h sums", "1h", ("time", "id"), None),
    "met_1h": ("SMHI weather at GbgA, hourly (instantaneous / as in the file; see each variable)", "1h", ("time", "id"), None),
    "network_reference_5min": ("network reference series: radar over the link domain + City gauges, state and event "
                               "id per 5 min bin", "5min", ("time",), None),
    "link_reference_state_5min": ("per-link reference state (wet / dry / unknown) and reasons, from references only",
                                  "5min", ("time", "cml_id"), (8928, 64)),
}

# ------------------------------------------------------------------ column descriptions
COLUMNS = {
    "cml_id": "link id (raw metadata 'Link')", "sublink_id": "sub-link id '<cml_id>_<1|2>'",
    "sublink": "sub-link number (1: site_0 -> site_1, 2: site_1 -> site_0)", "raw_sublink_id": "raw 'Sublink' number",
    "direction": "raw direction code A / B", "tx_site": "transmitting end (0 = site_0, 1 = site_1)",
    "rx_site": "receiving end (0 / 1)", "tx_site_id": "site id of the transmitter", "rx_site_id": "site id of the receiver",
    "gauge_id": "gauge id (City station name or 'GbgA')", "id": "gauge id", "event_id": "event id 'E<start>'",
    "run_id": "wet-run id '<event_id>_Rnn'", "site_id": "site (mast) id", "hop": "unordered pair of site ids; parallel "
    "links share it", "length_m": "path length (m, EPSG:3006)", "link_length_m": "path length (m)",
    "frequency_ghz": "carrier frequency (GHz)", "polarization": "'v' or 'h'", "length_class": "<1 km, 1-3 km, >3 km",
    "frequency_band": "<25 GHz or >=25 GHz (sub-link 1)", "stratum": "length class x frequency band",
    "excluded": "link excluded from the split (both sub-links unusable at step 3)",
    "split_link": "frozen link split: train / test / excluded", "split_event": "frozen event split: train / val / test",
    "site_0_id": "site id of end 0", "site_1_id": "site id of end 1", "n_cells": "radar cells crossed by the path",
    "n_radar_cells": "radar cells crossed by the path", "single_radar_cell": "the whole path lies in one radar cell",
    "lon": "longitude (deg, WGS84)", "lat": "latitude (deg, WGS84)", "x": "x (m, EPSG:3006)", "y": "y (m, EPSG:3006)",
    "x_3006_m": "x (m, EPSG:3006)", "y_3006_m": "y (m, EPSG:3006)", "n_links": "links with an end at the site",
    "n_ends": "link ends merged into the site", "max_member_offset_m": "largest offset of a merged end from the site (m)",
    "network": "'city' or 'smhi'", "gauge_type": "weighing or tipping_bucket", "quantization_mm": "gauge resolution (mm)",
    "height_m": "gauge height above ground (m)", "interval_s": "native interval (s)", "n_valid": "valid native steps",
    "completeness_pct": "valid share of native steps (%)", "total_mm": "L1 total over JJA (mm)",
    "n_negative": "negative amounts in the raw file", "jja_total_mm": "JJA total (mm)",
    "wet_fraction": "share of native steps with at least one quantum", "iy": "radar row index", "ix": "radar column index",
    "distance_to_cell_centre_m": "gauge distance to its radar cell centre (m)",
    "radar_pixel_jja_total_mm": "JJA radar total at the gauge pixel (mm)", "radar_over_gauge_ratio": "radar / gauge JJA total",
    "start": "start (UTC, interval start)", "end": "end (UTC, exclusive)", "duration_h": "duration (h)",
    "duration_min": "duration (min)", "gap_before_h": "non-wet time back to the previous event, capped at 48 h",
    "dry_before_h": "certified-dry hours in that gap", "gap_after_h": "non-wet time to the next event, capped at 48 h",
    "dry_after_h": "certified-dry hours in that gap", "merged_from": "wet runs merged into the event",
    "touches_period_edge": "less than one min_dry_gap from the period start / end",
    "radar_peak_domain_rain_mmh": "largest domain-mean radar rate (mm/h)", "radar_mean_domain_rain_mmh": "mean domain-mean "
    "radar rate (mm/h)", "radar_max_wet_area_fraction": "largest wet share of the domain cells",
    "radar_domain_total_mm": "domain-mean radar accumulation (mm averaged over the 720 domain cells)",
    "city_total_mean_mm": "mean of the City gauge totals (mm)", "city_total_max_mm": "largest City gauge total (mm)",
    "city_max_5min_rate_mmh": "largest 5 min City gauge rate (mm/h)", "city_max_5min_gauge": "gauge of that maximum",
    "city_max_1min_rate_mmh": "largest 1 min City gauge rate (mm/h)", "city_max_1min_gauge": "gauge of that maximum",
    "smhi_total_mm_info": "SMHI GbgA total over the overlapping 15 min bins (information only)",
    "n_unknown_bins": "5 min bins with state 'unknown'", "radar_missing_fraction": "share of bins without valid radar",
    "network_gap_minutes": "minutes of network gaps in the event", "n_weather_hours": "hourly weather values used",
    "month": "month of the start", "intensity_class": "tercile of the max 5 min City rate over events",
    "significant": "event (or run) passes the significance rule", "sig_reason": "which rule made it significant",
    "rain_class": "none_drizzle / significant / heavy", "n_runs": "wet runs in the event", "wet_core_start": "start of the "
    "first significant run", "wet_core_end": "end of the last significant run", "wet_core_hours": "wet core duration (h)",
    "gap_to_previous_run_h": "time since the previous run of the same event (h)", "significant_run": "run is significant",
    "radar_valid_fraction": "mean share of domain cells with radar data", "max_fraction_nan": "largest NaN share of alive "
    "sub-links", "n_alive": "alive sub-links that day",
    "distance_to_path_m": "gauge distance to the path segment (m)", "position_along_path": "0 = site_0 end, 1 = site_1 end",
    "distance_to_midpoint_m": "gauge distance to the path midpoint (m)", "distance_to_nearest_end_m": "to the nearest end (m)",
    "status": "quality status good / degraded / unusable (step 7 rules)", "status_reason": "rules that fired",
    "status_step3": "status of the step 3 QC (kept for reference)", "status_best": "better status of the two sub-links",
    "has_gauge_within_1km": "a gauge within 1 km of the path", "has_reference_gauge": "a gauge within 1 km of the path",
    "nearest_gauge_id": "nearest gauge at any distance", "nearest_gauge_network": "its network",
    "nearest_gauge_distance_m": "its distance to the path (m)", "radar_total_mm": "radar along the path, total in the event "
    "(mm)", "radar_peak_mmh": "largest 5 min R_mean in the event (mm/h)", "radar_wet_minutes": "minutes with R_mean >= "
    "0.1 mm/h", "radar_mean_wet_fraction": "mean wet share of the path", "radar_max_R_max_mmh": "largest R_max (mm/h)",
    "radar_mean_valid_weight_fraction": "mean share of the path with radar data", "gauge_network": "network of the matched "
    "gauge", "gauge_distance_to_path_m": "matched gauge distance to the path (m)", "gauge_position_along_path": "matched "
    "gauge position (0..1)", "gauge_total_mm": "matched gauge total in the event (mm)", "gauge_max_rate_mmh": "matched "
    "gauge max rate (mm/h)", "gauge_max_rate_interval": "interval of that rate (5min City, 15min SMHI)",
    "n_gauges_within_1km": "gauges within 1 km of the path", "cml_ok_event": "a sub-link has >= 80 % valid tl_mean bins in "
    "the event", "cml_ok_before": "... in the 3 h before", "cml_ok_after": "... in the 3 h after",
    "n_tsl_step_in_event": "10 s TSL steps in the event (both sub-links)", "within_500m": "gauge within 500 m",
    "within_1000m": "gauge within 1 km", "within_2000m": "gauge within 2 km",
    "rsl_min_step_db": "smallest RSL step (dB)", "rsl_n_distinct": "distinct RSL values", "tsl_n_distinct_values":
    "distinct TSL values", "tsl_top_share": "share of TSL at its most common value", "n_valid_tsl": "valid TSL samples",
    "tsl_n_steps": "TSL changes >= 1 dB", "tsl_step_fraction": "share of valid samples with a TSL change",
    "tsl_steps_per_day": "TSL changes per day of data", "tsl_behaviour": "no_tsl / toggling / constant / rare_steps / "
    "frequent_steps", "tsl_steps_in_events": "TSL steps inside events", "tsl_steps_in_dry": "TSL steps in dry periods",
    "tsl_step_rate_events_per_h": "TSL step rate inside events (per h)", "tsl_step_rate_dry_per_h": "rate in dry periods",
    "tsl_follows_fades_descriptive": "DESCRIPTIVE ONLY: TSL steps concentrate in rain", "share_not_blocked": "share of "
    "samples without a blocking flag", "flag_share_blocking_present": "share of samples present but blocking-flagged",
    "first_valid_day": "first day with valid tl", "last_valid_day": "last day with valid tl",
    "n_dry_days": "days with >= 1 h of valid dry tl", "n_dry_bins": "5 min bins in dry periods", "n_valid_tl":
    "valid tl_mean bins", "valid_fraction": "valid share of the dry bins", "dry_baseline_range_db": "p90 - p10 of the "
    "dry-period medians of tl (dB)", "n_dry_periods_used": "dry periods with >= 50 % valid bins",
    "dry_tl_10s_noise_median": "median 5 min tl_std in dry periods (10 s noise, dB)", "dry_rsl_10s_noise_median":
    "median 5 min rsl_std in dry periods (dB)", "dry_tl_median": "median dry tl (dB)", "dry_tl_std": "std of the dry 5 min "
    "tl means (dB)", "dry_tl_scaled_mad": "1.4826 x MAD of the dry tl means (dB)", "dry_rsl_std": "std of the dry rsl means",
    "dry_tl_diurnal_range": "max - min of the 24 hourly dry tl medians (dB)", "dry_tl_daily_drift_std": "std of the daily "
    "dry tl medians (dB)", "dry_tl_std_max": "larger dry tl noise of the two sub-links", "intra_bin": "",
    "n_neighbours_within_1km": "other links whose path is within 1 km", "n_neighbours_within_5km": "within 5 km",
    "completeness_tl_min": "worse sub-link tl completeness", "completeness_tl_max": "better sub-link tl completeness",
    "share_not_blocked_min": "worse sub-link not-blocked share", "frequency_ghz_sublink_1": "frequency of sub-link 1 (GHz)",
    "frequency_ghz_sublink_2": "frequency of sub-link 2 (GHz)", "site_0_antenna_diameter_m": "antenna diameter at site 0 (m)",
    "site_1_antenna_diameter_m": "antenna diameter at site 1 (m)",
}
PATTERNS = [
    (r"^city_total_(\w+)_mm$", "total of City gauge {0} in the event (mm)"),
    (r"^(temperature|humidity|wind_speed|wind_gust)_(min|max|mean)$", "{1} of the hourly SMHI {0} around the event"),
    (r"^completeness_(rsl|tsl|tl)$", "JJA share of samples with valid {0}"),
    (r"^completeness_tl_(jun|jul|aug)$", "share of samples with valid tl in {0}"),
    (r"^completeness_(jun|jul|aug)$", "share of valid native steps in {0}"),
    (r"^flag_share_(\w+)$", "JJA share of samples with the {0} flag"),
    (r"^valid_fraction_(event|before|after)_(sublink_[12]|min)$", "valid share of tl_mean bins {0} ({1})"),
    (r"^frac_(\w+?)_event_(sublink_[12]|max)$", "share of 5 min bins in the event whose flags_or holds {0} ({1})"),
    (r"^n_tsl_step_in_event_(sublink_[12])$", "10 s TSL steps in the event ({0})"),
    (r"^status_(sublink_[12])$", "quality status of {0}"), (r"^tsl_behaviour_(sublink_[12])$", "TSL behaviour of {0}"),
    (r"^n_gauges_within_(\d+)m$", "gauges within {0} m of the path"),
    (r"^nearest_(city|smhi)_gauge_id$", "nearest {0} gauge"), (r"^nearest_(city|smhi)_gauge_distance_m$", "its distance (m)"),
    (r"^(tl|rsl)_(median|mean|std|scaled_mad|iqr|p1|p99)$", "{1} of the dry 5 min {0} means (dB)"),
    (r"^(tl|rsl)_intra_bin_std_median$", "median 5 min {0}_std in dry periods (dB)"),
    (r"^tl_diurnal_range$", "max - min of the 24 hourly dry tl medians (dB)"),
    (r"^tl_daily_median_drift_std$", "std of the daily dry tl medians (dB)"),
]


def column_description(table: str, col: str) -> str:
    """Description of a table column (explicit entry, else a name pattern, else empty)."""
    if col in COLUMNS and COLUMNS[col]:
        return COLUMNS[col]
    for pat, text in PATTERNS:
        m = re.match(pat, col)
        if m:
            return text.format(*m.groups())
    return ""


def column_unit(col: str) -> str:
    """Unit of a table column, from its name (ambiguous suffixes are listed explicitly first)."""
    explicit = {"duration_min": "min", "network_gap_minutes": "min", "radar_wet_minutes": "min",
                "tsl_step_rate_events_per_h": "1/h", "tsl_step_rate_dry_per_h": "1/h", "wet_core_hours": "h",
                "temperature_min": "degC", "temperature_max": "degC", "temperature_mean": "degC", "humidity_min": "%",
                "humidity_max": "%", "humidity_mean": "%", "lon": "degrees_east", "lat": "degrees_north",
                "x": "m", "y": "m", "tsl_steps_per_day": "1/day"}
    if col in explicit:
        return explicit[col]
    if col.startswith("wind_"):
        return "m s-1"
    for suf, u in (("_mmh", "mm/h"), ("_mm", "mm"), ("_m", "m"), ("_db", "dB"), ("_ghz", "GHz"), ("_h", "h"),
                   ("_pct", "%"), ("_km2", "km2")):
        if col.endswith(suf):
            return u
    if col.startswith(("dry_tl", "dry_rsl", "tl_", "rsl_")) and not col.endswith(("_days", "_used", "_bins")):
        return "dB"
    return ""


# Descriptions of time-series variables whose source files carry units but no long_name.
VARIABLES = {
    "radar_wet_area_fraction": "share of the valid domain cells with floored radar rate >= 0.1 mm/h",
    "radar_valid_cell_fraction": "share of the 720 domain cells with radar data",
    "radar_n_wet_cells": "number of wet domain cells",
    "domain_mean_rain_mmh": "mean floored radar rate over the domain cells (mm/h)",
    "max_city_rate_mmh": "largest 5 min City gauge rate (mm/h)",
    "n_city_gauges_valid": "City gauges with a valid 5 min amount",
    "n_city_gauges_wet": "City gauges with >= 0.1 mm in the 5 min bin",
}


def variable_description(name: str, attrs: dict) -> str:
    """Description of a time-series variable: its long_name, else VARIABLES."""
    return attrs.get("long_name") or VARIABLES.get(name, "")
