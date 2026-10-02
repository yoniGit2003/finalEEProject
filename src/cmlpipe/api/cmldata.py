"""CMLData: read a cmlpipe release (tables + Zarr time series) with a small, documented API.

    from cmlpipe.api import CMLData
    d = CMLData.open("data/release/openmrg_v1")
    links = d.links(length_km=(1, 3), frequency_band=">=25GHz", status="good", split="train", has_gauge_within_m=1000)
    ev = d.events(significant=True, split="train")
    ds = d.cml(links.cml_id[:5], "2015-07-28", "2015-07-29", "5min")           # xarray, dims (time, sublink_id)
    df = d.frame(10227, "2015-07-28 12:00", "2015-07-28 20:00", res="5min")    # one aligned wide DataFrame
    w = d.event_window("E20150727_2140", cml_id=10227)                          # frame + event / run / dry markers
    d.describe("tl_mean")

Conventions: times are UTC, tz-naive, and a value at t covers [t, t + res) (interval start). `start` / `end` accept
strings or datetimes (tz-aware ones are converted to UTC) and select [start, end). Ids may be cml_ids (int) or
sublink_ids ("<cml_id>_<1|2>"). Tables come back as pandas DataFrames, time series as xarray Datasets (lazy Zarr,
only the requested chunks are read; no dask). Stores are opened once per CMLData object.

QC of the 10 s data (cml(..., res="10s", qc=...)):
  "mask"  rsl / tsl / tl set to NaN where the L2 flags hold a blocking flag (rsl: the rsl set, tsl: the tsl set,
          tl: their union; the same sets the L3 statistics use; listed in catalog["qc"])
  "raw"   the values as recorded (no masking, no flags)
  "flags" the raw values plus the uint16 `flags` bitmask
At 1 min and coarser the statistics already exclude blocked samples (the L3 rule); "mask" and "raw" return them as
stored and "flags" adds `flags_or` (bitwise OR of the flags in the bin).
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from cmlpipe.qc.flags import mask_of
from cmlpipe.release import zarrstore

RES_ORDER = ["10s", "1min", "5min", "15min", "1h"]
STATE_NAMES = {0: "dry", 1: "wet", -1: "unknown"}


class CMLDataError(ValueError):
    """Bad request: unknown id or resolution, empty or too large window."""


def _ts(t) -> pd.Timestamp:
    ts = pd.Timestamp(t)
    return ts.tz_convert("UTC").tz_localize(None) if ts.tzinfo is not None else ts


class CMLData:
    """Access to one release folder (see module docstring)."""

    def __init__(self, path: Path, catalog: dict):
        self.path = Path(path)
        self.catalog = catalog
        self._stores: dict[str, xr.Dataset] = {}

    # ------------------------------------------------------------------ opening
    @classmethod
    def open(cls, release_path) -> "CMLData":
        """Open a release folder (reads catalog.json only; tables and stores are opened on first use)."""
        p = Path(release_path)
        cat = p / "catalog.json"
        if not cat.exists():
            raise FileNotFoundError(f"{cat} not found: not a cmlpipe release folder")
        return cls(p, json.loads(cat.read_text(encoding="utf-8")))

    def __repr__(self) -> str:
        c = self.catalog
        return f"CMLData({c['name']} {c['version']}, {len(c['tables'])} tables, {len(c['timeseries'])} time series)"

    @lru_cache(maxsize=None)  # noqa: B019 - one object per release; caching per object is intended
    def table(self, name: str) -> pd.DataFrame:
        """Any release table by name (catalog['tables']), as a DataFrame."""
        if name not in self.catalog["tables"]:
            raise CMLDataError(f"unknown table {name!r}; tables: {sorted(self.catalog['tables'])}")
        return pd.read_parquet(self.path / self.catalog["tables"][name]["paths"]["parquet"])

    def store(self, name: str) -> xr.Dataset:
        """A time-series store (lazy, opened once)."""
        if name not in self.catalog["timeseries"]:
            raise CMLDataError(f"unknown time series {name!r}; available: {sorted(self.catalog['timeseries'])}")
        if name not in self._stores:
            self._stores[name] = zarrstore.open_store(self.path / self.catalog["timeseries"][name]["path"])
        return self._stores[name]

    # ------------------------------------------------------------------ tables
    @staticmethod
    def _filter(df: pd.DataFrame, filters: dict, special: dict) -> pd.DataFrame:
        m = np.ones(len(df), bool)
        for k, v in filters.items():
            if v is None:
                continue
            if k in special:
                m &= np.asarray(special[k](df, v), bool)
            elif k in df.columns:
                m &= df[k].isin(v).to_numpy() if isinstance(v, (list, tuple, set, np.ndarray, pd.Series)) else (df[k] == v).to_numpy()
            else:
                raise CMLDataError(f"unknown filter {k!r}; use a column name or one of {sorted(special)}")
        return df[m].reset_index(drop=True)

    def links(self, **filters) -> pd.DataFrame:
        """Links (one row per cml_id, columns: see describe / DATA_DICTIONARY 'links').

        Filters (all optional, combined with AND): length_km=(min, max); frequency_band="<25GHz" | ">=25GHz";
        status="good" | "degraded" | "unusable" (the better of the two sub-links, status_best); split="train" |
        "test" | "excluded"; has_gauge_within_m=250 | 500 | 1000; tsl_behaviour=<class of either sub-link>;
        cml_ids=[...]; or any column name = value (or list of values).
        """
        def gauge(df, v):
            col = f"n_gauges_within_{int(v)}m"
            if col not in df:
                raise CMLDataError(f"has_gauge_within_m must be one of 250, 500, 1000 (got {v})")
            return df[col] > 0
        special = {"length_km": lambda df, v: df.length_m.between(1000 * v[0], 1000 * v[1]),
                   "status": lambda df, v: df.status_best == v, "split": lambda df, v: df.split_link == v,
                   "has_gauge_within_m": gauge, "cml_ids": lambda df, v: df.cml_id.isin(list(v)),
                   "tsl_behaviour": lambda df, v: (df.tsl_behaviour_sublink_1 == v) | (df.tsl_behaviour_sublink_2 == v)}
        return self._filter(self.table("links"), filters, special)

    def sublinks(self, **filters) -> pd.DataFrame:
        """Sub-links (one row per sublink_id). Filters: length_km, frequency_band, status (of THIS sub-link), split,
        tsl_behaviour, cml_ids, or any column = value."""
        special = {"length_km": lambda df, v: df.length_m.between(1000 * v[0], 1000 * v[1]),
                   "split": lambda df, v: df.split_link == v, "cml_ids": lambda df, v: df.cml_id.isin(list(v))}
        return self._filter(self.table("sublinks"), filters, special)

    def gauges(self) -> pd.DataFrame:
        """Gauges (gauge_id, network, position, completeness, JJA totals, radar pixel comparison)."""
        return self.table("gauges")

    def events(self, significant: bool | None = None, split: str | None = None, month: int | None = None,
               rain_class: str | None = None) -> pd.DataFrame:
        """Events joined with their significance (columns of 'events' plus significant, rain_class, n_runs, wet core)."""
        sg = self.table("event_significance")[["event_id", "significant", "sig_reason", "rain_class", "n_runs",
                                               "wet_core_start", "wet_core_end", "wet_core_hours"]]
        df = self.table("events").merge(sg, on="event_id")
        return self._filter(df, {"significant": significant, "split_event": split, "month": month,
                                 "rain_class": rain_class}, {})

    def event_pairs(self, significant: bool | None = True, split=None, usable_only: bool = True,
                    max_flatline_share: float | None = 0.5, min_valid_fraction: float | None = 0.8,
                    return_dropped: bool = False):
        """Event x link pairs (event_links + event_links_flags_extra + event_significance) for modelling.

        Selection (rows outside it are not returned at all): significant=True | False | None (all events);
        split="train" | "val" | "test" (the event split) or a tuple (event_split, link_split), e.g. ("test", "train")
        for unseen events on training links; None in either place = no selection.
        Quality filters (None or False disables one); a pair that fails any of them is dropped:
          unusable   usable_only and no sub-link has status good / degraded
          flatline   frac_flatline_event_max > max_flatline_share (rsl frozen for most of the event)
          low_valid  the best valid share of tl_mean bins in the event over the usable sub-links < min_valid_fraction
        Returns the kept pairs; with return_dropped=True returns (kept, dropped), where dropped has a `drop_reason`
        column listing every failed filter ("flatline; low_valid").
        """
        fx = self.table("event_links_flags_extra").drop(columns=["split_event", "split_link"])
        sg = self.table("event_significance")[["event_id", "significant", "sig_reason", "rain_class",
                                               "wet_core_start", "wet_core_end", "wet_core_hours"]]
        df = self.table("event_links").merge(fx, on=["event_id", "cml_id"], how="left").merge(sg, on="event_id")

        # selection
        ev_split, link_split = split if isinstance(split, tuple) else (split, None)
        df = self._filter(df, {"significant": significant, "split_event": ev_split, "split_link": link_split}, {})

        # quality filters: one boolean "fails" column per filter
        usable = {s: df[f"status_sublink_{s}"].isin(["good", "degraded"]) for s in (1, 2)}
        fails = {}
        if usable_only:
            fails["unusable"] = ~(usable[1] | usable[2])
        if max_flatline_share is not None:
            fails["flatline"] = df.frac_flatline_event_max > max_flatline_share
        if min_valid_fraction is not None:
            # valid share of the better usable sub-link (an unusable sub-link does not count)
            vf = pd.concat([df[f"valid_fraction_event_sublink_{s}"].where(usable[s] | (not usable_only))
                            for s in (1, 2)], axis=1).max(axis=1)
            fails["low_valid"] = ~(vf >= min_valid_fraction)  # NaN (no usable sub-link) also fails
        failed = {name: f.to_numpy(bool) for name, f in fails.items()}
        reason = pd.Series(["; ".join(n for n, f in failed.items() if f[i]) for i in range(len(df))],
                           index=df.index, dtype=object)
        drop = (reason != "").to_numpy()
        kept = df[~drop].reset_index(drop=True)
        if not return_dropped:
            return kept
        return kept, df[drop].assign(drop_reason=reason[drop]).reset_index(drop=True)

    def event_runs(self, event_id: str) -> pd.DataFrame:
        """The wet runs of one event (run_id, start, end, totals, significance)."""
        self._check_event(event_id)
        r = self.table("event_runs")
        return r[r.event_id == event_id].reset_index(drop=True)

    def dry_periods(self) -> pd.DataFrame:
        """Certified dry periods (start, end, duration_h, radar_valid_fraction)."""
        return self.table("dry_periods")

    def status(self) -> pd.DataFrame:
        """Quality status per sub-link (status, reason and the metrics behind it)."""
        return self.table("quality_status")

    # ------------------------------------------------------------------ helpers
    def _check_event(self, event_id: str) -> pd.Series:
        ev = self.table("events").set_index("event_id")
        if event_id not in ev.index:
            raise CMLDataError(f"unknown event_id {event_id!r}")
        return ev.loc[event_id]

    def _sublink_ids(self, ids) -> list[str]:
        known = self.table("sublinks").sublink_id
        out = []
        for i in np.atleast_1d(ids if not isinstance(ids, (str, int, np.integer)) else [ids]):
            s = str(i)
            if "_" in s:
                if s not in set(known):
                    raise CMLDataError(f"unknown sublink_id {s!r}")
                out.append(s)
            else:
                subs = known[known.str.startswith(f"{int(s)}_")].tolist()
                if not subs:
                    raise CMLDataError(f"unknown cml_id {s}")
                out += subs
        return out

    def _cml_ids(self, ids) -> list[int]:
        known = set(self.table("links").cml_id)
        out = []
        for i in np.atleast_1d(ids if not isinstance(ids, (str, int, np.integer)) else [ids]):
            c = int(str(i).split("_")[0])
            if c not in known:
                raise CMLDataError(f"unknown cml_id {c}")
            if c not in out:
                out.append(c)
        return out

    def _window(self, start, end) -> tuple[pd.Timestamp, pd.Timestamp]:
        t0, t1 = _ts(start), _ts(end)
        if t1 <= t0:
            raise CMLDataError(f"empty window [{t0}, {t1})")
        p = self.catalog["period"]
        if t0 < pd.Timestamp(p["start"]) or t1 > pd.Timestamp(p["end_exclusive"]):
            raise CMLDataError(f"window [{t0}, {t1}) is outside the period [{p['start']}, {p['end_exclusive']})")
        return t0, t1

    @staticmethod
    def _slice(ds: xr.Dataset, t0, t1) -> xr.Dataset:
        return ds.sel(time=slice(t0, t1 - pd.Timedelta("1ns")))

    # ------------------------------------------------------------------ time series
    def cml(self, ids, start, end, res: str = "5min", variables=None, qc: str = "mask", allow_large: bool = False) -> xr.Dataset:
        """CML time series, dims (time, sublink_id), with sub-link metadata as coordinates.

        res: "10s" | "1min" | "5min" | "15min" | "1h". variables: list (default all). qc: see module docstring.
        10 s guard: more than 7 days, or more than 50 sub-links for more than 2 days, raises unless allow_large=True.
        """
        if res not in RES_ORDER:
            raise CMLDataError(f"unknown resolution {res!r}; available: {RES_ORDER}")
        if qc not in ("mask", "raw", "flags"):
            raise CMLDataError("qc must be 'mask', 'raw' or 'flags'")
        subs = self._sublink_ids(ids)
        t0, t1 = self._window(start, end)
        if res == "10s" and not allow_large:
            days = (t1 - t0) / pd.Timedelta("1D")
            if days > 7 or (len(subs) > 50 and days > 2):
                raise CMLDataError(f"10 s request of {len(subs)} sub-links x {days:.1f} days is large; "
                                   "pass allow_large=True (or use a coarser resolution)")
        ds = self._slice(self.store(f"cml_{res}"), t0, t1).sel(sublink_id=subs)
        flag_var = "flags" if res == "10s" else "flags_or"
        names = list(variables) if variables is not None else [v for v in ds.data_vars if v != flag_var]
        unknown = [v for v in names if v not in ds.data_vars]
        if unknown:
            raise CMLDataError(f"unknown variables {unknown}; available: {list(ds.data_vars)}")
        need = names + ([flag_var] if (qc == "flags" or (res == "10s" and qc == "mask")) and flag_var not in names else [])
        out = ds[need].load()
        if res == "10s" and qc == "mask":
            q = self.catalog["qc"]
            fl = out["flags"].values
            sets = {"rsl": mask_of(q["blocking_flags_rsl"]), "tsl": mask_of(q["blocking_flags_tsl"])}
            sets["tl"] = sets["rsl"] | sets["tsl"]
            for v in ("rsl", "tsl", "tl"):
                if v in out:
                    out[v] = out[v].where((fl & sets[v]) == 0)
            if "flags" not in names:
                out = out.drop_vars("flags")
            out.attrs["qc"] = "mask: samples with a blocking flag set to NaN (rsl / tsl sets; tl = union)"
        else:
            out.attrs["qc"] = qc
        return out

    def radar_along(self, cml_ids, start, end, res: str = "5min") -> xr.Dataset:
        """Radar along the link paths (R_mean, R_std, R_max, wet_fraction, valid_weight_fraction, n_valid_stamps),
        dims (time, cml_id). Available at 5min / 15min / 1h only (5 min native, never interpolated)."""
        name = f"radar_along_links_{res}"
        if name not in self.catalog["timeseries"]:
            raise CMLDataError(f"radar along links is not available at {res} (5min, 15min, 1h; never interpolated)")
        t0, t1 = self._window(start, end)
        return self._slice(self.store(name), t0, t1).sel(cml_id=self._cml_ids(cml_ids)).load()

    def gauge(self, gauge_ids, start, end, res: str = "5min") -> xr.Dataset:
        """Gauge amounts / rates, dims (time, gauge_id). City: 1min..1h; SMHI (GbgA): 15min, 1h (never disaggregated)."""
        t0, t1 = self._window(start, end)
        ids = [str(g) for g in np.atleast_1d(gauge_ids)]
        known = self.gauges().set_index("gauge_id").network
        parts = []
        for net in ("city", "smhi"):
            sel = [g for g in ids if g in known.index and known[g] == net]
            if not sel:
                continue
            name = f"gauges_{net}_{res}"
            if name not in self.catalog["timeseries"]:
                raise CMLDataError(f"{net} gauges are not available at {res}")
            ds = self._slice(self.store(name), t0, t1).sel(gauge_id=sel).load()
            parts.append(ds.drop_vars([c for c in ds.coords if c not in ("time", "gauge_id")]))
        bad = [g for g in ids if g not in known.index]
        if bad:
            raise CMLDataError(f"unknown gauge_id {bad}")
        return xr.concat(parts, "gauge_id", join="outer") if len(parts) > 1 else parts[0]

    def met(self, start, end) -> xr.Dataset:
        """SMHI hourly weather at GbgA (stamps are the instants / interval ends given in each variable's attrs)."""
        t0, t1 = self._window(start, end)
        return self._slice(self.store("met_1h"), t0, t1).load()

    def link_state(self, cml_ids, start, end) -> xr.Dataset:
        """Per-link reference state (5 min): state (0 dry, 1 wet, -1 unknown) and reasons, dims (time, cml_id)."""
        t0, t1 = self._window(start, end)
        return self._slice(self.store("link_reference_state_5min"), t0, t1).sel(cml_id=self._cml_ids(cml_ids)).load()

    def network_state(self, start, end) -> xr.Dataset:
        """Network reference series (5 min): radar wet area, gauges, state, event_id."""
        t0, t1 = self._window(start, end)
        return self._slice(self.store("network_reference_5min"), t0, t1).load()

    # ------------------------------------------------------------------ wide frames
    @staticmethod
    def _to_res(s: pd.Series, index: pd.DatetimeIndex, how: str) -> pd.Series:
        """A 5 min state series on another grid: finer -> forward-fill within the 5 min bin; coarser -> 'state' rule
        (wet if any bin wet, dry if all dry, else unknown) or the max."""
        step = index[1] - index[0] if len(index) > 1 else pd.Timedelta("5min")
        if step <= pd.Timedelta("5min"):
            return s.reindex(index.floor("5min")).set_axis(index)
        g = s.groupby(s.index.floor(step))
        if how == "state":
            key = s.index.floor(step)
            any_wet = (s == 1).groupby(key).any()
            all_dry = (s == 0).groupby(key).all()
            agg = pd.Series(np.where(any_wet, 1, np.where(all_dry, 0, -1)), index=any_wet.index)
        else:
            agg = g.max()
        return agg.reindex(index)

    def frame(self, cml_id, start, end, res: str = "5min", include=("cml", "flags", "radar", "gauge", "met", "state"),
              max_gauge_dist_m: float = 1000, qc: str = "mask") -> pd.DataFrame:
        """ONE aligned wide DataFrame for one link on the `res` grid of [start, end).

        Columns: cml statistics per sub-link with suffix _s1 / _s2 (at 10 s: tsl, rsl, tl); flags_or_s1/_s2 (or
        flags_s1/_s2 at 10 s); radar_* (along the path; only at 5min / 15min / 1h); gauge_* of the nearest gauge within
        max_gauge_dist_m of the path (from link_gauge_matches; none -> no columns); met_* hourly SMHI values
        forward-filled within the hour, with met_filled = True where the value is carried from the hour stamp;
        link_state / network_state (5 min states mapped to res; see _to_res). Notes on left-out sources are in
        df.attrs["notes"], the chosen gauge in df.attrs["gauge"].
        """
        c = self._cml_ids(cml_id)[0]
        t0, t1 = self._window(start, end)
        index = pd.date_range(t0, t1 - pd.Timedelta(res), freq=res)
        df = pd.DataFrame(index=index)
        df.index.name = "time"
        notes = []
        if "cml" in include:
            ds = self.cml(c, t0, t1, res, qc=qc)  # values with the chosen QC (default: masked at 10 s)
            for sid in ds.sublink_id.values:
                suf = "_s" + str(sid).split("_")[1]
                for v in ds.data_vars:
                    if v not in ("flags", "flags_or"):
                        df[v + suf] = ds[v].sel(sublink_id=sid).to_series().reindex(index)
        if "flags" in include:
            fv = "flags" if res == "10s" else "flags_or"
            fl = self._slice(self.store(f"cml_{res}"), t0, t1).sel(sublink_id=self._sublink_ids(c))[fv].load()
            for sid in fl.sublink_id.values:
                df[f"{fv}_s" + str(sid).split("_")[1]] = fl.sel(sublink_id=sid).to_series().reindex(index)
        if "radar" in include:
            if f"radar_along_links_{res}" in self.catalog["timeseries"]:
                r = self.radar_along(c, t0, t1, res).sel(cml_id=c)
                for v in r.data_vars:
                    df["radar_" + v] = r[v].to_series().reindex(index)
            else:
                notes.append(f"radar: not available at {res} (5 min native, never interpolated)")
        if "gauge" in include:
            m = self.table("link_gauge_matches")
            m = m[(m.cml_id == c) & (m.distance_to_path_m <= max_gauge_dist_m)].sort_values("distance_to_path_m")
            chosen = None
            for g in m.itertuples():
                if f"gauges_{g.network}_{res}" in self.catalog["timeseries"]:
                    chosen = g
                    break
                notes.append(f"gauge {g.gauge_id} ({g.network}) not available at {res}")
            if chosen is not None:
                gd = self.gauge(chosen.gauge_id, t0, t1, res).sel(gauge_id=chosen.gauge_id)
                for v in gd.data_vars:
                    df["gauge_" + v] = gd[v].to_series().reindex(index)
                df.attrs["gauge"] = {"gauge_id": chosen.gauge_id, "network": chosen.network,
                                     "distance_to_path_m": float(chosen.distance_to_path_m),
                                     "position_along_path": float(chosen.position_along_path)}
            else:
                notes.append(f"no gauge within {max_gauge_dist_m:g} m available at {res}")
        if "met" in include:
            mt = self.met(t0.floor("1h"), min(t1.ceil("1h") + pd.Timedelta("1h"), pd.Timestamp(self.catalog["period"]["end_exclusive"])))
            hour = index.floor("1h")
            for v in mt.data_vars:
                if not v.startswith("qc_"):
                    s = mt[v].isel(gauge_id=0).to_series() if "gauge_id" in mt[v].dims else mt[v].to_series()
                    df["met_" + v] = s.reindex(hour).to_numpy()
            df["met_filled"] = index != hour  # True: carried forward from the hourly stamp within the hour
        if "state" in include:
            ls = self.link_state(c, t0.floor("5min"), t1.ceil("5min")).state.sel(cml_id=c).to_series()
            ns = self.network_state(t0.floor("5min"), t1.ceil("5min")).state.to_series()
            # network reference codes: 0 dry, 1 wet, 2 unknown -> same codes as the link state (-1 unknown)
            ns = ns.replace({2: -1})
            df["link_state"] = self._to_res(ls, index, "state").astype("Int8")
            df["network_state"] = self._to_res(ns, index, "state").astype("Int8")
        df.attrs.update({"cml_id": c, "res": res, "time_convention": "interval_start (UTC)", "notes": notes})
        return df

    def event_window(self, event_id: str, cml_id=None, pad: str = "12h", res: str = "5min", **frame_kw) -> pd.DataFrame:
        """frame() of [event start - pad, event end + pad) with markers: in_event, run_id ('' outside runs),
        in_wet_core, in_dry_period. Without cml_id: the network state with the same markers."""
        ev = self._check_event(event_id)
        p = self.catalog["period"]
        t0 = max(pd.Timestamp(ev.start) - pd.Timedelta(pad), pd.Timestamp(p["start"]))
        t1 = min(pd.Timestamp(ev.end) + pd.Timedelta(pad), pd.Timestamp(p["end_exclusive"]))
        if cml_id is not None:
            df = self.frame(cml_id, t0, t1, res, **frame_kw)
        else:
            ns = self.network_state(t0.floor("5min"), t1.ceil("5min")).to_dataframe()
            index = pd.date_range(t0, t1 - pd.Timedelta(res), freq=res)
            df = ns.reindex(index.floor("5min")).set_axis(index) if pd.Timedelta(res) <= pd.Timedelta("5min") \
                else ns.groupby(ns.index.floor(res)).first().reindex(index)
            df.index.name = "time"
        idx = df.index
        df["in_event"] = (idx >= pd.Timestamp(ev.start)) & (idx < pd.Timestamp(ev.end))
        runs = self.event_runs(event_id)
        rid = np.full(len(idx), "", dtype=object)
        for r in runs.itertuples():
            rid[(idx >= r.start) & (idx < r.end)] = r.run_id
        df["run_id"] = rid
        sg = self.table("event_significance").set_index("event_id").loc[event_id]
        df["in_wet_core"] = ((idx >= sg.wet_core_start) & (idx < sg.wet_core_end)) if pd.notna(sg.wet_core_start) else False
        dp = self.dry_periods()
        dry = np.zeros(len(idx), bool)
        for r in dp.itertuples():
            dry |= (idx >= r.start) & (idx < r.end)
        df["in_dry_period"] = dry
        df.attrs.update({"event_id": event_id, "pad": pad})
        return df

    # ------------------------------------------------------------------ dictionary
    def describe(self, variable: str) -> pd.DataFrame:
        """Where `variable` appears (tables and time series) with its unit and description (from the catalog)."""
        rows = []
        for t, e in self.catalog["tables"].items():
            for c in e["columns"]:
                if c["name"] == variable:
                    rows.append({"product": t, "kind": "table", "dtype": c["dtype"], "unit": c["unit"],
                                 "description": c["description"]})
        for n, e in self.catalog["timeseries"].items():
            if variable in e["variables"]:
                v = e["variables"][variable]
                rows.append({"product": n, "kind": "timeseries", "dtype": v["dtype"], "unit": v["units"],
                             "description": v["description"]})
        if not rows:
            raise CMLDataError(f"{variable!r} is not a column or variable of this release")
        return pd.DataFrame(rows)
