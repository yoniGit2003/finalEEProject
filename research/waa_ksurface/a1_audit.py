"""STATUS (Detachment correction, 2026-11-11): part (b), the reproduction of the Part 1 length test, stands as a one-off
data-loading sanity check, diagnostic only. Nothing later builds on it.

Step A1 of the k-surface branch: data audit and reproduction gate (TRAIN only).

    python -m research.waa_ksurface.a1_audit

(a) rebuilds the Part 1 TRAIN dataset with the Part 1 code and checks the counts (206 links, 9 780 link-events);
(b) reruns the Part 1 length test on it and compares every intercept I and slope s with the stored Part 1 table;
(c) counts the TRAIN data (bands, length classes, rain classes on short links, the gauge-input set, g_e);
(d) locates the stored comparison fits (no refit);
(e) applies the 32 GHz rule (a band without an eligible length-test cell is excluded from the branch).
A failed gate in (a) or (b) stops the script with exit code 1 (the rule: STOP, no patch).
"""
import json
import sys

import numpy as np
import pandas as pd

from research.common.report import md_table, rel
from research.common.save import save
from research.wet_antenna import length_test as LTL
from research.waa_ksurface import data
from research.waa_ksurface import paths as P

C = data.branch_config()["a1"]
LT = data.LT


def length_class(L):
    """The length class name of every link length (config a1.length_classes_km)."""
    out = np.empty(len(L), object)
    for lo, hi in C["length_classes_km"]:
        out[(L >= lo) & (L < hi)] = f"L < {hi:g} km" if lo == 0 else (f"L > {lo:g} km" if hi > 50 else f"{lo:g}-{hi:g} km")
    return out


def length_test(D, le, variant):
    """The Part 1 length test (S3.1 settings, 1000 event-bootstrap draws) on one rain input."""
    t = LTL.steady_table(D, le, LT)
    return LTL.length_fits(t, int(D["le_event"].max()) + 1, variant, LT)


def compare(fits):
    """Every stored TRAIN cell of Part 1 against the rebuilt cell (the stored table repeats the observed I, s per model)."""
    ref = pd.read_csv(P.P1_LENGTH_TEST)
    ref = ref[ref.set == "TRAIN"].drop_duplicates(["input", "band", "class"])
    rows = []
    for r in ref.itertuples():
        m = fits[(fits.variant == r.input) & (fits.band == r.band) & (fits["class"] == r._5)]
        new = m.iloc[0] if len(m) else None
        row = {"input": r.input, "band": r.band, "class": r._5, "links_under_1km_part1": r.links_under_1km, "eligible_part1": r.eligible}
        for name, col in (("I", "I_obs"), ("I_lo", "I_lo"), ("I_hi", "I_hi"), ("s", "s_obs"), ("s_lo", "s_lo"), ("s_hi", "s_hi")):
            row[f"{name}_part1"] = getattr(r, col)
            row[f"{name}_rebuilt"] = float(new[name]) if new is not None else np.nan
        row["links_under_1km_rebuilt"] = int(new.links_under_1km) if new is not None else -1
        row["dI"], row["ds"] = abs(row["I_rebuilt"] - row["I_part1"]), abs(row["s_rebuilt"] - row["s_part1"])
        row["d_interval_ends_max"] = max(abs(row[f"{k}_rebuilt"] - row[f"{k}_part1"]) for k in ("I_lo", "I_hi", "s_lo", "s_hi"))
        rows.append(row)
    return pd.DataFrame(rows)


def counts(D1, le):
    """Link-events, links and minutes per band and length class (rain = radar / g_e, the Part 1 input V1)."""
    n5 = D1["le_n"] // 5
    valid = D1["w5"] > 0                                        # the 5 min bins with a valid A and radar
    rainy = valid & (D1["r5"] >= C["rainy_mmh"])
    le = le.assign(length_class=length_class(le.L_km.to_numpy()), valid_min=5 * np.add.reduceat(valid, D1["le_boff5"]),
                   rainy_min=5 * np.add.reduceat(rainy, D1["le_boff5"]))
    tab = le.groupby(["band_name", "length_class"]).agg(link_events=("cml_id", "size"), links=("cml_id", "nunique"), events=("event_id", "nunique"),
                                                         valid_minutes=("valid_min", "sum"), rainy_minutes=("rainy_min", "sum")).reset_index()
    # minutes per rain class on the links below 1 km (the Part 1 classes and the two classes A3a may add)
    short = np.repeat((le.L_km < 1.0).to_numpy(), n5)
    band = np.repeat(le.band_name.to_numpy(), n5)
    rows = []
    for lo, hi in [(0.0, 0.1)] + [tuple(c) for c in C["extra_classes_mmh"]] + [tuple(c) for c in LT["rain_classes_mmh"]] + [(20.0, 1e9)]:
        m = short & valid & (D1["r5"] >= lo) & (D1["r5"] < hi)
        rows.append({"rain_class_mmh": f"{lo:g}-{hi:g}" if hi < 1e8 else f">= {lo:g}", **{b: int(5 * (m & (band == b)).sum()) for b in D1["bands"]},
                     "all bands": int(5 * m.sum())})
    return tab, pd.DataFrame(rows)


def gauge_set(le):
    """The gauge-input set (link-events whose link has a 1 min City gauge within the Part 1 distance) and the g_e use."""
    g = data.g_table()
    has = le.gauge_id != ""
    sh = has & (le.L_km < 1.0)
    rows = [{"what": f"gauge-input set (City gauge within {data.DATA['gauge_max_dist_m']:g} m of the path): links", "count": int(le[has].cml_id.nunique())},
            {"what": "gauge-input set: link-events", "count": int(has.sum())},
            {"what": "gauge-input set: events", "count": int(le[has].event_id.nunique())},
            {"what": "gauge-input set, L < 1 km: links", "count": int(le[sh].cml_id.nunique())},
            {"what": "gauge-input set, L < 1 km: link-events", "count": int(sh.sum())},
            {"what": "TRAIN events with link-events", "count": int(le.event_id.nunique())},
            {"what": "events with a gauge ratio (g_e used)", "count": int(g.path_used.sum())},
            {"what": "events with g_e = 1 (fallback)", "count": int((~g.path_used).sum())},
            {"what": "link-events in events with a gauge ratio", "count": int(le.event_id.isin(g[g.path_used].event_id).sum())},
            {"what": "link-events in events with g_e = 1", "count": int((~le.event_id.isin(g[g.path_used].event_id)).sum())}]
    return pd.DataFrame(rows)


def stored_fits():
    """Where the TRAIN fits of the comparison models are, and whether they load (nothing is refitted)."""
    rows = []
    try:
        lit = json.loads(P.P1_LITERATURE_FITS.read_text(encoding="utf-8"))["fits"]
        for name in ("Schleiss b", "Pastorek b"):
            f = lit[name]
            rows.append({"model": name, "path": rel(P.P1_LITERATURE_FITS), "loads": True,
                         "content": f"shape {f['shape']}, amp per band, c_band, a_link of {len(f['a_link_per_link'])} links"})
    except Exception as e:
        rows.append({"model": "Schleiss b / Pastorek b", "path": rel(P.P1_LITERATURE_FITS), "loads": False, "content": f"{type(e).__name__}: {e}"})
    try:
        eq = json.loads(P.P1_EQUATION.read_text(encoding="utf-8"))
        rows.append({"model": "Part 1 equation", "path": rel(P.P1_EQUATION), "loads": True,
                     "content": f"every number of the equation, a_link of {len(eq['a_link']['per_link'])} links"})
    except Exception as e:
        rows.append({"model": "Part 1 equation", "path": rel(P.P1_EQUATION), "loads": False, "content": f"{type(e).__name__}: {e}"})
    fr = json.loads(P.P1_S8_FREEZE.read_text(encoding="utf-8"))["comparison_fits"]
    for key, name in (("B2", "no memory (B2 of S8)"), ("B0", "no WAA (B0 of S8)")):
        f = P.REPO / fr[key]["file"]
        rows.append({"model": name, "path": fr[key]["file"], "loads": bool(f.exists()),
                     "content": "pickle of the S8 fit" if f.exists() else "MISSING: the Part 1 fit cache is not on this machine; only its sha256 and its error numbers are committed"})
    return pd.DataFrame(rows)


def main():
    stop = []
    # ---- (a) the dataset
    D0, D1, le, drops, BANDS = data.train()
    n_links, n_le = int(le.cml_id.nunique()), len(le)
    ok_a = (n_links == C["expect"]["links"]) and (n_le == C["expect"]["link_events"])
    print(f"A1a: {n_links} links, {n_le} TRAIN link-events -> {'PASS' if ok_a else 'STOP'}", flush=True)
    funnel = pd.DataFrame([{"stage": k, "count": v} for k, v in drops.items()])
    save(funnel, P.A1_FUNNEL, index=False)
    if not ok_a:
        stop.append("A1a")
    # ---- (b) the length test
    cmp_ = fits = None
    if ok_a:
        fits = pd.concat([length_test(D0, le, "V0"), length_test(D1, le, "V1")], ignore_index=True)
        cmp_ = compare(fits)
        save(cmp_.round(5), P.A1_LENGTH_CMP, index=False)
        worst = float(np.nanmax(np.r_[cmp_.dI, cmp_.ds]))
        ok_b = bool(len(cmp_) and cmp_[["dI", "ds"]].notna().all().all() and worst <= C["tolerance"])
        print(f"A1b: {len(cmp_)} stored cells, largest |dI| {cmp_.dI.max():.5f} dB, largest |ds| {cmp_.ds.max():.5f} dB/km -> {'PASS' if ok_b else 'STOP'}", flush=True)
        if not ok_b:
            stop.append("A1b")
    # ---- (c) counts, (d) stored fits, (e) the 32 GHz rule
    tab, short_min = counts(D1, le)
    gs = gauge_set(le)
    sf = stored_fits()
    for t, f in ((tab, P.A1_COUNTS), (short_min, P.A1_SHORT_MINUTES), (gs, P.A1_GAUGE), (sf, P.A1_FITS)):
        save(t, f, index=False)
    elig = None
    if fits is not None:
        v1 = fits[(fits.variant == "V1") & fits.I.notna()]
        elig = v1.assign(eligible=v1.links_under_1km >= LT["verdict"]["min_short_links"]).groupby("band").eligible.sum().reindex(BANDS).fillna(0).astype(int)
        print("A1e: eligible cells per band (V1):", elig.to_dict(), flush=True)

    # ---- the report
    a = []
    a.append("# A1: data audit and reproduction gate (TRAIN only)\n")
    a.append(f"Written by `research/waa_ksurface/a1_audit.py`. VAL and TEST attenuation was not opened.\n")
    a.append(f"## (a) The rebuilt dataset: {'PASS' if ok_a else 'STOP'}\n")
    a.append(f"{n_links} links and {n_le} TRAIN link-events ({le.event_id.nunique()} events); the gate needs {C['expect']['links']} and {C['expect']['link_events']}. "
             f"Built with `system_id.build_dataset` called as in `system_id_s3.dataset`. The funnel (`{rel(P.A1_FUNNEL)}`):\n")
    a.append(md_table(funnel))
    if cmp_ is not None:
        a.append(f"\n## (b) The Part 1 length test reproduced: {'PASS' if 'A1b' not in stop else 'STOP'}\n")
        a.append(f"{len(cmp_)} stored TRAIN cells (inputs: {", ".join(sorted(cmp_.input.unique()))}) of `{rel(P.P1_LENGTH_TEST)}`. Largest difference: I {cmp_.dI.max():.5f} dB, s {cmp_.ds.max():.5f} dB/km "
                 f"(allowed {C['tolerance']}); interval ends {cmp_.d_interval_ends_max.max():.5f}. Table: `{rel(P.A1_LENGTH_CMP)}`.\n")
        a.append(md_table(cmp_[["input", "band", "class", "I_part1", "I_rebuilt", "s_part1", "s_rebuilt", "dI", "ds", "links_under_1km_rebuilt", "eligible_part1"]].round(4)))
        a.append(f"""
**The Part 1 definition of a steady-rain 5 min value** (config `wet_antenna.system_id_s3.length_test`, code `length_test.steady_table`):
a 5 min bin of a link-event with a valid A (at least 3 of its 5 minutes) and a valid radar path-mean R, whose R lies in a rain
class ({', '.join(f'{lo:g}-{hi:g}' for lo, hi in LT['rain_classes_mmh'])} mm/h, [lo, hi)), where (1) the {LT['wet_before_min']} min before it
({LT['wet_before_min'] // 5} previous bins) are all wet (R >= {LT['wet_threshold_mmh']} mm/h), and (2) the mean R of the {LT['prev_window_min']} min before it
({LT['prev_window_min'] // 5} previous bins, all with radar) is within {LT['prev_ratio'][0]}-{LT['prev_ratio'][1]} times the bin's R. Per link-event and class the kept bins are
reduced to one value: the median of (A - c) and the median of R.

**The Part 1 definition of an eligible cell** (a cell = band x rain class): the cell is fitted if it has at least
{LT['min_link_events']} link-event values and 3 distinct link lengths (line (A - c) = I + s L, Huber k = {LT['huber_k']}); it is ELIGIBLE if it also has at least
{LT['verdict']['min_short_links']} distinct links shorter than {LT['verdict']['short_km']:g} km.
""")
    a.append("\n## (c) The TRAIN data\n")
    a.append(f"Rain = radar / g_e (V1). A valid minute has A and radar; a rainy minute has R >= {C['rainy_mmh']} mm/h. (`{rel(P.A1_COUNTS)}`)\n")
    a.append(md_table(tab))
    a.append(f"\nValid minutes per rain class on the links below 1 km (`{rel(P.A1_SHORT_MINUTES)}`):\n")
    a.append(md_table(short_min))
    a.append(f"\nThe gauge-input set and the g_e use (`{rel(P.A1_GAUGE)}`):\n")
    a.append(md_table(gs))
    a.append(f"\n## (d) The stored comparison fits (nothing refitted; `{rel(P.A1_FITS)}`)\n")
    a.append(md_table(sf))
    if elig is not None:
        a.append("\n## (e) The 32 GHz rule\n")
        a.append("Eligible cells per band (V1, Part 1 classes): " + ", ".join(f"{b}: {n}" for b, n in elig.items()) + ".\n")
        out = [b for b, n in elig.items() if n == 0]
        a.append(("**Excluded from this branch for every model: " + ", ".join(out) + "** (no eligible length-test cell, so F and P cannot be measured)." if out
                  else "Every band has at least one eligible cell: no band is excluded."))
    if stop:
        a.append(f"\n## STOP\n\nFailed gate: {', '.join(stop)}. Nothing after it was run.")
    save("\n".join(a) + "\n", P.A1_REPORT)
    print(f"report: {rel(P.A1_REPORT)}", flush=True)
    if stop:
        sys.exit(1)


if __name__ == "__main__":
    main()
