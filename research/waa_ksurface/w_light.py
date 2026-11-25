"""Light-rain overestimate: a diagnostic (decisions.md, "Light-rain overestimate: diagnostic", logged before any code was run).

    python -m research.waa_ksurface.w_light

Reads ONLY the saved held-out predictions of the WAA-only comparison (research/.cache/waa_ksurface/w_cv/preds.npz, path A:
total and WAA part of the 6 models), the measured 5 min A - c, the radar rain input R = radar / g_e and the City gauge series
of the TRAIN setting (w_cv.setting()). No model is fitted; VAL and TEST are not opened.

Checks: (1) residual length test, median residual per link-event = a + b L (robust line), per model and rain class;
(2) tail split of the light-rain residual (after heavy rain or not), ours and Part 1; (3) median residual per model and class;
(4) gauge check on short links with a City gauge; (5) one figure. Intervals: event bootstrap (the 14 events resampled).
"""
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import data, paths as P, w_cv as W

CL = data.branch_config()["w_light"]
MODELS = W.MODELS
WAA_MODELS = [m for m in MODELS if m != "no WAA"]
CLASSES = [tuple(c) for c in CL["classes_mmh"]]                # half-open [lo, hi) in mm/h
CNAME = [f"{lo:g}-{hi:g}" for lo, hi in CLASSES]
LIGHT = [i for i, (lo, hi) in enumerate(CLASSES) if hi <= CL["light_max_mmh"]]
MARK = {"ours": "o", "Part 1 equation": "s", "Schleiss": "^", "Pastorek": "D", "no memory": "v", "no WAA": "x"}
LS = {"ours": "-", "Part 1 equation": "--", "Schleiss": "-.", "Pastorek": ":", "no memory": (0, (5, 1, 1, 1)), "no WAA": (0, (1, 1))}


def log(s):
    print(f"{time.strftime('%H:%M:%S')} {s}", flush=True)


# =================================================================================================== weighted statistics
def wmedian(v, w):
    """Median of v with integer weights w (the event-bootstrap multiplicity of each value). NaN if no weight."""
    o = np.argsort(v, kind="stable")
    cw = np.cumsum(w[o])
    if cw.size == 0 or cw[-1] <= 0:
        return np.nan
    return float(v[o][np.searchsorted(cw, 0.5 * cw[-1])])


def boot_median(v, ev, draws):
    """(median, lo, hi) of v; the interval from the same event draws everywhere (draws: reps x events)."""
    o = np.argsort(v, kind="stable")
    vs, es = v[o], ev[o]
    out = np.empty(len(draws))
    for i, c in enumerate(draws):
        cw = np.cumsum(c[es])
        out[i] = vs[np.searchsorted(cw, 0.5 * cw[-1])] if cw.size and cw[-1] > 0 else np.nan
    lo, hi = np.nanpercentile(out, CL["interval"])
    return (float(np.median(v)) if v.size else np.nan), float(lo), float(hi), out


def huber_line(x, y, w):
    """Robust line y = a + b x: Huber M-estimate by IRLS, tuning 1.345 x MAD scale, with case weights w."""
    X = np.c_[np.ones_like(x), x]
    keep = w > 0
    X, y, w = X[keep], y[keep], w[keep].astype(float)
    if len(y) < 3 or np.ptp(X[:, 1]) == 0:
        return np.nan, np.nan
    beta = np.linalg.lstsq(X * np.sqrt(w)[:, None], y * np.sqrt(w), rcond=None)[0]
    for _ in range(CL["huber_iter"]):
        r = y - X @ beta
        s = 1.4826 * np.median(np.abs(r - np.median(r))) or 1e-9
        k = CL["huber_c"] * s
        u = np.where(np.abs(r) <= k, 1.0, k / np.maximum(np.abs(r), 1e-12))
        sw = np.sqrt(w * u)
        nb = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
        if np.allclose(nb, beta, atol=1e-8):
            beta = nb
            break
        beta = nb
    return float(beta[0]), float(beta[1])


def ci(vals):
    lo, hi = np.nanpercentile(vals, CL["interval"])
    return float(lo), float(hi)


def fmt(x, lo, hi, d=3):
    return f"{x:+.{d}f} [{lo:+.{d}f}, {hi:+.{d}f}]"


def above(lo):
    return lo > 0


# =================================================================================================== the data
def load():
    st = W.setting()
    z = np.load(W.CACHE / "preds.npz")
    assert np.allclose(z["target"], st["target"], equal_nan=True), "preds.npz does not line up with the TRAIN setting"
    D, le = st["D"], st["le"]
    n5 = (D["le_n"] // 5).astype(int)
    o5 = np.r_[0, np.cumsum(n5)]
    R5 = W.to5(st["R1"])
    g1 = np.concatenate([D["rain_gauge"][o + W.PADM:o + W.PADM + k] for o, k in zip(D["le_roff"], D["le_n"])]).reshape(-1, 5)
    G5 = np.where(np.isfinite(g1).all(axis=1), g1.mean(axis=1), np.nan)
    tot = {m: z[f"total|A|{m}"] for m in MODELS}
    waa = {m: z[f"waa|A|{m}"] for m in MODELS}
    res = {m: tot[m] - st["target"] for m in MODELS}
    # previous 12 bins of the same window: heavy rain seen? full history?
    heavy_prev = np.zeros(len(R5), bool)
    full_hist = np.zeros(len(R5), bool)
    nb = CL["tail_bins"]
    for e in range(len(le)):
        h = (R5[o5[e]:o5[e + 1]] >= CL["heavy_mmh"]).astype(int)
        cs = np.r_[0, np.cumsum(h)]
        k = np.arange(len(h))
        heavy_prev[o5[e]:o5[e + 1]] = (cs[k] - cs[np.maximum(k - nb, 0)]) > 0
        full_hist[o5[e]:o5[e + 1]] = k >= nb
    cls = np.full(len(R5), -1)
    for i, (lo, hi) in enumerate(CLASSES):
        cls[(R5 >= lo) & (R5 < hi)] = i
    rng = np.random.default_rng(CL["seed"])
    nev = len(st["events"])
    draws = np.stack([np.bincount(rng.integers(0, nev, nev), minlength=nev) for _ in range(CL["boot_reps"])])
    log(f"setting: {len(le)} link-events, {nev} events, {int(st['scored'].sum())} scored bins; "
        f"{int((st['scored'] & (cls >= 0)).sum())} scored bins during rain in the 5 classes")
    return {"st": st, "le": le, "L": st["L"], "ev_le": st["ev_le"], "bin_le": st["bin_le"], "ev_bin": st["ev_bin"], "sc": st["scored"],
            "R5": R5, "G5": G5, "cls": cls, "res": res, "waa": waa, "heavy_prev": heavy_prev, "full_hist": full_hist, "draws": draws}


# =================================================================================================== the checks
def check1(X):
    """Per model and class: per link-event median residual, then a + b L (robust), event bootstrap."""
    rows = []
    for m in MODELS:
        for i, c in enumerate(CNAME):
            mask = X["sc"] & (X["cls"] == i)
            le_i = X["bin_le"][mask]
            r = X["res"][m][mask]
            df = pd.DataFrame({"le": le_i, "r": r}).groupby("le").r.agg(["median", "size"])
            df = df[df["size"] >= CL["min_bins_per_le"]]
            idx = df.index.to_numpy()
            x, y, ev = X["L"][idx], df["median"].to_numpy(), X["ev_le"][idx]
            a, b = huber_line(x, y, np.ones(len(y)))
            boot = np.array([huber_line(x, y, c_[ev]) for c_ in X["draws"]])
            alo, ahi = ci(boot[:, 0])
            blo, bhi = ci(boot[:, 1])
            if alo > 0 and bhi > 0 and blo > 0:
                read = "both"
            elif alo > 0 and not blo > 0:
                read = "antenna"
            elif blo > 0 and not alo > 0:
                read = "path"
            else:
                read = "neither"
            rows.append({"model": m, "class (mm/h)": c, "link-events": len(y), "a (dB)": a, "a lo": alo, "a hi": ahi,
                         "b (dB/km)": b, "b lo": blo, "b hi": bhi, "reading": read})
        log(f"   (1) {m}: done")
    return pd.DataFrame(rows)


def check2(X):
    light = X["sc"] & np.isin(X["cls"], LIGHT)
    groups = {"(i) within 60 min after R >= 2": light & X["heavy_prev"],
              "(ii) no R >= 2 in the previous 60 min": light & ~X["heavy_prev"] & X["full_hist"]}
    rows, boots = [], {}
    w = X["waa"]["ours"]
    for g, mask in groups.items():
        ev = X["ev_bin"][mask]
        row = {"group": g, "bins": int(mask.sum()), "minutes": int(mask.sum()) * 5}
        for m in ("ours", "Part 1 equation"):
            med, lo, hi, bt = boot_median(X["res"][m][mask], ev, X["draws"])
            row[f"{m} residual (dB)"], row[f"{m} lo"], row[f"{m} hi"] = med, lo, hi
            boots[(g, m)] = bt
        med, lo, hi, _ = boot_median(w[mask], ev, X["draws"])
        row["our WAA part (dB)"], row["WAA lo"], row["WAA hi"] = med, lo, hi
        low = w[mask] < CL["w_low_db"]
        sh = np.array([np.sum(c[ev] * low) / max(np.sum(c[ev]), 1) for c in X["draws"]])
        row["minutes W < 0.3 dB"] = int(low.sum()) * 5
        row["share W < 0.3 dB"], (row["share lo"], row["share hi"]) = float(low.mean()) if low.size else np.nan, ci(sh)
        rows.append(row)
    gi, gii = list(groups)
    diff = {}
    for m in ("ours", "Part 1 equation"):
        d = boots[(gi, m)] - boots[(gii, m)]
        diff[m] = (rows[0][f"{m} residual (dB)"] - rows[1][f"{m} residual (dB)"],) + ci(d)
    return pd.DataFrame(rows), diff


def check3(X):
    rows = []
    for i, c in enumerate(CNAME):
        mask = X["sc"] & (X["cls"] == i)
        ev = X["ev_bin"][mask]
        for m in MODELS:
            med, lo, hi, _ = boot_median(X["res"][m][mask], ev, X["draws"])
            rows.append({"class (mm/h)": c, "model": m, "bins": int(mask.sum()), "median residual (dB)": med, "lo": lo, "hi": hi,
                         "overestimate": bool(lo > 0)})
    df = pd.DataFrame(rows)
    read = {}
    for i, c in enumerate(CNAME):
        o = df[df["class (mm/h)"] == c].set_index("model").overestimate
        others = sum(o[m] for m in WAA_MODELS if m != "ours")
        if all(o[m] for m in MODELS):
            read[c] = "shared by all"
        elif o["ours"] and others >= 3 and not o["no WAA"]:
            read[c] = "shared by the WAA models"
        elif o["ours"] and others + int(o["no WAA"]) <= 1:
            read[c] = "ours only"
        else:
            read[c] = "mixed"
    return df, read


def check4(X):
    le_short_g = ((X["le"].L_km < CL["short_km"]) & (X["le"].gauge_id != "")).to_numpy()
    mask = X["sc"] & np.isin(X["cls"], LIGHT) & le_short_g[X["bin_le"]] & np.isfinite(X["G5"])
    ev = X["ev_bin"][mask]
    ro = X["G5"][mask] < CL["gauge_dry_mmh"]                     # radar > 0.1 in every light bin by construction
    share = float(ro.mean()) if ro.size else np.nan
    sh = np.array([np.sum(c[ev] * ro) / max(np.sum(c[ev]), 1) for c in X["draws"]])
    rows = []
    for m in MODELS:
        r = X["res"][m][mask]
        a = boot_median(r[ro], ev[ro], X["draws"])
        b = boot_median(r[~ro], ev[~ro], X["draws"])
        d = a[3] - b[3]
        rows.append({"model": m, "radar-only residual (dB)": a[0], "ro lo": a[1], "ro hi": a[2], "others residual (dB)": b[0],
                     "ot lo": b[1], "ot hi": b[2], "difference (dB)": a[0] - b[0], "diff lo": ci(d)[0], "diff hi": ci(d)[1]})
    info = {"link-events": int(len(np.unique(X["bin_le"][mask]))), "links": int(X["le"].cml_id[np.unique(X["bin_le"][mask])].nunique()),
            "bins": int(mask.sum()), "radar-only bins": int(ro.sum()), "share": share, "share lo": ci(sh)[0], "share hi": ci(sh)[1]}
    return pd.DataFrame(rows), info


# =================================================================================================== the figure
def figure(c1, c3):
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.6), sharex=True)
    x = np.arange(len(CNAME))
    for j, m in enumerate(MODELS):
        d = c3[c3.model == m].set_index("class (mm/h)").loc[CNAME]
        a1.plot(x, d["median residual (dB)"], color=W.COLOR[m], ls=LS[m], marker=MARK[m], ms=6, lw=1.8, label=m)
        e = c1[c1.model == m].set_index("class (mm/h)").loc[CNAME]
        xo = x + (j - 2.5) * 0.11
        a2.errorbar(xo, e["a (dB)"], yerr=[e["a (dB)"] - e["a lo"], e["a hi"] - e["a (dB)"]], fmt=MARK[m], color=W.COLOR[m], ms=5,
                    capsize=2, lw=1.2, label=m)
    for ax in (a1, a2):
        ax.axhline(0, color="#555555", lw=0.8)
        ax.set_xticks(x, CNAME)
        ax.set_xlabel("radar rain input class R (mm/h)", fontsize=9)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.grid(axis="y", color="#ececec", lw=0.6)
        ax.tick_params(labelsize=8)
    a1.set_ylabel("median residual, model - measured (dB)", fontsize=9)
    a1.set_title("residual per rain class (all scored bins, held-out)", fontsize=9, loc="left")
    a2.set_ylabel("antenna share a (dB), 90 % event bootstrap", fontsize=9)
    a2.set_title("robust line per link-event: median residual = a + b L", fontsize=9, loc="left")
    a1.legend(fontsize=7, frameon=False)
    fig.suptitle("Light-rain overestimate (TRAIN, leave-one-event-out predictions of the WAA-only comparison, path A)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, P.WL_FIGURE.with_suffix(".png"), dpi=300)
    save(fig, P.WL_FIGURE.with_suffix(".pdf"))
    plt.close(fig)


# =================================================================================================== report
def main():
    t0 = time.time()
    X = load()
    c1 = check1(X)
    log(f"=== (1) residual length test done ({(time.time() - t0) / 60:.1f} min) ===")
    c2, diff = check2(X)
    c3, read3 = check3(X)
    c4, info4 = check4(X)
    log(f"=== (2)-(4) done ({(time.time() - t0) / 60:.1f} min) ===")
    figure(c1, c3)
    save(c1.round(4), P.WL_LENGTH, index=False)
    save(c2.round(4), P.WL_TAIL, index=False)
    save(c3.round(4), P.WL_CLASSES, index=False)
    save(c4.round(4), P.WL_GAUGE, index=False)

    t1 = c1.assign(**{"a [90 %]": [fmt(r["a (dB)"], r["a lo"], r["a hi"]) for _, r in c1.iterrows()],
                      "b [90 %]": [fmt(r["b (dB/km)"], r["b lo"], r["b hi"], 4) for _, r in c1.iterrows()]})[
        ["model", "class (mm/h)", "link-events", "a [90 %]", "b [90 %]", "reading"]]
    t2 = pd.DataFrame([{"group": r.group, "minutes": r.minutes,
                        "ours residual": fmt(r["ours residual (dB)"], r["ours lo"], r["ours hi"]),
                        "Part 1 residual": fmt(r["Part 1 equation residual (dB)"], r["Part 1 equation lo"], r["Part 1 equation hi"]),
                        "our WAA part": fmt(r["our WAA part (dB)"], r["WAA lo"], r["WAA hi"]),
                        "minutes with W < 0.3 dB": f"{r['minutes W < 0.3 dB']} ({100 * r['share W < 0.3 dB']:.0f} % [{100 * r['share lo']:.0f}, {100 * r['share hi']:.0f}])"}
                       for _, r in c2.iterrows()])
    t3 = c3.pivot(index="model", columns="class (mm/h)", values="median residual (dB)").loc[MODELS, CNAME]
    t3s = t3.copy().astype(object)
    for _, r in c3.iterrows():
        t3s.loc[r.model, r["class (mm/h)"]] = fmt(r["median residual (dB)"], r.lo, r.hi) + (" *" if r.overestimate else "")
    t4 = pd.DataFrame([{"model": r.model, "radar-only": fmt(r["radar-only residual (dB)"], r["ro lo"], r["ro hi"]),
                        "others": fmt(r["others residual (dB)"], r["ot lo"], r["ot hi"]),
                        "difference": fmt(r["difference (dB)"], r["diff lo"], r["diff hi"])} for _, r in c4.iterrows()])

    # readings (the logged rules)
    ours1 = c1[c1.model == "ours"].set_index("class (mm/h)")
    r1 = {c: ours1.reading[c] for c in CNAME}
    d_ours = diff["ours"]
    tail = d_ours[1] > 0
    ii = c2.iloc[1]
    o4 = c4.set_index("model").loc["ours"]
    inp = o4["diff lo"] > 0
    S1 = ("Check 1 (ours, light classes " + ", ".join(f"{c}: {r1[c]}" for c in [CNAME[i] for i in LIGHT]) +
          f"; heavier: " + ", ".join(f"{c}: {r1[c]}" for c in [CNAME[i] for i in range(len(CNAME)) if i not in LIGHT]) + ").")
    S2 = (f"Check 2: our (i) - (ii) = {fmt(*d_ours)} dB -> {'tail (carried over from heavy rain: antenna drying)' if tail else 'not a tail'}; "
          f"group (ii) residual of ours {fmt(ii['ours residual (dB)'], ii['ours lo'], ii['ours hi'])} dB -> "
          f"{'the overestimate is also present without heavy rain before' if ii['ours lo'] > 0 else 'no overestimate without heavy rain before'}.")
    S3 = "Check 3: " + ", ".join(f"{CNAME[i]}: {read3[CNAME[i]]}" for i in LIGHT) + "."
    S4 = (f"Check 4: radar-only share {100 * info4['share']:.0f} % [{100 * info4['share lo']:.0f}, {100 * info4['share hi']:.0f}]; our difference "
          f"{fmt(o4['difference (dB)'], o4['diff lo'], o4['diff hi'])} dB -> {'input (the overestimate sits where the gauge sees no rain)' if inp else 'not input by this rule'}.")
    L = ["# Light-rain overestimate: a diagnostic\n",
         "Written by `research/waa_ksurface/w_light.py` (decisions.md, \"Light-rain overestimate: diagnostic\", logged before any code was run). "
         "Saved held-out predictions of the WAA-only comparison only (path A), TRAIN (14 learning events, 4563 link-events); no model fitted; "
         f"VAL and TEST not opened. Residual = model - measured (5 min A - c), scored bins. Intervals: event bootstrap, {CL['boot_reps']} draws, "
         f"seed {CL['seed']}, 90 %.\n",
         "## 1. Residual length test (per link-event median residual = a + b L, robust line)\n", md_table(t1), "",
         "## 2. Tail split (light rain, 0.1-1 mm/h)\n", md_table(t2), "",
         f"(i) - (ii): ours {fmt(*diff['ours'])} dB; Part 1 {fmt(*diff['Part 1 equation'])} dB.\n",
         "## 3. Median residual per model and class (* = interval above 0)\n", md_table(t3s.reset_index()), "",
         "Reading per class: " + "; ".join(f"{c}: {read3[c]}" for c in CNAME) + ".\n",
         "## 4. Gauge check (links < 1 km with a City gauge within 1 km; light-rain bins with a valid gauge 5 min rate)\n",
         f"{info4['links']} links, {info4['link-events']} link-events, {info4['bins']} bins; radar-only (gauge < {CL['gauge_dry_mmh']} mm/h): "
         f"{info4['radar-only bins']} bins = {100 * info4['share']:.1f} % [{100 * info4['share lo']:.1f}, {100 * info4['share hi']:.1f}].\n",
         md_table(t4), "",
         "## 5. Figure\n", f"`{rel(P.WL_FIGURE.with_suffix('.png'))}`\n",
         "## Readings (the logged rules)\n", "- " + S1, "- " + S2, "- " + S3, "- " + S4, "",
         f"Run time {(time.time() - t0) / 60:.1f} min."]
    save("\n".join(L) + "\n", P.WL_REPORT)
    log(f"report: {rel(P.WL_REPORT)}")
    log(f"=== light-rain diagnostic ended ({(time.time() - t0) / 60:.1f} min) ===")


if __name__ == "__main__":
    main()
