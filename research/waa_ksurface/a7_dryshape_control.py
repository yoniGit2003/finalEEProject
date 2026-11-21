"""Control for the within-window ratio (decisions.md, "Control for the within-window ratio", logged before the run).

Can baseline offsets plus the window selection alone make plain exponential drying (a K0 truth) look like the
measured slowing (ratio rate(late) / rate(early) = 0.23)? Fake data of the step 3c / A7 fake world, no fit at all;
the windows and the ratio are exactly those of the diagnostic (a7_dryshape.windows, a7_dryshape.ratios).

    python -m research.waa_ksurface.a7_dryshape_control

Truths: K0 with the two A7 K0 pairs (1/k_dry 8.88 min with 1/k_wet 1.74 min; 13.22 min with 10.92 min); each with
the offset error on (sd 0.136 dB per link-event, as step 3c) and off (the same dataset, offset sd 0).
"""
import numpy as np
import pandas as pd

from research.common.report import md_table
from research.common.save import save
from research.waa_ksurface import a7_dryshape as ds, a7_review as rv, data, fake6, s1_separability as s1, s3_recovery as s3
from research.waa_ksurface import paths as P
from research.wet_antenna.system_id import PADM

CC = data.branch_config()["a7c"]
log = rv.log


def fake(S, truth, offsets):
    """The step 3c fake dataset of a K0 truth (as s3c_a7_variant.run_one), with the offset error on or off."""
    D, le, learn = S["D"], S["le"], S["learn"]
    spec = {"kind": "const", "k_wet": 1.0 / truth["tau_wet_min"], "k_dry": 1.0 / truth["tau_dry_min"]}
    sc = np.repeat(np.where(learn, 1.0, s1.S1["radar_scale_added"]), D["le_n"] + 2 * PADM)
    Dt = dict(D)
    Dt["rain"] = np.asarray(D["rain"], float) / sc
    keep = fake6.A6["noise"]
    fake6.A6["noise"] = {**keep, "offset_sd_db": keep["offset_sd_db"] if offsets else 0.0}     # same random numbers either way
    try:
        Dk, _ = fake6.make(Dt, le, S["n_cells"], spec, S["curves"], S["iF"], S["iP"], truth["seed"], S["step"], S["dist"],
                           ghost_pool=np.sort(le.event_id[learn].unique()))
    finally:
        fake6.A6["noise"] = keep
    Dk["rain"] = Dk["rain"] * sc
    return Dk


def ratio_rows(S, Dk, name, offsets):
    use = S["learn"].astype(float)
    le_w, idx, target = ds.windows(Dk, fake6.window(Dk, np.asarray(Dk["rain"], float)), use)
    Y = target[idx]
    _, _, r = ds.ratios(Y)
    ev = S["le"].event_id.to_numpy()[le_w]
    ok = np.isfinite(r)
    rows = []
    cls = ds.C["class_db"]
    for c, m in [("all", ok), (f">= {cls:g} dB", ok & (Y[:, 0] >= cls)), (f"< {cls:g} dB", ok & (Y[:, 0] < cls))]:
        row = {"truth": name, "offsets": "on" if offsets else "off", "class": c, "windows kept (A - c > 0.15)": int(len(le_w)),
               "windows with ratio": int(m.sum()), "events": int(len(np.unique(ev[m])))}
        if m.sum():
            rr = r[m]
            lo, hi = ds.boot(ev[m], [lambda q: np.median(rr[q])], ds.C["boot_reps"], ds.C["boot_seed"])
            row.update({"median ratio": float(np.median(rr)), "5 %": float(lo[0]), "95 %": float(hi[0])})
        rows.append(row)
    log(f"   (control) {name}, offsets {'on' if offsets else 'off'}: {len(le_w)} windows kept, {int(ok.sum())} with a ratio, "
        f"median {np.median(r[ok]) if ok.any() else np.nan:.3f}")
    return rows


def reading(T, real):
    """The pre-registered reading per truth (offsets on)."""
    out = []
    for name, est in (("G", "gauge"), ("R", "radar")):
        t = T[(T.truth == name) & (T.offsets == "on")].set_index("class")
        a = t.loc["all"]
        contains1 = bool(a["5 %"] <= 1 <= a["95 %"])
        below = {c: bool(real[(est, c)] < t.loc[c, "5 %"]) for c in t.index if c != "all"}
        ne = contains1 or all(below.values())
        out.append({"truth": name, "fake interval (all) contains 1": contains1,
                    **{f"real {real[(est, c)]:.3f} below fake interval ({c})": v for c, v in below.items()},
                    "real 0.227 below fake interval (all)": bool(0.227 < a["5 %"]),
                    "reading": "not explained by offset and selection" if ne else "explained: stop to rethink"})
    return pd.DataFrame(out)


def main():
    S = s3.setting()
    rows = []
    for name, truth in CC["truths"].items():
        for offsets in (True, False):
            rows += ratio_rows(S, fake(S, truth, offsets), name, offsets)
    T = pd.DataFrame(rows)
    real_sum = pd.read_csv(P.A7D_SUM)
    real = {(r["k_wet estimate"], r["class"]): r["ratio measured"] for _, r in real_sum.iterrows()}
    Rd = reading(T, real)
    save(T, P.A7C_TABLE, index=False)
    out = T.copy()
    out["median ratio (90 %)"] = [f"{m:.3f} ({a:.3f} to {b:.3f})" if np.isfinite(m) else "" for m, a, b in
                                  zip(out.get("median ratio", np.nan), out.get("5 %", np.nan), out.get("95 %", np.nan))]
    out = out.drop(columns=[c for c in ("median ratio", "5 %", "95 %") if c in out])
    L = ["# Control for the within-window ratio\n",
         "Written by `research/waa_ksurface/a7_dryshape_control.py` (decisions.md, \"Control for the within-window ratio\", logged before the run). "
         "Fake data only (the step 3c / A7 fake world), no fit. Windows and ratio exactly as in `a7_dryshape`.\n",
         "Truths (K0): " + "; ".join(f"{n}: 1/k_dry {t['tau_dry_min']} min, 1/k_wet {t['tau_wet_min']} min, seed {t['seed']}" for n, t in CC["truths"].items()) + ".\n",
         "## Ratios\n", md_table(out), "", "## Reading (pre-registered; offsets on)\n", md_table(Rd), ""]
    save("\n".join(L) + "\n", P.A7C_REPORT)
    log(f"report: {P.A7C_REPORT}")


if __name__ == "__main__":
    main()
