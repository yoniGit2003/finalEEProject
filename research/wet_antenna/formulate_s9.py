"""wet_antenna step S9: formulate Part 1. No new analysis, no refitting, no new reading of TEST.

Writes, from the frozen S8b closed-F model (the S6 plug-back fit "F and P", fitted on TRAIN before the S8 freeze):
  - outputs/products/waa_equation_final.json   the official equation: form, every number with its 90 % interval, fixed
                                                dynamics, lag, rain input, valid range, bands, not identified, the fit's sha256
  - outputs/products/waa_equation_reference.json  3 link-events of the S8 TEST cache (already read once) with the model's 1 min
                                                W and path term, the reference of the unit test of waa_equation.py
  - the tables behind the figures of docs/WAA_research_overview.pdf (made by docs/WAA_research_overview_figures/figs.py): an
    example event, the memory curve, the closed F with its band, the dry-path curves (and the three boxes, kept as a table).
  - outputs/reports/formulate_s9.md, the report of this step.
The check of the overview and the document are written by hand (outputs/reports/formulate_s9_check.md, docs/WAA_research_overview.tex).
Settings: config `wet_antenna.formulate_s9`. Run from the repository root:
    python -m research.wet_antenna.formulate_s9
"""
import hashlib
import json
import pickle
import sys
import time

import numpy as np
import pandas as pd

from research.common.report import rel
from research.common.save import save
from research.wet_antenna import length_test as LTL
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.event_factors import rescaled
from research.wet_antenna.final_model import DATA, FIX, WA
from research.wet_antenna.validation import LB, dry_windows, lclass, load, radar_max, w_minutes
from research.wet_antenna.waa_equation import WAAEquation

S = WA["formulate_s9"]
PADM = sid.PADM


def kernel_parts(D, le, fz, idx):
    """The fit's own 1 min W and path term (c = 0) of the link-events idx (system_id kernel, agg = 1), in the order of the dataset."""
    idx = sorted(idx)
    Ds = sid.subset(D, np.isin(np.arange(len(le)), idx), renumber_events=True)
    les = le.iloc[idx].reset_index(drop=True)
    D0 = dict(Ds)
    D0["le_c"] = np.zeros(len(les))
    Pm = sid.Problem2(D0, "M4b", fz["M"], agg=1)
    x, _ = LTL.x_for(Pm, fz, LTL.link_ids(Ds, les))
    xW, xP = x.copy(), x.copy()
    xW[Pm.names == Pm.pname] = 0.0
    xP[Pm.names == Pm.fname] = 0.0
    return Ds, les, Pm.predict(xW, fz["d"]), Pm.predict(xP, fz["d"])


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    clock = time.time()
    d, BANDS, sets, gmap, mods, free = load()
    Dt1, let = sets["TRAIN"]
    C, z, fc = mods["closed F"]["C"], np.asarray(mods["closed F"]["z"]), mods["closed F"]["fz"]
    nb = len(BANDS)
    fit_file = P.S6_CACHE / "plug_shared_FP_F_and_P.pkl"
    fit_sha = hashlib.sha256(fit_file.read_bytes()).hexdigest()
    tr = pd.read_csv(P.S8B_TRAIN).set_index("quantity")

    # ------------------------------------------------------------------------------------------- the equation file
    def num(v, q):
        return {"value": float(v), "lo": float(tr.boot_lo[q]), "hi": float(tr.boot_hi[q])}
    _, _, a = C.values(z)
    ids = LTL.link_ids(Dt1, let)
    eq = {"name": "WAA equation of Part 1 (closed form adopted in S8b)",
          "form": {"A": "A(t) - c = W_f(t) + W_s(t) + L c_band R(t - d)",
                   "W_f": "dW_f/dt = (phi u - W_f) / tau_f, tau_f = tau_f_wet while W_f rises (phi u > W_f), tau_f_dry while it falls",
                   "W_s": "dW_s/dt = ((1 - phi) u - W_s) / tau_s",
                   "u": "u = a_link A_band (R(t - d) / 5)^b",
                   "update": "once per minute: W <- target + (W - target) exp(-1 / tau); the output after the update"},
          "rain_input": "R = radar path-mean rain (release radar_along_links_5min R_mean, held over its 5 minutes) / g_e; g_e = radar / City-gauge rain of "
                        "the event over its gauge links (City gauge within 1 km of the path), used if the gauge rain is >= 1.0 mm per gauge link and the "
                        "ratio is above 0, else g_e = 1 (S4.3, S8)",
          "offset_c": "c = the median of A over the dry minutes 120 to 30 min before the link's first wet bin (>= 30 such minutes, else 0); A = TL - "
                      "baseline (baseline v2: the median dry TL of the 24 h before the event)",
          "bands_ghz": {b: list(DATA["bands_ghz"][b]) for b in BANDS},
          "A_band_db": {b: num(z[i], f"A[{b}]") for i, b in enumerate(BANDS)},
          "b": num(z[nb], "b"),
          "c_band_db_km_per_mmh": {b: num(z[nb + 1 + i] / 5, f"c[{b}]") for i, b in enumerate(BANDS)},
          "a_link": {"p10": num(np.exp(np.percentile(np.log(a), 10)), "a_p10"), "p90": num(np.exp(np.percentile(np.log(a), 90)), "a_p90"),
                     "per_link": {str(k): float(v) for k, v in zip(ids, a)}, "new_link": 1.0},
          "dynamics_fixed": {k: float(v) for k, v in FIX.items()},
          "lag_min": int(C.d),
          "valid_mm_h": [0.0, 40.0],
          "tabulated_knots_mmh": [float(k) for k in C.kn],
          "intervals": "90 % of 40 cold event-bootstrap draws on TRAIN (S8b)",
          "not_identified": ["phi and tau_s (fixed, the slow store is not identified: S2b, S5a)", "the path term below 1 mm/h (S6)",
                             "the level of F apart from a_link (the ridge prior sets the median a_link to about 1: S6)"],
          "checks": {"fake_data": "9 of 9 numbers identified (S8b)", "TEST": "P2-P5 pass, P1 fails (4 of 8 cells: the WAA level in light rain), S8 and S8b",
                     "TEST_rmse_db": 0.6947},
          "fit_file": rel(fit_file), "fit_sha256": fit_sha,
          "provenance": "fitted on TRAIN in S6 (2026-10-18, before the S8 freeze); adopted by the S8b rule (2026-10-19)"}
    save(json.dumps(eq, indent=1), P.S9_EQUATION)
    eq_sha = hashlib.sha256(P.S9_EQUATION.read_bytes()).hexdigest()
    print(f"equation written, sha256 {eq_sha}", flush=True)

    # ------------------------------------------------------------------------------------------- the unit-test reference (S8 TEST cache)
    le_tr = pickle.loads((P.S8_CACHE / "le_test_train.pkl").read_bytes())[0].reset_index(drop=True)
    Dt = sid.load_data(P.S8_CACHE / "data_test_train.npz", BANDS)
    gt = pd.read_csv(P.S8_G_E)
    gm = {e: 1.0 for e in le_tr.event_id}
    gm.update({r.event_id: (r.g_path if r.path_used else 1.0) for r in gt.itertuples()})
    D1 = rescaled(Dt, le_tr, gm)
    picks = []
    for (lo, hi), _ in LB:
        for e in range(len(le_tr)):
            if lo <= le_tr.L_km[e] < hi:
                o, n = D1["le_roff"][e] + PADM, D1["le_n"][e]
                if (D1["rain"][o:o + n] >= 1.0).sum() >= 30:                 # the config rule: >= 30 min of V1 >= 1 mm/h
                    picks.append(e)
                    break
    Ds, les, Wk, Pk = kernel_parts(D1, le_tr, fc, picks)
    ref = {"what": "3 TEST link-events (S8 cache, V1 with the fixed g_e): the S8b model's 1 min W and path term from the fit's kernel", "lag_min": int(C.d),
           "pad_min": PADM, "tolerance_db": S["unit_test"]["tolerance_db"], "equation_sha256": eq_sha, "link_events": []}
    for k in range(len(les)):
        o, n, o1 = Ds["le_roff"][k], Ds["le_n"][k], Ds["le_boff1"][k]
        cid = les.cml_id[k]
        ref["link_events"].append({"event_id": les.event_id[k], "sublink_id": les.sublink_id[k], "L_km": float(les.L_km[k]), "band": les.band_name[k],
                                   "a_link": float(fc["a"].get(cid, 1.0)), "rain_mm_h": Ds["rain"][o + PADM:o + PADM + n].tolist(),
                                   "rain_before": Ds["rain"][o + PADM - C.d:o + PADM].tolist(), "W": Wk[o1:o1 + n].tolist(), "path": Pk[o1:o1 + n].tolist()})
    save(json.dumps(ref), P.S9_REFERENCE)
    weq = WAAEquation(eq)
    err = max(float(np.max(np.abs(weq.run(r_["rain_mm_h"], r_["L_km"], r_["band"], r_["a_link"], rain_before=r_["rain_before"])["W"] - np.asarray(r_["W"]))))
              for r_ in ref["link_events"])
    print(f"reference written; waa_equation.py reproduces W within {err:.2e} dB", flush=True)

    # ------------------------------------------------------------------------------------------- figure tables
    # example event (TRAIN): the event with the most link-events with V1 >= 10 mm/h; a short and a long 28 GHz link-event
    rmax = np.array([np.nanmax(np.r_[Dt1["r5"][Dt1["le_boff5"][e]:Dt1["le_boff5"][e] + Dt1["le_n"][e] // 5], 0]) for e in range(len(let))])
    heavy = pd.Series(rmax >= 10).groupby(let.event_id).sum()
    ev = heavy.idxmax()
    rsum = np.array([Dt1["rain"][Dt1["le_roff"][e] + PADM:Dt1["le_roff"][e] + PADM + Dt1["le_n"][e]].sum() for e in range(len(let))])
    cand = let[(let.event_id == ev) & (let.band_name == "28 GHz")].copy()
    cand["rsum"] = rsum[cand.index]
    picks2 = [int(cand[cand.L_km < 1].rsum.idxmax()), int(cand[cand.L_km > 4].rsum.idxmax())]
    Ds, les, Wk, Pk = kernel_parts(Dt1, let, fc, picks2)
    rows = []
    for k in range(len(les)):
        lab = "short" if les.L_km[k] < 1 else "long"
        o, n, o1 = Ds["le_roff"][k], Ds["le_n"][k], Ds["le_boff1"][k]
        y = np.where(Ds["w1"][o1:o1 + n] > 0, Ds["y1"][o1:o1 + n] - Ds["le_c"][k], np.nan)
        for t in range(n):
            rows.append({"link": lab, "L_km": float(les.L_km[k]), "event_id": ev, "sublink_id": les.sublink_id[k], "minute": t,
                         "rain_mm_h": float(Ds["rain"][o + PADM + t]), "A_minus_c": float(y[t]), "W": float(Wk[o1 + t]), "path": float(Pk[o1 + t])})
    save(pd.DataFrame(rows).round(4), P.S9_EXAMPLE, index=False)
    # the memory curve: 60 min of 5 mm/h on a 28 GHz link (a_link = 1)
    rain = np.zeros(180)
    rain[10:70] = 5.0
    out = weq.run(rain, 0.0, "28 GHz")
    save(pd.DataFrame({"minute": np.arange(180), "rain_mm_h": rain, "W_fast": out["W_fast"], "W_slow": out["W_slow"], "W": out["W"]}).round(5),
         P.S9_MEMORY, index=False)
    # the closed F and its band (the 40 S8b cold draws)
    R = np.round(np.arange(0, 40.01, 0.25), 2)
    draws = [np.asarray(pickle.loads(f.read_bytes())["z"]) for f in sorted(P.S8B_CACHE.glob("boot_*.pkl"))]
    rows = []
    for i, b in enumerate(BANDS):
        cur = np.array([zz[i] * np.power(R / 5, zz[nb]) for zz in draws])
        lo, hi = np.percentile(cur, [5, 95], axis=0)
        for k, r_ in enumerate(R):
            rows.append({"band": b, "R": r_, "F": z[i] * (r_ / 5) ** z[nb], "lo": lo[k], "hi": hi[k], "P": z[nb + 1 + i] / 5 * r_})
    save(pd.DataFrame(rows).round(5), P.S9_CLOSED_BAND, index=False)
    # the dry-path curves (S7 windows, TRAIN + VAL, the S8b model)
    Rm, T5, col = radar_max(d, sets)
    W1 = {s_: {"closed F": w_minutes(D, le, fc)} for s_, (D, le) in sets.items()}
    win, curves = dry_windows(sets, gmap, Rm, T5, col, W1)
    obs = np.array([c_["obs"] for c_ in curves])
    mod = np.array([c_["closed F"] for c_ in curves])
    nmin = win.dry_minutes.to_numpy()
    rows = []
    for _, nm in LB:
        sel = lclass(win.L_km.to_numpy()) == nm
        for t in range(60):
            m = sel & (nmin > t)
            if m.sum() >= 3:
                rows.append({"length": nm, "minute": t, "windows": int(m.sum()), "measured": float(np.nanmedian(obs[m, t])), "model": float(np.nanmedian(mod[m, t]))})
    save(pd.DataFrame(rows).round(4), P.S9_DRY_CURVES, index=False)
    # the three boxes after the S8 decision (S6b's table was written while M1 was the candidate by the rules; the decision of
    # 2026-10-19 kept the gauge timing and moved M1 to box B), with sizes and sources
    boxes = [
        ("A", "antenna level rises with the rain, F = A_band (R / 5)^b", "about 1 dB at 0.75 mm/h, 2.4-2.8 dB at 5 mm/h", "S6, S8b"),
        ("A", "one shape of F for all bands, one amplitude per band", "A_band 2.36-2.77 dB, intervals overlap", "S6, S8b"),
        ("A", "memory: a fast store wets in 2.6 min and dries in 10.3 min", "72 % of the water", "S2b, S3.4"),
        ("A", "a slow store holds the rest for hours (its numbers fixed)", "28 % of the water, 240 min", "S2b, S5a"),
        ("A", "path loss per km is a straight line, P = c_band R", "0.105-0.194 dB/km per mm/h", "S6, S8b"),
        ("A", "one multiplier per link", "80 % of links 0.77-1.24", "S6, S8b"),
        ("A", "W does not depend on the link length", "rho +0.05, equal dry-path error, TEST P3", "S7, S8"),
        ("B", "dew wets the same radomes without rain", "about 1.7 dB for about 4 h on humid nights", "S7.4"),
        ("B", "the antenna holds less water at the end of rain than the model", "0.2-0.3 dB, 0.5 dB after heavy rain", "S7.3a, S6b"),
        ("B", "a faster effective drying (3.7 min) removes that, but unbinds the path part", "0.17 dB at 5 min on VAL", "S6b"),
        ("B", "drying is faster with wind", "about 42 % per m/s, fails two rules", "S7.5, S6b"),
        ("C", "the WAA level in light rain varies between storms", "TEST 1.35-1.55 dB against 0.83-1.11 dB", "S8 P1"),
        ("C", "why some links wet more than others", "no link property explains it", "S7.1, S7.5"),
        ("C", "late water on long links after rain", "about 0.05 dB at 40 min", "S7.3a"),
        ("C", "the shape of F above 10 mm/h", "0.4 % of the minutes", "S6, S8b"),
        ("C", "the path loss below 1 mm/h", "< 0.1 dB/km, interval at 0", "S6"),
        ("C", "the rest of the error during rain, within storms", "61-64 % of the squared error", "S6b, S8"),
    ]
    save(pd.DataFrame(boxes, columns=["box", "effect", "size", "source"]), P.S9_BOXES, index=False)
    print(f"figure tables written ({len(win)} dry-path windows)", flush=True)

    # ------------------------------------------------------------------------------------------- the report
    errp = max(float(np.max(np.abs(weq.run(r_["rain_mm_h"], r_["L_km"], r_["band"], r_["a_link"], rain_before=r_["rain_before"])["path"]
                                   - np.asarray(r_["path"])))) for r_ in ref["link_events"])
    verdicts = pd.Series([ln.split("|")[4].strip() for ln in P.S9_CHECK.read_text(encoding="utf-8").splitlines()
                          if ln.startswith("| ") and ln.split("|")[1].strip().isdigit()]).value_counts()
    save(report(eq, eq_sha, fit_sha, ref, err, errp, verdicts, len(win)), P.S9_REPORT)
    print(f"report written; done in {(time.time() - clock) / 60:.1f} min")


def report(eq, eq_sha, fit_sha, ref, err, errp, verdicts, n_win):
    """The S9 report: what was written, the corrections of the overview, the open items."""
    A, c = eq["A_band_db"], eq["c_band_db_km_per_mmh"]
    f = lambda q: f"{q['value']:.3g} [{q['lo']:.3g}, {q['hi']:.3g}]"                          # noqa: E731
    les = "\n".join(f"| {r['event_id']} | {r['sublink_id']} | {r['L_km']:.2f} | {r['band']} | {r['a_link']:.3f} | {len(r['W'])} |"
                    for r in ref["link_events"])
    vt = ", ".join(f"{k} {v}" for k, v in verdicts.items())
    return f"""# wet_antenna step S9: formulate Part 1 (the WAA equation)

Step S9 (2026-10-19), settings `wet_antenna.formulate_s9` (written before the run). No new analysis, no refitting, no new reading
of TEST. **Part 1 is closed**: [FINAL.md](../../FINAL.md).

## 1. The official equation

`{rel(P.S9_EQUATION)}`, sha256 `{eq_sha}`. It is the closed-F model of S8b (the S6 plug-back fit "F and P", fitted on TRAIN
before the S8 freeze, `{eq['fit_file']}`, sha256 `{fit_sha}`).

```
A(t) - c = W_f + W_s + L c_band R(t - {eq['lag_min']} min)
dW_f/dt = ({eq['dynamics_fixed']['phi']:.2f} u - W_f) / tau_f,  tau_f = {eq['dynamics_fixed']['tau_f_wet']:.1f} min (wetting), {eq['dynamics_fixed']['tau_f_dry']:.1f} min (drying)
dW_s/dt = ({1 - eq['dynamics_fixed']['phi']:.2f} u - W_s) / {eq['dynamics_fixed']['tau_s']:.0f} min
u = a_link A_band (R / 5)^b
```

| number (90 % interval, 40 cold event-bootstrap draws) | 28 GHz | 32 GHz | 38 GHz |
|---|---|---|---|
| A_band (dB) | {f(A['28 GHz'])} | {f(A['32 GHz'])} | {f(A['38 GHz'])} |
| c_band (dB/km per mm/h) | {f(c['28 GHz'])} | {f(c['32 GHz'])} | {f(c['38 GHz'])} |

b = {f(eq['b'])}. a_link 10 % / 90 %: {f(eq['a_link']['p10'])} / {f(eq['a_link']['p90'])} (per link in the file, 1 for a new link).
Valid {eq['valid_mm_h'][0]:.0f}-{eq['valid_mm_h'][1]:.0f} mm/h. Not identified: {'. '.join(eq['not_identified'])}.

## 2. The code and its unit test

`research/wet_antenna/waa_equation.py` (`WAAEquation`) takes a rain series (mm/h, one value per minute), L, the band and a_link,
and returns W_f, W_s, W and the path term. It steps the stores as the fit does (`system_id.kernel2`) and evaluates F on the fit's
61 knots by default. The reference `{rel(P.S9_REFERENCE)}` holds 3 TEST link-events of the S8 cache (already read once; the first
of each length class with >= 30 min of V1 >= 1 mm/h) with the fit's own 1 min W and path term:

| event | sub-link | L (km) | band | a_link | minutes |
|---|---|---|---|---|---|
{les}

`waa_equation.py` reproduces them within {err:.1e} dB (W) and {errp:.1e} dB (path term), against the tolerance
{ref['tolerance_db']:g} dB. Test: `python -m pytest research/wet_antenna/tests`.

## 3. The overview checked against the repo

`{rel(P.S9_CHECK)}` (written by hand). The overview `docs/WAA_research_overview.pdf` had no LaTeX source in the repo, so
`docs/WAA_research_overview.tex` was rebuilt from the PDF. Verdicts of its {int(verdicts.sum())} statements: {vt}. The repo wins.
The changes, in short:

1. The pipeline averages to 1 and 5 minute values. The lower bands were dropped because all their links are 4 km or longer.
2. TEST is read (S8), with 3355 + 1010 link-events. TRAIN holds 53 events, 51 of them with link-events.
3. c is defined precisely. "No bend of A against L" became "a quadratic term in 4 of 15 rain groups".
4. The length test uses per-link medians in steady rain and a robust line, and its figure now shows the real TRAIN data.
5. The VAL timing is similar (2.0 / 8.8 min), not the same as TRAIN (2.6 / 10.3 min).
6. One shape for the three bands was chosen on unseen events. The training events alone prefer separate shapes.
7. The ITU exponent is 0.86-0.92. a_link is 0.77-1.24 for the final equation.
8. The fake recovery: 62 of 62 within 20 %, 61 within 10 %, and the 9 final numbers within 20 %.
9. How choices were judged: structure on left-out events, the closed form by AICc on the TRAIN points.
10. "A band difference up to 25 % cannot be excluded" was unsupported and was removed. P4 compares each class with the overall
    TEST / VAL ratio. The light-rain levels are 1.35-1.55 dB against 0.83-1.11 dB.
11. The two evaluation errors of S8 are named. The cause of the remaining error is marked as an interpretation.
12. The lab result on runoff now matches what the paper reports, and the title of reference [2] is corrected (Tiede, Chwala,
    Siart and Eibert 2025).
13. The next steps are R1, then Part 2.
14. Figures 1, 2, 4, 5 and 6 are real data. Figure 3 (the two stores) is an illustration from the equation and is captioned so.

## 4. The documents

- `docs/WAA_research_overview.tex` / `.pdf`: the one WAA document. It keeps the structure, wording and look of the first
  overview, with only the corrections above. S9 first wrote a short and a full (school) version in `docs/overview/`; the learners
  then merged them into this one document and deleted the folder (2026-10-19).
- A4, 11 pt, numbered sections and equations, no semicolons, no bold in the text. Built with
  `python docs/build_research_log.py WAA_research_overview` (no LaTeX error, no overfull box).
- Figures `docs/WAA_research_overview_figures/` (6), drawn by its `figs.py` from saved tables only (S3, S6, S7 and the tables of
  this step: `formulate_s9_example_event.csv`, `_memory.csv`, `_closed_f_band.csv`, `_dry_curves.csv` ({n_win} dry-path
  windows)). Black and grey only, with line styles, markers and text labels.

## 5. Open items

- R1, the next step: does the WAA correction improve link rain estimates against the gauges? It needs its own written plan and
  a decision about TEST, which was already read once.
- Then Part 2, the path law f(L, R̄, Z), with this equation removing W first.
- Box C stays open: the WAA level in light rain varies between storms (TEST P1), a_link is not explained by any link property,
  F above 10 mm/h, P below 1 mm/h, and most of the error during rain within storms.
- A small oddity, reported, not changed: the bootstrap interval of the 90th percentile of a_link
  ({eq['a_link']['p90']['lo']:.2f}-{eq['a_link']['p90']['hi']:.2f}) does not contain the fitted value ({eq['a_link']['p90']['value']:.2f}).
  The cold draws spread the link multipliers more than the full fit (S8b table).
- The fit caches (`research/.cache/wet_antenna/`) are still to be published now that the WAA stage is finished (ROADMAP part A).
"""


if __name__ == "__main__":
    main()
