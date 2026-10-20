"""Library of steps S7 and S6b: the frozen S6 models on the TRAIN / VAL V1 datasets, the dry-path windows and the event-bootstrap
helpers. Moved unchanged from validate_s7.py (S6b uses them too); dry_windows takes the models from its W1 argument.
Settings: config `wet_antenna.validate_s7`.
"""
import copy
import math
import pickle

import numpy as np
import pandas as pd

from cmlpipe.api import CMLData
from research.wet_antenna import length_test as LTL
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.final_model import CFG, DATA, M, SHAPES, fit_dict, length_bins, load_sets, predict, s6_models

LB = length_bins()                                                # ((lo, hi), name) of the three length classes
MODELS = ("candidate", "closed F")
PADM = sid.PADM


def lclass(L):
    out = np.empty(len(L), object)
    for (lo, hi), nm in LB:
        out[(L >= lo) & (L < hi)] = nm
    return out


def need(n):
    """Cells needed of n (decision 2026-10-19: at least 4 of 6, k of n = ceil(2 n / 3))."""
    return int(math.ceil(2 * n / 3 - 1e-9))


def ev_draws(n_ev, reps, seed):
    rng = np.random.default_rng(seed)
    return np.stack([np.bincount(rng.integers(0, n_ev, n_ev), minlength=n_ev) for _ in range(reps)]).astype(float)


def boot_items(ev, reps, seed):
    """For items with event codes ev (0..E-1): a list of index arrays, one per event-bootstrap draw."""
    codes, inv = np.unique(ev, return_inverse=True)
    groups = [np.nonzero(inv == i)[0] for i in range(len(codes))]
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(reps):
        pick = rng.integers(0, len(codes), len(codes))
        out.append(np.concatenate([groups[i] for i in pick]) if len(pick) else np.zeros(0, int))
    return out


def ci(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return (float(np.percentile(v, 5)), float(np.percentile(v, 95))) if len(v) else (np.nan, np.nan)


def fmt(v, lo=None, hi=None, nd=2):
    if v is None or not np.isfinite(v):
        return "n/a"
    s = f"{v:.{nd}f}"
    if lo is not None and np.isfinite(lo):
        s += f" [{lo:.{nd}f}, {hi:.{nd}f}]"
    return s


def load():
    d = CMLData.open(P.REPO / "data" / "release" / CFG["release"])
    BANDS = list(sid.merged_bands(d, DATA)[0])
    let = pickle.loads((P.SID3_CACHE / "le_train_all.pkl").read_bytes())[0].reset_index(drop=True)
    lev = pickle.loads((P.SID3_CACHE / "le_val_all.pkl").read_bytes())[0].reset_index(drop=True)
    g4 = pd.read_csv(P.COND4C_G)
    gmap = {s: dict(zip(g4[g4.set == s].event_id, g4[g4.set == s].g_path)) for s in ("TRAIN", "VAL")}
    Dt1 = sid.load_data(P.S5A_CACHE / "data_train_v1.npz", BANDS)
    _, Dv = load_sets(BANDS, let, lev, gmap)
    mods, _ = s6_models(Dt1, let, BANDS)
    free = fit_dict(pickle.loads((P.S5A_CACHE / "fit_M4b-B.pkl").read_bytes()), "M4b-B", Dt1, let)
    sets = {"TRAIN": (Dt1, let), "VAL": (Dv["V1"], lev)}
    return d, BANDS, sets, gmap, mods, free


def predictions(sets, mods, free):
    """5 min predictions of both frozen models (and the S5a free fit) on TRAIN and VAL V1."""
    out = {}
    for s_, (D, le) in sets.items():
        out[s_] = {m: predict(D, le, mods[m]["fz"], mods[m]["fz"]["M"])[0] for m in MODELS}
        out[s_]["S5a free (band shapes)"] = predict(D, le, free)[0]
    return out


def w_minutes(D, le, fz):
    """The model's W alone at 1 min (path term 0, c = 0) for every link-event of a dataset."""
    D0 = dict(D)
    D0["le_c"] = np.zeros(len(le))
    Pm = sid.Problem2(D0, "M4b", fz["M"], agg=1)
    x, _ = LTL.x_for(Pm, fz, LTL.link_ids(D, le))
    x[Pm.names == Pm.pname] = 0.0
    return Pm.predict(x, fz["d"])


def curve_at(C_, z, R):
    """F and P of a frozen model at rain rates R (per band; linear between its knots) and a_link."""
    F, Pv, a = C_.values(z)
    return np.array([np.interp(R, C_.kn, F[b]) for b in range(C_.nb)]), np.array([np.interp(R, C_.kn, Pv[b]) for b in range(C_.nb)]), a


def overlap(a, b):
    return bool(max(a[0], b[0]) <= min(a[1], b[1]))


def radar_max(d, sets, blank_test=True):
    """R_max of the path (every pixel) at the 5 min bins of the TRAIN / VAL windows and the 6 h before each link's first wet bin; all
    other rows NaN (TEST event times are never selected)."""
    st = d.store("radar_along_links_5min")
    T5 = pd.DatetimeIndex(st.time.values)
    col = {int(c): j for j, c in enumerate(st.cml_id.values)}
    need_ = set()
    for D, le in sets.values():
        b0 = T5.searchsorted(pd.to_datetime(le.window_start).to_numpy())
        fw = T5.searchsorted(pd.to_datetime(le.first_wet).to_numpy())
        for a, n, f in zip(b0, le.n_bins.to_numpy(), fw):
            need_.update(range(a, a + n))
            need_.update(range(max(f - 72, 0), f))
    ev = CMLData.open(P.REPO / "data" / "release" / CFG["release"]).events()
    test = np.zeros(len(T5), bool)
    for e in (ev[ev.split_event == "test"].itertuples() if blank_test else []):     # step S8 reads the TEST windows
        test[T5.searchsorted(e.start):T5.searchsorted(e.end)] = True
    idx = np.array(sorted(i for i in need_ if not test[i]))
    Rm = np.full((len(T5), len(col)), np.nan, np.float32)
    Rm[idx] = st.R_max.isel(time=idx).load().to_numpy()
    return Rm, T5, col


def dry_windows(sets, gmap, Rm, T5, col, W1):
    rows, curves = [], []
    for s_, (D, le) in sets.items():
        b0s = T5.searchsorted(pd.to_datetime(le.window_start).to_numpy())
        for e, r in enumerate(le.itertuples()):
            n5, o1, roff = D["le_n"][e] // 5, D["le_boff1"][e], D["le_roff"][e]
            n = D["le_n"][e]
            rmax = Rm[b0s[e]:b0s[e] + n5, col[r.cml_id]].astype(float)
            g = float(gmap[s_].get(r.event_id, 1.0))
            gau = D["rain_gauge"][roff + PADM:roff + PADM + n]
            g5 = gau.reshape(n5, 5)
            gauge_ok = np.where(np.isfinite(g5), g5 == 0, True).all(axis=1)
            dry = np.isfinite(rmax) & (rmax <= 0.05) & (rmax / g <= 0.05) & gauge_ok
            rain1 = D["rain"][roff + PADM:roff + PADM + n]
            for j in range(1, n5 - 1):
                if not (dry[j] and dry[j + 1] and not dry[j - 1]):
                    continue
                m0 = j * 5
                before = rain1[max(m0 - 60, 0):m0]
                if (before >= 2.0).sum() < 10:
                    continue
                k = j
                while k < n5 and dry[k]:
                    k += 1
                nmin = min(5 * (k - j), 60)
                yo = np.full(60, np.nan)
                seg = slice(o1 + m0, o1 + m0 + nmin)
                yv = D["y1"][seg] - D["le_c"][e]
                yv[D["w1"][seg] == 0] = np.nan
                yo[:nmin] = yv
                wm = {m: np.r_[W1[s_][m][seg], np.full(60 - nmin, np.nan)] for m in W1[s_]}
                rows.append({"set": s_, "le_idx": e, "event_id": r.event_id, "cml_id": r.cml_id, "sublink_id": r.sublink_id, "L_km": r.L_km,
                             "band": r.band_name, "length": lclass(np.array([r.L_km]))[0], "t0": pd.Timestamp(T5[b0s[e] + j]), "dry_minutes": nmin,
                             "rain_ge2_min_before": int((before >= 2).sum()), "heavy": bool((before >= 10).sum() >= 10), "g_e": g,
                             "W0_obs": float(np.nanmean(yo[:2])), "W0_candidate": float(wm["candidate"][0]) if "candidate" in wm else np.nan})
                curves.append({"obs": yo, **wm})
    return pd.DataFrame(rows), curves


def fit_tau(curve):
    t = np.arange(31)
    m = np.asarray(curve[:31], float)
    ok = np.isfinite(m) & (m > 0.05)
    if ok.sum() < 5:
        return np.nan
    slope = np.polyfit(t[ok], np.log(m[ok]), 1)[0]
    return float(-1 / slope) if slope < 0 else np.nan


def norm_median(arr, nmin, w0):
    """Median over windows of W(t) / W0 for t = 0-59 (windows still dry at t)."""
    out = np.full(60, np.nan)
    rel_ = arr / w0[:, None]
    for t in range(60):
        v = rel_[nmin > t, t]
        v = v[np.isfinite(v)]
        if len(v) >= 3:
            out[t] = np.median(v)
    return out


def error_budget(D, le, pred):
    """The squared error of 5 min A split into dry minutes, during rain (between / within events) and the first hour after the rain
    (moved unchanged from return_s6b.py; S8 uses it too)."""
    w = D["w5"] > 0
    r5 = np.nan_to_num(D["r5"], nan=-1)
    aft = D["after5"]
    res = D["y5"] - pred
    rain = w & (r5 >= 0.1)
    after = w & (r5 < 0.1) & (aft >= 0) & (aft < 60)
    dry = w & ~rain & ~after
    tot = float((res[w] ** 2).sum())
    rows = []
    for nm, g in (("dry minutes (noise floor)", dry), ("during rain (V1 >= 0.1 mm/h)", rain), ("after rain (dry path, first 60 min)", after)):
        sse = float((res[g] ** 2).sum())
        rows.append({"part": nm, "bins": int(g.sum()), "rmse_db": float(np.sqrt(sse / g.sum())), "share_of_squared_error": sse / tot})
    ev = D["le_event"][D["bin_le"]]
    n_e = np.bincount(ev[rain], minlength=ev.max() + 1).astype(float)
    s_e = np.bincount(ev[rain], weights=res[rain], minlength=ev.max() + 1)
    between = float(np.sum(np.where(n_e > 0, s_e ** 2 / np.maximum(n_e, 1), 0.0)))
    sse_r = float((res[rain] ** 2).sum())
    rows.append({"part": "  during rain: between events (the event's mean residual; the rain reference)", "bins": int(rain.sum()),
                 "rmse_db": float(np.sqrt(between / rain.sum())), "share_of_squared_error": between / tot})
    rows.append({"part": "  during rain: within events", "bins": int(rain.sum()), "rmse_db": float(np.sqrt((sse_r - between) / rain.sum())),
                 "share_of_squared_error": (sse_r - between) / tot})
    rows.append({"part": "all", "bins": int(w.sum()), "rmse_db": float(np.sqrt(tot / w.sum())), "share_of_squared_error": 1.0})
    return pd.DataFrame(rows)


# =================================================================================================== the S8 comparison models (moved unchanged from test_s8.py; S8b uses them too)
def step_g(R):
    """B1: 0 below 0.1 mm/h, 1 from 0.1001 mm/h on (a near-step on the knots)."""
    return np.clip((np.asarray(R, float) - 0.1) / 1e-4, 0.0, 1.0)


def comp_setup(name, nb):
    """(spec, M of the model, the static problem name) of a comparison model."""
    lin = {**SHAPES["linear"]}
    Mx = copy.deepcopy(M)
    if name == "B1":
        Mx["knots_mmh"] = [0.0, 0.1, 0.1001] + [k for k in M["knots_mmh"] if k > 0.1001]
        F = ("form", {"g": step_g, "lo": [], "hi": [], "starts": [[]], "text": "w0 (R > 0.1 mm/h)", "names": []}, True)
    elif name == "B0":
        F = ("fixed", np.zeros((nb, len(Mx["knots_mmh"]))))
    else:
        F = "shared"
    return {"F": F, "P": ("form", lin, True)}, Mx
