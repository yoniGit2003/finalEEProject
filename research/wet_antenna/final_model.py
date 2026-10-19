"""Library of steps S5a and S6: the model setup of S5a, a fitted model on another dataset (errors by length class, the
length test through the model), and the closed forms of F(R) and P(R) with the plug-back model.

Moved unchanged from final_s5a.py. Settings: config `wet_antenna.final_s5a` (the closed-form candidates, the AICc rule,
the length-test share) and the S2 / S2b / S3 sections it builds on.
"""
import copy
import json
import pickle

import numpy as np
import pandas as pd
from scipy.optimize import least_squares, minimize

from research.common.config import load_config
from research.wet_antenna import length_test as LTL
from research.wet_antenna import paths as P
from research.wet_antenna import system_id as sid
from research.wet_antenna.event_factors import rescaled

CFG = load_config()
WA = CFG["wet_antenna"]
S2, S2B, S3, S3B, C = WA["system_id_s2"], WA["system_id_s2b"], WA["system_id_s3"], WA["system_id_s3b"], WA["final_s5a"]
DATA, CHK, FK = S2["data"], S2["checks"], S2["fake"]
M = copy.deepcopy(S2["model"])                                # the S2b model settings (as S3b), plus the M4c bounds
M["increments"], M["maxiter"] = S2B["model"]["increments"], dict(S2B["model"]["maxiter"])
M["bounds"] = {**M["bounds"], "increment": S2B["model"]["increment_bounds"], "tau_s_min": S2B["model"]["tau_s_m4b_min"],
               "tau_s_split": C["tau_s_split_bounds"]}
J = {**S2["judge"], "order": S2B["judge"]["order"]}
LT = S3["length_test"]
KNOTS = np.asarray(M["knots_mmh"], float)
KEY_R = S2["length_rule"]["key_rain_mmh"]
FITS = C["fits"]
CF = C["closed_forms"]
FIX = dict(C["fits"]["M4b-B"]["fixed"])          # tau_f_wet 2.6, tau_f_dry 10.3, phi 0.72, tau_s 240
INC_HI = S2B["model"]["increment_bounds"][1]


# =================================================================================================== the model and its starts
def plan_start(Pm, key, itu_inc, slow=None):
    """A start of the S2b start plan ("g1+P0" = generic start g1, P = 0 x the ITU size: a STARTING value only)."""
    g, p = key.split("+")
    s = dict(M["starts"][g])
    if slow:
        s.update(slow)
    return Pm.start(s, Pinc=S2B["model"]["p_starts"][p] * np.asarray(itu_inc))


def problem(D, name):
    Pm = sid.Problem2(D, FITS[name]["model"], M)
    sid.fix_params(Pm, FITS[name]["fixed"])
    return Pm


# =================================================================================================== a fitted model on a dataset
def fit_dict(res, name, Dt, let):
    return {"x": res["all"]["x"], "d": int(res["all"]["d"]), "a": dict(zip(LTL.link_ids(Dt, let), res["a_link"])), "model": FITS[name]["model"],
            "params": res["params"]}


def predict(D, le, f, Mx=None):
    Pm = sid.Problem2(D, f["model"], Mx or M)
    x, n_new = LTL.x_for(Pm, f, LTL.link_ids(D, le))
    return Pm.predict(x, f["d"]), n_new


def length_bins():
    out = []
    for lo, hi in C["length_bins_km"]:
        out.append(((lo, hi), f"L < {hi:g} km" if lo == 0 else (f"L > {lo:g} km" if hi > 50 else f"{lo:g}-{hi:g} km")))
    return out


def val_table(Dv, lev, fz, inp):
    """VAL errors per fit and length class (and all lengths) for one input."""
    rows = []
    sets = [("all lengths", np.ones(len(lev), bool))] + [(nm, ((lev.L_km >= lo) & (lev.L_km < hi)).to_numpy()) for (lo, hi), nm in length_bins()]
    for k, f in fz.items():
        pr, n_new = predict(Dv, lev, f, f.get("M"))
        for nm, keep in sets:
            D_ = sid.subset(Dv, keep) if nm != "all lengths" else Dv
            p_ = pr if nm == "all lengths" else predict(D_, lev[keep].reset_index(drop=True), f, f.get("M"))[0]
            rows.append({"input": inp, "fit": k, "VAL set": nm, "link_events": int(keep.sum()), "links_a_link_1": n_new, **LTL.val_metrics(D_, p_, J)})
    return pd.DataFrame(rows)


def through(D, le, fz, set_, inp):
    """The length test through the model (S3b check 1): observed I, s with the event bootstrap, each fit's own line on the
    same bins; per cell both inside? (eligible = >= min_short_links links under 1 km)."""
    preds = {k: predict(D, le, f, f.get("M"))[0] for k, f in fz.items()}
    t = LTL.steady_table(D, le, LT, preds)
    n_ev = int(D["le_event"].max()) + 1
    obs = LTL.length_fits(t, n_ev, set_, LT)
    rows = []
    for k in fz:
        mod = LTL.length_fits(t, n_ev, set_, LT, ycol=f"pred_{k}", boot=False).set_index(["band", "class"])
        for r in obs.dropna(subset=["I"]).to_dict("records"):
            key = (r["band"], r["class"])
            if key not in mod.index or not np.isfinite(mod.loc[key, "I"]):
                continue
            q = mod.loc[key]
            rows.append({"input": inp, "set": set_, "fit": k, "band": r["band"], "class": r["class"], "R_cell": r["R_cell"], "links_under_1km": r["links_under_1km"],
                         "I_obs": r["I"], "I_lo": r["I_lo"], "I_hi": r["I_hi"], "I_model": q.I, "s_obs": r["s"], "s_lo": r["s_lo"], "s_hi": r["s_hi"], "s_model": q.s,
                         "eligible": bool(r["links_under_1km"] >= LT["verdict"]["min_short_links"]),
                         "I_inside": bool(r["I_lo"] <= q.I <= r["I_hi"]), "s_inside": bool(r["s_lo"] <= q.s <= r["s_hi"])})
    return pd.DataFrame(rows)


def through_verdict(cells):
    rows = []
    for (inp, s, f), g in cells.groupby(["input", "set", "fit"], sort=False):
        e = g[g.eligible]
        both = e.I_inside & e.s_inside
        share = float(both.mean()) if len(e) else np.nan
        rows.append({"input": inp, "set": s, "fit": f, "eligible_cells": len(e), "I_inside": int(e.I_inside.sum()), "s_inside": int(e.s_inside.sum()),
                     "both_inside": int(both.sum()), "share": share, "I_ratio_median": float(np.median(e.I_obs / e.I_model)) if len(e) else np.nan,
                     "s_ratio_median": float(np.median(e.s_obs / e.s_model)) if len(e) else np.nan,
                     "passes": bool(len(e) and share >= C["length_share"])})
    return pd.DataFrame(rows)


def load_sets(BANDS, let, lev, gmap):
    Dt = sid.load_data(P.SID3_CACHE / "data_train_all.npz", BANDS)
    Dv = sid.load_data(P.SID3_CACHE / "data_val_all.npz", BANDS)
    return Dt, {"V0": Dv, "V1": rescaled(Dv, lev, gmap["VAL"])}


# =================================================================================================== closed forms
FORMS = {   # name: (function of R and the parameters, parameter names, lower bounds, upper bounds, start values to try)
    "exp_sat": (lambda R, a, r0: a * (1 - np.exp(-R / r0)), ["F_max", "R0"], [1e-4, 1e-3], [100, 200], [[3, 1], [3, 3], [2, 0.3]]),
    "hyperbolic": (lambda R, a, r0: a * R / (R + r0), ["F_max", "R0"], [1e-4, 1e-3], [100, 200], [[3, 1], [4, 3], [2, 0.3]]),
    "tanh_sat": (lambda R, a, r0: a * np.tanh(R / r0), ["F_max", "R0"], [1e-4, 1e-3], [100, 200], [[3, 1], [3, 3], [2, 0.3]]),
    "capped_power": (lambda R, c, b, a: np.minimum(c * np.power(R, b), a), ["c", "b", "F_max"], [1e-4, 0.2, 1e-4], [100, 3, 100],
                     [[1.5, 0.5, 3], [1, 0.3, 3.5], [2, 1, 2.5]]),
    "power": (lambda R, c, b: c * np.power(R, b), ["c", "b"], [1e-5, 0.2], [100, 3], [[0.1, 1], [0.2, 0.8], [0.05, 1.3]]),
    "linear": (lambda R, c: c * R, ["c"], [1e-5], [100], [[0.1], [0.3]]),
    "quadratic": (lambda R, c, d: c * R + d * R ** 2, ["c", "d"], [1e-5, 0.0], [100, 10], [[0.1, 0.001], [0.2, 0.0]]),
    "shifted_power": (lambda R, c, r1, b: c * (np.power(R + r1, b) - np.power(r1, b)), ["c", "R1", "b"], [1e-5, 1e-3, 0.2], [100, 200, 3],
                      [[0.1, 1, 1.2], [0.2, 0.5, 0.8], [0.05, 3, 1.5]]),
}
FORM_TEXT = {"exp_sat": "F_max (1 - exp(-R / R0))", "hyperbolic": "F_max R / (R + R0)", "tanh_sat": "F_max tanh(R / R0)", "capped_power": "min(c R^b, F_max)",
             "power": "c R^b", "linear": "c R", "quadratic": "c R + d R^2", "shifted_power": "c ((R + R1)^b - R1^b)"}


def wls(form, R, y, se, starts=None):
    fn, names, lo, hi, st = FORMS[form]
    best = None
    for s0 in (starts or st):
        s0 = np.clip(s0, np.array(lo) * 1.0001 + 1e-9, np.array(hi) * 0.9999)
        try:
            r = least_squares(lambda th: (fn(R, *th) - y) / se, s0, bounds=(lo, hi), x_scale="jac", max_nfev=4000)
        except ValueError:
            continue
        if best is None or r.cost < best.cost:
            best = r
    return best


def aicc(wrss, n, k):
    if n - k - 1 <= 0 or wrss <= 0:
        return np.inf
    return n * np.log(wrss / n) + 2 * k + 2 * k * (k + 1) / (n - k - 1)


def pooled_fit(form, Rb, yb, seb, fb, variant):
    """Over the bands: (i) "one curve" for all bands; (ii) "amp linear in f": the first (amplitude) parameter
    A30 (1 + kappa (f - 30) / 10), the shape parameters shared."""
    fn, names, lo, hi, st = FORMS[form]
    R, y, se, f = np.concatenate(Rb), np.concatenate(yb), np.concatenate(seb), np.concatenate(fb)
    if variant == "one curve":
        r = wls(form, R, y, se)
        return r, names, len(names), y, se
    lo2, hi2 = list(lo) + [-3.0], list(hi) + [3.0]

    def res(th):
        amp = th[0] * (1 + th[-1] * (f - 30) / 10)
        return (fn(R, amp, *th[1:-1]) - y) / se
    best = None
    for s0 in st:
        s0 = np.r_[np.clip(s0, np.array(lo) * 1.0001 + 1e-9, np.array(hi) * 0.9999), 0.0]
        try:
            r = least_squares(res, s0, bounds=(lo2, hi2), x_scale="jac", max_nfev=4000)
        except ValueError:
            continue
        if best is None or r.cost < best.cost:
            best = r
    return best, names + ["kappa"], len(names) + 1, y, se


def choose(scores):
    """Lowest AICc; a candidate within aicc_simpler_within of it with fewer parameters is preferred."""
    s = scores[np.isfinite(scores.aicc)].sort_values("aicc")
    if not len(s):
        return None
    best = s.iloc[0]
    near = s[(s.aicc <= best.aicc + CF["aicc_simpler_within"]) & (s.k < best.k)]
    return (near.sort_values(["k", "aicc"]).iloc[0] if len(near) else best).form


class ClosedForm:
    """The kept model with F_band and P_band given by closed forms on a dense knot grid; the closed-form parameters and
    a_link are fitted, the time constants, phi and the lag are FIXED at the kept model's values."""

    def __init__(self, D, model, fixed_x_src, P_src, forms, d):
        kd = CF["plug_back_knots"]
        self.M = copy.deepcopy(M)
        self.M["knots_mmh"] = [0.0] + list(np.geomspace(kd["lo"], kd["hi"], kd["n"]))
        self.knots = np.asarray(self.M["knots_mmh"])
        self.Pm = sid.Problem2(D, model, self.M)
        self.x0 = sid.x_into(self.Pm, P_src, fixed_x_src, curves=False)   # dyn and a_link from the kept fit (curves set by build)
        self.forms, self.d = forms, d                        # forms = {(kind, band index): (form name, theta)}
        self.keys = list(forms)
        self.sizes = [len(FORMS[forms[k][0]][1]) for k in self.keys]
        self.lo = np.concatenate([FORMS[forms[k][0]][2] for k in self.keys])
        self.hi = np.concatenate([FORMS[forms[k][0]][3] for k in self.keys])
        self.ia = self.Pm.names == "loga"

    def incs(self, kind, b, th):
        v = FORMS[self.forms[(kind, b)][0]][0](self.knots, *th)
        return np.diff(v)

    def build(self, z):
        x = self.x0.copy()
        th_all, la = z[:len(self.lo)], z[len(self.lo):]
        Finc, Pinc = x[self.Pm.names == self.Pm.fname].reshape(self.Pm.nb, -1), x[self.Pm.names == self.Pm.pname].reshape(self.Pm.nb, -1)
        o = 0
        for (kind, b), n in zip(self.keys, self.sizes):
            (Finc if kind == "F" else Pinc)[b] = self.incs(kind, b, th_all[o:o + n])
            o += n
        x[self.Pm.names == self.Pm.fname], x[self.Pm.names == self.Pm.pname] = Finc.ravel(), Pinc.ravel()
        x[self.ia] = la
        return x

    def fun(self, z):
        x = self.build(z)
        loss, g = self.Pm.run(x, self.d, np.ones(self.Pm.ne))
        gF, gP = g[self.Pm.names == self.Pm.fname].reshape(self.Pm.nb, -1), g[self.Pm.names == self.Pm.pname].reshape(self.Pm.nb, -1)
        gz = np.zeros(len(z))
        o = 0
        th_all = z[:len(self.lo)]
        for (kind, b), n in zip(self.keys, self.sizes):
            th = th_all[o:o + n]
            gi = (gF if kind == "F" else gP)[b]
            for j in range(n):                                # the chain rule through the formula (central differences)
                h = 1e-6 * max(abs(th[j]), 1e-3)
                tp, tm = th.copy(), th.copy()
                tp[j] += h
                tm[j] -= h
                gz[o + j] = gi @ ((self.incs(kind, b, tp) - self.incs(kind, b, tm)) / (2 * h))
            o += n
        gz[len(self.lo):] = g[self.ia]
        return loss, gz

    def fit(self, maxiter):
        th0 = np.concatenate([self.forms[k][1] for k in self.keys])
        z0 = np.r_[np.clip(th0, self.lo, self.hi), self.x0[self.ia]]
        bounds = list(zip(self.lo, self.hi)) + [tuple(b) for b, a in zip(self.Pm.bounds, self.ia) if a]
        r = minimize(self.fun, z0, jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": maxiter, "maxfun": int(maxiter * 1.5), "ftol": 1e-10, "gtol": 1e-6})
        self.z, self.result = r.x, r
        return r


# =================================================================================================== parameters, knots, closed-form choice
def param_table(res, boots, s3b_params):
    rows = []
    for k, r in res.items():
        b = boots.get(k, {})
        for name, v in r["params"].items():
            arr = np.array([dr["params"].get(name, np.nan) for dr in b.values()], float)
            lo, hi = np.nanpercentile(arr, J["interval"]) if np.isfinite(arr).any() else (np.nan, np.nan)
            rows.append({"fit": k, "parameter": name, "value": v, "lo": lo, "hi": hi, "se": float(np.nanstd(arr, ddof=1)) if np.isfinite(arr).sum() > 1 else np.nan,
                         "draws": int(np.isfinite(arr).sum()), "S3b": s3b_params.get(name, np.nan) if k in ("M4b-B", "M4c") else np.nan})
    return pd.DataFrame(rows)


def knot_table(params, kept, BANDS, Dt1, let):
    """F and P of the kept model at the knots with their bootstrap SE (and the floor of the weights)."""
    rain_max = np.array([Dt1["rain"][Dt1["le_roff"][e] + sid.PADM:Dt1["le_roff"][e] + sid.PADM + Dt1["le_n"][e]].max() for e in range(len(let))])
    n_above = int((rain_max > 20).sum())
    knots = [k for k in CF["knots_mmh"] if not (k == 40 and n_above < CF["drop_40_if_fewer_than"])]
    pv = params[params.fit == kept].set_index("parameter")
    fl = CF["se_floor"]
    freq = let.groupby("band_name").frequency_ghz.median()
    rows = []
    for kind in ("F", "P"):
        for b in BANDS:
            for k in knots:
                r = pv.loc[f"{kind}[{b}]({k:g})"]
                floor = max(fl["rel"] * abs(r.value), fl["F_db"] if kind == "F" else fl["P_db_km"])
                rows.append({"kind": kind, "band": b, "f_ghz": float(freq[b]), "R": k, "value": r.value, "lo": r.lo, "hi": r.hi, "se_boot": r.se,
                             "se_used": max(r.se if np.isfinite(r.se) else 0.0, floor)})
    return pd.DataFrame(rows), n_above, knots


def fit_forms(kt):
    """Every candidate per band and pooled (two variants), with AICc; the chosen form per band and pooled."""
    rows, chosen = [], {}
    for kind, names in (("F", CF["F"]), ("P", CF["P"])):
        sub = kt[kt.kind == kind]
        for b, gb in sub.groupby("band", sort=False):
            R, y, se = gb.R.to_numpy(float), gb.value.to_numpy(float), gb.se_used.to_numpy(float)
            for form in names:
                r = wls(form, R, y, se)
                k = len(FORMS[form][1])
                wrss = float(2 * r.cost) if r is not None else np.inf
                rows.append({"kind": kind, "scope": b, "form": form, "formula": FORM_TEXT[form], "k": k, "n": len(R), "wrss": wrss, "aicc": aicc(wrss, len(R), k),
                             "params": json.dumps(dict(zip(FORMS[form][1], np.round(r.x, 5).tolist()))) if r is not None else "",
                             "max_rel_dev": float(np.max(np.abs(FORMS[form][0](R, *r.x) - y) / np.maximum(np.abs(y), 1e-9))) if r is not None else np.nan})
            sc = pd.DataFrame([q for q in rows if q["kind"] == kind and q["scope"] == b])
            chosen[(kind, b)] = choose(sc)
        groups = [g_ for _, g_ in sub.groupby("band", sort=False)]
        for variant in ("one curve", "amp linear in f"):
            for form in names:
                r, pn, k, y, se = pooled_fit(form, [g_.R.to_numpy(float) for g_ in groups], [g_.value.to_numpy(float) for g_ in groups],
                                             [g_.se_used.to_numpy(float) for g_ in groups], [g_.f_ghz.to_numpy(float) for g_ in groups], variant)
                wrss = float(2 * r.cost) if r is not None else np.inf
                rows.append({"kind": kind, "scope": f"pooled: {variant}", "form": form, "formula": FORM_TEXT[form], "k": k, "n": len(y), "wrss": wrss,
                             "aicc": aicc(wrss, len(y), k), "params": json.dumps(dict(zip(pn, np.round(r.x, 5).tolist()))) if r is not None else "", "max_rel_dev": np.nan})
            sc = pd.DataFrame([q for q in rows if q["kind"] == kind and q["scope"] == f"pooled: {variant}"])
            chosen[(kind, f"pooled: {variant}")] = choose(sc)
    t = pd.DataFrame(rows)
    t["chosen"] = [chosen.get((r.kind, r.scope)) == r.form for r in t.itertuples()]
    t["delta_aicc"] = t.aicc - t.groupby(["kind", "scope"]).aicc.transform("min")
    return t, chosen


def form_params(forms, kind, scope, form):
    """The fitted parameters of one candidate (in the order of FORMS)."""
    row = forms[(forms.kind == kind) & (forms.scope == scope) & (forms.form == form)].iloc[0]
    return np.array(list(json.loads(row.params).values()), float)


def form_intervals(kt, chosen, boots_kept, BANDS):
    """The chosen per-band form fitted to the knot curves of every cold bootstrap draw (same weights)."""
    rows = []
    for (kind, b), form in chosen.items():
        if b not in BANDS or form is None:
            continue
        gb = kt[(kt.kind == kind) & (kt.band == b)]
        R, se = gb.R.to_numpy(float), gb.se_used.to_numpy(float)
        main = wls(form, R, gb.value.to_numpy(float), se)
        draws = []
        for dr in boots_kept.values():
            y = np.array([dr["params"].get(f"{kind}[{b}]({k:g})", np.nan) for k in R])
            if np.isfinite(y).all():
                r = wls(form, R, y, se, starts=[main.x])
                if r is not None:
                    draws.append(r.x)
        draws = np.array(draws)
        for i, nm in enumerate(FORMS[form][1]):
            lo, hi = np.nanpercentile(draws[:, i], J["interval"]) if len(draws) else (np.nan, np.nan)
            rows.append({"kind": kind, "band": b, "form": form, "parameter": nm, "value": main.x[i], "lo": lo, "hi": hi, "draws": len(draws)})
    return pd.DataFrame(rows)


# =================================================================================================== step S6: the parametrisations (moved unchanged from tighten_s6.py)
class Curves:
    """M4b with the dynamics and the lag FIXED; F and P each given by a parametrisation; a_link free (the S3b ridge).

    spec[kind] (kind "F" or "P") is one of
      "free"                  the knot increments per band (as Problem2)
      "shared"                increments = amp_band * s, s = [softmax(0, v) (sums to 1 up to 5 mm/h), w (above 5, >= 0)]
      ("form", f, shared)     values = amp_band * g(knots; theta), g a closed-form shape with g(0) = 0, g(5) = 1; theta
                              shared over the bands or one set per band
      ("fixed", values)       the values at the knots (nb x nk+1), not fitted
    pen = profile penalties [(q, v, lam)], q(Finc, Pinc, la) -> (value, dF, dP, dla)."""

    def __init__(self, D, spec, x_ref, d, Mx=None, fixed=None, wle=None, Pm=None, free=()):
        self.M = Mx or M
        self.Pm = Pm if Pm is not None else sid.Problem2(D, "M4b", self.M)   # step S6b: a problem with a wind-dependent drying
        self.kn = np.asarray(self.M["knots_mmh"], float)
        self.nb, self.nk, self.nl = self.Pm.nb, self.Pm.nk, self.Pm.nl
        self.k5 = int(np.argmin(np.abs(self.kn[1:] - 5.0))) + 1
        self.x0 = np.asarray(x_ref, float).copy()
        for p, v in {**FIX, **(fixed or {})}.items():
            self.x0[self.Pm.names == p] = np.log(v / (1 - v)) if p == "phi" else np.log(v)
        self.d, self.spec, self.pen = int(d), spec, []
        self.wle = np.ones(self.Pm.ne) if wle is None else wle
        self.ia = self.Pm.names == "loga"
        self.blocks, lo, hi = [], [], []
        for kind in ("F", "P"):
            sp = spec[kind]
            if sp == "free":
                n, l_, h_ = self.nb * self.nk, [0.0] * (self.nb * self.nk), [INC_HI] * (self.nb * self.nk)
            elif sp == "shared":
                n = self.nb + (self.k5 - 1) + (self.nk - self.k5)
                l_ = [0.0] * self.nb + [-12.0] * (self.k5 - 1) + [0.0] * (self.nk - self.k5)
                h_ = [INC_HI] * self.nb + [12.0] * (self.k5 - 1) + [INC_HI] * (self.nk - self.k5)
            elif sp[0] == "form":
                f, shared = sp[1], sp[2]
                nt = len(f["lo"]) * (1 if shared else self.nb)
                n, l_, h_ = self.nb + nt, [0.0] * self.nb + list(f["lo"]) * (nt // max(len(f["lo"]), 1)), [100.0] * self.nb + list(f["hi"]) * (nt // max(len(f["lo"]), 1))
            else:
                n, l_, h_ = 0, [], []
            self.blocks.append((kind, len(lo), n))
            lo += l_
            hi += h_
        self.nc = len(lo)
        lo += [b[0] for b, a in zip(self.Pm.bounds, self.ia) if a]
        hi += [b[1] for b, a in zip(self.Pm.bounds, self.ia) if a]
        self.free, self.nf = list(free), len(free)             # step S6b: dynamic parameters fitted too (after a_link in z)
        for p in self.free:
            bd = self.Pm.bounds[int(np.nonzero(self.Pm.names == p)[0][0])]
            lo.append(bd[0])
            hi.append(bd[1])
        self.lo, self.hi = np.array(lo), np.array(hi)

    # --------------------------------------------------------------------------------------------- z <-> curves
    def incs(self, kind, zk):
        """(increments nb x nk, a function grad(g_inc) -> g_z)."""
        sp, nb, nk, k5 = self.spec[kind], self.nb, self.nk, self.k5
        if sp == "free":
            return zk.reshape(nb, nk), lambda g: g.ravel()
        if sp == "shared":
            amp, v, w = zk[:nb], zk[nb:nb + k5 - 1], zk[nb + k5 - 1:]
            e = np.exp(np.r_[0.0, v] - np.max(np.r_[0.0, v]))
            u = e / e.sum()
            s = np.r_[u, w]

            def grad(g):
                gamp = g @ s
                gs = amp @ g
                gu = gs[:k5]
                gv = u[1:] * (gu[1:] - gu @ u)
                return np.r_[gamp, gv, gs[k5:]]
            return np.outer(amp, s), grad
        if sp[0] == "fixed":
            return np.diff(np.asarray(sp[1], float), axis=1), lambda g: np.zeros(0)
        f, shared = sp[1], sp[2]
        amp, th = zk[:nb], zk[nb:]
        nth = len(f["lo"])
        ths = [th[:nth] if shared else th[b * nth:(b + 1) * nth] for b in range(nb)]
        gf = [shape_g(f, b) for b in range(nb)]
        G = np.array([np.diff(gf[b](self.kn, *t)) for b, t in enumerate(ths)])
        inc = amp[:, None] * G

        def grad(g):
            gamp = (g * G).sum(axis=1)
            gth = np.zeros(len(th))
            for b in range(nb):
                t = ths[b]
                for j in range(nth):
                    h = 1e-6 * max(abs(t[j]), 1e-3)
                    tp, tm = t.copy(), t.copy()
                    tp[j] += h
                    tm[j] -= h
                    dG = (np.diff(gf[b](self.kn, *tp)) - np.diff(gf[b](self.kn, *tm))) / (2 * h)
                    gth[(0 if shared else b * nth) + j] += amp[b] * (g[b] @ dG)
            return np.r_[gamp, gth]
        return inc, grad

    def build(self, z):
        x = self.x0.copy()
        out = {}
        for kind, o, n in self.blocks:
            inc, _ = self.incs(kind, z[o:o + n])
            out[kind] = inc
        x[self.Pm.names == self.Pm.fname] = out["F"].ravel()
        x[self.Pm.names == self.Pm.pname] = out["P"].ravel()
        x[self.ia] = z[self.nc:self.nc + self.nl]
        for i, p in enumerate(self.free):
            x[self.Pm.names == p] = z[self.nc + self.nl + i]
        return x

    def fun(self, z):
        x = self.build(z)
        loss, g = self.Pm.run(x, self.d, self.wle)
        gF = g[self.Pm.names == self.Pm.fname].reshape(self.nb, self.nk)
        gP = g[self.Pm.names == self.Pm.pname].reshape(self.nb, self.nk)
        gla = g[self.ia].copy()
        Finc, Pinc, la = x[self.Pm.names == self.Pm.fname].reshape(self.nb, self.nk), x[self.Pm.names == self.Pm.pname].reshape(self.nb, self.nk), z[self.nc:self.nc + self.nl]
        for q, v, lam in self.pen:
            val, dF, dP, dla = q(Finc, Pinc, la)
            loss += lam * (val - v) ** 2
            gF, gP, gla = gF + 2 * lam * (val - v) * dF, gP + 2 * lam * (val - v) * dP, gla + 2 * lam * (val - v) * dla
        gz = np.zeros(len(z))
        for kind, o, n in self.blocks:
            _, gr = self.incs(kind, z[o:o + n])
            gz[o:o + n] = gr(gF if kind == "F" else gP)
        gz[self.nc:self.nc + self.nl] = gla
        for i, p in enumerate(self.free):
            gz[self.nc + self.nl + i] = g[self.Pm.names == p].sum()
        return loss, gz

    def fit(self, z0, maxiter):
        z0 = np.clip(z0, self.lo, self.hi)
        r = minimize(self.fun, z0, jac=True, method="L-BFGS-B", bounds=list(zip(self.lo, self.hi)),
                     options={"maxiter": maxiter, "maxfun": int(maxiter * 1.5), "ftol": 1e-10, "gtol": 1e-6})
        pen = self.pen
        self.pen = []
        loss = float(self.Pm.run(self.build(r.x), self.d, self.wle, want_grad=False)[0])       # without the penalty
        self.pen = pen
        return {"z": r.x, "loss": loss, "loss_pen": float(r.fun), "nit": int(r.nit), "converged": bool(r.success)}

    # --------------------------------------------------------------------------------------------- projections and readouts
    def z_from_x(self, x):
        """A Problem2 x (the free layout, the same knots) carried into this parametrisation."""
        z = []
        for kind in ("F", "P"):
            inc = x[self.Pm.names == (self.Pm.fname if kind == "F" else self.Pm.pname)].reshape(self.nb, self.nk)
            sp = self.spec[kind]
            if sp == "free":
                z.append(inc.ravel())
            elif sp == "shared":
                amp = np.maximum(inc[:, :self.k5].sum(axis=1), 1e-3)
                s = (inc / amp[:, None]).mean(axis=0)
                u = np.maximum(s[:self.k5], 1e-5)
                u = u / u.sum()
                z.append(np.r_[amp, np.clip(np.log(u[1:] / u[0]), -12, 12), np.clip(s[self.k5:], 0, INC_HI)])
            else:
                raise ValueError("z_from_x: only free / shared")
        return np.r_[np.concatenate(z), x[self.ia], [x[self.Pm.names == p][0] for p in self.free]]

    def values(self, z):
        """F and P at the knots (nb x nk+1, 0 at R = 0) and a_link."""
        x = self.build(z)
        F = np.c_[np.zeros(self.nb), np.cumsum(x[self.Pm.names == self.Pm.fname].reshape(self.nb, self.nk), axis=1)]
        Pv = np.c_[np.zeros(self.nb), np.cumsum(x[self.Pm.names == self.Pm.pname].reshape(self.nb, self.nk), axis=1)]
        return F, Pv, np.exp(z[self.nc:self.nc + self.nl])

    def fz(self, z, ids):
        x = self.build(z)
        return {"x": x, "d": self.d, "a": dict(zip(ids, np.exp(z[self.nc:self.nc + self.nl]))), "model": "M4b", "M": self.M}


def shape_g(f, b):
    """The shape function of band b (a per-band mix uses its own form with the first parameters)."""
    if f.get("per_band"):
        sh = f["per_band"][b]
        n = len(sh["lo"])
        return lambda R, *t: np.asarray(sh["g"](R, *t[:n]), float) * np.ones_like(R)
    return lambda R, *t: np.asarray(f["g"](R, *t), float) * np.ones_like(R)


def readout(C_, z, bands):
    """The readable numbers of one fit: F, P at the knots per band, A / G / B / H, a_link percentiles."""
    F, Pv, a = C_.values(z)
    out = {}
    for b, nm in enumerate(bands):
        for j, k in enumerate(C_.kn[1:], 1):
            out[f"F[{nm}]({k:g})"] = float(F[b, j])
            out[f"P[{nm}]({k:g})"] = float(Pv[b, j])
    k5 = C_.k5
    for kind, V, A_, G_ in (("F", F, "A", "G"), ("P", Pv, "B", "H")):
        if C_.spec[kind] == "shared":
            for b, nm in enumerate(bands):
                out[f"{A_}[{nm}]"] = float(V[b, k5])
            for j, k in enumerate(C_.kn[1:], 1):
                out[f"{G_}({k:g})"] = float(V[0, j] / V[0, k5]) if V[0, k5] > 0 else np.nan
    la = np.log(a)
    out.update({"a_p10": float(np.exp(np.percentile(la, 10))), "a_p90": float(np.exp(np.percentile(la, 90))), "a_log_sd": float(la.std())})
    return out


# =================================================================================================== closed-form shapes
def _sat5(R, r0):
    return (R / (R + r0)) / (5 / (5 + r0))


def _exp5(R, r0):
    return (1 - np.exp(-R / r0)) / (1 - np.exp(-5 / r0))


def _capped(R, b, m):
    return np.minimum(np.power(R / 5, b), m) / min(1.0, m)


SHAPES = {   # name: g(R, *theta) with g(0) = 0, g(5) = 1; theta bounds and starts
    "sat5": {"kind": "F", "g": _sat5, "lo": [1e-3], "hi": [200.0], "starts": [[1.0], [3.0], [10.0]], "text": "F5 [R / (R + R0)] / [5 / (5 + R0)]", "names": ["R0"]},
    "exp5": {"kind": "F", "g": _exp5, "lo": [1e-2], "hi": [200.0], "starts": [[1.0], [3.0], [10.0]], "text": "F5 [1 - exp(-R / R0)] / [1 - exp(-5 / R0)]", "names": ["R0"]},
    "capped": {"kind": "F", "g": _capped, "lo": [0.05, 1e-3], "hi": [3.0, 100.0], "starts": [[0.5, 1.5], [0.4, 2.0], [0.7, 1.2]],
               "text": "F5 min((R / 5)^b, m) / min(1, m)", "names": ["b", "m"]},
    "linear": {"kind": "P", "g": lambda R: R / 5, "lo": [], "hi": [], "starts": [[]], "text": "P5 R / 5", "names": []},
    "power": {"kind": "P", "g": lambda R, b: np.power(R / 5, b), "lo": [0.2], "hi": [3.0], "starts": [[1.0], [0.8], [1.3]], "text": "P5 (R / 5)^b", "names": ["b"]},
    "quadratic": {"kind": "P", "g": lambda R, d: (R + d * R ** 2) / (5 + 25 * d), "lo": [0.0], "hi": [10.0], "starts": [[0.0], [0.05]],
                  "text": "P5 (R + d R^2) / (5 + 25 d)", "names": ["d"]},
}


def pysr_shape(expr_text, var_f=False):
    """A PySR equation as a shape: its numeric constants become parameters; g(R) = (h(R) - h(0)) / (h(5) - h(0)).
    Returns a SHAPES-like dict or None (no constant / not usable)."""
    import sympy as sp
    R, f = sp.symbols("R f")
    try:
        ex = sp.sympify(expr_text, locals={"R": R, "f": f})
    except Exception:                                        # noqa: BLE001
        return None
    consts = sorted({a for a in ex.atoms(sp.Float)}, key=lambda a: float(a))
    if var_f and f not in ex.free_symbols:
        return None
    if not var_f and f in ex.free_symbols:
        return None
    th = sp.symbols(f"t0:{len(consts)}") if consts else ()
    ex2 = ex.xreplace(dict(zip(consts, th)))
    args = (R, f) + tuple(th) if var_f else (R,) + tuple(th)
    h = sp.lambdify(args, ex2, "numpy")
    v0 = [float(c) for c in consts]
    lo = [min(c - 10 * abs(c) - 1, -1e-6) if c < 0 else max(c / 1e3, 1e-9) if c > 0 else -1 for c in v0]
    hi = [max(c + 10 * abs(c) + 1, 1e-6) if c > 0 else min(c / 1e3, -1e-9) if c < 0 else 1 for c in v0]
    if var_f:
        def g(Rv, fv, *t):
            with np.errstate(all="ignore"):
                Rv = np.asarray(Rv, float)
                return np.asarray(h(Rv, fv, *t), float) - np.asarray(h(0 * Rv, fv, *t), float)
    else:
        def g(Rv, *t):
            with np.errstate(all="ignore"):
                h0, h5 = np.asarray(h(np.array([0.0, 5.0]), *t), float) * np.ones(2)
                return (np.asarray(h(np.asarray(Rv, float), *t), float) * np.ones_like(Rv) - h0) / (h5 - h0)
    return {"g": g, "lo": lo, "hi": hi, "starts": [v0], "text": f"PySR: {expr_text}", "names": [str(t) for t in th], "structure": str(ex2)}


def usable(shape, theta):
    """>= 0 and non-decreasing on 0-40 mm/h (the rule of config closed_forms.every_form)."""
    r = np.linspace(0, 40, 801)
    with np.errstate(all="ignore"):
        v = np.asarray(shape["g"](r, *theta), float) * np.ones_like(r)
    return bool(np.isfinite(v).all() and (v >= -1e-9).all() and (np.diff(v) >= -1e-9).all())


# =================================================================================================== step S7: the frozen S6 models
def s6_models(Dt1, let, bands):
    """The two frozen descriptions of the S6 candidate, rebuilt from the S6 cache (read only, nothing refitted):
      "candidate"  the S6 plug-back "P only": F = the kept shared free knots (fixed), P = c_band R, a_link refitted
      "closed F"   the S6 plug-back "F and P": F = A_band g(R) with the PySR shape chosen in S6 (about (R / 5)^0.47), P = c_band R
    Both on the dense plug-back knot grid, dynamics and lag fixed. Returns {name: {"fz", "C", "z", "spec"}} and the kept
    shared fit (C, z); fz runs on any dataset with predict() (a_link of a link not in TRAIN = 1)."""
    r5 = pickle.loads((P.S5A_CACHE / "fit_M4b-B.pkl").read_bytes())
    x_free, d_free = np.asarray(r5["all"]["x"]), int(r5["all"]["d"])
    shared = {"F": "shared", "P": "shared"}
    Ck = Curves(Dt1, shared, x_free, d_free)
    zk = np.asarray(pickle.loads((P.S6_CACHE / "shared_FP.pkl").read_bytes())["z"])
    S6C = WA["tighten_s6"]["closed_forms"]["plug_back_knots"]
    Md = copy.deepcopy(M)
    Md["knots_mmh"] = [0.0] + list(np.geomspace(S6C["lo"], S6C["hi"], S6C["n"]))
    kn_d = np.asarray(Md["knots_mmh"])
    Pd, Psrc = sid.Problem2(Dt1, "M4b", Md), sid.Problem2(Dt1, "M4b", M)
    x_ref = sid.x_into(Pd, Psrc, Ck.build(zk), curves=False)
    F_k, _, _ = Ck.values(zk)
    forms = pd.read_csv(P.S6_FORMS)
    eq = forms[(forms.kind == "F") & forms.chosen].formula.iloc[0]
    assert eq.startswith("PySR: "), eq
    shF = {**pysr_shape(eq[len("PySR: "):]), "kind": "F"}
    lin = {**SHAPES["linear"]}
    ids = LTL.link_ids(Dt1, let)
    out = {}
    for name, cache, specF in (("candidate", "plug_shared_FP_P_only", ("fixed", np.array([np.interp(kn_d, Ck.kn, F_k[b]) for b in range(len(bands))]))),
                               ("closed F", "plug_shared_FP_F_and_P", ("form", shF, True))):
        spec = {"F": specF, "P": ("form", lin, True)}
        Cp = Curves(Dt1, spec, x_ref, Ck.d, Mx=Md)
        z = np.asarray(pickle.loads((P.S6_CACHE / f"{cache}.pkl").read_bytes())["z"])
        out[name] = {"fz": Cp.fz(z, ids), "C": Cp, "z": z, "spec": spec, "F_text": "free shared knots" if name == "candidate" else eq}
    return out, (Ck, zk)
