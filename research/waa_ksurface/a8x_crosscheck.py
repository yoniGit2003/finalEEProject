"""Naming cross-check: two independent tools on the step 8 input (decisions.md, "Naming cross-check", logged before the run).

    python -m research.waa_ksurface.a8x_crosscheck sparse      tool 1: best-subset sparse regression (63 subsets of 7 terms)
    python -m research.waa_ksurface.a8x_crosscheck bms         tool 2: the Bayesian Machine Scientist, one run per k_wet estimate
    python -m research.waa_ksurface.a8x_crosscheck report      agreement, the plug-back queue (reads the step 8 cache), the report

Same input and scoring as step 8 (1)-(2): the 30 A4 drying points, density weights, y = ln k_dry of the A7 K1 fit at both
estimates; every formula refitted by `a8_forms.fit_curve` (weighted, ln k) and ranked by `a8_forms.rank`. Light processes:
below-normal priority; a process stops itself if free memory falls below 1.5 GB. The step 8 run is not touched.
"""
import itertools
import os
import pickle
import re
import sys
import threading
import time

import numpy as np
import pandas as pd
import psutil
import sympy as sp

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a8_forms as F8, a8_name as A8
from research.waa_ksurface import paths as P

C = F8.data.branch_config()["a8x"]
EST = A8.EST
CACHE = P.CACHE / "a8x"
BMS_DIR = CACHE / "bms"
GB = 1024 ** 3
W = F8.W
LIB = {"1": sp.Integer(1), "W": W, "W^2": W ** 2, "1/W": 1 / W, "sqrt(W)": sp.sqrt(W), "exp(-W)": sp.exp(-W), "W/(1+W)": W / (1 + W)}


def log(s):
    print(f"{time.strftime('%H:%M:%S')} {s}", flush=True)


def light():
    """Below-normal priority, and a watch that stops this process below the memory floor."""
    psutil.Process().nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)

    def watch():
        while True:
            free = psutil.virtual_memory().available / GB
            if free < C["min_free_gb"]:
                log(f"=== free memory {free:.2f} GB below {C['min_free_gb']} GB: stopped ===")
                os._exit(3)
            time.sleep(5)
    threading.Thread(target=watch, daemon=True).start()


def scored(name, form, pts, source, extra=None, starts=None):
    """One candidate row: refitted at both estimates by the step 8 objective, its AICc."""
    row = {"name": name, "source": source, "p": form.p, "form": form.text, "struct": A8.structure(str(form.expr)) if form.p else str(form.expr), **(extra or {})}
    th = {}
    for e in EST:
        th[e], wrss = F8.fit_curve(form, pts["W"], pts["y"][e], pts["w"], starts=(starts or {}).get(e))
        row[f"wrss_{e}"], row[f"aicc_{e}"], row[f"formula_{e}"] = wrss, F8.aicc(wrss, len(pts["W"]), form.p), form.formula(th[e])
    return row, th


# =================================================================================================== tool 1
def sparse():
    pts = A8.points_part()
    rows, fits = [], {}
    for r in (1, 2, 3):
        for terms in itertools.combinations(LIB, r):
            cs = sp.symbols(f"c0:{r}")
            expr = sum(c * LIB[t] for c, t in zip(cs, terms))
            f = F8.Form("sparse: " + " + ".join(terms), expr, cs, [np.zeros(r)], "k = " + " + ".join(f"c{i}*{t}" for i, t in enumerate(terms)), source="sparse")
            Th = np.column_stack([np.asarray(sp.lambdify(W, LIB[t], "numpy")(pts["W"]), float) * np.ones(len(pts["W"])) for t in terms])
            starts = {}
            for e in EST:                                        # weighted linear LS of k, weights w / k^2 (first order of the ln k loss)
                k = np.exp(pts["y"][e])
                sw = np.sqrt(pts["w"]) / k
                starts[e] = [np.linalg.lstsq(Th * sw[:, None], k * sw, rcond=None)[0]]
            row, th = scored(f.name, f, pts, "sparse", {"terms": " + ".join(terms)}, starts)
            rows.append(row)
            fits[f.name] = {"terms": terms, "theta": th}
    t, order = F8.rank(pd.DataFrame(rows))
    log(f"sparse: chosen {order[0]}")
    return {"table": t, "order": order, "fits": fits}


# =================================================================================================== tool 2
def bms_expr_to_form(name, expr):
    """A BMS expression (parameters _a0_, ...) as a Form; its parameters are those that appear."""
    names = sorted(set(re.findall(r"_a\d+_", expr)), key=lambda s: int(s[2:-1]))
    txt = expr
    for n_ in names:
        txt = txt.replace(n_, n_.strip("_"))
    syms = sp.symbols(" ".join(n_.strip("_") for n_ in names)) if names else ()
    syms = list(syms) if isinstance(syms, (tuple, list)) else [syms]
    ex = sp.sympify(txt, locals={"W": W, **{str(s): s for s in syms}})
    syms = [s for s in syms if s in ex.free_symbols]
    starts = [np.ones(len(syms)), np.full(len(syms), 0.1), np.full(len(syms), -1.0), np.full(len(syms), 0.5)]
    return F8.Form(name, ex, syms, starts, f"k = {sp.sstr(ex)}", source="BMS")


def bms():
    sys.path.insert(0, str(BMS_DIR))
    from bms import BMS                                        # noqa: E402  (the authors' code, pinned commit)
    pts = A8.points_part()
    x = pd.DataFrame({"W": pts["W"]})
    out = {}
    for e in EST:
        t0 = time.time()
        y = pd.Series(np.exp(pts["y"][e]))
        m = BMS(x, y, npar=C["bms"]["npar"], prior_par=str(BMS_DIR / "prior_win" / C["bms"]["prior_file"]), ops=C["bms"]["ops"],
                max_size=C["bms"]["max_size"], seed=C["bms"]["seed"], quiet=True)
        tr = m.sample(nsteps=C["bms"]["burnin"] + C["bms"]["steps"], quiet=True)
        post = tr.burn(C["bms"]["burnin"]).model_posterior()
        post["probability"] = post["count"] / post["count"].sum()
        top = post.sort_values("count", ascending=False).head(C["bms"]["top"])
        out[e] = {"top": top.reset_index(), "mdl": str(m.mdl_model), "mdl_E": float(m.mdl), "samples": int(post["count"].sum()),
                  "minutes": (time.time() - t0) / 60}
        log(f"bms {e}: {out[e]['samples']} samples after burn-in, {(time.time() - t0) / 60:.1f} min; MDL model {out[e]['mdl']}")
    return out


# =================================================================================================== report
def step8_reference():
    """(name, row, which) of step 8's form: the accepted one if step 8 has finished with one, else its chosen form."""
    rk = A8.load("rank")
    t = rk["table"].set_index("name")
    acc = A8.accepted() if (A8.CACHE / "plug_0.pkl").exists() else None
    name = acc["name"] if acc else rk["order"][0]
    return name, t.loc[name], "accepted" if acc else ("chosen (tried, not accepted)" if (A8.CACHE / "plug_0.pkl").exists() else "chosen (plug-back not finished)")


def report():
    sp_ = pickle.loads((CACHE / "sparse.pkl").read_bytes())
    bm = pickle.loads((CACHE / "bms.pkl").read_bytes()) if (CACHE / "bms.pkl").exists() else None
    pts = A8.points_part()
    L = ["# Naming cross-check: sparse regression and the Bayesian Machine Scientist\n",
         "Written by `research/waa_ksurface/a8x_crosscheck.py` (decisions.md, \"Naming cross-check\", logged before the run; disclosure there: the step 8 "
         "ranking had been seen, its plug-back and bootstrap not). Input and scoring as step 8 (1)-(2). VAL and TEST were not opened.\n"]
    st = sp_["table"]
    save(st, P.TABLES / "a8x_crosscheck_sparse.csv", index=False)
    ch = st[st.chosen].iloc[0]
    cols = ["name", "p", "aicc_gauge", "aicc_radar", "daicc_gauge", "daicc_radar", "dmax", "chosen", "formula_radar"]
    L += ["## Tool 1: sparse regression (best subset, 63 subsets)\n", f"Chosen: **{ch['name']}** (p = {ch.p}); AICc gauge {ch.aicc_gauge:.2f}, radar "
          f"{ch.aicc_radar:.2f}; {ch.formula_gauge} (gauge); {ch.formula_radar} (radar).\n", "The 10 lowest Delta_max:\n",
          md_table(st.sort_values("dmax").head(10)[cols].round(3)), ""]
    rows = [ch.to_dict()]
    best_b = None
    if bm is not None:
        L += ["## Tool 2: the Bayesian Machine Scientist (bitbucket.org/rguimera/machine-scientist, commit f4cf1e8)\n",
              f"Prior `{C['bms']['prior_file']}` (nv = 1, np = {C['bms']['npar']}), operators {C['bms']['ops']}, max_size {C['bms']['max_size']}, seed "
              f"{C['bms']['seed']}, {C['bms']['steps']} steps after {C['bms']['burnin']} burn-in. Probability = share of the samples after burn-in (BMS's "
              "own unweighted likelihood on k). Each formula refitted and scored by the step 8 objective.\n"]
        brows = []
        for e in EST:
            for i, r in bm[e]["top"].iterrows():
                f = bms_expr_to_form(f"bms_{e}_{i + 1}", r["expr"])
                row, _ = scored(f.name, f, pts, "BMS", {"bms_estimate": e, "rank": i + 1, "probability": float(r["probability"]), "bms_expr": r["expr"]})
                brows.append(row)
        bt = pd.DataFrame(brows)
        bt, _ = F8.rank(bt)
        save(bt, P.TABLES / "a8x_crosscheck_bms.csv", index=False)
        L += [md_table(bt[["bms_estimate", "rank", "probability", "bms_expr", "p", "aicc_gauge", "aicc_radar", "dmax", "chosen"]].round(4)), "",
              f"MDL models: gauge `{bm['gauge']['mdl']}`, radar `{bm['radar']['mdl']}`.\n"]
        best_b = bt[bt.chosen].iloc[0]
        rows.append(best_b.to_dict())
    # agreement and the plug-back queue
    name8, r8, which = step8_reference()
    L += ["## Agreement and the plug-back queue\n", f"Step 8's form: **{name8}** ({which}); AICc gauge {r8.aicc_gauge:.2f}, radar {r8.aicc_radar:.2f}, p = {r8.p}.\n"]

    def agree(a, b):
        same = bool(a["struct"] == b["struct"])
        close = bool(abs(a["aicc_gauge"] - b["aicc_gauge"]) <= 2 and abs(a["aicc_radar"] - b["aicc_radar"]) <= 2)
        return same or close, same, (a["aicc_gauge"] - b["aicc_gauge"], a["aicc_radar"] - b["aicc_radar"])
    ag = []
    pairs = [("sparse", ch, "step 8", r8)] + ([("BMS", best_b, "step 8", r8), ("sparse", ch, "BMS", best_b)] if best_b is not None else [])
    for na, a, nb, b in pairs:
        ok, same, d = agree(a, b)
        ag.append({"pair": f"{na} vs {nb}", "same family": same, "dAICc gauge": d[0], "dAICc radar": d[1], "agree": ok})
    L += [md_table(pd.DataFrame(ag).round(2)), ""]
    allc = pd.concat([st] + ([bt] if bm is not None else []), ignore_index=True)
    q = allc[(allc.p <= r8.p) & (allc.aicc_gauge < r8.aicc_gauge - 2) & (allc.aicc_radar < r8.aicc_radar - 2)]
    save(q, P.TABLES / "a8x_crosscheck_queue.csv", index=False)
    L += [f"Queued for plug-back (p <= {r8.p} and AICc lower by more than 2 at both estimates): **{len(q)}**"
          + (": " + ", ".join(q.name) if len(q) else "") + ". They run after step 8 has finished and at least 4 GB is free.\n"]
    save("\n".join(L) + "\n", P.REPORTS / "a8x_crosscheck.md")
    log("report written")


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    part = sys.argv[1]
    if part == "report":
        return report()
    light()
    free = psutil.virtual_memory().available / GB
    log(f"=== {part} started; free memory {free:.1f} GB ===")
    if free < C["min_free_gb"]:
        log("=== not started: below the memory floor ===")
        sys.exit(3)
    r = sparse() if part == "sparse" else bms()
    (CACHE / f"{part}.pkl").write_bytes(pickle.dumps(r))
    log(f"=== {part} done; peak memory {psutil.Process().memory_info().peak_wset / GB:.2f} GB ===")


if __name__ == "__main__":
    main()
