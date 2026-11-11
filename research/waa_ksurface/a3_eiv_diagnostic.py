"""SUPERSEDED (Detachment correction, 2026-11-11): the diagnostic of the superseded step A3d. Kept for the record.

The one diagnostic after the failed errors-in-variables check of step A3d (nothing is patched with it).

    python -m research.waa_ksurface.a3_eiv_diagnostic

The same fake check as `a3_length_test.eiv_check`, with the rain error scaled by 0, 0.25, 0.5 and 1. If the errors of the
length test vanish without a rain error, the cause of the failure is the rain error itself.
"""
import numpy as np
import pandas as pd

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import a3_length_test as A3
from research.waa_ksurface import data, lt
from research.waa_ksurface import paths as P

SCALES = [0.0, 0.25, 0.5, 1.0]
REPS = 5
OUT = P.TABLES / "a3_eiv_diagnostic.csv"
REPORT = P.REPORTS / "a3_eiv_diagnostic.md"


def main():
    _, D1, le, _ = data.branch_train()
    cells = pd.read_csv(P.A3_CELLS)
    noise = pd.read_csv(P.A3_EIV_NOISE).iloc[0].to_dict()
    A3.C["eiv"]["reps"] = REPS                                  # fewer replicates (in memory only); same seeds as the first replicates of A3d
    parts = []
    for sc in SCALES:
        print(f"rain error x {sc}", flush=True)
        tab, _ = A3.eiv_check(D1, le, cells, {**noise, "sigma_1": sc * noise["sigma_1"]})
        parts.append(tab.assign(scale=sc, sigma_1=sc * noise["sigma_1"]))
    out = pd.concat(parts, ignore_index=True)
    save(out.round(5), OUT, index=False)
    E = A3.C["eiv"]["pass"]
    summ = out.groupby("scale").agg(sigma_1=("sigma_1", "first"), largest_I_err_db=("I_err", lambda v: float(v.abs().max())),
                                    largest_s_err_pct=("s_err_rel", lambda v: float(100 * v.abs().max())),
                                    cells_I_ok=("I_ok", "sum"), cells_s_ok=("s_ok", "sum"), cells=("I_ok", "size")).reset_index()
    wide = out.pivot_table(index=["band", "class"], columns="scale", values="I_err", sort=False).round(3)
    wide.columns = [f"I_err at x{c:g}" for c in wide.columns]
    a = ["# A3d diagnostic: the fake check against the size of the rain error\n",
         f"Written by `research/waa_ksurface/a3_eiv_diagnostic.py` after the STOP of A3d. The fake check of A3d with the rain error scaled "
         f"(σ = scale × {noise['sigma_1']:.3f} / √n_pix), {REPS} replicates per scale. Pass limits: {E['intercept_db']} dB and {100 * E['slope_rel']:.0f} %. "
         f"Nothing is patched with this table. (`{rel(OUT)}`)\n",
         md_table(summ.round(3)), "\nIntercept error (dB) per cell:\n", md_table(wide.reset_index())]
    save("\n".join(a) + "\n", REPORT)
    print(summ.round(3).to_string())
    print(wide.to_string())


if __name__ == "__main__":
    main()
