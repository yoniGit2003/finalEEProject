"""The one diagnostic after the failed fake test of the bin rules (step E3). Nothing is patched with it.

    python -m research.waa_ksurface.e3_diagnostic

The fake test of `e3_bin_rule.py` (no rain error) for both rules, with the memory of the store on / off and the noise on /
off. With neither memory nor noise only the method itself is left (a rain class holds a range of R, the cell sits at one
R_cell). The table shows which ingredient makes the intercept too high and the light-rain slopes wrong.
"""
import numpy as np
import pandas as pd

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import data, fakeworld as fw, ltest
from research.waa_ksurface import paths as P

LT, FK = ltest.config(), fw.config()
OUT = P.TABLES / "e3_diagnostic.csv"
SUMMARY = P.TABLES / "e3_diagnostic_summary.csv"
REPORT = P.REPORTS / "e3_diagnostic.md"


def main():
    D, le, n_cells = data.gauge_ratio_train()
    s = pd.read_csv(P.E3_SIZES).set_index("quantity").value
    sizes = {k: float(s[k]) for k in ("white_sd_db", "offset_sd_db", "link_sd_db")}
    n_short = int(float(s["n_short"]))
    parts, summ = [], []
    for name, rule in LT["rules"].items():
        for memory in (True, False):
            for noise in (True, False):
                tabs = []
                for rep in range(FK["reps"]):
                    Dk = fw.make(D, le, n_cells, sizes, 0.0, FK, FK["seed"] + rep, memory=memory, with_noise=noise)
                    tabs.append(ltest.cell_fits(ltest.steady_values(Dk, le, LT, rule), LT, n_short)[0])
                tab, n_fail, ok = fw.judge(tabs, FK)
                tab["s_err_abs"] = tab.s - tab.P_true
                parts.append(tab.assign(rule=name, memory=memory, noise=noise))
                summ.append({"rule": name, "memory": memory, "noise": noise, "cells": len(tab), "cells_outside": n_fail,
                             "I_err_median": float(tab.I_err.median()), "I_err_max_abs": float(tab.I_err.abs().max()), "I_outside": int((~tab.I_ok).sum()),
                             "s_outside": int((~tab.s_ok).sum()), "s_err_abs_max": float(tab.s_err_abs.abs().max()), "passes": ok})
                print(summ[-1], flush=True)
    out, summ = pd.concat(parts, ignore_index=True), pd.DataFrame(summ)
    save(out.round(5), OUT, index=False)
    save(summ.round(4), SUMMARY, index=False)
    wide = out[out.rule == "R1"].pivot_table(index=["band", "class"], columns=["memory", "noise"], values="I_err", sort=False).round(3)
    wide.columns = [f"I_err memory={'on' if m else 'off'} noise={'on' if n else 'off'}" for m, n in wide.columns]
    a = ["# E3 diagnostic: which ingredient of the fake world breaks the bin rule\n",
         f"Written by `research/waa_ksurface/e3_diagnostic.py` after the STOP of E3. No rain error in any run. {FK['reps']} replicates per row. "
         f"Nothing is patched with this table. (`{rel(SUMMARY)}`, `{rel(OUT)}`)\n",
         md_table(summ.round(3)), "\nIntercept error (dB) per cell, rule R1:\n", md_table(wide.reset_index())]
    save("\n".join(a) + "\n", REPORT)
    print(wide.to_string())


if __name__ == "__main__":
    main()
