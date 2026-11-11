"""Step E3 of the k-surface branch: the bin definitions of the length test, derived here and tried on fake data.

    python -m research.waa_ksurface.e3_bin_rule

What "steady rain" must mean so that the intercept of the length test is the steady WAA F(R): the rain has been near the
bin's own R for a while (so W has reached its level), with a band that is symmetric in log R (so the bin is not entered
systematically from heavier rain). Two rules are tried in order (R1: 20 min within a factor 2; R2: 30 min within a factor
1.5), in a fake world that has a known F and P, a store with memory and the measured noise, but NO rain error: the rule is
judged on its own. The script also measures the sizes the fake world needs (noise on dry minutes, spread between links) and
derives n_short, the number of short links a cell needs.
"""
import sys

import numpy as np
import pandas as pd

from research.common.report import md_table, rel
from research.common.save import save
from research.waa_ksurface import data, fakeworld as fw, ltest
from research.waa_ksurface import paths as P

LT, FK = ltest.config(), fw.config()


def link_spread(D, le, rule):
    """The between-link sd of the short links' residuals from the plain cell lines (real TRAIN data, gauge-ratio events):
    per link and storm the mean residual over the rain classes, then variance components over links."""
    t = ltest.steady_values(D, le, LT, rule)
    cells, _ = ltest.cell_fits(t, LT, LT["n_short_floor"])
    r = ltest.residuals(t, cells)
    r = r[(r.L_km < LT["short_km"]) & r.resid.notna()]
    per = r.groupby(["cml_id", "event_id"]).resid.mean().reset_index()
    between, within, n_links, n_values = ltest.variance_components(per.resid.to_numpy(), per.cml_id.to_numpy())
    return cells, {"link_sd_db": between, "storm_sd_db": within, "short_links": n_links, "link_storm_values": n_values}


def fake_test(D, le, n_cells, sizes, rule, n_short):
    """The plain length test in the fake world without rain error (memory and noise only)."""
    tabs = []
    for rep in range(FK["reps"]):
        Dk = fw.make(D, le, n_cells, sizes, 0.0, FK, FK["seed"] + rep)
        tabs.append(ltest.cell_fits(ltest.steady_values(Dk, le, LT, rule), LT, n_short)[0])
    return fw.judge(tabs, FK)


def main():
    d = data.open_release()
    D, le, n_cells = data.gauge_ratio_train(d)
    D_all, _, le_all, _, _ = data.train(d)
    print(f"data: {len(le)} TRAIN link-events of {le.cml_id.nunique()} links in {le.event_id.nunique()} events with a gauge ratio (28 and 38 GHz)", flush=True)

    # ---- sizes measured on real data
    sizes = fw.noise_sizes(D, le)
    real_cells, spread = link_spread(D, le, LT["rules"]["R1"])
    sizes.update(spread)
    n_short = int(max(LT["n_short_floor"], np.ceil((spread["link_sd_db"] / LT["tolerance_db"]) ** 2)))
    print(f"sizes: {sizes}; n_short = {n_short}", flush=True)
    short32 = int(le_all[(le_all.band_name == "32 GHz") & (le_all.L_km < LT["short_km"])].cml_id.nunique())

    # ---- the fake test of the rules, in order
    parts, adopted, verdicts = [], None, {}
    real_parts = [real_cells.assign(rule="R1")]
    for name, rule in LT["rules"].items():
        tab, n_fail, ok = fake_test(D, le, n_cells, sizes, rule, n_short)
        verdicts[name] = (len(tab), n_fail, ok)
        parts.append(tab.assign(rule=name))
        print(f"rule {name} {rule}: {len(tab)} cells judged, {n_fail} outside -> {'PASS' if ok else 'fail'}\n"
              + tab[["band", "class", "R_cell", "F_true", "I", "I_err", "s_err_rel", "links_short", "I_ok", "s_ok"]].round(3).to_string(), flush=True)
        if name != "R1":
            real_parts.append(ltest.cell_fits(ltest.steady_values(D, le, LT, rule), LT, n_short)[0].assign(rule=name))
        if ok and adopted is None:
            adopted = name
            break
    fake = pd.concat(parts, ignore_index=True)
    save(fake.round(5), P.E3_FAKE, index=False)
    save(pd.concat(real_parts, ignore_index=True)[["rule", "band", "class", "link_events", "links", "events", "bins", "links_short", "links_long", "fitted"]],
         P.E3_REAL_CELLS, index=False)
    out = pd.DataFrame([{"quantity": k, "value": v} for k, v in {**sizes, "n_short": n_short, "adopted_rule": adopted or "none",
                                                                 "links_32GHz_below_1km": short32}.items()])
    save(out, P.E3_SIZES, index=False)

    a = ["# E3: the bin definitions of the length test, and their fake test\n",
         f"Written by `research/waa_ksurface/e3_bin_rule.py`. {len(le)} TRAIN link-events, {le.cml_id.nunique()} links, {le.event_id.nunique()} events with a gauge ratio "
         "(28 and 38 GHz). The definitions and their reasons are in `docs/waa_ksurface/decisions.md` (step E3). VAL and TEST were not opened.\n",
         "## Sizes measured on the real TRAIN data\n",
         md_table(pd.DataFrame([
             {"what": "white noise of a 5 min value (dry offset window)", "value": round(sizes["white_sd_db"], 3), "unit": "dB", "from": f"{sizes['white_link_events']} link-events"},
             {"what": "offset error of a link-event (dry 30 min before the rain)", "value": round(sizes["offset_sd_db"], 3), "unit": "dB", "from": f"{sizes['offset_link_events']} link-events"},
             {"what": "between-link sd of the short links' residuals", "value": round(sizes["link_sd_db"], 3), "unit": "dB", "from": f"{sizes['short_links']} links, {sizes['link_storm_values']} link-storm values"},
             {"what": "within-link sd (storm to storm)", "value": round(sizes["storm_sd_db"], 3), "unit": "dB", "from": "same"}])),
         f"\n**n_short = max({LT['n_short_floor']}, ceil(({sizes['link_sd_db']:.3f} / {LT['tolerance_db']})²)) = {n_short}** distinct links below {LT['short_km']:g} km, "
         f"and as many of {LT['short_km']:g} km or more, for a cell to be eligible.\n",
         f"32 GHz has {short32} links below {LT['short_km']:g} km, so by this rule it has no eligible cell: **32 GHz stays excluded from the branch for every model.**\n",
         "## The fake test of the rules (memory and noise, no rain error)\n",
         f"Truth: F = {FK['truth']['f_top_db']:g}·R / (R + {FK['truth']['f_half_mmh']:g}) dB, P = 0.1·R / 0.2·R dB/km, one store with {FK['truth']['tau_wet_min']:g} min (wetting) and "
         f"{FK['truth']['tau_dry_min']:g} min (drying). {FK['reps']} replicates, medians. Pass: intercept within {FK['pass']['intercept_db']} dB and slope within "
         f"{100 * FK['pass']['slope_rel']:.0f} % in all judged cells but at most {FK['pass']['cells_may_fail']}.\n"]
    for name, (n, n_fail, ok) in verdicts.items():
        r = LT["rules"][name]
        a.append(f"### Rule {name}: {r['t_pre_min']} min within a factor {r['q']:g} — {n} cells judged, {n_fail} outside: {'PASS' if ok else 'FAIL'}\n")
        a.append(md_table(fake[fake.rule == name][["band", "class", "R_cell", "F_true", "I", "I_err", "P_true", "s", "s_err_rel", "links_short", "I_ok", "s_ok"]].round(3)))
        a.append("")
    a.append(f"**Adopted rule: {adopted or 'none (STOP)'}.** (`{rel(P.E3_FAKE)}`, `{rel(P.E3_SIZES)}`)\n")
    save("\n".join(a) + "\n", P.E3_REPORT)
    print(f"adopted: {adopted}; report: {rel(P.E3_REPORT)}", flush=True)
    if adopted is None:
        sys.exit(1)


if __name__ == "__main__":
    main()
