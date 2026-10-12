"""The datasets the fitting steps work on: the noise tables, tail sets on disk, event bootstraps.

Used by the method tests on simulated tails (steps 3b-1b to 3b-4) and by the fit of the real tails (step 4).
Settings come from research/config.yaml: `wet_antenna.simulator` and `wet_antenna.stress_pooled_fit.t1_noise`.
"""
import numpy as np
import pandas as pd

from research.common.config import load_config
from research.wet_antenna import paths as P
from research.wet_antenna import simulator as ws

CFG = load_config()
SIM, S1C = CFG["wet_antenna"]["simulator"], CFG["wet_antenna"]["stress_pooled_fit"]
STEP, SPM, LAGS = SIM["rsl_step_db"], SIM["samples_per_min"], SIM["noise_acf_lags_min"]


def load_noise():
    """Noise of the FIT (3a-fix, one AR(1) per sub-link) and the offset prior per sub-link.
    Returns (the noise table, {sub-link: (sigma_white, sigma_slow, tau_noise)}, {sub-link: offset std})."""
    nz = pd.read_csv(P.SIM_NOISE)
    noise = {r.sublink_id: (r.sigma_white_db, r.sigma_slow_db, r.tau_noise_min) for r in nz.itertuples()}
    offset_sd = pd.read_csv(P.BASELINE_UNCERTAINTY).set_index("sublink_id").offset_std.to_dict()
    return nz, noise, offset_sd


def load_tails(sim):
    """A tail set stored as 1 min rows (parquet file `sim`) as a list of tails, the input of pooled_fit.pack."""
    rows = pd.read_parquet(sim).sort_values(["tail_id", "minutes_since_t0_trace"])
    return [{"tail_id": t, "sublink_id": g.sublink_id.iloc[0], "cml_id": g.cml_id.iloc[0],
             "event_copy": int(g.event_copy.iloc[0]), "A_start": g.A_start_trace_db.iloc[0],
             "t": g.minutes_since_t0_trace.to_numpy(), "a": g.A_db.to_numpy(float)} for t, g in rows.groupby("tail_id")]


def bootstrap_skeleton(skel, rng):
    """Draw events with replacement and take all their tails; new tail ids, event_copy = draw number."""
    evs = np.sort(skel.event_id.unique())
    parts = []
    for k, e in enumerate(rng.choice(evs, len(evs), replace=True)):
        g = skel[skel.event_id == e].copy()
        g["orig_tail_id"] = g.tail_id
        g["tail_id"] = (k + 1) * 100000 + g.tail_id
        g["event_copy"] = k
        parts.append(g)
    return pd.concat(parts, ignore_index=True)


def shape(df, a0, minutes):
    """Median over tails of A / A_start at k minutes (mean of A over [k, k + 5)), as in step 3b-1b."""
    out = {}
    for k in minutes:
        v = df[(df.minutes_since_t0_trace >= k) & (df.minutes_since_t0_trace < k + 5)].groupby("tail_id").A_db.mean()
        out[k] = float((v / v.index.map(a0)).median())
    return out


def t1_noise(nz, block):
    """T1 noise: per sub-link the same white sigma and the same TOTAL slow variance as the 3a-fix noise (median
    0.15^2 dB^2), but split over two AR(1) wanders (tau_fast, tau_slow). The share of the fast part is one number for
    the whole network, chosen on a grid so that dry simulated series match the real control-tail statistics
    (within-tail std, ACF at 1 / 10 / 60 min). Returns (share, check table, per-share scores)."""
    T = S1C["t1_noise"]
    rng = np.random.default_rng(T["seed"])
    nb = T["blocks"]
    eps = rng.standard_normal((nb, block, SPM))                    # the same unit random numbers for every sub-link
    ph = rng.uniform(0, STEP, (nb, 1, 1))
    unit = [np.stack([ws.ar1(rng.standard_normal(block), 1.0, tau, rng.standard_normal()) for _ in range(nb)])
            for tau in (T["tau_fast_min"], T["tau_slow_min"])]
    fz = nz[nz.source == "fit"]
    cols = ["within_std_db", "on_level_share"] + [f"acf_{k}" for k in LAGS]
    scale = {"within_std_db": SIM["noise_weights"]["std_db"], **{f"acf_{k}": SIM["noise_weights"]["acf"] for k in LAGS}}
    sims, scores = {}, []
    for w in T["fast_share_grid"]:
        slow = (np.sqrt(w) * unit[0] + np.sqrt(1 - w) * unit[1])[:, :, None]
        rows = []
        for r in fz.itertuples():
            q = ws.quantize(r.sigma_slow_db * slow + r.sigma_white_db * eps, ph, STEP)
            rows.append(ws.dry_stats(q.mean(axis=2), q, LAGS))
        sims[w] = pd.DataFrame(rows, index=fz.index)
        err = {c: float((sims[w][c] - fz[c]).abs().median()) for c in scale}
        scores.append({"fast_share": w, **{f"abs_err_{c}": round(v, 4) for c, v in err.items()},
                       "score": round(sum((err[c] / scale[c]) ** 2 for c in scale), 4)})
    scores = pd.DataFrame(scores)
    w = float(scores.loc[scores.score.idxmin(), "fast_share"])
    chk = pd.DataFrame([{"statistic": c, "real_median": fz[c].median(), "sim_3afix_median": fz[f"sim_{c}"].median(),
                         "sim_T1_median": sims[w][c].median(),
                         "abs_err_3afix": (fz[f"sim_{c}"] - fz[c]).abs().median(),
                         "abs_err_T1": (sims[w][c] - fz[c]).abs().median()} for c in cols]).round(4)
    return w, chk, scores
