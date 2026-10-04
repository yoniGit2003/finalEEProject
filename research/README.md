# research/

The research part of the project. *[Placeholder: supervisor, dates, links to the report.]*

## The question

Rain along a microwave link attenuates its signal. The measured attenuation is A(t) = A_rain(t) + δ(t), where
A_rain is the loss caused by rain along the path and δ(t) collects everything else, mainly wet-antenna attenuation
(WAA): water on the antenna covers, which builds up in rain and dries afterwards. The common model is
A_rain ≈ a·R̄^b·L, which is exact only for uniform rain. We want to find, from data, (1) a dynamic model of WAA and
(2) the true dependence on the link, A_rain = a·R̄^b·f(L, R̄, Z), where uniform rain gives f = L. WAA must be removed
before learning f, because a constant WAA looks like f = L + L0. The aim is readable equations that keep the known
physics (gray-box), not only good predictions. Background: [`../docs/proposal.md`](../docs/proposal.md).

## Plan

1. **Simulator** with known ground truth: rain fields with a controlled cell size, the exact path integral for every
   link, a known WAA model, noise and quantization (RSL 0.3 dB, TSL 1 dB). Every method is tested here first.
2. **Stage 1, WAA.** Use the drying tail after rain ends (path attenuation is zero there, so what is left is WAA).
   Fit the literature models [1], [9]; learn a dynamic model with weak-form SINDy [10], rain as input; compare links.
3. **Stage 2, path law.** Compare the exact integral with a·R̄^b·L in simulation, learn f(L, R̄, Z) with PySR [7],
   then use real data with radar along the path as the reference (WAA removed first).
4. **Rain estimation and evaluation.** Invert the combined model on events not used for training; normalized bias
   and RMSE against gauges and radar by link length; false rain after events; comparison with the power law, the
   literature WAA models and an RNN baseline [11]; stability of the learned terms (bootstrap).

Links must never be chosen by how well the power law fits them (that would assume the answer).

## Getting the data

Research code uses the data only through `cmlpipe.api.CMLData` and the release pinned in [`config.yaml`](config.yaml):

```python
import yaml
from cmlpipe.api import CMLData
d = CMLData.open(f"data/release/{yaml.safe_load(open('research/config.yaml'))['release']}")   # from the repo root
```

Allowed imports from the pipeline: `cmlpipe.api` and `cmlpipe.viz` only, and the pipeline never imports research
(checked by `python pipeline/tools/check_boundaries.py`). What the data contains and how to get it:
[`../pipeline/docs/DATA_GUIDE.md`](../pipeline/docs/DATA_GUIDE.md).

## Planned layout

```
simulator/    rain fields, path integral, WAA, noise and quantization
waa/          Stage 1: drying tails, literature fits, weak-form SINDy
path_law/     Stage 2: f(L, R̄, Z) with PySR
evaluation/   rain estimation, metrics, baselines
notebooks/    exploration (results that matter go into scripts)
figures/      figures for the report
```

## Conventions

- Notation identical to the proposal: A(t), A_rain, δ(t), R̄, L, a, b, f(L, R̄, Z), Z.
- References [n] are those of the proposal; check claims against them, never invent citations.
- Prefer scripts that can be rerun over one-off notebooks; record the release version in every result.

## Status and next steps *(auto, from [`ROADMAP.md`](ROADMAP.md))*

<!-- AUTO:progress -->
Status: **in progress (not started in code yet)**; 0 / 18 roadmap items done (tick items in [`ROADMAP.md`](ROADMAP.md)).

- 2. Simulator (known ground truth): 0 / 4
- 3. Stage 1: WAA: 0 / 5
- 4. Stage 2: path law f(L, R̄, Z): 0 / 3
- 5. Rain estimation and evaluation: 0 / 6

Next steps:
- [ ] Rain fields with controlled cell size
- [ ] Exact path integral of a·R(x,t)^b along each link
- [ ] Known WAA model (wetting and drying)
- [ ] Noise and quantization (RSL 0.3 dB, TSL 1 dB)
- [ ] Find links with a drying tail after rain ends
- [ ] Fit literature WAA models [1], [9] per link
- ... and 12 more in the roadmap
<!-- /AUTO:progress -->
