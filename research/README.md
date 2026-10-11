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
What has been implemented so far, at the level of ideas: [`../docs/research_log.pdf`](../docs/research_log.pdf).

## Plan

1. **Baseline and A(t).** Estimate the dry-weather level of every sub-link from the certified dry periods of the
   release and with the minimum-attenuation method [8], then A(t) = TSL(t) - RSL(t) - baseline. This comes first:
   every later step uses A(t). The baseline must not follow the drying tail after rain, or the WAA disappears into it;
   what is left of the baseline in dry weather is part of δ(t).
2. **Stage 1, WAA.** Use the drying tail after rain ends (path attenuation is zero there, so what is left is WAA).
   Fit the literature models [1], [9]; learn a dynamic model with weak-form SINDy [10], rain as input; compare links.
3. **Stage 2, path law.** Learn f(L, R̄, Z) with PySR [7], with radar along the path as the reference (WAA removed
   first).
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

## Layout

```
config.yaml        all research settings (one section per module) and the pinned release
ROADMAP.md         research progress (the status block below is generated from it)
baseline_estm/     dry-weather baseline per sub-link and A(t)          (done: method v2 frozen)
WAA_estm/          Stage 1: drying tails (frozen tail set v1), then weak-form SINDy  (in progress)
path_law/          Stage 2: f(L, R̄, Z) with PySR                                  (planned)
evaluation/        rain estimation, metrics, comparisons                          (planned)
common/            shared helpers: outputs.py (save_output), tidy_outputs.py (not a module)
explorations/      playground: try ideas here; results that matter move into a module
```

Every module has the same layout (`README.md`, `METHOD.md`, `CATALOG.md`, `<module>_paths.py`, code, `validation/`,
`outputs/{products,tables,figures,reports}/`, git-ignored `.cache/`). The rules: [`../ORGANIZATION.md`](../ORGANIZATION.md);
what every output file means: the module's `CATALOG.md` (index: [`../CATALOG.md`](../CATALOG.md)).
Every output is written with `save_output(obj, path, level, purpose)` from `common/outputs.py`; the level (key /
product / diag / scratch, recorded in `outputs/MANIFEST.csv`) decides what is committed and what /tidy moves to the module's local `legacy/` folder
(ORGANIZATION.md section 5). Decisions from deleted step reports are kept in the module's `DECISIONS.md`.

## Conventions

- Notation identical to the proposal: A(t), A_rain, δ(t), R̄, L, a, b, f(L, R̄, Z), Z.
- References [n] are those of the proposal; check claims against them, never invent citations.
- Prefer scripts that can be rerun over one-off notebooks; record the release version in every result.

## Status and next steps *(auto, from [`ROADMAP.md`](ROADMAP.md))*

<!-- AUTO:progress -->
Status: **in progress (baseline frozen; Stage 1 WAA: drying tails frozen, fitting tools tested on simulated tails)**; 4 / 16 roadmap items done (tick items in [`ROADMAP.md`](ROADMAP.md)).

- 1. Baseline and attenuation A(t): 3 / 4
- 2. Stage 1: WAA: 1 / 5
- 3. Stage 2: path law f(L, R̄, Z): 0 / 1
- 4. Rain estimation and evaluation: 0 / 6

Next steps:
- [ ] Baseline per sub-link from the certified dry periods (constant / rolling) and the dynamic minimum-attenuation method [8]
- [ ] Fit literature WAA models [1], [9] per link
- [ ] Learn a dynamic WAA model with weak-form SINDy (rain as input)
- [ ] Compare the learned models across links (same form, different coefficients?)
- [ ] Drying time from SINDy vs drying time measured directly from tails
- [ ] Learn f(L, R̄, Z) with PySR, radar along the path as reference (WAA removed first)
- ... and 6 more in the roadmap
<!-- /AUTO:progress -->
