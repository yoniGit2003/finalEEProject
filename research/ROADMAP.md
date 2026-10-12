# Research roadmap

What is done and what is left in the research (research/), following the work scope of `docs/proposal.md`
(section 6). References [n] are those of the proposal. Tick an item by changing `- [ ]` to `- [x]`: the status block
of `research/README.md` and the one-line summary in the root `README.md` are generated from this file on every commit
(`pipeline/tools/update_readme.py`), so this is the only place where research progress is edited by hand.

Status: in progress (baseline frozen; Stage 1 WAA: drying tails frozen, fitting tools tested on simulated tails, first fit of the real tails done: no drying law identified yet)

## 1. Baseline and attenuation A(t)
- [ ] Baseline per sub-link from the certified dry periods (constant / rolling) and the dynamic minimum-attenuation method [8]
- [x] Check that the baseline does not follow the drying tail (otherwise WAA ends up in the baseline)
- [x] Baseline residual in dry periods: size and spread per sub-link (the noise floor of δ(t))
- [x] A(t) = TSL(t) - RSL(t) - baseline per sub-link, as the input of the next stages

## 2. Stage 1: WAA
- [x] Find links with a drying tail after rain ends
- [ ] Fit literature WAA models [1], [9] per link
- [ ] Learn a dynamic WAA model with weak-form SINDy (rain as input)
- [ ] Compare the learned models across links (same form, different coefficients?)
- [ ] Drying time from SINDy vs drying time measured directly from tails

## 3. Stage 2: path law f(L, R̄, Z)
- [ ] Learn f(L, R̄, Z) with PySR, radar along the path as reference (WAA removed first)

## 4. Rain estimation and evaluation
- [ ] Invert the combined model (A_rain + WAA) to get rain
- [ ] Normalized bias and normalized RMSE vs gauges and radar, split by link length
- [ ] Define and compute the false-rain-after-event metric
- [ ] Compare with power law, literature WAA models and the RNN baseline [11]
- [ ] Bootstrap stability of learned terms across links and events
- [ ] Final report
