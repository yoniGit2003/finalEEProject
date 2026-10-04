# Research roadmap

What is done and what is left in the research (research/), following the work scope of `docs/proposal.md`
(section 6). References [n] are those of the proposal. Tick an item by changing `- [ ]` to `- [x]`: the status block
of `research/README.md` and the one-line summary in the root `README.md` are generated from this file on every commit
(`pipeline/tools/update_readme.py`), so this is the only place where research progress is edited by hand.

Status: in progress (not started in code yet)

## 2. Simulator (known ground truth)
- [ ] Rain fields with controlled cell size
- [ ] Exact path integral of a·R(x,t)^b along each link
- [ ] Known WAA model (wetting and drying)
- [ ] Noise and quantization (RSL 0.3 dB, TSL 1 dB)

## 3. Stage 1: WAA
- [ ] Find links with a drying tail after rain ends
- [ ] Fit literature WAA models [1], [9] per link
- [ ] Learn a dynamic WAA model with weak-form SINDy (rain as input), first in simulation
- [ ] Compare the learned models across links (same form, different coefficients?)
- [ ] Drying time from SINDy vs drying time measured directly from tails

## 4. Stage 2: path law f(L, R̄, Z)
- [ ] Simulation: exact integral vs the linear law a·R̄^b·L
- [ ] Learn f(L, R̄, Z) with PySR in simulation
- [ ] Real data with radar along the path as reference (WAA removed first)

## 5. Rain estimation and evaluation
- [ ] Invert the combined model (A_rain + WAA) to get rain
- [ ] Normalized bias and normalized RMSE vs gauges and radar, split by link length
- [ ] Define and compute the false-rain-after-event metric
- [ ] Compare with power law, literature WAA models and the RNN baseline [11]
- [ ] Bootstrap stability of learned terms across links and events
- [ ] Final report
