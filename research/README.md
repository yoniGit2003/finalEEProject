# research/: the research code

The question, in the proposal's notation: the measured attenuation is A(t) = A_rain(t) + δ(t), where A_rain is the
loss caused by rain along the path and δ(t) collects everything else, mainly wet-antenna attenuation (WAA). The
common model is A_rain ≈ a·R̄^b·L, exact only for uniform rain. We look for (1) a dynamic model of WAA and (2) the
real dependence on the link, A_rain = a·R̄^b·f(L, R̄, Z), where uniform rain gives f = L. WAA must be removed before
f is learned, because a constant WAA looks like f = L + L0. Background: [`../docs/proposal.md`](../docs/proposal.md);
results so far: [`../docs/research_log.pdf`](../docs/research_log.pdf); plan: [`../docs/ROADMAP.md`](../docs/ROADMAP.md).

## How to run

From the repository root (setup and data: [`../docs/REPRODUCE.md`](../docs/REPRODUCE.md)):

```
python -m research                    # the menu: every step in order, what it does, how long it takes
python -m research run <step>         # one step;  `run all` = everything in order
python -m research.wet_antenna.extract_tails      # any script can also be run directly
```

## Layout

```
config.yaml      every research setting (one section per module) and the pinned data release
steps.yaml       every research step in run order (the list `python -m research` runs)
__main__.py      the `python -m research` command
common/          small shared helpers: save.py (save a figure / table / text), report.py, plotting.py, config.py
baseline/        the dry-weather baseline of every sub-link and A(t)                 (done: method v2 frozen)
wet_antenna/     Stage 1: drying tails (frozen tail set v1) and the drying law        (in progress)
explorations/    playground: local only, not in git
.cache/          saved fits and logs: not in git (python pipeline/get_caches.py downloads ours)
```

Planned modules: `path_law/` (Stage 2: f(L, R̄, Z) with PySR) and `evaluation/` (rain estimates, metrics).

## A module

```
<module>/
  README.md     what it computes, its status, how to run it, what each file does
  METHOD.md     the method, the decisions of each step with their numbers, open issues
  paths.py      every file location of the module (scripts never build paths themselves)
  *.py          library modules and entry scripts (each entry script has a main() and is listed in steps.yaml)
  outputs/      products/ (data other steps read), tables/, figures/, reports/ (one report per script)
```

Everything in `outputs/` is committed, except the few big files named in the root `.gitignore`.

## Rules

- **Data only through the release.** Research code reads data with `cmlpipe.api.CMLData` and the release pinned in
  `config.yaml`; it may import only `cmlpipe.api` and `cmlpipe.viz` (checked by
  `python pipeline/tools/check_boundaries.py`).
- **Settings in `config.yaml`, locations in `paths.py`**, never scattered in scripts. Decision rules are written
  into the config before a run.
- **Never select links or tails by how well a model fits them** (that would assume the answer).
- Notation as in the proposal: A(t), A_rain, δ(t), R̄, L, a, b, f(L, R̄, Z). References [n] are the proposal's.
- Figures must not rely on colour alone: use line styles, hatching and labels.
- A new script: write outputs with `from research.common.save import save`, give it a `main()`, add it to
  `steps.yaml`, then run `python -m research docs`.
