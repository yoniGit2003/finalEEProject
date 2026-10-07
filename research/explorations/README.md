# explorations/: the playground

Try ideas here. Anything goes inside this folder: scripts, notebooks, CSVs, figures, next to each other. Nothing here
is reorganized by `/tidy`, and nothing here needs a CATALOG entry.

Light conventions (they make the next step easier):

- Start each script with a docstring: the question it asks, and what it writes.
- Read data only through `cmlpipe.api.CMLData` and the release in `research/config.yaml` (as in the modules).
- Put caches in `explorations/.cache/` (git-ignored), not next to the scripts.

**Graduating.** When a result matters for the report or for another module, it moves into a research module
(`research/<module>/`) and follows that module's layout: the script goes next to the module's code, its outputs go to
`outputs/`, it reads locations from the module's paths file, and every output gets a `CATALOG.md` entry. The old
exploration file stays here (or is deleted by you); modules never read files from `explorations/`.

Current contents: `dry_excursions.py` -> `dry_excursions.csv`; `plot_link_baseline.py` -> `link_10200_days.png`,
`link_10200_event.png`.
