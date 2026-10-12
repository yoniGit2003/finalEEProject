# explorations/: the playground

A place to try ideas without rules. Everything here except this README is **local only** (ignored by git), so it
never reaches GitHub and nothing else in the repository depends on it.

- Read data through `cmlpipe.api.CMLData`, like the modules do.
- When an idea works, move it into a module: copy the script into `research/<module>/`, take file locations from the
  module's `paths.py`, write to its `outputs/` with `research.common.save.save`, and add the script to
  `research/steps.yaml`. Modules never read files from `explorations/`.
- The scripts that were here before the October 2026 reorganization import the old module names
  (`baseline_estm`, `WAA_estm`); change the imports to `research.baseline` / `research.wet_antenna` to run them again.
