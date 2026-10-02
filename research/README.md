# research/

Research code goes here (baselines, WAA, models, notebooks, figures, thesis). It may only use the data through
`cmlpipe.api.CMLData` and a pinned release version (`research/config.yaml`: `release: openmrg_v1`). It must not
import cmlpipe internals and cmlpipe must never import research.

`cmlpipe.viz` (plotting helpers) may also be imported. The rule is checked by `python tools/check_boundaries.py`
(also run by the test suite).

```python
import yaml
from cmlpipe.api import CMLData
cfg = yaml.safe_load(open("research/config.yaml"))
d = CMLData.open(f"data/release/{cfg['release']}")
```
