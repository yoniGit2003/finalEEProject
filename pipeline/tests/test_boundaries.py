"""The data layer / research boundary (tools/check_boundaries.py) holds, and the checker catches violations."""

import sys

from conftest import PIPELINE_ROOT

sys.path.insert(0, str(PIPELINE_ROOT / "tools"))
import check_boundaries as cb  # noqa: E402


def test_repository_respects_the_boundary():
    assert cb.violations() == []


def test_checker_catches_violations(tmp_path):
    (tmp_path / "research").mkdir()
    (tmp_path / "pipeline" / "src").mkdir(parents=True)
    (tmp_path / "research" / "ok.py").write_text("from cmlpipe.api import CMLData\nfrom cmlpipe import viz\n")
    (tmp_path / "research" / "bad.py").write_text("from cmlpipe.qc import flags\nimport cmlpipe\n")
    (tmp_path / "pipeline" / "src" / "bad.py").write_text("import research.models\n")
    bad = cb.violations(tmp_path)
    assert len(bad) == 3 and all("ok.py" not in b for b in bad)
