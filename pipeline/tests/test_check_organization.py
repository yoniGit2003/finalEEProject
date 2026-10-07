"""The repository follows ORGANIZATION.md (tools/check_organization.py), and the checker catches violations."""

import sys

from conftest import PIPELINE_ROOT

sys.path.insert(0, str(PIPELINE_ROOT / "tools"))
import check_organization as co  # noqa: E402


def test_repository_is_organized():
    assert co.check(co.tracked()) == []


def test_checker_catches_violations(tmp_path):
    good = {
        "research/m/README.md": "# m\n",
        "research/m/CATALOG.md": "`a.csv` and `c5_<id>.png`\n",
        "research/m/m_paths.py": 'OUT = "outputs"\n',
        "research/m/run.py": "import m_paths as P\n",
        "research/m/outputs/tables/a.csv": "x\n",
        "research/m/outputs/figures/c5_123.png": "",
    }
    bad = {
        "research/m/validation/README.md": "",                      # second README in a module
        "research/m/results.csv": "",                               # output outside outputs/
        "research/m/outputs/tables/b_v2.csv": "",                   # versioned name + no catalog entry
        "research/m/make.py": 'OUT = HERE / "outputs"\n',           # hard-coded output folder
        "research/m/archive/old/x.csv": "",                         # archive without NOTE.md
        "research/m/.gitignore": "",                                # nested .gitignore
        "notes.txt": "see [x](missing.md)\n",                       # not an allowed root file
        "docs/a.md": "see [x](missing.md)\n",                       # broken link
    }
    for name, text in {**good, **bad}.items():
        (tmp_path / name).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / name).write_text(text)
    assert co.check(sorted(good), tmp_path) == []
    found = "\n".join(co.check(sorted({**good, **bad}), tmp_path))
    for what in ("only one README", "results.csv: output file outside", "b_v2.csv: versioned name",
                 "b_v2.csv: no entry", "make.py:1: hard-coded", "missing NOTE.md", "nested .gitignore",
                 "notes.txt: not an allowed root file", "docs/a.md: broken link"):
        assert what in found, what
