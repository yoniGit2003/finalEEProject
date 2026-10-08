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


def test_repository_outputs_follow_the_levels():
    assert co.check_outputs(co.tracked()) == []


def test_output_checks_catch_violations(tmp_path):
    m = tmp_path / "research" / "m"
    files = {
        "m_paths.py": "",
        "outputs/MANIFEST.csv": "path,level,purpose,script,commit,date\n"
                                "outputs/figures/a.png,diag,x,s.py,abc,2026-01-01\n"
                                "outputs/tables/b.csv,bogus,x,s.py,abc,2026-01-01\n"
                                + "".join(f"outputs/figures/k{i}.png,key,x,s.py,abc,2026-01-01\n" for i in range(7)),
        "outputs/TIDY_LOG.md": "",
        "outputs/figures/a.png": "",
        "outputs/tables/b.csv": "",
        "outputs/tables/c_v2.csv": "",                              # not in the MANIFEST + versioned name
    }
    for name, text in files.items():
        (m / name).parent.mkdir(parents=True, exist_ok=True)
        (m / name).write_text(text)
    found = "\n".join(co.check_outputs(["research/m/outputs/figures/a.png"], tmp_path))
    for what in ("c_v2.csv: not in outputs/MANIFEST.csv", "c_v2.csv: versioned name", "b.csv: level 'bogus'",
                 "7 key outputs", "a.png: level diag but tracked"):
        assert what in found, what
    assert "TIDY_LOG.md" not in found and "MANIFEST.csv: not in" not in found
