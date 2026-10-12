"""research/steps.yaml (the list `python -m research` runs) matches the scripts, and docs/REPRODUCE.md is up to date.

Nothing from research/ is imported here (the pipeline never imports research): the files are read as text and
`python -m research docs --check` is run as a separate process.
"""

import ast
import subprocess
import sys

import yaml

from conftest import PIPELINE_ROOT

REPO = PIPELINE_ROOT.parent
RESEARCH = REPO / "research"
CONFIG = yaml.safe_load((RESEARCH / "steps.yaml").read_text(encoding="utf-8"))
STEPS = CONFIG["steps"]


def module_file(module):
    return REPO / (module.replace(".", "/") + ".py")


def entry_scripts():
    """Every research script that can be run: it has a main() function and an `if __name__ == "__main__"` block."""
    out = []
    for folder in ("baseline", "wet_antenna"):
        for f in sorted((RESEARCH / folder).rglob("*.py")):
            tree = ast.parse(f.read_text(encoding="utf-8"))
            has_main = any(isinstance(n, ast.FunctionDef) and n.name == "main" for n in tree.body)
            has_guard = any(isinstance(n, ast.If) and "__main__" in ast.unparse(n.test) for n in tree.body)
            if has_main and has_guard:
                out.append(".".join(f.relative_to(REPO).with_suffix("").parts))
    return out


def test_step_names_are_unique_and_complete():
    names = [s["name"] for s in STEPS]
    assert len(names) == len(set(names))
    for s in STEPS:
        assert s.get("module") and s.get("what") and s.get("time"), s["name"]


def test_every_listed_command_exists():
    for s in STEPS:
        assert module_file(s["module"]).exists(), f"{s['name']}: {s['module']} not found"
    for module in CONFIG.get("run_by_other_steps", {}):
        assert module_file(module).exists(), f"run_by_other_steps: {module} not found"


def test_every_entry_script_is_listed():
    listed = {s["module"] for s in STEPS} | set(CONFIG.get("run_by_other_steps", {}))
    missing = [m for m in entry_scripts() if m not in listed]
    assert not missing, f"add these scripts to research/steps.yaml: {missing}"


def test_reproduce_page_is_up_to_date():
    r = subprocess.run([sys.executable, "-m", "research", "docs", "--check"], cwd=REPO, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
