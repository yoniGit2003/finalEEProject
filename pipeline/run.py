"""DEPRECATED shim: use `python -m cmlpipe run <stage>` (see src/cmlpipe/cli.py).

    python pipeline/run.py --config configs/openmrg.yaml --steps events,split   # still works, prints a note
"""
import sys
from pathlib import Path

if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # works without `pip install -e .`
    from cmlpipe.cli import main
    print("note: pipeline/run.py is deprecated; use `python -m cmlpipe run <stage> [--all --force --workers N]`",
          file=sys.stderr)
    sys.exit(main(sys.argv[1:]))
