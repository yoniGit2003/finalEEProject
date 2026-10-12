"""Get the fit caches of the research scripts, so that the long fits (about 9 hours) do not have to be repeated.

    python pipeline/get_caches.py          # download the zip from the GitHub release and unpack it
    python pipeline/get_caches.py --pack   # (maintainers) make a new zip from the caches on this computer

Why: the wet-antenna scripts that test the fitting method save every finished fit in research/.cache/wet_antenna/.
With these files a script reruns in seconds to minutes and reproduces the committed tables; without them it fits
everything again (run times: docs/REPRODUCE.md). The caches are not in git (about 400 MB), so they are published as
one zip attached to a GitHub release of this repository.

Download: first a plain download; if the repository is private that is refused, and the GitHub command line tool
is used instead (`gh auth login` once, https://cli.github.com). Files that already exist are not overwritten.

Publishing a new zip after a new fitting step (maintainers):
    python pipeline/get_caches.py --pack                       # writes the zip next to the repository, prints its sha256
    gh release create <TAG> <zip> --title "Fit caches <date>" --notes "Fit caches of research/wet_antenna"
then set TAG and SHA256 below to the new values and commit.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CACHE = REPO / "research" / ".cache"
GITHUB = "yoniGit2003/finalEEProject"
TAG = "fit-caches-2026-10-12"                               # the GitHub release that holds the zip
ASSET = "fit_caches.zip"
SHA256 = ""                                                 # of the published zip (printed by --pack); "" = not checked
# what goes into the zip: the folders with cached fits (one per script). Logs, PySR working files and loaded-data
# caches are left out: the scripts rebuild them.
FOLDERS = ["wet_antenna/check_per_tail_fit", "wet_antenna/check_pooled_fit", "wet_antenna/stress_pooled_fit",
           "wet_antenna/check_sindy", "wet_antenna/check_free_curve", "wet_antenna/confirm_free_curve",
           "wet_antenna/check_discovery_method", "wet_antenna/check_quantization", "wet_antenna/fit_real_drying"]


def sha256_of(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def pack():
    """Make the zip from the caches on this computer and say how to publish it."""
    out = REPO.parent / f"{REPO.name}_{ASSET}"
    files = [f for folder in FOLDERS for f in sorted((CACHE / folder).rglob("*"))
             if f.is_file() and f.suffix not in (".tmp", ".log")]
    if not files:
        raise SystemExit(f"no cached fits found under {CACHE}")
    print(f"packing {len(files)} files into {out} ...")
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for f in files:
            z.write(f, f.relative_to(CACHE).as_posix())
    print(f"{out}: {out.stat().st_size / 1e6:.0f} MB\nsha256: {sha256_of(out)}\n\nPublish it with:\n"
          f'  gh release create {TAG} "{out}#{ASSET}" --repo {GITHUB} --title "Fit caches" '
          f'--notes "Fit caches of research/wet_antenna (see pipeline/get_caches.py)"\n'
          "then put the tag and the sha256 into TAG and SHA256 at the top of pipeline/get_caches.py.")


def download(target):
    """Download the zip to `target`: plain HTTPS first, then the `gh` tool (needed for a private repository)."""
    url = f"https://github.com/{GITHUB}/releases/download/{TAG}/{ASSET}"
    part = target.with_suffix(".part")
    try:
        print(f"downloading {url}")
        with urllib.request.urlopen(url, timeout=60) as response, open(part, "wb") as out:
            shutil.copyfileobj(response, out, 1 << 20)
        part.replace(target)
        return
    except Exception as e:                                  # noqa: BLE001  (404 for a private repository, network errors)
        part.unlink(missing_ok=True)
        print(f"   plain download failed ({type(e).__name__}: {e}); trying the GitHub command line tool")
    if shutil.which("gh") is None:
        raise SystemExit("The repository is private or the download is blocked, and `gh` is not installed.\n"
                         "Install it from https://cli.github.com, run `gh auth login`, then run this script again.")
    subprocess.run(["gh", "release", "download", TAG, "--repo", GITHUB, "--pattern", ASSET, "--dir",
                    str(target.parent), "--clobber"], check=True)


def unpack(zip_path):
    """Unpack into research/.cache/; a file that already exists is kept (your own fits are never overwritten)."""
    new = kept = 0
    with zipfile.ZipFile(zip_path) as z:
        for name in z.namelist():
            if name.endswith("/"):
                continue
            out = CACHE / name
            if out.exists():
                kept += 1
                continue
            out.parent.mkdir(parents=True, exist_ok=True)
            with z.open(name) as src, open(out, "wb") as dst:
                shutil.copyfileobj(src, dst)
            new += 1
    print(f"unpacked {new} files into {CACHE} ({kept} already there, kept)")


def main():
    parser = argparse.ArgumentParser(description="Download (or pack) the fit caches of the research scripts.")
    parser.add_argument("--pack", action="store_true", help="make a new zip from the caches on this computer")
    args = parser.parse_args()
    if args.pack:
        return pack()
    CACHE.mkdir(parents=True, exist_ok=True)
    zip_path = CACHE / ASSET
    if not zip_path.exists():
        download(zip_path)
    if SHA256 and sha256_of(zip_path) != SHA256:
        raise SystemExit(f"{zip_path}: wrong sha256 (incomplete download?). Delete the file and run again.")
    unpack(zip_path)
    zip_path.unlink()


if __name__ == "__main__":
    sys.exit(main())
