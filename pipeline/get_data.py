"""Get the data: download the raw OpenMRG dataset from Zenodo, then build the data release with the pipeline.

    python pipeline/get_data.py                  # download (318 MB) + build everything with 4 processes
    python pipeline/get_data.py --download-only  # only the raw data
    python pipeline/get_data.py --workers 2      # fewer processes for the build

What it does, in order:
  1. If data/raw/openmrg/ already holds the raw files, the download is skipped (the raw data is never changed).
     Otherwise OpenMRG.zip (v1.1, https://doi.org/10.5281/zenodo.7107689, CC BY-SA 4.0) is downloaded, its md5 is
     compared with the one Zenodo publishes, and it is unpacked into data/raw/openmrg/.
  2. `python -m cmlpipe run all` builds the work files (data/work/openmrg/) and the release
     (data/release/openmrg_v1/). The build is long (the recorded stages add up to more than half an hour; it was
     never timed in one go). Stages that are already up to date are skipped, so the command can be repeated.

If Zenodo refuses the download (it sometimes blocks a network for a while), download OpenMRG.zip in a browser from
the DOI above, put it in data/raw/ and run this script again: it will use that file.
Only the Python standard library is used for the download; the build needs the conda environment (see README.md).
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
import time
import urllib.request
import zipfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw" / "openmrg"
ZIP = REPO / "data" / "raw" / "OpenMRG.zip"
URL = "https://zenodo.org/api/records/7107689/files/OpenMRG.zip/content"
MD5 = "4e4cdf6ced4d45a084c348f245b8d690"                   # published by Zenodo for OpenMRG.zip (317 994 487 bytes)
# the files the pipeline reads (pipeline/configs/openmrg.yaml), relative to data/raw/openmrg/
NEEDED = ["cml/cml.nc", "cml/cml_metadata.csv", "gauges/city/CityGauges-2015JJA.csv",
          "gauges/city/CityGauges-metadata.csv", "gauges/smhi/GbgA-71420-2015JJA.csv",
          "gauges/smhi/GbgA-71420-metadata.csv", "radar/radar.nc"]


def missing_raw(folder=None):
    """The needed raw files that are not in `folder` (default: data/raw/openmrg)."""
    folder = folder or RAW
    return [f for f in NEEDED if not (folder / f).exists()]


def md5_of(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def download(target=None):
    """Download OpenMRG.zip to `target` (through a .part file, so an interrupted download is never used)."""
    target = target or ZIP
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_suffix(".part")
    print(f"downloading {URL}\n   -> {target}")
    request = urllib.request.Request(URL, headers={"User-Agent": "finalEEProject get_data.py"})
    with urllib.request.urlopen(request, timeout=60) as response, open(part, "wb") as out:
        total = int(response.headers.get("Content-Length") or 0)
        done = 0
        while True:
            block = response.read(1 << 20)
            if not block:
                break
            out.write(block)
            done += len(block)
            if total:
                print(f"\r   {done / 1e6:6.0f} / {total / 1e6:.0f} MB", end="", flush=True)
    print()
    part.replace(target)


def unpack(zip_path=None, folder=None):
    """Unpack the archive into `folder`. The archive keeps everything under one top folder; that folder is dropped,
    so that cml/cml.nc ends up directly in `folder`."""
    zip_path, folder = zip_path or ZIP, folder or RAW
    with zipfile.ZipFile(zip_path) as z:
        names = [n for n in z.namelist() if not n.endswith("/")]
        anchor = next((n for n in names if n.endswith("cml/cml.nc")), None)
        if anchor is None:
            raise SystemExit(f"{zip_path} does not look like OpenMRG.zip (no cml/cml.nc inside)")
        prefix = anchor[: -len("cml/cml.nc")]               # e.g. "OpenMRG/"
        print(f"unpacking {len(names)} files into {folder}")
        for name in names:
            if not name.startswith(prefix):
                continue
            out = folder / name[len(prefix):]
            out.parent.mkdir(parents=True, exist_ok=True)
            with z.open(name) as src, open(out, "wb") as dst:
                shutil.copyfileobj(src, dst)


def get_raw():
    """Step 1: make sure data/raw/openmrg/ holds the raw files."""
    if not missing_raw():
        print(f"raw data found in {RAW}: download skipped")
        return
    if RAW.exists() and any(RAW.iterdir()):
        raise SystemExit(f"{RAW} exists but these files are missing: {missing_raw()}\n"
                         "The raw folder is never changed by this script: fix or remove it by hand, then run again.")
    if not ZIP.exists():
        error = None
        for attempt in (1, 2, 3):                           # Zenodo sometimes answers with a time-out: try again
            try:
                download()
                error = None
                break
            except Exception as e:                          # noqa: BLE001  (network errors of any kind)
                print(f"   attempt {attempt} failed: {type(e).__name__}: {e}")
                error = e
                time.sleep(20)
        if error is not None:
            raise SystemExit(f"download failed: {type(error).__name__}: {error}\nDownload OpenMRG.zip in a browser "
                             f"from https://doi.org/10.5281/zenodo.7107689, put it at {ZIP} and run this script again.")
    found = md5_of(ZIP)
    if found != MD5:
        raise SystemExit(f"{ZIP}: md5 is {found}, expected {MD5} (incomplete download?). Delete the file and retry.")
    print("md5 ok")
    unpack()
    if missing_raw():
        raise SystemExit(f"after unpacking, these files are still missing: {missing_raw()}")
    ZIP.unlink()                                            # the unpacked files are what the pipeline reads
    print(f"raw data ready in {RAW}")


def build(workers):
    """Step 2: run the whole pipeline (work files and the release). Up-to-date stages are skipped."""
    cmd = [sys.executable, "-m", "cmlpipe", "--workers", str(workers), "run", "all"]
    print("running:", " ".join(cmd))
    subprocess.run(cmd, cwd=REPO, check=True)


def main():
    parser = argparse.ArgumentParser(description="Download raw OpenMRG from Zenodo and build the data release.")
    parser.add_argument("--download-only", action="store_true", help="only download and unpack the raw data")
    parser.add_argument("--workers", type=int, default=4, help="processes for the build (default 4)")
    args = parser.parse_args()
    get_raw()
    if not args.download_only:
        build(args.workers)
        print("done: the release is in data/release/ (open it with cmlpipe.api.CMLData.open)")


if __name__ == "__main__":
    main()
