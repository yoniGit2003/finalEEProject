"""Data access.

CMLData (cmldata) reads a release folder: tables, Zarr time series, aligned per-link frames, event windows.
access.load_link_window and explore work on the work folder (data/work/openmrg) during development.
"""
from cmlpipe.api.cmldata import CMLData, CMLDataError  # noqa: F401
