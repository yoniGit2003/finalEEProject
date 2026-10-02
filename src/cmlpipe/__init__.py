"""cmlpipe: a data pipeline for commercial microwave link (CML) research (OpenMRG, Gothenburg, JJA 2015).

The pipeline turns the raw CML, gauge, weather and radar files into documented, versioned products: L1 (common
format), L2 (QC flags), L3 (aggregates at 1 min - 1 h, radar along links), events / dry periods / split, an
inventory, and a read-only release with a Python access API. Definitions: docs/methods.md.

Subpackages
-----------
io         L1 schema + validators, raw-dataset adapters, ingest stages, the NetCDF reader
qc         flag bits (QCFlag), per-sample QC rules, qc stages
geometry   projections, radar grid, path / cell weights, gauge and neighbour matching
aggregate  L3 reductions and stages
events     rain events, dry periods, frozen split, per-link reference state
inventory  significance, wet runs, inventory tables, quality status
viz        figure style, quicklooks, link cards, inventory figures, run report
release    product metadata and the release builder
api        CMLData (read a release), work-folder helpers
Top level: config (YAML, stage hashes), provenance, parallel (worker processes), logging_utils, cli.
"""

__version__ = "0.2.0"
