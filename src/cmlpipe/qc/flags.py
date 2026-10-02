"""The L2 QC flag bits: defined ONCE here (QCFlag) and imported everywhere else.

One uint16 bitmask per 10 s sample, dims (time, sublink_id, cml_id) = the L1 dims. The file documents itself with
the CF attributes flag_masks / flag_meanings. Bits 8-15 are reserved: level_shift and pair_mismatch are deliberately
NOT flags (not calibrated), so they must stay 0. The rules and thresholds behind each bit: cmlpipe.qc.rules and
docs/methods.md.
"""

from __future__ import annotations

from enum import IntFlag

import numpy as np


class QCFlag(IntFlag):
    """QC flag bits of the L2 product (value = 1 << bit)."""

    missing_rsl = 1 << 0    # rsl is NaN
    missing_tsl = 1 << 1    # tsl is NaN
    out_of_range = 1 << 2   # rsl or tsl finite but outside the configured physical range
    rx_floor = 1 << 3       # receiver clipping: rsl held at the sub-link's floor level
    spike = 1 << 4          # isolated short outlier in rsl
    tsl_step = 1 << 5       # tsl changed by >= tsl_step min_db versus the previous valid sample
    flatline = 1 << 6       # rsl exactly constant for >= flatline min_duration
    network_gap = 1 << 7    # sample is NaN at a time when many alive sub-links are NaN (collector outage)

    @property
    def bit(self) -> int:
        """Bit position (0-7)."""
        return int(self).bit_length() - 1


# name -> bit position, in bit order (the form the code and the file attributes use)
FLAGS = {f.name: f.bit for f in QCFlag}
FLAG_RESERVED_BITS = tuple(range(8, 16))
FLAGS_VAR = "flags"


def mask_of(names) -> np.uint16:
    """uint16 mask with the bits of the named flags set (e.g. a blocking set)."""
    unknown = [n for n in names if n not in FLAGS]
    if unknown:
        raise ValueError(f"unknown flag names {unknown}; known: {list(FLAGS)}")
    return np.uint16(sum(int(QCFlag[n]) for n in names))


def decode_flags(flags) -> dict:
    """uint16 bitmask -> {flag name: boolean array}."""
    flags = np.asarray(flags)
    return {f.name: (flags & np.uint16(int(f))) != 0 for f in QCFlag}


def flag_attrs() -> dict:
    """CF attributes of the flags variable (flag_masks / flag_meanings)."""
    return {
        "long_name": "CML QC flags (bitmask)",
        "flag_masks": np.array([int(f) for f in QCFlag], dtype="uint16"),
        "flag_meanings": " ".join(FLAGS),
        "comment": "bits 8-15 reserved (level_shift and pair_mismatch are deliberately not flags yet); "
                   "decode with cmlpipe.qc.flags.decode_flags",
    }
