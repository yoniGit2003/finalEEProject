"""Quality control (L2): flag bits (flags.QCFlag), the per-sample rules (rules) and the qc stages (stage).

`from cmlpipe import qc` gives the rule functions directly (qc.flag_spikes, ...), as before the reorganisation."""
from cmlpipe.qc.flags import FLAGS, QCFlag, decode_flags, flag_attrs, mask_of  # noqa: F401
from cmlpipe.qc.rules import *  # noqa: F401,F403
