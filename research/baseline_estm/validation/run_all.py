"""Run the whole baseline validation (R1) in order; runtimes go to out/runtimes.csv.

python research/baseline_estm/validation/run_all.py          (sample from research/config.yaml, normally FULL)
BV_SAMPLE=DEV python research/baseline_estm/validation/run_all.py   (quick run on the 24 DEV sub-links -> out_dev/)
"""
import c1_fake_events, c2_after_window, c3_daily_cycle, c4_jumps, c5_flagged, c6_sensitivity, r1b_event_level
from common import Timer

if __name__ == "__main__":
    with Timer("total"):
        for m in (c1_fake_events, c2_after_window, c3_daily_cycle, c4_jumps, c5_flagged, c6_sensitivity,
                  r1b_event_level):
            m.main()
