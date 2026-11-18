"""Recompute all five certificates and compare them with the accompanying tables."""

import json
from pathlib import Path
import subprocess
import sys


CHECKS = Path(__file__).resolve().parent
CERTIFICATES = (
    ("derivative_seed.py", "derivative-seed.json", ("--cross-check",)),
    ("explicit_threshold.py", "explicit-threshold.json", ()),
    ("prime_bridge.py", "prime-bridge.json", ()),
    ("initial_intervals.py", "initial-intervals.json", ()),
    ("interval_replay.py", "interval-replay.json", ()),
)


def main():
    for script, table, options in CERTIFICATES:
        actual = subprocess.check_output(
            [sys.executable, "-B", script, *options], cwd=CHECKS, text=True
        )
        if json.loads(actual) != json.loads((CHECKS / table).read_text()):
            raise SystemExit(f"Recomputed output differs from {table}")
        print(f"Verified {table}", flush=True)


if __name__ == "__main__":
    main()
