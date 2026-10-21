#!/usr/bin/env python3
"""Positive-gain stress test / GI invariance check.

Thin wrapper. The implementation is `gi_sampen.gain_stress`; all arguments are passed
through unchanged. Run from anywhere once the package is installed
(`python -m pip install -e .`).
"""
import sys

from gi_sampen.gain_stress import main

if __name__ == "__main__":
    sys.exit(main())
