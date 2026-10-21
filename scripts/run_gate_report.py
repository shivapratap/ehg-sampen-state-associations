#!/usr/bin/env python3
"""Apply the pre-registered decision rule; writes gate_summary.csv.

Thin wrapper. The implementation is `gi_sampen.gate`; all arguments are passed
through unchanged. Run from anywhere once the package is installed
(`python -m pip install -e .`).
"""
import sys

from gi_sampen.gate import main

if __name__ == "__main__":
    sys.exit(main())
