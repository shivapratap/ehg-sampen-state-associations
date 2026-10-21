#!/usr/bin/env python3
"""AvgSampEn gain-sensitivity report.

Thin wrapper. The implementation is `gi_sampen.report_avg_gain`; all arguments are passed
through unchanged. Run from anywhere once the package is installed
(`python -m pip install -e .`).
"""
import sys

from gi_sampen.report_avg_gain import main

if __name__ == "__main__":
    sys.exit(main())
