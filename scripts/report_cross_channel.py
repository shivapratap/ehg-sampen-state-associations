#!/usr/bin/env python3
"""Cross-channel replication report (Table 1 source).

Thin wrapper. The implementation is `gi_sampen.report_cross_channel`; all arguments are passed
through unchanged. Run from anywhere once the package is installed
(`python -m pip install -e .`).
"""
import sys

from gi_sampen.report_cross_channel import main

if __name__ == "__main__":
    sys.exit(main())
