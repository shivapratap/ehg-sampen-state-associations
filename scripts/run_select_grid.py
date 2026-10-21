#!/usr/bin/env python3
"""Stage A1: fix the dimensionless tolerance grid from baseline windows.

Thin wrapper. The implementation is `gi_sampen.fixed_grid`; all arguments are passed
through unchanged. Run from anywhere once the package is installed
(`python -m pip install -e .`).
"""
import sys

from gi_sampen.fixed_grid import main

if __name__ == "__main__":
    sys.exit(main())
