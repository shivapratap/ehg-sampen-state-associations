#!/usr/bin/env python3
"""Build Figure 1 and its source CSVs.

Thin wrapper. The implementation is `gi_sampen.figure1`; all arguments are passed
through unchanged. Run from anywhere once the package is installed
(`python -m pip install -e .`).
"""
import sys

from gi_sampen.figure1 import main

if __name__ == "__main__":
    sys.exit(main())
