#!/usr/bin/env python3
"""Stage A2: GI feature extraction over the analysis units.

Thin wrapper. The implementation is `gi_sampen.extract`; all arguments are passed
through unchanged. Run from anywhere once the package is installed
(`python -m pip install -e .`).
"""
import sys

from gi_sampen.extract import main

if __name__ == "__main__":
    sys.exit(main())
