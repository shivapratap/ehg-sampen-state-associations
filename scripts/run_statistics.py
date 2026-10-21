#!/usr/bin/env python3
"""Mixed models and BH-FDR within the GI family.

Thin wrapper. The implementation is `gi_sampen.statistics`; all arguments are passed
through unchanged. Run from anywhere once the package is installed
(`python -m pip install -e .`).
"""
import sys

from gi_sampen.statistics import main

if __name__ == "__main__":
    sys.exit(main())
