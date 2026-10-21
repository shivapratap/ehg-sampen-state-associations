#!/usr/bin/env python3
"""GI resolution ablation across the decimal ladder.

Thin wrapper. The implementation is `gi_sampen.resolution`; all arguments are passed
through unchanged. Run from anywhere once the package is installed
(`python -m pip install -e .`).
"""
import sys

from gi_sampen.resolution import main

if __name__ == "__main__":
    sys.exit(main())
