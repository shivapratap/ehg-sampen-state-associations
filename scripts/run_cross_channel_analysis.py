#!/usr/bin/env python3
"""EHG9-EHG12 cross-channel GI feature extraction.

Thin wrapper. The implementation is `gi_sampen.cross_channel`; all arguments are passed
through unchanged. Run from anywhere once the package is installed
(`python -m pip install -e .`).
"""
import sys

from gi_sampen.cross_channel import main

if __name__ == "__main__":
    sys.exit(main())
