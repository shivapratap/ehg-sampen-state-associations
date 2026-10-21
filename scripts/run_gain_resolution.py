#!/usr/bin/env python3
"""Quantization-resolution sensitivity under gain.

Thin wrapper. The implementation is `gi_sampen.gain_resolution`; all arguments are passed
through unchanged. Run from anywhere once the package is installed
(`python -m pip install -e .`).
"""
import sys

from gi_sampen.gain_resolution import main

if __name__ == "__main__":
    sys.exit(main())
