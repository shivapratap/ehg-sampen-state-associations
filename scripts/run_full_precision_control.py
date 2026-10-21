#!/usr/bin/env python3
"""Unrounded / full-precision control (mechanism isolation).

Thin wrapper. The implementation is `gi_sampen.full_precision`; all arguments are passed
through unchanged. Run from anywhere once the package is installed
(`python -m pip install -e .`).
"""
import sys

from gi_sampen.full_precision import main

if __name__ == "__main__":
    sys.exit(main())
