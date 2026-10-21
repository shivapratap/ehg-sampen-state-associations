#!/usr/bin/env python3
"""Amplitude collinearity audit.

Thin wrapper. The implementation is `gi_sampen.amplitude`; all arguments are passed
through unchanged. Run from anywhere once the package is installed
(`python -m pip install -e .`).
"""
import sys

from gi_sampen.amplitude import main

if __name__ == "__main__":
    sys.exit(main())
