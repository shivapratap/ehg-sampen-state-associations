#!/usr/bin/env python3
"""Standard finite-resolution parity check.

Recomputes the six profile summaries on real analysis windows and compares
them against the cached values produced by the prior study's feature engine.
Defaults to 200 windows, the number reported in the manuscript; pass
``--parity N`` to change it.

Thin wrapper around ``gi_sampen.extract``.
"""
import sys

from gi_sampen.extract import main

if __name__ == "__main__":
    if not any(arg.startswith("--parity") for arg in sys.argv[1:]):
        sys.argv.extend(["--parity", "200"])
    sys.exit(main())
