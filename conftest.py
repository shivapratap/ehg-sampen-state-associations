"""Make the src/ layout importable without requiring an editable install.

Tests still work if the package IS installed (`pip install -e .`); this is
only a fallback so `pytest` runs out of the box on a fresh checkout.
"""
import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parent / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))
