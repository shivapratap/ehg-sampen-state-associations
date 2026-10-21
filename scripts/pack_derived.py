#!/usr/bin/env python3
"""MAINTAINER TOOL: pack the derived intermediate tables into data/derived/ehg_derived_inputs.zip.

The zip holds everything the statistical stages need so that nobody has to download the raw
recordings or run the 30-minute-to-2-hour extraction stages. Contents are exactly the files
listed in DERIVED below (paths relative to the repository root). A SHA-256 of the zip and of
every member is written to data/derived/MANIFEST.sha256.
"""
import hashlib, sys, zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from restore_derived import DERIVED_DIRS, DERIVED_FILES  # noqa: E402

OUT = ROOT / "data" / "derived"
ZIP = OUT / "ehg_derived_inputs.zip"
sha = lambda b: hashlib.sha256(b).hexdigest()

members = []
for d in DERIVED_DIRS:
    members += sorted(p for p in (ROOT / d).rglob("*") if p.is_file() and "__pycache__" not in p.parts)
members += [ROOT / f for f in DERIVED_FILES]
with zipfile.ZipFile(ZIP, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    for p in members:
        z.write(p, p.relative_to(ROOT).as_posix())
lines = [f"{sha(p.read_bytes())}  {p.relative_to(ROOT).as_posix()}" for p in members]
(OUT / "MANIFEST.sha256").write_text(f"{sha(ZIP.read_bytes())}  {ZIP.name}\n" + "\n".join(lines) + "\n")
print(f"{len(members)} files -> {ZIP} ({ZIP.stat().st_size/1e6:.1f} MB)")
