"""Trivial cross-cutting hashing utilities shared across pipeline stages.

Kept deliberately tiny and stage-agnostic -- everything that is actual
pipeline *logic* (auditing, windowing, feature extraction, etc.) lives in
the numbered stage scripts themselves, not here.
"""
from __future__ import annotations

import hashlib
from pathlib import Path


def sha256_of_file(path: Path, chunk_size: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(chunk_size), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_of_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def hash_directory_listing(directory: Path) -> str:
    """Hash of (filename, size) pairs for every file directly under
    `directory`, sorted -- a cheap way to detect "the directory changed"
    without hashing every byte of a multi-GB raw dataset on every run.
    """
    entries = sorted(
        (p.name, p.stat().st_size) for p in directory.iterdir() if p.is_file()
    )
    payload = "\n".join(f"{name}:{size}" for name, size in entries).encode("utf-8")
    return sha256_of_bytes(payload)
