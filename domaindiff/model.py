"""
Module: model
Purpose: Internal hit record and configuration (defaults documented in README).
Author: Dr. Arli Aditya Parikesit
Date: 2026
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass, field, fields
from pathlib import Path


class InputError(Exception):
    """Raised for any problem with user-supplied files. Message is user-facing."""


@dataclass(frozen=True)
class Hit:
    """One domain hit. Coordinates are 1-based, inclusive, on the protein."""

    start: int
    end: int
    acc: str
    evalue: float | None = None


# Run: protein ID -> list of hits on that protein.
Run = dict[str, list[Hit]]


@dataclass(frozen=True)
class Config:
    boundary_tolerance: int = 5
    evalue_cutoff: float = 1e-5
    cooccurrence_log2_threshold: float = 1.0
    min_protein_overlap: float = 0.5
    overlap_rule: str = "keep"  # "keep" | "best_evalue"
    interproscan_analyses: list[str] = field(default_factory=lambda: ["Pfam"])
    # Analyses whose InterProScan column 9 is an E-value; others (e.g. ProSiteProfiles, HAMAP) hold scores.
    interproscan_evalue_analyses: list[str] = field(default_factory=lambda: ["Pfam"])
    strip_accession_version: bool = True
    max_file_mb: int = 2048
    max_rows: int = 10_000_000

    def __post_init__(self) -> None:
        if self.overlap_rule not in ("keep", "best_evalue"):
            raise InputError("config: overlap_rule must be 'keep' or 'best_evalue'")
        if self.boundary_tolerance < 0 or self.max_file_mb <= 0 or self.max_rows <= 0:
            raise InputError("config: tolerance must be >= 0, limits must be > 0")
        if not 0 <= self.min_protein_overlap <= 1:
            raise InputError("config: min_protein_overlap must be in [0, 1]")
        if self.evalue_cutoff < 0 or self.cooccurrence_log2_threshold < 0:
            raise InputError("config: evalue_cutoff and log2 threshold must be >= 0")


def load_config(path: Path | None) -> Config:
    if path is None:
        return Config()
    try:
        data = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as e:
        raise InputError(f"config {path}: {e}") from e
    known = {f.name: f for f in fields(Config)}
    unknown = set(data) - set(known)
    if unknown:
        raise InputError(f"config: unknown keys {sorted(unknown)}")
    defaults = Config()
    for k, v in data.items():
        expected = type(getattr(defaults, k))
        ok = isinstance(v, expected) or (expected is float and isinstance(v, int))
        if not ok or isinstance(v, bool) != (expected is bool):
            raise InputError(f"config: {k} must be {expected.__name__}")
        if expected is list and not all(isinstance(x, str) for x in v):
            raise InputError(f"config: {k} must be a list of strings")
    return Config(**data)
