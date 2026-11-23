"""
Module: loaders
Purpose: Parse InterProScan TSV and HMMER --domtblout into one internal Run table.
Author: Dr. Arli Aditya Parikesit
Date: 2026
Parameters:
    - fmt: "interproscan" | "hmmscan" | "hmmsearch"
References:
    - InterProScan 5.78-109.0 / 6.0.3.0 TSV output docs (ebi-pf-team/interproscan-docs, interproscan6)
    - HMMER 3.4 user guide, tabular output (EddyRivasLab/hmmer, documentation/userguide/tabular.tex)
"""

from __future__ import annotations

import math
import re
import stat
from collections.abc import Iterator
from pathlib import Path

from .model import Config, Hit, InputError, Run

FORMATS = ("interproscan", "hmmscan", "hmmsearch")
MAX_FIELD_LEN = 1000
# Version suffix on HMM accessions, e.g. PF00069.29 or TIGR00001.1. Accessions with
# other punctuation (G3DSA:3.30.200.20) never match, so they are left untouched.
_VERSIONED = re.compile(r"^([A-Za-z]+\d+)\.\d+$")
_PRINTABLE = re.compile(r"^[\x21-\x7e]+$")  # visible ASCII, no spaces or control chars


def _lines(path: Path, cfg: Config) -> Iterator[tuple[int, str]]:
    """Yield (line number, line) after size, type and encoding checks."""
    try:
        st = path.stat()
    except OSError as e:
        raise InputError(f"{path}: {e.strerror}") from e
    if not stat.S_ISREG(st.st_mode):
        raise InputError(f"{path}: not a regular file")
    if st.st_size > cfg.max_file_mb * 1024 * 1024:
        raise InputError(f"{path}: larger than max_file_mb={cfg.max_file_mb}")
    with path.open("rb") as fh:
        if fh.read(2) == b"\x1f\x8b":
            raise InputError(f"{path}: gzip input is not supported; decompress first")
        fh.seek(0)
        for n, raw in enumerate(fh, 1):
            if n > cfg.max_rows:
                raise InputError(f"{path}: more than max_rows={cfg.max_rows} lines")
            try:
                yield n, raw.decode("utf-8").rstrip("\r\n")
            except UnicodeDecodeError as e:
                raise InputError(f"{path}:{n}: not valid UTF-8") from e


def _ident(value: str, what: str, where: str) -> str:
    if len(value) > MAX_FIELD_LEN or not _PRINTABLE.match(value):
        raise InputError(f"{where}: {what} must be 1-{MAX_FIELD_LEN} printable non-space ASCII chars")
    return value


def _coords(s: str, e: str, where: str) -> tuple[int, int]:
    try:
        start, end = int(s), int(e)
    except ValueError as err:
        raise InputError(f"{where}: non-integer coordinates {s!r}, {e!r}") from err
    if start < 1 or end < start:
        raise InputError(f"{where}: invalid coordinates start={start} end={end}")
    return start, end


def _evalue(s: str, where: str) -> float | None:
    if s == "-":
        return None
    try:
        v = float(s)
    except ValueError as err:
        raise InputError(f"{where}: non-numeric score {s!r}") from err
    if math.isnan(v) or v < 0:
        raise InputError(f"{where}: invalid score {s!r}")
    return v


def _norm_acc(acc: str, cfg: Config) -> str:
    if cfg.strip_accession_version:
        m = _VERSIONED.match(acc)
        if m:
            return m.group(1)
    return acc


def load(path: Path, fmt: str, cfg: Config, mapping: dict[str, str] | None = None) -> Run:
    """Read one annotation file. `mapping` (old -> new accession) is applied after version stripping."""
    if fmt not in FORMATS:
        raise InputError(f"unknown format {fmt!r}; choose from {FORMATS}")
    run: Run = {}
    analyses = set(cfg.interproscan_analyses)
    skipped: set[str] = set()
    for n, line in _lines(path, cfg):
        where = f"{path}:{n}"
        if not line.strip() or (fmt != "interproscan" and line.startswith("#")):
            continue
        if fmt == "interproscan":
            cols = line.split("\t")
            if len(cols) not in (11, 13, 14, 15):
                raise InputError(f"{where}: expected 11-15 tab-separated columns, got {len(cols)}")
            if "all" not in analyses and cols[3] not in analyses:
                skipped.add(cols[3][:50])
                continue
            protein, acc = cols[0], cols[4]
            start, end = _coords(cols[6], cols[7], where)
            ev = _evalue(cols[8], where) if cols[3] in cfg.interproscan_evalue_analyses else None
        else:
            cols = line.split(maxsplit=22)
            if len(cols) < 22:
                raise InputError(f"{where}: expected >= 22 whitespace-separated fields, got {len(cols)}")
            if fmt == "hmmscan":  # target = profile, query = protein sequence
                protein, acc = cols[3], cols[1] if cols[1] != "-" else cols[0]
            else:  # hmmsearch: target = protein sequence, query = profile
                protein, acc = cols[0], cols[4] if cols[4] != "-" else cols[3]
            start, end = _coords(cols[19], cols[20], where)  # envelope coordinates
            seq_len = cols[5] if fmt == "hmmscan" else cols[2]  # qlen / tlen = protein length
            if not seq_len.isdigit() or end > int(seq_len):
                raise InputError(f"{where}: envelope end {end} exceeds protein length {seq_len}; wrong --format?")
            ev = _evalue(cols[12], where)  # independent E-value of this domain
        protein = _ident(protein, "protein ID", where)
        acc = _norm_acc(_ident(acc, "domain accession", where), cfg)
        if mapping:
            acc = mapping.get(acc, acc)
        run.setdefault(protein, []).append(Hit(start, end, acc, ev))
    if not run and skipped:
        found = ", ".join(sorted(skipped))
        raise InputError(f"{path}: no rows match interproscan_analyses={sorted(analyses)}; file has: {found}")
    return run


def load_mapping(path: Path, cfg: Config) -> dict[str, str]:
    """Two-column TSV: accession in run A <TAB> accession in run B. '#' lines ignored."""
    mapping: dict[str, str] = {}
    for n, line in _lines(path, cfg):
        if not line.strip() or line.startswith("#"):
            continue
        cols = line.split("\t")
        where = f"{path}:{n}"
        if len(cols) != 2:
            raise InputError(f"{where}: mapping needs exactly 2 tab-separated columns")
        old, new = (_norm_acc(_ident(c.strip(), "accession", where), cfg) for c in cols)
        if mapping.get(old, new) != new:
            raise InputError(f"{where}: {old} mapped twice to different accessions")
        mapping[old] = new
    return mapping
