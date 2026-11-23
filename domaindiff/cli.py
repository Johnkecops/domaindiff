"""
Module: cli
Purpose: Command-line entry point: domaindiff RUN_A RUN_B --format FMT.
Author: Dr. Arli Aditya Parikesit
Date: 2026
Parameters:
    - run_a, run_b: annotation tables of the same proteome (same protein IDs)
    - --format: interproscan | hmmscan | hmmsearch (applies to both runs)
    - --mapping: optional TSV of run-A accession -> run-B accession (renames, merges)
    - --config: optional TOML file overriding defaults (see README)
    - --output-format: text | markdown | json | tsv
    - -o/--output: write to this path instead of stdout; refuses to overwrite without --force
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .checks import FAIL, compare
from .loaders import FORMATS, load, load_mapping
from .model import InputError, load_config
from .report import OUTPUT_FORMATS, render

EXIT_OK, EXIT_INPUT_ERROR, EXIT_INCOMPATIBLE = 0, 2, 3


def _printable(s: str) -> str:
    """Replace control characters (terminal escape sequences, newlines) in file names."""
    return "".join(c if c.isprintable() else "?" for c in s)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="domaindiff", description="Compare two protein domain annotation runs.")
    ap.add_argument("run_a", type=Path)
    ap.add_argument("run_b", type=Path)
    ap.add_argument("--format", required=True, choices=FORMATS, help="input format of both runs")
    ap.add_argument("--mapping", type=Path, help="TSV: run-A accession <TAB> run-B accession")
    ap.add_argument("--config", type=Path, help="TOML config file")
    ap.add_argument("--output-format", default="text", choices=OUTPUT_FORMATS)
    ap.add_argument("-o", "--output", type=Path, help="output file (default: stdout)")
    ap.add_argument("--force", action="store_true", help="overwrite --output if it exists")
    ap.add_argument("--version", action="version", version=f"domaindiff {__version__}")
    args = ap.parse_args(argv)

    try:
        if args.output is not None:
            inputs = [p for p in (args.run_a, args.run_b, args.mapping, args.config) if p is not None]
            if any(args.output.resolve() == p.resolve() for p in inputs):
                raise InputError("--output must not be one of the input files")
        cfg = load_config(args.config)
        mapping = load_mapping(args.mapping, cfg) if args.mapping else None
        a = load(args.run_a, args.format, cfg, mapping)
        b = load(args.run_b, args.format, cfg)
        results = compare(a, b, cfg, len(mapping) if mapping is not None else None)
    except InputError as e:
        print(f"domaindiff: error: {_printable(str(e))}", file=sys.stderr)
        return EXIT_INPUT_ERROR

    meta = {
        "version": __version__,
        "run_a": _printable(args.run_a.name),
        "run_b": _printable(args.run_b.name),
        "format": args.format,
        "mapping": _printable(args.mapping.name) if args.mapping else "none",
        "overlap_rule": cfg.overlap_rule,
        "interproscan_analyses": ",".join(cfg.interproscan_analyses),
        "coordinates": "1-based inclusive",
    }
    text = render(results, meta, args.output_format)
    if args.output is None:
        sys.stdout.write(text)
    else:
        try:
            with args.output.open("w" if args.force else "x", encoding="utf-8") as fh:
                fh.write(text)
        except FileExistsError:
            print(f"domaindiff: error: {_printable(str(args.output))} exists; use --force", file=sys.stderr)
            return EXIT_INPUT_ERROR
        except OSError as e:
            print(f"domaindiff: error: cannot write {_printable(str(args.output))}: {e.strerror}", file=sys.stderr)
            return EXIT_INPUT_ERROR
    return EXIT_INCOMPATIBLE if results[0].status == FAIL else EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
