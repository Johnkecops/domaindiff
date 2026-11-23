#!/usr/bin/env python3
"""
Module: pfam_dead_to_mapping
Purpose: Convert a Pfam-A.dead file into a domaindiff --mapping TSV (dead accession -> forward accession).
Author: Dr. Arli Aditya Parikesit
Date: 2026
Parameters:
    - dead: path to Pfam-A.dead (plain text or .gz), Stockholm-style records ending with '//'
    - out: output TSV path (written with exclusive create)
Notes:
    Only records with exactly one forward (FW) accession are written. Families killed without a
    forward target, or split into several targets, are counted and skipped, because a one-to-one
    rename cannot represent them.
References:
    - Pfam release files, EBI FTP: https://ftp.ebi.ac.uk/pub/databases/Pfam/
"""

import gzip
import re
import sys
from pathlib import Path

ACC = re.compile(r"PF\d{5}")


def records(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    rec: dict[str, str] = {}
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("//"):
                yield rec
                rec = {}
            elif line.startswith("#=GF "):
                tag, _, value = line[5:].partition(" ")
                rec[tag] = (rec.get(tag, "") + " " + value.strip()).strip()


def main(dead: str, out: str) -> None:
    written = no_target = multi = 0
    with open(out, "x", encoding="utf-8") as fh:
        fh.write("# Pfam dead family accession\tforward accession (from Pfam-A.dead FW)\n")
        for rec in records(Path(dead)):
            acc = ACC.search(rec.get("AC", ""))
            targets = ACC.findall(rec.get("FW", ""))
            if not acc:
                continue
            if not targets:
                no_target += 1
            elif len(set(targets)) > 1:
                multi += 1
            else:
                fh.write(f"{acc.group()}\t{targets[0]}\n")
                written += 1
    print(f"written={written} no_forward={no_target} multiple_forward={multi}", file=sys.stderr)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
