"""
Module: report
Purpose: Render check results as text, Markdown, JSON or TSV, escaping untrusted input.
Author: Dr. Arli Aditya Parikesit
Date: 2026
"""

from __future__ import annotations

import csv
import io
import json
import re
from dataclasses import asdict
from typing import Any

from .checks import CheckResult

OUTPUT_FORMATS = ("text", "markdown", "json", "tsv")
_MD_SPECIAL = re.compile(r"([\\`*_{}\[\]()#+\-.!|<>~&])")


def neutralize_cell(value: str) -> str:
    """Stop spreadsheets treating a cell as a formula (OWASP CSV injection guidance)."""
    return "'" + value if value[:1] in ("=", "+", "-", "@", "\t", "\r") else value


def md_escape(value: str) -> str:
    return _MD_SPECIAL.sub(r"\\\1", value)


def render(results: list[CheckResult], meta: dict[str, Any], fmt: str) -> str:
    if fmt == "json":
        return json.dumps({"meta": meta, "checks": [asdict(r) for r in results]}, indent=2) + "\n"
    if fmt == "tsv":
        buf = io.StringIO()
        w = csv.writer(buf, delimiter="\t", lineterminator="\n")
        w.writerow(["check", "status", "protein", "change", "detail"])
        for r in results:
            info = r.reason or "; ".join(f"{k}={v}" for k, v in r.summary.items())
            w.writerow([neutralize_cell(c) for c in (r.id, r.status, "", "check_status", info)])
            for row in r.rows:
                w.writerow([neutralize_cell(c) for c in (r.id, r.status, *row)])
        return buf.getvalue()
    if fmt == "markdown":
        out = ["# domaindiff report", ""]
        out += [f"- {md_escape(k)}: {md_escape(str(v))}" for k, v in meta.items()]
        for r in results:
            out += ["", f"## {md_escape(r.title)}: {r.status}", ""]
            if r.reason:
                out += [md_escape(r.reason), ""]
            out += [f"- {md_escape(k)}: {md_escape(str(v))}" for k, v in r.summary.items()]
            if r.rows:
                out += ["", "| protein | change | detail |", "|---|---|---|"]
                out += ["| " + " | ".join(md_escape(c) for c in row) + " |" for row in r.rows]
        return "\n".join(out) + "\n"
    out = ["domaindiff report"] + [f"  {k}: {v}" for k, v in meta.items()]
    for r in results:
        out += ["", f"[{r.status.upper()}] {r.title}"]
        if r.reason:
            out.append(f"  note: {r.reason}")
        out += [f"  {k}: {v}" for k, v in r.summary.items()]
        out += ["    " + "\t".join(c or "*" for c in row) for row in r.rows]
    return "\n".join(out) + "\n"
