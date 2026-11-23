#!/usr/bin/env python3
"""
Module: run_case_study
Purpose: Compare M. tuberculosis H37Rv Pfam 35.0 vs 39.0 hmmsearch runs with domaindiff; write tables and Figure 2.
Author: Dr. Arli Aditya Parikesit
Date: 2026
Parameters:
    - data/raw/mtb_pfam35.domtbl, data/raw/mtb_pfam39.domtbl: hmmsearch --cut_ga --domtblout output
    - data/raw/pfam39_dead_mapping.tsv: from scripts/pfam_dead_to_mapping.py
Outputs:
    - results/tables/case_checks.json, results/tables/case_<variant>.tsv, results/figures/fig2_case_study.png
References:
    - Eddy SR (2011) PLoS Comput Biol 7:e1002195; Mistry J et al. (2021) NAR 49:D412
"""

import json
import os
import subprocess  # noqa: S404  (runs domaindiff CLI on local files only)
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW, TAB, FIG = ROOT / "data/raw", ROOT / "results/tables", ROOT / "results/figures"
A, B, MAP = RAW / "mtb_pfam35.domtbl", RAW / "mtb_pfam39.domtbl", RAW / "pfam39_dead_mapping.tsv"


def run(extra: list[str], fmt: str) -> str:
    cmd = [sys.executable, "-m", "domaindiff", str(A), str(B), "--format", "hmmsearch", "--output-format", fmt]
    res = subprocess.run(cmd + extra, cwd=ROOT, capture_output=True, text=True, check=False)  # noqa: S603
    if res.returncode not in (0, 3):
        sys.exit(res.stderr)
    return res.stdout


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    cfg = TAB / "best_evalue.toml"
    cfg.write_text('overlap_rule = "best_evalue"\n', encoding="utf-8")
    variants = {
        "unmapped": [],
        "mapped": ["--mapping", str(MAP)],
        "mapped_best": ["--mapping", str(MAP), "--config", str(cfg)],
    }
    checks = {}
    for name, extra in variants.items():
        checks[name] = json.loads(run(extra, "json"))["checks"]
        (TAB / f"case_{name}.tsv").write_text(run(extra, "tsv"), encoding="utf-8")
    (TAB / "case_checks.json").write_text(json.dumps(checks, indent=1), encoding="utf-8")

    os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "mpl"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def get(name: str, cid: str) -> dict:
        return next(c for c in checks[name] if c["id"] == cid)["summary"]

    classes = ["insertion", "deletion", "substitution", "reordering", "mixed"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.0), dpi=300)
    x = range(len(classes))
    for off, name, col in ((-0.2, "unmapped", "#9db4cf"), (0.2, "mapped", "#2f5f8f")):
        vals = [get(name, "architecture")[k] for k in classes]
        ax1.bar(
            [i + off for i in x],
            vals,
            width=0.4,
            color=col,
            label="no mapping" if name == "unmapped" else "with mapping",
        )
        for i, v in zip(x, vals, strict=True):
            ax1.text(i + off, v, str(v), ha="center", va="bottom", fontsize=5.5)
    ax1.set_xticks(list(x), classes, rotation=30, fontsize=7)
    ax1.set_ylabel("Proteins", fontsize=8)
    ax1.set_title("Architecture change class", fontsize=8.5)
    ax1.legend(fontsize=6.5, frameon=False)
    ax1.tick_params(labelsize=7)
    labels = ["only in Pfam 35.0", "only in Pfam 39.0"]
    for off, name, col in ((-0.2, "unmapped", "#9db4cf"), (0.2, "mapped", "#2f5f8f")):
        t = get(name, "turnover")
        vals = [t["only_in_a"], t["only_in_b"]]
        ax2.bar([i + off for i in range(2)], vals, width=0.4, color=col)
        for i, v in enumerate(vals):
            ax2.text(i + off, v, str(v), ha="center", va="bottom", fontsize=6)
    ax2.set_xticks([0, 1], labels, fontsize=7)
    ax2.set_ylabel("Accessions", fontsize=8)
    ax2.set_title("Accession turnover", fontsize=8.5)
    ax2.tick_params(labelsize=7)
    for ax in (ax1, ax2):
        ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_case_study.png")
    print("ok")


if __name__ == "__main__":
    main()
