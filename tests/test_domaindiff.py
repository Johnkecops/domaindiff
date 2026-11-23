"""
Module: test_domaindiff
Purpose: Known-answer, invariance, hostile-input and scale tests. All data is synthetic.
Author: Dr. Arli Aditya Parikesit
Date: 2026
"""

from __future__ import annotations

import json
import random
import time
from collections import Counter
from pathlib import Path

import pytest

from domaindiff.checks import (
    FAIL,
    NOT_EVALUATED,
    PASS,
    WARN,
    classify,
    compare,
    resolve_overlaps,
)
from domaindiff.cli import main
from domaindiff.loaders import load, load_mapping
from domaindiff.model import Config, Hit, InputError, Run, load_config
from domaindiff.report import md_escape, neutralize_cell, render

CFG = Config()


def ips(p: str, acc: str, s: int, e: int, ev: str = "1e-20", analysis: str = "Pfam") -> str:
    return f"{p}\tabc123\t500\t{analysis}\t{acc}\tdesc\t{s}\t{e}\t{ev}\tT\t23-11-2026"


def write(tmp: Path, name: str, lines: list[str] | str) -> Path:
    path = tmp / name
    path.write_text(lines if isinstance(lines, str) else "\n".join(lines) + "\n", encoding="utf-8")
    return path


# Synthetic pair with known differences:
#  P1: PF00002 removed            -> lost, deletion
#  P2: PF00003 end 200 -> 230     -> boundary shift (tolerance 5)
#  P3: PF00004 / PF00005 swapped  -> reordering
#  P4: PF00006 renamed PF00099    -> lost+gained without mapping, nothing with mapping
#  P5: unchanged
RUN_A = [
    ips("P1", "PF00001", 10, 100),
    ips("P1", "PF00002", 150, 250),
    ips("P2", "PF00003", 20, 200),
    ips("P3", "PF00004", 10, 100),
    ips("P3", "PF00005", 150, 250),
    ips("P4", "PF00006", 30, 130),
    ips("P5", "PF00001", 10, 100),
]
RUN_B = [
    ips("P1", "PF00001", 10, 100),
    ips("P2", "PF00003", 20, 230),
    ips("P3", "PF00005", 10, 100),
    ips("P3", "PF00004", 150, 250),
    ips("P4", "PF00099", 30, 130),
    ips("P5", "PF00001", 10, 100),
]


def run_pair(tmp: Path, a: list[str] = RUN_A, b: list[str] = RUN_B, mapping: dict[str, str] | None = None):
    ra = load(write(tmp, "a.tsv", a), "interproscan", CFG, mapping)
    rb = load(write(tmp, "b.tsv", b), "interproscan", CFG)
    return {r.id: r for r in compare(ra, rb, CFG, len(mapping) if mapping is not None else None)}


def test_known_answers(tmp_path: Path) -> None:
    r = run_pair(tmp_path)
    assert r["compatibility"].status == PASS
    assert set(r["gain_loss"].rows) == {
        ("P1", "lost", "PF00002 x1"),
        ("P4", "lost", "PF00006 x1"),
        ("P4", "gained", "PF00099 x1"),
    }
    assert r["boundary_shifts"].rows == [
        ("P2", "shifted", "PF00003 20-200 -> 20-230"),
        ("P3", "shifted", "PF00004 10-100 -> 150-250"),  # a swap moves both domains
        ("P3", "shifted", "PF00005 150-250 -> 10-100"),
    ]
    arch = {row[0]: row[1] for row in r["architecture"].rows}
    assert arch == {"P1": "deletion", "P3": "reordering", "P4": "substitution"}
    assert set(r["turnover"].rows) == {
        ("", "unresolved_only_in_a", "PF00002"),  # P1 held the only copy
        ("", "unresolved_only_in_a", "PF00006"),
        ("", "unresolved_only_in_b", "PF00099"),
    }
    assert "unresolved" in r["turnover"].reason
    # PF00001+PF00002 co-occur on P1 in A only; PF00004+PF00005 on P3 in both.
    assert r["cooccurrence"].rows == [("", "only_in_a", "PF00001 + PF00002 n=1")]
    assert r["summary"].summary["distinct_architectures_a"] == 5
    assert r["threshold"].status == PASS  # all E-values 1e-20, filter changes nothing


def test_mapping_removes_rename(tmp_path: Path) -> None:
    mp = load_mapping(write(tmp_path, "map.tsv", ["# old\tnew", "PF00006.3\tPF00099"]), CFG)
    assert mp == {"PF00006": "PF00099"}  # version suffix stripped
    r = run_pair(tmp_path, mapping=mp)
    assert not any(row[0] == "P4" for row in r["gain_loss"].rows)
    assert r["turnover"].rows == [("", "only_in_a", "PF00002")]  # real loss, not labelled unresolved


def test_identity(tmp_path: Path) -> None:
    for res in run_pair(tmp_path, RUN_A, RUN_A).values():
        assert res.status == PASS, res.id
        assert res.rows == [] or res.id == "summary"


def test_symmetry(tmp_path: Path) -> None:
    fwd, rev = run_pair(tmp_path, RUN_A, RUN_B), run_pair(tmp_path, RUN_B, RUN_A)
    swap = {"gained": "lost", "lost": "gained", "insertion": "deletion", "deletion": "insertion"}
    assert sorted((p, swap.get(c, c), d) for p, c, d in fwd["gain_loss"].rows) == sorted(rev["gain_loss"].rows)
    for (pf, cf, df), (pr, cr, dr) in zip(fwd["architecture"].rows, rev["architecture"].rows, strict=True):
        x, y = df.split(" -> ")
        assert (pf, swap.get(cf, cf), f"{y} -> {x}") == (pr, cr, dr)
    assert len(fwd["boundary_shifts"].rows) == len(rev["boundary_shifts"].rows)
    s = fwd["compatibility"].summary
    assert s["jaccard"] == rev["compatibility"].summary["jaccard"]


def test_order_invariance(tmp_path: Path) -> None:
    base = run_pair(tmp_path)
    rng = random.Random(1)
    for _ in range(5):
        a, b = RUN_A[:], RUN_B[:]
        rng.shuffle(a)
        rng.shuffle(b)
        shuffled = run_pair(tmp_path, a, b)
        assert {k: (v.status, v.summary, v.rows) for k, v in shuffled.items()} == {
            k: (v.status, v.summary, v.rows) for k, v in base.items()
        }


def test_independent_cross_check() -> None:
    """Recompute gain/loss and pair counts with plain set logic on a random run."""
    rng = random.Random(7)
    accs = [f"PF{i:05d}" for i in range(12)]

    def rand_run() -> Run:
        run: Run = {}
        for p in range(60):
            k = rng.randint(1, 4)
            run[f"P{p}"] = [Hit(s * 100 + 1, s * 100 + 80, rng.choice(accs), 1e-10) for s in range(k)]
        return run

    a, b = rand_run(), rand_run()
    r = {x.id: x for x in compare(a, b, CFG)}
    gained = sum(max(0, [h.acc for h in b[p]].count(acc) - [h.acc for h in a[p]].count(acc)) for p in a for acc in accs)
    assert r["gain_loss"].summary["domains_gained"] == gained

    def pairs(run: Run) -> Counter[tuple[str, str]]:
        c: Counter[tuple[str, str]] = Counter()
        for hits in run.values():
            ds = sorted({h.acc for h in hits})
            c.update((x, y) for i, x in enumerate(ds) for y in ds[i + 1 :])
        return c

    pa, pb = pairs(a), pairs(b)
    assert r["cooccurrence"].summary["only_in_a"] == len(set(pa) - set(pb))
    assert r["cooccurrence"].summary["pairs_b"] == len(pb)


def test_incompatible_runs_stop(tmp_path: Path) -> None:
    r = run_pair(tmp_path, RUN_A, [ips("X1", "PF00001", 1, 50)])
    assert r["compatibility"].status == FAIL
    assert all(v.status == NOT_EVALUATED for k, v in r.items() if k != "compatibility")


def test_classify() -> None:
    assert classify(("A",), ("A", "B")) == "insertion"
    assert classify(("A", "B", "C"), ("A", "C")) == "deletion"
    assert classify(("A", "B"), ("A", "C")) == "substitution"
    assert classify(("A", "B"), ("B", "A")) == "reordering"
    assert classify(("A", "B"), ("C", "D", "E")) == "mixed"
    assert classify((), ("A",)) == "insertion"


def test_overlap_rule() -> None:
    hits = [
        Hit(10, 100, "PF1", 1e-5),
        Hit(50, 150, "PF2", 1e-30),
        Hit(160, 200, "PF3", 1e-3),
        Hit(200, 220, "PF4", 1.0),
    ]
    assert [h.acc for h in resolve_overlaps(hits, "keep")] == ["PF1", "PF2", "PF3", "PF4"]
    # PF4 shares residue 200 with PF3: touching counts as overlap
    assert [h.acc for h in resolve_overlaps(hits, "best_evalue")] == ["PF2", "PF3"]
    with pytest.raises(InputError, match="best_evalue"):
        resolve_overlaps([Hit(1, 10, "PF1", None)], "best_evalue")


def test_repeats_matched_by_position(tmp_path: Path) -> None:
    a = [ips("P", "PF1", 1, 50), ips("P", "PF1", 100, 150), ips("P", "PF1", 200, 250)]
    b = [ips("P", "PF1", 100, 150), ips("P", "PF1", 200, 250)]
    r = run_pair(tmp_path, a, b)
    assert r["boundary_shifts"].rows == []  # lost copy is not misreported as shifts
    assert r["gain_loss"].rows == [("P", "lost", "PF1 x1")]


def test_threshold_sensitivity(tmp_path: Path) -> None:
    a = [ips("P1", "PF1", 1, 50), ips("P1", "PF2", 60, 120, "0.01")]
    b = [ips("P1", "PF1", 1, 50)]
    r = run_pair(tmp_path, a, b)
    assert r["threshold"].summary == {"differences_before": 2, "differences_after": 0, "change": -2}
    assert r["threshold"].status == WARN
    r = run_pair(tmp_path, [ips("P1", "PF1", 1, 50, "-")], b)
    assert r["threshold"].status == NOT_EVALUATED


def test_hmmer_formats(tmp_path: Path) -> None:
    head = "# target name  accession  ...\n"
    scan = "PK  PF00069.29  264  sp|P1|X  -  400  1e-50  170.0  0.1  1  1  1e-52  1e-49  169.0  0.1  1  264  20  280  18  282  0.95  Protein kinase domain\n"  # noqa: E501
    search = "sp|P1|X  -  400  PK  PF00069.29  264  1e-50  170.0  0.1  1  1  1e-52  1e-49  169.0  0.1  1  264  20  280  18  282  0.95  some protein\n"  # noqa: E501
    for fmt, body in (("hmmscan", scan), ("hmmsearch", search)):
        run = load(write(tmp_path, f"{fmt}.txt", head + body), fmt, CFG)
        assert run == {"sp|P1|X": [Hit(18, 282, "PF00069", 1e-49)]}


def test_interproscan_analysis_filter_and_columns(tmp_path: Path) -> None:
    lines = [ips("P1", "PF1", 1, 50), ips("P1", "SM00220", 1, 50, analysis="SMART")]
    lines.append(ips("P1", "G3DSA:3.30.200.20", 1, 50, analysis="Gene3D") + "\tIPR000719\tProtein kinase\t-\t-")
    assert {h.acc for h in load(write(tmp_path, "x.tsv", lines), "interproscan", CFG)["P1"]} == {"PF1"}
    allcfg = Config(interproscan_analyses=["all"])
    accs = {h.acc for h in load(write(tmp_path, "x.tsv", lines), "interproscan", allcfg)["P1"]}
    assert "G3DSA:3.30.200.20" in accs  # punctuated accession not version-stripped


@pytest.mark.parametrize(
    "content",
    [
        "",
        "\n\n",
        "P1\tonly\tthree\n",
        ips("P1", "PF1", "x", 10),  # type: ignore[arg-type]
        ips("P1", "PF1", 50, 10),
        ips("P1", "PF1", -5, 10),
        ips("P1", "PF1", 0, 10),
        ips("P" * 5000, "PF1", 1, 10),
        ips("P1\x1b[31m", "PF1", 1, 10),
        ips("P1", "PF1", 1, 10, "nan"),
        ips("P1", "PF1", 1, 10, "abc"),
        ips("P1", "PF1", 1, 10) + "\n" + "PK PF1 264 P1 - 400 " + "1 " * 16,  # mixed formats
    ],
)
def test_hostile_inputs_clean(tmp_path: Path, content: str) -> None:
    path = write(tmp_path, "bad.tsv", content)
    try:
        run = load(path, "interproscan", CFG)
    except InputError:
        return
    assert run == {}  # only empty inputs may load without error


def test_non_utf8_gzip_and_nonregular(tmp_path: Path) -> None:
    bad = tmp_path / "bad.tsv"
    bad.write_bytes(b"P1\xff\xfe\n")
    with pytest.raises(InputError, match="UTF-8"):
        load(bad, "interproscan", CFG)
    bad.write_bytes(b"\x1f\x8b\x08\x00")
    with pytest.raises(InputError, match="gzip"):
        load(bad, "interproscan", CFG)
    with pytest.raises(InputError, match="regular file"):
        load(tmp_path, "interproscan", CFG)
    with pytest.raises(InputError, match="max_rows"):
        load(write(tmp_path, "m.tsv", [ips("P1", "PF1", 1, 10)] * 3), "interproscan", Config(max_rows=2))


def test_config(tmp_path: Path) -> None:
    assert load_config(write(tmp_path, "c.toml", "boundary_tolerance = 10\nevalue_cutoff = 1")).boundary_tolerance == 10
    for bad in ("nope = 1", "boundary_tolerance = 'x'", "overlap_rule = 'random'", "strip_accession_version = 1"):
        with pytest.raises(InputError):
            load_config(write(tmp_path, "c.toml", bad))


def test_output_escaping() -> None:
    assert neutralize_cell("=HYPERLINK(1)") == "'=HYPERLINK(1)"
    assert neutralize_cell("+1") == "'+1" and neutralize_cell("@x") == "'@x" and neutralize_cell("-1") == "'-1"
    assert neutralize_cell("P1") == "P1"
    assert md_escape("<script>|x") == "\\<script\\>\\|x"


def test_cli_all_formats_and_overwrite(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    a, b = write(tmp_path, "a.tsv", RUN_A), write(tmp_path, "b.tsv", RUN_B + [ips("=cmd", "PF1", 1, 10)])
    for fmt in ("text", "markdown", "json", "tsv"):
        assert main([str(a), str(b), "--format", "interproscan", "--output-format", fmt]) == 0
        out = capsys.readouterr().out
        assert "PF00099" in out
        if fmt == "tsv":
            assert "'=cmd" in out and "\t=cmd" not in out
        if fmt == "json":
            assert json.loads(out)["checks"][0]["status"] == "warn"
    out_path = tmp_path / "r.txt"
    assert main([str(a), str(b), "--format", "interproscan", "-o", str(out_path)]) == 0
    assert main([str(a), str(b), "--format", "interproscan", "-o", str(out_path)]) == 2
    assert main([str(a), str(b), "--format", "interproscan", "-o", str(out_path), "--force"]) == 0
    assert main([str(a), str(b), "--format", "interproscan", "-o", str(a), "--force"]) == 2
    assert a.read_text(encoding="utf-8").startswith("P1")  # input never overwritten
    assert main([str(a), str(tmp_path / "missing"), "--format", "interproscan"]) == 2
    assert main([str(a), str(write(tmp_path, "x.tsv", [ips("Z", "PF1", 1, 5)])), "--format", "interproscan"]) == 3


def test_render_text_has_every_check(tmp_path: Path) -> None:
    results = list(run_pair(tmp_path).values())
    text = render(results, {"version": "t"}, "text")
    assert all(r.title in text for r in results)


def test_scale() -> None:
    """Synthetic proteome-sized run (20,000 proteins, ~3 domains each). Prints timing; no hard target."""
    rng = random.Random(3)
    accs = [f"PF{i:05d}" for i in range(3000)]

    def rand_run() -> Run:
        return {
            f"P{p}": [Hit(s * 150 + 1, s * 150 + 120, rng.choice(accs), 1e-10) for s in range(rng.randint(1, 5))]
            for p in range(20_000)
        }

    a, b = rand_run(), rand_run()
    t = time.perf_counter()
    results = compare(a, b, CFG)
    elapsed = time.perf_counter() - t
    print(f"scale: 20,000 proteins compared in {elapsed:.2f}s")
    assert len(results) == 8


def test_score_columns_not_used_as_evalue(tmp_path: Path) -> None:
    cfg = Config(interproscan_analyses=["ProSiteProfiles"])
    run = load(write(tmp_path, "p.tsv", [ips("P1", "PS50011", 1, 50, "45.2", "ProSiteProfiles")]), "interproscan", cfg)
    assert run["P1"][0].evalue is None


def test_analysis_filter_matching_nothing_is_an_error(tmp_path: Path) -> None:
    with pytest.raises(InputError, match="SMART"):
        load(write(tmp_path, "s.tsv", [ips("P1", "SM00220", 1, 50, analysis="SMART")]), "interproscan", CFG)


def test_wrong_hmmer_flavour_rejected(tmp_path: Path) -> None:
    search = "sp|P1|X  -  400  PK  PF00069.29  264  1e-50  170.0  0.1  1  1  1e-52  1e-49  169.0  0.1  1  264  20  280  18  282  0.95  x\n"  # noqa: E501
    with pytest.raises(InputError, match="wrong --format"):
        load(write(tmp_path, "s.txt", search), "hmmscan", CFG)


def test_tsv_reports_every_check_status(tmp_path: Path) -> None:
    results = list(run_pair(tmp_path, RUN_A, [ips("X1", "PF00001", 1, 50)]).values())
    tsv = render(results, {}, "tsv")
    assert tsv.count("check_status") == 8 and "not evaluated" in tsv


def test_boundary_tolerance_is_exclusive(tmp_path: Path) -> None:
    a = [ips("P1", "PF1", 10, 100), ips("P1", "PF2", 200, 300)]
    b = [ips("P1", "PF1", 15, 100), ips("P1", "PF2", 206, 300)]  # 5 = tolerance, 6 > tolerance
    assert [r[2] for r in run_pair(tmp_path, a, b)["boundary_shifts"].rows] == ["PF2 200-300 -> 206-300"]


def test_jaccard_equal_to_threshold_passes() -> None:
    a: Run = {"P1": [Hit(1, 10, "PF1")], "P2": [Hit(1, 10, "PF1")]}
    b: Run = {"P1": [Hit(1, 10, "PF1")]}
    r = compare(a, b, Config(min_protein_overlap=0.5))
    assert r[0].status == WARN and r[0].summary["jaccard"] == 0.5


def test_cooccurrence_only_in_b_and_count_change() -> None:
    a: Run = {f"P{i}": [Hit(1, 10, "PF1"), Hit(20, 30, "PF2")] for i in range(2)}
    b: Run = {f"P{i}": [Hit(1, 10, "PF1"), Hit(20, 30, "PF2"), Hit(40, 50, "PF3")] for i in range(5)}
    a.update({f"P{i}": [Hit(1, 10, "PF1")] for i in range(2, 5)})
    co = {x.id: x for x in compare(a, b, CFG)}["cooccurrence"]
    rows = co.rows
    assert co.summary["count_changed"] == 1 and co.summary["only_in_b"] == 2
    assert ("", "count_changed", "PF1 + PF2 2 -> 5") in rows  # log2(6/3) = 1 >= 1
    assert ("", "only_in_b", "PF1 + PF3 n=5") in rows


def test_distinct_architectures_counted() -> None:
    a: Run = {"P1": [Hit(1, 10, "PF1")], "P2": [Hit(1, 10, "PF1")], "P3": [Hit(1, 10, "PF2")]}
    r = {x.id: x for x in compare(a, a, CFG)}["summary"]
    assert r.summary["distinct_architectures_a"] == 2 and r.summary["proteins_with_domains_a"] == 3
