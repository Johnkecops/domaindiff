"""
Module: checks
Purpose: One pure function per comparison check: two runs in, CheckResult out.
Author: Dr. Arli Aditya Parikesit
Date: 2026
Parameters:
    - see Config in model.py; definitions of domain, architecture and co-occurrence in README.md
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field
from itertools import combinations

from .model import Config, Hit, InputError, Run

PASS, WARN, FAIL, NOT_EVALUATED = "pass", "warn", "fail", "not evaluated"
Row = tuple[str, str, str]  # (protein or "" for proteome-level, change, detail)


@dataclass
class CheckResult:
    id: str
    title: str
    status: str
    reason: str = ""
    summary: dict[str, int | float | str] = field(default_factory=dict)
    rows: list[Row] = field(default_factory=list)


def _key(h: Hit) -> tuple[int, int, str]:
    return (h.start, h.end, h.acc)


def resolve_overlaps(hits: list[Hit], rule: str) -> list[Hit]:
    """Return hits sorted by (start, end, accession) after applying the overlap rule.

    keep: all hits are kept, overlapping or not.
    best_evalue: greedily keep the hit with the lowest E-value, drop any hit sharing
        at least one residue with an already-kept hit. Requires an E-value on every hit.
    """
    hits = sorted(hits, key=_key)
    if rule == "keep":
        return hits
    if any(h.evalue is None for h in hits):
        raise InputError("overlap_rule 'best_evalue' needs an E-value on every hit; use 'keep'")
    kept: list[Hit] = []
    ranked = sorted(hits, key=lambda h: (h.evalue or 0.0, *_key(h)))
    # ponytail: O(n^2) per protein; fine for real domain counts, interval tree if needed
    for h in ranked:
        if all(h.end < k.start or h.start > k.end for k in kept):
            kept.append(h)
    return sorted(kept, key=_key)


def prepare(run: Run, cfg: Config, cutoff: float | None = None) -> Run:
    """Optional E-value filter, then overlap rule. Proteins left with no hits are dropped."""
    out: Run = {}
    for p, hits in run.items():
        if cutoff is not None:
            hits = [h for h in hits if h.evalue is not None and h.evalue <= cutoff]
        if hits:
            out[p] = resolve_overlaps(hits, cfg.overlap_rule)
    return out


def architecture(hits: list[Hit]) -> tuple[str, ...]:
    return tuple(h.acc for h in hits)  # hits already sorted by prepare()


def fmt_arch(arch: tuple[str, ...]) -> str:
    return "; ".join(arch) if arch else "(none)"


# 1. Input compatibility ---------------------------------------------------------


def compatibility(a: Run, b: Run, cfg: Config) -> CheckResult:
    pa, pb = set(a), set(b)
    union, shared = pa | pb, pa & pb
    r = CheckResult("compatibility", "Input compatibility (shared protein IDs)", PASS)
    if not union:
        r.status, r.reason = FAIL, "neither run contains any domain hit"
        return r
    jac = len(shared) / len(union)
    r.summary = {
        "proteins_a": len(pa),
        "proteins_b": len(pb),
        "shared": len(shared),
        "only_a": len(pa - pb),
        "only_b": len(pb - pa),
        "jaccard": round(jac, 4),
        "min_protein_overlap": cfg.min_protein_overlap,
    }
    if jac < cfg.min_protein_overlap:
        r.status = FAIL
        r.reason = (
            f"protein ID overlap {jac:.3f} < {cfg.min_protein_overlap}; runs are not comparable "
            "under version 1 rules (different gene models need sequence matching)"
        )
    elif pa != pb:
        r.status = WARN
        r.reason = "some proteins have hits in one run only; treated as having no domains in the other"
    return r


# 2. Gain and loss ----------------------------------------------------------------


def gain_loss(a: Run, b: Run) -> CheckResult:
    rows: list[Row] = []
    gained = lost = 0
    for p in sorted(set(a) | set(b)):
        ca = Counter(h.acc for h in a.get(p, []))
        cb = Counter(h.acc for h in b.get(p, []))
        for acc, n in sorted((cb - ca).items()):
            rows.append((p, "gained", f"{acc} x{n}"))
            gained += n
        for acc, n in sorted((ca - cb).items()):
            rows.append((p, "lost", f"{acc} x{n}"))
            lost += n
    return CheckResult(
        "gain_loss",
        "Domain gain and loss per protein",
        WARN if rows else PASS,
        summary={"proteins_changed": len({r[0] for r in rows}), "domains_gained": gained, "domains_lost": lost},
        rows=rows,
    )


# 3. Boundary shifts --------------------------------------------------------------


def _match(xs: list[Hit], ys: list[Hit]) -> list[tuple[Hit, Hit]]:
    """Greedy minimum-distance pairing of copies of one accession on one protein."""
    # ponytail: greedy, not optimal assignment; O(n*m) per accession. Hungarian if repeats mispair.
    cand = sorted(
        (abs(x.start - y.start) + abs(x.end - y.end), i, j) for i, x in enumerate(xs) for j, y in enumerate(ys)
    )
    used_x: set[int] = set()
    used_y: set[int] = set()
    pairs = []
    for _, i, j in cand:
        if i not in used_x and j not in used_y:
            used_x.add(i)
            used_y.add(j)
            pairs.append((xs[i], ys[j]))
    return pairs


def boundary_shifts(a: Run, b: Run, cfg: Config) -> CheckResult:
    rows: list[Row] = []
    tol = cfg.boundary_tolerance
    for p in sorted(set(a) & set(b)):
        accs = {h.acc for h in a[p]} & {h.acc for h in b[p]}
        for acc in sorted(accs):
            xs = [h for h in a[p] if h.acc == acc]  # already sorted by coordinates
            ys = [h for h in b[p] if h.acc == acc]
            for x, y in sorted(_match(xs, ys), key=lambda xy: (_key(xy[0]), _key(xy[1]))):
                if abs(x.start - y.start) > tol or abs(x.end - y.end) > tol:
                    rows.append((p, "shifted", f"{acc} {x.start}-{x.end} -> {y.start}-{y.end}"))
    return CheckResult(
        "boundary_shifts",
        f"Boundary shifts beyond {tol} residues",
        WARN if rows else PASS,
        summary={"shifted_domains": len(rows), "proteins_affected": len({r[0] for r in rows}), "tolerance": tol},
        rows=rows,
    )


# 4. Architecture changes ---------------------------------------------------------


def _is_subsequence(short: tuple[str, ...], long: tuple[str, ...]) -> bool:
    it = iter(long)
    return all(any(x == y for y in it) for x in short)


def classify(x: tuple[str, ...], y: tuple[str, ...]) -> str:
    if Counter(x) == Counter(y):
        return "reordering"
    if len(y) > len(x) and _is_subsequence(x, y):
        return "insertion"
    if len(y) < len(x) and _is_subsequence(y, x):
        return "deletion"
    if len(x) == len(y):
        return "substitution"
    return "mixed"


def architecture_changes(a: Run, b: Run) -> CheckResult:
    rows: list[Row] = []
    for p in sorted(set(a) | set(b)):
        x, y = architecture(a.get(p, [])), architecture(b.get(p, []))
        if x != y:
            rows.append((p, classify(x, y), f"{fmt_arch(x)} -> {fmt_arch(y)}"))
    counts = Counter(r[1] for r in rows)
    summary: dict[str, int | float | str] = {"proteins_changed": len(rows)}
    summary.update({k: counts.get(k, 0) for k in ("insertion", "deletion", "substitution", "reordering", "mixed")})
    return CheckResult("architecture", "Architecture changes", WARN if rows else PASS, summary=summary, rows=rows)


# 5. Accession turnover -----------------------------------------------------------


def accession_turnover(a: Run, b: Run, mapping_size: int | None) -> CheckResult:
    acc_a = {h.acc for hs in a.values() for h in hs}
    acc_b = {h.acc for hs in b.values() for h in hs}
    label = "" if mapping_size is not None else "unresolved_"
    rows: list[Row] = [("", f"{label}only_in_a", acc) for acc in sorted(acc_a - acc_b)]
    rows += [("", f"{label}only_in_b", acc) for acc in sorted(acc_b - acc_a)]
    r = CheckResult(
        "turnover",
        "Accession turnover across the proteome",
        WARN if rows else PASS,
        summary={
            "only_in_a": len(acc_a - acc_b),
            "only_in_b": len(acc_b - acc_a),
            "mapping_entries": mapping_size if mapping_size is not None else "none supplied",
        },
        rows=rows,
    )
    if rows and mapping_size is None:
        r.reason = "no accession mapping supplied: turnover is unresolved and may include renamed or merged entries"
    return r


# 6. Co-occurrence pairs ----------------------------------------------------------


def pair_counts(run: Run) -> Counter[tuple[str, str]]:
    """Number of proteins on which each unordered pair of distinct accessions co-occurs."""
    c: Counter[tuple[str, str]] = Counter()
    for hits in run.values():
        # ponytail: k^2 in distinct accessions per protein; k is small in real proteomes
        c.update(combinations(sorted({h.acc for h in hits}), 2))
    return c


def cooccurrence(a: Run, b: Run, cfg: Config) -> CheckResult:
    ca, cb = pair_counts(a), pair_counts(b)
    rows: list[Row] = []
    strong = 0
    for pair in sorted(set(ca) | set(cb)):
        na, nb = ca[pair], cb[pair]
        name = f"{pair[0]} + {pair[1]}"
        if nb == 0:
            rows.append(("", "only_in_a", f"{name} n={na}"))
        elif na == 0:
            rows.append(("", "only_in_b", f"{name} n={nb}"))
        elif na != nb and abs(math.log2((nb + 1) / (na + 1))) >= cfg.cooccurrence_log2_threshold:
            rows.append(("", "count_changed", f"{name} {na} -> {nb}"))
            strong += 1
    return CheckResult(
        "cooccurrence",
        "Co-occurrence pair changes",
        WARN if rows else PASS,
        summary={
            "pairs_a": len(ca),
            "pairs_b": len(cb),
            "only_in_a": sum(1 for p in ca if p not in cb),
            "only_in_b": sum(1 for p in cb if p not in ca),
            "count_changed": strong,
            "log2_threshold": cfg.cooccurrence_log2_threshold,
        },
        rows=rows,
    )


# 7. Proteome summary -------------------------------------------------------------


def proteome_stats(run: Run) -> dict[str, int]:
    return {
        "proteins_with_domains": sum(1 for hs in run.values() if hs),
        "distinct_domains": len({h.acc for hs in run.values() for h in hs}),
        "distinct_architectures": len({architecture(hs) for hs in run.values() if hs}),
    }


def proteome_summary(a: Run, b: Run) -> CheckResult:
    sa, sb = proteome_stats(a), proteome_stats(b)
    rows: list[Row] = [("", k, f"{sa[k]} -> {sb[k]} ({sb[k] - sa[k]:+d})") for k in sa]
    summary: dict[str, int | float | str] = {f"{k}_a": v for k, v in sa.items()}
    summary.update({f"{k}_b": v for k, v in sb.items()})
    return CheckResult("summary", "Proteome-level summary", PASS if sa == sb else WARN, summary=summary, rows=rows)


# 8. Threshold sensitivity --------------------------------------------------------


def _n_differences(a: Run, b: Run) -> int:
    gl = gain_loss(a, b).summary
    return int(gl["domains_gained"]) + int(gl["domains_lost"]) + len(architecture_changes(a, b).rows)


def threshold_sensitivity(raw_a: Run, raw_b: Run, cfg: Config) -> CheckResult:
    r = CheckResult("threshold", f"Threshold sensitivity at E-value <= {cfg.evalue_cutoff:g}", NOT_EVALUATED)
    hits = [h for run in (raw_a, raw_b) for hs in run.values() for h in hs]
    if not hits or any(h.evalue is None for h in hits):
        r.reason = "not every hit has an E-value in both runs"
        return r
    before = _n_differences(prepare(raw_a, cfg), prepare(raw_b, cfg))
    after = _n_differences(prepare(raw_a, cfg, cfg.evalue_cutoff), prepare(raw_b, cfg, cfg.evalue_cutoff))
    r.status = WARN if after != before else PASS
    r.summary = {"differences_before": before, "differences_after": after, "change": after - before}
    r.reason = "differences = domains gained + lost + proteins with changed architecture"
    return r


# Driver --------------------------------------------------------------------------


def compare(raw_a: Run, raw_b: Run, cfg: Config, mapping_size: int | None = None) -> list[CheckResult]:
    a, b = prepare(raw_a, cfg), prepare(raw_b, cfg)
    compat = compatibility(a, b, cfg)
    later = [
        ("gain_loss", "Domain gain and loss per protein"),
        ("boundary_shifts", "Boundary shifts"),
        ("architecture", "Architecture changes"),
        ("turnover", "Accession turnover across the proteome"),
        ("cooccurrence", "Co-occurrence pair changes"),
        ("summary", "Proteome-level summary"),
        ("threshold", "Threshold sensitivity"),
    ]
    if compat.status == FAIL:
        return [compat] + [CheckResult(i, t, NOT_EVALUATED, "input compatibility failed") for i, t in later]
    return [
        compat,
        gain_loss(a, b),
        boundary_shifts(a, b, cfg),
        architecture_changes(a, b),
        accession_turnover(a, b, mapping_size),
        cooccurrence(a, b, cfg),
        proteome_summary(a, b),
        threshold_sensitivity(raw_a, raw_b, cfg),
    ]
