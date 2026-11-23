# domaindiff

Compare two protein domain annotation runs of the same proteome and report what changed: domains gained or lost, boundary shifts, architecture changes, accession turnover and co-occurrence pairs that appear or disappear.

The question it answers: how much of a domain-level result depends on the database release or on search parameters, before any evolutionary conclusion is drawn from it. It reports differences. It does not decide which run is correct.

Version 0.1.0. Python 3.11+, standard library only, no runtime dependencies.

## Install

```bash
git clone <this repository>
cd domaindiff
python3 -m pip install .
```

## Usage

```bash
domaindiff RUN_A RUN_B --format {interproscan,hmmscan,hmmsearch} \
    [--mapping MAP.tsv] [--config domaindiff.toml] \
    [--output-format {text,markdown,json,tsv}] [-o OUT [--force]]
```

Example with the synthetic files in `examples/`:

```bash
domaindiff examples/run_a.tsv examples/run_b.tsv --format interproscan --mapping examples/mapping.tsv
```

Exit codes: `0` comparison completed (with or without differences), `2` input or output error, `3` runs not comparable (protein ID overlap below `min_protein_overlap`).

`-o` refuses to overwrite an existing file unless `--force` is given, and never writes over an input file.

## Supported inputs

Both runs must use the same format.

| `--format` | Tool output | Protein | Domain accession | Coordinates | E-value |
|---|---|---|---|---|---|
| `interproscan` | InterProScan TSV (checked against 5.78-109.0 and 6.0.3.0 docs). 11, 13, 14 or 15 tab-separated columns, no header. | col 1 | col 5 (signature accession) | cols 7-8 | col 9 (`-` = missing) |
| `hmmscan` | HMMER 3.4 `hmmscan --domtblout` | query name (field 4) | target accession (field 2), or target name if `-` | envelope, fields 20-21 | i-Evalue, field 13 |
| `hmmsearch` | HMMER 3.4 `hmmsearch --domtblout` | target name (field 1) | query accession (field 5), or query name if `-` | envelope, fields 20-21 | i-Evalue, field 13 |

Coordinates are taken as 1-based and inclusive in every format. InterProScan 6 documents this explicitly; the InterProScan 5 docs and the HMMER guide do not state it outright, but their numbering (1..L) implies it.

For InterProScan input only the analyses listed in `interproscan_analyses` are read (default: Pfam). Mixing member databases puts redundant, overlapping signatures (Pfam, SMART, Gene3D, PANTHER) into one architecture, which is rarely what you want. Set it to `["all"]` to read everything.

Column 9 of InterProScan output is an E-value for some member databases and a score (higher is better) for others, such as ProSiteProfiles and HAMAP. It is read as an E-value only for analyses listed in `interproscan_evalue_analyses` (default: Pfam). For all other analyses the E-value is treated as missing, so check 8 reports `not evaluated` and `best_evalue` refuses to run.

If no row matches `interproscan_analyses`, the tool stops with an error that lists the analyses found in the file. Analysis names are case-sensitive.

For HMMER input, a hit whose envelope ends past the protein length (qlen for hmmscan, tlen for hmmsearch) is rejected, because that usually means the wrong `--format` was given. HMMER domain tables are read unfiltered. Run HMMER with `--cut_ga` so hits are filtered at Pfam gathering thresholds, which is how InterProScan filters Pfam.

Gzip input is rejected; decompress first.

## Definitions

These are the definitions the code implements. They are scientific choices and affect results.

- **Domain hit**: one row of input after filtering: protein ID, accession, start, end, optional E-value.
- **Accession**: the signature accession. A trailing version on HMM accessions (`PF00069.29` becomes `PF00069`) is removed by default, so a model rebuild between releases is not counted as turnover. Accessions containing other punctuation, such as `G3DSA:3.30.200.20`, are left unchanged.
- **Overlap rule**: applied to both runs before every check.
  - `keep` (default): every hit is kept, overlapping or not. The tool reports what the annotation program reported.
  - `best_evalue`: hits are taken in order of increasing E-value; a hit is dropped if it shares at least one residue with a hit already kept. Every hit must have an E-value, or the tool stops with an error. This resembles the clan-based "best match wins" overlap resolution of pfam_scan.pl (Finn et al. 2006, Nucleic Acids Res 34:D247, doi:10.1093/nar/gkj149) but is not clan-aware: it resolves overlaps between any two families.
- **Architecture**: the accessions of a protein's hits ordered N- to C-terminal by start, then end, then accession. Repeats are kept (`A; A; B`).
- **Co-occurrence pair**: an unordered pair of two *distinct* accessions found on the same protein. Each protein counts once per pair, no matter how many copies. A domain with itself is not a pair.
- **Protein missing from a run**: InterProScan and HMMER only report proteins with hits, so a protein present in one run only is treated as having no domains in the other.

## Checks

Each check reports `pass` (no differences), `warn` (differences found), or `not evaluated` with a reason. `warn` means the runs differ. It does not mean either run is wrong. Only the compatibility check can `fail`. TSV output has one `check_status` row per check, followed by that check's detail rows, so a `not evaluated` check is never silently missing.

1. **Input compatibility.** Jaccard overlap of protein IDs with at least one hit. Below `min_protein_overlap` the tool stops and marks every other check `not evaluated`.
2. **Gain and loss.** Per protein, the multiset difference of accessions, so losing one copy of a repeat counts.
3. **Boundary shifts.** Copies of the same accession on the same protein are paired greedily by smallest total coordinate distance; a pair is reported if start or end moved by more than `boundary_tolerance`. Unpaired copies show up as gain or loss, not as shifts. When two domains swap places, both are also reported as shifted.
4. **Architecture changes.** Classified as `insertion` (old architecture is a subsequence of the new, longer one), `deletion` (the reverse), `reordering` (same multiset, different order), `substitution` (same length, different content), or `mixed` (anything else).
5. **Accession turnover.** Accessions found anywhere in one run and nowhere in the other. Without `--mapping`, rows are labelled `unresolved_only_in_a/b` because they may include renamed or merged entries.
6. **Co-occurrence pairs.** Pairs found in one run only, and pairs whose protein count changes by at least `cooccurrence_log2_threshold` on a log2 scale with a pseudocount of 1: |log2((nB+1)/(nA+1))|.
7. **Proteome summary.** Proteins with any domain, distinct accessions, distinct architectures, per run and the difference.
8. **Threshold sensitivity.** If every hit in both runs has an E-value, both runs are filtered at `evalue_cutoff` and the number of differences (domains gained + lost + proteins with changed architecture) is recomputed. A drop means part of the difference sits near the threshold. This assumes the score column is an E-value, which holds for Pfam in InterProScan but not for every member database.

## Accession mapping

`--mapping` takes a two-column TSV, `run A accession<TAB>run B accession`, with `#` comment lines allowed. It is applied to run A after version stripping, before any check, so a renamed or merged family is not counted as gain, loss or turnover. Only one hop is applied (no chaining). You supply the file; the tool downloads nothing. For Pfam, the `Pfam-A.dead` file on the EBI FTP site has a forward-to (`FW`) field for many retired families that can be turned into this format. We have not tested that conversion.

## Configuration

Optional TOML file (`--config`). Unknown keys are an error. Defaults:

| Key | Default | Why |
|---|---|---|
| `boundary_tolerance` | `5` | Residues. A rebuilt HMM often moves envelope edges by a few residues; smaller moves rarely matter for architecture. A heuristic, not derived from measured data. |
| `evalue_cutoff` | `1e-5` | Only used by check 8. A widespread convention, not a recommendation. Pfam itself uses curated gathering thresholds. |
| `cooccurrence_log2_threshold` | `1.0` | Two-fold change in protein count (pseudocount 1). |
| `min_protein_overlap` | `0.5` | Jaccard below this suggests different gene models or proteomes, which version 1 cannot compare. A heuristic. |
| `overlap_rule` | `"keep"` | Report what the annotation tool reported; see Definitions. |
| `interproscan_analyses` | `["Pfam"]` | See Supported inputs. |
| `interproscan_evalue_analyses` | `["Pfam"]` | Analyses whose column 9 is an E-value. Pfam is the only one checked so far; add others after checking their documentation. |
| `strip_accession_version` | `true` | See Definitions. |
| `max_file_mb` | `2048` | Input size cap. |
| `max_rows` | `10000000` | Input line cap. |

## Limits

- Version 1 compares runs on the **same protein IDs**. Runs built on different gene models need sequence-based protein matching and are out of scope.
- No cross-species comparison.
- The overlap rule is not clan-aware. Clan-based resolution (as in pfam_scan.pl) would need the user to supply `Pfam-A.clans.tsv`.
- Boundary-shift pairing of repeated domains is greedy, not an optimal assignment. With many close repeats a copy can be paired with the wrong partner.
- Mixing input formats (one InterProScan run, one HMMER run) is not supported.
- No network access, no telemetry, no execution of external programs. It never runs InterProScan or HMMER.

## Open questions for domain review

These were raised in code review. They are scientific choices, so they are left for the maintainer to decide.

- **Compatibility conflates loss with mismatch.** Jaccard overlap uses only proteins with hits. If many proteins lose every domain in run B, the runs are called "not comparable" when this is really a large loss. Possible fixes: overlap = shared / smaller run, or an optional list of all proteome IDs.
- **Many-to-one merges in the mapping.** If PF1 and PF2 both map to PF3 and run B has one PF3 spanning both, the report shows a lost PF3, a boundary shift and a deletion. Possible fixes: merge adjacent identical hits after mapping, or label such rows "merge".
- **Classification choices.** A swap shows up as two boundary shifts. `A; B -> B; C` is labelled substitution. Repeats of one domain do not count as a co-occurrence pair.
- **HMMER filtering.** domtblout hits are used unfiltered (see Supported inputs).

## Prior art

A search on 23 November 2026 (web, OpenAlex, arXiv; a GitHub code search was rate-limited after one query; Bioconda and PyPI were not searched directly) found no tool that diffs domain annotations or architectures between two runs of the same proteome. Related tools compare proteins, not runs:

- DAMA (Bernardes et al. 2015, doi:10.1093/bioinformatics/btv582) resolves one protein's architecture using co-occurrence.
- CODD (Terrapon et al. 2009, doi:10.1093/bioinformatics/btp560) certifies weak domains by co-occurrence.
- DoMosaics (Moore et al. 2013, doi:10.1093/bioinformatics/btt640) visualises and compares domain arrangements.
- RADS/RAMPAGE (Terrapon et al. 2013, doi:10.1093/bioinformatics/btt379) searches by arrangement alignment.
- DOGMA (Dohmen et al. 2016, doi:10.1093/bioinformatics/btw231) scores proteome completeness from arrangements.
- FAS (Dosch et al. 2023, doi:10.1093/bioinformatics/btad226) scores feature-architecture similarity of protein pairs.

Outcome: build, keeping to the narrow niche of run-to-run comparison.

## Development

```bash
python3 -m venv .venv && .venv/bin/pip install --require-hashes -r requirements-dev.txt
.venv/bin/python -m pytest
.venv/bin/ruff check . && .venv/bin/ruff format --check . && .venv/bin/mypy
```

All test data is synthetic. See `SECURITY.md` for the audit record.

## Licence

MIT (see `LICENSE`).

AI Assistance Disclaimer: This codebase was developed with the assistance of Claude Code. While the AI provided code generation, debugging, and structural support, the human developer maintains full responsibility for reviewing, testing, and maintaining all content and functionality
