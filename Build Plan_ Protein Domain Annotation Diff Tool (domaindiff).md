# **Build Plan: Protein Domain Annotation Diff Tool**

Prepared 20 November 2026\. Working name: **domaindiff**. Target machine: MacBook Pro M3. Built with Claude Code.

## **1\. Summary**

Build a Python command-line tool that compares two protein domain annotation runs of the same proteome and reports what changed: domains gained or lost, altered domain architectures, and co-occurrence pairs that appear or disappear. The purpose is to let a researcher see how much of a domain-level result depends on the database release or the gene models, before drawing evolutionary conclusions from it.  
The tool reads two tables and writes one report. It runs no annotation software, makes no network calls, and executes nothing from its input.  
Release is blocked by two gates: a comprehensive code review (section 7\) and a security audit (section 8). Both must pass before the repository goes public.  
Effort estimate: two to three weekends. This is my estimate, not a measured figure. It holds only if version 1 is limited to runs that share protein identifiers (section 2). Matching proteins across different gene models would raise the effort noticeably.

## **2\. Scope**

**In scope**

> * Input: two domain annotation tables for the same species, in a common tabular format (for example InterProScan TSV or HMMER domain table output). Each row: protein ID, domain accession, start, end, score or E-value.  
> * Version 1 comparison case: **same protein set, different domain database release or different search parameters.** Protein IDs match exactly between runs.  
> * Output: a result per check (pass, warn, not evaluated), a table of changed proteins, and a summary. Formats: text, Markdown, JSON, TSV.

**Out of scope for version 1**

> * Running InterProScan, HMMER or any annotation tool.  
> * Comparing runs built on different gene models, where protein IDs differ. This needs sequence-based matching and is deferred to version 2\.  
> * Cross-species comparison.  
> * Deciding which annotation is correct. The tool reports differences; the scientist interprets them.  
> * Downloading database releases or accession mapping files. Users supply any mapping file themselves, which keeps the tool offline.

## **3\. Checks in version 1**

> 1. **Input compatibility.** Do the two runs cover the same protein IDs? Report the overlap. If overlap is low, stop and say the runs are not comparable under version 1 rules.  
> 2. **Domain gain and loss per protein.** Domains present in one run and absent in the other, per protein.  
> 3. **Boundary shifts.** Same domain on the same protein with start or end moved beyond a configurable tolerance.  
> 4. **Architecture changes.** The ordered list of domains on a protein differs between runs. Classified as insertion, deletion, substitution or reordering.  
> 5. **Accession turnover.** Domain accessions that exist in only one run across the whole proteome, which usually indicates entries added, retired or merged between releases. If the user supplies a mapping file of renamed or merged accessions, apply it first so renames are not counted as biological change.  
> 6. **Co-occurrence pair changes.** Domain pairs that co-occur on at least one protein in one run but in none in the other, and pairs whose co-occurrence count changes strongly.  
> 7. **Proteome-level summary.** Total proteins with any domain, total distinct domains, total distinct architectures, in each run and the difference.  
> 8. **Threshold sensitivity (optional).** If scores or E-values are present, report how many differences disappear when both runs are filtered at the same cut-off. This separates real annotation change from threshold effects.

## **4\. Design**

> * **Language:** Python 3.11 or later.  
> * **Dependencies:** pandas only, if possible.  
> * **Structure:** a loader per input format that converts to one internal table; one pure function per check (two tables in, result out); an emitter. Same shape as designlint, so project scaffolding, CI and audit setup can be reused.  
> * **Overlapping domains:** a stated, documented rule for how overlapping hits on one protein are ordered when building an architecture. This choice affects results, so it must be explicit and configurable, not buried in code.  
> * **Configuration:** boundary tolerance, E-value cut-off and co-occurrence change threshold in one config file. Defaults documented with their reasoning.  
> * **Unknown situations:** a check that cannot run (for example, no scores in the input) reports "not evaluated" with the reason.

## **5\. Phases**

> 1. **Prior-art search (half a day, do this first).** Search GitHub, Bioconda, Bioconductor and the literature for tools that diff domain annotations or architectures across database versions. Outcome decides whether to build, contribute, or narrow the niche.  
> 2. **Format survey and setup (one day).** Confirm the exact column layout of each supported input format against the official documentation of the current tool versions. Create the project folder outside \~/Downloads, initialize git, CLAUDE.md, licence, pyproject.toml, pre-commit hooks.  
> 3. **Loaders and compatibility check (one to two days).** Parse both formats into the internal table. Test on truncated, malformed and mixed-format files.  
> 4. **Core checks (two to three days).** Gain and loss, boundary shifts, architecture changes, accession turnover, co-occurrence pairs, summary.  
> 5. **Tests (one to two days).** See section 6\.  
> 6. **Code review gate (half a day).** See section 7\.  
> 7. **Security audit gate (half a day).** See section 8\.  
> 8. **Documentation and release (half a day).** README with install, usage, meaning of each check, and limits. Tag version 0.1.0. Optional Zenodo archive for a DOI.  
> 9. **Case study (optional, for a paper).** Annotate one or a few public proteomes against two database releases and report how many architectures and co-occurrence pairs change. This is what would turn the tool into a publishable result. Running the annotation is the expensive step; whether a full proteome scan is practical on the laptop depends on the tool, the database and available RAM, which I have not checked. Precomputed public annotations may be an alternative.

## **6\. Testing**

> * **Synthetic pairs with known answers.** Small hand-built annotation tables where the exact differences are known: one domain removed, one boundary moved, two domains swapped, one accession renamed. Each check must report exactly those differences.  
> * **Identity test.** A run compared with itself must report zero differences.  
> * **Symmetry test.** Swapping the two inputs must swap gains and losses and change nothing else.  
> * **Order invariance.** Shuffling rows in either input must not change any result.  
> * **Independent cross-check.** Recompute gain, loss and pair counts for a small real example with a short separate script (or by hand) and confirm agreement.  
> * **Hostile-input tests.** Empty file, header only, wrong column count, non-numeric coordinates, start greater than end, negative positions, very long IDs, non-UTF-8 bytes, cells starting with \= \+ \- @. Expected: clean report or clean error, never a crash.  
> * **Scale test.** A synthetic table at realistic proteome size, to confirm memory and run time are acceptable on the laptop. No target number is set here because I have not measured one.  
> * **Coverage target:** 90% of lines in the checks module. A working target, not a standard.

## **7\. Comprehensive code review (release gate 1\)**

**Process**

> 1. Automated pass: ruff (lint and format) and mypy (types) with zero errors.  
> 2. Claude review in a fresh session that did not write the code, against the checklist below.  
> 3. Human domain review by you. The central questions are whether the architecture definition, the overlap rule and the co-occurrence counting are biologically sensible. This is your area of expertise and Claude cannot be the final authority on it.  
> 4. Optional second reviewer who uses InterProScan or Pfam regularly, to catch format assumptions that are wrong in practice.  
> 5. Every finding is logged as a GitHub issue and either fixed or closed with a written reason.

**Checklist**

> * Format correctness: do the loaders match the documented column layouts of the stated tool versions? Are those versions named?  
> * Definitions: are "domain", "architecture" and "co-occurrence" defined once, in the documentation, and implemented exactly as defined?  
> * Overlap handling: is the rule explicit, tested, and configurable?  
> * Rename handling: are renamed or merged accessions kept separate from real gains and losses when a mapping is supplied, and clearly labelled as "unresolved turnover" when it is not?  
> * Honesty of output: "not evaluated" instead of a false pass; differences described as differences, not as errors in either run.  
> * Coordinates: is the convention (1-based, inclusive) stated and consistent across formats? Off-by-one errors are the classic bug here.  
> * Edge cases: proteins with no domains, repeated domains, one protein only, runs with no shared proteins.  
> * Performance: no accidental quadratic loops over proteins or pairs.  
> * Simplicity: no features beyond version 1 scope.  
> * Tests: would they fail if a check were wrong?  
> * Documentation: limits section is accurate, including that version 1 cannot compare different gene models.

## **8\. Security audit (release gate 2\)**

**Threat model.** The tool reads untrusted tables (from collaborators, students, public repositories) and writes reports that may be opened in a spreadsheet or rendered as Markdown or HTML. The main risks are (a) malicious or malformed content in an input affecting the tool or the machine, (b) report output carrying that content into another program, and (c) unpublished genome annotation leaking through reports, logs or test fixtures.  
**Audit items**

> * **No execution paths.** Search for subprocess, os.system, eval, exec, shell=True and pandas eval or query on input-derived strings. Expected count in the released package: zero. The tool must not call InterProScan or HMMER itself.  
> * **Safe parsing.** Config loaded with a safe loader. No pickle. Plain-text tabular input only in version 1\.  
> * **Compressed input.** If gzip input is supported, cap the decompressed size to prevent decompression bombs. Simplest option: do not support compressed input in version 1\.  
> * **Resource limits.** Cap file size and row count with clear errors.  
> * **Spreadsheet formula injection.** In TSV or CSV output, neutralize cells beginning with \= \+ \- @, since protein IDs and descriptions come from untrusted input.  
> * **Markup injection.** Escape IDs and descriptions echoed into Markdown or HTML reports.  
> * **File handling.** Reads only the two files given (plus an optional mapping file); writes only to standard output or an explicit path; no overwrite without a flag; no paths taken from file contents.  
> * **No network, no telemetry.** Confirm and state it in the README.  
> * **Data confidentiality.** Test fixtures must be synthetic or from public proteomes with the source recorded. No unpublished genome annotation from collaborators in the repository or its history.  
> * **Third-party data licences.** If any real database-derived file is included as a fixture, confirm its licence allows redistribution. I have not checked the current licence terms of the relevant databases; verify before including anything.  
> * **Dependencies.** pip-audit against the lock file; versions pinned with hashes; dependency count minimal.  
> * **Static analysis.** Run bandit; triage every finding.  
> * **Secrets.** Run a secret scanner (for example gitleaks) over the full git history before going public.  
> * **CI hardening.** GitHub Actions with read-only default permissions, actions pinned to commit hashes, no secrets exposed to pull requests from forks.  
> * **Release integrity.** Two-factor authentication on GitHub and PyPI; prefer PyPI trusted publishing over a stored token.  
> * **Claude Code hygiene.** Never run Claude Code in \~/Downloads. Do not let it read CLAUDE.md or AGENTS.md from downloaded material. Keep downloaded annotation files in a separate folder and treat them as untrusted data.

**Output of the audit:** a SECURITY.md listing each item, the result, the date, and how to report a vulnerability.  
**Limit of this audit:** it is a self-audit with automated tools and an AI reviewer, not an independent penetration test. Proportionate for an offline tool of this size, but the README should not claim more than was done.

## **9\. Release criteria**

> * Prior-art search done and its outcome written in the README.  
> * All version 1 checks implemented, each with a synthetic pair of known answer.  
> * Identity, symmetry and order-invariance tests pass.  
> * Loaders verified against the official format documentation of named tool versions.  
> * Test suite passes in CI on macOS and Linux.  
> * Code review gate closed, including your domain review of the definitions.  
> * Security audit gate closed: SECURITY.md complete, zero execution paths, pip-audit and bandit clean or triaged.

## **10\. Risks and open questions**

> * **Novelty is unverified.** I have not searched for existing tools. General annotation-comparison tools exist; I am not certain whether one does domain-architecture and co-occurrence diffs across database versions. Phase 1 exists for this reason.  
> * **Input formats may differ from my assumptions.** The column layouts described here are from general knowledge and may have changed between tool versions. Phase 2 verifies them against current documentation.  
> * **Renames versus real change.** Without a mapping of renamed or merged accessions, turnover is overstated. Where to obtain such a mapping for a given pair of releases, and in what format, is an open question I have not verified.  
> * **Overlap rule is a scientific choice.** Different rules give different architectures. The default needs your judgement and a citation if one exists.  
> * **Scope creep to different gene models.** This is the most interesting use case and the one closest to your dissertation, but it requires protein matching by sequence and would move the project above designlint in complexity. Keep it for version 2\.  
> * **Publication.** The tool alone is likely a software note. The optional case study (phase 9\) is what would support a stronger paper.  
> * **Licence.** Not yet chosen. MIT is a common simple default; confirm against i3L policy.