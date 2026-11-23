# Security

## Reporting a vulnerability

Please report suspected vulnerabilities privately through GitHub's "Report a vulnerability" (private security advisory) on this repository, not as a public issue. Include the input that triggers the problem and the domaindiff version.

## Threat model

domaindiff reads untrusted tables (from collaborators, students or public repositories) and writes reports that may be opened in a spreadsheet or rendered as Markdown. The risks considered are:

- malicious or malformed input affecting the tool or the machine;
- report output carrying input content into another program;
- unpublished annotation leaking through the repository.

## Audit record, 23 November 2026, version 0.1.0

This is a self-audit with automated tools and an AI reviewer. It is not an independent penetration test.

| Item | Result |
|---|---|
| No execution paths | Pass. grep for subprocess, os.system, eval, exec, shell=True, pickle, pandas query, and network modules in `domaindiff/`: 0 matches. The tool never runs InterProScan or HMMER. |
| Safe parsing | Pass. Config is read with stdlib `tomllib` (data only), and unknown keys and wrong types are rejected. No pickle. Plain-text tables only. |
| Compressed input | Pass. Gzip is not supported; gzip magic bytes are rejected with an error. |
| Resource limits | Pass. `max_file_mb` (default 2048) and `max_rows` (default 10,000,000). Only regular files are accepted, so FIFOs and devices are refused. IDs are capped at 1000 characters. |
| Input validation | Pass. Protein IDs and accessions must be printable non-space ASCII, which rejects terminal escape sequences. Coordinates must be integers with 1 <= start <= end. Non-UTF-8 bytes, NaN and non-numeric scores raise clean errors. Tested in `tests/test_domaindiff.py`. |
| Spreadsheet formula injection | Pass. TSV cells beginning with `= + - @`, tab or CR are prefixed with `'`. Tested. |
| Markup injection | Pass. Markdown output backslash-escapes ASCII punctuation, including `< > | [ ]`. File names printed in reports have control characters replaced. Tested. |
| File handling | Pass. The tool reads only the two inputs plus the optional mapping and config files. It writes to stdout or `-o`. `-o` uses exclusive create unless `--force` is given, and it never writes to an input path. No path is taken from file contents. |
| No network, no telemetry | Pass. There are no network imports, and this is stated in the README. |
| Data confidentiality | Pass. All fixtures and examples are synthetic (hand-built accessions and protein IDs). No collaborator data is included. |
| Third-party data licences | Not applicable. No database-derived file is included. |
| Dependencies | Pass. There are no runtime dependencies. Dev tools are pinned with hashes in `requirements-dev.txt` (uv pip compile --generate-hashes). pip-audit 2.10.1 found no known vulnerabilities. |
| Static analysis | Pass. bandit 1.9.4 found no issues. ruff 0.17.6 (including the flake8-bandit `S` rules) and mypy 2.4.1 `--strict` report no errors. |
| Secrets | Pass. gitleaks 8.30.1 over the working tree found no leaks. It also runs as a pre-commit hook. Rerun `gitleaks git .` over the full history before making the repository public. |
| CI hardening | Configured. The workflow has `permissions: contents: read`. actions/checkout and actions/setup-python are pinned to commit SHAs (verified with `git ls-remote`). `persist-credentials: false` is set and no secrets are used. Not yet exercised on GitHub. |
| Release integrity | Open. Enable 2FA on GitHub and PyPI, and use PyPI trusted publishing if the package is published. The maintainer has to do this. |
| Claude Code hygiene | Pass. The project lives outside `~/Downloads`, and no downloaded material is in the repository. |

## Code review gate (same date)

An independent AI review in a separate session found seven defects. All were fixed, each with a regression test:

1. Score columns (ProSiteProfiles, HAMAP) were read as E-values.
2. TSV output dropped check statuses.
3. An analysis filter that matched nothing gave a misleading error.
4. Weak tests let seven single-line mutants survive. Six are now caught; the seventh changed how hits without an E-value were ranked, and `best_evalue` now refuses such hits instead.
5. The wrong HMMER flavour was accepted silently.
6. File names were printed without escaping.
7. This file was missing.

The scientific judgement calls it raised are listed under "Open questions for domain review" in the README. They still need the maintainer's domain review.
