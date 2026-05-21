# Security Policy

## Supported Versions

`impopk` follows [Semantic Versioning](https://semver.org/). Security fixes
land on the latest `0.x` minor; we do not patch older minor lines.

| Version | Supported |
|---------|-----------|
| 0.2.x   | ✅        |
| < 0.2   | ❌        |

## Reporting a vulnerability

Please report security issues privately. Two acceptable channels:

1. GitHub: open a private security advisory via
   [Security → Advisories → "Report a vulnerability"](https://github.com/MarsicoFL/IMPOPk/security/advisories/new).
2. Email the maintainer at the address listed on the repository owner's
   GitHub profile, with `[impopk security]` in the subject.

Please include:

- Affected version (commit hash if from `main`, or release tag).
- A minimal reproducer.
- The impact you observed and any suspected attack surface.

We will acknowledge receipt within 72 hours and aim to publish a fix —
or a clear non-issue determination — within 14 days. Fixes ship as a
patch release on the supported minor line; the advisory is published
on GitHub once a tagged release is available.

## Scope

`impopk` reads genomic data files (PAF, AGC, TSVs, BEDs) and a small set
of CLI flags. The expected threat model is:

- Untrusted input files. Parsing must not panic, segfault, or
  arbitrary-execute on malformed input. Resource exhaustion (memory or
  CPU blow-up on adversarial inputs) is in-scope.
- Untrusted command-line arguments. The CLI must reject malformed input
  cleanly with a non-zero exit code.

Out of scope:

- Bugs in upstream crates (report to those projects; we will pull in
  patched versions in our next release).
- Bugs in optional external tools (`impg`, `agc`) when invoked as a
  subprocess.
- Cryptographic claims: `impopk` does no cryptography.

## Audit cadence

The dependency tree is checked against the
[RustSec Advisory Database](https://rustsec.org/) via `cargo audit` at
each release; the latest baseline shipped with v0.2.2 reported 0
vulnerabilities across 66 transitive dependencies.
