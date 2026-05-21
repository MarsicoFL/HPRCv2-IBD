# Contributing to impopk

Thanks for considering a contribution. This document covers the practical
details of getting a change in. For security-sensitive issues see
[`SECURITY.md`](SECURITY.md) instead — do not file a public issue.

## Quick checklist

- `cargo build --release --workspace` succeeds.
- `cargo test --workspace` passes.
- `cargo clippy --workspace -- -D warnings` is clean on `stable`.
- The pre-existing tests still pass; new tests cover any new public surface.
- `CHANGELOG.md` has an entry under the next unreleased section, written in
  user-visible language (no internal codenames, no review-process metadata).
- New public items have rustdoc.

CI runs the same checks on `stable` and on the declared MSRV (`1.74` as of
v0.2.x); both must pass before a PR can merge.

## Working tree layout

```
src/
├── common/        impopk-common — shared types
├── ibs-cli/       ibs binary    — wraps `impg similarity`
├── ibd-cli/       ibd binary    — 2-state IBD HMM
├── ancestry-cli/  ancestry      — N-state local-ancestry HMM
├── jacquard-cli/  jacquard      — 9 Δ-coefficient kinship
└── argraph-cli/   argraph       — experimental ARG inference (v0.1)
data/examples/     small precomputed inputs + ready-to-run recipes
paper/             LaTeX manuscript bundled with the release
```

## Local development

```bash
# Build everything
cargo build --release --workspace

# Run a single crate's tests (fast feedback)
cargo test -p impopk-ancestry-cli

# Run the lint gate
cargo clippy --workspace -- -D warnings

# Run all four example smoke tests
bash data/examples/ibd/run.sh
bash data/examples/ancestry/run.sh
bash data/examples/ancestry/run_weights.sh
bash data/examples/pedigree/run.sh
```

`CARGO_BUILD_JOBS=2` is recommended on memory-constrained machines —
the workspace builds cleanly with two compile jobs.

## Cargo.lock format

The committed `Cargo.lock` is at format **version 3** so the MSRV-1.74 CI
job can parse it. Newer local Cargo (≥ 1.78) silently rewrites the file
to version 4 on `cargo update` and similar commands. If your change adds
or updates a dependency:

1. Run `cargo update -p <changed-crate>` as usual.
2. Re-pin the lock: `sed -i 's/^version = 4$/version = 3/' Cargo.lock`.
3. Verify CI passes both `Test (stable)` and `Test (1.74)`.

Alternatively, propose bumping the workspace MSRV in the same PR.

## Style

- Code that you author should be rustfmt-clean against the committed
  `rustfmt.toml`. We do not enforce `cargo fmt --check` in CI yet because
  the existing tree carries unrelated drift; please at least format the
  files you touch.
- Public items must have rustdoc. The `window_weights` module is the
  reference: see `src/ancestry-cli/src/window_weights.rs`.
- New CLI flags follow `--kebab-case`. Document them in their `#[arg(...)]`
  doc-comment so they appear in `--help` automatically.

## Tests

Three layers:

1. **Module-internal `#[cfg(test)] mod tests`** — fast, in-process, no I/O.
2. **`tests/` integration tests** — black-box use of the public API,
   typically with `tempfile`.
3. **`data/examples/*/run.sh`** — end-to-end smoke against the built
   binary. Exercised in CI by the `ancestry-determinism` job.

A change that adds a public function should add at least one test at
layer 1; a CLI flag should add at least one at layer 3.

## Commit style

- One logical change per commit.
- Imperative subject ≤ 72 chars, no leading namespace prefix the reader
  can't grep on (e.g. `ancestry:`, `ibd:`, `ci:`, `docs:`, `chore:`).
- Wrap body at 72 columns; explain *why* in user-visible language.
- No AI attribution lines.

## Releasing

- Bump version in workspace `Cargo.toml` and every crate `Cargo.toml`.
- Re-pin `Cargo.lock` to version 3 (see above).
- Update `CHANGELOG.md`: move pending notes under the new release header.
- Tag `vX.Y.Z`. The merge to `main` and the tag push are independent.
