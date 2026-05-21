<!-- One-line summary on the title; details below. -->

## Summary
<!-- what changes and why, in user-visible language. -->

## Verification
- [ ] `cargo build --release --workspace` clean
- [ ] `cargo test --workspace` passing
- [ ] `cargo clippy --workspace -- -D warnings` clean on stable
- [ ] New public surface has rustdoc (see `window_weights` for the bar)
- [ ] `CHANGELOG.md` updated with a user-visible entry
- [ ] `Cargo.lock` still at format v3 (see `CONTRIBUTING.md` for the re-pin step if Cargo rewrote it)

## Related
<!-- linked issues, prior PRs, paper / data sources. -->
