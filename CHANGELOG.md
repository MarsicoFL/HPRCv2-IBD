# Changelog

All notable changes to `impopk` are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.1] — 2026-05-20

### Added

- `try_apply_window_weights(...) -> Result<()>` — non-panicking variant of
  `apply_window_weights` for downstream library users that cannot guarantee
  matched-length inputs. The panicking version is retained for callers that
  do.
- `WindowWeights::has(chrom, start, end)` and
  `WindowWeights::build_vector(observations) -> (Vec<f64>, usize)` — the
  latter returns both the per-window weight vector and the count of windows
  that fell back to the default 1.0 weight.
- CLI: when more than 10% of windows in the input similarity file have no
  matching row in `--window-weights`, `ancestry` now emits a one-shot warning
  that suggests a coordinate-system or grid mismatch (previously a silent
  no-op).
- CI: `data/examples/ancestry/run_weights.sh` smoke test asserts that
  `--window-weights` with all-1.0 weights is bit-equivalent (sorted) to the
  no-flag baseline, and that a low-coverage block runs cleanly. Wired into
  the `ancestry-determinism` workflow.
- README: inline one-liner that builds the weights TSV directly from
  `impg depth --combined-output`, removing the dependency on
  `research/bubble_v2/` for the documented user path.

### Changed

- `WindowWeights` internal layout: nested `HashMap<chrom, HashMap<(start,end), w>>`
  instead of a single `HashMap<(String,u64,u64), w>`. `get()` no longer
  allocates a `String` per lookup. No behavioural change.
- Doc tightening: `WeightMode::Mult` semantics (post-FB normalization makes
  it behave similarly to `Interp` at the same `w`; documented explicitly).
- `.gitignore`: the `paper/manuscript.pdf` + `paper/figures/*.pdf` overrides
  were replaced by a `/*.pdf` repo-root-only ignore, so tracked PDFs in any
  subdirectory no longer need explicit allowlist entries.

### Tests

- 16 unit tests in `window_weights::tests` (3 new: `try_apply_returns_error_on_length_mismatch`,
  `has_distinguishes_present_from_default`, `lookup_does_not_allocate_string`).
- 6 integration tests in `window_weights_integration_tests.rs` unchanged.
- 1 end-to-end CI test exercising the actual binary with the new flag.
- Full `ancestry-cli` test suite: 737 unit + integration, all passing.

### Internal review

All 12 issues from the v0.2.0 SE review are addressed (medium #1–#4 in this
release; minor #5–#9 and cosmetic #10–#12 either resolved or documented).

## [0.2.0] — 2026-05-20

### Added

- **`ancestry`: per-window support weighting for structurally variable regions.**
  Two new CLI flags accept a per-window confidence table derived from
  pangenome depth and shrink the per-window log-emission toward a uniform
  prior when alignment support is low. Targets centromeres and reference-absent
  structural haplotypes where a uniform-window HMM otherwise propagates noisy
  identity calls.
  - `--window-weights FILE` — TSV with `chrom, start, end, weight ∈ [0, 1]`.
    Windows not listed default to weight 1.0. Compatible with the output of
    `research/bubble_v2/01_centromere_fix/window_capping/make_capped_windows.py`.
  - `--weight-mode interp|mult` (default `interp`) — interpolation toward
    uniform `log(1/K)` (default) or multiplicative shrinkage.
- `WindowWeights`, `WeightMode`, `apply_window_weights`,
  `weights_for_observations`, `HasWindowKey` exported from
  `impopk_ancestry_cli` for downstream Rust users.
- Paper: new "Structural-haplotype-aware ancestry at the TAS2R cluster"
  Results subsection and accompanying figure (`paper/figures/fig_tas2r_bubble`).
  Methods section gains a paragraph describing the `(d_t, w_t)` weighting and
  the `--window-weights` interface.
- `research/bubble_v2/` directory: hypotheses, plan, state, validation tooling,
  the 2×2 deconvolution experiment runner, and the chr12 bubble-stats /
  SD-overlap / biology-recovery outputs.

### Tests

- 13 unit tests for `window_weights` (parsing, modes, edge cases).
- 6 integration tests including the centromere analog (adversarial 10-window
  block is normally followed by Viterbi; down-weighting it makes the path
  monotone).
- Full `ancestry-cli` test suite: 737 unit + integration, all passing.
- Clippy clean on the release profile.

### Notes

- `--window-weights` is opt-in: omitting it is a strict no-op.
- The transform is applied at the end of the emission pipeline, after every
  existing transform (label smoothing, ZCA whitening, kurtosis weighting, …)
  and immediately before forward-backward / Viterbi.
- Combining `--window-weights` with `--mask-bed` issues a warning. The two
  are intended as alternatives: a masked window already has uniform emissions
  and weighting is a no-op there.

## [0.1.0] — Initial release

- `ibs`, `ibd`, `ancestry`, `jacquard` Rust binaries.
- HMM-based local ancestry, IBD detection, kinship estimation from
  pangenome-derived pairwise identity, without phased VCFs.
- Validated against RFMix and hap-ibd on chr12 simulations, CEPH 1463
  platinum pedigree, bovine super-pangenome, and BXD recombinant inbred
  mice.
