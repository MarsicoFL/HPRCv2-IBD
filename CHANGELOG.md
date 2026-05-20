# Changelog

All notable changes to `impopk` are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
