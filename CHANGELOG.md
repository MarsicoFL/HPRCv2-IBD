# Changelog

All notable changes to `impopk` are documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and the project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.2] — 2026-05-20

### Changed

- Internal cleanup of in-package documentation: comments and docstrings no
  longer reference the project's internal development branch. The feature
  set and behaviour are unchanged.
- README and CLI flag descriptions no longer link to research artifacts
  outside the package; the awk one-liner that builds the weights TSV from
  `impg depth --combined-output` is the canonical user path.
- Module-level docs use generic depth examples (e.g. "a window where 90% of
  the panel aligns") instead of figures tied to a specific reference panel.

## [0.2.1] — 2026-05-20

### Added

- `try_apply_window_weights(...) -> Result<()>` — non-panicking variant of
  `apply_window_weights` for downstream library users that cannot guarantee
  matched-length inputs. The panicking version is retained for callers that
  already validate inputs.
- `WindowWeights::has(chrom, start, end)` distinguishes "present in the
  table with weight 1.0" from "not present, defaulting to 1.0".
- `WindowWeights::build_vector(observations) -> (Vec<f64>, usize)` returns
  both the per-window weight vector and the count of windows that fell back
  to the default.
- When more than 10% of windows in the input similarity file have no
  matching row in `--window-weights`, `ancestry` emits a one-shot warning
  suggesting a coordinate-system or grid mismatch.
- End-to-end smoke test `data/examples/ancestry/run_weights.sh` exercising
  the binary with the new flag.

### Changed

- `WindowWeights` internal layout: nested `HashMap<chrom, HashMap<(start, end), weight>>`
  instead of a single `HashMap<(String, u64, u64), weight>`. Lookups no
  longer allocate a `String`.
- `WeightMode::Mult` docstring clarified: after the forward-backward
  normalization the effect is similar to `Interp` at the same `w`, with
  slightly less pull toward the prior.
- `.gitignore`: `paper/manuscript.pdf` and `paper/figures/*.pdf` overrides
  removed in favour of `/*.pdf` (repo-root only); tracked PDFs in any
  subdirectory are now picked up by default.

## [0.2.0] — 2026-05-20

### Added

- `ancestry --window-weights FILE` and `--weight-mode {interp|mult}`:
  per-window confidence weights from a TSV (`chrom, start, end, weight ∈ [0, 1]`)
  shrink the per-window log-emission toward a uniform prior when alignment
  support is low. Targets centromeres and reference-absent structural
  haplotypes where a uniform-window HMM otherwise propagates noisy
  identity calls.
- `WindowWeights`, `WeightMode`, `apply_window_weights`,
  `weights_for_observations`, `HasWindowKey` exported from
  `impopk_ancestry_cli` for downstream Rust users.
- Paper: new "Structural-haplotype-aware ancestry at the TAS2R cluster"
  Results subsection and accompanying figure (`paper/figures/fig_tas2r_bubble`).
  Methods section gains a paragraph describing the depth-derived weighting
  and the `--window-weights` interface.

### Notes

- `--window-weights` is opt-in: omitting it is a strict no-op.
- The transform is applied at the end of the emission pipeline, after every
  existing transform (label smoothing, ZCA whitening, kurtosis weighting, …)
  and immediately before forward-backward / Viterbi.
- Combining `--window-weights` with `--mask-bed` issues a warning. The two
  are intended as alternatives: a masked window already has uniform
  emissions and weighting is a no-op there.

## [0.1.0] — Initial release

- `ibs`, `ibd`, `ancestry`, `jacquard` Rust binaries.
- HMM-based local ancestry, IBD detection, kinship estimation from
  pangenome-derived pairwise identity, without phased VCFs.
- Validated against RFMix and hap-ibd on chr12 simulations, the CEPH 1463
  platinum pedigree, a bovine super-pangenome, and BXD recombinant inbred
  mice.
