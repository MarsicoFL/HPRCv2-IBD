# Changelog

## [0.2.2] — 2026-05-20

- Cleanup of in-package docs and comments.
- README example for the depth-derived weights TSV.

## [0.2.1] — 2026-05-20

- `try_apply_window_weights` (non-panicking variant of `apply_window_weights`).
- `WindowWeights::has` and `WindowWeights::build_vector(observations) -> (Vec<f64>, usize)`.
- Warning when `--window-weights` matches < 90% of input windows.
- Nested `HashMap` layout: `WindowWeights::get` no longer allocates.
- End-to-end smoke test `data/examples/ancestry/run_weights.sh`.

## [0.2.0] — 2026-05-20

- `ancestry --window-weights FILE --weight-mode {interp,mult}`:
  per-window confidence weights from a TSV (`chrom, start, end, weight ∈ [0, 1]`)
  shrink the log-emission toward uniform where alignment support is low.
- Public exports: `WindowWeights`, `WeightMode`, `apply_window_weights`,
  `weights_for_observations`, `HasWindowKey`.
- Paper: TAS2R case-study section + figure.

`--window-weights` is opt-in; omitting it is a no-op.

## [0.1.0]

- `ibs`, `ibd`, `ancestry`, `jacquard` binaries: HMM-based IBD detection,
  local ancestry, kinship from pangenome-derived pairwise identity.
- Validated against RFMix and hap-ibd on chr12 simulations, the CEPH 1463
  platinum pedigree, a bovine super-pangenome, and BXD recombinant inbred
  mice.
