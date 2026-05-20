# STATE — bubble_v2

**Last updated:** 2026-05-20

## Current iteration

Iteration 1 — implementation & experimental setup complete; data run pending.

## Status

| Step                                       | Status                                   |
|--------------------------------------------|------------------------------------------|
| Directory skeleton                         | done                                     |
| README / PLAN / STATE                      | done                                     |
| Locate n=50 chimera inputs                 | done                                     |
| Bubble stats @ 1000 bp / 200 bp            | done                                     |
| SD overlap validation                      | done                                     |
| Biology recovery check                     | done                                     |
| (A) window capping script                  | done                                     |
| (B) support weighting — design doc         | done                                     |
| (B) support weighting — Rust impl          | done (branch `bubble_v2` of IMPOPk)      |
| (B) unit tests                             | 13 passing                               |
| (B) integration tests                      | 6 passing                                |
| (B) full ancestry-cli test suite           | 737 unit + many integration; all green   |
| (B) clippy                                 | clean                                    |
| (B) end-to-end CLI smoke test              | passes (uniform 20-window toy)           |
| Auto-deconvolution 2×2 experiment runner   | written (`experiments/deconv_2x2/`)      |
| n=50 benchmark with E4 = A+B               | ready to run via `deconv_2x2/run.sh`     |
| Higher-resolution impg depth runs          | not needed                               |
| (H3) pop_coverage feature                  | blocked on n=50 results                  |
| Cross-chrom replication                    | blocked on n=50 results                  |
| IMPOPk sync report                         | pending (scope shifted, see below)       |

## Headline findings (iteration 0)

1. At `--min-interval-len 1000`, chr12 has 55 bubbles. The 10–50 kb tier is
   75 % SD-supported. ≥ 50 kb tier is all in the centromere and 0 % SD-supported
   (expected — centromere is satellite, not SD).
2. **TAS2R is the only genuinely bubble-driven publication case.** SH2B3/ATXN2,
   ANKS1B, TBX5/TBX3, SETD1B have **zero detected bubbles** at any tested
   resolution. Those "wins" in the v1 technical note are auto-deconvolution
   artifacts, not bubble structure.
3. The centromere (36.4 – 37.2 Mb) is **over-detected**: 5 bubbles at 1000 bp
   covering 768 kb, depth dropping to 22/234. At 200 bp it's 750 bubbles
   covering 771 kb, depth as low as 11/234. The signal is real but noisy.

## Iteration 1 deliverables

- `src/ancestry-cli/src/window_weights.rs` — module with `WindowWeights`,
  `WeightMode`, `apply_window_weights`, `weights_for_observations`, plus the
  `HasWindowKey` trait (implemented for `AncestryObservation` in `hmm.rs`).
- `src/ancestry-cli/src/main.rs` — CLI flags `--window-weights FILE` and
  `--weight-mode interp|mult`, wired into both the pass-1 and pass-2
  emission pipelines.
- `src/ancestry-cli/tests/window_weights_integration_tests.rs` — 6 property
  tests including the "noisy block becomes monotone" centromere analog.
- `experiments/deconv_2x2/` — runner + scorer for the 2×2 factorial that
  isolates the bubble effect from the auto-deconvolution effect.

## Next step

Execute the 2×2 factorial via `experiments/deconv_2x2/run.sh`. Two `impg
similarity` passes (uniform + bubble-aware windows) are the slow step; four
ancestry decodings and scoring are fast. Output → `out/concordance_2x2.tsv`.

## Decisions waiting on Franco

- IMPOPk sync report scope. With `bubble_v2` already pushed to IMPOPk, the
  remaining question is whether to backfill the missing cycles-82-92 tests
  from HPRCv2-IBD into IMPOPk main. Not blocking.
- Are these defaults right for the 2×2 run, or should the optional features
  (`--auto-configure`, `--estimate-params`, mask BED) be threaded through
  consistently?

The biology-recovery table in `04_validation/FINDINGS.md` is the suggested
unit test: bubble_v2 must win at TAS2R, tie everywhere bubbleless, and
recover the centromere to ≥ 85 %.
