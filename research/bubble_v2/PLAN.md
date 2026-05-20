# PLAN — bubble_v2

## Experiment matrix

Every row reuses the chr12 n=50 AMR 70 % AFR / 30 % EUR chimera benchmark
from `bubble_aware/`. Identical IBS source (`impg similarity` on combined
PAF), identical ref panel of 46 haplotypes, identical scoring script.

| ID  | Pipeline                                  | Centromere %     | chr12 mean %     | TAS2R %       |
|-----|-------------------------------------------|------------------|------------------|----------------|
| E0  | uniform 10 kb auto (v1 baseline)          | ≈ 92             | 90.50            | 80 (failing)   |
| E1  | bubble-aware auto (v1 published)          | **53.70**        | 93.55            | 100            |
| E2  | bubble-aware + (A) capping @ 50 kb        | **? — H1 test**  | ?                | should hold 100 |
| E3  | bubble-aware + (B) support weighting      | **? — H2 test**  | ?                | should hold 100 |
| E4  | bubble-aware + (A) + (B) combined         | **? — H1+H2**    | ?                | should hold 100 |
| E5  | E4 + (H3) pop_coverage feature            | ?                | ?                | should improve  |
| E6  | E4 across resolutions {1k,500,200,100,50} | sweep            | sweep            | sweep          |

## Success criteria

- **Minimum useful (H1+H2):** centromere ≥ 85 %, TAS2R unchanged at 100 %,
  chr12 mean ≥ 95 % (catches up to RFMix at 95.33).
- **Paper-grade contribution:** chr12 mean ≥ 96 % (parity with FLARE), zero
  TAS2R regression, centromere ≥ 90 %.
- **Novel finding (H3):** at least one locus outside chr12 where pop_coverage
  catches an ancestry call no other method makes (analog of TAS2R but
  independently discovered).

## Order of execution

Three things can run in parallel:

1. **Validation pipeline (no code changes).** Compute bubble stats + SD/SV
   overlap on existing 1000 bp and 200 bp depth runs. Result: confidence
   that detected bubbles correspond to known SVs.
2. **(A) window capping.** Pure Python change in adaptive-windows builder.
   Output is a new BED. Can run end-to-end ancestry on it the same day.
3. **(B) support weighting.** Rust change in `ancestry-cli` to accept a
   per-window weights file and multiply log-emissions by `w_t`. Then plumb
   the depth column from impg depth output through.

(A) ships faster (Python only); (B) is the bigger change. Validation goes
ahead in parallel with no risk.

## Validation criteria

For the detected-bubbles-are-real claim:

1. **SD overlap (Vollger 2025).** ≥ 80 % of bubbles ≥ 10 kb overlap an SD
   interval in `data/masks/chm13v2.0_SD.bed`.
2. **SV overlap (HPRC Liao 2023, vcfbub).** Bubbles 200 bp – 10 kb enriched
   for SV calls. Target Jaccard ≥ 0.5 at matching length scales.
3. **Cross-resolution stability.** Bubbles ≥ 1 kb detected at 200 bp should
   be detected at 1000 bp (containment ≥ 0.8). Smaller bubbles can drop.
4. **Biology recovery.** TAS2R, KRT, SH2B3/ATXN2, SETD1B all detected as
   bubbles at every resolution.

If 1 or 4 fail, the bubble detection itself is suspect — stop and inspect.

## Open questions for Franco

- Does the chr12 ground-truth tracts file from v1 live in `bubble_aware/`
  or `paper_v3_submission/`? Need exact path to symlink into 00_inputs/.
- HPRC SV catalog (Liao 2023) — is the VCF already downloaded somewhere,
  or do we need to fetch it from the HPRC FTP?
- For (B) Rust change: is it OK to add a flag to ancestry-cli on a feature
  branch, or do you want a parallel binary `ancestry-bubble` so the main CLI
  stays stable?
