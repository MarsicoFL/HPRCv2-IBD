# Validation findings — iteration 0

**Date:** 2026-05-20
**Inputs:** `impg depth` runs from `bubble_aware/` at `--min-interval-len`
1000 bp and 200 bp; Vollger-style SD catalog at `data/masks/chm13v2.0_SD.bed`
(83,516 intervals).

## 1. Bubble counts and length distribution (chr12, ref-anchored)

| Resolution | Bubbles | Bubble territory | Median len | Max len      | ≥10 kb | ≥50 kb |
|------------|---------|------------------|------------|--------------|--------|--------|
| 1000 bp    |    55   | 951 kb (0.71 %)  | 2,778 bp   | **456,957 bp** |  7    |  3    |
| 200 bp     | 1,551   | 1,850 kb (1.39 %) |   678 bp   |  21,620 bp   |  9    |  0    |

Key observation: increasing detection resolution **fragments the centromere**
from a few mega-bubbles into hundreds of small ones. At 200 bp no single
bubble exceeds 50 kb. The total bubble territory roughly doubles, but most
of the new bubbles are < 1 kb.

## 2. SD overlap by length tier (Vollger catalog)

| Resolution | Tier         | n      | SD-supported | % SD |
|------------|--------------|--------|--------------|------|
| 1000 bp    | ≥ 50 kb      |     3  |     0        |  0   |
| 1000 bp    | 10–50 kb     |     4  |     3        | 75   |
| 1000 bp    | 1–10 kb      |    48  |    14        | 29   |
| 200 bp     | 10–50 kb     |     9  |     2        | 22   |
| 200 bp     | 1–10 kb      |   541  |    25        |  4.6 |
| 200 bp     | < 1 kb       | 1,001  |    40        |  4.0 |
| **both**   | **centromere(36–38 Mb)** | many | **0** |  **0** |

Take-aways:

- At 1000 bp the 10–50 kb tier is **75 % SD-supported** — strong validation
  that bubbles at the right length scale correspond to known structural
  duplications.
- The centromere bubbles have **0 % SD overlap**. This is expected: the
  centromere is alpha-satellite / HSAT, not segmental duplication. The
  Vollger catalog is not the correct reference for centromere validation —
  a T2T-CHM13 censat track would be.
- Below 1 kb, SD overlap collapses to ≈ 4 %. Most sub-1 kb "bubbles" are
  wfmash alignment noise, not biology.

**Implication:** use `--min-interval-len 1000` and filter to bubbles
≥ 1 kb. That leaves 55 chr12 bubbles. Going to 200 bp produces 1,551 but
the additional 1,500 are mostly noise.

## 3. Biology recovery — which named loci are bubble-driven?

| Locus                  | Span (kb) | 1000 bp     | 200 bp        | Verdict        |
|------------------------|-----------|-------------|---------------|----------------|
| TAS2R full             |   200     | 1 (4.5 kb)  | 10 (20 kb)    | ✅ bubble-driven |
| **TAS2R core**         |     6     | 1 (4.5 kb, depth 218) | 6 (17 kb) | ✅ bubble-driven |
| **Centromere alpha**   |   800     | 5 (768 kb)  | 750 (771 kb)  | ⚠️ over-detected |
| KRT type-II            |   200     | 2 (12 kb)   | 4 (12 kb)     | ✅ bubble-driven |
| Subtelomere SD         |   500     | 2 (7 kb)    | 9 (12 kb)     | ✅ bubble-driven |
| ANKS1B q23.1           | 2,180     | **0**       | **0**         | ❌ NOT bubble-driven |
| SH2B3 / ATXN2          |   340     | **0**       | **0**         | ❌ NOT bubble-driven |
| TBX5 / TBX3            |   460     | **0**       | **0**         | ❌ NOT bubble-driven |
| SETD1B q24.32          | 1,530     | **0**       | **0**         | ❌ NOT bubble-driven |

This re-reads four of the five extended-scan "wins" from the v1 technical
note: SH2B3/ATXN2, TBX5/TBX3, ANKS1B, SETD1B all have **zero detected
bubbles** at either resolution. The technical note already acknowledged
this for SH2B3 in passing — those wins come from the auto-deconvolution
toggle, not from local bubble structure.

**Only TAS2R survives as a genuinely bubble-driven publication case.**

## 4. Validation criteria for bubble_v2

The biology-recovery table doubles as a unit test for any new pipeline:

| Locus class           | Expected behavior                                |
|-----------------------|--------------------------------------------------|
| Bubble present (TAS2R, KRT, Subtelo) | bubble_v2 wins or ties vs uniform |
| Bubble noisy (centromere) | bubble_v2 (with A+B) recovers to ≥ 85 %     |
| No bubble (SH2B3, ANKS1B, …) | bubble_v2 indistinguishable from uniform |

If bubble_v2 wins at SH2B3 despite zero detected bubbles, something
non-bubble-related is doing the work (likely auto-deconvolution interaction).

## 5. Open questions

- The centromere validation needs T2T-CHM13 censat / HSAT annotation.
  Available at <https://github.com/marbl/CHM13/issues/47> or
  the T2T-CHM13 v2.0 release. Worth checking what's in `data/masks/`.
- Cross-resolution stability (Jaccard of bubble territory at 1000 bp vs
  200 bp) not yet computed but expected to be high for ≥ 1 kb bubbles.
- HPRC SV catalog (Liao 2023) overlap pending — useful for the 1–10 kb
  tier where Vollger SDs alone give only 29 % support.
