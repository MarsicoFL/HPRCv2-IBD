# bubble_v2 — bubble-aware ancestry, round 2

Continuation of `bubble_aware/` (one-off, closed 2026-04-23) as a structured
research line aimed at three concrete gaps the technical note left open.

## Why round 2

The chr12 v1 result (`bubble_aware/manuscript/technical_note.md`):

- Bubble-aware impop_k auto: **93.55 %** mean concordance (n=50 AMR chimeras),
  +3.06 pp over uniform-10 kb auto (CI [+2.05, +4.05], 39/50 wins).
- TAS2R core (chr12:10.95–10.98 Mb): **9/9 vs RFMix 0/9** on AFR donor HG01960
  — the publication case.
- Centromere (chr12:36.4–37.2 Mb): **53.7 %** vs RFMix/impop_opt ≈ 92 %.
  Five atomic bubbles of 100–450 kb let alpha-satellite noise dominate the
  posterior across 762 kb.

Round 2 attacks the centromere failure without sacrificing TAS2R, and adds a
novel angle (continuous depth as ancestry-informative feature) that Shuo/Erik
have not yet considered.

## Hypotheses, in order

1. **H1 — Window capping (A).** Subdividing any bubble window ≥ 50 kb into
   50-kb sub-windows turns the centromere into ≈ 15 observations instead of
   5. Single bad sub-windows no longer dominate.

2. **H2 — Support weighting (B).** Emission weight `w = mean_depth / max_depth`
   per window. Centromere bubbles drop to weight ≈ 0.09 (depth 22/234), HMM
   defaults to flanking interpolation in those windows.

3. **H1 + H2 compose.** Each addresses a different failure mode (too few
   observations vs untrusted observations). The expected effect is multiplicative
   recovery on the centromere.

4. **H3 — Continuous depth as population-informative signal.** For each window,
   `pop_coverage_k = (# samples of pop k aligning) / |pop k|` is a K-vector.
   At TAS2R the AFR-specific structural haplotype produces `pop_coverage` ≈
   `(1.0, 0.5, …)` — that pattern alone is ancestry-informative even before
   looking at identity. No VCF method can access this signal.

## Directory layout

```
research/bubble_v2/
├── README.md            # this file
├── PLAN.md              # detailed experiments, success criteria
├── STATE.md             # current iteration / next step
├── 00_inputs/           # symlinks to PAF, AGC, sample lists, ref panel, GT tracts
├── 01_centromere_fix/
│   ├── window_capping/      # (A) — adaptive windows with bubble size cap
│   └── support_weighting/   # (B) — Rust change in ancestry-cli emissions
├── 02_depth_continuous/     # (H3) per-window pop_coverage feature
├── 03_resolution_sweep/     # impg depth at --min-interval-len {1000,500,200,100,50}
├── 04_validation/           # bubble stats + SD/SV catalog overlap + biology check
└── 05_cross_chrom/          # chr1, chr2 (LCT), chr20, chr22 replication
```

## Inputs

Reuses everything that already worked in v1:

- impg fork (Shuo's `fix/depth-review-p0-p1`): already compiled in
  `bubble_aware/tools/impg_shuo/target/release/impg`.
- n=50 AMR chimera benchmark: chr12 PAF, 46-haplotype ref panel, ground-truth
  tracts. Paths documented in `00_inputs/` (symlinks, not copies).

## Reproduction

Each experiment directory has its own runner (`run.sh`) and writes outputs
inside its own directory. The top-level state lives in `STATE.md`.

## What is and is not in scope

In scope:
- Algorithmic changes to ancestry-cli (per-window emission weights).
- Bubble window construction logic (capping, depth-driven splits).
- Validation against Vollger 2025 SD catalog and HPRC SV catalog (Liao 2023).

Out of scope:
- Touching `bubble_aware/` (frozen historical record).
- IBS extraction beyond `impg similarity` — `ibs-from-paf` is retired.
- Cross-population transfer / cross-chromosome parameter sharing (separate
  research line in `docs/research/`).
