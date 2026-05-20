# 2×2 factorial: bubble-aware × auto-deconvolution

## Why this experiment

`research/bubble_v2/04_validation/FINDINGS.md` showed that 4 of 5 supposedly
bubble-driven wins from the v1 technical note have **zero detected bubbles**
at any tested resolution (SH2B3/ATXN2, ANKS1B, TBX5/TBX3, SETD1B). The own
note already admits this in passing for SH2B3 — the win comes from how
bubble-aware windows interact with auto-deconvolution, not from bubble
structure.

This 2×2 factorial isolates the two contributions:

|                         | deconv OFF | deconv ON (default) |
|-------------------------|------------|---------------------|
| Uniform 10 kb           | cell U-OFF | cell U-ON           |
| Bubble-aware capped+weighted | cell B-OFF | cell B-ON      |

## What we expect to learn

The four cells let us read three independent effects:

1. **Bubble vs uniform (deconv fixed):** `B-ON − U-ON` and `B-OFF − U-OFF`.
   At loci where bubbles matter (TAS2R, centromere), bubble-aware should win
   under both deconv settings.

2. **Deconv vs no-deconv (windowing fixed):** `U-ON − U-OFF` and
   `B-ON − B-OFF`. At loci where deconv matters (the 4 zero-bubble loci),
   deconv-OFF should outperform deconv-ON regardless of windowing.

3. **Interaction:** `(B-ON − B-OFF) − (U-ON − U-OFF)`. If non-zero, bubble
   structure changes how deconvolution behaves — interesting on its own.

## Inputs

Reuses `00_inputs/` from `research/bubble_v2/` (which itself reuses
`paper_v2_submission/experiments/sim_expansion/`):

- `combined.paf` — 183 MB chimera PAF
- `tracts_n50.tsv` — ground-truth tracts for 50 chimeras
- `populations.tsv` — AFR/EUR sample map (46-haplotype ref panel)
- `queries.txt` / `subset.txt` — 50 chimera haplotype IDs

Plus two derived window BEDs:

- `windows/uniform.bed` — 13,328 uniform 10 kb windows (baseline)
- `windows/bubble_capped.bed` — 13,286 windows (uniform full + capped
  bubbles, max 50 kb), with parallel `weights.tsv`

## Cost

- 2× IBS via `impg similarity` (slow step, minutes-to-hours per window set).
- 4× ancestry decoding (fast, seconds per query × 50 queries × 4 = ~minutes).
- 1× scoring (Python, seconds).

The IBS computation is the bottleneck. The bubble-aware variant has slightly
fewer windows (13,286 vs 13,328) so it's marginally cheaper, not more.

## Reproduction

```bash
cd research/bubble_v2/experiments/deconv_2x2
./run.sh
```

The runner is idempotent: existing IBS / ancestry outputs are not recomputed.
To force a clean run, `rm -rf out/`.

## Output

A `concordance_2x2.tsv` with one row per locus and one column per cell.
Targets from `FINDINGS.md`:

| Locus class           | Expected outcome                                          |
|-----------------------|-----------------------------------------------------------|
| TAS2R_core            | B-OFF and B-ON both win over U-OFF and U-ON (bubble effect) |
| Centromere            | B-OFF ≫ U-OFF (capping+weighting recovery)               |
| SH2B3/ANKS1B/TBX5/SETD1B | (U-ON, B-ON) ≈ each other; (U-OFF, B-OFF) ≈ each other; pair-vs-pair gap = deconv effect, not bubble effect |
| WHOLE_CHR12 mean      | B-OFF + B-ON should beat both U cells if the centromere fix delivers |

## Status

`run.sh` is **written, not executed** in this checkpoint. Execution is
deferred to a session where the slow IBS step can run in the foreground
(or backgrounded via Bash run_in_background). See parent `STATE.md`.
