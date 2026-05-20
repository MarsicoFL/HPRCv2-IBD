# (A) Window capping — implementation

## What changed vs v1

`make_capped_windows.py` differs from `bubble_aware/.../make_adaptive_windows.py`
in three ways:

1. **Bubble size cap.** Any bubble window with length ≥ `--max-bubble-bp`
   (default 50 kb) is subdivided into uniform sub-windows. Sub-windows
   inherit the parent bubble's depth and support fraction; only the
   `window_type` field is tagged `bubble_capped_*` so downstream code can
   distinguish.
2. **Per-window weights file.** A parallel TSV with one row per emitted
   window: `chrom, start, end, weight = depth / max_depth`. Full-depth
   windows get weight `1.0`. Bubbles inherit `depth / 234`. This file is
   the contract with the (B) support-weighted emission Rust change.
3. **Per-bubble parent ID.** Sub-windows from the same parent bubble share
   a `parent_bubble_id` (e.g. `b3`) so post-decoding analysis can re-merge
   them when needed.

## Default output (chr12, cap = 50 kb)

```
$ python3 make_capped_windows.py \
    --bubbles-bed bubble_aware/.../chr12_bubbles.bed \
    --grid-step 10000 --max-bubble-bp 50000 \
    --chrom-rename 'CHM13#0#chr12->chr12' \
    --out windows_cap50k.bed \
    --out-detailed windows_cap50k.detailed.tsv \
    --out-weights windows_cap50k.weights.tsv
```

| Metric                         | v1 uncapped | (A) cap 50 kb |
|--------------------------------|-------------|---------------|
| Full-depth windows             | 13,258      | 13,258        |
| Atomic bubble windows          | 14          | 11            |
| Subdivided bubble windows      |  —          | 17            |
| Bubbles capped (parent count)  |  —          |  3            |
| Total windows                  | 13,272      | 13,286        |

The 3 capped bubbles are exactly the 3 ≥ 50 kb intervals from the chr12
bubble stats table — all in the pericentromere (36.4 – 37.2 Mb). The
457 kb mega-bubble becomes 10 sub-windows of ≈ 45 kb each; the 210 kb
follow-up becomes 5; the 67 kb after that becomes 2. Total 17 new
sub-windows replace 3 atomic ones in the failure region.

## Parameter sweep (not run yet)

`--max-bubble-bp ∈ {25_000, 50_000, 100_000}` is cheap to sweep — pure
Python over an existing BED. Will run alongside the (B) Rust change once
that ships.

## What this script does NOT do

- It does **not** change the IBS computation. The downstream
  `impg similarity` / pipeline still computes one identity score per
  emitted window. Capped sub-windows inside a single bubble will receive
  *correlated but not identical* identity scores, because impg integrates
  over each sub-window's coordinates independently.
- It does **not** apply weights to anything. The weights TSV is just
  data; nothing in ancestry-cli reads it yet. That is (B).
