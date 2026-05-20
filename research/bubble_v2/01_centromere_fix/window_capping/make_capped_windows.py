#!/usr/bin/env python3
"""
Build adaptive windows with a bubble-size cap.

Differences vs v1's make_adaptive_windows.py:

  1. Bubble windows larger than --max-bubble-bp are subdivided into uniform
     sub-windows of --max-bubble-bp. Default 50 kb. Targets the centromere
     case where 5 atomic bubbles of 100–450 kb dominated the HMM posterior.

  2. Emits a parallel window-weights TSV (`<out>.weights.tsv`) with one row
     per window containing `weight = depth / max_depth`. This is what the
     (B) support-weighted emission step in ancestry-cli will consume.

Input: chr12_bubbles.bed from combined_to_bed.py (CHM13#0#chr12, start, end,
       window_type, depth, missing_count, ...).

Outputs:
  --out             standard BED (chrom, start, end) for IBS step
  --out-detailed    TSV with (chrom, start, end, window_type, depth,
                              support_frac, parent_bubble_id)
  --out-weights     TSV with (chrom, start, end, weight) — for (B)
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path


def load_bubbles_bed(path: Path):
    """Read combined_to_bed.py output, yielding (chrom, start, end, wtype, depth)."""
    with path.open() as fh:
        header = fh.readline()
        assert header.startswith("#"), f"expected '#' header, got: {header!r}"
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            yield parts[0], int(parts[1]), int(parts[2]), parts[3], int(parts[4])


def merge_adjacent_full(rows):
    """Merge adjacent 'full' rows into single intervals so the re-tiling
    produces a clean uniform grid (not a sequence of inherited fragments).
    """
    out = []
    for chrom, start, end, wtype, depth in rows:
        if out and out[-1][0] == chrom and out[-1][2] == start \
                and out[-1][3] == "full" and wtype == "full":
            prev = out[-1]
            out[-1] = (prev[0], prev[1], end, "full", prev[4])
        else:
            out.append((chrom, start, end, wtype, depth))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bubbles-bed", required=True, type=Path,
                    help="combined_to_bed.py output (chr12_bubbles.bed)")
    ap.add_argument("--grid-step", type=int, default=10_000,
                    help="window size in full-depth regions")
    ap.add_argument("--max-bubble-bp", type=int, default=50_000,
                    help="cap on atomic bubble window size; larger bubbles are "
                         "subdivided into sub-windows of this size")
    ap.add_argument("--n-total", type=int, default=234,
                    help="total sample count = max possible depth")
    ap.add_argument("--chrom-rename", default=None,
                    help="rename chromosome, e.g. 'CHM13#0#chr12->chr12'")
    ap.add_argument("--out", required=True, type=Path,
                    help="BED (chrom, start, end) — for IBS step")
    ap.add_argument("--out-detailed", required=True, type=Path,
                    help="detailed TSV with window_type, depth, support_frac")
    ap.add_argument("--out-weights", required=True, type=Path,
                    help="weights TSV (chrom, start, end, weight) for (B)")
    args = ap.parse_args()

    rows = merge_adjacent_full(load_bubbles_bed(args.bubbles_bed))

    rename_src, rename_dst = None, None
    if args.chrom_rename:
        if "->" in args.chrom_rename:
            rename_src, rename_dst = args.chrom_rename.split("->")
        else:
            rename_dst = args.chrom_rename

    def chrom_out(c):
        if rename_dst is None:
            return c
        if rename_src is None or c == rename_src:
            return rename_dst
        return c

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out_detailed.parent.mkdir(parents=True, exist_ok=True)
    args.out_weights.parent.mkdir(parents=True, exist_ok=True)

    n_full = 0
    n_bubble_atomic = 0
    n_bubble_subdiv = 0
    n_bubbles_capped = 0
    bp_full = 0
    bp_bubble = 0

    with args.out.open("w") as bed, \
         args.out_detailed.open("w") as det, \
         args.out_weights.open("w") as wts:

        det.write("chrom\tstart\tend\twindow_type\tdepth\tsupport_frac\t"
                  "parent_bubble_id\n")
        wts.write("chrom\tstart\tend\tweight\n")

        bubble_id = 0
        for chrom, b_start, b_end, wtype, depth in rows:
            out_chrom = chrom_out(chrom)

            if wtype == "full":
                # Re-tile full regions at grid_step
                pos = b_start
                while pos < b_end:
                    w_end = min(pos + args.grid_step, b_end)
                    bed.write(f"{out_chrom}\t{pos}\t{w_end}\n")
                    det.write(f"{out_chrom}\t{pos}\t{w_end}\tfull\t{depth}\t"
                              f"{depth / args.n_total:.4f}\t-\n")
                    wts.write(f"{out_chrom}\t{pos}\t{w_end}\t1.0000\n")
                    bp_full += w_end - pos
                    pos = w_end
                    n_full += 1
                continue

            # Bubble case
            length = b_end - b_start
            bubble_id += 1
            weight = depth / args.n_total
            if length <= args.max_bubble_bp:
                bed.write(f"{out_chrom}\t{b_start}\t{b_end}\n")
                det.write(f"{out_chrom}\t{b_start}\t{b_end}\tbubble_{wtype}\t"
                          f"{depth}\t{weight:.4f}\tb{bubble_id}\n")
                wts.write(f"{out_chrom}\t{b_start}\t{b_end}\t{weight:.4f}\n")
                bp_bubble += length
                n_bubble_atomic += 1
            else:
                # Subdivide into floor(length / max_bubble_bp) ceil-rounded sub-windows
                n_sub = math.ceil(length / args.max_bubble_bp)
                sub_size = math.ceil(length / n_sub)
                pos = b_start
                while pos < b_end:
                    w_end = min(pos + sub_size, b_end)
                    bed.write(f"{out_chrom}\t{pos}\t{w_end}\n")
                    det.write(f"{out_chrom}\t{pos}\t{w_end}\tbubble_capped_{wtype}\t"
                              f"{depth}\t{weight:.4f}\tb{bubble_id}\n")
                    wts.write(f"{out_chrom}\t{pos}\t{w_end}\t{weight:.4f}\n")
                    bp_bubble += w_end - pos
                    pos = w_end
                    n_bubble_subdiv += 1
                n_bubbles_capped += 1

    print(
        f"Full-depth windows:      {n_full:,} ({bp_full:,} bp)",
        f"Atomic bubble windows:   {n_bubble_atomic:,}",
        f"Subdivided bubble wins:  {n_bubble_subdiv:,} (from {n_bubbles_capped} large bubbles)",
        f"Total bubble bp:         {bp_bubble:,}",
        f"Total windows:           {n_full + n_bubble_atomic + n_bubble_subdiv:,}",
        sep="\n", file=sys.stderr,
    )


if __name__ == "__main__":
    main()
