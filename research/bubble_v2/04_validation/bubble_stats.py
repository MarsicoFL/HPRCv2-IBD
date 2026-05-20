#!/usr/bin/env python3
"""
Compute bubble statistics from impg depth combined.bed output.

The combined.bed format from Shuo's fork is:
    seq_name <TAB> length <TAB> depth <TAB> samples_csv

It is sequential along the reference: start/end are reconstructed by
accumulating `length`. A bubble is any interval with depth < max_depth.
"""
from __future__ import annotations

import argparse
import csv
import statistics
import sys
from pathlib import Path


def parse_combined(path: Path):
    """Yield (chrom, start, end, depth, n_samples) for each row in the
    sequential combined.bed produced by `impg depth --combined-output`.
    """
    cursor = {}
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            length = int(parts[1])
            depth = int(parts[2])
            samples = parts[3].split(",") if len(parts) > 3 else []
            start = cursor.get(chrom, 0)
            end = start + length
            cursor[chrom] = end
            yield chrom, start, end, depth, len(samples)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined", required=True, type=Path,
                    help="impg depth combined.bed (sequential, lengths in col 2)")
    ap.add_argument("--out-bed", required=True, type=Path,
                    help="standard BED with start/end + depth")
    ap.add_argument("--out-stats", required=True, type=Path,
                    help="summary stats TSV")
    ap.add_argument("--label", default="UNKNOWN",
                    help="label for the stats row (e.g., '1000bp')")
    args = ap.parse_args()

    intervals = list(parse_combined(args.combined))
    if not intervals:
        sys.exit("no intervals parsed")

    max_depth = max(r[3] for r in intervals)
    total_bp = sum(end - start for _, start, end, _, _ in intervals)

    args.out_bed.parent.mkdir(parents=True, exist_ok=True)
    with args.out_bed.open("w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["#chrom", "start", "end", "depth", "n_samples",
                    "is_bubble", "frac_present"])
        for chrom, start, end, depth, n in intervals:
            is_bubble = int(depth < max_depth)
            frac = depth / max_depth
            w.writerow([chrom, start, end, depth, n, is_bubble, f"{frac:.4f}"])

    bubbles = [(c, s, e, d) for c, s, e, d, _ in intervals if d < max_depth]
    lengths = [e - s for _, s, e, _ in bubbles]
    fullcov_bp = total_bp - sum(lengths)

    def pct(values, p):
        if not values:
            return 0
        s = sorted(values)
        k = max(0, min(len(s) - 1, int(round(p / 100 * (len(s) - 1)))))
        return s[k]

    def count_ge(values, t):
        return sum(1 for v in values if v >= t)

    row = {
        "label": args.label,
        "total_intervals": len(intervals),
        "max_depth": max_depth,
        "total_bp": total_bp,
        "fullcov_bp": fullcov_bp,
        "fullcov_pct": f"{100 * fullcov_bp / total_bp:.4f}",
        "n_bubbles": len(bubbles),
        "bubble_bp": sum(lengths),
        "bubble_pct": f"{100 * sum(lengths) / total_bp:.4f}",
        "len_median": int(statistics.median(lengths)) if lengths else 0,
        "len_mean": int(statistics.mean(lengths)) if lengths else 0,
        "len_q95": pct(lengths, 95),
        "len_max": max(lengths) if lengths else 0,
        "n_ge_200bp": count_ge(lengths, 200),
        "n_ge_1kb": count_ge(lengths, 1_000),
        "n_ge_10kb": count_ge(lengths, 10_000),
        "n_ge_50kb": count_ge(lengths, 50_000),
        "n_ge_100kb": count_ge(lengths, 100_000),
    }

    args.out_stats.parent.mkdir(parents=True, exist_ok=True)
    write_header = not args.out_stats.exists()
    with args.out_stats.open("a") as fh:
        w = csv.writer(fh, delimiter="\t")
        if write_header:
            w.writerow(row.keys())
        w.writerow(row.values())

    print(f"[{args.label}] {row['n_bubbles']} bubbles, "
          f"{row['bubble_pct']}% of {row['total_bp']/1e6:.2f} Mb, "
          f"median {row['len_median']} bp, max {row['len_max']} bp, "
          f">=50kb: {row['n_ge_50kb']}", file=sys.stderr)


if __name__ == "__main__":
    main()
