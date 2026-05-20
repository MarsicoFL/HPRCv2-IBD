#!/usr/bin/env python3
"""
Overlap detected bubbles vs Vollger-style SD catalog (BED format).

For each bubble:
  - is it covered (≥ frac) by any SD interval?
  - which length tier does it fall into?

Reports SD-overlap rate by length tier and the centromere-localized bubbles.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def load_bubbles_bed(path: Path):
    """Read the bubble BED produced by bubble_stats.py (chrom start end depth ...)."""
    out = []
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            depth = int(parts[3])
            is_bubble = parts[5] == "1"
            if is_bubble:
                out.append((chrom, start, end, depth))
    return out


def load_sd_bed(path: Path):
    """Read SD intervals from a BED file, indexed by chromosome."""
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("track") or not line.strip():
                continue
            parts = line.split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            by_chrom[chrom].append((start, end))
    for c in by_chrom:
        by_chrom[c].sort()
    return by_chrom


def overlap_bp(a_start: int, a_end: int, intervals: list[tuple[int, int]]) -> int:
    """Total bp of `a` covered by any interval in sorted `intervals`."""
    total = 0
    for s, e in intervals:
        if e <= a_start:
            continue
        if s >= a_end:
            break
        total += min(e, a_end) - max(s, a_start)
    return total


def normalize_chrom(name: str) -> str:
    """Normalize CHM13#0#chr12 → chr12; pass through otherwise."""
    if "#" in name:
        return name.split("#")[-1]
    return name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bubbles", required=True, type=Path)
    ap.add_argument("--sd-bed", required=True, type=Path)
    ap.add_argument("--label", default="UNKNOWN")
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="bubble counted as SD-supported if SD covers ≥ this fraction")
    ap.add_argument("--out-detail", required=True, type=Path,
                    help="per-bubble TSV with SD overlap")
    ap.add_argument("--out-summary", required=True, type=Path)
    args = ap.parse_args()

    bubbles = load_bubbles_bed(args.bubbles)
    sds = load_sd_bed(args.sd_bed)

    # Normalize bubble chrom names to match SD catalog
    norm_bubbles = [(normalize_chrom(c), s, e, d) for c, s, e, d in bubbles]

    detail_rows = []
    tiers: dict[str, dict[str, int]] = defaultdict(lambda: {"total": 0, "supported": 0})
    centromere_supported = {"total": 0, "supported": 0}

    for chrom, s, e, depth in norm_bubbles:
        length = e - s
        intervals = sds.get(chrom, [])
        cov = overlap_bp(s, e, intervals)
        frac = cov / length if length > 0 else 0.0
        supported = frac >= args.min_frac

        tier = (
            "≥50kb" if length >= 50_000 else
            "10–50kb" if length >= 10_000 else
            "1–10kb" if length >= 1_000 else
            "<1kb"
        )
        tiers[tier]["total"] += 1
        tiers[tier]["supported"] += int(supported)

        in_centromere = (chrom == "chr12") and (s < 38_000_000) and (e > 35_000_000)
        if in_centromere:
            centromere_supported["total"] += 1
            centromere_supported["supported"] += int(supported)

        detail_rows.append([
            args.label, chrom, s, e, length, depth,
            cov, f"{frac:.4f}", int(supported), int(in_centromere), tier,
        ])

    args.out_detail.parent.mkdir(parents=True, exist_ok=True)
    write_header = not args.out_detail.exists()
    with args.out_detail.open("a") as fh:
        w = csv.writer(fh, delimiter="\t")
        if write_header:
            w.writerow(["label", "chrom", "start", "end", "length", "depth",
                        "sd_overlap_bp", "sd_overlap_frac",
                        "sd_supported", "in_centromere", "tier"])
        w.writerows(detail_rows)

    write_header = not args.out_summary.exists()
    with args.out_summary.open("a") as fh:
        w = csv.writer(fh, delimiter="\t")
        if write_header:
            w.writerow(["label", "tier", "n_bubbles", "n_sd_supported",
                        "sd_supported_pct"])
        for tier in ["≥50kb", "10–50kb", "1–10kb", "<1kb"]:
            t = tiers[tier]
            pct = 100 * t["supported"] / t["total"] if t["total"] else 0
            w.writerow([args.label, tier, t["total"], t["supported"], f"{pct:.1f}"])
        c = centromere_supported
        pct = 100 * c["supported"] / c["total"] if c["total"] else 0
        w.writerow([args.label, "centromere(36-38Mb)",
                    c["total"], c["supported"], f"{pct:.1f}"])

    print(f"[{args.label}] tiers:", dict(tiers),
          f"centromere: {centromere_supported}")


if __name__ == "__main__":
    main()
