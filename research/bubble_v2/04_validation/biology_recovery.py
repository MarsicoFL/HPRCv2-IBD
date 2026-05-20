#!/usr/bin/env python3
"""
Check that biologically interesting loci on chr12 are detected as bubbles.

Loci from bubble_aware/manuscript/technical_note.md (Table 2 + section 8):
- TAS2R bitter-receptor cluster      ~10.95 Mb
- TAS2R core (paper case)            10.95–10.98 Mb
- Centromere alpha-satellite + HSAT  36.4–37.16 Mb (5 atomic intervals)
- KRT type-II keratin cluster        58.09 Mb
- SH2B3/ATXN2 (EUR-selection locus)  109.56–109.90 Mb (340 kb)
- SETD1B / SBNO1 / KDM2B             122.70–124.23 Mb (1.53 Mb)
- 12q24.31 subtelomere SD            131.70 Mb
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


LOCI = [
    ("TAS2R_full",       "chr12", 10_800_000, 11_000_000),
    ("TAS2R_core",       "chr12", 10_952_000, 10_958_000),
    ("Centromere_alpha", "chr12", 36_400_000, 37_200_000),
    ("KRT_typeII",       "chr12", 58_000_000, 58_200_000),
    ("ANKS1B_q23.1",     "chr12", 98_870_000, 101_050_000),
    ("SH2B3_ATXN2",      "chr12", 109_560_000, 109_900_000),
    ("TBX5_TBX3",        "chr12", 113_080_000, 113_540_000),
    ("SETD1B_q24.32",    "chr12", 122_700_000, 124_230_000),
    ("Subtelomere_SD",   "chr12", 131_500_000, 132_000_000),
]


def load_bubbles(path: Path):
    out = []
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            chrom = parts[0].split("#")[-1] if "#" in parts[0] else parts[0]
            start, end = int(parts[1]), int(parts[2])
            if parts[5] == "1":
                out.append((chrom, start, end, int(parts[3])))
    return out


def overlapping(bubbles, locus_chrom, locus_start, locus_end):
    out = []
    for c, s, e, d in bubbles:
        if c != locus_chrom:
            continue
        if e <= locus_start or s >= locus_end:
            continue
        out.append((s, e, d))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bubbles", action="append", required=True,
                    help="format: label=path/to/bubbles.bed (repeatable)")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    sources = []
    for spec in args.bubbles:
        label, path = spec.split("=", 1)
        sources.append((label, load_bubbles(Path(path))))

    rows = [["locus", "chrom", "start", "end", "length_kb"]
            + [s[0] + "_n" for s in sources]
            + [s[0] + "_total_bp" for s in sources]
            + [s[0] + "_min_depth" for s in sources]]

    for name, chrom, s, e in LOCI:
        row = [name, chrom, s, e, f"{(e - s) / 1000:.1f}"]
        n_counts = []
        bp_counts = []
        depths = []
        for label, bubbles in sources:
            hits = overlapping(bubbles, chrom, s, e)
            n_counts.append(len(hits))
            bp_counts.append(sum(he - hs for hs, he, _ in hits))
            depths.append(min((d for _, _, d in hits), default=-1))
        row.extend(n_counts + bp_counts + depths)
        rows.append(row)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerows(rows)

    for row in rows:
        print("\t".join(str(x) for x in row))


if __name__ == "__main__":
    main()
