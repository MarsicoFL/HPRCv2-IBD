#!/usr/bin/env python3
"""
Score the 2x2 deconv experiment.

For each (pipeline, deconv) cell, compute per-locus concordance against
the n=50 ground-truth tracts. Loci of interest are the ones audited in
research/bubble_v2/04_validation/FINDINGS.md:

- TAS2R_core    (bubble-driven win in v1)
- Centromere    (failure mode bubble_v2 targets)
- SH2B3_ATXN2   (claimed win — zero bubbles)
- ANKS1B        (claimed win — zero bubbles)
- TBX5_TBX3     (claimed win — zero bubbles)
- SETD1B        (claimed win — zero bubbles)
- WHOLE_CHR12   (chr12 mean)
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


LOCI = [
    ("TAS2R_full",      "chr12",  10_800_000, 11_000_000),
    ("TAS2R_core",      "chr12",  10_952_000, 10_958_000),
    ("Centromere",      "chr12",  36_400_000, 37_200_000),
    ("KRT_typeII",      "chr12",  58_000_000, 58_200_000),
    ("ANKS1B_q23.1",    "chr12",  98_870_000, 101_050_000),
    ("SH2B3_ATXN2",     "chr12", 109_560_000, 109_900_000),
    ("TBX5_TBX3",       "chr12", 113_080_000, 113_540_000),
    ("SETD1B_q24.32",   "chr12", 122_700_000, 124_230_000),
    ("WHOLE_CHR12",     "chr12",           0, 133_275_309),
]


def parse_run_arg(spec: str) -> tuple[str, Path]:
    label, path = spec.split("=", 1)
    return label, Path(path)


def load_tracts(path: Path):
    """Ground-truth tracts: chim_id, chrom, start, end, ancestry."""
    out = defaultdict(list)
    with path.open() as fh:
        header = next(fh).rstrip().split("\t")
        ci = {h: i for i, h in enumerate(header)}
        for line in fh:
            parts = line.rstrip().split("\t")
            chim = parts[ci["chim_id"]]
            chrom = parts[ci["chrom"]]
            s = int(parts[ci["start"]])
            e = int(parts[ci["end"]])
            anc = parts[ci["ancestry"]]
            out[chim].append((chrom, s, e, anc))
    for c in out:
        out[c].sort(key=lambda r: (r[0], r[1]))
    return dict(out)


def load_decoded(path: Path):
    """Decoded segments from ancestry-cli output."""
    out = defaultdict(list)
    with path.open() as fh:
        header = next(fh).rstrip().split("\t")
        ci = {h: i for i, h in enumerate(header)}
        for line in fh:
            parts = line.rstrip().split("\t")
            if not parts or not parts[0]:
                continue
            chrom = parts[ci["chrom"]]
            s = int(parts[ci["start"]])
            e = int(parts[ci["end"]])
            sample = parts[ci["sample"]]
            anc = parts[ci["ancestry"]]
            out[sample].append((chrom, s, e, anc))
    return dict(out)


def query_overlap(segments, chrom, qs, qe):
    """Yield (overlap_start, overlap_end, ancestry) for each segment that
    intersects [qs, qe) on chrom."""
    for c, s, e, anc in segments:
        if c != chrom:
            continue
        if e <= qs or s >= qe:
            continue
        yield max(s, qs), min(e, qe), anc


def chim_to_sample(chim_id: str) -> str:
    # ground truth uses chim_id like "0", decoded uses sample like "CHIM_00#1"
    return f"CHIM_{int(chim_id):02d}#1"


def score_locus(gt_tracts, decoded, chrom, lstart, lend):
    """Total bp of agreement and total covered bp across all chimeras."""
    agree_bp = 0
    total_bp = 0
    for chim, gt in gt_tracts.items():
        sample = chim_to_sample(chim)
        if sample not in decoded:
            continue
        dec = decoded[sample]

        # Per-base voting: at each base in [lstart, lend), check truth and call.
        # Implementation: walk merged interval cuts to avoid O(L) loop.
        cuts = {lstart, lend}
        for c, s, e, _ in gt:
            if c == chrom and s < lend and e > lstart:
                cuts.add(max(s, lstart)); cuts.add(min(e, lend))
        for c, s, e, _ in dec:
            if c == chrom and s < lend and e > lstart:
                cuts.add(max(s, lstart)); cuts.add(min(e, lend))
        sorted_cuts = sorted(c for c in cuts if lstart <= c <= lend)

        for a, b in zip(sorted_cuts, sorted_cuts[1:]):
            gt_anc = None
            dec_anc = None
            for c, s, e, anc in gt:
                if c == chrom and s <= a and e >= b:
                    gt_anc = anc; break
            for c, s, e, anc in dec:
                if c == chrom and s <= a and e >= b:
                    dec_anc = anc; break
            if gt_anc is None or dec_anc is None:
                continue
            total_bp += (b - a)
            if gt_anc == dec_anc:
                agree_bp += (b - a)
    return agree_bp, total_bp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gt", required=True, type=Path)
    ap.add_argument("--runs", nargs="+", required=True,
                    help="label=path repeated for each pipeline cell")
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    gt = load_tracts(args.gt)
    runs = [parse_run_arg(spec) for spec in args.runs]
    decoded_by_label = {label: load_decoded(p) for label, p in runs}

    rows = [["locus", "chrom", "start", "end"] + [r[0] for r in runs]]
    for name, chrom, s, e in LOCI:
        cells = []
        for label, _ in runs:
            agree, total = score_locus(gt, decoded_by_label[label], chrom, s, e)
            pct = 100 * agree / total if total else float("nan")
            cells.append(f"{pct:.2f}")
        rows.append([name, chrom, str(s), str(e), *cells])

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerows(rows)

    for r in rows:
        print("\t".join(r))


if __name__ == "__main__":
    main()
