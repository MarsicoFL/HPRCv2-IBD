#!/bin/bash
# 2x2 factorial: {uniform 10 kb, bubble-aware capped+weighted} × {deconv ON, deconv OFF}
#
# Goal: isolate the bubble-vs-deconv contribution to chr12 ancestry concordance
# on the n=50 AMR chimera benchmark. The bubble_v2 validation work showed that
# 4 of 5 supposedly-bubble-driven wins (SH2B3/ANKS1B/TBX5/SETD1B) have ZERO
# detected bubbles, so the "win" must come from something else. The auto-config
# default toggles deconvolution; this experiment teases the two apart.
#
# Inputs assumed at $INPUTS_DIR (defaults to ../../../HPRCv2-IBD/research/bubble_v2/00_inputs/):
#   combined.paf           — chimera PAF
#   tracts_n50.tsv         — ground-truth tracts
#   populations.tsv        — AFR/EUR/AMR populations
#   queries.txt            — 50 chimera haplotypes
#   subset.txt             — same as queries
#   chm13_SD.bed           — for validation
#   impg_shuo              — Shuo's impg fork (already built)
#
# Outputs in $OUT_DIR (default ./out):
#   ibs_uniform.tsv, ibs_bubble.tsv
#   ancestry_<pipeline>_<deconv>.tsv  for 4 combinations
#   concordance_2x2.tsv               aggregated per-locus and overall
set -euo pipefail

INPUTS_DIR=${INPUTS_DIR:-../../../../../../HPRCv2-IBD/research/bubble_v2/00_inputs}
OUT_DIR=${OUT_DIR:-out}
ANCESTRY=${ANCESTRY:-../../../target/release/ancestry}
IBS=${IBS:-../../../target/release/ibs}
THREADS=${THREADS:-8}

mkdir -p "$OUT_DIR"

# ---------------------------------------------------------------------------
# Step 1 — Build the two window BEDs (uniform + capped)
# ---------------------------------------------------------------------------
WINDOWS_DIR=$OUT_DIR/windows
mkdir -p "$WINDOWS_DIR"

# Uniform 10kb grid — matches v1 baseline
python3 - <<'EOF' >"$WINDOWS_DIR/uniform.bed"
for s in range(0, 133_275_309, 10_000):
    e = min(s + 10_000, 133_275_309)
    print(f"chr1\t{s}\t{e}".replace("chr1", "chr12"))
EOF

# Capped bubble-aware (50kb cap) — produced by make_capped_windows.py
# Expects bubbles.bed already in INPUTS_DIR (from earlier validation work)
python3 ../../01_centromere_fix/window_capping/make_capped_windows.py \
    --bubbles-bed "$INPUTS_DIR/../../../../bubble_aware/experiments/02_adaptive_windows/chr12_bubbles.bed" \
    --grid-step 10000 --max-bubble-bp 50000 \
    --chrom-rename "CHM13#0#chr12->chr12" \
    --out "$WINDOWS_DIR/bubble_capped.bed" \
    --out-detailed "$WINDOWS_DIR/bubble_capped.detailed.tsv" \
    --out-weights "$WINDOWS_DIR/bubble_capped.weights.tsv"

# ---------------------------------------------------------------------------
# Step 2 — IBS for both window sets (the slow part)
# ---------------------------------------------------------------------------
# Uses impg similarity via the `ibs` wrapper. Each ~minutes-to-hours depending
# on machine. Skipped if outputs exist.
for kind in uniform bubble_capped; do
    OUT_IBS="$OUT_DIR/ibs_${kind}.tsv"
    if [[ -s "$OUT_IBS" ]]; then
        echo "[skip] $OUT_IBS exists"
        continue
    fi
    echo "[ibs] computing identity for $kind windows..."
    # See bubble_aware/README.md §4 for the v1 invocation that produced the
    # paired identity TSV. This experiment intentionally uses `impg similarity`
    # (not the retired ibs-from-paf) so the result is reproducible from the
    # public branch.
    $IBS \
        --alignment "$INPUTS_DIR/combined.paf" \
        --bed "$WINDOWS_DIR/${kind}.bed" \
        --subset-sequence-list "$INPUTS_DIR/subset.txt" \
        --output "$OUT_IBS" \
        -t $THREADS
done

# ---------------------------------------------------------------------------
# Step 3 — Run 4 ancestry decodings (the 2x2 factorial)
# ---------------------------------------------------------------------------
run_ancestry() {
    local label=$1     # uniform_deconvOFF, etc.
    local sim=$2
    local extra=$3     # extra CLI args
    local out="$OUT_DIR/ancestry_${label}.tsv"
    if [[ -s "$out" ]]; then echo "[skip] $out"; return; fi
    echo "[ancestry] $label..."
    $ANCESTRY \
        --similarity-file "$sim" \
        --populations "$INPUTS_DIR/populations.tsv" \
        --query-samples "$INPUTS_DIR/queries.txt" \
        --region chr12 --region-length 133275309 --window-size 10000 \
        --estimate-params \
        $extra \
        --output "$out" \
        -t $THREADS
}

# Uniform × deconv OFF
run_ancestry uniform_deconvOFF \
    "$OUT_DIR/ibs_uniform.tsv" \
    ""

# Uniform × deconv ON
run_ancestry uniform_deconvON \
    "$OUT_DIR/ibs_uniform.tsv" \
    "--auto-configure"

# Bubble-aware capped × deconv OFF, with support weighting (B)
run_ancestry bubble_deconvOFF \
    "$OUT_DIR/ibs_bubble_capped.tsv" \
    "--window-weights $WINDOWS_DIR/bubble_capped.weights.tsv --weight-mode interp"

# Bubble-aware capped × deconv ON, with support weighting (B)
run_ancestry bubble_deconvON \
    "$OUT_DIR/ibs_bubble_capped.tsv" \
    "--window-weights $WINDOWS_DIR/bubble_capped.weights.tsv --weight-mode interp --auto-configure"

# ---------------------------------------------------------------------------
# Step 4 — Per-locus concordance scoring
# ---------------------------------------------------------------------------
python3 score_2x2.py \
    --gt "$INPUTS_DIR/tracts_n50.tsv" \
    --runs \
        uniform_deconvOFF="$OUT_DIR/ancestry_uniform_deconvOFF.tsv" \
        uniform_deconvON="$OUT_DIR/ancestry_uniform_deconvON.tsv" \
        bubble_deconvOFF="$OUT_DIR/ancestry_bubble_deconvOFF.tsv" \
        bubble_deconvON="$OUT_DIR/ancestry_bubble_deconvON.tsv" \
    --out "$OUT_DIR/concordance_2x2.tsv"

column -t "$OUT_DIR/concordance_2x2.tsv"
