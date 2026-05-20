#!/usr/bin/env bash
# Example: ancestry with per-window support weighting (v0.2.0+).
#
# Demonstrates the --window-weights / --weight-mode pair on the same chr12
# example used by run.sh. Three runs:
#
#   uniform     — weights = 1.0 everywhere (must equal no-flag baseline)
#   lowcov      — a 200 kb block (15.15–15.17 Mb) gets weight 0.1, simulating
#                 a region with structurally degraded alignment support
#   no_flag     — control, no --window-weights at all
#
# Exit codes:
#   - 0 if all three runs complete and uniform matches no_flag bit-for-bit
#     (proves the flag is a strict no-op when weights are all 1.0)
#   - non-zero otherwise
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output_weights

ANC_BIN=${ANC_BIN:-../../../target/release/ancestry}

"$ANC_BIN" \
    --similarity-file  input/ibs.tsv \
    --window-size      10000 \
    --populations      input/populations.tsv \
    --query-samples    input/queries.txt \
    --emission-model   max \
    --estimate-params \
    --threads          4 \
    --output           output_weights/ancestry_no_flag.tsv

"$ANC_BIN" \
    --similarity-file  input/ibs.tsv \
    --window-size      10000 \
    --populations      input/populations.tsv \
    --query-samples    input/queries.txt \
    --emission-model   max \
    --estimate-params \
    --window-weights   input/weights_uniform.tsv \
    --weight-mode      interp \
    --threads          4 \
    --output           output_weights/ancestry_uniform_weights.tsv

"$ANC_BIN" \
    --similarity-file  input/ibs.tsv \
    --window-size      10000 \
    --populations      input/populations.tsv \
    --query-samples    input/queries.txt \
    --emission-model   max \
    --estimate-params \
    --window-weights   input/weights_lowcov_block.tsv \
    --weight-mode      interp \
    --threads          4 \
    --output           output_weights/ancestry_lowcov_block.tsv

# Property 1: weights-all-1.0 must be bit-identical to the no-flag baseline
# (modulo HashMap iteration order: sort first).
if ! diff \
    <(sort output_weights/ancestry_no_flag.tsv) \
    <(sort output_weights/ancestry_uniform_weights.tsv); then
    echo "FAIL: --window-weights with all-1.0 weights changed the output" >&2
    exit 1
fi

# Property 2: the low-cov block run completes and produces well-formed output.
# We don't assert a specific decoded ancestry — chr12:15.15–15.17 Mb identity
# may already decisively call one population, and this 200 kb is much smaller
# than the chimera benchmark. The smoke check is just that the binary handles
# the file end-to-end.
test -s output_weights/ancestry_lowcov_block.tsv

echo "OK: --window-weights is a strict no-op at w=1.0, and runs cleanly at w<1.0"
