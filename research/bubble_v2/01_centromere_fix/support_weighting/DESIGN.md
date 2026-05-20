# (B) Support-weighted emission — design (not yet implemented)

## What it does

For each window `t`, the emission log-probability for state `k` is
currently computed from the per-population identity statistics (max,
mean, median, or top-K). Under (B), that log-emission is scaled by a
per-window confidence weight `w_t ∈ (0, 1]`:

```
log_e'(t, k) = w_t * log_e(t, k) + (1 - w_t) * log_uniform
```

This is a convex interpolation between the data-driven emission and a
uniform emission. When `w_t = 1` (full-depth window) the behavior is
unchanged. When `w_t = 0.1` (centromere window, depth 22/234) the
emission is 90 % uniform and the HMM is forced to lean on transitions —
i.e. flanking-window interpolation — for the posterior in that window.

This is equivalent in spirit to the v1 ad-hoc fix of `--mask-bed`, but
continuous: a window can be partially trusted instead of binary on/off.

## CLI contract

Add to `ancestry-cli`:

```
--window-weights FILE      TSV with columns: chrom, start, end, weight (0..1)
                           Windows in the input similarity file that have no
                           matching row in this TSV default to weight 1.0.
--weight-mode interp|mult  How to apply the weight to the log-emission.
                              interp (default): convex blend with log(1/K)
                              mult:            log_e' = w * log_e
                                               (sharper down-weighting)
```

The weight TSV is exactly what `make_capped_windows.py` produces. No new
file format.

## Where the change lives

The emission computation in `src/ancestry-cli/src/hmm.rs` builds a
`log_e[t][k]` matrix from the IBS observations for each query. The
support-weighting step inserts there, **after** all the existing
transforms (label smoothing, ZCA whitening, kurtosis weighting, …) and
**before** the forward-backward / Viterbi step.

Sketch:

```rust
fn apply_window_weights(
    log_e: &mut [Vec<f64>],
    weights: &[f64],
    mode: WeightMode,
    k: usize,
) {
    let log_uniform = -(k as f64).ln();
    for (t, row) in log_e.iter_mut().enumerate() {
        let w = weights[t];
        for v in row.iter_mut() {
            *v = match mode {
                WeightMode::Interp => w * *v + (1.0 - w) * log_uniform,
                WeightMode::Mult   => w * *v,
            };
        }
    }
}
```

The `weights` slice is built once per query from the `--window-weights`
file by looking up each window in the input similarity TSV.

## Why interpolation, not multiplication

`mode=mult` shrinks all emissions uniformly when `w` is low. That makes
the HMM less confident overall but does not concentrate posterior on
flanking interpolation — it just flattens the path likelihood. The
correct semantics for "I don't trust this window" is "use the prior /
flanking evidence", which is what interpolating toward uniform achieves.
`mode=mult` is included as a comparison knob, not the default.

## Composability with existing features

(B) sits at the end of the emission pipeline so it composes with every
existing transform (whitening, label smoothing, etc.). The one tricky
case is `--mask-bed`: a masked window already sets all `log_e[t][k]` to
`log_uniform`, so applying (B) on top is a no-op. We will therefore
expect `--mask-bed` and `--window-weights` to be used as alternatives,
not together. The CLI should warn (not fail) if both are passed.

## Validation harness

A unit test that:

1. Builds a 100-window toy HMM with K=2, all log_e from a single
   informative population in the first 50 windows and the other in the
   last 50 — clean step function.
2. Drops a 5-window block in the middle with `w = 0.1`, all other
   weights = 1.
3. Asserts that the Viterbi path in those 5 windows follows the
   transition prior, not the local (noisy) emission.

## What does NOT change

- The IBS computation. `ibs` / `impg similarity` are untouched.
- The HMM topology, transition prior, Baum-Welch loop.
- Default behavior: omitting `--window-weights` is a no-op.

## Estimated work

- Code: ~80 lines in `hmm.rs` + ~30 lines of CLI wiring in `main.rs`.
- Tests: 3–4 unit tests in `tests/`, ~150 lines.
- Integration test: one end-to-end on a tiny synthetic case (~200 windows).

Should fit in a single focused session.
