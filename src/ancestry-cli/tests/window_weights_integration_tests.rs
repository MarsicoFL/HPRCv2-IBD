//! End-to-end behavioural tests for per-window support weighting.
//!
//! Property under test: if you take a noiseless K-state HMM run and drop a
//! contiguous block of windows to low weight (analog of a low-coverage
//! pangenome region), Viterbi+forward-backward should ignore the local
//! emission and follow the transition prior across that block. That is the
//! whole point of (B) — making centromeric / low-support windows
//! non-informative without dropping them.

use std::collections::HashMap;

use impopk_ancestry_cli::{
    apply_window_weights, forward_backward_from_log_emissions,
    precompute_log_emissions, viterbi_from_log_emissions, weights_for_observations,
    AncestralPopulation, AncestryHmmParams, AncestryObservation, WeightMode, WindowWeights,
};

fn make_pop(name: &str, haplotypes: &[&str]) -> AncestralPopulation {
    AncestralPopulation {
        name: name.to_string(),
        haplotypes: haplotypes.iter().map(|s| s.to_string()).collect(),
    }
}

/// Build a controlled HMM with K=2 and a single ground-truth ancestry (A) for
/// all windows. The adversarial `noisy_block` window range gets emissions
/// pointing at the wrong state (B). This is the centromere analog: a
/// contiguous patch of misleading identity in an otherwise uniform region.
fn build_hmm_scenario(
    n_windows: usize,
    noisy_block: std::ops::Range<usize>,
) -> (Vec<AncestryObservation>, AncestryHmmParams) {
    let pop_a = make_pop("A", &["A_h1", "A_h2"]);
    let pop_b = make_pop("B", &["B_h1", "B_h2"]);
    let params = AncestryHmmParams::new(vec![pop_a.clone(), pop_b.clone()], 0.001);

    let mut observations = Vec::with_capacity(n_windows);
    let window_size = 1_000u64;

    for t in 0..n_windows {
        let mut sims = HashMap::new();
        let (high, low) = if noisy_block.contains(&t) {
            (&pop_b, &pop_a) // adversarial flip
        } else {
            (&pop_a, &pop_b) // ground truth = A
        };
        for h in &high.haplotypes {
            sims.insert(h.clone(), 0.99);
        }
        for h in &low.haplotypes {
            sims.insert(h.clone(), 0.40);
        }
        observations.push(AncestryObservation {
            chrom: "chr1".to_string(),
            start: (t as u64) * window_size,
            end: ((t as u64) + 1) * window_size,
            sample: "Q1".to_string(),
            similarities: sims,
            coverage_ratios: None,
            haplotype_consistency_bonus: None,
        });
    }

    (observations, params)
}

fn count_state_transitions(states: &[usize]) -> usize {
    states.windows(2).filter(|w| w[0] != w[1]).count()
}

#[test]
fn noisy_block_misleads_unweighted_viterbi() {
    // Sanity: without weighting, a 10-window adversarial block flips the
    // Viterbi path inside the block (truth = A everywhere, baseline picks B
    // for the block). This is the failure mode (B) targets.
    let n = 60;
    let noisy = 25..35usize;
    let (obs, params) = build_hmm_scenario(n, noisy.clone());

    let log_e = precompute_log_emissions(&obs, &params);
    let states = viterbi_from_log_emissions(&log_e, &params);

    let wrong_in_block = (noisy.clone()).filter(|&t| states[t] != 0).count();
    assert!(
        wrong_in_block >= 5,
        "expected adversarial block to flip ≥5 windows without weighting; \
         got {wrong_in_block}. States in block: {:?}",
        &states[noisy]
    );
}

#[test]
fn low_weights_let_viterbi_follow_transition_prior() {
    // The property test: when the noisy block is down-weighted, the HMM
    // ignores the adversarial signal and stays in the surrounding (correct)
    // state across the block.
    let n = 60;
    let noisy = 25..35usize;
    let (obs, params) = build_hmm_scenario(n, noisy.clone());

    let mut log_e = precompute_log_emissions(&obs, &params);
    let mut weights = vec![1.0_f64; n];
    for t in noisy.clone() {
        weights[t] = 0.05;
    }
    apply_window_weights(&mut log_e, &weights, WeightMode::Interp);
    let states = viterbi_from_log_emissions(&log_e, &params);

    // With effective weights, the entire path should remain in state A.
    let wrong_in_block = (noisy.clone()).filter(|&t| states[t] != 0).count();
    assert_eq!(
        wrong_in_block, 0,
        "weighted decoder should stay in state A across the block; \
         got states {:?}",
        &states[noisy.clone()]
    );

    let block_transitions = count_state_transitions(&states[noisy.clone()]);
    assert_eq!(block_transitions, 0, "block should be monotone");

    // Outside the block, flanks should remain correct (A).
    assert!(states[..noisy.start].iter().all(|&s| s == 0));
    assert!(states[noisy.end..].iter().all(|&s| s == 0));
}

#[test]
fn empty_weight_file_is_noop_for_inference() {
    // With an empty WindowWeights, weights are all 1.0 → apply_window_weights
    // is a no-op and inference matches the unweighted baseline exactly.
    let n = 40;
    let (obs, params) = build_hmm_scenario(n, 0..0);
    let log_e_baseline = precompute_log_emissions(&obs, &params);
    let states_baseline = viterbi_from_log_emissions(&log_e_baseline, &params);
    let post_baseline = forward_backward_from_log_emissions(&log_e_baseline, &params);

    let weights = WindowWeights::empty();
    let ws = weights_for_observations(&weights, &obs);
    let mut log_e = log_e_baseline.clone();
    apply_window_weights(&mut log_e, &ws, WeightMode::Interp);
    let states = viterbi_from_log_emissions(&log_e, &params);
    let post = forward_backward_from_log_emissions(&log_e, &params);

    assert_eq!(states, states_baseline);
    for (a, b) in post.iter().zip(post_baseline.iter()) {
        for (x, y) in a.iter().zip(b.iter()) {
            assert!(
                (x - y).abs() < 1e-12,
                "posteriors differ: {x} vs {y}"
            );
        }
    }
}

#[test]
fn weights_lookup_uses_window_coordinates() {
    // Verify that a WindowWeights loaded from a TSV is correctly matched to
    // the (chrom, start, end) of each observation.
    use std::io::Write;
    let n = 30;
    let (obs, _params) = build_hmm_scenario(n, 0..0);

    // Build a weights TSV for windows 5..10 only.
    let mut tmp = tempfile::NamedTempFile::new().unwrap();
    writeln!(tmp, "chrom\tstart\tend\tweight").unwrap();
    for t in 5..10usize {
        let start = (t as u64) * 1_000;
        let end = start + 1_000;
        writeln!(tmp, "chr1\t{start}\t{end}\t0.10").unwrap();
    }
    tmp.flush().unwrap();

    let ww = WindowWeights::load(tmp.path()).unwrap();
    let ws = weights_for_observations(&ww, &obs);

    assert_eq!(ws.len(), n);
    for (t, &w) in ws.iter().enumerate() {
        if (5..10).contains(&t) {
            assert!((w - 0.10).abs() < 1e-12, "window {t} weight {w}");
        } else {
            assert_eq!(w, 1.0, "window {t} should default to 1.0, got {w}");
        }
    }
}

#[test]
fn mult_mode_shrinks_evidence_uniformly() {
    // mult-mode at w=0.5 should halve every log-emission entry in the block;
    // shape is preserved (argmax unchanged) but absolute magnitudes shrink.
    let n = 20;
    let (obs, params) = build_hmm_scenario(n, 0..0);

    let log_e_base = precompute_log_emissions(&obs, &params);

    let block = 5..15usize;
    let mut weights = vec![1.0_f64; n];
    for t in block.clone() {
        weights[t] = 0.5;
    }
    let mut log_e = log_e_base.clone();
    apply_window_weights(&mut log_e, &weights, WeightMode::Mult);

    for t in block {
        for (k, &base) in log_e_base[t].iter().enumerate() {
            if base.is_finite() {
                assert!(
                    (log_e[t][k] - 0.5 * base).abs() < 1e-12,
                    "mult mode mismatch at t={t}, k={k}: {base} -> {}",
                    log_e[t][k]
                );
            }
        }
    }
}

#[test]
fn weights_apply_only_to_listed_windows() {
    // Confirm that windows outside the weights file are not perturbed at all.
    let n = 20;
    let (obs, params) = build_hmm_scenario(n, 0..0);
    let log_e_base = precompute_log_emissions(&obs, &params);

    let mut weights = vec![1.0_f64; n];
    for t in 5..15usize {
        weights[t] = 0.3;
    }
    let mut log_e = log_e_base.clone();
    apply_window_weights(&mut log_e, &weights, WeightMode::Interp);

    // Outside the block: untouched bit-for-bit.
    for t in (0..5).chain(15..n) {
        for (k, &base) in log_e_base[t].iter().enumerate() {
            assert_eq!(
                log_e[t][k], base,
                "window {t} outside weight block was modified"
            );
        }
    }
}
