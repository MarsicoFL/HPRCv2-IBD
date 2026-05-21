//! Micro-benchmarks for the `window_weights` module.
//!
//! Realistic shape:
//!   - 13 500 windows on chr12 (10 kb grid)
//!   - 3 ancestral populations (K = 3)
//!   - 100 windows down-weighted (centromere-style block)
//!
//! Run with `cargo bench -p impopk-ancestry-cli`.

use std::io::Write;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use impopk_ancestry_cli::{
    apply_window_weights, weights_for_observations, HasWindowKey, WeightMode, WindowWeights,
};

const N_WINDOWS: usize = 13_500;
const WINDOW_SIZE: u64 = 10_000;
const N_STATES: usize = 3;
const N_DOWNWEIGHTED: usize = 100;

#[derive(Clone)]
struct Win {
    chrom: String,
    start: u64,
    end: u64,
}
impl HasWindowKey for Win {
    fn chrom(&self) -> &str {
        &self.chrom
    }
    fn start(&self) -> u64 {
        self.start
    }
    fn end(&self) -> u64 {
        self.end
    }
}

fn synthetic_windows() -> Vec<Win> {
    (0..N_WINDOWS as u64)
        .map(|i| Win {
            chrom: "chr12".to_string(),
            start: i * WINDOW_SIZE,
            end: (i + 1) * WINDOW_SIZE,
        })
        .collect()
}

fn synthetic_weights_file() -> tempfile::NamedTempFile {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    writeln!(f, "chrom\tstart\tend\tweight").unwrap();
    // First 100 windows down-weighted (centromere-style)
    for i in 0..N_DOWNWEIGHTED as u64 {
        let s = i * WINDOW_SIZE;
        writeln!(f, "chr12\t{}\t{}\t0.0940", s, s + WINDOW_SIZE).unwrap();
    }
    // Rest at full trust (uniform 1.0)
    for i in N_DOWNWEIGHTED as u64..N_WINDOWS as u64 {
        let s = i * WINDOW_SIZE;
        writeln!(f, "chr12\t{}\t{}\t1.0000", s, s + WINDOW_SIZE).unwrap();
    }
    f.flush().unwrap();
    f
}

fn bench_load(c: &mut Criterion) {
    let f = synthetic_weights_file();
    c.bench_function("WindowWeights::load 13.5k rows", |b| {
        b.iter(|| black_box(WindowWeights::load(f.path()).unwrap()));
    });
}

fn bench_get(c: &mut Criterion) {
    let f = synthetic_weights_file();
    let w = WindowWeights::load(f.path()).unwrap();
    c.bench_function("WindowWeights::get 13.5k mixed hits", |b| {
        b.iter(|| {
            let mut acc = 0.0_f64;
            for i in 0..N_WINDOWS as u64 {
                acc += w.get("chr12", i * WINDOW_SIZE, (i + 1) * WINDOW_SIZE);
            }
            black_box(acc)
        });
    });
}

fn bench_build_vector(c: &mut Criterion) {
    let f = synthetic_weights_file();
    let w = WindowWeights::load(f.path()).unwrap();
    let obs = synthetic_windows();
    c.bench_function("weights_for_observations 13.5k", |b| {
        b.iter(|| black_box(weights_for_observations(&w, &obs)));
    });
}

fn bench_apply_interp(c: &mut Criterion) {
    // log-emission matrix: 13.5k × 3, all finite
    let make = || vec![vec![-2.0_f64, -1.5, -1.0]; N_WINDOWS];
    let weights: Vec<f64> = (0..N_WINDOWS)
        .map(|i| if i < N_DOWNWEIGHTED { 0.094 } else { 1.0 })
        .collect();
    c.bench_function("apply_window_weights interp 13.5k×3", |b| {
        b.iter_batched(
            make,
            |mut log_e| {
                apply_window_weights(&mut log_e, &weights, WeightMode::Interp);
                black_box(log_e)
            },
            criterion::BatchSize::SmallInput,
        );
    });
}

fn bench_apply_mult(c: &mut Criterion) {
    let make = || vec![vec![-2.0_f64, -1.5, -1.0]; N_WINDOWS];
    let weights: Vec<f64> = (0..N_WINDOWS)
        .map(|i| if i < N_DOWNWEIGHTED { 0.094 } else { 1.0 })
        .collect();
    c.bench_function("apply_window_weights mult 13.5k×3", |b| {
        b.iter_batched(
            make,
            |mut log_e| {
                apply_window_weights(&mut log_e, &weights, WeightMode::Mult);
                black_box(log_e)
            },
            criterion::BatchSize::SmallInput,
        );
    });
}

criterion_group!(
    benches,
    bench_load,
    bench_get,
    bench_build_vector,
    bench_apply_interp,
    bench_apply_mult,
);
criterion_main!(benches);
