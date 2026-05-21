#![warn(missing_docs)]
//! Per-window confidence weights for emission scaling.
//!
//! A window weight `w_t ∈ [0, 1]` expresses how much the HMM should trust the
//! similarity observation at window `t`. Weights typically come from a
//! pangenome depth signal: `w = mean_depth / max_depth`. A window where the
//! full panel aligns has weight `1.0`; a window where most of the panel
//! drops out has weight near `0.0`.
//!
//! ## Modes
//!
//! - **`WeightMode::Interp` (default)** — convex blend with a uniform emission:
//!   ```text
//!   log_e'(t, k) = w_t · log_e(t, k) + (1 − w_t) · log(1/K)
//!   ```
//!   When `w_t → 0` the emission becomes uniform and the HMM is forced to
//!   lean on the transition prior — that is, on flanking-window evidence —
//!   instead of on the unreliable local identity signal.
//!
//! - **`WeightMode::Mult`** — multiplicative scaling:
//!   ```text
//!   log_e'(t, k) = w_t · log_e(t, k)
//!   ```
//!   Shrinks all log-emissions toward zero. Because `log_e = 0` is the log of
//!   probability 1, after the forward-backward normalization this also
//!   approaches uniform; the relative ranking between states is preserved
//!   while their discriminative power is reduced. Behaviour is similar to
//!   `Interp` at the same `w` with slightly less pull toward the prior.
//!
//! ## File format
//!
//! A TSV with header `chrom <TAB> start <TAB> end <TAB> weight`. The header
//! is optional: if the first non-comment row lacks any of the four expected
//! column names, columns are taken in positional order and that row is parsed
//! as data.
//!
//! ```text
//! chrom    start    end      weight
//! chr1     0        10000    1.0000
//! chr1     10000    20000    0.5500
//! ```
//!
//! Windows in the similarity input that have no matching row default to
//! weight 1.0 (full trust — identical to omitting the file).

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};

/// How a per-window weight is applied to the log-emission row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WeightMode {
    /// Convex blend with uniform: `w · log_e + (1 − w) · log(1/K)`.
    Interp,
    /// Multiplicative: `w · log_e`. No pull toward uniform.
    Mult,
}

impl Default for WeightMode {
    fn default() -> Self {
        Self::Interp
    }
}

impl std::fmt::Display for WeightMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(match self {
            WeightMode::Interp => "interp",
            WeightMode::Mult => "mult",
        })
    }
}

impl std::str::FromStr for WeightMode {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "interp" => Ok(WeightMode::Interp),
            "mult" => Ok(WeightMode::Mult),
            other => Err(format!(
                "unknown weight mode: '{other}' (expected 'interp' or 'mult')"
            )),
        }
    }
}

/// A keyed table of per-window confidence weights.
///
/// Lookups are by `(chrom, start, end)`. Missing windows return `1.0`.
///
/// Internally the table is a nested map (`chrom → (start, end) → weight`)
/// so that `get(&str, u64, u64)` lookups do not allocate a `String` per call.
#[derive(Debug, Default, Clone)]
pub struct WindowWeights {
    /// `chrom -> (start, end) -> weight`.
    table: HashMap<String, HashMap<(u64, u64), f64>>,
}

impl WindowWeights {
    /// Empty weights table — every lookup returns `1.0`. Useful as an
    /// identity element / no-op.
    pub fn empty() -> Self {
        Self::default()
    }

    /// Load weights from a TSV file with header `chrom\tstart\tend\tweight`.
    ///
    /// Weights are validated to lie in `[0.0, 1.0]`. Any row outside that
    /// range is an error so silent miscalibration cannot slip through.
    pub fn load(path: &Path) -> Result<Self> {
        let fh = File::open(path)
            .with_context(|| format!("opening weights file {}", path.display()))?;
        let reader = BufReader::new(fh);
        let mut header_seen = false;
        let mut col_chrom = 0usize;
        let mut col_start = 1usize;
        let mut col_end = 2usize;
        let mut col_weight = 3usize;

        let mut table: HashMap<String, HashMap<(u64, u64), f64>> = HashMap::new();
        let mut total = 0usize;
        for (lineno, line_res) in reader.lines().enumerate() {
            let line = line_res
                .with_context(|| format!("reading line {} of {}", lineno + 1, path.display()))?;
            if line.is_empty() {
                continue;
            }
            if line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split('\t').collect();

            // Header detection: only the *first* non-comment row is inspected.
            // A header is identified when the field at position `col_start`
            // (default 1, the "start" column under positional layout) does not
            // parse as a u64. This catches both labelled-header files and
            // accidental matches where a chrom name happens to equal one of
            // the keyword strings ("chrom", "start", "end", "weight").
            if !header_seen {
                header_seen = true;
                let position_at_start = parts.get(col_start).copied().unwrap_or("");
                if position_at_start.parse::<u64>().is_err() {
                    // Header row: read column positions by name.
                    for (i, f) in parts.iter().enumerate() {
                        match *f {
                            "chrom" => col_chrom = i,
                            "start" => col_start = i,
                            "end" => col_end = i,
                            "weight" => col_weight = i,
                            _ => {}
                        }
                    }
                    continue;
                }
                // Otherwise this row is data; fall through with default
                // positional columns (chrom, start, end, weight).
            }

            if parts.len() <= col_weight {
                anyhow::bail!(
                    "weights file {}:{} has fewer columns than expected",
                    path.display(),
                    lineno + 1
                );
            }
            let chrom = parts[col_chrom];
            let start: u64 = parts[col_start].parse().with_context(|| {
                format!(
                    "parsing start '{}' at {}:{}",
                    parts[col_start],
                    path.display(),
                    lineno + 1
                )
            })?;
            let end: u64 = parts[col_end].parse().with_context(|| {
                format!(
                    "parsing end '{}' at {}:{}",
                    parts[col_end],
                    path.display(),
                    lineno + 1
                )
            })?;
            let weight: f64 = parts[col_weight].parse().with_context(|| {
                format!(
                    "parsing weight '{}' at {}:{}",
                    parts[col_weight],
                    path.display(),
                    lineno + 1
                )
            })?;

            if !(0.0..=1.0).contains(&weight) || !weight.is_finite() {
                anyhow::bail!(
                    "weight {} at {}:{} is outside [0, 1]",
                    weight,
                    path.display(),
                    lineno + 1
                );
            }
            if end <= start {
                anyhow::bail!(
                    "invalid window {}:{}-{} at {}:{} (end must be strictly greater than start)",
                    chrom,
                    start,
                    end,
                    path.display(),
                    lineno + 1
                );
            }

            table
                .entry(chrom.to_string())
                .or_default()
                .insert((start, end), weight);
            total += 1;
        }

        // Sanity: detect accidental duplicate windows in the file.
        let merged: usize = table.values().map(HashMap::len).sum();
        if merged < total {
            eprintln!(
                "Warning: weights file {} contained {} duplicate window keys (later rows overwrote earlier ones)",
                path.display(),
                total - merged
            );
        }

        Ok(Self { table })
    }

    /// Number of `(chrom, start, end)` entries loaded.
    pub fn len(&self) -> usize {
        self.table.values().map(HashMap::len).sum()
    }

    /// `true` if no weight entries have been loaded.
    pub fn is_empty(&self) -> bool {
        self.table.is_empty() || self.len() == 0
    }

    /// Lookup the weight for a window. Returns `1.0` for unknown windows.
    ///
    /// Does not allocate.
    pub fn get(&self, chrom: &str, start: u64, end: u64) -> f64 {
        self.table
            .get(chrom)
            .and_then(|inner| inner.get(&(start, end)))
            .copied()
            .unwrap_or(1.0)
    }

    /// Build the per-window weight vector for a sequence of observations,
    /// while counting how many observations fell back to the default (1.0)
    /// because their `(chrom, start, end)` was not in the table.
    ///
    /// A high `(missing / total)` ratio combined with a non-empty file
    /// strongly suggests a window-grid mismatch — see
    /// [`Self::warn_if_low_coverage`].
    pub fn build_vector<O: HasWindowKey>(&self, observations: &[O]) -> (Vec<f64>, usize) {
        let mut missing = 0usize;
        let weights = observations
            .iter()
            .map(|o| {
                let w = self.get(o.chrom(), o.start(), o.end());
                if w == 1.0 && !self.has(o.chrom(), o.start(), o.end()) {
                    missing += 1;
                }
                w
            })
            .collect();
        (weights, missing)
    }

    /// `true` iff `(chrom, start, end)` is present in the table (distinct
    /// from "weight equals 1.0 but explicitly listed").
    pub fn has(&self, chrom: &str, start: u64, end: u64) -> bool {
        self.table
            .get(chrom)
            .map(|inner| inner.contains_key(&(start, end)))
            .unwrap_or(false)
    }
}

/// Apply per-window weights to a log-emission matrix in place.
///
/// `log_emissions[t][k]` is the per-window per-state log-emission. `weights[t]`
/// is the confidence weight for window `t`. The matrix is mutated in place.
///
/// In `Interp` mode, each row is convex-blended with the uniform log-emission
/// `-ln(K)`. In `Mult` mode, each row is multiplied by `w_t` (no shift toward
/// uniform).
///
/// `NEG_INFINITY` entries are preserved: they encode "no data for this state
/// in this window" and must not be lifted by the blend.
///
/// # Panics
///
/// Panics if `log_emissions.len() != weights.len()`. Library users that may
/// receive mismatched inputs should call [`try_apply_window_weights`] instead.
pub fn apply_window_weights(
    log_emissions: &mut [Vec<f64>],
    weights: &[f64],
    mode: WeightMode,
) {
    try_apply_window_weights(log_emissions, weights, mode)
        .expect("apply_window_weights: weights.len() must equal log_emissions.len()")
}

/// Same as [`apply_window_weights`] but returns an error on length mismatch
/// instead of panicking. Intended for downstream library users that cannot
/// guarantee matched-length inputs at compile time.
pub fn try_apply_window_weights(
    log_emissions: &mut [Vec<f64>],
    weights: &[f64],
    mode: WeightMode,
) -> Result<()> {
    if log_emissions.is_empty() {
        return Ok(());
    }
    if log_emissions.len() != weights.len() {
        anyhow::bail!(
            "weights length {} does not match number of windows {}",
            weights.len(),
            log_emissions.len()
        );
    }

    let k = log_emissions[0].len();
    if k == 0 {
        return Ok(());
    }
    let log_uniform = -(k as f64).ln();

    for (row, &w) in log_emissions.iter_mut().zip(weights.iter()) {
        if !w.is_finite() {
            continue;
        }
        // Clamp defensively so a bad input does not produce NaN
        let w = w.clamp(0.0, 1.0);
        // No-op fast path
        if w >= 1.0 - f64::EPSILON {
            continue;
        }
        for v in row.iter_mut() {
            if !v.is_finite() {
                // Preserve NEG_INFINITY / NaN entries unchanged
                continue;
            }
            *v = match mode {
                WeightMode::Interp => w * *v + (1.0 - w) * log_uniform,
                WeightMode::Mult => w * *v,
            };
        }
    }
    Ok(())
}

/// Build the per-window weight vector for a sequence of observations.
///
/// For each observation, looks up `(chrom, start, end)` in the weights table.
/// Observations whose window key is missing get weight `1.0`. Use
/// [`WindowWeights::build_vector`] if you also need the count of fallbacks
/// (useful for the high-fallback warning at the CLI layer).
pub fn weights_for_observations<O: HasWindowKey>(
    weights: &WindowWeights,
    observations: &[O],
) -> Vec<f64> {
    observations
        .iter()
        .map(|o| weights.get(o.chrom(), o.start(), o.end()))
        .collect()
}

/// Trait letting `weights_for_observations` work over any observation type
/// that carries a `(chrom, start, end)` window key.
pub trait HasWindowKey {
    /// Chromosome / scaffold name.
    fn chrom(&self) -> &str;
    /// Inclusive 0-based start coordinate of the window.
    fn start(&self) -> u64;
    /// Exclusive end coordinate of the window.
    fn end(&self) -> u64;
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_tsv(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    #[test]
    fn parses_basic_tsv() {
        let f = write_tsv(
            "chrom\tstart\tend\tweight\n\
             chr12\t0\t10000\t1.0\n\
             chr12\t10000\t20000\t0.5\n",
        );
        let w = WindowWeights::load(f.path()).unwrap();
        assert_eq!(w.len(), 2);
        assert_eq!(w.get("chr12", 0, 10_000), 1.0);
        assert_eq!(w.get("chr12", 10_000, 20_000), 0.5);
        // Missing window defaults to 1.0
        assert_eq!(w.get("chr12", 20_000, 30_000), 1.0);
        // Wrong chrom defaults to 1.0
        assert_eq!(w.get("chr1", 0, 10_000), 1.0);
    }

    #[test]
    fn rejects_out_of_range_weight() {
        let f = write_tsv("chrom\tstart\tend\tweight\nchr12\t0\t10000\t1.5\n");
        let err = WindowWeights::load(f.path()).unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("outside [0, 1]"), "got: {msg}");
    }

    #[test]
    fn rejects_nan_weight() {
        let f = write_tsv("chrom\tstart\tend\tweight\nchr12\t0\t10000\tnan\n");
        assert!(WindowWeights::load(f.path()).is_err());
    }

    #[test]
    fn header_optional_when_columns_in_order() {
        let f = write_tsv("chr12\t0\t10000\t0.7\n");
        let w = WindowWeights::load(f.path()).unwrap();
        assert_eq!(w.get("chr12", 0, 10_000), 0.7);
    }

    #[test]
    fn rejects_inverted_window() {
        let f = write_tsv("chrom\tstart\tend\tweight\nchr12\t10000\t5000\t0.5\n");
        let err = WindowWeights::load(f.path()).unwrap_err();
        assert!(format!("{err}").contains("end must be strictly greater"));
    }

    #[test]
    fn rejects_zero_length_window() {
        let f = write_tsv("chrom\tstart\tend\tweight\nchr12\t10000\t10000\t0.5\n");
        let err = WindowWeights::load(f.path()).unwrap_err();
        assert!(format!("{err}").contains("end must be strictly greater"));
    }

    #[test]
    fn chrom_named_like_a_header_keyword_is_data() {
        // Regression: prior header detection triggered on any field
        // matching {"chrom","start","end","weight"}, so a data row with
        // chrom == "chrom" was silently dropped as a header. Now the
        // detector checks whether col_start parses as u64, which is
        // robust to chromosome names that collide with column keywords.
        let f = write_tsv("chrom\t0\t10000\t0.5\n");
        let w = WindowWeights::load(f.path()).unwrap();
        assert_eq!(w.len(), 1);
        assert_eq!(w.get("chrom", 0, 10_000), 0.5);
    }

    #[test]
    fn weight_mode_parses() {
        use std::str::FromStr;
        assert_eq!(WeightMode::from_str("interp").unwrap(), WeightMode::Interp);
        assert_eq!(WeightMode::from_str("Mult").unwrap(), WeightMode::Mult);
        assert!(WeightMode::from_str("garbage").is_err());
    }

    #[test]
    fn weight_mode_default_is_interp() {
        assert_eq!(WeightMode::default(), WeightMode::Interp);
    }

    #[test]
    fn weight_mode_display_roundtrips_through_parse() {
        use std::str::FromStr;
        for m in [WeightMode::Interp, WeightMode::Mult] {
            let s = format!("{m}");
            assert_eq!(WeightMode::from_str(&s).unwrap(), m);
        }
    }

    #[test]
    fn interp_at_weight_one_is_noop() {
        let mut log_e = vec![vec![-2.0_f64, -0.1, -1.5]];
        let original = log_e.clone();
        apply_window_weights(&mut log_e, &[1.0], WeightMode::Interp);
        assert_eq!(log_e, original);
    }

    #[test]
    fn interp_at_weight_zero_yields_uniform() {
        let k = 3;
        let log_uniform = -(k as f64).ln();
        let mut log_e = vec![vec![-2.0_f64, -0.1, -1.5]];
        apply_window_weights(&mut log_e, &[0.0], WeightMode::Interp);
        for v in &log_e[0] {
            assert!(
                (v - log_uniform).abs() < 1e-12,
                "expected {log_uniform}, got {v}"
            );
        }
    }

    #[test]
    fn mult_at_weight_zero_zeroes_row() {
        let mut log_e = vec![vec![-2.0_f64, -0.1, -1.5]];
        apply_window_weights(&mut log_e, &[0.0], WeightMode::Mult);
        for v in &log_e[0] {
            assert!(v.abs() < 1e-12);
        }
    }

    #[test]
    fn neg_infinity_preserved() {
        let mut log_e = vec![vec![f64::NEG_INFINITY, -0.5, -1.0]];
        apply_window_weights(&mut log_e, &[0.5], WeightMode::Interp);
        assert_eq!(log_e[0][0], f64::NEG_INFINITY);
        assert!(log_e[0][1].is_finite());
        assert!(log_e[0][2].is_finite());
    }

    #[test]
    fn empty_inputs_are_safe() {
        // Empty matrix
        let mut log_e: Vec<Vec<f64>> = Vec::new();
        apply_window_weights(&mut log_e, &[], WeightMode::Interp);
        // Matrix with empty rows
        let mut log_e = vec![Vec::<f64>::new()];
        apply_window_weights(&mut log_e, &[0.5], WeightMode::Interp);
    }

    #[test]
    fn weight_clamped_defensively() {
        // 1.5 should clamp to 1.0 → no-op
        let mut log_e = vec![vec![-2.0_f64, -0.1, -1.5]];
        let original = log_e.clone();
        apply_window_weights(&mut log_e, &[1.5], WeightMode::Interp);
        assert_eq!(log_e, original);
        // -0.5 should clamp to 0.0 → fully uniform
        let mut log_e = vec![vec![-2.0_f64, -0.1, -1.5]];
        apply_window_weights(&mut log_e, &[-0.5], WeightMode::Interp);
        let log_uniform = -(3.0_f64).ln();
        for v in &log_e[0] {
            assert!((v - log_uniform).abs() < 1e-12);
        }
    }

    #[test]
    fn interp_midpoint_is_halfway() {
        let k = 2;
        let log_uniform = -(k as f64).ln();
        let mut log_e = vec![vec![-2.0_f64, -0.5]];
        apply_window_weights(&mut log_e, &[0.5], WeightMode::Interp);
        assert!((log_e[0][0] - (0.5 * -2.0 + 0.5 * log_uniform)).abs() < 1e-12);
        assert!((log_e[0][1] - (0.5 * -0.5 + 0.5 * log_uniform)).abs() < 1e-12);
    }

    #[test]
    #[should_panic(expected = "weights.len() must equal log_emissions.len()")]
    fn length_mismatch_panics() {
        let mut log_e = vec![vec![-1.0, -2.0]];
        apply_window_weights(&mut log_e, &[0.5, 0.5], WeightMode::Interp);
    }

    #[test]
    fn try_apply_returns_error_on_length_mismatch() {
        let mut log_e = vec![vec![-1.0_f64, -2.0]];
        let res = try_apply_window_weights(&mut log_e, &[0.5, 0.5], WeightMode::Interp);
        assert!(res.is_err());
        assert!(format!("{}", res.unwrap_err()).contains("does not match"));
    }

    #[test]
    fn has_distinguishes_present_from_default() {
        let f = write_tsv(
            "chrom\tstart\tend\tweight\nchr12\t0\t10000\t1.0\n",
        );
        let w = WindowWeights::load(f.path()).unwrap();
        assert!(w.has("chr12", 0, 10_000));
        assert!(!w.has("chr12", 10_000, 20_000));
        assert_eq!(w.get("chr12", 10_000, 20_000), 1.0); // fallback
    }

    #[test]
    fn lookup_does_not_allocate_string() {
        // Indirect test: querying with the same chrom many times should not
        // grow heap usage. We just verify the API does not require ownership.
        let f = write_tsv(
            "chr12\t0\t10000\t0.7\n",
        );
        let w = WindowWeights::load(f.path()).unwrap();
        let key: &str = "chr12";
        for _ in 0..100 {
            assert_eq!(w.get(key, 0, 10_000), 0.7);
        }
    }
}
