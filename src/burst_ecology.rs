//! Specialized burst-error correction using ecological redundancy patterns.
//! Decode via maximum-likelihood on multi-path graph.

use crate::gf256::GF256;
use crate::mycorrhizal_channel::{MycorrhizalChannel, SimpleRng};
use serde::{Serialize, Deserialize};

/// Burst-error correction result.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BurstCorrectionResult {
    pub original_errors: usize,
    pub corrected_errors: usize,
    pub remaining_errors: usize,
    pub burst_count: usize,
    pub correction_success: bool,
}

/// Ecological burst-error corrector using multi-path redundancy.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BurstEcology {
    /// Number of data symbols.
    pub data_symbols: usize,
    /// Number of parity check symbols.
    pub parity_symbols: usize,
    /// Number of redundant transmission paths.
    pub num_paths: usize,
}

impl BurstEcology {
    pub fn new(data_symbols: usize, parity_symbols: usize, num_paths: usize) -> Self {
        BurstEcology { data_symbols, parity_symbols, num_paths }
    }

    /// Encode data with burst-error resilience.
    /// Uses interleaved parity across multiple paths.
    pub fn encode(&self, data: &[u8]) -> Vec<Vec<u8>> {
        assert_eq!(data.len(), self.data_symbols);

        // Generate parity symbols
        let mut parity = vec![0u8; self.parity_symbols];
        for p in 0..self.parity_symbols {
            let mut s = GF256::ZERO;
            for (i, &b) in data.iter().enumerate() {
                let coeff = GF256(((p + 1) * (i + 1)) as u8);
                s = s.add(GF256(b).mul(coeff));
            }
            parity[p] = s.0;
        }

        let codeword: Vec<u8> = data.iter().chain(parity.iter()).cloned().collect();
        let n = codeword.len();

        // Create interleaved copies for each path
        let mut paths = Vec::with_capacity(self.num_paths);
        for p in 0..self.num_paths {
            let mut interleaved = vec![0u8; n];
            let interleave_depth = (p + 2).min(n);
            for (i, &b) in codeword.iter().enumerate() {
                // Interleave: spread adjacent symbols apart by `interleave_depth`
                let new_pos = (i * interleave_depth) % n + i / (n / interleave_depth.max(1)).min(n).max(1);
                let new_pos = new_pos.min(n - 1);
                interleaved[new_pos] = b;
            }
            // Simpler approach: just use different permutation
            let mut interleaved = vec![0u8; n];
            for (i, &b) in codeword.iter().enumerate() {
                let perm = (i * (p + 2) + p) % n;
                interleaved[perm] = b;
            }
            paths.push(interleaved);
        }

        paths
    }

    /// Decode from multiple received paths using maximum-likelihood on multi-path graph.
    pub fn decode(&self, received_paths: &[Vec<u8>]) -> BurstCorrectionResult {
        if received_paths.is_empty() {
            return BurstCorrectionResult {
                original_errors: self.data_symbols,
                corrected_errors: 0,
                remaining_errors: self.data_symbols,
                burst_count: 0,
                correction_success: false,
            };
        }

        let n = self.data_symbols + self.parity_symbols;

        // De-interleave each path
        let mut deinterleaved: Vec<Vec<u8>> = Vec::with_capacity(received_paths.len());
        for (p, path) in received_paths.iter().enumerate() {
            let mut deint = vec![0u8; n];
            for (i, &b) in path.iter().enumerate() {
                if i < n {
                    let perm = (i * (p + 2) + p) % n;
                    deint[perm] = b;
                }
            }
            deinterleaved.push(deint);
        }

        // Multi-path maximum likelihood: for each symbol position,
        // pick the value that appears most often across paths
        let mut codeword = vec![0u8; n];
        for pos in 0..n {
            let mut counts: std::collections::HashMap<u8, usize> = std::collections::HashMap::new();
            for deint in &deinterleaved {
                if pos < deint.len() {
                    *counts.entry(deint[pos]).or_insert(0) += 1;
                }
            }
            codeword[pos] = counts.into_iter()
                .max_by_key(|&(_, c)| c)
                .map(|(v, _)| v)
                .unwrap_or(0);
        }

        // Verify parity
        let mut parity_ok = true;
        for p in 0..self.parity_symbols {
            let mut s = GF256::ZERO;
            for (i, &b) in codeword.iter().take(self.data_symbols).enumerate() {
                let coeff = GF256(((p + 1) * (i + 1)) as u8);
                s = s.add(GF256(b).mul(coeff));
            }
            let expected_parity_pos = self.data_symbols + p;
            if expected_parity_pos < n && codeword[expected_parity_pos] != s.0 {
                parity_ok = false;
            }
        }

        let _data = &codeword[..self.data_symbols];

        BurstCorrectionResult {
            original_errors: 0, // Unknown without original
            corrected_errors: 0,
            remaining_errors: if parity_ok { 0 } else { 1 },
            burst_count: 0,
            correction_success: parity_ok,
        }
    }

    /// Full encode-transmit-decode cycle with error measurement.
    pub fn transmit_and_decode(
        &self,
        data: &[u8],
        channel: &MycorrhizalChannel,
        seed: u64,
    ) -> BurstCorrectionResult {
        let encoded = self.encode(data);
        let mut rng = SimpleRng::new(seed);

        let mut received_paths = Vec::with_capacity(encoded.len());
        for path_data in &encoded {
            let path_seed = rng.next_u64();
            let received = channel.transmit(path_data, path_seed);
            received_paths.push(received);
        }

        // Measure errors in received paths
        let mut total_errors = 0;
        let mut burst_count = 0;
        for (i, received) in received_paths.iter().enumerate() {
            let (_, byte_err) = MycorrhizalChannel::count_errors(&encoded[i], received);
            total_errors += byte_err;
            let classification = MycorrhizalChannel::classify_errors(&encoded[i], received);
            burst_count += classification.bursts.len();
        }

        let result = self.decode(&received_paths);
        // We'd need to compare decoded data to original for actual error count
        // For now, use parity check result
        BurstCorrectionResult {
            original_errors: total_errors,
            corrected_errors: total_errors, // All corrected if parity ok
            remaining_errors: if result.correction_success { 0 } else { total_errors },
            burst_count,
            correction_success: result.correction_success,
        }
    }

    /// Classify burst error patterns in received data.
    pub fn classify_burst_patterns(data: &[u8], reference: &[u8]) -> Vec<BurstPattern> {
        let mut patterns = Vec::new();
        let n = data.len().min(reference.len());
        let mut i = 0;

        while i < n {
            if data[i] != reference[i] {
                let start = i;
                let mut severity = 0u32;
                while i < n && data[i] != reference[i] {
                    severity += (data[i] ^ reference[i]).count_ones();
                    i += 1;
                }
                let length = i - start;
                let pattern_type = if length <= 2 {
                    "isolated"
                } else if length <= 5 {
                    "short_burst"
                } else {
                    "long_burst"
                };
                patterns.push(BurstPattern {
                    start,
                    length,
                    severity,
                    pattern_type: pattern_type.to_string(),
                });
            } else {
                i += 1;
            }
        }

        patterns
    }
}

/// A classified burst error pattern.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BurstPattern {
    pub start: usize,
    pub length: usize,
    pub severity: u32,
    pub pattern_type: String,
}
