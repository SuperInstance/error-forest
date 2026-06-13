//! Error-correcting codes derived from plant signaling strategies.
//! Multi-path redundancy encoding, chemical parity, root-echo retransmission.

use crate::gf256::GF256;
use crate::mycorrhizal_channel::MycorrhizalChannel;
use serde::{Serialize, Deserialize};

/// An error-correcting code inspired by plant phytochemical signaling.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PhytoCode {
    /// Number of data symbols per codeword.
    pub data_symbols: usize,
    /// Number of parity symbols.
    pub parity_symbols: usize,
    /// Number of redundant transmission paths.
    pub redundancy_paths: usize,
}

impl PhytoCode {
    pub fn new(data_symbols: usize, parity_symbols: usize, redundancy_paths: usize) -> Self {
        PhytoCode { data_symbols, parity_symbols, redundancy_paths }
    }

    /// Total codeword length.
    pub fn codeword_len(&self) -> usize {
        self.data_symbols + self.parity_symbols
    }

    /// Maximum correctable errors (parity-based).
    pub fn max_correctable(&self) -> usize {
        self.parity_symbols / 2
    }

    /// Vandermonde coefficient for parity p at data position i.
    fn coeff(p: usize, i: usize) -> GF256 {
        if i == 0 { GF256::ONE } else { GF256(2).pow(((p + 1) * i) as u32) }
    }

    /// Encode data using phytochemical parity strategy.
    /// Generates parity symbols via GF(256) Vandermonde-style encoding.
    pub fn encode(&self, data: &[u8]) -> Vec<u8> {
        assert_eq!(data.len(), self.data_symbols);

        let data_gf: Vec<GF256> = data.iter().map(|&b| GF256(b)).collect();
        let mut parity = Vec::with_capacity(self.parity_symbols);

        for p in 0..self.parity_symbols {
            let mut symbol = GF256::ZERO;
            for (i, &d) in data_gf.iter().enumerate() {
                symbol = symbol.add(d.mul(Self::coeff(p, i)));
            }
            parity.push(symbol.0);
        }

        let mut codeword = data.to_vec();
        codeword.extend_from_slice(&parity);
        codeword
    }

    /// Encode with multi-path redundancy — generates `redundancy_paths + 1` copies
    /// of the codeword (base + redundant), each with a different linear transformation.
    pub fn encode_multipath(&self, data: &[u8]) -> Vec<Vec<u8>> {
        let base = self.encode(data);
        let mut paths = vec![base.clone()];

        for p in 1..=self.redundancy_paths {
            let mut transformed = Vec::with_capacity(base.len());
            for (i, &b) in base.iter().enumerate() {
                // Apply a simple linear transformation per path (skip zero coefficient)
                let shift = ((p * 7 + i * 13) % 254 + 1) as u8; // 1..=254
                let coeff = GF256(shift);
                transformed.push(GF256(b).mul(coeff).0);
            }
            paths.push(transformed);
        }

        paths
    }

    /// Compute syndrome for parity check p: data parity minus stored parity.
    fn compute_syndrome(&self, received: &[u8], p: usize) -> GF256 {
        let mut s = GF256::ZERO;
        for (i, &b) in received.iter().take(self.data_symbols).enumerate() {
            s = s.add(GF256(b).mul(Self::coeff(p, i)));
        }
        // Subtract stored parity
        if self.data_symbols + p < received.len() {
            s = s.sub(GF256(received[self.data_symbols + p]));
        }
        s
    }

    /// Decode from a single (possibly corrupted) codeword.
    /// Uses syndrome-based error detection and correction.
    pub fn decode(&self, received: &[u8]) -> Result<Vec<u8>, String> {
        if received.len() < self.data_symbols {
            return Err("received data too short".into());
        }

        // Compute syndromes: S_p = computed_parity_from_data - stored_parity
        let syndromes: Vec<GF256> = (0..self.parity_symbols)
            .map(|p| self.compute_syndrome(received, p))
            .collect();

        // If all syndromes are zero, no errors
        if syndromes.iter().all(|&s| s == GF256::ZERO) {
            return Ok(received[..self.data_symbols].to_vec());
        }

        // Try single-error correction in data region
        // If error at data position j with value e:
        //   S_p = e * coeff(p, j) for all p
        //   => e = S_0 / coeff(0, j) = S_0 (since coeff(0,j) = alpha^j for p=0)
        // Wait: coeff(p, j) = alpha^((p+1)*j)
        // So e = S_0 / alpha^((0+1)*j) = S_0 / alpha^j

        for j in 0..self.data_symbols {
            let c0 = Self::coeff(0, j);
            if c0 == GF256::ZERO || syndromes[0] == GF256::ZERO {
                continue;
            }
            let error_val = syndromes[0].div(c0);

            let mut consistent = true;
            for p in 1..self.parity_symbols {
                let expected = error_val.mul(Self::coeff(p, j));
                if expected != syndromes[p] {
                    consistent = false;
                    break;
                }
            }

            if consistent {
                let mut corrected = received.to_vec();
                corrected[j] = GF256(corrected[j]).sub(error_val).0;
                // Verify correction
                let all_zero = (0..self.parity_symbols)
                    .all(|p| self.compute_syndrome(&corrected, p) == GF256::ZERO);
                if all_zero {
                    return Ok(corrected[..self.data_symbols].to_vec());
                }
            }
        }

        // Try error in parity region — just return data as-is if parity is wrong
        // Check if data alone produces consistent syndromes with any single parity error
        for j in 0..self.parity_symbols {
            let parity_pos = self.data_symbols + j;
            if parity_pos >= received.len() {
                continue;
            }
            // An error in parity j only affects syndrome j
            let all_others_zero = syndromes.iter().enumerate()
                .all(|(p, &s)| p == j || s == GF256::ZERO);
            if all_others_zero && syndromes[j] != GF256::ZERO {
                // Single parity error, data is fine
                return Ok(received[..self.data_symbols].to_vec());
            }
        }

        Err("could not correct errors".into())
    }

    /// Decode from multi-path redundancy — use all received copies to recover data.
    pub fn decode_multipath(&self, received_paths: &[Vec<u8>]) -> Result<Vec<u8>, String> {
        if received_paths.is_empty() {
            return Err("no received paths".into());
        }

        // Try decoding each path individually
        for (p, path) in received_paths.iter().enumerate() {
            let untransformed: Vec<u8> = if p == 0 {
                path.clone()
            } else {
                path.iter().enumerate().map(|(i, &b)| {
                    let shift = ((p * 7 + i * 13) % 254 + 1) as u8;
                    GF256(b).div(GF256(shift)).0
                }).collect()
            };

            if let Ok(data) = self.decode(&untransformed) {
                return Ok(data);
            }
        }

        // Majority vote fallback
        let n = self.data_symbols;
        let mut result = vec![0u8; n];
        for pos in 0..n {
            let mut counts: std::collections::HashMap<u8, usize> = std::collections::HashMap::new();
            for (p, path) in received_paths.iter().enumerate() {
                if pos < path.len() {
                    let val = if p == 0 {
                        path[pos]
                    } else {
                        let shift = ((p * 7 + pos * 13) % 254 + 1) as u8;
                        GF256(path[pos]).div(GF256(shift)).0
                    };
                    *counts.entry(val).or_insert(0) += 1;
                }
            }
            result[pos] = counts.into_iter()
                .max_by_key(|&(_, c)| c)
                .map(|(v, _)| v)
                .unwrap_or(0);
        }

        Ok(result)
    }

    /// Compare PhytoCode burst error resilience to naive repetition coding.
    pub fn compare_to_repetition(
        &self,
        data: &[u8],
        channel: &MycorrhizalChannel,
        seed: u64,
    ) -> (usize, usize) {
        let codeword = self.encode(data);

        let received = channel.transmit(&codeword, seed);
        let phyto_errors = match self.decode(&received) {
            Ok(decoded) => data.iter().zip(decoded.iter())
                .map(|(&a, &b)| if a != b { 1 } else { 0 })
                .sum(),
            Err(_) => data.len(),
        };

        let copies = self.redundancy_paths + 1;
        let mut rng = crate::mycorrhizal_channel::SimpleRng::new(seed + 1);
        let repetition_votes: Vec<Vec<u8>> = (0..copies)
            .map(|_| channel.transmit(data, rng.next_u64()))
            .collect();

        let n = data.len();
        let rep_result: Vec<u8> = (0..n).map(|pos| {
            let mut counts: std::collections::HashMap<u8, usize> = std::collections::HashMap::new();
            for vote in &repetition_votes {
                if pos < vote.len() {
                    *counts.entry(vote[pos]).or_insert(0) += 1;
                }
            }
            counts.into_iter()
                .max_by_key(|&(_, c)| c)
                .map(|(v, _)| v)
                .unwrap_or(0)
        }).collect();

        let rep_errors: usize = data.iter().zip(rep_result.iter())
            .map(|(&a, &b)| if a != b { 1 } else { 0 })
            .sum();

        (phyto_errors, rep_errors)
    }
}
