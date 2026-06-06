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

    /// Encode data using phytochemical parity strategy.
    /// Generates parity symbols via GF(256) Vandermonde-style encoding.
    pub fn encode(&self, data: &[u8]) -> Vec<u8> {
        assert_eq!(data.len(), self.data_symbols);

        let data_gf: Vec<GF256> = data.iter().map(|&b| GF256(b)).collect();
        let mut parity = Vec::with_capacity(self.parity_symbols);

        for p in 0..self.parity_symbols {
            let mut symbol = GF256::ZERO;
            for (i, &d) in data_gf.iter().enumerate() {
                // Vandermonde row: alpha^(p * i) where alpha is the generator
                let coeff = if i == 0 {
                    GF256::ONE
                } else {
                    GF256(2).pow(((p + 1) * i) as u32)
                };
                symbol = symbol.add(d.mul(coeff));
            }
            parity.push(symbol.0);
        }

        let mut codeword = data.to_vec();
        codeword.extend_from_slice(&parity);
        codeword
    }

    /// Encode with multi-path redundancy — generates `redundancy_paths` copies
    /// of the codeword, each with a different linear transformation.
    pub fn encode_multipath(&self, data: &[u8]) -> Vec<Vec<u8>> {
        let base = self.encode(data);
        let mut paths = vec![base.clone()];

        for p in 1..self.redundancy_paths {
            let mut transformed = Vec::with_capacity(base.len());
            for (i, &b) in base.iter().enumerate() {
                // Apply a simple linear transformation per path
                let shift = ((p * 7 + i * 13) % 255) as u8;
                let coeff = GF256(shift);
                transformed.push(GF256(b).mul(coeff).0);
            }
            paths.push(transformed);
        }

        paths
    }

    /// Decode from a single (possibly corrupted) codeword.
    /// Uses syndrome-based error detection and correction.
    pub fn decode(&self, received: &[u8]) -> Result<Vec<u8>, String> {
        if received.len() < self.data_symbols {
            return Err("received data too short".into());
        }

        let rx: Vec<GF256> = received.iter().map(|&b| GF256(b)).collect();

        // Compute syndromes
        let mut syndromes = Vec::with_capacity(self.parity_symbols);
        for p in 0..self.parity_symbols {
            let mut s = GF256::ZERO;
            for (i, &r) in rx.iter().enumerate() {
                let coeff = if i == 0 { GF256::ONE } else { GF256(2).pow(((p + 1) * i) as u32) };
                s = s.add(r.mul(coeff));
            }
            syndromes.push(s);
        }

        // If all syndromes are zero, no errors
        if syndromes.iter().all(|&s| s == GF256::ZERO) {
            return Ok(received[..self.data_symbols].to_vec());
        }

        // Simple error correction: try single-error correction
        // For position j: S_p = e * alpha^((p+1)*j) for all p
        // => S_0 / S_p = alpha^((p+1)*j - j) = alpha^(p*j) ... 
        // Try each position
        let n = rx.len().min(self.codeword_len());
        for j in 0..n {
            if syndromes[0] == GF256::ZERO {
                continue;
            }
            // Expected ratio: S_p / S_0 should equal alpha^(p * j) for single error at position j
            let mut consistent = true;
            let alpha_j = if j == 0 { GF256::ONE } else { GF256(2).pow(j as u32) };
            let error_val = syndromes[0].mul(alpha_j.inv());

            for p in 1..self.parity_symbols {
                let expected = error_val.mul(GF256(2).pow(((p + 1) * j) as u32));
                if expected != syndromes[p] {
                    consistent = false;
                    break;
                }
            }

            if consistent && error_val != GF256::ZERO {
                let mut corrected = received.to_vec();
                corrected[j] = GF256(corrected[j]).sub(error_val).0;
                // Verify correction
                let mut ok = true;
                let corrected_gf: Vec<GF256> = corrected.iter().map(|&b| GF256(b)).collect();
                for p in 0..self.parity_symbols {
                    let mut s = GF256::ZERO;
                    for (i, &r) in corrected_gf.iter().enumerate() {
                        let coeff = if i == 0 { GF256::ONE } else { GF256(2).pow(((p + 1) * i) as u32) };
                        s = s.add(r.mul(coeff));
                    }
                    if s != GF256::ZERO {
                        ok = false;
                        break;
                    }
                }
                if ok {
                    return Ok(corrected[..self.data_symbols].to_vec());
                }
            }
        }

        Err("could not correct errors".into())
    }

    /// Decode from multi-path redundancy — use all received copies to recover data.
    /// Falls back to single-path decode if only one copy is usable.
    pub fn decode_multipath(&self, received_paths: &[Vec<u8>]) -> Result<Vec<u8>, String> {
        if received_paths.is_empty() {
            return Err("no received paths".into());
        }

        // Try decoding each path individually
        for (p, path) in received_paths.iter().enumerate() {
            // Inverse transform for paths > 0
            let untransformed: Vec<u8> = if p == 0 {
                path.clone()
            } else {
                path.iter().enumerate().map(|(i, &b)| {
                    let shift = ((p * 7 + i * 13) % 255) as u8;
                    let coeff = GF256(shift);
                    GF256(b).div(coeff).0
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
                        let shift = ((p * 7 + pos * 13) % 255) as u8;
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
    /// Returns (phyto_errors, repetition_errors) after channel transmission.
    pub fn compare_to_repetition(
        &self,
        data: &[u8],
        channel: &MycorrhizalChannel,
        seed: u64,
    ) -> (usize, usize) {
        let codeword = self.encode(data);

        // PhytoCode: transmit codeword, decode
        let received = channel.transmit(&codeword, seed);
        let phyto_result = self.decode(&received);
        let phyto_errors = match phyto_result {
            Ok(decoded) => data.iter().zip(decoded.iter())
                .map(|(&a, &b)| if a != b { 1 } else { 0 })
                .sum(),
            Err(_) => data.len(),
        };

        // Naive repetition: send data `redundancy_paths + 1` times, majority vote
        let copies = self.redundancy_paths + 1;
        let mut rng = crate::mycorrhizal_channel::SimpleRng::new(seed + 1);
        let mut repetition_votes: Vec<Vec<u8>> = Vec::new();
        for _ in 0..copies {
            let received = channel.transmit(data, rng.next_u64());
            repetition_votes.push(received);
        }

        // Majority vote
        let n = data.len();
        let mut rep_result = vec![0u8; n];
        for pos in 0..n {
            let mut counts: std::collections::HashMap<u8, usize> = std::collections::HashMap::new();
            for vote in &repetition_votes {
                if pos < vote.len() {
                    *counts.entry(vote[pos]).or_insert(0) += 1;
                }
            }
            rep_result[pos] = counts.into_iter()
                .max_by_key(|&(_, c)| c)
                .map(|(v, _)| v)
                .unwrap_or(0);
        }

        let rep_errors: usize = data.iter().zip(rep_result.iter())
            .map(|(&a, &b)| if a != b { 1 } else { 0 })
            .sum();

        (phyto_errors, rep_errors)
    }
}
