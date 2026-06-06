//! Reed-Solomon codes over GF(256) for comparison with ecological codes.

use crate::gf256::GF256;
use serde::{Serialize, Deserialize};

/// A Reed-Solomon encoder/decoder over GF(256).
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ReedSolomon {
    /// Number of data symbols.
    pub data_symbols: usize,
    /// Number of parity symbols (2t — can correct up to t symbol errors).
    pub parity_symbols: usize,
}

impl ReedSolomon {
    pub fn new(data_symbols: usize, parity_symbols: usize) -> Self {
        ReedSolomon { data_symbols, parity_symbols }
    }

    /// Maximum correctable symbol errors.
    pub fn t(&self) -> usize {
        self.parity_symbols / 2
    }

    /// Total codeword length.
    pub fn n(&self) -> usize {
        self.data_symbols + self.parity_symbols
    }

    /// Generate the generator polynomial g(x) = prod(x - alpha^i) for i=0..2t-1
    fn generator_poly(&self) -> Vec<GF256> {
        let mut g = vec![GF256::ONE];
        for i in 0..self.parity_symbols {
            let root = GF256(2).pow(i as u32);
            let factor = vec![GF256::ONE, root];
            g = GF256::poly_mul(&g, &factor);
        }
        g
    }

    /// Encode data into a codeword with parity appended.
    pub fn encode(&self, data: &[u8]) -> Vec<u8> {
        assert_eq!(data.len(), self.data_symbols, "data length mismatch");

        let g = self.generator_poly();
        // Pad data polynomial with parity_symbols zeros
        let mut padded: Vec<GF256> = data.iter().map(|&b| GF256(b)).collect();
        padded.extend(vec![GF256::ZERO; self.parity_symbols]);

        // Compute remainder
        let (_, remainder) = GF256::poly_div(&padded, &g);

        // Codeword = data || remainder
        let mut codeword = data.to_vec();
        let rem_start = remainder.len().saturating_sub(self.parity_symbols);
        for i in rem_start..remainder.len() {
            codeword.push(remainder[i].0);
        }
        // Pad if remainder is shorter
        while codeword.len() < self.n() {
            codeword.push(0);
        }
        codeword
    }

    /// Compute syndromes from received codeword.
    fn syndromes(&self, received: &[GF256]) -> Vec<GF256> {
        let mut synd = Vec::with_capacity(self.parity_symbols);
        for i in 0..self.parity_symbols {
            let alpha_i = GF256(2).pow(i as u32);
            let val = GF256::eval_poly(received, alpha_i);
            synd.push(val);
        }
        synd
    }

    /// Decode received codeword. Returns Ok(decoded_data) or Err if uncorrectable.
    pub fn decode(&self, received: &[u8]) -> Result<Vec<u8>, String> {
        if received.len() != self.n() {
            return Err("received length mismatch".into());
        }

        let rx: Vec<GF256> = received.iter().map(|&b| GF256(b)).collect();
        let synd = self.syndromes(&rx);

        // Check if all syndromes are zero
        if synd.iter().all(|&s| s == GF256::ZERO) {
            return Ok(received[..self.data_symbols].to_vec());
        }

        // Berlekamp-Massey to find error locator polynomial
        let error_loc = self.berlekamp_massey(&synd);

        // Chien search for error positions
        let n = self.n();
        let mut error_positions = Vec::new();
        for i in 0..n {
            // Evaluate error locator at alpha^(-i)
            let alpha_neg_i = if i == 0 { GF256::ONE } else { GF256(2).pow((255 - (i % 255)) as u32) };
            let val = GF256::eval_poly(&error_loc, alpha_neg_i);
            if val == GF256::ZERO {
                error_positions.push(i);
            }
        }

        if error_positions.len() != error_loc.len() - 1 {
            return Err("could not find all error positions".into());
        }

        // Forney algorithm to find error magnitudes
        let corrections = self.forney(&synd, &error_loc, &error_positions);

        // Apply corrections
        let mut corrected = rx.clone();
        for (&pos, &correction) in error_positions.iter().zip(corrections.iter()) {
            corrected[pos] = corrected[pos].sub(correction);
        }

        // Verify
        let new_synd = self.syndromes(&corrected);
        if new_synd.iter().all(|&s| s == GF256::ZERO) {
            Ok(corrected[..self.data_symbols].iter().map(|g| g.0).collect())
        } else {
            Err("correction failed — too many errors".into())
        }
    }

    /// Berlekamp-Massey algorithm.
    fn berlekamp_massey(&self, syndromes: &[GF256]) -> Vec<GF256> {
        let n = self.parity_symbols;
        let mut c = vec![GF256::ONE];
        let mut b = vec![GF256::ONE];
        let mut l = 0usize;
        let mut m = 1usize;
        let mut beta = GF256::ONE;

        for n_iter in 0..n {
            let d = (0..=l)
                .map(|i| c[i].mul(syndromes[n_iter - i]))
                .fold(GF256::ZERO, |a, v| a.add(v));

            if d == GF256::ZERO {
                m += 1;
            } else if 2 * l <= n_iter {
                let t = c.clone();
                let coeff = d.mul(beta.inv());
                while c.len() < b.len() + m {
                    c.push(GF256::ZERO);
                }
                for (i, &bi) in b.iter().enumerate() {
                    c[i + m] = c[i + m].sub(coeff.mul(bi));
                }
                l = n_iter + 1 - l;
                b = t;
                beta = d;
                m = 1;
            } else {
                let coeff = d.mul(beta.inv());
                while c.len() < b.len() + m {
                    c.push(GF256::ZERO);
                }
                for (i, &bi) in b.iter().enumerate() {
                    c[i + m] = c[i + m].sub(coeff.mul(bi));
                }
                m += 1;
            }
        }
        c
    }

    /// Forney algorithm for error magnitudes.
    fn forney(&self, syndromes: &[GF256], error_loc: &[GF256], positions: &[usize]) -> Vec<GF256> {
        // Compute error evaluator polynomial: S(x) * sigma(x) mod x^2t
        let omega = {
            let s_poly: Vec<GF256> = syndromes.iter().cloned().collect();
            let product = GF256::poly_mul(&s_poly, error_loc);
            let two_t = self.parity_symbols;
            if product.len() > two_t {
                product[..two_t].to_vec()
            } else {
                product
            }
        };

        // Formal derivative of sigma
        let sigma_prime: Vec<GF256> = error_loc.iter()
            .enumerate()
            .filter_map(|(i, &c)| {
                if i % 2 == 1 { Some(c) } else if i + 1 < error_loc.len() { None } else { None }
            })
            .enumerate()
            .map(|(i, c)| {
                let coeff = GF256((i as u8 * 2 + 1).min(254));
                c.mul(coeff) // simplified: for GF(2^m), formal derivative keeps only odd-indexed terms
            })
            .collect();

        positions.iter().map(|&pos| {
            let alpha_pos = if pos == 0 { GF256::ONE } else { GF256(2).pow(pos as u32) };
            let alpha_inv = alpha_pos.inv();

            let omega_val = GF256::eval_poly(&omega, alpha_inv);
            let sigma_prime_val = GF256::eval_poly(&sigma_prime, alpha_inv);

            if sigma_prime_val == GF256::ZERO {
                GF256::ZERO // Can't determine, return 0
            } else {
                // e_i = - alpha_i * Omega(alpha_i^-1) / sigma'(alpha_i^-1)
                alpha_pos.mul(omega_val).div(sigma_prime_val)
            }
        }).collect()
    }
}
