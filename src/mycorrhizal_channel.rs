//! Noisy multi-path biological channels with realistic noise profiles.
//! Models burst errors, asymmetric attenuation, and interference patterns
//! observed in mycorrhizal fungal networks.

use serde::{Serialize, Deserialize};
use std::collections::HashSet;

/// Noise profile describing channel characteristics.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NoiseProfile {
    /// Probability that a burst error occurs at any given symbol.
    pub burst_probability: f64,
    /// Length of burst errors (consecutive corrupted symbols).
    pub burst_length: usize,
    /// Signal attenuation factor (0.0 = total loss, 1.0 = no attenuation).
    pub attenuation: f64,
    /// Per-symbol random error probability (non-burst).
    pub random_error_rate: f64,
}

impl Default for NoiseProfile {
    fn default() -> Self {
        NoiseProfile {
            burst_probability: 0.05,
            burst_length: 4,
            attenuation: 0.9,
            random_error_rate: 0.01,
        }
    }
}

impl NoiseProfile {
    /// Create a low-noise profile (clean forest).
    pub fn low_noise() -> Self {
        NoiseProfile {
            burst_probability: 0.01,
            burst_length: 2,
            attenuation: 0.95,
            random_error_rate: 0.005,
        }
    }

    /// Create a high-noise profile (disturbed ecosystem).
    pub fn high_noise() -> Self {
        NoiseProfile {
            burst_probability: 0.15,
            burst_length: 8,
            attenuation: 0.6,
            random_error_rate: 0.05,
        }
    }

    /// Create a burst-dominant profile.
    pub fn burst_dominant() -> Self {
        NoiseProfile {
            burst_probability: 0.1,
            burst_length: 10,
            attenuation: 0.85,
            random_error_rate: 0.001,
        }
    }
}

/// A node in the mycorrhizal network.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Node {
    pub id: usize,
    pub x: f64,
    pub y: f64,
    pub health: f64, // 0.0–1.0, affects signal quality
}

/// A multi-path biological channel modeling mycorrhizal networks.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MycorrhizalChannel {
    pub nodes: usize,
    pub paths: Vec<Vec<usize>>,
    pub noise_profile: NoiseProfile,
    pub nodes_data: Vec<Node>,
}

impl MycorrhizalChannel {
    /// Create a new channel with given nodes and connectivity.
    pub fn new(nodes: usize, noise_profile: NoiseProfile) -> Self {
        let mut paths = Vec::new();
        let nodes_data: Vec<Node> = (0..nodes)
            .map(|id| Node {
                id,
                x: (id as f64 * 17.3) % 100.0,
                y: (id as f64 * 23.7) % 100.0,
                health: 1.0,
            })
            .collect();

        // Create multi-path connectivity — each node connects to 2-3 neighbors
        for i in 0..nodes {
            let j = (i + 1) % nodes;
            paths.push(vec![i, j]);
            if i + 2 < nodes {
                paths.push(vec![i, i + 2]);
            }
        }

        MycorrhizalChannel {
            nodes,
            paths,
            noise_profile,
            nodes_data,
        }
    }

    /// Add a specific path between nodes.
    pub fn add_path(&mut self, path: Vec<usize>) {
        self.paths.push(path);
    }

    /// Get all paths from source to destination (simple BFS).
    pub fn find_paths(&self, source: usize, dest: usize, max_paths: usize) -> Vec<Vec<usize>> {
        let mut adjacency: Vec<HashSet<usize>> = vec![HashSet::new(); self.nodes];
        for path in &self.paths {
            for window in path.windows(2) {
                adjacency[window[0]].insert(window[1]);
                adjacency[window[1]].insert(window[0]);
            }
        }

        let mut results = Vec::new();
        let mut queue = vec![vec![source]];
        let mut visited_paths = HashSet::new();

        while let Some(current_path) = queue.pop() {
            if results.len() >= max_paths {
                break;
            }
            let last = *current_path.last().unwrap();
            if last == dest {
                let key: Vec<usize> = current_path.clone();
                if visited_paths.insert(key) {
                    results.push(current_path.clone());
                }
                continue;
            }
            if current_path.len() > self.nodes {
                continue;
            }
            for &neighbor in &adjacency[last] {
                if !current_path.contains(&neighbor) {
                    let mut new_path = current_path.clone();
                    new_path.push(neighbor);
                    queue.push(new_path);
                }
            }
        }
        results
    }

    /// Transmit data through the channel, applying noise.
    /// Returns the received data with errors applied.
    pub fn transmit(&self, data: &[u8], seed: u64) -> Vec<u8> {
        let mut result = data.to_vec();
        let mut rng = SimpleRng::new(seed);

        let n = data.len();

        // Apply burst errors
        let mut i = 0;
        while i < n {
            if rng.next_f64() < self.noise_profile.burst_probability {
                let burst_start = i;
                let burst_end = (i + self.noise_profile.burst_length).min(n);
                for j in burst_start..burst_end {
                    result[j] ^= rng.next_u8() & 0xFF;
                }
                i = burst_end;
            } else {
                i += 1;
            }
        }

        // Apply random errors
        for byte in result.iter_mut() {
            if rng.next_f64() < self.noise_profile.random_error_rate {
                *byte ^= rng.next_u8() & 0xFF;
            }
        }

        // Apply attenuation — scale values toward zero
        let att = self.noise_profile.attenuation;
        if att < 1.0 {
            for byte in result.iter_mut() {
                let f = (*byte as f64) * att;
                *byte = f as u8;
            }
        }

        result
    }

    /// Transmit through a specific path, where path health affects error rate.
    pub fn transmit_along_path(&self, data: &[u8], path: &[usize], seed: u64) -> Vec<u8> {
        let path_health: f64 = path.iter()
            .filter(|&&n| n < self.nodes_data.len())
            .map(|&n| self.nodes_data[n].health)
            .product::<f64>()
            .powf(1.0 / path.len().max(1) as f64);

        let mut noisy_profile = self.noise_profile.clone();
        noisy_profile.burst_probability /= path_health.max(0.1);
        noisy_profile.random_error_rate /= path_health.max(0.1);
        noisy_profile.attenuation = self.noise_profile.attenuation * path_health;

        let channel = MycorrhizalChannel {
            nodes: self.nodes,
            paths: self.paths.clone(),
            noise_profile: noisy_profile,
            nodes_data: self.nodes_data.clone(),
        };
        channel.transmit(data, seed)
    }

    /// Multi-path transmission: send data along all available paths, merge by majority vote.
    pub fn transmit_multipath(&self, data: &[u8], source: usize, dest: usize, seed: u64) -> Vec<u8> {
        let paths = self.find_paths(source, dest, 5);
        if paths.is_empty() {
            // Direct transmission if no path found
            return self.transmit(data, seed);
        }

        let mut rng = SimpleRng::new(seed);
        let mut received: Vec<Vec<u8>> = Vec::new();

        for (_i, path) in paths.iter().enumerate() {
            let path_seed = rng.next_u64();
            let received_data = self.transmit_along_path(data, path, path_seed);
            received.push(received_data);
        }

        // Majority vote per byte position
        let n = data.len();
        let mut result = vec![0u8; n];
        for pos in 0..n {
            let mut counts: std::collections::HashMap<u8, usize> = std::collections::HashMap::new();
            for recv in &received {
                if pos < recv.len() {
                    *counts.entry(recv[pos]).or_insert(0) += 1;
                }
            }
            result[pos] = counts.into_iter()
                .max_by_key(|&(_, c)| c)
                .map(|(v, _)| v)
                .unwrap_or(0);
        }
        result
    }

    /// Count bit errors between original and received data.
    pub fn count_errors(original: &[u8], received: &[u8]) -> (usize, usize) {
        let n = original.len().min(received.len());
        let mut bit_errors = 0;
        let mut byte_errors = 0;
        for i in 0..n {
            if original[i] != received[i] {
                byte_errors += 1;
                bit_errors += (original[i] ^ received[i]).count_ones() as usize;
            }
        }
        (bit_errors, byte_errors)
    }

    /// Detect burst error patterns in corrupted data.
    pub fn classify_errors(original: &[u8], received: &[u8]) -> ErrorClassification {
        let n = original.len().min(received.len());
        let mut bursts = Vec::new();
        let mut random_errors = 0;
        let mut i = 0;

        while i < n {
            if original[i] != received[i] {
                let start = i;
                let mut length = 0;
                while i < n && original[i] != received[i] {
                    length += 1;
                    i += 1;
                }
                if length >= 3 {
                    bursts.push(BurstError { start, length });
                } else {
                    random_errors += length;
                }
            } else {
                i += 1;
            }
        }

        let total_errors: usize = bursts.iter().map(|b| b.length).sum::<usize>() + random_errors;
        ErrorClassification {
            bursts,
            random_errors,
            total_errors,
        }
    }
}

/// Classification of error patterns.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ErrorClassification {
    pub bursts: Vec<BurstError>,
    pub random_errors: usize,
    pub total_errors: usize,
}

/// A burst error descriptor.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BurstError {
    pub start: usize,
    pub length: usize,
}

/// Simple deterministic PRNG for reproducible noise generation.
#[derive(Clone, Debug)]
pub struct SimpleRng {
    state: u64,
}

impl SimpleRng {
    pub fn new(seed: u64) -> Self {
        SimpleRng { state: if seed == 0 { 1 } else { seed } }
    }

    pub fn next_u64(&mut self) -> u64 {
        // xorshift64
        self.state ^= self.state << 13;
        self.state ^= self.state >> 7;
        self.state ^= self.state << 17;
        self.state
    }

    pub fn next_u8(&mut self) -> u8 {
        self.next_u64() as u8
    }

    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() & 0x1FFFFF) as f64 / 0x1FFFFF as f64
    }
}
