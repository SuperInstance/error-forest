//! Syndrome decoding generalized to hub-and-spoke networks.
//! Mother trees act as parity nodes. Detect which "tree" (node) is compromised.

use serde::{Serialize, Deserialize};
use crate::gf256::GF256;

/// Result of syndrome decoding.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SyndromeResult {
    /// Position of detected error (None if no error detected).
    pub error_position: Option<usize>,
    /// Whether correction was successful.
    pub corrected: bool,
    /// The corrected value at the error position.
    pub corrected_value: Option<u8>,
    /// Confidence in the detection (0.0–1.0).
    pub confidence: f64,
}

/// Hub-and-spoke network with mother tree as parity node.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HubTree {
    /// Number of spoke nodes (data nodes).
    pub spokes: usize,
    /// Hub parity symbols per spoke group.
    pub hub_parity: usize,
    /// Connectivity matrix: which hubs connect to which spokes.
    pub connectivity: Vec<Vec<usize>>,
    /// Current node health status (0.0 = dead, 1.0 = healthy).
    pub node_health: Vec<f64>,
}

impl HubTree {
    /// Create a hub-and-spoke network with given spoke count.
    pub fn new(spokes: usize, hub_parity: usize) -> Self {
        // Create groups of spokes, each group shares a hub (parity node)
        let group_size = hub_parity + 1;
        let num_hubs = (spokes + group_size - 1) / group_size;

        let mut connectivity = Vec::with_capacity(num_hubs);
        for h in 0..num_hubs {
            let mut group = Vec::new();
            for s in 0..=hub_parity {
                let spoke_idx = h * group_size + s;
                if spoke_idx < spokes {
                    group.push(spoke_idx);
                }
            }
            connectivity.push(group);
        }

        let node_health = vec![1.0; spokes + num_hubs];

        HubTree {
            spokes,
            hub_parity,
            connectivity,
            node_health,
        }
    }

    /// Number of hub (parity) nodes.
    pub fn num_hubs(&self) -> usize {
        self.connectivity.len()
    }

    /// Total nodes (spokes + hubs).
    pub fn total_nodes(&self) -> usize {
        self.spokes + self.num_hubs()
    }

    /// Compute parity for a hub given spoke data.
    pub fn compute_hub_parity(&self, hub_idx: usize, spoke_data: &[u8]) -> Vec<u8> {
        let group = &self.connectivity[hub_idx.min(self.connectivity.len() - 1)];
        let mut parity = vec![0u8; self.hub_parity];

        for (p, parity_val) in parity.iter_mut().enumerate() {
            let mut s = GF256::ZERO;
            for &spoke_idx in group {
                if spoke_idx < spoke_data.len() {
                    let coeff = GF256(((p + 1) * (spoke_idx + 1)) as u8);
                    s = s.add(GF256(spoke_data[spoke_idx]).mul(coeff));
                }
            }
            *parity_val = s.0;
        }

        parity
    }

    /// Compute syndromes for all hub groups.
    /// Returns syndrome values per hub. Non-zero means error in that group.
    pub fn compute_syndromes(&self, spoke_data: &[u8], stored_parity: &[Vec<u8>]) -> Vec<Vec<u8>> {
        let mut syndromes = Vec::with_capacity(self.num_hubs());

        for h in 0..self.num_hubs() {
            let current_parity = self.compute_hub_parity(h, spoke_data);
            let stored = stored_parity.get(h).cloned().unwrap_or_default();

            let syndrome: Vec<u8> = current_parity.iter()
                .zip(stored.iter())
                .map(|(&a, &b)| GF256(a).sub(GF256(b)).0)
                .collect();

            syndromes.push(syndrome);
        }

        syndromes
    }

    /// Detect which node is compromised using syndrome analysis.
    pub fn detect_failed_node(&self, spoke_data: &[u8], stored_parity: &[Vec<u8>]) -> SyndromeResult {
        let syndromes = self.compute_syndromes(spoke_data, stored_parity);

        // Find which hub group has non-zero syndromes
        let mut failed_groups: Vec<usize> = Vec::new();
        for (h, syndrome) in syndromes.iter().enumerate() {
            if syndrome.iter().any(|&s| s != 0) {
                failed_groups.push(h);
            }
        }

        if failed_groups.is_empty() {
            return SyndromeResult {
                error_position: None,
                corrected: true,
                corrected_value: None,
                confidence: 1.0,
            };
        }

        // If exactly one group failed, find the specific spoke
        if failed_groups.len() == 1 {
            let group_idx = failed_groups[0];
            let group = &self.connectivity[group_idx];
            let syndrome = &syndromes[group_idx];

            // Try each spoke in the group
            for &spoke_idx in group {
                if spoke_idx >= spoke_data.len() {
                    continue;
                }
                // Check if zeroing out this spoke would fix the syndrome
                // Compute what the parity should be if this spoke is the error
                let mut test_data = spoke_data.to_vec();
                test_data[spoke_idx] = 0;

                // Reconstruct: what value at spoke_idx would zero the syndrome?
                // S_p = e * coeff(p, spoke_idx) for all p
                // e = S_0 / coeff(0, spoke_idx)
                let coeff_0 = GF256(1 * (spoke_idx + 1) as u8);
                if coeff_0 == GF256::ZERO || syndrome.is_empty() {
                    continue;
                }
                let error_val = GF256(syndrome[0]).div(coeff_0);

                // Verify against all parity syndromes
                let mut consistent = true;
                for (p, &s) in syndrome.iter().enumerate() {
                    let coeff = GF256(((p + 1) * (spoke_idx + 1)) as u8);
                    let expected = error_val.mul(coeff);
                    if expected != GF256(s) {
                        consistent = false;
                        break;
                    }
                }

                if consistent {
                    let corrected_val = GF256(spoke_data[spoke_idx]).sub(error_val).0;
                    return SyndromeResult {
                        error_position: Some(spoke_idx),
                        corrected: true,
                        corrected_value: Some(corrected_val),
                        confidence: 1.0 / group.len() as f64,
                    };
                }
            }

            // Couldn't pinpoint the spoke
            return SyndromeResult {
                error_position: Some(self.spokes + group_idx), // Hub itself
                corrected: false,
                corrected_value: None,
                confidence: 0.5,
            };
        }

        // Multiple groups failed — likely a shared node or multiple errors
        SyndromeResult {
            error_position: None,
            corrected: false,
            corrected_value: None,
            confidence: 0.1,
        }
    }

    /// Full encoding: compute and store parity for all hubs.
    pub fn encode(&self, spoke_data: &[u8]) -> Vec<Vec<u8>> {
        (0..self.num_hubs())
            .map(|h| self.compute_hub_parity(h, spoke_data))
            .collect()
    }

    /// Detect and correct a single-node failure.
    pub fn decode(&self, spoke_data: &[u8], stored_parity: &[Vec<u8>]) -> (Vec<u8>, SyndromeResult) {
        let result = self.detect_failed_node(spoke_data, stored_parity);
        let mut corrected_data = spoke_data.to_vec();

        if let (Some(pos), Some(val)) = (result.error_position, result.corrected_value) {
            if pos < corrected_data.len() {
                corrected_data[pos] = val;
            }
        }

        (corrected_data, result)
    }

    /// Simulate node failure by corrupting data at given positions.
    pub fn simulate_failures(&mut self, positions: &[usize]) {
        for &pos in positions {
            if pos < self.node_health.len() {
                self.node_health[pos] = 0.0;
            }
        }
    }
}
