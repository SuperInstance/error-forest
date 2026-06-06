//! Gossip protocol inspired by spore dispersal patterns.
//! Provable eventual consistency with ecological convergence guarantees.

use serde::{Serialize, Deserialize};
use std::collections::{HashMap, HashSet};

/// Configuration for the spore gossip protocol.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SporeConfig {
    /// Number of peers to forward to per round (fanout).
    pub fanout: usize,
    /// Time-to-live (max hops before spore dies).
    pub ttl: usize,
    /// Target convergence ratio (0.0–1.0).
    pub convergence_target: f64,
}

impl Default for SporeConfig {
    fn default() -> Self {
        SporeConfig {
            fanout: 3,
            ttl: 10,
            convergence_target: 0.99,
        }
    }
}

/// A spore carrying data through the network.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Spore {
    pub id: u64,
    pub source: usize,
    pub data: Vec<u8>,
    pub ttl: usize,
    pub visited: HashSet<usize>,
    pub generation: usize,
}

/// State of a node in the gossip network.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NodeState {
    pub id: usize,
    pub data: HashMap<u64, Vec<u8>>,
    pub received_generations: HashMap<u64, usize>,
    pub neighbors: Vec<usize>,
}

/// Result of gossip convergence analysis.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ConvergenceResult {
    /// Number of rounds to reach convergence target.
    pub rounds: usize,
    /// Final convergence ratio.
    pub convergence_ratio: f64,
    /// Messages sent total.
    pub messages_sent: usize,
    /// Whether convergence target was reached.
    pub converged: bool,
    /// Per-round convergence ratios.
    pub convergence_history: Vec<f64>,
}

/// The spore gossip network.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SporeNetwork {
    pub config: SporeConfig,
    pub nodes: Vec<NodeState>,
    pub total_spores_originated: u64,
}

impl SporeNetwork {
    /// Create a new spore gossip network.
    pub fn new(num_nodes: usize, config: SporeConfig) -> Self {
        let nodes = (0..num_nodes)
            .map(|id| {
                // Create a connected graph: each node connects to its neighbors
                let mut neighbors = Vec::new();
                for f in 1..=config.fanout {
                    let neighbor = (id + f) % num_nodes;
                    if neighbor != id {
                        neighbors.push(neighbor);
                    }
                    if id >= f {
                        neighbors.push(id - f);
                    }
                }
                neighbors.sort();
                neighbors.dedup();

                NodeState {
                    id,
                    data: HashMap::new(),
                    received_generations: HashMap::new(),
                    neighbors,
                }
            })
            .collect();

        SporeNetwork {
            config,
            nodes,
            total_spores_originated: 0,
        }
    }

    /// Originate a spore from a source node.
    pub fn originate(&mut self, source: usize, data: Vec<u8>) -> Spore {
        let id = self.total_spores_originated;
        self.total_spores_originated += 1;

        // Source immediately stores the data
        self.nodes[source].data.insert(id, data.clone());
        self.nodes[source].received_generations.insert(id, 0);

        Spore {
            id,
            source,
            data,
            ttl: self.config.ttl,
            visited: vec![source].into_iter().collect(),
            generation: 0,
        }
    }

    /// Run one round of gossip dissemination.
    /// Returns the number of messages sent this round.
    pub fn gossip_round(&mut self, pending: &mut Vec<Spore>) -> usize {
        let mut new_pending = Vec::new();
        let mut messages_sent = 0;

        for mut spore in pending.drain(..) {
            if spore.ttl == 0 {
                continue;
            }

            let current_node = *spore.visited.iter().last().unwrap_or(&spore.source);
            if current_node >= self.nodes.len() {
                continue;
            }

            let neighbors = self.nodes[current_node].neighbors.clone();
            let mut rng = SimpleRng2::new(spore.id.wrapping_add(spore.generation as u64));

            // Select fanout neighbors that haven't been visited
            let candidates: Vec<usize> = neighbors.iter()
                .filter(|&n| !spore.visited.contains(n))
                .cloned()
                .collect();

            let selected: Vec<usize> = if candidates.len() <= self.config.fanout {
                candidates
            } else {
                // Deterministic pseudo-random selection
                let mut indices: Vec<usize> = candidates;
                // Fisher-Yates shuffle with deterministic rng
                for i in (1..indices.len()).rev() {
                    let j = (rng.next() as usize) % (i + 1);
                    indices.swap(i, j);
                }
                indices.into_iter().take(self.config.fanout).collect()
            };

            for neighbor in selected {
                messages_sent += 1;

                // Deliver data to neighbor
                let next_gen = spore.generation + 1;
                let prev_gen = self.nodes[neighbor].received_generations.get(&spore.id).copied();

                if prev_gen.unwrap_or(usize::MAX) > next_gen || !self.nodes[neighbor].data.contains_key(&spore.id) {
                    self.nodes[neighbor].data.insert(spore.id, spore.data.clone());
                    self.nodes[neighbor].received_generations.insert(spore.id, next_gen);

                    let mut new_spore = spore.clone();
                    new_spore.visited.insert(neighbor);
                    new_spore.ttl -= 1;
                    new_spore.generation = next_gen;
                    new_pending.push(new_spore);
                }
            }
        }

        *pending = new_pending;
        messages_sent
    }

    /// Run gossip until convergence target is reached or max rounds exceeded.
    pub fn run_until_convergence(&mut self, spores: Vec<Spore>, max_rounds: usize) -> ConvergenceResult {
        let mut pending = spores;
        let mut rounds = 0;
        let mut messages_sent = 0;
        let mut convergence_history = Vec::new();

        for round in 0..max_rounds {
            if pending.is_empty() {
                break;
            }

            let msg = self.gossip_round(&mut pending);
            messages_sent += msg;
            rounds = round + 1;

            let convergence = self.compute_convergence();
            convergence_history.push(convergence);

            if convergence >= self.config.convergence_target {
                return ConvergenceResult {
                    rounds,
                    convergence_ratio: convergence,
                    messages_sent,
                    converged: true,
                    convergence_history,
                };
            }
        }

        let convergence = self.compute_convergence();
        ConvergenceResult {
            rounds,
            convergence_ratio: convergence,
            messages_sent,
            converged: convergence >= self.config.convergence_target,
            convergence_history,
        }
    }

    /// Compute current convergence ratio across all nodes.
    pub fn compute_convergence(&self) -> f64 {
        if self.nodes.is_empty() {
            return 1.0;
        }

        // All data IDs that exist in the network
        let all_ids: HashSet<u64> = self.nodes.iter()
            .flat_map(|n| n.data.keys().copied())
            .collect();

        if all_ids.is_empty() {
            return 1.0;
        }

        let total_possible = self.nodes.len() * all_ids.len();
        if total_possible == 0 {
            return 1.0;
        }

        let total_have: usize = self.nodes.iter()
            .map(|n| all_ids.intersection(&n.data.keys().copied().collect::<HashSet<_>>()).count())
            .sum();

        total_have as f64 / total_possible as f64
    }

    /// Theoretical bound on convergence rounds.
    /// Based on epidemic model: O(log(N) / log(fanout + 1)).
    pub fn theoretical_convergence_bound(&self) -> usize {
        let n = self.nodes.len() as f64;
        let f = (self.config.fanout + 1) as f64;
        (n.ln() / f.ln()).ceil() as usize + self.config.ttl.min(3)
    }
}

/// Simple deterministic RNG for gossip selection.
#[derive(Clone, Debug)]
struct SimpleRng2 {
    state: u64,
}

impl SimpleRng2 {
    fn new(seed: u64) -> Self {
        SimpleRng2 { state: if seed == 0 { 1 } else { seed } }
    }

    fn next(&mut self) -> u64 {
        self.state ^= self.state << 13;
        self.state ^= self.state >> 7;
        self.state ^= self.state << 17;
        self.state
    }
}
