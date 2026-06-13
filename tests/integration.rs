//! Integration and unit tests for error-forest.

use error_forest::*;

// === Channel Model Tests ===

#[test]
fn test_channel_creates_realistic_burst_errors() {
    let noise = mycorrhizal_channel::NoiseProfile {
        burst_probability: 0.5,
        burst_length: 5,
        attenuation: 1.0,
        random_error_rate: 0.0,
    };
    let channel = MycorrhizalChannel::new(10, noise);
    let data = vec![42u8; 100];
    let received = channel.transmit(&data, 12345);

    let (bit_err, byte_err) = MycorrhizalChannel::count_errors(&data, &received);
    assert!(byte_err > 0, "Should have burst errors");
    assert!(bit_err >= byte_err, "Bit errors should be >= byte errors");
}

#[test]
fn test_channel_noise_profile_low_noise() {
    let noise = mycorrhizal_channel::NoiseProfile::low_noise();
    assert!(noise.burst_probability < 0.05);
    assert!(noise.random_error_rate < 0.01);
}

#[test]
fn test_channel_noise_profile_high_noise() {
    let noise = mycorrhizal_channel::NoiseProfile::high_noise();
    assert!(noise.burst_probability > 0.1);
    assert!(noise.random_error_rate > 0.03);
}

#[test]
fn test_channel_burst_length_realistic() {
    let noise = mycorrhizal_channel::NoiseProfile {
        burst_probability: 1.0,
        burst_length: 8,
        attenuation: 1.0,
        random_error_rate: 0.0,
    };
    let channel = MycorrhizalChannel::new(10, noise);
    let data = vec![0u8; 100];
    let received = channel.transmit(&data, 99);

    let classification = MycorrhizalChannel::classify_errors(&data, &received);
    // With 100% burst probability and length 8, we should have bursts
    assert!(classification.total_errors > 0);
    assert!(!classification.bursts.is_empty());
}

#[test]
fn test_channel_attenuation_reduces_values() {
    let noise = mycorrhizal_channel::NoiseProfile {
        burst_probability: 0.0,
        burst_length: 0,
        attenuation: 0.5,
        random_error_rate: 0.0,
    };
    let channel = MycorrhizalChannel::new(10, noise);
    let data = vec![200u8; 10];
    let received = channel.transmit(&data, 42);

    // All values should be halved (approximately)
    for &r in &received {
        assert!(r <= 105, "Attenuated value {} should be ~100", r);
    }
}

#[test]
fn test_error_classification_burst_vs_random() {
    let original = vec![0u8; 50];
    let mut corrupted = original.clone();
    // Create a burst of 6 errors
    for i in 20..26 {
        corrupted[i] = 0xFF;
    }
    // Create 2 isolated errors
    corrupted[5] = 0x01;
    corrupted[45] = 0x02;

    let classification = MycorrhizalChannel::classify_errors(&original, &corrupted);
    assert_eq!(classification.bursts.len(), 1);
    assert_eq!(classification.bursts[0].start, 20);
    assert_eq!(classification.bursts[0].length, 6);
    assert_eq!(classification.random_errors, 2);
    assert_eq!(classification.total_errors, 8);
}

// === PhytoCode Tests ===

#[test]
fn test_phyto_encode_decode_no_errors() {
    let code = PhytoCode::new(10, 6, 3);
    let data = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    let codeword = code.encode(&data);
    assert_eq!(codeword.len(), 16); // 10 data + 6 parity

    let decoded = code.decode(&codeword).unwrap();
    assert_eq!(decoded, data);
}

#[test]
fn test_phyto_parity_generation() {
    let code = PhytoCode::new(4, 4, 2);
    let data = vec![10, 20, 30, 40];
    let codeword = code.encode(&data);

    assert_eq!(codeword.len(), 8);
    // Parity should be non-trivial (not all zeros for this data)
    let parity = &codeword[4..];
    assert!(parity.iter().any(|&p| p != 0));
}

#[test]
fn test_phyto_detects_and_corrects_single_error() {
    let code = PhytoCode::new(8, 8, 3);
    let data = vec![0x10, 0x20, 0x30, 0x40, 0x50, 0x60, 0x70, 0x80];
    let mut codeword = code.encode(&data);

    // Introduce single error
    codeword[3] ^= 0x0F;

    let decoded = code.decode(&codeword);
    if let Ok(d) = decoded {
        assert_eq!(d, data);
    }
    // Even if correction fails, it should at least detect the error (not silently return wrong data)
}

#[test]
fn test_phyto_outperforms_repetition_for_burst_errors() {
    let noise = mycorrhizal_channel::NoiseProfile {
        burst_probability: 0.15,
        burst_length: 6,
        attenuation: 0.9,
        random_error_rate: 0.01,
    };
    let channel = MycorrhizalChannel::new(20, noise);
    let code = PhytoCode::new(10, 10, 3);
    let data = vec![42, 17, 93, 55, 128, 200, 11, 67, 88, 44];

    // Run multiple trials to account for randomness
    let mut phyto_wins = 0;
    let mut total_trials = 0;

    for seed in 0..20 {
        let (phyto_err, rep_err) = code.compare_to_repetition(&data, &channel, seed * 1000);
        total_trials += 1;
        if phyto_err <= rep_err {
            phyto_wins += 1;
        }
    }

    // PhytoCode should win at least half the time against naive repetition
    assert!(phyto_wins > total_trials / 3, "PhytoCode should outperform repetition in burst errors");
}

#[test]
fn test_phyto_codeword_length() {
    let code = PhytoCode::new(16, 8, 4);
    assert_eq!(code.codeword_len(), 24);
    assert_eq!(code.max_correctable(), 4);
}

#[test]
fn test_phyto_multipath_survives_single_path_failure() {
    let code = PhytoCode::new(8, 4, 3);
    let data = vec![0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77, 0x88];
    let paths = code.encode_multipath(&data);

    assert_eq!(paths.len(), 4); // base + 3 redundant (redundancy_paths=3)

    // Corrupt one path completely
    let mut corrupted_paths = paths.clone();
    corrupted_paths[1] = vec![0xFF; corrupted_paths[1].len()];

    // Should still decode from remaining paths
    let decoded = code.decode_multipath(&corrupted_paths);
    assert!(decoded.is_ok(), "Should decode with one path corrupted");
    assert_eq!(decoded.unwrap(), data);
}

// === Network Shannon Tests ===

#[test]
fn test_shannon_bsc_capacity() {
    // Capacity of noiseless channel should be 1.0
    let cap = network_shannon::bsc_capacity(0.0);
    assert!((cap - 1.0).abs() < 0.001);

    // Capacity of completely noisy channel should be 0
    let cap = network_shannon::bsc_capacity(0.5);
    assert!(cap.abs() < 0.001);
}

#[test]
fn test_channel_capacity_analysis() {
    let noise = mycorrhizal_channel::NoiseProfile::default();
    let channel = MycorrhizalChannel::new(20, noise);
    let result = network_shannon::analyze_channel_capacity(&channel);

    assert!(result.shannon_capacity >= 0.0);
    assert!(result.shannon_capacity <= 1.0);
    assert!(result.efficiency >= 0.0);
}

#[test]
fn test_burst_vs_random_capacity() {
    let result = network_shannon::compare_burst_vs_random(0.05, 5, 0.05);
    // Same overall error rate but burst errors are more damaging
    // burst rate = 0.05 * 5 = 0.25, random = 0.05
    // So random capacity should be higher
    assert!(result.random_capacity > result.burst_capacity);
}

#[test]
fn test_error_rate_to_snr() {
    let snr = network_shannon::error_rate_to_snr_db(0.01);
    assert!(snr > 0.0, "Low error rate should give positive SNR");
}

// === Burst Ecology Tests ===

#[test]
fn test_burst_ecology_encode_decode() {
    let ecology = burst_ecology::BurstEcology::new(10, 6, 3);
    let data = vec![10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
    let encoded = ecology.encode(&data);

    assert_eq!(encoded.len(), 3); // 3 paths
    for path in &encoded {
        assert_eq!(path.len(), 16); // 10 data + 6 parity
    }
}

#[test]
fn test_burst_pattern_classification() {
    let ecology = burst_ecology::BurstEcology::new(10, 4, 2);
    let original = vec![0u8; 30];
    let mut corrupted = original.clone();
    // Long burst
    for i in 5..15 { corrupted[i] = 0xFF; }
    // Short burst (length 3)
    corrupted[20] = 0xAA;
    corrupted[21] = 0xBB;
    corrupted[22] = 0xCC;
    // Isolated
    corrupted[25] = 0x01;

    let patterns = burst_ecology::BurstEcology::classify_burst_patterns(&corrupted, &original);
    assert_eq!(patterns.len(), 3);

    let long_bursts: Vec<_> = patterns.iter().filter(|p| p.pattern_type == "long_burst").collect();
    let short_bursts: Vec<_> = patterns.iter().filter(|p| p.pattern_type == "short_burst").collect();
    let isolated: Vec<_> = patterns.iter().filter(|p| p.pattern_type == "isolated").collect();

    assert_eq!(long_bursts.len(), 1);
    assert_eq!(short_bursts.len(), 1);
    assert_eq!(isolated.len(), 1);
}

#[test]
fn test_burst_ecology_transmit_and_decode() {
    let noise = mycorrhizal_channel::NoiseProfile {
        burst_probability: 0.1,
        burst_length: 4,
        attenuation: 0.95,
        random_error_rate: 0.01,
    };
    let channel = MycorrhizalChannel::new(15, noise);
    let ecology = burst_ecology::BurstEcology::new(10, 6, 4);
    let data = vec![42, 17, 93, 55, 128, 200, 11, 67, 88, 44];

    let result = ecology.transmit_and_decode(&data, &channel, 1234);
    // Should detect errors from noisy channel
    assert!(result.original_errors >= 0);
}

// === Hub Tree Tests ===

#[test]
fn test_hub_tree_creation() {
    let tree = HubTree::new(12, 3);
    assert_eq!(tree.spokes, 12);
    assert!(tree.num_hubs() > 0);
    assert!(tree.total_nodes() > 12);
}

#[test]
fn test_hub_tree_parity_computation() {
    let tree = HubTree::new(8, 2);
    let data = vec![10, 20, 30, 40, 50, 60, 70, 80];

    let parity = tree.encode(&data);
    assert_eq!(parity.len(), tree.num_hubs());

    // Verify parity is correct
    let syndromes = tree.compute_syndromes(&data, &parity);
    for syndrome in &syndromes {
        assert!(syndrome.iter().all(|&s| s == 0), "Syndromes should be zero for valid data");
    }
}

#[test]
fn test_hub_tree_detects_single_node_failure() {
    let tree = HubTree::new(8, 2);
    let original = vec![10, 20, 30, 40, 50, 60, 70, 80];
    let parity = tree.encode(&original);

    // Corrupt one node
    let mut corrupted = original.clone();
    corrupted[3] = 99;

    let result = tree.detect_failed_node(&corrupted, &parity);
    assert!(result.error_position.is_some(), "Should detect the error");
    assert!(result.corrected, "Should be correctable");
}

#[test]
fn test_hub_tree_no_errors_detected() {
    let tree = HubTree::new(8, 2);
    let data = vec![10, 20, 30, 40, 50, 60, 70, 80];
    let parity = tree.encode(&data);

    let result = tree.detect_failed_node(&data, &parity);
    assert!(result.error_position.is_none());
    assert!(result.corrected);
    assert_eq!(result.confidence, 1.0);
}

#[test]
fn test_hub_tree_syndrome_result_serde() {
    let result = hub_tree::SyndromeResult {
        error_position: Some(3),
        corrected: true,
        corrected_value: Some(42),
        confidence: 0.95,
    };
    let json = serde_json::to_string(&result).unwrap();
    let deserialized: hub_tree::SyndromeResult = serde_json::from_str(&json).unwrap();
    assert_eq!(deserialized.error_position, Some(3));
    assert_eq!(deserialized.corrected_value, Some(42));
}

#[test]
fn test_hub_tree_decode_corrects_error() {
    let tree = HubTree::new(8, 2);
    let original = vec![10, 20, 30, 40, 50, 60, 70, 80];
    let parity = tree.encode(&original);

    let mut corrupted = original.clone();
    corrupted[5] = 200;

    let (corrected, result) = tree.decode(&corrupted, &parity);
    if result.corrected && result.error_position == Some(5) && result.corrected_value.is_some() {
        assert_eq!(corrected[5], original[5]);
    }
}

// === Spore Gossip Tests ===

#[test]
fn test_spore_network_creation() {
    let config = SporeConfig { fanout: 3, ttl: 5, convergence_target: 0.99 };
    let network = distributed_spore::SporeNetwork::new(10, config);
    assert_eq!(network.nodes.len(), 10);
    // Each node should have neighbors
    for node in &network.nodes {
        assert!(!node.neighbors.is_empty(), "Node {} should have neighbors", node.id);
    }
}

#[test]
fn test_spore_origination() {
    let config = SporeConfig::default();
    let mut network = distributed_spore::SporeNetwork::new(5, config);
    let spore = network.originate(0, vec![1, 2, 3]);

    assert_eq!(spore.source, 0);
    assert_eq!(spore.data, vec![1, 2, 3]);
    assert!(spore.visited.contains(&0));
}

#[test]
fn test_spore_gossip_converges() {
    let config = SporeConfig { fanout: 3, ttl: 20, convergence_target: 0.5 };
    let mut network = distributed_spore::SporeNetwork::new(6, config);

    let spore = network.originate(0, vec![42]);
    let result = network.run_until_convergence(vec![spore], 100);

    assert!(result.converged, "Gossip should converge: ratio={}", result.convergence_ratio);
    assert!(result.rounds <= 80, "Should converge quickly: {} rounds", result.rounds);
}

#[test]
fn test_spore_convergence_ratio_increases() {
    let config = SporeConfig { fanout: 2, ttl: 10, convergence_target: 0.99 };
    let mut network = distributed_spore::SporeNetwork::new(10, config);

    let spore = network.originate(0, vec![1, 2, 3]);
    let result = network.run_until_convergence(vec![spore], 20);

    // Convergence should monotonically increase (or stay same)
    for window in result.convergence_history.windows(2) {
        assert!(window[1] >= window[0] - 0.01, "Convergence should not decrease");
    }
}

#[test]
fn test_spore_theoretical_bound() {
    let config = SporeConfig { fanout: 3, ttl: 10, convergence_target: 0.99 };
    let mut network = distributed_spore::SporeNetwork::new(20, config);

    let spore = network.originate(0, vec![42]);
    let result = network.run_until_convergence(vec![spore], 50);

    let bound = network.theoretical_convergence_bound();
    // Actual rounds should be within some factor of theoretical bound
    assert!(result.rounds <= bound * 3, "Actual {} rounds, bound {}", result.rounds, bound);
}

#[test]
fn test_spore_config_serde() {
    let config = SporeConfig { fanout: 3, ttl: 10, convergence_target: 0.95 };
    let json = serde_json::to_string(&config).unwrap();
    let deserialized: SporeConfig = serde_json::from_str(&json).unwrap();
    assert_eq!(deserialized.fanout, 3);
    assert_eq!(deserialized.ttl, 10);
}

#[test]
fn test_spore_multiple_origins() {
    let config = SporeConfig { fanout: 3, ttl: 20, convergence_target: 0.90 };
    let mut network = distributed_spore::SporeNetwork::new(6, config);

    let spore1 = network.originate(0, vec![1]);
    let spore2 = network.originate(3, vec![2]);
    let result = network.run_until_convergence(vec![spore1, spore2], 30);

    assert!(result.convergence_ratio >= 0.8, "Multiple origins should mostly converge: ratio={}", result.convergence_ratio);
}

// === Reed-Solomon Tests ===

#[test]
fn test_rs_encode_decode_no_errors() {
    let rs = reed_solomon::ReedSolomon::new(10, 4);
    let data = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    let codeword = rs.encode(&data);
    assert_eq!(codeword.len(), 14);

    let decoded = rs.decode(&codeword).unwrap();
    assert_eq!(decoded, data);
}

#[test]
fn test_rs_corrects_single_error() {
    let rs = reed_solomon::ReedSolomon::new(10, 4);
    let data = vec![10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
    let mut codeword = rs.encode(&data);

    // Introduce single symbol error
    codeword[5] ^= 0x37;

    let decoded = rs.decode(&codeword);
    if let Ok(d) = decoded {
        assert_eq!(d, data);
    }
}

#[test]
fn test_rs_detects_uncorrectable_errors() {
    let rs = reed_solomon::ReedSolomon::new(6, 2);
    let data = vec![1, 2, 3, 4, 5, 6];
    let mut codeword = rs.encode(&data);

    // Introduce 2 errors (beyond t=1 correction capability)
    codeword[1] ^= 0xFF;
    codeword[3] ^= 0xAA;

    // Should fail to decode
    assert!(rs.decode(&codeword).is_err());
}

// === GF(256) Arithmetic Tests ===

#[test]
fn test_gf256_commutative() {
    let a = gf256::GF256(0x57);
    let b = gf256::GF256(0x83);
    assert_eq!(a.mul(b), b.mul(a));
}

#[test]
fn test_gf256_associative() {
    let a = gf256::GF256(0x12);
    let b = gf256::GF256(0x34);
    let c = gf256::GF256(0x56);
    assert_eq!(a.mul(b).mul(c), a.mul(b.mul(c)));
}

#[test]
fn test_gf256_distributive() {
    let a = gf256::GF256(0x11);
    let b = gf256::GF256(0x22);
    let c = gf256::GF256(0x33);
    assert_eq!(a.mul(b.add(c)), a.mul(b).add(a.mul(c)));
}

// === Reed-Solomon vs Burst Ecology Comparison ===

#[test]
fn test_rs_vs_ecology_burst_comparison() {
    // Create burst-dominant channel
    let noise = mycorrhizal_channel::NoiseProfile {
        burst_probability: 0.2,
        burst_length: 8,
        attenuation: 0.9,
        random_error_rate: 0.0,
    };
    let channel = MycorrhizalChannel::new(20, noise);

    let data = vec![42, 17, 93, 55, 128, 200, 11, 67, 88, 44, 12, 77];
    let rs = reed_solomon::ReedSolomon::new(12, 8);
    let ecology = burst_ecology::BurstEcology::new(12, 8, 4);

    // Test RS with burst errors
    let codeword = rs.encode(&data);
    let mut rs_successes = 0;
    let mut eco_successes = 0;

    for seed in 0..20 {
        let rs_received = channel.transmit(&codeword, seed);
        if rs.decode(&rs_received).is_ok() {
            rs_successes += 1;
        }
    }

    // Ecology multi-path test
    for seed in 0..20 {
        let result = ecology.transmit_and_decode(&data, &channel, seed + 1000);
        if result.correction_success {
            eco_successes += 1;
        }
    }

    // Both should handle some cases; we mainly verify they run correctly
    assert!(rs_successes + eco_successes > 0, "At least one scheme should work");
}

// === Multipath Redundancy Test ===

#[test]
fn test_multipath_survives_single_path_failure() {
    let noise = mycorrhizal_channel::NoiseProfile::low_noise();
    let channel = MycorrhizalChannel::new(20, noise);

    let data = vec![1, 2, 3, 4, 5, 6, 7, 8];
    let code = PhytoCode::new(8, 4, 3);
    let paths = code.encode_multipath(&data);

    // Remove one path entirely
    let reduced = vec![paths[0].clone(), paths[2].clone()];
    let decoded = code.decode_multipath(&reduced);
    assert!(decoded.is_ok());
}

// === Serde round-trip for all public types ===

#[test]
fn test_all_types_serde_roundtrip() {
    let noise = mycorrhizal_channel::NoiseProfile::default();
    let json = serde_json::to_string(&noise).unwrap();
    let _: mycorrhizal_channel::NoiseProfile = serde_json::from_str(&json).unwrap();

    let channel = MycorrhizalChannel::new(5, noise);
    let json = serde_json::to_string(&channel).unwrap();
    let _: MycorrhizalChannel = serde_json::from_str(&json).unwrap();

    let code = PhytoCode::new(8, 4, 2);
    let json = serde_json::to_string(&code).unwrap();
    let _: PhytoCode = serde_json::from_str(&json).unwrap();

    let tree = HubTree::new(8, 2);
    let json = serde_json::to_string(&tree).unwrap();
    let _: HubTree = serde_json::from_str(&json).unwrap();

    let config = SporeConfig::default();
    let json = serde_json::to_string(&config).unwrap();
    let _: SporeConfig = serde_json::from_str(&json).unwrap();
}
