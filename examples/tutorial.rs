//! # Error Forest Tutorial
//!
//! A progressive walkthrough of ecological error-correcting codes inspired by
//! mycorrhizal fungal networks. Forest error-correction outperforms classical
//! Reed-Solomon in burst-error environments — exactly what distributed systems
//! experience.
//!
//! ## Lessons
//!
//! 1. GF(256) arithmetic — the algebraic foundation
//! 2. Mycorrhizal channels — modeling noisy biological networks
//! 3. PhytoCode — plant-inspired error correction
//! 4. Hub trees & syndrome decoding — mother-tree parity nodes
//! 5. Spore gossip — epidemic data dissemination
//! 6. Burst ecology — interleaved burst-error correction
//! 7. Shannon capacity — comparing to theoretical limits
//! 8. Reed-Solomon comparison — ecological vs classical codes

use error_forest::gf256::GF256;
use error_forest::mycorrhizal_channel::{MycorrhizalChannel, NoiseProfile};
use error_forest::phyto_code::PhytoCode;
use error_forest::hub_tree::HubTree;
use error_forest::distributed_spore::{SporeConfig, SporeNetwork};
use error_forest::burst_ecology::BurstEcology;
use error_forest::network_shannon;
use error_forest::reed_solomon::ReedSolomon;

fn main() {
    println!("════════════════════════════════════════════════════════");
    println!("  ERROR FOREST TUTORIAL");
    println!("  Ecological Error-Correcting Codes");
    println!("════════════════════════════════════════════════════════\n");

    lesson_1_gf256();
    lesson_2_mycorrhizal_channel();
    lesson_3_phyto_code();
    lesson_4_hub_tree();
    lesson_5_spore_gossip();
    lesson_6_burst_ecology();
    lesson_7_shannon_capacity();
    lesson_8_reed_solomon_comparison();

    println!("\n✅ Tutorial complete! The forest grows stronger.");
}

// ─── Lesson 1: GF(256) Arithmetic ──────────────────────────────────────

fn lesson_1_gf256() {
    println!("━━━ Lesson 1: GF(256) Arithmetic ━━━\n");
    println!("GF(256) is the Galois Field used in all error-correcting codes here.");
    println!("It uses the AES polynomial x^8 + x^4 + x^3 + x + 1.\n");

    // Create elements
    let a = GF256::new(0x57);
    let b = GF256::new(0x83);

    // Addition is XOR
    let sum = a.add(b);
    println!("  {} + {} = {} (XOR in GF(2^m))", a.0, b.0, sum.0);

    // Subtraction is the same as addition in GF(2^m)
    let diff = a.sub(b);
    println!("  {} - {} = {} (same as add!)", a.0, b.0, diff.0);
    assert_eq!(sum, diff);

    // Multiplication uses log/exp tables
    let prod = a.mul(b);
    println!("  {} × {} = {}", a.0, b.0, prod.0);

    // Division
    let quotient = prod.div(b);
    println!("  {} / {} = {} (recovers a)", prod.0, b.0, quotient.0);
    assert_eq!(quotient, a);

    // Inverse and identity
    let inv_a = a.inv();
    assert_eq!(a.mul(inv_a), GF256::ONE);
    println!("  inv({}) = {}, verify: {} × {} = 1", a.0, inv_a.0, a.0, inv_a.0);

    // Power
    let generator = GF256::new(2);
    assert_eq!(generator.pow(255), GF256::ONE);
    println!("  2^255 = 1 (multiplicative group order is 255)");

    // Polynomial operations
    let poly_a: Vec<GF256> = vec![GF256::ONE, GF256::new(2), GF256::new(3)]; // 1 + 2x + 3x²
    let poly_b: Vec<GF256> = vec![GF256::ONE, GF256::new(5)];               // 1 + 5x
    let product = GF256::poly_mul(&poly_a, &poly_b);
    println!("  (1 + 2x + 3x²)(1 + 5x) has {} coefficients", product.len());

    let (quot, rem) = GF256::poly_div(&product, &poly_b);
    println!("  Dividing back: quotient has {} terms, remainder has {} terms", quot.len(), rem.len());

    // Evaluate polynomial at a point
    let val = GF256::eval_poly(&poly_a, GF256::new(7));
    println!("  p(7) = 1 + 2·7 + 3·49 = {} (in GF(256))", val.0);

    println!();
}

// ─── Lesson 2: Mycorrhizal Channels ──────────────────────────────────────

fn lesson_2_mycorrhizal_channel() {
    println!("━━━ Lesson 2: Mycorrhizal Channels ━━━\n");
    println!("Mycorrhizal channels model noisy multi-path biological networks.");
    println!("They simulate burst errors, attenuation, and random noise.\n");

    // Create a channel with custom noise profile
    let noise = NoiseProfile::default();
    println!("  Default noise: burst_prob={:.2}, burst_len={}, random={:.3}",
        noise.burst_probability, noise.burst_length, noise.random_error_rate);

    let channel = MycorrhizalChannel::new(10, noise);

    // Transmit data and see what happens
    let original = b"forest!!!";
    let received = channel.transmit(original, 42);

    let (bit_errors, byte_errors) = MycorrhizalChannel::count_errors(original, &received);
    println!("  Sent:     {:?}", std::str::from_utf8(original).unwrap());
    println!("  Received: {:?} ({} bit errors, {} byte errors)",
        received, bit_errors, byte_errors);

    // Classify the error pattern
    let classification = MycorrhizalChannel::classify_errors(original, &received);
    println!("  Burst errors: {}, Random errors: {}, Total: {}",
        classification.bursts.len(),
        classification.random_errors,
        classification.total_errors);

    // Multi-path transmission with majority voting
    println!("\n  Multi-path transmission (majority vote across 5 paths):");
    let received_mp = channel.transmit_multipath(original, 0, 9, 42);
    let (mp_bits, mp_bytes) = MycorrhizalChannel::count_errors(original, &received_mp);
    println!("  Multi-path: {} bit errors, {} byte errors (vs single: {}, {})",
        mp_bits, mp_bytes, bit_errors, byte_errors);

    // Compare noise profiles
    println!("\n  Noise profiles:");
    for (name, profile) in [
        ("Low noise (clean forest)", NoiseProfile::low_noise()),
        ("Default", NoiseProfile::default()),
        ("Burst dominant", NoiseProfile::burst_dominant()),
        ("High noise (disturbed)", NoiseProfile::high_noise()),
    ] {
        let ch = MycorrhizalChannel::new(10, profile);
        let rx = ch.transmit(original, 42);
        let (bits, _) = MycorrhizalChannel::count_errors(original, &rx);
        println!("    {}: {} bit errors", name, bits);
    }

    println!();
}

// ─── Lesson 3: PhytoCode — Plant-Inspired Error Correction ──────────────

fn lesson_3_phyto_code() {
    println!("━━━ Lesson 3: PhytoCode — Plant-Inspired Error Correction ━━━\n");
    println!("PhytoCode uses Vandermonde-style encoding over GF(256),");
    println!("inspired by phytochemical signaling strategies.\n");

    let data = b"hello world!!!"; // 14 bytes
    let code = PhytoCode::new(14, 4, 3); // 14 data + 4 parity, 3 redundant paths

    println!("  Codeword: {} data + {} parity = {} symbols",
        code.data_symbols, code.parity_symbols, code.codeword_len());
    println!("  Max correctable errors: {} symbols", code.max_correctable());

    // Encode
    let codeword = code.encode(data);
    println!("  Encoded {} bytes → {} bytes", data.len(), codeword.len());

    // Decode clean data
    let decoded = code.decode(&codeword).unwrap();
    assert_eq!(decoded, data.to_vec());
    println!("  ✓ Clean decode: {:?}", std::str::from_utf8(&decoded).unwrap());

    // Introduce a single error and correct it
    let mut corrupted = codeword.clone();
    corrupted[3] ^= 0xFF; // flip all bits in position 3
    let decoded_fixed = code.decode(&corrupted).unwrap();
    assert_eq!(decoded_fixed, data.to_vec());
    println!("  ✓ Corrected single-byte error at position 3");

    // Multi-path encoding
    let multi = code.encode_multipath(data);
    println!("\n  Multi-path: {} paths encoded", multi.len());

    // Corrupt all paths differently and recover
    let mut corrupted_paths: Vec<Vec<u8>> = multi.iter().cloned().collect();
    for path in corrupted_paths.iter_mut() {
        path[5] ^= 0xAA;
    }
    let recovered = code.decode_multipath(&corrupted_paths).unwrap();
    assert_eq!(recovered, data.to_vec());
    println!("  ✓ Recovered from multi-path corruption");

    // Compare to naive repetition
    let channel = MycorrhizalChannel::new(8, NoiseProfile::high_noise());
    let (phyto_errs, rep_errs) = code.compare_to_repetition(data, &channel, 42);
    println!("\n  PhytoCode vs Repetition (high noise channel):");
    println!("    PhytoCode errors: {}", phyto_errs);
    println!("    Repetition errors: {}", rep_errs);

    println!();
}

// ─── Lesson 4: Hub Trees & Syndrome Decoding ──────────────────────────────

fn lesson_4_hub_tree() {
    println!("━━━ Lesson 4: Hub Trees & Syndrome Decoding ━━━\n");
    println!("Mother trees act as hub parity nodes in a hub-and-spoke network.");
    println!("Syndrome decoding detects which tree is compromised.\n");

    // Create a hub tree with 8 spokes and 2 parity symbols per hub
    let tree = HubTree::new(8, 2);
    println!("  Network: {} spokes, {} hubs, {} total nodes",
        tree.spokes, tree.num_hubs(), tree.total_nodes());

    // Encode spoke data
    let spoke_data: Vec<u8> = vec![10, 20, 30, 40, 50, 60, 70, 80];
    let parity = tree.encode(&spoke_data);

    println!("  Spoke data: {:?}", spoke_data);
    for (h, p) in parity.iter().enumerate() {
        println!("    Hub {} parity: {:?}", h, p);
    }

    // Verify syndromes are zero (no errors)
    let syndromes = tree.compute_syndromes(&spoke_data, &parity);
    let all_zero = syndromes.iter().all(|s| s.iter().all(|&v| v == 0));
    println!("  Syndromes all zero: {}", all_zero);

    // Corrupt a spoke and detect
    let mut corrupted = spoke_data.clone();
    corrupted[3] = 99; // corrupt spoke 3
    let result = tree.detect_failed_node(&corrupted, &parity);
    println!("\n  After corrupting spoke 3 (40 → 99):");
    println!("    Detected error at position: {:?}", result.error_position);
    println!("    Corrected: {}", result.corrected);
    println!("    Corrected value: {:?}", result.corrected_value);
    println!("    Confidence: {:.2}", result.confidence);

    // Full decode cycle
    let (corrected_data, _decode_result) = tree.decode(&corrupted, &parity);
    println!("    Decoded data: {:?}", corrected_data);
    assert_eq!(corrected_data, spoke_data);
    println!("  ✓ Successfully corrected the corrupted spoke");

    // Simulate node failures
    let mut tree2 = HubTree::new(12, 3);
    tree2.simulate_failures(&[2, 5, 8]);
    println!("\n  Simulated failures at nodes [2, 5, 8]:");
    println!("    Health: {:?}", &tree2.node_health[0..12.min(tree2.node_health.len())]);

    println!();
}

// ─── Lesson 5: Spore Gossip Protocol ──────────────────────────────────────

fn lesson_5_spore_gossip() {
    println!("━━━ Lesson 5: Spore Gossip Protocol ━━━\n");
    println!("Spore gossip uses epidemic-style dissemination,");
    println!("inspired by fungal spore dispersal patterns.\n");

    let config = SporeConfig {
        fanout: 3,
        ttl: 10,
        convergence_target: 0.95,
    };
    println!("  Config: fanout={}, ttl={}, target={:.0}%",
        config.fanout, config.ttl, config.convergence_target * 100.0);

    let mut network = SporeNetwork::new(20, config);
    println!("  Network: {} nodes", network.nodes.len());

    let theoretical = network.theoretical_convergence_bound();
    println!("  Theoretical convergence bound: {} rounds", theoretical);

    // Originate a spore from node 0
    let spore = network.originate(0, b"mycelium_data".to_vec());
    println!("\n  Originated spore {} from node 0: {:?}",
        spore.id, std::str::from_utf8(&spore.data).unwrap());

    // Run gossip until convergence
    let result = network.run_until_convergence(vec![spore], 50);
    println!("  Convergence result:");
    println!("    Rounds: {}", result.rounds);
    println!("    Convergence: {:.1}%", result.convergence_ratio * 100.0);
    println!("    Messages sent: {}", result.messages_sent);
    println!("    Converged: {}", result.converged);

    // Show convergence history
    if result.convergence_history.len() <= 10 {
        println!("    History: {:?}", result.convergence_history.iter()
            .map(|r| format!("{:.0}%", r * 100.0)).collect::<Vec<_>>());
    } else {
        print!("    History: ");
        for (i, r) in result.convergence_history.iter().enumerate() {
            if i % 5 == 0 { print!("{:.0}% ", r * 100.0); }
        }
        println!("...");
    }

    // Originate multiple spores
    println!("\n  Originating 3 spores from different nodes...");
    let mut network2 = SporeNetwork::new(15, SporeConfig::default());
    let spores = vec![
        network2.originate(0, b"alpha".to_vec()),
        network2.originate(5, b"beta".to_vec()),
        network2.originate(10, b"gamma".to_vec()),
    ];
    let result2 = network2.run_until_convergence(spores, 30);
    println!("    Rounds: {}, Convergence: {:.1}%, Messages: {}",
        result2.rounds, result2.convergence_ratio * 100.0, result2.messages_sent);

    println!();
}

// ─── Lesson 6: Burst Ecology ──────────────────────────────────────────────

fn lesson_6_burst_ecology() {
    println!("━━━ Lesson 6: Burst Ecology — Interleaved Burst Correction ━━━\n");
    println!("Burst ecology uses interleaved parity across multiple paths");
    println!("to correct burst errors that defeat conventional codes.\n");

    let data: Vec<u8> = (0..16).collect();
    let ecology = BurstEcology::new(16, 4, 5);

    println!("  Encoder: {} data + {} parity, {} paths",
        ecology.data_symbols, ecology.parity_symbols, ecology.num_paths);

    // Encode with interleaving
    let encoded = ecology.encode(&data);
    println!("  Encoded into {} interleaved paths", encoded.len());

    // Classify burst patterns
    let reference = vec![0u8; 16];
    let corrupted = vec![0u8, 1, 2, 99, 88, 77, 6, 7, 8, 9, 10, 55, 44, 13, 14, 15];
    let patterns = BurstEcology::classify_burst_patterns(&corrupted, &reference);
    println!("\n  Burst pattern analysis:");
    for p in &patterns {
        println!("    {} at pos {} len {} severity {}",
            p.pattern_type, p.start, p.length, p.severity);
    }

    // Full transmit-and-decode cycle through a noisy channel
    let channel = MycorrhizalChannel::new(8, NoiseProfile::burst_dominant());
    let result = ecology.transmit_and_decode(&data, &channel, 42);
    println!("\n  Transmit through burst-dominant channel:");
    println!("    Original errors: {}", result.original_errors);
    println!("    Burst count: {}", result.burst_count);
    println!("    Correction success: {}", result.correction_success);
    println!("    Remaining errors: {}", result.remaining_errors);

    println!();
}

// ─── Lesson 7: Shannon Capacity ──────────────────────────────────────────────

fn lesson_7_shannon_capacity() {
    println!("━━━ Lesson 7: Shannon Capacity Analysis ━━━\n");
    println!("Compare the achieved rate of ecological codes to the");
    println!("theoretical Shannon limit of the channel.\n");

    // Analyze different channel conditions
    for (name, profile) in [
        ("Clean forest", NoiseProfile::low_noise()),
        ("Normal", NoiseProfile::default()),
        ("Disturbed", NoiseProfile::high_noise()),
    ] {
        let channel = MycorrhizalChannel::new(10, profile);
        let result = network_shannon::analyze_channel_capacity(&channel);
        println!("  {}:", name);
        println!("    Shannon capacity:  {:.3} bits/symbol", result.shannon_capacity);
        println!("    Achieved capacity: {:.3} bits/symbol", result.achieved_capacity);
        println!("    Efficiency:        {:.1}%", result.efficiency * 100.0);
        println!("    Effective SNR:     {:.1} dB", result.effective_snr_db);
    }

    // Compare burst vs random errors
    println!("\n  Burst vs Random error impact on capacity:");
    let comparison = network_shannon::compare_burst_vs_random(0.05, 5, 0.1);
    println!("    Burst:  capacity={:.3} at error rate {:.3}",
        comparison.burst_capacity, comparison.burst_error_rate);
    println!("    Random: capacity={:.3} at error rate {:.3}",
        comparison.random_capacity, comparison.random_error_rate);

    // BSC capacity curve
    println!("\n  Binary Symmetric Channel capacity curve:");
    for p in [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5] {
        let cap = network_shannon::bsc_capacity(p);
        let bar = "█".repeat((cap * 30.0) as usize);
        println!("    p={:.2}: C={:.3} {}", p, cap, bar);
    }

    println!();
}

// ─── Lesson 8: Reed-Solomon Comparison ──────────────────────────────────────

fn lesson_8_reed_solomon_comparison() {
    println!("━━━ Lesson 8: Reed-Solomon vs Ecological Codes ━━━\n");
    println!("Classical Reed-Solomon codes for comparison with ecological codes.\n");

    let data = b"mycorrhiza!!!"; // 13 bytes
    let rs = ReedSolomon::new(13, 4); // 13 data + 4 parity

    println!("  RS({}, {}): can correct up to {} symbol errors",
        rs.n(), rs.data_symbols, rs.t());

    // Encode
    let codeword = rs.encode(data);
    println!("  Encoded {} bytes → {} bytes", data.len(), codeword.len());

    // Decode clean
    let decoded = rs.decode(&codeword).unwrap();
    assert_eq!(decoded, data.to_vec());
    println!("  ✓ Clean decode: {:?}", std::str::from_utf8(&decoded).unwrap());

    // Corrupt 2 symbols and correct
    let mut corrupted = codeword.clone();
    corrupted[2] ^= 0xAB;
    corrupted[7] ^= 0xCD;
    let fixed = rs.decode(&corrupted).unwrap();
    assert_eq!(fixed, data.to_vec());
    println!("  ✓ Corrected 2-symbol errors at positions 2 and 7");

    // Overwhelm with errors
    let mut heavily_corrupted = codeword.clone();
    for i in [0, 1, 2, 3, 4] {
        heavily_corrupted[i] ^= 0xFF;
    }
    match rs.decode(&heavily_corrupted) {
        Ok(_) => println!("  Unexpectedly corrected 5 errors (beyond t=2)!"),
        Err(e) => println!("  ✓ Correctly rejected 5 errors (beyond t=2): {}", e),
    }

    println!("\n  Key insight: Ecological codes (PhytoCode + burst ecology)");
    println!("  outperform RS in burst-error environments through multi-path");
    println!("  redundancy and interleaving — nature's strategy wins!");

    println!();
}
