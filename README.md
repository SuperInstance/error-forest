# 🌲 error-forest

**Ecological signaling networks modeled as error-correcting codes.**

Mycorrhizal fungal networks — the underground internet connecting forest trees — have evolved sophisticated strategies for reliable communication over noisy, burst-error-prone channels. This library translates those biological strategies into formally analyzable error-correcting codes that **outperform Reed-Solomon in burst-error environments** — exactly what distributed systems experience.

## Why Forests?

Mother trees distribute nutrients and chemical signals through fungal networks that span hundreds of meters. These networks face:
- **Burst errors** from root damage, drought, and chemical interference
- **Multi-path fading** as signals traverse different fungal hyphae
- **Asymmetric attenuation** from varying soil conditions
- **Node failures** when trees die or connections sever

Yet forests maintain remarkably reliable information transfer. The strategies they've evolved over 400 million years map directly to problems in distributed systems.

## Modules

### `mycorrhizal_channel` — Biological Channel Model

Multi-path noisy channels with realistic noise profiles:

```rust
use error_forest::{MycorrhizalChannel, mycorrhizal_channel::NoiseProfile};

let noise = NoiseProfile {
    burst_probability: 0.08,
    burst_length: 6,
    attenuation: 0.9,
    random_error_rate: 0.01,
};
let channel = MycorrhizalChannel::new(20, noise);

let data = vec![1, 2, 3, 4, 5, 6, 7, 8];
let received = channel.transmit(&data, 42);
```

### `phyto_code` — Plant Signaling Codes

Error-correcting codes using Vandermonde-style parity across multiple transmission paths:

```rust
use error_forest::PhytoCode;

let code = PhytoCode::new(8, 6, 3); // 8 data, 6 parity, 3 redundant paths
let data = vec![0x10, 0x20, 0x30, 0x40, 0x50, 0x60, 0x70, 0x80];

// Single-path encode/decode
let codeword = code.encode(&data);
let decoded = code.decode(&codeword).unwrap();

// Multi-path with redundancy
let paths = code.encode_multipath(&data);
let recovered = code.decode_multipath(&paths).unwrap();
```

### `network_shannon` — Channel Capacity Analysis

Compute channel capacity and compare to the Shannon limit:

```rust
use error_forest::{MycorrhizalChannel, network_shannon};

let channel = MycorrhizalChannel::new(20, Default::default());
let result = network_shannon::analyze_channel_capacity(&channel);
println!("Shannon capacity: {:.3} bits/symbol", result.shannon_capacity);
println!("Achieved: {:.3} bits/symbol ({:.1}% efficiency)", 
    result.achieved_capacity, result.efficiency * 100.0);
```

### `burst_ecology` — Ecological Burst-Error Correction

Specialized burst-error correction using interleaved multi-path redundancy:

```rust
use error_forest::burst_ecology::BurstEcology;

let ecology = BurstEcology::new(16, 8, 4);
let data = vec![0u8; 16]; // your data
let encoded = ecology.encode(&data); // 4 interleaved paths
let result = ecology.decode(&encoded);
```

### `hub_tree` — Syndrome Decoding for Hub Networks

Mother trees as parity nodes in a hub-and-spoke topology:

```rust
use error_forest::HubTree;

let tree = HubTree::new(12, 2); // 12 spokes, groups of 3 with parity
let data = vec![10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120];
let parity = tree.encode(&data);

// Detect and locate failures
let result = tree.detect_failed_node(&data, &parity);
if let Some(pos) = result.error_position {
    println!("Node {} compromised!", pos);
}
```

### `distributed_spore` — Spore Gossip Protocol

Eventual consistency with ecological convergence guarantees:

```rust
use error_forest::{SporeConfig, distributed_spore::SporeNetwork};

let config = SporeConfig { fanout: 3, ttl: 10, convergence_target: 0.99 };
let mut network = SporeNetwork::new(20, config);

let spore = network.originate(0, vec![42]);
let result = network.run_until_convergence(vec![spore], 30);
println!("Converged in {} rounds ({} messages)", result.rounds, result.messages_sent);
```

### `reed_solomon` — Reed-Solomon for Comparison

Standard RS codes over GF(256) to benchmark against ecological approaches:

```rust
use error_forest::reed_solomon::ReedSolomon;

let rs = ReedSolomon::new(10, 4); // Corrects up to 2 symbol errors
let data = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
let codeword = rs.encode(&data);
let decoded = rs.decode(&codeword).unwrap();
```

## Core Types

| Type | Description |
|------|-------------|
| `MycorrhizalChannel` | Multi-path noisy channel with burst error profiles |
| `NoiseProfile` | Burst probability, length, attenuation, random error rate |
| `PhytoCode` | Vandermonde-parity error-correcting code |
| `HubTree` | Hub-and-spoke syndrome decoder |
| `SyndromeResult` | Error position, correction status, confidence |
| `SporeConfig` | Gossip fanout, TTL, convergence target |

All public types derive `Serialize`/`Deserialize` via serde.

## Performance Characteristics

- **Burst errors**: Ecological codes (PhytoCode + BurstEcology) outperform Reed-Solomon when burst length exceeds the RS error-correction capacity, thanks to interleaved multi-path redundancy.
- **Multi-path**: Survives single-path failures entirely — data recovers from remaining paths.
- **Convergence**: Spore gossip reaches 99% convergence in O(log N / log(fanout+1)) rounds, matching epidemic spreading models.
- **Detection**: Hub-tree syndrome decoding identifies compromised nodes in O(spokes × parity) time.

## Architecture

```
error-forest/
├── src/
│   ├── lib.rs                 # Library root & re-exports
│   ├── gf256.rs               # GF(2^8) arithmetic (primitive polynomial 0x11D)
│   ├── mycorrhizal_channel.rs # Channel model with burst noise
│   ├── phyto_code.rs          # Vandermonde parity codes
│   ├── network_shannon.rs     # Shannon capacity analysis
│   ├── burst_ecology.rs       # Interleaved burst correction
│   ├── hub_tree.rs            # Hub-spoke syndrome decoding
│   ├── distributed_spore.rs   # Gossip protocol
│   └── reed_solomon.rs        # RS(255, k) for comparison
└── tests/
    └── integration.rs         # 41 integration tests
```

## Dependencies

- `serde` — Serialization for all public types
- No other external dependencies. GF(256) arithmetic, Reed-Solomon, and all algorithms are implemented from scratch.

## Testing

```bash
cargo test          # 46 tests (5 unit + 41 integration)
cargo test -- --nocapture  # with output
```

## License

MIT
