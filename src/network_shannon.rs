//! Channel capacity computation for mycorrhizal-style networks.
//! Compare achieved rates to theoretical Shannon limit.

use crate::mycorrhizal_channel::{MycorrhizalChannel, NoiseProfile};
use serde::{Serialize, Deserialize};

/// Shannon channel capacity result.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CapacityResult {
    /// Theoretical Shannon capacity (bits per symbol).
    pub shannon_capacity: f64,
    /// Achieved capacity by the ecological code (bits per symbol).
    pub achieved_capacity: f64,
    /// Channel efficiency (achieved / shannon).
    pub efficiency: f64,
    /// Effective SNR of the channel.
    pub effective_snr_db: f64,
}

/// Compute binary entropy H(p).
fn binary_entropy(p: f64) -> f64 {
    if p <= 0.0 || p >= 1.0 {
        return 0.0;
    }
    -p * p.log2() - (1.0 - p) * (1.0 - p).log2()
}

/// Shannon capacity for a binary symmetric channel with error probability p.
pub fn bsc_capacity(error_prob: f64) -> f64 {
    1.0 - binary_entropy(error_prob)
}

/// Compute effective error probability from a noise profile.
pub fn effective_error_rate(noise: &NoiseProfile) -> f64 {
    // Combine burst and random error rates
    let burst_contribution = noise.burst_probability * noise.burst_length as f64;
    let total_rate = noise.random_error_rate + burst_contribution;
    total_rate.min(1.0)
}

/// Compute SNR from error probability (BSC model: SNR relates to crossover probability).
pub fn error_rate_to_snr_db(error_prob: f64) -> f64 {
    if error_prob <= 0.0 {
        return f64::INFINITY;
    }
    if error_prob >= 0.5 {
        return -f64::INFINITY;
    }
    // Approximate: Q^{-1}(p) ≈ sqrt(2) * erfc^{-1}(2p)
    // Simplified: SNR ≈ (1-2p)^2 / (4p(1-p))
    let snr_linear = (1.0 - 2.0 * error_prob).powi(2) / (4.0 * error_prob * (1.0 - error_prob));
    10.0 * snr_linear.log10()
}

/// Full capacity analysis for a mycorrhizal channel.
pub fn analyze_channel_capacity(channel: &MycorrhizalChannel) -> CapacityResult {
    let error_rate = effective_error_rate(&channel.noise_profile);
    let shannon_cap = bsc_capacity(error_rate);

    // Estimate achieved capacity by simulation
    let achieved = estimate_achieved_capacity(channel, 1000, 42);

    let efficiency = if shannon_cap > 0.0 {
        achieved / shannon_cap
    } else {
        0.0
    };

    CapacityResult {
        shannon_capacity: shannon_cap,
        achieved_capacity: achieved,
        efficiency: efficiency,
        effective_snr_db: error_rate_to_snr_db(error_rate),
    }
}

/// Estimate achieved capacity by simulation.
/// Transmit random data, measure successful decode rate.
pub fn estimate_achieved_capacity(
    channel: &MycorrhizalChannel,
    num_trials: usize,
    seed: u64,
) -> f64 {
    let mut rng = crate::mycorrhizal_channel::SimpleRng::new(seed);
    let data_len = 16;
    let mut total_bits = 0usize;
    let mut error_bits = 0usize;

    for _trial in 0..num_trials {
        // Generate random data
        let data: Vec<u8> = (0..data_len).map(|_| rng.next_u8()).collect();
        let trial_seed = rng.next_u64();

        let received = channel.transmit(&data, trial_seed);

        let (bit_err, _) = MycorrhizalChannel::count_errors(&data, &received);
        total_bits += data_len * 8;
        error_bits += bit_err;
    }

    let ber = error_bits as f64 / total_bits as f64;
    // Capacity = 1 - H(BER) for BSC
    bsc_capacity(ber.min(0.5))
}

/// Compare capacity under burst vs random error models.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct BurstVsRandom {
    pub burst_capacity: f64,
    pub random_capacity: f64,
    pub burst_error_rate: f64,
    pub random_error_rate: f64,
}

/// Analyze how burst errors affect capacity differently than random errors.
pub fn compare_burst_vs_random(
    burst_prob: f64,
    burst_length: usize,
    random_prob: f64,
) -> BurstVsRandom {
    let burst_noise = NoiseProfile {
        burst_probability: burst_prob,
        burst_length,
        attenuation: 1.0,
        random_error_rate: 0.0,
    };
    let random_noise = NoiseProfile {
        burst_probability: 0.0,
        burst_length: 0,
        attenuation: 1.0,
        random_error_rate: random_prob,
    };

    let burst_rate = effective_error_rate(&burst_noise);
    let random_rate = effective_error_rate(&random_noise);

    BurstVsRandom {
        burst_capacity: bsc_capacity(burst_rate),
        random_capacity: bsc_capacity(random_rate),
        burst_error_rate: burst_rate,
        random_error_rate: random_rate,
    }
}
