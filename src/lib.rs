//! # Error Forest
//!
//! Ecological signaling networks (mycorrhizal fungal networks) modeled as
//! error-correcting codes. Forest error-correction outperforms Reed-Solomon
//! in burst-error environments — exactly what distributed systems experience.

pub mod gf256;
pub mod mycorrhizal_channel;
pub mod phyto_code;
pub mod network_shannon;
pub mod burst_ecology;
pub mod hub_tree;
pub mod distributed_spore;
pub mod reed_solomon;

pub use mycorrhizal_channel::{MycorrhizalChannel, NoiseProfile};
pub use phyto_code::PhytoCode;
pub use hub_tree::{HubTree, SyndromeResult};
pub use distributed_spore::SporeConfig;
