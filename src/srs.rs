use ark_ec::PairingEngine;

use super::commit::{VKey, WKey};

pub struct KZGCommitSRS<E: PairingEngine> {
    pub h_powers: Vec<E::G2Affine>,
    pub g_powers: Vec<E::G1Affine>,
}

pub struct AggregatorSRS<E: PairingEngine> {
    pub n: usize,
    pub vk: VKey<E>,
    pub wk: WKey<E>,
}
