use ark_ec::pairing::Pairing;

use super::commit::{VKey, WKey};

pub struct KZGCommitSRS<E: Pairing> {
    pub g_powers: Vec<E::G1Affine>,
    pub h_powers: Vec<E::G2Affine>,
}

pub struct AggregatorSRS<E: Pairing> {
    pub n: usize,
    pub vk: VKey<E>,
    pub wk: WKey<E>,
}

pub struct KZGVerifySRS<E: Pairing> {
    pub g1_base: E::G1Affine,
    pub g2_base: E::G2Affine,
    pub g: E::G1Affine,
    pub h: E::G2Affine,
}
