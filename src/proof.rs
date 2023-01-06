use ark_ec::pairing::Pairing;
use ark_serialize::{
    Read, Write,
    CanonicalSerialize, CanonicalDeserialize, SerializationError,
};

use crate::{
    commit::{wire::*, permutation::*, lookup::*},
    transcript::{TranscriptProtocol, TranscriptAppender},
    kzg::KZGOpening,
};

#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Eq, PartialEq)]
pub struct PublicProof<E: PairingEngine> {
    pub q_lo_commit: E::G1Affine,
    pub q_mid_commit: E::G1Affine,
    pub q_hi_commit: E::G1Affine,
    
    // Evaluations
    pub sigma1: E::Fr,
    pub sigma2: E::Fr,
    pub t_tag: E::Fr,
    pub t: E::Fr,
    pub t_next: E::Fr,
}

#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Eq, PartialEq)]
pub struct IndividualProof<E: PairingEngine> {
    pub a_commit: E::G1Affine,
    pub b_commit: E::G1Affine,
    pub c_commit: E::G1Affine,
    pub f_commit: E::G1Affine,
    pub z1_commit: E::G1Affine,
    pub z2_commit: E::G1Affine,
    pub h1_commit: E::G1Affine,
    pub h2_commit: E::G1Affine,
    pub aw_opening: E::G1Affine,
    pub saw_opening: E::G1Affine,

    // Evaluations
    pub a: E::Fr,
    pub b: E::Fr,
    pub c: E::Fr,
    pub f: E::Fr,
    pub z1_next: E::Fr,
    pub z2_next: E::Fr,
    pub h1_next: E::Fr,
    pub h2: E::Fr,
}

pub struct PrimitiveProof<E: PairingEngine> {
    pub wire_comms: WirePrimitiveCommitments<E>,
    pub perm_comms: PermPrimitiveCommitments<E>,
    pub lookup_comms: LookupPrimitiveCommitments<E>,
    pub wire_ips: WirePrimitiveInnerProducts<E>,
    pub perm_ips: PermPrimitiveInnerProducts<E>,
    pub lookup_ips: LookupPrimitiveInnerProducts<E>,
}

#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Eq, PartialEq)]
pub struct GipaFinalValue<E: PairingEngine> {
    // check with wk
    pub final_a_r: E::Fr,
    pub final_b_r: E::Fr,
    pub final_c_r: E::Fr,
    pub final_ab_r: E::Fr,
    pub final_ac_r: E::Fr,
    pub final_bc_r: E::Fr,
    pub final_abc_r: E::Fr,
    pub final_f_r: E::Fr,
    pub final_z2_next_r: E::Fr,
    pub final_h1_next_z2_next_r: E::Fr,
    pub final_h2_z2_next_r: E::Fr,
    // check with vk
    pub final_r: E::Fr,
    pub final_a: E::Fr,
    pub final_b: E::Fr,
    pub final_a_g: E::G1Affine,
    pub final_b_g: E::G1Affine,
    pub final_c_g: E::G1Affine,
    pub final_z1_next: E::Fr,
    pub final_h1_next: E::Fr,
    pub final_h2: E::Fr,
    pub final_z1_g: E::G1Affine,
    pub final_z2_g: E::G1Affine,
    pub final_h1_g: E::G1Affine,
    pub final_f_g: E::G1Affine,
    pub final_h2_g: E::G1Affine,
    pub final_vk: E::G2Affine,
    pub final_wk: E::G1Affine,
}

impl<E, T> TranscriptAppender<E, T> for GipaFinalValue<E>
where
    E: PairingEngine,
    T: TranscriptProtocol<E>,
{
    fn append_in_transcript(&self, transcript: &mut T) {
        transcript.append_fr("final_a_r", &self.final_a_r);
        transcript.append_fr("final_b_r", &self.final_b_r);
        transcript.append_fr("final_c_r", &self.final_c_r);
        transcript.append_fr("final_ab_r", &self.final_ab_r);
        transcript.append_fr("final_ac_r", &self.final_ac_r);
        transcript.append_fr("final_bc_r", &self.final_bc_r);
        transcript.append_fr("final_abc_r", &self.final_abc_r);
        transcript.append_fr("final_f_r", &self.final_f_r);
        transcript.append_fr("final_z2_next_r", &self.final_z2_next_r);
        transcript.append_fr("final_h1_next_z2_next_r", &self.final_h1_next_z2_next_r);
        transcript.append_fr("final_h2_z2_next_r", &self.final_h2_z2_next_r);
        transcript.append_fr("final_r", &self.final_r);
        transcript.append_fr("final_a", &self.final_a);
        transcript.append_fr("final_b", &self.final_b);
        transcript.append_g1("final_a_g", &self.final_a_g);
        transcript.append_g1("final_b_g", &self.final_b_g);
        transcript.append_g1("final_c_g", &self.final_c_g);
        transcript.append_fr("final_z1_next", &self.final_z1_next);
        transcript.append_fr("final_h1_next", &self.final_h1_next);
        transcript.append_fr("final_h2", &self.final_h2);
        transcript.append_g1("final_z1_g", &self.final_z1_g);
        transcript.append_g1("final_z2_g", &self.final_z2_g);
        transcript.append_g1("final_h1_g", &self.final_h1_g);
        transcript.append_g1("final_f_g", &self.final_f_g);
        transcript.append_g1("final_h2_g", &self.final_h2_g);
        transcript.append_g2("final_vk", &self.final_vk);
        transcript.append_g1("final_wk", &self.final_wk);
    }
}

pub struct GipaProof<E: PairingEngine> {
    pub wire_comms: Vec<(WireGipaCommitments<E>, WireGipaCommitments<E>)>,
    pub perm_comms: Vec<(PermGipaCommitments<E>, PermGipaCommitments<E>)>,
    pub lookup_comms: Vec<(LookupGipaCommitments<E>, LookupGipaCommitments<E>)>,
    pub final_value: GipaFinalValue<E>,
    pub vk_opening: KZGOpening<E::G2Affine>,
    pub wk_opening: KZGOpening<E::G1Affine>,
}

pub struct AggregateProof<E: PairingEngine> {
    pub n: u64,
    pub prim_proof: PrimitiveProof<E>,
    pub gipa_proof: GipaProof<E>,
}