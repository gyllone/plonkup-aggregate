use ark_ec::pairing::Pairing;
use ark_serialize::{CanonicalSerialize, CanonicalDeserialize};

use crate::{
    commit::{wire::*, permutation::*, lookup::*},
    transcript::{TranscriptProtocol, TranscriptAppender},
    kzg::KZGOpening,
};

#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Eq, PartialEq)]
pub struct PublicProof<E: Pairing> {
    pub q_lo_commit: E::G1Affine,
    pub q_mid_commit: E::G1Affine,
    pub q_hi_commit: E::G1Affine,
    
    // Evaluations
    pub sigma1: E::ScalarField,
    pub sigma2: E::ScalarField,
    pub t_tag: E::ScalarField,
    pub t: E::ScalarField,
    pub t_next: E::ScalarField,
}

#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Eq, PartialEq)]
pub struct IndividualProof<E: Pairing> {
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
    pub a: E::ScalarField,
    pub b: E::ScalarField,
    pub c: E::ScalarField,
    pub f: E::ScalarField,
    pub z1_next: E::ScalarField,
    pub z2_next: E::ScalarField,
    pub h1_next: E::ScalarField,
    pub h2: E::ScalarField,
}

pub struct PrimitiveProof<E: Pairing> {
    pub wire_comms: WirePrimitiveCommitments<E>,
    pub perm_comms: PermPrimitiveCommitments<E>,
    pub lookup_comms: LookupPrimitiveCommitments<E>,
    pub wire_ips: WirePrimitiveInnerProducts<E>,
    pub perm_ips: PermPrimitiveInnerProducts<E>,
    pub lookup_ips: LookupPrimitiveInnerProducts<E>,
}

#[derive(CanonicalDeserialize, CanonicalSerialize, derivative::Derivative)]
#[derivative(Clone, Debug, Eq, PartialEq)]
pub struct GipaFinalValue<E: Pairing> {
    // check with wk
    pub a_r: E::ScalarField,
    pub b_r: E::ScalarField,
    pub c_r: E::ScalarField,
    pub ab_r: E::ScalarField,
    pub ac_r: E::ScalarField,
    pub bc_r: E::ScalarField,
    pub abc_r: E::ScalarField,
    pub f_r: E::ScalarField,
    pub z2_next_r: E::ScalarField,
    pub h1_next_z2_next_r: E::ScalarField,
    pub h2_z2_next_r: E::ScalarField,
    // check with vk
    pub a: E::ScalarField,
    pub b: E::ScalarField,
    pub a_g: E::G1Affine,
    pub b_g: E::G1Affine,
    pub c_g: E::G1Affine,
    pub z1_next: E::ScalarField,
    pub h1_next: E::ScalarField,
    pub h2: E::ScalarField,
    pub z1_g: E::G1Affine,
    pub z2_g: E::G1Affine,
    pub h1_g: E::G1Affine,
    pub f_g: E::G1Affine,
    pub h2_g: E::G1Affine,
    pub vk: E::G2Affine,
    pub wk: E::G1Affine,
}

impl<E, T> TranscriptAppender<E, T> for GipaFinalValue<E>
where
    E: Pairing,
    T: TranscriptProtocol<E>,
{
    fn append_in_transcript(&self, transcript: &mut T) {
        transcript.append_fr("a_r", &self.a_r);
        transcript.append_fr("b_r", &self.b_r);
        transcript.append_fr("c_r", &self.c_r);
        transcript.append_fr("ab_r", &self.ab_r);
        transcript.append_fr("ac_r", &self.ac_r);
        transcript.append_fr("bc_r", &self.bc_r);
        transcript.append_fr("abc_r", &self.abc_r);
        transcript.append_fr("f_r", &self.f_r);
        transcript.append_fr("z2_next_r", &self.z2_next_r);
        transcript.append_fr("h1_next_z2_next_r", &self.h1_next_z2_next_r);
        transcript.append_fr("h2_z2_next_r", &self.h2_z2_next_r);
        transcript.append_fr("a", &self.a);
        transcript.append_fr("b", &self.b);
        transcript.append_g1("a_g", &self.a_g);
        transcript.append_g1("b_g", &self.b_g);
        transcript.append_g1("c_g", &self.c_g);
        transcript.append_fr("z1_next", &self.z1_next);
        transcript.append_fr("h1_next", &self.h1_next);
        transcript.append_fr("h2", &self.h2);
        transcript.append_g1("z1_g", &self.z1_g);
        transcript.append_g1("z2_g", &self.z2_g);
        transcript.append_g1("h1_g", &self.h1_g);
        transcript.append_g1("f_g", &self.f_g);
        transcript.append_g1("h2_g", &self.h2_g);
        transcript.append_g2("vk", &self.vk);
        transcript.append_g1("wk", &self.wk);
    }
}

pub struct GipaProof<E: Pairing> {
    pub wire_comms: Vec<(WireGipaCommitments<E>, WireGipaCommitments<E>)>,
    pub perm_comms: Vec<(PermGipaCommitments<E>, PermGipaCommitments<E>)>,
    pub lookup_comms: Vec<(LookupGipaCommitments<E>, LookupGipaCommitments<E>)>,
    pub final_value: GipaFinalValue<E>,
    pub vk_opening: KZGOpening<E::G2Affine>,
    pub wk_opening: KZGOpening<E::G1Affine>,
}

pub struct AggregateProof<E: Pairing> {
    pub n: u64,
    pub prim_proof: PrimitiveProof<E>,
    pub gipa_proof: GipaProof<E>,
}