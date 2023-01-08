use ark_ec::pairing::Pairing;
use ark_ff::Field;
use itertools::izip;
use num_traits::One;

use crate::{
    error::Error,
    commit::{wire::*, permutation::*, lookup::*},
    proof::*,
    srs::KZGVerifySRS,
    kzg::{verify_kzg_opening_on_g1, verify_kzg_opening_on_g2},
    transcript::{TranscriptProtocol, TranscriptAppender},
    util::poly_eval_from_transcript,
};

pub fn verify_aggregated_proof<E, T>(
    vs: &KZGVerifySRS<E>,
    proof: AggregateProof<E>,
    transcript: &mut T,
) -> Result<(), Error>
where
    E: Pairing,
    T: TranscriptProtocol<E>,
{
    // Primitive part

    let prim_proof = proof.prim_proof;

    prim_proof.wire_comms.append_in_transcript(transcript);
    prim_proof.perm_comms.append_in_transcript(transcript);
    prim_proof.lookup_comms.append_in_transcript(transcript);

    // Derive a random scalar to perform a linear combination of proofs
    let r = transcript.challenge_fr("r");

    prim_proof.wire_ips.append_in_transcript(transcript);
    prim_proof.perm_ips.append_in_transcript(transcript);
    prim_proof.lookup_ips.append_in_transcript(transcript);

    let final_wc = WireFinalCommitments::from_primitive(
        prim_proof.wire_comms,
        prim_proof.wire_ips,
    );
    let final_pc = PermFinalCommitments::from_primitive(
        prim_proof.perm_comms,
        prim_proof.perm_ips,
    );
    let final_lc = LookupFinalCommitments::from_primitive(
        prim_proof.lookup_comms,
        prim_proof.lookup_ips,
    );

    // GIPA part

    verify_gipa(
        vs,
        final_wc,
        final_pc,
        final_lc,
        proof.gipa_proof,
        transcript,
        r,
    )
}

fn verify_gipa<E, T>(
    vs: &KZGVerifySRS<E>,
    mut final_wc: WireFinalCommitments<E>,
    mut final_pc: PermFinalCommitments<E>,
    mut final_lc: LookupFinalCommitments<E>,
    proof: GipaProof<E>,
    transcript: &mut T,
    r: E::ScalarField,
) -> Result<(), Error>
where
    E: Pairing,
    T: TranscriptProtocol<E>,
{
    let mut challenges: Vec<E::ScalarField> = Vec::new();
    let mut challenges_inv: Vec<E::ScalarField> = Vec::new();

    for (wc, pc, lc) in izip!(
        proof.wire_comms,
        proof.perm_comms,
        proof.lookup_comms,
    ) {
        // Fiat-Shamir challenge
        wc.0.append_in_transcript(transcript);
        wc.1.append_in_transcript(transcript);
        pc.0.append_in_transcript(transcript);
        pc.1.append_in_transcript(transcript);
        lc.0.append_in_transcript(transcript);
        lc.1.append_in_transcript(transcript);
        let x = transcript.challenge_fr("x");
        let x_inv = x.inverse().unwrap();

        final_wc.merge(wc.0, wc.1, &x, &x_inv);
        final_pc.merge(pc.0, pc.1, &x, &x_inv);
        final_lc.merge(lc.0, lc.1, &x, &x_inv);

        challenges.push(x);
        challenges_inv.push(x_inv);
    }

    challenges.reverse();
    challenges_inv.reverse();

    // Fiat-Shamir challenge
    proof.final_value.append_in_transcript(transcript);
    // KZG challenge point
    let z = transcript.challenge_fr("z");

    // Check vk
    verify_vk(
        vs,
        &challenges,
        proof.final_value.vk,
        proof.vk_opening,
        z,
    )?;

    // Check wk
    verify_wk(
        vs,
        &challenges_inv,
        proof.final_value.wk,
        proof.wk_opening,
        z,
        r.inverse().unwrap(),
    )?;

    // Final gipa check
    let final_one = poly_eval_from_transcript(
        &challenges,
        E::ScalarField::one(),
        E::ScalarField::one(),
    );
    let final_r = poly_eval_from_transcript(
        &challenges_inv,
        E::ScalarField::one(),
        r,
    );

    final_wc.check(&proof.final_value, &final_one, &final_r);
    final_pc.check(&proof.final_value, &final_one, &final_r);
    final_lc.check(&proof.final_value, &final_one, &final_r);

    Ok(())
}

fn verify_vk<E: Pairing>(
    vs: &KZGVerifySRS<E>,
    challenges: &[E::ScalarField],
    final_vk: E::G2Affine,
    vk_opening: E::G2Affine,
    z: E::ScalarField,
) -> Result<(), Error> {
    let poly_eval = poly_eval_from_transcript(
        challenges,
        z,
        E::ScalarField::one(),
    );

    verify_kzg_opening_on_g2::<E>(
        vs.g1_base,
        vs.g2_base,
        vs.g,
        final_vk,
        vk_opening,
        poly_eval,
        z,
    )
}

fn verify_wk<E: Pairing>(
    vs: &KZGVerifySRS<E>,
    challenges_inv: &[E::ScalarField],
    final_wk: E::G1Affine,
    wk_opening: E::G1Affine,
    z: E::ScalarField,
    r_inv: E::ScalarField,
) -> Result<(), Error> {
    let poly_eval = poly_eval_from_transcript(
        challenges_inv,
        z,
        r_inv,
    );

    verify_kzg_opening_on_g1::<E>(
        vs.g1_base,
        vs.g2_base,
        vs.h,
        final_wk,
        wk_opening,
        poly_eval,
        z,
    )
}
