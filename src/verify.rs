use ark_ec::PairingEngine;
use ark_ff::Field;
use itertools::izip;
use num_traits::Zero;

use crate::error::Error;
use super::{
    proof::*,
    transcript::{TranscriptProtocol, TranscriptAppender},
};

fn verify_gipa<E, T>(
    proof: GipaProof<E>,
    transcript: &mut T,
    r_inv: E::Fr,
) -> Result<(), Error>
where
    E: PairingEngine,
    T: TranscriptProtocol<E>,
{
    let wire_comms = proof.wire_comms;
    let perm_comms = proof.perm_comms;
    let lookup_comms = proof.lookup_comms;
    let mut final_value = proof.final_value;

    // GIPA

    let mut challenges: Vec<E::Fr> = Vec::new();
    let mut challenges_inv: Vec<E::Fr> = Vec::new();

    // for (wc, pc, lc) in izip!(
    //     wire_comms.iter(),
    //     perm_comms.iter(),
    //     lookup_comms.iter(),
    // ) {
    //     // Fiat-Shamir challenge
    //     wc.0.append_in_transcript(transcript);
    //     wc.1.append_in_transcript(transcript);
    //     pc.0.append_in_transcript(transcript);
    //     pc.1.append_in_transcript(transcript);
    //     lc.0.append_in_transcript(transcript);
    //     lc.1.append_in_transcript(transcript);
    //     let x = transcript.challenge_fr("x");
    //     let x_inv = x.inverse().unwrap();

    //     challenges.push(x);
    //     challenges_inv.push(x_inv);
    // }

    let mut vk_comm = E::Fr::zero();
    let mut wk_comm = E::Fr::zero();
    let mut ip = E::Fr::zero();

    for (wc, pc, lc) in izip!(
        wire_comms,
        perm_comms,
        lookup_comms,
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




        challenges.push(x);
        challenges_inv.push(x_inv);
    }



    Ok(())
}