use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_poly::{UVPolynomial, univariate::DensePolynomial};
use itertools::Itertools;
use num_traits::One;

use crate::{
    error::Error,
    commit::{VKey, WKey, wire::*, permutation::*, lookup::*},
    transcript::{TranscriptProtocol, TranscriptAppender},
    kzg::{create_kzg_opening, KZGOpening},
    srs::{KZGCommitSRS, AggregatorSRS},
};

pub fn aggregate_proofs<E, T, I>(
    ks: &KZGCommitSRS<E>,
    aggs: &AggregatorSRS<E>,
    proofs: Vec<IndividualProof<E>>,
    instances: I,
    transcript: &mut T,
) -> Result<AggregateProof<E>, Error>
where
    E: PairingEngine,
    T: TranscriptProtocol<E>,
    I: IntoIterator<Item = Vec<E::Fr>>,
{
    let n = proofs.len();
    assert!(n.is_power_of_two(), "length of proofs must be power of 2");
    assert_eq!(n, aggs.n);

    let (evaluations, commits): (Vec<_>, Vec<_>) =
        proofs.into_iter().map(|proof| {
            (
                (
                    proof.a,
                    proof.b,
                    proof.c,
                    proof.f,
                    proof.z1_next,
                    proof.z2_next,
                    proof.h1_next,
                    proof.h2,
                ),
                (
                    proof.a_commit,
                    proof.b_commit,
                    proof.c_commit,
                    proof.z1_commit,
                    proof.z2_commit,
                    proof.h1_commit,
                    proof.h2_commit,
                    proof.f_commit,
                )
            )
        }).unzip();

    let (a, b, c, f, z1_next, z2_next, h1_next, h2):
        (Vec<_>, Vec<_>, Vec<_>, Vec<_>, Vec<_>, Vec<_>, Vec<_>, Vec<_>) =
            evaluations.into_iter().map(|evaluation| {
                (
                    evaluation.0,
                    evaluation.1,
                    evaluation.2,
                    evaluation.3,
                    evaluation.4,
                    evaluation.5,
                    evaluation.6,
                    evaluation.7,
                )
            }).multiunzip();
    
    let (a_g, b_g, c_g, z1_g, z2_g, h1_g, f_g, h2_g):
        (Vec<_>, Vec<_>, Vec<_>, Vec<_>, Vec<_>, Vec<_>, Vec<_>, Vec<_>) =
            commits.into_iter().map(|commit| {
                (
                    commit.0,
                    commit.1,
                    commit.2,
                    commit.3,
                    commit.4,
                    commit.5,
                    commit.6,
                    commit.7,
                )
            }).multiunzip();

    let ab = scalar_mul(&a, &b);
    let ac = scalar_mul(&a, &c);
    let bc = scalar_mul(&b, &c);
    let abc = scalar_mul(&ab, &c);
    let h1_next_z2_next = scalar_mul(&h1_next, &z2_next);
    let h2_z2_next = scalar_mul(&h2, &z2_next);

    // Primitive Part

    // Wire commit
    let wire_comms = wire_primitive_commit(
        &a,
        &b,
        &c,
        &ab,
        &ac,
        &bc,
        &abc,
        &a_g,
        &b_g,
        &c_g,
        &aggs.vk,
        &aggs.wk,
    );
    wire_comms.append_in_transcript(transcript);
    
    // Permutation commit
    let perm_comms = perm_primitive_commit(
        &z1_next,
        &z1_g,
        &aggs.vk,
    );
    perm_comms.append_in_transcript(transcript);

    // Lookup commit
    let lookup_comms = lookup_primitive_commit(
        &h1_next,
        &h2,
        &f,
        &z2_next,
        &h1_next_z2_next,
        &h2_z2_next,
        &z2_g,
        &h1_g,
        &f_g,
        &h2_g,
        &aggs.vk,
        &aggs.wk,
    );
    lookup_comms.append_in_transcript(transcript);

    // Derive a random scalar to perform a linear combination of proofs
    let r = transcript.challenge_fr("r");
    let r_inv = r.inverse().unwrap();

    // 1,r, r^2, r^3, r^4 ...
    let r_vec = structured_scalar_power(n, &r);
    // 1,r^-1, r^-2, r^-3
    let r_inv_vec: Vec<_> = cfg_iter!(&r_vec)
        .map(|ri| ri.inverse().unwrap())
        .collect();

    // Shift vectors by r
    let a_r = scalar_mul(&a, &r_vec);
    let b_r = scalar_mul(&b, &r_vec);
    let c_r = scalar_mul(&c, &r_vec);
    let ab_r = cover_scalar_mul(ab, &r_vec);
    let ac_r = cover_scalar_mul(ac, &r_vec);
    let bc_r = cover_scalar_mul(bc, &r_vec);
    let abc_r = cover_scalar_mul(abc, &r_vec);
    let f_r = cover_scalar_mul(f, &r_vec);
    let z2_next_r = cover_scalar_mul(z2_next, &r_vec);
    let h1_next_z2_next_r = cover_scalar_mul(h1_next_z2_next, &r_vec);
    let h2_z2_next_r = cover_scalar_mul(h2_z2_next, &r_vec);

    // Wire inner product
    let wire_ips = wire_primitive_inner_product(
        &r_vec,
        &a_r,
        &b_r,
        &c_r,
        &ab_r,
        &ac_r,
        &bc_r,
        &abc_r,
        &a_g,
        &b_g,
        &c_g,
    );
    wire_ips.append_in_transcript(transcript);

    // Permutation inner product
    let perm_ips = perm_primitive_inner_product(
        &r_vec,
        &a_r,
        &b_r,
        &c_r,
        &ab_r,
        &ac_r,
        &bc_r,
        &abc_r,
        &z1_next,
        &z1_g,
    );
    perm_ips.append_in_transcript(transcript);

    // Lookup inner product
    let lookup_ips = lookup_primitive_inner_product(
        &r_vec,
        &h1_next,
        &h2,
        &f_r,
        &z2_next_r,
        &h1_next_z2_next_r,
        &h2_z2_next_r,
        &z2_g,
        &h1_g,
        &f_g,
        &h2_g,
    );
    lookup_ips.append_in_transcript(transcript);

    let prim_proof = PrimitiveProof {
        wire_comms,
        perm_comms,
        lookup_comms,
        wire_ips,
        perm_ips,
        lookup_ips,
    };

    // Gipa Part

    let vk = aggs.vk.clone();
    // Shift wk by r_inv
    let wk_r_inv = aggs.wk.scale(&r_inv_vec);

    let gipa_proof = prove_gipa(
        ks,
        n,
        vk,
        wk_r_inv,
        r_vec,
        a,
        b,
        z1_next,
        h1_next,
        h2,
        a_r,
        b_r,
        c_r,
        ab_r,
        ac_r,
        bc_r,
        abc_r,
        f_r,
        z2_next_r,
        h1_next_z2_next_r,
        h2_z2_next_r,
        a_g,
        b_g,
        c_g,
        z1_g,
        z2_g,
        h1_g,
        f_g,
        h2_g,
        transcript,
        r_inv,
    )?;

    Ok(AggregateProof {
        n: n as u64,
        prim_proof,
        gipa_proof,
    })
}

fn prove_gipa<E, T>(
    ks: &KZGCommitSRS<E>,
    mut n: usize,
    mut vk: VKey<E>,
    mut wk: WKey<E>,
    mut r: Vec<E::Fr>,
    mut a: Vec<E::Fr>,
    mut b: Vec<E::Fr>,
    mut z1_next: Vec<E::Fr>,
    mut h1_next: Vec<E::Fr>,
    mut h2: Vec<E::Fr>,
    mut a_r: Vec<E::Fr>,
    mut b_r: Vec<E::Fr>,
    mut c_r: Vec<E::Fr>,
    mut ab_r: Vec<E::Fr>,
    mut ac_r: Vec<E::Fr>,
    mut bc_r: Vec<E::Fr>,
    mut abc_r: Vec<E::Fr>,
    mut f_r: Vec<E::Fr>,
    mut z2_next_r: Vec<E::Fr>,
    mut h1_next_z2_next_r: Vec<E::Fr>,
    mut h2_z2_next_r: Vec<E::Fr>,
    mut a_g: Vec<E::G1Affine>,
    mut b_g: Vec<E::G1Affine>,
    mut c_g: Vec<E::G1Affine>,
    mut z1_g: Vec<E::G1Affine>,
    mut z2_g: Vec<E::G1Affine>,
    mut h1_g: Vec<E::G1Affine>,
    mut f_g: Vec<E::G1Affine>,
    mut h2_g: Vec<E::G1Affine>,
    transcript: &mut T,
    r_inv: E::Fr,
) -> Result<GipaProof<E>, Error>
where
    E: PairingEngine,
    T: TranscriptProtocol<E>,
{
    // storing the values for including in the proof
    let mut wire_comms = Vec::new();
    let mut perm_comms = Vec::new();
    let mut lookup_comms = Vec::new();
    let mut challenges: Vec<E::Fr> = Vec::new();
    let mut challenges_inv: Vec<E::Fr> = Vec::new();

    // GIPA
    while n > 1 {
        // recursive step
        let split = n / 2;
        
        let wire_comm = wire_gipa_commit(
            split,
            &r,
            &a,
            &b,
            &a_r,
            &b_r,
            &c_r,
            &ab_r,
            &ac_r,
            &bc_r,
            &abc_r,
            &a_g,
            &b_g,
            &c_g,
            &vk,
            &wk,
        );

        let perm_comm = permutation_gipa_commit(
            split,
            &r,
            &a_r,
            &b_r,
            &c_r,
            &ab_r,
            &ac_r,
            &bc_r,
            &abc_r,
            &z1_next,
            &z1_g,
            &vk,
        );

        let lookup_comm = lookup_gipa_commit(
            split,
            &r,
            &h1_next,
            &h2,
            &f_r,
            &z2_next_r,
            &h1_next_z2_next_r,
            &h2_z2_next_r,
            &z2_g,
            &h1_g,
            &f_g,
            &h2_g,
            &vk,
            &wk,
        );

        // Fiat-Shamir challenge
        wire_comm.0.append_in_transcript(transcript);
        wire_comm.1.append_in_transcript(transcript);
        perm_comm.0.append_in_transcript(transcript);
        perm_comm.1.append_in_transcript(transcript);
        lookup_comm.0.append_in_transcript(transcript);
        lookup_comm.1.append_in_transcript(transcript);
        let x = transcript.challenge_fr("x");
        let x_inv = x.inverse().unwrap();

        // Compress and update vectors for next recursion
        compress_scalars(&mut a_r, split, &x_inv);
        compress_scalars(&mut b_r, split, &x_inv);
        compress_scalars(&mut c_r, split, &x_inv);
        compress_scalars(&mut ab_r, split, &x_inv);
        compress_scalars(&mut ac_r, split, &x_inv);
        compress_scalars(&mut bc_r, split, &x_inv);
        compress_scalars(&mut abc_r, split, &x_inv);
        compress_scalars(&mut f_r, split, &x_inv);
        compress_scalars(&mut z2_next_r, split, &x_inv);
        compress_scalars(&mut h1_next_z2_next_r, split, &x_inv);
        compress_scalars(&mut h2_z2_next_r, split, &x_inv);

        compress_scalars(&mut r, split, &x);
        compress_scalars(&mut a, split, &x);
        compress_scalars(&mut b, split, &x);
        compress_affines(&mut a_g, split, &x);
        compress_affines(&mut b_g, split, &x);
        compress_affines(&mut c_g, split, &x);
        compress_scalars(&mut z1_next, split, &x);
        compress_scalars(&mut h1_next, split, &x);
        compress_scalars(&mut h2, split, &x);
        compress_affines(&mut z1_g, split, &x);
        compress_affines(&mut z2_g, split, &x);
        compress_affines(&mut h1_g, split, &x);
        compress_affines(&mut f_g, split, &x);
        compress_affines(&mut h2_g, split, &x);

        vk.compress(split, &x_inv);
        wk.compress(split, &x);

        wire_comms.push(wire_comm);
        perm_comms.push(perm_comm);
        lookup_comms.push(lookup_comm);
        challenges.push(x);
        challenges_inv.push(x_inv);

        n = split;
    }

    let final_value = GipaFinalValue {
        final_a_r: a_r[0],
        final_b_r: b_r[0],
        final_c_r: c_r[0],
        final_ab_r: ab_r[0],
        final_bc_r: bc_r[0],
        final_ac_r: ac_r[0],
        final_abc_r: abc_r[0],
        final_f_r: f_r[0],
        final_z2_next_r: z2_next_r[0],
        final_h1_next_z2_next_r: h1_next_z2_next_r[0],
        final_h2_z2_next_r: h2_z2_next_r[0],
        final_r: r[0],
        final_a: a[0],
        final_b: b[0],
        final_a_g: a_g[0],
        final_b_g: b_g[0],
        final_c_g: c_g[0],
        final_z1_next: z1_next[0],
        final_h1_next: h1_next[0],
        final_h2: h2[0],
        final_z1_g: z1_g[0],
        final_z2_g: z2_g[0],
        final_h1_g: h1_g[0],
        final_f_g: f_g[0],
        final_h2_g: h2_g[0],
        final_vk: vk[0],
        final_wk: wk[0],
    };
    
    // Fiat-Shamir challenge
    final_value.append_in_transcript(transcript);
    // KZG challenge point
    let z = transcript.challenge_fr("z");
    
    // Prove final commitment keys are wellformed
    // we reverse the transcript so the polynomial in kzg opening is constructed
    // correctly - the formula indicates x_{l-j}. Also for deriving KZG
    // challenge point, input must be the last challenge.
    challenges.reverse();
    challenges_inv.reverse();

    let vk_opening = prove_commitment_key(
        &ks.h_powers,
        &challenges_inv,
        z,
        E::Fr::one(),
    )?;

    let wk_opening = prove_commitment_key(
        &ks.g_powers,
        &challenges,
        z,
        r_inv,
    )?;

    Ok(GipaProof {
        wire_comms,
        perm_comms,
        lookup_comms,
        final_value,
        vk_opening,
        wk_opening,
    })
}

fn prove_commitment_key<G: AffineCurve>(
    powers: &[G],
    challenges: &[G::ScalarField],
    point: G::ScalarField,
    shift: G::ScalarField,
) -> Result<KZGOpening<G>, Error> {
    let v_coeffs = poly_coeffs_from_transcript(challenges, shift);
    let v_poly = DensePolynomial::from_coefficients_vec(v_coeffs);

    create_kzg_opening(powers, &v_poly, point)
}

#[inline]
fn scalar_mul<F: Field>(x: &[F], y: &[F]) -> Vec<F> {
    x.iter().zip(y.iter()).map(|(x, y)| *x * y).collect()
}

#[inline]
fn cover_scalar_mul<F: Field>(x: Vec<F>, y: &[F]) -> Vec<F> {
    x.into_iter().zip(y.iter()).map(|(x, y)| x * y).collect()
}
