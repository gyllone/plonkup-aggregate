use std::ops::AddAssign;
use ark_ec::pairing::{Pairing, PairingOutput};

use crate::{transcript::{TranscriptProtocol, TranscriptAppender}, proof::GipaFinalValue};
use super::*;

pub struct LookupPrimitiveCommitments<E: Pairing> {
    // commit with vk
    pub h1_next_vk: E::G2Affine,
    pub h2_vk: E::G2Affine,
    pub z2_g_vk: PairingOutput<E>,
    pub h1_g_vk: PairingOutput<E>,
    pub f_g_vk: PairingOutput<E>,
    pub h2_g_vk: PairingOutput<E>,
    // commit with wk
    pub f_wk: E::G1Affine,
    pub z2_next_wk: E::G1Affine,
    pub h1_next_z2_next_wk: E::G1Affine,
    pub h2_z2_next_wk: E::G1Affine,
}

impl<E, T> TranscriptAppender<E, T> for LookupPrimitiveCommitments<E>
where
    E: Pairing,
    T: TranscriptProtocol<E>,
{
    fn append_in_transcript(&self, transcript: &mut T) {
        transcript.append_g2("h1_next_vk", &self.h1_next_vk);
        transcript.append_g2("h2_vk", &self.h2_vk);
        transcript.append_gt("z2_g_vk", &self.z2_g_vk);
        transcript.append_gt("h1_g_vk", &self.h1_g_vk);
        transcript.append_gt("f_g_vk", &self.f_g_vk);
        transcript.append_gt("h2_g_vk", &self.h2_g_vk);
        transcript.append_g1("f_wk", &self.f_wk);
        transcript.append_g1("z2_next_wk", &self.z2_next_wk);
        transcript.append_g1("h1_next_z2_next_wk", &self.h1_next_z2_next_wk);
        transcript.append_g1("h2_z2_next_wk", &self.h2_z2_next_wk);
    }
}

pub(crate) fn lookup_primitive_commit<E: Pairing>(
    h1_next: &[E::ScalarField],
    h2: &[E::ScalarField],
    f: &[E::ScalarField],
    z2_next: &[E::ScalarField],
    h1_next_z2_next: &[E::ScalarField],
    h2_z2_next: &[E::ScalarField],
    z2_g: &[E::G1Affine],
    h1_g: &[E::G1Affine],
    f_g: &[E::G1Affine],
    h2_g: &[E::G1Affine],
    vk: &VKey<E>,
    wk: &WKey<E>,
) -> LookupPrimitiveCommitments<E> {
    LookupPrimitiveCommitments {
        h1_next_vk: multiexponentiation(h1_next, vk),
        h2_vk: multiexponentiation(h2, vk),
        z2_g_vk: pairing::<E>(z2_g, vk),
        h1_g_vk: pairing::<E>(h1_g, vk),
        f_g_vk: pairing::<E>(f_g, vk),
        h2_g_vk: pairing::<E>(h2_g, vk),
        f_wk: multiexponentiation(f, wk),
        z2_next_wk: multiexponentiation(z2_next, wk),
        h1_next_z2_next_wk: multiexponentiation(h1_next_z2_next, wk),
        h2_z2_next_wk: multiexponentiation(h2_z2_next, wk),
    }
}

pub struct LookupPrimitiveInnerProducts<E: Pairing> {
    pub z2_g_r: E::G1Affine,
    pub z2_g_f_r: E::G1Affine,
    pub h1_g_z2_next_r: E::G1Affine,
    pub h1_g_h1_next_z2_next_r: E::G1Affine,
    pub h1_g_h2_z2_next_r: E::G1Affine,
    pub f_r: E::ScalarField,
    pub z2_next_r: E::ScalarField,
    pub h1_next_r: E::ScalarField,
    pub h2_r: E::ScalarField,
    pub h1_next_z2_next_r: E::ScalarField,
    pub h2_z2_next_r: E::ScalarField,
    pub f_g_r: E::G1Affine,
    pub h2_g_r: E::G1Affine,
}

impl<E, T> TranscriptAppender<E, T> for LookupPrimitiveInnerProducts<E>
where
    E: Pairing,
    T: TranscriptProtocol<E>,
{
    fn append_in_transcript(&self, transcript: &mut T) {
        transcript.append_g1("z2_g_r", &self.z2_g_r);
        transcript.append_g1("z2_g_f_r", &self.z2_g_f_r);
        transcript.append_g1("h1_g_z2_next_r", &self.h1_g_z2_next_r);
        transcript.append_g1("h1_g_h1_next_z2_next_r", &self.h1_g_h1_next_z2_next_r);
        transcript.append_g1("h1_g_h2_z2_next_r", &self.h1_g_h2_z2_next_r);
        transcript.append_fr("f_r", &self.f_r);
        transcript.append_fr("z2_next_r", &self.z2_next_r);
        transcript.append_fr("h1_next_r", &self.h1_next_r);
        transcript.append_fr("h2_r", &self.h2_r);
        transcript.append_fr("h1_next_z2_next_r", &self.h1_next_z2_next_r);
        transcript.append_fr("h2_z2_next_r", &self.h2_z2_next_r);
        transcript.append_g1("f_g_r", &self.f_g_r);
        transcript.append_g1("h2_g_r", &self.h2_g_r);
    }
}

pub(crate) fn lookup_primitive_inner_product<E: Pairing>(
    r: &[E::ScalarField],
    h1_next: &[E::ScalarField],
    h2: &[E::ScalarField],
    f_r: &[E::ScalarField], // scaled f(z) * r^i
    z2_next_r: &[E::ScalarField], // scaled z2(ωz) * r^i
    h1_next_z2_next_r: &[E::ScalarField], // scaled h1(ωz)z2(ωz) * r^i
    h2_z2_next_r: &[E::ScalarField], // scaled h2(z)z2(ωz) * r^i
    z2_g: &[E::G1Affine],
    h1_g: &[E::G1Affine],
    f_g: &[E::G1Affine],
    h2_g: &[E::G1Affine],
) -> LookupPrimitiveInnerProducts<E> {
    LookupPrimitiveInnerProducts {
        z2_g_r: multiexponentiation(r, z2_g),
        z2_g_f_r: multiexponentiation(f_r, z2_g),
        h1_g_z2_next_r: multiexponentiation(z2_next_r, h1_g),
        h1_g_h1_next_z2_next_r: multiexponentiation(h1_next_z2_next_r, h1_g),
        h1_g_h2_z2_next_r: multiexponentiation(h2_z2_next_r, h1_g),
        f_r: f_r.iter().sum(),
        z2_next_r: z2_next_r.iter().sum(),
        h1_next_r: scalars_inner_product(h1_next, r),
        h2_r: scalars_inner_product(h2, r),
        h1_next_z2_next_r: h1_next_z2_next_r.iter().sum(),
        h2_z2_next_r: h2_z2_next_r.iter().sum(),
        f_g_r: multiexponentiation(r, f_g),
        h2_g_r: multiexponentiation(r, h2_g),
    }
}

pub struct LookupGipaCommitments<E: Pairing> {
    // commit with vk
    pub h1_next_vk_comm: E::G2Affine,
    pub h2_vk_comm: E::G2Affine,
    pub z2_g_vk_comm: PairingOutput<E>,
    pub h1_g_vk_comm: PairingOutput<E>,
    pub f_g_vk_comm: PairingOutput<E>,
    pub h2_g_vk_comm: PairingOutput<E>,
    // commit with wk
    pub f_r_wk_comm: E::G1Affine,
    pub z2_next_r_wk_comm: E::G1Affine,
    pub h1_next_z2_next_r_wk_comm: E::G1Affine,
    pub h2_z2_next_r_wk_comm: E::G1Affine,
    // inner product
    pub z2_g_r_ip: E::G1Affine,
    pub z2_g_f_r_ip: E::G1Affine,
    pub h1_g_z2_next_r_ip: E::G1Affine,
    pub h1_g_h1_next_z2_next_r_ip: E::G1Affine,
    pub h1_g_h2_z2_next_r_ip: E::G1Affine,
    pub f_r_ip: E::ScalarField,
    pub z2_next_r_ip: E::ScalarField,
    pub h1_next_r_ip: E::ScalarField,
    pub h2_r_ip: E::ScalarField,
    pub h1_next_z2_next_r_ip: E::ScalarField,
    pub h2_z2_next_r_ip: E::ScalarField,
    pub f_g_r_ip: E::G1Affine,
    pub h2_g_r_ip: E::G1Affine,
}

impl<E, T> TranscriptAppender<E, T> for LookupGipaCommitments<E>
where
    E: Pairing,
    T: TranscriptProtocol<E>,
{
    fn append_in_transcript(&self, transcript: &mut T) {
        transcript.append_g2("h1_next_vk_comm", &self.h1_next_vk_comm);
        transcript.append_g2("h2_vk_comm", &self.h2_vk_comm);
        transcript.append_gt("z2_g_vk_comm", &self.z2_g_vk_comm);
        transcript.append_gt("h1_g_vk_comm", &self.h1_g_vk_comm);
        transcript.append_gt("f_g_vk_comm", &self.f_g_vk_comm);
        transcript.append_gt("h2_g_vk_comm", &self.h2_g_vk_comm);
        transcript.append_g1("f_r_wk_comm", &self.f_r_wk_comm);
        transcript.append_g1("z2_next_r_wk_comm", &self.z2_next_r_wk_comm);
        transcript.append_g1("h1_next_z2_next_r_wk_comm", &self.h1_next_z2_next_r_wk_comm);
        transcript.append_g1("h2_z2_next_r_wk_comm", &self.h2_z2_next_r_wk_comm);
        transcript.append_g1("z2_g_r_ip", &self.z2_g_r_ip);
        transcript.append_g1("z2_g_f_r_ip", &self.z2_g_f_r_ip);
        transcript.append_g1("h1_g_z2_next_r_ip", &self.h1_g_z2_next_r_ip);
        transcript.append_g1("h1_g_h1_next_z2_next_r_ip", &self.h1_g_h1_next_z2_next_r_ip);
        transcript.append_g1("h1_g_h2_z2_next_r_ip", &self.h1_g_h2_z2_next_r_ip);
        transcript.append_fr("f_r_ip", &self.f_r_ip);
        transcript.append_fr("z2_next_r_ip", &self.z2_next_r_ip);
        transcript.append_fr("h1_next_r_ip", &self.h1_next_r_ip);
        transcript.append_fr("h2_r_ip", &self.h2_r_ip);
        transcript.append_fr("h1_next_z2_next_r_ip", &self.h1_next_z2_next_r_ip);
        transcript.append_fr("h2_z2_next_r_ip", &self.h2_z2_next_r_ip);
        transcript.append_g1("f_g_r_ip", &self.f_g_r_ip);
        transcript.append_g1("h2_g_r_ip", &self.h2_g_r_ip);
    }
}

pub(crate) fn lookup_gipa_commit<E: Pairing>(
    split: usize,
    r: &[E::ScalarField],
    h1_next: &[E::ScalarField],
    h2: &[E::ScalarField],
    f_r: &[E::ScalarField], // scaled f(z) * r^i
    z2_next_r: &[E::ScalarField], // scaled z2(ωz) * r^i
    h1_next_z2_next_r: &[E::ScalarField], // scaled h1(ωz)z2(ωz) * r^i
    h2_z2_next_r: &[E::ScalarField], // scaled h2(z)z2(ωz) * r^i
    z2_g: &[E::G1Affine],
    h1_g: &[E::G1Affine],
    f_g: &[E::G1Affine],
    h2_g: &[E::G1Affine],
    vk: &VKey<E>,
    wk: &WKey<E>, // scaled key w^r^-1
) -> (LookupGipaCommitments<E>, LookupGipaCommitments<E>) {
    let (r_left, r_right) = r.split_at(split);
    let (h1_next_left, h1_next_right) = h1_next.split_at(split);
    let (h2_left, h2_right) = h2.split_at(split);

    let (f_r_left, f_r_right) = f_r.split_at(split);
    let (z2_next_r_left, z2_next_r_right) =
        z2_next_r.split_at(split);
    let (h1_next_z2_next_r_left, h1_next_z2_next_r_right) =
        h1_next_z2_next_r.split_at(split);
    let (h2_z2_next_r_left, h2_z2_next_r_right) =
        h2_z2_next_r.split_at(split);

    let (z2_g_left, z2_g_right) = z2_g.split_at(split);
    let (h1_g_left, h1_g_right) = h1_g.split_at(split);
    let (h2_g_left, h2_g_right) = h2_g.split_at(split);
    let (f_g_left, f_g_right) = f_g.split_at(split);

    let (vk_left, vk_right) = vk.split_at(split);
    let (wk_left, wk_right) = wk.split_at(split);

    // Commit with vk
    // ----------------------------------------------------------------
    //   (h1_ω, vk)       |                       |
    let h1_next_vk_comm_left =
        multiexponentiation(h1_next_right, vk_left);
    let h1_next_vk_comm_right =
        multiexponentiation(h1_next_left, vk_right);
    // ----------------------------------------------------------------
    //   (h2, vk)         |                       |
    let h2_vk_comm_left =
        multiexponentiation(h2_right, vk_left);
    let h2_vk_comm_right =
        multiexponentiation(h2_left, vk_right);
    // ----------------------------------------------------------------
    //   (z2_g, vk)       |                       |
    let z2_g_vk_comm_left = pairing::<E>(z2_g_right, vk_left);
    let z2_g_vk_comm_right = pairing::<E>(z2_g_left, vk_right);
    // ----------------------------------------------------------------
    //   (h1_g, vk)       |                       |
    let h1_g_vk_comm_left = pairing::<E>(h1_g_right, vk_left);
    let h1_g_vk_comm_right = pairing::<E>(h1_g_left, vk_right);
    // ----------------------------------------------------------------
    //   (f_g, vk)        |                       |
    let f_g_vk_comm_left = pairing::<E>(f_g_right, vk_left);
    let f_g_vk_comm_right = pairing::<E>(f_g_left, vk_right);
    // ----------------------------------------------------------------
    //   (h2_g, vk)       |                       |
    let h2_g_vk_comm_left = pairing::<E>(h2_g_right, vk_left);
    let h2_g_vk_comm_right = pairing::<E>(h2_g_left, vk_right);

    // Commit with wk
    // ----------------------------------------------------------------
    //                    |       (f * r, wk)     |
    let f_r_wk_comm_left =
        multiexponentiation(f_r_left, wk_right);
    let f_r_wk_comm_right =
        multiexponentiation(f_r_right, wk_left);
    // ----------------------------------------------------------------
    //                    |     (z2_ω * r, wk)    |
    let z2_next_r_wk_comm_left =
        multiexponentiation(z2_next_r_left, wk_right);
    let z2_next_r_wk_comm_right =
        multiexponentiation(z2_next_r_right, wk_left);
    // ----------------------------------------------------------------
    //                    |  (h1_ω*z2_ω * r, wk)  |
    let h1_next_z2_next_r_wk_comm_left =
        multiexponentiation(h1_next_z2_next_r_left, wk_right);
    let h1_next_z2_next_r_wk_comm_right =
        multiexponentiation(h1_next_z2_next_r_right, wk_left);
    // ----------------------------------------------------------------
    //                    |   (h2*z2_ω * r, wk)   |
    let h2_z2_next_r_wk_comm_left =
        multiexponentiation(h2_z2_next_r_left, wk_right);
    let h2_z2_next_r_wk_comm_right =
        multiexponentiation(h2_z2_next_r_right, wk_left);

    // Inner Product
    // ----------------------------------------------------------------
    //                    |                       |  (z2_g, r)
    let z2_g_r_ip_left =
        multiexponentiation(r_left, z2_g_right);
    let z2_g_r_ip_right =
        multiexponentiation(r_right, z2_g_left);
    // ----------------------------------------------------------------
    //                    |                       |  (z2_g, f * r)
    let z2_g_f_r_ip_left =
        multiexponentiation(f_r_left, z2_g_right);
    let z2_g_f_r_ip_right =
        multiexponentiation(f_r_right, z2_g_left);
    // ----------------------------------------------------------------
    //                    |                       |  (h1_g, z2_ω * r)
    let h1_g_z2_next_r_ip_left =
        multiexponentiation(z2_next_r_left, h1_g_right);
    let h1_g_z2_next_r_ip_right =
        multiexponentiation(z2_next_r_right, h1_g_left);
    // ----------------------------------------------------------------
    //                    |                       |  (h1_g, h1_ω*z2_ω * r)
    let h1_g_h1_next_z2_next_r_ip_left =
        multiexponentiation(h1_next_z2_next_r_left, h1_g_right);
    let h1_g_h1_next_z2_next_r_ip_right =
        multiexponentiation(h1_next_z2_next_r_right, h1_g_left);
    // ----------------------------------------------------------------
    //                    |                       |  (h1_g, h2*z2_ω * r)
    let h1_g_h2_z2_next_r_ip_left =
        multiexponentiation(h2_z2_next_r_left, h1_g_right);
    let h1_g_h2_z2_next_r_ip_right =
        multiexponentiation(h2_z2_next_r_right, h1_g_left);
    // ----------------------------------------------------------------
    //                    |                       |  (1, f * r)
    let f_r_ip_left = f_r_left.iter().sum();
    let f_r_ip_right = f_r_right.iter().sum();
    // ----------------------------------------------------------------
    //                    |                       |  (1, z2_ω * r)
    let z2_next_r_ip_left = z2_next_r_left.iter().sum();
    let z2_next_r_ip_right = z2_next_r_right.iter().sum();
    // ---------------------------------------------------------------
    //                    |                       |  (h1_ω, r)
    let h1_next_r_ip_left =
        scalars_inner_product(h1_next_right, r_left);
    let h1_next_r_ip_right =
        scalars_inner_product(h1_next_left, r_right);
    // ---------------------------------------------------------------
    //                    |                       |  (h2, r)
    let h2_r_ip_left =
        scalars_inner_product(h2_right, r_left);
    let h2_r_ip_right =
        scalars_inner_product(h2_left, r_right);
    // ----------------------------------------------------------------
    //                    |                       |  (1, h1_ω*z2_ω * r) == (h1_ω, z2_ω * r)
    let h1_next_z2_next_r_ip_left = h1_next_z2_next_r_left.iter().sum();
    let h1_next_z2_next_r_ip_right = h1_next_z2_next_r_right.iter().sum();
    // ----------------------------------------------------------------
    //                    |                       |  (1, h2*z2_ω * r) == (h2, z2_ω * r)
    let h2_z2_next_r_ip_left = h2_z2_next_r_left.iter().sum();
    let h2_z2_next_r_ip_right = h2_z2_next_r_right.iter().sum();
    // ----------------------------------------------------------------
    //                    |                       |  (f_g, r)
    let f_g_r_ip_left = multiexponentiation(r_left, f_g_right);
    let f_g_r_ip_right = multiexponentiation(r_right, f_g_left);
    // ----------------------------------------------------------------
    //                    |                       |  (h2_g, r)
    let h2_g_r_ip_left = multiexponentiation(r_left, h2_g_right);
    let h2_g_r_ip_right = multiexponentiation(r_right, h2_g_left);

    let left_comms = LookupGipaCommitments {
        h1_next_vk_comm: h1_next_vk_comm_left,
        h2_vk_comm: h2_vk_comm_left,
        z2_g_vk_comm: z2_g_vk_comm_left,
        h1_g_vk_comm: h1_g_vk_comm_left,
        f_g_vk_comm: f_g_vk_comm_left,
        h2_g_vk_comm: h2_g_vk_comm_left,
        f_r_wk_comm: f_r_wk_comm_left,
        z2_next_r_wk_comm: z2_next_r_wk_comm_left, 
        h1_next_z2_next_r_wk_comm: h1_next_z2_next_r_wk_comm_left,
        h2_z2_next_r_wk_comm: h2_z2_next_r_wk_comm_left,
        z2_g_r_ip: z2_g_r_ip_left,
        z2_g_f_r_ip: z2_g_f_r_ip_left,
        h1_g_z2_next_r_ip: h1_g_z2_next_r_ip_left,
        h1_g_h1_next_z2_next_r_ip: h1_g_h1_next_z2_next_r_ip_left,
        h1_g_h2_z2_next_r_ip: h1_g_h2_z2_next_r_ip_left,
        f_r_ip: f_r_ip_left,
        z2_next_r_ip: z2_next_r_ip_left,
        h1_next_r_ip: h1_next_r_ip_left,
        h2_r_ip: h2_r_ip_left,
        h1_next_z2_next_r_ip: h1_next_z2_next_r_ip_left,
        h2_z2_next_r_ip: h2_z2_next_r_ip_left,
        f_g_r_ip: f_g_r_ip_left,
        h2_g_r_ip: h2_g_r_ip_left,
    };
    let right_comms = LookupGipaCommitments {
        h1_next_vk_comm: h1_next_vk_comm_right,
        h2_vk_comm: h2_vk_comm_right,
        z2_g_vk_comm: z2_g_vk_comm_right,
        h1_g_vk_comm: h1_g_vk_comm_right,
        f_g_vk_comm: f_g_vk_comm_right,
        h2_g_vk_comm: h2_g_vk_comm_right,
        f_r_wk_comm: f_r_wk_comm_right,
        z2_next_r_wk_comm: z2_next_r_wk_comm_right, 
        h1_next_z2_next_r_wk_comm: h1_next_z2_next_r_wk_comm_right,
        h2_z2_next_r_wk_comm: h2_z2_next_r_wk_comm_right,
        z2_g_r_ip: z2_g_r_ip_right,
        z2_g_f_r_ip: z2_g_f_r_ip_right,
        h1_g_z2_next_r_ip: h1_g_z2_next_r_ip_right,
        h1_g_h1_next_z2_next_r_ip: h1_g_h1_next_z2_next_r_ip_right,
        h1_g_h2_z2_next_r_ip: h1_g_h2_z2_next_r_ip_right,
        f_r_ip: f_r_ip_right,
        z2_next_r_ip: z2_next_r_ip_right,
        h1_next_r_ip: h1_next_r_ip_right,
        h2_r_ip: h2_r_ip_right,
        h1_next_z2_next_r_ip: h1_next_z2_next_r_ip_right,
        h2_z2_next_r_ip: h2_z2_next_r_ip_right,
        f_g_r_ip: f_g_r_ip_right,
        h2_g_r_ip: h2_g_r_ip_right,
    };

    (left_comms, right_comms)
}

pub struct LookupFinalCommitments<E: Pairing> {
    // commit with vk
    h1_next_vk_comm: E::G2,
    h2_vk_comm: E::G2,
    z2_g_vk_comm: PairingOutput<E>,
    h1_g_vk_comm: PairingOutput<E>,
    f_g_vk_comm: PairingOutput<E>,
    h2_g_vk_comm: PairingOutput<E>,
    // commit with wk
    f_r_wk_comm: E::G1,
    z2_next_r_wk_comm: E::G1,
    h1_next_z2_next_r_wk_comm: E::G1,
    h2_z2_next_r_wk_comm: E::G1,
    // inner product
    z2_g_r_ip: E::G1,
    z2_g_f_r_ip: E::G1,
    h1_g_z2_next_r_ip: E::G1,
    h1_g_h1_next_z2_next_r_ip: E::G1,
    h1_g_h2_z2_next_r_ip: E::G1,
    f_r_ip: E::ScalarField,
    z2_next_r_ip: E::ScalarField,
    h1_next_r_ip: E::ScalarField,
    h2_r_ip: E::ScalarField,
    h1_next_z2_next_r_ip: E::ScalarField,
    h2_z2_next_r_ip: E::ScalarField,
    f_g_r_ip: E::G1,
    h2_g_r_ip: E::G1,
}

impl<E: Pairing> LookupFinalCommitments<E> {
    pub(crate) fn from_primitive(
        comm: LookupPrimitiveCommitments<E>,
        ip: LookupPrimitiveInnerProducts<E>
    ) -> Self {
        Self {
            h1_next_vk_comm: comm.h1_next_vk.into(),
            h2_vk_comm: comm.h2_vk.into(),
            z2_g_vk_comm: comm.z2_g_vk,
            h1_g_vk_comm: comm.h1_g_vk,
            f_g_vk_comm: comm.f_g_vk,
            h2_g_vk_comm: comm.h2_g_vk,
            f_r_wk_comm: comm.f_wk.into(),
            z2_next_r_wk_comm: comm.z2_next_wk.into(),
            h1_next_z2_next_r_wk_comm: comm.h1_next_z2_next_wk.into(),
            h2_z2_next_r_wk_comm: comm.h2_z2_next_wk.into(),
            z2_g_r_ip: ip.z2_g_r.into(),
            z2_g_f_r_ip: ip.z2_g_f_r.into(),
            h1_g_z2_next_r_ip: ip.h1_g_z2_next_r.into(),
            h1_g_h1_next_z2_next_r_ip: ip.h1_g_h1_next_z2_next_r.into(),
            h1_g_h2_z2_next_r_ip: ip.h1_g_h2_z2_next_r.into(),
            f_r_ip: ip.f_r,
            z2_next_r_ip: ip.z2_next_r,
            h1_next_r_ip: ip.h1_next_r,
            h2_r_ip: ip.h2_r,
            h1_next_z2_next_r_ip: ip.h1_next_z2_next_r,
            h2_z2_next_r_ip: ip.h2_z2_next_r,
            f_g_r_ip: ip.f_g_r.into(),
            h2_g_r_ip: ip.h2_g_r.into(),
        }
    }

    pub(crate) fn merge(
        &mut self,
        left: LookupGipaCommitments<E>,
        right: LookupGipaCommitments<E>,
        x: &E::ScalarField,
        x_inv: &E::ScalarField,
    ) {
        self.h1_next_vk_comm.add_assign(left.h1_next_vk_comm * x + (right.h1_next_vk_comm * x_inv));
        self.h2_vk_comm.add_assign(left.h2_vk_comm * x + (right.h2_vk_comm * x_inv));
        self.z2_g_vk_comm.add_assign(left.z2_g_vk_comm * x + (right.z2_g_vk_comm * x_inv));
        self.h1_g_vk_comm.add_assign(left.h1_g_vk_comm * x + (right.h1_g_vk_comm * x_inv));
        self.f_g_vk_comm.add_assign(left.f_g_vk_comm * x + (right.f_g_vk_comm * x_inv));
        self.h2_g_vk_comm.add_assign(left.h2_g_vk_comm * x + (right.h2_g_vk_comm * x_inv));
        self.f_r_wk_comm.add_assign(left.f_r_wk_comm * x + (right.f_r_wk_comm * x_inv));
        self.z2_next_r_wk_comm.add_assign(left.z2_next_r_wk_comm * x + (right.z2_next_r_wk_comm * x_inv));
        self.h1_next_z2_next_r_wk_comm.add_assign(
            left.h1_next_z2_next_r_wk_comm * x + (right.h1_next_z2_next_r_wk_comm * x_inv),
        );
        self.h2_z2_next_r_wk_comm.add_assign(
            left.h2_z2_next_r_wk_comm * x + (right.h2_z2_next_r_wk_comm * x_inv),
        );
        self.z2_g_r_ip.add_assign(left.z2_g_r_ip * x + (right.z2_g_r_ip * x_inv));
        self.z2_g_f_r_ip.add_assign(left.z2_g_f_r_ip * x + (right.z2_g_f_r_ip * x_inv));
        self.h1_g_z2_next_r_ip.add_assign(
            left.h1_g_z2_next_r_ip * x + (right.h1_g_z2_next_r_ip * x_inv),
        );
        self.h1_g_h1_next_z2_next_r_ip.add_assign(
            left.h1_g_h1_next_z2_next_r_ip * x + (right.h1_g_h1_next_z2_next_r_ip * x_inv),
        );
        self.h1_g_h2_z2_next_r_ip.add_assign(
            left.h1_g_h2_z2_next_r_ip * x + (right.h1_g_h2_z2_next_r_ip * x_inv),
        );
        self.f_r_ip.add_assign(left.f_r_ip * x + (right.f_r_ip * x_inv));
        self.z2_next_r_ip.add_assign(left.z2_next_r_ip * x + (right.z2_next_r_ip * x_inv));
        self.h1_next_r_ip.add_assign(left.h1_next_r_ip * x + (right.h1_next_r_ip * x_inv));
        self.h2_r_ip.add_assign(left.h2_r_ip * x + (right.h2_r_ip * x_inv));
        self.h1_next_z2_next_r_ip.add_assign(
            left.h1_next_z2_next_r_ip * x + (right.h1_next_z2_next_r_ip * x_inv),
        );
        self.h2_z2_next_r_ip.add_assign(
            left.h2_z2_next_r_ip * x + (right.h2_z2_next_r_ip * x_inv),
        );
        self.f_g_r_ip.add_assign(left.f_g_r_ip * x + (right.f_g_r_ip * x_inv));
        self.h2_g_r_ip.add_assign(left.h2_g_r_ip * x + (right.h2_g_r_ip * x_inv));
    }

    pub(crate) fn check(
        self,
        fv: &GipaFinalValue<E>,
        f1: &E::ScalarField,
        fr: &E::ScalarField,
    ) -> bool {
        // Commit with vk
        // (h1_ω, vk)
        if fv.vk * fv.h1_next != self.h1_next_vk_comm {
        }
        // (h2, vk)
        if fv.vk * fv.h2 != self.h2_vk_comm {
        }
        // (z2_g, vk)
        if E::pairing(fv.z2_g, fv.vk) != self.z2_g_vk_comm {
        }
        // (h1_g, vk)
        if E::pairing(fv.h1_g, fv.vk) != self.h1_g_vk_comm {
        }
        // (f_g, vk)
        if E::pairing(fv.f_g, fv.vk) != self.f_g_vk_comm {
        }
        // (h2_g, vk)
        if E::pairing(fv.h2_g, fv.vk) != self.h2_g_vk_comm {
        }

        // Commit with wk
        // (f * r, wk)
        if fv.wk * fv.f_r == self.f_r_wk_comm {
        }
        // (z2_ω * r, wk)
        if fv.wk * fv.z2_next_r == self.z2_next_r_wk_comm {
        }
        // (h1_ω*z2_ω * r, wk)
        if fv.wk * fv.h1_next_z2_next_r == self.h1_next_z2_next_r_wk_comm {
        }
        // (h2*z2_ω * r, wk)
        if fv.wk * fv.h2_z2_next_r == self.h2_z2_next_r_wk_comm {
        }

        // Inner product
        // (z2_g, r)
        if fv.z2_g * fr == self.z2_g_r_ip {
        }
        // (z2_g, f * r)
        if fv.z2_g * fr == self.z2_g_f_r_ip {
        }
        // (h1_g, z2_ω * r)
        if fv.h1_g * fv.z2_next_r == self.h1_g_z2_next_r_ip {
        }
        // (h1_g, h1_ω*z2_ω * r)
        if fv.h1_g * fv.h1_next_z2_next_r == self.h1_g_h1_next_z2_next_r_ip {
        }
        // (h1_g, h2*z2_ω * r)
        if fv.h1_g * fv.h2_z2_next_r == self.h1_g_h2_z2_next_r_ip {
        }
        // (1, f * r)
        if fv.f_r * f1 == self.f_r_ip {
        }
        // (1, z2_ω * r)
        if fv.z2_next_r * f1 == self.z2_next_r_ip {
        }
        // (h1_ω, r)
        if fv.h1_next * fr == self.h1_next_r_ip {
        }
        // (h2, r)
        if fv.h2 * fr == self.h2_r_ip {
        }
        // (1, h1_ω*z2_ω * r) == (h1_ω, z2_ω * r)
        if fv.h1_next_z2_next_r * f1 == self.h1_next_z2_next_r_ip
            || fv.h1_next * fv.z2_next_r == self.h1_next_z2_next_r_ip {
        }
        // (1, h2*z2_ω * r) == (h2, z2_ω * r)
        if fv.h2_z2_next_r * f1 == self.h2_z2_next_r_ip
            || fv.h2 * fv.z2_next_r == self.h2_z2_next_r_ip {
        }
        // (f_g, r)
        if fv.f_g * fr == self.f_g_r_ip {
        }
        // (h2_g, r)
        if fv.h2_g * fr == self.h2_g_r_ip {
        }

        true
    }
}
