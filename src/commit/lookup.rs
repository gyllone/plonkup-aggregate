use ark_ec::PairingEngine;

use crate::transcript::{TranscriptProtocol, TranscriptAppender};
use super::*;

pub struct LookupPrimitiveCommitments<E: PairingEngine> {
    // commit with vk
    pub h1_next_vk_comm: E::G2Affine,
    pub h2_vk_comm: E::G2Affine,
    pub z2_g_vk_comm: E::Fqk,
    pub h1_g_vk_comm: E::Fqk,
    pub f_g_vk_comm: E::Fqk,
    pub h2_g_vk_comm: E::Fqk,
    // commit with wk
    pub f_wk_comm: E::G1Affine,
    pub z2_next_wk_comm: E::G1Affine,
    pub h1_next_z2_next_wk_comm: E::G1Affine,
    pub h2_z2_next_wk_comm: E::G1Affine,
}

impl<E, T> TranscriptAppender<E, T> for LookupPrimitiveCommitments<E>
where
    E: PairingEngine,
    T: TranscriptProtocol<E>,
{
    fn append_in_transcript(&self, transcript: &mut T) {
        transcript.append_g2("h1_next_vk_comm", &self.h1_next_vk_comm);
        transcript.append_g2("h2_vk_comm", &self.h2_vk_comm);
        transcript.append_gt("z2_g_vk_comm", &self.z2_g_vk_comm);
        transcript.append_gt("h1_g_vk_comm", &self.h1_g_vk_comm);
        transcript.append_gt("f_g_vk_comm", &self.f_g_vk_comm);
        transcript.append_gt("h2_g_vk_comm", &self.h2_g_vk_comm);
        transcript.append_g1("f_wk_comm", &self.f_wk_comm);
        transcript.append_g1("z2_next_wk_comm", &self.z2_next_wk_comm);
        transcript.append_g1("h1_next_z2_next_wk_comm", &self.h1_next_z2_next_wk_comm);
        transcript.append_g1("h2_z2_next_wk_comm", &self.h2_z2_next_wk_comm);
    }
}

pub struct LookupPrimitiveInnerProducts<E: PairingEngine> {
    pub z2_g_r: E::G1Affine,
    pub z2_g_f_r: E::G1Affine,
    pub h1_g_z2_next_r: E::G1Affine,
    pub h1_g_h1_next_z2_next_r: E::G1Affine,
    pub h1_g_h2_z2_next_r: E::G1Affine,
    pub f_r: E::Fr,
    pub z2_next_r: E::Fr,
    pub h1_next_r: E::Fr,
    pub h2_r: E::Fr,
    pub h1_next_z2_next_r: E::Fr,
    pub h2_z2_next_r: E::Fr,
    pub f_g_r: E::G1Affine,
    pub h2_g_r: E::G1Affine,
}

impl<E, T> TranscriptAppender<E, T> for LookupPrimitiveInnerProducts<E>
where
    E: PairingEngine,
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

pub struct LookupGipaCommitments<E: PairingEngine> {
    // commit with vk
    pub h1_next_vk_comm: E::G2Affine,
    pub h2_vk_comm: E::G2Affine,
    pub z2_g_vk_comm: E::Fqk,
    pub h1_g_vk_comm: E::Fqk,
    pub f_g_vk_comm: E::Fqk,
    pub h2_g_vk_comm: E::Fqk,
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
    pub f_r_ip: E::Fr,
    pub z2_next_r_ip: E::Fr,
    pub h1_next_r_ip: E::Fr,
    pub h2_r_ip: E::Fr,
    pub h1_next_z2_next_r_ip: E::Fr,
    pub h2_z2_next_r_ip: E::Fr,
    pub f_g_r_ip: E::G1Affine,
    pub h2_g_r_ip: E::G1Affine,
}

impl<E, T> TranscriptAppender<E, T> for LookupGipaCommitments<E>
where
    E: PairingEngine,
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

pub(crate) fn lookup_primitive_commit<E: PairingEngine>(
    h1_next: &[E::Fr],
    h2: &[E::Fr],
    f: &[E::Fr],
    z2_next: &[E::Fr],
    h1_next_z2_next: &[E::Fr],
    h2_z2_next: &[E::Fr],
    z2_g: &[E::G1Affine],
    h1_g: &[E::G1Affine],
    f_g: &[E::G1Affine],
    h2_g: &[E::G1Affine],
    vk: &VKey<E>,
    wk: &WKey<E>,
) -> LookupPrimitiveCommitments<E> {
    // commit with vk
    let h1_next_vk_comm = multiexponentiation(h1_next, vk);
    let h2_vk_comm = multiexponentiation(h2, vk);
    let z2_g_vk_comm = pairing::<E>(z2_g, vk);
    let h1_g_vk_comm = pairing::<E>(h1_g, vk);
    let f_g_vk_comm = pairing::<E>(f_g, vk);
    let h2_g_vk_comm = pairing::<E>(h2_g, vk);

    // commit with wk
    let f_wk_comm = multiexponentiation(f, wk);
    let z2_next_wk_comm = multiexponentiation(z2_next, wk);
    let h1_next_z2_next_wk_comm = multiexponentiation(h1_next_z2_next, wk);
    let h2_z2_next_wk_comm = multiexponentiation(h2_z2_next, wk);

    LookupPrimitiveCommitments {
        h1_next_vk_comm,
        h2_vk_comm,
        z2_g_vk_comm,
        h1_g_vk_comm,
        f_g_vk_comm,
        h2_g_vk_comm,
        f_wk_comm,
        z2_next_wk_comm,
        h1_next_z2_next_wk_comm,
        h2_z2_next_wk_comm,
    }
}

pub(crate) fn lookup_primitive_inner_product<E: PairingEngine>(
    r: &[E::Fr],
    h1_next: &[E::Fr],
    h2: &[E::Fr],
    f_r: &[E::Fr], // scaled f(z) * r^i
    z2_next_r: &[E::Fr], // scaled z2(ωz) * r^i
    h1_next_z2_next_r: &[E::Fr], // scaled h1(ωz)z2(ωz) * r^i
    h2_z2_next_r: &[E::Fr], // scaled h2(z)z2(ωz) * r^i
    z2_g: &[E::G1Affine],
    h1_g: &[E::G1Affine],
    f_g: &[E::G1Affine],
    h2_g: &[E::G1Affine],
) -> LookupPrimitiveInnerProducts<E> {
    let z2_g_r = multiexponentiation(r, z2_g);
    let z2_g_f_r = multiexponentiation(f_r, z2_g);
    let h1_g_z2_next_r = multiexponentiation(z2_next_r, h1_g);
    let h1_g_h1_next_z2_next_r = multiexponentiation(h1_next_z2_next_r, h1_g);
    let h1_g_h2_z2_next_r = multiexponentiation(h2_z2_next_r, h1_g);
    let f_r = f_r.iter().sum();
    let z2_next_r = z2_next_r.iter().sum();
    let h1_next_r = scalars_inner_product(h1_next, r);
    let h2_r = scalars_inner_product(h2, r);
    let h1_next_z2_next_r = h1_next_z2_next_r.iter().sum();
    let h2_z2_next_r = h2_z2_next_r.iter().sum();
    let f_g_r = multiexponentiation(r, f_g);
    let h2_g_r = multiexponentiation(r, h2_g);

    LookupPrimitiveInnerProducts {
        z2_g_r,
        z2_g_f_r,
        h1_g_z2_next_r,
        h1_g_h1_next_z2_next_r,
        h1_g_h2_z2_next_r,
        f_r,
        z2_next_r,
        h1_next_r,
        h2_r,
        h1_next_z2_next_r,
        h2_z2_next_r,
        f_g_r,
        h2_g_r,
    }
}

pub(crate) fn lookup_gipa_commit<E: PairingEngine>(
    split: usize,
    r: &[E::Fr],
    h1_next: &[E::Fr],
    h2: &[E::Fr],
    f_r: &[E::Fr], // scaled f(z) * r^i
    z2_next_r: &[E::Fr], // scaled z2(ωz) * r^i
    h1_next_z2_next_r: &[E::Fr], // scaled h1(ωz)z2(ωz) * r^i
    h2_z2_next_r: &[E::Fr], // scaled h2(z)z2(ωz) * r^i
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
