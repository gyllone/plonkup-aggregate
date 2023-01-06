use ark_ec::PairingEngine;

use crate::transcript::{TranscriptProtocol, TranscriptAppender};
use super::*;

pub struct PermPrimitiveCommitments<E: PairingEngine> {
    // commit with vk
    pub z1_g_vk_comm: E::Fqk,
    pub z1_next_vk_comm: E::G2Affine,
}

impl<E, T> TranscriptAppender<E, T> for PermPrimitiveCommitments<E>
where
    E: PairingEngine,
    T: TranscriptProtocol<E>,
{
    fn append_in_transcript(&self, transcript: &mut T) {
        transcript.append_gt("z1_g_vk_comm", &self.z1_g_vk_comm);
        transcript.append_g2("z1_next_vk_comm", &self.z1_next_vk_comm);
    }
}

pub struct PermPrimitiveInnerProducts<E: PairingEngine> {
    pub z1_g_r: E::G1Affine,
    pub z1_g_a_r: E::G1Affine,
    pub z1_g_b_r: E::G1Affine,
    pub z1_g_c_r: E::G1Affine,
    pub z1_g_ab_r: E::G1Affine,
    pub z1_g_ac_r: E::G1Affine,
    pub z1_g_bc_r: E::G1Affine,
    pub z1_g_abc_r: E::G1Affine,
    pub z1_next_r: E::Fr,
    pub z1_next_a_r: E::Fr,
    pub z1_next_b_r: E::Fr,
    pub z1_next_c_r: E::Fr,
    pub z1_next_ab_r: E::Fr,
    pub z1_next_ac_r: E::Fr,
    pub z1_next_bc_r: E::Fr,
    pub z1_next_abc_r: E::Fr,
}

impl<E, T> TranscriptAppender<E, T> for PermPrimitiveInnerProducts<E>
where
    E: PairingEngine,
    T: TranscriptProtocol<E>,
{
    fn append_in_transcript(&self, transcript: &mut T) {
        transcript.append_g1("z1_g_r", &self.z1_g_r);
        transcript.append_g1("z1_g_a_r", &self.z1_g_a_r);
        transcript.append_g1("z1_g_b_r", &self.z1_g_b_r);
        transcript.append_g1("z1_g_c_r", &self.z1_g_c_r);
        transcript.append_g1("z1_g_ab_r", &self.z1_g_ab_r);
        transcript.append_g1("z1_g_ac_r", &self.z1_g_ac_r);
        transcript.append_g1("z1_g_bc_r", &self.z1_g_bc_r);
        transcript.append_g1("z1_g_abc_r", &self.z1_g_abc_r);
        transcript.append_fr("z1_next_r", &self.z1_next_r);
        transcript.append_fr("z1_next_a_r", &self.z1_next_a_r);
        transcript.append_fr("z1_next_b_r", &self.z1_next_b_r);
        transcript.append_fr("z1_next_c_r", &self.z1_next_c_r);
        transcript.append_fr("z1_next_ab_r", &self.z1_next_ab_r);
        transcript.append_fr("z1_next_ac_r", &self.z1_next_ac_r);
        transcript.append_fr("z1_next_bc_r", &self.z1_next_bc_r);
        transcript.append_fr("z1_next_abc_r", &self.z1_next_abc_r);
    }
}

pub struct PermGipaCommitments<E: PairingEngine> {
    // commit with vk
    pub z1_g_vk_comm: E::Fqk,
    pub z1_next_vk_comm: E::G2Affine,
    // inner product
    pub z1_g_r_ip: E::G1Affine,
    pub z1_g_a_r_ip: E::G1Affine,
    pub z1_g_b_r_ip: E::G1Affine,
    pub z1_g_c_r_ip: E::G1Affine,
    pub z1_g_ab_r_ip: E::G1Affine,
    pub z1_g_ac_r_ip: E::G1Affine,
    pub z1_g_bc_r_ip: E::G1Affine,
    pub z1_g_abc_r_ip: E::G1Affine,
    pub z1_next_r_ip: E::Fr,
    pub z1_next_a_r_ip: E::Fr,
    pub z1_next_b_r_ip: E::Fr,
    pub z1_next_c_r_ip: E::Fr,
    pub z1_next_ab_r_ip: E::Fr,
    pub z1_next_ac_r_ip: E::Fr,
    pub z1_next_bc_r_ip: E::Fr,
    pub z1_next_abc_r_ip: E::Fr,
}

impl<E, T> TranscriptAppender<E, T> for PermGipaCommitments<E>
where
    E: PairingEngine,
    T: TranscriptProtocol<E>,
{
    fn append_in_transcript(&self, transcript: &mut T) {
        transcript.append_gt("z1_g_vk_comm", &self.z1_g_vk_comm);
        transcript.append_g2("z1_next_vk_comm", &self.z1_next_vk_comm);
        transcript.append_g1("z1_g_r_ip", &self.z1_g_r_ip);
        transcript.append_g1("z1_g_a_r_ip", &self.z1_g_a_r_ip);
        transcript.append_g1("z1_g_b_r_ip", &self.z1_g_b_r_ip);
        transcript.append_g1("z1_g_c_r_ip", &self.z1_g_c_r_ip);
        transcript.append_g1("z1_g_ab_r_ip", &self.z1_g_ab_r_ip);
        transcript.append_g1("z1_g_ac_r_ip", &self.z1_g_ac_r_ip);
        transcript.append_g1("z1_g_bc_r_ip", &self.z1_g_bc_r_ip);
        transcript.append_g1("z1_g_abc_r_ip", &self.z1_g_abc_r_ip);
        transcript.append_fr("z1_next_r_ip", &self.z1_next_r_ip);
        transcript.append_fr("z1_next_a_r_ip", &self.z1_next_a_r_ip);
        transcript.append_fr("z1_next_b_r_ip", &self.z1_next_b_r_ip);
        transcript.append_fr("z1_next_c_r_ip", &self.z1_next_c_r_ip);
        transcript.append_fr("z1_next_ab_r_ip", &self.z1_next_ab_r_ip);
        transcript.append_fr("z1_next_ac_r_ip", &self.z1_next_ac_r_ip);
        transcript.append_fr("z1_next_bc_r_ip", &self.z1_next_bc_r_ip);
        transcript.append_fr("z1_next_abc_r_ip", &self.z1_next_abc_r_ip);
    }
}

pub(crate) fn perm_primitive_commit<E: PairingEngine>(
    z1_next: &[E::Fr],
    z1_g: &[E::G1Affine],
    vk: &VKey<E>,
) -> PermPrimitiveCommitments<E> {
    // commit with vk
    let z1_g_vk_comm = pairing::<E>(z1_g, vk);
    let z1_next_vk_comm = multiexponentiation(z1_next, vk);

    PermPrimitiveCommitments {
        z1_g_vk_comm,
        z1_next_vk_comm,
    }
}

pub(crate) fn perm_primitive_inner_product<E: PairingEngine>(
    r: &[E::Fr],
    a_r: &[E::Fr], // scaled a(z) * r
    b_r: &[E::Fr], // scaled b(z) * r
    c_r: &[E::Fr], // scaled c(z) * r
    ab_r: &[E::Fr], // scaled a(z)b(z) * r
    ac_r: &[E::Fr], // scaled a(z)c(z) * r
    bc_r: &[E::Fr], // scaled b(z)c(z) * r
    abc_r: &[E::Fr], // scaled a(z)b(z)c(z) * r
    z1_next: &[E::Fr],
    z1_g: &[E::G1Affine],
) -> PermPrimitiveInnerProducts<E> {
    let z1_g_r = multiexponentiation(r, z1_g);
    let z1_g_a_r = multiexponentiation(a_r, z1_g);
    let z1_g_b_r = multiexponentiation(b_r, z1_g);
    let z1_g_c_r = multiexponentiation(c_r, z1_g);
    let z1_g_ab_r = multiexponentiation(ab_r, z1_g);
    let z1_g_ac_r = multiexponentiation(ac_r, z1_g);
    let z1_g_bc_r = multiexponentiation(bc_r, z1_g);
    let z1_g_abc_r = multiexponentiation(abc_r, z1_g);
    let z1_next_r = scalars_inner_product(z1_next, r);
    let z1_next_a_r = scalars_inner_product(z1_next, a_r);
    let z1_next_b_r = scalars_inner_product(z1_next, b_r);
    let z1_next_c_r = scalars_inner_product(z1_next, c_r);
    let z1_next_ab_r = scalars_inner_product(z1_next, ab_r);
    let z1_next_ac_r = scalars_inner_product(z1_next, ac_r);
    let z1_next_bc_r = scalars_inner_product(z1_next, bc_r);
    let z1_next_abc_r = scalars_inner_product(z1_next, abc_r);

    PermPrimitiveInnerProducts {
        z1_g_r,
        z1_g_a_r,
        z1_g_b_r,
        z1_g_c_r,
        z1_g_ab_r,
        z1_g_ac_r,
        z1_g_bc_r,
        z1_g_abc_r,
        z1_next_r,
        z1_next_a_r,
        z1_next_b_r,
        z1_next_c_r,
        z1_next_ab_r,
        z1_next_ac_r,
        z1_next_bc_r,
        z1_next_abc_r,
    }
}

pub(crate) fn permutation_gipa_commit<E>(
    split: usize,
    r: &[E::Fr],
    a_r: &[E::Fr], // scaled a(z) * r
    b_r: &[E::Fr], // scaled b(z) * r
    c_r: &[E::Fr], // scaled c(z) * r
    ab_r: &[E::Fr], // scaled a(z)b(z) * r
    ac_r: &[E::Fr], // scaled a(z)c(z) * r
    bc_r: &[E::Fr], // scaled b(z)c(z) * r
    abc_r: &[E::Fr], // scaled a(z)b(z)c(z) * r
    z1_next: &[E::Fr],
    z1_g: &[E::G1Affine],
    vk: &VKey<E>,
) -> (PermGipaCommitments<E>, PermGipaCommitments<E>)
where
    E: PairingEngine,
{
    let (r_left, r_right) = r.split_at(split);

    let (a_r_left, a_r_right) = a_r.split_at(split);
    let (b_r_left, b_r_right) = b_r.split_at(split);
    let (c_r_left, c_r_right) = c_r.split_at(split);
    let (ab_r_left, ab_r_right) = ab_r.split_at(split);
    let (ac_r_left, ac_r_right) = ac_r.split_at(split);
    let (bc_r_left, bc_r_right) = bc_r.split_at(split);
    let (abc_r_left, abc_r_right) = abc_r.split_at(split);
    let (z1_next_left, z1_next_right) = z1_next.split_at(split);

    let (z1_g_left, z1_g_right) = z1_g.split_at(split);

    let (vk_left, vk_right) = vk.split_at(split);

    // Commit with vk
    // -------------------------------------------------------
    // (z1_g, vk)      |                  |
    let z1_g_vk_comm_left = pairing::<E>(z1_g_right, vk_left);
    let z1_g_vk_comm_right = pairing::<E>(z1_g_left, vk_right);
    // -------------------------------------------------------
    // (z1_next, vk)   |                  |
    let z1_next_vk_comm_left =
        multiexponentiation(z1_next_right, vk_left);
    let z1_next_vk_comm_right =
        multiexponentiation(z1_next_left, vk_right);

    // Inner product
    // -------------------------------------------------------
    //                 |                  |  (z1_g, r)
    let z1_g_r_ip_left =
        multiexponentiation(r_left, z1_g_right);
    let z1_g_r_ip_right =
        multiexponentiation(r_right, z1_g_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_g, a * r)
    let z1_g_a_r_ip_left =
        multiexponentiation(a_r_left, z1_g_right);
    let z1_g_a_r_ip_right =
        multiexponentiation(a_r_right, z1_g_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_g, b * r)
    let z1_g_b_r_ip_left =
        multiexponentiation(b_r_left, z1_g_right);
    let z1_g_b_r_ip_right =
        multiexponentiation(b_r_right, z1_g_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_g, c * r)
    let z1_g_c_r_ip_left =
        multiexponentiation(c_r_left, z1_g_right);
    let z1_g_c_r_ip_right =
        multiexponentiation(c_r_right, z1_g_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_g, ab * r)
    let z1_g_ab_r_ip_left =
        multiexponentiation(ab_r_left, z1_g_right);
    let z1_g_ab_r_ip_right =
        multiexponentiation(ab_r_right, z1_g_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_g, ac * r)
    let z1_g_ac_r_ip_left =
        multiexponentiation(ac_r_left, z1_g_right);
    let z1_g_ac_r_ip_right =
        multiexponentiation(ac_r_right, z1_g_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_g, bc * r)
    let z1_g_bc_r_ip_left =
        multiexponentiation(bc_r_left, z1_g_right);
    let z1_g_bc_r_ip_right =
        multiexponentiation(bc_r_right, z1_g_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_g, abc * r)
    let z1_g_abc_r_ip_left =
        multiexponentiation(abc_r_left, z1_g_right);
    let z1_g_abc_r_ip_right =
        multiexponentiation(abc_r_right, z1_g_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_next, r)
    let z1_next_r_ip_left =
        scalars_inner_product(r_left, z1_next_right);
    let z1_next_r_ip_right =
        scalars_inner_product(r_right, z1_next_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_next, a * r)
    let z1_next_a_r_ip_left =
        scalars_inner_product(a_r_left, z1_next_right);
    let z1_next_a_r_ip_right =
        scalars_inner_product(a_r_right, z1_next_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_next, b * r)
    let z1_next_b_r_ip_left =
        scalars_inner_product(b_r_left, z1_next_right);
    let z1_next_b_r_ip_right =
        scalars_inner_product(b_r_right, z1_next_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_next, c * r)
    let z1_next_c_r_ip_left =
        scalars_inner_product(c_r_left, z1_next_right);
    let z1_next_c_r_ip_right =
        scalars_inner_product(c_r_right, z1_next_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_next, ab * r)
    let z1_next_ab_r_ip_left =
        scalars_inner_product(ab_r_left, z1_next_right);
    let z1_next_ab_r_ip_right =
        scalars_inner_product(ab_r_right, z1_next_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_next, ac * r)
    let z1_next_ac_r_ip_left =
        scalars_inner_product(ac_r_left, z1_next_right);
    let z1_next_ac_r_ip_right =
        scalars_inner_product(ac_r_right, z1_next_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_next, bc * r)
    let z1_next_bc_r_ip_left =
        scalars_inner_product(bc_r_left, z1_next_right);
    let z1_next_bc_r_ip_right =
        scalars_inner_product(bc_r_right, z1_next_left);
    // -------------------------------------------------------
    //                 |                  |  (z1_next, abc * r)
    let z1_next_abc_r_ip_left =
        scalars_inner_product(abc_r_left, z1_next_right);
    let z1_next_abc_r_ip_right =
        scalars_inner_product(abc_r_right, z1_next_left);

    let left_comms = PermGipaCommitments {
        z1_g_vk_comm: z1_g_vk_comm_left,
        z1_next_vk_comm: z1_next_vk_comm_left,
        z1_g_r_ip: z1_g_r_ip_left,
        z1_g_a_r_ip: z1_g_a_r_ip_left,
        z1_g_b_r_ip: z1_g_b_r_ip_left,
        z1_g_c_r_ip: z1_g_c_r_ip_left,
        z1_g_ab_r_ip: z1_g_ab_r_ip_left,
        z1_g_ac_r_ip: z1_g_ac_r_ip_left,
        z1_g_bc_r_ip: z1_g_bc_r_ip_left,
        z1_g_abc_r_ip: z1_g_abc_r_ip_left,
        z1_next_r_ip: z1_next_r_ip_left,
        z1_next_a_r_ip: z1_next_a_r_ip_left,
        z1_next_b_r_ip: z1_next_b_r_ip_left,
        z1_next_c_r_ip: z1_next_c_r_ip_left,
        z1_next_ab_r_ip: z1_next_ab_r_ip_left,
        z1_next_ac_r_ip: z1_next_ac_r_ip_left,
        z1_next_bc_r_ip: z1_next_bc_r_ip_left,
        z1_next_abc_r_ip: z1_next_abc_r_ip_left,
    };
    let right_comms = PermGipaCommitments {
        z1_g_vk_comm: z1_g_vk_comm_right,
        z1_next_vk_comm: z1_next_vk_comm_right,
        z1_g_r_ip: z1_g_r_ip_right,
        z1_g_a_r_ip: z1_g_a_r_ip_right,
        z1_g_b_r_ip: z1_g_b_r_ip_right,
        z1_g_c_r_ip: z1_g_c_r_ip_right,
        z1_g_ab_r_ip: z1_g_ab_r_ip_right,
        z1_g_ac_r_ip: z1_g_ac_r_ip_right,
        z1_g_bc_r_ip: z1_g_bc_r_ip_right,
        z1_g_abc_r_ip: z1_g_abc_r_ip_right,
        z1_next_r_ip: z1_next_r_ip_right,
        z1_next_a_r_ip: z1_next_a_r_ip_right,
        z1_next_b_r_ip: z1_next_b_r_ip_right,
        z1_next_c_r_ip: z1_next_c_r_ip_right,
        z1_next_ab_r_ip: z1_next_ab_r_ip_right,
        z1_next_ac_r_ip: z1_next_ac_r_ip_right,
        z1_next_bc_r_ip: z1_next_bc_r_ip_right,
        z1_next_abc_r_ip: z1_next_abc_r_ip_right,
    };

    (left_comms, right_comms)
}
