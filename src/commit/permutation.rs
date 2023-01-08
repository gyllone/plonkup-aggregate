use std::ops::AddAssign;
use ark_ec::pairing::{Pairing, PairingOutput};

use crate::{
    transcript::{TranscriptProtocol, TranscriptAppender},
    proof::GipaFinalValue,
};
use super::*;

pub struct PermPrimitiveCommitments<E: Pairing> {
    // commit with vk
    pub z1_g_vk: PairingOutput<E>,
    pub z1_next_vk: E::G2Affine,
}

impl<E, T> TranscriptAppender<E, T> for PermPrimitiveCommitments<E>
where
    E: Pairing,
    T: TranscriptProtocol<E>,
{
    fn append_in_transcript(&self, transcript: &mut T) {
        transcript.append_gt("z1_g_vk", &self.z1_g_vk);
        transcript.append_g2("z1_next_vk", &self.z1_next_vk);
    }
}

pub(crate) fn perm_primitive_commit<E: Pairing>(
    z1_next: &[E::ScalarField],
    z1_g: &[E::G1Affine],
    vk: &VKey<E>,
) -> PermPrimitiveCommitments<E> {
    PermPrimitiveCommitments {
        z1_g_vk: pairing::<E>(z1_g, vk),
        z1_next_vk: multiexponentiation(z1_next, vk),
    }
}

pub struct PermPrimitiveInnerProducts<E: Pairing> {
    pub z1_g_r: E::G1Affine,
    pub z1_g_a_r: E::G1Affine,
    pub z1_g_b_r: E::G1Affine,
    pub z1_g_c_r: E::G1Affine,
    pub z1_g_ab_r: E::G1Affine,
    pub z1_g_ac_r: E::G1Affine,
    pub z1_g_bc_r: E::G1Affine,
    pub z1_g_abc_r: E::G1Affine,
    pub z1_next_r: E::ScalarField,
    pub z1_next_a_r: E::ScalarField,
    pub z1_next_b_r: E::ScalarField,
    pub z1_next_c_r: E::ScalarField,
    pub z1_next_ab_r: E::ScalarField,
    pub z1_next_ac_r: E::ScalarField,
    pub z1_next_bc_r: E::ScalarField,
    pub z1_next_abc_r: E::ScalarField,
}

pub(crate) fn perm_primitive_inner_product<E: Pairing>(
    r: &[E::ScalarField],
    a_r: &[E::ScalarField], // scaled a(z) * r
    b_r: &[E::ScalarField], // scaled b(z) * r
    c_r: &[E::ScalarField], // scaled c(z) * r
    ab_r: &[E::ScalarField], // scaled a(z)b(z) * r
    ac_r: &[E::ScalarField], // scaled a(z)c(z) * r
    bc_r: &[E::ScalarField], // scaled b(z)c(z) * r
    abc_r: &[E::ScalarField], // scaled a(z)b(z)c(z) * r
    z1_next: &[E::ScalarField],
    z1_g: &[E::G1Affine],
) -> PermPrimitiveInnerProducts<E> {
    PermPrimitiveInnerProducts {
        z1_g_r: multiexponentiation(r, z1_g),
        z1_g_a_r: multiexponentiation(a_r, z1_g),
        z1_g_b_r: multiexponentiation(b_r, z1_g),
        z1_g_c_r: multiexponentiation(c_r, z1_g),
        z1_g_ab_r: multiexponentiation(ab_r, z1_g),
        z1_g_ac_r: multiexponentiation(ac_r, z1_g),
        z1_g_bc_r: multiexponentiation(bc_r, z1_g),
        z1_g_abc_r: multiexponentiation(abc_r, z1_g),
        z1_next_r: scalars_inner_product(z1_next, r),
        z1_next_a_r: scalars_inner_product(z1_next, a_r),
        z1_next_b_r: scalars_inner_product(z1_next, b_r),
        z1_next_c_r: scalars_inner_product(z1_next, c_r),
        z1_next_ab_r: scalars_inner_product(z1_next, ab_r),
        z1_next_ac_r: scalars_inner_product(z1_next, ac_r),
        z1_next_bc_r: scalars_inner_product(z1_next, bc_r),
        z1_next_abc_r: scalars_inner_product(z1_next, abc_r),
    }
}

impl<E, T> TranscriptAppender<E, T> for PermPrimitiveInnerProducts<E>
where
    E: Pairing,
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

pub struct PermGipaCommitments<E: Pairing> {
    // commit with vk
    pub z1_g_vk_comm: PairingOutput<E>,
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
    pub z1_next_r_ip: E::ScalarField,
    pub z1_next_a_r_ip: E::ScalarField,
    pub z1_next_b_r_ip: E::ScalarField,
    pub z1_next_c_r_ip: E::ScalarField,
    pub z1_next_ab_r_ip: E::ScalarField,
    pub z1_next_ac_r_ip: E::ScalarField,
    pub z1_next_bc_r_ip: E::ScalarField,
    pub z1_next_abc_r_ip: E::ScalarField,
}

impl<E, T> TranscriptAppender<E, T> for PermGipaCommitments<E>
where
    E: Pairing,
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

pub(crate) fn permutation_gipa_commit<E>(
    split: usize,
    r: &[E::ScalarField],
    a_r: &[E::ScalarField], // scaled a(z) * r
    b_r: &[E::ScalarField], // scaled b(z) * r
    c_r: &[E::ScalarField], // scaled c(z) * r
    ab_r: &[E::ScalarField], // scaled a(z)b(z) * r
    ac_r: &[E::ScalarField], // scaled a(z)c(z) * r
    bc_r: &[E::ScalarField], // scaled b(z)c(z) * r
    abc_r: &[E::ScalarField], // scaled a(z)b(z)c(z) * r
    z1_next: &[E::ScalarField],
    z1_g: &[E::G1Affine],
    vk: &VKey<E>,
) -> (PermGipaCommitments<E>, PermGipaCommitments<E>)
where
    E: Pairing,
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

pub(crate) struct PermFinalCommitments<E: Pairing> {
    // commit with vk
    z1_g_vk_comm: PairingOutput<E>,
    z1_next_vk_comm: E::G2,
    // inner product
    z1_g_r_ip: E::G1,
    z1_g_a_r_ip: E::G1,
    z1_g_b_r_ip: E::G1,
    z1_g_c_r_ip: E::G1,
    z1_g_ab_r_ip: E::G1,
    z1_g_ac_r_ip: E::G1,
    z1_g_bc_r_ip: E::G1,
    z1_g_abc_r_ip: E::G1,
    z1_next_r_ip: E::ScalarField,
    z1_next_a_r_ip: E::ScalarField,
    z1_next_b_r_ip: E::ScalarField,
    z1_next_c_r_ip: E::ScalarField,
    z1_next_ab_r_ip: E::ScalarField,
    z1_next_ac_r_ip: E::ScalarField,
    z1_next_bc_r_ip: E::ScalarField,
    z1_next_abc_r_ip: E::ScalarField,
}

impl<E: Pairing> PermFinalCommitments<E> {
    pub(crate) fn from_primitive(
        comm: PermPrimitiveCommitments<E>,
        ip: PermPrimitiveInnerProducts<E>,
    ) -> Self {
        Self {
            z1_g_vk_comm: comm.z1_g_vk,
            z1_next_vk_comm: comm.z1_next_vk.into(),
            z1_g_r_ip: ip.z1_g_r.into(),
            z1_g_a_r_ip: ip.z1_g_a_r.into(),
            z1_g_b_r_ip: ip.z1_g_b_r.into(),
            z1_g_c_r_ip: ip.z1_g_c_r.into(),
            z1_g_ab_r_ip: ip.z1_g_ab_r.into(),
            z1_g_ac_r_ip: ip.z1_g_ac_r.into(),
            z1_g_bc_r_ip: ip.z1_g_bc_r.into(),
            z1_g_abc_r_ip: ip.z1_g_abc_r.into(),
            z1_next_r_ip: ip.z1_next_r,
            z1_next_a_r_ip: ip.z1_next_a_r,
            z1_next_b_r_ip: ip.z1_next_b_r,
            z1_next_c_r_ip: ip.z1_next_c_r,
            z1_next_ab_r_ip: ip.z1_next_ab_r,
            z1_next_ac_r_ip: ip.z1_next_ac_r,
            z1_next_bc_r_ip: ip.z1_next_bc_r,
            z1_next_abc_r_ip: ip.z1_next_abc_r,
        }
    }

    pub(crate) fn merge(
        &mut self,
        left: PermGipaCommitments<E>,
        right: PermGipaCommitments<E>,
        x: &E::ScalarField,
        x_inv: &E::ScalarField,
    ) {
        self.z1_g_vk_comm.add_assign(left.z1_g_vk_comm * x + (right.z1_g_vk_comm * x_inv));
        self.z1_next_vk_comm.add_assign(left.z1_next_vk_comm * x + (right.z1_next_vk_comm * x_inv));
        self.z1_g_r_ip.add_assign(left.z1_g_r_ip * x + (right.z1_g_r_ip * x_inv));
        self.z1_g_a_r_ip.add_assign(left.z1_g_a_r_ip * x + (right.z1_g_a_r_ip * x_inv));
        self.z1_g_b_r_ip.add_assign(left.z1_g_b_r_ip * x + (right.z1_g_b_r_ip * x_inv));
        self.z1_g_c_r_ip.add_assign(left.z1_g_c_r_ip * x + (right.z1_g_c_r_ip * x_inv));
        self.z1_g_ab_r_ip.add_assign(left.z1_g_ab_r_ip * x + (right.z1_g_ab_r_ip * x_inv));
        self.z1_g_ac_r_ip.add_assign(left.z1_g_ac_r_ip * x + (right.z1_g_ac_r_ip * x_inv));
        self.z1_g_bc_r_ip.add_assign(left.z1_g_bc_r_ip * x + (right.z1_g_bc_r_ip * x_inv));
        self.z1_g_abc_r_ip.add_assign(left.z1_g_abc_r_ip * x + (right.z1_g_abc_r_ip * x_inv));
        self.z1_next_r_ip.add_assign(left.z1_next_r_ip * x + (right.z1_next_r_ip * x_inv));
        self.z1_next_a_r_ip.add_assign(left.z1_next_a_r_ip * x + (right.z1_next_a_r_ip * x_inv));
        self.z1_next_b_r_ip.add_assign(left.z1_next_b_r_ip * x + (right.z1_next_b_r_ip * x_inv));
        self.z1_next_c_r_ip.add_assign(left.z1_next_c_r_ip * x + (right.z1_next_c_r_ip * x_inv));
        self.z1_next_ab_r_ip.add_assign(left.z1_next_ab_r_ip * x + (right.z1_next_ab_r_ip * x_inv));
        self.z1_next_ac_r_ip.add_assign(left.z1_next_ac_r_ip * x + (right.z1_next_ac_r_ip * x_inv));
        self.z1_next_bc_r_ip.add_assign(left.z1_next_bc_r_ip * x + (right.z1_next_bc_r_ip * x_inv));
        self.z1_next_abc_r_ip.add_assign(left.z1_next_abc_r_ip * x + (right.z1_next_abc_r_ip * x_inv));
    }

    pub(crate) fn check(
        self,
        fv: &GipaFinalValue<E>,
        f1: &E::ScalarField,
        fr: &E::ScalarField,
    ) -> bool {
        // Check with vk
        // (z1_g, vk)
        if E::pairing(fv.z1_g, fv.vk) != self.z1_g_vk_comm {
        }
        // (z1_next, vk)
        if fv.vk * fv.z1_next != self.z1_next_vk_comm {
        }
        
        // Check inner product
        // (z1_g, r)
        if fv.z1_g * fr != self.z1_g_r_ip {
        }
        // (z1_g, a_r)
        if fv.z1_g * fv.a_r != self.z1_g_a_r_ip {
        }
        // (z1_g, b_r)
        if fv.z1_g * fv.b_r != self.z1_g_b_r_ip {
        }
        // (z1_g, c_r)
        if fv.z1_g * fv.c_r != self.z1_g_c_r_ip {
        }
        // (z1_g, ab_r)
        if fv.z1_g * fv.ab_r != self.z1_g_ab_r_ip {
        }
        // (z1_g, ac_r)
        if fv.z1_g * fv.ac_r != self.z1_g_ac_r_ip {
        }
        // (z1_g, bc_r)
        if fv.z1_g * fv.bc_r != self.z1_g_bc_r_ip {
        }
        // (z1_g, abc_r)
        if fv.z1_g * fv.abc_r != self.z1_g_abc_r_ip {
        }
        // (z1_next, r)
        if fv.z1_next * fr != self.z1_next_r_ip {
        }
        // (z1_next, a_r)
        if fv.z1_next * fv.a_r != self.z1_next_a_r_ip {
        }
        // (z1_next, b_r)
        if fv.z1_next * fv.b_r != self.z1_next_b_r_ip {
        }
        // (z1_next, c_r)
        if fv.z1_next * fv.c_r != self.z1_next_c_r_ip {
        }
        // (z1_next, ab_r)
        if fv.z1_next * fv.ab_r != self.z1_next_ab_r_ip {
        }
        // (z1_next, ac_r)
        if fv.z1_next * fv.ac_r != self.z1_next_ac_r_ip {
        }
        // (z1_next, bc_r)
        if fv.z1_next * fv.bc_r != self.z1_next_bc_r_ip {
        }
        // (z1_next, abc_r)
        if fv.z1_next * fv.abc_r != self.z1_next_abc_r_ip {
        }

        true
    }
}
