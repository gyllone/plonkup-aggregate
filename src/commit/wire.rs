use std::ops::AddAssign;
use ark_ec::pairing::Pairing;

use crate::{transcript::{TranscriptProtocol, TranscriptAppender}, proof::GipaFinalValue};
use super::*;

pub struct WirePrimitiveCommitments<E: Pairing> {
    // commit with vk
    pub a_vk: E::G2Affine,
    pub b_vk: E::G2Affine,
    pub a_g_vk: PairingOutput<E>,
    pub b_g_vk: PairingOutput<E>,
    pub c_g_vk: PairingOutput<E>,
    // commit with wk
    pub a_wk: E::G1Affine,
    pub b_wk: E::G1Affine,
    pub c_wk: E::G1Affine,
    pub ab_wk: E::G1Affine,
    pub ac_wk: E::G1Affine,
    pub bc_wk: E::G1Affine,
    pub abc_wk: E::G1Affine,
}

impl<E, T> TranscriptAppender<E, T> for WirePrimitiveCommitments<E>
where
    E: Pairing,
    T: TranscriptProtocol<E>,
{
    fn append_in_transcript(&self, transcript: &mut T) {
        transcript.append_g2("a_vk", &self.a_vk);
        transcript.append_g2("b_vk", &self.b_vk);
        transcript.append_gt("a_g_vk", &self.a_g_vk);
        transcript.append_gt("b_g_vk", &self.b_g_vk);
        transcript.append_gt("c_g_vk", &self.c_g_vk);
        transcript.append_g1("a_wk", &self.a_wk);
        transcript.append_g1("b_wk", &self.b_wk);
        transcript.append_g1("c_wk", &self.c_wk);
        transcript.append_g1("ab_wk", &self.ab_wk);
        transcript.append_g1("ac_wk", &self.ac_wk);
        transcript.append_g1("bc_wk", &self.bc_wk);
        transcript.append_g1("abc_wk", &self.abc_wk);
    }
}

pub(crate) fn wire_primitive_commit<E: Pairing>(
    a: &[E::ScalarField],
    b: &[E::ScalarField],
    c: &[E::ScalarField],
    ab: &[E::ScalarField],
    ac: &[E::ScalarField],
    bc: &[E::ScalarField],
    abc: &[E::ScalarField],
    a_g: &[E::G1Affine],
    b_g: &[E::G1Affine],
    c_g: &[E::G1Affine],
    vk: &VKey<E>,
    wk: &WKey<E>,
) -> WirePrimitiveCommitments<E> {
    WirePrimitiveCommitments {
        a_vk: multiexponentiation(a, vk),
        b_vk: multiexponentiation(b, vk),
        a_g_vk: pairing::<E>(a_g, vk),
        b_g_vk: pairing::<E>(b_g, vk),
        c_g_vk: pairing::<E>(c_g, vk),
        a_wk: multiexponentiation(a, wk),
        b_wk: multiexponentiation(b, wk),
        c_wk: multiexponentiation(c, wk),
        ab_wk: multiexponentiation(ab, wk),
        ac_wk: multiexponentiation(ac, wk),
        bc_wk: multiexponentiation(bc, wk),
        abc_wk: multiexponentiation(abc, wk),
    }
}

pub struct WirePrimitiveInnerProducts<E: Pairing> {
    pub a_r: E::ScalarField,
    pub b_r: E::ScalarField,
    pub c_r: E::ScalarField,
    pub ab_r: E::ScalarField,
    pub ac_r: E::ScalarField,
    pub bc_r: E::ScalarField,
    pub abc_r: E::ScalarField,
    pub a_g_r: E::G1Affine,
    pub b_g_r: E::G1Affine,
    pub c_g_r: E::G1Affine,
}

impl<E, T> TranscriptAppender<E, T> for WirePrimitiveInnerProducts<E>
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
        transcript.append_g1("a_g_r", &self.a_g_r);
        transcript.append_g1("b_g_r", &self.b_g_r);
        transcript.append_g1("c_g_r", &self.c_g_r);
    }
}

pub(crate) fn wire_primitive_inner_product<E: Pairing>(
    r: &[E::ScalarField],
    a_r: &[E::ScalarField], // scaled a(z) * r
    b_r: &[E::ScalarField], // scaled b(z) * r
    c_r: &[E::ScalarField], // scaled c(z) * r
    ab_r: &[E::ScalarField], // scaled a(z)b(z) * r
    ac_r: &[E::ScalarField], // scaled a(z)c(z) * r
    bc_r: &[E::ScalarField], // scaled b(z)c(z) * r
    abc_r: &[E::ScalarField], // scaled a(z)b(z)c(z) * r
    a_g: &[E::G1Affine],
    b_g: &[E::G1Affine],
    c_g: &[E::G1Affine],
) -> WirePrimitiveInnerProducts<E> {
    WirePrimitiveInnerProducts {
        a_r: a_r.iter().sum(),
        b_r: b_r.iter().sum(),
        c_r: c_r.iter().sum(),
        ab_r: ab_r.iter().sum(),
        ac_r: ac_r.iter().sum(),
        bc_r: bc_r.iter().sum(),
        abc_r: abc_r.iter().sum(),
        a_g_r: multiexponentiation(r, a_g),
        b_g_r: multiexponentiation(r, b_g),
        c_g_r: multiexponentiation(r, c_g),
    }
}

pub struct WireGipaCommitments<E: Pairing> {
    // commit with vk
    pub a_vk_comm: E::G2Affine,
    pub b_vk_comm: E::G2Affine,
    pub a_g_vk_comm: PairingOutput<E>,
    pub b_g_vk_comm: PairingOutput<E>,
    pub c_g_vk_comm: PairingOutput<E>,
    // commit with wk^r_inv
    pub a_r_wk_comm: E::G1Affine,
    pub b_r_wk_comm: E::G1Affine,
    pub c_r_wk_comm: E::G1Affine,
    pub ab_r_wk_comm: E::G1Affine,
    pub ac_r_wk_comm: E::G1Affine,
    pub bc_r_wk_comm: E::G1Affine,
    pub abc_r_wk_comm: E::G1Affine,
    // inner product
    pub a_r_ip: E::ScalarField,
    pub b_r_ip: E::ScalarField,
    pub c_r_ip: E::ScalarField,
    pub ab_r_ip: E::ScalarField,
    pub ac_r_ip: E::ScalarField,
    pub bc_r_ip: E::ScalarField,
    pub abc_r_ip: E::ScalarField,
    pub a_g_r_ip: E::G1Affine,
    pub b_g_r_ip: E::G1Affine,
    pub c_g_r_ip: E::G1Affine,
}

impl<E, T> TranscriptAppender<E, T> for WireGipaCommitments<E>
where
    E: Pairing,
    T: TranscriptProtocol<E>,
{
    fn append_in_transcript(&self, transcript: &mut T) {
        transcript.append_g2("a_vk_comm", &self.a_vk_comm);
        transcript.append_g2("b_vk_comm", &self.b_vk_comm);
        transcript.append_gt("a_g_vk_comm", &self.a_g_vk_comm);
        transcript.append_gt("b_g_vk_comm", &self.b_g_vk_comm);
        transcript.append_gt("c_g_vk_comm", &self.c_g_vk_comm);
        transcript.append_g1("a_r_wk_comm", &self.a_r_wk_comm);
        transcript.append_g1("b_r_wk_comm", &self.b_r_wk_comm);
        transcript.append_g1("c_r_wk_comm", &self.c_r_wk_comm);
        transcript.append_g1("ab_r_wk_comm", &self.ab_r_wk_comm);
        transcript.append_g1("ac_r_wk_comm", &self.ac_r_wk_comm);
        transcript.append_g1("bc_r_wk_comm", &self.bc_r_wk_comm);
        transcript.append_g1("abc_r_wk_comm", &self.abc_r_wk_comm);
        transcript.append_fr("a_r_ip", &self.a_r_ip);
        transcript.append_fr("b_r_ip", &self.b_r_ip);
        transcript.append_fr("c_r_ip", &self.c_r_ip);
        transcript.append_fr("ab_r_ip", &self.ab_r_ip);
        transcript.append_fr("ac_r_ip", &self.ac_r_ip);
        transcript.append_fr("bc_r_ip", &self.bc_r_ip);
        transcript.append_fr("abc_r_ip", &self.abc_r_ip);
        transcript.append_g1("a_g_r_ip", &self.a_g_r_ip);
        transcript.append_g1("b_g_r_ip", &self.b_g_r_ip);
        transcript.append_g1("c_g_r_ip", &self.c_g_r_ip);
    }
}

pub(crate) fn wire_gipa_commit<E: Pairing>(
    split: usize,
    r: &[E::ScalarField],
    a: &[E::ScalarField],
    b: &[E::ScalarField],
    a_r: &[E::ScalarField], // scaled a(z) * r
    b_r: &[E::ScalarField], // scaled b(z) * r
    c_r: &[E::ScalarField], // scaled c(z) * r
    ab_r: &[E::ScalarField], // scaled a(z)b(z) * r
    ac_r: &[E::ScalarField], // scaled a(z)c(z) * r
    bc_r: &[E::ScalarField], // scaled b(z)c(z) * r
    abc_r: &[E::ScalarField], // scaled a(z)b(z)c(z) * r
    a_g: &[E::G1Affine],
    b_g: &[E::G1Affine],
    c_g: &[E::G1Affine],
    vk: &VKey<E>,
    wk: &WKey<E>, // scaled key w^r_inv
) -> (WireGipaCommitments<E>, WireGipaCommitments<E>) {
    let (r_left, r_right) = r.split_at(split);
    let (a_left, a_right) = a.split_at(split);
    let (b_left, b_right) = b.split_at(split);
    
    let (a_r_left, a_r_right) = a_r.split_at(split);
    let (b_r_left, b_r_right) = b_r.split_at(split);
    let (c_r_left, c_r_right) = c_r.split_at(split);
    let (ab_r_left, ab_r_right) = ab_r.split_at(split);
    let (ac_r_left, ac_r_right) = ac_r.split_at(split);
    let (bc_r_left, bc_r_right) = bc_r.split_at(split);
    let (abc_r_left, abc_r_right) = abc_r.split_at(split);

    let (a_g_left, a_g_right) = a_g.split_at(split);
    let (b_g_left, b_g_right) = b_g.split_at(split);
    let (c_g_left, c_g_right) = c_g.split_at(split);

    let (vk_left, vk_right) = vk.split_at(split);
    let (wk_left, wk_right) = wk.split_at(split);

    // Commit with vk
    // -------------------------------------------------------
    // (a, vk)         |                  |
    let a_vk_comm_left = multiexponentiation(a_right, vk_left);
    let a_vk_comm_right = multiexponentiation(a_left, vk_right);
    // (b, vk)         |                  |
    let b_vk_comm_left = multiexponentiation(b_right, vk_left);
    let b_vk_comm_right = multiexponentiation(b_left, vk_right);
    // (a_g, vk)       |                  |
    let a_g_vk_comm_left = pairing::<E>(a_g_right, vk_left);
    let a_g_vk_comm_right = pairing::<E>(a_g_left, vk_right);
    // (b_g, vk)       |                  |
    let b_g_vk_comm_left = pairing::<E>(b_g_right, vk_left);
    let b_g_vk_comm_right = pairing::<E>(b_g_left, vk_right);
    // (c_g, vk)       |                  |
    let c_g_vk_comm_left = pairing::<E>(c_g_right, vk_left);
    let c_g_vk_comm_right = pairing::<E>(c_g_left, vk_right);

    // Commit with wk
    // -------------------------------------------------------
    //                 |   (a * r, wk)    |
    let a_r_wk_comm_left =
        multiexponentiation(a_r_left, wk_right);
    let a_r_wk_comm_right =
        multiexponentiation(a_r_right, wk_left);
    // -------------------------------------------------------
    //                 |   (b * r, wk)    |
    let b_r_wk_comm_left =
        multiexponentiation(b_r_left, wk_right);
    let b_r_wk_comm_right =
        multiexponentiation(b_r_right, wk_left);
    // -------------------------------------------------------
    //                 |   (c * r, wk)    |
    let c_r_wk_comm_left =
        multiexponentiation(c_r_left, wk_right);
    let c_r_wk_comm_right =
        multiexponentiation(c_r_right, wk_left);
    // -------------------------------------------------------
    //                 |   (ab * r, wk)   |
    let ab_r_wk_comm_left =
        multiexponentiation(ab_r_left, wk_right);
    let ab_r_wk_comm_right =
        multiexponentiation(ab_r_right, wk_left);
    // -------------------------------------------------------
    //                 |   (ac * r, wk)   |
    let ac_r_wk_comm_left =
        multiexponentiation(ac_r_left, wk_right);
    let ac_r_wk_comm_right =
        multiexponentiation(ac_r_right, wk_left);
    // -------------------------------------------------------
    //                 |   (bc * r, wk)   |  
    let bc_r_wk_comm_left =
        multiexponentiation(bc_r_left, wk_right);
    let bc_r_wk_comm_right =
        multiexponentiation(bc_r_right, wk_left);
    // -------------------------------------------------------
    //                 |   (abc * r, wk)  |
    let abc_r_wk_comm_left =
        multiexponentiation(abc_r_left, wk_right);
    let abc_r_wk_comm_right =
        multiexponentiation(abc_r_right, wk_left);

    // Inner product
    // -------------------------------------------------------
    //                 |                  |   (a, r) == (1, a * r)
    let a_r_ip_left = a_r_left.iter().sum();
    let a_r_ip_right = a_r_right.iter().sum();
    // -------------------------------------------------------
    //                 |                  |   (b, r) == (1, b * r)
    let b_r_ip_left = b_r_left.iter().sum();
    let b_r_ip_right = b_r_right.iter().sum();
    // -------------------------------------------------------
    //                 |                  |   (1, c * r)
    let c_r_ip_left = c_r_left.iter().sum();
    let c_r_ip_right = c_r_right.iter().sum();
    // -------------------------------------------------------
    //                 |                  |   (1, ab * r) == (a, b * r)
    let ab_r_ip_left = ab_r_left.iter().sum();
    let ab_r_ip_right = ab_r_right.iter().sum();
    // -------------------------------------------------------
    //                 |                  |  (1, ac * r) == (a, c * r)
    let ac_r_ip_left = ac_r_left.iter().sum();
    let ac_r_ip_right = ac_r_right.iter().sum();
    // -------------------------------------------------------
    //                 |                  |  (1, bc * r) == (b, c * r)
    let bc_r_ip_left = bc_r_left.iter().sum();
    let bc_r_ip_right = bc_r_right.iter().sum();
    // -------------------------------------------------------
    //                 |                  |  (1, abc * r) == (a, bc * r)
    let abc_r_ip_left = abc_r_left.iter().sum();
    let abc_r_ip_right = abc_r_right.iter().sum();
    // -------------------------------------------------------
    //                 |                  |   (a_g, r)
    let a_g_r_ip_left =
        multiexponentiation(r_left, a_g_right);
    let a_g_r_ip_right =
        multiexponentiation(r_right, a_g_left);
    // -------------------------------------------------------
    //                 |                  |   (b_g, r)
    let b_g_r_ip_left =
        multiexponentiation(r_left, b_g_right);
    let b_g_r_ip_right =
        multiexponentiation(r_right, b_g_left);
    // -------------------------------------------------------
    //                 |                  |   (c_g, r)
    let c_g_r_ip_left =
        multiexponentiation(r_left, c_g_right);
    let c_g_r_ip_right =
        multiexponentiation(r_right, c_g_left);

    let left_comms = WireGipaCommitments {
        a_vk_comm: a_vk_comm_left,
        b_vk_comm: b_vk_comm_left,
        a_g_vk_comm: a_g_vk_comm_left,
        b_g_vk_comm: b_g_vk_comm_left,
        c_g_vk_comm: c_g_vk_comm_left,
        a_r_wk_comm: a_r_wk_comm_left,
        b_r_wk_comm: b_r_wk_comm_left,
        c_r_wk_comm: c_r_wk_comm_left,
        ab_r_wk_comm: ab_r_wk_comm_left,
        ac_r_wk_comm: ac_r_wk_comm_left,
        bc_r_wk_comm: bc_r_wk_comm_left,
        abc_r_wk_comm: abc_r_wk_comm_left,
        a_r_ip: a_r_ip_left,
        b_r_ip: b_r_ip_left,
        c_r_ip: c_r_ip_left,
        ab_r_ip: ab_r_ip_left,
        ac_r_ip: ac_r_ip_left,
        bc_r_ip: bc_r_ip_left,
        abc_r_ip: abc_r_ip_left,
        a_g_r_ip: a_g_r_ip_left,
        b_g_r_ip: b_g_r_ip_left,
        c_g_r_ip: c_g_r_ip_left,
    };
    let right_comms = WireGipaCommitments {
        a_vk_comm: a_vk_comm_right,
        b_vk_comm: b_vk_comm_right,
        a_g_vk_comm: a_g_vk_comm_right,
        b_g_vk_comm: b_g_vk_comm_right,
        c_g_vk_comm: c_g_vk_comm_right,
        a_r_wk_comm: a_r_wk_comm_right,
        b_r_wk_comm: b_r_wk_comm_right,
        c_r_wk_comm: c_r_wk_comm_right,
        ab_r_wk_comm: ab_r_wk_comm_right,
        ac_r_wk_comm: ac_r_wk_comm_right,
        bc_r_wk_comm: bc_r_wk_comm_right,
        abc_r_wk_comm: abc_r_wk_comm_right,
        a_r_ip: a_r_ip_right,
        b_r_ip: b_r_ip_right,
        c_r_ip: c_r_ip_right,
        ab_r_ip: ab_r_ip_right,
        ac_r_ip: ac_r_ip_right,
        bc_r_ip: bc_r_ip_right,
        abc_r_ip: abc_r_ip_right,
        a_g_r_ip: a_g_r_ip_right,
        b_g_r_ip: b_g_r_ip_right,
        c_g_r_ip: c_g_r_ip_right,
    };

    (left_comms, right_comms)
}

pub(crate) struct WireFinalCommitments<E: Pairing> {
    // commit with vk
    a_vk_comm: E::G2,
    b_vk_comm: E::G2,
    a_g_vk_comm: PairingOutput<E>,
    b_g_vk_comm: PairingOutput<E>,
    c_g_vk_comm: PairingOutput<E>,
    // commit with wk
    a_r_wk_comm: E::G1,
    b_r_wk_comm: E::G1,
    c_r_wk_comm: E::G1,
    ab_r_wk_comm: E::G1,
    ac_r_wk_comm: E::G1,
    bc_r_wk_comm: E::G1,
    abc_r_wk_comm: E::G1,
    // inner product
    a_r_ip: E::ScalarField,
    b_r_ip: E::ScalarField,
    c_r_ip: E::ScalarField,
    ab_r_ip: E::ScalarField,
    ac_r_ip: E::ScalarField,
    bc_r_ip: E::ScalarField,
    abc_r_ip: E::ScalarField,
    a_g_r_ip: E::G1,
    b_g_r_ip: E::G1,
    c_g_r_ip: E::G1,
}

impl<E: Pairing> WireFinalCommitments<E> {
    pub(crate) fn from_primitive(
        comm: WirePrimitiveCommitments<E>,
        ip: WirePrimitiveInnerProducts<E>,
    ) -> Self {
        Self {
            a_vk_comm: comm.a_vk.into(),
            b_vk_comm: comm.b_vk.into(),
            a_g_vk_comm: comm.a_g_vk,
            b_g_vk_comm: comm.b_g_vk,
            c_g_vk_comm: comm.c_g_vk,
            a_r_wk_comm: comm.a_wk.into(),
            b_r_wk_comm: comm.b_wk.into(),
            c_r_wk_comm: comm.c_wk.into(),
            ab_r_wk_comm: comm.ab_wk.into(),
            ac_r_wk_comm: comm.ac_wk.into(),
            bc_r_wk_comm: comm.bc_wk.into(),
            abc_r_wk_comm: comm.abc_wk.into(),
            a_r_ip: ip.a_r,
            b_r_ip: ip.b_r,
            c_r_ip: ip.c_r,
            ab_r_ip: ip.ab_r,
            ac_r_ip: ip.ac_r,
            bc_r_ip: ip.bc_r,
            abc_r_ip: ip.abc_r,
            a_g_r_ip: ip.a_g_r.into(),
            b_g_r_ip: ip.b_g_r.into(),
            c_g_r_ip: ip.c_g_r.into(),
        }
    }

    pub(crate) fn merge(
        &mut self,
        left: WireGipaCommitments<E>,
        right: WireGipaCommitments<E>,
        x: &E::ScalarField,
        x_inv: &E::ScalarField,
    ) {
        self.a_vk_comm.add_assign(left.a_vk_comm * x + (right.a_vk_comm * x_inv));
        self.b_vk_comm.add_assign(left.b_vk_comm * x + (right.b_vk_comm * x_inv));
        self.a_g_vk_comm.add_assign(left.a_g_vk_comm * x + (right.a_g_vk_comm * x_inv));
        self.b_g_vk_comm.add_assign(left.b_g_vk_comm * x + (right.b_g_vk_comm * x_inv));
        self.c_g_vk_comm.add_assign(left.c_g_vk_comm * x + (right.c_g_vk_comm * x_inv));
        self.a_r_wk_comm.add_assign(left.a_r_wk_comm * x + (right.a_r_wk_comm * x_inv));
        self.b_r_wk_comm.add_assign(left.b_r_wk_comm * x + (right.b_r_wk_comm * x_inv));
        self.c_r_wk_comm.add_assign(left.c_r_wk_comm * x + (right.c_r_wk_comm * x_inv));
        self.ab_r_wk_comm.add_assign(left.ab_r_wk_comm * x + (right.ab_r_wk_comm * x_inv));
        self.ac_r_wk_comm.add_assign(left.ac_r_wk_comm * x + (right.ac_r_wk_comm * x_inv));
        self.bc_r_wk_comm.add_assign(left.bc_r_wk_comm * x + (right.bc_r_wk_comm * x_inv));
        self.abc_r_wk_comm.add_assign(left.abc_r_wk_comm * x + (right.abc_r_wk_comm * x_inv));
        self.a_r_ip.add_assign(left.a_r_ip * x + (right.a_r_ip * x_inv));
        self.b_r_ip.add_assign(left.b_r_ip * x + (right.b_r_ip * x_inv));
        self.c_r_ip.add_assign(left.c_r_ip * x + (right.c_r_ip * x_inv));
        self.ab_r_ip.add_assign(left.ab_r_ip * x + (right.ab_r_ip * x_inv));
        self.ac_r_ip.add_assign(left.ac_r_ip * x + (right.ac_r_ip * x_inv));
        self.bc_r_ip.add_assign(left.bc_r_ip * x + (right.bc_r_ip * x_inv));
        self.abc_r_ip.add_assign(left.abc_r_ip * x + (right.abc_r_ip * x_inv));
        self.a_g_r_ip.add_assign(left.a_g_r_ip * x + (right.a_g_r_ip * x_inv));
        self.b_g_r_ip.add_assign(left.b_g_r_ip * x + (right.b_g_r_ip * x_inv));
        self.c_g_r_ip.add_assign(left.c_g_r_ip * x + (right.c_g_r_ip * x_inv));
    }

    pub(crate) fn check(
        self,
        fv: &GipaFinalValue<E>,
        f1: &E::ScalarField,
        fr: &E::ScalarField,
    ) -> bool {
        // Check with vk
        // (a, vk)
        if fv.vk * fv.a != self.a_vk_comm {
        }
        // (b, vk)
        if fv.vk * fv.b != self.b_vk_comm {
        }
        // (a_g, vk)
        if E::pairing(fv.a_g, fv.vk) != self.a_g_vk_comm {
        }
        // (b_g, vk)
        if E::pairing(fv.b_g, fv.vk) != self.b_g_vk_comm {
        }
        // (c_g, vk)
        if E::pairing(fv.c_g, fv.vk) != self.c_g_vk_comm {
        }

        // Check with wk
        // (a * r, wk)
        if fv.wk * fv.a_r != self.a_r_wk_comm {
        }
        // (b * r, wk)
        if fv.wk * fv.b_r != self.b_r_wk_comm {
        }
        // (c * r, wk)
        if fv.wk * fv.c_r != self.c_r_wk_comm {
        }
        // (ab * r, wk)
        if fv.wk * fv.ab_r != self.ab_r_wk_comm {
        }
        // (ac * r, wk)
        if fv.wk * fv.ac_r != self.ac_r_wk_comm {
        }
        // (bc * r, wk)
        if fv.wk * fv.bc_r != self.bc_r_wk_comm {
        }
        // (abc * r, wk)
        if fv.wk * fv.abc_r != self.abc_r_wk_comm {
        }

        // Check inner product
        // (a, r) == (1, a * r)
        if fv.a * fr != self.a_r_ip || fv.a_r * f1 != self.a_r_ip {
        }
        // (b, r) == (1, b * r)
        if fv.b * fr != self.b_r_ip || fv.b_r * f1 != self.b_r_ip {
        }
        // (1, c * r)
        if fv.c_r * f1 != self.c_r_ip {
        }
        // (1, ab * r) == (a, b * r)
        if fv.ab_r * f1 != self.ab_r_ip || fv.a * fv.b_r != self.ab_r_ip {
        }
        // (1, ac * r) == (a, c * r)
        if fv.ac_r * f1 != self.ac_r_ip || fv.a * fv.c_r != self.ac_r_ip {
        }
        // (1, bc * r) == (b, c * r)
        if fv.bc_r * f1 != self.bc_r_ip || fv.b * fv.c_r != self.bc_r_ip {
        }
        // (1, abc * r) == (a, bc * r)
        if fv.abc_r * f1 != self.abc_r_ip || fv.a * fv.bc_r != self.abc_r_ip {
        }
        // (a_g, r)
        if fv.a_g * fr != self.a_g_r_ip {
        }
        // (b_g, r)
        if fv.b_g * fr != self.b_g_r_ip {
        }
        // (c_g, r)
        if fv.c_g * fr != self.c_g_r_ip {
        }

        true
    }
}
