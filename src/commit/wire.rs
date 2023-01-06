use ark_ec::{AffineCurve, ProjectiveCurve, PairingEngine};
use ark_ff::{Field, PrimeField};

use crate::transcript::{TranscriptProtocol, TranscriptAppender};
use super::*;

pub struct WirePrimitiveCommitments<E: PairingEngine> {
    // commit with vk
    pub a_vk_comm: E::G2Affine,
    pub b_vk_comm: E::G2Affine,
    pub a_g_vk_comm: E::Fqk,
    pub b_g_vk_comm: E::Fqk,
    pub c_g_vk_comm: E::Fqk,
    // commit with wk
    pub a_wk_comm: E::G1Affine,
    pub b_wk_comm: E::G1Affine,
    pub c_wk_comm: E::G1Affine,
    pub ab_wk_comm: E::G1Affine,
    pub ac_wk_comm: E::G1Affine,
    pub bc_wk_comm: E::G1Affine,
    pub abc_wk_comm: E::G1Affine,
}

impl<E, T> TranscriptAppender<E, T> for WirePrimitiveCommitments<E>
where
    E: PairingEngine,
    T: TranscriptProtocol<E>,
{
    fn append_in_transcript(&self, transcript: &mut T) {
        transcript.append_g2("a_vk_comm", &self.a_vk_comm);
        transcript.append_g2("b_vk_comm", &self.b_vk_comm);
        transcript.append_gt("a_g_vk_comm", &self.a_g_vk_comm);
        transcript.append_gt("b_g_vk_comm", &self.b_g_vk_comm);
        transcript.append_gt("c_g_vk_comm", &self.c_g_vk_comm);
        transcript.append_g1("a_r_wk_comm", &self.a_wk_comm);
        transcript.append_g1("b_r_wk_comm", &self.b_wk_comm);
        transcript.append_g1("c_r_wk_comm", &self.c_wk_comm);
        transcript.append_g1("ab_r_wk_comm", &self.ab_wk_comm);
        transcript.append_g1("ac_r_wk_comm", &self.ac_wk_comm);
        transcript.append_g1("bc_r_wk_comm", &self.bc_wk_comm);
        transcript.append_g1("abc_r_wk_comm", &self.abc_wk_comm);
    }
}

pub struct WirePrimitiveInnerProducts<E: PairingEngine> {
    pub a_r: E::Fr,
    pub b_r: E::Fr,
    pub c_r: E::Fr,
    pub ab_r: E::Fr,
    pub ac_r: E::Fr,
    pub bc_r: E::Fr,
    pub abc_r: E::Fr,
    pub a_g_r: E::G1Affine,
    pub b_g_r: E::G1Affine,
    pub c_g_r: E::G1Affine,
}

impl<E, T> TranscriptAppender<E, T> for WirePrimitiveInnerProducts<E>
where
    E: PairingEngine,
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

pub struct WireGipaCommitments<E: PairingEngine> {
    // commit with vk
    pub a_vk_comm: E::G2Affine,
    pub b_vk_comm: E::G2Affine,
    pub a_g_vk_comm: E::Fqk,
    pub b_g_vk_comm: E::Fqk,
    pub c_g_vk_comm: E::Fqk,
    // commit with wk^r_inv
    pub a_r_wk_comm: E::G1Affine,
    pub b_r_wk_comm: E::G1Affine,
    pub c_r_wk_comm: E::G1Affine,
    pub ab_r_wk_comm: E::G1Affine,
    pub ac_r_wk_comm: E::G1Affine,
    pub bc_r_wk_comm: E::G1Affine,
    pub abc_r_wk_comm: E::G1Affine,
    // inner product
    pub a_r_ip: E::Fr,
    pub b_r_ip: E::Fr,
    pub c_r_ip: E::Fr,
    pub ab_r_ip: E::Fr,
    pub ac_r_ip: E::Fr,
    pub bc_r_ip: E::Fr,
    pub abc_r_ip: E::Fr,
    pub a_g_r_ip: E::G1Affine,
    pub b_g_r_ip: E::G1Affine,
    pub c_g_r_ip: E::G1Affine,
}

impl<E: PairingEngine> WireGipaCommitments<E> {
    pub(crate) fn scale(&self, shift: E::Fr) -> Self {
        let shift_repr = shift.into_repr();

        let a_vk_comm = self.a_vk_comm.mul(shift).into_affine();
        let b_vk_comm = self.b_vk_comm.mul(shift).into_affine();
        let a_g_vk_comm = self.a_g_vk_comm.pow(&shift_repr);
        let b_g_vk_comm = self.b_g_vk_comm.pow(&shift_repr);
        let c_g_vk_comm = self.c_g_vk_comm.pow(&shift_repr);

        let a_r_wk_comm = self.a_r_wk_comm.mul(shift).into_affine();
        let b_r_wk_comm = self.b_r_wk_comm.mul(shift).into_affine();
        let c_r_wk_comm = self.c_r_wk_comm.mul(shift).into_affine();
        let ab_r_wk_comm = self.ab_r_wk_comm.mul(shift).into_affine();
        let ac_r_wk_comm = self.ac_r_wk_comm.mul(shift).into_affine();
        let bc_r_wk_comm = self.bc_r_wk_comm.mul(shift).into_affine();
        let abc_r_wk_comm = self.abc_r_wk_comm.mul(shift).into_affine();

        let a_r_ip = self.a_r_ip * shift;
        

    }
}

impl<E, T> TranscriptAppender<E, T> for WireGipaCommitments<E>
where
    E: PairingEngine,
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

pub(crate) fn wire_primitive_commit<E: PairingEngine>(
    a: &[E::Fr],
    b: &[E::Fr],
    c: &[E::Fr],
    ab: &[E::Fr],
    ac: &[E::Fr],
    bc: &[E::Fr],
    abc: &[E::Fr],
    a_g: &[E::G1Affine],
    b_g: &[E::G1Affine],
    c_g: &[E::G1Affine],
    vk: &VKey<E>,
    wk: &WKey<E>,
) -> WirePrimitiveCommitments<E> {
    // commit with vk
    let a_vk_comm = multiexponentiation(a, vk);
    let b_vk_comm = multiexponentiation(b, vk);
    let a_g_vk_comm = pairing::<E>(a_g, vk);
    let b_g_vk_comm = pairing::<E>(b_g, vk);
    let c_g_vk_comm = pairing::<E>(c_g, vk);

    // commit with wk
    let a_wk_comm = multiexponentiation(a, wk);
    let b_wk_comm = multiexponentiation(b, wk);
    let c_wk_comm = multiexponentiation(c, wk);
    let ab_wk_comm = multiexponentiation(ab, wk);
    let ac_wk_comm = multiexponentiation(ac, wk);
    let bc_wk_comm = multiexponentiation(bc, wk);
    let abc_wk_comm = multiexponentiation(abc, wk);

    WirePrimitiveCommitments {
        a_vk_comm,
        b_vk_comm,
        a_g_vk_comm,
        b_g_vk_comm,
        c_g_vk_comm,
        a_wk_comm,
        b_wk_comm,
        c_wk_comm,
        ab_wk_comm,
        ac_wk_comm,
        bc_wk_comm,
        abc_wk_comm,
    }
}

pub(crate) fn wire_primitive_inner_product<E: PairingEngine>(
    r: &[E::Fr],
    a_r: &[E::Fr], // scaled a(z) * r
    b_r: &[E::Fr], // scaled b(z) * r
    c_r: &[E::Fr], // scaled c(z) * r
    ab_r: &[E::Fr], // scaled a(z)b(z) * r
    ac_r: &[E::Fr], // scaled a(z)c(z) * r
    bc_r: &[E::Fr], // scaled b(z)c(z) * r
    abc_r: &[E::Fr], // scaled a(z)b(z)c(z) * r
    a_g: &[E::G1Affine],
    b_g: &[E::G1Affine],
    c_g: &[E::G1Affine],
) -> WirePrimitiveInnerProducts<E> {
    let a_r = a_r.iter().sum();
    let b_r = b_r.iter().sum();
    let c_r = c_r.iter().sum();
    let ab_r = ab_r.iter().sum();
    let ac_r = ac_r.iter().sum();
    let bc_r = bc_r.iter().sum();
    let abc_r = abc_r.iter().sum();
    let a_g_r = multiexponentiation(r, a_g);
    let b_g_r = multiexponentiation(r, b_g);
    let c_g_r = multiexponentiation(r, c_g);

    WirePrimitiveInnerProducts {
        a_r,
        b_r,
        c_r,
        ab_r,
        ac_r,
        bc_r,
        abc_r,
        a_g_r,
        b_g_r,
        c_g_r,
    }
}

pub(crate) fn wire_gipa_commit<E>(
    split: usize,
    r: &[E::Fr],
    a: &[E::Fr],
    b: &[E::Fr],
    a_r: &[E::Fr], // scaled a(z) * r
    b_r: &[E::Fr], // scaled b(z) * r
    c_r: &[E::Fr], // scaled c(z) * r
    ab_r: &[E::Fr], // scaled a(z)b(z) * r
    ac_r: &[E::Fr], // scaled a(z)c(z) * r
    bc_r: &[E::Fr], // scaled b(z)c(z) * r
    abc_r: &[E::Fr], // scaled a(z)b(z)c(z) * r
    a_g: &[E::G1Affine],
    b_g: &[E::G1Affine],
    c_g: &[E::G1Affine],
    vk: &VKey<E>,
    wk: &WKey<E>, // scaled key w^r_inv
) -> (WireGipaCommitments<E>, WireGipaCommitments<E>)
where
    E: PairingEngine,
{
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
