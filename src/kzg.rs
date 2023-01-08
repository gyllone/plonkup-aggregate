use ark_ec::{VariableBaseMSM, AffineRepr, CurveGroup, pairing::Pairing};
use ark_poly::{univariate::DensePolynomial, Polynomial, DenseUVPolynomial};
use num_traits::{One, Zero};

use crate::error::Error;

pub type KZGOpening<G> = G;

pub(crate) fn create_kzg_opening<G: AffineRepr>(
    powers: &[G],
    poly: &DensePolynomial<G::ScalarField>,
    point: G::ScalarField,
) -> Result<KZGOpening<G>, Error>
where
    G::Group: VariableBaseMSM<MulBase = G>,
{
    if poly.degree() + 1 > powers.len() {
        return Err(Error::PolynomialDegreeTooLarge);
    }

    // (p(x) - p(z)) / (x - z)
    // Observe that this quotient does not change with z because
    // p(z) is the remainder term. We can therefore omit p(z) when computing the quotient.
    let divisor_poly =
        DensePolynomial::from_coefficients_vec(vec![-point, G::ScalarField::one()]);
    let quotient_poly = poly / &divisor_poly;

    let coeffs = quotient_poly.coeffs();
    let mut num_leading_zeros = 0;
    while num_leading_zeros < coeffs.len() && coeffs[num_leading_zeros].is_zero() {
        num_leading_zeros += 1;
    }
    let res = G::Group::msm_unchecked(
        &powers[num_leading_zeros..],
        &coeffs[num_leading_zeros..],
    );
    
    Ok(res.into())
}

pub(crate) fn verify_kzg_opening_on_g1<E: Pairing>(
    g1_base: E::G1Affine,
    g2_base: E::G2Affine,
    tau: E::G2Affine,
    comm: E::G1Affine,
    pi: E::G1Affine,
    poly_eval: E::ScalarField,
    point: E::ScalarField,
) -> Result<(), Error> {
    // $ e([h(z)]_1, [\tau]_2-z\cdot[1]_2)=e([f(\tau)]_1-f(z)\cdot[1]_1, [1]_2) $
    let left = E::pairing(pi, tau + (g2_base * -point));
    let right = E::pairing(comm + (g1_base * -poly_eval), g2_base);

    if left != right {
        return Err(Error::PairingCheckFailure);
    } else {
        Ok(())
    }
}

pub(crate) fn verify_kzg_opening_on_g2<E: Pairing>(
    g1_base: E::G1Affine,
    g2_base: E::G2Affine,
    tau: E::G1Affine,
    comm: E::G2Affine,
    pi: E::G2Affine,
    poly_eval: E::ScalarField,
    point: E::ScalarField,
) -> Result<(), Error> {
    // $ e([\tau]_1-z\cdot[1]_1,[h(z)]_2)=e([1]_1, [f(\tau)]_2-f(z)\cdot[1]_2) $
    let left = E::pairing(tau + (g1_base * -point), pi);
    let right = E::pairing(g1_base, comm + (g2_base * -poly_eval));

    if left != right {
        return Err(Error::PairingCheckFailure);
    } else {
        Ok(())
    }
}
