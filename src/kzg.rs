use ark_std::cfg_iter;
use ark_ec::{AffineCurve, ProjectiveCurve, msm::VariableBaseMSM};
use ark_ff::PrimeField;
use ark_poly::{UVPolynomial, univariate::DensePolynomial, Polynomial};
use num_traits::{One, Zero};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::error::Error;

pub type KZGOpening<G> = G;

pub(super) fn create_kzg_opening<G: AffineCurve>(
    powers: &[G],
    poly: &DensePolynomial<G::ScalarField>,
    point: G::ScalarField,
) -> Result<KZGOpening<G>, Error> {
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
    let coeffs: Vec<_> = cfg_iter!(&coeffs[num_leading_zeros..])
        .map(|s| s.into_repr())
        .collect();

    let w = VariableBaseMSM::multi_scalar_mul(
        &powers[num_leading_zeros..],
        &coeffs,
    );
    
    Ok(w.into_affine())
}