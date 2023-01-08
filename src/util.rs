
//! PLONKUP aggeration system

use ark_std::{cfg_iter, cfg_iter_mut};
use ark_ec::AffineRepr;
use ark_ff::Field;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Returns the vector used for the linear combination fo the inner pairing product
/// between A and B for the Groth16 aggregation: A^r * B. It is required as it
/// is not enough to simply prove the ipp of A*B, we need a random linear
/// combination of those.
pub(crate) fn structured_scalar_power<F: Field>(num: usize, s: &F) -> Vec<F> {
    let mut powers = vec![F::one()];
    for i in 1..num {
        powers.push(powers[i - 1] * s);
    }
    powers
}

pub(crate) fn compress_affines<G: AffineRepr>(
    vec: &mut Vec<G>,
    split: usize,
    scale: &G::ScalarField,
) {
    let (left, right) = vec.split_at_mut(split);
    assert_eq!(left.len(), right.len());
    
    cfg_iter_mut!(left)
        .zip(cfg_iter!(right))
        .for_each(|(left, right)| {
            *left = (*right * scale + *left).into();
        });

    vec.truncate(split);
}

pub(crate) fn compress_scalars<F: Field>(
    vec: &mut Vec<F>,
    split: usize,
    scale: &F,
) {
    let (left, right) = vec.split_at_mut(split);
    assert_eq!(left.len(), right.len());

    left.iter_mut()
        .zip(right.iter_mut())
        .for_each(|(l, r)| {
            r.mul_assign(scale);
            l.add_assign(*r);
        });

    vec.truncate(split);
}

/// It returns the evaluation of the polynomial $\prod (1 + x_{l-j}(rX)^{2j}$ at
/// the point z, where transcript contains the reversed order of all challenges (the x).
/// THe challenges must be in reversed order for the correct evaluation of the
/// polynomial in O(logn)
pub(crate) fn poly_eval_from_transcript<F: Field>(
    transcript: &[F],
    mut point: F,
    r: F,
) -> F {
    let mut res = F::one();
    for (i, x) in transcript.iter().enumerate() {
        if i == 0 {
            point.mul_assign(r);
        } else {
            point = point.square();
        }
        res.mul_assign(point * x + F::one());
    }

    res
}

// This method expects the coefficients in reverse order so transcript[i] =
// x_{l-j}.
// f(X) = \prod_{i=0}^{l-1} \left ( 1 + x_{l-i} (rX)^{2^i} \right )
pub(crate) fn poly_coeffs_from_transcript<F: Field>(
    transcript: &[F],
    mut r: F,
) -> Vec<F> {
    let mut coefficients = Vec::new();
    for (i, x) in transcript.iter().enumerate() {
        if i == 0 {
            coefficients.push(F::one());
            coefficients.push(r * x);
        } else {
            r = r.square();
            let shift = r * x;
            for j in 0..coefficients.len() {
                coefficients.push(coefficients[j] * shift);
            }
        }
    }

    coefficients
}
