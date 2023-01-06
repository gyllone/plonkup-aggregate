pub mod wire;
pub mod permutation;
pub mod lookup;

use std::ops::{Deref, DerefMut};
use ark_std::{cfg_iter, cfg_iter_mut};
use ark_ec::{AffineCurve, ProjectiveCurve, PairingEngine, msm::VariableBaseMSM};
use ark_ff::{Field, PrimeField};
#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub type VKey<E> = Key<<E as PairingEngine>::G2Affine>;

pub type WKey<E> = Key<<E as PairingEngine>::G1Affine>;

/// Key is a generic commitment key that is instanciated with basis and powers.
#[derive(Debug, Clone)]
pub struct Key<G: AffineCurve>(Vec<G>);

impl<G: AffineCurve> Key<G> {
    pub(crate) fn scale(&self, s: &[G::ScalarField]) -> Self {
        assert_eq!(self.len(), s.len());

        let k: Vec<_> = cfg_iter!(&self)
            .zip(cfg_iter!(s))
            .map(|(g, s)| g.mul(*s).into_affine())
            .collect();

        Self(k)
    }

    pub(crate) fn compress(
        &mut self,
        split: usize,
        scale: &G::ScalarField,
    ) {
        let (left, right) = self.0.split_at_mut(split);
        assert_eq!(left.len(), right.len());

        cfg_iter_mut!(left)
            .zip(cfg_iter_mut!(right))
            .for_each(|(left, right)| {
                let mut r = right.mul(*scale);
                r.add_assign_mixed(left);
                *left = r.into_affine();
            });

        self.0.truncate(split);
    }
}

impl<G: AffineCurve> Deref for Key<G> {
    type Target = [G];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<G: AffineCurve> DerefMut for Key<G> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

///
fn pairing<E: PairingEngine>(
    g1: &[E::G1Affine],
    g2: &[E::G2Affine],
) -> E::Fqk {
    let pairs: Vec<_> = cfg_iter!(g1)
        .zip(cfg_iter!(g2))
        .map(|(g1, g2)| {
            (
                E::G1Prepared::from(*g1),
                E::G2Prepared::from(*g2),
            )
        })
        .collect();

    E::product_of_pairings(&pairs)
}

/// Need to be optimized with precomputed tables
fn multiexponentiation<G: AffineCurve>(
    scalars: &[G::ScalarField],
    bases: &[G],
) -> G {
    assert_eq!(scalars.len(), bases.len());

    let scalars: Vec<_> = cfg_iter!(scalars).map(|s| s.into_repr()).collect();
    let result = VariableBaseMSM::multi_scalar_mul(bases, &scalars);

    result.into_affine()
}

///
fn scalars_inner_product<F: Field>(
    scalars_a: &[F],
    scalars_b: &[F],
) -> F {
    assert_eq!(scalars_a.len(), scalars_b.len());

    cfg_iter!(scalars_a)
        .zip(cfg_iter!(scalars_b))
        .map(|(a, b)| *a * b)
        .sum()
}