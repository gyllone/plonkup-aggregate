pub mod wire;
pub mod permutation;
pub mod lookup;

use std::ops::{Deref, DerefMut};
use ark_std::{cfg_iter, cfg_iter_mut};
use ark_ec::{VariableBaseMSM, AffineRepr, pairing::{Pairing, PairingOutput}};
use ark_ff::Field;
#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub type VKey<E> = Key<<E as Pairing>::G2Affine>;

pub type WKey<E> = Key<<E as Pairing>::G1Affine>;

/// Key is a generic commitment key that is instanciated with basis and powers.
#[derive(Debug, Clone)]
pub struct Key<G: AffineRepr>(Vec<G>);

impl<G: AffineRepr> Key<G> {
    pub(crate) fn scale(&self, s: &[G::ScalarField]) -> Self {
        assert_eq!(self.len(), s.len());

        let k: Vec<_> = cfg_iter!(&self)
            .zip(cfg_iter!(s))
            .map(|(g, s)| (*g * s).into())
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
            .zip(cfg_iter!(right))
            .for_each(|(left, right)| {
                *left = (*right * scale + *left).into();
            });

        self.0.truncate(split);
    }
}

impl<G: AffineRepr> Deref for Key<G> {
    type Target = [G];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<G: AffineRepr> DerefMut for Key<G> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

///
fn pairing<E: Pairing>(
    g1: &[E::G1Affine],
    g2: &[E::G2Affine],
) -> PairingOutput<E> {
    let a: Vec<_> = cfg_iter!(g1)
        .map(|g1| E::G1Prepared::from(*g1))
        .collect();
    let b: Vec<_> = cfg_iter!(g2)
        .map(|g2| E::G2Prepared::from(*g2))
        .collect();

    E::multi_pairing(a, b)
}

/// Need to be optimized with precomputed tables
fn multiexponentiation<G: AffineRepr>(
    scalars: &[G::ScalarField],
    bases: &[G],
) -> G
where
    G::Group: VariableBaseMSM<MulBase = G>,
{
    assert_eq!(scalars.len(), bases.len());

    G::Group::msm_unchecked(bases, scalars).into()
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
