use ark_ec::pairing::{Pairing, PairingOutput};

/// Transcript adds an abstraction over the Merlin transcript
/// For convenience
pub trait TranscriptProtocol<E: Pairing> {
    ///
    fn new(label: &'static str) -> Self;

    ///
    fn append_fr(&mut self, label: &'static str, item: &E::ScalarField);

    ///
    fn append_g1(&mut self, label: &'static str, item: &E::G1Affine);

    ///
    fn append_g2(&mut self, label: &'static str, item: &E::G2Affine);

    ///
    fn append_gt(&mut self, label: &'static str, item: &PairingOutput<E>);

    /// Compute a `label`ed challenge variable.
    fn challenge_fr(&mut self, label: &'static str) -> E::ScalarField;
}

///
pub trait TranscriptAppender<E, T>
where
    E: Pairing,
    T: TranscriptProtocol<E>,
{
    ///
    fn append_in_transcript(&self, transcript: &mut T);
}
