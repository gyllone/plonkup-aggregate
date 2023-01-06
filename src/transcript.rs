use ark_ec::PairingEngine;

/// Transcript adds an abstraction over the Merlin transcript
/// For convenience
pub trait TranscriptProtocol<E: PairingEngine> {
    ///
    fn new(label: &'static str) -> Self;

    ///
    fn append_fr(&mut self, label: &'static str, item: &E::Fr);

    ///
    fn append_g1(&mut self, label: &'static str, item: &E::G1Affine);

    ///
    fn append_g2(&mut self, label: &'static str, item: &E::G2Affine);

    ///
    fn append_gt(&mut self, label: &'static str, item: &E::Fqk);

    /// Compute a `label`ed challenge variable.
    fn challenge_fr(&mut self, label: &'static str) -> E::Fr;
}

///
pub trait TranscriptAppender<E, T>
where
    E: PairingEngine,
    T: TranscriptProtocol<E>,
{
    ///
    fn append_in_transcript(&self, transcript: &mut T);
}
