
// We have quite some functions that require quite some args by it's nature.
// It can be refactored but for now, we avoid these warns.
#![allow(clippy::too_many_arguments)]
// #![deny(missing_docs)]

mod util;
mod kzg;
mod commit;

pub mod srs;
pub mod error;
pub mod proof;
pub mod prove;
pub mod verify;
pub mod transcript;
