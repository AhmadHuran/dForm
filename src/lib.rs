pub mod data;
pub mod psp;
pub mod linalg;

use special::Gamma;

pub use crate::{
    data::*,
    psp::PseudoGTH,
    linalg::*,
};

pub trait ELogGamma {
    type Out;
    fn elgamma(&self) -> Self::Out;
}

impl ELogGamma for f64 {
    type Out = f64;
    fn elgamma(&self) -> Self::Out {
        self.ln_gamma().0.exp()
    }
}
