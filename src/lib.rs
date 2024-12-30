//! Provides various methods for evaluating both spectral domain and spatial domain Green's
//! functions for time-harmonic electromagnetic problems.

use num::complex::Complex64;

pub mod material;
pub mod planar;
pub mod telegrapher;

pub type TETM = (Complex64, Complex64);

/// Spectral domain functions with possible tail arguments. Mainly used for Green's functions in
/// spectral domain, i.e. $\tilde{G}(k_p, \rho, z, z_p)$ $\tilde{G}: \C\times \R^3 \to \C$
pub trait SpFn: Sized {
    type InputTail;
    type Output;

    type CurryFOut: SpFn;
    type CurrySOut: SpFn;
    type CurryOOut: SpFn;

    /// Evaluate function
    fn eval(&self, kp: Complex64, tail: Self::InputTail) -> Self::Output;

    /// Optional functions for currying
    fn curry_freq(&self, _: f64) -> Option<Self::CurryFOut> { Option::None }
    fn curry_src<F: SpFn>(&self) -> Option<Self::CurrySOut> { Option::None }
    fn curry_obs<F: SpFn>(&self) -> Option<Self::CurryOOut> { Option::None }
}

impl SpFn for () {
    type InputTail = ();
    type Output = ();
    type CurryFOut= ();
    type CurrySOut = ();
    type CurryOOut = ();

    fn eval(&self, _: Complex64, _: Self::InputTail) -> Self::Output {
        ()
    }

}

/// Marker trait for spectral domain variable dependency
pub trait SpDep {}

/// Marker type for frequency dependency
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct FreqDep;
/// Marker type for source vertical position $z_p$ dependency
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct SrcZpDep;
/// Marker type for observer vertical position $z$ dependency
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct ObsZDep;
/// Marker type for planar distance $\rho$ dependency
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct RhoDep;

impl SpDep for FreqDep {}
impl SpDep for SrcZpDep {}
impl SpDep for ObsZDep {}
impl SpDep for () {}
