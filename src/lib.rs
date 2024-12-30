//! Provides various methods for evaluating both spectral domain and spatial domain Green's
//! functions for time-harmonic electromagnetic problems.

use num::complex::Complex64;

pub mod material;
pub mod planar;
pub mod telegrapher;

/// Spectral domain functions with possible tail arguments. Mainly used for Green's functions in
/// spectral domain, i.e. $\tilde{G}(k_p, \rho, z, z_p)$ $\tilde{G}: \C\times \R^3 \to \C$
pub trait SpectralFn: Sized {
    type InputTail;
    type Output;

    type CurryFOut: SpectralFn;
    type CurrySOut: SpectralFn;
    type CurryOOut: SpectralFn;

    /// Evaluate function
    fn eval(&self, kp: Complex64, tail: Self::InputTail) -> Self::Output;

    /// Optional functions for currying
    fn curry_freq(&self, _: f64) -> Option<Self::CurryFOut> { Option::None }
    fn curry_src<F: SpectralFn>(&self) -> Option<Self::CurrySOut> { Option::None }
    fn curry_obs<F: SpectralFn>(&self) -> Option<Self::CurryOOut> { Option::None }
}

impl SpectralFn for () {
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
pub trait SpectralDependency {}

/// Marker type for frequency dependency
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct FrequencyDep;
/// Marker type for source vertical position $z_p$ dependency
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct SrcZDep;
/// Marker type for observer vertical position $z$ dependency
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct ObsZDep;
/// Marker type for planar distance $\rho$ dependency
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct RhoDep;

impl SpectralDependency for FrequencyDep {}
impl SpectralDependency for SrcZDep {}
impl SpectralDependency for ObsZDep {}
impl SpectralDependency for () {}
