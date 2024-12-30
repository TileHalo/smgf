//! Linear Isotropic and spatially dispersioneless material definition and helpers

use std::{f64::consts::PI, marker::PhantomData};

use num::complex::Complex64;
use scilib::constant::{EPSILON_0, MU_0};

use crate::{FreqDep, SpDep, TETM};


/// Sp domain propagation constant
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpPropagation<D: SpDep>{
    /// Frequency
    pub f: f64,
    mat: Material,
    dep: PhantomData<D>
}

/// Sp domain impedance
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpImpedance<D: SpDep>(SpPropagation<D>);

/// Sp domain reflection
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpReflection<D: SpDep, U: SpDep>(SpImpedance<D>, SpImpedance<U>);

impl<D: SpDep> SpPropagation<D> {
    pub fn new() -> Self {
        SpPropagation{f: 0.0, mat: Material::new(),dep: PhantomData }
    }
    fn with_material(mat: Material) -> Self {
        SpPropagation{f: 0.0, mat, dep: PhantomData }
    }
}

impl<D: SpDep> From<Material> for SpPropagation<D> {
    fn from(value: Material) -> Self {
        Self::with_material(value)
    }
}

impl crate::SpFn for SpPropagation<FreqDep> {
    type InputTail = f64;
    type Output = Complex64;
    type CurryFOut= SpPropagation<()>;
    type CurrySOut = ();
    type CurryOOut = ();

    fn eval(&self, kp: Complex64, tail: Self::InputTail) -> Self::Output {
        Complex64::sqrt(self.mat.cplx_prop_const(tail).powi(2) - kp.powi(2))
    }
    fn curry_freq(&self, f: f64) -> Option<Self::CurryFOut> {
        Some(SpPropagation{ f, mat: self.mat, dep: PhantomData})
    }
}

impl crate::SpFn for SpPropagation<()> {
    type InputTail = ();
    type Output = Complex64;
    type CurryFOut= ();
    type CurrySOut = ();
    type CurryOOut = ();

    fn eval(&self, kp: Complex64, _: Self::InputTail) -> Self::Output {
        Complex64::sqrt(self.mat.cplx_prop_const(self.f).powi(2) - kp.powi(2))
    }
}


impl<D: SpDep> SpImpedance<D> {
    pub fn new() -> Self {
        SpImpedance(SpPropagation::<D>::new())
    }
    fn with_material(mat: Material) -> Self {
        SpImpedance(SpPropagation::<D>::with_material(mat))
    }
}

impl<D: SpDep> From<Material> for SpImpedance<D> {
    fn from(value: Material) -> Self {
        Self::with_material(value)
    }
}

impl crate::SpFn for SpImpedance<FreqDep> {
    type InputTail = f64;
    type Output = TETM;
    type CurryFOut= SpImpedance<()>;
    type CurrySOut = ();
    type CurryOOut = ();

    fn eval(&self, kp: Complex64, f: Self::InputTail) -> Self::Output {
        let kz = self.0.mat.sp_cplx_wavenumber();
        let kzh = self.0.mat.sp_cplx_wavenumber();
        let (den, nume) = (
            EPSILON_0 * self.0.mat.cplx_rel_permittivity(f) * 2.0 * PI * f,
            MU_0 * self.0.mat.mu_r * 2.0 * PI * f,
        );
        let ze = kz.eval(kp, f) / den;
        let zh = nume / kzh.eval(kp, f);
        (ze, zh)
    }

    fn curry_freq(&self, f: f64) -> Option<Self::CurryFOut> {
        let mut b = Self::CurryFOut::from(self.0.mat);
        b.0.f = f;
        Some(b)
    }
}

impl crate::SpFn for SpImpedance<()> {
    type InputTail = ();
    type Output = TETM;
    type CurryFOut= ();
    type CurrySOut = ();
    type CurryOOut = ();

    fn eval(&self, kp: Complex64, _: Self::InputTail) -> Self::Output {
        let bulvaani = SpImpedance::<FreqDep>::from(self.0.mat);

        bulvaani.eval(kp, self.0.f)
    }
}

impl<D: SpDep, U: SpDep> SpReflection<D, U> {
    pub fn new() -> Self {
        SpReflection(SpImpedance::<D>::new(),SpImpedance::<U>::new())
    }
    fn with_materials(a: Material, b: Material) -> Self {
        SpReflection(SpImpedance::<D>::from(a), SpImpedance::<U>::from(b))
    }
}

impl<D: SpDep, U: SpDep> From<(Material, Material)> for SpReflection<D, U> {
    fn from((a, b): (Material, Material)) -> Self {
        Self::with_materials(a, b)
    }
}

impl crate::SpFn for SpReflection<FreqDep, FreqDep> {
    type InputTail = f64;
    type Output = TETM;
    type CurryFOut= SpReflection<(), ()>;
    type CurrySOut = ();
    type CurryOOut = ();

    fn eval(&self, kp: Complex64, f: Self::InputTail) -> Self::Output {

        let (ze1, zh1) = self.0.eval(kp, f);
        let (ze2, zh2) = self.1.eval(kp, f);

        ((ze2 - ze1)/(ze1 + ze2),
        (zh2 - zh1)/(zh1 + zh2))
    }

    fn curry_freq(&self, f: f64) -> Option<Self::CurryFOut> {
        let mut b = Self::CurryFOut::from((self.0.0.mat, self.1.0.mat));
        b.0.0.f = f;
        b.1.0.f = f;
        Some(b)
    }
}

impl crate::SpFn for SpReflection<(), ()> {
    type InputTail = ();
    type Output = TETM;
    type CurryFOut= ();
    type CurrySOut = ();
    type CurryOOut = ();

    fn eval(&self, kp: Complex64, _: Self::InputTail) -> Self::Output {
        let bulvaani = SpReflection::<FreqDep, FreqDep>::from((self.1.0.mat, self.1.0.mat));

        bulvaani.eval(kp, self.0.0.f)
    }
}


/// Various (electric) loss definitions for media.
/// In all documentation we have for complex permittivity $\varepsilon = \varepsilon' - j\varepsilon''$
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LossFormulation {
    /// Loss tangent $\tan\delta, that is $\tan\delta = \frac{\varepsilon''}{\varepsilon'}$.
    LossTangent(f64),
    /// Conductivity $\sigma$, that is $\varepsilon'' = j\frac{\sigma}{\omega}$ and thus lossess
    /// are dependent on frequency
    Conductivity(f64),
    /// Material is PEC
    PEC,
    /// Material is lossless dielectric, that is $\varepsilon'' = 0$.
    Dielectric,
}

/// A generic linear isotropic material
///
/// ## Usage:
/// Create either vacuum (constant `Material::VACUUM` or `Material::new()`), metal (`Material::metal()`)
/// or dielectric (`Material::dielectric()`).
///
/// Various properties can then be calculated:
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Material {
    /// Relative permittivity, $\varepsilon' = \varepsilon_0\varepsilon_r$
    pub eps_r: f64,
    /// Relative permeability, $\varepsilon' = \varepsilon_0\varepsilon_r$
    pub mu_r: f64,
    /// Material losses
    pub lossess: LossFormulation,
}

impl Material {
    /// Vacuum material
    pub const VACUUM: Self = Self {
        eps_r: 1.0,
        mu_r: 1.0,
        lossess: LossFormulation::Dielectric,
    };
    /// PEC
    pub const PEC: Self = Self {
        eps_r: 1.0,
        mu_r: 1.0,
        lossess: LossFormulation::PEC,
    };
    /// Create a new material with properties of vacuum
    pub fn new() -> Self {
        Self::VACUUM
    }
    /// Create a new metal (i.e. dielectric material with conductivity $\sigma$)
    pub fn metal(sigma: f64) -> Self {
        Self {
            eps_r: 1.0,
            mu_r: 1.0,
            lossess: LossFormulation::Conductivity(sigma),
        }
    }

    /// Create a new dielectric with loss tangent. Use `Some(val)` for lossy materials and `None`
    /// for lossless
    pub fn dielectric(eps_r: f64, tand: Option<f64>) -> Self {
        let lossess = match tand {
            Some(d) => LossFormulation::LossTangent(d),
            None => LossFormulation::Dielectric,
        };
        Self {
            eps_r,
            mu_r: 1.0,
            lossess,
        }
    }

    /// Calculate wavelength from frequency `f` (in Hz).
    pub fn wavelength(&self, f: f64) -> f64 {
        match self.lossess {
            LossFormulation::PEC => 0.0,
            _ => (scilib::constant::C / f64::sqrt(self.mu_r * self.eps_r)) / f,
        }
    }

    fn cplx_rel_permittivity(&self, f: f64) -> Complex64 {
        match self.lossess {
            LossFormulation::LossTangent(tand) => self.eps_r - Complex64::I * tand * self.eps_r,
            LossFormulation::Conductivity(sigma) => {
                self.eps_r - Complex64::I * sigma / (EPSILON_0 * 2.0 * PI * f)
            }
            LossFormulation::Dielectric => self.eps_r - Complex64::ZERO,
            LossFormulation::PEC => Complex64::I * f64::INFINITY,
        }
    }

    /// Complex relative permittivity
    pub fn cplx_permittivity(&self, f: f64) -> Complex64 {
        EPSILON_0 * self.cplx_rel_permittivity(f)
    }

    /// Complex wavenumber $\gamma = \alpha + j\beta$
    pub fn cplx_prop_const(&self, f: f64) -> Complex64 {
        2.0 * PI / (scilib::constant::C / f)
            * Complex64::sqrt(self.cplx_rel_permittivity(f) * self.mu_r)
    }

    /// Attenuation constant $\alpha$
    pub fn attenuation_const(&self, f: f64) -> f64 {
        self.cplx_prop_const(f).re
    }

    /// phase constant $\alpha$
    pub fn phase_const(&self, f: f64) -> f64 {
        self.cplx_prop_const(f).im
    }

    /// Complex wave impedance
    pub fn cplx_impedance(&self, f: f64) -> Complex64 {
        match self.lossess {
            LossFormulation::LossTangent(_) | LossFormulation::Conductivity(_) => {
                Complex64::I * MU_0 * 2.0 * PI * f / self.cplx_prop_const(f)
            }
            LossFormulation::Dielectric => {
                (f64::sqrt(self.mu_r / self.eps_r) * scilib::constant::Z_0).into()
            }
            LossFormulation::PEC => 0.0.into(),
        }
    }

    /// Complex reflection coefficient
    pub fn reflection_coeff(&self, b: Material, f: f64) -> Complex64 {
        let (z1, z2) = (self.cplx_impedance(f), b.cplx_impedance(f));

        (z2 - z1)/(z1 + z2)
    }

    /// Sp complex wavenumber
    pub fn sp_cplx_wavenumber(&self) -> SpPropagation<FreqDep>{
        SpPropagation::from(*self)
    }

    /// Return complex impedances in spectral domain, for electric field and magnetic field
    /// respectively
    pub fn sp_cplx_impedance(&self) -> SpImpedance<FreqDep> {
        SpImpedance::from(*self)
    }

    /// Sp reflection coefficient
    pub fn sp_reflection_coeff(&self, b: Material) -> SpReflection<FreqDep, FreqDep> {
        SpReflection::from((*self, b))
    }
}

impl Default for Material {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use assert_approx_eq::assert_approx_eq;
    use scilib::constant::C;
    use crate::*;

    use super::*;

    #[test]
    fn vacuum_test() {
        let vcm = Material::VACUUM;

        let freq = [1.0, 1.0e3, 1e6, 1e9, 5.0e9, 1.5e9, 5.0e12, 7.5e18];
        let ress = [
            C / 1.0,
            C / 1.0e3,
            C / 1e6,
            C / 1e9,
            C / 5.0e9,
            C / 1.5e9,
            C / 5.0e12,
            C / 7.5e18,
        ];

        assert_eq!(vcm.eps_r, 1.0);
        assert_eq!(vcm.mu_r, 1.0);

        for (f, res) in freq.iter().zip(ress.iter()) {
            // Wavelength
            let ans = vcm.wavelength(*f);
            let err = ans * 1.0e-15; // Due to exploding values force error into
            assert_approx_eq!(ans, res, err);

            // permittivity
            let ans = vcm.cplx_permittivity(*f);
            let err = ans.norm() * 1.0e-15; // Due to exploding values force error into
            assert_approx_eq!(ans.re, EPSILON_0, err);
            assert_approx_eq!(ans.im, 0.0, err);

            // wavenumber
            let ans = vcm.cplx_prop_const(*f);
            let err = ans.norm() * 1.0e-15; // Due to exploding values force error into
            assert_approx_eq!(ans.re, 2.0 * PI / res, err);
            assert_approx_eq!(ans.im, 0.0, err);

            let k = vcm.sp_cplx_wavenumber().curry_freq(*f).unwrap();
            let z = vcm.sp_cplx_impedance().curry_freq(*f).unwrap();
            let kpp: &[Complex64] = &[
                0.0.into(),
                1.0.into(),
                200.0.into(),
                1.0e9.into(),
                1.0e20.into(),
            ];
            for kp in kpp {
                let an = k.eval(*kp, ());
                let r = Complex64::sqrt((2.0 * PI / res).powi(2) - kp.powi(2));
                assert_approx_eq!(an.re, r.re, err);
                assert_approx_eq!(an.im, r.im, err);

                let (ze, zh) = z.eval(*kp, ());
                let re = r / (2.0 * PI * f * EPSILON_0);
                let eerr = an.norm() * 1.0e-5; // Due to exploding values force error into
                assert_approx_eq!(ze.re, re.re, eerr);
                assert_approx_eq!(ze.im, re.im, eerr);
                                               // something decent
                let rh = (2.0 * PI * f * MU_0) / r;
                let herr = an.norm() * 1.0e-5;
                assert_approx_eq!(zh.re, rh.re, herr);
                assert_approx_eq!(zh.im, rh.im, herr);
            }
        }
    }

    // Check up to permittivity, others should be fine
    #[test]
    fn metal_test() {
        let metal = Material::metal(1.0e8);

        let freq = [1.0, 1.0e3, 1e6, 1e9, 5.0e9, 1.5e9, 5.0e12, 7.5e18];
        let ress = [
            C / 1.0,
            C / 1.0e3,
            C / 1e6,
            C / 1e9,
            C / 5.0e9,
            C / 1.5e9,
            C / 5.0e12,
            C / 7.5e18,
        ];

        assert_eq!(metal.eps_r, 1.0);
        assert_eq!(metal.mu_r, 1.0);

        for (f, res) in freq.iter().zip(ress.iter()) {
            // Wavelength
            let ans = metal.wavelength(*f);
            let err = ans * 1.0e-15;
            assert_approx_eq!(ans, res, err);

            // permittivity
            let ans = metal.cplx_permittivity(*f);
            let err = ans.norm() * 1.0e-15;
            assert_approx_eq!(ans.re, EPSILON_0, err);
            assert_approx_eq!(ans.im, -1.0e8 / (2.0 * PI * f), err);
        }
    }

    #[test]
    fn lossless_dielectric_test() {
        let dielectric = Material::dielectric(4.0, None);

        let freq = [1.0, 1.0e3, 1e6, 1e9, 5.0e9, 1.5e9, 5.0e12, 7.5e18];
        let ress = [
            0.5 * C / 1.0,
            0.5 * C / 1.0e3,
            0.5 * C / 1e6,
            0.5 * C / 1e9,
            0.5 * C / 5.0e9,
            0.5 * C / 1.5e9,
            0.5 * C / 5.0e12,
            0.5 * C / 7.5e18,
        ];

        assert_eq!(dielectric.eps_r, 4.0);
        assert_eq!(dielectric.mu_r, 1.0);

        for (f, res) in freq.iter().zip(ress.iter()) {
            // Wavelength
            let ans = dielectric.wavelength(*f);
            let err = ans * 1.0e-15;
            assert_approx_eq!(ans, res, err);

            // permittivity
            let ans = dielectric.cplx_permittivity(*f);
            let err = ans.norm() * 1.0e-15;
            assert_approx_eq!(ans.re, 4.0 * EPSILON_0, err);
            assert_eq!(ans.im, 0.0);

            // wavenumber
            let ans = dielectric.cplx_prop_const(*f);
            let err = ans.norm() * 1.0e-15;
            assert_approx_eq!(ans.re, 2.0 * PI / res, err);
            assert_eq!(ans.im, 0.0);

            let k = dielectric.sp_cplx_wavenumber().curry_freq(*f).unwrap();
            let z = dielectric.sp_cplx_impedance().curry_freq(*f).unwrap();
            let kpp: &[Complex64] = &[
                0.0.into(),
                1.0.into(),
                200.0.into(),
                1.0e9.into(),
                1.0e20.into(),
            ];
            for kp in kpp {
                let an = k.eval(*kp, ());
                let r = Complex64::sqrt((2.0 * PI / res).powi(2) - kp.powi(2));
                let err = an.norm() * 1.0e-15;
                assert_approx_eq!(an.re, r.re, err);
                assert_approx_eq!(an.im, r.im, err);

                let (ze, zh) = z.eval(*kp, ());
                let rr = r / (2.0 * PI * f * EPSILON_0 * 4.0);
                let err = ze.norm() * 1.0e-15; // Due to exploding values force error into
                                               // something decent
                assert_approx_eq!(ze.re, rr.re, err);
                assert_approx_eq!(ze.im, rr.im, err);

                let rr = (2.0 * PI * f * MU_0) / r;
                let err = zh.norm() * 1.0e-10;
                assert_approx_eq!(zh.re, rr.re, err);
                assert_approx_eq!(zh.im, rr.im, err);
            }
        }
    }

    #[test]
    fn lossy_dielectric_test() {
        let dielectric = Material::dielectric(4.0, Some(0.15));

        let freq = [1.0, 1.0e3, 1e6, 1e9, 5.0e9, 1.5e9, 5.0e12, 7.5e18];
        let ress = [
            0.5 * C / 1.0,
            0.5 * C / 1.0e3,
            0.5 * C / 1e6,
            0.5 * C / 1e9,
            0.5 * C / 5.0e9,
            0.5 * C / 1.5e9,
            0.5 * C / 5.0e12,
            0.5 * C / 7.5e18,
        ];

        assert_eq!(dielectric.eps_r, 4.0);
        assert_eq!(dielectric.mu_r, 1.0);

        for (f, res) in freq.iter().zip(ress.iter()) {
            // Wavelength
            let ans = dielectric.wavelength(*f);
            let err = ans * 1.0e-12;
            assert_approx_eq!(ans, res, err);

            // permittivity
            let ans = dielectric.cplx_permittivity(*f);
            let err = ans.norm() * 1.0e-12;
            assert_approx_eq!(ans.re, 4.0 * EPSILON_0, err);
            assert_approx_eq!(ans.im, -4.0 * 0.15 * EPSILON_0, err);
        }
    }
}
