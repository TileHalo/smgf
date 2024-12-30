//! Planar stratified media

use std::marker::PhantomData;

use num::complex::Complex64;

use crate::{material, FreqDep, SpDep, SpFn, TETM};

/// Generalized reflection coefficient $$\Gamma_{p,r} = \frac{\Gamma_{p+1, p} +
/// Gamma_{p+1,r}e^{-2jk_{z,p+1}(z_p - z_{p + 1})}}{1 + \Gamma_{p+1, p}
/// Gamma_{p+1,r}e^{-2jk_{z,p+1}(z_p - z_{p + 1})}}$$

#[derive(Debug, Clone, Default)]
pub struct SpGeneralizedReflection<D: SpDep> {
    f: f64,
    pub lprev: Option<Box<SpGeneralizedReflection<D>>>,
    pub rprev: Option<Box<SpGeneralizedReflection<D>>>,
    a: PlanarLayer,
    b: PlanarLayer,
    dep: PhantomData<D>
}

impl<D: SpDep> SpGeneralizedReflection<D> {

    fn eval_refl(&self, kp: Complex64, f: f64) -> (TETM, (Complex64,Complex64)) {
        let ((_, _), (rhe, rhh)) = match self.rprev {
           Some(ref r)  => r.eval_refl(kp, f),
            None => ((Complex64::ZERO, Complex64::ZERO),(Complex64::ZERO, Complex64::ZERO)),
        };
        let ((_, _), (lhe, lhh)) = match self.lprev {
           Some(ref l)  => l.eval_refl(kp, f),
            None => ((Complex64::ZERO, Complex64::ZERO), (Complex64::ZERO, Complex64::ZERO))
        };

        let ka = self.a.material.sp_cplx_wavenumber();
        let kb = self.a.material.sp_cplx_wavenumber();

        let (gre, grh) = self.a.material.sp_reflection_coeff(self.b.material).eval(kp, f);
        let (gle, glh) = self.b.material.sp_reflection_coeff(self.a.material).eval(kp, f);
        let lexp = Complex64::exp(-2.0*Complex64::I*ka.eval(kp, f)*self.a.height);
        let rexp = Complex64::exp(-2.0*Complex64::I*kb.eval(kp, f)*self.a.height);

        (((gle + lhe*lexp)/(1.0 + gle*lhe*lexp), (glh + lhh*lexp)/(1.0 + glh*lhh*lexp)),
        ((gre + rhe*rexp)/(1.0 + gre*rhe*rexp), (grh + rhh*rexp)/(1.0 + grh*rhh*rexp)))
    }
    pub fn new() -> Self {
        SpGeneralizedReflection {
            f: 0.0,
            lprev: None,
            rprev: None,
            a: PlanarLayer::default(),
            b: PlanarLayer::default(),
            dep: PhantomData,
        }
    }

    pub fn from_prev(
        lprev: Option<Box<SpGeneralizedReflection<D>>>,
        rprev: Option<Box<SpGeneralizedReflection<D>>>,
        a: PlanarLayer,
        b: PlanarLayer,
    ) -> Self {
        Self { f: 0.0, lprev, rprev, a, b, dep: PhantomData }
    }

}

impl<'a> SpFn for SpGeneralizedReflection<FreqDep> {
    type InputTail = f64;
    type Output = (TETM, TETM);
    type CurryFOut = SpGeneralizedReflection<()>;
    type CurrySOut = ();
    type CurryOOut = ();

    fn eval(&self, kp: Complex64, f: Self::InputTail) -> Self::Output {
        self.eval_refl(kp, f)
    }
}

impl<'a> SpFn for SpGeneralizedReflection<()> {
    type InputTail = f64;
    type Output = (TETM, TETM);
    type CurryFOut = ();
    type CurrySOut = ();
    type CurryOOut = ();

    fn eval(&self, kp: Complex64, _: Self::InputTail) -> Self::Output {
        self.eval_refl(kp, self.f)
    }

}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlanarLayer {
    pub height: f64,
    pub material: material::Material,
}

impl PlanarLayer {
    pub fn new() -> Self {
        PlanarLayer {
            height: 0.0,
            material: material::Material::VACUUM,
        }
    }

    pub fn with_material(height: f64, material: material::Material) -> Self {
        PlanarLayer { height, material }
    }
}

impl Default for PlanarLayer {
    fn default() -> Self {
        Self::new()
    }
}

pub struct PlanarStratifiedMedia<const N: usize> {
    pub layers: [PlanarLayer; N],
}

impl<const N: usize> PlanarStratifiedMedia<N> {
    pub fn new() -> Self {
        PlanarStratifiedMedia {
            layers: [PlanarLayer::default(); N],
        }
    }

    pub fn reflections<'a>(&self) -> [SpGeneralizedReflection<FreqDep>;N] {
        let mut refl = core::array::from_fn(|_| SpGeneralizedReflection::new());

        for (i, r) in refl.iter_mut().enumerate() {
            *r = SpGeneralizedReflection::from_prev(None, None, self.layers[i], self.layers[i + 1]);
        }
        let mut ran = 0..N;
        while let (Some(i), Some(j)) = (ran.next(), ran.next_back()) {
            if i != 0 {
                let p = Box::new(refl[i - 1].clone());
                refl[i].lprev = Some(p);
            }
            if j != N {
                let p = Box::new(refl[j + 1].clone());
                refl[j].rprev = Some(p);
            }
        }
        refl
    }
}

impl<const N: usize> Default for PlanarStratifiedMedia<N> {
    fn default() -> Self {
        Self::new()
    }
}
