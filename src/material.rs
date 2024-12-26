
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
    Dielectric
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
    pub lossess: LossFormulation
}

impl Material {
    pub const VACUUM: Self = Self { eps_r: 1.0, mu_r: 1.0, lossess:  LossFormulation::Dielectric  };
    pub const PEC: Self = Self { eps_r: 1.0, mu_r: 1.0, lossess:  LossFormulation::PEC  };
    /// Create a new material with properties of vacuum
    pub fn new() -> Self {
        Self::VACUUM
    }
    /// Create a new metal (i.e. dielectric material with conductivity $\sigma$)
    pub fn metal(sigma: f64) -> Self {
        Self { eps_r: 1.0, mu_r: 1.0, lossess:  LossFormulation::Conductivity(sigma)  }
    }

    /// Create a new dielectric with loss tangent. Use `Some(val)` for lossy materials and `None`
    /// for lossless
    pub fn dielectric(eps_r: f64, tand: Option<f64>) -> Self {
        let lossess = match tand {
            Some(d) => LossFormulation::LossTangent(d),
            None => LossFormulation::Dielectric
        };
        Self { eps_r, mu_r: 1.0, lossess}
    }

    /// Calculate wavelength from frequency `f` (in Hz).
    pub fn wavelength(&self, f: f64) -> f64 {
        (scilib::constant::C/f64::sqrt(self.mu_r*self.eps_r))/f
    }
}

impl Default for Material {
    fn default() -> Self {
        Self::new()
    }
}
