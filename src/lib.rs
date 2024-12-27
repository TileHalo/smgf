use num::complex::Complex64;

pub mod material;
pub mod planar;
pub mod telegrapher;

macro_rules! zip {
    ($x: expr) => ($x);
    ($x: expr, $($y: expr), +) => (
        $x.into_iter().zip(
            zip!($($y), +))
    )
}

/// Curries the spectral frequency dependent Green's function with frequency
pub fn curry_frequency(
    fdsgf: impl Fn(Complex64, f64, f64, f64, f64) -> Complex64,
    f: f64,
) -> impl Fn(Complex64, f64, f64, f64) -> Complex64 {
    move |kp, p, z, zp| fdsgf(kp, p, z, zp, f)
}
/// Curries the spectral frequency dependent Green's function with frequencies in the iterator
pub fn curry_frequency_vect<I>(
    fdsgf: impl Fn(Complex64, f64, f64, f64, f64) -> Complex64 + Clone,
    freqs: I,
) -> impl Iterator
where
    I: IntoIterator<Item = f64>,
{
    let mut v = Vec::new();
    for f in freqs {
        let func = fdsgf.clone();
        v.push(curry_frequency(func, f));
    }
    v.into_iter()
}

/// Curries the spectral Green's function with source height $z_p$
pub fn curry_zp(
    fdsgf: impl Fn(Complex64, f64, f64, f64) -> Complex64,
    zp: f64,
) -> impl Fn(Complex64, f64, f64) -> Complex64 {
    move |kp, p, z| fdsgf(kp, p, z, zp)
}
/// Curries the spectral Green's function with source heights $z_p$ from iterator
pub fn curry_zp_vect<I>(
    fdsgf: impl Fn(Complex64, f64, f64, f64) -> Complex64 + Clone,
    zps: I,
) -> impl Iterator
where
    I: IntoIterator<Item = f64>,
{
    let mut v = Vec::new();
    for f in zps {
        let func = fdsgf.clone();
        v.push(curry_zp(func, f));
    }
    v.into_iter()
}

/// Curries the spectral Green's function with observer height $z$
pub fn curry_z(
    fdsgf: impl Fn(Complex64, f64, f64, f64) -> Complex64,
    z: f64,
) -> impl Fn(Complex64, f64, f64) -> Complex64 {
    move |kp, p, zp| fdsgf(kp, p, z, zp)
}
/// Curries the spectral Green's function with observer heights $z$ from iterator
pub fn curry_z_vect<I>(
    fdsgf: impl Fn(Complex64, f64, f64, f64) -> Complex64 + Clone,
    zs: I,
) -> impl Iterator
where
    I: IntoIterator<Item = f64>,
{
    let mut v = Vec::new();
    for f in zs {
        let func = fdsgf.clone();
        v.push(curry_z(func, f));
    }
    v.into_iter()
}

/// Curries the spectral Green's function with observer height $z$ and source heigh $z_p$
pub fn curry_z_zp(
    fdsgf: impl Fn(Complex64, f64, f64, f64) -> Complex64,
    z: f64,
    zp: f64,
) -> impl Fn(Complex64, f64) -> Complex64 {
    move |kp, p| fdsgf(kp, p, z, zp)
}
/// Curries the spectral Green's function with source heights $zp_p$ from iterator
pub fn curry_z_zp_vect<I, J>(
    fdsgf: impl Fn(Complex64, f64, f64, f64) -> Complex64 + Clone,
    zs: I,
    zps: J,
) -> impl Iterator
where
    I: IntoIterator<Item = f64>,
    J: IntoIterator<Item = f64>,
{
    let mut v = Vec::new();
    for (z, zp) in zs.into_iter().zip(zps.into_iter()) {
        let func = fdsgf.clone();
        v.push(curry_z_zp(func, z, zp));
    }
    v.into_iter()
}

/// Curries the spectral Green's function with observer height $z$ and source heigh $z_p$
pub fn curry_spatial(
    fdsgf: impl Fn(Complex64, f64, f64, f64) -> Complex64,
    z: f64,
    zp: f64,
    p: f64,
) -> impl Fn(Complex64) -> Complex64 {
    move |kp| fdsgf(kp, p, z, zp)
}
/// Curries the spectral Green's function with source heights $zp_p$ from iterator
pub fn curry_spatial_vect<I, J, U>(
    fdsgf: impl Fn(Complex64, f64, f64, f64) -> Complex64 + Clone,
    zs: I,
    zps: J,
    ps: U,
) -> impl Iterator
where
    I: IntoIterator<Item = f64>,
    J: IntoIterator<Item = f64>,
    U: IntoIterator<Item = f64>,
{
    let mut v = Vec::new();
    for (z, (zp, p)) in zip!(zs, zps, ps) {
        let func = fdsgf.clone();
        v.push(curry_spatial(func, z, zp, p));
    }
    v.into_iter()
}
